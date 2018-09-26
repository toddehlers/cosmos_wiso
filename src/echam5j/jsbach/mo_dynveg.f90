MODULE mo_dynveg
!
! Computation of pft fractions per gridbox - scheme adopted from lpj
!

  USE mo_jsbach_grid,      ONLY: grid_type, domain_type
  USE mo_kind,             ONLY: dp
  USE mo_exception,        ONLY: finish, message, message_text
  USE mo_linked_list,      ONLY: t_stream
  USE mo_mpi,              ONLY: p_parallel_io, p_io, p_parallel, p_bcast
  USE mo_climbuf,          ONLY: climbuf
  USE mo_netCDF,           ONLY: file_info, nf_max_name
  USE mo_io,               ONLY: IO_READ
  USE mo_land_surface,     ONLY: scale_cover_fract, fract_small

  IMPLICIT NONE

  ! === BEGIN OF PUBLIC PART ======================================================================================================

  TYPE dynveg_type  !! contains the state variables
     INTEGER          :: ntiles
     REAL(dp),POINTER :: bio_exist(:,:) 
     REAL(dp),POINTER :: litter_moisture(:)
     REAL(dp),POINTER :: burned_fpc(:,:)
     REAL(dp),POINTER :: damaged_fpc(:,:)    
     REAL(dp),POINTER :: act_fpc(:,:) 
     REAL(dp),POINTER :: pot_fpc(:,:) 
     REAL(dp),POINTER :: bare_fpc(:) 
     REAL(dp),POINTER :: desert_fpc(:)
     REAL(dp),POINTER :: carbon_2_LeafLitterPool(:)
     REAL(dp),POINTER :: carbon_2_WoodLitterPool(:)
     REAL(dp),POINTER :: carbon_2_slowSoilPool_damage(:)
     REAL(dp),POINTER :: carbon_2_slowSoilPool_fire(:)
     REAL(dp),POINTER :: carbon_2_atmos(:)
  END TYPE dynveg_type
  TYPE(dynveg_type), SAVE :: dynveg

  TYPE dynveg_params_type
     INTEGER           :: nlct                         !! number of landcover types found in dynveg landcover library file
     LOGICAL,  POINTER :: dynamic_PFT(:)               !! indicates those PFTs which shall take part in the vegetation dynamics
     LOGICAL,  POINTER :: woody_PFT(:)                 !! indicates those PFTs which are of woody type (in contrast to grasses)
     REAL(dp), POINTER :: bclimit_min_cold_mmtemp(:)   !! PFT-specific minimum coldest monthly mean temperature
     REAL(dp), POINTER :: bclimit_max_cold_mmtemp(:)   !! PFT-specific maximum coldest monthly mean temperature
     REAL(dp), POINTER :: bclimit_max_warm_mmtemp(:)   !! PFT-specific upper limit of warmest-month temperature
     REAL(dp), POINTER :: bclimit_min_temprange(:)     !! PFT-specific 20-year average min warmest - coldest month temperature range
     REAL(dp), POINTER :: bclimit_min_gdd(:)           !! PFT-specific minimum growing degree days (at or above 5 deg C)
     REAL(dp), POINTER :: gdd_base(:)                  !! PFT-specific GDD base
     REAL(dp), POINTER :: upper_tlim(:)                !! PFT-specific base to calculate GDD_upper_tlim 
     REAL(dp), POINTER :: tau_pft(:)                   !! PFT-specific time scale
  END TYPE dynveg_params_type
  TYPE(dynveg_params_type), SAVE :: dynveg_params

  TYPE dynveg_options_type
     LOGICAL                :: dynveg_all         !! competition between all PFTs, not just between woody PFTs
     LOGICAL                :: dynveg_feedback    !! activate climate feedback of the dynamic vegetation
     LOGICAL                :: init_running_means !! initialize running means of the climate buffer
     LOGICAL                :: read_climbuf       !! Read long term climate variables from file
     LOGICAL                :: read_fpc           !! Read fractional plant cover from file
     CHARACTER(nf_max_name) :: climbuf_file_name
     CHARACTER(nf_max_name) :: fpc_file_name
  END TYPE dynveg_options_type
  TYPE(dynveg_options_type), SAVE :: dynveg_options

  PUBLIC :: dynveg_type
  PUBLIC :: dynveg_params_type
  PUBLIC :: dynveg_options_type
  PUBLIC :: update_dynveg 
  PUBLIC :: init_dynveg
  PUBLIC :: config_dynveg
  PUBLIC :: fpc_daily
  PUBLIC :: potential_tree_fpc
  PUBLIC :: desert_fraction
  PUBLIC :: fpc_to_cover_fract
  PUBLIC :: cover_fract_to_fpc
  PUBLIC :: scale_fpc

  PUBLIC :: frac_wood_2_litter_wind, frac_wood_2_litter_fire, frac_wood_2_atmos

  REAL(dp), PARAMETER :: frac_wood_2_litter_wind = 0.5_dp   !! fraction of wood transferred to the litter pool by wind break
  REAL(dp), PARAMETER :: frac_wood_2_litter_fire = 0.25_dp  !! fraction of wood transferred to the litter pool by fire
  REAL(dp), PARAMETER :: frac_wood_2_atmos     = 0.25_dp    !! fraction of wood carbon emitted to the atmosphere by fire

  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ==============================================================================
  PRIVATE

  REAL(dp), PARAMETER :: sum_npp_min           = 1.e-12_dp  !! minimum total npp for tree cover, kg C/m2
  REAL(dp), PARAMETER :: act_fpc_min           = 0.005_dp   !! minimum actual FPC to be considered to calculate potential FPC
  REAL(dp), PARAMETER :: tree_fpc_max          = 1.0_dp     !! maximum FPC of woody PFTs
  REAL(dp), PARAMETER :: tolerance             = 1.e-9_dp   !! margin by which the sum of actual FPC may exceed 1
  REAL(dp), PARAMETER :: accelerate_vegetation = 1.0_dp     !! factor by which vegetation dynamics are accelerated (should be 1 for production runs) 
  REAL(dp), PARAMETER :: npp_nonlinearity      = 1.5_dp     !! parameter controlling the non-linearity in the dynamic equation with respect to NPP
  REAL(dp), PARAMETER :: desert_extend         = 0.65_dp    !! parameter controlling the extend of desert (the lower the value the more desert)
  REAL(dp), PARAMETER :: desert_margin         = 1.5_dp     !! parameter controlling the transition from vegetated land to desert
                                                            !! (the higher the value the sharper the transition)
  REAL(dp), PARAMETER :: tau_desert            = 50._dp     !! time constant by which the extension of deserts is adapted
  REAL(dp), PARAMETER :: life_to_estab         = 20._dp     !! ratio of life-time of plants to their establishment time scale
                                                            !! (only used, when competiton of woody types and grass is active)
  REAL(dp), PARAMETER :: fire_litter_threshold = 16.67_dp   !! minimal amount of litter [mol(C)/m^2(grid box)] for fire
  REAL(dp), PARAMETER :: fire_rel_hum_threshold= 70._dp     !! maximal relative humidity for fire
  REAL(dp), PARAMETER :: fire_minimum_woody    = 0.002_dp   !! minimal fraction of act_fpc of woody PFT to be burned each year
  REAL(dp), PARAMETER :: fire_minimum_grass    = 0.006_dp   !! minimal fraction of act_fpc of grass PFT to be burned each year
  REAL(dp), PARAMETER :: fire_tau_woody        = 7.0_dp     !! return period of fire for woody PFT [year] at 0 % relative humidity
  REAL(dp), PARAMETER :: fire_tau_grass        = 2.0_dp     !! return period of fire for grass PFT [year] at 0 % relative humidity
  REAL(dp), PARAMETER :: wind_threshold        = 2.2_dp     !! factor by which the previous day maximum wind speed must be larger than
                                                            !! the climatological daily maximum wind speed to allow any wind damage
  REAL(dp), PARAMETER :: wind_damage_scale     = 4.e-03_dp  !! scaling factor for wind damage

  REAL(dp), POINTER :: max_green_bio(:,:)            !! maximum value of green biomass within a year
  REAL(dp), POINTER :: sum_green_bio_memory(:)       !! vegetated fraction calculated from green biomass

  LOGICAL, SAVE      :: module_initialized = .false. !! Signifies whether module has been initialized
  INTEGER, SAVE      :: ntiles
  INTEGER, SAVE      :: nsoil

  TYPE(t_stream),  POINTER :: IO_dynveg    !! Memory stream for dynveg model state

  TYPE(file_info),   SAVE  :: fpc_file               !! Input file for FPC
  REAL(dp), POINTER, SAVE  :: init_act_fpc(:,:)      !! initial values for act_fpc if read from file
  REAL(dp), POINTER, SAVE  :: init_bare_fpc(:)       !! initial values for bare_fpc if read from file
  REAL(dp), POINTER, SAVE  :: init_desert_fpc(:)     !! initial values for desert_fpc if read from file

! !VERSION CONTROL:
  CHARACTER(len=*), PARAMETER :: version = '$ID$'


! !GLOBAL VARIABLES:

! initialization of bioclimatic limits from LPJ lookup table
! initial values for FPC: from original vegetation map of JSBACH
 
!      6  flammability threshold
!      8  fire resistance nindex
!
!     BIOCLIMATIC LIMITS
!
!     28 minimum coldest monthly mean temperature
!     29 maximum coldest monthly mean temperature
!     30 minimum growing degree days (at or above 5 deg C)
!     31 upper limit of temperature of the warmest month 
!     32 lower limit of growth efficiency (g/m2)
!     
!     PARAMETERS ADDED LATER
!
!     33 GDD base
!     34 20-year average min warmest - coldest month temperature range
!     10 wind damage threshold (m/s)

! Plant Functional Types (PFT) currently regarded in this scheme (and correspondence to LPJ PFTs):

! JSBACH: 1. tropical broadleaved evergreen forest       - LPJ type 1
! JSBACH: 2. tropical deciduous broadleaved forest       - LPJ type 2
! JSBACH: 3. temper./boreal evergreen forest             - LPJ types 4 (extended cold tolerance), 3 and 6
! JSBACH: 4. temper./boreal deciduous forest             - LPJ types 5, 7  and 8
! JSBACH: 5. raingreen shrubs                            - LPJ types 1, 3, 4
! JSBACH: 6. cold shrubs                                 - LPJ type 8 (tundra)
! JSBACH: 7. C3 perennial grass	                         - LPJ type 9
! JSBACH: 8. C4 perennial grass                          - LPJ type 10 
! JSBACH: 9. crops
! JSBACH: 10. pasture

CONTAINS

  ! --- config_dynveg ------------------------------------------------------------------
  !
  ! Reads in parameters from lctlib-file

  SUBROUTINE config_dynveg(dynveg_params, dynveg_options)

    USE mo_namelist,         ONLY: position_nml, POSITIONED
    USE mo_filename,         ONLY: find_next_free_unit
    USE mo_util_string,      ONLY: tolower
    USE mo_jsbach,           ONLY: debug, nml_unit
    USE mo_doctor,           ONLY: nout

    TYPE(dynveg_params_type),  INTENT(out) :: dynveg_params
    TYPE(dynveg_options_type), INTENT(out) :: dynveg_options
 
    !! --- parameters

    character(len=2),parameter  :: blank_set = " "//achar(9) !! the blank characters: BLANK and TAB

    !! --- locals

    INTEGER                :: UnitNo,read_status
    INTEGER                :: pos,pos_comment,length
    CHARACTER(len=256)     :: line
    CHARACTER(len=30)      :: key
    INTEGER                :: nlct              !! Number of landcover types found in dynveg landcover library file
    INTEGER, ALLOCATABLE   :: itmp(:)           !! temporary array used for input of logicals

    LOGICAL                :: exists_bclimit_min_cold_mmtemp, exists_bclimit_max_cold_mmtemp, exists_bclimit_max_warm_mmtemp
    LOGICAL                :: exists_bclimit_min_temprange, exists_bclimit_min_gdd, exists_gdd_base, exists_upper_tlim
    LOGICAL                :: exists_dynamic_PFT, exists_woody_PFT, exists_tau_pft


    !! Namelist Parameters
    CHARACTER(NF_MAX_NAME) :: lctlib_dynveg      !! lctlib for the dynamic vegetation
    CHARACTER(NF_MAX_NAME) :: fpc_file_name      !! file name for initial PFTs
    CHARACTER(NF_MAX_NAME) :: climbuf_file_name  !! file name for initial climate data
    LOGICAL                :: dynveg_all         !! competition between all PFTs, not just between woody PFTs
    LOGICAL                :: dynveg_feedback    !! activate climate feedback of the dynamic vegetation
    LOGICAL                :: init_running_means !! initialize running means of the climate buffer
    LOGICAL                :: read_climbuf       !! Read long term climate variables from file
    LOGICAL                :: read_fpc           !! Read fractional plant cover from file

    INCLUDE 'dynveg_ctl.inc'


    !! Read namelist dynveg_ctl

    IF (p_parallel_io) THEN

       ! define default values
       lctlib_dynveg = 'lctlib.def'
       read_fpc = .FALSE. 
       fpc_file_name = 'fpc.nc'
       dynveg_all = .FALSE.
       dynveg_feedback = .TRUE.
       init_running_means = .FALSE.
       read_climbuf = .FALSE.
       climbuf_file_name = 'climbuf.nc'
       
       CALL position_nml ('DYNVEG_CTL', status=read_status)
       SELECT CASE (read_status)
       CASE (POSITIONED)
          READ (nml_unit, dynveg_ctl)
          CALL message('config_dynveg', 'Namelist DYNVEG_CTL: ')
          WRITE(nout, dynveg_ctl)
       END SELECT
    ENDIF
    IF (p_parallel_io) THEN

       CALL message('config_dynveg','Land cover library file: '//TRIM(lctlib_dynveg))

       dynveg_options%read_fpc = read_fpc
       fpc_file%file_name = fpc_file_name
       IF (read_fpc) THEN
          CALL message('config_dynveg','Fractional plant cover read from external file')
          CALL message('config_dynveg','File containing initial values for fpc: '//TRIM(fpc_file_name))
       ENDIF

       dynveg_options%dynveg_all = dynveg_all
       IF (dynveg_all) THEN
          CALL message('config_dynveg', 'Competition between woody types and grass active')
       ENDIF

       dynveg_options%dynveg_feedback = dynveg_feedback
       IF (dynveg_feedback) THEN
          CALL message('config_dynveg', 'Dynamic vegetation feedback active')
       ENDIF

       dynveg_options%init_running_means = init_running_means
       IF (init_running_means) THEN
          CALL message('config_dynveg', 'Dynveg running means are initialised')
       ENDIF

       dynveg_options%read_climbuf = read_climbuf
       dynveg_options%climbuf_file_name = climbuf_file_name
       IF (read_climbuf) THEN
          CALL message('config_dynveg','Climate buffer read from external file')
          CALL message('config_dynveg','File containing initial values for climbuf: '//TRIM(climbuf_file_name))
       ENDIF
    ENDIF
    IF (p_parallel) THEN
       CALL p_bcast(dynveg_options%read_fpc, p_io)
       CALL p_bcast(dynveg_options%read_climbuf, p_io)
       CALL p_bcast(dynveg_options%dynveg_all, p_io)
       CALL p_bcast(dynveg_options%dynveg_feedback, p_io)
       CALL p_bcast(dynveg_options%init_running_means, p_io)
    ENDIF


    ! Read basic parameters from the lctlib file

    IF (p_parallel_io) THEN

       !! open parameter file
       
       UnitNo = find_next_free_unit(50,100)

       OPEN(unit=UnitNo, file=lctlib_dynveg, form='FORMATTED', status='OLD', iostat=read_status)
       IF (read_status /= 0) CALL finish('config_dynveg','Error when opening dynveg library file '//trim(lctlib_dynveg))

       ! --- find keyword NLCT in landcover library file
       ! ---- Routine adopted from LCTLIB Routine ------

       DO
          READ(UnitNo,'(A256)', IOSTAT=read_status) line
          IF (read_status /= 0) THEN
             CALL finish('config_dynveg',&
               'Read error, probably because keyword NLCT is missing in dynveg land cover library file '//TRIM(lctlib_dynveg))
          ENDIF

          ! Skip comment lines

          pos_comment = SCAN(line,"#")                 ! Start position of comment, 0 if none present
          IF (pos_comment > 1) THEN
             line = TRIM(ADJUSTL(line(1:pos_comment))) ! Disregard everything to the right of "#", ..
             ! .. adjust to the left and trim to the right
          ELSE IF (pos_comment == 1) THEN
             CYCLE                                     ! Disregard whole line
          ELSE
             line = TRIM(ADJUSTL(line))                ! Only disregard blanks to the left and right
          ENDIF
          length=LEN_TRIM(line)
          IF(length== 0) CYCLE                         ! Line is empty
          
          pos = SCAN(line,blank_set)                   ! Position of first blank character
          IF(pos == 0)   CALL finish('config_dynveg',"Wrong syntax in dynveg definitions")
          
          READ(line(1:pos-1),'(A)') key
          IF (debug) CALL message('config_dynveg - reading',TRIM(key))
          
          IF(tolower(TRIM(key)) == 'nlct') THEN
             READ(line(pos:length),*,IOSTAT=read_status) nlct
             IF (read_status /= 0) THEN
                CALL finish('config_dynveg',&
                     'Could not read number of landcover types (keyword: NLCT) from '//TRIM(lctlib_dynveg))
             ENDIF
             EXIT ! found number of landcover types in landcover library file --- continue after loop
          ENDIF
       END DO

    ENDIF

    IF (p_parallel) CALL p_bcast(nlct, p_io)
    dynveg_params%nlct = nlct

    !! --- Allocate memory for PFT-specific parameters

    ALLOCATE(dynveg_params%dynamic_PFT(1:nlct))
    exists_dynamic_PFT=.FALSE.
    ALLOCATE(dynveg_params%woody_PFT(1:nlct))
    exists_woody_PFT=.FALSE.
    ALLOCATE(dynveg_params%bclimit_min_cold_mmtemp(1:nlct))
    exists_bclimit_min_cold_mmtemp = .FALSE.
    ALLOCATE(dynveg_params%bclimit_max_cold_mmtemp(1:nlct))
    exists_bclimit_max_cold_mmtemp = .FALSE.
    ALLOCATE(dynveg_params%bclimit_max_warm_mmtemp(1:nlct))
    exists_bclimit_max_warm_mmtemp = .FALSE.
    ALLOCATE(dynveg_params%bclimit_min_temprange(1:nlct))
    exists_bclimit_min_temprange = .FALSE.
    ALLOCATE(dynveg_params%bclimit_min_gdd(1:nlct))
    exists_bclimit_min_gdd = .FALSE.
    ALLOCATE(dynveg_params%gdd_base(1:nlct))
    exists_gdd_base = .FALSE.
    ALLOCATE(dynveg_params%upper_tlim(1:nlct))
    exists_upper_tlim = .FALSE.
    ALLOCATE(dynveg_params%tau_pft(1:nlct))
    exists_tau_pft = .FALSE.


    !! --- Read the PFT specific parameters from the library file

    IF (p_parallel_io) THEN

       ALLOCATE(itmp(1:nlct))  !! Temporary memory

       REWIND(unit=UnitNo) ! Go back to beginning of dynamic landcover library file

       DO
          READ(UnitNo,'(A256)', IOSTAT=read_status) line
          IF (read_status /= 0) EXIT                   !! Finished reading

          ! Look for comment
          pos_comment = SCAN(line,"#")                 !! Start position of comment, 0 if none present
          IF (pos_comment > 1) THEN
             line = TRIM(ADJUSTL(line(1:pos_comment))) !! Disregard everything to the right of "#", ..
             ! .. adjust to the left and trim to the right
          ELSE IF (pos_comment == 1) THEN
             CYCLE                                     !! Disregard whole line
          ELSE
             line = TRIM(ADJUSTL(line))                !! Only disregard blanks to the left and right
          ENDIF
          length=LEN_TRIM(line)
          IF (length== 0) CYCLE                        !! Line is empty

          pos = SCAN(line,blank_set)                   !! Position of first blank character
          IF (pos == 0) CALL finish('config_dynveg',"Wrong syntax in dynveg lctlib-file "//trim(lctlib_dynveg))

          READ(line(1:pos-1),'(A)') key

          IF (tolower(TRIM(key)) .EQ. "nlct") CYCLE    !! nlct already read above
          IF (debug) CALL message('config_dynveg',TRIM(key))
!!$          SELECT CASE (tolower(TRIM(key)))

          IF (tolower(TRIM(key)) == 'dynamic_pft') THEN  !! Indicates those PFTs that shall take part in vegetation dynamics
             READ(line(pos:length),*) itmp(1:nlct)
             WHERE(itmp(:) == 0) 
                dynveg_params%dynamic_PFT(:) = .FALSE.
             ELSEWHERE
                dynveg_params%dynamic_PFT(:) = .TRUE.
             END WHERE
             exists_dynamic_PFT = .TRUE.

          ELSEIF (tolower(TRIM(key)) == 'woody_pft') THEN       !! Indicates woody PFTs
             READ(line(pos:length),*)  itmp(1:nlct)
             WHERE(itmp(:) == 0) 
                dynveg_params%woody_PFT(:) = .FALSE.
             ELSEWHERE
                dynveg_params%woody_PFT (:) = .TRUE.
             END WHERE
             exists_woody_PFT = .TRUE.

          ELSEIF (tolower(TRIM(key)) == 'bclimit_min_cold_mmtemp') THEN   !! PFT-specific minimum coldest monthly mean temperature
             READ(line(pos:length),*) dynveg_params%bclimit_min_cold_mmtemp(1:nlct)
             exists_bclimit_min_cold_mmtemp = .TRUE.

          ELSEIF (tolower(TRIM(key)) == 'bclimit_max_cold_mmtemp') THEN   !! PFT-specific maximum coldest monthly mean temperature
             READ(line(pos:length),*) dynveg_params%bclimit_max_cold_mmtemp(1:nlct)
             exists_bclimit_max_cold_mmtemp = .TRUE.

          ELSEIF (tolower(TRIM(key)) == 'bclimit_max_warm_mmtemp') THEN   !! PFT-specific maximum warmest monthly mean temperature
             READ(line(pos:length),*) dynveg_params%bclimit_max_warm_mmtemp(1:nlct)
             exists_bclimit_max_warm_mmtemp = .TRUE.

          ELSEIF (tolower(TRIM(key)) == 'bclimit_min_temprange') THEN     !! PFT-specific minimum temperature range between warmest
             READ(line(pos:length),*) dynveg_params%bclimit_min_temprange(1:nlct) !!      and coldest month (20-year average)
             exists_bclimit_min_temprange = .TRUE.

          ELSEIF (tolower(TRIM(key)) == 'bclimit_min_gdd') THEN           !! PFT-specific minimum growing degree days
             READ(line(pos:length),*) dynveg_params%bclimit_min_gdd(1:nlct)      !!       (at or above 5 deg C)
             exists_bclimit_min_gdd = .TRUE.

          ELSEIF (tolower(TRIM(key)) == 'gdd_base') THEN                  !! PFT-specific GDD base
             READ(line(pos:length),*) dynveg_params%gdd_base(1:nlct)
             exists_gdd_base = .TRUE.

          ELSEIF (tolower(TRIM(key)) == 'upper_tlim') THEN                !! PFT-specific upper temperature limit
             READ(line(pos:length),*) dynveg_params%upper_tlim(1:nlct)           !!       (used to calculate gdd_upper_tlim)
             exists_upper_tlim = .TRUE.

          ELSEIF (tolower(TRIM(key)) == 'tau_pft') THEN                   !! PFT-specific time_scale
             READ(line(pos:length),*) dynveg_params%tau_pft(1:nlct)
             exists_tau_pft = .TRUE.

          ENDIF
       END DO

       DEALLOCATE(itmp) !! release temporary memory

       !! Check whether all keywords have been found

       IF(.NOT. exists_dynamic_PFT) &
            CALL finish('config_dynveg','No data for DYNAMIC_PFT found in '//TRIM(lctlib_dynveg))
       IF(.NOT. exists_woody_PFT) &
            CALL finish('config_dynveg','No data for WOODY_PFT found in '//TRIM(lctlib_dynveg))
       IF(.NOT. exists_bclimit_min_cold_mmtemp) &
            CALL finish('config_dynveg','No data for BCLIMIT_MIN_COLD_MMTEMP found in '//TRIM(lctlib_dynveg))
       IF(.NOT. exists_bclimit_max_cold_mmtemp) &
            CALL finish('config_dynveg','No data for BCLIMIT_MAX_COLD_MMTEMP found in '//TRIM(lctlib_dynveg))
       IF(.NOT. exists_bclimit_max_warm_mmtemp) &
            CALL finish('config_dynveg','No data for BCLIMIT_MAX_WARM_MMTEMP found in '//TRIM(lctlib_dynveg))
       IF(.NOT. exists_bclimit_min_temprange ) &
            CALL finish('config_dynveg','No data for BCLIMIT_MIN_TEMPRANGE found in '//TRIM(lctlib_dynveg))
       IF(.NOT.exists_bclimit_min_gdd ) &
            CALL finish('config_dynveg','No data for BCLIMIT_MIN_GDD found in '//TRIM(lctlib_dynveg))
       IF(.NOT. exists_gdd_base) &
            CALL finish('config_dynveg','No data for GDD_BASE found in '//TRIM(lctlib_dynveg))
       IF(.NOT. exists_upper_tlim) &
            CALL finish('config_dynveg','No data for UPPER_TLIM found in '//TRIM(lctlib_dynveg))
       IF(.NOT. exists_tau_pft) &
            CALL finish('config_dynveg','No data for TAU_PFT found in '//TRIM(lctlib_dynveg))

    ENDIF

    IF (p_parallel) THEN
       CALL p_bcast(dynveg_params%dynamic_pft(:), p_io)
       CALL p_bcast(dynveg_params%woody_pft(:), p_io)
       CALL p_bcast(dynveg_params%bclimit_min_cold_mmtemp(:), p_io)
       CALL p_bcast(dynveg_params%bclimit_max_cold_mmtemp(:), p_io)
       CALL p_bcast(dynveg_params%bclimit_max_warm_mmtemp(:), p_io)
       CALL p_bcast(dynveg_params%bclimit_min_temprange(:), p_io)
       CALL p_bcast(dynveg_params%bclimit_min_gdd(:), p_io)
       CALL p_bcast(dynveg_params%gdd_base(:), p_io)
       CALL p_bcast(dynveg_params%upper_tlim(:), p_io)
       CALL p_bcast(dynveg_params%tau_pft(:), p_io)
    ENDIF

  END SUBROUTINE config_dynveg

  ! --- init_dynveg --------------------------------------------------------------------------------------------------------
  !
  ! Initialises this dynveg module. In particular the structure "dynveg" is initialised,
  ! this means, that the initialisation of dynveg is done here.
  
  SUBROUTINE init_dynveg(grid, domain, no_tiles, no_soil, lrestart, fileformat, stream)

    USE mo_linked_list,            ONLY: LAND, TILES
    USE mo_memory_base,            ONLY: new_stream,default_stream_setting, &
                                         add =>add_stream_element
    USE mo_netcdf,                 ONLY: max_dim_name, io_inq_varid, io_get_var_double
    USE mo_io,                     ONLY: io_open, io_close
    USE mo_transpose,              ONLY: scatter_gp
    USE mo_decomposition,          ONLY: global_decomposition
    USE mo_climbuf,                ONLY: init_climbuf
    USE mo_temp,                   ONLY: zreal2d, zreal3d, zzreal2d, zreal2d_ptr
    USE mo_jsbach,                 ONLY: missing_value

    TYPE(grid_type),   INTENT(in)        :: grid
    TYPE(domain_type), INTENT(in)        :: domain
    LOGICAL,           INTENT(in)        :: lrestart      ! true for restarted runs (throughout the whole run)
    INTEGER,           INTENT(in)        :: no_tiles      ! maximum number of living pfts in a grid cell
    INTEGER,           INTENT(in)        :: no_soil       ! number of soil levels
    INTEGER,           INTENT(in)        :: fileformat    ! output file format (grib/netcdf)
    TYPE(t_stream), POINTER, OPTIONAL    :: stream

    ! local variables

    INTEGER                     :: IO_file_id, IO_var_id, i
    INTEGER                     :: g_nland, l_nland
    INTEGER                     :: dim1p(1), dim1(1)
    INTEGER                     :: dim2p(2), dim2(2)
    CHARACTER(max_dim_name) :: dim1n(1), dim2n(2)
    
    !! Check initialization

    IF (module_initialized)  CALL finish("init_dynveg","Attempt to initialize twice.")

    nsoil=no_soil
    ntiles=no_tiles


    !! Read in namelist and lctlib-file for dynamical vegetation

    CALL config_dynveg(dynveg_params, dynveg_options)

    !! Stream definition
    
    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'dynveg', filetype=fileformat)
          ! Set default stream options
          CALL default_stream_setting(stream, table=180, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
       IO_dynveg => stream
    ELSE
       ! Add new stream
       CALL new_stream(IO_dynveg, 'dynveg', filetype=fileformat)
       ! Set default stream options
       CALL default_stream_setting(IO_dynveg, table=180, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF

    !! Initialisation of the climate buffer

    CALL init_climbuf(grid, domain, ntiles, lrestart, fileformat, nsoil, dynveg_options%read_climbuf, &
         dynveg_options%climbuf_file_name, IO_dynveg)

    g_nland = grid%nland
    l_nland = domain%nland

    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    dim2p = (/ l_nland, ntiles /)
    dim2  = (/ g_nland, ntiles /)
    dim2n(1) = 'landpoint'
    dim2n(2) = 'tiles'

    ! --- Define variables as stream elements
    CALL add(IO_dynveg,'burned_fpc'      ,dynveg%burned_fpc,  longname='burned area fraction',   units='', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=32, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_dynveg,'damaged_fpc'     ,dynveg%damaged_fpc, longname='damaged area fraction',  units='', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=33, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_dynveg,'pot_fpc'         ,dynveg%pot_fpc,     longname='potential plant cover fraction', units='', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=38, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_dynveg,'bio_exist'       ,dynveg%bio_exist,   longname='bio exist',              units='', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=39, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_dynveg,'act_fpc'         ,dynveg%act_fpc,     longname='fractional plant cover', units='', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=31, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_dynveg,'max_green_bio'   ,max_green_bio,      longname='maximum amount of green biomass in the year', &
             units='mol(C) m-2(canopy)', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=36, lmiss=.TRUE., missval=missing_value,   lpost=.false.)
    CALL add(IO_dynveg,'sum_green_bio_memory',sum_green_bio_memory,longname='vegetated fraction calculated from green biomass', &
             units='mol(C) m-2(canopy)', &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=40, lmiss=.TRUE., missval=missing_value,   lpost=.true.)
    CALL add(IO_dynveg,'desert_fpc'      ,dynveg%desert_fpc,  longname='desert fraction',        units='', &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=34, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_dynveg,'bare_fpc'        ,dynveg%bare_fpc,    longname='bare soil fraction',     units='', &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=35, lmiss=.TRUE., missval=missing_value)
    IF (dynveg_options%dynveg_feedback) THEN
       CALL add(IO_dynveg,'carbon_2_LeafLitterPool', dynveg%carbon_2_LeafLitterPool, lpost=.false., &
                longname='C transferred to leaf litter by vegetation dynamics', &
                units='mol(C)/m^2(grid box)', ldims=dim1p, gdims=dim1, dimnames=dim1n, &
                code=200, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_dynveg,'carbon_2_WoodLitterPool', dynveg%carbon_2_WoodLitterPool, lpost=.false., &
                longname='C transferred to wood litter by vegetation dynamics and fire', &
                units='mol(C)/m^2(grid box)', ldims=dim1p, gdims=dim1, dimnames=dim1n, &
                code=201, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_dynveg,'carbon_2_slowSoilPool_damage', dynveg%carbon_2_slowSoilPool_damage, lpost=.false., &
                longname='C transferred to slow soil pool by damage', &
                units='mol(C)/m^2(grid box)', ldims=dim1p, gdims=dim1, dimnames=dim1n, &
                code=202, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_dynveg,'carbon_2_slowSoilPool_fire', dynveg%carbon_2_slowSoilPool_fire, lpost=.false., &
                longname='C transferred to slow soil pool by fire', &
                units='mol(C)/m^2(grid box)', ldims=dim1p, gdims=dim1, dimnames=dim1n, &
                code=203, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_dynveg,'carbon_2_atmos', dynveg%carbon_2_atmos, &
                longname='CO2 emitted to the atmosphere by fire', &
                units='mol(C)/m^2(grid box)', ldims=dim1p, gdims=dim1, dimnames=dim1n, &
                code=204, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    ENDIF
    IF (dynveg_options%read_fpc) THEN

       ALLOCATE(init_act_fpc(l_nland,ntiles))
       ALLOCATE(init_bare_fpc(l_nland))
       ALLOCATE(init_desert_fpc(l_nland))
       ALLOCATE(zzreal2d(domain%ndim,domain%nblocks))

       IF (p_parallel_io) THEN
          fpc_file%opened = .FALSE.
          CALL IO_open(TRIM(fpc_file%file_name), fpc_file, IO_READ)
          IO_file_id = fpc_file%file_id

          ALLOCATE(zreal3d(grid%nlon,grid%nlat,ntiles))
          ALLOCATE(zreal2d(grid%nlon,grid%nlat))
       ENDIF

       !! --- Get initial values for act_fpc (and allocate further memory)

       IF (p_parallel_io) THEN 
          CALL IO_inq_varid(IO_file_id, 'act_fpc', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_gp(zreal2d_ptr, zzreal2d, global_decomposition)
          init_act_fpc(:,i) = PACK(zzreal2d, MASK=domain%mask)
       END DO

       !! --- Get bare_fpc

       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'bare_fpc', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal2d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       IF (p_parallel_io) zreal2d_ptr => zreal2d(:,:)
       CALL scatter_gp(zreal2d_ptr, zzreal2d, global_decomposition)
       init_bare_fpc(:) = PACK(zzreal2d, MASK=domain%mask)

       !! --- Get desert_fpc

       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'desert_fpc', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal2d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       IF (p_parallel_io) zreal2d_ptr => zreal2d(:,:)
       CALL scatter_gp(zreal2d_ptr, zzreal2d, global_decomposition)
       init_desert_fpc(:) = PACK(zzreal2d, MASK=domain%mask)

       !! --- Finish

       IF (p_parallel_io) THEN
          CALL IO_close(fpc_file)
          DEALLOCATE(zreal3d)
          DEALLOCATE(zreal2d)
       ENDIF
       DEALLOCATE(zzreal2d)

    ENDIF

    IF (lrestart) RETURN

    ! -- Initialisation at beginning of an experiment

    dynveg%bio_exist = 1._dp
    dynveg%burned_fpc = 0._dp
    dynveg%damaged_fpc = 0._dp
    dynveg%pot_fpc = 0._dp
    max_green_bio = 0._dp
    sum_green_bio_memory = 0._dp
    IF (dynveg_options%dynveg_feedback) THEN
       dynveg%carbon_2_LeafLitterPool = 0._dp
       dynveg%carbon_2_WoodLitterPool = 0._dp
       dynveg%carbon_2_slowSoilPool_damage = 0._dp
       dynveg%carbon_2_slowSoilPool_fire = 0._dp
       dynveg%carbon_2_atmos = 0._dp
    ENDIF

    CALL message('init_dynveg','fpc initialised')

  END SUBROUTINE init_dynveg

  ! --- update_dynveg------------------------------------------------------------------------------------------------------

  SUBROUTINE update_dynveg(lstart, lresume, nidx, kidx0, kidx1, &
                          temp_air, rel_hum_air, temp_surf, wind10, &
                          precip, rel_soil_moisture, & 
                          npp_rate, sla, lai, veg_fract_correction, &
                          glacier, is_present, read_cover_fract, &
                          cpool_green, cpool_reserve, &
                          cpool_woods, cpool_litter_leaf, &
                          cpool_litter_wood, cpool_fast, &
                          cpool_slow, rock_fract, &
                          cover_fract, veg_ratio_max, co2_emission)

! !DESCRIPTION
! 
! the module simulates vegetation dynamics based on averaged annual npp values
!
! Input:
!   new_year:  flag for the end of the year (.TRUE. at the end of the year).
!   npp_yr_sum     :  annual NPP (g C/m2/yr) for a given PFT in a given grid
!                     cell
!   wind_speed     :  extreme wind speed (m/s)
!   litter_moisture:  daily litter moisture (%) in a given grid cell for the 
!                     fire module 
!   cpool_litter    :  amount of above-ground litter [mol(C)/m^2(grid box)] for the fire 
!                     module
!   bioclimate_ann :  annually averaged climatology
!
! In-/output:
!   act_fpc        :  fractional projective cover (between 0 and 1) of each PFT with respect
!                     to the part of the grid box, which is vegetated and not excluded for agriculture 
!
! Output:
!   burned_fpc     :  burned area (between 0 and 1) for a given PFT in a given 
!                     grid cell
!   damaged_fpc    :  wind-damaged area (between 0 and 1) for a given PFT in a
!                     given grid cell
!   pot_fpc        :  --
!
!------------------------------------------------------------------------------
    USE mo_time_control,  ONLY: get_time_step, get_date_components, get_year_day,   &
                                start_date, current_date, previous_date, delta_time
    USE mo_climbuf,       ONLY: update_climbuf
    USE mo_cbal_cpools,   ONLY: relocate_carbon, relocate_carbon_desert, relocate_carbon_fire, relocate_carbon_damage
    USE mo_utils,         ONLY: average_tiles
    use mo_jsbach_constants, ONLY: molarMassCO2_kg

! input parameters

    LOGICAL,   INTENT(in)    :: lstart                 ! .TRUE. at the beginning of an experiment
    LOGICAL,   INTENT(in)    :: lresume                ! .TRUE. at the beginning of restarted runs
    INTEGER,   INTENT(in)    :: nidx                   ! Vector lenght
    INTEGER,   INTENT(in)    :: kidx0, kidx1           ! first and last index in the global array
    REAL(dp),  INTENT(in)    :: temp_air(:)            ! air temperature of the lowest atmosphere level
    REAL(dp),  INTENT(in)    :: rel_hum_air(:)         ! relative air humidity of the lowest atmosphere level
    REAL(dp),  INTENT(in)    :: temp_surf(:,:)         ! surface temperature
    REAL(dp),  INTENT(in)    :: wind10(:)              ! 10 m wind speed [m/s]
    REAL(dp),  INTENT(in)    :: precip(:)              ! precipitation as rain
    REAL(dp),  INTENT(in)    :: rel_soil_moisture(:,:,:)  ! relative soil moisture
    REAL(dp),  INTENT(in)    :: npp_rate(:,:)          ! Rate of NPP [mol(C)/m^2(canopy) s]
    REAL(dp),  INTENT(in)    :: sla(:,:)               ! specific leaf area
    REAL(dp),  INTENT(in)    :: lai(:,:)               ! leaf area index
    REAL(dp),  INTENT(in)    :: veg_fract_correction(:,:) ! corrects vegetated area 1-exp(-LAI_max/2)
    LOGICAL,   INTENT(in)    :: glacier(:)             ! glacier mask
    LOGICAL,   INTENT(in)    :: is_present(:,:)        ! determines, if tile is handled anyway
    LOGICAL,   INTENT(in)    :: read_cover_fract       ! read cover fractions from initial rather than restart file

! in- and output parameters
    REAL(dp),  INTENT(inout) :: cpool_green(:,:)       ! green biomass [mol(C)/m^2(canopy)]
    REAL(dp),  INTENT(inout) :: cpool_reserve(:,:)     ! sugars, starches etc. [mol(C)/m^2(canopy)]
    REAL(dp),  INTENT(inout) :: cpool_woods(:,:)       ! wood biomass [mol(C)/m^2(canopy)]
    REAL(dp),  INTENT(inout) :: cpool_litter_leaf(:,:) ! leaf litter biomass [mol(C)/m^2(canopy)]
    REAL(dp),  INTENT(inout) :: cpool_litter_wood(:,:) ! woody litter biomass [mol(C)/m^2(canopy)]
    REAL(dp),  INTENT(inout) :: cpool_fast(:,:)        ! biomass in the fast soil pool [mol(C)/m^2(canopy)]
    REAL(dp),  INTENT(inout) :: cpool_slow(:,:)        ! biomass in the slow soil pool [mol(C)/m^2(canopy)]
    REAL(dp),  INTENT(inout) :: rock_fract(:)          ! fraction of area not vegetated at all 
                                                       !   (lakes, rocks, streets, ...)
    REAL(dp),  INTENT(inout) :: cover_fract(:,:)       ! fractional cover of the tiles
    REAL(dp),  INTENT(inout) :: veg_ratio_max(:)       ! fractional cover of vegetated area

! output parameters
    REAL(dp),  INTENT(out) :: co2_emission(:)          ! co2-flux to the atmosphere due to fires [kg(CO2)/m^2(gridbox)s]

! local variables

    INTEGER    :: istep                               ! current time step (since initialization)
    INTEGER    :: year, year_at_prev_ts, start_year
    INTEGER    :: day, day_at_prev_ts
    INTEGER    :: time_steps_per_day                  ! number of time steps per day
    INTEGER    :: time_steps_in_last_year             ! number of time steps within ot previous year
    INTEGER    :: pft
    LOGICAL    :: new_year, new_day
    REAL(dp)   :: cpool_litter(nidx)                  ! biomass of litter [mol(C)/m^2(grid box)] 
    REAL(dp)   :: veg_ratio_max_old(nidx)          ! fractional cover of vegetated area of last year
    REAL(dp)   :: cover_fract_old(nidx,ntiles)     ! JSBACH cover fractions of last year
    REAL(dp)   :: act_fpc_old(nidx,ntiles)         ! fractional projective cover of last day

!-- initialise FPC

    IF (lstart .OR. lresume) THEN
       IF (dynveg_options%read_fpc) THEN
          CALL message('update_dynveg','Initial values of FPC read from file' )
          dynveg%act_fpc(kidx0:kidx1,:) = init_act_fpc(kidx0:kidx1,:)
          dynveg%bare_fpc(kidx0:kidx1) = init_bare_fpc(kidx0:kidx1)
          dynveg%desert_fpc(kidx0:kidx1) = init_desert_fpc(kidx0:kidx1)
       ELSE IF (.NOT. lresume .OR. read_cover_fract) THEN
!
!      -- calculate initial fpc from cover_fract of the jsbach initial data file
!
          CALL message('update_dynveg','FPC initialized' )
          dynveg%act_fpc(kidx0:kidx1,:) = 0._dp
          dynveg%bare_fpc(kidx0:kidx1) = 0._dp
          dynveg%desert_fpc(kidx0:kidx1) = 0._dp
          CALL scale_cover_fract(nidx, ntiles, glacier(:), cover_fract(:,:))
          CALL cover_fract_to_fpc(nidx, ntiles, cover_fract(:,:), veg_ratio_max(:),   &
               rock_fract(:), glacier(:), dynveg_params%dynamic_pft(:), &
               dynveg%act_fpc(kidx0:kidx1,1:ntiles), dynveg%bare_fpc(kidx0:kidx1),  &
               dynveg%desert_fpc(kidx0:kidx1))
       ENDIF

       WHERE (glacier(:))
          rock_fract(:) = 1._dp
       END WHERE
       DO pft=1,ntiles
          WHERE (rock_fract(:) >= 1._dp - EPSILON(1._dp))
             dynveg%act_fpc(kidx0:kidx1,pft) = 0._dp
          END WHERE
       END DO

!   -- make sure, sum of fpc plus bare fraction is one in each grid cell

       CALL scale_fpc(nidx, ntiles, glacier(:), dynveg_params%dynamic_pft(:), &
            dynveg%act_fpc(kidx0:kidx1,:), dynveg%bare_fpc(kidx0:kidx1))

       IF (dynveg_options%read_fpc .and. dynveg_options%dynveg_feedback) THEN
!
!   -- call fpc_to_cover_fract to get consistent veg_ratio_max
!
          CALL fpc_to_cover_fract(nidx, ntiles, dynveg%act_fpc(kidx0:kidx1,:), dynveg%bare_fpc(kidx0:kidx1),  &
                                  dynveg%desert_fpc(kidx0:kidx1), rock_fract(:), glacier(:),                  &
                                  dynveg_params%woody_pft(:), dynveg_params%dynamic_pft(:), cover_fract(:,:), &
                                  veg_ratio_max(:))
       ENDIF
    ENDIF

!-- update climate buffer

    CALL update_climbuf(nidx, kidx0, kidx1, &
         dynveg_params%gdd_base(:), dynveg_params%upper_tlim(1:ntiles), & 
         dynveg_options%init_running_means, &
         temp_air(:), temp_surf(:,:), wind10(:), &
         precip, rel_soil_moisture(:,1:nsoil,1:ntiles), & 
         npp_rate(:,1:ntiles), rel_hum_air(:))

!-- find out time step and whether new day or new year just started

    istep=get_time_step()

    call get_date_components(current_date, YEAR=year, DAY=day)
    call get_date_components(previous_date, YEAR=year_at_prev_ts, DAY=day_at_prev_ts)  
    call get_date_components(start_date, YEAR=start_year)  
    new_year = (year /= year_at_prev_ts)
    new_day = (day /= day_at_prev_ts)
    time_steps_per_day = 86400/INT(delta_time)
    time_steps_in_last_year = INT(get_year_day(previous_date)) * time_steps_per_day

! -- annual calculations

    IF (new_year) THEN
       IF (istep > time_steps_per_day) THEN

!      -- calculation of bioclimatic limits
          CALL bioclim_limits(climbuf%min_mmtemp20(kidx0:kidx1), climbuf%max_mmtemp20(kidx0:kidx1), &
                              climbuf%prev_year_gdd(kidx0:kidx1,:), dynveg%bio_exist(kidx0:kidx1,:))
    
!      -- potential FPC
          CALL potential_tree_fpc(nidx, ntiles, dynveg_options%dynveg_all, &
                                  dynveg_params%woody_pft(:), dynveg_params%dynamic_pft(:), &
                                  climbuf%ave_npp5(kidx0:kidx1,:), dynveg%bio_exist(kidx0:kidx1,:), &
                                  dynveg%act_fpc(kidx0:kidx1,:), dynveg%pot_fpc(kidx0:kidx1,:))

!      -- calculation of desert extend
          CALL desert_fraction(dynveg_options%init_running_means, nidx, ntiles, &
                               dynveg_params%woody_pft(:), dynveg_params%dynamic_pft(:), &
                               dynveg%act_fpc(kidx0:kidx1,:), dynveg%bare_fpc(kidx0:kidx1), sla(:,:), &
                               max_green_bio(kidx0:kidx1,:), sum_green_bio_memory(kidx0:kidx1), &
                               dynveg%desert_fpc(kidx0:kidx1))

!      -- reset yearly accumulated variables          
          max_green_bio(kidx0:kidx1,:) = 0._dp
       ENDIF

       IF (dynveg_options%dynveg_feedback) THEN

!      -- conversion of FPCs to JSBACH cover fractions and new vegetated fraction (veg_ratio_max)
          veg_ratio_max_old(:) = veg_ratio_max(:)
          cover_fract_old(:,:) = cover_fract(:,:)
          CALL fpc_to_cover_fract(nidx, ntiles, dynveg%act_fpc(kidx0:kidx1,:), dynveg%bare_fpc(kidx0:kidx1),  &
                                  dynveg%desert_fpc(kidx0:kidx1), rock_fract(:), glacier(:), &
                                  dynveg_params%woody_pft(:), dynveg_params%dynamic_pft(:), cover_fract(:,:), &
                                  veg_ratio_max(:))

!      -- relocate carbon accouting for the change in the cover fraction of each PFT
          CALL relocate_carbon(cover_fract_old(:,:), cover_fract(:,:), veg_fract_correction(:,:), fract_small,   &
                               Cpool_green(:,:), Cpool_woods(:,:), Cpool_reserve(:,:),                           &
                               Cpool_litter_leaf(:,:), Cpool_litter_wood(:,:), Cpool_fast(:,:), Cpool_slow(:,:))

!      -- scaling cpools to account for change in vegetated fraction in order to conserve the total mass of carbon
          CALL relocate_carbon_desert(nidx, ntiles, veg_ratio_max(:), veg_ratio_max_old(:), fract_small,           &
                                      lai(:,:), sla(:,:), Cpool_green(:,:), Cpool_woods(:,:), Cpool_reserve(:,:),  &
                                      Cpool_litter_leaf(:,:), Cpool_litter_wood(:,:), Cpool_fast(:,:), Cpool_slow(:,:))

       ENDIF

    ENDIF

! -- daily calculations

    IF (new_day) THEN

!      -- determine maximum of green biomass for the current year
       DO pft = 1,ntiles
          IF (dynveg_params%dynamic_pft(pft)) THEN
             max_green_bio(kidx0:kidx1,pft) = MAX(max_green_bio(kidx0:kidx1,pft),cpool_green(:,pft))
          ENDIF
       END DO

!      -- determine grid box average litter
       IF (dynveg_options%dynveg_feedback) THEN

          CALL average_tiles((cpool_litter_leaf(1:nidx,1:ntiles) + cpool_litter_wood(1:nidx,1:ntiles)) *                 &
                             veg_fract_correction(1:nidx,1:ntiles) * SPREAD(veg_ratio_max(1:nidx),NCOPIES=ntiles,DIM=2), &
                             is_present(1:nidx,1:ntiles), cover_fract(1:nidx,1:ntiles), cpool_litter(1:nidx))

          act_fpc_old(:,:) = dynveg%act_fpc(kidx0:kidx1,:)

       ELSE

          cpool_litter(:) = 0._dp

       ENDIF

!      -- calculation of FPC (competition of vegetation types and disturbances (fire, wind break))

       CALL fpc_daily(nidx, ntiles, dynveg_options%dynveg_all, dynveg_params%woody_pft(:), dynveg_params%dynamic_pft(:), &
                      dynveg_params%tau_pft(:), climbuf%ave_npp5(kidx0:kidx1,:), dynveg%bio_exist(kidx0:kidx1,:), &
                      dynveg%pot_fpc(kidx0:kidx1,:), climbuf%rel_hum_air(kidx0:kidx1),                &
                      climbuf%prev_day_max_wind10(kidx0:kidx1), climbuf%max_wind10(kidx0:kidx1),      &
                      cpool_litter(1:nidx),                                                           &
                      dynveg%act_fpc(kidx0:kidx1,:), dynveg%burned_fpc(kidx0:kidx1,:),                &
                      dynveg%damaged_fpc(kidx0:kidx1,:), dynveg%bare_fpc(kidx0:kidx1))

!      -- rescaling of actual fpc and bare fpc
 
       CALL scale_fpc(nidx, ntiles, glacier(:), dynveg_params%dynamic_pft(:), dynveg%act_fpc(kidx0:kidx1,:), &
            dynveg%bare_fpc(kidx0:kidx1))

       IF (dynveg_options%dynveg_feedback) THEN

          dynveg%carbon_2_WoodLitterPool(kidx0:kidx1) = 0._dp

!      -- changes in carbon pools caused by wind break
          CALL relocate_carbon_damage(nidx, ntiles, act_fpc_old(:,:), dynveg%damaged_fpc(kidx0:kidx1,:),   &
                                      cover_fract(:,:), fract_small,                                       &
                                      veg_fract_correction(:,:),  frac_wood_2_litter_wind,                 &
                                      Cpool_green(:,:), Cpool_woods(:,:), Cpool_reserve(:,:),              &
                                      Cpool_litter_leaf(:,:), Cpool_litter_wood(:,:), Cpool_slow(:,:),     &
                                      dynveg%carbon_2_slowSoilPool_damage(kidx0:kidx1),                    &
                                      dynveg%carbon_2_LeafLitterPool(kidx0:kidx1), &
                                      dynveg%carbon_2_WoodLitterPool(kidx0:kidx1))

!      -- changes in carbon pools caused by fire
          CALL relocate_carbon_fire(nidx, ntiles, act_fpc_old(:,:), dynveg%burned_fpc(kidx0:kidx1,:),               &
                                    cover_fract(:,:), fract_small,                                                  &
                                    veg_fract_correction(:,:), frac_wood_2_atmos, frac_wood_2_litter_fire,          &
                                    Cpool_green(:,:), Cpool_woods(:,:), Cpool_reserve(:,:),                         &
                                    Cpool_litter_leaf(:,:), Cpool_litter_wood(:,:), Cpool_slow(:,:),                &
                                    dynveg%carbon_2_WoodLitterPool(kidx0:kidx1),                                    &
                                    dynveg%carbon_2_slowSoilPool_fire(kidx0:kidx1), dynveg%carbon_2_atmos(kidx0:kidx1))

!      -- scale fluxes from vegetated area to whole grid box
          dynveg%carbon_2_slowSoilPool_damage(kidx0:kidx1) = dynveg%carbon_2_slowSoilPool_damage(kidx0:kidx1) * &
                                                             veg_ratio_max(1:nidx)
          dynveg%carbon_2_LeafLitterPool(kidx0:kidx1) = dynveg%carbon_2_LeafLitterPool(kidx0:kidx1) * veg_ratio_max(1:nidx)
          dynveg%carbon_2_WoodLitterPool(kidx0:kidx1) = dynveg%carbon_2_WoodLitterPool(kidx0:kidx1) * veg_ratio_max(1:nidx)
          dynveg%carbon_2_slowSoilPool_fire(kidx0:kidx1) = dynveg%carbon_2_slowSoilPool_fire(kidx0:kidx1) * veg_ratio_max(1:nidx)
          dynveg%carbon_2_atmos(kidx0:kidx1) = dynveg%carbon_2_atmos(kidx0:kidx1) * veg_ratio_max(1:nidx)
       ENDIF

    ENDIF

!-- Calculate CO2 flux to the atmosphere (updated once a day)

    IF (dynveg_options%dynveg_feedback) THEN
       CO2_emission(:) = dynveg%carbon_2_atmos(kidx0:kidx1) * molarMassCO2_kg/86400._dp
    ELSE
       CO2_emission(:) = 0._dp
    END IF

  END SUBROUTINE update_dynveg

!------------------------------------------------------------------------------
!
! !IROUTINE potential_tree_fpc
!
! !SUBROUTINE INTERFACE:

  SUBROUTINE potential_tree_fpc(nidx, ntiles, dynveg_all, woody_pft, dynamic_pft, npp_ave, &
                                bio_exist, act_fpc, pot_fpc)

! !DESCRIPTION:
!
! subroutine calculates potential FPC (in absence of disturbances) based on NPP for each PFT
!
! called once a year!
!------------------------------------------------------------------------------

! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: nidx            ! vector length
    INTEGER,  INTENT(in)  :: ntiles
    LOGICAL,  INTENT(in)  :: dynveg_all      ! flag to activate competition between woody types and grass
    LOGICAL,  INTENT(in)  :: woody_pft(:)
    LOGICAL,  INTENT(in)  :: dynamic_pft(:)
    REAL(dp), INTENT(in)  :: bio_exist(:,:)
    REAL(dp), INTENT(in)  :: npp_ave(:,:)    ! NPP averaged
    REAL(dp), INTENT(in)  :: act_fpc(:,:)

! !OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(inout) :: pot_fpc(:,:)

! !LOCAL VARIABLES:
!
    INTEGER           :: pft
    REAL(dp)              :: sum_npp(nidx)
      
!------------------------------------------------------------------------------
!
!-- summing up NPP
!
    sum_npp(:)=0._dp
    DO pft=1, ntiles
       IF (dynveg_all) THEN
          IF (dynamic_pft(pft)) THEN
             WHERE (npp_ave(:,pft) > sum_npp_min .AND. bio_exist(:,pft) > 0.5_dp)
                sum_npp(:) = sum_npp(:) + npp_ave(:,pft) ** npp_nonlinearity * MAX(act_fpc_min,act_fpc(:,pft))
             END WHERE
          ENDIF
       ELSE
          IF (woody_pft(pft) .AND. dynamic_pft(pft)) THEN
             WHERE (npp_ave(:,pft) > sum_npp_min .AND. bio_exist(:,pft) > 0.5_dp)
                sum_npp(:) = sum_npp(:) + npp_ave(:,pft) ** npp_nonlinearity * MAX(act_fpc_min,act_fpc(:,pft))
             END WHERE
          ENDIF
       ENDIF
    ENDDO
    DO pft=1, ntiles
       IF (dynveg_all) THEN
          IF (dynamic_pft(pft)) THEN
             WHERE (npp_ave(:,pft) > sum_npp_min .AND. sum_npp(:) > sum_npp_min .AND. bio_exist(:,pft) > 0.5_dp) 
                pot_fpc(:,pft) = npp_ave(:,pft) ** npp_nonlinearity / sum_npp(:) * MAX(act_fpc_min,act_fpc(:,pft))
             ELSEWHERE
                pot_fpc(:,pft) = 0._dp
             END WHERE
          ELSE
             pot_fpc(:,pft) = 0._dp
          ENDIF
       ELSE
          IF (woody_pft(pft) .AND. dynamic_pft(pft)) THEN
             WHERE (npp_ave(:,pft) > sum_npp_min .AND. sum_npp(:) > sum_npp_min .AND. bio_exist(:,pft) > 0.5_dp)
                pot_fpc(:,pft) = npp_ave(:,pft) ** npp_nonlinearity / sum_npp(:) * &
                                 tree_fpc_max * MAX(act_fpc_min,act_fpc(:,pft))
             ELSEWHERE
                pot_fpc(:,pft) = 0._dp
             END WHERE
          ELSE
             pot_fpc(:,pft) = 0._dp
          ENDIF
       ENDIF
    ENDDO

  END SUBROUTINE potential_tree_fpc

!------------------------------------------------------------------------------
!
! !IROUTINE desert_fraction
!
! !SUBROUTINE INTERFACE:

  SUBROUTINE desert_fraction(init_running_means, nidx, ntiles,       &
                             woody_pft, dynamic_pft,              &
                             act_fpc, bare_fpc, sla,              &
                             max_green_bio, sum_green_bio_memory, &
                             desert_fpc)
! !DESCRIPTION:
!
!------------------------------------------------------------------------------


! !INPUT PARAMETERS:
!
    LOGICAL,  INTENT(in)  :: init_running_means
    INTEGER,  INTENT(in)  :: nidx            ! vector length
    INTEGER,  INTENT(in)  :: ntiles
    LOGICAL,  INTENT(in)  :: woody_pft(:)
    LOGICAL,  INTENT(in)  :: dynamic_pft(:)
    REAL(dp), INTENT(in)  :: act_fpc(:,:)
    REAL(dp), INTENT(in)  :: bare_fpc(:)
    REAL(dp), INTENT(in)  :: sla(:,:)
    REAL(dp), INTENT(in)  :: max_green_bio(:,:)

! !OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(inout) :: sum_green_bio_memory(:)
    REAL(dp), INTENT(inout) :: desert_fpc(:)

! !LOCAL VARIABLES:
!
    INTEGER           :: pft
    REAL(dp)          :: sum_green_bio(nidx)
    REAL(dp)          :: sum_act_fpc(nidx)
    REAL(dp)          :: sum_grass_fpc(nidx)
    INTEGER           :: n_grass_pft
      
!------------------------------------------------------------------------------

!
! calculating desert fraction from green biomass
!
    sum_green_bio(:) = 0._dp
    sum_act_fpc(:) = 0._dp
    sum_grass_fpc(:) = 0._dp
    n_grass_pft = 0

    DO pft = 1,ntiles
       IF (.NOT. woody_pft(pft) .AND. dynamic_pft(pft)) THEN
          sum_grass_fpc(:) = sum_grass_fpc(:) + act_fpc(:,pft)
          n_grass_pft = n_grass_pft + 1
       END IF
    END DO

    DO pft = 1,ntiles
       IF (woody_pft(pft) .AND. dynamic_pft(pft)) THEN
          sum_green_bio(:) = sum_green_bio(:) + MAX(0._dp, act_fpc(:,pft) * (1._dp -             &
            exp(-desert_extend * ((max_green_bio(:,pft) * sla(:,pft)) ** desert_margin) / &
            (3._dp ** (desert_margin - 1._dp)))))
          sum_act_fpc(:) = sum_act_fpc(:) + act_fpc(:,pft)
       ELSEIF (.NOT. woody_pft(pft) .AND. dynamic_pft(pft)) THEN
          WHERE (sum_grass_fpc(:) > EPSILON(1._dp))
             sum_green_bio(:) = sum_green_bio(:) + &
               MAX(0._dp, act_fpc(:,pft) * (1._dp + bare_fpc(:) / sum_grass_fpc(:)) * (1._dp -      &
               exp(-desert_extend * ((max_green_bio(:,pft) * sla(:,pft)) ** desert_margin) / &
               (3._dp ** (desert_margin - 1._dp)))))
             sum_act_fpc(:) = sum_act_fpc(:) + act_fpc(:,pft) * (1._dp + bare_fpc(:) / sum_grass_fpc(:))
          ELSEWHERE
             sum_green_bio(:) = sum_green_bio(:) + &
               MAX(0._dp, (bare_fpc(:) / REAL(n_grass_pft)) * (1._dp -                              &
               exp(-desert_extend * ((max_green_bio(:,pft) * sla(:,pft)) ** desert_margin) / &
               (3._dp ** (desert_margin - 1._dp)))))
             sum_act_fpc(:) = sum_act_fpc(:) + (bare_fpc(:) / REAL(n_grass_pft))
          END WHERE
       ENDIF
    ENDDO

    WHERE (sum_act_fpc(:) > EPSILON(1._dp)) 
       sum_green_bio(:) = sum_green_bio(:) / sum_act_fpc(:) 
    END WHERE

    IF (init_running_means) THEN 
       sum_green_bio_memory(:) = sum_green_bio(:)
    ELSE
       sum_green_bio_memory(:) = (sum_green_bio_memory(:) * (tau_desert - 1._dp) + sum_green_bio(:)) / tau_desert
    END IF

    desert_fpc(:) = MIN(1._dp - 2._dp * EPSILON(1._dp),MAX(0._dp,1._dp - sum_green_bio_memory(:)))

  END SUBROUTINE desert_fraction

!------------------------------------------------------------------------------
!
! !IROUTINE: fpc_daily
!
! !SUBROUTINE INTERFACE:
!
  SUBROUTINE fpc_daily(nidx, ntiles, dynveg_all, woody_pft, dynamic_pft, tau_pft, &
                       npp_ave, bio_exist, pot_fpc, rel_hum_air, &
                       wind_prev_day, wind_clim, cpool_litter,     &
                       act_fpc, burned_fpc, damaged_fpc, bare_fpc)
!
! !DESCRIPTION:
!
! Calculation of daily FPC for trees, shrubs, and herbaceous plants
!
! called: each day
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in) :: nidx            ! vector length
    INTEGER,  INTENT(in) :: ntiles
    LOGICAL,  INTENT(in) :: dynveg_all      ! flag to activate competition between woody types and grass
    LOGICAL,  INTENT(in) :: woody_pft(:)
    LOGICAL,  INTENT(in) :: dynamic_pft(:)
    REAL(dp), INTENT(in) :: tau_pft(:)
    REAL(dp), INTENT(in) :: npp_ave(:,:)
    REAL(dp), INTENT(in) :: bio_exist(:,:)
    REAL(dp), INTENT(in) :: pot_fpc(:,:)
    REAL(dp), INTENT(in) :: rel_hum_air(:)
    REAL(dp), INTENT(in) :: wind_prev_day(:)
    REAL(dp), INTENT(in) :: wind_clim(:)
    REAL(dp), INTENT(in) :: cpool_litter(:)

! !INPUT/OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(inout) :: act_fpc(:,:)

! !OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(inout) :: burned_fpc(:,:)
    REAL(dp), INTENT(inout) :: damaged_fpc(:,:)
    REAL(dp), INTENT(inout) :: bare_fpc(:)

! !LOCAL VARIABLES:
!  
    INTEGER     :: pft
    REAL(dp)    :: woody_act_fpc(nidx)
    REAL(dp)    :: grass_act_fpc(nidx)
    REAL(dp)    :: total_act_fpc(nidx)
    REAL(dp)    :: delta_time_yr
    REAL(dp)    :: grass_sum(nidx)
    REAL(dp)    :: non_woody_fpc(nidx)
    REAL(dp)    :: woody_estab_fpc(nidx)

!------------------------------------------------------------------------------
! Initialisations
! ---------------
    woody_act_fpc(:)=0._dp
    grass_act_fpc(:)=0._dp
    delta_time_yr = 1._dp / 365._dp        ! time step in years
    grass_sum(:)=0._dp

    burned_fpc(:,:)=0._dp
    damaged_fpc(:,:)=0._dp

    woody_estab_fpc(:)=0._dp   
    non_woody_fpc(:)=bare_fpc(:)
    DO pft=1, ntiles
       IF (.NOT. dynamic_pft(pft)) act_fpc(:,pft) = 0._dp
       IF (.NOT. woody_pft(pft) .AND. dynamic_pft(pft)) non_woody_fpc(:) = non_woody_fpc(:) + act_fpc(:,pft)
    ENDDO

    DO pft=1, ntiles
       IF (woody_pft(pft) .AND. dynamic_pft(pft)) THEN
!
!-- area_burned and wind_damage
!
          burned_fpc(:,pft) = fire_minimum_woody * delta_time_yr * act_fpc(:,pft)
          WHERE (rel_hum_air(:) < fire_rel_hum_threshold .and. cpool_litter(:) > fire_litter_threshold)
             burned_fpc(:,pft) = (fire_minimum_woody + ((fire_rel_hum_threshold - rel_hum_air(:)) / &
                                 (fire_tau_woody * fire_rel_hum_threshold))) * delta_time_yr * act_fpc(:,pft)
          END WHERE
          WHERE (wind_prev_day(:) > wind_clim(:) * wind_threshold)
             damaged_fpc(:,pft) = MIN(MAX(0._dp,act_fpc(:,pft) - fract_small), &
                                  act_fpc(:,pft) * delta_time_yr * wind_damage_scale * wind_prev_day(:) ** 3._dp / wind_clim(:))
          END WHERE
! 
!-- dynamic equation for act_fpc, woody types
!
          IF (dynveg_all) THEN
             act_fpc(:,pft) = act_fpc(:,pft) &
                + accelerate_vegetation * ((pot_fpc(:,pft) * non_woody_fpc(:) - act_fpc(:,pft) / life_to_estab) &
                / tau_pft(pft) * delta_time_yr &
                - burned_fpc(:,pft) - damaged_fpc(:,pft))
             woody_estab_fpc(:) = woody_estab_fpc(:) + accelerate_vegetation * (pot_fpc(:,pft) * non_woody_fpc(:) &
                / tau_pft(pft) * delta_time_yr)
          ELSE
             act_fpc(:,pft) = act_fpc(:,pft) &
                + accelerate_vegetation * ((pot_fpc(:,pft) - act_fpc(:,pft)) / tau_pft(pft) * delta_time_yr &
                - burned_fpc(:,pft) - damaged_fpc(:,pft))
          ENDIF

          WHERE (act_fpc(:,pft) < 0._dp)
             act_fpc(:,pft) = 0._dp
          END WHERE

          woody_act_fpc(:) = woody_act_fpc(:) + act_fpc(:,pft)
       ENDIF
    ENDDO
!
!-- dynamic equation for act_fpc, grass
!
    IF (dynveg_all) THEN
       DO pft = 1,ntiles
          IF (.NOT. woody_pft(pft) .AND. dynamic_pft(pft)) THEN
             WHERE (non_woody_fpc(:) > EPSILON(1._dp))
                act_fpc(:,pft) = MAX(0._dp,act_fpc(:,pft) * (1._dp - woody_estab_fpc(:) / non_woody_fpc(:)))
             ENDWHERE
          ENDIF
       ENDDO
    ENDIF

    DO pft = 1,ntiles
       IF (.NOT. woody_pft(pft) .AND. dynamic_pft(pft)) THEN 
          WHERE (bio_exist(:,pft) > 0.5_dp .AND. npp_ave(:,pft) > EPSILON(1._dp))
             grass_sum(:) = grass_sum(:) + npp_ave(:,pft) * MAX(act_fpc_min,act_fpc(:,pft))
             grass_act_fpc(:) = grass_act_fpc(:) + MAX(act_fpc_min,act_fpc(:,pft))
          END WHERE
       ENDIF
    ENDDO

    DO pft = 1,ntiles
       IF (.NOT. woody_pft(pft) .AND. dynamic_pft(pft)) THEN 

!-- area_burned of grasses, pasture and crops
!
          burned_fpc(:,pft) = fire_minimum_grass * delta_time_yr * act_fpc(:,pft)
          WHERE (rel_hum_air(:) < fire_rel_hum_threshold .and. cpool_litter(:) > fire_litter_threshold)
            burned_fpc(:,pft) = (fire_minimum_grass + ((fire_rel_hum_threshold - rel_hum_air(:)) / &
                                (fire_tau_grass * fire_rel_hum_threshold))) * delta_time_yr * act_fpc(:,pft)
          END WHERE

          IF (dynveg_all) THEN
             WHERE (bio_exist(:,pft) > 0.5_dp .AND. npp_ave(:,pft) > EPSILON(1._dp))
                act_fpc(:,pft) = act_fpc(:,pft) &
                   + accelerate_vegetation * ((pot_fpc(:,pft) * bare_fpc(:)) &
                   / tau_pft(pft) * delta_time_yr &
                   - burned_fpc(:,pft))
             ELSEWHERE
                act_fpc(:,pft) = act_fpc(:,pft) &
                   + accelerate_vegetation * ((pot_fpc(:,pft) * bare_fpc(:) - act_fpc(:,pft) / life_to_estab) &
                   / tau_pft(pft) * delta_time_yr &
                   - burned_fpc(:,pft))
             END WHERE
          ELSE
             WHERE (bio_exist(:,pft) > 0.5_dp .AND. npp_ave(:,pft) > EPSILON(1._dp))  
                act_fpc(:,pft) = act_fpc(:,pft) - burned_fpc(:,pft)  &
                   + (1._dp-woody_act_fpc(:)-grass_act_fpc(:))/grass_sum(:)*npp_ave(:,pft)*MAX(act_fpc_min,act_fpc(:,pft)) &
                   / tau_pft(pft) * delta_time_yr
             ELSEWHERE
                act_fpc(:,pft) =act_fpc(:,pft) - burned_fpc(:,pft) &
                   -act_fpc(:,pft) / tau_pft(pft) * delta_time_yr
             END WHERE
          ENDIF
       ENDIF

! -sigma_fpc(pft) ? stochastic disturbance term to be used in the future

       WHERE (act_fpc(:,pft) < 0)
          act_fpc(:,pft) = 0._dp
       END WHERE

    ENDDO

!-- calculate total area fraction covered by dynamic vegetation

    total_act_fpc(:)=0._dp
    DO pft= 1,ntiles 
       IF (dynamic_pft(pft)) THEN 
          total_act_fpc(:) = total_act_fpc(:) + act_fpc(:,pft)
       ENDIF
    ENDDO

!-- calculation of bare soil fraction
    
    bare_fpc(:) = MAX(0._dp,1._dp - total_act_fpc(:))

  END SUBROUTINE fpc_daily

!------------------------------------------------------------------------------
!
! !IROUTINE: bioclim_limits
!
! !SUBROUTINE INTERFACE:
!
  SUBROUTINE bioclim_limits (min_mmtemp20, max_mmtemp20, prev_year_gdd, bio_exist)
!
! !DESCRIPTION:
!
! calculation of bioclimatic limits for each PFT based on updated
! LPJ lookup-table
!
! Limits based on 20-year running averages of coldest-month mean
! temperature and growing degree days (5 degree base, except larches
! (2 degrees)), and minimal temperature range required (larches).
! For SURVIVAL, coldest month temperature and GDD should be
! at least as high as PFT-specific limits.
! For REGENERATION, PFT must be able to survive AND coldest month
! temperature should be no higher than a PFT-specific limit.
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(in)  :: min_mmtemp20(:)  ! Minimum monthly mean temp. (20yr climatology)
    REAL(dp), INTENT(in)  :: max_mmtemp20(:)  ! Maximum monthly mean temp. (20yr climatology)
    REAL(dp), INTENT(in)  :: prev_year_gdd(:,:) ! GDD of previous year

! ! OUTPUT PARAMETERS: 
!
    REAL(dp), INTENT(inout) :: bio_exist(:,:)

! !LOCAL VARIABLES:
!
    INTEGER   :: pft
    REAL(dp)  :: tcmin         !PFT-specific minimum coldest-month temperature limit
    REAL(dp)  :: tcmax         !PFT-specific maximum coldest-month temperature limit
    REAL(dp)  :: gddmin        !PFT-specific minimum GDD
    REAL(dp)  :: twmax         !PFT-specific maximum warmest-month temperature limit
    REAL(dp)  :: min_temprange !PFT-specific minimum difference of 20-year average warmest 
                                  !                    minus coldest month temperature
!------------------------------------------------------------------------------
    DO pft=1,ntiles
       IF (dynveg_params%dynamic_pft(pft)) THEN

          tcmin=dynveg_params%bclimit_min_cold_mmtemp(pft)
          tcmax=dynveg_params%bclimit_max_cold_mmtemp(pft)
          twmax=dynveg_params%bclimit_max_warm_mmtemp(pft)
          min_temprange=dynveg_params%bclimit_min_temprange(pft)
          gddmin=dynveg_params%bclimit_min_gdd(pft)

          bio_exist(:,pft) = 1._dp
          WHERE (min_mmtemp20(:) < tcmin .OR. min_mmtemp20(:) > tcmax   &
               .OR. prev_year_gdd(:,pft) < gddmin .OR. max_mmtemp20(:) > twmax   &
               .OR. (max_mmtemp20(:) - min_mmtemp20) < min_temprange)
             bio_exist(:,pft) = 0._dp
          END WHERE

       ENDIF
    ENDDO
      
  END SUBROUTINE bioclim_limits



  SUBROUTINE fpc_to_cover_fract (nidx, ntiles, act_fpc, bare_fpc, desert_fpc, rock_fract, &
                                 glacier, woody_pft, dynamic_pft, cover_fract, veg_ratio_max)
!
! !DESCRIPTION:
!
! convert fractional plant cover calculated within dynveg to jsbach
! cover fractions
!
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: nidx
    INTEGER,  INTENT(in)  :: ntiles
    REAL(dp), INTENT(in)  :: act_fpc(:,:)  ! actual fpc in dynveg
    REAL(dp), INTENT(in)  :: bare_fpc(:)   ! dynveg bare soil fraction
    REAL(dp), INTENT(in)  :: desert_fpc(:) ! dynveg desert fraction
    REAL(dp), INTENT(in)  :: rock_fract(:) ! area fraction inappropriate for
                                           !    vegetation (rocks, lakes, ...)
    LOGICAL,  INTENT(in)  :: glacier(:)    ! logical glacier mask
    LOGICAL,  INTENT(in)  :: woody_pft(:)  ! PFTs that build wood
    LOGICAL,  INTENT(in)  :: dynamic_pft(:)! PFTs taking part in competition
!
! !IN- and OUTPUT PARAMETERS: 
!
    REAL(dp), INTENT(inout) :: cover_fract(:,:)  ! cover fraction within jsbach
!
! !OUTPUT PARAMETERS:
! 
    REAL(dp), INTENT(out) :: veg_ratio_max(:) ! maximum vegetated area
                                                !    fraction in jsbach

! !LOCAL VARIABLES:
!
    INTEGER   :: pft, i
    INTEGER   :: n_grass_pft             ! number of grass PFT taking part in the competition
    REAL(dp)  :: sum_grass_fpc(nidx)     ! part of the vegetated area covered by grass
    REAL(dp)  :: excluded_fract(nidx)    ! area fraction not considered in dynveg
!------------------------------------------------------------------------------

!-- Find out fraction not regarded in dynveg and number of dynamic PFTs

    n_grass_pft = 0
    excluded_fract(:) = 0._dp
    sum_grass_fpc(:) = 0._dp

    DO pft=1,ntiles
       IF ( .NOT. dynamic_pft(pft)) excluded_fract(:) = excluded_fract(:) + cover_fract(:,pft)
       IF ( .NOT. woody_pft(pft) .AND. dynamic_pft(pft)) THEN
          n_grass_pft = n_grass_pft + 1
          sum_grass_fpc(:) = sum_grass_fpc(:) + act_fpc(:,pft)
       ENDIF
    ENDDO

!-- calculate new veg_ratio_max

    DO i = 1,nidx
       IF (.NOT. glacier(i)) THEN
          veg_ratio_max(i) = MAX(fract_small, 1._dp - (desert_fpc(i) + rock_fract(i)))
       ELSE
          veg_ratio_max(i) = 0._dp
       ENDIF
    END DO

!-- calculate new cover_fractions

    DO pft=1,ntiles

       IF (woody_pft(pft) .AND. dynamic_pft(pft)) THEN
          DO i = 1,nidx
             IF (.NOT. glacier(i)) THEN
                cover_fract(i,pft) = act_fpc(i,pft) * (1._dp - excluded_fract(i))
             ENDIF
          END DO
       ELSEIF ( .NOT. woody_pft(pft) .AND. dynamic_pft(pft)) THEN
          DO i = 1,nidx
             IF (.NOT. glacier(i)) THEN
                cover_fract(i,pft) = act_fpc(i,pft) * (1._dp + bare_fpc(i) / sum_grass_fpc(i)) * &
                                     (1._dp - excluded_fract(i))
             ENDIF
          END DO
       ENDIF

    END DO

    CALL scale_cover_fract (nidx, ntiles, glacier(:), cover_fract(:,:))

    IF (ANY(SUM(cover_fract(:,:),DIM=2) > 1._dp + ntiles*EPSILON(1._dp)) .OR. &
        ANY(SUM(cover_fract(:,:),DIM=2) < 1._dp - ntiles*EPSILON(1._dp))) THEN
       WRITE (message_text,*) 'sum of cover_fract /= 1: ', &
            MINVAL(SUM(cover_fract(:,:),DIM=2)), MAXVAL(SUM(cover_fract(:,:),DIM=2)),  MAXLOC(SUM(cover_fract(:,:),DIM=2))
       CALL finish ('fpc_to_cover_fract', message_text)
    END IF
    DO pft = 1,ntiles
       IF (dynamic_pft(pft)) THEN
          IF (ANY(.NOT. glacier(:) .AND. cover_fract(:,pft) < fract_small)) THEN
             WRITE (message_text,*) 'cover_fract too small: ', MINVAL(cover_fract(:,pft))
             CALL finish ('fpc_to_cover_fract', message_text)
          END IF
       END IF
    END DO

  END SUBROUTINE fpc_to_cover_fract
!------------------------------------------------------------------------------

  SUBROUTINE cover_fract_to_fpc (nidx, ntiles, cover_fract, veg_ratio_max, &
                                 rock_fract, glacier, dynamic_pft, act_fpc, &
                                 bare_fpc, desert_fpc)
!
! !DESCRIPTION:
!
! convert jsbach cover fractions (cover_fract) to dynveg FPCs
!
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: nidx
    INTEGER,  INTENT(in)  :: ntiles
    REAL(dp), INTENT(in)  :: cover_fract(:,:) ! cover fraction within jsbach
    REAL(dp), INTENT(in)  :: veg_ratio_max(:) ! maximum vegetated area
                                              !    fraction in jsbach
    REAL(dp), INTENT(in)  :: rock_fract(:) ! area fraction inappropriate for
                                           !    vegetation (rocks, lakes, ...)
    LOGICAL,  INTENT(in)  :: glacier(:)    ! logical glacier mask
    LOGICAL,  INTENT(in)  :: dynamic_pft(:)! PFTs taking part in competition
!
! !IN- and OUTPUT PARAMETERS:
! 
    REAL(dp), INTENT(inout) :: act_fpc(:,:)  ! actual fpc in dynveg
    REAL(dp), INTENT(inout) :: bare_fpc(:)   ! fraction of bare ground
    REAL(dp), INTENT(inout) :: desert_fpc(:) ! dynveg desert fraction

! !LOCAL VARIABLES:
!
    INTEGER   :: pft
    REAL(dp)  :: effective_cover_fract(nidx,ntiles) ! actual vegetation cover fraction
    REAL(dp)  :: excluded_fract(nidx)    ! area fraction not considered in dynveg
    REAL(dp)  :: desert_fract(nidx)      ! fraction of desert in whole grid box
!------------------------------------------------------------------------------


!-- Vegetation cover fractions with respect to the whole grid cell

    DO pft=1,ntiles
       effective_cover_fract(:,pft) = cover_fract(:,pft) * veg_ratio_max(:)
    ENDDO

!-- Find out area fraction not regarded in dynveg

    excluded_fract(:) = rock_fract(:)
    DO pft=1,ntiles
       IF ( .NOT. dynamic_pft(pft)) THEN
          excluded_fract(:) = excluded_fract(:) + effective_cover_fract(:,pft)
       ENDIF
    ENDDO

!-- Calculate desert fraction

    desert_fract(:) = (1._dp - veg_ratio_max(:)) - rock_fract(:)

!-- calculate FPCs

    WHERE (glacier(:) .OR. excluded_fract(:) >= 1._dp - EPSILON(1._dp) &
         .OR. desert_fpc(:) >= 1._dp - EPSILON(1._dp))
       bare_fpc(:) = 1._dp
       desert_fpc(:) = 1._dp
    ELSEWHERE
       bare_fpc(:) = 0._dp
       desert_fpc(:) = desert_fract(:) / (1._dp - excluded_fract(:))
    END WHERE

    DO pft=1,ntiles
       IF (dynamic_pft(pft)) THEN
          WHERE (glacier(:) .OR. excluded_fract(:) >= 1._dp - EPSILON(1._dp) &
               .OR. desert_fpc(:) >= 1._dp - EPSILON(1._dp))
             act_fpc(:,pft) = 0._dp
          ELSEWHERE
             act_fpc(:,pft) = effective_cover_fract(:,pft) / &
                  (1._dp - desert_fpc(:)) / (1._dp - excluded_fract(:))   
          END WHERE
       END IF
    END DO

!-- Rescaling to assure that the sum of act_fpc + bare soil is exactly one
    CALL scale_fpc(nidx, ntiles, glacier(:), dynamic_pft(:), act_fpc(:,:), bare_fpc(:))

    IF (ANY(SUM(act_fpc(:,:),DIM=2)+bare_fpc(:) > 1._dp+ntiles*EPSILON(1._dp)) .OR. &
        ANY(SUM(act_fpc(:,:),DIM=2)+bare_fpc(:) < 1._dp-ntiles*EPSILON(1._dp))) THEN
       WRITE (message_text,*) 'sum of act_fpc and bare_fpc /= 1: ', &
            MINVAL(SUM(act_fpc(:,:),DIM=2)+bare_fpc(:))
       CALL finish ('cover_fract_to_fpc',message_text)
    END IF
    DO pft = 1,ntiles
       IF (ANY(.NOT. glacier(:) .AND. act_fpc(:,pft) < fract_small)) THEN
          WRITE (message_text,*) 'act_fpc too small: ', MINVAL(act_fpc(:,pft))
          CALL finish ('cover_fract_to_fpc',message_text)
       END IF
    END DO

  END SUBROUTINE cover_fract_to_fpc
!------------------------------------------------------------------------------

  SUBROUTINE scale_fpc (nidx, ntiles, glacier, dynamic_pft, act_fpc, bare_fpc)
!
! !DESCRIPTION:
!
! Rescaling of fractional plant coverage to assure that
!  - the sum of act_fpc and bare soil is exactly one
!  - all tiles have a minimum vegetated fraction
!
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: nidx           ! vector length
    INTEGER,  INTENT(in)  :: ntiles         ! number of tiles
    LOGICAL,  INTENT(in)  :: glacier(:)     ! logical glacier mask
    LOGICAL,  INTENT(in)  :: dynamic_pft(:) ! PFTs taking part in competition
!
! !IN- and OUTPUT PARAMETERS:
! 
    REAL(dp), INTENT(inout) :: act_fpc(:,:) ! actual fpc in dynveg
    REAL(dp), INTENT(inout) :: bare_fpc(:)  ! fraction of bare ground

! !LOCAL VARIABLES:
!
    INTEGER   :: pft
    INTEGER   :: nsparce(nidx)  ! number of PFTs with less then fract_small vegetation
    REAL(dp)  :: sum_fpc(nidx)  ! sum of the different cover fractions
    REAL(dp)  :: excess(nidx)   ! extra fraction, that needs to be redistibuted
!------------------------------------------------------------------------------

    sum_fpc(:) = bare_fpc(:)
    nsparce(:) = 0
    excess(:) = 0._dp
    DO pft = 1,ntiles
       IF (dynamic_pft(pft)) THEN
          WHERE (act_fpc(:,pft) > fract_small)
             sum_fpc(:) = sum_fpc(:) + act_fpc(:,pft)
          ELSEWHERE
             nsparce(:) = nsparce(:) + 1
             act_fpc(:,pft) = fract_small
          END WHERE
       ELSE
          act_fpc(:,pft) = 0._dp
       END IF
    END DO
    sum_fpc(:) = sum_fpc(:) + REAL(nsparce(:))*fract_small

    WHERE (sum_fpc <= 1._dp - EPSILON(1._dp))
       bare_fpc(:) = bare_fpc(:) + (1._dp - sum_fpc(:))
    ELSEWHERE (sum_fpc >= 1._dp + EPSILON(1._dp))
       bare_fpc(:) = bare_fpc(:) / sum_fpc(:)
    END WHERE

    DO pft = 1,ntiles
       IF (dynamic_pft(pft)) THEN
          WHERE (sum_fpc >= 1._dp + EPSILON(1._dp) .AND. act_fpc(:,pft) > fract_small*sum_fpc(:))
             act_fpc(:,pft) = act_fpc(:,pft) / sum_fpc(:)
          ELSEWHERE (sum_fpc >= 1._dp + EPSILON(1._dp) .AND. act_fpc(:,pft) > fract_small)
             excess(:) = fract_small - act_fpc(:,pft)/sum_fpc(:)
             act_fpc(:,pft) = fract_small
          ELSEWHERE (sum_fpc >= 1._dp + EPSILON(1._dp))
             act_fpc(:,pft) = fract_small
          END WHERE
       END IF
    END DO
    DO pft = 1,ntiles
       IF (dynamic_pft(pft)) THEN
          WHERE (excess(:) > EPSILON(1._dp) .AND. act_fpc(:,pft) > (2*fract_small - act_fpc(:,pft)/sum_fpc(:)))
             act_fpc(:,pft) = act_fpc(:,pft) - (fract_small - act_fpc(:,pft)/sum_fpc(:))
             excess(:) = 0._dp
          END WHERE
       END IF
    END DO


    WHERE (glacier(:))
       bare_fpc(:) = 1._dp
    END WHERE
    DO pft = 1,ntiles
       WHERE (glacier(:))
          act_fpc(:,pft) = 0._dp
       END WHERE
    END DO

  END SUBROUTINE scale_fpc

END MODULE mo_dynveg
