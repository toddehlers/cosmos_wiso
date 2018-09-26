module mo_cbal_bethy
!
! Computation of the net primary productivity (NPP) as in BETHY. 
!

  USE mo_jsbach_grid,      ONLY: grid_type, domain_type, kstart, kend, nidx
  USE mo_kind,             ONLY: dp
  USE mo_exception,        ONLY: finish
  USE mo_linked_list,      ONLY: t_stream
  USE mo_netCDF,           ONLY: FILE_INFO
  USE mo_io,               ONLY: IO_open, IO_READ, IO_close
  USE mo_mpi,              ONLY: p_parallel_io, p_io, p_pe, p_parallel, p_bcast
  USE mo_jsbach_constants, ONLY: molarMassCO2_kg
  USE mo_time_control,     ONLY: delta_time    !! time step in seconds

  IMPLICIT NONE

  ! === BEGIN OF PUBLIC PART =======================================================================================================


  TYPE cbalance_type  ! contains the state variables of the cbalance module (dims: domain%nland x vegetation%ntiles)
     INTEGER          :: ntiles
     REAL(dp),POINTER :: Cpool_green(:,:)    !! C-pool for leaves, fine roots, vegetative organs and other green (living) parts 
                                             !! .. of vegetation [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: Cpool_reserve(:,:)  !! C-pool for carbohydrate reserve (sugars, starches) that allows plants to survive 
                                             !! .. bad times[mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: Cpool_woods(:,:)    !! C-pool for stems, thick roots and other (dead) structural material of living 
                                             !! .. plants [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: Cpool_litter_leaf(:,:) !! C-pool for litter originating from leaves [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: Cpool_litter_wood(:,:) !! C-pool for litter originating from woody parts of the plants [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: Cpool_fast(:,:)     !! C-pool for below ground organic material that quickly decomposes partly[mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: Cpool_slow(:,:)     !! C-pool for below ground organic material (coming from the fast pool) [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: NPP_Rate(:,:)       !! The instantaneous NPP rate [mol(C)/m^2(canopy) s]
     REAL(dp),POINTER :: NPP_Rate_acc(:,:)   !! averaged NPP rate [mol(C)/m^2(canopy) s]
     REAL(dp),POINTER :: LAI_sum(:,:)        !! used to accumulate LAI over a day
     REAL(dp),POINTER :: NPP_sum(:,:)        !! used to accumulated NPP-Rate over a day 
     REAL(dp),POINTER :: GPP_sum(:,:)        !! used to accumulated GPP-Rate over a day     
     REAL(dp),POINTER :: topSoilTemp_sum(:)  !! used to accumulated upper layer soil temperature over a day
     REAL(dp),POINTER :: alpha_sum(:,:)      !! used to accumulated water stress coefficient alpha over a day
     
     REAL(dp),POINTER :: soil_respiration(:,:)    !! mean daily rate of heterotrophic (soil) respiration [mol(CO2)/m^2(ground)] 
     REAL(dp),POINTER :: NPP_flux_correction(:,:) !! Daily updated flux correction from yesterdays carbon balance [mol(CO2)/m^2(canopy) s]
  END TYPE cbalance_type

  TYPE cbalance_diag_type  !! contains diagnostic variables  (1.dim: domain%nland, 2.dim: vegetation%ntiles)
     REAL(dp),POINTER :: boxC_green(:,:)              !! As Cpool_green() but in kg and relative to grid box area [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxC_reserve(:,:)            !! As Cpool_reserve() but in kg and relative to grid box area [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxC_woods(:,:)              !! As Cpool_woods() but in kg and relative to grid box area [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxC_litter_leaf(:,:)        !! As Cpool_litter_leaves() but in kg and relative to grid box area [kg(C)/m^2(grid box)]
     REAL(dp),POINTER :: boxC_litter_wood(:,:)        !! As Cpool_litter_wood() but in kg and relative to grid box area [kg(C)/m^2(grid box)]
     REAL(dp),POINTER :: boxC_fast(:,:)               !! As Cpool_fast() but in kg and relative to grid box area [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxC_slow(:,:)               !! As Cpool_slow() but in kg and relative to grid box area [kg(C)/m^2(grid box)     
     REAL(dp),POINTER :: LAI_previousDayMean(:,:)     !! mean value of LAI the day before yesterday
     REAL(dp),POINTER :: LAI_yDayMean(:,:)            !! mean value of LAI yesterday (from  LAI_sum())
     REAL(dp),POINTER :: NPP_yDayMean(:,:)            !! mean value of NPP-Rate yesterday (from NPP_sum()) [mol(CO2)/(m^2(canopy) s)]
     REAL(dp),POINTER :: GPP_yDayMean(:,:)            !! mean value of GPP-Rate yesterday (from GPP_sum()) [mol(CO2)/(m^2(canopy) s)]
     REAL(dp),POINTER :: topSoilTemp_yDayMean(:)      !! mean value of upper layer soil temperature yesterday (from topSoilTemp_sum()) [K]
     REAL(dp),POINTER :: alpha_yDayMean(:,:)          !! mean value of water stress coefficient alpha yesterday (from alpha_sum())

     REAL(dp),POINTER :: box_soil_respiration(:,:)    !! Daily updated soil respiration rate [mol(CO2)/m^2(grid box) s]
     REAL(dp),POINTER :: box_NPP_yDayMean(:,:)        !! Daily updated mean NPP rate [mol(CO2)/m^2(grid box) s]
     REAL(dp),POINTER :: box_NPP_flux_correction(:,:) !! Daily updated flux correction for NPP [mol(CO2)/m^2(grid box) s]
     REAL(dp),POINTER :: box_GPP_yDayMean(:,:)        !! Daily updated mean GPP rate [mol(CO2)/m^2(grid box) s] 
     
     REAL(dp),POINTER :: litter_flux(:,:)             !! Total litter carbon entering the soil pools from green and wood pool plus
                                                      !! ... excess NPP [mol(C)/m^2(canopy) s]
     REAL(dp),POINTER :: box_litter_flux(:,:)         !! Same as litter_flux, but relative to grid box area [mol(C)/m^2(grid box) s]
     REAL(dp),POINTER :: box_Cpools_total(:)          !! Sum of all carbon pools [mol(C)/m^2(grid box)]

  END TYPE cbalance_diag_type

  PUBLIC :: cbalance_type
  PUBLIC :: cbalance_diag_type
  PUBLIC :: update_cbalance_bethy      ! This subroutine computes NPP
  PUBLIC :: init_cbalance_bethy        ! This subroutine initializes this module
  PUBLIC :: NPP_rate_bethy

  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ===============================================================================

  PRIVATE 

  !! parameters

  INTEGER, PARAMETER              :: sec_per_day         = 86400  !! seconds per day

  !! module variables

  LOGICAL, SAVE                   :: module_initialized = .FALSE. !! Signifies whether module has been initialized
  INTEGER, SAVE                   :: ntiles
  INTEGER, SAVE                   :: time_steps_per_day           !! Number of time steps per day

  TYPE(cbalance_diag_type), SAVE  :: cbalance_diag

  TYPE(t_stream), POINTER, SAVE   :: IO_cbalance  ! Memory stream for cbalance model state
  TYPE(t_stream), POINTER, SAVE   :: IO_diag      ! Memory stream for cbalance diagnostic output

  LOGICAL                         :: read_cpools = .FALSE.
  TYPE(FILE_INFO),   SAVE         :: cpool_file   ! Input file for C pools
  REAL(dp), POINTER, SAVE         :: init_Cpool_green(:,:),init_Cpool_woods(:,:),init_Cpool_reserve(:,:)
  REAL(dp), POINTER, SAVE         :: init_Cpool_litter_leaf(:,:),init_Cpool_litter_wood(:,:)
  REAL(dp), POINTER, SAVE         :: init_Cpool_fast(:,:),init_Cpool_slow(:,:)

CONTAINS

  ! --- init_cbalance_bethy() -----------------------------------------------------------------------------------------------------
  !
  ! Initializes this cbalance module. In particular the structure "cbalance" is initialized. 
  !
  SUBROUTINE init_cbalance_bethy(grid, domain, cbalance, fileformat, useDynveg, diag_stream, stream )

    USE mo_linked_list,            ONLY: LAND, TILES, NETCDF
    USE mo_grib,                   ONLY: land_table
    USE mo_memory_base,            ONLY: new_stream,default_stream_setting, &
                                         add =>add_stream_element
    USE mo_decomposition,          ONLY: global_decomposition
    USE mo_transpose,              ONLY: scatter_gp
    USE mo_netcdf,                 ONLY: NF_MAX_NAME,MAX_DIM_NAME,io_inq_varid,io_get_var_double
    USE mo_exception,              ONLY: message
    USE mo_time_event,             ONLY: io_time_event
    USE mo_cbal_cpools,            ONLY: printCpoolParameters
    USE mo_jsbach,                 ONLY: missing_value, nml_unit
    USE mo_temp                         ! Provides temporary arrays
    USE mo_namelist,               ONLY: position_nml, POSITIONED
    USE mo_doctor,                 ONLY: nout

    TYPE(grid_type),    INTENT(in)     :: grid
    TYPE(domain_type),  INTENT(in)     :: domain
    TYPE(cbalance_type),intent(inout)  :: cbalance     !! Pointer to cbalance structure (allocated in the calling routine)
                                                       !! .. that shall be initialized
    INTEGER,            INTENT(in)     :: fileformat   !! output file format (grib/netcdf)
    LOGICAL,            INTENT(in)     :: useDynveg    !! Experiment with dynamic vegetation 
    TYPE(t_stream), POINTER            :: diag_stream
    TYPE(t_stream), POINTER, OPTIONAL  :: stream
    
    ! local variables

    CHARACTER(NF_MAX_NAME) :: cpool_file_name
    INTEGER :: IO_file_id, IO_var_id

    INTEGER :: g_nland, l_nland
    INTEGER                     :: dim1p(1), dim1(1)
    INTEGER                     :: dim3p(2), dim3(2)
    CHARACTER(LEN=max_dim_name) :: dim1n(1), dim3n(2)

    INTEGER :: i, read_status

    INCLUDE 'cbalance_ctl.inc'
    
    !! Check initialization

    if(module_initialized)  call finish("init_cbalance_bethy()","Attempt to initialize twice.")

    !! Printout parameters 

    call printCpoolParameters

    !! set parameters

    ntiles = cbalance%ntiles

    IF (ASSOCIATED(diag_stream)) THEN
       IO_diag => diag_stream
    ELSE
       CALL finish('init_cbalance_bethy', 'Diagnostic stream not present')
    END IF

    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'cbalance', filetype=fileformat, interval=io_time_event(1,'day','first',0))
          ! Set default stream options
          CALL default_stream_setting(stream, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
       IO_cbalance => stream
    ELSE
       ! Add new stream
       CALL new_stream(IO_cbalance, 'cbalance', filetype=fileformat, interval=io_time_event(1,'day','first',0))
       ! Set default stream options
       CALL default_stream_setting(IO_cbalance, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF

    !! Read namelist cbalance_ctl

    IF (p_parallel_io) THEN

       ! define default values
       read_cpools = .FALSE.
       cpool_file_name = 'cpools.nc'

       CALL position_nml ('CBALANCE_CTL', status=read_status)
       SELECT CASE (read_status)
       CASE (POSITIONED)
          READ (nml_unit, cbalance_ctl)
          CALL message('init_cbalance_bethy', 'Namelist CBALANCE_CTL: ')
          WRITE(nout, cbalance_ctl)
       END SELECT
    ENDIF

    IF (p_parallel_io) THEN
       IF (read_cpools) THEN
          CALL message('init_cbalance_bethy',' Reading initial values of the carbon pools from file '//TRIM(cpool_file_name))
       END IF
    END IF
    IF (p_parallel) THEN
       CALL p_bcast(read_cpools, p_io)
    END IF

    ! --- compute the number of time steps per day

    if( mod(86400,int(delta_time)) .ne. 0) then ! For computing the mean day temperature (see below) it is assumed that each day ..
                                                ! is computed with a fixed number of time steps. Therefore the program is  ..
                                                ! stopped wheen day length is not an integer multiple of the time step.
       call finish("init_cbalance_bethy","ERROR: Day length is not an integer multiple of the time step!")
    else 
       time_steps_per_day = 86400/int(delta_time)
    end if

    g_nland = grid%nland
    l_nland = domain%nland

    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    dim3p = (/ l_nland, ntiles /)
    dim3  = (/ g_nland, ntiles /)
    dim3n(1) = 'landpoint'
    dim3n(2) = 'tiles'

    ! --- Define state variables as stream elements

    CALL add(IO_cbalance,'Cpool_green'         ,cbalance%Cpool_green, lpost=.NOT.useDynveg .OR. fileformat==NETCDF, &
             longname='C-Pool for Green Parts of Vegetation',     units='mol(C) m-2(canopy)', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=150,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_woods'         ,cbalance%Cpool_woods, lpost=.NOT.useDynveg .OR. fileformat==NETCDF, &
             longname='C-Pool for Structural Material of Plants', units='mol(C) m-2(canopy)', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=151,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_reserve'     ,cbalance%Cpool_reserve, lpost=.NOT.useDynveg .OR. fileformat==NETCDF, &
             longname='C-Pool for reserve carbohydrates (starches, sugars)', units='mol(C) m-2(canopy)',&
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=152,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_litter_leaf' ,cbalance%Cpool_litter_leaf, lpost=.NOT.useDynveg .OR. fileformat==NETCDF, &
             longname='C-Pool for leaf litter',     units='mol(C) m-2(canopy)', contnorest=.TRUE.,      &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=177,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_litter_wood' ,cbalance%Cpool_litter_wood, lpost=.NOT.useDynveg .OR. fileformat==NETCDF, &
             longname='C-Pool for woody litter',    units='mol(C) m-2(canopy)', contnorest=.TRUE.,      &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=178,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_fast'          ,cbalance%Cpool_fast,  lpost=.NOT.useDynveg .OR. fileformat==NETCDF, &
             longname='C-Pool for fastly respirated soil organic material', units='mol(C) m-2(canopy)', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=153,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_slow'          ,cbalance%Cpool_slow, lpost=.NOT.useDynveg .OR. fileformat==NETCDF, &
             longname='C-Pool for slowly respirated soil organic material', units='mol(C) m-2(canopy)', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=154,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'NPP_Rate'            ,cbalance%NPP_Rate,     lpost=.FALSE.,                   &
             longname='Net Primary Production Rate', units='mol(CO2) m-2(canopy) s-1',                  &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=155,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'NPP_Rate_acc'        ,cbalance%NPP_Rate_acc, laccu=.TRUE., contnorest=.TRUE., &
             longname='Net Primary Production Rate (avg)', units='mol(CO2) m-2(canopy) s-1',            &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=158,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'LAI_sum'             ,cbalance%LAI_sum,             ldims=dim3p, gdims=dim3,  &
             dimnames=dim3n, code=156, lrerun=.TRUE., lpost=.FALSE.)
    CALL add(IO_cbalance,'NPP_sum'             ,cbalance%NPP_sum,             ldims=dim3p, gdims=dim3,  &
             dimnames=dim3n, code=157, lrerun=.TRUE., lpost=.FALSE.)
    CALL add(IO_cbalance,'GPP_sum'             ,cbalance%GPP_sum,             ldims=dim3p, gdims=dim3,  &
             dimnames=dim3n, code=158, lrerun=.TRUE., lpost=.FALSE.)
    CALL add(IO_cbalance,'topSoilTemp_sum'     ,cbalance%topSoilTemp_sum,     ldims=dim1p, gdims=dim1,  &
             dimnames=dim1n, code=159, lrerun=.TRUE., lpost=.FALSE.)
    CALL add(IO_cbalance,'alpha_sum'           ,cbalance%alpha_sum,           ldims=dim3p, gdims=dim3,  &
             dimnames=dim3n, code=160, lrerun=.TRUE., lpost=.FALSE.)
    CALL add(IO_cbalance,'soil_respiration'   ,cbalance%soil_respiration, lpost=.NOT.useDynveg .OR. fileformat==NETCDF, &
             longname='Soil respiration', units='mol(C) m-2(canopy)',                                           &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=156, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'NPP_flux_correction' ,cbalance%NPP_flux_correction,                                   &
             longname='Flux correction for NPP',   units='mol(CO2) m-2(canopy) s-1',                            &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=157, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)

    ! --- Define diagnostic variables as stream elements

    CALL add(IO_cbalance,'boxC_green',cbalance_diag%boxC_green,                                                 &
             longname='C-Pool for Green Parts of Vegetation',     units='mol(C) m-2(grid box)',                 &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=160, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'boxC_woods',cbalance_diag%boxC_woods,                                                 &
             longname='C-Pool for Structural Material of Plants', units='mol(C) m-2(grid box)',                 &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=161, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'boxC_reserve',cbalance_diag%boxC_reserve,                                             &
             longname='C-Pool for reserve carbohydrates (starches, sugars)', units='mol(C) m-2(grid box)',      &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=162, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'boxC_litter_leaf',cbalance_diag%boxC_litter_leaf,                                     &
             longname='C-Pool for leaf litter', units='mol(C) m-2(grid box)',                                   &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=179, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'boxC_litter_wood',cbalance_diag%boxC_litter_wood,                                     &
             longname='C-Pool for woody litter', units='mol(C) m-2(grid box)',                                  &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=159, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'boxC_fast',cbalance_diag%boxC_fast,                                                   &
             longname='C-Pool for fastly respirated soil organic material', units='mol(C) m-2(grid box)',       &
             ldims=dim3p,gdims= dim3, dimnames=dim3n, code=163, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value) 
    CALL add(IO_cbalance,'boxC_slow',cbalance_diag%boxC_slow,                                                   &
             longname='C-Pool for slowly respirated soil organic material',  units='mol(C) m-2(grid box)',      &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=164, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)

    CALL add(IO_cbalance,'LAI_previousDayMean',cbalance_diag%LAI_previousDayMean,              &
             longname='Mean Leaf Area Index of the Day Before Previous Day', units='', lpost=.FALSE., &
             ldims=dim3p,gdims= dim3, dimnames=dim3n, code=180, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'LAI_yDayMean',cbalance_diag%LAI_yDayMean,                            &
             longname='Mean Leaf Area Index of the Previous Day', units='', contnorest=.TRUE., &
             ldims=dim3p,gdims= dim3, dimnames=dim3n, code=165, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'NPP_yDayMean',cbalance_diag%NPP_yDayMean,                            &
             longname='Mean NPP Rate of the Previous Day',   units='mol(CO2) m-2(canopy) s-1', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=166, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'GPP_yDayMean',cbalance_diag%GPP_yDayMean,                            &
             longname='Mean GPP Rate of the Previous Day',   units='mol(CO2) m-2(canopy) s-1', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=167, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'topSoilTemp_yDayMean',cbalance_diag%topSoilTemp_yDayMean,            & 
             longname='Previous Day Mean Temperature of the Uppermost Soil Layer', units='K',  &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=168, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'alpha_yDayMean',cbalance_diag%alpha_yDayMean,                        &
             longname='Previous Day Mean Value of the Water Stress Coefficient', units='',     &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=169, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)


    CALL add(IO_cbalance,'box_soil_respiration',cbalance_diag%box_soil_respiration,                         &
             longname='Soil respiration', units='mol(C) m-2(grid box) s-1',                                     &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=170, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_NPP_yDayMean',cbalance_diag%box_NPP_yDayMean,                                 &
             longname='Mean NPP Rate of the Previous Day',   units='mol(CO2) m-2(grid box) s-1',            &
             ldims=dim3p,gdims= dim3, dimnames=dim3n, code=171, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_NPP_flux_correction',cbalance_diag%box_NPP_flux_correction,                   &
             longname='Flux correction for NPP',   units='mol(CO2) m-2(grid box) s-1',                      &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=172, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_GPP_yDayMean',cbalance_diag%box_GPP_yDayMean,                                 &
             longname='Mean GPP Rate of the Previous Day',   units='mol(CO2) m-2(grid box) s-1',            &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=173, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'litter_flux',cbalance_diag%litter_flux, lpost=.NOT.useDynveg .OR. fileformat==NETCDF, &
             longname='Total litter flux entering the soil pools', units='mol(CO2) m-2(canopy) s-1',   &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=174, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_litter_flux',cbalance_diag%box_litter_flux,                                     &
             longname='Total litter flux entering the soil pools', units='mol(CO2) m-2(grid box) s-1', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=175, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_Cpools_total',cbalance_diag%box_Cpools_total,                                 &
             longname='Sum of carbon from all carbon pools',   units='mol(CO2) m-2(grid box)',              &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=176, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)

    ! --- initializations -------------------------------

    IF (read_cpools) THEN

       CALL message('init_cbalance_bethy','Reading initial carbon pools from file')

       ALLOCATE(init_Cpool_green(l_nland,ntiles))
       ALLOCATE(init_Cpool_woods(l_nland,ntiles))
       ALLOCATE(init_Cpool_reserve(l_nland,ntiles))
       ALLOCATE(init_Cpool_litter_leaf(l_nland,ntiles))
       ALLOCATE(init_Cpool_litter_wood(l_nland,ntiles))
       ALLOCATE(init_Cpool_fast(l_nland,ntiles))
       ALLOCATE(init_Cpool_slow(l_nland,ntiles))
       
       ALLOCATE(zreal2d(domain%ndim,domain%nblocks))

       !! --- Get green pool (and allocate further memory)

       IF (p_parallel_io) THEN
          cpool_file%opened = .FALSE.
          CALL IO_open(TRIM(cpool_file_name), cpool_file, IO_READ)
          IO_file_id = cpool_file%file_id

          ALLOCATE(zreal3d(grid%nlon,grid%nlat,ntiles))

          CALL IO_inq_varid(IO_file_id, 'Cpool_green', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       END IF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
          init_Cpool_green(:,i) = PACK(zreal2d, MASK=domain%mask)
       END DO

       !! --- Get wood pool

       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'Cpool_woods', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       END IF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
          init_Cpool_woods(:,i) = PACK(zreal2d, MASK=domain%mask)
       END DO

       !! --- Get reserve pool

       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'Cpool_reserve', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       END IF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
          init_Cpool_reserve(:,i) = PACK(zreal2d, MASK=domain%mask)
       END DO

       !! --- Get leaf litter pool

       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'Cpool_litter_leaf', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       END IF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
          init_Cpool_litter_leaf(:,i) = PACK(zreal2d, MASK=domain%mask)
       END DO

       !! --- Get wood litter pool

       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'Cpool_litter_wood', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       END IF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
          init_Cpool_litter_wood(:,i) = PACK(zreal2d, MASK=domain%mask)
       END DO

       !! --- Get fast pool

       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'Cpool_fast', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
       END IF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
          init_Cpool_fast(:,i) = PACK(zreal2d, MASK=domain%mask)
       END DO

       !! --- Get slow pool

       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'Cpool_slow', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
       END IF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
          init_Cpool_slow(:,i) = PACK(zreal2d, MASK=domain%mask)
       END DO

       !! --- Finish

       IF (p_parallel_io) THEN
          CALL IO_close(cpool_file)
          DEALLOCATE(zreal3d)
       END IF
       DEALLOCATE(zreal2d)

       IF (ANY(ABS(init_Cpool_green) > 10000._dp)) THEN
          CALL message('init_cbalance_bethy','Found cpool values in ini file with absolute value > 10000')
          CALL message('                   ','  This is likely due to using a cpool file that was')
          CALL message('                   ','  generated with a different land-sea mask.')
          CALL message('                   ','  Check your setup!')
          CALL finish('init_cbalance_bethy','')
       END IF

    ELSE
       cbalance%Cpool_green(:,:)   = 0.0_dp ![mol[C]/m2]
       cbalance%Cpool_woods(:,:)   = 0.0_dp
       cbalance%Cpool_reserve(:,:) = 0.0_dp
       cbalance%Cpool_litter_leaf(:,:) = 60.0_dp
       cbalance%Cpool_litter_wood(:,:) = 180.0_dp
       cbalance%Cpool_fast(:,:)        = 60.0_dp
       cbalance%Cpool_slow(:,:)        = 2400.0_dp
       CALL message('init_cbalance_bethy','Carbon pools initialized')
    END IF

    cbalance%NPP_Rate = 0.0_dp
    cbalance%NPP_Rate_acc = 0.0_dp
    cbalance%soil_respiration = 0.0_dp
    cbalance%NPP_flux_correction = 0.0_dp

    cbalance_diag%boxC_green = 0.0_dp
    cbalance_diag%boxC_woods = 0.0_dp
    cbalance_diag%boxC_reserve = 0.0_dp
    cbalance_diag%boxC_litter_leaf = 0.0_dp
    cbalance_diag%boxC_litter_wood = 0.0_dp
    cbalance_diag%boxC_fast = 0.0_dp
    cbalance_diag%boxC_slow = 0.0_dp
    cbalance_diag%LAI_yDayMean = 0.0_dp
    cbalance_diag%LAI_previousDayMean = 0.0_dp
    cbalance_diag%NPP_yDayMean = 0.0_dp
    cbalance_diag%GPP_yDayMean = 0.0_dp
    cbalance_diag%topSoilTemp_yDayMean = 0.0_dp
    cbalance_diag%alpha_yDayMean = 0.0_dp
    cbalance_diag%box_soil_respiration = 0.0_dp
    cbalance_diag%box_NPP_yDayMean = 0.0_dp
    cbalance_diag%box_NPP_flux_correction = 0.0_dp
    cbalance_diag%box_GPP_yDayMean = 0.0_dp
    cbalance_diag%litter_flux = 0.0_dp
    cbalance_diag%box_litter_flux = 0.0_dp
    cbalance_diag%box_Cpools_total = 0.0_dp

  END SUBROUTINE init_cbalance_bethy

  ! --- update_cbalance_bethy()  --------------------------------------------------------------------------------------------------

  SUBROUTINE update_cbalance_bethy(domain, &
                                   cbalance, &
                                   lctlib, &
                                   surface, &
                                   grossAssimilation, &
                                   darkRespiration, &
                                   topSoilTemp, &
                                   alpha, &
                                   LAI, &
                                   veg_fract_correction, &
                                   net_CO2_flux) 
    
    USE mo_jsbach,                 ONLY: debug
    USE mo_time_control,           ONLY: lstart,lresume
    USE mo_land_surface,           ONLY: land_surface_type
    USE mo_jsbach_lctlib,          ONLY: lctlib_type
    USE mo_exception,              ONLY: message
    USE mo_time_control,           ONLY: get_date_components, current_date, previous_date
    USE mo_cbal_cpools,            ONLY: update_Cpools
#if defined (__SX__) && defined (_OPENMP)
    USE omp_lib,                   ONLY: omp_get_thread_num
#endif

    TYPE(domain_type),      INTENT(in)    :: domain  ! Information on the relation between global and local grid point indexing
    TYPE(cbalance_type),    INTENT(inout) :: cbalance
    TYPE(lctlib_type),      INTENT(in)    :: lctlib
    TYPE(land_surface_type),INTENT(in)    :: surface

    REAL(dp),INTENT(in)               :: grossAssimilation(:,:)    !! Gross primary product (where photorespiration has already  been ..
    !                                                              !! .. accounted for) [mol(CO2)/(m^2 s)]
    REAL(dp),INTENT(in)               :: darkRespiration(:,:)      !! dark respiration of leaves [mol(CO2)/(m^2 s)]
    REAL(dp),INTENT(in)               :: topSoilTemp(:)            !! As temperature of the fast pool the temperature in the top soil layer (on all tiles identical)
                                                                   !! .. is taken [Kelvin]
    REAL(dp),INTENT(in)               :: alpha(:,:)                !! relative soil humidity factor for heterotrophic respiration
    REAL(dp),INTENT(in)               :: LAI(:,:)                  !! leaf area index of (current time step)
    REAL(dp),INTENT(in)               :: veg_fract_correction(:,:) !! Correction factor for cover fractions 1-exp(-LAI_max/2)
    REAL(dp),INTENT(out)              :: net_CO2_flux(:)           !! grid cell average of net CO2-flux between biosphere and atmosphere [kg(CO2)/(m^2(ground) s)]
                                                                   !! .. (positive for emission to atmosphere, negative for absorption by biosphere)
    ! local variables

    integer :: itile
    integer :: kidx0,kidx1
    integer :: day_of_month,day_of_month_at_prev_ts

    real(dp) :: frac_npp_2_woodPool(1:nidx,1:ntiles)
    real(dp) :: frac_npp_2_reservePool(1:nidx,1:ntiles)
    real(dp) :: tau_Cpool_litter_leaf(1:nidx,1:ntiles)
    real(dp) :: tau_Cpool_litter_wood(1:nidx,1:ntiles)
    real(dp) :: LAI_shed_constant(1:nidx,1:ntiles)
    real(dp) :: frac_C_fast2atmos(1:nidx,1:ntiles)
    real(dp) :: Max_C_content_woods(1:nidx,1:ntiles)
    real(dp) :: specificLeafArea_C(nidx,1:ntiles)
    real(dp) :: areaWeightingFactor(1:nidx,1:ntiles)

    logical  :: first_ts_of_day  !! .true.: current time step is first time step of day; .false. : not first time step of day

#if defined (__SX__) && defined (_OPENMP)
    INTEGER :: tid
#endif

#if defined (__SX__) && defined (_OPENMP)
    tid = omp_get_thread_num()
#endif

    ! --- set block range indices

    kidx0 = kstart    ! first value of kindex() (the index of the first element of a block in the processor domain)
    kidx1 = kend      ! last value of kindex()  (the index of the last element of a block in the processor domain)

    ! initializations

    net_CO2_flux(1:nidx) = 0._dp

    IF (read_cpools .AND. (lstart .OR. lresume)) THEN
       cbalance%Cpool_green(kidx0:kidx1,:)   = init_Cpool_green(kidx0:kidx1,:)
       cbalance%Cpool_woods(kidx0:kidx1,:)   = init_Cpool_woods(kidx0:kidx1,:)
       cbalance%Cpool_reserve(kidx0:kidx1,:) = init_Cpool_reserve(kidx0:kidx1,:)
       cbalance%Cpool_litter_leaf(kidx0:kidx1,:) = init_Cpool_litter_leaf(kidx0:kidx1,:)
       cbalance%Cpool_litter_wood(kidx0:kidx1,:) = init_Cpool_litter_wood(kidx0:kidx1,:)
       cbalance%Cpool_fast(kidx0:kidx1,:)    = init_Cpool_fast(kidx0:kidx1,:)
       cbalance%Cpool_slow(kidx0:kidx1,:)    = init_Cpool_slow(kidx0:kidx1,:)
#if defined (__SX__) && defined (_OPENMP)
       WRITE(0,*) 'Cpools overwritten with initial values for kidx0,kidx1,nland on PE/thread: ',kidx0,kidx1,domain%nland,p_pe,tid
#else
       WRITE(0,*) 'Cpools overwritten with initial values for kidx0,kidx1,nland on PE: ',kidx0,kidx1,domain%nland,p_pe
#endif
    END IF

    ! --- update first_ts_of_day

    CALL get_date_components(current_date,DAY=day_of_month)
    CALL get_date_components(previous_date, DAY=day_of_month_at_prev_ts)
    first_ts_of_day = ( day_of_month /= day_of_month_at_prev_ts )
    IF (debug .AND. first_ts_of_day) CALL message('update_cbalance_bethy','First time step of new day')

    ! --- go ...

    ! compute net primary production rate

    cbalance%NPP_rate(kidx0:kidx1,:) = NPP_rate_bethy(grossAssimilation(1:nidx,:),darkRespiration(1:nidx,:))

    cbalance%NPP_rate_acc(kidx0:kidx1,:) = cbalance%NPP_rate_acc(kidx0:kidx1,:) + cbalance%NPP_rate(kidx0:kidx1,:) * delta_time

    ! Prepare area weighting factor to rescale from 1/[m^2(canopy)] to 1/[m^2(grid box)]
    areaWeightingFactor(:,:) = veg_fract_correction(:,:) * surface%cover_fract(kidx0:kidx1,:) &
                               * SPREAD(surface%veg_ratio_max(kidx0:kidx1),DIM=2,NCOPIES=ntiles)


    ! +++++

    IF( .NOT. first_ts_of_day .OR. lstart) THEN ! perform daily sums ---------------------------------------------------

       cbalance%LAI_sum(kidx0:kidx1,:)          = cbalance%LAI_sum(kidx0:kidx1,:)         + lai(1:nidx,:)
       cbalance%NPP_sum(kidx0:kidx1,:)          = cbalance%NPP_sum(kidx0:kidx1,:)         + cbalance%NPP_rate(kidx0:kidx1,:)
       cbalance%GPP_sum(kidx0:kidx1,:)          = cbalance%GPP_sum(kidx0:kidx1,:)         + grossAssimilation(1:nidx,:)
       cbalance%topSoilTemp_sum(kidx0:kidx1)    = cbalance%topSoilTemp_sum(kidx0:kidx1)   + topSoilTemp(1:nidx)
       cbalance%alpha_sum(kidx0:kidx1,:)        = cbalance%alpha_sum(kidx0:kidx1,:)       + alpha(1:nidx,:)

    ELSE ! A new day begins ==> perform carbon balance -------------------------------------------------------------------
       
      ! Compute previous days means 

      cbalance_diag%LAI_previousDayMean(kidx0:kidx1,:)  = cbalance_diag%LAI_yDayMean(kidx0:kidx1,:)
      cbalance_diag%LAI_yDayMean(kidx0:kidx1,:)         = cbalance%LAI_sum(kidx0:kidx1,:)/time_steps_per_day
      cbalance_diag%NPP_yDayMean(kidx0:kidx1,:)         = cbalance%NPP_sum(kidx0:kidx1,:)/time_steps_per_day
      cbalance_diag%GPP_yDayMean(kidx0:kidx1,:)         = cbalance%GPP_sum(kidx0:kidx1,:)/time_steps_per_day
      cbalance_diag%topSoilTemp_yDayMean(kidx0:kidx1)   = cbalance%topSoilTemp_sum(kidx0:kidx1)/time_steps_per_day
      cbalance_diag%alpha_yDayMean(kidx0:kidx1,:)       = cbalance%alpha_sum(kidx0:kidx1,:)/time_steps_per_day

      ! Restart summing of previous days values

      cbalance%LAI_sum(kidx0:kidx1,:)             = lai(1:nidx,:)
      cbalance%NPP_sum(kidx0:kidx1,:)             = cbalance%NPP_Rate(kidx0:kidx1,:)
      cbalance%GPP_sum(kidx0:kidx1,:)             = grossAssimilation(1:nidx,:)
      cbalance%topSoilTemp_sum(kidx0:kidx1)       = topSoilTemp(1:nidx)
      cbalance%alpha_sum(kidx0:kidx1,:)           = alpha(1:nidx,:)

      ! Construct arrays with land cover type (PFT) properties (SLA, max. wood mass, etc.)

      DO itile=1,ntiles
         frac_npp_2_woodPool(1:nidx,itile)    = lctlib%frac_npp_2_woodPool(surface%cover_type(kidx0:kidx1,itile))
         frac_npp_2_reservePool(1:nidx,itile) = lctlib%frac_npp_2_reservePool(surface%cover_type(kidx0:kidx1,itile))
         tau_Cpool_litter_leaf(1:nidx,itile)  = lctlib%tau_Cpool_litter_leaf(surface%cover_type(kidx0:kidx1,itile))
         tau_Cpool_litter_wood(1:nidx,itile)  = lctlib%tau_Cpool_litter_wood(surface%cover_type(kidx0:kidx1,itile))
         LAI_shed_constant(1:nidx,itile)      = lctlib%LAI_shed_constant(surface%cover_type(kidx0:kidx1,itile))
         frac_C_fast2atmos(1:nidx,itile)      = lctlib%frac_C_fast2atmos(surface%cover_type(kidx0:kidx1,itile))
         Max_C_content_woods(1:nidx,itile)    = lctlib%Max_C_content_woods(surface%cover_type(kidx0:kidx1,itile))
         specificLeafArea_C(1:nidx,itile)     = lctlib%specificLeafArea_C(surface%cover_type(kidx0:kidx1,itile))
      END DO

      ! update of carbon pools and compute soil respiration rate

      CALL update_Cpools(cbalance_diag%LAI_yDayMean(kidx0:kidx1,:),          &
                         cbalance_diag%LAI_previousDayMean(kidx0:kidx1,:),   &
                         cbalance_diag%NPP_yDayMean(kidx0:kidx1,:),          &
                         SPREAD(cbalance_diag%topSoilTemp_yDayMean(kidx0:kidx1),DIM=2,NCOPIES=ntiles),  &
                         SPREAD(cbalance_diag%alpha_yDayMean(kidx0:kidx1,1),DIM=2,NCOPIES=ntiles),      &
                         frac_npp_2_woodPool(1:nidx,:),                      &
                         frac_npp_2_reservePool(1:nidx,:),                   &
                         tau_Cpool_litter_leaf(1:nidx,:),                    &
                         tau_Cpool_litter_wood(1:nidx,:),                    &
                         LAI_shed_constant(1:nidx,:),                        &
                         frac_C_fast2atmos(1:nidx,:),                        &
                         Max_C_content_woods(1:nidx,:),                      &
                         specificLeafArea_C(1:nidx,:),                       &
                         cbalance%Cpool_green(kidx0:kidx1,:),                &
                         cbalance%Cpool_woods(kidx0:kidx1,:),                &
                         cbalance%Cpool_reserve(kidx0:kidx1,:),              &
                         cbalance%Cpool_litter_leaf(kidx0:kidx1,:),          &
                         cbalance%Cpool_litter_wood(kidx0:kidx1,:),          &
                         cbalance%Cpool_fast(kidx0:kidx1,:),                 &
                         cbalance%Cpool_slow(kidx0:kidx1,:),                 &
                         cbalance%soil_respiration(kidx0:kidx1,:),           &
                         cbalance%NPP_flux_correction(kidx0:kidx1,:),        &
                         cbalance_diag%litter_flux(kidx0:kidx1,:)            &
                        )

      ! Compute carbon contents per grid box area by weighting pools with fractions of grid box covered by vegetation
      ! -------------------------------------------------------------------------------------------------------------
      
      cbalance_diag%boxC_green(kidx0:kidx1,:)   = cbalance%Cpool_green(kidx0:kidx1,:)   * areaWeightingFactor(:,:)
      cbalance_diag%boxC_reserve(kidx0:kidx1,:) = cbalance%Cpool_reserve(kidx0:kidx1,:) * areaWeightingFactor(:,:)
      cbalance_diag%boxC_woods(kidx0:kidx1,:)   = cbalance%Cpool_woods(kidx0:kidx1,:)   * areaWeightingFactor(:,:)
      cbalance_diag%boxC_litter_leaf(kidx0:kidx1,:) = cbalance%Cpool_litter_leaf(kidx0:kidx1,:) * areaWeightingFactor(:,:)
      cbalance_diag%boxC_litter_wood(kidx0:kidx1,:) = cbalance%Cpool_litter_wood(kidx0:kidx1,:) * areaWeightingFactor(:,:)
      cbalance_diag%boxC_fast(kidx0:kidx1,:)    = cbalance%Cpool_fast(kidx0:kidx1,:)    * areaWeightingFactor(:,:)
      cbalance_diag%boxC_slow(kidx0:kidx1,:)    = cbalance%Cpool_slow(kidx0:kidx1,:)    * areaWeightingFactor(:,:)

#ifndef __PGI
      cbalance_diag%box_Cpools_total(kidx0:kidx1) = SUM(    cbalance_diag%boxC_green(kidx0:kidx1,:)   &
                                                          + cbalance_diag%boxC_woods(kidx0:kidx1,:)   &
                                                          + cbalance_diag%boxC_reserve(kidx0:kidx1,:) &
                                                          + cbalance_diag%boxC_litter_leaf(kidx0:kidx1,:) &
                                                          + cbalance_diag%boxC_litter_wood(kidx0:kidx1,:) &
                                                          + cbalance_diag%boxC_fast(kidx0:kidx1,:)    &
                                                          + cbalance_diag%boxC_slow(kidx0:kidx1,:)    &
                                                        ,DIM=2)
#else
      !
      ! Differnt notation to allow compilation with the PGI compiler (pgf95 6.1-1) on Linux
      ! 
      cbalance_diag%box_Cpools_total(kidx0:kidx1) = SUM(cbalance_diag%boxC_green(kidx0:kidx1,:),DIM=2)   &
                                                  + SUM(cbalance_diag%boxC_woods(kidx0:kidx1,:),DIM=2)   &
                                                  + SUM(cbalance_diag%boxC_reserve(kidx0:kidx1,:),DIM=2) &
                                                  + SUM(cbalance_diag%boxC_litter_leaf(kidx0:kidx1,:),DIM=2) &
                                                  + SUM(cbalance_diag%boxC_litter_wood(kidx0:kidx1,:),DIM=2) &
                                                  + SUM(cbalance_diag%boxC_fast(kidx0:kidx1,:),DIM=2)    &
                                                  + SUM(cbalance_diag%boxC_slow(kidx0:kidx1,:),DIM=2)
#endif

      ! Other diagnostics
      ! -----------------

      cbalance_diag%box_soil_respiration(kidx0:kidx1,:)   = cbalance%soil_respiration(kidx0:kidx1,:)    * areaWeightingFactor(:,:)
      cbalance_diag%box_NPP_yDayMean(kidx0:kidx1,:)       = cbalance_diag%NPP_yDayMean(kidx0:kidx1,:)   * areaWeightingFactor(:,:)
      cbalance_diag%box_NPP_flux_correction(kidx0:kidx1,:)= cbalance%NPP_flux_correction(kidx0:kidx1,:) * areaWeightingFactor(:,:)
      cbalance_diag%box_GPP_yDayMean(kidx0:kidx1,:)       = cbalance_diag%GPP_yDayMean(kidx0:kidx1,:)   * areaWeightingFactor(:,:)
      cbalance_diag%box_litter_flux(kidx0:kidx1,:)        = cbalance_diag%litter_flux(kidx0:kidx1,:)    * areaWeightingFactor(:,:)

   end if
   
   ! Compute net CO2-flux exchanged with atmosphere at each time step
   !-----------------------------------------------------------------
   ! Note: carbon loss of biosphere means a positive C-flux to atmosphere (i.e. NEP and net CO2-flux have opposite signs)

   net_CO2_flux(1:nidx) = - sum(areaWeightingFactor(1:nidx,:)                     & !! Minus because of convention (atmosphere gain is positive)
                                 * (    cbalance%NPP_Rate(kidx0:kidx1,:)          & !! .. actual NPP rate in this time step
                                      + cbalance%NPP_flux_correction(kidx0:kidx1,:)    & !! .. flux corrected NPP rate 
                                      + cbalance%soil_respiration(kidx0:kidx1,:)  & !! ..
                                   ),DIM=2)                                       & !! ..
                            * molarMassCO2_kg                                       !! .. times conversion factor from mol to kg CO2

 END SUBROUTINE update_cbalance_bethy


  ! --- NPP_rate_bethy() ----------------------------------------------------------------------------------------------------------

  ! Following BETHY the net primary productivity is estimated from dark respiration (R_d) and  gross primary production (GPP). 
  ! More precisely:  
  !
  !    (1) NPP = GPP - R_m - R_g,
  !
  ! where "R_g" is the growth respiration and "R_m" the maintenance respiration. Maintenance respiration can be estimated from
  ! dark respiration by
  !
  !                 R_d
  !    (2) R_m = ----------             (Eq. (128) in Knorr (in molar units))
  !              f_aut_leaf
  !
  ! where "f_aut_leaf" is the leaf fraction of plant-total (autotrophic) respiration. Note that "GPP" is here the Farquar
  ! productivity, where it is already accounted for photorespiration.  To estimate NPP from (1) it therefore remains to 
  ! determine the growth respiration. There is no growth respiration for NPP<0 and is otherwise a certain fraction of NPP:
  !
  !               / (cCost -1) NPP  for NPP>0
  !    (3) R_g = <
  !               \ 0               otherwise,
  !
  ! where "cCost" are the relative costs (measured in carbon) to produce 1 unit of carbon. Entering this into (1) one thus finds
  !
  !
  !                      cCost -1
  !    (4) R_g = MAX(0,  --------- (GPP -R_m)) 
  !                        cCost
  !
  ! See: W. Knorr, "Satellite Remote Sensing and Modelling of the Global CO2 Exchange of Land Vegetation: A Synthesis Study",
  !      Examensarbeit Nr. 49, (MPI for Meterology, Hamburg, 1998).
  !

  elemental pure function NPP_rate_bethy(grossAssimilation,darkRespiration)
    real(dp),intent(in) :: grossAssimilation ! Gross primary product (where photorespiration has already  ..
                                                                 ! .. been accounted for) [mol(CO2)/(m^2 s)]
    real(dp),intent(in) :: darkRespiration   ! dark respiration of leaves [mol(CO2)/(m^2 s)]
    real(dp)            :: NPP_rate_bethy    ! net primary production rate [mol(CO2)/(m^2 s)]

    ! locals

    REAL(dp), parameter :: f_aut_leaf             = 0.40_dp      ! leaf fraction of plant-total (autotrophic) respiration
    REAL(dp), parameter :: cCost                  = 1.25_dp      ! relative costs (measured in carbon) to produce 1 unit of carbon

    real(dp),parameter   :: hlp =  (cCost -1._dp)/cCost
    real(dp)             :: maintenanceRespiration ! R_m
    real(dp)             :: growthRespiration      ! R_g

    ! Go ...
    
    maintenanceRespiration = darkRespiration/f_aut_leaf
    growthRespiration = max(0._dp,hlp*(grossAssimilation - maintenanceRespiration))
    NPP_rate_bethy = grossAssimilation - maintenanceRespiration - growthRespiration

  end function NPP_rate_bethy

end module mo_cbal_bethy

