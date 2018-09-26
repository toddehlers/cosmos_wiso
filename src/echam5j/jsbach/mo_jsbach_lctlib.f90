MODULE mo_jsbach_lctlib
  !
  ! This module provides the data from the landcover library file, which contains additional information on the 
  ! landcover data used in a particular run of JSBACH.
  ! The name of the landcover library file is specified in the JSBACH configuration file (run.def) under the 
  ! keyword "LCT_FILE".
  !
  ! Authors: Reiner Schnur, 7/2003
  !          Christian H. Reick, 8/2003
  !
  USE mo_mpi,       ONLY: p_parallel, p_parallel_io, p_bcast, p_io
  USE mo_kind,      ONLY: dp
  USE mo_exception, ONLY: message, finish
  USE mo_netcdf,    ONLY: nf_max_name
  USE mo_jsbach_constants, ONLY : C3,C4

  IMPLICIT NONE

  ! === BEGIN OF PUBLIC PART =======================================================================================================

  PUBLIC :: init_lctlib ! This subroutine initializes lctlib, i.e. it allocates memory and reads in the data


  ! --- lctlib ---------------------------------------------------------------------------------------------------------------------

  TYPE lctlib_type
     ! --- general parameters ------------------------------------------------------------------------------------------------------
     INTEGER           :: nlct                     !! Number of landcover types
     INTEGER           :: npft                     !! Number of vegetation types, including crop
     INTEGER, POINTER  :: LctClass(:)              !! Landcover class identifier
     CHARACTER(LEN=16), POINTER  :: LctName(:)     !! Names of landcover types
     LOGICAL, POINTER  :: NaturalVegFlag(:)        !! "true" for natural vegetation landcover types
     LOGICAL, POINTER  :: ForestFlag(:)            !! "ture" for forest landcover types
     LOGICAL, POINTER  :: CropFlag(:)              !! "true" for agricultural landcover types
     LOGICAL, POINTER  :: LakeFlag(:)              !! "true" for a "Vegetation" class that signifies lake or ocean
     LOGICAL, POINTER  :: GlacierFlag(:)           !! "true" for a land surface class that signifies glaciers
     LOGICAL, POINTER  :: BareSoilFlag(:)          !! "true" if not vegetation, lake or glacier
     ! --- Albedo -------------------------------------------------------------------------------------------------------------------
     REAL(dp), POINTER :: AlbedoSnowVisMin(:)      !! Minimum snow albedo in the visible range
     REAL(dp), POINTER :: AlbedoSnowVisMax(:)      !! Maximum snow albedo in the visible range
     REAL(dp), POINTER :: AlbedoSnowNirMin(:)      !! Minimum snow albedo in the NIR range
     REAL(dp), POINTER :: AlbedoSnowNirMax(:)      !! Maximum snow albedo in the NIR range
     REAL(dp), POINTER :: AlbedoSnowMin(:)         !! Minimum snow albedo
     REAL(dp), POINTER :: AlbedoSnowMax(:)         !! Maximum snow albedo
     REAL(dp), POINTER :: AlbedoCanopyVIS(:)       !! Albedo of the canopy (vegetation) in the visible range
     REAL(dp), POINTER :: AlbedoCanopyNIR(:)       !! Albedo of the canopy (vegetation) in the NIR range
     ! --- parameters used by bethy ------------------------------------------------------------------------------------------------
     LOGICAL, POINTER  :: NitrogenScalingFlag(:)   !! Indicates, whether nitrogen scaling shall be applied to that vegetation type
     LOGICAL, POINTER  :: C4flag(:)                !! Photosynthetic pathway: C4=.true or C3=.false.
     REAL(dp), POINTER :: CarboxRate(:)            !! Maximum carboxilation rate at 25 degrees Celsius [1.E-6 * Mol(CO2)/m^2/s] ...
                                                   !! ... (Table 2.6 in Knorr)
     REAL(dp), POINTER :: ETransport(:)            !! Maximum electron transport rate at 25 degrees Celsius [1.E-6 * Mol/m^2/s] ...
                                                   !! ... (Table 2.6 in Knorr)
     REAL(dp), POINTER :: VegHeight(:)             !! Typical height of the vegetation classes [m]
     REAL(dp), POINTER :: VegRoughness(:)          !! Typical roughness length of the vegetation classes [m]
     ! --- Parameters used in Phenology scheme -------------------------------------------------------------------------------------
     INTEGER, POINTER  :: PhenologyType(:)         !! Phenology type (only for natural vegetation) ...
                                                   !! ... (none: 0; evergreen: 1; summergreen: 2; raingreen: 3; grasses: 4)
     REAL(dp), POINTER :: MaxLAI(:)                !! Upper LAI boundary for LoGoP-Scheme (phenology) in [m2/m2]
     REAL(dp), POINTER :: specificLeafArea_C(:)    !! Carbon content per leaf area in [mol(Carbon)/m^2(leaf)]
     ! --- Parameters used in the carbon balance model -----------------------------------------------------------------------------
     REAL(dp), POINTER :: frac_npp_2_woodPool(:)   !! Fraction of NPP that is maximally put into the wood pool of the carbon balance model
     REAL(dp), POINTER :: frac_npp_2_reservePool(:)!! Fraction of NPP that is optimally put into the reserve pool of the carbon balance model
     REAL(dp), POINTER :: tau_Cpool_litter_leaf(:) !! Time constant by which leaf litter pool is depreciated  [days]
     REAL(dp), POINTER :: tau_Cpool_litter_wood(:) !! Time constant by which woody litter pool is depreciated  [days]
     REAL(dp), POINTER :: LAI_shed_constant(:)     !! 
     REAL(dp), POINTER :: frac_C_fast2atmos(:)     !! Fraction of heterotrophic loss of fast pool emitted to the atmosphere (rest enters slow pool)
     REAL(dp), POINTER :: Max_C_content_woods(:)   !! Maximum carbon content in the woody parts of plants [mol(C)/m^2(canopy)] (for forests
                                                   !! .. this is closely related to the yield, usually measured in m^3(wood)/hectar)
     ! --- Other parameters --------------------------------------------------------------------------------------------------------
     REAL(dp), POINTER :: CanopyResistanceMin(:)  !! Minimum canopy resistance (optional, only used in VIC scheme
                                                  !! and if BETHY is not used)
     REAL(dp), POINTER :: StemArea(:)             !! Area of stems and branches of woody plants
  END TYPE lctlib_type
  PUBLIC :: lctlib_type                            

  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ===============================================================================

  PRIVATE 


  ! Save nlct and npft for convenience
  INTEGER, SAVE             :: nlct ! The number of land cover types
  INTEGER, SAVE             :: npft ! Number of plant functional types, including crop

  ! === Private declarations =======================================================================================================

CONTAINS

  SUBROUTINE init_lctlib(lctlib_file_name, lctlib)
    !
    ! Get landcover library, i.e. lookup table for landcover types
    !
    ! The structure of the lookup table in the input file is as follows. It can contain only three types of lines in arbitrary
    ! order:
    !
    !     Comment lines: These contain before a '#' only blanks
    !       Blank lines: These contain only blanks
    !        Data lines: The first nonblank characters form a keyword which is follwed by the data (separated by blanks)
    !
    ! Keywords are the same as the names of the lctlib components, with one exception:
    !
    !        NLCT : This is the keyword after which the number of landcover types, i.e. the number of data columns
    !               in the file is expected. This number has to be identical with the number of landcover types
    !               from the landcover data file (keyword: LCTLIB_FILE in the JSBACH configuration file).
    !               
    ! After all other keywords NLCT columns of data are expected.
    !
    ! Example:
    !            # ---- LANDCOVER LIBRARY FILE -------
    !            NLCT 3
    !            LctClass               3 5 9
    !            # Phenology types: 0= none 1=summergreen, 2=evergreen, 3=raingreen, 4=grasses
    !            PhenologyType          2 2 4
    !            # C4flag: 0=C3, 1=C4
    !            C4flag                 0 1 1
    !
    ! The first string on each line (except first line and comment lines) must correspond to the name of the lctlib component
    ! (case doesn't matter)
    ! For lctlib components of type LOGICAL use 0/1 in the file to indicate .FALSE./.TRUE.
    ! For lctlib components of type INTEGER the numbers on the line must be integer values
    ! For lctlib components of type REAL the numbers on the line can be either integer or floating point
    ! Comments can appear on any line, everything to the right of and including "#" is disregarded
    !
    ! The file can contain more keywords than needed --- therefore the same file can be used by several model components.
    !
    USE mo_util_string, ONLY: tolower

    CHARACTER(nf_max_name), INTENT(in) :: lctlib_file_name
    TYPE(lctlib_type),    INTENT(out) :: lctlib ! LCT library to initialize

    INTEGER, PARAMETER            :: lctlib_file_unit = 66

    CHARACTER(len=30)  :: key
    CHARACTER(len=256) :: line
    INTEGER            :: pos_comment, read_status
    integer            :: pos,length
    character(len=2)   :: blank_set = " "//achar(9) !! the blank characters: BLANK and TAB
    integer,allocatable:: itmp(:)                   !! temporary array used for input of logicals

    LOGICAL            :: exist_LctClass               = .FALSE.
    LOGICAL            :: exist_LctName                = .FALSE.
    LOGICAL            :: exist_NaturalVegFlag         = .FALSE.
    LOGICAL            :: exist_ForestFlag             = .FALSE.
    LOGICAL            :: exist_CropFlag               = .FALSE.
    LOGICAL            :: exist_LakeFlag               = .FALSE.
    LOGICAL            :: exist_GlacierFlag            = .FALSE.
    LOGICAL            :: exist_BareSoilFlag           = .FALSE.
    LOGICAL            :: exist_AlbedoSnowVisMin       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowVisMax       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowNirMin       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowNirMax       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowMin          = .FALSE.
    LOGICAL            :: exist_AlbedoSnowMax          = .FALSE.
    LOGICAL            :: exist_AlbedoCanopyVIS        = .FALSE.
    LOGICAL            :: exist_AlbedoCanopyNIR        = .FALSE.
    LOGICAL            :: exist_nitrogenScalingFlag    = .FALSE.
    LOGICAL            :: exist_C4flag                 = .FALSE.
    LOGICAL            :: exist_CarboxRate             = .FALSE.
    LOGICAL            :: exist_ETransport             = .FALSE.
    LOGICAL            :: exist_VegHeight              = .FALSE.
    LOGICAL            :: exist_VegRoughness           = .FALSE.
    LOGICAL            :: exist_PhenologyType          = .FALSE.
    LOGICAL            :: exist_CanopyResistanceMin    = .FALSE.
    LOGICAL            :: exist_MaxLAI                 = .FALSE.
    LOGICAL            :: exist_StemArea               = .FALSE.
    LOGICAL            :: exist_specificLeafArea_C     = .FALSE.
    LOGICAL            :: exist_frac_npp_2_woodPool    = .FALSE.
    LOGICAL            :: exist_frac_npp_2_reservePool = .FALSE.
    LOGICAL            :: exist_tau_Cpool_litter_leaf  = .FALSE.
    LOGICAL            :: exist_tau_Cpool_litter_wood  = .FALSE.
    LOGICAL            :: exist_LAI_shed_constant      = .FALSE.
    LOGICAL            :: exist_frac_C_fast2atmos      = .FALSE.
    LOGICAL            :: exist_Max_C_content_woods    = .FALSE.
    integer            :: i

    !nlct = lctlib%nlct

    ! Determine name of landcover library file

    ! Read lctlib data and eventually allocate also memory for lctlib data on io-processor

    IF (p_parallel_io) THEN

       OPEN(unit=lctlib_file_unit, file=lctlib_file_name, form='FORMATTED', status='OLD', iostat=read_status)
       IF (read_status /= 0) CALL finish('init_lctlib','Error opening landcover library file')

       ! --- find keyword NLCT in landcover library file

       DO
          READ(lctlib_file_unit,'(A256)', IOSTAT=read_status) line
          IF (read_status /= 0) THEN
             CALL finish('init_lctlib','No keyword NLCT found in land cover library file '//TRIM(lctlib_file_name))
          END IF
          ! Look for comment
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
          IF(length== 0) CYCLE                 ! Line is empty
          
          pos = SCAN(line,blank_set)                   ! Position of first blank character
          IF(pos == 0)   CALL finish('init_lctlib',"Wrong syntax in lctlib definitions")
          
          READ(line(1:pos-1),'(A)') key
          
          IF(tolower(TRIM(key)) == 'nlct') THEN
             READ(line(pos:length),*,IOSTAT=read_status) nlct
             IF (read_status /= 0) THEN
                CALL finish('init_lctlib','Could not read number of landcover types (keyword: NLCT) from '//TRIM(lctlib_file_name))
             END IF
             EXIT ! found number of landcover types in landcover library file --- continue after loop
          END IF
       END DO

    END IF

    IF (p_parallel) CALL p_bcast(nlct, p_io)
    lctlib%nlct = nlct

    ALLOCATE(lctlib%LctClass(nlct))
    ALLOCATE(lctlib%LctName(nlct))
    ALLOCATE(lctlib%NaturalVegFlag(nlct))
    ALLOCATE(lctlib%ForestFlag(nlct))
    ALLOCATE(lctlib%CropFlag(nlct))
    ALLOCATE(lctlib%LakeFlag(nlct))
    ALLOCATE(lctlib%GlacierFlag(nlct))
    ALLOCATE(lctlib%BareSoilFlag(nlct))
    ALLOCATE(lctlib%AlbedoSnowVisMin(nlct))
    ALLOCATE(lctlib%AlbedoSnowVisMax(nlct))
    ALLOCATE(lctlib%AlbedoSnowNirMin(nlct))
    ALLOCATE(lctlib%AlbedoSnowNirMax(nlct))
    ALLOCATE(lctlib%AlbedoSnowMin(nlct))
    ALLOCATE(lctlib%AlbedoSnowMax(nlct))
    ALLOCATE(lctlib%AlbedoCanopyVIS(nlct))
    ALLOCATE(lctlib%AlbedoCanopyNIR(nlct))
    ALLOCATE(lctLib%nitrogenScalingFlag(nlct))
    ALLOCATE(lctlib%C4flag(nlct))
    ALLOCATE(lctlib%CarboxRate(nlct))
    ALLOCATE(lctlib%ETransport(nlct))
    ALLOCATE(lctlib%VegHeight(nlct))
    ALLOCATE(lctlib%VegRoughness(nlct))
    ALLOCATE(lctlib%PhenologyType(nlct))
    ALLOCATE(lctlib%CanopyResistanceMin(nlct))
    ALLOCATE(lctlib%MaxLAI(nlct))
    ALLOCATE(lctlib%StemArea(nlct))
    ALLOCATE(lctlib%specificLeafArea_C(nlct))
    ALLOCATE(lctlib%frac_npp_2_woodPool(nlct))
    ALLOCATE(lctlib%frac_npp_2_reservePool(nlct))
    ALLOCATE(lctlib%tau_Cpool_litter_leaf(nlct))
    ALLOCATE(lctlib%tau_Cpool_litter_wood(nlct))
    ALLOCATE(lctlib%LAI_shed_constant(nlct))
    ALLOCATE(lctlib%frac_C_fast2atmos(nlct))
    ALLOCATE(lctlib%Max_C_content_woods(nlct))
    IF (p_parallel_io) THEN

       REWIND(unit=lctlib_file_unit) ! Go back to beginning of landcover library file

       ALLOCATE(itmp(nlct)) 

       DO
          READ(lctlib_file_unit,'(A256)', IOSTAT=read_status) line
          IF (read_status /= 0) EXIT                   ! Finished reading

          ! Look for comment
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
          IF(length== 0) CYCLE                 ! Line is empty
          
          pos = SCAN(line,blank_set)                   ! Position of first blank character
          IF(pos == 0)   CALL finish('init_lctlib',"Wrong syntax in lctlib definitions")
          
          READ(line(1:pos-1),'(A)') key
          key = tolower(key)

          IF(TRIM(key) .EQ. "nlct") CYCLE ! nlct already read above

          SELECT CASE (TRIM(key))
          CASE ('lctclass')                                           ! Landcover class (including Plant Functional Types)
             READ(line(pos:length),*) lctlib%LctClass(1:nlct)
             exist_LctClass = .TRUE.
          CASE ('lctname')                                           ! Name of landcover type 
             READ(line(pos:length),*) lctlib%LctName(1:nlct)
             exist_LctName = .TRUE.
          CASE ('naturalvegflag')                                     ! Indicates natural vegetation as landcover
             READ(line(pos:length),*) itmp(1:nlct)
             WHERE(itmp(:) == 0) 
                lctlib%naturalVegFlag(:) = .FALSE.
             ELSEWHERE
                lctlib%naturalVegFlag(:) = .TRUE.
             END WHERE
             exist_NaturalVegFlag = .TRUE.
          CASE ('forestflag')                                     ! Indicates natural vegetation as laaandcover
             READ(line(pos:length),*) itmp(1:nlct)
             WHERE(itmp(:) == 0) 
                lctlib%ForestFlag(:) = .FALSE.
             ELSEWHERE
                lctlib%ForestFlag(:) = .TRUE.
             END WHERE
             exist_ForestFlag = .TRUE.
          CASE ('cropflag')                                     ! Indicates natural vegetation as landcover
             READ(line(pos:length),*) itmp(1:nlct)
             WHERE(itmp(:) == 0) 
                lctlib%cropFlag(:) = .FALSE.
             ELSEWHERE
                lctlib%cropFlag(:) = .TRUE.
             END WHERE
             exist_CropFlag = .TRUE.
          CASE ('lakeflag')                                           ! Indicates whether vegetation class signifies lake or ocean
             READ(line(pos:length),*) itmp(1:nlct)
             WHERE(itmp(:) == 0) 
                lctlib%LakeFlag(:) = .FALSE.
             ELSEWHERE
                lctlib%LakeFlag(:) = .TRUE.
             END WHERE
             exist_LakeFlag = .TRUE.
          CASE ('glacierflag')                                           ! Indicates whether vegetation class signifies glacier
             READ(line(pos:length),*) itmp(1:nlct)
             WHERE(itmp(:) == 0) 
                lctlib%GlacierFlag(:) = .FALSE.
             ELSEWHERE
                lctlib%GlacierFlag(:) = .TRUE.
             END WHERE
             exist_GlacierFlag = .TRUE.
          CASE ('baresoilflag')                                           ! Indicates whether vegetation class signifies bare soil
             READ(line(pos:length),*) itmp(1:nlct)
             WHERE(itmp(:) == 0) 
                lctlib%BareSoilFlag(:) = .FALSE.
             ELSEWHERE
                lctlib%BareSoilFlag(:) = .TRUE.
             END WHERE
             exist_BareSoilFlag = .TRUE.
          CASE ('albedosnowvismin')
             READ(line(pos:length),*) lctlib%AlbedoSnowVisMin(1:nlct)
             exist_AlbedoSnowVisMin = .TRUE.
          CASE ('albedosnowvismax')
             READ(line(pos:length),*) lctlib%AlbedoSnowVisMax(1:nlct)
             exist_AlbedoSnowVisMax = .TRUE.
          CASE ('albedosnownirmin')
             READ(line(pos:length),*) lctlib%AlbedoSnowNirMin(1:nlct)
             exist_AlbedoSnowNirMin = .TRUE.
          CASE ('albedosnownirmax')
             READ(line(pos:length),*) lctlib%AlbedoSnowNirMax(1:nlct)
             exist_AlbedoSnowNirMax = .TRUE.
          CASE ('albedosnowmin')
             READ(line(pos:length),*) lctlib%AlbedoSnowMin(1:nlct)
             exist_AlbedoSnowMin = .TRUE.
          CASE ('albedosnowmax')
             READ(line(pos:length),*) lctlib%AlbedoSnowMax(1:nlct)
             exist_AlbedoSnowMax = .TRUE.
          CASE ('albedocanopyvis')
             READ(line(pos:length),*) lctlib%AlbedoCanopyVIS(1:nlct)
             exist_AlbedoCanopyVIS = .TRUE.
          CASE ('albedocanopynir')
             READ(line(pos:length),*) lctlib%AlbedoCanopyNIR(1:nlct)
             exist_AlbedoCanopyNIR = .TRUE.
          CASE ('nitrogenscalingflag')        !Whether nitrogen scaling should be accounted for (.false. input as 0 and .true. as 1)
             READ(line(pos:length),*) itmp(1:nlct)
             WHERE(itmp(:) == 0) 
                lctlib%nitrogenScalingFlag(:) = .FALSE.
             ELSEWHERE
                lctlib%nitrogenScalingFlag(:) = .TRUE.
             END WHERE
             exist_nitrogenScalingFlag  = .TRUE.
          CASE ('c4flag')                              ! Photosynthetic pathway (C4: .true.; C3: .false.)
             READ(line(pos:length),*) itmp(1:nlct)
             WHERE(itmp(:) == 1) 
   !             lctlib%nitrogenScalingFlag(:) = C4
                lctlib%c4flag(:) = .TRUE.
             ELSEWHERE
   !             lctlib%nitrogenScalingFlag(:) = C3
                lctlib%c4flag(:) = .FALSE.
             END WHERE
             exist_C4flag = .TRUE.
          CASE ('carboxrate')                                         ! Carbox rate
             READ(line(pos:length),*) lctlib%CarboxRate(1:nlct)
             exist_CarboxRate = .TRUE.
          CASE ('etransport')                                         ! E-transport
             READ(line(pos:length),*) lctLib%ETransport(1:nlct)
           exist_ETransport = .true.
          CASE ('vegheight')                                          ! typical vegetation height
             READ(line(pos:length),*) lctLib%VegHeight(1:nlct)
             exist_VegHeight = .true.
          CASE ('vegroughness')                                       ! typical vegetation roughness length
             READ(line(pos:length),*) lctLib%VegRoughness(1:nlct)
             exist_VegRoughness = .true.
          CASE ('phenologytype')                                      ! Phenology type
             READ(line(pos:length),*) lctlib%PhenologyType(1:nlct)
             exist_PhenologyType = .TRUE.
          CASE ('canopyresistancemin')                                ! Minimum canopy resistance
             READ(line(pos:length),*) lctlib%CanopyResistanceMin(1:nlct)
             exist_CanopyResistanceMin = .TRUE.
         CASE ('maxlai')                                              ! Maximum LAI for LogroP
             READ(line(pos:length),*) lctlib%MaxLAI(1:nlct)
             exist_MaxLAI = .TRUE.
         CASE ('stemarea')                                            ! Area of stems and branches
             READ(line(pos:length),*) lctlib%StemArea(1:nlct)
             exist_StemArea = .TRUE.
         CASE ('specificleafarea_c')                                  ! Carbon content per leaf area for LogroP
             READ(line(pos:length),*) lctlib%specificLeafArea_C(1:nlct)
             exist_specificLeafArea_C = .TRUE.
         CASE ('frac_npp_2_woodpool')                                 ! fraction of NPP directed to wood pool (see cbalance)
             READ(line(pos:length),*) lctlib%frac_npp_2_woodPool(1:nlct)
             exist_frac_npp_2_woodPool = .TRUE.
         CASE ('frac_npp_2_reservepool')                              ! fraction of NPP directed to reserve pool (see cbalance)
             READ(line(pos:length),*) lctlib%frac_npp_2_reservePool(1:nlct)
             exist_frac_npp_2_reservePool = .TRUE.
         CASE ('tau_cpool_litter_leaf')                               ! Time constant by which leaf litter pool is depreciated  [days]
             READ(line(pos:length),*) lctlib%tau_Cpool_litter_leaf(1:nlct)
             exist_tau_Cpool_litter_leaf = .TRUE.
         CASE ('tau_cpool_litter_wood')                               ! Time constant by which woody litter pool is depreciated  [days]
             READ(line(pos:length),*) lctlib%tau_Cpool_litter_wood(1:nlct)
             exist_tau_Cpool_litter_wood = .TRUE.
         CASE ('lai_shed_constant')                                   ! Time constant by which leaves are constantly shedded [days-1]
             READ(line(pos:length),*) lctlib%LAI_shed_constant(1:nlct)
             exist_LAI_shed_constant = .TRUE.
         CASE ('frac_c_fast2atmos')                                   ! fraction of loss of fast pool emitted to atmosphere (see cbalance)
             READ(line(pos:length),*) lctlib%frac_C_fast2atmos(1:nlct)
             exist_frac_C_fast2atmos = .TRUE.
         CASE ('max_c_content_woods')                                 ! Maximum carbon content of wood pool (see cbalance)
             READ(line(pos:length),*) lctlib%Max_C_content_woods(1:nlct)
             exist_Max_C_content_woods = .TRUE.
          CASE default
             ! nothing to do
          END SELECT
       END DO
       DEALLOCATE(itmp)

       IF(.NOT. exist_LctClass) &
            CALL finish('init_lctlib','No data for LctClass found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_LctName) &
!            CALL finish('init_lctlib','No data for LctName found in '//TRIM(lctlib_file_name))
            ! Set default
            lctlib%LctName = "undefined"
       IF(.NOT. exist_NaturalVegFlag) &
            CALL finish('init_lctlib','No data for NaturalVegFlag found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_ForestFlag) &
            CALL finish('init_lctlib','No data for ForestFlag found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_CropFlag) &
            CALL finish('init_lctlib','No data for CropFlag found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_LakeFlag) &
            CALL finish('init_lctlib','No data for LakeFlag found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_GlacierFlag) &
            CALL finish('init_lctlib','No data for GlacierFlag found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_BareSoilFlag) &
            CALL finish('init_lctlib','No data for BareSoilFlag found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowVisMin) &
            CALL finish('init_lctlib','No data for AlbedoSnowVisMin found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowVisMax) &
            CALL finish('init_lctlib','No data for AlbedoSnowVisMax found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowNirMin) &
            CALL finish('init_lctlib','No data for AlbedoSnowNirMin found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowNirMax) &
            CALL finish('init_lctlib','No data for AlbedoSnowNirMax found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowMin) &
            CALL finish('init_lctlib','No data for AlbedoSnowMin found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowMax) &
            CALL finish('init_lctlib','No data for AlbedoSnowMax found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoCanopyVIS) &
            CALL finish('init_lctlib','No data for AlbedoCanopyVIS found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoCanopyNIR) &
            CALL finish('init_lctlib','No data for albedoCanopyNIR found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_nitrogenScalingFlag) &
            CALL finish('init_lctLib','No data for NitrogenScalingFlag found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_C4flag) &
            CALL finish('init_lctLib','No data for C4flag found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_CarboxRate) &
            CALL finish('init_lctLib','No data for CarboxRate found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_ETransport) &
            CALL finish('init_lctLib','No data for ETransport found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_VegHeight) &
            CALL finish('init_lctLib','No data for VegHeight found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_VegRoughness) &
            CALL finish('init_lctLib','No data for VegRoughness found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_PhenologyType) &
            CALL finish('init_lctlib','No data for PhenologyType found in '//TRIM(lctlib_file_name))
       IF (.NOT. exist_CanopyResistanceMin) &
            lctlib%CanopyResistanceMin = -1.0_dp
       IF(.NOT. exist_MaxLAI) &
            CALL finish('init_lctlib','No data for MaxLAI found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_StemArea) &
            CALL finish('init_lctlib','No data for StemArea found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_specificLeafArea_C) &
            CALL finish('init_lctlib','No data for specificLeafArea_C found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_frac_npp_2_woodPool) &
            CALL finish('init_lctlib','No data for frac_npp_2_woodPool found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_frac_npp_2_reservePool) &
            CALL finish('init_lctlib','No data for frac_npp_2_reservePool found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_tau_Cpool_litter_leaf) &
            CALL finish('init_lctlib','No data for tau_Cpool_litter_leaf found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_tau_Cpool_litter_wood) &
            CALL finish('init_lctlib','No data for tau_Cpool_litter_wood found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_LAI_shed_constant) &
            CALL finish('init_lctlib','No data for LAI_shed_constant found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_frac_C_fast2atmos) &
            CALL finish('init_lctlib','No data for frac_C_fast2atmos found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_Max_C_content_woods) &
            CALL finish('init_lctlib','No data for Max_C_content_woods found in '//TRIM(lctlib_file_name))
       CLOSE(unit=lctlib_file_unit) !! Closing access to lctlib file

       npft = 0
       !------------------------------------------------------------------------------------------
       ! kalle, 151004
       ! Set carbox rate to  [Mol(CO2)/m^2/s] 
       ! SET e-transport rate to [Mol(CO2)/m^2/s]       
       DO i=1,nlct
         lctLib%ETransport(i)=lctLib%ETransport(i) * 1.e-06_dp
         lctlib%CarboxRate(i)=lctlib%CarboxRate(i) * 1.e-06_dp
       END DO
       !------------------------------------------------------------------------------------------
       DO i=1,nlct
          IF(lctlib%naturalVegFlag(i) .OR. lctlib%cropFlag(i)) THEN

             npft = npft+1 ! count PFTs

!!$             ! Check consistency of PFT flags
!!$             IF(lctlib%naturalVegFlag(i) .AND. lctlib%cropFlag(i)) THEN
!!$                CALL finish('init_lctlib',&
!!$                     'Inconsistency (NaturalVegFlag and CropFlag) in file '//trim(lctlib_file_name))
!!$             END IF
          
             ! Check consistency with lake flag
             IF(lctlib%LakeFlag(i)) THEN
                CALL finish('init_lctlib',&
                     'Inconsistent entries for PFT flags and LakeFlag in landcover library file '//trim(lctlib_file_name))
             END IF

             ! Check consistency with glacier flag
             IF(lctlib%GlacierFlag(i)) THEN
                CALL finish('init_lctlib',&
                     'Inconsistent entries for PFT flags and GlacierFlag in landcover library file '//trim(lctlib_file_name))
             END IF

             ! Check consistency with PhenologyType
             IF(lctlib%naturalVegFlag(i) .AND. & 
                  (lctlib%PhenologyType(i) .le.0 .or. lctlib%PhenologyType(i) > 5) ) then
                CALL finish('init_lctlib',&
                     'Inconsistent entries for NaturalVegFlag and PhenologyType in landcover library file '//trim(lctlib_file_name))
             END IF

          END IF
       END DO

    END IF


    IF(p_parallel) THEN
       CALL p_bcast(lctlib%LctClass, p_io)
       DO i=1,nlct
          CALL p_bcast(lctlib%LctName(i), p_io)
       ENDDO
       CALL p_bcast(lctlib%NaturalVegFlag,p_io)
       CALL p_bcast(lctlib%ForestFlag,p_io)
       CALL p_bcast(lctlib%CropFlag,p_io)
       CALL p_bcast(lctlib%LakeFlag, p_io)
       CALL p_bcast(lctlib%GlacierFlag, p_io)
       CALL p_bcast(lctlib%AlbedoSnowVisMin, p_io)
       CALL p_bcast(lctlib%AlbedoSnowVisMax, p_io)
       CALL p_bcast(lctlib%AlbedoSnowNirMin, p_io)
       CALL p_bcast(lctlib%AlbedoSnowNirMax, p_io)
       CALL p_bcast(lctlib%AlbedoSnowMin, p_io)
       CALL p_bcast(lctlib%AlbedoSnowMax, p_io)
       CALL p_bcast(lctlib%AlbedoCanopyVIS, p_io)
       CALL p_bcast(lctlib%AlbedoCanopyNIR, p_io)
       CALL p_bcast(lctlib%NitrogenScalingFlag,p_io)
       CALL p_bcast(lctlib%C4flag, p_io)
       CALL p_bcast(lctlib%CarboxRate, p_io)
       CALL p_bcast(lctlib%ETransport, p_io)
       CALL p_bcast(lctlib%VegHeight, p_io)
       CALL p_bcast(lctlib%VegRoughness, p_io)
       CALL p_bcast(lctlib%PhenologyType, p_io)
       CALL p_bcast(lctlib%CanopyResistanceMin, p_io)
       CALL p_bcast(lctlib%MaxLAI, p_io)
       CALL p_bcast(lctlib%StemArea, p_io)
       CALL p_bcast(lctlib%specificLeafArea_C, p_io)
       CALL p_bcast(lctlib%frac_npp_2_woodPool, p_io)
       CALL p_bcast(lctlib%frac_npp_2_reservePool, p_io)
       CALL p_bcast(lctlib%tau_Cpool_litter_leaf, p_io)
       CALL p_bcast(lctlib%tau_Cpool_litter_wood, p_io)
       CALL p_bcast(lctlib%LAI_shed_constant, p_io)
       CALL p_bcast(lctlib%frac_C_fast2atmos, p_io)
       CALL p_bcast(lctlib%Max_C_content_woods, p_io)
       CALL p_bcast(npft, p_io)
    END IF

    IF(npft>nlct) CALL finish('init_lctlib','Christian did bad programming!')

    lctlib%npft = npft

  END SUBROUTINE init_lctlib
  
END MODULE mo_jsbach_lctlib
