!+ Type definitions and common functions for JSBACH
! 
MODULE mo_jsbach

  ! 
  ! Description: 
  !   <Say what this module is for> 
  ! 
  ! Current Code Owner: jsbach_admin
  ! 
  ! History: 
  !  
  ! Version   Date     Comment 
  ! -------   ----     ------- 
  ! 4.0.3     01/06/28 Original Fortran 90 code. Reiner Schnur
  ! 
  ! Code Description: 
  !   Language:           Fortran 90. 
  !   Software Standards: "European Standards for Writing and  
  !     Documenting Exchangeable Fortran 90 Code". 
  ! 
  ! Modules used: 
  ! 
  USE mo_time_control, ONLY: lresume, dt_start, dt_resume, dt_stop, delta_time, no_cycles, &
                             no_days, no_steps, putdata, putrerun, trigfiles, ldebugev, &
                             p_bcast_event
  USE mo_time_event,   ONLY: io_time_event,                                                &
                             TRIG_LAST, TRIG_FIRST, TRIG_NONE,                             &
                             TIME_INC_SECONDS, TIME_INC_MINUTES,                           &
                             TIME_INC_HOURS , TIME_INC_DAYS,                               &
                             TIME_INC_MONTHS, TIME_INC_YEARS
  USE mo_control,      ONLY: ltimer, ldebugio
  USE mo_filename,     ONLY: name_limit, path_limit, out_expname, GRIB, NETCDF
  USE mo_netcdf,       ONLY: nf_fill_real, nf_max_name
  USE mo_mpi,          ONLY: p_parallel, p_parallel_io, p_bcast, p_io
  USE mo_kind,         ONLY: dp
  USE mo_exception,    ONLY: finish, message, message_text, int2string, real2string
  USE mo_doctor,       ONLY: nout

  IMPLICIT NONE

  ! Global (i.e. public) Declarations: 
  ! Global Type Definitions: 

  !! Structure to hold options that define the model structure and switch on/off
  !! certain model operations. They are defined from namelist parameters. 
  TYPE options_type
     INTEGER      :: ntiles                 !! Number of sub-grid tiles on land surface
                                            !!   Has to be equal to the <ntiles> dimension in jsbach initial file 
                                            !!   and smaller than number of land cover types <nlct>
     LOGICAL      :: Standalone             !! Type of model run (standalone or coupled)
     CHARACTER(name_limit) :: Experiment    !! Label for experiment
     CHARACTER(len=10) :: Coupling          !! Type of coupling between LSS and atmosphere 
                                            !! (Choices: "explicit","semi","implicit")
     !
     CHARACTER(len=10) :: LSS               !! Which LSS to use?
     LOGICAL      :: UseVic                 !! Use VIC model?
     LOGICAL      :: UseEchamLand           !! Use old ECHAM land surface scheme?
     LOGICAL      :: UseBethyLand           !! Use old BETHY land surface scheme?
     LOGICAL      :: UseBethy               !! Use BETHY model?
     LOGICAL      :: UsePhenology           !! Use phenology model?
     LOGICAL      :: UseAlbedo              !! Use albedo model?
     LOGICAL      :: UseDynveg              !! Use dynamic vegetation DYNVEG?
     LOGICAL      :: LCChange               !! Read landcover maps from external files?
     LOGICAL      :: ReadCoverFract         !! Read cover fractions from initial file instead of restart file in restarted runs
     !
     INTEGER      :: FileType               !! File type for output files (NETCDF or GRIB)
     LOGICAL      :: OutputModelState       !! Output full model state stream?
     CHARACTER(path_limit) :: RestartPrefix !! Prefix for construction of restart files
     LOGICAL      :: Timer                  !! Use Timer?
     LOGICAL      :: ResetTime              !! Allow model to overwrite start of restart files with time of first call?
     CHARACTER(nf_max_name)  :: GridFile    !! File containing grid information
     CHARACTER(nf_max_name)  :: LctlibFile  !! Name of the land cover library file
     CHARACTER(nf_max_name)  :: VegFile     !! File containig initial data for the vegetation 
     CHARACTER(nf_max_name)  :: SurfFile    !! File containig initial data for the surface
     CHARACTER(nf_max_name)  :: SoilFile    !! File containig initial data for the soil
     !
     ! The following measurement heights are only used by the standalone model for the forcing data. In coupled mode, they
     ! are not changed from their default values of 0. They are added to the geopotential height passed to the interface and
     ! describing the height of the lowest atmospheric level where wind, temperature and humidity are taken from in the coupled 
     ! mode. In order for this to work in standalone mode, the interface has to be passed the geopotential of the surface
     ! elevation, i.e. elevation * Gravity. That is:
     ! Coupled mode: Pass geopotential of lowest atmosphere level and set the following three values to zero
     ! Standalone:   Pass geopotential of elevation to interface and set the three values to there heights above surface
     REAL(dp)     :: HeightTemperature      !! Height above surface at which temperature measurements were taken [m]
     REAL(dp)     :: HeightWind             !! Height above surface at which wind measurements were taken [m]
     REAL(dp)     :: HeightHumidity         !! Height above surface at which humidity measurements were taken
     !
  END TYPE options_type

  !! Global Scalars: 
  LOGICAL, SAVE  :: jsbach_FirstCall   = .TRUE. !! True only on first call
  LOGICAL, SAVE  :: jsbach_NewTimeStep = .TRUE. !! TRUE iff start of new time step
  LOGICAL, SAVE  :: debug                       !! Generate debug output in standard output?
  LOGICAL, SAVE  :: test_stream                 !! Write test stream in mo_test module?
  LOGICAL, SAVE  :: lpost_echam                 !! Write JSBACH output variables even if they are part of the ECHAM output stream
  REAL(dp), SAVE :: missing_value               !! Missing value used in model output
  INTEGER, SAVE  :: nml_unit                    !! Unit of the jsbach namelist

  LOGICAL :: module_configured  = .FALSE.

CONTAINS 

  !!+ Initialize configuration options and parameters for JSBACH
  SUBROUTINE jsbach_config(options)

    USE mo_util_string, ONLY: tolower, toupper
    USE mo_filename,    ONLY: out_filetype, find_next_free_unit
    USE mo_namelist,    ONLY: POSITIONED, open_nml, position_nml

    !! Description:
    !! This subroutine reads the namelists jsbach_ctl (and jsbalone_ctl for stand alone jsbach runs)
    !! into the global structure <options>

    !! Called from the interface routine *jsbach_init*

    TYPE(options_type), INTENT(inout) :: options

    !! Namelist parameters

    INTEGER                 :: ntiles           ! number of tiles
    LOGICAL                 :: standalone       ! Type of model run 
                                                !    true: stand-alone jsbach run
                                                !    false: jsbach driven by an atmosphere model
    CHARACTER(len=10)       :: coupling         ! Type of coupling: explicit, semi, implicit
    CHARACTER(len=10)       :: LSS              ! Land surface sceme: ECHAM
    LOGICAL                 :: use_bethy        ! Use BETHY model (photosynthesis, respiration)
    LOGICAL                 :: use_phenology    ! Calculate LAI using the phenology module
    LOGICAL                 :: use_albedo       ! Calculate albedo depending on vegetaion
    LOGICAL                 :: use_dynveg       ! Use the dynamic vegetation module
    LOGICAL                 :: lcc              ! Read land cover maps from external files?
    CHARACTER(len=10)       :: file_type        ! Output file format
    LOGICAL                 :: out_state        ! write the jsbach stream
    LOGICAL                 :: test_stream      ! Additional stream for model testing
    LOGICAL                 :: read_cover_fract ! read cover fractions from initial rather than restart file
    CHARACTER(nf_max_name)  :: grid_file        ! File containing grid information
    CHARACTER(nf_max_name)  :: lctlib_file      ! Name of the land cover library file
    CHARACTER(nf_max_name)  :: veg_file         ! File containig initial data for the vegetation 
    CHARACTER(nf_max_name)  :: surf_file        ! File containig initial data for the surface
    CHARACTER(nf_max_name)  :: soil_file        ! File containig initial data for the surface

    !! other local parameters
    INTEGER                 :: read_status
    

    INCLUDE 'jsbach_ctl.inc'
    INCLUDE 'jsbalone_ctl.inc'


    IF (p_parallel_io) THEN

       !! Open the jsbach namelist file

       nml_unit = find_next_free_unit (51,100)
       CALL open_nml ('namelist.jsbach', unit=nml_unit)


       !! Read namelist jsbach_ctl

       !!  define default values
       ntiles = -1
       standalone = .TRUE.
       coupling = 'implicit'
       lss = 'ECHAM'
       use_bethy = .FALSE.
       use_phenology = .FALSE.
       use_albedo = .FALSE.
       use_dynveg = .FALSE.
       lcc = .FALSE.
       file_type = 'GRIB'
       out_state = .TRUE.
       lpost_echam = .FALSE.
       missing_value = NF_FILL_REAL
       test_stream = .FALSE.
       debug = .FALSE.
       grid_file = 'jsbach.nc'
       lctlib_file = 'lctlib.def'
       veg_file = 'jsbach.nc'
       surf_file = 'jsbach.nc'
       soil_file = 'jsbach.nc'
       read_cover_fract = .FALSE.
       
       CALL position_nml ('JSBACH_CTL', status=read_status)
       SELECT CASE (read_status)
       CASE (POSITIONED)
          READ (nml_unit, jsbach_ctl)
          CALL message('jsbach_config', 'Namelist JSBACH_CTL: ')
          WRITE(nout, jsbach_ctl)
       END SELECT

       options%Standalone = standalone
       WRITE(*,*) 'Type of model run: standalone = ', options%Standalone
       options%Coupling = TRIM(coupling)
       WRITE(*,*) 'Type of coupling: ',options%Coupling
          
       options%LSS = TRIM(lss)
       WRITE(*,*) 'Using LSS: ', options%LSS
       options%UseVic = .FALSE.
       options%UseEchamLand = .FALSE.
       options%UseBethyLand = .FALSE.
       SELECT CASE(tolower(options%LSS))
       CASE ('vic')
          options%UseVic = .TRUE.
       CASE ('echam')
          options%UseEchamLand = .TRUE.
       CASE ('bethy')
          options%UseBethyLand = .TRUE.
       CASE default
          CALL finish('jsbach_config','You have to use some land surface scheme')
       END SELECT

       options%UseBethy = use_bethy
       WRITE (message_text,*) 'Using BETHY model: ', options%UseBethy
       CALL message('jsbach_config', message_text)

       options%UsePhenology = use_phenology
       WRITE (message_text,*) 'Using phenology model: ', options%UsePhenology
       CALL message('jsbach_config', message_text)

       IF (options%UsePhenology .AND. .NOT. options%UseBethy) &
            CALL finish('jsbach_config','Phenology model can only be used together with BETHY')

       options%UseAlbedo = use_albedo
       WRITE (message_text,*) 'Using albedo model: ', options%UseAlbedo
       CALL message('jsbach_config', message_text)

       options%UseDynveg = use_dynveg
       WRITE (message_text,*) 'Using dynamic vegetation: ', options%UseDynveg
       CALL message('jsbach_config', message_text)

       options%LCChange = lcc
       WRITE (message_text,*) 'Read maps for land cover change: ', options%LCChange
       CALL message('jsbach_config', message_text)

       options%ReadCoverFract = read_cover_fract
       WRITE (message_text,*) 'Read cover fractions from initial rather than restart file: ', read_cover_fract
       CALL message('jsbach_config', message_text)

       options%ntiles = ntiles
       WRITE (message_text,*) 'Number of tiles: ', ntiles
       CALL message('jsbach_config', message_text)

       file_type = toupper(file_type)
       SELECT CASE(TRIM(file_type))
       CASE ('GRIB')
          options%FileType = GRIB
       CASE ('NETCDF')
          options%FileType = NETCDF
       CASE default
          CALL finish('jsbach_config','Wrong FILE_TYPE')
       END SELECT
       WRITE (message_text,*) 'JSBACH output file type: ', file_type
       CALL message('jsbach_config', message_text)

       IF (options%Standalone) out_filetype = options%FileType

       options%OutputModelState = out_state
       WRITE (message_text,*) 'Output full model state', out_state
       CALL message('jsbach_config', message_text)

       IF (.NOT. lpost_echam) THEN
          CALL message('jsbach_config',' Output variables defined in ECHAM and JSBACH streams written by ECHAM only.')
       END IF
       
       WRITE (message_text,*) 'Fill value for output files', missing_value
       CALL message('jsbach_config', message_text)

       WRITE (message_text,*) 'Test Stream: ', test_stream
       CALL message('jsbach_config', message_text)

       options%GridFile = grid_file
       WRITE (message_text,*) 'Grid information is read from file: ', grid_file
       CALL message('jsbach_config', message_text)

       options%VegFile = veg_file
       WRITE (message_text,*) 'Vegetation data is read from file: ', veg_file
       CALL message('jsbach_config', message_text)

       options%SurfFile = surf_file
       WRITE (message_text,*) 'Surface data is read from file: ', surf_file
       CALL message('jsbach_config', message_text)

       options%SoilFile = soil_file
       WRITE (message_text,*) 'Soil data is read from file: ', soil_file
       CALL message('jsbach_config', message_text)

       options%LctlibFile = lctlib_file
       WRITE (message_text,*) 'Land cover type library: ', lctlib_file
       CALL message('jsbach_config', message_text)

       WRITE (message_text,*) 'DEBUG:', debug
       CALL message('jsbach_config', message_text)
    ENDIF

    ! Broadcast to processors
    IF (p_parallel) THEN
       CALL p_bcast(options%ntiles, p_io)
       CALL p_bcast(options%Standalone, p_io)
       CALL p_bcast(options%Coupling, p_io)
       CALL p_bcast(options%LSS, p_io)
       CALL p_bcast(options%UseVic, p_io)
       CALL p_bcast(options%UseEchamLand, p_io)
       CALL p_bcast(options%UseBethyLand, p_io)
       CALL p_bcast(options%UseBethy, p_io)
       CALL p_bcast(options%UsePhenology, p_io)
       CALL p_bcast(options%UseAlbedo, p_io)
       CALL p_bcast(options%UseDynveg, p_io)
       CALL p_bcast(options%LCChange, p_io)
       CALL p_bcast(options%ReadCoverFract, p_io)
       CALL p_bcast(options%FileType, p_io)
       CALL p_bcast(options%OutputModelState, p_io)
       CALL p_bcast(lpost_echam, p_io)
       CALL p_bcast(missing_value, p_io)
       CALL p_bcast(test_stream, p_io)
       CALL p_bcast(debug, p_io)
       CALL p_bcast(options%GridFile, p_io)
       CALL p_bcast(options%VegFile, p_io)
       CALL p_bcast(options%SurfFile, p_io)
       CALL p_bcast(options%SoilFile, p_io)
       CALL p_bcast(options%LctlibFile, p_io)
    ENDIF

    IF (options%standalone) THEN
       IF (p_parallel_io) THEN

          !! Read namelist jsbalone_ctl

          !!  define default values
          out_expname = 'xxxxxx'
          lresume = .FALSE.
          no_cycles = 1
          dt_start(:) = 0
          dt_resume(:) = 0
          dt_stop(:) = 0
          no_days = -1
          no_steps = -1
          delta_time = 0
          ltimer = .FALSE.
          putdata  = io_time_event(1, TIME_INC_DAYS, TRIG_FIRST, 0)
          putrerun = io_time_event(1, TIME_INC_MONTHS, TRIG_LAST, 0)
          trigfiles = io_time_event(1, TIME_INC_MONTHS, TRIG_FIRST, 0)

          CALL position_nml ('JSBALONE_CTL', status=read_status)
          SELECT CASE (read_status)
          CASE (POSITIONED)
             READ (nml_unit, jsbalone_ctl)
             CALL message('jsbach_config', 'Namelist JSBALONE_CTL: ')
             WRITE(nout, jsbalone_ctl)
          END SELECT

          WRITE (message_text,*) 'Experiment name', TRIM(out_expname)
          CALL message('jsbach_config', message_text)
          options%Experiment = out_expname
          options%RestartPrefix = 'rerun_' // TRIM(out_expname) ! This is the echam prefix (hardcoded)
          
          WRITE (message_text,*) 'Initializing model from restart files', lresume
          CALL message('jsbach_config', message_text)

          WRITE (message_text,*) 'Number of restart cycles', no_cycles
          CALL message('jsbach_config', message_text)

          WRITE (message_text,*) 'dt_start:', dt_start
          CALL message('jsbach_config', message_text)

          WRITE (message_text,*) 'dt_stop:', dt_stop
          CALL message('jsbach_config', message_text)

          WRITE (message_text,*) 'restart periode:', putrerun
          CALL message('jsbach_config', message_text)

          WRITE (message_text,*) 'output periode:', putdata
          CALL message('jsbach_config', message_text)

          IF (no_days /= -1) THEN
             WRITE (message_text,*) 'Number of days of this run: ', no_days
             CALL message('jsbach_config', message_text)
          END IF

          IF (no_steps /= -1) THEN
             WRITE (message_text,*) 'Number of time steps in this run: ', no_steps
             CALL message('jsbach_config', message_text)
          END IF

          WRITE (message_text,*) 'Time step in seconds: ', delta_time
          CALL message('jsbach_config', message_text)

          options%Timer = ltimer
          WRITE (message_text,*) 'Check performances using Timer:', ltimer
          CALL message('jsbach_config', message_text)

          ldebugio = debug
          ldebugev = debug
       ENDIF

       IF (p_parallel) THEN
          CALL p_bcast(options%Experiment, p_io)
          CALL p_bcast(options%RestartPrefix, p_io)
          CALL p_bcast(options%Timer, p_io)
       ! The following are used from ECHAM5 modules and are (partly) duplicated in options, but need
       ! to be broadcast if running in standalone mode so that the modules can use them
          CALL p_bcast(lresume, p_io)
          CALL p_bcast(out_expname, p_io)
          CALL p_bcast(no_cycles, p_io)
          CALL p_bcast(dt_start, p_io)
          CALL p_bcast(dt_resume, p_io)
          CALL p_bcast(dt_stop, p_io)
          CALL p_bcast(no_days, p_io)
          CALL p_bcast(no_steps, p_io)
          CALL p_bcast(delta_time, p_io)
          CALL p_bcast(ltimer, p_io)
          CALL p_bcast(ldebugio, p_io)
          CALL p_bcast(ldebugev, p_io)
          CALL p_bcast_event(putdata, p_io)
          CALL p_bcast_event(trigfiles, p_io)
          CALL p_bcast_event(putrerun, p_io)
       ENDIF
    END IF  ! standalone

    module_configured = .TRUE.

  END SUBROUTINE jsbach_config
  !
  !----------------------------------------------------------------------------
  
END MODULE mo_jsbach

!Local Variables:
!mode: f90
!fill-column: 100
!End:
