!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                            !!
!!                          D R I V E R                       !!
!!                                                            !!
!!             for running the carbon pool model offline      !!
!!                                                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program cbalone_driver

  USE mo_mpi,                 ONLY: p_start
  USE mo_jsbach_lctlib,       ONLY: init_lctlib,lctlib_type
  USE mo_cbal_cpools,         ONLY: update_Cpools,relocate_carbon
  USE mo_kind,                ONLY: dp 
  USE mo_cbalone_io,          ONLY: getDailyData,writeSingleTimeStep,read_landcover_fractions
  USE mo_cbalone_memory,      ONLY: cbal_offline_type,grid_offline_type,vegetation_offline_type
  USE mo_cbalone_memory,      ONLY: landcover_change_type
  USE mo_cbalone_memory,      ONLY: initGrid,initVegetation,initializePools,initCbalance,initLandCoverChange
  USE mo_cbalone_memory,      ONLY: dynveg_offline_type, dynveg_clim_type
  USE mo_cbalone_dynveg,      ONLY: init_offline_dynveg, update_offline_dynveg
  USE mo_data_scaling,        ONLY: data_scaling_type,initDataScaling,scaleDailyData
  USE mo_exception,           ONLY: int2string,finish
  USE mo_util_string,         ONLY: tolower
  USE mo_netcdf,              ONLY: nf_max_name
  USE mo_jsbach_constants,    ONLY: molarMassCO2_kg
  USE mo_dynveg,              ONLY: config_dynveg, fpc_daily, potential_tree_fpc, desert_fraction, &
                                    fpc_to_cover_fract, dynveg_params_type, dynveg_options_type
  USE mo_filename,            ONLY: find_next_free_unit
  USE mo_namelist,            ONLY: open_nml, position_nml, POSITIONED
#ifdef NAG
  USE f90_unix_io,            ONLY: flush
#endif

  implicit none

  !! === PARAMETERS ====================================================================================
  
  ! --- Standardnames for input files ----

  character(len=*),parameter :: iniFileVegetation =  "jsbach.nc" !! File with information on grid and initial landcover
  character(nf_max_name),parameter :: lctlib_file = "lctlib.def"  !! File name of the Land cover library
  character(len=*),parameter :: outStemName =        "Cbalone"   !! Output file names are constructed as: outStemName.xxxxxx.yyyymm.nc
                                                                 !! .. where xxxxxx is the repetion number

  ! --- Parameters for landcover change

  real(dp),parameter :: frac_wood_2_atmos    = 0.6 !! Fraction of carbon from the wood pool to be released into the atmosphere on human landcover change
  real(dp),parameter :: frac_green_2_atmos   = 0.6 !! Fraction of carbon from the green pool to be released into the atmosphere on human landcover change
  real(dp),parameter :: frac_reserve_2_atmos = 0.6 !! Fraction of carbon from the reserve pool to be released into the atmosphere on human landcover chang

  ! --- Elements of filename for driver data (NPP,LAI,alpha,temperature)
  ! Note: It is assumed that the filenames have the form  experiment_yyyymm.specificName.nc
  !       The infornmation on experiment and specific name are taken from run.def

  character(len=128) :: experiment
  character(len=128) :: specificName

  ! --- start and end of integration (read in from namelist)

  integer   :: climate_yearstart   !! first year from external (climate) data
  integer   :: climate_yearend     !! last year  from external (climate) data
  integer   :: run_year_first !! First absolute year of run 
  integer   :: run_year_last  !! Last absolute year of run

  character(len=5) :: out_interval !! output interval

  character(len=1024) :: driver_data_path !! path where the driver data (LAI,NPP,alpha,temperature) are found

  ! --- Control of initialization of carbon pools (options read in from run.def)

  logical             :: read_cpools         !! true: Cpools are initialised from 'cpool_file_name' 
  character(len=1024) :: cpool_file_name     !! file name, including path, with initial data of the carbon pools

  ! --- Control of dynamic vegetation

  logical             :: use_dynveg           !! true if dynamic vegetation module is used
 
  ! --- structures -----------

  TYPE(lctlib_type)             :: lctlib
  TYPE(grid_offline_type)       :: grid
  TYPE(cbal_offline_type)       :: cbalance
  TYPE(vegetation_offline_type) :: vegetation
  TYPE(landcover_change_type)   :: LC_change
  TYPE(data_scaling_type)       :: data_scaling
  TYPE(dynveg_clim_type)        :: dynveg_clim
  TYPE(dynveg_offline_type)     :: dynveg
  TYPE(dynveg_params_type)      :: dynveg_params
  TYPE(dynveg_options_type)     :: dynveg_options

  ! --- allocatables

  real(dp),allocatable,dimension(:,:)  ::  specificLeafArea_C
  real(dp),allocatable,dimension(:,:)  ::  veg_fract_correction
  real(dp),allocatable,dimension(:,:)  ::  frac_npp_2_woodPool
  real(dp),allocatable,dimension(:,:)  ::  frac_npp_2_reservePool
  real(dp),allocatable,dimension(:,:)  ::  tau_Cpool_litter_leaf
  real(dp),allocatable,dimension(:,:)  ::  tau_Cpool_litter_wood
  real(dp),allocatable,dimension(:,:)  ::  LAI_shed_constant
  real(dp),allocatable,dimension(:,:)  ::  frac_C_fast2atmos
  real(dp),allocatable,dimension(:,:)  ::  Max_C_content_woods
  real(dp),allocatable,dimension(:,:)  ::  coverFractBox
  real(dp),allocatable,dimension(:,:)  ::  coverFract_new           !! The new landcover fractions read in each year
  real(dp),allocatable,dimension(:,:)  ::  coverFract_interpolated  !! The landcover fractions interpolated throughout the year
  real(dp),allocatable,dimension(:)    ::  carbon_2_atmos,carbon_2_fastSoilPool,carbon_2_slowSoilPool
  real(dp),allocatable,dimension(:)    ::  rock_fract    !! area fraction not suitable for vegetation (lakes, rocks, ...)
  real(dp),allocatable,dimension(:)    ::  veg_ratio_max_old
  real(dp),allocatable,dimension(:)    ::  veg_ratio_old_new
  real(dp),allocatable,dimension(:,:)  ::  cover_fract_old

  ! --- other locals

  integer             ::  itile,ihlp
  integer             ::  day,nday
  integer             ::  climate_year,month,number_of_days_in_year
  character(len=128)  ::  outFile
  character(len=128)  ::  inFile
  integer             ::  digits
  character(len=16)   ::  sHlp,sOut
  integer             ::  run_year !! counts years including repeats
  integer             ::  dayInYear !! actual day number in year, where Jan, 1, has day number 1
  integer             ::  noOfDaysInThisYear !! Total number of days in a year
  integer             ::  outputInterval
  character(len=1024) ::  LCC_fileName
  character(len=1024) ::  CO2_file,climatology_diff_file
  character(len=1024) ::  fpc_file_name
  integer             ::  ref_year_past,ref_year_recent !! Only used when CO2-scaling is switched on
  LOGICAL             ::  new_year, first_year, second_year
  LOGICAL             ::  read_fpc
  LOGICAL             ::  input_scaling
  LOGICAL             ::  lcc
  REAL(dp)            ::  step
  INTEGER             ::  read_status   ! error status reading the namelist
  INTEGER             ::  nml_unit      ! logical unit for the namelist


  ! --- constants ----

  integer,parameter   :: sec_per_day         = 86400
  integer,parameter   :: INTERVAL_DAY=1,INTERVAL_MONTH=2,INTERVAL_YEAR=3



  INCLUDE 'cbalone_ctl.inc'

  !! === GO =====================================================================

  !! --- initialization of MPI
  call p_start

  !! --- read namelist cbalone_ctl
  driver_data_path = ''
  experiment = ''
  specificName = ''
  climate_yearstart = -9999
  climate_yearend = -9999
  run_year_first = -9999
  run_year_last = -9999
  out_interval = 'YEAR'
  input_scaling = .FALSE.
  co2_file = 'co2.nc'
  climatology_diff_file = 'climdiff.nc'
  ref_year_past = -9999
  ref_year_recent = -9999
  use_dynveg = .FALSE.
  lcc = .FALSE.
  read_cpools=.FALSE.      !! false: Carbon pools are initialized with zero
                           !! true: Carbon pools initialized from an external file
  cpool_file_name='cpools.nc'
  read_fpc = .FALSE.
  fpc_file_name = 'fpc.nc'

  nml_unit = find_next_free_unit (51,100)
  CALL open_nml ('namelist.cbalone', unit=nml_unit)
  CALL position_nml ('CBALONE_CTL', status=read_status)
  SELECT CASE (read_status)
  CASE (POSITIONED)
     READ (nml_unit, cbalone_ctl)
  END SELECT

  IF (TRIM(driver_data_path) == '') &
       STOP 'No entry for DRIVER_DATA_PATH found in namelist.'
  IF (TRIM(experiment) == '') &
       STOP 'No entry for keyword EXPERIMENT found in namelist.'
  IF (TRIM(specificName) == '') &
       STOP 'No entry for keyword SPECIFICNAME found in namelist.'
  IF (climate_yearstart > climate_yearend) &
       STOP 'Namelist paramter climate_yearstart > climate_yearend.'
  IF (run_year_first > run_year_last) &
       STOP 'Namelist paramter run_year_first > run_year_last.'

  out_interval=toLower(out_interval)
  select case(TRIM(out_interval))
     case('day')
        outputInterval=INTERVAL_DAY
     case('month')
        outputInterval=INTERVAL_MONTH
     case('year')
        outputInterval=INTERVAL_YEAR
     case default
        stop 'No valid value for keyword OUT_INTERVAL found in namelist.'
  end select

  write(*,*) " === Cbalone info =============="
  write(*,*) " Climate data: First year: ",climate_yearstart
  write(*,*) " Climate data: Last  year: ",climate_yearend
  write(*,*) " First year of run : ",run_year_first
  write(*,*) " Last  year of run : ",run_year_last
  write(*,*) " Output interval   : ",trim(out_interval)
  write(*,*) " ==============================="

  IF (read_cpools) THEN
     WRITE(*,*) "Init file for carbon pools: ",TRIM(cpool_file_name)
  END IF

  IF (use_dynveg) THEN
     WRITE(*,*) "Dynamic vegetation module is active."
     IF (read_fpc) THEN
        WRITE(*,*) "Reading initial vegetation cover fractions from "//fpc_file_name
     END IF
  END IF

  LC_change%do_landcover_change = lcc

  data_scaling%do_input_scaling = input_scaling
  if (data_scaling%do_input_scaling) then
     write(*,*) "This is a run with climate input scaling"
     write(*,*) "CO2 development taken from file: ",trim(CO2_file)
     write(*,*) "Climatological differences of driver data are taken from file: ",TRIM(climatology_diff_file)
     write(*,*) "For input data scaling reference years ",ref_year_past," and ",ref_year_recent," are used."
  else
     write(*,*) "No input scaling performed in this run"
  end if

  !! --- initialization of lctlib

  call init_lctlib(lctlib_file, lctlib)
  write(*,*) "LCT-lib initialized"

  !! --- initialization of grid

  call initGrid(iniFileVegetation,grid)
  write(*,*) "Grid initialized."

  ! --- initialization of vegetation specific values
  call initVegetation(iniFileVegetation,grid,vegetation)
  write(*,*) "Vegetation initialized."

  ! --- init memory for cbalance structure

  call initCbalance(grid,vegetation%ntiles,cbalance)
  write(*,*) "Cbalance initialized."

  ! --- initialization of Cpools
  if (read_cpools) then
     call initializePools(grid,cbalance,vegetation,trim(cpool_file_name))
  else
     call initializePools(grid,cbalance,vegetation)
  end if

  ! --- initialization of memory for lancover change emissions
  if(LC_change%do_landcover_change) then
     call initLandCoverChange(grid,LC_change)

     !! --- printout parameters 
     
     write(*,*) "=== Parameters Landcover Change ============"
     write(*,*) "   frac_wood_2_atmos=", frac_wood_2_atmos
     write(*,*) "  frac_green_2_atmos=", frac_green_2_atmos
     write(*,*) "frac_reserve_2_atmos=",frac_reserve_2_atmos
     write(*,*) "============================================"

  end if

  ! --- initialise memory for dynamic vegetation
  IF (use_dynveg) THEN
     allocate(rock_fract(grid%nland))
     rock_fract(:) = 0._dp
     CALL config_dynveg(dynveg_params, dynveg_options)
     CALL init_offline_dynveg(grid, vegetation, dynveg, &
          dynveg_clim, rock_fract, read_fpc, fpc_file_name, dynveg_params, dynveg_options)
  END IF

  ! --- initialization of memory for input data scaling
  if(data_scaling%do_input_scaling) then
     call initDataScaling(CO2_file,climatology_diff_file,grid,vegetation%ntiles,ref_year_past,ref_year_recent, &
          run_year_first,run_year_last,data_scaling)
  end if

  !! Allocate local memory

  allocate(frac_npp_2_woodPool(1:grid%nland,vegetation%ntiles))
  allocate(frac_npp_2_reservePool(1:grid%nland,vegetation%ntiles))
  allocate(tau_Cpool_litter_leaf(1:grid%nland,vegetation%ntiles))
  allocate(tau_Cpool_litter_wood(1:grid%nland,vegetation%ntiles))
  allocate(LAI_shed_constant(1:grid%nland,vegetation%ntiles))
  allocate(frac_C_fast2atmos(1:grid%nland,vegetation%ntiles))
  allocate(Max_C_content_woods(1:grid%nland,vegetation%ntiles))
  allocate(specificLeafArea_C(1:grid%nland,vegetation%ntiles))
  allocate(veg_fract_correction(1:grid%nland,vegetation%ntiles))
  allocate(coverFractBox(1:grid%nland,vegetation%ntiles))
  allocate(carbon_2_atmos(1:grid%nland))
  allocate(carbon_2_fastSoilPool(1:grid%nland))
  allocate(carbon_2_slowSoilPool(1:grid%nland))

  allocate(coverFract_new(1:grid%nland,vegetation%ntiles))
  allocate(coverFract_interpolated(1:grid%nland,vegetation%ntiles))

  do itile=1,vegetation%ntiles
    frac_npp_2_woodPool(1:grid%nland,itile)    = lctLib%frac_npp_2_woodPool(vegetation%coverType(1:grid%nland,itile))
    frac_npp_2_reservePool(1:grid%nland,itile) = lctLib%frac_npp_2_reservePool(vegetation%coverType(1:grid%nland,itile))
    tau_Cpool_litter_leaf(1:grid%nland,itile)  = lctLib%tau_Cpool_litter_leaf(vegetation%coverType(1:grid%nland,itile))
    tau_Cpool_litter_wood(1:grid%nland,itile)  = lctLib%tau_Cpool_litter_wood(vegetation%coverType(1:grid%nland,itile))
    LAI_shed_constant(1:grid%nland,itile)      = lctLib%LAI_shed_constant(vegetation%coverType(1:grid%nland,itile))
    frac_C_fast2atmos(1:grid%nland,itile)      = lctLib%frac_C_fast2atmos(vegetation%coverType(1:grid%nland,itile))
    Max_C_content_woods(1:grid%nland,itile)    = lctLib%Max_C_content_woods(vegetation%coverType(1:grid%nland,itile))
    specificLeafArea_C(1:grid%nland,itile)     = lctLib%specificLeafArea_C(vegetation%coverType(1:grid%nland,itile))
    veg_fract_correction(1:grid%nland,itile)   = 1.0_dp - exp(-lctLib%MaxLAI(vegetation%coverType(1:grid%nland,itile))/2.0_dp)
  end do

  !! initialize fields

  cbalance%box_NEP_wholeRun(:) = 0.0_dp

  !! prepare for averages
  
  call nullifyFields(cbalance,LC_change)

  ! --- read in land cover fractions in case of landcover change

  if(LC_change%do_landcover_change) then !! Note: Here the landcover map of the year before the first year is needed!
     digits = int(log10(real(run_year_first-1)+epsilon(1.0_dp)))
     if(digits>4) then
        write(*,*) "ERROR (1): Incomplete programming: "//&
             "Year numbers (="//trim(int2string(run_year_first-1))//&
             ") with more than 4 digits cannot be handled for getting landcover fraction files" 
        stop
     end if
     LCC_fileName = "cover_fract."//"0000"//".nc"
     LCC_fileName(16-digits:16) = adjustl(int2string(run_year_first-1))
     call read_landcover_fractions(LCC_fileName,grid,vegetation%ntiles,vegetation%coverFract) !! Overwrite cover fractions read in from jsbach init file
  end if

  ! --- start of repeat loop

  run_year=run_year_first - 1

  RepeatLoop: do

     ! --- loop over all years of the data

     do climate_year=climate_yearstart,climate_yearend

        run_year = run_year+1 !! update year counter
        new_year = .true.

        write(*,*) "run_year/climate_year=",run_year,"/",climate_year

        !! preparations

        noOfDaysInThisYear = Get_GregorianYearLength(climate_year)

        !! Determine name of output file

        if(outputInterval /= INTERVAL_DAY) then
           outfile=''
           outFile=adjustl(trim(outStemName))
           ihlp=len(adjustl(trim(outStemName)))
           write(outFile(ihlp+1:ihlp+7),'(".000000")')
           digits = int(log10(real(run_year)+epsilon(1.0_dp)))
           write(sHlp,'(I6)') run_year
           outFile(ihlp+7-digits:ihlp+7) =  adjustl(trim(sHlp))
           outFile = trim(outFile)//"."//trim(int2string(climate_year))//".nc"
        end if

        !! loop over all months

        number_of_days_in_year = 0
        do month=1,12

           write(*,*) "run_year,month=",run_year,month

           !! Construct name of output file
           if(outputInterval == INTERVAL_DAY) then
              outfile=''
              outFile=adjustl(trim(outStemName))
              ihlp=len(adjustl(trim(outStemName)))
              write(outFile(ihlp+1:ihlp+7),'(".000000")')
              digits = int(log10(real(run_year)+epsilon(1.0_dp)))
              write(sHlp,'(I6)') run_year
              outFile(ihlp+7-digits:ihlp+7) =  adjustl(trim(sHlp))
              if(month<10) then
                 outFile = trim(outFile)//"0"//trim(int2string(month))//"."//trim(int2string(climate_year))//".nc"
              else
                 outFile = trim(outFile)//trim(int2string(month))//"."//trim(int2string(climate_year))//".nc"
              end if
           end if

           !! Read in new landcover fractions at begin of each January
           
           if(LC_change%do_landcover_change .and. month == 1) then
              digits = int(log10(real(run_year)+epsilon(1.0_dp)))
              if(digits>4) then
                 write(*,*) "ERROR (2): Incomplete programming: "//&
                      "Year numbers (="//trim(int2string(run_year))//&
                      ")  with more than 4 digits cannot be handled for getting landcover fraction files" 
                 stop
              end if
              LCC_fileName = "cover_fract."//"0000"//".nc"
              LCC_fileName(16-digits:16) = adjustl(int2string(run_year))
              call read_landcover_fractions(LCC_fileName,grid,vegetation%ntiles,coverFract_new)
           end if

           ! --- get data for a full month

           if(month<10) then
              write(sHlp(1:2),'("0",I1)') month
           else 
              write(sHlp(1:2),'(I2)') month
           end if
           if(climate_year >999) then
              inFile = trim(experiment)//trim(int2string(climate_year))//sHlp(1:2)//trim(specificName)//".nc"
           else 
              inFile = trim(experiment)//"0"//trim(int2string(climate_year))//sHlp(1:2)//trim(specificName)//".nc"
           end if
           call getDailyData(trim(driver_data_path)//trim(inFile),grid, use_dynveg, cbalance, dynveg_clim, &
                            vegetation%ntiles, nday)

           ! --- scale the data according to CO2 concentration

           if(data_scaling%do_input_scaling) call scaleDailyData(nday,data_scaling,run_year,month,cbalance)

           !! --- start day loop

           do day=1,nday

              ! --- do landcover change
              
              if(LC_change%do_landcover_change) then
                 dayInYear = getYearDay(climate_year,month,day)
                 coverFract_interpolated(:,:) =   vegetation%coverFract(:,:)                                                  &
                      + (coverFract_new(:,:) - vegetation%coverFract(:,:))/real(1 + noOfDaysInThisYear- dayInYear,dp)

                 call relocate_carbon(vegetation%coverFract(:,:),coverFract_interpolated(:,:),                           &
                      veg_fract_correction(:,:), epsilon(1._dp),                                                         &
                      cbalance%Cpool_green(:,:), cbalance%Cpool_woods(:,:), cbalance%Cpool_reserve(:,:),                 &
                      cbalance%Cpool_litter_leaf(:,:), cbalance%Cpool_litter_wood(:,:),                                  &
                      cbalance%Cpool_fast(:,:), cbalance%Cpool_slow(:,:),                                                &
                      frac_wood_2_atmos = frac_wood_2_atmos,                                                             &
                      frac_green_2_atmos = frac_green_2_atmos,                                                           &
                      frac_reserve_2_atmos = frac_reserve_2_atmos,                                                       &
                      carbon_2_atmos = carbon_2_atmos(:),                                                                &
                      carbon_2_fastSoilPool = carbon_2_fastSoilPool(:),                                                  &
                      carbon_2_slowSoilPool = carbon_2_slowSoilPool(:))

                 vegetation%coverFract(:,:) = coverFract_interpolated(:,:)

                 LC_change%LCC_sum_box_C2atmos(:)=     &
                      LC_change%LCC_sum_box_C2atmos(:)    + vegetation%vegRatioMax(:)*carbon_2_atmos(:)
                 LC_change%LCC_sum_box_C2fastSoilPool(:) =     &
                      LC_change%LCC_sum_box_C2fastSoilPool(:)     + vegetation%vegRatioMax(:)*carbon_2_fastSoilPool(:)
                 LC_change%LCC_sum_box_C2slowSoilPool(:) =     &
                      LC_change%LCC_sum_box_C2slowSoilPool(:)     + vegetation%vegRatioMax(:)*carbon_2_slowSoilPool(:)
                 LC_change%LCC_emissions_wholeRun(:) = &
                      LC_change%LCC_emissions_wholeRun(:) + vegetation%vegRatioMax(:)*carbon_2_atmos(:)*molarMassCO2_kg
              end if

              !! vegetation cover fractions may have changed, so redetermine fractions in box for later use here

              coverFractBox(:,:) = vegetation%coverFract(:,:) * veg_fract_correction(:,:)      &
                                   * spread(vegetation%vegRatioMax(:),DIM=2,NCOPIES=vegetation%ntiles)

              ! update Cpools
              if (day > 1) then
                cbalance%LAI_previousDayMean(1:grid%nland,1:vegetation%ntiles) = &
                  cbalance%LAI_yDaymean(1:grid%nland,1:vegetation%ntiles,day-1)
              end if

              call update_Cpools(cbalance%LAI_yDaymean(1:grid%nland,1:vegetation%ntiles,day),         &
                   cbalance%LAI_previousDayMean(1:grid%nland,1:vegetation%ntiles),      &
                   cbalance%NPP_yDayMean(1:grid%nland,1:vegetation%ntiles,day),         &
                   cbalance%topSoilTemp_yDayMean(1:grid%nland,1:vegetation%ntiles,day), &
                   cbalance%alpha_yDayMean(1:grid%nland,1:vegetation%ntiles,day),       &
                   frac_npp_2_woodPool(1:grid%nland,1:vegetation%ntiles),               &
                   frac_npp_2_reservePool(1:grid%nland,1:vegetation%ntiles),            &
                   tau_Cpool_litter_leaf(1:grid%nland,1:vegetation%ntiles),             &
                   tau_Cpool_litter_wood(1:grid%nland,1:vegetation%ntiles),             &
                   LAI_shed_constant(1:grid%nland,1:vegetation%ntiles),                 &
                   frac_C_fast2atmos(1:grid%nland,1:vegetation%ntiles),                 &
                   Max_C_content_woods(1:grid%nland,1:vegetation%ntiles),               &
                   specificLeafArea_C(1:grid%nland,1:vegetation%ntiles),                &
                   cbalance%Cpool_green(1:grid%nland,1:vegetation%ntiles),              &
                   cbalance%Cpool_woods(1:grid%nland,1:vegetation%ntiles),              &
                   cbalance%Cpool_reserve(1:grid%nland,1:vegetation%ntiles),            &
                   cbalance%Cpool_litter_leaf(1:grid%nland,1:vegetation%ntiles),        &
                   cbalance%Cpool_litter_wood(1:grid%nland,1:vegetation%ntiles),        &
                   cbalance%Cpool_fast(1:grid%nland,1:vegetation%ntiles),               &
                   cbalance%Cpool_slow(1:grid%nland,1:vegetation%ntiles),               &
                   cbalance%soil_respiration(1:grid%nland,1:vegetation%ntiles),         &
                   cbalance%NPP_flux_correction(1:grid%nland,1:vegetation%ntiles),      &
                   cbalance%litter_flux(1:grid%nland,1:vegetation%ntiles)                &
                   )
              if (day == nday) then
                cbalance%LAI_previousDayMean(1:grid%nland,1:vegetation%ntiles) = &
                  cbalance%LAI_yDaymean(1:grid%nland,1:vegetation%ntiles,day)                 
              endif

             cbalance%avg_Cpool_green(:,:)   = cbalance%avg_Cpool_green(:,:)  + coverFractBox(:,:)*cbalance%Cpool_green(:,:)
             cbalance%avg_Cpool_woods(:,:)   = cbalance%avg_Cpool_woods(:,:)  + coverFractBox(:,:)*cbalance%Cpool_woods(:,:)
             cbalance%avg_Cpool_reserve(:,:) = cbalance%avg_Cpool_reserve(:,:)+ coverFractBox(:,:)*cbalance%Cpool_reserve(:,:)
             cbalance%avg_Cpool_litter_leaf(:,:) = cbalance%avg_Cpool_litter_leaf(:,:)+ &
                                                   coverFractBox(:,:)*cbalance%Cpool_litter_leaf(:,:)
             cbalance%avg_Cpool_litter_wood(:,:) = cbalance%avg_Cpool_litter_wood(:,:)+ &
                                                   coverFractBox(:,:)*cbalance%Cpool_litter_wood(:,:)
             cbalance%avg_Cpool_fast(:,:)    = cbalance%avg_Cpool_fast(:,:)   + coverFractBox(:,:)*cbalance%Cpool_fast(:,:)
             cbalance%avg_Cpool_slow(:,:)    = cbalance%avg_Cpool_slow(:,:)   + coverFractBox(:,:)*cbalance%Cpool_slow(:,:)

#ifndef __PGI
             cbalance%box_Cpools_total(:)    = sum(coverFractBox(:,:)                                                  &
                                                     * (   cbalance%Cpool_green(:,:)   + cbalance%Cpool_woods(:,:)   &
                                                   + cbalance%Cpool_litter_leaf(:,:) + cbalance%Cpool_litter_wood(:,:) &
                                                         + cbalance%Cpool_reserve(:,:) + cbalance%Cpool_fast(:,:)    &
                                                         + cbalance%Cpool_slow(:,:) )                                &
                                                   ,DIM=2)
#else
             !
             ! Different notation to allow compilation with the PGI compiler (pgf95 6.1-1)
             !
             cbalance%box_Cpools_total(:) = 0._dp
             do itile=1,vegetation%ntiles
                cbalance%box_Cpools_total(:) =  cbalance%box_Cpools_total(:) + coverFractBox(:,itile) &
                                               * ( cbalance%Cpool_green(:,itile) + cbalance%Cpool_woods(:,itile) &
                                                 + cbalance%Cpool_litter_leaf(:,itile) + cbalance%Cpool_litter_wood(:,itile) &
                                                 + cbalance%Cpool_reserve(:,itile) + cbalance%Cpool_fast(:,itile) &
                                                 + cbalance%Cpool_slow(:,itile) )
             end do
#endif
             cbalance%box_NEP_wholeRun(:) =                                                                                     &
                     cbalance%box_NEP_wholeRun(:)                                                                               &
                   + sum( coverFractBox(:,:)                                                                                    &
                         * (cbalance%NPP_yDayMean(:,:,day) + cbalance%NPP_flux_correction(:,:) + cbalance%soil_respiration(:,:)) &
                         ,DIM=2) * sec_per_day

             cbalance%avg_soil_respiration(:,:)    = cbalance%avg_soil_respiration(:,:)                             &
                                                        + coverFractBox(:,:) * cbalance%soil_respiration(:,:)
             cbalance%avg_NPP_yDayMean(:,:)        = cbalance%avg_NPP_yDayMean(:,:)                                 &
                                                        + coverFractBox(:,:) * cbalance%NPP_yDayMean(:,:,day)
             cbalance%avg_NPP_flux_correction(:,:) = cbalance%avg_NPP_flux_correction(:,:)                          &
                                                        + coverFractBox(:,:) * cbalance%NPP_flux_correction(:,:)
             cbalance%avg_litter_flux(:,:)          = cbalance%avg_litter_flux(:,:)                                   &
                                                        + coverFractBox(:,:) * cbalance%litter_flux(:,:)
             cbalance%avg_box_NEP(:) = cbalance%avg_box_NEP(:)                                                                   &
                   + sum(coverFractBox(:,:)                                                                                      &
                         * (cbalance%NPP_yDayMean(:,:,day) + cbalance%NPP_flux_correction(:,:) + cbalance%soil_respiration(:,:)) &
                         ,DIM=2) * molarMassCO2_kg
            !
            ! --- update dynamic vegetation
            !
             IF (use_dynveg) THEN
                CALL update_offline_dynveg(day, run_year, run_year_first, new_year, &
                     grid, vegetation, dynveg, dynveg_clim, cbalance, dynveg_params, dynveg_options, &
                     specificLeafArea_C(:,:), rock_fract(:), veg_fract_correction(:,:))
             END IF

             number_of_days_in_year = number_of_days_in_year + 1

             if(outputInterval == INTERVAL_DAY) then
                step = real(run_year + real(month-1)/12._dp + real(day-1)/real(nday)/12._dp)
                call writeSingleTimeStep(trim(outFile), grid, cbalance, use_dynveg, dynveg_options%dynveg_feedback, &
                                         dynveg, vegetation%ntiles, LC_change, step)
                call nullifyFields(cbalance,LC_change)
             end if

             new_year = .false.
          end do ! end of day loop

          if(outputInterval == INTERVAL_DAY) then !! close the output file
             call writeSingleTimeStep("", grid, cbalance, use_dynveg, dynveg_options%dynveg_feedback,&
                                      dynveg, 0, LC_change, 0.0_dp, .true.) !! close the output file
          end if

          if(outputInterval == INTERVAL_MONTH) then
             call averageFields(cbalance,nday)
             step = real(run_year + real(2*month-1)/24._dp)
             call writeSingleTimeStep(trim(outFile), grid, cbalance, use_dynveg, dynveg_options%dynveg_feedback, &
                                      dynveg, vegetation%ntiles, LC_change, step)
             call nullifyFields(cbalance,LC_change)
          end if

       end do ! end month loop

       if(outputInterval == INTERVAL_MONTH) then !! close the output file
          call writeSingleTimeStep("",grid,cbalance,use_dynveg, dynveg_options%dynveg_feedback, dynveg, 0, &
                                   LC_change,0.0_dp,.true.) !! close the output file
       end if

       if(outputInterval == INTERVAL_YEAR) then
          call averageFields(cbalance,number_of_days_in_year)
          step = real(run_year)
          call writeSingleTimeStep(trim(outFile), grid, cbalance, use_dynveg, dynveg_options%dynveg_feedback, &
                                   dynveg, vegetation%ntiles, LC_change, step)
          call nullifyFields(cbalance,LC_change)
          call writeSingleTimeStep("", grid, cbalance, use_dynveg, dynveg_options%dynveg_feedback, dynveg, 0, &
                                   LC_change, 0.0_dp, .true.) !! close the output file
       end if

       !! Check for end of run

       if(run_year == run_year_last) then
          exit RepeatLoop
       end if

       !! Force writing to stdout

       call flush(6)

     end do               ! end loop over all climate data years
     
  end do RepeatLoop                            ! end repeat loop

  write(*,*) "Normal End!"
  

contains

  subroutine averageFields(cbalance,noOfTimeSteps)
    USE mo_cbalone_memory, ONLY: cbal_offline_type
    TYPE(cbal_offline_type),intent(inout) :: cbalance
    integer,intent(in) :: noOfTimeSteps

    cbalance%avg_Cpool_green(:,:) = cbalance%avg_Cpool_green(:,:)/real(noOfTimeSteps)
    cbalance%avg_Cpool_woods(:,:) = cbalance%avg_Cpool_woods(:,:)/real(noOfTimeSteps)
    cbalance%avg_Cpool_reserve(:,:) = cbalance%avg_Cpool_reserve(:,:)/real(noOfTimeSteps)
    cbalance%avg_Cpool_litter_leaf(:,:) = cbalance%avg_Cpool_litter_leaf(:,:)/real(noOfTimeSteps)
    cbalance%avg_Cpool_litter_wood(:,:) = cbalance%avg_Cpool_litter_wood(:,:)/real(noOfTimeSteps)
    cbalance%avg_Cpool_fast(:,:) = cbalance%avg_Cpool_fast(:,:)/real(noOfTimeSteps)
    cbalance%avg_Cpool_slow(:,:) = cbalance%avg_Cpool_slow(:,:)/real(noOfTimeSteps)

    cbalance%avg_soil_respiration(:,:)    = cbalance%avg_soil_respiration(:,:)/real(noOfTimeSteps)
    cbalance%avg_NPP_yDayMean(:,:)        = cbalance%avg_NPP_yDayMean(:,:)/real(noOfTimeSteps)
    cbalance%avg_NPP_flux_correction(:,:) = cbalance%avg_NPP_flux_correction(:,:)/real(noOfTimeSteps)
    cbalance%avg_litter_flux(:,:)         = cbalance%avg_litter_flux(:,:)/real(noOfTimeSteps)
    cbalance%avg_box_NEP(:)               = cbalance%avg_box_NEP(:)/real(noOfTimeSteps)

  end subroutine averageFields

  subroutine nullifyFields(cbalance,LC_change)
    USE mo_cbalone_memory, ONLY: cbal_offline_type
    TYPE(cbal_offline_type),intent(inout)     :: cbalance
    TYPE(landcover_change_type),intent(inout) :: LC_change

    cbalance%avg_Cpool_green(:,:)   = 0._dp
    cbalance%avg_Cpool_woods(:,:)   = 0._dp
    cbalance%avg_Cpool_reserve(:,:) = 0._dp
    cbalance%avg_Cpool_litter_leaf(:,:)   = 0._dp
    cbalance%avg_Cpool_litter_wood(:,:)   = 0._dp
    cbalance%avg_Cpool_fast(:,:)    = 0._dp
    cbalance%avg_Cpool_slow(:,:)    = 0._dp

    cbalance%avg_soil_respiration(:,:)   = 0._dp
    cbalance%avg_NPP_yDayMean(:,:)       = 0._dp
    cbalance%avg_NPP_flux_correction(:,:)= 0._dp
    cbalance%avg_litter_flux(:,:)        = 0._dp
    cbalance%avg_box_NEP(:)              = 0._dp

    if(LC_change%do_landcover_change) then
       LC_change%LCC_sum_box_C2atmos(:) = 0.0_dp
       LC_change%LCC_sum_box_C2fastSoilPool(:)  = 0.0_dp
       LC_change%LCC_sum_box_C2slowSoilPool(:)  = 0.0_dp
    end if

  end subroutine nullifyFields

  !! --- getYearDay() -----------------------------
  !!
  !! Returns the number of a day in a year, when counts starts with 1 at first of january
  !!
  integer function getYearDay(year,month,day)
    integer,intent(in) :: year
    integer,intent(in) :: month
    integer,intent(in) :: day

    logical :: leapYear

    integer,parameter :: lastDayOfPrevMonth_normal(1:12)   =(/0,31,59,90,120,151,181,212,243,273,304,334/)
    integer,parameter :: lastDayOfPrevMonth_leapYear(1:12) =(/0,31,60,91,121,152,182,213,244,274,305,335/)

    leapYear = (MOD(year,4)==0 .AND. MOD(year,100)/=0) .OR. MOD(year,400)==0 
    if(leapYear) then
       getYearDay = lastDayOfPrevMonth_leapYear(month) + day
    else
       getYearDay = lastDayOfPrevMonth_normal(month) + day
    end if

  end function getYearDay

  !! --- get_GregorianYearLength() -----------------------------
  !!
  !! Returns the number days in a year according to the (proleptic) gregorian calender 
  !!
  integer function get_GregorianYearLength(year)
    integer,intent(in) :: year

    if ( (MOD(year,4)==0 .AND. MOD(year,100)/=0) .OR. MOD(year,400)==0 ) then
       get_gregorianYearLength = 366
    else
       get_gregorianYearLength = 365
    end if

  end function get_gregorianYearLength

end program cbalone_driver
