! ==================================================================================================================================
!
!                                               DRIVER OF THE STANDALONE-JSBACH
!
! ==================================================================================================================================

PROGRAM jsbalone_driver

  USE mo_kind,             ONLY : dp
  USE mo_jsbach,           ONLY : options_type, jsbach_FirstCall
  USE mo_mpi,              ONLY : p_pe, p_io, p_parallel,p_start, p_stop, p_bcast,p_parallel_io  
  USE mo_time_control,     ONLY : init_manager, init_times, lfirst_cycle, init_events, &
                                  construct_events, get_time_step, putdata, putrerun, trigfiles, &
                                  p_bcast_event, time_set, time_reset, l_putrerun, lstop, &
                                  l_trigfiles, current_date
  USE mo_timer,            ONLY : init_timer, cleanup_timer
  USE mo_doctor,           ONLY : nout
  USE mo_jsbalone_forcing, ONLY : init_forcing,update_forcing,finish_forcing
  USE mo_jsbalone_forcing, ONLY : forcing_type                  ! This structure contains all the forcing fields
  USE mo_zenith,           ONLY : compute_declination, get_zenith_angle
  USE mo_grib,             ONLY : open_output_streams, close_output_streams, out_streams
  USE mo_io,               ONLY : write_streams
  USE mo_jsbach_grid,      ONLY : grid_type, domain_type
  USE mo_exception,        ONLY : finish, message, message_text, int2string
  USE mo_jsbach_interface, ONLY : jsbach_inter_1d, jsbach_init    ! The interface to the land surface model
  USE mo_jsbach_constants, ONLY : Gravity
  USE mo_kind,             ONLY : dp
  USE mo_jsbalone,         ONLY : init_driving, update_driving, drive_type
  USE mo_machine,          ONLY : machine_setup
!!$  USE mo_jsbach_bvocemis, ONLY: bvoc_type, init_bvoc
  USE mo_param_switches,   ONLY : lsurf

  implicit none

  ! Local structures to get grid and domain description back from the interface (jsbach_init)
  TYPE(grid_type)   :: grid
  TYPE(domain_type) :: domain

  ! Local structure to get JSBACH options from interface (jsbach_init)
  TYPE(options_type) :: options

  ! Local structure to hold forcing fields
  TYPE(forcing_type) :: forcing

  ! Local structure for driving field (derived vars from forcing fields)
  TYPE(drive_type)      :: driving

  ! Local structure for bvoc emisions from mo_jsbach_bvocemis (adapted from megan)
!!$  TYPE(bvoc_type)    :: bvoc

  INTEGER, EXTERNAL :: util_cputime !  External functions
  REAL(dp), EXTERNAL:: util_walltime
  REAL(dp):: zutime, zstime, zrtime, zwtime

  INTEGER           :: istep,status,jsec_add
  INTEGER           :: ntile
  INTEGER           :: fileformat     ! output file format (grib/netcdf)
  REAL(dp) :: time_add
  REAL(dp),DIMENSION(:),ALLOCATABLE :: cos_zenith, evapo_act, evapo_pot
  REAL(dp) :: declination

  !! only for testing BEGIN
  real(dp),ALLOCATABLE :: hlpField1(:)
  !! only for testing END


  !=================================================================================================
  ! Start MPI
  CALL p_start
  !
  !=================================================================================================
  ! Initialize wallclock timer
!$OMP PARALLEL
!$OMP MASTER
  zwtime = util_walltime()
!$OMP END MASTER
!$OMP END PARALLEL
  !=================================================================================================
  ! Initialize JSBACH (read configuration, set dates, compute decomposition, 
  !                    read restart or ini files for grid, soil, vegetation parameters etc.)
  CALL machine_setup

  CALL jsbach_init(grid, domain, options, ntile)

  CALL init_forcing(grid, domain, options, forcing)

  CALL init_driving(grid, domain, options, forcing, driving, ntile)

!!$  CALL init_bvoc(grid, domain, options, bvoc , ntile) 

  lsurf = .true.

  !
  !=================================================================================================
  ! Allocation of fields needed to exchange data with interface

    ALLOCATE(cos_zenith(domain%nland), STAT=status)
    if(status /= 0) call finish("jsbach_driver", "ERROR: Field cos_zenith() could not be allocated.")
    ALLOCATE(evapo_act(domain%nland),STAT=status)
    if(status /= 0) call finish("jsbach_driver", "ERROR: Field forcingFields%air_temp() could not be allocated.")
    ALLOCATE(evapo_pot(domain%nland),STAT=status)
    if(status /= 0) call finish("jsbach_driver", "ERROR: Field forcingFields%air_temp() could not be allocated.")

    evapo_pot(:)= 0.5_dp !! ONLY FOR TEST-PURPOSES!
    evapo_act(:)= 0.5_dp !! ONLY FOR TEST-PURPOSES!

  ! Main loop

  main_loop:DO

     ! evaluate events at begin of time step

     CALL time_set

     IF (l_trigfiles) THEN
        CALL close_output_streams
        CALL open_output_streams
     ENDIF

     ! compute cosine of zenith angle
     jsec_add = 0
     CALL compute_declination(jsec_add, declination)
     time_add = 0._dp
     cos_zenith = get_zenith_angle(domain%nland, domain, time_add, declination)

     ! generate forcing data for particular time step 
     CALL update_forcing(grid, domain,cos_zenith, declination, evapo_act, evapo_pot, forcing)

     CALL out_streams
     IF (l_putrerun) THEN
        CALL write_streams
     ENDIF

!TESTING
!forcing%air_temp(1:domain%nland) = -20.0_dp
!forcing%rad_PAR_down(1:domain%nland) = 0.1_dp
!forcing%rad_lw_down(1:domain%nland) = 0.1_dp
!forcing%rad_sw_down(1:domain%nland) = 50.0_dp

     CALL update_driving( domain%nland, &
                          domain%nland, &
                          driving, &
                          domain%elev(1:domain%nland) * Gravity, &
                          forcing%wind_speed(1:domain%nland), &
                          forcing%wind_speed(1:domain%nland), & ! should be 10 meter wind
                          !forcing%wind_u(1:domain%nland), &
                          !forcing%wind_v(1:domain%nland), &
                          forcing%air_temp(1:domain%nland)+273.15_dp,   &
                          forcing%spec_humidity(1:domain%nland), &
                          forcing%precip_rain(1:domain%nland), &
                          forcing%precip_snow(1:domain%nland), &
                          forcing%rad_lw_down(1:domain%nland), &   
                          forcing%rad_PAR_down(1:domain%nland), & ! should be sw_vis_net
                          forcing%frac_PAR_diffuse(1:domain%nland), & !should be sw_vis_frac_diffuse
                          forcing%rad_NIR_down(1:domain%nland), & ! should be sw_nir_net
                          forcing%rad_NIR_down(1:domain%nland), & ! should be sw_nir_frac_diffuse (not used in jsbach yet)
                          forcing%air_pressure(1:domain%nland), &
                          cos_zenith(1:domain%nland), &
                          forcing%CO2_concentr(1:domain%nland), &
                          options%HeightHumidity, &               ! Attention: Is used for 'LayerHeight'-> Sechiba uses 2m
                          jsbach_FirstCall, & 
                          driving%sensible_heat(1:domain%nland), &
                          driving%evap_act(1:domain%nland), &
                          driving%soil_temp(1:domain%nland), &
                          domain%elev(1:domain%nland) * Gravity &
                          )
     ! call JSBACH-interface to land surface

     call jsbach_inter_1d(domain%nland, &
                          domain%nland, &
                          domain%elev(1:domain%nland) * Gravity, &
                          forcing%wind_speed(1:domain%nland), &
                          forcing%wind_speed(1:domain%nland), & ! should be 10 meter wind
                          !forcing%wind_u(1:domain%nland), &
                          !forcing%wind_v(1:domain%nland), &
                          forcing%air_temp(1:domain%nland)+273.15_dp,   &
                          forcing%spec_humidity(1:domain%nland), &
                          forcing%precip_rain(1:domain%nland), &
                          forcing%precip_snow(1:domain%nland), &
                          forcing%rad_lw_down(1:domain%nland), &   
                          forcing%rad_PAR_down(1:domain%nland), & ! should be sw_vis_net
                          forcing%frac_PAR_diffuse(1:domain%nland), & !should be sw_vis_frac_diffuse
                          forcing%rad_NIR_down(1:domain%nland), & ! should be sw_nir_net
                          forcing%rad_NIR_down(1:domain%nland), & ! should be sw_nir_frac_diffuse (not used in jsbach yet)
                          forcing%air_pressure(1:domain%nland), &
                          cos_zenith(1:domain%nland), &
                          declination,                &
                          forcing%CO2_concentr(1:domain%nland), &
                          evap_act=driving%evap_act(1:domain%nland), & 
                          evap_pot=driving%evap_pot(1:domain%nland), & !output from jsbach ?!?
                          sensible=driving%sensible_heat(1:domain%nland), &
                          temp_soil_new=driving%soil_temp(1:domain%nland), &
                          etacoef = driving%etacoef(1:domain%nland)  , &
                          etbcoef = driving%etbcoef(1:domain%nland)  , &
                          eqacoef = driving%eqacoef(1:domain%nland)  , &
                          eqbcoef = driving%eqbcoef(1:domain%nland),  &
                          echam_zchl = driving%zchl(1:domain%nland) ,&
                          cdrag = driving%cdrag(1:domain%nland), &
                          latent = driving%latent(1:domain%nland) &                          
                          )


     CALL time_reset 

     ! check for last time step (then lstop=.true.) and eventually exit the main loop

     IF (lstop) THEN 
!        IF (p_parallel) THEN
           istep = get_time_step() ! get the number of the time step
           CALL message('jsbach_driver','Step '//int2string(istep)//'completed.')
!           WRITE (nout,'(a,i7,a,i4,a)') &
!                'Step ', istep, 'on PE ', p_pe, ' completed.'
!        ELSE
!           WRITE (nout,'(a,i7,a)') 'Step ', istep, ' completed.'
!        END IF

        EXIT main_loop   
     END IF
      
  END DO main_loop

  IF (lstop) THEN
     CALL close_output_streams
  END IF

  ! Cleanup timer
  ! CALL cleanup_timer

  ! Deallocate memory of forcing fields

  call finish_forcing(forcing)

  DEALLOCATE(cos_zenith, evapo_pot, evapo_act)

!$OMP PARALLEL
!$OMP MASTER
  status = util_cputime(zutime, zstime)
!$OMP END MASTER
!$OMP END PARALLEL
  IF (status == -1) THEN
     CALL message('stepon','Cannot determine used CPU time')
  ELSE
!$OMP PARALLEL
!$OMP MASTER
     zwtime = util_walltime()
!$OMP END MASTER
!$OMP END PARALLEL
     zrtime = (zutime+zstime)/zwtime
     CALL message ('', '')
     WRITE (message_text,'(a,f10.2,a)') ' Wallclock        : ', zwtime, ' s'
     CALL message('',message_text)
     WRITE (message_text,'(a,f10.2,a)') ' CPU-time (user)  : ', zutime, ' s'
     CALL message('',message_text)
     WRITE (message_text,'(a,f10.2,a)') ' CPU-time (system): ', zstime, ' s'
     CALL message('',message_text)
     WRITE (message_text,'(a,f10.2,a)') ' Ratio            : ', 100*zrtime, ' %'
     CALL message('',message_text)
     CALL message ('', '')
  END IF

  ! Stop MPI
  CALL p_stop

END PROGRAM jsbalone_driver
