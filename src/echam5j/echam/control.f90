#if defined(__uxp__) || defined(__SX__) || defined(ES)
#define FAST_AND_DIRTY 1
#endif

SUBROUTINE control

  ! Description:
  !
  ! Control routine for the model.
  !
  ! Method:
  !
  ! This subroutine controls the running of the model.
  ! 
  ! *control* is called from the main program (*master*).
  !
  ! Externals:
  !   *initialize*    called to initialize modules.
  !   *stepon*        controls time integration
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! R.G and M.J, ECMWF, December 1982, changed
  ! U. Schlese, MPI, August 1989, new structure
  ! U. Schlese, DKRZ, September 1994, interface *drive* added
  ! U. Schlese, DKRZ, January 1995, reading of optional files added
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schlese, DKRZ and M. Esch, MPI, July 1999, modifications for ECHAM5
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! U. Schlese, DKRZ, November 1999, "drive" removed, "stepon" called directly
  ! U. Schlese, DKRZ, December 1999, modifications for coupling
  ! M.A. Giorgetta, MPI, May 2000, ozone initialization removed, see setrad
  ! S. Legutke, MPI M&D, July 00, modifications for coupling interface
  ! I. Kirchner, MPI, December 2000, date/time control/nudging
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! I. Kirchner, MPI, Aug 2002, tendency diagnostics revision
  ! A. Rhodin, DWD, June 2002, call subroutines cleanup_...
  ! M. Esch, MPI, September 2002, modifications for mixed layer ocean
  ! L. Kornblueh, MPI, October 2003, added setup for AMIP2 global diagnostics
  ! U. Schulzweida, MPI, March 2007, added daily SST and SIC support
  ! 
  ! for more details see file AUTHORS
  ! 

  USE mo_sst,             ONLY: readsst, readice, readflux, &
                                readdailysst, readdailyice
  USE mo_control,         ONLY: lcouple, lnudge, ltdiag, lhd, &
                                lnmi, lmlo, numfl1, numfl2,   &
                                ldiagamip, lso4, ldailysst
  USE mo_machine,         ONLY: machine_setup
  USE mo_legendre,        ONLY: inileg, cleanup_legendre
  USE mo_tracer,          ONLY: reademi, cleanstatr
  USE mo_nmi,             ONLY: NMI_Init, NMI_Close
  USE mo_time_control,    ONLY: lstart, lfirst_cycle, lresume,         &
                                construct_events, lstop, l_diagrad
  USE mo_clim,            ONLY: readtslclim, readvltclim, readvgratclim
  USE mo_nudging_init,    ONLY: NudgingInit, NDG_CLEAN_MEM, NDG_CLOSE
  USE mo_diag_tendency,   ONLY: DIAG_Init, IDIAG_INIT, IDIAG_FREE
  USE mo_diag_amip2,      ONLY: init_amip2_diag
  USE mo_couple,          ONLY: couple_init
  USE mo_column,          ONLY: resetcolumn
  USE mo_grib,            ONLY: init_grib, cleanup_grib
  USE mo_memory_streams,  ONLY: init_memory, free_memory
  USE mo_geoloc,          ONLY: init_geoloc, cleanup_geoloc
  USE mo_advection,       ONLY: iadvec, tpcore, semi_lagrangian,       &
                                spitfire
  USE mo_semi_lagrangian, ONLY: init_semi_lagrangian,                  &
                                cleanup_semi_lagrangian
  USE mo_spitfire,        ONLY: init_spitfire
  USE mo_tpcore,          ONLY: init_tpcore, cleanup_tpcore
  USE mo_timer,           ONLY: init_timer, cleanup_timer
  USE mo_aero_tanre,      ONLY: init_aero, cleanup_aero
  USE mo_radiation,       ONLY: iaero
  USE mo_decomposition,   ONLY: ldc=>local_decomposition,              &
                                cleanup_decomposition
  USE mo_exception,       ONLY: message_text, message, finish
  USE mo_hydrology,       ONLY: init_hydrology, cleanup_hydrology
  USE mo_scan_buffer,     ONLY: cleanup_scanbuffer
  USE mo_o3clim,          ONLY: cleanup_o3clim
  USE mo_clim,            ONLY: cleanup_clim
  USE mo_sst,             ONLY: cleanup_sst
  USE mo_tmp_buffer,      ONLY: cleanup_tmp_buffer
  USE mo_gaussgrid,       ONLY: cleanup_gaussgrid
#ifdef FFT991
  USE mo_fft991,          ONLY: cleanup_fft991
#else
  USE mo_fft992,          ONLY: cleanup_fft992
#endif
  USE m_alloc_mods,       ONLY: dealloc_mods
  USE mo_netcdf,          ONLY: cleanup_netcdf
  USE mo_io,              ONLY: cleanup_io
  USE mo_doctor,          ONLY: nerr, nout
  USE mo_field,           ONLY: field1, field2
  USE mo_so4,             ONLY: cleanup_so4
  USE mo_surface,         ONLY: init_surface
  USE mo_hd_highres_io,   ONLY: hd_highres_open, &
                                hd_highres_close
!---wiso-code
  USE mo_memory_wiso,     ONLY: cleanup_scanbuffer_wiso
  USE mo_wiso,            ONLY: lwiso
!---wiso-code-end

#ifdef FAST_AND_DIRTY
  USE mo_memory_f,        ONLY: resort_restart_read_memory_f
#endif

  IMPLICIT NONE

!  External subroutines 
  EXTERNAL :: stepon, initialize, inipost,                                &
              labrun, readfld, inhysi, ioinitial, scan2

!  Executable statements 

!-- Default advection scheme

  iadvec = tpcore

!--  Print machine specific values

  CALL machine_setup

  CALL construct_events

  CALL init_timer

!--  Initialize modules and parallel decomposition

  CALL initialize

!-- Initialize time independent surface parameters 

  CALL init_geoloc

  IF ( iaero == 2 .OR. iaero == 3 .OR. iaero == 4 ) THEN
    CALL init_aero
  END IF

  IF ( .NOT. ldc% lreg ) THEN
    IF ( iadvec .EQ. semi_lagrangian .OR.                              &
         iadvec .EQ. spitfire        .OR.                              &
         l_diagrad                   .OR.                              &
         lnudge                            ) THEN
      WRITE(message_text,*) "NPROMA does not work with:"
      CALL message('control', TRIM(message_text))
      WRITE(message_text,*) "  - semi_lagrangian advection"
      CALL message('control', TRIM(message_text))
      WRITE(message_text,*) "  - spitfire advection"
      CALL message('control', TRIM(message_text))
      WRITE(message_text,*) "  - l_diagrad"
      CALL message('control', TRIM(message_text))
      WRITE(message_text,*) "  - nudging"
      CALL message('control', TRIM(message_text))
      CALL finish('control','NPROMA not implemented')
    END IF
  END IF

!-- Initialize tendency diagnostics
 
  IF (ltdiag) CALL DIAG_Init(IDIAG_INIT)

!--  Initialize memory

  CALL init_memory

  !-- Preset values needed in the advection scheme

  SELECT CASE (iadvec)
  CASE (semi_lagrangian)
    CALL init_semi_lagrangian
  CASE (spitfire)
    CALL init_spitfire
  CASE (tpcore)
    CALL init_tpcore
  END SELECT

!
!--  Initialize HD Model
!
  IF (lhd) THEN
    write (0,*) 'HD model selected.'
    CALL init_hydrology
!hd-hi switched off HD high resolution output
!hd-hi    CALL hd_highres_open
  ENDIF

!--  Compute *Legendre polynomials and parameters needed
!    for the *Legendre transforms.

  CALL inileg

!---wiso-code

! read climate land surface temperatures before calling ioinitial 
! to enable temperature-depended initialisation of water isotope surface fields
! (initialisation is done in *ioinitial*) 

  CALL readtslclim

!---wiso-code-end

  IF (lresume) THEN
     CALL iorestart
  ELSE IF (lstart) THEN
     CALL ioinitial
  END IF

  CALL init_surface

  CALL inhysi

#ifdef FAST_AND_DIRTY  
  CALL resort_restart_read_memory_f
#endif
  IF (lstart) CALL scan2

  IF (lnmi .AND. lfirst_cycle) CALL NMI_Init(lnudge)


!-- Initialize postprocessing

  CALL inipost

  IF (ldiagamip) CALL init_amip2_diag

  CALL init_grib

!-- Read optional sst-file if not in coupled mode

  IF (.NOT. lcouple) THEN
    IF ( ldailysst ) THEN
      CALL readdailysst
    ELSE
      CALL readsst
    END IF
  END IF
  
! --  read optional sea ice if not in coupled mode

  IF (.NOT. lcouple) THEN
    IF ( ldailysst ) THEN
      CALL readdailyice
    ELSE
      CALL readice
    END IF
  END IF
#ifdef __cpl_fluxes4
  IF (lcouple ) CALL readice
#endif

! --  read optional flux correction if in mixed layer mode

  IF (lmlo) CALL readflux

!---wiso-code

! --  read optional isotope values of ocean surface waters if not in coupled mode

  IF (lwiso .AND. (.NOT. lcouple)) THEN
    CALL readwisosw_d
  END IF

! Comment out call of readtslclim at this point (routine has already been called above)

!-- Read  climate land surface temperatures

!  CALL readtslclim

!---wiso-code-end

!-- Read climate leaf-area index and climate vegetation ratio

  CALL readvltclim

  CALL readvgratclim

!-- Read optional fields

  CALL reademi

  CALL readfld

!
!---  control a coupled run
!        -initialise  data exchange with coupler
!        -check whether Echam control parameters are consistent  
!                                          with coupling system 

  IF (lcouple .AND. lfirst_cycle) THEN 
    CALL couple_init(nout,nerr)
  END IF

  CALL labrun

!-- Start time integration

  CALL stepon 

!-- Clean up

  CALL resetcolumn
  IF (iadvec == tpcore) CALL cleanup_tpcore
  IF (iadvec == semi_lagrangian) CALL cleanup_semi_lagrangian

  CALL cleanup_timer

  IF (lhd) THEN
    CALL cleanup_hydrology
!hd-hi switched off HD high resolution output
!hd-hi    CALL hd_highres_close
  ENDIF

  IF (numfl2 > 0) DEALLOCATE (field2)
  IF (numfl1 > 0) DEALLOCATE (field1)

  CALL free_memory
  CALL cleanup_scanbuffer
!---wiso-code
  IF (lwiso) THEN
    CALL cleanup_scanbuffer_wiso
  END IF
!---wiso-code-end
  CALL cleanup_o3clim
  CALL cleanup_clim
  CALL cleanup_sst
  CALL cleanup_grib
  CALL cleanup_legendre
  CALL cleanup_aero
  CALL cleanup_tmp_buffer
  CALL cleanup_gaussgrid
#ifdef FFT991
  CALL cleanup_fft991
#else
  CALL cleanup_fft992
#endif
  CALL cleanup_decomposition
  CALL dealloc_mods
  CALL cleanup_netcdf
  CALL cleanup_io
  CALL cleanup_geoloc
  CALL cleanstatr
  IF(lso4) CALL cleanup_so4
  IF (lnudge) THEN
     CALL NudgingInit(NDG_CLEAN_MEM)
     IF (lstop) CALL NudgingInit(NDG_CLOSE)
  END IF
  
  IF (lnmi .AND. lstop) CALL NMI_Close

  IF (ltdiag)  CALL DIAG_Init(IDIAG_FREE)
 
END SUBROUTINE control
