SUBROUTINE inictl

  ! Description:
  !
  ! Preset constants in mo_control.
  !
  ! Method:
  !
  ! Calculate space for memory manager
  !
  ! *inictl* is called from *initialize*.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, August 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! H.-S. Bauer, MPI, Jul 1998, changed
  ! I. Kirchner, MPI, August 1998, tendency diagnostics
  ! L. Kornblueh, MPI, April 1998, added NWP forecast mode
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! M. Esch, MPI, July 1999, remove nudging, nmi switches
  ! M. Esch, MPI, July 1999, modifications for ECHAM5
  ! U. Schlese, DKRZ, December 1999, modifications for coupling
  ! I. Kirchner, MPI, October 2000, revision, time control
  ! L. Kornblueh, MPI, January 2001, revision of time control
  ! I. Kirchner, MPI, March 2001, revision
  ! A. Rhodin, MPI, June 2001, g3x,g4x fields removed
  ! L. Kornblueh, MPI, October 2001, added missing broadcast of nsub in runctl
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! M. Esch, MPI, September 2002, switch for mixed layer ocean
  ! L. Kornblueh, MPI, April 2003, switch for port test added
  ! U. Schulzweida, MPI, March 2007, added daily SST and SIC support
  ! S. Lorenz, MPI, November 2007, added switch for volcanic forcing
  !
  ! for more details see file AUTHORS
  !

  USE mo_exception,     ONLY: finish, message, message_text
  USE mo_mpi,           ONLY: p_io, p_parallel, p_parallel_io, p_bcast, p_pe
  USE mo_doctor,        ONLY: nin, nout, nerr
  USE mo_control,       ONLY: lamip, lcouple, lipcc, ldebug, lmidatm, lmlo,   &
                              lnmi, lnudge, lnwp, ltdiag, lhd, lhd_que, lso4, &
                              ldiagamip, lprint_m0, ldebugio, ldebugmem,      &
                              ldebughd, loldrerun, ltimer, ltctest,           &
                              ngl, nlev, nproca, nprocb, nproma,              &
                              numfl1, numfl2, nhd_diag, nflu, na2stre,        &
                              nisp, nigp, ndiahdf, nist,  nice, na2stat,      &
                              nhf1, nhg1, nhg2, nhg3, nhgl1, nfl1, nfl2,      &
                              njin, njout, subjob_cmd, ldailysst,             &
                              stdout_redir, stderr_redir, lsolc, lreff,       &
                              lco2, lco2_flxcor, lco2_2perc, ico2_emis,       &
                              lco2_nudge, co2_nudge_vmr, co2_nudge_tau,       &
                              lslm,l_volc

  USE mo_port_test,     ONLY: lport
  USE mo_namelist,      ONLY: open_nml, position_nml, nnml, POSITIONED
  USE mo_time_base,     ONLY: set_calendar_type, JULIAN, CYL360
  USE mo_time_control,  ONLY: delta_time, no_days, no_steps, no_cycles,       &
                              dt_start, dt_stop, dt_resume,                   &
                              labort, l_orbvsop87,                            &
                              putdata, putrerun, putocean, getocean,          &
                              puthd, gethd,                                   &
                              trigjob, nsub, subflag,                         &
                              p_bcast_event, NSUB_MAX, trigfiles,             &
                              lresume, ldebugev, lfirst_cycle
  USE mo_advection,     ONLY: iadvec
  USE mo_filename,      ONLY: out_datapath, out_expname,                      &
                              out_filetype, trac_filetype, rerun_filetype

  IMPLICIT NONE

  ! Local scalars: 

  LOGICAL :: ly360 = .FALSE.

  INTEGER ::  i

  INTEGER           :: ierr  ! error return value from position_nml
  CHARACTER(len=16) :: fname ! filename for output redirection

  INCLUDE 'runctl.inc'

  ! Executable statements 

  ! 1. Preset constants

!-- 2. Read namelist runctl

  IF (p_parallel_io) THEN

     ! This is the first namelist to read.
     ! Initialize namelist module.

     CALL open_nml ('namelist.echam', unit=nin)

     ! read namelist.

     CALL position_nml ('RUNCTL', status=ierr)
     SELECT CASE (ierr)
     CASE (POSITIONED)
       trac_filetype = 0
       READ (nnml, runctl)
       IF(trac_filetype == 0) trac_filetype = out_filetype
     END SELECT
  ENDIF
  IF (p_parallel) THEN
     CALL p_bcast (nisp, p_io)
     CALL p_bcast (nigp, p_io)
     CALL p_bcast (ndiahdf, p_io)
     CALL p_bcast (nist, p_io)
     CALL p_bcast (nice, p_io)
     CALL p_bcast (nflu, p_io)
     CALL p_bcast (na2stat, p_io)
     CALL p_bcast (na2stre, p_io)
     CALL p_bcast (nhf1, p_io)
     CALL p_bcast (nhg1, p_io)
     CALL p_bcast (nhg2, p_io)
     CALL p_bcast (nhg3, p_io)
     CALL p_bcast (nhgl1, p_io)
     CALL p_bcast (nfl1, p_io)
     CALL p_bcast (nfl2, p_io)
     CALL p_bcast (njin, p_io)
     CALL p_bcast (njout, p_io)

     CALL p_bcast (ltimer, p_io)
     CALL p_bcast (ldebugio, p_io)
     CALL p_bcast (ldebugmem, p_io)
     CALL p_bcast (ldebughd, p_io)
     CALL p_bcast (ldebugev, p_io)
     CALL p_bcast (loldrerun, p_io)
     CALL p_bcast (ltctest, p_io)

     CALL p_bcast (lresume, p_io)

     CALL p_bcast (subjob_cmd, p_io)

     CALL p_bcast (out_datapath, p_io)
     CALL p_bcast (out_expname,  p_io)
     CALL p_bcast (rerun_filetype, p_io)
     CALL p_bcast (out_filetype, p_io)
     CALL p_bcast (trac_filetype, p_io)

     CALL p_bcast (stdout_redir, p_io)
     CALL p_bcast (stderr_redir, p_io)

     CALL p_bcast (lprint_m0, p_io)
     CALL p_bcast (numfl1, p_io)
     CALL p_bcast (numfl2, p_io)
     CALL p_bcast (ldebug, p_io)
     CALL p_bcast (ldailysst, p_io)
     CALL p_bcast (lamip, p_io)
     CALL p_bcast (ldiagamip, p_io)
     CALL p_bcast (labort, p_io)
     CALL p_bcast (lnwp, p_io)
     CALL p_bcast (lnudge, p_io)
     CALL p_bcast (lmidatm, p_io)
     CALL p_bcast (lmlo, p_io)
     CALL p_bcast (iadvec, p_io)
     CALL p_bcast (lnmi, p_io)
     CALL p_bcast (ltdiag, p_io)
     CALL p_bcast (lport, p_io)
     CALL p_bcast (nproma, p_io)
     CALL p_bcast (nproca, p_io)
     CALL p_bcast (nprocb, p_io)
     CALL p_bcast (lcouple, p_io)
     CALL p_bcast (lipcc, p_io)
     CALL p_bcast (lso4, p_io)
     CALL p_bcast (lsolc, p_io)
     CALL p_bcast (lreff, p_io)
     CALL p_bcast (lco2, p_io)
     CALL p_bcast (lco2_flxcor, p_io)
     CALL p_bcast (lco2_2perc, p_io)
     CALL p_bcast (ico2_emis, p_io)
     CALL p_bcast (lco2_nudge, p_io)
     CALL p_bcast (co2_nudge_vmr, p_io)
     CALL p_bcast (co2_nudge_tau, p_io)
     CALL p_bcast (lslm, p_io)

     CALL p_bcast (l_orbvsop87, p_io)
     CALL p_bcast (ly360, p_io)
     CALL p_bcast (delta_time, p_io)
        
     CALL p_bcast (dt_start, p_io)
     CALL p_bcast (dt_resume, p_io)
     CALL p_bcast (dt_stop, p_io)

     CALL p_bcast (no_days, p_io)
     CALL p_bcast (no_cycles, p_io)
     CALL p_bcast (no_steps, p_io)

     CALL p_bcast (nsub, p_io)

     CALL p_bcast_event(putdata, p_io)
     CALL p_bcast_event(trigfiles, p_io)
     CALL p_bcast_event(putrerun, p_io)

     nsub = MIN(MAX(nsub,0),NSUB_MAX)
     CALL p_bcast(subflag, p_io)
     DO i=1,nsub
        CALL p_bcast_event(trigjob(i), p_io)
     END DO
     CALL p_bcast_event(putocean, p_io)
     CALL p_bcast_event(getocean, p_io)
     CALL p_bcast_event(puthd, p_io)
     CALL p_bcast_event(gethd, p_io)

     CALL p_bcast (lhd, p_io)
     CALL p_bcast (lhd_que, p_io)
     CALL p_bcast (nhd_diag, p_io)

     CALL p_bcast (l_volc, p_io)

  ENDIF

  ! reset lresume for the second rerun cycle during an initial run
  IF (.NOT. lfirst_cycle) lresume = .TRUE.

  ! redirect output

  IF ( stderr_redir == p_pe+1                         .or. & 
     (-stderr_redir /= p_pe+1 .and. stderr_redir < 0)) THEN
    write (fname,'("echam5_pe",i5.5,".e")') p_pe
    write (nerr,*) ' stderr redirected to file ',fname,' on PE ',p_pe
    close (nerr)
    open  (nerr,file=fname)
    write (nerr,*) ' stderr redirected to file ',fname,' on PE ',p_pe
  ENDIF

  IF ( stdout_redir == p_pe+1                         .or. & 
     (-stdout_redir /= p_pe+1 .and. stdout_redir < 0)) THEN
    write (fname,'("echam5_pe",i5.5,".o")') p_pe
    write (nout,*) ' stdout redirected to file ',fname,' on PE ',p_pe
    close (nout)
    open  (nout,file=fname)
    write (nout,*) ' stdout redirected to file ',fname,' on PE ',p_pe
  ENDIF

  ! check consistency of orbit and year length

  IF (p_parallel_io) THEN
    IF (l_orbvsop87 .AND. ly360) THEN
      CALL finish('inictl', &
           ' ly360=.TRUE. cannot run with real orbit (l_orbvsop87=.TRUE.).')
    ENDIF
  ENDIF
  IF (ly360) THEN
    CALL set_calendar_type(CYL360)
  ENDIF

  IF ( nproma < 0 ) nproma = 0  ! default is lon/lat ordering

!lk  ! correct number of rerun cycles
!lk  no_cycles = MAX(no_cycles,1)
  IF (no_cycles < 1) THEN
    no_cycles = 0
    CALL message('inictl', &
         ' Debug mode - DT_STOP is finishing this specific run.')
  END IF

  ! correct maximal number of subjobs
  nsub = MIN(MAX(nsub,0),NSUB_MAX)

  IF (lamip .AND. ldailysst) THEN
    CALL finish('inictl', &
           ' ldailysst=.TRUE. does not work with lamip=.TRUE.')
  END IF

!-- 4.2 Write *runctl.*

   IF (.NOT. p_parallel) THEN 
     WRITE (nout,runctl)
   END IF

END SUBROUTINE inictl
