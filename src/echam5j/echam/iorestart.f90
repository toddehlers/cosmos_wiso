SUBROUTINE iorestart

  ! Description:
  !
  ! Reads netCDF history files for a resumed run.
  !
  ! Method:
  !
  ! *iorestart* positions data sets at the beginning of a rerun,
  ! writing data description records, and setting up necessary work
  ! files.
  !
  ! Information is written to the data description records of
  ! appropriate files, and work files are written if necessary.
  !
  !
  ! Authors:
  !
  ! L. Kornblueh, MPI, May 1999, f90 rewrite
  ! U. Schulzweida, MPI, May 1999, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_tracer,        ONLY: xtini
  USE mo_jsbach_interface, ONLY: jsbach_restart
  USE mo_doctor,        ONLY: nout
  USE mo_io,            ONLY: IO_read_streams
  USE mo_memory_gl,     ONLY: xt
  USE mo_memory_g1a,    ONLY: xtm1
  USE mo_control,       ONLY: ltimer
  USE mo_timer,         ONLY: timer_start, timer_stop, timer_netcdf
  USE mo_mpi,           ONLY: p_pe, p_io

!---wiso-code
  USE mo_wiso,          ONLY: lwiso, lwiso_rerun
  USE mo_memory_base,   ONLY: ostreams, nstreams
  USE mo_exception,     ONLY: message_text, message, finish
!---wiso-code-end

  IMPLICIT NONE

!---wiso-code
  INTEGER :: i
!---wiso-code-end
  
  !  Executable statements 

  ! Restart from history files

!---wiso-code

  ! check if water isotope rerun files shall be used; 
  ! if not, set lrerun flag of ostream to .false. 
  
  IF (lwiso) THEN
    IF ( .NOT. lwiso_rerun) THEN
      DO i = 1,nstreams
         IF ( (ostreams(i)%name == 'wiso') .OR. (ostreams(i)%name == 'js_wiso') .OR. (ostreams(i)%name == 'sf_wiso') ) THEN
            ostreams(i)%lrerun = .FALSE.
         END IF
      END DO
    END IF
  END IF

!---wiso-code-end

  ! Read all restart files

  IF (ltimer) CALL timer_start(timer_netcdf)
  CALL IO_read_streams
  IF (ltimer) CALL timer_stop(timer_netcdf)

!---wiso-code

  ! check if *wiso* rerun files shall be used; 
  ! if not, initialize water isotope values
  ! and set back lrerun flag of wiso ostream to .true. for next iteration 
  
  IF (lwiso) THEN
    IF ( .NOT. lwiso_rerun) THEN
      CALL iniwiso
      DO i = 1,nstreams
         IF ( (ostreams(i)%name == 'wiso') .OR. (ostreams(i)%name == 'js_wiso') .OR. (ostreams(i)%name == 'sf_wiso') ) THEN
            ostreams(i)%lrerun = .TRUE.
         END IF
      END DO
    ELSE
      IF (p_pe == p_io) WRITE (nout,'(a,/,/)') 'Water isotopes resumed from *wiso* rerun file.'
    END IF
  END IF

!---wiso-code-end

  CALL xtini (xt, xtm1)

  CALL jsbach_restart

  IF (p_pe == p_io) &
       WRITE (nout,'(/)')

END SUBROUTINE iorestart
