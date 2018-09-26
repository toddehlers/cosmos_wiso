
SUBROUTINE datim(kdate,ktime)

  ! Description:
  !
  ! Return date and time in ECMWF format.
  !
  ! Method:
  !
  ! *datim* returns the date and time in the standard *ECMWF*
  ! format, as used by the operational system.
  !
  ! *kdate*   - date returned as *yyyymmdd* (integer).
  ! *ktime*   - time returned as *hhmmss* (integer).
  !
  ! Authors:
  !
  ! J. K. Gibson, ECMWF, February 1983, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kdate, ktime

  !  Local scalars: 
  CHARACTER  (8) :: ydate
  CHARACTER (10) :: ytime


  !  Executable statements 

  ! Obtain date and time

  CALL date_and_time(ydate, ytime)

  READ (ydate,'(I8)') kdate
  READ (ytime,'(I6)') ktime

  RETURN
END SUBROUTINE datim
