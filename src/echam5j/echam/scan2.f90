SUBROUTINE scan2

  ! Description:
  !
  ! 2nd scan over the latitude lines controlling the inverse Legendre
  !
  ! Method:
  !
  ! This subroutine scans over the latitude lines to perform
  ! inverse *legendre transforms(*lti*).
  !
  ! *scan2* is called from *stepon*
  !
  ! The spectral components are located in long term storage
  ! arrays. *posts2* allocate and release buffers,
  ! and perform other tasks associated with the i/o.
  !
  ! Externals:
  ! *lti*       inverse legendre transforms
  ! *posts2*    input/output f-buffer
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, January 1995, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_call_trans,    ONLY: spectral_to_legendre

  IMPLICIT NONE

  !  External subroutines 

  EXTERNAL lti

  !  Executable statements 

!-- Transpose: spectral space -> Legendre space

  call spectral_to_legendre

!-- Inverse *Legendre transforms

  CALL lti

END SUBROUTINE scan2
