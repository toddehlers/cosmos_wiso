SUBROUTINE presf(pf,kdimp,ph,klen)

  ! Description:
  !
  ! Compute full-level pressures from half-level values.
  !
  ! Method:
  !
  ! Full-level pressures are defined as the arithmetic
  ! average of the two adjoining half-level pressures.
  !
  ! *presf* is called from *physc*. 
  ! Parameters are:
  !    *pf*        *computed full-level pressures.
  !    *ph*        *half-level pressures.
  !    *kdimp*     *first dimension of 2-d arrays *pf* and *ph*
  !    *klen*      *number of points for which computation is
  !                 performed.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, May 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,    ONLY: dp
  USE mo_control, ONLY: nlev

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER ,INTENT(in) :: kdimp, klen

  !  Array arguments 
  REAL(dp),INTENT(in) :: ph(kdimp, *)
  REAL(dp),INTENT(inout):: pf(kdimp, *)

  !  Local scalars: 
  INTEGER :: jlev, jlon


  !  Executable statements 

!-- 1. Compute full-level pressure values

  DO jlev = 1, nlev
    DO jlon = 1, klen
      pf(jlon,jlev) = (ph(jlon,jlev)+ph(jlon,jlev+1))*.5_dp
    END DO
  END DO

  RETURN
END SUBROUTINE presf
