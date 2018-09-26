SUBROUTINE geopot(phi,ptv,plnpr,palpha,phis,kdim,klen)

  ! Description:
  !
  ! Calculates full- or half-level geopotentials.
  !
  ! Method:
  !
  ! Integrate the hydrostatic equation in the vertical
  ! to obtain full- or half-level values of the geopotential.
  !
  ! *geopot* is called during the calculation of adiabatic
  ! tendencies, prior to the calculation of the physical paramet-
  ! erizations, and in the post-processing. 
  !
  ! Parameters are
  !   *phi*         *computed geopotential.
  !   *ptv*         *virtual temperature.
  !   *plnpr*       *logarithm of ratio of pressures, computed
  !                  by *auxhyb*.
  !   *palpha*      *for full-level values use *alpha* as computed by *auxhyb*.
  !                  Set *alpha* to zero for values at half-levels below
  !                  the full levels.
  !   *phis*        *surface geopotential.
  !   *kdim*        *first dimension of 2-d arrays *phi,*
  !                  *ptv,* *plnpr,* and *palpha.*
  !   *klen*        *number of points for which calculation is performed.
  !
  ! Required constants are obtained from module *mo_hyb*. 
  ! The latter should have been initialized by a call of subroutine *inihyb*.
  !
  ! Results are computed for *klen* consecutive points for *nlev* levels.
  !
  ! The choice of full- or half-level values is determined
  ! by the specification of the input array *alpha.*
  !
  ! External documentation of the model equations and the
  ! organization of the vertical calculation.
  !
  ! Authors:
  !
  ! A. J. Simmons, ECMWF, January 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,      ONLY: dp
  USE mo_control,   ONLY: nlev
  USE mo_hyb,       ONLY: nlevm1

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER ,INTENT(in) :: kdim, klen

  !  Array arguments 
  REAL(dp) ,INTENT(in)  :: palpha(kdim,*), phis(*), plnpr(kdim,*)      &
                         , ptv(kdim,*)
  REAL(dp) ,INTENT(inout) :: phi(kdim,*)

  !  Local scalars: 
  INTEGER :: ikp, jk, jl


  !  Executable statements 

!-- 1. Integrate hydrostatic equation

  DO jl = 1, klen
    phi(jl,nlev) = palpha(jl,nlev)*ptv(jl,nlev) + phis(jl)
  END DO

  DO jk = nlevm1, 1, -1
    ikp = jk + 1

    DO jl = 1, klen
      phi(jl,jk) = palpha(jl,jk)*ptv(jl,jk) +                          &
                                       (plnpr(jl,ikp)-palpha(jl,ikp))* &
                                         ptv(jl,ikp) + phi(jl,ikp)

    END DO
  END DO

  RETURN
END SUBROUTINE geopot
