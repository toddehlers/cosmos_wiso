MODULE mo_semi_impl

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  ! ---------------------------------------------------------------
  !
  ! module *mo_semi_impl* - quantities needed for the semi-implicit scheme
  !
  ! ---------------------------------------------------------------

  REAL(dp) :: betadt        !  =0 : explicit scheme(for *d*,*t*,*alps*).
                            !  =1.:semi implicit scheme.
  REAL(dp) :: betazq        !  =0 : explicit scheme(for *v0*,*q*).
                            !  =1.:semi implicit scheme.
  REAL(dp) :: apr           !  reference surf. pressure semi-implicit scheme.
  REAL(dp) :: tr            !  reference temperature semi-implicit scheme.
  REAL(dp) :: ulat          !  latitude of maximum !u!+!v! (real winds).
  REAL(dp) :: uvmax         !  max(!u!+!v!) (real winds).
  REAL(dp) :: ulm           !  linearization wind profile on the latitude
                            !  and at the level of maximum !u!+!v!.
  REAL(dp) :: vcrit         !  critical velocity above which
                            !  horizontal diffusion is enhanced for
                            !  t63 with dt=20min.
  REAL(dp) :: vcheck        !  threshold value for check of high windspeed
  REAL(dp) :: hdamp         !  damping factor for strong
                            !  stratospheric damping.
  REAL(dp), ALLOCATABLE :: vmax(:)

  INTEGER :: nulev          !  *level of maximum !u!+!v! (real winds).

  REAL(dp) :: eps  = 0.1_dp !   time filtering coefficient.

END MODULE mo_semi_impl
