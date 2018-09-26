MODULE mod_smooth
!
! -- smooth.h   01-10-95   Version 2.0   Author: Laurent Terray
!
!@
!@  Contents : variables controlling how smoothing is handled
!@  --------
!@
!@ Common to all N* variables :
!@       indices for the gcm-1 grid of regions corresponding
!@       to the boundaries of the gcm-2 domain.
!@
!@ -- nsltb : southern outermost point of the smoothing region
!@
!@ -- nslte : southern innermost point of the smoothing region
!@
!@ -- nnltb : northern outermost point of the smoothing region
!@
!@ -- nnlte : northern innermost point of the smoothing region
!@
!@ -- qalfa : control weights in the transition zone (north and south)
!@
!@ -- qbeta : control weights in the transition zone (west and east)
!@
!@ -- nliss : width of the smoothing region (west and east boundaries)
!@
!@ -- nwlgmx : longitude index of the furthest east point on the western  
!@             boundary
!@
!@ -- nelgmx : longitude indexof the furthest west point on the eastern  
!@             boundary
!@
! -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  INTEGER(kind=ip_intwp_p), PARAMETER :: nsltb = 18
  INTEGER(kind=ip_intwp_p), PARAMETER :: nslte = 26
  INTEGER(kind=ip_intwp_p), PARAMETER :: nnltb = 49
  INTEGER(kind=ip_intwp_p), PARAMETER :: nnlte = 41
  INTEGER(kind=ip_intwp_p), PARAMETER :: qalfa = 0.125
  INTEGER(kind=ip_intwp_p), PARAMETER :: qbeta = 0.25
  INTEGER(kind=ip_intwp_p), PARAMETER :: nliss = 5
  INTEGER(kind=ip_intwp_p), PARAMETER :: nelgmx = 229
  INTEGER(kind=ip_intwp_p), PARAMETER :: nwlgmx = 57
!
END MODULE mod_smooth
! -------------------------------------------------------------------
