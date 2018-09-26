MODULE mod_nproc
!
! -- nproc.h   97-09-03   Version 2.2   Author: S. Valcke
!    *******
!              18-03-02   Version 2.5   changed in module
!@
!@  Contents :
!@  --------
!@
!@ -- nproc  : process ID of each model -used in PIPE, SIPC and GMEM
!@
!     -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: nproc

!     -------------------------------------------------------------------

END MODULE mod_nproc
