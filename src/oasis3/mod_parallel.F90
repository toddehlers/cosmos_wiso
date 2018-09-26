MODULE mod_parallel
!
! -- parallel.h   11-09-95   Version 2.0   Author: Laurent Terray
!    **********
!                 18-03-02   Version 2.5   changed in module
!@
!@  Contents : variables related to parallel data decomposition
!@  --------
!@
!@ -- nparal : parallel decomposition description (see CLIM manual) (2D)
!@
!@ -- cparal : type of parallel decomposition (idem) (1D)
!@
!     -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  INTEGER (kind=ip_intwp_p), DIMENSION(:,:), ALLOCATABLE :: nparal
!
  CHARACTER(len=8),DIMENSION(:),ALLOCATABLE :: cparal
!
!     -------------------------------------------------------------------
END MODULE mod_parallel



