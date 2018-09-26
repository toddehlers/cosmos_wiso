MODULE mod_rainbow
!
! -- rainbow.h   05-08-96   Version 2.1   Author: Laurent Terray
!    *********
!@
!@  Contents : variables and arrays related to mapping and subgrid interpolator
!@  --------
!@
!@ Mapping --->>>
!@
!@ -- amapp : weight for each gcm-1 mesh proportional to overlapped area (1D)
!@
!@ -- nmapp : neighbors adress on gcm-1 grid for a given gcm-2 grid point (1D)
!@
!@ -- lmapp : I/O initialization flag for each field
!@
!@ Subgrid --->>>
!@
!@ -- asubg : weight for each gcm-2 mesh proportional to overlapped area (1D)
!@
!@ -- nsubg : neighbors adress on gcm-2 grid for a given gcm-1 grid point (1D)
!@
!@ -- lsubg : I/O initialization flag for each field
!@
!     -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  USE mod_parameter
!
  INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: nmapp, nsubg
!
  REAL (kind=ip_realwp_p), DIMENSION(:),ALLOCATABLE :: amapp, asubg
!
  LOGICAL,DIMENSION(:),ALLOCATABLE :: lmapp, lsubg
!
!     -------------------------------------------------------------------
!
END MODULE mod_rainbow
