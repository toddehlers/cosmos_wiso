MODULE mod_memory
!
! -- memory.h   23-08-95   Version 2.0   Author: Laurent Terray
!    ********   
!               18-03-98   Version 2.5   changed in module
!@
!@  Contents : variables related to pseudo-dynamic memory allocation
!@  --------
!@
!@ -- nsizold : size for each field sub-array (before interpolation) (1D)
!@
!@ -- nsiznew : size for each field sub-array (after interpolation) (1D)
!@
!@ -- nadrold : adress for each field sub-array (before interpolation) (1D)
!@
!@ -- nadrold_grid : adress for each grid sub-array (before interpolation) (1D)
!@
!@ -- nadrnew : adress for each field sub-array (after interpolation) (1D)
!@
!@ -- nadrnew_grid : adress for each grid sub-array (after interpolation) (1D)
!@
!@ -- mskold : macro array for masks (before interpolation) (1D)
!@
!@ -- msknew : macro array for masks (after interpolation) (1D) 
!@
!@ -- nwork : integer work array (1D)
!@
!@ -- fldold : macro array for fields (before interpolation) (1D)
!@
!@ -- fldnew : macro array for fields (after interpolation) (1D)
!@
!@ -- xgrold : macro array for longitudes (before interpolation) (1D)
!@
!@ -- xgrnew : macro array for longitudes (after interpolation) (1D)
!@
!@ -- ygrold : macro array for latitudes (before interpolation) (1D)
!@
!@ -- ygrnew : macro array for latitudes (after interpolation) (1D)
!@
!@ -- surold : macro array for mesh surfaces (before interpolation) (1D) 
!@
!@ -- surnew : macro array for mesh surfaces (after interpolation) (1D)
!@
!@ -- work : real work array (1D)
!@
!     -------------------------------------------------------------------   
!
  USE mod_kinds_oasis
  USE mod_parameter
!
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: nsizold, nsiznew
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: nadrold, nadrnew
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: nadrold_grid
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: nadrnew_grid
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: mskold, msknew
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: nwork
!
  REAL (kind=ip_realwp_p), DIMENSION(:), ALLOCATABLE :: fldold, fldnew
  REAL (kind=ip_realwp_p), DIMENSION(:), ALLOCATABLE :: xgrold, xgrnew
  REAL (kind=ip_realwp_p), DIMENSION(:), ALLOCATABLE :: ygrold, ygrnew
  REAL (kind=ip_realwp_p), DIMENSION(:), ALLOCATABLE :: surold, surnew
  REAL (kind=ip_realwp_p), DIMENSION(:), ALLOCATABLE :: work
!
!     -------------------------------------------------------------------   
!
END MODULE mod_memory



