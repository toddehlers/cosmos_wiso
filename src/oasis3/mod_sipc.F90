MODULE mod_sipc
!
! -- sipc.h   97-08-11   Version 2.2   Author: S.Valcke,A.Piacentini
!    ******
!             18-03-02   Version 2.5   changed in module
!@
!@  Contents : variables describing pools (set of shared memory segments) 
!@  --------
!@
!@ -- mpoolinit(r/w) : handles associated to model pools for passing initial 
!@                     information(1 for reading and 1 for writing)
!@
!@ -- mpoolidin : handles associated to pools used to pass field to oasis
!@
!@ -- mpoolidou : handles associated to pools used to pass field from oasis
!@ 
!     -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: mpoolidin, mpoolidou
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: mpoolinitr, mpoolinitw
!
!     -------------------------------------------------------------------
END MODULE mod_sipc
