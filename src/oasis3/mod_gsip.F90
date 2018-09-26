MODULE mod_gsip
!
! -- mod_gsip  
!    ********
!             18-10-04    S. Valcke    Created      
!@
!@  Contents : variables Gossip channels
!@  --------
!@
!@ -- iga_gsipw : handles associated to model channel for writing 
!@                info and fields to each process of each models 
!@ -- iga_gsipr : handles associated to model channel for reading 
!@               info and fields from each process of each models
!     -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: iga_gsipw, iga_gsipr
!
!     -------------------------------------------------------------------
END MODULE mod_gsip
