MODULE mod_coast
!
! -- coast.h   01-11-95   Version 2.0   Author: Laurent Terray
!    *******
!@
!@  Contents : variables and arrays related to the coast correction
!@  --------
!@
!@ Anaism --->>>
!@
!@    Rappel:
!@ -- nmesh : number of ocean gcm-1 squares overlapped by a given gcm-2 
!@            square (1D)
!@    End of rappel
!@
!@ -- ncoast : number of coastal ocean squares on the gcm-2 grid for which
!@              there are no underlying ocean squares on the gcm-1 grid
!@
!@ -- npcoast : coast mismatch data array (2D)
!@              (n,1) --> 1D index of nth point described above
!@              (n,2) --> number of neighbours suitable for extrapolating SST
!@              (n,3-6) --> 1D indices of neighbours of nth point
!@ 
!@ -- nfcoast : flag to perform coast mismatch correction
!@
!@ -- lcoast  : initialization flag
!@
! -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  INTEGER (kind=ip_intwp_p), DIMENSION(:,:), ALLOCATABLE :: npcoast
  INTEGER (kind=ip_intwp_p) :: ncoast, nfcoast
!
  LOGICAL :: lcoast
!
! -------------------------------------------------------------------
!
END MODULE mod_coast
