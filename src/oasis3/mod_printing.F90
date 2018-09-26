MODULE mod_printing
!
! -- printing.h   30-04-99   Version 2.3   Author: Sophie Valcke
!    **********  
!@
!@  Contents : variables related to printing level in cplout
!@  --------
!@
!@ -- nlogprt : printing level in output file cplout: 0 = no printing
!@      1 = main routines and field names when treated, 2 = complete output
!@
! -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  INTEGER(kind=ip_intwp_p) :: nlogprt
!
END MODULE mod_printing
! -------------------------------------------------------------------




