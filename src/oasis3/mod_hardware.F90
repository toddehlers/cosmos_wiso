MODULE mod_hardware
!
! -- hardware.h   25-07-95   Version 2.0   Author: Laurent Terray
!    **********   24-12-97   Version 2.2   addition of nsigxxx
!                 23-10-98   Version 2.3   addition of ctype (J.Latour F.S.E.)
!@
!@  Contents : variables related to machine and message passing type
!@  --------
!@
!@ -- cchan : type of message passing used to run OASIS (PIPE, SVIPC, CLIM)
!@
! -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  CHARACTER(len=4) :: cchan
!
END MODULE mod_hardware
!
! -------------------------------------------------------------------
