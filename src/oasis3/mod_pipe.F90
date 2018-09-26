MODULE mod_pipe
!
! -- pipe.h   25-09-95   Version 2.0   Author: Laurent Terray
!    ******   20-09-97   Version 2.2   Mods: suppress nproc (S. Valcke)
!             18-03-02   Version 2.5   changed in module
!@
!@  Contents : variables describing pipes (FIFO) used for message passing
!@  --------
!@
!@ -- cprnam : files associated to reading model pipes (1D)
!@
!@ -- cpwnam : files associated to writing model pipes (1D)
!@
!     -------------------------------------------------------------------
!
  CHARACTER(len=8), DIMENSION(:), ALLOCATABLE :: cprnam, cpwnam
!
!     -------------------------------------------------------------------
END MODULE mod_pipe
