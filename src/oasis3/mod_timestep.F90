MODULE mod_timestep
!
! -- timestep.h   29-08-95   Version 2.0beta   Author: Laurent Terray
!    **********   01-02-96   Version 2.0 : addition of nitfn (L. Terray)
!                 18-03-02   Version 2.5   changed in module
!@
!@  Contents : variables controlling timestepping
!@  --------
!@
!@   OASIS information :
!@   -----------------
!@
!@ -- ntime : total simulated time (in seconds)
!@
!@ -- niter : number of iterations
!@
!@ -- nitfn : last iteration number
!@
!@ -- nstep : timestep value (in seconds)
!@
!@   Remote models information :
!@   -------------------------
!@
!@ -- mstep : number of time steps (1D)
!@
!@ -- mfcpl : frequency of coupling in timesteps (1D)
!@
!@ -- mdt   : length of timesteps in seconds (1D)
!@
!     -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  INTEGER (kind=ip_intwp_p) :: ntime, niter, nitfn, nstep
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: mstep, mfcpl, mdt
!
!     -------------------------------------------------------------------
!
END MODULE mod_timestep



