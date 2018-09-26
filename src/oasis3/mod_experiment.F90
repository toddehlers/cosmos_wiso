MODULE mod_experiment
!
! -- experiment.h   26-08-95   Version 2.0   Author: Laurent Terray
!    ************   14-12-97   Version 2.2   addition of info mode (lmodinf)
!                   14-03-99   Version 2.3   CHARACTER*4 cjobnam (S. Valcke)
!                   07-09-00   Version 2.5   INTEGER nbcplproc (S. Valcke)
!                                            INTEGER nbtotproc
!                                            CHARACTER*80 cmpiarg
!                   18-03-02   Version 2.5   changed in module    
!@
!@  !ontents : variables related to the simulation being performed
!@  --------
!@
!@ -- cjobnam : experiment name
!@
!@ -- cmodnam : models name 
!@
!@ -- nmodel : number of models being coupled
!@
!@ -- nmseq : number of sequential models
!@
!@ -- lmodinf : information mode (extended header)
!@
!@ -- nbcplproc: number of model processes implied in the coupling
!@
!@ -- nbtotproc: total number of model processes
!@
!@ -- cmpiarg: launching argument for each model
!
! -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  INTEGER (kind=ip_intwp_p) :: nmseq
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: nbcplproc, nbtotproc
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: iga_unitmod
!
  CHARACTER(len=4) :: cjobnam
  CHARACTER(len=6), DIMENSION(:), ALLOCATABLE :: cmodnam
  CHARACTER(len=80), DIMENSION(:), ALLOCATABLE :: cmpiarg
!
  LOGICAL :: lmodinf
!
!     -------------------------------------------------------------------
!
END MODULE mod_experiment



