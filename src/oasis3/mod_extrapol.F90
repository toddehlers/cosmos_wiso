MODULE mod_extrapol
!
! -- extrapol.h   16-12-97   Version 2.2   Author: Laurent Terray
!    **********   30-03-99   Version 2.3   READ/WRITE flag, FILE and dataset
!                                           index for NINENN weights
!                 18-03-02   Version 2.5   changed in module
!@
!@  Contents : variables and arrays related to extrapolation
!@  --------
!@
!@ WEIGHT --->>>
!@
!@ -- aextra : weight for each gcm-1 mesh proportional to overlapped area (1D)
!@
!@ -- nextra : neighbors adress on gcm-1 grid for a given gcm-2 grid point (1D)
!@
!@ -- lextra : I/O initialization flag for each field
!@
!@ NINENN --->>>
!@
!@ -- niwtn : flag to read/write EXTRAP/NINENN parameters (1D)
!@
!@ -- niwtng : flag to read/write EXTRAP/NINENN parameters when extrap
!@             is called by GLORED (1D)
!@
!@ -- cwninenn : file name for NINENN parameter FILE
!@
!@ -- nninnfl : flag to identify different EXTRAP/NINENN parameter sets 
!@              within all NINENN/EXTRAP analyses (1D)
!@
!@ -- nninnflg : flag to identify different EXTRAP/NINENN parameter sets 
!@               WHEN extrap is called by GLORED within all NINENN/EXTRAP
!@               analyses (1D)
!@
!@ COMMON --->>>
!@
!@ -- lweight : flag indicating IF EXTRAP/NINENN parameter sets have 
!@              already been calculated or read (.TRUE.) or not (.FALSE.)
!@
!@ -- lextrapdone : logical indicating if EXTRAP analysis has been done 
!@                  for field
! -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  USE mod_parameter
!
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: nextra
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: niwtn, nninnfl 
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: niwtng, nninnflg
  INTEGER (kind=ip_intwp_p), DIMENSION(:,:), ALLOCATABLE :: iincre
  INTEGER (kind=ip_intwp_p), DIMENSION(:,:,:), ALLOCATABLE :: iaddress
!
  REAL (kind=ip_realwp_p), DIMENSION(:), ALLOCATABLE :: aextra
  REAL (kind=ip_realwp_p), DIMENSION(:,:,:), ALLOCATABLE :: zweights
!
  LOGICAL, DIMENSION(:), ALLOCATABLE :: lextra, lweight
  LOGICAL, DIMENSION(:), ALLOCATABLE :: lextrapdone
!
  CHARACTER(len=8), PARAMETER ::  cwninenn = 'nweights'
!
!
! -------------------------------------------------------------------
!
END MODULE mod_extrapol


