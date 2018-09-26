MODULE mod_anais
!
! -- anais.h   01-11-95   Version 2.0   Author: Laurent Terray
!    *******             
!              31-08-96   Version 2.1 : add new fields for anais package
!
!              18-03-02   Version 2.5 : changed in module
!@
!@  Contents : variables and arrays related to the ANAIS interpolator
!@  --------
!@
!@ Anaism --->>>
!@
!@ -- nmesh : number of ocean gcm-1 squares overlapped by a given gcm-2 
!@            square (1D)
!@
!@ -- amint : weight for each gcm-1 mesh proportional to overlapped area (1D)
!@
!@ -- nmint : neighbors adress on gcm-1 grid for a given gcm-2 grid point (1D)
!@
!@ -- naismfl : flag to identify different ANAISM parameter sets (1D)
!@
!@ -- naismvoi : maximum number of overlapped neighbors (1D)
!@
!@ -- niwtm : flag to read/write ANAISM parameters (1D)
!@
!@ -- cwanaism : file name for ANAISM parameter file  
!@
!@ Anaisg --->>>
!@
!@ -- agint : weight for each gcm-1 grid point with gaussian distribution (1D)
!@
!@ -- ngint : neighbors adress on gcm-1 grid for a given gcm-2 grid point (1D)
!@
!@ -- varmul : variance multiplicator
!@
!@ -- naisgfl : flag to identify different ANAISG parameter sets (1D)
!@
!@ -- naisgvoi :  maximum number of used neighbors (1D)
!@
!@ -- niwtg : flag to read/write ANAISG parameters (1D)
!@
!@ -- cwanaisg : file name for ANAISG parameter file
!@
!@ Common --->>>
!@
!@ -- linit : I/O initialization flag for each field
!@
!@ -- cnaisout : ANAIS output file name
!@
!     -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  USE mod_parameter
!
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: ngint,nmint,nmesh
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: naismfl, naisgfl
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: naismvoi, naisgvoi
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: niwtm, niwtg
!
  REAL (kind=ip_realwp_p), DIMENSION(:), ALLOCATABLE :: agint, amint, varmul
!
  LOGICAL, DIMENSION(:), ALLOCATABLE :: linit
!
  CHARACTER(len=8), PARAMETER :: cwanaisg = 'gweights' 
  CHARACTER(len=8), PARAMETER :: cwanaism = 'mweights'
  CHARACTER(len=8), PARAMETER :: cnaisout = 'anaisout'
!
!   -------------------------------------------------------------------
END MODULE mod_anais
