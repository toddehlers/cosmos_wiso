MODULE mod_comclim
!
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*    =comclim.h=  CLIM 1.1 internal include file
!                  Coupling Library for Interfacing Models
!
!     - - - - - - - - - - - - - - - - - - - - - - - - - - -
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
  USE mod_kinds_oasis
  USE mod_clim
!
!-----Unit number for trace file
!
  INTEGER(kind=ip_intwp_p)	:: nulprt
!
!-----Models descriptors
!
!     ncplprocs: total number of processes involved in the coupling
!                counting oasis process (1)
  INTEGER(kind=ip_intwp_p) :: ncplprocs, mynum, mytid
  INTEGER(kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: modtid
!  INTEGER*4 :: mynummod
!
!-----Minimum exchange frequency, total time of the simulation, total 
!     number of fields exchanged in the simulation and local communicator of 
!     each component.
! 
  INTEGER(kind=ip_intwp_p) :: ig_frqmin, ig_ntime, ig_clim_nfield, ig_local_comm
!
!-----Ports descriptors
!
  INTEGER(kind=ip_intwp_p) :: nports 
  INTEGER(kind=ip_intwp_p), DIMENSION(:,:), ALLOCATABLE :: myport, irport
  INTEGER(kind=ip_intwp_p), DIMENSION(:,:), ALLOCATABLE :: mydist, irdist
!
!
!-----Links descriptors
!
  INTEGER(kind=ip_intwp_p) :: nlinks 
  INTEGER(kind=ip_intwp_p), DIMENSION(:,:), ALLOCATABLE :: mylink
!
!-----Data encoding
!
  INTEGER(kind=ip_intwp_p) :: mycode
  INTEGER(kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: ncode
!
!-----Inquiry descriptors
!
  INTEGER(kind=ip_intwp_p) :: nbsend, nbrecv
!
!-----Time out stuff
!
  INTEGER(kind=ip_intwp_p) :: ntiret, ntiogp, ntiout
!
!-----Timesteps descriptors
!
  INTEGER(kind=ip_intwp_p) :: mystep, mystdt, myfcpl
  INTEGER(kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: ig_nstep
  INTEGER(kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: nstdt
  INTEGER(kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: nfcpl 
!
!-----Character strings
!
  CHARACTER*8 :: cgroup
  CHARACTER*12 :: cnaprt
  CHARACTER*32, DIMENSION(:), ALLOCATABLE :: cports, clrport
  CHARACTER*32 :: cmynam
  CHARACTER*32, DIMENSION(:), ALLOCATABLE :: cnames  
  CHARACTER*4 :: ctype
!
!-----Packing area for Export and Import
!
  REAL(kind=ip_realwp_p), DIMENSION(:), ALLOCATABLE :: pkwork  
!
!
!-----CLIM/MPI2 specific variable
!-----mpi_rank :  Rank of process in MPI communicator
!-----mpi_size :  Size of communicator
!-----mpi_err :  Error code of MPI calls
!-----mpi_comm :  Communicator handle

  INTEGER(kind=ip_intwp_p) :: mpi_rank, mpi_size, mpi_err, mpi_comm

!
!-----Buffer for MPI_BSend
!
  REAL(kind=ip_double_p), DIMENSION(:), ALLOCATABLE :: dg_bufsend
!
!
#endif
END MODULE mod_comclim


