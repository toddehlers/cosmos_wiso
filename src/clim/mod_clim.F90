MODULE mod_clim
!
! -- clim.h   18-08-95   Version 2.0   Author: Laurent Terray
!    ******
!             26-10-99   Version 2.4   Jean Latour (F.S.E.) MPI-2 support
!@
! -- mod_clim.f90  12-06-02  Version 2.5   A. Caubel 
!    ************
!
!@  Contents : variables related to the CLIM library
!@  --------
!@ For complete definition, see the CLIM manual
!@
!
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
!
  USE mod_kinds_oasis
!
  INTEGER(kind=ip_intwp_p) :: CLIM_MaxTag
!
!-----Parameter sizes
!
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Void = 0  
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_MaxSegments = 338 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_ParSize = 2*CLIM_MaxSegments+2 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Clength = 32 
!
!-----Dimension of buffer for packing / unpacking messages with MPI
!
  INTEGER(kind=ip_intwp_p) :: ig_maxtype, ig_CLIMmax
!
!-----Logical for buffered send
!
  LOGICAL :: lg_bsend
!
!-----Ports status
!
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_In      = 1 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Out     = 0 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_InOut   = 2 
!
!----Field local transformation
!
  INTEGER(kind=ip_intwp_p), PARAMETER :: ip_instant = 1
  INTEGER(kind=ip_intwp_p), PARAMETER :: ip_average = 2
  INTEGER(kind=ip_intwp_p), PARAMETER :: ip_accumul = 3
  INTEGER(kind=ip_intwp_p), PARAMETER :: ip_min = 4
  INTEGER(kind=ip_intwp_p), PARAMETER :: ip_max = 5
!
!-----Parallel distribution
!
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Strategy = 1 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Segments = 2 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Serial   = 0 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Apple    = 1 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Box      = 2 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Orange   = 3 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Offset   = 2 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Length   = 3 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_SizeX    = 3 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_SizeY    = 4 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_LdX      = 5 
!
!-----Datatypes
!
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Integer = 1 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Real    = 4  
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Double  = 8 
!
!-----Quit parameters
!
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_ContPvm = 0 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_StopPvm = 1 
!
!-----Error Codes
!
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_MaxCodes  = -22 
!
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Ok	 = 0 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_FastExit  = -1 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_BadName   = -2 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_BadPort   = -3 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_BadType   = -4 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_DoubleDef = -5 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_NotStep   = -6 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_IncStep   = -7 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_IncSize   = -8 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_NotClim   = -9 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_TimeOut   = -10 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Pvm       = -11 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_FirstCall = -12 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_PbRoute   = -13 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Group     = -14 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_BadTaskId = -15 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_NoTask    = -16 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_InitBuff  = -17 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Pack      = -18 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Unpack    = -19 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Down      = -20 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_PvmExit   = -21 
  INTEGER(kind=ip_intwp_p), PARAMETER :: CLIM_Mpi       = -22 
!
#endif
END MODULE mod_clim

