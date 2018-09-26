MODULE mod_string
!
! -- string.h   26-07-95   Version 2.0   Author: Laurent Terray
!    *******    25-09-96   Version 2.1   Addition of extra time step (nfend)
!                                        and field integral flag (nintflx)
!               12-11-97   Version 2.2   Addition of nmxdel
!               05-09-00   Version 2.5   Remove integral flag nintflx   
!               18-03-02   Version 2.5   changed in module
!@
!@  Contents : set of self consistent strings for each field
!@  --------
!@
!@ First line of SSCS :
!@ ------------------
!@
!@ -- cnaminp : symbolic name of input indirect fields 
!@
!@ -- cnamout : symbolic name of output indirect fields
!@
!@ -- cg_input_field : symbolic name of all input fields
!@
!@ -- cg_output_field : symbolic name of all output fields 
!@
!@ -- numlab : label number of exchanged fields indirectly
!@
!@ -- ig_numlab : label number of all exchanged fields
!@
!@ -- nfexch : coupling frequency of exchanged indirect fields (in seconds)
!@
!@ -- ig_freq : coupling frequency of all exchanged fields (in seconds)
!@
!@ -- ig_ntrans : number of analysis performed for each field
!@
!@ -- cficinp : restart file name of indirect input fields 
!@
!@ -- cg_restart_file : restart file name of all input fields  
!@
!@ -- ig_nbr_rstfile : number of different restart file
!@
!@ -- ig_no_rstfile : restart file number
!@
!@ -- cg_name_rstfile : name of restart file corresponding to ig_no_rstfile 
!@
!@ -- cficout : file name for output field in cases 'PIPE' or 'NONE'
!@
!@ -- cg_input_file : file name for input file in cases 'IGNOUT' or 'EXPOUT'
!@
!@ -- ig_lag : lag of exported field
!@
!@ -- nluninp : logical unit for input field
!@
!@ -- nlunout : logical unit for output field
!@
!@ -- cstate : field I/O status of indirect fields
!@
!@ -- ig_total_state : field I/O status of all the fields
!@
!@ -- lg_state : LOGICAL indicating if the field goes through Oasis (.true.)
!@               or not (.false.)
!@
!@ -- ifrqmin : minimum exchange frequency of indirect fields
!@
!@ -- ig_total_frqmin :  minimum exchange frequency of all the fields
!@
!@ -- ig_number_field : Rank of the field going through oasis in total number 
!@    of fields 
!@
!@ -- ig_local_trans : Local transformation of fields 
!@
!@ Second line of SSCS :
!@ -------------------
!@
!@ -- nlonbf : number of longitudes for initial fields
!@
!@ -- nlatbf : number of latitudes for initial fields
!@
!@ -- nlonaf : number of longitudes for initial fields
!@
!@ -- nlataf : number of latitudes for initial fields
!@
!@ -- cficbf : root name of specific files for initial fields
!@
!@ -- cficaf : root name of specific files for final fields
!@
!@ -- ig_grid_nbrbf : grid number associated to a field before interpolation
!@ 
!@ -- ig_grid_nbraf : grid number associated to a field after interpolation
!@
!@ -- nseqn  : sequential index of indirect fields
!@
!@ -- nlagn  : lag index of indirect fields
!@
!@ -- ig_total_nseqn : sequential index of all the fields
!@
!@ -- ig_invert: index=1 if INVERT is used for that field
!@
!@ -- ig_reverse: index=1 if REVERSE is used for that field
!@   
!@ Third line of SSCS :
!@ ------------------
!@
!@ -- canal : names of analysis performed for each field 
!@
!@ Relation between symbolic name of the field and port number of 
!@ CLIM library :
!@  
!@ -- ig_portin_id : port number of input fields
!@
!@ -- ig_portout_id : port number of output fields
!@
!  -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: numlab, ig_numlab, nfexch, ig_ntrans
  INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: ig_total_ntrans
  INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: nluinp, nluout, nlonbf
  INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: nlatbf, nlonaf, nlataf
  INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: nseqn, ig_total_nseqn
  INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: ig_freq, ig_lag, nlagn 
  INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: ig_invert, ig_reverse
  INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: ig_number_field, ig_no_rstfile
  INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: ig_total_state, ig_local_trans
  INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: ig_portin_id, ig_portout_id
  INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: ig_grid_nbrbf, ig_grid_nbraf
  INTEGER (kind=ip_intwp_p)                          :: ifrqmin, ig_nbr_rstfile
  INTEGER (kind=ip_intwp_p)                          :: ig_total_frqmin
  
  
  LOGICAL,DIMENSION(:),ALLOCATABLE :: lg_state
!
  CHARACTER(len=8),DIMENSION(:),ALLOCATABLE :: cnaminp, cnamout
  CHARACTER(len=8),DIMENSION(:,:),ALLOCATABLE :: canal 
  CHARACTER(len=8),DIMENSION(:),ALLOCATABLE :: cg_name_rstfile, cg_restart_file
  CHARACTER(len=8),DIMENSION(:),ALLOCATABLE :: cficinp, cficout
  CHARACTER(len=32),DIMENSION(:),ALLOCATABLE :: cg_input_file
  CHARACTER(len=8),DIMENSION(:),ALLOCATABLE :: cg_input_field, cg_output_field
  CHARACTER(len=8),DIMENSION(:),ALLOCATABLE :: cficbf, cficaf, cstate
  CHARACTER(len=4),DIMENSION(:),ALLOCATABLE :: cga_locatorbf, cga_locatoraf

!     -------------------------------------------------------------------
END MODULE mod_string
