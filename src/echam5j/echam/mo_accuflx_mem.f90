MODULE mo_accuflx_mem
  !
  ! Authors:
  ! -------
  ! Philip Stier, MPI-MET                                                       2001/2002

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream , HYBRID, SURFACE, BELOWSUR, &
                              ABOVESUR10, NETCDF, GRIB, HYBRID_H
  USE mo_memory_base,   ONLY: new_stream, add_stream_element,       &
                              default_stream_setting, add_stream_reference, AUTO
  USE mo_netCDF,        ONLY: max_dim_name
  USE mo_parameters,    ONLY: jpgrnd
  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_control,       ONLY: nlon, nlev, ngl  ! size of global arrays
  USE mo_time_event,    ONLY: io_time_event
  USE mo_memory_g3b,    ONLY: netht_sw, netht_lw, emter1, trsol1,  &
                              emtef01, trsof01, emtef1, trsof1

  IMPLICIT NONE

  PRIVATE

  !--- Service routines: ----------------------------------------------------------------

  PUBLIC :: construct_stream_accuflx            ! construct the accuflx stream

  !--- Declarations for stream accuflx: ----------------------------------------------------

  TYPE (t_stream), PUBLIC, POINTER :: accuflx
  
  !TYPE(io_time_event)      :: input_out=io_time_event(1,'months','last',0) ! Output interval 
  TYPE(io_time_event),SAVE  :: input_out=io_time_event(1,'months','last',0) ! Output interval 

  REAL(dp), PUBLIC, POINTER :: d_aflx_sw(:,:,:)
  REAL(dp), PUBLIC, POINTER :: d_aflx_lw(:,:,:)
  REAL(dp), PUBLIC, POINTER :: d_aflx_swc(:,:,:)
  REAL(dp), PUBLIC, POINTER :: d_aflx_lwc(:,:,:)
  
  !--- Declare dimensions for the streams: ----------------------------------------------------------------------

  INTEGER :: lnlon,   lnlev,  lngl ! size of local arrays on PE
  INTEGER :: nlevp1
  INTEGER :: lnlevp1
  INTEGER :: dim1(2), dim1p(2)
  INTEGER :: dim2(2), dim2p(2)
  INTEGER :: dim3(3), dim3p(3)
  CHARACTER (max_dim_name) :: dim1n(2), dim2n(2), dim3n(3)
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE construct_stream_accuflx

    !--- Set dimensions to be used in the declarations:

    lnlon=lc%nproma
    lnlev=lc%nlev
    lngl =lc%ngpblks

    nlevp1  = nlev  + 1
    lnlevp1 = lnlev + 1

    dim1p = (/ lnlon, lnlev  /)
    dim1  = (/  nlon,  nlev  /)
    dim1n = (/  "lon","lev" /)

    dim2p = (/ lnlon, lngl  /)
    dim2  = (/  nlon,  ngl  /)
    dim2n = (/  "lon","lat" /)
    
    dim3p = (/ lnlon,  lnlevp1, lngl  /)
    dim3  = (/  nlon,   nlevp1,  ngl  /)
    dim3n = (/  "lon ","ilev","lat "/)
    
    
    !--- 1) Construct the accuflx stream: ------------------------------------------------------------------------------

    CALL new_stream (accuflx,'accuflx',filetype=GRIB,lrerun=.true.)

    !--- Add standard fields for post-processing:

    CALL add_stream_reference (accuflx, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (accuflx, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (accuflx, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (accuflx, 'gboxarea','geoloc',lpost=.TRUE.)

    !--- 2) Add stream elements: ------------------------------------------------------------------------------------

    CALL default_stream_setting (accuflx, units = 'Wm**-2',    &
                                      lrerun   = .TRUE. ,     &
                                      contnorest=.TRUE. ,    &
                                      laccu    = .TRUE. ,     &
                                      lpost    = .TRUE. ,     &
                                      gdims    = dim3 ,       &
                                      dimnames = dim3n ,      &
                                      leveltype = HYBRID_H,   &
                                      table    = 199,         &
                                      code     = AUTO         )

    !--- 2.1) flux anomalies, perturbed minus unperturbed
    
    CALL add_stream_element (accuflx, 'd_aflx_sw', d_aflx_sw, dim3p, code=73,  &
                             longname='Accumulated SW flux anomalies - all sky ')
    CALL add_stream_element (accuflx, 'd_aflx_lw', d_aflx_lw, dim3p, code=74,  &
                             longname='Accumulated LW flux anomalies - all sky')
    CALL add_stream_element (accuflx, 'd_aflx_swc', d_aflx_swc, dim3p, code=75,  &
                             longname='Accumulated SW flux anomalies - clear')
    CALL add_stream_element (accuflx, 'd_aflx_lwc', d_aflx_lwc, dim3p, code=76,  &
                             longname='Accumulated LW flux anomalies - clear')
    !
    !--- 2.2) heating rates
    !
    CALL add_stream_element (accuflx, 'netht_sw', netht_sw, code=77,  &
                             longname='Net SW heating rates', units='K/s')
    CALL add_stream_element (accuflx, 'netht_lw', netht_lw, code=78,  &
                             longname='Net LW heating rates', units='K/s')

    !
    !--- 2.3) flux variables of second call of radiation not written per default
    !
    CALL add_stream_element (accuflx,'emter1',  emter1,  lpost=.FALSE.)
    CALL add_stream_element (accuflx,'trsol1',  trsol1,  lpost=.FALSE.)
    CALL add_stream_element (accuflx,'emtef01', emtef01, lpost=.FALSE.)
    CALL add_stream_element (accuflx,'trsof01', trsof01, lpost=.FALSE.)
    CALL add_stream_element (accuflx,'emtef1',  emtef1,  lpost=.FALSE.)
    CALL add_stream_element (accuflx,'trsof1',  trsof1,  lpost=.FALSE.)
    
  END SUBROUTINE construct_stream_accuflx
  

!---------------------------------------------------------------------------------------
END MODULE mo_accuflx_mem
