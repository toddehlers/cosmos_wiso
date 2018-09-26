SUBROUTINE readwisosw_d

! Routine to read in isotope values of ocean surface waters
! M. Werner, AWI, 2012


  USE mo_kind,          ONLY: dp
  USE mo_io,            ONLY: ini_spec, ini_surf, io_file_id           &
                            , io_read, io_var_id, io_get_var_double    &
                            , io_open_unit, io_inq_varid, io_close	   &
                            , ini_wisosw

  USE mo_control,       ONLY: ldebugio
  USE mo_mpi,           ONLY: p_io, p_pe
  USE mo_doctor,        ONLY: nerr, nout
  USE mo_memory_base,   ONLY: get_stream_element, memory_info          &
                            , get_stream_element_info
  USE mo_decomposition, ONLY: global_decomposition
  USE mo_transpose,     ONLY: scatter_gp
  USE mo_filename,      ONLY: NETCDF  
  USE mo_wiso,          ONLY: nwiso, nisw
  USE mo_memory_wiso,   ONLY: wiso
  
  IMPLICIT NONE

  TYPE (memory_info) :: info

  REAL(dp), POINTER  :: zin(:,:,:), zin2(:,:,:), zptr(:,:,:)

  INTEGER  :: jt


! read in water isotope ocean surface fields
  IF (nwiso > 0) THEN
     IF (p_pe == p_io) THEN

        WRITE(nout,'(/,A,I2)') ' Read surf seawater wiso climatology from unit: ', nisw

        ini_wisosw%format = NETCDF
        CALL IO_open_unit(nisw, ini_wisosw, IO_READ)
        IO_file_id = ini_wisosw%file_id

     ENDIF

     IF (p_pe == p_io) THEN

        IF (ldebugio) WRITE(nerr,*) 'IO_initial : read ', 'wisosw_d'

        CALL get_stream_element_info (wiso, 'wisosw_d', info)

        ALLOCATE (zin(info%gdim(1), info%gdim(3), info%gdim(2)))
        ALLOCATE (zin2(info%gdim(1), info%gdim(2), info%gdim(3)))

        CALL IO_INQ_VARID (IO_file_id, 'wisosw_d', IO_var_id)

        CALL IO_GET_VAR_DOUBLE (IO_file_id, IO_var_id, zin(:,:,:))

        DO jt = 1, nwiso
           zin2(:,jt,:) = zin(1:info%gdim(1),1:info%gdim(3),jt) 
        END DO

     ENDIF
       
     CALL get_stream_element (wiso, 'wisosw_d', zptr)

     CALL scatter_gp (zin2, zptr, global_decomposition)

     IF (p_pe == p_io)  DEALLOCATE (zin)     
     IF (p_pe == p_io)  DEALLOCATE (zin2)     
     
     IF (p_pe == p_io) CALL IO_close (ini_wisosw)
        
  ENDIF ! nwiso>0

END SUBROUTINE readwisosw_d
