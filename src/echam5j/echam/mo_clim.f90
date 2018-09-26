MODULE mo_clim

  ! U. Schulzweida, MPI, May 2002, blocking (nproma)

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  REAL(dp), ALLOCATABLE :: vltclim(:,:,:)    ! (nlon,ngl,0:13)
  REAL(dp), ALLOCATABLE :: vgratclim(:,:,:)  ! (nlon,ngl,0:13)
  REAL(dp), ALLOCATABLE :: tslclim(:,:,:)    ! (nlon,ngl,0:13)

CONTAINS

  SUBROUTINE readtslclim

    USE mo_doctor,        ONLY: nout
    USE mo_control,       ONLY: ngl, nlon, ntslcl
    USE mo_exception,     ONLY: finish
    USE mo_mpi,           ONLY: p_pe, p_io   
    USE mo_decomposition, ONLY: lc => local_decomposition, global_decomposition
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_io

    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_tslclim(:,:,:)

    INTEGER    :: i
    INTEGER    :: io_ngl, io_nlon


    ! Allocate memory for tslclim per PE

    IF (.NOT. ALLOCATED(tslclim)) ALLOCATE (tslclim(lc%nproma, lc%ngpblks, 0:13))

    IF (p_pe == p_io) THEN

       WRITE (nout, '(/)')
       WRITE (nout,*) 'Reading tslclim from', ' unit ', ntslcl

       ! Open file

       yearnc%format = NETCDF
       CALL IO_open_unit (ntslcl, yearnc, IO_READ)
       IO_file_id = yearnc%file_id

       ! Check resolution

       CALL IO_inq_dimid  (IO_file_id, 'lat', io_var_id)
       CALL IO_inq_dimlen (IO_file_id, io_var_id, io_ngl)
       CALL IO_inq_dimid  (IO_file_id, 'lon', io_var_id)
       CALL IO_inq_dimlen (IO_file_id, io_var_id, io_nlon)

       IF (io_nlon /= nlon .OR. io_ngl /= ngl) THEN
          WRITE(nerr,*) 'readtslclim: unexpected resolution ', io_nlon, io_ngl
          CALL finish ('readtslclim', 'unexpected resolution')
       END IF

       ! Allocate memory for tslclim global field

       ALLOCATE (zin(lc%nlon, lc%nlat, 0:13))

       ! Read data

       CALL IO_inq_varid (IO_file_id, 'TSLCLIM', io_var_id)
       CALL IO_get_var_double (IO_file_id, io_var_id, zin(:,:,1:12))

       ! Close file

       CALL IO_close(yearnc)

       zin(:,:,0)  = zin(:,:,12)
       zin(:,:,13) = zin(:,:,1)

       zin(:,:,:)  = MAX(zin(:,:,:), 0.001_dp)

    END IF

    NULLIFY (gl_tslclim)
    DO i = 0, 13
       IF (p_pe == p_io) gl_tslclim => zin(:,:,i:i)
       CALL scatter_gp (gl_tslclim, tslclim(:,:,i:i), global_decomposition)
    END DO

    IF (p_pe == p_io) THEN
       DEALLOCATE (zin)
    END IF

    RETURN
  END SUBROUTINE readtslclim

  SUBROUTINE readvltclim

    USE mo_doctor,        ONLY: nout
    USE mo_control,       ONLY: ngl, nlon, nvltcl
    USE mo_exception,     ONLY: finish
    USE mo_mpi,           ONLY: p_pe, p_io   
    USE mo_decomposition, ONLY: lc => local_decomposition, global_decomposition
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_io

    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_vltclim(:,:,:)

    INTEGER    :: i
    INTEGER    :: io_ngl, io_nlon


    ! Allocate memory for vltclim per PE

    IF (.NOT. ALLOCATED(vltclim)) ALLOCATE (vltclim(lc%nproma, lc%ngpblks, 0:13))

    IF (p_pe == p_io) THEN

       WRITE (nout, '(/)')
       WRITE (nout,*) 'Reading vltclim from', ' unit ', nvltcl

       ! Open file

       CALL IO_open_unit (nvltcl, yearnc, IO_READ)
       IO_file_id = yearnc%file_id

       ! Check resolution

       CALL IO_inq_dimid  (IO_file_id, 'lat', io_var_id)
       CALL IO_inq_dimlen (IO_file_id, io_var_id, io_ngl)
       CALL IO_inq_dimid  (IO_file_id, 'lon', io_var_id)
       CALL IO_inq_dimlen (IO_file_id, io_var_id, io_nlon)

       IF (io_nlon /= nlon .OR. io_ngl /= ngl) THEN
          WRITE(nerr,*) 'readvltclim: unexpected resolution ', io_nlon, io_ngl
          CALL finish ('readvltclim', 'unexpected resolution')
       END IF

       ! Allocate memory for vltclim global field

       ALLOCATE (zin(lc%nlon, lc%nlat, 0:13))

       ! Read data

       CALL IO_inq_varid (IO_file_id, 'VLTCLIM', io_var_id)
       CALL IO_get_var_double (IO_file_id, io_var_id, zin(:,:,1:12))

       ! Close file

       CALL IO_close(yearnc)

       zin(:,:,0)  = zin(:,:,12)
       zin(:,:,13) = zin(:,:,1)

       zin(:,:,:)  = MAX(zin(:,:,:), 0.001_dp)

    END IF

    NULLIFY (gl_vltclim)
    DO i = 0, 13
       IF (p_pe == p_io) gl_vltclim => zin(:,:,i:i)
       CALL scatter_gp (gl_vltclim, vltclim(:,:,i:i), global_decomposition)
    END DO

    IF (p_pe == p_io) THEN
       DEALLOCATE (zin)
    END IF

    RETURN
  END SUBROUTINE readvltclim

  SUBROUTINE readvgratclim

    USE mo_doctor,        ONLY: nout
    USE mo_control,       ONLY: ngl, nlon, nvgratcl
    USE mo_exception,     ONLY: finish
    USE mo_mpi,           ONLY: p_pe, p_io   
    USE mo_decomposition, ONLY: lc => local_decomposition, global_decomposition
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_io

    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_vgratclim(:,:,:)

    INTEGER    :: i
    INTEGER    :: io_ngl, io_nlon


    ! Allocate memory for vgratclim per PE

    IF (.NOT. ALLOCATED(vgratclim)) ALLOCATE (vgratclim(lc%nproma, lc%ngpblks, 0:13))

    IF (p_pe == p_io) THEN

       WRITE (nout, '(/)')
       WRITE (nout,*) 'Reading vgratclim from', ' unit ', nvgratcl

       ! Open file

       CALL IO_open_unit (nvgratcl, yearnc, IO_READ)
       IO_file_id = yearnc%file_id

       ! Check resolution

       CALL IO_inq_dimid  (IO_file_id, 'lat', io_var_id)
       CALL IO_inq_dimlen (IO_file_id, io_var_id, io_ngl)
       CALL IO_inq_dimid  (IO_file_id, 'lon', io_var_id)
       CALL IO_inq_dimlen (IO_file_id, io_var_id, io_nlon)

       IF (io_nlon /= nlon .OR. io_ngl /= ngl) THEN
          WRITE(nerr,*) 'readvgratclim: unexpected resolution ', io_nlon, io_ngl
          CALL finish ('readvgratclim', 'unexpected resolution')
       END IF

       ! Allocate memory for vgratclim global field

       ALLOCATE (zin(lc%nlon, lc%nlat, 0:13))

       ! Read data

       CALL IO_inq_varid (IO_file_id, 'VGRATCLIM', io_var_id)
       CALL IO_get_var_double (IO_file_id, io_var_id, zin(:,:,1:12))

       ! Close file

       CALL IO_close(yearnc)

       zin(:,:,0)  = zin(:,:,12)
       zin(:,:,13) = zin(:,:,1)

    END IF

    NULLIFY (gl_vgratclim)
    DO i = 0, 13
       IF (p_pe == p_io) gl_vgratclim => zin(:,:,i:i)
       CALL scatter_gp (gl_vgratclim, vgratclim(:,:,i:i), global_decomposition)
    END DO

    IF (p_pe == p_io) THEN
       DEALLOCATE (zin)
    END IF

    RETURN
  END SUBROUTINE readvgratclim
!------------------------------------------------------------------------------
  subroutine cleanup_clim
    !
    ! cleanup module variables
    !
    if (allocated (vltclim))   deallocate (vltclim)
    if (allocated (vgratclim)) deallocate (vgratclim)
    if (allocated (tslclim))   deallocate (tslclim)
  end subroutine cleanup_clim

END MODULE mo_clim
