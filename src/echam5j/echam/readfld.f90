SUBROUTINE readfld

  ! Description:
  !
  ! *readfld* reads optional global fields
  !           to be used in various subroutines
  !
  ! Method:
  !
  ! Memory is allocated from the heap and pointers
  ! for the start addresses are stored in modules.
  ! (See *xtemiss* for the use of these pointers)
  !
  ! *readfld* is called from *control*
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, January 1995, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, Oct 1999, netCDF version
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
  USE mo_field,         ONLY: field1, field2
  USE mo_control,       ONLY: ngl, nlon, numfl1, numfl2,    &
                              nfl1, nfl2
  USE mo_time_control,  ONLY: lstart
  USE mo_decomposition, ONLY: dc    => local_decomposition, &
                              gl_dc => global_decomposition
  USE mo_transpose,     ONLY: scatter_gp
  USE mo_test_trans,    ONLY: test_gridpoint

  IMPLICIT NONE

  REAL(dp), POINTER :: ftmp (:,:,:) ! global read buffer

  !  Executable statements 

!-- 1. Allocate memory and read fields

  IF (numfl1==0 .AND. numfl2==0) RETURN

!-- 1.1   Allocate memory for global fields

  ! Fields to be read at nstep=0 only

  IF (numfl1>0 .AND. lstart) THEN
    ALLOCATE (field1(dc%nglon, numfl1, dc%nglat))
  END IF

  ! Fields to be read at nstep=nresum (including nstep=0)

  IF (numfl2>0) THEN
    ALLOCATE (field2(dc%nglon, numfl2, dc%nglat))
  END IF

  ! 1.2 Open file(s) and read

!-- Fields to be read at nstep=0 only

  IF (numfl1>0) THEN

    CALL readfield(nfl1, numfl1, ftmp)

    CALL scatter_gp(ftmp, field1, gl_dc)
    IF (ASSOCIATED(ftmp)) DEALLOCATE (ftmp)
    CALL test_gridpoint (field1, 'readfld: field1')

  END IF

  ! Fields to be read at nstep=nresum (including nstep=0)

  IF (numfl2>0) THEN

    CALL readfield(nfl2, numfl2, ftmp)

    CALL scatter_gp(ftmp, field2, gl_dc)
    IF (ASSOCIATED(ftmp)) DEALLOCATE (ftmp)
    CALL test_gridpoint (field2, 'readfld: field2')

  END IF

  RETURN

CONTAINS

  SUBROUTINE readfield(nfl, numfl, ftmp)

    USE mo_kind,          ONLY: dp
    USE mo_exception,     ONLY: finish
    USE mo_mpi,           ONLY: p_pe, p_io
    USE mo_doctor,        ONLY: nout, nerr
    USE mo_io,            ONLY: IO_inq_dimid, IO_inq_dimlen, IO_file_id, &
                                io_var_id, IO_inq_varid, IO_get_vara_double, &
                                ini_field, IO_READ, IO_open_unit, IO_close

    REAL(dp), POINTER :: ftmp (:,:,:) ! global read buffer

    !  Local scalars: 
    INTEGER       :: start(3), COUNT(3)
    INTEGER       :: io_ngl, io_nlon, io_nfld
    INTEGER       :: js
    INTEGER       :: nfl, numfl

    NULLIFY(ftmp)
    IF(p_pe == p_io) THEN

      WRITE (nout,*)
      IF (numfl == 1) THEN
        WRITE (nout,*) 'Reading ', numfl, ' field from', ' unit ', nfl
      ELSE
        WRITE (nout,*) 'Reading ', numfl, ' fields from', ' unit ', nfl
      END IF

      CALL IO_open_unit (nfl, ini_field, IO_READ)
      IO_file_id = ini_field%file_id

      ! Check resolution
      CALL IO_inq_dimid  (IO_file_id, 'lat', io_var_id)
      CALL IO_inq_dimlen (IO_file_id, io_var_id, io_ngl)
      CALL IO_inq_dimid  (IO_file_id, 'lon', io_var_id)
      CALL IO_inq_dimlen (IO_file_id, io_var_id, io_nlon)
      CALL IO_inq_dimid  (IO_file_id, 'nfield', io_var_id)
      CALL IO_inq_dimlen (IO_file_id, io_var_id, io_nfld)

      IF (io_nlon /= nlon .OR. io_ngl /= ngl) THEN
         WRITE(nerr,*) 'readfield: unexpected resolution ', io_nlon, io_ngl
         CALL finish ('readfield', 'unexpected resolution')
      END IF

      IF (numfl > io_nfld) THEN
         WRITE(nerr,*) 'readfield: only ', io_nfld, ' field(s) available '
         CALL finish ('readfield', 'too much fields selected')
      END IF

      ! Read field
      ALLOCATE (ftmp(nlon,numfl,ngl)); ftmp(:,:,:) = 0._dp

      CALL IO_inq_varid (IO_file_id, 'FIELD', io_var_id)

      DO js = 1, numfl
        COUNT(:) = (/ nlon, ngl, 1 /)
        start(:) = (/ 1, 1, js /)
        CALL IO_get_vara_double (IO_file_id, io_var_id, start, count, ftmp(:,js,:))
      END DO

      CALL IO_close (ini_field)
    ENDIF

  END SUBROUTINE readfield

END SUBROUTINE readfld
