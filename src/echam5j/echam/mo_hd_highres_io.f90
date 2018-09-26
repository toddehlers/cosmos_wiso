MODULE mo_hd_highres_io

  USE mo_kind,         ONLY: dp
  use mo_exception,    ONLY: message_text, message, finish
  use mo_time_control, ONLY: get_date_components, next_date 
  use mo_filename,     ONLY: compose_filenames, standard_grib_file

  IMPLICIT NONE

  PRIVATE
  
  INCLUDE 'netcdf.inc'

  ! netCDF id
  INTEGER ::  ncid

  ! dimension ids

  INTEGER ::  lon_dim
  INTEGER ::  lat_dim
  INTEGER ::  time_dim

  ! dimension lengths

  INTEGER, PARAMETER ::  lon_len = 720
  INTEGER, PARAMETER ::  lat_len = 360
  INTEGER, PARAMETER ::  time_len = NF_UNLIMITED

  ! variable ids

  INTEGER ::  lon_id
  INTEGER ::  lat_id
  INTEGER ::  time_id

  INTEGER ::  friv_id

  ! variable shapes

  INTEGER :: dim1(1)
  INTEGER :: dim3(3)

  ! grid basic parameter

  REAL(dp), PARAMETER :: florg = -180.0_dp
  REAL(dp), PARAMETER :: fborg =   90.0_dp
  REAL(dp), PARAMETER :: fscal =    0.5_dp

  INTEGER, SAVE :: icount = 0

  PUBLIC :: hd_highres_open
  PUBLIC :: hd_highres_close
  PUBLIC :: hd_highres_write

CONTAINS

  SUBROUTINE hd_highres_open

    ! data variables

    REAL(dp) :: lon(lon_len)
    REAL(dp) :: lat(lat_len)

    INTEGER :: jb, jl

    CHARACTER(len=256) :: hd_highres_filename

    CALL compose_filenames
    hd_highres_filename = TRIM(standard_grib_file)//'_hd_higres.nc'

    write (0,*) 'HD model high resolution river discharge output: ', TRIM(hd_highres_filename)

    ! enter define mode

    CALL check_nf(nf_create(TRIM(hd_highres_filename), NF_CLOBBER, ncid))

    WRITE (message_text,*)                                    &
         'HD model high resolution river discharge output: ', &
         TRIM(hd_highres_filename)
    CALL message('hd_highres_open', message_text)


    ! define dimensions

    CALL check_nf(nf_def_dim(ncid, 'lon', lon_len, lon_dim))
    CALL check_nf(nf_def_dim(ncid, 'lat', lat_len, lat_dim))
    CALL check_nf(nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim))

    ! define variables

    dim1(1) = lon_dim
    CALL check_nf(nf_def_var(ncid, 'lon', NF_DOUBLE, 1, dim1, lon_id))
    dim1(1) = lat_dim
    CALL check_nf(nf_def_var(ncid, 'lat', NF_DOUBLE, 1, dim1, lat_id))
    dim1(1) = time_dim
    CALL check_nf(nf_def_var(ncid, 'time', NF_DOUBLE, 1, dim1, time_id))

    dim3(1:3) = (/ lon_dim, lat_dim, time_dim /)
    CALL check_nf(nf_def_var(ncid, 'friv', NF_REAL, 3, dim3, friv_id))

    ! assign attributes
    CALL check_nf(nf_put_att_text(ncid, lon_id, 'long_name', 9, 'longitude'))
    CALL check_nf(nf_put_att_text(ncid, lon_id, 'units', 12, 'degrees_east'))
    CALL check_nf(nf_put_att_text(ncid, lon_id, 'standard_name', 9, 'longitude'))

    CALL check_nf(nf_put_att_text(ncid, lat_id, 'long_name', 8, 'latitude'))
    CALL check_nf(nf_put_att_text(ncid, lat_id, 'units', 13, 'degrees_north'))
    CALL check_nf(nf_put_att_text(ncid, lat_id, 'standard_name', 8, 'latitude'))

    CALL check_nf(nf_put_att_text(ncid, time_id, 'units', 16, 'day as %Y%m%d.%f'))
    CALL check_nf(nf_put_att_text(ncid, time_id, 'calendar',19, 'proleptic_gregorian'))

    CALL check_nf(nf_put_att_text(ncid, friv_id, 'long_name', 15, 'river discharge'))
    CALL check_nf(nf_put_att_text(ncid, friv_id, 'standard_name', 45, 'water_volume_transport_into_ocean_from_rivers'))
    CALL check_nf(nf_put_att_text(ncid, friv_id, 'units', 6, 'm3 s-1'))

    CALL check_nf(nf_put_att_text(ncid, NF_GLOBAL, 'conventions', 6, 'CF-1.0'))
    CALL check_nf(nf_put_att_text(ncid, NF_GLOBAL, 'source', 15, 'ECHAM5/HD model'))
    CALL check_nf(nf_put_att_text(ncid, NF_GLOBAL, 'institution', 36, 'Max-Planck-Institute for Meteorology'))

    ! leave define mode
    CALL check_nf(nf_enddef(ncid))

    DO jb = 1, lat_len
      lat(jb) = fborg-REAL(jb,dp)*fscal+0.5_dp*fscal
    ENDDO
    DO jl = 1, lon_len
      lon(jl) = REAL(jl,dp)*fscal+florg-0.5_dp*fscal
      IF (lon(jl) >= 180.0_dp) lon(jl) = lon(jl)-360.0_dp
    ENDDO

    CALL check_nf(nf_put_var_double(ncid, lat_id, lat))
    CALL check_nf(nf_put_var_double(ncid, lon_id, lon))

  END SUBROUTINE hd_highres_open

  SUBROUTINE hd_highres_close

    CALL check_nf(nf_close(ncid))

  END SUBROUTINE hd_highres_close
    
  SUBROUTINE hd_highres_write(friv)
       
    REAL(dp), INTENT(in) :: friv(:,:) 

    ! starts and counts for array sections of record variables

    INTEGER ::  tstart(1), tcount(1)
    INTEGER ::  fstart(3), fcount(3)

    INTEGER :: year, month, day, hour, minute, second
    REAL(dp) :: yyyymmdd

    icount = icount+1

    CALL get_date_components(next_date, year, month, day, hour, minute, second)
    yyyymmdd = ABS(year)*10000+month*100+day              &
             + (hour*3600+minute*60+second)/86400._dp
    IF (year < 0) yyyymmdd = -yyyymmdd

    tstart(1) = icount
    tcount(1) = 1
    CALL check_nf(nf_put_vara_double(ncid, time_id, tstart, tcount, yyyymmdd))

    fstart(1:3) = (/ 1, 1, icount /)
    fcount(1:3) = (/ lon_len, lat_len, 1 /)
    CALL check_nf(nf_put_vara_double(ncid, friv_id, fstart, fcount, friv))

  END SUBROUTINE hd_highres_write
  
  SUBROUTINE check_nf(iret)

    INTEGER, INTENT(in) :: iret
    
    IF (iret /= NF_NOERR) THEN
      write (message_text,'(a)') nf_strerror(iret)
      call finish ('mo_hd_highres_io', TRIM(message_text))
    ENDIF

  END SUBROUTINE check_nf

END MODULE mo_hd_highres_io
