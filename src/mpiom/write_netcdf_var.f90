      SUBROUTINE WRITE_NETCDF_VAR(ncid,desc,arr,klev,time)

!**************************************************************************
!
! Gathers a global variable from all PEs and writes it to a NETCDF file
!
! The NETCDF File is only accessed by p_io
!
!**************************************************************************

      USE mo_param1
      USE mo_parallel
      
      implicit none

      INTEGER ncid, klev, time, ndims
      CHARACTER (LEN=*) desc
      REAL arr(ie,je,klev)

      REAL arr_g(ie_g,je_g,klev)
      CHARACTER (LEN=80) err_text

      INCLUDE 'netcdf.inc'
      INTEGER ncstat,ncvarid,k
      INTEGER, ALLOCATABLE	:: start(:),count(:)

! Gather variable

      DO k=1,klev
         CALL gather_arr(arr(:,:,k),arr_g(:,:,k),p_io)
      ENDDO

      if (klev==1.and.time==0) then
        ndims=2
      else if (klev==1.or.time==0) then
        ndims=3
      else
        ndims=4
      endif
      ALLOCATE(start(ndims))
      ALLOCATE(count(ndims))

      start(1)=1; count(1)=ie_g
      start(2)=1; count(2)=je_g
      if (klev.gt.1.or.time.gt.0) then
        if (klev.gt.1.and.time.gt.0) then
          start(3)=1; count(3)=klev
          start(4)=time; count(4)=1
        else if (klev.gt.1.and.time==0) then
          start(3)=1; count(3)=klev
        else
          start(3)=time; count(3)=1
        endif
      endif

! Write NETCDF data

      IF(p_pe==p_io) THEN
        err_text = 'NETDCF: Problem with varid for '//desc
        ncstat = NF_INQ_VARID(ncid,desc,ncvarid)
        IF ( ncstat .NE. NF_NOERR ) CALL stop_all(err_text)
        err_text = 'NETDCF: Problem with put for '//desc
        ncstat = NF_PUT_VARA_DOUBLE(ncid,ncvarid,start,count,arr_g)
        IF ( ncstat .NE. NF_NOERR ) CALL stop_all(err_text)
        err_text = 'NETDCF: Problem with flush for '//desc
        ncstat = NF_SYNC(ncid)
        IF ( ncstat .NE. NF_NOERR ) CALL stop_all(err_text)
      ENDIF

      END
