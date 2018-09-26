      SUBROUTINE SAVE_TIMESER_ADD

!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/save_timeser_add.f90,v $\\

!
!****************************************************************
!
!**** *SAVE_TIMSER_ADD*
!

!
!     Purpose
!     -------
!     - save add time series.
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *SAVE_TIMESER_ADD*
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *PARAM1_ADD.h* - declaration of ocean/sediment tracer.
!     *COMMON*     *COMMO1_ADD.h* - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS_ADD.h*  - std I/O logical units.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd dimension) [m].
!
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************
!ik added trap depths for time series
!ik changed output to ASCII formatted file
!ik added budget output to be written at the end of addout (being too lazy
!ik to write a seperate routine)

      USE mo_timeser_add
      USE mo_control_add
      USE mo_contra
      use mo_param1_add 

      use mo_mpi 
      use mo_parallel     

implicit none
      
      INTEGER :: i,n,l,l1,l2,l3
      REAL :: ts1_time(lents1)
      REAL :: ts1recv(nvarts1_add,nelets1,lents1)
      CHARACTER*33 ystring

#ifdef PNETCDF
      INCLUDE 'netcdf.inc'

      INTEGER :: ncid,ncstat,ncvarid
      INTEGER :: nvarts1id, nelets1id, nctimeid
      INTEGER :: ncdims(3)

      WRITE(io_stdo_add,*)'Save Timeser ADD in NetCDF Format'
#endif

      
! For a parallel run, sum up ts1 on p_io
! Note that this code also works for nonparallel runs,
! it just does nothing in this case!

      ! All non-I/O-PEs just send to p_io
      IF(p_pe/=p_io .AND. p_pe<nprocxy) CALL p_send(ts1,p_io,123)

      ! The I/O-PE sums up
      IF(p_pe==p_io) THEN
        DO n=0,nprocxy-1
          IF(n/=p_io) THEN
            CALL p_recv(ts1recv,n,123)
            ts1(:,:,:) = ts1(:,:,:) + ts1recv(:,:,:)
          ENDIF
        ENDDO
      ENDIF

      IF(p_pe==p_io) THEN ! Output is only done by I/O-PE




#ifdef PNETCDF

!
! Calculate time
!

      ystring(1:33)='seconds since 0001-01-01 23:59:59'

      if (addstartyear .ge. 1000) then
        WRITE(ystring(15:18),'(I4)') addstartyear
      elseif (addstartyear .ge. 100) then
        WRITE(ystring(16:18),'(I3)') addstartyear   
      elseif (addstartyear .ge. 10) then
        WRITE(ystring(17:18),'(I2)') addstartyear     
      else
        WRITE(ystring(18:18),'(I1)') addstartyear    
      endif
      
      if (addstartmonth .ge. 10) then
        WRITE(ystring(20:21),'(I2)') addstartmonth
      else
        WRITE(ystring(21:21),'(I1)') addstartmonth
      endif

      if (addstartday .ge. 10) then
        WRITE(ystring(23:24),'(I2)') addstartday
      else
        WRITE(ystring(24:24),'(I1)') addstartday
      endif
      
      WRITE(io_stdo_add,*) 'Write timeser_add in ', ystring   

      do i=1,lents1
         ts1_time(i) = dt_add*nfreqts1/2. + dt_add*nfreqts1*(i-1)
!      WRITE(io_stdo_bgc,*)'TEST:', ts1_time(i)/86400
      enddo
      


!
! Open netCDF data file
!
      ncstat = NF_CREATE('timeser_add.nc',NF_CLOBBER, ncid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('SAVE_TIMESER_ADD: Problem with netCDF1')

      ! Define global attributes
      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'title'             &
     &,21, 'marine add timeseries')
      IF ( ncstat .NE. NF_NOERR ) call stop_all ('SAVE_TIMESER_ADD: Problem with netCDF2')

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'source'             &
     &,7, 'ADDCONTRA')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF3')

!      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'author'             &
!     &,42, 'Patrick Wetzel, <patrick.wetzel@dkrz.de>')
!      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_BGC: Problem with netCDF4')

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'institution'             &
     &,45, 'Max-Planck-Institute for Meteorology & AWI, Hamburg & Bremenhaven')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF5')

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'Conventions'             &
     &,12, 'non-standard')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('write_all_to_ncdf.f90: Problem with netCDF6')

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'History'             &
     &,7, 'created')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('write_all_to_ncdf.f90: Problem with netCDF7')

      WRITE(io_stdo_add,*)'Global attributes of Timeser ADD'
!
! Define dimension
! ---------------------------------------------------------------------------------
      ncstat = NF_DEF_DIM(ncid, 'dim_nvarts1',nvarts1_add,nvarts1id)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF11')

      ncstat = NF_DEF_DIM(ncid, 'dim_nelets1',nelets1,nelets1id)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF12')

!      ncstat = NF_DEF_DIM(ncid, 'dim_lents1',lents1,lents1id)
!      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCD13')


      ncstat = NF_DEF_DIM(ncid, 'time',lents1,nctimeid)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF14')

      WRITE(io_stdo_add,*)'Define dimensions of Timeser ADD'

       
! Define grid variables
! ---------------------------------------------------------------------------------


      ncdims(1) = nctimeid
      
      ncstat = NF_DEF_VAR(ncid,'time',NF_REAL,1,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF21')

      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name',4, 'time')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF22')

      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',33,ystring)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF22a')

      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'calendar',9, 'gregorian')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF22b')


      WRITE(io_stdo_add,*)'Define time variable of Timeser ADD'

!               time:long_name = "time" ;
!               time:units = "day since 1999-10-01 00:00:00" ;
!               time:calendar = "gregorian" ; 

! Define data variables
! ---------------------------------------------------------------------------------

      ncdims(1) = nvarts1id
      ncdims(2) = nelets1id
      ncdims(3) = nctimeid
      
      
      ncstat = NF_DEF_VAR(ncid,'timeseries',NF_REAL,3,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF31')

      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',4, 'none')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF32')

      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                  &
     &,15, 'timeseries data')
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF33')

      ncstat = NF_PUT_ATT_REAL(ncid,ncvarid,'missing_value',NF_FLOAT,1,99999.)
 
      WRITE(io_stdo_add,*)'Define timeseries variable of Timeser ADD'

! Define data variables
! ---------------------------------------------------------------------------------

      ncdims(1) = nelets1id
      ncdims(2) = nctimeid
      
      CALL NETCDF_DEF_VARSG(ncid,5,'h2o18',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**3',5, 'H2O18',rmasko,40,io_stdo_add)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 200)

      CALL NETCDF_DEF_VARSG(ncid,5,'hDo16',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**3',5, 'HDO16',rmasko,42,io_stdo_add)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 202)
      
      CALL NETCDF_DEF_VARSG(ncid,5,'h2o16',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**3',5, 'H2O16',rmasko,41,io_stdo_add)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 201)                

      CALL NETCDF_DEF_VARSG(ncid,9,'h2o18_ice',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**3',9, 'H2O18_ICE',rmasko,43,io_stdo_add)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 203)

      CALL NETCDF_DEF_VARSG(ncid,9,'hDo16_ice',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**3',9, 'HDO16',rmasko,44,io_stdo_add)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 202)
      
      CALL NETCDF_DEF_VARSG(ncid,9,'h2o16_ice',2,ncdims,ncvarid,                      &
     &   9,'kmol/m**3',9, 'H2O16_ICE',rmasko,45,io_stdo_add)
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, 201)                

     
     
! ---------------------------------------------------------------------------------
! END Define variables
      ncstat = NF_ENDDEF(ncid)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF40')

      WRITE(io_stdo_add,*)'end define of Timeser ADD'

!
! Set fill mode
! ---------------------------------------------------------------------------------
      ncstat = NF_SET_FILL(ncid,NF_NOFILL,ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF41')
  
      WRITE(io_stdo_add,*)'fill mode of Timeser BGC'

! Save grid variables
! ---------------------------------------------------------------------------------

      ncstat = NF_INQ_VARID(ncid,'time',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF42')

      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1_time)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF43')

      WRITE(io_stdo_add,*)'fill time of Timeser ADD'
      

! Save timeseries data
!---------------------------------------------------------------------
      ncstat = NF_INQ_VARID(ncid,'timeseries',ncvarid )
 !      WRITE(0,*) 'timeseries', NF_STRERROR(ncstat)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF44')

      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF45')

      WRITE(io_stdo_add,*)'fill timeseries of Timeser ADD'

!
! Write timeser data 
!
!-----------------------------------------------------------------------

      ncstat = NF_INQ_VARID(ncid,'h2o18',ncvarid )
 !     WRITE(0,*) 'h2O18', NF_STRERROR(ncstat)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF44')
      ncstat = NF_INQ_VARID(ncid,'hDo16',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF44')
      ncstat = NF_INQ_VARID(ncid,'h2o16',ncvarid)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF44')
      ncstat = NF_INQ_VARID(ncid,'h2o18_ice',ncvarid)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF44')
      ncstat = NF_INQ_VARID(ncid,'hDo16_ice',ncvarid)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF44')
      ncstat = NF_INQ_VARID(ncid,'h2o16_ice',ncvarid)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF44')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsh2o18,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF45')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itshDo16,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF45')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsh2o16,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF45') 
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsh2o18_ice,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF45')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itshDo16_ice,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF45')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,ts1(itsh2o16_ice,:,:))
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF45') 

!
! Close File

      ncstat = NF_CLOSE(ncid)
      IF ( ncstat .NE. NF_NOERR ) CALL stop_all('SAVE_TIMESER_ADD: Problem with netCDF50')

      WRITE(io_stdo_add,*)'Timeseries ADD written to "timeser_add.nc"'


#else
      IF( lents1 .LT. lts1 ) THEN
         WRITE(io_stdo_add,*)                                          &
     &   'Warning : too many samples in time series 1!!!'
      ELSE
         io_timeser_add = io_stdi_add
         WRITE(io_stdo_add,*)'Opening timeser_add (',io_timeser_add,').'

         CLOSE(io_timeser_add)
         OPEN(io_timeser_add,FILE='timeser_add',STATUS='UNKNOWN'       &
     &       ,FORM='FORMATTED')
         WRITE(io_timeser_add,*) 'time steps, number of variables,     & 
     &                            number of stations '
         WRITE(io_timeser_add,*) lts1,nvarts1_add,nelets1
         WRITE(io_timeser_add,*) (rlonts1(l),l=1,nelets1)
         WRITE(io_timeser_add,*) (rlatts1(l),l=1,nelets1)
 
         DO l2=1,nelets1 !outer vertical loop is station
         DO l3=1,lts1 !inner vertical loop is time
           WRITE(io_timeser_add,433)(ts1(l1,l2,l3),l1=1,nvarts1_add) !horizontal loop is variable type
         ENDDO
         ENDDO

         CLOSE(io_timeser_add)
433      FORMAT(19e15.6) 

      ENDIF
#endif
      ENDIF ! p_pe==p_io

         WRITE(io_stdo_add,*)'Memory deallocation ts1 ...'
         WRITE(io_stdo_add,*)'No. of variables: ',nvarts1_add
         WRITE(io_stdo_add,*)'No. of stations : ',nelets1
         WRITE(io_stdo_add,*)'No. of time steps: ',lts1
         DEALLOCATE (ts1)

      !CALL p_barrier ! Not really necessary, just to wait for the I/O-PE

      RETURN
      END
