      SUBROUTINE AUFW_ADD(kpie,kpje,kpke,pddpo,pgila,pgiph,ptiestu   &
     &                    ,kplyear,kplmon,kplday,kpldtoce)

!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/aufw_bgc.f90,v $\\
!$Revision: 1.3.2.1.16.1.4.2 $\\
!$Date: 2006/04/06 10:01:07 $\\

!****************************************************************
!
!**** *AUFW_ADD* - write marine add restart data.
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     
!
!     Purpose
!     -------
!     Write restart data for continuation of interrupted integration.
!
!     Method
!     -------
!     The add data are written to an extra file, other than the ocean data.
!     The netCDF version also writes grid description variables.
!     The netCDF file is selfdescribing and can be used for 
!     visualization (e.g. STOMPP, grads). Before the data are saved, all
!     values at dry cells are set to rmissing (e.g.rmasko).
!     The time stamp of the add restart file (idate) is taken from the
!     ocean time stamp through the SBR parameter list. The only time
!     control variable proper to the add is the time step number (idate(5)). 
!     It can differ from that of the ocean (idate(4)) by the difference
!     of the offsets of restart files.

!       changed : add uses ocean time
!
!**   Interface.
!     ----------
!
!     *CALL*       *AUFW_ADD(kpie,kpje,kpke,pddpo,pgila,pgiph,ptiestu)
!
!     *COMMON*     *MO_PARAM1_ADD* - declaration of ocean/sediment tracer.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd dimension) [m].
!     *REAL*    *pgila*   - geographical longitude of grid points [degree E].
!     *REAL*    *pgiph*   - geographical latitude  of grid points [degree N].
!     *REAL*    *ptiestu* - depth of layers [m].
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************

!ik nocectra is the number of all ADD element (size of ocectra(,,,l))



      USE mo_contra
      USE mo_control_add
      use mo_param1_add 

      use mo_commo1, ONLY: ie_g,je_g,ldays,lmonts,lyears
      use mo_parallel

      implicit none


      INTEGER  kpie,kpje,kpke,kpicycli
      INTEGER  kplyear,kplmon,kplday,kpldtoce
      REAL pddpo(kpie,kpje,kpke)    
      REAL pgila(kpie*2,kpje*2)
      REAL pgiph(kpie*2,kpje*2)
      REAL ptiestu(kpke)
      
      CHARACTER(LEN=80) err_text
      INTEGER  i,j,k,l,kmon,jj,ii,kk,kt,ll
      integer icode(1)
      INTEGER idate(5)
      REAL pi,rad,radi,rmissing
#ifdef PNETCDF
      INCLUDE 'netcdf.inc'
      INTEGER ncid,ncvarid,ncstat,ncoldmod,ncdims(4)                  &
     &       ,nclatid,nclonid,nclevid,nclev1id                        &
     &       ,nctraid,nccheid,ncksid,ncsedid                          &
     &       ,ncmonid,nstart2(2),ncount2(2),nstride2(2)

      INTEGER timeid, START_1D,COUNT_1D
      REAL timeaufw
      REAL zfield(kpie,kpje)


      pi        = 4.*ATAN(1.)
      rad       = pi/180.
      radi      = 1./rad

#endif
!
! Check of fields written to restart file.
!
      CALL CHCK_ADD(io_stdo_add,kpicycli,                                 &
     &'Check values of ocean addtional conservetive tracer written to restart file :',  &
     & kpie,kpje,kpke,pddpo)

      idate(1) = kplyear
      idate(2) = kplmon
      idate(3) = kplday
      idate(4) = kpldtoce
      idate(5) = ldtadd
      WRITE(io_stdo_add,*) ' '
      WRITE(io_stdo_add,*) 'Writing restart file at date :'              &
     &,'YY=',idate(1),' MM=',idate(2),' day=',idate(3)
      WRITE(io_stdo_add,*) 'Ocean model step number is ',idate(4)
      WRITE(io_stdo_add,*) 'Add   model step number is ',idate(5)
!
!  Masking aqueous sea water tracer.
!

      rmissing = rmasko
     


      DO l=1,nocectra
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
      IF(pddpo(i,j,k) .LT. 0.5) THEN
         ocectra(i,j,k,l)=rmasko
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO l=1,nicectra
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
      IF(pddpo(i,j,k) .LT. 0.5) THEN
         icectra(i,j,k,l)=rmasko
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO

!





#ifdef PNETCDF
!
! Open netCDF data file
!
      IF(p_pe==p_io) THEN

      ncstat = NF_CREATE('restartw_add.nc',NF_CLOBBER, ncid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF1')

!
! Define dimension
!
      ncstat = NF_DEF_DIM(ncid, 'lon', ie_g, nclonid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF2')

      ncstat = NF_DEF_DIM(ncid, 'lat', je_g, nclatid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF3')

      ncstat = NF_DEF_DIM(ncid, 'depth', kpke, nclevid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF4')

      ncstat = NF_DEF_DIM(ncid, 'ntra', nocectra, nctraid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF5')

      ncstat = NF_DEF_DIM(ncid, 'nche', 8, nccheid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF6')


     
! Paddy:
      ncstat = NF_DEF_DIM(ncid, 'nmon',12 , ncmonid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF8a')
     
      ncstat = NF_DEF_DIM(ncid, 'lev1', 1, nclev1id)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF8b')


      ncstat = NF_DEF_DIM(ncid, 'time',NF_UNLIMITED, timeid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('AUFW: Problem with netCDF7'      )
!

!
! Define global attributes
!
      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'title'               &
     &,35, 'restart data for marine add modules') 
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF9')

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'history'             &
     &,35, 'restart data for marine add  modules')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF9')

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'conventions'         &
     &,6, 'COARDS')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF9')

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'source'              &
     &,24, 'Marine add  model output MPI-OM/grob')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF10')

      ncstat = NF_PUT_ATT_INT(ncid, NF_GLOBAL, 'date', NF_INT, 5, idate)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF11')

!
! Define variables : grid
!
      ncdims(1) = nclonid
      ncdims(2) = nclatid

      ncstat = NF_DEF_VAR(ncid,'scal_lon',NF_DOUBLE,2,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF13')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',8, 'degree E')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF14')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,34, '2-d longitude of scalar grid cells')
      icode=1
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF15')


      ncstat = NF_DEF_VAR(ncid,'scal_lat',NF_DOUBLE,2,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF16')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',8, 'degree N')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF17')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,33, '2-d latitude of scalar grid cells')
      icode=2
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF18')

      ncstat = NF_DEF_VAR(ncid,'scal_wdep',NF_DOUBLE,2,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF16a')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',5, 'meter')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF17a')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,37, '2-d water depth at scalar grid points')
      icode=3
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF18a')




      ncdims(1) = nclevid

      ncstat = NF_DEF_VAR(ncid,'depth',NF_DOUBLE,1,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF19')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',5, 'meter')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF20')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,32, '1-d layer depths of ocean tracer')
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('AUFW: Problem with netCDF21')

      ncdims(1) = timeid

      CALL NETCDF_DEF_VARDB(ncid,4,'time',1,ncdims,ncvarid,          &
     &   16,'day as %Y%m%d.%f',4,'time',rmasko,31,io_stdo_add)



     
 

!
! Define variables : advected ocean tracer
!
      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = nclevid


     

      CALL NETCDF_DEF_VARDB(ncid,5,'h2o18',3,ncdims,ncvarid,           &
     &   9,'kmol/m**3',5,'H2O18',rmissing,25,io_stdo_add)
      icode=4
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
    
      CALL NETCDF_DEF_VARDB(ncid,5,'hDo16',3,ncdims,ncvarid,           &
     &   9,'kmol/m**3',5,'HDO16',rmissing,27,io_stdo_add)
      icode=11
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      
      CALL NETCDF_DEF_VARDB(ncid,5,'h2o16',3,ncdims,ncvarid,           &
     &   9,'kmol/m**3',5,'H2O16',rmissing,26,io_stdo_add)
      icode=7
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode) 

      CALL NETCDF_DEF_VARDB(ncid,9,'h2o18_ice',3,ncdims,ncvarid,           &
     &   9,'kmol/m**3',9,'H2O18_ICE',rmissing,28,io_stdo_add)
      icode=18
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
    
      CALL NETCDF_DEF_VARDB(ncid,9,'hDo16_ice',3,ncdims,ncvarid,           &
     &   9,'kmol/m**3',9,'HDO16_ICE',rmissing,27,io_stdo_add)
      icode=19
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode)
      
      CALL NETCDF_DEF_VARDB(ncid,9,'h2o16_ice',3,ncdims,ncvarid,           &
     &   9,'kmol/m**3',9,'H2O16_ICE',rmissing,26,io_stdo_add)
      icode=20
      ncstat = NF_PUT_ATT_INT(ncid,ncvarid, 'code', NF_INT, 1, icode) 

   

!
      ncstat = NF_ENDDEF(ncid)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF00')


!
! Set fill mode
!
      ncstat = NF_SET_FILL(ncid,NF_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF97')

      ENDIF 

!
! Write grid describing data
!

      nstart2(1) = 1
      nstart2(2) = 1
      ncount2(1) = kpie
      ncount2(2) = kpje
      nstride2(1) = 2
      nstride2(2) = 2
      DO j=1,ncount2(2)
         jj=nstart2(2)+(j-1)*nstride2(2)
         DO i=1,ncount2(1)
            ii=nstart2(1)+(i-1)*nstride2(1)
            zfield(i,j) = pgila(ii,jj)*radi
         ENDDO
      ENDDO

      CALL write_netcdf_var(ncid,'scal_lon',zfield(1,1),1,0)

               
      DO j=1,ncount2(2)
         jj=nstart2(2)+(j-1)*nstride2(2)
         DO i=1,ncount2(1)
            ii=nstart2(1)+(i-1)*nstride2(1)
            zfield(i,j) = pgiph(ii,jj)*radi
         ENDDO
      ENDDO

     CALL write_netcdf_var(ncid,'scal_lat', zfield(1,1),1,0)

      DO j=1,kpje
         DO i=1,kpie
            zfield(i,j) = 0.0
         ENDDO
      ENDDO
      DO k=1,kpke
         DO j=1,kpje
            DO i=1,kpie
               zfield(i,j) = zfield(i,j) + pddpo(i,j,k)
            ENDDO
         ENDDO
      ENDDO
             
      CALL write_netcdf_var(ncid,'scal_wdep',zfield(1,1),1,0)



      IF(p_pe==p_io) THEN
        ncstat = NF_INQ_VARID(ncid,'depth',ncvarid )
        IF ( ncstat .NE. NF_NOERR ) CALL stop_all('AUFW: netCDF102')
        ncstat = NF_PUT_VAR_DOUBLE (ncid,ncvarid,ptiestu(1) )
        IF ( ncstat .NE. NF_NOERR ) CALL stop_all('AUFW: netCDF103')




       timeaufw=LDAYS+LMONTS*100+LYEARS*10000
       START_1D =1
       COUNT_1D =1


       ncstat = NF_INQ_VARID(ncid,'time',ncvarid )
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF103d')
       ncstat = NF_PUT_VARA_DOUBLE (ncid,ncvarid,start_1d,count_1d,timeaufw)
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('AUFW: Problem with netCDF103d')


       
      ENDIF







!
! Write restart data : ocean aquateous tracer
!

#ifdef uwe_restart
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie

! replace rmasko by start value so restart file can be used for newly introduced wet points (e.g., lakes)
      IF(pddpo(i,j,k) .LT. 0.5) THEN       ! dry points
         ocectra(i,j,k,ih2o18)= 2.0052  ! [kg/m3]
         ocectra(i,j,k,ihDo16)= 0.15576  ! [kg/m3]
	     ocectra(i,j,k,ih2o16)= 1000 ! [kg/m3]
      ENDIF

      ENDDO
      ENDDO
      ENDDO
#endif/*uwe_restart*/

      ll=size(ocectra(1,1,1,:))
      do l = 1, ll
       call bounds_exch('p',ocectra(:,:,:,l),'in aufw 1')
      enddo

      ll=size(icectra(1,1,1,:))
      do l = 1, ll
       call bounds_exch('p',icectra(:,:,:,l),'in aufw 1')
      enddo


      CALL write_netcdf_var(ncid,'h2o18',ocectra(1,1,1,ih2o18),kpke,0)
      CALL write_netcdf_var(ncid,'hDo16',ocectra(1,1,1,ihDo16),kpke,0)
      CALL write_netcdf_var(ncid,'h2o16',ocectra(1,1,1,ih2o16),kpke,0)
      CALL write_netcdf_var(ncid,'h2o18_ice',icectra(1,1,1,ih2o18_ice),kpke,0)
      CALL write_netcdf_var(ncid,'hDo16_ice',icectra(1,1,1,ihDo16_ice),kpke,0)
      CALL write_netcdf_var(ncid,'h2o16_ice',icectra(1,1,1,ih2o16_ice),kpke,0)     


#else         ! no netcdf
! does this work with the MPI anymore - c.f. read_netcdf_var ??
      OPEN(io_rsto_add,FILE='restart_add',STATUS='UNKNOWN'            &
     &               ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      WRITE(io_rsto_add)                                              &
     &             (((ocectra(i,j,k,ih2o18),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((ocectra(i,j,k,ihDo16),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((ocectra(i,j,k,ih2o16),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((icectra(i,j,k,ih2o18_ice),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((icectra(i,j,k,ihDo16_ice),i=1,kpie),j=1,kpje),k=1,kpke)&
     &            ,(((icectra(i,j,k,ih2o16_ice),i=1,kpie),j=1,kpje),k=1,kpke)
      
     

      REWIND io_rsto_add

      CLOSE(io_rsto_add)
#endif

      WRITE(io_stdo_add,*) 'End of AUFW_ADD'
      WRITE(io_stdo_add,*) '***************'

      RETURN
      END
