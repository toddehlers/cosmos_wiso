      SUBROUTINE OPEN_ADDMEAN_3D(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,   &
     &        pgila,pgiph,ptiestu,kplyear,kplmon,kplday,kpldtoce)

!****************************************************************
!
!**** OPEN_ADDMEAN_3D* - open netcdf files for 3-dimensional add mean data.
!
!    
!
!     Purpose
!     -------
!     open netcdf files for 3-dimensional add mean data
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *OPEN_ADDMEAN_3D(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,
!                  pgila,pgiph,ptiestu,kplyear,kplmon,kplday,kpldtoce) *
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st REAL :: of model grid.
!     *INTEGER* *kpje*    - 2nd REAL :: of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) REAL :: of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd REAL ::) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st REAL ::) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd REAL ::) [m].
!     *REAL*    *pgila*   - geographical longitude of grid points [degree E].
!     *REAL*    *pgiph*   - geographical latitude  of grid points [degree N].
!     *REAL*    *ptiestu* - depth of layers [m].
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************

      USE mo_addmean
      USE mo_contra
      use mo_param1_add 

      USE mo_control_add

      USE mo_commo1, ONLY: ie_g, je_g
      USE mo_parallel
 
      implicit none
      INTEGER :: kpie,kpje,kpke,i,j,k,ii,jj
      INTEGER :: kplyear,kplmon,kplday,kpldtoce
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL :: pgila(kpie*2,kpje*2)
      REAL :: pgiph(kpie*2,kpje*2)
      REAL :: ptiestu(kpke)


      INCLUDE 'netcdf.inc'
      INTEGER :: ncvarid,ncstat,ncoldmod,ncdims(4)
      INTEGER :: nclatid,nclonid,nclevid
      INTEGER :: timeid,nstart2(2),ncount2(2),nstride2(2)
      INTEGER :: icode(1)
      INTEGER :: idate(5)
      INTEGER :: chunksize,initialsize

      REAL :: zfield(kpie,kpje)

      REAL :: pi,rad,radi

!-----------------------------------------------------------------------

      pi        = 4.*ATAN(1.)
      rad       = pi/180.
      radi      = 1./rad


      idate(1) = kplyear
      idate(2) = kplmon
      idate(3) = kplday
      idate(4) = kpldtoce
      idate(5) = ldtadd

      IF (p_pe==p_io) THEN
      WRITE(io_stdo_add,*) ' '
      WRITE(io_stdo_add,*) 'Creating 3D addmean data at date :'    &
     &  ,'YY=',idate(1),' MM=',idate(2),' day=',idate(3)
      WRITE(io_stdo_add,*) 'Ocean model step number is ',idate(4)
      WRITE(io_stdo_add,*) 'Add   model step number is ',idate(5)


!
! Open netCDF data file
!
!-----------------------------------------------------------------------
      
!     WRITE(0,*) 'ncstat ', ncstat
      ncstat = NF_CREATE('addmean_3d.nc',NF_CLOBBER, nc_3d_id_add)
!     WRITE(0,*) 'ncstat ', ncstat, ' NF_NOERR ', NF_NOERR
!     WRITE(0,*)  NF_STRERROR(ncstat)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF1')

      chunksize = 1024*1024*32 ! 32 MB --> man 3f netcdf
      initialsize = 1024*1024*32 ! 32 MB --> man 3f netcdf
      

!
! Define dimension
!
!-----------------------------------------------------------------------

      ncstat = NF_DEF_DIM(nc_3d_id_add, 'lon', ie_g, nclonid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF2')

      ncstat = NF_DEF_DIM(nc_3d_id_add, 'lat', je_g, nclatid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF3')

      ncstat = NF_DEF_DIM(nc_3d_id_add, 'depth', kpke, nclevid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF4')

      ncstat = NF_DEF_DIM(nc_3d_id_add, 'time',NF_UNLIMITED, timeid)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF7'      )
!
! Define global attributes
!
!-----------------------------------------------------------------------
 
      ncstat = NF_PUT_ATT_TEXT(nc_3d_id_add, NF_GLOBAL,'title'             &
     &,42, 'Yearly mean output for marine add modules')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF8')

      ncstat = NF_PUT_ATT_TEXT(nc_3d_id_add, NF_GLOBAL,'history'           &
     &,42, 'Yearly mean output for marine add modules')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF9')

      ncstat = NF_PUT_ATT_TEXT(nc_3d_id_add, NF_GLOBAL,'conventions'       &
     &,6, 'COARDS')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF10')

      ncstat = NF_PUT_ATT_TEXT(nc_3d_id_add, NF_GLOBAL,'source'            &
     &,24, 'Marine add model output HOPC68/grob')
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF11')

      ncstat = NF_PUT_ATT_INT(nc_3d_id_add, NF_GLOBAL, 'date', NF_INT, 5, idate)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF11')



!-----------------------------------------------------------------------

      ncdims(1) = timeid

      CALL NETCDF_DEF_VARDB(nc_3d_id_add,4,'time',1,ncdims,ncvarid,          &
     &   16,'day as %Y%m%d.%f',4,'time',rmasko,31,io_stdo_add)

!
! Define variables : grid
!
!-----------------------------------------------------------------------

      ncdims(1) = nclonid
      ncdims(2) = nclatid

      CALL NETCDF_DEF_VARSG(nc_3d_id_add,8,'scal_lon',2,ncdims,ncvarid,          &
     &   8,'degree E',34,'2-d longitude of scalar grid cells',rmasko,20,io_stdo_add) 
      icode=1
      ncstat = NF_PUT_ATT_INT(nc_3d_id_add,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_3d_id_add,8,'scal_lat',2,ncdims,ncvarid,          &
     &   8,'degree N',33,'2-d latitude of scalar grid cells',rmasko,21,io_stdo_add) 
      icode=2
      ncstat = NF_PUT_ATT_INT(nc_3d_id_add,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_3d_id_add,9,'scal_wdep',2,ncdims,ncvarid,          &
     &   5,'meter',37,'2-d water depth at scalar grid points',rmasko,22,io_stdo_add) 
      icode=3
      ncstat = NF_PUT_ATT_INT(nc_3d_id_add,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_3d_id_add,6,'size_x',2,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in x-direction',rmasko,23,io_stdo_add) 
      icode=4
      ncstat = NF_PUT_ATT_INT(nc_3d_id_add,ncvarid, 'code', NF_INT, 1, icode)

      CALL NETCDF_DEF_VARSG(nc_3d_id_add,6,'size_y',2,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in y-direction',rmasko,24,io_stdo_add) 
      icode=5
      ncstat = NF_PUT_ATT_INT(nc_3d_id_add,ncvarid, 'code', NF_INT, 1, icode)

      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = nclevid

      CALL NETCDF_DEF_VARSG(nc_3d_id_add,6,'size_z',3,ncdims,ncvarid,          &
     &   5,'meter',33,'size of grid cells in z-direction',rmasko,25,io_stdo_add)  
      icode=6
      ncstat = NF_PUT_ATT_INT(nc_3d_id_add,ncvarid, 'code', NF_INT, 1, icode)          
      
      ncdims(1) = nclevid

      CALL NETCDF_DEF_VARSG(nc_3d_id_add,5,'depth',1,ncdims,ncvarid,          &
     &   5,'meter',32,'1-d layer depths of ocean tracer',rmasko,26,io_stdo_add)  


!
! Define variables : total add data (js: what means 'total' here?)
!
!-----------------------------------------------------------------------

      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = nclevid
      ncdims(4) = timeid




      
     
        CALL NETCDF_DEF_VARSG(nc_3d_id_add,7,'h2o18_t',4,ncdims,ncvarid,    & 
     &    9,'kmol/m**3',19,'h2o18 concentration',rmasko,46,io_stdo_add)         
      icode=12
      ncstat = NF_PUT_ATT_INT(nc_3d_id_add,ncvarid, 'code', NF_INT, 1, icode)    

        CALL NETCDF_DEF_VARSG(nc_3d_id_add,7,'hDo16_t',4,ncdims,ncvarid,    & 
     &    9,'kmol/m**3',19,'hDo16 concentration',rmasko,48,io_stdo_add)         
      icode=14
      ncstat = NF_PUT_ATT_INT(nc_3d_id_add,ncvarid, 'code', NF_INT, 1, icode)
        
        CALL NETCDF_DEF_VARSG(nc_3d_id_add,7,'h2o16_t',4,ncdims,ncvarid,     &
     &    9,'kmol/m**3',19,'h2o16 concentration',rmasko,47,io_stdo_add)
      icode=13
      ncstat = NF_PUT_ATT_INT(nc_3d_id_add,ncvarid, 'code', NF_INT, 1, icode)
        
        CALL NETCDF_DEF_VARSG(nc_3d_id_add,9,'h2o18_ice',4,ncdims,ncvarid,    & 
     &    9,'kmol/m**3',23,'h2o18_ice concentration',rmasko,49,io_stdo_add)         
      icode=15
      ncstat = NF_PUT_ATT_INT(nc_3d_id_add,ncvarid, 'code', NF_INT, 1, icode)    

        CALL NETCDF_DEF_VARSG(nc_3d_id_add,9,'hDo16_ice',4,ncdims,ncvarid,    & 
     &    9,'kmol/m**3',23,'hDo16_ice concentration',rmasko,50,io_stdo_add)         
      icode=16
      ncstat = NF_PUT_ATT_INT(nc_3d_id_add,ncvarid, 'code', NF_INT, 1, icode)
        
        CALL NETCDF_DEF_VARSG(nc_3d_id_add,9,'h2o16_ice',4,ncdims,ncvarid,     &
     &    9,'kmol/m**3',23,'h2o16_ice concentration',rmasko,51,io_stdo_add)
      icode=17
      ncstat = NF_PUT_ATT_INT(nc_3d_id_add,ncvarid, 'code', NF_INT, 1, icode)

!
! END Define variables
!
!-----------------------------------------------------------------------

      ncstat = NF_ENDDEF(nc_3d_id_add)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF00')
!
! Set fill mode
!
!-----------------------------------------------------------------------

      ncstat = NF_SET_FILL(nc_3d_id_add,NF_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF97')

      END IF ! p_pe==p_io

!
! Write grid describing data
!
!-----------------------------------------------------------------------

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
      CALL write_netcdf_var(nc_3d_id_add,'scal_lon',zfield(1,1),1,0)
               
      DO j=1,ncount2(2)
         jj=nstart2(2)+(j-1)*nstride2(2)
         DO i=1,ncount2(1)
            ii=nstart2(1)+(i-1)*nstride2(1)
            zfield(i,j) = pgiph(ii,jj)*radi
         ENDDO
      ENDDO
      CALL write_netcdf_var(nc_3d_id_add,'scal_lat', zfield(1,1),1,0)

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
      CALL write_netcdf_var(nc_3d_id_add,'scal_wdep',zfield(1,1),1,0)

      IF(p_pe==p_io) THEN 
        ncstat = NF_INQ_VARID(nc_3d_id_add,'depth',ncvarid )
        IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF102a')
        ncstat = NF_PUT_VAR_DOUBLE (nc_3d_id_add,ncvarid,ptiestu(1) )
        IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF103a')
      ENDIF

      CALL write_netcdf_var(nc_3d_id_add,'size_x',pdlxp(1,1),1,0)
      CALL write_netcdf_var(nc_3d_id_add,'size_y',pdlyp(1,1),1,0)
      CALL write_netcdf_var(nc_3d_id_add,'size_z',pddpo(1,1,1),kpke,0)

!
! Close File
!
!-----------------------------------------------------------------------
!      IF(p_pe==p_io) THEN
!        ncstat = NF_CLOSE(nc_3d_id_add)
!        IF ( ncstat .NE. NF_NOERR ) call stop_all('OPEN_ADDMEAN_3D: Problem with netCDF200')
!      ENDIF
!

      RETURN
      END
