      SUBROUTINE WRITE_ADDMEAN_3D(kpie,kpje,kpke,pddpo)
      
!****************************************************************
!
!**** *WRITE_BGCMEAN_2D* - calculate and write 3-dimensional add mean data.
!

!     Purpose
!     -------
!     Write add mean data.
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *WRITE_ADDMEAN_3D(kpie,kpje,kpke,pddpo)*
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
      use mo_commo1, only: lyears,lmonts,ldays

      USE mo_control_add

      use mo_parallel
 
      implicit none
      INTEGER :: kpie,kpje,kpke,i,j,k,l,nk
      REAL :: pddpo(kpie,kpje,kpke)
      real :: time3d

      INTEGER :: START_3D(4),COUNT_3D(4) ! size of internal array 
      INTEGER :: START_1D,COUNT_1D


      INCLUDE 'netcdf.inc'
      INTEGER :: ncvarid,ncstat,ncoldmod


!-----------------------------------------------------------------------
      START_1D =meantime_2d_add
     
      COUNT_1D =1

      START_3D(1) =1  
      START_3D(2) =1
      START_3D(3) =1
      START_3D(4) =meantime_3d_add
      
      COUNT_3D(1) =kpie  
      COUNT_3D(2) =kpje
      COUNT_3D(3) =kpke
      COUNT_3D(4) =1      
      
!-----------------------------------------------------------------------

 
      WRITE(io_stdo_add,*) 'Writing 3D addmean data at step:', meantime_3d_add

!
!  Masking addmean data.
!
!-----------------------------------------------------------------------

      DO l=1,naddt3d 
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie      
        IF(pddpo(i,j,k) .LT. 0.5) THEN
          addt3d(i,j,k,l)=rmasko
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO l=1,naddice 
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie      
        IF(pddpo(i,j,k) .LT. 0.5) THEN
          addice(i,j,k,l)=rmasko
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
        

!
! Set fill mode
!
!-----------------------------------------------------------------------
      if (p_pe==p_io) then

      ncstat = NF_SET_FILL(nc_3d_id_add,NF_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('ADDMEAN: Problem with netCDF97')          

       time3d=LDAYS+LMONTS*100+LYEARS*10000


       ncstat = NF_INQ_VARID(nc_3d_id_add,'time',ncvarid )
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('WRITE_ADDMEAN_3D: Problem with netCDF103d')
       ncstat = NF_PUT_VARA_DOUBLE (nc_3d_id_add,ncvarid,start_1d,count_1d,time3d)
       IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('WRITE_ADDMEAN_3D: Problem with netCDF103d')

      end if
!
! Write bgcmean data : total 3-D add data
!
!-----------------------------------------------------------------------
      WRITE(0,*) 'Max of addt3d ', MAXVAL(addt3d(:,:,:,jh2o18_t))
      nk = kpke
      CALL write_netcdf_var(nc_3d_id_add,'h2o18_t',addt3d(1,1,1,jh2o18_t),nk,meantime_3d_add)
      CALL write_netcdf_var(nc_3d_id_add,'hDo16_t',addt3d(1,1,1,jhDo16_t),nk,meantime_3d_add) 
      CALL write_netcdf_var(nc_3d_id_add,'h2o16_t',addt3d(1,1,1,jh2o16_t),nk,meantime_3d_add)
      CALL write_netcdf_var(nc_3d_id_add,'h2o18_ice',addice(1,1,1,jh2o18_ice),nk,meantime_3d_add)
      CALL write_netcdf_var(nc_3d_id_add,'hDo16_ice',addice(1,1,1,jhDo16_ice),nk,meantime_3d_add) 
      CALL write_netcdf_var(nc_3d_id_add,'h2o16_ice',addice(1,1,1,jh2o16_ice),nk,meantime_3d_add)
     
!
! Reset mean fields
!
!-----------------------------------------------------------------------
!
      DO l=1,naddt3d
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
          addt3d(i,j,k,l) = 0.	
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO l=1,naddice
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
          addice(i,j,k,l) = 0.	
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      
      RETURN
      END
