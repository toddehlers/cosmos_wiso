      SUBROUTINE AUFR_ADD(kpie,kpje,kpke,pddpo,kplyear,kplmon     &
     &                    ,kplday,kpldtoce)

!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/aufr_bgc.f90,v $\\
!$Revision: 1.3.2.1.16.1.4.2 $\\
!$Date: 2006/04/06 10:01:07 $\\

!****************************************************************
!
!**** *AUFR_ADD* - reads marine add restart data.
!
!     
!
!     Purpose
!     -------
!     Read restart data to continue an interrupted integration.
!
!     Method
!     -------
!     The add data are read from a file other than the ocean data.
!     The time stamp of the add restart file (idate) is specified from the
!     ocean time stamp through the SBR parameter list of AUFW_ADD. The only 
!     time control variable proper to the add is the time step number 
!     (idate(5)). It can differ from that of the ocean (idate(4)) by the 
!     difference of the offsets of restart files.
!
!**   Interface.
!     ----------
!
!     *CALL*       *AUFR_ADD(kpie,kpje,kpke,pddpo
!                            ,kplyear,kplmon,kplday,kpldtoce)*
!
!     *COMMON*     *MO_PARAM1_ADD* - declaration of ocean/sediment tracer.
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
!     *INTEGER* *kplyear*   - year  in ocean restart date
!     *INTEGER* *kplmon*  - month in ocean restart date
!     *INTEGER* *kplday*    - day   in ocean restart date
!     *INTEGER* *kpldtoce*  - step  in ocean restart date
!
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************
!ik nocetra is the number of all ADD elements (size of ocetra(,,,l))
!ik nosedi is the number of all elements interacting with the sediment

      USE mo_contra
      USE mo_control_add
      use mo_param1_add 

      use mo_mpi
      use mo_parallel
      
      implicit none
      
      INTEGER  kpie,kpje,kpke
      INTEGER  kplyear,kplmon,kplday,kpldtoce
      REAL pddpo(kpie,kpje,kpke)
      INTEGER  i,j,k,l,kmon

      INTEGER idate(5)

      INTEGER :: restyear            !  year of restart file
      INTEGER :: restmonth           !  month of restart file
      INTEGER :: restday             !  day of restart file
      INTEGER :: restdtoce           !  time step number from add ocean file

#ifdef PNETCDF                
      INCLUDE 'netcdf.inc'
      INTEGER ncid,ncstat,ncvarid
!
! Open netCDF data file
!
!js 2 lines from emr deactivated for production runs

      IF(p_pe==p_io) THEN
        ncstat = NF_OPEN('restartr_add.nc',NF_NOWRITE, ncid)
        IF ( ncstat .NE. NF_NOERR ) CALL STOP_ALL('AUFR: Problem with netCDF1')

!
! Read restart data : date
!

        ncstat = NF_GET_ATT_INT(ncid, NF_GLOBAL,'date', idate)
        IF ( ncstat .NE. NF_NOERR )                                      &
        &  CALL STOP_ALL('AUFR: Problem reading date of restart file')
        ENDIF
      CALL p_bcast(idate,p_io)
      restyear  = idate(1)
      restmonth = idate(2)
      restday   = idate(3)
      restdtoce = idate(4)
      ldtadd = idate(5)
      WRITE(io_stdo_add,*) ' '
      WRITE(io_stdo_add,*) 'Date of add restart file : '
      WRITE(io_stdo_add,*) ' year  = ',restyear
      WRITE(io_stdo_add,*) ' month = ',restmonth
      WRITE(io_stdo_add,*) ' day   = ',restday
      WRITE(io_stdo_add,*) ' oce_step = ',restdtoce
      WRITE(io_stdo_add,*) ' add_step = ',ldtadd
      WRITE(io_stdo_add,*) ' '

!
! Compare with date read from ocean restart file
!
! As the ocean is already in its first step, its counter has 
! gone up one step already for the year and month. The ocean day
! counter is still at its restart date. Therefore:



! Memorize ocean start date :   
      addstartyear  = kplyear
      addstartmonth = kplmon
      addstartday   = kplday

      IF ( kplyear  .NE. restyear  ) THEN
         WRITE(io_stdo_add,*)                                     &
     &   'WARNING: restart years in oce/add are not the same : '  &
     &   ,kplyear,'/',restyear,' !!!'
      ENDIF

      IF ( kplmon .NE. restmonth ) THEN
         WRITE(io_stdo_add,*)                                     &
     &   'WARNING: restart months in oce/add are not the same : '   &
     &   ,kplmon,'/',restmonth,' !!!'
!         STOP 'Stop : restart months in oce/add are not the same.'
      ENDIF

      IF ( kplday   .NE. restday   ) THEN
         WRITE(io_stdo_add,*)                                     &
     &   'WARNING: restart days in oce/add are not the same : '   &
     &   ,kplday,'/',restday,' !!!'
!         STOP 'Stop : restart days in oce/add are not the same.'
      ENDIF 

      IF ( kpldtoce .NE. ldtadd   ) THEN
         WRITE(io_stdo_add,*)                                       & 
     &   'WARNING: restart step no.  in oce/add are not the same : '&
     &   ,kpldtoce,'/',ldtadd,' !!!'
      ENDIF

!
! Read restart data : ocean aqueous tracer
!                
      CALL read_netcdf_var(ncid,'h2o18',ocectra(1,1,1,ih2o18),kpke)
      CALL read_netcdf_var(ncid,'hDo16',ocectra(1,1,1,ihDo16),kpke)
      CALL read_netcdf_var(ncid,'h2o16',ocectra(1,1,1,ih2o16),kpke)
      CALL read_netcdf_var(ncid,'h2o18_ice',icectra(1,1,1,ih2o18_ice),kpke)
      CALL read_netcdf_var(ncid,'hDo16_ice',icectra(1,1,1,ihDo16_ice),kpke)
      CALL read_netcdf_var(ncid,'h2o16_ice',icectra(1,1,1,ih2o16_ice),kpke)

     


!
!Check aquateous restart data for topography
! 

      DO i    =1,kpie
      DO j    =1,kpje
      DO k    =1,kpke      
      DO l    =1,nocectra 
         IF (pddpo(i,j,k) .le. 0.5 ) THEN
            IF (kchck_add == 1) THEN
               IF ( ocectra(i,j,k,l) .NE. rmasko ) THEN
                  WRITE(io_stdo_add,*) 'ocectra not properly masked at :'  &
     &                              ,i,j,k,l,ocectra(i,j,k,l)
               END IF
            END IF
         ELSE
            IF ( ocectra(i,j,k,l) .EQ. rmasko ) THEN
               WRITE(io_stdo_add,*) 'land mask values at wet points :'  &
     &                             ,i,j,k,l,ocectra(i,j,k,l)

               IF (pddpo(i-1,j,1) .gt.0.5) ocectra(i,j,k,l)=ocectra(i-1,j,1,l)
               IF (pddpo(i+1,j,1) .gt.0.5) ocectra(i,j,k,l)=ocectra(i+1,j,1,l)
               IF (pddpo(i,j+1,1) .gt.0.5) ocectra(i,j,k,l)=ocectra(i,j+1,1,l)
               IF (pddpo(i,j-1,1) .gt.0.5) ocectra(i,j,k,l)=ocectra(i,j-1,1,l)
            END IF
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO


      DO i    =1,kpie
      DO j    =1,kpje
      DO k    =1,kpke      
      DO l    =1,nicectra 
         IF (pddpo(i,j,k) .le. 0.5 ) THEN
            IF (kchck_add == 1) THEN
               IF ( icectra(i,j,k,l) .NE. rmasko ) THEN
                  WRITE(io_stdo_add,*) 'icectra not properly masked at :'  &
     &                              ,i,j,k,l,icectra(i,j,k,l)
               END IF
            END IF
         ELSE
            IF ( icectra(i,j,k,l) .EQ. rmasko ) THEN
               WRITE(io_stdo_add,*) 'land mask values at wet points :'  &
     &                             ,i,j,k,l,icectra(i,j,k,l)

               IF (pddpo(i-1,j,1) .gt.0.5) icectra(i,j,k,l)=icectra(i-1,j,1,l)
               IF (pddpo(i+1,j,1) .gt.0.5) icectra(i,j,k,l)=icectra(i+1,j,1,l)
               IF (pddpo(i,j+1,1) .gt.0.5) icectra(i,j,k,l)=icectra(i,j+1,1,l)
               IF (pddpo(i,j-1,1) .gt.0.5) icectra(i,j,k,l)=icectra(i,j-1,1,l)
            END IF
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO



#else  !! Attention !! Update required !!!!!

     ! Attention - this doesn't work any more for MPI
     ! One big FORTRAN READ like here is a very, very bad idea
     ! for parallel programs !!!!!

#ifndef NOMPI
     call stop_all("we can't read in add restart data with MPI and no NETCDF now")
! would need to read these in on p_io and scatter if we were actually i
! going to do this - c.f. read_netcdf_var.F90
#else

      OPEN(io_rsti_add,FILE='restart_add',STATUS='UNKNOWN'              &
     &               ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      READ(io_rsti_add)                                                 &
     &            (((ocectra(i,j,k,ih2o18),i=1,kpie),j=1,kpje),k=1,kpke)&
     &           ,(((ocectra(i,j,k,ihDo16),i=1,kpie),j=1,kpje),k=1,kpke)&
     &           ,(((ocectra(i,j,k,ih2o16),i=1,kpie),j=1,kpie),k=1,kpke)&
     &           ,(((icectra(i,j,k,ih2o18_ice),i=1,kpie),j=1,kpie),k=1,kpke)&
     &           ,(((icectra(i,j,k,ihDo16_ice),i=1,kpie),j=1,kpie),k=1,kpke)&
     &           ,(((icectra(i,j,k,ih2o16_ice),i=1,kpie),j=1,kpie),k=1,kpke)  

      CLOSE (io_rsti_add)
#endif

#endif
!
!     
!  Masking aqueous sea water tracer.
!
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


      return
!     
!  Restrict to positive values (until error is found only !!!!) js: which error?
!
      DO l=1,nocectra
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
      IF(pddpo(i,j,k) .GT. 0.5) THEN
	 ocectra(i,j,k,l)=MAX(ocectra(i,j,k,l),0.)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO l=1,nicectra
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
      IF(pddpo(i,j,k) .GT. 0.5) THEN
	 icectra(i,j,k,l)=MAX(icectra(i,j,k,l),0.)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END
