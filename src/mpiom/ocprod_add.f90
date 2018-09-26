      SUBROUTINE OCPROD_ADD(kpie,kpje,kpke,pddpo,                     &
     &                    pdlxp,pdlyp,pdpio,kplmon)

!$Source: /scratch/local1/m212047/patrick/SRC_MPI/src_hamocc/RCS/ocprod.f90,v $\\
!$Revision: 1.1 $\\
!$Date: 2005/01/28 08:37:45 $\\

!**********************************************************************
!
!**** *OCPROD_ADD* - .
!
!    
!     Purpose
!     -------
!     compute 3d mean field
!     Method:
!     ------
!     kchck=1 can be used to check max/min of bgc arrays on wet/dry cells.
!     Note: prosil is used only for k=1,2. It is adressed, however, for
!           k=1,4 in loop 100. To save memory,  js:note seems dated???
!
!     _ant fields are natural PLUS anthropogenic (not anthropogenic only!!!)
!
!     *CALL*       *OCPROD*
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *ptho*    - potential temperature [deg C].
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st dimension) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd dimension) [m].
!     *REAL*    *pdpio*   - inverse thickness of grid cell (3rd dimension)[m].
!
!     Externals
!     ---------
!     .
!**********************************************************************
! nocectra is the number of all ADD elements (size of ocetra(,,,l))

      USE mo_timeser_add
      USE mo_contra
     
      use mo_param1_add

      USE mo_control_add
      USE mo_addmean

      use mo_parallel

implicit none

      INTEGER :: kplmon,kpie,kpje,kpke
      INTEGER :: i,j,k,l
      
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: pdpio(kpie,kpje,kpke)
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
     
      
      REAL ::  volcell
      
     
      
!

      call maschk_add(kpie,kpje,kpke,0)

	    
! sum-up 3d tracers, averaged annually in avrg_addmean_3d.f90
! these go into addmean_3d files
!      WRITE(0,*) 'Max of pddpo ', MAXVAL(pddpo(:,:,:))

      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
         IF(pddpo(i,j,k).GT.0.5) THEN
	 
             addt3d(i,j,k,jh2o18_t) =  		                       &
     &       addt3d(i,j,k,jh2o18_t) + ocectra(i,j,k,ih2o18)
             addt3d(i,j,k,jhDo16_t) =  		                       &
     &       addt3d(i,j,k,jhDo16_t) + ocectra(i,j,k,ihDo16)
             addt3d(i,j,k,jh2o16_t) =  		                       &
     &       addt3d(i,j,k,jh2o16_t) + ocectra(i,j,k,ih2o16)
             addice(i,j,k,jh2o18_ice) =  		                       &
     &       addice(i,j,k,jh2o18_ice) + icectra(i,j,k,ih2o18_ice)
             addice(i,j,k,jhDo16_ice) =  		                       &
     &       addice(i,j,k,jhDo16_ice) + icectra(i,j,k,ihDo16_ice)
             addice(i,j,k,jh2o16_ice) =  		                       &
     &       addice(i,j,k,jh2o16_ice) + icectra(i,j,k,ih2o16_ice)
            
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      call maschk_add(kpie,kpje,kpke,-4)

!
! Sampling timeseries-1 : global inventory (this will be stored in same file as station data with index 0)
!

      DO k=1,kpke
      DO j=2,kpje-1
         DO i=2,kpie-1
            IF(pddpo(i,j,k).GT.0.5) THEN

               volcell = pdlxp(i,j)*pdlyp(i,j)*pddpo(i,j,k)*1.e-6          ! why *1.e-6? pdlxp [m], pdlyp [m], pddpo [m]
                                                                           ! does not fit unit in output file [kmol/m3]
                                                                           ! unit in output file should be [10^6 kmol P /day]

               ts1(itsh2o18,1,lts1) = ts1(itsh2o18,1,lts1)           &   ! lts1 counts time steps for sampling
     &                               + ocectra(i,j,k,ih2o18)*volcell

               ts1(itshDo16,1,lts1) = ts1(itshDo16,1,lts1)           &   
     &                               + ocectra(i,j,k,ihDo16)*volcell
             
               ts1(itsh2o16,1,lts1) = ts1(itsh2o16,1,lts1)           & 
     &                               + ocectra(i,j,k,ih2o16)*volcell
               ts1(itsh2o18_ice,1,lts1) = ts1(itsh2o18_ice,1,lts1)           &   ! lts1 counts time steps for sampling
     &                               + icectra(i,j,k,ih2o18_ice)*volcell

               ts1(itshDo16_ice,1,lts1) = ts1(itshDo16_ice,1,lts1)           &   
     &                               + icectra(i,j,k,ihDo16_ice)*volcell
             
               ts1(itsh2o16_ice,1,lts1) = ts1(itsh2o16_ice,1,lts1)           & 
     &                               + icectra(i,j,k,ih2o16_ice)*volcell
               
               
            ENDIF
         ENDDO
      ENDDO
      ENDDO
!
      call maschk_add(kpie,kpje,kpke,-6)

! Check maximum/minimum values in wet/dry cells.
!
      IF( kchck_add .EQ. 1) CALL CHCK_ADD(io_stdo_add,icycliadd,           &
     &'Check values of ocean tracer at exit from SBR OCPROD_ADD :',        &
     & kpie,kpje,kpke,pddpo)


      RETURN
      END
