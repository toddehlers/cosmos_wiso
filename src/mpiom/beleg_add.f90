      SUBROUTINE BELEG_ADD(kpie,kpje,kpke,pddpo)

!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/beleg_bgc.f90,v $\\
!$Revision: 1.2.2.1.16.1.2.1.2.1 $\\
!$Date: 2006/04/03 11:27:49 $\\

!****************************************************************
!
!**** *BELEG_ADD* - initialize add variables.
!
!     
!     
!     Purpose
!     -------
!     - set start values for  add variables.
!
!     Method
!     -------
!     - add variables are initialized. They might be overwritten
!       if read from restart by call of AUFR_ADD.
!     - physical constants are defined
!     - fields used for budget calculations (should be in extra SBR!)
!       and as boundary conditions are set to zero.
!     
!
!**   Interface.
!     ----------
!
!     *CALL*       *BELEG_ADD*
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *PARAM1_ADD.h* - declaration of ocean/sediment tracer.
!     *COMMON*     *COMMO1_ADD.h* - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS_ADD.h*  - std I/O logical units.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd dimension) [m].
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************

      USE mo_contra
      USE mo_control_add
      USE mo_addmean
      use mo_param1_add 
      
      USE MO_COMMO1

      use mo_parallel



implicit none      
            
      

      INTEGER :: i,j,k,l,ii,jj,m,kpie,kpje,kpke,kmon
      INTEGER :: IF1,IF2,IF3,IF4
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: dummy(kpie,kpje)
      REAL :: north, south
      


      REAL :: xpi,rad,radi,rmissing

      xpi       = 4.*ATAN(1.)
      rad       = xpi/180.
      radi      = 1./rad


!      WRITE(0,*) 'Start of beleg_add'
!
!  Initialize overall time step counter.
!
      ldtadd = 0
!
!  Initialize 3D time step counter.
!
     
      meancnt_add_3D = 0

      
      WRITE(0,*) 'Start of initializtion first layer'
      
!
!  Initial values for aquatic (advected) ocean tracers (for new start, start with standard mean)
! 
       
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
 	
     
      IF(pddpo(i,j,k) .GT. 0.5) THEN     ! wet points
 
         ocectra(i,j,k,ih2o18)= 2.0052           ! [kg/m3]
         ocectra(i,j,k,ihDo16)= 0.15576          ! [kg/m3]
         ocectra(i,j,k,ih2o16)= 1000.0           ! [kg/m3]
         icectra(i,j,k,ih2o18_ice)= 0.
         icectra(i,j,k,ihDo16_ice)= 0.
         icectra(i,j,k,ih2o16_ice)= 0.
 
      ELSE                              ! dry points
         ocectra(i,j,k,ih2o18)= rmasko
         ocectra(i,j,k,ihDo16)= rmasko
         ocectra(i,j,k,ih2o16)= rmasko

      ENDIF
      
      


      ENDDO
      ENDDO
      ENDDO


!  Initial values for aquatic (advected) ocean tracers (for new start, start with noaa set base)
 

!      OPEN(211, FILE='/home/sx8/xxu/mpiom-1.3.0-O18_icetest/setup/mk_phc/GR60 GR60L20_INID18O_PHC', &
!      &     ACCESS='SEQUENTIAL', FORM='UNFORMATTED')       

!      DO k=1,kpke
!      READ(211)IF1,IF2,IF3,IF4
!      READ(211)((dummy(ii,jj),ii=1,kpie),jj=kpje,1,-1)
!      DO j=1,kpje
!      DO i=1,kpie
      
     
!      IF(pddpo(i,j,k) .GT. 0.5) THEN     ! wet points
! 
!          ocectra(i,j,k,ih2o18)= dummy(i,j)        ! [kg/m3]
! 	  ocectra(i,j,k,ih2o16)= 1000.             ! [kg/m3]
!      ELSE                              ! dry points
!          ocectra(i,j,k,ih2o18)= rmasko
!	  ocectra(i,j,k,ih2o16)= rmasko
!
!      ENDIF


!      ENDDO
!      ENDDO
!      ENDDO



     WRITE(0,*) 'End of initial ocectra'
!
      
!
! Values for addmean
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
