      SUBROUTINE CHCK_ADD(kunit,kpicycli,ystring,kpie,kpje,kpke,pddpo)

!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/chck_bgc.f90,v $\\
!$Revision: 1.2.2.1.4.1.2.2.4.1.2.2.2.3.2.1 $\\
!$Name: mpiom_1_2_0 $\\
!
!*********************************************************************
!
!**** *CHCK_ADD* - check fields
!
!     
!     
!     Purpose
!     -------
!     - check max/min values on wet/dry cells.
!
!     Method
!     -------
!     -
!
!**   Interface.
!     ----------
!
!     *CALL*       *CHCK_ADD(kunit,kpicycli,ystring,kpie,kpje,kpke,pddpo)*
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *PARAM1_ADD.h* - declaration of ocean/sediment tracer.
!     *COMMON*     *COMMO1_ADD.h* - ocean/sediment tracer arrays.
!
!     called by:
!
!               aufw_add
!               ocprod
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kunit*    - stdout logical unit.
!     *INTEGER* *kpicycli* - flag for cyclicity.
!     *CHARACTER* *ystring*- .
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!
!     Externals
!     ---------
!     EXTR
!
!***********************************************************************

      USE mo_contra

      use mo_param1_add 

      USE mo_control_add

      USE mo_param1, ONLY: ie_g

implicit none

      CHARACTER*(*) ystring
      INTEGER :: kunit,kpicycli,kpie,kpje,kpke,i,j,k,l
      REAL ::  pddpo(kpie,kpje,kpke)

      WRITE(kunit,*) ystring

                         
! check for invalid values of oceanic tracers in dry and wet cells
!
      DO l=1,nocectra
         WRITE(kunit,*)                                                &
     &   'Check values of ocean additional conservitive tracer no. ',l,' :'
      DO k=1,kpke,5
         WRITE(kunit,*)'      layer ',k,' :'
         CALL EXTR(kpie,kpje,ocectra(1,1,k,l),pddpo(1,1,k),rmasko,kunit)
      ENDDO
      ENDDO

      DO l=1,nicectra
         WRITE(kunit,*)                                                &
     &   'Check values of ocean ice additional conservitive tracer no. ',l,' :'
      DO k=1,kpke,5
         WRITE(kunit,*)'      layer ',k,' :'
         CALL EXTR(kpie,kpje,icectra(1,1,k,l),pddpo(1,1,k),rmasko,kunit)
      ENDDO
      ENDDO

!                     
!                        
! check for cyclicity: ocean addtional conservitive tracer.
!
      DO l=1,nocectra
      DO k=1,kpke
      DO j=1,kpje
      IF (ABS(ocectra(1,j,k,l)-ocectra(kpie-1,j,k,l)).GT.1.e-15 .OR.     &
     &    ABS(ocectra(2,j,k,l)-ocectra(kpie  ,j,k,l)).GT.1.e-15     ) THEN
          WRITE(kunit,*)                                                &
     &   'Ocean additional conservertive tracer no. ',l,' is not cyclic at j,k=',j,k,' : ',     &
     &    ocectra(1,j,k,l),ocectra(kpie-1,j,k,l),                        &
     &    ocectra(2,j,k,l),ocectra(kpie  ,j,k,l)
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      DO l=1,nicectra
      DO k=1,kpke
      DO j=1,kpje
      IF (ABS(icectra(1,j,k,l)-icectra(kpie-1,j,k,l)).GT.1.e-15 .OR.     &
     &    ABS(icectra(2,j,k,l)-icectra(kpie  ,j,k,l)).GT.1.e-15     ) THEN
          WRITE(kunit,*)                                                &
     &   'Ocean additional conservertive tracer no. ',l,' is not cyclic at j,k=',j,k,' : ',     &
     &    icectra(1,j,k,l),icectra(kpie-1,j,k,l),                        &
     &    icectra(2,j,k,l),icectra(kpie  ,j,k,l)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      


      RETURN
      END
