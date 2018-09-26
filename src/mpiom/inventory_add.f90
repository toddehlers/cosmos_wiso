      SUBROUTINE INVENTORY_ADD(kpie,kpje,kpke)
!*******************************************************************
!
!**** *INVENTORY_ADD* - calculate the ADD inventory.
!
!    
!     
!     Purpose
!     -------
!     - calculate the ADD inventory.
!
!     Method
!     -------
!     -
!
!**   Interface.
!     ----------
!
!     *CALL*       *INVENTORY_ADD*
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!
!     Externals
!     ---------
!     GLOBAL_SUM (../src_oce/mo_parallel.f90)
!
!**********************************************************************
    
      USE mo_contra
      USE mo_control_add
      USE mo_addmean
      use mo_param1_add 
      
      USE MO_COMMO1
      USE MO_COMMOAU1

      USE mo_parallel

      implicit none

      INTEGER :: kpie,kpje,kpke,i,j,k,l

     
      REAL :: zocectrato(nocectra),zicectrato(nicectra)
  
      REAL :: ztotvol,ztotarea,vol
  
  
      

!
!  oceanic tracers
!----------------------------------------------------------------------

      ztotvol=0.
      k=1
      DO j=2,kpje-1
      DO i=2,kpie-1
         ztotvol=ztotvol+WETO(I,J,K)*dlxp(i,j)*dlyp(i,j)*(DDPO(i,j,k) &
     &                     +ZO(I,J)-SICTHO(I,J)*RHOICWA               &
     &                     -SICSNO(I,J)*RHOSNWA)	 
      ENDDO
      ENDDO     
      
      
      DO k=2,kpke
      DO j=2,kpje-1
      DO i=2,kpie-1
         ztotvol=ztotvol+WETO(I,J,K)*dlxp(i,j)*dlyp(i,j)*DDPO(i,j,k) 
      ENDDO
      ENDDO
      ENDDO

      CALL global_sum(ztotvol)

      DO l=1,nocectra
        zocectrato(l) = 0.0
        k=1
        DO j=2,kpje-1
        DO i=2,kpie-1
          vol    = dlxp(i,j)*dlyp(i,j)*(DDPO(i,j,k)                   &
     &                   +ZO(I,J)-SICTHO(I,J)*RHOICWA                 &
     &                   -SICSNO(I,J)*RHOSNWA)
          zocectrato(l)=zocectrato(l)+WETO(I,J,K)*ocectra(i,j,k,l)*vol
        ENDDO
        ENDDO      
      
        DO k=2,kpke
        DO j=2,kpje-1
        DO i=2,kpie-1
          vol    = dlxp(i,j)*dlyp(i,j)*DDPO(i,j,k)
          zocectrato(l)=zocectrato(l)+WETO(I,J,K)*ocectra(i,j,k,l)*vol
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      CALL global_sum(zocectrato)
    
      DO l=1,nicectra
        zicectrato(l) = 0.0
        k=1
        DO j=2,kpje-1
        DO i=2,kpie-1
          vol    = dlxp(i,j)*dlyp(i,j)*(DDPO(i,j,k)                   &
     &                   +ZO(I,J)-SICTHO(I,J)*RHOICWA                 &
     &                   -SICSNO(I,J)*RHOSNWA)
          zicectrato(l)=zicectrato(l)+WETO(I,J,K)*icectra(i,j,k,l)*vol
        ENDDO
        ENDDO      
      
        DO k=2,kpke
        DO j=2,kpje-1
        DO i=2,kpie-1
          vol    = dlxp(i,j)*dlyp(i,j)*DDPO(i,j,k)
          zicectrato(l)=zicectrato(l)+WETO(I,J,K)*icectra(i,j,k,l)*vol
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      CALL global_sum(zicectrato)

      WRITE(io_stdo_add,*) ' '
      WRITE(io_stdo_add,*) 'Global inventory of ocean additional coonservative tracers'
      WRITE(io_stdo_add,*) '------------------------------------------'
      WRITE(io_stdo_add,*) ' '
      WRITE(io_stdo_add,*) '       total[kmol]  concentration[kmol/m^3]'
      WRITE(io_stdo_add,*) ' '
      DO l=1,nocectra
      WRITE(io_stdo_add,*) 'No. ',l,' ocectra ',zocectrato(l),zocectrato(l)/ztotvol
      WRITE(io_stdo_add,*) 'No. ',l,' icectra ',zicectrato(l),zicectrato(l)/ztotvol
      WRITE(io_stdo_add,*) 'No. ',l,' ocectra+icectra ',(zocectrato(l)+zicectrato(l)),(zocectrato(l)+zicectrato(l))/ztotvol
      ENDDO
      
   

      RETURN
      END
