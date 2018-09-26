      SUBROUTINE maschk_add(kpie,kpje,kpke,imas)
!*******************************************************************
!
!**** *INVENTORY_ADD* - calculate the ADD inventory.
!
!    
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
!     none.
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

      INTEGER :: kpie,kpje,kpke,i,j,k,l,imas
      REAL :: zocectrato(nocectra),zicectrato(nicectra)
      REAL :: ztotvol,ztotarea,vol

      return
      i=25
      j=14
      k=5
      l=1
      write(0,*)'mas',imas,ocectra(i,j,k+1,l),ocectra(i,j,k,l)
   

 
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
    
      WRITE(io_stdo_add,*) ' '
      WRITE(io_stdo_add,*) 'Global inventory of advected ocean tracers'
      WRITE(io_stdo_add,*) '------------------------------------------'
      WRITE(io_stdo_add,*) ' '
      WRITE(io_stdo_add,*) '       total[kg]  concentration[kg/m^3]'
      WRITE(io_stdo_add,*) ' '
      DO l=1,nocectra
      WRITE(io_stdo_add,*) 'No. ',l,    &
     &         zocectrato(l),zocectrato(l)/ztotvol
      ENDDO
       

      RETURN
      END
