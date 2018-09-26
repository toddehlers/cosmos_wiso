      SUBROUTINE OCSCHEP
!
!     HORIZONTAL DIFFUSION OF MOMENTUM
!     
!     BY ERNST MAIER-REIMER
!     MODIFIED BY UWE MIKOLAJEWICZ 2/99
!           MAKE AUSF 2D, ADAPT TO SPATIALLY VARIABLE GRID
!           AUS TIME CONSTANT
!
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_UNITS
      IMPLICIT NONE

      INTEGER i,j,k
      REAL P11(IE,JE),P22(IE,JE),P12(IE,JE),AUSF(IE,JE),AUSFPSI(IE,JE)
      REAL DLXXYY, UX, UY, VX, VY

!
!UWE  AUS DIMENSIONLESS MOMENTUM DIFFUSION COEFFICIENT
!
!     AUSF     2D-MOMENTUM DIFFUSION COEFFICENT P-POINTS      [M2/S]
!     AUSFPSI  THE SAME FOR PSI-POINTS
!
!     SPATIAL DEPENDENCE ACCORDING TO BRYAN ET AL. 1975
!
      DO 2712 J=1,JE
        DO 2712 I=1,IE
          P11(I,J)=0.
          P12(I,J)=0.
          P22(I,J)=0.
          UDIV(I,J)=0.
          AUSFPSI(I,J)=0.
          DLXXYY=MAX(DLXP(I,J),DLYP(I,J))
          AUSF(I,J)=AUS                                                 &
     &       *DLXXYY**2
 2712 CONTINUE
!
      DO J=2,JE-1
       DO I=2,IE-1
        AUSFPSI(I,J)=0.25*(AUSF(I,J)+AUSF(I+1,J)                        &
     &                     +AUSF(I,J+1)+AUSF(I+1,J+1))
       ENDDO
      ENDDO
!
      CALL bounds_exch('s',AUSFPSI,'ocschep 1')
!
      DO 1 K=1,KE
      DO 1 J=2,JE1
      DO 1 I=2,IE1
! ZUSATZTERM KANN EVTL. AUF 0 GESETZT WERDEN, DURCHSTROEME!!!!
!
#ifdef SCHEP
      UK1O(I,J,K)=UKO(I,J,K)*AMSUO(I,J,K)                               &
     &      -2.*(UKO(I,J-1,K)+UKO(I,J+1,K))*(1.-AMSUO(I,J,K))
#else
      UK1O(I,J,K)=UKO(I,J,K)*AMSUO(I,J,K)
#endif /*SCHEP*/
#ifdef SCHEP
      VK1E(I,J,K)=VKE(I,J,K)*AMSUE(I,J,K)                               &
     &      -2.*(VKE(I-1,J,K)+VKE(I+1,J,K))*(1.-AMSUE(I,J,K))
#else
      VK1E(I,J,K)=VKE(I,J,K)*AMSUE(I,J,K)
#endif /*SCHEP*/
1     CONTINUE
!
      CALL bounds_exch('u',UK1O,'ocschep 2')
      CALL bounds_exch('v',VK1E,'ocschep 3')
!
      DO 8921 K=1,KE
      DO 8821 J=2,JE1
      DO 8821 I=2,IE1
      UX=(DLYU(I,J)*DDUO(I,J,K)*UKO(I,J,K)                              &
     &     -DLYU(I-1,J)*DDUO(I-1,J,K)*UKO(I-1,J,K))                     &
     &      *DPIO(I,J,K)/(DLXP(I,J)*DLYP(I,J))
!     DV/DY AUF P-PUNKT
      VY=(DDUE(I,J-1,K)*DLXV(I,J-1)*VKE(I,J-1,K)                        &
     &    -DDUE(I,J,K)*DLXV(I,J)*VKE(I,J,K))                            &
     &     *DPIO(I,J,K)/(DLXP(I,J)*DLYP(I,J))  

!     DV/DX AUF PSI-PUNKT
      VX=(VK1E(I+1,J,K)-VK1E(I,J,K))*2./(DLXU(I,J)+DLXU(I,J+1))
!     DU/DY AUF PSI-PUNKT
      UY=(UK1O(I,J,K)-UK1O(I,J+1,K))*2./                                &
     &                (DLXV(I,J)+DLXV(I+1,J))
      P11(I,J)=AUSF(I,J)*DDPO(I,J,K)*(UX-VY)
      P22(I,J)=-P11(I,J)
      P12(I,J)=AUSFPSI(I,J)*DDPSIO(I,J,K)*(UY+VX)
      UDIV(I,J)=UX+VY
8821  CONTINUE

      CALL bounds_exch('p',P11,'ocschep 4')
      CALL bounds_exch('p',P22,'ocschep 5')
      CALL bounds_exch('s',P12,'ocschep 6')
      CALL bounds_exch('p',UDIV,'ocschep 7')
!
      DO 8822 J=2,JE1
      DO 8822 I=2,IE1
!
!replace dwio
!     UKO(I,J,K)=UOO(I,J,K)+                                            &
!    & DWIO(I,J,K)*(DTDXUO(I,J)*(P11(I+1,J)-P11(I,J))                   &
!    & +(DPYO(I,J-1)*P12(I,J-1)-DPYO(I,J)*P12(I,J)))
      UKO(I,J,K)=UOO(I,J,K)+                                            &
     & (amsuo(i,j,k)/(almzer+dduo(i,j,k)))                              &
     &     *(dtdxuo(i,j)*(p11(i+1,j)-p11(i,j))                          &
     & +(dpyo(i,j-1)*p12(i,j-1)-dpyo(i,j)*p12(i,j)))

!
!replace dwie
!     VKE(I,J,K)=VOE(I,J,K)                                             &
!    & + DWIE(I,J,K)*((DTDXUE(I,J)*P12(I,J)-DTDXUE(I-1,J)*P12(I-1,J))   &
!    & + DPYE(I,J)*(P22(I,J)-P22(I,J+1)))
      VKE(I,J,K)=VOE(I,J,K)                                             &
     & + (amsue(i,j,k)/(almzer+ddue(i,j,k)))                            &
     &           *((dtdxue(i,j)*p12(i,j)-dtdxue(i-1,j)*p12(i-1,j))      &
     & + dpye(i,j)*(p22(i,j)-p22(i,j+1)))

8822  CONTINUE
!
      DO 8823 J=2,JE1
      DO 8823 I=2,IE1
      UOO(I,J,K)=UKO(I,J,K)
      VOE(I,J,K)=VKE(I,J,K)
!
8823  CONTINUE
8921  CONTINUE

      CALL bounds_exch('u',UKO,'ocschep 8')
      CALL bounds_exch('v',VKE,'ocschep 9')
      CALL bounds_exch('u',UOO,'ocschep 10')
      CALL bounds_exch('v',VOE,'ocschep 11')

      RETURN
      END
