      SUBROUTINE BELEG
!****************************************************************
!
!     BBBBB   EEEEEE  L       EEEEEE  GGGGG
!     B    B  E       L       E       G
!     BBBBB   EEEEE   L       EEEEE   GGGGG
!     B    B  E       L       E       G    G
!     BBBBB   EEEEEE  LLLLLL  EEEEEE  GGGGGG
!
!
!*****************************************************************
!
!***********************************************************************
!     BELEG : READS AND COMPUTES START VALUES FOR DIFFERENT FIELDS
!             INCLUDING  A.)       ORIENTATION OF GRID POINTS ON EARTH
!                        B.)       DISTANCES BETWEEN POINTS
!                        C.)       LAYER DEPTH CONFIGURATION
!                        D.)       GRID CONFIGURATION
!
!---------------------------------------------------------------
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_UNITS
      IMPLICIT NONE

      REAL DMAXO, DMINO, DELZ,  RHO
      INTEGER I, J, K, II, JJ
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     LONGITUDE AND LATITUDE OF GRID POINTS
!
!=========LAENGE UND BREITE
!     ODD FIELDS J=1,3,5,...,JTO-1
!
      DMAXO=-1.E20
      DMINO=1.E20
!
      DO 152 J=1,JE1
!
      DO 151 I=1,IE
!
!----------------------------------------------------------------
!     360 DEGREE CROSS OVER CHECK
!
!----------------------------------------------------------------
!
      II = I + p_ioff
      JJ = J + p_joff
      IF(DLXU(I,J).LT.1.)WRITE(IO_STDOUT,*) 'DLXU ',II,JJ,DLXU(I,J)
      IF(DLXV(I,J).LT.1.)WRITE(IO_STDOUT,*) 'DLXV ',II,JJ,DLXV(I,J)
      IF(DLXP(I,J).LT.1.)WRITE(IO_STDOUT,*) 'DLXP ',II,JJ,DLXP(I,J)
      IF(DLYU(I,J).LT.1.)WRITE(IO_STDOUT,*) 'DLYU ',II,JJ,DLYU(I,J)
      IF(DLYV(I,J).LT.1.)WRITE(IO_STDOUT,*) 'DLYV ',II,JJ,DLYV(I,J)
      IF(DLYP(I,J).LT.1.)WRITE(IO_STDOUT,*) 'DLYP ',II,JJ,DLYP(I,J)
!
      DTDXUO(I,J)=DT/DLXU(I,J)
      DTDXPO(I,J)=DT/DLXP(I,J)
!
      DTDYO(I,J)=DT/DLYP(I,J)
      DPYO(I,J)=DT/DLYP(I,J)
!
151   CONTINUE
!
!     COMPUTE DISTANCES BETWEEN ZETA UND VECTOR POINTS
!
!     EVEN FIELDS J=2,4,6,...,JTO
!
      DO 153 I=1,IE
!
      II = I + p_ioff
      JJ = J + p_joff
      DELZ=DLXV(I,J)
      IF(DELZ.LE.0.)WRITE(IO_STDOUT,*) 'NULL ',II,JJ
!
      DMAXO=MAX(DMAXO,DLXV(I,J))
      DMINO=MIN(DMINO,DLXV(I,J))
      DTDXUE(I,J) = DT / DLXV(I,J)
      DTDXPE(I,J) = DT / DLXV(I,J)
      DPYE(I,J)  = DT/DLYV(I,J)
!
153   CONTINUE
 152  CONTINUE
!
      CALL global_min(DMINO)
      CALL global_max(DMAXO)
      WRITE(IO_STDOUT,15290)DMINO/1000.,DMAXO/1000.
15290 FORMAT(' RESOLUTION KM  MIN. : ',F10.3,'  MAX. : ',F10.3,' ODD ')
      WRITE(IO_STDOUT,*) (DLXP(5,J),J=1,JE)
!
!-----------------------------------------------------------------------
!
!     LAYER DEPTH CONFIGURATION  (IS ACTUALLY DECOFTED IN MPGR PROGRAM
!                              BY : DATA DZW/..../ WHICH GIVES THE
!                              LAYER THICKNESSES STARTING FROM THE TOP)
!
!                               TIESTU(K) : DEPTH OF U-POINT IN LAYER K
!                               TIESTW(K) : DEPTH OF W-POINT OF THE
!                                           UPPER BOUNDARY OF LAYER K
      TIESTU(1) = 0.5 * DZW(1)
      TIESTW(1) = 0.0
!
      DO 44 K=1,KE
      DWI(K)       = 1. / DZW(K)
      TIESTW(K+1)  = TIESTW(K) + DZW(K)
   44 CONTINUE
!
!-------------------------
!     CALCULATION OF VECTOR POINT DISTANCES DZ(1,..,KEP)
!
      DO 4544 K=2,KE
      TIESTU(K) = 0.5 * ( TIESTW(K+1) + TIESTW(K) )
4544  CONTINUE
       TIESTU(KEP)=9000.
!
      DZ(1) = TIESTU(1)
      DO 4545 K=2,KEP
4545  DZ(K) = TIESTU(K) - TIESTU(K-1)
!
      DO 4504 K=1,KE
      WRITE(IO_STDOUT,4546)K,TIESTW(K),DZW(K)
      WRITE(IO_STDOUT,4547)K,TIESTU(K),DZ (K)
4504  CONTINUE
!
      WRITE(IO_STDOUT,4546)KEP,TIESTW(K),ZERO
      WRITE(IO_STDOUT,4547)KEP,TIESTU(K),DZ (K)
!
4546  FORMAT(' LAYER ',I2,' W-POINT DEPTH ',F6.0,30X,' THICKNESS : ',   &
     &F6.1)
4547  FORMAT(' LAYER ',I2,' U-POINT DEPTH ',15X,F6.0,' DISTANCE  : ',   &
     &F6.1)
!
!--------------------------- END OF LAYER CONFIGURATION ----------------
!
      WRITE(IO_STDOUT,*)' REFERENCE STRATIFICATION : '
!
      DO 701 K=1,KE
!
         DI(K)    = 1. / DZ(K)
         PREFF(K) = 0.1026 * TIESTU(K)
         ROREF(K) = RHO ( SAF(K),TAF(K),PREFF(K) )
!
      WRITE(IO_STDOUT,70199)K,ROREF(K),SAF(K),TAF(K),PREFF(K)
70199 FORMAT(' LAYER ',I3,' RHO : ',F12.4,' S,T,P ',3F10.5)
!
      WRITE(IO_STDOUT,*)IE_G,JE_G,KE
      DO 701 J=1,JE
      DO 701 I=1,IE
!
      STABIO(I,J,K)=ZERO
      
      DDUE(I,J,K)=ZERO
      DDUO(I,J,K)=ZERO
!
      DDPO(I,J,K)=ZERO
      DPIO(I,J,K)=ZERO
!
      WETO(I,J,K)  = ZERO
      AMSUE(I,J,K) = ZERO
      AMSUO(I,J,K) = ZERO
!
      PO(I,J,K)=ZERO
!
      DVO(I,J,K)=ZERO
!
      VKE(I,J,K) = ZERO
      UKO(I,J,K) = ZERO
!
      UK1O(I,J,K)=ZERO
!
      VK1E(I,J,K)=ZERO
!
      WO(I,J,K)=ZERO
!
      T1O(I,J,K)=ZERO
!
      S1O(I,J,K)=ZERO
!
!......................................
!     TEMPERATURE AND SALINITY FIELDS  SET TO REFERENCE STRATIFICATION
!                 TO START FROM A DIFFERENT STRATIFICATION SPECIFY
!                 SBR STATUS0 AND CALL STATUS0
!......................................
!
      THO(I,J,K)=TAF(K)
!
      SAO(I,J,K)=SAF(K)
!
701   CONTINUE
!
      DO 702 I=1,ILT
!
      B(I)=ZERO
      X(I)=ZERO
702   CONTINUE
!
! FOR LAYER KEP   :
!
      DI(KEP)=1./DZ(KE)
!
      DO 344 J=1,JE
      DO 344 I=1,IE
!
      WO(I,J,KEP)=ZERO
      DVO(I,J,KEP)=ZERO
!
344    CONTINUE
!
      DO 7896 J=1,JE
      DO 7896 I=1,IE
!
      DVO(I,J,1)=ZERO
!
      HEATO(I,J)=ZERO
      EMINPO(I,J)=ZERO
!
      V1E(I,J)=ZERO
!
      U1O(I,J)=ZERO
!
      Z1O(I,J)=ZERO
!
      ZO(I,J)=ZERO
!
7896  CONTINUE
!
      RETURN
      END
