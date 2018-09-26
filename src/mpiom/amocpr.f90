      SUBROUTINE AMOCPR
#ifdef AMOCEMR
      USE MO_PARAM1
      USE MO_COMMO1, ONLY: GIPH_G, GILA_G
      USE MO_PARALLEL
      USE MO_UNITS
      PARAMETER(NSOSTA=6,IPRO=100)
!
       COMMON/AMOC/ A(NSOSTA,NSOSTA,3,3),                               &
     & BX(NSOSTA,NSOSTA,3,3),VEC(NSOSTA,NSOSTA,3)                       &
     & ,C(NSOSTA,NSOSTA,3,3)                                            &
     & ,ISOSTA(NSOSTA),JSOSTA(NSOSTA),ARCL(NSOSTA,NSOSTA)               &
     & ,XSP(IE_G,JE_G),YSP(IE_G,JE_G),ZSP(IE_G,JE_G)                    &
     & ,IRAY(IPRO,NSOSTA,NSOSTA),JRAY(IPRO,NSOSTA,NSOSTA),              &
     & BERAY(IPRO,NSOSTA,NSOSTA),DERAY(IPRO,NSOSTA,NSOSTA)              &
     &, TRATI(NSOSTA,NSOSTA),SEGMA(NSOSTA,NSOSTA)                       &
     & ,XRAY(IPRO,NSOSTA,NSOSTA),YRAY(IPRO,NSOSTA,NSOSTA)
!C       DATA ISOSTA/73,100,99,129,120,140/
!C       DATA JSOSTA/34,12,63,19,68,26/
       DATA ISOSTA/89,74,87,85,63,62/
       DATA JSOSTA/ 8,26,32,58,36,70/
       DO 1 J=1,JE_G
       DO 1 I=1,IE_G
       ZSP(I,J)=SIN(GIPH_G(2*I,2*J))
       XSP(I,J)=COS(GIPH_G(2*I,2*J))*COS(GILA_G(2*I,2*J))
       YSP(I,J)=COS(GIPH_G(2*I,2*J))*SIN(GILA_G(2*I,2*J))
1      CONTINUE
       IF (ICYCLI.EQ.1) THEN
         XSP(1,:) = XSP(IE_G-1,:)
         YSP(1,:) = YSP(IE_G-1,:)
         ZSP(1,:) = ZSP(IE_G-1,:)
         XSP(IE_G,:) = XSP(2,:)
         YSP(IE_G,:) = YSP(2,:)
         ZSP(IE_G,:) = ZSP(2,:)
       ENDIF
!
       DO 3000 I=2,NSOSTA
       DO 3000 J=1,I-1
!
!      PHI=90.+GRIDSQ-FLOAT(J)*GRIDSH
!      LAMBDA=I*GRIDS+PHI-90
       PH1=GIPH_G(2*ISOSTA(I),2*JSOSTA(I))
       PH2=GIPH_G(2*ISOSTA(J),2*JSOSTA(J))
       AL1=GILA_G(2*ISOSTA(I),2*JSOSTA(I))
       AL2=GILA_G(2*ISOSTA(J),2*JSOSTA(J))
!
!       WRITE(IO_STDOUT,*) ' GRIX ??? ', GRIX
!       WRITE(IO_STDOUT,*) ' GRADE ',AAL1,APH1,AAL2,APH2
       X1=COS(PH1)*COS(AL1)
       Y1=COS(PH1)*SIN(AL1)
       Z1=SIN(PH1)
       X2=COS(PH2)*COS(AL2)
       Y2=COS(PH2)*SIN(AL2)
       Z2=SIN(PH2)
!
       XA=Y1*Z2-Y2*Z1
       YA=Z1*X2-Z2*X1
       ZA=X1*Y2-X2*Y1
       WRITE(IO_STDOUT,*) ' X1 ', X1,Y1,Z1
       WRITE(IO_STDOUT,*) ' X2', X2,Y2,Z2
       BETR=SQRT(XA**2+YA**2+ZA**2)
       XA=XA/BETR
       YA=YA/BETR
       ZA=ZA/BETR
       XR=Y1*ZA-YA*Z1
       YR=Z1*XA-ZA*X1
       ZR=X1*YA-XA*Y1
!
       WRITE(IO_STDOUT,*)' ACHS ', XA,YA,ZA
       WRITE(IO_STDOUT,*) 'XR ',XR,YR,ZR
       S1=X1**2 + Y1**2 + Z1**2
       S2=X2**2 + Y2**2 + Z2**2
       S3=XA**2 + YA**2 + ZA**2
       S4=XR**2 + YR**2 + ZR**2
       S5=X1*XA + Y1*YA + Z1*ZA
       S6=X1*XR + Y1*YR + Z1*ZR
       S7=XA*XR + YA*YR + ZA*ZR
       WRITE(IO_STDOUT,*)S1,S2,S3,S4,S5,S6,S7
       A(I,J,1,1)=X1
       A(I,J,2,1)=XR
       A(I,J,3,1)=XA
       A(I,J,1,2)=Y1
       A(I,J,2,2)=YR
       A(I,J,3,2)=YA
       A(I,J,1,3)=Z1
       A(I,J,2,3)=ZR
       A(I,J,3,3)=ZA
!       WRITE(IO_STDOUT,*) ' ROTATION MATRIX '
!       WRITE(IO_STDOUT,600) A
!       WRITE(IO_STDOUT,*) ' VECTORBASIS '
!       WRITE(IO_STDOUT,600)X1,Y1,Z1,X2,Y2,Z2,XA,YA,ZA
       DEPH=ACOS(X1*X2+Y1*Y2+Z1*Z2)
!
        WRITE(IO_STDOUT,*) ' DEPH ', DEPH
        DDPH=DEPH*0.01
       SEGMA(I,J)=RADIUS*DDPH
       DO 3000 JJJ=1,IPRO
!
       WIN=DDPH*(FLOAT(JJJ)-0.5)
!
       BX(I,J,1,1)=COS(WIN)
       BX(I,J,2,2)=BX(I,J,1,1)
       BX(I,J,2,1)=-SIN(WIN)
       BX(I,J,1,2)=-BX(I,J,2,1)
       BX(I,J,3,3)=1.
       DO 2 M=1,3
       VEC(I,J,M)=0.
       DO 2 K=1,3
       C(I,J,K,M)=0.
       DO 3 L=1,3
       C(I,J,K,M)=C(I,J,K,M)+A(I,J,L,M)*BX(I,J,L,K)
!
3      CONTINUE
2      CONTINUE
!       WRITE(IO_STDOUT,600) (C(I,J,1,L),L=1,3)
!
       XGL=C(I,J,1,1)
       YGL=C(I,J,1,2)
       ZGL=C(I,J,1,3)    
       AL=ATAN2(YGL,XGL)
       PH=ASIN(ZGL)
       DISS=15.
       DO 1317 JJ=1,JE_G-1
       DO 1317 II=1,IE_G-1
       IF(ZSP(II,JJ).LT.0.85) GO TO 1317
       YM=0.25*(YSP(II,JJ)+YSP(II+1,JJ)+YSP(II,JJ+1)+YSP(II+1,JJ+1))
       XM=0.25*(XSP(II,JJ)+XSP(II+1,JJ)+XSP(II,JJ+1)+XSP(II+1,JJ+1))
       DIF=ABS(XGL-XM)+ABS(YGL-YM)
       IF(DIF.GT.DISS) GO TO 1317
!
       DISS=DIF
       XT=XGL-XSP(II,JJ)
       YT=YGL-YSP(II,JJ)
       BETR=SQRT(XT**2+YT**2)
       BE=1.
       DE=1.
!
       X1=XSP(II+1,JJ)-XSP(II,JJ)
       Y1=YSP(II+1,JJ)-YSP(II,JJ)
       X2=XSP(II,JJ+1)-XSP(II,JJ)
       Y2=YSP(II,JJ+1)-YSP(II,JJ)
!
        DETI=1./(X1*Y2-X2*Y1)
        AL=DETI*(Y2*XT-X2*YT)
        BE=DETI*(-Y1*XT+X1*YT)
!
       IRAY(JJJ,I,J)=II
       JRAY(JJJ,I,J)=JJ
       BERAY(JJJ,I,J)=1.-AL
       DERAY(JJJ,I,J)=1.-BE    
       XRAY(JJJ,I,J)=II+AL
       YRAY(JJJ,I,J)=JJ+BE
1317   CONTINUE
!       WRITE(IO_STDOUT,*)' WODENN?',I,J
!     1 ,IRAY(JJJ,I,J),JRAY(JJJ,I,J)
!     1 ,XRAY(JJJ,I,J),YRAY(JJJ,I,J)
!       WRITE(IO_STDOUT,*)' ABER ',XGL,YGL,TISS
!       WRITE(IO_STDOUT,*) 'WIEDENN ',X1,Y1,XT,YT
3000   CONTINUE
       IF(p_pe==p_io) THEN
         OPEN(IO_OU_RAY,FILE='RAY', FORM='UNFORMATTED')
         WRITE(IO_OU_RAY)XRAY,YRAY
         CLOSE(IO_OU_RAY)
       ENDIF
#endif /*AMOCEMR*/
       RETURN
       END
