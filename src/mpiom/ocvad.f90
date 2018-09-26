SUBROUTINE OCVAD(TRF)
  !OtB_OLDMOMADV ifndef OLDMOMADV
  USE MO_PARAM1
  USE MO_PARALLEL
  USE MO_COMMO1
  USE MO_UNITS
  !
  REAL TRP(IE,JE,KEP),TRM(IE,JE,0:KE)
  REAL TRF(IE,JE,KE),TRVOL(IE,JE,KE)
  REAL TR1(IE,JE,KE)
  REAL WTP(IE,JE,KEP),WTM(IE,JE,0:KE)
  
!$OMP PARALLEL PRIVATE(i,j,k,klo,kup,abl,zwabl,wup,wlo,up,um,scha,sch,schi,rn,uos,uwe,vsu,vno)
!
!$OMP DO
  DO K=1,KE
    DO J=1,JE
      DO I=1,IE
        TRVOL(I,J,K)=DLXV(I,J)*DLYV(I,J)*DDUE(I,J,K)
        TRP(I,J,K)=0.
        TRM(I,J,K)=0.
        TR1(I,J,K)=TRF(I,J,K)
        WTP(I,J,K)=0.
        WTM(I,J,K)=0.
      ENDDO
    ENDDO
  ENDDO
  !
!$OMP DO
  DO J=1,JE
    DO I=1,IE
      TRP(I,J,KEP)=0.
      TRM(I,J,0)=0.
      WTP(I,J,KEP)=0.
      WTM(I,J,0)=0.
    ENDDO
  ENDDO
  !
!$OMP DO
  DO K=1,KE
    KLO=MIN(K+1,KE)
    KUP=MAX(K-1,1)
    DO J=2,JE1
      DO I=2,IE1
        ABL=ABS(TRF(I,J,KUP)-TRF(I,J,KLO))
        ZWABL=ABS(TRF(I,J,KUP)+TRF(I,J,KLO)-2.*TRF(I,J,K))
        WUP=0.5*(WO(I,J,K)+WO(I,J+1,K))
        WLO=0.5*(WO(I,J,K+1)+WO(I,J+1,K+1))
        UP=0.5*DT*DLYV(I,J)*DLXV(I,J)*(WUP+ABS(WUP))
        UM=0.5*DT*DLYV(I,J)*DLXV(I,J)*(ABS(WLO)-WLO)
        !
        SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))
        !
        SCH=MIN(1.,SCHA*TRVOL(I,J,K)/(UP+UM+1.E-20))
        SCHI=1.-SCH
        TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J,KUP)))
        TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J,KLO)))
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      ENDDO
    ENDDO
  ENDDO
  !
!$OMP DO
  DO J=2,JE1
    DO I=1,IE1
      TRP(I,J,KEP)=0.
      WTP(I,J,KEP)=0.
      TRM(I,J,0)=0.
      WTM(I,J,0)=0.
      WTP(I,J,1)=0.
    ENDDO
  ENDDO
  !
  ! RJ: no boundary exchange necessary here, since TRP, WTP, TRM, WTM
  !     are used only in the inner domain in the following loop
!$OMP DO
  DO K=1,KE
    DO J=2,JE1
      DO I=2,IE1
        IF(AMSUE(I,J,K).GT.0.5) THEN
          !RJ        CONT=TRVOL(I,J,K)*TRF(I,J,K)
          !
          RN=TRVOL(I,J,K)+WTP(I,J,K+1)-WTP(I,J,K)-WTM(I,J,K)+WTM(I,J,K-1)
          TR1(I,J,K)=(TRF(I,J,K)*TRVOL(I,J,K)+TRP(I,J,K+1)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I,J,K-1))/RN
          TRVOL(I,J,K)=RN
        ENDIF
      ENDDO
    ENDDO
    DO J=2,JE1
      DO I=2,IE1
        TRF(I,J,K)=TR1(I,J,K)
      ENDDO
    ENDDO
  ENDDO
  !
!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('v',TR1,'ocvad 1')
  CALL bounds_exch('v+',TRVOL,'ocvad 2')
  CALL bounds_exch('v',TRF,'ocvad 3')
!$OMP END SINGLE
!#endif
  !
!$OMP DO
  DO K=1,KE
    DO J=2,JE1
      DO I=2,IE1
        ABL=ABS(TRF(I+1,J,K)-TRF(I-1,J,K))
        ZWABL=ABS(TRF(I+1,J,K)+TRF(I-1,J,K)-2.*TRF(I,J,K))
        UOS=0.5*(UKO(I,J,K)*DDUO(I,J,K)+UKO(I,J+1,K)*DDUO(I,J+1,K))
!EMR    UWE=0.5*(UKO(I-1,J,K)*DDUO(I,J,K)+UKO(I,J+1,K)*DDUO(I,J+1,K))
        UWE=0.5*(UKO(I-1,J,K)*DDUO(I-1,J,K)+UKO(I-1,J+1,K)*DDUO(I-1,J+1,K))
        UP=0.5*DT*DLYpsi(I,J)*(UOS+ABS(UOS))
        UM=0.5*DT*(ABS(UWE)-UWE)*DLYpsi(I-1,J)
        SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))
        !
        SCH=MIN(1.,SCHA*TRVOL(I,J,K)/(UP+UM+1.E-20))
        SCHI=1.-SCH
        !
        TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I+1,J,K)))
        TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I-1,J,K)))
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      ENDDO
    ENDDO
  ENDDO
  !
!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('uu',TRP,'ocvad 4')
  CALL bounds_exch('uu',TRM,'ocvad 5')
  CALL bounds_exch('uu',WTP,'ocvad 6')
  CALL bounds_exch('uu',WTM,'ocvad 7')
!$OMP END SINGLE
!#endif
  !
!$OMP DO
  DO K=1,KE
    DO J=2,JE1
      DO I=2,IE1
        IF(AMSUE(I,J,K).GT.0.5) THEN
          !
          RN=TRVOL(I,J,K)+WTP(I-1,J,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I+1,J,K)
          TR1(I,J,K)=(TRF(I,J,K)*TRVOL(I,J,K)+TRP(I-1,J,K)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I+1,J,K))/RN
          TRVOL(I,J,K)=RN
        ENDIF
      ENDDO
    ENDDO
    DO J=2,JE1
      DO I=2,IE1
        TRF(I,J,K)=TR1(I,J,K)
      ENDDO
    ENDDO
  ENDDO
  !
!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('v',TR1,'ocvad 8')
  CALL bounds_exch('v+',TRVOL,'ocvad 9')
  CALL bounds_exch('v',TRF,'ocvad 10')
!$OMP END SINGLE
!#endif
  !
!$OMP DO
  DO K=1,KE
    DO J=2,JE1
      DO I=2,IE1
        ABL=ABS(TRF(I,J-1,K)-TRF(I,J+1,K))
        ZWABL=ABS(TRF(I,J-1,K)+TRF(I,J+1,K)-2.*TRF(I,J,K))
        VNO=0.5*(VKE(I,J,K)*DDUE(I,J,K)+VKE(I,J-1,K)*DDUE(I,J-1,K))
        VSU=0.5*(VKE(I,J,K)*DDUE(I,J,K)+VKE(I,J+1,K)*DDUE(I,J+1,K))       
        UP=0.5*DT*DLXp(I,J)*(VNO+ABS(VNO))
        UM=0.5*DT*DLXp(I,J+1)*(ABS(VSU)-VSU)

        !
        SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))
        !
        SCH=MIN(1.,SCHA*TRVOL(I,J,K)/(UP+UM+1.E-20))
        SCHI=1.-SCH
        TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J-1,K)))
        TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J+1,K)))
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      ENDDO
    ENDDO
  ENDDO
  !
!$OMP SINGLE
!#ifdef bounds_exch_save


  CALL bounds_exch('v+',TRP,'ocvad 11')
  CALL bounds_exch('v+',TRM,'ocvad 12')
  CALL bounds_exch('v+',WTP,'ocvad 13')
  CALL bounds_exch('v+',WTM,'ocvad 14')

!#endif

#ifdef bounds_exch_tp
    if(p_joff.eq.0) then
    do k=1,ke
    do i=2,ie1
    tf=trp(i,1,k)
    wf=wtp(i,1,k)
!    wtp(i,1,k)=wtm(i,1,k)
!    trp(i,1,k)=-trm(i,1,k)
    wtm(i,1,k)=wf
    trm(i,1,k)=-tf
    enddo
    enddo
    endif
#endif



!$OMP END SINGLE
  !
!$OMP DO
  DO K=1,KE
    DO J=2,JE1
      DO I=2,IE1
        IF(AMSUE(I,J,K).GT.0.5) THEN
          !
          RN=TRVOL(I,J,K)+WTP(I,J+1,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I,J-1,K)
          TR1(I,J,K)=(TRF(I,J,K)*TRVOL(I,J,K)+TRP(I,J+1,K)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I,J-1,K))/RN
          TRVOL(I,J,K)=RN
        ENDIF
      ENDDO
    ENDDO
    DO J=2,JE1
      DO I=2,IE1
        TRF(I,J,K)=TR1(I,J,K)
      ENDDO
    ENDDO
  ENDDO
  !
!$OMP END PARALLEL
  !OtB_OLDMOMADV endif /*OLDMOMADV*/
!#ifdef bounds_exch_save
  CALL bounds_exch('v',TRF,'ocvad 15')
!#endif
END SUBROUTINE OCVAD
