      SUBROUTINE CALCGMVIS

      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMOAU1
      USE MO_UNITS
!
!-----------------------------------------------------------------------
!
!       CALCULATION OF RICHARDSON NUMBER DEPENDENT
!       COEFFICIENTS AFTER VISBECK ET AL. JPO 27, 1997
!
!
!         RINUM : ( (G/RHO)*D(RHO)/DZ ) / ( (D(VELOCITY)/DZ)**2 )
!                  RICHARDSON NUMBER (EVEN,ODD ==> RINUME,RINUMO)
!         GOR   : G/RHO
!
!
!---------------------------------------------------------------------
!
      DIMENSION VISRICH_M(IE,JE),VISROSSBY_M(IE,JE)
!=======================================================================
!
      VISAL=0.005
      VISMAX= 2.0*CAH00 
      VISMIN=25.0
!:: SET UPPER LIMIT FOR DEPTH INTEGRATED RICHARDSON NUMBER
      DEPUP=20.0
!:: SET LOWER LIMIT FOR DEPTH INTEGRATED RICHARDSON NUMBER
      DEPDOWN=1000.0
       KVISD=1
       KVISU=1
      DO k=1,ke
       IF (DEPUP .GE. TIESTU(K)) KVISU= K
       IF (DEPDOWN .GE. TIESTU(K)) KVISD= K
      ENDDO
!
!    
      WRITE(IO_STDOUT,*) ' CALCULATING VISBECK COEFFICIENT BETWEEN'
      WRITE(IO_STDOUT,*) ' K= ',KVISU,' AND K= ',KVISD              
      GOR=G/RHOWAT
!
!--------------------------------------------------------------------
      DO I=1,IE
       DO J=1,JE
        VISRICH_M(I,J)=0.0
        VISROSSBY_M(I,J)=0.0
        THELP(I,J)=0.0
        SHELP(I,J)=0.0
        BOLX(I,J)=MAX(BOLX(I,J),VISMIN)
        BOLY(I,J)=MAX(BOLY(I,J),VISMIN)
       ENDDO
      ENDDO
!
! CALCULATE VERTICALLY AVERAGED RICHARDSEN NR. FOR
! CALCULATION OF VISBECK COEFFICIENTS
! AFTER GENT ET AL. J. CLIM. 2002
!
      DO J=2,JE1
       DO I=2,IE-1
        KDOWN=MIN(KBOT(I,J),KVISD)
        IF (KDOWN .GT. KVISU) THEN
        DO K=KVISU,KDOWN 
         KO=MAX(K-1,1)
!
         DUDO=ALMZER + WETO(I,J,K) * DI(K)**2                           &
     & * (   ( UKO(I-1,J,K) - UKO(I-1,J,KO) )**2                        &
     &     + ( UKO(I,J,K)   - UKO(I,J,KO)   )**2                        &
     &     + ( VKE(I,J-1,K) - VKE(I,J-1,KO) )**2                        &
     &     + ( VKE(I,J,K)   - VKE(I,J,KO)   )**2   ) *0.5
!
         HHO=WETO(I,J,K)*(AMSUO(I,J,K)+AMSUO(I-1,J,K)+AMSUE(I,J-1,K)    &
     & +AMSUE(I,J,K))*0.25
!
         RINUMO=HHO*MAX(GOR*STABIO(I,J,K)/DUDO,0.)
! INTEGRATE VERTICALLY
           VISRICH_M(I,J)=VISRICH_M(I,J)+RINUMO*DZ(K)
           VISN=SQRT(HHO*MAX(GOR*STABIO(I,J,K),0.))
!:: USE THELP AS TEMP FIELD TO SUM N
           THELP(I,J)=THELP(I,J)+VISN*DZ(K)        
!:: USE SHELP AS TEMP FIELD TO SUM EFFECTIVE DEPTH RANGE
           SHELP(I,J)=SHELP(I,J)+DZ(K)
!
        ENDDO
       ENDIF
!
       ENDDO
      ENDDO
!
!OtB_GMVISETAL_TEST ifdef GMVISETAL_TEST
!OtB_GMVISETAL_TEST      ii1=ldt
!OtB_GMVISETAL_TEST      ii4=ie*je
!OtB_GMVISETAL_TEST      ii3=-100
!OtB_GMVISETAL_TEST      ii2=251
!OtB_GMVISETAL_TEST      WRITE(161) ii1,ii2,ii3,ii4
!OtB_GMVISETAL_TEST      WRITE(161) ((VISRICH_M(I,J)/(SHELP(I,J)+ALMZER),I=1,IE),J=1,JE)
!OtB_GMVISETAL_TEST      ii2=252
!OtB_GMVISETAL_TEST      WRITE(161) ii1,ii2,ii3,ii4
!OtB_GMVISETAL_TEST      WRITE(161) ((THELP(I,J)/(SHELP(I,J)+ALMZER),I=1,IE),J=1,JE)
!OtB_GMVISETAL_TEST      ii2=253
!OtB_GMVISETAL_TEST      WRITE(161) ii1,ii2,ii3,ii4
!OtB_GMVISETAL_TEST      WRITE(161) ((SHELP(I,J),I=1,IE),J=1,JE)
!OtB_GMVISETAL_TEST      igm=56
!OtB_GMVISETAL_TEST      jgm=24
!OtB_GMVISETAL_TEST      WRITE(*,*) 'gmvis at 56,24'
!OtB_GMVISETAL_TEST      WRITE(*,*) 'depto',depto(igm,jgm)
!OtB_GMVISETAL_TEST      WRITE(*,*) visrich_m(igm,jgm),thelp(igm,jgm),shelp(igm,jgm)
!OtB_GMVISETAL_TEST endif /*GMVISETAL_TEST*/
!
!:: CALCULATE VISBECK ET AL. DIFFUSION COEFFICIENTS FOR GM
        DO I=2,IE1
         DO J=2,JE1
          IF (WETO(I,J,1) .GT. 0.5) THEN
!:: CALCULATE ROSSBY RADIUS
          CORIO=ABS(0.25*(FTWOU(I,J)+FTWOU(I-1,J)+                      &
     &                   FTWOV(I,J)+FTWOV(I,J-1)))
!:: LIMIT F TO 2.5 degree off equator
          CORIO=MAX(2.5e-6,CORIO)
          VISROSSBY_M(I,J)=THELP(I,J)/CORIO
          VISROSSBY_M(I,J)=MIN(2.5e6,VISROSSBY_M(I,J))
!:: CALCULATE RICHARDSON_NUMBER
          RICHI=VISRICH_M(I,J)/(SHELP(I,J)+ALMZER)
          VISLX=MAX(VISROSSBY_M(I,J),DLXP(I,J))
          VISLY=MAX(VISROSSBY_M(I,J),DLYP(I,J))
          VISFAC=CORIO*VISAL/SQRT(RICHI+ALMZER)
          IF (SHELP(I,J) .LT. ALMZER) VISFAC=0.0
! STORE ROSSBY NR. FOR DIAGNOSTICS
          VISROSSBY_M(I,J)=VISLX                        
          THELP(I,J)=VISLX**2*VISFAC             
          SHELP(I,J)=VISLY**2*VISFAC             
          THELP(I,J)=MAX(VISMIN,THELP(I,J))
          THELP(I,J)=MIN(VISMAX,THELP(I,J))
          SHELP(I,J)=MAX(VISMIN,SHELP(I,J))
          SHELP(I,J)=MIN(VISMAX,SHELP(I,J))
          ELSE
           THELP(I,J)=VISMIN
           SHELP(I,J)=VISMIN
          ENDIF
        ENDDO
       ENDDO
!
       CALL bounds_exch('p',THELP,'calcgmvis 10')
       CALL bounds_exch('p',SHELP,'calcgmvis 10')
       CALL bounds_exch('p',VISRICH_M,'calcgmvis 10')
       CALL bounds_exch('p',VISROSSBY_M,'calcgmvis 10')
!
! TAKE 0.5 times old day bolx/boly
       DO I=1,IE1
        DO J=1,JE1
         BOLX(I,J)=0.5*BOLX(I,J)+0.25*(THELP(I,J)+THELP(I+1,J))
         BOLY(I,J)=0.5*BOLY(I,J)+0.25*(SHELP(I,J)+SHELP(I,J+1))
        ENDDO
       ENDDO

       CALL bounds_exch('u',BOLX,'calcgmvis 10')
       CALL bounds_exch('v',BOLY,'calcgmvis 10')

!OtB_GMVISETAL_TEST ifdef GMVISETAL_TEST
!OtB_GMVISETAL_TEST      ii1=ldt
!OtB_GMVISETAL_TEST      ii4=ie*je
!OtB_GMVISETAL_TEST      ii3=-100
!OtB_GMVISETAL_TEST      ii2=254
!OtB_GMVISETAL_TEST      WRITE(161) ii1,ii2,ii3,ii4
!OtB_GMVISETAL_TEST      WRITE(161) ((THELP(I,J),I=1,IE),J=1,JE)
!OtB_GMVISETAL_TEST      ii2=255
!OtB_GMVISETAL_TEST      WRITE(161) ii1,ii2,ii3,ii4
!OtB_GMVISETAL_TEST      WRITE(161) ((SHELP(I,J),I=1,IE),J=1,JE)
!OtB_GMVISETAL_TEST      ii2=256
!OtB_GMVISETAL_TEST      WRITE(161) ii1,ii2,ii3,ii4
!OtB_GMVISETAL_TEST      WRITE(161) ((VISROSSBY_M(I,J),I=1,IE),J=1,JE)
!OtB_GMVISETAL_TEST endif /*GMVISETAL_TEST*/
!
!=====================================================================

      RETURN
      END
