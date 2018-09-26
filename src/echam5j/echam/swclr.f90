SUBROUTINE SWCLR &
  &( KIDIA , KFDIA , KBDIM  , KLEV  , KAER  ,KNEWAER, KAERH, KNU &
  &, PAER  , PSWEXT, PSWSSA, PSWASY &
  &, PALBP , PDSIG , PRAYL , PSEC &
  &, PCGAZ , PPIZAZ, PRAY1 , PRAY2 , PREFZ , PRJ  &
  &, PRK   , PRMU0 , PTAUAZ, PTRA1 , PTRA2 , PTRCLR &
  &, PDIFFUSE, PTAUAZ_N, PPIZAZ_N, PCGAZ_N) !new set of optical properties

!**** *SWCLR* - CLEAR-SKY COLUMN COMPUTATIONS

!     PURPOSE.
!     --------
!           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
!     CLEAR-SKY COLUMN

!**   INTERFACE.
!     ----------

!          *SWCLR* IS CALLED EITHER FROM *SW1S*
!                                OR FROM *SWNI*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------


!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 94-11-15
!        Modified : 96-03-19 JJM-PhD (loop 107 in absence of aerosols)
!        JJMorcrette 990128 : sunshine duration
!        99-05-25   JJMorcrette    Revised aerosols
!        June 2005: CCagnazzo from 4 to 6 bands, following 26Rcycle Morcrette
!        Code
   
   
!     ------------------------------------------------------------------


USE MO_KIND       , ONLY : DP

USE MO_AERO_TANRE , ONLY : RTAUA =>TAUA  ,RPIZA =>PIZA  ,RCGA =>CGA
USE MO_AERO_GADS  , ONLY : RTAUAN=>TAUAN ,RPIZAN=>PIZAN ,RCGAN=>CGAN
USE MO_RADIATION  , ONLY : NDFAER
USE MO_SW         , ONLY : NOVLP ,NSW
USE MO_CONTROL    , ONLY : L_VOLC

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER :: KAER
INTEGER :: KNEWAER
INTEGER :: KFDIA
INTEGER :: KIDIA
INTEGER :: KLEV
INTEGER :: KBDIM
INTEGER :: KNU



!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL(DP):: PAER(KBDIM,KLEV,KAER+KNEWAER), PALBP(KBDIM,NSW)&
  &,  PDSIG(KBDIM,KLEV)&
  &,  PRAYL(KBDIM)&
  &,  PSEC(KBDIM)

INTEGER :: KAERH(KBDIM,KLEV)

REAL(DP)::&
     &PCGAZ(KBDIM,KLEV)     &
  &,  PPIZAZ(KBDIM,KLEV)&
  &,  PRAY1(KBDIM,KLEV+1)  , PRAY2(KBDIM,KLEV+1)&
  &,  PREFZ(KBDIM,2,KLEV+1), PRJ(KBDIM,6,KLEV+1)&
  &,  PRK(KBDIM,6,KLEV+1)  , PRMU0(KBDIM,KLEV+1)&
  &,  PTAUAZ(KBDIM,KLEV)&
  &,  PTRA1(KBDIM,KLEV+1)  , PTRA2(KBDIM,KLEV+1)&
  &,  PTRCLR(KBDIM)        , PDIFFUSE(KBDIM)

!volcanic optical properties
REAL(DP) :: PSWEXT(KBDIM,KLEV,KNU), PSWSSA(KBDIM,KLEV,KNU) &
  &,        PSWASY(KBDIM,KLEV,KNU)

!new set of optical properties without the delta transformation
REAL(DP):: PTAUAZ_N(KBDIM,KLEV), PPIZAZ_N(KBDIM,KLEV), PCGAZ_N(KBDIM,KLEV)

!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL(DP):: ZC0I(KBDIM,KLEV+1)&
  &,  ZC0I2(KBDIM), ZCLEAR2(KBDIM) &
  &,  ZCLE0(KBDIM,KLEV), ZCLEAR(KBDIM) &
  &,  ZR21(KBDIM)&
  &,  ZR23(KBDIM) , ZSS0(KBDIM) , ZSCAT(KBDIM)&
  &,  ZTR(KBDIM,2,KLEV+1)

REAL(DP):: ZGG_1(KBDIM), ZTO1_1(KBDIM), ZW_1(KBDIM), ZREF(KBDIM)&
  &,       PR1(KBDIM), PT1(KBDIM), PR2(KBDIM), PT2(KBDIM)       &
  &,       ZRMUZ1(KBDIM), ZRMUZ(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IKL, JA, JAE, JAJ, JK, JKL, JKLP1, JKM1, JL, ICAE, IH

!     LOCAL REAL SCALARS
REAL(DP):: ZBMU0, ZBMU1, ZCORAE, ZDEN, ZDEN1, ZFACOA,&
          &ZFF, ZGAP, ZGAR, ZMU1, ZMUE, ZRATIO, ZRE11, &
          &ZTO, ZTRAY, ZWW, PRMU1, ZMUETO1
!     ------------------------------------------------------------------

!*         1.    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH
!                --------------------------------------------


DO JK = 1 , KLEV+1
  DO JA = 1 , 6
    DO JL = KIDIA,KFDIA
      PRJ(JL,JA,JK) = 0._DP
      PRK(JL,JA,JK) = 0._DP
    ENDDO
  ENDDO
ENDDO

! ------   NB: 'PAER' AEROSOLS ARE ENTERED FROM TOP TO BOTTOM

DO JK = 1 , KLEV
  IKL=KLEV+1-JK
  DO JL = KIDIA,KFDIA
    PCGAZ(JL,JK) = 0._DP
    PPIZAZ(JL,JK) =  0._DP
    PTAUAZ(JL,JK) = 0._DP

    PCGAZ_N(JL,JK) = 0._DP
    PPIZAZ_N(JL,JK) =  0._DP
    PTAUAZ_N(JL,JK) = 0._DP

  ENDDO
  DO JAE=1,KAER
    DO JL = KIDIA,KFDIA
      PTAUAZ(JL,JK)=PTAUAZ(JL,JK)+PAER(JL, IKL ,JAE)*RTAUA(KNU,JAE)
      PPIZAZ(JL,JK)=PPIZAZ(JL,JK)+PAER(JL, IKL ,JAE)&
       &* RTAUA(KNU,JAE)*RPIZA(KNU,JAE)
      PCGAZ(JL,JK) =  PCGAZ(JL,JK) +PAER(JL, IKL ,JAE)&
       &* RTAUA(KNU,JAE)*RPIZA(KNU,JAE)*RCGA(KNU,JAE)
    ENDDO
  ENDDO
  DO JAE=KAER+1,KAER+KNEWAER
    ICAE=NDFAER(JAE-KAER)
    DO JL = KIDIA,KFDIA
      IH=KAERH(JL,IKL)
      PTAUAZ(JL,JK)=PTAUAZ(JL,JK)+PAER(JL, IKL ,JAE)*RTAUAN(IH,KNU,ICAE)
      PPIZAZ(JL,JK)=PPIZAZ(JL,JK)+PAER(JL, IKL ,JAE)&
       &* RTAUAN(IH,KNU,ICAE)*RPIZAN(IH,KNU,ICAE)
      PCGAZ(JL,JK) =  PCGAZ(JL,JK) +PAER(JL, IKL ,JAE)&
       &* RTAUAN(IH,KNU,ICAE)*RPIZAN(IH,KNU,ICAE)*RCGAN(IH,KNU,ICAE)
    ENDDO
  ENDDO

!--MANU-- volcanic aerosols
 IF (L_VOLC) THEN
  DO JL = KIDIA, KFDIA
   PTAUAZ(JL,JK) = PTAUAZ(JL,JK) + PSWEXT(JL,IKL,KNU)
   PPIZAZ(JL,JK) = PPIZAZ(JL,JK) + PSWEXT(JL,IKL,KNU)*PSWSSA(JL,IKL,KNU)
   PCGAZ(JL,JK)  = PCGAZ(JL,JK)  + PSWEXT(JL,IKL,KNU)*PSWSSA(JL,IKL,KNU)  &
                                 * PSWASY(JL,IKL,KNU)
  ENDDO
 ENDIF
!---------------------------

  DO JL = KIDIA,KFDIA
    IF (KAER+KNEWAER /= 0 .AND. PPIZAZ(JL,JK)>EPSILON(1._DP)) THEN
      PCGAZ(JL,JK)=PCGAZ(JL,JK)/PPIZAZ(JL,JK)
      PPIZAZ(JL,JK)=PPIZAZ(JL,JK)/PTAUAZ(JL,JK)
      ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
      ZRATIO = ZTRAY / (ZTRAY + PTAUAZ(JL,JK))

!new set of optical properties without delta transformation

      PTAUAZ_N(JL,JK)= ZTRAY + PTAUAZ(JL,JK)
      PCGAZ_N(JL,JK) = (PTAUAZ(JL,JK)*PPIZAZ(JL,JK)*PCGAZ(JL,JK))&
       &/(ZTRAY+PTAUAZ(JL,JK)*PPIZAZ(JL,JK))
      PPIZAZ_N(JL,JK)= (ZTRAY+PTAUAZ(JL,JK)*PPIZAZ(JL,JK))&
       &/PTAUAZ_N(JL,JK)

!changes made to "pcgaz" in the set of transformed equations
!      ZGAR = PCGAZ(JL,JK)
!      ZFF = ZGAR * ZGAR
!      PTAUAZ(JL,JK)=ZTRAY+PTAUAZ(JL,JK)*(1._DP-PPIZAZ(JL,JK)*ZFF)
!      PCGAZ(JL,JK) = ZGAR * (1._DP - ZRATIO) / (1._DP + ZGAR)
!      PPIZAZ(JL,JK) =ZRATIO+(1._DP-ZRATIO)*PPIZAZ(JL,JK)*(1._DP-ZFF)&
!       &/ (1._DP - PPIZAZ(JL,JK) * ZFF)

!mixture of Rayleigh scatt and aerosols combined first by a weighted average
!and then a delta transformation done to the mixture
      ZFF=PCGAZ_N(JL,JK)*PCGAZ_N(JL,JK) 
      PTAUAZ(JL,JK)=PTAUAZ_N(JL,JK)*(1._DP-PPIZAZ_N(JL,JK)*ZFF)
      PCGAZ(JL,JK)=PCGAZ_N(JL,JK)/(1._DP+PCGAZ_N(JL,JK))
      PPIZAZ(JL,JK) =(PPIZAZ_N(JL,JK)*(1._DP-ZFF))&
       &/ (1._DP - PPIZAZ_N(JL,JK) * ZFF)

    ELSE
      ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
      PTAUAZ(JL,JK) = ZTRAY
      PCGAZ(JL,JK) = 0._DP
      PPIZAZ(JL,JK) = 1._DP-EPSILON(1._DP)
    ENDIF
  ENDDO
ENDDO

!*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
!                ----------------------------------------------


DO JL = KIDIA,KFDIA
  ZR23(JL) = 0._DP
  ZC0I(JL,KLEV+1) = 0._DP
  ZCLEAR(JL) = 1._DP
  ZCLEAR2(JL) = 1._DP
  ZSCAT(JL) = 0._DP
ENDDO

JK = 1
JKL = KLEV+1 - JK
JKLP1 = JKL + 1
DO JL = KIDIA,KFDIA

!repitition of applying the delta transformation to ptauaz

!  ZFACOA = 1._DP - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)

  ZFACOA = 1._DP 
  ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
  ZR21(JL) = EXP(-ZCORAE   )
  ZSS0(JL) = 1._DP-ZR21(JL)
  ZCLE0(JL,JKL) = ZSS0(JL)

  IF (NOVLP == 1) THEN
!* maximum-random      
    ZCLEAR(JL) = ZCLEAR(JL)&
     &*(1._DP-MAX(ZSS0(JL),ZSCAT(JL)))&
     &/(1._DP-MIN(ZSCAT(JL),1._DP-EPSILON(1._DP)))
    ZCLEAR2(JL) = ZCLEAR2(JL)*(1._DP-ZSS0(JL))
    ZC0I(JL,JKL) = 1._DP - ZCLEAR(JL)
    ZSCAT(JL) = ZSS0(JL)
  ELSEIF (NOVLP == 2) THEN
!* maximum
    ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
    ZCLEAR2(JL) = ZCLEAR2(JL)*(1._DP-ZSS0(JL))
    ZC0I(JL,JKL) = ZSCAT(JL)
  ELSEIF (NOVLP == 3) THEN
!* random
    ZCLEAR(JL)=ZCLEAR(JL)*(1._DP-ZSS0(JL))
    ZSCAT(JL) = 1._DP - ZCLEAR(JL)
    ZC0I(JL,JKL) = ZSCAT(JL)
  ENDIF
ENDDO

DO JK = 2 , KLEV
  JKL = KLEV+1 - JK
  JKLP1 = JKL + 1
  DO JL = KIDIA,KFDIA

!repitition of applying the delta transformation to ptauaz

!    ZFACOA = 1._DP - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
    
    ZFACOA = 1._DP 
    ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
    ZR21(JL) = EXP(-ZCORAE   )
    ZSS0(JL) = 1._DP-ZR21(JL)
    ZCLE0(JL,JKL) = ZSS0(JL)
    ZCLEAR2(JL) = ZCLEAR2(JL)*(1._DP-ZSS0(JL))

    IF (NOVLP == 1) THEN
!* maximum-random      
      ZCLEAR(JL) = ZCLEAR(JL)&
       &*(1._DP-MAX(ZSS0(JL),ZSCAT(JL)))&
       &/(1._DP-MIN(ZSCAT(JL),1._DP-EPSILON(1._DP)))
      ZC0I(JL,JKL) = 1._DP - ZCLEAR(JL)
      ZSCAT(JL) = ZSS0(JL)
      IF (JKL == 1) ZC0I2(JL) = 1._DP - ZCLEAR2(JL)
    ELSEIF (NOVLP == 2) THEN
!* maximum
      ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
      ZC0I(JL,JKL) = ZSCAT(JL)
      IF (JKL == 1) ZC0I2(JL) = 1._DP-ZCLEAR2(JL)
    ELSEIF (NOVLP == 3) THEN
!* random
      ZCLEAR(JL)=ZCLEAR(JL)*(1._DP-ZSS0(JL))
      ZSCAT(JL) = 1._DP - ZCLEAR(JL)
      ZC0I(JL,JKL) = ZSCAT(JL)
      IF (JKL == 1) ZC0I2(JL) = 1._DP - ZCLEAR(JL)
    ENDIF
  ENDDO
ENDDO
!     ------------------------------------------------------------------

!*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
!                -----------------------------------------------


DO JL = KIDIA,KFDIA
  PRAY1(JL,KLEV+1) = 0._DP
  PRAY2(JL,KLEV+1) = 0._DP
  PREFZ(JL,2,1) = PALBP(JL,KNU)
  PREFZ(JL,1,1) = PALBP(JL,KNU)
  PTRA1(JL,KLEV+1) = 1._DP
  PTRA2(JL,KLEV+1) = 1._DP
ENDDO

DO JK = 2 , KLEV+1
  JKM1 = JK-1
  DO JL = KIDIA,KFDIA
   PR1(JL) = 0._DP
   PT1(JL) = 0._DP
   PR2(JL) = 0._DP
   PT2(JL) = 0._DP

!     ------------------------------------------------------------------

!*         3.1  EQUIVALENT ZENITH ANGLE
!               -----------------------


    ZMUE = (1._DP-ZC0I(JL,JK)) * PSEC(JL)+ ZC0I(JL,JK) * 1.66_DP
    PRMU0(JL,JK) = 1._DP/ZMUE
    PRMU1 = 0.5_DP

!     ------------------------------------------------------------------

!*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
!               ----------------------------------------------------

!original
!
!    ZGAP = PCGAZ(JL,JKM1)
!    ZBMU0 = 0.5_DP - 0.75_DP * ZGAP / ZMUE
!    ZWW = PPIZAZ(JL,JKM1)
!    ZTO = PTAUAZ(JL,JKM1)
!    ZDEN = 1._DP + (1._DP - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE &
!     &+ (1-ZWW) * (1._DP - ZWW +2._DP*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE
!    PRAY1(JL,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE / ZDEN
!    PTRA1(JL,JKM1) = 1._DP / ZDEN
!
!
!    ZMU1 = 0.5_DP
!    ZBMU1 = 0.5_DP - 0.75_DP * ZGAP * ZMU1
!    ZDEN1= 1._DP + (1._DP - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1 &
!     &+ (1-ZWW) * (1._DP - ZWW +2._DP*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1
!    PRAY2(JL,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1
!    PTRA2(JL,JKM1) = 1._DP / ZDEN1

!MANU--modified

    ZREF(JL) = PREFZ(JL,1,JKM1)
    ZRMUZ(JL) = PRMU0(JL,JK)
    ZRMUZ1(JL) = PRMU1
    ZGG_1(JL)=PCGAZ_N(JL,JKM1)
    ZTO1_1(JL)=PTAUAZ_N(JL,JKM1)
    ZW_1(JL) = PPIZAZ_N(JL,JKM1)

 ENDDO

!pr1 and pt1 are reflectivity and transmissivity without reflection using the 
!equivalent zenith angle, prmu0

  CALL SWDE_WR ( KIDIA, KFDIA , KBDIM &
   &, ZGG_1  , ZREF  , ZRMUZ , ZTO1_1 , ZW_1 &
   &, PR1 , PT1 )

!pr2 and pt2 are reflectivity and transmissivity without reflection using the 
!fixed zenith angle, prmu1

 CALL SWDE_WR ( KIDIA, KFDIA , KBDIM &
   &, ZGG_1  , ZREF  , ZRMUZ1 , ZTO1_1 , ZW_1 &
   &, PR2 , PT2 )

 DO JL = KIDIA,KFDIA
   
   PRAY1(JL,JKM1) = PR1(JL)
   PTRA1(JL,JKM1) = PT1(JL)

   PRAY2(JL,JKM1) = PR2(JL)
   PTRA2(JL,JKM1) = PT2(JL)

    PREFZ(JL,1,JK) = (PRAY1(JL,JKM1)&
     &+ PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)&
     &* PTRA2(JL,JKM1)&
     &/ (1._DP-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))

    ZTR(JL,1,JKM1) = (PTRA1(JL,JKM1)&
     &/ (1._DP-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))

    PREFZ(JL,2,JK) = (PRAY1(JL,JKM1)&
     &+ PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)&
     &* PTRA2(JL,JKM1) )

    ZTR(JL,2,JKM1) = PTRA1(JL,JKM1)

  ENDDO
ENDDO

DO JL = KIDIA,KFDIA
  ZMUE = (1._DP-ZC0I(JL,1))*PSEC(JL)+ZC0I(JL,1)*1.66_DP
  PRMU0(JL,1)=1._DP/ZMUE
  PTRCLR(JL)=1._DP-ZC0I(JL,1)
  ZMUETO1 = 1._DP / (1._DP - ZC0I2(JL) * 0.5_DP)
  PDIFFUSE(JL) = ZC0I2(JL) * 0.5_DP * ZMUETO1
ENDDO


!     ------------------------------------------------------------------

!*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
!                 -------------------------------------------------


!!IF (KNU == 1) THEN
IF (KNU <= 3) THEN    !! new
  JAJ = 2
  DO JL = KIDIA,KFDIA
    PRJ(JL,JAJ,KLEV+1) = 1._DP
    PRK(JL,JAJ,KLEV+1) = PREFZ(JL, 1,KLEV+1)
  ENDDO

  DO JK = 1 , KLEV
    JKL = KLEV+1 - JK
    JKLP1 = JKL + 1
    DO JL = KIDIA,KFDIA
      ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,  1,JKL)
      PRJ(JL,JAJ,JKL) = ZRE11
      PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,  1,JKL)
    ENDDO
  ENDDO

ELSE

  DO JAJ = 1 , 2
    DO JL = KIDIA,KFDIA
      PRJ(JL,JAJ,KLEV+1) = 1._DP
      PRK(JL,JAJ,KLEV+1) = PREFZ(JL,JAJ,KLEV+1)
    ENDDO

    DO JK = 1 , KLEV
      JKL = KLEV+1 - JK
      JKLP1 = JKL + 1
      DO JL = KIDIA,KFDIA
        ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,JAJ,JKL)
        PRJ(JL,JAJ,JKL) = ZRE11
        PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,JAJ,JKL)
      ENDDO
    ENDDO
  ENDDO

ENDIF
!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SWCLR
