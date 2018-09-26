SUBROUTINE SWNI &
 &( KIDIA , KFDIA , KBDIM  , KLEV , KAER, KNEWAER, KAERH, KNU &
 &, PAER  , PSWEXT, PSWSSA, PSWASY &! sw PADS optical prop. volc
 &, PAKI  , PALBD , PALBP, PCG &
 &, PCLD  , PCLEAR &
 &, PDSIG , POMEGA, POZ   , PRMU , PSEC, PTAU &
 &, PUD   , PWV   , PQS &
 &, PFDOWN, PFUP  , PCDOWN, PCUP , PSUDU2 &
 &, P_NIR, P_NIR_DIFFUSE &
 &)

!**** *SWNI* - SHORTWAVE RADIATION, NEAR-INFRARED SPECTRAL INTERVALS

!     PURPOSE.
!     --------

!          COMPUTES THE SHORTWAVE RADIATION FLUXES IN THE NEAR-INFRARED 
!     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

!**   INTERFACE.
!     ----------

!          *SWNI* IS CALLED FROM *SW*.


!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING TO
!     CONTINUUM SCATTERING
!          2. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING FOR
!     A GREY MOLECULAR ABSORPTION
!          3. LAPLACE TRANSFORM ON THE PREVIOUS TO GET EFFECTIVE AMOUNTS
!     OF ABSORBERS
!          4. APPLY H2O AND U.M.G. TRANSMISSION FUNCTIONS
!          5. MULTIPLY BY OZONE TRANSMISSION FUNCTION

!     EXTERNALS.
!     ----------

!          *SWCLR*, *SWR*, *SWDE*, *SWTT*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!        94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
!        95-12-07   J.-J. MORCRETTE    NEAR-INFRARED SW
!        990128     JJMorcrette        Sunshine duration
!        99-05-25   JJMorcrette        Revised aerosols

!     ------------------------------------------------------------------


USE MO_KIND  , ONLY : DP

USE MO_SW    , ONLY : NSW      ,RRAY     ,RSUN     ,RSWCE    ,RSWCP

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER :: KAER
INTEGER :: KNEWAER
INTEGER :: KFDIA
INTEGER :: KIDIA
INTEGER :: KLEV
INTEGER :: KBDIM
INTEGER :: KNU


!#include "yoeaer.h"
!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL(DP):: PAER(KBDIM,KLEV,KAER+KNEWAER)    , PAKI(KBDIM,2,NSW)&
  &,  PALBD(KBDIM,NSW)      , PALBP(KBDIM,NSW)&
  &,  PCG(KBDIM,NSW,KLEV)   , PCLD(KBDIM,KLEV)&
  &,  PCLEAR(KBDIM)         , PDSIG(KBDIM,KLEV)&
  &,  POMEGA(KBDIM,NSW,KLEV), POZ(KBDIM,KLEV)&
  &,  PQS(KBDIM,KLEV)&
  &,  PRMU(KBDIM)           , PSEC(KBDIM)&
  &,  PTAU(KBDIM,NSW,KLEV)  , PUD(KBDIM,5,KLEV+1)&
  &,  PWV(KBDIM,KLEV)


INTEGER :: KAERH(KBDIM,KLEV)

REAL(DP):: PFDOWN(KBDIM,KLEV+1)  , PFUP(KBDIM,KLEV+1)&
  &,  PCDOWN(KBDIM,KLEV+1)  , PCUP(KBDIM,KLEV+1)&
  &,  PSUDU2(KBDIM)
! volc:
REAL(DP) :: PSWEXT(KBDIM,KLEV,NSW), PSWSSA(KBDIM,KLEV,NSW) &
  &,        PSWASY(KBDIM,KLEV,NSW)

REAL(DP) ::  P_NIR(KBDIM), P_NIR_DIFFUSE(KBDIM)

!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

INTEGER :: IIND2(2), IIND3(3)
REAL(DP):: ZCGAZ(KBDIM,KLEV)  , ZDIFF(KBDIM)         , ZDIRF(KBDIM)&
  &,  ZFD(KBDIM,KLEV+1)  &
  &,  ZG(KBDIM)          , ZGG(KBDIM)&
  &,  ZPIZAZ(KBDIM,KLEV)&
  &,  ZRAYL(KBDIM)       , ZRAY1(KBDIM,KLEV+1)  , ZRAY2(KBDIM,KLEV+1)&
  &,  ZREF(KBDIM)        , ZREFZ(KBDIM,2,KLEV+1)&
  &,  ZRE1(KBDIM)        , ZRE2(KBDIM)&
  &,  ZRJ(KBDIM,6,KLEV+1), ZRJ0(KBDIM,6,KLEV+1)&
  &,  ZRK(KBDIM,6,KLEV+1), ZRK0(KBDIM,6,KLEV+1)&
  &,  ZRL(KBDIM,8)&
  &,  ZRMUE(KBDIM,KLEV+1), ZRMU0(KBDIM,KLEV+1)  , ZRMUZ(KBDIM)&
  &,  ZRNEB(KBDIM)       , ZRUEF(KBDIM,8)       , ZR1(KBDIM) &
  &,  ZR2(KBDIM,2)       , ZR3(KBDIM,3)         , ZR4(KBDIM)&
  &,  ZR21(KBDIM)        , ZR22(KBDIM)&
  &,  ZS(KBDIM)&
  &,  ZTAUAZ(KBDIM,KLEV) , ZTO1(KBDIM)          , ZTR(KBDIM,2,KLEV+1)&
  &,  ZTRA1(KBDIM,KLEV+1), ZTRA2(KBDIM,KLEV+1)&
  &,  ZTRCLD(KBDIM)      , ZTRCLR(KBDIM)&
  &,  ZTR1(KBDIM)        , ZTR2(KBDIM)&
  &,  ZW(KBDIM)          , ZW1(KBDIM)           , ZW2(KBDIM,2)&
  &,  ZW3(KBDIM,3)       , ZW4(KBDIM)           , ZW5(KBDIM) &
  &,  Z_FRACT_DIFFUSE_CLEAR(KBDIM), Z_FRACT_DIFFUSE_CLOUDY(KBDIM)

REAL(DP) :: ZWC1(KBDIM), ZWC4(KBDIM), ZWC5(KBDIM),ZRC4(KBDIM), &
  &         ZRC1(KBDIM)
REAL(DP) :: ZTAUAZ_N(KBDIM,KLEV), ZPIZAZ_N(KBDIM,KLEV), ZCGAZ_N(KBDIM,KLEV)
REAL(DP) :: ZFDOWNC(KBDIM,KLEV+1)  , ZFUPC(KBDIM,KLEV+1)


!     LOCAL INTEGER SCALARS
INTEGER :: IABS, IKL, IKM1, JABS, JAJ, JAJP, JK, JKKI,&
             &JKKP4, JKL, JKLP1, JKM1, JL, JN, JN2J, JREF

!     LOCAL REAL SCALARS
REAL(DP):: ZAA, ZBB, ZCNEB, ZRE11, ZRKI, ZRMUM1, ZWH2O
REAL(DP):: ZCHKS, ZCHKG, ZARGJ, ZARGK, ZARGJ0, ZARGK0

REAL(DP) :: ZEPS
!*         0.        Initializations
!                       ---------------
ZEPS=EPSILON(1.0_DP)

!     ------------------------------------------------------------------

!*         1.     NEAR-INFRARED SPECTRAL INTERVAL (0.68-4.00 MICRON)
!                 --------------------------------------------------



!*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
!                 -----------------------------------------


DO JL = KIDIA,KFDIA
  ZRMUM1 = 1._DP - PRMU(JL)
  ZRAYL(JL) =  RRAY(KNU,1) + ZRMUM1   * (RRAY(KNU,2) + ZRMUM1 &
   &* (RRAY(KNU,3) + ZRMUM1   * (RRAY(KNU,4) + ZRMUM1 &
   &* (RRAY(KNU,5) + ZRMUM1   *  RRAY(KNU,6)     ))))
  ZRAYL(JL) =  MAX(ZRAYL(JL),0._DP)
ENDDO


!     ------------------------------------------------------------------

!*         2.    CONTINUUM SCATTERING CALCULATIONS
!                ---------------------------------


!*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
!                --------------------------------


CALL SWCLR &
  &( KIDIA , KFDIA , KBDIM ,  KLEV , KAER ,KNEWAER, KAERH, KNU &
  &, PAER  , PSWEXT, PSWSSA, PSWASY &
  &, PALBP , PDSIG , ZRAYL,  PSEC &
  &, ZCGAZ , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0 &
  &, ZRK0  , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2, ZTRCLR &
  &, Z_FRACT_DIFFUSE_CLEAR &
  &, ZTAUAZ_N, ZPIZAZ_N, ZCGAZ_N)
!-MANU--additional variables

!*         2.1.A   SCATTERING CALCULATIONS WITH GREY MOLECULAR ABSORPTION
!                ------------------------------------------------------


JN = 2

DO JABS=1,2


!*         A.1  SURFACE CONDITIONS
!               ------------------


  DO JL = KIDIA,KFDIA
    ZREFZ(JL,2,1) = PALBD(JL,KNU)
    ZREFZ(JL,1,1) = PALBD(JL,KNU)
  ENDDO


!*         A.2  INTRODUCING CLOUD EFFECTS
!               -------------------------


  DO JK = 2 , KLEV+1
    JKM1 = JK - 1
    IKL=KLEV+1-JKM1
    DO JL = KIDIA,KFDIA
        ZAA=PUD(JL,JABS,JKM1)
        ZBB=ZAA
      ZRKI = PAKI(JL,JABS,KNU)
      ZCHKS= MIN(200._DP,ZRKI * ZAA * 1.66_DP)
      ZS(JL) = EXP(-ZCHKS)
      ZCHKG= MIN(200._DP,ZRKI * ZAA / ZRMU0(JL,JK))
      ZG(JL) = EXP(-ZCHKG)
      ZTR1(JL) = 0._DP
      ZRE1(JL) = 0._DP
      ZTR2(JL) = 0._DP
      ZRE2(JL) = 0._DP

!combining aerosol+gaseous absorption

!original ECHAM5
!      ZW(JL)= 1._DP
!      ZTO1(JL) = ZTAUAZ(JL,JKM1) / ZPIZAZ(JL,JKM1)&
!      &+ ZBB * ZRKI
!      ZR21(JL) =  ZTAUAZ(JL,JKM1)
!      ZR22(JL) = 0._DP
!      ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)&
!       &+ (1._DP - ZR22(JL)) * ZCGAZ(JL,JKM1)
!      ZW(JL) = ZR21(JL) / ZTO1(JL)

!-MANU--modified
      ZTO1(JL) =  ZTAUAZ_N(JL,JKM1) &
          &      +ZBB*ZRKI
      ZW(JL)   =  ZTAUAZ_N(JL,JKM1)*ZPIZAZ_N(JL,JKM1)
      ZGG(JL)  =  ZTAUAZ_N(JL,JKM1)*ZPIZAZ_N(JL,JKM1)*ZCGAZ_N(JL,JKM1)

!-asymmetry parameter

       IF(ZW(JL)>ZEPS) THEN
       ZGG(JL)  = ZGG(JL)/ZW(JL)
       ELSE
       ZGG(JL)  = 0.0_DP
       ENDIF

!-single scattering albedo

       IF(ZTO1(JL)>ZEPS) THEN
       ZW(JL)   = ZW(JL)/ZTO1(JL)
       ELSE
       ZW(JL)   = 1.0_DP
       ENDIF

       ZW(JL)=MIN(ZW(JL),1.0_DP-ZEPS)
!---

      ZREF(JL) = ZREFZ(JL,1,JKM1)
      ZRMUZ(JL) = ZRMU0(JL,JK)

    ENDDO


    CALL SWDE ( KIDIA, KFDIA, KBDIM &
     &, ZGG  , ZREF , ZRMUZ, ZTO1, ZW &
     &, ZRE1 , ZRE2 , ZTR1 , ZTR2     )

    DO JL = KIDIA,KFDIA

      ZREFZ(JL,2,JK) = (ZRAY1(JL,JKM1)&
       &+ ZREFZ(JL,2,JKM1) * ZTRA1(JL,JKM1)&
       &* ZTRA2(JL,JKM1)) * ZG(JL) * ZS(JL)

      ZTR(JL,2,JKM1)=ZTRA1(JL,JKM1) * ZG(JL)

      ZREFZ(JL,1,JK)=(ZRAY1(JL,JKM1)&
       &+ZREFZ(JL,1,JKM1)*ZTRA1(JL,JKM1)*ZTRA2(JL,JKM1)&
       &/(1._DP-ZRAY2(JL,JKM1)*ZREFZ(JL,1,JKM1)))*ZG(JL)*ZS(JL)

      ZTR(JL,1,JKM1)=(ZTRA1(JL,JKM1)/(1._DP-ZRAY2(JL,JKM1)&
       &* ZREFZ(JL,1,JKM1)))&
       &* ZG(JL)
    ENDDO
  ENDDO


!*         3.3  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
!               -------------------------------------------------


  DO JREF=1,2

    JN = JN + 1
    

    DO JL = KIDIA,KFDIA
      ZRJ0(JL,JN,KLEV+1) = 1._DP
      ZRK0(JL,JN,KLEV+1) = ZREFZ(JL,JREF,KLEV+1)
    ENDDO

    DO JK = 1 , KLEV
      JKL = KLEV+1 - JK
      JKLP1 = JKL + 1
      DO JL = KIDIA,KFDIA
        ZRE11 = ZRJ0(JL,JN,JKLP1) * ZTR(JL,JREF,JKL)
        ZRJ0(JL,JN,JKL) = ZRE11
        ZRK0(JL,JN,JKL) = ZRE11 * ZREFZ(JL,JREF,JKL)
      ENDDO
    ENDDO
  ENDDO
ENDDO


!     ------------------------------------------------------------------

!*         4.    INVERT GREY AND CONTINUUM FLUXES
!                --------------------------------



!*         4.1   UPWARD (ZRK) AND DOWNWARD (ZRJ) PSEUDO-FLUXES
!                ---------------------------------------------



DO JK = 1 , KLEV+1
  DO JAJ = 1 , 5 , 2
    JAJP = JAJ + 1
    DO JL = KIDIA,KFDIA
      ZRJ0(JL,JAJ,JK)=        ZRJ0(JL,JAJ,JK) - ZRJ0(JL,JAJP,JK)
      ZRK0(JL,JAJ,JK)=        ZRK0(JL,JAJ,JK) - ZRK0(JL,JAJP,JK)
      ZRJ0(JL,JAJ,JK)= MAX( ZRJ0(JL,JAJ,JK) , EPSILON(1._DP) )
      ZRK0(JL,JAJ,JK)= MAX( ZRK0(JL,JAJ,JK) , EPSILON(1._DP) )
    ENDDO
  ENDDO
ENDDO

DO JK = 1 , KLEV+1
  DO JAJ = 2 , 6 , 2
    DO JL = KIDIA,KFDIA
      ZRJ0(JL,JAJ,JK)= MAX( ZRJ0(JL,JAJ,JK) , EPSILON(1._DP) )
      ZRK0(JL,JAJ,JK)= MAX( ZRK0(JL,JAJ,JK) , EPSILON(1._DP) )
    ENDDO
  ENDDO
ENDDO
!
!*         4.2    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE
!                 ---------------------------------------------


DO JK = 1 , KLEV+1
  JKKI = 1
  DO JAJ = 1 , 2
    IIND2(1)=JAJ
    IIND2(2)=JAJ
    DO JN = 1 , 2
    JN2J = JN + 2 * JAJ
      JKKP4 = JKKI + 4

!*         4.2.1  EFFECTIVE ABSORBER AMOUNTS
!                 --------------------------


      DO JL = KIDIA,KFDIA
        ZARGJ0     = MAX( ZRJ0(JL,JN,JK) / ZRJ0(JL,JN2J,JK),1._DP)
        ZARGK0     = MAX( ZRK0(JL,JN,JK) / ZRK0(JL,JN2J,JK),1._DP)
        ZW2(JL,1) = LOG( ZARGJ0 )/ PAKI(JL,JAJ,KNU)
        ZW2(JL,2) = LOG( ZARGK0 )/ PAKI(JL,JAJ,KNU)
      ENDDO
!
!*         4.2.2  TRANSMISSION FUNCTION
!                 ---------------------


      CALL SWTT1 ( KIDIA,KFDIA,KBDIM, KNU, 2, IIND2 &
       &, ZW2 &
       &, ZR2                              )

      DO JL = KIDIA,KFDIA
        ZRL(JL,JKKI) = ZR2(JL,1)
        ZRUEF(JL,JKKI) = ZW2(JL,1)
        ZRL(JL,JKKP4) = ZR2(JL,2)
        ZRUEF(JL,JKKP4) = ZW2(JL,2)
      ENDDO

      JKKI=JKKI+1
    ENDDO
  ENDDO

!*         4.3    UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION
!                 ------------------------------------------------------


  DO JL = KIDIA,KFDIA
    ZFDOWNC(JL,JK) = ZRJ0(JL,1,JK) * ZRL(JL,1) * ZRL(JL,3)&
     &+ ZRJ0(JL,2,JK) * ZRL(JL,2) * ZRL(JL,4)
    ZFUPC(JL,JK)   = ZRK0(JL,1,JK) * ZRL(JL,5) * ZRL(JL,7)&
     &+ ZRK0(JL,2,JK) * ZRL(JL,6) * ZRL(JL,8)
  ENDDO
ENDDO




!*         2.2   CLOUDY FRACTION OF THE COLUMN
!                -----------------------------


CALL SWR &
  &( KIDIA , KFDIA , KBDIM , KLEV  , KNU &
  &, PALBD , PCG   , PCLD , POMEGA, PSEC , PTAU &
  &, ZCGAZ , ZPIZAZ, ZRAY1, ZRAY2 , ZREFZ, ZRJ  , ZRK, ZRMUE &
  &, ZTAUAZ, ZTRA1 , ZTRA2, ZTRCLD &
  &, Z_FRACT_DIFFUSE_CLOUDY &
  &, ZTAUAZ_N, ZPIZAZ_N, ZCGAZ_N)
!-MANU--additional variables

!     ------------------------------------------------------------------

!*         3.    SCATTERING CALCULATIONS WITH GREY MOLECULAR ABSORPTION
!                ------------------------------------------------------


JN = 2

DO JABS=1,2


!*         3.1  SURFACE CONDITIONS
!               ------------------


  DO JL = KIDIA,KFDIA
    ZREFZ(JL,2,1) = PALBD(JL,KNU)
    ZREFZ(JL,1,1) = PALBD(JL,KNU)
  ENDDO


!*         3.2  INTRODUCING CLOUD EFFECTS
!               -------------------------

  DO JK = 2 , KLEV+1
    JKM1 = JK - 1
    IKL=KLEV+1-JKM1
    DO JL = KIDIA,KFDIA

      ZRNEB(JL) = PCLD(JL,JKM1)

      IF (JABS == 1.AND. ZRNEB(JL) > 2._DP*EPSILON(1._DP)) THEN
        ZWH2O=MAX(PWV(JL,IKL),EPSILON(1._DP))
        ZCNEB=MAX(EPSILON(1._DP),MIN(ZRNEB(JL),1._DP-EPSILON(1._DP)))
        ZBB=PUD(JL,JABS,JKM1)*PQS(JL,IKL)/ZWH2O
        ZAA=MAX((PUD(JL,JABS,JKM1)-ZCNEB*ZBB)/(1._DP-ZCNEB),EPSILON(1._DP))
      ELSE
        ZAA=PUD(JL,JABS,JKM1)
        ZBB=ZAA
      ENDIF
      ZRKI = PAKI(JL,JABS,KNU)
      ZCHKS= MIN(200._DP,ZRKI * ZAA * 1.66_DP)
      ZS(JL) = EXP(-ZCHKS)
      ZCHKG= MIN(200._DP,ZRKI * ZAA / ZRMUE(JL,JK))
      ZG(JL) = EXP(-ZCHKG)
      ZTR1(JL) = 0._DP
      ZRE1(JL) = 0._DP
      ZTR2(JL) = 0._DP
      ZRE2(JL) = 0._DP

!combining aerosols+clouds+gaseous absorption

!original ECHAM5
!      ZW(JL)= POMEGA(JL,KNU,JKM1)
!      ZTO1(JL) = PTAU(JL,KNU,JKM1) / ZW(JL)&
!      &+ ZTAUAZ(JL,JKM1) / ZPIZAZ(JL,JKM1) &
!      &+ ZBB * ZRKI
!      ZR21(JL) = PTAU(JL,KNU,JKM1) + ZTAUAZ(JL,JKM1)
!      ZR22(JL) = PTAU(JL,KNU,JKM1) / ZR21(JL)
!      ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)&
!       &+ (1._DP - ZR22(JL)) * ZCGAZ(JL,JKM1)
!      ZW(JL) = ZR21(JL) / ZTO1(JL)

!-MANU--modified
      ZTO1(JL) =  ZTAUAZ_N(JL,JKM1) &
           &      +PTAU(JL,KNU,JKM1) &
           &      +ZBB*ZRKI
       ZW(JL)   =  ZTAUAZ_N(JL,JKM1)*ZPIZAZ_N(JL,JKM1) &
           &      +PTAU(JL,KNU,JKM1)*POMEGA(JL,KNU,JKM1)
       ZGG(JL)  =  ZTAUAZ_N(JL,JKM1)*ZPIZAZ_N(JL,JKM1)*ZCGAZ_N(JL,JKM1)&
           &      +PTAU(JL,KNU,JKM1)*POMEGA(JL,KNU,JKM1)*PCG(JL,KNU,JKM1)

       IF(ZW(JL)>ZEPS) THEN
       ZGG(JL)  = ZGG(JL)/ZW(JL)
       ELSE
       ZGG(JL)  = 0.0_DP
       ENDIF

       IF(ZTO1(JL)>ZEPS) THEN
       ZW(JL)   = ZW(JL)/ZTO1(JL)
       ELSE
       ZW(JL)=1.0_DP
       ENDIF

       ZW(JL)=MIN(ZW(JL),1.0_DP-ZEPS)

!----
      ZREF(JL) = ZREFZ(JL,1,JKM1)
      ZRMUZ(JL) = ZRMUE(JL,JK)
    ENDDO

    CALL SWDE ( KIDIA, KFDIA, KBDIM &
     &, ZGG  , ZREF , ZRMUZ, ZTO1, ZW &
     &, ZRE1 , ZRE2 , ZTR1 , ZTR2     )

    DO JL = KIDIA,KFDIA


      ZREFZ(JL,2,JK) = (1._DP-ZRNEB(JL)) * (ZRAY1(JL,JKM1)&
       &+ ZREFZ(JL,2,JKM1) * ZTRA1(JL,JKM1)&
       &* ZTRA2(JL,JKM1) ) * ZG(JL) * ZS(JL)&
       &+ ZRNEB(JL) * ZRE1(JL)

      ZTR(JL,2,JKM1)=ZRNEB(JL)*ZTR1(JL)&
       &+ (ZTRA1(JL,JKM1)) * ZG(JL) * (1._DP-ZRNEB(JL))

 
      ZREFZ(JL,1,JK)=(1._DP-ZRNEB(JL))*(ZRAY1(JL,JKM1)&
       &+ZREFZ(JL,1,JKM1)*ZTRA1(JL,JKM1)*ZTRA2(JL,JKM1)&
       &/(1._DP-ZRAY2(JL,JKM1)*ZREFZ(JL,1,JKM1)))*ZG(JL)*ZS(JL)&
       &+ ZRNEB(JL) * ZRE2(JL)

      ZTR(JL,1,JKM1)= ZRNEB(JL)*ZTR2(JL)&
       &+ (ZTRA1(JL,JKM1)/(1._DP-ZRAY2(JL,JKM1)&
       &* ZREFZ(JL,1,JKM1)))&
       &* ZG(JL) * (1._DP -ZRNEB(JL))

    ENDDO
  ENDDO


!*         3.3  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
!               -------------------------------------------------


  DO JREF=1,2

    JN = JN + 1
    

    DO JL = KIDIA,KFDIA
      ZRJ(JL,JN,KLEV+1) = 1._DP
      ZRK(JL,JN,KLEV+1) = ZREFZ(JL,JREF,KLEV+1)
    ENDDO

    DO JK = 1 , KLEV
      JKL = KLEV+1 - JK
      JKLP1 = JKL + 1
      DO JL = KIDIA,KFDIA
        ZRE11 = ZRJ(JL,JN,JKLP1) * ZTR(JL,JREF,JKL)
        ZRJ(JL,JN,JKL) = ZRE11
        ZRK(JL,JN,JKL) = ZRE11 * ZREFZ(JL,JREF,JKL)
      ENDDO
    ENDDO
  ENDDO
ENDDO


!     ------------------------------------------------------------------

!*         4.    INVERT GREY AND CONTINUUM FLUXES
!                --------------------------------



!*         4.1   UPWARD (ZRK) AND DOWNWARD (ZRJ) PSEUDO-FLUXES
!                ---------------------------------------------



DO JK = 1 , KLEV+1
  DO JAJ = 1 , 5 , 2
    JAJP = JAJ + 1
    DO JL = KIDIA,KFDIA
      ZRJ(JL,JAJ,JK)=        ZRJ(JL,JAJ,JK) - ZRJ(JL,JAJP,JK)
      ZRK(JL,JAJ,JK)=        ZRK(JL,JAJ,JK) - ZRK(JL,JAJP,JK)
      ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , EPSILON(1._DP) )
      ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , EPSILON(1._DP) )
    ENDDO
  ENDDO
ENDDO

DO JK = 1 , KLEV+1
  DO JAJ = 2 , 6 , 2
    DO JL = KIDIA,KFDIA
      ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , EPSILON(1._DP) )
      ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , EPSILON(1._DP) )
    ENDDO
  ENDDO
ENDDO
!
!*         4.2    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE
!                 ---------------------------------------------


DO JK = 1 , KLEV+1
  JKKI = 1
  DO JAJ = 1 , 2
    IIND2(1)=JAJ
    IIND2(2)=JAJ
    DO JN = 1 , 2
    JN2J = JN + 2 * JAJ
      JKKP4 = JKKI + 4

!*         4.2.1  EFFECTIVE ABSORBER AMOUNTS
!                 --------------------------


      DO JL = KIDIA,KFDIA
        ZARGJ     = MAX( ZRJ(JL,JN,JK) / ZRJ(JL,JN2J,JK),1._DP)
        ZARGK     = MAX( ZRK(JL,JN,JK) / ZRK(JL,JN2J,JK),1._DP)
        ZW2(JL,1) = LOG( ZARGJ )/ PAKI(JL,JAJ,KNU)
        ZW2(JL,2) = LOG( ZARGK )/ PAKI(JL,JAJ,KNU)
      ENDDO
!
!*         4.2.2  TRANSMISSION FUNCTION
!                 ---------------------


      CALL SWTT1 ( KIDIA,KFDIA,KBDIM, KNU, 2, IIND2 &
       &, ZW2 &
       &, ZR2                              )

      DO JL = KIDIA,KFDIA
        ZRL(JL,JKKI) = ZR2(JL,1)
        ZRUEF(JL,JKKI) = ZW2(JL,1)
        ZRL(JL,JKKP4) = ZR2(JL,2)
        ZRUEF(JL,JKKP4) = ZW2(JL,2)
      ENDDO

      JKKI=JKKI+1
    ENDDO
  ENDDO

!*         4.3    UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION
!                 ------------------------------------------------------


  DO JL = KIDIA,KFDIA
    PFDOWN(JL,JK) = ZRJ(JL,1,JK) * ZRL(JL,1) * ZRL(JL,3)&
     &+ ZRJ(JL,2,JK) * ZRL(JL,2) * ZRL(JL,4)
    PFUP(JL,JK)   = ZRK(JL,1,JK) * ZRL(JL,5) * ZRL(JL,7)&
     &+ ZRK(JL,2,JK) * ZRL(JL,6) * ZRL(JL,8)
  ENDDO
ENDDO



!     ------------------------------------------------------------------

!*         5.    MOLECULAR ABSORPTION ON CLEAR-SKY FLUXES
!                ----------------------------------------



!*         5.1   DOWNWARD FLUXES
!                ---------------


!JAJ = 2
IIND3(1)=1
IIND3(2)=2
IIND3(3)=3

DO JL = KIDIA,KFDIA
  ZW3(JL,1)=0._DP
  ZW3(JL,2)=0._DP
  ZW3(JL,3)=0._DP
  ZW4(JL)  =0._DP
  ZW5(JL)  =0._DP
  ZR4(JL)  =1._DP
! ZFD(JL,KLEV+1)= PFDOWNC(JL,KLEV+1)
ENDDO
DO JK = 1 , KLEV
  IKL = KLEV+1-JK
  DO JL = KIDIA,KFDIA
    ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
    ZW3(JL,2)=ZW3(JL,2)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
    ZW3(JL,3)=ZW3(JL,3)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
    ZW4(JL)  =ZW4(JL)  +PUD(JL,4,IKL)/ZRMU0(JL,IKL)
    ZW5(JL)  =ZW5(JL)  +PUD(JL,5,IKL)/ZRMU0(JL,IKL)
  ENDDO

  CALL SWTT1 ( KIDIA,KFDIA,KBDIM, KNU, 3, IIND3 &
   &, ZW3 &
   &, ZR3                              )

  DO JL = KIDIA,KFDIA
    ZR4(JL) = EXP(-RSWCE(KNU)*ZW4(JL)-RSWCP(KNU)*ZW5(JL))
!    ZFD(JL,IKL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)* PFDOWNC(JL,IKL)
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  ZDIFF(JL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)*ZTRCLD(JL)
  ZDIRF(JL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)*ZTRCLR(JL)
  PSUDU2(JL) = ((1._DP-PCLEAR(JL)) * ZDIFF(JL)&
   &+PCLEAR(JL) * ZDIRF(JL)) * RSUN(KNU)
ENDDO



!*         6.     INTRODUCTION OF OZONE AND H2O CONTINUUM ABSORPTION
!                 --------------------------------------------------

IABS=3

!*         6.1    DOWNWARD FLUXES
!                 ---------------

DO JL = KIDIA,KFDIA
  ZW1(JL)=0._DP
  ZW4(JL)=0._DP
  ZW5(JL)=0._DP
  ZR1(JL)=0._DP
  ZWC1(JL)=0._DP
  ZWC4(JL)=0._DP
  ZWC5(JL)=0._DP
  ZRC1(JL)=0._DP
  PFDOWN(JL,KLEV+1) = ((1._DP-PCLEAR(JL))*PFDOWN(JL,KLEV+1)&
   &+ PCLEAR(JL) * ZFDOWNC(JL,KLEV+1)) * RSUN(KNU)
  PCDOWN(JL,KLEV+1) = ZFDOWNC(JL,KLEV+1) * RSUN(KNU)
ENDDO

DO JK = 1 , KLEV
  IKL=KLEV+1-JK
  DO JL = KIDIA,KFDIA
    ZW1(JL) = ZW1(JL)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
    ZW4(JL) = ZW4(JL)+PUD(JL,4,IKL)/ZRMUE(JL,IKL)
    ZW5(JL) = ZW5(JL)+PUD(JL,5,IKL)/ZRMUE(JL,IKL)
    ZR4(JL) = EXP(-RSWCE(KNU)*ZW4(JL)-RSWCP(KNU)*ZW5(JL))

!gases absorption for clear sky fluxes

    ZWC1(JL) = ZWC1(JL)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
    ZWC4(JL) = ZWC4(JL)+PUD(JL,4,IKL)/ZRMU0(JL,IKL)
    ZWC5(JL) = ZWC5(JL)+PUD(JL,5,IKL)/ZRMU0(JL,IKL)
    ZRC4(JL) = EXP(-RSWCE(KNU)*ZWC4(JL)-RSWCP(KNU)*ZWC5(JL))
  ENDDO

  CALL SWTT ( KIDIA,KFDIA,KBDIM, KNU, IABS, ZW1, ZR1 )

!gases absorption for clear sky fluxes

  CALL SWTT ( KIDIA,KFDIA,KBDIM, KNU, IABS, ZWC1, ZRC1 )

 IF (IKL == 1) THEN
  DO JL = KIDIA,KFDIA
      P_NIR_DIFFUSE(JL) = ((1._DP-PCLEAR(JL))*ZR1(JL)*ZR4(JL)                       &
       & * PFDOWN(JL,IKL)*Z_FRACT_DIFFUSE_CLOUDY(JL)                                &
       & + PCLEAR(JL)*ZRC1(JL)*ZRC4(JL)*ZFDOWNC(JL,IKL)*Z_FRACT_DIFFUSE_CLEAR(JL))  &
       & * RSUN(KNU)
  ENDDO
 END IF
  DO JL = KIDIA,KFDIA
    PFDOWN(JL,IKL) = ((1._DP-PCLEAR(JL))*ZR1(JL)*ZR4(JL)*PFDOWN(JL,IKL)&
     &+PCLEAR(JL)*ZRC1(JL)*ZRC4(JL)*ZFDOWNC(JL,IKL)) * RSUN(KNU)
    PCDOWN(JL,IKL) = ZRC1(JL)*ZRC4(JL)*ZFDOWNC(JL,IKL) * RSUN(KNU)
  ENDDO
 IF (IKL == 1) THEN
  DO JL = KIDIA,KFDIA
      P_NIR_DIFFUSE(JL) = P_NIR_DIFFUSE(JL)/(PFDOWN(JL,IKL)+EPSILON(1._DP))
  ENDDO
 END IF
ENDDO


!*         6.2    UPWARD FLUXES
!                 -------------

DO JL = KIDIA,KFDIA
  PFUP(JL,1) = ((1._DP-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,1)&
   &+PCLEAR(JL)*ZRC1(JL)*ZRC4(JL)*ZFDOWNC(JL,1)*PALBP(JL,KNU)) * RSUN(KNU)
  P_NIR(JL) = PFDOWN(JL,1) - PFUP(JL,1)
  PCUP(JL,1) = ZRC1(JL)*ZRC4(JL)*ZFDOWNC(JL,1)*PALBP(JL,KNU)* RSUN(KNU)
ENDDO

DO JK = 2 , KLEV+1
  IKM1=JK-1
  DO JL = KIDIA,KFDIA
    ZW1(JL) = ZW1(JL)+POZ(JL  ,IKM1)*1.66_DP
    ZW4(JL) = ZW4(JL)+PUD(JL,4,IKM1)*1.66_DP
    ZW5(JL) = ZW5(JL)+PUD(JL,5,IKM1)*1.66_DP
    ZR4(JL) = EXP(-RSWCE(KNU)*ZW4(JL)-RSWCP(KNU)*ZW5(JL))
  ENDDO

  CALL SWTT ( KIDIA,KFDIA,KBDIM, KNU, IABS, ZW1, ZR1 )

  DO JL = KIDIA,KFDIA
    PFUP(JL,JK) = ((1._DP-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,JK)&
     &+PCLEAR(JL)*ZR1(JL)*ZR4(JL) * ZFUPC(JL,JK)) * RSUN(KNU)
    PCUP(JL,JK) = ZR1(JL)*ZR4(JL) * ZFUPC(JL,JK) * RSUN(KNU)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SWNI
