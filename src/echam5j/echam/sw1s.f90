SUBROUTINE SW1S &
 &( KIDIA , KFDIA , KBDIM , KLEV , KAER , KNEWAER, KAERH, KNU &
 &, PAER  , PSWEXT, PSWSSA, PSWASY      & !sw PADS optical prop. -MANU
 &, PALBD , PALBP, PCG  , PCLD , PCLEAR &
 &, PDSIG , POMEGA, POZ  , PRMU , PSEC , PTAU  , PUD  &
 &, PFD   , PFU   , PCD  , PCU  , PSUDU1 &
 &, P_VIS , P_VIS_FRACT_DIFFUSE &
 &)

!**** *SW1S* - SHORTWAVE RADIATION, FIRST SPECTRAL INTERVAL

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
!     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

!**   INTERFACE.
!     ----------

!          *SW1S* IS CALLED FROM *SW*.


!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES QUANTITIES FOR THE CLEAR-SKY FRACTION OF THE
!     COLUMN
!          2. COMPUTES UPWARD AND DOWNWARD FLUXES CORRESPONDING TO
!     CONTINUUM SCATTERING
!          3. MULTIPLY BY OZONE TRANSMISSION FUNCTION

!     EXTERNALS.
!     ----------

!          *SWCLR*, *SWR*, *SWTT*, *SWUVO3*

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
!        96-01-15   J.-J. MORCRETTE    SW in nsw SPECTRAL INTERVALS 
!        990128     JJMorcrette        sunshine duration
!        99-05-25   JJMorcrette        Revised aerosols
!        June 2005:  CCagnazzo 6 spectral bands from Morcrette code CY26R3

!     ------------------------------------------------------------------


USE MO_KIND  , ONLY : DP

USE MO_SW    , ONLY : NSW,    RRAY     ,RSUN


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

REAL(DP):: PAER(KBDIM,KLEV,KAER+KNEWAER)&
  &,  PALBD(KBDIM,NSW)      , PALBP(KBDIM,NSW)&
  &,  PCG(KBDIM,NSW,KLEV)   , PCLD(KBDIM,KLEV) &
  &,  PCLEAR(KBDIM)&
  &,  PDSIG(KBDIM,KLEV)&
  &,  POMEGA(KBDIM,NSW,KLEV), POZ(KBDIM,KLEV)&
  &,  PRMU(KBDIM)           , PSEC(KBDIM)&
  &,  PTAU(KBDIM,NSW,KLEV)  , PUD(KBDIM,5,KLEV+1)

INTEGER :: KAERH(KBDIM,KLEV)

REAL(DP):: PFD(KBDIM,KLEV+1)     , PFU(KBDIM,KLEV+1)&
  &,  PCD(KBDIM,KLEV+1)     , PCU(KBDIM,KLEV+1)&
  &,  PSUDU1(KBDIM)

!---MANU-- DECLARE
REAL(DP) :: PSWEXT(KBDIM,KLEV,NSW), PSWSSA(KBDIM,KLEV,NSW) &
  &, PSWASY(KBDIM,KLEV,NSW)

REAL(DP)::  P_VIS(KBDIM), P_VIS_FRACT_DIFFUSE(KBDIM)
!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

!!INTEGER :: IIND(6)
INTEGER :: IIND(4)

REAL(DP):: ZCGAZ(KBDIM,KLEV)&
  &,  ZDIFF(KBDIM)        , ZDIRF(KBDIM)        &
  &,  ZDIFT(KBDIM)        , ZDIRT(KBDIM)        &
  &,  ZPIZAZ(KBDIM,KLEV)&
  &,  ZRAYL(KBDIM), ZRAY1(KBDIM,KLEV+1), ZRAY2(KBDIM,KLEV+1)&
  &,  ZREFZ(KBDIM,2,KLEV+1)&
  &,  ZRJ(KBDIM,6,KLEV+1), ZRJ0(KBDIM,6,KLEV+1)&
  &,  ZRK(KBDIM,6,KLEV+1), ZRK0(KBDIM,6,KLEV+1)&
  &,  ZRMUE(KBDIM,KLEV+1), ZRMU0(KBDIM,KLEV+1)&
  &,  ZR(KBDIM,4)&			    ! new
!  &,  ZR(KBDIM,6)&			    ! old
  &,  ZTAUAZ(KBDIM,KLEV)&
  &,  ZTRA1(KBDIM,KLEV+1), ZTRA2(KBDIM,KLEV+1)&
  &,  ZTRCLD(KBDIM)      , ZTRCLR(KBDIM)&
!  &,  ZW(KBDIM,6)                          ! old
  &,  ZW(KBDIM,4),  ZO(KBDIM,2) ,ZT(KBDIM,2) & ! new
  &,  Z_FRACT_DIFFUSE_CLEAR(KBDIM), Z_FRACT_DIFFUSE_CLOUDY(KBDIM)
   
!-MANU--declaring the variables
REAL(DP)::ZTAUAZ_N(KBDIM,KLEV),ZPIZAZ_N(KBDIM,KLEV),ZCGAZ_N(KBDIM,KLEV)

!     LOCAL INTEGER SCALARS
INTEGER :: IKL, IKM1, JAJ, JK, JL


!     ------------------------------------------------------------------

!*         1.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)
!                 ----------------------- ------------------



!*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
!                 -----------------------------------------


DO JL = KIDIA,KFDIA
  ZRAYL(JL) =  RRAY(KNU,1) + PRMU(JL) * (RRAY(KNU,2) + PRMU(JL)&
   &* (RRAY(KNU,3) + PRMU(JL) * (RRAY(KNU,4) + PRMU(JL)&
   &* (RRAY(KNU,5) + PRMU(JL) *  RRAY(KNU,6)       ))))
ENDDO


!     ------------------------------------------------------------------

!*         2.    CONTINUUM SCATTERING CALCULATIONS
!                ---------------------------------


!*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
!                --------------------------------


CALL SWCLR &
  &( KIDIA  , KFDIA , KBDIM  , KLEV , KAER ,KNEWAER, KAERH, KNU &
  &, PAER   , PSWEXT, PSWSSA, PSWASY &
  &, PALBP , PDSIG , ZRAYL, PSEC &
  &, ZCGAZ  , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0 &
  &, ZRK0   , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2, ZTRCLR &
  &, Z_FRACT_DIFFUSE_CLEAR &
  &, ZTAUAZ_N, ZPIZAZ_N, ZCGAZ_N)


!*         2.2   CLOUDY FRACTION OF THE COLUMN
!                -----------------------------


CALL SWR &
  &( KIDIA ,KFDIA ,KBDIM  ,KLEV  , KNU &
  &, PALBD ,PCG   ,PCLD  ,POMEGA, PSEC , PTAU &
  &, ZCGAZ ,ZPIZAZ,ZRAY1 ,ZRAY2 , ZREFZ, ZRJ  ,ZRK , ZRMUE &
  &, ZTAUAZ,ZTRA1 ,ZTRA2 ,ZTRCLD &
  &, Z_FRACT_DIFFUSE_CLOUDY &
  &, ZTAUAZ_N, ZPIZAZ_N, ZCGAZ_N)


!     ------------------------------------------------------------------
!*         3.2   SIX SPECTRAL INTERVALS
!                ----------------------

!!!*         3.    OZONE ABSORPTION
!                ----------------


!!IIND(1)=1
!!IIND(2)=2
!!IIND(3)=3
!!IIND(4)=1
!!IIND(5)=2
!!IIND(6)=3

IIND(1)=1
IIND(2)=2
IIND(3)=1
IIND(4)=2



!*         3.1   DOWNWARD FLUXES
!                ---------------


JAJ = 2

DO JL = KIDIA,KFDIA
  ZW(JL,1)=0._DP
  ZW(JL,2)=0._DP
  ZW(JL,3)=0._DP
  ZW(JL,4)=0._DP

  ZO(JL,1)=0._DP  !! new
  ZO(JL,2)=0._DP  !! new
!!  ZW(JL,5)=0._DP
!!  ZW(JL,6)=0._DP
  PFD(JL,KLEV+1)=((1._DP-PCLEAR(JL))*ZRJ(JL,JAJ,KLEV+1)&
   &+ PCLEAR(JL) *ZRJ0(JL,JAJ,KLEV+1)) * RSUN(KNU)
  PCD(JL,KLEV+1)= ZRJ0(JL,JAJ,KLEV+1) * RSUN(KNU)
ENDDO
DO JK = 1 , KLEV
  IKL = KLEV+1-JK
  DO JL = KIDIA,KFDIA
    ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKL)/ZRMUE(JL,IKL)
    ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKL)/ZRMUE(JL,IKL)
    ZW(JL,3)=ZW(JL,3)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)   !! new
    ZW(JL,4)=ZW(JL,4)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)   !! new

!!    ZW(JL,3)=ZW(JL,3)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
!!    ZW(JL,4)=ZW(JL,4)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
!!    ZW(JL,5)=ZW(JL,5)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
!!    ZW(JL,6)=ZW(JL,6)+POZ(JL,  IKL)/ZRMU0(JL,IKL)

      ZO(JL,1)=ZO(JL,1)+POZ(JL,  IKL)/ZRMUE(JL,IKL)  !! new
      ZO(JL,2)=ZO(JL,2)+POZ(JL,  IKL)/ZRMU0(JL,IKL)  !! new
  ENDDO

!!  CALL SWTT1 ( KIDIA, KFDIA, KBDIM, KNU, 6 &
  CALL SWTT1 ( KIDIA, KFDIA, KBDIM, KNU, 4 &     !! new
   &, IIND &
   &, ZW  &
   &, ZR                          )

!! This is new:
   CALL SWUVO3 ( KIDIA, KFDIA, KBDIM, KNU, 2 &
     &, ZO  &
     &, ZT  &
     & )


  DO JL = KIDIA,KFDIA
!!    ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZR(JL,3)*ZRJ(JL,JAJ,IKL)
!!    ZDIRF(JL) = ZR(JL,4)*ZR(JL,5)*ZR(JL,6)*ZRJ0(JL,JAJ,IKL)
    ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZT(JL,1)*ZRJ(JL,JAJ,IKL)  !! new
    ZDIRF(JL) = ZR(JL,3)*ZR(JL,4)*ZT(JL,2)*ZRJ0(JL,JAJ,IKL) !! new
    PFD(JL,IKL) = ((1._DP-PCLEAR(JL)) * ZDIFF(JL)&
     &+PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
    PCD(JL,IKL) = ZDIRF(JL) * RSUN(KNU)
  ENDDO
    IF (IKL == 1) THEN
  DO JL = KIDIA,KFDIA
      P_VIS_FRACT_DIFFUSE(JL) = ((1._DP-PCLEAR(JL)) * ZDIFF(JL)        &
          * Z_FRACT_DIFFUSE_CLOUDY(JL) + PCLEAR(JL) * ZDIRF(JL)        &
          * Z_FRACT_DIFFUSE_CLEAR(JL)) * RSUN(KNU)                     &
          / (PFD(JL,IKL) + EPSILON(1._DP))
  ENDDO
    END IF
ENDDO

DO JL=KIDIA,KFDIA
!!  ZDIFT(JL) = ZR(JL,1)*ZR(JL,2)*ZR(JL,3)*ZTRCLD(JL)
!!  ZDIRT(JL) = ZR(JL,4)*ZR(JL,5)*ZR(JL,6)*ZTRCLR(JL)
  ZDIFT(JL) = ZR(JL,1)*ZR(JL,2)*ZT(JL,1)*ZTRCLD(JL)  !! new
  ZDIRT(JL) = ZR(JL,3)*ZR(JL,4)*ZT(JL,2)*ZTRCLR(JL)  !! new
  PSUDU1(JL) = ((1._DP-PCLEAR(JL)) * ZDIFT(JL)&
   &+PCLEAR(JL) * ZDIRT(JL)) * RSUN(KNU)
ENDDO


!*         3.2   UPWARD FLUXES
!                -------------


DO JL = KIDIA,KFDIA
  PFU(JL,1) = ((1._DP-PCLEAR(JL))*ZDIFF(JL)*PALBD(JL,KNU)&
   &+ PCLEAR(JL) *ZDIRF(JL)*PALBP(JL,KNU))&
   &* RSUN(KNU)
  P_VIS(JL) = PFD(JL,1) - PFU(JL,1)
  PCU(JL,1) = ZDIRF(JL) * PALBP(JL,KNU) * RSUN(KNU)
ENDDO

DO JK = 2 , KLEV+1
  IKM1=JK-1
  DO JL = KIDIA,KFDIA
    ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKM1)*1.66_DP
    ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKM1)*1.66_DP
    ZW(JL,3)=ZW(JL,3)+PUD(JL,1,IKM1)*1.66_DP
    ZW(JL,4)=ZW(JL,4)+PUD(JL,2,IKM1)*1.66_DP
!!    ZW(JL,5)=ZW(JL,5)+PUD(JL,2,IKM1)*1.66_DP
!!    ZW(JL,6)=ZW(JL,6)+POZ(JL,  IKM1)*1.66_DP
    
    ZO(JL,1)=ZO(JL,1)+POZ(JL,  IKM1)*1.66_DP   !! new
    ZO(JL,2)=ZO(JL,2)+POZ(JL,  IKM1)*1.66_DP   !! new
  ENDDO

!!  CALL SWTT1 ( KIDIA, KFDIA, KBDIM, KNU, 6 &
  CALL SWTT1 ( KIDIA, KFDIA, KBDIM, KNU, 4 &     !! new
   &, IIND &
   &, ZW  &
   &, ZR                          )

!! This is new:
  CALL SWUVO3 ( KIDIA, KFDIA, KBDIM, KNU, 2 &
    &, ZO  &
    &, ZT  &
    & )

  DO JL = KIDIA,KFDIA
!!    ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZR(JL,3)*ZRK(JL,JAJ,JK)
!!    ZDIRF(JL) = ZR(JL,4)*ZR(JL,5)*ZR(JL,6)*ZRK0(JL,JAJ,JK)
    ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZT(JL,1)*ZRK(JL,JAJ,JK)
    ZDIRF(JL) = ZR(JL,3)*ZR(JL,4)*ZT(JL,2)*ZRK0(JL,JAJ,JK)
    PFU(JL,JK) = ((1._DP-PCLEAR(JL)) * ZDIFF(JL)&
     &+PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
    PCU(JL,JK) = ZDIRF(JL) * RSUN(KNU)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SW1S
