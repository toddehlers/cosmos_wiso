!OPTIONS XOPT(HSFUN)
SUBROUTINE SWDE_WR &
  &( KIDIA, KFDIA, KBDIM &
  &, PGG  , PREF , PRMUZ, PTO1, PW &
  &, PRE1 , PTR1 &
  &)

!**** *SWDE_WR* - DELTA-EDDINGTON IN A CLOUDY LAYER

!     PURPOSE.
!     --------
!           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLOUDY
!     LAYER USING THE DELTA-EDDINGTON'S APPROXIMATION.

!**   INTERFACE.
!     ----------
!          *SWDE* IS CALLED BY *SWR*, *SWNI*


!        EXPLICIT ARGUMENTS :
!        --------------------
! PGG    : (KBDIM)             ; ASSYMETRY FACTOR
! PREF   : (KBDIM)             ; REFLECTIVITY OF THE UNDERLYING LAYER
! PRMUZ  : (KBDIM)             ; COSINE OF SOLAR ZENITH ANGLE
! PTO1   : (KBDIM)             ; OPTICAL THICKNESS
! PW     : (KBDIM)             ; SINGLE SCATTERING ALBEDO
!     ==== OUTPUTS ===
! PRE1   : (KBDIM)             ; LAYER REFLECTIVITY ASSUMING NO
!                             ; REFLECTION FROM UNDERLYING LAYER
! PTR1   : (KBDIM)             ; LAYER TRANSMISSIVITY ASSUMING NO
!                             ; REFLECTION FROM UNDERLYING LAYER

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          STANDARD DELTA-EDDINGTON LAYER CALCULATIONS.

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 88-12-15
!                   96-05-30 Michel Deque (security in EXP()) 
    
!     ------------------------------------------------------------------


!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------


USE MO_KIND  , ONLY : DP

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER :: KFDIA
INTEGER :: KIDIA
INTEGER :: KBDIM

REAL(DP):: PGG(KBDIM),PREF(KBDIM),PRMUZ(KBDIM),PTO1(KBDIM),PW(KBDIM)
REAL(DP):: PRE1(KBDIM),PTR1(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: JL

!     LOCAL REAL SCALARS
REAL(DP):: ZA11, ZA12, ZA13, ZA21, ZA22, ZA23, ZALPHA,&
          &ZAM2B, ZAP2B, ZARG, ZARG2, &
          &ZBETA, ZC1A, ZC2A, ZDENA,  &
          &ZDT, ZEXKM, ZEXKP, ZEXMU0, ZFF, ZGP, ZRI0A, &
          &ZRI0B, ZRI1A, ZRI1B, &
          &ZRK, ZRM2, ZRP, ZTOP, ZWCP, ZWM, ZX1, &
          &ZX2, ZXM2P, ZXP2P



!     ------------------------------------------------------------------

!*         1.      DELTA-EDDINGTON CALCULATIONS


DO JL   =   KIDIA,KFDIA

!*         1.1     SET UP THE DELTA-MODIFIED PARAMETERS


  ZFF = PGG(JL)*PGG(JL)
  ZGP = PGG(JL)/(1._DP+PGG(JL))
  ZTOP = (1._DP- PW(JL) * ZFF) * PTO1(JL)
  ZWCP = (1-ZFF)* PW(JL) /(1._DP- PW(JL) * ZFF)
  ZDT = 2._DP/3._DP
  ZX1 = 1._DP-ZWCP*ZGP
  ZWM = 1._DP-ZWCP
  ZRM2 =  PRMUZ(JL) * PRMUZ(JL)
  ZRK = SQRT(3._DP*ZWM*ZX1)
  ZX2 = 4._DP*(1._DP-ZRK*ZRK*ZRM2)
  ZRP=ZRK/ZX1
  ZALPHA = 3._DP*ZWCP*ZRM2*(1._DP+ZGP*ZWM)/ZX2
  ZBETA = 3._DP*ZWCP* PRMUZ(JL) *(1._DP+3._DP*ZGP*ZRM2*ZWM)/ZX2
!     ZARG=MIN(ZTOP/PRMUZ(JL),200.)
  ZARG=MAX(-200._DP,MIN(ZTOP/PRMUZ(JL),200._DP))
  ZEXMU0=EXP(-ZARG)
  ZARG2=MIN(ZRK*ZTOP,200._DP)
  ZEXKP=EXP(ZARG2)
  ZEXKM = 1._DP/ZEXKP
  ZXP2P = 1._DP+ZDT*ZRP
  ZXM2P = 1._DP-ZDT*ZRP
  ZAP2B = ZALPHA+ZDT*ZBETA
  ZAM2B = ZALPHA-ZDT*ZBETA

!*         1.2     WITHOUT REFLECTION FROM THE UNDERLYING LAYER


  ZA11 = ZXP2P
  ZA12 = ZXM2P
  ZA13 = ZAP2B
  ZA22 = ZXP2P*ZEXKP
  ZA21 = ZXM2P*ZEXKM
  ZA23 = ZAM2B*ZEXMU0
  ZDENA = ZA11 * ZA22 - ZA21 * ZA12
  ZC1A = (ZA22*ZA13-ZA12*ZA23)/ZDENA
  ZC2A = (ZA11*ZA23-ZA21*ZA13)/ZDENA
  ZRI0A = ZC1A+ZC2A-ZALPHA
  ZRI1A = ZRP*(ZC1A-ZC2A)-ZBETA
  PRE1(JL) = (ZRI0A-ZDT*ZRI1A)/ PRMUZ(JL)
  ZRI0B = ZC1A*ZEXKM+ZC2A*ZEXKP-ZALPHA*ZEXMU0
  ZRI1B = ZRP*(ZC1A*ZEXKM-ZC2A*ZEXKP)-ZBETA*ZEXMU0
  PTR1(JL) = ZEXMU0+(ZRI0B+ZDT*ZRI1B)/ PRMUZ(JL)

ENDDO
RETURN
END SUBROUTINE SWDE_WR
