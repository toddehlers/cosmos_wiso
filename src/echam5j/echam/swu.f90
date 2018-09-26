!OPTIONS XOPT(HSFUN)
SUBROUTINE SWU &
  &( KIDIA, KFDIA , KBDIM  , KLEV &
  &, PSCT , PCARDI, PPMB , PPSOL, PRMU0, PTAVE, PWV &
  &, PAKI , PDSIG , PFACT, PRMU , PSEC , PUD &
  &)

!**** *SWU* - SHORTWAVE RADIATION, ABSORBER AMOUNTS

!     PURPOSE.
!     --------
!           COMPUTES THE ABSORBER AMOUNTS USED IN SHORTWAVE RADIATION
!     CALCULATIONS

!**   INTERFACE.
!     ----------
!          *SWU* IS CALLED BY *SW*


!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES ABSORBER AMOUNTS WITH TEMPERATURE AND PRESSURE
!     SCALING.

!     EXTERNALS.
!     ----------

!          *SWTT*

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
!        June 2005: CCagnazzo from 4 to 6 bands, following 26Rcycle Morcrette Code

!     ------------------------------------------------------------------


USE MO_KIND  , ONLY : DP

USE MO_SW    , ONLY : NSW      ,&
            &         RPDH1    ,RPDU1    ,RPNH     ,RPNU     ,&
            &         RTDH2O   ,RTDUMG   ,RTH2O    ,RTUMG


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER :: KFDIA
INTEGER :: KIDIA
INTEGER :: KLEV
INTEGER :: KBDIM

!     DUMMY REAL SCALARS
REAL(DP):: PSCT



!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL(DP):: PPMB(KBDIM,KLEV+1), PPSOL(KBDIM)&
  &,  PRMU0(KBDIM)      , PTAVE(KBDIM,KLEV) , PWV(KBDIM,KLEV) &
  &,  PCARDI(KBDIM,KLEV)

REAL(DP):: PAKI(KBDIM,2,NSW)&
  &,  PDSIG(KBDIM,KLEV) , PFACT(KBDIM)      , PRMU(KBDIM)&
  &,  PSEC(KBDIM)       , PUD(KBDIM,5,KLEV+1)

!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

INTEGER :: IIND(2)
REAL(DP):: ZN175(KBDIM), ZN190(KBDIM), ZO175(KBDIM)&
  &,  ZO190(KBDIM), ZSIGN(KBDIM)&
  &,  ZR(KBDIM,2) , ZSIGO(KBDIM), ZUD(KBDIM,2)

!     LOCAL INTEGER SCALARS
INTEGER :: JA, JK, JKL, JKLP1, JKP1, JL, JNU

!     LOCAL REAL SCALARS
REAL(DP):: ZDSCO2, ZDSH2O, ZFPPW, ZRTH, ZRTU, ZWH2O


!     ------------------------------------------------------------------

!*         1.     COMPUTES AMOUNTS OF ABSORBERS
!                 -----------------------------


IIND(1)=1
IIND(2)=2


!*         1.1    INITIALIZES QUANTITIES
!                 ----------------------


DO JL = KIDIA,KFDIA
  PUD(JL,1,KLEV+1)=0._DP
  PUD(JL,2,KLEV+1)=0._DP
  PUD(JL,3,KLEV+1)=0._DP
  PUD(JL,4,KLEV+1)=0._DP
  PUD(JL,5,KLEV+1)=0._DP
  PFACT(JL)= PRMU0(JL) * PSCT
  PRMU(JL)=SQRT(1224._DP* PRMU0(JL) * PRMU0(JL) + 1._DP) / 35._DP
  PSEC(JL)=1._DP/PRMU(JL)
ENDDO

!*          1.3    AMOUNTS OF ABSORBERS
!                  --------------------


DO JL= KIDIA,KFDIA
  ZUD(JL,1) = 0._DP
  ZUD(JL,2) = 0._DP
  ZO175(JL) = PPSOL(JL)** RPDU1
  ZO190(JL) = PPSOL(JL)** RPDH1
  ZSIGO(JL) = PPSOL(JL)
ENDDO

DO JK = 1 , KLEV
  JKP1 = JK + 1
  JKL = KLEV+1 - JK
  JKLP1 = JKL+1
  DO JL = KIDIA,KFDIA
    ZRTH=(RTH2O/PTAVE(JL,JK))**RTDH2O
    ZRTU=(RTUMG/PTAVE(JL,JK))**RTDUMG
    ZWH2O = MAX (PWV(JL,JKL) , EPSILON(1._DP) )
    ZSIGN(JL) = 100._DP * PPMB(JL,JKP1)
    PDSIG(JL,JK) = (ZSIGO(JL) - ZSIGN(JL))/PPSOL(JL)
    ZN175(JL) = ZSIGN(JL) ** RPDU1
    ZN190(JL) = ZSIGN(JL) ** RPDH1
    ZDSCO2 = ZO175(JL) - ZN175(JL)
    ZDSH2O = ZO190(JL) - ZN190(JL)
    PUD(JL,1,JK) = RPNH * ZDSH2O * ZWH2O  * ZRTH
    PUD(JL,2,JK) = RPNU * ZDSCO2 * PCARDI(JL,JKL) * ZRTU
    ZFPPW=1.6078_DP*ZWH2O/(1._DP+0.608_DP*ZWH2O)
    PUD(JL,4,JK)=PUD(JL,1,JK)*ZFPPW
    PUD(JL,5,JK)=PUD(JL,1,JK)*(1._DP-ZFPPW)
    ZUD(JL,1) = ZUD(JL,1) + PUD(JL,1,JK)
    ZUD(JL,2) = ZUD(JL,2) + PUD(JL,2,JK)
    ZSIGO(JL) = ZSIGN(JL)
    ZO175(JL) = ZN175(JL)
    ZO190(JL) = ZN190(JL)
  ENDDO
ENDDO


!*         1.4    COMPUTES CLEAR-SKY GREY ABSORPTION COEFFICIENTS
!                 -----------------------------------------------


DO JA = 1,2
  DO JL = KIDIA,KFDIA
    ZUD(JL,JA) = ZUD(JL,JA) * PSEC(JL)
  ENDDO
ENDDO

!! DO JNU= 2,NSW
DO JNU=4,NSW   !! new

  CALL SWTT1 ( KIDIA,KFDIA,KBDIM, JNU, 2, IIND &
   &, ZUD &
   &, ZR                            )

  DO JA = 1,2
    DO JL = KIDIA,KFDIA
      PAKI(JL,JA,JNU) = -LOG( ZR(JL,JA) ) / ZUD(JL,JA)
    ENDDO
  ENDDO
ENDDO


!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SWU
