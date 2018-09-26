#if defined (__SX__) || defined (ES) || defined(_UNICOSMP)
#  define UNROLLNGX
#endif
!******************************************************************************
!                                                                             *
!                  Optical depths developed for the                           *
!                                                                             *
!                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
!                                                                             *
!            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
!                        840 MEMORIAL DRIVE                                   *
!                        CAMBRIDGE, MA 02139                                  *
!                                                                             *
!                           ELI J. MLAWER                                     *
!                         STEVEN J. TAUBMAN                                   *
!                         SHEPARD A. CLOUGH                                   *
!                                                                             *
!                       email:  mlawer@aer.com                                *
!                                                                             *
!        The authors wish to acknowledge the contributions of the             *
!        following people:  Patrick D. Brown, Michael J. Iacono,              *
!        Ronald E. Farren, Luke Chen, Robert Bergstrom.                       *
!                                                                             *
!******************************************************************************
! Modified by:                                                                *
!      JJ Morcrette 980714 ECMWF      for use on ECMWF's Fujitsu VPP770       *
!         Reformatted for F90 by JJMorcrette, ECMWF                           * 
!         - replacing COMMONs by MODULEs                                      *
!         - changing labelled to unlabelled DO loops                          *
!         - creating set-up routines for all block data statements            *
!         - reorganizing the parameter statements                             * 
!         - passing KLEV as argument                                          *
!         - suppressing some equivalencing                                    *
!                                                                             *
!      D Salmond    9907   ECMWF      Speed-up modifications                  *
!      D Salmond    000515 ECMWF      Speed-up modifications                  *
!******************************************************************************
!     TAUMOL                                                                  *
!                                                                             *
!     This file contains the subroutines TAUGBn (where n goes from            *
!     1 to 16).  TAUGBn calculates the optical depths and Planck fractions    *
!     per g-value and layer for band n.                                       *
!                                                                             *
!  Output:  optical depths (unitless)                                         *
!           fractions needed to compute Planck functions at every layer       *
!               and g-value                                                   *
!                                                                             *
!     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
!     COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
!                                                                             *
!  Input                                                                      *
!                                                                             *
!     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
!     COMMON /PRECISE/  ONEMINUS                                              *
!     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
!    &                  PZ(0:MXLAY),TZ(0:MXLAY),TBOUND                        *
!     COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,                              *
!    &                  COLH2O(MXLAY),COLCO2(MXLAY),                          *
!    &                  COLO3(MXLAY),COLN2O(MXLAY),COLCH4(MXLAY),             *
!    &                  COLO2(MXLAY),CO2MULT(MXLAY)                           *
!     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
!    &                  FAC10(MXLAY),FAC11(MXLAY)                             *
!     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)                        *
!     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
!                                                                             *
!     Description:                                                            *
!     NG(IBAND) - number of g-values in band IBAND                            *
!     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
!                   atmospheres that are stored for band IBAND per            *
!                   pressure level and temperature.  Each of these            *
!                   atmospheres has different relative amounts of the         *
!                   key species for the band (i.e. different binary           *
!                   species parameters).                                      *
!     NSPB(IBAND) - same for upper atmosphere                                 *
!     ONEMINUS - since problems are caused in some cases by interpolation     *
!                parameters equal to or greater than 1, for these cases       *
!                these parameters are set to this value, slightly < 1.        *
!     PAVEL - layer pressures (mb)                                            *
!     TAVEL - layer temperatures (degrees K)                                  *
!     PZ - level pressures (mb)                                               *
!     TZ - level temperatures (degrees K)                                     *
!     LAYTROP - layer at which switch is made from one combination of         *
!               key species to another                                        *
!     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
!               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
!               respectively (molecules/cm**2)                                *
!     CO2MULT - for bands in which carbon dioxide is implemented as a         *
!               trace species, this is the factor used to multiply the        *
!               band's average CO2 absorption coefficient to get the added    *
!               contribution to the optical depth relative to 355 ppm.        *
!     FACij(LAY) - for layer LAY, these are factors that are needed to        *
!                  compute the interpolation factors that multiply the        *
!                  appropriate reference k-values.  A value of 0 (1) for      *
!                  i,j indicates that the corresponding factor multiplies     *
!                  reference k-value for the lower (higher) of the two        *
!                  appropriate temperatures, and altitudes, respectively.     *
!     JP - the index of the lower (in altitude) of the two appropriate        *
!          reference pressure levels needed for interpolation                 *
!     JT, JT1 - the indices of the lower of the two appropriate reference     *
!               temperatures needed for interpolation (for pressure           *
!               levels JP and JP+1, respectively)                             *
!     SELFFAC - scale factor needed to water vapor self-continuum, equals     *
!               (water vapor density)/(atmospheric density at 296K and        *
!               1013 mb)                                                      *
!     SELFFRAC - factor needed for temperature interpolation of reference     *
!                water vapor self-continuum data                              *
!     INDSELF - index of the lower of the two appropriate reference           *
!               temperatures needed for the self-continuum interpolation      *
!                                                                             *
!  Data input                                                                 *
!     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG) *
!        (note:  n is the band number)                                        *
!                                                                             *
!     Description:                                                            *
!     KA - k-values for low reference atmospheres (no water vapor             *
!          self-continuum) (units: cm**2/molecule)                            *
!     KB - k-values for high reference atmospheres (all sources)              *
!          (units: cm**2/molecule)                                            *
!     SELFREF - k-values for water vapor self-continuum for reference         *
!               atmospheres (used below LAYTROP)                              *
!               (units: cm**2/molecule)                                       *
!                                                                             *
!     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
!     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
!                                                                             *
!******************************************************************************


SUBROUTINE RRTM_TAUMOL1 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,&
  &COLH2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.
!     Revised by Michael J. Iacono, Atmospheric & Environmental Research.

!     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
 
! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,NG1
USE MO_RRTWN   , ONLY : NSPA   ,NSPB
USE MO_RRTA1   , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
             &FORREF   ,SELFREF 

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)
REAL(DP):: FORFAC(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)

INTEGER :: IND0(KBDIM),IND1(KBDIM),INDS(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, LAY, IPLON
INTEGER :: IXP, IXC0


!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  Below LAYTROP, the water vapor self-continuum 
!     is interpolated (in temperature) separately.  

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(1) + 1
      IND1(IPLON) =  (JP(IPLON,LAY)   *5+(JT1(IPLON,LAY)-1))*NSPA(1) + 1
      INDS(IPLON) = INDSELF(IPLON,LAY)
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG1
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=8
!OCL  UNROLL(8)
!DIR$ UNROLL 8
      DO IG = 1, NG1
#endif
        TAU (IPLON,IG,LAY) = COLH2O(IPLON,LAY) *              &
             (FAC00(IPLON,LAY) * ABSA(IND0(IPLON)  ,IG) +     &
              FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+1,IG) +     &
              FAC01(IPLON,LAY) * ABSA(IND1(IPLON)  ,IG) +     &
              FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+1,IG) +     &
            SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) +   &
           SELFFRAC(IPLON,LAY) *                              &
           (SELFREF(INDS(IPLON)+1,IG) - SELFREF(INDS(IPLON),IG)))  &
             + FORFAC(IPLON,LAY) * FORREF(IG) )               &
             + TAUAERL(IPLON,LAY,1)
        PFRAC(IPLON,IG,LAY) = FRACREFA(IG)
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      IND0(IPLON)  = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(1) + 1
      IND1(IPLON)  = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(1) + 1
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG1
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=8
!OCL  UNROLL(8)
!DIR$ UNROLL 8
      DO IG = 1, NG1
#endif
        TAU (IPLON,IG,LAY) = COLH2O(IPLON,LAY) *          &
             (FAC00(IPLON,LAY) * ABSB(IND0(IPLON)  ,IG) + &
              FAC10(IPLON,LAY) * ABSB(IND0(IPLON)+1,IG) + &
              FAC01(IPLON,LAY) * ABSB(IND1(IPLON)  ,IG) + &
              FAC11(IPLON,LAY) * ABSB(IND1(IPLON)+1,IG)   &
           + FORFAC(IPLON,LAY) * FORREF(IG) )             &
          + TAUAERL(IPLON,LAY,1)
        PFRAC(IPLON,IG,LAY) = FRACREFB(IG)
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL1
!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL2 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,COLDRY,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,&
  &COLH2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up
!     JJMorcrette 2000-07-14 bugfix


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,NG2   ,NGS1
USE MO_RRTWN   , ONLY : NSPA   ,NSPB
USE MO_RRTA2   , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
             &FORREF   ,SELFREF , REFPARAM
USE MO_RRTBG2  , ONLY : CORR1  ,CORR2

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

REAL(DP):: COLDRY(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)
REAL(DP):: FORFAC(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


REAL(DP):: FRACINT(KBDIM)

INTEGER :: INDEX(KBDIM)

INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
INTEGER :: IFPX(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IFP, IFRAC, IG, JFRAC, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: FP, H2OPARAM, WATER


!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  Below LAYTROP, the water vapor self-continuum is 
!     interpolated (in temperature) separately.

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!CDIR NODEP
!OCL NOVREC
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      WATER = 1.E20_dp * COLH2O(IPLON,LAY) / COLDRY(IPLON,LAY)
      H2OPARAM = WATER / (WATER +.002_dp)

      IF (H2OPARAM >= REFPARAM(2)) THEN
        INDEX(IPLON) = 2
      ELSE
!DIR$ UNROLL 11
!CDIR UNROLL=11
!OCL  UNROLL(11)
        DO JFRAC = 2, 12
          IF (H2OPARAM < REFPARAM(JFRAC)) THEN
            INDEX(IPLON) = JFRAC + 1
          END IF
        ENDDO
      ENDIF

      !---- JJM_000714
      IFRAC = INDEX(IPLON)
      FRACINT(IPLON) = (H2OPARAM-REFPARAM(IFRAC)) / &
           (REFPARAM(IFRAC-1)-REFPARAM(IFRAC))
    ENDDO

!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)

      FP = FAC11(IPLON,LAY) + FAC01(IPLON,LAY)
      IFP = 2.E2_dp*FP + 0.5_dp

      !---MI 981104        
      IFPX(IPLON) = MAX(0,IFP)

      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(2) + 1
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(2) + 1
      INDS(IPLON) = INDSELF(IPLON,LAY)
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG2
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      IFRAC = INDEX(IPLON)
      IFP   = IFPX(IPLON)
#ifdef UNROLLNGX
!DIR$ UNROLL 14
!CDIR UNROLL=14
!OCL  UNROLL(14)
      DO IG = 1, NG2
#endif
        TAU (IPLON,NGS1+IG,LAY) = COLH2O(IPLON,LAY) *                     &
             (CORR2(IFP) * (FAC00(IPLON,LAY) * ABSA(IND0(IPLON)  ,IG)  +  &
                            FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+1,IG)) +  &
              CORR1(IFP) * (FAC01(IPLON,LAY) * ABSA(IND1(IPLON)  ,IG)  +  &
                            FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+1,IG)) +  &
                          SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) + &
                         SELFFRAC(IPLON,LAY) *                            &
                         (SELFREF(INDS(IPLON)+1,IG) - SELFREF(INDS(IPLON),IG))) &
                       + FORFAC(IPLON,LAY) * FORREF(IG) )                 &
                       + TAUAERL(IPLON,LAY,2)
        PFRAC(IPLON,NGS1+IG,LAY) = FRACREFA(IG,IFRAC) + FRACINT(IPLON) * &
             (FRACREFA(IG,IFRAC-1) - FRACREFA(IG,IFRAC))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      FP = FAC11(IPLON,LAY) + FAC01(IPLON,LAY)
      IFP = 2.E2_dp*FP + 0.5_dp

      !---MI 981104        
      IFPX(IPLON) = MAX(0,IFP)

      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(2) + 1
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(2) + 1
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG2
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      IFP   = IFPX(IPLON)
#ifdef UNROLLNGX
!DIR$ UNROLL 14
!CDIR UNROLL=14
!OCL  UNROLL(14)
      DO IG = 1, NG2
#endif
        TAU (IPLON,NGS1+IG,LAY) = COLH2O(IPLON,LAY) *                    &
             (CORR2(IFP) * (FAC00(IPLON,LAY) * ABSB(IND0(IPLON)  ,IG)  + &
                            FAC10(IPLON,LAY) * ABSB(IND0(IPLON)+1,IG)) + &
              CORR1(IFP) * (FAC01(IPLON,LAY) * ABSB(IND1(IPLON)  ,IG)  + &
                            FAC11(IPLON,LAY) * ABSB(IND1(IPLON)+1,IG))   &
                         + FORFAC(IPLON,LAY) * FORREF(IG) )              &
                        + TAUAERL(IPLON,LAY,2)
        PFRAC(IPLON,NGS1+IG,LAY) = FRACREFB(IG)
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL2
!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL3 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCO2,COLN2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)

! Modifications
!
!     D Salmond 2000-05-15 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,NG3   ,NGS2
USE MO_RRTWN   , ONLY : NSPA   ,NSPB
USE MO_RRTA3   , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
             &FORREF   ,SELFREF , ABSN2OA ,&
             &ABSN2OB  ,ETAREF ,H2OREF ,N2OREF  , CO2REF  ,&
             &STRRAT

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)
REAL(DP):: FORFAC(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLCO2(KBDIM,KLEV)
REAL(DP):: COLN2O(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM), N2OMULT(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, NS, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: COLREF1, COLREF2, CURRN2O,                      &
           FP, FS, RATIO, SPECMULT, SPECPARM, WCOMB1,      &
           WCOMB2


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!CDIR NODEP
!DIR$ CONCURRENT
!OCL NOVREC
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT*COLCO2(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 8._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IF (JS  ==  8) THEN
        IF (FS  >=  0.9_dp) THEN
          JS = 9
          FS = 10._dp * (FS - 0.9_dp)
        ELSE
          FS = FS / 0.9_dp
        ENDIF
      ENDIF

      NS = JS + INT(FS + 0.5_dp)
      FP = FAC01(IPLON,LAY) + FAC11(IPLON,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(3) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(3) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      COLREF1 = N2OREF(JP(IPLON,LAY))
      COLREF2 = N2OREF(JP(IPLON,LAY)+1)
      IF (NS  ==  10) THEN
        WCOMB1 = 1.0_dp/H2OREF(JP(IPLON,LAY))
        WCOMB2 = 1.0_dp/H2OREF(JP(IPLON,LAY)+1)
      ELSE
        WCOMB1 = (1.0_dp-ETAREF(NS))/(STRRAT * CO2REF(JP(IPLON,LAY)))
        WCOMB2 = (1.0_dp-ETAREF(NS))/(STRRAT * CO2REF(JP(IPLON,LAY)+1))
      ENDIF
      RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
      CURRN2O = SPECCOMB(IPLON) * RATIO
      N2OMULT(IPLON) = COLN2O(IPLON,LAY) - CURRN2O
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG3
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!DIR$ UNROLL 16
!CDIR UNROLL=16
!OCL  UNROLL(16)
      DO IG = 1, NG3
#endif
        TAU (IPLON,NGS2+IG,LAY) = SPECCOMB(IPLON) *                      &
             ((1. - FS) *(FAC00(IPLON,LAY) * ABSA(IND0(IPLON)   ,IG) +   &
                          FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+10,IG) +   &
                          FAC01(IPLON,LAY) * ABSA(IND1(IPLON)   ,IG) +   &
                          FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+10,IG))+   &
                    FS * (FAC00(IPLON,LAY) * ABSA(IND0(IPLON)+ 1,IG) +   &
                          FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+11,IG) +   &
                          FAC01(IPLON,LAY) * ABSA(IND1(IPLON)+ 1,IG) +   &
                          FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+11,IG))) + &
                          COLH2O(IPLON,LAY) *                            &
              SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) +            &
              SELFFRAC(IPLON,LAY) *                                      &
              (SELFREF(INDS(IPLON)+1,IG) - SELFREF(INDS(IPLON),IG))      &
              + FORFAC(IPLON,LAY) * FORREF(IG) )                         &
              + N2OMULT(IPLON) * ABSN2OA(IG)                             &
              + TAUAERL(IPLON,LAY,3)
        PFRAC(IPLON,NGS2+IG,LAY) = FRACREFA(IG,JS) + FS *         &
              (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
!CDIR NODEP
!DIR$ CONCURRENT
!OCL NOVREC
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT*COLCO2(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 4._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      NS = JS + INT(FS + 0.5_dp)
      FP = FAC01(IPLON,LAY) + FAC11(IPLON,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(3) + JS
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(3) + JS
      COLREF1 = N2OREF(JP(IPLON,LAY))
      COLREF2 = N2OREF(JP(IPLON,LAY)+1)
      IF (NS  ==  5) THEN
        WCOMB1 = 1.0_dp/H2OREF(JP(IPLON,LAY))
        WCOMB2 = 1.0_dp/H2OREF(JP(IPLON,LAY)+1)
      ELSE
        WCOMB1 = (1.0_dp-ETAREF(NS))/(STRRAT * CO2REF(JP(IPLON,LAY)))
        WCOMB2 = (1.0_dp-ETAREF(NS))/(STRRAT * CO2REF(JP(IPLON,LAY)+1))
      ENDIF
      RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
      CURRN2O = SPECCOMB(IPLON) * RATIO
      N2OMULT(IPLON) = COLN2O(IPLON,LAY) - CURRN2O
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG3
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!DIR$ UNROLL 16
!CDIR UNROLL=16
!OCL  UNROLL(16)
      DO IG = 1, NG3
#endif
        TAU (IPLON,NGS2+IG,LAY) = SPECCOMB(IPLON) *                    &
             ((1. - FS) *(FAC00(IPLON,LAY) * ABSB(IND0(IPLON)  ,IG) +  &
                          FAC10(IPLON,LAY) * ABSB(IND0(IPLON)+5,IG) +  &
                          FAC01(IPLON,LAY) * ABSB(IND1(IPLON)  ,IG) +  &
                          FAC11(IPLON,LAY) * ABSB(IND1(IPLON)+5,IG))+  &
                    FS * (FAC01(IPLON,LAY) * ABSB(IND1(IPLON)+1,IG) +  &
                          FAC10(IPLON,LAY) * ABSB(IND0(IPLON)+6,IG) +  &
                          FAC00(IPLON,LAY) * ABSB(IND0(IPLON)+1,IG) +  &
                          FAC11(IPLON,LAY) * ABSB(IND1(IPLON)+6,IG)))  &
                       + COLH2O(IPLON,LAY)*FORFAC(IPLON,LAY)*FORREF(IG)  &
                       + N2OMULT(IPLON) * ABSN2OB(IG)                  &
                       + TAUAERL(IPLON,LAY,3)
        PFRAC(IPLON,NGS2+IG,LAY) = FRACREFB(IG,JS) + FS *       &
             (FRACREFB(IG,JS+1) - FRACREFB(IG,JS))
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL3
!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL4 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCO2,COLO3,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)

! Modifications
!
!     D Salmond 2000-05-15 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,NG4   ,NGS3
USE MO_RRTWN   , ONLY : NSPA   ,NSPB
USE MO_RRTA4   , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
             &SELFREF,STRRAT1,STRRAT2

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLCO2(KBDIM,KLEV)
REAL(DP):: COLO3(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: FS, SPECMULT, SPECPARM


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately. 
 
  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT1*COLCO2(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 8._dp * (SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(4) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(4) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG4
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!DIR$ UNROLL 14
!CDIR UNROLL=14
!OCL  UNROLL(14)
      DO IG = 1, NG4
#endif
        TAU (IPLON,NGS3+IG,LAY) = SPECCOMB(IPLON) *                     &
             ((1. - FS)*(FAC00(IPLON,LAY) * ABSA(IND0(IPLON)   ,IG)  +  &
                         FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+ 9,IG)  +  &
                         FAC01(IPLON,LAY) * ABSA(IND1(IPLON)   ,IG)  +  &
                         FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+ 9,IG)) +  &
                   FS * (FAC01(IPLON,LAY) * ABSA(IND1(IPLON)+ 1,IG)  +  &
                         FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+10,IG)  +  &
                         FAC00(IPLON,LAY) * ABSA(IND0(IPLON)+ 1,IG)  +  &
                         FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+10,IG))) + &
                        COLH2O(IPLON,LAY) *                             &
                       SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) +  &
                       SELFFRAC(IPLON,LAY) *                            &
                       (SELFREF(INDS(IPLON)+1,IG) - SELFREF(INDS(IPLON),IG)))  &
                      + TAUAERL(IPLON,LAY,4)
        PFRAC(IPLON,NGS3+IG,LAY) = FRACREFA(IG,JS) + FS *        &
             (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      SPECCOMB(IPLON) = COLO3(IPLON,LAY) + STRRAT2*COLCO2(IPLON,LAY)
      SPECPARM = COLO3(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 4._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IF (JS  >  1) THEN
        JS = JS + 1
        ELSEIF (FS  >=  0.0024_dp) THEN
        JS = 2
        FS = (FS - 0.0024_dp)/0.9976_dp
      ELSE
        JS = 1
        FS = FS/0.0024_dp
      ENDIF
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(4) + JS
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(4) + JS
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG4
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!DIR$ UNROLL 14
!CDIR UNROLL=14
!OCL  UNROLL(14)
      DO IG = 1, NG4
#endif
        TAU (IPLON,NGS3+IG,LAY) = SPECCOMB(IPLON) *                   &
             ((1. - FS)*(FAC00(IPLON,LAY) * ABSB(IND0(IPLON)  ,IG) +  &
                         FAC10(IPLON,LAY) * ABSB(IND0(IPLON)+6,IG) +  &
                         FAC01(IPLON,LAY) * ABSB(IND1(IPLON)  ,IG) +  &
                         FAC11(IPLON,LAY) * ABSB(IND1(IPLON)+6,IG)) + &
                   FS * (FAC10(IPLON,LAY) * ABSB(IND0(IPLON)+7,IG) +  &
                         FAC01(IPLON,LAY) * ABSB(IND1(IPLON)+1,IG) +  &
                         FAC00(IPLON,LAY) * ABSB(IND0(IPLON)+1,IG) +  &
                         FAC11(IPLON,LAY) * ABSB(IND1(IPLON)+7,IG)))  &
                     + TAUAERL(IPLON,LAY,4)
        PFRAC(IPLON,NGS3+IG,LAY) = FRACREFB(IG,JS) + FS *      &
              (FRACREFB(IG,JS+1) - FRACREFB(IG,JS))
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL4
!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL5 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,WX,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCO2, COLO3,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)

! Modifications
!
!     D Salmond 2000-05-15 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,JPXSEC ,NG5   ,NGS4
USE MO_RRTWN   , ONLY : NSPA   ,NSPB
USE MO_RRTA5   , ONLY : ABSA   ,ABSB   ,CCL4   , FRACREFA, FRACREFB,&
             &SELFREF,STRRAT1,STRRAT2

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

REAL(DP):: WX(KBDIM,JPXSEC,KLEV)
!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLCO2(KBDIM,KLEV)
REAL(DP):: COLO3(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: FS, SPECMULT, SPECPARM


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT1*COLCO2(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 8._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(5) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(5) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG5
#endif
!!! !OCL NOVREC
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!DIR$ UNROLL 16
!CDIR UNROLL=16
!OCL  UNROLL(16)
      DO IG = 1, NG5
#endif
        TAU (IPLON,NGS4+IG,LAY) = SPECCOMB(IPLON) *                   &
          ((1. - FS)*(FAC00(IPLON,LAY) * ABSA(IND0(IPLON),IG) +       &
                      FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+9,IG) +     &
                      FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+9,IG) +     &
                      FAC01(IPLON,LAY) * ABSA(IND1(IPLON),IG)) +      &
                 FS* (FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+10,IG) +    &
                      FAC01(IPLON,LAY) * ABSA(IND1(IPLON)+1,IG) +     &
                      FAC00(IPLON,LAY) * ABSA(IND0(IPLON)+1,IG) +     &
                      FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+10,IG)) ) + &
                   COLH2O(IPLON,LAY) *                                &
                   SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) +    &
                   SELFFRAC(IPLON,LAY) *                              &
                  (SELFREF(INDS(IPLON)+1,IG) - SELFREF(INDS(IPLON),IG)))     &
                + WX(IPLON,1,LAY) * CCL4(IG)                          &
                + TAUAERL(IPLON,LAY,5)
        PFRAC(IPLON,NGS4+IG,LAY) = FRACREFA(IG,JS) + FS * &
             (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      SPECCOMB(IPLON) = COLO3(IPLON,LAY) + STRRAT2*COLCO2(IPLON,LAY)
      SPECPARM = COLO3(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 4._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(5) + JS
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(5) + JS
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG5
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!DIR$ UNROLL 16
!CDIR UNROLL=16
!OCL  UNROLL(16)
      DO IG = 1, NG5
#endif
        TAU (IPLON,NGS4+IG,LAY) = SPECCOMB(IPLON) *               &
        ((1. - FS)*(FAC00(IPLON,LAY) * ABSB(IND0(IPLON),IG) +     &
                    FAC10(IPLON,LAY) * ABSB(IND0(IPLON)+5,IG) +   &
                    FAC01(IPLON,LAY) * ABSB(IND1(IPLON),IG) +     &
                    FAC11(IPLON,LAY) * ABSB(IND1(IPLON)+5,IG) ) + &
              FS * (FAC01(IPLON,LAY) * ABSB(IND1(IPLON)+1,IG) +   &
                    FAC10(IPLON,LAY) * ABSB(IND0(IPLON)+6,IG) +   &
                    FAC00(IPLON,LAY) * ABSB(IND0(IPLON)+1,IG) +   &
                    FAC11(IPLON,LAY) * ABSB(IND1(IPLON)+6,IG)))   &
             + WX(IPLON,1,LAY) * CCL4(IG)                         &
             + TAUAERL(IPLON,LAY,5)
        PFRAC(IPLON,NGS4+IG,LAY) = FRACREFB(IG,JS) + FS *  &
             (FRACREFB(IG,JS+1) - FRACREFB(IG,JS))
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL5
!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL6 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,WX,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
  &COLH2O,CO2MULT,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,JPXSEC ,NG6   ,NGS5
USE MO_RRTWN   , ONLY : NSPA
USE MO_RRTA6   , ONLY : ABSA   ,ABSCO2 ,CFC11ADJ , CFC12  ,&
             &FRACREFA, SELFREF

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

REAL(DP):: WX(KBDIM,JPXSEC,KLEV)
!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: CO2MULT(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, LAY, IPLON
INTEGER :: IXP, IXC0


!     Compute the optical depth by interpolating in ln(pressure) and
!     temperature. The water vapor self-continuum is interpolated
!     (in temperature) separately.  

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(6) + 1
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(6) + 1
      INDS(IPLON) = INDSELF(IPLON,LAY)
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG6
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 8
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG6
#endif
        TAU (IPLON,NGS5+IG,LAY) = COLH2O(IPLON,LAY) *        &
             (FAC00(IPLON,LAY) * ABSA(IND0(IPLON)  ,IG) +    &
              FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+1,IG) +    &
              FAC01(IPLON,LAY) * ABSA(IND1(IPLON)  ,IG) +    &
              FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+1,IG) +    &
             SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) + &
             SELFFRAC(IPLON,LAY)*                            &
             (SELFREF(INDS(IPLON)+1,IG)-SELFREF(INDS(IPLON),IG))))  &
             + WX(IPLON,2,LAY) * CFC11ADJ(IG)                &
             + WX(IPLON,3,LAY) * CFC12(IG)                   &
             + CO2MULT(IPLON,LAY) * ABSCO2(IG)               &
             + TAUAERL(IPLON,LAY,6)
        PFRAC(IPLON,NGS5+IG,LAY) = FRACREFA(IG)
      ENDDO
    ENDDO

    !     Nothing important goes on above LAYTROP in this band.
    IXC0 = KPROMA - IXC0
#ifndef UNROLLNGX
    DO IG = 1, NG6
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 8
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG6
#endif
        TAU (IPLON,NGS5+IG,LAY) = 0.0_dp      &
             + WX(IPLON,2,LAY) * CFC11ADJ(IG) &
             + WX(IPLON,3,LAY) * CFC12(IG)    &
             + TAUAERL(IPLON,LAY,6)
        PFRAC(IPLON,NGS5+IG,LAY) = FRACREFA(IG)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL6
!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL7 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLO3,CO2MULT,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)

! Modifications
!
!     D Salmond 2000-05-15 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,NG7   ,NGS6
USE MO_RRTWN   , ONLY : NSPA   ,NSPB
USE MO_RRTA7   , ONLY : ABSA   ,ABSB   ,ABSCO2 ,FRACREFA ,FRACREFB,&
             &SELFREF,STRRAT

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLO3(KBDIM,KLEV)
REAL(DP):: CO2MULT(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: FS, SPECMULT, SPECPARM


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.
  
  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT*COLO3(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 8._dp*SPECPARM
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(7) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(7) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG7
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!DIR$ UNROLL 12
!CDIR UNROLL=12
!OCL  UNROLL(12)
      DO IG = 1, NG7
#endif
        TAU (IPLON,NGS6+IG,LAY) = SPECCOMB(IPLON) *                     &
             ((1. - FS)*(FAC00(IPLON,LAY) * ABSA(IND0(IPLON),IG) +      &
                         FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+9,IG) +    &
                         FAC01(IPLON,LAY) * ABSA(IND1(IPLON),IG) +      &
                         FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+9,IG) )+   &
                   FS * (FAC01(IPLON,LAY) * ABSA(IND1(IPLON)+1,IG) +    &
                         FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+10,IG) +   &
                         FAC00(IPLON,LAY) * ABSA(IND0(IPLON)+1,IG) +    &
                         FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+10,IG))) + &
                        COLH2O(IPLON,LAY) *                             &
                       SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) +  &
                      SELFFRAC(IPLON,LAY) *                             &
                     (SELFREF(INDS(IPLON)+1,IG) - SELFREF(INDS(IPLON),IG)))    &
                    + CO2MULT(IPLON,LAY) * ABSCO2(IG)                   &
                    + TAUAERL(IPLON,LAY,7)
        PFRAC(IPLON,NGS6+IG,LAY) = FRACREFA(IG,JS) + FS * &
             (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(7) + 1
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(7) + 1
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG7
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 12
!CDIR UNROLL=12
!OCL  UNROLL(12)
      DO IG = 1, NG7
#endif
        TAU (IPLON,NGS6+IG,LAY) = COLO3(IPLON,LAY) *        &
             (FAC00(IPLON,LAY) * ABSB(IND0(IPLON)  ,IG) +   &
              FAC10(IPLON,LAY) * ABSB(IND0(IPLON)+1,IG) +   &
              FAC01(IPLON,LAY) * ABSB(IND1(IPLON)  ,IG) +   &
              FAC11(IPLON,LAY) * ABSB(IND1(IPLON)+1,IG))    &
             + CO2MULT(IPLON,LAY) * ABSCO2(IG)              &
             + TAUAERL(IPLON,LAY,7)
        PFRAC(IPLON,NGS6+IG,LAY) = FRACREFB(IG)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL7
!*******************************************************************************
SUBROUTINE RRTM_TAUMOL8 (KPROMA,KBDIM,KLEV,IXS,IXLOS,IXHIGS,TAU,WX,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
  &COLH2O,COLO3,COLN2O,CO2MULT,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,JPXSEC ,NG8   ,NGS7
USE MO_RRTWN   , ONLY : NSPA   ,NSPB
USE MO_RRTA8   , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
             &SELFREF ,ABSCO2A ,ABSCO2B,&
             &ABSN2OA , ABSN2OB,CFC12  ,CFC22ADJ, H2OREF  ,&
             &N2OREF  , O3REF

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXS(KLEV), IXLOS(KBDIM,KLEV), IXHIGS(KBDIM,KLEV)

REAL(DP):: WX(KBDIM,JPXSEC,KLEV)
!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLO3(KBDIM,KLEV)
REAL(DP):: COLN2O(KBDIM,KLEV)
REAL(DP):: CO2MULT(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: COLREF1, COLREF2, CURRN2O, FP, RATIO, WCOMB1, WCOMB2
REAL(DP):: N2OMULT(KBDIM)


!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  

  DO LAY = 1, KLEV
    IXC0 = IXS(LAY)

!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOS(IXP,LAY)
      FP = FAC01(IPLON,LAY) + FAC11(IPLON,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT(IPLON,LAY)-1))*NSPA(8) + 1
      IND1(IPLON) = (JP(IPLON,LAY)*5+(JT1(IPLON,LAY)-1))*NSPA(8) + 1
      INDS(IPLON) = INDSELF(IPLON,LAY)
      COLREF1 = N2OREF(JP(IPLON,LAY))
      COLREF2 = N2OREF(JP(IPLON,LAY)+1)
      WCOMB1 = 1.0_dp/H2OREF(JP(IPLON,LAY))
      WCOMB2 = 1.0_dp/H2OREF(JP(IPLON,LAY)+1)
      RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
      CURRN2O = COLH2O(IPLON,LAY) * RATIO
      N2OMULT(IPLON) = COLN2O(IPLON,LAY) - CURRN2O
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG8
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOS(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 8
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG8
#endif
        TAU (IPLON,NGS7+IG,LAY) = COLH2O(IPLON,LAY) *         &
             (FAC00(IPLON,LAY) * ABSA(IND0(IPLON)  ,IG) +     &
              FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+1,IG) +     &
              FAC01(IPLON,LAY) * ABSA(IND1(IPLON)  ,IG) +     &
              FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+1,IG) +     &
             SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) +  &
             SELFFRAC(IPLON,LAY) *                            &
             (SELFREF(INDS(IPLON)+1,IG) - SELFREF(INDS(IPLON),IG)))) &
             + WX(IPLON,3,LAY) * CFC12(IG)                    &
             + WX(IPLON,4,LAY) * CFC22ADJ(IG)                 &
             + CO2MULT(IPLON,LAY) * ABSCO2A(IG)               &
             + N2OMULT(IPLON) * ABSN2OA(IG)                   &
             + TAUAERL(IPLON,LAY,8)
        PFRAC(IPLON,NGS7+IG,LAY) = FRACREFA(IG)
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGS(IXP,LAY)
      FP = FAC01(IPLON,LAY) + FAC11(IPLON,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-7)*5+(JT (IPLON,LAY)-1))*NSPB(8) + 1
      IND1(IPLON) = ((JP(IPLON,LAY)-6)*5+(JT1(IPLON,LAY)-1))*NSPB(8) + 1
      COLREF1 = N2OREF(JP(IPLON,LAY))
      COLREF2 = N2OREF(JP(IPLON,LAY)+1)
      WCOMB1 = 1.0_dp/O3REF(JP(IPLON,LAY))
      WCOMB2 = 1.0_dp/O3REF(JP(IPLON,LAY)+1)
      RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
      CURRN2O = COLO3(IPLON,LAY) * RATIO
      N2OMULT(IPLON) = COLN2O(IPLON,LAY) - CURRN2O
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG8
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGS(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 8
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG8
#endif
        TAU (IPLON,NGS7+IG,LAY) = COLO3(IPLON,LAY) *        &
             (FAC00(IPLON,LAY) * ABSB(IND0(IPLON)  ,IG) +   &
              FAC10(IPLON,LAY) * ABSB(IND0(IPLON)+1,IG) +   &
              FAC01(IPLON,LAY) * ABSB(IND1(IPLON)  ,IG) +   &
              FAC11(IPLON,LAY) * ABSB(IND1(IPLON)+1,IG))    &
             + WX(IPLON,3,LAY) * CFC12(IG)                  &
             + WX(IPLON,4,LAY) * CFC22ADJ(IG)               &
             + CO2MULT(IPLON,LAY) * ABSCO2B(IG)             &
             + N2OMULT(IPLON) * ABSN2OB(IG)                 &
             + TAUAERL(IPLON,LAY,8)
        PFRAC(IPLON,NGS7+IG,LAY) = FRACREFB(IG)
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL8
!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL9 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLN2O,COLCH4,LAYSWTCH,LAYLOW,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,NG9   ,NGS8
USE MO_RRTWN   , ONLY : NSPA   ,NSPB
USE MO_RRTA9   , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
             &SELFREF , ABSN2O ,CH4REF ,&
             &ETAREF  , H2OREF ,N2OREF ,STRRAT

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLN2O(KBDIM,KLEV)
REAL(DP):: COLCH4(KBDIM,KLEV)
INTEGER :: LAYSWTCH(KBDIM)
INTEGER :: LAYLOW(KBDIM)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM), N2OMULT(KBDIM), FFRAC(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, IOFF(KBDIM), JS, LAY, NS, IPLON
INTEGER :: IXP, IXC0
INTEGER :: JFRAC

!     LOCAL REAL SCALARS
REAL(DP):: COLREF1, COLREF2, CURRN2O, &
           FP, FS, RATIO, SPECMULT, SPECPARM, WCOMB1, &
           WCOMB2

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.
  

  IOFF(:) = 0

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!CDIR NODEP
!DIR$ CONCURRENT
!OCL NOVREC
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT*COLCH4(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 8._dp * (SPECPARM)
      JS = 1 + INT(SPECMULT)
      JFRAC = JS
      FS = MOD(SPECMULT,1.0_dp)
      FFRAC(IPLON) = FS
      IF (JS  ==  8) THEN
        IF (FS .LE. 0.68_dp) THEN
          FS = FS/0.68_dp
        ELSEIF (FS  <=  0.92_dp) THEN
          JS = JS + 1
          FS = (FS-0.68_dp)/0.24_dp
        ELSE
          JS = JS + 2
          FS = (FS-0.92_dp)/0.08_dp
        ENDIF
      ELSEIF (JS == 9) THEN
        JS = 10
        FS = 1.0_dp
        JFRAC = 8
        FFRAC(IPLON) = 1.0_dp
      ENDIF
      FP = FAC01(IPLON,LAY) + FAC11(IPLON,LAY)
      NS = JS + INT(FS + 0.5_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(9) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(9) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      IF (LAY  ==  LAYLOW(IPLON)) IOFF(IPLON) = NG9
      IF (LAY  ==  LAYSWTCH(IPLON)) IOFF(IPLON) = 2*NG9
      COLREF1 = N2OREF(JP(IPLON,LAY))
      COLREF2 = N2OREF(JP(IPLON,LAY)+1)
      IF (NS  ==  11) THEN
        WCOMB1 = 1.0_dp/H2OREF(JP(IPLON,LAY))
        WCOMB2 = 1.0_dp/H2OREF(JP(IPLON,LAY)+1)
      ELSE
        WCOMB1 = (1.0_dp-ETAREF(NS))/(STRRAT * CH4REF(JP(IPLON,LAY)))
        WCOMB2 = (1.0_dp-ETAREF(NS))/(STRRAT * CH4REF(JP(IPLON,LAY)+1))
      ENDIF
      RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
      CURRN2O = SPECCOMB(IPLON) * RATIO
      N2OMULT(IPLON) = COLN2O(IPLON,LAY) - CURRN2O
      ZFS(IPLON) = FS
!    ENDDO

!#ifndef UNROLLNGX
!    DO IG = 1, NG9
!#endif
!    DO IXP = 1, IXC0
!      IPLON = IXLOW(IXP,LAY)
!      FS = ZFS(IPLON)
!#ifdef UNROLLNGX
!DIR$ UNROLL 12
!CDIR UNROLL=12
!OCL  UNROLL(12)
      DO IG = 1, NG9
!#endif
        TAU (IPLON,NGS8+IG,LAY) = SPECCOMB(IPLON) *                       &
           ((1. - FS)*(FAC00(IPLON,LAY) * ABSA(IND0(IPLON)   ,IG) +       &
                       FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+11,IG) +       &
                       FAC01(IPLON,LAY) * ABSA(IND1(IPLON)   ,IG) +       &
                       FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+11,IG)) +      &
               FS* (   FAC00(IPLON,LAY) * ABSA(IND0(IPLON)+ 1,IG) +       &
                       FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+12,IG) +       &
                       FAC01(IPLON,LAY) * ABSA(IND1(IPLON)+ 1,IG) +       &
                       FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+12,IG))) +     &
                      COLH2O(IPLON,LAY) *                                 &
                     SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) +      &
                    SELFFRAC(IPLON,LAY) *                                 &
                   (SELFREF(INDS(IPLON)+1,IG) - SELFREF(INDS(IPLON),IG))) &
                 + N2OMULT(IPLON) * ABSN2O(IG+IOFF(IPLON))                &
                 + TAUAERL(IPLON,LAY,9)
        PFRAC(IPLON,NGS8+IG,LAY) = FRACREFA(IG,JFRAC) + FFRAC(IPLON) *&
             (FRACREFA(IG,JFRAC+1) - FRACREFA(IG,JFRAC))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(9) + 1
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(9) + 1
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG9
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 12
!CDIR UNROLL=12
!OCL  UNROLL(12)
      DO IG = 1, NG9
#endif
        TAU (IPLON,NGS8+IG,LAY) = COLCH4(IPLON,LAY) *     &
             (FAC00(IPLON,LAY) * ABSB(IND0(IPLON)  ,IG) + &
              FAC10(IPLON,LAY) * ABSB(IND0(IPLON)+1,IG) + &
              FAC01(IPLON,LAY) * ABSB(IND1(IPLON)  ,IG) + &
              FAC11(IPLON,LAY) * ABSB(IND1(IPLON)+1,IG))  &
             + TAUAERL(IPLON,LAY,9)
        PFRAC(IPLON,NGS8+IG,LAY) = FRACREFB(IG)
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL9
!*******************************************************************************
SUBROUTINE RRTM_TAUMOL10 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
  &COLH2O,PFRAC)

!     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,NG10   ,NGS9
USE MO_RRTWN   , ONLY : NSPA   ,NSPB
USE MO_RRTA10  , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IND0(KBDIM), IND1(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, LAY, IPLON
INTEGER :: IXP, IXC0


!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(10) + 1
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(10) + 1
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG10
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 6
!CDIR UNROLL=6
!OCL  UNROLL(6)
      DO IG = 1, NG10
#endif
        TAU (IPLON,NGS9+IG,LAY) = COLH2O(IPLON,LAY) *     &
             (FAC00(IPLON,LAY) * ABSA(IND0(IPLON)  ,IG) + &
              FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+1,IG) + &
              FAC01(IPLON,LAY) * ABSA(IND1(IPLON)  ,IG) + &
              FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+1,IG))  &
             + TAUAERL(IPLON,LAY,10)
        PFRAC(IPLON,NGS9+IG,LAY) = FRACREFA(IG)
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(10) + 1
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(10) + 1
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG10
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 6
!CDIR UNROLL=6
!OCL  UNROLL(6)
      DO IG = 1, NG10
#endif
        TAU (IPLON,NGS9+IG,LAY) = COLH2O(IPLON,LAY) *     &
             (FAC00(IPLON,LAY) * ABSB(IND0(IPLON)  ,IG) + &
              FAC10(IPLON,LAY) * ABSB(IND0(IPLON)+1,IG) + &
              FAC01(IPLON,LAY) * ABSB(IND1(IPLON)  ,IG) + &
              FAC11(IPLON,LAY) * ABSB(IND1(IPLON)+1,IG))  &
             + TAUAERL(IPLON,LAY,10)
        PFRAC(IPLON,NGS9+IG,LAY) = FRACREFB(IG)
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL10
!******************************************************************************
SUBROUTINE RRTM_TAUMOL11 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
  &COLH2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,NG11  ,NGS10
USE MO_RRTWN   , ONLY : NSPA   ,NSPB
USE MO_RRTA11  , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
            &SELFREF

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, LAY, IPLON
INTEGER :: IXP, IXC0


!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  Below LAYTROP, the water vapor self-continuum 
!     is interpolated (in temperature) separately.
  
  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT( IPLON,LAY)-1))*NSPA(11) + 1
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(11) + 1
      INDS(IPLON) = INDSELF(IPLON,LAY)
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG11
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 8
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG11
#endif
        TAU (IPLON,NGS10+IG,LAY) = COLH2O(IPLON,LAY) *        &
             (FAC00(IPLON,LAY) * ABSA(IND0(IPLON)  ,IG) +     &
              FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+1,IG) +     &
              FAC01(IPLON,LAY) * ABSA(IND1(IPLON)  ,IG) +     &
              FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+1,IG) +     &
             SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) +  &
             SELFFRAC(IPLON,LAY) *                            &
             (SELFREF(INDS(IPLON)+1,IG) - SELFREF(INDS(IPLON),IG)))) &
             + TAUAERL(IPLON,LAY,11)
        PFRAC(IPLON,NGS10+IG,LAY) = FRACREFA(IG)
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(11) + 1
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(11) + 1
   ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG11
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 8
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG11
#endif
        TAU (IPLON,NGS10+IG,LAY) = COLH2O(IPLON,LAY) *     &
             (FAC00(IPLON,LAY) * ABSB(IND0(IPLON)  ,IG) +  &
              FAC10(IPLON,LAY) * ABSB(IND0(IPLON)+1,IG) +  &
              FAC01(IPLON,LAY) * ABSB(IND1(IPLON)  ,IG) +  &
              FAC11(IPLON,LAY) * ABSB(IND1(IPLON)+1,IG))   &
             + TAUAERL(IPLON,LAY,11)
        PFRAC(IPLON,NGS10+IG,LAY) = FRACREFB(IG)
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL11
!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL12 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCO2,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND  ,JPGPT  ,NG12 ,NGS11
USE MO_RRTWN   , ONLY : NSPA
USE MO_RRTA12  , ONLY : ABSA   ,FRACREFA,SELFREF,STRRAT

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLCO2(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL (dp) :: ZFS(KBDIM), SPECCOMB(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL (dp) :: FS, SPECMULT, SPECPARM


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT*COLCO2(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 8._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(12) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(12) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG12
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!DIR$ UNROLL 8
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG12
#endif
        TAU (IPLON,NGS11+IG,LAY) = SPECCOMB(IPLON) *                    &
             ((1. - FS)*(FAC00(IPLON,LAY) * ABSA(IND0(IPLON)  ,IG) +    &
                         FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+9,IG) +    &
                         FAC01(IPLON,LAY) * ABSA(IND1(IPLON)  ,IG) +    &
                         FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+9,IG)) +   &
               FS *(     FAC01(IPLON,LAY) * ABSA(IND1(IPLON)+ 1,IG) +   &
                         FAC00(IPLON,LAY) * ABSA(IND0(IPLON)+ 1,IG) +   &
                         FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+10,IG) +   &
                         FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+10,IG))) + &
                        COLH2O(IPLON,LAY) *                             &
                       SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) +  &
                      SELFFRAC(IPLON,LAY) *                             &
                     (SELFREF(INDS(IPLON)+1,IG) - SELFREF(INDS(IPLON),IG)))  &
                    + TAUAERL(IPLON,LAY,12)
        PFRAC(IPLON,NGS11+IG,LAY) = FRACREFA(IG,JS) + FS * &
             (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
#ifndef UNROLLNGX
    DO IG = 1, NG12
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 8
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG12
#endif
        TAU  (IPLON,NGS11+IG,LAY) = TAUAERL(IPLON,LAY,12)
        PFRAC(IPLON,NGS11+IG,LAY) = 0.0_dp
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL12
!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL13 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLN2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 13:  2080-2250 cm-1 (low - H2O,N2O; high - nothing)

! Modifications
!
!     D Salmond 2000-05-15 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,NG13  ,NGS12
USE MO_RRTWN   , ONLY : NSPA
USE MO_RRTA13  , ONLY : ABSA   ,FRACREFA,SELFREF,STRRAT

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLN2O(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: FS, SPECMULT, SPECPARM


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately. 
 
  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT*COLN2O(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 8._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(13) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(13) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG13
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!DIR$ UNROLL 4
!CDIR UNROLL=4
!OCL  UNROLL(4)
      DO IG = 1, NG13
#endif
        TAU (IPLON,NGS12+IG,LAY) = SPECCOMB(IPLON) *                     &
             ((1. - FS)*(FAC00(IPLON,LAY) * ABSA(IND0(IPLON),IG) +       &
                         FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+9,IG) +     &
                         FAC01(IPLON,LAY) * ABSA(IND1(IPLON),IG) +       &
                         FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+9,IG)) +    &
                 FS* (   FAC01(IPLON,LAY) * ABSA(IND1(IPLON)+1,IG) +     &
                         FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+10,IG) +    &
                         FAC00(IPLON,LAY) * ABSA(IND0(IPLON)+1,IG) +     &
                         FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+10,IG)) ) + &
                        COLH2O(IPLON,LAY) *                              &
                       SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) +   &
                      SELFFRAC(IPLON,LAY) *                              &
                     (SELFREF(INDS(IPLON)+1,IG) - SELFREF(INDS(IPLON),IG))) &
                    + TAUAERL(IPLON,LAY,13)
        PFRAC(IPLON,NGS12+IG,LAY) = FRACREFA(IG,JS) + FS * &
             (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
#ifndef UNROLLNGX
    DO IG = 1, NG13
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 4
!CDIR UNROLL=4
!OCL  UNROLL(4)
      DO IG = 1, NG13
#endif
        TAU  (IPLON,NGS12+IG,LAY) = TAUAERL(IPLON,LAY,13)
        PFRAC(IPLON,NGS12+IG,LAY) = 0.0_dp
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL13
!******************************************************************************
SUBROUTINE RRTM_TAUMOL14 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
  &COLCO2,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)

! Modifications
!
!     D Salmond 1999-07-14 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,NGS13, NG14
USE MO_RRTWN   , ONLY : NSPA   ,NSPB
USE MO_RRTA14  , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
               &SELFREF

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PROFDATA             
REAL(DP):: COLCO2(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, LAY, IPLON
INTEGER :: IXP, IXC0


!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  Below LAYTROP, the water vapor self-continuum 
!     is interpolated (in temperature) separately.  

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(14) + 1
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(14) + 1
      INDS(IPLON) = INDSELF(IPLON,LAY)
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG14
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 2
!CDIR UNROLL=2
!OCL  UNROLL(2)
      DO IG = 1, NG14
#endif
        TAU (IPLON,NGS13+IG,LAY) = COLCO2(IPLON,LAY) *         &
             (FAC00(IPLON,LAY) * ABSA(IND0(IPLON)  ,IG) +      &
              FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+1,IG) +      &
              FAC01(IPLON,LAY) * ABSA(IND1(IPLON)  ,IG) +      &
              FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+1,IG) +      &
            SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) +    &
           SELFFRAC(IPLON,LAY) *                               &
           (SELFREF(INDS(IPLON)+1,IG) - SELFREF(INDS(IPLON),IG)))) &
          + TAUAERL(IPLON,LAY,14)
        PFRAC(IPLON,NGS13+IG,LAY) = FRACREFA(IG)
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(14) + 1
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(14) + 1
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG14
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 2
!CDIR UNROLL=2
!OCL  UNROLL(2)
      DO IG = 1, NG14
#endif
        TAU (IPLON,NGS13+IG,LAY) = COLCO2(IPLON,LAY) *     &
             (FAC00(IPLON,LAY) * ABSB(IND0(IPLON)  ,IG) +  &
              FAC10(IPLON,LAY) * ABSB(IND0(IPLON)+1,IG) +  &
              FAC01(IPLON,LAY) * ABSB(IND1(IPLON)  ,IG) +  &
              FAC11(IPLON,LAY) * ABSB(IND1(IPLON)+1,IG))   &
          + TAUAERL(IPLON,LAY,14)
        PFRAC(IPLON,NGS13+IG,LAY) = FRACREFB(IG)
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL14
!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL15 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCO2,COLN2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; high - nothing)

! Modifications
!
!     D Salmond 1999-07-14 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,NGS14, NG15
USE MO_RRTWN   , ONLY : NSPA
USE MO_RRTA15  , ONLY : ABSA   ,FRACREFA,SELFREF,STRRAT


IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLCO2(KBDIM,KLEV)
REAL(DP):: COLN2O(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: FS, SPECMULT, SPECPARM


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately. 
 
  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLN2O(IPLON,LAY) + STRRAT*COLCO2(IPLON,LAY)
      SPECPARM = COLN2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(SPECPARM,ONEMINUS)
      SPECMULT = 8._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(15) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(15) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG15
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!DIR$ UNROLL 2
!CDIR UNROLL=2
!OCL  UNROLL(2)
      DO IG = 1, NG15
#endif
        TAU (IPLON,NGS14+IG,LAY) = SPECCOMB(IPLON) *                    &
             ((1. - FS)*(FAC00(IPLON,LAY) * ABSA(IND0(IPLON),IG) +      &
                         FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+9,IG) +    &
                         FAC01(IPLON,LAY) * ABSA(IND1(IPLON),IG) +      &
                         FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+9,IG)) +   &
                  FS *  (FAC01(IPLON,LAY) * ABSA(IND1(IPLON)+1,IG) +    &
                         FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+10,IG) +   &
                         FAC00(IPLON,LAY) * ABSA(IND0(IPLON)+1,IG) +    &
                         FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+10,IG))) + &
                        COLH2O(IPLON,LAY) *                             &
                       SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) +  &
                      SELFFRAC(IPLON,LAY) *                             &
                     (SELFREF(INDS(IPLON)+1,IG) - SELFREF(INDS(IPLON),IG)))    &
                    + TAUAERL(IPLON,LAY,15)
        PFRAC(IPLON,NGS14+IG,LAY) = FRACREFA(IG,JS) + FS * &
             (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
#ifndef UNROLLNGX
    DO IG = 1, NG15
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 2
!CDIR UNROLL=2
!OCL  UNROLL(2)
      DO IG = 1, NG15
#endif
        TAU  (IPLON,NGS14+IG,LAY) = TAUAERL(IPLON,LAY,15)
        PFRAC(IPLON,NGS14+IG,LAY) = 0.0_dp
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL15
!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL16 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCH4,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)

! Modifications
!
!     D Salmond 1999-07-14 speed-up


USE MO_KIND    , ONLY : DP
USE MO_PARRRTM , ONLY : JPBAND ,JPGPT  ,NGS15
USE MO_RRTWN   , ONLY : NSPA
USE MO_RRTA16  , ONLY : ABSA   ,FRACREFA,SELFREF,STRRAT


IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLCH4(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, IPLON
INTEGER :: IXP, IXC0
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM)

!     LOCAL REAL SCALARS
REAL(DP):: FS, SPECMULT, SPECPARM


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately. 
 
  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT*COLCH4(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(SPECPARM,ONEMINUS)
      SPECMULT = 8._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+( JT(IPLON,LAY)-1))*NSPA(16) + JS
      IND1(IPLON) = (JP(IPLON,LAY)*5    +(JT1(IPLON,LAY)-1))*NSPA(16) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
!    DO IG = 1, NG16
    DO IG = 1, 2
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!DIR$ UNROLL 2
!CDIR UNROLL=2
!OCL  UNROLL(2)
!      DO IG = 1, NG16
      DO IG = 1, 2
#endif
        TAU (IPLON,NGS15+IG,LAY) = SPECCOMB(IPLON) *                  &
           ((1. - FS)*(FAC00(IPLON,LAY) * ABSA(IND0(IPLON),IG) +      &
                       FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+9,IG) +    &
                       FAC01(IPLON,LAY) * ABSA(IND1(IPLON),IG) +      &
                       FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+9,IG)) +   &
               FS * (  FAC01(IPLON,LAY) * ABSA(IND1(IPLON)+1,IG) +    &
                       FAC10(IPLON,LAY) * ABSA(IND0(IPLON)+10,IG) +   &
                       FAC00(IPLON,LAY) * ABSA(IND0(IPLON)+1,IG) +    &
                       FAC11(IPLON,LAY) * ABSA(IND1(IPLON)+10,IG))) + &
                      COLH2O(IPLON,LAY) *                             &
                     SELFFAC(IPLON,LAY) * (SELFREF(INDS(IPLON),IG) +  &
                    SELFFRAC(IPLON,LAY) *                             &
                    (SELFREF(INDS(IPLON)+1,IG) - SELFREF(INDS(IPLON),IG)))   &
                  + TAUAERL(IPLON,LAY,16)
        PFRAC(IPLON,NGS15+IG,LAY) = FRACREFA(IG,JS) + FS * &
             (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
#ifndef UNROLLNGX
!    DO IG = 1, NG16
    DO IG = 1, 2
#endif
!DIR$ CONCURRENT
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!DIR$ UNROLL 2
!CDIR UNROLL=2
!OCL  UNROLL(2)
!      DO IG = 1, NG16
      DO IG = 1, 2
#endif
        TAU  (IPLON,NGS15+IG,LAY) = TAUAERL(IPLON,LAY,16)
        PFRAC(IPLON,NGS15+IG,LAY) = 0.0_dp
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE RRTM_TAUMOL16
