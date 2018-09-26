SUBROUTINE RRTM_SETCOEF_140GP (KPROMA,KBDIM,KLEV,COLDRY,WKL &
 &, FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1 &
 &, COLH2O,COLCO2,COLO3,COLN2O,COLCH4,CO2MULT &
 &, LAYTROP,LAYSWTCH,LAYLOW,PAVEL,TAVEL,SELFFAC,SELFFRAC,INDSELF)

!     Reformatted for F90 by JJMorcrette, ECMWF, 980714
!     Include longitude loop, jan2001, Marco A. Giorgetta, MPI

!     Purpose:  For a given atmosphere, calculate the indices and
!     fractions related to the pressure and temperature interpolations.
!     Also calculate the values of the integrated Planck functions 
!     for each band at the level and layer temperatures.

USE MO_KIND   , ONLY : DP
USE MO_PARRRTM, ONLY : JPINPX
USE MO_RRTRF  , ONLY : PREFLOG   ,TREF

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV

REAL(DP):: COLDRY(KBDIM,KLEV)
REAL(DP):: WKL(KBDIM,JPINPX,KLEV)

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
REAL(DP):: COLCO2(KBDIM,KLEV)
REAL(DP):: COLO3 (KBDIM,KLEV)
REAL(DP):: COLN2O(KBDIM,KLEV)
REAL(DP):: COLCH4(KBDIM,KLEV)
!!$REAL(DP):: COLO2 (KBDIM,KLEV)
REAL(DP):: CO2MULT(KBDIM,KLEV)
INTEGER :: LAYTROP(KBDIM)
INTEGER :: LAYSWTCH(KBDIM)
INTEGER :: LAYLOW(KBDIM)

!- from PROFILE             
REAL(DP):: PAVEL(KBDIM,KLEV)
REAL(DP):: TAVEL(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)


  ! The following source code follows that provided by Uwe Schulzweida,
  ! as given in rrtm_setcoef_140gp.f90@@/main/clean/2,
  ! except that COLO2 is not computed.
  ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  !     LOCAL INTEGER SCALARS
  INTEGER :: JP1, LAY, IPLON

  !     LOCAL REAL SCALARS
  REAL(dp) :: CO2REG, COMPFP, FACTOR, FP, FT, FT1, PLOG, SCALEFAC, STPFAC, WATER


  STPFAC = 296._dp/1013._dp

  LAYTROP(:)  = 0
  LAYSWTCH(:) = 0
  LAYLOW(:)   = 0


  DO LAY = 1, KLEV
    DO IPLON = 1, KPROMA
      !        Find the two reference pressures on either side of the
      !        layer pressure.  Store them in JP and JP1.  Store in FP the
      !        fraction of the difference (in ln(pressure)) between these
      !        two values that the layer pressure lies.
      PLOG = LOG(PAVEL(IPLON,LAY))
      JP(IPLON,LAY) = INT(36._dp - 5*(PLOG+0.04_dp))
      IF (JP(IPLON,LAY)  <  1) THEN
        JP(IPLON,LAY) = 1
      ELSE IF (JP(IPLON,LAY)  >  58) THEN
        JP(IPLON,LAY) = 58
      ENDIF
      JP1 = JP(IPLON,LAY) + 1
      FP = 5._dp * (PREFLOG(JP(IPLON,LAY)) - PLOG)

      !        Determine, for each reference pressure (JP and JP1), which
      !        reference temperature (these are different for each  
      !        reference pressure) is nearest the layer temperature but does
      !        not exceed it.  Store these indices in JT and JT1, resp.
      !        Store in FT (resp. FT1) the fraction of the way between JT
      !        (JT1) and the next highest reference temperature that the 
      !        layer temperature falls.
      JT(IPLON,LAY) = INT(3._dp + (TAVEL(IPLON,LAY)-TREF(JP(IPLON,LAY)))/15._dp)
      IF (JT(IPLON,LAY)  <  1) THEN
        JT(IPLON,LAY) = 1
      ELSE IF (JT(IPLON,LAY)  >  4) THEN
        JT(IPLON,LAY) = 4
      ENDIF
      FT = ((TAVEL(IPLON,LAY)-TREF(JP(IPLON,LAY)))/15._dp) - REAL(JT(IPLON,LAY)-3)
      JT1(IPLON,LAY) = INT(3._dp + (TAVEL(IPLON,LAY)-TREF(JP1))/15._dp)
      IF (JT1(IPLON,LAY)  <  1) THEN
        JT1(IPLON,LAY) = 1
      ELSE IF (JT1(IPLON,LAY)  >  4) THEN
        JT1(IPLON,LAY) = 4
      ENDIF
      FT1 = ((TAVEL(IPLON,LAY)-TREF(JP1))/15._dp) - REAL(JT1(IPLON,LAY)-3)

      WATER = WKL(IPLON,1,LAY)/COLDRY(IPLON,LAY)
      SCALEFAC = PAVEL(IPLON,LAY) * STPFAC / TAVEL(IPLON,LAY)

      !        If the pressure is less than ~100mb, perform a different
      !        set of species interpolations.
      !         IF (PLOG .LE. 4.56) GO TO 5300
      !--------------------------------------         
      IF (PLOG  >  4.56_dp) THEN
        LAYTROP(IPLON) =  LAYTROP(IPLON) + 1
        !        For one band, the "switch" occurs at ~300 mb. 
        IF (PLOG  >=  5.76_dp) LAYSWTCH(IPLON) = LAYSWTCH(IPLON) + 1
        IF (PLOG  >=  6.62_dp) LAYLOW(IPLON) = LAYLOW(IPLON) + 1

        FORFAC(IPLON,LAY) = SCALEFAC / (1.0_dp+WATER)

        !        Set up factors needed to separately include the water vapor
        !        self-continuum in the calculation of absorption coefficient.
        !C           SELFFAC(IPLON,LAY) = WATER * SCALEFAC / (1.+WATER)
        SELFFAC(IPLON,LAY) = WATER * FORFAC(IPLON,LAY)
        FACTOR = (TAVEL(IPLON,LAY)-188.0_dp)/7.2_dp
        INDSELF(IPLON,LAY) = MIN(9, MAX(1, INT(FACTOR)-7))
        SELFFRAC(IPLON,LAY) = FACTOR - REAL(INDSELF(IPLON,LAY) + 7,dp)

        !        Calculate needed column amounts.
        COLH2O(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,1,LAY)
        COLCO2(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,2,LAY)
        COLO3(IPLON,LAY)  = 1.E-20_dp * WKL(IPLON,3,LAY)
        COLN2O(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,4,LAY)
        COLCH4(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,6,LAY)
!!$        COLO2(IPLON,LAY)  = 1.E-20_dp * WKL(IPLON,7,LAY)
        IF (COLCO2(IPLON,LAY)  ==  0.0_dp) COLCO2(IPLON,LAY) = 1.E-32_dp * COLDRY(IPLON,LAY)
        IF (COLN2O(IPLON,LAY)  ==  0.0_dp) COLN2O(IPLON,LAY) = 1.E-32_dp * COLDRY(IPLON,LAY)
        IF (COLCH4(IPLON,LAY)  ==  0.0_dp) COLCH4(IPLON,LAY) = 1.E-32_dp * COLDRY(IPLON,LAY)
        !        Using E = 1334.2 cm-1.
        CO2REG = 3.55E-24_dp * COLDRY(IPLON,LAY)
        CO2MULT(IPLON,LAY)= (COLCO2(IPLON,LAY) - CO2REG) *&
             272.63_dp*EXP(-1919.4_dp/TAVEL(IPLON,LAY))/(8.7604E-4_dp*TAVEL(IPLON,LAY))
        !         GO TO 5400
        !------------------
      ELSE
        !        Above LAYTROP.
        ! 5300    CONTINUE

        !        Calculate needed column amounts.
        FORFAC(IPLON,LAY) = SCALEFAC / (1.0_dp+WATER)

        COLH2O(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,1,LAY)
        COLCO2(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,2,LAY)
        COLO3(IPLON,LAY)  = 1.E-20_dp * WKL(IPLON,3,LAY)
        COLN2O(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,4,LAY)
        COLCH4(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,6,LAY)
!!$        COLO2(IPLON,LAY)  = 1.E-20_dp * WKL(IPLON,7,LAY)
        IF (COLCO2(IPLON,LAY)  ==  0.0_dp) COLCO2(IPLON,LAY) = 1.E-32_dp * COLDRY(IPLON,LAY)
        IF (COLN2O(IPLON,LAY)  ==  0.0_dp) COLN2O(IPLON,LAY) = 1.E-32_dp * COLDRY(IPLON,LAY)
        IF (COLCH4(IPLON,LAY)  ==  0.0_dp) COLCH4(IPLON,LAY) = 1.E-32_dp * COLDRY(IPLON,LAY)
        CO2REG = 3.55E-24_dp * COLDRY(IPLON,LAY)
        CO2MULT(IPLON,LAY)= (COLCO2(IPLON,LAY) - CO2REG) *&
             272.63_dp*EXP(-1919.4_dp/TAVEL(IPLON,LAY))/(8.7604E-4_dp*TAVEL(IPLON,LAY))
        !----------------     
      ENDIF
      ! 5400    CONTINUE

      !        We have now isolated the layer ln pressure and temperature,
      !        between two reference pressures and two reference temperatures 
      !        (for each reference pressure).  We multiply the pressure 
      !        fraction FP with the appropriate temperature fractions to get 
      !        the factors that will be needed for the interpolation that yields
      !        the optical depths (performed in routines TAUGBn for band n).

      COMPFP = 1.0_dp - FP
      FAC10(IPLON,LAY) = COMPFP * FT
      FAC00(IPLON,LAY) = COMPFP * (1.0_dp - FT)
      FAC11(IPLON,LAY) = FP * FT1
      FAC01(IPLON,LAY) = FP * (1.0_dp - FT1)

    ENDDO
  ENDDO

  DO IPLON = 1, KPROMA
    ! MT 981104 
    !-- Set LAYLOW for profiles with surface pressure less than 750 hPa. 
    IF (LAYLOW(IPLON) == 0) LAYLOW(IPLON)=1
  ENDDO

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ! End of source code from Uwe Schulzweidas version
  ! rrtm_setcoef_140gp.f90@@/main/clean/2

RETURN
END SUBROUTINE RRTM_SETCOEF_140GP
