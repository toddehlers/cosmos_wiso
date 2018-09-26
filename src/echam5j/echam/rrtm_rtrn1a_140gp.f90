SUBROUTINE RRTM_RTRN1A_140GP (KPROMA,KBDIM,KLEV,ISTART,IEND,ICLDLYR,CLDFRAC,TAUCLD,ABSS1 &
  &, OD,TAUSF1,TOTDFLUC,TOTDFLUX,TOTUFLUC,TOTUFLUX &
  &, TAVEL,TZ,TBOUND,PFRAC,SEMISS,SEMISLW)

!     Reformatted for F90 by JJMorcrette, ECMWF, 980714
!     Speed-up by D.Salmond, ECMWF, 9907
!     Bug-fix by M.J. Iacono, AER, Inc., 9911
!     Bug-fix by JJMorcrette, ECMWF, 991209 (RAT1, RAT2 initialization)
!     Speed-up by D. Salmond, ECMWF, 9912
!     Bug-fix by JJMorcrette, ECMWF, 0005 (extrapolation T<160K)
!     Speed-up by D. Salmond, ECMWF, 000515
!     Speed-up by Uwe Schulzweida, MPIFMET, 2001/01
!     Change data i/o from modules to argument lists,
!        and include new speed-up of Uwe Schulzweida
!                     March 2001, Marco A. Giorgetta
!     U. Schulzweida, MPI, May 2002, blocking (nproma)

!-* This program calculates the upward fluxes, downward fluxes,
!   and heating rates for an arbitrary atmosphere.  The input to
!   this program is the atmospheric profile and all Planck function
!   information.  First-order "numerical" quadrature is used for the 
!   angle integration, i.e. only one exponential is computed per layer
!   per g-value per band.  Cloud overlap is treated with a generalized
!   maximum/random method in which adjacent cloud layers are treated
!   with maximum overlap, and non-adjacent cloud groups are treated
!   with random overlap.  For adjacent cloud layers, cloud information
!   is carried from the previous two layers.


  USE MO_KIND   , ONLY : DP          ! Kind parameter for REAL precision
  USE MO_PARRRTM, ONLY : JPBAND , &  ! Number of longwave spectral bands
       &                 JPGPT       ! Total number of g-point subintervals
  USE MO_RRTAB  , ONLY : BPADE       ! Pade constant
  USE MO_RRTWN  , ONLY : TOTPLNK, &  ! Planck function for spectral bands
       &                 DELWAVE     ! Band width
  USE MO_RRTFTR , ONLY : NGB         ! Band index of g-points


  IMPLICIT NONE

  ! Argument list
  ! =============

  ! Input
  ! -----
  INTEGER :: KPROMA                    ! Number of columns
  INTEGER :: KBDIM                     ! first dimension of 2-d arrays
  INTEGER :: KLEV                      ! Number of layers
  INTEGER :: ISTART                    ! Index of 1st band for Planck emission
  INTEGER :: IEND                      ! Index of last band for Planck emission
  INTEGER :: ICLDLYR(KBDIM,KLEV)       ! Cloud indicator
  REAL(DP):: CLDFRAC(KBDIM,KLEV)       ! Cloud fraction
  REAL(DP):: TAUCLD(KBDIM,KLEV,JPBAND) ! Spectral cloud optical thickness
  REAL(DP):: ABSS1 (KBDIM,JPGPT*KLEV)  ! 
  REAL(DP):: OD    (KBDIM,JPGPT,KLEV)  ! Clear-sky optical thickness
  REAL(DP):: TAUSF1(KBDIM,JPGPT*KLEV)  ! 
  REAL(DP):: TAVEL(KBDIM,KLEV)         ! Layer temperature
  REAL(DP):: TZ(KBDIM,0:KLEV)          ! Level temperature
  REAL(DP):: TBOUND(KBDIM)             ! Surface temperature
  REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)   ! Planck function fractions
  REAL(DP):: SEMISS(KBDIM,JPBAND)      ! Surface spectral emissivity
  REAL(DP):: SEMISLW(KBDIM)            ! Surface emissivity

  ! Output
  ! ------
  REAL(DP):: TOTDFLUC(KBDIM,0:KLEV)    ! Clear-sky downward longwave flux
  REAL(DP):: TOTDFLUX(KBDIM,0:KLEV)    ! Downward longwave flux
  REAL(DP):: TOTUFLUC(KBDIM,0:KLEV)    ! Clear-sky upward longwave flux
  REAL(DP):: TOTUFLUX(KBDIM,0:KLEV)    ! Upward longwave flux


  !--------------------------------------------------------------------------
  !
  ! Maximum/Random cloud overlap variables
  ! for upward radiaitve transfer
  !  FACCLR2  fraction of clear radiance from previous layer that needs to 
  !           be switched to cloudy stream
  !  FACCLR1  fraction of the radiance that had been switched in the previous
  !           layer from cloudy to clear that needs to be switched back to
  !           cloudy in the current layer
  !  FACCLD2  fraction of cloudy radiance from previous layer that needs to 
  !           be switched to clear stream
  !  FACCLD1  fraction of the radiance that had been switched in the previous
  !           layer from clear to cloudy that needs to be switched back to
  !           clear in the current layer
  ! for downward radiaitve transfer
  !  FACCLR2D fraction of clear radiance from previous layer that needs to 
  !           be switched to cloudy stream
  !  FACCLR1D fraction of the radiance that had been switched in the previous
  !           layer from cloudy to clear that needs to be switched back to
  !           cloudy in the current layer
  !  FACCLD2D fraction of cloudy radiance from previous layer that needs to 
  !           be switched to clear stream
  !  FACCLD1D fraction of the radiance that had been switched in the previous
  !           layer from clear to cloudy that needs to be switched back to
  !           clear in the current layer
  !
  !--------------------------------------------------------------------------
  

  ! The following source code follows that provided by Uwe Schulzweida,
  ! as given in rrtm_rtrn1a_140gp.f90@@/main/clean/2
  ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  INTEGER :: INDLAY(KPROMA,KLEV),INDLEV(KPROMA,0:KLEV)

  REAL(dp) :: BBU1(KPROMA,JPGPT*KLEV),BBUTOT1(KPROMA,JPGPT*KLEV)
  REAL(dp) :: TLAYFRAC(KPROMA,KLEV),TLEVFRAC(KPROMA,0:KLEV)
  REAL(dp) :: BGLEV(KPROMA,JPGPT)
  REAL(dp) :: PLVL(KPROMA,JPBAND,0:KLEV),PLAY(KPROMA,JPBAND,0:KLEV),WTNUM(3)
  REAL(dp) :: SEMIS(KPROMA,JPGPT),RADUEMIT(KPROMA,JPGPT)

  REAL(dp) :: RADCLRU1(KPROMA,JPGPT) ,RADCLRD1(KPROMA,JPGPT)
  REAL(dp) :: RADLU1(KPROMA,JPGPT)   ,RADLD1(KPROMA,JPGPT)
  REAL(dp) :: TRNCLD(KPROMA,JPBAND,KLEV)
  REAL(dp) :: ABSCLDNW
  REAL(dp) :: ATOT1(KPROMA,JPGPT*KLEV)

  REAL(dp) :: SURFEMIS(KPROMA,JPBAND),PLNKEMIT(KPROMA,JPBAND)

  !******
  REAL(dp) :: clrradu(KPROMA,jpgpt),cldradu(KPROMA,jpgpt),oldcld
  REAL(dp) :: oldclr,rad(KPROMA,jpgpt)
  REAL(dp) :: faccld1(KPROMA,klev+1),faccld2(KPROMA,klev+1)
  REAL(dp) :: facclr1(KPROMA,klev+1),facclr2(KPROMA,klev+1)
  REAL(dp) :: faccmb1(KPROMA,klev+1),faccmb2(KPROMA,klev+1)
  REAL(dp) :: faccld1d(KPROMA,0:klev),faccld2d(KPROMA,0:klev)
  REAL(dp) :: facclr1d(KPROMA,0:klev),facclr2d(KPROMA,0:klev)
  REAL(dp) :: faccmb1d(KPROMA,0:klev),faccmb2d(KPROMA,0:klev)
  REAL(dp) :: clrradd(KPROMA,jpgpt),cldradd(KPROMA,jpgpt)
  INTEGER :: istcld(KPROMA,klev+1),istcldd(KPROMA,0:klev)

  !     LOCAL INTEGER SCALARS
  INTEGER :: ICCLD, ICCLDL(KPROMA), IX      
  INTEGER :: IPLON
  INTEGER :: IBAND, ICLDDN(KPROMA), IENT, INDBOUND(KPROMA), INDEX, IPR, LAY, LEV, NBI

  !     LOCAL REAL SCALARS
  REAL(dp) :: BBD(KPROMA,jpgpt),DELBGDN(KPROMA,jpgpt)
  REAL(dp) :: DELBGUP(KPROMA,jpgpt),BGLAY(KPROMA,jpgpt)
  REAL(dp) :: BBDTOT, CLDSRC, DBDTLAY, DBDTLEV,                    &
              DRAD1(KPROMA), DRADCL1(KPROMA), FACTOT1,           &
              FMAX, FMIN, GASSRC, ODSM, PLANKBND, RADCLD,                      &
              RADD, RADMOD, RAT1(KPROMA), RAT2(KPROMA), SUMPL(KPROMA),               &
              SUMPLEM(KPROMA), TBNDFRAC(KPROMA), TRNS, TTOT, URAD1(KPROMA), URADCL1(KPROMA)
  REAL(dp) :: CLDRADD_L,CLRRADD_L,RAD_L
  !******


  WTNUM(1) = 0.5_dp
  WTNUM(2) = 0.0_dp
  WTNUM(3) = 0.0_dp

  DO IPLON = 1, KPROMA
    !-start JJM_000511
    IF (TBOUND(IPLON) < 339._dp .AND. TBOUND(IPLON) >= 160._dp ) THEN
      INDBOUND(IPLON) = TBOUND(IPLON) - 159._dp
      TBNDFRAC(IPLON) = TBOUND(IPLON) - INT(TBOUND(IPLON))
    ELSE IF (TBOUND(IPLON) >= 339._dp ) THEN
      INDBOUND(IPLON) = 180
      TBNDFRAC(IPLON) = TBOUND(IPLON) - 339._dp
    ELSE IF (TBOUND(IPLON) < 160._dp ) THEN
      INDBOUND(IPLON) = 1
      TBNDFRAC(IPLON) = TBOUND(IPLON) - 160._dp
    ENDIF
    !-end JJM_000511
  END DO

  TOTUFLUC(1:KPROMA,:) = 0.0_dp
  TOTDFLUC(1:KPROMA,:) = 0.0_dp
  TOTUFLUX(1:KPROMA,:) = 0.0_dp
  TOTDFLUX(1:KPROMA,:) = 0.0_dp

  TRNCLD(1:KPROMA,:,:) = 0.0_dp

  DO LAY = 0, KLEV
    DO IPLON = 1, KPROMA
      !-start JJM_000511
      IF (TZ(IPLON,LAY) < 339._dp .AND. TZ(IPLON,LAY) >= 160._dp ) THEN
        INDLEV(IPLON,LAY) = TZ(IPLON,LAY) - 159._dp
        TLEVFRAC(IPLON,LAY) = TZ(IPLON,LAY) - INT(TZ(IPLON,LAY))
      ELSE IF (TZ(IPLON,LAY) >= 339._dp ) THEN
        INDLEV(IPLON,LAY) = 180
        TLEVFRAC(IPLON,LAY) = TZ(IPLON,LAY) - 339._dp
      ELSE IF (TZ(IPLON,LAY) < 160._dp ) THEN
        INDLEV(IPLON,LAY) = 1
        TLEVFRAC(IPLON,LAY) = TZ(IPLON,LAY) - 160._dp
      ENDIF
      !-end JJM_000511
    ENDDO
  ENDDO

  !- jjm_991209
  FACCLD1(1:KPROMA,:) = 0.0_dp
  FACCLD2(1:KPROMA,:) = 0.0_dp
  FACCLR1(1:KPROMA,:) = 0.0_dp
  FACCLR2(1:KPROMA,:) = 0.0_dp
  FACCMB1(1:KPROMA,:) = 0.0_dp
  FACCMB2(1:KPROMA,:) = 0.0_dp
  FACCLD1D(1:KPROMA,:)  = 0.0_dp
  FACCLD2D(1:KPROMA,:)  = 0.0_dp
  FACCLR1D(1:KPROMA,:)  = 0.0_dp
  FACCLR2D(1:KPROMA,:)  = 0.0_dp
  FACCMB1D(1:KPROMA,:)  = 0.0_dp
  FACCMB2D(1:KPROMA,:)  = 0.0_dp

  RAT1(1:KPROMA)    = 0.0_dp
  RAT2(1:KPROMA)    = 0.0_dp

  !- jjm_991209

  SUMPL(1:KPROMA)   = 0.0_dp
  SUMPLEM(1:KPROMA) = 0.0_dp

  ISTCLD(1:KPROMA,1)     = 1
  ISTCLDD(1:KPROMA,KLEV) = 1

  DO LEV = 1, KLEV
    DO IPLON = 1, KPROMA
      !-- DS_000515
      !-start JJM_000511
      IF (TAVEL(IPLON,LEV) < 339._dp .AND. TAVEL(IPLON,LEV) >= 160._dp ) THEN
        INDLAY(IPLON,LEV) = TAVEL(IPLON,LEV) - 159._dp
        TLAYFRAC(IPLON,LEV) = TAVEL(IPLON,LEV) - INT(TAVEL(IPLON,LEV))
      ELSE IF (TAVEL(IPLON,LEV) >= 339._dp ) THEN
        INDLAY(IPLON,LEV) = 180
        TLAYFRAC(IPLON,LEV) = TAVEL(IPLON,LEV) - 339._dp
      ELSE IF (TAVEL(IPLON,LEV) < 160._dp ) THEN
        INDLAY(IPLON,LEV) = 1
        TLAYFRAC(IPLON,LEV) = TAVEL(IPLON,LEV) - 160._dp
      ENDIF
      !-end JJM_000511
      !-- DS_000515
    END DO
  END DO

  DO LEV = 1, KLEV
    DO IPLON = 1, KPROMA
      IF (ICLDLYR(IPLON,LEV) == 1) THEN

        !mji    
        ISTCLD(IPLON,LEV+1) = 0

        IF (LEV  ==  KLEV) THEN
          FACCLD1(IPLON,LEV+1) = 0.0_dp
          FACCLD2(IPLON,LEV+1) = 0.0_dp
          FACCLR1(IPLON,LEV+1) = 0.0_dp
          FACCLR2(IPLON,LEV+1) = 0.0_dp
          FACCMB1(IPLON,LEV+1) = 0.0_dp
          FACCMB2(IPLON,LEV+1) = 0.0_dp
          !mji      ISTCLD(IPLON,LEV+1) = 0

        ELSE IF (CLDFRAC(IPLON,LEV+1)  >=  CLDFRAC(IPLON,LEV)) THEN
          FACCLD1(IPLON,LEV+1) = 0.0_dp
          FACCLD2(IPLON,LEV+1) = 0.0_dp


          IF (ISTCLD(IPLON,LEV)  ==  1) THEN
            !mji        ISTCLD(IPLON,LEV+1) = 0
            FACCLR1(IPLON,LEV+1) = 0.0_dp
            !mji        
            FACCLR2(IPLON,LEV+1) = 0.0_dp
            IF (CLDFRAC(IPLON,LEV) < 1.0_dp) THEN
              FACCLR2(IPLON,LEV+1) = (CLDFRAC(IPLON,LEV+1)-CLDFRAC(IPLON,LEV)) &
                                   / (1.0_dp-CLDFRAC(IPLON,LEV))
            END IF
          ELSE
            FMAX = MAX(CLDFRAC(IPLON,LEV),CLDFRAC(IPLON,LEV-1))
            !mji
            IF (CLDFRAC(IPLON,LEV+1)  >  FMAX) THEN
              FACCLR1(IPLON,LEV+1) = RAT2(IPLON)
              FACCLR2(IPLON,LEV+1) = (CLDFRAC(IPLON,LEV+1)-FMAX)/(1.0_dp-FMAX)
              !mji          
            ELSE IF (CLDFRAC(IPLON,LEV+1) < FMAX) THEN
              FACCLR1(IPLON,LEV+1) = (CLDFRAC(IPLON,LEV+1)-CLDFRAC(IPLON,LEV)) &
                                   / (CLDFRAC(IPLON,LEV-1)-CLDFRAC(IPLON,LEV))
              FACCLR2(IPLON,LEV+1) = 0.0_dp
              !mji
            ELSE
              FACCLR1(IPLON,LEV+1) = RAT2(IPLON)
              FACCLR2(IPLON,LEV+1) = 0.0_dp
            ENDIF
          ENDIF
          IF (FACCLR1(IPLON,LEV+1) > 0.0_dp .OR. FACCLR2(IPLON,LEV+1) > 0.0_dp) THEN
            RAT1(IPLON) = 1.0_dp
            RAT2(IPLON) = 0.0_dp
          ENDIF
        ELSE
          FACCLR1(IPLON,LEV+1) = 0.0_dp
          FACCLR2(IPLON,LEV+1) = 0.0_dp
          IF (ISTCLD(IPLON,LEV)  ==  1) THEN
            !mji        ISTCLD(IPLON,LEV+1) = 0
            FACCLD1(IPLON,LEV+1) = 0.0_dp
            FACCLD2(IPLON,LEV+1) = (CLDFRAC(IPLON,LEV)-CLDFRAC(IPLON,LEV+1)) &
                                 /  CLDFRAC(IPLON,LEV)
          ELSE
            FMIN = MIN(CLDFRAC(IPLON,LEV),CLDFRAC(IPLON,LEV-1))
            IF (CLDFRAC(IPLON,LEV+1)  <=  FMIN) THEN
              FACCLD1(IPLON,LEV+1) = RAT1(IPLON)
              FACCLD2(IPLON,LEV+1) = (FMIN-CLDFRAC(IPLON,LEV+1))/FMIN
            ELSE
              FACCLD1(IPLON,LEV+1) = (CLDFRAC(IPLON,LEV)-CLDFRAC(IPLON,LEV+1)) &
                                   / (CLDFRAC(IPLON,LEV)-FMIN)
              FACCLD2(IPLON,LEV+1) = 0.0_dp
            ENDIF
          ENDIF
          IF (FACCLD1(IPLON,LEV+1) > 0.0_dp .OR. FACCLD2(IPLON,LEV+1) > 0.0_dp) THEN
            RAT1(IPLON) = 0.0_dp
            RAT2(IPLON) = 1.0_dp
          ENDIF
        ENDIF
        !fcc
        IF (LEV == 1) THEN
          FACCMB1(IPLON,LEV+1) = 0._dp
          FACCMB2(IPLON,LEV+1) = FACCLD1(IPLON,LEV+1) * FACCLR2(IPLON,LEV)
        ELSE
          FACCMB1(IPLON,LEV+1) = FACCLR1(IPLON,LEV+1) * FACCLD2(IPLON,LEV) * &
                                 CLDFRAC(IPLON,LEV-1)
          FACCMB2(IPLON,LEV+1) = FACCLD1(IPLON,LEV+1) * FACCLR2(IPLON,LEV) * &
                                 (1.0_dp - CLDFRAC(IPLON,LEV-1)) 
        ENDIF
        !end fcc
      ELSE
        !-- DS_000515
        ISTCLD(IPLON,LEV+1) = 1
      ENDIF
    ENDDO
  ENDDO

  !- jjm_991209
  RAT1(1:KPROMA)    = 0.0_dp
  RAT2(1:KPROMA)    = 0.0_dp
  !- jjm_991209

  DO LEV = KLEV, 1, -1
    DO IPLON = 1, KPROMA
      IF (ICLDLYR(IPLON,LEV) == 1) THEN
        !mji
        ISTCLDD(IPLON,LEV-1) = 0  
        IF (LEV  ==  1) THEN
          FACCLD1D(IPLON,LEV-1) = 0.0_dp
          FACCLD2D(IPLON,LEV-1) = 0.0_dp
          FACCLR1D(IPLON,LEV-1) = 0.0_dp
          FACCLR2D(IPLON,LEV-1) = 0.0_dp
          FACCMB1D(IPLON,LEV-1) = 0.0_dp
          FACCMB2D(IPLON,LEV-1) = 0.0_dp
          !mji      ISTCLDD(IPLON,LEV-1) = 0.0_dp
          ELSEIF (CLDFRAC(IPLON,LEV-1)  >=  CLDFRAC(IPLON,LEV)) THEN
          FACCLD1D(IPLON,LEV-1) = 0.0_dp
          FACCLD2D(IPLON,LEV-1) = 0.0_dp
          IF (ISTCLDD(IPLON,LEV)  ==  1) THEN
            !mji        ISTCLDD(IPLON,LEV-1) = 0
            FACCLR1D(IPLON,LEV-1) = 0.0_dp
            FACCLR2D(IPLON,LEV-1) = 0.0_dp
            IF (CLDFRAC(IPLON,LEV) < 1.0_dp) THEN
              FACCLR2D(IPLON,LEV-1) = (CLDFRAC(IPLON,LEV-1)-CLDFRAC(IPLON,LEV))&
                                    / (1.0_dp-CLDFRAC(IPLON,LEV))
            END IF
          ELSE
            FMAX = MAX(CLDFRAC(IPLON,LEV),CLDFRAC(IPLON,LEV+1))
            !mji
            IF (CLDFRAC(IPLON,LEV-1)  >  FMAX) THEN
              FACCLR1D(IPLON,LEV-1) = RAT2(IPLON)
              FACCLR2D(IPLON,LEV-1) = (CLDFRAC(IPLON,LEV-1)-FMAX)/(1.0_dp-FMAX)
              !mji
            ELSE IF (CLDFRAC(IPLON,LEV-1) < FMAX) THEN
              FACCLR1D(IPLON,LEV-1) = (CLDFRAC(IPLON,LEV-1)-CLDFRAC(IPLON,LEV))&
                                    / (CLDFRAC(IPLON,LEV+1)-CLDFRAC(IPLON,LEV))
              FACCLR2D(IPLON,LEV-1) = 0.0_dp
              !mji
            ELSE          
              FACCLR1D(IPLON,LEV-1) = RAT2(IPLON)
              FACCLR2D(IPLON,LEV-1) = 0.0_dp
            ENDIF
          ENDIF
          IF (FACCLR1D(IPLON,LEV-1) > 0.0_dp .OR. FACCLR2D(IPLON,LEV-1) > 0.0_dp)THEN
            RAT1(IPLON) = 1.0_dp
            RAT2(IPLON) = 0.0_dp
          ENDIF
        ELSE
          FACCLR1D(IPLON,LEV-1) = 0.0_dp
          FACCLR2D(IPLON,LEV-1) = 0.0_dp
          IF (ISTCLDD(IPLON,LEV)  ==  1) THEN
            !mji        ISTCLDD(IPLON,LEV-1) = 0
            FACCLD1D(IPLON,LEV-1) = 0.0_dp
            FACCLD2D(IPLON,LEV-1) = (CLDFRAC(IPLON,LEV)-CLDFRAC(IPLON,LEV-1)) &
                                  /  CLDFRAC(IPLON,LEV)
          ELSE
            FMIN = MIN(CLDFRAC(IPLON,LEV),CLDFRAC(IPLON,LEV+1))
            IF (CLDFRAC(IPLON,LEV-1)  <=  FMIN) THEN
              FACCLD1D(IPLON,LEV-1) = RAT1(IPLON)
              FACCLD2D(IPLON,LEV-1) = (FMIN-CLDFRAC(IPLON,LEV-1))/FMIN
            ELSE
              FACCLD1D(IPLON,LEV-1) = (CLDFRAC(IPLON,LEV)-CLDFRAC(IPLON,LEV-1))&
                                    / (CLDFRAC(IPLON,LEV)-FMIN)
              FACCLD2D(IPLON,LEV-1) = 0.0_dp
            ENDIF
          ENDIF
          IF (FACCLD1D(IPLON,LEV-1) > 0.0_dp .OR. FACCLD2D(IPLON,LEV-1) > 0.0_dp)THEN
            RAT1(IPLON) = 0.0_dp
            RAT2(IPLON) = 1.0_dp
          ENDIF
        ENDIF
        !mag
        IF (LEV == KLEV) THEN
          FACCMB1D(IPLON,LEV-1) = 0._dp
          FACCMB2D(IPLON,LEV-1) = FACCLD1D(IPLON,LEV-1) * FACCLR2D(IPLON,LEV)
        ELSE
          FACCMB1D(IPLON,LEV-1) = FACCLR1D(IPLON,LEV-1) * FACCLD2D(IPLON,LEV) &
                                * CLDFRAC(IPLON,LEV+1)
          FACCMB2D(IPLON,LEV-1) = FACCLD1D(IPLON,LEV-1) * FACCLR2D(IPLON,LEV) &
                                * (1.0_dp - CLDFRAC(IPLON,LEV+1))
        ENDIF
        !end mag
      ELSE
        ISTCLDD(IPLON,LEV-1) = 1
      ENDIF
    ENDDO
  ENDDO

  !- Loop over frequency bands.

  DO IBAND = ISTART, IEND
    DO IPLON = 1, KPROMA
      DBDTLEV  = TOTPLNK(INDBOUND(IPLON)+1,IBAND)-TOTPLNK(INDBOUND(IPLON),IBAND)
      PLANKBND = DELWAVE(IBAND) &
               * (TOTPLNK(INDBOUND(IPLON),IBAND) + TBNDFRAC(IPLON) * DBDTLEV)
      DBDTLEV  = TOTPLNK(INDLEV(IPLON,0)+1,IBAND)-TOTPLNK(INDLEV(IPLON,0),IBAND)
      PLVL(IPLON,IBAND,0) = DELWAVE(IBAND) &
           * (TOTPLNK(INDLEV(IPLON,0),IBAND) + TLEVFRAC(IPLON,0)*DBDTLEV)

      SURFEMIS(IPLON,IBAND) = SEMISS(IPLON,IBAND)
      PLNKEMIT(IPLON,IBAND) = SURFEMIS(IPLON,IBAND) * PLANKBND
      SUMPLEM(IPLON)  = SUMPLEM(IPLON) + PLNKEMIT(IPLON,IBAND)
      SUMPL(IPLON)    = SUMPL(IPLON)   + PLANKBND
      !--DS
    ENDDO
  ENDDO

  !---
  DO IBAND = ISTART, IEND
    DO LEV = 1, KLEV
      DO IPLON = 1, KPROMA
        !----              
        !- Calculate the integrated Planck functions for at the
        !  level and layer temperatures.
        !  Compute cloud transmittance for cloudy layers.
        DBDTLEV = TOTPLNK(INDLEV(IPLON,LEV)+1,IBAND) &
                - TOTPLNK(INDLEV(IPLON,LEV),IBAND)
        DBDTLAY = TOTPLNK(INDLAY(IPLON,LEV)+1,IBAND) &
                - TOTPLNK(INDLAY(IPLON,LEV),IBAND)
        PLAY(IPLON,IBAND,LEV) = DELWAVE(IBAND) &
             * (TOTPLNK(INDLAY(IPLON,LEV),IBAND)+TLAYFRAC(IPLON,LEV)*DBDTLAY)
        PLVL(IPLON,IBAND,LEV) = DELWAVE(IBAND) &
             * (TOTPLNK(INDLEV(IPLON,LEV),IBAND)+TLEVFRAC(IPLON,LEV)*DBDTLEV)
        IF (ICLDLYR(IPLON,LEV) > 0) THEN
          TRNCLD(IPLON,IBAND,LEV) = EXP(-TAUCLD(IPLON,LEV,IBAND))
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  SEMISLW(1:KPROMA) = SUMPLEM(1:KPROMA) / SUMPL(1:KPROMA)

  !- Initialize for radiative transfer.
  DO IPR = 1, JPGPT
    NBI = NGB(IPR)
    DO IPLON = 1, KPROMA
      RADCLRD1(IPLON,IPR) = 0.0_dp
      RADLD1(IPLON,IPR)   = 0.0_dp
      SEMIS(IPLON,IPR)    = SURFEMIS(IPLON,NBI)
      RADUEMIT(IPLON,IPR) = PFRAC(IPLON,IPR,1) * PLNKEMIT(IPLON,NBI)
      BGLEV(IPLON,IPR)    = PFRAC(IPLON,IPR,KLEV) * PLVL(IPLON,NBI,KLEV)
    ENDDO
  ENDDO

  !- Downward radiative transfer.
  !  *** DRAD1 holds summed radiance for total sky stream
  !  *** DRADCL1 holds summed radiance for clear sky stream

  ICLDDN(1:KPROMA) = 0

  DO LEV = KLEV, 1, -1
    DRAD1(1:KPROMA)   = 0.0_dp
    DRADCL1(1:KPROMA) = 0.0_dp

    ICCLD = 0
    DO IPLON = 1, KPROMA
      IF (ICLDLYR(IPLON,LEV) == 1) THEN
        ICCLD = ICCLD + 1
        ICCLDL(ICCLD) = IPLON
      ENDIF
    ENDDO

    IENT = JPGPT * (LEV-1)
    DO IPR = 1, JPGPT
      NBI = NGB(IPR)
      INDEX = IENT + IPR

      DO IPLON = 1, KPROMA

        BGLAY(IPLON,IPR) = PFRAC(IPLON,IPR,LEV) * PLAY(IPLON,NBI,LEV)
        !----            
        DELBGUP(IPLON,IPR)     = BGLEV(IPLON,IPR) - BGLAY(IPLON,IPR)
        BBU1(IPLON,INDEX) = BGLAY(IPLON,IPR) + TAUSF1(IPLON,INDEX) * DELBGUP(IPLON,IPR)
        !--DS            
        BGLEV(IPLON,IPR) = PFRAC(IPLON,IPR,LEV) * PLVL(IPLON,NBI,LEV-1)
        !----            
        DELBGDN(IPLON,IPR) = BGLEV(IPLON,IPR) - BGLAY(IPLON,IPR)
        BBD(IPLON,IPR) = BGLAY(IPLON,IPR) + TAUSF1(IPLON,INDEX) * DELBGDN(IPLON,IPR)
      ENDDO
    ENDDO

    DO IPR = 1, JPGPT
      NBI = NGB(IPR)
      INDEX = IENT + IPR

!CDIR NODEP
!OCL NOVREC
!DIR$ CONCURRENT
      DO IX = 1, ICCLD
        IPLON = ICCLDL(IX)

        !  *** Cloudy layer
        ICLDDN(IPLON) = 1
          
        !- total-sky downward flux          
        ODSM = OD(IPLON,IPR,LEV) + TAUCLD(IPLON,LEV,NBI)
        FACTOT1 = ODSM / (BPADE + ODSM)
        BBUTOT1(IPLON,INDEX) = BGLAY(IPLON,IPR) + FACTOT1 * DELBGUP(IPLON,IPR)
        ABSCLDNW = 1.0_dp - TRNCLD(IPLON,NBI,LEV)
        ATOT1(IPLON,INDEX) = ABSS1(IPLON,INDEX) + ABSCLDNW &
                           - ABSS1(IPLON,INDEX) * ABSCLDNW
        BBDTOT = BGLAY(IPLON,IPR) + FACTOT1 * DELBGDN(IPLON,IPR)
        GASSRC = BBD(IPLON,IPR) * ABSS1(IPLON,INDEX)
        !***
        IF (ISTCLDD(IPLON,LEV)  ==  1) THEN
          CLDRADD_L = CLDFRAC(IPLON,LEV) * RADLD1(IPLON,IPR)
          CLRRADD_L = RADLD1(IPLON,IPR) - CLDRADD_L
          RAD_L = 0.0_dp
        ELSE
          CLDRADD_L = CLDRADD(IPLON,IPR)
          CLRRADD_L = CLRRADD(IPLON,IPR)
          RAD_L = RAD(IPLON,IPR)
        ENDIF
        TTOT = 1.0_dp - ATOT1(IPLON,INDEX)
        CLDSRC = BBDTOT * ATOT1(IPLON,INDEX)

        ! Separate RT equations for clear and cloudy streams      
        CLDRADD(IPLON,IPR) = CLDRADD_L * TTOT &
                           + CLDFRAC(IPLON,LEV) * CLDSRC
        CLRRADD(IPLON,IPR) = CLRRADD_L * &
                             (1.0_dp - ABSS1(IPLON,INDEX)) + &
                             (1.0_dp - CLDFRAC(IPLON,LEV)) * GASSRC

        !  Total sky downward radiance
        RADLD1(IPLON,IPR) = CLDRADD(IPLON,IPR) + CLRRADD(IPLON,IPR)
        DRAD1(IPLON) = DRAD1(IPLON) + RADLD1(IPLON,IPR)

        !- clear-sky downward flux          
        RADCLRD1(IPLON,IPR) = RADCLRD1(IPLON,IPR) &
                            + (BBD(IPLON,IPR)-RADCLRD1(IPLON,IPR))*ABSS1(IPLON,INDEX)
        DRADCL1(IPLON) = DRADCL1(IPLON) + RADCLRD1(IPLON,IPR)

        !* Code to account for maximum/random overlap:
        !   Performs RT on the radiance most recently switched between clear and
        !   cloudy streams
        RADMOD = RAD_L * (FACCLR1D(IPLON,LEV-1) * &
               (1.0_dp-ABSS1(IPLON,INDEX)) +                 &
               FACCLD1D(IPLON,LEV-1) *  TTOT) -              &
               FACCMB1D(IPLON,LEV-1) * GASSRC +              &
               FACCMB2D(IPLON,LEV-1) * CLDSRC

        !   Computes what the clear and cloudy streams would have been had no
        !   radiance been switched       
        OLDCLD = CLDRADD(IPLON,IPR) - RADMOD
        OLDCLR = CLRRADD(IPLON,IPR) + RADMOD

        !   Computes the radiance to be switched between clear and cloudy.
        RAD(IPLON,IPR) = -RADMOD + FACCLR2D(IPLON,LEV-1)*OLDCLR - &
                          FACCLD2D(IPLON,LEV-1)*OLDCLD
        CLDRADD(IPLON,IPR) = CLDRADD(IPLON,IPR) + RAD(IPLON,IPR)
        CLRRADD(IPLON,IPR) = CLRRADD(IPLON,IPR) - RAD(IPLON,IPR)

      ENDDO
    ENDDO

    DO IPR = 1, JPGPT
      NBI = NGB(IPR)
      INDEX = IENT + IPR


      DO IPLON=1,KPROMA
        IF (ICLDLYR(IPLON,LEV) .NE. 1) THEN
          !  *** Clear layer
          !  *** DRAD1 holds summed radiance for total sky stream
          !  *** DRADCL1 holds summed radiance for clear sky stream

          !- total-sky downward flux          
          RADLD1(IPLON,IPR) = RADLD1(IPLON,IPR) &
                            + (BBD(IPLON,IPR)-RADLD1(IPLON,IPR))*ABSS1(IPLON,INDEX)
          DRAD1(IPLON) = DRAD1(IPLON) + RADLD1(IPLON,IPR)
          !- clear-sky downward flux          
          !-  Set clear sky stream to total sky stream as long as layers
          !-  remain clear.  Streams diverge when a cloud is reached.
          IF (ICLDDN(IPLON) == 1) THEN
            RADCLRD1(IPLON,IPR) = RADCLRD1(IPLON,IPR) &
                                + (BBD(IPLON,IPR)-RADCLRD1(IPLON,IPR))*ABSS1(IPLON,INDEX)
            DRADCL1(IPLON) = DRADCL1(IPLON) + RADCLRD1(IPLON,IPR)
          ELSE
            RADCLRD1(IPLON,IPR) = RADLD1(IPLON,IPR)
            DRADCL1(IPLON) = DRAD1(IPLON)
          ENDIF
          ENDIF

      ENDDO
    ENDDO

    TOTDFLUC(1:KPROMA,LEV-1) = DRADCL1(1:KPROMA) * WTNUM(1)
    TOTDFLUX(1:KPROMA,LEV-1) = DRAD1(1:KPROMA)   * WTNUM(1)

  ENDDO


  ! Spectral reflectivity and reflectance
  ! Includes the contribution of spectrally varying longwave emissivity 
  ! and reflection from the surface to the upward radiative transfer.
  ! Note: Spectral and Lambertian reflections are identical for the one
  ! angle flux integration used here.

  URAD1(1:KPROMA)   = 0.0_dp
  URADCL1(1:KPROMA) = 0.0_dp

  !IF (IREFLECT  ==  0) THEN
  !- Lambertian reflection.
  DO IPR = 1, JPGPT
    DO IPLON = 1, KPROMA
      ! Clear-sky radiance
      !      RADCLD = 2.0_dp * (RADCLRD1(IPLON,IPR) * WTNUM(1) )
      RADCLD = RADCLRD1(IPLON,IPR)
      RADCLRU1(IPLON,IPR) = RADUEMIT(IPLON,IPR) &
                          + (1.0_dp - SEMIS(IPLON,IPR)) * RADCLD
      URADCL1(IPLON) = URADCL1(IPLON) + RADCLRU1(IPLON,IPR)

      ! Total sky radiance
      !      RADD = 2.0_dp * (RADLD1(IPLON,IPR) * WTNUM(1) )
      RADD = RADLD1(IPLON,IPR)
      RADLU1(IPLON,IPR) = RADUEMIT(IPLON,IPR) &
                        + (1.0_dp - SEMIS(IPLON,IPR)) * RADD
      URAD1(IPLON) = URAD1(IPLON) + RADLU1(IPLON,IPR)
    ENDDO
  ENDDO
  TOTUFLUC(1:KPROMA,0) = URADCL1(1:KPROMA) * 0.5_dp
  TOTUFLUX(1:KPROMA,0) = URAD1(1:KPROMA) * 0.5_dp
  !ELSE
  !- Specular reflection.
  !  DO IPR = 1, JPGPT
  !    DO IPLON = 1, KPROMA
  !      RADCLU = RADUEMIT(IPLON,IPR)
  !      RADCLRU1(IPLON,IPR) = RADCLU + (1.0_dp - SEMIS(IPLON,IPR)) * RADCLRD1(IPLON,IPR)
  !      URADCL1(IPLON) = URADCL1(IPLON) + RADCLRU1(IPLON,IPR)
  !
  !      RADU = RADUEMIT(IPLON,IPR)
  !      RADLU1(IPLON,IPR) = RADU + (1.0_dp - SEMIS(IPLON,IPR)) * RADLD1(IPLON,IPR)
  !      URAD1(IPLON) = URAD1(IPLON) + RADLU1(IPLON,IPR)
  !    ENDDO
  !  ENDDO
  !  TOTUFLUC(1:KPROMA,0) = URADCL1(1:KPROMA) * WTNUM(1)
  !  TOTUFLUX(1:KPROMA,0) = URAD1(1:KPROMA)   * WTNUM(1)
  !ENDIF


  !- Upward radiative transfer.
  !- *** URAD1 holds the summed radiance for total sky stream
  !- *** URADCL1 holds the summed radiance for clear sky stream
  DO LEV = 1, KLEV

    URAD1(1:KPROMA)   = 0.0_dp
    URADCL1(1:KPROMA) = 0.0_dp

    ICCLD = 0
    DO IPLON = 1, KPROMA
      IF (ICLDLYR(IPLON,LEV) == 1) THEN
        ICCLD = ICCLD + 1
        ICCLDL(ICCLD) = IPLON
      ENDIF
    END DO

    IENT = JPGPT * (LEV-1)
    DO IPR = 1, JPGPT
      INDEX = IENT + IPR

!CDIR NODEP
!OCL NOVREC
!DIR$ CONCURRENT
      DO IX = 1, ICCLD
        IPLON = ICCLDL(IX)
        !- *** Cloudy layer
        !- total-sky upward flux          
        GASSRC = BBU1(IPLON,INDEX) * ABSS1(IPLON,INDEX)

        !- If first cloudy layer in sequence, split up radiance into clear and
        !    cloudy streams depending on cloud fraction
        IF (ISTCLD(IPLON,LEV)  ==  1) THEN
          CLDRADU(IPLON,IPR) = CLDFRAC(IPLON,LEV) * RADLU1(IPLON,IPR)
          CLRRADU(IPLON,IPR) = RADLU1(IPLON,IPR) - CLDRADU(IPLON,IPR)
          RAD(IPLON,IPR) = 0.0_dp
        ENDIF
        TTOT = 1.0_dp - ATOT1(IPLON,INDEX)
        TRNS = 1.0_dp - ABSS1(IPLON,INDEX)
        CLDSRC = BBUTOT1(IPLON,INDEX) * ATOT1(IPLON,INDEX)

        !- Separate RT equations for clear and cloudy streams      
        CLDRADU(IPLON,IPR) = CLDRADU(IPLON,IPR) * TTOT &
                           + CLDFRAC(IPLON,LEV) * CLDSRC
        CLRRADU(IPLON,IPR) = CLRRADU(IPLON,IPR) * TRNS &
                           + (1.0_dp - CLDFRAC(IPLON,LEV)) * GASSRC

        !- total sky upward flux
        RADLU1(IPLON,IPR) = CLDRADU(IPLON,IPR) + CLRRADU(IPLON,IPR)
        URAD1(IPLON) = URAD1(IPLON) + RADLU1(IPLON,IPR)

        !* Code to account for maximum/random overlap:
        !   Performs RT on the radiance most recently switched between clear and
        !   cloudy streams
        RADMOD = RAD(IPLON,IPR) * (FACCLR1(IPLON,LEV+1) * TRNS + &
               FACCLD1(IPLON,LEV+1) *  TTOT) - &
               FACCMB1(IPLON,LEV+1) * GASSRC + &
               FACCMB2(IPLON,LEV+1) * CLDSRC

        !   Computes what the clear and cloudy streams would have been had no
        !   radiance been switched       
        OLDCLD = CLDRADU(IPLON,IPR) - RADMOD
        OLDCLR = CLRRADU(IPLON,IPR) + RADMOD

        !   Computes the radiance to be switched between clear and cloudy.
        RAD(IPLON,IPR) = -RADMOD + FACCLR2(IPLON,LEV+1)*OLDCLR - &
                          FACCLD2(IPLON,LEV+1)*OLDCLD
        CLDRADU(IPLON,IPR) = CLDRADU(IPLON,IPR) + RAD(IPLON,IPR)
        CLRRADU(IPLON,IPR) = CLRRADU(IPLON,IPR) - RAD(IPLON,IPR)

      ENDDO
    ENDDO

    DO IPR = 1, JPGPT
      INDEX = IENT + IPR

      DO IPLON = 1, KPROMA
        IF (ICLDLYR(IPLON,LEV) .NE. 1) THEN

          !- *** Clear layer
          !- total-sky upward flux          
          RADLU1(IPLON,IPR) = RADLU1(IPLON,IPR)+(BBU1(IPLON,INDEX) &
                            - RADLU1(IPLON,IPR))*ABSS1(IPLON,INDEX)
          URAD1(IPLON) = URAD1(IPLON) + RADLU1(IPLON,IPR)

        ENDIF

      ENDDO
    ENDDO

!CDIR UNROLL=4
    DO IPR = 1, JPGPT
      DO IPLON = 1, KPROMA
      INDEX = IENT + IPR

        !- clear-sky upward flux
        !   Upward clear and total sky streams must be separate because surface
        !   reflectance is different for each.
        RADCLRU1(IPLON,IPR) = RADCLRU1(IPLON,IPR)+(BBU1(IPLON,INDEX) &
                            - RADCLRU1(IPLON,IPR))*ABSS1(IPLON,INDEX)
        URADCL1(IPLON) = URADCL1(IPLON) + RADCLRU1(IPLON,IPR)

      ENDDO

    ENDDO

    TOTUFLUC(1:KPROMA,LEV) = URADCL1(1:KPROMA) * WTNUM(1)
    TOTUFLUX(1:KPROMA,LEV) = URAD1(1:KPROMA)   * WTNUM(1)

  END DO

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ! End of source code from Uwe Schulzweidas version
  ! rrtm_rtrn1a_140gp.f90@@/main/clean/2

RETURN
END SUBROUTINE RRTM_RTRN1A_140GP
