MODULE mo_ssodrag

  ! Description:
  !
  ! Set up parameters for gravity wave drag calculations
  !
  ! Authors:
  !           Martin Miller, ECMWF, Jan 1990
  !           Francois Lott, LMD,   Jul 1999  
  !           Elisa Manzini, MPI,   Aug 2000
  !
  ! References: 
  !     Lott, 1999: Alleviation of stationary biases in a GCM through...
  !                 Monthly Weather Review, 127, pp 788-801.

  USE mo_kind,      ONLY: dp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE

  ! nombre de vrais traceurs
  INTEGER, PARAMETER :: nqmx=2
  INTEGER, PARAMETER :: nbtr=nqmx-2+1/(nqmx-1)

  INTEGER :: nktopg ! Security value for blocked flow level
  INTEGER :: ntop = 1   ! An estimate to qualify the upper levels of
  !                       the model where one wants to impose strees
  !                       profiles
  INTEGER ::  nstra ! no documentation
  !
  ! Parameters depending on model resolution
  !
  REAL(dp) :: gpicmea   ! (PEAK-mean) threshold for activation of scheme
  REAL(dp) :: gstd      ! Standard deviation threshold for activation of scheme
  REAL(dp) :: gkdrag    ! Gravity wave drag coefficient (G in (3),      LOTT 1999)
  REAL(dp) :: gkwake    ! Bluff-body drag coefficient for low level wake 
  !                    (Cd in (2), LOTT 1999)

  !      SET_UP THE "TUNABLE PARAMETERS" OF THE VARIOUS SSO SCHEMES

  REAL(dp), PARAMETER :: gfrcrit = 0.5_dp  ! Critical Non-dimensional mountain 
  !                                    Height (HNC in (1), LOTT 1999)
  REAL(dp), PARAMETER :: grcrit  = 0.25_dp ! Critical Richardson Number (Ric, End of
  !                                   first column p791 of LOTT 1999)
  REAL(dp), PARAMETER :: gklift  = 0.00_dp ! Mountain Lift coefficient
  !                                   (Cl in (4),     LOTT 1999)
  REAL(dp), PARAMETER :: grahilo = 1.00_dp ! Set-up the trapped waves fraction
  !                                  (Beta , End of first column,  LOTT 1999)
  REAL(dp), PARAMETER :: ghmax   = 10000.0_dp  ! Not used
  REAL(dp), PARAMETER :: gvcrit  = 0.1_dp     ! no documentation

  !       SET_UP  VALUES OF SECURITY PARAMETERS

  REAL(dp), PARAMETER :: gsigcr = 0.80_dp    ! Security value for blocked flow depth
  REAL(dp), PARAMETER :: gssec  = 0.0001_dp  ! Security min value for low-level B-V 
  !                                     frequency
  REAL(dp), PARAMETER :: gtsec  = 0.00001_dp ! Security min value for anisotropy 
  !                                     and GW stress.
  REAL(dp), PARAMETER :: gvsec  = 0.10_dp    ! Security min value for ulow

CONTAINS
  !======================================================================
  SUBROUTINE sugwd(klon,klev)

  USE mo_control, ONLY: nvclev, vct, nn

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: klon, klev

  ! local scalar
  INTEGER :: jk
  REAL(dp)    :: zstra, zsigt, zpm1r, zpr

  !          SET THE VALUES OF THE PARAMETERS
  !
  IF (nn == 21) THEN
    gpicmea = 800.0_dp
    gstd    = 400.0_dp
    gkdrag  = 0.8_dp
    gkwake  = 0.8_dp
  ELSE IF(nn == 31) THEN
    gpicmea = 500.0_dp
    gstd    = 200.0_dp
    gkdrag  = 0.8_dp
    gkwake  = 0.8_dp
  ELSE IF(nn == 42) THEN
    gpicmea = 400.0_dp
    gstd    = 100.0_dp
    gkdrag  = 0.8_dp
    gkwake  = 0.8_dp
  ELSE IF(nn == 63) THEN
    gpicmea = 300.0_dp
    gstd    =  75.0_dp
    gkdrag  = 0.8_dp
    gkwake  = 0.8_dp
  ELSE IF(nn == 85) THEN
    gpicmea = 200.0_dp
    gstd    =  50.0_dp
    gkdrag  = 0.8_dp
    gkwake  = 0.8_dp
  ELSE IF(nn == 106) THEN
    gpicmea = 150.0_dp
    gstd    =  30.0_dp
    gkdrag  = 0.8_dp
    gkwake  = 0.8_dp
  ELSE IF(nn == 159) THEN
    gpicmea = 125.0_dp
    gstd    =  20.0_dp
    gkdrag  = 0.8_dp
    gkwake  = 0.8_dp
  ELSE IF(nn == 213) THEN
    gpicmea = 125.0_dp
    gstd    =  20.0_dp
    gkdrag  = 0.8_dp
    gkwake  = 0.8_dp
  ELSE IF(nn == 255) THEN
    gpicmea = 100.0_dp
    gstd    =  10.0_dp
    gkdrag  = 0.8_dp
    gkwake  = 0.8_dp
  ELSE IF(nn == 319) THEN
    gpicmea = 100.0_dp
    gstd    =  10.0_dp
    gkdrag  = 0.8_dp
    gkwake  = 0.8_dp
  ELSE IF(nn == 511) THEN
    gpicmea = 100.0_dp
    gstd    =  10.0_dp
    gkdrag  = 0.8_dp
    gkwake  = 0.8_dp
  ELSE
    CALL finish ('mo_ssodrag', 'Truncation not supported.')
  ENDIF

  ! PRINT *,' DANS SUGWD NLEV=',klev

  !
  zpr=80000.0_dp
  zstra=0.001_dp
  zsigt=0.94_dp
  !old  ZSIGT=0.85_dp
  !
  DO 110 jk=klev,1,-1
    zpm1r = 0.5_dp*(vct(jk)+vct(jk+1)+zpr*(vct(nvclev+jk)+vct(nvclev+jk+1)))
    zpm1r = zpm1r/zpr

    IF(zpm1r >= zsigt)THEN
      nktopg=jk
    ENDIF
    IF(zpm1r >= zstra)THEN
      nstra=jk-1
    ENDIF
110 END DO

!    PRINT *,' DANS SUGWD nktopg=', nktopg
!    PRINT *,' DANS SUGWD nstra=', nstra
!    PRINT *,' DANS SUGWD ntop=', ntop
    !


!    WRITE(unit=6,fmt='('' *** SSO essential constants ***'')')
!    WRITE(unit=6,fmt='('' *** SPECIFIED IN SUGWD ***'')')
!    WRITE(unit=6,fmt='('' Gravity wave ct '',E13.7,'' '')')gkdrag
!    WRITE(unit=6,fmt='('' Trapped/total wave dag '',E13.7,'' '')')    &
!         grahilo
!    WRITE(unit=6,fmt='('' Critical Richardson   = '',E13.7,'' '')')   &
!         grcrit
!    WRITE(unit=6,fmt='('' Critical Froude'',e13.7)') gfrcrit
!    WRITE(unit=6,fmt='('' Low level Wake bluff cte'',e13.7)') gkwake
!    WRITE(unit=6,fmt='('' Low level lift  cte'',e13.7)') gklift

  END SUBROUTINE sugwd
  !======================================================================
END MODULE mo_ssodrag
