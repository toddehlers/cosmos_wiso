MODULE mo_orbit

  USE mo_kind, ONLY: dp
  USE mo_constants, ONLY: pi => api
  USE mo_exception, ONLY : finish

  IMPLICIT NONE

  PUBLIC   :: orbit
  INTERFACE orbit
     MODULE PROCEDURE orbit_vsop87
     MODULE PROCEDURE orbit_pcmdi
  END INTERFACE
  PRIVATE  :: orbit_vsop87
  PRIVATE  :: orbit_pcmdi

  REAL(dp), PARAMETER :: dtor = pi/180_dp
  REAL(dp), PARAMETER :: stor = dtor/3600_dp
  REAL(dp), PARAMETER :: sectot = 1.0_dp/(86400.0_dp*36525.0_dp)

  INTEGER, PARAMETER :: deltattable(187) = (/ &
       1240, 1150, 1060, 980, 910, 850, 790, 740, 700, 650, &
       620,  580,  550, 530, 500, 480, 460, 440, 420, 400, & 
       370,  350,  330, 310, 280, 260, 240, 220, 200, 180, & 
       160,  140,  130, 120, 110, 100,  90,  90,  90,  90, & 
       90,   90,   90,  90, 100, 100, 100, 100, 100, 110, & 
       110,  110,  110, 110, 110, 110, 110, 120, 120, 120, &
       120,  120,  130, 130, 130, 130, 140, 140, 140, 150, &
       150,  150,  150, 160, 160, 160, 160, 160, 170, 170, &
       170,  170,  170, 170, 170, 170, 160, 160, 150, 140, &
       137,  131,  127, 125, 125, 125, 125, 125, 125, 123, &
       120,  114,  106,  96,  86,  75,  66,  60,  57,  56, & 
       57,   59,   62,  65,  68,  71,  73,  75,  77,  78, &
       79,   75,   64,  54,  29,  16, -10, -27, -36, -47, &
       -54,  -52,  -55, -56, -58, -59, -62, -64, -61, -47, &
       -27,    0,   26,  54,  77, 105, 134, 160, 182, 202, &
       212,  224,  235, 239, 243, 240, 239, 239, 237, 240, &
       243,  253,  262, 273, 282, 291, 300, 307, 314, 322, &
       331,  340,  350, 365, 383, 402, 422, 445, 465, 485, &
       505,  522,  538, 549, 558, 569, 580 /)

  REAL(dp), PARAMETER ::  oblpoly(11) = (/ &
       84381.448_dp, -4680.93_dp, -1.55_dp, 1999.25_dp, -51.38_dp, &
         -249.67_dp,  -39.050_dp,  7.12_dp,   27.87_dp,   5.79_dp, &
            2.45_dp /)

  REAL(dp), PARAMETER :: pom(4)    = (/ &
       125.0445222_dp,  -1934.1362608_dp,  0.00207833_dp,  2.220e-6_dp /)
  REAL(dp), PARAMETER :: pmmoon(4) = (/ &
       357.5277233_dp,  35999.0503400_dp, -0.00016030_dp, -3.330e-6_dp /)
  REAL(dp), PARAMETER :: pmsun(4)  = (/ &
       134.9629814_dp, 477198.8673981_dp,  0.00869720_dp,  1.778e-5_dp /)
  REAL(dp), PARAMETER :: pf(4)     = (/ &
       93.2719103_dp, 483202.0175381_dp, -0.00368250_dp,  3.056e-6_dp /)
  REAL(dp), PARAMETER :: pd(4)     = (/ &
       297.8503631_dp, 445267.1114800_dp, -0.00191420_dp,  5.278e-6_dp /)

  REAL(dp), PARAMETER :: abconst = 20.49552_dp*stor

  INTEGER, PARAMETER :: fk5system = 2

  REAL(dp) ::  declination    ! declination of the sun

CONTAINS

  SUBROUTINE sunlonandecc (t, lon, e)

    REAL(dp), INTENT(in)  :: t
    REAL(dp), INTENT(out) :: lon, e

    REAL(dp) :: l, m, c, e2

    l = (280.46646_dp+t*(36000.76983_dp+t*0.0003032_dp))*dtor;
    m = (357.52910_dp+t*(35999.05028_dp-t*0.0001561_dp))*dtor;
    e = 0.016708617_dp-t*(0.000042040_dp+t*0.0000001236_dp);
    e2 = e*e
    ! equation of the center, in terms of e and m
    c = e*(2-0.25_dp*e2)*SIN(m)+1.25_dp*e2*SIN(2*m) &
         +1.0833333333_dp*e*e2*SIN(3*m)
    lon = l+c

  END SUBROUTINE sunlonandecc

  SUBROUTINE aberration (t, obl, ra, decl, system)

    REAL(dp), INTENT(in)  :: t, obl
    REAL(dp), INTENT(inout) :: ra, decl
    INTEGER, OPTIONAL, INTENT(in) :: system

    REAL(dp) :: cosobl, sinobl
    REAL(dp) :: lon, e, pir, tmp
    REAL(dp) :: cosra, sinra
    REAL(dp) :: cosdecl, sindecl
    REAL(dp) :: coslon, sinlon

    cosobl = COS(obl)
    sinobl = SIN(obl)
    cosra = COS(ra)
    sinra = SIN(ra)
    cosdecl = COS(decl)
    sindecl = SIN(decl)

    CALL sunlonandecc (t, lon, e)

    coslon = COS (lon)
    sinlon = SIN (lon)

    IF (PRESENT(system)) THEN  ! FK5 - include the e-terms.
      pir = (102.93735_dp+t*(1.71954_dp+t*0.00046_dp))*dtor
      coslon = coslon-e*COS(pir)
      sinlon = sinlon-e*SIN(pir)
    END IF

    ra   = ra-abconst*(cosra*coslon*cosobl+sinra*sinlon)/cosdecl
    tmp  = coslon*(sinobl*cosdecl-sinra*sindecl*cosobl)+cosra*sindecl*sinlon
    decl = decl-abconst*tmp

  END SUBROUTINE aberration

  SUBROUTINE ecltoequ(l, b, obl, ra, decl)

    REAL(dp), INTENT(in)  :: l, b, obl
    REAL(dp), INTENT(out) :: ra, decl

    REAL(dp) :: sinobl, cosobl
    REAL(dp) :: sinl, cosl
    REAL(dp) :: sinb, cosb

    sinobl = SIN(obl)
    cosobl = COS(obl)
    sinl = SIN(l)
    cosl = COS(l)
    sinb = SIN(b)
    cosb = COS(b)
    ra = ATAN2(cosb*sinl*cosobl-sinb*sinobl,cosb*cosl)
    IF (ra < 0) ra = ra+2*pi
    decl = ASIN(sinb*cosobl+cosb*sinobl*sinl)

  END SUBROUTINE ecltoequ

  SUBROUTINE sphtorect(s, r)

    REAL(dp), INTENT(in)  :: s(3)
    REAL(dp), INTENT(out) :: r(3)

    r(1) = s(3)*COS(s(1))*COS(s(2))
    r(2) = s(3)*SIN(s(1))*COS(s(2))
    r(3) = s(3)*SIN(s(2))

  END SUBROUTINE sphtorect

  SUBROUTINE recttosph (r, s)

    REAL(dp), INTENT(in)  :: r(3)
    REAL(dp), INTENT(out) :: s(3)

    s(1) = ATAN2(r(2),r(1))
    IF (s(1) < 0.0_dp) s(1) = s(1)+2*pi
    s(3) = SQRT(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
    s(2) = ASIN(r(3)/s(3))

  END SUBROUTINE recttosph

  SUBROUTINE heliotogeo(shelio, searth, sgeo)

    REAL(dp), INTENT(in)  :: shelio(3), searth(3)
    REAL(dp), INTENT(out) :: sgeo(3) 

    REAL(dp) :: r(3), rearth(3)
    INTEGER :: i

    CALL sphtorect (shelio, r)
    CALL sphtorect (searth, rearth)
    DO i = 1, 3
      r(i) = r(i)-rearth(i)
    END DO
    CALL recttosph (r, sgeo)

  END SUBROUTINE heliotogeo

  FUNCTION siderealtime(t) RESULT (sidtime)

    REAL(dp) :: sidtime

    REAL(dp), INTENT(in) :: t

    REAL(dp) :: theta

    theta = t*(360.98564736629_dp*36525 &
           +t*(0.000387933_dp-t/38710000))
    sidtime = MODULO (((280.46061837_dp+theta)*dtor), 2*pi)

  END FUNCTION siderealtime

  FUNCTION evalpoly (p, n, x) RESULT (poly)

    REAL(dp) :: poly

    INTEGER, INTENT(in) :: n
    REAL(dp), INTENT(in) :: p(n), x

    INTEGER :: i

    poly = p(n)
    DO i = n-1, 1, -1
      poly = poly*x+p(i)
    END DO

  END FUNCTION evalpoly

  SUBROUTINE nutationconst(t, nutlon, nutobl)

    REAL(dp), INTENT(in)  :: t
    REAL(dp), INTENT(out) :: nutlon, nutobl

    REAL(dp) :: om, msun, mmoon, f, d

    REAL(dp) :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11

    om    = MODULO ((evalpoly (pom,    4, t)*dtor),2*pi)
    msun  = MODULO ((evalpoly (pmmoon, 4, t)*dtor),2*pi)
    mmoon = MODULO ((evalpoly (pmsun,  4, t)*dtor),2*pi)
    f     = MODULO ((evalpoly (pf,     4, t)*dtor),2*pi)
    d     = MODULO ((evalpoly (pd,     4, t)*dtor),2*pi)

    n1 = (-0.01742_dp * t - 17.1996_dp) * SIN(om)            &
         + (-1.3187_dp - 0.00016_dp * t) * SIN(2*f-2*d+2*om) &
         + (-0.2274_dp - 0.00002_dp * t) * SIN(2*f+2*om)     &
         + (0.2062_dp + 0.00002_dp * t) * SIN(2*om)          &
         + (0.1426_dp - 0.00034_dp * t) * SIN(msun)          &
         + (0.0712_dp + 0.00001_dp * t) * SIN(mmoon)
    n2 = + (-0.0517_dp + 0.00012_dp * t) * SIN(msun+2*f-2*d+2*om) &
         + (-0.0386_dp - 0.00004_dp * t) * SIN(2*f+om)            &
         -0.0301_dp * SIN(mmoon+2*f+2*om)                         &
         + (0.0217_dp - 0.00005_dp * t) * SIN(-msun+2*f-2*d+2*om) &
         -0.0158_dp * SIN(mmoon-2*d)                              &
         + (0.0129_dp + 0.00001_dp * t) * SIN(2*f-2*d+om)
    n3 = +0.0123_dp * SIN(-mmoon+2*f+2*om)                &
         +0.0063_dp * SIN(2*d)                            &
         + (0.0063_dp + 0.00001_dp * t) * SIN(mmoon+om)   &
         -0.0059_dp * SIN(-mmoon+2*f+2*d+2*om)            &
         + (-0.0058_dp - 0.00001_dp * t) * SIN(-mmoon+om) &
         -0.0051_dp * SIN(mmoon+2*f+om)
    n4 = +0.0048_dp * SIN(2*mmoon-2*d)      &
         +0.0046_dp * SIN(-2*mmoon+2*f+om)  &
         -0.0038_dp * SIN(2*f+2*d+2*om)     &
         -0.0031_dp * SIN(2*mmoon+2*f+2*om) &
         +0.0029_dp * SIN(2*mmoon)          &
         +0.0029_dp * SIN(mmoon+2*f-2*d+2*om)
    n5 = +0.0026_dp * SIN(2*f)                        &
         -0.0022_dp * SIN(2*f-2*d)                    &
         +0.0021_dp * SIN(-mmoon+2*f+om)              &
         + (0.0017_dp - 0.00001_dp * t) * SIN(2*msun) &
         +0.0016_dp * SIN(-mmoon+2*d+om)              &
         + (-0.0016_dp + 0.00001_dp * t) * SIN(2*msun+2*f-2*d+2*om)
    n6 = -0.0015_dp * SIN(msun+om)           &
         -0.0013_dp * SIN(mmoon-2*d+om)      &
         -0.0012_dp * SIN(-msun+om)          &
         +0.0011_dp * SIN(2*mmoon-2*f)       &
         -0.0010_dp * SIN(-mmoon+2*f+2*d+om) &
         -0.0008_dp * SIN(mmoon+2*f+2*d+2*om)
    n7 = +0.0007_dp * SIN(msun+2*f+2*om)  &
         -0.0007_dp * SIN(mmoon+msun-2*d) &
         -0.0007_dp * SIN(-msun+2*f+2*om) &
         -0.0007_dp * SIN(2*f+2*d+om)     &
         +0.0006_dp * SIN(mmoon+2*d)      &
         +0.0006_dp * SIN(2*mmoon+2*f-2*d+2*om)
    n8 = +0.0006_dp * SIN(mmoon+2*f-2*d+om) &
         -0.0006_dp * SIN(-2*mmoon+2*d+om)  &
         -0.0006_dp * SIN(2*d+om)           &
         +0.0005_dp * SIN(mmoon-msun)       &
         -0.0005_dp * SIN(-msun+2*f-2*d+om) &
         -0.0005_dp * SIN(-2*d+om)
    n9 = -0.0005_dp * SIN(2*mmoon+2*f+om)      &
         -0.0003_dp * SIN(-2*mmoon+2*f+2*om)   &
         +0.0004_dp * SIN(2*mmoon-2*d+om)      &
         +0.0004_dp * SIN(msun+2*f-2*d+om)     &
         -0.0003_dp * SIN(mmoon-msun+2*f+2*om) &
         -0.0003_dp * SIN(-mmoon-msun+2*f+2*d+2*om)
    n10 =-0.0003_dp * SIN(3*mmoon+2*f+2*om)   &
         -0.0003_dp * SIN(-msun+2*f+2*d+2*om) &
         -0.0003_dp * SIN(mmoon-msun-d)       &
         -0.0004_dp * SIN(mmoon-d)            &
         -0.0004_dp * SIN(msun-2*d)           &
         +0.0004_dp * SIN(mmoon-2*f)
    n11 =-0.0004_dp * SIN(d)          &
         -0.0003_dp * SIN(mmoon+msun) &
         +0.0003_dp * SIN(mmoon+2*f)

    nutlon = n1+n2+n3+n4+n5+n6+n7+n8+n9+n10+n11

    n1 =   (0.00089_dp * t + 9.2025_dp) * COS(om)           &
         + (0.5736_dp - 0.00031_dp * t) * COS(2*f-2*d+2*om) &
         + (0.0977_dp - 0.00005_dp * t) * COS(2*f+2*om)     &
         + (-0.0895_dp + 0.00005_dp * t) * COS(2*om)        &
         + (0.0054_dp - 0.00001_dp * t) * COS(msun)         &
         -0.0007_dp * COS(mmoon)
    n2 = + (0.0224_dp - 0.00006_dp * t) * COS(msun+2*f-2*d+2*om)   &
         +0.0200_dp * COS(2*f+om)                               &
         + (0.0129_dp - 0.00001_dp * t) * COS(mmoon+2*f+2*om)      &
         + (-0.0095_dp + 0.00003_dp * t) * COS(-msun+2*f-2*d+2*om) &
         -0.0070_dp * COS(2*f-2*d+om)                           &
         -0.0053_dp * COS(-mmoon+2*f+2*om)
    n3 = -0.0033_dp * COS(mmoon+om)            &
         +0.0026_dp * COS(-mmoon+2*f+2*d+2*om) &
         +0.0032_dp * COS(-mmoon+om)           &
         +0.0027_dp * COS(mmoon+2*f+om)        &
         -0.0024_dp * COS(-2*mmoon+2*f+om)     &
         +0.0016_dp * COS(2*f+2*d+2*om)
    n4 = +0.0013_dp * COS(2*mmoon+2*f+2*om)    &
         -0.0012_dp * COS(mmoon+2*f-2*d+2*om)  &
         -0.0010_dp * COS(-mmoon+2*f+om)       &
         -0.0008_dp * COS(-mmoon+2*d+om)       &
         +0.0007_dp * COS(2*msun+2*f-2*d+2*om) &
         +0.0009_dp * COS(msun+om)
    n5 = +0.0007_dp * COS(mmoon-2*d+om)       &
         +0.0006_dp * COS(-msun+om)           &
         +0.0005_dp * COS(-mmoon+2*f+2*d+om)  &
         +0.0003_dp * COS(mmoon+2*f+2*d+2*om) &
         -0.0003_dp * COS(msun+2*f+2*om)      &
         +0.0003_dp * COS(-msun+2*f+2*om)
    n6 = +0.0003_dp * COS(2*f+2*d+om)           &
         -0.0003_dp * COS(2*mmoon+2*f-2*d+2*om) &
         -0.0003_dp * COS(mmoon+2*f-2*d+om)     &
         +0.0003_dp * COS(-2*mmoon+2*d+om)      &
         +0.0003_dp * COS(2*d+om)               &
         +0.0003_dp * COS(-msun+2*f-2*d+om)
    n7 = +0.0003_dp * COS(-2*d+om) &
         +0.0003_dp * COS(2*mmoon+2*f+om)

    nutobl = n1+n2+n3+n4+n5+n6+n7

    nutlon = nutlon*stor
    nutobl = nutobl*stor

  END SUBROUTINE nutationconst

  FUNCTION obliquity(t) RESULT(obl)

    REAL(dp) :: obl
    REAL(dp), INTENT(in) :: t

    REAL(dp) :: u

    u = 0.01_dp*t

    obl = evalpoly(oblpoly, 11, u)*stor

  END FUNCTION obliquity

  FUNCTION approxdeltat(t) RESULT (deltat)

    REAL(dp) :: deltat
    REAL(dp), INTENT(in) :: t

    REAL(dp) :: y
    INTEGER :: ind

    y = 2000+t*100
    IF (y > 2000.0_dp) THEN
      deltat = 102.3_dp+t*(123.5_dp+t*32.5_dp)
    ELSE 
      IF (y < 1620_dp) THEN
        IF (y < 948_dp) THEN
          deltat =  2715.6_dp+t*(573.36_dp+t*46.5_dp)
        ELSE 
          deltat = 50.6_dp+t*(67.5_dp+t*22.5_dp)
        END IF
      ELSE
        ! interpolate from the above table
        ind = INT((y-1620)/2)+1
        IF (ind > 186) ind = 186
        y = y/2-ind-810
        deltat = 0.1_dp*(deltattable(ind)+(deltattable(ind+1)-deltattable(ind))*y)
      END IF
    END IF

  END FUNCTION approxdeltat

  FUNCTION jdtot (jd) RESULT (t)

    ! convert Julian Day to centuries since J2000.0.

    REAL(dp) :: t ! the T value corresponding to the Julian Day

    REAL(dp), INTENT(in) :: jd

    t = (jd-2451545.0_dp)/36525.0_dp ! the Julian Day to convert

  END FUNCTION jdtot

!----------------------------------------------------------------------

  SUBROUTINE orbit_vsop87 (jde, ra, dec, dis, gha)

    USE mo_vsop87, ONLY: earth_position

    REAL(dp), INTENT(in)  :: jde
    REAL(dp), INTENT(out) :: ra, dec, dis, gha

    REAL(dp) :: t
    REAL(dp) :: obl
    REAL(dp) :: nutlon, nutobl
    REAL(dp) :: l, b, r
    REAL(dp) :: shelio(3), searth(3), sgeo(3)

    ! Preliminary calculations 

    t = jdtot(jde)
    obl = obliquity(t)

    CALL nutationconst (t, nutlon, nutobl)

    ! Main calculations

    CALL earth_position (t, l, b, r)

    shelio(:) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
    searth(:) = (/ l , b, r /)

    CALL heliotogeo (shelio, searth, sgeo)

    CALL ecltoequ (sgeo(1)+nutlon, sgeo(2), obl+nutobl, ra, dec)
    CALL aberration (t, obl, ra, dec, system=fk5system)

    dis  = sgeo(3)
    gha  = siderealtime(t)
    declination = dec

  END SUBROUTINE orbit_vsop87

!----------------------------------------------------------------------

  SUBROUTINE orbit_pcmdi (pvetim, pdisse, pdec, pra)
    !
    ! Description:
    !
    ! Computes the solar position.
    ! 
    ! Method:
    !
    ! This routine computes three orbital parameters depending on the time
    ! of the day as well as of the year (both in radians). Basic equations 
    ! used are the Kepler equation for the eccentric anomaly, and Lacaille's
    ! formula.
    ! 
    ! Input argument:
    !
    ! pvetim - time of the year from the vernal equinox in radians     
    !
    ! Output arguments:
    !
    ! pdisse - radius length in AU 
    !          (ratio of the solar constant to its annual mean) 
    ! pdec   - declination of Sun
    ! pra    - right ascension of Sun  
    !
    ! references:
    !
    ! Monin, A. S.: An Introduction to the Theory of Climate
    !               D. Reidel Publishing Company, Dordrecht, 1986 (pp 10-12).
    ! Meeus, J.: Astronomische Algorithmen, 2ed
    !            Johann Ambrosius Barth, Leipzig, 1994 (pp 199-222).
    !
    ! S. J. Lorenz, Uni Bremen, July 1996, original version  
    ! S. J. Lorenz, Uni Bremen, July 1998, changed
    ! U. Schlese, DKRZ, September 1998, changed
    ! L. Kornblueh, MPI, December 1998, f90 rewrite
    ! L. Kornblueh, MPI, February 2003, precision changes and proper commenting
    !

    USE mo_constants, ONLY : api
    USE mo_radiation, ONLY : cecc, cobld, clonp

    REAL(dp), INTENT(in)  :: pvetim
    REAL(dp), INTENT(out) :: pdisse, pdec, pra

    ! changed from 1.0e-6_dp, due to stability problems in the solution
    ! explained in Meeus.

    REAL(dp), PARAMETER :: ceps = 1.0e-9_dp        

    REAL(dp) :: zoblr, zlonpr, zsqecc, &
         zeve, ztlonpr, zm, zeold, zenew, zeps, zcose,  &
         zdisse, znu, zlambda, zsinde, zdecli

    REAL(dp) :: za, zb, zs, zs0, zsqrt, zz

    INTEGER :: iter

    ! conversion of inclination (obliquity) from degrees to rad

    zoblr = cobld*api/180.0_dp

    ! conversion of longitude of Perihelion from vernal equinox fromm degrees
    ! to rad

    zlonpr = clonp*api/180.0_dp

    ! intermediate variables   

    zsqecc = SQRT((1.0_dp+cecc)/(1.0_dp-cecc))

    ! calculation of eccentric anomaly of vernal equinox (Lacaille's formula)

    zeve = 2*ATAN(TAN(0.5_dp*zlonpr)/zsqecc)

    ! calculation of true anomaly of vernal equinox (Kepler)

    ztlonpr = zeve-cecc*SIN(zeve)

    ! use Newtons method for determing the eccentric anomaly for the 
    ! actual true anomaly. 

    ! true anomaly

    zm = pvetim-ztlonpr

    zeold = zm/(1.0_dp-cecc)
    zenew = zm

    ! for the iteration a first guess of the eccentric anomly is required
    ! for some special cases the original assumption does not converge.
    ! Following Meeus (pp. 213-214) this is covered by the following 
    ! calculation.

    IF ( cecc > 0.975_dp ) THEN  
      IF ( ABS(zm) < 0.52359_dp) THEN            ! M < ~30 deg in radians
        za = (1.0_dp-cecc)/(4.0_dp*cecc+0.5_dp)
        zb = zm/(8.0_dp*cecc+1.0_dp)
        zsqrt = SIGN(SQRT(zb*zb+za*za*za), zb)
        zz = (ABS(zb+zsqrt))**1.5_dp
        zz = SIGN(zz,zb+zsqrt) 
        zs0 = zz-0.5_dp*za
        zs = zs0-(0.078_dp*zs0**5)/(1.0_dp+cecc)
        
        zenew = zm+cecc*(3.0_dp*zs-4.0_dp*zs*zs*zs)
      END IF
    END IF
    
    ! do the Newton iterations

    iter = 0

    DO
      zeps = zeold-zenew

      IF (iter >= 25) THEN
        CALL finish('orbit_pcmdi','Eccentric anomaly not found!')
      END IF

      IF (ABS(zeps) < ceps) EXIT
      
      iter = iter+1
      
      zeold = zenew

      ! the f(x), and f'(x) is a bit to much simplified - please 
      ! take a look in the ECHAM5 documantation for explanation.

      zcose = COS(zenew)
      zenew = (zm+cecc*(SIN(zenew)-zenew*zcose))/(1.0_dp-cecc*zcose)
    END DO

    ! calculation of distance Earth-Sun in AU

    zdisse = (1.0_dp/(1.0_dp-cecc*COS(zenew)))**2

    ! Calculation of the true anomaly

    znu     = 2.0_dp*ATAN(zsqecc*TAN(zenew*0.5_dp))

    zlambda = znu+zlonpr
    zsinde  = SIN(zoblr)*SIN(zlambda)

    ! Calculation of the declination.

    zdecli  = ASIN(zsinde)

    ! finalize return values

    pdisse  = zdisse
    pdec    = zdecli
    pra     = ATAN2(COS(zoblr)*SIN(zlambda), COS(zlambda))

    declination = pdec

  END SUBROUTINE orbit_pcmdi

END MODULE mo_orbit
