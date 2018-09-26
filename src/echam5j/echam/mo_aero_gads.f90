MODULE mo_aero_gads

!- Description:
!
!  This module provides the optical properties of the 12 tabulated
!  GADS aerosols for 8 different humidity ranges. The GADS aerosols are:
!
!    Aerosol-components:             radius (micron)        density
!                                    r(0%)   r(80%)  sigma  (g/cm3)
!    1.  water-insoluble       inso  0.471   0.471   2.51    2.0
!    2.  water-soluble         waso  0.0212  0.031   2.24    1.8
!    3.  soot                  soot  0.0118  0.0118  2.00    2.3(1.0)
!    4.  sea-salt (nuc)        ssnm  0.03    0.052   2.03    2.2
!    5.  sea-salt (acc)        ssam  0.209   0.416   2.03    2.2
!    6.  sea-salt (coa)        sscm  1.75    3.49    2.03    2.2
!    7.  mineral (nuc)         minm  0.07    0.07    1.95    2.6
!    8.  mineral (acc)         miam  0.39    0.39    2.00    2.6
!    9.  mineral (coa)         micm  1.90    1.90    2.15    2.6
!    10. mineral transported   mitr  0.5     0.5     2.20    2.6
!    11. sulphate droplets     suso  0.0695  0.118   2.03    1.7
!    12. ammonium sulfate      sulf  0.065   0.097   1.50    1.8
!
!    Coefficients are given for rel. humidity:
!    0%, 50%, 70%, 80%, 90%, 95%, 98%, 99%
!
!  The GADS aerosol types are selected by the given array ndfaer,
!  see comments below in su_aero_gads for more details.
!
!  This module contains:
!
!  A) SW and LW optical properties for all 12 GADS aerosol types for
!     8 different relative humidity ranges
!  B) the subroutine su_aero_gads to initialise the variables of A)
!  C) the function aero_gads to get the lon-lev distribution of
!     the GADS aerosols (in mass mixing ratio) from the xtm1 array
!     of mo_memory_g1a.
!
!  This module replaces the following parts of the ECHAM4.f90 model:
!   - mo_aerosols, GADS part
!   - subroutine suaerx
!
!- Author:
!
!  Marco Giorgetta, MPI, April 2000

!=======================================================================

  USE mo_kind,       ONLY: dp
  USE mo_memory_g1a ,ONLY: xtm1
  USE mo_constants  ,ONLy: g
  USE mo_radiation  ,ONLY: newaer,ndfaer

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tauan,pizan,cgan,caern, &
            su_aero_gads,aero_gads,aero_gads4,fcvaer

!=======================================================================
!
! B)

  ! replaces ECHAM4 mo_aerosols
  ! shortwave
  REAL(dp) :: tauan(8,4,12)  !  extinction cross section in cm2/particle
  REAL(dp) :: pizan(8,4,12)  !  single scattering albedo
  REAL(dp) :: cgan(8,4,12)   !  assymetry factor
  ! longwave:
  REAL(dp) :: caern(8,5,12)  !  extinction cross section in cm2/particle
  REAL(dp) :: fcvaer(12)     !  conversion factor in part/g

!=======================================================================

  CONTAINS

!=======================================================================
!
! C)

SUBROUTINE su_aero_gads  ! former suaerx

  ! Description:
  ! Initialize aerosols, the module that contains the
  ! radiative characteristics of the aerosols.
  !
  ! Method:
  !
  ! Implicit arguments:
  ! module aerosols
  !
  ! Aerosol-components:             radius (micron)        density
  !                                 r(0%)   r(80%)  sigma  (g/cm3)
  ! 1.  water-insoluble       inso  0.471   0.471   2.51    2.0
  ! 2.  water-soluble         waso  0.0212  0.031   2.24    1.8
  ! 3.  soot                  soot  0.0118  0.0118  2.00    2.3(1.0)
  ! 4.  sea-salt (nuc)        ssnm  0.03    0.052   2.03    2.2
  ! 5.  sea-salt (acc)        ssam  0.209   0.416   2.03    2.2
  ! 6.  sea-salt (coa)        sscm  1.75    3.49    2.03    2.2
  ! 7.  mineral (nuc)         minm  0.07    0.07    1.95    2.6
  ! 8.  mineral (acc)         miam  0.39    0.39    2.00    2.6
  ! 9.  mineral (coa)         micm  1.90    1.90    2.15    2.6
  ! 10. mineral transported   mitr  0.5     0.5     2.20    2.6
  ! 11. sulphate droplets     suso  0.0695  0.118   2.03    1.7
  ! 12. ammonium sulfate      sulf  0.065   0.097   1.50    1.8
  !
  ! Coefficients are given for rel. humidity:
  ! 0%, 50%, 70%, 80%, 90%, 95%, 98%, 99%
  !
  ! Extinction/absorption cross sections in cm2/particle
  !
  ! Source: global aerosol dataset (gads),
  ! koepke, hess, schult and shettle (1996),
  ! to appear in "theoretical and applied climate"
  !
  ! Conversion from mass mixing ratio to aerosol amount (part/g)
  ! n(k) = fcvaer * (dp(k)/(10*g)) * r(k)
  ! n(k) in part/cm2 in layer k
  ! dp(k)/(10*g) in g/cm2 (computed in lwu)
  ! r(k) in kg/kg (mass mixing ratio) computed in echam4
  ! fcvaer in part/g
  !
  ! fcvaer=1.e12/[(4*pi/3)*(rad**3)*exp(9*(ln(sigma)**2)/2)*rho]
  !
  ! Where:
  ! rad = geometric mean radius in mu (log-normal distribution)
  ! sigma = standard deviation (log-normal distribution)
  ! rho = specific density in g/cm3
  !
  ! ndfaer(12) : aerosol definition array
  ! contains the number of aerosol component, corresponding
  ! with the gads number. element k=(1,.,12) corresponds
  ! with the aerosol mixing ratio written in paer(.,5+k)
  ! newaer = total number of new aerosol components, computed
  ! within echam4 and written in paer, element 6,.,5+newaer
  ! newaer=0 -> no additional aerosols
  !
  ! Reference:
  ! knmi/mpi documentation =to be done=
  !
  ! Authors:
  !
  ! R. Van Dorland, KNMI, May 1995, original source
  ! R. Van Dorland, KNMI, April 1996, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! M.A. Giorgetta, MPI, may 2000, arrays of SW optical properties expanded to 4 bands
  ! 
  ! for more details see file AUTHORS
  !

  !  Local scalars: 
  REAL(dp) :: pi, radiu3, rlnsig, rlnsig2, volpar, zepext
  INTEGER  :: ih, in, ja
  LOGICAL  :: loa

  !  Local arrays: 
  REAL(dp) :: absswa(8,4,12)  !!$, abslwa(8,5,12)
  REAL(dp) :: asyswa(8,4,12)  !!$, asylwa(8,5,12)
  REAL(dp) :: extswa(8,4,12)     , extlwa(8,5,12)
  REAL(dp) :: rhoaer(12), rmnaer(12), sigaer(12)

  !  Intrinsic functions 
  INTRINSIC ATAN, EXP, LOG

  !  Data statements 

!-- 1. Optical parameters

!-- 1.1 RH = 0%
!       =======

!-- 1.1.1   Shortwave coefficients

  ! extswa = sw-extinction cross section (cm2/part)

  ! Short wave interval 1
  DATA (extswa(1,1,ja),ja=1,12)/0.8403E-07_dp, 0.5129E-10_dp, 0.7433E-11_dp, &
        0.4737E-10_dp, 0.9989E-08_dp, 0.5493E-06_dp, 0.8058E-09_dp, 0.3044E-07_dp, 0.7749E-06_dp, &
        0.5735E-07_dp, 0.8252E-09_dp, 0.1449E-09_dp/
  ! Short wave interval 2
  DATA (extswa(1,2,ja),ja=1,12)/0.8013E-07_dp, 0.7632E-11_dp, 0.1674E-11_dp, &
        0.4936E-11_dp, 0.6247E-08_dp, 0.6236E-06_dp, 0.1849E-09_dp, 0.2824E-07_dp, 0.8592E-06_dp, &
        0.6126E-07_dp, 0.1990E-09_dp, 0.2000E-10_dp/

  ! absswa = sw-absorption cross section (cm2/part)

  ! Short wave interval 1
  DATA (absswa(1,1,ja),ja=1,12)/0.2441E-07_dp, 0.1830E-11_dp, 0.5538E-11_dp, &
        0.1587E-15_dp, 0.9777E-13_dp, 0.3668E-10_dp, 0.5926E-10_dp, 0.6167E-08_dp, 0.2935E-06_dp, &
        0.1412E-07_dp, 0.7231E-16_dp, 0.5392E-11_dp/
  ! Short wave interval 2
  DATA (absswa(1,2,ja),ja=1,12)/0.1315E-07_dp, 0.8173E-12_dp, 0.1520E-11_dp, &
        0.4646E-13_dp, 0.4271E-10_dp, 0.2224E-07_dp, 0.8144E-11_dp, 0.2177E-08_dp, 0.2069E-06_dp, &
        0.6006E-08_dp, 0.2291E-10_dp, 0.3586E-11_dp/

  ! asyswa = sw-asymmetry factor

  ! Short wave interval 1
  DATA (asyswa(1,1,ja),ja=1,12)/0.8431E+00_dp, 0.6269E+00_dp, 0.3703E+00_dp, &
        0.6116E+00_dp, 0.7003E+00_dp, 0.7893E+00_dp, 0.6792E+00_dp, 0.7817E+00_dp, 0.8549E+00_dp, &
        0.8245E+00_dp, 0.7177E+00_dp, 0.5283E+00_dp/
  ! Short wave interval 2
  DATA (asyswa(1,2,ja),ja=1,12)/0.8406E+00_dp, 0.4442E+00_dp, 0.1477E+00_dp, &
        0.3534E+00_dp, 0.6712E+00_dp, 0.7549E+00_dp, 0.4914E+00_dp, 0.6984E+00_dp, 0.8182E+00_dp, &
        0.7206E+00_dp, 0.5404E+00_dp, 0.2325E+00_dp/

!-- 1.1.2 Longwave coefficients

  ! extlwa = lw-extinction cross section (cm2/part)

  ! Long wave interval 1+6
  DATA (extlwa(1,1,ja),ja=1,12)/0.4582E-07_dp, 0.4628E-12_dp, 0.3579E-12_dp, &
        0.2972E-12_dp, 0.1175E-08_dp, 0.6199E-06_dp, 0.1202E-10_dp, 0.1052E-07_dp, 0.8878E-06_dp, &
        0.3358E-07_dp, 0.4457E-10_dp, 0.0000E+00_dp/
  ! Long wave interval 2
  DATA (extlwa(1,2,ja),ja=1,12)/0.4273E-07_dp, 0.5035E-12_dp, 0.1328E-12_dp, &
        0.2735E-12_dp, 0.2310E-09_dp, 0.4774E-06_dp, 0.2067E-10_dp, 0.6756E-08_dp, 0.8534E-06_dp, &
        0.2240E-07_dp, 0.1582E-10_dp, 0.0000E+00_dp/
  ! Long wave interval 3
  DATA (extlwa(1,3,ja),ja=1,12)/0.5310E-07_dp, 0.8876E-12_dp, 0.2192E-12_dp, &
        0.1566E-12_dp, 0.4019E-09_dp, 0.6050E-06_dp, 0.2845E-10_dp, 0.9278E-08_dp, 0.8141E-06_dp, &
        0.2911E-07_dp, 0.7337E-10_dp, 0.0000E+00_dp/
  ! Long wave interval 4
  DATA (extlwa(1,4,ja),ja=1,12)/0.7324E-07_dp, 0.1048E-11_dp, 0.2257E-12_dp, &
        0.1529E-12_dp, 0.5115E-09_dp, 0.6716E-06_dp, 0.3168E-10_dp, 0.1820E-07_dp, 0.9749E-06_dp, &
        0.5227E-07_dp, 0.7312E-10_dp, 0.0000E+00_dp/
  ! Long wave interval 5
  DATA (extlwa(1,5,ja),ja=1,12)/0.4131E-07_dp, 0.5159E-12_dp, 0.8653E-13_dp, &
        0.4815E-12_dp, 0.2991E-09_dp, 0.4515E-06_dp, 0.1552E-10_dp, 0.6227E-08_dp, 0.9415E-06_dp, &
        0.2294E-07_dp, 0.5965E-11_dp, 0.0000E+00_dp/

!!$  ! abslwa = lw-absorption cross section (cm2/part)
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (abslwa(1,1,ja),ja=1,12)/0.1297E-07_dp, 0.3678E-12_dp, 0.3569E-12_dp, &
!!$        0.2541E-12_dp, 0.1415E-09_dp, 0.8443E-07_dp, 0.8333E-11_dp, 0.2282E-08_dp, 0.2920E-06_dp, &
!!$        0.6815E-08_dp, 0.4033E-10_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (abslwa(1,2,ja),ja=1,12)/0.1635E-07_dp, 0.4998E-12_dp, 0.1327E-12_dp, &
!!$        0.2728E-12_dp, 0.1525E-09_dp, 0.1225E-06_dp, 0.2050E-10_dp, 0.5019E-08_dp, 0.4544E-06_dp, &
!!$        0.1502E-07_dp, 0.1555E-10_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (abslwa(1,3,ja),ja=1,12)/0.2137E-07_dp, 0.8576E-12_dp, 0.2191E-12_dp, &
!!$        0.1522E-12_dp, 0.9758E-10_dp, 0.7808E-07_dp, 0.2789E-10_dp, 0.6268E-08_dp, 0.4322E-06_dp, &
!!$        0.1761E-07_dp, 0.7134E-10_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (abslwa(1,4,ja),ja=1,12)/0.3043E-07_dp, 0.9997E-12_dp, 0.2256E-12_dp, &
!!$        0.1473E-12_dp, 0.9979E-10_dp, 0.8620E-07_dp, 0.2925E-10_dp, 0.1026E-07_dp, 0.4731E-06_dp, &
!!$        0.2799E-07_dp, 0.7053E-10_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (abslwa(1,5,ja),ja=1,12)/0.2166E-07_dp, 0.5143E-12_dp, 0.8652E-13_dp, &
!!$        0.4812E-12_dp, 0.2576E-09_dp, 0.1962E-06_dp, 0.1542E-10_dp, 0.4599E-08_dp, 0.4941E-06_dp, &
!!$        0.1528E-07_dp, 0.5896E-11_dp, 0.0000E+00_dp/
!!$
!!$  ! asylwa = lw-asymmetry factor
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (asylwa(1,1,ja),ja=1,12)/0.7771E+00_dp, 0.1893E+00_dp, 0.2827E-01_dp, &
!!$        0.1017E+00_dp, 0.4712E+00_dp, 0.6762E+00_dp, 0.2002E+00_dp, 0.5606E+00_dp, 0.7375E+00_dp, &
!!$        0.5866E+00_dp, 0.2038E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (asylwa(1,2,ja),ja=1,12)/0.6524E+00_dp, 0.8294E-01_dp, 0.5490E-02_dp, &
!!$        0.2416E-01_dp, 0.3209E+00_dp, 0.6895E+00_dp, 0.6468E-01_dp, 0.3716E+00_dp, 0.7186E+00_dp, &
!!$        0.3967E+00_dp, 0.4908E-01_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (asylwa(1,3,ja),ja=1,12)/0.7344E+00_dp, 0.1379E+00_dp, 0.1266E-01_dp, &
!!$        0.5255E-01_dp, 0.4315E+00_dp, 0.7167E+00_dp, 0.1135E+00_dp, 0.4950E+00_dp, 0.8095E+00_dp, &
!!$        0.5643E+00_dp, 0.1008E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (asylwa(1,4,ja),ja=1,12)/0.6692E+00_dp, 0.1597E+00_dp, 0.1326E-01_dp, &
!!$        0.5806E-01_dp, 0.4368E+00_dp, 0.6749E+00_dp, 0.1527E+00_dp, 0.3926E+00_dp, 0.7079E+00_dp, &
!!$        0.4551E+00_dp, 0.1153E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (asylwa(1,5,ja),ja=1,12)/0.5493E+00_dp, 0.5586E-01_dp, 0.2660E-02_dp, &
!!$        0.1324E-01_dp, 0.2328E+00_dp, 0.6017E+00_dp, 0.4756E-01_dp, 0.2371E+00_dp, 0.5267E+00_dp, &
!!$        0.2313E+00_dp, 0.2430E-01_dp, 0.0000E+00_dp/

!-- 1.2 RH = 50%
!       ========

!-- 1.2.1 Shortwave coefficients

  ! extswa = sw-extinction cross section (cm2/part)

  ! Short wave interval 1
  DATA (extswa(2,1,ja),ja=1,12)/0.8403E-07_dp, 0.8274E-10_dp, 0.7433E-11_dp, &
&      0.1658E-09_dp, 0.2449E-07_dp, 0.1398E-05_dp, 0.8058E-09_dp, 0.3044E-07_dp, 0.7749E-06_dp, &
&      0.5735E-07_dp, 0.1832E-08_dp, 0.1671E-09_dp/
  ! Short wave interval 2
  DATA (extswa(2,2,ja),ja=1,12)/0.8013E-07_dp, 0.1353E-10_dp, 0.1674E-11_dp, &
&      0.2448E-10_dp, 0.1875E-07_dp, 0.1534E-05_dp, 0.1849E-09_dp, 0.2824E-07_dp, 0.8592E-06_dp, &
&      0.6126E-07_dp, 0.5293E-09_dp, 0.2361E-10_dp/

  ! absswa = sw-absorption cross section (cm2/part)

  ! Short wave interval 1
  DATA (absswa(2,1,ja),ja=1,12)/0.2441E-07_dp, 0.1805E-11_dp, 0.5538E-11_dp, &
&      0.1602E-15_dp, 0.8013E-13_dp, 0.3480E-10_dp, 0.5926E-10_dp, 0.6167E-08_dp, 0.2935E-06_dp, &
&      0.1412E-07_dp, 0.1515E-15_dp, 0.5416E-11_dp/
  ! Short wave interval 2
  DATA (absswa(2,2,ja),ja=1,12)/0.1315E-07_dp, 0.1448E-11_dp, 0.1520E-11_dp, &
&      0.2142E-11_dp, 0.1030E-08_dp, 0.1706E-06_dp, 0.8144E-11_dp, 0.2177E-08_dp, 0.2069E-06_dp, &
&      0.6006E-08_dp, 0.4720E-10_dp, 0.3753E-11_dp/

  ! asyswa = sw-asymmetry factor

  ! Short wave interval 1
  DATA (asyswa(2,1,ja),ja=1,12)/0.8431E+00_dp, 0.6838E+00_dp, 0.3703E+00_dp, &
&      0.7181E+00_dp, 0.7782E+00_dp, 0.8065E+00_dp, 0.6792E+00_dp, 0.7817E+00_dp, 0.8549E+00_dp, &
&      0.8245E+00_dp, 0.7695E+00_dp, 0.5503E+00_dp/
  ! Short wave interval 2
  DATA (asyswa(2,2,ja),ja=1,12)/0.8406E+00_dp, 0.4922E+00_dp, 0.1477E+00_dp, &
&      0.4613E+00_dp, 0.7595E+00_dp, 0.8360E+00_dp, 0.4914E+00_dp, 0.6984E+00_dp, 0.8182E+00_dp, &
&      0.7206E+00_dp, 0.6182E+00_dp, 0.2488E+00_dp/

!-- 1.2.2 Longwave coefficients

  ! extlwa = lw-extinction cross section (cm2/part)

  ! Long wave interval 1+6
  DATA (extlwa(2,1,ja),ja=1,12)/0.4582E-07_dp, 0.9593E-12_dp, 0.3579E-12_dp, &
&      0.1898E-11_dp, 0.5181E-08_dp, 0.1592E-05_dp, 0.1202E-10_dp, 0.1052E-07_dp, 0.8878E-06_dp, &
&      0.3358E-07_dp, 0.7401E-10_dp, 0.0000E+00_dp/
  ! Long wave interval 2
  DATA (extlwa(2,2,ja),ja=1,12)/0.4273E-07_dp, 0.2065E-11_dp, 0.1328E-12_dp, &
&      0.6208E-11_dp, 0.3555E-08_dp, 0.1340E-05_dp, 0.2067E-10_dp, 0.6756E-08_dp, 0.8534E-06_dp, &
&      0.2240E-07_dp, 0.7754E-10_dp, 0.0000E+00_dp/
  ! Long wave interval 3
  DATA (extlwa(2,3,ja),ja=1,12)/0.5310E-07_dp, 0.1635E-11_dp, 0.2192E-12_dp, &
&      0.2298E-11_dp, 0.2165E-08_dp, 0.1400E-05_dp, 0.2845E-10_dp, 0.9278E-08_dp, 0.8141E-06_dp, &
&      0.2911E-07_dp, 0.1054E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 4
  DATA (extlwa(2,4,ja),ja=1,12)/0.7324E-07_dp, 0.1784E-11_dp, 0.2257E-12_dp, &
&      0.1406E-11_dp, 0.1943E-08_dp, 0.1559E-05_dp, 0.3168E-10_dp, 0.1820E-07_dp, 0.9749E-06_dp, &
&      0.5227E-07_dp, 0.1073E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 5
  DATA (extlwa(2,5,ja),ja=1,12)/0.4131E-07_dp, 0.1425E-11_dp, 0.8653E-13_dp, &
&      0.3925E-11_dp, 0.2425E-08_dp, 0.1390E-05_dp, 0.1552E-10_dp, 0.6227E-08_dp, 0.9415E-06_dp, &
&      0.2294E-07_dp, 0.4072E-10_dp, 0.0000E+00_dp/

!!$  ! abslwa = lw-absorption cross section (cm2/part)
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (abslwa(2,1,ja),ja=1,12)/0.1297E-07_dp, 0.7242E-12_dp, 0.3569E-12_dp, &
!!$&      0.1517E-11_dp, 0.9292E-09_dp, 0.3974E-06_dp, 0.8333E-11_dp, 0.2282E-08_dp, 0.2920E-06_dp, &
!!$&      0.6815E-08_dp, 0.5515E-10_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (abslwa(2,2,ja),ja=1,12)/0.1635E-07_dp, 0.2057E-11_dp, 0.1327E-12_dp, &
!!$&      0.6199E-11_dp, 0.3144E-08_dp, 0.7801E-06_dp, 0.2050E-10_dp, 0.5019E-08_dp, 0.4544E-06_dp, &
!!$&      0.1502E-07_dp, 0.7667E-10_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (abslwa(2,3,ja),ja=1,12)/0.2137E-07_dp, 0.1585E-11_dp, 0.2191E-12_dp, &
!!$&      0.2271E-11_dp, 0.1234E-08_dp, 0.4693E-06_dp, 0.2789E-10_dp, 0.6268E-08_dp, 0.4322E-06_dp, &
!!$&      0.1761E-07_dp, 0.1021E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (abslwa(2,4,ja),ja=1,12)/0.3043E-07_dp, 0.1704E-11_dp, 0.2256E-12_dp, &
!!$&      0.1375E-11_dp, 0.8316E-09_dp, 0.4074E-06_dp, 0.2925E-10_dp, 0.1026E-07_dp, 0.4731E-06_dp, &
!!$&      0.2799E-07_dp, 0.1029E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (abslwa(2,5,ja),ja=1,12)/0.2166E-07_dp, 0.1422E-11_dp, 0.8652E-13_dp, &
!!$&      0.3921E-11_dp, 0.2133E-08_dp, 0.7712E-06_dp, 0.1542E-10_dp, 0.4599E-08_dp, 0.4941E-06_dp, &
!!$&      0.1528E-07_dp, 0.4043E-10_dp, 0.0000E+00_dp/
!!$
!!$  ! asylwa = lw-asymmetry factor
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (asylwa(2,1,ja),ja=1,12)/0.7771E+00_dp, 0.2211E+00_dp, 0.2827E-01_dp, &
!!$&      0.1628E+00_dp, 0.5745E+00_dp, 0.7846E+00_dp, 0.2002E+00_dp, 0.5606E+00_dp, 0.7375E+00_dp, &
!!$&      0.5866E+00_dp, 0.2601E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (asylwa(2,2,ja),ja=1,12)/0.6524E+00_dp, 0.8795E-01_dp, 0.5490E-02_dp, &
!!$&      0.4291E-01_dp, 0.3705E+00_dp, 0.8112E+00_dp, 0.6468E-01_dp, 0.3716E+00_dp, 0.7186E+00_dp, &
!!$&      0.3967E+00_dp, 0.5209E-01_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (asylwa(2,3,ja),ja=1,12)/0.7344E+00_dp, 0.1540E+00_dp, 0.1266E-01_dp, &
!!$&      0.8921E-01_dp, 0.5479E+00_dp, 0.8700E+00_dp, 0.1135E+00_dp, 0.4950E+00_dp, 0.8095E+00_dp, &
!!$&      0.5643E+00_dp, 0.1229E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (asylwa(2,4,ja),ja=1,12)/0.6692E+00_dp, 0.1716E+00_dp, 0.1326E-01_dp, &
!!$&      0.9582E-01_dp, 0.5682E+00_dp, 0.8538E+00_dp, 0.1527E+00_dp, 0.3926E+00_dp, 0.7079E+00_dp, &
!!$&      0.4551E+00_dp, 0.1344E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (asylwa(2,5,ja),ja=1,12)/0.5493E+00_dp, 0.6162E-01_dp, 0.2660E-02_dp, &
!!$&      0.2622E-01_dp, 0.2907E+00_dp, 0.7183E+00_dp, 0.4756E-01_dp, 0.2371E+00_dp, 0.5267E+00_dp, &
!!$&      0.2313E+00_dp, 0.2648E-01_dp, 0.0000E+00_dp/

!-- 1.3 RH = 70%
!       ========

!-- 1.3.1 Shortwave coefficients

  ! extswa = sw-extinction cross section (cm2/part)

  ! Short wave interval 1
  DATA (extswa(3,1,ja),ja=1,12)/0.8403E-07_dp, 0.1020E-09_dp, 0.7433E-11_dp, &
&      0.2307E-09_dp, 0.3053E-07_dp, 0.1766E-05_dp, 0.8058E-09_dp, 0.3044E-07_dp, 0.7749E-06_dp, &
&      0.5735E-07_dp, 0.2339E-08_dp, 0.2728E-09_dp/
  ! Short wave interval 2
  DATA (extswa(3,2,ja),ja=1,12)/0.8013E-07_dp, 0.1740E-10_dp, 0.1674E-11_dp, &
&      0.3689E-10_dp, 0.2479E-07_dp, 0.1921E-05_dp, 0.1849E-09_dp, 0.2824E-07_dp, 0.8592E-06_dp, &
&      0.6126E-07_dp, 0.7293E-09_dp, 0.4174E-10_dp/

  ! absswa = sw-absorption cross section (cm2/part)

  ! Short wave interval 1
  DATA (absswa(3,1,ja),ja=1,12)/0.2441E-07_dp, 0.1805E-11_dp, 0.5538E-11_dp, &
&      0.1652E-15_dp, 0.8282E-13_dp, 0.3621E-10_dp, 0.5926E-10_dp, 0.6167E-08_dp, 0.2935E-06_dp, &
&      0.1412E-07_dp, 0.2015E-15_dp, 0.5514E-11_dp/
  ! Short wave interval 2
  DATA (absswa(3,2,ja),ja=1,12)/0.1315E-07_dp, 0.1830E-11_dp, 0.1520E-11_dp, &
&      0.3299E-11_dp, 0.1465E-08_dp, 0.2247E-06_dp, 0.8144E-11_dp, 0.2177E-08_dp, 0.2069E-06_dp, &
&      0.6006E-08_dp, 0.6110E-10_dp, 0.4339E-11_dp/

  ! asyswa = sw-asymmetry factor

  ! Short wave interval 1
  DATA (asyswa(3,1,ja),ja=1,12)/0.8431E+00_dp, 0.7014E+00_dp, 0.3703E+00_dp, &
&      0.7348E+00_dp, 0.7869E+00_dp, 0.7982E+00_dp, 0.6792E+00_dp, 0.7817E+00_dp, 0.8549E+00_dp, &
&      0.8245E+00_dp, 0.7779E+00_dp, 0.6151E+00_dp/
  ! Short wave interval 2
  DATA (asyswa(3,2,ja),ja=1,12)/0.8406E+00_dp, 0.5106E+00_dp, 0.1477E+00_dp, &
&      0.4875E+00_dp, 0.7710E+00_dp, 0.8439E+00_dp, 0.4914E+00_dp, 0.6984E+00_dp, 0.8182E+00_dp, &
&      0.7206E+00_dp, 0.6471E+00_dp, 0.3058E+00_dp/

!-- 1.3.2 Longwave coefficients

  ! extlwa = lw-extinction cross section (cm2/part)

  ! Long wave interval 1+6
  DATA (extlwa(3,1,ja),ja=1,12)/0.4582E-07_dp, 0.1286E-11_dp, 0.3579E-12_dp, &
&      0.2861E-11_dp, 0.7574E-08_dp, 0.2019E-05_dp, 0.1202E-10_dp, 0.1052E-07_dp, 0.8878E-06_dp, &
&      0.3358E-07_dp, 0.9770E-10_dp, 0.0000E+00_dp/
  ! Long wave interval 2
  DATA (extlwa(3,2,ja),ja=1,12)/0.4273E-07_dp, 0.3081E-11_dp, 0.1328E-12_dp, &
&      0.9478E-11_dp, 0.5423E-08_dp, 0.1747E-05_dp, 0.2067E-10_dp, 0.6756E-08_dp, 0.8534E-06_dp, &
&      0.2240E-07_dp, 0.1158E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 3
  DATA (extlwa(3,3,ja),ja=1,12)/0.5310E-07_dp, 0.2068E-11_dp, 0.2192E-12_dp, &
&      0.3481E-11_dp, 0.3214E-08_dp, 0.1783E-05_dp, 0.2845E-10_dp, 0.9278E-08_dp, 0.8141E-06_dp, &
&      0.2911E-07_dp, 0.1231E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 4
  DATA (extlwa(3,4,ja),ja=1,12)/0.7324E-07_dp, 0.2133E-11_dp, 0.2257E-12_dp, &
&      0.2098E-11_dp, 0.2820E-08_dp, 0.1976E-05_dp, 0.3168E-10_dp, 0.1820E-07_dp, 0.9749E-06_dp, &
&      0.5227E-07_dp, 0.1218E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 5
  DATA (extlwa(3,5,ja),ja=1,12)/0.4131E-07_dp, 0.2011E-11_dp, 0.8653E-13_dp, &
&      0.5806E-11_dp, 0.3685E-08_dp, 0.1828E-05_dp, 0.1552E-10_dp, 0.6227E-08_dp, 0.9415E-06_dp, &
&      0.2294E-07_dp, 0.6254E-10_dp, 0.0000E+00_dp/

!!$  ! abslwa = lw-absorption cross section (cm2/part)
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (abslwa(3,1,ja),ja=1,12)/0.1297E-07_dp, 0.9434E-12_dp, 0.3569E-12_dp, &
!!$&      0.2197E-11_dp, 0.1368E-08_dp, 0.5401E-06_dp, 0.8333E-11_dp, 0.2282E-08_dp, 0.2920E-06_dp, &
!!$&      0.6815E-08_dp, 0.6466E-10_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (abslwa(3,2,ja),ja=1,12)/0.1635E-07_dp, 0.3069E-11_dp, 0.1327E-12_dp, &
!!$&      0.9461E-11_dp, 0.4747E-08_dp, 0.1016E-05_dp, 0.2050E-10_dp, 0.5019E-08_dp, 0.4544E-06_dp, &
!!$&      0.1502E-07_dp, 0.1140E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (abslwa(3,3,ja),ja=1,12)/0.2137E-07_dp, 0.2005E-11_dp, 0.2191E-12_dp, &
!!$&      0.3434E-11_dp, 0.1853E-08_dp, 0.6340E-06_dp, 0.2789E-10_dp, 0.6268E-08_dp, 0.4322E-06_dp, &
!!$&      0.1761E-07_dp, 0.1180E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (abslwa(3,4,ja),ja=1,12)/0.3043E-07_dp, 0.2035E-11_dp, 0.2256E-12_dp, &
!!$&      0.2049E-11_dp, 0.1243E-08_dp, 0.5563E-06_dp, 0.2925E-10_dp, 0.1026E-07_dp, 0.4731E-06_dp, &
!!$&      0.2799E-07_dp, 0.1151E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (abslwa(3,5,ja),ja=1,12)/0.2166E-07_dp, 0.2005E-11_dp, 0.8652E-13_dp, &
!!$&      0.5799E-11_dp, 0.3196E-08_dp, 0.1016E-05_dp, 0.1542E-10_dp, 0.4599E-08_dp, 0.4941E-06_dp, &
!!$&      0.1528E-07_dp, 0.6184E-10_dp, 0.0000E+00_dp/
!!$
!!$  ! asylwa = lw-asymmetry factor
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (asylwa(3,1,ja),ja=1,12)/0.7771E+00_dp, 0.2350E+00_dp, 0.2827E-01_dp, &
!!$&      0.1814E+00_dp, 0.5963E+00_dp, 0.8002E+00_dp, 0.2002E+00_dp, 0.5606E+00_dp, 0.7375E+00_dp, &
!!$&      0.5866E+00_dp, 0.3208E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (asylwa(3,2,ja),ja=1,12)/0.6524E+00_dp, 0.9228E-01_dp, 0.5490E-02_dp, &
!!$&      0.4970E-01_dp, 0.3929E+00_dp, 0.8266E+00_dp, 0.6468E-01_dp, 0.3716E+00_dp, 0.7186E+00_dp, &
!!$&      0.3967E+00_dp, 0.9611E-01_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (asylwa(3,3,ja),ja=1,12)/0.7344E+00_dp, 0.1630E+00_dp, 0.1266E-01_dp, &
!!$&      0.1021E+00_dp, 0.5780E+00_dp, 0.8854E+00_dp, 0.1135E+00_dp, 0.4950E+00_dp, 0.8095E+00_dp, &
!!$&      0.5643E+00_dp, 0.1983E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (asylwa(3,4,ja),ja=1,12)/0.6692E+00_dp, 0.1796E+00_dp, 0.1326E-01_dp, &
!!$&      0.1095E+00_dp, 0.6000E+00_dp, 0.8729E+00_dp, 0.1527E+00_dp, 0.3926E+00_dp, 0.7079E+00_dp, &
!!$&      0.4551E+00_dp, 0.2156E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (asylwa(3,5,ja),ja=1,12)/0.5493E+00_dp, 0.6556E-01_dp, 0.2660E-02_dp, &
!!$&      0.3114E-01_dp, 0.3119E+00_dp, 0.7400E+00_dp, 0.4756E-01_dp, 0.2371E+00_dp, 0.5267E+00_dp, &
!!$&      0.2313E+00_dp, 0.5150E-01_dp, 0.0000E+00_dp/

!-- 1.4 RH = 80%
!       ========

!-- 1.4.1 Shortwave coefficients

  ! extswa = sw-extinction cross section (cm2/part)

  ! Short wave interval 1
  DATA (extswa(4,1,ja),ja=1,12)/0.8403E-07_dp, 0.1226E-09_dp, 0.7433E-11_dp, &
&      0.2989E-09_dp, 0.3643E-07_dp, 0.2135E-05_dp, 0.8058E-09_dp, 0.3044E-07_dp, 0.7749E-06_dp, &
&      0.5735E-07_dp, 0.2812E-08_dp, 0.4309E-09_dp/
  ! Short wave interval 2
  DATA (extswa(4,2,ja),ja=1,12)/0.8013E-07_dp, 0.2170E-10_dp, 0.1674E-11_dp, &
&      0.5086E-10_dp, 0.3101E-07_dp, 0.2309E-05_dp, 0.1849E-09_dp, 0.2824E-07_dp, 0.8592E-06_dp, &
&      0.6126E-07_dp, 0.9283E-09_dp, 0.7163E-10_dp/

  ! absswa = sw-absorption cross section (cm2/part)

  ! Short wave interval 1
  DATA (absswa(4,1,ja),ja=1,12)/0.2441E-07_dp, 0.1806E-11_dp, 0.5538E-11_dp, &
&      0.1698E-15_dp, 0.8187E-13_dp, 0.3751E-10_dp, 0.5926E-10_dp, 0.6167E-08_dp, 0.2935E-06_dp, &
&      0.1412E-07_dp, 0.2469E-15_dp, 0.5616E-11_dp/
  ! Short wave interval 2
  DATA (absswa(4,2,ja),ja=1,12)/0.1315E-07_dp, 0.2236E-11_dp, 0.1520E-11_dp, &
&      0.4527E-11_dp, 0.1917E-08_dp, 0.2803E-06_dp, 0.8144E-11_dp, 0.2177E-08_dp, 0.2069E-06_dp, &
&      0.6006E-08_dp, 0.7438E-10_dp, 0.5004E-11_dp/

  ! asyswa = sw-asymmetry factor

  ! Short wave interval 1
  DATA (asyswa(4,1,ja),ja=1,12)/0.8431E+00_dp, 0.7144E+00_dp, 0.3703E+00_dp, &
&      0.7461E+00_dp, 0.7934E+00_dp, 0.7894E+00_dp, 0.6792E+00_dp, 0.7817E+00_dp, 0.8549E+00_dp, &
&      0.8245E+00_dp, 0.7831E+00_dp, 0.6638E+00_dp/
  ! Short wave interval 2
  DATA (asyswa(4,2,ja),ja=1,12)/0.8406E+00_dp, 0.5257E+00_dp, 0.1477E+00_dp, &
&      0.5077E+00_dp, 0.7786E+00_dp, 0.8493E+00_dp, 0.4914E+00_dp, 0.6984E+00_dp, 0.8182E+00_dp, &
&      0.7206E+00_dp, 0.6619E+00_dp, 0.3599E+00_dp/

!-- 1.4.2 Longwave coefficients

  ! extlwa = lw-extinction cross section (cm2/part)

  ! Long wave interval 1+6
  DATA (extlwa(4,1,ja),ja=1,12)/0.4582E-07_dp, 0.1652E-11_dp, 0.3579E-12_dp, &
&      0.3943E-11_dp, 0.1028E-07_dp, 0.2448E-05_dp, 0.1202E-10_dp, 0.1052E-07_dp, 0.8878E-06_dp, &
&      0.3358E-07_dp, 0.1212E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 2
  DATA (extlwa(4,2,ja),ja=1,12)/0.4273E-07_dp, 0.4182E-11_dp, 0.1328E-12_dp, &
&      0.1294E-10_dp, 0.7479E-08_dp, 0.2160E-05_dp, 0.2067E-10_dp, 0.6756E-08_dp, 0.8534E-06_dp, &
&      0.2240E-07_dp, 0.1534E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 3
  DATA (extlwa(4,3,ja),ja=1,12)/0.5310E-07_dp, 0.2517E-11_dp, 0.2192E-12_dp, &
&      0.4737E-11_dp, 0.4409E-08_dp, 0.2179E-05_dp, 0.2845E-10_dp, 0.9278E-08_dp, 0.8141E-06_dp, &
&      0.2911E-07_dp, 0.1390E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 4
  DATA (extlwa(4,4,ja),ja=1,12)/0.7324E-07_dp, 0.2466E-11_dp, 0.2257E-12_dp, &
&      0.2836E-11_dp, 0.3835E-08_dp, 0.2408E-05_dp, 0.3168E-10_dp, 0.1820E-07_dp, 0.9749E-06_dp, &
&      0.5227E-07_dp, 0.1334E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 5
  DATA (extlwa(4,5,ja),ja=1,12)/0.4131E-07_dp, 0.2645E-11_dp, 0.8653E-13_dp, &
&      0.7794E-11_dp, 0.5110E-08_dp, 0.2276E-05_dp, 0.1552E-10_dp, 0.6227E-08_dp, 0.9415E-06_dp, &
&      0.2294E-07_dp, 0.8419E-10_dp, 0.0000E+00_dp/

!!$  ! abslwa = lw-absorption cross section (cm2/part)
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (abslwa(4,1,ja),ja=1,12)/0.1297E-07_dp, 0.1178E-11_dp, 0.3569E-12_dp, &
!!$&      0.2917E-11_dp, 0.1855E-08_dp, 0.6885E-06_dp, 0.8333E-11_dp, 0.2282E-08_dp, 0.2920E-06_dp, &
!!$&      0.6815E-08_dp, 0.7343E-10_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (abslwa(4,2,ja),ja=1,12)/0.1635E-07_dp, 0.4165E-11_dp, 0.1327E-12_dp, &
!!$&      0.1291E-10_dp, 0.6476E-08_dp, 0.1249E-05_dp, 0.2050E-10_dp, 0.5019E-08_dp, 0.4544E-06_dp, &
!!$&      0.1502E-07_dp, 0.1507E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (abslwa(4,3,ja),ja=1,12)/0.2137E-07_dp, 0.2440E-11_dp, 0.2191E-12_dp, &
!!$&      0.4664E-11_dp, 0.2529E-08_dp, 0.8038E-06_dp, 0.2789E-10_dp, 0.6268E-08_dp, 0.4322E-06_dp, &
!!$&      0.1761E-07_dp, 0.1323E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (abslwa(4,4,ja),ja=1,12)/0.3043E-07_dp, 0.2350E-11_dp, 0.2256E-12_dp, &
!!$&      0.2761E-11_dp, 0.1699E-08_dp, 0.7126E-06_dp, 0.2925E-10_dp, 0.1026E-07_dp, 0.4731E-06_dp, &
!!$&      0.2799E-07_dp, 0.1249E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (abslwa(4,5,ja),ja=1,12)/0.2166E-07_dp, 0.2637E-11_dp, 0.8652E-13_dp, &
!!$&      0.7783E-11_dp, 0.4373E-08_dp, 0.1262E-05_dp, 0.1542E-10_dp, 0.4599E-08_dp, 0.4941E-06_dp, &
!!$&      0.1528E-07_dp, 0.8310E-10_dp, 0.0000E+00_dp/
!!$
!!$  ! asylwa = lw-asymmetry factor
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (asylwa(4,1,ja),ja=1,12)/0.7771E+00_dp, 0.2471E+00_dp, 0.2827E-01_dp, &
!!$&      0.1968E+00_dp, 0.6129E+00_dp, 0.8115E+00_dp, 0.2002E+00_dp, 0.5606E+00_dp, 0.7375E+00_dp, &
!!$&      0.5866E+00_dp, 0.3416E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (asylwa(4,2,ja),ja=1,12)/0.6524E+00_dp, 0.9676E-01_dp, 0.5490E-02_dp, &
!!$&      0.5568E-01_dp, 0.4127E+00_dp, 0.8383E+00_dp, 0.6468E-01_dp, 0.3716E+00_dp, 0.7186E+00_dp, &
!!$&      0.3967E+00_dp, 0.1099E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (asylwa(4,3,ja),ja=1,12)/0.7344E+00_dp, 0.1717E+00_dp, 0.1266E-01_dp, &
!!$&      0.1133E+00_dp, 0.6021E+00_dp, 0.8955E+00_dp, 0.1135E+00_dp, 0.4950E+00_dp, 0.8095E+00_dp, &
!!$&      0.5643E+00_dp, 0.2212E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (asylwa(4,4,ja),ja=1,12)/0.6692E+00_dp, 0.1878E+00_dp, 0.1326E-01_dp, &
!!$&      0.1214E+00_dp, 0.6249E+00_dp, 0.8851E+00_dp, 0.1527E+00_dp, 0.3926E+00_dp, 0.7079E+00_dp, &
!!$&      0.4551E+00_dp, 0.2390E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (asylwa(4,5,ja),ja=1,12)/0.5493E+00_dp, 0.6950E-01_dp, 0.2660E-02_dp, &
!!$&      0.3550E-01_dp, 0.3302E+00_dp, 0.7565E+00_dp, 0.4756E-01_dp, 0.2371E+00_dp, 0.5267E+00_dp, &
!!$&      0.2313E+00_dp, 0.6024E-01_dp, 0.0000E+00_dp/

!-- 1.5 RH = 90%
!       ========

!-- 1.5.1 Shortwave coefficients

  ! extswa = sw-extinction cross section (cm2/part)

  ! Short wave interval 1
  DATA (extswa(5,1,ja),ja=1,12)/0.8403E-07_dp, 0.1726E-09_dp, 0.7433E-11_dp, &
&      0.4774E-09_dp, 0.5069E-07_dp, 0.3040E-05_dp, 0.8058E-09_dp, 0.3044E-07_dp, 0.7749E-06_dp, &
&      0.5735E-07_dp, 0.3826E-08_dp, 0.7078E-09_dp/
  ! Short wave interval 2
  DATA (extswa(5,2,ja),ja=1,12)/0.8013E-07_dp, 0.3283E-10_dp, 0.1674E-11_dp, &
&      0.9127E-10_dp, 0.4679E-07_dp, 0.3253E-05_dp, 0.1849E-09_dp, 0.2824E-07_dp, 0.8592E-06_dp, &
&      0.6126E-07_dp, 0.1397E-08_dp, 0.1304E-09_dp/

  ! absswa = sw-absorption cross section (cm2/part)

  ! Short wave interval 1
  DATA (absswa(5,1,ja),ja=1,12)/0.2441E-07_dp, 0.1817E-11_dp, 0.5538E-11_dp, &
&      0.1826E-15_dp, 0.8520E-13_dp, 0.3997E-10_dp, 0.5926E-10_dp, 0.6167E-08_dp, 0.2935E-06_dp, &
&      0.1412E-07_dp, 0.3613E-15_dp, 0.5734E-11_dp/
  ! Short wave interval 2
  DATA (absswa(5,2,ja),ja=1,12)/0.1315E-07_dp, 0.3229E-11_dp, 0.1520E-11_dp, &
&      0.7821E-11_dp, 0.3089E-08_dp, 0.4195E-06_dp, 0.8144E-11_dp, 0.2177E-08_dp, 0.2069E-06_dp, &
&      0.6006E-08_dp, 0.1048E-09_dp, 0.5992E-11_dp/

  ! asyswa = sw-asymmetry factor

  ! Short wave interval 1
  DATA (asyswa(5,1,ja),ja=1,12)/0.8431E+00_dp, 0.7340E+00_dp, 0.3703E+00_dp, &
&      0.7628E+00_dp, 0.8031E+00_dp, 0.7686E+00_dp, 0.6792E+00_dp, 0.7817E+00_dp, 0.8549E+00_dp, &
&      0.8245E+00_dp, 0.7890E+00_dp, 0.7068E+00_dp/
  ! Short wave interval 2
  DATA (asyswa(5,2,ja),ja=1,12)/0.8406E+00_dp, 0.5519E+00_dp, 0.1477E+00_dp, &
&      0.5435E+00_dp, 0.7898E+00_dp, 0.8566E+00_dp, 0.4914E+00_dp, 0.6984E+00_dp, 0.8182E+00_dp, &
&      0.7206E+00_dp, 0.6846E+00_dp, 0.4193E+00_dp/

!-- 1.5.2 Longwave coefficients

  ! extlwa = lw-extinction cross section (cm2/part)

  ! Long wave interval 1+6
  DATA (extlwa(5,1,ja),ja=1,12)/0.4582E-07_dp, 0.2614E-11_dp, 0.3579E-12_dp, &
&      0.7126E-11_dp, 0.1804E-07_dp, 0.3496E-05_dp, 0.1202E-10_dp, 0.1052E-07_dp, 0.8878E-06_dp, &
&      0.3358E-07_dp, 0.1824E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 2
  DATA (extlwa(5,2,ja),ja=1,12)/0.4273E-07_dp, 0.6921E-11_dp, 0.1328E-12_dp, &
&      0.2222E-10_dp, 0.1322E-07_dp, 0.3184E-05_dp, 0.2067E-10_dp, 0.6756E-08_dp, 0.8534E-06_dp, &
&      0.2240E-07_dp, 0.2425E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 3
  DATA (extlwa(5,3,ja),ja=1,12)/0.5310E-07_dp, 0.3593E-11_dp, 0.2192E-12_dp, &
&      0.8119E-11_dp, 0.7913E-08_dp, 0.3177E-05_dp, 0.2845E-10_dp, 0.9278E-08_dp, 0.8141E-06_dp, &
&      0.2911E-07_dp, 0.1762E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 4
  DATA (extlwa(5,4,ja),ja=1,12)/0.7324E-07_dp, 0.3204E-11_dp, 0.2257E-12_dp, &
&      0.4837E-11_dp, 0.6887E-08_dp, 0.3495E-05_dp, 0.3168E-10_dp, 0.1820E-07_dp, 0.9749E-06_dp, &
&      0.5227E-07_dp, 0.1595E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 5
  DATA (extlwa(5,5,ja),ja=1,12)/0.4131E-07_dp, 0.4221E-11_dp, 0.8653E-13_dp, &
&      0.1312E-10_dp, 0.9252E-08_dp, 0.3392E-05_dp, 0.1552E-10_dp, 0.6227E-08_dp, 0.9415E-06_dp, &
&      0.2294E-07_dp, 0.1358E-09_dp, 0.0000E+00_dp/

!!$  ! abslwa = lw-absorption cross section (cm2/part)
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (abslwa(5,1,ja),ja=1,12)/0.1297E-07_dp, 0.1757E-11_dp, 0.3569E-12_dp, &
!!$&      0.4859E-11_dp, 0.3235E-08_dp, 0.1067E-05_dp, 0.8333E-11_dp, 0.2282E-08_dp, 0.2920E-06_dp, &
!!$&      0.6815E-08_dp, 0.9403E-10_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (abslwa(5,2,ja),ja=1,12)/0.1635E-07_dp, 0.6889E-11_dp, 0.1327E-12_dp, &
!!$&      0.2215E-10_dp, 0.1117E-07_dp, 0.1808E-05_dp, 0.2050E-10_dp, 0.5019E-08_dp, 0.4544E-06_dp, &
!!$&      0.1502E-07_dp, 0.2371E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (abslwa(5,3,ja),ja=1,12)/0.2137E-07_dp, 0.3477E-11_dp, 0.2191E-12_dp, &
!!$&      0.7958E-11_dp, 0.4402E-08_dp, 0.1234E-05_dp, 0.2789E-10_dp, 0.6268E-08_dp, 0.4322E-06_dp, &
!!$&      0.1761E-07_dp, 0.1650E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (abslwa(5,4,ja),ja=1,12)/0.3043E-07_dp, 0.3039E-11_dp, 0.2256E-12_dp, &
!!$&      0.4674E-11_dp, 0.2991E-08_dp, 0.1115E-05_dp, 0.2925E-10_dp, 0.1026E-07_dp, 0.4731E-06_dp, &
!!$&      0.2799E-07_dp, 0.1461E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (abslwa(5,5,ja),ja=1,12)/0.2166E-07_dp, 0.4205E-11_dp, 0.8652E-13_dp, &
!!$&      0.1309E-10_dp, 0.7684E-08_dp, 0.1861E-05_dp, 0.1542E-10_dp, 0.4599E-08_dp, 0.4941E-06_dp, &
!!$&      0.1528E-07_dp, 0.1335E-09_dp, 0.0000E+00_dp/
!!$
!!$  ! asylwa = lw-asymmetry factor
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (asylwa(5,1,ja),ja=1,12)/0.7771E+00_dp, 0.2699E+00_dp, 0.2827E-01_dp, &
!!$&      0.2268E+00_dp, 0.6413E+00_dp, 0.8302E+00_dp, 0.2002E+00_dp, 0.5606E+00_dp, 0.7375E+00_dp, &
!!$&      0.5866E+00_dp, 0.3757E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (asylwa(5,2,ja),ja=1,12)/0.6524E+00_dp, 0.1065E+00_dp, 0.5490E-02_dp, &
!!$&      0.6830E-01_dp, 0.4523E+00_dp, 0.8579E+00_dp, 0.6468E-01_dp, 0.3716E+00_dp, 0.7186E+00_dp, &
!!$&      0.3967E+00_dp, 0.1352E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (asylwa(5,3,ja),ja=1,12)/0.7344E+00_dp, 0.1897E+00_dp, 0.1266E-01_dp, &
!!$&      0.1364E+00_dp, 0.6457E+00_dp, 0.9105E+00_dp, 0.1135E+00_dp, 0.4950E+00_dp, 0.8095E+00_dp, &
!!$&      0.5643E+00_dp, 0.2609E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (asylwa(5,4,ja),ja=1,12)/0.6692E+00_dp, 0.2054E+00_dp, 0.1326E-01_dp, &
!!$&      0.1460E+00_dp, 0.6689E+00_dp, 0.9022E+00_dp, 0.1527E+00_dp, 0.3926E+00_dp, 0.7079E+00_dp, &
!!$&      0.4551E+00_dp, 0.2796E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (asylwa(5,5,ja),ja=1,12)/0.5493E+00_dp, 0.7794E-01_dp, 0.2660E-02_dp, &
!!$&      0.4494E-01_dp, 0.3667E+00_dp, 0.7844E+00_dp, 0.4756E-01_dp, 0.2371E+00_dp, 0.5267E+00_dp, &
!!$&      0.2313E+00_dp, 0.7760E-01_dp, 0.0000E+00_dp/

!-- 1.6 RH=95%
!       ======

!-- 1.6.1 Shortwave coefficients

  ! extswa = sw-extinction cross section (cm2/part)

  ! Short wave interval 1
  DATA (extswa(6,1,ja),ja=1,12)/0.8403E-07_dp, 0.2505E-09_dp, 0.7433E-11_dp, &
&      0.7658E-09_dp, 0.7332E-07_dp, 0.4530E-05_dp, 0.8058E-09_dp, 0.3044E-07_dp, 0.7749E-06_dp, &
&      0.5735E-07_dp, 0.5421E-08_dp, 0.1265E-08_dp/
  ! Short wave interval 2
  DATA (extswa(6,2,ja),ja=1,12)/0.8013E-07_dp, 0.5168E-10_dp, 0.1674E-11_dp, &
&      0.1659E-09_dp, 0.7301E-07_dp, 0.4800E-05_dp, 0.1849E-09_dp, 0.2824E-07_dp, 0.8592E-06_dp, &
&      0.6126E-07_dp, 0.2238E-08_dp, 0.2678E-09_dp/

  ! absswa = sw-absorption cross section (cm2/part)

  ! Short wave interval 1
  DATA (absswa(6,1,ja),ja=1,12)/0.2441E-07_dp, 0.1832E-11_dp, 0.5538E-11_dp, &
&      0.2036E-15_dp, 0.1003E-12_dp, 0.4488E-10_dp, 0.5926E-10_dp, 0.6167E-08_dp, 0.2935E-06_dp, &
&      0.1412E-07_dp, 0.5622E-15_dp, 0.5873E-11_dp/
  ! SHORT WAVE INTERVAL 2
  DATA (absswa(6,2,ja),ja=1,12)/0.1315E-07_dp, 0.4798E-11_dp, 0.1520E-11_dp, &
&      0.1343E-10_dp, 0.5117E-08_dp, 0.6597E-06_dp, 0.8144E-11_dp, 0.2177E-08_dp, 0.2069E-06_dp, &
&      0.6006E-08_dp, 0.1583E-09_dp, 0.7776E-11_dp/

  ! asyswa = sw-asymmetry factor

  ! Short wave interval 1
  DATA (asyswa(6,1,ja),ja=1,12)/0.8431E+00_dp, 0.7505E+00_dp, 0.3703E+00_dp, &
&      0.7760E+00_dp, 0.8121E+00_dp, 0.7394E+00_dp, 0.6792E+00_dp, 0.7817E+00_dp, 0.8549E+00_dp, &
&      0.8245E+00_dp, 0.7934E+00_dp, 0.7466E+00_dp/
  ! Short wave interval 2
  DATA (asyswa(6,2,ja),ja=1,12)/0.8406E+00_dp, 0.5784E+00_dp, 0.1477E+00_dp, &
&      0.5787E+00_dp, 0.7989E+00_dp, 0.8602E+00_dp, 0.4914E+00_dp, 0.6984E+00_dp, 0.8182E+00_dp, &
&      0.7206E+00_dp, 0.7076E+00_dp, 0.4889E+00_dp/

!-- 1.6.2 Longwave coefficients

  ! extlwa = lw-extinction cross section (cm2/part)

  ! Long wave interval 1+6
  DATA (extlwa(6,1,ja),ja=1,12)/0.4582E-07_dp, 0.4298E-11_dp, 0.3579E-12_dp, &
&      0.1331E-10_dp, 0.3312E-07_dp, 0.5204E-05_dp, 0.1202E-10_dp, 0.1052E-07_dp, 0.8878E-06_dp, &
&      0.3358E-07_dp, 0.3094E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 2
  DATA (extlwa(6,2,ja),ja=1,12)/0.4273E-07_dp, 0.1131E-10_dp, 0.1328E-12_dp, &
&      0.3802E-10_dp, 0.2414E-07_dp, 0.4877E-05_dp, 0.2067E-10_dp, 0.6756E-08_dp, 0.8534E-06_dp, &
&      0.2240E-07_dp, 0.4059E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 3
  DATA (extlwa(6,3,ja),ja=1,12)/0.5310E-07_dp, 0.5273E-11_dp, 0.2192E-12_dp, &
&      0.1394E-10_dp, 0.1505E-07_dp, 0.4847E-05_dp, 0.2845E-10_dp, 0.9278E-08_dp, 0.8141E-06_dp, &
&      0.2911E-07_dp, 0.2451E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 4
  DATA (extlwa(6,4,ja),ja=1,12)/0.7324E-07_dp, 0.4294E-11_dp, 0.2257E-12_dp, &
&      0.8314E-11_dp, 0.1331E-07_dp, 0.5306E-05_dp, 0.3168E-10_dp, 0.1820E-07_dp, 0.9749E-06_dp, &
&      0.5227E-07_dp, 0.2073E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 5
  DATA (extlwa(6,5,ja),ja=1,12)/0.4131E-07_dp, 0.6747E-11_dp, 0.8653E-13_dp, &
&      0.2221E-10_dp, 0.1756E-07_dp, 0.5240E-05_dp, 0.1552E-10_dp, 0.6227E-08_dp, 0.9415E-06_dp, &
&      0.2294E-07_dp, 0.2317E-09_dp, 0.0000E+00_dp/

!!$  ! abslwa = lw-absorption cross section (cm2/part)
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (abslwa(6,1,ja),ja=1,12)/0.1297E-07_dp, 0.2685E-11_dp, 0.3569E-12_dp, &
!!$&      0.8194E-11_dp, 0.5916E-08_dp, 0.1719E-05_dp, 0.8333E-11_dp, 0.2282E-08_dp, 0.2920E-06_dp, &
!!$&      0.6815E-08_dp, 0.1315E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (abslwa(6,2,ja),ja=1,12)/0.1635E-07_dp, 0.1124E-10_dp, 0.1327E-12_dp, &
!!$&      0.3784E-10_dp, 0.1972E-07_dp, 0.2704E-05_dp, 0.2050E-10_dp, 0.5019E-08_dp, 0.4544E-06_dp, &
!!$&      0.1502E-07_dp, 0.3941E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (abslwa(6,3,ja),ja=1,12)/0.2137E-07_dp, 0.5085E-11_dp, 0.2191E-12_dp, &
!!$&      0.1356E-10_dp, 0.7942E-08_dp, 0.1970E-05_dp, 0.2789E-10_dp, 0.6268E-08_dp, 0.4322E-06_dp, &
!!$&      0.1761E-07_dp, 0.2233E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (abslwa(6,4,ja),ja=1,12)/0.3043E-07_dp, 0.4042E-11_dp, 0.2256E-12_dp, &
!!$&      0.7939E-11_dp, 0.5508E-08_dp, 0.1819E-05_dp, 0.2925E-10_dp, 0.1026E-07_dp, 0.4731E-06_dp, &
!!$&      0.2799E-07_dp, 0.1827E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (abslwa(6,5,ja),ja=1,12)/0.2166E-07_dp, 0.6716E-11_dp, 0.8652E-13_dp, &
!!$&      0.2214E-10_dp, 0.1403E-07_dp, 0.2833E-05_dp, 0.1542E-10_dp, 0.4599E-08_dp, 0.4941E-06_dp, &
!!$&      0.1528E-07_dp, 0.2262E-09_dp, 0.0000E+00_dp/
!!$
!!$  ! asylwa = lw-asymmetry factor
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (asylwa(6,1,ja),ja=1,12)/0.7771E+00_dp, 0.2953E+00_dp, 0.2827E-01_dp, &
!!$&      0.2601E+00_dp, 0.6695E+00_dp, 0.8495E+00_dp, 0.2002E+00_dp, 0.5606E+00_dp, 0.7375E+00_dp, &
!!$&      0.5866E+00_dp, 0.4137E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (asylwa(6,2,ja),ja=1,12)/0.6524E+00_dp, 0.1188E+00_dp, 0.5490E-02_dp, &
!!$&      0.8378E-01_dp, 0.4989E+00_dp, 0.8774E+00_dp, 0.6468E-01_dp, 0.3716E+00_dp, 0.7186E+00_dp, &
!!$&      0.3967E+00_dp, 0.1672E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (asylwa(6,3,ja),ja=1,12)/0.7344E+00_dp, 0.2116E+00_dp, 0.1266E-01_dp, &
!!$&      0.1640E+00_dp, 0.6911E+00_dp, 0.9236E+00_dp, 0.1135E+00_dp, 0.4950E+00_dp, 0.8095E+00_dp, &
!!$&      0.5643E+00_dp, 0.3083E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (asylwa(6,4,ja),ja=1,12)/0.6692E+00_dp, 0.2276E+00_dp, 0.1326E-01_dp, &
!!$&      0.1755E+00_dp, 0.7133E+00_dp, 0.9163E+00_dp, 0.1527E+00_dp, 0.3926E+00_dp, 0.7079E+00_dp, &
!!$&      0.4551E+00_dp, 0.3285E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (asylwa(6,5,ja),ja=1,12)/0.5493E+00_dp, 0.8860E-01_dp, 0.2660E-02_dp, &
!!$&      0.5682E-01_dp, 0.4096E+00_dp, 0.8124E+00_dp, 0.4756E-01_dp, 0.2371E+00_dp, 0.5267E+00_dp, &
!!$&      0.2313E+00_dp, 0.1021E+00_dp, 0.0000E+00_dp/

!-- 1.7 RH = 98%
  ! ========

!-- 1.7.1 Shortwave coefficients

  ! extswa = sw-extinction cross section (cm2/part)

  ! Short wave interval 1
  DATA (extswa(7,1,ja),ja=1,12)/0.8403E-07_dp, 0.4058E-09_dp, 0.7433E-11_dp, &
&      0.1379E-08_dp, 0.1247E-06_dp, 0.8043E-05_dp, 0.8058E-09_dp, 0.3044E-07_dp, 0.7749E-06_dp, &
&      0.5735E-07_dp, 0.8467E-08_dp, 0.3109E-08_dp/
  ! Short wave interval 2
  DATA (extswa(7,2,ja),ja=1,12)/0.8013E-07_dp, 0.9357E-10_dp, 0.1674E-11_dp, &
&      0.3537E-09_dp, 0.1339E-06_dp, 0.8430E-05_dp, 0.1849E-09_dp, 0.2824E-07_dp, 0.8592E-06_dp, &
&      0.6126E-07_dp, 0.4131E-08_dp, 0.8590E-09_dp/

  ! absswa = sw-absorption cross section (cm2/part)

  ! Short wave interval 1
  DATA (absswa(7,1,ja),ja=1,12)/0.2441E-07_dp, 0.1858E-11_dp, 0.5538E-11_dp, &
&      0.2519E-15_dp, 0.1371E-12_dp, 0.6279E-10_dp, 0.5926E-10_dp, 0.6167E-08_dp, 0.2935E-06_dp, &
&      0.1412E-07_dp, 0.1026E-14_dp, 0.6085E-11_dp/
  ! Short wave interval 2
  DATA (absswa(7,2,ja),ja=1,12)/0.1315E-07_dp, 0.8036E-11_dp, 0.1520E-11_dp, &
&      0.2645E-10_dp, 0.1019E-07_dp, 0.1257E-05_dp, 0.8144E-11_dp, 0.2177E-08_dp, 0.2069E-06_dp, &
&      0.6006E-08_dp, 0.2774E-09_dp, 0.1346E-10_dp/

  ! asyswa = sw-asymmetry factor

  ! Short wave interval 1
  DATA (asyswa(7,1,ja),ja=1,12)/0.8431E+00_dp, 0.7666E+00_dp, 0.3703E+00_dp, &
&      0.7877E+00_dp, 0.8232E+00_dp, 0.6873E+00_dp, 0.6792E+00_dp, 0.7817E+00_dp, 0.8549E+00_dp, &
&      0.8245E+00_dp, 0.7963E+00_dp, 0.7896E+00_dp/
  ! Short wave interval 2
  DATA (asyswa(7,2,ja),ja=1,12)/0.8406E+00_dp, 0.6101E+00_dp, 0.1477E+00_dp, &
&      0.6212E+00_dp, 0.8089E+00_dp, 0.8574E+00_dp, 0.4914E+00_dp, 0.6984E+00_dp, 0.8182E+00_dp, &
&      0.7206E+00_dp, 0.7333E+00_dp, 0.5940E+00_dp/

!-- 1.7.2 Longwave coefficients

  ! extlwa = lw-extinction cross section (cm2/part)

  ! Long wave interval 1+6
  DATA (extlwa(7,1,ja),ja=1,12)/0.4582E-07_dp, 0.8276E-11_dp, 0.3579E-12_dp, &
&      0.3055E-10_dp, 0.7533E-07_dp, 0.9172E-05_dp, 0.1202E-10_dp, 0.1052E-07_dp, 0.8878E-06_dp, &
&      0.3358E-07_dp, 0.6583E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 2
  DATA (extlwa(7,2,ja),ja=1,12)/0.4273E-07_dp, 0.2050E-10_dp, 0.1328E-12_dp, &
&      0.7506E-10_dp, 0.5463E-07_dp, 0.8857E-05_dp, 0.2067E-10_dp, 0.6756E-08_dp, 0.8534E-06_dp, &
&      0.2240E-07_dp, 0.7933E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 3
  DATA (extlwa(7,3,ja),ja=1,12)/0.5310E-07_dp, 0.8769E-11_dp, 0.2192E-12_dp, &
&      0.2783E-10_dp, 0.3689E-07_dp, 0.8788E-05_dp, 0.2845E-10_dp, 0.9278E-08_dp, 0.8141E-06_dp, &
&      0.2911E-07_dp, 0.4146E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 4
  DATA (extlwa(7,4,ja),ja=1,12)/0.7324E-07_dp, 0.6501E-11_dp, 0.2257E-12_dp, &
&      0.1676E-10_dp, 0.3388E-07_dp, 0.9513E-05_dp, 0.3168E-10_dp, 0.1820E-07_dp, 0.9749E-06_dp, &
&      0.5227E-07_dp, 0.3269E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 5
  DATA (extlwa(7,5,ja),ja=1,12)/0.4131E-07_dp, 0.1204E-10_dp, 0.8653E-13_dp, &
&      0.4356E-10_dp, 0.4252E-07_dp, 0.9559E-05_dp, 0.1552E-10_dp, 0.6227E-08_dp, 0.9415E-06_dp, &
&      0.2294E-07_dp, 0.4637E-09_dp, 0.0000E+00_dp/

!!$  ! abslwa = lw-absorption cross section (cm2/part)
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (abslwa(7,1,ja),ja=1,12)/0.1297E-07_dp, 0.4635E-11_dp, 0.3569E-12_dp, &
!!$&      0.1611E-10_dp, 0.1367E-07_dp, 0.3328E-05_dp, 0.8333E-11_dp, 0.2282E-08_dp, 0.2920E-06_dp, &
!!$&      0.6815E-08_dp, 0.2204E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (abslwa(7,2,ja),ja=1,12)/0.1635E-07_dp, 0.2033E-10_dp, 0.1327E-12_dp, &
!!$&      0.7446E-10_dp, 0.4233E-07_dp, 0.4747E-05_dp, 0.2050E-10_dp, 0.5019E-08_dp, 0.4544E-06_dp, &
!!$&      0.1502E-07_dp, 0.7600E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (abslwa(7,3,ja),ja=1,12)/0.2137E-07_dp, 0.8380E-11_dp, 0.2191E-12_dp, &
!!$&      0.2665E-10_dp, 0.1784E-07_dp, 0.3773E-05_dp, 0.2789E-10_dp, 0.6268E-08_dp, 0.4322E-06_dp, &
!!$&      0.1761E-07_dp, 0.3584E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (abslwa(7,4,ja),ja=1,12)/0.3043E-07_dp, 0.6029E-11_dp, 0.2256E-12_dp, &
!!$&      0.1562E-10_dp, 0.1283E-07_dp, 0.3580E-05_dp, 0.2925E-10_dp, 0.1026E-07_dp, 0.4731E-06_dp, &
!!$&      0.2799E-07_dp, 0.2671E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (abslwa(7,5,ja),ja=1,12)/0.2166E-07_dp, 0.1196E-10_dp, 0.8652E-13_dp, &
!!$&      0.4329E-10_dp, 0.3193E-07_dp, 0.5058E-05_dp, 0.1542E-10_dp, 0.4599E-08_dp, 0.4941E-06_dp, &
!!$&      0.1528E-07_dp, 0.4464E-09_dp, 0.0000E+00_dp/
!!$
!!$  ! asylwa = lw-asymmetry factor
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (asylwa(7,1,ja),ja=1,12)/0.7771E+00_dp, 0.3290E+00_dp, 0.2827E-01_dp, &
!!$&      0.3061E+00_dp, 0.7039E+00_dp, 0.8752E+00_dp, 0.2002E+00_dp, 0.5606E+00_dp, 0.7375E+00_dp, &
!!$&      0.5866E+00_dp, 0.4644E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (asylwa(7,2,ja),ja=1,12)/0.6524E+00_dp, 0.1376E+00_dp, 0.5490E-02_dp, &
!!$&      0.1078E+00_dp, 0.5665E+00_dp, 0.9005E+00_dp, 0.6468E-01_dp, 0.3716E+00_dp, 0.7186E+00_dp, &
!!$&      0.3967E+00_dp, 0.2183E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (asylwa(7,3,ja),ja=1,12)/0.7344E+00_dp, 0.2433E+00_dp, 0.1266E-01_dp, &
!!$&      0.2053E+00_dp, 0.7484E+00_dp, 0.9390E+00_dp, 0.1135E+00_dp, 0.4950E+00_dp, 0.8095E+00_dp, &
!!$&      0.5643E+00_dp, 0.3774E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (asylwa(7,4,ja),ja=1,12)/0.6692E+00_dp, 0.2605E+00_dp, 0.1326E-01_dp, &
!!$&      0.2196E+00_dp, 0.7673E+00_dp, 0.9318E+00_dp, 0.1527E+00_dp, 0.3926E+00_dp, 0.7079E+00_dp, &
!!$&      0.4551E+00_dp, 0.4004E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (asylwa(7,5,ja),ja=1,12)/0.5493E+00_dp, 0.1048E+00_dp, 0.2660E-02_dp, &
!!$&      0.7572E-01_dp, 0.4727E+00_dp, 0.8465E+00_dp, 0.4756E-01_dp, 0.2371E+00_dp, 0.5267E+00_dp, &
!!$&      0.2313E+00_dp, 0.1501E+00_dp, 0.0000E+00_dp/

!-- 1.8 RH = 99%
  ! ========

!-- 1.8.1 Shortwave coefficients

  ! extswa = sw-extinction cross section (cm2/part)

  ! Short wave interval 1
  DATA (extswa(8,1,ja),ja=1,12)/0.8403E-07_dp, 0.5535E-09_dp, 0.7433E-11_dp, &
&      0.1983E-08_dp, 0.1884E-06_dp, 0.1262E-04_dp, 0.8058E-09_dp, 0.3044E-07_dp, 0.7749E-06_dp, &
&      0.5735E-07_dp, 0.1194E-07_dp, 0.5980E-08_dp/
  ! Short wave interval 2
  DATA (extswa(8,2,ja),ja=1,12)/0.8013E-07_dp, 0.1377E-09_dp, 0.1674E-11_dp, &
&      0.5700E-09_dp, 0.2096E-06_dp, 0.1316E-04_dp, 0.1849E-09_dp, 0.2824E-07_dp, 0.8592E-06_dp, &
&      0.6126E-07_dp, 0.6646E-08_dp, 0.2113E-08_dp/

  ! absswa = sw-absorption cross section (cm2/part)

  ! Short wave interval 1
  DATA (absswa(8,1,ja),ja=1,12)/0.2441E-07_dp, 0.1877E-11_dp, 0.5538E-11_dp, &
&      0.3045E-15_dp, 0.1893E-12_dp, 0.9280E-10_dp, 0.5926E-10_dp, 0.6167E-08_dp, 0.2935E-06_dp, &
&      0.1412E-07_dp, 0.1741E-14_dp, 0.6232E-11_dp/
  ! Short wave interval 2
  DATA (absswa(8,2,ja),ja=1,12)/0.1315E-07_dp, 0.1126E-10_dp, 0.1520E-11_dp, &
&      0.4069E-10_dp, 0.1704E-07_dp, 0.2086E-05_dp, 0.8144E-11_dp, 0.2177E-08_dp, 0.2069E-06_dp, &
&      0.6006E-08_dp, 0.4364E-09_dp, 0.2317E-10_dp/

  ! asyswa = sw-asymmetry factor

  ! Short wave interval 1
  DATA (asyswa(8,1,ja),ja=1,12)/0.8431E+00_dp, 0.7745E+00_dp, 0.3703E+00_dp, &
&      0.7928E+00_dp, 0.8306E+00_dp, 0.6428E+00_dp, 0.6792E+00_dp, 0.7817E+00_dp, 0.8549E+00_dp, &
&      0.8245E+00_dp, 0.7971E+00_dp, 0.8071E+00_dp/
  ! Short wave interval 2
  DATA (asyswa(8,2,ja),ja=1,12)/0.8406E+00_dp, 0.6293E+00_dp, 0.1477E+00_dp, &
&      0.6464E+00_dp, 0.8154E+00_dp, 0.8488E+00_dp, 0.4914E+00_dp, 0.6984E+00_dp, 0.8182E+00_dp, &
&      0.7206E+00_dp, 0.7503E+00_dp, 0.6659E+00_dp/

!-- 2.82 Longwave coefficients

  ! extlwa = lw-extinction cross section (cm2/part)

  ! Long wave interval 1+6
  DATA (extlwa(8,1,ja),ja=1,12)/0.4582E-07_dp, 0.1277E-10_dp, 0.3579E-12_dp, &
&      0.5291E-10_dp, 0.1367E-06_dp, 0.1429E-04_dp, 0.1202E-10_dp, 0.1052E-07_dp, 0.8878E-06_dp, &
&      0.3358E-07_dp, 0.1224E-08_dp, 0.0000E+00_dp/
  ! Long wave interval 2
  DATA (extlwa(8,2,ja),ja=1,12)/0.4273E-07_dp, 0.2976E-10_dp, 0.1328E-12_dp, &
&      0.1162E-09_dp, 0.9982E-07_dp, 0.1402E-04_dp, 0.2067E-10_dp, 0.6756E-08_dp, 0.8534E-06_dp, &
&      0.2240E-07_dp, 0.1347E-08_dp, 0.0000E+00_dp/
  ! Long wave interval 3
  DATA (extlwa(8,3,ja),ja=1,12)/0.5310E-07_dp, 0.1232E-10_dp, 0.2192E-12_dp, &
&      0.4355E-10_dp, 0.7212E-07_dp, 0.1389E-04_dp, 0.2845E-10_dp, 0.9278E-08_dp, 0.8141E-06_dp, &
&      0.2911E-07_dp, 0.6695E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 4
  DATA (extlwa(8,4,ja),ja=1,12)/0.7324E-07_dp, 0.8735E-11_dp, 0.2257E-12_dp, &
&      0.2651E-10_dp, 0.6846E-07_dp, 0.1485E-04_dp, 0.3168E-10_dp, 0.1820E-07_dp, 0.9749E-06_dp, &
&      0.5227E-07_dp, 0.5123E-09_dp, 0.0000E+00_dp/
  ! Long wave interval 5
  DATA (extlwa(8,5,ja),ja=1,12)/0.4131E-07_dp, 0.1740E-10_dp, 0.8653E-13_dp, &
&      0.6734E-10_dp, 0.8208E-07_dp, 0.1510E-04_dp, 0.1552E-10_dp, 0.6227E-08_dp, 0.9415E-06_dp, &
&      0.2294E-07_dp, 0.8039E-09_dp, 0.0000E+00_dp/

!!$  ! abslwa = lw-absorption cross section (cm2/part)
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (abslwa(8,1,ja),ja=1,12)/0.1297E-07_dp, 0.6616E-11_dp, 0.3569E-12_dp, &
!!$&      0.2500E-10_dp, 0.2564E-07_dp, 0.5502E-05_dp, 0.8333E-11_dp, 0.2282E-08_dp, 0.2920E-06_dp, &
!!$&      0.6815E-08_dp, 0.3481E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (abslwa(8,2,ja),ja=1,12)/0.1635E-07_dp, 0.2946E-10_dp, 0.1327E-12_dp, &
!!$&      0.1149E-09_dp, 0.7405E-07_dp, 0.7330E-05_dp, 0.2050E-10_dp, 0.5019E-08_dp, 0.4544E-06_dp, &
!!$&      0.1502E-07_dp, 0.1272E-08_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (abslwa(8,3,ja),ja=1,12)/0.2137E-07_dp, 0.1167E-10_dp, 0.2191E-12_dp, &
!!$&      0.4114E-10_dp, 0.3265E-07_dp, 0.6188E-05_dp, 0.2789E-10_dp, 0.6268E-08_dp, 0.4322E-06_dp, &
!!$&      0.1761E-07_dp, 0.5483E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (abslwa(8,4,ja),ja=1,12)/0.3043E-07_dp, 0.7993E-11_dp, 0.2256E-12_dp, &
!!$&      0.2418E-10_dp, 0.2423E-07_dp, 0.5975E-05_dp, 0.2925E-10_dp, 0.1026E-07_dp, 0.4731E-06_dp, &
!!$&      0.2799E-07_dp, 0.3868E-09_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (abslwa(8,5,ja),ja=1,12)/0.2166E-07_dp, 0.1725E-10_dp, 0.8652E-13_dp, &
!!$&      0.6676E-10_dp, 0.5867E-07_dp, 0.7862E-05_dp, 0.1542E-10_dp, 0.4599E-08_dp, 0.4941E-06_dp, &
!!$&      0.1528E-07_dp, 0.7620E-09_dp, 0.0000E+00_dp/
!!$
!!$  ! asylwa = lw-asymmetry factor
!!$
!!$  ! Long wave interval 1+6
!!$  DATA (asylwa(8,1,ja),ja=1,12)/0.7771E+00_dp, 0.3514E+00_dp, 0.2827E-01_dp, &
!!$&      0.3371E+00_dp, 0.7266E+00_dp, 0.8936E+00_dp, 0.2002E+00_dp, 0.5606E+00_dp, 0.7375E+00_dp, &
!!$&      0.5866E+00_dp, 0.5032E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 2
!!$  DATA (asylwa(8,2,ja),ja=1,12)/0.6524E+00_dp, 0.1514E+00_dp, 0.5490E-02_dp, &
!!$&      0.1258E+00_dp, 0.6173E+00_dp, 0.9149E+00_dp, 0.6468E-01_dp, 0.3716E+00_dp, 0.7186E+00_dp, &
!!$&      0.3967E+00_dp, 0.2588E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 3
!!$  DATA (asylwa(8,3,ja),ja=1,12)/0.7344E+00_dp, 0.2657E+00_dp, 0.1266E-01_dp, &
!!$&      0.2353E+00_dp, 0.7859E+00_dp, 0.9501E+00_dp, 0.1135E+00_dp, 0.4950E+00_dp, 0.8095E+00_dp, &
!!$&      0.5643E+00_dp, 0.4312E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 4
!!$  DATA (asylwa(8,4,ja),ja=1,12)/0.6692E+00_dp, 0.2840E+00_dp, 0.1326E-01_dp, &
!!$&      0.2516E+00_dp, 0.8015E+00_dp, 0.9430E+00_dp, 0.1527E+00_dp, 0.3926E+00_dp, 0.7079E+00_dp, &
!!$&      0.4551E+00_dp, 0.4564E+00_dp, 0.0000E+00_dp/
!!$  ! Long wave interval 5
!!$  DATA (asylwa(8,5,ja),ja=1,12)/0.5493E+00_dp, 0.1167E+00_dp, 0.2660E-02_dp, &
!!$&      0.9024E-01_dp, 0.5216E+00_dp, 0.8689E+00_dp, 0.4756E-01_dp, 0.2371E+00_dp, 0.5267E+00_dp, &
!!$&      0.2313E+00_dp, 0.1908E+00_dp, 0.0000E+00_dp/

!-- 2. Data radius, sigma and rho

!-- 2.1 Geometric mean radius (micron)
!       log-normal distribution, rel. hum = 0%

  DATA (rmnaer(ja),ja=1,12)/0.4710_dp, 0.0212_dp, 0.0118_dp, 0.0300_dp, 0.2090_dp, 1.7500_dp, &
&      0.0700_dp, 0.3900_dp, 1.9000_dp, 0.5000_dp, 0.0695_dp, 0.065_dp/

!-- 2.2 Standard deviation log-normal distribution

  DATA (sigaer(ja),ja=1,12)/2.51_dp, 2.24_dp, 2.00_dp, 2.03_dp, 2.03_dp, 2.03_dp, 1.95_dp, 2.00_dp, &
&      2.15_dp, 2.20_dp, 2.03_dp, 1.50_dp/

!-- 2.3 Specific densities in g/cm3 of dry aerosol

  DATA (rhoaer(ja),ja=1,12)/2.0_dp, 1.8_dp, 2.3_dp, 2.2_dp, 2.2_dp, 2.2_dp, 2.6_dp, 2.6_dp, 2.6_dp, 2.6_dp, &
&      1.7_dp, 1.8_dp/


  !  Executable statements 

!-- 3. Computation conversion factor (part/g)

  pi = 4._dp*ATAN(1._dp)
  DO ja = 1, 12
    radiu3 = rmnaer(ja)*rmnaer(ja)*rmnaer(ja)
    rlnsig = LOG(sigaer(ja))
    rlnsig2 = rlnsig*rlnsig
    volpar = 4._dp*pi*radiu3/3._dp*EXP(9._dp*rlnsig2/2._dp)
    fcvaer(ja) = 1.E+12_dp/(volpar*rhoaer(ja))
  END DO


!-- 3.5 Copy SW opyical properties for old band 2 for near-IR to new bands 3 and 4
!       which together with the new band 2 cover the old band 2.

  DO in = 3,4
     absswa(:,in,:)=absswa(:,2,:)
     asyswa(:,in,:)=asyswa(:,2,:)
     extswa(:,in,:)=extswa(:,2,:)
  ENDDO


!-- 4. Computation single scattering albedo

  ! Humidity classes: ih=1,8

  zepext = 1.E-20_dp

  DO ih = 1, 8
    DO ja = 1, 12
      DO in = 1, 4
        loa = extswa(ih,in,ja) > 0._dp
        IF (loa) THEN
          tauan(ih,in,ja) = extswa(ih,in,ja)
          pizan(ih,in,ja) = 1._dp - absswa(ih,in,ja)/extswa(ih,in,ja)
          cgan(ih,in,ja) = asyswa(ih,in,ja)
        ELSE
          tauan(ih,in,ja) = zepext
          pizan(ih,in,ja) = 0.5_dp
          cgan(ih,in,ja) = 0.5_dp
        END IF
      END DO
      DO in = 1, 5
        caern(ih,in,ja) = extlwa(ih,in,ja)
      END DO
    END DO
  END DO

!-- 5. Choose additional aerosol components in *radctl*

  ! Example: soot in paer(,,6)
  !          waso in paer(,,7)
  !          suso in paer(,,8)
  !
  !      -> Then array ndfaer should be written as: ndfaer=3,2,11
  !      -> in subroutine setrad newaer gets the value 3


END SUBROUTINE su_aero_gads

!=======================================================================

FUNCTION aero_gads(jrow,pdp)

!- Author:
!
!  Marco Giorgetta, MPI, April 2000

  INTEGER, INTENT(in)                                 :: jrow
  REAL(dp), INTENT(in), DIMENSION(:,:)                :: pdp
  REAL(dp), DIMENSION(SIZE(pdp,1),SIZE(pdp,2),newaer) :: aero_gads

  INTEGER :: klon,klev,jaer

  klon=SIZE(pdp,1)
  klev=SIZE(pdp,2)

  do jaer=1,newaer
     aero_gads(:,:,jaer)=xtm1(1:klon,1:klev,jaer,jrow) &
          &           *fcvaer(ndfaer(jaer))*pdp/10._dp/g
  enddo

END FUNCTION aero_gads

!=======================================================================

FUNCTION aero_gads4(jrow,pdp)

!- Author:
!
!  Marco Giorgetta, MPI, April 2000
!  Uli Schlese    , MPI, March 2004

  USE mo_memory_g3b        ,ONLY: so4nat, so4all

  INTEGER, INTENT(in)                                 :: jrow
  REAL(dp), INTENT(in), DIMENSION(:,:)                :: pdp
  REAL(dp), DIMENSION(SIZE(pdp,1),SIZE(pdp,2),newaer) :: aero_gads4

  INTEGER :: klon,klev,jaer
  REAL(dp)    :: zconv

  klon=SIZE(pdp,1)
  klev=SIZE(pdp,2)
!
! conversion factor (kg(sulfur)/kg(air)) --> kg(sulfate)/kg(air) = 3.
! mass scattering efficiency for dry sulfate = 3.41 m**2/g
! mass scattering efficiency for ammonium sulfate = 5. m**2/g
! Feichter et al. (1996), Clim. Dyn.
!

  zconv=3._dp*5._dp/3.41_dp

  do jaer=1,newaer
     aero_gads4(:,:,jaer)=(so4all(1:klon,1:klev,jrow)-      &
                            so4nat(1:klon,1:klev,jrow))      &
                         *fcvaer(ndfaer(jaer))*pdp/10._dp/g*zconv
  enddo

END FUNCTION aero_gads4

!=======================================================================

END MODULE mo_aero_gads
