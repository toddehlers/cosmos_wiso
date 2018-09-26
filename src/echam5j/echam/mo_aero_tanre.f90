MODULE mo_aero_tanre

!- Description:
!
!  This module contains a global aerosol climatology following Tanre(). 
!
!  This climatology distinguishes spacial distributions of sea, land,
!  urban, and desert aerosols. Further the climatology contains constant
!  background aerosols of tropospheric, stratospheric and volcanic type.
!
!  These aerosols are mixed together to 5 aerosol mixtures:
!    1: land + desert + tropospheric background
!    2: sea
!    3: urban
!    4: volcanic background
!    5: stratospheric background
!
!  There is no annual cycle in the aerosol distribution.
!
!  This module contains:
!
!  A) the constant naer = 5 = number of aerosol mixture types
!  B) spectral coefficients for the horizontal spatial distribution
!     of sea, land, urban and desert aerosols
!  C) variables for the vertical distribution of sea, land, urban
!     and desert aerosols
!  D) SW and LW optical properties for naer aerosol mixtures
!  E) the subroutine su_aero_tanre to initialize the variables of C)
!  F) the subroutine aero_tanre to extract a longitude height section of
!     the naer aerosol mixtures
!  G) the subroutine cleanup_aero_ to deallocate module variables
!
!  A,D,E, and G are public
!
!  This module replaces the following parts of the ECHAM4.f90 model:
!   - mo_rad1
!   - mo_rad2
!   - mo_aerosols, Tanre part
!   - subroutine aerosol
!   - subroutine aerdis
!   - the computation of the longitude height distribution in radint
!
!- Author:
!
!  Marco Giorgetta, MPI, April 1999
!  U. Schulzweida,  MPI, May 2002, blocking (nproma)

!=======================================================================

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: naer, &
            taua,piza,cga,caer, &
            su_aero_tanre,aero_tanre, init_aero, cleanup_aero

!=======================================================================
!
! A)

  INTEGER, PARAMETER :: naer=5 ! number of aerosol mixture types
                               ! 1: land + desert + tropospheric background
                               ! 2: sea
                               ! 3: urban
                               ! 4: volcanic background
                               ! 5: stratospheric background

!=======================================================================
!
! B)

  ! Horizontal distribution of sea, land, urban and desert aerosols
  ! specified by T10 truncated spectral representations.
  !
  ! replaces ECHAM4 module mo_rad1

  ! sea aerosol cosine coefficients
  REAL(dp), PARAMETER, DIMENSION(66) :: caesc= &
     (/ 0.6688E+00_dp, -0.1172E+00_dp, -0.1013E+00_dp,  0.1636E-01_dp, -0.3699E-01_dp, &
        0.1775E-01_dp, -0.9635E-02_dp,  0.1290E-02_dp,  0.4681E-04_dp, -0.9106E-04_dp, &
        0.9355E-04_dp, -0.7076E-01_dp, -0.1782E-01_dp,  0.1856E-01_dp,  0.1372E-01_dp, &
        0.8210E-04_dp,  0.2149E-02_dp,  0.4856E-03_dp,  0.2231E-03_dp,  0.1824E-03_dp, &
        0.1960E-05_dp,  0.2057E-01_dp,  0.2703E-01_dp,  0.2424E-01_dp,  0.9716E-02_dp, &
        0.1312E-02_dp, -0.8846E-03_dp, -0.3347E-03_dp,  0.6231E-04_dp,  0.6397E-04_dp, &
       -0.3341E-02_dp, -0.1295E-01_dp, -0.4598E-02_dp,  0.3242E-03_dp,  0.8122E-03_dp, &
       -0.2975E-03_dp, -0.7757E-04_dp,  0.7793E-04_dp,  0.4455E-02_dp, -0.1584E-01_dp, &
       -0.2551E-02_dp,  0.1174E-02_dp,  0.1335E-04_dp,  0.5112E-04_dp,  0.5605E-04_dp, &
        0.7412E-04_dp,  0.1857E-02_dp, -0.1917E-03_dp,  0.4460E-03_dp,  0.1767E-04_dp, &
       -0.5281E-04_dp, -0.5043E-03_dp,  0.2467E-03_dp, -0.2497E-03_dp, -0.2377E-04_dp, &
       -0.3954E-04_dp,  0.2666E-03_dp, -0.8186E-03_dp, -0.1441E-03_dp, -0.1904E-04_dp, &
        0.3337E-03_dp, -0.1696E-03_dp, -0.2503E-04_dp,  0.1239E-03_dp, -0.9983E-04_dp, &
       -0.5283E-04_dp/)

  ! sea aerosol sine coefficients
  REAL(dp), PARAMETER, DIMENSION(55) :: caess= &
    (/ -0.3374E-01_dp, -0.3247E-01_dp, -0.1012E-01_dp,  0.6002E-02_dp,  0.5190E-02_dp, &
        0.7784E-03_dp, -0.1090E-02_dp,  0.3294E-03_dp,  0.1719E-03_dp, -0.5866E-05_dp, &
       -0.4124E-03_dp, -0.3742E-01_dp, -0.5054E-02_dp,  0.3430E-02_dp,  0.5513E-03_dp, &
       -0.6235E-03_dp,  0.2892E-03_dp, -0.9730E-04_dp,  0.7078E-04_dp, -0.3300E-01_dp, &
        0.5104E-03_dp, -0.2156E-02_dp, -0.3194E-02_dp, -0.5079E-03_dp, -0.5517E-03_dp, &
        0.4632E-04_dp,  0.5369E-04_dp, -0.2731E-01_dp,  0.5126E-02_dp,  0.2241E-02_dp, &
       -0.5789E-03_dp, -0.3048E-03_dp, -0.1774E-03_dp,  0.1946E-05_dp, -0.8247E-02_dp, &
        0.2338E-02_dp,  0.1021E-02_dp,  0.1575E-04_dp,  0.2612E-05_dp,  0.1995E-04_dp, &
       -0.1319E-02_dp,  0.1384E-02_dp, -0.4159E-03_dp, -0.2337E-03_dp,  0.5764E-04_dp, &
        0.1495E-02_dp, -0.3727E-03_dp,  0.6075E-04_dp, -0.4642E-04_dp,  0.5368E-03_dp, &
       -0.7619E-04_dp,  0.3774E-04_dp,  0.1206E-03_dp, -0.4104E-06_dp,  0.2158E-04_dp/)

  ! land aerosol cosine coefficients
  REAL(dp), PARAMETER, DIMENSION(66) :: caelc= &
    (/  0.1542E+00_dp,  0.8245E-01_dp, -0.1879E-03_dp,  0.4864E-02_dp, -0.5527E-02_dp, &
       -0.7966E-02_dp, -0.2683E-02_dp, -0.2011E-02_dp, -0.8889E-03_dp, -0.1058E-03_dp, &
       -0.1614E-04_dp,  0.4206E-01_dp,  0.1912E-01_dp, -0.9476E-02_dp, -0.6780E-02_dp, &
        0.1767E-03_dp, -0.5422E-03_dp, -0.7753E-03_dp, -0.2106E-03_dp, -0.9870E-04_dp, &
       -0.1721E-04_dp, -0.9536E-02_dp, -0.9580E-02_dp, -0.1050E-01_dp, -0.5747E-02_dp, &
       -0.1282E-02_dp,  0.2248E-03_dp,  0.1694E-03_dp, -0.4782E-04_dp, -0.2441E-04_dp, &
        0.5781E-03_dp,  0.6212E-02_dp,  0.1921E-02_dp, -0.1102E-02_dp, -0.8145E-03_dp, &
        0.2497E-03_dp,  0.1539E-03_dp, -0.2538E-04_dp, -0.3993E-02_dp,  0.9777E-02_dp, &
        0.4837E-03_dp, -0.1304E-02_dp,  0.2417E-04_dp, -0.1370E-04_dp, -0.3731E-05_dp, &
        0.1922E-02_dp, -0.5167E-03_dp,  0.4295E-03_dp, -0.1888E-03_dp,  0.2427E-04_dp, &
        0.4012E-04_dp,  0.1529E-02_dp, -0.2120E-03_dp,  0.8166E-04_dp,  0.2579E-04_dp, &
        0.3488E-04_dp,  0.2140E-03_dp,  0.2274E-03_dp, -0.3447E-05_dp, -0.1075E-04_dp, &
       -0.1018E-03_dp,  0.2864E-04_dp,  0.3442E-04_dp, -0.1002E-03_dp,  0.7117E-04_dp, &
        0.2045E-04_dp/)

  ! land aerosol sine coefficients
  REAL(dp), PARAMETER, DIMENSION(55) :: caels= &
    (/  0.1637E-01_dp,  0.1935E-01_dp,  0.1080E-01_dp,  0.2784E-02_dp,  0.1606E-03_dp, &
        0.1860E-02_dp,  0.1263E-02_dp, -0.2707E-03_dp, -0.2290E-03_dp, -0.9761E-05_dp, &
       -0.7317E-02_dp,  0.2465E-01_dp,  0.6799E-02_dp, -0.1913E-02_dp,  0.1382E-02_dp, &
        0.6691E-03_dp,  0.1414E-03_dp,  0.3527E-04_dp, -0.5210E-04_dp,  0.1873E-01_dp, &
        0.2977E-02_dp,  0.4650E-02_dp,  0.2509E-02_dp,  0.3680E-03_dp,  0.1481E-03_dp, &
       -0.6594E-04_dp, -0.5634E-04_dp,  0.1592E-01_dp, -0.1875E-02_dp, -0.1093E-02_dp, &
        0.3022E-03_dp,  0.2625E-03_dp,  0.3252E-04_dp, -0.3803E-04_dp,  0.4218E-02_dp, &
       -0.1843E-02_dp, -0.1351E-02_dp, -0.2952E-03_dp, -0.8171E-05_dp, -0.1473E-04_dp, &
        0.9076E-03_dp, -0.1057E-02_dp,  0.2676E-03_dp,  0.1307E-03_dp, -0.3628E-04_dp, &
       -0.9158E-03_dp,  0.4335E-03_dp,  0.2927E-04_dp,  0.6602E-04_dp, -0.3570E-03_dp, &
        0.5760E-04_dp, -0.3465E-04_dp, -0.8535E-04_dp, -0.2011E-04_dp,  0.6612E-06_dp/)

  ! urban aerosol cosine coefficients
  REAL(dp), PARAMETER, DIMENSION(66) :: caeuc= &
    (/  0.8005E-01_dp,  0.7095E-01_dp,  0.2014E-01_dp, -0.1412E-01_dp, -0.2425E-01_dp, &
       -0.1332E-01_dp, -0.2904E-02_dp,  0.5068E-03_dp,  0.9369E-03_dp,  0.4114E-03_dp, &
        0.7549E-04_dp,  0.1922E-01_dp,  0.2534E-01_dp,  0.2088E-01_dp,  0.1064E-01_dp, &
        0.1063E-02_dp, -0.2526E-02_dp, -0.2091E-02_dp, -0.9660E-03_dp, -0.2030E-03_dp, &
        0.3865E-04_dp, -0.9900E-02_dp, -0.5964E-02_dp,  0.2223E-02_dp,  0.4941E-02_dp, &
        0.3277E-02_dp,  0.1038E-02_dp, -0.1480E-03_dp, -0.2844E-03_dp, -0.1208E-03_dp, &
        0.3999E-02_dp,  0.6282E-02_dp,  0.2813E-02_dp,  0.1475E-02_dp,  0.4571E-03_dp, &
       -0.1349E-03_dp, -0.9011E-04_dp, -0.1936E-04_dp,  0.1994E-02_dp,  0.3540E-02_dp, &
        0.8837E-03_dp,  0.1992E-03_dp,  0.3092E-04_dp, -0.7979E-04_dp, -0.2664E-04_dp, &
       -0.5006E-04_dp,  0.6447E-03_dp,  0.5550E-03_dp,  0.1197E-03_dp,  0.6657E-04_dp, &
        0.1488E-04_dp, -0.9141E-04_dp, -0.2896E-03_dp, -0.1561E-03_dp, -0.6524E-04_dp, &
       -0.1559E-04_dp, -0.1082E-03_dp, -0.4126E-03_dp, -0.1732E-03_dp, -0.8286E-04_dp, &
       -0.1993E-04_dp,  0.3850E-04_dp,  0.2870E-04_dp,  0.4493E-04_dp,  0.4721E-04_dp, &
        0.1338E-04_dp/)

  ! urban aerosol sine coefficients
  REAL(dp), PARAMETER, DIMENSION(55) :: caeus= &
    (/  0.6646E-02_dp,  0.8373E-02_dp,  0.5463E-02_dp,  0.4554E-02_dp,  0.3301E-02_dp, &
        0.5725E-03_dp, -0.7482E-03_dp, -0.6222E-03_dp, -0.2603E-03_dp, -0.5127E-04_dp, &
       -0.3849E-04_dp,  0.9741E-02_dp,  0.8190E-02_dp,  0.5712E-02_dp,  0.3039E-02_dp, &
        0.5290E-03_dp, -0.2044E-03_dp, -0.2309E-03_dp, -0.1160E-03_dp,  0.9160E-02_dp, &
        0.1286E-01_dp,  0.1170E-01_dp,  0.5491E-02_dp,  0.1393E-02_dp, -0.6288E-04_dp, &
       -0.2715E-03_dp, -0.1047E-03_dp,  0.4873E-02_dp,  0.3545E-02_dp,  0.3069E-02_dp, &
        0.1819E-02_dp,  0.6947E-03_dp,  0.1416E-03_dp, -0.1538E-04_dp, -0.4351E-03_dp, &
       -0.1907E-02_dp, -0.5774E-03_dp, -0.2247E-03_dp,  0.5345E-04_dp,  0.9052E-04_dp, &
       -0.3972E-04_dp, -0.9665E-04_dp,  0.7912E-04_dp, -0.1094E-04_dp, -0.6776E-05_dp, &
        0.2724E-03_dp,  0.1973E-03_dp,  0.6837E-04_dp,  0.4313E-04_dp, -0.7174E-05_dp, &
        0.8527E-05_dp, -0.2160E-05_dp, -0.7852E-04_dp,  0.3453E-06_dp, -0.2402E-05_dp/)

  ! desert aerosol cosine coefficients
  REAL(dp), PARAMETER, DIMENSION(66) :: caedc= &
    (/  0.2840E-01_dp,  0.1775E-01_dp, -0.1069E-01_dp, -0.1553E-01_dp, -0.3299E-02_dp, &
        0.3583E-02_dp,  0.2274E-02_dp,  0.5767E-04_dp, -0.3678E-03_dp, -0.1050E-03_dp, &
        0.2133E-04_dp,  0.2326E-01_dp,  0.1566E-01_dp, -0.3130E-02_dp, -0.8253E-02_dp, &
       -0.2615E-02_dp,  0.1247E-02_dp,  0.1059E-02_dp,  0.1196E-03_dp, -0.1303E-03_dp, &
       -0.5094E-04_dp,  0.1185E-01_dp,  0.7238E-02_dp, -0.1562E-02_dp, -0.3665E-02_dp, &
       -0.1182E-02_dp,  0.4678E-03_dp,  0.4448E-03_dp,  0.8307E-04_dp, -0.3468E-04_dp, &
        0.5273E-02_dp,  0.3037E-02_dp, -0.4014E-03_dp, -0.1202E-02_dp, -0.4647E-03_dp, &
        0.5148E-04_dp,  0.1014E-03_dp,  0.2996E-04_dp,  0.2505E-02_dp,  0.1495E-02_dp, &
        0.2438E-03_dp, -0.1223E-03_dp, -0.7669E-04_dp, -0.1638E-04_dp,  0.1869E-05_dp, &
        0.1094E-02_dp,  0.6131E-03_dp,  0.1508E-03_dp,  0.1765E-04_dp,  0.1360E-05_dp, &
       -0.7998E-06_dp,  0.4475E-03_dp,  0.2737E-03_dp,  0.6430E-04_dp, -0.6759E-05_dp, &
       -0.6761E-05_dp,  0.1992E-03_dp,  0.1531E-03_dp,  0.4828E-04_dp,  0.5103E-06_dp, &
        0.7454E-04_dp,  0.5917E-04_dp,  0.2152E-04_dp,  0.9300E-05_dp,  0.9790E-05_dp, &
       -0.8853E-05_dp/)

  ! desert aerosol sine coefficients
  REAL(dp), PARAMETER, DIMENSION(55) :: caeds= &
    (/  0.9815E-02_dp,  0.8436E-02_dp,  0.1087E-02_dp, -0.2717E-02_dp, -0.1755E-02_dp, &
       -0.1559E-03_dp,  0.2367E-03_dp,  0.8808E-04_dp,  0.2001E-05_dp, -0.1244E-05_dp, &
        0.1041E-01_dp,  0.8039E-02_dp,  0.1005E-02_dp, -0.1981E-02_dp, -0.1090E-02_dp, &
        0.1595E-05_dp,  0.1787E-03_dp,  0.4644E-04_dp, -0.1052E-04_dp,  0.6593E-02_dp, &
        0.3983E-02_dp, -0.1527E-03_dp, -0.1235E-02_dp, -0.5078E-03_dp,  0.3649E-04_dp, &
        0.1005E-03_dp,  0.3182E-04_dp,  0.3225E-02_dp,  0.1672E-02_dp, -0.7752E-04_dp, &
       -0.4312E-03_dp, -0.1872E-03_dp, -0.1666E-04_dp,  0.1872E-04_dp,  0.1133E-02_dp, &
        0.5643E-03_dp,  0.7747E-04_dp, -0.2980E-04_dp, -0.2092E-04_dp, -0.8590E-05_dp, &
        0.2988E-03_dp,  0.6714E-04_dp, -0.6249E-05_dp,  0.1052E-04_dp,  0.8790E-05_dp, &
        0.1569E-03_dp, -0.1175E-04_dp, -0.3033E-04_dp, -0.9777E-06_dp,  0.1101E-03_dp, &
        0.6827E-05_dp, -0.1023E-04_dp,  0.4231E-04_dp,  0.4905E-05_dp,  0.6229E-05_dp/)


!=======================================================================
!
! C)

  ! Vertical distribution of sea, land, urban and desert aerosols, 
  ! and background aerosols of tropospheric, stratospheric and 
  ! volcanic type.
  !
  ! replaces ECHAM4 module mo_rad2

  REAL(dp)             :: caeops    ! for optical properties of sea aerosols
  REAL(dp)             :: caeopl    ! for optical properties of land aerosols
  REAL(dp)             :: caeopu    ! for optical properties of urban aerosols
  REAL(dp)             :: caeopd    ! for optical properties of desert aerosols
  REAL(dp),ALLOCATABLE :: cvdaes(:) ! for vertical distribution of sea aerosols
  REAL(dp),ALLOCATABLE :: cvdael(:) ! for vertical distribution of land aerosols
  REAL(dp),ALLOCATABLE :: cvdaeu(:) ! for vertical distribution of urban aerosols
  REAL(dp),ALLOCATABLE :: cvdaed(:) ! for vertical distribution of desert aerosols
  REAL(dp)             :: ctrbga    ! for tropospheric background aerosol
  REAL(dp)             :: cvobga    ! for volcanic background aerosol
  REAL(dp)             :: cstbga    ! for stratospheric background aerosol
  REAL(dp)             :: ctrpt     ! temperature exponent for strat. aerosol
  REAL(dp)             :: caeadk(3) ! for moisture adsorption of aerosols
  REAL(dp)             :: caeadm    ! for moisture adsorption of aerosols


!=======================================================================
!
! D)

  ! Optical properties of sea, land, urban and desert aerosols in
  ! 4 SW  bands and 5 LW bands (as in ECMWF)
  !
  ! replaces ECHAM4 module mo_aerosols

  !  shortwave
  !  index 1: SW band
  !  index 2: aerosol mixture type

  ! normalized optical thickness at 0.55 micron
  REAL(dp),DIMENSION(6,naer) :: taua=RESHAPE((/&
         0.730719_dp, 0.730719_dp, 0.730719_dp, 0.730719_dp, 0.730719_dp, 0.730719_dp,&
         0.912819_dp, 0.912819_dp, 0.912819_dp, 0.912819_dp, 0.912819_dp, 0.912819_dp,&
         0.725059_dp, 0.725059_dp, 0.725059_dp, 0.725059_dp, 0.725059_dp, 0.725059_dp,&
         0.745405_dp, 0.745405_dp, 0.745405_dp, 0.745405_dp, 0.745405_dp, 0.745405_dp,&
         0.682188_dp, 0.682188_dp, 0.682188_dp, 0.682188_dp, 0.682188_dp, 0.682188_dp /)&
        ,SHAPE=(/6,naer/))

  ! single scattering albedo
  REAL(dp),DIMENSION(6,naer) :: piza=RESHAPE((/&
         0.872212_dp, 0.872212_dp, 0.872212_dp, 0.872212_dp, 0.872212_dp, 0.872212_dp,&
         0.982545_dp, 0.982545_dp, 0.982545_dp, 0.982545_dp, 0.982545_dp, 0.982545_dp,&
         0.623143_dp, 0.623143_dp, 0.623143_dp, 0.623143_dp, 0.623143_dp, 0.623143_dp,&
         0.944887_dp, 0.944887_dp, 0.944887_dp, 0.944887_dp, 0.944887_dp, 0.944887_dp,&
         0.997975_dp, 0.997975_dp, 0.997975_dp, 0.997975_dp, 0.997975_dp, 0.997975_dp /)&
        ,SHAPE=(/6,naer/))

  ! asymmetry factor
  REAL(dp),DIMENSION(6,naer) :: cga=RESHAPE((/&
         0.647596_dp, 0.647596_dp, 0.647596_dp, 0.647596_dp, 0.647596_dp, 0.647596_dp,&
         0.739002_dp, 0.739002_dp, 0.739002_dp, 0.739002_dp, 0.739002_dp, 0.739002_dp,&
         0.580845_dp, 0.580845_dp, 0.580845_dp, 0.580845_dp, 0.580845_dp, 0.580845_dp,&
         0.662657_dp, 0.662657_dp, 0.662657_dp, 0.662657_dp, 0.662657_dp, 0.662657_dp,&
         0.624246_dp, 0.624246_dp, 0.624246_dp, 0.624246_dp, 0.624246_dp, 0.624246_dp /)&
        ,SHAPE=(/6,naer/))

  !  longwave
  !  index 1: LW band
  !  index 2: aerosol mixture type

  ! absorption coefficients
  REAL(dp),DIMENSION(5,naer) :: caer=RESHAPE((/&
         0.038520_dp, 0.037196_dp, 0.040532_dp, 0.054934_dp, 0.038520_dp ,&
         0.12613_dp , 0.18313_dp , 0.10357_dp , 0.064106_dp, 0.126130_dp ,&
         0.012579_dp, 0.013649_dp, 0.018652_dp, 0.025181_dp, 0.012579_dp ,&
         0.011890_dp, 0.016142_dp, 0.021105_dp, 0.028908_dp, 0.011890_dp ,&
         0.013792_dp, 0.026810_dp, 0.052203_dp, 0.066338_dp, 0.013792_dp /)&
        ,SHAPE=(/5,naer/))

!=======================================================================

  CONTAINS

!=======================================================================
!
! E)

  SUBROUTINE su_aero_tanre(petah)

  ! Description:
  !
  ! parameters for the vertical distributions of aerosols.
  !
  ! Method:
  !
  ! This routine computes the values *cvdaen* (*n=*s,*l,*u or *d
  ! for sea,land,urban or desert) of a surface-normalised vertical
  ! distribution of aerosols' optical dephts from the argument *petah*
  ! (vertical coordinate) at *klevp1* levels. It also sets values for
  ! non-geographically weighted total optical depths (at 0.55 e-06
  ! wave-length) *caeopn* for the same four types and similear optical
  ! dephts divided by pressure for background well-mixed aerosols
  ! of three types *cmnbga* (*mn*=*tr*,*vo* or *st* for tropospheric,
  ! volcanic (stratospheric ashes) or stratospheric (sulfuric type)).
  ! It finally set values for the power to be applied to a temperature
  ! ratio smaller than one in order to obtain an index one in the
  ! stratosphere and zero in the troposphere with a relatively smooth
  ! transition (*ctrpt*), as well as for adsorption coefficients for
  ! water to the three type of tropospheric aerosols (*caeadk*) with
  ! a minimum value (in the whole atmosphere) for the sum of the
  ! products of *caeadk* by the optical depths divided by pressure
  ! thickness: *caeadm*.
  !
  ! *su_aero_tanre* is called from *physc*.
  ! there is one dummy argument:
  ! 
  !   *petah*  is the vertical coordinate array of size nlev+1.
  !
  ! *su_aero_tanre initializes :
  !   *cvdaen* (*n=*s,*l,*u or*d) are the normalised vertical distributions.
  !   *cmnbga* (*mn*=*tr*,*vo* or *st*) 
  !            are the background optical depths divided by pressure.
  !   *caeopn* (*n=*s,*l,*u or *d) are the total optical dephts 
  !            for the vertically varying aerosols.
  !   *ctrpt*  is the temperature exponent for the stratospheric definition.
  !   *caeadk* (1,2,3) and
  !   *caeadm* are the constants for the definition of the quantity
  !   of water vapour that will be adsorbed to the dry aerosols to form
  !   moist aerosols.
  !
  ! straightforward, equivalent heigths are given in meters (8434
  ! for the atmosphere) and tropospheric and stratospheric pressure
  ! boundary values are set at 101325 and 19330 *pascal.
  !
  ! Authors:
  !
  ! J. F. Geleyn, ECMWF, November 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! M. A. Giorgetta, MPI, April 1999, modified for new module mo_aero_tanre
  ! 
  ! for more details see file AUTHORS
  !

  !  Array argument
  REAL(dp), INTENT(in) :: petah(:)

  !  Module scalar variables
  !  REAL :: caeopd, caeopl, caeops, caeopu, &
  !          cstbga, ctrbga, cvobga, caeadm, ctrpt
  !
  !  Module array variables
  !  REAL, DIMENSION(nlevp1) :: cvdaed, cvdael, cvdaes, cvdaeu
  !  REAL, DIMENSION(3)      :: caeadk(3)

  !  Local: 
  INTEGER :: klevp1, jk
  REAL(dp) :: zhsd, zhsl, zhss, zhsu

  !  Intrinsic functions 
  INTRINSIC MAX,SIZE

  !  Executable statements 

  !  get vertical extension
  klevp1 = SIZE(petah)

  !  and allocate vertical arrays
  IF (.NOT. ALLOCATED(cvdaes)) ALLOCATE(cvdaes(klevp1))
  IF (.NOT. ALLOCATED(cvdael)) ALLOCATE(cvdael(klevp1))
  IF (.NOT. ALLOCATED(cvdaeu)) ALLOCATE(cvdaeu(klevp1))
  IF (.NOT. ALLOCATED(cvdaed)) ALLOCATE(cvdaed(klevp1))

  !  Computations

  zhss = MAX(1._dp,8434._dp/1000._dp)
  zhsl = MAX(1._dp,8434._dp/1000._dp)
  zhsu = MAX(1._dp,8434._dp/1000._dp)
  zhsd = MAX(1._dp,8434._dp/3000._dp)
  cvdaes(1) = 0._dp
  cvdael(1) = 0._dp
  cvdaeu(1) = 0._dp
  cvdaed(1) = 0._dp
  IF (petah(1)/=0._dp) THEN
    cvdaes(1) = petah(1)**zhss
    cvdael(1) = petah(1)**zhsl
    cvdaeu(1) = petah(1)**zhsu
    cvdaed(1) = petah(1)**zhsd
  END IF
  DO jk = 2, klevp1
    cvdaes(jk) = petah(jk)**zhss
    cvdael(jk) = petah(jk)**zhsl
    cvdaeu(jk) = petah(jk)**zhsu
    cvdaed(jk) = petah(jk)**zhsd
  END DO
  ctrbga = 0.03_dp/(101325._dp-19330._dp)
  cvobga = 0.007_dp/19330._dp
!  cvobga =1.E-45_dp
  cstbga = 0.045_dp/19330._dp
!  cstbga =1.E-45_dp
  caeops = 0.05_dp
  caeopl = 0.2_dp
  caeopu = 0.1_dp
  caeopd = 1.9_dp
  ctrpt  = 30._dp
  caeadk(1) =  0.3876E-03_dp
  caeadk(2) =  0.6693E-02_dp
  caeadk(3) =  0.8563E-03_dp

  caeadm = 2.6E-10_dp

END SUBROUTINE su_aero_tanre

!=======================================================================
!
! F)

FUNCTION aero_tanre(krow,kproma,pdp,ppf,pth)

  ! follows ECHAM4 radint

  USE mo_geoloc,         ONLY: zaes_x, zael_x, zaeu_x, zaed_x

  INTEGER                             :: krow,kproma,klev

  REAL(dp), INTENT(in) , DIMENSION(:,:)   :: pdp,ppf,pth
  REAL(dp), DIMENSION(kproma,SIZE(pdp,2),naer) :: aero_tanre

  REAL(dp), DIMENSION(kproma,SIZE(pdp,2),naer) :: zaero

  REAL(dp), DIMENSION(kproma)        :: zaetrn, zaetro

  REAL(dp)    :: zaeqsn, zaeqln, zaequn, zaeqdn, &
                 zaeqso, zaeqlo, zaequo, zaeqdo, &
                 zaetr

  INTEGER :: jl,jk

  INTRINSIC MAX, MIN, SQRT


  klev = SIZE(pdp,2)

  ! output: aerosol longitude height sections in aero_tanre

  DO jl = 1, kproma
    zaetro(jl) = 1._dp
  END DO
  DO jk = 1, klev
    DO jl = 1, kproma
      zaeqso = caeops*zaes_x(jl,krow)*cvdaes(jk)
      zaeqsn = caeops*zaes_x(jl,krow)*cvdaes(jk+1)
      zaeqlo = caeopl*zael_x(jl,krow)*cvdael(jk)
      zaeqln = caeopl*zael_x(jl,krow)*cvdael(jk+1)
      zaequo = caeopu*zaeu_x(jl,krow)*cvdaeu(jk)
      zaequn = caeopu*zaeu_x(jl,krow)*cvdaeu(jk+1)
      zaeqdo = caeopd*zaed_x(jl,krow)*cvdaed(jk)
      zaeqdn = caeopd*zaed_x(jl,krow)*cvdaed(jk+1)
      IF (ppf(jl,jk)<999._dp) THEN
         zaetr= 1._dp ! above 10 hPa 
      ELSE
         zaetrn(jl) = zaetro(jl)*(MIN(1.0_dp,pth(jl,jk)/pth(jl,jk+1)))**ctrpt
         zaetr= SQRT(zaetrn(jl)*zaetro(jl))
         zaetro(jl) = zaetrn(jl)
      END IF
      zaero(jl,jk,1) = (1._dp-zaetr)*(ctrbga*pdp(jl,jk)+zaeqln-zaeqlo+zaeqdn-zaeqdo)
      zaero(jl,jk,2) = (1._dp-zaetr)*(zaeqsn-zaeqso)
      zaero(jl,jk,3) = (1._dp-zaetr)*(zaequn-zaequo)
      zaero(jl,jk,4) = zaetr*cvobga*pdp(jl,jk)
      zaero(jl,jk,5) = zaetr*cstbga*pdp(jl,jk)
    END DO
  END DO

  ! optical thickness is not negative
  aero_tanre(:,:,:)=MAX(zaero(:,:,1:naer),0.0_dp) 

END FUNCTION aero_tanre

SUBROUTINE init_aero

  USE mo_gaussgrid,      ONLY: coslon, sinlon, gl_twomu
  USE mo_geoloc,         ONLY: zaes_x, zael_x, zaeu_x, zaed_x
  USE mo_decomposition,  ONLY: ldc => local_decomposition
  USE mo_transpose,      ONLY: reorder

  REAL(dp), DIMENSION(21) :: zfaes,zfael,zfaeu,zfaed
  REAL(dp), DIMENSION(66) :: zalp

  REAL(dp)    :: zsin
  REAL(dp)    :: zcos1, zcos2, zcos3, zcos4, zcos5,  &
                 zcos6, zcos7, zcos8, zcos9, zcos10, &
                 zsin1, zsin2, zsin3, zsin4, zsin5,  &
                 zsin6, zsin7, zsin8, zsin9, zsin10
  REAL(dp)    :: zaes(ldc% nglon, ldc% nglat)
  REAL(dp)    :: zael(ldc% nglon, ldc% nglat)
  REAL(dp)    :: zaeu(ldc% nglon, ldc% nglat)
  REAL(dp)    :: zaed(ldc% nglon, ldc% nglat)

  INTEGER :: jmm,imm,imnc,imns,jnn,jl

  !  Local array bounds
  INTEGER :: nglon, nglat, nlof
  INTEGER :: jlat, jglat

  !  External subroutines
  EXTERNAL  legtri  ! for the Legendre transform at a given latitude


  !  Executable statements

  !  Local array bounds

  nglon = ldc% nglon      ! local number of longitudes
  nglat = ldc% nglat      ! local number of latitudes

  DO jlat = 1, nglat

    jglat = ldc% glat(jlat)
    nlof  = ldc% glon(jlat) ! longitude offset to global field
     
    zsin = .5_dp*gl_twomu(jglat)

    ! Prepare Legendre transform at given latitute

    CALL legtri(zsin,11,zalp)

    ! Legendre transform

    zfaes(:) = 0._dp
    zfael(:) = 0._dp
    zfaeu(:) = 0._dp
    zfaed(:) = 0._dp
    imm  = 0
    imnc = 0
    imns = 0
    DO jmm = 1, 11
      imm = imm + 1
      DO jnn = jmm, 11
        imnc = imnc + 1
        zfaes(imm) = zfaes(imm) + zalp(imnc)*caesc(imnc)
        zfael(imm) = zfael(imm) + zalp(imnc)*caelc(imnc)
        zfaeu(imm) = zfaeu(imm) + zalp(imnc)*caeuc(imnc)
        zfaed(imm) = zfaed(imm) + zalp(imnc)*caedc(imnc)
      END DO
      IF (jmm/=1) THEN
        imm = imm + 1
        DO jnn = jmm, 11
          imns = imns + 1
          zfaes(imm) = zfaes(imm) + zalp(imns+11)*caess(imns)
          zfael(imm) = zfael(imm) + zalp(imns+11)*caels(imns)
          zfaeu(imm) = zfaeu(imm) + zalp(imns+11)*caeus(imns)
          zfaed(imm) = zfaed(imm) + zalp(imns+11)*caeds(imns)
        END DO
      END IF
    END DO

    ! Fourier transform

    DO jl = 1, nglon

      zcos1 = coslon(jl+nlof)
      zsin1 = sinlon(jl+nlof)
      zcos2 = zcos1*zcos1 - zsin1*zsin1
      zsin2 = zsin1*zcos1 + zcos1*zsin1
      zcos3 = zcos2*zcos1 - zsin2*zsin1
      zsin3 = zsin2*zcos1 + zcos2*zsin1
      zcos4 = zcos3*zcos1 - zsin3*zsin1
      zsin4 = zsin3*zcos1 + zcos3*zsin1
      zcos5 = zcos4*zcos1 - zsin4*zsin1
      zsin5 = zsin4*zcos1 + zcos4*zsin1
      zcos6 = zcos5*zcos1 - zsin5*zsin1
      zsin6 = zsin5*zcos1 + zcos5*zsin1
      zcos7 = zcos6*zcos1 - zsin6*zsin1
      zsin7 = zsin6*zcos1 + zcos6*zsin1
      zcos8 = zcos7*zcos1 - zsin7*zsin1
      zsin8 = zsin7*zcos1 + zcos7*zsin1
      zcos9 = zcos8*zcos1 - zsin8*zsin1
      zsin9 = zsin8*zcos1 + zcos8*zsin1
      zcos10= zcos9*zcos1 - zsin9*zsin1
      zsin10= zsin9*zcos1 + zcos9*zsin1
    
      zaes(jl,jlat) = zfaes(1) + 2._dp*(zfaes(2)*zcos1+zfaes(3)*zsin1+zfaes(4)*zcos2+   &
                      zfaes(5)*zsin2+zfaes(6)*zcos3+zfaes(7)*zsin3+zfaes(8)*zcos4+      &
                      zfaes(9)*zsin4+zfaes(10)*zcos5+zfaes(11)*zsin5+zfaes(12)*zcos6+   &
                      zfaes(13)*zsin6+zfaes(14)*zcos7+zfaes(15)*zsin7+zfaes(16)*zcos8+  &
                      zfaes(17)*zsin8+zfaes(18)*zcos9+zfaes(19)*zsin9+zfaes(20)*zcos10+ &
                      zfaes(21)*zsin10)
      zael(jl,jlat) = zfael(1) + 2._dp*(zfael(2)*zcos1+zfael(3)*zsin1+zfael(4)*zcos2+   &
                      zfael(5)*zsin2+zfael(6)*zcos3+zfael(7)*zsin3+zfael(8)*zcos4+      &
                      zfael(9)*zsin4+zfael(10)*zcos5+zfael(11)*zsin5+zfael(12)*zcos6+   &
                      zfael(13)*zsin6+zfael(14)*zcos7+zfael(15)*zsin7+zfael(16)*zcos8+  &
                      zfael(17)*zsin8+zfael(18)*zcos9+zfael(19)*zsin9+zfael(20)*zcos10+ &
                      zfael(21)*zsin10)
      zaeu(jl,jlat) = zfaeu(1) + 2._dp*(zfaeu(2)*zcos1+zfaeu(3)*zsin1+zfaeu(4)*zcos2+   &
                      zfaeu(5)*zsin2+zfaeu(6)*zcos3+zfaeu(7)*zsin3+zfaeu(8)*zcos4+      &
                      zfaeu(9)*zsin4+zfaeu(10)*zcos5+zfaeu(11)*zsin5+zfaeu(12)*zcos6+   &
                      zfaeu(13)*zsin6+zfaeu(14)*zcos7+zfaeu(15)*zsin7+zfaeu(16)*zcos8+  &
                      zfaeu(17)*zsin8+zfaeu(18)*zcos9+zfaeu(19)*zsin9+zfaeu(20)*zcos10+ &
                      zfaeu(21)*zsin10)
      zaed(jl,jlat) = zfaed(1) + 2._dp*(zfaed(2)*zcos1+zfaed(3)*zsin1+zfaed(4)*zcos2+   &
                      zfaed(5)*zsin2+zfaed(6)*zcos3+zfaed(7)*zsin3+zfaed(8)*zcos4+      &
                      zfaed(9)*zsin4+zfaed(10)*zcos5+zfaed(11)*zsin5+zfaed(12)*zcos6+   &
                      zfaed(13)*zsin6+zfaed(14)*zcos7+zfaed(15)*zsin7+zfaed(16)*zcos8+  &
                      zfaed(17)*zsin8+zfaed(18)*zcos9+zfaed(19)*zsin9+zfaed(20)*zcos10+ &
                      zfaed(21)*zsin10)

    END DO

  END DO

  CALL reorder (zaes_x, zaes)
  CALL reorder (zael_x, zael)
  CALL reorder (zaeu_x, zaeu)
  CALL reorder (zaed_x, zaed)

END SUBROUTINE init_aero
!------------------------------------------------------------------------------
  SUBROUTINE cleanup_aero
    !
    ! deallocate module variables
    !
    IF (ALLOCATED(cvdaes)) DEALLOCATE(cvdaes)
    IF (ALLOCATED(cvdael)) DEALLOCATE(cvdael)
    IF (ALLOCATED(cvdaeu)) DEALLOCATE(cvdaeu)
    IF (ALLOCATED(cvdaed)) DEALLOCATE(cvdaed)
  END SUBROUTINE cleanup_aero
!=======================================================================

END MODULE mo_aero_tanre
