MODULE mo_vsop87

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  REAL(dp), PARAMETER :: pi   = 3.14159265358979323846_dp
  REAL(dp), PARAMETER :: dtor = pi/180_dp
  REAL(dp), PARAMETER :: stor = dtor/3600_dp

  TYPE terms 
    REAL(dp) :: A
    REAL(dp) :: B
    REAL(dp) :: C
  END TYPE terms

  ! Data given are based on VSOP87D

  TYPE (terms), PARAMETER :: L1(64) = (/ &
       terms(1.753470456730000e+00_dp, 0.00000000000e+00_dp, 0.000000000000000e+00_dp), &
       terms(3.341656456000000e-02_dp, 4.66925680417e+00_dp, 6.283075849991400e+03_dp), &
       terms(3.489427500000000e-04_dp, 4.62610241759e+00_dp, 1.256615169998280e+04_dp), &
       terms(3.497056000000000e-05_dp, 2.74411800971e+00_dp, 5.753384884896800e+03_dp), &
       terms(3.417571000000000e-05_dp, 2.82886579606e+00_dp, 3.523118349000000e+00_dp), &
       terms(3.135896000000000e-05_dp, 3.62767041758e+00_dp, 7.771377146812050e+04_dp), &
       terms(2.676218000000000e-05_dp, 4.41808351397e+00_dp, 7.860419392439200e+03_dp), &
       terms(2.342687000000000e-05_dp, 6.13516237631e+00_dp, 3.930209696219600e+03_dp), &
       terms(1.324292000000000e-05_dp, 7.42463563520e-01_dp, 1.150676976979360e+04_dp), &
       terms(1.273166000000000e-05_dp, 2.03709655772e+00_dp, 5.296909650946000e+02_dp), &
       terms(1.199167000000000e-05_dp, 1.10962944315e+00_dp, 1.577343542447800e+03_dp), &
       terms(9.902500000000000e-06_dp, 5.23268129594e+00_dp, 5.884926846583200e+03_dp), &
       terms(9.018550000000000e-06_dp, 2.04505443513e+00_dp, 2.629831979980000e+01_dp), &
       terms(8.572229999999999e-06_dp, 3.50849156957e+00_dp, 3.981490034082000e+02_dp), &
       terms(7.797859999999999e-06_dp, 1.17882652114e+00_dp, 5.223693919802200e+03_dp), &
       terms(7.531410000000000e-06_dp, 2.53339053818e+00_dp, 5.507553238667400e+03_dp), &
       terms(5.052640000000000e-06_dp, 4.58292563052e+00_dp, 1.884922754997420e+04_dp), &
       terms(4.923790000000000e-06_dp, 4.20506639861e+00_dp, 7.755226113240000e+02_dp), &
       terms(3.566550000000000e-06_dp, 2.91954116867e+00_dp, 6.731030280000000e-02_dp), &
       terms(3.170870000000000e-06_dp, 5.84901952218e+00_dp, 1.179062908865880e+04_dp), &
       terms(2.841250000000000e-06_dp, 1.89869034186e+00_dp, 7.962980068163999e+02_dp), &
       terms(2.710390000000000e-06_dp, 3.14886076490e-01_dp, 1.097707880469900e+04_dp), &
       terms(2.428100000000000e-06_dp, 3.44811409060e-01_dp, 5.486777843175000e+03_dp), &
       terms(2.061600000000000e-06_dp, 4.80646606059e+00_dp, 2.544314419883400e+03_dp), &
       terms(2.053850000000000e-06_dp, 1.86947813692e+00_dp, 5.573142801433100e+03_dp), &
       terms(2.022610000000000e-06_dp, 2.45767795458e+00_dp, 6.069776754553400e+03_dp), &
       terms(1.555160000000000e-06_dp, 8.33060738070e-01_dp, 2.132990954380000e+02_dp), &
       terms(1.322120000000000e-06_dp, 3.41118275555e+00_dp, 2.942463423291600e+03_dp), &
       terms(1.261840000000000e-06_dp, 1.08302630210e+00_dp, 2.077539549240000e+01_dp), &
       terms(1.151320000000000e-06_dp, 6.45449116830e-01_dp, 9.803210682000000e-01_dp), &
       terms(1.028510000000000e-06_dp, 6.35998467270e-01_dp, 4.694002954707600e+03_dp), &
       terms(1.018950000000000e-06_dp, 9.75692218240e-01_dp, 1.572083878487840e+04_dp), &
       terms(1.017240000000000e-06_dp, 4.26679821365e+00_dp, 7.113547000800000e+00_dp), &
       terms(9.920600000000000e-07_dp, 6.20992940258e+00_dp, 2.146165416475200e+03_dp), &
       terms(9.760700000000001e-07_dp, 6.81012722700e-01_dp, 1.554203994342000e+02_dp), &
       terms(8.580300000000000e-07_dp, 5.98322631256e+00_dp, 1.610006857376741e+05_dp), &
       terms(8.512800000000000e-07_dp, 1.29870743025e+00_dp, 6.275962302990600e+03_dp), &
       terms(8.471100000000000e-07_dp, 3.67080093025e+00_dp, 7.143069561812909e+04_dp), &
       terms(7.963700000000000e-07_dp, 1.80791330700e+00_dp, 1.726015465469040e+04_dp), &
       terms(7.875600000000000e-07_dp, 3.03698313141e+00_dp, 1.203646073488820e+04_dp), &
       terms(7.465100000000000e-07_dp, 1.75508916159e+00_dp, 5.088628839766800e+03_dp), &
       terms(7.387400000000000e-07_dp, 3.50319443167e+00_dp, 3.154687084895600e+03_dp), &
       terms(7.354700000000000e-07_dp, 4.67926565481e+00_dp, 8.018209311238001e+02_dp), &
       terms(6.962700000000000e-07_dp, 8.32975969660e-01_dp, 9.437762934887000e+03_dp), &
       terms(6.244899999999999e-07_dp, 3.97763880587e+00_dp, 8.827390269874801e+03_dp), &
       terms(6.114800000000000e-07_dp, 1.81839811024e+00_dp, 7.084896781115200e+03_dp), &
       terms(5.696300000000000e-07_dp, 2.78430398043e+00_dp, 6.286598968340400e+03_dp), &
       terms(5.611600000000000e-07_dp, 4.38694880779e+00_dp, 1.414349524243060e+04_dp), &
       terms(5.557700000000000e-07_dp, 3.47006009062e+00_dp, 6.279552731642400e+03_dp), &
       terms(5.199200000000000e-07_dp, 1.89149458340e-01_dp, 1.213955350910680e+04_dp), &
       terms(5.160500000000000e-07_dp, 1.33282746983e+00_dp, 1.748016413067000e+03_dp), &
       terms(5.114500000000000e-07_dp, 2.83068645010e-01_dp, 5.856477659115400e+03_dp), &
       terms(4.900000000000000e-07_dp, 4.87350650330e-01_dp, 1.194447010224600e+03_dp), &
       terms(4.103600000000000e-07_dp, 5.36817351402e+00_dp, 8.429241266466601e+03_dp), &
       terms(4.093800000000000e-07_dp, 2.39850881707e+00_dp, 1.965104848109800e+04_dp), &
       terms(3.920000000000000e-07_dp, 6.16832995016e+00_dp, 1.044738783960440e+04_dp), &
       terms(3.677000000000000e-07_dp, 6.04133859347e+00_dp, 1.021328554621100e+04_dp), &
       terms(3.659600000000000e-07_dp, 2.56955238628e+00_dp, 1.059381930189200e+03_dp), &
       terms(3.595400000000000e-07_dp, 1.70876111898e+00_dp, 2.352866153771800e+03_dp), &
       terms(3.556600000000000e-07_dp, 1.77597314691e+00_dp, 6.812766815086000e+03_dp), &
       terms(3.329100000000000e-07_dp, 5.93094994590e-01_dp, 1.778984561978500e+04_dp), &
       terms(3.041200000000000e-07_dp, 4.42944641350e-01_dp, 8.399684731811189e+04_dp), &
       terms(3.004700000000000e-07_dp, 2.73975123935e+00_dp, 1.349867409658800e+03_dp), &
       terms(2.535200000000000e-07_dp, 3.16470953405e+00_dp, 4.690479836358600e+03_dp) /)

  TYPE (terms), PARAMETER :: L2(34) = (/ &
       terms(6.283319667474910e+03_dp, 0.00000000000e+00_dp, 0.000000000000000e+00_dp), &
       terms(2.060588630000000e-03_dp, 2.67823455584e+00_dp, 6.283075849991400e+03_dp), &
       terms(4.303430000000000e-05_dp, 2.63512650414e+00_dp, 1.256615169998280e+04_dp), &
       terms(4.252640000000000e-06_dp, 1.59046980729e+00_dp, 3.523118349000000e+00_dp), &
       terms(1.192610000000000e-06_dp, 5.79557487799e+00_dp, 2.629831979980000e+01_dp), &
       terms(1.089770000000000e-06_dp, 2.96618001993e+00_dp, 1.577343542447800e+03_dp), &
       terms(9.347800000000000e-07_dp, 2.59212835365e+00_dp, 1.884922754997420e+04_dp), &
       terms(7.212200000000000e-07_dp, 1.13846158196e+00_dp, 5.296909650946000e+02_dp), &
       terms(6.776800000000000e-07_dp, 1.87472304791e+00_dp, 3.981490034082000e+02_dp), &
       terms(6.732700000000000e-07_dp, 4.40918235168e+00_dp, 5.507553238667400e+03_dp), &
       terms(5.902700000000000e-07_dp, 2.88797038460e+00_dp, 5.223693919802200e+03_dp), &
       terms(5.597600000000000e-07_dp, 2.17471680261e+00_dp, 1.554203994342000e+02_dp), &
       terms(4.540700000000000e-07_dp, 3.98030798050e-01_dp, 7.962980068163999e+02_dp), &
       terms(3.636900000000000e-07_dp, 4.66247398350e-01_dp, 7.755226113240000e+02_dp), &
       terms(2.895800000000000e-07_dp, 2.64707383882e+00_dp, 7.113547000800000e+00_dp), &
       terms(2.084400000000000e-07_dp, 5.34138275149e+00_dp, 9.803210682000000e-01_dp), &
       terms(1.909700000000000e-07_dp, 1.84628332577e+00_dp, 5.486777843175000e+03_dp), &
       terms(1.850800000000000e-07_dp, 4.96855124577e+00_dp, 2.132990954380000e+02_dp), &
       terms(1.729300000000000e-07_dp, 2.99116864949e+00_dp, 6.275962302990600e+03_dp), &
       terms(1.623300000000000e-07_dp, 3.21648304700e-02_dp, 2.544314419883400e+03_dp), &
       terms(1.583200000000000e-07_dp, 1.43049285325e+00_dp, 2.146165416475200e+03_dp), &
       terms(1.461500000000000e-07_dp, 1.20532366323e+00_dp, 1.097707880469900e+04_dp), &
       terms(1.246100000000000e-07_dp, 2.83432285512e+00_dp, 1.748016413067000e+03_dp), &
       terms(1.187700000000000e-07_dp, 3.25804815607e+00_dp, 5.088628839766800e+03_dp), &
       terms(1.180800000000000e-07_dp, 5.27379790480e+00_dp, 1.194447010224600e+03_dp), &
       terms(1.151400000000000e-07_dp, 2.07502418155e+00_dp, 4.694002954707600e+03_dp), &
       terms(1.064100000000000e-07_dp, 7.66141992020e-01_dp, 5.535694028424000e+02_dp), &
       terms(9.969000000000000e-08_dp, 1.30262991097e+00_dp, 6.286598968340400e+03_dp), &
       terms(9.720999999999999e-08_dp, 4.23925472239e+00_dp, 1.349867409658800e+03_dp), &
       terms(9.452000000000000e-08_dp, 2.69957062864e+00_dp, 2.427286039740000e+02_dp), &
       terms(8.577000000000001e-08_dp, 5.64475868067e+00_dp, 9.517184062506000e+02_dp), &
       terms(7.576000000000000e-08_dp, 5.30062664886e+00_dp, 2.352866153771800e+03_dp), &
       terms(6.385000000000001e-08_dp, 2.65033984967e+00_dp, 9.437762934887000e+03_dp), &
       terms(6.101000000000000e-08_dp, 4.66632584188e+00_dp, 4.690479836358600e+03_dp) /)

  TYPE (terms), PARAMETER :: L3(20) = (/ &
       terms(5.291887000000000e-04_dp, 0.00000000000e+00_dp, 0.000000000000000e+00_dp), &
       terms(8.719837000000000e-05_dp, 1.07209665242e+00_dp, 6.283075849991400e+03_dp), &
       terms(3.091250000000000e-06_dp, 8.67288188320e-01_dp, 1.256615169998280e+04_dp), &
       terms(2.733900000000000e-07_dp, 5.29787169100e-02_dp, 3.523118349000000e+00_dp), &
       terms(1.633400000000000e-07_dp, 5.18826691036e+00_dp, 2.629831979980000e+01_dp), &
       terms(1.575200000000000e-07_dp, 3.68457889430e+00_dp, 1.554203994342000e+02_dp), &
       terms(9.541000000000001e-08_dp, 7.57422976750e-01_dp, 1.884922754997420e+04_dp), &
       terms(8.937000000000000e-08_dp, 2.05705419118e+00_dp, 7.771377146812050e+04_dp), &
       terms(6.952000000000000e-08_dp, 8.26733054100e-01_dp, 7.755226113240000e+02_dp), &
       terms(5.064000000000000e-08_dp, 4.66284525271e+00_dp, 1.577343542447800e+03_dp), &
       terms(4.061000000000000e-08_dp, 1.03057162962e+00_dp, 7.113547000800000e+00_dp), &
       terms(3.810000000000000e-08_dp, 3.44050803490e+00_dp, 5.573142801433100e+03_dp), &
       terms(3.463000000000000e-08_dp, 5.14074632811e+00_dp, 7.962980068163999e+02_dp), &
       terms(3.169000000000000e-08_dp, 6.05291851171e+00_dp, 5.507553238667400e+03_dp), &
       terms(3.020000000000000e-08_dp, 1.19246506441e+00_dp, 2.427286039740000e+02_dp), &
       terms(2.886000000000000e-08_dp, 6.11652627155e+00_dp, 5.296909650946000e+02_dp), &
       terms(2.714000000000000e-08_dp, 3.06378810250e-01_dp, 3.981490034082000e+02_dp), &
       terms(2.538000000000000e-08_dp, 2.27992810679e+00_dp, 5.535694028424000e+02_dp), &
       terms(2.371000000000000e-08_dp, 4.38118838167e+00_dp, 5.223693919802200e+03_dp), &
       terms(2.079000000000000e-08_dp, 3.75435330484e+00_dp, 9.803210682000000e-01_dp) /)

  TYPE (terms), PARAMETER :: L4(7) = (/ &
       terms(2.892260000000000e-06_dp, 5.84384198723e+00_dp, 6.283075849991400e+03_dp), &
       terms(3.495500000000000e-07_dp, 0.00000000000e+00_dp, 0.000000000000000e+00_dp), &
       terms(1.681900000000000e-07_dp, 5.48766912348e+00_dp, 1.256615169998280e+04_dp), &
       terms(2.962000000000000e-08_dp, 5.19577265202e+00_dp, 1.554203994342000e+02_dp), &
       terms(1.288000000000000e-08_dp, 4.72200252235e+00_dp, 3.523118349000000e+00_dp), &
       terms(7.140000000000000e-09_dp, 5.30045809128e+00_dp, 1.884922754997420e+04_dp), &
       terms(6.350000000000000e-09_dp, 5.96925937141e+00_dp, 2.427286039740000e+02_dp) /)

  TYPE (terms), PARAMETER :: L5(3) = (/ &
       terms(1.140840000000000e-06_dp, 3.14159265359e+00_dp, 0.000000000000000e+00_dp), &
       terms(7.717000000000000e-08_dp, 4.13446589358e+00_dp, 6.283075849991400e+03_dp), &
       terms(7.650000000000001e-09_dp, 3.83803776214e+00_dp, 1.256615169998280e+04_dp) /)

  TYPE (terms), PARAMETER :: L6(1) = (/ &
       terms(8.780000000000000e-09_dp, 3.14159265359e+00_dp, 0.000000000000000e+00_dp) /)

  TYPE (terms), PARAMETER :: B1(5) = (/ &
       terms(2.796200000000000e-06_dp, 3.19870156017e+00_dp, 8.433466158130829e+04_dp), &
       terms(1.016430000000000e-06_dp, 5.42248619256e+00_dp, 5.507553238667400e+03_dp), &
       terms(8.044500000000000e-07_dp, 3.88013204458e+00_dp, 5.223693919802200e+03_dp), &
       terms(4.380600000000000e-07_dp, 3.70444689758e+00_dp, 2.352866153771800e+03_dp), &
       terms(3.193300000000000e-07_dp, 4.00026369781e+00_dp, 1.577343542447800e+03_dp) /)

  TYPE (terms), PARAMETER :: B2(2) = (/ &
       terms(9.029999999999999e-08_dp, 3.89729061890e+00_dp, 5.507553238667400e+03_dp), &
       terms(6.177000000000000e-08_dp, 1.73038850355e+00_dp, 5.223693919802200e+03_dp) /)

  TYPE (terms), PARAMETER :: R1(40) = (/ &
       terms(1.000139887990000e+00_dp, 0.00000000000e+00_dp, 0.000000000000000e+00_dp), &
       terms(1.670699626000000e-02_dp, 3.09846350771e+00_dp, 6.283075849991400e+03_dp), &
       terms(1.395602300000000e-04_dp, 3.05524609620e+00_dp, 1.256615169998280e+04_dp), &
       terms(3.083720000000000e-05_dp, 5.19846674381e+00_dp, 7.771377146812050e+04_dp), &
       terms(1.628461000000000e-05_dp, 1.17387749012e+00_dp, 5.753384884896800e+03_dp), &
       terms(1.575568000000000e-05_dp, 2.84685245825e+00_dp, 7.860419392439200e+03_dp), &
       terms(9.247990000000000e-06_dp, 5.45292234084e+00_dp, 1.150676976979360e+04_dp), &
       terms(5.424440000000000e-06_dp, 4.56409149777e+00_dp, 3.930209696219600e+03_dp), &
       terms(4.721100000000000e-06_dp, 3.66100022149e+00_dp, 5.884926846583200e+03_dp), &
       terms(3.459830000000000e-06_dp, 9.63686176870e-01_dp, 5.507553238667400e+03_dp), &
       terms(3.287800000000000e-06_dp, 5.89983646482e+00_dp, 5.223693919802200e+03_dp), &
       terms(3.067840000000000e-06_dp, 2.98671395120e-01_dp, 5.573142801433100e+03_dp), &
       terms(2.431890000000000e-06_dp, 4.27349536153e+00_dp, 1.179062908865880e+04_dp), &
       terms(2.118290000000000e-06_dp, 5.84714540314e+00_dp, 1.577343542447800e+03_dp), &
       terms(1.857520000000000e-06_dp, 5.02194447178e+00_dp, 1.097707880469900e+04_dp), &
       terms(1.748440000000000e-06_dp, 3.01193636534e+00_dp, 1.884922754997420e+04_dp), &
       terms(1.098350000000000e-06_dp, 5.05510636285e+00_dp, 5.486777843175000e+03_dp), &
       terms(9.831599999999999e-07_dp, 8.86813112770e-01_dp, 6.069776754553400e+03_dp), &
       terms(8.649900000000000e-07_dp, 5.68959778254e+00_dp, 1.572083878487840e+04_dp), &
       terms(8.582500000000000e-07_dp, 1.27083733351e+00_dp, 1.610006857376741e+05_dp), &
       terms(6.490300000000000e-07_dp, 2.72506137870e-01_dp, 1.726015465469040e+04_dp), &
       terms(6.291600000000000e-07_dp, 9.21771088320e-01_dp, 5.296909650946000e+02_dp), &
       terms(5.705600000000000e-07_dp, 2.01374292014e+00_dp, 8.399684731811189e+04_dp), &
       terms(5.573600000000000e-07_dp, 5.24159798933e+00_dp, 7.143069561812909e+04_dp), &
       terms(4.938400000000000e-07_dp, 3.24501240359e+00_dp, 2.544314419883400e+03_dp), &
       terms(4.696300000000000e-07_dp, 2.57805070386e+00_dp, 7.755226113240000e+02_dp), &
       terms(4.466100000000000e-07_dp, 5.53715807302e+00_dp, 9.437762934887000e+03_dp), &
       terms(4.251500000000000e-07_dp, 6.01110242003e+00_dp, 6.275962302990600e+03_dp), &
       terms(3.896800000000000e-07_dp, 5.36071738169e+00_dp, 4.694002954707600e+03_dp), &
       terms(3.824500000000000e-07_dp, 2.39255343974e+00_dp, 8.827390269874801e+03_dp), &
       terms(3.749000000000000e-07_dp, 8.29529223320e-01_dp, 1.965104848109800e+04_dp), &
       terms(3.695700000000000e-07_dp, 4.90107591914e+00_dp, 1.213955350910680e+04_dp), &
       terms(3.566000000000000e-07_dp, 1.67468058995e+00_dp, 1.203646073488820e+04_dp), &
       terms(3.453700000000000e-07_dp, 1.84270693282e+00_dp, 2.942463423291600e+03_dp), &
       terms(3.319300000000000e-07_dp, 2.43703000980e-01_dp, 7.084896781115200e+03_dp), &
       terms(3.192100000000000e-07_dp, 1.83682297810e-01_dp, 5.088628839766800e+03_dp), &
       terms(3.184600000000000e-07_dp, 1.77775642085e+00_dp, 3.981490034082000e+02_dp), &
       terms(2.846400000000000e-07_dp, 1.21344868176e+00_dp, 6.286598968340400e+03_dp), &
       terms(2.779300000000000e-07_dp, 1.89934330904e+00_dp, 6.279552731642400e+03_dp), &
       terms(2.627500000000000e-07_dp, 4.58896850401e+00_dp, 1.044738783960440e+04_dp) /)

  TYPE (terms), PARAMETER :: R2(10) = (/ &
       terms(1.030186080000000e-03_dp, 1.10748969588e+00_dp, 6.283075849991400e+03_dp), &
       terms(1.721238000000000e-05_dp, 1.06442301418e+00_dp, 1.256615169998280e+04_dp), &
       terms(7.022150000000000e-06_dp, 3.14159265359e+00_dp, 0.000000000000000e+00_dp), &
       terms(3.234600000000000e-07_dp, 1.02169059149e+00_dp, 1.884922754997420e+04_dp), &
       terms(3.079900000000000e-07_dp, 2.84353804832e+00_dp, 5.507553238667400e+03_dp), &
       terms(2.497100000000000e-07_dp, 1.31906709482e+00_dp, 5.223693919802200e+03_dp), &
       terms(1.848500000000000e-07_dp, 1.42429748614e+00_dp, 1.577343542447800e+03_dp), &
       terms(1.007800000000000e-07_dp, 5.91378194648e+00_dp, 1.097707880469900e+04_dp), &
       terms(8.654000000000001e-08_dp, 1.42046854427e+00_dp, 6.275962302990600e+03_dp), &
       terms(8.634000000000000e-08_dp, 2.71461506020e-01_dp, 5.486777843175000e+03_dp) /)

  TYPE (terms), PARAMETER :: R3(6) = (/ &
       terms(4.359385000000000e-05_dp, 5.78455133738e+00_dp, 6.283075849991400e+03_dp), &
       terms(1.236330000000000e-06_dp, 5.57934722157e+00_dp, 1.256615169998280e+04_dp), &
       terms(1.234100000000000e-07_dp, 3.14159265359e+00_dp, 0.000000000000000e+00_dp), &
       terms(8.792000000000000e-08_dp, 3.62777733395e+00_dp, 7.771377146812050e+04_dp), &
       terms(5.689000000000000e-08_dp, 1.86958905084e+00_dp, 5.573142801433100e+03_dp), &
       terms(3.301000000000000e-08_dp, 5.47027913302e+00_dp, 1.884922754997420e+04_dp) /)

  TYPE (terms), PARAMETER :: R4(2) = (/ &
       terms(1.445950000000000e-06_dp, 4.27319435148e+00_dp, 6.283075849991400e+03_dp), &
       terms(6.729000000000000e-08_dp, 3.91697608662e+00_dp, 1.256615169998280e+04_dp) /)

  TYPE (terms), PARAMETER :: R5(1) = (/ &
       terms(3.858000000000000e-08_dp, 2.56384387339e+00_dp, 6.283075849991400e+03_dp) /)

CONTAINS

  FUNCTION sum_vsop87 (t, term, order)

    REAL(dp) :: sum_vsop87

    REAL(dp), INTENT(in) :: t
    INTEGER, INTENT(in) :: order
    TYPE (terms), INTENT(in) :: term(:)

    REAL(dp) :: td10, tn
    INTEGER :: k

    td10 = t*0.1_dp
    sum_vsop87 = 0.0_dp

    IF (.NOT.(order == 1 .AND. td10 == 0.0_dp)) THEN
      tn = td10**(order-1)
    ELSE
      tn = 1.0_dp
    ENDIF

    DO k = 1, UBOUND(term,1)
      sum_vsop87 = sum_vsop87+term(k)%A*COS(term(k)%B+term(k)%C*td10)
    END DO
    sum_vsop87 = sum_vsop87*tn

  END FUNCTION sum_vsop87

  SUBROUTINE earth_position (t, l, b, r)

    REAL(dp), INTENT(in)  :: t  ! number of centuries since J2000

    ! Spherical coordinates, FK5

    REAL(dp), INTENT(out) :: l  ! longitude [radians]
    REAL(dp), INTENT(out) :: b  ! latitude  [radians]
    REAL(dp), INTENT(out) :: r  ! radius    [AU] (AU - Astronomical Units)

    REAL(dp) :: ld

    ! Calculate L, B, R

    l   = 0.0_dp; b   = 0.0_dp; r   = 0.0_dp; 

    ! Calculate L

    l = l+sum_vsop87 (t, L1, 1)
    l = l+sum_vsop87 (t, L2, 2)
    l = l+sum_vsop87 (t, L3, 3)
    l = l+sum_vsop87 (t, L4, 4)
    l = l+sum_vsop87 (t, L5, 5)
    l = l+sum_vsop87 (t, L6, 6)

    l = MODULO(l, 2.0_dp*pi)

    ! Calculate B

    b = b+sum_vsop87 (t, B1, 1)
    b = b+sum_vsop87 (t, B2, 2)

    ! Calculate R

    r = r+sum_vsop87 (t, R1, 1)
    r = r+sum_vsop87 (t, R2, 2)
    r = r+sum_vsop87 (t, R3, 3)
    r = r+sum_vsop87 (t, R4, 4)
    r = r+sum_vsop87 (t, R5, 5)

    ! Convert from Dynamic to FK5 equator & ecliptic (need extra factor 0.1
    ! for t in the equation for ld compared to Meeus original development). 

    ld = l-0.1_dp*t*(1.397_dp+0.000031_dp*t)*dtor;
    l  = l+(-0.09033_dp+0.03916_dp*TAN(b)*(COS(ld)+SIN(ld)))*stor;
    b  = b+0.03916_dp*(COS(ld)-SIN(ld))*stor;

  END SUBROUTINE earth_position

END MODULE mo_vsop87