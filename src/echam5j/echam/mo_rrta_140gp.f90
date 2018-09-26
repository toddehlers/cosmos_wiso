MODULE MO_RRTA1


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA1* - RRTM COEFFICIENTS FOR INTERVAL 1
!     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG1  = 8

REAL(DP):: FRACREFA(NG1)  , FRACREFB(NG1)
REAL(DP):: KA(5,13,NG1)   , ABSA(65,NG1)
REAL(DP):: KB(5,13:59,NG1), ABSB(235,NG1)
REAL(DP):: SELFREF(10,NG1), FORREF(NG1)
!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(ABSB,256)
!CDIR DUPLICATE(SELFREF,256)


EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! FORREF  : REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE MO_RRTA1
MODULE MO_RRTA2


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA2* - RRTM COEFFICIENTS FOR INTERVAL 2
!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG2  = 14

!     The ith set of reference fractions are from the ith reference
!     pressure level.
REAL(DP):: FRACREFA(NG2,13), FRACREFB(NG2), REFPARAM(13)
REAL(DP):: KA(5,13,NG2)   , ABSA(65,NG2)
REAL(DP):: KB(5,13:59,NG2), ABSB(235,NG2)
REAL(DP):: SELFREF(10,NG2), FORREF(NG2)

!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(ABSB,256)
!CDIR DUPLICATE(SELFREF,256)
!CDIR DUPLICATE(FRACREFA,256)


EQUIVALENCE (KA(1,1,1),ABSA(1,1)),(KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! REFPARAM: REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL
! FORREF  : REAL     
!     -----------------------------------------------------------------
END MODULE MO_RRTA2


MODULE MO_RRTA3


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA3* - RRTM COEFFICIENTS FOR INTERVAL 3
!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG3  = 16

REAL(DP):: FRACREFA(NG3,10) ,FRACREFB(NG3,5)

REAL(DP), DIMENSION(16) :: FORREF
REAL(DP), DIMENSION(16) :: ABSN2OA
REAL(DP), DIMENSION(16) :: ABSN2OB
REAL(DP), DIMENSION(10) :: ETAREF
REAL(DP), DIMENSION(59) :: H2OREF
REAL(DP), DIMENSION(59) :: N2OREF
REAL(DP), DIMENSION(59) :: CO2REF

REAL(DP):: KA(10,5,13,NG3)  ,ABSA(650,NG3)
REAL(DP):: KB(5,5,13:59,NG3),ABSB(1175,NG3)
REAL(DP):: SELFREF(10,NG3)
REAL(DP):: STRRAT

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,1,13,1),ABSB(1,1))

!     ------------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! ABSN2OA : REAL
! ABSN2OB : REAL
! CO2REF  : REAL
! ETAREF  : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! H2OREF  : REAL
! KA      : REAL     
! KB      : REAL     
! N2OREF  : REAL
! SELFREF : REAL     
! STRRAT  : REAL
!     -----------------------------------------------------------------
!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(ABSB,256)
!CDIR DUPLICATE(SELFREF,256)
!CDIR DUPLICATE(FRACREFA,256)
!CDIR DUPLICATE(FRACREFB,256)
!CDIR DUPLICATE(H2OREF,256)
!CDIR DUPLICATE(N2OREF,256)
!CDIR DUPLICATE(CO2REF,256)
!CDIR DUPLICATE(ETAREF,256)

END MODULE MO_RRTA3
MODULE MO_RRTA4


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA4* - RRTM COEFFICIENTS FOR INTERVAL 5
!     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG4  = 14

REAL(DP):: FRACREFA(NG4,9)  ,FRACREFB(NG4,6)
REAL(DP):: KA(9,5,13,NG4)   ,ABSA(585,NG4)
REAL(DP):: KB(6,5,13:59,NG4),ABSB(1410,NG4)
REAL(DP):: SELFREF(10,NG4)
REAL(DP):: STRRAT1
REAL(DP):: STRRAT2

!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(ABSB,256)
!CDIR DUPLICATE(SELFREF,256)
!CDIR DUPLICATE(FRACREFA,256)
!CDIR DUPLICATE(FRACREFB,256)

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,1,13,1),ABSB(1,1))

!     ------------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE MO_RRTA4
MODULE MO_RRTA5


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA5* - RRTM COEFFICIENTS FOR INTERVAL 5
!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG5  = 16

REAL(DP):: FRACREFA(NG5,9) ,FRACREFB(NG5,5)

REAL(DP), DIMENSION(NG5) :: CCL4

REAL(DP):: KA(9,5,13,NG5)   ,ABSA(585,NG5)
REAL(DP):: KB(5,5,13:59,NG5),ABSB(1175,NG5)
REAL(DP):: SELFREF(10,NG5)
REAL(DP):: STRRAT1
REAL(DP):: STRRAT2

!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(ABSB,256)
!CDIR DUPLICATE(SELFREF,256)
!CDIR DUPLICATE(FRACREFA,256)
!CDIR DUPLICATE(FRACREFB,256)

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! CCL4    : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
! STRRAT1 : REAL
! STRRAT2 : REAL    
!     -----------------------------------------------------------------
END MODULE MO_RRTA5
MODULE MO_RRTA6


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA6* - RRTM COEFFICIENTS FOR INTERVAL 6
!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG6  = 8

REAL(DP), DIMENSION(NG6) :: FRACREFA

REAL(DP), DIMENSION(NG6) :: CFC11ADJ
REAL(DP), DIMENSION(NG6) :: CFC12
REAL(DP), DIMENSION(NG6) :: ABSCO2

REAL(DP):: KA(5,13,NG6),ABSA(65,NG6)
REAL(DP):: SELFREF(10,NG6)

!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(SELFREF,256)

EQUIVALENCE (KA(1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSCO2  : REAL 
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! CFC11ADJ: REAL
! CFC12   : REAL
! KA      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE MO_RRTA6
MODULE MO_RRTA7


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA7* - RRTM COEFFICIENTS FOR INTERVAL 7
!     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG7  = 12

REAL(DP):: FRACREFA(NG7,9)

REAL(DP), DIMENSION(NG7) :: FRACREFB
REAL(DP), DIMENSION(NG7) :: ABSCO2
REAL(DP):: KA(9,5,13,NG7) ,ABSA(585,NG7)
REAL(DP):: KB(5,13:59,NG7),ABSB(235,NG7)
REAL(DP):: SELFREF(10,NG7)

!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(ABSB,256)
!CDIR DUPLICATE(SELFREF,256)
!CDIR DUPLICATE(FRACREFA,256)
!CDIR DUPLICATE(FRACREFB,256)

REAL(DP):: STRRAT

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL 
! ABSCO2  : REAL
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL  
! STRRAT  : REAL   
!     -----------------------------------------------------------------
END MODULE MO_RRTA7
MODULE MO_RRTA8


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA8* - RRTM COEFFICIENTS FOR INTERVAL 8
!     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG8  = 8

REAL(DP), DIMENSION(NG8) :: FRACREFA
REAL(DP), DIMENSION(NG8) :: FRACREFB
REAL(DP), DIMENSION(NG8) :: CFC12
REAL(DP), DIMENSION(NG8) :: CFC22ADJ
REAL(DP), DIMENSION(NG8) :: ABSCO2A
REAL(DP), DIMENSION(NG8) :: ABSCO2B
REAL(DP), DIMENSION(NG8) :: ABSN2OA
REAL(DP), DIMENSION(NG8) :: ABSN2OB
REAL(DP), DIMENSION(59)  :: H2OREF
REAL(DP), DIMENSION(59)  :: N2OREF
REAL(DP), DIMENSION(59)  :: O3REF

REAL(DP):: KA(5,7,NG8)    ,ABSA(35,NG8)
REAL(DP):: KB(5,7:59,NG8) ,ABSB(265,NG8)
REAL(DP):: SELFREF(10,NG8)

!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(ABSB,256)
!CDIR DUPLICATE(SELFREF,256)
!CDIR DUPLICATE(H2OREF,256)
!CDIR DUPLICATE(N2OREF,256)
!CDIR DUPLICATE(O3REF,256)
EQUIVALENCE (KA(1,1,1),ABSA(1,1)),(KB(1,7,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! ABSCO2A : REAL     
! ABSCO2B : REAL     
! ABSN2OA : REAL     
! ABSN2OB : REAL 
! CFC12   : REAL     
! CFC22ADJ: REAL     
! FRACREFA: REAL    
! FRACREFB: REAL    
! H2OREF  : REAL    
! KA      : REAL     
! KB      : REAL     
! N2OREF  : REAL    
! O3REF   : REAL    
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE MO_RRTA8
MODULE MO_RRTA9


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA9* - RRTM COEFFICIENTS FOR INTERVAL 9
!     BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG9  = 12

REAL(DP):: FRACREFA(NG9,9)

REAL(DP), DIMENSION(NG9) :: FRACREFB
REAL(DP), DIMENSION(13) :: N2OREF
REAL(DP), DIMENSION(13) :: H2OREF
REAL(DP), DIMENSION(13) :: CH4REF
REAL(DP), DIMENSION(11) :: ETAREF
! 36 = 3*NG9      
REAL(DP), DIMENSION(36) :: ABSN2O

REAL(DP):: KA(11,5,13,NG9) ,ABSA(715,NG9)
REAL(DP):: KB(5,13:59,NG9) ,ABSB(235,NG9)
REAL(DP):: SELFREF(10,NG9)
REAL(DP):: STRRAT

!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(ABSB,256)
!CDIR DUPLICATE(SELFREF,256)
!CDIR DUPLICATE(ABSN2O,256)
!CDIR DUPLICATE(FRACREFA,256)
!CDIR DUPLICATE(N2OREF,256)
!CDIR DUPLICATE(H2OREF,256)
!CDIR DUPLICATE(CH4REF,256)
!CDIR DUPLICATE(ETAREF,256)


EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! ABSN2O  : REAL    
! CH4REF  : REAL
! ETAREF  : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! H2OREF  : REAL
! N2OREF  : REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
! STRRAT  : REAL
!     -----------------------------------------------------------------
END MODULE MO_RRTA9
MODULE MO_RRTA10


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA14* - RRTM COEFFICIENTS FOR INTERVAL 10
!     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG10 = 6

REAL(DP), DIMENSION(NG10) :: FRACREFA
REAL(DP), DIMENSION(NG10) :: FRACREFB

REAL(DP):: KA(5,13,NG10)   , ABSA(65,NG10)
REAL(DP):: KB(5,13:59,NG10), ABSB(235,NG10)

!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(ABSB,256)

EQUIVALENCE (KA(1,1,1),ABSA(1,1)),(KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
!     -----------------------------------------------------------------
END MODULE MO_RRTA10
MODULE MO_RRTA11


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA11* - RRTM COEFFICIENTS FOR INTERVAL 11
!     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG11 = 8

REAL(DP), DIMENSION(NG11) :: FRACREFA
REAL(DP), DIMENSION(NG11) :: FRACREFB

REAL(DP):: KA(5,13,NG11)   , ABSA(65,NG11)
REAL(DP):: KB(5,13:59,NG11), ABSB(235,NG11)
REAL(DP):: SELFREF(10,NG11)

!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(ABSB,256)
!CDIR DUPLICATE(SELFREF,256)

EQUIVALENCE (KA(1,1,1),ABSA(1,1)),(KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE MO_RRTA11
MODULE MO_RRTA12


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA12* - RRTM COEFFICIENTS FOR INTERVAL 12
!     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG12 = 8

REAL(DP):: FRACREFA(NG12,9)
REAL(DP):: KA(9,5,13,NG12) ,ABSA(585,NG12)
REAL(DP):: SELFREF(10,NG12)

!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(SELFREF,256)
!CDIR DUPLICATE(FRACREFA,256)

REAL(DP):: STRRAT

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL
! STRRAT1 : REAL     
!     -----------------------------------------------------------------
END MODULE MO_RRTA12
MODULE MO_RRTA13


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA13* - RRTM COEFFICIENTS FOR INTERVAL 13
!     BAND 13:  2080-2250 cm-1 (low - H2O,N2O; high - nothing)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG13 = 4

REAL(DP):: FRACREFA(NG13,9)

REAL(DP):: KA(9,5,13,NG13) ,ABSA(585,NG13)
REAL(DP):: SELFREF(10,NG13)
REAL(DP):: STRRAT

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1))

!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(SELFREF,256)
!CDIR DUPLICATE(FRACREFA,256)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL
! STRRAT1 : REAL     
!     -----------------------------------------------------------------
END MODULE MO_RRTA13
MODULE MO_RRTA14


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA14* - RRTM COEFFICIENTS FOR INTERVAL 14
!     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG14 = 2

REAL(DP), DIMENSION(NG14) :: FRACREFA
REAL(DP), DIMENSION(NG14) :: FRACREFB

REAL(DP):: KA(5,13,NG14)   ,ABSA(65,NG14)
REAL(DP):: KB(5,13:59,NG14),ABSB(235,NG14)
REAL(DP):: SELFREF(10,NG14)

!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(ABSB,256)
!CDIR DUPLICATE(SELFREF,256)

EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE MO_RRTA14


MODULE MO_RRTA15


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA15* - RRTM COEFFICIENTS FOR INTERVAL 15
!     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; high - nothing)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG15 = 2

REAL(DP):: FRACREFA(NG15,9)

REAL(DP):: KA(9,5,13,NG15) ,ABSA(585,NG15)
REAL(DP):: SELFREF(10,NG15)
REAL(DP):: STRRAT

!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(SELFREF,256)
!CDIR DUPLICATE(FRACREFA,256)

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL 
! STRRAT  : REAL    
!     -----------------------------------------------------------------
END MODULE MO_RRTA15
MODULE MO_RRTA16


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTA16* - RRTM COEFFICIENTS FOR INTERVAL 16
!     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)
!     -----------------------------------------------------------------

INTEGER, PARAMETER :: NG16 = 2

REAL(DP):: FRACREFA(NG16,9)

REAL(DP):: KA(9,5,13,NG16) ,ABSA(585,NG16)
REAL(DP):: SELFREF(10,NG16)
REAL(DP):: STRRAT

!CDIR DUPLICATE(ABSA,256)
!CDIR DUPLICATE(SELFREF,256)
!CDIR DUPLICATE(FRACREFA,256)

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL     
! STRRAT  : REAL
!     -----------------------------------------------------------------
END MODULE MO_RRTA16
