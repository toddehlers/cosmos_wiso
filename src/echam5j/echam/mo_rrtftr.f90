MODULE MO_RRTFTR

USE MO_KIND   , ONLY : DP
USE MO_PARRRTM, ONLY : JPBAND, JPGPT, JPG


IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

!    -------------------------------------------------------------------

INTEGER :: NGC(JPBAND)
INTEGER :: NGS(JPBAND)
INTEGER :: NGN(JPGPT)
INTEGER :: NGB(JPGPT)

INTEGER :: NGM(JPG*JPBAND)
REAL(DP):: WT(JPG)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
!  NGC   : INTEGER :
!  NGS   : INTEGER :
!  NGN   : INTEGER :
!  NGB   : INTEGER :
!  NGM   : INTEGER :
!  WT    : REAL    :
!    -------------------------------------------------------------------
END MODULE MO_RRTFTR
