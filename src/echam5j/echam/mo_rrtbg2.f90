MODULE MO_RRTBG2

USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

!    -------------------------------------------------------------------

REAL(DP):: CORR1(0:200)
REAL(DP):: CORR2(0:200)

!CDIR DUPLICATE(CORR1,256)
!CDIR DUPLICATE(CORR2,256)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! CORR1  :  REAL   : 
! CORR2  :  REAL   :
!    -------------------------------------------------------------------
END MODULE MO_RRTBG2
