MODULE MO_RRTAB


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

!    -------------------------------------------------------------------

REAL(DP), DIMENSION(0:5000) :: TRANS
!CDIR DUPLICATE(TRANS,256)

REAL(DP):: BPADE

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! TRANS  :  REAL    
! BPADE  :  REAL     
!     -----------------------------------------------------------------
END MODULE MO_RRTAB
