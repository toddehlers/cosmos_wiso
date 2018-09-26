MODULE MO_RRTRF


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *MO_RRTRF* - RRTM REFERENCE ATMOSPHERE
!     -----------------------------------------------------------------

REAL(DP), DIMENSION(59) :: PREF
REAL(DP), DIMENSION(59) :: PREFLOG
REAL(DP), DIMENSION(59) :: TREF

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! PREF   :  REAL    
! PREFLOG: REAL
! TREF   : REAL
!     -----------------------------------------------------------------
END MODULE MO_RRTRF
