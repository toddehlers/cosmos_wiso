SUBROUTINE SURRTAB

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** AER'S RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!     -----------------------------------------------------------------

!     mag, MPI, 25 February 2000: comment added

USE MO_KIND  , ONLY : DP

!   Output
USE MO_RRTAB , ONLY : TRANS, BPADE

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER :: ITR

!     LOCAL REAL SCALARS
REAL(DP) :: ZTAU, ZTFN


BPADE=1._DP/0.278_DP
TRANS(0)   =1._DP
TRANS(5000)=0._DP
DO ITR=1,4999
  ZTFN=REAL(ITR,dp)/5000._DP
  ZTAU=BPADE*ZTFN/(1._DP-ZTFN)
  TRANS(ITR)=EXP(-ZTAU)
ENDDO

!CDIR DU_UPDATE(TRANS)

!     -----------------------------------------------------------------

RETURN
END SUBROUTINE SURRTAB
