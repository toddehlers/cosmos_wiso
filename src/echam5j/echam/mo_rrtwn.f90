MODULE MO_RRTWN


USE MO_KIND  , ONLY : DP

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

INTEGER , DIMENSION(16) :: NG
INTEGER , DIMENSION(16) :: NSPA
INTEGER , DIMENSION(16) :: NSPB

REAL(DP), DIMENSION(16) :: WAVENUM1
REAL(DP), DIMENSION(16) :: WAVENUM2
REAL(DP), DIMENSION(16) :: DELWAVE

REAL(DP), DIMENSION(181,16) :: TOTPLNK
REAL(DP), DIMENSION(181)    :: TOTPLK16

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----   : ----    : -------
!  NG     : INTEGER : Number of k-coefficients in spectral intervals
!  NSPA   : INTEGER :
!  NSPB   : INTEGER :
! WAVENUM1: REAL    : Lower wavenumber spectral limit
! WAVENUM2: REAL    : Higher wavenumber spectral limit
! DELWAVE : REAL    : Spectral interval width
! TOTPLNK : REAL    :
! TOTPLK16: REAL    :
!     -----------------------------------------------------------------
END MODULE MO_RRTWN
