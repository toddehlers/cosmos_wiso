      SUBROUTINE izero (ka, kna)
C****
C               ********************************
C               * OASIS SUBROUTINE  -  LEVEL T *
C               * ----------------     ------- *
C               ********************************
C
C**** *izero*  - Utility routine
C
C     Purpose:
C     -------
C     Zero integer array from element 1 to element kna
C
C**   Interface:
C     ---------
C       *CALL*  *izero (ka, kna)*
C
C     Input:
C     -----
C                ka     : array to be zeroed (integer 1D)
C                kna    : array dimension (integer)
C
C     Output:
C     ------
C                ka     : array filled up with zeros (integer 1D)
C
C     Workspace:
C     ---------
C     None
C
C     Externals:
C     ---------
C     None
C
C     Reference:
C     ---------
C     See OASIS manual (1995)
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       2.0       L. Terray      95/09/01  created
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------------------- Include files ---------------------------
C
      USE mod_kinds_oasis
      USE mod_unit
C
C* ---------------------------- Argument declarations -------------------
C
      INTEGER (kind=ip_intwp_p) ka(kna)
C
C* ---------------------------- Local declarations ----------------------
C
      INTEGER (kind=ip_intwp_p) inzero
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Zero the array 
C        ---------------
C
      inzero = 0
      DO 110 ja = 1, kna
        ka(ja) = inzero
  110 CONTINUE
C
C
C*    2. End of routine
C        --------------
C
      RETURN
      END
