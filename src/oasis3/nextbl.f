      FUNCTION nextbl (cdstr, kistr, knstr)
C****
C               ******************************
C               * OASIS FUNCTION  -  LEVEL T *
C               * --------------     ------- *
C               ******************************
C
C**** *nextbl*  - Search function
C
C     Purpose:
C     -------
C     Find the first blank in a character string
C
C**   Interface:
C     ---------
C       *ii =*  *nextbl (cdstr, kistr, knstr)*
C
C     Input:
C     -----
C                cdstr : string to be searched (char string)
C                kistr : initial search position within the string (integer)
C                knstr : final search position within the string (integer)
C
C     Output:
C     ------
C     None
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
      USE mod_unit
C
C* ---------------------------- Argument declarations -------------------
C
      CHARACTER*1 cdstr
      DIMENSION cdstr(knstr)
C
C* ---------------------------- Local declarations -------------------
C
      CHARACTER (len=1), SAVE :: clblank = ' '
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C
C*    1. Find first blank character
C        --------------------------
C
      DO 110 ji = kistr, knstr
        idum = ji
        IF (cdstr(ji) .EQ. clblank) GO TO 120
  110 CONTINUE
  120 CONTINUE
      nextbl = idum
      IF (idum .GE. knstr) nextbl = -1
C
C
C*    2. End of function
C        ---------------
C
      RETURN
      END
