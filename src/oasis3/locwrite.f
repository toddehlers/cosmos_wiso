      SUBROUTINE locwrite (cdfldn, pfield, kdimax, knulre, kflgre)
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 0 *
C               * -------------     ------- *
C               *****************************
C
C**** *locwrite*  - Write binary field on unit knulre
C
C     Purpose:
C     -------
C     Write string cdfldn and array pfield on unit knulre
C
C**   Interface:
C     ---------
C       *CALL*  *locwrite (cdfldn, pfield, kdimax, knulre, kflgre)*
C
C     Input:
C     -----
C                cdfldn : character string locator
C                kdimax : dimension of field to be written 
C                knulre : logical unit to be written
C                pfield : field array (real 1D) 
C
C     Output:
C     ------
C                kflgre : error status flag
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
C       2.3       S. Valcke      99/04/30  added: printing levels
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------------------- Include files ---------------------------
C
      USE mod_kinds_oasis
      USE mod_unit
      USE mod_printing
C
C* ---------------------------- Argument declarations -------------------
C
      REAL (kind=ip_realwp_p) pfield(kdimax)
      CHARACTER*8 cdfldn
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Initialization
C        --------------
C
      IF (nlogprt .GE. 2) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '           ROUTINE locwrite  -  Level 0'
          WRITE (UNIT = nulou,FMT = *) 
     $    '           ****************     *******'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = 1001) knulre
          WRITE (UNIT = nulou,FMT = *) ' '
      ENDIF
C
C* Formats
C
 1001 FORMAT(5X,' Write binary file connected to unit = ',I3)
C
C     2. Find field in file
C        ------------------
C
C* Write string
      WRITE (UNIT = knulre, ERR = 210) cdfldn
C* Write associated field
      WRITE (UNIT = knulre, ERR = 210) pfield
C* Writing done and ok
      kflgre = 0
      GO TO 220
C* Problem in Writing
 210  kflgre = 1
 220  CONTINUE
C
C
C*    3. End of routine
C        --------------
C
      IF (nlogprt .GE. 2) THEN
          WRITE (UNIT = nulou,FMT = *) 
     $    '          --------- End of routine locwrite ---------'
          WRITE (UNIT = nulou,FMT = *) ' '
          CALL FLUSH (nulou)
      ENDIF
      RETURN
      END


