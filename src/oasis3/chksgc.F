      SUBROUTINE chksgc
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL C *
C               * -------------     ------- *
C               *****************************
C
C**** *chksgc*  - Initialize signal handling
C
C     Purpose:
C     -------
C     Initialize signal handling related to fork stuff and fpe trapping
C
C**   Interface:
C     ---------
C       *CALL*  *chksgc*
C
C     Input:
C     -----
C     None
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
C     fsigctl, signal, ferror, getfpe
C
C     Reference:
C     ---------
C     See OASIS manual (1997)
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       2.2       L. Terray      97/12/24  created
C       2.3       S. Valcke      99/04/30  added: printing levels
C       2.3       L. Terray      99/09/15  added: GMEM branch
C       2.5       S. Valcke      2K/09/04  Remove cmach
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------------------- Include files ---------------------------
C
      USE mod_unit
      USE mod_hardware
      USE mod_printing
C
C* ---------------------------- External declarations -------------------
C
      EXTERNAL ferror, getfpe
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C
C*    1. Initialization
C        --------------
C
      IF (nlogprt .GE. 2) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '           ROUTINE chksgc  -  Level C'
          WRITE (UNIT = nulou,FMT = *) 
     $    '           **************     *******'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' Initialize signal handling '
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
      ENDIF
C
C
C*    2. Signal handling initialization
C        ------------------------------
C 1- Insure proper termination in case of child process death
C 2- Catch floating point error within OASIS
C
C* Case PIPE-CRAY
C
#ifdef use_comm_PIPE
          CALL fsigctl ('REGISTER', 'SIGCLD', ferror)
          CALL fsigctl ('REGISTER', 'SIGFPE', getfpe)
C
C* Case SIPC,GMEM-IEEE
C
#elif defined use_comm_SIPC || defined use_comm_GMEM
          CALL fpecatch

C
C* Child process death trapping only for SIPC (taken care thru MPI for GMEM)
C
#ifdef use_comm_SIPC
             CALL cldcatch
#endif
C
C* Others
C
#else
          IF (nlogprt .GE. 2) THEN 
              WRITE (UNIT = nulou,FMT = *) 
     $            ' CLIM or exotic cases ===>>> No signal handling' 
          ENDIF
#endif
C
C
C*    3. End of routine
C        --------------
C
      IF (nlogprt .GE. 2) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '          --------- End of routine chksgc ---------'
          CALL FLUSH (nulou)
      ENDIF
      RETURN
      END
