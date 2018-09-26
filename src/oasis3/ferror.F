      SUBROUTINE ferror
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL C *
C               * -------------     ------- *
C               *****************************
C
C**** *ferror*  - Signal handler
C
C     Purpose:
C     -------
C     ferror is executed each time signal sigcld is caught.
C     It stops execution of run due to child process death.
C
C**   Interface:
C     ---------
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
C     wait, waitcld
C
C     Reference:
C     ---------
C     See OASIS manual (1997)
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       1.0       L. Terray      94/01/01  created
C       2.0       L. Terray      95/12/20  modified: new structure
C       2.2       L. Terray      97/10/10  added: Include nproc.h
C       2.2       S. Valcke, L.T 97/12/24  added: call wait and waitcld 
C       2.3       L. Terray      99/09/15  added: GMEM branch
C       2.5       S. Valcke      2K/09/04  Remove cmach
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------- Include files and use of modules---------------------------
C
      USE mod_kinds_oasis
      USE mod_experiment
      USE mod_parameter
      USE mod_nproc
      USE mod_unit
      USE mod_hardware
C
C* ---------------------------- Local declarations ----------------------
C
      INTEGER (kind=ip_intwp_p) wait
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Signal sigcld has been caught; check which process is dead
C        ----------------------------------------------------------
C
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' ROUTINE ferror  -  Level C'
      WRITE (UNIT = nulou,FMT = *) ' **************     *******'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' Catch child process death'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
#ifdef use_comm_PIPE 
          iwone = wait(isone)
#elif defined use_comm_SIPC || defined use_comm_GMEM
          CALL waitcld (isone, iwone)
#endif
      DO 110 jm = 1, ig_nmodel
        IF (iwone .EQ. nproc(jm)) THEN 
            WRITE (UNIT = nulou,FMT = *) 
     $          ' WARNING: sigcld has been caught by coupler'
            WRITE (UNIT = nulou,FMT = *) 
     $          ' =======  ======                    ======='
            WRITE (UNIT = nulou,FMT = *) 
     $          ' model ', cmodnam(jm), 
     $          ' process pid = ',nproc(jm),' is dead'
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) 
     $          ' We STOP the simulation now !!!'
            WRITE (UNIT = nulou,FMT = *) ' '
        ENDIF 
 110  CONTINUE 
      CALL FLUSH (nulou)
C
C
C*    2. End of simulation if abnormal termination
C        -----------------------------------------
C
      CALL HALTE ('STOP in ferror')
C
C
C*    3. End of routine
C        --------------
C
      RETURN
      END
      SUBROUTINE fignore
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL C *
C               * -------------     ------- *
C               *****************************
C
C**** *fignore*  - Signal handler
C
C     Purpose:
C     -------
C     fignore is executed each time signal sigcld is caught.
C**   Interface:
C     ---------
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
C
C     Reference:
C     ---------
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       2.5       S. Valcke      01/08/02  Created
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------- Include files and use of modules-------------------   
C
      USE mod_unit
C
C* ---------------------------- Local declarations ----------------------
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' ROUTINE fignore -  Level C'
      WRITE (UNIT = nulou,FMT = *) ' ***************    *******'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
C
C*    End of routine
C     --------------
C
      RETURN
      END
