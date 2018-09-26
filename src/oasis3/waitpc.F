      SUBROUTINE waitpc
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL C *
C               * -------------     ------- *
C               *****************************
C
C**** *waitpc*  - Insure proper termination of simulation
C
C     Purpose:
C     -------
C     Wait for gcm's termination so that child processes end up
C     before their parent.
C     Use system call wait to pause calling process until
C     completion of child processes.
C
C**   Interface:
C     ---------
C       *CALL*  *waitpc*
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
C     wait, waitcld, SIPC_End, CLIM_Quit
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
C       2.0       L. Terray      95/09/01  modified: new structure
C       2.2       S. Valcke      97/09/03  added: introduction of SVIPC
C       2.2       L. Terray      97/12/24  added: change waitsipc routine
C       2.3       S. Valcke      99/04/30  added: printing levels
C       2.3       L. Terray      99/09/15  added: GMEM branch. replace
C                                                 SIPC_Destroy by SIPC_End
C       2.5       S. Valcke      2K/09/04  Remove cmach
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------- Include files and USE of modules---------------------------
C
      USE mod_kinds_oasis
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
      USE mod_clim
#endif
      USE mod_parameter
      USE mod_sipc
      USE mod_experiment
      USE mod_pipe
      USE mod_nproc
      USE mod_unit
      USE mod_hardware
      USE mod_printing
C
C* ---------------------------- Local declarations ----------------------
C
      INTEGER (kind=ip_intwp_p) wait
      INTEGER (kind=ip_intwp_p) iwone
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
     $    '           ROUTINE waitpc  -  Level C'
          WRITE (UNIT = nulou,FMT = *) 
     $    '           **************     *******'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' Wait until gcm completion '
          WRITE (UNIT = nulou,FMT = *) ' Exit if interpolator only '
          WRITE (UNIT = nulou,FMT = *) ' '
      ENDIF
C
C
C*    2. Wait until gcm's completion  
C        ---------------------------
C
      iwone = 0
      IF (cchan .EQ. 'PIPE' .OR. cchan .EQ. 'SIPC'
     $    .OR. cchan .EQ. 'GMEM') THEN 
          DO 210 jm = 1, ig_nmodel
#ifdef use_comm_PIPE            
            iwone = wait(isone)
#elif defined use_comm_SIPC || defined use_comm_GMEM   
            CALL waitcld(isone,iwone)
#endif
            IF (iwone .eq. -1) THEN
                WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
                WRITE (UNIT = nulou,FMT = *) 
     $              ' ===>>> Model numbered = ', jm
                WRITE (UNIT = nulou,FMT = *) 
     $              ' ======         ====        ====     '
                WRITE (UNIT = nulou,FMT = *) 
     $          ' was not an running child when waitpc was called '
                CALL FLUSH(nulou)
            ELSE
                DO 220 jk = 1, ig_nmodel
                  IF (iwone .EQ. nproc(jk)) THEN
                      IF (nlogprt .GE. 1) THEN 
                          CALL prcout
     $                    ('End of simulation for model named', 
     $                     cmodnam(jk), 1)
                          CALL prtout
     $                    ('Exit of process nproc =', nproc(jk), 1)
                          WRITE (UNIT = nulou,FMT = *) 
     $                    ' coupler waits until end of other gcms '
                      ENDIF
                      GO TO 210
                  ENDIF 
 220            CONTINUE 
                WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
                WRITE (UNIT = nulou,FMT = *) 
     $              ' ===>>> : first wait return code is = '
     $              ,iwone
                WRITE (UNIT = nulou,FMT = *) 
     $              ' ======         ====        ====      '
                WRITE (UNIT = nulou,FMT = *) ' '
                WRITE (UNIT = nulou,FMT = *) 
     $              ' Either gcm is probably multitasked '
                WRITE (UNIT = nulou,FMT = *) 
     $              ' The leaving process with pid = ', iwone
                WRITE (UNIT = nulou,FMT = *) 
     $              ' is different from the initial process values  '
                WRITE (UNIT = nulou,FMT = *) 
     $              ' We go on til normal termination !!!'
             ENDIF
 210      CONTINUE
C
C* Destroy shared-memory pools used for exchanging fields and
C  pool for reading/writing initial infos from/to models
C
#if defined use_comm_SIPC || use_comm_GMEM
          CALL SIPC_End(ig_nmodel)
#endif
C
       ENDIF
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
       CALL CLIM_Quit( infos)
       IF (infos .NE. CLIM_Ok) THEN
          CALL prtout
     $         ('An error occured while leaving CLIM. Error =',
     $         infos, 2)
       ENDIF
#endif
C
C
C
C*    3. End of program
C        --------------
C
      IF (nlogprt .GE. 2) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '          --------- End of routine waitpc ---------'
          CALL FLUSH (nulou)
      ENDIF
      RETURN
      END

