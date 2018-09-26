      SUBROUTINE chkpar
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 0 *
C               * -------------     ------- *
C               *****************************
C
C**** *chkpar*  - Parameter checking
C
C     Purpose:
C     -------
C     Checks option compatibility between oasis and remote models
C     as well as some basic dimension checks
C
C**   Interface:
C     ---------
C       *CALL*  *chkpar*
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
C
C     Externals:
C     ---------
C     imaxim
C
C     Reference:
C     ---------
C     See OASIS manual (1995) 
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       2.0       L. Terray      95/10/01  created
C       2.2       L. Terray      97/11/13  added: check nitfn if cchan
C                                                 equals NONE
C       2.3       S. Valcke      99/04/30  added: printing levels
C       2.3       S. Valcke      99/10/12  check if nninnfl and nninnflg
C                                          are different
C       2.4       S. Valcke      2K/02/03  check CLIM_MaxPOrt and jpfield
C       2.4       S. Valcke      2K/02/03  check CLIM_MaxSegment and jpparal
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* --------------- Include files and use of modules---------------------------
C
      USE mod_kinds_oasis
      USE mod_parameter
      USE mod_string
      USE mod_extrapol
      USE mod_experiment
      USE mod_timestep
      USE mod_unit
      USE mod_hardware
      USE mod_printing
C
C* ---------------------------- Local declarations ----------------------
C
      INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: itime, 
     $    iexch
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Allocations and Initializations
C        -------------------------------
C
      IF (nlogprt .GE. 1) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '           ROUTINE chkpar  -  Level 0'
          WRITE (UNIT = nulou,FMT = *) 
     $    '           **************     *******'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' control of run time options'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
      ENDIF
C
C* In the CLIM CASE, mstep, mfcpl and mdt are not defined anymore in inicmc.f.
C
      IF (cchan.ne.'MPI1'.and.cchan.ne.'MPI2'.and.cchan.ne.'GSIP') then
C
C* Array allocation
C
          ALLOCATE (itime(ig_nmodel))
          ALLOCATE (iexch(ig_nmodel))
C
C* Initialization of remote models time variables.
C
          DO 110 jm = 1, ig_nmodel
            itime(jm) = mstep(jm) * mdt(jm)
            iexch(jm) = mfcpl(jm) * mdt(jm)
 110      CONTINUE

C
C*    2. Basic checks
C        ------------
C
C*    3. Time compatibility checks
C        -------------------------
C
          DO 310 jm = 1, ig_nmodel
            WRITE (UNIT = nulou,FMT = *)'itime ',itime(jm)
            WRITE (UNIT = nulou,FMT = *)'ntime ',ntime
            IF (itime(jm) .NE. ntime) THEN
                WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
                WRITE (UNIT = nulou,FMT = *) 
     $              ' ===>>> : total runtime for model',jm
                WRITE (UNIT = nulou,FMT = *)
     $              '          incompatible with coupler total time '
                WRITE (UNIT = nulou,FMT = *) 
     $              ' ======         =======              '
                WRITE (UNIT = nulou,FMT = *) ' '
                WRITE (UNIT = nulou,FMT = *) 
     $              ' The current simulation uses asynchroneous 
     $              coupling'
                WRITE (UNIT = nulou,FMT = *) 
     $              '                             -------------'
                WRITE (UNIT = nulou,FMT = *) ' '
            ENDIF
C
C* The exchange frequency data sent from remote models is the
C  smallest frequency used in a given model. This is a crude 
C  check of the compatibility of coupling options.
C
            WRITE (UNIT = nulou,FMT = *)'iexch ',iexch(jm)
            WRITE (UNIT = nulou,FMT = *)'nstep ',nstep
            IF (MOD(iexch(jm),nstep) .NE. 0) THEN
                WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
                WRITE (UNIT = nulou,FMT = *) 
     $              ' ===>>> : exchange time incompatible with model',
     $              jm
                WRITE (UNIT = nulou,FMT = *) 
     $              ' ======            =====          '
                WRITE (UNIT = nulou,FMT = *) ' '
                WRITE (UNIT = nulou,FMT = *) 
     $              ' We STOP        !!! rerun with new  parameters'
                WRITE (UNIT = nulou,FMT = *) 
     $              '                               ---------------'
                WRITE (UNIT = nulou,FMT = *) ' '
                CALL HALTE ('STOP in chkpar')
            ENDIF
 310      CONTINUE 
          DEALLOCATE (itime)
          DEALLOCATE (iexch)
      ENDIF
C
C
C*    4. End of routine
C        --------------
C
      IF (nlogprt .GE. 1) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '          ---------- End of routine chkpar --------'
          CALL FLUSH (nulou)
      ENDIF
      RETURN
      END


