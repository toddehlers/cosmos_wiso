      SUBROUTINE initim
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 0 *
C               * -------------     ------- *
C               *****************************
C
C**** *initim*  - Initialize time
C
C     Purpose:
C     -------
C     Initialize date at run starting date and find coupler timestep
C     value
C
C**   Interface:
C     ---------
C       *CALL*  *initim*
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
C       1.0       L. Terray      94/01/01  created
C       2.0beta   L. Terray      95/08/30  modified : new structure
C       2.0       L. Terray      96/02/01  modified : get final iteration
C       2.3       S. Valcke      99/04/30  added: printing levels
C       2.5       S. Valcke      01/12/10  Explicit intialization of nddeb
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------- Include files and USE of modules---------------------------
C
      USE mod_parameter
      USE mod_string
      USE mod_timestep 
      USE mod_unit
      USE mod_calendar
      USE mod_printing
C
C* ---------------------------- Local functions -------------------------
C
C* - Time functions
C    Date is under form: ccaammdd
C            ndd   --->>> extract dd from ccaammdd
C            nmm   --->>> extract mm from ccaammdd
C            nccaa --->>> extract ccaa from ccaammdd
C            ncth  --->>> turn seconds into hours
C
      ndd(kidat) = MOD(kidat,100)
      nmm(kidat) = MOD((kidat - ndd(kidat))/100,100)
      nccaa(kidat) = kidat / 10000
      ncth(ksec) = ksec / 3600
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Initializations
C         ---------------
C     
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) 
     $           '           ROUTINE initim  -  Level 0'
            WRITE (UNIT = nulou,FMT = *) 
     $           '           **************     *******'
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) 
     $           ' Determine time iteration parameters '
            WRITE (UNIT = nulou,FMT = *) 
     $           ' and initialize simulation date '
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
C     
C* First, get initial date of this run
C
         ig_date(:) = 0
         njini = ndd(ndate)
         ig_date(3) = njini
         nmini = nmm(ndate)
         ig_date(2) = nmini
         naini = nccaa(ndate)
         ig_date(1) = naini
C     
C* Print initial date of this run
C
         WRITE (UNIT = nulou,FMT = 1001) njini,  nmini, naini
         WRITE (UNIT = nulou,FMT = *) ' '
C
C* Following lines are called only if one field (at least) goes through Oasis
C
         IF (lg_oasis_field) THEN
C     
C* Get time iteration parameters
C     
C* Get minimum exchange frequency
C     
            ifrqmin = iminim (nfexch, ig_nfield)
C     
C* Get greatest common divisor for all frequencies
C     
            nstep = idivmax (nfexch, ig_nfield, ifrqmin)
C     
C* Get total number of iterations
C     
            niter = ntime / nstep
C     
C* Get final iteration number
C     
            nitfn = niter - 1
C     
         IF (nlogprt .GE. 1) THEN
C*    Print number of iterations and timestep
C     
            WRITE (UNIT = nulou,FMT = *) 
     $           '      Number of iterations = ', niter
            WRITE (UNIT = nulou,FMT = *) '     ====================   '
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) 
     $           '      Value of timestep = ', nstep
            WRITE (UNIT = nulou,FMT = *) '     =================   '
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
C     
C* Get beginning date of the whole coupled simulation
C     
         nddeb = 0
         njdeb = ndd(nddeb)
         nmdeb = nmm(nddeb)
         nadeb = nccaa(nddeb)
      ENDIF
C     
C* Formats
C     
 1001    FORMAT(5X,'  Run initial date : ',I2,' - ',I2,' - ',I4,
     $        /,5X,'  ================ ')
C     
C     
C*    2. End of routine
C     --------------
C     
         IF (nlogprt .GE. 1) THEN 
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) 
     $           '          --------- End of routine initim ---------'
            CALL FLUSH (nulou)
         ENDIF
      RETURN
      END
