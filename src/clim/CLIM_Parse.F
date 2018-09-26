      SUBROUTINE CLIM_Parse (cdarg, cdspawn, kargs, karmax)
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL T *
C               * -------------     ------- *
C               *****************************
C
C**** *parse*  - Parsing routine
C
C     Purpose:
C     -------
C     parse string cdstr into cdspawn. Blanks are the separators between args.
C     Each argument is defined as a continuous set of non-blanks characters
C     the first argument is put into cdspawn(1), the second in cdspawn(2) ...
C
C**   Interface:
C     ---------
C       *CALL*  *CLIM_Parse (cdarg, cdspawn, kargs, karmax)*
C
C     Input:
C     -----
C                cdarg   : line to be parsed (char string, including blanks)
C                karmax  : maximum number of arguments to be found (parameter)
C
C     Output:
C     ------
C                cdspawn : extracted character strings (possibly blank string)
C                karg    : number of arguments found (possibly 0)
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
C     See OASIS manual (2000) 
C 
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       1.0       J. Latour     00/09/15   created
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------------------- Include files ---------------------------
C
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
      USE mod_kinds_oasis
      USE mod_clim
      USE mod_comclim
C
C* ---------------------------- Argument declarations -------------------
C
      CHARACTER*80 cdarg
      INTEGER (kind=ip_intwp_p)      kargs, karmax
      CHARACTER*24 cdspawn(karmax)
C
C* ---------------------------- Local declarations -------------------
C
      CHARACTER*1  cdstr
      DIMENSION    cdstr(80)
      LOGICAL lnewarg
      INTEGER (kind=ip_intwp_p), PARAMETER :: maxchar=80 
      INTEGER (kind=ip_intwp_p), PARAMETER :: maxarg=24
      INTEGER (kind=ip_intwp_p) ii, ji, nulin
C
C* - Parameter declarations MUST match dimension declarations of cdstr
C* - and cdspawn
C
      CHARACTER (len=1), PARAMETER :: clblank = ' '
      CHARACTER (len=1), PARAMETER :: clcmt = '#'
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Initialize arrays
C        -----------------
C
C* - Transfer cdarg into cdstr
C
      DO 100 ji = 1,80
        cdstr(ji) = cdarg(ji:ji)
  100 CONTINUE
C
C* - Skip line if it is a comment
C     
      nulin=99
 110  IF (cdstr(1) .NE. clcmt) GO TO 130
      READ (UNIT = nulin, FMT = 1001) cdarg
      DO 120 ji = 1, 80
         cdstr(ji) = cdarg(ji:ji)
 120  CONTINUE
      GO TO 110
 130  CONTINUE 
 1001 FORMAT(A80)
C     
C* - Fill cdspawn with blanks and kargs with 0
C
      DO 210 ji = 1, karmax
        cdspawn(ji) = clblank
  210 CONTINUE 

      kargs = 0
      lnewarg = .TRUE.
C
C* - 2. Do the extraction work
C       ----------------------
C* - Check every character in cdstr, and fill cdspawn with contiguous
C* - strings of nonblank characters : the arguments.
C     
      DO 220 ji = 1, maxchar 
         IF ( cdstr(ji) .NE. clblank ) THEN
            IF ( lnewarg ) THEN
               kargs = kargs + 1 
               IF ( kargs .GT. karmax ) THEN
                  WRITE (UNIT = nulprt,FMT = *)
     $            ' Too many arguments in list passed to model'
                  WRITE (UNIT = nulprt,FMT = *)
     $            ' We STOP!!! Check the file namcouple'
                  CALL flush(nulprt)
                  CALL halte('STOP in CLIM_Parse.f')
               ENDIF
               ii = 1 
               cdspawn(kargs)(ii:ii) = cdstr(ji)
               lnewarg = .FALSE.
            ELSE
               ii = ii+1
               IF ( ii .GT. maxarg ) THEN
                  WRITE (UNIT = nulprt,FMT = *)
     $            ' Too many characters in argument passed to model'
                  WRITE (UNIT = nulprt,FMT = *)
     $            ' We STOP!!! Check the file namcouple'
                  CALL flush(nulprt)
                  CALL halte('STOP in CLIM_Parse.f')
               ENDIF
               cdspawn(kargs)(ii:ii) = cdstr(ji) 
            ENDIF
         ELSE
            lnewarg = .TRUE.
         ENDIF
  220 CONTINUE
C
C*    3. End of routine
C        --------------
C
#endif
      RETURN
      END
