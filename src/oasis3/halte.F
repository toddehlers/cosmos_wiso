      SUBROUTINE halte (cdtext)
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL C *
C               * -------------     ------- *
C               *****************************
C
C**** *halte*  - Abort the program
C
C     Purpose:
C     -------
C     Print an error message to standard output and abort the coupler
C
C**   Interface:
C     ---------
C       *CALL*  *halte (cdtext)*
C
C     Input:
C     -----
C                cdtext   : character string to be printed
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
C     See OASIS 2.2 manual (1997) 
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       2.2       S. Valcke      97/11/18  created
C       2.3       L. Terray      99/02/24  modified: X format for NEC
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------------------- Include files ---------------------------
C
      USE mod_unit
      USE mod_hardware
#if defined use_comm_MPI1 || defined use_comm_MPI2 || defined use_comm_GSIP || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_GSIP && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
      USE mod_comclim
#endif
C
C* ---------------------------- Argument declarations ----------------------
C
      CHARACTER*(*) cdtext
C
C* ---------------------------- Local declarations ----------------------
C
      CHARACTER cpbase
      CHARACTER*10 cprpt, cpdots
      CHARACTER*69 cline
      PARAMETER ( cpbase = '-' )
      PARAMETER ( cprpt = '* ===>>> :' )
      PARAMETER ( cpdots = '  ------  ' )
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Print text 
C        ----------
C
          cline = ' '
          ilen = len(cdtext)
          DO 110 jl = 1, ilen
            cline(jl:jl) = cpbase
 110      CONTINUE
          WRITE(UNIT = nulou,FMT='(/,A,1X,A)') cpdots, cline
          WRITE(UNIT = nulou,FMT='(/,A,1X,A,/)') cprpt, cdtext
          WRITE(UNIT = nulou,FMT='(A,1X,A,/)')cpdots, cline
C
C
C*    2. Close the coupler output
C        ------------------------
C
          CLOSE(nulou)
C
C
C*    3. Abort the coupler
C        -----------------
C
#if defined use_comm_MPI1 || defined use_comm_MPI2
      IF (cchan .eq. 'MPI2' .or. cchan .EQ. 'MPI1' ) THEN
          CALL MPI_ABORT (mpi_comm, 0, ierror)
          IF (ierror /= CLIM_Ok) THEN
              WRITE (*,'(a)') ' MPI_ABORT failed'
              WRITE (*,'(a,i4)') ' Error =  ', ierror
              CALL abort
          ENDIF
      ENDIF
#endif
#if defined use_comm_GSIP || defined use_comm_SIPC || defined use_comm_GMEM || defined use_comm_PIPE || defined use_comm_NONE
      CALL abort
#endif
      RETURN
      END
