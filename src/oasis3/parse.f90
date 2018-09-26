SUBROUTINE parse (cdone, cdtwo, knumb, klen, kleng)
  USE mod_kinds_oasis
!****
!               *****************************
!               * OASIS ROUTINE  -  LEVEL T *
!               * -------------     ------- *
!               *****************************
!
!**** *parse*  - Parsing routine
!
!     Purpose:
!     -------
!     Find the knumb'th string in cdone and put it in cdtwo.
!     A string is defined as a continuous set of non-blanks characters
!
!**   Interface:
!     ---------
!       *CALL*  *parse (cdone, cdtwo, knumb, klen, kleng)*
!
!     Input:
!     -----
!                cdone : line to be parsed (char string)
!                knumb : rank within the line of the extracted string (integer)
!                klen  : length of the input line (integer)
!
!     Output:
!     ------
!                cdtwo : extracted character string (char string)
!                kleng : length of the extracted string (integer)
!
!     Workspace:
!     ---------
!     None
!
!     Externals:
!     ---------
!
!     Reference:
!     ---------
!     See OASIS manual (1995)
!
!     History:
!     -------
!       Version   Programmer     Date      Description
!       -------   ----------     ----      -----------  
!       2.0       L. Terray      95/09/01  created
!                 O. Marti     2000/11/08  simplify by using F90 
!                                          CHARACTER functions
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!* ---------------------------- Include files ---------------------------
!
  USE mod_unit
!
!* ---------------------------- Argument declarations -------------------
!
  INTEGER (kind=ip_intwp_p), INTENT ( in) :: knumb, klen
  CHARACTER (len=klen), INTENT ( inout) :: cdone 
  CHARACTER (len=klen), INTENT ( out) :: cdtwo
  INTEGER (kind=ip_intwp_p), INTENT ( out) :: kleng
!
!* ---------------------------- Local declarations -------------------
!
  CHARACTER (len=klen) :: clline
  CHARACTER (len=klen) :: clwork
  CHARACTER (len=1), SAVE :: clblank = ' ', clcmt = '#'
!
!* ---------------------------- Poema verses ----------------------------
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!*    1. Skip line if it is a comment
!        ----------------------------
!
100 IF (cdone(1:1) .NE. clcmt) GO TO 120
  READ (UNIT = nulin, FMT = 1001) clline 
  cdone(1:klen) = clline(1:klen)
  GO TO 100
120 CONTINUE 
1001 FORMAT(A80)
!
!
!*    2. Do the extraction job
!        ---------------------
!
!* - Fill cdtwo with blanks
!
  cdtwo = clblank
!
!* Fill temporary string and remove leading blanks
!
  clwork = ADJUSTL ( cdone)
!
!* - If there are no more characters, kleng=-1
!
  IF ( LEN_TRIM ( clwork) .LE. 0) THEN
      kleng = -1
      RETURN
  END IF
!
!* - If this is the one we're looking for, skip
!    otherwise go knumb-1 more sets of characters
!
  IF (knumb .GE. 2) THEN
      DO jl = 1, knumb-1
        ii = INDEX ( clwork, clblank) - 1
        clwork ( 1:ii) = clblank
        clwork = ADJUSTL ( clwork)
!
!* - If there are no more characters, kleng=-1
!
        IF (LEN_TRIM ( clwork) .LE. 0) THEN
            kleng = -1
            RETURN
        END IF
      END DO
  END IF
!
!* - Find the length of this set of characters
!
  kleng = INDEX ( clwork, clblank) - 1
!
!* - Copy to cdtwo
!
  cdtwo ( 1:kleng) = clwork ( 1: kleng)
!
!
!*    3. End of routine
!        --------------
!
  RETURN
END SUBROUTINE parse




