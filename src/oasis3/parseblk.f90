SUBROUTINE parseblk (cdone, cdtwo, knumb, klen, kleng)
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
!     Get the rest of the line starting at the knumb'th string.
!     A string is defined as a continuous set of non-blanks characters
!
!**   Interface:
!     ---------
!       *CALL*  *parseblk (cdone, cdtwo, knumb, klen, kleng)*
!
!     Input:
!     -----
!                cdone : line to be parsed (char string)
!                knumb : rank within the line of the starting string (integer)
!                klen  : length of the input line (integer)
!
!     Output:
!     ------
!                cdtwo : extracted rest of line, including blanks (char string)
!                kleng : length of the extracted string (integer)
!
!     Workspace:
!     ---------
!     None
!
!     Externals:
!     ---------
!
!     History:
!     -------
!       Version   Programmer     Date      Description
!       -------   ----------     ----      -----------  
!       2.5       S. Valcke      00/09/08  Adapted from parse.f
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
  INTEGER (kind=ip_intwp_p) :: il, kleng_aux
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
  il = INDEX ( cdone, clblank)
  kleng_aux = 1
  IF (INDEX ( cdone, clblank).EQ.1) THEN
      DO WHILE (cdone(il+1:il+1).EQ.clblank)
        kleng_aux = kleng_aux +1
        il = il+1
        IF (il+1.GT.klen) GO TO 130
      ENDDO
  ENDIF
130 CONTINUE
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
        il = ii + 1 
        DO WHILE (clwork(il:il).EQ.clblank)
          kleng_aux = kleng_aux +1
          il = il + 1
          IF (il.GT.klen) GO TO 140
        ENDDO
140 CONTINUE
        kleng_aux = kleng_aux + ii
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
!* - Find the length of the rest of the line
!
  kleng = klen - kleng_aux
!
!* - Copy to cdtwo
!
  cdtwo ( 1:kleng) = clwork ( 1: kleng)
!
!*    3. End of routine
!        --------------
!
  RETURN
END SUBROUTINE parseblk

