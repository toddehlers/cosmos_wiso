SUBROUTINE skip (cd_one, id_len)
!
!**** SKIP
!
!     Purpose:
!       Skip line if it is a comment
!
!     Interface: 
!       Call skip (cl_one)
!
!     Method:
!       Read the first caracter of the line and skip line if 
!       it is a comment
!
!     External:
!       none
!
!     Files:
!       none
!   
!     References:
!
!     History:
!     --------
!       Version   Programmer     Date        Description
!       ------------------------------------------------
!       2.5       A.Caubel       2002/04/04  created
!
!*-----------------------------------------------------------------------
!
!** + DECLARATIONS
!
!
!** ++ Include files
!
  USE mod_kinds_oasis
  USE mod_unit
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: ib,id_len
  CHARACTER(len=80) :: cl_line
  CHARACTER(len=1), DIMENSION(id_len) :: cd_one
  CHARACTER(len=1) :: cl_two
!
!*-----------------------------------------------------------------------
!
  cl_two='#'
100 IF (cd_one(1) .NE. cl_two) GO TO 120
  READ (UNIT = nulin, FMT = 1001) cl_line
  DO ib = 1,id_len
    cd_one(ib) = cl_line(ib:ib)
  END DO
  GO TO 100
120 CONTINUE 
1001 FORMAT(A80)
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE skip
!
!*========================================================================
