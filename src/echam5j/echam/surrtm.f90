SUBROUTINE surrtm

  !=============================================================================
  !
  !- Description:
  !
  !   Initialisation of RRTM modules
  !
  !   M.A. Giorgetta, MPI, March 2000
  !
  !=============================================================================

  IMPLICIT NONE

  EXTERNAL surrtab,surrtpk,surrtrf,surrtftr,surrtbg2,surrta

  ! init RRTM modules
  ! -----------------
  CALL surrtab        ! mo_rrtab
  CALL surrtpk        ! mo_rrtwn
  CALL surrtrf        ! mo_rrtrf
  CALL surrtftr       ! mo_rrtftr
  CALL surrtbg2       ! mo_rrtbg2
  CALL surrta         ! mo_rrtaN (N=1:16) from file 

END SUBROUTINE surrtm
