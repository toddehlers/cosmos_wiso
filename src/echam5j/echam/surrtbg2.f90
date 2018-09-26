SUBROUTINE surrtbg2

  !=============================================================================
  !
  !- Description:
  !
  !   Calculate lookup tables in mo_rrtbg2 for functions needed in routine taumol
  !
  !   M.A. Giorgetta, MPI, June 2000
  !
  !=============================================================================

  USE MO_KIND  , ONLY: dp

  ! Output
  USE mo_rrtbg2, ONLY: corr1, corr2

  IMPLICIT NONE

  INTEGER :: i
  REAL(dp):: fp, rtfp

  corr1(0)   = 1._dp
  corr1(200) = 1._dp
  corr2(0)   = 1._dp
  corr2(200) = 1._dp

  DO i = 1,199
     fp       = 0.005_dp*REAL(i,dp)
     rtfp     = SQRT(fp)
     corr1(i) = rtfp/fp
     corr2(i) = (1._dp-rtfp)/(1._dp-fp)
  ENDDO

!CDIR DU_UPDATE(CORR1)
!CDIR DU_UPDATE(CORR2)

END SUBROUTINE SURRTBG2
