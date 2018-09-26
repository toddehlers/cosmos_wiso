MODULE mo_physc1

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module mo_physc1 constants to communicate between the main program
  !                    and the radiation subroutines.
  !
  ! ----------------------------------------------------------------
  
  !   ratio atmospheric height/radius of the earth.
  REAL(dp), PARAMETER :: crae = 0.1277e-02_dp   
  REAL(dp) :: cdisse      ! solar constant normalised by its annual mean.
  REAL(dp) :: cdissem
  REAL(dp) :: czen1       ! orbital parameters.
  REAL(dp) :: czen2
  REAL(dp) :: czen3
  REAL(dp) :: czen1m
  REAL(dp) :: czen2m
  REAL(dp) :: czen3m
  REAL(dp) :: cosrad(128)  ! cos of longitudes used in solange.
  REAL(dp) :: sinrad(128)  ! sin of longitudes used in solange.

END MODULE mo_physc1
