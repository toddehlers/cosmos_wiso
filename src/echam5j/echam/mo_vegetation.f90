MODULE mo_vegetation
  
  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  ! _______________________________________________________________________
  !
  ! module *mo_vegetation* variables for computation of evapotranspiration
  !
  !
  !      Version 1  C.A Blondin  2/12/86  ECMWF
  !
  ! _______________________________________________________________________

  REAL(dp) :: cva     !   constant to define the stomatal resistance
  REAL(dp) :: cvb     !   constant to define the stomatal resistance
  REAL(dp) :: cvc     !   minimum stomatal resistance
  REAL(dp) :: cvbc    !   cvb*cvc
  REAL(dp) :: cvxpklt !   exp(cvklt)
  REAL(dp) :: cvxmklt !   exp(-cvklt)
  REAL(dp) :: cvkc    !   cvk*cvc
  REAL(dp) :: cvabc   !   (cva+cvbc)/cvc
  REAL(dp) :: cvroots !   percentage of roots in the 1st soil layer
  REAL(dp) :: cvrootd !   percentage of roots in the 2nd soil layer
  REAL(dp) :: cvrootc !   percentage of roots in the 3rd soil layer
  REAL(dp) :: cvrad   !   fraction of the net s.w radiation contributing to p.a.r
  REAL(dp) :: cvinter !   efficency of interception of precipitation as rain
  REAL(dp) :: cvk

END MODULE mo_vegetation
