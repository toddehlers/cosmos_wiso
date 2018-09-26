SUBROUTINE update_cover_fract(kdim, ntiles, ntiles_lct, veg_ratio, fract)

  USE mo_interpo, ONLY: wgt1, wgt2, nmw1, nmw2
  USE mo_radiation, ONLY: nmonth
  USE mo_exception, ONLY: finish, message
  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  INTEGER, INTENT(in) :: kdim, ntiles, ntiles_lct
  REAL(dp), INTENT(in)  :: veg_ratio(kdim,0:13)
  REAL(dp), INTENT(inout) :: fract(kdim,ntiles)

  INTEGER :: i

!  IF (ntiles_lct /= 3) CALL message('update_cover_fract','This only works for 3 tiles (ECHAM compatibility')
  IF (ntiles_lct /= 2) CALL message('update_cover_fract','Watch out for jsbach_interface and the correct usage of cover_fract !')
  IF (ntiles_lct /= 2) CALL message('update_cover_fract','REWRITTEN FOR 2TILE')

  IF (nmonth == 0) THEN
     ! Interpolation in time
     fract(:,1) = wgt1*veg_ratio(:,nmw1) + wgt2*veg_ratio(:,nmw2)
  ELSE
     ! Perpetual month
     fract(:,1) = veg_ratio(:,nmonth)
  END IF
  IF (ntiles_lct == 3) THEN
     WHERE (fract(:,3) >= 1.) fract(:,1) = 0._dp
     fract(:,2) = MAX(0._dp, 1._dp - fract(:,1) - fract(:,3))
     
     IF (ntiles /= ntiles_lct) THEN
        DO i=1,ntiles/ntiles_lct
           fract(:,ntiles_lct*(i-1)+1:ntiles_lct*i) = fract(:,1:ntiles_lct)
        END DO
     END IF
  ENDIF

END SUBROUTINE update_cover_fract
