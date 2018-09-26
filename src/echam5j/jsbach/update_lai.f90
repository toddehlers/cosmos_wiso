SUBROUTINE update_lai(kdim, ntiles, lai_clim, lai)

  USE mo_interpo, ONLY: wgt1, wgt2, nmw1, nmw2
  USE mo_radiation, ONLY: nmonth
  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  INTEGER, INTENT(in) :: kdim, ntiles
  REAL(dp), INTENT(in)  :: lai_clim(kdim,ntiles,0:13)
  REAL(dp), INTENT(out) :: lai(kdim,ntiles)

  IF (nmonth == 0) THEN
     ! Interpolation in time
     lai(:,:) = wgt1*lai_clim(:,:,nmw1) + wgt2*lai_clim(:,:,nmw2)
  ELSE
     ! Perpetual month
     lai(:,:) = lai_clim(:,:,nmonth)
  END IF

END SUBROUTINE update_lai
