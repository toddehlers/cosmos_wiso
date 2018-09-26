MODULE mo_hdiff

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  ! ---------------------------------------------------------------
  !
  ! module *mo_hdiff* - coefficients for horizontal diffusion.
  !
  ! ---------------------------------------------------------------

  REAL(dp)    :: damhih      !  for extra diffusion in middle atmosphere
  REAL(dp)    :: dampth
  REAL(dp)    :: difvo       !  coefficient for vorticity.
  REAL(dp)    :: difd        !  coefficient for divergence.
  REAL(dp)    :: dift        !  coefficient for temperature.
  REAL(dp), ALLOCATABLE    :: diftcor(:) !  correction profile for temperature
  REAL(dp)    :: enstdif     !  *factor by which stratospheric
                             !  horizontal diffusion is increased from one
                             !  level to next level above.
  INTEGER :: nlvstd1         !  *last (uppermost) layer at which
                             !  stratospheric horizontal diffusion is
                             !  enhanced.
  INTEGER :: nlvstd2         !  *first (lowest) layer at which
                             !  stratospheric horizontal diffusion is
                             !  enhanced.
  LOGICAL :: ldiahdf         !  .true. for statistics of horizontal diffusion

END MODULE mo_hdiff
