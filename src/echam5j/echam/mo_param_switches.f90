MODULE mo_param_switches

  IMPLICIT NONE

  ! M.A. Giorgetta, March 2000, lrad added
  !
  ! ----------------------------------------------------------------
  !
  ! module *mo_param_switches* switches related to the parameterisations of
  !                            diabatic processes except for radiation.
  !
  ! ----------------------------------------------------------------

  LOGICAL :: lphys    !   *true for parameterisation of diabatic processes.
  LOGICAL :: lrad     !   *true for radiation.
  LOGICAL :: lvdiff   !   *true for vertical diffusion.
  LOGICAL :: lcond    !   *true for large scale condensation scheme.
  LOGICAL :: lsurf    !   *true for surface exchanges.
  LOGICAL :: lcover   !   *true for prognostic cloud cover scheme
  LOGICAL :: lconv    !   *true to allow convection
  INTEGER :: iconv    !   *1,2,3 for different convection schemes
  LOGICAL :: lgwdrag  !   *true for gravity wave drag scheme
  LOGICAL :: lice     !   *true for sea-ice temperature calculation

END MODULE mo_param_switches
