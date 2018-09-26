MODULE mo_forecast_switches

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_forecast_switches* switches related to the dynamics
  !                         and the general control of the forecast.
  !
  ! ----------------------------------------------------------------

  LOGICAL :: lsimdt   !  *true for semi implicit time scheme for
                      !  divergence,temperature and surface pressure
                      !  equations.
  LOGICAL :: lsimzq   !  *true for semi implicit time scheme for
                      !  vorticity and humidity equations.
  LOGICAL :: lvtmpc1  !  *true for virtual temperature.
                      !  i.e:*rv*.ne.*rd*.
  LOGICAL :: lvtmpc2  !  *true for influence of humidity on *cp*.
                      !  i.e:*cpv*.ne.*cpd*.
  LOGICAL :: lumax    !  *true to compute and print information on
                      !  maximum wind.
  LOGICAL :: lzondia

END MODULE mo_forecast_switches
