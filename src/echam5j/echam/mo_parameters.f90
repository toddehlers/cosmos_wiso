MODULE mo_parameters

  IMPLICIT NONE

  ! parameters controlling array sizes.

  INTEGER, PARAMETER :: jpm     = 106 ! max zonal wave number
  INTEGER, PARAMETER :: jpn     = jpm ! max meridional wave number for m=0
  INTEGER, PARAMETER :: jpk     = jpm ! max meridional wave number
  INTEGER, PARAMETER :: jpnlev  = 999 ! number of vertical levels.
                                      ! for ntrn only, f90 namelist restriction
  INTEGER, PARAMETER :: jpgrnd  = 5   ! number of soil layers
                                      
  INTEGER, PARAMETER :: jptrac  = 100 ! maximum number of prognostic tracers
  INTEGER, PARAMETER :: jps     = 3   ! basic Spitfire variables without tracers
                                      ! i.e. humidity, cloud water and cloud ice
  INTEGER, PARAMETER :: jpmp1=jpm+1

END MODULE mo_parameters
