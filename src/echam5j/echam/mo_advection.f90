MODULE mo_advection

  ! L. Kornblueh, MPI, August 2001, initial coding

  IMPLICIT NONE

  PUBLIC

  ! defined advection schemes for passive tracer (in advection sense)

  INTEGER, PARAMETER :: no_advection    = 0
  INTEGER, PARAMETER :: semi_lagrangian = 1
  INTEGER, PARAMETER :: spitfire        = 2
  INTEGER, PARAMETER :: tpcore          = 3

  ! select variable, set in namelist runctl 

  INTEGER :: iadvec

  ! general dimension parameters to be set in the init section of
  ! the diferent advection modules 

  INTEGER :: nalatd, nalond, nalev, nacnst, nalat, nalon

END MODULE mo_advection
