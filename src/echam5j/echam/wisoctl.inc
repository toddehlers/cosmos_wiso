!* Namelist WISOCTL for control parameters for water isotopes
!
! name        type     purpose
! ----        ----     -------
! LWISO       LOGICAL  switch water isotope diagnostics on/off
! LWISO_RERUN LOGICAL  switch for initial/resumed water isotope values
! NWISO       INTEGER  number of water isotope tracers
! NWISOTYP    INTEGER  type of water isotope tracer (1=H2-16O, 2=H2-18O, 3=HDO)
!
NAMELIST /WISOCTL/ LWISO, LWISO_RERUN, NWISO, NWISOTYP

