MODULE mod_unit
!
! -- unit.h   23-08-95   Version 2.0   Author: Laurent Terray
!    ******   20-06-97   Version 2.2   Mods: add nudum (S. Valcke)
!             30-03-99   Version 2.3   Mods: add nulgn (S. Valcke)
!             23-03-01   Version 2.5   Mods: add lncdfxxx (S. Valcke)
!@
!@  Contents : unit numbers
!@  --------
!@
!@ -- nulin : logical unit for coupler input
!@
!@ -- nulou : logical unit for coupler output
!@
!@ -- nulgr : logical unit for gcm's grids
!@
!@ -- nulma : logical unit for gcm's masks
!@
!@ -- nulsa : logical unit for gcm's surfaces
!@
!@ -- nultr : trace file for CLIM and PVM messages
!@
!@ -- nulcc : file for ANAISM data
!@
!@ -- nulgg : file for ANAISG data
!@
!@ -- nulgn : file for NINENN weights and addresses
!@
!@ -- nulan : file for ANAIS(M-G) specific output
!@
!@ -- nulrd : file for reduced grid masks
!@
!@ -- nudum : dummy file signaling that SVIPC pools for exchange of 
!@            initial infos are OPENED
!@
! -------------------------------------------------------------------   
!
  USE mod_kinds_oasis
  INTEGER(kind=ip_intwp_p) ::  nulin, nulou, nulgr, nulma, nulsu, nultr, nulcc 
  INTEGER(kind=ip_intwp_p) ::  nulgg, nulan, nulrd, nudum, nulgn 
!
END MODULE mod_unit
! -------------------------------------------------------------------   











