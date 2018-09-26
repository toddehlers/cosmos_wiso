MODULE mo_post

  USE mo_kind,        ONLY: dp
  USE mo_parameters

  IMPLICIT NONE

  !
  ! module *mo_post* - control variables for postprocessing
  !
  !      e. kirk     uni hamburg     2-mar-89
  !

  REAL(dp) :: setval(256)

  INTEGER, ALLOCATABLE :: npplev(:,:)  !   selected levels
  INTEGER :: npbits(256)
  INTEGER :: nlevg3(256)

  INTEGER :: nlalo                !   nlon * ngl
  INTEGER :: nung4x
  INTEGER :: ncdmin
  INTEGER :: ncdmax

  LOGICAL :: laccu(256)
  LOGICAL :: lppspe
  LOGICAL :: lppd
  LOGICAL :: lppvo
  LOGICAL :: lppt
  LOGICAL :: lppp
  LOGICAL :: lppq
  LOGICAL :: lppxl
  LOGICAL :: lppxi

  CHARACTER (8) ::  yn(256)         !   array names

END MODULE mo_post
