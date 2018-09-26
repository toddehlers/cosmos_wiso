MODULE mo_gl1

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

! global buffered variables

  REAL(dp), TARGET, ALLOCATABLE :: q(:,:)
  REAL(dp), TARGET, ALLOCATABLE :: x(:,:)
  REAL(dp), TARGET, ALLOCATABLE :: xt(:,:,:)

END MODULE mo_gl1
