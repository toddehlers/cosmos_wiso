MODULE mo_field

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  REAL(dp), ALLOCATABLE :: field1(:,:,:)               ! (nlon,numfl1,ngl)
  REAL(dp), ALLOCATABLE :: field2(:,:,:)               ! (nlon,numfl2,ngl)

END MODULE mo_field
