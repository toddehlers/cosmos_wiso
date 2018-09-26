MODULE mo_base_types
  USE mo_kind,     ONLY:dp
  TYPE base_pointer
     REAL(dp), POINTER :: vp(:)
     CHARACTER (8) :: cp
     INTEGER :: np
  END TYPE base_pointer

  TYPE pointer_set
     TYPE (base_pointer), POINTER :: set(:)
     INTEGER :: nsets
  END TYPE pointer_set

END MODULE mo_base_types
