      MODULE MO_COMMO3
      USE MO_PARAM1
!===> COMMO3.h
      LOGICAL, POINTER :: LTLEV(:,:,:), LSLEV(:,:,:)
      REAL, POINTER :: SLEVI(:,:,:), TLEVI(:,:,:)
!<===END COMMO3

      CONTAINS

      SUBROUTINE alloc_mem_commo3

      ALLOCATE(SLEVI(IE,JE,KE),TLEVI(IE,JE,KE))
      ALLOCATE(LSLEV(IE,JE,KE),LTLEV(IE,JE,KE))

      END SUBROUTINE alloc_mem_commo3
      END MODULE MO_COMMO3
