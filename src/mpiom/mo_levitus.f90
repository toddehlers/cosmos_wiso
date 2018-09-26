MODULE mo_levitus
  IMPLICIT NONE
  LOGICAL, ALLOCATABLE ::  LTLEV(:,:,:), LSLEV(:,:,:)
  REAL, ALLOCATABLE ::  SLEVI(:,:,:), TLEVI(:,:,:)
CONTAINS
  SUBROUTINE init_levitus (ie, je, ke)
    INTEGER, INTENT(in) :: ie, je, ke
    ALLOCATE(SLEVI(IE,JE,KE),TLEVI(IE,JE,KE))
    ALLOCATE(LSLEV(IE,JE,KE),LTLEV(IE,JE,KE))
    SLEVI(:,:,:)=0.0
    TLEVI(:,:,:)=0.0
    LSLEV(:,:,:)=.FALSE.
    LTLEV(:,:,:)=.FALSE.
END SUBROUTINE init_levitus
  SUBROUTINE free_levitus
   DEALLOCATE(SLEVI,TLEVI)
   DEALLOCATE(LSLEV,LTLEV)
END SUBROUTINE free_levitus
END MODULE mo_levitus
