      MODULE MO_COMMO2
      USE MO_PARAM1

      IMPLICIT NONE

!==> COMMO2
#ifndef __coupled
!UWE***FORCING FIELDS
      REAL, POINTER :: TXO1(:,:),TXO2(:,:),TYE1(:,:),                   &
     &   TYE2(:,:),FCLOU1(:,:),FSWR1(:,:),FU101(:,:),                   &
     &   FPREC1(:,:),FTDEW1(:,:),FSLP1(:,:),                            &
     &   FCLOU2(:,:),FSWR2(:,:),FU102(:,:),                             &
     &   FPREC2(:,:),FTDEW2(:,:),                                       &
     &   TAFO1(:,:),TAFO2(:,:),FSLP2(:,:)
#endif /*__coupled*/
!<==END COMMO2

      CONTAINS

      SUBROUTINE alloc_mem_commo2

#ifndef __coupled
      ALLOCATE(TXO1(IE,JE),TXO2(IE,JE),TYE1(IE,JE),                     &
     &   TYE2(IE,JE),FCLOU1(IE,JE),FSWR1(IE,JE),FU101(IE,JE),           &
     &   FPREC1(IE,JE),FTDEW1(IE,JE),FSLP1(IE,JE),                      &
     &   FCLOU2(IE,JE),FSWR2(IE,JE),FU102(IE,JE),                       &
     &   FPREC2(IE,JE),FTDEW2(IE,JE),                                   &
     &   TAFO1(IE,JE),TAFO2(IE,JE),FSLP2(IE,JE))
#endif /*__coupled*/
      END SUBROUTINE alloc_mem_commo2

      END MODULE MO_COMMO2
