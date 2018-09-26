!-------------------------------------------------------------------
! BOP
!
! !MODULE: mod_kinds_oasis
MODULE mod_kinds_oasis
!
! !USES:
!
! !PUBLIC TYPES:
IMPLICIT NONE
SAVE
!
! !PUBLIC MEMBER FUNCTIONS:
!
! !PUBLIC DATA MEMBERS:
!
! !PARAMETERS:
  INTEGER, PARAMETER :: ip_single_p = SELECTED_REAL_KIND(6,37)
  INTEGER, PARAMETER :: ip_double_p = SELECTED_REAL_KIND(12,307)
#ifdef use_realtype_single
  INTEGER, PARAMETER :: ip_realwp_p = ip_single_p
  LOGICAL, PARAMETER :: ll_single = .TRUE.
#else
  INTEGER, PARAMETER :: ip_realwp_p = ip_double_p
  LOGICAL, PARAMETER :: ll_single = .FALSE. 
#endif
  INTEGER, PARAMETER :: ip_i2_p = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: ip_i4_p = SELECTED_INT_KIND(9)
#ifdef SX
  INTEGER, PARAMETER :: ip_i8_p = SELECTED_INT_KIND(15)
#else
  INTEGER, PARAMETER :: ip_i8_p = SELECTED_INT_KIND(18)
#endif
  INTEGER, PARAMETER :: ip_intwp_p = ip_i4_p 
!
! !DESCRIPTION:
! This modules contains the parameters defining the precision used for 
! real and integer variables
!   
  END MODULE mod_kinds_oasis

