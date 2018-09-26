MODULE mo_util_string
  !----------------------------------------------
  ! This module holds string conversion utilities
  !----------------------------------------------
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tolower        ! Conversion   : 'ABCXYZ' -> 'abcxyz'   
  PUBLIC :: toupper        ! Conversion   : 'abcxyz' -> 'ABCXYZ'
  PUBLIC :: char2          ! Conversion   : INTEGER  -> CHAR (LEN=2)
  PUBLIC :: separator      ! Format string: (/"-----...-----"/)
  !-----------------
  ! module variables
  !-----------------
#ifdef __PGI
  CHARACTER(len=*) ,PARAMETER :: separator = '(78("-"))'
#else
  CHARACTER(len=*) ,PARAMETER :: separator = '("'//REPEAT('-',78)//'")'
#endif
!==============================================================================
CONTAINS
!==============================================================================
  FUNCTION tolower (upper)
    !-----------------------------------
    ! Conversion: Uppercase -> Lowercase
    !-----------------------------------
    CHARACTER(LEN=*)              ,INTENT(in) :: upper
    CHARACTER(LEN=LEN_TRIM(upper))            :: tolower

    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('a')-ICHAR('A')

    DO i=1,LEN_TRIM(upper)
      IF (ICHAR(upper(i:i)) >= ICHAR('A') .AND. &
          ICHAR(upper(i:i)) <= ICHAR('Z')) THEN
        tolower(i:i) = CHAR( ICHAR(upper(i:i)) + idel )
      ELSE
        tolower(i:i) = upper(i:i)
      END IF
    END DO

  END FUNCTION tolower
!------------------------------------------------------------------------------
  FUNCTION toupper (lower)
    !-----------------------------------
    ! Conversion: Lowercase -> Uppercase
    !-----------------------------------
    CHARACTER(LEN=*)              ,INTENT(in) :: lower
    CHARACTER(LEN=LEN_TRIM(lower))            :: toupper

    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('A')-ICHAR('a')

    DO i=1,LEN_TRIM(lower)
      IF (ICHAR(lower(i:i)) >= ICHAR('a') .AND. &
          ICHAR(lower(i:i)) <= ICHAR('z')) THEN
        toupper(i:i) = CHAR( ICHAR(lower(i:i)) + idel )
      ELSE
        toupper(i:i) = lower(i:i)
      END IF
    END DO

  END FUNCTION toupper
!------------------------------------------------------------------------------
  FUNCTION char2 (i, zero)
    !----------------------------------------
    ! Conversion: INTEGER -> CHARACTER(LEN=2)
    !----------------------------------------
    CHARACTER(LEN=2)                       :: char2 ! result
    INTEGER          ,INTENT(in)           :: i     ! argument
    CHARACTER        ,INTENT(in) ,OPTIONAL :: zero  ! padding instead of '0'

    INTEGER ,PARAMETER :: i0 = ICHAR ('0')

    IF (i>99 .OR. i<0) THEN
      char2 = '**'
    ELSE
      char2(1:1) = CHAR(    i/10  + i0)
      char2(2:2) = CHAR(MOD(i,10) + i0)
    ENDIF

    IF(PRESENT(zero)) THEN
      IF(char2(1:1) == '0') char2(1:1) = zero
      IF(char2(2:2) == '0') char2(2:2) = zero
    ENDIF
  END FUNCTION char2
!------------------------------------------------------------------------------
END MODULE mo_util_string
