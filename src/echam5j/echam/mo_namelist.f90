MODULE mo_namelist

! position_nml - position namelist file for reading
!
! Author:
!
! L. Kornblueh, MPI, March 2001, original source
! L. Kornblueh, MPI, January 2003, added error message for OPEN
! L. Kornblueh, MPI, April 2005, added close_nml
!

  USE mo_util_string,   ONLY: tolower
  USE mo_exception,     ONLY: message, finish, message_text

  IMPLICIT NONE

  ! Public entities

  PRIVATE
  PUBLIC :: position_nml           ! function  : position namelist file
  PUBLIC :: open_nml               ! subroutine: open namelist file
  PUBLIC :: close_nml              ! subroutine: close namelist file
  PUBLIC :: nnml                   ! default namelist unit
  PUBLIC :: POSITIONED, MISSING,  &! return values from position_nml
            LENGTH_ERROR, READ_ERROR
  
  ! return values of function 'position_nml':

  INTEGER, PARAMETER :: POSITIONED   =  0 ! file pointer set to namelist group
  INTEGER, PARAMETER :: MISSING      =  1 ! namelist group is missing
  INTEGER, PARAMETER :: LENGTH_ERROR =  2 
  INTEGER, PARAMETER :: READ_ERROR   =  3

  INTEGER, SAVE      :: nnml         = -1 ! default namelist unit

CONTAINS
!==============================================================================
  SUBROUTINE open_nml (file, unit)

    CHARACTER(len=*) ,INTENT(in)           :: file
    INTEGER          ,INTENT(in)           :: unit
    INTEGER :: istat

    nnml = unit
    WRITE (message_text,*) ' This is namelist ' //file
    CALL message('',message_text)
    OPEN (nnml, file=file, iostat=istat, status='old', action='read', &
         delim='apostrophe')

    IF (istat /= 0) THEN
      CALL finish ('open_nml','Could not open '//TRIM(file))
    END IF

  END SUBROUTINE open_nml  
!==============================================================================
  SUBROUTINE close_nml ()

    INTEGER :: istat

    IF (nnml == -1) THEN
      CALL finish ('open_nml','No namelist file opened.')
    END IF

    CLOSE (nnml, IOSTAT=istat)

    IF (istat /= 0) THEN
      CALL finish ('open_nml','Could not close namelist file.')
    END IF

    nnml = -1

  END SUBROUTINE close_nml
!==============================================================================
  SUBROUTINE position_nml (name, unit, REWIND, status)
    
    ! position_nml - position namelist file for reading
    !
    ! Purpose:
    !
    ! To position namelist file at correct place for reading
    ! namelist /name/ (case independent). 
    !

    CHARACTER(len=*), INTENT(in)            :: name   ! namelist group name
    INTEGER,          INTENT(in)  ,OPTIONAL :: unit   ! file unit number
    LOGICAL,          INTENT(in)  ,OPTIONAL :: REWIND ! default: true
    INTEGER,          INTENT(out) ,OPTIONAL :: status ! error return value

    CHARACTER(len=256) :: yline    ! line read
    CHARACTER(len=256) :: test     ! uppercase namelist group name
    INTEGER            :: stat     ! local copy of status variable
    INTEGER            :: ios      ! status variable from read operation
    LOGICAL            :: lrew     ! local copy of rewind flag
    INTEGER            :: iunit    ! local copy of unit number
    INTEGER            :: len_name ! length of requested namelist group name
    CHARACTER          :: ytest    ! character to test for delimiter
    CHARACTER(len=12)  :: code     ! error code printed
    INTEGER            :: ind      ! index from index routine
    INTEGER            :: indc     ! index of comment character (!)

    lrew  = .TRUE.   ; IF (PRESENT(REWIND)) lrew  = REWIND
    iunit =  nnml    ; IF (PRESENT(unit  )) iunit = unit   
    stat  =  MISSING
    code  = 'MISSING'

    len_name = LEN_TRIM (name)

    IF (len_name > LEN(test)) THEN
      stat =  LENGTH_ERROR
      code = 'LENGTH_ERROR'
    ENDIF

    test = tolower (name)

    ! Reposition file at beginning:

    IF (lrew) REWIND (iunit)

    ! Search start of namelist
    
    DO
      IF (stat /= MISSING) EXIT

      yline = ' '
    
      READ (iunit,'(a)',IOSTAT=ios) yline
      IF (ios < 0) THEN
        EXIT  ! MISSING
      ELSE IF (ios > 0) THEN
        stat =  READ_ERROR
        code = 'READ_ERROR'
        EXIT
      END IF

      yline = tolower (yline)
    
      ind = INDEX(yline,'&'//TRIM(test))
    
      IF (ind == 0) CYCLE
      
      indc = INDEX(yline,'!')

      IF (indc > 0 .AND. indc < ind) CYCLE

      ! test for delimiter
    
      ytest = yline(ind+len_name+1:ind+len_name+1)

      IF ( (LGE(ytest,'0') .AND. LLE(ytest,'9')) .OR. &
           (LGE(ytest,'a') .AND. LLE(ytest,'z')) .OR. &
            ytest == '_'                         .OR. &
           (LGE(ytest,'A') .AND. LLE(ytest,'Z'))) THEN
        CYCLE
      ELSE 
        stat = POSITIONED
        BACKSPACE (iunit)
        EXIT
      END IF
    ENDDO
    
    IF (PRESENT (status)) status = stat
      SELECT CASE (stat)
      CASE (POSITIONED)
        RETURN
      CASE (MISSING)
        IF (PRESENT (status)) RETURN
      END SELECT

    CALL finish ('position_nml','namelist /'//TRIM(test)//'/ '//code)

  END SUBROUTINE position_nml
!==============================================================================
END MODULE mo_namelist
