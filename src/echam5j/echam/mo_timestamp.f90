MODULE mo_timestamp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_timestamp

CONTAINS

  SUBROUTINE write_timestamp(application)
    
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: application

    CHARACTER(len= 8) :: date
    CHARACTER(len=10) :: time
    CHARACTER(len=10) :: iso_date
    CHARACTER(len=12) :: new_time

    CALL date_and_time(date, time)

    ! reformat date to ISO-8601

    WRITE(iso_date,'(a4,a1,a2,a1,a2)') &
         date(1:4), '-', date(5:6), '-', date(7:8)

    ! reformat time string

    WRITE (new_time,'(a2,a1,a2,a1,a2,a1,a3)') &
         time(1:2), ':', time(3:4), ':', time(5:6), '.', time(8:10)

    ! write out ...

    IF (PRESENT(application)) THEN
      WRITE (0,*) TRIM(application), ': ', iso_date, ' ',new_time
    ELSE
      WRITE (0,*) iso_date, ' ',new_time
    ENDIF

  END SUBROUTINE write_timestamp

END MODULE mo_timestamp
