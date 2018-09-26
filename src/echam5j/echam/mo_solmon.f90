MODULE mo_solmon

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_solmon
  PUBLIC :: interpolate_solmon

  PUBLIC :: reff_stratm
  PUBLIC :: solcm
  PUBLIC :: solcm0, reff_stratm0

  REAL(dp) :: reff_stratm(12)
  REAL(dp) :: solcm(12)
  REAL(dp) :: solcm0, reff_stratm0

CONTAINS

  SUBROUTINE init_solmon

    ! Description:
    !
    ! Read 12 month of solc and reff
    !
    ! Method:
    !
    ! Time series are read from file solmon
    ! 
    ! Authors:
    !
    ! M. Esch, MPI, Jan 2005, original source
    !
    ! for more details see file AUTHORS
    !

    USE mo_control,      ONLY: lcouple, lipcc, lsolc, lreff
    USE mo_doctor,       ONLY: nout
    USE mo_mpi,          ONLY: p_parallel_io, p_io, p_bcast, p_parallel
    USE mo_time_control, ONLY: get_date_components, next_date, lstart
    USE mo_filename,     ONLY: find_next_free_unit

    INTEGER :: pyr, pmo                  ! p* present date/time
    INTEGER :: mm, nn, i, nsolunit

    IF (p_parallel_io) THEN 
      nsolunit = find_next_free_unit (51,100)

      OPEN(nsolunit,file='solmon',status='old',form='formatted',action='read')
      CALL get_date_components(next_date,pyr,pmo)
      ! WRITE(nout,*) 'next_date,pyr,pmo ', next_date,pyr,pmo
      ! set to default values for NAG ompiler
      solcm(12) = 1365._dp
      IF(lcouple .OR. lipcc) solcm(12) = 1367._dp
      reff_stratm(12) = 0.2_dp

      DO
        solcm0 = solcm(12)
        reff_stratm0 = reff_stratm(12)
        DO i = 1,12
          READ(nsolunit,*) mm,nn,solcm(i),reff_stratm(i)
        END DO
        IF (nn == pyr) EXIT
      END DO

      DO i=1,12
        WRITE(nout,*) 'SOLC and REFF_STRAT set at time to ', &
             i,nn,solcm(i),reff_stratm(i)
      END DO
      WRITE(nout,*) 'READ from unit', nsolunit
    ENDIF

    IF (p_parallel) THEN
      CALL p_bcast (solcm0, p_io)
      CALL p_bcast (solcm, p_io)
      CALL p_bcast (reff_stratm0, p_io)
      CALL p_bcast (reff_stratm, p_io)
    ENDIF

  END SUBROUTINE init_solmon

  SUBROUTINE interpolate_solmon

    USE mo_exception,       ONLY: message, finish, message_text
    USE mo_time_conversion, ONLY: time_native, TC_get, TC_convert
    USE mo_time_control,    ONLY: current_date, get_time_step, lresume
    USE mo_mpi,             ONLY: p_parallel_io
    USE mo_radiation,       ONLY: solc, reff_strat
    USE mo_control,         ONLY: lsolc, lreff

    TYPE(time_native) :: solmon_date 

    INTEGER :: yr, mo, dy, hr, mn, se

    ! Description:
    !

    CALL TC_convert(current_date, solmon_date) 
    CALL TC_get (solmon_date, yr, mo, dy, hr, mn, se)

    IF(lsolc) THEN
      IF(lresume .AND. mo==12) THEN
        solc=solcm0
      ELSE
        solc=solcm(mo)
      END IF
    END IF
    !
    reff_strat=0.2_dp
    IF(lreff) THEN
      IF(lresume .AND. mo==12) THEN
        reff_strat=reff_stratm0
      ELSE
        reff_strat=reff_stratm(mo)
      END IF
    END IF

    WRITE (message_text,'(a,i8,a,f9.4,a,f9.4)') &
         'nstep: ', get_time_step(), ' SOLC = ', solc, &
         ' REFF_STRAT = ', reff_strat
    CALL message('', TRIM(message_text))

  END SUBROUTINE interpolate_solmon

END MODULE mo_solmon
