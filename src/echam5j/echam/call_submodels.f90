!
! This file holds the calls to the entries of the submodels 
! attached to ECHAM.
!
! Change this file to attach submodels to ECHAM.
!
!==============================================================================
  SUBROUTINE call_init_submodels
  !
  ! This routine is called from 'initialize'
  !
  ! new modules must be introduced from here (call new_submodel)
  !
    USE mo_jsbach_interface, ONLY: jsbach_init
    USE mo_co2,              ONLY: init_submodel_co2

    CALL jsbach_init
    CALL init_submodel_co2

  END SUBROUTINE call_init_submodels 
!==============================================================================
  SUBROUTINE call_request_tracer
  !
  ! This routine is called from 'initrac', module 'mo_tracer'
  ! 
  ! tracer request routines are called from here
  !
  ! Template for a tracer request routine:
  !
  !    SUBROUTINE request_tracer
  !      USE mo_tracer, ONLY: new_tracer
  !      IMPLICIT NONE
  !      CALL new_tracer('tracername1' ,'modulename', ...)
  !      CALL new_tracer('tracername2' ,'modulename', ...)
  !    END SUBROUTINE request_tracer
  !
    USE mo_sub_nml, ONLY: request_tracer_nml, set_tracer_nml
    USE mo_co2,     ONLY: init_co2
    IMPLICIT NONE
    !
    ! routine to request tracers via namelist group /NEW_TRACER/
    !
    !CALL request_tracer_nml
    !CALL set_tracer_nml

    CALL init_co2

  END SUBROUTINE call_request_tracer
!==============================================================================
  SUBROUTINE call_init_submodel_memory
  !
  ! This routine is called from subroutine 'init_memory', module
  ! 'mo_memory_streams'. Routines are called from here to allocate memory
  ! and to define output streams.
  !
  !  USE mo_jsbach_io, ONLY: jsbach_init_memory
    USE mo_accuflx_mem,  ONLY: construct_stream_accuflx
    USE mo_radiation,    ONLY: ldblrad

    IMPLICIT NONE

    IF (ldblrad) THEN
       CALL construct_stream_accuflx
    ENDIF

  !  CALL jsbach_init_memory

  END SUBROUTINE call_init_submodel_memory
!------------------------------------------------------------------------------
  SUBROUTINE call_free_submodel_memory
  !
  ! Routines are called from here to deallocate memory
  !
  !  USE mo_jsbach_io, ONLY: jsbach_free_memory

    IMPLICIT NONE

  ! CALL jsbach_free_memory

  END SUBROUTINE call_free_submodel_memory
!==============================================================================
  SUBROUTINE call_chem1 (kproma, kbdim, klev, klevp1, ktrac, paphp1, papp1, &
                         ptte, ptm1, pxtm1, pxtte)

    !
    ! This subroutine calls the different parts of the chemistry
    ! module.
    !
    ! It is called by xtdriver (module mo_tracer) which is called by 'physc'
    !
    USE mo_kind,      ONLY : dp
    USE mo_sub_echam, ONLY : radionucl_sink
    IMPLICIT NONE

    ! Scalar arguments
    INTEGER :: kproma, kbdim, klev, klevp1, ktrac

    ! Array arguments
    REAL(dp):: paphp1(kbdim,klevp1), papp1(kbdim,klev), ptte(kbdim,klev),  &
               ptm1(kbdim,klev)
    REAL(dp):: pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac)

    !
    ! commonly used processes calculated by ECHAM:
    !
    ! EXAMPLE:
!!$    CALL radionucl_sink (kproma, kbdim, klev, pxtm1, pxtte)

  END SUBROUTINE call_chem1
!------------------------------------------------------------------------------
  SUBROUTINE call_chem2 (kproma, kbdim, klev, klevp1, ktrac, paphp1, papp1, &
                         ptte, ptm1, pxtm1, pxtte)
    !
    ! similar to call_chem1, but called from other place within physc
    USE mo_kind,      ONLY : dp
    USE mo_sub_echam, ONLY : radionucl_sink
    IMPLICIT NONE

    ! Scalar arguments
    INTEGER :: kproma, kbdim, klev, klevp1, ktrac

    ! Array arguments
    REAL(dp):: paphp1(kbdim,klevp1), papp1(kbdim,klev), ptte(kbdim,klev),  &
               ptm1(kbdim,klev)
    REAL(dp):: pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac)

  END SUBROUTINE call_chem2
!==============================================================================
  SUBROUTINE call_diagn
  !
  ! This routine is called by xtdiagn (mo_tracer) called at the end of physc
  !  physc. Diagnostics at the end of a time step may be called from here
  !
    IMPLICIT NONE
 END SUBROUTINE call_diagn
!==============================================================================
  SUBROUTINE call_init_tracers
    USE mo_co2, ONLY: reference_co2
  !  
  ! This routine is called from xtini to allow initialization of tracer
  ! variables besides initialization with constant values or from the
  ! rerun file.
  !
  ! Base the decision whether to initialize on the following conditions:
  !
  ! trlist% ti(jt)% init   == 0       (not initialized so far)
  ! trlist% ti(jt)% ninit: >= INITIAL (initialisation requested)
  !
  ! Set trlist%ti(jt) to a value =/0 afterwards.
  !

    CALL reference_co2

  END SUBROUTINE call_init_tracers
!==============================================================================
  SUBROUTINE call_read_bcond
  !
  ! This routine is called from 'stepon'. Routines which read forcing
  ! data may be called from here.
  !
    USE mo_time_control,    ONLY: current_date, get_date_components
    USE mo_co2,            ONLY: read_co2_emission

    IMPLICIT NONE

    INTEGER :: imonth

    CALL get_date_components (current_date, month=imonth)

    CALL read_co2_emission

  END SUBROUTINE call_read_bcond
!==============================================================================
  SUBROUTINE call_chem_bcond (kproma, kbdim, klev, pxtems, pxtte, jrow)

    !
    ! This is the interface routine between ECHAM's standard tracer treatment
    ! and the chemistry routines.
    !
    ! This routine is called from vdiff
    !
    ! M. Schultz and Hans-Stefan Bauer, MPI Hamburg, November 2000
    !

    USE mo_kind,            ONLY: dp
    USE mo_tracer,         ONLY: ntrac
    USE mo_co2,           ONLY: co2_emissions

    IMPLICIT NONE
    
    ! Arguments

    INTEGER, INTENT(in)    :: kproma, kbdim, klev, jrow
    REAL(dp),INTENT(inout) :: pxtte(kbdim,klev,ntrac), &
                              pxtems(kbdim,ntrac)

    CALL co2_emissions(kproma, kbdim, ntrac, pxtems, jrow)

  END SUBROUTINE call_chem_bcond
!==============================================================================
