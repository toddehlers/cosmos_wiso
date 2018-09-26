MODULE mo_sub_echam
!-----------------------------------------------------------------------
! This module holds the simple submodule provided by echam. It allows to
! request tracers via the namelist and provides routines for commonly
! used processes. All these routines are called from routines within
! call_submodels.f90
!
! Authors: 
! A.Rhodin MPI/DWD April 2002 original code
!------------------------------------------------------------------------
  !-------------
  ! Modules used
  !-------------
  USE mo_kind,   ONLY: dp
  USE mo_tracer, ONLY: trlist
  IMPLICIT NONE
  !----------------
  ! Public entities
  !----------------
  PUBLIC :: radionucl_sink       ! apply exponential decay
  PRIVATE

CONTAINS
!==============================================================================
  !==========
  ! Processes
  !==========
!------------------------------------------------------------------------------
  SUBROUTINE radionucl_sink (kproma, kbdim, klev, pxtm1, pxtte)
    !------------------------------------------------------------
    ! Description:
    !
    ! Calculates the decrease of tracer concentration for a given 
    ! exponential decay time
    !
    ! Method:
    !
    ! The mass mixing-ratio of tracers is multiplied with
    ! exp(time-step/decay-time).
    ! This routine could also be used for emission or sink
    ! above the surface.
    !
    ! *radionucl_sink* is called from *physc*
    !
    ! Authors:
    !
    ! J. Feichter, MI, August 1991, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    !------------------------------------------------------------
    USE mo_time_control,  ONLY: time_step_len
    !----------
    ! Arguments
    !----------
    INTEGER ,INTENT(in)    :: kproma, kbdim, klev             ! dimensions
    REAL(dp),INTENT(in)    :: pxtm1(kbdim,klev,trlist% ntrac) ! concentr. t-dt
    REAL(dp),INTENT(inout) :: pxtte(kbdim,klev,trlist% ntrac) ! tendency
    !----------------
    ! Local variables
    !----------------
    INTEGER :: jt                 ! tracer index
    REAL(dp):: zxtp1(kproma,klev) ! tracer concentration at t+dt
    REAL(dp):: ztmst              ! time step
    REAL(dp):: zqtmst             ! 1 / time step
    !----------------------
    ! Executable statements
    !----------------------
    ztmst  = time_step_len
    zqtmst = 1.0_dp/ztmst

    DO jt = 1, trlist% ntrac
       IF (trlist% ti(jt)% tdecay /= 0._dp) THEN
          zxtp1(:,:)           = pxtm1(1:kproma,:,jt) + pxtte(1:kproma,:,jt) * ztmst
          zxtp1(:,:)           = EXP (-ztmst/trlist% ti(jt)% tdecay) * zxtp1(:,:)
          pxtte(1:kproma,:,jt) = (zxtp1(:,:)-pxtm1(1:kproma,:,jt))*zqtmst
       END IF
    END DO

  END SUBROUTINE radionucl_sink
!==============================================================================
END MODULE mo_sub_echam
