MODULE mo_tmp_buffer

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  ! Matrix for Helmholtz-equation

  REAL(dp), ALLOCATABLE :: cn(:,:,:)   ! setdyn.f  (nlev,nlev,nkp1) -

CONTAINS
  !------------------------------------------------------------------------------
  SUBROUTINE cleanup_tmp_buffer
    !
    ! deallocate module variables
    !
    IF (ALLOCATED(cn)) DEALLOCATE (cn)
  END SUBROUTINE cleanup_tmp_buffer
  !------------------------------------------------------------------------------
END MODULE mo_tmp_buffer
