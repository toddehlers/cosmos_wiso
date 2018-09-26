MODULE mo_diag_radiation

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_diag_radiation
  PUBLIC :: print_diag_radiation

  PUBLIC :: diag_rad1, diag_rad2

  REAL(dp), ALLOCATABLE :: diag_rad1(:,:)
  REAL(dp)              :: diag_rad2(7)

CONTAINS

  SUBROUTINE init_diag_radiation

    USE mo_control, ONLY: nlevp1

    LOGICAL, SAVE :: not_used = .TRUE.
    
    IF (not_used) THEN

      ALLOCATE (diag_rad1(nlevp1,4))

      not_used = .FALSE.

    END IF

    diag_rad1(1:nlevp1,1:4) = 0.0_dp
    diag_rad2(1:7)          = 0.0_dp

  END SUBROUTINE init_diag_radiation

  SUBROUTINE print_diag_radiation

    ! Description:
    !
    ! Complets and prints diagnostics for radiation.
    !
    ! Method:
    !
    ! print_diag_radiation is called from scan1
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, June 1983, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! I. Kirchner, MPI, November 2000, date/time control
    ! L. Kornblueh, MPI, October 2003, packed in module and parallelized
    ! 

    USE mo_control, ONLY: nlevp1, nlev
    USE mo_doctor,  ONLY: nout
    USE mo_mpi,     ONLY: p_field_sum, p_parallel_io

    !  Local scalars: 
    REAL(dp) :: zalbpr, zdift1, zdift2, zdift3, zdift4, zdift5
    INTEGER :: jk

    !  Local arrays: 
    REAL(dp) :: zdiat(nlevp1), zdiag1(nlevp1,4), zdiag2(7)

    !  Executable statements 

    zdiag1(:,1) = p_field_sum(diag_rad1(:,1))
    zdiag1(:,2) = p_field_sum(diag_rad1(:,2))
    zdiag1(:,3) = p_field_sum(diag_rad1(:,3))
    zdiag1(:,4) = p_field_sum(diag_rad1(:,4))

    zdiag2(:)   = p_field_sum(diag_rad2(:))

    IF (p_parallel_io) THEN

      !-- 1. Complete diagnostics

      zdiat(1:nlevp1) = zdiag1(1:nlevp1,3) + zdiag1(1:nlevp1,4)
      zdift1 = zdiag2(1) - zdiag2(2)
      zdift2 = zdiag2(4) - zdiag2(5)
      zdift3 = zdiag2(6) - zdiag2(7)
      zdift4 = zdift1 - zdiag2(3)
      zdift5 = zdift2 - zdift3
      IF (zdiag2(1) > 0.01_dp) THEN
        zalbpr = zdiag2(2)/zdiag2(1)*100.0_dp
      ELSE
        zalbpr = 0.0_dp
      END IF
      WRITE (nout, '(/,a,/,a,/)') &
           ' Global Radiation horizontally averaged:', &
           '   (Fluxes - T=Top, B=Bottom, D=Down, U=Up, S=Solar, L=Thermal, N=Net)'
      WRITE (nout, '(a)') ' Temperature  [K]'
      WRITE (nout, '(10f6.1)') zdiag1(1:nlev,1)
      WRITE (nout, '(a)') ' Cloud Cover     '
      WRITE (nout, '(10f6.1)') zdiag1(1:nlev,2)
      WRITE (nout, '(a)') ' Shortwave heating rates [K/day]  '
      WRITE (nout, '(10f6.1)') zdiag1(1:nlev,3)
      WRITE (nout, '(a)') ' Longwave heating rates [K/day]  '
      WRITE (nout, '(10f6.1)') zdiag1(1:nlev,4)
      WRITE (nout, '(a)') ' Net heating rates [K/day]  ' 
      WRITE (nout, '(10f6.1)') zdiat(1:nlev)
      
      WRITE (nout, '(a)') ' Fluxes:'
      WRITE (nout, '(3(a,f6.1))') &
           '   TSD = ', zdiag2(1), ' TSU = ', zdiag2(2), ' TLU = ', zdiag2(3)
      WRITE (nout, '(4(a,f6.1))') &
           '   BSD = ', zdiag2(4),  ' BSU = ', zdiag2(5), &
           ' BLU = ', zdiag2(6), ' BLD = ', zdiag2(7)
      WRITE (nout, '(5(a,f6.1))') '   TS  = ', zdift1, &
           ' BS  = ',zdift2 ,' BL  = ', zdift3, ' TND = ', zdift4,&
           ' BND = ', zdift5
      WRITE (nout, '(a,f5.1,a)') '   Planetary Albedo = ', zalbpr, ' %'
      WRITE (nout, '(2/)')

    END IF

  END SUBROUTINE print_diag_radiation

END MODULE mo_diag_radiation
