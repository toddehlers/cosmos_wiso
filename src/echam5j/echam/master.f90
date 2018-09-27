PROGRAM master

  !----------------------------------------------------------------------------
  !
  ! Copyright 2000-2005 by Max Planck Institute for Meteorology
  !
  ! This software can be used only under the conditions of the 
  ! "MPI-M Software Licence Agreement", which must be signed 
  ! by each user.
  !
  ! The "MPI-M Software Licence Agreement" and the information on
  ! the distribution of MPI-M models are given on:
  !
  ! http://www.mpimet.mpg.de/en/extra/models/distribution/index.php
  !
  ! The ECHAM5-webpage can be found at:
  !
  ! http://www.mpimet.mpg.de/en/extra/models/echam/echam5.php
  !
  !----------------------------------------------------------------------------
  !
  ! Call the control subroutine (*control*).
  !
  ! Externals:
  !
  ! *control*   called to control the run.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 2000, date/time control
  ! S. Legutke, MPI M&D, Juli 2001, redirect stdout for coupling
  ! 
  ! for more details see file AUTHORS

  USE mo_kind,         ONLY: dp
  USE mo_doctor,       ONLY: nout, nin
  USE mo_mpi,          ONLY: p_start, p_stop, p_pe, p_io
  USE mo_time_control, ONLY: lbreak, lstop, labort
  USE mo_exception,    ONLY: finish, message
  USE mo_util_string,  ONLY: separator
#if defined (__oasis)
  USE mo_couple,       ONLY: couple_quit
#endif
#ifdef _OPENMP
  USE omp_lib,         ONLY: omp_set_dynamic, omp_set_num_threads
#endif

  IMPLICIT NONE

  !  External functions 
  REAL(dp), EXTERNAL :: util_walltime

  !  External subroutines 
  EXTERNAL control
#ifdef __XT3__
  EXTERNAL :: util_base_iobuf
#endif

  REAL(dp) :: zwtime

#ifdef __XT3__
  ! set buffer for stderr and stdout

  CALL util_base_iobuf
#endif

#ifdef _OPENMP
  CALL omp_set_dynamic(.FALSE.)
  CALL omp_set_num_threads(1)
#endif

  ! Initialize wallclock timer

!$OMP PARALLEL
!$OMP MASTER
  zwtime = util_walltime()
!$OMP END MASTER
!$OMP END PARALLEL

  ! Start MPI

  CALL p_start


  IF (p_pe == p_io) THEN
#if defined (__prism) 
!--  Redirect standard output file to atmout if coupled run.
 
        OPEN (UNIT=nout,FILE='atmout',STATUS='UNKNOWN',FORM ='FORMATTED')
        WRITE(nout,*)' '
        WRITE(nout,*)'Atmosphere standard output is assigned to file atmout.'
        WRITE(nout,*)' '
#endif

  ! Print version

        
     WRITE (nout,separator)

     WRITE (nout,'(/,a,/,a,/)')                                        &
          '  ECHAM         - Release  5.4.01 ',                        &
          '  Copyright by Max-Planck-Institute for Meteorology, 2007', &
          '  Read master.f90 and licence.pdf before using ECHAM5'

     WRITE (nout,*)
     
!---wiso-code

     WRITE (nout,'(a)') '  ++ Added water isotope code, Version 2.3'
     WRITE (nout,'(a)') '  ++ M. Werner, P. Langebroek, B. Haese'
     WRITE (nout,'(a)') '  ++ AWI Bremerhaven'
     
!---wiso-code-end

     WRITE (nout,separator)
     WRITE (nout,*)
  END IF

  DO                               ! Loop over rerun cycles
     CALL control
     IF (lbreak .OR. lstop) EXIT
     CALL message('','Start next rerun cycle.')
  END DO

 CALL p_stop                      ! Stop MPI

  IF (lstop) THEN
    CALL message('','Experiment finished.')
    IF (labort) CALL finish('master','Run terminated.',1)
  ELSE
    CALL message('','Experiment stopped.')
  END IF

END PROGRAM master
