SUBROUTINE maxwind(ulz,vmaxz)

  ! Description:
  !
  ! Computes maximum winds for horizontal diffusion
  ! and diagnostics.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, April 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, December 2000, time control
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
  USE mo_control,       ONLY: nlev
  USE mo_semi_impl,     ONLY: nulev, uvmax, vcheck, vmax 
  USE mo_doctor,        ONLY: nout
  USE mo_decomposition, ONLY: dc => local_decomposition
  USE mo_global_op,     ONLY: maxval_latit
  USE mo_time_control,  ONLY: get_time_step

  IMPLICIT NONE

  !  Array arguments 
  REAL(dp) :: ulz  (           nlev  , dc% nglat)
  REAL(dp) :: vmaxz(           nlev  , dc% nglat)

  ! Local arrays
  INTEGER :: nmaxlev(1)   

  !  Local scalars: 
  INTEGER ::  istep


  !  Executable statements 

!-- 1. Maximum level winds for diffusion

  vmax(:) = maxval_latit (vmaxz(:,:))

  istep = get_time_step()

!-- 2. Maximum wind for diagnostics

  nmaxlev = MAXLOC(vmax)                    ! vertical index of max
  nulev   = nmaxlev(1)                      ! for postatd
  uvmax   = vmax(nulev)                     ! max value (absolute)
!  jlat = ilmax(nulev)                      ! latitude index of max
!  ulat = ASIN(twomu(jlat)*0.5_dp)*180._dp/api    ! for postatd
!  ulm = ulz(nulev,jlat)/sqcst(jlat)        ! for postatd

  ! Check for high windspeeds

  IF (uvmax > vcheck) THEN
    WRITE (nout,'(a,f4.0,a)') ' WARNING! high wind speed: ',uvmax,' m/s'
!    WRITE (nout,*) ' Level: ', nulev, ' Latitude: ', jlat, ' NSTEP= ', istep
    WRITE (nout,*) ' Level: ', nulev, ' NSTEP= ', istep

  END IF

END SUBROUTINE maxwind
