SUBROUTINE setdyn

  ! Description:
  !
  ! Preset and modify constants in dynamics,initialisation
  ! and general purposes modules.
  !
  ! Method:
  !
  ! *setdyn* is called from *initialize*.
  !
  ! Externals.:
  !
  ! *inicon*        preset constants in *mo_constants*.
  ! *inictl*        preset constants in *mo_control*.
  ! *init_fft992*   preset constants in *mo_fft992*.
  ! *inigau*        preset constants in *mo_gaussgrid*.
  ! *inihyb*        preset constants in *mo_hyb*.
  ! *inhysi*        modify constants in *mo_hyb*.
  ! *helmo*         compute matrix for the *helmoltz equation.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, December 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, November 1998, nudging
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! I. Kirchner, MPI, December 2000, time control
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,             ONLY: dp
  USE mo_exception,        ONLY: message, finish
  USE mo_tmp_buffer,       ONLY: cn
  USE mo_mpi,              ONLY: p_parallel_io, p_parallel, p_io,      &
                                 p_bcast, p_pe
  USE mo_doctor,           ONLY: nout
  USE mo_control,          ONLY: nn, nlon, nlev, nlevp1, nm, nk, nkp1, &
                                 nvclev, vct, lmidatm
  USE mo_truncation,       ONLY: ntrk, ntrm, ntrn, scpar
  USE mo_semi_impl,        ONLY: vcrit, vcheck, hdamp, eps, apr, tr,   &
                                 betadt, betazq
  USE mo_hdiff,            ONLY: nlvstd1, nlvstd2, ldiahdf, enstdif,   &
                                 dampth, damhih, diftcor
  USE mo_hyb,              ONLY: apsurf, inihyb
  USE mo_forecast_switches,ONLY: lumax, lzondia, lsimdt, lsimzq,       &
                                 lvtmpc1, lvtmpc2
  USE mo_constants,        ONLY: a, vtmpc1, vtmpc2
  USE mo_upper_sponge,     ONLY: spdrag, enspodi, nlvspd1, nlvspd2
#ifdef FFT991
  USE mo_fft991,           ONLY: init_fft991
#else
  USE mo_fft992,           ONLY: init_fft992
#endif
  USE mo_gaussgrid,        ONLY: inigau ! module subroutine
  USE mo_time_control,     ONLY: lstart, diagvert, diagdyn,            &
                                 delta_time, p_bcast_event
  USE mo_namelist,         ONLY: position_nml, nnml, POSITIONED

  IMPLICIT NONE

  !  Local scalars: 
  REAL(dp):: zaa, zp0, zt, zt0
  INTEGER :: idt, jk, jl
  INTEGER :: ierr         ! error return value from position_nml

  !  Local arrays: 
  REAL(dp):: zp0a(1), zpf(nlev), zph(nlevp1)

  !  External subroutines 

  EXTERNAL helmo, inhysi, pres, presf, sudif

  INCLUDE 'dynctl.inc'

  !  Executable statements 

!-- 1. Preset constants in modules

!-- 1.2 Preset constants in *mo_forecast_switches*

  lsimdt = .TRUE.
  lsimzq = .TRUE.
  lvtmpc1 = .TRUE.
  lvtmpc2 = .TRUE.
  lumax = .FALSE.
  lzondia = .FALSE.
  ldiahdf = .FALSE.

!-- 1.3 Preset constants for fft

#ifdef FFT991
  CALL init_fft991(nlon)
#else
  CALL init_fft992(nlon)
#endif

!-- 1.4 Preset constants in *mo_gaussgrid*

  CALL inigau

!-- 1.5 Preset constants in *mo_hdiff*

  damhih = 100._dp
  enstdif = 1._dp
  nlvstd1 = 1
  nlvstd2 = 1
  spdrag  = 0.0_dp
  enspodi = 1._dp
  nlvspd1 = 1
  nlvspd2 = 1
  IF (nn==511) THEN
    dampth = 1._dp
  ELSE IF (nn==319) THEN
    dampth = 1._dp
  ELSE IF (nn==255) THEN
    dampth = 1._dp
  ELSE IF (nn==213) THEN
    dampth = 2._dp
  ELSE IF (nn==159) THEN
    dampth = 2._dp
  ELSE IF (nn==106) THEN
    dampth = 3._dp
  ELSE IF (nn==85) THEN
    dampth = 5._dp
  ELSE IF (nn==63) THEN
    dampth = 7._dp
  ELSE IF (nn==42) THEN
    dampth = 9._dp
  ELSE IF (nn==31) THEN
    dampth = 12._dp
    IF(nlev==11) dampth = 15._dp
  ELSE IF (nn==21) THEN
    damhih = 1000._dp
    IF (lmidatm) THEN
      dampth = 192._dp
    ELSE
      dampth = 6._dp
      IF(nlev==11) dampth = 15._dp
    ENDIF
  ELSE
    CALL message('setdyn',' This model resolution is not supported')
    CALL finish('setdyn','Run terminated.')
  END IF

!-- 1.5.1 Preset constants in *mo_diff*

  CALL sudif

!-- 1.6 Preset constants in *mo_truncation*

  CALL scpar (nm, nn, nk)

  IF (nlev/=19) THEN
    CALL message('setdyn','check setting of stratospheric diffusion')
  END IF

  IF (nm==nn .AND. nm==nk) THEN
    IF (nm>=106 .AND. nlev>=19) THEN
      ntrn(1) = nm - 24
      ntrn(2) = nm - 22
      ntrn(3) = nm - 20
      ntrn(4) = nm - 18
      ntrn(5) = nm - 16
      ntrn(6) = nm - 13
      ntrn(7) = nm - 10
      ntrn(8) = nm - 6
      ntrn(9) = nm - 3
      ntrn(10) = nm - 1
!!$    ELSE IF (nm>=106 .AND. nlev>=16) THEN
!!$      ntrn(1) = nm - 22
!!$      ntrn(2) = nm - 18
!!$      ntrn(3) = nm - 14
!!$      ntrn(4) = nm - 10
!!$      ntrn(5) = nm - 6
!!$      ntrn(6) = nm - 3
!!$      ntrn(7) = nm - 1
    ELSE IF (nm==63 .AND. nlev>=19) THEN
      ntrn(1) = 56
      ntrn(2) = 56
      ntrn(3) = 57
      ntrn(4) = 57
      ntrn(5) = 58
      ntrn(6) = 58
      ntrn(7) = 59
      ntrn(8) = 60
      ntrn(9) = 61
      ntrn(10) = 62
!!$    ELSE IF (nm==63 .AND. nlev>=16) THEN
!!$      ntrn(1) = 56
!!$      ntrn(2) = 57
!!$      ntrn(3) = 58
!!$      ntrn(4) = 59
!!$      ntrn(5) = 60
!!$      ntrn(6) = 61
!!$      ntrn(7) = 62
    END IF

  END IF

!-- 1.7 Preset constants in *mo_semi_impl*

  apr = 80000._dp
  tr = 300._dp
  IF (nn==106) THEN
    vcrit = 68._dp
  ELSE
    vcrit = 85._dp
  END IF
  hdamp = 1._dp
  vcheck = 200._dp

!-- 1.8 Preset constants in *mo_hyb*

  CALL inihyb

!-- 1.9 Compute icao-based correction profile for
!       temperature diffusion

  zp0 = 101320._dp
  zp0a(1) = zp0
  zt0 = 288._dp
  zaa = 1._dp/5.256_dp

  CALL pres(zph,1,zp0a,1)
  CALL presf(zpf,1,zph,1)

  DO jk = 1, nlev
     zt = zt0*((zpf(jk)/zp0)**zaa)
     IF (zt>216.5_dp) THEN
        diftcor(jk) = 0.5_dp*(vct(nvclev+jk)+vct(nvclev+jk+1))*zaa*zt*zp0/zpf(jk)
     ELSE
        diftcor(jk) = 0._dp
     END IF
  END DO

  !-- 2. Read namelists

  IF (p_parallel_io) THEN
     CALL position_nml ('DYNCTL', status=ierr)
     SELECT CASE (ierr)
     CASE (POSITIONED)
       READ (nnml, dynctl)
     END SELECT
  ENDIF
  IF (p_parallel) THEN
     CALL p_bcast_event (diagdyn,  p_io)
     CALL p_bcast_event (diagvert, p_io)
     CALL p_bcast (ntrn, p_io)
     CALL p_bcast (nlvstd1, p_io)
     CALL p_bcast (nlvstd2, p_io)
     CALL p_bcast (lumax, p_io)
     CALL p_bcast (lzondia, p_io)
     CALL p_bcast (ldiahdf, p_io)
     CALL p_bcast (vcrit, p_io)
     CALL p_bcast (hdamp, p_io)
     CALL p_bcast (enstdif, p_io)
     CALL p_bcast (apsurf, p_io)
     CALL p_bcast (vcheck, p_io)
     CALL p_bcast (eps, p_io)
     CALL p_bcast (dampth, p_io)
     CALL p_bcast (damhih, p_io)
     CALL p_bcast (spdrag, p_io)
     CALL p_bcast (enspodi, p_io)
     CALL p_bcast (nlvspd1, p_io)
     CALL p_bcast (nlvspd2, p_io)
  ENDIF

  !-- 3. Modify and derive some constants

  !-- 3.1 Modify constants in *mo_constants*

  IF ( .NOT. lvtmpc1) vtmpc1 = 0._dp
  IF ( .NOT. lvtmpc2) vtmpc2 = 0._dp

!-- 3.2 Modify constants in *mo_truncation*

  IF (nm==nn .AND. nm==nk) THEN

    DO jl = 1, nlev
      ntrm(jl) = ntrn(jl)
      ntrk(jl) = ntrn(jl)
    END DO

  ELSE

    DO jl = 1, nlev
      ntrn(jl) = nn
    END DO

  END IF

!-- 3.4 Modify constants in *mo_semi_impl*

  IF (lmidatm) THEN
    ! vcrit is defined as the courant velocity: a/(nn*delta_time)
    ! multiplied by the truncation nn. it is used in hdiff.
    vcrit = a/delta_time
  ELSE
    ! The following statement is based on the empirical
    ! constatation that a 20 minutes time step was always safe
    ! for t63 if the wind does not exceed 85(vcrit)m/s.
    vcrit = vcrit*1200._dp*63._dp/delta_time
  ENDIF

  IF (p_pe == p_io) THEN
     WRITE (nout, '(/,a,e12.3,/,a,e12.3,/)') &
          ' Damping factor for strong stratospheric damping      (hdamp) = ', &
          hdamp, &
          ' Critical velocity for horizontal diffusion enhancing (vcrit) = ', &
          vcrit
  END IF


  IF (lsimdt) THEN
    betadt = 0.75_dp
  ELSE
    betadt = 0._dp
  END IF

  IF (lsimzq) THEN
    betazq = 1._dp
  ELSE
    betazq = 0._dp
  END IF

!-- 3.5 Modify constants in *mo_hyb*

  CALL inhysi

!-- 3.6 Compute matrix for the *helmoltz equation

  IF (.NOT. ALLOCATED(cn)) ALLOCATE (cn(nlev,nlev,nkp1))
  idt = 2
  IF (lstart) idt = 1
  CALL helmo(idt)

!-- 4. Write namelist values

  IF (.NOT. p_parallel) THEN
    WRITE (nout,dynctl)
  END IF

  RETURN
END SUBROUTINE setdyn
