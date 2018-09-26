MODULE mo_midatm

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

CONTAINS

  SUBROUTINE gwspectrum ( krow, kproma, kbdim, klev  &
       ,paphm1,      papm1                       &
       ,ptm1,        pum1,        pvm1           &
       ,paprflux                                 &
       ,ptte,        pvol,        pvom   )


    !
    ! Description:
    ! 
    !   Hines parameterization from ccc/mam (Hines, 1997a,b):       
    !   physical tendencies of the prognostic variables u,v 
    !   due to vertical transports by a broad band spectrum
    !   of gravity waves. 
    !
    !   Note that diffusion coefficient and  heating rate 
    !   only calculated if iheatcal = 1.
    !
    !        *gwspectrum* is called from *physc*.
    !
    !  Authors:
    !    
    !   c. mclandress   ists   august 1995
    !   n. mcfarlane    cccma  may 1995
    !   m. charron      mpi    2000-2001
    !   e. manzini      mpi    february 2002 (re-write, based on cccgwd)                 
    !   h. schmidt      mpi    march 2003 
    !
    !
    USE mo_constants,         ONLY: cpd, rd, rhoh2o, g
    USE mo_gwspectrum,        ONLY: lextro, lfront, lozpr, iheatcal, &
                                    emiss_lev, rmscon, kstar, naz,  & 
                                    pcrit, pcons 

    ! scalar argument with intent(IN)
    INTEGER ,INTENT(in) ::  krow, kproma, kbdim, klev

    !  Array arguments with intent(IN):
    ! Input 1D
    REAL(dp) ,INTENT(IN) :: paprflux(kbdim)     ! precipitation flux
    ! Input 2D
    REAL(dp) ,INTENT(IN) :: paphm1(kbdim,klev+1) ! half level pressure (t-dt)
    REAL(dp) ,INTENT(IN) :: papm1(kbdim,klev)    ! full level pressure (t-dt)
    REAL(dp) ,INTENT(IN) :: ptm1(kbdim,klev)     ! temperature (t-dt)
    REAL(dp) ,INTENT(IN) :: pum1(kbdim,klev)     ! zonal wind (t-dt)
    REAL(dp) ,INTENT(IN) :: pvm1(kbdim,klev)     ! meridional wind (t-dt)

    !  Array arguments with intent(InOut):
    ! - input/output 2d
    REAL(dp) ,INTENT(INOUT) :: ptte(kbdim,klev)  ! tendency of temperature
    REAL(dp) ,INTENT(INOUT) :: pvol(kbdim,klev)  ! tendency of meridional wind
    REAL(dp) ,INTENT(INOUT) :: pvom(kbdim,klev)  ! tendency of zonal wind

    !  Local arrays for ccc/mam hines gwd scheme:

    ! Important local parameter (passed to all subroutines):
    INTEGER, PARAMETER :: nazmth = 8  ! max azimuth array dimension size 

    REAL(dp) :: pressg(kproma)  ! Surface pressure (pascal)  
    REAL(dp) :: zpr(kproma)     ! precipitation (check dims: echam5 change)

    ! * Vertical positioning arrays and work arrays:                       
    REAL(dp) :: sgj(kproma,klev), shj(kproma,klev), shxkj(kproma,klev)
    REAL(dp) ::  dsgj(kproma,klev), dttdsf(kproma)

    REAL(dp) :: utendgw(kproma,klev) ! zonal tend, gravity wave spectrum (m/s^2)
    REAL(dp) :: vtendgw(kproma,klev) ! merid tend, gravity wave spectrum (m/s^2)
    REAL(dp) :: ttendgw(kproma,klev) ! temperature tend, gravity wave spectrum (K/s)
    REAL(dp) :: diffco(kproma,klev)  ! diffusion coefficient (m^2/s) 

    REAL(dp) ::  flux_u(kproma,klev) ! zonal momentum flux (pascals)
    REAL(dp) ::  flux_v(kproma,klev) ! meridional momentum flux (pascals) 

    REAL(dp) :: uhs(kproma,klev)      ! zonal wind (m/s), input for hines param
    REAL(dp) :: vhs(kproma,klev)      ! merid wind (m/s), input for hines param
    REAL(dp) :: bvfreq(kproma,klev)   ! background brunt vassala frequency (rad/s)
    REAL(dp) :: density(kproma,klev)  ! background density (kg/m^3)
    REAL(dp) :: visc_mol(kproma,klev) ! molecular viscosity (m^2/s) 
    REAL(dp) :: alt(kproma,klev)      ! background altitude (m)

    REAL(dp) :: rmswind(kproma)        ! rms gravity wave  wind, lowest level (m/s) 
    REAL(dp) :: anis(kproma,nazmth)    ! anisotropy factor (sum over azimuths = 1) 
    REAL(dp) :: k_alpha(kproma,nazmth) ! horizontal wavenumber of each azimuth (1/m)
    LOGICAL :: lorms(kproma)       ! .true. for rmswind /=0 at launching level 

    REAL(dp) :: m_alpha(kproma,klev,nazmth) ! cutoff vertical wavenumber (1/m)
    REAL(dp) :: mmin_alpha(kproma,nazmth)   ! minumum value of m_alpha
    REAL(dp) :: sigma_t(kproma,klev)        ! total rms gw wind (m/s)

    ! gw variances from orographic sources (for coupling to a orogwd)
    REAL(dp) :: sigsqmcw(kproma,klev,nazmth), sigmatm(kproma,klev)

    !
    ! Local scalars:
    INTEGER  :: jk, jl
    INTEGER  :: levbot     ! gravity wave spectrum lowest level
    REAL(dp)     :: rgocp, zpcons, hscal, ratio

    IF (lextro) THEN

       !
       !--  Initialize the ccc/mam hines gwd scheme
       !

       utendgw(:,:) = 0.0_dp
       vtendgw(:,:) = 0.0_dp
       ttendgw(:,:) = 0.0_dp

       diffco(:,:) = 0.0_dp

       flux_u(:,:) = 0.0_dp
       flux_v(:,:) = 0.0_dp 

       uhs(:,:) = 0.0_dp
       vhs(:,:) = 0.0_dp

       ! Wind variances form orographic gravity waves
       ! Note: the code is NOT fully implemeted for this case!

       sigsqmcw(:,:,:) = 0.0_dp
       sigmatm(:,:)    = 0.0_dp

       ! precipitation (check the units!):
       zpcons = (1000.0_dp*86400.0_dp)/rhoh2o
       zpr(1:kproma)=zpcons*paprflux(1:kproma)

       rgocp=rd/cpd

       ! Vertical positioning arrays: 
       DO jk=1,klev
          DO jl=1,kproma
             shj(jl,jk)=papm1(jl,jk)/paphm1(jl,klev+1)
             sgj(jl,jk)=papm1(jl,jk)/paphm1(jl,klev+1)
             dsgj(jl,jk)=(paphm1(jl,jk+1)-paphm1(jl,jk))/paphm1(jl,klev+1)
             shxkj(jl,jk)=(papm1(jl,jk)/paphm1(jl,klev+1))**rgocp 
          END DO
       END DO

       ! Surface pressure: 
       DO jl=1,kproma
          pressg(jl)=paphm1(jl,klev+1)
       END DO

       !
       !     * calculate b v frequency at all points
       !     * and smooth bvfreq.

       DO jk=2,klev
         DO jl=1,kproma
           dttdsf(jl)=(ptm1(jl,jk)/shxkj(jl,jk)-ptm1(jl,jk-1)/shxkj(jl,jk-1)) &
                     /(shj(jl,jk)-shj(jl,jk-1))
           dttdsf(jl)=MIN(dttdsf(jl), -5.0_dp/sgj(jl,jk))
           dttdsf(jl)=dttdsf(jl)*(sgj(jl,jk)**rgocp)

           bvfreq(jl,jk)=SQRT(-dttdsf(jl)*sgj(jl,jk)/rd)*g/ptm1(jl,jk)
         END DO
       END DO

       bvfreq(:,1) = bvfreq(:,2)

       DO jk=2,klev
          DO jl=1,kproma
             ratio=5.0_dp*LOG(sgj(jl,jk)/sgj(jl,jk-1))
             bvfreq(jl,jk) = (bvfreq(jl,jk-1) + ratio*bvfreq(jl,jk))/(1.0_dp+ratio)
          END DO
       END DO

       !     * altitude and density at bottom.

       alt(:,klev) = 0.0_dp

       DO jl=1,kproma
          hscal = rd * ptm1(jl,klev) / g
          density(jl,klev) = sgj(jl,klev) * pressg(jl) / (g*hscal)
       END DO

       !     * altitude and density at remaining levels.

       DO jk=klev-1,1,-1
          DO jl=1,kproma
             hscal = rd * ptm1(jl,jk) / g
             alt(jl,jk) = alt(jl,jk+1) + hscal * dsgj(jl,jk) / sgj(jl,jk)
             density(jl,jk) = sgj(jl,jk) * pressg(jl) / (g*hscal)
          END DO
       END DO

       !
       !     * set molecular viscosity to a very small value.
       !     * if the model top is greater than 100 km then the actual
       !     * viscosity coefficient could be specified here.


       visc_mol(:,:) = 1.5e-5_dp


       ! use single value for azimuthal-dependent horizontal wavenumber:
       ! kstar = (old latitudinal dependence, introduce here if necessary)

       k_alpha(:,:) = kstar

       !     * defile bottom launch level (emission level of gws)

       levbot = klev-emiss_lev   

       !     * initialize switch for column calculation

       lorms(:) = .FALSE.

       !     * background wind minus value at bottom launch level.

       DO jk=1,levbot
         DO jl=1,kproma 
           uhs(jl,jk) = pum1(jl,jk) - pum1(jl,levbot)
           vhs(jl,jk) = pvm1(jl,jk) - pvm1(jl,levbot)
         END DO
       END DO

       !     * specify root mean square wind at bottom launch level.

       DO jl=1,kproma 
          rmswind(jl) = rmscon
          anis(jl,:)   = 1.0_dp/REAL(naz,dp)
       END DO

       !     * gravity waves from fronts:
       IF (lfront) THEN
          CALL gw_fronts(krow, kproma, nazmth, rmswind, anis)
       ENDIF

       !     * modulation by precipitation:
       IF (lozpr) THEN
          DO jl=1,kproma 
             IF (zpr(jl) .GT. pcrit) THEN
                rmswind(jl) = rmscon + ( (zpr(jl)-pcrit)/zpr(jl) )*pcons
             ENDIF
          END DO
       ENDIF

       DO jl=1,kproma 
          IF (rmswind(jl) .GT. 0.0_dp) THEN
             lorms(jl) = .TRUE.
          ENDIF
       END DO

       !
       !     * calculate gw tendencies (note that diffusion coefficient and
       !     * heating rate only calculated if iheatcal = 1).
       !
       CALL hines_extro ( kproma, klev, nazmth,                          & 
                          utendgw, vtendgw, ttendgw, diffco,             &
                          flux_u, flux_v,                                & 
                          uhs, vhs, bvfreq, density, visc_mol, alt,      & 
                          rmswind, anis, k_alpha, sigsqmcw,              &
                          m_alpha,  mmin_alpha ,sigma_t, sigmatm,        & 
                          levbot, lorms)


       !   update tendencies: 
       !
       DO jk=1, klev
          DO jl=1,kproma
             pvom(jl,jk) = pvom(jl,jk) + utendgw(jl,jk)
             pvol(jl,jk) = pvol(jl,jk) + vtendgw(jl,jk)
          END DO
       END DO
       !
       IF ( iheatcal == 1 ) THEN  
          DO jk=1, klev
             DO jl=1,kproma   
                ptte(jl,jk) = ptte(jl,jk) + ttendgw(jl,jk)
             END DO
          END DO
       ENDIF
       !

    ENDIF
    !     * end of hines calculations.

    !-----------------------------------------------------------------------
  END SUBROUTINE gwspectrum

  SUBROUTINE hines_extro ( nlons, nlevs, nazmth,                          &
                           drag_u, drag_v, heat, diffco, flux_u, flux_v,  &
                           vel_u, vel_v, bvfreq, density, visc_mol, alt,  &
                           rmswind, anis, k_alpha, sigsqmcw,               &
                           m_alpha,  mmin_alpha, sigma_t, sigmatm,        &
                           lev2,  lorms)
    !
    !  main routine for hines' "extrowave" gravity wave parameterization based
    !  on hines' doppler spread theory. this routine calculates zonal
    !  and meridional components of gravity wave drag, heating rates
    !  and diffusion coefficient on a longitude by altitude grid.
    !  no "mythical" lower boundary region calculation is made. 
    !
    !  aug. 13/95 - c. mclandress
    !  sept. /95  - n. mcfarlane
    !  1995- 2002 - e. manzini
    !
    !  modifications:
    !
    !  output arguements:
    !
    !     * drag_u = zonal component of gravity wave drag (m/s^2).
    !     * drag_v = meridional component of gravity wave drag (m/s^2).
    !     * heat   = gravity wave heating (k/sec).
    !     * diffco = diffusion coefficient (m^2/sec)
    !     * flux_u = zonal component of vertical momentum flux (pascals)
    !     * flux_v = meridional component of vertical momentum flux (pascals)
    !
    !  input arguements:
    !
    !     * vel_u      = background zonal wind component (m/s).
    !     * vel_v      = background meridional wind component (m/s).
    !     * bvfreq     = background brunt vassala frequency (radians/sec).
    !     * density    = background density (kg/m^3) 
    !     * visc_mol   = molecular viscosity (m^2/s)
    !     * alt        = altitude of momentum, density, buoyancy levels (m)
    !     *              (note: levels ordered so that alt(i,1) > alt(i,2), etc.)
    !     * rmswind   = root mean square gravity wave wind at lowest level (m/s).
    !     * anis      = anisotropy factor (sum over azimuths is one)
    !     * lorms     = .true. for drag computation (column selector)
    !     * k_alpha    = horizontal wavenumber of each azimuth (1/m).
    !     * lev2       = index of last level (eg bottom) for drag calculation 
    !     *              (i.e., lev1 < lev2 <= nlevs).
    !     * nlons      = number of longitudes.
    !     * nlevs      = number of vertical levels.
    !     * nazmth     = azimuthal array dimension (nazmth >= naz).
    ! 
    !  ouput diagnostics:
    !     * m_alpha      = cutoff vertical wavenumber (1/m).
    !     * mmin_alpha   = minimum value of cutoff wavenumber.
    !     * sigma_t      = total rms horizontal wind (m/s).

    !  work arrays.
    !

    !     * v_alpha      = wind component at each azimuth (m/s) and if iheatcal=1
    !     *                holds vertical derivative of cutoff wavenumber.
    !     * sigma_alpha  = total rms wind in each azimuth (m/s).
    !     * ak_alpha     = spectral amplitude factor at each azimuth 
    !     *                (i.e.,{ajkj}) in m^4/s^2.
    !     * densb        = background density at bottom level.
    !     * bvfb         = buoyancy frequency at bottom level and
    !     *                work array for icutoff = 1.
    !
    !     * losigma_t    =  .true. for total sigma not zero

    USE mo_gwspectrum,     ONLY: iheatcal,                        &
                                 kstar, m_min,                    &
                                 naz, slope, f1, f2, f3, f5, f6,  &
                                 icutoff, alt_cutoff, smco, nsmax

    INTEGER :: nlons, nlevs, nazmth, lev2

    REAL(dp):: drag_u(nlons,nlevs),   drag_v(nlons,nlevs) 
    REAL(dp):: heat(nlons,nlevs),     diffco(nlons,nlevs)
    REAL(dp):: flux_u(nlons,nlevs),   flux_v(nlons,nlevs)
    REAL(dp):: flux(nlons,nlevs,nazmth)
    REAL(dp):: vel_u(nlons,nlevs),    vel_v(nlons,nlevs)
    REAL(dp):: bvfreq(nlons,nlevs),   density(nlons,nlevs)
    REAL(dp):: visc_mol(nlons,nlevs), alt(nlons,nlevs)
    REAL(dp):: rmswind(nlons),      bvfb(nlons),   densb(nlons)
    REAL(dp):: anis(nlons,nazmth)
    REAL(dp):: sigma_t(nlons,nlevs), sigsqmcw(nlons,nlevs,nazmth)
    REAL(dp):: sigma_alpha(nlons,nlevs,nazmth), sigmatm(nlons,nlevs)

    REAL(dp):: m_alpha(nlons,nlevs,nazmth), v_alpha(nlons,nlevs,nazmth)
    REAL(dp):: ak_alpha(nlons,nazmth),      k_alpha(nlons,nazmth)
    REAL(dp):: mmin_alpha(nlons,nazmth)    
    REAL(dp):: smoothr1(nlons,nlevs), smoothr2(nlons,nlevs)

    LOGICAL :: lorms(nlons), losigma_t(nlons,nlevs)
    !
    !  internal variables.
    !
    INTEGER :: i, n, l, lev1, il1, il2, iprint

    !----------------------------------------------------------------------- 
    !

    ! range of longitude index:
    il1 = 1      
    il2 = nlons     

    lev1=1              ! top level index

    iprint = 0       !     * iprint     = 1 to print out various arrays.

    !
    !  buoyancy and density at bottom level.
    !
    DO i = il1,il2
       bvfb(i)  = bvfreq(i,lev2)
       densb(i) = density(i,lev2)
    END DO
    !
    !  initialize some variables
    !
    DO n = 1,naz
       DO l=lev1,lev2
          DO i=il1,il2
             m_alpha(i,l,n) =  m_min
          END DO
       END DO
    END DO
    !
    !  compute azimuthal wind components from zonal and meridional winds.
    !
    CALL hines_wind ( v_alpha,   & 
                      vel_u, vel_v, naz,   &
                      il1, il2, lev1, lev2, nlons, nlevs, nazmth )
    !
    !  calculate cutoff vertical wavenumber and velocity variances.
    !
    CALL hines_wavnum ( m_alpha, sigma_t, sigma_alpha, ak_alpha,   &
                        mmin_alpha, losigma_t,                     &
                        v_alpha, visc_mol, density, densb,         &
                        bvfreq, bvfb, rmswind, anis, lorms,        &
                        sigsqmcw, sigmatm,                         &
                        il1, il2, lev1, lev2, nlons, nlevs, nazmth)
    !
    !  smooth cutoff wavenumbers and total rms velocity in the vertical 
    !  direction nsmax times, using flux_u as temporary work array.
    !   
    IF (nsmax.GT.0)  THEN
       DO n = 1,naz
          DO l=lev1,lev2
             DO i=il1,il2
                smoothr1(i,l) = m_alpha(i,l,n)
             END DO
          END DO
          CALL vert_smooth (smoothr1,smoothr2, smco, nsmax,  &
                            il1, il2, lev1, lev2, nlons, nlevs )
          DO l=lev1,lev2
             DO i=il1,il2
                m_alpha(i,l,n) = smoothr1(i,l)
             END DO
          END DO
       END DO
       CALL vert_smooth ( sigma_t, smoothr2, smco, nsmax, &
                          il1, il2, lev1, lev2, nlons, nlevs )
    END IF
    !
    !  calculate zonal and meridional components of the
    !  momentum flux and drag.
    !
    CALL hines_flux ( flux_u, flux_v, flux, drag_u, drag_v,       &
                      alt, density, densb,                        &
                      m_alpha,  ak_alpha, k_alpha,                &
                      m_min, slope, naz,                          &
                      il1, il2, lev1, lev2, nlons, nlevs, nazmth, &
                      lorms )
    !
    !  cutoff drag above alt_cutoff, using bvfb as temporary work array.
    !
    IF (icutoff.EQ.1)  THEN	
       CALL hines_exp ( drag_u, bvfb, alt, alt_cutoff,  &
                        il1, il2, lev1, lev2, nlons, nlevs )
       CALL hines_exp ( drag_v, bvfb, alt, alt_cutoff,  &
                        il1, il2, lev1, lev2, nlons, nlevs )
    END IF
    !
    !  print out various arrays for diagnostic purposes.
    !
    IF (iprint.EQ.1)  THEN
       CALL hines_print ( flux_u, flux_v, drag_u, drag_v, alt,     &
                          sigma_t, sigma_alpha, v_alpha, m_alpha,  &
                          1, 1, il1, il2, lev1, lev2,              &
                          naz, nlons, nlevs, nazmth)
    END IF
    !
    !  if not calculating heating rate and diffusion coefficient then finished.
    !
    IF (iheatcal.NE.1)  RETURN
    !
    !
    !  heating rate and diffusion coefficient.

    CALL hines_heat ( heat, diffco,                                 &
                      alt, bvfreq, density, sigma_t, sigma_alpha,   &
                      flux, visc_mol, kstar, f1, f2, f3, f5, f6,    &
                      naz, il1, il2, lev1, lev2, nlons, nlevs,      &
                      nazmth, losigma_t )
    !

    !  finished.
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_extro

  SUBROUTINE hines_wavnum ( m_alpha, sigma_t, sigma_alpha, ak_alpha,     &
                            mmin_alpha, losigma_t,                       &
                            v_alpha, visc_mol, density, densb,           &
                            bvfreq, bvfb, rms_wind, anis, lorms,         &
                            sigsqmcw, sigmatm,                           &
                            il1, il2, levtop, levbot, nlons, nlevs, nazmth)
    !
    !  this routine calculates the cutoff vertical wavenumber and velocity
    !  variances on a longitude by altitude grid for the hines' doppler 
    !  spread gravity wave drag parameterization scheme.
    !  note: (1) only values of four or eight can be used for # azimuths (naz).
    !        (2) only values of 1.0, 1.5 or 2.0 can be used for slope (slope). 
    !        (3) if m_min not zero, only slope=1. can be used. 
    !
    !  aug. 10/95 - c. mclandress
    !  2000-2001  - m. charron
    !  2002       - e. manzini
    !
    !  output arguements:
    !
    !     * m_alpha      = cutoff wavenumber at each azimuth (1/m).
    !     * sigma_t      = total rms horizontal wind (m/s).
    !     * sigma_alpha  = total rms wind in each azimuth (m/s).
    !     * ak_alpha     = spectral amplitude factor at each azimuth 
    !     *                (i.e.,{ajkj}) in m^4/s^2.
    !     * losigma_t    =  .true. for total sigma not zero
    !     * mmin_alpha   = minimum value of cutoff wavenumber.


    !  input arguements:
    !
    !     * v_alpha  = wind component at each azimuth (m/s). 
    !     * visc_mol = molecular viscosity (m^2/s)
    !     * density  = background density (kg/m^3).
    !     * densb    = background density at model bottom (kg/m^3).
    !     * bvfreq   = background brunt vassala frequency (radians/sec).
    !     * bvfb     = background brunt vassala frequency at model bottom.
    !     * rms_wind = root mean square gravity wave wind at lowest level (m/s).
    !     * anis      = anisotropy factor (sum over azimuths is one)
    !     * lorms       = .true. for drag computation at lowest level 

    !     * levbot   = index of lowest vertical level.
    !     * levtop   = index of highest vertical level 
    !     *            (note: if levtop < levbot then level index 
    !     *             increases from top down).
    !     * il1      = first longitudinal index to use (il1 >= 1).
    !     * il2      = last longitudinal index to use (il1 <= il2 <= nlons).

    !     * nlons    = number of longitudes.
    !     * nlevs    = number of vertical levels.
    !     * nazmth   = azimuthal array dimension (nazmth >= naz).
    !

    !
    !  work arrays:
    !
    !     * i_alpha    = hines' integral at a single level.
    !
    !     * do_alpha = .true. for the azimuths and longitudes for
    !                   which to continue to compute the drag above  
    !                   the lowest level


    USE mo_doctor,       ONLY: nerr
    USE mo_exception,    ONLY: finish 
    USE mo_gwspectrum,   ONLY: kstar, m_min, slope, f1, f2, f3, naz

    IMPLICIT NONE

    INTEGER :: il1, il2, levtop, levbot, nlons, nlevs, nazmth
    REAL(dp):: m_alpha(nlons,nlevs,nazmth)
    REAL(dp):: sigma_alpha(nlons,nlevs,nazmth)
    REAL(dp):: sigalpmc(nlons,nlevs,nazmth)
    REAL(dp):: sigsqh_alpha(nlons,nlevs,nazmth)
    REAL(dp):: sigsqmcw(nlons,nlevs,nazmth)
    REAL(dp):: sigma_t(nlons,nlevs)
    REAL(dp):: sigmatm(nlons,nlevs)
    REAL(dp):: ak_alpha(nlons,nazmth)
    REAL(dp):: v_alpha(nlons,nlevs,nazmth)
    REAL(dp):: visc_mol(nlons,nlevs)
    REAL(dp):: f2mod(nlons,nlevs)
    REAL(dp):: density(nlons,nlevs),  densb(nlons)
    REAL(dp):: bvfreq(nlons,nlevs),   bvfb(nlons),  rms_wind(nlons)
    REAL(dp):: anis(nlons,nazmth) 
    REAL(dp):: i_alpha(nlons,nazmth), mmin_alpha(nlons,nazmth)

    LOGICAL :: lorms(nlons), losigma_t(nlons,nlevs), do_alpha(nlons,nazmth)
    !
    ! internal variables.
    !
    INTEGER :: i, l, n, istart, lend, lincr, lbelow

    REAL(dp):: m_sub_m_turb, m_sub_m_mol, m_trial, mmsq
    REAL(dp):: visc, visc_min, sp1, f2mfac

    REAL(dp):: n_over_m(nlons), sigfac(nlons)
    !-----------------------------------------------------------------------     
    !

    visc_min = 1.e-10_dp

    sp1 = slope + 1.0_dp
    mmsq = m_min**2

    !
    !  indices of levels to process.
    !
    IF (levbot > levtop)  THEN
       istart = levbot - 1     
       lend   = levtop         
       lincr  = -1
    ELSE
       WRITE (nerr,*) ' Error: level index not increasing downward '
       CALL finish('hines_wavnum','Run terminated')
    END IF


    !   initialize logical flags and arrays
    DO l=1,nlevs
       losigma_t(:,l) = lorms(:)
    ENDDO
    DO n=1,nazmth
       do_alpha(:,n) = lorms(:)
    ENDDO

    sigsqh_alpha(:,:,:) = 0.0_dp
    i_alpha(:,:) = 0.0_dp

    !
    ! calculate azimuthal variances at bottom level using anisotropy factor
    !
    DO n = 1,naz
       DO i = il1,il2
          sigsqh_alpha(i,levbot,n) = anis(i,n)* rms_wind(i)**2
       END DO
    END DO
    !
    !  velocity variances at bottom level.
    !
    CALL hines_sigma ( sigma_t, sigma_alpha,     &
                       sigsqh_alpha, naz, levbot,     &
                       il1, il2, nlons, nlevs, nazmth)

    CALL hines_sigma ( sigmatm, sigalpmc,     &
                       sigsqmcw, naz, levbot,     &
                       il1, il2, nlons, nlevs, nazmth)
    !
    !  calculate cutoff wavenumber and spectral amplitude factor 
    !  at bottom level where it is assumed that background winds vanish
    !  and also initialize minimum value of cutoff wavnumber.
    !
    IF ( ABS(slope-1.0_dp) < EPSILON(1.0_dp) ) THEN
       DO n = 1,naz
          DO i = il1,il2
             IF (lorms(i)) THEN
                m_alpha(i,levbot,n) =  bvfb(i) /    &
                     ( f1 * sigma_alpha(i,levbot,n)    &
                     + f2 * sigma_t(i,levbot) )
                ak_alpha(i,n)   = 2.0_dp * sigsqh_alpha(i,levbot,n)    &
                     / ( m_alpha(i,levbot,n)**2 - mmsq )
                mmin_alpha(i,n) = m_alpha(i,levbot,n)
             ENDIF
          END DO
       END DO
    ELSE
       DO n = 1,naz
          DO i = il1,il2
             IF (lorms(i)) THEN
                m_alpha(i,levbot,n) =  bvfb(i) /    & 
                     ( f1 * sigma_alpha(i,levbot,n)    & 
                     + f2 * sigma_t(i,levbot) )
                ak_alpha(i,n)   = sigsqh_alpha(i,levbot,n)    & 
                     / ( m_alpha(i,levbot,n)**sp1 / sp1 )
                mmin_alpha(i,n) = m_alpha(i,levbot,n)
             ENDIF
          END DO
       END DO
    ENDIF
    !
    !  calculate quantities from the bottom upwards, 
    !  starting one level above bottom.
    !
    DO l = istart,lend,lincr
       !
       !  level beneath present level.
       !
       lbelow = l - lincr 
       !
       !  calculate n/m_m where m_m is maximum permissible value of the vertical
       !  wavenumber (i.e., m > m_m are obliterated) and n is buoyancy frequency.
       !  m_m is taken as the smaller of the instability-induced 
       !  wavenumber (m_sub_m_turb) and that imposed by molecular viscosity
       !  (m_sub_m_mol). since variance at this level is not yet known
       !  use value at level below.
       !

!OCL VCT(MASK)
       DO i = il1,il2
          IF (losigma_t(i,lbelow))   THEN

             f2mfac=sigmatm(i,lbelow)**2
             f2mod(i,lbelow) =1.0_dp+ 2.0_dp*f2mfac  &
                  / ( f2mfac+sigma_t(i,lbelow)**2 )

             visc = MAX ( visc_mol(i,l), visc_min )
             m_sub_m_turb = bvfreq(i,l)   &
                  / ( f2 *f2mod(i,lbelow)*sigma_t(i,lbelow))
             m_sub_m_mol = (bvfreq(i,l)*kstar/visc)**0.33333333_dp/f3

             IF (m_sub_m_turb < m_sub_m_mol)  THEN
                n_over_m(i) = f2 *f2mod(i,lbelow)*sigma_t(i,lbelow)
             ELSE
                n_over_m(i) = bvfreq(i,l) / m_sub_m_mol 
             END IF

          ENDIF
       END DO
       !
       !  calculate cutoff wavenumber at this level.
       !
       DO n = 1,naz
          DO i = il1,il2
             IF ( do_alpha(i,n) .AND. losigma_t(i,lbelow) ) THEN
                !
                !  calculate trial value (variance at this level is not yet known:
                !  use value at level below). if trial value negative or larger 
                !  minimum value (not permitted) then set it to minimum value. 
                !                                                                      
                m_trial = bvfb(i) / ( f1 * ( sigma_alpha(i,lbelow,n)+   & 
                     sigalpmc(i,lbelow,n)) + n_over_m(i) + v_alpha(i,l,n) )

                IF (m_trial <= 0.0_dp .OR. m_trial > mmin_alpha(i,n))  THEN
                   m_trial = mmin_alpha(i,n)
                END IF
                m_alpha(i,l,n) = m_trial

                !  do not permit cutoff wavenumber to be less than minimum  value.

                IF (m_alpha(i,l,n) < m_min) THEN
                   m_alpha(i,l,n) = m_min
                ENDIF
                !
                !  reset minimum value of cutoff wavenumber if necessary.
                !
                IF (m_alpha(i,l,n) < mmin_alpha(i,n))  THEN
                   mmin_alpha(i,n) = m_alpha(i,l,n)
                END IF
             ELSE

                m_alpha(i,l,n) = m_min

             ENDIF
          END DO
       END DO
       !
       !  calculate the hines integral at this level.
       !
       CALL hines_intgrl ( i_alpha,                                     &
                           v_alpha, m_alpha, bvfb, m_min, slope, naz,   &
                           l, il1, il2, nlons, nlevs, nazmth,           &
                           lorms, do_alpha )

       !
       !  calculate the velocity variances at this level.
       !
       DO i = il1,il2
          sigfac(i) = densb(i) / density(i,l) * bvfreq(i,l) / bvfb(i) 
       END DO
       DO n = 1,naz
          DO i = il1,il2
             sigsqh_alpha(i,l,n) = sigfac(i) * ak_alpha(i,n) * i_alpha(i,n)
          END DO
       END DO
       CALL hines_sigma ( sigma_t, sigma_alpha, sigsqh_alpha, naz, l, &
                          il1, il2, nlons, nlevs, nazmth )

       CALL hines_sigma ( sigmatm, sigalpmc, sigsqmcw, naz, l,   &
                          il1, il2, nlons, nlevs, nazmth )

       !
       !  if total rms wind zero (no more drag) then set drag to false
       !
       DO i=il1,il2
          IF ( sigma_t(i,l) < EPSILON(1.0_dp)) THEN
             losigma_t(i,l) = .FALSE.
          ENDIF
       ENDDO
       !
       !  end of level loop.
       !
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_wavnum

  SUBROUTINE hines_wind (v_alpha,vel_u,vel_v,  &
                         naz,il1,il2,lev1,lev2,nlons,nlevs,nazmth)
    !
    !  this routine calculates the azimuthal horizontal background wind components 
    !  on a longitude by altitude grid for the case of 4 or 8 azimuths for
    !  the hines' doppler spread gwd parameterization scheme.
    !
    !  aug. 7/95 - c. mclandress
    !
    !  output arguement:
    !
    !     * v_alpha   = background wind component at each azimuth (m/s). 
    !     *             (note: first azimuth is in eastward direction
    !     *              and rotate in counterclockwise direction.)
    !
    !  input arguements:
    !
    !     * vel_u     = background zonal wind component (m/s).
    !     * vel_v     = background meridional wind component (m/s).
    !     * naz       = actual number of horizontal azimuths used (must be 4 or 8).
    !     * il1       = first longitudinal index to use (il1 >= 1).
    !     * il2       = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1      = first altitude level to use (lev1 >=1). 
    !     * lev2      = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons     = number of longitudes.
    !     * nlevs     = number of vertical levels.
    !     * nazmth    = azimuthal array dimension (nazmth >= naz).
    !
    !  constants in data statements.
    !
    !     * cos45 = cosine of 45 degrees. 		
    !     * umin  = minimum allowable value for zonal or meridional 
    !     *         wind component (m/s).
    !
    !  subroutine arguements.

    INTEGER  :: naz, il1, il2, lev1, lev2
    INTEGER  :: nlons, nlevs, nazmth
    REAL(dp) :: v_alpha(nlons,nlevs,nazmth)
    REAL(dp) ::  vel_u(nlons,nlevs), vel_v(nlons,nlevs)
    !
    !  internal variables.
    !
    INTEGER  :: i, l
    REAL(dp) :: u, v, cos45, umin
    !----------------------------------------------------------------------- 

    cos45 = 0.7071068_dp 
    umin  = 0.001_dp 

    SELECT CASE (naz) 
       !
    CASE(4)  !  case with 4 azimuths.

       DO l = lev1,lev2
!CDIR NODEP
          DO i = il1,il2
             u = vel_u(i,l)
             v = vel_v(i,l)
             IF (ABS(u) .LT. umin)  u = umin 
             IF (ABS(v) .LT. umin)  v = umin 
             v_alpha(i,l,1) = u 
             v_alpha(i,l,2) = v
             v_alpha(i,l,3) = - u
             v_alpha(i,l,4) = - v
          END DO
       END DO
       !
    CASE (8)   !  case with 8 azimuths.
       !
       DO l = lev1,lev2
!CDIR NODEP
          DO i = il1,il2
             u = vel_u(i,l)
             v = vel_v(i,l)
             IF (ABS(u) .LT. umin)  u = umin  
             IF (ABS(v) .LT. umin)  v = umin  
             v_alpha(i,l,1) = u 
             v_alpha(i,l,2) = cos45 * ( v + u )
             v_alpha(i,l,3) = v
             v_alpha(i,l,4) = cos45 * ( v - u )
             v_alpha(i,l,5) = - u
             v_alpha(i,l,6) = - v_alpha(i,l,2)
             v_alpha(i,l,7) = - v
             v_alpha(i,l,8) = - v_alpha(i,l,4)
          END DO
       END DO
       !
    END SELECT
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_wind

  SUBROUTINE hines_flux ( flux_u, flux_v, flux, drag_u, drag_v,        & 
                          alt, density, densb,                         &
                          m_alpha, ak_alpha, k_alpha,                  &
                          m_min, slope, naz,                           &
                          il1, il2, lev1, lev2, nlons, nlevs, nazmth,  &
                          lorms )
    !
    !  calculate zonal and meridional components of the vertical flux 
    !  of horizontal momentum and corresponding wave drag (force per unit mass)
    !  on a longitude by altitude grid for the hines' doppler spread 
    !  gwd parameterization scheme.
    !  note: only 4 or 8 azimuths can be used.
    !
    !  aug. 6/95 - c. mclandress
    !       2001 - m. charron
    !
    !  output arguements:
    !
    !     * flux_u = zonal component of vertical momentum flux (pascals)
    !     * flux_v = meridional component of vertical momentum flux (pascals)
    !     * drag_u = zonal component of drag (m/s^2).
    !     * drag_v = meridional component of drag (m/s^2).
    !
    !  input arguements:
    !
    !     * alt       = altitudes (m).
    !     * density   = background density (kg/m^3).
    !     * densb     = background density at bottom level (kg/m^3).
    !     * m_alpha   = cutoff vertical wavenumber (1/m).
    !     * ak_alpha  = spectral amplitude factor (i.e., {ajkj} in m^4/s^2).
    !     * k_alpha   = horizontal wavenumber (1/m).
    !     * slope     = slope of incident vertical wavenumber spectrum.
    !     * m_min     = minimum allowable cutoff wavenumber (1/m)
    !     *             for spectral slope of one.
    !     * naz       = actual number of horizontal azimuths used (must be 4 or 8).
    !     * il1       = first longitudinal index to use (il1 >= 1).
    !     * il2       = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1      = first altitude level to use (lev1 >=1). 
    !     * lev2      = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons     = number of longitudes.
    !     * nlevs     = number of vertical levels.
    !     * nazmth    = azimuthal array dimension (nazmth >= naz).
    !     * lorms     = .true. for drag computation (column selector)
    !
    !  constant in data statement.
    !
    !     * cos45 = cosine of 45 degrees. 		
    !
    !  subroutine arguements.
    !

    INTEGER :: naz, il1, il2, lev1, lev2, lev2p
    INTEGER :: nlons, nlevs, nazmth
    REAL(dp) ::  slope, m_min
    REAL(dp) ::  flux_u(nlons,nlevs), flux_v(nlons,nlevs)
    REAL(dp) ::  flux(nlons,nlevs,nazmth)
    REAL(dp) ::  drag_u(nlons,nlevs), drag_v(nlons,nlevs)
    REAL(dp) ::  alt(nlons,nlevs),    density(nlons,nlevs), densb(nlons)
    REAL(dp) ::  m_alpha(nlons,nlevs,nazmth)
    REAL(dp) ::  ak_alpha(nlons,nazmth), k_alpha(nlons,nazmth)

    LOGICAL :: lorms(nlons)
    !
    !  internal variables.
    !
    INTEGER :: i, l, lev1p, lev2m, k
    REAL(dp) ::  cos45, dendz, dendz2
    !-----------------------------------------------------------------------
    cos45 = 0.7071068_dp   
    !
    lev1p = lev1 + 1
    lev2m = lev2 - 1
    lev2p = lev2 + 1
    !
    !  sum over azimuths for case where slope = 1.
    !
    IF ( ABS(slope-1.0_dp) .LT. EPSILON(1.0_dp) )  THEN
       !
       !  case with 4 azimuths.
       !
       IF (naz.EQ.4)  THEN
          DO l = lev1,lev2
             DO i = il1,il2
                flux(i,l,:) = ak_alpha(i,:)*k_alpha(i,:)*(m_alpha(i,l,:)-m_min)
                flux_u(i,l) = flux(i,l,1) - flux(i,l,3)
                flux_v(i,l) = flux(i,l,2) - flux(i,l,4)
             END DO
          END DO
       END IF
       !
       !  case with 8 azimuths.
       !
       IF (naz.EQ.8)  THEN
          DO l = lev1,lev2
             DO k = 1, nazmth
                DO i = il1,il2
                  flux(i,l,k) = ak_alpha(i,k)*k_alpha(i,k)*(m_alpha(i,l,k)-m_min)
                END DO
             END DO
             DO i = il1,il2
                flux_u(i,l) = flux(i,l,1) - flux(i,l,5) + cos45 *     &
                     ( flux(i,l,2) - flux(i,l,4) - flux(i,l,6) + flux(i,l,8) )
                flux_v(i,l) = flux(i,l,3) - flux(i,l,7) + cos45 *     &
                     ( flux(i,l,2) + flux(i,l,4) - flux(i,l,6) - flux(i,l,8) )
             END DO
          END DO
       END IF

    END IF
    !
    !  sum over azimuths for case where slope not equal to 1.
    !
    IF ( ABS(slope-1.0_dp) .GT. EPSILON(1.0_dp) )  THEN
       !
       !  case with 4 azimuths.
       !
       IF (naz.EQ.4)  THEN
          DO l = lev1,lev2
             DO i = il1,il2
                flux(i,l,:) = ak_alpha(i,:)*k_alpha(i,:)*m_alpha(i,l,:)**slope
                flux_u(i,l) = flux(i,l,1) - flux(i,l,3)
                flux_v(i,l) = flux(i,l,2) - flux(i,l,4)
             END DO
          END DO
       END IF
       !
       !  case with 8 azimuths.
       !
       IF (naz.EQ.8)  THEN
          DO l = lev1,lev2
             DO k = 1, nazmth
                DO i = il1,il2
                  flux(i,l,k) = ak_alpha(i,k)*k_alpha(i,k)*m_alpha(i,l,k)**slope
                END DO
             END DO
             DO i = il1,il2
                flux_u(i,l) = flux(i,l,1) - flux(i,l,5) + cos45 *     &
                     ( flux(i,l,2) - flux(i,l,4) - flux(i,l,6) + flux(i,l,8) )
                flux_v(i,l) = flux(i,l,3) - flux(i,l,7) + cos45 *     &
                     ( flux(i,l,2) + flux(i,l,4) - flux(i,l,6) - flux(i,l,8) )
             END DO
          END DO
       END IF

    END IF
    !
    !  calculate flux from sum.
    !
    DO l = lev1,lev2
       DO i = il1,il2
          flux_u(i,l) = flux_u(i,l) * densb(i) / slope
          flux_v(i,l) = flux_v(i,l) * densb(i) / slope
       END DO
       DO k = 1, nazmth
          DO i = il1,il2
            flux(i,l,k) = flux(i,l,k) * densb(i) / slope
          END DO
       END DO
    END DO
    !
    !  calculate drag at intermediate levels
    !      
    DO l = lev1p,lev2m
       DO i = il1,il2
          IF (lorms(i)) THEN
             dendz2 = density(i,l) * ( alt(i,l-1) - alt(i,l) ) 
             drag_u(i,l) = - ( flux_u(i,l-1) - flux_u(i,l) ) / dendz2
             drag_v(i,l) = - ( flux_v(i,l-1) - flux_v(i,l) ) / dendz2       
          ENDIF
       END DO
    END DO
    !
    !  calculate drag at intermediate levels using centered differences (not used) 
    !ccc       dendz2 = density(i,l) * ( alt(i,l+1) - alt(i,l-1) )
    !ccc       drag_u(i,l) = - ( flux_u(i,l+1) - flux_u(i,l-1) ) / dendz2
    !ccc       drag_v(i,l) = - ( flux_v(i,l+1) - flux_v(i,l-1) ) / dendz2


    !  drag at first and last levels using one-side differences.
    ! 
    DO i = il1,il2
       IF (lorms(i)) THEN
          dendz = density(i,lev1) * ( alt(i,lev1) - alt(i,lev1p) ) 
          drag_u(i,lev1) =  flux_u(i,lev1)  / dendz
          drag_v(i,lev1) =  flux_v(i,lev1)  / dendz
       ENDIF
    END DO
    DO i = il1,il2
       IF (lorms(i)) THEN
          dendz = density(i,lev2) * ( alt(i,lev2m) - alt(i,lev2) )
          drag_u(i,lev2) = - ( flux_u(i,lev2m) - flux_u(i,lev2) ) / dendz
          drag_v(i,lev2) = - ( flux_v(i,lev2m) - flux_v(i,lev2) ) / dendz
       ENDIF
    END DO
    IF (nlevs .GT. lev2) THEN
       DO i = il1,il2
          IF (lorms(i)) THEN
             dendz = density(i,lev2p) * ( alt(i,lev2) - alt(i,lev2p) )
             drag_u(i,lev2p) = -  flux_u(i,lev2)  / dendz
             drag_v(i,lev2p) = - flux_v(i,lev2)  / dendz
          ENDIF
       END DO
    ENDIF
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_flux

  SUBROUTINE hines_heat ( heat, diffco,                                 &
                          alt, bvfreq, density, sigma_t, sigma_alpha,   &
                          flux, visc_mol, kstar, f1, f2, f3, f5, f6,    &
                          naz, il1, il2, lev1, lev2, nlons, nlevs,      & 
                          nazmth, losigma_t )
    !
    !  this routine calculates the gravity wave induced heating and 
    !  diffusion coefficient on a longitude by altitude grid for  
    !  the hines' doppler spread gravity wave drag parameterization scheme.
    !
    !  This routine can be used for nonzero minimum cutoff wavenumber (m_min)
    !  only in the case of spectral slope=1, in which case m_min is not needed
    !  since its vertical derivative is zero.
    !
    !  aug. 6/95 - c. mclandress
    !  2001      - m. charron  
    !
    !  output arguements:
    !
    !     * heat   = gravity wave heating (k/sec).
    !     * diffco = diffusion coefficient (m^2/sec)
    !
    !  input arguements:
    !
    !
    !     * bvfreq      = background brunt vassala frequency (rad/sec).
    !     * density     = background density (kg/m^3).
    !
    !     * sigma_t     = total rms horizontal wind (m/s).
    !
    !     * visc_mol    = molecular viscosity (m^2/s).
    !     * kstar       = typical gravity wave horizontal wavenumber (1/m).
    !     * slope       = slope of incident vertical wavenumber spectrum.
    !     * f1,f2,f3,f5,f6 = hines's fudge factors.
    !
    !     * il1         = first longitudinal index to use (il1 >= 1).
    !     * il2         = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1        = first altitude level to use (lev1 >=1). 
    !     * lev2        = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons       = number of longitudes.
    !     * nlevs       = number of vertical levels.
    !     * nazmth      = azimuthal array dimension (nazmth >= naz).
    !     * losigma_t   = .true. for total sigma not zero
    !
    USE mo_constants,    ONLY: cpd

    IMPLICIT NONE

    INTEGER ::  naz, il1, il2, lev1, lev2, nlons, nlevs, nazmth
    REAL(dp)::  kstar, f1, f2, f3, f5, f6
    REAL(dp)::  heat(nlons,nlevs), diffco(nlons,nlevs)
    REAL(dp)::  alt(nlons,nlevs), bvfreq(nlons,nlevs), density(nlons,nlevs) 
    REAL(dp)::  sigma_t(nlons,nlevs),  sigma_alpha(nlons,nlevs,nazmth)
    REAL(dp)::  flux(nlons,nlevs,nazmth), visc_mol(nlons,nlevs)
    LOGICAL ::  losigma_t(nlons,nlevs)
    !
    ! internal variables.
    !
    INTEGER  :: i, l, n, lev1p, lev2m
    REAL(dp) :: m_sub_m_turb, m_sub_m_mol, m_sub_m, heatng, dendz2
    REAL(dp) :: visc, visc_min

    REAL(dp) :: dfdz(nlons,nlevs,nazmth)
    !-----------------------------------------------------------------------   

    visc_min = 1.e-10_dp 

    lev1p = lev1 + 1
    lev2m = lev2 - 1   

    DO l = lev1p,lev2m
       DO i = il1,il2
          IF (losigma_t(i,l)) THEN
             dendz2 = density(i,l) * ( alt(i,l-1) - alt(i,l) )
             visc    = MAX ( visc_mol(i,l), visc_min )
             m_sub_m_turb = bvfreq(i,l) / ( f2 * sigma_t(i,l) )
             m_sub_m_mol  = (bvfreq(i,l)*kstar/visc)**0.33333333_dp/f3
             m_sub_m      = MIN ( m_sub_m_turb, m_sub_m_mol )
!CDIR UNROLL=8
             dfdz(i,l,:) = ( flux(i,l-1,:) - flux(i,l,:) ) / dendz2 &
                  & * ( f1*sigma_alpha(i,l,:) + bvfreq(i,l)/m_sub_m )
          ENDIF
       END DO
    END DO

    DO i = il1,il2
       IF (losigma_t(i,lev1)) THEN
          dendz2 = density(i,lev1) * ( alt(i,lev1) - alt(i,lev1p) )
          visc    = MAX ( visc_mol(i,lev1), visc_min )
          m_sub_m_turb = bvfreq(i,lev1) / ( f2 * sigma_t(i,lev1) )
          m_sub_m_mol  = (bvfreq(i,lev1)*kstar/visc)**0.33333333_dp/f3
          m_sub_m      = MIN ( m_sub_m_turb, m_sub_m_mol )
!CDIR UNROLL=8
          dfdz(i,lev1,:) = -flux(i,lev1,:) / dendz2 &
               & * ( f1*sigma_alpha(i,lev1,:) + bvfreq(i,lev1)/m_sub_m )
       ENDIF
    END DO

    DO i = il1,il2
       IF (losigma_t(i,lev2)) THEN
          dendz2 = density(i,lev2) * ( alt(i,lev2m) - alt(i,lev2) )
          visc    = MAX ( visc_mol(i,lev2), visc_min )
          m_sub_m_turb = bvfreq(i,lev2) / ( f2 * sigma_t(i,lev2) )
          m_sub_m_mol  = (bvfreq(i,lev2)*kstar/visc)**0.33333333_dp/f3
          m_sub_m      = MIN ( m_sub_m_turb, m_sub_m_mol )
!CDIR UNROLL=8
          dfdz(i,lev2,:) = ( flux(i,lev2m,:) - flux(i,lev2,:) ) / dendz2 &
               & * ( f1*sigma_alpha(i,lev2,:) + bvfreq(i,lev2)/m_sub_m )
       ENDIF
    END DO
    !
    !  heating and diffusion.

    !
    !  maximum permissible value of cutoff wavenumber is the smaller 
    !  of the instability-induced wavenumber (m_sub_m_turb) and 
    !  that imposed by molecular viscosity (m_sub_m_mol).
    !
    !
    DO l = lev1,lev2
       DO i = il1,il2
          IF (losigma_t(i,l)) THEN
             visc    = MAX ( visc_mol(i,l), visc_min )
             m_sub_m_turb = bvfreq(i,l) / ( f2 * sigma_t(i,l) )
             m_sub_m_mol  = (bvfreq(i,l)*kstar/visc)**0.33333333_dp/f3
             m_sub_m      = MIN ( m_sub_m_turb, m_sub_m_mol )
             heatng = 0.0_dp
             DO n=1,naz
                heatng = heatng - f5 * dfdz(i,l,n)
             ENDDO
             diffco(i,l) = f6 * heatng**0.33333333_dp / m_sub_m**1.33333333_dp
             heat(i,l)   = heatng / cpd
          ENDIF
       END DO
    END DO

    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_heat

  SUBROUTINE hines_sigma (sigma_t,sigma_alpha,sigsqh_alpha,  &
                          naz,lev,il1,il2,nlons,nlevs,nazmth)
    !
    !  this routine calculates the total rms and azimuthal rms horizontal 
    !  velocities at a given level on a longitude by altitude grid for 
    !  the hines' doppler spread gwd parameterization scheme.
    !  note: only four or eight azimuths can be used.
    !
    !  aug. 7/95 - c. mclandress
    !
    !  output arguements:
    !
    !     * sigma_t      = total rms horizontal wind (m/s).
    !     * sigma_alpha  = total rms wind in each azimuth (m/s).
    !
    !  input arguements:
    !
    !     * sigsqh_alpha = portion of wind variance from waves having wave
    !     *                normals in the alpha azimuth (m/s).
    !     * naz       = actual number of horizontal azimuths used (must be 4 or 8).
    !     * lev       = altitude level to process.
    !     * il1       = first longitudinal index to use (il1 >= 1).
    !     * il2       = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * nlons     = number of longitudes.
    !     * nlevs     = number of vertical levels.
    !     * nazmth    = azimuthal array dimension (nazmth >= naz).
    !
    !  subroutine arguements.
    !

    INTEGER  :: lev, naz, il1, il2
    INTEGER  :: nlons, nlevs, nazmth
    REAL(dp) ::  sigma_t(nlons,nlevs)
    REAL(dp) :: sigma_alpha(nlons,nlevs,nazmth)
    REAL(dp) :: sigsqh_alpha(nlons,nlevs,nazmth)
    !
    !  internal variables.
    !
    INTEGER  :: i, n
    REAL(dp) :: sum_even, sum_odd 
    !-----------------------------------------------------------------------     
    !
    !  calculate azimuthal rms velocity for the 4 azimuth case.
    !
    IF (naz.EQ.4)  THEN
!CDIR NODEP
       DO i = il1,il2
          sigma_alpha(i,lev,1) = SQRT(sigsqh_alpha(i,lev,1)+sigsqh_alpha(i,lev,3))
          sigma_alpha(i,lev,2) = SQRT(sigsqh_alpha(i,lev,2)+sigsqh_alpha(i,lev,4))
          sigma_alpha(i,lev,3) = sigma_alpha(i,lev,1)
          sigma_alpha(i,lev,4) = sigma_alpha(i,lev,2)
       END DO
    END IF
    !
    !  calculate azimuthal rms velocity for the 8 azimuth case.
    !
    IF (naz.EQ.8)  THEN
!CDIR NODEP
       DO i = il1,il2
          sum_odd  = ( sigsqh_alpha(i,lev,1) + sigsqh_alpha(i,lev,3)   &
               + sigsqh_alpha(i,lev,5) + sigsqh_alpha(i,lev,7) ) * 0.5_dp
          sum_even = ( sigsqh_alpha(i,lev,2) + sigsqh_alpha(i,lev,4)   &
               + sigsqh_alpha(i,lev,6) + sigsqh_alpha(i,lev,8) ) * 0.5_dp
          sigma_alpha(i,lev,1) = SQRT( sigsqh_alpha(i,lev,1)   &
               + sigsqh_alpha(i,lev,5) + sum_even )
          sigma_alpha(i,lev,2) = SQRT( sigsqh_alpha(i,lev,2)   &
               + sigsqh_alpha(i,lev,6) + sum_odd )
          sigma_alpha(i,lev,3) = SQRT( sigsqh_alpha(i,lev,3)   &
               + sigsqh_alpha(i,lev,7) + sum_even )
          sigma_alpha(i,lev,4) = SQRT( sigsqh_alpha(i,lev,4)   &
               + sigsqh_alpha(i,lev,8) + sum_odd )
          sigma_alpha(i,lev,5) = sigma_alpha(i,lev,1)
          sigma_alpha(i,lev,6) = sigma_alpha(i,lev,2)
          sigma_alpha(i,lev,7) = sigma_alpha(i,lev,3)
          sigma_alpha(i,lev,8) = sigma_alpha(i,lev,4)
       END DO
    END IF
    !
    !  calculate total rms velocity.
    !
    DO i = il1,il2
       sigma_t(i,lev) = 0.0_dp
    END DO
    DO n = 1,naz
       DO i = il1,il2
          sigma_t(i,lev) = sigma_t(i,lev) + sigsqh_alpha(i,lev,n)
       END DO
    END DO
    DO i = il1,il2
       sigma_t(i,lev) = SQRT( sigma_t(i,lev) )
    END DO
    !
    !-----------------------------------------------------------------------     
  END SUBROUTINE hines_sigma

  SUBROUTINE hines_intgrl (i_alpha,                                     &
                           v_alpha, m_alpha, bvfb, m_min, slope, naz,   &
                           lev, il1, il2, nlons, nlevs, nazmth,         &
                           lorms, do_alpha)
    !
    !  this routine calculates the vertical wavenumber integral
    !  for a single vertical level at each azimuth on a longitude grid
    !  for the hines' doppler spread gwd parameterization scheme.
    !  note: (1) only spectral slopes of 1, 1.5 or 2 are permitted.
    !        (2) the integral is written in terms of the product qm
    !            which by construction is always less than 1. series
    !            solutions are used for small |qm| and analytical solutions
    !            for remaining values.
    !
    !  aug. 8/95 - c. mclandress
    !  2001      - m. charron
    !  2003      - l. kornblueh
    !
    !  output arguement:
    !
    !     * i_alpha = hines' integral.
    !
    !  input arguements:
    !
    !     * v_alpha = azimuthal wind component (m/s). 
    !     * m_alpha = azimuthal cutoff vertical wavenumber (1/m).
    !     * bvfb    = background brunt vassala frequency at model bottom.
    !     * m_min   = minimum allowable cutoff vertical wavenumber (1/m)
    !     *           for spectral slope of one.
    !     * slope   = slope of initial vertical wavenumber spectrum 
    !     *           (must use slope = 1., 1.5 or 2.)
    !     * naz     = actual number of horizontal azimuths used.
    !     * lev     = altitude level to process.
    !     * il1     = first longitudinal index to use (il1 >= 1).
    !     * il2     = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * nlons   = number of longitudes.
    !     * nlevs   = number of vertical levels.
    !     * nazmth  = azimuthal array dimension (nazmth >= naz).
    !     * lorms     = .true. for drag computation (column selector)

    !
    !  constants in data statements:
    !
    !     * qmin = minimum value of q_alpha (avoids indeterminant form of integral)
    !     * qm_min = minimum value of q_alpha * m_alpha (used to avoid numerical
    !     *          problems).
    !

    USE mo_doctor,     ONLY: nout, nerr
    USE mo_exception,  ONLY: finish 

    INTEGER  :: lev, naz, il1, il2, nlons, nlevs, nazmth
    REAL(dp) :: i_alpha(nlons,nazmth)
    REAL(dp) :: v_alpha(nlons,nlevs,nazmth)
    REAL(dp) :: m_alpha(nlons,nlevs,nazmth)
    REAL(dp) :: bvfb(nlons), rbvfb(nlons), slope, m_min

    LOGICAL :: lorms(nlons), do_alpha(nlons,nazmth)
    LOGICAL :: lerror(nlons)
    !
    !  internal variables.
    !
    INTEGER  :: i, n
    REAL(dp) :: q_alpha, qm, qmm, sqrtqm, q_min, qm_min

    !  variables for sparse vector optimization
    INTEGER :: ic, ixi(nlons*nazmth), ix, ixnaz(nlons*nazmth)

    !-----------------------------------------------------------------------
    !
    !  initialize local scalar and arrays

    q_min = 1.0_dp 
    qm_min = 0.01_dp 

    DO i = il1,il2
      rbvfb(i)=1.0_dp/bvfb(i)
    ENDDO
    !
    !  for integer value slope = 1.
    !
    IF ( ABS(slope-1.0_dp) < EPSILON(1.0_dp) )  THEN
       ic = 0
       DO n = 1,naz
!OCL VCT(MASK)
          DO i = il1,il2
             IF (lorms(i)) THEN
                ! 
                IF (m_alpha(i,lev,n) > m_min) THEN
                   q_alpha = v_alpha(i,lev,n) * rbvfb(i)
                   qm      = q_alpha * m_alpha(i,lev,n)
                   qmm     = q_alpha * m_min
                   !
                   !  if |qm| is small then use first 4 terms series of taylor
                   !  series expansion of integral in order to avoid 
                   !  indeterminate form of integral,
                   !  otherwise use analytical form of integral.
                   !
                   IF ( ABS(q_alpha) .LT. q_min .OR. ABS(qm).LT. qm_min)  THEN
                      ! taylor series expansion is a very rare event.
                      ! do sparse processing separately
                      ic = ic+1
                      ixi(ic) = i
                      ixnaz(ic) = n
                   ELSE
                      i_alpha(i,n) = - ( LOG(1.0_dp-qm) - LOG(1.0_dp-qmm) + qm - qmm) / q_alpha**2
                   END IF
                   !
                   !  If i_alpha negative due to round off error, set it to zero
                   !
                   i_alpha(i,n) = MAX( i_alpha(i,n) , 0.0_dp )
                ELSE
                   i_alpha(i,n) = 0.0_dp
                   do_alpha(i,n) = .FALSE.
                ENDIF
                !
             ENDIF
          END DO
       END DO
       ! taylor series expansion is a very rare event.
       ! do sparse processing here separately 
!CDIR NODEP
       DO ix =1, ic
          n = ixnaz(ix)
          i = ixi(ix)
          q_alpha = v_alpha(i,lev,n) * rbvfb(i)
          qm      = q_alpha * m_alpha(i,lev,n)
          qmm     = q_alpha * m_min
          IF ( ABS(q_alpha) < EPSILON(1.0_dp) )  THEN
             i_alpha(i,n) = ( m_alpha(i,lev,n)**2  - m_min**2 ) * 0.5_dp
          ELSE
             i_alpha(i,n) = ( qm**2 * 0.5_dp + qm**3  / 3.0_dp  + qm**4 * 0.25_dp + qm**5  * 0.2_dp   &
                           - qmm**2 * 0.5_dp - qmm**3 / 3.0_dp - qmm**4 * 0.25_dp - qmm**5 * 0.2_dp ) &
                 / q_alpha**2
          END IF
          i_alpha(i,n) = MAX( i_alpha(i,n) , 0.0_dp )
       END DO
    END IF
    !
    !  for integer value slope = 2.
    !
    IF ( ABS(slope-2.0_dp) < EPSILON(1.0_dp) )  THEN
       ic = 0
       DO n = 1,naz
!OCL VCT(MASK)
          DO i = il1,il2
             IF ( lorms(i) ) THEN
                !
                q_alpha = v_alpha(i,lev,n) * rbvfb(i)
                qm = q_alpha * m_alpha(i,lev,n)
                !
                !  if |qm| is small then use first 4 terms series of taylor 
                !  series expansion of integral in order to avoid 
                !  indeterminate form of integral,
                !  otherwise use analytical form of integral.
                !
                IF ( ABS(q_alpha) .LT. q_min .OR. ABS(qm) .LT. qm_min)  THEN  
                   ! taylor series expansion is a very rare event.
                   ! do sparse processing separately
                   ic = ic+1
                   ixi(ic) = i
                   ixnaz(ic) = n
                ELSE
                   i_alpha(i,n) = - ( LOG(1.0_dp-qm) + qm + qm**2 * 0.5_dp)    &
                        / q_alpha**3
                ENDIF
                !
             ENDIF
          END DO
       END DO
       ! taylor series expansion is a very rare event.
       ! do sparse processing here separately 
!CDIR NODEP
       DO ix = 1, ic
          n = ixnaz(ix)
          i = ixi(ix)
          q_alpha = v_alpha(i,lev,n) * rbvfb(i)
          qm = q_alpha * m_alpha(i,lev,n)
          IF ( ABS(q_alpha) < EPSILON(1.0_dp) )  THEN
             i_alpha(i,n) = m_alpha(i,lev,n)**3 / 3.0_dp
          ELSE
             i_alpha(i,n) = ( qm**3/3._dp + qm**4/4._dp + qm**5/5._dp    &
               + qm**6/6._dp ) / q_alpha**3
           END IF
       END DO
    END IF
    !
    !  for real value slope = 1.5
    !
    IF ( ABS(slope-1.5_dp) < EPSILON(1.0_dp) )  THEN
       ic = 0 
       DO n = 1,naz
!OCL VCT(MASK)
          DO i = il1,il2
             IF ( lorms(i) ) THEN
                !
                q_alpha = v_alpha(i,lev,n) * rbvfb(i)
                qm = q_alpha * m_alpha(i,lev,n)       
                !
                !  if |qm| is small then use first 4 terms series of taylor 
                !  series expansion of integral in order to avoid 
                !  indeterminate form of integral,
                !  otherwise use analytical form of integral.
                !
                IF (ABS(q_alpha) .LT. q_min .OR. ABS(qm) .LT. qm_min)  THEN  
                   ! taylor series expansion is a very rare event.
                   ! do sparse processing separately
                   ic = ic+1
                   ixi(ic) = i
                   ixnaz(ic) = n
                ELSE
                   qm     = ABS(qm)
                   sqrtqm = SQRT(qm)
                   IF (q_alpha .GE. 0.0_dp)  THEN
                      i_alpha(i,n) = ( LOG( (1.0_dp+sqrtqm)/(1.0_dp-sqrtqm) )  &
                           -2.0_dp*sqrtqm*(1.0_dp+qm/3.0_dp) ) / q_alpha**2.5_dp
                   ELSE
                      i_alpha(i,n) = 2.0_dp * ( ATAN(sqrtqm) + sqrtqm*(qm/3.0_dp-1.0_dp) ) &
                           / ABS(q_alpha)**2.5_dp
                   ENDIF
                ENDIF
                !
             ENDIF
          END DO
       END DO
       ! taylor series expansion is a very rare event.
       ! do sparse processing here separately 
!CDIR NODEP
       DO ix = 1, ic
          n = ixnaz(ix)
          i = ixi(ix)
          q_alpha = v_alpha(i,lev,n) * rbvfb(i)
          qm = q_alpha * m_alpha(i,lev,n)   
          IF ( ABS(q_alpha) < EPSILON(1.0_dp) )  THEN
             i_alpha(i,n) = m_alpha(i,lev,n)**2.5_dp / 2.5_dp
          ELSE
             i_alpha(i,n) = ( qm/2.5_dp + qm**2/3.5_dp   &
                  + qm**3/4.5_dp + qm**4/5.5_dp )   &
                  * m_alpha(i,lev,n)**1.5_dp / q_alpha
          END IF
       ENDDO
    END IF
    !
    !  if integral is negative (which in principal should not happen) then
    !  print a message and some info since execution will abort when calculating
    !  the variances.
    !
    DO n = 1,naz
       lerror(:) = .FALSE.
       DO i = il1, il2
          IF (i_alpha(i,n) < 0.0_dp)  THEN
             lerror(i) = .TRUE.
             EXIT 
          END IF
       END DO

       IF (ANY(lerror)) THEN
          WRITE (nout,*) 
          WRITE (nout,*) '******************************'
          WRITE (nout,*) 'hines integral i_alpha < 0 '
          WRITE (nout,*) '  longitude i=',i
          WRITE (nout,*) '  azimuth   n=',n
          WRITE (nout,*) '  level   lev=',lev
          WRITE (nout,*) '  i_alpha =',i_alpha(i,n)
          WRITE (nout,*) '  v_alpha =',v_alpha(i,lev,n)
          WRITE (nout,*) '  m_alpha =',m_alpha(i,lev,n)
          WRITE (nout,*) '  q_alpha =',v_alpha(i,lev,n)*rbvfb(i)
          WRITE (nout,*) '  qm      =',v_alpha(i,lev,n)*rbvfb(i)*m_alpha(i,lev,n)
          WRITE (nout,*) '******************************'
          WRITE (nerr,*) ' Error: Hines i_alpha integral is negative  '
          CALL finish(' hines_intgrl','Run terminated')
       END IF
  
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_intgrl

  SUBROUTINE hines_print (flux_u, flux_v, drag_u, drag_v, alt, sigma_t,    &
                          sigma_alpha, v_alpha, m_alpha,                   &
                          iu_print, iv_print,                              &
                          ilprt1, ilprt2, levprt1, levprt2, naz,           &
                          nlons, nlevs, nazmth)
    !
    !  print out altitude profiles of various quantities from
    !  hines' doppler spread gravity wave drag parameterization scheme.
    !  (note: only for naz = 4 or 8). 
    !
    !  aug. 8/95 - c. mclandress
    !
    !  input arguements:
    !
    !     * iu_print = 1 to print out values in east-west direction.
    !     * iv_print = 1 to print out values in north-south direction.

    !     * ilprt1   = first longitudinal index to print.
    !     * ilprt2   = last longitudinal index to print.
    !     * levprt1  = first altitude level to print.
    !     * levprt2  = last altitude level to print.
    !
    USE mo_doctor, ONLY: nout

    INTEGER  :: naz, ilprt1, ilprt2, levprt1, levprt2
    INTEGER  :: nlons, nlevs, nazmth
    INTEGER  :: iu_print, iv_print
    REAL(dp) :: flux_u(nlons,nlevs), flux_v(nlons,nlevs)
    REAL(dp) :: drag_u(nlons,nlevs), drag_v(nlons,nlevs)
    REAL(dp) :: alt(nlons,nlevs), sigma_t(nlons,nlevs)
    REAL(dp) :: sigma_alpha(nlons,nlevs,nazmth)
    REAL(dp) :: v_alpha(nlons,nlevs,nazmth), m_alpha(nlons,nlevs,nazmth)
    !
    !  internal variables.
    !
    INTEGER :: n_east, n_west, n_north, n_south
    INTEGER ::  i, l
    !-----------------------------------------------------------------------
    !
    !  azimuthal indices of cardinal directions.
    !
    n_east = 1
    IF (naz.EQ.4)  THEN
       n_west  = 3       
       n_north = 2
       n_south = 4       
    ELSE IF (naz.EQ.8)  THEN
       n_west  = 5       
       n_north = 3
       n_south = 7       
    END IF
    !
    !  print out values for range of longitudes.
    !
    DO i = ilprt1,ilprt2
       !
       !  print east-west wind, sigmas, cutoff wavenumbers, flux and drag.
       !
       IF (iu_print.EQ.1)  THEN
          WRITE (nout,*) 
          WRITE (nout,'(a,i3)') 'hines gw (east-west) at longitude i =',i
          WRITE (nout,6005) 
6005      FORMAT (15x,' u ',2x,'sig_e',2x,'sig_t',3x,'m_e',  &
               &            4x,'m_w',4x,'fluxu',5x,'gwdu')
          DO l = levprt1,levprt2
             WRITE (nout,6701) alt(i,l)/1.e3_dp, v_alpha(i,l,n_east),   &
                  &                          sigma_alpha(i,l,n_east), sigma_t(i,l),  &
                  &                          m_alpha(i,l,n_east)*1.e3_dp,   &
                  &                          m_alpha(i,l,n_west)*1.e3_dp,  &
                  &                          flux_u(i,l)*1.e5_dp, drag_u(i,l)*24.0_dp*3600.0_dp
          END DO
6701      FORMAT (' z=',f7.2,1x,3f7.1,2f7.3,f9.4,f9.3)
       END IF
       !
       !  print north-south winds, sigmas, cutoff wavenumbers, flux and drag.
       !
       IF (iv_print.EQ.1)  THEN
          WRITE(nout,*) 
          WRITE(nout,'(a,i3)') 'hines gw (north-south) at longitude i =',i
          WRITE(nout,6006) 
6006      FORMAT (15x,' v ',2x,'sig_n',2x,'sig_t',3x,'m_n',   &
               &            4x,'m_s',4x,'fluxv',5x,'gwdv')
          DO l = levprt1,levprt2
             WRITE (nout,6701) alt(i,l)/1.e3_dp, v_alpha(i,l,n_north),    &
                  &                          sigma_alpha(i,l,n_north), sigma_t(i,l),   &
                  &                          m_alpha(i,l,n_north)*1.e3_dp,    &
                  &                          m_alpha(i,l,n_south)*1.e3_dp,   &
                  &                          flux_v(i,l)*1.e5_dp, drag_v(i,l)*24.0_dp*3600.0_dp
          END DO
       END IF
       !
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_print

  SUBROUTINE hines_exp (darr, data_zmax, alt, alt_exp,     &
                        il1, il2, lev1, lev2, nlons, nlevs)
    !
    !  this routine exponentially damps a longitude by altitude array 
    !  of darr above a specified altitude.
    !
    !  aug. 13/95 - c. mclandress
    !
    !  output arguements:
    !
    !     * darr = modified data array.
    !
    !  input arguements:
    !
    !     * darr    = original data array.
    !     * alt     = altitudes.
    !     * alt_exp = altitude above which exponential decay applied.

    !     * il1     = first longitudinal index to use (il1 >= 1).
    !     * il2     = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1    = first altitude level to use (lev1 >=1). 
    !     * lev2    = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons   = number of longitudes.
    !     * nlevs   = number of vertical
    !
    !  input work arrays:
    !
    !     * data_zmax = data values just above altitude alt_exp.
    !

    USE mo_doctor,       ONLY: nerr
    USE mo_exception,    ONLY: finish 

    INTEGER  :: il1, il2, lev1, lev2, nlons, nlevs
    REAL(dp) :: alt_exp
    REAL(dp) :: darr(nlons,nlevs), data_zmax(nlons), alt(nlons,nlevs)
    !
    ! internal variables.
    !
    INTEGER  :: levbot, levtop, lincr, i, l
    REAL(dp) :: hscale
    !-----------------------------------------------------------------------     
    
    hscale = 5.e3_dp 

    !  index of lowest altitude level (bottom of drag calculation).
    !
    levbot = lev2
    levtop = lev1
    lincr  = 1
    IF (levbot > levtop)  THEN
       levbot = lev1
       levtop = lev2
       lincr  = -1
    ELSE
       WRITE (nerr,*) ' Error: level index not increasing downward '
       CALL finish('hines_exp','Run terminated')
    END IF
    !
    !  data values at first level above alt_exp.
    !
    DO i = il1,il2
       DO l = levtop,levbot,lincr
          IF (alt(i,l) .GE. alt_exp)  THEN
             data_zmax(i) = darr(i,l) 
          END IF
       END DO
    END DO
    !
    !  exponentially damp field above alt_exp to model top at l=1.
    !
    DO l = 1,lev2 
       DO i = il1,il2
          IF (alt(i,l) .GE. alt_exp)  THEN
             darr(i,l) = data_zmax(i) * EXP( (alt_exp-alt(i,l))/hscale )
          END IF
       END DO
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_exp

  SUBROUTINE vert_smooth (darr, work, coeff, nsmooth,         &
                          il1, il2, lev1, lev2, nlons, nlevs)
    !
    !  smooth a longitude by altitude array in the vertical over a
    !  specified number of levels using a three point smoother. 
    !
    !  note: input array darr is modified on output!
    !
    !  aug. 3/95 - c. mclandress
    !
    !  output arguement:
    !
    !     * darr    = smoothed array (on output).
    !
    !  input arguements:
    !
    !     * darr    = unsmoothed array of data (on input).
    !     * work    = work array of same dimension as darr.
    !     * coeff   = smoothing coefficient for a 1:coeff:1 stencil.
    !     *           (e.g., coeff = 2 will result in a smoother which
    !     *           weights the level l gridpoint by two and the two 
    !     *           adjecent levels (l+1 and l-1) by one).
    !     * nsmooth = number of times to smooth in vertical.
    !     *           (e.g., nsmooth=1 means smoothed only once, 
    !     *           nsmooth=2 means smoothing repeated twice, etc.)
    !     * il1     = first longitudinal index to use (il1 >= 1).
    !     * il2     = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1    = first altitude level to use (lev1 >=1). 
    !     * lev2    = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons   = number of longitudes.
    !     * nlevs   = number of vertical levels.
    !
    !  subroutine arguements.
    !
    INTEGER  :: nsmooth, il1, il2, lev1, lev2, nlons, nlevs
    REAL(dp) :: coeff
    REAL(dp) :: darr(nlons,nlevs), work(nlons,nlevs)
    !
    !  internal variables.
    !
    INTEGER  :: i, l, ns, lev1p, lev2m
    REAL(dp) :: sum_wts
    !-----------------------------------------------------------------------     
    !
    !  calculate sum of weights.
    !
    sum_wts = coeff + 2.0_dp
    !
    lev1p = lev1 + 1
    lev2m = lev2 - 1
    !
    !  smooth nsmooth times
    !
    DO ns = 1,nsmooth
       !
       !  copy darr into work array.
       !
       DO l = lev1,lev2
          DO i = il1,il2
             work(i,l) = darr(i,l)
          END DO
       END DO
       !
       !  smooth array work in vertical direction and put into darr.
       !
       DO l = lev1p,lev2m
          DO i = il1,il2
             darr(i,l) = (work(i,l+1)+coeff*work(i,l)+work(i,l-1) ) / sum_wts 
          END DO
       END DO
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vert_smooth

  SUBROUTINE  gw_fronts(krow, kproma, nazmth, rmswind, ani)
    !
    !
    !  may 22/2000 - m. charron
    !
    !  output arguements:
    !
    !     * rmswind    = root mean square gravity wave wind at lowest level (m/s).
    !     * ani        = anisotropy factor (sum over azimuths is one)
    !
    !  input arguments:
    !
    !     * klon     = number of longitudes 
    !     * nazmth     = azimuthal array dimension (nazmth >= naz).

    USE mo_gwspectrum,     ONLY: rmscon, rms_front, front_thres, naz
    !
    !  subroutine arguements.
    !
    INTEGER,  INTENT(IN)  :: krow, kproma, nazmth
    REAL(dp), INTENT(OUT) :: rmswind(kproma)
    REAL(dp), INTENT(OUT) :: ani(kproma,nazmth)
    !
    !  internal variables.
    !
    INTEGER  :: jl, jdir
    INTEGER  :: opp(nazmth)
    REAL(dp) :: angle(kproma), gen(kproma)


    IF ( naz == 8 ) THEN
       opp(1)=5
       opp(2)=6
       opp(3)=7
       opp(4)=8
       opp(5)=1
       opp(6)=2
       opp(7)=3
       opp(8)=4
    ELSE IF ( naz == 4 ) THEN
       opp(1)=3
       opp(2)=4
       opp(3)=1
       opp(4)=2
    END IF

    CALL calculate_gen ( krow, kproma, gen, angle )

    DO jl=1,kproma
       IF ( gen(jl) >= front_thres*2.7777778E-14_dp ) THEN
          rmswind(jl) = rms_front
          jdir=INT( angle(jl)/360.0_dp*REAL(naz,dp) ) + 1
          ani(jl,    :    ) = 0.0_dp
          ani(jl,    jdir ) = 0.5_dp
          ani(jl,opp(jdir)) = 0.5_dp
       ELSE
          rmswind(jl) = rmscon
          ani(jl,:)   = 1.0_dp/REAL(naz,dp)
       END IF
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE gw_fronts

  SUBROUTINE calculate_gen(krow, il2, gen, angle)
    !
    ! determine the value of the frontogenesis function
    !
    !  may 22/2000   - m. charron - first version
    !  march 29/2001 - m. charron - reorganized as a module
    !
    !  output arguement:
    !
    !     * gen    = frontogenesis function in (K/m)^2/hour
    !     * angle  = orientation (in degrees) of the gradiant of temp
    !
    !  input arguements:
    !
    !     * il2     = last longitudinal index to use (1 <= il2 <= nlons).
    
    USE mo_gaussgrid,  ONLY: gl_twomu, gl_sqcst
    USE mo_control,    ONLY: nvclev, vct, nlev
    USE mo_memory_g1a, ONLY: dalpslm1, dalpsmm1, vom1, dm1, alpsm1, tm1
    USE mo_memory_g2a, ONLY: dtlm1, dtmm1, um1, vm1, dudlm1, dvdlm1
    USE mo_gwspectrum, ONLY: emiss_lev, naz
    USE mo_geoloc,     ONLY: ilat

    !
    !  subroutine arguements.
    !
    INTEGER,                   INTENT(IN)  :: krow, il2
    REAL(dp), DIMENSION(il2),  INTENT(OUT) :: gen, angle
    !
    !  internal variables.
    !
    INTEGER                  :: jl, jglat, jrow, iplev
    REAL(dp), DIMENSION(il2) :: dalpsdl, dalpsdm, dudl, dvdl, ps  , vort, div
    REAL(dp), DIMENSION(il2) :: dtdl   , dtdm   , uu  , vv  , uup1, uum1, vvp1
    REAL(dp), DIMENSION(il2) :: vvm1   , tt     , ttp1, ttm1
    REAL(dp), DIMENSION(il2) :: mu, cstmu, zrcst
    REAL(dp)                 :: za, zb, zap1, zbp1
    REAL(dp)                 :: zaplus, zaminus, zbplus, zbminus
    REAL(dp)                 :: term_1, term_2, term_3l, term_3m, term_4l, term_4m
    REAL(dp)                 :: term_5, term_6, term_7 , term_8
    REAL(dp),PARAMETER       :: a=6371000.0_dp, pr=100000.0_dp, kappa=2.0_dp/7.0_dp, pi=3.141592654_dp

!------------------------------------------------------------------------------

    jrow=krow
    iplev=nlev-emiss_lev

    ps     (:)=EXP(alpsm1(:,jrow))
    za        =vct     (iplev)
    zap1      =vct     (iplev+1)
    zb        =vct     (iplev+nvclev)
    zbp1      =vct     (iplev+nvclev+1)
    zaplus    =zap1+za
    zaminus   =zap1-za
    zbplus    =zbp1+zb
    zbminus   =zbp1-zb
    DO jl=1,il2
      jglat     =ilat(jl,jrow)
      mu(jl)    =0.5_dp*gl_twomu(jglat)
      zrcst(jl) =1.0_dp/gl_sqcst(jglat)
    ENDDO
    cstmu  (:)=1.0_dp-mu(:)*mu(:)
    dtdl   (:)=a*cstmu(:)*dtlm1 (:,iplev,jrow)
    dtdm   (:)=a*dtmm1 (:,iplev,jrow)
    dalpsdl(:)=a*cstmu(:)*dalpslm1(:,jrow)
    dalpsdm(:)=a*dalpsmm1(:,jrow)

    uu  (:) = um1 (:,iplev  ,jrow)
    vv  (:) = vm1 (:,iplev  ,jrow)
    uup1(:) = um1 (:,iplev+1,jrow)
    uum1(:) = um1 (:,iplev-1,jrow)
    vvp1(:) = vm1 (:,iplev+1,jrow)
    vvm1(:) = vm1 (:,iplev-1,jrow)
    tt  (:) = tm1 (:,iplev  ,jrow)
    ttp1(:) = tm1 (:,iplev+1,jrow)
    ttm1(:) = tm1 (:,iplev-1,jrow)
    vort(:) = vom1(:,iplev  ,jrow)
    div (:) = dm1 (:,iplev  ,jrow)

    dudl(:) = dudlm1(:,iplev  ,jrow)*zrcst(:)
    dvdl(:) = dvdlm1(:,iplev  ,jrow)*zrcst(:)

    DO jl=1,il2
      term_1  = (2.0_dp*pr/(zaplus+ps(jl)*zbplus))**kappa
      term_2  = 1.0_dp/(zaminus+ps(jl)*zbminus)
      term_3l = ps(jl)*dalpsdl(jl)*zbplus*term_2
      term_3m = ps(jl)*dalpsdm(jl)*zbplus*term_2
      term_4l = ( dtdl(jl) - 0.25_dp*(ttp1(jl)-ttm1(jl)) * term_3l ) * term_1
      term_4m = ( dtdm(jl) - 0.25_dp*(ttp1(jl)-ttm1(jl)) * term_3m ) * term_1
      term_5  = 0.25_dp * (uup1(jl) - uum1(jl))
      term_6  = 0.25_dp * (vvp1(jl) - vvm1(jl))
      term_7  = 1.0_dp/(a*a*cstmu(jl))
      term_8  = a * term_7
      gen(jl) = -term_7*term_4l**2                                             &
           *(term_8*(dudl(jl)-term_5*term_3l)-mu(jl)*vv(jl)*term_8)            &
           -cstmu(jl)/(a*a) * term_4m**2                                       &
           *(1.0_dp/a*(a*div(jl)-dudl(jl)/cstmu(jl)-term_6*term_3m)            &
           +mu(jl)*vv(jl)*term_8)-term_7* term_4l*term_4m                      &
           *( 1.0_dp/a*(dvdl(jl)-term_6*term_3l)                               &
           +cstmu(jl)/a*(dvdl(jl)/cstmu(jl)-a*vort(jl)-term_5*term_3m)         &
           +2.0_dp*mu(jl)*uu(jl)/a)
      angle(jl) = ATAN2(cstmu(jl)*dtdm(jl),dtdl(jl))*180.0_dp/pi+180.0_dp/REAL(naz,dp)
      IF ( angle(jl) < 0.0_dp ) angle(jl) = angle(jl) + 360.0_dp
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE calculate_gen

END MODULE mo_midatm
