#if defined (NAG)
#define ARGCHECK 1
#endif
!OCL NOALIAS

SUBROUTINE physc(krow, kglat)
!
! Description:
!
! Controls the calls to the various physical subroutines.
!
!
!  *physc* is called from *gpc*.
!
!  Externals:
!
!  *geopot*    computes full level geopotentials.
!  *pres*      computes half level pressure.
!  *presf*     computes full level pressure.
!  *radiation* controls radiation computations.
!  *radheat*   computes radiation tendencies.
!  *vdiff*     computes vertical exchange by turbulence.
!  *ssodrag*   computes gravity wave drag.
!  *cucall*    controls mass-flux scheme.
!  *cloud*     computes large scale water phase changes and cloud cover.
!  *ml_ocean*  computes mixed layer ocean
!
!
!  Authors:
!
!  M. Jarraud, ECMWF, January 1982, original source
!
!     Modifications.
!     --------------
!  M.A. Giorgetta, MPI-Hamburg, May 2000, modified for ECHAM5
!  A.Tompkins,     MPI-Hamburg, June  2000, cloud cover scheme
!  U. Schlese,     MPI-Hamburg, July  2000, calling sequence changed
!  I. Kirchner,    MPI, July 2000, tendency diagnostics revision
!  L. Kornblueh,   MPI, August 2001, changes for different advection
!                       schemes, adapt polefilter call 
!  U. Schulzweida, MPI, May 2002, blocking (nproma)
!  I. Kirchner,    MPI, August 2002, nmi revision/extension
!  U. Schlese, M. Esch, MPI, September 2002, mixed layer ocean
!  L. Kornblueh,   MPI, August 2004, new tropopause calculation
!  S.J. Lorenz, M. Esch, MPI, Januar 2008, volcanic/dblrad
!
!
USE mo_kind,              ONLY: dp
USE mo_exception,         ONLY: message, finish, int2string
USE mo_memory_g1a,        ONLY: xlm1, xim1, tm1, qm1, alpsm1, xtm1
USE mo_memory_g2a,        ONLY: vm1, um1
USE mo_memory_g3a
USE mo_memory_g3b
!---wiso-code
USE mo_memory_wiso
!---wiso-code-end
USE mo_control,           ONLY: ltdiag, lcouple, lmidatm, lhd, lmlo,   &
                                nlev, nlevp1, ldiagamip, lso4, ltimer, &
                                lco2_nudge, lco2_flxcor, lco2_2perc,   &
                                l_volc                
USE mo_hyb,               ONLY: delb, nlevm1
USE mo_param_switches,    ONLY: lcond, lconv, lsurf, lgwdrag, lrad
USE mo_constants,         ONLY: cpd, vtmpc1, vtmpc2, g, tmelt
USE mo_radiation,         ONLY: solc, ldblrad
USE mo_cloud,             ONLY: ctaus, ctaul, ctauk
USE mo_hydrology,         ONLY: hydrology_collect
USE mo_scan_buffer,       ONLY: vo, vol, vom, qte, xlte, xite, tte, xtte, &
                                alnpr, alpha, alpste, vervel
USE mo_physc1,            ONLY: cdisse
USE mo_physc2,            ONLY: cqsncr, cwlmax
USE mo_tracer,            ONLY: ntrac, xtdriver1, xtdriver2, xtdiagn
!
USE mo_midatm,            ONLY: gwspectrum
USE mo_ssortns,           ONLY: ssodrag
!
USE mo_decomposition,     ONLY: ldc => local_decomposition
USE mo_diag_tendency,     ONLY: pdiga
USE mo_column,            ONLY: get_col_pol, lcotra
USE mo_geoloc,            ONLY: amu0_x, rdayl_x, sqcst_2d, &
                                philat_2d, budw_2d, twomu_2d
!
USE mo_time_control,      ONLY: lstart, lresume, delta_time, l_trigrad,&
                                lfirst_day, time_step_len, l_getocean
USE mo_advection
!---wiso-code
!USE mo_spitfire,          ONLY: pole_filter
USE mo_spitfire,          ONLY: pole_filter_wiso
!---wiso-code

USE mo_timer,             ONLY: timer_start, timer_stop, timer_radiation, timer_cloud
!
USE mo_nmi,               ONLY: nmi_phase, NMI_ACCU, NMI_USE_AVG,      &
                                dh_t, dh_m, dh_l, buf_t, buf_m, buf_l, &
                                lnmi_run, lnmi_cloud
USE mo_diag_amip2,        ONLY: collect_amip2_diag 
USE mo_diag_radiation,    ONLY: diag_rad1, diag_rad2
USE mo_tropopause,        ONLY: WMO_tropopause

!---wiso-code

USE mo_wiso,              ONLY: lwiso, nwiso

!---wiso-code-end

! JSBACH interface add-ons
USE mo_surface_memory,    ONLY: jrsfl, jrsfc, jssfl, jssfc, ztrfli,    &
                                zsofli, 			       &
!---wiso-code
                                jwisorsfl, jwisorsfc, jwisossfl, jwisossfc
!---wiso-code-end
USE mo_surface_memory,    ONLY: jsswnir, jsswdifnir, jsswvis, jsswdifvis, &
                                jsswniracc, jsswvisacc, jsswdifniracc, jsswdifvisacc

USE mo_co2,               ONLY: l_co2flxcorr, co2m1, co2atmos, co2flux_cpl, &
                                co2_burden, co2_burden_old, co2_flux_corr_acc, &
                                co2_burden_corr_acc1, co2_flux_corr_acc1, co2_flux_corr_acc1_old, &
                                co2_burden_corr_acc2, &
                                co2_flux, co2_flux_land, co2_flux_ocean, &
                                co2_emission, co2_flux_total, diag_co2, ico2idx, &
                                co2_nudge_mmr, co2_nudge_relax, co2_burden_nudge, co2_burden_nudge_acc
#ifdef __cpl_co2
USE mo_co2,               ONLY: co2_flux_atmosphere_ocean
#endif
USE mo_test
USE mo_jsbach, ONLY: test_stream

IMPLICIT NONE
!
INTEGER :: krow, kglat
! 
!  Local array bounds
INTEGER :: nglon, nproma, nbdim
!
!  Local scalars:
REAL(dp) :: zcst, zrcst, ztwodt, zprat,  zcdnc,  zn1,  zn2                 &
       , zepsec, zsigfac, zsigh, zsn_mm, zqdp, zco2eps, zco2tend, zco2p1, zco2_nudge_corr
INTEGER :: jlev, jl, kfdia, kidia, ktdia, jk, nexp
INTEGER :: jglat
LOGICAL :: loconv, locond

!---wiso-code

!  Local scalars - water isotopes

INTEGER :: jt

!---wiso-code-end
! 
!  Local arrays:
REAL(dp) ::  zdadc(ldc%nproma), zdpsdt(ldc%nproma), zgeo(ldc%nproma)          &
        ,ztvm1(ldc%nproma,nlev), ztslnew(ldc%nproma)                     &
! zbetaa: qt distribution minimum in beta
! zbetab: qt distribution maximum in beta
! zvdiffp:  dq/dt from vdiff scheme needed
! zhmixtau:  timescale of mixing for horizontal eddies
! zvmixtau:  timescale of mixing for vertical turbulence
        ,zbetaa(ldc%nproma,nlev)                                        &
        ,zbetab(ldc%nproma,nlev)                                        &
        ,zvdiffp(ldc%nproma,nlev)                                       &
        ,zhmixtau(ldc%nproma,nlev)                                      &
        ,zvmixtau(ldc%nproma,nlev)                                      &
        ,zi0(ldc%nproma)
REAL(dp) ::  zqtec(ldc%nproma,nlev), zbetass(ldc%nproma,nlev)
REAL(dp) ::  zqtold(ldc%nproma), zflux_corr_acc(ldc%nproma)
!
!---wiso-code

REAL(dp) ::  zwisoqtec(ldc%nproma,nlev,nwiso)
REAL(dp) ::  zwisoqtold(ldc%nproma,nwiso)

!---wiso-code-end

LOGICAL :: lonorth(ldc%nproma)

!  Surface fluxes over land/water/ice:
REAL(dp) :: zhfsl(ldc%nproma),  zhfsw(ldc%nproma),  zhfsi(ldc%nproma),  &
        zhflw(ldc%nproma),    &
        zevapw(ldc%nproma), zevapi(ldc%nproma), &
        ztrflw(ldc%nproma),                     &
        zsofll(ldc%nproma), zsoflw(ldc%nproma)

!---wiso-code
!  Surface fluxes over water/ice - water isotopes
! ATTENTION - In JSBACH the field zwisoevapl - surface evaporation (accumulated) over land is deleted!
!        zwisoevapw - surface evaporation (accumulated) over water
!        zwisoevapi - surface evaporation (accumulated) over ice
! wiso-code ECHAM5 stand alone:
!        REAL(dp) :: zwisoevapl(ldc%nproma,nwiso),zwisoevapw(ldc%nproma,nwiso),zwisoevapi(ldc%nproma,nwiso)
REAL(dp) :: zwisoevapw(ldc%nproma,nwiso),zwisoevapi(ldc%nproma,nwiso)
!---wiso-code-end

REAL(dp) ::  zfrl(ldc%nproma),  zfrw(ldc%nproma),   zfri(ldc%nproma)          &
        ,zcvsi(ldc%nproma), zcvsc(ldc%nproma),  zwlmx(ldc%nproma)         &
        ,zcvs(ldc%nproma),  zcvw(ldc%nproma)
REAL(dp) :: zteffl4(ldc%nproma),ztsnew(ldc%nproma), zradtemp_old(ldc%nproma)
REAL(dp) :: zraind(ldc%nproma), zsnowd(ldc%nproma)

!---wiso-code
REAL(dp) :: zwisoraind(ldc%nproma,nwiso), zwisosnowd(ldc%nproma,nwiso)
!---wiso-code-end

!  Local arrays for the HD-model and glacier calving model
REAL(dp) :: zros_hd(ldc%nproma), zdrain_hd(ldc%nproma)
REAL(dp) :: zalac(ldc%nproma)

INTEGER :: ilab(ldc%nproma,nlev), itype(ldc%nproma)
INTEGER :: invb(ldc%nproma)
!
INTEGER :: itrpwmo(ldc%nproma), itrpwmop1(ldc%nproma)
!
!---wiso-code

!  Local arrays for the HD-model and glacier calving model - water isotopes
REAL(dp) :: zwisoros_hd(ldc%nproma,nwiso), zwisodrain_hd(ldc%nproma,nwiso)
REAL(dp) :: zwisoalac(ldc%nproma,nwiso)

!---wiso-code-end

!
!    Arrays internal to physics
REAL(dp) :: qhfla(ldc%nproma), evapot(ldc%nproma), zprecip(ldc%nproma)
!
!---wiso-code

REAL(dp) :: wisoqhfla(ldc%nproma,nwiso), wisoevapot(ldc%nproma,nwiso), zwisoprecip(ldc%nproma,nwiso)

!---wiso-code-end

REAL(dp) :: zdtime

! Local arrays
LOGICAL  :: loland(ldc%nproma), loglac(ldc%nproma)
REAL(dp) :: geom1(ldc%nproma,nlev)
REAL(dp) :: aphm1(ldc%nproma,nlevp1), apm1(ldc%nproma,nlev) 
REAL(dp) :: aphp1(ldc%nproma,nlevp1), app1(ldc%nproma,nlev)
REAL(dp) :: rsfc(ldc%nproma), ssfc(ldc%nproma)
REAL(dp) :: rsfl(ldc%nproma), ssfl(ldc%nproma)

!---wiso-code

! Local arrays - water isotopes
REAL(dp) :: wisorsfc(ldc%nproma,nwiso), wisossfc(ldc%nproma,nwiso)
REAL(dp) :: wisorsfl(ldc%nproma,nwiso), wisossfl(ldc%nproma,nwiso)

!---wiso-code-end
! 
!  External subroutines
EXTERNAL geopot, pres, presf, radiation, vdiff, cloud,                 &
         cucall, radheat, collect
! 
!  Intrinsic functions
INTRINSIC EXP
!
!  Local array bounds
  jglat = kglat           ! global continuous latitude index
  nglon = ldc% nglon      ! local number of longitudes

  nbdim = ldc% nproma

  IF ( krow == ldc% ngpblks ) THEN
    nproma = ldc% npromz
  ELSE
    nproma = ldc% nproma
  END IF
!
!*    COMPUTATIONAL CONSTANTS.
!     ------------- ----------
!
  zdtime = delta_time
  ztwodt=time_step_len
  zepsec=1.E-12_dp
  zsigfac=0.15_dp
!
!     ----------------------------------------------------------------
!
!*        2.    ALLOCATE STORAGE 
!               -------- ------- 
!
200 CONTINUE
!
!     ------------------------------------------------------------
!
!*        3.    COMPUTE SOME FIELDS NEEDED BY THE PHYSICAL ROUTINES.
!               ------- ---- ------ ------ -- --- -------- ---------
!

300 CONTINUE
!
!*        3.0   Set start conditions for JSBACH interface variables
!
! INITIALIZE SURFACE PARAMETERS WHICH COME FROM MO_SURFACE IN THE 
!  NEXT TIME STEP
  IF (lstart) THEN
     IF (krow == 1) CALL message('physc','Set rain and snow to ZERO for JSBACH (first time step)')
     jrsfl(:,krow) = 0.0_dp  ! convectiv, large scale
     jrsfc(:,krow) = 0.0_dp  ! snow and rain fall
     jssfl(:,krow) = 0.0_dp  ! with units in
     jssfc(:,krow) = 0.0_dp  ! kg/m2/s, set for initial values

 !---wiso-code
    IF (lwiso) THEN

    DO jt=1,nwiso
        jwisorsfl(:,jt,krow) = 0.0_dp  ! large scale rain fall in kg/m2/s, set for initial values - water isotopes
        jwisorsfc(:,jt,krow) = 0.0_dp  ! convectiv rain fall in kg/m2/s, set for initial values - water isotopes
        jwisossfl(:,jt,krow) = 0.0_dp  ! large scale snow fall in kg/m2/s, set for initial values - water isotopes
        jwisossfc(:,jt,krow) = 0.0_dp  ! convectiv snow fall in kg/m2/s, set for initial values - water isotopes
     END DO
     
     END IF
 !---wiso-code-end
  END IF

  zraind = 0._dp
  zsnowd = 0._dp
  zraind(1:nproma) = jrsfl(1:nproma,krow) + jrsfc(1:nproma,krow)
  zsnowd(1:nproma) = jssfl(1:nproma,krow) + jssfc(1:nproma,krow)

!---wiso-code
 IF (lwiso) THEN

 zwisoraind = 0._dp
 zwisosnowd = 0._dp
 DO jt=1,nwiso
  zwisoraind(1:nproma,jt) = jwisorsfl(1:nproma,jt,krow) + jwisorsfc(1:nproma,jt,krow)
  zwisosnowd(1:nproma,jt) = jwisossfl(1:nproma,jt,krow) + jwisossfc(1:nproma,jt,krow)
 END DO
 
 END IF
!---wiso-code-end

!
!*        3.1   COMPUTE VIRTUAL TEMPERATURE AT T-DT AND SET *ZGEO* TO 0.
!
310 CONTINUE
!
  ztvm1(1:nproma,:) = tm1(1:nproma,:,krow)*(1._dp+vtmpc1*qm1(1:nproma,:,krow) &
                      -(xlm1(1:nproma,:,krow)+xim1(1:nproma,:,krow)))
!
  zgeo(1:nproma)=0._dp
!
!*        3.2   COMPUTE (PHI-PHIS) AT T-DT USING LN(P) AT T.
!
320 CONTINUE
!
  CALL geopot(geom1,ztvm1,alnpr(:,:,krow),alpha(:,:,krow),zgeo,nbdim,nproma)
!
!
!*        3.3    COMPUTE PRESSURE AT FULL AND HALF LEVELS AT T-DT.
!
330 CONTINUE
!
  aphm1(1:nproma,nlevp1)=EXP(alpsm1(1:nproma,krow))
!
  CALL pres(aphm1,nbdim,aphm1(1,nlevp1),nproma)
!
  CALL presf(apm1,nbdim,aphm1,nproma)
!
!*        3.4   COMPUTE REAL WINDS AND WIND TENDENCIES.
!
340 CONTINUE
!
  DO jlev = 1, nlev
!DIR$ CONCURRENT
    DO jl = 1, nproma
      zrcst=1._dp/sqcst_2d(jl,krow)
      um1(jl,jlev,krow) = zrcst*um1(jl,jlev,krow)
      vm1(jl,jlev,krow) = zrcst*vm1(jl,jlev,krow)
      vol(jl,jlev,krow) = zrcst*vol(jl,jlev,krow)
      vom(jl,jlev,krow) = zrcst*vom(jl,jlev,krow)
    END DO
  END DO
!
! Horizontal wind shear for horizontal mixing of variance
!
  DO jlev=1,nlev
     DO jl=1,nproma
        zhmixtau(jl,jlev)=ctauk*ABS(vo(jl,jlev,krow)) 
        zhmixtau(jl,jlev)=MIN(ctaus,MAX(ctaul,zhmixtau(jl,jlev)))
        zvmixtau(jl,jlev)=0._dp
     END DO
  END DO
!
! ---------------------------!!!!!-----------------------------------
!   no moisture computations at the first day to prevent instability
!   at the beginning caused by unbalanced initial fields  !!!!
!
!     IF(lfirst_day) THEN
!       locond=.false.
!       loconv=.false.
!     ELSE
       locond=lcond
       loconv=lconv
!     END IF

     IF (lnmi_run) THEN
       ! control cloud parameterisation in nmi initialization mode
       locond=lnmi_cloud
       loconv=lnmi_cloud
       SELECT CASE(nmi_phase)
       CASE(NMI_ACCU)     ! store tendencies
!DIR$ CONCURRENT
         buf_t(:,:) = tte(:,:,krow)
!DIR$ CONCURRENT
         buf_m(:,:) = vom(:,:,krow)
!DIR$ CONCURRENT
         buf_l(:,:) = vol(:,:,krow)
       CASE(NMI_USE_AVG)  ! store tendencies
!DIR$ CONCURRENT
         buf_t(:,:) = tte(:,:,krow)
!DIR$ CONCURRENT
         buf_m(:,:) = vom(:,:,krow)
!DIR$ CONCURRENT
         buf_l(:,:) = vol(:,:,krow)
       CASE default
       END SELECT
     END IF
!
! ------------------------------------------------------------------
!
!*        3.5   ESTIMATE ADIABATIC CONVERSION OF POTENTIAL ENERGY.
!
350 CONTINUE
!
  zdadc(1:nproma)=0._dp
!
  DO 351 jl=1,nproma
     zdpsdt(jl)=aphm1(jl,nlevp1)*alpste(jl,krow)
     zdadc(jl)=zdadc(jl)+geospm(jl,krow)*zdpsdt(jl)
351 END DO
!
  DO 353 jlev=1,nlev
     DO 352 jl=1,nproma
        zdadc(jl)=zdadc(jl)+(1._dp+vtmpc2*qm1(jl,jlev,krow))*cpd*   &
                (tte(jl,jlev,krow)*(aphm1(jl,jlev+1)-aphm1(jl,jlev)) &
                        +tm1(jl,jlev,krow)*delb(jlev)*zdpsdt(jl))
352  END DO
353 END DO
!
!
!*        3.6   COMPUTE LOGICAL MASK FOR LAND AND GLACIER.
!
360 CONTINUE
!
  DO 365 jl=1,nproma
     loland(jl)=slm(jl,krow).GT.0._dp
     loglac(jl)=loland(jl).AND.glac(jl,krow).GT.0._dp
     lonorth(jl)=philat_2d(jl,krow).GT.0 ! true in northern hemisphere
365 END DO
!
!       3.7 Weighting factors for fractional surface coverage
!           Accumulate ice portion for diagnostics
!
!DIR$ CONCURRENT
   DO jl=1,nproma
      zfrl(jl)=slm(jl,krow)
      zfrw(jl)=(1._dp-zfrl(jl))*(1._dp-seaice(jl,krow))
      zfri(jl)=1._dp-zfrl(jl)-zfrw(jl) 
!!$      friac(jl,krow)=friac(jl,krow)+zdtime*zfri(jl)
!! friac now directly imported from mo_surface into g3b stream
   END DO
  IF(lcouple) THEN
    DO jl=1,nproma
       IF(slf(jl,krow).GT.1.0_dp-zepsec) THEN
         tsi(jl,krow)=tmelt
         tsw(jl,krow)=tmelt
       END IF
    END DO
  ENDIF
!
!      3.8  Skin reservoir, wet skin fraction and snow cover
!           (bare land, canopy, lake ice)
!
   DO jl=1,nproma
     IF (.NOT.loglac(jl)) THEN
        zwlmx(jl)=cwlmax*(1._dp+vlt(jl,krow))
        zcvw(jl)=MIN(wl(jl,krow)/zwlmx(jl),1.0_dp)
        zsn_mm=1000._dp*sn(jl,krow)
        zsigh=SQRT(zsn_mm/(zsn_mm+zepsec+zsigfac*orostd(jl,krow)))
        zcvs(jl)=cqsncr*TANH(zsn_mm/10._dp)*zsigh
        zcvsc(jl)=MIN(1._dp,snc(jl,krow)/(zwlmx(jl)-cwlmax+EPSILON(1._dp)))
        IF (zcvs(jl).LT.EPSILON(1._dp) .AND. zcvsc(jl).GE.EPSILON(1._dp)) THEN
           zcvs(jl)=zcvsc(jl)
        END IF
     ELSE
        zwlmx(jl)=0._dp
        zcvw(jl)=0._dp
        zcvs(jl)=1._dp
        zcvsc(jl)=0._dp
     END IF

   END DO
!
!---wiso-code
 IF (lwiso) THEN

! set snow depth on glaciers to zero for all non-glacier points

 DO jl=1,nproma
   IF (.NOT.loglac(jl)) THEN
     snglac(jl,krow) = 0._dp
     wisosnglac(jl,:,krow) = 0._dp
   END IF
 END DO
 
 END IF
!---wiso-code-end
!
370 CONTINUE
!
!
!*         3.8    SET LOOP  VALUES FOR PHYSICAL PARAMETERIZATIONS
!
  kidia=1
  kfdia=nproma
  ktdia=1
!
!         3.9   DETERMINE TROPOPAUSE HEIGHT AND MASS BUDGETS
!
!
  CALL WMO_tropopause (nproma, nbdim, nlev,                 &
                       tm1(:,:,krow),  apm1, tropo(:,krow), &
                       itrpwmo, itrpwmop1)
!
!
!*    3.12 INITIALISATION OF CLOUD DROPLET NUMBER CONCENTRATION 
!          (1/M**3) USED IN RADLSW AND CLOUD 
!
  IF (lstart) THEN
     DO 4103 jk=ktdia,nlev        
!DIR$ CONCURRENT
        DO 4102 jl=kidia,kfdia       
           nexp=2
           zprat=(MIN(8._dp,80000._dp/apm1(jl,jk)))**nexp
           IF (loland(jl).AND.(.NOT.loglac(jl))) THEN
              zn1= 50._dp
              zn2=220._dp
           ELSE 
              zn1= 50._dp
              zn2= 80._dp
           ENDIF
           IF (apm1(jl,jk).LT.80000._dp) THEN
              zcdnc=1.e6_dp*(zn1+(zn2-zn1)*(EXP(1._dp-zprat)))
           ELSE
              zcdnc=zn2*1.e6_dp
           ENDIF
           acdnc(jl,jk,krow)=zcdnc
           acdncm(jl,jk,krow)=zcdnc
4102    END DO
4103 END DO
  ENDIF
!
     DO 4104 jl=kidia,kfdia       
           itype(jl)=NINT(rtype(jl,krow))
4104 END DO
!
     DO 4106 jk=ktdia,nlev
        DO 4105 jl=kidia,kfdia
           zqtec(jl,jk)=0.0_dp
4105    END DO
4106 END DO
!
!---wiso-code
  IF (lwiso) THEN

       DO jt=1,nwiso
        DO jk=ktdia,nlev
         DO jl=kidia,kfdia
             zwisoqtec(jl,jk,jt)=0.0_dp
         END DO
        END DO
       END DO

  END IF
!---wiso-code-end

!
!*        3.13   DIAGNOSE CURRENT CLOUD COVER
!
  IF(locond) THEN

     IF (ltimer) CALL timer_start(timer_cloud)

#ifdef ARGCHECK
     CALL cover( nproma, nbdim, ktdia, nlev, nlevp1                   &
               , itype,             zfrw                              &
               , invb,              rintop(:,krow)                    &
               , aphm1,             apm1                              &
               , qm1(:,:,krow),     tm1(:,:,krow)                     &
               , xlm1(:,:,krow),    xim1(:,:,krow)                    &
               , vervel(:,:,krow)                                     &
               , xvar(:,:,krow),    xskew(:,:,krow)                   &
               , aclc(:,:,krow)                                       &
               , zbetaa,            zbetab                            &
               , zbetass                                              &
                )
#else
     CALL cover( nproma, nbdim, ktdia, nlev, nlevp1                   &
               , itype,             zfrw                              &
               , invb,              rintop(1,krow)                    &
               , aphm1,             apm1                              &
               , qm1(1,1,krow),     tm1(1,1,krow)                     &
               , xlm1(1,1,krow),    xim1(1,1,krow)                    &
               , vervel(1,1,krow)                                     &
               , xvar(1,1,krow),    xskew(1,1,krow)                   &
               , aclc(1,1,krow)                                       &
               , zbetaa,            zbetab                            &
               , zbetass                                              &
                )
#endif
     IF (ltimer) CALL timer_stop(timer_cloud)

  ENDIF
!
!*        4.    RADIATION PARAMETERISATION.
!               --------- -----------------
!
400 CONTINUE
!
!
410 CONTINUE
!

  IF (lrad) THEN

  IF (l_trigrad) THEN

    IF (ltimer) CALL timer_start(timer_radiation)

#ifdef ARGCHECK
    CALL radiation(nproma,nbdim,nlev,nlevp1,krow,jglat,         &
                   twomu_2d(:,krow),budw_2d(:,krow),            &
                   aphm1,apm1,tm1(:,:,krow),                    &
                   co2m1(:,:,krow),                             &
                   qm1(:,:,krow),                               &
                   acdnc(:,:,krow),aclc(:,:,krow),              &
                   xlm1(:,:,krow),xim1(:,:,krow),               &
                   loland,loglac,                               &
                   albedo(:,krow),                              &
                   albedo_vis(:,krow),albedo_nir(:,krow),       &
                   radtemp(:,krow),                             &
                   ao3(:,:,krow),                               &
                   aclcv(:,krow),                               &
                   emter(:,:,krow),trsol(:,:,krow),             &
                   emtef(:,:,krow),trsof(:,:,krow),             &
                   emtef0(:,:,krow),trsof0(:,:,krow),           &
                   so4nat(:,:,krow),so4all(:,:,krow),           &
                   diag_rad1,diag_rad2,                         &
                   sswnir(:,krow),sswdifnir(:,krow),            &
                   sswvis(:,krow),sswdifvis(:,krow))
#else
    CALL radiation(nproma,nbdim,nlev,nlevp1,krow,jglat,         &
                   twomu_2d(1,krow),budw_2d(1,krow),            &
                   aphm1,apm1,tm1(1,1,krow),                    &
                   co2m1(1,1,krow),                             &
                   qm1(1,1,krow),                               &
                   acdnc(1,1,krow),aclc(1,1,krow),              &
                   xlm1(1,1,krow),xim1(1,1,krow),               &
                   loland,loglac,                               &
                   albedo(1,krow),                              &
                   albedo_vis(1,krow),albedo_nir(1,krow),       &
                   radtemp(1,krow),                             &
                   ao3(1,1,krow),                               &
                   aclcv(1,krow),                               &
                   emter(1,1,krow),trsol(1,1,krow),             &
                   emtef(1,1,krow),trsof(1,1,krow),             &
                   emtef0(1,1,krow),trsof0(1,1,krow),           &
                   so4nat(1,1,krow),so4all(1,1,krow),           &
                   diag_rad1,diag_rad2,                         &
                   sswnir(1,krow),sswdifnir(1,krow),            &
                   sswvis(1,krow),sswdifvis(1,krow))
#endif

!
!   Second call of radiation computation
!     for diagnose of instantaneous fluxes (switch dblrad)
!     double radiation currently used for volcanic forcing only
!
    IF (ldblrad .AND. l_volc) THEN
!
!   Second call of radiation computation without volcanic forcing
       l_volc = .FALSE.

#ifdef ARGCHECK
    CALL radiation(nproma,nbdim,nlev,nlevp1,krow,jglat,         &
                   twomu_2d(:,krow),budw_2d(:,krow),            &
                   aphm1,apm1,tm1(:,:,krow),                    &
                   co2m1(:,:,krow),                             &
                   qm1(:,:,krow),                               &
                   acdnc(:,:,krow),aclc(:,:,krow),              &
                   xlm1(:,:,krow),xim1(:,:,krow),               &
                   loland,loglac,                               &
                   albedo(:,krow),                              &
                   albedo_vis(:,krow),albedo_nir(:,krow),       &
                   radtemp(:,krow),                             &
                   ao3(:,:,krow),                               &
                   aclcv(:,krow),                               &
                   emter1(:,:,krow),trsol1(:,:,krow),           &
                   emtef1(:,:,krow),trsof1(:,:,krow),           &
                   emtef01(:,:,krow),trsof01(:,:,krow),         &
                   so4nat(:,:,krow),so4all(:,:,krow),           &
                   diag_rad1,diag_rad2,                         &
                   sswnir(:,krow),sswdifnir(:,krow),            &
                   sswvis(:,krow),sswdifvis(:,krow))
#else
    CALL radiation(nproma,nbdim,nlev,nlevp1,krow,jglat,         &
                   twomu_2d(1,krow),budw_2d(1,krow),            &
                   aphm1,apm1,tm1(1,1,krow),                    &
                   co2m1(1,1,krow),                             &
                   qm1(1,1,krow),                               &
                   acdnc(1,1,krow),aclc(1,1,krow),              &
                   xlm1(1,1,krow),xim1(1,1,krow),               &
                   loland,loglac,                               &
                   albedo(1,krow),                              &
                   albedo_vis(1,krow),albedo_nir(1,krow),       &
                   radtemp(1,krow),                             &
                   ao3(1,1,krow),                               &
                   aclcv(1,krow),                               &
                   emter1(1,1,krow),trsol1(1,1,krow),           &
                   emtef1(1,1,krow),trsof1(1,1,krow),           &
                   emtef01(1,1,krow),trsof01(1,1,krow),         &
                   so4nat(1,1,krow),so4all(1,1,krow),           &
                   diag_rad1,diag_rad2,                         &
                   sswnir(1,krow),sswdifnir(1,krow),            &
                   sswvis(1,krow),sswdifvis(1,krow))

#endif
!
!   Switch back to use of volcanic forcing
       l_volc = .TRUE.

    ENDIF ! ldblrad
!
    IF (ltimer) CALL timer_stop(timer_radiation)

  END IF

  ELSE
     ! lrad=.FALSE.
     ! --> no radiative effect of the atmosphere or the surface
     emter(:,:,krow)=0._dp
     trsol(:,:,krow)=1._dp
     emtef(:,:,krow)=0._dp
     trsof(:,:,krow)=1._dp
     emtef0(:,:,krow)=0._dp
     trsof0(:,:,krow)=1._dp
     ! --> no radiative heating or cooling of 
     !     (1) the surface, as computed in vdiff, or
     !     (2) the atmosphere, as computed in radheat
  END IF

!
! ----------------------------------------------------------------------
!
!       Update solar incidence *zi0* 
!     
!
  IF (lrad) THEN
     DO  jl=1,nproma
        zi0(jl)=cdisse*solc*amu0_x(jl,krow)*rdayl_x(jl,krow)
     END DO
  ELSE
     DO  jl=1,nproma
        zi0(jl)=0._dp
     END DO
  END IF

!
!       Compute diffuse, direct, NIR and visible fluxes for JSBACH
!
     DO  jl=1,nproma
        jsswnir(jl,krow) = zi0(jl) * sswnir(jl,krow)
        jsswdifnir(jl,krow) = sswdifnir(jl,krow)
        jsswvis(jl,krow) = zi0(jl) * sswvis(jl,krow)
        jsswdifvis(jl,krow) = sswdifvis(jl,krow)
        jsswniracc(jl,krow) = jsswniracc(jl,krow) + jsswnir(jl,krow) * zdtime
        jsswdifniracc(jl,krow) = jsswdifniracc(jl,krow) + jsswdifnir(jl,krow)  * jsswnir(jl,krow) * zdtime
        jsswvisacc(jl,krow) = jsswvisacc(jl,krow) + jsswvis(jl,krow) * zdtime
        jsswdifvisacc(jl,krow) = jsswdifvisacc(jl,krow) + jsswdifvis(jl,krow) * jsswvis(jl,krow) * zdtime
     END DO

!
!     ------------------------------------------------------------
!
!*              VERTICAL EXCHANGE OF U,V,T,Q BY TURBULENCE.
!               -------- -------- -- - - - - -- -----------
!
!
!      COMPUTE PRESSURE AT FULL AND HALF LEVELS AT T+DT.
!
  DO 522 jl=1,nproma
     aphp1(jl,nlevp1)=EXP(alpsm1(jl,krow)+ztwodt*alpste(jl,krow))
522 END DO
!
  CALL pres(aphp1,nbdim,aphp1(1,nlevp1),nproma)
!
  CALL presf(app1,nbdim,aphp1,nproma)
!
  IF (ltdiag) THEN
! prepare next fields for VDIFF
    pdiga(1:nglon,:, 3,krow) = pdiga(1:nglon,:, 3,krow) - vom(1:nproma,:,krow)
    pdiga(1:nglon,:, 8,krow) = pdiga(1:nglon,:, 8,krow) - vol(1:nproma,:,krow)
    pdiga(1:nglon,:,16,krow) = pdiga(1:nglon,:,16,krow) - tte(1:nproma,:,krow)
  ENDIF
!
  ! Compute CO2 burden and correction term
         
  co2_burden_old(1:nproma,krow) = co2_burden(1:nproma,krow)
  co2_burden(1:nproma,krow) = 0._dp
  DO jk=ktdia,nlev
     co2_burden(1:nproma,krow) = co2_burden(1:nproma,krow) &
             + co2m1(1:nproma,jk,krow) * (aphm1(1:nproma,jk+1)-aphm1(1:nproma,jk))/g
  END DO

  IF (ico2idx > 0 .AND. lco2_flxcor) THEN
     IF (ALL(co2_burden_old(1:nproma,krow) < 0._dp)) THEN
        IF (krow == 1) THEN
           CALL message('physc','co2_burden < 0')
           CALL message('physc','   assuming that it was not present in restart file,')
           CALL message('physc','   setting correction in CO2 burden to zero')
           CALL message('physc','   (this should only happen on the first time step of the first')
           CALL message('physc','    restart after introduction of flux correction, or for lstart=true)')
        END IF
        IF (ANY(co2_flux_total(1:nproma,krow) /= 0._dp)) THEN
           CALL finish('physc','co2_flux_total should also be zero in this case!')
        END IF
        IF (ANY(co2_flux_corr_acc1_old /= 0._dp)) THEN
           CALL finish('physc','co2_flux_corr_acc1_old should also be zero in this case!')
        END IF        
        co2_burden_old(1:nproma,krow) = co2_burden(1:nproma,krow)
     END IF

     IF (l_co2flxcorr .OR. lstart) THEN
        co2_burden_corr_acc1(1:nproma,krow) = 0._dp ! set to zero after coupling
        IF (krow == 1) CALL message('physc','Resetting acc. CO2 burden correction')
     ELSE
        co2_flux_corr_acc1_old(1:nproma,krow) = co2_flux_corr_acc1(1:nproma,krow)
     END IF

     ! Compute flux correction field from (old burden + fluxes + global mean flux corr.) - (new burden).
     ! This is the instant. flux correction multiplied by zdtime, i.e. it is an accumulated value
     ! Uses old co2_flux_total which is only updated in vdiff
     zflux_corr_acc = 0.0_dp
     zflux_corr_acc(1:nproma) = &
          (co2_burden_old(1:nproma,krow) + (co2_flux_total(1:nproma,krow)+co2_flux_corr_acc1_old(1:nproma,krow))*zdtime) &
         - co2_burden(1:nproma,krow)
!     IF (lco2_nudge) THEN
!        NOTE: This doesn't work as intended!!!
!        ! Compute flux correction from burden in the absence of nudging; otherwise, the flux
!        ! correction would falsely "correct", i.e. counteract, the nudging term, which is an
!        ! additional, artificial source/sink of CO2 into the system
!        ! (co2_burden in the previous equation contains the CO2 added/subtracted by nudging
!        ! in the last time step, while co2_burden_old doesn't contain this source/sink.
!        zflux_corr_acc(1:nproma) = zflux_corr_acc(1:nproma) + co2_burden_nudge(1:nproma,krow)
!     END IF

     ! Accumulate flux correction field for output interval
     ! (co2_flux_corr_acc has laccu=.TRUE., i.e. is divided by interval length before output)
     co2_flux_corr_acc(1:nproma,krow) = co2_flux_corr_acc(1:nproma,krow) + zflux_corr_acc(1:nproma)

     ! Accumulate CO2 burden correction for daily correction update
     ! (co2_burden_corr_acc1 has laccu=.FALSE. and is reset above to zero at correction time step,
     !  but not divided by interval length)
     co2_burden_corr_acc1(1:nproma,krow) = co2_burden_corr_acc1(1:nproma,krow) + zflux_corr_acc(1:nproma)
  END IF

#ifdef __cpl_co2
  CALL co2_flux_atmosphere_ocean(krow, nproma, nlev)
#endif

#ifdef ARGCHECK
  IF (lwiso) THEN 
    CALL vdiff (nproma, nbdim, ktdia, nlev, nlevm1, nlevp1, ntrac     &
            , krow                                                    &
!
            , xtm1(:,:,:,krow)                                        &   ! contains CO2 concentration
!
            , qm1(:,:,krow),     tm1(:,:,krow),     um1(:,:,krow)     &
            , vm1(:,:,krow),     xlm1(:,:,krow),    xim1(:,:,krow)    &
            , xvar(:,:,krow)                                          &
!
            , ahfl(:,krow),      ahfs(:,krow),      az0(:,krow)       &
            , dew2(:,krow),      evap(:,krow),      forest(:,krow)    &
            , temp2(:,krow),     t2max(:,krow)                        &
            , t2min(:,krow),     wind10w(:,krow),   vdis(:,krow)      &
            , u10(:,krow),       v10(:,krow),       ustr(:,krow)      &
            , vstr(:,krow),      wimax(:,krow),     wind10(:,krow)    &
            , wsmx(:,krow),      vlt(:,krow)                          &
            , vgrat(:,krow)                                           &
            , tsl(:,krow),       tsw(:,krow),       tsi(:,krow)       &
            , ocu(:,krow),       ocv(:,krow)                          &
            , az0l(:,krow),      az0w(:,krow),      az0i(:,krow)      &
            , zhfsl,             zhfsw,             zhfsi             &
            , zhflw                                                   &
            , zevapw,            zevapi                               &
            , ahfslac(:,krow),   ahfswac(:,krow),   ahfsiac(:,krow)   &
            , ahfllac(:,krow),   ahflwac(:,krow),   ahfliac(:,krow)   &
            , evaplac(:,krow),   evapwac(:,krow),   evapiac(:,krow)   &
            , ustrl(:,krow),     ustrw(:,krow),     ustri(:,krow)     &
            , vstrl(:,krow),     vstrw(:,krow),     vstri(:,krow)     &
            , sn(:,krow),        snc(:,krow),       tslm1(:,krow)     &
            , ws(:,krow),        albedo(:,krow),    albedo_vis(:,krow) &
            , albedo_nir(:,krow),alsol(:,krow)                        &
!
            , tke(:,:,krow),     tkem1(:,:,krow),   tkem(:,:,krow)    &
            , aclc(:,:,krow),    emter(:,:,krow)                      &
!
            , aphm1,             apm1,              geom1             &
            , ztvm1,             zvdiffp,           zvmixtau          &
!
            , zcvs,              zcvw                                 &
            , qhfla,             evapot                               &
            , ztslnew,           zwlmx                                &
            , zfrl,              zfrw,              zfri              &
            , loland,            loglac                               &
!
            , xtte(:,:,:,krow)                                        &
!
            , vol(:,:,krow),     vom(:,:,krow),     qte(:,:,krow)     &
            , tte(:,:,krow),     xlte(:,:,krow),    xite(:,:,krow)    &
            , siced(:,krow)                                           &
! Rain and snow over last time step
            , zraind(:),            zsnowd(:)                          &
            , jsswnir(:,krow), jsswdifnir(:,krow)                      &
            , jsswvis(:,krow), jsswdifvis(:,krow)                      &
            , zi0,                    trsol(:,:,krow)                  &
            , ztrflw,           ztrfli(:,krow)      &
            , zsofll,            zsoflw,           zsofli(:,krow)      &
            , trfllac(:,krow),   trflwac(:,krow),  trfliac(:,krow)     &
            , sofllac(:,krow),   soflwac(:,krow),  sofliac(:,krow)     &
            , alsoi(:,krow),    alsow(:,krow),     radtemp(:,krow)     &
            , sni(:,krow)                                              &
            , ahfice(:,krow)                                           &
            , fluxres(:,krow)                                          &
            , qres(:,krow)                                             &
            , ahfcon(:,krow)                                           &
            , ahfres(:,krow)                                           &
            , zteffl4 , ztsnew , tsurf(:,krow)                         &
            , seaice(:,krow)                                           &
            , zradtemp_old, runoff(1:nbdim,krow), drain(1:nbdim,krow)  &
            , zalac, zros_hd, zdrain_hd                                &
            , co2_flux_ocean(1:nbdim,krow), co2_flux_land(1:nbdim,krow)&
            , co2_flux(1:nbdim,krow), co2_flux_corr_acc1(1:nbdim,krow) &
!---wiso-code
            , lwiso, nwiso                                            &
            , wisoqm1(:,:,:,krow), wisoxlm1(:,:,:,krow), wisoxim1(:,:,:,krow)     &
            , wisoevap(:,:,krow)                                                  &
            , zwisoevapw, zwisoevapi                                              &
            , wisoevaplac(:,:,krow), wisoevapwac(:,:,krow), wisoevapiac(:,:,krow) &
            , wisosn(:,:,krow), wisows(:,:,krow)                                  &
            , wisosw_d(:,:,krow)                                                  &       
            , wisoqhfla, wisoevapot                                               &
            , wisoqte(:,:,:,krow), wisoxlte(:,:,:,krow), wisoxite(:,:,:,krow)     &
! Rain and snow over last time step
            , zwisoraind(:,:), zwisosnowd(:,:)                                    &
            , wisorunoff(:,:,krow), wisodrain(:,:,krow)                           &
            , zwisoalac, zwisoros_hd, zwisodrain_hd)
!---wiso-code-end
  ELSE
    CALL vdiff (nproma, nbdim, ktdia, nlev, nlevm1, nlevp1, ntrac     &
            , krow                                                    &
!
            , xtm1(:,:,:,krow)                                        &   ! contains CO2 concentration
!
            , qm1(:,:,krow),     tm1(:,:,krow),     um1(:,:,krow)     &
            , vm1(:,:,krow),     xlm1(:,:,krow),    xim1(:,:,krow)    &
            , xvar(:,:,krow)                                          &
!
            , ahfl(:,krow),      ahfs(:,krow),      az0(:,krow)       &
            , dew2(:,krow),      evap(:,krow),      forest(:,krow)    &
            , temp2(:,krow),     t2max(:,krow)                        &
            , t2min(:,krow),     wind10w(:,krow),   vdis(:,krow)      &
            , u10(:,krow),       v10(:,krow),       ustr(:,krow)      &
            , vstr(:,krow),      wimax(:,krow),     wind10(:,krow)    &
            , wsmx(:,krow),      vlt(:,krow)                          &
            , vgrat(:,krow)                                           &
            , tsl(:,krow),       tsw(:,krow),       tsi(:,krow)       &
            , ocu(:,krow),       ocv(:,krow)                          &
            , az0l(:,krow),      az0w(:,krow),      az0i(:,krow)      &
            , zhfsl,             zhfsw,             zhfsi             &
            , zhflw                                                   &
            , zevapw,            zevapi                               &
            , ahfslac(:,krow),   ahfswac(:,krow),   ahfsiac(:,krow)   &
            , ahfllac(:,krow),   ahflwac(:,krow),   ahfliac(:,krow)   &
            , evaplac(:,krow),   evapwac(:,krow),   evapiac(:,krow)   &
            , ustrl(:,krow),     ustrw(:,krow),     ustri(:,krow)     &
            , vstrl(:,krow),     vstrw(:,krow),     vstri(:,krow)     &
            , sn(:,krow),        snc(:,krow),       tslm1(:,krow)     &
            , ws(:,krow),        albedo(:,krow),    albedo_vis(:,krow) &
            , albedo_nir(:,krow),alsol(:,krow)                        &
!
            , tke(:,:,krow),     tkem1(:,:,krow),   tkem(:,:,krow)    &
            , aclc(:,:,krow),    emter(:,:,krow)                      &
!
            , aphm1,             apm1,              geom1             &
            , ztvm1,             zvdiffp,           zvmixtau          &
!
            , zcvs,              zcvw                                 &
            , qhfla,             evapot                               &
            , ztslnew,           zwlmx                                &
            , zfrl,              zfrw,              zfri              &
            , loland,            loglac                               &
!
            , xtte(:,:,:,krow)                                        &
!
            , vol(:,:,krow),     vom(:,:,krow),     qte(:,:,krow)     &
            , tte(:,:,krow),     xlte(:,:,krow),    xite(:,:,krow)    &
            , siced(:,krow)                                           &
! Rain and snow over last time step
            , zraind(:),            zsnowd(:)                          &
            , jsswnir(:,krow), jsswdifnir(:,krow)                      &
            , jsswvis(:,krow), jsswdifvis(:,krow)                      &
            , zi0,                    trsol(:,:,krow)                  &
            , ztrflw,           ztrfli(:,krow)      &
            , zsofll,            zsoflw,           zsofli(:,krow)      &
            , trfllac(:,krow),   trflwac(:,krow),  trfliac(:,krow)     &
            , sofllac(:,krow),   soflwac(:,krow),  sofliac(:,krow)     &
            , alsoi(:,krow),    alsow(:,krow),     radtemp(:,krow)     &
            , sni(:,krow)                                              &
            , ahfice(:,krow)                                           &
            , fluxres(:,krow)                                          &
            , qres(:,krow)                                             &
            , ahfcon(:,krow)                                           &
            , ahfres(:,krow)                                           &
            , zteffl4 , ztsnew , tsurf(:,krow)                         &
            , seaice(:,krow)                                           &
            , zradtemp_old, runoff(1:nbdim,krow), drain(1:nbdim,krow)  &
            , zalac, zros_hd, zdrain_hd                                &
            , co2_flux_ocean(1:nbdim,krow), co2_flux_land(1:nbdim,krow)&
            , co2_flux(1:nbdim,krow), co2_flux_corr_acc1(1:nbdim,krow) &
!---wiso-code
            , lwiso, nwiso)
!---wiso-code-end
  END IF
#else
  IF (lwiso) THEN 
    CALL vdiff (nproma, nbdim, ktdia, nlev, nlevm1, nlevp1, ntrac     &
            , krow                                                    &
!
            , xtm1(:,:,:,krow)                                        &
!
            , qm1(1,1,krow),     tm1(1,1,krow),     um1(1,1,krow)     &
            , vm1(1,1,krow),     xlm1(1,1,krow),    xim1(1,1,krow)    &
            , xvar(1,1,krow)                                          &
!
            , ahfl(1,krow),      ahfs(1,krow),      az0(1,krow)       &
            , dew2(1,krow),      evap(1,krow),      forest(1,krow)    &
            , temp2(1,krow),     t2max(1,krow)                        &
            , t2min(1,krow),     wind10w(1,krow),   vdis(1,krow)      &
            , u10(1,krow),       v10(1,krow),       ustr(1,krow)      &
            , vstr(1,krow),      wimax(1,krow),     wind10(1,krow)    &
            , wsmx(1,krow),      vlt(1,krow)                          &
            , vgrat(1,krow)                                           &
            , tsl(1,krow),       tsw(1,krow),       tsi(1,krow)       &
            , ocu(1,krow),       ocv(1,krow)                          &
            , az0l(1,krow),      az0w(1,krow),      az0i(1,krow)      &
            , zhfsl,             zhfsw,             zhfsi             &
            , zhflw                                                   &
            , zevapw,            zevapi                               &
            , ahfslac(1,krow),   ahfswac(1,krow),   ahfsiac(1,krow)   &
            , ahfllac(1,krow),   ahflwac(1,krow),   ahfliac(1,krow)   &
            , evaplac(1,krow),   evapwac(1,krow),   evapiac(1,krow)   &
            , ustrl(1,krow),     ustrw(1,krow),     ustri(1,krow)     &
            , vstrl(1,krow),     vstrw(1,krow),     vstri(1,krow)     &
            , sn(1,krow),        snc(1,krow),       tslm1(1,krow)     &
            , ws(1,krow),        albedo(1,krow),    albedo_vis(1,krow) &
            , albedo_nir(1,krow),alsol(1,krow)                        &
!
            , tke(1,1,krow),     tkem1(1,1,krow),   tkem(1,1,krow)    &
            , aclc(1,1,krow),    emter(1,1,krow)                      &
!
            , aphm1,             apm1,              geom1             &
            , ztvm1,             zvdiffp,           zvmixtau          &
!
            , zcvs,              zcvw                                 &
            , qhfla,             evapot                               &
            , ztslnew,           zwlmx                                &
            , zfrl,              zfrw,              zfri              &
            , loland,            loglac                               &
!
            , xtte(1,1,1,krow)                                        &
!
            , vol(1,1,krow),     vom(1,1,krow),     qte(1,1,krow)     &
            , tte(1,1,krow),     xlte(1,1,krow),    xite(1,1,krow)    &
            , siced(1,krow)                                           &
! Rain and snow over last time step
            , zraind(1),            zsnowd(1)                          &
            , jsswnir(1,krow), jsswdifnir(1,krow)                      &
            , jsswvis(1,krow), jsswdifvis(1,krow)                      &
            , zi0,                    trsol(1,1,krow)                  &
            , ztrflw,           ztrfli(1,krow)                         &
            , zsofll,            zsoflw,           zsofli(1,krow)      &
            , trfllac(1,krow),   trflwac(1,krow),  trfliac(1,krow)     &
            , sofllac(1,krow),   soflwac(1,krow),  sofliac(1,krow)     &
            , alsoi(1,krow),    alsow(1,krow),     radtemp(1,krow)     &
            , sni(1,krow)                                              &
            , ahfice(1,krow)                                           &
            , fluxres(1,krow)                                          &
            , qres(1,krow)                                             &
            , ahfcon(1,krow)                                           &
            , ahfres(1,krow)                                           &
            , zteffl4 , ztsnew , tsurf(1,krow)                         &
            , seaice(1,krow)                                           &
            , zradtemp_old, runoff(1,krow), drain(1,krow)              &
            , zalac, zros_hd, zdrain_hd                                &
            , co2_flux_ocean(1,krow), co2_flux_land(1,krow)            &
            , co2_flux(1,krow), co2_flux_corr_acc1(1,krow)             &
!---wiso-code
            , lwiso, nwiso                                            &
            , wisoqm1(1,1,1,krow), wisoxlm1(1,1,1,krow), wisoxim1(1,1,1,krow)     &
            , wisoevap(1,1,krow)                                                  &
            , zwisoevapw, zwisoevapi                                              &
            , wisoevaplac(1,1,krow), wisoevapwac(1,1,krow), wisoevapiac(1,1,krow) &
            , wisosn(1,1,krow), wisows(1,1,krow)                                  &
            , wisosw_d(1,1,krow)                                                  &
            , wisoqhfla, wisoevapot                                               &
            , wisoqte(1,1,1,krow), wisoxlte(1,1,1,krow), wisoxite(1,1,1,krow)     &
            , zwisoraind(1,1), zwisosnowd(1,1)                                    &
            , wisorunoff(1,1,krow), wisodrain(1,1,krow)                           &
            , zwisoalac, zwisoros_hd, zwisodrain_hd)
!---wiso-code-end
  ELSE
    CALL vdiff (nproma, nbdim, ktdia, nlev, nlevm1, nlevp1, ntrac     &
            , krow                                                    &
!
            , xtm1(:,:,:,krow)                                        &
!
            , qm1(1,1,krow),     tm1(1,1,krow),     um1(1,1,krow)     &
            , vm1(1,1,krow),     xlm1(1,1,krow),    xim1(1,1,krow)    &
            , xvar(1,1,krow)                                          &
!
            , ahfl(1,krow),      ahfs(1,krow),      az0(1,krow)       &
            , dew2(1,krow),      evap(1,krow),      forest(1,krow)    &
            , temp2(1,krow),     t2max(1,krow)                        &
            , t2min(1,krow),     wind10w(1,krow),   vdis(1,krow)      &
            , u10(1,krow),       v10(1,krow),       ustr(1,krow)      &
            , vstr(1,krow),      wimax(1,krow),     wind10(1,krow)    &
            , wsmx(1,krow),      vlt(1,krow)                          &
            , vgrat(1,krow)                                           &
            , tsl(1,krow),       tsw(1,krow),       tsi(1,krow)       &
            , ocu(1,krow),       ocv(1,krow)                          &
            , az0l(1,krow),      az0w(1,krow),      az0i(1,krow)      &
            , zhfsl,             zhfsw,             zhfsi             &
            , zhflw                                                   &
            , zevapw,            zevapi                               &
            , ahfslac(1,krow),   ahfswac(1,krow),   ahfsiac(1,krow)   &
            , ahfllac(1,krow),   ahflwac(1,krow),   ahfliac(1,krow)   &
            , evaplac(1,krow),   evapwac(1,krow),   evapiac(1,krow)   &
            , ustrl(1,krow),     ustrw(1,krow),     ustri(1,krow)     &
            , vstrl(1,krow),     vstrw(1,krow),     vstri(1,krow)     &
            , sn(1,krow),        snc(1,krow),       tslm1(1,krow)     &
            , ws(1,krow),        albedo(1,krow),    albedo_vis(1,krow) &
            , albedo_nir(1,krow),alsol(1,krow)                        &
!
            , tke(1,1,krow),     tkem1(1,1,krow),   tkem(1,1,krow)    &
            , aclc(1,1,krow),    emter(1,1,krow)                      &
!
            , aphm1,             apm1,              geom1             &
            , ztvm1,             zvdiffp,           zvmixtau          &
!
            , zcvs,              zcvw                                 &
            , qhfla,             evapot                               &
            , ztslnew,           zwlmx                                &
            , zfrl,              zfrw,              zfri              &
            , loland,            loglac                               &
!
            , xtte(1,1,1,krow)                                        &
!
            , vol(1,1,krow),     vom(1,1,krow),     qte(1,1,krow)     &
            , tte(1,1,krow),     xlte(1,1,krow),    xite(1,1,krow)    &
            , siced(1,krow)                                           &
! Rain and snow over last time step
            , zraind(1),            zsnowd(1)                          &
            , jsswnir(1,krow), jsswdifnir(1,krow)                      &
            , jsswvis(1,krow), jsswdifvis(1,krow)                      &
            , zi0,                    trsol(1,1,krow)                  &
            , ztrflw,           ztrfli(1,krow)                         &
            , zsofll,            zsoflw,           zsofli(1,krow)      &
            , trfllac(1,krow),   trflwac(1,krow),  trfliac(1,krow)     &
            , sofllac(1,krow),   soflwac(1,krow),  sofliac(1,krow)     &
            , alsoi(1,krow),    alsow(1,krow),     radtemp(1,krow)     &
            , sni(1,krow)                                              &
            , ahfice(1,krow)                                           &
            , fluxres(1,krow)                                          &
            , qres(1,krow)                                             &
            , ahfcon(1,krow)                                           &
            , ahfres(1,krow)                                           &
            , zteffl4 , ztsnew , tsurf(1,krow)                         &
            , seaice(1,krow)                                           &
            , zradtemp_old, runoff(1,krow), drain(1,krow)              &
            , zalac, zros_hd, zdrain_hd                                &
            , co2_flux_ocean(1,krow), co2_flux_land(1,krow)            &
            , co2_flux(1,krow), co2_flux_corr_acc1(1,krow)             &
!---wiso-code
            , lwiso, nwiso)
!---wiso-code-end
  END IF

#endif

!
!
  IF (ltdiag) THEN
! store VDIFF increment
    pdiga(1:nglon,:, 3,krow) = pdiga(1:nglon,:, 3,krow) + vom(1:nproma,:,krow)
    pdiga(1:nglon,:, 8,krow) = pdiga(1:nglon,:, 8,krow) + vol(1:nproma,:,krow)
    pdiga(1:nglon,:,16,krow) = pdiga(1:nglon,:,16,krow) + tte(1:nproma,:,krow)

! prepare next fields for RADHEAT
    pdiga(1:nglon,:,15,krow) = pdiga(1:nglon,:,15,krow) - tte(1:nproma,:,krow)
  ENDIF
!
  IF (ntrac.GT.0) THEN
     CALL xtdriver1 (nproma, nbdim, nlev, nlevp1, ntrac, aphp1(:,:),   &
                     app1(:,:), tte(:,:,krow), tm1(:,:,krow),          &
                     xtm1(:,:,:,krow), xtte(:,:,:,krow))
  ENDIF
!
!------------------------------------------------------------------------------
!
!*            ADD RADIATION TENDENCIES EVERY TIME STEP.
!
420 CONTINUE
!
  IF (lrad) THEN
#ifdef ARGCHECK
     IF (.NOT. ldblrad) THEN
        ALLOCATE (emter1(nbdim,nlevp1,krow))
        emter1(:,:,:)=0._dp
        ALLOCATE (trsol1(nbdim,nlevp1,krow))
        trsol1(:,:,:)=0._dp
        ALLOCATE (emtef01(nbdim,nlevp1,krow))
        emtef01(:,:,:)=0._dp
        ALLOCATE (trsof01(nbdim,nlevp1,krow))
        trsof01(:,:,:)=0._dp
        ALLOCATE (netht_sw(nbdim,nlevp1,krow))
        netht_sw(:,:,:)=0._dp
        ALLOCATE (netht_lw(nbdim,nlevp1,krow))
        netht_lw(:,:,:)=0._dp
     END IF
  CALL radheat (nproma, nbdim, nlev, nlevp1,                           &
          krow,                                                        &
          ldblrad,                                                     &
          zi0,                                                         &
          tm1(:,:,krow)  ,   qm1(:,:,krow),                            &
          trsof(:,:,krow),   trsol(:,:,krow),                          &
          emtef(:,:,krow),   emter(:,:,krow),                          &
          emtef0(:,:,krow),  trsof0(:,:,krow),                         &
          emter1(:,:,krow),  trsol1(:,:,krow),                         &
          emtef01(:,:,krow),  trsof01(:,:,krow),                       &
          netht_sw(:,:,krow), netht_lw(:,:,krow),                      &
          srad0(:,krow),     srads(:,krow),                            &
          sradl(:,krow),     srafl(:,krow),                            &
          srad0u(:,krow),    sradsu(:,krow),                           &
          sraf0(:,krow),     srafs(:,krow),                            &
          srad0d(:,krow),                                              &
          trad0(:,krow),     trads(:,krow),                            &
          tradl(:,krow),     trafl(:,krow),                            &
          traf0(:,krow),     trafs(:,krow),                            &
          tradsu(:,krow),                                              &
          albedo(:,krow),                                              &
          aphm1,             apm1,                                     &
          ztslnew,           tte(:,:,krow),                            &
          zradtemp_old,      ztsnew  )
     IF (.NOT. ldblrad) DEALLOCATE (emter1,trsol1,emtef01,trsof01)
     IF (.NOT. ldblrad) DEALLOCATE (netht_sw,netht_lw)
#else
 CALL radheat (nproma, nbdim, nlev, nlevp1,                           &
          krow,                                                        &
          ldblrad,                                                     &
          zi0,                                                         &
          tm1(1,1,krow)  ,   qm1(1,1,krow),                            &
          trsof(1,1,krow),   trsol(1,1,krow),                          &
          emtef(1,1,krow),   emter(1,1,krow),                          &
          emtef0(1,1,krow),  trsof0(1,1,krow),                         &
          emter1(1,1,krow),  trsol1(1,1,krow),                         &
          emtef01(1,1,krow),  trsof01(1,1,krow),                       &
          netht_sw(1,1,krow), netht_lw(1,1,krow),                      &
          srad0(1,krow),     srads(1,krow),                            &
          sradl(1,krow),     srafl(1,krow),                            &
          srad0u(1,krow),    sradsu(1,krow),                           &
          sraf0(1,krow),     srafs(1,krow),                            &
          srad0d(1,krow),                                              &
          trad0(1,krow),     trads(1,krow),                            &
          tradl(1,krow),     trafl(1,krow),                            &
          traf0(1,krow),     trafs(1,krow),                            &
          tradsu(1,krow),                                              &
          albedo(1,krow),                                              &
          aphm1,             apm1,                                     &
          ztslnew,           tte(1,1,krow),                            &
          zradtemp_old,      ztsnew  )
#endif
  ELSE
     ! lrad=.FALSE.
     ! --> no radiative effect at the surface
     ztrflw(:)=0._dp
     ztrfli(:,krow)=0._dp
     zsofll(:)=0._dp
     zsoflw(:)=0._dp
     zsofli(:,krow)=0._dp
     ! --> lake, ml_ocean, licetemp and sicetemp get zero fluxes 
  END IF
!
  IF (ltdiag) THEN
! store RADHEAT increment
    pdiga(1:nglon,:,15,krow) = pdiga(1:nglon,:,15,krow) + tte(1:nproma,:,krow)

! prepare next fields for GWDRAG
    pdiga(1:nglon,:, 4,krow) = pdiga(1:nglon,:, 4,krow) - vom(1:nproma,:,krow)
    pdiga(1:nglon,:, 9,krow) = pdiga(1:nglon,:, 9,krow) - vol(1:nproma,:,krow)
    pdiga(1:nglon,:,17,krow) = pdiga(1:nglon,:,17,krow) - tte(1:nproma,:,krow)
  ENDIF
!
!
!     ------------------------------------------------------------
!
!        ***  GRAVITY WAVE DRAG PARAMETERISATION  ***
!
  IF (lmidatm) THEN
    IF (lresume) THEN
       aprflux(1:nproma,krow) = 0._dp
       aprfluxm(1:nproma,krow) = 0._dp
    ENDIF
 

    CALL gwspectrum ( krow, nproma,  nbdim,  nlev,                             &
                  aphm1,            apm1,                                      &
                  tm1(:,:,krow),    um1(:,:,krow),   vm1(:,:,krow),            &
                  aprflux(:,krow),                                             &
                  tte(:,:,krow),    vol(:,:,krow),   vom(:,:,krow) )
!
  END IF

  IF (lgwdrag) THEN

    CALL ssodrag ( nproma,  nbdim,  nlev,                                      &
                  aphm1,            apm1,            geom1,                    &
                  tm1(:,:,krow),    um1(:,:,krow),   vm1(:,:,krow),            &
                  oromea(:,krow),   orostd(:,krow),  orosig(:,krow),           &
                  orogam(:,krow),                                              &
                  orothe(:,krow),   oropic(:,krow),  oroval(:,krow),           &
                  ustrgw(:,krow),   vstrgw(:,krow),  vdisgw(:,krow),           &
                  tte(:,:,krow),    vol(:,:,krow),   vom(:,:,krow) )

  END IF

  IF (ltdiag) THEN
! store GWDRAG increment
    pdiga(1:nglon,:, 4,krow) = pdiga(1:nglon,:, 4,krow) + vom(1:nproma,:,krow)
    pdiga(1:nglon,:, 9,krow) = pdiga(1:nglon,:, 9,krow) + vol(1:nproma,:,krow)
    pdiga(1:nglon,:,17,krow) = pdiga(1:nglon,:,17,krow) + tte(1:nproma,:,krow)

! prepare next fields for CUCALL
    pdiga(1:nglon,:, 5,krow) = pdiga(1:nglon,:, 5,krow) - vom(1:nproma,:,krow)
    pdiga(1:nglon,:,10,krow) = pdiga(1:nglon,:,10,krow) - vol(1:nproma,:,krow)
    pdiga(1:nglon,:,18,krow) = pdiga(1:nglon,:,18,krow) - tte(1:nproma,:,krow)
  ENDIF
!
!     ------------------------------------------------------------------
!
!*        6.    CONVECTION PARAMETERISATION.
!               ---------- -----------------
!
600 CONTINUE
!
!
!*        6.3    COMPUTE *T* AND *Q* TENDENCIES BY MOIST CONVECTION.
!*                AND SHALLOW CONVECTION.
!
630 CONTINUE
!
!
   itype(1:nproma)=0  !!!!!
!
!
!
!*         6.3.1   INITIALIZE ARRAYS FOR CONVECTIVE PRECIPITATION
!*                 AND COPY ARRAYS FOR CONVECTIVE CLOUD PARAMETERS
!*                 -----------------------------------------------
!
631 CONTINUE
!
!
   xtec(1:nproma,:,krow)=0._dp  !!!!! 
!
  DO 632 jl=1,nproma
     rsfc(jl)=0._dp
     ssfc(jl)=0._dp
632 END DO

!---wiso-code
  IF (lwiso) THEN

   wisoxtec(1:nproma,:,:,krow)=0._dp

  DO jt=1,nwiso   
   DO jl=1,nproma
     wisorsfc(jl,jt)=0._dp
     wisossfc(jl,jt)=0._dp
   END DO
  END DO

  END IF
!---wiso-code-end
!
!*         6.3.3   CALL SUBROUTINE CUCALL FOR CUMULUS PARAMETERIZATION
!                  ---------------------------------------------------
!
633 CONTINUE
!
!
  IF (loconv) THEN
!
#ifdef ARGCHECK
  IF (lwiso) THEN 
     CALL cucall(nproma, nbdim, nlev, nlevp1, nlevm1, ilab,            &
                ntrac,                                                 &
                xtm1(:,:,:,krow), xtte(:,:,:,krow),                    &
                tm1(:,:,krow),    qm1(:,:,krow),    um1(:,:,krow),     &
                vm1(:,:,krow),    xlm1(:,:,krow),   xim1(:,:,krow),    &
                tte(:,:,krow),    qte(:,:,krow),    vom(:,:,krow),     &
                vol(:,:,krow),    xlte(:,:,krow),   xite(:,:,krow),    &
                vervel(:,:,krow), xtec(:,:,krow),                      &
                zqtec,            qhfla,                               &
                app1,             aphp1,            geom1,             &
                rsfc,             ssfc,                                &
                aprc(:,krow),     aprs(:,krow),                        &
                itype,            loland,                              &
                topmax(:,krow),                                        &                
!---wiso-code
                lwiso, nwiso,                                          &
                wisoqm1(:,:,:,krow),                                   &
                wisoxlm1(:,:,:,krow), wisoxim1(:,:,:,krow),            &
                wisoqte(:,:,:,krow),                                   &
                wisoxlte(:,:,:,krow), wisoxite(:,:,:,krow),            &
                wisoxtec(:,:,:,krow),                                  &
                zwisoqtec,                                             &
                wisorsfc,             wisossfc,                        &
                wisoaprc(:,:,krow),   wisoaprs(:,:,krow))
!---wiso-code-end
  ELSE
     CALL cucall(nproma, nbdim, nlev, nlevp1, nlevm1, ilab,            &
                ntrac,                                                 &
                xtm1(:,:,:,krow), xtte(:,:,:,krow),                    &
                tm1(:,:,krow),    qm1(:,:,krow),    um1(:,:,krow),     &
                vm1(:,:,krow),    xlm1(:,:,krow),   xim1(:,:,krow),    &
                tte(:,:,krow),    qte(:,:,krow),    vom(:,:,krow),     &
                vol(:,:,krow),    xlte(:,:,krow),   xite(:,:,krow),    &
                vervel(:,:,krow), xtec(:,:,krow),                      &
                zqtec,            qhfla,                               &
                app1,             aphp1,            geom1,             &
                rsfc,             ssfc,                                &
                aprc(:,krow),     aprs(:,krow),                        &
                itype,            loland,                              &
                topmax(:,krow),                                        &
!---wiso-code
                lwiso, nwiso)
!---wiso-code-end
  END IF
#else
  IF (lwiso) THEN 
     CALL cucall(nproma, nbdim, nlev, nlevp1, nlevm1, ilab,            &
                ntrac,                                                 &
                xtm1(1,1,1,krow), xtte(1,1,1,krow),                    &
                tm1(1,1,krow),    qm1(1,1,krow),    um1(1,1,krow),     &
                vm1(1,1,krow),    xlm1(1,1,krow),   xim1(1,1,krow),    &
                tte(1,1,krow),    qte(1,1,krow),    vom(1,1,krow),     &
                vol(1,1,krow),    xlte(1,1,krow),   xite(1,1,krow),    &
                vervel(1,1,krow), xtec(1,1,krow),                      &
                zqtec,            qhfla,                               &
                app1,             aphp1,            geom1,             &
                rsfc,             ssfc,                                &
                aprc(1,krow),     aprs(1,krow),                        &
                itype,            loland,                              &
                topmax(1,krow),                                        &
!---wiso-code
                lwiso, nwiso,                                          &
                wisoqm1(1,1,1,krow),                                   &
                wisoxlm1(1,1,1,krow), wisoxim1(1,1,1,krow),            &
                wisoqte(1,1,1,krow),                                   &
                wisoxlte(1,1,1,krow), wisoxite(1,1,1,krow),            &
                wisoxtec(1,1,1,krow),                                  &
                zwisoqtec,                                             &
                wisorsfc,             wisossfc,                        &
                wisoaprc(1,1,krow),   wisoaprs(1,1,krow))
!---wiso-code-end
  ELSE
     CALL cucall(nproma, nbdim, nlev, nlevp1, nlevm1, ilab,            &
                ntrac,                                                 &
                xtm1(1,1,1,krow), xtte(1,1,1,krow),                    &
                tm1(1,1,krow),    qm1(1,1,krow),    um1(1,1,krow),     &
                vm1(1,1,krow),    xlm1(1,1,krow),   xim1(1,1,krow),    &
                tte(1,1,krow),    qte(1,1,krow),    vom(1,1,krow),     &
                vol(1,1,krow),    xlte(1,1,krow),   xite(1,1,krow),    &
                vervel(1,1,krow), xtec(1,1,krow),                      &
                zqtec,            qhfla,                               &
                app1,             aphp1,            geom1,             &
                rsfc,             ssfc,                                &
                aprc(1,krow),     aprs(1,krow),                        &
                itype,            loland,                              &
                topmax(1,krow),                                        &
!---wiso-code
                lwiso, nwiso)
!---wiso-code-end
  END IF
#endif
!
     DO jl=kidia,kfdia       
       rtype(jl,krow)=REAL(itype(jl),dp)
     END DO
!
  ELSE
!       NECESSARY COMPUTATIONS IF MASSFLUX IS BY-PASSED
!
     ilab(1:nproma,1:nlev)=0
!
  ENDIF
!
  IF (ltdiag) THEN
! store CUCALL increment (massflux)
    pdiga(1:nglon,:, 5,krow) = pdiga(1:nglon,:, 5,krow) + vom(1:nproma,:,krow)
    pdiga(1:nglon,:,10,krow) = pdiga(1:nglon,:,10,krow) + vol(1:nproma,:,krow)
    pdiga(1:nglon,:,18,krow) = pdiga(1:nglon,:,18,krow) + tte(1:nproma,:,krow)

! prepare next fields for COND
    pdiga(1:nglon,:,19,krow) = pdiga(1:nglon,:,19,krow) - tte(1:nproma,:,krow)
  ENDIF
!
!     ------------------------------------------------------------
!
!*       7.    LARGE SCALE CONDENSATION.
!              ----- ----- -------------
!
700 CONTINUE
!
!
  IF(locond) THEN
!
!     POLE FILTER FOR TENDENCIES
!
     IF (iadvec == spitfire .and. .not. ldc%col_1d)                    & 
!---wiso-code
! replace original routine pole_filter by new one including water isotopes
! (new routine is defined in mo_spitfire)
!       CALL pole_filter(tte(:,:,krow), qte(:,:,krow), krow)
       CALL pole_filter_wiso(tte(:,:,krow), qte(:,:,krow), wisoqte(:,:,:,krow), krow)

! attention - water isotopes not included in column-model code, yet!       
     IF (lcotra) CALL get_col_pol (tte(:,:,krow), qte(:,:,krow), krow)
!---wiso-code-end

     IF (ltimer) CALL timer_start(timer_cloud)
#ifdef ARGCHECK
  IF (lwiso) THEN
  
     CALL cloud( nproma, nbdim, ktdia, nlev, nlevp1                   &
               , aphm1,             vervel(:,:,krow)                  &
               , apm1,              app1,             acdnc(:,:,krow) &
               , qm1(:,:,krow),     tm1(:,:,krow),    ztvm1           &
               , xlm1(:,:,krow),    xim1(:,:,krow),   xtec(:,:,krow)  &
               , xvar(:,:,krow),    xskew(:,:,krow),  zqtec           &
               , zbetaa,            zbetab                            &
               , zvdiffp,           zhmixtau,         zvmixtau        &
               , geom1,             zbetass                           &
               , invb                                                 &
               , aclc(:,:,krow),    aclcac(:,:,krow)                  &
               , relhum(:,:,krow)                                     &
               , aclcov(:,krow),    aprl(:,krow),     qvi(:,krow)     &
               , xlvi(:,krow),      xivi(:,krow)                      &
               , ssfl,              rsfl                              &
               , qte(:,:,krow),     tte(:,:,krow)                     &
               , xlte(:,:,krow),    xite(:,:,krow)                    &
               , aprs(:,krow)                                         &
!---wiso-code
               , lwiso, nwiso                                         &
               , wisoqm1(:,:,:,krow)                                  &
               , wisoxlm1(:,:,:,krow), wisoxim1(:,:,:,krow)           &
               , wisoxtec(1,1,1,krow)                                 &
               , zwisoqtec                                            &
               , wisoaprl(:,:,krow),   wisoqvi(:,:,krow)              &
               , wisoxlvi(:,:,krow),   wisoxivi(:,:,krow)             &
               , wisossfl,             wisorsfl                       &               
               , wisoqte(:,:,:,krow)                                  &
               , wisoxlte(:,:,:,krow), wisoxite(:,:,:,krow)           &
               , wisoaprs(:,:,krow))

  ELSE

     CALL cloud( nproma, nbdim, ktdia, nlev, nlevp1                   &
               , aphm1,             vervel(:,:,krow)                  &
               , apm1,              app1,             acdnc(:,:,krow) &
               , qm1(:,:,krow),     tm1(:,:,krow),    ztvm1           &
               , xlm1(:,:,krow),    xim1(:,:,krow),   xtec(:,:,krow)  &
               , xvar(:,:,krow),    xskew(:,:,krow),  zqtec           &
               , zbetaa,            zbetab                            &
               , zvdiffp,           zhmixtau,         zvmixtau        &
               , geom1,             zbetass                           &
               , invb                                                 &
               , aclc(:,:,krow),    aclcac(:,:,krow)                  &
               , relhum(:,:,krow)                                     &
               , aclcov(:,krow),    aprl(:,krow),     qvi(:,krow)     &
               , xlvi(:,krow),      xivi(:,krow)                      &
               , ssfl,              rsfl                              &
               , qte(:,:,krow),     tte(:,:,krow)                     &
               , xlte(:,:,krow),    xite(:,:,krow)                    &
               , aprs(:,krow)                                         &
!---wiso-code
               , lwiso, nwiso)

  END IF
  
!---wiso-code-end
#else

  IF (lwiso) THEN

     CALL cloud( nproma, nbdim, ktdia, nlev, nlevp1                   &
               , aphm1,             vervel(1,1,krow)                  &
               , apm1,              app1,             acdnc(1,1,krow) &
               , qm1(1,1,krow),     tm1(1,1,krow),    ztvm1           &
               , xlm1(1,1,krow),    xim1(1,1,krow),   xtec(1,1,krow)  &
               , xvar(1,1,krow),    xskew(1,1,krow),  zqtec           &
               , zbetaa,            zbetab                            &
               , zvdiffp,           zhmixtau,         zvmixtau        &
               , geom1,             zbetass                           &
               , invb                                                 &
               , aclc(1,1,krow),    aclcac(1,1,krow)                  &
               , relhum(1,1,krow)                                     &
               , aclcov(1,krow),    aprl(1,krow),     qvi(1,krow)     &
               , xlvi(1,krow),      xivi(1,krow)                      &
               , ssfl,              rsfl                              &
               , qte(1,1,krow),     tte(1,1,krow)                     &
               , xlte(1,1,krow),    xite(1,1,krow)                    &
               , aprs(1,krow)                                         &
!---wiso-code
               , lwiso, nwiso                                         &
               , wisoqm1(1,1,1,krow)                                  &
               , wisoxlm1(1,1,1,krow), wisoxim1(1,1,1,krow)           &
               , wisoxtec(1,1,1,krow)                                 &
               , zwisoqtec                                            &
               , wisoaprl(1,1,krow),   wisoqvi(1,1,krow)              &
               , wisoxlvi(1,1,krow),   wisoxivi(1,1,krow)             &
               , wisossfl,             wisorsfl                       &               
               , wisoqte(1,1,1,krow)                                  &
               , wisoxlte(1,1,1,krow), wisoxite(1,1,1,krow)           &
               , wisoaprs(1,1,krow))
!---wiso-code-end

  ELSE
  
     CALL cloud( nproma, nbdim, ktdia, nlev, nlevp1                   &
               , aphm1,             vervel(1,1,krow)                  &
               , apm1,              app1,             acdnc(1,1,krow) &
               , qm1(1,1,krow),     tm1(1,1,krow),    ztvm1           &
               , xlm1(1,1,krow),    xim1(1,1,krow),   xtec(1,1,krow)  &
               , xvar(1,1,krow),    xskew(1,1,krow),  zqtec           &
               , zbetaa,            zbetab                            &
               , zvdiffp,           zhmixtau,         zvmixtau        &
               , geom1,             zbetass                           &
               , invb                                                 &
               , aclc(1,1,krow),    aclcac(1,1,krow)                  &
               , relhum(1,1,krow)                                     &
               , aclcov(1,krow),    aprl(1,krow),     qvi(1,krow)     &
               , xlvi(1,krow),      xivi(1,krow)                      &
               , ssfl,              rsfl                              &
               , qte(1,1,krow),     tte(1,1,krow)                     &
               , xlte(1,1,krow),    xite(1,1,krow)                    &
               , aprs(1,krow)                                         &
!---wiso-code
               , lwiso, nwiso)
!---wiso-code-end
  
  END IF
#endif
     IF (ltimer) CALL timer_stop(timer_cloud)

  ELSE
!
!              NECESSARY COMPUTATIONS IF *CLOUD* IS BY-PASSED.
!
    ssfl(1:nproma) = 0._dp
    rsfl(1:nproma) = 0._dp

!---wiso-code
  IF (lwiso) THEN

    wisossfl(1:nproma,1:nwiso) = 0._dp
    wisorsfl(1:nproma,1:nwiso) = 0._dp

  END IF
!---wiso-code-end
!
    aclc(1:nproma,:,krow) = 0._dp
!
  ENDIF
!
! store COND increment
  IF (ltdiag) &
       pdiga(1:nglon,:, 19,krow) = pdiga(1:nglon,:, 19,krow) + tte(1:nproma,:,krow)
 
  IF (lmidatm) THEN
    DO jl=1,nproma
       aprflux(jl,krow)=rsfl(jl)+ssfl(jl)+rsfc(jl)+ssfc(jl)
    END DO
  END IF
  DO jl=1,nproma
    zprecip(jl)=rsfl(jl)+ssfl(jl)+rsfc(jl)+ssfc(jl)
  END DO

!---wiso-code
  IF (lwiso) THEN

  DO jt=1,nwiso
    DO jl=1,nproma
      zwisoprecip(jl,jt)=wisorsfl(jl,jt)+wisossfl(jl,jt)+wisorsfc(jl,jt)+wisossfc(jl,jt)
    END DO
  END DO

  END IF
!---wiso-code-end
!
!     ----------------------------------------------------------------
!
!*        8.    Computation of new surface values over land points
!
!
  IF (lsurf) THEN   

!
!
!*   10.1  Compute mixed layer ocean physics
!
! JSBACH - deleted 300904 (kalle)
! ATTENTION-290904 (kalle.)
! SUBROUTINE ML_OCEAN is VOID !! (?)
! ATTENTION: fluxres comes from licetemp and NOT from sicetemp
! This ist important for the calculation row within mo_surface !
! Here sicetemp is calculated BEFORE ml_ocean is called from
! physc. WORKAROUND would be the call to ml_ocean within
! mo_surface between licetemp and sicetemp.
! ATTENTION
!!$    IF (lmlo) THEN
!!$      CALL ml_ocean ( nproma                                          &
!!$            , slm(1,krow)                                             &
!!$            , lonorth                                                 &
!!$            , seaice(1,krow),    siced(1,krow),    alake(1,krow)      &
!!$            , tsi(1,krow),       tsw(1,krow)                          &
!!$            , zhflw,             zhfsw,            fluxres(1,krow)    &
!!$            , ztrflw,            zsoflw                               &
!!$            , amlcorr(1,krow),   amlcorac(1,krow), amlheatac(1,krow)  &
!!$            , zevapi,            sni(1,krow),      zcvsi              &
!!$            , ahfres(1,krow),    zfri           )
!!$    END IF 


  END IF

    IF (lnmi_run) THEN
      SELECT CASE(nmi_phase)
      CASE(NMI_ACCU)    ! accumulate tendencies
!DIR$ CONCURRENT
        dh_t(:,:,krow) = dh_t(:,:,krow) + tte(:,:,krow) - buf_t(:,:)
!DIR$ CONCURRENT
        dh_m(:,:,krow) = dh_m(:,:,krow) + vom(:,:,krow) - buf_m(:,:)
!DIR$ CONCURRENT
        dh_l(:,:,krow) = dh_l(:,:,krow) + vol(:,:,krow) - buf_l(:,:)

      CASE(NMI_USE_AVG) ! prepare/reset tendencies
!DIR$ CONCURRENT
        tte(:,:,krow) = buf_t(:,:) + dh_t(:,:,krow)
!DIR$ CONCURRENT
        vom(:,:,krow) = buf_m(:,:) + dh_m(:,:,krow)
!DIR$ CONCURRENT
        vol(:,:,krow) = buf_l(:,:) + dh_l(:,:,krow)

      END SELECT
    END IF

!---wiso-code
  IF (lwiso) THEN
  
! final adjustment of water isotope values of vapour, cloud water and cloud ice

#ifdef ARGCHECK
    CALL physcadjwiso (nproma,nbdim,nlev,nwiso,   &
                       qm1(:,:,krow),       xlm1(:,:,krow),       xim1(:,:,krow),        &
                       qte(:,:,krow),       xlte(:,:,krow),       xite(:,:,krow),        &
                       wisoqm1(:,:,:,krow), wisoxlm1(:,:,:,krow), wisoxim1(:,:,:,krow),  &
                       wisoqte(:,:,:,krow), wisoxlte(:,:,:,krow), wisoxite(:,:,:,krow))
#else
    CALL physcadjwiso (nproma,nbdim,nlev,nwiso,   &
                       qm1(1,1,krow),       xlm1(1,1,krow),       xim1(1,1,krow),        &
                       qte(1,1,krow),       xlte(1,1,krow),       xite(1,1,krow),        &
                       wisoqm1(1,1,1,krow), wisoxlm1(1,1,1,krow), wisoxim1(1,1,1,krow),  &
                       wisoqte(1,1,1,krow), wisoxlte(1,1,1,krow), wisoxite(1,1,1,krow))
#endif

  END IF
!---wiso-code-end

!
!      p-e budget correction for coupling
!
      DO jl=1,nproma
        zqtold(jl)=qtnew(jl,krow)
        qtnew(jl,krow)=0._dp
      END DO
      DO jk=ktdia,nlev
!DIR$ CONCURRENT
        DO jl=1,nproma
          qtnew(jl,krow)=qtnew(jl,krow)+(qm1(jl,jk,krow)              &
                         +xlm1(jl,jk,krow)+xim1(jl,jk,krow))          &
                         *(aphm1(jl,jk+1)-aphm1(jl,jk))/g
        END DO
      END DO

!---wiso-code
  IF (lwiso) THEN

!   p-e budget correction for coupling
!
      DO jt=1,nwiso
        DO jl=1,nproma
          zwisoqtold(jl,jt)=wisoqtnew(jl,jt,krow)
          wisoqtnew(jl,jt,krow)=0._dp
        END DO
      END DO
      DO jt=1,nwiso
        DO jk=ktdia,nlev
!DIR$ CONCURRENT
          DO jl=1,nproma
            wisoqtnew(jl,jt,krow)=wisoqtnew(jl,jt,krow)+(wisoqm1(jl,jk,jt,krow)              &
                                  +wisoxlm1(jl,jk,jt,krow)+wisoxim1(jl,jk,jt,krow))          &
                                  *(aphm1(jl,jk+1)-aphm1(jl,jk))/g
          END DO
        END DO
      END DO

  END IF
!---wiso-code-end

!     
!     Accumulate p-e correction for standard diagnostics
!
!DIR$ CONCURRENT
      DO jl=1,nproma
        apmeb(jl,krow)=apmeb(jl,krow)-(qtnew(jl,krow)-zqtold(jl))     &
                       -(zprecip(jl)+qhfla(jl))*zdtime
      END DO

!---wiso-code
  IF (lwiso) THEN

!     Accumulate p-e correction for standard diagnostics
!
      DO jt=1,nwiso
!DIR$ CONCURRENT
        DO jl=1,nproma
          wisoapmeb(jl,jt,krow)=wisoapmeb(jl,jt,krow)-(wisoqtnew(jl,jt,krow)-zwisoqtold(jl,jt))     &
                         -(zwisoprecip(jl,jt)+wisoqhfla(jl,jt))*zdtime
        END DO
      END DO

  END IF
!---wiso-code-end
!
!     Vertical integral of anthropogenic sulfur burden 
!
      IF(lso4) THEN
        DO jk=ktdia,nlev
          DO jl=1,nproma
            abso4(jl,krow)=abso4(jl,krow)                               &
                    +(so4all(jl,jk,krow)-so4nat(jl,jk,krow))            &
                    *(aphm1(jl,jk+1)-aphm1(jl,jk))/g*zdtime
          END DO
        END DO
      ENDIF
!
!
!!$ HILFESTELLUNG ZUR KOPPLUNG: kalle, 020904.
!!$      zhflw(jl) = ocean%latent_heat_flux_inst(1:nproma,krow)
!!$      zhfsw(jl) = ocean%sensible_heat_flux_inst(1:nproma,krow) 
!!$      ahfice(jl,krow) = ice%ahfice(1:nproma,krow) 
!!$      ztrflw(jl) = ocean%trflw(1:nproma,krow)
!!$      zsoflw(jl)= ocean%soflw(1:nproma,krow) 
!!$      qres(jl,krow) = ice%qres(1:nproma,krow) 
!!$      zevapw(jl) = ocean%evaporation_inst(1:nproma, krow)
!!$      zevapi(jl)= ice%evaporation_inst(1:nproma, krow)
!!$      ustrw(jl,krow) = ocean%u_stress(1:nproma,krow) 
!!$      vstrw(jl,krow) = ocean%v_stress(1:nproma,krow) 
!!$      

      IF (lcouple) THEN
         IF (l_getocean) THEN
            apmebco(:,krow) = 0._dp ! set to zero after coupling
            rain(:,krow) = 0._dp    !       "         "
!---wiso-code
        IF (lwiso) THEN
          wisoapmebco(:,:,krow) = 0._dp ! set to zero after coupling
          wisorain(:,:,krow) = 0._dp    !       "         "
        END IF
!---wiso-code-end
         END IF

         ! Accumulate p-e correction for coupling
!DIR$ CONCURRENT
         DO jl=1,nproma
            apmebco(jl,krow)=apmebco(jl,krow)-(qtnew(jl,krow)-zqtold(jl))  &
                                 -(zprecip(jl)+qhfla(jl))*zdtime
            rain(jl,krow)=rain(jl,krow)+(rsfl(jl)+rsfc(jl))*zdtime
         END DO

!---wiso-code
    IF (lwiso) THEN

!     Accumulate p-e correction for coupling
      DO jt=1,nwiso
!DIR$ CONCURRENT
        DO jl=1,nproma
          wisoapmebco(jl,jt,krow)=wisoapmebco(jl,jt,krow)-(wisoqtnew(jl,jt,krow)-zwisoqtold(jl,jt)) &
                         -(zwisoprecip(jl,jt)+wisoqhfla(jl,jt))*zdtime
          wisorain(jl,jt,krow)=wisorain(jl,jt,krow)+(wisorsfl(jl,jt)+wisorsfc(jl,jt))*zdtime
        END DO
      END DO

    END IF
!---wiso-code-end

         !
         ! collect data needed as input for the ocean model
         !
#ifdef ARGCHECK

  IF (lwiso) THEN

   CALL collect ( nproma                                               &
        , zhflw,          zhfsw,         ahfice(:,krow)                &
        , ztrflw,         zsoflw                                       &
        , qres(:,krow),   zevapw,        zevapi                        &
        , ustrw(:,krow),  vstrw(:,krow), ustri(:,krow), vstri(:,krow)  &
        , alake(:,krow),  slf(:,krow),   seaice(:,krow)                &
        , wind10w(:,krow), co2m1(:,nlev,krow), co2_flux_ocean(:,krow)  &
        , awhea(:,krow),  awsol(:,krow), awfre(:,krow), awust(:,krow)  &
        , awvst(:,krow),  aicon(:,krow), aiqre(:,krow), aifre(:,krow)  &
        , aiust(:,krow),  aivst(:,krow), awsta(:,krow)                 &
        , co2atmos(:,krow), co2flux_cpl(:,krow)                        &
        , rsfc,           ssfc,          rsfl,          ssfl           &
!---wiso-code
        , lwiso, nwiso                                                 &  
        , zwisoevapw,     zwisoevapi                                   &
        , wisoawfre(:,:,krow)                                          &
        , wisoaifre(:,:,krow)                                          &
        , wisorsfc,      wisossfc,       wisorsfl,      wisossfl)
 !---wiso-code-end
       
  ELSE
  
   CALL collect ( nproma                                               &
        , zhflw,          zhfsw,         ahfice(:,krow)                &
        , ztrflw,         zsoflw                                       &
        , qres(:,krow),   zevapw,        zevapi                        &
        , ustrw(:,krow),  vstrw(:,krow), ustri(:,krow), vstri(:,krow)  &
        , alake(:,krow),  slf(:,krow),   seaice(:,krow)                &
        , wind10w(:,krow), co2m1(:,nlev,krow), co2_flux_ocean(:,krow)  &
        , awhea(:,krow),  awsol(:,krow), awfre(:,krow), awust(:,krow)  &
        , awvst(:,krow),  aicon(:,krow), aiqre(:,krow), aifre(:,krow)  &
        , aiust(:,krow),  aivst(:,krow), awsta(:,krow)                 &
        , co2atmos(:,krow), co2flux_cpl(:,krow)                        &
        , rsfc,           ssfc,          rsfl,          ssfl           &
!---wiso-code
        , lwiso, nwiso)
!---wiso-code-end

  END IF
#else

  IF (lwiso) THEN

   CALL collect ( nproma                                               &
        , zhflw,          zhfsw,         ahfice(1,krow)                &
        , ztrflw,         zsoflw                                       &
        , qres(1,krow),   zevapw,        zevapi                        &
        , ustrw(1,krow),  vstrw(1,krow), ustri(1,krow), vstri(1,krow)  &
        , alake(1,krow),  slf(1,krow),   seaice(1,krow)                &
        , wind10w(1,krow), co2m1(1,nlev,krow), co2_flux_ocean(1,krow)  &
        , awhea(1,krow),  awsol(1,krow), awfre(1,krow), awust(1,krow)  &
        , awvst(1,krow),  aicon(1,krow), aiqre(1,krow), aifre(1,krow)  &
        , aiust(1,krow),  aivst(1,krow), awsta(1,krow)                 &
        , co2atmos(1,krow), co2flux_cpl(1,krow)                        &
        , rsfc,           ssfc,          rsfl,          ssfl           &
!---wiso-code
        , lwiso,         nwiso                                         &  
        , zwisoevapw,    zwisoevapi                                    &
        , wisoawfre(1,1,krow)                                          &
        , wisoaifre(1,1,krow)                                          &
        , wisorsfc,      wisossfc,       wisorsfl,      wisossfl)
!---wiso-code-end

  ELSE
  
   CALL collect ( nproma                                               &
        , zhflw,          zhfsw,         ahfice(1,krow)                &
        , ztrflw,         zsoflw                                       &
        , qres(1,krow),   zevapw,        zevapi                        &
        , ustrw(1,krow),  vstrw(1,krow), ustri(1,krow), vstri(1,krow)  &
        , alake(1,krow),  slf(1,krow),   seaice(1,krow)                &
        , wind10w(1,krow), co2m1(1,nlev,krow), co2_flux_ocean(1,krow)  &
        , awhea(1,krow),  awsol(1,krow), awfre(1,krow), awust(1,krow)  &
        , awvst(1,krow),  aicon(1,krow), aiqre(1,krow), aifre(1,krow)  &
        , aiust(1,krow),  aivst(1,krow), awsta(1,krow)                 &
        , co2atmos(1,krow), co2flux_cpl(1,krow)                        &
        , rsfc,           ssfc,          rsfl,          ssfl           &
!---wiso-code
        , lwiso,         nwiso)
!---wiso-code-end
  
  END IF
#endif
        afre_residual(:,krow) = awfre(:,krow) + aifre(:,krow)
!---wiso-code
      IF (lwiso) THEN
        DO jt=1,nwiso
          wisoafre_residual(:,:,krow) = wisoawfre(:,:,krow) + wisoaifre(:,:,krow)
        END DO
     END IF
!---wiso-code-end
      END IF
      !
      !     collect data needed as input for the HD-model 
      !     and glacier calving model
      !    (includes diagnostics of discharge and calving on atmospheric grid)
      !
      IF (lhd) THEN

 IF (lwiso) THEN

     CALL hydrology_collect ( nproma                                   &
                            , aros(:,krow),    adrain(:,krow)          &
                            , apmecal(:,krow)                          &
                            , disch(:,krow),   runtoc(:,krow)          &
                            , zros_hd,         zdrain_hd               &
                            , zalac                                    &
!---wiso-code
                            , lwiso, nwiso                             &  
                            , wisoaros(:,:,krow), wisoadrain(:,:,krow) &
                            , wisoapmecal(:,:,krow)                    &
                            , wisodisch(:,:,krow),wisoruntoc(:,:,krow) &
                            , zwisoros_hd, zwisodrain_hd               &
                            , zwisoalac)
!---wiso-code-end

 ELSE
 
     CALL hydrology_collect ( nproma                                   &
                            , aros(:,krow),    adrain(:,krow)          &
                            , apmecal(:,krow)                          &
                            , disch(:,krow),   runtoc(:,krow)          &
                            , zros_hd,         zdrain_hd               &
                            , zalac                                    &
!---wiso-code
                            , lwiso, nwiso)
!---wiso-code-end

 END IF

      END IF

      IF (ico2idx > 0) THEN
         ! Check for negative CO2 concentration
         IF (MINVAL(xtm1(1:nproma,ktdia:nlev,ico2idx,krow)) <= 0._dp) THEN
            CALL finish('physc','Negative CO2 concentration. Stop.')         
         END IF
         ! Limit change in CO2 concentration to +/- 2% of absolute concentration over ocean
         IF (lco2_2perc) THEN
            DO jk=ktdia,nlev
               DO jl=1,nproma
                  IF (.NOT. loland(jl)) THEN
                     zco2eps = 0.02_dp * xtm1(jl,jk,ico2idx,krow) / ztwodt
                     zco2tend = xtte(jl,jk,ico2idx,krow)
                     xtte(jl,jk,ico2idx,krow) = MIN(MAX(xtte(jl,jk,ico2idx,krow), -zco2eps), zco2eps)
                     IF (xtte(jl,jk,ico2idx,krow) /= zco2tend) THEN
                        CALL message('physc','Limit change in CO2 concentration, jl,jk,krow = '// &
                                     int2string(jl)//', '//int2string(jk)//', '//int2string(krow))
                     END IF
                     co2_burden_corr_acc2(jl,krow) = co2_burden_corr_acc2(jl,krow) &
                                      + (xtte(jl,jk,ico2idx,krow) - zco2tend) * ztwodt &
                                      * zdtime * (aphm1(jl,jk+1)-aphm1(jl,jk)) / g
                  END IF
               END DO
            END DO
         END IF

         ! CO2 nudging
         IF (lco2_nudge) THEN
            co2_burden_nudge(:,:) = 0._dp
            DO jk=ktdia,nlev-2
               DO jl=1,nproma
                  ! Current concentration
                  zco2p1 = xtm1(jl,jk,ico2idx,krow) + xtte(jl,jk,ico2idx,krow) * ztwodt
                  ! Compute nudging correction in concentration (mass mixing ratio)
                  zco2_nudge_corr = co2_nudge_relax * (co2_nudge_mmr - zco2p1)
                  ! Add nudging term (correction) to old concentration => new concentration
                  zco2p1 = zco2p1 + zco2_nudge_corr
                  ! Compute new tendencies from difference between new and old concentration
                  xtte(jl,jk,ico2idx,krow) = (zco2p1 - xtm1(jl,jk,ico2idx,krow)) / ztwodt
                  ! Compute change in burden due to nudging, not accumulated
                  co2_burden_nudge(jl,krow) = co2_burden_nudge(jl,krow) &
                       + zco2_nudge_corr  * (aphm1(jl,jk+1)-aphm1(jl,jk))/g
               END DO
            END DO
            ! Compute change in burden due to nudging, accumulated
            co2_burden_nudge_acc(1:nproma,krow) = co2_burden_nudge_acc(1:nproma,krow) &
                 + zdtime * co2_burden_nudge(1:nproma,krow)
         END IF

      END IF

      
!
! allow submodels to calculate some diagnostics
!

  IF (ntrac>0) THEN
    CALL xtdriver2 (nproma, nbdim, nlev, nlevp1, ntrac, aphp1(:,:),    &
     app1(:,:), tte(:,:,krow), tm1(:,:,krow), xtm1(:,:,:,krow), xtte(:,:,:,krow))
    CALL xtdiagn
  END IF

  ! Accumulate CO2 fluxes and burden
  CALL diag_co2(krow, nproma, nbdim, nlevp1, aphm1(:,:), zfrl(:), zfrw(:), zfri(:))

!
! daily block statistics - AMIP2 global diagnostics
!  
  IF(ldiagamip) THEN
    CALL collect_amip2_diag(nproma,nbdim,nlev,nlevp1,krow                    &
       ,tm1(:,:,krow),um1(:,:,krow),vm1(:,:,krow),aphm1(:,:),geospm(:,krow)  &
       ,ustr(:,krow),ustrgw(:,krow),ustrm(:,krow),ustrgwm(:,krow)       &
       ,tslm1(:,krow),seaicem(:,krow),loland)
  END IF
!
! Save large-scale/convective rain/snow rate for next time step
!
  jrsfl(1:nproma,krow) = rsfl(1:nproma)
  jrsfc(1:nproma,krow) = rsfc(1:nproma)
  jssfl(1:nproma,krow) = ssfl(1:nproma)
  jssfc(1:nproma,krow) = ssfc(1:nproma)

!---wiso-code
 IF (lwiso) THEN

 DO jt=1,nwiso
    jwisorsfl(1:nproma,jt,krow) = wisorsfl(1:nproma,jt)
    jwisorsfc(1:nproma,jt,krow) = wisorsfc(1:nproma,jt)
    jwisossfl(1:nproma,jt,krow) = wisossfl(1:nproma,jt)
    jwisossfc(1:nproma,jt,krow) = wisossfc(1:nproma,jt)
 END DO
 
 END IF
!---wiso-code-end 

!      
!*        9.    RESTORE WINDS AND WIND TENDENCIES.
!
900 CONTINUE
!
!
!
!*        9.2   RESTORE WINDS AND WIND TENDENCIES.
!
920 CONTINUE
!
  DO jlev = 1, nlev
!DIR$ CONCURRENT
    DO jl = 1, nproma
      zcst = sqcst_2d(jl,krow)
      um1(jl,jlev,krow) = zcst*um1(jl,jlev,krow)
      vm1(jl,jlev,krow) = zcst*vm1(jl,jlev,krow)
      vol(jl,jlev,krow) = zcst*vol(jl,jlev,krow)
      vom(jl,jlev,krow) = zcst*vom(jl,jlev,krow)
    END DO
  END DO

!
!     ------------------------------------------------------------
!
!*       10.    RELEASE SPACE.
!               ------- ------
!
1200 CONTINUE
!
!     ------------------------------------------------------------
!
  RETURN
END SUBROUTINE physc
