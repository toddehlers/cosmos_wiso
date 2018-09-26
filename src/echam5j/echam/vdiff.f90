SUBROUTINE vdiff ( kproma, kbdim, ktdia, klev, klevm1, klevp1, ktrac   &
         , krow                                                        &
!-----------------------------------------------------------------------
! - 3D from mo_memory_g1a
         , pxtm1                                                       &
! - 2D from mo_memory_g1a
         , pqm1,       ptm1,       pum1                                &
         , pvm1,       pxlm1,      pxim1                               &
         , pxvar                                                       &
! - 1D 
         , pahfl,      pahfs,      paz0                                &
         , pdew2,      pevap,      pforest                             &
         , ptemp2,     pt2max                                          &
         , pt2min,     pwind10w,   pvdis                               &
         , pu10,       pv10,       pustr                               &
         , pvstr,      pwimax,     pwind10                             &
         , pwsmx,      pvlt                                            &
         , pvgrat                                                      &
         , ptsl,       ptsw,       ptsi                                &
         , pocu,       pocv                                            &
         , paz0l,      paz0w,      paz0i                               &
         , pahfsl,     pahfsw,     pahfsi                              &
         , pahflw                                                      &
         , pevapw,     pevapi                                          &
         , pahfslac,   pahfswac,   pahfsiac                            &
         , pahfllac,   pahflwac,   pahfliac                            &
         , pevaplac,   pevapwac,   pevapiac                            &
         , pustrl,     pustrw,     pustri                              &
         , pvstrl,     pvstrw,     pvstri                              &
         , psn,        psnc,       ptslm1                              &
         , pws,        palbedo,    palbedo_vis                         &
         , palbedo_nir,palsol                                          &
! - 2D from mo_memory_g3b
         , ptke,       ptkem1,     ptkem                               &
         , paclc,      pemter                                          &
! - 2D within physics only
         , paphm1,     papm1,      pgeom1                              &
         , ptvm1,      pvdiffp,    pvmixtau                            &
! - 1D within physics only.
         , pcvs,       pcvw                                            &
         , pqhfla,     pevapot                                         &
         , ptslnew,    pwlmx                                           &
         , pfrl,       pfrw,       pfri                                &
         , lpland,     lpglac                                          &
! - Tendencies
! - 3D
         , pxtte                                                       &
! - 2D
         , pvol,       pvom,       pqte                                &
         , ptte,       pxlte,      pxite                               &
         , psiced                                                      &
! - jsbach addons for whole surface calculations
         , praind,      psnowd                                         &
         , jsswnir,  jsswdifnir                                        &
         , jsswvis,  jsswdifvis                                        &
         , pi0 , ptrsol                                                &
         , ptrflw,      ptrfli                                         &
         , psofll,       psoflw,      psofli                           &
         , ptrfllac,     ptrflwac,    ptrfliac                         &
         , psofllac,     psoflwac,    psofliac                         &
         , palsoi,       palsow,      zth_new                          &
         , psni                                                        &
         , pahfice                                                     &
         , pfluxres                                                    &
         , pqres                                                       &
         , pahfcon                                                     &
         , pahfres                                                     &
         , pzteffl4                                                    &
         , pztsnew, ptsurf                                             &
         , pseaice, pradtemp_old, prunoff, pdrain                      &
         , palac, pros_hd, pdrain_hd                                   &
         , pco2fluxo, pco2fluxl, pco2flux, pco2flux_corr               &                
!---wiso-code

         , lwiso, kwiso                                                &
! - 2D
         , pwisoqm1,   pwisoxlm1,  pwisoxim1                           &
! - 1D
         , pwisoevap                                                   &
         , pwisoevapw,  pwisoevapi                                     &
         , pwisoevaplac, pwisoevapwac, pwisoevapiac                    &
         , pwisosn,     pwisows                                        &
         , pwisosw_d                                                   &
! - 1D within physics only.
         , pwisoqhfla,  pwisoevapot                                    &
! - Tendencies
! - 2D
         , pwisoqte,    pwisoxlte,  pwisoxite                          &
! - jsbach addons for whole surface calculations
         , pwisoraind,  pwisosnowd                                     &
         , pwisorunoff, pwisodrain                                     &
         , pwisoalac,   pwisoros_hd, pwisodrain_hd)
!---wiso-code-end

!
!**** *vdiff* - does the vertical exchange of u, v, t, q, xl, xi and xt
!               by turbulence.
!
!
!     Subject.
!
!       This routine computes the physical tendencies of the seven
!   prognostic variables u, v, t, q, xl, xi and xt due to the vertical
!   exchange by turbulent (= non-moist convective) processes.
!   These tendencies are obtained as the difference between
!   the result of an implicit time-step starting from values at t-1
!   and these t-1 (???) values.
!   all the diagnostic computations (exchange coefficients, ...) are
!   done from the t-1 values. As a by-product the roughness length
!   over sea is updated accordingly to the *charnock formula. heat and
!   moisture surface fluxes and their derivatives against ts and ws,
!   later to be used for soil processes treatment, are also
!   computed as well as a stability value to be used as a diagnostic
!   of the depth of the well mixed layer in convective computations.
!
!**   Interface.
!     ----------
!
!          *vdiff* is called from *physc*.
!
!      Arguments.
!      ----------
!
!  - 3d from mo_memory_g1a
!
!  pxtm1    : tracer variables (t-dt)
!
!  - 2d from mo_memory_g1a
!
!  pqm1     : humidity (t-dt)
!  ptm1     : temperature (t-dt)
!  pum1     : zonal wind (t-dt)
!  pvm1     : meridional wind (t-dt)
!  pxlm1    : cloud water (t-dt)
!  pxim1    : cloud ice (t-dt)
!  pxvar    : distribution width (b-a) (t-dt)

! - 2d from mo_memory_g3
!
!  paclc    : cloud cover
!
!  ptke     : turbulent kinetic energy at t+dt (unfiltered)
!  ptkem    :            "             at t    (unfiltered)
!  ptkem1   :            "             at t-dt   (filtered)
!
!  - 1d from mo_memory_g3
!
!  ptsl     : surface temperature over land
!  ptsw     :             "       over water
!  ptsi     :             "       over ice
!
!  pocu     : ocean u-velocity
!  pocv     : ocean v-velocity
!
!  pahfs    : surface sensible heat flux (accumulated)
!  pahfsl   :             "              over land
!  pahfsw   :             "              over water
!  pahfsi   :             "              over ice
!
!  pahfl    : surface latent heat flux   (accumulated)
!  pahfll   :             "              over land
!  pahflw   :             "              over water
!  pahfli   :             "              over ice
!
!  pevap    : surface evaporation (accumulated)
!  pevapl   :             "        over land
!  pevapw   :             "        over water
!  pevapi   :             "        over ice
!
!  paz0     : roughness length
!  paz0l    :      "            over land
!  paz0w    :      "            over water
!  paz0i    :      "            over ice
!
!  pustr    : u-stress (accumulated)
!  pustrl   :     "     over land
!  pustrw   :     "     over sea
!  pustri   :     "     over ice
!
!  pvstr    : v-stress (accumulated)
!  pvstrl   :     "     over land
!  pvstrw   :     "     over water
!  pvstri   :     "     over ice
!
!  pdew2    : dew point temperature at 2 meter
!  peforest : forest coverage
!  psn      : snow depth
!  psnc     : snow depth on canopy
!  ptemp2   : temperature at 2 meter
!  ptsm1    : surface temperature (t-dt)
!  pt2max   : maximum temp. at 2 m between output intervals
!  pt2min   : minimun temp. at 2 m between output intervals
!  pwind10w : 10m wind over water
!  pu10     : u-wind at 10 meter
!  pv10     : v-wind at 10 meter
!  pwind10  : wind speed at 10 meter (accumulated)
!  pwimax   : maximum windspeed at 10 m. between output intervals
!  pvdis    : boundary layer dissipation (accumulated)
!  pws      : surface soil wetness
!  pwsmx    : field capacity of soil
!  pwlmx    : skin reservoir
!  pvlt     : leaf area index
!  pvgrat   : vegetation ratio
!
! - 2d within physics only
!
!  paphm1   : half level pressure (t-dt)
!  papm1    : full level pressure (t-dt)
!  ptvm1    : virtual temperature at t-dt
!  pvdiffp  : rate of change of qv due to vdiff routine for cover
!  pvmixtau : vdiff mixing timescale for variance and skewness
!
! - 1d within physics only
!
!  pgeom1   : geopotential above surface (t-dt)
!  pqhfla   : moisture flux at the surface
!  pevapot  : potential evaporation
!  pcvs     : fractional snow cover (defined in *physc*)
!  pcvw     : wet skin fraction
!  ktropo   : tropopause index
!  lpland   : land-sea flag
!
!        Tendencies
!
!  - 3d
!
!  pxtte    : tendencies of tracer variables
!
!  - 2d
!  pvol     : tendency of meridional wind
!  pvom     : tendency of zonal wind
!  pqte     : tendency of humidity
!  ptte     : tendency of temperature
!  pxlte    : tendency of cloud water
!  pxite    : tendency of cloud ice
!
!
!     Method.
!     -------
!
!        First an auxialiary variable cp(q)t+gz is created on which
!   the vertical diffusion process will work like on u,v and q. then
!   along the vertical and at the surface, exchange coefficients (with
!   the dimension of a pressure thickness) are computed for momentum
!   and for heat (sensible plus latent). the letters m and h are used
!   to distinguish them. the diffusioncoefficents depend on the
!   turbulent kinetic energy (tke) calculated by an additional
!   prognostic equation, which considers advection of tke.
!        In the second part of the routine the implicit linear
!   systems for u,v first and t,q second are solved by a *gaussian
!   elimination back-substitution method. for t and q the lower
!   boundary condition depends on the surface state.
!   for tke the lower boundary condition depends on the square of
!   the frictional velocity.
!   over land, two different regimes of evaporation prevail:
!   a stomatal resistance dependent one over the vegetated part
!   and a soil relative humidity dependent one over the
!   bare soil part of the grid mesh.
!   potential evaporation takes place over the sea, the snow
!   covered part and the liquid water covered part of the
!   grid mesh as well as in case of dew deposition.
!        Finally one returns to the variable temperature to compute
!   its tendency and the later is modified by the dissipation's effect
!   (one assumes no storage in the turbulent kinetic energy range) and
!   the effect of moisture diffusion on cp. z0 is updated and the
!   surface fluxes of t and q and their derivatives are prepared and
!   stored like the difference between the implicitely obtained
!   cp(q)t+gz and cp(q)t at the surface.
!
!
!     Reference.
!
!          See vertical diffusion's part of the model's documentation
!     for details about the mathematics of this routine.
!
!     Authors.
!
!     u. schlese     dkrz-hamburg  feb-93
!       modified     e. roeckner  - 1994
!
!     j.-p. schulz   mpi - 1997 : implementation of implicit
!                                 coupling between land surface
!                                 and atmosphere.
!     m. esch, mpi, june 1999, echam5-modifications
!
!
!     based  on  original ecmwf version by j.f. geleyn  - 1982
!                              modified by c.b. blondin - 1986
!                                          h. feichter  - 1991
!                                          s. brinkop   - 1992
!                                          m. claussen  - 1993
USE mo_kind,             ONLY: dp
USE mo_geoloc,           ONLY: coriol_2d, amu0_x
USE mo_param_switches,   ONLY: lvdiff
USE mo_control,          ONLY: vct, nvclev
USE mo_physc2,           ONLY: clam, ckap, cb, cc, cchar, cvdifts, cfreec, cgam
USE mo_constants,        ONLY: vtmpc1, cpd, rd, g, vtmpc2, tmelt, alv, als
USE mo_vegetation,       ONLY: cva, cvb, cvc, cvbc, cvk, cvkc, cvabc,  &
                               cvrad
USE mo_tracer,           ONLY: trlist
USE mo_convect_tables,   ONLY: tlucua, jptlucu1, jptlucu2, &
                               lookuperror, lookupoverflow
USE mo_radiation,        ONLY: cemiss, co2mmr, ico2
USE mo_greenhouse_gases, ONLY: ghg_co2mmr
USE mo_exception,        ONLY: finish, message
USE mo_time_control,     ONLY: delta_time, lstart, time_step_len
USE mo_semi_impl,        ONLY: eps
USE mo_physc1,           ONLY: crae
USE mo_memory_g3b,       ONLY: tkevn, cdum, udif, vdif, ustarm
USE mo_surface
USE mo_hyb,              ONLY: apsurf
USE mo_co2,              ONLY: ico2idx  ! Index of CO2 in tracer list

!---wiso-code
USE mo_wiso,             ONLY: talphal1, talphal2, talphal3, talphal4, &
                               talphas1, talphas2, talphas3, talphas4, &
                               tkinsl,   tkinfa1,  tkinfa2,  tkinfa3,  &
                               toce,     tdifrel,                      &
                               tnat,     tsatbase, tsatfac,            &
                               lwisofracl, lwisokinl,                  &
                               cwisomin, cwisosec
!---wiso-code-end
!
IMPLICIT NONE

  integer, parameter :: xp = selected_real_kind(30,200)

!
! Variables declared
!
  INTEGER :: it, itop, itopp1, jk, jl, jt, klev, klevm1, klevp1
  INTEGER :: ibl, iblm1, krow, iblmax, iblmin, ktop
  INTEGER :: kproma, kbdim, ktdia, ktrac, xjrow
  REAL(dp)    :: z1dgam, z2geomf, zalf, zalh2
  REAL(dp)    :: zb, zbet, zblend, zbuoy, zc
  REAL(dp)    :: zchar, zchneu, zcons11, zcons12
  REAL(dp)    :: zcons13, zcons14, zcons15, zcons16, zcons17, zcons18
  REAL(dp)    :: zcons2, zcons23, zcons25, zcons29, zcons3, zcons30, zcons5
  REAL(dp)    :: zcons6, zcons8, zcons9, zcor, zcpd, zda1
  REAL(dp)    :: zdisx, zdisxt, zdivv, zdivv1, zdisc, zdisl
  REAL(dp)    :: zdqdt, zdqtot, zds, zdtdt, zdudt, zdus1
  REAL(dp)    :: zdus2, zdvdt, zdximdt, zdxlmdt, zdz
  REAL(dp)    :: zepcor, zepdu2, zepevap, zephum, zeps, zepsec, zepshr
  REAL(dp)    :: zepsr, zepz0o, zepzzo, zes, zfac, zfox
  REAL(dp)    :: zfreec, zfux, zgam, zh1, zh2, zhexp, zhtq
  REAL(dp)    :: zhuv, zkap, zkappa, zktest, zlam, zm1, zm2
  REAL(dp)    :: zm4, zmix, zmult1, zmult2, zmult3, zmult4
  REAL(dp)    :: zmult5, zplmax, zplmin, zqddif, zqdp
  REAL(dp)    :: zqtmit, zrd, zrdrv, zri, zrvrd, zsdep1
  REAL(dp)    :: zsdep2, zsh, zshear, zshn, zsm, zsmn
  REAL(dp)    :: zteldif, ztest, ztkemin, ztkesq
  REAL(dp)    :: ztmst, ztpfac1, ztpfac2, ztpfac3, ztpfac4
  REAL(dp)    :: zucf, zustf, zusus1, zva, zvabc, zvb, zvbc, zvc
  REAL(dp)    :: zvk, zvkc, zvrad, zwstf, zz2geo, zzb, zztvm, zzzlam
  REAL(dp)    :: zwslim, zdtime, zrhodz,  ztvh,  zdisv, zfav

!---wiso-code

  LOGICAL  :: lwiso
  INTEGER  :: kwiso

  REAL(dp)    :: zwisodqdt
  REAL(dp)    :: zwisodximdt, zwisodxlmdt

!---wiso-code-end

!     ------------------------------------------------------------------
!
  REAL(dp) ::                                                             &
       pxtm1(kbdim,klev,ktrac)
!
  REAL(dp) ::                                                             &
       pqm1(kbdim,klev),     ptm1(kbdim,klev),     pum1(kbdim,klev)       &
      ,pvm1(kbdim,klev),     pxlm1(kbdim,klev),    pxim1(kbdim,klev)      &
      ,pxvar(kbdim,klev)
  REAL(dp) ::                                                             &
       pahfl(kbdim),         pahfs(kbdim),         paz0(kbdim)            &
      ,pdew2(kbdim),         pevap(kbdim),         pforest(kbdim)         &
      ,ptemp2(kbdim),        pt2max(kbdim)                                &
      ,pt2min(kbdim),        pwind10w(kbdim),      pvdis(kbdim)           &
      ,pu10(kbdim),          pv10(kbdim),          pustr(kbdim)           &
      ,pvstr(kbdim),         pwimax(kbdim),        pwind10(kbdim)         &
      ,pwsmx(kbdim),         pvlt(kbdim)                                  &
      ,pvgrat(kbdim)                                                      &
      ,ptsl(kbdim),          ptsw(kbdim),          ptsi(kbdim)            &
      ,pocu(kbdim),          pocv(kbdim)                                  &
      ,paz0l(kbdim),         paz0w(kbdim),         paz0i(kbdim)           &
      ,pahfsl(kbdim),        pahfsw(kbdim),        pahfsi(kbdim)          &
      ,pahflw(kbdim)                                                      &
      ,pevapw(kbdim),        pevapi(kbdim)                                &
      ,pahfslac(kbdim),      pahfswac(kbdim),      pahfsiac(kbdim)        &
      ,pahfllac(kbdim),      pahflwac(kbdim),      pahfliac(kbdim)        &
      ,pevaplac(kbdim),      pevapwac(kbdim),      pevapiac(kbdim)        &
      ,pustrl(kbdim),        pustrw(kbdim),        pustri(kbdim)          &
      ,pvstrl(kbdim),        pvstrw(kbdim),        pvstri(kbdim)          &
      ,psn(kbdim),           psnc(kbdim),          ptslm1(kbdim)          &
      ,pws(kbdim),           palbedo(kbdim),       palbedo_vis(kbdim)     &
      ,palbedo_nir(kbdim),   palsol(kbdim)
  REAL(dp) ::                                                             &
       ptke(kbdim,klev),     ptkem1(kbdim,klev),   ptkem(kbdim,klev)      &
      ,paclc(kbdim,klev),    pemter(kbdim,klevp1)
  REAL(dp) ::                                                             &
       paphm1(kbdim,klevp1), papm1(kbdim,klev),    pgeom1(kbdim,klev)     &
      ,ptvm1(kbdim,klev),    pvdiffp(kbdim,klev),  pvmixtau(kbdim,klev)
  REAL(dp) ::                                                             &
       pcvs(kbdim),          pcvw(kbdim)                                  &
      ,pqhfla(kbdim),        pevapot(kbdim)                               &
      ,ptslnew(kbdim),       pwlmx(kbdim)
  LOGICAL ::                lpland(kbdim),        lpglac(kbdim)
  REAL(dp) ::                                                             &
       pfrl(kbdim),          pfrw(kbdim),          pfri(kbdim)
!
  REAL(dp) ::                                                             &
       pxtte(kbdim,klev,ktrac)
!
  REAL(dp) ::                                                          &
       pvol(kbdim,klev),     pvom(kbdim,klev),     pqte(kbdim,klev)    &
      ,ptte(kbdim,klev),     pxlte(kbdim,klev),    pxite(kbdim,klev)   &
      ,psiced(kbdim)                                                   &
      ,praind(kbdim),      psnowd(kbdim)                         & ! rain and snow total [kg/m2*s] over last time step
      ,jsswnir(kbdim),         jsswdifnir(kbdim)                       &
      ,jsswvis(kbdim),         jsswdifvis(kbdim)                       &
      ,pi0(kbdim)           , ptrsol(kbdim,klevp1)                     &
      , ptrflw(kbdim),          ptrfli(kbdim)    &
      , psofll(kbdim),         psoflw(kbdim),          psofli(kbdim)    &
      , ptrfllac(kbdim),       ptrflwac(kbdim),        ptrfliac(kbdim)  &
      , psofllac(kbdim),       psoflwac(kbdim),        psofliac(kbdim)  &
      , palsoi(kbdim),         palsow(kbdim),          pseaice(kbdim)   &
      , pseaice_new(kbdim),    psiced_new(kbdim)
      

  REAL(dp)  ::                                                          &
       psni(kbdim),           pahfice(kbdim),          pfluxres(kbdim), &
       pqres(kbdim),          pahfcon(kbdim),          pahfres(kbdim),  &
       pzti(kbdim),                                                     &
       pzteffl4(kbdim),       pztsnew(kbdim),          ptsurf(kbdim),   &
       prunoff(kbdim),        pdrain(kbdim),           ptsw_new(kbdim), &
       pradtemp_old(kbdim)
  REAL(dp)  :: palac(kbdim), pros_hd(kbdim), pdrain_hd(kbdim)
  REAL(dp)  :: pco2fluxo(kbdim), pco2fluxl(kbdim), pco2flux(kbdim), pco2flux_corr(kbdim)

!---wiso-code
  REAL(dp), OPTIONAL ::                                                 &
       pwisoqm1(kbdim,klev,kwiso), pwisoxlm1(kbdim,klev,kwiso),         &
       pwisoxim1(kbdim,klev,kwiso),                                     &
       pwisoevap(kbdim,kwiso),                                          &
       pwisoevapw(kbdim,kwiso), pwisoevapi(kbdim,kwiso),                &
       pwisoevaplac(kbdim,kwiso), pwisoevapwac(kbdim,kwiso),            &
       pwisoevapiac(kbdim,kwiso),                                       &
       pwisosn(kbdim,kwiso), pwisows(kbdim,kwiso),                      &
       pwisosw_d(kbdim,kwiso),                                          &
       pwisoqhfla(kbdim,kwiso), pwisoevapot(kbdim,kwiso),               &
       pwisoqte(kbdim,klev,kwiso), pwisoxlte(kbdim,klev,kwiso),         &
       pwisoxite(kbdim,klev,kwiso),                                     &
       pwisoraind(kbdim,kwiso), pwisosnowd(kbdim,kwiso),                &
       pwisorunoff(kbdim,kwiso), pwisodrain(kbdim,kwiso),               &
       pwisoalac(kbdim,kwiso), pwisoros_hd(kbdim,kwiso), pwisodrain_hd(kbdim,kwiso)
       
!---wiso-code-end

!!$      ,pswnet(kbdim), pswdown(kbdim)

!
!     local variables
!
  LOGICAL  :: lo
!
  REAL(dp) :: zxtdif(kbdim,klev,ktrac)
  REAL(dp) :: zdxtdt(kbdim,klev), zdxtdtsum, zdxtdtmean
  REAL(dp) :: zxtems(kbdim,ktrac)
  REAL(dp) :: zcfm(kbdim,klev),    zdis(kbdim,klev)                       &
         ,zcfh(kbdim,klev),    zcptgz(kbdim,klev),   zebsm(kbdim,klev)    &
         ,zudif(kbdim,klev),   zvdif(kbdim,klev)                          &
         ,ztcoe(kbdim)
  REAL(dp) :: ztdif(kbdim,klev)                                           &
         ,zqdif(kbdim,klev),   zebsh(kbdim,klev),    zvidis(kbdim)        &
         ,z1mxtm1(kbdim)
  REAL(dp) :: zhdyn(kbdim),        zteta1(kbdim,klev),   zlteta1(kbdim,klev)  &
         ,ztvir1(kbdim,klev),  zhh(kbdim,klevm1),    zqss(kbdim,klev)     &
         ,zxldif(kbdim,klev),  zxidif(kbdim,klev),   zedif(kbdim,klev)    &
         ,ztkevn(kbdim,klev),  zx(kbdim,klev)                             &
         ,zvardif(kbdim,klev)
  REAL(dp) :: zqssm(kbdim,klevm1),ztmitte(kbdim,klevm1),ztvirmit(kbdim,klevm1)&
         ,zfaxen(kbdim,klevm1),zfaxe(kbdim,klev),    zccover(kbdim,klevm1)&
         ,zlwcmit(kbdim,klevm1),ztemit(kbdim,klevm1),zqmit(kbdim,klevm1)  &
         ,zcdum(kbdim,klev),   zcfv(kbdim,klev),     zebsv(kbdim,klev)
!
  REAL(dp) :: zcrae, zmu0(kbdim)

  REAL(dp) ::  zqshear(kbdim,klev),zrho(kbdim,klevp1),zqflux(kbdim,klevp1)    &
          ,zvarpr(kbdim,klevp1)
!
  INTEGER  :: ihpblc(kbdim),   ihpbld(kbdim), ihpbl(kbdim)
  REAL(dp) :: zghabl(kbdim)
  REAL(dp) :: zph(klevp1), zp(klev), zh(klev)

!---wiso-code
  REAL(dp) :: zwisoqdif(kbdim,klev,kwiso)
  REAL(dp) :: zwisoxldif(kbdim,klev,kwiso), zwisoxidif(kbdim,klev,kwiso)
  REAL(dp) :: zwisokinw(kbdim,kwiso),  &                ! kinetic fraction factor at air/sea/(land) interface over ocean
              zwisokinl(kbdim,kwiso),  &                ! kinetic fraction factor at air/sea/(land) interface over land
              zwspeed, zwspeedmin,     &                ! wind speed and minimum wind speed for kin. fraction
              zwisosurfmin, zwisoroundoff
!---wiso-code-end

!
!     ------------------------------------------------------------------
!
!     THE FOLLOWING VARIABLES ARE NEEDED FOR THE SUBROUTINES SURFTEMP
!

  REAL(dp) ::    ztdif_new(kbdim), zqdif_new(kbdim)
  REAL(dp) ::    zqsurf_new(kbdim), ztvh_new(kbdim), zth_new(kbdim)
  REAL(dp) ::    ztte_corr(kbdim)

  REAL(dp) ::    zco2(kbdim)     ! CO2 concentration [kg/kg] in lowest level

!---wiso-code
  REAL(dp) ::    zwisoqdif_new(kbdim,kwiso)
!---wiso-code-end

!  Executable statements 

  lookupoverflow = .FALSE.

!
!     ------------------------------------------------------------------
!
!*    PHYSICAL CONSTANTS.
!
!          *ZLAM* IS THE ASYMPTOTIC MIXING LENGTH FOR MOMENTUM EXCHANGE,
!     *ZKAP* IS THE VON KARMAN CONSTANT, *ZB*, *ZC* AND *ZD* ARE SOME
!     CONSTANTS FOR THE FORMULAE ABOUT STABILITY DEPENDENCY RESPECTIVELY
!     NEAR THE NEUTRAL CASE, IN THE UNSTABLE CASE AND IN THE STABLE
!     CASE AND *ZCHAR* IS THE CONSTANT OF THE *CHARNOCK FORMULA.
!     *ZQWSSAT* AND *ZQSNCR* ARE THE INVERSES OF CRITICAL VALUES FOR
!     SOIL WATER AND SNOW DEPTH THAT ARE USED IN THE COMPUTATION OF THE
!     EVAPOTRANSPIRATION'S EFFICIENCY.
!
  xjrow=krow
  zlam=clam
  zkap=ckap
  zb=cb
  zc=cc
  zchar=cchar
  zva=cva
  zvb=cvb
  zvc=cvc
  zvbc=cvbc
  zvk=cvk
  zvkc=cvkc
  zvabc=cvabc
  zvrad=cvrad
  zrvrd=vtmpc1+1._dp
  zrdrv=1._dp/zrvrd
  zcpd=cpd
  zrd=rd
  zkappa=zrd/zcpd
!
!*      PARAMETERS FOR BOUNDARY LAYER DIAGNOSTICS
!
  zhuv=10._dp*g
  zhtq=2._dp*g
  zephum=5.e-2_dp
!
!*    SECURITY PARAMETERS.
!
  zepdu2=1.0_dp
  zepshr=1.e-5_dp
  zepzzo=1.5e-05_dp
  zepz0o=2._dp
  zepcor=5.e-05_dp
  ztkemin=1.e-10_dp
  zepsr=1.e-10_dp
  zepevap=1.e-10_dp
  zepsec=1.e-2_dp
!
!*    COMPUTATIONAL CONSTANTS.
!
  zdtime = delta_time
  ztmst  = time_step_len
  ztpfac1=cvdifts
  ztpfac2=1._dp/ztpfac1
  ztpfac3=1._dp-ztpfac2
  ztpfac4=1._dp+ztpfac3
  zzzlam=1._dp
  zcons2=0.5_dp*zkap/g
  zcons3=zlam
  zcons5=3._dp*zb*zc*g**2
  zcons6=1._dp/3._dp
  zcons8=2._dp*zb
  zcons9=3._dp*zb
  zcons11=3._dp*zb*zc
  zcons12=ztpfac1*ztmst*g/rd
  zcons13=1._dp/ztmst
  zcons14=zchar*rd/(g**2*ztmst)
  zcons15=1._dp/(g*ztmst)
  zcons16=cpd*vtmpc2
  zcons18=ztpfac1*ztmst*g**2
  zcons17=1._dp/zkap**2
  zcons25=zcons2/zcons3
  zcons29=ztpfac1*ztmst
  zcons30=1._dp/(ztpfac1*ztmst*g)
  zplmax=0.75_dp
  zplmin=0.35_dp
  zwslim=zplmin
  zblend=100._dp
  zchneu=.3_dp
  zfreec=cfreec
  zgam=cgam
  z1dgam=1._dp/zgam
  zh1= 2.22_dp
  zh2= 0.22_dp
  zm1= 1.24_dp
  zm2= 2.37_dp
  zm4= 3.69_dp
  zshn=zh1*zh2*SQRT(2._dp)
  zsmn=zshn*zm1*zm2/zm4
  zda1=1._dp/zsmn**3
  zustf=1._dp/zsmn**2
  zwstf=0.2_dp
!
  itop=1
  itopp1=itop+1
      
!---wiso-code
  IF (lwiso) THEN

! Computational Constants

! Minimum Windspeed for calculation of kinetic fractionation effects

  zwspeedmin=7._dp

! Security Paramters for calculation of delta values

  zwisosurfmin=1.e-12_dp
  zwisoroundoff=1.e-8_dp
  
  ENDIF
!---wiso-code-end
  
  ibl = klev-1            ! CO2 flux is distributed into lowest layers klev to ibl
  iblm1=ibl-1

!
! ------------------------------------------------------------
!

  ! CO2 mixing in PBL.
  ! Compute lowest (klev-2, approx. 500m in L19/31) and highest (highest below 2km) level ktop.
  ! CO2 is ideally mixed between bottom and ktop. 
  iblmin = klev-2

!-- half level pressure values, assuming 101320. Pa surface pressure

  DO jk=1,klevp1
    zph(jk)=vct(jk)+vct(jk+nvclev)*101320.0_dp
  END DO
!
! -- full level pressure
!
  DO jk = 1, klev
    zp(jk)=(zph(jk)+zph(jk+1))*0.5_dp
  END DO
!
  DO jk = 1, klev
    zh(jk)=(zph(klevp1)-zp(jk))/(g*1.25_dp)
  END DO
!
! -- search for highest level below 2000m
!
  DO jk = 1, klev
    iblmax=jk
    IF(zh(jk).LT.2000.0_dp) EXIT
  END DO
!
!     ------------------------------------------------------------------
!
!*         2.     NEW THERMODYNAMIC VARIABLE AND BOUNDARY CONDITIONS.
!
!*         2.1     REPLACE T BY CP(Q)*T+GZ IN THE ATMOSPHERE.
!
210 CONTINUE
  DO 212 jk=ktdia,klev
     DO 211 jl=1,kproma
        zx(jl,jk)=pxlm1(jl,jk)+pxim1(jl,jk)   ! total cloud water
        zcptgz(jl,jk)=pgeom1(jl,jk)+ptm1(jl,jk)                        &
                                   *cpd*(1._dp+vtmpc2*pqm1(jl,jk))
        zteta1(jl,jk)=ptm1(jl,jk)*(100000._dp/papm1(jl,jk))**zkappa
        ztvir1(jl,jk)=zteta1(jl,jk)*(1._dp+vtmpc1*pqm1(jl,jk)-zx(jl,jk))
        lo=ptm1(jl,jk).GE.tmelt
        zfaxe(jl,jk)=MERGE(alv,als,lo)
        zbet=zfaxe(jl,jk)/zcpd
        zusus1=zbet*zteta1(jl,jk)/ptm1(jl,jk)*zx(jl,jk)
        zlteta1(jl,jk)=zteta1(jl,jk)-zusus1
        it = NINT(ptm1(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/papm1(jl,jk)
        zes=MIN(zes,0.5_dp)
        zqss(jl,jk)=zes/(1._dp-vtmpc1*zes)
211  END DO
212 END DO

  IF (lookupoverflow) CALL lookuperror ('vdiff (1)   ')

  DO 214 jk=ktdia,klevm1
     DO 213 jl=1,kproma
        zhh(jl,jk)=(pgeom1(jl,jk)-pgeom1(jl,jk+1))
        zsdep1=(paphm1(jl,jk)-paphm1(jl,jk+1))                         &
              /(paphm1(jl,jk)-paphm1(jl,jk+2))
        zsdep2=(paphm1(jl,jk+1)-paphm1(jl,jk+2))                       &
              /(paphm1(jl,jk)  -paphm1(jl,jk+2))
        zqssm(jl,jk)=zsdep1*zqss(jl,jk)+zsdep2*zqss(jl,jk+1)
        ztmitte(jl,jk)=zsdep1*ptm1(jl,jk)+zsdep2*ptm1(jl,jk+1)
        ztvirmit(jl,jk)=zsdep1*ztvir1(jl,jk)+zsdep2*ztvir1(jl,jk+1)
        zfaxen(jl,jk)=zsdep1*zfaxe(jl,jk)+zsdep2*zfaxe(jl,jk+1)
        zlwcmit(jl,jk)=zsdep1*zx(jl,jk)+zsdep2*zx(jl,jk+1)
        zqmit(jl,jk)=zsdep1*pqm1(jl,jk)+zsdep2*pqm1(jl,jk+1)
        ztemit(jl,jk)=zsdep1*zteta1(jl,jk)+zsdep2*zteta1(jl,jk+1)
        zccover(jl,jk)=paclc(jl,jk)*zsdep1+paclc(jl,jk+1)*zsdep2
213  END DO
214 END DO

  IF (lvdiff) THEN

!
!!------------------- MO_SURFACE Coupling-----------------
! INITIALISE udif and vdif which previously came from momentum
! transfer und is then used for wind stress calculations
! wind stress is now calculated before momentum transfer
!- kalle 010904
!-----------------------------------------------------------
     IF (lstart) THEN

        DO jl=1,kproma
           udif(jl,xjrow)=ztpfac2*pum1(jl,klev)
           vdif(jl,xjrow)=ztpfac2*pvm1(jl,klev)
           ustarm(jl,xjrow)=1._dp
        END DO
        
     END IF
! JSBACH-end udif, vdif--------------------------------------

!
! Compute planetary boundary layer extension
!
! JSBACH note: in standard ECHAM5, ustarm is computed before zhdyn for the current time step.
!     Here, ustarm comes from the call to mo_surface at the previous timestep. 
!     But this should have only a minor effect on ihpbl and ghabl.
!
    DO jl = 1,kproma
       zcor=MAX(ABS(coriol_2d(jl,krow)),zepcor)
       zhdyn(jl)=MIN(pgeom1(jl,1)/g,zchneu*ustarm(jl,xjrow)/zcor)       
       ihpblc(jl)=klev
       ihpbld(jl)=klev
    END DO

     DO jk=klevm1,1,-1
        DO jl=1,kproma
           zds=zcptgz(jl,jk)-zcptgz(jl,klev)
           zdz=pgeom1(jl,jk)/g-zhdyn(jl)
           ihpblc(jl)=MERGE(jk,ihpblc(jl),                             &
                            ihpblc(jl).EQ.klev.AND.zds.GT.0._dp)
           ihpbld(jl)=MERGE(jk,ihpbld(jl),                             &
                            ihpbld(jl).EQ.klev.AND.zdz.GE.0._dp)
     END DO
  END DO
!
     DO  jl=1,kproma
        ihpbl(jl)=MIN(ihpblc(jl),ihpbld(jl))
        zghabl(jl)=MIN(50000._dp,pgeom1(jl,ihpbl(jl)))
     END DO
!

!
!     ==================================================================
!
!*       3.5   Vertical loop: Computation of basic quantities:
!              wind shear, buoyancy, Ri-number, mixing length
!
     DO 372 jk=ktdia,klevm1
        DO 361 jl=1,kproma
           zqtmit=zlwcmit(jl,jk)+zqmit(jl,jk)
           zfux=zfaxen(jl,jk)/(zcpd*ztmitte(jl,jk))
           zfox=zfaxen(jl,jk)/(zrd*ztmitte(jl,jk))
           zmult1=1._dp+vtmpc1*zqtmit
           zmult2=zfux*zmult1-zrvrd
           zmult3=zrdrv*zfox*zqssm(jl,jk)                              &
                  /(1._dp+zrdrv*zfux*zfox*zqssm(jl,jk))
           zmult5=zmult1-zmult2*zmult3
           zmult4=zfux*zmult5-1._dp
           zdus1=zccover(jl,jk)*zmult5+(1._dp-zccover(jl,jk))*zmult1
           zdus2=zccover(jl,jk)*zmult4+(1._dp-zccover(jl,jk))*vtmpc1
           zteldif=(zlteta1(jl,jk)-zlteta1(jl,jk+1))/zhh(jl,jk)*g
           zdqtot=(pqm1(jl,jk)+zx(jl,jk))-(pqm1(jl,jk+1)+zx(jl,jk+1))
           zqddif=zdqtot/zhh(jl,jk)*g
           zqshear(jl,jk)=zqddif !store for variance production
           zbuoy=(zteldif*zdus1+ztemit(jl,jk)*zdus2*zqddif)            &
                 *g/ztvirmit(jl,jk)
           zdivv=(pum1(jl,jk)-pum1(jl,jk+1))**2
           zdivv1=(pvm1(jl,jk)-pvm1(jl,jk+1))**2
           zshear=(zdivv+zdivv1)*(g/zhh(jl,jk))**2
           zri=zbuoy/MAX(zshear,zepshr)
!
!      ASYMPTOTIC MIXING LENGTH FOR MOMENTUM AND
!      HEAT (ZLAM) ABOVE THE PBL AS A FUNCTION OF HEIGHT
!      ACCORDING TO HOLTSLAG AND BOVILLE (1992), J. CLIMATE.
!
           zhexp=EXP(1._dp-pgeom1(jl,jk)/pgeom1(jl,ihpbl(jl)))
           zlam=zzzlam+(zcons3-zzzlam)*zhexp
           IF(jk.GE.ihpbl(jl)) THEN
              zcons23=zcons25
           ELSE
              zcons23=zcons2/zlam
           END IF
!
!     MIXING LENGTH (BLACKADAR) + STABILITY DEPENDENT FUNCTION
!
           z2geomf=pgeom1(jl,jk)+pgeom1(jl,jk+1)
           zz2geo=zcons2*z2geomf
           zmix=zz2geo/(1._dp+zcons23*z2geomf)
!
!      STABILITY FUNCTIONS (LOUIS, 1979)
!
           IF(zri.LT.0._dp) THEN
              zalh2=zmix*zmix
              zucf=1._dp/                                                 &
                   (1._dp+zcons5*zalh2*SQRT(ABS(zri)*(((pgeom1(jl,jk)     &
                       /pgeom1(jl,jk+1))**zcons6-1._dp)/(pgeom1(jl,jk)    &
                       -pgeom1(jl,jk+1)))**3/(pgeom1(jl,jk+1))))
              zsh=zshn*(1._dp-zcons9*zri*zucf)*zmix
              zsm=zsmn*(1._dp-zcons8*zri*zucf)*zmix
           ELSE
              zsh=zshn/(1._dp+zcons8*zri*SQRT(1._dp+zri))*zmix
              zsm=zsmn/(1._dp+zcons8*zri/SQRT(1._dp+zri))*zmix
           END IF
!
!       Dimensionless coefficients multiplied by pressure
!            thicknesses for momentum and heat exchange
!
           zzb=zshear*zsm-zbuoy*zsh
           zdisl=zda1*zmix/ztmst
           zktest=1._dp+(zzb*ztmst+SQRT(ptkem1(jl,jk))*2._dp)/zdisl
           IF (zktest.LE.1._dp) THEN
              ztkevn(jl,jk)=ztkemin
           ELSE
              ztkevn(jl,jk)=MAX(ztkemin,(zdisl*(SQRT(zktest)-1._dp))**2)
           END IF
           IF(lstart) THEN
              ptkem1(jl,jk)=ztkevn(jl,jk)
              ptkem(jl,jk)=ztkevn(jl,jk)
           END IF
           ztkesq=SQRT(MAX(ztkemin,ptkem1(jl,jk)))
           zztvm=(ptvm1(jl,jk)+ptvm1(jl,jk+1))*0.5_dp
           zalf=paphm1(jl,jk+1)/(zztvm*zhh(jl,jk)*zrd)
           zcfm(jl,jk)=zsm*ztkesq*zcons18*zalf
           zcfh(jl,jk)=zsh*ztkesq*zcons18*zalf
           zcfv(jl,jk)=0.5_dp*zcfh(jl,jk)
           zcdum(jl,jk)=zcfm(jl,jk)/ztkesq*SQRT(ztkevn(jl,jk))
361     END DO
372  END DO

!---wiso-code
  IF (lwiso) THEN

! calculation of the kinetic fractionation factor *zwisokinw* over open water
! kinetic fractionation according to formula of brutsaert, 1975 (see also: g. hoffmann, dissertation, p.21)
     DO jt=1,kwiso
        DO jl=1,kproma
!          absolute value of windspeed
           zwspeed=sqrt(pum1(jl,klev)**2+pvm1(jl,klev)**2)
           IF (zwspeed.le.zwspeedmin) THEN
              zwisokinw(jl,jt)=tkinsl(jt)
           ELSE
              zwisokinw(jl,jt)=tkinfa1(jt)*zwspeed+tkinfa2(jt)
           ENDIF
           zwisokinw(jl,jt)=1._dp-zwisokinw(jl,jt)
        END DO
     END DO

! calculation of the kinetic fractionation factor *zwisokinl* over land surfaces
  IF (lwisokinl .EQ. 0) THEN
! no kinetic fractionation over land surface
     DO jt=1,kwiso
        DO jl=1,kproma
           zwisokinl(jl,jt)=1._dp
        END DO
     END DO
  ENDIF

  IF (lwisokinl .EQ. 1) THEN
! kinetic fractionation according to formula of brutsaert, 1975 (see also: g. hoffmann, dissertation, p.21)
     DO jt=1,kwiso
        DO jl=1,kproma
!          absolute value of windspeed
           zwspeed=sqrt(pum1(jl,klev)**2+pvm1(jl,klev)**2)
!          kin. fractionation according to formula of brutsaert, 1975
!          (see also: g. hoffmann, dissertation, p.21)
           IF (zwspeed.le.zwspeedmin) THEN
              zwisokinl(jl,jt)=tkinsl(jt)
           ELSE
              zwisokinl(jl,jt)=tkinfa1(jt)*zwspeed+tkinfa2(jt)
           ENDIF
           zwisokinl(jl,jt)=1._dp-zwisokinl(jl,jt)
        END DO
     END DO
  ENDIF

  IF (lwisokinl .EQ. 2) THEN
! kinetic fractionation over land surface according to Mathieu & Bariac (1996)
! with n-power coefficient (n=0.67) as suggested by Riley et al. (2002)
     DO jt=1,kwiso
        DO jl=1,kproma
           zwisokinl(jl,jt)=1._dp/(tdifrel(jt)**0.67_dp)
        END DO
     END DO
  ENDIF
    
     
  END IF
!---wiso-code-end

!     ------------------------------------------------------------------
!
!*         5.     DIFFUSION IMPLICIT COMPUTATIONS FOR HEAT (S.+L.).
!
500  CONTINUE
     DO 502 jk=1,klev
        DO 501 jl=1,kproma
           ztdif(jl,jk)=0._dp
           zqdif(jl,jk)=0._dp
           zxldif(jl,jk)=0._dp
           zxidif(jl,jk)=0._dp
           zvardif(jl,jk)=0._dp
501     END DO
502  END DO

!---wiso-code
  IF (lwiso) THEN

    DO jt=1,kwiso
       DO jk=1,klev
          DO jl=1,kproma
             zwisoqdif (jl,jk,jt)=0._dp
             zwisoxldif(jl,jk,jt)=0._dp
             zwisoxidif(jl,jk,jt)=0._dp
          END DO
       END DO
    END DO

  END IF     
!---wiso-code-end

!
!
!*         5.1     SETTING OF RIGHT HAND SIDES.
!
510  CONTINUE
     DO 512 jk=itop,klev
        DO 511 jl=1,kproma
           ztdif(jl,jk)=ztpfac2*zcptgz(jl,jk)
           zqdif(jl,jk)=ztpfac2*pqm1(jl,jk)
           zxldif(jl,jk)=ztpfac2*pxlm1(jl,jk)
           zxidif(jl,jk)=ztpfac2*pxim1(jl,jk)
           zvardif(jl,jk)=ztpfac2*pxvar(jl,jk)
511     END DO
512  END DO

!---wiso-code
  IF (lwiso) THEN

     DO jt=1,kwiso
       DO jk=itop,klev
            DO jl=1,kproma
               zwisoqdif (jl,jk,jt)=ztpfac2*pwisoqm1 (jl,jk,jt)
               zwisoxldif(jl,jk,jt)=ztpfac2*pwisoxlm1(jl,jk,jt)
               zwisoxidif(jl,jk,jt)=ztpfac2*pwisoxim1(jl,jk,jt)
             END DO
        END DO
      END DO
     
  END IF
!---wiso-code-end

!
!
!*         5.2     TOP LAYER ELIMINATION.
!
520  CONTINUE
!
     DO 521 jl=1,kproma
        zqdp=1._dp/(paphm1(jl,itopp1)-paphm1(jl,itop))
        zdisc=1._dp/(1._dp+zcfh(jl,itop)*zqdp)
        zebsh(jl,itop)=zdisc*(zcfh(jl,itop)*zqdp)
        zdisv=1._dp/(1._dp+zcfv(jl,itop)*zqdp)
        zebsv(jl,itop)=zdisv*(zcfv(jl,itop)*zqdp)
        ztdif(jl,itop)=zdisc*ztdif(jl,itop)
        zqdif(jl,itop)=zdisc*zqdif(jl,itop)
        zxldif(jl,itop)=zdisc*zxldif(jl,itop)
        zxidif(jl,itop)=zdisc*zxidif(jl,itop)
        zvardif(jl,itop)=zdisv*zvardif(jl,itop)
521  END DO

!---wiso-code
  IF (lwiso) THEN

     DO jt=1,kwiso
        DO jl=1,kproma
           zqdp=1._dp/(paphm1(jl,itopp1)-paphm1(jl,itop))
           zdisc=1._dp/(1._dp+zcfh(jl,itop)*zqdp)
           zwisoqdif (jl,itop,jt)=zdisc*zwisoqdif (jl,itop,jt)
           zwisoxldif(jl,itop,jt)=zdisc*zwisoxldif(jl,itop,jt)
           zwisoxidif(jl,itop,jt)=zdisc*zwisoxidif(jl,itop,jt)
        END DO
     END DO
     
  END IF
!---wiso-code-end

!
!
!*         5.3     ELIMINATION FOR MIDDLE LAYERS.
!
530  CONTINUE
!
     DO 532 jk=itopp1,klevm1
        DO 531 jl=1,kproma
           zqdp=1._dp/(paphm1(jl,jk+1)-paphm1(jl,jk))
           zfac=zcfh(jl,jk-1)*zqdp
           zdisc=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,jk-1))+zcfh(jl,jk)*zqdp)
           zebsh(jl,jk)=zdisc*(zcfh(jl,jk)*zqdp)
           zfav=zcfv(jl,jk-1)*zqdp
           zdisv=1._dp/(1._dp+zfav*(1._dp-zebsv(jl,jk-1))+zcfv(jl,jk)*zqdp)
           zebsv(jl,jk)=zdisv*(zcfv(jl,jk)*zqdp)
           ztdif(jl,jk)=zdisc*(ztdif(jl,jk)+zfac*ztdif(jl,jk-1))
           zqdif(jl,jk)=zdisc*(zqdif(jl,jk)+zfac*zqdif(jl,jk-1))
           zxldif(jl,jk)=zdisc*(zxldif(jl,jk)+zfac*zxldif(jl,jk-1))
           zxidif(jl,jk)=zdisc*(zxidif(jl,jk)+zfac*zxidif(jl,jk-1))
           zvardif(jl,jk)=zdisv*(zvardif(jl,jk)+zfav*zvardif(jl,jk-1))
531     END DO
532  END DO
!
!---wiso-code
  IF (lwiso) THEN

    DO jt=1,kwiso
       DO jk=itopp1,klevm1
          DO jl=1,kproma
           zqdp=1._dp/(paphm1(jl,jk+1)-paphm1(jl,jk))
           zfac=zcfh(jl,jk-1)*zqdp
           zdisc=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,jk-1))              &
                                                    +zcfh(jl,jk)*zqdp)
           zwisoqdif (jl,jk,jt) = zdisc*(zwisoqdif (jl,jk,jt)+zfac*zwisoqdif (jl,jk-1,jt))
           zwisoxldif(jl,jk,jt) = zdisc*(zwisoxldif(jl,jk,jt)+zfac*zwisoxldif(jl,jk-1,jt))
           zwisoxidif(jl,jk,jt) = zdisc*(zwisoxidif(jl,jk,jt)+zfac*zwisoxidif(jl,jk-1,jt))
          END DO
       END DO
    END DO
    
  END IF
!---wiso-code-end

!
!
!*         5.4     BOTTOM LAYER ELIMINATION.

     DO 541 jl=1,kproma
        zqdp=1._dp/(paphm1(jl,klevp1)-paphm1(jl,klev))
        zfac=zcfh(jl,klevm1)*zqdp
        zdisx=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1)))
        zfav=zcfv(jl,klevm1)*zqdp
        zdisv=1._dp/(1._dp+zfav*(1._dp-zebsv(jl,klevm1)))
        zxldif(jl,klev)=zdisx*(zxldif(jl,klev)+zfac*zxldif(jl,klevm1))
        zxidif(jl,klev)=zdisx*(zxidif(jl,klev)+zfac*zxidif(jl,klevm1))
        zvardif(jl,klev)=zdisv*(zvardif(jl,klev)+zfav*zvardif(jl,klevm1))

541  END DO
!
!---wiso-code
  IF (lwiso) THEN

    DO jt=1,kwiso
      DO jl=1,kproma
        zqdp=1._dp/(paphm1(jl,klevp1)-paphm1(jl,klev))
        zfac=zcfh(jl,klevm1)*zqdp
        zdisx=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1)))

        zwisoxldif(jl,klev,jt)=zdisx*(zwisoxldif(jl,klev,jt)+zfac*zwisoxldif(jl,klevm1,jt))
        zwisoxidif(jl,klev,jt)=zdisx*(zwisoxidif(jl,klev,jt)+zfac*zwisoxidif(jl,klevm1,jt))
      END DO
    END DO
    
  END IF
!---wiso-code-end

!
! solar zenith angle
     zcrae = crae*(crae+2._dp)
     zmu0(1:kproma) = crae/(SQRT(amu0_x(1:kproma,krow)**2+zcrae)-amu0_x(1:kproma,krow))

     IF (ico2idx > 0) THEN
        IF (trlist%ti(ico2idx)%nvdiff /= 1) &
             CALL finish('vdiff','Need %nvdiff = 1 for CO2 tracer')
        zco2(1:kproma) = pxtm1(1:kproma,klev,ico2idx)
!!$     ELSE
!!$        CALL finish('vdiff','CO2 not found in tracer list')
     ELSE IF (ico2 == 2) THEN
        zco2(1:kproma) = co2mmr
     ELSE IF (ico2 == 4) THEN
        zco2(1:kproma) = ghg_co2mmr
     ELSE
        CALL finish('vdiff','co2: this "ico2" is not supported')
     END IF

  IF (lwiso) THEN
     CALL update_surface(                 & !! logistic parameters
       kproma,                            & !! length of vector
       krow,                              &
       klev,                              & !! lowest level number
       klevp1,                            & !! lln plus 1
       klevm1,                            & !! lln minus 1
!!---INPUT------------------!! Atmospheric conditions lowest level INPUT
       pxlm1(1:kproma,klev),              & !! liquid clouds (water)
       pxim1(1:kproma,klev),              & !! frozen clouds (ice)
       pgeom1(1:kproma,klev),             & !! geopotential above surface
       ptm1(1:kproma,klev),               & !! temperature (t-dt)
       pqm1(1:kproma,klev),               & !! humidity (t-dt)
       papm1(1:kproma,klev),              & !! full level pressure (means middle of the layer)
       paphm1(1:kproma,1:klevp1),         & !! half level pressures (bottom of the layer)
       pum1(1:kproma,klev),               & !! wind_u
       pvm1(1:kproma,klev),               & !! wind_v
       praind(1:kproma),                  & !! rain (convective and large scale) over last time step
       psnowd(1:kproma),                  & !! snow (convective and large scale) over last time step
       pemter(1:kproma,klevp1),           & !! longwave down
       jsswvis(1:kproma),                 & !! net surface visible
       jsswdifvis(1:kproma),              & !! fraction of diffuse visible
       jsswnir(1:kproma),                 & !! net surface near infrared
       jsswdifnir(1:kproma),              & !! fraction of diffuse near infrared
       zmu0(1:kproma),                    & !! solar zenith angle
       paclc(1:kproma,klev),              & !! cloud cover
       ptsw(1:kproma),                    & !! ocean_temp read from sst field or ocean model
       pocu(1:kproma),                    & !! ocean_u_velocity
       pocv(1:kproma),                    & !! ocean_v_velocity
       psiced(1:kproma),                  & !! seaice depth from clsst
       zcfh(1:kproma,1:klev),             & !! turb diffusion coeff for RM-scheme
       zebsh(1:kproma,1:klev),            & !! coeff from diffusion scheme for RM
       zqdif(1:kproma,1:klev),            & !! humid. coeff from diffusion scheme
       ztdif(1:kproma,1:klev),            & !! temp. coeff from diffusion scheme
       udif(1:kproma,xjrow),              & !! wind speed coeff from momentum diffusion
       vdif(1:kproma,xjrow),              & !! wind speed coeff ..for wind stress calc.
       zghabl(1:kproma),                  & !! Geopotential of PBL extension 
       pi0(1:kproma),                     & !! solar incidence factor from radiation scheme
       ptrsol(1:kproma,klevp1),           & !! solar radiation
       zco2(1:kproma),                    & !! CO2 concentration in lowest level
!! - OUTPUT FROM SURFACE SCHEME ----------------------------------------
       palbedo(1:kproma),       palbedo_vis(1:kproma),                 &
       palbedo_nir(1:kproma),                                          &
       ptrflw(1:kproma),        ptrfli(1:kproma),                      &
       psofll(1:kproma),        psoflw(1:kproma),                      &
       psofli(1:kproma) ,       ptrfllac(1:kproma),                    &
       ptrflwac(1:kproma),      ptrfliac(1:kproma),                    &
       psofllac(1:kproma),      psoflwac(1:kproma),                    &
       psofliac(1:kproma),      palsol(1:kproma),                      &
       palsoi(1:kproma),        palsow(1:kproma),                      &
       ustarm(1:kproma,xjrow),  cdum(1:kproma,xjrow),                  & 
       tkevn(1:kproma,xjrow),   ztdif_new(1:kproma),                   &
       zqdif_new(1:kproma),     ztvh_new(1:kproma),                    &
       zqsurf_new(1:kproma),    zth_new(1:kproma),                     &
       pwind10w(1:kproma) ,     pu10(1:kproma),                       &
       pv10(1:kproma) ,         pwimax(1:kproma),                      &
       pwind10(1:kproma),       pdew2(1:kproma),                       &
       ptemp2(1:kproma),        pt2max(1:kproma),                      &
       pt2min(1:kproma),        pevaplac(1:kproma),                    &
       pevapwac(1:kproma),      pevapiac(1:kproma),                    &
       pevap(1:kproma),         pahfllac(1:kproma),                    &
       pahflwac(1:kproma),      pahfliac(1:kproma),                    &
       pahfl(1:kproma),         pahfslac(1:kproma),                    &
       pahfswac(1:kproma),      pahfsiac(1:kproma),                    &
       pahfs(1:kproma),         pqhfla(1:kproma),                      &
       pevapw(1:kproma),                      &
       pevapi(1:kproma),        pahfsl(1:kproma),                      &
       pahfsw(1:kproma),        pahfsi(1:kproma),                      &
       pevapot(1:kproma),                                              &
       pahflw(1:kproma),                                               &
       psni(1:kproma),          pahfice(1:kproma),                     &  !Note: psni is INOUT (comes already from ocean)
       pfluxres(1:kproma),      pqres(1:kproma),                       &
       pahfcon(1:kproma),       pahfres(1:kproma),                     &
       ptsi(1:kproma),          ptslnew(1:kproma),                     &
       pzti(1:kproma),          pzteffl4(1:kproma),                    &
       pztsnew(1:kproma),       ptsurf(1:kproma),                      &
       paz0w(1:kproma),         paz0i(1:kproma),                       &
       paz0l(1:kproma),         paz0(1:kproma),                        &
       pustrl(1:kproma),        pvstrl(1:kproma),                      &
       pustrw(1:kproma),        pvstrw(1:kproma),                      &
       pustri(1:kproma),        pvstri(1:kproma),                      &
       prunoff(1:kproma),       pdrain(1:kproma),                      &
       pustr(1:kproma),         pvstr(1:kproma),                       &
       ztte_corr(1:kproma),     ptsw_new(1:kproma),                    &
       pseaice(1:kproma),       pseaice_new(1:kproma),                 &
       psiced_new(1:kproma),    pradtemp_old(1:kproma),                &
       palac(1:kproma), pros_hd(1:kproma), pdrain_hd(1:kproma),        &
       pco2fluxo(1:kproma), pco2fluxl(1:kproma),pco2flux(1:kproma),    &
!---wiso-code
       ! Input
       lwiso, kwiso, lwisofracl,                                       &
       pwisoqm1(1:kproma,klev,1:kwiso),                                &
       pwisoraind(1:kproma,1:kwiso), pwisosnowd(1:kproma,1:kwiso),     &
       zwisoqdif(1:kproma,1:klev,1:kwiso),                             &
       zwisokinw(1:kproma,1:kwiso), zwisokinl(1:kproma,1:kwiso),       &	
       pwisosw_d(1:kproma,1:kwiso),                                    &
       ! Output
       zwisoqdif_new(1:kproma,1:kwiso),                                &		
       pwisoevaplac(1:kproma,1:kwiso), pwisoevapwac(1:kproma,1:kwiso), &	
       pwisoevapiac(1:kproma,1:kwiso), pwisoevap(1:kproma,1:kwiso),    &		
       pwisoqhfla(1:kproma,1:kwiso),                                   &	
       pwisoevapw(1:kproma,1:kwiso),                                   &	
       pwisoevapi(1:kproma,1:kwiso), pwisoevapot(1:kproma,1:kwiso),    &
       pwisorunoff(1:kproma,1:kwiso), pwisodrain(1:kproma,1:kwiso),    &
       pwisoalac(1:kproma,1:kwiso), pwisoros_hd(1:kproma,1:kwiso),     &
       pwisodrain_hd(1:kproma,1:kwiso))
!---wiso-code-end
  ELSE
     CALL update_surface(                 & !! logistic parameters
       kproma,                            & !! length of vector
       krow,                              &
       klev,                              & !! lowest level number
       klevp1,                            & !! lln plus 1
       klevm1,                            & !! lln minus 1
!!---INPUT------------------!! Atmospheric conditions lowest level INPUT
       pxlm1(1:kproma,klev),              & !! liquid clouds (water)
       pxim1(1:kproma,klev),              & !! frozen clouds (ice)
       pgeom1(1:kproma,klev),             & !! geopotential above surface
       ptm1(1:kproma,klev),               & !! temperature (t-dt)
       pqm1(1:kproma,klev),               & !! humidity (t-dt)
       papm1(1:kproma,klev),              & !! full level pressure (means middle of the layer)
       paphm1(1:kproma,1:klevp1),         & !! half level pressures (bottom of the layer)
       pum1(1:kproma,klev),               & !! wind_u
       pvm1(1:kproma,klev),               & !! wind_v
       praind(1:kproma),                  & !! rain (convective and large scale) over last time step
       psnowd(1:kproma),                  & !! snow (convective and large scale) over last time step
       pemter(1:kproma,klevp1),           & !! longwave down
       jsswvis(1:kproma),                 & !! net surface visible
       jsswdifvis(1:kproma),              & !! fraction of diffuse visible
       jsswnir(1:kproma),                 & !! net surface near infrared
       jsswdifnir(1:kproma),              & !! fraction of diffuse near infrared
       zmu0(1:kproma),                    & !! solar zenith angle
       paclc(1:kproma,klev),              & !! cloud cover
       ptsw(1:kproma),                    & !! ocean_temp read from sst field or ocean model
       pocu(1:kproma),                    & !! ocean_u_velocity
       pocv(1:kproma),                    & !! ocean_v_velocity
       psiced(1:kproma),                  & !! seaice depth from clsst
       zcfh(1:kproma,1:klev),             & !! turb diffusion coeff for RM-scheme
       zebsh(1:kproma,1:klev),            & !! coeff from diffusion scheme for RM
       zqdif(1:kproma,1:klev),            & !! humid. coeff from diffusion scheme
       ztdif(1:kproma,1:klev),            & !! temp. coeff from diffusion scheme
       udif(1:kproma,xjrow),              & !! wind speed coeff from momentum diffusion
       vdif(1:kproma,xjrow),              & !! wind speed coeff ..for wind stress calc.
       zghabl(1:kproma),                  & !! Geopotential of PBL extension 
       pi0(1:kproma),                     & !! solar incidence factor from radiation scheme
       ptrsol(1:kproma,klevp1),           & !! solar radiation
       zco2(1:kproma),                    & !! CO2 concentration in lowest level
!! - OUTPUT FROM SURFACE SCHEME ----------------------------------------
       palbedo(1:kproma),       palbedo_vis(1:kproma),                 &
       palbedo_nir(1:kproma),                                          &
       ptrflw(1:kproma),        ptrfli(1:kproma),                      &
       psofll(1:kproma),        psoflw(1:kproma),                      &
       psofli(1:kproma) ,       ptrfllac(1:kproma),                    &
       ptrflwac(1:kproma),      ptrfliac(1:kproma),                    &
       psofllac(1:kproma),      psoflwac(1:kproma),                    &
       psofliac(1:kproma),      palsol(1:kproma),                      &
       palsoi(1:kproma),        palsow(1:kproma),                      &
       ustarm(1:kproma,xjrow),  cdum(1:kproma,xjrow),                  & 
       tkevn(1:kproma,xjrow),   ztdif_new(1:kproma),                   &
       zqdif_new(1:kproma),     ztvh_new(1:kproma),                    &
       zqsurf_new(1:kproma),    zth_new(1:kproma),                     &
       pwind10w(1:kproma) ,     pu10(1:kproma),                       &
       pv10(1:kproma) ,         pwimax(1:kproma),                      &
       pwind10(1:kproma),       pdew2(1:kproma),                       &
       ptemp2(1:kproma),        pt2max(1:kproma),                      &
       pt2min(1:kproma),        pevaplac(1:kproma),                    &
       pevapwac(1:kproma),      pevapiac(1:kproma),                    &
       pevap(1:kproma),         pahfllac(1:kproma),                    &
       pahflwac(1:kproma),      pahfliac(1:kproma),                    &
       pahfl(1:kproma),         pahfslac(1:kproma),                    &
       pahfswac(1:kproma),      pahfsiac(1:kproma),                    &
       pahfs(1:kproma),         pqhfla(1:kproma),                      &
       pevapw(1:kproma),                      &
       pevapi(1:kproma),        pahfsl(1:kproma),                      &
       pahfsw(1:kproma),        pahfsi(1:kproma),                      &
       pevapot(1:kproma),                                              &
       pahflw(1:kproma),                                               &
       psni(1:kproma),          pahfice(1:kproma),                     &  !Note: psni is INOUT (comes already from ocean)
       pfluxres(1:kproma),      pqres(1:kproma),                       &
       pahfcon(1:kproma),       pahfres(1:kproma),                     &
       ptsi(1:kproma),          ptslnew(1:kproma),                     &
       pzti(1:kproma),          pzteffl4(1:kproma),                    &
       pztsnew(1:kproma),       ptsurf(1:kproma),                      &
       paz0w(1:kproma),         paz0i(1:kproma),                       &
       paz0l(1:kproma),         paz0(1:kproma),                        &
       pustrl(1:kproma),        pvstrl(1:kproma),                      &
       pustrw(1:kproma),        pvstrw(1:kproma),                      &
       pustri(1:kproma),        pvstri(1:kproma),                      &
       prunoff(1:kproma),       pdrain(1:kproma),                      &
       pustr(1:kproma),         pvstr(1:kproma),                       &
       ztte_corr(1:kproma),     ptsw_new(1:kproma),                    &
       pseaice(1:kproma),       pseaice_new(1:kproma),                 &
       psiced_new(1:kproma),    pradtemp_old(1:kproma),                &
       palac(1:kproma), pros_hd(1:kproma), pdrain_hd(1:kproma),        &
       pco2fluxo(1:kproma), pco2fluxl(1:kproma),pco2flux(1:kproma),    &
!---wiso-code
       ! Input
       lwiso, kwiso)
!---wiso-code-end
  END IF
!
!    INITIALIZE SURFACE EMISSION FOR TRACERS
!
     IF(ktrac.GT.0) THEN

        DO jt=1,ktrac
           DO 3230 jl=1,kproma
              zxtems(jl,jt)=0._dp
3230       END DO
        END DO

        ! Copy CO2 flux (from ocean and from land) into emissions field
        IF (ico2idx > 0) zxtems(1:kproma,ico2idx) = pco2flux(1:kproma)
        !
        !     SURFACE EMISSIONS AND DRY DEPOSITION
        !
!!$        IF(lemis) THEN
!!$           DO 324 jl=1,kproma
!!$              z1mxtm1(jl)=papm1(jl,klev)                               &
!!$                   /(ptm1(jl,klev)*rd*(1._dp+vtmpc1*pqm1(jl,klev)))
!!$324        END DO
!!$!
!!$           CALL xtemiss (kproma,   klev,     irow,   cvdifts,  zdtime,   &
!!$                         pxtm1,  zxtems,   z1mxtm1,                    &
!!$                         lpland, pforest,  psn)
!!$        END IF

        ! Add emissions to CO2 flux (also computes total flux, including emissions but without flux correction)
        CALL call_chem_bcond(kproma, kbdim, klev, zxtems, pxtte, xjrow)

        ! Add CO2 flux correction
        IF (ico2idx > 0) zxtems(1:kproma,ico2idx) = zxtems(1:kproma,ico2idx) + pco2flux_corr(1:kproma)

     END IF

     DO 505 jt=1,trlist% ntrac
        IF (trlist% ti(jt)% nvdiff /= 1) CYCLE
        DO 504 jk=1,klev
           DO 503 jl=1,kproma
              zxtdif(jl,jk,jt)=0._dp
503        END DO
504     END DO

        DO 516 jk=itop,iblm1
!!$        DO 516 jk=itop,klev
           DO 514 jl=1,kproma
              zxtdif(jl,jk,jt)=ztpfac2*pxtm1(jl,jk,jt)
514        END DO
516     END DO

        DO jk=ibl,klev
           DO jl=1,kproma
              zqdp = 1._dp/(paphm1(jl,klevp1)-paphm1(jl,ibl))
              zxtdif(jl,jk,jt) = ztpfac2 * pxtm1(jl,jk,jt)   &
                       + ztmst * g * zqdp * zxtems(jl,jt)
           END DO
        END DO

        DO 526 jl=1,kproma
           zqdp=1._dp/(paphm1(jl,itopp1)-paphm1(jl,itop))
           zdisc=1._dp/(1._dp+zcfh(jl,itop)*zqdp)
           zxtdif(jl,itop,jt)=zdisc*zxtdif(jl,itop,jt)
526     END DO

        DO 536 jk=itopp1,klevm1
           DO 534 jl=1,kproma
              zqdp=1._dp/(paphm1(jl,jk+1)-paphm1(jl,jk))
              zfac=zcfh(jl,jk-1)*zqdp
              zdisc=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,jk-1))+zcfh(jl,jk)*zqdp)
              zxtdif(jl,jk,jt)=zdisc *        &
                              (zxtdif(jl,jk,jt) + zfac*zxtdif(jl,jk-1,jt))
534        END DO
536     END DO

        DO 543 jl=1,kproma
           zqdp=1._dp/(paphm1(jl,klevp1)-paphm1(jl,klev))
           zfac=zcfh(jl,klevm1)*zqdp
           zdisxt=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1)))
           zxtdif(jl,klev,jt)=zdisxt * (zxtdif(jl,klev,jt)     &
!!$                                        + ztmst * g * zqdp * zxtems(jl,jt) &
                                        + zfac*zxtdif(jl,klevm1,jt)        &
                                       )
543     END DO

505  END DO


!--------------------
! Hand it over ! begin
!--------------------

     DO jl=1,kproma
        zcdum(jl,klev)  = cdum(jl,xjrow)
        zcfm(jl,klev)   = cdum(jl,xjrow)
        ztkevn(jl,klev) = tkevn(jl,xjrow)
        ztdif(jl,klev)  = ztdif_new(jl)
        zqdif(jl,klev)  = zqdif_new(jl)
        ptsw(jl)        = ptsw_new(jl)
        pseaice(jl)     = pseaice_new(jl)
        psiced(jl)      = psiced_new(jl)
     END DO
     IF(lstart) THEN
        DO jl=1,kproma
           ptkem1(jl,klev)=ztkevn(jl,klev)
           ptkem(jl,klev)=ztkevn(jl,klev)
        END DO
     END IF

!---wiso-code
  IF (lwiso) THEN

     DO jt=1,kwiso
         DO jl=1,kproma
             zwisoqdif(jl,klev,jt)  = zwisoqdif_new(jl,jt)
         END DO
     END DO

  END IF
!---wiso-code-end 

!--------------------
! Hand it over ! end
!--------------------

!
!*         5.5     BACK-SUBSTITUTION.
!
550  CONTINUE
!
     DO 552 jk=klevm1,itop,-1
        DO 551 jl=1,kproma
           ztdif(jl,jk)=ztdif(jl,jk)+zebsh(jl,jk)*ztdif(jl,jk+1)
           zqdif(jl,jk)=zqdif(jl,jk)+zebsh(jl,jk)*zqdif(jl,jk+1)
           zxldif(jl,jk)=zxldif(jl,jk)+zebsh(jl,jk)*zxldif(jl,jk+1)
           zxidif(jl,jk)=zxidif(jl,jk)+zebsh(jl,jk)*zxidif(jl,jk+1)
           zvardif(jl,jk)=zvardif(jl,jk)+zebsv(jl,jk)*zvardif(jl,jk+1)
551     END DO
552  END DO

!---wiso-code
  IF (lwiso) THEN

     DO jt=1,kwiso
        DO jk=klevm1,itop,-1
           DO jl=1,kproma
              zwisoqdif (jl,jk,jt)=zwisoqdif (jl,jk,jt)+zebsh(jl,jk)*zwisoqdif (jl,jk+1,jt)
              zwisoxldif(jl,jk,jt)=zwisoxldif(jl,jk,jt)+zebsh(jl,jk)*zwisoxldif(jl,jk+1,jt)
              zwisoxidif(jl,jk,jt)=zwisoxidif(jl,jk,jt)+zebsh(jl,jk)*zwisoxidif(jl,jk+1,jt)
           ENDDO
        ENDDO
     ENDDO
     
  END IF
!---wiso-code-end

     DO 558 jt=1,trlist% ntrac
        IF (trlist% ti(jt)% nvdiff /= 1) CYCLE
        DO 556 jk=klevm1,itop,-1
           DO 554 jl=1,kproma
              zxtdif(jl,jk,jt)=zxtdif(jl,jk,jt)+                       &
                               zebsh(jl,jk)*zxtdif(jl,jk+1,jt)
554        END DO
556     END DO
558  END DO

!
!     ==================================================================
!
!*       3.8        DIFFUSION IMPLICIT COMPUTATIONS FOR TKE
!
     DO 380 jk=ktdia,klev
        DO 381 jl=1,kproma
           zedif(jl,jk)=ztpfac2*ztkevn(jl,jk)
381     END DO
380  END DO
!
     DO 385 jl=1,kproma
        ztcoe(jl)=(zcdum(jl,itop)+zcdum(jl,itopp1))*0.5_dp
        zqdp=1._dp/(papm1(jl,itopp1)-papm1(jl,itop))
        zdisc=1._dp/(1._dp+(zcdum(jl,itop)+zcdum(jl,itopp1))*0.5_dp*zqdp)
        zebsm(jl,itop)=zdisc*(zcdum(jl,itop)+zcdum(jl,itopp1))*0.5_dp*zqdp
        zedif(jl,itop)=zdisc*zedif(jl,itop)
385  END DO
!
     DO 386 jk=itopp1,klev-2
        DO 387 jl=1,kproma
           zqdp=1._dp/(papm1(jl,jk+1)-papm1(jl,jk))
           zfac=ztcoe(jl)*zqdp
           ztcoe(jl)=(zcdum(jl,jk+1)+zcdum(jl,jk))*0.5_dp
           zdisc=1._dp/(1._dp+zfac*(1._dp-zebsm(jl,jk-1))+(zcdum(jl,jk+1)+      &
                     zcdum(jl,jk))*0.5_dp*zqdp)
           zebsm(jl,jk)=zdisc*(zcdum(jl,jk+1)+zcdum(jl,jk))*0.5_dp*zqdp
           zedif(jl,jk)=zdisc*(zedif(jl,jk)+zfac*zedif(jl,jk-1))
387     END DO
386  END DO
!
     DO 390 jl=1,kproma
        zqdp=1._dp/(papm1(jl,klev)-papm1(jl,klevm1))
        zfac=ztcoe(jl)*zqdp
        ztcoe(jl)=(zcdum(jl,klev)+zcdum(jl,klevm1))*0.5_dp
        zdisc=1._dp/(1._dp+zfac*(1._dp-zebsm(jl,klev-2))+(zcdum(jl,klev)+       &
                  zcdum(jl,klevm1))*0.5_dp*zqdp)
        zedif(jl,klevm1)=zdisc*((zcdum(jl,klev)+zcdum(jl,klevm1))      &
                         *0.5_dp*zqdp*zedif(jl,klev)+zedif(jl,klevm1)     &
                         +zfac*zedif(jl,klev-2))
390  END DO
!
     DO 392 jk=klev-2,itop,-1
        DO 393 jl=1,kproma
           zedif(jl,jk)=zedif(jl,jk)+zebsm(jl,jk)*zedif(jl,jk+1)
393     END DO
392  END DO
!
!*    TIME INTEGRATION OF TURBULENT KINETIC ENERGY AND CHECK
!
     DO 394 jk=itop,klev
        ztest=0._dp
        DO 395 jl=1,kproma
           ptke(jl,jk)=zedif(jl,jk)+ztpfac3*ztkevn(jl,jk)
           ztest=ztest+MERGE(1._dp,0._dp,ptke(jl,jk)<0._dp)
395     END DO
        IF(ztest.NE.0._dp) CALL finish('vdiff','TKE IS NEGATIVE')
394  END DO
!
!*    TIME FILTER FOR TURBULENT KINETIC ENERGY
!
     IF(.NOT.lstart) THEN
       zeps=eps
     ELSE
       zeps=0._dp
     END IF
     DO 397 jk=ktdia,klev
       DO 396 jl=1,kproma
         ptkem1(jl,jk)=                                               &
          ptkem(jl,jk)+zeps*(ptkem1(jl,jk)-2._dp*ptkem(jl,jk)+ptke(jl,jk))
         ptkem(jl,jk)=ptke(jl,jk)
396     END DO
397  END DO
!
!     ------------------------------------------------------------------
!
!*         4.     DIFFUSION IMPLICIT COMPUTATIONS FOR MOMENTUM.
!
!*         4.1     SETTING OF RIGHT HAND SIDES.
!
!!$!411, 412
   DO jk=itop,klevm1
      DO  jl=1,kproma
         zudif(jl,jk)=ztpfac2*pum1(jl,jk)
         zvdif(jl,jk)=ztpfac2*pvm1(jl,jk)
      END DO
   END DO
         
   DO 413 jl=1,kproma
         zudif(jl,klev)=ztpfac2*(pum1(jl,klev)-   &
                       pocu(jl)*(1._dp-pfrl(jl)))
         zvdif(jl,klev)=ztpfac2*(pvm1(jl,klev)-   &
                       pocv(jl)*(1._dp-pfrl(jl)))
413 ENDDO
  
!*         4.2     TOP LAYER ELIMINATION.
!421
     DO jl=1,kproma
        zqdp=1._dp/(paphm1(jl,itopp1)-paphm1(jl,itop))
        zdisc=1._dp/(1._dp+zcfm(jl,itop)*zqdp)
        zebsm(jl,itop)=zdisc*(zcfm(jl,itop)*zqdp)
        zudif(jl,itop)=zdisc*zudif(jl,itop)
        zvdif(jl,itop)=zdisc*zvdif(jl,itop)
  END DO
!
!*         4.3     ELIMINATION FOR MIDDLE LAYERS.
!432, 431
     DO  jk=itopp1,klevm1
        DO  jl=1,kproma
           zqdp=1._dp/(paphm1(jl,jk+1)-paphm1(jl,jk))
           zfac=zcfm(jl,jk-1)*zqdp
           zdisc=1._dp/(1._dp+zfac*(1._dp-zebsm(jl,jk-1))+zcfm(jl,jk)*zqdp)
           zebsm(jl,jk)=zdisc*(zcfm(jl,jk)*zqdp)
           zudif(jl,jk)=zdisc*(zudif(jl,jk)+zfac*zudif(jl,jk-1))
           zvdif(jl,jk)=zdisc*(zvdif(jl,jk)+zfac*zvdif(jl,jk-1))
     END DO
  END DO
!
!*         4.4     BOTTOM LAYER ELIMINATION.
!441
     DO jl=1,kproma
        zqdp=1._dp/(paphm1(jl,klevp1)-paphm1(jl,klev))
        zfac=zcfm(jl,klevm1)*zqdp
        zdisc=1._dp/(1._dp+zfac*(1._dp-zebsm(jl,klevm1))+zcfm(jl,klev)*zqdp)
        zudif(jl,klev)=zdisc*(zudif(jl,klev)+zfac*zudif(jl,klevm1))
        zvdif(jl,klev)=zdisc*(zvdif(jl,klev)+zfac*zvdif(jl,klevm1))
  END DO
!
!*         4.5     BACK-SUBSTITUTION.
!452, 451
     DO  jk=klevm1,itop,-1
        DO  jl=1,kproma
           zudif(jl,jk)=zudif(jl,jk)+zebsm(jl,jk)*zudif(jl,jk+1)
           zvdif(jl,jk)=zvdif(jl,jk)+zebsm(jl,jk)*zvdif(jl,jk+1)
     END DO
  END DO
!
!*         4.6     INCREMENTATION OF U AND V TENDENCIES AND STORAGE OF
!*                 THE DISSIPATION.
!461
     DO  jl=1,kproma
        zvidis(jl)=0._dp
  END DO
!471, 462
     DO  jk=itop,klev
        DO  jl=1,kproma
           zdudt=(zudif(jl,jk)-ztpfac2*pum1(jl,jk))*zcons13
           pvom(jl,jk)=pvom(jl,jk)+zdudt
           zdvdt=(zvdif(jl,jk)-ztpfac2*pvm1(jl,jk))*zcons13
           pvol(jl,jk)=pvol(jl,jk)+zdvdt
           zdis(jl,jk)=0.5_dp*((ztpfac2*pum1(jl,jk)-zudif(jl,jk))         &
                           *(ztpfac4*pum1(jl,jk)+zudif(jl,jk))         &
                           +(ztpfac2*pvm1(jl,jk)-zvdif(jl,jk))         &
                           *(ztpfac4*pvm1(jl,jk)+zvdif(jl,jk)))
           zvidis(jl)=zvidis(jl)+                                      &
                      zdis(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
     END DO
  END DO
    
! HAND IT OVER (JSBACH - MO_SURFACE----------------------

     DO jl=1,kproma
        udif(jl,xjrow) = zudif(jl,klev)
        vdif(jl,xjrow) = zvdif(jl,klev)
! BLM dissipation -----------------------------------------
        pvdis(jl)=pvdis(jl)+zdtime*zcons15*zvidis(jl)

     END DO

!---------------------------------------------------------------------------------------------------

!
!*         5.6     INCREMENTATION OF T AND Q TENDENCIES.
!
560  CONTINUE
!
     DO 571 jk=itop,klev
        DO 561 jl=1,kproma
           zqdif(jl,jk)=zqdif(jl,jk)+ztpfac3*pqm1(jl,jk)
           zdqdt=(zqdif(jl,jk)-pqm1(jl,jk))*zcons13
           pqte(jl,jk)=pqte(jl,jk)+zdqdt
           ztdif(jl,jk)=ztdif(jl,jk)+ztpfac3*zcptgz(jl,jk)
           zdtdt=((ztdif(jl,jk)+zdis(jl,jk)-pgeom1(jl,jk))             &
                 /(cpd*(1._dp+vtmpc2*zqdif(jl,jk)))-ptm1(jl,jk))*zcons13
           ptte(jl,jk)=ptte(jl,jk)+zdtdt
           zxldif(jl,jk)=zxldif(jl,jk)+ztpfac3*pxlm1(jl,jk)
           zxidif(jl,jk)=zxidif(jl,jk)+ztpfac3*pxim1(jl,jk)
           zdxlmdt=(zxldif(jl,jk)-pxlm1(jl,jk))*zcons13
           zdximdt=(zxidif(jl,jk)-pxim1(jl,jk))*zcons13
           pxlte(jl,jk)=pxlte(jl,jk)+zdxlmdt
           pxite(jl,jk)=pxite(jl,jk)+zdximdt
           pxvar(jl,jk)=zvardif(jl,jk)+ztpfac3*pxvar(jl,jk)
           pvdiffp(jl,jk)=zdqdt+zdxlmdt+zdximdt !store for production
561     END DO
571  END DO

!---wiso-code
  IF (lwiso) THEN

!    Incrementation Of Q Tendencies - Water Isotopes

     DO jt=1,kwiso
        DO jk=itop,klev
           DO jl=1,kproma
              zwisoqdif(jl,jk,jt)=zwisoqdif(jl,jk,jt)+ztpfac3*pwisoqm1(jl,jk,jt)
              zwisodqdt=(zwisoqdif(jl,jk,jt)-pwisoqm1(jl,jk,jt))*zcons13
              pwisoqte(jl,jk,jt)=pwisoqte(jl,jk,jt)+zwisodqdt
              zwisoxldif(jl,jk,jt)=zwisoxldif(jl,jk,jt)+ztpfac3*pwisoxlm1(jl,jk,jt)
              zwisoxidif(jl,jk,jt)=zwisoxidif(jl,jk,jt)+ztpfac3*pwisoxim1(jl,jk,jt)
              zwisodxlmdt=(zwisoxldif(jl,jk,jt)-pwisoxlm1(jl,jk,jt))*zcons13
              zwisodximdt=(zwisoxidif(jl,jk,jt)-pwisoxim1(jl,jk,jt))*zcons13
              pwisoxlte(jl,jk,jt)=pwisoxlte(jl,jk,jt)+zwisodxlmdt
              pwisoxite(jl,jk,jt)=pwisoxite(jl,jk,jt)+zwisodximdt
           ENDDO
        ENDDO
      ENDDO

  END IF      
!---wiso-code-end

     ! Correction of tte for snow melt
     DO jl = 1,kproma
        ptte(jl,klev)=ptte(jl,klev)-ztte_corr(jl)
     END DO
!
     zdxtdt = 0._dp
     IF (trlist% anyvdiff /= 0) THEN
        DO jt=1,trlist% ntrac
           IF (trlist% ti(jt)% nvdiff /= 1) CYCLE
           DO jk=itop,klev
              DO jl=1,kproma
                 zxtdif(jl,jk,jt)=zxtdif(jl,jk,jt)+                    &
                                  ztpfac3*pxtm1(jl,jk,jt)
                 zdxtdt(jl,jk) = zxtdif(jl,jk,jt)-pxtm1(jl,jk,jt)
              END DO
           END DO
           IF (jt == ico2idx) THEN             ! ideal mixing in the PBL for CO2
              DO jl=1,kproma
                 ktop = MIN(MAX(ihpbl(jl),iblmax),iblmin) 
                 zdxtdtsum = 0._dp
                 DO jk=ktop,klev
                    zdxtdtsum = zdxtdtsum + zdxtdt(jl,jk) * (paphm1(jl,jk+1)-paphm1(jl,jk))
                 END DO
                 zdxtdtmean = zdxtdtsum / (paphm1(jl,klevp1)-paphm1(jl,ktop))
                 DO jk=ktop,klev
                    zdxtdt(jl,jk) = zdxtdtmean
                 END DO
              END DO
           END IF
           pxtte(1:kproma,1:klev,jt) = pxtte(1:kproma,1:klev,jt) + zdxtdt(1:kproma,1:klev) * zcons13
        END DO
     END IF
!
! back out moisture flux
!
    DO jl=1,kproma
       zqflux(jl,itop)=0._dp 
       zvarpr(jl,itop)=0._dp
       zrho(jl,klevp1)=paphm1(jl,klevp1)/(rd*ztvh_new(jl))        !surf density
       zdqtot=(pqm1(jl,klev)+zx(jl,klev))- zqsurf_new(jl)
       zqshear(jl,klev)=zdqtot*g/pgeom1(jl,klev) !qshear in lev 19

    ENDDO !jl

    DO jk=itop+1,klevp1
       DO jl=1,kproma
          IF (jk<klevp1) THEN
             ztvh=(ptm1(jl,jk)*(1._dp+vtmpc1*pqm1(jl,jk)-zx(jl,jk))      &
                         +ptm1(jl,jk-1)                               &
                            *(1._dp+vtmpc1*pqm1(jl,jk-1)-zx(jl,jk-1)))/2._dp
             zrho(jl,jk)=paphm1(jl,jk)/(rd*ztvh)
          ENDIF
          zrhodz=-(paphm1(jl,jk)-paphm1(jl,jk-1))/g
          zqflux(jl,jk)=zrhodz*pvdiffp(jl,jk-1)+zqflux(jl,jk-1)
          zvarpr(jl,jk)=zqshear(jl,jk-1)*zqflux(jl,jk)/zrho(jl,jk)
       ENDDO !jl
    ENDDO !jk


    DO jk=itop,klev
       DO jl=1,kproma
          pvdiffp(jl,jk)=(zvarpr(jl,jk)+zvarpr(jl,jk+1))/2._dp
          zhexp=EXP(1._dp-pgeom1(jl,jk)/pgeom1(jl,ihpbl(jl)))
          zlam=zzzlam+(zcons3-zzzlam)*zhexp
          IF(jk.GE.ihpbl(jl)) THEN
             zcons23=zcons25
          ELSE
             zcons23=zcons2/zlam
          END IF
          z2geomf=2._dp*pgeom1(jl,jk)
          zz2geo=zcons2*z2geomf
          zmix=zz2geo/(1._dp+zcons23*z2geomf)
          IF(jk.EQ.1) THEN
             ztkesq=SQRT(MAX(ztkemin,ptkem1(jl,1)))
          ELSE
             ztkesq=SQRT(MAX(ztkemin,0.5_dp*(ptkem1(jl,jk-1)             &
                                           +ptkem1(jl,jk))))
          END IF
          pvmixtau(jl,jk)=ztkesq/(zmix*zda1)
       ENDDO !jl
    ENDDO !jk
!
!     ------------------------------------------------------------------
!
!*         6.     NECESSARY COMPUTATIONS IF SUBROUTINE IS BY-PASSED.
!
600  CONTINUE
!
  ELSE
     DO  601 jl=1,kproma
        pevapot(jl)=0._dp
        pqhfla(jl)=0._dp
601  END DO

!---wiso-code
  IF (lwiso) THEN

!   back out moisture flux - water isotopes

     DO jt=1,kwiso
      DO jl=1,kproma
        pwisoevapot(jl,jt)=0._dp
        pwisoqhfla(jl,jt)=0._dp
      END DO
     END DO

  END IF
!---wiso-code-end

  END IF

!
!     ------------------------------------------------------------------
!
  RETURN
END SUBROUTINE vdiff
