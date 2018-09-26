!*    namelist *COLUMNCTL* controls column model
!
!     A Rhodin MPI Hamburg Nov 1999
!
!     Namelist parameters for single column model (modules/mo_column)
!
NAMELIST /columnctl/ &
!
! name        type       purpose
! ----        ----       -------
  comode,  &! CHARACTER  specify the mode of the column model run:
            !   (len=8)                                         (default='')
            !            ''        3d model run
            !            'traject' write trajectory (run 3d model)
            !            'resid'   write residui
            !            'force'   force column model to full model traject.
            !            'add'     add residui to column model trajectory
            !            'free'    run the column model without forcing
  lrewind, &! LOGICAL    rewind forcing and residui file before reading
            !            for constant forcing               (default=.false.)
  nforce,  &! INTEGER    time steps to force (for add,free) (default=1)
  lvctch,  &! LOGICAL    read vct table from forcing file   (default=.false.)
            !
            !            specify the location(s) of the column model run:
  lat_1d,  &! INTEGER    global latitude  index (traject':up to 100 indices)
  lon_1d,  &! INTEGER    global longitude index (traject':up to 100 indices)
            !            or specify a region ('traject' mode only)
  lat_1db, &! INTEGER    global latitude  indices (lower, upper bound)
  lon_1db, &! INTEGER    global longitude indices (lower, upper bound)
            !
            !            specify the pressure gradient ('free' mode only)
  ug_1d,   &! REAL       geostrophic wind component              (default=0.)
  vg_1d,   &! REAL       geostrophic wind component              (default=0.)
  k_vg,    &! INTEGER    take vg above this level from 'traject' (default=0)
  lfres,   &! LOGICAL    T for forcing in resid-mode             (default=T)
  uvfac,   &! LOGICAL    /=0 if u,v are scaled by sqcst(irow)
            !
            !            flags for column model integration:
  int_t,   &! INTEGER    temperature                              (default=1)
  int_ps,  &! INTEGER    log surface pressure                     (default=1)
  int_uv,  &! INTEGER    wind fields                              (default=1)
  int_vtr, &! INTEGER    vertical transport of q,xl,xi,..,tracers (default=1)
            !
            !            flags for application of forcing fields
            !
  for_t,   &! INTEGER    temperature                        (default=1)
  for_ps,  &! INTEGER    log surface pressure               (default=1)
  for_uv,  &! INTEGER    wind fields                        (default=1)
  for_q,   &! INTEGER    specific humidity                  (default=0)
  for_x,   &! INTEGER    cloud water content                (default=0)
  for_xt,  &! INTEGER    tracer                             (default=0)
            !
  for_test,&! INTEGER    test switches                      (default= 0)
  for_def, &! INTEGER    default value                      (default= 1)
  for_ext, &! INTEGER    external forcing switches          (default= 0)
  for_all, &! INTEGER    default for flags not activated    (default= 0)
            !
  for_d,   &! INTEGER    divergence                         (default=for_def)
  for_vo,  &! INTEGER    vorticity                          (default=for_def)
  for_dt,  &! INTEGER    horizontal temperature gradient    (default=for_def)
  for_dp,  &! INTEGER    log surface pressure gradient      (default=for_def)
  for_tr,  &! INTEGER    qe,xe,xte tendencies by transport  (default=for_def)
  for_trp, &! INTEGER    surf.press. modific. by transport  (default=for_def)
  for_trt, &! INTEGER    xte tendencies by transport        (default=for_def)
  for_zm,  &! INTEGER    zonal means (ul, u0, du0)          (default=for_def)
  for_pol, &! INTEGER    polfilter tendencies in physc      (default=for_def)
            !
  for_tte, &! INTEGER    temperature tendency term                (default=0)
  for_uvt, &! INTEGER    wind tendency terms                      (default=0)
  for_pst, &! INTEGER    surface pressure tendency                (default=0)
  for_pt,  &! INTEGER    vertical velocity                        (default=0)
  for_sst, &! REAL       prescribe SeaSurfaceTemperature          (default=0)
  chk_q_0, &! INTEGER    check that specific hum. is >= 0         (default=0)
            !
  cnudguv, &! REAL(nlev) nudging weight for u,v                   (default=0)
  cnudgt,  &! REAL(nlev) nudging weight for temperature           (default=0)
  cnudgp,  &! REAL       nudging weight for log sfc pressure      (default=0)
  cnudgq,  &! REAL(nlev) nudging weight for specific humidity     (default=0)
  cnudgx,  &! REAL(nlev) nudging weight for cloud water + ice     (default=0)
  cnudgxt   ! REAL(nlev) nudging weight for tracers               (default=0)
