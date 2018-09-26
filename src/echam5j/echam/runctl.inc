! runctl.inc
! 
!     ------------------------------------------------------------------
NAMELIST /runctl/ &
  nisp,           &! logical unit for initial spectral fields
  nigp,           &! logical unit for initial grid point fields
  ndiahdf,        &! logical unit for hdiff diagnostics.
  nist,           &! logical unit for surf.temp. file
  nice,           &! logical unit for amip2 sea ice file
  na2stat,        &! logical units for amip2 global statistics
  nhf1,           &! logical units for fourier history files
  nhg1,           &! 
  nhg2,           &! 
  nhg3,           &! 
  nhgl1,          &! logical unit for grid point slt work file
  nfl1,           &! logical unit for optional file read at nstep=0
  nfl2,           &! logical unit for optional file read at resume time
  njin,           &! logical unit for "jobn" input file
  njout,          &! logical unit for "subjobn" output file
  ldebugio,       &! switch debug I/O
  ldebugmem,      &! switch debug memory
  lresume,        &! switch resume model from restart files
  ldebugev,       &! switch debugging events
  ldebughd,       &! switch debugging hd model
  subjob_cmd,     &! subjob command
  loldrerun,      &! switch old rerun-file name-convention (default=F)
  ltimer,         &! switch use timer
  ltctest,        &! switch to test time control
  out_datapath,   &! name of data path with ending '/'
  out_expname,    &! prefix of file names with experiment name
  out_filetype,   &! output file type (GRIB=1, NetCDF=2)
  trac_filetype,  &! output file type (GRIB=1, NetCDF=2)
  stdout_redir,   &! standard output redirection flag
  stderr_redir,   &! standard error  redirection flag
  rerun_filetype, &! rerun file type (NetCDF=2, HDF=3)
  lnwp,           &! switch NWP mode on/off
  lmidatm,        &! switch middle atmosphere model version on/off
  lmlo,           &! switch mixed layer ocean on/off
  lcouple,        &! switch ocean coupling on/off
  lipcc,          &! switch for using IPCC parameters
  lamip,          &! switch AMIP sst on/off
  ldiagamip,      &! switch AMIP diagnostics on/off
  ldebug,         &! switch mass fixer diagnostics on/off
  lport,          &! switch on random perturbation for port test
  nproma,         &! working dimension for grid-point computations
  nproca,         &! number of processors for set A
  nprocb,         &! number of processors for set B
  ltdiag,         &! switch tendency diagnostic on/off -> mo_diag_tendency
  lnmi,           &! switch NMI on/off      -> mo_nmi
  lnudge,         &! switch nudging on/off  -> mo_nudging
  l_orbvsop87,    &! use orbit functions from vsop87 (T) or pcmdi (AMIP) (F)
  ly360,          &! use 360 days year
  numfl1,         &! number of optional fields read at LSTART
  numfl2,         &! number of optional fields read at LSTART/LRESUME
  trigfiles,      &! compose new filenames
  putdata,        &! PUT model DATA to file
  subflag,        &! bind output stream subjobs
  nsub,           &! number of subjobs
  trigjob,        &! subjob submit interval
  labort,         &! stop job query after last rerun cycle
  dt_start,       &! external defined start date
  dt_resume,      &! external defined restart date
  dt_stop,        &! external defined stop date
  delta_time,     &! external defined time step
  no_cycles,      &! stop after number of rerun cycles
  no_steps,       &! stop after number of steps
  no_days,        &! stop after number of days
  putocean,       &! receive ocean data coupling intervals
  getocean,       &! send atmospheric data coupling intervals
  putrerun,       &! store-rerun-files interval
  lhd,            &! for hydrologic discharge model (HD model)
  puthd,          &! coupling interval - put data to   HD-model
  gethd,          &! coupling interval - get data from HD-model
  lhd_que,        &! for comments output from HD model
  nhd_diag,       &! number of region for HD model diagnostics
  iadvec,         &! selects advection scheme: 
                   ! 0 - NONE, 1 - SLT, 2 - SPITFIRE, 3 - FFSL
  lso4,           &! switch for sulfate
  lco2,           &! switch for interactive CO2 (as tracer)
  lco2_flxcor,    &! switch for CO2 flux correction (if lco2=.TRUE.)
  lco2_2perc,     &! switch for workaround to limit CO2 change to 2% of concentration
  ico2_emis,      &! Switch for CO2 emissions
  lco2_nudge,     &! Switch for CO2 nudging
  co2_nudge_vmr,  &! Target value for CO2 nudging (ppm)
  co2_nudge_tau,  &! Relaxation time for CO2 nudging (seconds)
  lsolc,          &! switch for variable solar constant
  lreff,          &! switch for effective radius (volcanic contr. in the
                   ! stratosphere
  lslm,           &! additional array for 1/0 SLM in initial file
  l_volc           ! volcanic aerosols on/off
!     ------------------------------------------------------------------


