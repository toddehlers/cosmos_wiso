MODULE mo_control

  ! Control variables for model housekeeping.
  !
  ! U. Schlese, DKRZ, December 1994
  ! A. Rhodin, MPI, January 1999,
  !      Subroutine m_control renamed to alloc_mods and moved from
  !      module mo_control to module m_alloc_mods.
  ! L. Kornblueh, MPI, June 1999,
  !      added nproca and nprocb for driving the parallel decomposition
  ! M. Esch, MPI, June 1999, ECHAM5-modifications
  ! I. Kirchner, MPI, December 2000, time control
  ! I. Kirchner, MPI, March 2001, revision
  ! M. Esch, MPI, September 2002, add switch for mixed layer ocean
  ! M. Esch, MPI, November  2003, add switch for scenario runs
  ! ------------------------------------------------------------------

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  REAL(dp), POINTER :: vct(:)=>NULL() ! vertical coefficients table.

  INTEGER, SAVE :: nproma   = -1  !   working dimension for grid-point computations
  INTEGER, SAVE :: nproca   = 1   !   number of processors in set A
  INTEGER, SAVE :: nprocb   = 1   !   number of processors in set A
  INTEGER :: nm             !   max zonal wave number.
  INTEGER :: nn             !   max meridional wave number for m=0.
  INTEGER :: nk             !   max meridional wave number.
  INTEGER :: ngl            !   number of gaussian latitudes.
  INTEGER :: nlon           !   max number of points on each latitude line.
  INTEGER :: nlev           !   number of vertical levels.
  INTEGER :: nmp1           !   max zonal wave number + 1.
  INTEGER :: nnp1           !   max meridional wave number + 1.
  INTEGER :: nkp1
  INTEGER :: n2mp1          !   2 * (max zonal wave number + 1).
  INTEGER :: n4mp1          !   4 * (max zonal wave number + 1).
  INTEGER :: nlp2           !   max number of points per latitude line + 2.
  INTEGER :: nlevp1         !   *nlev+1.
  INTEGER :: nsp            !   number of spectral coefficients.
  INTEGER :: n2sp           !   2*number of spectral coefficients.
  INTEGER :: nhgl           !   (number of gaussian latitudes)/2.
  INTEGER :: nscan          !   current scan number.

  INTEGER :: nspace1        !   memory manager space for use of root task
  INTEGER :: nspace2        !   memory manager space for use of subtasks

  INTEGER :: nspadd         !   memory manager space increase

  INTEGER :: nvclev         !   number of levels with vertical coefficients.
  INTEGER, SAVE :: numfl1  = 0    !   number of optional fields read at nstep=0
  INTEGER, SAVE :: numfl2  = 0    !   number of optional fields read at nstep=nresum

  LOGICAL, SAVE :: ldebug    = .FALSE. !   .true. for mass fixer diagnostics
  LOGICAL, SAVE :: ldailysst = .FALSE. !   .true. for using daily SST and SIC
  LOGICAL, SAVE :: lamip     = .FALSE. !   .true. for using variable sst
  LOGICAL, SAVE :: ldiagamip = .FALSE. !   .true. for AMIP diagnostics
  LOGICAL, SAVE :: lcouple   = .FALSE. !   .true. for a coupled run
  LOGICAL, SAVE :: lipcc     = .FALSE. !   .true. for run using IPCC parameters
  LOGICAL, SAVE :: lnwp      = .FALSE. !   .false. for climate mode .true. for NWP mode
  LOGICAL, SAVE :: lnudge    = .FALSE. !   .true. for Nudging mode
  LOGICAL, SAVE :: lmidatm   = .FALSE. !   .true. for middle atmosphere model version
  LOGICAL, SAVE :: lmlo      = .FALSE. !   .true. for mixed layer ocean

  LOGICAL, SAVE :: lprint_m0 = .FALSE. !   .false. for printing t(m0) and timestep time in stepon 
  LOGICAL, SAVE :: lnmi      = .FALSE. !   .true. normal mode initialisation
  LOGICAL, SAVE :: ltdiag    = .FALSE. !   .true. run with additional tendency diagnostics
  LOGICAL, SAVE :: lspit     = .TRUE.  !   .true. (default) for spitfire switched on
  LOGICAL, SAVE :: lcolumn   = .FALSE. !   .true. column model mode
  LOGICAL, SAVE :: lvctch    = .FALSE. !   .true. if column model has changed vct
  LOGICAL, SAVE :: lslm      = .FALSE. !   .true. if additional 1/0 SLM is included in initial file

  INTEGER, SAVE :: nhd_diag  = 0       !   number of region for HD model diagnostics
  LOGICAL, SAVE :: lso4      = .FALSE. !   switch for so4 (sulfate aerosol)
  LOGICAL, SAVE :: lco2      = .FALSE. !   switch for interactive CO2 transport (as tracer)
  LOGICAL, SAVE :: lco2_flxcor = .TRUE.!   switch for CO2 flux correction (if lco2=.TRUE.)
  LOGICAL, SAVE :: lco2_2perc = .FALSE. !  switch for workaround to limit CO2 change to 2% of concentration
  INTEGER, SAVE :: ico2_emis = 0       !   Should CO2 emissions be used?
  LOGICAL, SAVE :: lsolc     = .FALSE. !   switch for variable solar constant
  LOGICAL, SAVE :: lreff     = .FALSE. !   switch for effective radius (volcanic contr. in the
                                       !   stratosphere
  LOGICAL,  SAVE :: lco2_nudge = .FALSE.          ! Nudging of CO2
  REAL(dp), SAVE :: co2_nudge_vmr = 278.0E-06_dp  ! Target value for CO2 nudging
  REAL(dp), SAVE :: co2_nudge_tau = 0._dp         ! tau for co2 nudging

  LOGICAL, SAVE :: l_volc    = .FALSE. !   switch for volcanic forcing

  LOGICAL :: lhd      = .FALSE.  !   .true. for hydrologic discharge model
  LOGICAL :: lhd_que  = .FALSE.  !   .true. for comments output from HD model


  ! Spectral and grid point initial files.
  INTEGER :: nisp  = 23   !  *nisp*      logical unit for initial spectral fields.
  INTEGER :: nigp  = 24   !  *nigp*      logical unit for initial grid point fields.

  ! Climate sea surface temperature and sea ice annual cycle file
  INTEGER :: nist  = 20   !  *nist*      logical unit for surf.temp. file
  INTEGER :: nice  = 96   !  *nice*      logical unit for amip ice file

  ! Climate flux correction file
  INTEGER :: nflu  = 42   !  *nflu*      logical unit for flux correction file

  ! Climate leaf area index and vegetation ratio annual cycle file
  INTEGER :: nvltcl   = 90 !  *nvltcl*    logical unit for climate leaf area index
  INTEGER :: nvgratcl = 91 !  *nvgratcl*  logical unit for climate vegetation ratio

  ! Climate land surface temperature annual cycle file
  INTEGER :: ntslcl   = 92 !  *ntslcl*    logical unit for climate land surface temperature

  ! optional files
  INTEGER :: ndiahdf = 10 !  *ndiahdf*   logical unit for hdiff diagnostics.
  INTEGER :: nfl1    = 13 !  *nfl1*      logical unit for optional file read at nstep=0
  INTEGER :: nfl2    = 14 !  *nfl2*      logical unit for optional file read at resume time
  INTEGER :: na2stat = 77 !  logical unit for AMIP2 global statistics (formatted)
  INTEGER :: na2stre = 78 !  logical unit for AMIP2 global statistics rerun file

  ! History files.
  INTEGER :: nhf1  = 31   !  *nhf1*
  INTEGER :: nhgl1 = 32   !  *nhgl1*     logical unit for grid point slt work file
  INTEGER :: nhg1  = 35   !  *nhg1*
  INTEGER :: nhg2  = 36   !   *to*       logical units for grid point history files.
  INTEGER :: nhg3  = 37   !

  INTEGER :: nsurface = 79   !  surface  log unit for boundary atm-surface jsbach

  ! Subjob files "jobn" and "subjobn"
  INTEGER :: njin   = 30  !  *njin*      logical unit for "jobn" input file
  INTEGER :: njout  = 39  !  *njout*     logical unit for "subjobn" output file

  ! Switches
  LOGICAL :: ltimer    = .FALSE. !  *ltimer*    *true* to use timer
  LOGICAL :: ldebugio  = .FALSE. !  *ldebugio*  *true* to debug IO
  LOGICAL :: ldebugmem = .FALSE. !  *ldebugmem* *true* to debug memory
  LOGICAL :: ldebughd  = .FALSE. !  *ldebughd*  *true* to debug hd model
  LOGICAL :: loldrerun = .FALSE. !  *loldrerun* *true* for old filenames 
  LOGICAL :: ltctest   = .FALSE. !  *ltctest*   *true* to test time control

  ! Special variables
  CHARACTER(256) :: subjob_cmd = 'qsub'  ! *subjob_cmd*  command for submit jobs

  ! output redirection
  INTEGER :: stdout_redir = 0
  INTEGER :: stderr_redir = 0

END MODULE mo_control
