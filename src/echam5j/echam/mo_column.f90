MODULE mo_column
!
! This module gathers all routines and variables required to run
! ECHAM in column model mode.
!
! Working principle: 
! All terms which cannot be calculated by the column model are stored 
! on a file by the 3D model and read in later by the column model.
!
! Authors:
!
! A. Rhodin,   MPI, November 1999, original source
! I. Kirchner, MPI, December 2000, time control update
! A. Rhodin,   DWD, March    2003, nproma blocking
!
  !
  ! modules used
  !

  USE mo_kind,          ONLY: dp
#ifdef NAG
  USE f90_unix,         ONLY: flush
#endif
  USE mo_decomposition, ONLY: gdc => global_decomposition, &
                              ldc => local_decomposition
  USE mo_transpose,     ONLY: tag_gather_gp
  USE mo_mpi,           ONLY: p_pe, p_io, p_nprocs,        & 
                              p_bcast, p_send, p_recv,     &
                              p_parallel, p_parallel_io
  USE mo_doctor,        ONLY: nout, nerr
  USE mo_exception,     ONLY: finish
  USE mo_control,       ONLY: vct, nvclev, nlevp1
  USE mo_time_control,  ONLY: lresume, lstart, &
                              delta_time, get_time_step, time_step_len
  USE mo_namelist,      ONLY: position_nml, nnml, POSITIONED
  USE mo_advection,     ONLY: iadvec
  USE mo_filename,      ONLY: find_next_free_unit
  !
  ! Some more modules are accessed by module subroutines to store/restore
  ! module variables:
  !
  ! mo_scan_buffer, mo_memory_g1a, mo_memory_g2a, mo_gaussgrid
  !
  IMPLICIT NONE
  !
  ! public entities
  !
  PRIVATE
  PUBLIC :: lcotra       ! logical   : true if column model or 'traject' run
  PUBLIC :: inicolumn    ! subroutine: read column model namelist
  PUBLIC :: setcolumn    ! subroutine: set module variables within this module
  PUBLIC :: resetcolumn  ! subroutine: deinitialize column model
  PUBLIC :: get_col_ffti ! subroutine: get information for column model
  PUBLIC :: get_col_dyn  ! subroutine: get information for column model
  PUBLIC :: get_col_tran ! subroutine: get information for column model
  PUBLIC :: get_col_pt   ! subroutine: get information for column model
  PUBLIC :: get_col_pol  ! subroutine: get information for column model
  PUBLIC :: cal_col_expl ! subroutine: explicit integration for column model
  PUBLIC :: write_column ! subroutine: write field to plot file
  PUBLIC :: write_column1

  PUBLIC :: lat_1d       ! global latitude  index of column
  PUBLIC :: lon_1d       ! global longitude index of column
  PUBLIC :: lat_1db      ! global latitude  index (edges of area)
  PUBLIC :: lon_1db      ! global longitude index (edges of area)

  PUBLIC :: sst_1d       ! prescribed SST
  PUBLIC :: for_pst      ! surface pressure tendency forcing flag
  PUBLIC :: for_pt       ! pressure tendency forcing flag
  PUBLIC :: pst          ! surface pressure tendency prescribed
  PUBLIC :: pt_1d        ! pressure tendency prescribed
  PUBLIC :: int_vtr      ! flag: vertical transport of q,xl,xi,..,tracers
  !
  ! module variables (namelist columnctl)
  ! 
  INTEGER,PARAMETER:: mc        = 100! max.number of columns in namelist
  INTEGER,PARAMETER:: ml        = 100! max.number of levels (nudging coeffs.)

  CHARACTER(len=8) :: comode    = '' ! mode
  INTEGER          :: lat_1d(mc)= -1 ! global latitude index of column
  INTEGER          :: lon_1d(mc)= -1 ! global longitude index of column
  INTEGER          :: lat_1db(2)= -1 ! global latitude index of column (area)
  INTEGER          :: lon_1db(2)= -1 ! global longitude index of column (area)
  INTEGER          :: nforce    =  1 ! force until timestep nforce
                                     ! flags for column model integration
  LOGICAL          :: lrewind=.FALSE.! rewind before reading each timestep
  INTEGER          :: int_t     = 1  ! temperature
  INTEGER          :: int_ps    = 1  ! log surface pressure
  INTEGER          :: int_uv    = 1  ! wind fields
  INTEGER          :: int_vtr   = 1  ! vertical transport of q,xl,xi,..,xt_
                                     ! flags for column model forcing
  INTEGER          :: for_t     = 1  ! temperature
  INTEGER          :: for_ps    = 1  ! log surface pressure
  INTEGER          :: for_uv    = 1  ! wind fields
  INTEGER          :: for_q     = 2  ! specific humidity
  INTEGER          :: for_x     = 0  ! cloud water content
  INTEGER          :: for_xt    = 0  ! tracer
  INTEGER          :: for_def   = 1  ! default value
  INTEGER          :: for_test  = -1 ! test switch
  INTEGER          :: for_ext   = -1 ! external forcing switch
  INTEGER          :: for_tr    = 1  ! qe,xe,xte tendencies (due to transport) 
  INTEGER          :: for_trp   = 1  ! surf. pressure modifications(transport)
  INTEGER          :: for_trt   = 1  ! xte tendencies (due to transport)
  INTEGER          :: for_d     = 1  ! divergence
  INTEGER          :: for_vo    = 1  ! vorticity
  INTEGER          :: for_dt    = 1  ! temperature gradient
  INTEGER          :: for_dp    = 1  ! log surface pressure gradient
  INTEGER          :: for_zm    = 0  ! zonal means calculated in spectral space
  INTEGER          :: for_pol   = 1  ! polfilter tendencies in physc
  INTEGER          :: for_pst   = 0  ! surface pressure tendency [Pa/s]
  INTEGER          :: for_pt    = 0  ! pressure tendency [Pa/s]
  INTEGER          :: for_sst   = 0  ! Sea Surface Temperature
  INTEGER          :: for_tte   = 0  ! Temperature tendency
  INTEGER          :: for_uvt   = 0  ! Wind tendencies
  REAL(dp)         :: ug_1d     = 0._dp ! geostrophic wind component
  REAL(dp)         :: vg_1d     = 0._dp ! geostrophic wind component
  INTEGER          :: uvfac     =-1  ! 0 if u,v are not scaled by gl_sqcst
  INTEGER          :: k_vg      = 0  ! estimate geostr. wind above this level
  INTEGER          :: chk_q_0   = 0  ! check that specific humidity is >= 0
  LOGICAL          :: lfres = .TRUE. ! T to force in 'resid' mode
  INTEGER          :: for_all   = 0  ! default settings for remaining variables

  REAL(dp)         :: cnudguv(ml)= 0._dp ! nudging weigth for u,v
  REAL(dp)         :: cnudgt (ml)= 0._dp ! nudging weigth for temperature
  REAL(dp)         :: cnudgp     = 0._dp ! nudging weigth for log sfc pressure
  REAL(dp)         :: cnudgq (ml)= 0._dp ! nudging weigth for water vapour
  REAL(dp)         :: cnudgx (ml)= 0._dp ! nudging weigth for liq. water + ice
  REAL(dp)         :: cnudgxt(ml)= 0._dp ! nudging weigth for tracer
  !
  ! other module variables
  ! 
  LOGICAL          :: lcotra=.FALSE.! column model or 'traject' run
  LOGICAL          :: lcnudge=.FALSE.! nudging applied
  INTEGER          :: p_col    = -1 ! PE running the column model
  LOGICAL          :: lfull =.FALSE.! full model run
  INTEGER          :: ncol          ! number of SCM columns (>=1 for 'traject')
  INTEGER          :: nlat          ! number of SCM latitudes
  INTEGER          :: nlon          ! number of SCM longitudes
  INTEGER          :: nlcol         ! number of SCM columns on this pe
  INTEGER          :: iread         ! index of column in input files
  INTEGER,ALLOCATABLE :: pes (:)    ! processor corresponding to SCM locations
  INTEGER,ALLOCATABLE :: lats(:)    ! SCM global latitude  indices
  INTEGER,ALLOCATABLE :: lons(:)    ! SCM global longitude indices
  INTEGER,ALLOCATABLE :: lrow(:)    ! SCM local  latitude  indices
  INTEGER,ALLOCATABLE :: lcol(:)    ! SCM local  longitude indices
  INTEGER          :: unit_f1d = -1 ! fortran unit for forcing file
  INTEGER          :: unit_r1d = -1 ! fortran unit for residuum file
  INTEGER          :: unit_b1d = -1 ! fortran unit for results before correct. 
  INTEGER          :: unit_a1d = -1 ! fortran unit for results after correct.
  LOGICAL          :: fapp  =.FALSE.! forcing applied in this step
  REAL(dp),ALLOCATABLE :: ug    (:,:,:) ! geostrophic wind component profile
  REAL(dp),ALLOCATABLE :: vg    (:,:,:) ! geostrophic wind component profile
  REAL(dp),ALLOCATABLE :: pst   (:,:)   ! surface pressure tendency
  REAL(dp),ALLOCATABLE :: pt_1d (:,:,:) ! pressure tendency [Pa/s]
  REAL(dp),ALLOCATABLE :: sst_1d(:,:)   ! prescribed Sea Surface Temperature
  REAL(dp),ALLOCATABLE :: tte_3d(:,:,:) ! temperat. tendency (polefilter, physc)
  REAL(dp),ALLOCATABLE :: qte_3d(:,:,:) ! humidity  tendency (polefilter, physc)
  !
  ! interfaces
  !
  INTERFACE write_column            ! write field x to plot file:
    MODULE PROCEDURE write_column3  ! ( x(lon,lev,lat) ,name)
    MODULE PROCEDURE write_column2  ! ( x(lon,    lat) ,name)
    MODULE PROCEDURE write_column3r ! ( x(lon,lev)     ,name ,jlat)
    MODULE PROCEDURE write_column2r ! ( x(lon)         ,name ,jlat)
  END INTERFACE

CONTAINS
!==============================================================================
! Initialization routines
!==============================================================================
  SUBROUTINE inicolumn (lcolumn, lvctch, nlev, nlevp1, nvclev, vct)
  LOGICAL ,INTENT (out)   :: lcolumn ! .true. if a column model is run
  LOGICAL ,INTENT (out)   :: lvctch  ! .true. if vct table is changed
  INTEGER ,INTENT (inout) :: nlev    ! levels (changed if lvctch == .true.)
  INTEGER ,INTENT (inout) :: nlevp1  ! levels (changed if lvctch == .true.)
  INTEGER ,INTENT (inout) :: nvclev  ! levels (changed if lvctch == .true.)
  REAL(dp)    ,POINTER        :: vct (:) ! hybrid level coefficients     ( '' )
    REAL(dp) ,POINTER :: x (:)
    INTEGER       :: ierr
    !
    ! Initialize module variables from namelist.
    !
    INCLUDE 'columnctl.inc'
    !
    ! default values
    !
    lat_1d  = -1
    lon_1d  = -1
    lat_1db = -1
    lon_1db = -1
    nforce  =  1
    comode  = ''
    lcolumn = .FALSE.
    lrewind = .FALSE.
    lvctch  = .FALSE.
    int_t   =  1
    int_ps  =  1
    int_uv  =  1
    int_vtr =  1
    for_t   =  1
    for_ps  =  1
    for_uv  =  1
    for_q   =  2
    for_x   =  0
    for_xt  =  0
    for_def =  1
    for_test= -1
    for_ext = -1
    for_d   = -1
    for_vo  = -1
    for_dt  = -1
    for_dp  = -1
    for_tr  = -1
    for_trp = -1
    for_trt = -1
    for_zm  = -1
    for_pol = -1
    for_tte = -1
    for_uvt = -1
    for_pst =  0
    for_pt  =  0
    for_sst =  0
    ug_1d   =  0._dp
    vg_1d   =  0._dp
    uvfac   =  -1
    k_vg    =  0
    lfres   = .TRUE.
    chk_q_0 =  0
    cnudguv =  0._dp
    cnudgt  =  0._dp
    cnudgp  =  0._dp
    cnudgq  =  0._dp
    cnudgx  =  0._dp
    cnudgxt =  0._dp
    !
    ! read namelist
    !
    IF (p_parallel_io) THEN
      CALL position_nml ('COLUMNCTL', status=ierr)
      SELECT CASE (ierr)
      CASE (POSITIONED)
        READ (nnml, columnctl)
      END SELECT
    ENDIF
    IF (p_parallel) THEN
      CALL p_bcast (comode,  p_io)
    ENDIF
    IF (comode == '') RETURN 
    IF (p_parallel) THEN
      CALL p_bcast (lat_1d,  p_io)
      CALL p_bcast (lon_1d,  p_io)
      CALL p_bcast (lat_1db, p_io)
      CALL p_bcast (lon_1db, p_io)
      CALL p_bcast (nforce,  p_io)
      CALL p_bcast (lrewind, p_io)
      CALL p_bcast (lvctch,  p_io)
      CALL p_bcast (int_t,   p_io)
      CALL p_bcast (int_ps,  p_io)
      CALL p_bcast (int_uv,  p_io)
      CALL p_bcast (int_vtr, p_io)
      CALL p_bcast (for_t,   p_io)
      CALL p_bcast (for_ps,  p_io)
      CALL p_bcast (for_uv,  p_io)
      CALL p_bcast (for_q,   p_io)
      CALL p_bcast (for_x,   p_io)
      CALL p_bcast (for_xt,  p_io)
      CALL p_bcast (for_def, p_io)
      CALL p_bcast (for_test,p_io)
      CALL p_bcast (for_ext, p_io)
      CALL p_bcast (for_d,   p_io)
      CALL p_bcast (for_vo,  p_io)
      CALL p_bcast (for_dt,  p_io)
      CALL p_bcast (for_dp,  p_io)
      CALL p_bcast (for_tr,  p_io)
      CALL p_bcast (for_trt, p_io)
      CALL p_bcast (for_trp, p_io)
      CALL p_bcast (for_zm,  p_io)
      CALL p_bcast (for_pol, p_io)
      CALL p_bcast (for_tte, p_io)
      CALL p_bcast (for_uvt, p_io)
      CALL p_bcast (for_pst, p_io)
      CALL p_bcast (for_pt,  p_io)
      CALL p_bcast (for_sst, p_io)
      CALL p_bcast (ug_1d,   p_io)
      CALL p_bcast (vg_1d,   p_io)
      CALL p_bcast (uvfac,   p_io)
      CALL p_bcast (k_vg    ,p_io)
      CALL p_bcast (lfres   ,p_io)
      CALL p_bcast (chk_q_0 ,p_io)
      CALL p_bcast (cnudguv ,p_io)
      CALL p_bcast (cnudgt  ,p_io)
      CALL p_bcast (cnudgp  ,p_io)
      CALL p_bcast (cnudgq  ,p_io)
      CALL p_bcast (cnudgx  ,p_io)
      CALL p_bcast (cnudgxt ,p_io)
      CALL p_bcast (for_all ,p_io)
    ENDIF
    !
    ! set flags for individual variables
    ! for_test: old forcing from 3D ECHAM
    !
    IF (for_test== 1) THEN
      for_def = 0
      for_ext = 0
    ENDIF
    IF (for_d  ==-1) for_d   = for_test
    IF (for_vo ==-1) for_vo  = for_test
    IF (for_dt ==-1) for_dt  = for_test
    IF (for_dp ==-1) for_dp  = for_test
    IF (for_tr ==-1) for_tr  = for_test
    IF (for_trt==-1) for_trt = for_test
    IF (for_trp==-1) for_trp = for_test
    IF (for_zm ==-1) for_zm  = for_test
    IF (for_pol==-1) for_pol = for_test
    IF (uvfac  ==-1) uvfac   = for_test
    !
    ! for_ext: forcing from external data (Analyses etc)
    !
    IF (for_ext== 1) THEN
      for_def = 0
      IF (for_q  ==-1) for_q   = for_ext
      IF (for_x  ==-1) for_x   = for_ext
      IF (for_d  ==-1) for_d   = for_ext
      IF (for_pst==-1) for_pst = for_ext
      IF (for_sst==-1) for_sst = for_ext
      IF (for_tr ==-1) for_tr  = for_ext
      IF (for_trt==-1) for_trt = for_ext
      IF (for_tte==-1) for_tte = for_ext
      IF (uvfac  ==-1) uvfac   = 0
    ENDIF
    !
    ! for_def: forcing from 3D ECHAM
    !
    IF (for_d  ==-1) for_d   = for_def
    IF (for_vo ==-1) for_vo  = 0
    IF (for_dt ==-1) for_dt  = 0
    IF (for_dp ==-1) for_dp  = for_def
    IF (for_tr ==-1) for_tr  = for_def
    IF (for_trp==-1) for_trp = for_def
    IF (for_trt==-1) for_trt = for_def
    IF (for_zm ==-1) for_zm  = 0
    IF (for_pol==-1) for_pol = for_def
    IF (for_tte==-1) for_tte = for_def
    IF (for_uvt==-1) for_uvt = for_def
    IF (uvfac  ==-1) uvfac   = for_def
    !
    ! for_all: default for all forcing flags not activated so far
    !
    IF (for_q  == 0) for_q   = for_all
    IF (for_x  == 0) for_x   = for_all
    IF (for_xt == 0) for_xt  = for_all
    IF (for_d  == 0) for_d   = for_all
    IF (for_vo == 0) for_vo  = for_all
    IF (for_dt == 0) for_dt  = for_all
    IF (for_dp == 0) for_dp  = for_all
    IF (for_tr == 0) for_tr  = for_all
    IF (for_trp== 0) for_trp = for_all
    IF (for_trt== 0) for_trt = for_all
    IF (for_zm == 0) for_zm  = for_all
    IF (for_pol== 0) for_pol = for_all
    IF (for_tte== 0) for_tte = for_all
    IF (for_uvt== 0) for_uvt = for_all
    !
    ! is nudging applied to any variable?
    !
    lcnudge = ANY(cnudguv/=0._dp) .OR. ANY(cnudgt/=0._dp) .OR.     cnudgp /=0._dp &
        .OR. ANY(cnudgq /=0._dp) .OR. ANY(cnudgx/=0._dp) .OR. ANY(cnudgxt/=0._dp)
    !
    ! calculate ncol,nlat,nlon
    ! (number of SCM columns, SCM latitudes, SCM longitudes)
    !
    ncol = 0
    IF (ANY(lat_1d/=-1).OR.ANY(lon_1d/=-1)) THEN  ! individual locations:
      lat_1db    = -1                             ! do not specify areas
      lon_1db    = -1                             ! do not specify areas
      nlat       = COUNT(lat_1d/=-1)              ! SCM latitudes
      nlon       = COUNT(lon_1d/=-1)              ! SCM longitudes
      IF(nlat/=nlon)CALL finish('inicolumn','error in lat_1db/lon_1db')
      ncol=nlat                                   ! SCM columns
    ELSE                                          !
      IF (lat_1db(2)==-1) lat_1db(2) = lat_1db(1) ! edges of areas specified
      nlat = MAX(0,lat_1db(2)-lat_1db(1)+1)       ! SCM latitudes
      IF (ANY(lat_1db==-1)) nlat = 0
      IF (lon_1db(2)==-1) lon_1db(2) = lon_1db(1)
      nlon = MAX(0,lon_1db(2)-lon_1db(1)+1)       ! SCM longitudes
      IF (ANY(lon_1db==-1)) nlon = 0
      ncol  = nlat * nlon                         ! SCM columns
    ENDIF
    !
    ! check value of comode, number of SCM columns
    !
    SELECT CASE (comode)
    CASE ('')
    CASE ('traject')
      IF(ncol==0) CALL finish('inicolumn','location of columns not specified') 
    CASE ('resid','force','add','free')
      IF(ncol/=1) CALL finish('inicolumn',&
                              'exactly one column location must be specified')
    CASE default
      CALL finish('inicolumn','comode='//comode)
    END SELECT
    !
    ! set return arguments
    !
    lcotra  =  comode /= ''                            ! traject or SCM run
    lcolumn = (comode /= '' .AND. comode /= 'traject') ! SCM run
    lvctch  = lvctch .AND. lcolumn                     ! VCT changed by SCM
    IF(lcolumn .AND. iadvec==3) iadvec=0               ! swich off FFSL
    !
    ! change vertical levels
    !
    IF(lvctch) THEN
      NULLIFY (x)
      IF (p_pe==p_io) THEN
        !
        ! read AK,BK coefficients from forcing-file
        !
        unit_f1d = find_next_free_unit (51,100)
        OPEN (unit_f1d,file='forcing' ,form='unformatted')
        CALL read_colx (x, 'AK', unit_f1d)
        nvclev = SIZE(x)
        DEALLOCATE (vct)
        ALLOCATE (vct (2*nvclev))
        vct(1:nvclev) = x
        CALL read_colx (x, 'BK', unit_f1d)
        vct(nvclev+1:nvclev+nvclev) = x
        DEALLOCATE (x)
        CLOSE (unit_f1d); unit_f1d = -1
      ENDIF
      !
      ! broadcast changed nvclev,vct; recalculate nlev,nlevp1
      !
      CALL p_bcast (nvclev, p_io)
      nlev   = nvclev-1
      nlevp1 = nvclev
      IF (p_pe/=p_io) THEN
        DEALLOCATE (vct)
        ALLOCATE (vct (2*nvclev))
      ENDIF
      CALL p_bcast (vct, p_io)            
      !
      ! printout changed vct
      !
      IF (p_pe==p_io) THEN
        WRITE(nout,*)
        WRITE(nout,'(a)') REPEAT('-',72)
        WRITE(nout,'(a)')' Column model has changed hybrid levels:'
        WRITE(nout,'(a)')' this option does not work with restart files.'
        WRITE(nout,'(a)')' this option does not work with nudging.'
        WRITE(nout,*)
        WRITE(nout,'(a,i4)')  '      nlev      ='  ,nlev
        WRITE(nout,'(a,i4)')  '      nlevp1    ='  ,nlevp1
        WRITE(nout,'(a,i4)')  '      nvclev    ='  ,nvclev
        WRITE(nout,*)
        IF(lresume) &
         CALL finish('inicolumn','lvctch=.true. conflicts with lresume=.true.')
      ENDIF
    ENDIF
  END SUBROUTINE inicolumn
!------------------------------------------------------------------------------
  SUBROUTINE setcolumn
  !
  ! Set module variables within this module.
  ! This routine must be called after the parallel decomposition has
  ! been set (after init_decomposition)
  !
    INTEGER :: p, i, j, k, l
    INTEGER :: lon, lat
    INTEGER :: tlat(ncol)
    INTEGER :: tlon(ncol)
    LOGICAL :: llful

    lfull    = .FALSE.
    p_col    = -1
    !
    ! allocate arrays for global indices
    !
    IF(ldc% pe==p_io) THEN
      !
      ! allocate on I/O processor
      !
      ALLOCATE (lats  (ncol))  ;lats  = -1
      ALLOCATE (lons  (ncol))  ;lons  = -1
      ALLOCATE (pes   (ncol))  ;pes   = -1
    ELSE
      !
      ! zero length on non I/O processor
      !
      ALLOCATE (lats  (0))
      ALLOCATE (lons  (0))
      ALLOCATE (pes   (0))
    ENDIF
    iread = 1
    !
    ! loop over pe's
    !
    DO p=1,SIZE(gdc)
      IF(gdc(p)% col_1d)THEN
        !
        ! SCM processor found:
        ! set p_col (SCM processor id, lat/lon indices)
        !
        IF(ncol>1) CALL finish('setcolumn','more than one gp in column model')
        p_col = gdc(p)% pe
        IF(ldc% pe == gdc(p)% pe) THEN
          !
          ! if I am the SCM processor, set llat_1d = 1, llon_1d
          ! (local lat/lon indices)
          !
          nlcol    = 1
          ALLOCATE (lrow (nlcol)) ;lrow = 1
          ALLOCATE (lcol (nlcol)) ;lcol = 1
          IF(ldc% pe==p_io) THEN
            lats = lat_1d(1)
            lons = lon_1d(1)
          ENDIF
        ENDIF
      ENDIF
      IF(.NOT.gdc(p)% col_1d) THEN
        !
        ! for SCM processors (i.e. full model run in traject mode):
        ! set lfull, llat_1d, llon_1d
        !
        l = 0
        i=0; j=0; k=0
        !
        ! loop over SCM columns
        !
        llful=.FALSE.
        DO
         IF (lat_1db(1)/=-1) THEN
           k=MAX(lon_1db(1),k)
           j=MAX(lat_1db(1)-1,j)
           j=j+1
           IF(j>lat_1db(2)) THEN
             j=lat_1db(1)
             k=k+1
             IF(k>lon_1db(2)) EXIT
           ENDIF
           lon=k
           lat=j
         ELSE
           k=k+1
           IF(k>ncol) EXIT
           IF(lon_1d(k)==-1) CYCLE
           lat=lat_1d(k)
           lon=lon_1d(k)
         ENDIF
         i=i+1
         !
         ! check for valid lon/lat indices
         !
         IF(lat<1)        CALL finish('setcolumn','SCM latitude  index < 1')
         IF(lon<1)        CALL finish('setcolumn','SCM longitude index < 1')
         IF(lat>ldc%nlat) CALL finish('setcolumn','SCM latitude  index > nlat')
         IF(lon>ldc%nlon) CALL finish('setcolumn','SCM longitude index > nlon')
         !
         ! determine local lat/lon indices on northern hemisphere
         !
         IF(gdc(p)%glats(1)<=lat .AND. gdc(p)% glate(1)>=lat .AND. &
            gdc(p)%glons(1)<=lon .AND. gdc(p)% glone(1)>=lon) THEN
           lfull   = .TRUE.
           llful   = .TRUE.
           l   = l+1
           tlat(l) = lat - gdc(p)% glats(1) + 1
           tlon(l) = lon - gdc(p)% glons(1) + 1
         ENDIF
         !
         ! determine local lat/lon indices on southern hemisphere
         !
         IF(gdc(p)%glats(2)<=lat .AND. gdc(p)% glate(2)>=lat .AND. &
            gdc(p)%glons(2)<=          gdc(p)% glone(2)      .AND. &
            gdc(p)%glons(2)<=lon .AND. gdc(p)% glone(2)>=lon) THEN
           lfull   = .TRUE.
           llful   = .TRUE.
           l   = l+1
           tlat(l) = lat - gdc(p)% glats(2) + 1 + gdc(p)% nglh(1)
           tlon(l) = lon - gdc(p)% glons(2) + 1
         ENDIF
         !
         ! local lat/lon indices on southern hemisphere, rotated coordinates
         !
         IF(gdc(p)%glats(2)<=lat .AND. gdc(p)% glate(2)>=lat .AND. &
            gdc(p)%glons(2)>           gdc(p)% glone(2)      .AND. &
           (gdc(p)%glons(2)<=lon .OR.  gdc(p)% glone(2)>=lon))THEN
           lfull   = .TRUE.
           llful   = .TRUE.
           l   = l+1
           tlat(l) = lat - gdc(p)% glats(2) + 1 + gdc(p)% nglh(1)
           tlon(l) = lon - gdc(p)% glons(2) + 1
           IF (tlon(l)<1) tlon(l) = tlon(l) + gdc(p)% nlon
         ENDIF
         !
         ! on I/O processor: store global indices, 
         !                   store id's of processors with SCM columns
         !
         IF(ldc% pe==p_io .AND. llful) THEN
           lats(i)  = lat
           lons(i)  = lon
           pes (i)  = gdc(p)% pe
         ENDIF
        END DO ! SCM column loop
        !
        ! store local indices
        !
        IF (ldc% pe==gdc(p)% pe) THEN
          nlcol = l
          ALLOCATE (lrow (nlcol))
          ALLOCATE (lcol (nlcol))
          !
          ! for NPROMA blocking (nproma/=nglon) recalculate indices
          !
          DO l = 1, nlcol
            i = (tlon(l)-1) + (ldc% nglon) * (tlat(l)-1)
            lcol(l) = MOD (i,ldc% nproma) + 1
            lrow(l) =      i/ldc% nproma  + 1
          END DO
        ENDIF
      ENDIF
    END DO ! PE loop
    !
    ! allocate arrays not allocated so far, check parameters, set defaults
    !
    IF (.NOT. ALLOCATED (lrow)) ALLOCATE (lrow (0))
    IF (.NOT. ALLOCATED (lcol)) ALLOCATE (lcol (0))
    IF (for_pst > 0) THEN
      IF (.NOT.ALLOCATED(pst)) ALLOCATE (pst (ldc%nglon,ldc%nglat)); pst = 0._dp
    ENDIF
    IF (for_pt > 0) THEN
      IF (.NOT.ALLOCATED(pt_1d)) &
        ALLOCATE (pt_1d(ldc%nglon,ldc%nlev+1,ldc%nglat))
      pt_1d = 0._dp
    ENDIF
    IF (for_sst > 0) THEN
      IF (.NOT.ALLOCATED(sst_1d)) ALLOCATE (sst_1d (ldc%nglon, ldc%nglat))
      sst_1d = 0._dp
    ENDIF
    IF (for_pol > 0) THEN
      IF (.NOT.ALLOCATED(tte_3d)) THEN
        ALLOCATE (tte_3d (ldc%nglon,ldc%nlev,ldc%nglat))
        ALLOCATE (qte_3d (ldc%nglon,ldc%nlev,ldc%nglat))
        tte_3d = 0._dp
        qte_3d = 0._dp
      ENDIF
    ENDIF
    IF (ldc% col_1d) THEN
      ALLOCATE (ug     (ldc%nglon, ldc%nlev, ldc%nglat)); ug     = ug_1d
      ALLOCATE (vg     (ldc%nglon, ldc%nlev, ldc%nglat)); vg     = vg_1d
      k_vg  = MIN (k_vg , ldc%nlev) 
    ENDIF
    !
    ! printout
    !
    IF (p_col == -1 .AND. .NOT. lfull) RETURN
    IF (p_pe == p_io) THEN
      WRITE(nout,*)
      WRITE(nout,'(a)') REPEAT('-',72)
      WRITE(nout,'(a)')' Column model initialised:'
      WRITE(nout,*)
      WRITE(nout,'(a,1x,a)')'      comode    =', comode
      WRITE(nout,'(a,i4)')  '      p_col     =', p_col
      WRITE(nout,'(a,l4)')  '      lfull     =', lfull
      WRITE(nout,*)
      WRITE(nout,'(a,i4)')  '      ncol      =', ncol
      WRITE(nout,'(a,i4,7i7/(14x,8i7))') &
                            '      jlat      =', lats
      WRITE(nout,'(a,i4,7i7/(14x,8i7))') &
                            '      ilon      =', lons
      WRITE(nout,'(a,8f7.2/(17x,8f7.2))') &
                            '      longitude :',360._dp/ldc%nlon*(lons-1)
      WRITE(nout,'(a,8f7.2/(17x,8f7.2))') &
                            '      approx.lat:',90._dp-180._dp*(lats-0.5_dp)/ldc%nlat
      WRITE(nout,*)
      WRITE(nout,'(a,l4)')  '      lrewind   ='  ,lrewind
      WRITE(nout,*)
      WRITE(nout,'(a,f7.2)')'      ug_1d     =', ug_1d
      WRITE(nout,'(a,f7.2)')'      vg_1d     =', vg_1d
      WRITE(nout,'(a,i4)')  '      uvfac     =', uvfac
      WRITE(nout,'(a,i4)')  '      k_vg      =', k_vg
      WRITE(nout,'(a,l4)')  '      lfres     =', lfres
      WRITE(nout,'(a,i4)')  '      chk_q_0   =', chk_q_0 
      WRITE(nout,*)
      WRITE(nout,'(a,i4)')  '      for_def   =', for_def
      WRITE(nout,'(a,i4)')  '      for_test  =', for_test
      WRITE(nout,'(a,i4)')  '      for_ext   =', for_ext
      WRITE(nout,'(a,i4)')  '      for_all   =', for_all
      WRITE(nout,*)
      WRITE(nout,'(a,i4)')  '      int_t     =', int_t
      WRITE(nout,'(a,i4)')  '      int_ps    =', int_ps
      WRITE(nout,'(a,i4)')  '      int_uv    =', int_uv
      WRITE(nout,'(a,i4)')  '      int_vtr   =', int_vtr
      WRITE(nout,*)
      WRITE(nout,'(a,i4)')  '      for_t     =', for_t
      WRITE(nout,'(a,i4)')  '      for_ps    =', for_ps
      WRITE(nout,'(a,i4)')  '      for_uv    =', for_uv
      WRITE(nout,'(a,i4)')  '      for_q     =', for_q
      WRITE(nout,'(a,i4)')  '      for_x     =', for_x
      WRITE(nout,'(a,i4)')  '      for_xt    =', for_xt
      WRITE(nout,*)
      WRITE(nout,'(a,i4)')  '      for_tr    =', for_tr
      WRITE(nout,'(a,i4)')  '      for_trp   =', for_trp
      WRITE(nout,'(a,i4)')  '      for_trt   =', for_trt
      WRITE(nout,'(a,i4)')  '      for_d     =', for_d
      WRITE(nout,'(a,i4)')  '      for_vo    =', for_vo
      WRITE(nout,'(a,i4)')  '      for_dt    =', for_dt
      WRITE(nout,'(a,i4)')  '      for_dp    =', for_dp
      WRITE(nout,'(a,i4)')  '      for_zm    =', for_zm
      WRITE(nout,'(a,i4)')  '      for_tte   =', for_tte
      WRITE(nout,'(a,i4)')  '      for_uvt   =', for_uvt
      WRITE(nout,'(a,i4)')  '      for_pst   =', for_pst
      WRITE(nout,'(a,i4)')  '      for_pt    =', for_pt
      WRITE(nout,'(a,i4)')  '      for_pol   =', for_pol
      WRITE(nout,'(a,i4)')  '      for_sst   =', for_sst
      WRITE(nout,*)
      WRITE(nout,'(a,l4)')  '      lcnudge   =', lcnudge
      DO i=ml,1,-1
      IF(cnudguv(i)==0.0_dp .AND. cnudgt (i)==0.0_dp .AND. cnudgq(i)==0.0_dp .AND. &
         cnudgx (i)==0.0_dp .AND. cnudgxt(i)==0.0_dp .AND.       i > 1       ) CYCLE
      WRITE(nout,'(a,(10f5.2)/(17x,10f5.2))')'      cnudguv   =', cnudguv(1:i)
      WRITE(nout,'(a,(10f5.2)/(17x,10f5.2))')'      cnudgt    =', cnudgt (1:i)
      WRITE(nout,'(a,   f5.2              )')'      cnudgp    =', cnudgp
      WRITE(nout,'(a,(10f5.2)/(17x,10f5.2))')'      cnudgq    =', cnudgq (1:i)
      WRITE(nout,'(a,(10f5.2)/(17x,10f5.2))')'      cnudgx    =', cnudgx (1:i)
      WRITE(nout,'(a,(10f5.2)/(17x,10f5.2))')'      cnudgxt   =', cnudgxt(1:i)
      EXIT
      END DO
      WRITE(nout,'(a)') REPEAT('-',72)
      WRITE(nout,*)
      !
      ! open files
      !
      IF (p_col/=-1 .OR. lfull) &
        unit_f1d = find_next_free_unit (51,100)
        OPEN (unit_f1d,file='forcing' ,form='unformatted')
      IF (p_col/=-1) THEN
        unit_r1d = find_next_free_unit (52,100)
        unit_b1d = find_next_free_unit (53,100)
        unit_a1d = find_next_free_unit (54,100)
        OPEN (unit_r1d,file='residui' ,form='unformatted')
        OPEN (unit_b1d,file='result1' ,form='unformatted')
        OPEN (unit_a1d,file='result2' ,form='unformatted')
      ENDIF
    ENDIF
    !
    ! security check
    !
    IF (p_col /= -1 .AND. ldc% nproca * ldc% nprocb /= 1) &
      CALL finish ('inicolumn','nproca,nprocb must be 1 in SCM run')
    IF (p_col /= -1 .AND. p_nprocs /= 1) &
      CALL finish ('inicolumn','number of processors must be 1 in SCM run')
  END SUBROUTINE setcolumn
!------------------------------------------------------------------------------
  SUBROUTINE resetcolumn
    !
    ! clean up module variables
    !
    comode = ' '
    lat_1d =  0
    lon_1d =  0
    lfull  = .FALSE.
    p_col  = -1
    ncol   =  0
    nlcol  =  0
    !
    ! close output files
    !
    IF (p_pe == p_io) THEN
      IF(unit_f1d>0) CLOSE (unit_f1d); unit_f1d = -1
      IF(unit_r1d>0) CLOSE (unit_r1d); unit_r1d = -1
      IF(unit_b1d>0) CLOSE (unit_b1d); unit_b1d = -1
      IF(unit_a1d>0) CLOSE (unit_a1d); unit_a1d = -1
    ENDIF
    !
    ! deallocate
    !
    IF (ldc% col_1d) DEALLOCATE (ug, vg)
    IF (ALLOCATED (pst)   ) DEALLOCATE (pst)
    IF (ALLOCATED (pt_1d )) DEALLOCATE (pt_1d)
    IF (ALLOCATED (sst_1d)) DEALLOCATE (sst_1d)
    IF (ALLOCATED (tte_3d)) DEALLOCATE (tte_3d)
    IF (ALLOCATED (qte_3d)) DEALLOCATE (qte_3d)
    IF (ALLOCATED (lats)  ) DEALLOCATE (lats)
    IF (ALLOCATED (lons)  ) DEALLOCATE (lons)
    IF (ALLOCATED (lrow)  ) DEALLOCATE (lrow)
    IF (ALLOCATED (lcol)  ) DEALLOCATE (lcol)
    IF (ALLOCATED (pes)   ) DEALLOCATE (pes)
  END SUBROUTINE resetcolumn
!==============================================================================
! Routines to be called during grid point calculations (scan1sl)
!==============================================================================
  SUBROUTINE get_col_tran
  USE mo_scan_buffer,    ONLY: alps, qte, xlte, xite, xtte, ul
  USE mo_memory_g1a,     ONLY: alpsm1
  !
  ! Prepare or forcing fields for column model.
  ! To be called from scan1sl after the semilagrangian transport
  ! (after CALL slt2)
  !
    INTEGER               :: i
    CHARACTER(len=6)      :: name = 'XTEnnn'
    IF (comode/=' ') THEN
      !
      ! fetch from file or other PEs
      !
      IF(for_zm >0) CALL process_colzm(ul,      'UL'      ,for_zm ==1)
      IF(for_tr >0) CALL process_col3 (qte,     'QE'      ,for_tr ==1)
      IF(for_tr >0) CALL process_col3 (xlte,    'XLTE'    ,for_tr ==1)
      IF(for_tr >0) CALL process_col3 (xite,    'XITE'    ,for_tr ==1)
      IF(for_trp>0) CALL process_col2 (alps,    'ALPS_TR' ,for_trp==1)
      IF(for_trp>0) CALL process_col2 (alpsm1,  'ALPSM1'  ,for_trp==1)
      IF(for_trt>0) THEN
        DO i=1,SIZE(xtte,3)
          WRITE (name(4:6),'(i3.3)') i
            CALL process_col3 (xtte(:,:,i,:), name ,for_trt==1)
        END DO
      ENDIF
    ENDIF
  END SUBROUTINE get_col_tran
!------------------------------------------------------------------------------
  SUBROUTINE get_col_dyn
  USE mo_scan_buffer,  ONLY: vom, vol, tte ! modified tendencies
  !
  ! Prepare or use forcing fields for column model.
  ! To be called from scan1sl after CALL dyn
  !
    IF (comode/=' ') THEN
      !
      ! fetch from other PEs
      !
      IF(for_uvt>0) CALL process_col3 (vom,     'VOM'     ,for_uvt==1)
      IF(for_uvt>0) CALL process_col3 (vol,     'VOL'     ,for_uvt==1)
      IF(for_tte>0) CALL process_col3 (tte,     'TTE'     ,for_tte==1)
    ENDIF
  END SUBROUTINE get_col_dyn
!------------------------------------------------------------------------------
  SUBROUTINE get_col_ffti
  USE mo_scan_buffer, ONLY: alps, d, dalpsl, dalpsm,  &
                            dtl, dtm, du0, t, u0, &
                            u, v, vo
  USE mo_memory_gl,   ONLY: q, xl, xi, xt
  USE mo_memory_g1a,  ONLY: alpsm1
  USE mo_gaussgrid,   ONLY: gl_sqcst
  !
  ! Prepare or get forcing fields for column model.
  ! To be called from scan1sl after the inverse fft 
  ! (after CALL fourier_to_gridpoint).
  !    
    CHARACTER(len=5) :: name = 'XTnnn'
    INTEGER          :: i, jglat
    REAL(dp)             :: fak
    !
    ! fetch from other PEs
    !
    IF (comode/=' ') THEN
      !
      ! scale to real winds
      !
      IF(ldc%col_1d .AND. uvfac ==0) THEN
        jglat = ldc%glat (lrow(1))              ! global index north -> south
        fak  = 1._dp/gl_sqcst(jglat)
        u(1,:,1) = fak * u(1,:,1)
        v(1,:,1) = fak * v(1,:,1)
      ENDIF
      !
      ! flush files
      !
      IF (p_pe == p_io) THEN
        IF (lfull)            CALL flush (unit_f1d)
        IF (p_col/=-1) THEN
          IF(comode=='resid') CALL flush (unit_r1d)
          CALL flush (unit_b1d)
          CALL flush (unit_a1d)
        ENDIF
      ENDIF
      !
      ! rewind forcing, residui files
      !
      IF (lrewind .AND. p_pe == p_io .AND. p_col/=-1 .AND. &
         (comode=='force' .OR. comode=='add')) THEN
        REWIND (unit_f1d)
        REWIND (unit_r1d)
      ENDIF
      !
      ! write nstep, latitude/longitude index
      !
      CALL process_nstep_lonlat
      !
      ! prognostic variables in column model
      !
      IF (for_ps >0) CALL process_col2 (alps,'ALPS' ,for_ps ==1,cnudgp )
      IF (for_ps >0) CALL process_col2 (alps,'PSREF',.FALSE.,diagn=.TRUE.)
      IF (for_pst>0) pst = 0._dp
      IF (for_pst>0) CALL process_col2 (  pst, 'PS+' ,for_pst==1)
      IF (for_t  >0) CALL process_col3 (t, 'T'   ,for_t  ==1,cnudgt)
      IF (k_vg>0 .AND. ldc%col_1d) THEN
        IF (for_uv >0) CALL process_col3 ( u, 'U',for_uv ==1,cnudguv,fo=ug)
        IF (for_uv >0) CALL process_col3 ( v, 'V',for_uv ==1,cnudguv,fo=vg)
        ug (:, k_vg +1:,:) = SPREAD (ug (:, k_vg ,:), 2, ldc%nlev-k_vg )
        vg (:, k_vg +1:,:) = SPREAD (vg (:, k_vg ,:), 2, ldc%nlev-k_vg )
      ELSE
        IF (for_uv >0) CALL process_col3 ( u, 'U',for_uv ==1,cnudguv)
        IF (for_uv >0) CALL process_col3 ( v, 'V',for_uv ==1,cnudguv)
      ENDIF
      IF (for_q  >0) CALL process_col3 (  q     , 'Q',for_q  ==1,cnudgq )
      IF (for_x  >0) CALL process_col3 (  xl    ,'XL',for_x  ==1,cnudgx )
      IF (for_x  >0) CALL process_col3 (  xi    ,'XI',for_x  ==1,cnudgx )
      IF (for_xt >0) THEN
        DO i=1,SIZE(xt,3)
          WRITE (name(3:5),'(i3.3)') i
          CALL process_col3 (xt(:,:,i,:), name ,for_xt==1, cnudgxt)
        END DO
      ENDIF
      !
      ! safeguard for negative moisture
      !
      IF (chk_q_0 == 1) THEN
        IF (p_pe == p_io .AND. ANY(q <0._dp)) &
          WRITE(0,*) ' WARNING !!! Q  < 0 :',MINLOC(q ),MINVAL(q )
        q =MAX(q ,0.0_dp)
      ENDIF
      !
      ! surface pressure tendency
      !
      IF (for_pst>0) THEN
        WHERE(pst> 10000._dp)
          pst = (pst - EXP(alpsm1)) / time_step_len
        ENDWHERE
      ENDIF
      !
      ! rescale winds
      !
      IF(ldc%col_1d .AND. uvfac ==0) THEN
        jglat = ldc% glat (lrow(1))              ! global index north -> south
        fak   = gl_sqcst(jglat)
        u(1,:,1) = fak * u(1,:,1)
        v(1,:,1) = fak * v(1,:,1)
      ENDIF
      !
      ! diagnostic variables in column model
      !
      IF (for_dt>0) CALL process_col3 (   dtm, 'DTM'    ,for_dt==1)
      IF (for_dt>0) CALL process_col3 (   dtl, 'DTL'    ,for_dt==1)
      IF (for_vo>0) CALL process_col3 (    vo, 'VO'     ,for_vo==1)
      IF (for_d >0) CALL process_col3 (     d, 'D'      ,for_d ==1)
      IF (for_dp>0) CALL process_col2 (dalpsl, 'DALPSL' ,for_dp==1)
      IF (for_dp>0) CALL process_col2 (dalpsm, 'DALPSM' ,for_dp==1)
      !
      ! zonal means
      !
      IF (for_zm>0) CALL process_colzm(  u0, 'U0'     ,for_zm==1)
!     IF (for_zm>0) CALL process_colzm(  ul, 'UL'     ,for_zm==1)
      IF (for_zm>0) CALL process_colzm( du0, 'DU0'    ,for_zm==1)
      !
      ! SST
      !
      IF (for_sst>0) CALL process_col2(  sst_1d, 'SST'    ,for_sst==1)
    ENDIF
  END SUBROUTINE get_col_ffti
!------------------------------------------------------------------------------
  SUBROUTINE get_col_pt (zsdiv, d, zvgrps, delp, aph, jrow)
! USE mo_hyb, ONLY:         delb
  REAL(dp)    ,INTENT(inout) :: zsdiv (:,:)   ! pressure tendency
  REAL(dp)    ,INTENT(inout) :: d     (:,:,:) ! divergence
  REAL(dp)    ,INTENT(in)    :: zvgrps(:,:)   ! v * grad(ps)
  REAL(dp)    ,INTENT(in)    :: delp  (:,:)   ! pressure difference across layers
  REAL(dp)    ,INTENT(in)    :: aph   (:,:)   ! pressure
  INTEGER ,INTENT(in)    :: jrow          ! local row index
    !
    ! correct pressure tendencies in dyn (within latitude loop)
    !
    INTEGER :: jl, jk, ikp
!   REAL(dp)    :: zdelb
    REAL(dp)    :: zpsc (ldc% nglon)
    !
    !-- prescribe vertical velocity (Pa/s) in SCM run
    !   this procedure is exact only for zero surface pressure gradient
    !
    IF (comode/=' ') THEN
      IF(for_pt > 0) THEN
        pt_1d(:,:,jrow) = 0._dp
        IF (ldc%col_1d) THEN
          !
          ! SCM run on this PE, routine is called once
          !
          CALL process_col3 (pt_1d,'PT',for_pt==1)
          DO jk = 1, ldc% nlev
            ikp = jk + 1
!           zdelb = delb(jk)
            DO jl = 1, ldc% nglon
              zsdiv(jl,ikp) = pt_1d(jl,ikp,jrow)
              d(jl,jk,jrow) = (zsdiv(jl,ikp) - zsdiv(jl,jk))/delp(jl,jk)
!                               - zdelb*zvgrps(jl,jk)
            END DO
          END DO
        ELSE
          !
          ! 3D model run on this PE, routine is called for each latitude
          ! gather field:
          !
          IF(jrow == 1) CALL process_col3 (pt_1d,'PT',for_pt==1)
        END IF
      END IF
      !
      !-- prescribe surface pressure tendency in SCM run
      !   this procedure is exact only for zero surface pressure gradient
      !
      IF (ldc%col_1d .AND. for_pst == 1) THEN
        zpsc(:) = (pst(:,jrow) - zsdiv(:,ldc% nlev+1)) / aph(:,ldc% nlev+1)
        DO jk = 1, ldc% nlev
          ikp = jk + 1
!         zdelb = delb(jk)
          DO jl = 1, ldc% nglon
            zsdiv(jl,ikp) = zsdiv(jl,ikp) + aph(jl,ikp) * zpsc (jl)
            d(jl,jk,jrow) = (zsdiv(jl,ikp) - zsdiv(jl,jk))/delp(jl,jk)
!                             - zdelb*zvgrps(jl,jk)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE get_col_pt
!------------------------------------------------------------------------------
  SUBROUTINE get_col_pol (tte, qte, jlat)
  REAL(dp)    ,INTENT(inout) :: tte(:,:) ! temperature tendency (within row loop)
  REAL(dp)    ,INTENT(inout) :: qte(:,:) ! moisture    tendency (within row loop)
  INTEGER ,INTENT(in)    :: jlat     ! local row index
    !
    ! correct for polefilter in physc (within latitude loop)
    !
    IF (comode/=' ') THEN
      IF(for_pol > 0) THEN
        tte_3d(:,:,jlat) = tte; qte_3d(:,:,jlat) = qte
        IF (ldc% col_1d) THEN
          !
          ! SCM run on this PE, routine is called once
          !
          CALL process_col3 (tte_3d,'TTE_POL',for_pol==1)
          CALL process_col3 (qte_3d,'QTE_POL',for_pol==1)
          tte = tte_3d(:,:,jlat); qte = qte_3d(:,:,jlat)
        ELSE
          !
          ! 3D model run on this PE, routine is called for each latitude
          ! gather field:
          !
          !
          ! send field on last latitude cycle
          !
          IF((comode=='traject'.AND.jlat == ldc% nglat  ).OR.  &
             (comode/='traject'.AND.jlat == MAXVAL(lats))) THEN
            CALL process_col3 (tte_3d,'TTE_POL',for_pol==1)
            CALL process_col3 (qte_3d,'QTE_POL',for_pol==1)
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE get_col_pol
!------------------------------------------------------------------------------
  SUBROUTINE cal_col_expl (ztodt, jlat)
  !
  ! Perform explicit time integration in a column model run 
  ! for variables treated implicitly in the full model.
  ! To be called from scan1sl before the call to 'si1'
  !
  USE mo_scan_buffer, ONLY: u, v, t, alps,   &! new values
                            d, vo, dtm, dtl, &
                            dalpsl, dalpsm,          &
                            vol, vom, tte, alpste
  USE mo_memory_g2a,  ONLY: um1, vm1                          ! old values
  USE mo_memory_g1a,  ONLY: alpsm1, tm1                       ! old values
  USE mo_gaussgrid,   ONLY: gl_coriol, gl_sqcst               ! Coriolis param.
  REAL(dp)    ,INTENT(in) :: ztodt                                ! time step
  INTEGER ,INTENT(in) :: jlat                                 ! latit. index
    INTEGER :: jglat       ! globel indices: north -> south
    INTEGER :: ngl         ! number of gaussian latitudes
    REAL(dp)    :: zcorio      ! coriolis parameter
    REAL(dp)    :: zcst        ! scale faktor for winds outside 'physc'
    IF(ldc% col_1d) THEN
      !
      ! get coriolis parameter
      !
      ngl   = ldc% nlat                      ! number of gaussian latitudes
      jglat = ldc% glat(jlat)                ! global continuous north -> south
      zcorio = gl_coriol(jglat)
      zcst   = gl_sqcst (jglat)
      !
      ! explicit integration for prognostic variables
      !
      IF (ug_1d/=0._dp .OR. vg_1d/=0._dp .OR. k_vg/=0) THEN
!       vom = vom + zcorio * (v(:,:,jlat) - zcst*vg(:,:,jlat))
!       vol = vol - zcorio * (u(:,:,jlat) - zcst*ug(:,:,jlat))
        vom(:,:,jlat) = vom(:,:,jlat) + zcorio * (- zcst*vg(:,:,jlat))
        vol(:,:,jlat) = vol(:,:,jlat) - zcorio * (- zcst*ug(:,:,jlat))
      ENDIF
      IF (int_uv==1) u    (:,:,jlat) = um1   (:,:,jlat) + ztodt * vom   (:,:,jlat)
      IF (int_uv==1) v    (:,:,jlat) = vm1   (:,:,jlat) + ztodt * vol   (:,:,jlat)
      IF (int_t ==1) t    (:,:,jlat) = tm1   (:,:,jlat) + ztodt * tte   (:,:,jlat)
      IF (int_ps==1) alps (:,  jlat) = alpsm1(:,  jlat) + ztodt * alpste(:,  jlat)
      !
      ! zero prescribed terms if flag == 0
      !
      IF (for_d ==0) d      (:,:,jlat) = 0._dp
      IF (for_vo==0) vo     (:,:,jlat) = 0._dp
      IF (for_dt==0) dtm    (:,:,jlat) = 0._dp
      IF (for_dt==0) dtl    (:,:,jlat) = 0._dp
      IF (for_dp==0) dalpsl (:,  jlat) = 0._dp
      IF (for_dp==0) dalpsm (:,  jlat) = 0._dp
      nforce = nforce - 1
    ENDIF
  END SUBROUTINE cal_col_expl
!==============================================================================
! Routines to be called for diagnostics
!==============================================================================
  SUBROUTINE write_column3 (x, name)
  REAL(dp)             ,INTENT(inout) :: x (:,:,:)
  CHARACTER(len=*) ,INTENT(in)    :: name
      CALL process_col3 (x, name, .FALSE.)
  END SUBROUTINE write_column3
!------------------------------------------------------------------------------
  SUBROUTINE write_column2 (x, name)
  REAL(dp)             ,INTENT(inout) :: x (:,:)
  CHARACTER(len=*) ,INTENT(in)    :: name
      CALL process_col2 (x, name, .FALSE.)
  END SUBROUTINE write_column2
!------------------------------------------------------------------------------
  SUBROUTINE write_column1 (x, name, jlat)
  REAL(dp)             ,INTENT(in) :: x (:)
  CHARACTER(len=*) ,INTENT(in) :: name
  INTEGER          ,INTENT(in) :: jlat
    REAL(dp) :: tmp (ldc%nglon,SIZE(x),ldc%nglat)
    IF (ncol==1) THEN
      IF (ANY (jlat==lrow)) THEN
        tmp = 0._dp
        tmp (1,:,jlat) = x
        CALL process_col3 (tmp, name, .FALSE.)
      ENDIF
    ENDIF
  END SUBROUTINE write_column1
!------------------------------------------------------------------------------
  SUBROUTINE write_column3r (x, name, jlat)
  REAL(dp)             ,INTENT(in) :: x (:,:)
  CHARACTER(len=*) ,INTENT(in) :: name
  INTEGER          ,INTENT(in) :: jlat
    REAL(dp) :: tmp (SIZE(x,1),SIZE(x,2),ldc%nglat)
    IF (ncol==1) THEN
      IF (ANY (jlat==lrow)) THEN
        tmp = 0._dp
        tmp (:,:,jlat) = x
        CALL process_col3 (tmp, name, .FALSE.)
      ENDIF
    ENDIF
  END SUBROUTINE write_column3r
!------------------------------------------------------------------------------
  SUBROUTINE write_column2r (x, name, jlat)
  REAL(dp)             ,INTENT(in) :: x (:)
  CHARACTER(len=*) ,INTENT(in) :: name
  INTEGER          ,INTENT(in) :: jlat
    REAL(dp) :: tmp (SIZE(x),ldc%nglat)
    IF (ncol==1) THEN
      IF (ANY (jlat==lrow)) THEN
        tmp = 0._dp
        tmp (:,jlat) = x
        CALL process_col2 (tmp, name, .FALSE.)
      ENDIF
    ENDIF
  END SUBROUTINE write_column2r
!==============================================================================
! Private routines within this module
!==============================================================================
!------------------------------------------------------------------------------
  SUBROUTINE gather_col3 (co, gp, p_dest)
  REAL(dp)    ,INTENT(in)           :: gp (:,:,:) ! local field (lon ,lev,lat)
  REAL(dp)    ,POINTER              :: co (:,:,:) ! column      (ncol,lev, 1)
  INTEGER ,INTENT(in)           :: p_dest     ! destination PE
  !
  ! get a column from the full model (distributed domain)
  !
    REAL(dp)    :: tmp (ncol,SIZE(gp,2))
    INTEGER :: i, j, pe

    IF (p_pe == p_dest) THEN
      IF(SIZE(co,1)/=ncol.OR.SIZE(co,3)/=1) THEN
        WRITE(nerr,*)'gather_col3: p_pe, p_dest, lfull=',p_pe, p_dest, lfull
        CALL finish ('gather_col3: shape(co)/=(/ncol,:,1/)')
      ENDIF

!      IF (lfull /= p_dest) THEN
!        CALL p_recv (co, lfull, tag_gather_gp)
!      ELSE
!        DO i=1,nlcol
!          co(i,:,1) = gp(lcol(i),:,lrow(i))
!        END DO
!      ENDIF 
!    ELSE IF (p_pe == lfull) THEN
!      DO i=1,nlcol
!        tmp(i,:) = gp(lcol(i),:,lrow(i))
!      END DO
!      CALL p_send (tmp, p_dest, tag_gather_gp)
!    ENDIF

      DO pe = 0, p_nprocs-1
        IF (ANY(pes(:ncol)==pe)) THEN
          IF(pe==p_dest) THEN
            DO i=1,nlcol
              tmp(i,:) = gp(lcol(i),:,lrow(i))
            END DO
          ELSE
            CALL p_recv (tmp, pe, tag_gather_gp)
          ENDIF
          j=0
          DO i=1,ncol
            IF (pes(i)==pe) THEN
              j=j+1
              co(i,:,1) = tmp(j,:)
            ENDIF
          END DO
        ENDIF
      END DO
    ELSE IF (nlcol>0) THEN
      DO i=1,nlcol
        tmp(i,:) = gp(lcol(i),:,lrow(i))
      END DO
      CALL p_send (tmp, p_dest, tag_gather_gp)
    ENDIF

  END SUBROUTINE gather_col3
!------------------------------------------------------------------------------
  SUBROUTINE gather_col2 (co, gp, p_dest)
  REAL(dp)    ,INTENT(in)           :: gp (:,:) ! local field (lon ,lat)
  REAL(dp)    ,POINTER              :: co (:,:) ! column      (ncol, 1 )
  INTEGER ,INTENT(in)           :: p_dest   ! destination PE
  !
  ! get surface field from the full model (distributed domain)
  !
    REAL(dp)    :: tmp (ncol)
    INTEGER :: i, j, pe
    IF (p_pe == p_dest) THEN
      IF(SIZE(co,1)/=ncol.OR.SIZE(co,2)/=1) &
        CALL finish ('gather_col2: shape(co)/=(/ncol,1/)')
      DO pe = 0, p_nprocs-1
        IF (ANY(pes(:ncol)==pe)) THEN
          IF(pe==p_dest) THEN
            DO i=1,nlcol
              tmp(i) = gp(lcol(i),lrow(i))
            END DO
          ELSE
            CALL p_recv (tmp, pe, tag_gather_gp)
          ENDIF
          j=0
          DO i=1,ncol
            IF (pes(i)==pe) THEN
              j=j+1
              co(i,1) = tmp(j)
            ENDIF
          END DO
        ENDIF
      END DO
    ELSE IF (nlcol>0) THEN
      DO i=1,nlcol
        tmp(i) = gp(lcol(i),lrow(i))
      END DO
      CALL p_send (tmp, p_dest, tag_gather_gp)
    ENDIF
  END SUBROUTINE gather_col2
!------------------------------------------------------------------------------
  SUBROUTINE gather_colzm (co, gp, p_dest)
  REAL(dp)    ,INTENT(in)           :: gp (:,:) ! local field (lev ,lat)
  REAL(dp)    ,POINTER              :: co (:,:) ! column      (ncol,lev)
  INTEGER ,INTENT(in)           :: p_dest   ! destination PE
  !
  ! get zonal mean from the full model (distributed domain)
  !
    REAL(dp)    :: tmp (ncol,SIZE(gp,1))
    INTEGER :: i, j, pe
    IF (p_pe == p_dest) THEN
      IF(SIZE(co,1)/=ncol) CALL finish ('gather_col2: size(co,1)/=ncol')
      DO pe = 0, p_nprocs-1
        IF (ANY(pes(:ncol)==pe)) THEN
          IF(pe==p_dest) THEN
            DO i=1,nlcol
              tmp(i,:) = gp(:,lrow(i))
            END DO
          ELSE
            CALL p_recv (tmp, pe, tag_gather_gp)
          ENDIF
          j=0
          DO i=1,ncol
            IF (pes(i)==pe) THEN
              j=j+1
              co(i,:) = tmp(j,:)
            ENDIF
          END DO
        ENDIF
      END DO
    ELSE IF (nlcol>0) THEN
      DO i=1,nlcol
        tmp(i,:) = gp(:,lrow(i))
      END DO
      CALL p_send (tmp, p_dest, tag_gather_gp)
    ENDIF
  END SUBROUTINE gather_colzm
!------------------------------------------------------------------------------
  SUBROUTINE process_nstep_lonlat
  USE mo_gaussgrid,   ONLY: gl_sqcst
    !
    ! write time step, longitude/latitude index, cvt table
    !
    REAL(dp) ,ALLOCATABLE :: ijlonlat (:,:) ! output buffer for lon/lat indices
    REAL(dp) ,ALLOCATABLE :: facuv    (:)   ! factor for u,v
    REAL(dp)              :: rstep          ! output buffer for time step
    CHARACTER(len=16) :: nam16          ! output buffer for field identifier
    INTEGER           :: iu(4)          ! output units to process
    INTEGER           :: i              ! do loop index
    INTEGER           :: jglat          ! latitude index
    INTEGER           :: nstep

    IF (p_pe==p_io) THEN
      !
      ! determine output units
      !
      iu = 0
      SELECT CASE (comode)
      CASE ('')
      CASE ('traject','test')
        iu(1) = unit_f1d
      CASE ('resid')
        IF (p_col/=-1) THEN
          iu(2) = unit_r1d
          iu(3) = unit_b1d
          iu(4) = unit_a1d
        ENDIF
      CASE ('force','add','free')
        IF (p_col/=-1) THEN
          iu(3) = unit_b1d
          iu(4) = unit_a1d
        ENDIF
      CASE default
      END SELECT
      !
      ! write
      !
      ALLOCATE (ijlonlat (ncol,2))
      ALLOCATE (facuv    (ncol  ))
      ijlonlat (:,1) = lons
      ijlonlat (:,2) = lats
      nstep = get_time_step()
      rstep = nstep
      IF(uvfac==0) THEN
        facuv=1._dp
      ELSE
        DO i=1,ncol
          jglat = lats(i)                           ! global index n->s
          facuv(i)= gl_sqcst(jglat)
        END DO
      END IF
      DO i = 1, SIZE(iu)
        IF (iu(i)==0) CYCLE
        nam16 = 'NSTEP'
        WRITE (iu(i)) nam16, 1, 1, 0, 0
        WRITE (iu(i)) rstep
        IF (.NOT.(lresume.OR.lstart)) CYCLE
        nam16 = 'ILON_JLAT' 
        WRITE (iu(i)) nam16, 2, SHAPE(ijlonlat), 0
        WRITE (iu(i)) ijlonlat
        nam16 = 'DTIME'
        WRITE (iu(i)) nam16, 2, 1, 1, 0
        WRITE (iu(i)) delta_time
        nam16 = 'AK'
        WRITE (iu(i)) nam16, 2, 1, nlevp1, 0
        WRITE (iu(i)) vct(1:nlevp1)
        nam16 = 'BK'
        WRITE (iu(i)) nam16, 2, 1, nlevp1, 0
        WRITE (iu(i)) vct(nvclev+1:nvclev+nlevp1)
        nam16 = 'FACUV'
        WRITE (iu(i)) nam16, 2, ncol, 1, 0
        WRITE (iu(i)) facuv
      END DO
      DEALLOCATE (ijlonlat,facuv)
    ENDIF
  END SUBROUTINE process_nstep_lonlat
!------------------------------------------------------------------------------
  SUBROUTINE process_col3 (gp, name, lforce, wnudge, fo)
  REAL(dp)              ,INTENT(inout) :: gp (:,:,:)
  CHARACTER (len=*) ,INTENT(in)    :: name
  LOGICAL           ,INTENT(in)    :: lforce
  REAL(dp)    ,OPTIONAL ,INTENT(in)    :: wnudge(:)
  REAL(dp)    ,OPTIONAL ,INTENT(out)   :: fo (:,:,:)
    REAL(dp), POINTER      :: f (:,:,:)
    REAL(dp), POINTER      :: r (:,:,:)
    REAL(dp), POINTER      :: m (:,:,:)
    REAL(dp), POINTER      :: b (:,:,:)
    INTEGER            :: k
    CHARACTER (len=16) :: nam16
    nam16 = name
    !
    ! allocate temporaries
    !
    NULLIFY (f,r,m)
    IF (p_pe==p_io) THEN
      ALLOCATE (f(ncol,SIZE(gp,2),1))
      ALLOCATE (r(ncol,SIZE(gp,2),1))
      ALLOCATE (m(ncol,SIZE(gp,2),1))
      ALLOCATE (b(ncol,SIZE(gp,2),1))
      f = 0._dp
      r = 0._dp
    ENDIF
    !
    ! get forcing if full model is running
    !
    IF (lfull) CALL gather_col3 (f,gp,p_io)
    !
    ! read forcing
    !
    IF (lforce) THEN
      IF(p_pe==p_io) THEN
        SELECT CASE (comode)
        CASE ('force','resid')
          CALL read_col3 (f, name, unit_f1d)
        CASE ('add')
          CALL read_col3 (f, name, unit_f1d)
          CALL read_col3 (r, name, unit_r1d)
        CASE ('free')
          IF (nforce > 0) CALL read_col3 (f, name, unit_f1d)
        END SELECT
      ENDIF
    ENDIF
    !
    ! get column model result if column model is running
    !
    IF (p_col>-1.AND.(p_pe==p_io.OR.p_pe==p_col)) THEN
      IF(p_pe==p_io) THEN
        IF(p_io==p_col) THEN
          m = gp
        ELSE
          CALL p_recv (m ,p_col, tag_gather_gp)
        ENDIF
        b = m
      ENDIF
      IF(p_pe==p_col .AND. p_io/=p_col) CALL p_send (gp, p_io, tag_gather_gp)
    ENDIF
    !
    ! process forcing
    !
    fapp = .FALSE.
    IF (lforce) THEN
      SELECT CASE (comode)
      CASE ('resid')
        IF(p_pe==p_io) THEN
          r = f - m
!         r = r + (f - (m + r))
          m = m + r
          IF(lfres) m = f
        ENDIF
        fapp = .TRUE.
      CASE ('free')
      CASE ('add')
        IF(p_pe==p_io) THEN
          m = m + r
          IF(lcnudge .AND. PRESENT(wnudge)) THEN
            DO k=1,SIZE(m,2)
              m(:,k,:) = wnudge(k) * f(:,k,:) + (1._dp-wnudge(k)) * m(:,k,:)
            END DO
          END IF
        END IF
      CASE default
        IF(p_pe==p_io) m = f
        fapp = .TRUE.
      END SELECT
      IF (nforce > 0.AND. p_pe==p_io) m = f
      IF (nforce > 0) fapp = .TRUE.
    ENDIF
    !
    ! write forcing, residui and results
    !
    IF(p_pe==p_io) THEN
      SELECT CASE (comode)
      CASE ('traject','test')
        WRITE (unit_f1d) nam16, 3, SHAPE(f)
        WRITE (unit_f1d) f
      END SELECT
      SELECT CASE (comode)
      CASE ('resid','test')
        WRITE (unit_r1d) nam16, 3, SHAPE(r)
        WRITE (unit_r1d) r
      END SELECT
      IF (p_col>-1) THEN
        WRITE (unit_b1d) nam16, 3, SHAPE(r)
        WRITE (unit_b1d) b
        WRITE (unit_a1d) nam16, 3, SHAPE(r)
        WRITE (unit_a1d) m
      ENDIF
    ENDIF
    !
    ! apply forcing
    !
    IF (lforce) THEN
      SELECT CASE (comode)
      CASE ('add','traject','resid','force','free')
        IF(p_pe==p_io) THEN
          IF(p_io==p_col) THEN
            gp = m
            IF(PRESENT(fo)) fo = f 
          ELSE
            CALL p_send (m ,p_col, tag_gather_gp)
            IF(PRESENT(fo)) CALL p_send (f ,p_col, tag_gather_gp)
          ENDIF
        ENDIF
        IF(p_pe==p_col .AND. p_io/=p_col) THEN
          CALL p_recv (gp, p_io, tag_gather_gp)
          IF(PRESENT(fo)) CALL p_recv (fo, p_io, tag_gather_gp)
        END IF
      END SELECT
    ENDIF
    !
    ! deallocate temporaries
    !
    IF (p_pe==p_io)  DEALLOCATE (f,r,m,b)
  END SUBROUTINE process_col3
!------------------------------------------------------------------------------
  SUBROUTINE process_col2 (gp, name, lforce, wnudge, diagn)
  REAL(dp)              ,INTENT(inout) :: gp (:,:)
  CHARACTER (len=*) ,INTENT(in)    :: name
  LOGICAL           ,INTENT(in)    :: lforce
  REAL(dp)    ,OPTIONAL ,INTENT(in)    :: wnudge
  LOGICAL ,OPTIONAL ,INTENT(in)    :: diagn
    REAL(dp), POINTER      :: f (:,:)
    REAL(dp), POINTER      :: r (:,:)
    REAL(dp), POINTER      :: m (:,:)
    REAL(dp), POINTER      :: b (:,:)
    CHARACTER (len=16) :: nam16
    LOGICAL            :: ldiagn
    nam16  = name
    ldiagn = .FALSE.; IF (PRESENT(diagn)) ldiagn = diagn
    !
    ! allocate temporaries
    !
    NULLIFY (f,r,m)
    IF (p_pe==p_io) THEN
      ALLOCATE (f(ncol,1))
      ALLOCATE (r(ncol,1))
      ALLOCATE (m(ncol,1))
      ALLOCATE (b(ncol,1))
      f = 0._dp
      r = 0._dp
    ENDIF
    !
    ! get forcing if full model is running
    !
    IF (lfull) CALL gather_col2 (f,gp,p_io)
    !
    ! read forcing
    !
    IF (lforce) THEN
      IF(p_pe==p_io) THEN
        SELECT CASE (comode)
        CASE ('force','resid')
          CALL read_col2 (f, name, unit_f1d)
        CASE ('add')
          CALL read_col2 (f, name, unit_f1d)
          CALL read_col2 (r, name, unit_r1d)
        CASE ('free')
          IF (nforce > 0) CALL read_col2 (f, name, unit_f1d)
        END SELECT
      ENDIF
    ENDIF
    !
    ! get column model result if column model is running
    !
    IF (p_col>-1.AND.(p_pe==p_io.OR.p_pe==p_col)) THEN
      IF(p_pe==p_io) THEN
        IF(p_io==p_col) THEN
          m = gp
        ELSE
          CALL p_recv (m ,p_col, tag_gather_gp)
        ENDIF
        b = m
      ENDIF
      IF(p_pe==p_col .AND. p_io/=p_col) CALL p_send (gp, p_io, tag_gather_gp)
    ENDIF
    !
    ! process forcing
    !
    IF (lforce) THEN
      IF(p_pe==p_io) THEN
        SELECT CASE (comode)
        CASE ('resid')
          r = f - m
!         r = r + (f - (m + r))
          m = m + r
          IF(lfres) m = f
        CASE ('free')
        CASE ('add')
          m = m + r
          IF(lcnudge .AND. PRESENT(wnudge)) THEN
            m(:,:) = wnudge * f(:,:) + (1._dp-wnudge) * m(:,:)
          END IF
        CASE default
          m = f
        END SELECT
        IF (nforce > 0) m = f
      ENDIF
    ENDIF
    !
    ! diagnostic output
    !
    IF (ldiagn .AND. p_pe==p_io .AND. comode=='resid') r = m
    !
    ! write forcing
    !
    IF(p_pe==p_io) THEN
      SELECT CASE (comode)
      CASE ('traject','test')
        WRITE (unit_f1d) nam16, 2, SHAPE(f), 0
        WRITE (unit_f1d) f
      END SELECT
      SELECT CASE (comode)
      CASE ('resid','test')
        WRITE (unit_r1d) nam16, 2, SHAPE(r), 0
        WRITE (unit_r1d) r
      END SELECT
      IF (p_col>-1) THEN
        WRITE (unit_a1d) nam16, 2, SHAPE(r), 0
        WRITE (unit_a1d) m
        WRITE (unit_b1d) nam16, 2, SHAPE(r), 0
        WRITE (unit_b1d) b
      ENDIF
    ENDIF
    !
    ! apply forcing
    !
    IF (lforce) THEN
      SELECT CASE (comode)
      CASE ('add','traject','resid','force','free')
        IF(p_pe==p_io) THEN
          IF(p_io==p_col) THEN
            gp = m
          ELSE
            CALL p_send (m ,p_col, tag_gather_gp)
          ENDIF
        ENDIF
        IF(p_pe==p_col .AND. p_io/=p_col) CALL p_recv (gp, p_io, tag_gather_gp)
      END SELECT
    ENDIF
    !
    ! deallocate temporaries
    !
    IF (p_pe==p_io)  DEALLOCATE (f,r,m,b)
  END SUBROUTINE process_col2
!------------------------------------------------------------------------------
  SUBROUTINE process_colzm (gp, name, lforce)
  REAL(dp)              ,INTENT(inout) :: gp (:,:) ! zonal mean (lev,lat)
  CHARACTER (len=*) ,INTENT(in)    :: name     ! name of field
  LOGICAL           ,INTENT(in)    :: lforce   ! modify field ?
    REAL(dp), POINTER      :: f (:,:)
    REAL(dp), POINTER      :: r (:,:)
    REAL(dp), POINTER      :: m (:,:)
    REAL(dp), POINTER      :: b (:,:)
    CHARACTER (len=16) :: nam16
    nam16 = name
    !
    ! allocate temporaries
    !
    NULLIFY (f,r,m)
    IF (p_pe==p_io) THEN
      ALLOCATE (f(ncol,SIZE(gp,1))) ! (ncol,lev)
      ALLOCATE (r(ncol,SIZE(gp,1))) !
      ALLOCATE (m(ncol,SIZE(gp,1))) !
      ALLOCATE (b(ncol,SIZE(gp,1))) !
      f = 0._dp
      r = 0._dp
    ENDIF
    !
    ! get forcing if full model is running
    !
    IF (lfull) CALL gather_colzm(f,gp,p_io)
    !
    ! read forcing
    !
    IF (lforce) THEN
      IF(p_pe==p_io) THEN
        SELECT CASE (comode)
        CASE ('force','resid')
          CALL read_col2 (f, name, unit_f1d)
        CASE ('add')
          CALL read_col2 (f, name, unit_f1d)
          CALL read_col2 (r, name, unit_r1d)
        CASE ('free')
          IF (nforce > 0) CALL read_col2 (f, name, unit_f1d)
        END SELECT
      ENDIF
    ENDIF
    !
    ! get column model result if column model is running
    !
    IF (p_col>-1.AND.(p_pe==p_io.OR.p_pe==p_col)) THEN
      IF(p_pe==p_io) THEN
        IF(p_io==p_col) THEN
          m(1,:) = gp (:,1)
        ELSE
          CALL p_recv (m ,p_col, tag_gather_gp)
        ENDIF
        b = m
      ENDIF
      IF(p_pe==p_col .AND. p_io/=p_col) CALL p_send (gp, p_io, tag_gather_gp)
    ENDIF
    !
    ! process forcing
    !
    IF (lforce) THEN
      IF(p_pe==p_io) THEN
        SELECT CASE (comode)
        CASE ('resid')
          r = f - m
!         r = r + (f - (m + r))
          m = m + r
          IF(lfres) m = f
        CASE ('free')
        CASE ('add')
          m = m + r
        CASE default
          m = f
        END SELECT
        IF (nforce > 0) m = f
      ENDIF
    ENDIF
    !
    ! write forcing
    !
    IF(p_pe==p_io) THEN
      SELECT CASE (comode)
      CASE ('traject','test')
        WRITE (unit_f1d) nam16, 2, SHAPE(f), 0
        WRITE (unit_f1d) f
      END SELECT
      SELECT CASE (comode)
      CASE ('resid','test')
        WRITE (unit_r1d) nam16, 2, SHAPE(r), 0
        WRITE (unit_r1d) r
      END SELECT
      IF (p_col>-1) THEN
        WRITE (unit_a1d) nam16, 2, SHAPE(r), 0
        WRITE (unit_a1d) m
        WRITE (unit_b1d) nam16, 2, SHAPE(r), 0
        WRITE (unit_b1d) b
      ENDIF
    ENDIF
    !
    ! apply forcing
    !
    IF (lforce) THEN
      SELECT CASE (comode)
      CASE ('add','traject','resid','force','free')
        IF(p_pe==p_io) THEN
          IF(p_io==p_col) THEN
            gp(:,1) = m (1,:)
          ELSE
            CALL p_send (m ,p_col, tag_gather_gp)
          ENDIF
        ENDIF
        IF(p_pe==p_col .AND. p_io/=p_col) CALL p_recv (gp, p_io, tag_gather_gp)
      END SELECT
    ENDIF
    !
    ! deallocate temporaries
    !
    IF (p_pe==p_io)  DEALLOCATE (f,r,m,b)
  END SUBROUTINE process_colzm
!------------------------------------------------------------------------------
  SUBROUTINE read_col3 (x, name, iunit)
  REAL(dp)             ,INTENT(out) :: x(:,:,:)
  CHARACTER(len=*) ,INTENT(in)  :: name 
  INTEGER          ,INTENT(in)  :: iunit
    CHARACTER(len=16) :: nam16
    INTEGER           :: rank, shap(3), ios
    REAL(dp) ,ALLOCATABLE :: z(:,:,:)
    INTEGER           :: ireadold
l1: DO
      READ (iunit, iostat=ios) nam16, rank, shap
      IF (ios/=0) EXIT
      IF (nam16=='ILON_JLAT' .AND. iunit==unit_f1d) THEN
        ALLOCATE (z(shap(1),shap(2),1))
        READ (iunit) z
        ireadold = iread
        DO iread = 1,shap(1)
          IF (z(iread,1,1)==lon_1d(1).AND.z(iread,2,1)==lat_1d(1)) THEN
            IF (ireadold /= iread)                           &
              WRITE(nout,*)'read_col2: iread set to ',iread, &
              ',ilon,jlat=',lon_1d(1),lat_1d(1)
            DEALLOCATE (z)
            CYCLE l1
          ENDIF
        END DO
        WRITE(nerr,*)'read_col2: ilon,jlat=',lon_1d(1),lat_1d(1)
        CALL finish('read_col2','lat/lon-index not found in forcing file')
      ENDIF
      ALLOCATE (z(shap(1),shap(2),shap(3)))
      READ (iunit) z
      IF (nam16==name) THEN
        IF (rank /= 3 .OR. ANY (shap(2:3)/=SHAPE(x(1,:,:)))) &
          CALL finish ('read_col3', 'invalid rank/shape of field '//name)
        IF (shap(1)>1) THEN
          x = z (iread:iread,:,:)
        ELSE
          x = z
        ENDIF
        DEALLOCATE (z)
        RETURN
      ENDIF
      DEALLOCATE (z)
    END DO l1
    CALL finish ('read_col3', 'field not found: '//name)
  END SUBROUTINE read_col3
!------------------------------------------------------------------------------
  SUBROUTINE read_col2 (x, name, iunit)
  REAL(dp)             ,INTENT(out) :: x(:,:)
  CHARACTER(len=*) ,INTENT(in)  :: name
  INTEGER          ,INTENT(in)  :: iunit
    CHARACTER(len=16) :: nam16
    INTEGER           :: rank, shap(3), ios
    REAL(dp) ,ALLOCATABLE :: z(:,:)
l1: DO
      READ (iunit, iostat=ios) nam16, rank, shap
      IF (ios/=0) EXIT
      ALLOCATE (z(shap(1),shap(2)))
      READ (iunit) z
      IF (nam16=='ILON_JLAT' .AND. iunit==unit_f1d) THEN
        DO iread = 1,shap(1)
          IF (z(iread,1)==lon_1d(1).AND.z(iread,2)==lat_1d(1)) THEN
            WRITE(nout,*)'read_col2: iread set to ',iread,&
                                      ',ilon,jlat=',lon_1d(1),lat_1d(1)
            DEALLOCATE (z)
            CYCLE l1
          ENDIF
        END DO
        WRITE(nerr,*)'read_col2: ilon,jlat=',lon_1d(1),lat_1d(1)
        CALL finish('read_col2','lat/lon-index not found in forcing file')
      ENDIF
      IF (nam16==name) THEN
        IF (rank /= 2 .OR. shap(2)/=SIZE(x,2)) THEN
          CALL finish ('read_col2', 'invalid rank/shape of field '//name)
        ENDIF
        IF (shap(1)>1) THEN
          x = z (iread:iread,:)
        ELSE
          x = z
        ENDIF
        DEALLOCATE (z)
        RETURN
      ENDIF
      DEALLOCATE (z)
    END DO l1
    CALL finish ('read_col2', 'field not found: '//name)
  END SUBROUTINE read_col2
!------------------------------------------------------------------------------
  SUBROUTINE read_colx (x, name, iunit)
  REAL(dp)             ,POINTER    :: x(:)
  CHARACTER(len=*) ,INTENT(in) :: name
  INTEGER          ,INTENT(in) :: iunit
    CHARACTER(len=16) :: nam16
    INTEGER           :: rank, shap(3), ios
l1: DO
      IF(ASSOCIATED(x)) DEALLOCATE (x)
      READ (iunit, iostat=ios) nam16, rank, shap
      IF (ios/=0) EXIT
      IF (nam16==name) THEN
        ALLOCATE (x(PRODUCT(shap(1:rank))))
        READ (iunit) x
        RETURN
      ELSE
        READ (iunit)
      ENDIF
    END DO l1
    CALL finish ('read_colx', 'field not found: '//name)      
  END SUBROUTINE read_colx
!------------------------------------------------------------------------------
END MODULE mo_column
