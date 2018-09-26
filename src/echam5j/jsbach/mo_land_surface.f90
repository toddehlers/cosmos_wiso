MODULE mo_land_surface

  USE mo_jsbach,        ONLY: debug
  USE mo_jsbach_grid,   ONLY: grid_type, domain_type, kstart, kend, nidx
  USE mo_jsbach_lctlib, ONLY: lctlib_type

  USE mo_netcdf,        ONLY: FILE_INFO
  USE mo_linked_list,   ONLY: t_stream
  USE mo_mpi,           ONLY: p_parallel, p_parallel_io, p_bcast, p_io, p_pe
  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: message, message_text, finish, int2string

  IMPLICIT NONE

  PUBLIC ::  init_land_surface
  PUBLIC ::  init_albedo
  PUBLIC ::  update_land_surface_fast
  PUBLIC ::  update_albedo
  PUBLIC ::  land_surface_diagnostics
  PUBLIC ::  scale_cover_fract

  TYPE land_surface_type
     INTEGER              :: ntiles_lct
     INTEGER              :: ntiles
     INTEGER, POINTER, DIMENSION(:,:)     :: & !! (nland, ntiles_lct)
          cover_type                           !! Index into LctLibrary
     REAL(dp),    POINTER, DIMENSION(:,:) :: & !! (nland, ntiles)
          cover_type_real,                   & !! cover_type converted to REAL for netCDF files
          cover_fract,                       & !! Fraction of coverage for each land cover type
          albedo,                            & !! Surface albedo
          albedo_vis,                        & !! Surface albedo in the visible range
          albedo_nir,                        & !! Surface albedo in the NIR range
          veg_ratio_actual_per_tile            !! This is the actual fraction of a tile covered by a vegetation canopy, i.e. this value
                                               !! includes the actual leaf area index (i.e. foliar projective cover)

     REAL(dp),    POINTER, DIMENSION(:,:) :: & !! (nland, month)
          veg_ratio
     REAL(dp), POINTER, DIMENSION(:)      :: & !! (nland)
          veg_ratio_max,                     & !! maximal fraction of grid box covered by leaves, assuming infinite LAI
          rock_fract,                        & !! fraction of grid box not suitable for plants
          elevation,                         & !! Mean land grid cell elevation (m)
          oro_std_dev,                       & !! Standard deviation of orography
          forest_fract,                      & !! Forest fraction (TEMPORARY!)
          area,                              & !! Area of land fraction of grid cell (m^2)
          swdown_acc,                        & !! accumulated solar downward radiation
          swdown_reflect_acc                   !! accumulated reflected solar radiation
     LOGICAL, POINTER, DIMENSION(:,:) ::     & !! (nland,  ntiles)
          is_bare_soil,     &     !! Tile is bare soil (but not glacier)?
          is_vegetation,    &     !! Tile is vegetation (natural or crop)?
          is_forest,        &     !! Tile is forest
          is_glacier,       &     !! Tile is glacier?
          is_lake,          &     !! Tile is lake?
          is_present              !! Is tile present, i.e. should it be handled by the land model?
  END TYPE land_surface_type
  PUBLIC ::  land_surface_type

  TYPE land_surface_diag_type
     REAL(dp), POINTER, DIMENSION(:)     :: & !! (nland)
          albedo
  END TYPE land_surface_diag_type

  TYPE albedo_options_type
     LOGICAL :: useAlbedocanopy
  END TYPE albedo_options_type
  TYPE(albedo_options_type), SAVE :: albedo_options

  TYPE albedo_params_type      !vg: Some arrays of the land_surface_type should go in here
     REAL(dp) :: dummy             !vg: just to have an entry to enable compilation 
  END TYPE albedo_params_type
  TYPE(albedo_params_type), SAVE :: albedo_params

  REAL(dp), PARAMETER :: fract_small = 1.e-10_dp       !! very small fraction (e.g. minimum value of cover_fract)
  PUBLIC :: fract_small

  PRIVATE

  INTEGER, SAVE :: nlct = -1        !! Number of land cover types
  INTEGER, SAVE :: ntiles_lct = -1  !! Maximum number of land cover types actually used (<= nclt)
  INTEGER, SAVE :: ntiles = -1      !! Number of tiles (== ntiles_lct if only the tiling structure according to
                                    !! land cover types is used, but can be larger if other sub-tiling structures
                                    !! are used, e.g. snowbands, distributed precipitation.

  TYPE(t_stream), POINTER :: IO_land_surface
  TYPE(t_stream), POINTER :: IO_diag
  TYPE(FILE_INFO), SAVE   :: land_surface_file

  TYPE(land_surface_diag_type), SAVE    :: land_surface_diag

  LOGICAL, SAVE :: module_configured  = .FALSE.
  LOGICAL, SAVE :: module_initialized = .FALSE.

  REAL(dp), POINTER, SAVE :: init_cover_fract(:,:)
  REAL(dp), POINTER, SAVE :: init_veg_ratio_max(:)
  PUBLIC :: init_cover_fract, init_veg_ratio_max

  REAL(dp), PARAMETER :: SkyViewFactor = 1.0_dp        !! Constant in albedo calculation
  REAL(dp), PARAMETER :: AlbedoCanopySnow = 0.25_dp    !! Albedo of snow covered canopy
  REAL(dp), PARAMETER :: AlbedoGlacierVisMin = 0.78_dp !! Albedo of glacier in the visible range at the melting point
  REAL(dp), PARAMETER :: AlbedoGlacierVisMax = 0.9_dp  !! Albedo of glacier in the visible range at hard frost
  REAL(dp), PARAMETER :: AlbedoGlacierNirMin = 0.44_dp !! Albedo of glacier in the NIR range at at the melting point
  REAL(dp), PARAMETER :: AlbedoGlacierNirMax = 0.8_dp  !! Albedo of glacier in the NIR range at hard frost

CONTAINS

  SUBROUTINE land_surface_init_io(ntiles, nlct_help, IO_file_name)

  ! land_surface_init_io is called from init_land_surface

    USE mo_netCDF, ONLY : add_dim, IO_inq_dimid, IO_inq_dimlen, NF_MAX_NAME
    USE mo_io,  ONLY: IO_open, IO_READ, IO_close
    USE mo_linked_list, ONLY : TILES

    INTEGER, INTENT(in)  :: ntiles
    INTEGER, INTENT(in)  :: nlct_help
    CHARACTER(NF_MAX_NAME), INTENT(in) :: IO_file_name

    INTEGER :: IO_file_id, IO_dim_id
    INTEGER :: ntiles_help
    REAL(dp)  :: tile_values(ntiles)
    INTEGER :: i

    IF (p_parallel_io) THEN


       ! --- get number of landcover types actually used in grid boxes

       land_surface_file%opened = .FALSE.
       CALL IO_open(TRIM(IO_file_name), land_surface_file, IO_READ)
       IO_file_id = land_surface_file%file_id
       CALL IO_inq_dimid(IO_file_id, 'lct', IO_dim_id)
       CALL IO_inq_dimlen(IO_file_id, IO_dim_id, nlct)
       CALL IO_inq_dimid(IO_file_id, 'ntiles', IO_dim_id)
       CALL IO_inq_dimlen(IO_file_id, IO_dim_id, ntiles_help)

       CALL IO_close(land_surface_file)

       IF (nlct_help   /= nlct)   CALL finish('land_surface_init_io', 'nlct from IO_file differs from the one from LctLibrary')
       IF (ntiles_help /= ntiles) CALL finish('land_surface_init_io', 'ntiles from IO_file differs from the one in run.def')

    END IF

    IF (p_parallel) THEN
       CALL p_bcast(nlct, p_io)
    END IF

    IF (debug) CALL message('land_surface_init_io','Adding dimensions')
    CALL add_dim("lct",nlct,"available land cover types")
    ntiles_lct = ntiles   !vg:  Needs to be adapted if full flexibility for ntiles and ntiles_lct is needed 
    CALL add_dim("tiles_lct", ntiles_lct, "land cover type")

    DO i=1, ntiles
      tile_values(i)=float(i)
    END DO
    CALL add_dim ("tiles", ntiles, longname="land surface tile", units="", value=tile_values(:), levtyp=70, indx=TILES)

  END SUBROUTINE land_surface_init_io

  !=================================================================================================
  SUBROUTINE land_surface_init_memory(g_nland, l_nland, ntiles_help, land_surface, useDynveg, fileformat, diag_stream, stream)

  ! land_surface_init_memory is called from init_land_surface

    USE mo_jsbach,      ONLY : missing_value
    USE mo_linked_list, ONLY : LAND, TILES
    USE mo_memory_base, ONLY : new_stream, default_stream_setting, &
                               add =>add_stream_element
    USE mo_netCDF,      ONLY : max_dim_name
    USE mo_grib,        ONLY : land_table

    INTEGER,                 INTENT(in)    :: g_nland, l_nland, ntiles_help
    LOGICAL,                 INTENT(in)    :: useDynveg
    TYPE(land_surface_type), INTENT(inout) :: land_surface
    INTEGER,                 INTENT(in)    :: fileformat     ! output file format (grib/netcdf)            
    TYPE(t_stream),          POINTER       :: diag_stream
    TYPE(t_stream), POINTER, OPTIONAL :: stream

    INTEGER                     :: dim1p(1), dim1(1)
    INTEGER                     :: dim2p(2), dim2(2)
    CHARACTER(LEN=max_dim_name) :: dim1n(1), dim2n(2)

    IF (debug) CALL message('land_surface_init_memory','Entering ...')

    IF (ASSOCIATED(diag_stream)) THEN
       IO_diag => diag_stream
    ELSE
       CALL finish('land_surface_init_memory', 'Diagnostic stream not present')
    END IF

    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'land surface', filetype=fileformat)
          ! Set default stream options
          CALL default_stream_setting(stream, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
       IO_land_surface => stream
    ELSE
       ! Add new stream
       CALL new_stream(IO_land_surface, 'land surface', filetype=fileformat)
       ! Set default stream options
       CALL default_stream_setting(IO_land_surface, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF


    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    dim2p = (/ l_nland, ntiles_help /)
    dim2  = (/ g_nland, ntiles_help /)
    dim2n(1) = 'landpoint'
    dim2n(2) = 'tiles'

    CALL add(IO_land_surface, 'cover_type',        land_surface%cover_type_real, longname='Land Cover', &
             units='',     ldims=dim2p, gdims=dim2, dimnames=dim2n, code=11, lrerun=.FALSE., lpost=.FALSE.)
    CALL add(IO_land_surface, 'cover_fract',       land_surface%cover_fract,     longname='Land Cover Fraction', &
             units='',     ldims=dim2p, gdims=dim2, dimnames=dim2n, code=12, lrerun=.TRUE.,  lmiss=.TRUE., missval=missing_value)
    CALL add(IO_land_surface, 'albedo' ,           land_surface%albedo,          longname='Land Surface Albedo', &
             units='',     ldims=dim2p, gdims=dim2, dimnames=dim2n, code=13,  lpost=.FALSE.)
    CALL add(IO_land_surface, 'albedo_vis',        land_surface%albedo_vis,      longname='Surface Albedo in the Visible Range',  &
             units='',     ldims=dim2p, gdims=dim2, dimnames=dim2n, code=14,  lmiss=.TRUE., missval=missing_value)
    CALL add(IO_land_surface, 'albedo_nir' ,       land_surface%albedo_nir,      longname='Surface Albedo in the NIR', &
             units='',     ldims=dim2p, gdims=dim2, dimnames=dim2n, code=15,  lmiss=.TRUE., missval=missing_value)
    CALL add(IO_land_surface, 'elevation',         land_surface%elevation,       longname='Surface Altitude', &
             units='m',    ldims=dim1p, gdims=dim1, dimnames=dim1n, code=16,  lpost=.FALSE.)
    CALL add(IO_land_surface, 'orography_std_dev', land_surface%oro_std_dev,     longname='Standard Deviation of the Orography', &
             units='m',    ldims=dim1p, gdims=dim1, dimnames=dim1n, code=17,  lpost=.FALSE.)
    CALL add(IO_land_surface, 'forest_fract',      land_surface%forest_fract,    longname='Forest Fraction', &
             units='',     ldims=dim1p, gdims=dim1, dimnames=dim1n, code=18,  lpost=.FALSE.)
    CALL add(IO_land_surface, 'veg_ratio_actual_per_tile',land_surface%veg_ratio_actual_per_tile,            &
             longname='fraction of tile in vegetated part of grid box', units='',                            &
             ldims=dim2p,gdims=dim2,dimnames=dim2n,code=19, lpost=.FALSE.)
    CALL add(IO_land_surface, 'veg_ratio_max',     land_surface%veg_ratio_max,   longname='Maximum Vegetation Fraction', &
             units='',     ldims=dim1p, gdims=dim1, dimnames=dim1n, code=20, lpost=useDynveg)
    IF (useDynveg) THEN
       CALL add(IO_land_surface, 'rock_fract',     land_surface%rock_fract,      longname='Rock Fraction', &
                units='',  ldims=dim1p, gdims=dim1, dimnames=dim1n, code=23, lpost=.FALSE., contnorest=.TRUE.)
    ENDIF
    CALL add(IO_land_surface, 'swdown_acc',        land_surface%swdown_acc,      longname='Surface Downwelling Solar Radiation', &
             units='W m-2',ldims=dim1p, gdims=dim1, dimnames=dim1n, code=21, laccu=.true., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_land_surface, 'swdown_reflect_acc',land_surface%swdown_reflect_acc,longname='Surface Upwelling Solar Radiation', &
             units='W m-2',ldims=dim1p, gdims=dim1, dimnames=dim1n, code=22, laccu=.true., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag,         'albedo',            land_surface_diag%albedo,     longname='Land Surface Albedo', &
             units='',     ldims=dim1p, gdims=dim1, dimnames=dim1n, code=13, laccu=.FALSE., lpost=.FALSE., &
             lmiss=.TRUE., missval=missing_value)

    if (debug) CALL message('land_surface_init_memory','    exit')

    ALLOCATE(land_surface%cover_type(l_nland,ntiles_help))

  END SUBROUTINE land_surface_init_memory

  !=================================================================================================
  SUBROUTINE init_land_surface(grid, domain, lctlib, land_surface, surf_file, useDynveg, isRestart, &
       fileformat, read_cover_fract, IO_diag_stream, IO_stream)

    USE mo_decomposition, ONLY: global_decomposition
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_io,            ONLY: IO_open, IO_READ
    USE mo_netcdf,        ONLY: io_inq_dimid, io_inq_dimlen, io_inq_varid, io_get_var_double, nf_max_name
    Use mo_temp                          ! Provides temporary arrays

    TYPE(grid_type),   INTENT(in)        :: grid
    TYPE(domain_type), INTENT(in)        :: domain
    TYPE(lctlib_type), INTENT(in)        :: lctlib
    TYPE(land_surface_type), INTENT(inout) :: land_surface
    CHARACTER(nf_max_name), INTENT(in)   :: surf_file
    LOGICAL,           INTENT(in)        :: useDynveg
    LOGICAL,           INTENT(in)        :: isRestart
    INTEGER,           INTENT(in)        :: fileformat      ! output file format (grib/netcdf)
    LOGICAL,           INTENT(in)        :: read_cover_fract ! read cover fractions from initial
                                                             ! instead of restart file
    TYPE(t_stream),    POINTER           :: IO_diag_stream
    TYPE(t_stream),    POINTER, OPTIONAL :: IO_stream

    TYPE(FILE_INFO) :: IO_file
    INTEGER  :: IO_file_id, IO_var_id, IO_dim_id
    INTEGER  :: i, j, ilct, znlon, znlat
    REAL(dp) :: k
    integer  :: status

    if(debug) call message("init_land_surface","Start initialization of mo_land_surface")

    ntiles = land_surface%ntiles  ! land_surface%ntiles is set in calling jsbach_init

    call land_surface_init_io(ntiles, lctlib%nlct, surf_file)  ! Sets nlct  and test whether lctlib%nlct is consistent with IO file

    ! --- generate land surface stream

    CALL land_surface_init_memory(grid%nland, domain%nland, ntiles, land_surface, useDynveg, fileformat, &
         IO_diag_stream, stream=IO_stream)

    IF (.NOT. ASSOCIATED(IO_land_surface)) &
         CALL finish('init_land_surface','No memory stream for land surface')

    IF (p_parallel_io) THEN

       ! Open ini file
       CALL message('init_land_surface','Reading land surface fields from '//TRIM(land_surface_file%file_name))
       IO_file%opened = .FALSE.
       CALL IO_open(TRIM(land_surface_file%file_name), IO_file, IO_READ)
       IO_file_id = IO_file%file_id

       ! Check resolution
       CALL IO_inq_dimid  (IO_file_id, 'lat', IO_dim_id)
       CALL IO_inq_dimlen (IO_file_id, IO_dim_id, znlat)
       CALL IO_inq_dimid  (IO_file_id, 'lon', IO_dim_id)
       CALL IO_inq_dimlen (IO_file_id, IO_dim_id, znlon)

       IF (znlon /= grid%nlon .OR. znlat /= grid%nlat) THEN
          CALL finish('init_land_surface', 'Unexpected resolution:'//int2string(znlon)//' '//int2string(znlat))
       ENDIF

    ENDIF

    ! Temporary storage for local domain fields
    ALLOCATE(zreal2d(domain%ndim,domain%nblocks),STAT=status)
    if(status .ne. 0) call finish('init_land_surface','Allocation failure (1)')

    ! Land surface cover types
    IF (p_parallel_io) THEN
       ALLOCATE(zreal3d(grid%nlon,grid%nlat,ntiles_lct),STAT=status)
       if(status .ne. 0) call finish('init_land_surface','Allocation failure (2)')
       CALL IO_inq_varid(IO_file_id, 'cover_type', IO_var_id)
       CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
    END IF
    NULLIFY(zreal2d_ptr)
    DO i=1,ntiles_lct
       IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
       CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
       land_surface%cover_type_real(:,i) = PACK(zreal2d, MASK=domain%mask)
    END DO

    ! Distribute cover types to additional sub-tiling structure
    IF (ntiles /= ntiles_lct) THEN
       DO i = 1,ntiles/ntiles_lct
          DO j = 1,ntiles_lct
             land_surface%cover_type_real(:,(i-1)*ntiles_lct+j ) = land_surface%cover_type_real(:,j)
          END DO
       END DO
    END IF

    ! Convert cover_type from REAL to INTEGER
    land_surface%cover_type = NINT(land_surface%cover_type_real)

!! -- Read in landcover fractions from the initial file. 
!!    In restarted runs they are read from restart file, unless read_cover_fract is true (namelist jsbach_ctl)

    IF (.NOT. isRestart .OR. read_cover_fract ) THEN

       ALLOCATE(init_cover_fract(domain%nland,ntiles))

       IF (p_parallel_io) THEN
          zreal3d = 0._dp
          k = 0._dp
          DO j=1,grid%nlat
             DO i=1,grid%nlon
                IF (grid%mask(i,j)) THEN
                   k = k + 1._dp
                   zreal3d(i,j,1) = k
                END IF
             END DO
          END DO
       END IF

       ! Land surface cover fractions
       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'cover_fract', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          ! File may contain fraction of <nlct> land cover types as fraction of total grid cell (including ocean/sea ice).
          ! Scale surface coverage to fractions of total land cover. This also adjusts total cover if it is greater than 1.
          ALLOCATE(zzreal2d(grid%nlon, grid%nlat))
          zzreal2d = SUM(zreal3d, DIM=3)   ! Total land cover (possibly as fraction of grid cell)
          IF (ANY(grid%mask .AND. zzreal2d < 1.e-5_dp)) &
               CALL finish('init_land_surface', 'Land cover inconsistent with global land-sea mask')
          DO i=1,ntiles_lct
             WHERE (grid%mask) 
                zreal3d(:,:,i) = zreal3d(:,:,i) / zzreal2d(:,:)
             END WHERE
          END DO
          DEALLOCATE(zzreal2d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles_lct
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
          init_cover_fract(:,i) = PACK(zreal2d, MASK=domain%mask)
       ENDDO

       ! Distribute land cover fractions over ntiles (at the moment evenly over the additional tiling dimensions)
       IF (ntiles /= ntiles_lct) THEN
          init_cover_fract(:,1:ntiles_lct) = init_cover_fract(:,1:ntiles_lct) / (ntiles/ntiles_lct)
          DO i=2,ntiles/ntiles_lct
             init_cover_fract(:,(i-1)*ntiles_lct+1:i*ntiles_lct) = init_cover_fract(:,1:ntiles_lct)
          END DO
       END IF

       ! With Dynveg: Also read veg_ratio_max from initial file for consistency reasons
       IF (useDynveg) THEN
          ALLOCATE(init_veg_ratio_max(domain%nland))

          ! Read veg_ratio_max

          NULLIFY(zreal2d_ptr)
          IF (p_parallel_io) THEN
             ALLOCATE(zreal2d_ptr(grid%nlon,grid%nlat))
             CALL IO_inq_varid(IO_file_id, 'veg_ratio_max', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
          ENDIF
          CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
          init_veg_ratio_max(:) = PACK(zreal2d, MASK=domain%mask)
       END IF

    END IF
       
    IF (p_parallel_io) DEALLOCATE(zreal3d)


    ! ECHAM5 compatibility: read monthly vegetation ratio and use this as fraction x of first tile (vegetation), and
    ! use 1-x as fraction of second tile (bare soil). Only the third tile is kept from previous read (glacier)
    ALLOCATE(land_surface%veg_ratio(domain%nland, 0:13))
    IF (p_parallel_io) THEN
       ALLOCATE(zreal3d(grid%nlon,grid%nlat,12), STAT=status)
       if(status .ne. 0) call finish('init_land_surface','Allocation failure')
       CALL IO_inq_varid(IO_file_id, 'veg_fract', IO_var_id)
       CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
    END IF
    NULLIFY(zreal2d_ptr)
    DO i=1,12
       IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
       CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
       land_surface%veg_ratio(:,i) = PACK(zreal2d, MASK=domain%mask)
    END DO
    land_surface%veg_ratio(:,0) = land_surface%veg_ratio(:,12)
    land_surface%veg_ratio(:,13) = land_surface%veg_ratio(:,1)

    ! cover_fract is updated each time step from mo_jsbach_interface

    ! Set convenience logical fields <is_vegetation>, <is_forest>, <is_lake>, <is_glacier>, <is_bare_soil>
    ALLOCATE(land_surface%is_bare_soil(domain%nland, ntiles), &
             land_surface%is_vegetation(domain%nland, ntiles), &
             land_surface%is_forest(domain%nland, ntiles), &
             land_surface%is_lake(domain%nland, ntiles), &
             land_surface%is_glacier(domain%nland, ntiles), &
             land_surface%is_present(domain%nland, ntiles) )

    land_surface%is_bare_soil   = .FALSE.
    land_surface%is_vegetation  = .FALSE.
    land_surface%is_forest      = .FALSE.
    land_surface%is_lake        = .FALSE.
    land_surface%is_glacier     = .FALSE.
    land_surface%is_present     = .FALSE.

!!$    DO i=1,ntiles
       DO ilct=1,nlct
          WHERE (land_surface%cover_type == ilct)
             land_surface%is_vegetation = lctlib%NaturalVegFlag(ilct) .OR. lctlib%CropFlag(ilct)
             land_surface%is_forest     = lctlib%ForestFlag    (ilct)
             land_surface%is_lake       = lctlib%LakeFlag      (ilct)
             land_surface%is_glacier    = lctlib%GlacierFlag   (ilct)
             land_surface%is_bare_soil  = lctlib%BareSoilFlag  (ilct)
          END WHERE
       END DO
!!$    END DO

!!$    land_surface%is_vegetation = land_surface%is_vegetation .AND. land_surface%cover_fract > 0._dp
!!$    land_surface%is_forest     = land_surface%is_forest     .AND. land_surface%cover_fract > 0._dp
!!$    land_surface%is_lake       = land_surface%is_lake       .AND. land_surface%cover_fract > 0._dp
!!$    land_surface%is_glacier    = land_surface%is_glacier    .AND. land_surface%cover_fract > 0._dp
!!$    land_surface%is_bare_soil  = land_surface%is_bare_soil  .AND. land_surface%cover_fract > 0._dp
    ! Logical mask for tiles that should be processed by the land model
    land_surface%is_present    = land_surface%is_vegetation .OR. land_surface%is_bare_soil .OR. &
                                 land_surface%is_lake       .OR. land_surface%is_glacier

    ! If this is a restart run the model state variables are read in jsbach_init (Standalone) or the GCM
    ! by calling io_read_streams. We can therefore exit now:
    IF (isRestart) THEN
       NULLIFY(zreal2d_ptr)
       DEALLOCATE(zreal2d)
       IF (p_parallel_io) DEALLOCATE(zreal3d)
       module_initialized = .TRUE.
       RETURN
    ENDIF
    
    ! If this is not a restart run we continue and get the land surface model state from ini file

    land_surface%albedo = 0.0_dp
    land_surface%albedo_vis = 0.0_dp
    land_surface%albedo_nir = 0.0_dp
    land_surface_diag%albedo = 0.0_dp
    land_surface%swdown_acc = 0.0_dp
    land_surface%swdown_reflect_acc = 0.0_dp

    ! Elevation
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) THEN
       ALLOCATE(zreal2d_ptr(grid%nlon,grid%nlat))
       CALL io_inq_varid (IO_file_id, 'elevation', IO_var_id)
       CALL io_get_var_double (IO_file_id, IO_var_id, zreal2d_ptr)
    ENDIF
    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
    land_surface%elevation = PACK(zreal2d, MASK=domain%mask)
    
    ! Standard deviation of orography
    IF (p_parallel_io) THEN
       CALL io_inq_varid(IO_file_id, 'orography_std_dev', IO_var_id)
       CALL io_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
    END IF
    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
    land_surface%oro_std_dev = PACK(zreal2d, MASK=domain%mask)

    ! Forest fraction (temporary until JSBACH is verified against ECHAM5)
    IF (p_parallel_io) THEN
       CALL io_inq_varid(IO_file_id, 'forest_fract', IO_var_id)
       CALL io_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
    END IF
    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
    land_surface%forest_fract = PACK(zreal2d, MASK=domain%mask)

    ! maximal fraction of each grid box covered by leaves, assuming infinite LAI
    IF (p_parallel_io) THEN
       CALL io_inq_varid(IO_file_id, 'veg_ratio_max', IO_var_id)
       CALL io_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
    END IF
    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
    land_surface%veg_ratio_max = PACK(zreal2d, MASK=domain%mask)

    ! Grid cell area
    !ALLOCATE(grid%area(grid%nland))
    !IF (p_parallel_io) THEN
    !   CALL io_inq_varid (IO_file_id, 'area', IO_var_id)
    !   CALL io_get_var_double (IO_file_id, IO_var_id, zreal2d_ptr)
    !ENDIF
    !CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
    !grid%area = PACK(zreal2d, MASK=grid%mask_land)

    IF (p_parallel_io) THEN
       DEALLOCATE(zreal2d_ptr)
       DEALLOCATE(zreal3d)
    END IF
    NULLIFY(zreal2d_ptr)
    DEALLOCATE(zreal2d)

    IF (useDynveg) land_surface%rock_fract(:) = 0._dp

  END SUBROUTINE init_land_surface

  SUBROUTINE update_land_surface_fast(kidx, lctlib, land_surface, &
          surface_temperature, snow_fract, background_albedo, canopy_snow_fract, lai )

! This routine is coded for ECHAM5 compatibility!

    USE mo_constants,     ONLY: tmelt

    INTEGER, INTENT(in) :: kidx
    TYPE(lctlib_type),       INTENT(in)    :: lctlib
    TYPE(land_surface_type), INTENT(inout) :: land_surface
    REAL(dp), INTENT(in), DIMENSION(kidx, ntiles) ::   &
         surface_temperature, &
         snow_fract,          &                    ! Fraction of snow covered ground
         background_albedo,   &                    ! background albedo
         canopy_snow_fract,   &                    ! Fraction of snow covered canopy
         lai

    ! Local variables
    REAL(dp)    :: min_temp_snow_albedo             ! Temperature threshold below which maximum snow albedo is used
    REAL(dp)    :: snow_albedo_ground(kidx,ntiles)  ! Temperature dependend snow albedo over ground
    REAL(dp)    :: sky_view_fract(kidx,ntiles)      ! Fraction of bare ground below canopy for albedo calculation
    REAL(dp)    :: forest_fract(kidx,ntiles)        ! Factor for albedo from canopy (uses forest fraction and sky view factor)
    INTEGER :: kidx0, kidx1, i, j, ilct, itile
!!$    INTEGER :: ctype(kidx)

    kidx0 = kstart
    kidx1 = kend

    min_temp_snow_albedo = tmelt - 5.0_dp
    snow_albedo_ground = 0._dp

    do itile=1,ntiles
      do i=1,kidx
        j=kidx0+i-1
        ilct = land_surface%cover_type(j,itile)
        ! Temperature dependend snow albedo over ground
        IF (surface_temperature(i,itile) >= tmelt) THEN
          snow_albedo_ground(i,itile) = lctlib%AlbedoSnowMin(ilct)
        ELSE IF (surface_temperature(i,itile) < min_temp_snow_albedo) THEN
          snow_albedo_ground(i,itile) = lctlib%AlbedoSnowMax(ilct)
        ELSE
         snow_albedo_ground(i,itile) = lctlib%AlbedoSnowMin(ilct) + &
                     (tmelt - surface_temperature(i,itile)) * (lctlib%AlbedoSnowMax(ilct) - lctlib%AlbedoSnowMin(ilct)) / &
                     (tmelt - min_temp_snow_albedo)
        END IF
      END DO
    ENDDO

! The following loop gives identical results with or without  openmp
!
!!$    do itile=1,ntiles
!!$       ctype = land_surface%cover_type(kidx0:kidx1,itile)
!!$       DO ilct=1,nlct
!!$          WHERE (ctype == ilct)
!!$             ! Temperature dependend snow albedo over ground
!!$             WHERE (surface_temperature(:,itile) >= tmelt)
!!$                snow_albedo_ground(:,itile) = lctlib%MinSnowAlbedo(ilct)
!!$             ELSEWHERE (surface_temperature(:,itile) < min_temp_snow_albedo)
!!$                snow_albedo_ground(:,itile) = lctlib%MaxSnowAlbedo(ilct)
!!$             ELSEWHERE
!!$                snow_albedo_ground(:,itile) = lctlib%MinSnowAlbedo(ilct) + &
!!$                     (tmelt - surface_temperature(:,itile)) * (lctlib%MaxSnowAlbedo(ilct) - lctlib%MinSnowAlbedo(ilct)) / &
!!$                                                  (tmelt - min_temp_snow_albedo)
!!$             END WHERE
!!$          END WHERE
!!$       END DO
!!$    END DO

! The following loop gives different results with and without  openmp
!
!!$    DO ilct=1,nlct
!!$       WHERE (land_surface%cover_type(kidx0:kidx1,:) == ilct)
!!$          ! Temperature dependend snow albedo over ground
!!$          WHERE (surface_temperature >= tmelt)
!!$             snow_albedo_ground = lctlib%MinSnowAlbedo(ilct)
!!$          ELSEWHERE (surface_temperature < min_temp_snow_albedo)
!!$             snow_albedo_ground = lctlib%MaxSnowAlbedo(ilct)
!!$          ELSEWHERE
!!$             snow_albedo_ground = lctlib%MinSnowAlbedo(ilct) + &
!!$                  (tmelt - surface_temperature) * (lctlib%MaxSnowAlbedo(ilct) - lctlib%MinSnowAlbedo(ilct)) / &
!!$                                                  (tmelt - min_temp_snow_albedo)
!!$          END WHERE
!!$       END WHERE
!!$    END DO

    sky_view_fract(:,:) = EXP(-SkyViewFactor * MAX(lai(:,:),2.0_dp)) ! the value 2.0 should be replaced by StemArea

    ! calculate fraction for which albedo is computed from canopy
    forest_fract(:,:) = SPREAD(land_surface%forest_fract(kidx0:kidx1), DIM=2, NCOPIES=ntiles) * (1._dp - sky_view_fract(:,:))

    WHERE (land_surface%is_glacier(kidx0:kidx1,:))
       land_surface%albedo(kidx0:kidx1,:) = snow_albedo_ground
    ELSEWHERE
       ! albedo = weighted mean of albedo of ground below canopy and albedo of canopy
       land_surface%albedo(kidx0:kidx1,:) = &
            MAX(((1._dp - forest_fract) * (snow_fract * snow_albedo_ground + (1._dp - snow_fract) * background_albedo) + &
            forest_fract * (canopy_snow_fract * AlbedoCanopySnow + (1._dp - canopy_snow_fract) * background_albedo)),    &
            background_albedo)
    END WHERE

  END SUBROUTINE update_land_surface_fast

  !=================================================================================================

  SUBROUTINE init_albedo ! vg: hier her gehoert einiges aus update_land_surface!

    !! Read namelist and lctlib-file for albedo

    CALL config_albedo(albedo_params, albedo_options)
  END SUBROUTINE init_albedo

  !=================================================================================================

  SUBROUTINE config_albedo(albedo_params, albedo_options)

    USE mo_namelist,         ONLY: position_nml, POSITIONED
    USE mo_jsbach,           ONLY: nml_unit
    USE mo_doctor,           ONLY: nout

    TYPE(albedo_params_type), INTENT(inout) :: albedo_params
    TYPE(albedo_options_type), INTENT(inout) :: albedo_options

    !! Locals
    INTEGER :: read_status

    !! Namelist Parameters
    LOGICAL :: use_albedocanopy   !! true: use map of canopy albedo
                                  !! false: use PFT specific albedo values 
    INCLUDE 'albedo_ctl.inc'

    !! Read namelist albedo_ctl

    IF (p_parallel_io) THEN

       ! define default values
       use_albedocanopy = .FALSE.
       
       CALL position_nml ('ALBEDO_CTL', status=read_status)
       SELECT CASE (read_status)
       CASE (POSITIONED)
          READ (nml_unit, albedo_ctl)
          CALL message('config_albedo', 'Namelist ALBEDO_CTL: ')
          WRITE(nout, albedo_ctl)
       END SELECT
    END IF

    IF (p_parallel_io) THEN
       WRITE (message_text,*) 'use_albedocanopy: ', use_albedocanopy
       CALL message('config_albedo', message_text)

       albedo_options%UseAlbedocanopy = use_albedocanopy
    END IF

    IF (p_parallel) THEN
       CALL p_bcast(albedo_options%UseAlbedocanopy, p_io)
    END IF
 
  END SUBROUTINE config_albedo

  !=================================================================================================

  SUBROUTINE update_albedo(lctlib,                                   &
          l_inquire, l_trigrad, l_trigradm1, l_start,                &
          l_standalone,                                              &
          nidx, ntiles, cover_type, is_glacier, is_forest,           &
          veg_ratio_max, radiation_net_vis, radiation_net_nir,       &
          surface_temperature, snow_fract,                           &
          background_albedo_soil_vis, background_albedo_soil_nir,    &
          background_albedo_veg_vis, background_albedo_veg_nir,      &
          lai, canopy_snow_fract,                                    &
          albedo_vis, albedo_nir, albedo)

    USE mo_constants,     ONLY : tmelt
#if defined (__SX__) && defined (_OPENMP)
    USE omp_lib,          ONLY: omp_get_thread_num, omp_get_num_threads
#endif

    TYPE(lctlib_type),       INTENT(in)    :: lctlib
    LOGICAL, INTENT(in)  ::  l_inquire              ! trigger model initialisation
    LOGICAL, INTENT(in)  ::  l_trigrad              ! trigger full radiation time step
    LOGICAL, INTENT(in)  ::  l_trigradm1            ! trigger one time step before full radiation time step
    LOGICAL, INTENT(in)  ::  l_start                ! trigger first time step after model initialisation
    LOGICAL, INTENT(in)  ::  l_standalone           ! model runs without atmosphere (w/o echam e.g. and has no pre-init and no triger)
    INTEGER, INTENT(in)  ::  nidx                   ! domain
    INTEGER, INTENT(in)  ::  ntiles                 ! number of tiles
    INTEGER, INTENT(in), DIMENSION(nidx,ntiles) ::          &
         cover_type                                 ! number of plant functional type (PFT)
    LOGICAL, INTENT(in), DIMENSION(nidx,ntiles) ::          &
         is_glacier,                   &            ! glacier flag
         is_forest                                  ! forest flag
    REAL(dp), INTENT(in), DIMENSION(nidx)  ::                   &
         veg_ratio_max,                &            ! maximal fraction of the grid boy covered by vegetation
         radiation_net_vis,            &            ! net solar radiation in the visible range [W/m2]
         radiation_net_nir                          ! net solar radiation in the NIR range [W/m2]
    REAL(dp), INTENT(in), DIMENSION(nidx,ntiles)  ::            &
         surface_temperature,          &            !
         snow_fract,                   &            ! fraction of snow covered ground
         background_albedo_soil_vis,   &            ! background albedo vegetation NIR
         background_albedo_soil_nir,   &            ! background albedo soil NIR
         background_albedo_veg_vis,    &            ! background albedo vegetation visible
         background_albedo_veg_nir,    &            ! background albedo soil visible
         lai,                          &            ! leaf area index
         canopy_snow_fract                          ! fraction of snow covered canopy (forest)
    REAL(dp), INTENT(inout), DIMENSION(nidx,ntiles)  ::         &
         albedo_vis,                   &            ! albedo of the visible range
         albedo_nir,                   &            ! albedo of the NIR range
         albedo                                     ! albedo of the whole solar range

    ! Local variables
    REAL(dp)  ::  background_albedo_vis                ! albedo without snow in the visible range
    REAL(dp)  ::  background_albedo_nir                ! albedo without snow in the NIR range
    REAL(dp)  ::  background_albedo_canopy_vis         ! albedo without snow of the canopy in the visible range
    REAL(dp)  ::  background_albedo_canopy_nir         ! albedo without snow of the canopy in the NIR range
    REAL(dp)  ::  fraction_down_vis(nidx,ntiles)       ! fraction of solar downward radiation in the visible range
    REAL(dp)  ::  min_temp_snow_albedo                 ! temperature threshold below which maximum snow albedo is used
    REAL(dp)  ::  snow_albedo_soil_vis(nidx,ntiles)    ! albedo(temp.) of snow (covering the soil) in the visible range
    REAL(dp)  ::  snow_albedo_soil_nir(nidx,ntiles)    ! albedo(temp.) of snow (covering the soil) in the NIR range
    REAL(dp)  ::  sky_view_fract                       ! fraction of bare ground below canopy for albedo calculation
    REAL(dp)  ::  sky_view_fract_stem                  ! fraction added to snow covered canopy due to stem area
    LOGICAL  ::  l_rad(nidx)                           ! flag to indicate if it is day or night
    INTEGER  ::  i,itile, tid, nt

#if defined (__SX__) && defined (_OPENMP)
    IF (debug) THEN
       tid = omp_get_thread_num()
       nt = omp_get_num_threads()
       CALL message('update_albedo', 'OpenMP thread #'//int2string(tid)//' of '//int2string(nt)//' nidx: '//int2string(nidx))
    END IF
#endif
    min_temp_snow_albedo = tmelt - 5.0_dp
    snow_albedo_soil_vis(:,:) = 0.0_dp
    snow_albedo_soil_nir(:,:) = 0.0_dp

    IF (l_standalone .AND. l_start) THEN
       albedo(:,:) = 0.2_dp
       albedo_vis(:,:) = 0.1_dp
       albedo_nir(:,:) = 0.3_dp
    END IF

    l_rad(:) = .FALSE.
    IF (l_inquire) THEN
      l_rad(:) = .TRUE.
    ELSE
      DO i = 1,nidx
        IF (radiation_net_vis(i) + radiation_net_nir(i) > 1.0e-09_dp) l_rad(i) = .TRUE.
      END DO
    END IF

    ! calculate albedo of snow
    IF (l_standalone .OR. l_trigradm1 .OR. l_inquire) THEN
    DO itile = 1,ntiles
      DO i = 1,nidx
        IF (surface_temperature(i,itile) >= tmelt) THEN
          snow_albedo_soil_vis(i,itile) = lctlib%AlbedoSnowVisMin(cover_type(i,itile))
          snow_albedo_soil_nir(i,itile) = lctlib%AlbedoSnowNirMin(cover_type(i,itile))
        ELSE IF (surface_temperature(i,itile) < min_temp_snow_albedo) THEN
          snow_albedo_soil_vis(i,itile) = lctlib%AlbedoSnowVisMax(cover_type(i,itile))
          snow_albedo_soil_nir(i,itile) = lctlib%AlbedoSnowNirMax(cover_type(i,itile))
        ELSE
          snow_albedo_soil_vis(i,itile) = lctlib%AlbedoSnowVisMin(cover_type(i,itile)) +             &
            (tmelt - surface_temperature(i,itile)) * (lctlib%AlbedoSnowVisMax(cover_type(i,itile)) - &
            lctlib%AlbedoSnowVisMin(cover_type(i,itile))) / (tmelt - min_temp_snow_albedo)
          snow_albedo_soil_nir(i,itile) = lctlib%AlbedoSnowNirMin(cover_type(i,itile)) +             &
            (tmelt - surface_temperature(i,itile)) * (lctlib%AlbedoSnowNirMax(cover_type(i,itile)) - &
            lctlib%AlbedoSnowNirMin(cover_type(i,itile))) / (tmelt - min_temp_snow_albedo)
        END IF
      END DO
    END DO

    ! calculate new glacier albedo (visible and NIR) only at grid points with solar radiation or at model initialisation
    DO itile = 1,ntiles
      DO i = 1,nidx
        IF (l_rad(i)) THEN
          IF (surface_temperature(i,itile) >= tmelt) THEN
            albedo_vis(i,itile) = AlbedoGlacierVisMin
            albedo_nir(i,itile) = AlbedoGlacierNirMin          
          ELSE IF (surface_temperature(i,itile) < min_temp_snow_albedo) THEN
            albedo_vis(i,itile) = AlbedoGlacierVisMax
            albedo_nir(i,itile) = AlbedoGlacierNirMax
          ELSE
            albedo_vis(i,itile) = AlbedoGlacierVisMin + &
              (tmelt - surface_temperature(i,itile)) * (AlbedoGlacierVisMax - AlbedoGlacierVisMin) / &
              (tmelt - min_temp_snow_albedo)
            albedo_nir(i,itile) = AlbedoGlacierNirMin + &
              (tmelt - surface_temperature(i,itile)) * (AlbedoGlacierNirMax - AlbedoGlacierNirMin) / &
              (tmelt - min_temp_snow_albedo)
          END IF
        END IF
      END DO
    END DO
    END IF
    ! diagnose albedo of glaciers (whole solar spectral range) only at grid points with solar radiation or at model initialisation
    IF (l_standalone .OR. l_trigrad .AND. .NOT. l_inquire) THEN
    DO itile = 1,ntiles
      DO i = 1,nidx
        IF (l_rad(i) .AND. is_glacier(i,itile)) THEN
          fraction_down_vis(i,itile) = (radiation_net_vis(i) / (1.0_dp - albedo_vis(i,itile))) /   &
                                       ((radiation_net_vis(i) / (1.0_dp - albedo_vis(i,itile))) +  &
                                       (radiation_net_nir(i) / (1.0_dp - albedo_nir(i,itile))))
          albedo(i,itile) = fraction_down_vis(i,itile) * albedo_vis(i,itile) +  &
                            (1.0_dp - fraction_down_vis(i,itile)) * albedo_nir(i,itile)
        END IF
      END DO
    END DO
    END IF
    IF (l_inquire) THEN
    DO itile = 1,ntiles
      DO i = 1,nidx
        IF (is_glacier(i,itile)) THEN
          albedo(i,itile) = 0.5_dp * albedo_vis(i,itile) + 0.5_dp * albedo_nir(i,itile)
        END IF
      END DO
    END DO
    END IF
    ! calculate new albedo (visible and NIR) only at grid points with solar radiation or at model initialisation
    IF (l_standalone .OR. l_trigradm1 .OR. l_inquire) THEN
    DO itile = 1,ntiles
      DO i = 1,nidx
        ! albedo of the canopy (vegetation)
        IF (albedo_options%useAlbedocanopy) THEN
          background_albedo_canopy_vis = background_albedo_veg_vis(i,itile)
          background_albedo_canopy_nir = background_albedo_veg_nir(i,itile)
        ELSE
          background_albedo_canopy_vis = lctlib%AlbedoCanopyVIS(cover_type(i,itile))
          background_albedo_canopy_nir = lctlib%AlbedoCanopyNIR(cover_type(i,itile))
        END IF
        ! fraction for which albedo is computed from canopy
        sky_view_fract = 1.0_dp - (veg_ratio_max(i) * &
                         (1.0_dp - EXP(-0.5_dp * lai(i,itile))))
        ! fraction which is added to the (snow covered) canopy fraction due to stem area
        sky_view_fract_stem = sky_view_fract - (1.0_dp - (veg_ratio_max(i) * &
                              (1.0_dp - EXP(-0.5_dp * (lai(i,itile) + lctlib%StemArea(cover_type(i,itile)))))))

        IF (l_rad(i) .AND. .NOT. is_glacier(i,itile)) THEN
          background_albedo_vis = (1.0_dp - sky_view_fract) * background_albedo_canopy_vis + &
                                  sky_view_fract * background_albedo_soil_vis(i,itile)
          background_albedo_nir = (1.0_dp - sky_view_fract) * background_albedo_canopy_nir + &
                                  sky_view_fract * background_albedo_soil_nir(i,itile)
          ! albedo of forests = weighted mean of albedo of ground below canopy and albedo of canopy
          IF (is_forest(i,itile)) THEN
            albedo_vis(i,itile) =                                                                                 &
              MAX((sky_view_fract - sky_view_fract_stem) * (snow_fract(i,itile) * snow_albedo_soil_vis(i,itile) + &
              (1.0_dp - snow_fract(i,itile)) * background_albedo_soil_vis(i,itile)) +                             &
              sky_view_fract_stem  * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                             &
              (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_soil_vis(i,itile)) +                      &
              (1.0_dp - sky_view_fract) * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                        &
              (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_canopy_vis),                              &
              background_albedo_vis)
            albedo_nir(i,itile) =                                                                                 &
              MAX((sky_view_fract - sky_view_fract_stem) * (snow_fract(i,itile) * snow_albedo_soil_nir(i,itile) + &
              (1.0_dp - snow_fract(i,itile)) * background_albedo_soil_nir(i,itile)) +                             &
              sky_view_fract_stem * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                              &
              (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_soil_nir(i,itile)) +                      &
              (1.0_dp - sky_view_fract) * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                        &
              (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_canopy_nir),                              &
              background_albedo_nir)
          ELSE
            albedo_vis(i,itile) =                                              &
              MAX(snow_fract(i,itile) * snow_albedo_soil_vis(i,itile) +        &
              (1.0_dp - snow_fract(i,itile)) * background_albedo_vis,          &
              background_albedo_vis)
            albedo_nir(i,itile) =                                              &
              MAX(snow_fract(i,itile) * snow_albedo_soil_nir(i,itile) +        &
              (1.0_dp - snow_fract(i,itile)) * background_albedo_nir,          &
              background_albedo_nir)
          END IF
        END IF
      END DO
    END DO
    END IF
    ! diagnose albedo (whole solar spectral range) only at grid points with solar radiation or at model initialisation
    IF (l_standalone .OR. l_trigrad .AND. .NOT. l_inquire) THEN
    DO itile = 1,ntiles
      DO i = 1,nidx
        IF (l_rad(i) .AND. .NOT. is_glacier(i,itile)) THEN
          fraction_down_vis(i,itile) = (radiation_net_vis(i) / (1.0_dp - albedo_vis(i,itile))) /  &
                                       ((radiation_net_vis(i) / (1.0_dp - albedo_vis(i,itile))) + &
                                       (radiation_net_nir(i) / (1.0_dp - albedo_nir(i,itile))))
          albedo(i,itile) = fraction_down_vis(i,itile) * albedo_vis(i,itile) +                    &
                            (1.0_dp - fraction_down_vis(i,itile)) * albedo_nir(i,itile)
        END IF
      END DO
    END DO
    END IF
    IF (l_inquire) THEN
    DO itile = 1,ntiles
      DO i = 1,nidx
        IF (.NOT. is_glacier(i,itile)) THEN
          albedo(i,itile) = 0.5_dp * albedo_vis(i,itile) + 0.5_dp * albedo_nir(i,itile)
        END IF
      END DO
    END DO
    END IF

  END SUBROUTINE update_albedo

  SUBROUTINE land_surface_diagnostics(surface)

    USE mo_utils,        ONLY: average_tiles

    TYPE(land_surface_type), INTENT(in) :: surface

    !Local variables
    LOGICAL  :: mask(nidx,surface%ntiles)
    REAL(dp)     :: fract(nidx,surface%ntiles)
    INTEGER  :: ntiles
    INTEGER  :: kidx0, kidx1

    ntiles = surface%ntiles
    kidx0   = kstart
    kidx1   = kend

    ! Compute grid box averages
    mask  = surface%is_present(kidx0:kidx1,1:ntiles)
    fract = surface%cover_fract(kidx0:kidx1,1:ntiles)

    CALL average_tiles(surface%albedo(kidx0:kidx1,1:ntiles), mask, fract, land_surface_diag%albedo(kidx0:kidx1))

  END SUBROUTINE land_surface_diagnostics

!------------------------------------------------------------------------------
  SUBROUTINE scale_cover_fract (nidx, ntiles, glacier, cover_fract)
!
! !DESCRIPTION:
!
! Rescaling of cover fractions to assure that
!  - the sum of cover fractions is one
!  - all non-glacier tiles have a minimum vegetated fraction
!  - all glacier points have cover fraction of one on the first tile
!
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: nidx               ! vector length
    INTEGER,  INTENT(in)  :: ntiles             ! number of tiles
    LOGICAL,  INTENT(in)  :: glacier(:)         ! logical glacier mask
!
! !IN- and OUTPUT PARAMETERS:
! 
    REAL(dp), INTENT(inout) :: cover_fract(:,:) ! vegetated fraction

! !LOCAL VARIABLES:
!
    INTEGER   :: pft
    INTEGER   :: nsparce(nidx)     ! number of PFTs with less then fract_small vegetation
    REAL(dp)  :: sum_fract(nidx)   ! sum of all cover fractions
!------------------------------------------------------------------------------

    sum_fract(:) = 0._dp
    nsparce(:) = ntiles
    DO pft = 1,ntiles
       WHERE (cover_fract(:,pft) > fract_small)
          sum_fract(:) = sum_fract(:) + cover_fract(:,pft)
          nsparce(:) = nsparce(:) - 1
       END WHERE
    END DO
    DO pft = 1,ntiles
       WHERE (cover_fract(:,pft) > fract_small)
          cover_fract(:,pft) = cover_fract(:,pft) &
                              + cover_fract(:,pft)/sum_fract(:) * (1 - (sum_fract(:) + nsparce(:)*fract_small))
       ELSEWHERE
          cover_fract(:,pft) = fract_small
       END WHERE
    END DO

    WHERE (glacier(:))
       cover_fract(:,1) = 1._dp
    END WHERE
    DO pft = 2,ntiles
       WHERE (glacier(:))
          cover_fract(:,pft) = 0._dp
       END WHERE
    END DO

  END SUBROUTINE scale_cover_fract

END MODULE mo_land_surface
