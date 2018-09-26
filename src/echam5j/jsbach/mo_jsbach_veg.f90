Module mo_jsbach_veg
  !
  ! This module serves to
  ! (i) provide and initialize the infrastructure needed for the communication between the various submodels, and
  ! (ii) to call the initialization routines of the submodels.
  !
  ! I (CHR) would prefer to put the whole code of this module into mo_jsbach_interface, because it is actually
  ! an initialization of the interface. This would allow to make the fields (streams) initialized her PRIVATE instead of PUBLIC ---
  ! I am a real adict of data encapsulation!!
  !
  USE mo_netCDF, ONLY: FILE_INFO, NF_MAX_NAME
  USE mo_jsbach, ONLY: debug
  USE mo_jsbach_grid, ONLY: kstart, kend, nidx
  USE mo_linked_list, ONLY: t_stream
  USE mo_mpi, ONLY: p_parallel, p_parallel_io, p_bcast, p_io, p_pe
  USE mo_doctor, ONLY: nout, nerr
  Use mo_kind, Only : dp 
  Use mo_exception, Only : message, finish, int2string

!---wiso-code
  USE mo_wiso, ONLY: lwiso, nwiso
!---wiso-code-end

  IMPLICIT NONE

  ! === BEGIN OF PUBLIC PART =======================================================================================================

  ! --- public subroutines

  PUBLIC :: vegetation_type, init_vegetation, veg_diagnostics

  ! --- public parameters and variables

  ! --- public fields used for communication between submodels in the jsbach-interface AND ONLY IN THE JSBACH-INTERFACE!!!!!

  TYPE vegetation_type
     INTEGER              :: ntiles
     INTEGER              :: nroot_zones
     INTEGER, POINTER, DIMENSION(:,:)     :: & !! (nland, ntiles)
          veg_type                                   !! Index into LctLibrary
     REAL(dp),    POINTER, DIMENSION(:,:) :: & !! (nland, ntiles)
          canopy_conductance_bethy,          &       !! Canopy resistance
          canopy_conductance_limited,        &       !! Canopy resistance
          lai,                               &       !! Leaf Area Index
          lai_logrop,                        &       !! Leaf Area Index
          lai_max,                           &       !! Projected annual maximum LAI for each phenology type ...
          veg_fract_correction,              &       !! Correction factor for cover fraction 1-exp(-LAI_max/2)
          snow_depth_canopy,                 &       !! Snow depth on canopy [m]
          snow_fract_canopy                          !! Fraction of snow covered canopy
     REAL(dp),  POINTER, DIMENSION(:,:,:) :: & !! (nland, ntiles, 0:13)
          lai_clim                                   !! Climatological LAI
     REAL(dp),  POINTER, DIMENSION(:,:,:) :: & !! (nland, nroot_zones, ntiles)
          root_depth,                        &       !! Depth of root zones
          root_fract                                 !! Fraction of roots within each root zone
!---wiso-code
     REAL(dp),  POINTER, DIMENSION(:,:,:) :: &   !! (nland, nwiso, ntiles)
          wiso_snow_depth_canopy                 !! Snow depth on canopy [m]
!---wiso-code-end
  
  END TYPE vegetation_type

  TYPE vegetation_diag_type
     REAL(dp), POINTER, DIMENSION(:) :: &
          canopy_conductance_bethy,     &
          canopy_conductance_limited,   &
          lai,                          &
          lai_max,                      &
          lai_logrop,                   &
          snow_depth_canopy,            &
          snow_fract_canopy
!---wiso-code
     REAL(dp),  POINTER, DIMENSION(:,:) :: &
          wiso_snow_depth_canopy
!---wiso-code-end
  END TYPE vegetation_diag_type


  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ===============================================================================

  PRIVATE 

  ! === Private declarations =======================================================================================================

  INTEGER, SAVE :: nlct        = -1 !! Number of land cover types (including all PFTs)
!!$  INTEGER, SAVE :: ntiles      = -1 !! Number of tiles
  INTEGER, SAVE :: nroot_zones = -1 !! Number of root zones

  LOGICAL, SAVE   :: vegdist_changed = .FALSE. 
                                 !! Set to .TRUE. in the dynamic vegetation module
                                 !! if the distribution of vegetation types has
                                 !! changed in which case the global vegetation
                                 !! arrays have to be reinitialized.
 
  TYPE(t_stream),  POINTER, SAVE  :: IO_veg
  TYPE(t_stream),  POINTER, SAVE  :: IO_diag      !! Memory stream for diagnostic output

!---wiso-code
  TYPE(t_stream),  POINTER, SAVE  :: IO_veg_wiso  !! 
  TYPE(t_stream),  POINTER, SAVE  :: IO_diag_wiso !! Memory stream for soil diagnostic output - water isotopes
!---wiso-code-end

  TYPE(FILE_INFO), SAVE       :: veg_file
  TYPE(vegetation_diag_type), SAVE    :: veg_diag

Contains 

  !
  !=================================================================================================
  SUBROUTINE veg_init_io(IO_file_name)

    USE mo_netCDF,      ONLY: IO_inq_dimid, IO_inq_dimlen, add_dim
    USE mo_io,          ONLY: IO_open, IO_READ, IO_close
    USE mo_linked_list, ONLY: ROOTZONES

    CHARACTER(NF_MAX_NAME), INTENT(in) :: IO_file_name
    INTEGER :: IO_file_id, IO_dim_id
    INTEGER :: i    
    REAL(dp), ALLOCATABLE  :: root_values(:)


    IF (p_parallel_io) THEN
       ! Get number of root zones
!!$       call message('veg_init_io','Reading number of root zones from '//TRIM(IO_file_name))
       veg_file%opened = .FALSE.
       CALL IO_open(TRIM(IO_file_name), veg_file, IO_READ)
       IO_file_id = veg_file%file_id
!!$       CALL IO_inq_dimid(IO_file_id, 'rzone', IO_dim_id)
!!$       CALL IO_inq_dimlen(IO_file_id, IO_dim_id, nroot_zones)
       CALL IO_close(veg_file)
       nroot_zones = 1 !! ONLY FOR TESTING!!!
!!$       call message('veg_init_io',&
!!$            'WARNING: For TESTING number of root zones artificially preset to '//int2string(nroot_zones)//' !!')

    END IF

    IF (p_parallel) THEN
       CALL p_bcast(nroot_zones, p_io)
    ENDIF

    ALLOCATE (root_values(nroot_zones))
    DO i=1,nroot_zones 
       root_values(i) = float(i)
    END DO
    CALL add_dim ("root_zone", nroot_zones, longname="root zone", units="", value=root_values(:), levtyp=72, indx=ROOTZONES)
    DEALLOCATE (root_values)

  END SUBROUTINE veg_init_io
  !
  !=================================================================================================
  SUBROUTINE veg_init_memory(g_nland, l_nland, ntiles, veg, fileformat, diag_stream, diag_wiso_stream, stream, wiso_stream)

    USE mo_jsbach,      ONLY : missing_value
    USE mo_linked_list, ONLY : LAND, TILES, ROOTZONES
    USE mo_memory_base, ONLY : new_stream, default_stream_setting, &
                               add =>add_stream_element
    USE mo_netCDF,      ONLY : max_dim_name
    USE mo_grib,        ONLY : land_table
!---wiso-code
    USE mo_wiso, ONLY: lwiso, nwiso
!---wiso-code

    INTEGER,               INTENT(in)    :: g_nland, l_nland, ntiles
    TYPE(vegetation_type), INTENT(inout) :: veg
    INTEGER,               INTENT(in)    :: fileformat            ! output file format (grib/netcdf)
    TYPE(t_stream), POINTER              :: diag_stream
    TYPE(t_stream), POINTER, OPTIONAL    :: stream
!---wiso-code
    TYPE(t_stream), POINTER, OPTIONAL    :: diag_wiso_stream
    TYPE(t_stream), POINTER, OPTIONAL    :: wiso_stream
!---wiso-code-end

    INTEGER                     :: dim1p(1), dim1(1)
    INTEGER                     :: dim2p(2), dim2(2)
    INTEGER                     :: dim3p(3), dim3(3)
    CHARACTER(LEN=max_dim_name) :: dim1n(1), dim2n(2), dim3n(3)
!---wiso-code
    INTEGER                     :: dim4p(3), dim4(3)
    INTEGER                     :: dim5p(2), dim5(2)
    CHARACTER(LEN=max_dim_name) :: dim4n(3), dim5n(2) 
!---wiso-code-end

    integer                     :: status

    if (debug) CALL message('veg_init_memory','Entering ...')

    IF (ASSOCIATED(diag_stream)) THEN
       IO_diag => diag_stream
    ELSE
       CALL finish('veg_init_memory', 'Diagnostic stream not present')
    END IF

    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'veg', filetype=fileformat)
          ! Set default stream options
          CALL default_stream_setting(stream, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
       IO_veg => stream
    ELSE
       ! Add new stream
       CALL new_stream(IO_veg, 'veg', filetype=fileformat)
       ! Set default stream options
       CALL default_stream_setting(IO_veg, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF

!---wiso-code
    IF (lwiso) THEN

      IF (ASSOCIATED(diag_wiso_stream)) THEN
        IO_diag_wiso => diag_wiso_stream
      ELSE
        CALL finish('veg_init_memory', 'Isotope diagnostic stream not present')
      END IF

      IF (PRESENT(wiso_stream)) THEN
        IF (.NOT.ASSOCIATED(wiso_stream)) THEN
          ! Add new stream
          CALL new_stream(wiso_stream, 'veg_wiso', filetype = fileformat)
          ! Set default stream options
          CALL default_stream_setting(wiso_stream, table = land_table, repr = LAND, lpost = .TRUE., lrerun = .TRUE., &
                                      leveltype=TILES)
        ENDIF
        IO_veg_wiso => wiso_stream
      ELSE
        ! Add new stream
        CALL new_stream(IO_veg_wiso, 'veg_wiso', filetype = fileformat)
        ! Set default stream options
        CALL default_stream_setting(IO_veg_wiso, table = land_table, repr = LAND, lpost = .TRUE., lrerun = .TRUE., &
                                    leveltype = TILES)
      ENDIF
        
    ENDIF
!---wiso-code-end


    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    dim2p = (/ l_nland, ntiles /)
    dim2  = (/ g_nland, ntiles /)
    dim2n(1) = 'landpoint'
    dim2n(2) = 'tiles'

    dim3p = (/ l_nland, nroot_zones, ntiles /)
    dim3  = (/ g_nland, nroot_zones, ntiles /)
    dim3n(1) = 'landpoint'
    dim3n(2) = 'root_zone'
    dim3n(3) = 'tiles'

!---wiso-code
    dim4p = (/ l_nland, nwiso, ntiles /)
    dim4  = (/ g_nland, nwiso, ntiles /)    
    dim4n(1) = 'landpoint'
    dim4n(2) = 'isotopetyp'
    dim4n(3) = 'tiles'

    dim5p = (/ l_nland, nwiso /)
    dim5  = (/ g_nland, nwiso /)    
    dim5n(1) = 'landpoint'
    dim5n(2) = 'isotopetyp'
!---wiso-code-end

    CALL add(IO_veg, 'canopy_cond_bethy',  veg%canopy_conductance_bethy,       longname='Canopy Conductance BETHY', &
             units='m/s', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=105, lpost=.FALSE.)
    CALL add(IO_veg, 'canopy_cond_limited',veg%canopy_conductance_limited,     longname='Water Limited Canopy Conductance', &
             units='m/s', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=106, lpost=.FALSE.)
    CALL add(IO_veg, 'lai',                veg%lai,                            longname='Leaf Area Index', &
             units='',    ldims=dim2p, gdims=dim2, dimnames=dim2n, code=107, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_veg, 'lai_max',            veg%lai_max,                        longname='Maximum Annual Leaf Area Index', &
             units='',    ldims=dim2p, gdims=dim2, dimnames=dim2n, code=108, lpost=.FALSE.)
    CALL add(IO_veg, 'veg_fract_correction', veg%veg_fract_correction, &
             longname='Correction factor for cover fraction 1-exp(-LAI_max/2)', &
             units='',    ldims=dim2p, gdims=dim2, dimnames=dim2n, code=116, lpost=.FALSE.)
    CALL add(IO_veg, 'snow_depth_canopy',  veg%snow_depth_canopy,              longname='Snow Depth on Canopy', &
             units='m',   ldims=dim2p, gdims=dim2, dimnames=dim2n, code=109, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_veg, 'snow_fract_canopy',  veg%snow_fract_canopy,              longname='Snow Fraction on Canopy', &
             units='',    ldims=dim2p, gdims=dim2, dimnames=dim2n, code=110, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_veg, 'root_depth',         veg%root_depth,                     longname='Depth of Root Zone', &
              units='m',  ldims=dim3p, gdims=dim3, dimnames=dim3n, code=113, lpost=.FALSE., leveltype=ROOTZONES)
    CALL add(IO_veg, 'root_fraction',      veg%root_fract,                     longname='Fraction of Roots within each Root Zone', &
             units='',    ldims=dim3p, gdims=dim3, dimnames=dim3n, code=114, lpost=.FALSE., leveltype=ROOTZONES)
    CALL add(IO_veg, 'lai_logrop',         veg%lai_logrop,                     longname='Leaf Area Index from LOGROP', &
             units='',    ldims=dim2p, gdims=dim2, dimnames=dim2n, code=115, lpost=.FALSE.)

!---wiso-code
    IF (lwiso) THEN
       CALL add(IO_veg_wiso, 'wiso_snow_depth_canopy', veg%wiso_snow_depth_canopy, &
       longname='Snow Depth on Canopy - Water Isotopes', units='m', &
       ldims=dim4p, gdims=dim4, dimnames=dim4n, code=120, lmiss=.TRUE., missval=missing_value)
    END IF
!---wiso-code

    CALL add(IO_diag, 'canopy_cond_bethy',       veg_diag%canopy_conductance_bethy,  longname='Canopy Conductance BETHY', &
         units='m/s',     ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=105, lpost=.FALSE., &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'canopy_cond_limited',     veg_diag%canopy_conductance_limited, longname='Water Limited Canopy Conductance', &
         units='m/s',     ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=106, lpost=.FALSE., &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'lai',                     veg_diag%lai,                       longname='Leaf Area Index', &
         units='',        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=107, &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'lai_max',                 veg_diag%lai_max,                   longname='Maximum Annual Leaf Area Index', &
         units='',        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=108, lpost=.FALSE., &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'lai_logrop',              veg_diag%lai_logrop,                longname='Leaf Area Index from LOGROP', &
         units='',        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=115, lpost=.FALSE., &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'snow_depth_canopy',       veg_diag%snow_depth_canopy,         longname='Snow Depth on Canopy', &
         units='m',       ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=109, &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'snow_fract_canopy',       veg_diag%snow_fract_canopy,         longname='Snow Fraction on Canopy', &
         units='',        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=110, &
         lmiss=.TRUE., missval=missing_value)

!---wiso-code
    IF (lwiso) THEN
       CALL add(IO_diag_wiso, 'wiso_snow_depth_canopy', veg_diag%wiso_snow_depth_canopy, &
       longname='Snow Depth on Canopy - Water Isotopes', units='m', &
       ldims=dim5p, gdims=dim5, dimnames=dim5n, laccu=.FALSE., code=120, &
       lmiss=.TRUE., missval=missing_value)
    END IF
!---wiso-code

    if (debug) CALL message('veg_init_memory','    exit')

  END SUBROUTINE veg_init_memory
  !
  !=================================================================================================
  SUBROUTINE init_vegetation(grid, domain, land_surface, lctlib, vegetation, IO_file_name, &
       isRestart, fileformat, IO_diag_stream, IO_stream, IO_diag_wiso_stream, IO_wiso_stream)

    USE mo_jsbach_grid,   ONLY: grid_type, domain_type
    USE mo_land_surface,  ONLY: land_surface_type
    USE mo_jsbach_lctlib, ONLY: lctlib_type

    USE mo_decomposition, ONLY: global_decomposition
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_netcdf,        ONLY: nf_max_name, io_inq_varid, io_get_var_double
    USE mo_io,            ONLY: IO_open, IO_READ, IO_close
    Use mo_temp                          ! Provides temporary arrays
    USE mo_utils,         ONLY: get_current_month
!---wiso-code
    USE mo_time_control,  ONLY: lresume
    USE mo_wiso,          ONLY: lwiso, lwiso_rerun, nwiso, tnat
!---wiso-code

    TYPE(grid_type),         INTENT(in)    :: grid
    TYPE(domain_type),       INTENT(in)    :: domain
    TYPE(land_surface_type), INTENT(in)    :: land_surface
    TYPE(lctlib_type),       INTENT(in)    :: lctlib
    TYPE(vegetation_type),   INTENT(inout) :: vegetation
    CHARACTER(nf_max_name),  INTENT(in)    :: IO_file_name
    LOGICAL,                 INTENT(in)    :: isRestart
    INTEGER,                 INTENT(in)    :: fileformat   ! output file format (grib/netcdf) 
    TYPE(t_stream), POINTER           :: IO_diag_stream
    TYPE(t_stream), POINTER, OPTIONAL :: IO_stream
!---wiso-code
    TYPE(t_stream), POINTER, OPTIONAL  :: IO_diag_wiso_stream
    TYPE(t_stream), POINTER, OPTIONAL  :: IO_wiso_stream
!---wiso-code-end

    TYPE(FILE_INFO) :: IO_file
    INTEGER :: ntiles
    INTEGER :: IO_file_id, IO_var_id, IO_dim_id
    INTEGER :: i, j
    integer :: status,month
!---wiso-code
    INTEGER :: jt
!---wiso-code-end

    if(debug) call message("init_vegetation","Start initialization of mo_jsbach_veg")

    ntiles = land_surface%ntiles
    vegetation%ntiles = ntiles

    call veg_init_io(IO_file_name)  ! Sets nroot_zones
    vegetation%nroot_zones = nroot_zones

    ! --- generate vegetation stream

!---wiso-code
    IF (lwiso) THEN
       call veg_init_memory(grid%nland, domain%nland, ntiles, vegetation, fileformat, IO_diag_stream, IO_diag_wiso_stream, &
                            stream=IO_stream, wiso_stream=IO_wiso_stream)
    ELSE
       call veg_init_memory(grid%nland, domain%nland, ntiles, vegetation, fileformat, IO_diag_stream, stream=IO_stream)
    ENDIF
!---wiso-code-end

    IF (.NOT. ASSOCIATED(IO_veg)) &
         CALL finish('init_vegetation','No memory stream for vegetation')

    IF (p_parallel_io) THEN

       ! Open ini file
       call message('init_vegetation','Reading vegetation fields from '//TRIM(veg_file%file_name))
       IO_file%opened = .FALSE.
       CALL IO_open(TRIM(veg_file%file_name), IO_file, IO_READ)
       IO_file_id = IO_file%file_id

    ENDIF

    ! Temporary storage for local domain fields
    ALLOCATE(zreal2d(domain%ndim,domain%nblocks),STAT=status)
    if(status .ne. 0) call finish('init_vegetation','Allocation failure (1)')

    ! Vegetation cover
    ALLOCATE(vegetation%veg_type(domain%nland,ntiles))
    CALL message('init_vegetation', 'Setting veg_type')
    vegetation%veg_type = MERGE(land_surface%cover_type, -1, land_surface%is_vegetation)
    
call message("init_vegetation","CHECK -1")

	    CALL message('init_vegetation', 'veg_type set')

    !! Read in global LAI field for each calendar month
    ! Attention: The LAI read in here distinguishes not between different LCTs! Therefore here the total LAI
    !            in a grid cell is distributed to the various LCTs in proportion to their presence in the
    !            the grid cell ("vegetation_cover"). Note that this procedure assumes consistency between coverage
    !            data and LAI-data!! IT WOULD BE MUCH BETTER TO HAVE LCT-SPECIFIC DATA OF LAI!!! 

    
call message("init_vegetation","CHECK 0")

    ALLOCATE(vegetation%lai_clim(domain%nland,ntiles,0:13))
    IF (p_parallel_io) THEN
       ALLOCATE(zreal3d(grid%nlon,grid%nlat,12),STAT=status)
       if(status .ne. 0) call finish('init_vegetation','Allocation failure')
       CALL IO_inq_varid(IO_file_id, 'lai_clim', IO_var_id)
       CALL io_get_var_double(IO_file_id, IO_var_id, zreal3d)
!!$       month = get_current_month()
!!$       CALL message('init_vegetation', '   ... reading climatological LAI for month '//int2string(month))
    ENDIF
    
    NULLIFY(zreal2d_ptr)
    DO i=1,12
       IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
       CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
       ! set clim. LAI to total LAI
       vegetation%lai_clim(:,:,i) = MAX(SPREAD(PACK(zreal2d, MASK=domain%mask),DIM=2, NCOPIES=ntiles),0.001_dp)
    END DO

    DO j=1,12
       vegetation%lai_clim(:,:,j) = MERGE(vegetation%lai_clim(:,:,j),0._dp,land_surface%is_vegetation)
    END DO
    vegetation%lai_clim(:,:,0) = vegetation%lai_clim(:,:,12)
    vegetation%lai_clim(:,:,13) = vegetation%lai_clim(:,:,1)

    
    IF (p_parallel_io) DEALLOCATE(zreal3d)

!---wiso-code
    IF (lwiso) THEN

    ! For a restart without previous isotope diagnostics: 
    ! Initialize wiso_snow_depth_canopy from default JSBACH snow_depth_canopy field
    ! (assume at start a delta value of zero permill) 
    
      IF (lresume .AND. (.NOT. lwiso_rerun)) THEN
        DO jt=1,nwiso
          DO i=1,ntiles
            vegetation%wiso_snow_depth_canopy(:,jt,i) = tnat(jt)*vegetation%snow_depth_canopy(:,i)
          END DO
        END DO
      END IF

    END IF
!---wiso-code-end

    ! If this is a restart run the model state variables are read in jsbach_init (Standalone) or the GCM
    ! by calling io_read_streams. We can therefore exit now:
    IF (isRestart) THEN
       DEALLOCATE(zreal2d)
       IF (p_parallel_io) CALL IO_close(IO_file)
       RETURN
    ENDIF
    ! ----------------------------------------------------------------------------------------------------------------------
    ! If this is not a restart run we continue and get the vegetation model state from ini file
    !

    ! Root Depth
!!$    IF (p_parallel_io) THEN
!!$       ALLOCATE(zreal4d(nlon,nlat,nroot_zones,nlct))
     !!  CALL IO_inq_varid(IO_file_id, 'root_depth', IO_var_id)
     !!  CALL IO_get_var_double(IO_file_id, IO_var_id, zreal4d)
!!$       zreal4d = 75.
!!$    ENDIF
!!$    NULLIFY(zreal2d_ptr)
!!$    DO j=1,nlct
!!$       DO i=1,nroot_zones
!!$          IF (p_parallel_io) zreal2d_ptr => zreal4d(:,:,i,j)
!!$          CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
!!$          root_depth(:,i,j) = PACK(zreal2d, MASK=grid%mask_land)
!!$       ENDDO
!!$    ENDDO
    vegetation%root_depth = 1._dp  !! FOR TESTING (one root zone with same depth as the one soil layer!)
!!$
!!$    IF (ANY(SUM(root_depth,DIM=2) <= 0)) &
!!$         CALL finish('JSBACH - init_vegegation',&
!!$                     'Root zone depths must sum to a value > 0')
!!$    
!!$    ! Root Fraction
!!$    IF (p_parallel_io) THEN
     !!  CALL IO_inq_varid(IO_file_id, 'root_fract', IO_var_id)
     !!  CALL IO_get_var_double(IO_file_id, IO_var_id, zreal4d)
!!$       zreal4d = 1.
!!$    ENDIF
!!$    NULLIFY(zreal2d_ptr)
!!$    DO j=1,nlct
!!$       DO i=1,nroot_zones
!!$          IF (p_parallel_io) zreal2d_ptr => zreal4d(:,:,i,j)
!!$          CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
!!$          root_fract(:,i,j) = PACK(zreal2d, MASK=grid%mask_land)
!!$       ENDDO
!!$    ENDDO
!!$    IF (p_parallel_io) DEALLOCATE(zreal4d)
!!$    DEALLOCATE(zreal2d)
!!$
!!$    ALLOCATE(zreal2d(grid%nland,nlct))
!!$    zreal2d = SUM(root_fract,DIM=2)
!!$    IF (ANY(zreal2d /= 1.)) THEN
!!$       CALL message('JSBACH - init_vegetation',&
!!$                    'Sum of root zone fractions not one, normalizing ...')
!!$       FORALL(j=1:nroot_zones) &
!!$            root_fract(:,j,:) = root_fract(:,j,:) / zreal2d
!!$    ENDIF

    DEALLOCATE(zreal2d)

    IF (p_parallel_io) CALL IO_close(IO_file)

    vegetation%snow_fract_canopy = 0._dp
    vegetation%snow_depth_canopy = 0._dp

!---wiso-code
    IF (lwiso) THEN
      vegetation%wiso_snow_depth_canopy = 0._dp  
    END IF
!---wiso-code-end

    vegetation%lai = 0._dp

    vegetation%lai_max = 0._dp
    vegetation%veg_fract_correction = 1.0_dp
    
    DO i=1,ntiles
       WHERE (land_surface%is_vegetation(:,i))
          vegetation%lai_max(:,i) = &
!!$               MAX(MIN(MAXVAL(vegetation%lai_clim(:,i,:),DIM=2), &
!!$               lctlib%MaxLAI(land_surface%cover_type(:,i)) - 0.05_dp), 0.001_dp)
               lctlib%MaxLAI(land_surface%cover_type(:,i))
          WHERE (lctlib%PhenologyType(land_surface%cover_type(:,i)) == 5)
             vegetation%veg_fract_correction(:,i) = 1.0_dp - exp(-vegetation%lai_max(:,i) / 3.0_dp) ! Crops more clumpy ..
          ELSEWHERE
             vegetation%veg_fract_correction(:,i) = 1.0_dp - exp(-vegetation%lai_max(:,i) / 2.0_dp) ! .. than natural vegetation.
          END WHERE
       END WHERE
    END DO

    veg_diag%canopy_conductance_bethy = 0._dp
    veg_diag%canopy_conductance_limited = 0._dp
    veg_diag%lai = 0._dp
    veg_diag%lai_max = 0._dp
    veg_diag%lai_logrop = 0._dp
    veg_diag%snow_depth_canopy = 0._dp
    veg_diag%snow_fract_canopy = 0._dp

!---wiso-code
    IF (lwiso) THEN
      veg_diag%wiso_snow_depth_canopy = 0._dp  
    END IF
!---wiso-code-end

    if(debug) call message("init_vegetation","Initialization of mo_jsbach_veg finished.")

  END SUBROUTINE init_vegetation

  !
  !=================================================================================================
  SUBROUTINE veg_diagnostics(surface, veg)

    USE mo_land_surface, ONLY: land_surface_type
    USE mo_utils,        ONLY: average_tiles
!---wiso-code
    USE mo_wiso,         ONLY: lwiso, nwiso
!---wiso-code-end

    TYPE(land_surface_type), INTENT(in) :: surface
    TYPE(vegetation_type),   INTENT(inout) :: veg

    !Local variables
    INTEGER  :: ntiles
    INTEGER  :: kidx0, kidx1
!---wiso-code
    INTEGER  :: jt
!---wiso-code-end

    ntiles  = veg%ntiles
    kidx0   = kstart
    kidx1   = kend

    CALL average_tiles(veg%canopy_conductance_bethy(kidx0:kidx1,1:ntiles), &
         surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
         veg_diag%canopy_conductance_bethy(kidx0:kidx1))
    CALL average_tiles(veg%canopy_conductance_limited(kidx0:kidx1,1:ntiles), &
         surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
         veg_diag%canopy_conductance_limited(kidx0:kidx1))
    CALL average_tiles(veg%lai(kidx0:kidx1,1:ntiles), &
         surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
         veg_diag%lai(kidx0:kidx1))
    CALL average_tiles(veg%lai_max(kidx0:kidx1,1:ntiles), &
         surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
         veg_diag%lai_max(kidx0:kidx1))
    CALL average_tiles(veg%lai_logrop(kidx0:kidx1,1:ntiles), &
         surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
         veg_diag%lai_logrop(kidx0:kidx1))
    CALL average_tiles(veg%snow_depth_canopy(kidx0:kidx1,1:ntiles), &
         surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
         veg_diag%snow_depth_canopy(kidx0:kidx1))
    CALL average_tiles(veg%snow_fract_canopy(kidx0:kidx1,1:ntiles), &
         surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
         veg_diag%snow_fract_canopy(kidx0:kidx1))

!---wiso-code
    IF (lwiso) THEN
      DO jt=1,nwiso
        CALL average_tiles(veg%wiso_snow_depth_canopy(kidx0:kidx1,jt,1:ntiles), &
           surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
           veg_diag%wiso_snow_depth_canopy(kidx0:kidx1,jt))
      END DO
    END IF
!---wiso-code-end

  END SUBROUTINE veg_diagnostics

End module mo_jsbach_veg

!Local Variables:
!mode: f90
!fill-column: 100
!End:
