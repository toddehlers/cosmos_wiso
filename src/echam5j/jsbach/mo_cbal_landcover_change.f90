module mo_cbal_landcover_change

  !! Christian Reick, 2006-08-15

  use mo_kind,          only: dp
  use mo_mpi,           only: p_parallel_io
  use mo_exception,     only: finish, message,int2string
  use mo_linked_list,   only: t_stream

  implicit none

  type landcover_change_type
     real(dp),pointer :: LCC_sum_box_C2atmos(:)      !! Carbon released to atmosphere by landcover change
                                                     !! .. (summed through whole output interval) [mol(C)/m^2(grid box)]
     real(dp),pointer :: LCC_sum_box_C2fastSoilPool(:) !! Carbon released from green and reserve pool to fast soil pool by landcover change
                                                       !! .. (summed through whole output interval) [mol(C)/m^2(grid box)]
     real(dp),pointer :: LCC_sum_box_C2slowSoilPool(:) !! Carbon released from wood pool to slow soil pool by landcover change
                                                       !! .. (summed through whole output interval) [mol(C)/m^2(grid box)]
     real(dp),pointer :: LCC_box_C2A_wholeRun(:)     !! Amount of landcover change emissions summed over the whole run [mol(C)/m^2(grid box)]
     real(dp),pointer :: LCC_coverFract_target(:,:)  !! Fraction of vegetated part of a gridbox that due to landcover change should be ..
                                                     !! .. reached at the end of the year. (nland x ntiles)
  end type landcover_change_type

  private 

  public :: init_landcover_change
  public :: do_landcover_change

  !! --- Parameters ---------------------------------------------------------------

  real(dp),parameter :: frac_wood_2_atmos    = 0.5 !! Fraction of carbon from the wood pool to be released into the atmosphere on human landcover change
  real(dp),parameter :: frac_green_2_atmos   = 0.5 !! Fraction of carbon from the green pool to be released into the atmosphere on human landcover change 
  real(dp),parameter :: frac_reserve_2_atmos = 0.5 !! Fraction of carbon from the reserve pool to be released into the atmosphere on human landcover change

  character(len=*),parameter :: coverFractVarName = "cover_fract" !! name of variable for landcover fractions in netcdf input file


  !! --- private variables ---------------------------------------------------------

  TYPE(landcover_change_type), SAVE :: landcover_change
  TYPE(t_stream), POINTER, SAVE     :: LCC_stream   !! Memory stream for model state
  TYPE(t_stream), POINTER           :: IO_diag      !! Memory stream for diagnostic output

CONTAINS


  !! --- init_landcover_change() ----------------------------------------------------

  SUBROUTINE init_landcover_change(g_nland, l_nland, ntiles, isRestart, fileformat, diag_stream, stream)

    USE mo_land_surface, ONLY: land_surface_type
    USE mo_linked_list,  ONLY: LAND, TILES
    USE mo_grib,         ONLY: land_table
    USE mo_netCDF,       ONLY: max_dim_name
    USE mo_time_control, ONLY: get_date_components
    USE mo_memory_base,  ONLY: new_stream,default_stream_setting, add =>add_stream_element
    USE mo_exception,    ONLY: message,real2string

    INTEGER, INTENT(in)               :: g_nland, l_nland
    INTEGER, INTENT(in)               :: ntiles          !! number of tiles
    LOGICAL, INTENT(in)               :: isRestart       !! restart flag
    INTEGER, INTENT(in)               :: fileformat      !! output file format (grib/netcdf)
    TYPE(t_stream), POINTER           :: diag_stream
    TYPE(t_stream), POINTER, OPTIONAL :: stream          !! stream to which local stream can optionally be associated

    ! --- local variables

    INTEGER                     :: dim1p(1), dim1(1)
    INTEGER                     :: dim3p(2), dim3(2)
    CHARACTER(LEN=max_dim_name) :: dim1n(1), dim3n(2)

    !! --- printout parameters 

    call message("init_landcover_change()","=============================================================")
    call message("init_landcover_change()","   frac_wood_2_atmos="//trim(real2string(frac_wood_2_atmos)))
    call message("init_landcover_change()","  frac_green_2_atmos="//trim(real2string(frac_green_2_atmos)))
    call message("init_landcover_change()","frac_reserve_2_atmos="//trim(real2string(frac_reserve_2_atmos)))
    call message("init_landcover_change()","=============================================================")

    IF (ASSOCIATED(diag_stream)) THEN
       IO_diag => diag_stream
    ELSE
       CALL finish('init_landcover_change', 'Diagnostic stream not present')
    END IF

    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(LCC_Stream, 'landCoverChange', filetype=fileformat)
          ! Set default stream options
          CALL default_stream_setting(stream, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       END IF
       LCC_Stream => stream
    ELSE
       ! Add new stream
       CALL new_stream(LCC_Stream, 'landCoverChange', filetype=fileformat)
       ! Set default stream options
       CALL default_stream_setting(LCC_Stream, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    END IF

    !! --- add stream elements

    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    dim3p = (/ l_nland, ntiles /)
    dim3  = (/ g_nland, ntiles /)
    dim3n(1) = 'landpoint'
    dim3n(2) = 'tiles'

    CALL add(LCC_stream, 'LCC_box_C2A_wholeRun', landcover_change%LCC_box_C2A_wholeRun,                          &
             longname="Carbon released to atmosphere by landcover change summed through whole run",        &
             units='mol(C) m-2(grid box)',                                                                     &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=200, lrerun=.TRUE., contnorest=.TRUE., lpost=.TRUE.)

    CALL add(LCC_Stream, 'LCC_coverFract_target', landcover_change%LCC_coverFract_target,                        &
             longname="Vegetation cover fractions that should be reached at the end of the year as a result of landcover change.", &
             units='',                                                                     &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=201, lrerun=.TRUE., contnorest=.TRUE., lpost=.true.)

    CALL add(IO_diag, 'LCC_sum_box_C2atmos', landcover_change%LCC_sum_box_C2atmos,                            &
             longname="Sum of carbon emitted to atmosphere from land cover change",                            &
             units='mol(C) m-2(vegetated area)',                                                               &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=161, lpost=.TRUE.)
    CALL add(IO_diag, 'LCC_sum_box_C2fastSoilPool', landcover_change%LCC_sum_box_C2fastSoilPool,              &
             longname="Sum of green and reserve carbon relocated by landcover change to fast soil pool", &
             units='mol(C) m-2(vegetated area)',                                                               &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=162, lpost=.TRUE.)
    CALL add(IO_diag, 'LCC_sum_box_C2slowSoilPool', landcover_change%LCC_sum_box_C2slowSoilPool,              &
             longname="Wood carbon relocated by landcover change to slow soil pool", &
             units='mol(C) m-2(vegetated area)',                                                               &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=163, lpost=.TRUE.)


    landcover_change%LCC_sum_box_C2atmos        = 0.0_dp
    landcover_change%LCC_sum_box_C2fastSoilPool = 0.0_dp
    landcover_change%LCC_sum_box_C2slowSoilPool = 0.0_dp

    IF(isRestart ) THEN
       landcover_change%LCC_coverFract_target      = 0.0_dp
       landcover_change%LCC_box_C2A_wholeRun       = 0.0_dp
    END IF

  END SUBROUTINE init_landcover_change


  !! --- do_landcover_change() ----------------------------------------------------
  !!
  !! This routine:
  !!    1. Replaces the old map of land cover fractions by a new map that is read from file
  !!    2. Performs the necessary changes in the carbon pools that go along
  !!       with landcover changes
  !!
  !! NOTE: This routine assumes that the restart files have already been
  !! read, because they contain the old landcover distribution that is needed
  !! to derive the changes in the carbon pools
  !!
  subroutine do_landcover_change(grid,domain,ntiles,is_glacier,cbalance,veg_ratio_max,veg_fract_correction,&
                                 landcover_fract_current,CO2_emissions)

    use mo_jsbach_grid,      only: grid_type,domain_type,kstart,kend,nidx
    use mo_cbal_bethy,       only: cbalance_type
    use mo_cbal_cpools,      only: relocate_carbon
    USE mo_time_control,     ONLY: get_date_components, current_date, previous_date, get_year_day
    use mo_jsbach_constants, only: molarMassCO2_kg
    use mo_time_conversion,  only: year_len
    use mo_time_control,     only: l_trigfiles
    USE mo_time_control,     ONLY: delta_time    !! time step in seconds
    USE mo_land_surface,     ONLY: scale_cover_fract

    integer,                intent(in)    :: ntiles                       !! number of tiles of the land points
    type(grid_type),        intent(in)    :: grid
    type(domain_type),      intent(in)    :: domain
    logical,                intent(in)    :: is_glacier(:)                !! logical glacier mask
    type(cbalance_type),    intent(inout) :: cbalance                     !! the carbon pools that will be changed in this call
    real(dp),               intent(in)    :: veg_ratio_max(:)             !! maximum fraction of grid cell that can be covered by vegetation
    real(dp),               intent(in)    :: veg_fract_correction(:,:)    !! Correction factor for cover fractions 1-exp(-LAI_max/2)
    real(dp),               intent(inout) :: landcover_fract_current(:,:) !! On call the old landcover fractions, on return the new ones
    real(dp),               intent(out)   :: CO2_emissions(:)             !! The CO2-emissions from landcover change in kg(CO2)/m^2 in this timestep

    !! locals

    character(len=1024) :: fileName
    integer             :: current_year,previous_year,dayNo_in_current_year
    integer             :: current_day,previous_day
    real(dp)            :: C2atmos(1:nidx)
    real(dp)            :: C2slowSoilPool(1:nidx)
    real(dp)            :: C2fastSoilPool(1:nidx)
    real(dp)            :: landcover_fract_new(1:nidx,1:ntiles)
    integer             :: noOfDays_in_current_year

    !! Preparations

    call get_date_components(current_date,YEAR=current_year,DAY=current_day)
    call get_date_components(previous_date,YEAR=previous_year,DAY=previous_day)

    !! -- Restart summations after a restart

    if(l_trigfiles) then  !! Restart summations with begin of a new output interval
       landcover_change%LCC_sum_box_C2atmos(kstart:kend) = 0.0_dp
       landcover_change%LCC_sum_box_C2fastSoilPool(kstart:kend) = 0.0_dp
       landcover_change%LCC_sum_box_C2slowSoilPool(kstart:kend) = 0.0_dp
    endif

    !! --- If at the first time step of a new year or first time step of new run ==> read in new landcover map
    !!     (run might start in the middle of a year)


    if(current_year /= previous_year) then

       !! Generic name of land surface file with new landcover
       IF (current_year < 1000 .AND. current_year >= 100) THEN
          fileName = "cover_fract.0"//TRIM(int2string(current_year))//".nc"
       ELSE IF (current_year >= 1000 .AND. current_year < 10000) THEN
          fileName = "cover_fract."//TRIM(int2string(current_year))//".nc"
       ELSE
          CALL finish('do_landcover_change','Only years between 100 and 9999 supperted currently')
       END IF

       call read_landcover_fractions(fileName,coverFractVarName,grid,domain,ntiles) !! Updates LCC_coverFract_target(:,:) from external file
       CALL scale_cover_fract (nidx, ntiles, is_glacier(:), landcover_fract_current(:,:))

    end if

    !! --- Check for exiting the routine without landcover change

    if(current_day == previous_day) then !! We are not at the first time step of a day ..
       CO2_emissions(:) = 0.0_dp         !! .. so we do no landcover change and have no CO2-emissions
       return                            !! .. and exit the routine
    end if

    !! --- determine new cover fractions (we are at the first time step of a day)

    noOfDays_in_current_year = year_len(current_date)
    dayNo_in_current_year = floor(get_year_day(current_date))
    landcover_fract_new(1:nidx,:) =                                                                                 &
                 landcover_fract_current(:,:)                                                                       &
                 + ( landcover_change%LCC_coverFract_target(kstart:kend,:) - landcover_fract_current(1:nidx,:) )    &
                 / real(1 + noOfDays_in_current_year - dayNo_in_current_year)

    CALL scale_cover_fract (nidx, ntiles, is_glacier(:), landcover_fract_current(:,:))

    
    !! --- Do the changes in the carbon pools and compute amount of carbon to be released to atmosphere

    call relocate_carbon(landcover_fract_current(1:nidx,1:ntiles),landcover_fract_new(1:nidx,1:ntiles), &
                         veg_fract_correction(1:nidx,1:ntiles), epsilon(1._dp),                         &
                         cbalance%Cpool_green(kstart:kend,1:ntiles),                                    &
                         cbalance%Cpool_woods(kstart:kend,1:ntiles),                                    &
                         cbalance%Cpool_reserve(kstart:kend,1:ntiles),                                  &
                         cbalance%Cpool_litter_leaf(kstart:kend,1:ntiles),                              &
                         cbalance%Cpool_litter_wood(kstart:kend,1:ntiles),                              &
                         cbalance%Cpool_fast(kstart:kend,1:ntiles),                                     &
                         cbalance%Cpool_slow(kstart:kend,1:ntiles),                                     &
                         frac_wood_2_atmos = frac_wood_2_atmos,                                         &
                         frac_green_2_atmos = frac_green_2_atmos,                                       &
                         frac_reserve_2_atmos = frac_reserve_2_atmos,                                   &
                         carbon_2_atmos = C2atmos(1:nidx),                                              &
                         carbon_2_fastSoilPool = C2fastSoilPool(1:nidx),                                &
                         carbon_2_slowSoilPool = C2slowSoilPool(1:nidx))

    !! Sum landcover change emissions

    landcover_change%LCC_sum_box_C2atmos(kstart:kend) = &
                                    landcover_change%LCC_sum_box_C2atmos(kstart:kend) + veg_ratio_max(1:nidx)*C2atmos(1:nidx)
    landcover_change%LCC_sum_box_C2fastSoilPool(kstart:kend) = &
                      landcover_change%LCC_sum_box_C2fastSoilPool(kstart:kend) + veg_ratio_max(1:nidx)*C2fastSoilPool(1:nidx)
    landcover_change%LCC_sum_box_C2slowSoilPool(kstart:kend) = &
                      landcover_change%LCC_sum_box_C2slowSoilPool(kstart:kend) + veg_ratio_max(1:nidx)*C2slowSoilPool(1:nidx)

    !! Sum landcover change emissions through whole run

    landcover_change%LCC_box_C2A_wholeRun(kstart:kend) = &
                                   landcover_change%LCC_box_C2A_wholeRun(kstart:kend) + veg_ratio_max(1:nidx)*C2atmos(1:nidx)

    !! Compute emissions from landcover change: 
    !! .. i.e. conversion from [mol(C)/m^2(vegetated area)]during whole timestep to [kg(CO2)/m^2(grid box) s]

    CO2_emissions(:) = C2atmos(1:nidx) * veg_ratio_max(:) * molarMassCO2_kg/delta_time
                                   
    !! replace the old by the new landcover fractions

    landcover_fract_current(1:nidx,:) = landcover_fract_new(1:nidx,:)
    
  end subroutine do_landcover_change

  !! --- read_landcover_fractions() ----------------------------------------------------
  !!
  !! This routine reads (when necessary) the new landcover fractions and updates the module variable landcover_fract_new(:,:)
  !!
  subroutine read_landcover_fractions(fileName,varName,grid,domain,ntiles)
    use mo_jsbach_grid,      only: grid_type,domain_type
    use mo_exception,        only: finish, message,int2string
    use mo_netcdf,           only: FILE_INFO,IO_inq_dimid,IO_inq_dimlen,IO_inq_varid,io_get_var_double
    use mo_decomposition,    only: global_decomposition
    use mo_transpose,        only: scatter_gp
    use mo_io,               only: IO_open, IO_READ, IO_close
    use mo_jsbach,           only: jsbach_NewTimeStep
    use mo_jsbach_constants, only: molarMassCO2_kg
    USE mo_temp

    character(len=*),  intent(in)   :: fileName   !! Name of netcdf input file from which the new landcover fractions shall be read in
    character(len=*),  intent(in)   :: varName    !! Name of variable of the new landcover fractions in netcdf input file
    type(grid_type),   intent(in)   :: grid
    type(domain_type), intent(in)   :: domain
    integer,           intent(in)   :: ntiles   

    integer                     :: IO_file_id, IO_var_id, IO_dim_id
    type(FILE_INFO)             :: IO_file
    integer                     :: znlon, znlat,zntiles
    integer                     :: i,status


    if(jsbach_NewTimeStep) then
       if (p_parallel_io) then
          ! Open ini file
          call message('read_landcover_fractions','Reading new land surface fields from '//trim(filename))
          IO_file%opened = .false.
          call IO_open(trim(fileName), IO_file, IO_READ)
          IO_file_id = IO_file%file_id

          ! Check resolution
          call IO_inq_dimid  (IO_file_id, 'lat', IO_dim_id)
          call IO_inq_dimlen (IO_file_id, IO_dim_id, znlat)
          call IO_inq_dimid  (IO_file_id, 'lon', IO_dim_id)
          call IO_inq_dimlen (IO_file_id, IO_dim_id, znlon)
          call IO_inq_dimid  (IO_file_id, 'ntiles', IO_dim_id)
          call IO_inq_dimlen (IO_file_id, IO_dim_id, zntiles)

          if (znlon /= grid%nlon .or. znlat /= grid%nlat .or. zntiles /= ntiles) then
             call finish('read_landcover_fractions()', 'Unexpected grid resolution:'//int2string(znlon)//'x'//int2string(znlat))
          endif
          if (zntiles /= ntiles) then
             call finish('read_landcover_fractions()', 'Unexpected number of tiles:'//int2string(zntiles))
          endif


          !! allocate temporary memory

          allocate(zreal3d(grid%nlon,grid%nlat,ntiles),STAT=status)
          if(status .ne. 0) call finish('read_landcover_fractions()','Allocation failure (1)')

          !! read cover fractions

          call IO_inq_varid(IO_file_id, trim(varName), IO_var_id)
          call IO_get_var_double(IO_file_id, IO_var_id, zreal3d)

          !! Rescale cover fractions
          !! Reason: the cover fractions just read in are cover fractions with respect to the area of the
          !! whole grid cell. But what we need is the cover fraction with respect to the area of the grid 
          !! cell covered by vegetation. Therefore:
          !! Scale coverage to fractions of total vegetation land cover.
          ALLOCATE(zreal2d(grid%nlon, grid%nlat), STAT=status)
          if(status .ne. 0) call finish('read_landcover_fractions()','Allocation failure (2)')
          zreal2d = sum(zreal3d, DIM=3)   ! Total fraction of grid cell covered by vegetation
          if (any(grid%mask(:,:) .and. zreal2d(:,:) < 1.e-5_dp)) &
               call finish('read_landcover_fractions()', &
               'Land cover fractions from file '//trim(filename)//' inconsistent with global land-sea mask')
          do i=1,ntiles
             where (grid%mask) 
                zreal3d(:,:,i) = zreal3d(:,:,i) / zreal2d(:,:) !! scaling to grid cell area covered by vegetation
             end where
          end do
          deallocate(zreal2d)
          call IO_close(IO_file)
       endif !! end parallel_io

       !! Bring cover fractions to the other processors

       allocate(zreal2d(domain%ndim,domain%nblocks),STAT=status)
       if(status .ne. 0) call finish('read_landcover_fractions()','Allocation failure (3)')
       nullify(zreal2d_ptr)
       do i=1,ntiles
          if (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          call scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
          landcover_change%LCC_coverFract_target(:,i) = pack(zreal2d, MASK=domain%mask)
       enddo
    
       !! Free temporary memory
       if (p_parallel_io) deallocate(zreal3d)
       deallocate(zreal2d)

    end if !! jsbach_NewTimeStep

  end subroutine read_landcover_fractions
  
end module mo_cbal_landcover_change
