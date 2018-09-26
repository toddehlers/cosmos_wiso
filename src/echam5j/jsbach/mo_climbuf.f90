!-----------------------------------------------------------------------------------------------------------------------------------
!
! Module to calculate the climate variables
!
!-----------------------------------------------------------------------------------------------------------------------------------

MODULE mo_climbuf

!-----------------------------------------------------------------------------------------------------------------------------------
 
  USE mo_jsbach_grid,  ONLY: grid_type, domain_type
  USE mo_kind,         ONLY: dp 
  USE mo_jsbach,       ONLY: debug, missing_value
  USE mo_exception,    ONLY: finish, message
  USE mo_netCDF,       ONLY: FILE_INFO, NF_MAX_NAME
  USE mo_io,           ONLY: IO_open, IO_READ, IO_close
  USE mo_mpi,          ONLY: p_parallel_io
  USE mo_time_event,   ONLY: io_time_event

  IMPLICIT NONE  

  ! --- public subroutines ---

  PUBLIC :: init_climbuf             ! Allocates and initializes memory, sets standard parameters
  PUBLIC :: update_climbuf           ! Collects weather data for climatology
  PUBLIC :: climbuf_type

  TYPE climbuf_type                  ! contains the state variables
     REAL(dp), POINTER :: min_mmtemp20(:)
     REAL(dp), POINTER :: max_mmtemp20(:)
     REAL(dp), POINTER :: prev_year_gdd(:,:)
     REAL(dp), POINTER :: gdd_upper_tlim(:,:)
     REAL(dp), POINTER :: prev_year_npp(:,:)
     REAL(dp), POINTER :: ave_npp5(:,:)
     REAL(dp), POINTER :: prev_year_precip(:)
     REAL(dp), POINTER :: prev_year_soiltemp(:,:)
     REAL(dp), POINTER :: prev_year_soilmoist(:,:)
     REAL(dp), POINTER :: max_wind10(:)
     REAL(dp), POINTER :: prev_day_max_wind10(:)
     REAL(dp), POINTER :: rel_hum_air(:)
  end TYPE climbuf_type

  TYPE(climbuf_type), PUBLIC   :: climbuf


  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ===============================================================================

  PRIVATE                             ! Make ALL following objects private

  REAL(dp), PARAMETER :: persist_rel_hum = 0.95_dp   ! factor concerning the smoothing of relative air humidity
  REAL(dp), PARAMETER :: persist_wind10 = 0.9995_dp  ! factor concerning the smoothing of daily maximum wind speed

  INTEGER, SAVE       :: ntiles       ! number of tiles
  INTEGER, SAVE       :: nsoil        ! number of soil layers vertical
  ! --- private fields (state variables and fields for intra-module communication)
  INTEGER, SAVE       :: time_steps_per_day   ! number of time steps per day

  REAL(dp), POINTER :: temp_sum_month(:)
  REAL(dp), POINTER :: temp_sum_day(:)
  REAL(dp), POINTER :: surftemp_sum(:,:)
  REAL(dp), POINTER :: rel_soilmoist_sum(:,:)
  REAL(dp), POINTER :: max_wind10_act(:)
  REAL(dp), POINTER :: min_mmtemp_of_yr(:)
  REAL(dp), POINTER :: max_mmtemp_of_yr(:)
  REAL(dp), POINTER :: ann_npp_sum(:,:)
  REAL(dp), POINTER :: ann_precip_sum(:)
  REAL(dp), POINTER :: ann_gdd_sum(:,:)
  REAL(dp), POINTER :: gdd_upper_tlim_sum(:,:)


CONTAINS

  ! --- init_climbuf() -------------------------------------------------------------------------------------------------------------
  !
  ! This routine initializes the climate buffer module. It has to be called before the first time step.

  SUBROUTINE init_climbuf (grid, domain, notiles, isRestart, fileformat, nosoillayers, read_climbuf, climbuf_file_name, stream)

    USE mo_time_control,           ONLY: delta_time !! time step in seconds
    USE mo_memory_base,            ONLY: new_stream,default_stream_setting, &
                                         add =>add_stream_element
    USE mo_linked_list,            ONLY: t_stream, LAND, TILES
    USE mo_netCDF,                 ONLY: max_dim_name, io_inq_varid, io_get_var_double
    USE mo_temp                                                             ! Provided temporary arrays
    USE mo_decomposition,          ONLY: global_decomposition
    USE mo_transpose,              ONLY: scatter_gp


    TYPE(grid_type), INTENT(in)          :: grid
    TYPE(domain_type), INTENT(in)        :: domain
    INTEGER, INTENT(in)                  :: notiles, nosoillayers         
    LOGICAL, INTENT(in)                  :: isRestart
    INTEGER, INTENT(in)                  :: fileformat           ! output file format (grib/netcdf)
    LOGICAL, INTENT(in)                  :: read_climbuf
    CHARACTER(NF_MAX_NAME), INTENT(in)   :: climbuf_file_name
    TYPE(t_stream), POINTER, OPTIONAL    :: stream

    ! local variables

    TYPE(t_stream),POINTER      :: climbufStream

    INTEGER :: g_nland, l_nland
    INTEGER                     :: dim1p(1), dim1(1)
    INTEGER                     :: dim3p(2), dim3(2)
    CHARACTER(LEN=max_dim_name) :: dim3n(2), dim1n(1)   
    
    REAL(dp), POINTER, SAVE :: init_climbuf_data(:)
    TYPE(FILE_INFO)         :: climbuf_file    !! Input file for climate buffer
    INTEGER :: IO_file_id, IO_var_id



    IF (debug) CALL message('init_climbuf','Start initialization of climate buffer')
    
    g_nland        = grid%nland
    l_nland        = domain%nland
    nsoil          = nosoillayers
    ntiles         = notiles
    
    ! --- compute the number of time steps per day

    IF( MOD(86400,INT(delta_time)) .NE. 0) THEN ! For computing the mean day temperature (see below) it is assumed that each day ..
                                                ! is computed with a fixed number of time steps. Therefore the program is  ..
                                                ! stopped wheen day length is not an integer multiple of the time step.
       CALL finish("init_climbuf","ERROR: Day length is not an integer multiple of the time step!")
    ELSE 
       time_steps_per_day = 86400/INT(delta_time)
    END IF

    ! --- Open new stream for climate data
    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'climbuf', filetype=fileformat, interval=io_time_event(1,'days','first',0))
          ! Set default stream options
          CALL default_stream_setting(stream, table  = 198, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
        climbufStream => stream
    ELSE
       ! Add new stream
       CALL new_stream(climbufStream, 'climbuf', filetype=fileformat, interval=io_time_event(1,'days','first',0))
       ! Set default stream options
       CALL default_stream_setting(climbufStream, table  = 198, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF

    ! --- Define state variables as stream elements

    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)
 
    dim3p = (/ l_nland, ntiles /)
    dim3  = (/ g_nland, ntiles /)
    dim3n(1) = 'landpoint'
    dim3n(2) = 'tiles'

    CALL add(climbufStream,'temp_sum_day',         temp_sum_day, dim1p, dim1, dimnames=dim1n, laccu=.false., code=1, &
             lmiss=.TRUE., missval=missing_value,  lpost=.false.)
    CALL add(climbufStream,'temp_sum_month',     temp_sum_month, dim1p, dim1, dimnames=dim1n, laccu=.false., code=2, &
             lmiss=.TRUE., missval=missing_value,  lpost=.false.)
    CALL add(climbufStream,'min_mmtemp_of_yr', min_mmtemp_of_yr, dim1p, dim1, dimnames=dim1n, laccu=.false., code=3, &
             lmiss=.TRUE., missval=missing_value,  lpost=.false.)
    CALL add(climbufStream,'max_mmtemp_of_yr', max_mmtemp_of_yr, dim1p, dim1, dimnames=dim1n, laccu=.false., code=4, &
             lmiss=.TRUE., missval=missing_value,  lpost=.false.)
    CALL add(climbufStream,'min_mmtemp20', climbuf%min_mmtemp20, dim1p, dim1, dimnames=dim1n, laccu=.false., code=5, &
             lmiss=.TRUE., missval=missing_value,  lpost=.false.)
    CALL add(climbufStream,'max_mmtemp20', climbuf%max_mmtemp20, dim1p, dim1, dimnames=dim1n, laccu=.false., code=6, &
             lmiss=.TRUE., missval=missing_value,  lpost=.false.)
    CALL add(climbufStream,'surftemp_sum',         surftemp_sum, dim3p, dim3, dimnames=dim3n, laccu=.false., code=7, &
             lmiss=.TRUE., missval=missing_value,  lpost=.false.)
    CALL add(climbufStream,'prev_year_soiltemp',  &
                                     climbuf%prev_year_soiltemp, dim3p, dim3, dimnames=dim3n, laccu=.false., code=8, &
             lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(climbufStream,'rel_soilmoist_sum',rel_soilmoist_sum,dim3p, dim3, dimnames=dim3n, laccu=.false., code=9, &
             lmiss=.TRUE., missval=missing_value,  lpost=.false.)
    CALL add(climbufStream,'prev_year_soilmoist', &
                                    climbuf%prev_year_soilmoist, dim3p, dim3, dimnames=dim3n, laccu=.false., code=10, &
             lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(climbufStream,'ann_precip_sum',     ann_precip_sum, dim1p, dim1, dimnames=dim1n, laccu=.false., code=11, &
             lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(climbufStream,'prev_year_precip',    &
                                       climbuf%prev_year_precip, dim1p, dim1, dimnames=dim1n, laccu=.false., code=12, &
             lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(climbufStream,'ann_gdd_sum',           ann_gdd_sum, dim3p, dim3, dimnames=dim3n, laccu=.false., code=13, &
             lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(climbufStream,'prev_year_gdd',climbuf%prev_year_gdd,dim3p, dim3, dimnames=dim3n, laccu=.false., code=14, &
             lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(climbufStream,'gdd_upper_tlim_sum',gdd_upper_tlim_sum,dim3p,dim3,dimnames=dim3n, laccu=.false., code=15, &
             lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(climbufStream,'gdd_upper_tlim',climbuf%gdd_upper_tlim,dim3p,dim3,dimnames=dim3n, laccu=.false., code=16, &
             lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(climbufStream,'ann_npp_sum',           ann_npp_sum, dim3p, dim3, dimnames=dim3n, laccu=.false., code=17, &
             lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(climbufStream,'prev_year_npp',climbuf%prev_year_npp,dim3p, dim3, dimnames=dim3n, laccu=.false., code=18, &
             lmiss=.TRUE., missval=missing_value, lpost=.true.)
    CALL add(climbufStream,'ave_npp5',         climbuf%ave_npp5, dim3p, dim3, dimnames=dim3n, laccu=.false., code=19, &
             lmiss=.TRUE., missval=missing_value, lpost=.true.)
    CALL add(climbufStream,'max_wind10_act',     max_wind10_act, dim1p, dim1, dimnames=dim1n, laccu=.false., code=22, &
             lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(climbufStream,'max_wind10',     climbuf%max_wind10, dim1p, dim1, dimnames=dim1n, laccu=.false., code=23, &
             lmiss=.TRUE., missval=missing_value)
    CALL add(climbufStream,'prev_day_max_wind10', &
                                    climbuf%prev_day_max_wind10, dim1p, dim1, dimnames=dim1n, laccu=.false., code=20, &
             lmiss=.TRUE., missval=missing_value)
    CALL add(climbufStream,'rel_hum_air',   climbuf%rel_hum_air, dim1p, dim1, dimnames=dim1n, laccu=.false., code=21, &
             lmiss=.TRUE., missval=missing_value)

    ! If this is a restart run we can exit now since the model state variables are read from restart file

    IF (isRestart) THEN
       IF (debug) THEN
          CALL message('init_climbuf()','State variables read from restart file')
       END IF
       RETURN
    ENDIF

    ! --- initializations -------------------------------

    climbuf%rel_hum_air(:) = 50._dp

    IF (read_climbuf) THEN

       CALL message('init_climbuf','Reading initial climate information from file')

       ALLOCATE(init_climbuf_data(l_nland))
       
       ALLOCATE(zreal2d(domain%ndim,domain%nblocks))
       ALLOCATE(zzreal2d(grid%nlon,grid%nlat))

       IF (p_parallel_io) THEN
          climbuf_file%opened = .FALSE.
          CALL IO_open(TRIM(climbuf_file_name), climbuf_file, IO_READ)
          IO_file_id = climbuf_file%file_id

          CALL IO_inq_varid(IO_file_id, 'tsurf', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
       END IF

       NULLIFY(zreal2d_ptr)
       IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
       CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
       init_climbuf_data(:) = PACK(zreal2d, MASK=domain%mask)
       climbuf%min_mmtemp20 = init_climbuf_data
       climbuf%max_mmtemp20 = init_climbuf_data
       min_mmtemp_of_yr = init_climbuf_data
       max_mmtemp_of_yr = init_climbuf_data

       IF (p_parallel_io) THEN
          CALL IO_close(climbuf_file)
       END IF
       DEALLOCATE(zreal2d)
       DEALLOCATE(zzreal2d)

       IF (debug) THEN
          CALL message('init_climbuf()','Read climbuf_fields: Initialization of ClimateBuffer finished.')
       END IF

    ELSE
       min_mmtemp_of_yr=1000._dp
       max_mmtemp_of_yr=-1000._dp
       climbuf%min_mmtemp20 = 0._dp
       climbuf%max_mmtemp20 = 0._dp
    END IF

    IF (isRestart) RETURN

    ! -- Initialisation at beginning of an experiment

    temp_sum_day  = 0._dp
    temp_sum_month = 0._dp
    surftemp_sum = 0._dp
    climbuf%prev_year_soiltemp = 0._dp
    rel_soilmoist_sum = 0._dp
    climbuf%prev_year_soilmoist = 0._dp
    ann_precip_sum = 0._dp
    climbuf%prev_year_precip = 0._dp
    ann_gdd_sum = 0._dp
    climbuf%prev_year_gdd = 0._dp
    gdd_upper_tlim_sum = 0._dp
    climbuf%gdd_upper_tlim = 0._dp
    ann_npp_sum = 0._dp
    climbuf%prev_year_npp = 0._dp
    climbuf%ave_npp5 = 0._dp
    max_wind10_act = 0._dp
    climbuf%prev_day_max_wind10 = 0._dp
    climbuf%max_wind10 = 15._dp

    CALL message('init_climbuf','climate buffer initialised')

  END SUBROUTINE init_climbuf

!-----------------------------------------------------------------------------------------------------------------------------------
! Update the climate buffer
!
  SUBROUTINE update_climbuf(kidx, kidx0, kidx1, &
       gdd_base, upper_tlim, init_running_means, &
       temp_air, &
       temp_surf, wind10, &
       precip, rel_soil_water, &
       npp_rate, rel_hum_air)

!-----------------------------------------------------------------------------------------------------------------------------------
   USE mo_jsbach,                 ONLY: debug
   USE mo_exception,              ONLY: message
   USE mo_time_control,           ONLY: get_year_day, get_date_components, current_date, previous_date, start_date, &
                                        lstart, get_time_step, delta_time

! INTENT(IN)
!------------
   INTEGER,                   INTENT(in)  ::  kidx, kidx0, kidx1
   REAL(dp), DIMENSION(:),    INTENT(in)  ::  gdd_base                   ! base temperature to calculate growing degree days
   REAL(dp), DIMENSION(:),    INTENT(in)  ::  upper_tlim                 ! base temperature to calculate gdd_upper_tlim
   LOGICAL,                   INTENT(in)  ::  init_running_means         ! flag to initialize runing means 
                                                                         !     (after at least one complete year)
   REAL(dp), DIMENSION(:),    INTENT(in)  ::  temp_air                   ! air temperature of the lowest level [degC]
   REAL(dp), DIMENSION(:,:),  INTENT(in)  ::  temp_surf                  ! surface temperature [degC]
   REAL(dp), DIMENSION(:),    INTENT(in)  ::  wind10                     ! 10m wind speed [m/s]
   REAL(dp), DIMENSION(:),    INTENT(in)  ::  precip                     ! precipitation [m/s]
   REAL(dp), DIMENSION(:,:,:),INTENT(in)  ::  rel_soil_water             ! relative soil water content
   REAL(dp), DIMENSION(:,:),  INTENT(in)  ::  npp_rate                   ! NPP rate [mol(C)/m^2(canopy) s]
   REAL(dp), DIMENSION(:),    INTENT(in)  ::  rel_hum_air                ! relative air humidity

! LOCAL VARIABLES
!----------------
   INTEGER :: day_of_month,  day_of_month_at_prev_ts
   INTEGER :: month_of_year, month_of_year_at_prev_ts
   INTEGER :: year, year_at_prev_ts, start_year
   INTEGER :: istep, time_steps_in_last_year
   LOGICAL :: new_day, new_month, new_year, first_day
   REAL(dp):: month_mean_temp(kidx), prev_day_mean_temp(kidx)
   REAL(dp):: average_filling(kidx,ntiles)
   INTEGER :: pft, ns

!-----------------------------------------------------------------------------------------------------------------------------------
   IF (debug) CALL message ('update climbuf','starting routine')

   ! --- New day, new month, new year ?

   CALL get_date_components(current_date, YEAR=year, MONTH=month_of_year, DAY=day_of_month)
   CALL get_date_components(previous_date, YEAR=year_at_prev_ts, &
                            MONTH=month_of_year_at_prev_ts, DAY=day_of_month_at_prev_ts)
   CALL get_date_components(start_date, YEAR=start_year)
   istep=get_time_step()
   
   first_day = ( istep < time_steps_per_day)
   new_day   = ( day_of_month  /=  day_of_month_at_prev_ts)
   new_month = ( month_of_year /= month_of_year_at_prev_ts)
   new_year  = ( year          /=          year_at_prev_ts)

   time_steps_in_last_year = INT(get_year_day(previous_date)) * time_steps_per_day

   IF (debug .AND. new_day)   CALL message('update climbuf','First time step of new day')
   IF (debug .AND. new_month) CALL message('update climbuf','First time step of new month')
   IF (debug .AND. new_year)  CALL message('update climbuf','First time step of new year')

   !   -------------------
   ! -- Climate Variables
   !   -------------------

   ! --- Variables based on air temperature
   !
   IF (.NOT. new_day .OR. lstart) THEN
      temp_sum_day(kidx0:kidx1) = temp_sum_day(kidx0:kidx1) + temp_air(:)
   ELSE
      IF (first_day) THEN
         prev_day_mean_temp(:) = temp_sum_day(kidx0:kidx1) / istep
      ELSE
         prev_day_mean_temp(:) = temp_sum_day(kidx0:kidx1) / time_steps_per_day
      ENDIF
      temp_sum_day(kidx0:kidx1) = temp_air(:)

   ! -- Growing Degree Days

      DO pft = 1,ntiles
         ann_gdd_sum(kidx0:kidx1,pft) = &
              ann_gdd_sum(kidx0:kidx1,pft) + MAX(prev_day_mean_temp(:) - gdd_base(pft), 0._dp)
         gdd_upper_tlim_sum(kidx0:kidx1,pft) = &
              gdd_upper_tlim_sum(kidx0:kidx1,pft) + MAX(prev_day_mean_temp(:) - upper_tlim(pft), 0._dp)
      END DO

      IF (new_year) THEN

      ! --- provide previous year GDD and reset for current year
         DO pft = 1, ntiles
            climbuf%prev_year_gdd(kidx0:kidx1,pft) = ann_gdd_sum(kidx0:kidx1,pft)
            ann_gdd_sum(kidx0:kidx1,pft) =  0._dp
            climbuf%gdd_upper_tlim(kidx0:kidx1,pft) = gdd_upper_tlim_sum(kidx0:kidx1,pft)
            gdd_upper_tlim_sum(kidx0:kidx1,pft) = 0._dp
         END DO
      ENDIF

   !
   ! --- Minimum/Maximum monthly mean temperature (20 year climatology) 
   !
      temp_sum_month(kidx0:kidx1) = temp_sum_month(kidx0:kidx1) + prev_day_mean_temp(:)

      IF (new_month) THEN

      ! --- provide previous month mean temperature and reset for current month
         month_mean_temp(:) = temp_sum_month(kidx0:kidx1) / day_of_month_at_prev_ts      
         temp_sum_month(kidx0:kidx1) = 0._dp

      ! --- calculate coldest/warmest month of the year so far
         min_mmtemp_of_yr(kidx0:kidx1) = min(min_mmtemp_of_yr(kidx0:kidx1), month_mean_temp(:))
         max_mmtemp_of_yr(kidx0:kidx1) = max(max_mmtemp_of_yr(kidx0:kidx1), month_mean_temp(:))

         IF (new_year) THEN

         ! --- bulid kind of 20yr running mean of coldest and warmest monthly temperature
            IF (init_running_means) THEN
               IF (debug) CALL message ('update climbuf','initialisation of mmtemp20')
               climbuf%min_mmtemp20(kidx0:kidx1) = min_mmtemp_of_yr(kidx0:kidx1)
               climbuf%max_mmtemp20(kidx0:kidx1) = max_mmtemp_of_yr(kidx0:kidx1)
            ELSE
               climbuf%min_mmtemp20(kidx0:kidx1) = (climbuf%min_mmtemp20(kidx0:kidx1)*19._dp + min_mmtemp_of_yr(kidx0:kidx1)) &
                                                    /20._dp
               climbuf%max_mmtemp20(kidx0:kidx1) = (climbuf%max_mmtemp20(kidx0:kidx1)*19._dp + max_mmtemp_of_yr(kidx0:kidx1)) &
                                                    /20._dp
            ENDIF

         ! --- reset for current year
            min_mmtemp_of_yr(kidx0:kidx1) = 1000._dp
            max_mmtemp_of_yr(kidx0:kidx1) = -1000._dp

         END IF

      END IF
   END IF

   ! --- Maximum daily wind speed and relative air humidity smoothed in time

   IF (.NOT. new_day .OR. lstart) THEN
      WHERE (wind10(:) > max_wind10_act(kidx0:kidx1))
         max_wind10_act(kidx0:kidx1) = wind10(:)
      END WHERE
   ELSE
      climbuf%prev_day_max_wind10(kidx0:kidx1) = max_wind10_act(kidx0:kidx1)
      climbuf%max_wind10(kidx0:kidx1) = climbuf%max_wind10(kidx0:kidx1) * persist_wind10 + &
                                        max_wind10_act(kidx0:kidx1) * (1._dp - persist_wind10)
      max_wind10_act(kidx0:kidx1) = wind10(:)
      climbuf%rel_hum_air(kidx0:kidx1) = climbuf%rel_hum_air(kidx0:kidx1) * persist_rel_hum + &
                                         MIN(rel_hum_air(1:kidx),100._dp) * (1._dp - persist_rel_hum)
   END IF

   ! Average soil moisture integrated column
   average_filling(:,:) = 0._dp
   DO ns = 1,nsoil
      average_filling(:,:) = average_filling(:,:) + rel_soil_water(:,ns,:)
   END DO
   average_filling(:,:) = average_filling(:,:) / REAL(nsoil,dp)

   IF (.NOT. new_year .OR. lstart) THEN
      ann_npp_sum(kidx0:kidx1,:)  = ann_npp_sum(kidx0:kidx1,:) + NPP_Rate(1:kidx,:) * delta_time
      ann_precip_sum(kidx0:kidx1) = ann_precip_sum(kidx0:kidx1) + precip(1:kidx) * delta_time
      surftemp_sum(kidx0:kidx1,:)  = surftemp_sum(kidx0:kidx1,:) + temp_surf(:,:)
      rel_soilmoist_sum(kidx0:kidx1,:) = rel_soilmoist_sum(kidx0:kidx1,:) + average_filling(:,:)
   ELSE

      ! --- provide previous year NPP and reset for current year

      climbuf%prev_year_npp(kidx0:kidx1,:) = ann_npp_sum(kidx0:kidx1,:)
      ann_npp_sum(kidx0:kidx1,:) = NPP_Rate(1:kidx,:) * delta_time

      ! --- bulid kind of 5yr running mean of NPP
 
      IF (init_running_means) THEN
         climbuf%ave_npp5(kidx0:kidx1,:) = climbuf%prev_year_npp(kidx0:kidx1,:)
      ELSE
         climbuf%ave_npp5(kidx0:kidx1,:) = (climbuf%ave_npp5(kidx0:kidx1,:)*4._dp + climbuf%prev_year_npp(kidx0:kidx1,:)) / 5._dp
      ENDIF

      ! --- provide previous year precipitation and reset for current year

      climbuf%prev_year_precip(kidx0:kidx1) = ann_precip_sum(kidx0:kidx1)
      ann_precip_sum(kidx0:kidx1) = precip(:) * delta_time

      ! --- provide previous year soil temperature and reset for current year

      climbuf%prev_year_soiltemp(kidx0:kidx1,:) = surftemp_sum(kidx0:kidx1,:) / time_steps_in_last_year
      surftemp_sum(kidx0:kidx1,:) = temp_surf(:,:)

      ! --- provide previous year soil moisture and reset for current year

      climbuf%prev_year_soilmoist(kidx0:kidx1,:) = rel_soilmoist_sum(kidx0:kidx1,:) / time_steps_in_last_year
      rel_soilmoist_sum(kidx0:kidx1,:) = average_filling(:,:)

   END IF

   IF (debug) CALL message ('update climbuf','end of routine')

  END SUBROUTINE update_climbuf

END MODULE mo_climbuf
