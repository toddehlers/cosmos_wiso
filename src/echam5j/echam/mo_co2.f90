MODULE mo_co2

  ! L. Kornblueh, MPI, November 2004, implementation
  ! R. Schnur, MPI, November 2004, added some diagnostics and new stream variables for land/ocean fluxes
  ! R. Schnur, MPI, December 2004, added routines to use anthropogenic CO2 emissions

  USE mo_kind,          ONLY: dp
  USE mo_tracer,        ONLY: new_submodel, new_tracer, get_tracer, &
                              OFF, ON, RESTART, CONSTANT
  USE mo_radiation,     ONLY: co2mmr, ico2
  USE mo_memory_base,   ONLY: new_stream, default_stream_setting,   &
                              add_stream_element,                   &
                              GRIB, NETCDF, GAUSSIAN
  USE mo_linked_list,   ONLY: t_stream
  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_mpi,           ONLY: p_parallel_io, p_io, p_bcast
  USE mo_exception,     ONLY: message_text, message, finish, real2string, int2string
  USE mo_control,       ONLY: lcouple, ico2_emis, lco2, lco2_flxcor, lco2_2perc, &
                              lco2_nudge, co2_nudge_vmr, co2_nudge_tau
  USE mo_constants,     ONLY: g, amco2, amd   ! accelaration of gravity [m s-2], molecular weight of CO2 [g mol-1],
                                              ! molecular weight of dry air [g mol-1]
  USE mo_filename,      ONLY: out_filetype

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_submodel_co2, init_co2, reference_co2, diag_co2, read_co2_emission, &
            co2_emissions, co2_flux_correction
#ifdef __cpl_co2
  PUBLIC :: co2_flux_atmosphere_ocean
#endif

  REAL(dp), POINTER, PUBLIC :: co2(:,:,:)
  REAL(dp), POINTER, PUBLIC :: co2m1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: co2atmos(:,:)  ! CO2 concentration in lowest level, accumulated for coupling (collect.f90)
  REAL(dp), POINTER, PUBLIC :: co2flux_cpl(:,:) ! upward ocean CO2 flux, accumulated for coupling (collect.f90)
#ifdef __cpl_co2
  REAL(dp), POINTER, PUBLIC :: co2trans(:,:)  ! CO2 ocean-atmosphere transfer velocity (not including wind speed)
  REAL(dp), POINTER, PUBLIC :: co2ocean(:,:)  ! CO2 concentration of the ocean
#endif
  REAL(dp), POINTER, PUBLIC :: co2_flux(:,:)
  REAL(dp), POINTER, PUBLIC :: co2_flux_land(:,:)
  REAL(dp), POINTER, PUBLIC :: co2_flux_ocean(:,:)
  REAL(dp), POINTER, PUBLIC :: co2_flux_acc(:,:)
  REAL(dp), POINTER, PUBLIC :: co2_flux_land_acc(:,:)
  REAL(dp), POINTER, PUBLIC :: co2_flux_ocean_acc(:,:)
  REAL(dp), POINTER, PUBLIC :: co2_flux_anthro_acc(:,:)
  REAL(dp), POINTER, PUBLIC :: co2_burden(:,:)
  REAL(dp), POINTER, PUBLIC :: co2_burden_old(:,:)
  REAL(dp), POINTER, PUBLIC :: co2_burden_acc(:,:)
  REAL(dp), POINTER, PUBLIC :: co2_flux_corr_acc(:,:)        ! CO2 flux correction, accumulated for output time step
  REAL(dp), POINTER, PUBLIC :: co2_burden_corr_acc1(:,:)     ! CO2 burden correction, accumulated for daily correction time step
  ! Note that co2_flux_corr_acc1 and co2_flux_corr_acc1_old are globally constant, but are kept in spatial fields because they
  ! need to be written/read to/from the restart files.
  REAL(dp), POINTER, PUBLIC :: co2_flux_corr_acc1(:,:)       ! CO2 flux correction, accumulated for daily correction time step
  REAL(dp), POINTER, PUBLIC :: co2_flux_corr_acc1_old(:,:)   ! CO2 flux correction, accumulated for daily correction time step
  REAL(dp), POINTER, PUBLIC :: co2_burden_corr_acc2(:,:)     ! CO2 burden correction in tendencies, accumulated for output time step
  REAL(dp), POINTER, PUBLIC :: co2_burden_nudge(:,:)         ! Correction term by CO2 nudging
  REAL(dp), POINTER, PUBLIC :: co2_burden_nudge_acc(:,:)     ! Correction term by CO2 nudging, accumulated for output time step
  REAL(dp), POINTER, PUBLIC :: co2_emission(:,:)
  REAL(dp), POINTER, PUBLIC :: co2_emission_acc(:,:)
  REAL(dp), POINTER, PUBLIC :: co2_flux_total(:,:)      ! Flux + flux correction + emissions
!!$  REAL(dp), POINTER         :: co2_emission_in(:,:), &   ! CO2 emission field from netCDF file for current year
!!$                               co2_emission_inm1, &  !                 "                       previous year
!!$                               co2_emission_inp1     !                 "                       following year
  REAL(dp), PUBLIC :: co2_nudge_mmr     = 0._dp
  REAL(dp), PUBLIC :: co2_nudge_relax   = 0._dp

  INTEGER :: emis_no_years = -1
  REAL(dp) :: emis_base_year
  REAL(dp), ALLOCATABLE :: emis_years(:)
  INTEGER :: emis_id = -1  ! NetCDF file id

  INTEGER, PUBLIC :: ico2idx
  TYPE(t_stream), POINTER :: xtco2

  LOGICAL, SAVE, PUBLIC :: l_co2flxcorr = .FALSE.

CONTAINS

  SUBROUTINE init_submodel_co2

    INTEGER :: co2_submodel_id
!    INTEGER :: i, j, ierr 
    INTEGER :: nlon, nglat, lnlon, lnglat 

    LOGICAL :: lpost = .TRUE.
    LOGICAL :: lrerun = .TRUE.

    IF (ico2 == 1 .AND. .NOT. lco2) THEN
       lco2 = .TRUE.
       CALL message('mo_co2','Switching on CO2 tracer module (since ico2=1)')
    END IF

    IF (.NOT. lco2) THEN
!!$       lpost = .FALSE.
!!$       lrerun = .FALSE.
       CALL message('mo_co2', 'Not using CO2 tracer module')
    END IF

    IF (lco2 .AND. .NOT. lco2_flxcor) THEN
       CALL message('mo_co2', 'CO2 flux correction turned off by namelist switch!')
    END IF

    IF (lco2 .AND. lco2_nudge .AND. lco2_flxcor) THEN
       CALL message('mo_co2', 'CO2 flux correction turned off while doing nudging')
       lco2_flxcor = .FALSE.
    END IF

    IF (lco2) THEN
       IF (lco2_2perc) THEN
          CALL message('mo_co2', 'CO2 2% limitation in physc switched on')
       ELSE
          CALL message('mo_co2', 'CO2 2% limitation in physc switched off')
       END IF
    END IF

    nlon   = lc%nlon
    nglat  = lc%nlat
    lnlon  = lc%nproma
    lnglat = lc%ngpblks

    CALL new_stream (xtco2, 'co2', filetype = out_filetype, &
         post_suf = '_co2', rest_suf = '_co2') 

    CALL default_stream_setting (xtco2,           &
         table = 191, bits = 16, repr = GAUSSIAN, &
         lrerun = lrerun, lpost = lpost)

    IF (lco2) CALL new_submodel ('CO2', co2_submodel_id)

    ! add extra fields required for usage with interactive carbon cycle

    ! Instantaneous CO2 flux from the surface (ocean+land) into the atmosphere
    ! This will be updated every time step for the land surface, but only
    ! at each coupling time step for the ocean
    IF (lcouple) THEN
       CALL add_stream_element (xtco2, 'co2atmos', co2atmos, &
            (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
            code = 16, longname = 'CO2 concentration at surface for ocean coupling',            & 
            units = 'kg kg-1', laccu=.FALSE., lpost=.FALSE., &    ! is accumulated in collect at coupling interval
            lrerun=lrerun, contnorest=.TRUE.) 
    END IF
    CALL add_stream_element (xtco2, 'co2_flux_inst', co2_flux, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
         code = 2, bits=24,                               & 
         longname = 'upward surface CO2 flux',            & 
         units = 'kg m-2 s-1', laccu=.FALSE.) 
    CALL add_stream_element (xtco2, 'co2_flux_l_inst', co2_flux_land, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
         code = 3, bits=24,                               & 
         longname = 'upward land CO2 flux',            & 
         units = 'kg m-2 s-1', laccu=.FALSE., lpost=.FALSE.) 
    CALL add_stream_element (xtco2, 'co2_flux_o_inst', co2_flux_ocean, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
         code = 4, bits=24,                               & 
         longname = 'upward ocean CO2 flux',            & 
         units = 'kg m-2 s-1', laccu=.FALSE., lpost=.FALSE.)
   CALL add_stream_element (xtco2, 'co2flux_cpl', co2flux_cpl,        &
         (/ lnlon, lnglat /), (/ nlon, nglat /),          &
         code = 19, & 
         longname = 'upward ocean CO2 flux for ocean coupling', & 
         units = 'kg m-2 s-1', laccu=.FALSE., lpost=.FALSE., & ! is accumulated in collect at coupling interval
         lrerun=lrerun, contnorest=.TRUE.)
#ifdef __cpl_co2
    CALL add_stream_element (xtco2, 'co2trans', co2trans,  &
         (/ lnlon, lnglat /), (/ nlon, nglat /),           & 
         code = 18, bits=24,                               & 
         longname = 'transfer velocity of ocean/atmosphere CO2 flux',    & 
         units = '10-9 mol s m-4', laccu=.FALSE., lpost=.TRUE.)
    CALL add_stream_element (xtco2, 'co2ocean', co2ocean,  &
         (/ lnlon, lnglat /), (/ nlon, nglat /),           & 
         code = 17,bits=24,                                & 
         longname = 'pco2 of surface ocean water',         & 
         units = 'ppm CO2', laccu=.FALSE., lpost=.TRUE.)
#endif
    CALL add_stream_element (xtco2, 'co2_flux_total_inst', co2_flux_total, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),           & 
         code = 15, bits=24,                               & 
         longname = 'total upward surface CO2 flux',            & 
         units = 'kg m-2 s-1', laccu=.FALSE., lpost=lpost .AND. lco2, contnorest=.TRUE.) 

    CALL add_stream_element (xtco2, 'co2_burden_inst', co2_burden, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                  &
         code=14, &
         longname = 'CO2 content',                                &
         units = 'kg m-2', laccu=.FALSE., lrerun=lrerun, contnorest=.TRUE.) 

    CALL add_stream_element (xtco2, 'co2_burden_old', co2_burden_old, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                  &
         longname = 'old CO2 content',                                &
         units = 'kg m-2', laccu=.FALSE., lpost=.FALSE., lrerun=lrerun, contnorest=.TRUE.) 

    IF (lco2_nudge) THEN
       CALL add_stream_element(xtco2, 'co2_burden_nudge_inst', co2_burden_nudge, &
            (/ lnlon, lnglat /), (/ nlon, nglat /),                  &
            longname = 'CO2 nudging correction',                                &
            units = 'kg m-2', laccu=.FALSE., lpost=.FALSE., lrerun=lrerun, contnorest=.TRUE.) 
    END IF
         
    IF (ico2_emis > 0) THEN
       CALL add_stream_element (xtco2, 'co2_emis_inst', co2_emission, &
            (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
            code = 12, bits=24, & 
            longname = 'upward CO2 emission',            & 
            units = 'kg m-2 s-1', laccu=.FALSE., lpost=.FALSE., lrerun=lrerun, contnorest=.TRUE.)
       ! lrerun must be true because emissions are read only at the start of a new year

!!$       CALL add_stream_element(xtco2, 'co2_emis_in', co2_emission_in, &
!!$            (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
!!$            laccu=.FALSE., lpost=.FALSE., lrerun=.FALSE.)
!!$       CALL add_stream_element(xtco2, 'co2_emis_inm1', co2_emission_inm1, &
!!$            (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
!!$            laccu=.FALSE., lpost=.FALSE., lrerun=.FALSE.)
!!$       CALL add_stream_element(xtco2, 'co2_emis_inp1', co2_emission_inp1, &
!!$            (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
!!$            laccu=.FALSE., lpost=.FALSE., lrerun=.FALSE.)
    END IF

    ! Accumulated CO2 flux
    CALL add_stream_element (xtco2, 'co2_flux', co2_flux_acc, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                  & 
         code = 5, bits=24,                                       & 
         longname = 'upward surface CO2 flux (acc.)',             & 
         units = 'kg m-2 s-1', laccu=.TRUE., lrerun=lrerun, contnorest=.TRUE.) 
    CALL add_stream_element (xtco2, 'co2_flx_land', co2_flux_land_acc, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
         code = 6, bits=24,                               & 
         longname = 'upward land CO2 flux (acc.)',            & 
         units = 'kg m-2 s-1', laccu=.TRUE., lrerun=lrerun, contnorest=.TRUE.) 
    CALL add_stream_element (xtco2, 'co2_flx_ocean', co2_flux_ocean_acc, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
         code = 7, bits=24,                               & 
         longname = 'upward ocean CO2 flux (acc.)',            & 
         units = 'kg m-2 s-1', laccu=.TRUE., lpost=lpost.AND.lcouple,lrerun=lrerun, contnorest=.TRUE.) 
    IF (ico2_emis > 0) THEN
       CALL add_stream_element (xtco2, 'co2_emis', co2_emission_acc, &
            (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
            code = 11, bits=24,                              & 
            longname = 'upward CO2 emission (acc.)',            & 
            units = 'kg m-2 s-1', laccu=.TRUE., lrerun=lrerun, contnorest=.TRUE.) 
    END IF

    CALL add_stream_element (xtco2, 'co2_flx_anthro', co2_flux_anthro_acc, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                           &
         code=20, bits=24,                                                 &
         longname = 'upward anthropogenic CO2 flux (acc.)',                &
         units = 'kg m-2 s-1', laccu=.TRUE., lpost=lpost.AND.lcouple, lrerun=lrerun, contnorest=.TRUE.) 

    ! Accumulated vertical integral of CO2 concentration
    CALL add_stream_element (xtco2, 'co2_burden', co2_burden_acc, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                  &
         code = 8,                                                &
         longname = 'CO2 content (acc.)',                                &
         units = 'kg m-2', laccu=.TRUE., lpost=.TRUE., lrerun=lrerun, contnorest=.TRUE.) 
    CALL add_stream_element (xtco2, 'co2_flux_corr', co2_flux_corr_acc, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
         code = 10,                                        &
         longname = 'CO2 flux correction (acc.)',       & 
         units = 'kg m-2 s-1', laccu=.TRUE., lpost=lpost .AND. lco2, &
         lrerun=lrerun .AND. lco2, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_burden_corr_acc1', co2_burden_corr_acc1, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
         units = 'kg m-2', laccu=.FALSE., lpost=.FALSE., &   ! Is accumulated manually for daily correction time step
         lrerun=lrerun .AND. lco2, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_flux_corr_acc1', co2_flux_corr_acc1, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
         units = 'kg m-2', laccu=.FALSE., lpost=.FALSE., &
         lrerun=lrerun .AND. lco2, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_flux_corr_acc1_old', co2_flux_corr_acc1_old, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
         units = 'kg m-2', laccu=.FALSE., lpost=.FALSE., &
         lrerun=lrerun .AND. lco2, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_burden_corr_acc2', co2_burden_corr_acc2, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),          & 
         code = 9,                                        &
         longname = 'CO2 correction in tendencies (acc.)',       & 
         units = 'kg m-2', laccu=.TRUE., lpost=lpost .AND. lco2, &
         lrerun=.FALSE., contnorest=.TRUE.)
    IF (lco2_nudge) THEN
       CALL add_stream_element(xtco2, 'co2_burden_nudge_acc', co2_burden_nudge_acc, &
            (/ lnlon, lnglat /), (/ nlon, nglat /),                  &
            code = 13, &
            longname = 'CO2 nudging correction (acc.)',                                &
            units = 'kg m-2', laccu=.TRUE., lpost=lpost .AND. lco2, &
            lrerun=lrerun .AND. lco2, contnorest=.TRUE.) 
    END IF

    IF (lcouple) co2atmos = 0._dp
    co2flux_cpl = 0._dp
    co2_flux = 0._dp
    co2_flux_land = 0._dp
    co2_flux_ocean = 0._dp
    co2_flux_acc = 0._dp
    co2_flux_land_acc = 0._dp
    co2_flux_ocean_acc = 0._dp
    co2_flux_anthro_acc = 0._dp
    co2_flux_total = 0._dp
    co2_burden = -1._dp
    co2_burden_old = 0._dp
    co2_burden_acc = 0._dp
    co2_flux_corr_acc = 0._dp
    co2_burden_corr_acc1 = 0._dp
    co2_flux_corr_acc1 = 0._dp
    co2_flux_corr_acc1_old = 0._dp
    co2_burden_corr_acc2 = 0._dp
    IF (lco2_nudge) THEN
       co2_burden_nudge = 0._dp
       co2_burden_nudge_acc = 0._dp
    END IF
    IF (ico2_emis > 0) THEN
       co2_emission = 0._dp
       co2_emission_acc = 0._dp
    END IF

  END SUBROUTINE init_submodel_co2

  SUBROUTINE init_co2

    USE mo_time_control, ONLY: time_step_len, lresume

    INTEGER :: ierr

    IF (lco2) THEN
       ! add 3d tracer field reference from mo_memory_gl
       CALL new_tracer ('CO2',                        &
            'CO2', idx = ico2idx,                     &
            longname = 'mass fraction of CO2 in air', &
            units = 'kg kg-1',                        &
            code = 1, table = 191, bits = 16,         &
            ninit = RESTART+CONSTANT, vini = co2mmr,  &
            nwrite = ON, nrerun = ON, nsemis = OFF,   &
            ierr = ierr)
       IF (p_parallel_io .AND. .NOT. lresume) THEN
          WRITE (message_text,'(a,e12.6)') &
               'Initialized global CO2 mass mixing ratio to ', co2mmr
          CALL message('init_co2', TRIM(message_text))
       END IF

       IF (lco2_nudge) THEN

          co2_nudge_mmr=co2_nudge_vmr*amco2/amd
          CALL message('init_co2', &
                       'Nudging CO2 concentration to target vmr/mmr '//real2string(co2_nudge_vmr)//&
                                                                 ', '//real2string(co2_nudge_mmr))
          IF (co2_nudge_tau < EPSILON(1._dp)) &
               CALL finish('init_co2', 'Please specify relaxation time for CO2 nudging')
          co2_nudge_relax = MAX(0._dp, MIN(time_step_len/co2_nudge_tau,1._dp))
          CALL message('init_co2', &
                       'Relaxation time for CO2 nudging: '//real2string(co2_nudge_tau))

       END IF

    ELSE
       ! make sure co2m1 is allocated even though it's never used (otherwise, a bus error will
       ! occur when subroutines are called with co2m1 in the parameter list (e.g. radiation)
       CALL add_stream_element(xtco2, 'mass fraction of CO2 in air', co2m1, &
            (/ lc%nproma, lc%nlev, lc%ngpblks /), (/ lc%nlon, lc%nlev, lc%nlat /), &
            lrerun=.FALSE., lpost=.FALSE.)
       co2m1 = co2mmr
       ico2idx = -1
    END IF

  END SUBROUTINE init_co2

  SUBROUTINE reference_co2
    
    INTEGER :: ierr

    IF (.NOT. lco2) RETURN

    CALL get_tracer ('CO2', idx = ico2idx, pxt = co2, pxtm1 = co2m1, &
         ierr = ierr) 

  END SUBROUTINE reference_co2

  SUBROUTINE read_co2_emission

    ! Read annual carbon emissions from netCDF file and convert to CO2 

    USE mo_time_conversion, ONLY: time_native, &
                                  TC_get, TC_convert, &
                                  year_len, day_in_year, &
                                  print_date
    USE mo_time_control, ONLY: previous_date, current_date, get_time_step, lstart, lresume
    USE mo_transpose, ONLY: scatter_gp
    USE mo_decomposition, ONLY: gl_dc => global_decomposition
    USE mo_control, ONLY: nlon, ngl
    USE mo_netcdfstream, ONLY: nf

    REAL(dp), PARAMETER :: amc = 12.0107_dp   ! Molecular weight of C [g mol-1]

    INTEGER :: nvarid, ndimid
    INTEGER :: zcount(3), zstart(3)
    REAL(dp), POINTER :: ztmp1(:,:)

!    INTEGER :: i
    TYPE(time_native) :: emis_date
!    INTEGER :: iyear, iyearm, iyearp
    INTEGER :: iyear
    INTEGER :: yr, mo, dy, hr, mn, se, yrm
!    INTEGER :: day, second

    INCLUDE 'netcdf.inc'

    IF (.NOT. lco2 .OR. ico2_emis < 1) THEN
       !CALL message('mo_co2','Not reading CO2 emissions')
       RETURN
    END IF

    IF (emis_no_years < 0) THEN

       IF (p_parallel_io) THEN
          ! Select scenario
          IF (ico2_emis == 1) CALL nf(nf_open('carbon_emissions_A1B.nc', nf_nowrite, emis_id))
          IF (ico2_emis == 2) CALL nf(nf_open('carbon_emissions_B1.nc', nf_nowrite, emis_id))
          IF (ico2_emis == 3) CALL nf(nf_open('carbon_emissions_A2.nc', nf_nowrite, emis_id))

          CALL nf(nf_inq_dimid(emis_id, 'time', ndimid))
          CALL nf(nf_inq_dimlen(emis_id, ndimid, emis_no_years))
       END IF
       CALL p_bcast(emis_no_years, p_io)
       ALLOCATE(emis_years(emis_no_years))
       IF (p_parallel_io) THEN
          CALL nf(nf_inq_varid(emis_id, 'time', nvarid))
          CALL nf(nf_get_var_double(emis_id, nvarid, emis_years))
       ENDIF

       CALL p_bcast(emis_years, p_io)
       emis_base_year = emis_years(1)
       
    END IF

    CALL TC_convert(previous_date, emis_date) 
    CALL TC_get (emis_date, yrm, mo, dy, hr, mn, se)
    CALL TC_convert(current_date, emis_date) 
    CALL TC_get (emis_date, yr, mo, dy, hr, mn, se)

    IF (yrm == yr) RETURN      ! same year, nothing to do
!!$    IF (yrm == yr .AND. .NOT. (lstart .OR. lresume)) RETURN      ! same year, nothing to do

!!$    ! We're either at the beginning of a new year or this is an initial or resumed run.
!!$    ! The emission fields for yr-1,yr,yr+1 have therefore to be read from the netCDF file.

    iyear = yr-emis_base_year+1   ! set right index to access in emission fields

    IF (yr < emis_base_year .OR. iyear > emis_no_years) THEN
       CALL finish('read_co2_emission', 'Year not present in emissions file')
    END IF

    co2_emission = 0._dp

    IF (p_parallel_io) THEN
       ALLOCATE(ztmp1(nlon,ngl))
       zstart(:) = (/ 1,1,iyear /)
       zcount(:) = (/ nlon,ngl,1 /)
       CALL nf(nf_inq_varid(emis_id, 'carbon_emission', nvarid))
       CALL nf(nf_get_vara_double(emis_id, nvarid, zstart, zcount, ztmp1))
    END IF
    CALL scatter_gp(ztmp1, co2_emission, gl_dc)
!!$    CALL scatter_gp(ztmp1, co2_emission_in, gl_dc)
!!$
!!$    IF (yr > emis_base_year) THEN
!!$       IF (p_parallel_io) THEN
!!$          zstart(:) = (/ 1,1,iyear-1 /)
!!$          CALL nf(nf_inq_varid(emis_id, 'carbon_emission', nvarid))
!!$          CALL nf(nf_get_vara_double(emis_id, nvarid, zstart, zcount, ztmp1))
!!$       END IF
!!$       CALL scatter_gp(ztmp1, co2_emission_inm1, gl_dc)
!!$    ELSE
!!$       co2_emission_inm1 = co2_emission_in
!!$    END IF
!!$
!!$    IF (iyear < emis_no_years) THEN
!!$       IF (p_parallel_io) THEN
!!$          zstart(:) = (/ 1,1,iyear-1 /)
!!$          CALL nf(nf_inq_varid(emis_id, 'carbon_emission', nvarid))
!!$          CALL nf(nf_get_vara_double(emis_id, nvarid, zstart, zcount, ztmp1))
!!$       END IF
!!$       CALL scatter_gp(ztmp1, co2_emission_inp1, gl_dc)
!!$    ELSE
!!$       co2_emission_inp1 = co2_emission_in
!!$    END IF

    IF (p_parallel_io) THEN
       DEALLOCATE(ztmp1)
!!$       WRITE(message_text,*) 'nstep: ', get_time_step(), ' Read carbon emissions for years: ', yr-1,yr,yr+1
       WRITE(message_text,*) 'nstep: ', get_time_step(), ' Read carbon emissions for year: ', yr
       CALL message('', TRIM(message_text))
    END IF
    
    ! Convert carbon emission in g m2-1 s-1 to CO2 emission in kg m2-1 s-1
    co2_emission   = co2_emission   * amco2 / (amc * 1000._dp)
!!$    co2_emission_in   = co2_emission_in   * amco2 / (amc * 1000._dp)
!!$    co2_emission_inm1 = co2_emission_inm1 * amco2 / (amc * 1000._dp)
!!$    co2_emission_inp1 = co2_emission_inp1 * amco2 / (amc * 1000._dp)

  END SUBROUTINE read_co2_emission

  SUBROUTINE co2_emissions(kproma, kbdim, ktrac, pxtems, jrow)

    INTEGER, INTENT(in) :: kproma, kbdim, ktrac, jrow
    REAL(dp), INTENT(inout) :: pxtems(kbdim,ktrac)

    IF (.NOT. lco2) RETURN

    IF (ico2_emis > 0) THEN
       pxtems(1:kproma, ico2idx) = pxtems(1:kproma,ico2idx) + co2_emission(1:kproma,jrow)
    END IF

    co2_flux_total(1:kproma,jrow) = pxtems(1:kproma,ico2idx)
!!$    co2_flux_total(1:kproma,jrow) = 0._dp

  END SUBROUTINE co2_emissions

#ifdef __cpl_co2
  SUBROUTINE co2_flux_atmosphere_ocean(jrow, kproma, klev)

    ! Calculation of the CO2 flux between atmosphere and ocean [kg m-2 s-1]

    USE mo_memory_g3b,   ONLY: wind10w

    INTEGER, INTENT(in) :: jrow, kproma, klev

    ! The range (Max-Min) of co2trans is smaller than 1.e-14. Array with a range smaller 1.e-12 are set
    ! to zero in grib output. To avoid this problem the co2trans was multiplied by 1.e6 (mo_couple).
    
    co2_flux_ocean(1:kproma,jrow) = amco2 * co2trans(1:kproma,jrow) * wind10w(1:kproma,jrow)**2 * &
                                    (co2ocean(1:kproma,jrow) - (co2m1(1:kproma,klev,jrow) * &
                                    1.0E+06_dp * amd / amco2)) * 1.e-6

  END SUBROUTINE co2_flux_atmosphere_ocean
#endif

  SUBROUTINE diag_co2(jrow, kproma, kbdim, klevp1, paphm1, pfrl, pfrw, pfri)

    USE mo_time_control, ONLY: delta_time
    USE mo_mpi,          ONLY: p_pe, p_io

    ! Scalar arguments
    INTEGER, INTENT(in) :: jrow, kproma, kbdim, klevp1

    ! Array arguments
    REAL(dp), INTENT(in) :: paphm1(kbdim,klevp1)
    REAL(dp), INTENT(in) :: pfrl(kbdim)          ! Land fraction
    REAL(dp), INTENT(in) :: pfrw(kbdim)          ! Ocean fraction
    REAL(dp), INTENT(in) :: pfri(kbdim)          ! Sea ice fraction

    INTEGER :: jk

    co2_flux_acc(1:kproma,jrow)       = co2_flux_acc(1:kproma,jrow)                                  &
                                        + co2_flux(1:kproma,jrow) * delta_time
    co2_flux_land_acc(1:kproma,jrow)  = co2_flux_land_acc(1:kproma,jrow)                             &
                                        + pfrl(1:kproma) * co2_flux_land(1:kproma,jrow) * delta_time
    co2_flux_ocean_acc(1:kproma,jrow) = co2_flux_ocean_acc(1:kproma,jrow)                            &
                                        + (pfrw(1:kproma)+pfri(1:kproma))                            &
                                          * co2_flux_ocean(1:kproma,jrow) * delta_time
    co2_flux_anthro_acc(1:kproma,jrow) = co2_flux_anthro_acc(1:kproma,jrow)                          &
                                         + (co2_burden(1:kproma,jrow)-co2_burden_old(1:kproma,jrow)) &
                                         - co2_flux(1:kproma,jrow) * delta_time
    IF (ico2_emis > 0) THEN
       co2_emission_acc(1:kproma,jrow)   = co2_emission_acc(1:kproma,jrow)                              &
                                           + co2_emission(1:kproma,jrow) * delta_time
    END IF
    
    co2_burden_acc(1:kproma,jrow) = co2_burden_acc(1:kproma,jrow)                                   &
                                    + co2_burden(1:kproma,jrow) * delta_time


  END SUBROUTINE diag_co2

  SUBROUTINE co2_flux_correction

    USE mo_control,       ONLY: ngl, nlon
    USE mo_decomposition, ONLY: gl_dc => global_decomposition
    USE mo_transpose,     ONLY: gather_gp, scatter_gp
    USE mo_gaussgrid,     ONLY: gl_budw
    USE mo_mpi,           ONLY: p_pe, p_io, p_bcast
    USE mo_time_control,  ONLY: get_date_components, current_date, previous_date, next_date, delta_time, &
                                lresume, l_putocean, get_time_step

    IMPLICIT NONE

    ! local arrays
    REAL(dp), POINTER :: zco2_burden_corr_glob(:,:)
    REAL(dp)          :: ztime = 0.0_dp
    INTEGER           :: day_of_month, day_of_month_at_prev_ts
    REAL(dp)          :: ztemp(ngl)

    REAL(dp)     :: zco2_burden_corr_mean, zco2_flux_corr_mean !, zco2_flux_corr

    INTEGER :: jlat, istep

    IF (.NOT. lco2) RETURN

    ! Set values to zero if CO2 flux correction not requested (e.g. in case of CO2 nudging)
    IF (.NOT. lco2_flxcor) THEN
       co2_flux_corr_acc1(:,:) = 0._dp
       co2_flux_corr_acc1_old(:,:) = 0._dp
       RETURN
    END IF

    CALL get_date_components(next_date,DAY=day_of_month)
    CALL get_date_components(current_date, DAY=day_of_month_at_prev_ts)
    l_co2flxcorr = day_of_month /= day_of_month_at_prev_ts

    ztime = 86400._dp
    istep = get_time_step()

!!$    CALL message('co2_flux_correction','istep = '//TRIM(int2string(istep))//', ztime = '//&
!!$                                        TRIM(real2string(ztime))//')')

    IF (.NOT. l_co2flxcorr) THEN
!!$       CALL message('co2_flux_correction','Not doing flux correction')
       RETURN
    END IF

    CALL message('co2_flux_correction','Computing global mean CO2 flux correction')

    IF (p_pe == p_io) THEN
       ALLOCATE(zco2_burden_corr_glob(nlon,ngl))
    END IF

!!$    co2_flux_corr_old = co2_flux_corr
    co2_flux_corr_acc1_old(:,:) = co2_flux_corr_acc1(:,:)

    CALL gather_gp (zco2_burden_corr_glob, co2_burden_corr_acc1, gl_dc)

    ! Global mean
    IF (p_pe == p_io) THEN
       DO jlat=1,ngl
          ztemp(jlat) = SUM(zco2_burden_corr_glob(1:nlon,jlat)) * gl_budw(jlat)
       END DO
       zco2_burden_corr_mean = SUM(ztemp)

       ! Apply global CO2 burden correction to CO2 flux
       zco2_flux_corr_mean = zco2_burden_corr_mean / ztime
       WRITE(message_text,*) 'CO2 flux correction (global): ',zco2_flux_corr_mean, ' kg m-2 s-1'
       CALL message('', message_text)

    END IF

    CALL p_bcast(zco2_flux_corr_mean, p_io)
    co2_flux_corr_acc1(:,:) = zco2_flux_corr_mean

    IF (p_pe == p_io) THEN
       DEALLOCATE(zco2_burden_corr_glob)
    END IF

    RETURN

  END SUBROUTINE co2_flux_correction


END MODULE mo_co2
