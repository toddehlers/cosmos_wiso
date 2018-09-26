!+ Definitions for soil parameters
! 
Module mo_soil

    !
    ! Description:
    !   Variables that hold soil parameters for JSBACH
    !
    ! Current Code Owner: jsbach_admin
    !
    ! History:
    !
    ! Version   Date        Comment
    ! -------   ----        -------
    ! 0.1       2001/07/01  Original code. Reiner Schnur
    !
    ! Code Description:
    !   Language:           Fortran 90.
    !   Software Standards: "European Standards for Writing and
    !     Documenting Exchangeable Fortran 90 Code".
    !
    ! Modules used:
    !
    USE mo_netCDF, ONLY: FILE_INFO, NF_MAX_NAME, BELOWSUR, SOILLEV
    USE mo_jsbach, ONLY: debug, test_stream
    USE mo_jsbach_grid, ONLY: kstart, kend, nidx
    USE mo_linked_list, ONLY: t_stream
    USE mo_mpi, ONLY: p_parallel, p_parallel_io, p_bcast, p_io, p_pe, p_nprocs
    USE mo_doctor, ONLY: nout, nerr
    USE mo_kind, ONLY: dp
    USE mo_filename, ONLY: path_limit
    USE mo_exception, ONLY: finish, message, message_text
    USE mo_test, ONLY: test

!---wiso-code
    USE mo_wiso, ONLY: lwiso, nwiso
!---wiso-code-end

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: soil_type, soil_param_type, init_soil, update_soil, soil_diagnostics, get_soil_diag , get_soil_diag_wiso

    ! Global (i.e. public) Declarations:
    ! Global Type Definitions:
    !! Structure to hold global soil parameters for each land cell and soil layer.
    TYPE soil_param_type
        !
        ! Variables describing soil characteristics and being read from initialization files
        !
        REAL(dp), POINTER, DIMENSION(:) :: & ! (nland)
        Ds, & !! Fraction of maximum subsurface flow rate
        Dsmax, & !! Maximum subsurface flow rate [mm/day]
        Ws, & !! Fraction of maximum soil moisture where non-linear baseflow occurs
        B_infilt, & !! Infiltration parameter
        AvgTemp, & !! Average soil temperature used as the bottom
        !!    boundary for soil heat flux solutions
        Roughness, & !! Surface roughness of bare soil
        VolHeatCapacity, & !! Volumetric heat capacity of soil [J/m**3/K]
        ThermalDiffusivity !! Thermal diffusivity of soil [m**2/s]

        ! The dimension of the following pointer is: (nland,nsoil)
        REAL(dp), POINTER :: Ksat(:,:) !! Saturated hydraulic conductivity [mm/day]
        REAL(dp), POINTER :: Wpwp(:,:) !! Soil moisture content at permanent wilting point [m]
        REAL(dp), POINTER :: Depth(:,:) !! Depth of soil layer [m]
        REAL(dp), POINTER :: InitMoisture(:,:) !! Initial layer moisture level [m]
        REAL(dp), POINTER :: MaxMoisture(:,:) !! Maximum moisture content per layer (field capacity) [m]
        REAL(dp), POINTER :: BulkDensity(:,:) !! Bulk density of soil layer [kg/m^3]
        REAL(dp), POINTER :: SoilDensity(:,:) !! Soil particle density, normally 2685 kg/m^3
        REAL(dp), POINTER :: Porosity(:,:) !! Soil Porosity [fraction]
        REAL(dp), POINTER :: MaxInfiltration(:,:) !! Maximum infiltration rate
    END TYPE soil_param_type
    !
    !----------------------------------------------------------------------------------------------------------
    ! Soil state variables
    TYPE soil_type
        INTEGER :: nsoil !! Number of soil layers
        INTEGER :: ntsoil !! Number of soil layers for soil temperature calculation
        INTEGER :: ntiles !! Number of soil tiles
        REAL(dp), POINTER, DIMENSION(:) :: & !! (nland)
        csat, & !!
        cair, & !!
        relative_humidity_air, & !! Relative humidity of the air in the lowest layer of the atmosphere
        csat_transpiration

!---wiso-code
        REAL(dp), POINTER, DIMENSION(:,:) :: & !! (nland,nwiso)
        wiso_csat,               & !!
        wiso_cair,               & !!
        wiso_csat_fra,           & !!
        wiso_cair_fra,           & 
        wiso_csat_transpiration, &
        wiso_csat_transp_fra,    &
        q_wiso_Acoef_new,       &
        q_wiso_Bcoef_new
!---wiso-code-end

        REAL(dp), POINTER, DIMENSION(:,:) :: & !! (nland, ntiles)
        skin_reservoir, & !! Water content of skin reservoir
        albedo, & !! background albedo (ECHAM albedo scheme) or = albedo_soil_vis (JSBACH albedo scheme)
        albedo_vegetation_vis, & !! albedo of vegetation in the visible range
        albedo_vegetation_nir, & !! albedo of vegetation in the NIR range
        albedo_soil_vis, & !! albedo of soil in the visible range
        albedo_soil_nir, & !! albedo of soil in the NIR range
        snow, & !! Snow on ground [m water equivalent]
!---wiso-code
        snowglac, & !! Snow on ground on glaciers [m water equivalent]
!---wiso-code-end
        snow_fract, & !! Fraction of snow covered ground
        snow_acc, & !! Snow accumulation (budget) at non-glacier points [kg/(m^2 s)]
        snow_melt_acc, & !! Snow melt (accumulated) [kg/(m^2 s)]
        runoff_acc, & !! Total runoff (surface runoff + drainage) at non-glacier points (accumulated)
        drainage_acc, & !! Drainage at non-glacier points (accumulated) [kg/(m^2 s)]
        glacier_depth, & !! Glacier depth (including snow) [m water equivalent]
        glacier_precip_minus_evap_acc, & !! Precipitation minus evaporation for glaciers (accumulated)
        glacier_runoff_acc, & !! Glacier runoff (rain+snow/ice melt) (accumulated) [kg/(m^2 s)]
        surface_temperature, & !! Temperature of the surface
        surface_temperature_old, & !! Temperature of the surface at time step t-dt
        surface_temperature_unfiltered, & !! Temperature of the surface at time step t+dt (unfiltered)
        radiative_temperature, & !! Temp for radiation derived from dry_satic_energy_new
        sat_surface_specific_humidity, & !! Saturated surface specific humidity
        evapotranspiration, & !! Evaporation and transpiration (for time step)
        evapotranspiration_acc, & !! Evaporation and transpiration (accumulated)
        evaporation_pot, & !! Potential evaporation (for time step)
        transpiration, & !! Transpiration (for time step)
        transpiration_acc, & !! Transpiration (accumulated)
        sensible_heat_flux, & !! Sensible heat flux (for time step) [W/m^2]
        sensible_heat_acc, & !!                  (accumulated)
        latent_heat_flux, & !! Latent heat flux (for time step) [W/m^2]
        latent_heat_acc, & !!                  (accumulated)
        ground_heat_flux, & !! Ground heat flux (for time step) [W/m^2]
        ground_heat_flux_acc, & !!                  (accumulated)
        heat_capacity, & !! Surface heat capacity [J/m**2/K]
        dry_static_energy_new

!---wiso-code
        REAL(dp), POINTER, DIMENSION(:,:) :: & !! (nland, ntiles)
        remember_cveg, &
        remember_canopy, &
        remember_rel_hum, &
        remember_echamzchl, &
        remember_cfrac, &
        remember_snfrac, &
        remember_wlfrac, &
        remember_wind, &
        remember_airm, &
        remember_satsphum
        REAL(dp), POINTER, DIMENSION(:,:,:) :: & !! (nland,nwiso,ntiles)
        wiso_skin_reservoir, & !! Water content of skin reservoir - Isotopes
        wiso_snow, & !! Snow on ground [m water equivalent] - Isotopes
        wiso_snowglac, & !! Snow on ground on glaciers [m water equivalent] - Isotopes
        wiso_snow_acc, & !! Snow accumulation (budget) at non- points [kg/(m^2 s)] - Isotopes
        wiso_snow_melt_acc, & !! Snow melt (accumulated) [kg/(m^2 s)] - Isotopes
        wiso_runoff_acc, & !! Total runoff (surface runoff + drainage) at non-glacier points (accumulated) - Isotopes
        wiso_drainage_acc, & !! Drainage at non-glacier points (accumulated) [kg/(m^2 s)] - Isotopes
        wiso_glacier_depth, & !! Glacier depth (including snow) [m water equivalent] - Isotopes
        wiso_glac_p_minus_e_acc, & !! - Precipitation minus evaporation for glaciers (accumulated) - Isotopes
        wiso_glacier_runoff_acc, & !! Glacier runoff (rain+snow/ice melt) (accumulated) [kg/(m^2 s)] - Isotopes
        wiso_sat_surface_specific_hum, & !! Saturated surface specific humidity - Isotopes
        wiso_evapotranspiration, & !! Evaporation and transpiration (for time step) - Isotopes
        wiso_evapotranspiration_acc, & !! Evaporation and transpiration (accumulated) - Isotopes
        wiso_evaporation_pot, & !! Potential evaporation (for time step) - Isotopes
        wiso_transpiration, & !! Transpiration (for time step) - Isotopes
        wiso_transpiration_acc !! - Transpiration (accumulated) - Isotopes
!---wiso-code-end

        REAL(dp), POINTER, DIMENSION(:,:,:) :: & !! (nland, nsoil, ntiles)
        moisture, & !! Moisture content of the unfrozen sub-layer
        relative_moisture, & !! relative water content of the soil layers to max capacity
        root_fract !! Fraction of roots within each soil layer

!---wiso-code
    !    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: & !! (nland, nsoils, nwiso, ntile)
    !   Vorerst nur drei Dimensionen, da sonst alle Speicher umgeschrieben werden muessen
        REAL(dp), POINTER, DIMENSION(:,:,:) :: & !! (nland, nwiso, ntiles)
        wiso_moisture !! Moisture content of the unfrozen sub-layer - Isotopes
!---wiso-code-end

        REAL(dp), POINTER, DIMENSION(:,:,:) :: & !! (nland, ntsoil, ntiles)
        soil_temperature, & !! Soil temperature
        c_soil_temperature, & !! Coefficient to compute soil_temperature
        d_soil_temperature !!                 "
    END TYPE soil_type
    !
    !----------------------------------------------------------------------------------------------------------
    TYPE soil_diag_type
        REAL(dp), POINTER, DIMENSION(:) :: &
        surface_temperature, &
        surface_radiative_temp, &
        sat_surface_specific_humidity, &
        ground_heat_flux_acc, &
        heat_capacity, &
        evapotranspiration_acc, &
        transpiration_acc, &
        sensible_heat_acc, &
        latent_heat_acc, &
        albedo, &
        skin_reservoir, &
        snow, &
!---wiso-code
        snowglac, &
!---wiso-code-end
        snow_acc, &
        snow_fract, &
        snow_melt_acc, &
        runoff_acc, &
        drainage_acc, &
        glacier_depth, &
        glacier_precip_minus_evap_acc, &
        glacier_runoff_acc

        REAL(dp), POINTER, DIMENSION(:,:) :: &
        moisture, &
        soil_temperature

!---wiso-code
        REAL(dp), POINTER, DIMENSION(:,:) :: &
        wiso_sat_surface_specific_hum, &
        wiso_evapotranspiration_acc, &
        wiso_transpiration_acc, &
        wiso_skin_reservoir, &
        wiso_snow, &
        wiso_snowglac, &
        wiso_snow_acc, &
        wiso_snow_melt_acc, &
        wiso_runoff_acc, &
        wiso_drainage_acc, &
        wiso_glacier_depth, &
        wiso_glac_p_minus_e_acc, &
        wiso_glacier_runoff_acc, &
        wiso_moisture
!---wiso-code-end

    END TYPE soil_diag_type
    !
    !----------------------------------------------------------------------------------------------------------
    TYPE soil_options_type
        REAL(dp) :: SkinReservoirMax !! Maximum content of skin reservoir over bare soil (set to 0. to disable)
        REAL(dp) :: MoistureFractCritical !! Fractional soil moisture content at the critical point (fraction of maximum moisture)
        REAL(dp) :: MoistureFractWilting !! Fractional soil moisture content at the wilting point (fraction of maximum moisture)
        REAL(dp) :: MaxMoistureLimit !! Upper limit for maximum soil moisture content
        REAL(dp) :: CriticalSnowDepth !! Critical snow depth for correction of new surface temperature [m water equivalent]
    END TYPE soil_options_type
    !
    !----------------------------------------------------------------------------------------------------------
    INTEGER, SAVE :: nsoil = -1 !! Number of soil layers (for water balance)
    INTEGER, SAVE :: ntsoil = -1 !! Number of soil layers (for energy balance)
    LOGICAL, SAVE :: flag_soil_moisture = .true. !! Flag eventually indicating, that negatives values of soil moisture
    !! already occurred

    !  TYPE(soil_param_type), SAVE :: soil_param !! Global soil parameters

    TYPE(t_stream), POINTER :: IO_soil !! Memory stream for soil model state
    TYPE(t_stream), POINTER :: IO_diag !! Memory stream for soil diagnostic output
!---wiso-code
    TYPE(t_stream), POINTER :: IO_soil_wiso !! Memory stream for soil model state - water isotopes
    TYPE(t_stream), POINTER :: IO_diag_wiso !! Memory stream for soil diagnostic output - water isotopes
!---wiso-code-end
    TYPE(FILE_INFO), SAVE :: soil_file !! Input file for soil parameters

    TYPE(soil_diag_type), SAVE :: soil_diag
    TYPE(soil_options_type), SAVE :: soil_options

    ! Parameters used for snow cover fraction
    REAL(dp), PARAMETER :: &
    zepsec = 1.E-12_dp, &
    zsigfac = 0.15_dp, &
    zqsncr = 0.95_dp ! inverse of equivalent water height when snow is considered to completely cover the ground

CONTAINS
    !
    !=================================================================================================
    !
    SUBROUTINE config_soil

        USE mo_jsbach, ONLY: nml_unit
        USE mo_namelist, ONLY: position_nml, POSITIONED

        !! locals
        INTEGER :: read_status

        !! Namelist Parameters
        REAL(dp) :: skin_res_max !! Maximum content of skin reservoir of bare soil [m]
        REAL(dp) :: moist_crit_fract !! Critical value of soil moisutre
        REAL(dp) :: moist_wilt_fract !! Soil moisture content at permanent wilting point
        REAL(dp) :: moist_max_limit !! Upper limit for maximum soil moisture content
        REAL(dp) :: crit_snow_depth !! Critical snow depth for correction of surface temperature for melt

        INCLUDE 'soil_ctl.inc'

        !! Read namelist soil_ctl

        IF (p_parallel_io) THEN

            ! define default values
            skin_res_max = 2.E-4_dp
            moist_crit_fract = 0.75_dp
            moist_wilt_fract = 0.35_dp
            moist_max_limit = -1.0_dp
            crit_snow_depth = 5.85036E-3_dp

            CALL position_nml('SOIL_CTL', status = read_status)
            SELECT CASE (read_status)
            CASE (POSITIONED)
                READ (nml_unit, soil_ctl)
                CALL message('config_soil', 'Namelist SOIL_CTL: ')
                WRITE(nout, soil_ctl)
            END SELECT
        ENDIF

        IF (p_parallel_io) THEN
            soil_options % SkinReservoirMax = skin_res_max
            WRITE(message_text, *) 'Maximum content of (bare soil) skin reservoir: ', soil_options % SkinReservoirMax
            CALL message('config_soil', message_text)

            soil_options % MoistureFractCritical = moist_crit_fract
            WRITE(message_text, *) 'Fractional soil moisture at critical point: ', soil_options % MoistureFractCritical
            CALL message('config_soil', message_text)

            soil_options % MoistureFractWilting = moist_wilt_fract
            WRITE(message_text, *) 'Fractional soil moisture at permanent wilting point: ', soil_options % MoistureFractWilting
            CALL message('config_soil', message_text)

            soil_options % MaxMoistureLimit = moist_max_limit
            WRITE(message_text, *) 'Upper limit for maximum soil moisture content: ', soil_options % MaxMoistureLimit
            CALL message('config_soil', message_text)

            soil_options % CriticalSnowDepth = crit_snow_depth
            WRITE(message_text, *) 'Critical snow depth for correction of surface temperature for melt: ', &
                                    soil_options % CriticalSnowDepth
            CALL message('config_soil', message_text)
        END IF

        IF (p_parallel) THEN
            CALL p_bcast(soil_options % SkinReservoirMax, p_io)
            CALL p_bcast(soil_options % MoistureFractCritical, p_io)
            CALL p_bcast(soil_options % MoistureFractWilting, p_io)
            CALL p_bcast(soil_options % MaxMoistureLimit, p_io)
            CALL p_bcast(soil_options % CriticalSnowDepth, p_io)
        END IF

    END SUBROUTINE config_soil
    !
    !=================================================================================================
    !
    SUBROUTINE soil_init_io(IO_file_name)

        USE mo_netCDF, ONLY: add_dim, IO_inq_dimid, IO_inq_dimlen
        USE mo_io, ONLY: IO_open, IO_READ, IO_close

        CHARACTER(NF_MAX_NAME), INTENT(in) :: IO_file_name
        INTEGER :: IO_file_id, IO_dim_id

        IF (debug) CALL message('soil_init_io', '')

        IF (p_parallel_io) THEN

            ! Get number of soil layers
            !       WRITE(nout, *) 'Reading number of soil layers from', TRIM(IO_file_name)
            soil_file % opened = .FALSE.
            CALL IO_open(TRIM(IO_file_name), soil_file, IO_READ)
            IO_file_id = soil_file % file_id
            !       CALL IO_inq_dimid(IO_file_id, 'soil_layer', IO_dim_id)
            !       CALL IO_inq_dimlen(IO_file_id, IO_dim_id, nsoil)
            CALL IO_close(soil_file)
            nsoil = 1 !! FOR TESTING ONLY
            ntsoil = 5 !! FOR TESTING ONLY

            WRITE(nout, *) 'soil_init_io - Number of soil layers (water) :', nsoil
            WRITE(nout, *) 'soil_init_io - Number of soil layers (energy):', ntsoil

        ENDIF

        IF (p_parallel) THEN
            CALL p_bcast(nsoil, p_io)
            CALL p_bcast(ntsoil, p_io)
        ENDIF

        ! Add dimensions
        IF (debug) CALL message('soil_init_io', 'Adding dimensions')

        CALL add_dim("soil_layer", nsoil, longname = 'soil levels (water)', levtyp = 71, &
        units = 'cm', value = (/1._dp/), indx = SOILLEV)
        CALL add_dim("soil_layer_temp", ntsoil, longname = 'soil levels (energy)', levtyp = 111, &
        units = 'cm', value = (/3._dp, 19._dp, 78._dp, 268._dp, 698._dp/), indx = BELOWSUR)

    END SUBROUTINE soil_init_io
    !
    !=================================================================================================

    SUBROUTINE soil_init_memory(g_nland, l_nland, ntiles, useDynveg, soil, fileformat, diag_stream, diag_wiso_stream, &
                                stream, wiso_stream)

        USE mo_jsbach, ONLY: missing_value, lpost_echam
        USE mo_linked_list, ONLY: NETCDF, LAND, TILES
        USE mo_memory_base, ONLY: new_stream, default_stream_setting, &
        add => add_stream_element
        USE mo_netCDF, ONLY: max_dim_name
        USE mo_grib, ONLY: land_table

!---wiso-code
        USE mo_wiso, ONLY: lwiso, nwiso
!---wiso-code-end

        INTEGER, INTENT(in) :: g_nland, l_nland, ntiles
        INTEGER, INTENT(in) :: fileformat
        LOGICAL, INTENT(in) :: useDynveg

        TYPE(soil_type), INTENT(inout) :: soil
        TYPE(t_stream), POINTER :: diag_stream
        TYPE(t_stream), POINTER, OPTIONAL :: stream
!---wiso-code
        TYPE(t_stream), POINTER, OPTIONAL :: diag_wiso_stream
        TYPE(t_stream), POINTER, OPTIONAL :: wiso_stream
!---wiso-code-end

        INTEGER :: dim1p(1), dim1(1)
        INTEGER :: dim2p(2), dim2(2)
        INTEGER :: dim3p(2), dim3(2)
        INTEGER :: dim4p(3), dim4(3)
        INTEGER :: dim5p(3), dim5(3)
        INTEGER :: dim6p(2), dim6(2)
        CHARACTER(LEN = max_dim_name) :: dim1n(1), dim2n(2), dim3n(2), dim4n(3), dim5n(3), dim6n(2)
!---wiso-code
	INTEGER :: iwiso
        INTEGER :: dim7p(2), dim7(2)
        INTEGER :: dim8p(3), dim8(3)
        INTEGER :: dim9p(4), dim9(4)
        CHARACTER(LEN = max_dim_name) :: dim7n(2), dim8n(3), dim9n(4)
!---wiso-code-end

        IF (ntiles < 1) CALL finish('soil_init_memory', 'tiles not initialized')
        IF (nsoil < 1) CALL finish('soil_init_memory', 'soil layers not initialized (water balance)')
        IF (ntsoil < 1) CALL finish('soil_init_memory', 'soil layers not initialized (energy balance)')

        IF (ASSOCIATED(diag_stream)) THEN
            IO_diag => diag_stream
        ELSE
            CALL finish('soil_init_memory', 'Diagnostic stream not present')
        END IF

        IF (PRESENT(stream)) THEN
            IF (.NOT.ASSOCIATED(stream)) THEN
                ! Add new stream
                CALL new_stream(stream, 'soil', filetype = fileformat)
                ! Set default stream options
                CALL default_stream_setting(stream, table = land_table, repr = LAND, lpost = .TRUE., lrerun = .TRUE.)
            ENDIF
            IO_soil => stream
        ELSE
            ! Add new stream
            CALL new_stream(IO_soil, 'soil', filetype = fileformat)
            ! Set default stream options
            CALL default_stream_setting(IO_soil, table = land_table, repr = LAND, lpost = .TRUE., lrerun = .TRUE., &
                                        leveltype = TILES)
        ENDIF

!---wiso-code
        IF (lwiso) THEN

        IF (ASSOCIATED(diag_wiso_stream)) THEN
            IO_diag_wiso => diag_wiso_stream
        ELSE
            CALL finish('soil_init_memory', 'Isotope diagnostic stream not present')
        END IF

        IF (PRESENT(wiso_stream)) THEN
            IF (.NOT.ASSOCIATED(wiso_stream)) THEN
                ! Add new stream
                CALL new_stream(wiso_stream, 'soil_wiso', filetype = fileformat)
                ! Set default stream options
                CALL default_stream_setting(wiso_stream, table = land_table, repr = LAND, lpost = .TRUE., lrerun = .TRUE.)
            ENDIF
            IO_soil_wiso => wiso_stream
        ELSE
            ! Add new stream
            CALL new_stream(IO_soil_wiso, 'soil_wiso', filetype = fileformat)
            ! Set default stream options
            CALL default_stream_setting(IO_soil_wiso, table = land_table, repr = LAND, lpost = .TRUE., lrerun = .TRUE., &
                                        leveltype = TILES)
        ENDIF
        
        ENDIF
!---wiso-code-end

        dim1p = (/ l_nland /)
        dim1 = (/ g_nland /)
        dim1n = (/ 'landpoint' /)

        dim2p = (/ l_nland, nsoil /)
        dim2 = (/ g_nland, nsoil /)
        dim2n(1) = 'landpoint'
        dim2n(2) = 'soil_layer'

        dim3p = (/ l_nland, ntiles /)
        dim3 = (/ g_nland, ntiles /)
        dim3n(1) = 'landpoint'
        dim3n(2) = 'tiles'

        dim4p = (/ l_nland, nsoil, ntiles /)
        dim4 = (/ g_nland, nsoil, ntiles /)
        dim4n(1) = 'landpoint'
        dim4n(2) = 'soil_layer'
        dim4n(3) = 'tiles'

        dim5p = (/ l_nland, ntsoil, ntiles /)
        dim5 = (/ g_nland, ntsoil, ntiles /)
        dim5n(1) = 'landpoint'
        dim5n(2) = 'soil_layer_temp'
        dim5n(3) = 'tiles'

        dim6p = (/ l_nland, ntsoil /)
        dim6 = (/ g_nland, ntsoil /)
        dim6n(1) = 'landpoint'
        dim6n(2) = 'soil_layer_temp'

!---wiso-code
        dim7p = (/ l_nland, nwiso /)
        dim7  = (/ g_nland, nwiso /)
        dim7n(1) = 'landpoint'
        dim7n(2) = 'isotopetyp'

        dim8p = (/ l_nland, nwiso, ntiles  /)
        dim8  = (/ g_nland, nwiso, ntiles  /)
        dim8n(1) = 'landpoint'
        dim8n(2) = 'isotopetyp'
        dim8n(3) = 'tiles'

! Attention - nsoil is not coded yet, because nsoil=1 in this echam&jsbach version
!    dim9p = (/ l_nland, nsoil, ntiles, nwiso  /)
!    dim9  = (/ g_nland, nsoil, ntiles, nwiso  /)
!    dim9n(1) = 'landpoint'
!    dim9n(2) = 'soil_layer'
!    dim9n(3) = 'tiles'
!    dim9n(4) = 'isotopetyp'
!---wiso-code-end

        CALL add(IO_soil, 'surface_temperature', soil % surface_temperature, longname = 'Land Surface Temperature', &
        units = 'K', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 35, lpost = .FALSE.)
        CALL add(IO_soil, 'surface_temperature_old', soil % surface_temperature_old, &
        ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 36, lpost = .FALSE.)
        CALL add(IO_soil, 'surface_temp_unfiltered', soil % surface_temperature_unfiltered, &
        ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 37, lpost = .FALSE.)
        CALL add(IO_soil, 'surface_radiative_temp', soil % radiative_temperature, longname = 'Surface Radiative Temperature', &
        ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 38, lpost = .FALSE.)
        CALL add(IO_soil, 'surface_qsat', soil % sat_surface_specific_humidity, &
        longname = 'Soil Specific Humidity at Saturation', &
        units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 39, lpost = .FALSE.)
        CALL add(IO_soil, 'ground_heat_flux_inst', soil % ground_heat_flux, longname = 'Ground Heat Flux (inst.)', &
        units = 'W m-2', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 40, lpost = .FALSE.)
        CALL add(IO_soil, 'ground_heat_flux', soil % ground_heat_flux_acc, longname = 'Ground Heat Flux (ave.)', &
        units = 'W m-2', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 41, laccu = .TRUE., lpost = .FALSE.)
        CALL add(IO_soil, 'heat_capacity', soil % heat_capacity, longname = 'Heat Capacity of the Uppermost Soil Layer', &
        units = 'J m-2 K-1', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 42, lpost = .FALSE.)
        CALL add(IO_soil, 'evapotranspiration_inst', soil % evapotranspiration, longname = 'Evapotranspiration (inst.)', &
        units = 'kg m-2 s-1', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 43, lpost = .FALSE.)
        CALL add(IO_soil, 'evapotranspiration', soil % evapotranspiration_acc, longname = 'Evapotranspiration (avg)', &
        units = 'kg m-2 s-1', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 44, laccu = .TRUE., lpost = .FALSE.)
        CALL add(IO_soil, 'evaporation_pot_inst', soil % evaporation_pot, longname = 'Potential Evaporation', &
        units = 'kg m-2 s-1', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 45, lpost = .FALSE.)
        CALL add(IO_soil, 'transpiration_inst', soil % transpiration, longname = 'Transpiration (inst.)', &
        units = 'kg m-2 s-1', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 75, lpost = .FALSE.)
        CALL add(IO_soil, 'transpiration', soil % transpiration_acc, longname = 'Transpiration (avg)', &
        units = 'kg m-2 s-1', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 76, laccu = .TRUE., lpost = .FALSE.)
        CALL add(IO_soil, 'sensible_heat_flux_inst', soil % sensible_heat_flux, longname = 'Sensible Heat Flux (inst.)', &
        units = 'W m-2', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 46, lpost = .FALSE.)
        CALL add(IO_soil, 'sensible_heat_flux', soil % sensible_heat_acc, longname = 'Sensible Heat Flux (avg)', &
        units = 'W m-2', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 47, laccu = .TRUE., lpost = .FALSE.)
        CALL add(IO_soil, 'latent_heat_flux_inst', soil % latent_heat_flux, longname = 'Latent Heat Flux (inst.)', &
        units = 'W m-2', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 48, lpost = .FALSE.)
        CALL add(IO_soil, 'latent_heat_flux', soil % latent_heat_acc, longname = 'Latent Heat Flux (avg)', &
        units = 'W m-2', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 49, laccu = .TRUE., lpost = .FALSE.)
        CALL add(IO_soil, 'soil_albedo', soil % albedo, longname = 'Soil Albedo', &
        units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 50, lpost = .FALSE., lrerun = .FALSE.)
        CALL add(IO_soil, 'albedo_veg_vis', soil % albedo_vegetation_vis, &
        ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 51, lpost = .FALSE., lrerun = .FALSE.)
        CALL add(IO_soil, 'albedo_veg_nir', soil % albedo_vegetation_nir, &
        ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 52, lpost = .FALSE., lrerun = .FALSE.)
        CALL add(IO_soil, 'albedo_soil_vis', soil % albedo_soil_vis, &
        ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 53, lpost = .FALSE., lrerun = .FALSE.)
        CALL add(IO_soil, 'albedo_soil_nir', soil % albedo_soil_nir, &
        ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 54, lpost = .FALSE., lrerun = .FALSE.)
        CALL add(IO_soil, 'skin_reservoir', soil % skin_reservoir, longname = 'Water Content of the Skin Reservoir', &
        units = 'm', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 55, lpost = .TRUE., &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_soil, 'soil_moisture', soil % moisture, longname = 'Water Content of the Soil Layers', &
        units = 'm', ldims = dim4p, gdims = dim4, dimnames = dim4n, code = 57, lpost = .FALSE., leveltype = SOILLEV)
        CALL add(IO_soil, 'rel_soil_moisture', soil % relative_moisture, longname = 'Relative Water Content of the Soil Layers', &
        units = '', ldims = dim4p, gdims = dim4, dimnames = dim4n, code = 58, lpost = .FALSE., leveltype = SOILLEV)
        CALL add(IO_soil, 'snow', soil % snow, longname = 'Snow Depth in Water Equivalent', &
        units = 'm', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 59, lpost = .FALSE.)
        CALL add(IO_soil, 'snow_fract', soil % snow_fract, longname = 'Snow Cover Fraction', &
        units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 60, lpost = .FALSE.)
        CALL add(IO_soil, 'snow_acc', soil % snow_acc, longname = 'Snowfall (avg)', &
        units = 'kg m-2 s-1', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 61, laccu = .TRUE., lpost = .FALSE.)
        CALL add(IO_soil, 'snow_melt', soil % snow_melt_acc, longname = 'Snow Melt Flux', &
        units = 'kg m-2 s-1', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 62, laccu = .TRUE., lpost = .FALSE.)
        CALL add(IO_soil, 'runoff', soil % runoff_acc, longname = 'Surface Runoff and Drainage', &
        units = 'kg m-2 s-1', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 63, laccu = .TRUE., lpost = .FALSE.)
        CALL add(IO_soil, 'drainage', soil % drainage_acc, longname = 'Drainage', &
        units = 'kg m-2 s-1', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 64, laccu = .TRUE., lpost = .FALSE.)
        CALL add(IO_soil, 'glacier_depth', soil % glacier_depth, longname = 'Glacier Depth', &
        units = 'm', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 65, lpost = .FALSE.)
        CALL add(IO_soil, 'glacier_precip_minus_evap', soil % glacier_precip_minus_evap_acc, &
        longname = 'Precipitation minus Evaporation over Glaciers', &
        units = 'm', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 66, laccu = .TRUE., lpost = .FALSE.)
        CALL add(IO_soil, 'glacier_runoff', soil % glacier_runoff_acc, longname = 'Glacier Runoff', &
        units = 'kg m-2 s-1', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 67, laccu = .TRUE., lpost = .FALSE.)
        CALL add(IO_soil, 'soil_temperature', soil % soil_temperature, longname = 'Soil Temperature', &
        units = 'K', ldims = dim5p, gdims = dim5, dimnames = dim5n, code = 68, lpost = .FALSE., leveltype = BELOWSUR)
        CALL add(IO_soil, 'c_soil_temp_coef', soil % c_soil_temperature, &
        ldims = dim5p, gdims = dim5, dimnames = dim5n, code = 69, lpost = .FALSE., leveltype = BELOWSUR)
        CALL add(IO_soil, 'd_soil_temp_coef', soil % d_soil_temperature, &
        ldims = dim5p, gdims = dim5, dimnames = dim5n, code = 70, lpost = .FALSE., leveltype = BELOWSUR)
        CALL add(IO_soil, 'dry_static_energy_new', soil % dry_static_energy_new, &
        ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 71, lpost = .FALSE.)

        CALL add(IO_soil, 'csat', soil % csat, dim1p, dim1, dimnames = dim1n, code = 72, lpost = .FALSE.)
        CALL add(IO_soil, 'cair', soil % cair, dim1p, dim1, dimnames = dim1n, code = 73, lpost = .FALSE.)
        IF (useDynveg) THEN
            CALL add(IO_soil, 'relative_humidity_air_inst', soil % relative_humidity_air, longname = 'Relative Humidity', &
            units = '', ldims = dim1p, gdims = dim1, dimnames = dim1n, code = 56, lpost = .TRUE., &
            lmiss = .TRUE., missval = missing_value)
        ENDIF
        CALL add(IO_soil, 'csat_transpiration', soil % csat_transpiration, dim1p, dim1, dimnames = dim1n, code = 74, &
        lpost = .FALSE.)

!---wiso-code
        IF (lwiso) THEN

           CALL add(IO_soil_wiso, 'snowglac', soil % snowglac, &
                longname = 'Snow Depth on Glaciers in Water Equivalent', &
                units = 'm', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 77, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_surface_qsat', soil % wiso_sat_surface_specific_hum, &
                longname = 'Soil Specific Humidity at Saturation - Water Isotopes',          &
                units = '', ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 80, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_evapotranspiration_inst', soil % wiso_evapotranspiration, &
                longname = 'Evapotranspiration (inst.) - Water Isotopes', &
                units = 'kg m-2 s-1', ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 81, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_evapotranspiration', soil % wiso_evapotranspiration_acc, &
                longname = 'Evapotranspiration (avg) - Water Isotopes', units = 'kg m-2 s-1', &
                ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 82, laccu = .TRUE., lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_evaporation_pot_inst', soil % wiso_evaporation_pot, &
                longname = 'Potential Evaporation - Water Isotopes', units = 'kg m-2 s-1', &
                ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 83, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_transpiration_inst', soil % wiso_transpiration, &
                longname = 'Transpiration (inst.) - Water Isotopes', &
                units = 'kg m-2 s-1', ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 84, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_transpiration', soil % wiso_transpiration_acc, &
                longname = 'Transpiration(avg) - Water Isotopes', units = 'kg m-2 s-1', &
                ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 85, laccu = .TRUE., lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_skin_reservoir', soil % wiso_skin_reservoir, &
                longname = 'Water Content of the Skin Reservoir - Water Isotopes', &
                units = 'm', ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 86, lpost = .TRUE., &
                lmiss = .TRUE., missval = missing_value)
           CALL add(IO_soil_wiso, 'wiso_soil_moisture', soil % wiso_moisture, &
                longname = 'Water Content of the Soil Layers - Water Isotopes', &
                units = 'm', ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 87, lpost = .FALSE., leveltype = SOILLEV)
           CALL add(IO_soil_wiso, 'wiso_snow', soil % wiso_snow, longname = 'Snow Depth in Water Equivalent - Water Isotopes', &
                units = 'm', ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 89, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_snowglac', soil % wiso_snowglac, &
                longname = 'Snow Depth on Glaciers in Water Equivalent - Water Isotopes', &
                units = 'm', ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 78, lpost = .FALSE.) 
           CALL add(IO_soil_wiso, 'wiso_snow_acc', soil % wiso_snow_acc, longname = 'Snowfall (avg) - Water Isotopes', &
                units = 'kg m-2 s-1', ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 91, laccu = .TRUE., lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_snow_melt', soil % wiso_snow_melt_acc, longname = 'Snow Melt Flux - Water Isotopes', &
                units = 'kg m-2 s-1', ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 92, laccu = .TRUE., lpost = .FALSE.) 
           CALL add(IO_soil_wiso, 'wiso_runoff', soil % wiso_runoff_acc, longname = 'Surface Runoff and Drainage - Water Isotopes',&
                units = 'kg m-2 s-1', ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 93, laccu = .TRUE., lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_drainage', soil % wiso_drainage_acc, longname = 'Drainage - Water Isotopes', &
                units = 'kg m-2 s-1', ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 94, laccu = .TRUE., lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_glacier_depth', soil % wiso_glacier_depth, longname = 'Glacier Depth - Water Isotopes', &
                units = 'm', ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 95, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_glac_p_minus_e', soil % wiso_glac_p_minus_e_acc, &
               longname = 'Precipitation minus Evaporation over Glaciers - Water Isotopes', &
               units = 'm', ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 96, laccu = .TRUE., lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_glacier_runoff', soil % wiso_glacier_runoff_acc, &
                longname = 'Glacier Runoff - Water Isotopes', &
                units = 'kg m-2 s-1', ldims = dim8p, gdims = dim8, dimnames = dim8n, code = 97, laccu = .TRUE., lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_csat', soil % wiso_csat, ldims = dim7p, gdims = dim7, dimnames = dim7n, &
                code = 98, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_cair', soil % wiso_cair, ldims = dim7p, gdims = dim7, dimnames = dim7n, &
                code = 99, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_csat_transpiration', soil % wiso_csat_transpiration, ldims = dim7p, gdims = dim7, &
                dimnames = dim7n, code = 100, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_csat_fra', soil % wiso_csat_fra, ldims = dim7p, gdims = dim7, dimnames = dim7n, &
                code = 122, lpost = .FALSE., lrerun = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_cair_fra', soil % wiso_cair_fra, ldims = dim7p, gdims = dim7, dimnames = dim7n, &
                code = 123, lpost = .FALSE., lrerun = .FALSE.)
           CALL add(IO_soil_wiso, 'wiso_csat_transp_fra', soil % wiso_csat_transp_fra, ldims = dim7p, gdims = dim7, &
                dimnames = dim7n, code = 124, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'remember_cveg', soil % remember_cveg, longname = 'Remember coverfraction Veg', &
                units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 127, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'remember_canopy', soil % remember_canopy, longname = 'Remember canopy res.', &
                units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 128, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'remember_rel_hum', soil % remember_rel_hum, longname = 'Remember relative humidity', &
                units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 129, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'remember_echamzchl', soil % remember_echamzchl, longname = 'Remember echam zchl', &
                units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 130, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'remember_cfrac', soil % remember_cfrac, longname = 'Remember cfrac', &
                units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 131, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'remember_snfrac', soil % remember_snfrac, longname = 'Remember snfrac', &
                units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 132, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'remember_wlfrac', soil % remember_wlfrac, longname = 'Remember wlfrac', &
                units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 133, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'remember_wind', soil % remember_wind, longname = 'Remember windspeed', &
                units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 134, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'remember_airm', soil % remember_airm, longname = 'Remember air moisture', &
                units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 135, lpost = .FALSE.) 
           CALL add(IO_soil_wiso, 'remember_satsphum', soil % remember_satsphum, longname = 'Remember satturation specific hum', &
                units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 136, lpost = .FALSE.)
           CALL add(IO_soil_wiso, 'q_wiso_Acoef_new', soil % q_wiso_Acoef_new, longname = 'q_wiso_Acoef_new', &
                units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 137, lpost = .FALSE., lrerun=.FALSE.) 
           CALL add(IO_soil_wiso, 'q_wiso_Bcoef_new', soil % q_wiso_Bcoef_new, longname = 'q_wiso_Bcoef_new', &
                units = '', ldims = dim3p, gdims = dim3, dimnames = dim3n, code = 138, lpost = .FALSE., lrerun=.FALSE.) 
        END IF
!---wiso-code-end

        CALL add(IO_diag, 'surface_temperature', soil_diag % surface_temperature, longname = 'Land Surface Temperature', &
        units = 'K', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .FALSE., code = 35, lpost = lpost_echam, &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'surface_radiative_temp', soil_diag % surface_radiative_temp, &
        longname = 'Surface Radiative Temperature', &
        units = 'K', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .FALSE., code = 36, lpost = .FALSE., &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'surface_qsat', soil_diag % sat_surface_specific_humidity, &
        longname = 'Soil Specific Humidity at Saturation', &
        units = 'kg/kg', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .FALSE., code = 39, lpost = .FALSE., &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'ground_heat_flux', soil_diag % ground_heat_flux_acc, longname = 'Ground Heat Flux (avg)', &
        units = 'W/m^2', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 40, &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'heat_capacity', soil_diag % heat_capacity, longname = 'Heat Capacity of the Uppermost Soil Layer', &
        units = 'J m-2 K-1', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 42, lpost = .FALSE., &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'evapotranspiration', soil_diag % evapotranspiration_acc, longname = 'Evapotranspiration (avg)', &
        units = 'kg/m^2s', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 44, &
        lmiss = .TRUE., missval = missing_value, bits = 24)
        CALL add(IO_diag, 'transpiration', soil_diag % transpiration_acc, longname = 'Transpiration (avg)', &
        units = 'kg/m^2s', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 76, lpost = .TRUE., &
        lmiss = .TRUE., missval = missing_value, bits = 24)
        CALL add(IO_diag, 'sensible_heat_flx', soil_diag % sensible_heat_acc, longname = 'Sensible Heat Flux (avg)', &
        units = 'W/m^2', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 47, lpost = lpost_echam, &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'latent_heat_flx', soil_diag % latent_heat_acc, longname = 'Latent Heat Flux (avg)', &
        units = 'W/m^2', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 49, lpost = lpost_echam, &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'soil_albedo', soil_diag % albedo, longname = 'Soil Albedo', &
        units = '', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .FALSE., code = 50, lpost = .FALSE., &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'soil_moisture', soil_diag % moisture, longname = 'Water Content of the Soil', &
        units = 'm', ldims = dim2p, gdims = dim2, dimnames = dim2n, laccu = .TRUE., code = 57, lpost = lpost_echam, &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'skin_reservoir', soil_diag % skin_reservoir, longname = 'Water Content of the Skin Reservoir', &
        units = 'm', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 55, lpost = lpost_echam, &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'snow', soil_diag % snow, longname = 'Snow Depth in Water Equivalent', &
        units = 'm', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 59, lpost = lpost_echam, &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'snow_acc', soil_diag % snow_acc, longname = 'Snow Depth (acc.)', &
        units = 'kg/m^2', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 61, lpost = lpost_echam, &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'snow_fract', soil_diag % snow_fract, longname = 'Snow Cover Fraction', &
        units = '', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .FALSE., code = 60, &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'snow_melt', soil_diag % snow_melt_acc, longname = 'Snow Melt', &
        units = 'kg/m^2s', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 62, lpost = lpost_echam, &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'runoff', soil_diag % runoff_acc, longname = 'Surface Runoff and Drainage', &
        units = 'kg/m^2s', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 63, lpost = lpost_echam, &
        lmiss = .TRUE., missval = missing_value, bits = 24)
        CALL add(IO_diag, 'drainage', soil_diag % drainage_acc, longname = 'Drainage', &
        units = 'kg/m^2s', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 64, lpost = lpost_echam, &
        lmiss = .TRUE., missval = missing_value, bits = 24)
        CALL add(IO_diag, 'glacier_depth', soil_diag % glacier_depth, longname = 'Glacier Depth', &
        units = 'm', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 65, lpost = lpost_echam, &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'glacier_precip_minus_evap', soil_diag % glacier_precip_minus_evap_acc, &
        longname = 'Precipitation minus Evaporation over Glaciers', &
        units = 'kg/m^2s', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 66, lpost = lpost_echam, &
        lmiss = .TRUE., missval = missing_value, bits = 24)
        CALL add(IO_diag, 'glacier_runoff', soil_diag % glacier_runoff_acc, longname = 'Glacier Runoff', &
        units = 'kg/m^2s', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 67, &
        lmiss = .TRUE., missval = missing_value)
        CALL add(IO_diag, 'soil_temperature', soil_diag % soil_temperature, longname = 'Soil Temperature', &
        units = 'K', ldims = dim6p, gdims = dim6, dimnames = dim6n, laccu = .FALSE., code = 68, &
        lmiss = .TRUE., missval = missing_value)

!---wiso-code
        IF (lwiso) THEN
           CALL add(IO_diag_wiso, 'snowglac', soil_diag % snowglac, longname = 'Snow Depth on Glaciers in Water Equivalent', &
                units = 'm', ldims = dim1p, gdims = dim1, dimnames = dim1n, laccu = .TRUE., code = 77, lpost = lpost_echam, &
                lmiss = .TRUE., missval = missing_value)
           CALL add(IO_diag_wiso, 'wiso_surface_qsat', soil_diag % wiso_sat_surface_specific_hum, &
                longname = 'Soil Specific Humidity at Saturation - Water Isotopes - Water Isotopes', units = 'kg/kg', &
                ldims = dim7p, gdims = dim7, dimnames = dim7n, laccu = .FALSE., code = 80, lpost = .FALSE., &
                lmiss = .TRUE., missval = missing_value)
           CALL add(IO_diag_wiso, 'wiso_evapotranspiration', soil_diag % wiso_evapotranspiration_acc, &
                longname = 'Evapotranspiration (avg) - Water Isotopes', units = 'kg m-2 s-1', &
                ldims = dim7p, gdims = dim7, dimnames = dim7n, laccu = .TRUE., code = 82, &
                lmiss = .TRUE., missval = missing_value, bits = 24)
           CALL add(IO_diag_wiso, 'wiso_transpiration', soil_diag % wiso_transpiration_acc, &
                longname = 'Transpiration(avg) - Water Isotopes', units = 'kg m-2 s-1', &
                ldims = dim7p, gdims = dim7, dimnames = dim7n, laccu = .TRUE., code = 85, lpost=.TRUE., &
                lmiss = .TRUE., missval = missing_value, bits = 24)
           CALL add(IO_diag_wiso, 'wiso_soil_moisture', soil_diag % wiso_moisture, &
                longname = 'Water Content of the Soil - Water Isotopes', units = 'm', &
                ldims = dim7p, gdims = dim7, dimnames = dim7n, laccu = .TRUE., code = 87, lpost = lpost_echam, &
                lmiss = .TRUE., missval = missing_value)
           CALL add(IO_diag_wiso, 'wiso_skin_reservoir', soil_diag % wiso_skin_reservoir, &
                longname = 'Water Content of the Skin Reservoir - Water Isotopes', units = 'm', &
                ldims = dim7p, gdims = dim7, dimnames = dim7n, laccu = .TRUE., code = 86, lpost = lpost_echam, &
                lmiss = .TRUE., missval = missing_value)
           CALL add(IO_diag_wiso, 'wiso_snow', soil_diag % wiso_snow, &
                longname = 'Snow Depth in Water Equivalent - Water Isotopes', units = 'm', &
                ldims = dim7, gdims = dim7, dimnames = dim7n, laccu = .TRUE., code = 89, lpost = lpost_echam, &
                lmiss = .TRUE., missval = missing_value)
           CALL add(IO_diag_wiso, 'wiso_snowglac', soil_diag % wiso_snowglac, &
                longname = 'Snow Depth on Glaciers - Water Isotopes', units = 'm', &
                ldims = dim7p, gdims = dim7, dimnames = dim7n, laccu = .TRUE., code = 78, lpost = lpost_echam, &
                lmiss = .TRUE., missval = missing_value)
           CALL add(IO_diag_wiso, 'wiso_snow_acc', soil_diag % wiso_snow_acc, &
                longname = 'Snow Depth (acc.) - Water Isotopes', units = 'kg m-2', &
                ldims = dim7p, gdims = dim7, dimnames = dim7n, laccu = .TRUE., code = 91, lpost = lpost_echam, &
                lmiss = .TRUE., missval = missing_value)
           CALL add(IO_diag_wiso, 'wiso_snow_melt', soil_diag % wiso_snow_melt_acc, &
                longname = 'Snow Melt - Water Isotopes', units = 'kg m-2 s-1', &
                ldims = dim7p, gdims = dim7, dimnames = dim7n, laccu = .TRUE., code = 92, lpost = lpost_echam, &
                lmiss = .TRUE., missval = missing_value)
           CALL add(IO_diag_wiso, 'wiso_runoff', soil_diag % wiso_runoff_acc, &
                longname = 'Surface Runoff and Drainage - Water Isotopes', units = 'kg m-2 s-1', &
                ldims = dim7p, gdims = dim7, dimnames = dim7n, laccu = .TRUE., code = 93, lpost = lpost_echam, &
                lmiss = .TRUE., missval = missing_value, bits = 24)
           CALL add(IO_diag_wiso, 'wiso_drainage', soil_diag % wiso_drainage_acc, &
                longname = 'Drainage - Water Isotopes', units = 'kg m-2 s-1', &
                ldims = dim7p, gdims = dim7, dimnames = dim7n, laccu = .TRUE., code = 94, lpost = lpost_echam, &
                lmiss = .TRUE., missval = missing_value, bits = 24)
           CALL add(IO_diag_wiso, 'wiso_glacier_depth', soil_diag % wiso_glacier_depth, &
                longname = 'Glacier Depth - Water Isotopes', units = 'm', &
                ldims = dim7p, gdims = dim7, dimnames = dim7n, laccu = .TRUE., code = 95, lpost = lpost_echam, &
                lmiss = .TRUE., missval = missing_value)
           CALL add(IO_diag_wiso, 'wiso_glac_p_minus_e', soil_diag % wiso_glac_p_minus_e_acc, &
                longname = 'Precipitation minus Evaporation over Glaciers - Water Isotopes', units = 'kg m-2 s-1', &
                ldims = dim7p, gdims = dim7, dimnames = dim7n, laccu = .TRUE., code = 96, lpost = lpost_echam, &
                lmiss = .TRUE., missval = missing_value, bits = 24)
           CALL add(IO_diag_wiso, 'wiso_glacier_runoff', soil_diag%wiso_glacier_runoff_acc, &
                longname = 'Glacier Runoff - Water Isotopes', units = 'kg m-2 s-1', &
                ldims = dim7p, gdims = dim7, dimnames = dim7n, laccu = .TRUE., code = 97, &
                lmiss = .TRUE., missval = missing_value)
        END IF
!---wiso-code

    END SUBROUTINE soil_init_memory
  !
  !=================================================================================================
  !
  SUBROUTINE init_soil(grid, domain, surface, soil_param, soil, IO_file_name, isRestart, useAlbedo, &
       useDynveg, fileformat, IO_diag_stream, IO_stream, IO_diag_wiso_stream, IO_wiso_stream)

    ! Initialize soil

    ! Called from jsbach_init

    USE mo_land_surface, ONLY: land_surface_type
    USE mo_jsbach_grid, ONLY: grid_type, domain_type
    USE mo_netcdf, ONLY: NF_NOERR, nf_inq_varid, nf_get_var_double, &
                         IO_inq_dimid, IO_inq_dimlen, IO_inq_varid, IO_get_var_double, nf_max_name
    USE mo_decomposition, ONLY: global_decomposition
    USE mo_transpose, ONLY: scatter_gp
    USE mo_io,    ONLY : IO_open, IO_READ, IO_close
    USE mo_temp                         ! Provided temporary arrays

!---wiso-code
    USE mo_time_control, ONLY: lresume
    USE mo_constants,    ONLY: tmelt
    USE mo_wiso,         ONLY: lwiso, lwiso_rerun, nwiso, tnat, twisosur1, twisosur2, nwisotyp, snglacmx
!---wiso-code-end

    TYPE(grid_type),         INTENT(in)    :: grid
    TYPE(domain_type),       INTENT(in)    :: domain
    TYPE(land_surface_type), INTENT(in)    :: surface
    TYPE(soil_param_type),   INTENT(inout) :: soil_param
    TYPE(soil_type),         INTENT(inout) :: soil
    CHARACTER(nf_max_name),  INTENT(in)    :: IO_file_name
    LOGICAL,                 INTENT(in)    :: isRestart
    LOGICAL,                 INTENT(in)    :: UseAlbedo
    LOGICAL,                 INTENT(in)    :: UseDynveg
    INTEGER,                 INTENT(in)    :: fileformat       ! output file format (grib/netcdf)
    TYPE(t_stream), POINTER            :: IO_diag_stream
    TYPE(t_stream), POINTER, OPTIONAL  :: IO_stream
!---wiso-code
    TYPE(t_stream), POINTER, OPTIONAL  :: IO_diag_wiso_stream
    TYPE(t_stream), POINTER, OPTIONAL  :: IO_wiso_stream
!---wiso-code-end
    INTEGER :: ntiles
    TYPE(FILE_INFO) :: IO_file
    INTEGER :: IO_file_id, IO_var_id, IO_dim_id
    INTEGER :: i, j, znlon, znlat, status
!---wiso-code
    INTEGER :: jt, itile, jl
!---wiso-code-end

    REAL(dp), POINTER :: fao(:)            !! FAO soil type flag [0...5]
    REAL(dp), ALLOCATABLE :: surface_temperature_clim(:,:)  !! Climatological values of surface temperature

!---wiso-code
    REAL(dp), ALLOCATABLE :: zwisosurf(:,:)
!---wiso-code-end

    IF (debug) CALL message('init_soil','')

    CALL config_soil

    ntiles = soil%ntiles

    CALL soil_init_io (IO_file_name)
    soil%nsoil  = nsoil
    soil%ntsoil = ntsoil

!---wiso-code
    IF (lwiso) THEN
       CALL soil_init_memory(grid%nland, domain%nland, ntiles, useDynveg, soil, fileformat, &
            IO_diag_stream, IO_diag_wiso_stream, stream=IO_stream, wiso_stream=IO_wiso_stream)
    ELSE
       CALL soil_init_memory(grid%nland, domain%nland, ntiles, useDynveg, soil, fileformat, &
            IO_diag_stream, stream=IO_stream)
    ENDIF
!---wiso-code-end

    IF (.NOT. ASSOCIATED(IO_soil)) &
         CALL finish('init_soil', 'No memory stream for soil')

    IF (p_parallel_io) THEN

       ! Open ini file
       WRITE(nout, *) 'Reading soil parameters from ', TRIM(soil_file%file_name)
       IO_file%opened = .FALSE.
       CALL IO_open(TRIM(soil_file%file_name), IO_file, IO_READ)
       IO_file_id = IO_file%file_id

       ! Check resolution
       CALL IO_inq_dimid  (IO_file_id, 'lat', IO_dim_id)
       CALL IO_inq_dimlen (IO_file_id, IO_dim_id, znlat)
       CALL IO_inq_dimid  (IO_file_id, 'lon', IO_dim_id)
       CALL IO_inq_dimlen (IO_file_id, IO_dim_id, znlon)

       IF (znlon /= grid%nlon .OR. znlat /= grid%nlat) THEN
          WRITE(nerr,*) 'Unexpected resolution:', znlon, znlat
          CALL finish('init_soil', 'Unexpected resolution')
       ENDIF

    ENDIF

    ! Read global soil parameters

    ! (nland,nsoil) variables

    ! Temporary storage for local domain fields
    ALLOCATE(zreal2d(domain%ndim,domain%nblocks))

    ! Soil layer depths
    ALLOCATE(soil_param%Depth(domain%nland, nsoil))
    IF (p_parallel_io) THEN
       ALLOCATE(zreal3d(grid%nlon,grid%nlat,nsoil))
!       CALL IO_inq_varid(IO_file_id, 'depth', IO_var_id)
!       CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       zreal3d = 1._dp  !! FOR TESTING
    ENDIF
    NULLIFY(zreal2d_ptr)
    DO i=1,nsoil
       IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
       CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
       soil_param%Depth(:,i) = PACK(zreal2d, MASK=domain%mask)
    ENDDO

    ! Field capacity
    ALLOCATE(soil_param%MaxMoisture(domain%nland,nsoil))
    IF (p_parallel_io) THEN
       CALL IO_inq_varid(IO_file_id, 'maxmoist', IO_var_id)
       CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
    END IF
    NULLIFY(zreal2d_ptr)
    DO i=1,nsoil
       IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
       CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
       soil_param%MaxMoisture(:,i) = PACK(zreal2d, MASK=domain%mask)
    ENDDO

    IF (soil_options%MaxMoistureLimit > 0.0_dp) THEN
       soil_param%MaxMoisture(:,:) = MIN(soil_param%MaxMoisture(:,:), soil_options%MaxMoistureLimit)
    END IF

    ! Ksat
!    ALLOCATE(soil_param%Ksat(grid%nland, nsoil))
!    IF (p_parallel_io) THEN
!       CALL IO_inq_varid(IO_file_id, 'ksat', IO_var_id)
!       CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
!    ENDIF
!    NULLIFY(zreal2d_ptr)
!    DO i=1,nsoil
!       IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
!       CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
!       soil_param%Ksat(:,i) = PACK(zreal2d, MASK=grid%mask_land)
!    ENDDO

    ! Initial soil moisture
    ALLOCATE(soil_param%InitMoisture(domain%nland,nsoil))
    IF (p_parallel_io) THEN
!       IF (nsoil == 1) THEN  !! ECHAM compatibility
          CALL IO_inq_varid(IO_file_id, 'init_moist', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
!       ELSE
!       zreal3d = 0.5_dp !! FOR TESTING ONLY
    ENDIF
    NULLIFY(zreal2d_ptr)
    DO i=1,nsoil
       IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
       CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
       soil_param%InitMoisture(:,i) = PACK(zreal2d, MASK=domain%mask)
    ENDDO

    IF (p_parallel_io) DEALLOCATE(zreal3d)

    ! (nland) variables

    ! Get FAO soil types
    ALLOCATE(fao(domain%nland))
    IF (p_parallel_io) THEN
       ALLOCATE(zzreal2d(grid%nlon,grid%nlat))
       CALL io_inq_varid(IO_file_id, 'fao', IO_var_id)
       CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
    END IF
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
    fao = PACK(zreal2d, MASK=domain%mask)

    ! Derive soil properties from FAO soil types
    ALLOCATE(soil_param%VolHeatCapacity(domain%nland))
    ALLOCATE(soil_param%ThermalDiffusivity(domain%nland))
    WHERE(NINT(fao) == 0) 
       soil_param%VolHeatCapacity = 2.25E+06_dp
       soil_param%ThermalDiffusivity = 7.4E-7_dp
    ELSEWHERE(NINT(fao) == 1) 
       soil_param%VolHeatCapacity = 1.93E+06_dp
       soil_param%ThermalDiffusivity = 8.7E-7_dp
    ELSEWHERE(NINT(fao) == 2)
       soil_param%VolHeatCapacity = 2.10E+06_dp
       soil_param%ThermalDiffusivity = 8.0E-7_dp
    ELSEWHERE(NINT(fao) == 3)
       soil_param%VolHeatCapacity = 2.25E+06_dp
       soil_param%ThermalDiffusivity = 7.4E-7_dp
    ELSEWHERE(NINT(fao) == 4)
       soil_param%VolHeatCapacity = 2.36E+06_dp
       soil_param%ThermalDiffusivity = 7.1E-7_dp
    ELSEWHERE(NINT(fao) == 5)
       soil_param%VolHeatCapacity = 2.48E+06_dp
       soil_param%ThermalDiffusivity = 6.7E-7_dp
    ELSEWHERE
       soil_param%VolHeatCapacity = -1.0_dp
       soil_param%ThermalDiffusivity = -1.0_dp
    END WHERE
    IF (ANY(soil_param%VolHeatCapacity < 0.0_dp)) &
         CALL finish('init_soil','FAO soil types must be between 0 and 5')
    DEALLOCATE(fao)

    ! Surface roughness length
    ALLOCATE(soil_param%Roughness(domain%nland))
    IF (p_parallel_io) THEN
       CALL IO_inq_varid(IO_file_id, 'roughness_length_oro', IO_var_id)
       CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
    ENDIF
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
    soil_param%Roughness = PACK(zreal2d, MASK=domain%mask)

    ! Ds
    ALLOCATE(soil_param%Ds(domain%nland))
!!$    IF (p_parallel_io) THEN
!!$       CALL IO_inq_varid(IO_file_id, 'ds', IO_var_id)
!!$       CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
!!$    ENDIF
!!$    NULLIFY(zreal2d_ptr)
!!$    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
!!$    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
!!$    soil_param%Ds = PACK(zreal2d, MASK=grid%mask_land)

    ! Average temperature at damping depth
    ALLOCATE(soil_param%AvgTemp(domain%nland))
!!$    IF (p_parallel_io) THEN
!!$       CALL IO_inq_varid(IO_file_id, 'avg_temp', IO_var_id)
!!$       CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
!!$    ENDIF
!!$    NULLIFY(zreal2d_ptr)
!!$    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
!!$    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
!!$    soil_param%AvgTemp = PACK(zreal2d, MASK=grid%mask_land)

    DEALLOCATE(zreal2d)
!!$    IF (p_parallel_io) THEN
!!$       DEALLOCATE(zzreal2d)
!!$       CALL IO_close(IO_file)
!!$    ENDIF

    soil_diag%surface_temperature = 0._dp
    soil_diag%surface_radiative_temp = 0._dp
    soil_diag%sat_surface_specific_humidity = 0._dp
    soil_diag%ground_heat_flux_acc = 0._dp
    soil_diag%heat_capacity = 0._dp
    soil_diag%evapotranspiration_acc = 0._dp
    soil_diag%transpiration_acc = 0._dp
    soil_diag%sensible_heat_acc = 0._dp
    soil_diag%latent_heat_acc = 0._dp
    soil_diag%albedo = 0._dp
    soil_diag%moisture = 0._dp
    soil_diag%skin_reservoir = 0._dp
    soil_diag%snow = 0._dp
!---wiso-code
    IF (lwiso) THEN
      soil_diag%snowglac = 0._dp
    END IF
!---wiso-code-end
    soil_diag%snow_acc = 0._dp
    soil_diag%snow_fract = 0._dp
    soil_diag%snow_melt_acc = 0._dp
    soil_diag%runoff_acc = 0._dp
    soil_diag%drainage_acc = 0._dp
    soil_diag%glacier_depth = 0._dp
    soil_diag%glacier_precip_minus_evap_acc = 0._dp
    soil_diag%glacier_runoff_acc = 0._dp
    soil_diag%soil_temperature = 0._dp

!---wiso-code
    IF (lwiso) THEN

    soil_diag%wiso_evapotranspiration_acc = 0._dp
    soil_diag%wiso_transpiration_acc = 0._dp
    soil_diag%wiso_sat_surface_specific_hum = 0._dp
    soil_diag%wiso_moisture = 0._dp
    soil_diag%wiso_skin_reservoir = 0._dp
    soil_diag%wiso_snow = 0._dp
    soil_diag%wiso_snowglac = 0._dp
    soil_diag%wiso_snow_acc = 0._dp
    soil_diag%wiso_snow_melt_acc = 0._dp
    soil_diag%wiso_runoff_acc = 0._dp
    soil_diag%wiso_drainage_acc = 0._dp
    soil_diag%wiso_glacier_depth = 0._dp
    soil_diag%wiso_glac_p_minus_e_acc = 0._dp
    soil_diag%wiso_glacier_runoff_acc = 0._dp
    
    END IF
!---wiso-code-end

    ALLOCATE(zreal2d(domain%ndim,domain%nblocks))

    IF (.NOT. UseAlbedo) THEN
    ! background albedo
    IF (p_parallel_io) THEN
       status = nf_inq_varid(IO_file_id,'albedo', IO_var_id)
       IF (status /= NF_NOERR) THEN
          CALL message('init_soil', 'Initializing background albedo')
          zzreal2d(:,:) = 0.0_dp
       ELSE
          status = nf_get_var_double(IO_file_id, IO_var_id, zzreal2d)
       ENDIF
    ENDIF
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
!!$    FORALL (i=1:ntiles) &
    DO i=1,ntiles
       soil%albedo(:,i) = PACK(zreal2d, MASK=domain%mask)   
    END DO
    ELSE
    ! albedo of vegetation and soil for NIR and visible range
    IF (p_parallel_io) THEN
       status = nf_inq_varid(IO_file_id,'albedo_soil_vis', IO_var_id)
       IF (status /= NF_NOERR) THEN
          CALL message('init_soil', 'Initializing albedo soil visible')
          zzreal2d(:,:) = 0.0_dp
       ELSE
          status = nf_get_var_double(IO_file_id, IO_var_id, zzreal2d)
       ENDIF
    ENDIF
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
    DO i=1,ntiles
       soil%albedo_soil_vis(:,i) = PACK(zreal2d, MASK=domain%mask)
       soil%albedo(:,i) = soil%albedo_soil_vis(:,i)
    END DO

    IF (p_parallel_io) THEN
       status = nf_inq_varid(IO_file_id,'albedo_soil_nir', IO_var_id)
       IF (status /= NF_NOERR) THEN
          CALL message('init_soil', 'Initializing albedo soil NIR')
          zzreal2d(:,:) = 0.0_dp
       ELSE
          status = nf_get_var_double(IO_file_id, IO_var_id, zzreal2d)
       ENDIF
    ENDIF
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
    DO i=1,ntiles
       soil%albedo_soil_nir(:,i) = PACK(zreal2d, MASK=domain%mask)
    END DO

    IF (p_parallel_io) THEN
       status = nf_inq_varid(IO_file_id,'albedo_veg_vis', IO_var_id)
       IF (status /= NF_NOERR) THEN
          CALL message('init_soil', 'Initializing albedo vegetation visible')
          zzreal2d(:,:) = 0.0_dp
       ELSE
          status = nf_get_var_double(IO_file_id, IO_var_id, zzreal2d)
       ENDIF
    ENDIF
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
    DO i=1,ntiles
       soil%albedo_vegetation_vis(:,i) = PACK(zreal2d, MASK=domain%mask)   
    END DO

    IF (p_parallel_io) THEN
       status = nf_inq_varid(IO_file_id,'albedo_veg_nir', IO_var_id)
       IF (status /= NF_NOERR) THEN
          CALL message('init_soil', 'Initializing albedo vegetation NIR')
          zzreal2d(:,:) = 0.0_dp
       ELSE
          status = nf_get_var_double(IO_file_id, IO_var_id, zzreal2d)
       ENDIF
    ENDIF
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
    DO i=1,ntiles
       soil%albedo_vegetation_nir(:,i) = PACK(zreal2d, MASK=domain%mask)   
    END DO
    END IF

!---wiso-code
    IF (lwiso) THEN

    ! For a restart without previous isotope diagnostics: 
    ! Initialize *wiso* variables from default JSBACH restart fields
    ! (assume at start a delta value of zero permill) 
    
    IF (lresume .AND. (.NOT. lwiso_rerun)) THEN
      DO jt=1,nwiso
        DO i=1,ntiles
          soil%snowglac(:,i) = snglacmx
          soil%wiso_snowglac(:,jt,i) = tnat(jt)*soil%snowglac(:,i)
          soil%wiso_evapotranspiration(:,jt,i) = tnat(jt)*soil%evapotranspiration(:,i)
          soil%wiso_evapotranspiration_acc(:,jt,i) = tnat(jt)*soil%evapotranspiration_acc(:,i)
          soil%wiso_evaporation_pot(:,jt,i) = tnat(jt)*soil%evaporation_pot(:,i)
          soil%wiso_transpiration(:,jt,i) = tnat(jt)*soil%transpiration(:,i)
          soil%wiso_transpiration_acc(:,jt,i) = tnat(jt)*soil%transpiration_acc(:,i)
          soil%wiso_skin_reservoir(:,jt,i) = tnat(jt)*soil%skin_reservoir(:,i)
          soil%wiso_moisture(:,jt,i) = tnat(jt)*soil%moisture(:,1,i)
          soil%wiso_snow(:,jt,i) = tnat(jt)*soil%snow(:,i)
          soil%wiso_snow_acc(:,jt,i) = tnat(jt)*soil%snow_acc(:,i)
          soil%wiso_snow_melt_acc(:,jt,i) = tnat(jt)*soil%snow_melt_acc(:,i)
          soil%wiso_runoff_acc(:,jt,i) = tnat(jt)*soil%runoff_acc(:,i)
          soil%wiso_drainage_acc(:,jt,i) = tnat(jt)*soil%drainage_acc(:,i)
          soil%wiso_glacier_depth(:,jt,i) = tnat(jt)*soil%glacier_depth(:,i)
          soil%wiso_glac_p_minus_e_acc(:,jt,i) = tnat(jt)*soil%glacier_precip_minus_evap_acc(:,i)
          soil%wiso_glacier_runoff_acc(:,jt,i) = tnat(jt)*soil%glacier_runoff_acc(:,i)
          
          soil%wiso_csat(:,jt) = tnat(jt)*soil%csat(:)
          soil%wiso_cair(:,jt) = tnat(jt)*soil%cair(:) 
          soil%wiso_csat_transpiration(:,jt) = tnat(jt)*soil%csat_transpiration(:)          
          soil%wiso_csat_transp_fra(:,jt) = 0._dp
          
          soil%remember_cveg(:,i) = 0._dp
          soil%remember_canopy(:,i) = 0._dp
          soil%remember_rel_hum(:,i) = 0._dp
          soil%remember_echamzchl(:,i) = 0._dp
          soil%remember_cfrac(:,i) = 0._dp
          soil%remember_snfrac(:,i) = 0._dp
          soil%remember_wlfrac(:,i) = 0._dp
          soil%remember_wind(:,i) = 0._dp
          soil%remember_airm(:,i) = 0._dp
          soil%remember_satsphum(:,i) = 0._dp
        END DO
      END DO
    END IF
    
    END IF
!---wiso-code-end

    ! If this is a restart run we can exit now since the model state variables are read from restart file
    IF (isRestart) THEN
       DEALLOCATE(zreal2d)
       IF (p_parallel_io) DEALLOCATE(zzreal2d)
       RETURN
    ENDIF

    ! If this is not a restart run initialize soil model state from values read from initialization file

    ! Set initial value for soil moisture
!!$    FORALL (i=1:nsoil, j=1:ntiles) &
    DO i=1,nsoil
       DO j=1,ntiles
          soil%moisture(:,i,j) = soil_param%InitMoisture(:,i)
       END DO
    END DO

    ! Initial snow depth
    IF (p_parallel_io) THEN
       status = nf_inq_varid(IO_file_id, 'snow', IO_var_id)
       IF (status /= NF_NOERR) THEN
          CALL message('init_soil', 'Initializing snow to zero')
          zzreal2d(:,:) = 0.0_dp
       ELSE
          status = nf_get_var_double(IO_file_id, IO_var_id, zzreal2d)
       ENDIF
    ENDIF
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
!!$    FORALL (i=1:ntiles) &
    DO i=1,ntiles
       soil%snow(:,i) = PACK(zreal2d, MASK=domain%mask)
       ! Roesch et al, 2002, Climate Dynamics (compare with snow_fract calculation in update_soil)
       WHERE (surface%is_glacier(:,i))
         soil%snow_fract(:,i) = 1.0_dp
       ELSEWHERE
         soil%snow_fract(:,i) = zqsncr * TANH(soil%snow(:,i) * 100._dp) *            &
            SQRT(soil%snow(:,i) * 1000._dp / (soil%snow(:,i) * 1000._dp + zepsec +   &
            zsigfac * surface%oro_std_dev(:)))
       ENDWHERE
    END DO

!---wiso-code
    IF (lwiso) THEN

    ! Initial snow depth on glaciers
    IF (p_parallel_io) zzreal2d(:,:) = snglacmx
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
!!$    FORALL (i=1:ntiles) &
    DO i=1,ntiles
       soil%snowglac(:,i) = PACK(zreal2d, MASK=domain%mask)
    END DO

    END IF
!---wiso-code-end

    ! Calculate fraction of roots in each soil layer
!!$    FORALL(i=1:grid%nland,j=1:npft, vegetation_cover(i,j) > 0.0_dp) &
!!$       soil_root_fract(i,:,j) = calc_root_fractions(root_depth(i,:,j), &
!!$                                                    root_fract(i,:,j), &
!!$                                                    soil_param%Depth(i,:))
!!$    IF (ANY(soil_root_fract(:,1,:)==0.0_dp)) THEN
!!$       CALL finish('JSBACH - calc_root_fractions', &
!!$                   'Root fractions sum equals zero')
!!$    END IF

    ! Climatological values of surface temperature
    ALLOCATE (surface_temperature_clim(domain%nland,12))

    IF (p_parallel_io) THEN
       ALLOCATE(zreal3d(grid%nlon,grid%nlat,12))
       CALL IO_inq_varid(IO_file_id, 'surf_temp', IO_var_id)
       CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
    ENDIF
    NULLIFY(zreal2d_ptr)
    DO i=1,12
       IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
       CALL scatter_gp(zreal2d_ptr, zreal2d, global_decomposition)
       surface_temperature_clim(:,i) = PACK(zreal2d, MASK=domain%mask)
    ENDDO
    IF (p_parallel_io) DEALLOCATE(zreal3d)

    ! Initialize soil and surface temperature
    soil%skin_reservoir = 0._dp

    CALL ini_soil_temp(domain%nland, ntsoil, &
         surface_temperature_clim(:,:), &
         soil%soil_temperature(:,:,1), &
         soil%surface_temperature(:,1))

    DEALLOCATE (surface_temperature_clim)

    IF (ntiles > 1) THEN
      DO i=2,ntiles
        soil%soil_temperature(:,:,i) = soil%soil_temperature(:,:,1)
        soil%surface_temperature(:,i) = soil%surface_temperature(:,1)
      END DO
    ENDIF
    soil%surface_temperature_old        = soil%surface_temperature
    soil%surface_temperature_unfiltered = soil%surface_temperature

    soil%latent_heat_flux = 0._dp
    soil%latent_heat_acc = 0._dp
    soil%sensible_heat_flux = 0._dp
    soil%sensible_heat_acc = 0._dp

    soil%heat_capacity = 0._dp
    soil%ground_heat_flux = 0._dp
    soil%ground_heat_flux_acc = 0._dp

    soil%snow_acc = 0._dp
    soil%snow_melt_acc = 0._dp
    soil%glacier_depth = 0._dp
    soil%glacier_precip_minus_evap_acc = 0._dp
    soil%runoff_acc = 0._dp
    soil%drainage_acc = 0._dp
    soil%glacier_runoff_acc = 0._dp

!---wiso-code
    IF (lwiso) THEN
    
    ! zwisosurf = inital_deviation_1 from SMOW at surface * surf.temp. - initial_deviation_2
    ALLOCATE (zwisosurf(domain%nland, nwiso))
    
    DO jt=1,nwiso
       zwisosurf(:,jt)=twisosur1(jt)*(soil%surface_temperature(:,1) - tmelt) - twisosur2(jt)
       IF (nwisotyp(jt).eq.2) THEN  ! if tracer = 18o: zwisosurf = min (-4./1000., zwisosurf)
           zwisosurf(:,jt)=min( -4._dp/1000._dp,zwisosurf(:,jt))
       ELSE IF (nwisotyp(jt).eq.3) THEN  ! if tracer = hdo: zwisosurf = min (-22./1000., zwisosurf)
           zwisosurf(:,jt)=min(-22._dp/1000._dp,zwisosurf(:,jt))
       ENDIF
    END DO

    DO jt=1,nwiso
       DO i=1,ntiles
         ! set tracer reservoirs = nat.isotope-ratio * corresponding wetness reservoir * (1.+zwisosurf)
         soil%wiso_moisture(:,jt,i) = soil % moisture(:,1,i) * tnat(jt) * (1._dp + zwisosurf(:,jt))
         soil%wiso_snow(:,jt,i)     = soil % snow(:,i) * tnat(jt)  * (1._dp + zwisosurf(:,jt))
         soil%wiso_snowglac(:,jt,i) = soil % snowglac(:,i) * tnat(jt)  * (1._dp + zwisosurf(:,jt))
         
         ! maximum ws value for any tracer: wsmx (isotope-independent)
         soil%wiso_moisture(:,jt,i) = MIN(soil%wiso_moisture(:,jt,i),soil_param%MaxMoisture(:,1)) 
        END DO
     END DO

    soil%wiso_skin_reservoir = 0._dp
    soil%wiso_snow_acc = 0._dp
    soil%wiso_snow_melt_acc = 0._dp
    soil%wiso_glacier_depth = 0._dp
    soil%wiso_glac_p_minus_e_acc = 0._dp
    soil%wiso_runoff_acc = 0._dp
    soil%wiso_drainage_acc = 0._dp
    soil%wiso_glacier_runoff_acc = 0._dp
    
    DEALLOCATE (zwisosurf)

    END IF
!---wiso-code-end

    IF (useDynveg) soil%relative_humidity_air = 0._dp

    IF (p_parallel_io) THEN
       DEALLOCATE(zzreal2d)
    END IF
    DEALLOCATE(zreal2d)

  END SUBROUTINE init_soil

  SUBROUTINE update_soil(kidx, surface, soil, soil_param, useDynveg, canopy_conductance_max, root_depth,  &
        lai, canopy_snow, canopy_snow_fract,                                                                &
        cdrag, t_Acoef, t_Bcoef, q_Acoef, q_Bcoef, air_temperature, air_moisture,                           &
        surface_pressure, windspeed, wind10, rad_longwave_down, rad_shortwave_net, rain, snow,              &
        cair, csat, p_echam_zchl, zhsoil_avg, tte_corr_avg, canopy_conductance_limited,                     &
        glac_runoff_evap, surf_runoff_hd, drainage_hd,                                                      &
!---wiso-code
        lwisofracl,                                                                                         &
        wiso_canopy_snow, q_wiso_Acoef, q_wiso_Bcoef,                                                       &
        wiso_air_moisture, wisokinl,                                                                        &
        wiso_nenner,                                                                                        &
        wiso_helpqdif,                                                                                      &
        wiso_rain, wiso_snow, wiso_cair, wiso_csat, wiso_cair_fra, wiso_csat_fra,                           &
        q_wiso_Acoef_new, q_wiso_Bcoef_new,                                                                 &
        wiso_glac_runoff_evap, wiso_surf_runoff_hd, wiso_drainage_hd)
!---wiso-code-end

        USE mo_land_surface,     ONLY: land_surface_type
        USE mo_jsbach_constants, ONLY: Gravity, RhoH2O, SpecificHeatDryAirConstPressure,                 &
                                       SpecificHeatVaporConstPressure, Emissivity, StefanBoltzmann,      &
                                       LatentHeatVaporization, LatentHeatSublimation, GasConstantDryAir, &
                                       GasConstantWaterVapor
        USE mo_constants,        ONLY: tmelt
        USE mo_atmosphere,       ONLY: sat_specific_humidity
        USE mo_utils,            ONLY: average_tiles

#ifndef STANDALONE
        USE mo_physc2,           ONLY: cvdifts ! Factor for time step weighting in ECHAM5
#endif
        USE mo_time_control,     ONLY: lstart, delta_time, time_step_len
        USE mo_semi_impl,        ONLY: eps
        USE mo_param_switches,   ONLY: lsurf

!---wiso-code
        USE mo_time_control,     ONLY: lresume
        USE mo_wiso,             ONLY: lwiso, lwiso_rerun, nwiso, tnat, twisorhoh2o, snglacmx, &
                                       talphal1, talphal2, talphal3, talphal4, &
                                       talphas1, talphas2, talphas3, talphas4, &
                                       nwisotyp,                     &
                                       tsatbase, tsatfac, tdifrel,   &
                                       tkinfa1, tkinfa2, tkinsl,     &
                                       cwisomin, cwisosec 
!---wiso-code-end

        INTEGER,                 INTENT(in)     :: kidx

        TYPE(land_surface_type), INTENT(in)     :: surface
        TYPE(soil_param_type),   INTENT(in)     :: soil_param

        TYPE(soil_type),         INTENT(inout)  :: soil
        
        LOGICAL,                 INTENT(in)     :: useDynveg

        REAL(dp), DIMENSION(:,:), INTENT(in) :: &                           ! Dimension (nidx,ntiles)
                                                lai                         ! leaf aria index
        REAL(dp), DIMENSION(:,:), INTENT(in) :: &                           ! Dimension (nidx,ntiles)
                                                canopy_conductance_max      ! Unstressed canopy resistance


        REAL(dp), DIMENSION(:),   INTENT(in) :: &                           ! Dimension nidx
                                                cdrag, &
                                                t_Acoef, t_Bcoef, &         ! land%zetnl, land%zftnl - Richtmyer-Morten-Coef.
                                                q_Acoef, q_Bcoef, &         ! land%zeqnl, land%zfqnl - Richtmyer-Morten-Coef.
                                                air_temperature, &          ! Temperature at lowest atmospheric level
                                                air_moisture, &             ! Specific humidity at lowest atmospheric level
                                                surface_pressure, &         ! Surface pressure
                                                windspeed, &                ! Wind speed
                                                wind10, &                   ! 10m wind speed
                                                rad_longwave_down, &        ! Longwave radiation down
                                                rad_shortwave_net, &        ! Net shortwave radiation
                                                rain, snow                  ! Rain and snow [m water equivalent]
                                                
!---wiso-code
        INTEGER,                  INTENT(in), OPTIONAL :: lwisofracl
        REAL(dp), DIMENSION(:),   INTENT(in), OPTIONAL :: &                           ! Dimension nidx
                                                wiso_helpqdif

        REAL(dp), DIMENSION(:,:), INTENT(in), OPTIONAL :: &                              ! Dimension nidx, nwiso
                                                wiso_air_moisture,          &  ! humidity at lowest level - Isotopes
                                                q_wiso_Acoef, q_wiso_Bcoef, &  ! land%zeqnl, land%zfqnl - Ri.-Mor.-Coef.
                                                wiso_nenner,                & 
                                                wisokinl,                   &  ! kinetic fractionation factor over land
                                                wiso_rain, wiso_snow           ! Rain and snow [m water equivalent]
!---wiso-code-end

        REAL(dp), DIMENSION(:),   INTENT(in) :: p_echam_zchl
        REAL(dp), INTENT(in) ::                 root_depth(:,:,:)

        REAL(dp), DIMENSION(:,:), INTENT(inout) :: &                        ! Dimension (nidx,ntiles)
                                                canopy_snow, &              ! Snow depth in canopy
                                                canopy_snow_fract
                                                
!---wiso-code
        REAL(dp), DIMENSION(:,:,:), INTENT(inout), OPTIONAL :: &                      ! Dimension (nidx,nwiso,ntiles)
                                                wiso_canopy_snow            ! Snow depth in canopy - isotopes
!---wiso-code-end

        REAL(dp), DIMENSION(:,:), INTENT(out) :: &                          ! Dimension (nidx,ntiles)
                                                canopy_conductance_limited  ! resistance lim. by water avaiability (water stress)
        
        REAL(dp), DIMENSION(:),   INTENT(out) :: &                          ! Dimension nidx
                                                cair, csat, zhsoil_avg, tte_corr_avg, &
                                                ! INPUT for HD-Model in coupled ocean model
                                                glac_runoff_evap, surf_runoff_hd, drainage_hd
                                         
!---wiso-code
        REAL(dp), DIMENSION(:,:), INTENT(out), OPTIONAL :: &                          ! Dimension (nidx,nwiso)
                                                 wiso_cair, wiso_csat, &
                                                 wiso_cair_fra, wiso_csat_fra, &
                                                 q_wiso_Acoef_new, q_wiso_Bcoef_new, &
                                                ! INPUT for HD-Model in coupled ocean model
                                                wiso_glac_runoff_evap, wiso_surf_runoff_hd, wiso_drainage_hd
!---wiso-code-end                                         

        ! Local variables
        REAL(dp), DIMENSION(kidx) :: &
                           net_radiation, & ! Net radiation at surface
                           dry_static_energy, & ! Surface dry static energy (= C_p * T )
                           dry_static_energy_new, & ! New dry static energy
                           surface_qsat_new, & ! New surface saturated specific humidity
                           air_qsat, & ! Saturated specific humidity at lowest atmospheric level
                           evapotranspiration_no_snow_skin, & ! Evapotranspiration without that from snow and the skin reservoir
                           zdqsl, & ! Sensitivity of saturated surface specific humidity wrt temperature
                           zcpq, & ! Conversion factor for humidity from dry static energy
                           snow_avg, & ! Snow [m water equivalent] averaged over all tiles
!---wiso-code
                           snowglac_avg, & ! Snow on glaciers [m water equivalent] averaged over all tiles
!---wiso-code-end
                           melt_water_excess, & ! water from snow melting which exceeds skin reservoir and infiltrates soil
                           csat_transpiration   ! fraction of grid box that contributes to transpiration
                                                ! (considered to be completely wet)
                                   
        REAL(dp), DIMENSION(kidx) :: &
                           sat_surf_specific_hum_avg, &
                           ground_heat_flux_avg, heat_capacity_avg, surface_temperature_avg, &
                           soil_moisture_avg

        REAL(dp), DIMENSION(kidx, soil % ntiles) :: &
                           skin_reservoir_max, &
                           wet_skin_fract, &
                           glacier_precip_minus_evap, & ! P-E for glaciers [m]
                           surface_runoff, drainage, & ! Surface runoff and drainage for HD model [m]
                           qsat_fact, qair_fact, & ! Factors for implicit coupling
                           air_dry_static_energy_new, &
                           air_moisture_new, &
                           evaporation_pot, & ! Potential evaporation
                           canopy_resistance, & ! Water-limited canopy resistance
                           soil_moisture_root, & ! Soil moisture in root zone
                           soil_moisture_root_max, & ! Field capacity in root zone
                           water_stress_factor, & ! Water stress factor (=1: no stress, =0: infinite stress)
                           relative_humidity, & ! Relative humidity (Eq. 3.3.2.9 ECHAM3 Manual)
                           zhsoil, tte_corr, &
                           sat_surface_specific_hum_old, &
                           qsat_veg, qair_veg, &
                           csat_tiles, cair_tiles, &
                           qsat_transpiration, csat_transpiration_tiles

!---wiso-code
        REAL(dp), DIMENSION(kidx,nwiso) :: wiso_evapotrans_no_snow_skin,  &
                                           wiso_snow_avg, 	              &
                                           wiso_snowglac_avg, 	          &
                                           wiso_csat_transpiration, 	  &
                                           wiso_csat_transp_fra, 	      &
                                           wiso_soil_moisture_avg,        &
                                           wiso_melt_water_excess,        &    
                                           zdelta_air

        REAL(dp), DIMENSION(kidx, nwiso, soil % ntiles) ::           &
                           wiso_glac_p_minus_e,           & ! P-E for glaciers [m]
                           wiso_surface_runoff, wiso_drainage,       & ! Surface runoff and drainage for HD model [m]
                           qwisosat_fact, qwisoair_fact,             & ! Factors for implicit coupling
                           qwisosat_fact_fra, qwisoair_fact_fra,     & 
                           wiso_air_moisture_new,                    &
                           wiso_evaporation_pot,                     & ! Potential evaporation
                           qwisosat_veg, qwisoair_veg,               &
                           qwisosat_veg_fra, qwisoair_veg_fra,       &
                           wiso_csat_tiles, wiso_cair_tiles,           &
                           wiso_csat_tiles_fra, wiso_cair_tiles_fra,   &
                           qwisosat_transpiration, wiso_csat_transpiration_tiles, &
                           qwisosat_transp_fra, wiso_csat_transp_tiles_fra, &
                           zdelta_snow, zdelta_soil, zdelta_skin, zdelta_air_tiles, &  
                           zwisofracl, zwisofraci,                   &
                           wisoqklevl_tiles

        REAL(dp), DIMENSION(kidx, soil % ntiles) :: zsatval

        REAL(dp) :: zwisosurfmin, zwisoroundoff, zwspeedmin
        REAL(dp) :: zdelta, zwisodenom     
!---wiso-code-end                   

        LOGICAL, DIMENSION(kidx, soil % ntiles) :: &
                           soil_mask ! True if not glacier and not lake

        REAL(dp) :: zcons15, zcons30, ztpfac2, ztpfac3, vtmpc2

    INTEGER :: ntiles, kidx0, kidx1, itile, nroot_zones, i, j, jt, jl

#ifdef STANDALONE
    REAL(dp), PARAMETER :: cvdifts = 1.0_dp
#endif

    EXTERNAL update_surfacetemp

    ntiles = soil%ntiles
    kidx0 = kstart
    kidx1 = kend
    nroot_zones = SIZE(root_depth,2)

    zcons15 = 1._dp / (Gravity*time_step_len)
    zcons30 = 1._dp / (cvdifts*Gravity*time_step_len)
    ztpfac2 = 1._dp / cvdifts
    ztpfac3 = 1._dp - ztpfac2
    vtmpc2 = SpecificHeatVaporConstPressure / SpecificHeatDryAirConstPressure - 1._dp

    soil_mask = .NOT. surface%is_glacier(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) &
         .AND. surface%is_present(kidx0:kidx1,:)

!---wiso-code
    IF (lwiso) THEN

    ! Minimum Windspeed for calculation of kinetic fractionation effects
    zwspeedmin=7._dp
    ! Securityparameters for calculation of delta values
    zwisosurfmin   = 1.e-12_dp
    zwisoroundoff  = 1.e-8_dp
    
    END IF
!---wiso-code-end

    !----------------------------------------------------------------------------------------

    IF (lstart) THEN
       csat(1:nidx) = .5_dp
       cair(1:nidx) = .5_dp
       csat_transpiration(1:nidx) = .5_dp

!---wiso-code
    IF (lwiso) THEN

    ! calculate the isotope ratios for snow-, soil- and skinreservoir
    DO jt=1,nwiso
       wiso_csat(:,jt)      = .5_dp * tnat(jt)
       wiso_cair(:,jt)      = .5_dp * tnat(jt) 
       wiso_csat_transpiration(:,jt) = .5_dp * tnat(jt)
       wiso_csat_fra(:,jt)           = 0._dp
       wiso_cair_fra(:,jt)           = 0._dp
       wiso_csat_transp_fra(:,jt)    = 0._dp
    END DO
    
    END IF
!---wiso-code-end

    ELSE
       csat(1:nidx) = soil % csat(kidx0:kidx1)
       cair(1:nidx) = soil % cair(kidx0:kidx1)
       csat_transpiration(1:nidx) = soil % csat_transpiration(kidx0:kidx1)

!---wiso-code
    IF (lwiso) THEN

    ! For the calculation of wiso_csat, wiso_cair qm1 and wisoqm1 are needed. So calculation for isotopes moved from the end
    ! of previous timestep to beginning of actuell timestep!     

    qwisosat_fact(:,:,:)          = 0._dp
    qwisoair_fact(:,:,:)          = 0._dp
    qwisosat_veg(:,:,:)           = 0._dp
    qwisoair_veg(:,:,:) 	      = 0._dp
    qwisosat_fact_fra(:,:,:)      = 0._dp
    qwisoair_fact_fra(:,:,:)      = 0._dp
    qwisosat_veg_fra(:,:,:)       = 0._dp
    qwisoair_veg_fra(:,:,:) 	  = 0._dp
    qwisosat_transpiration(:,:,:) = 0._dp
    qwisosat_transp_fra(:,:,:)    = 0._dp
 
    DO jt=1,nwiso
       zdelta_air(:,jt)=tnat(jt)
       WHERE ((wiso_air_moisture(1:nidx,jt).gt.zwisosurfmin).and.(air_moisture(1:nidx).gt.zwisosurfmin))      
          zdelta_air(:,jt)=(wiso_air_moisture(1:nidx,jt) / air_moisture(1:nidx))
       END WHERE 
       ! cut off rounding errors
       WHERE (abs(1._dp-zdelta_air(:,jt)).lt.zwisoroundoff) 
          zdelta_air(:,jt)=1_dp
       END WHERE
       !END IF
       DO itile=1,ntiles
          DO i=1,nidx
             j = kidx0 + i - 1

             ! calculate isotope ratio of snow reservoir (distinguish between land and glaciers)
             zdelta_snow(i,jt,itile)=tnat(jt)
             IF (.NOT. surface%is_glacier(j,itile)) THEN
               IF ((soil%wiso_snow(j,jt,itile).gt.zwisosurfmin).and.(soil%snow(j,itile).gt.zwisosurfmin))         &
                  zdelta_snow(i,jt,itile)=(soil%wiso_snow(j,jt,itile)/soil%snow(j,itile))
             ELSE
               IF ((soil%wiso_snowglac(j,jt,itile).gt.zwisosurfmin).and.(soil%snowglac(j,itile).gt.zwisosurfmin)) &
                  zdelta_snow(i,jt,itile)=(soil%wiso_snowglac(j,jt,itile)/soil%snowglac(j,itile))
             END IF
             ! cut off rounding errors
             IF (abs(1._dp-zdelta_snow(i,jt,itile)).lt.zwisoroundoff) zdelta_snow(i,jt,itile)=1._dp

             ! calculate isotope ratio of soil reservoir
             zdelta_soil(i,jt,itile)=tnat(jt)
             IF ((soil%wiso_moisture(j,jt,itile).gt.zwisosurfmin).and.(soil%moisture(j,1,itile).gt.zwisosurfmin))&
                  zdelta_soil(i,jt,itile)=(soil%wiso_moisture(j,jt,itile)/soil%moisture(j,1,itile))
             ! cut off rounding errors
             IF (abs(1._dp-zdelta_soil(i,jt,itile)).lt.zwisoroundoff) zdelta_soil(i,jt,itile)=1._dp

             ! calculate isotope ratio of skin reservoir
             zdelta_skin(i,jt,itile)=tnat(jt)
             IF ((soil%wiso_skin_reservoir(j,jt,itile).gt.zwisosurfmin)                                          &
                  .and.(soil%skin_reservoir(j,itile).gt.zwisosurfmin))                                           &
                 zdelta_skin(i,jt,itile)=(soil%wiso_skin_reservoir(j,jt,itile)/soil%skin_reservoir(j,itile))
             ! cut off rounding errors
             IF (abs(1._dp-zdelta_skin(i,jt,itile)).lt.zwisoroundoff) zdelta_skin(i,jt,itile)=1._dp

             ! calculate isotope ratio of air_moisture
             zdelta_air_tiles(i,jt,itile)=zdelta_air(i,jt)
          END DO
       END DO
    END DO

    ! This loop should be equivalent to the first part of loop 331 of ECHAM5
    ! At this part zwisowet is calculated, multiplied by isotope ratios
    DO jt = 1, nwiso
       DO itile = 1, ntiles
          DO i = 1, nidx
             j = kidx0 + i - 1
             IF (lwisofracl==2) THEN ! fractionation of E and T
             IF (surface % is_present(j, itile)) THEN
                IF (soil % moisture(j, 1, itile)>                                                                   &
                   (soil_options % MoistureFractWilting * soil_param % MaxMoisture(j, 1))) THEN

                      qwisosat_veg(i,jt,itile) = soil%remember_snfrac(j,itile) * zdelta_snow(i,jt,itile)            &
                      + (1._dp - soil%remember_snfrac(j,itile)) * ( soil%remember_wlfrac(j,itile) * zdelta_skin(i,jt,itile) )
                      qwisoair_veg(i,jt,itile) = qwisosat_veg(i,jt,itile)
                      qwisosat_transpiration(i,jt,itile) = 0._dp
                      qwisosat_veg_fra(i,jt,itile) = (1._dp - soil%remember_snfrac(j,itile))                        &
                                                   * (1._dp - soil%remember_wlfrac(j,itile))                        &
                         / ( 1._dp + soil%remember_echamzchl(j,itile) * soil%remember_canopy(j,itile)                &
                                                                     * MAX(1.0_dp,soil%remember_wind(j,itile)) )
                      qwisoair_veg_fra(i,jt,itile) = qwisosat_veg_fra(i,jt,itile) 
                      qwisosat_transp_fra(i,jt,itile) = qwisosat_veg_fra(i,jt,itile)

                ELSE

                      qwisosat_veg(i,jt,itile) = soil%remember_snfrac(j,itile) * zdelta_snow(i,jt,itile)              &
                      + (1._dp - soil%remember_snfrac(j,itile)) *soil%remember_wlfrac(j,itile) * zdelta_skin(i,jt,itile)
                      qwisosat_transpiration(i,jt,itile) = 0._dp
                      qwisoair_veg(i,jt,itile) = qwisosat_veg(i,jt,itile)
                      qwisoair_veg_fra(i,jt,itile) = 0._dp 
                      qwisosat_veg_fra(i,jt,itile) = 0._dp
                      qwisosat_transp_fra(i,jt,itile) = 0._dp

                END IF
             END IF
             ELSE
             IF (surface % is_present(j, itile)) THEN
                IF (soil % moisture(j, 1, itile)>                                                                   &
                   (soil_options % MoistureFractWilting * soil_param % MaxMoisture(j, 1))) THEN

                      qwisosat_veg(i,jt,itile) = soil%remember_snfrac(j,itile) * zdelta_snow(i,jt,itile)            &
                      + (1._dp - soil%remember_snfrac(j,itile)) * (soil%remember_wlfrac(j,itile) * zdelta_skin(i,jt,itile)&
                      + (1._dp - soil%remember_wlfrac(j,itile)) * zdelta_soil(i,jt,itile)                           &
                         / (1._dp + soil%remember_echamzchl(j,itile) * soil%remember_canopy(j,itile)                &
                                                                     * MAX(1.0_dp,soil%remember_wind(j,itile))))
                      qwisosat_transpiration(i,jt,itile) = (1._dp - soil%remember_snfrac(j,itile))                  &
                      * (1._dp - soil%remember_wlfrac(j,itile)) * zdelta_soil(i,jt,itile)                           &
                         / (1._dp + soil%remember_echamzchl(j,itile) * soil%remember_canopy(j,itile)                &
                                                                     * MAX(1.0_dp,soil%remember_wind(j,itile)))
                      qwisoair_veg(i,jt,itile) = qwisosat_veg(i,jt,itile)
                      qwisoair_veg_fra(i,jt,itile) = 0._dp 
                      qwisosat_veg_fra(i,jt,itile) = 0._dp
                      qwisosat_transp_fra(i,jt,itile) = 0._dp

                ELSE

                      qwisosat_veg(i,jt,itile) = soil%remember_snfrac(j,itile) * zdelta_snow(i,jt,itile)              &
                      + (1._dp - soil%remember_snfrac(j,itile)) *soil%remember_wlfrac(j,itile) * zdelta_skin(i,jt,itile)
                      qwisosat_transpiration(i,jt,itile) = 0._dp
                      qwisoair_veg(i,jt,itile) = qwisosat_veg(i,jt,itile)
                      qwisoair_veg_fra(i,jt,itile) = 0._dp 
                      qwisosat_veg_fra(i,jt,itile) = 0._dp
                      qwisosat_transp_fra(i,jt,itile) = 0._dp

                END IF
             END IF
             END IF

            ! This loop should be equivalent to the second part of loop 331 of ECHAM5
            ! At this part is the first calculation of zcsat and zcair
             IF (lwisofracl==0) THEN 
             IF (soil%remember_satsphum(j,itile) .NE. 0._dp) THEN
               IF (soil%remember_rel_hum(j,itile) > soil % remember_airm(j,itile) / &
                 soil%remember_satsphum(j,itile) .AND. soil%remember_rel_hum(j,itile) > 1.e-10_dp ) THEN
 
                qwisosat_fact(i,jt,itile)=soil%remember_snfrac(j,itile)*zdelta_snow(i,jt,itile)                            &
                   +(1._dp-soil%remember_snfrac(j,itile))*( soil%remember_wlfrac(j,itile)*zdelta_skin(i,jt,itile)          &
                   +(1._dp-soil%remember_wlfrac(j,itile))*zdelta_soil(i,jt,itile)*soil%remember_rel_hum(j,itile) )
                qwisoair_fact(i,jt,itile) = soil%remember_snfrac(j,itile)*zdelta_snow(i,jt,itile)                          &
                   +(1._dp-soil%remember_snfrac(j,itile))*(soil%remember_wlfrac(j,itile)*zdelta_skin(i,jt,itile)           &
                                                         +(1._dp-soil%remember_wlfrac(j,itile))*zdelta_soil(i,jt,itile))
                qwisosat_fact_fra(i,jt,itile) = 0._dp 
                qwisoair_fact_fra(i,jt,itile) = 0._dp 

               ELSE IF (surface%is_present(j,itile)) THEN
        
                qwisosat_fact(i,jt,itile) = soil%remember_snfrac(j,itile) * zdelta_snow(i,jt,itile)                        &
                   + (1._dp - soil%remember_snfrac(j,itile)) * soil%remember_wlfrac(j,itile) * zdelta_skin(i,jt,itile)
                qwisoair_fact(i,jt,itile) = qwisosat_fact(i,jt,itile)
                qwisosat_fact_fra(i,jt,itile) = 0._dp 
                qwisoair_fact_fra(i,jt,itile) = 0._dp

               END IF
             END IF
             ELSE
             IF (soil%remember_satsphum(j,itile) .NE. 0._dp) THEN
               IF (soil%remember_rel_hum(j,itile) > soil % remember_airm(j,itile) / &
                 soil%remember_satsphum(j,itile) .AND. soil%remember_rel_hum(j,itile) > 1.e-10_dp ) THEN

                qwisosat_fact(i,jt,itile)= soil%remember_snfrac(j,itile)*zdelta_snow(i,jt,itile)                           &
                   +(1._dp-soil%remember_snfrac(j,itile))*( soil%remember_wlfrac(j,itile)*zdelta_skin(i,jt,itile) )
                qwisoair_fact(i,jt,itile) = soil%remember_snfrac(j,itile)*zdelta_snow(i,jt,itile)                          &
                   +(1._dp-soil%remember_snfrac(j,itile))*(soil%remember_wlfrac(j,itile)*zdelta_skin(i,jt,itile))
                qwisosat_fact_fra(i,jt,itile) = wisokinl(i,jt) * (1._dp-soil%remember_snfrac(j,itile))                     &
                                              * (1._dp-soil%remember_wlfrac(j,itile)) * soil%remember_rel_hum(j,itile)
                qwisoair_fact_fra(i,jt,itile) = wisokinl(i,jt) * (1._dp-soil%remember_snfrac(j,itile))                     &
                                              * (1._dp-soil%remember_wlfrac(j,itile))

             ELSE IF (surface%is_present(j,itile)) THEN
       
                qwisosat_fact(i,jt,itile) = soil%remember_snfrac(j,itile) * zdelta_snow(i,jt,itile)                        &
                   + (1._dp - soil%remember_snfrac(j,itile)) * soil%remember_wlfrac(j,itile) * zdelta_skin(i,jt,itile)
                qwisoair_fact(i,jt,itile) = qwisosat_fact(i,jt,itile)
                qwisosat_fact_fra(i,jt,itile) = 0._dp 
                qwisoair_fact_fra(i,jt,itile) = 0._dp

               END IF
             END IF
             END IF

             !dew formation over soil only (not above glacier)
             IF(soil % remember_airm(j,itile)>soil%remember_satsphum(j,itile)) THEN
                IF(surface%is_glacier(j,1)) THEN

                   qwisosat_fact(i,jt,itile)=soil%remember_snfrac(j,itile)*zdelta_snow(i,jt,itile)                          &
                      +(1._dp-soil%remember_snfrac(j,itile))*( soil%remember_wlfrac(j,itile)*zdelta_skin(i,jt,itile)        &
                      +(1._dp-soil%remember_wlfrac(j,itile))*zdelta_soil(i,jt,itile)*soil%remember_rel_hum(i,itile) )
                   qwisoair_fact(i,jt,itile) = soil%remember_snfrac(j,itile)*zdelta_snow(i,jt,itile)                        &
                      +(1._dp-soil%remember_snfrac(j,itile))*(soil%remember_wlfrac(j,itile)*zdelta_skin(i,jt,itile)         &
                                           +(1._dp-soil%remember_wlfrac(j,itile))*zdelta_soil(i,jt,itile))
                   qwisosat_fact_fra(i,jt,itile) = 0._dp 
                   qwisoair_fact_fra(i,jt,itile) = 0._dp 

                ELSE

                   qwisosat_fact(i,jt,itile) = 0._dp 
                   qwisoair_fact(i,jt,itile) = 0._dp          
                   qwisosat_fact_fra(i,jt,itile) = 1._dp 
                   qwisoair_fact_fra(i,jt,itile) = 1._dp

                END IF
             END IF
          END DO
       END DO
    END DO

    ! This loop should be equivalent to the third part of loop 331 of ECHAM5
    ! At this part are csat and cair recalculated with depending at pvgrat
    DO jt = 1, nwiso
       DO itile = 1, ntiles
          DO i = 1, nidx
             j = kidx0 + i - 1
             wiso_csat_tiles(i,jt,itile) = soil%remember_cveg(j,itile) * qwisosat_veg(i,jt,itile)   &
                      + (1._dp - soil%remember_cveg(j,itile)) * qwisosat_fact(i,jt,itile)
             wiso_cair_tiles(i,jt,itile) = soil%remember_cveg(j,itile) * qwisoair_veg(i,jt,itile)   &
                      + (1._dp - soil%remember_cveg(j,itile)) * qwisoair_fact(i,jt,itile)
             wiso_csat_transpiration_tiles(i,jt,itile) = soil%remember_cveg(j,itile) * qwisosat_transpiration(i,jt,itile)
             wiso_csat_tiles_fra(i,jt,itile) = soil%remember_cveg(j,itile) * qwisosat_veg_fra(i,jt,itile) &
                      + (1._dp - soil%remember_cveg(j,itile)) * qwisosat_fact_fra(i,jt,itile)
             wiso_cair_tiles_fra(i,jt,itile) = soil%remember_cveg(j,itile) * qwisoair_veg_fra(i,jt,itile) &
                      + (1._dp - soil%remember_cveg(j,itile)) * qwisoair_fact_fra(i,jt,itile)
             wiso_csat_transp_tiles_fra(i,jt,itile) = soil%remember_cveg(j,itile) * qwisosat_transp_fra(i,jt,itile)
          END DO
       END DO
    END DO

    ! ECHAM5 compatibility: one value for whole grid box
    DO jt=1,nwiso
       CALL average_tiles(wiso_csat_tiles(1:nidx,jt,:), surface % is_present(kidx0:kidx1,:) .AND. .NOT.               &
                      surface % is_lake(kidx0:kidx1,:), soil % remember_cfrac(kidx0:kidx1,:), wiso_csat(1:nidx,jt))
       CALL average_tiles(wiso_cair_tiles(1:nidx,jt,:), surface % is_present(kidx0:kidx1,:) .AND. .NOT.               &
                      surface % is_lake(kidx0:kidx1,:), soil % remember_cfrac(kidx0:kidx1,:), wiso_cair(1:nidx,jt))
       CALL average_tiles(wiso_csat_transpiration_tiles(1:nidx,jt,:), surface % is_present(kidx0:kidx1,:) .AND. .NOT. &
                         surface % is_lake(kidx0:kidx1,:), soil % remember_cfrac(kidx0:kidx1,:), wiso_csat_transpiration(1:nidx,jt))
       CALL average_tiles(wiso_csat_tiles_fra(1:nidx,jt,:), surface % is_present(kidx0:kidx1,:) .AND. .NOT.           &
                      surface % is_lake(kidx0:kidx1,:), soil % remember_cfrac(kidx0:kidx1,:), wiso_csat_fra(1:nidx,jt))
       CALL average_tiles(wiso_cair_tiles_fra(1:nidx,jt,:), surface % is_present(kidx0:kidx1,:) .AND. .NOT.           &
                      surface % is_lake(kidx0:kidx1,:), soil % remember_cfrac(kidx0:kidx1,:), wiso_cair_fra(1:nidx,jt))
       CALL average_tiles(wiso_csat_transp_tiles_fra(1:nidx,jt,:), surface % is_present(kidx0:kidx1,:) .AND. .NOT.    &
                         surface % is_lake(kidx0:kidx1,:), soil % remember_cfrac(kidx0:kidx1,:), wiso_csat_transp_fra(1:nidx,jt))
    END DO

    DO jt = 1, nwiso
       DO i = 1, nidx
          j = kidx0 + i - 1
          soil%wiso_csat(j,jt)            = wiso_csat(i,jt)
          soil%wiso_cair(j,jt)            = wiso_cair(i,jt)
          soil%wiso_csat_transpiration(j,jt) = wiso_csat_transpiration(i,jt)
          soil%wiso_csat_fra(j,jt)        = wiso_csat_fra(i,jt)
          soil%wiso_cair_fra(j,jt)        = wiso_cair_fra(i,jt)
          soil%wiso_csat_transp_fra(j,jt) = wiso_csat_transp_fra(i,jt)
       END DO
    END DO
    
    END IF
!---wiso-code-end

    END IF

    !----------------------------------------------------------------------------------------
    ! Soil moisture in root zone

    DO itile=1,ntiles
       soil_moisture_root(:,itile) = &
            calc_moist_root_zone(soil%moisture(kidx0:kidx1,:,itile), &
                                 soil_param%Depth(kidx0:kidx1,:), root_depth(:,:,itile))
       soil_moisture_root_max(:,itile) = &
            calc_moist_root_zone(soil_param%MaxMoisture(kidx0:kidx1,:), &
                                soil_param%Depth(kidx0:kidx1,:), root_depth(:,:,itile))
    END DO
    soil_moisture_root = MERGE(soil_moisture_root, 0._dp, surface%is_vegetation(kidx0:kidx1,:))
    soil_moisture_root_max = MERGE(soil_moisture_root_max, 0._dp, surface%is_vegetation(kidx0:kidx1,:))

    !------------------------------------------------------------------------------------------

    ! Surface dry static energy 
    !    ... (loop 331 in vdiff)

    DO itile=1,ntiles
       soil%sat_surface_specific_humidity(kidx0:kidx1,itile) = &
            sat_specific_humidity(soil%surface_temperature(kidx0:kidx1,itile), surface_pressure(1:nidx))
       IF (ANY(soil%sat_surface_specific_humidity(kidx0:kidx1,itile) > 0.99_dp*HUGE(1._dp))) THEN
          CALL message('sat_specific_humidity', 'lookup table overflow')
       ENDIF
       sat_surface_specific_hum_old(1:nidx,itile) = soil%sat_surface_specific_humidity(kidx0:kidx1,itile)
    END DO

    CALL average_tiles(soil%surface_temperature(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
                       .AND. .NOT. surface%is_lake(kidx0:kidx1,:), &
                       surface%cover_fract(kidx0:kidx1,:), surface_temperature_avg(:))
    CALL average_tiles(soil%sat_surface_specific_humidity(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. surface%is_lake(kidx0:kidx1,:), &
         surface%cover_fract(kidx0:kidx1,:),sat_surf_specific_hum_avg(:))

    ! Modify minimum canopy resistance (no water stress) according to water limitation in the soil root zone 
    canopy_resistance = 1.e20_dp
    water_stress_factor = 0._dp
    canopy_conductance_limited = 1._dp/1.e20_dp

    DO itile=1,ntiles
    DO i=1,nidx
    j=kidx0+i-1
    IF (surface%is_vegetation(j,itile)) THEN
       water_stress_factor(i,itile) = calc_water_stress_factor(soil_moisture_root(i,itile), soil_moisture_root_max(i,itile), &
                                             soil_options%MoistureFractCritical, soil_options%MoistureFractWilting)
       IF (water_stress_factor(i,itile) > EPSILON(1._dp) .AND. canopy_conductance_max(i,itile) > EPSILON(1._dp) .AND. &
            air_moisture(i) .LE. sat_surf_specific_hum_avg(i)) THEN
          canopy_resistance(i,itile) = 1._dp / (canopy_conductance_max(i,itile) * water_stress_factor(i,itile) + 1.e-20_dp)
       ELSE
          canopy_resistance(i,itile) = 1.e20_dp
       END IF
       canopy_conductance_limited(i,itile) = 1._dp/MAX(canopy_resistance(i,itile),1.e-20_dp)
    END IF
    END DO
    END DO

    ! Sensitivity of saturated surface specific humidity to temperature
    zdqsl = (sat_specific_humidity(surface_temperature_avg(:) + 0.001_dp, surface_pressure(1:nidx)) - &
             sat_surf_specific_hum_avg(:)) * 1000._dp

    ! Compute ECHAM-type skin reservoir and snow fraction (bare soil and soil under canopy)
    skin_reservoir_max(:,:) = MERGE(soil_options%SkinReservoirMax * (1._dp + lai(:,:) *               &
                              SPREAD(surface%veg_ratio_max(kidx0:kidx1),DIM=2,ncopies=ntiles)), 0._dp, &
                              .NOT. surface%is_glacier(kidx0:kidx1,:))

    evaporation_pot = 0._dp
    do itile=1,ntiles
    do i=1,nidx
    j=kidx0+i-1
    IF (skin_reservoir_max(i,itile) > EPSILON(1._dp)) THEN
       canopy_snow_fract(i,itile) = MIN(1._dp, canopy_snow(i,itile) / skin_reservoir_max(i,itile))
    ELSE
       canopy_snow_fract(i,itile) = 0._dp
    END IF

    IF (soil_mask(i,itile)) THEN
       wet_skin_fract(i,itile) = MERGE(MIN(1._dp, soil%skin_reservoir(j,itile) / skin_reservoir_max(i,itile)), &
            0.0_dp, soil_options%SkinReservoirMax > EPSILON(1._dp))
      
       ! Roesch et al, 2002, Climate Dynamics
       soil%snow_fract(j,itile) = zqsncr * TANH(soil%snow(j,itile) * 100._dp) * &
            SQRT(soil%snow(j,itile) * 1000._dp / (soil%snow(j,itile) * 1000._dp + zepsec + &
            zsigfac * surface%oro_std_dev(j)))
       IF (soil%snow_fract(j,itile) .LT. EPSILON(1._dp) .AND. canopy_snow_fract(i,itile) .GE. EPSILON(1._dp)) &
            soil%snow_fract(j,itile) = canopy_snow_fract(i,itile)

       ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than 
       ! equivalent snow water content from soil and canopy; same for skin reservoir
       ! Potential evaporation using old values of air and surface humidity
       evaporation_pot(i,itile) = zcons30 * cdrag(i) * &
            (soil%sat_surface_specific_humidity(j,itile) - air_moisture(i))
       IF (soil%snow_fract(j,itile) > 0._dp) THEN
          soil%snow_fract(j,itile) = soil%snow_fract(j,itile) / &
               MAX(1._dp, soil%snow_fract(j,itile) * evaporation_pot(i,itile) * delta_time / &
               (RhoH2O*(soil%snow(j,itile) + canopy_snow(i,itile))))
       END IF
       IF (wet_skin_fract(i,itile) > 0._dp) THEN
          wet_skin_fract(i,itile) = wet_skin_fract(i,itile) / &
               MAX(1._dp, (1._dp - soil%snow_fract(j,itile)) * evaporation_pot(i,itile) * delta_time / &
               (RhoH2O * MAX(EPSILON(1._dp),soil%skin_reservoir(j,itile))))
       END IF
       
    ELSE IF (surface%is_glacier(j,itile)) THEN
       wet_skin_fract(i,itile) = 0.0_dp
       soil%snow_fract(j,itile) = 1.0_dp
    ELSE
       wet_skin_fract(i,itile) = 0.0_dp
       soil%snow_fract(j,itile) = 0.0_dp
    END IF
    END DO
    END do

    !----------------------------------------------------------------------------------------------------------------------

    ! Note: at the moment, surface_temperature and sat_surface_specific_humidity are the same for all tiles in a grid box
    dry_static_energy(:) = surface_temperature_avg(:) * SpecificHeatDryAirConstPressure * &
         (1._dp+ vtmpc2 * ( csat(:) * sat_surf_specific_hum_avg(:) + (1._dp - cair(:)) * air_moisture(:)))
    
    ! Conversion factor for humidity from dry static energy
    zcpq(:) = dry_static_energy(:) / surface_temperature_avg(:)

    net_radiation(:) = rad_shortwave_net(:) + &
         rad_longwave_down(:) - Emissivity * StefanBoltzmann * surface_temperature_avg(:)**4

    ! Compute new surface temperature and moisture
    dry_static_energy_new = 0._dp
    surface_qsat_new      = 0._dp

    CALL average_tiles(soil%ground_heat_flux(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. surface%is_lake(kidx0:kidx1,:), &
         surface%cover_fract(kidx0:kidx1,:),ground_heat_flux_avg(:))
    CALL average_tiles(soil%heat_capacity(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
         surface%cover_fract(kidx0:kidx1,:),heat_capacity_avg(:))

    IF (.NOT. lstart) THEN
       CALL update_surfacetemp(nidx, zcpq(:),                                         &
            t_Acoef(:), t_Bcoef(:), q_Acoef(:), q_Bcoef(:),                                  &
            dry_static_energy(:), sat_surf_specific_hum_avg(:), zdqsl(:),            &
            net_radiation(:), ground_heat_flux_avg(:),                                       &
            zcons30*cdrag(:), cair(:), csat(:),                                          &
            SUM(surface%cover_fract(kidx0:kidx1,:) * soil%snow_fract(kidx0:kidx1,:), DIM=2), &
            heat_capacity_avg(:),                                                            &
            dry_static_energy_new(:), surface_qsat_new(:))
    ELSE
       dry_static_energy_new(:) = dry_static_energy(:)
       surface_qsat_new(:) = sat_surf_specific_hum_avg(:)
    ENDIF

    ! New land surface temperature and surface saturated humidity
    DO itile=1,ntiles
       IF (.NOT. lstart) THEN
          soil%sat_surface_specific_humidity(kidx0:kidx1,itile) = surface_qsat_new(:)
       END IF
       soil%dry_static_energy_new(kidx0:kidx1,itile) = dry_static_energy_new(:)
    END DO

    ! New unfiltered land surface temperature at t+dt (\tilde(X)^(t+1))
    soil%surface_temperature_unfiltered(kidx0:kidx1,:) = &
         SPREAD((ztpfac2 * dry_static_energy_new(:) + ztpfac3 * dry_static_energy(:)) / zcpq, NCOPIES=ntiles, DIM=2)

    ! Correction for snowmelt
    CALL average_tiles(soil%snow(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
         surface%cover_fract(kidx0:kidx1,:), snow_avg(:))
    WHERE (snow_avg(:) > soil_options%CriticalSnowDepth .OR. ANY(surface%is_glacier(kidx0:kidx1,:), DIM=2)) &
       dry_static_energy_new(:) = (MIN(soil%surface_temperature_unfiltered(kidx0:kidx1,1),tmelt)                &
       * zcpq(:) - ztpfac3 * dry_static_energy(:)) / ztpfac2

!!$    FORALL(itile=1:ntiles)
    DO itile=1,ntiles
       ! Compute temperature and moisture at lowest atmospheric level by back-substitution
       air_dry_static_energy_new(:,itile) = t_Acoef(:) * dry_static_energy_new(:) + t_Bcoef(:)
       air_moisture_new(:,itile) = q_Acoef(:) * soil%sat_surface_specific_humidity(kidx0:kidx1,itile) + q_Bcoef(:)
    END DO
!!$    END FORALL

!---wiso-code
    IF (lwiso) THEN

    DO jt=1,nwiso
       DO itile = 1,ntiles  
          zwisofracl(1:nidx,jt,itile) = exp( talphal1(jt) / ( soil%surface_temperature(kidx0:kidx1,itile)**2._dp )          &
                                      + talphal2(jt) / ( soil%surface_temperature(kidx0:kidx1,itile) ) + talphal3(jt) )
          zwisofraci(1:nidx,jt,itile) = exp( talphas1(jt) / ( soil%surface_temperature(kidx0:kidx1,itile)**2._dp )          &
                                      + talphas2(jt) / ( soil%surface_temperature(kidx0:kidx1,itile) ) + talphas3(jt) )
          ! effective fractionation over ice if necessary (isotopes only)
          WHERE (nwisotyp(jt).NE.1.AND.soil%surface_temperature(kidx0:kidx1,itile).LT.tmelt)
             zsatval(1:nidx,itile) = tsatbase - tsatfac * (soil%surface_temperature(kidx0:kidx1,itile) - tmelt)
             zwisofraci(1:nidx,jt,itile) = zwisofraci(1:nidx,jt,itile) * ( zsatval(1:nidx,itile) &
                                         / (1._dp + zwisofraci(1:nidx,jt,itile) * (zsatval(1:nidx,itile)-1._dp)*tdifrel(jt) ) )
          END WHERE
       END DO
    END DO

    DO jt=1,nwiso
       DO itile = 1,ntiles
          DO i=1,nidx
             j = kidx0 + i - 1
             IF (soil % remember_airm(j,itile)>soil%remember_satsphum(j,itile).AND.soil_mask(i,itile)) THEN
                ! if new surface air moisture is greater then vapour
                IF (soil%sat_surface_specific_humidity(j,itile).GT.air_moisture(i)) THEN
                   ! all vapour condensates and no fractionation occurs
                   soil%wiso_sat_surface_specific_hum(j,jt,itile) = &
                        zdelta_air_tiles(i,jt,itile) * soil%sat_surface_specific_humidity(j,itile) 
                ELSE
                   IF (soil%dry_static_energy_new(j,itile).GT.tmelt) THEN
                   ! dew formation - fractionation as a closed system
                      zwisodenom = air_moisture(i)+(zwisofraci(i,jt,itile)-1._dp)*soil%sat_surface_specific_humidity(j,itile)
                      zdelta=tnat(jt)
                      IF (ABS(zwisodenom).GT.cwisomin) zdelta=zwisofracl(i,jt,itile)*wiso_air_moisture(i,jt)/zwisodenom
                      IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                      soil%wiso_sat_surface_specific_hum(j,jt,itile) = &
                                    zdelta * soil%sat_surface_specific_humidity(j,itile)
                   ELSE
                   ! rime formation - fractionation as a open system
                      IF (ABS(air_moisture(i)).GT.cwisomin) THEN
                         zdelta=1._dp-(soil%sat_surface_specific_humidity(j,itile)/air_moisture(i))
                         soil%wiso_sat_surface_specific_hum(j,jt,itile) = &
                                    wiso_air_moisture(i,jt) * (1._dp - zdelta**zwisofraci(i,jt,itile))
                      ELSE
                         zdelta=tnat(jt)
                         soil%wiso_sat_surface_specific_hum(j,jt,itile) = &
                                    zdelta * soil%sat_surface_specific_humidity(j,itile)
                      END IF
                   END IF
                END IF
             ELSE
                zwisofracl(i,jt,itile) = zdelta_soil(i,jt,itile) / zwisofracl(i,jt,itile)
                soil%wiso_sat_surface_specific_hum(j,jt,itile) = &
                     zwisofracl(i,jt,itile) * soil%sat_surface_specific_humidity(j,itile)
             END IF
          END DO
       END DO
    END DO

!   include wiso_cair_frac and kin in Richtmeyer-Morton 
    DO jt=1,nwiso
       q_wiso_Acoef_new(1:nidx,jt) = 1._dp / ( wiso_nenner(1:nidx,jt) + &
                                              wiso_cair_fra(1:nidx,jt) * q_wiso_Acoef(1:nidx,jt) ) &
                                   * q_wiso_Acoef(1:nidx,jt)
       q_wiso_Bcoef_new(1:nidx,jt) = 1._dp / ( wiso_nenner(1:nidx,jt) + &
                                              wiso_cair_fra(1:nidx,jt) * q_wiso_Acoef(1:nidx,jt) ) &
                                   * q_wiso_Bcoef(1:nidx,jt)
    END DO

    ! Compute moisture at lowest atmospheric level by back-substitution - water isotope
    ! Calculation is analog at zwisoqklevl at update land
    DO jt=1,nwiso
       DO itile = 1,ntiles
          wiso_air_moisture_new(1:nidx,jt,itile) = ( q_wiso_Acoef_new(1:nidx,jt) * &
                 ( wiso_csat(1:nidx,jt) * soil%sat_surface_specific_humidity(kidx0:kidx1,itile) &
                   - wiso_cair(1:nidx,jt) * (air_moisture_new(1:nidx,itile) * wiso_helpqdif(1:nidx)) * cvdifts ) ) &
                 + ( q_wiso_Acoef_new(1:nidx,jt) * wiso_csat_fra(1:nidx,jt) &
                                         * soil%wiso_sat_surface_specific_hum(kidx0:kidx1,jt,itile) ) & 
              + q_wiso_Bcoef_new(1:nidx,jt)
       END DO
    END DO
    
    END IF
!---wiso-code-end

    DO itile=1,ntiles
       soil%radiative_temperature(kidx0:kidx1,itile) = dry_static_energy_new(:) / zcpq(:)
    END DO

    ! Compute sensible heat flux
!!$    FORALL(itile=1:ntiles)
    DO itile=1,ntiles
       soil%sensible_heat_flux(kidx0:kidx1,itile) = zcons30 * cdrag(:)                                               &
       * (air_dry_static_energy_new(:,itile) - dry_static_energy_new(:) / zcpq(:) * SpecificHeatDryAirConstPressure  &
       * (1._dp + vtmpc2 * (air_moisture(:) + cair(:) * (air_moisture_new(:,itile) - air_moisture(:)))))
    END DO
!!$    END FORALL

    ! Compute evaporation and latent heat flux
!!$    FORALL(itile=1:ntiles)
    DO itile=1,ntiles
       ! Evapotranspiration
       soil%evapotranspiration(kidx0:kidx1,itile) = zcons30 * cdrag(:)  &
       * (cair(:) * air_moisture_new(:,itile) - csat(:)  * soil%sat_surface_specific_humidity(kidx0:kidx1,itile))
       ! Transpiration
       soil%transpiration(kidx0:kidx1,itile) = zcons30 * cdrag(:)  &
       * csat_transpiration(:) * (air_moisture_new(:,itile) - soil%sat_surface_specific_humidity(kidx0:kidx1,itile))
       ! Potential evaporation
       soil%evaporation_pot(kidx0:kidx1,itile) = zcons30 * cdrag(:)  &
       * (air_moisture_new(:,itile) - soil%sat_surface_specific_humidity(kidx0:kidx1,itile))
    END DO
!!$    END FORALL

!---wiso-code
    IF (lwiso) THEN

    ! Compute evaporation - water isotopes
    DO jt=1,nwiso
       DO itile = 1,ntiles
          DO i=1,nidx
             j=kidx0+i-1
             ! Evapotranspiration
             soil % wiso_evapotranspiration(j,jt,itile) = zcons30 * cdrag(i) * ( &
                          ( wiso_cair(i,jt) * air_moisture_new(i,itile) &
                            - wiso_csat(i,jt) * soil%sat_surface_specific_humidity(j,itile) ) &
                        + ( wiso_cair_fra(i,jt) * wiso_air_moisture_new(i,jt,itile) &
                            - wiso_csat_fra(i,jt) * soil % wiso_sat_surface_specific_hum(j,jt,itile) ) )
             ! Transpiration
             soil % wiso_transpiration(j,jt,itile) = zcons30 * cdrag(i) * ( &
                          ( wiso_csat_transpiration(i,jt) * ( air_moisture_new(i,itile) &
                                                            - soil%sat_surface_specific_humidity(j,itile) ) ) &
                        + ( wiso_csat_transp_fra(i,jt) * ( wiso_air_moisture_new(i,jt,itile) &
                                                           - soil % wiso_sat_surface_specific_hum(j,jt,itile) ) ) )
             ! Potential evaporation
             soil % wiso_evaporation_pot(j,jt,itile) = soil%evaporation_pot(j,itile)
          END DO
       END DO
    END DO
    
    END IF
!---wiso-code-end

    ! Latent heat flux
    soil%latent_heat_flux(kidx0:kidx1,:) = LatentHeatVaporization  * soil%evapotranspiration(kidx0:kidx1,:) &
    + (LatentHeatSublimation - LatentHeatVaporization) * soil%snow_fract(kidx0:kidx1,:) * soil%evaporation_pot(kidx0:kidx1,:)

    !
    !--------------------------------------------------------------------------------------------------------
    ! 
    tte_corr = 0._dp
    glacier_precip_minus_evap = 0._dp
    melt_water_excess = 0._dp
    surface_runoff = 0._dp
    drainage = 0._dp

!---wiso-code
    IF (lwiso) THEN

    wiso_glac_p_minus_e = 0._dp
    wiso_melt_water_excess = 0._dp
    wiso_surface_runoff = 0._dp
    wiso_drainage = 0._dp
    
    END IF
!---wiso-code-end

    IF (lwiso) THEN
    
    DO itile=1,ntiles
       CALL update_surf_down(&
                ! kbdim
                nidx, &
                ! pqm1
                air_moisture(1:nidx), &
                ! ptsl, pspeed10, ptm1
                soil % surface_temperature_unfiltered(kidx0:kidx1, itile), wind10(1:nidx), air_temperature(1:nidx), &
                ! pwl, pcvw, pwlmx
                soil % skin_reservoir(kidx0:kidx1, itile), wet_skin_fract(1:nidx, itile), skin_reservoir_max(1:nidx, itile), &
                ! psn, pcvs
                soil % snow(kidx0:kidx1, itile), soil % snow_fract(kidx0:kidx1, itile), &
                ! psnc, pgld
                canopy_snow(1:nidx, itile), soil % glacier_depth(kidx0:kidx1, itile), &
                ! pgrndcapc
                soil % heat_capacity(kidx0:kidx1, itile), &
                ! pevapl, pevapot
                soil % evapotranspiration(kidx0:kidx1, itile), soil % evaporation_pot(kidx0:kidx1, itile), &
                ! pevwsd
                evapotranspiration_no_snow_skin(1:nidx), &
                ! prain, psnow
                rain, snow, &
                ! lmask
                surface % is_present(kidx0:kidx1, itile), &
                ! lpglac
                surface % is_glacier(kidx0:kidx1, itile), &
                ! palac
                glacier_precip_minus_evap(:, itile), &
                ! psnacl, psnmel
                soil % snow_acc(kidx0:kidx1, itile), soil % snow_melt_acc(kidx0:kidx1, itile), &
                ! progl, pmlres
                soil % glacier_runoff_acc(kidx0:kidx1, itile), melt_water_excess(1:nidx), &
                tte_corr(1:nidx, itile), &
!---wiso-code
                lwisofracl, &
                ! pwisosn
                soil % wiso_snow(kidx0:kidx1,1:nwiso,itile), &
                ! pwisosnc, pwisogld
                wiso_canopy_snow(1:nidx,1:nwiso,itile), soil % wiso_glacier_depth(kidx0:kidx1,1:nwiso,itile), &
                ! psnglac, pwisosnglac
                soil % snowglac(kidx0:kidx1, itile), soil % wiso_snowglac(kidx0:kidx1,1:nwiso,itile), &
                ! pwisosnmel
                soil % wiso_snow_melt_acc(kidx0:kidx1,1:nwiso,itile), &
                ! pwisowl
                soil % wiso_skin_reservoir(kidx0:kidx1,1:nwiso,itile), &
                ! pwisoevapl, 
                soil % wiso_evapotranspiration(kidx0:kidx1,1:nwiso,itile), &
                ! pwisoevapot
                soil % wiso_evaporation_pot(kidx0:kidx1,1:nwiso,itile), &
                ! pwisoevwsd
                wiso_evapotrans_no_snow_skin(1:nidx,1:nwiso), &
                ! pwisorain, pwisosnow, pwisorogl
                wiso_rain, wiso_snow, soil % wiso_glacier_runoff_acc(kidx0:kidx1,1:nwiso,itile), &
                ! pwisoalac, pwisosnacl
                wiso_glac_p_minus_e(:,1:nwiso,itile), soil % wiso_snow_acc(kidx0:kidx1,1:nwiso,itile), &
                ! pwisomlres
                wiso_melt_water_excess(1:nidx,1:nwiso))
!---wiso-code-end

       CALL update_soiltemp(nidx,soil%ntsoil                                                               &
                         , soil%surface_temperature_unfiltered(kidx0:kidx1,itile),soil%snow(kidx0:kidx1,itile) &
                         , soil_param%ThermalDiffusivity(kidx0:kidx1)                                      &
                         , soil_param%VolHeatCapacity(kidx0:kidx1)                                         &
                         , soil%c_soil_temperature(kidx0:kidx1,:,itile)                                    &
                         , soil%d_soil_temperature(kidx0:kidx1,:,itile)                                    &
                         , soil%soil_temperature(kidx0:kidx1,:,itile)                                      &
                         , soil%heat_capacity(kidx0:kidx1,itile),soil%ground_heat_flux(kidx0:kidx1,itile)  &
                         , surface%is_present(kidx0:kidx1,itile)                                           &
                         , surface%is_glacier(kidx0:kidx1,itile))
   
       CALL update_surf_up(&
                nidx, itile, &
                ! pwl, pcvw,  pwlmx
                soil % skin_reservoir(kidx0:kidx1, itile), wet_skin_fract(:, itile), skin_reservoir_max(:, itile), &
                ! pws, pwsmx
                soil % moisture(kidx0:kidx1, 1, itile), soil_param % MaxMoisture(kidx0:kidx1, 1), &
                ! ptsoil, pcvs
                soil % soil_temperature(kidx0:kidx1, 1, itile), soil % snow_fract(kidx0:kidx1, itile), &
                ! porostd
                surface % oro_std_dev(kidx0:kidx1), &
                ! pevapl, pevapot
                soil % evapotranspiration(kidx0:kidx1, itile), soil % evaporation_pot(kidx0:kidx1, itile), &
                ! pevwsd
                evapotranspiration_no_snow_skin(1:nidx), &
                ! prain
                rain(1:nidx), &
                ! lmask, lpglac
                surface % is_present(kidx0:kidx1, itile), surface % is_glacier(kidx0:kidx1, itile), &
                ! prunoff, pdrain
                soil % runoff_acc(kidx0:kidx1, itile), soil % drainage_acc(kidx0:kidx1, itile), &
                ! pros_hd, pdrain_hd, pmlres
                surface_runoff(1:nidx, itile), drainage(1:nidx, itile), melt_water_excess(1:nidx), &
!---wiso-code
                lwisofracl, &
                ! pwisows, pwisowl
                soil%wiso_moisture(kidx0:kidx1,1:nwiso,itile), soil%wiso_skin_reservoir(kidx0:kidx1,1:nwiso,itile), &
                ! pwisorain
                wiso_rain(1:nidx,1:nwiso), &
                ! pwisorunoff,     pwisodrain
                soil%wiso_runoff_acc(kidx0:kidx1,1:nwiso,itile), soil%wiso_drainage_acc(kidx0:kidx1,1:nwiso,itile), &
                ! pwisoevapot,     pwisoevwsd
                soil % wiso_evaporation_pot(kidx0:kidx1,1:nwiso,itile), wiso_evapotrans_no_snow_skin(1:nidx,1:nwiso), &
                ! pwisoros_hd,     pwisodrain_hd
                wiso_surface_runoff(1:nidx,1:nwiso,itile), wiso_drainage(1:nidx,1:nwiso,itile), &
                ! pwisomlres
                wiso_melt_water_excess(1:nidx,1:nwiso))
!---wiso-code-end
    END DO
    
    ELSE
    
    DO itile=1,ntiles
       CALL update_surf_down(&
                ! kbdim
                nidx, &
                ! pqm1
                air_moisture(1:nidx), &
                ! ptsl, pspeed10, ptm1
                soil % surface_temperature_unfiltered(kidx0:kidx1, itile), wind10(1:nidx), air_temperature(1:nidx), &
                ! pwl, pcvw, pwlmx
                soil % skin_reservoir(kidx0:kidx1, itile), wet_skin_fract(1:nidx, itile), skin_reservoir_max(1:nidx, itile), &
                ! psn, pcvs
                soil % snow(kidx0:kidx1, itile), soil % snow_fract(kidx0:kidx1, itile), &
                ! psnc, pgld
                canopy_snow(1:nidx, itile), soil % glacier_depth(kidx0:kidx1, itile), &
                ! pgrndcapc
                soil % heat_capacity(kidx0:kidx1, itile), &
                ! pevapl, pevapot
                soil % evapotranspiration(kidx0:kidx1, itile), soil % evaporation_pot(kidx0:kidx1, itile), &
                ! pevwsd
                evapotranspiration_no_snow_skin(1:nidx), &
                ! prain, psnow
                rain, snow, &
                ! lmask
                surface % is_present(kidx0:kidx1, itile), &
                ! lpglac
                surface % is_glacier(kidx0:kidx1, itile), &
                ! palac
                glacier_precip_minus_evap(:, itile), &
                ! psnacl, psnmel
                soil % snow_acc(kidx0:kidx1, itile), soil % snow_melt_acc(kidx0:kidx1, itile), &
                ! progl, pmlres
                soil % glacier_runoff_acc(kidx0:kidx1, itile), melt_water_excess(1:nidx), &
                tte_corr(1:nidx, itile))

       CALL update_soiltemp(nidx,soil%ntsoil                                                               &
                         , soil%surface_temperature_unfiltered(kidx0:kidx1,itile),soil%snow(kidx0:kidx1,itile) &
                         , soil_param%ThermalDiffusivity(kidx0:kidx1)                                      &
                         , soil_param%VolHeatCapacity(kidx0:kidx1)                                         &
                         , soil%c_soil_temperature(kidx0:kidx1,:,itile)                                    &
                         , soil%d_soil_temperature(kidx0:kidx1,:,itile)                                    &
                         , soil%soil_temperature(kidx0:kidx1,:,itile)                                      &
                         , soil%heat_capacity(kidx0:kidx1,itile),soil%ground_heat_flux(kidx0:kidx1,itile)  &
                         , surface%is_present(kidx0:kidx1,itile)                                           &
                         , surface%is_glacier(kidx0:kidx1,itile))
   
       CALL update_surf_up(&
                nidx, itile, &
                ! pwl, pcvw,  pwlmx
                soil % skin_reservoir(kidx0:kidx1, itile), wet_skin_fract(:, itile), skin_reservoir_max(:, itile), &
                ! pws, pwsmx
                soil % moisture(kidx0:kidx1, 1, itile), soil_param % MaxMoisture(kidx0:kidx1, 1), &
                ! ptsoil, pcvs
                soil % soil_temperature(kidx0:kidx1, 1, itile), soil % snow_fract(kidx0:kidx1, itile), &
                ! porostd
                surface % oro_std_dev(kidx0:kidx1), &
                ! pevapl, pevapot
                soil % evapotranspiration(kidx0:kidx1, itile), soil % evaporation_pot(kidx0:kidx1, itile), &
                ! pevwsd
                evapotranspiration_no_snow_skin(1:nidx), &
                ! prain
                rain(1:nidx), &
                ! lmask, lpglac
                surface % is_present(kidx0:kidx1, itile), surface % is_glacier(kidx0:kidx1, itile), &
                ! prunoff, pdrain
                soil % runoff_acc(kidx0:kidx1, itile), soil % drainage_acc(kidx0:kidx1, itile), &
                ! pros_hd, pdrain_hd, pmlres
                surface_runoff(1:nidx, itile), drainage(1:nidx, itile), melt_water_excess(1:nidx))
    END DO
    
    END IF ! lwiso


!---wiso-code
    IF (lwiso) THEN

    ! For a restart without previous isotope diagnostics: 
    ! Initialize *wiso* variables from default JSBACH restart fields
    ! (assume at start a delta value of zero permill) 
    
    IF (lresume .AND. (.NOT. lwiso_rerun)) THEN
      DO jt=1,nwiso
        DO itile=1,ntiles
          soil % snowglac(kidx0:kidx1, itile) = snglacmx
          soil%wiso_snowglac(kidx0:kidx1,jt,itile) = tnat(jt)*soil%snowglac(kidx0:kidx1,itile)
          soil%wiso_evapotranspiration(kidx0:kidx1,jt,itile) = tnat(jt)*soil%evapotranspiration(kidx0:kidx1,itile)
          soil%wiso_evapotranspiration_acc(kidx0:kidx1,jt,itile) = tnat(jt)*soil%evapotranspiration_acc(kidx0:kidx1,itile)
          soil%wiso_evaporation_pot(kidx0:kidx1,jt,itile) = tnat(jt)*soil%evaporation_pot(kidx0:kidx1,itile)
          soil%wiso_transpiration(kidx0:kidx1,jt,itile) = tnat(jt)*soil%transpiration(kidx0:kidx1,itile)
          soil%wiso_transpiration_acc(kidx0:kidx1,jt,itile) = tnat(jt)*soil%transpiration_acc(kidx0:kidx1,itile)
          soil%wiso_skin_reservoir(kidx0:kidx1,jt,itile) = tnat(jt)*soil%skin_reservoir(kidx0:kidx1,itile)
          soil%wiso_moisture(kidx0:kidx1,jt,itile) = tnat(jt)*soil%moisture(kidx0:kidx1,1,itile)
          soil%wiso_snow(kidx0:kidx1,jt,itile) = tnat(jt)*soil%snow(kidx0:kidx1,itile)
          soil%wiso_snow_acc(kidx0:kidx1,jt,itile) = tnat(jt)*soil%snow_acc(kidx0:kidx1,itile)
          soil%wiso_snow_melt_acc(kidx0:kidx1,jt,itile) = tnat(jt)*soil%snow_melt_acc(kidx0:kidx1,itile)
          soil%wiso_runoff_acc(kidx0:kidx1,jt,itile) = tnat(jt)*soil%runoff_acc(kidx0:kidx1,itile)
          soil%wiso_drainage_acc(kidx0:kidx1,jt,itile) = tnat(jt)*soil%drainage_acc(kidx0:kidx1,itile)
          soil%wiso_glacier_depth(kidx0:kidx1,jt,itile) = tnat(jt)*soil%glacier_depth(kidx0:kidx1,itile)
          soil%wiso_glac_p_minus_e_acc(kidx0:kidx1,jt,itile) = tnat(jt)*soil%glacier_precip_minus_evap_acc(kidx0:kidx1,itile)
          soil%wiso_glacier_runoff_acc(kidx0:kidx1,jt,itile) = tnat(jt)*soil%glacier_runoff_acc(kidx0:kidx1,itile)

          wiso_canopy_snow(:,jt,itile) = tnat(jt) * canopy_snow(:,itile)
          wiso_glac_p_minus_e(:,jt,itile) = tnat(jt) * glacier_precip_minus_evap(:,itile)
          wiso_surface_runoff(:,jt,itile) = tnat(jt) * surface_runoff(:,itile)
          wiso_drainage(:,jt,itile) = tnat(jt) * drainage(:,itile)
        END DO    
      END DO
    END IF

    END IF
!---wiso-code-end

    CALL average_tiles(tte_corr, surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
         surface%cover_fract(kidx0:kidx1,:), tte_corr_avg)
    CALL average_tiles(soil%moisture(kidx0:kidx1,1,:), surface%is_present(kidx0:kidx1,:) &
                      .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                      soil_moisture_avg)
    CALL average_tiles(soil%snow(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:) &
                      .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                      snow_avg)
    CALL average_tiles(soil%surface_temperature_unfiltered(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:) &
                      .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
                      surface_temperature_avg)

!---wiso-code
    IF (lwiso) THEN

     !--- calculate average over the tiles at the gridboxes at soil%moisture and soil%snow
     DO jt=1,nwiso
        CALL average_tiles(soil % wiso_moisture(kidx0:kidx1,jt,:), surface % is_present(kidx0:kidx1,:)          &
                           .AND. .NOT.surface % is_lake(kidx0:kidx1,:), surface % cover_fract(kidx0:kidx1,:),   &
                           wiso_soil_moisture_avg(:,jt))
        CALL average_tiles(soil % wiso_snow(kidx0:kidx1,jt,:), surface % is_present(kidx0:kidx1,:)              &
                           .AND. .NOT.surface % is_lake(kidx0:kidx1,:), surface % cover_fract(kidx0:kidx1,:),   &
                           wiso_snow_avg(:,jt))
     END DO
     !--- calculate average over the tiles at the gridboxes at soil%snowglac
     CALL average_tiles(soil%snowglac(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:)                         &
                      .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:),           &
                      snowglac_avg)
     DO jt=1,nwiso
        CALL average_tiles(soil % wiso_snowglac(kidx0:kidx1,jt,:), surface % is_present(kidx0:kidx1,:)          &
                           .AND. .NOT.surface % is_lake(kidx0:kidx1,:), surface % cover_fract(kidx0:kidx1,:),   &
                           wiso_snowglac_avg(:,jt))
     END DO
    
    END IF
!---wiso-code-end

    ! Test that soil moisture has still a positive value
    IF (MINVAL(soil_moisture_avg(:)) < 0._dp .AND. flag_soil_moisture) THEN
      CALL message('update_soil','Warning: Negative soil moisture occurred!')
      flag_soil_moisture = .false.
    END IF
    soil_moisture_avg(:) = MAX(EPSILON(1._dp),soil_moisture_avg(:))
    DO itile=1,ntiles
      soil%moisture(kidx0:kidx1,1,itile) = soil_moisture_avg(:)
      soil%snow(kidx0:kidx1,itile) = snow_avg(:)
      soil%surface_temperature_unfiltered(kidx0:kidx1,itile) = surface_temperature_avg(:)
      soil%soil_temperature(kidx0:kidx1,1,itile) = surface_temperature_avg(:)
    END DO

!---wiso-code
    IF (lwiso) THEN

    !--- act. values
    DO jt=1,nwiso
       wiso_soil_moisture_avg(:,jt) = MAX(EPSILON(1._dp),wiso_soil_moisture_avg(:,jt))
       DO itile = 1, ntiles
          DO i=1,nidx
             j=kidx0+i-1
             soil % wiso_moisture(j, jt, itile) = wiso_soil_moisture_avg(i,jt)
             soil % wiso_snow(j, jt, itile) = wiso_snow_avg(i,jt)
          END DO
       END DO
    END DO

    DO itile=1,ntiles
      soil%snowglac(kidx0:kidx1,itile) = snowglac_avg(:)
    END DO
    DO jt=1,nwiso
       DO itile = 1, ntiles
          DO i=1,nidx
             j=kidx0+i-1
             soil % wiso_snowglac(j, jt, itile) = wiso_snowglac_avg(i,jt)
          END DO
       END DO
    END DO

    END IF
!---wiso-code-end

        ! Time filter for surface temperature
    IF (lsurf) THEN
       IF (.NOT. lstart) THEN
!          WHERE (surface%is_present(kidx0:kidx1,:))
             soil%surface_temperature(kidx0:kidx1,:) = soil%surface_temperature_old(kidx0:kidx1,:) + &
                  eps * (soil%surface_temperature(kidx0:kidx1,:)                                     &
                  - 2._dp * soil%surface_temperature_old(kidx0:kidx1,:)                                 &
                  + soil%surface_temperature_unfiltered(kidx0:kidx1,:))             
             soil%surface_temperature_old(kidx0:kidx1,:) = soil%surface_temperature_unfiltered(kidx0:kidx1,:)
!          END WHERE
       ELSE
          soil%surface_temperature(kidx0:kidx1,:) = soil%surface_temperature_old(kidx0:kidx1,:)
       END IF
    ELSE
       soil%surface_temperature_old(kidx0:kidx1,:) = soil%surface_temperature(kidx0:kidx1,:)
       soil%surface_temperature_unfiltered(kidx0:kidx1,:) = soil%surface_temperature(kidx0:kidx1,:)
    END IF

    !--------------------------------------------------------------------------------------------------------
    ! Re-calculate qsat,qair and zhsoil with updated moisture parameters due to time mismatch with old echam
    !--------------------------------------------------------------------------------------------------------

    DO itile=1,ntiles
    DO i=1,nidx
    j=kidx0+i-1
    IF (soil_mask(i,itile)) THEN
       wet_skin_fract(i,itile) = MERGE(MIN(1._dp, soil%skin_reservoir(j,itile) / skin_reservoir_max(i,itile)), &
            0.0_dp, soil_options%SkinReservoirMax > EPSILON(1._dp))
      
       ! Roesch et al, 2002, Climate Dynamics
       soil%snow_fract(j,itile) = zqsncr * TANH(soil%snow(j,itile) * 100._dp) * &
            SQRT(soil%snow(j,itile)*1000. / (soil%snow(j,itile)*1000._dp + zepsec + &
            zsigfac * surface%oro_std_dev(j)))

       IF (soil%snow_fract(j,itile) .LT. EPSILON(1._dp) .AND. canopy_snow_fract(i,itile) .GE. EPSILON(1._dp)) &
            soil%snow_fract(j,itile) = canopy_snow_fract(i,itile)
       ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than
       ! equivalent snow water content from soil and canopy; same for skin reservoir
       ! Potential evaporation using old values of air and surface humidity
       evaporation_pot(i,itile) = zcons30 * cdrag(i) * &
            (sat_surface_specific_hum_old(i,itile) - air_moisture(i))
       IF (soil%snow_fract(j,itile) > 0._dp .AND. soil%snow(j,itile) + canopy_snow(i,itile) > 0._dp) THEN
          soil%snow_fract(j,itile) = soil%snow_fract(j,itile) / &
               MAX(1._dp, soil%snow_fract(j,itile) * evaporation_pot(i,itile) * delta_time / &
               (RhoH2O*(soil%snow(j,itile) + canopy_snow(i,itile))))
       END IF

       IF (wet_skin_fract(i,itile) > 0._dp) THEN
          wet_skin_fract(i,itile) = wet_skin_fract(i,itile) / &
               MAX(1._dp, (1._dp - soil%snow_fract(j,itile)) * evaporation_pot(i,itile) * delta_time / &
               (RhoH2O * MAX(EPSILON(1._dp),soil%skin_reservoir(j,itile))))
       END IF
       
    ELSE IF (surface%is_glacier(j,itile)) THEN
       wet_skin_fract(i,itile) = 0.0_dp
       soil%snow_fract(j,itile) = 1.0_dp
    ELSE
       wet_skin_fract(i,itile) = 0.0_dp
       soil%snow_fract(j,itile) = 0.0_dp
    END IF
    END DO
    END DO

    relative_humidity = 0.0_dp

    ! Calculate relative humidity from water content in first soil layer
    relative_humidity = &
         calc_relative_humidity(soil%moisture(kidx0:kidx1,1,:), &
         SPREAD(soil_param%MaxMoisture(kidx0:kidx1,1), NCOPIES=ntiles, DIM=2), &
         soil_options%MoistureFractWilting)

    !------------------------------------------------------------------------------------------------------

    qsat_fact = 0._dp
    qair_fact = 0._dp
    qsat_veg(:,:) = 0._dp
    qair_veg(:,:) = 0._dp
    qsat_transpiration(:,:) = 0._dp
    zhsoil = 0._dp

    ! Entspricht der Berechnung von 'zwet' in Loop 331 im ECHAM5 Mode
    DO itile=1,ntiles
    DO i=1,nidx
    j=kidx0+i-1  
    IF (surface%is_present(j,itile)) THEN

    ! Neuberechnung wenn es sich um Land mit Bedingung Feuchte > Feuchtekapazitaet in allen (Blatt-)Schichten
    ! qsat_veg = schnebedeckter Teil + wasserbedeckter Teil bzw. 'skin-reservoir' + Teil der neu durchs Blattdach kommt
    ! qsat_transpiration = max. was neu an zutranspiriert und liegen bleibt
 
       IF (soil%moisture(j,1,itile) > & 
            (soil_options%MoistureFractWilting * soil_param%MaxMoisture(j,1))) THEN
          qsat_veg(i,itile) = soil%snow_fract(j,itile) + (1._dp - soil%snow_fract(j,itile)) *   &
                              (wet_skin_fract(i,itile) + (1._dp - wet_skin_fract(i,itile)) /    &
                              (1._dp + p_echam_zchl(i) * canopy_resistance(i,itile) *           &
                              MAX(1.0_dp,windspeed(i))))
          qsat_transpiration(i,itile) = (1._dp - soil%snow_fract(j,itile)) *                    &
                                        (1._dp - wet_skin_fract(i,itile)) /                     &
                                        (1._dp + p_echam_zchl(i) * canopy_resistance(i,itile) * &
                                        MAX(1.0_dp,windspeed(i)))
       ELSE

    ! Neuberechnung wenn es sich um Land mit Bedingung Feuchte <= Feuchtekapazitaet in allen (Blatt-)Schichten
    ! qsat_veg= schnebedeckter Teil + wasserbedeckter Teil bzw. 'skin-reservoir'

          qsat_veg(i,itile) = soil%snow_fract(j,itile) + (1._dp - soil%snow_fract(j,itile)) *   &
                              wet_skin_fract(i,itile)
          qsat_transpiration(i,itile) = 0._dp
       END IF
       qair_veg(i,itile) = qsat_veg(i,itile)
    END IF
    END DO
    END DO
 
    ! Falls 'Saettigungsverhaeltnis (Luft mit Wasser)' > 'Feuchte / Saettigungsfeuchte'
    ! qair_fact = 1
    ! qsat_fact = schnebedeckter Teil + 'skin-reservoir' + nicht wasserbedeckter Teil * 'Saettigungsverhaeltnis'

    WHERE (relative_humidity(1:nidx,:) > SPREAD(air_moisture(1:nidx), NCOPIES=ntiles, DIM=2) / &
         SPREAD(sat_surf_specific_hum_avg(1:nidx), NCOPIES=ntiles, DIM=2) .AND. &
         relative_humidity(1:nidx,:) > 1.e-10_dp )

       qsat_fact = soil%snow_fract(kidx0:kidx1,:) +                                   & 
                   (1._dp - soil%snow_fract(kidx0:kidx1,:)) *                         &
                   (wet_skin_fract(1:nidx,:) + (1._dp - wet_skin_fract(1:nidx,:)) *   &
                   relative_humidity(1:nidx,:))
       qair_fact = 1._dp

    ! Falls 'Saettigungsverhaeltnis (Luft mit Wasser)' <= 'Feuchte / Saettigungsfeuchte'
    ! qair_fact = qsat_fact
    ! qsat_fact = schnebedeckter Teil + 'skin-reservoir' 

    ELSEWHERE (surface%is_present(kidx0:kidx1,:))

       qsat_fact = soil%snow_fract(kidx0:kidx1,:) + &
                   (1._dp - soil%snow_fract(kidx0:kidx1,:)) * wet_skin_fract(1:nidx,:)
       qair_fact = qsat_fact

    END WHERE

    ! Falls Feuchte > 'Saettigungsfeuchte', ist beides 1!

    WHERE(SPREAD(air_moisture(1:nidx), NCOPIES=ntiles, DIM=2) > &
         SPREAD(sat_surf_specific_hum_avg(1:nidx), NCOPIES=ntiles, DIM=2))
       qsat_fact = 1._dp
       qair_fact = 1._dp
    END WHERE
  
    ! An Stellem mit Vegetation wird q..._veg verwendet, sonst q..._fract

    csat_tiles(1:nidx,:)=Surface%veg_ratio_actual_per_tile(kidx0:kidx1,:) * qsat_veg(1:nidx,:) +            &
                         (1._dp - Surface%veg_ratio_actual_per_tile(kidx0:kidx1,:)) * qsat_fact(1:nidx,:)
    cair_tiles(1:nidx,:)=Surface%veg_ratio_actual_per_tile(kidx0:kidx1,:) * qair_veg(1:nidx,:) +            &
                         (1._dp - Surface%veg_ratio_actual_per_tile(kidx0:kidx1,:)) * qair_fact(1:nidx,:)
    csat_transpiration_tiles(1:nidx,:) = Surface%veg_ratio_actual_per_tile(kidx0:kidx1,:) * qsat_transpiration(1:nidx,:)

    ! ECHAM5 compatibility: one surface temperature for whole grid box

    CALL average_tiles(csat_tiles(1:nidx,:), surface%is_present(kidx0:kidx1,:) .AND. .NOT.                &
                       surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), csat(1:nidx))
    CALL average_tiles(cair_tiles(1:nidx,:), surface%is_present(kidx0:kidx1,:) .AND. .NOT.                &
                       surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), cair(1:nidx))
    CALL average_tiles(csat_transpiration_tiles(1:nidx,:), surface%is_present(kidx0:kidx1,:) .AND. .NOT.  &
                       surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), csat_transpiration(1:nidx))

    WHERE (SPREAD(air_moisture(1:nidx), NCOPIES=ntiles, DIM=2) < sat_surface_specific_hum_old(1:nidx,:))
       zhsoil = soil%snow_fract(kidx0:kidx1,:) +                                 & 
                (1._dp - soil%snow_fract(kidx0:kidx1,:)) *                       &
                (wet_skin_fract(1:nidx,:) + (1._dp - wet_skin_fract(1:nidx,:)) * &
                relative_humidity(1:nidx,:))
    ELSEWHERE
       zhsoil = 1._dp
    END WHERE
    CALL average_tiles(zhsoil, surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
         surface%cover_fract(kidx0:kidx1,:), zhsoil_avg)

    soil%csat(kidx0:kidx1) = csat(1:nidx)
    soil%cair(kidx0:kidx1) = cair(1:nidx)
    soil%csat_transpiration(kidx0:kidx1) = csat_transpiration(1:nidx)

!---wiso-code
    IF (lwiso) THEN

    ! Land surface saturation specific humidity - water isotopes
    ! calulate equilibrium and kinetic fractionation factor
    DO itile = 1,ntiles
       soil%remember_cveg(kidx0:kidx1,itile)      = Surface%veg_ratio_actual_per_tile(kidx0:kidx1,itile)
       soil%remember_canopy(kidx0:kidx1,itile)    = canopy_resistance(1:nidx,itile)
       soil%remember_rel_hum(kidx0:kidx1,itile)   = relative_humidity(1:nidx,itile)
       soil%remember_echamzchl(kidx0:kidx1,itile) = p_echam_zchl(1:nidx)
       soil%remember_cfrac(kidx0:kidx1,itile)     = surface%cover_fract(kidx0:kidx1,itile)
       soil%remember_snfrac(kidx0:kidx1,itile)    = soil%snow_fract(kidx0:kidx1,itile)
       soil%remember_wlfrac(kidx0:kidx1,itile)    = wet_skin_fract(1:nidx,itile)
       soil%remember_wind(kidx0:kidx1,itile)      = windspeed(1:nidx)
       soil%remember_airm(kidx0:kidx1,itile)      = air_moisture(1:nidx)
       soil%remember_satsphum(kidx0:kidx1,itile)  = sat_surf_specific_hum_avg(1:nidx)
    END DO
    
    END IF
!---wiso-code-end

    !------------------------------------------------------------------------------------------------------------
    ! END RECALC QSAT QAir and zhsoil
    !------------------------------------------------------------------------------------------------------------
    ! Note: glacier_precip_minus_evap is in m water equivalent; to convert to kg/m^2s multiply by RhoH2O/delta_time,
    !       and for accumulation this is again multiplied by delta_time, i.e. just multiply by RhoH2O
    ! Accumulate P-E for glaciers
    WHERE (surface%is_present(kidx0:kidx1,:))
       soil%glacier_precip_minus_evap_acc(kidx0:kidx1,:) = soil%glacier_precip_minus_evap_acc(kidx0:kidx1,:) + &
                                                           glacier_precip_minus_evap * RhoH2O
    END WHERE

    ! INPUT for HD-Model in coupled case:----------
    glac_runoff_evap = 0._dp
    surf_runoff_hd = 0._dp
    drainage_hd = 0._dp

    CALL average_tiles(glacier_precip_minus_evap(1:nidx,1:ntiles), surface%is_present(kidx0:kidx1,1:ntiles), & 
                       surface%cover_fract(kidx0:kidx1,1:ntiles), glac_runoff_evap(1:nidx))
    CALL average_tiles(surface_runoff(1:nidx,1:ntiles), surface%is_present(kidx0:kidx1,1:ntiles), &
                       surface%cover_fract(kidx0:kidx1,1:ntiles), surf_runoff_hd(1:nidx))
    CALL average_tiles(drainage(1:nidx,1:ntiles), surface%is_present(kidx0:kidx1,1:ntiles), &
                       surface%cover_fract(kidx0:kidx1,1:ntiles), drainage_hd(1:nidx))

!---wiso-code
    IF (lwiso) THEN

    DO jt=1,nwiso
      WHERE (surface%is_present(kidx0:kidx1,:))
        soil%wiso_glac_p_minus_e_acc(kidx0:kidx1,jt,:) = soil%wiso_glac_p_minus_e_acc(kidx0:kidx1,jt,:) + &
                                                           wiso_glac_p_minus_e(:,jt,:) * twisorhoh2o(jt)
      END WHERE
    END DO


    ! INPUT for HD-Model in coupled case:----------
    wiso_glac_runoff_evap = 0._dp
    wiso_surf_runoff_hd = 0._dp
    wiso_drainage_hd = 0._dp

    DO jt=1,nwiso
      CALL average_tiles(wiso_glac_p_minus_e(1:nidx,jt,1:ntiles), surface%is_present(kidx0:kidx1,1:ntiles), & 
                       surface%cover_fract(kidx0:kidx1,1:ntiles), wiso_glac_runoff_evap(1:nidx,jt))
      CALL average_tiles(wiso_surface_runoff(1:nidx,jt,1:ntiles), surface%is_present(kidx0:kidx1,1:ntiles), &
                       surface%cover_fract(kidx0:kidx1,1:ntiles), wiso_surf_runoff_hd(1:nidx,jt))
      CALL average_tiles(wiso_drainage(1:nidx,jt,1:ntiles), surface%is_present(kidx0:kidx1,1:ntiles), &
                       surface%cover_fract(kidx0:kidx1,1:ntiles), wiso_drainage_hd(1:nidx,jt))
    END DO
    
    END IF
!---wiso-code-end

    IF (useDynveg) THEN
       air_qsat(1:nidx) = sat_specific_humidity(surface_temperature_avg(1:nidx), surface_pressure(1:nidx))
       IF (ANY(air_qsat(1:nidx) > 0.99_dp*HUGE(1._dp))) &
          CALL message('sat_specific_humidity', 'lookup table overflow')
       soil%relative_humidity_air(kidx0:kidx1) = 100._dp * air_moisture(1:nidx) * &
         ((1._dp - GasConstantDryAir / GasConstantWaterVapor) * &
         air_qsat(1:nidx) + (GasConstantDryAir / GasConstantWaterVapor)) / &
         (air_qsat(1:nidx) * ((1._dp - GasConstantDryAir / GasConstantWaterVapor) * air_moisture(1:nidx) + &
         (GasConstantDryAir / GasConstantWaterVapor)))
    END IF

  END SUBROUTINE update_soil

  FUNCTION calc_moist_root_zone(soil_moisture, soil_depth, root_depth) RESULT(root_moisture)

    ! Calculate moisture in root zone
    ! This assumes for now that depths of root zones are equal to depths of soil moisture layers, 
    ! except the last root zone

    REAL(dp), INTENT(in) :: soil_moisture(:,:)    ! nland x nsoil
    REAL(dp), INTENT(in) :: soil_depth(:,:)       ! nland x nsoil
    REAL(dp), INTENT(in) :: root_depth(:,:)       ! nland x nroot_zones
    REAL(dp)             :: root_moisture(SIZE(soil_moisture,1))

    INTEGER :: n_soil, n_root

    n_soil = SIZE(soil_moisture,2)
    n_root = SIZE(root_depth,2)

    root_moisture(:) = SUM(soil_moisture(:,1:n_root-1), DIM=2)
    root_moisture(:) = root_moisture(:) + soil_moisture(:,n_root) * root_depth(:,n_root) / soil_depth(:,n_root)

  END FUNCTION calc_moist_root_zone

  ELEMENTAL FUNCTION calc_water_stress_factor(moisture, moisture_max, fract_critical, fract_wilting) RESULT(stress)

    REAL(dp), INTENT(in)  :: moisture
    REAL(dp), INTENT(in)  :: moisture_max
    REAL(dp), INTENT(in)  :: fract_critical, fract_wilting
    REAL(dp)              :: stress

    REAL(dp) :: moisture_critical, moisture_wilting

    moisture_critical = fract_critical * moisture_max
    moisture_wilting  = fract_wilting  * moisture_max

    stress = MAX(0._dp, MIN(1._dp, (moisture - moisture_wilting) / (moisture_critical - moisture_wilting) ))

  END FUNCTION calc_water_stress_factor

  FUNCTION calc_relative_humidity(moisture, moisture_max, fract_wilting) RESULT(rel_humidity)

    USE mo_constants,          ONLY: api

    REAL(dp), INTENT(in)  :: moisture(:,:)
    REAL(dp), INTENT(in)  :: moisture_max(:,:)
    REAL(dp), INTENT(in)  :: fract_wilting
    REAL(dp)              :: rel_humidity(SIZE(moisture,1),SIZE(moisture,2))

    REAL(dp) moisture_limit(SIZE(moisture,1),SIZE(moisture,2)), evap_stop(SIZE(moisture,1),SIZE(moisture,2))

    evap_stop      = MIN(0.1_dp, moisture_max)
    moisture_limit = moisture_max - evap_stop

    WHERE (moisture > moisture_limit .AND. moisture > fract_wilting*moisture_max)
       rel_humidity = 0.5_dp * (1._dp - COS((moisture-moisture_limit) * api / evap_stop))
    ELSEWHERE
       rel_humidity = 0._dp
    END WHERE

  END FUNCTION calc_relative_humidity

  SUBROUTINE soil_diagnostics(surface, soil)

    USE mo_land_surface, ONLY: land_surface_type
    USE mo_utils,        ONLY: average_tiles
    USE mo_time_control, ONLY: delta_time, l_putdata
!---wiso-code
    USE mo_wiso,         ONLY: lwiso, nwiso
!---wiso-code-end

    TYPE(land_surface_type), INTENT(in) :: surface
    TYPE(soil_type),         INTENT(inout) :: soil

    !Local variables
    LOGICAL  :: mask(nidx,soil%ntiles)
    REAL(dp) :: fract(nidx,soil%ntiles)
    INTEGER  :: itile, ntiles, i, jt
    INTEGER  :: kidx0, kidx1

    ntiles  = soil%ntiles
    kidx0   = kstart
    kidx1   = kend

    ! Accumulate variables
    WHERE (surface%is_present(kidx0:kidx1,:))
       ! Ground heat flux
       soil%ground_heat_flux_acc(kidx0:kidx1,:) = soil%ground_heat_flux_acc(kidx0:kidx1,:) + &
            soil%ground_heat_flux(kidx0:kidx1,:) * delta_time
       ! Evapotranspiration, transpiration and latent/sensible heat fluxes
       soil%evapotranspiration_acc(kidx0:kidx1,:) = soil%evapotranspiration_acc(kidx0:kidx1,:) + &
            soil%evapotranspiration(kidx0:kidx1,:) * delta_time
       soil%transpiration_acc(kidx0:kidx1,:) = soil%transpiration_acc(kidx0:kidx1,:) + &
            soil%transpiration(kidx0:kidx1,:) * delta_time
       soil%latent_heat_acc(kidx0:kidx1,:) = soil%latent_heat_acc(kidx0:kidx1,:) + &
            soil%latent_heat_flux(kidx0:kidx1,:) * delta_time
       soil%sensible_heat_acc(kidx0:kidx1,:) = soil%sensible_heat_acc(kidx0:kidx1,:) + &
            soil%sensible_heat_flux(kidx0:kidx1,:) * delta_time
    END WHERE

 !---wiso-code
    IF (lwiso) THEN

   ! Accumulate variables - Water Isotopes
    DO jt=1,nwiso
       WHERE (surface%is_present(kidx0:kidx1,:))
          ! Evapotranspiration, transpiration and latent/sensible heat fluxes
          soil%wiso_evapotranspiration_acc(kidx0:kidx1,jt,:) = soil%wiso_evapotranspiration_acc(kidx0:kidx1,jt,:) + &
               soil%wiso_evapotranspiration(kidx0:kidx1,jt,:) * delta_time
          soil%wiso_transpiration_acc(kidx0:kidx1,jt,:) = soil%wiso_transpiration_acc(kidx0:kidx1,jt,:) + &
               soil%wiso_transpiration(kidx0:kidx1,jt,:) * delta_time
       END WHERE
    END DO

    END IF
!---wiso-code-end

    ! Compute grid box averages
    mask  = surface%is_present(kidx0:kidx1,1:ntiles)
    fract = surface%cover_fract(kidx0:kidx1,1:ntiles)

    CALL average_tiles(soil%surface_temperature(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%surface_temperature(kidx0:kidx1))

    CALL average_tiles(soil%radiative_temperature(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%surface_radiative_temp(kidx0:kidx1))
    
    CALL average_tiles(soil%sat_surface_specific_humidity(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%sat_surface_specific_humidity(kidx0:kidx1))

    CALL average_tiles(soil%glacier_runoff_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%glacier_runoff_acc(kidx0:kidx1))

!!$    ! Variables that are not passed back to the atmosphere in each time step are only averaged at the output time step
!!$    IF (.NOT. l_putdata(IO_diag_stream%post_idx)) RETURN

    CALL average_tiles(soil%ground_heat_flux_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%ground_heat_flux_acc(kidx0:kidx1))

    CALL average_tiles(soil%heat_capacity(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%heat_capacity(kidx0:kidx1))

    CALL average_tiles(soil%evapotranspiration_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%evapotranspiration_acc(kidx0:kidx1))

    CALL average_tiles(soil%transpiration_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%transpiration_acc(kidx0:kidx1))

    CALL average_tiles(soil%sensible_heat_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%sensible_heat_acc(kidx0:kidx1))

    CALL average_tiles(soil%latent_heat_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%latent_heat_acc(kidx0:kidx1))

    CALL average_tiles(soil%albedo(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%albedo(kidx0:kidx1))

    CALL average_tiles(soil%skin_reservoir(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%skin_reservoir(kidx0:kidx1))

    CALL average_tiles(soil%snow(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%snow(kidx0:kidx1))

    CALL average_tiles(soil%snow_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%snow_acc(kidx0:kidx1))

    CALL average_tiles(soil%snow_fract(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%snow_fract(kidx0:kidx1))

    CALL average_tiles(soil%snow_melt_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%snow_melt_acc(kidx0:kidx1))

    CALL average_tiles(soil%runoff_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%runoff_acc(kidx0:kidx1))

    CALL average_tiles(soil%drainage_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%drainage_acc(kidx0:kidx1))

    CALL average_tiles(soil%glacier_depth(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%glacier_depth(kidx0:kidx1))

    CALL average_tiles(soil%glacier_precip_minus_evap_acc(kidx0:kidx1,1:ntiles), mask, fract, &
         soil_diag%glacier_precip_minus_evap_acc(kidx0:kidx1))

    DO i=1,nsoil
       CALL average_tiles(soil%moisture(kidx0:kidx1,i,1:ntiles), mask, fract, &
            soil_diag%moisture(kidx0:kidx1,i))
    END DO

    DO i=1,ntsoil
       CALL average_tiles(soil%soil_temperature(kidx0:kidx1,i,1:ntiles), mask, fract, &
            soil_diag%soil_temperature(kidx0:kidx1,i))
    END DO

!---wiso-code
    IF (lwiso) THEN

    ! Accumulate variables over tiles - Water Isotopes
    DO jt=1,nwiso
       CALL average_tiles(soil%wiso_evapotranspiration_acc(kidx0:kidx1,jt,1:ntiles), mask, fract, &
            soil_diag%wiso_evapotranspiration_acc(kidx0:kidx1,jt))
       CALL average_tiles(soil%wiso_transpiration_acc(kidx0:kidx1,jt,1:ntiles), mask, fract, &
            soil_diag%wiso_transpiration_acc(kidx0:kidx1,jt))
       CALL average_tiles(soil%wiso_skin_reservoir(kidx0:kidx1,jt,1:ntiles), mask, fract, &
            soil_diag%wiso_skin_reservoir(kidx0:kidx1,jt))
       CALL average_tiles(soil%wiso_snow(kidx0:kidx1,jt,1:ntiles), mask, fract, soil_diag%wiso_snow(kidx0:kidx1,jt))
       CALL average_tiles(soil%wiso_snow_acc(kidx0:kidx1,jt,1:ntiles), mask, fract, &
            soil_diag%wiso_snow_acc(kidx0:kidx1,jt))
       CALL average_tiles(soil%wiso_snow_melt_acc(kidx0:kidx1,jt,1:ntiles), mask, fract,&
            soil_diag%wiso_snow_melt_acc(kidx0:kidx1,jt))
       CALL average_tiles(soil%wiso_runoff_acc(kidx0:kidx1,jt,1:ntiles), mask, fract, &
            soil_diag%wiso_runoff_acc(kidx0:kidx1,jt))
       CALL average_tiles(soil%wiso_drainage_acc(kidx0:kidx1,jt,1:ntiles), mask, fract, &
            soil_diag%wiso_drainage_acc(kidx0:kidx1,jt))
       CALL average_tiles(soil%wiso_sat_surface_specific_hum(kidx0:kidx1,jt,1:ntiles), mask, fract, &
            soil_diag%wiso_sat_surface_specific_hum(kidx0:kidx1,jt))
       CALL average_tiles(soil%wiso_moisture(kidx0:kidx1,jt,1:ntiles), mask, fract, &
            soil_diag%wiso_moisture(kidx0:kidx1,jt))
       CALL average_tiles(soil%wiso_glacier_depth(kidx0:kidx1,jt,1:ntiles), mask, fract, &
            soil_diag%wiso_glacier_depth(kidx0:kidx1,jt))
       CALL average_tiles(soil%wiso_glacier_runoff_acc(kidx0:kidx1,jt,1:ntiles), mask, fract, &
            soil_diag%wiso_glacier_runoff_acc(kidx0:kidx1,jt))
       CALL average_tiles(soil%wiso_glac_p_minus_e_acc(kidx0:kidx1,jt,1:ntiles), mask, fract, &
            soil_diag%wiso_glac_p_minus_e_acc(kidx0:kidx1,jt))
    END DO
    
    END IF
!---wiso-code-end

  END SUBROUTINE soil_diagnostics

  FUNCTION get_soil_diag(kdim, what) RESULT(outfield)

    ! only 1-d at the moment

    INTEGER, INTENT(in) :: kdim
    CHARACTER(len=*), INTENT(in) :: what
    REAL(dp) :: outfield(kdim)

    INTEGER :: k0, k1

    k0 = kstart
    k1 = kend


    IF (kdim /= (k1-k0+1)) CALL finish('get_soil_diag','Wrong dimension')

    SELECT CASE(TRIM(what))
    CASE('surface_temperature')
       outfield(1:kdim) = soil_diag%surface_temperature(k0:k1)
    CASE('surface_radiative_temperature')
       outfield(1:kdim) = soil_diag%surface_radiative_temp(k0:k1)
    CASE('sat_surface_specific_humidity')
       outfield(1:kdim) = soil_diag%sat_surface_specific_humidity(k0:k1)
    CASE('runoff_acc')
       outfield(1:kdim) = soil_diag%runoff_acc(k0:k1)
    CASE('glacier_depth')
       outfield(1:kdim) = soil_diag%glacier_depth(k0:k1)
    CASE('snow_melt_acc')
       outfield(1:kdim) = soil_diag%snow_melt_acc(k0:k1)
    CASE('glacier_precip_minus_evap_acc')
       outfield(1:kdim) = soil_diag%glacier_precip_minus_evap_acc(k0:k1)
    CASE('snow_acc')
       outfield(1:kdim) = soil_diag%snow_acc(k0:k1)
    CASE('glacier_runoff_acc')
       outfield(1:kdim) = soil_diag%glacier_runoff_acc(k0:k1)
    CASE default
       CALL finish('get_soil_diag','Unknown request')
    END SELECT

  END FUNCTION get_soil_diag

!---wiso-code
  FUNCTION get_soil_diag_wiso(kdim, kwiso, what) RESULT(outfield)

    ! only 1-d at the moment

    INTEGER, INTENT(in) :: kdim, kwiso
    CHARACTER(len=*), INTENT(in) :: what
    REAL(dp) :: outfield(kdim,kwiso)

    INTEGER :: k0, k1, jt

    k0 = kstart
    k1 = kend


    IF (kdim /= (k1-k0+1)) CALL finish('get_soil_diag_wiso','Wrong dimension')

    DO jt=1,kwiso
       SELECT CASE(TRIM(what))
       CASE('wiso_sat_surface_specific_hum')
          outfield(1:kdim,jt) = soil_diag%wiso_sat_surface_specific_hum(k0:k1,jt)
       CASE('wiso_runoff_acc')
          outfield(1:kdim,jt) = soil_diag%wiso_runoff_acc(k0:k1,jt)
       CASE('wiso_glacier_depth')
          outfield(1:kdim,jt) = soil_diag%wiso_glacier_depth(k0:k1,jt)
       CASE('wiso_snow_melt_acc')
          outfield(1:kdim,jt) = soil_diag%wiso_snow_melt_acc(k0:k1,jt)
       CASE('wiso_glac_p_minus_e_acc')
          outfield(1:kdim,jt) = soil_diag%wiso_glac_p_minus_e_acc(k0:k1,jt)
       CASE('wiso_snow_acc')
          outfield(1:kdim,jt) = soil_diag%wiso_snow_acc(k0:k1,jt)
       CASE('wiso_glacier_runoff_acc')
          outfield(1:kdim,jt) = soil_diag%wiso_glacier_runoff_acc(k0:k1,jt)
       CASE default
          CALL finish('get_soil_diag_wiso','Unknown request')
       END SELECT
    END DO

  END FUNCTION get_soil_diag_wiso

!---wiso-code-end

  SUBROUTINE ini_soil_temp(klon,ntsoil,tslclim,tsoil,tsl)
    !
    ! Initialize soil and surface temperature 
    !
    ! Method:
    !
    ! Starting from the tslclim field temperatures are set in relation
    ! to the depth of the soil layer and position of the initial
    ! day in the annual cycle.
    !
    ! tsl is at 0.07 m
    ! thickness of layers 0.065, 0.254, 0.913, 2.902, 5.700 m (s. soiltemp)
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, April 2000, original source
    ! L. Dumenil, MPI, June 1989, original source of *initemp*
    !             which is now part of this routine
    ! I. Kirchner, MPI, November 2000, date/time control
    ! U. Schulzweida, MPI, May 2002, blocking (nproma)
    ! E. Roeckner, MPI, Sep 2002, initialization of mixed layer ocean
    ! T. Raddatz, MPI, Mai 2005, adaption to JSBACH
    ! 
    USE mo_constants,          ONLY: api
    USE mo_radiation,          ONLY: nmonth
    USE mo_time_control,       ONLY: get_date_components, ndaylen, start_date, get_month_len, get_year_len

    INTEGER, INTENT(in)  ::  klon                 !! number of land grid points
    INTEGER, INTENT(in)  ::  ntsoil               !! number of thermal soil layers
    REAL(dp), INTENT(in)     ::  tslclim(klon,12)     !! climatological surface temperature for each kalendar month
    REAL(dp), INTENT(out)    ::  tsoil(klon,ntsoil)   !! temperature of the soil
    REAL(dp), INTENT(out)    ::  tsl(klon)            !! surface temperature at time step t

    !  Local
    REAL(dp) :: zdmax(klon), zkap, zsqrt, zday, yearl
    REAL(dp) :: zd(ntsoil), znmea(klon), zrange(klon)
    INTEGER  ::  jg, jl, im
    INTEGER  ::  iyday
    INTEGER  ::  jmax(klon), nmomid(12), nmonthl(12)
    INTEGER  ::  yr, mo, dy, hr, mn, se 

    ! Computational constants

    zkap = 7.5E-7_dp

    !
    !    Date handling
    !

    ! get year, month, day of start date

    CALL get_date_components(start_date, yr, mo, dy, hr, mn, se)

    ! get length of the year

    yearl = get_year_len(yr)

    ! get length of month

    DO im = 1,12
       nmonthl(im) = REAL(get_month_len(yr,im),dp)
       nmomid(im)  = nmonthl(im) * 0.5_dp
    END DO

    ! calendar day of the start date

    IF (nmonth == 0) THEN
       ! annual cycle run
       iyday = dy
       DO im = 1, mo-1
          iyday = iyday + nmonthl(im)
       END DO
    ELSE
       ! perpetual months
       iyday = nmomid(nmonth)
       DO im = 1, nmonth-1
          iyday = iyday + nmonthl(im)
       END DO
    END IF
    zday = REAL(iyday,dp)
    zsqrt = SQRT(zkap * yearl * REAL(ndaylen,dp) / api)

    ! calendar month of maximum surface temperature

    jmax(:) = MAXLOC(tslclim(:,:),DIM=2)

    ! difference between months of maximum and minimum surface temperature

    zrange(:) = MAXVAL(tslclim(:,:),DIM=2) - MINVAL(tslclim(:,:),DIM=2)

    ! calendar day of maximum surface temperature 

    zdmax(:) = 0._dp
    DO jl= 1,klon
       DO im = 1,jmax(jl)
          zdmax(jl) = zdmax(jl) + REAL(nmonthl(im),dp)
       END DO
       zdmax(jl) = zdmax(jl) - REAL(nmomid(jmax(jl)),dp)
    END DO

    !
    !    Initialize soil temperatures
    !

    ! layer depths

    zd(1) =            0.5_dp * 0.065_dp - 0.07_dp
    zd(2) = 0.065_dp + 0.5_dp * 0.254_dp - 0.07_dp
    zd(3) = 0.065_dp + 0.254_dp + 0.5_dp * 0.913_dp - 0.07_dp
    zd(4) = 0.065_dp + 0.254_dp + 0.913_dp + 0.5_dp * 2.902_dp - 0.07_dp
    zd(5) = 0.065_dp + 0.254_dp + 0.913_dp + 2.902_dp + 0.5_dp * 5.7_dp - 0.07_dp

    ! calculate annual mean surface temperature

    znmea(:) = SUM(tslclim(:,:),DIM=2) / 12._dp

    ! set soil temperature of the five layers

    DO jg = 1,ntsoil
       tsoil(:,jg) = znmea(:) + 0.5_dp * zrange(:) * EXP(-zd(jg) / zsqrt)       &
            * COS(2._dp * api * (zday - zdmax(:)) / yearl - zd(jg) / zsqrt)
    END DO

    !  set all time levels of surface temperature to uppermost soil temp.

    tsl(:)   = tsoil(:,1)

    RETURN

  END SUBROUTINE ini_soil_temp

End module mo_soil

