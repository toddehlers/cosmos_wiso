MODULE mo_land_boundary

  USE mo_jsbach_grid, ONLY: kstart, kend, nidx
  USE mo_jsbach_lctlib, ONLY: lctlib_type
  USE mo_land_surface, ONLY: land_surface_type
  USE mo_soil, ONLY: soil_type
  USE mo_utils, ONLY: average_tiles
  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  REAL(dp), PARAMETER :: blending_height = 100._dp     ! blending height [m]
  REAL(dp), PARAMETER :: roughness_snow = 0.001_dp     ! roughness length of surfaces covered by snow [m]
  REAL(dp), PARAMETER :: roughness_bare = 0.005_dp     ! roughness length of bare land [m]

CONTAINS

  SUBROUTINE update_land_boundary_up(lctlib, land_surface, soil,                                            &
      roughness_oro, roughness_heat, roughness_momentum,                                                    &
      albedo, surface_temperature, sat_surface_specific_humidity, evapotranspiration, evaporation_pot,      &
      sensible_heat_flux_avg, latent_heat_flux_avg, radiative_temp_avg,                                     &
      dry_static_energy_new_avg, soil_moisture_avg, snow_depth_avg,                                         &
      runoff_avg, drainage_avg, skin_res_avg,                                                               &
!---wiso-code
      wiso_sat_surface_specific_hum, wiso_evapotranspiration, wiso_evaporation_pot,                         &
      wiso_soil_moisture_avg, wiso_snow_depth_avg,                                                          &
      wiso_runoff_avg,   wiso_drainage_avg, wiso_skin_res_avg)
!---wiso-code-end

    USE mo_soil, ONLY: get_soil_diag, &
!---wiso-code
                       get_soil_diag_wiso
    USE mo_wiso, ONLY: lwiso, nwiso, tnat
!---wiso-code-end

    TYPE(lctlib_type),       INTENT(in) :: lctlib
    TYPE(land_surface_type), INTENT(in) :: land_surface
    TYPE(soil_type),         INTENT(in) :: soil
    REAL(dp),                INTENT(in), DIMENSION(nidx) :: &
         roughness_oro                         !! Surface roughness due to orographie
    REAL(dp),                INTENT(out), DIMENSION(nidx) :: &
         roughness_heat,                    &  !! Surface roughness of heat
         roughness_momentum,                &  !! Surface roughness of momentum
         albedo,                            &  !! Land albedo for grid box
         surface_temperature,               &  !! Surface temperature for grid box
         sat_surface_specific_humidity,     &  !! Saturated surface specific humidity for grid box
         evapotranspiration,                &  !! Evapotranspiration for grid box
         evaporation_pot,                   &  !! Potential evaporation for grid box
         sensible_heat_flux_avg,            &  !! sensible heat flux from surface (W/m2)
         latent_heat_flux_avg,              &  !! latent heat flux from surface (W/m2)
         radiative_temp_avg,                &
         dry_static_energy_new_avg,         &
         soil_moisture_avg, snow_depth_avg, &
         runoff_avg, drainage_avg,          &
         skin_res_avg

!---wiso-code
    REAL(dp),                INTENT(out), DIMENSION(nidx,nwiso), OPTIONAL :: &
         wiso_skin_res_avg,                 & !! Water content of skin reservoir for grid box - Isotopes
         wiso_snow_depth_avg,               & !! Snow on ground for grid box - Isotopes
         wiso_drainage_avg,                 & !! Drainage at non-glacier for grid box - Isotopes
         wiso_evapotranspiration,           & !! Evapotranspiration for grid box - Isotopes
         wiso_evaporation_pot,              & !! Potential evaporation for grid box - Isotopes
         wiso_soil_moisture_avg,            & !! Moisture content of the unfrozen sub-layer for grid box - Isotopes
         wiso_sat_surface_specific_hum,     & !! Saturated surface specific humidity for grid box - Isotopes
         wiso_runoff_avg                      !! surface runoff - Isotopes
!---wiso-code-end

    REAL(dp)  ::  roughness_sum(nidx)
    INTEGER   ::  kidx0, kidx1, ntiles, i, j, k, jt, jl

    kidx0 = kstart
    kidx1 = kend
    ntiles = land_surface%ntiles
    
    ! Surface roughness
    DO i = 1,nidx
      j = kidx0 - 1 + i
      roughness_sum(i) = 0._dp
      DO k = 1,ntiles
        roughness_sum(i) = roughness_sum(i) + ((1._dp - land_surface%veg_ratio_max(j)) *      & ! bare land
                           land_surface%cover_fract(j,k) /                                    &
                           (LOG(blending_height / roughness_bare) ** 2))
        roughness_sum(i) = roughness_sum(i) + (land_surface%veg_ratio_max(j) *                & ! vegetated land
                           land_surface%cover_fract(j,k) /                                    &
                           (LOG(blending_height / lctlib%VegRoughness(land_surface%cover_type(j,k))) ** 2))
      END DO
      roughness_momentum(i) = SQRT((blending_height * EXP(-1._dp / SQRT(roughness_sum(i)))) ** 2 + roughness_oro(i) ** 2)
      roughness_sum(i) = 0._dp
      DO k = 1,ntiles
        roughness_sum(i) = roughness_sum(i) + (soil%snow_fract(j,k) *                         & ! snow covered land
                           land_surface%cover_fract(j,k) /                                    &
                           (LOG(blending_height / roughness_snow) ** 2))
        roughness_sum(i) = roughness_sum(i) + ((1._dp - soil%snow_fract(j,k)) *               & ! not snow covered land
                           land_surface%cover_fract(j,k) /                                    &
                           (LOG(blending_height / MIN(1._dp,roughness_momentum(i))) ** 2))
      END DO
      roughness_heat(i) = blending_height * EXP(-1._dp / SQRT(roughness_sum(i)))
    END DO


    surface_temperature(1:nidx) = get_soil_diag(nidx, 'surface_temperature')
    sat_surface_specific_humidity(1:nidx) = get_soil_diag(nidx, 'sat_surface_specific_humidity')

!---wiso-code
    IF (lwiso) THEN

    ! Land surface saturation specific humidity - water isotopes
    wiso_sat_surface_specific_hum(1:nidx,1:nwiso) = get_soil_diag_wiso(nidx,nwiso,'wiso_sat_surface_specific_hum')

    END IF
!---wiso-code-end

    radiative_temp_avg(1:nidx) = get_soil_diag(nidx, 'surface_radiative_temperature')
    runoff_avg(1:nidx) = get_soil_diag(nidx, 'runoff_acc')
!---wiso-code
    IF (lwiso) THEN

    wiso_runoff_avg(1:nidx,1:nwiso) = get_soil_diag_wiso(nidx,nwiso,'wiso_runoff_acc')

    END IF
!---wiso-code-end

    CALL average_tiles(land_surface%albedo(kidx0:kidx1,1:ntiles), land_surface%is_present(kidx0:kidx1,1:ntiles), &
                       land_surface%cover_fract(kidx0:kidx1,1:ntiles), albedo(:) )

    CALL average_tiles(soil%evapotranspiration(kidx0:kidx1,1:ntiles), land_surface%is_present(kidx0:kidx1,1:ntiles) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,1:ntiles), evapotranspiration(:))
    CALL average_tiles(soil%evaporation_pot(kidx0:kidx1,1:ntiles), land_surface%is_present(kidx0:kidx1,1:ntiles) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,1:ntiles), evaporation_pot(:))
    CALL average_tiles(soil%latent_heat_flux(kidx0:kidx1,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),latent_heat_flux_avg(:))
    CALL average_tiles(soil%sensible_heat_flux(kidx0:kidx1,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),sensible_heat_flux_avg(:))
    CALL average_tiles(soil%dry_static_energy_new(kidx0:kidx1,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),dry_static_energy_new_avg(:))
    CALL average_tiles(soil%moisture(kidx0:kidx1,1,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),soil_moisture_avg(:))
    CALL average_tiles(soil%snow(kidx0:kidx1,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),snow_depth_avg(:))
    CALL average_tiles(soil%drainage_acc(kidx0:kidx1,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),drainage_avg(:))
    CALL average_tiles(soil%skin_reservoir(kidx0:kidx1,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),skin_res_avg(:))

!---wiso-code
    IF (lwiso) THEN

    DO jt=1,nwiso
        CALL average_tiles(soil%wiso_skin_reservoir(kidx0:kidx1,jt,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),wiso_skin_res_avg(:,jt))
        CALL average_tiles(soil%wiso_snow(kidx0:kidx1,jt,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),wiso_snow_depth_avg(:,jt))
        CALL average_tiles(soil%wiso_drainage_acc(kidx0:kidx1,jt,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),wiso_drainage_avg(:,jt))
        CALL average_tiles(soil%wiso_evapotranspiration(kidx0:kidx1,jt,1:ntiles), land_surface%is_present(kidx0:kidx1,1:ntiles) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,1:ntiles), wiso_evapotranspiration(:,jt))
        CALL average_tiles(soil%wiso_evaporation_pot(kidx0:kidx1,jt,1:ntiles), land_surface%is_present(kidx0:kidx1,1:ntiles) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,1:ntiles), wiso_evaporation_pot(:,jt))
        CALL average_tiles(soil%wiso_moisture(kidx0:kidx1,jt,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),wiso_soil_moisture_avg(:,jt))
        CALL average_tiles(soil%wiso_runoff_acc(kidx0:kidx1,jt,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),wiso_runoff_avg(:,jt))
    END DO

    END IF
!---wiso-code-end

  END SUBROUTINE update_land_boundary_up

END MODULE mo_land_boundary
