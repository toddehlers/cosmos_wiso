!!+ Interface between main program and JSBACH
! 
MODULE mo_jsbach_interface

  ! 
  !! Description: 
  !!   This module defines the interface between a main program (climate model
  !!   or driver program).
  !!   The interface is independent of the forcing data which are provided from
  !!   the main program. It is specific to any one land surface scheme (JSBACH
  !!   in this case) and does all the soil_wetness conversions from the general grid-based
  !!   forcing input to the variables (structure, units, etc.) as needed by
  !!   the land surface scheme. In particular:
  !!   - Input fields are gathered to keep just land points
  !!   - Call the routines that make up the JSBACH land surface scheme
  !!   - Scatter output fields to complete global grids
  ! 
  !! Current Code Owner: jsbach_admin
  ! 
  !! History: 
  !  
  !! Version   Date       Comment 
  !! -------   ----       ------- 
  !! 0.1       2001/06/28 Original code (based on J. Polcher's
  !!                      routine intersurf.f90). Reiner Schnur
  ! 
  Use mo_jsbach,        Only : options_type,                        &
                               jsbach_FirstCall, jsbach_NewTimeStep,  &
                               debug, test_stream
  USE mo_soil,          ONLY : soil_type, soil_param_type, update_soil, soil_diagnostics, get_soil_diag, &
!---wiso-code
                               get_soil_diag_wiso
!---wiso-code-end
  USE mo_land_surface,  ONLY : land_surface_type, land_surface_diagnostics
  USE mo_jsbach_lctlib, ONLY : lctlib_type
  USE mo_jsbach_grid,   ONLY : grid_type, domain_type, kstart, kend, nidx, kindex
  Use mo_exception,     Only : finish, message,int2string,logical2string
  USE mo_doctor,        ONLY : nout
  Use mo_kind,          Only : dp
  USE mo_jsbach_veg,    ONLY : vegetation_type, veg_diagnostics
  USE mo_cbal_bethy,    Only : cbalance_type, init_cbalance_bethy, update_cbalance_bethy
  USE mo_bethy,         ONLY : bethy_type, init_bethy, update_bethy, bethy_diagnostics
  USE mo_util_string,   ONLY : tolower

  USE mo_mpi, ONLY : p_pe, p_io, p_parallel_io, p_parallel, p_bcast

  USE mo_linked_list, ONLY: t_stream
  
  USE mo_test, ONLY: init_test, init_test_2, test

  USE mo_time_control, ONLY : lresume
  USE mo_time_conversion, ONLY: time_days

  USE mo_dynveg, ONLY : init_dynveg, update_dynveg

#if defined (__SX__) && defined (_OPENMP)
  USE omp_lib,          ONLY: omp_get_thread_num, &
                              omp_get_num_threads
#endif

  IMPLICIT NONE

  PRIVATE

  ! Global (i.e. public) Declarations: 
  PUBLIC :: jsbach_interface, jsbach_inter_1d, update_current_call
  INTERFACE jsbach_interface
     MODULE PROCEDURE jsbach_inter_1d !, jsbach_inter_2d
  END INTERFACE


  TYPE land_type
     PRIVATE
     TYPE(grid_type)                :: Grid
     TYPE(domain_type)              :: Domain
     INTEGER                        :: nlct
     INTEGER                        :: npft
     INTEGER                        :: ntiles
     TYPE(land_surface_type)        :: Surface          ! State of land surface as seen by the atmosphere
     TYPE(lctlib_type)              :: LctLibrary       ! Lookup table of land cover characteristics (<nlct>)
     TYPE(soil_param_type)          :: SoilParam        ! Spatial distribution of soil characteristics
     TYPE(soil_type)                :: Soil             ! State of soil (<ntiles>)
     TYPE(vegetation_type)          :: Vegetation       ! State of vegetation (<ntiles>)
     TYPE(bethy_type)               :: Bethy            ! State of BETHY (including %Canopy for each canopy layer)
     TYPE(cbalance_type)            :: Cbalance         ! State of cbalance (NPP etc.)
     TYPE(t_stream), POINTER :: IO_stream
     TYPE(t_stream), POINTER :: IO_diag_stream
!---wiso-code
     TYPE(t_stream), POINTER :: IO_wiso_stream
     TYPE(t_stream), POINTER :: IO_diag_wiso_stream
!---wiso-code-end
     TYPE(t_stream), POINTER :: IO_veg
  END TYPE land_type

  PUBLIC :: jsbach_init, jsbach_restart

  TYPE(options_type), SAVE       :: theOptions
  TYPE(land_type), SAVE          :: theLand

  TYPE(time_days), SAVE          :: old_time_step      ! Used to check for jsbach_NewTimeStep
!$OMP THREADPRIVATE(old_time_step)


#ifdef STANDALONE
  EXTERNAL jsbalone_init_decomp
  EXTERNAL jsbalone_iniphy
#endif

CONTAINS 
  
  !!+ JSBACH interface routine for 1-d call
  SUBROUTINE jsbach_inter_1d ( &
       ! Input arguments
       kdim, &                            !! Length of vectors
       kland, &                           !! Number of land points in vectors ( = SUM(lmask) )
       ! First level conditions
       geopot, wind, wind10, temp_air, qair,  &
       ! Rain, snow, radiation, surface pressure
       precip_rain, precip_snow, &
       lwdown, &
       sw_vis_net, &                      !! net solar radiation in the visible band
       sw_vis_frac_diffuse, &             !! fraction of diffuse solar radiation contained in sw_vis_net
       sw_nir_net, &                      !! net solar radiation in the near infrared band
       sw_nir_frac_diffuse, &             !! fraction of diffuse solar radiation contained in sw_nir_net
       pressure, &
       czenith, &
       declination, &                     !! Solar declination 
       ! CO2 concentration
       CO2_concentration, &
       ! Output: Fluxes
       evap_act, evap_pot, sensible, latent, CO2_flux, &
       ! Surface temperatures and atm. humidity
       tsoil_rad, temp_soil_new, qsurf, &
       ! Surface radiative properties
       albedo, albedo_vis, albedo_nir, emis, &
       ! roughness length and drag coef.
       z0h, z0m, cdrag, &
       ! Variables for implicit coupling
       etAcoef, etBcoef, eqAcoef, eqBcoef, cair, csat, zhsoil, &
       echam_zchl, &
       mask_land, mode, surf_dry_static_energy, &
       kblock, soil_wetness, snow_depth, &
       runoff, drainage, skin_res, &
       tte_corr, &
       glac_runoff_evap, surf_runoff_hd , drainage_hd, &
       glacier_depth, snow_melt_acc, glacier_p_minus_e_acc, snow_acc, glacier_runoff_acc, &
       wsmax,                                                                               &
!---wiso-code
       ! Input - Isotopes
       lwisofracl,                                                                          &
       wiso_eqAcoef, wiso_eqBcoef, wiso_nenner,                                             & 
       wiso_helpqdif,                                                                       &
       wisoqair,                                                                            &
       wiso_precip_rain, wiso_precip_snow,                                                  &
       wisokinl,                                                                            &
       ! Output - Isotopes
       wisocair, wisocsat, wisocair_fra, wisocsat_fra,                                      &
       wiso_eqAcoef_new, wiso_eqBcoef_new,                                                  &
       wiso_evap_act, wiso_evap_pot,                                                        &
       wiso_qsurf,                                                                          &
       wiso_soil_wetness, wiso_snow_depth,                                                  &
       wiso_runoff, wiso_drainage, wiso_skin_res,                                           &
       wiso_glac_runoff_evap, wiso_surf_runoff_hd, wiso_drainage_hd,                        &
       wiso_glacier_depth, wiso_snow_melt_acc, wiso_glacier_p_minus_e_acc, wiso_snow_acc, wiso_glacier_runoff_acc       &
!---wiso-code-end
       )
 
    !! Description:

    USE mo_land_surface, ONLY: update_land_surface_fast, update_albedo, init_cover_fract, init_veg_ratio_max
    USE mo_land_boundary, ONLY: update_land_boundary_up
    USE mo_canopy, ONLY: unstressed_canopy_cond_par
    USE mo_jsbach_constants,    ONLY : RhoH2O
    USE mo_constants,    ONLY: tmelt

    USE mo_time_control, ONLY : delta_time, lstart, l_trigrad, l_trigradm1
    USE mo_phenology,    ONLY : update_phenology
    USE mo_utils,        ONLY : average_tiles
    USE mo_cbal_landcover_change, ONLY: do_landcover_change

!---wiso-code
    USE mo_wiso,                    ONLY : lwiso, nwiso
!---wiso-code-end

    ! Subroutine arguments
    INTEGER,            INTENT(in)    :: kdim                     !! Length of vectors (if call from ECHAM5 this
                                                                  !! will be nproma, or npromz for last block)
    INTEGER, OPTIONAL,  INTENT(in)    :: kblock                   !! Index of block  in domain that is to be
                                                                  !! processed. If missing: one call for whole domain
    INTEGER,            INTENT(in)    :: kland                    !! Number of land points in vectors
    LOGICAL, OPTIONAL,  INTENT(in)    :: mask_land(kdim)          !! Land-sea mask (land includes glaciers)
    CHARACTER(LEN=7), OPTIONAL, INTENT(in) :: mode
    !LOGICAL,            INTENT(in)    :: mask_glacier(kdim)      !! Mask for glacier boxes
    !LOGICAL,            INTENT(in)    :: mask_lake(kdim)         !! Mask for lakes
    REAL(dp), OPTIONAL, INTENT(in)    :: geopot(kdim)             !! Geopotential at lowest atmospheric level [m^2/s^2]
    !REAL,               INTENT(in)    :: uwind(kdim), vwind(kdim) !! Lowest level u/v wind [m/s]
    REAL(dp), OPTIONAL, INTENT(in)    :: wind(kdim)               !! Lowest level wind speed [m/s]
    REAL(dp), OPTIONAL, INTENT(in)    :: wind10(kdim)             !! 10m wind speed [m/s] (for update_surface_down)
    REAL(dp), OPTIONAL, INTENT(in)    :: temp_air(kdim)           !! Lowest level air temperature [Kelvin]
    REAL(dp), OPTIONAL, INTENT(in)    :: qair(kdim)               !! Lowest level specific humidity
    REAL(dp), OPTIONAL, INTENT(in)    :: etAcoef(kdim), &
                                         etBcoef(kdim), &
                                         eqAcoef(kdim), &
                                         eqBcoef(kdim)
    REAL(dp), OPTIONAL, INTENT(in)    :: cdrag(kdim)               !! Surface drag
    REAL(dp), OPTIONAL, INTENT(in)    :: precip_rain(kdim)         !! Precipitation as rain [kg/(m^2 s)]
    REAL(dp), OPTIONAL, INTENT(in)    :: precip_snow(kdim)         !! Precipitation as snow [kg/(m^2 s)]
    REAL(dp), OPTIONAL, INTENT(in)    :: CO2_concentration(kdim)   !! Atmospheric CO2 concentration [kg(CO2)/kg(air)]
    REAL(dp), OPTIONAL, INTENT(in)    :: lwdown(kdim)              !! Downward longwave flux
    !REAL(dp),           INTENT(in)    :: swdown(kdim)             !! Downward shortwave flux
    !REAL(dp),           INTENT(in)    :: swnet(kdim)              !! Net surface shortwave flux (for land part only!)
    REAL(dp), OPTIONAL, INTENT(in)    :: sw_vis_net(kdim)          !! net surface visible radiation [W/m^2]
    REAL(dp), OPTIONAL, INTENT(in)    :: sw_vis_frac_diffuse(kdim) !! fraction of diffuse radiation contained in sw_vis_net
    REAL(dp), OPTIONAL, INTENT(in)    :: sw_nir_net(kdim)          !! net surface NIR radiation [W/m^2]
    REAL(dp), OPTIONAL, INTENT(in)    :: sw_nir_frac_diffuse(kdim) !! fraction of diffuse radiation contained in sw_nir_net
    REAL(dp), OPTIONAL, INTENT(in)    :: czenith(kdim)             !! Cosine of solar zenith angle
    REAL(dp), OPTIONAL, INTENT(in)    :: declination               !! Solar declination
    REAL(dp), OPTIONAL, INTENT(in)    :: pressure(kdim)            !! Surface pressure

!---wiso-code
    INTEGER,  OPTIONAL, INTENT(in)    :: lwisofracl
    REAL(dp), OPTIONAL, INTENT(in)    :: wiso_helpqdif(kdim)
    REAL(dp), OPTIONAL, INTENT(in)    :: wiso_eqAcoef(kdim,nwiso), &
                                         wiso_eqBcoef(kdim,nwiso), &
                                         wiso_nenner(kdim,nwiso), &
                                         wisoqair(kdim,nwiso)
    REAL(dp), OPTIONAL, INTENT(in)    :: wiso_precip_rain(kdim,nwiso)   !! Precipitation as rain [kg/(m^2 s)] - Isotopes
    REAL(dp), OPTIONAL, INTENT(in)    :: wiso_precip_snow(kdim,nwiso)   !! Precipitation as snow [kg/(m^2 s)] - Isotopes
    REAL(dp), OPTIONAL, INTENT(in)    :: wisokinl(kdim,nwiso)           !! kinetic fractionation factor
!---wiso-code-end

    REAL(dp), OPTIONAL, INTENT(out)   :: evap_act(kdim)            !! Total of evaporation (actual) [kg/(m^2 s)]
    REAL(dp), OPTIONAL, INTENT(out)   :: evap_pot(kdim)            !! Potential evaporation [kg/(m^2 s)]
    REAL(dp), OPTIONAL, INTENT(out)   :: sensible(kdim)            !! Sensible heat flux [W/m^2]
    REAL(dp), OPTIONAL, INTENT(out)   :: latent(kdim)              !! Latent heat flux [W/m^2]
    REAL(dp), OPTIONAL, INTENT(out)   :: CO2_flux(kdim)            !! CO2 flux [kg(CO2)/m^2 s] in each time step
    REAL(dp), OPTIONAL, INTENT(out)   :: tsoil_rad(kdim)           !!
    REAL(dp), OPTIONAL, INTENT(out)   :: temp_soil_new(kdim)       !! New soil temperature
    REAL(dp), OPTIONAL, INTENT(out)   :: qsurf(kdim)               !! Surface specific humidity
    REAL(dp), OPTIONAL, INTENT(out)   :: albedo(kdim)              !! Albedo
    REAL(dp), OPTIONAL, INTENT(out)   :: albedo_vis(kdim)          !! Albedo of the visible range
    REAL(dp), OPTIONAL, INTENT(out)   :: albedo_nir(kdim)          !! Albedo of the NIR range
    REAL(dp), OPTIONAL, INTENT(out)   :: emis(kdim)                !! Emissivity
    REAL(dp), OPTIONAL, INTENT(out)   :: z0h(kdim)                 !! Surface roughness for heat
    REAL(dp), OPTIONAL, INTENT(out)   :: z0m(kdim)                 !! Surface roughness for momentum
    REAL(dp), OPTIONAL, INTENT(out)   :: runoff(kdim)              !! Surface runoff for g3b
    REAL(dp), OPTIONAL, INTENT(out)   :: drainage(kdim)            !! Surface drainage for g3b
    REAL(dp), OPTIONAL, INTENT(out)   :: skin_res(kdim)            !! Skin reservoir for g3b
    REAL(dp), OPTIONAL, INTENT(out)   :: tte_corr(kdim)            !! correction of tte without zdp (!)
    REAL(dp), OPTIONAL, INTENT(out)   :: snow_depth(kdim)          !! snow on soil acc. to sn in g3b
    REAL(dp), OPTIONAL, INTENT(out)   :: soil_wetness(kdim)        !! Soil moisture acc. to pws in echam for diagn - g3b
    REAL(dp), OPTIONAL, INTENT(out)   :: cair(kdim), csat(kdim)    !!
    REAL(dp), OPTIONAL, INTENT(out)   :: zhsoil(kdim)              !!
    REAL(dp), OPTIONAL, INTENT(out)   :: surf_dry_static_energy(kdim) !! qslnew in vdiff (zqsl*radiative_temp)
    REAL(dp), OPTIONAL, INTENT(in)    :: echam_zchl(kdim)          !!
    REAL(dp), OPTIONAL, INTENT(out)   :: glac_runoff_evap(kdim)    !! For HD-Model in ocean coupling
    REAL(dp), OPTIONAL, INTENT(out)   :: surf_runoff_hd(kdim)      !! For HD-Model in ocean coupling
    REAL(dp), OPTIONAL, INTENT(out)   :: drainage_hd(kdim)         !! For HD-Model in ocean coupling
    ! The following variables are passed back to echam just to be put into the standard echam output streams.
    ! This should be revised once the echam and jsbach output is reconciled.
    REAL(dp), OPTIONAL, INTENT(out), DIMENSION(kdim)   :: &
         glacier_depth,            &
         snow_melt_acc,            &
         glacier_p_minus_e_acc,    &
         snow_acc,                 &
         glacier_runoff_acc,       &
         wsmax

!---wiso-code
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_evap_act(kdim,nwiso)            !! Total of evaporation (actual) [kg/(m^2 s)]
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_evap_pot(kdim,nwiso)            !! Potential evaporation [kg/(m^2 s)]
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_qsurf(kdim,nwiso)               !! Surface specific humidity
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_runoff(kdim,nwiso)              !! Surface runoff for g3b
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_drainage(kdim,nwiso)            !! Surface drainage for g3b
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_skin_res(kdim,nwiso)            !! Skin reservoir for g3b
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_snow_depth(kdim,nwiso)          !! snow on soil acc. to sn in g3b
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_soil_wetness(kdim,nwiso)        !! Soil moisture acc. to pws in echam for diagn-g3b
    REAL(dp), OPTIONAL, INTENT(out)   :: wisocair(kdim,nwiso), wisocsat(kdim,nwiso)
    REAL(dp), OPTIONAL, INTENT(out)   :: wisocair_fra(kdim,nwiso), wisocsat_fra(kdim,nwiso)
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_eqAcoef_new(kdim,nwiso), wiso_eqBcoef_new(kdim,nwiso)
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_glac_runoff_evap(kdim,nwiso)    !! For HD-Model in ocean coupling
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_surf_runoff_hd(kdim,nwiso)      !! For HD-Model in ocean coupling
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_drainage_hd(kdim,nwiso)         !! For HD-Model in ocean coupling

    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_glacier_depth(kdim,nwiso)
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_snow_melt_acc(kdim,nwiso)
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_glacier_p_minus_e_acc(kdim,nwiso)
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_snow_acc(kdim,nwiso)
    REAL(dp), OPTIONAL, INTENT(out)   :: wiso_glacier_runoff_acc(kdim,nwiso)
!---wiso-code-end

    !! Local declarations for packed (gathered) input fields
    !LOGICAL, DIMENSION(kland) ::   zmask_glacier                       !! Packed mask for glaciers
    !LOGICAL, DIMENSION(kland) ::   zmask_land                          !! Packed mask for lakes
    REAL(dp), DIMENSION(kland) ::   zgeopot, zwind, zwind10, ztemp_air, zqair
    REAL(dp), DIMENSION(kland) ::   zetAcoef, zetBcoef, zeqAcoef, zeqBcoef
    REAL(dp), DIMENSION(kland) ::   zcdrag
    REAL(dp), DIMENSION(kland) ::   zprecip_rain, zprecip_snow         !! Packed precipitation
    REAL(dp), DIMENSION(kland) ::   zCO2                               !! Packed CO2 concentration
    REAL(dp), DIMENSION(kland) ::   zlwdown                            !! Packed longwave radiation
    REAL(dp), DIMENSION(kland) ::   zsw_vis_net                        !! Packed radiation from visible band [W/m^2]
    REAL(dp), DIMENSION(kland) ::   zsw_vis_frac_diffuse               !! Packed fraction of diffuse visible radiation
    REAL(dp), DIMENSION(kland) ::   zsw_nir_net                        !! Packed radiation from NIR band [W/m^2]
    REAL(dp), DIMENSION(kland) ::   zsw_nir_frac_diffuse               !! Packed fraction of diffuse NIR radiation
    REAL(dp), DIMENSION(kland) ::   zczenith                           !! Packed cosine of solar zenith angle
    REAL(dp), DIMENSION(kland) ::   zpressure                          !! Packed surface pressure
    REAL(dp), DIMENSION(kland) ::   zevap_act, zevap_pot               !! Packed actual evapotranspiration and potential evaporation
    REAL(dp), DIMENSION(kland) ::   zsensible, zlatent                 !! Packed surface fluxes
    REAL(dp), DIMENSION(kland) ::   zCO2_flux                          !! CO2 flux to the atmosphere [kg(CO2)/m^2 s]
    REAL(dp), DIMENSION(kland) ::   zCO2_flux_landcover_change         !! CO2 flux from landcover change     
    REAL(dp), DIMENSION(kland) ::   zCO2_flux_dynveg                   !! CO2 flux due to fires in the dynamic vegetation module   
    REAL(dp), DIMENSION(kland) ::   ztsoil_rad, ztemp_soil_new, zqsurf
    REAL(dp), DIMENSION(kland) ::   zalbedo                            !! Packed albedo
    REAL(dp), DIMENSION(kland) ::   zalbedo_vis                        !! Packed albedo of visible range
    REAL(dp), DIMENSION(kland) ::   zalbedo_nir                        !! Packed albedo of NIR range
    REAL(dp), DIMENSION(kland) ::   zemis                              !! Packed emissivity
    REAL(dp), DIMENSION(kland) ::   zz0h                               !! Packed surface roughness for heat
    REAL(dp), DIMENSION(kland) ::   zz0m                               !! Packed surface roughness for momentum
    REAL(dp), DIMENSION(kland) ::   zrunoff                            !! Packed surface runoff
    REAL(dp), DIMENSION(kland) ::   zdrainage                          !! Packed surface drainage
    REAL(dp), DIMENSION(kland) ::   zsnow_depth                        !! Packed soil snow
    REAL(dp), DIMENSION(kland) ::   zsoil_wetness                      !! PAcked soil moisture (acc. to pws in old echam)
    REAL(dp), DIMENSION(kland) ::   zcair, zcsat
    REAL(dp), DIMENSION(kland) ::   p_echam_zchl
    REAL(dp), DIMENSION(kland) ::   zskin_res

    !! Other local declarations
    LOGICAL,  DIMENSION(kdim)  ::   mask                               !! Land mask or all true's
    REAL(dp), DIMENSION(kland) ::   swdown                             !! Total shortwave radiation (visible + near infrared) [W/m^2]
    REAL(dp), DIMENSION(kland) ::   swnet                              !! Total net shortwave radiation [W/m^2]
    REAL(dp), DIMENSION(kland) ::   sw_par                             !! Radiation in PAR-band [W/m^2]
    REAL(dp), DIMENSION(kland) ::   sw_par_fract_direct                !! fraction of direct radiation in PAR-band
    REAL(dp), DIMENSION(kland) ::   zzhsoil                            !! relative hum of land surface acc to zhsoil in vdiff
    REAL(dp), DIMENSION(kland) ::   ztte_corr                          !! correction of tte multiplied by zqp (!)
    REAL(dp), DIMENSION(kland) ::   zdry_static_energy                 !! new surface dry static energy
    INTEGER,  DIMENSION(kland,theLand%ntiles) :: pheno_type            !! Mapping of land cover types to phenology types
    REAL(dp), DIMENSION(kland,theLand%ntiles) :: specificLeafArea_C 

    REAL(dp), DIMENSION(kland) :: zglac_runoff_evap
    REAL(dp), DIMENSION(kland) :: zsurf_runoff_hd
    REAL(dp), DIMENSION(kland) :: zdrainage_hd

!---wiso-code
    REAL(dp), DIMENSION(kland,nwiso) :: zwisokinl !! Packed kinetic fractionation factor
    REAL(dp), DIMENSION(kland)       :: zwiso_helpqdif

    REAL(dp), DIMENSION(kland,nwiso) :: zwisoeqAcoef, zwisoeqBcoef, zwisoqair
    REAL(dp), DIMENSION(kland,nwiso) :: zwiso_nenner
    REAL(dp), DIMENSION(kland,nwiso) :: zwisoqsurf
    REAL(dp), DIMENSION(kland,nwiso) :: zwisoprecip_rain, zwisoprecip_snow !! Packed precipitation
    REAL(dp), DIMENSION(kland,nwiso) :: zwisoevap_act, zwisoevap_pot !Packed actual evapotranspiration and potential evaporation
    REAL(dp), DIMENSION(kland,nwiso) :: zwisorunoff !! Packed surface runoff
    REAL(dp), DIMENSION(kland,nwiso) :: zwisodrainage !! Packed surface drainage
    REAL(dp), DIMENSION(kland,nwiso) :: zwisosnow_depth !! Packed soil snow
    REAL(dp), DIMENSION(kland,nwiso) :: zwisosoil_wetness !! Packed soil moisture (acc. to pws in old echam)
    REAL(dp), DIMENSION(kland,nwiso) :: zwisocair, zwisocsat
    REAL(dp), DIMENSION(kland,nwiso) :: zwisocair_fra, zwisocsat_fra
    REAL(dp), DIMENSION(kland,nwiso) :: zwiso_eqAcoef_new(kdim,nwiso), zwiso_eqBcoef_new(kdim,nwiso)
    REAL(dp), DIMENSION(kland,nwiso) :: zwisoskin_res

    REAL(dp), DIMENSION(kland,nwiso) :: zwisoglac_runoff_evap
    REAL(dp), DIMENSION(kland,nwiso) :: zwisosurf_runoff_hd
    REAL(dp), DIMENSION(kland,nwiso) :: zwisodrainage_hd

    REAL(dp), DIMENSION(kland,nwiso) :: zwisoglacierdepth, zwisosnowmeltacc, zwisoglacierpminuseacc, &
                                        zwisosnowacc, zwisoglacierrunoffacc  
    ! Variables for PACKing
    REAL(dp), DIMENSION(kdim)  :: help_wiso_eqAcoef
    REAL(dp), DIMENSION(kdim)  :: help_wiso_eqBcoef
    REAL(dp), DIMENSION(kdim)  :: help_wiso_nenner
    REAL(dp), DIMENSION(kdim)  :: help_wisoqair
    REAL(dp), DIMENSION(kdim)  :: help_wisoprecip_rain
    REAL(dp), DIMENSION(kdim)  :: help_wisoprecip_snow
    REAL(dp), DIMENSION(kdim)  :: help_wisokinl
    REAL(dp), DIMENSION(kland) :: help_zwisoeqAcoef
    REAL(dp), DIMENSION(kland) :: help_zwisoeqBcoef
    REAL(dp), DIMENSION(kland) :: help_zwiso_nenner
    REAL(dp), DIMENSION(kland) :: help_zwisoqair
    REAL(dp), DIMENSION(kland) :: help_zwisoprecip_rain
    REAL(dp), DIMENSION(kland) :: help_zwisoprecip_snow
    REAL(dp), DIMENSION(kland) :: help_zwisokinl

!---wiso-code-end

    LOGICAL :: l_inquire
    INTEGER :: ntiles, nsoil,isoil

    INTEGER                 :: i, kidx0, kidx1, itile, jt, jl, jk
    !DEBUG WISO
    REAL(dp)                :: zdelta

#if defined (__SX__) && defined (_OPENMP)
    INTEGER :: tid          ! OpenMP thread number
#endif

#if defined (__SX__) && defined (_OPENMP)
    IF (debug) THEN
       tid = omp_get_thread_num()
       CALL message('jsbach_interface', 'OpenMP thread #'//int2string(tid))
       CALL message('                ', '   kdim    = '//int2string(kdim))
       CALL message('                ', '   kland   = '//int2string(kland))
       CALL message('                ', '   kblock  = '//int2string(kblock))
    END IF
#endif

    IF (PRESENT(mask_land)) THEN
       IF (COUNT(mask_land(1:kdim)) /= kland) THEN
          CALL message('jsbach_interface ','kland, COUNT(mask_land) = '// &
                            int2string(kland)//', '//int2string(COUNT(mask_land(1:kdim))))
          CALL finish('jsbach_interface ','wrong mask_land')
       END IF
    END IF
    
    l_inquire = .FALSE.
    IF (PRESENT(mode)) THEN
       IF (TRIM(mode) == 'inquire' .OR. TRIM(mode) == 'enquire') THEN
          l_inquire = .TRUE.
       END IF
    END IF

    IF (kland == 0) THEN
       IF (PRESENT(evap_act)) evap_act                             = 0.0_dp
       IF (PRESENT(evap_pot)) evap_pot                             = 0.0_dp
       IF (PRESENT(sensible)) sensible                             = 0.0_dp
       IF (PRESENT(latent)) latent                                 = 0.0_dp
       IF (PRESENT(CO2_flux)) CO2_flux                             = 0.0_dp
!!$       IF (PRESENT(tsoil_rad)) tsoil_rad                           = 280.0_dp
!!$       IF (PRESENT(temp_soil_new)) temp_soil_new                   = 280.0_dp
       IF (PRESENT(tsoil_rad)) tsoil_rad                           = 0.0_dp
       IF (PRESENT(temp_soil_new)) temp_soil_new                   = 0.0_dp
       IF (PRESENT(qsurf)) qsurf                                   = 0.0_dp
       IF (PRESENT(albedo)) albedo                                 = 0.0_dp
       IF (PRESENT(albedo_vis)) albedo_vis                         = 0.0_dp
       IF (PRESENT(albedo_nir)) albedo_nir                         = 0.0_dp
       IF (PRESENT(emis)) emis                                     = 0.0_dp
       IF (PRESENT(z0h)) z0h                                       = 0.0_dp
       IF (PRESENT(z0m)) z0m                                       = 0.0_dp
       IF (PRESENT(cair)) cair                                     = 0.0_dp
       IF (PRESENT(csat)) csat                                     = 0.0_dp
       IF (PRESENT(zhsoil)) zhsoil                                 = 0.0_dp
!!$       IF (PRESENT(surf_dry_static_energy)) surf_dry_static_energy = 2.276e+05_dp
       IF (PRESENT(surf_dry_static_energy)) surf_dry_static_energy = 0.0_dp
       IF (PRESENT(soil_wetness)) soil_wetness                     = 0.0_dp
       IF (PRESENT(snow_depth)) snow_depth                         = 0.0_dp
       IF (PRESENT(runoff)) runoff                                 = 0.0_dp
       IF (PRESENT(drainage)) drainage                             = 0.0_dp
       IF (PRESENT(skin_res)) skin_res                             = 0.0_dp
       IF (PRESENT(tte_corr)) tte_corr                             = 0.0_dp
       IF (PRESENT(glac_runoff_evap ))glac_runoff_evap             = 0.0_dp
       IF (PRESENT(surf_runoff_hd )) surf_runoff_hd                = 0.0_dp
       IF (PRESENT(drainage_hd )) drainage_hd                      = 0.0_dp
       IF (PRESENT(glacier_depth)) glacier_depth                   = 0.0_dp
       IF (PRESENT(snow_melt_acc)) snow_melt_acc                   = 0.0_dp
       IF (PRESENT(glacier_p_minus_e_acc)) glacier_p_minus_e_acc   = 0.0_dp
       IF (PRESENT(snow_acc)) snow_acc                             = 0.0_dp
       IF (PRESENT(glacier_runoff_acc)) glacier_runoff_acc         = 0.0_dp
       IF (PRESENT(wsmax)) wsmax                                   = 0.0_dp

       RETURN
    END IF

!---wiso-code
    IF (lwiso) THEN

    IF (kland == 0) THEN
       IF (PRESENT(wiso_evap_act)) wiso_evap_act                   = 0.0_dp
       IF (PRESENT(wiso_evap_pot)) wiso_evap_pot                   = 0.0_dp
       IF (PRESENT(wiso_qsurf)) wiso_qsurf                         = 0.0_dp
       IF (PRESENT(wisocair)) wisocair                             = 0.0_dp
       IF (PRESENT(wisocsat)) wisocsat                             = 0.0_dp
       IF (PRESENT(wisocair_fra)) wisocair_fra                     = 0.0_dp
       IF (PRESENT(wisocsat_fra)) wisocsat_fra                     = 0.0_dp
       IF (PRESENT(wiso_eqAcoef_new)) wiso_eqAcoef_new             = 0.0_dp
       IF (PRESENT(wiso_eqBcoef_new)) wiso_eqBcoef_new             = 0.0_dp
       IF (PRESENT(wiso_soil_wetness)) wiso_soil_wetness           = 0.0_dp
       IF (PRESENT(wiso_snow_depth)) wiso_snow_depth               = 0.0_dp
       IF (PRESENT(wiso_runoff)) wiso_runoff                       = 0.0_dp
       IF (PRESENT(wiso_drainage)) wiso_drainage                   = 0.0_dp
       IF (PRESENT(wiso_skin_res)) wiso_skin_res                   = 0.0_dp
       IF (PRESENT(wiso_glac_runoff_evap )) wiso_glac_runoff_evap  = 0.0_dp
       IF (PRESENT(wiso_surf_runoff_hd )) wiso_surf_runoff_hd      = 0.0_dp
       IF (PRESENT(wiso_drainage_hd )) wiso_drainage_hd            = 0.0_dp
       IF (PRESENT(wiso_glacier_depth)) wiso_glacier_depth         = 0.0_dp
       IF (PRESENT(wiso_snow_melt_acc)) wiso_snow_melt_acc         = 0.0_dp
       IF (PRESENT(wiso_glacier_p_minus_e_acc)) wiso_glacier_p_minus_e_acc = 0.0_dp
       IF (PRESENT(wiso_snow_acc)) wiso_snow_acc                   = 0.0_dp
       IF (PRESENT(wiso_glacier_runoff_acc)) wiso_glacier_runoff_acc = 0.0_dp
       RETURN
    END IF
    
    END IF
!---wiso-code-end

    ! If first call do initializations (.g. post restart)
    IF ( jsbach_FirstCall ) THEN
       IF (debug) CALL message('jsbach_interface', ' First call to interface       - '//logical2string(jsbach_FirstCall))
       jsbach_FirstCall = .FALSE.
    ENDIF
    CALL update_current_call(kdim, kland, kblock=kblock, mask=mask_land)

    kidx0 = kstart
    kidx1 = kend

#if defined (__SX__) && defined (_OPENMP)
    IF (debug) THEN
       tid = omp_get_thread_num()
       CALL message('jsbach_interface', 'OpenMP thread #'//int2string(tid)//' '//int2string(kblock)//' '//int2string(kland)//&
            ' '//int2string(nidx)//' '//int2string(kidx0)//' '//int2string(kidx1))
    END IF
#endif

    IF (PRESENT(mask_land)) THEN
       mask = mask_land
    ELSE
       mask = .TRUE.
    ENDIF

    ntiles = theLand%ntiles
    nsoil = theLand%Soil%nsoil

    !Generate local packed forcing arrays for each domain (processor)
    !zmask_glacier = PACK(mask_glacier, MASK=mask)
    !zmask_lake   = PACK(mask_lake, MASK=mask)
    IF (PRESENT(geopot)) zgeopot      = PACK(geopot, MASK=mask)
    !zuwind       = PACK(uwind, MASK=mask)
    !zvwind       = PACK(vwind, MASK=mask)
    IF (PRESENT(wind)) zwind       = PACK(wind, MASK=mask)
    IF (PRESENT(wind10)) zwind10     = PACK(wind10, MASK=mask)
    IF (PRESENT(temp_air)) ztemp_air    = PACK(temp_air, MASK=mask)
    IF (PRESENT(qair)) zqair        = PACK(qair, MASK=mask)
    IF (PRESENT(cdrag))   zcdrag       = PACK(cdrag, MASK=mask)
    IF (PRESENT(etAcoef)) zetAcoef     = PACK(etAcoef, MASK=mask)
    IF (PRESENT(etBcoef)) zetBcoef     = PACK(etBcoef, MASK=mask)
    IF (PRESENT(eqAcoef)) zeqAcoef     = PACK(eqAcoef, MASK=mask)
    IF (PRESENT(eqBcoef)) zeqBcoef     = PACK(eqBcoef, MASK=mask)
    IF (PRESENT(precip_rain)) zprecip_rain = PACK(precip_rain, MASK=mask)
    IF (PRESENT(precip_snow)) zprecip_snow = PACK(precip_snow, MASK=mask)
    IF (PRESENT(lwdown)) zlwdown      = PACK(lwdown, MASK=mask)
!!$    zswnet       = PACK(swnet, MASK=mask)
!!$    zswdown      = PACK(swdown, MASK=mask)
    IF (PRESENT(sw_vis_net)) zsw_vis_net = PACK(sw_vis_net, MASK=mask)
    IF (PRESENT(sw_vis_frac_diffuse)) zsw_vis_frac_diffuse = PACK(sw_vis_frac_diffuse, MASK=mask)
    IF (PRESENT(sw_nir_net)) zsw_nir_net = PACK(sw_nir_net, MASK=mask)
    IF (PRESENT(sw_nir_frac_diffuse)) zsw_nir_frac_diffuse = PACK(sw_nir_frac_diffuse, MASK=mask)
    IF (PRESENT(czenith)) zczenith     = PACK(czenith, MASK=mask)
    IF (PRESENT(pressure)) zpressure    = PACK(pressure, MASK=mask)
    IF (PRESENT(CO2_concentration)) zCO2         = PACK(CO2_concentration, MASK=mask)
    IF (PRESENT(echam_zchl)) p_echam_zchl = PACK(echam_zchl, MASK=mask)

!---wiso-code
    IF (lwiso) THEN

    ! Generate local packed forcing arrays for each domain (processor) - Isotopes
    IF (PRESENT(wiso_helpqdif)) zwiso_helpqdif     = PACK(wiso_helpqdif, MASK=mask)
   
    IF (PRESENT(wiso_precip_rain)) THEN
      DO jt=1,nwiso
         help_wisoprecip_rain(:)    = wiso_precip_rain(:,jt)
         help_zwisoprecip_rain      = PACK(help_wisoprecip_rain, MASK=mask)
         zwisoprecip_rain(:,jt)     = help_zwisoprecip_rain(:)
      END DO
    END IF
    IF (PRESENT(wiso_precip_snow)) THEN
      DO jt=1,nwiso
         help_wisoprecip_snow(:)    = wiso_precip_snow(:,jt)
         help_zwisoprecip_snow      = PACK(help_wisoprecip_snow, MASK=mask)
         zwisoprecip_snow(:,jt)     = help_zwisoprecip_snow(:)
      END DO
    END IF
    IF (PRESENT(wiso_eqBcoef)) THEN
      DO jt=1,nwiso
         help_wiso_eqBcoef(:) = wiso_eqBcoef(:,jt)
         help_zwisoeqBcoef    = PACK(help_wiso_eqBcoef, MASK=mask)
         zwisoeqBcoef(:,jt)   = help_zwisoeqBcoef(:)
      END DO
    END IF 
    IF (PRESENT(wiso_eqAcoef)) THEN
      DO jt=1,nwiso
         help_wiso_eqAcoef(:) = wiso_eqAcoef(:,jt)
         help_zwisoeqAcoef    = PACK(help_wiso_eqAcoef, MASK=mask)
         zwisoeqAcoef(:,jt)   = help_zwisoeqAcoef(:)
      END DO
    END IF
    IF (PRESENT(wiso_nenner)) THEN
      DO jt=1,nwiso
         help_wiso_nenner(:)  = wiso_nenner(:,jt)
         help_zwiso_nenner    = PACK(help_wiso_nenner, MASK=mask)
         zwiso_nenner(:,jt)   = help_zwiso_nenner(:)
      END DO
    END IF 
    IF (PRESENT(wisoqair)) THEN
      DO jt=1,nwiso
         help_wisoqair(:)     = wisoqair(:,jt)
         help_zwisoqair       = PACK(help_wisoqair, MASK=mask)
         zwisoqair(:,jt)      = help_zwisoqair(:)
      END DO
    END IF 
    IF (PRESENT(wisokinl)) THEN
      DO jt=1,nwiso
         help_wisokinl(:)     = wisokinl(:,jt)
         help_zwisokinl       = PACK(help_wisokinl, MASK=mask)
         zwisokinl(:,jt)      = help_zwisokinl(:)
      END DO
    END IF
    
    END IF
!---wiso-code-end

    IF (lstart .OR. (lresume .AND. theOptions%ReadCoverFract)) THEN
       theLand%Surface%cover_fract(:,:) = init_cover_fract(:,:)
       IF (jsbach_NewTimeStep .AND. .NOT. lstart) CALL message('jsbach_interface','Using cover fractions from ini file')
       IF (theOptions%UseDynveg) THEN
         theLand%Surface%veg_ratio_max(:) = init_veg_ratio_max(:)
         IF (jsbach_NewTimeStep .AND. .NOT. lstart) CALL message('jsbach_interface','Using veg_ratio_max from ini file')
       END IF
    END IF

    !! --- Landcover change ---
    !! In case landcover change is driven by a sequence of external maps, read the appropriate map, replace the old one 
    !! (theLand%Surface%cover_fract), perform the necessary changes in the carbon pools and compute the CO2 emissions from
    !! the change of landcover 
    !! NOTE: landcover change has to be done before any other submodels of JSBACH are called to guarantee that during the
    !!       whole timestep the same landcover map is used! 
    !! NOTE: If the model is coupled to mpiom/hamocc, the prism setup starts the run one time step before the actual initial
    !!       date, usually the last time step in December of the previous year. Therefore, although do_landcover_change is
    !!       called, no land cover data are read (since current_date and previous_date belong to the same year). This is ok,
    !!       since for this first time step, the cover fractions from the ini file (in %cover_fract) can be used. However, it
    !!       seems to be important that the cover_fract in the ini file is the same as the one for the first year ... 
    !!       otherwise, the model crashes with negative CO2 concentrations, probably because the change in land cover is
    !!       too sudden, leading to strange effects in the CO2 transport.
    zCO2_flux_landcover_change(:) = 0._dp
    IF (theOptions%LCChange) THEN
       CALL do_landcover_change(theLand%Grid,theLand%Domain, theLand%ntiles,                           &
                                theLand%Surface%is_glacier(kidx0:kidx1,1),                             &
                                theLand%Cbalance,                                                      & !! <-- inout
                                theLand%Surface%veg_ratio_max(kidx0:kidx1),                            &
                                theLand%Vegetation%veg_fract_correction(kidx0:kidx1,1:theLand%ntiles), &
                                theLand%Surface%cover_fract(kidx0:kidx1,1:theLand%ntiles),             & !! <-- inout
                                zCO2_flux_landcover_change(1:nidx))                                      !! <-- out
    END IF

    DO itile=1,ntiles
       DO isoil=1,nsoil
          theLand%Soil%relative_moisture(kidx0:kidx1,isoil,itile) = &
               theLand%Soil%moisture(kidx0:kidx1,isoil,itile)/theLand%soilparam%MaxMoisture(kidx0:kidx1,1)
       END DO
       specificLeafArea_C(1:nidx,itile) = &
            theLand%LctLibrary%specificLeafArea_C(theLand%Surface%cover_type(kidx0:kidx1,itile))
    END DO

    IF (theOptions%UseDynveg .AND. .NOT. l_inquire) THEN
       IF (debug) CALL message('jsbach_interface_1d','Calling update_dynveg')
       CALL update_dynveg (lstart, lresume, nidx, kidx0, kidx1, &
           ztemp_air - tmelt, theLand%Soil%relative_humidity_air(kidx0:kidx1), &
           theLand%Soil%surface_temperature(kidx0:kidx1,:) - tmelt, zwind10, &
           zprecip_rain(:) + zprecip_snow(:), theLand%Soil%relative_moisture(kidx0:kidx1,1:nsoil,1:ntiles), & 
           theLand%cbalance%NPP_rate(kidx0:kidx1,1:ntiles), SpecificLeafArea_C(1:nidx,1:ntiles), &
           theLand%Vegetation%lai(kidx0:kidx1,1:ntiles), theLand%Vegetation%veg_fract_correction(kidx0:kidx1,1:ntiles), &
           theLand%Surface%is_glacier(kidx0:kidx1,1), theLand%Surface%is_present(kidx0:kidx1,1:ntiles), &
           theOptions%ReadCoverFract, &
           theLand%Cbalance%Cpool_green(kidx0:kidx1,1:ntiles), theLand%Cbalance%Cpool_reserve(kidx0:kidx1,1:ntiles), &
           theLand%Cbalance%Cpool_woods(kidx0:kidx1,1:ntiles), theLand%Cbalance%Cpool_litter_leaf(kidx0:kidx1,1:ntiles), &
           theLand%Cbalance%Cpool_litter_wood(kidx0:kidx1,1:ntiles), theLand%Cbalance%Cpool_fast(kidx0:kidx1,1:ntiles), &
           theLand%Cbalance%Cpool_slow(kidx0:kidx1,1:ntiles), theLand%Surface%rock_fract(kidx0:kidx1), &
           theLand%surface%cover_fract(kidx0:kidx1,1:ntiles), theLand%Surface%veg_ratio_max(kidx0:kidx1), &
           zCO2_flux_dynveg(:))
       IF (debug) CALL message('jsbach_interface_1d','Returned from update_dynveg')
    ENDIF

!!$    CALL message('jsbach_interface', 'Radiation Trigger l_trigrad      -'//logical2string(l_trigrad))
!!$    CALL message('jsbach_interface', 'Radiation Trigger l_trigradm1      -'//logical2string(l_trigradm1))


    IF (theOptions%UsePhenology) THEN
       DO itile=1,ntiles
          WHERE (theLand%Surface%is_vegetation(kidx0:kidx1,itile))
             pheno_type(1:nidx,itile) = theLand%LctLibrary%PhenologyType(theLand%Surface%cover_type(kidx0:kidx1,itile))
          ELSEWHERE
             pheno_type(1:nidx,itile) = 0
          END WHERE
       END DO
    END IF

    ! If first time step initialize lai with climatological value

    IF (lstart) THEN
#if defined (__SX__) && defined (_OPENMP)
       IF (debug) CALL message('jsbach_interface','Calling update_lai for start of model; thread '//int2string(tid))
#else
       IF (debug) CALL message('jsbach_interface','Calling update_lai for start of model')
#endif

       CALL update_lai(nidx, ntiles, theLand%Vegetation%lai_clim(kidx0:kidx1,1:ntiles,0:13), &
                                     theLand%Vegetation%lai     (kidx0:kidx1,1:ntiles))
       IF (theOptions%UsePhenology) THEN
          DO i=1,ntiles
             theLand%Vegetation%lai_logrop(kidx0:kidx1,i) = &
                  MIN(theLand%Vegetation%lai(kidx0:kidx1,i),theLand%Vegetation%lai_max(kidx0:kidx1,i))
          END DO
          theLand%Vegetation%lai(kidx0:kidx1,:) = theLand%Vegetation%lai_logrop(kidx0:kidx1,:)
          DO i=1,ntiles
             theLand%Surface%veg_ratio_actual_per_tile(kidx0:kidx1,i) = &
                  theLand%Surface%veg_ratio_max(kidx0:kidx1) * &
                  (1._dp - exp(-0.5_dp * theLand%Vegetation%lai(kidx0:kidx1,i)))
          END DO
       ELSE
          CALL update_cover_fract(nidx, ntiles, theLand%Surface%ntiles_lct, &
                  theLand%Surface%veg_ratio(kidx0:kidx1,0:13), &
                  theLand%Surface%veg_ratio_actual_per_tile(kidx0:kidx1,1:ntiles))
       END IF

    END IF

    ! Calculate albedo

    IF (.NOT. theOptions%UseAlbedo) THEN
      IF (l_inquire .OR. lstart) THEN
        CALL update_land_surface_fast(nidx, theLand%LctLibrary, theLand%Surface, &
             theLand%Soil%surface_temperature(kidx0:kidx1,:),                    &
             theLand%Soil%snow_fract         (kidx0:kidx1,:),                    &
             theLand%Soil%albedo             (kidx0:kidx1,:),                    &
             theLand%Vegetation%snow_fract_canopy(kidx0:kidx1,:),                &
             theLand%Vegetation%lai          (kidx0:kidx1,:))
        theLand%Surface%albedo_vis(kidx0:kidx1,:) = theLand%Surface%albedo(kidx0:kidx1,:)
        theLand%Surface%albedo_nir(kidx0:kidx1,:) = theLand%Surface%albedo(kidx0:kidx1,:)
      END IF
    ELSE
      IF (theOptions%Standalone .OR. l_inquire .OR. l_trigrad) THEN
#if defined (__SX__) && defined (_OPENMP)
        IF (debug) CALL message('jsbach_interface', 'update_albedo called l_inquire or l_trigrad; thread '//int2string(tid))
#else
        IF (debug) CALL message('jsbach_interface', 'update_albedo called l_inquire or l_trigrad')
#endif
        CALL update_albedo(theLand%LctLibrary,                    &
             l_inquire, l_trigrad, l_trigradm1, lstart,           &
             theOptions%Standalone,                               &
             nidx,                                                &
             theLand%Surface%ntiles,                              &
             theLand%Surface%cover_type          (kidx0:kidx1,:), &
             theLand%Surface%is_glacier          (kidx0:kidx1,:), &
             theLand%Surface%is_forest           (kidx0:kidx1,:), &
             theLand%Surface%veg_ratio_max       (kidx0:kidx1),   &
             zsw_vis_net                         (1:nidx),        &
             zsw_nir_net                         (1:nidx),        &
             theLand%Soil%surface_temperature    (kidx0:kidx1,:), &
             theLand%Soil%snow_fract             (kidx0:kidx1,:), &
             theLand%Soil%albedo_soil_vis        (kidx0:kidx1,:), &
             theLand%Soil%albedo_soil_nir        (kidx0:kidx1,:), &
             theLand%Soil%albedo_vegetation_vis  (kidx0:kidx1,:), &
             theLand%Soil%albedo_vegetation_nir  (kidx0:kidx1,:), &
             theLand%Vegetation%lai              (kidx0:kidx1,:), &
             theLand%Vegetation%snow_fract_canopy(kidx0:kidx1,:), &
             theLand%Surface%albedo_vis          (kidx0:kidx1,:), &
             theLand%Surface%albedo_nir          (kidx0:kidx1,:), &
             theLand%Surface%albedo              (kidx0:kidx1,:))
      END IF
    END IF

    ! Diagnose weighted albedo
    CALL average_tiles(theLand%Surface%albedo(kidx0:kidx1,1:ntiles), theLand%Surface%is_present(kidx0:kidx1,1:ntiles), &
                       theLand%Surface%cover_fract(kidx0:kidx1,1:ntiles), zalbedo(:))
    CALL average_tiles(theLand%Surface%albedo_vis(kidx0:kidx1,:), theLand%Surface%is_present(kidx0:kidx1,:), &
                       theLand%Surface%cover_fract(kidx0:kidx1,:), zalbedo_vis(:))
    CALL average_tiles(theLand%Surface%albedo_nir(kidx0:kidx1,:), theLand%Surface%is_present(kidx0:kidx1,:), &
                       theLand%Surface%cover_fract(kidx0:kidx1,:), zalbedo_nir(:))

    IF (l_inquire) THEN
       IF (kblock == 1) CALL message('jsbach_interface_1d','Inquiry mode')
       IF (PRESENT(albedo)) albedo         = UNPACK(zalbedo,mask,0.0_dp)
       IF (PRESENT(albedo_vis)) albedo_vis = UNPACK(zalbedo_vis,mask,0.0_dp)
       IF (PRESENT(albedo_nir)) albedo_nir = UNPACK(zalbedo_nir,mask,0.0_dp)
       ! kindex is allocated at the beginning of the interface through call to update_current_call
       DEALLOCATE(kindex)
       RETURN
    END IF

    ! Determine total solar shortwave radiation

    swnet(:) = zsw_vis_net(:) + zsw_nir_net(:)

    ! Determine downward shortwave radiation

    swdown(:) =  swnet(:) / (1._dp - zalbedo(:))

    ! Compute incoming photosynthetically active radiation (par) from radiation in visible band

    sw_par(:) =  zsw_vis_net(:) / (1._dp - zalbedo_vis(:))

    ! Compute fraction of direct part in PAR 

    sw_par_fract_direct(:) = 1._dp - zsw_vis_frac_diffuse(:)

    theLand%Surface%swdown_acc(kidx0:kidx1) = theLand%Surface%swdown_acc(kidx0:kidx1) + &
                                              (zsw_vis_net(:) / (1.0_dp - zalbedo_vis(:)) + &
                                              zsw_nir_net(:) / (1.0_dp - zalbedo_nir(:))) * delta_time
    theLand%Surface%swdown_reflect_acc(kidx0:kidx1) = theLand%Surface%swdown_reflect_acc(kidx0:kidx1) + &
                                              (zalbedo_vis(:) * zsw_vis_net(:) / (1.0_dp - zalbedo_vis(:)) + &
                                              zalbedo_nir(:) * zsw_nir_net(:) / (1.0_dp - zalbedo_nir(:))) * delta_time

    zCO2_flux = 0._dp

    
    ! Fist call to BETHY model (water unlimited case): 
    ! calculate unstressed stomatal conductance at prescribed leaf internal CO2 concentration
    IF (theOptions%UseBethy) THEN
       
       IF (debug) CALL message('jsbach_interface_1d','Calling bethy (first time)')

       CALL update_bethy(nidx, theLand%Domain, theLand%Surface%is_vegetation(kidx0:kidx1,1:ntiles), & 
            theLand%Bethy, &
            theLand%Surface%cover_type(kidx0:kidx1,:), &
            theLand%LctLibrary, &
            theOptions%UseAlbedo, &                                    !! Flag (true if interactive albedo is used)
            .FALSE.,                                               &   !! Flag for water limitation in photosynthesis
            zczenith             (:),                              &   !! Cosine of solar zenith angle
            declination             ,                              &   !! Solar declination
            theLand%Vegetation%lai(kidx0:kidx1,1:ntiles),          &   !! Leaf area index
            swdown               (:),                              &   !! Total downward shortwave radiation
            sw_par               (:),                              &   !! incoming shortwave radiation from PAR-band
            sw_par_fract_direct  (:),                              &   !! fraction of direct radiation in PAR-band
            zpressure            (:),                              &   !! Surface pressure
            ztemp_air            (:),                              &   !! canopy temperature (= air temperature ?)
            theLand%Soil%albedo(kidx0:kidx1,1:ntiles),             &   !! Soil albedo
            zCO2                 (:),                              &   !! CO2 concentration air
            theLand%Vegetation%canopy_conductance_bethy(kidx0:kidx1,1:ntiles), &   !! Unlimited canopy conductance (output)
            theLand%Vegetation%canopy_conductance_limited(kidx0:kidx1,1:ntiles) &   !! dummy
            )

    ELSE

       ! Compute canopy resistance without water stress

       IF (debug .AND. theOptions%UseEchamLand) CALL message('jsbach_interface_1d','Calling unstressed_canopy_cond_par')

       DO itile=1,ntiles
          IF (theOptions%UseEchamLand) THEN
             CALL unstressed_canopy_cond_par(theLand%Vegetation%lai(kidx0:kidx1,itile), zsw_vis_net(:), &
                  theLand%Vegetation%canopy_conductance_bethy(kidx0:kidx1,itile))
!!$          ELSE IF (theOptions%UseVic) THEN
!!$             CALL unstressed_canopy_resistance_rmin( &
!!$                  theLand%Vegetation% lai                (kidx0:kidx1,itile),                             &
!!$                  theLand%LctLibrary% CanopyResistanceMin(theLand%Surface%cover_type(kidx0:kidx1,itile)), &
!!$                  theLand%Vegetation% canopy_resistance  (kidx0:kidx1,itile))
          END IF
       END DO

       theLand%Vegetation%canopy_conductance_bethy(kidx0:kidx1,1:ntiles) = 1.e-5_dp   ! XXX

    END IF
    
    ! Note: theLand%Vegetation%canopy_resistance now contains the canopy resistance for the case that there's no water stress.
    ! It will be used to compute potential evapotranspiration (but not in the ECHAM formulation). Afterwards, it will be 
    ! corrected for water stress and used for the computation of actual evapotranspiration from vegetation.

    ! --- Run the land surface scheme ---------------------------------------------------------------------------------------

  IF (lwiso) THEN
    CALL update_soil(                                                                               &
                     ! kidx, surface, soil, soil_param, useDynveg
                     nidx, theLand%Surface, theLand%Soil, theLand%SoilParam, theOptions%UseDynveg,  &
                     ! canopy_conductance_max
                     theLand%Vegetation%canopy_conductance_bethy(kidx0:kidx1,:),                    &
                     ! root_depth
                     theLand%Vegetation%root_depth(kidx0:kidx1,:,:),                                &
                     ! lai
                     theLand%Vegetation%lai(kidx0:kidx1,:),                                         &
                     ! cancopy_snow
                     theLand%Vegetation%snow_depth_canopy(kidx0:kidx1,:),                           &
                     ! cancopy_snow_fract
                     theLand%Vegetation%snow_fract_canopy(kidx0:kidx1,:),                           &
                     ! cdrag, t_Acoef, t_Bcoef, q_Acoef, q_Bcoef, air_temperature, air_moisture
                     zcdrag, zetAcoef, zetBcoef, zeqAcoef, zeqBcoef, ztemp_air, zqair,              &
                     ! surface_pressure, windspeed, wind10
                     zpressure, zwind, zwind10,                                                     &
                     ! rad_longwave_down, rad_shortwave_net, rain, snow, cair, csat
                     zlwdown, swnet, zprecip_rain, zprecip_snow, zcair, zcsat,                      &
                     ! p_echam_zchl, zhsoil_avg, tte_corr_avg
                     p_echam_zchl, zzhsoil, ztte_corr,                                              &
                     ! canopy_conductance_limited
                     theLand%Vegetation%canopy_conductance_limited(kidx0:kidx1,:),                  &
                     ! glac_runoff_evap, surf_runoff_hd, drainage_hd
                     zglac_runoff_evap, zsurf_runoff_hd, zdrainage_hd,                              &
!---wiso-code
                     lwisofracl,                                                                    &
                     ! wiso_canopy_snow
                     theLand%Vegetation%wiso_snow_depth_canopy(kidx0:kidx1,:,:),                    &
                     ! q_wiso_Acoef, q_wiso_Bcoef, wiso_air_moisture, wisokinl
                     zwisoeqAcoef, zwisoeqBcoef, zwisoqair, zwisokinl,                              &
                     ! wiso_nenner
                     zwiso_nenner,                                                                  &
                     ! wiso_helpqdif
                     zwiso_helpqdif,                                                                &
                     ! wiso_rain, wiso_snow
                     zwisoprecip_rain, zwisoprecip_snow,                                            &
                     ! wiso_cair, wiso_csat, wisocair_fra, wisocsat_fra
                     zwisocair, zwisocsat, zwisocair_fra, zwisocsat_fra,                            &
                     ! wisoqklevl
                     zwiso_eqAcoef_new, zwiso_eqBcoef_new,       		                            &
                     ! glac_runoff_evap, surf_runoff_hd, drainage_hd
                     zwisoglac_runoff_evap, zwisosurf_runoff_hd, zwisodrainage_hd)
!---wiso-code-end
  ELSE
    CALL update_soil(                                                                               &
                     ! kidx, surface, soil, soil_param, useDynveg
                     nidx, theLand%Surface, theLand%Soil, theLand%SoilParam, theOptions%UseDynveg,  &
                     ! canopy_conductance_max
                     theLand%Vegetation%canopy_conductance_bethy(kidx0:kidx1,:),                    &
                     ! root_depth
                     theLand%Vegetation%root_depth(kidx0:kidx1,:,:),                                &
                     ! lai
                     theLand%Vegetation%lai(kidx0:kidx1,:),                                         &
                     ! cancopy_snow
                     theLand%Vegetation%snow_depth_canopy(kidx0:kidx1,:),                           &
                     ! cancopy_snow_fract
                     theLand%Vegetation%snow_fract_canopy(kidx0:kidx1,:),                           &
                     ! cdrag, t_Acoef, t_Bcoef, q_Acoef, q_Bcoef, air_temperature, air_moisture
                     zcdrag, zetAcoef, zetBcoef, zeqAcoef, zeqBcoef, ztemp_air, zqair,              &
                     ! surface_pressure, windspeed, wind10
                     zpressure, zwind, zwind10,                                                     &
                     ! rad_longwave_down, rad_shortwave_net, rain, snow, cair, csat
                     zlwdown, swnet, zprecip_rain, zprecip_snow, zcair, zcsat,                      &
                     ! p_echam_zchl, zhsoil_avg, tte_corr_avg
                     p_echam_zchl, zzhsoil, ztte_corr,                                              &
                     ! canopy_conductance_limited
                     theLand%Vegetation%canopy_conductance_limited(kidx0:kidx1,:),                  &
                     ! glac_runoff_evap, surf_runoff_hd, drainage_hd
                     zglac_runoff_evap, zsurf_runoff_hd, zdrainage_hd)
  END IF

!!$    WHERE (theLand%Vegetation%canopy_resistance(kidx0:kidx1,:) > 0.9_dp*HUGE(1._dp)) &
!!$         theLand%Vegetation%canopy_resistance(kidx0:kidx1,:) = 1.E20_dp


    ! Second call to BETHY model (water limited case): 
    ! calculate net assimilation at prescribed unstressed stomatal conductance

   IF (theOptions%UseBethy) THEN
      IF (debug) CALL message('jsbach_interface_1d','Calling bethy (second time)')
      CALL update_bethy(nidx, theLand%Domain, theLand%Surface%is_vegetation(kidx0:kidx1,1:ntiles), &
            theLand%Bethy, &
            theLand%Surface%cover_type(kidx0:kidx1,:), &
            theLand%LctLibrary, &
            theOptions%UseAlbedo, &                                    !! Flag (true if interactive albedo is used)
            .TRUE.,                                                &   !! Flag for water limited photosynthesis
            zczenith             (:),                              &   !! Cosine of solar zenith angle
            declination             ,                              &   !! Solar declination
            theLand%Vegetation%lai(kidx0:kidx1,1:ntiles),          &   !! Leaf area index
            swdown               (:),                              &   !! Total downward shortwave radiation
            sw_par               (:),                              &   !! incoming shortwave radiation from PAR-band
            sw_par_fract_direct  (:),                              &   !! fraction of direct radiation in PAR-band
            zpressure            (:),                              &   !! Surface pressure
            ztemp_air            (:),                              &   !! canopy temperature (= air temperature ?)
            theLand%Soil%albedo(kidx0:kidx1,1:ntiles),             &   !! Soil albedo
            zCO2                 (:),                              &   !! CO2 concentration air
            theLand%Vegetation%canopy_conductance_bethy(kidx0:kidx1,1:ntiles), &   !! Canopy conductance (output, 
                                               !! should be the same as limited conductance input on following line)
            theLand%Vegetation%canopy_conductance_limited(kidx0:kidx1,1:ntiles) &   !! Limited canopy conductance (input)
            )
      
      CALL update_cbalance_bethy(theLand%Domain, &
                                 theLand%Cbalance, &
                                 theLand%LctLibrary, &
                                 theLand%Surface, &
                                 theLand%Bethy%gross_assimilation(kidx0:kidx1,1:ntiles), &  ! input: GPP rel. to ground area
                                 theLand%Bethy%dark_respiration(kidx0:kidx1,1:ntiles), &  ! input: dark respiration rel. to ground area
                                 theLand%Soil%soil_temperature(kidx0:kidx1,3,1),  & ! input (third layer, only one tile because temperature is identical on all tiles)
                                 theLand%Soil%relative_moisture(kidx0:kidx1,1,1:ntiles),                  & ! input
                                 theLand%Vegetation%lai(kidx0:kidx1,1:ntiles),                            & ! Leaf area index (input)
                                 theLand%Vegetation%veg_fract_correction(kidx0:kidx1,1:ntiles),           & ! correction factor 1-exp(-LAI_max/2)
                                 zCO2_flux(1:nidx) )                                                        ! output: [kg(CO2)/m^2(gridbox) s]

   ELSE
      zCO2_flux(1:nidx) = 0._dp
   END IF

   ! Update phenology, i.e. recompute LAI

   IF (theOptions%UsePhenology) THEN

      IF (debug) CALL message('jsbach_interface_1d','Calling update_phenology')

      CALL update_phenology(nidx, &
           pheno_type(:, 1:ntiles), &
           theLand%Vegetation%lai_logrop(kidx0:kidx1, 1:ntiles), &
           theLand%Vegetation%lai_max(kidx0:kidx1, 1:ntiles), &
           (ztemp_air(:)-tmelt), &
           theLand%Soil%relative_moisture(kidx0:kidx1,1:nsoil,1:ntiles), &
           theLand%cbalance%NPP_rate(kidx0:kidx1,1:ntiles), &
           specificLeafArea_C(1:nidx,1:ntiles), &
           theLand%Domain%lat(kidx0:kidx1))

      IF (debug) CALL message('jsbach_interface_1d','Returned from update_phenology')

      theLand%Vegetation%lai(kidx0:kidx1,:) = theLand%Vegetation%lai_logrop(kidx0:kidx1,:)
      DO i = 1,ntiles
         theLand%Surface%veg_ratio_actual_per_tile(kidx0:kidx1,i) = &  
              theLand%Surface%veg_ratio_max(kidx0:kidx1) * &
              (1._dp - EXP(-0.5_dp * theLand%Vegetation%lai(kidx0:kidx1,i)))
      END DO

   ELSE

      IF (debug) CALL message('jsbach_interface_1d','Calling update_lai')

   END IF

   WHERE (theLand%Vegetation%lai_logrop(kidx0:kidx1,1:ntiles) < EPSILON(1._dp)) &
        theLand%Vegetation%lai_logrop(kidx0:kidx1,1:ntiles) = 0._dp
   WHERE (theLand%Vegetation%lai(kidx0:kidx1,1:ntiles) < EPSILON(1._dp)) &
        theLand%Vegetation%lai(kidx0:kidx1,1:ntiles) = 0._dp

   DO i=1,ntiles
      WHERE (theLand%Vegetation%lai_logrop(kidx0:kidx1,i) < 1.e-30_dp .AND. &
           theLand%Vegetation%lai_logrop(kidx0:kidx1,i) > 0._dp) &
           theLand%Vegetation%lai_logrop(kidx0:kidx1,i) = 1.e-30_dp
      WHERE (theLand%Vegetation%lai(kidx0:kidx1,i) < 1.e-30_dp .AND. theLand%Vegetation%lai(kidx0:kidx1,i) > 0._dp) &
           theLand%Vegetation%lai(kidx0:kidx1,i) = 1.e-30_dp
   END do

   ! Calculate albedo

   IF (.NOT. theOptions%UseAlbedo) THEN
      CALL update_land_surface_fast(nidx, theLand%LctLibrary, theLand%Surface, &
           theLand%Soil%surface_temperature(kidx0:kidx1,:),                    &
           theLand%Soil%snow_fract         (kidx0:kidx1,:),                    &
           theLand%Soil%albedo             (kidx0:kidx1,:),                    &
           theLand%Vegetation%snow_fract_canopy(kidx0:kidx1,:),                &
           theLand%Vegetation%lai          (kidx0:kidx1,:))
      theLand%Surface%albedo_vis(kidx0:kidx1,:) = theLand%Surface%albedo(kidx0:kidx1,:)
      theLand%Surface%albedo_nir(kidx0:kidx1,:) = theLand%Surface%albedo(kidx0:kidx1,:)
   ELSE
      IF (l_trigradm1 .AND. .NOT. lstart .AND. .NOT. theOptions%Standalone) THEN
         IF (debug) CALL message('jsbach_interface', 'update_albedo called l_trigradm1')
         CALL update_albedo(theLand%LctLibrary,                    &
              l_inquire, l_trigrad, l_trigradm1, lstart,           &
              theOptions%Standalone,                               &
              nidx,                                                &
              theLand%Surface%ntiles,                              &
              theLand%Surface%cover_type          (kidx0:kidx1,:), &
              theLand%Surface%is_glacier          (kidx0:kidx1,:), &
              theLand%Surface%is_forest           (kidx0:kidx1,:), &
              theLand%Surface%veg_ratio_max       (kidx0:kidx1),   &
              zsw_vis_net                         (1:nidx),        &
              zsw_nir_net                         (1:nidx),        &
              theLand%Soil%surface_temperature    (kidx0:kidx1,:), &
              theLand%Soil%snow_fract             (kidx0:kidx1,:), &
              theLand%Soil%albedo_soil_vis        (kidx0:kidx1,:), &
              theLand%Soil%albedo_soil_nir        (kidx0:kidx1,:), &
              theLand%Soil%albedo_vegetation_vis  (kidx0:kidx1,:), &
              theLand%Soil%albedo_vegetation_nir  (kidx0:kidx1,:), &
              theLand%Vegetation%lai              (kidx0:kidx1,:), &
              theLand%Vegetation%snow_fract_canopy(kidx0:kidx1,:), &
              theLand%Surface%albedo_vis          (kidx0:kidx1,:), &
              theLand%Surface%albedo_nir          (kidx0:kidx1,:), &
              theLand%Surface%albedo              (kidx0:kidx1,:))
      END IF
      CALL average_tiles(theLand%Surface%albedo_vis(kidx0:kidx1,:), theLand%Surface%is_present(kidx0:kidx1,:), &
           theLand%Surface%cover_fract(kidx0:kidx1,:), zalbedo_vis(:))
      CALL average_tiles(theLand%Surface%albedo_nir(kidx0:kidx1,:), theLand%Surface%is_present(kidx0:kidx1,:), &
           theLand%Surface%cover_fract(kidx0:kidx1,:), zalbedo_nir(:))
   END IF

   ! Compute diagnostic output
!!$   CALL land_surface_diagnostics(theLand%Surface)
   CALL soil_diagnostics(theLand%Surface, theLand%Soil)
   CALL veg_diagnostics (theLand%Surface, theLand%Vegetation)
   IF (theOptions%UseBethy) CALL bethy_diagnostics(theLand%Surface, theLand%Bethy)

   ! Compute roughness length and grid box areal average values
  IF (lwiso) THEN
   CALL update_land_boundary_up(theLand%LctLibrary, theLand%Surface, theLand%Soil,        &
        theLand%SoilParam%Roughness(kidx0:kidx1), zz0h, zz0m,                             &
        zalbedo, ztemp_soil_new, zqsurf, zevap_act, zevap_pot,                            &
        zsensible, zlatent, ztsoil_rad, zdry_static_energy, zsoil_wetness, zsnow_depth,   &
        zrunoff, zdrainage, zskin_res,							                          &
!---wiso-code
        zwisoqsurf, zwisoevap_act, zwisoevap_pot,                                         &
        zwisosoil_wetness,  zwisosnow_depth,                                              &
        zwisorunoff,      zwisodrainage,      zwisoskin_res)
!---wiso-code-end
  ELSE
   CALL update_land_boundary_up(theLand%LctLibrary, theLand%Surface, theLand%Soil,        &
        theLand%SoilParam%Roughness(kidx0:kidx1), zz0h, zz0m,                             &
        zalbedo, ztemp_soil_new, zqsurf, zevap_act, zevap_pot,                            &
        zsensible, zlatent, ztsoil_rad, zdry_static_energy, zsoil_wetness, zsnow_depth,   &
        zrunoff, zdrainage, zskin_res)
  END IF

   !! compute total CO2-flux by adding co2-flux from landcover change and fires
   IF (theOptions%LCChange) THEN
      zCO2_flux(:) = zCO2_flux(:) + zCO2_flux_landcover_change(:)
   END IF
   IF (theOptions%UseDynveg) THEN
      zCO2_flux(:) = zCO2_flux(:) + zCO2_flux_dynveg(:)
   END IF

   ! Unpack output fields
   IF (PRESENT(evap_act)) evap_act                             = UNPACK(zevap_act,mask,0.0_dp)
   IF (PRESENT(evap_pot)) evap_pot                             = UNPACK(zevap_pot,mask,0.0_dp)
   IF (PRESENT(sensible)) sensible                             = UNPACK(zsensible,mask,0.0_dp)
   IF (PRESENT(latent)) latent                                 = UNPACK(zlatent,mask,0.0_dp)
   IF (PRESENT(CO2_flux)) CO2_flux                             = UNPACK(zCO2_flux, mask, 0.0_dp)
   IF (PRESENT(tsoil_rad)) tsoil_rad                           = UNPACK(ztsoil_rad,mask,0.0_dp)
   IF (PRESENT(temp_soil_new)) temp_soil_new                   = UNPACK(ztemp_soil_new,mask,0.0_dp)
   IF (PRESENT(qsurf)) qsurf                                   = UNPACK(zqsurf,mask,0.0_dp)
   IF (PRESENT(albedo)) albedo                                 = UNPACK(zalbedo,mask,0.0_dp)
   IF (PRESENT(albedo_vis)) albedo_vis                         = UNPACK(zalbedo_vis,mask,0.0_dp)
   IF (PRESENT(albedo_nir)) albedo_nir                         = UNPACK(zalbedo_nir,mask,0.0_dp)
   IF (PRESENT(emis)) emis                                     = UNPACK(zemis,mask,0.0_dp)
   IF (PRESENT(z0h)) z0h                                       = UNPACK(zz0h,mask,0.0_dp)
   IF (PRESENT(z0m)) z0m                                       = UNPACK(zz0m,mask,0.0_dp)
   IF (PRESENT(cair)) cair                                     = UNPACK(zcair,mask,0.0_dp)
   IF (PRESENT(csat)) csat                                     = UNPACK(zcsat,mask,0.0_dp)
   IF (PRESENT(zhsoil)) zhsoil                                 = UNPACK(zzhsoil,mask,0.0_dp)
   IF (PRESENT(tsoil_rad)) tsoil_rad                           = UNPACK(ztsoil_rad,mask,0.0_dp)
   IF (PRESENT(surf_dry_static_energy)) surf_dry_static_energy = UNPACK(zdry_static_energy,mask,2.276e+05_dp)
   IF (PRESENT(soil_wetness)) soil_wetness                     = UNPACK(zsoil_wetness,mask,0.0_dp)
   IF (PRESENT(snow_depth)) snow_depth                         = UNPACK(zsnow_depth,mask,0.0_dp)
   IF (PRESENT(runoff)) runoff                                 = UNPACK(zrunoff,mask,0.0_dp)
   IF (PRESENT(drainage)) drainage                             = UNPACK(zdrainage,mask,0.0_dp)
   IF (PRESENT(skin_res)) skin_res                             = UNPACK(zskin_res,mask,0.0_dp)
   IF (PRESENT(tte_corr)) tte_corr                             = UNPACK(ztte_corr,mask,0.0_dp)

   IF (PRESENT(glac_runoff_evap ))glac_runoff_evap         = UNPACK(zglac_runoff_evap,mask,0.0_dp)
   IF (PRESENT(surf_runoff_hd )) surf_runoff_hd            = UNPACK(zsurf_runoff_hd,mask,0.0_dp)
   IF (PRESENT(drainage_hd )) drainage_hd                  = UNPACK(zdrainage_hd,mask,0.0_dp)
   IF (PRESENT(glacier_depth)) glacier_depth = UNPACK(get_soil_diag(nidx, 'glacier_depth'),mask,0.0_dp)
   IF (PRESENT(snow_melt_acc)) snow_melt_acc = UNPACK(get_soil_diag(nidx, 'snow_melt_acc'),mask,0.0_dp)
   IF (PRESENT(glacier_p_minus_e_acc)) glacier_p_minus_e_acc = &
                                       UNPACK(get_soil_diag(nidx, 'glacier_precip_minus_evap_acc'),mask,0.0_dp)
   IF (PRESENT(snow_acc)) snow_acc = UNPACK(get_soil_diag(nidx, 'snow_acc'),mask,0.0_dp)
   IF (PRESENT(wsmax)) wsmax = UNPACK(SUM(theLand%SoilParam%MaxMoisture(kidx0:kidx1,1:nsoil),DIM=2),mask,0.0_dp)
   IF (PRESENT(glacier_runoff_acc)) glacier_runoff_acc = &
                                       UNPACK(get_soil_diag(nidx, 'glacier_runoff_acc'),mask,0.0_dp)

!---wiso-code
    IF (lwiso) THEN

   zwisoglacierdepth(1:nidx,1:nwiso) = get_soil_diag_wiso(nidx,nwiso,'wiso_glacier_depth')
   zwisosnowmeltacc(1:nidx,1:nwiso) = get_soil_diag_wiso(nidx,nwiso,'wiso_snow_melt_acc')
   zwisoglacierpminuseacc(1:nidx,1:nwiso) = get_soil_diag_wiso(nidx,nwiso,'wiso_glac_p_minus_e_acc')
   zwisosnowacc(1:nidx,1:nwiso) = get_soil_diag_wiso(nidx,nwiso,'wiso_snow_acc')
   zwisoglacierrunoffacc(1:nidx,1:nwiso) = get_soil_diag_wiso(nidx,nwiso,'wiso_glacier_runoff_acc')

   ! Unpack outputfiles - Isotopes
   DO jt=1,nwiso
     IF (PRESENT(wiso_evap_act)) wiso_evap_act(:,jt)                  = UNPACK(zwisoevap_act(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_evap_pot)) wiso_evap_pot(:,jt)                  = UNPACK(zwisoevap_pot(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_qsurf)) wiso_qsurf(:,jt)                        = UNPACK(zwisoqsurf(:,jt),mask,0.0_dp)
     IF (PRESENT(wisocair)) wisocair(:,jt)                            = UNPACK(zwisocair(:,jt),mask,0.0_dp)
     IF (PRESENT(wisocsat)) wisocsat(:,jt)                            = UNPACK(zwisocsat(:,jt),mask,0.0_dp)
     IF (PRESENT(wisocair_fra)) wisocair_fra(:,jt)                    = UNPACK(zwisocair_fra(:,jt),mask,0.0_dp)
     IF (PRESENT(wisocsat_fra)) wisocsat_fra(:,jt)                    = UNPACK(zwisocsat_fra(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_eqAcoef_new)) wiso_eqAcoef_new(:,jt)            = UNPACK(zwiso_eqAcoef_new(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_eqBcoef_new)) wiso_eqBcoef_new(:,jt)            = UNPACK(zwiso_eqBcoef_new(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_soil_wetness)) wiso_soil_wetness(:,jt)          = UNPACK(zwisosoil_wetness(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_snow_depth)) wiso_snow_depth(:,jt)              = UNPACK(zwisosnow_depth(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_runoff)) wiso_runoff(:,jt)                      = UNPACK(zwisorunoff(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_drainage)) wiso_drainage(:,jt)                  = UNPACK(zwisodrainage(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_skin_res)) wiso_skin_res(:,jt)                  = UNPACK(zwisoskin_res(:,jt),mask,0.0_dp)

     IF (PRESENT(wiso_glac_runoff_evap)) wiso_glac_runoff_evap(:,jt)  = UNPACK(zwisoglac_runoff_evap(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_surf_runoff_hd)) wiso_surf_runoff_hd(:,jt)      = UNPACK(zwisosurf_runoff_hd(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_drainage_hd)) wiso_drainage_hd(:,jt)            = UNPACK(zwisodrainage_hd(:,jt),mask,0.0_dp)

     IF (PRESENT(wiso_glacier_depth)) wiso_glacier_depth(:,jt)        = UNPACK(zwisoglacierdepth(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_snow_melt_acc)) wiso_snow_melt_acc(:,jt)        = UNPACK(zwisosnowmeltacc(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_glacier_p_minus_e_acc)) wiso_glacier_p_minus_e_acc(:,jt) = UNPACK(zwisoglacierpminuseacc(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_snow_acc)) wiso_snow_acc(:,jt)                  = UNPACK(zwisosnowacc(:,jt),mask,0.0_dp)
     IF (PRESENT(wiso_glacier_runoff_acc)) wiso_glacier_runoff_acc(:,jt) = UNPACK(zwisoglacierrunoffacc(:,jt),mask,0.0_dp)
   END DO
   
   END IF
!---wiso-code-end

   ! kindex is allocated at the beginning of the interface through call to update_current_call
   DEALLOCATE(kindex)

 END SUBROUTINE jsbach_inter_1d

  SUBROUTINE jsbach_init(grid, domain, options, ntile)

    ! Initialization of the global grid, local PE domains, and all submodels of JSBACH

    USE mo_jsbach,        ONLY: jsbach_config
    USE mo_jsbach_grid,   ONLY: init_grid, init_domain
    USE mo_land_surface,  ONLY: init_land_surface
    USE mo_jsbach_lctlib, ONLY: init_lctlib
    USE mo_jsbach_veg,    ONLY: init_vegetation
    USE mo_soil,          ONLY: init_soil
    USE mo_phenology,     ONLY: init_phenology
    USE mo_time_control,  ONLY: lfirst_cycle, init_manager, init_times, construct_events, p_bcast_event, &
                                ec_manager_init, init_events, putdata, trigfiles, putrerun
    USE mo_time_event,    ONLY: io_time_event
    USE mo_linked_list,   ONLY: LAND, TILES
    USE mo_grib,          ONLY: init_grib, land_table
    USE mo_memory_base,   ONLY: new_stream, default_stream_setting
    USE mo_gaussgrid,     ONLY: inigau
    USE mo_lookup_tables, ONLY: init_lookup_table
    USE mo_convect_tables,ONLY: init_convect_tables
    USE mo_filename,      ONLY: GRIB


    USE mo_cbal_landcover_change, ONLY: init_landcover_change

!---wiso-code
    USE mo_wiso,          ONLY: lwiso
!---wiso-code-end

    TYPE(grid_type),    INTENT(out), OPTIONAL :: grid
    TYPE(domain_type),  INTENT(out), OPTIONAL :: domain
    TYPE(options_type), INTENT(out), OPTIONAL :: options
    INTEGER ,           INTENT(out), OPTIONAL :: ntile

    INTEGER :: IO_timestep, istep

    IF (debug) CALL message('jsbach_init','Begin')
    CALL jsbach_config(theOptions)

    CALL init_grid(theLand%Grid, &                                      ! Set global dimensions
                   theOptions%Standalone, lresume, theOptions%GridFile)

    IF (theOptions%Standalone) THEN
       CALL construct_events                                            ! Set up memory for event manager
       CALL get_dates(IO_timestep, istep)                               ! Get timestep
       CALL ec_manager_init(IO_timestep, istep)                         ! Initialize time manager
       CALL p_bcast_event(putdata, p_io)
       CALL p_bcast_event(trigfiles, p_io)
       CALL p_bcast_event(putrerun, p_io)
       IF (lfirst_cycle) CALL init_manager                              ! Re-initialize time manager (?)
       CALL init_times                                                  ! Initialize all dates
       CALL init_events                                                 ! Initialize event manager
       CALL jsbalone_init_decomp(theLand%Grid%nlon, theLand%Grid%nlat, & ! Initialize decomposition
            theLand%Grid%mask, theLand%Grid%nproca, theLand%Grid%nprocb, &
            theLand%Grid%npedim)
       CALL inigau                                                      ! Initialize gausgrid parameters 
       CALL jsbalone_iniphy                                             ! Physical params
       CALL init_lookup_table                                           ! lookup tables from echam
       CALL init_convect_tables                                         ! convection tables from echam
       IF (theOptions%FileType == GRIB) CALL init_grib                  ! Grib initialization
    ENDIF

    ! Add JSBACH I/O stream
    CALL new_stream(theLand%IO_stream, 'jsbach', filetype=theOptions%FileType, lpost=theOptions%OutputModelState, lrerun=.TRUE.)
    CALL default_stream_setting(theLand%IO_stream, lpost=theOptions%OutputModelState, lrerun=.TRUE., repr=LAND, table=land_table)

    CALL new_stream(theLand%IO_diag_stream, 'land', filetype=theOptions%FileType, lpost=.TRUE., lrerun=.FALSE.)
    CALL default_stream_setting(theLand%IO_diag_stream, lpost=.TRUE., lrerun=.FALSE., repr=LAND, table=land_table)

!---wiso-code
    IF (lwiso) THEN

      CALL new_stream(theLand%IO_wiso_stream, 'js_wiso', filetype=theOptions%FileType, lpost=theOptions%OutputModelState, &
                      lrerun=.TRUE.)
      CALL default_stream_setting(theLand%IO_wiso_stream, lpost=theOptions%OutputModelState, lrerun=.TRUE., repr=LAND, &
                      table=land_table)

      CALL new_stream(theLand%IO_diag_wiso_stream, 'la_wiso', filetype=theOptions%FileType, lpost=.TRUE., lrerun=.FALSE.)
      CALL default_stream_setting(theLand%IO_diag_wiso_stream, lpost=.TRUE., lrerun=.FALSE., repr=LAND, table=land_table)

    END IF
!---wiso-code-end

    CALL new_stream(theLand%IO_veg, 'veg', filetype=theOptions%FileType, lpost=.TRUE., lrerun=.TRUE., &
                    interval=io_time_event(1,'days','first',-43200))
    CALL default_stream_setting(theLand%IO_veg, lpost=.TRUE., lrerun=.TRUE., repr=LAND, table=land_table)

    CALL init_domain(theLand%Grid, theLand%Domain, theOptions%FileType, theLand%IO_stream)

    CALL init_lctlib(theOptions%LctlibFile, theLand%LctLibrary)       ! Gets number of land cover types

    !------------------------------------------------------------------------------------------------------
    ! Determine total number of tiles
    theLand%ntiles = theOptions%ntiles

    if (present(ntile)) ntile =  theLand%ntiles

    ! Here, other subroutines could be called which introduce additional sub-tiling structures 
    ! (e.g. snowbands, distributed precipitation). In this case, theLand%ntiles, theLand%Surface%ntiles,
    ! theLand%Soil%ntiles and theLand%Vegetation%ntiles would be updated here.
    ! ...
    !
    ! The number of soil tiles is set to the total number of land cover types (times additional sub-tiling)
    theLand%Surface%ntiles = theLand%ntiles
    theLand%Soil%ntiles = theLand%ntiles
    theLand%Vegetation%ntiles = theLand%ntiles
    IF (theOptions%UseBethy) THEN
       theLand%Cbalance%ntiles = theLand%ntiles
       theLand%Bethy%ntiles = theLand%ntiles
    ENDIF

    !------------------------------------------------------------------------------------------------------
    ! Now that the total number of tiles is determined, the land surface, vegetation and soil modules can be initialized.
    CALL init_land_surface(theLand%Grid, theLand%Domain, theLand%LctLibrary, theLand%Surface, theOptions%SurfFile,  &
                         theOptions%UseDynveg, lresume, theOptions%filetype, theOptions%ReadCoverFract, &
                         theLand%IO_diag_stream, theLand%IO_stream)
    theLand%Domain%elev => theLand%Surface%elevation   !! Note: when this is a restart run, theLand%Surface%elevation has been allocated but not initialized

    ! default_stream_settings needs to be called here again, as TILES (index into IO_dim_ids) has only been defined in
    ! the previous call to init_land_surface (through add_dim(..., levtyp=70, indx=TILES) in  land_surface_init_io).
    CALL default_stream_setting(theLand%IO_stream, leveltype=TILES)
    IF (test_stream) CALL init_test(theLand%Grid, theLand%Domain, theLand%ntiles)

!---wiso-code
    IF (lwiso) THEN
    
       CALL init_vegetation(theLand%Grid, theLand%Domain, theLand%Surface, theLand%LctLibrary, theLand%Vegetation, &
            theOptions%VegFile, lresume, theOptions%filetype, theLand%IO_diag_stream, theLand%IO_stream, &
            theLand%IO_diag_wiso_stream, theLand%IO_wiso_stream)
           
       CALL init_soil(theLand%Grid, theLand%Domain, theLand%Surface, theLand%SoilParam, theLand%Soil, theOptions%SoilFile, &
            lresume, theOptions%UseAlbedo, theOptions%UseDynveg, theOptions%filetype, theLand%IO_diag_stream, theLand%IO_stream, &
            theLand%IO_diag_wiso_stream, theLand%IO_wiso_stream)
    ELSE

      CALL init_vegetation(theLand%Grid, theLand%Domain, theLand%Surface, theLand%LctLibrary, theLand%Vegetation, &
           theOptions%VegFile, lresume, theOptions%filetype, theLand%IO_diag_stream, theLand%IO_stream)

      CALL init_soil(theLand%Grid, theLand%Domain, theLand%Surface, theLand%SoilParam, theLand%Soil, theOptions%SoilFile, &
           lresume, theOptions%UseAlbedo, theOptions%UseDynveg, theOptions%filetype, theLand%IO_diag_stream, theLand%IO_stream)

    END IF
!---wiso-code-end

    IF (theOptions%UsePhenology) &
         CALL init_phenology(theLand%Grid%nland, theLand%Domain%nland, theLand%Vegetation%ntiles,theLand%Soil%nsoil, &
                             theLand%Vegetation%lai_max, lresume, theOptions%filetype, theLand%IO_veg)

    IF (theOptions%LCChange) &
       CALL init_landcover_change(theLand%Grid%nland, theLand%Domain%nland, theLand%ntiles, &
                                  lresume, theOptions%filetype, theLand%IO_diag_stream, theLand%IO_stream)

    IF (theOptions%UseBethy) THEN
       CALL init_cbalance_bethy(theLand%Grid, theLand%Domain,theLand%Cbalance, theOptions%filetype, theOptions%UseDynveg, &
         theLand%IO_diag_stream, theLand%IO_veg)
       CALL init_bethy(theLand%Grid, theLand%Domain, theLand%Bethy, theOptions%filetype, &
         theLand%IO_diag_stream, theLand%IO_stream)
    ENDIF

    IF (theOptions%UseDynveg) THEN
       CALL init_dynveg(theLand%Grid, theLand%Domain, theLand%Vegetation%ntiles, theLand%Soil%nsoil, &
            lresume, theOptions%filetype, theLand%IO_veg)
    ENDIF

    !IF (options%UseVic) THEN 
    !   CALL vic_init
    !ELSE IF (options%UseEchamLand) THEN
    !   CALL echamland_init
    !ENDIF


    !
    !------------------------------------------------------------------------------------------------------
    ! Read restart files if this is a restarted run
    IF (lresume) THEN
       CALL jsbach_restart
    ENDIF

    IF (PRESENT(grid)) grid = theLand%Grid
    IF (PRESENT(domain)) domain = theLand%Domain
    IF (PRESENT(options)) options = theOptions

    IF (debug) CALL message('jsbach_init','End')

  END SUBROUTINE jsbach_init

  SUBROUTINE jsbach_restart

    USE mo_timer,            ONLY: timer_start, timer_stop, timer_netcdf
    USE mo_io,               ONLY: IO_read_streams

    ! Note: in coupled mode the restart streams are read from the GCM 
    !       (SUBROUTINE *iorestart* in ECHAM5, called from *control*)
    IF (.not. theOptions%Standalone) THEN !! case ECHAM + JSBACH

       ! Convert REAL netCDF representation of surface cover types to INTEGER
       theLand%Surface%cover_type = NINT(theLand%Surface%cover_type_real)

    else !! case JBACH ALONE

        !! -- Read restart files and set time
       IF (theOptions%Timer) CALL timer_start(timer_netcdf)
       CALL io_read_streams
       IF (theOptions%Timer) CALL timer_stop(timer_netcdf)
      
       ! Convert REAL netCDF representation of surface cover types to INTEGER
       theLand%Surface%cover_type = NINT(theLand%Surface%cover_type_real)

    END IF

  END SUBROUTINE jsbach_restart

  SUBROUTINE get_dates(IO_timestep, istep)

    ! For standalone mode only.
    ! Corresponds to ECHAM5 subroutine IO_init in mo_io.f90
    ! This routine doesn't need to be called from ECHAM in coupled mode.

    USE mo_io, ONLY: IO_READ, IO_open, IO_close
    USE mo_netCDF, ONLY: FILE_INFO, NF_GLOBAL, io_get_att_int
    USE mo_time_control, ONLY: INIT_STEP, delta_time, dt_start, &
                               lresume, resume_date, start_date, inp_convert_date, write_date
    USE mo_time_conversion, ONLY : time_native, TC_set, TC_convert

    INTEGER, INTENT(out) :: IO_timestep, istep

    TYPE(FILE_INFO) :: IO_file
    INTEGER :: IO_file_id
    INTEGER :: forecast_date, verification_date   ! YYYYMMDD 
    INTEGER :: forecast_time, verification_time   ! HHMMSS
    TYPE(time_native) :: date_nat

    IF (debug) CALL message('jsbach_init_io','BEGIN')

    ! Initialize time
    IF (p_parallel_io) THEN
       IF (lresume) THEN
          ! Get time information from restart file
          IO_file%opened = .FALSE.
          CALL IO_open(trim(theOptions%RestartPrefix)//'_jsbach', IO_file, IO_READ)
          IO_file_id = IO_file%file_id
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'nstep', istep)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'timestep', IO_timestep)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'fdate', forecast_date)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'ftime', forecast_time)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'vdate', verification_date)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'vtime', verification_time)
          CALL IO_close(IO_file)
       ELSE
          istep = INIT_STEP
          IO_timestep = delta_time
       ENDIF
    ENDIF
       
    IF (p_parallel) THEN
       CALL p_bcast(istep, p_io)
       CALL p_bcast(IO_timestep, p_io)
       IF (lresume) THEN
          CALL p_bcast(forecast_date, p_io)
          CALL p_bcast(forecast_time, p_io)
          CALL p_bcast(verification_date, p_io)
          CALL p_bcast(verification_time, p_io)
       ENDIF
    ENDIF
       
    IF (lresume) THEN
       CALL inp_convert_date(forecast_date, forecast_time, start_date)
       CALL inp_convert_date(verification_date, verification_time, resume_date)
    ELSE
       IF (SUM(dt_start(:)) /= 0) THEN
          CALL TC_set(dt_start(1),dt_start(2),dt_start(3), &
               dt_start(4),dt_start(5),dt_start(6), date_nat)
          CALL TC_convert(date_nat, start_date)
          resume_date = start_date
       ELSE
          CALL finish('jsbach_init_io', 'Start date not set')
       ENDIF
    ENDIF

    IF (p_parallel_io) THEN
       CALL write_date(start_date, 'Start date (initial/restart): ')
       CALL write_date(start_date, 'Resume date (initial/restart): ')
    ENDIF

    IF (debug) CALL message('get_dates','END')

  END SUBROUTINE get_dates

  SUBROUTINE update_current_call(kdim, kland, kblock, mask)

    USE mo_jsbach_grid,     ONLY: update_comm_to_echam5mods
    USE mo_time_control,    ONLY: current_date, lstart, lresume
    USE mo_time_conversion, ONLY: TC_set, TC_convert, time_native, OPERATOR(==)

    INTEGER,           INTENT(in) :: kdim
    INTEGER,           INTENT(in) :: kland
    INTEGER, OPTIONAL, INTENT(in) :: kblock
    LOGICAL, OPTIONAL, INTENT(in) :: mask(kdim)

!!$    INTEGER, ALLOCATABLE, TARGET, SAVE :: kindex_help(:)
    LOGICAL :: mask_help(kdim)    !! Land mask or all true's
    INTEGER,  ALLOCATABLE  :: zint1d(:)
    INTEGER,  ALLOCATABLE  :: zint2d(:,:)

    INTEGER :: i

    TYPE(time_native) :: date_nat

#if defined (__SX__) && defined (_OPENMP)
    INTEGER :: tid
    tid = omp_get_thread_num()
#endif

    IF (PRESENT(kblock)) THEN
       theLand%Domain%kblock = kblock
    ELSE
       theLand%Domain%kblock = 1
       IF (kdim /= theLand%Domain%nlon * theLand%Domain%nlat .AND. kdim /= theLand%Domain%nland) THEN
          call message('jsbach_inter_1d: kdim, nlon, nlat, nland for domain: ',&
               int2string(kdim)//' '//int2string(theLand%Domain%nlon)//' '//int2string(theLand%Domain%nlat)//' '// &
                                                                    int2string(theLand%Domain%nland))
          CALL finish('jsbach_inter_1d','Dimension in call inconsistent with domain')
       ENDIF
    ENDIF

    ! Is this the first call for a new time step on this PE or thread?
    IF (lstart .OR. lresume) THEN
       ! Initialize old_time_step. This is done here (and not e.g. in jsbach_init) so that old_time_step is initialized
       ! on all threads in OpenMP (old_time_step is threadprivate)
       ! Note: start date of experiment must be larger than Jan 1, 00 00:00
       CALL TC_set(0,1,1,0,0,0, date_nat)
       CALL TC_convert(date_nat, old_time_step)
    END IF

    IF (.NOT.(current_date == old_time_step)) THEN
       jsbach_NewTimeStep = .TRUE.
       old_time_step = current_date
    ELSE 
       jsbach_NewTimeStep = .FALSE.
    END IF

    IF (PRESENT(mask)) THEN
       mask_help = mask
    ELSE
       mask_help = .TRUE.
    ENDIF
       
    ! Index of current land boxes into packed vector of land boxes for each PE's domain
    nidx = kland
    ALLOCATE(kindex(kland))
!!$    IF (debug) &
!!$         CALL message('jsbach_inter_1d','Allocated kindex:'// &
!!$                                              int2string(SIZE(kindex))//'  (PE '//int2string(p_pe)//')')

    IF (ALLOCATED(zint1d)) CALL finish('jsbach_inter_1d', 'zint1d already allocated') 

    ALLOCATE(zint1d(theLand%Domain%nland)) ; 
    zint1d = (/ (i, i=1,theLand%Domain%nland) /)        ! Index of domain's land points
    IF (PRESENT(kblock)) THEN
       ALLOCATE(zint2d(theLand%Domain%ndim,theLand%Domain%nblocks))
       zint2d(:,:) = UNPACK(zint1d, MASK=theLand%Domain%mask, FIELD=0) 
                                              ! Location of all domain land points in (ndim,nblocks) grid
       kindex(1:kland) = PACK(zint2d(1:kdim,kblock), MASK=mask_help)
                                              ! Note: kdim can be smaller than grid%ndim for last block!
                                              ! Note: if no mask given, mask=.T. and kland=kdim
       DEALLOCATE(zint2d)
    ELSE   ! kblock not present, i.e. kdim = grid%ndim = grid%nlon * grid%nlat, i.e. whole domain in one call (block, nblocks=1)
       kindex(:) = zint1d(:)     ! Note: kland = grid%nland
    ENDIF
    kstart = kindex(1)
    kend = kindex(kland)
    IF (ANY(kindex == 0)) THEN
       CALL message('mo_jsbach_interface - update_current_call','nidx, kstart, kend = '// &
                                      int2string(nidx)//' '//int2string(kstart)//' '//int2string(kend))
       CALL finish('              -  ','Grid of calling programm inconsistent with JSBACH grid')
    ENDIF

    DEALLOCATE(zint1d)

    ! Is this the last call (block) for this time step?
!!$    IF (theLand%Domain%kblock == theLand%Domain%nblocks) THEN
!!$       theLand%Domain%LastBlock = .TRUE.
!!$       IF (kdim /= theLand%Domain%ndimz) &
!!$            CALL finish('jsbach_inter_1d','Wrong vector length in last block of domain')
!!$    ELSE IF (theLand%Domain%kend == theLand%Domain%nland) THEN
!!$       theLand%Domain%LastBlock = .TRUE.
!!$       IF (kdim /= theLand%Domain%ndim) &
!!$            CALL finish('jsbach_inter_1d','Wrong vector length in last block of domain')
!!$    ELSE
!!$       theLand%Domain%LastBlock = .FALSE.
!!$    ENDIF

!!$    IF (debug) THEN
!!$       CALL message('jsbach_interface','PE '//int2string(p_pe)//': ')
!!$                    ' First call for this time step - ',jsbach_NewTimeStep
!!$       print*,'               .', ' Last call for this time step  - ',theLand%Domain%LastBlock
!!$       IF (theLand%Domain%nblocks > 1) THEN
!!$          CALL message('               .', ' Processing block # '//int2string(kblock)//'('//&
!!$               int2string(theLand%Domain%nland)//', '//int2string(kland)//', '//int2string(kstart)//', '//int2string(kend))
!!$       ELSE
!!$          CALL message('               .', ' Processing whole domain in one call')
!!$       ENDIF
!!$    ENDIF

    IF (test_stream) CALL init_test_2(kstart,kend, kblock , &
         theLand%Surface%cover_fract(kstart:kend,:), &
         theLand%Surface%is_present(kstart:kend,:) .AND. &
         .NOT. theLand%Surface%is_lake(kstart:kend,:), theLand%Domain)

#if defined (__SX__) && defined (_OPENMP)
    IF (tid == 0) &
#endif
         CALL update_comm_to_echam5mods(theLand%Grid, theLand%Domain)

  END SUBROUTINE update_current_call

END MODULE mo_jsbach_interface


!Local Variables:
!mode: f90
!fill-column: 100
!End:
