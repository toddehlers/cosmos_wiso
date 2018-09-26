!!
!! This module belongs to the standalone version of JSBACH
!!
!+ Definition for offline driving secundary parameters (other than forcing [= weather])
!
Module mo_jsbalone

  !
  ! Description:
  !   Derives additional parameters for offline driving of jsbach
  !   like wind profile, stability params, Richtmeyr-Morton Coeffs
  !   ect. 
  !   Some additional parameters use values from the last time step
  !   
  !
  ! Current Code Owner:  jsbach_admin
  ! 
  ! History:
  !
  ! Version       Date                Comment
  ! -------       ------              -------
  ! 0.1           2005/03/22          Org. Code admin
  !
  !
  ! Modules used:
  !
  USE mo_kind,              ONLY: dp
  USE mo_jsbach,            ONLY: options_type
  USE mo_jsbach_grid,       ONLY: grid_type, domain_type
  USE mo_exception,         ONLY: finish, message, int2string,real2string
  USE mo_mpi,               ONLY: p_parallel_io, p_parallel, p_io, p_bcast
  USE mo_jsbalone_forcing,  ONLY: forcing_type
  USE mo_linked_list,       ONLY: t_stream
 
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Drive_type, init_driving, update_driving
  !
  !-----------------------------------------------------------------------------------------------------------------
  !
  TYPE Drive_type
     REAL(dp), POINTER, Dimension(:)  :: & !! Dimension (nland)
          Air_Temp_old , &                !! Air Temperature from last time step
          geopot, &
          wind, &
          wind10, &
          temp_air, &
          qair, &
          eair, &
          precip_rain, &
          precip_snow, &
          co2, &
          lwdown, &
          vis_net, &
          vis_frac, &
          nir_net, &
          nir_frac, &
          czenith, &
          pressure, &
          etacoef, &
          eqacoef, &
          etbcoef, &
          eqbcoef, &
          for_rau, &
          sensible_heat, &
          evap_act, &
          soil_temp, &
          evap_pot, &
          cdrag, &
          zchl, &
          latent, &
          test1, &
          test2, &
          test3, &
          test4, &
          test5, &
          test6, &
          test7, &
          test8, &
          test9, &
! DIAG TYPE : upper elements have to be cleared !!
          Air_Temp_old_acc , &                !! Air Temperature from last time step
          geopot_acc, &
          wind_acc, &
          wind10_acc, &
          temp_air_acc, &
          qair_acc, &
          eair_acc, &
          precip_rain_acc, &
          precip_snow_acc, &
          co2_acc, &
          lwdown_acc, &
          vis_net_acc, &
          vis_frac_acc, &
          nir_net_acc, &
          nir_frac_acc, &
          czenith_acc, &
          pressure_acc, &
          etacoef_acc, &
          eqacoef_acc, &
          etbcoef_acc, &
          eqbcoef_acc, &
          for_rau_acc, &
          sensible_heat_acc, &
          evap_act_acc, &
          soil_temp_acc, &
          evap_pot_acc, &
          cdrag_acc, &
          zchl_acc, &
          latent_acc, &
          test1_acc, &
          test2_acc, &
          test3_acc, &
          test4_acc, &
          test5_acc, &
          test6_acc, &
          test7_acc, &
          test8_acc, &
          test9_acc
  END TYPE Drive_type

  !
  !-----------------------------------------------------------------------------------------------------------------
  !
  LOGICAL, SAVE  :: module_configured = .FALSE. 
  LOGICAL, SAVE  :: module_initialized = .FALSE.

  TYPE(t_stream),  POINTER, SAVE :: IO_driving      !! Memory stream for driving variables in offline mode

  !
  !-----------------------------------------------------------------------------------------------------------------
  !
  CONTAINS
  !
  !-----------------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE config_driving

    module_configured = .true.
    
  END SUBROUTINE config_driving
  !
  !=================================================================================================================
  !
  SUBROUTINE driving_init_memory(g_nland, l_nland, ntiles, Drv_Var, stream)

    USE mo_linked_list, ONLY : LAND
    USE mo_memory_base, ONLY : add =>add_stream_element
    USE mo_netCDF,      ONLY : max_dim_name
    
    INTEGER         , INTENT(in)             :: g_nland, l_nland, ntiles
    TYPE(drive_type), INTENT(inout)          :: Drv_Var
    TYPE(t_stream)  , POINTER, OPTIONAL      :: stream

    INTEGER                     :: dim1p(1), dim1(1)
    CHARACTER(LEN=max_dim_name) :: dim1n(1)


    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    CALL add(stream, 'air_temp_old', Drv_Var%Air_Temp_old, dim1p, dim1, &
         laccu=.false., lpost=.false., dimnames=dim1n, code=1)
    CALL add(stream, 'geopot'      , Drv_Var%geopot, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'wind'        , Drv_Var%wind, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'wind10'      , Drv_Var%wind10, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'temp_air'    , Drv_Var%temp_air, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'qair'        , Drv_Var%qair, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'eair'        , Drv_Var%eair, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'precip_rain' , Drv_Var%precip_rain, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'precip_snow' , Drv_Var%precip_snow, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'co2'         , Drv_Var%co2, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'lwdown'      , Drv_Var%lwdown, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'vis_net'     , Drv_Var%vis_net, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'vis_frac'    , Drv_Var%vis_frac, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'nir_net'     , Drv_Var%nir_net, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'nir_frac'    , Drv_Var%nir_frac, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'czenith'     , Drv_Var%czenith, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'pressure'    , Drv_Var%pressure, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'etacoef'     , Drv_Var%etacoef, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'etbcoef'     , Drv_Var%etbcoef, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'eqacoef'     , Drv_Var%eqacoef, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'eqbcoef'     , Drv_Var%eqbcoef, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'for_rau'     , Drv_Var%for_rau , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'sensible_heat', Drv_Var%sensible_heat , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'evap_act'    , Drv_Var%evap_act , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'soil_temp'   , Drv_Var%soil_temp , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'evap_pot'    , Drv_Var%evap_pot , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'cdrag'       , Drv_Var%cdrag , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'zchl'        , Drv_Var%zchl  , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'latent'      , Drv_Var%latent  , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'test1'        , Drv_Var%test1 , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'test2'        , Drv_Var%test2 , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'test3'        , Drv_Var%test3 , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'test4'        , Drv_Var%test4 , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'test5'        , Drv_Var%test5 , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'test6'        , Drv_Var%test6 , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'test7'        , Drv_Var%test7 , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'test8'        , Drv_Var%test8 , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'test9'        , Drv_Var%test9 , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
!
!   DIAG TYPE: Upper elements have to be cleaned ! 
!   
    CALL add(stream, 'air_temp_old_acc', Drv_Var%Air_Temp_old_acc, dim1p, dim1, &
         laccu=.true., lpost=.true., dimnames=dim1n, code=1)
    CALL add(stream, 'geopot_acc'      , Drv_Var%geopot_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'wind_acc'        , Drv_Var%wind_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'wind10_acc'      , Drv_Var%wind10_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'temp_air_acc'    , Drv_Var%temp_air_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'qair_acc'        , Drv_Var%qair_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'eair_acc'        , Drv_Var%eair_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'precip_rain_acc' , Drv_Var%precip_rain_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'precip_snow_acc' , Drv_Var%precip_snow_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'co2_acc'         , Drv_Var%co2_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'lwdown_acc'      , Drv_Var%lwdown_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'vis_net_acc'     , Drv_Var%vis_net_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'vis_frac_acc'    , Drv_Var%vis_frac_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'nir_net_acc'     , Drv_Var%nir_net_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'nir_frac_acc'    , Drv_Var%nir_frac_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'czenith_acc'     , Drv_Var%czenith_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'pressure_acc'    , Drv_Var%pressure_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'etacoef_acc'     , Drv_Var%etacoef_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'etbcoef_acc'     , Drv_Var%etbcoef_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'eqacoef_acc'     , Drv_Var%eqacoef_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'eqbcoef_acc'     , Drv_Var%eqbcoef_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'for_rau_acc'     , Drv_Var%for_rau_acc , dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'sensible_heat_acc', Drv_Var%sensible_heat_acc , dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'evap_act_acc'    , Drv_Var%evap_act_acc , dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'soil_temp_acc'   , Drv_Var%soil_temp_acc , dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'evap_pot_acc'    , Drv_Var%evap_pot_acc , dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'cdrag_acc'       , Drv_Var%cdrag_acc , dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'zchl_acc'        , Drv_Var%zchl_acc  , dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'latent_acc'      , Drv_Var%latent_acc  , dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'test1_acc'        , Drv_Var%test1_acc , dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'test2_acc'        , Drv_Var%test2_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'test3_acc'        , Drv_Var%test3_acc , dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'test4_acc'        , Drv_Var%test4_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'test5_acc'        , Drv_Var%test5_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'test6_acc'        , Drv_Var%test6_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'test7_acc'        , Drv_Var%test7_acc , dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'test8_acc'        , Drv_Var%test8_acc , dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
    CALL add(stream, 'test9_acc'        , Drv_Var%test9_acc , dim1p, dim1, dimnames=dim1n, code=1, laccu=.true.)
!    CALL add(stream, ''        , Drv_Var% , dim1p, dim1, dimnames=dim1n, code=1)


  END SUBROUTINE driving_init_memory
  !
  !=================================================================================================================
  !
  SUBROUTINE init_driving(grid, domain, global_options, forcingFields, drv_var, ntiles)         

    USE mo_linked_list, ONLY : LAND
    USE mo_memory_base, ONLY : new_stream, default_stream_setting

    TYPE(grid_type)   ,    INTENT(in)   :: grid
    TYPE(domain_type) ,    INTENT(in)   :: domain
    TYPE(options_type),    INTENT(in)   :: global_options
    TYPE(forcing_type),    INTENT(in)   :: forcingFields   !! Forcing fieldes for first initialisation (new start) 
    INTEGER           ,    INTENT(in)   :: ntiles
    TYPE(drive_type)  ,    INTENT(inout):: Drv_Var

    CALL config_driving

    CALL new_stream(IO_driving, 'driving', filetype=global_options%filetype, lpost=.true., lrerun=.TRUE.)

    CALL default_stream_setting(IO_driving, lpost=.true., lrerun=.TRUE., repr=LAND)

    call driving_init_memory(grid%nland, domain%nland, ntiles, Drv_Var, stream=IO_driving)

    module_initialized = .TRUE.
    
  END SUBROUTINE init_driving
  !
  !=================================================================================================================
  !
  SUBROUTINE update_driving( &
       kdim, &                            
       kland, &                          
       drv_var, &                     
       geopot, wind, wind10, temp_air, qair_forcing,  &
       precip_rain, precip_snow, &
       lwdown, &
       sw_vis_net, &                       
       sw_vis_frac_diffuse, &              
       sw_nir_net, &                     
       sw_nir_frac_diffuse, &            
       pressure, &
       czenith, &
       CO2_concentration, &
       LayerHeight, &
       FirstCall, &
       fluxsens, &
       vevapp , &
       soil_temp_new, &
       pgeom &
       )

    USE mo_time_control,     ONLY: delta_time
    USE mo_jsbach_constants, Only: SpecificHeatDryAirConstPressure, GasConstantDryAir, Gravity
    USE mo_convect_tables,   ONLY: tlucua, jptlucu1, jptlucu2, lookuperror, lookupoverflow
 

    INTEGER,           INTENT(in)    :: kdim                      !! Length of vectors (if call from ECHAM5 this
    INTEGER,           INTENT(in)    :: kland                     !! Number of land points in vectors
    TYPE(drive_type), INTENT(inout)  :: drv_var                   !! driving variables 
    REAL(dp),              INTENT(in)    :: geopot(kdim)              !! Geopotential at lowest atmospheric level [m^2/s^2]
    REAL(dp),              INTENT(in)    :: wind(kdim)                !! Lowest level wind speed [m/s]
    REAL(dp),              INTENT(in)    :: wind10(kdim)              !! 10m wind speed [m/s] (for update_surface_down)
    REAL(dp),              INTENT(in)    :: temp_air(kdim)            !! Lowest level air temperature [Kelvin]
    REAL(dp),              INTENT(in)    :: qair_forcing(kdim)        !! Lowest level specific humidity
    REAL(dp),              INTENT(in)    :: precip_rain(kdim)         !! Precipitation as rain [kg/(m^2 s)]
    REAL(dp),              INTENT(in)    :: precip_snow(kdim)         !! Precipitation as snow [kg/(m^2 s)]
    REAL(dp),              INTENT(in)    :: CO2_concentration(kdim)   !! Atmospheric CO2 concentration [kg(CO2)/kg(air)]
    REAL(dp),              INTENT(in)    :: lwdown(kdim)              !! Downward longwave flux
    REAL(dp),              INTENT(in)    :: sw_vis_net(kdim)          !! solar radiation in the visible band reaching surface [W/m^2]
    REAL(dp),              INTENT(in)    :: sw_vis_frac_diffuse(kdim) !! fraction of diffuse radiation contained in sw_vis_net
    REAL(dp),              INTENT(in)    :: sw_nir_net(kdim)          !! solar radiation in the near infrared band reaching surface [W/m^2]
    REAL(dp),              INTENT(in)    :: sw_nir_frac_diffuse(kdim) !! fraction of diffuse radiation contained in sw_nir_net
    REAL(dp),              INTENT(in)    :: czenith(kdim)             !! Cosine of solar zenith angle
    REAL(dp),              INTENT(in)    :: pressure(kdim)            !! Surface pressure
    REAL(dp),              INTENT(in)    :: LAyerHeight               !! Defines lowest layer height-> where measurements are taken
    LOGICAL,           INTENT(in)    :: FirstCall                 !! 
    REAL(dp),              INTENT(inout)    :: vevapp(kdim)              !! total evaporation from LAST TIME STEP (calc. by jsabch)
    REAL(dp),              INTENT(inout)    :: fluxsens(kdim)            !! sensible heat flux from last time step (calc. by jsbach)
    REAL(dp),              INTENT(in)    :: soil_temp_new(kdim)       !! new soil temperature (calc. by jsbach)
    REAL(dp),              INTENT(in)    :: pgeom(kdim)               !! Geopotential
    REAL(dp), PARAMETER :: relaxcoef = 1._dp
    REAL(dp), PARAMETER :: ckap = 0.4_dp
    REAL(dp), DIMENSION(kland) :: for_qair, for_eair, old_qair, old_eair, qair_obs, eair_obs, for_rau
    REAL(dp), DIMENSION(kland) :: zlflu, zlev_vec, for_tair
    REAL(dp), DIMENSION(kland) :: zchnl, zcons, zgeom
    REAL(dp) ::  zkap, zcons12,i

    !-------------------------------------------------------------------------------------------------
    ! Write fields to output stream
    !-------------------------------------------------------------------------------------------------
    
    Drv_Var%geopot(1:kland) = geopot(:)
    Drv_Var%wind = wind
    Drv_Var%wind10 = wind10
    Drv_Var%temp_air = temp_air
    Drv_Var%qair = qair_forcing 
    Drv_Var%precip_rain = precip_rain 
    Drv_Var%precip_snow = precip_snow 
    Drv_Var%co2 = CO2_concentration
    Drv_Var%lwdown = lwdown 
    Drv_Var%vis_net = sw_vis_net 
    Drv_Var%vis_frac = sw_vis_frac_diffuse 
    Drv_Var%nir_net = sw_nir_net 
    Drv_Var%nir_frac = sw_nir_frac_diffuse 
    Drv_Var%czenith = czenith
    Drv_Var%pressure = pressure
    Drv_Var%cdrag(1:kland) = 1600._dp
    Drv_Var%zchl(1:kland) = 0.06_dp

    Drv_Var%geopot_acc(1:kland) = Drv_Var%geopot_acc(1:kland) + delta_time*geopot(:)
    Drv_Var%wind_acc = Drv_Var%wind_acc + delta_time*wind
    Drv_Var%wind10_acc = Drv_Var%wind10_acc + delta_time*wind10
    Drv_Var%temp_air_acc = Drv_Var%temp_air_acc + delta_time*temp_air
    Drv_Var%qair_acc = Drv_Var%qair_acc + delta_time*qair_forcing 
    Drv_Var%precip_rain_acc = Drv_Var%precip_rain_acc + delta_time*precip_rain 
    Drv_Var%precip_snow_acc = Drv_Var%precip_snow_acc + delta_time*precip_snow 
    Drv_Var%co2_acc = Drv_Var%co2_acc + delta_time*CO2_concentration
    Drv_Var%lwdown_acc = Drv_Var%lwdown_acc + delta_time*lwdown 
    Drv_Var%vis_net_acc= Drv_Var%vis_net_acc + delta_time*sw_vis_net 
    Drv_Var%vis_frac_acc = Drv_Var%vis_frac_acc + delta_time*sw_vis_frac_diffuse 
    Drv_Var%nir_net_acc = Drv_Var%nir_net_acc + delta_time*sw_nir_net 
    Drv_Var%nir_frac_acc = Drv_Var%nir_frac_acc + delta_time*sw_nir_frac_diffuse 
    Drv_Var%czenith_acc = Drv_Var%czenith_acc + delta_time*czenith
    Drv_Var%pressure_acc = Drv_Var%pressure_acc + delta_time*pressure
    Drv_Var%cdrag_acc(1:kland) = Drv_Var%cdrag_acc + delta_time*1600.
    Drv_Var%zchl_acc(1:kland) = Drv_Var%zchl_acc + delta_time*0.06
    Drv_Var%evap_act_acc(1:kland) = Drv_Var%evap_act_acc(1:kland) + delta_time *Drv_Var%evap_act
    Drv_Var%latent_acc(1:kland) = Drv_Var%latent_acc(1:kland) + delta_time * Drv_Var%latent
    Drv_Var%sensible_heat_acc(1:kland) = Drv_Var%sensible_heat_acc(1:kland) + delta_time * Drv_Var%sensible_heat
!    Drv_Var%_acc(1:kland) = Drv_Var%_acc(1:kland) + delta_time
!    Drv_Var%_acc(1:kland) = Drv_Var%_acc(1:kland) + delta_time
!    Drv_Var%_acc(1:kland) = Drv_Var%_acc(1:kland) + delta_time
!    Drv_Var%_acc(1:kland) = Drv_Var%_acc(1:kland) + delta_time

    !-------------------------------------------------------------------------------------------------
    ! Some Constants
    !-------------------------------------------------------------------------------------------------
    
    zcons12       = 1.5_dp * delta_time * 9.81_dp / 287.05_dp


    !-------------------------------------------------------------------------------------------------
    ! If First Call Initialise some fields
    !-------------------------------------------------------------------------------------------------
    IF(FirstCall) then
       old_qair(1:kland) = qair_forcing(1:kland)
       old_eair(1:kland) = SpecificHeatDryAirConstPressure * temp_air + Gravity * LayerHeight
       vevapp(1:kland)   = -5.0E-08_dp
       fluxsens(1:kland) = -100.0_dp
       Drv_Var%for_rau(1:kland) = pressure(:)/(GasConstantDryAir*temp_air(:))
    END IF
    !----------------------------------------------------------------------------------------
    ! Hand over
    !----------------------------------------------------------------------------------------
    qair_obs(1:kland)    = qair_forcing(1:kland)
    eair_obs(1:kland)    = SpecificHeatDryAirConstPressure * temp_air + Gravity * LayerHeight
!!$    qair_obs(1:kland)    = 0.002_dp
!!$    eair_obs(1:kland)    = SpecificHeatDryAirConstPressure * 263.0_dp + Gravity * 2.0_dp
    for_rau(1:kland)     = Drv_Var%for_rau(1:kland)
    old_qair(:)          = Drv_Var%qair(1:kland)
    old_eair(:)          = Drv_Var%eair(1:kland)
    zgeom                = pgeom
    !IF IMPLICITE Coupling:
    zlflu = delta_time / relaxcoef
    ! ELSE :
    ! zlflu = LayerHeight / 2._dp

    !------------------------------------------------------------------------------------
    ! Approximation of cdrag
    !------------------------------------------------------------------------------------
!    WHERE(zgeom(:) .lt. 1._dp) zgeom = 10000.0_dp
!       zchnl(:)          = (ckap / LOG(1._dp+ pgeom(:) / (9.81_dp * z0h(:))))**2
where(zgeom .lt. 1._dp) zgeom = 1000._dp

       zchnl(:)          = (ckap / LOG(1._dp+ zgeom(:) / (9.81_dp)))**2
       zcons(:)          = zcons12 * pressure(:) / temp_air(:)
       Drv_Var%cdrag(1:kland)  = zcons(:) * wind(:) * zchnl(:)
!do i=1,kland
!write(*,*)'fehler',i, Drv_Var%cdrag(i), zcons(i), wind(i),pressure(i), pgeom(i)
!end do
    !------------------------------------------------------------------------------------

    for_qair(:) = old_qair(:) &
         & +delta_time / zlflu * (relaxcoef * (qair_obs(:) - old_qair(:)) &
         &            + vevapp(:) / (delta_time * for_rau(:)))

    for_eair(:) = old_eair(:) &
         & + delta_time / zlflu * (relaxcoef * (eair_obs(:) - old_eair(:)) &
         & + fluxsens(:) / for_rau(:))

    ! TEST
    Drv_Var%test1(1:kland) = old_qair(:)
    Drv_Var%test2(1:kland) = old_eair(:)
    Drv_Var%test3(1:kland) = delta_time
    Drv_Var%test4(1:kland) = zlflu
    Drv_Var%test5(1:kland) = relaxcoef
    Drv_Var%test6(1:kland) = qair_obs(:)
    Drv_Var%test7(1:kland) = (qair_obs(:) - old_qair(:))
    Drv_Var%test8(1:kland) = vevapp(:)
    Drv_Var%test9(1:kland) =  vevapp(:) / (delta_time * for_rau(:))

    !TEST

    !------------------------------------------------------------------------------------------------------------------
    for_tair(:) = (for_eair(:) - Gravity * LayerHeight) / SpecificHeatDryAirConstPressure
    !-------
!    epot_sol(:,:) = cp_air*temp_sol_NEW(:)
    !-------


    for_rau(:) =  pressure(:)/(GasConstantDryAir*temp_air(:))
!!$    !-------


!!$    !---------------------------------------------------------------------------------------------------------------
!!$    ! Computation of Relaxation Coefficients
!!$    !---------------------------------------------------------------------------------------------------------------    
!!$
!!$    
!!$    
!!$
!!$    !---------------------------------------------------------------------------------------------------------------
!!$    ! Computation of Richtmeyr-Morton Coefficients
!!$    !---------------------------------------------------------------------------------------------------------------    

    Drv_Var%etacoef = 0.0_dp
    Drv_Var%eqacoef = 0.0_dp
!    Drv_Var%etbcoef(1:kland) = relaxcoef * (eair_obs(:) - old_eair(:))
!    Drv_Var%eqbcoef(1:kland) = relaxcoef * (qair_obs(:) - old_qair(:))
    Drv_Var%etbcoef(1:kland) = relaxcoef * (for_eair)
    Drv_Var%eqbcoef(1:kland) = relaxcoef * (for_qair)
     
    !---------------------------------------------------------------------------------------------------------------
    ! Save old Values of qair and eair
    !---------------------------------------------------------------------------------------------------------------    
    Drv_Var%qair(1:kland) = for_qair(:)
    Drv_Var%eair(1:kland) = for_eair(:)
    Drv_Var%for_rau(1:kland) = for_rau(:)

  END SUBROUTINE update_driving
  !
  !=================================================================================================================
  !
End Module mo_jsbalone
