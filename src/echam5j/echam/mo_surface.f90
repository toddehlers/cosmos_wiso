MODULE mo_surface
  !---------------------------------------------------------------------
  ! USE statements 
  !---------------------------------------------------------------------
  USE mo_kind,                 ONLY: dp
  USE mo_doctor,               ONLY: nout
  USE mo_surface_types
  USE mo_surface_memory,       ONLY: land
  USE mo_surface_memory,       ONLY: ocean
  USE mo_surface_memory,       ONLY: ice 
  USE mo_surface_memory,       ONLY: box
  USE mo_surface_memory,       ONLY: surface
  USE mo_jsbach_interface,     ONLY: jsbach_inter_1d
  !---------------------------------------------------------------------
  IMPLICIT NONE 
  ! 
  !---------------------------------------------------------------------
CONTAINS
  !---------------------------------------------------------------------
  ! SUBROUTINES: update_surface is one and only routine here 
  !---------------------------------------------------------------------
  SUBROUTINE update_surface (                         & ! logistic parameters
       current_nproma,                                & ! vector length
       jrow,                                          & ! krow
       levels,                                        & ! number of levels (lowest level)
       levels_plus_1,                                 & ! one up
       levels_minus_1,                                & ! one down
  !!-----------------------------!! atmospheric conditions lowest level-
       cloud_water,                                   & ! total amount of water in liquid droplets
       cloud_ice,                                     & ! total amount of water in cloud ice
       geopotential,                                  & ! geopot.
       atm_temp,                                      & ! atmosphere temperature (lowest layer)
       atm_spec_humidity,                             & ! specific humidity atmosphere (lowest layer)
       atm_full_lev_press,                            & ! full level pressure
       atm_half_lev_press,                            & ! half level pressure
       wind_u,                                        & ! wind speed zonal
       wind_v,                                        & ! wind speed meridional
       rain,                                          & ! total rain fall
       snow,                                          & ! total snow fall
       longwave_net,                                  & ! long wave net radiation
       sw_vis,                                        & ! net surface visible
       sw_vis_frac_diffuse,                           & ! fraction of diffuse visible
       sw_nir,                                        & ! net surface NIR
       sw_nir_frac_diffuse,                           & ! fraction of diffuse NIR
       cos_zenith_angle,                              & ! zenith angle for radiation between l_trigrad
       cloud_cover,                                   & ! integrated cloud cover [0,...,1]
       !!-----------------------------------!! -------------------------
       ocean_temp,                                    & ! SST is read by echam routines and handed over here
       ocean_u,                                       & ! ocean_u_velocity
       ocean_v,                                       & ! ocean_v_velocity
       ice_depth,                                     & ! ICE depth comes from clsst 
       !!-----------------------------------!! -------------------------
       zcfh,                                          & ! 
       zebsh,                                         & !
       zqdif_pre,                                     & !
       ztdif_pre,                                     & !
       zudif,                                         & !
       zvdif,                                         & !
       zghabl,                                        & !
       pi0,                                           & !       
       ptrsol,                                        & !
       pco2_concentration,                            & ! CO2 concentration in lowest level (mass mixing ratio)
       !!-------OUTPUT--------------------------------------------------
       palbedo,              palbedo_vis,             &
       palbedo_nir,                                   &
       ptrflw,               ptrfli,                  &
       psofll,               psoflw,                  &
       psofli,               ptrfllac,                &
       ptrflwac,             ptrfliac,                &
       psofllac,             psoflwac,                &
       psofliac,             palsol,                  &
       palsoi,               palsow ,                 &
       ustarm,               momentum_exchange_coeff, &
       tkevn_cond,           ztdif_new,               &
       zqdif_new,            ztvh,                    &
       zqsurf,               zth,                     &
       pwind10w ,            pu10 ,                   &
       pv10 ,                pwimax ,                 &
       pwind10,              pdew2,                   &
       ptemp2,               pt2max,                  &
       pt2min,               pevaplac,                &
       pevapwac,             pevapiac,                &
       pevap,                pahfllac,                &
       pahflwac,             pahfliac,                &
       pahfl,                pahfslac,                &
       pahfswac,             pahfsiac,                &
       pahfs,                pqhfla,                  &
       pevapw,                                        &
       pevapi,               pahfsl,                  &
       pahfsw,               pahfsi,                  &
       pevapot,                                       &
       pahflw,                                        &
       psni,                 pahfice,                 &
       pfluxres,             pqres,                   &
       pahfcon,              pahfres,                 &
       ptsi,                 ptslnew,                 &
       pzti,                 pzteffl4,                &
       pztsnew,              ptsurf,                  &
       paz0w,                paz0i,                   &
       paz0l,                paz0,                    &
       pustrl,               pvstrl,                  &
       pustrw,               pvstrw,                  &
       pustri,               pvstri,                  &
       prunoff,              pdrainage,               &
       pustr,                pvstr,                   &
       ptte_corr,            ptsw_new,                &
       pseaice,              pseaice_new,             &
       psiced_new,           pradtemp_old,            &
       palac, pros_hd, pdrain_hd,                     &
       pco2_flux_ocean, pco2_flux_land, pco2_flux,    &
!---wiso-code
       lwiso, kwiso, lwisofracl,                      &
       ! Input - ISOTOPES
       wiso_atm_spec_humidity,                        &
       wiso_rain,            wiso_snow,               &
       zwisoqdif_pre,        zwisokinw,   zwisokinl,  &
       pwisosw_d,                                     &
       ! Output - ISOTOPES
       zwisoqdif_new,                                 &
       pwisoevaplac,         pwisoevapwac,            &
       pwisoevapiac,         pwisoevap,               &
       pwisoqhfla,                                    &
       pwisoevapw,                                    &
       pwisoevapi,           pwisoevapot,             &
       pwisorunoff,          pwisodrainage,           &
       pwisoalac, pwisoros_hd, pwisodrain_hd          &
!---wiso-code-end
       )

  !---------------------------------------------------------------------
  ! USE statements 
  !---------------------------------------------------------------------
    USE mo_memory_g3b, ONLY:                    &
         slf,                                   &   !! Sea-land fraction (contains fraction of land [0,....,1]
         friac, tslm1, ws, sn, wl,              &
         gld, snmel, apmegl, snacl, rogl, wsmx
!---wiso-code
    USE mo_wiso, ONLY:                          &
         lwiso_rerun, tnat
    USE mo_memory_wiso, ONLY:                   &
         wisows, wisosn, wisowl,                &
         wisogld, wisosnmel, wisoapmegl,        &
         wisosnacl, wisorogl
    USE mo_time_control,         ONLY: lresume
!---wiso-code-end

    USE mo_surface_land
    USE mo_surface_ocean
    USE mo_surface_ice
    USE mo_surface_boundary
    USE mo_time_control,         ONLY: lstart, delta_time
    USE mo_orbit,                ONLY: declination

    !-------------------------------------------------------------------
    !                         :: NEW VARIABLE NAME          :: OLD ECHAM NAME AND DIMENSIONS
    !-------------------------------------------------------------------
    INTEGER, INTENT(in)       :: current_nproma             !! kbdim
    INTEGER, INTENT(IN)       :: jrow                       !! krow
    INTEGER, INTENT(IN)       :: levels                     !! klev
    INTEGER, INTENT(in)       :: levels_plus_1              !! klevp1
    INTEGER, INTENT(in)       :: levels_minus_1             !! klevm1
!---wiso-code
    LOGICAL, INTENT(IN)       :: lwiso
    INTEGER, INTENT(IN)       :: kwiso
    INTEGER, INTENT(IN), OPTIONAL :: lwisofracl
!---wiso-code-end
    REAL(dp),    INTENT(in)       :: geopotential(:)            !! pgeom1(kdim,klev)
    REAL(dp),    INTENT(in)       :: wind_u(:)                  !! pum1
    REAL(dp),    INTENT(in)       :: wind_v(:)                  !! pvm1
    REAL(dp),    INTENT(in)       :: atm_temp(:)                !! ptm1
    REAL(dp),    INTENT(in)       :: atm_spec_humidity(:)       !! pqm1
    REAL(dp),    INTENT(in)       :: rain(:)                    !! rsfl + rsfc from physc
    REAL(dp),    INTENT(in)       :: snow (:)                   !! ssfl + ssfc from physc
    REAL(dp),    INTENT(in)       :: longwave_net(:)            !! pemter
    REAL(dp),    INTENT(in)       :: sw_vis(:)                  !! net surface visible
    REAL(dp),    INTENT(in)       :: sw_vis_frac_diffuse(:)     !! fraction of diffuse visible
    REAL(dp),    INTENT(in)       :: sw_nir(:)                  !! net surface NIR
    REAL(dp),    INTENT(in)       :: sw_nir_frac_diffuse(:)     !! fraction of diffuse visible
    REAL(dp),    INTENT(in)       :: atm_full_lev_press(:)      !! papm1
    REAL(dp),    INTENT(in)       :: cos_zenith_angle(:)        !! zm0
    REAL(dp),    INTENT(in)       :: atm_half_lev_press(:,:)    !! paphm1(kbdim, klevp1)
    REAL(dp),    INTENT(in)       :: cloud_water(:)             !! pxlm1
    REAL(dp),    INTENT(in)       :: cloud_ice(:)               !! pxim1
    REAL(dp),    INTENT(in)       :: cloud_cover(:)             !! paclc
    REAL(dp),    INTENT(in)       :: ocean_temp(:)              !! ptsw
    REAL(dp),    INTENT(in)       :: ocean_u(:)                 !! ocu
    REAL(dp),    INTENT(in)       :: ocean_v(:)                 !! ocv
    REAL(dp),    INTENT(in)       :: zcfh(:,:)                  !! zcfh: dimless exch coeff. from diffusion scheme
    REAL(dp),    INTENT(in)       :: zebsh(:,:)                 !! zebsh from heat diffusion scheme
    REAL(dp),    INTENT(in)       :: zqdif_pre(:,:)             !! zqdif from heat diffusion scheme
    REAL(dp),    INTENT(in)       :: ztdif_pre(:,:)             !! ztdif from heat diffusion scheme
    REAL(dp),    INTENT(in)       :: zudif(:)                   !! zudif from momentum transfer scheme
    REAL(dp),    INTENT(in)       :: zvdif(:)                   !! zvdif from momentum transfer scheme
    REAL(dp),    INTENT(in)       :: zghabl(:)                  !! zghabl from pbl-height
    REAL(dp),    INTENT(in)       :: pi0(:)                     !! solar incidence factor for solar radiation
    REAL(dp),    INTENT(in)       :: ptrsol(:)                  !! ptrsol solar radiation
    REAL(dp),    INTENT(in)       :: ice_depth(:)               !! psiced comes from clsst
    REAL(dp),    INTENT(in)       :: pseaice(:)                 !! seaice comes from clsst
    REAL(dp),    INTENT(in)       :: pco2_concentration(:)      !! CO2 tracer
!---wiso-code
    REAL(dp),    INTENT(in), OPTIONAL :: wiso_atm_spec_humidity(:,:)!! pwisoqm1
    REAL(dp),    INTENT(in), OPTIONAL :: wiso_rain(:,:)             !! rsfl + rsfc from physc
    REAL(dp),    INTENT(in), OPTIONAL :: wiso_snow(:,:)             !! ssfl + ssfc from physc
    REAL(dp),    INTENT(in), OPTIONAL :: zwisoqdif_pre(:,:,:)       !! zqdif from heat diffusion scheme
    REAL(dp),    INTENT(in), OPTIONAL :: zwisokinw(:,:)             !! kin. fractionation factor over ocean water
    REAL(dp),    INTENT(in), OPTIONAL :: zwisokinl(:,:)             !! kin. fractionation factor over land
    REAL(dp),    INTENT(in), OPTIONAL :: pwisosw_d(:,:)             !! pointer of surface ocean water delta values
!---wiso-code-end
    !! - OUTPUT --------------------------------------------------------
    REAL(dp),    INTENT(out)      :: ustarm(:)                  !! ref. to zustarm in vdiff (341)
    REAL(dp),    INTENT(out)      :: momentum_exchange_coeff(:) !! zcdum for lowest layer-soil
    REAL(dp),    INTENT(out)      :: tkevn_cond(:)              !! ptkevn
    REAL(dp),    INTENT(out)      :: ztdif_new(:)               !! ztdif
    REAL(dp),    INTENT(out)      :: zqdif_new(:)               !! zqdif
    REAL(dp),    INTENT(out)      :: ztvh(:)                    !! ztvh
    REAL(dp),    INTENT(out)      :: zqsurf(:)                  !! zqsl*zcsat*pfrl+pfrw*zqsw+pfri*zqsi
    REAL(dp),    INTENT(out)      :: zth(:)
    REAL(dp),    INTENT(out)      :: pwind10w(:)
    REAL(dp),    INTENT(out)      :: pu10(:)
    REAL(dp),    INTENT(out)      :: pv10(:)
    REAL(dp),    INTENT(out)      :: pwimax(:)
    REAL(dp),    INTENT(out)      :: pwind10(:)
    REAL(dp),    INTENT(out)      :: pdew2(:)
    REAL(dp),    INTENT(out)      :: ptemp2(:)
    REAL(dp),    INTENT(out)      :: pt2max(:)
    REAL(dp),    INTENT(out)      :: pt2min(:)
    REAL(dp),    INTENT(out)      :: pevaplac(:)
    REAL(dp),    INTENT(out)      :: pevapwac(:)
    REAL(dp),    INTENT(out)      :: pevapiac(:)
    REAL(dp),    INTENT(out)      :: pevap(:)
    REAL(dp),    INTENT(out)      :: pahfllac(:)
    REAL(dp),    INTENT(out)      :: pahflwac(:)
    REAL(dp),    INTENT(out)      :: pahfliac(:)
    REAL(dp),    INTENT(out)      :: pahfl(:)
    REAL(dp),    INTENT(out)      :: pahfslac(:)
    REAL(dp),    INTENT(out)      :: pahfswac(:)
    REAL(dp),    INTENT(out)      :: pahfsiac(:)
    REAL(dp),    INTENT(out)      :: pahfs(:)
    REAL(dp),    INTENT(out)      :: pqhfla(:)
    REAL(dp),    INTENT(out)      :: pevapw(:)
    REAL(dp),    INTENT(out)      :: pevapi(:)
    REAL(dp),    INTENT(out)      :: pahfsl(:)
    REAL(dp),    INTENT(out)      :: pahfsw(:)
    REAL(dp),    INTENT(out)      :: pahfsi(:)
    REAL(dp),    INTENT(out)      :: pevapot(:)
    REAL(dp),    INTENT(out)      :: pahflw(:)
!!$    REAL(dp),    INTENT(out)      :: pahfli(:)
    REAL(dp),    INTENT(inout)    :: psni(:)
    REAL(dp),    INTENT(out)      :: pahfice(:)
    REAL(dp),    INTENT(out)      :: pfluxres(:)
    REAL(dp),    INTENT(out)      :: pqres(:)
    REAL(dp),    INTENT(out)      :: pahfcon(:)
    REAL(dp),    INTENT(out)      :: pahfres (:)
    REAL(dp),    INTENT(out)      :: ptsi (:)
    REAL(dp),    INTENT(out)      :: ptslnew(:)
    REAL(dp),    INTENT(out)      :: pzti(:)
    REAL(dp),    INTENT(out)      :: pzteffl4(:)
    REAL(dp),    INTENT(out)      :: pztsnew(:)
    REAL(dp),    INTENT(out)      :: ptsurf(:)
    REAL(dp),    INTENT(out)      :: paz0w(:)
    REAL(dp),    INTENT(out)      :: paz0i(:)
    REAL(dp),    INTENT(out)      :: paz0l(:)
    REAL(dp),    INTENT(out)      :: paz0(:)
    REAL(dp),    INTENT(out)      :: pustrl(:)
    REAL(dp),    INTENT(out)      :: pvstrl(:)
    REAL(dp),    INTENT(out)      :: pustrw(:)
    REAL(dp),    INTENT(out)      :: pvstrw(:)
    REAL(dp),    INTENT(out)      :: pustri(:)
    REAL(dp),    INTENT(out)      :: pvstri(:)
    REAL(dp),    INTENT(out)      :: pustr(:)
    REAL(dp),    INTENT(out)      :: pvstr(:)
    REAL(dp),    INTENT(out)      :: ptte_corr(:)
    REAL(dp),    INTENT(out)      :: ptsw_new(:)
    REAL(dp),    INTENT(out)      :: pseaice_new(:)
    REAL(dp),    INTENT(out)      :: psiced_new(:)
    REAL(dp),    INTENT(out)      :: prunoff(:)
    REAL(dp),    INTENT(out)      :: pdrainage(:)
    REAL(dp),    INTENT(out)      :: pradtemp_old(:)
    REAL(dp),    INTENT(out)      :: palbedo(:)
    REAL(dp),    INTENT(out)      :: palbedo_vis(:)
    REAL(dp),    INTENT(out)      :: palbedo_nir(:)
    REAL(dp),    INTENT(out)      :: palsol(:) 
    REAL(dp),    INTENT(out)      :: palsoi(:)  
    REAL(dp),    INTENT(out)      :: palsow(:)
!!$    REAL(dp),    INTENT(out)      :: ptrfll(:)
    REAL(dp),    INTENT(out)      :: ptrflw(:)
    REAL(dp),    INTENT(out)      :: ptrfli(:) 
    REAL(dp),    INTENT(out)      :: psofll(:)
    REAL(dp),    INTENT(out)      :: psoflw(:)
    REAL(dp),    INTENT(out)      :: psofli(:)               
    REAL(dp),    INTENT(out)      :: ptrfllac(:)
    REAL(dp),    INTENT(out)      :: ptrflwac(:)
    REAL(dp),    INTENT(out)      :: ptrfliac(:)             
    REAL(dp),    INTENT(out)      :: psofllac(:)
    REAL(dp),    INTENT(out)      :: psoflwac(:)
    REAL(dp),    INTENT(out)      :: psofliac(:)  
    REAL(dp),    INTENT(out)      :: palac(:), pros_hd(:), pdrain_hd(:)
    REAL(dp),    INTENT(in)       :: pco2_flux_ocean(:)         
    REAL(dp),    INTENT(out)      :: pco2_flux_land(:)         
    REAL(dp),    INTENT(out)      :: pco2_flux(:)

!---wiso-code
    REAL(dp),    INTENT(out), OPTIONAL :: zwisoqdif_new(:,:)               !! zqdif
    REAL(dp),    INTENT(out), OPTIONAL :: pwisoevaplac(:,:)
    REAL(dp),    INTENT(out), OPTIONAL :: pwisoevapwac(:,:)
    REAL(dp),    INTENT(out), OPTIONAL :: pwisoevapiac(:,:)
    REAL(dp),    INTENT(out), OPTIONAL :: pwisoevap(:,:)
    REAL(dp),    INTENT(out), OPTIONAL :: pwisoqhfla(:,:)
    REAL(dp),    INTENT(out), OPTIONAL :: pwisoevapw(:,:)
    REAL(dp),    INTENT(out), OPTIONAL :: pwisoevapi(:,:)
    REAL(dp),    INTENT(out), OPTIONAL :: pwisoevapot(:,:)
    REAL(dp),    INTENT(out), OPTIONAL :: pwisorunoff(:,:)
    REAL(dp),    INTENT(out), OPTIONAL :: pwisodrainage(:,:)
    REAL(dp),    INTENT(out), OPTIONAL :: pwisoalac(:,:) 
    REAL(dp),    INTENT(out), OPTIONAL :: pwisoros_hd(:,:)
    REAL(dp),    INTENT(out), OPTIONAL :: pwisodrain_hd(:,:)
!---wiso-code-end
             
    !-------------------------------------------------------------------
    ! LOCAL VARIABLES
    !-------------------------------------------------------------------
    REAL(dp)       :: zqdif(current_nproma, levels)
    REAL(dp)       :: ztdif(current_nproma, levels)
    REAL(dp)       :: radtemp(current_nproma)
    REAL(dp)       :: longwave_down(current_nproma)
    REAL(dp)       :: radtemp_old(current_nproma)
    REAL(dp)       :: soil_wetness(current_nproma)
    REAL(dp)       :: snow_depth(current_nproma)
    REAL(dp)       :: zrunoff(current_nproma)
    REAL(dp)       :: zdrainage(current_nproma)
    REAL(dp)       :: skin_reservoir(current_nproma)
    REAL(dp)       :: tte_corr(current_nproma)

!---wiso-code
    REAL(dp)       :: zwisoqdif(current_nproma, levels, 1:kwiso)
    REAL(dp)       :: wiso_soil_wetness(current_nproma,1:kwiso)
    REAL(dp)       :: wiso_snow_depth(current_nproma,1:kwiso)
    REAL(dp)       :: wiso_skin_reservoir(current_nproma,1:kwiso)
    REAL(dp)       :: zwisorunoff(current_nproma,1:kwiso)
    REAL(dp)       :: zwisodrainage(current_nproma,1:kwiso)
    REAL(dp)       :: zdelta
    INTEGER        :: jt,jl,jk
!---wiso-code-end
    !-------------------------------------------------------------------
    ! MASK FOR MASKOUT LAND,OCEAN,ICE
    !-------------------------------------------------------------------
    REAL(dp) ::  zero(current_nproma)
    INTEGER  ::  nproma
!#if (! defined (__prism)) && (defined (__PGI))
    INTEGER  :: np
!#endif
    !-------------------------------------------------------------------
    ! GET CURRENT VECTOR AND ROW DEFINITION OF DOMAIN 
    !-------------------------------------------------------------------
    nproma   = current_nproma
    zero(:)  = 0._dp

    !-------------------------------------------------------------------
    ! INITIALISATION NEW START OF MODEL 
    !-------------------------------------------------------------------
    IF (lstart) THEN
       ! New start of model
       land%wind_10_meter(1:nproma,jrow)=                              &
                                       SQRT(wind_u(:)**2 + wind_v(:)**2)
!!$    ELSE
!!$       ! From restart files
!!$       land%zhsoil(1:nproma,jrow)       = hsoil(1:nproma,jrow)      
!!$       land%zcair(1:nproma,jrow)        = cair(1:nproma,jrow) 
!!$       land%zcsat(1:nproma,jrow)        = csat(1:nproma,jrow) 
    END IF

    !-------------------------------------------------------------------
    ! PERMANENT INITIALISATION 
    !-------------------------------------------------------------------
    zqdif(:,:)                                  = zqdif_pre(:,:)
    ztdif(:,:)                                  = ztdif_pre(:,:)
!---wiso-code
    IF (lwiso) THEN

    zwisoqdif(:,:,:)                            = zwisoqdif_pre(:,:,:)
    
    END IF
!---wiso-code-end
    ocean%surface_temperature(1:nproma,jrow)    = ocean_temp(:)         ! comes from clsst every time step
    box%seaice(1:nproma,jrow)                   = pseaice(:)            ! comes from clsst every time step   
    ice%ice_depth(1:nproma,jrow)                = ice_depth(:)          ! comes from clsst every time step
    ice%snow_water_equivalent(1:nproma,jrow)    = psni(:)               ! comes from ocean every time step in coupled mode (only over ice)

    ! FRACTIONAL COVER and MASK DEFINITION
    ! Note: seaice cover changes over time, therefore this isn't in init_surface
    !-------------------------------------------------------------------

#if (! defined (__prism)) && (defined (__PGI))
    ! Bugfix for the pgi compiler (version 8.0-4).
    ! The bug is somehow depending on the land sea mask. The fix is needed
    ! only in uncoupled (cosmos-as) runs.
    DO np=1,nproma
       box%ocean_fract(np,jrow)          =                          &
                                    (1._dp-box%land_fract(np,jrow)) &
                                    * (1._dp-box%seaice(np,jrow)) 
       box%seaice_fract(np,jrow)         =                          &
                                     1._dp-box%land_fract(np,jrow)  &
                                     - box%ocean_fract(np,jrow)
    END DO
#else
    box%ocean_fract(1:nproma,jrow)          =                          &
                                    (1._dp-box%land_fract(1:nproma,jrow)) &
                                    * (1._dp-box%seaice(1:nproma,jrow)) 
    box%seaice_fract(1:nproma,jrow)         =                          &
                                     1._dp-box%land_fract(1:nproma,jrow)  &
                                     - box%ocean_fract(1:nproma,jrow)

#endif
    box%frac_ice_cover_acc(1:nproma,jrow)   =                          &
                        box%frac_ice_cover_acc(1:nproma,jrow)          &
                        + delta_time * box%seaice_fract(1:nproma,jrow)
    surface%is_seaice(1:nproma,jrow)        =                          &
                                 box%seaice_fract(1:nproma,jrow) > 0.0_dp &
                           .OR. box%lake_fract(1:nproma,jrow) .GE. 0.5_dp
    surface%is_ocean (1:nproma,jrow)        =                          &
                                 box%ocean_fract (1:nproma,jrow) > 0.0_dp &
                           .OR. box%lake_fract(1:nproma,jrow) .GE. 0.5_dp

    WHERE (.NOT. surface%is_land(1:nproma,jrow))  
       ice%zcvsi(1:nproma,jrow) =                                      &
                 TANH(ice%snow_water_equivalent(1:nproma,jrow) * 100._dp)
    ELSEWHERE
       ice%zcvsi(1:nproma,jrow) = 0._dp
    END WHERE
    !-------------------------------------------------------------------
    ! additional grid points for ocean coupling
    surface%is_ocean_for_coup(1:nproma,jrow) =                         &
         slf(1:nproma,jrow) .lt. 1.0_dp                                &
         .and.(.not. surface%is_ocean(1:nproma,jrow))                  
    surface%is_ice_for_coup(1:nproma,jrow)   =                         &
         (surface%is_ocean_for_coup(1:nproma,jrow)  .or.               &
          surface%is_ocean(1:nproma,jrow)         ) .and.              &
         (.not. surface%is_seaice(1:nproma,jrow))

    !-------------------------------------------------------------------
    ! long wave radiation
    CALL longwave_down_rad( nproma                                     &
                          , longwave_net(1:nproma)                     &
                          , land%surface_temperature(1:nproma,jrow)    &
                          , ocean%surface_temperature(1:nproma,jrow)   &
                          , ice%surface_temperature(1:nproma,jrow)     &
                          , box%land_fract(1:nproma,jrow)              &
                          , box%ocean_fract(1:nproma,jrow)             &
                          , box%seaice_fract(1:nproma,jrow)            & 
                          , longwave_down(1:nproma) )
         
    box%longwave_down_acc(1:nproma,jrow)   =                           &
                        box%longwave_down_acc(1:nproma,jrow)           &
                        + delta_time * longwave_down(1:nproma)

    ! ATMOSPHERIC LOWEST LEVEL PARAMERTERS
    !-------------------------------------------------------------------
    CALL atm_conditions( nproma                                        &
                       , cloud_water(:)                                &
                       , cloud_ice(:)                                  &
                       , geopotential(:)                               &
                       , atm_temp(:)                                   &
                       , atm_spec_humidity(:)                          &
                       , atm_full_lev_press(:)                         &
                       , box%atm_tot_cloud_water(1:nproma,jrow)        &
                       , box%atm_dry_stat_energy(1:nproma,jrow)        &
                       , box%atm_pot_temp(1:nproma,jrow)               &
                       , box%atm_vir_pot_temp(1:nproma,jrow)           &
                       , box%atm_lat_heat_fact(1:nproma,jrow)          &
                       , box%atm_liq_wat_pot_temp(1:nproma,jrow)       &
                       , box%atm_sat_spec_hum(1:nproma,jrow)          )

    ! calculate surface and lowest level parameters for grid 
    ! and grid fractions of land, water, ice
    !-------------------------------------------------------------------
    CALL precalc_land( nproma                                          &
                     , surface%is_land(1:nproma,jrow)                  &  !! in
                     , jrow                                            &
                     , land%surface_temperature(1:nproma,jrow)         &
                     , atm_half_lev_press(:,levels_plus_1)             &  !! in
                     , wind_u(:)                                       &
                     , wind_v(:)                                       &  !! in
                     , atm_spec_humidity(:)                            &
                     , box%atm_tot_cloud_water(1:nproma,jrow)          &  !! in
                     , box%atm_sat_spec_hum(1:nproma,jrow)             &
                     , box%atm_pot_temp(1:nproma,jrow)                 &  !! in
                     , box%atm_vir_pot_temp(1:nproma,jrow)             &
                     , box%atm_lat_heat_fact(1:nproma,jrow)            &  !! in
                     , cloud_cover(:)                                  &
                     , box%atm_liq_wat_pot_temp(1:nproma,jrow)         &  !! in
                     , geopotential(:)                                 &
                     , land%roughness_heat(1:nproma,jrow)              &  !! in
                     , land%roughness_momentum(1:nproma,jrow)          &  !! in
                     , atm_temp(:)                                     &
                     , zghabl(:)                                       &  !! in
                     , land%ZDQSL(1:nproma,jrow)                       &
                     , land%ZRIL(1:nproma,jrow)                        &  !! out, in  
                     , land%ZQSL(1:nproma,jrow)                        &
                     , land%ZCFNCL(1:nproma,jrow)                      &  !! out, out
                     , land%ZCHL(1:nproma,jrow)                        &
                     , land%ZCFHL(1:nproma,jrow)                       &  !! out  
                     , land%zbnl(1:nproma,jrow)                        &
                     , land%zbhnl(1:nproma,jrow)                       &  !! out
                     , land%zbml(1:nproma,jrow)                        &
                     , land%zbhl(1:nproma,jrow)                        &  !! out
                     , land%zustarl(1:nproma,jrow)                     &
                     , land%ztkevl(1:nproma,jrow)                      &  !! out
                     , land%zhsoil(1:nproma,jrow)                      &
                     , land%zcfml(1:nproma,jrow)                       &  !! in, out
                     , land%zustl(1:nproma,jrow)                       &
                     , land%zcptl(1:nproma,jrow)                       &
                     , land%zcpq(1:nproma,jrow)                        &  !! out
                     , land%zcair(1:nproma,jrow)                       &
                     , land%zcsat(1:nproma,jrow)                       &
                     )

  IF (lwiso) THEN
    CALL precalc_ocean( nproma                                         &
                      , (surface%is_ocean(1:nproma,jrow) .or.          &
                      surface%is_ocean_for_coup(1:nproma,jrow))        &
                      , ocean%surface_temperature(1:nproma,jrow)       &
                      , atm_half_lev_press(:,levels_plus_1)            &
                      , atm_spec_humidity(:)                           &
                      , box%atm_tot_cloud_water(1:nproma,jrow)         &
                      , atm_temp(:)                                    &
                      , box%atm_sat_spec_hum(1:nproma,jrow)            &
                      , box%atm_pot_temp(1:nproma,jrow)                &
                      , box%atm_vir_pot_temp(1:nproma,jrow)            &
                      , box%atm_lat_heat_fact(1:nproma,jrow)           &
                      , cloud_cover(:)                                 &
                      , box%atm_liq_wat_pot_temp(1:nproma,jrow)        &
                      , wind_u(:)                                      &
                      , wind_v(:)                                      &
                      , ocean_u(:)                                     &
                      , ocean_v(:)                                     &
                      , ocean%roughness(1:nproma,jrow)                 &
                      , geopotential(:)                                &
                      , zghabl(:)                                      &
                      , ocean%zqsw(1:nproma,jrow)                      &
                      , ocean%zcptw(1:nproma,jrow)                     &
                      , ocean%zriw(1:nproma,jrow)                      &
                      , ocean%zcfhw(1:nproma,jrow)                     &
                      , ocean%zchw(1:nproma,jrow)                      &
                      , ocean%zbnw(1:nproma,jrow)                      &
                      , ocean%zbmw(1:nproma,jrow)                      &
                      , ocean%zbhw(1:nproma,jrow)                      &
                      , ocean%zustarw(1:nproma,jrow)                   &
                      , ocean%ztkevw(1:nproma,jrow)                    &
                      , ocean%zcfmw(1:nproma,jrow)                     &
                      , ocean%zustw(1:nproma,jrow)                     &
!---wiso-code
                      , lwiso, kwiso                                   &  !! in
                      , pwisosw_d(1:nproma,1:kwiso)                    &  !! in
                      , ocean%zwisoqsw(1:nproma,1:kwiso,jrow))            !! out
!---wiso-code-end

    CALL precalc_ice( nproma                                           &
                    , (surface%is_seaice(1:nproma,jrow).or.            &
                    surface%is_ice_for_coup(1:nproma,jrow))            &
                    , ice%surface_temperature(1:nproma,jrow)           &
                    , atm_half_lev_press(:,levels_plus_1)              &
                    , atm_spec_humidity(:)                             &
                    , box%atm_tot_cloud_water(1:nproma,jrow)           &
                    , atm_temp(:)                                      &
                    , box%atm_sat_spec_hum(1:nproma,jrow)              &
                    , box%atm_pot_temp(1:nproma,jrow)                  &
                    , box%atm_vir_pot_temp(1:nproma,jrow)              &
                    , box%atm_lat_heat_fact(1:nproma,jrow)             &
                    , cloud_cover(:)                                   &
                    , box%atm_liq_wat_pot_temp(1:nproma,jrow)          &
                    , geopotential(:)                                  &
                    , wind_u(:)                                        &
                    , wind_v(:)                                        &
                    , ocean_u(:)                                       &
                    , ocean_v(:)                                       &
                    , ice%roughness(1:nproma,jrow)                     &
                    , zghabl(:)                                        &
                    , ice%zqsi(1:nproma,jrow)                          &
                    , ice%zcpti(1:nproma,jrow)                         &        
                    , ice%zrii(1:nproma,jrow)                          &
                    , ice%zcfmi(1:nproma,jrow)                         & 
                    , ice%zchi(1:nproma,jrow)                          &
                    , ice%zcfhi(1:nproma,jrow)                         & 
                    , ice%zbni(1:nproma,jrow)                          &
                    , ice%zbmi(1:nproma,jrow)                          & 
                    , ice%zbhi(1:nproma,jrow)                          &
                    , ice%zustari(1:nproma,jrow)                       & 
                    , ice%ztkevi(1:nproma,jrow)                        &
                    , ice%zusti(1:nproma,jrow)                         &
!---wiso-code
                    , lwiso, kwiso                                     &  !! in
                    , pwisosw_d(1:nproma,1:kwiso)                      &  !! in
                    , ice%zwisoqsi(1:nproma,1:kwiso,jrow))                !! out
!---wiso-code-end
  ELSE
    CALL precalc_ocean( nproma                                         &
                      , (surface%is_ocean(1:nproma,jrow) .or.          &
                      surface%is_ocean_for_coup(1:nproma,jrow))        &
                      , ocean%surface_temperature(1:nproma,jrow)       &
                      , atm_half_lev_press(:,levels_plus_1)            &
                      , atm_spec_humidity(:)                           &
                      , box%atm_tot_cloud_water(1:nproma,jrow)         &
                      , atm_temp(:)                                    &
                      , box%atm_sat_spec_hum(1:nproma,jrow)            &
                      , box%atm_pot_temp(1:nproma,jrow)                &
                      , box%atm_vir_pot_temp(1:nproma,jrow)            &
                      , box%atm_lat_heat_fact(1:nproma,jrow)           &
                      , cloud_cover(:)                                 &
                      , box%atm_liq_wat_pot_temp(1:nproma,jrow)        &
                      , wind_u(:)                                      &
                      , wind_v(:)                                      &
                      , ocean_u(:)                                     &
                      , ocean_v(:)                                     &
                      , ocean%roughness(1:nproma,jrow)                 &
                      , geopotential(:)                                &
                      , zghabl(:)                                      &
                      , ocean%zqsw(1:nproma,jrow)                      &
                      , ocean%zcptw(1:nproma,jrow)                     &
                      , ocean%zriw(1:nproma,jrow)                      &
                      , ocean%zcfhw(1:nproma,jrow)                     &
                      , ocean%zchw(1:nproma,jrow)                      &
                      , ocean%zbnw(1:nproma,jrow)                      &
                      , ocean%zbmw(1:nproma,jrow)                      &
                      , ocean%zbhw(1:nproma,jrow)                      &
                      , ocean%zustarw(1:nproma,jrow)                   &
                      , ocean%ztkevw(1:nproma,jrow)                    &
                      , ocean%zcfmw(1:nproma,jrow)                     &
                      , ocean%zustw(1:nproma,jrow)                     &
!---wiso-code
                      , lwiso, kwiso)                                     !! in
!---wiso-code-end
    
    CALL precalc_ice( nproma                                           &
                    , (surface%is_seaice(1:nproma,jrow).or.            &
                    surface%is_ice_for_coup(1:nproma,jrow))            &
                    , ice%surface_temperature(1:nproma,jrow)           &
                    , atm_half_lev_press(:,levels_plus_1)              &
                    , atm_spec_humidity(:)                             &
                    , box%atm_tot_cloud_water(1:nproma,jrow)           &
                    , atm_temp(:)                                      &
                    , box%atm_sat_spec_hum(1:nproma,jrow)              &
                    , box%atm_pot_temp(1:nproma,jrow)                  &
                    , box%atm_vir_pot_temp(1:nproma,jrow)              &
                    , box%atm_lat_heat_fact(1:nproma,jrow)             &
                    , cloud_cover(:)                                   &
                    , box%atm_liq_wat_pot_temp(1:nproma,jrow)          &
                    , geopotential(:)                                  &
                    , wind_u(:)                                        &
                    , wind_v(:)                                        &
                    , ocean_u(:)                                       &
                    , ocean_v(:)                                       &
                    , ice%roughness(1:nproma,jrow)                     &
                    , zghabl(:)                                        &
                    , ice%zqsi(1:nproma,jrow)                          &
                    , ice%zcpti(1:nproma,jrow)                         &        
                    , ice%zrii(1:nproma,jrow)                          &
                    , ice%zcfmi(1:nproma,jrow)                         & 
                    , ice%zchi(1:nproma,jrow)                          &
                    , ice%zcfhi(1:nproma,jrow)                         & 
                    , ice%zbni(1:nproma,jrow)                          &
                    , ice%zbmi(1:nproma,jrow)                          & 
                    , ice%zbhi(1:nproma,jrow)                          &
                    , ice%zustari(1:nproma,jrow)                       & 
                    , ice%ztkevi(1:nproma,jrow)                        &
                    , ice%zusti(1:nproma,jrow)                         &
!---wiso-code
                    , lwiso, kwiso)                                       !! in
!---wiso-code-end
  END IF

    !-------------------------------------------------------------------
    ! CALCULATE NEW ZO VALUES

    CALL update_z0_ocean( nproma                                       &
                        , levels_plus_1                                &
                      , (surface%is_ocean(1:nproma,jrow) .or.          &
                      surface%is_ocean_for_coup(1:nproma,jrow))        &
!                        , surface%is_ocean(1:nproma, jrow)             &
                        , ocean%zcfmw(1:nproma,jrow)                   &
                        , zudif(1:nproma)                              &
                        , zvdif(1:nproma)                              &
                        , atm_temp(1:nproma)                           &
                        , atm_spec_humidity(1:nproma)                  &
                        , box%atm_tot_cloud_water(1:nproma,jrow)       &
                        , atm_half_lev_press(1:nproma,:)               &
                        , ocean%roughness(1:nproma,jrow)              )
    
    CALL update_z0_ice( nproma                                         &
                    , (surface%is_seaice(1:nproma,jrow).or.            &
                    surface%is_ice_for_coup(1:nproma,jrow))            &
!                      , surface%is_seaice(1:nproma,jrow)               &
                      , ice%roughness(1:nproma,jrow)                  )
    !-------------------------------------------------------------------
    !   RICHTMEYR-MORTON COEFFICIENTS FOR DIFFERENT SURFACE FRACTIONS
    !
  IF (lwiso) THEN
    CALL RICHTMEYR_land( nproma                                        &
                       , surface%is_land(1:nproma,jrow)                &
                       , jrow                                          &
                       , levels, levels_plus_1                         &
                       , levels_minus_1                                &
                       , atm_half_lev_press(:,:)                       &
                       , zcfh(:,:)                                     &
                       , zebsh(:,:)                                    &
                       , zqdif(:,:)                                    & 
                       , ztdif(:,:)                                    &
                       , land%zcfhl(1:nproma, jrow)                    &
                       , land%zcair(1:nproma, jrow)                    &
                       , land%zcsat(1:nproma, jrow)                    &
                       , land%zetnl(1:nproma, jrow)                    &
                       , land%zftnl(1:nproma, jrow)                    &
                       , land%zeqnl(1:nproma, jrow)                    &
                       , land%zfqnl(1:nproma, jrow)                    &
!---wiso-code
                       , lwiso, kwiso                                  &
                       ! Inputvariables:
                       , box%land_fract(1:nproma,jrow)                 &  !Hilfsvar. fuer qdif in mo_soil
                       , box%ocean_fract(1:nproma,jrow)                &  !Hilfsvar. fuer qdif in mo_soil
                       , box%seaice_fract(1:nproma,jrow)               &  !Hilfsvar. fuer qdif in mo_soil
                       , ocean%zcfhw(1:nproma, jrow)                   &  !Hilfsvar. fuer qdif in mo_soil
                       , ice%zcfhi(1:nproma, jrow)                     &  !Hilfsvar. fuer qdif in mo_soil
                       , surface%is_ocean(1:nproma,jrow)               &  !Hilfsvar. fuer qdif in mo_soil
                       , surface%is_seaice(1:nproma,jrow)              &  !Hilfsvar. fuer qdif in mo_soil
                       , zwisoqdif(:,:,1:kwiso)                        &
                       ! Outputvariables: EN and FN coefficient of the Richtmyer-Morton-Scheme, water isotopes
                       , land%zwisoeqnl(1:nproma,1:kwiso,jrow)         &
                       , land%zwisofqnl(1:nproma,1:kwiso,jrow)         &
                       , land%zwiso_nenner(1:nproma,1:kwiso,jrow)      &
                       , land%zwiso_helpqdif(1:nproma,jrow))
!---wiso-code-end

    CALL RICHTMEYR_ocean( nproma                                       &
                        , (surface%is_ocean(1:nproma,jrow) .or.        &
                           surface%is_ocean_for_coup(1:nproma,jrow))   &
                        , levels, levels_plus_1                        &
                        , levels_minus_1                               &
                        , atm_half_lev_press(:,:)                      &
                        , zcfh(:,:)                                    &
                        , zebsh(:,:)                                   &
                        , zqdif(:,:)                                   &
                        , ztdif(:,:)                                   &
                        , ocean%zcfhw(1:nproma, jrow)                  &
                        , ocean%zetnw(1:nproma, jrow)                  &
                        , ocean%zftnw(1:nproma, jrow)                  &
                        , ocean%zeqnw(1:nproma, jrow)                  &
                        , ocean%zfqnw(1:nproma, jrow)                  &
!---wiso-code
                        , lwiso, kwiso                                 &
                        ! Inputvariables: kinetic fractionaiton factor, diffusion coefficient of ???humidity???
                        , zwisokinw(:,1:kwiso)                          &
                        , zwisoqdif(:,:,1:kwiso)                       &
                        ! Outputvariables: EN and FN coefficient of the Richtmyer-Morton-Scheme, water isotopes
                        , ocean%zwisoeqnw(1:nproma,1:kwiso,jrow)       &
                        , ocean%zwisofqnw(1:nproma,1:kwiso,jrow))
!---wiso-code-end
    
    CALL RICHTMEYR_ice( nproma                                         &
                      , (surface%is_seaice(1:nproma,jrow).or.          &
                        surface%is_ice_for_coup(1:nproma,jrow))        &
                      , levels, levels_plus_1                          &
                      , levels_minus_1                                 &
                      , atm_half_lev_press(:,:)                        &
                      , zcfh(:,:)                                      &
                      , zebsh(:,:)                                     &
                      , zqdif(:,:)                                     &
                      , ztdif(:,:)                                     &
                      , ice%zcfhi(1:nproma, jrow)                      &
                      , ice%zetni(1:nproma, jrow)                      &
                      , ice%zftni(1:nproma, jrow)                      &
                      , ice%zeqni(1:nproma, jrow)                      &
                      , ice%zfqni(1:nproma, jrow)                      &
!---wiso-code
                      , lwiso, kwiso                                   &
                      ! Inputvariables: diffusion coefficient of ???humidity???
                      , zwisoqdif(:,:,1:kwiso)                         &
                      ! Outputvariables: EN and FN coefficient of the Richtmyer-Morton-Scheme, water isotopes
                      , ice%zwisoeqni(1:nproma,1:kwiso,jrow)           &
                      , ice%zwisofqni(1:nproma,1:kwiso,jrow))
!---wiso-code-end
  ELSE
    CALL RICHTMEYR_land( nproma                                        &
                       , surface%is_land(1:nproma,jrow)                &
                       , jrow                                          &
                       , levels, levels_plus_1                         &
                       , levels_minus_1                                &
                       , atm_half_lev_press(:,:)                       &
                       , zcfh(:,:)                                     &
                       , zebsh(:,:)                                    &
                       , zqdif(:,:)                                    & 
                       , ztdif(:,:)                                    &
                       , land%zcfhl(1:nproma, jrow)                    &
                       , land%zcair(1:nproma, jrow)                    &
                       , land%zcsat(1:nproma, jrow)                    &
                       , land%zetnl(1:nproma, jrow)                    &
                       , land%zftnl(1:nproma, jrow)                    &
                       , land%zeqnl(1:nproma, jrow)                    &
                       , land%zfqnl(1:nproma, jrow)                    &
!---wiso-code
                       , lwiso, kwiso)
!---wiso-code-end

    CALL RICHTMEYR_ocean( nproma                                       &
                        , (surface%is_ocean(1:nproma,jrow) .or.        &
                           surface%is_ocean_for_coup(1:nproma,jrow))   &
                        , levels, levels_plus_1                        &
                        , levels_minus_1                               &
                        , atm_half_lev_press(:,:)                      &
                        , zcfh(:,:)                                    &
                        , zebsh(:,:)                                   &
                        , zqdif(:,:)                                   &
                        , ztdif(:,:)                                   &
                        , ocean%zcfhw(1:nproma, jrow)                  &
                        , ocean%zetnw(1:nproma, jrow)                  &
                        , ocean%zftnw(1:nproma, jrow)                  &
                        , ocean%zeqnw(1:nproma, jrow)                  &
                        , ocean%zfqnw(1:nproma, jrow)                  &
!---wiso-code
                        , lwiso, kwiso)
!---wiso-code-end
    
    CALL RICHTMEYR_ice( nproma                                         &
                      , (surface%is_seaice(1:nproma,jrow).or.          &
                        surface%is_ice_for_coup(1:nproma,jrow))        &
                      , levels, levels_plus_1                          &
                      , levels_minus_1                                 &
                      , atm_half_lev_press(:,:)                        &
                      , zcfh(:,:)                                      &
                      , zebsh(:,:)                                     &
                      , zqdif(:,:)                                     &
                      , ztdif(:,:)                                     &
                      , ice%zcfhi(1:nproma, jrow)                      &
                      , ice%zetni(1:nproma, jrow)                      &
                      , ice%zftni(1:nproma, jrow)                      &
                      , ice%zeqni(1:nproma, jrow)                      &
                      , ice%zfqni(1:nproma, jrow)                      &
!---wiso-code
                      , lwiso, kwiso)
!---wiso-code-end
  END IF

    !-------------------------------------------------------------------
    ! CALCULATE NEW ALBEDO VALUES
    CALL update_albedo_ocean((surface%is_ocean(1:nproma,jrow) .or.     &
                      surface%is_ocean_for_coup(1:nproma,jrow))        &
                      !surface%is_ocean(1:nproma,jrow)          &
                            , ocean%albedo(1:nproma,jrow)             )
    CALL update_albedo_ice((surface%is_seaice(1:nproma,jrow).or.       &
                    surface%is_ice_for_coup(1:nproma,jrow))            &
                    !surface%is_seaice(1:nproma,jrow)           &
                          , ice%surface_temperature(1:nproma,jrow)     &
                          , ice%snow_water_equivalent(1:nproma,jrow)   &
                          , ice%albedo(1:nproma,jrow)                 )
    ocean%albedo_vis(1:nproma,jrow) = ocean%albedo(1:nproma,jrow)
    ocean%albedo_nir(1:nproma,jrow) = ocean%albedo(1:nproma,jrow)
    ice%albedo_vis(1:nproma,jrow) = ice%albedo(1:nproma,jrow)
    ice%albedo_nir(1:nproma,jrow) = ice%albedo(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! Calculate windstress
    CALL update_stress_land( nproma                                    &
                           , surface%is_land(1:nproma,jrow)            &
                           , land%zcfml(1:nproma,jrow)                 &
                           , zudif(:)                                  &
                           , zvdif(:)                                  &
                           , land%u_stress(1:nproma,jrow)              &
                           , land%v_stress(1:nproma,jrow)             )
    CALL update_stress_ocean( nproma                                   &
                      , (surface%is_ocean(1:nproma,jrow) .or.          &
                      surface%is_ocean_for_coup(1:nproma,jrow))        &
!                            , surface%is_ocean(1:nproma,jrow)          &
                            , ocean%zcfmw(1:nproma,jrow)               &
                            , zudif(:)                                 &
                            , zvdif(:)                                 &
                            , ocean%u_stress(1:nproma,jrow)            &
                            , ocean%v_stress(1:nproma,jrow)           )
    CALL update_stress_ice( nproma                                     &
                    , (surface%is_seaice(1:nproma,jrow).or.            &
                    surface%is_ice_for_coup(1:nproma,jrow))            &
!                          , surface%is_seaice(1:nproma,jrow)           &
                          , ice%zcfmi(1:nproma,jrow)                   &
                          , zudif(:)                                   &
                          , zvdif(:)                                   &
                          , ice%u_stress(1:nproma,jrow)                &
                          , ice%v_stress(1:nproma,jrow)               )
    !-------------------------------------------------------------------
    ! Derive snow cover fraction from snow water content
    !
    ice%snow_cover_fract(1:nproma, jrow) =                             &
                TANH(ice%snow_water_equivalent(1:nproma, jrow) * 100._dp)    
    !-------------------------------------------------------------------
    ! PRECALC FOR SURFACE LAYER ELIMINATION
    
  IF (lwiso) THEN
    CALL update_ocean( nproma                                           &
                      , (surface%is_ocean(1:nproma,jrow) .or.           &
                      surface%is_ocean_for_coup(1:nproma,jrow))         &
                      , ocean%zetnw(1:nproma,jrow)                      &
                      , ocean%zftnw(1:nproma,jrow)                      &
                      , ocean%zeqnw(1:nproma,jrow)                      &
                      , ocean%zfqnw(1:nproma,jrow)                      &
                      , ocean%zcptw(1:nproma,jrow)                      &
                      , ocean%zqsw(1:nproma,jrow)                       &
                      , ocean%ztklevw(1:nproma,jrow)                    &
                      , ocean%zqklevw(1:nproma,jrow)                    &
!---wiso-code
                      , lwiso, kwiso                                    &
                      ! Input - calculated in subroutine "Richtmeyer-Ocean"
                      , ocean%zwisoeqnw(1:nproma,1:kwiso,jrow)          &
                      , ocean%zwisofqnw(1:nproma,1:kwiso,jrow)          &
                      ! Input - calculated in subroutine "precalc_ocean"
                      , zwisokinw(1:nproma,1:kwiso)                      &
                      , ocean%zwisoqsw(1:nproma,1:kwiso,jrow)           &
                      ! Output
                      , ocean%zwisoqklevw(1:nproma,1:kwiso,jrow))
!---wiso-code-end
    
    CALL update_ice( nproma                                            &
                    , (surface%is_seaice(1:nproma,jrow).or.            &
                    surface%is_ice_for_coup(1:nproma,jrow))            & 
                    , ice%zetni(1:nproma, jrow)                        &
                    , ice%zftni(1:nproma, jrow)                        &
                    , ice%zeqni(1:nproma, jrow)                        &
                    , ice%zfqni(1:nproma, jrow)                        &
                    , ice%zcpti(1:nproma, jrow)                        &
                    , ice%zqsi(1:nproma, jrow)                         &
                    , ice%ztklevi(1:nproma,jrow)                       &
                    , ice%zqklevi(1:nproma,jrow)                       &
!---wiso-code
                    , lwiso, kwiso                                     &
                    ! Input - calculated in subroutine "Richtmeyer-Ocean"
                    , ice%zwisoeqni(1:nproma,1:kwiso,jrow)             &
                    , ice%zwisofqni(1:nproma,1:kwiso,jrow)             &
                    ! Input - calculated in subroutine "precalc_ocean"
                    , ice%zwisoqsi(1:nproma,1:kwiso,jrow)              &
                    ! Output
                    , ice%zwisoqklevi(1:nproma,1:kwiso,jrow))
!---wiso-code-end
  ELSE
    CALL update_ocean( nproma                                           &
                      , (surface%is_ocean(1:nproma,jrow) .or.           &
                      surface%is_ocean_for_coup(1:nproma,jrow))         &
                      , ocean%zetnw(1:nproma,jrow)                      &
                      , ocean%zftnw(1:nproma,jrow)                      &
                      , ocean%zeqnw(1:nproma,jrow)                      &
                      , ocean%zfqnw(1:nproma,jrow)                      &
                      , ocean%zcptw(1:nproma,jrow)                      &
                      , ocean%zqsw(1:nproma,jrow)                       &
                      , ocean%ztklevw(1:nproma,jrow)                    &
                      , ocean%zqklevw(1:nproma,jrow)                    &
!---wiso-code
                      , lwiso, kwiso)
!---wiso-code-end
    
    CALL update_ice( nproma                                            &
                    , (surface%is_seaice(1:nproma,jrow).or.            &
                    surface%is_ice_for_coup(1:nproma,jrow))            & 
                    , ice%zetni(1:nproma, jrow)                        &
                    , ice%zftni(1:nproma, jrow)                        &
                    , ice%zeqni(1:nproma, jrow)                        &
                    , ice%zfqni(1:nproma, jrow)                        &
                    , ice%zcpti(1:nproma, jrow)                        &
                    , ice%zqsi(1:nproma, jrow)                         &
                    , ice%ztklevi(1:nproma,jrow)                       &
                    , ice%zqklevi(1:nproma,jrow)                       &
!---wiso-code
                    , lwiso, kwiso)
!---wiso-code-end
  END IF

!---wiso-code
    IF (lwiso) THEN

    land%zcair_old(1:nproma,jrow) = land%zcair(1:nproma,jrow)
    
    END IF
!---wiso-code-end

  IF (lwiso) THEN
    CALL jsbach_inter_1d( nproma,                                      &
                          COUNT(surface%is_land(1:nproma,jrow)),       &
         geopot         = geopotential,                                &
         wind           = SQRT(wind_u(1:nproma)**2 + wind_v(1:nproma)**2),     &
         wind10         = land%wind_10_meter(1:nproma,jrow),           & 
         ! NOTE: This is the old 10m wind speed, we don't
         ! have the new one yet (as is the case in ECHAM)
         temp_air       = atm_temp(1:nproma),                          &
         qair           = atm_spec_humidity(1:nproma),                 &
         precip_rain    = rain,                                        &
         precip_snow    = snow,                                        &
         lwdown         = longwave_down,                               &
         sw_vis_net     = sw_vis(:),                                   &
         sw_vis_frac_diffuse = sw_vis_frac_diffuse(:),                 &
         sw_nir_net     = sw_nir(:),                                   &
         sw_nir_frac_diffuse  = sw_nir_frac_diffuse(:),                &
         pressure       = atm_half_lev_press(1:nproma,levels_plus_1),  &
         czenith        = cos_zenith_angle(1:nproma),                  &
         declination    = declination,                                 &
         CO2_concentration = pco2_concentration(1:nproma),             &
         cdrag          = land%zcfhl(1:nproma,jrow),                   &
         etAcoef        = land%zetnl(1:nproma,jrow),                   &
         etBcoef        = land%zftnl(1:nproma,jrow),                   &
         eqAcoef        = land%zeqnl(1:nproma,jrow),                   &
         eqBcoef        = land%zfqnl(1:nproma,jrow),                   &
         cair           = land%zcair(1:nproma,jrow),                   &
         csat           = land%zcsat(1:nproma,jrow),                   &
         zhsoil         = land%zhsoil(1:nproma,jrow),                  &
         albedo         = land%albedo(1:nproma,jrow),                  &
         albedo_vis     = land%albedo_vis(1:nproma,jrow),              &
         albedo_nir     = land%albedo_nir(1:nproma,jrow),              &
         z0h            = land%roughness_heat(1:nproma,jrow),          &
         z0m            = land%roughness_momentum(1:nproma,jrow),      &
         evap_act       = land%evaporation_inst(1:nproma,jrow),        &
         evap_pot       = land%evaporation_pot(1:nproma,jrow),         &
         tsoil_rad      = land%surface_temperature_rad(1:nproma,jrow), &
         temp_soil_new  = land%surface_temperature_new(1:nproma,jrow), &
         qsurf          = land%surface_qsat_new(1:nproma,jrow),        &
         mask_land      = surface%is_land(1:nproma,jrow),              &
         kblock         = jrow,                                        &
         latent         = land%latent_heat_flux_inst(1:nproma,jrow),   &
         sensible       = land%sensible_heat_flux_inst(1:nproma,jrow), &
         echam_zchl     = land%zchl(1:nproma,jrow),                    &
         surf_dry_static_energy = land%dry_static_energy_new(1:nproma,jrow), &
         soil_wetness   = soil_wetness(1:nproma),                      &
         snow_depth     = snow_depth(1:nproma),                        &
         runoff         = zrunoff(1:nproma),                           &
         drainage       = zdrainage(1:nproma),                         &
         skin_res       = skin_reservoir(1:nproma),                    & 
         tte_corr       = tte_corr(1:nproma),                          &
         glac_runoff_evap = palac(1:nproma),                           &
         surf_runoff_hd = pros_hd(1:nproma),                           &
         drainage_hd    = pdrain_hd(1:nproma),                         &
         glacier_depth = gld(1:nproma,jrow),                           &
         snow_melt_acc = snmel(1:nproma,jrow),                         &
         glacier_p_minus_e_acc = apmegl(1:nproma,jrow),                &
         snow_acc      = snacl(1:nproma,jrow),                         &
         glacier_runoff_acc = rogl(1:nproma,jrow),                     &
         wsmax = wsmx(1:nproma,jrow),                                  &
         CO2_flux       = pco2_flux_land(1:nproma),                    &
!---wiso-code
         ! Input - Isotopes
         lwisofracl         = lwisofracl,                                         &
         wiso_eqAcoef       = land%zwisoeqnl(1:nproma,1:kwiso,jrow),              &
         wiso_eqBcoef       = land%zwisofqnl(1:nproma,1:kwiso,jrow),              &
         wiso_nenner        = land%zwiso_nenner(1:nproma,1:kwiso,jrow),           &
         wiso_helpqdif      = land%zwiso_helpqdif(1:nproma,jrow),                 &  
         wisoqair           = wiso_atm_spec_humidity(1:nproma,1:kwiso),           &
         wiso_precip_rain   = wiso_rain,                                          &  
         wiso_precip_snow   = wiso_snow,                                          &  
         wisokinl           = zwisokinl(1:nproma,1:kwiso),                        &
         ! Output - Isotopes
         wisocair           = land%zwisocair(1:nproma,1:kwiso,jrow),              &
         wisocsat           = land%zwisocsat(1:nproma,1:kwiso,jrow),              &
         wisocair_fra       = land%zwisocair_fra(1:nproma,1:kwiso,jrow),          &
         wisocsat_fra       = land%zwisocsat_fra(1:nproma,1:kwiso,jrow),          &
         wiso_eqAcoef_new   = land%zwisoeqnl_new(1:nproma,1:kwiso,jrow),          &
         wiso_eqBcoef_new   = land%zwisofqnl_new(1:nproma,1:kwiso,jrow),          &
         wiso_evap_act      = land%wiso_evaporation_inst(1:nproma,1:kwiso,jrow),  &
         wiso_evap_pot      = land%wiso_evaporation_pot(1:nproma,1:kwiso,jrow),   &
         wiso_qsurf         = land%wiso_surface_qsat_new(1:nproma,1:kwiso,jrow),  &
         wiso_soil_wetness  = wiso_soil_wetness(1:nproma,1:kwiso),                &
         wiso_snow_depth    = wiso_snow_depth(1:nproma,1:kwiso),                  &
         wiso_runoff        = zwisorunoff(1:nproma,1:kwiso),                      &
         wiso_drainage      = zwisodrainage(1:nproma,1:kwiso),                    &
         wiso_skin_res      = wiso_skin_reservoir(1:nproma,1:kwiso),              &
         wiso_glac_runoff_evap = pwisoalac(1:nproma,1:kwiso),                     &
         wiso_surf_runoff_hd = pwisoros_hd(1:nproma,1:kwiso),                     &
         wiso_drainage_hd    = pwisodrain_hd(1:nproma,1:kwiso),                   &
         wiso_glacier_depth = wisogld(1:nproma,1:kwiso,jrow),                     &
         wiso_snow_melt_acc = wisosnmel(1:nproma,1:kwiso,jrow),                   &
         wiso_glacier_p_minus_e_acc = wisoapmegl(1:nproma,1:kwiso,jrow),          &
         wiso_snow_acc      = wisosnacl(1:nproma,1:kwiso,jrow),                   &
         wiso_glacier_runoff_acc = wisorogl(1:nproma,1:kwiso,jrow))
!---wiso-code-end 
  ELSE
    CALL jsbach_inter_1d( nproma,                                      &
                          COUNT(surface%is_land(1:nproma,jrow)),       &
         geopot         = geopotential,                                &
         wind           = SQRT(wind_u(1:nproma)**2 + wind_v(1:nproma)**2),     &
         wind10         = land%wind_10_meter(1:nproma,jrow),           & 
         ! NOTE: This is the old 10m wind speed, we don't
         ! have the new one yet (as is the case in ECHAM)
         temp_air       = atm_temp(1:nproma),                          &
         qair           = atm_spec_humidity(1:nproma),                 &
         precip_rain    = rain,                                        &
         precip_snow    = snow,                                        &
         lwdown         = longwave_down,                               &
         sw_vis_net     = sw_vis(:),                                   &
         sw_vis_frac_diffuse = sw_vis_frac_diffuse(:),                 &
         sw_nir_net     = sw_nir(:),                                   &
         sw_nir_frac_diffuse  = sw_nir_frac_diffuse(:),                &
         pressure       = atm_half_lev_press(1:nproma,levels_plus_1),  &
         czenith        = cos_zenith_angle(1:nproma),                  &
         declination    = declination,                                 &
         CO2_concentration = pco2_concentration(1:nproma),             &
         cdrag          = land%zcfhl(1:nproma,jrow),                   &
         etAcoef        = land%zetnl(1:nproma,jrow),                   &
         etBcoef        = land%zftnl(1:nproma,jrow),                   &
         eqAcoef        = land%zeqnl(1:nproma,jrow),                   &
         eqBcoef        = land%zfqnl(1:nproma,jrow),                   &
         cair           = land%zcair(1:nproma,jrow),                   &
         csat           = land%zcsat(1:nproma,jrow),                   &
         zhsoil         = land%zhsoil(1:nproma,jrow),                  &
         albedo         = land%albedo(1:nproma,jrow),                  &
         albedo_vis     = land%albedo_vis(1:nproma,jrow),              &
         albedo_nir     = land%albedo_nir(1:nproma,jrow),              &
         z0h            = land%roughness_heat(1:nproma,jrow),          &
         z0m            = land%roughness_momentum(1:nproma,jrow),      &
         evap_act       = land%evaporation_inst(1:nproma,jrow),        &
         evap_pot       = land%evaporation_pot(1:nproma,jrow),         &
         tsoil_rad      = land%surface_temperature_rad(1:nproma,jrow), &
         temp_soil_new  = land%surface_temperature_new(1:nproma,jrow), &
         qsurf          = land%surface_qsat_new(1:nproma,jrow),        &
         mask_land      = surface%is_land(1:nproma,jrow),              &
         kblock         = jrow,                                        &
         latent         = land%latent_heat_flux_inst(1:nproma,jrow),   &
         sensible       = land%sensible_heat_flux_inst(1:nproma,jrow), &
         echam_zchl     = land%zchl(1:nproma,jrow),                    &
         surf_dry_static_energy = land%dry_static_energy_new(1:nproma,jrow), &
         soil_wetness   = soil_wetness(1:nproma),                      &
         snow_depth     = snow_depth(1:nproma),                        &
         runoff         = zrunoff(1:nproma),                           &
         drainage       = zdrainage(1:nproma),                         &
         skin_res       = skin_reservoir(1:nproma),                    & 
         tte_corr       = tte_corr(1:nproma),                          &
         glac_runoff_evap = palac(1:nproma),                           &
         surf_runoff_hd = pros_hd(1:nproma),                           &
         drainage_hd    = pdrain_hd(1:nproma),                         &
         glacier_depth = gld(1:nproma,jrow),                           &
         snow_melt_acc = snmel(1:nproma,jrow),                         &
         glacier_p_minus_e_acc = apmegl(1:nproma,jrow),                &
         snow_acc      = snacl(1:nproma,jrow),                         &
         glacier_runoff_acc = rogl(1:nproma,jrow),                     &
         wsmax = wsmx(1:nproma,jrow),                                  &
         CO2_flux       = pco2_flux_land(1:nproma))
  END IF

    CALL update_land(                                                  &
                      nproma                                           &
                    , surface%is_land(1:nproma, jrow)                  &
                    , land%zetnl(1:nproma, jrow)                       &
                    , land%zftnl(1:nproma, jrow)                       &
                    , land%zeqnl(1:nproma, jrow)                       &
                    , land%zfqnl(1:nproma, jrow)                       &
                    , land%dry_static_energy_new(1:nproma,jrow)        &
                    , land%surface_qsat_new(1:nproma,jrow)             &
                    , land%ztklevl(1:nproma,jrow)                      &
                    , land%zqklevl(1:nproma,jrow)                      &
                    )

    !-------------------------------------------------------------------
    ! Copy land albedo to g3b memory stream
    !    alsol(1:nproma,jrow) = land%albedo(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! update ztdif and zq dif due to new surface values for humidity and 
    ! temperature at the blending height (means lowest atmosphere level)

    CALL blend_zq_zt( nproma                                           &
                    , box%land_fract(1:nproma,jrow)                    &
                    , box%ocean_fract(1:nproma,jrow)                   &
                    , box%seaice_fract(1:nproma,jrow)                  &
                    , land%zcfhl(1:nproma,jrow)                        &
                    , ocean%zcfhw(1:nproma, jrow)                      &
                    , ice%zcfhi(1:nproma, jrow)                        &
                    , land%zcair_old(1:nproma,jrow)                    &  !zcair ref. to old value (!)
                    , land%ztklevl(1:nproma,jrow)                      &
                    , ocean%ztklevw(1:nproma,jrow)                     &
                    , ice%ztklevi(1:nproma,jrow)                       &
                    , land%zqklevl(1:nproma,jrow)                      &
                    , ocean%zqklevw(1:nproma,jrow)                     &
                    , ice%zqklevi(1:nproma,jrow)                       &
                    , ztdif(:,levels)                                  &
                    , zqdif(:,levels)                                  &
                    , surface%is_land(1:nproma,jrow)                   &
                    , surface%is_ocean(1:nproma,jrow)                  &
                    , surface%is_seaice(1:nproma,jrow)                 &
                    )

!---wiso-code
    IF (lwiso) THEN

    CALL update_land_wiso(                                                &
                    kwiso                                                 &
                    , nproma                                              &
                    , surface%is_land(1:nproma, jrow)                     &
                    ! Input - calculated in subroutine "Richtmeyer-Land"
                    , land%zwisoeqnl_new(1:nproma,1:kwiso,jrow)           &
                    , land%zwisofqnl_new(1:nproma,1:kwiso,jrow)           &
                    ! Input - calculated in subroutine "jsbach_inter_1d"
                    , land%surface_qsat_new(1:nproma,jrow)                &
                    , land%wiso_surface_qsat_new(1:nproma,1:kwiso,jrow)   &
                    , land%zwisocsat(1:nproma,1:kwiso,jrow)               &
                    , land%zwisocair(1:nproma,1:kwiso,jrow)               &
                    , land%zwisocsat_fra(1:nproma,1:kwiso,jrow)           &
                    , land%zwisocair_fra(1:nproma,1:kwiso,jrow)           &
                    ! Input - needed for calc. of zqdif
                    , zqdif(1:nproma,levels)                              &
                    ! Output
                    , land%zwisoqklevl(1:nproma,1:kwiso,jrow)             &
                   )
    !-------------------------------------------------------------------
    ! Copy land albedo to g3b memory stream
    !    alsol(1:nproma,jrow) = land%albedo(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! update ztdif and zq dif due to new surface values for humidity and 
    ! temperature at the blending height (means lowest atmosphere level)

    ! first timestep land%zcair_old(1:nproma,jrow) = 0
    DO jt=1,kwiso
      IF (lstart) THEN
         land%zwisocair(:,jt,jrow) = 0._dp
         land%zwisocair_fra(:,jt,jrow) = 0._dp
      END IF
    END DO

    CALL blend_zq_zt_wiso(                                             &
                      nproma                                           &
                    , box%land_fract(1:nproma,jrow)                    &
                    , box%ocean_fract(1:nproma,jrow)                   &
                    , box%seaice_fract(1:nproma,jrow)                  &
                    , land%zcfhl(1:nproma,jrow)                        &
                    , ocean%zcfhw(1:nproma, jrow)                      &
                    , ice%zcfhi(1:nproma, jrow)                        &
                    , surface%is_land(1:nproma,jrow)                   &
                    , surface%is_ocean(1:nproma,jrow)                  &
                    , surface%is_seaice(1:nproma,jrow)                 &
                    , land%zcair_old(1:nproma,jrow)                    &  !zcair ref. to old value (!)                
                    , land%zwisocair(1:nproma,1:kwiso,jrow)            &
                    , land%zwisocair_fra(1:nproma,1:kwiso,jrow)        &
                    , land%zwisoqklevl(1:nproma,1:kwiso,jrow)          &
                    , ocean%zwisoqklevw(1:nproma,1:kwiso,jrow)         &
                    , ice%zwisoqklevi(1:nproma,1:kwiso,jrow)           &
                    , zwisoqdif(:,levels,1:kwiso)                      &
                    )

    END IF
!---wiso-code-end

    ztdif_new(1:nproma) = ztdif(1:nproma,levels)
    zqdif_new(1:nproma) = zqdif(1:nproma,levels)
!---wiso-code
    IF (lwiso) THEN

    zwisoqdif_new(1:nproma,1:kwiso) = zwisoqdif(1:nproma,levels,1:kwiso)
    
    ! For a restart without previous isotope diagnostics: 
    ! Initialize zwisoqdif_new from default JSBACH zqdif_new field
    ! (assume at start a delta value of zero permill) 
    
    IF (lresume .AND. (.NOT. lwiso_rerun)) THEN
      DO jt=1,kwiso
        zwisoqdif_new(1:nproma,jt) = tnat(jt) * zqdif_new(1:nproma)        
      END DO
    END IF
    
    END IF
!---wiso-code-end

    !-------------------------------------------------------------------

  IF (lwiso) THEN
    CALL postproc_ocean( nproma                                             &
                       , levels                                             &
                       , levels_plus_1                                      &
                       , (surface%is_ocean(1:nproma,jrow) .or.              &
                          surface%is_ocean_for_coup(1:nproma,jrow))         &
                       , ocean%zcfhw(1:nproma, jrow)                        &
                       , ocean%zqsw(1:nproma, jrow)                         &
                       , zqdif(:,levels)                                    &
                       , atm_spec_humidity(:)                               &
                       , geopotential(:)                                    &
                       , ztdif(:,levels)                                    &
                       , box%atm_dry_stat_energy(1:nproma,jrow)             &
                       , ocean%zcptw(1:nproma, jrow)                        &
                       , ocean%surface_temperature(1:nproma,jrow)           &
                       , ocean%zbnw(1:nproma, jrow)                         &
                       , ocean%zbhw(1:nproma, jrow)                         &
                       , ocean%zriw(1:nproma, jrow)                         &
                       , atm_temp(:)                                        &
                       , atm_full_lev_press(:)                              &
                       , box%atm_tot_cloud_water(1:nproma, jrow)            &
                       , wind_u(:)                                          &
                       , wind_v(:)                                          &
                       , ocean_u(:)                                         &
                       , ocean_v(:)                                         &
                       , ocean%evaporation_inst(1:nproma, jrow)             &
                       , ocean%sensible_heat_flux_inst(1:nproma, jrow)      &
                       , ocean%latent_heat_flux_inst(1:nproma, jrow)        &
                       , ocean%wind_10_meter(1:nproma, jrow)                &
                       , ocean%temp_2_meter(1:nproma, jrow)                 &
                       , atm_half_lev_press(:,:)                            &
                       , ocean%zbmw(1:nproma,jrow)                          &
                       , ocean%dewpoint_2_meter(1:nproma,jrow)              &
                       , ocean%u_wind_10_meter(1:nproma,jrow)               &
                       , ocean%v_wind_10_meter(1:nproma,jrow)               &
!---wiso-code
                       ! Input
                       , lwiso, kwiso                                       &
                       , ocean%zwisoqsw(1:nproma,1:kwiso,jrow)              &
                       , zwisoqdif(1:nproma,levels,1:kwiso)                 &
                       , wiso_atm_spec_humidity(1:nproma,1:kwiso)           &
                       , zwisokinw(1:nproma,1:kwiso)                        &
                       ! Output
                       , ocean%wiso_evaporation_inst(1:nproma,1:kwiso,jrow))
!---wiso-code-end

    CALL postproc_ice( nproma                                          &  ! kdim
                     , levels                                          &  ! klev
                     , levels_plus_1                                   &  ! klevp1
                     , (surface%is_seaice(1:nproma,jrow).or.           &  ! mask
                        surface%is_ice_for_coup(1:nproma,jrow))        &
                     , ice%zcfhi(1:nproma, jrow)                       &  ! zcfhiX
                     , ice%zqsi(1:nproma, jrow)                        &  ! zqsi
                     , zqdif(:,levels)                                 &  ! zqdif
                     , atm_spec_humidity(:)                            &  ! pqm1
                     , geopotential(:)                                 &  ! pgeom1
                     , ztdif(:,levels)                                 &  ! ztdif
                     , box%atm_dry_stat_energy(1:nproma,jrow)          &  ! zcptgz
                     , ice%zcpti (1:nproma, jrow)                      &  ! zcpti
                     , ice%surface_temperature(1:nproma,jrow)          &  ! ptsi
                     , ice%zbni(1:nproma, jrow)                        &  ! zbni
                     , ice%zbhi (1:nproma, jrow)                       &  ! zbhi
                     , ice%zrii(1:nproma, jrow)                        &  ! zrii
                     , atm_temp(:)                                     &  ! ptm1
                     , atm_full_lev_press(:)                           &  ! papm1
                     , box%atm_tot_cloud_water(1:nproma, jrow)         &  ! zx 
                     , wind_u(:)                                       &  ! pum1 
                     , wind_v(:)                                       &  ! pvm1
                     , ice%evaporation_inst(1:nproma, jrow)            &  ! pevapi
                     , ice%sensible_heat_flux_inst(1:nproma, jrow)     &  ! pahfsi
                     , ice%latent_heat_flux_inst (1:nproma, jrow)      &  ! pahfli
                     , ice%wind_10_meter(1:nproma, jrow)               &  ! zspeedi
                     , ice%temp_2_meter(1:nproma, jrow)                &  ! zt2i
                     , atm_half_lev_press(:,:)                         &  ! paphm1
                     , ice%zbmi(1:nproma, jrow)                        &  ! zbmi
                     , ice%dewpoint_2_meter(1:nproma, jrow)            &  ! zdew2i
                     , ice%u_wind_10_meter(1:nproma,jrow)              &  ! zu10i
                     , ice%v_wind_10_meter(1:nproma,jrow)              &  ! zv10i
!---wiso-code
                     ! Input
                     , lwiso, kwiso                                    &  ! lwiso, kwiso
                     , ice%zwisoqsi(1:nproma,1:kwiso,jrow)             &  ! zwisoqsi
                     , zwisoqdif(1:nproma,levels,1:kwiso)              &  ! zwisoqdif
                     , wiso_atm_spec_humidity(1:nproma,1:kwiso)        &  ! pwisoqm1
                     ! Output
                     , ice%wiso_evaporation_inst(1:nproma,1:kwiso,jrow))  ! pwisoevapi
!---wiso-code-end

  ELSE
    CALL postproc_ocean( nproma                                             &
                       , levels                                             &
                       , levels_plus_1                                      &
                       , (surface%is_ocean(1:nproma,jrow) .or.              &
                          surface%is_ocean_for_coup(1:nproma,jrow))         &
                       , ocean%zcfhw(1:nproma, jrow)                        &
                       , ocean%zqsw(1:nproma, jrow)                         &
                       , zqdif(:,levels)                                    &
                       , atm_spec_humidity(:)                               &
                       , geopotential(:)                                    &
                       , ztdif(:,levels)                                    &
                       , box%atm_dry_stat_energy(1:nproma,jrow)             &
                       , ocean%zcptw(1:nproma, jrow)                        &
                       , ocean%surface_temperature(1:nproma,jrow)           &
                       , ocean%zbnw(1:nproma, jrow)                         &
                       , ocean%zbhw(1:nproma, jrow)                         &
                       , ocean%zriw(1:nproma, jrow)                         &
                       , atm_temp(:)                                        &
                       , atm_full_lev_press(:)                              &
                       , box%atm_tot_cloud_water(1:nproma, jrow)            &
                       , wind_u(:)                                          &
                       , wind_v(:)                                          &
                       , ocean_u(:)                                         &
                       , ocean_v(:)                                         &
                       , ocean%evaporation_inst(1:nproma, jrow)             &
                       , ocean%sensible_heat_flux_inst(1:nproma, jrow)      &
                       , ocean%latent_heat_flux_inst(1:nproma, jrow)        &
                       , ocean%wind_10_meter(1:nproma, jrow)                &
                       , ocean%temp_2_meter(1:nproma, jrow)                 &
                       , atm_half_lev_press(:,:)                            &
                       , ocean%zbmw(1:nproma,jrow)                          &
                       , ocean%dewpoint_2_meter(1:nproma,jrow)              &
                       , ocean%u_wind_10_meter(1:nproma,jrow)               &
                       , ocean%v_wind_10_meter(1:nproma,jrow)               &
!---wiso-code
                       ! Input
                       , lwiso, kwiso)
!---wiso-code-end

    CALL postproc_ice( nproma                                          &  ! kdim
                     , levels                                          &  ! klev
                     , levels_plus_1                                   &  ! klevp1
                     , (surface%is_seaice(1:nproma,jrow).or.           &  ! mask
                        surface%is_ice_for_coup(1:nproma,jrow))        &
                     , ice%zcfhi(1:nproma, jrow)                       &  ! zcfhiX
                     , ice%zqsi(1:nproma, jrow)                        &  ! zqsi
                     , zqdif(:,levels)                                 &  ! zqdif
                     , atm_spec_humidity(:)                            &  ! pqm1
                     , geopotential(:)                                 &  ! pgeom1
                     , ztdif(:,levels)                                 &  ! ztdif
                     , box%atm_dry_stat_energy(1:nproma,jrow)          &  ! zcptgz
                     , ice%zcpti (1:nproma, jrow)                      &  ! zcpti
                     , ice%surface_temperature(1:nproma,jrow)          &  ! ptsi
                     , ice%zbni(1:nproma, jrow)                        &  ! zbni
                     , ice%zbhi (1:nproma, jrow)                       &  ! zbhi
                     , ice%zrii(1:nproma, jrow)                        &  ! zrii
                     , atm_temp(:)                                     &  ! ptm1
                     , atm_full_lev_press(:)                           &  ! papm1
                     , box%atm_tot_cloud_water(1:nproma, jrow)         &  ! zx 
                     , wind_u(:)                                       &  ! pum1 
                     , wind_v(:)                                       &  ! pvm1
                     , ice%evaporation_inst(1:nproma, jrow)            &  ! pevapi
                     , ice%sensible_heat_flux_inst(1:nproma, jrow)     &  ! pahfsi
                     , ice%latent_heat_flux_inst (1:nproma, jrow)      &  ! pahfli
                     , ice%wind_10_meter(1:nproma, jrow)               &  ! zspeedi
                     , ice%temp_2_meter(1:nproma, jrow)                &  ! zt2i
                     , atm_half_lev_press(:,:)                         &  ! paphm1
                     , ice%zbmi(1:nproma, jrow)                        &  ! zbmi
                     , ice%dewpoint_2_meter(1:nproma, jrow)            &  ! zdew2i
                     , ice%u_wind_10_meter(1:nproma,jrow)              &  ! zu10i
                     , ice%v_wind_10_meter(1:nproma,jrow)              &  ! zv10i
!---wiso-code
                     ! Input
                     , lwiso, kwiso)
!---wiso-code-end

  END IF

    CALL postproc_land( nproma                                         &
                      , surface%is_land(1:nproma,jrow)                 &
                      , geopotential(:)                                &
                      , land%zbnl(1:nproma,jrow)                       &
                      , land%zbml(1:nproma,jrow)                       &
                      , wind_u(:)                                      &
                      , wind_v(:)                                      &
                      , land%ZRIL(1:nproma,jrow)                       & 
                      , land%wind_10_meter(1:nproma,jrow)              &
                      , land%zbhnl(1:nproma,jrow)                      &
                      , land%zbhl(1:nproma,jrow)                       &
                      , land%zcptl(1:nproma,jrow)                      &
                      , box%atm_dry_stat_energy(1:nproma,jrow)         &
                      , atm_spec_humidity(:)                           &
                      , land%temp_2_meter(1:nproma,jrow)               &
                      , atm_temp(:)                                    &
                      , atm_full_lev_press(:)                          &
                      , atm_half_lev_press(:,:)                        &
                      , levels_plus_1                                  &
                      , box%atm_tot_cloud_water(1:nproma, jrow)        &
                      , land%dewpoint_2_meter(1:nproma, jrow)          &
                      , land%u_wind_10_meter(1:nproma,jrow)            &
                      , land%v_wind_10_meter(1:nproma,jrow)           )         
    !-------------------------------------------------------------------
    ! Average albedo
!    IF(l_trigrad) THEN
       CALL surface_box_average( nproma                                &
                               , land%albedo (1:nproma,jrow)           &
                               , ice%albedo  (1:nproma,jrow)           &
                               , ocean%albedo(1:nproma,jrow)           &
                               , box%albedo(1:nproma,jrow)             &
                               , jrow                                  &
                               , surface%is_land(1:nproma,jrow)        &
                               , surface%is_ocean(1:nproma,jrow)       &
                               , surface%is_seaice(1:nproma,jrow)      &
                               , box%land_fract(1:nproma,jrow)         &
                               , box%ocean_fract(1:nproma,jrow)        &
                               , box%seaice_fract(1:nproma,jrow)      )
       CALL surface_box_average( nproma                                &
                               , land%albedo_vis(1:nproma,jrow)        &
                               , ice%albedo_vis(1:nproma,jrow)         &
                               , ocean%albedo_vis(1:nproma,jrow)       &
                               , box%albedo_vis(1:nproma,jrow)         &
                               , jrow                                  &
                               , surface%is_land(1:nproma,jrow)        &
                               , surface%is_ocean(1:nproma,jrow)       &
                               , surface%is_seaice(1:nproma,jrow)      &
                               , box%land_fract(1:nproma,jrow)         &
                               , box%ocean_fract(1:nproma,jrow)        &
                               , box%seaice_fract(1:nproma,jrow)      )
       CALL surface_box_average( nproma                                &
                               , land%albedo_nir(1:nproma,jrow)        &
                               , ice%albedo_nir(1:nproma,jrow)         &
                               , ocean%albedo_nir(1:nproma,jrow)       &
                               , box%albedo_nir(1:nproma,jrow)         &
                               , jrow                                  &
                               , surface%is_land(1:nproma,jrow)        &
                               , surface%is_ocean(1:nproma,jrow)       &
                               , surface%is_seaice(1:nproma,jrow)      &
                               , box%land_fract(1:nproma,jrow)         &
                               , box%ocean_fract(1:nproma,jrow)        &
                               , box%seaice_fract(1:nproma,jrow)      )
!    END IF
    !-------------------------------------------------------------------
    ! calculate surface lw and sw radiation for ice and lakeice

    CALL ice_rad( nproma                                               &
                    , (surface%is_seaice(1:nproma,jrow).or.            &
                    surface%is_ice_for_coup(1:nproma,jrow))            &
!                , surface%is_seaice(1:nproma,jrow)                     &
                , longwave_down(:)                                     &
                , ice%surface_temperature(1:nproma,jrow)               &
                , pi0(:)                                               &
                , ptrsol(:)                                            &
                , ice%albedo  (1:nproma,jrow)                          &
                , box%albedo(1:nproma,jrow)                            &
                , ice%sofli(1:nproma,jrow)                             &
                , ice%trfli(1:nproma,jrow)                            )
    CALL ocean_rad( nproma                                             &
                      , (surface%is_ocean(1:nproma,jrow) .or.          &
                      surface%is_ocean_for_coup(1:nproma,jrow))        &
!                  , surface%is_ocean(1:nproma,jrow)                    &
                  , longwave_down(:)                                   &
                  , ocean%surface_temperature(1:nproma,jrow)           &
                  , pi0(:)                                             &
                  , ptrsol(:)                                          &
                  , ocean%albedo  (1:nproma,jrow)                      &
                  , box%albedo(1:nproma,jrow)                          &
                  , ocean%soflw(1:nproma,jrow)                         &
                  , ocean%trflw(1:nproma,jrow)                        )
    CALL land_rad( nproma                                              &
                 , surface%is_land(1:nproma,jrow)                      &
                 , longwave_down(:)                                    &
                 , land%surface_temperature(1:nproma,jrow)             &
                 , pi0(:)                                              &
                 , ptrsol(:)                                           &
                 , land%albedo  (1:nproma,jrow)                        &
                 , box%albedo(1:nproma,jrow)                           &
                 , land%sofll(1:nproma,jrow)                           &
                 , land%trfll(1:nproma,jrow)                           &
                 , land%surface_temperature_rad(1:nproma,jrow)         &
                 , pzteffl4(1:nproma)                                 )
    land%sofllac(1:nproma,jrow)  =  land%sofllac(1:nproma,jrow)        &
             + box%land_fract(1:nproma,jrow) * delta_time              &
             * land%sofll(1:nproma,jrow)
    land%trfllac(1:nproma,jrow)  =  land%trfllac(1:nproma,jrow)        &
             + box%land_fract(1:nproma,jrow) * delta_time              &
             * land%trfll(1:nproma,jrow)
    ocean%soflwac(1:nproma,jrow) =  ocean%soflwac(1:nproma,jrow)       &
             + box%ocean_fract(1:nproma,jrow) * delta_time             &
             * ocean%soflw(1:nproma,jrow)
    ocean%trflwac(1:nproma,jrow) =  ocean%trflwac(1:nproma,jrow)       &
             + box%ocean_fract(1:nproma,jrow) * delta_time             &
             * ocean%trflw(1:nproma,jrow)
    ice%sofliac(1:nproma,jrow)   =  ice%sofliac(1:nproma,jrow)         &
             + box%seaice_fract(1:nproma,jrow) * delta_time            &
             * ice%sofli(1:nproma,jrow)
    ice%trfliac(1:nproma,jrow)   =  ice%trfliac(1:nproma,jrow)         &
             + box%seaice_fract(1:nproma,jrow) * delta_time            &
             * ice%trfli(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! Average roughness length
    CALL surface_box_average( nproma                                   &
                            , land%roughness_momentum(1:nproma,jrow)   &
                            , ice%roughness(1:nproma,jrow)             &
                            , ocean%roughness(1:nproma,jrow)           &
                            , box%roughness(1:nproma,jrow)             &
                            , jrow                                     &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )
    !-------------------------------------------------------------------
    ! Average momentum exchange coefficient
    CALL surface_box_average( nproma                                   &
                            , land%zcfml(1:nproma,jrow)                &
                            , ice%zcfmi(1:nproma,jrow)                 &
                            , ocean%zcfmw(1:nproma,jrow)               &
                            , box%momentum_ex_coef(1:nproma,jrow)      &
                            , jrow                                     &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         ) 
    momentum_exchange_coeff(1:nproma)    =                             &
                                   box%momentum_ex_coef(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! Average ustar
    CALL average_ustar( nproma                                         &
                      , levels_plus_1                                  & !! in
                      , box%land_fract(1:nproma,jrow)                  &
                      , box%ocean_fract(1:nproma,jrow)                 & !! in
                      , box%seaice_fract(1:nproma,jrow)                &
                      , land%zustl(1:nproma,jrow)                      & !! in
                      , ocean%zustw(1:nproma,jrow)                     &
                      , ice%zusti(1:nproma,jrow)                       & !! in
                      , atm_temp(:)                                    &
                      , atm_spec_humidity(:)                           & !! in
                      , box%atm_tot_cloud_water(1:nproma,jrow)         &
                      , atm_half_lev_press(:,:)                        & !! in
                      , box%ustar(1:nproma,jrow)                       &
                      , surface%is_land(1:nproma,jrow)                 & !! out, in
                      , surface%is_ocean(1:nproma,jrow)                &
                      , surface%is_seaice(1:nproma,jrow)              )
    ustarm(1:nproma)                   = box%ustar(1:nproma,jrow) 
    !-------------------------------------------------------------------
    ! Average TKE condition
    CALL surface_box_average( nproma                                   &
                            , land%ztkevl(1:nproma,jrow)               &
                            , ice%ztkevi(1:nproma,jrow)                &
                            , ocean%ztkevw(1:nproma,jrow)              &
                            , box%ztkevn(1:nproma,jrow)                &
                            , jrow                                     & !! out, in
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )
    tkevn_cond(:)              = box%ztkevn(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! Average wind stress
    CALL surface_box_average( nproma                                   &
                            , land%u_stress(1:nproma,jrow)             &
                            , ice%u_stress(1:nproma,jrow)              &
                            , ocean%u_stress(1:nproma,jrow)            &
                            , box%u_stress(1:nproma,jrow)              &
                            , jrow                                     & !! out, in
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )
    CALL surface_box_average( nproma                                   &
                            , land%v_stress(1:nproma,jrow)             &
                            , ice%v_stress(1:nproma,jrow)              &
                            , ocean%v_stress(1:nproma,jrow)            &
                            , box%v_stress(1:nproma,jrow)              &
                            , jrow                                     & !! out, in
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )
    box%u_stress_acc(1:nproma,jrow) =  box%u_stress_acc(1:nproma,jrow) &
                          + box%u_stress(1:nproma,jrow) * delta_time
    box%v_stress_acc(1:nproma,jrow) =  box%v_stress_acc(1:nproma,jrow) &
                          + box%v_stress(1:nproma,jrow) * delta_time
    !-------------------------------------------------------------------
    ! Average 10 meter wind speed
    CALL surface_box_average( nproma                                   &
                            , land%wind_10_meter(1:nproma,jrow)        &
                            , ice%wind_10_meter(1:nproma,jrow)         &
                            , ocean%wind_10_meter(1:nproma,jrow)       &
                            , box%wind_10_meter(1:nproma,jrow)         &
                            , jrow                                     & !! out, in
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )
    CALL surface_box_average( nproma                                   &
                            , land%u_wind_10_meter(1:nproma,jrow)      &
                            , ice%u_wind_10_meter(1:nproma,jrow)       &
                            , ocean%u_wind_10_meter(1:nproma,jrow)     &
                            , box%u_wind_10_meter(1:nproma,jrow)       &
                            , jrow                                     & !! out, in
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )
    CALL surface_box_average( nproma                                   &
                            , land%v_wind_10_meter(1:nproma,jrow)      &
                            , ice%v_wind_10_meter(1:nproma,jrow)       &
                            , ocean%v_wind_10_meter(1:nproma,jrow)     &
                            , box%v_wind_10_meter(1:nproma,jrow)       &
                            , jrow                                     & !! out, in
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )

    box%wind_10_meter_acc(1:nproma,jrow)    =                          &
                                  box%wind_10_meter_acc(1:nproma,jrow) &
                       + box%wind_10_meter(1:nproma,jrow) *delta_time
    box%wind_10_max(1:nproma,jrow)          =                          &
                                  MAX(box%wind_10_max(1:nproma,jrow) , &
                        box%wind_10_meter(1:nproma,jrow) )
    !-------------------------------------------------------------------
    ! Average moisture flux (evaporation)
    CALL surface_box_average( nproma                                   &
                            , land%evaporation_inst(1:nproma,jrow)     &
                            , ice%evaporation_inst(1:nproma,jrow)      &
                            , ocean%evaporation_inst(1:nproma,jrow)    &
                            , box%evaporation_inst(1:nproma,jrow)      &
                            , jrow                                     &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )

    ice%evaporation_acc(1:nproma,jrow) =                               &
                                 ice%evaporation_acc(1:nproma,jrow)    &
                       + MERGE((ice%evaporation_inst(1:nproma,jrow)    &
                       * box%seaice_fract(1:nproma,jrow)), zero,       &
                       surface%is_seaice(1:nproma,jrow)) * delta_time
    ocean%evaporation_acc(1:nproma,jrow) =                             &
                               ocean%evaporation_acc(1:nproma,jrow)    &
                     + MERGE((ocean%evaporation_inst(1:nproma,jrow)    &
                     * box%ocean_fract(1:nproma,jrow)), zero,          &
                     surface%is_ocean(1:nproma,jrow)) * delta_time
    land%evaporation_acc(1:nproma,jrow) =                              &
                               land%evaporation_acc(1:nproma,jrow)     &
                     + MERGE((land%evaporation_inst(1:nproma,jrow)     &
                     * box%land_fract(1:nproma,jrow)), zero,           &
                     surface%is_land(1:nproma,jrow)) *delta_time         
    box%evaporation_acc(1:nproma,jrow) =                               &
                               box%evaporation_acc(1:nproma,jrow)      &
                   + box%evaporation_inst(1:nproma,jrow) * delta_time

!---wiso-code
! Average moisture flux (evaporation) - Isotopes
    IF (lwiso) THEN
   
    CALL surface_box_avg_wiso( nproma                                                &
                             , land%wiso_evaporation_inst(1:nproma,1:kwiso,jrow)     &
                             , ice%wiso_evaporation_inst(1:nproma,1:kwiso,jrow)      &
                             , ocean%wiso_evaporation_inst(1:nproma,1:kwiso,jrow)    &
                             , box%wiso_evaporation_inst(1:nproma,1:kwiso,jrow)      &
                             , jrow, kwiso                                           &
                             , surface%is_land(1:nproma,jrow)                        &
                             , surface%is_ocean(1:nproma,jrow)                       &
                             , surface%is_seaice(1:nproma,jrow)                      &
                             , box%land_fract(1:nproma,jrow)                         &
                             , box%ocean_fract(1:nproma,jrow)                        &
                             , box%seaice_fract(1:nproma,jrow)                       )

    DO jt=1,kwiso
        ice%wiso_evaporation_acc(1:nproma,jt,jrow) = ice%wiso_evaporation_acc(1:nproma,jt,jrow)             &
                                                   + MERGE((ice%wiso_evaporation_inst(1:nproma,jt,jrow)     &
                                                   * box%seaice_fract(1:nproma,jrow)), zero,                &
                                                     surface%is_seaice(1:nproma,jrow)) * delta_time
        ocean%wiso_evaporation_acc(1:nproma,jt,jrow) = ocean%wiso_evaporation_acc(1:nproma,jt,jrow)         &
                                                     + MERGE((ocean%wiso_evaporation_inst(1:nproma,jt,jrow) &
                                                     * box%ocean_fract(1:nproma,jrow)), zero,               &
                                                       surface%is_ocean(1:nproma,jrow)) * delta_time
        land%wiso_evaporation_acc(1:nproma,jt,jrow) = land%wiso_evaporation_acc(1:nproma,jt,jrow)           &
                                                    + MERGE((land%wiso_evaporation_inst(1:nproma,jt,jrow)   &
                                                    * box%land_fract(1:nproma,jrow)), zero,                 &
                                                      surface%is_land(1:nproma,jrow)) *delta_time
        box%wiso_evaporation_acc(1:nproma,jt,jrow) = box%wiso_evaporation_acc(1:nproma,jt,jrow)             &
                                                   + box%wiso_evaporation_inst(1:nproma,jt,jrow) * delta_time
    END DO
    
    END IF
!---wiso-code-end

    !-------------------------------------------------------------------
    ! Average sensible heat flux
    CALL surface_box_average( nproma                                   &
                      , land%sensible_heat_flux_inst(1:nproma,jrow)    &
                      , ice%sensible_heat_flux_inst(1:nproma,jrow)     &
                      , ocean%sensible_heat_flux_inst(1:nproma,jrow)   &
                      , box%sensible_heat_flux_inst(1:nproma,jrow)     &
                      , jrow                                           &
                      , surface%is_land(1:nproma,jrow)                 &
                      , surface%is_ocean(1:nproma,jrow)                &
                      , surface%is_seaice(1:nproma,jrow)               &
                      , box%land_fract(1:nproma,jrow)                  &
                      , box%ocean_fract(1:nproma,jrow)                 &
                      , box%seaice_fract(1:nproma,jrow)               )

    land%sensible_flux_acc(1:nproma,jrow) =                            &
                                 land%sensible_flux_acc(1:nproma,jrow) & 
                + MERGE(land%sensible_heat_flux_inst(1:nproma,jrow),   &
                zero, surface%is_land(1:nproma,jrow))                  &
                * box%land_fract(1:nproma,jrow) * delta_time
    ocean%sensible_flux_acc(1:nproma,jrow) =                           &
                                ocean%sensible_flux_acc(1:nproma,jrow) & 
                + MERGE(ocean%sensible_heat_flux_inst(1:nproma,jrow),  &
                zero, surface%is_ocean(1:nproma,jrow))                 &
                * box%ocean_fract(1:nproma,jrow) * delta_time
    ice%sensible_flux_acc(1:nproma,jrow) =                             &
                                  ice%sensible_flux_acc(1:nproma,jrow) &
                + MERGE(ice%sensible_heat_flux_inst(1:nproma,jrow),    &
                zero,surface%is_seaice(1:nproma,jrow))                 &
                * box%seaice_fract(1:nproma,jrow) * delta_time
    box%sensible_heat_flux_acc(1:nproma,jrow) =                        &
               box%sensible_heat_flux_acc(1:nproma,jrow) +             &
               box%sensible_heat_flux_inst(1:nproma,jrow) * delta_time
    !-------------------------------------------------------------------
    ! Average latent heat flux
    CALL  surface_box_average( nproma                                  &
                      , land%latent_heat_flux_inst(1:nproma,jrow)      &
                      , ice%latent_heat_flux_inst(1:nproma,jrow)       &
                      , ocean%latent_heat_flux_inst(1:nproma,jrow)     &
                      , box%latent_heat_flux_inst(1:nproma,jrow)       &
                      , jrow                                           &
                      , surface%is_land(1:nproma,jrow)                 &
                      , surface%is_ocean(1:nproma,jrow)                &
                      , surface%is_seaice(1:nproma,jrow)               &
                      , box%land_fract(1:nproma,jrow)                  &
                      , box%ocean_fract(1:nproma,jrow)                 &
                      , box%seaice_fract(1:nproma,jrow)               )

    ice%latent_heat_flux_acc(1:nproma,jrow)     =                      &
                               ice%latent_heat_flux_acc(1:nproma,jrow) &
                     + MERGE((ice%latent_heat_flux_inst(1:nproma,jrow) &
                     * box%seaice_fract(1:nproma,jrow)), zero,         &
                     surface%is_seaice(1:nproma,jrow)) * delta_time
    ocean%latent_heat_flux_acc(1:nproma,jrow)   =                      &
                             ocean%latent_heat_flux_acc(1:nproma,jrow) &
                   + MERGE((ocean%latent_heat_flux_inst(1:nproma,jrow) &
                   * box%ocean_fract(1:nproma,jrow)), zero,            &
                   surface%is_ocean(1:nproma,jrow)) * delta_time
    land%latent_heat_flux_acc(1:nproma,jrow)    =                      &
                              land%latent_heat_flux_acc(1:nproma,jrow) & 
                    + MERGE((land%latent_heat_flux_inst(1:nproma,jrow) &
                    * box%land_fract(1:nproma,jrow)), zero,            &
                    surface%is_land(1:nproma,jrow)) * delta_time
    box%latent_heat_flux_acc(1:nproma,jrow) =                          &
                             box%latent_heat_flux_acc(1:nproma,jrow) + &
                 box%latent_heat_flux_inst(1:nproma,jrow) * delta_time
    !-------------------------------------------------------------------
    ! Average 2 meter temperature
    CALL  surface_box_average( nproma                                  &
                             , land%temp_2_meter(1:nproma,jrow)        &
                             , ice%temp_2_meter(1:nproma,jrow)         &
                             , ocean%temp_2_meter(1:nproma,jrow)       &
                             , box%temp_2_meter(1:nproma,jrow)         &
                             , jrow                                    &
                             , surface%is_land(1:nproma,jrow)          &
                             , surface%is_ocean(1:nproma,jrow)         &
                             , surface%is_seaice(1:nproma,jrow)        &
                             , box%land_fract(1:nproma,jrow)           &
                             , box%ocean_fract(1:nproma,jrow)          &
                             , box%seaice_fract(1:nproma,jrow)        )
    box%maxtemp_2_meter(1:nproma,jrow) =                               &
                           MAX( box%maxtemp_2_meter(1:nproma,jrow),    &
                                box%temp_2_meter(1:nproma,jrow)) 
    box%mintemp_2_meter(1:nproma,jrow) =                               &
                           MIN( box%mintemp_2_meter(1:nproma,jrow),    &
                                box%temp_2_meter(1:nproma,jrow)) 
    !-------------------------------------------------------------------
    ! Average 2 meter dew point
    CALL  surface_box_average( nproma                                  &
                             , land%dewpoint_2_meter(1:nproma,jrow)    &
                             , ice%dewpoint_2_meter(1:nproma,jrow)     &
                             , ocean%dewpoint_2_meter(1:nproma,jrow)   &
                             , box%dewpoint_2_meter(1:nproma,jrow)     &
                             , jrow                                    &
                             , surface%is_land(1:nproma,jrow)          &
                             , surface%is_ocean(1:nproma,jrow)         &
                             , surface%is_seaice(1:nproma,jrow)        &
                             , box%land_fract(1:nproma,jrow)           &
                             , box%ocean_fract(1:nproma,jrow)          &
                             , box%seaice_fract(1:nproma,jrow)        )
    !-------------------------------------------------------------------
    ! Average virtual surface temperature and surface humidity
    CALL  average_tvh_qsurf( nproma                                    &
                           , land%surface_temperature(1:nproma,jrow)   &
                           , land%zhsoil(1:nproma,jrow)                &
                           , ocean%surface_temperature(1:nproma,jrow)  &
                           , ice%surface_temperature(1:nproma,jrow)    &
                           , box%land_fract(1:nproma,jrow)             &
                           , box%ocean_fract(1:nproma,jrow)            &
                           , box%seaice_fract(1:nproma,jrow)           &
                           , land%zcsat(1:nproma, jrow)                &
                           , ocean%zqsw(1:nproma,jrow)                 &
                           , ice%zqsi(1:nproma,jrow)                   &
                           , box%surf_vir_temp(1:nproma,jrow)          &
                           , box%surface_humidity(1:nproma,jrow)       &
                           , land%zqsl(1:nproma,jrow)                  &
                           , surface%is_land(1:nproma,jrow)            &
                           , surface%is_ocean(1:nproma,jrow)           &
                           , surface%is_seaice(1:nproma,jrow)         )
    !-------------------------------------------------------------------
    ! LAKE PHYSICS (old lake routine, 
    ! up to now located in mo_surface_ocean, mo_surface_ice)

    CALL s_lake( nproma                                                &
               , box%seaice(1:nproma,jrow)                             &
               , ice%ice_depth(1:nproma,jrow)                          &
               , box%lake_fract(1:nproma,jrow)                         &
               , ice%surface_temperature(1:nproma,jrow)                &
               , ocean%surface_temperature(1:nproma,jrow)              &
               , ocean%latent_heat_flux_inst(1:nproma,jrow)            &
               , ocean%sensible_heat_flux_inst(1:nproma,jrow)          &
               , ocean%fluxres(1:nproma,jrow)                          &
               , ocean%trflw(1:nproma,jrow)                            &
               , ocean%soflw(1:nproma,jrow)                            &
               , ice%evaporation_inst(1:nproma,jrow)                   &
               , ice%snow_water_equivalent(1:nproma,jrow)              &
               , ice%zcvsi(1:nproma,jrow)                              &
               , ice%melting(1:nproma,jrow)                            &
               , box%seaice_fract(1:nproma,jrow)                      )

    CALL s_licetemp( nproma                                            &
                   , ice%ice_depth(1:nproma,jrow)                      &
                   , ice%snow_water_equivalent(1:nproma,jrow)          &
                   , box%lake_fract(1:nproma,jrow)                     &
                   , ice%surface_temperature(1:nproma,jrow)            &
                   , ice%trfli(1:nproma,jrow)                          &
                   , ice%sofli(1:nproma,jrow)                          &
                   , ice%ahfice(1:nproma,jrow)                         &
                   , ocean%fluxres(1:nproma,jrow)                      &
                   , ice%ahfcon(1:nproma,jrow)                         &
                   , ice%melting(1:nproma,jrow)                        &
                   , ice%evaporation_inst(1:nproma,jrow)               &
                   , snow(:)                                           &
                   , ice%sensible_heat_flux_inst(1:nproma,jrow)        &
                   , ice%latent_heat_flux_inst(1:nproma,jrow)          &
                   , ice%zcvsi(1:nproma,jrow)                          &
                   , box%seaice_fract(1:nproma,jrow)                  )

    CALL s_sicetemp( nproma                                            &
                   , ice%ice_depth(1:nproma,jrow)                      &
                   , ice%snow_water_equivalent(1:nproma,jrow)          &
                   , box%lake_fract(1:nproma,jrow)                     &
                   , slf(1:nproma,jrow)                                &
                   , ice%surface_temperature(1:nproma,jrow)            &
                   , ice%trfli(1:nproma,jrow)                          &
                   , ice%sofli(1:nproma,jrow)                          &
                   , ice%ahfice(1:nproma,jrow)                         &
                   , ocean%fluxres(1:nproma,jrow)                      &
                   , ice%qres(1:nproma,jrow)                           &
                   , ice%ahfcon(1:nproma,jrow)                         &
                   , ice%melting(1:nproma,jrow)                        &
                   , ice%sensible_heat_flux_inst(1:nproma,jrow)        &
                   , ice%latent_heat_flux_inst(1:nproma,jrow)          &
                   , box%seaice_fract(1:nproma,jrow)                   &
                    , (surface%is_seaice(1:nproma,jrow).or.            &
                    surface%is_ice_for_coup(1:nproma,jrow))            &
!                   , surface%is_seaice(1:nproma,jrow)                 
                    )

    !-------------------------------------------------------------------
    ! Average surface temperature
    CALL surface_box_average( nproma                                   &
                       , land%surface_temperature_new(1:nproma,jrow)   &
                       , ice%surface_temperature(1:nproma,jrow)        &
                       , ocean%surface_temperature(1:nproma,jrow)      &
                       , box%surface_temperature(1:nproma,jrow)        &
                       , jrow                                          & !! out, in
                       , surface%is_land(1:nproma,jrow)                &
                       , surface%is_ocean(1:nproma,jrow)               &
                       , surface%is_seaice(1:nproma,jrow)              &
                       , box%land_fract(1:nproma,jrow)                 &
                       , box%ocean_fract(1:nproma,jrow)                &
                       , box%seaice_fract(1:nproma,jrow)              )
    box%surface_temperature_acc(1:nproma,jrow) =                       &
                  box%surface_temperature_acc(1:nproma,jrow)           &
                + box%surface_temperature(1:nproma,jrow) * delta_time
    !-------------------------------------------------------------------
    ! Radiative Temperature for subroutine RADIATION (from physc)
    CALL surface_box_average( nproma                                   &
                    , (land%surface_temperature_new(1:nproma,jrow)**4) &
                    , (ice%surface_temperature(1:nproma,jrow)**4)      &
                    , (ocean%surface_temperature(1:nproma,jrow)**4)    &
                    , radtemp(1:nproma)                                &
                    , jrow                                             & !! out, in
                    , surface%is_land(1:nproma,jrow)                   &
                    , surface%is_ocean(1:nproma,jrow)                  &
                    , surface%is_seaice(1:nproma,jrow)                 &
                    , box%land_fract(1:nproma,jrow)                    &
                    , box%ocean_fract(1:nproma,jrow)                   &
                    , box%seaice_fract(1:nproma,jrow)                 )
    !-------------------------------------------------------------------
    ! Radiative Temperature for subroutine RADHEAT (from physc)
    CALL surface_box_average( nproma                                   &
                       , (land%surface_temperature(1:nproma,jrow)**4)  &
                       , (ice%surface_temperature(1:nproma,jrow)**4)   &
                       , (ocean%surface_temperature(1:nproma,jrow)**4) &
                       , radtemp_old(1:nproma)                         &
                       , jrow                                          & !! out, in
                       , surface%is_land(1:nproma,jrow)                &
                       , surface%is_ocean(1:nproma,jrow)               &
                       , surface%is_seaice(1:nproma,jrow)              &
                       , box%land_fract(1:nproma,jrow)                 &
                       , box%ocean_fract(1:nproma,jrow)                &
                       , box%seaice_fract(1:nproma,jrow)              )
    land%surface_temperature(1:nproma,jrow)              =             &
                            land%surface_temperature_new(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! MASK declaration for specific surface part
    ! (pingo masks for maskout)
    box%landpingo                                        = 0._dp
    box%oceanpingo                                       = 0._dp
    box%icepingo                                         = 0._dp
    WHERE(surface%is_land) box%landpingo                 = 1._dp
    WHERE(surface%is_ocean .or. surface%is_ocean_for_coup) box%oceanpingo = 1._dp
    WHERE(surface%is_seaice .or. surface%is_ice_for_coup) box%icepingo   = 1._dp
    !-------------------------------------------------------------------
    ! Variables for ECHAM (updates Parameters and old stream vars ect.
    zth(1:nproma)    = (radtemp(1:nproma) )**0.25_dp
    pradtemp_old(1:nproma) = (radtemp_old(1:nproma))**0.25_dp
    ztvh(1:nproma)   = box%surf_vir_temp(1:nproma,jrow)
    zqsurf(1:nproma) = box%surface_humidity(1:nproma,jrow) 
    ptsw_new(:)      = ocean%surface_temperature(1:nproma,jrow)
    palsol(:)        = land%albedo (1:nproma,jrow)
    palsoi(:)        = ice%albedo  (1:nproma,jrow)
    palsow(:)        = ocean%albedo(1:nproma,jrow)
    palbedo(:)       = box%albedo(1:nproma,jrow)
    palbedo_vis(:)   = box%albedo_vis(1:nproma,jrow)
    palbedo_nir(:)   = box%albedo_nir(1:nproma,jrow)
!!$    ptrfll(:)        = land%trfll(1:nproma,jrow)
    ptrflw(:)        = ocean%trflw(1:nproma,jrow)
    ptrfli(:)        = ice%trfli(1:nproma,jrow)
    psofll(:)        = land%sofll(1:nproma,jrow)
    psoflw(:)        = ocean%soflw(1:nproma,jrow)
    psofli(:)        = ice%sofli(1:nproma,jrow)   
    ptrfllac(:)      = land%trfllac(1:nproma,jrow)
    ptrflwac(:)      = ocean%trflwac(1:nproma,jrow)
    ptrfliac(:)      = ice%trfliac(1:nproma,jrow)
    psofllac(:)      = land%sofllac(1:nproma,jrow)
    psoflwac(:)      = ocean%soflwac(1:nproma,jrow)
    psofliac(:)      = ice%sofliac(1:nproma,jrow)
    pwind10w(:)      = ocean%wind_10_meter(1:nproma,jrow)
    pu10(:)          = box%u_wind_10_meter(1:nproma,jrow)   
    pv10(:)          = box%v_wind_10_meter(1:nproma,jrow)
    pwimax(:)        = box%wind_10_max(1:nproma,jrow)
    pwind10(:)       = box%wind_10_meter_acc(1:nproma,jrow)
    pdew2(:)         = box%dewpoint_2_meter(1:nproma,jrow)    
    ptemp2(:)        = box%temp_2_meter(1:nproma,jrow) 
    pt2max(:)        = box%maxtemp_2_meter(1:nproma,jrow)
    pt2min(:)        = box%mintemp_2_meter(1:nproma,jrow)
    pevaplac(:)      = land%evaporation_acc(1:nproma,jrow)
    pevapwac(:)      = ocean%evaporation_acc(1:nproma,jrow) 
    pevapiac(:)      = ice%evaporation_acc(1:nproma,jrow)
    pevap(:)         = box%evaporation_acc(1:nproma,jrow)
    pahfllac(:)      = land%latent_heat_flux_acc(1:nproma,jrow)  
    pahflwac(:)      = ocean%latent_heat_flux_acc(1:nproma,jrow)
    pahfliac(:)      = ice%latent_heat_flux_acc(1:nproma,jrow)   
    pahfl(:)         = box%latent_heat_flux_acc(1:nproma,jrow)
    pahfslac(:)      = land%sensible_flux_acc(1:nproma,jrow) 
    pahfswac(:)      = ocean%sensible_flux_acc(1:nproma,jrow) 
    pahfsiac(:)      = ice%sensible_flux_acc(1:nproma,jrow) 
    pahfs(:)         = box%sensible_heat_flux_acc(1:nproma,jrow) 
    pqhfla(:)        = box%evaporation_inst(1:nproma,jrow)  
!!$    pevapl(:)        = land%evaporation_inst(1:nproma,jrow)
    pevapw(:)        = ocean%evaporation_inst(1:nproma,jrow)  
    pevapi(:)        = ice%evaporation_inst(1:nproma,jrow) 
    pahfsl(:)        = land%sensible_heat_flux_inst(1:nproma,jrow) 
    pahfsw(:)        = ocean%sensible_heat_flux_inst(1:nproma,jrow)
    pahfsi(:)        = ice%sensible_heat_flux_inst(1:nproma,jrow)
    pevapot(:)       = land%evaporation_pot(1:nproma,jrow)
    pahflw(:)        = ocean%latent_heat_flux_inst(1:nproma,jrow)
!!$    pahfli(:)        = ice%latent_heat_flux_inst(1:nproma,jrow)
    psni(:)          = ice%snow_water_equivalent(1:nproma,jrow)
    pahfice(:)       = ice%ahfice(1:nproma,jrow)
    pfluxres(:)      = ocean%fluxres(1:nproma,jrow)
    pqres(:)         = ice%qres(1:nproma,jrow)
    pahfcon(:)       = ice%ahfcon(1:nproma,jrow)  
    pahfres(:)       = ice%melting(1:nproma,jrow)
    ptsi(:)          = ice%surface_temperature(1:nproma,jrow)
    ptslnew(:)       = land%surface_temperature_rad(1:nproma,jrow)
    pzti(:)          = radtemp(1:nproma)**0.25_dp
    pztsnew(:)       = ( box%land_fract  (1:nproma,jrow) * pzteffl4(1:nproma) + &
                         box%seaice_fract(1:nproma,jrow) * ice%surface_temperature(1:nproma,jrow)**4 + &
                         box%ocean_fract (1:nproma,jrow) * ocean%surface_temperature(1:nproma,jrow)**4 ) &
                         **0.25_dp
    ptsurf(:)        = box%surface_temperature_acc(1:nproma,jrow)
    paz0w(:)         = ocean%roughness(1:nproma,jrow)
    paz0i(:)         = ice%roughness(1:nproma,jrow)
    paz0l(:)         = land%roughness_momentum(1:nproma,jrow)
    paz0(:)          = box%roughness(1:nproma,jrow) 
    pustrl(:)        = land%u_stress(1:nproma,jrow)
    pvstrl(:)        = land%v_stress(1:nproma,jrow) 
    pustrw(:)        = ocean%u_stress(1:nproma,jrow)
    pvstrw(:)        = ocean%v_stress(1:nproma,jrow) 
    pustri(:)        = ice%u_stress(1:nproma,jrow)
    pvstri(:)        = ice%v_stress(1:nproma,jrow)
    prunoff(:)       = zrunoff(1:nproma)
    pdrainage(:)     = zdrainage(1:nproma)
    pustr(:)         = box%u_stress_acc(1:nproma,jrow)
    pvstr(:)         = box%v_stress_acc(1:nproma,jrow)
    ptte_corr(:)     = tte_corr(:)/                                    &
                        (atm_half_lev_press(:,levels_plus_1)           &
                        -atm_half_lev_press(:,levels))
    pseaice_new(:)   = box%seaice(1:nproma,jrow)
    psiced_new(:)    = ice%ice_depth(1:nproma,jrow)
    pco2_flux(:)     = MERGE(box%land_fract (1:nproma,jrow) * pco2_flux_land(1:nproma), &
                             zero, surface%is_land (1:nproma,jrow)) + &
                       MERGE((box%ocean_fract(1:nproma,jrow)+box%seaice_fract(1:nproma,jrow)) * pco2_flux_ocean(1:nproma),  &
                             zero, surface%is_ocean(1:nproma,jrow))

!---wiso-code
    IF (lwiso) THEN

    pwisoevaplac(:,:)   = land%wiso_evaporation_acc(1:nproma,1:kwiso,jrow)
    pwisoevapwac(:,:)   = ocean%wiso_evaporation_acc(1:nproma,1:kwiso,jrow)
    pwisoevapiac(:,:)   = ice%wiso_evaporation_acc(1:nproma,1:kwiso,jrow)
    pwisoevap(:,:)      = box%wiso_evaporation_acc(1:nproma,1:kwiso,jrow) 
    pwisoqhfla(:,:)     = box%wiso_evaporation_inst(1:nproma,1:kwiso,jrow) 
    pwisoevapw(:,:)     = ocean%wiso_evaporation_inst(1:nproma,1:kwiso,jrow)  
    pwisoevapi(:,:)     = ice%wiso_evaporation_inst(1:nproma,1:kwiso,jrow) 
    pwisoevapot(:,:)    = land%wiso_evaporation_pot(1:nproma,1:kwiso,jrow)
    pwisorunoff(:,:)    = zwisorunoff(1:nproma,1:kwiso)
    pwisodrainage(:,:)  = zwisodrainage(1:nproma,1:kwiso)
    
    END IF
!---wiso-code-end

    !-------------------------------------------------------------------
    ! Variables for g3b stream only for diagnose !
    friac(1:nproma,jrow)    = box%frac_ice_cover_acc(1:nproma,jrow)
    tslm1(1:nproma,jrow)    = land%surface_temperature_new(1:nproma,jrow) 
    ws(1:nproma,jrow)       = soil_wetness(1:nproma)
    sn(1:nproma,jrow)       = snow_depth(1:nproma)
    wl(1:nproma,jrow)       = skin_reservoir(1:nproma)

!---wiso-code
    IF (lwiso) THEN

    DO jt=1,kwiso
       wisows(1:nproma,jt,jrow)     = wiso_soil_wetness(1:nproma,jt)
       wisosn(1:nproma,jt,jrow)     = wiso_snow_depth(1:nproma,jt)
       wisowl(1:nproma,jt,jrow)     = wiso_skin_reservoir(1:nproma,jt)
    END DO
    
    END IF
!---wiso-code-end

  END SUBROUTINE update_surface

  SUBROUTINE init_surface

    USE mo_surface_ocean, ONLY: update_albedo_ocean
    USE mo_surface_ice, ONLY: init_albedo_ice, update_albedo_ice
    USE mo_surface_boundary, ONLY: surface_box_average
    USE mo_memory_g3b, ONLY: &
         slm,                &   !! Sea-land mask, currently read in ioinitial (contains 1 or 0)
         slf,                &
         alake,              &    !! lake faction 
         alsol, alsoi, alsow, albedo, albedo_vis, albedo_nir
    USE mo_decomposition,     ONLY: ldc => local_decomposition
    USE mo_control, ONLY: nlev 
    USE mo_time_control, ONLY: lresume, time_set
!---wiso-code
    USE mo_wiso, ONLY: lwiso, nwiso
!---wiso-code-end

    REAL(dp) :: rain(ldc%nproma), snow(ldc%nproma)
!---wiso-code
    REAL(dp) :: wiso_rain(ldc%nproma,nwiso), wiso_snow(ldc%nproma,nwiso)
!---wiso-code-end
    INTEGER :: nproma, jrow, jt, jk, jl

    CALL init_albedo_ice

    box%land_fract(:,:)             = slm(:,:)
    surface%is_land(:,:)            = slm(:,:) > 0.0_dp
    box%lake_fract(:,:)             = alake(:,:)

    surface%is_seaice(:,:)         = .FALSE.
    surface%is_ocean(:,:)          = .FALSE.
    surface%is_ice_for_coup(:,:)   = .FALSE.
    surface%is_ocean_for_coup(:,:) = .FALSE.

    IF (.NOT. lresume) THEN

       box%ocean_fract(:,:)            = 0.0_dp
       box%seaice_fract(:,:)           = 0.0_dp
       box%frac_ice_cover_acc(:,:)     = 0.0_dp

       CALL time_set     ! Computes initial time interpolation weights needed for update_lai in JSBACH inquire mode

       ! Compute surface albedo
       
       box%ocean_fract(:,:)          =                          &
            (1._dp-box%land_fract(:,:)) &
            * (1._dp-box%seaice(:,:)) 
       box%seaice_fract(:,:)         =                          &
            1._dp-box%land_fract(:,:)  &
            - box%ocean_fract(:,:)
       surface%is_seaice(:,:)        =                          &
            box%seaice_fract(:,:) > 0.0_dp &
            .OR. box%lake_fract(:,:) .GE. 0.5_dp
       surface%is_ocean (:,:)        =                          &
            box%ocean_fract (:,:) > 0.0_dp &
            .OR. box%lake_fract(:,:) .GE. 0.5_dp
       surface%is_ocean_for_coup(:,:) =                         &
            slf(:,:) .lt. 1.0_dp                                &
            .and.(.not. surface%is_ocean(:,:))                  
       surface%is_ice_for_coup(:,:)   =                         &
            (surface%is_ocean_for_coup(:,:)  .or.               &
            surface%is_ocean(:,:)         ) .and.               &
            (.not. surface%is_seaice(:,:))

       rain = 0._dp
       snow = 0._dp
       alsol = 0._dp
!---wiso-code
       IF (lwiso) THEN

       wiso_rain(:,:) = 0._dp
       wiso_snow(:,:) = 0._dp
       
       END IF
!---wiso-code-end

       DO jrow=1,ldc%ngpblks
          IF ( jrow == ldc% ngpblks ) THEN
             nproma = ldc% npromz
          ELSE
             nproma = ldc% nproma
          END IF

  IF (lwiso) THEN
          CALL jsbach_inter_1d( nproma,                            &
               COUNT(surface%is_land(1:nproma,jrow)),              &
               precip_rain    = rain(1:nproma),                    &
               precip_snow    = snow(1:nproma),                    &
               albedo         = land%albedo(1:nproma,jrow),        &
               albedo_vis     = land%albedo_vis(1:nproma,jrow),    &
               albedo_nir     = land%albedo_nir(1:nproma,jrow),    &
               mode           = 'inquire',                         &
               mask_land      = surface%is_land(1:nproma,jrow),    &
               kblock         = jrow,                              &
!---wiso-code                                                       
               wiso_precip_rain = wiso_rain(1:nproma,1:nwiso),     &
               wiso_precip_snow = wiso_snow(1:nproma,1:nwiso))
!---wiso-code-end
  ELSE
          CALL jsbach_inter_1d( nproma,                            &
               COUNT(surface%is_land(1:nproma,jrow)),              &
               precip_rain    = rain(1:nproma),                    &
               precip_snow    = snow(1:nproma),                    &
               albedo         = land%albedo(1:nproma,jrow),        &
               albedo_vis     = land%albedo_vis(1:nproma,jrow),    &
               albedo_nir     = land%albedo_nir(1:nproma,jrow),    &
               mode           = 'inquire',                         &
               mask_land      = surface%is_land(1:nproma,jrow),    &
               kblock         = jrow)
  END IF


          CALL update_albedo_ocean((surface%is_ocean(1:nproma,jrow) .or.     &
               surface%is_ocean_for_coup(1:nproma,jrow))        &
               , ocean%albedo(1:nproma,jrow))
          CALL update_albedo_ice((surface%is_seaice(1:nproma,jrow).or.       &
               surface%is_ice_for_coup(1:nproma,jrow))            &
               , ice%surface_temperature(1:nproma,jrow)     &
               , ice%snow_water_equivalent(1:nproma,jrow)   &
               , ice%albedo(1:nproma,jrow)                 )

          ocean%albedo_vis(1:nproma,jrow) = ocean%albedo(1:nproma,jrow)
          ocean%albedo_nir(1:nproma,jrow) = ocean%albedo(1:nproma,jrow)
          ice%albedo_vis(1:nproma,jrow) = ice%albedo(1:nproma,jrow)
          ice%albedo_nir(1:nproma,jrow) = ice%albedo(1:nproma,jrow)

          CALL surface_box_average( nproma             &
               , land%albedo(1:nproma,jrow)            &
               , ice%albedo(1:nproma,jrow)             &
               , ocean%albedo(1:nproma,jrow)           &
               , box%albedo(1:nproma,jrow)             &
               , jrow                                  &
               , surface%is_land(1:nproma,jrow)        &
               , surface%is_ocean(1:nproma,jrow)       &
               , surface%is_seaice(1:nproma,jrow)      &
               , box%land_fract(1:nproma,jrow)         &
               , box%ocean_fract(1:nproma,jrow)        &
               , box%seaice_fract(1:nproma,jrow)      )

          CALL surface_box_average( nproma             &
               , land%albedo_vis(1:nproma,jrow)        &
               , ice%albedo_vis(1:nproma,jrow)         &
               , ocean%albedo_vis(1:nproma,jrow)       &
               , box%albedo_vis(1:nproma,jrow)         &
               , jrow                                  &
               , surface%is_land(1:nproma,jrow)        &
               , surface%is_ocean(1:nproma,jrow)       &
               , surface%is_seaice(1:nproma,jrow)      &
               , box%land_fract(1:nproma,jrow)         &
               , box%ocean_fract(1:nproma,jrow)        &
               , box%seaice_fract(1:nproma,jrow)      )

          CALL surface_box_average( nproma             &
               , land%albedo_nir(1:nproma,jrow)        &
               , ice%albedo_nir(1:nproma,jrow)         &
               , ocean%albedo_nir(1:nproma,jrow)       &
               , box%albedo_nir(1:nproma,jrow)         &
               , jrow                                  &
               , surface%is_land(1:nproma,jrow)        &
               , surface%is_ocean(1:nproma,jrow)       &
               , surface%is_seaice(1:nproma,jrow)      &
               , box%land_fract(1:nproma,jrow)         &
               , box%ocean_fract(1:nproma,jrow)        &
               , box%seaice_fract(1:nproma,jrow)      )

          alsol(1:nproma,jrow)  = land%albedo (1:nproma,jrow)
          alsoi(1:nproma,jrow)  = ice%albedo  (1:nproma,jrow)
          alsow(1:nproma,jrow)  = ocean%albedo(1:nproma,jrow)
          albedo(1:nproma,jrow) = box%albedo(1:nproma,jrow)

       END DO
       land%roughness_momentum(:,:)   = 0.1_dp
       land%roughness_heat(:,:)       = 0.1_dp
       ocean%roughness(:,:)           = 0.001_dp
       ice%roughness(:,:)             = 0.001_dp

       alsol = MERGE(alsol, 0._dp, surface%is_land)

       ice%snow_water_equivalent(:,:) = 0.0_dp
       land%zhsoil(:,:)               = 0.0_dp
       land%zcair(:,:)                = 0.0_dp
       land%zcsat(:,:)                = 0.0_dp
!---wiso-code
       IF (lwiso) THEN

       land%zwisocair(:,:,:)      = 0.0_dp
       land%zwisocsat(:,:,:)      = 0.0_dp
       land%zwisocair_fra(:,:,:)  = 0.0_dp
       land%zwisocsat_fra(:,:,:)  = 0.0_dp
       
       END IF
!---wiso-code-end

    END IF

  END SUBROUTINE init_surface

END MODULE mo_surface
