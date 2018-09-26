module mo_cbal_cpools
    USE mo_kind,  ONLY : dp
    implicit none

    public :: update_Cpools
    public :: relocate_carbon
    public :: relocate_carbon_desert
    public :: relocate_carbon_fire
    public :: relocate_carbon_damage
    public :: printCpoolParameters

    ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ===============================================================================


    private

    REAL(dp), Parameter :: Kg_C_perMol            = 12.011e-3_dp ! self explaining, or ?
    REAL(dp), parameter :: Q10                    = 1.8_dp       ! Empirical parameter for temperature dependence of heterotrophic respiration rate
    REAL(dp), parameter :: kappa                  = 1.0_dp       ! Empirical parameter for soil moisture dependence of heterotrophic respiration rate
    REAL(dp), parameter :: tau_Cpool_reserve      = 1.5_dp *365._dp  ! time constant by which reserve pool is depreciated [days]
    REAL(dp), parameter :: tau_Cpool_wood         = 100._dp*365._dp  ! time constant by which wood pool is depreciated  [days]
    REAL(dp), parameter :: tau_Cpool_fast         = 1.5_dp *365._dp  ! time constant by which fast pool is depreciated  [days]
    REAL(dp), parameter :: tau_Cpool_slow         = 150._dp*365._dp  ! time constant by which slow pool is depreciated  [days]
    REAL(dp), parameter :: referenceTemp_Q10      = 273.15_dp    ! Reference temperature in the Q10 formula [Kelvin]
    REAL(dp), parameter :: greenC2leafC           = 2.0_dp       ! ratio of carbon green pool (leaves, fine roots, starches, sugars) to leaf carbon
    REAL(dp), parameter :: reserveC2leafC         = 1.0_dp       ! ratio of carbon in reserve pool (starches, sugars) to leaf carbon 
    REAL(dp), parameter :: frac_above_groundC     = 0.5_dp       ! fraction of wood allocated above ground
    REAL(dp), parameter :: frac_C_litter2atmos    = 0.1_dp       ! fraction of heterotrophic loss of litter emitted to the atmosphere (rest enters soil pool)
    REAL(dp), parameter :: sec_per_day            = 86400._dp    ! seconds per day
    REAL(dp), parameter :: alpha_critical         = 0.35         ! critical value of relative soil moisture below which heterotrophic respiration almost stops
    REAL(dp), parameter :: alpha_min              = 0.05         ! effective lowest value for soil moisture entering the heterotrophic respiration

contains

  ! --- printCpoolParameters() -----------------------------------------------------------------------------------------------------------
  !
  ! Prints out pthe parameters used for the carbon pool model
  !
  subroutine printCpoolParameters(outUnit)
    use mo_exception,  ONLY: message,real2string
    use mo_doctor, ONLY: nerr
    integer,optional,intent(in) :: outUnit

    integer :: iout

    if(present(outUnit)) then
       iout=outUnit
    else
       iout=nerr
    end if

    call message("","=== CpoolParameters =========================================",kout=iout)
    call message("","              Q10="//trim(real2string(Q10)),kout=iout)
    call message("","            kappa="//trim(real2string(kappa)),kout=iout)
    call message("","referenceTemp_Q10="//trim(real2string(referenceTemp_Q10))//" Kelvin",kout=iout)
    call message("","tau_Cpool_reserve="//trim(real2string(tau_Cpool_reserve/365.))//" years",kout=iout)
    call message("","   tau_Cpool_wood="//trim(real2string(tau_Cpool_wood/365.))//" years",kout=iout)
    call message("","   tau_Cpool_fast="//trim(real2string(tau_Cpool_fast/365.))//" years",kout=iout)
    call message("","   tau_Cpool_slow="//trim(real2string(tau_Cpool_slow/365.))//" years",kout=iout)
    call message("","     greenC2leafC="//trim(real2string(greenC2leafC)),kout=iout)
    call message("","   reserveC2leafC="//trim(real2string(reserveC2leafC)),kout=iout)
    call message(""," frac_above_groundC="//trim(real2string(frac_above_groundC)),kout=iout)
    call message("","frac_C_litter2atmos="//trim(real2string(frac_C_litter2atmos)),kout=iout)
    call message("","=============================================================",kout=iout)

  end subroutine printCpoolParameters


  ! --- update_Cpools() -----------------------------------------------------------------------------------------------------------
  !
  ! Updates carbon pools and computes soil respiration rate and Net Ecosystem Product [mol(C)/m^2 s].
  ! Has to be called once each day for each covertype, where a carbon pool exists.
  !
  elemental pure subroutine update_Cpools(LAI, LAI_previous, NPPrate, topSoilTemp, alpha,                &
                                          frac_npp_2_woodPool,frac_npp_2_reservePool,                    &
                                          tau_Cpool_litter_leaf,tau_Cpool_litter_wood,LAI_shed_constant, &
                                          frac_C_fast2atmos,Max_C_content_woods, specific_leaf_area_C,   &
                                          Cpool_green, Cpool_woods, Cpool_reserve,                       &
                                          Cpool_litter_leaf, Cpool_litter_wood, Cpool_fast, Cpool_slow,  &
                                          soilResp_rate,NPP_flux_correction,total_litter_flux)
    real(dp),intent(in)    :: LAI                    !! Yesterdays mean LAI
    real(dp),intent(in)    :: LAI_previous           !! The day before yesterdays mean LAI
    real(dp),intent(in)    :: NPPrate                !! Yesterdays mean NPPrate [mol(C)/m^2 s]
    real(dp),intent(in)    :: topSoilTemp            !! Yesterdays mean temperature of upper soil layer [degree Kelvin]
    real(dp),intent(in)    :: alpha                  !! Yesterdays mean water stress factor (between 0 and 1)
    real(dp),intent(in)    :: frac_npp_2_woodPool    !! Fraction of NPP to be put maximally into the green pool
    real(dp),intent(in)    :: frac_npp_2_reservePool !! Fraction of NPP to be put into the optimally into the reserve pool
    real(dp),intent(in)    :: tau_Cpool_litter_leaf  !! Time constant by which leaf litter pool is depreciated  [days]
    real(dp),intent(in)    :: tau_Cpool_litter_wood  !! Time constant by which woody litter pool is depreciated  [days]
    real(dp),intent(in)    :: LAI_shed_constant      !! Leaf shedding at a constant rate for evergreens etc.
    real(dp),intent(in)    :: frac_C_fast2atmos      !!Fraction of heterotrophic loss of fast pool emitted to the atmosphere (rest enters slow pool)
    real(dp),intent(in)    :: Max_C_content_woods    !! Maximum carbon content of wood pool (from lctLib)
    real(dp),intent(in)    :: specific_leaf_area_C   !! Specific leaf area (from lctLib) [m^2(leaf)/mol(Carbon)]
    real(dp),intent(inout) :: Cpool_green            !! Value of green carbon pool: on input last value; updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods            !! Value of wood carbon pool: on input last value; updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve          !! Value of reserve carbon pool: on input last value; updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_leaf      !! Value of leaf litter carbon pool: on input last value; updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_wood      !! Value of woody litter carbon pool: on input last value; updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_fast             !! Value of fast soil carbon pool: on input last value; updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_slow             !! Value of slow soil carbon pool: on input last value; updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(out)   :: soilResp_rate          !! Soil (=heterotrophic) respiration rate  [mol(C)/m^2 s] (Note: this is a loss, therefore nrgative)
    real(dp),intent(out)   :: NPP_flux_correction    !! Amount by which the NPP rate entering the routine has to be corrected. This correction arises
                                                     !! .. either because otherwise the reserve pool would get negative (positive correction), or the 
                                                     !! .. wood pool would exceed its maximum value (negative correction). [mol(C)/m^2 s]
    real(dp),intent(out)   :: total_litter_flux      !! Diagnostic of the total flux of litter from leaves and wood and excess carbon 
                                                     !! ... [mol(C)/m^2(canopy) s]

    ! locals

    real(dp) :: d_Cpool_fast,d_Cpool_slow
    real(dp) :: d_Cpool_litter_leaf,d_Cpool_litter_wood
    real(dp) :: Cpool_green_max,Cpool_reserve_optimal
    real(dp) :: NPP_2_greenPool,NPP_2_woodPool,NPP_2_reservePool
    real(dp) :: C_2_litter_leafPool,C_2_litter_woodPool
    real(dp) :: C_2_slowPool,C_2_fastPool
    real(dp) :: excess_carbon,missing_carbon
    real(dp) :: leaf_shedding
    real(dp) :: alpha_mod

    ! Preparations
    !-------------------------------------------

    alpha_mod = max(alpha_min,(alpha-alpha_critical)/(1._dp-alpha_critical)) !! modified soil moist. such that het. resp. is zero below alpha_critical

    ! Initializations
    !-------------------------------------------
    
    NPP_flux_correction = 0.0_dp
    Cpool_green_max = 0.0_dp
    Cpool_reserve_optimal = 0.0_dp
    total_litter_flux = 0.0_dp
    Cpool_green_max = greenC2leafC * LAI / specific_leaf_area_C
    Cpool_reserve_optimal = reserveC2leafC * LAI / specific_leaf_area_C

    ! Preliminary determination of distribution of NPP to the various pools (may be corrected afterwards)
    !----------------------------------------------------------------------------------------------------

    if(NPPrate < 0.0_dp) then                                               !! In case of negative NPP
       NPP_2_reservePool = NPPrate * sec_per_day                            !! .. the reserve pool has to carry the whole burden
       NPP_2_greenPool   = 0.0_dp
       NPP_2_woodPool    = 0.0_dp
    else                                                                    !! But for positive NPP it is distributed according to 
       NPP_2_woodPool    = frac_npp_2_woodPool * NPPrate * sec_per_day      !! .. predefined relative fraction
       NPP_2_reservePool = frac_npp_2_reservePool * NPPrate * sec_per_day
       NPP_2_greenPool   = (1.0_dp - frac_npp_2_woodPool - frac_npp_2_reservePool)*NPPrate * sec_per_day
    end if

    ! Update Wood Pool (structural carbon of living plants)
    !------------------------------------------------------

    if(NPPrate > 0.0_dp) then                                         !! Only if NPP>0 carbon enters the wood pool.
       Cpool_woods = Cpool_woods + NPP_2_woodPool                     !! Then it is attempted to put all NPP share into it
       excess_carbon = max(0.0_dp,Cpool_woods-Max_C_content_woods)    !! .. but thereby the pool may get too large
       Cpool_woods = Cpool_woods - excess_carbon                      !! .. so that the excess carbon has once more to be subtracted
       NPP_2_greenPool = NPP_2_greenPool + excess_carbon              !! .. and instead made available to the green pool
    end if
                                                                      !! Assuming that forests continously die at the inverse lifetime of trees,
    C_2_litter_woodPool = frac_above_groundC * Cpool_woods / tau_Cpool_wood    !! .. above ground wood is put into the woody litter pool
    C_2_slowPool = (1._dp - frac_above_groundC) * Cpool_woods / tau_Cpool_wood !! .... below ground wood is put into the slow soil pool
    Cpool_woods = Cpool_woods - Cpool_woods / tau_Cpool_wood                   !! .... and both are substracted from the wood pool

    ! Update Reserve Pool (sugars and starches) and determination how much will enter the Green Pool
    !-----------------------------------------------------------------------------------------------

                                                                      !! Some carbon of the reserve pool is always lost (fruits, barks,..)
    C_2_litter_leafPool = Cpool_reserve / tau_Cpool_reserve           !! .. and enters the leaf litter pool
    Cpool_reserve = Cpool_reserve - C_2_litter_leafPool               !! .... and is substracted from the reserve pool

    if(NPP_2_reservePool < 0.0_dp) then                               !! If NPP is negative
       Cpool_reserve = Cpool_reserve + NPP_2_reservePool              !! .. take it from the reserve pool. But
       if(Cpool_reserve < 0.0_dp) then                                !! .. if thereby the reserve pool gets negative
          NPP_flux_correction = -Cpool_reserve/sec_per_day            !! .... we give the deficit as a flux correction back to the calling routine 
          Cpool_reserve = 0.0_dp                                      !! .... and force the reserve pool to zero.
       end if                                                         !!
    else                                                              !! Otherwise 
       if(Cpool_reserve < Cpool_reserve_optimal) then                 !! .. If the reserve pool is smaller than optimal,
          Cpool_reserve = Cpool_reserve +  NPP_2_reservePool          !! .... it is filled by the available NPP
          excess_carbon =                                         &   !! .... Thereby it may happen that it gets larger than optimal
                    max(0.0_dp,Cpool_reserve - Cpool_reserve_optimal) !! .... so that there is excess carbon
          Cpool_reserve = Cpool_reserve - excess_carbon               !! .... that needs not be taken up
          NPP_2_greenPool = NPP_2_greenPool + excess_carbon           !! .... but can better be used to increase the green pool.
       else                                                           !! ...Otherwise (reserve pool is already larger than optimal)
          NPP_2_greenPool = NPP_2_greenPool + NPP_2_reservePool       !! .... all NPP is left for the green pool
       end if
    end if


    ! Update Green Pool (leaves and fine roots): 1. Leaf shedding
    !------------------------------------------------------------

    IF (LAI >= LAI_previous) THEN                                    !! IF LAI is increasing or constant the leaf shedding is given by
      leaf_shedding = LAI * LAI_shed_constant                        !! .. a constant loss of leaves for evergreens, raingreens, grasses.
    ELSE                                                             !! .... Otherwise
      leaf_shedding = LAI_previous - LAI                             !! .... the leaf shedding is given by the decrease in LAI.
    END IF

    Cpool_green = Cpool_green - leaf_shedding / specific_leaf_area_C !! The leaf shedding removes carbon from the green pool,
    C_2_litter_leafPool = C_2_litter_leafPool + &                    !! .. which
                          leaf_shedding / specific_leaf_area_C       !! .... is added to the leaf litter pool.
    missing_carbon = -min(0.0_dp,Cpool_green)                        !! Thereby the green pool may go below zero, so that there is carbon missing
    Cpool_green = Cpool_green + missing_carbon                       !! .. which is again added to the green pool and
    C_2_litter_leafPool = C_2_litter_leafPool - missing_carbon       !! .... removed from the leaf litter pool

    excess_carbon = max(0.0_dp,Cpool_green - Cpool_green_max)         !! If by a reduction of LAI the maximum value of the green pool
                                                                      !! .. is still smaller than the current value,
    C_2_fastPool = excess_carbon                                      !! .... there is excess carbon to be put into the fast soil pool
    Cpool_green = Cpool_green - excess_carbon                         !! .... and to be substracted from the green pool

    ! Update Green Pool (leaves and fine roots): 2. Growth. In case of too much NPP, try to put it into the reserve pool.
    !------------------------------------------------------------------------------------------------------------------

    if(NPP_2_greenPool > 0.0_dp) then                                 !! If there is NPP available for the green pool,
       Cpool_green = Cpool_green + NPP_2_greenPool                    !! .. it is filled by the available NPP.
       excess_carbon = max(0.0_dp,Cpool_green - Cpool_green_max)      !! .. Thereby it may get larger than appropriate for current LAI.
       Cpool_green = Cpool_green - excess_carbon                      !! .. This excess carbon needs not be taken up but
       if(Cpool_reserve < Cpool_reserve_optimal) then                 !! .. if the reserve pool is smaller than optimal
          Cpool_reserve = Cpool_reserve + excess_carbon               !! .... it is tried to put the carbon there.
          excess_carbon =                                          &  !! .... Thereby it may happen that
                    max(0.0_dp,Cpool_reserve - Cpool_reserve_optimal) !! .... the reserve pool increases beyond the optimal value. In that case
          Cpool_reserve = Cpool_reserve - excess_carbon               !! .... the excess carbon is once more removed from the reserve pool
       end if                                                         !! .... and instead
       C_2_fastPool = C_2_fastPool + excess_carbon                    !! .... put into the fast pool
       total_litter_flux = excess_carbon / sec_per_day                !! Remember the amount of NPP that could not be stored in the plants for  
    endif                                                             !! .. diagnostic output

    ! Update Leaf Litter Pool (litter from dead leaves, fruits, debris from bark)

    d_Cpool_litter_leaf = -MIN(alpha_mod**kappa * q10**((topSoilTemp - referenceTemp_Q10)/10._dp)  / & !! Decomposition of leaf litter
                          tau_Cpool_litter_leaf,1._dp) * Cpool_litter_leaf                         !! .. according to Q10-model
    Cpool_litter_leaf = Cpool_litter_leaf + C_2_litter_leafPool + d_Cpool_litter_leaf

    ! Update Wood Litter Pool (litter from all above ground dead woody material)

    d_Cpool_litter_wood = -MIN(q10**((topSoilTemp - referenceTemp_Q10)/10._dp)  / &                !! Decomposition of woody litter
                          tau_Cpool_litter_wood,1._dp) * Cpool_litter_wood                         !! .. according to Q10-model
                                                                                                   !! .. without soil moisture dependence
    Cpool_litter_wood = Cpool_litter_wood + C_2_litter_woodPool + d_Cpool_litter_wood

    ! Update Fast Pool (fast respiration pool from already destructured leafs, fine roots, fruits)

    d_Cpool_fast = -MIN(alpha_mod**kappa * q10**((topSoilTemp - referenceTemp_Q10)/10._dp)  /        & !! Fast soil respiration
                   tau_Cpool_fast,1._dp) * Cpool_fast                                              !! .. according to Q10-model
    Cpool_fast = Cpool_fast                                             &
                 - (1._dp - frac_C_litter2atmos) * d_Cpool_litter_leaf  &
                 + C_2_fastPool                                         &
                 + d_Cpool_fast
              
    ! Update Soil Pool (slow Carbon pool like organic C-materials in soil, dead steams and dead roots)
    !-------------------------------------------------------------------------------------------------
    d_Cpool_slow = -MIN(alpha_mod**kappa * q10**((topSoilTemp - referenceTemp_Q10)/10._dp)  / tau_Cpool_slow    & !! Slow soil respiration by funghi
                             ,1._dp)  * Cpool_slow                                                                 !! .. according to Q10-model
    Cpool_slow = Cpool_slow                                  &
                 - (1._dp - frac_C_litter2atmos) * d_Cpool_litter_wood  &
              - (1._dp-frac_C_fast2atmos) * d_Cpool_fast      &
              + C_2_slowPool                                 &
              + d_Cpool_slow

    ! Soil respiration 

    soilResp_rate = (frac_C_litter2atmos * d_Cpool_litter_leaf + & !! heterotrophic respiration from leaf litter decomposition
                     frac_C_litter2atmos * d_Cpool_litter_wood + & !! heterotrophic respiration from woody litter decomposition
                     frac_C_fast2atmos * d_Cpool_fast +          & !! heterotrophic soil respiration from fast pool
                     d_Cpool_slow) / sec_per_day                   !! heterotrophic soil respiration from slow pool

  end subroutine update_Cpools

  ! --- relocate_carbon() ---------------------------------------------------------------------------------------------------
  !
  ! The routine models the consequences of transient shifts in the cover fractions of the tiles for the carbon pools.
  ! This change in vegetation composition is either due to human land-cover change or due to vegetation dynamics.
  ! Concerning the carbon stored on those partitions of shrinking tiles that are lost to extending tiles during the time interval considered:
  ! In the case of landcover change the carbon from living plant pools (Cpool_green, Cpool_reserve and Cpool_woods)
  ! is relocated, partly by a release to the atmosphere (e.g. by fire stubbing) and the remaining carbon into soil pools.
  ! In the case of vegetation dynamics it is assumed that the carbon from living plant pools is already lost to other reservoirs
  ! (atmosphere, litter, soil) by the processes that removed the vegetation (general carbon flow reflecting minimum mortality, fire, wind break)

  subroutine relocate_carbon(cf_old, cf_new, veg_fract_correction, fract_small,               &
                             Cpool_green, Cpool_woods, Cpool_reserve,                         &
                             Cpool_litter_leaf, Cpool_litter_wood, Cpool_fast, Cpool_slow,    &
                             frac_wood_2_atmos, frac_green_2_atmos, frac_reserve_2_atmos,     &
                             carbon_2_atmos, carbon_2_fastSoilPool, carbon_2_slowSoilPool)

    use mo_exception, ONLY : finish

    real(dp),intent(in)    :: cf_old(:,:)                !! Cover fractions before landcover change or vegetation dynamics [] 
    real(dp),intent(in)    :: cf_new(:,:)                !! Cover fraction after landcover change or vegetation dynamics []
    real(dp),intent(in)    :: veg_fract_correction(:,:)  !! Correction factor for cover fractions 1-exp(-LAI_max/2) (accounts vor sparseness of vegetation)
    real(dp),intent(in)    :: fract_small                !! very small fraction (minimum of cf_new on non-glacier land points)
    real(dp),intent(inout) :: Cpool_green(:,:)           !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods(:,:)           !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve(:,:)         !! Value of reserve carbon pool [mol(C)/m^(canopy)2]
    real(dp),intent(inout) :: Cpool_litter_leaf(:,:)     !! Value of leaf litter carbon pool [mol(C)/m^(canopy)2]
    real(dp),intent(inout) :: Cpool_litter_wood(:,:)     !! Value of wood litter carbon pool [mol(C)/m^(canopy)2]
    real(dp),intent(inout) :: Cpool_fast(:,:)            !! Value of fast soil carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_slow(:,:)            !! Value of slow soil carbon pool [mol(C)/m^2(canopy)]
                                                         !! .. to slow soil pool [mol(C)/m^2(vegetated area)]
    real(dp),optional,intent(in)  :: frac_wood_2_atmos         !! Fraction of wood pool immediately emitted to atmosphere (e.g. by burning vegetation down)
    real(dp),optional,intent(in)  :: frac_green_2_atmos        !! Fraction of green pool immediately emitted to atmosphere (e.g. by burning vegetation down)
    real(dp),optional,intent(in)  :: frac_reserve_2_atmos      !! Fraction of reserve pool immediately emitted to atmosphere (e.g. by burning vegetation down)
    real(dp),optional,intent(out) :: carbon_2_atmos(:)         !! Amount of carbon directly emitted to atmosphere
                                                               !! ..in this timestep [mol(C)/m^2(vegetated area)]
    real(dp),optional,intent(out) :: carbon_2_fastSoilPool(:)  !! Amount of carbon relocated by land cover change from green and reserve pool
                                                               !! .. to fast soil pool [mol(C)/m^2(vegetated area)]
    real(dp),optional,intent(out) :: carbon_2_slowSoilPool(:)  !! Amount of carbon relocated by land cover change from wood pool
                                                               !! .. to slow soil pool [mol(C)/m^2(vegetated area)]

    !! locals

    real(dp) :: Cfast_2_fast(1:size(cf_old,DIM=1))                      !! Total carbon relocated from fast pools of shrinking to fast pools of extending tiles
    real(dp) :: Cslow_2_slow(1:size(cf_old,DIM=1))                      !! Total carbon relocated from slow pools of shrinking to slow pools of extending tiles
    real(dp) :: Clitter_leaf_2_litter_leaf(1:size(cf_old,DIM=1))        !! Total carbon relocated from leaf litter pools of shrinking to leaf litter pools of extending tiles
    real(dp) :: Clitter_wood_2_litter_wood(1:size(cf_old,DIM=1))        !! Total carbon relocated from wood litter pools of shrinking to wood litter pools of extending tiles
    logical  :: tiles_extend(1:size(cf_old,DIM=1),1:size(cf_old,DIM=2)) !! Logical mask indicating those tiles which are extending by landcover change
    real(dp) :: shrinkRatio (1:size(cf_old,DIM=1),1:size(cf_old,DIM=2)) !! Ratio old and new cover fractions of extending tiles
    real(dp) :: sum_cf_new_old(1:size(cf_old,DIM=1))                    !! Sum of cover fractions added to old cover fractions of extending tiles
    logical  :: land_cover_change                                       !! Flag: land cover change (.true.) or vegetation dynamics (.false.) 
    integer :: ntiles,i

    !! Avoid wrong use of the routine!
    !! This routine handles either land cover change or vegetation dynamics:
    !! For land cover change pass frac_wood_2_atmos, frac_green_2_atmos, frac_reserve_2_atmos, carbon_2_atmos,
    !! carbon_2_fastSoilPool and carbon_2_slowSoilPool to this routine.

    if (present(frac_wood_2_atmos) .or. present(frac_green_2_atmos) .or. present(frac_reserve_2_atmos) .or. &
           present(carbon_2_atmos) .or. present(carbon_2_fastSoilPool) .or. present(carbon_2_slowSoilPool)) then
       if (.not. present(frac_wood_2_atmos) .or. .not. present(frac_green_2_atmos) .or. .not. present(frac_reserve_2_atmos) .or. &
           .not. present(carbon_2_atmos) .or. .not. present(carbon_2_fastSoilPool) .or. .not. present(carbon_2_slowSoilPool)) &
          call finish('relocate_carbon','at least one variable missing to handle land cover change')
       land_cover_change=.true.
    else
       land_cover_change=.false.
    endif

    !! preparations

    ntiles = size(cf_old,DIM=2)

    !! determine extending and shrinking tiles

    where(cf_new(:,:) > cf_old(:,:))  
       tiles_extend(:,:) = .true.
    elsewhere 
       tiles_extend(:,:) = .false.
    end where

    !! determine amount of carbon released from living plant pools to atmosphere, litter and soil

    if (land_cover_change) then

       carbon_2_slowSoilPool(:) = SUM( (cf_old(:,:) - cf_new(:,:)) * veg_fract_correction(:,:) &  !! C to be relocated from wood pool of 
                                      * (1.0_dp - frac_wood_2_atmos) * Cpool_woods(:,:)        &  !! .. shrinking tiles to slow pool of extending tiles
                                      ,MASK=.not. tiles_extend(:,:),DIM=2 )

    carbon_2_atmos(:) = SUM( (cf_old(:,:) - cf_new(:,:)) * veg_fract_correction(:,:) &
                               * ( frac_green_2_atmos   * Cpool_green(:,:)           &
                               + frac_wood_2_atmos    * Cpool_woods(:,:)             &
                               + frac_reserve_2_atmos * Cpool_reserve(:,:))          &
                               ,MASK=.not. tiles_extend(:,:),DIM=2)

       carbon_2_fastSoilPool(:) = SUM( (cf_old(:,:) - cf_new(:,:)) * veg_fract_correction(:,:) &  !! C to be relocated from green, and reserve pool of 
                                      *( (1.0_dp - frac_green_2_atmos) * Cpool_green(:,:)      &  !! .. shrinking tiles to fast pool of extending tiles
                                      + (1.0_dp - frac_reserve_2_atmos) * Cpool_reserve(:,:))  &
                                      ,MASK=.not. tiles_extend(:,:),DIM=2 )

    else  !! rescale living plant pools of shrinking tiles 

       where(.NOT. tiles_extend(:,:) .AND. cf_new(:,:) > 0.5_dp*fract_small) !! cf_new should never be substantially smaller than
                                                                             !! small_fract on non-glacier land points
          shrinkRatio(:,:) = cf_old(:,:) / cf_new(:,:)
       elsewhere
          shrinkRatio(:,:) = 1._dp
       end where      

       Cpool_green(:,:)   = Cpool_green(:,:)   * shrinkRatio(:,:) 
       Cpool_reserve(:,:) = Cpool_reserve(:,:) * shrinkRatio(:,:)
       Cpool_woods(:,:)   = Cpool_woods(:,:)   * shrinkRatio(:,:)       

    end if

    Cfast_2_fast(:)   =  SUM( (cf_old(:,:) - cf_new(:,:)) * veg_fract_correction(:,:)  &             !! C to be relocated from fast pool of 
                            * Cpool_fast(:,:),MASK=.not. tiles_extend(:,:),DIM=2 )                   !! .. shrinking to extending tiles

    Cslow_2_slow(:)   =  SUM( (cf_old(:,:) - cf_new(:,:)) * veg_fract_correction(:,:)  &             !! C to be relocated from slow pool of
                            * Cpool_slow(:,:),MASK=.not. tiles_extend(:,:),DIM=2 )                   !! .. shrinking to extending tiles

    Clitter_leaf_2_litter_leaf(:) = SUM( (cf_old(:,:) - cf_new(:,:)) * veg_fract_correction(:,:)  &  !! C to be relocated from leaf litter pool of
                                       * Cpool_litter_leaf(:,:),MASK=.not. tiles_extend(:,:),DIM=2 ) !! .. shrinking to extending tiles 

    Clitter_wood_2_litter_wood(:) = SUM( (cf_old(:,:) - cf_new(:,:)) * veg_fract_correction(:,:)  &  !! C to be relocated from wood litter pool of
                                       * Cpool_litter_wood(:,:),MASK=.not. tiles_extend(:,:),DIM=2 ) !! .. shrinking to extending tiles 

    !! Lower down the carbon density [mol(C)/m^2(canopy)] for extending tiles

    where(tiles_extend(:,:) .AND. cf_new(:,:) > 0.5_dp*fract_small)
       shrinkRatio(:,:) = cf_old(:,:) / cf_new(:,:)
    elsewhere
       shrinkRatio(:,:) = 1._dp
    end where

    Cpool_green(:,:)       = Cpool_green(:,:)       * shrinkRatio(:,:) 
    Cpool_reserve(:,:)     = Cpool_reserve(:,:)     * shrinkRatio(:,:)
    Cpool_woods(:,:)       = Cpool_woods(:,:)       * shrinkRatio(:,:)
    Cpool_litter_leaf(:,:) = Cpool_litter_leaf(:,:) * shrinkRatio(:,:)
    Cpool_litter_wood(:,:) = Cpool_litter_wood(:,:) * shrinkRatio(:,:)
    Cpool_fast(:,:)        = Cpool_fast(:,:)        * shrinkRatio(:,:)
    Cpool_slow(:,:)        = Cpool_slow(:,:)        * shrinkRatio(:,:)

    !! sum of cover fractions added to cover fractions of extending tiles

    sum_cf_new_old(:) = 0._dp

    do i = 1,ntiles
       where(tiles_extend(:,i) .AND. (cf_new(:,i) - cf_old(:,i)) > EPSILON(1._dp))
          sum_cf_new_old(:) = sum_cf_new_old(:) + (cf_new(:,i) - cf_old(:,i)) * veg_fract_correction(:,i)
       end where
    end do

     !! Distribute freed carbon to extending tiles

     do i = 1,ntiles
        if (land_cover_change) then
           where(tiles_extend(:,i) .AND. (cf_new(:,i) - cf_old(:,i)) > EPSILON(1._dp))
              Cpool_litter_leaf(:,i) = Cpool_litter_leaf(:,i) + Clitter_leaf_2_litter_leaf(:) *   &
                                       (cf_new(:,i) - cf_old(:,i)) / (sum_cf_new_old(:) * cf_new(:,i))
              Cpool_litter_wood(:,i) = Cpool_litter_wood(:,i) + Clitter_wood_2_litter_wood(:) *   &
                                       (cf_new(:,i) - cf_old(:,i)) / (sum_cf_new_old(:) * cf_new(:,i))
              Cpool_fast(:,i) = Cpool_fast(:,i) + (carbon_2_fastSoilPool(:) + Cfast_2_fast(:)) *  &
                                (cf_new(:,i) - cf_old(:,i)) / (sum_cf_new_old(:) * cf_new(:,i))
              Cpool_slow(:,i) = Cpool_slow(:,i) + (carbon_2_slowSoilPool(:) + Cslow_2_slow(:)) *  &
                                (cf_new(:,i) - cf_old(:,i)) / (sum_cf_new_old(:) * cf_new(:,i))
           end where
        else
          where(tiles_extend(:,i) .AND. (cf_new(:,i) - cf_old(:,i)) > EPSILON(1._dp))
             Cpool_litter_leaf(:,i) = Cpool_litter_leaf(:,i) + Clitter_leaf_2_litter_leaf(:) * &
                                      (cf_new(:,i) - cf_old(:,i)) / (sum_cf_new_old(:) * cf_new(:,i))
             Cpool_litter_wood(:,i) = Cpool_litter_wood(:,i) + Clitter_wood_2_litter_wood(:) * &
                                      (cf_new(:,i) - cf_old(:,i)) / (sum_cf_new_old(:) * cf_new(:,i))
             Cpool_fast(:,i) = Cpool_fast(:,i) + Cfast_2_fast(:) *  &
                                      (cf_new(:,i) - cf_old(:,i)) / (sum_cf_new_old(:) * cf_new(:,i))
             Cpool_slow(:,i) = Cpool_slow(:,i) + Cslow_2_slow(:) *  &
                                      (cf_new(:,i) - cf_old(:,i)) / (sum_cf_new_old(:) * cf_new(:,i))
          end where
       end if
    end do

  end subroutine relocate_carbon


  ! --- relocate_carbon_desert() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of desert dynamics for the carbon pools. More precisely:
  ! It is assumed that no carbon fluxes have to be considered as these are already taken into account by the other routines in this module.
  ! So, the amount of carbon per unit area is only rescaled to conserv the carbon mass.
  ! The only exception is the green pool, which should never exceed its prescribed maximum carbon content.
  ! In this (very exceptional) case the excessive carbon is transferred to the leaf litter pool.

  pure subroutine relocate_carbon_desert(nidx, ntiles, veg_ratio_max, veg_ratio_max_old, fract_small,  &
                                         lai, sla, Cpool_green, Cpool_woods, Cpool_reserve,            &
                                         Cpool_litter_leaf, Cpool_litter_wood, Cpool_fast, Cpool_slow)

    integer,intent(in)     :: nidx                        !! Vector length
    integer,intent(in)     :: ntiles                      !! Number of tiles
    real(dp),intent(in)    :: veg_ratio_max(:)            !! Fractional cover of vegetated area
    real(dp),intent(in)    :: veg_ratio_max_old(:)        !! Fractional cover of vegetated area of last year
    real(dp),intent(in)    :: fract_small                 !! Minimum value of veg_ratio_max for non-glacier land points 
    real(dp),intent(in)    :: lai(:,:)                    !! leaf area index [m^2(leaf)/m^2(canopy)]
    real(dp),intent(in)    :: sla(:,:)                    !! specific leaf area [m^2(leaf)/mol(C)]
    real(dp),intent(inout) :: Cpool_green(:,:)            !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods(:,:)            !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve(:,:)          !! Value of reserve carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_leaf(:,:)      !! Value of leaf litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_wood(:,:)      !! Value of wood litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_fast(:,:)             !! Value of fast soil carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_slow(:,:)             !! Value of slow soil carbon pool [mol(C)/m^2(canopy)]

    !! locals

    real(dp)  ::  veg_ratio_old_new(nidx)          ! ratio of old fraction of vegetated land to new fraction of vegetated land
    real(dp)  ::  cpool_green_excess(nidx,ntiles)  ! 

    WHERE (veg_ratio_max(:) > 0.5_dp*fract_small)
       veg_ratio_old_new(:) = veg_ratio_max_old(:) / veg_ratio_max(:)
    ELSEWHERE
       veg_ratio_old_new(:) = 1._dp                ! never happens on non-glacier land points
    ENDWHERE

    cpool_green_excess(:,:) = MAX(0._dp,cpool_green(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2) - &
                              greenC2leafC * lai(:,:) / sla(:,:))
    cpool_green(:,:) = MIN(cpool_green(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2), greenC2leafC * lai(:,:) / sla(:,:))
    cpool_reserve(:,:) = cpool_reserve(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
    cpool_woods(:,:) = cpool_woods(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
    cpool_litter_leaf(:,:) = cpool_litter_leaf(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2) + cpool_green_excess(:,:)
    cpool_litter_wood(:,:) = cpool_litter_wood(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
    cpool_fast(:,:) = cpool_fast(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
    cpool_slow(:,:) = cpool_slow(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)

  end subroutine relocate_carbon_desert


  ! --- relocate_carbon_fire() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of fire for the carbon pools. More precisely:
  ! It is assumed that for the burned area the carbon from the litter pools (Cpool_litter_leaf, Cpool_litter_wood) is partly 
  ! released to the atmosphere and partly relocated into the woody litter pool
  ! and the carbon from the living plant pools (Cpool_green, Cpool_reserve and Cpool_woods) is partly released to the atmosphere and partly 
  ! relocated into the slow soil pool.
  !
  pure subroutine relocate_carbon_fire(nidx, ntiles, cf_act, cf_burned, cf, fract_small,             &
                                       veg_fract_correction, frac_wood_2_atmos, frac_wood_2_litter,  &
                                       Cpool_green, Cpool_woods, Cpool_reserve,                      &
                                       Cpool_litter_leaf, Cpool_litter_wood, Cpool_slow,             &
                                       carbon_2_WoodLitterPool,                                      &
                                       carbon_2_slowSoilPool, carbon_2_atmos)

    integer,intent(in)     :: nidx                        !! Vector length
    integer,intent(in)     :: ntiles                      !! Number of tiles
    real(dp),intent(in)    :: cf_act(:,:)                 !! Fractional projective cover of each PFT of the vegetated area (before fire)
    real(dp),intent(in)    :: cf_burned(:,:)              !! Fraction of the vegetated area of each tile burned till the last call of this routine
    real(dp),intent(in)    :: cf(:,:)                     !! Cover fractions
    real(dp),intent(in)    :: fract_small                 !! Minimum value for cf_act for non-glacier land points
    real(dp),intent(in)    :: veg_fract_correction(:,:)   !! Correction factor for cover fractions 1-exp(-LAI_max/2) (accounts vor sparseness of vegetation)
    real(dp),intent(in)    :: frac_wood_2_atmos           !! Fraction of wood pool immediately emitted to atmosphere
    real(dp),intent(in)    :: frac_wood_2_litter          !! Fraction of wood pool relocated to the woody litter pool
    real(dp),intent(inout) :: Cpool_green(:,:)            !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods(:,:)            !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve(:,:)          !! Value of reserve carbon pool [mol(C)/m^(canopy)2]
    real(dp),intent(inout) :: Cpool_litter_leaf(:,:)      !! Value of leaf litter carbon pool [mol(C)/m^(canopy)2]
    real(dp),intent(inout) :: Cpool_litter_wood(:,:)      !! Value of wood litter carbon pool [mol(C)/m^(canopy)2]
    real(dp),intent(inout) :: Cpool_slow(:,:)             !! Value of slow soil carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: carbon_2_WoodLitterPool(:)  !! Amount of carbon relocated by wind damage and fire
                                                          !! .. to the wood litter pool [mol(C)/m^2(vegetated area)]
    real(dp),intent(out)   :: carbon_2_slowSoilPool(:)    !! Amount of carbon relocated from wood pool by fire 
                                                          !! .. to slow soil pool [mol(C)/m^2(vegetated area)]
    real(dp),intent(out)   :: carbon_2_atmos(:)           !! Amount of carbon immediately released by fire
                                                          !! .. to the atmosphere [mol(C)/m^2(vegetated area)]
    !! locals

    real(dp)  ::  carbon_2_WoodLitterPool_tiled(nidx,ntiles)
    real(dp)  ::  carbon_2_slowSoilPool_tiled(nidx,ntiles)
    real(dp)  ::  carbon_2_atmos_tiled(nidx,ntiles)

    !! preparations

    carbon_2_WoodLitterPool_tiled(:,:) = 0._dp
    carbon_2_slowSoilPool_tiled(:,:) = 0._dp
    carbon_2_atmos_tiled(:,:) = 0._dp

    !! transfer some carbon from the wood pool to the slow soil pool and
    !! determine amount of carbon released from the wood pool to the woody litter pool as well as the slow soil pool and
    !! from living plant pools and litter pools to the atmosphere

    WHERE (cf_act(:,:)  > 0.5_dp*fract_small)  !! cf_act should never be substantially smaller than fract_small for 
                                               !! non-glacial points 
       Cpool_slow(:,:) = Cpool_slow(:,:) + Cpool_woods(:,:) * (1._dp - frac_wood_2_atmos - frac_wood_2_litter) * &
                         (cf_burned(:,:) / cf_act(:,:))
       carbon_2_WoodLitterPool_tiled(:,:) = Cpool_woods(:,:) * frac_wood_2_litter * &
                                            cf(:,:) * veg_fract_correction(:,:) * (cf_burned(:,:) / cf_act(:,:))
       carbon_2_slowSoilPool_tiled(:,:) = Cpool_woods(:,:) * (1._dp - frac_wood_2_atmos - frac_wood_2_litter) * &
                                          cf(:,:) * veg_fract_correction(:,:) * (cf_burned(:,:) / cf_act(:,:))
       carbon_2_atmos_tiled(:,:) = (Cpool_green(:,:) + Cpool_reserve(:,:) + Cpool_litter_leaf(:,:) + Cpool_litter_wood(:,:) + &
                                   (Cpool_woods(:,:) * frac_wood_2_atmos)) * &
                                   cf(:,:) * veg_fract_correction(:,:) * (cf_burned(:,:) / cf_act(:,:))
    ENDWHERE

    !! sum carbon fluxes for all tiles

    carbon_2_WoodLitterPool(:) = carbon_2_WoodLitterPool(:) + SUM(carbon_2_WoodLitterPool_tiled(:,:),DIM=2)
    carbon_2_slowSoilPool(:) = SUM(carbon_2_slowSoilPool_tiled(:,:),DIM=2)
    carbon_2_atmos(:) = SUM(carbon_2_atmos_tiled(:,:),DIM=2)

    !! lower down the carbon density of litter pools
    !! transfer carbon from the wood pool to the woody litter pool

    WHERE (cf_act(:,:)  >  0.5_dp*fract_small)
       Cpool_litter_leaf(:,:) = Cpool_litter_leaf(:,:) * (1._dp - (cf_burned(:,:)/ cf_act(:,:)))
       Cpool_litter_wood(:,:) = Cpool_litter_wood(:,:) * (1._dp - (cf_burned(:,:)/ cf_act(:,:)))
       Cpool_litter_wood(:,:) = Cpool_litter_wood(:,:) + Cpool_woods(:,:) * frac_wood_2_litter * (cf_burned(:,:) / cf_act(:,:))
    ENDWHERE

    !! lower down the carbon density of living plant pools

    WHERE (cf_act(:,:)  >  0.5_dp*fract_small)
       Cpool_green(:,:) = Cpool_green(:,:) * (1._dp - (cf_burned(:,:) / cf_act(:,:)))
       Cpool_woods(:,:) = Cpool_woods(:,:) * (1._dp - (cf_burned(:,:)/ cf_act(:,:)))
       Cpool_reserve(:,:) = Cpool_reserve(:,:) * (1._dp - (cf_burned(:,:)/ cf_act(:,:)))
     ENDWHERE

  end subroutine relocate_carbon_fire

  ! --- relocate_carbon_damage() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of damages to the vegetation (e.g. wind break) for the carbon pools. More precisely:
  ! It is assumed that for the damaged area the carbon from the living plant pools (Cpool_green, Cpool_reserve and Cpool_woods)
  ! is partly relocated into the litter pools and partly into the slow soil pool.
  !
  pure subroutine relocate_carbon_damage(nidx, ntiles, cf_act, cf_damaged, cf, fract_small, &
                                         veg_fract_correction, frac_wood_2_litter,          &
                                         Cpool_green, Cpool_woods, Cpool_reserve,           &
                                         Cpool_litter_leaf, Cpool_litter_wood, Cpool_slow,  &
                                         carbon_2_slowSoilPool, carbon_2_LeafLitterPool,    &
                                         carbon_2_WoodLitterPool)

    integer,intent(in)     :: nidx                        !! Vector length
    integer,intent(in)     :: ntiles                      !! Number of tiles
    real(dp),intent(in)    :: cf_act(:,:)                 !! Fractional projective cover of each PFT of the vegetated area (before damage)
    real(dp),intent(in)    :: cf_damaged(:,:)             !! Fraction of the vegetated area of each tile damaged till the last call of this routine
    real(dp),intent(in)    :: cf(:,:)                     !! Cover fractions
    real(dp),intent(in)    :: fract_small                 !! Small Fraction: Minimum value of cf_act for non-glacier land points
    real(dp),intent(in)    :: veg_fract_correction(:,:)   !! Correction factor for cover fractions 1-exp(-LAI_max/2) (accounts vor sparseness of vegetation)
    real(dp),intent(in)    :: frac_wood_2_litter          !! Fraction of wood carbon relocated from the wood pool to the the litter wood pool  
    real(dp),intent(inout) :: Cpool_green(:,:)            !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods(:,:)            !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve(:,:)          !! Value of reserve carbon pool [mol(C)/m^(canopy)2]
    real(dp),intent(inout) :: Cpool_litter_leaf(:,:)      !! Value of leaf litter carbon pool [mol(C)/m^(canopy)2]
    real(dp),intent(inout) :: Cpool_litter_wood(:,:)      !! Value of wood litter carbon pool [mol(C)/m^(canopy)2]
    real(dp),intent(inout) :: Cpool_slow(:,:)             !! Value of slow soil carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(out)   :: carbon_2_slowSoilPool(:)    !! Amount of carbon relocated by wind damage
                                                          !! .. to the slow soil pool [mol(C)/m^2(vegetated area)]
    real(dp),intent(out)   :: carbon_2_LeafLitterPool(:)  !! Amount of carbon relocated by wind damage
                                                          !! .. to the leaf litter pool [mol(C)/m^2(vegetated area)]
    real(dp),intent(inout) :: carbon_2_WoodLitterPool(:)  !! Amount of carbon relocated by wind damage and fire
                                                          !! .. to the wood litter pool [mol(C)/m^2(vegetated area)]
    !! locals

    real(dp)  ::  carbon_2_slowSoilPool_tiled(nidx,ntiles)
    real(dp)  ::  carbon_2_WoodLitterPool_tiled(nidx,ntiles)
    real(dp)  ::  carbon_2_LeafLitterPool_tiled(nidx,ntiles)

    !! preparations

    carbon_2_slowSoilPool_tiled(:,:) = 0._dp
    carbon_2_WoodLitterPool_tiled(:,:) = 0._dp
    carbon_2_LeafLitterPool_tiled(:,:) = 0._dp

    !! transfer carbon from the wood pool to the slow soil pool and and to the wood litter pool
    !! transfer carbon from the green and reserve pool to the leaf litter pool
    !! determine amount of carbon relocated

    WHERE (cf_act(:,:)  > 0.5_dp*fract_small) !! cf_act should never be substantially smaller than fract_small for 
                                              !! non-glacial points 
       Cpool_slow(:,:) = Cpool_slow(:,:) + Cpool_woods(:,:) * (1._dp - frac_wood_2_litter) * (cf_damaged(:,:) / cf_act(:,:))
       Cpool_litter_wood(:,:) = Cpool_litter_wood(:,:) + Cpool_woods(:,:) * frac_wood_2_litter * (cf_damaged(:,:) / cf_act(:,:))
       Cpool_litter_leaf(:,:) = Cpool_litter_leaf(:,:) + (Cpool_green(:,:) + Cpool_reserve(:,:)) * (cf_damaged(:,:) / cf_act(:,:))
       carbon_2_slowSoilPool_tiled(:,:) = Cpool_woods(:,:) * (1._dp - frac_wood_2_litter) * &
                                          cf(:,:) * veg_fract_correction(:,:) * (cf_damaged(:,:) / cf_act(:,:))
       carbon_2_WoodLitterPool_tiled(:,:) = Cpool_woods(:,:) * frac_wood_2_litter * &
                                            cf(:,:) * veg_fract_correction(:,:) * (cf_damaged(:,:) / cf_act(:,:))
       carbon_2_LeafLitterPool_tiled(:,:) = (Cpool_green(:,:) + Cpool_reserve(:,:)) * &
                                            cf(:,:) * veg_fract_correction(:,:) * (cf_damaged(:,:) / cf_act(:,:))
    ENDWHERE

    !! sum carbon fluxes for all tiles

    carbon_2_slowSoilPool(:) = SUM(carbon_2_slowSoilPool_tiled(:,:),DIM=2)
    carbon_2_WoodLitterPool(:) = carbon_2_WoodLitterPool(:) + SUM(carbon_2_WoodLitterPool_tiled(:,:),DIM=2)
    carbon_2_LeafLitterPool(:) = SUM(carbon_2_LeafLitterPool_tiled(:,:),DIM=2)

    !! lower down the carbon density of living plant pools

    WHERE (cf_act(:,:)  > 0.5_dp*fract_small)
       Cpool_green(:,:) = Cpool_green(:,:) * (1._dp - (cf_damaged(:,:) / cf_act(:,:)))
       Cpool_woods(:,:) = Cpool_woods(:,:) * (1._dp - (cf_damaged(:,:)/ cf_act(:,:)))
       Cpool_reserve(:,:) = Cpool_reserve(:,:) * (1._dp - (cf_damaged(:,:)/ cf_act(:,:)))
    ENDWHERE

  end subroutine relocate_carbon_damage

end module mo_cbal_cpools
