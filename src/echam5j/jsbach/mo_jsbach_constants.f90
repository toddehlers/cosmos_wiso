!+ Constants for JSBACH
! 
MODULE mo_jsbach_constants

  ! 
  !! Description: 
  !!   Declaration of constants used in JSBACH
  ! 
  !! Current Code Owner: jsbach_admin
  ! 
  !! History: 
  !  
  !! Version   Date        Comment 
  !! -------   ----        ------- 
  !! 0.1       2001/06/28  Original code. Reiner Schnur
  ! 
  ! Code Description: 
  !   Language:           Fortran 90. 
  !   Software Standards: "European Standards for Writing and  
  !     Documenting Exchangeable Fortran 90 Code". 
  ! 
  ! Modules used: 
  ! 
  USE mo_kind, ONLY : dp

  IMPLICIT NONE 

  !! Global (i.e. public) Declarations: 

  !! Global Parameters: 

  REAL(dp),    PARAMETER :: GramsPerKilogram = 1000._dp
  REAL(dp),    PARAMETER :: Kelvin = 273.15_dp                !! Conversion deg C to K
  REAL(dp),    PARAMETER :: Tmelt = 273.15_dp                 !! Melting temperature of snow/ice [K]
  ! Standard atmosphere: (Ref.: Fundamentals of Atmospheric Modeling, p26)
  REAL(dp),    PARAMETER :: LapseRate   = 6.5_dp              !! Temperature lapse rate of US Std. Atmos. [K/km]
  REAL(dp),    PARAMETER :: SeaLevelPressure = 1.01325E5_dp   !! Air pressure at sealevel of standard ICAO-atmosphere [N/m^2 = Pa
  REAL(dp),    PARAMETER :: SeaLevelTemperature = 288._dp     !! Temperature at sea level of standard atmosphere
  ! End standard atmosphere
  REAL(dp),    PARAMETER :: Gravity = 9.80665_dp              !! Effective gravity at surface of earth [m/s^2]

  REAL(dp),    PARAMETER :: UniversalGasConstant = 8.3145_dp  !! Universal gas constant [J / (mole * K)
  REAL(dp),    PARAMETER :: MolecularWeightDryAir = 28.966_dp !! Molecular weight of dry air [g / mole]
  REAL(dp),    PARAMETER :: MolecularWeightWater = 18.02_dp   !! Molecular weight of water [g / mole]
  REAL(dp),    PARAMETER :: GasConstantDryAir = UniversalGasConstant/MolecularWeightDryAir*GramsPerKilogram
                                                              !! Gas constant for dry air [J/(K*kg)], = 287.04
  REAL(dp),    PARAMETER :: GasConstantWaterVapor = UniversalGasConstant/MolecularWeightWater*GramsPerKilogram
                                                              !! Gas constant for water vapor [J/(K*kg)], = 461.4
  REAL(dp),    PARAMETER :: Eps        = GasConstantDryAir / GasConstantWaterVapor
  REAL(dp),    PARAMETER :: SpecificHeatDryAirConstPressure = 1005.46_dp
                                                              !! Specific heat of dry air at constant pressure
  REAL(dp),    PARAMETER :: SpecificHeatVaporConstPressure = 1869.46_dp
                                                              !! Specific heat of water vapor at constant pressure
  !REAL,    PARAMETER :: SpecificHeatMoistAir = 

  REAL(dp),    PARAMETER :: MaxSnowTemp =  0.5_dp             !! Maximum temperature at which snow can fall (C)  
  REAL(dp),    PARAMETER :: MinRainTemp = -0.5_dp             !! Minimum temperature at which rain can fall (C) 
  REAL(dp),    PARAMETER :: MinSnow     = 1.e-5_dp            !! Minimum amount of precip that is allowed to fall as snow (mm)
  REAL(dp),    PARAMETER :: RhoH2O      = 1000._dp            !! Density of liquid water (kg/m^3)
  REAL(dp),    PARAMETER :: EnvLapse    = -0.0006_dp          !! Environmental lapse rate for Penman evaporation [C/m]
  REAL(dp),    PARAMETER :: LatentHeatVaporization = 2.5008e6_dp  !! Latent heat of vaporization [J/kg]
  REAL(dp),    PARAMETER :: LatentHeatSublimation  = 2.8345e6_dp  !! Latent heat of sublimation [J/kg]
  REAL(dp),    PARAMETER :: LatentHeatFusion       = LatentHeatSublimation - LatentHeatVaporization  !! Latent heat of fusion [J/kg]

  REAL(dp),    PARAMETER :: StormThreshold = 0.001_dp         !! Precipitation threshold above which a new storm is declared

  REAL(dp),    PARAMETER :: StefanBoltzmann = 5.67e-8_dp      !! Stefan-Boltzmann constant [W/(m^2 K^4)] 
                                                              !! Fundamentals of Atmospheric Modeling gives 5.67051e-8
  REAL(dp),    PARAMETER :: vonKarman = 0.4_dp                !! von Karman constant for evapotranspiration

  REAL(dp),    PARAMETER :: Emissivity = 0.996_dp             !! Surface emissivity

  REAL(dp),    PARAMETER :: AlbedoDenseVegetation = 0.15_dp   !! Albedo of dense vegetation (Brutsaert, 1982, p136)
  REAL(dp),    PARAMETER :: AbsorbtivitySoilBelowCanopy = 0.05_dp !! Absorbtivity of soil below dense vegetation

  real(dp),    parameter :: molarMassCO2_kg     = 4.4E-2_dp   !! Mass of 1 mol CO2 in kg 


  ! --- abreviations -----------

  logical, PARAMETER, PUBLIC :: C4 = .true.     ! Constant for C4 photosynthetic pathway of vegetation types
  logical, PARAMETER, PUBLIC :: C3 = .false.    ! Constant for C3 photosynthetic pathway of vegetation types

END MODULE mo_jsbach_constants
