!*    NAMELIST *SOIL_CTL* To control jsbach soil module.
!                                                                 
!      
!   NAME               TYPE      PURPOSE                                       
!   ----               ----      -------                                       
!                                                                       
!  *SKIN_RES_MAX*      REAL    Maximum content of skin reservoir of bare soil [m]
!  *MOIST_CRIT_FRACT*  REAL    Critical value of soil moisutre above which 
!                              transpiration is not affected by the moisture 
!                              stress in the soil, expressed as fraction of the 
!                              maximum soil moisture content
!  *MOIST_WILT_FRACT*  REAL    Soil moisture content at permanent wilting point,
!                              expressed as fraction of maximum soil moisture 
!                              content
!  *MOIST_MAX_LIMIT*   REAL    Upper limit for maximum soil moisture content: 
!                              If positive, max_moisture from ini file is cut 
!                              off at this value.
!  *CRIT_SNOW_DEPTH*   REAL    Critical snow depth for correction of surface 
!                              temperature for melt [m water equivalent] 

NAMELIST /SOIL_CTL/ &
      SKIN_RES_MAX, MOIST_CRIT_FRACT, MOIST_WILT_FRACT, MOIST_MAX_LIMIT, &
      CRIT_SNOW_DEPTH
!                                                                       
!     ------------------------------------------------------------------
