!*    NAMELIST *ALBEDO_CTL* To control land surface albedo.
!                                                                 
!      
!      NAME     TYPE      PURPOSE                                       
!      ----     ----      -------                                       
!                                                                       
!     *USE_ALBEDOCANOPY*  LOGICAL  true: use map of canopy albedo 
!                                  false: use PFT specific albedo values 

NAMELIST /ALBEDO_CTL/ &
      USE_ALBEDOCANOPY
!                                                                       
!     ------------------------------------------------------------------
