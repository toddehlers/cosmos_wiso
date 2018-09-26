!*    NAMELIST *JSBACH_CTL* To control the land surface scheme JSBACH.
!                                                                 
!      
!      NAME     TYPE      PURPOSE                                       
!      ----     ----      -------                                       
!                                                                       
!     *NTILES*           INTEGER   Number of tiles defined on each grid cell
!     *STANDALONE*       LOGICAL   Type of model run 
!                                    standalone: stand-alone jsbach run
!                                    coupled: jsbach driven by an atmosphere model
!     *COUPLING*         CHARACTER Type of coupling: explicit, semi, implicit
!     *LSS*              CHARACTER Land surface sceme: ECHAM
!     *USE_BETHY*        LOGICAL   Use BETHY model (photosynthesis, respiration)
!     *USE_PHENOLOGY*    LOGICAL   Calculate LAI using the phenology module
!     *USE_ALBEDO*       LOGICAL   Calculate albedo depending on vegetaion
!     *USE_DYNVEG*       LOGICAL   Use the dynamic vegetation module
!     *LCC*              LOGICAL   Read land cover maps from external files?
!     *FILE_TYPE*        CHARACTER Output format: grib or netcdf
!     *OUT_STATE*        LOGICAL      write the jsbach stream
!     *LPOST_ECHAM*      LOGICAL   write jsbach output variables, even if they
!                                  are includeded in the echam output
!     *MISSING_VALUE*    REAL      Missing value for the output
!     *TEST_STREAM*      LOGICAL   Additional stream for model testing
!     *DEBUG*            LOGICAL   Additional output for debugging
!     *GRID_FILE*        CHARACTER File containing grid information
!     *LCTLIB_FILE*      CHARACTER Name of the land cover library file
!     *VEG_FILE*         CHARACTER File containig initial data for the vegetation 
!     *SURF_FILE*        CHARACTER File containig initial data for the surface
!     *SOIL_FILE*        CHARACTER File containig initial data for the soil
!     *READ_COVER_FRACT* LOGICAL   read cover fractions from file

NAMELIST /JSBACH_CTL/ &
      NTILES, STANDALONE, LSS, USE_BETHY, USE_PHENOLOGY, USE_ALBEDO, &
      USE_DYNVEG, LCC, FILE_TYPE, OUT_STATE, LPOST_ECHAM, &
      READ_COVER_FRACT, MISSING_VALUE, TEST_STREAM, DEBUG, GRID_FILE, &
      LCTLIB_FILE, VEG_FILE, SURF_FILE, SOIL_FILE
!                                                                       
!     ------------------------------------------------------------------
