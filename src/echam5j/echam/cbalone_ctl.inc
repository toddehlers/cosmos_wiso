!*    NAMELIST *CBALONE_CTL* To control (stand alone) carbon balance
!                                                                 
!      
!    NAME     TYPE      PURPOSE                                       
!    ----     ----      -------                                       
!                                                                       
!   *DRIVER_DATA_PATH*      CHARACTER    path to the driver data
!   *EXPERIMENT*            CHARACTER    name of the experiment providing driving data
!   *SPECIFICNAME*          CHARACTER    driving data identifier: filenames are
!                                        assumed to be EXPERIMENT_yyyymm.SPECIFICNAME.nc
!   *CLIMATE_YEARSTART*     INTEGER      first year of the driving data 
!   *CLIMATE_YEAREND*       INTEGER      last year of the driving data
!   *RUN_YEAR_FIRST*        INTEGER      first year of the cbalance experiment
!   *RUN_YEAR_LAST*         INTEGER      last year of the cbalance experiment
!   *OUT_INTERVAL*          CHARACTER    output interval: YEAR / MONTH / DAY
!   *INPUT_SCALING*         LOGICAL      true for CO2-scaling
!   *CO2_FILE*              CHARACTER    CO2 development file - only with input_scaling
!   *CLIMATOLOGY_DIFF_FILE* CHARACTER    Climatological differences of driver data
!                                         - only with input_scaling
!   *REF_YEAR_PAST*         INTEGER      reference years for input scaling
!   *REF_YEAR_RECENT*       INTEGER
! 
! !compare jsbach_ctl
!   *USE_DYNVEG*            LOGICAL      true to run the dynamic vegetation
!   *LCC*                   LOGICAL      Land cover change: read cover fractions
!                                        from external file
! !compare cbalance_ctl
!   *READ_CPOOLS*           LOGICAL      read inital values for carbon pools from file
!
! !compare dynveg_ctl
!   *READ_FPC*              LOGICAL      read initial FPCs from FPC_FILE_NAME 
!   *FPC_FILE_NAME*         CHARACTER    file containing initail values for FPC

NAMELIST /CBALONE_CTL/ &
      DRIVER_DATA_PATH, EXPERIMENT, SPECIFICNAME, CLIMATE_YEARSTART, &
      CLIMATE_YEAREND, RUN_YEAR_FIRST, RUN_YEAR_LAST, OUT_INTERVAL, &
      INPUT_SCALING, CO2_FILE, CLIMATOLOGY_DIFF_FILE, REF_YEAR_PAST, &
      REF_YEAR_RECENT, USE_DYNVEG, LCC, READ_CPOOLS, READ_FPC, FPC_FILE_NAME 
!                                                                       
!     ------------------------------------------------------------------
