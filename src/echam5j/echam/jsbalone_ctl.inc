!*    NAMELIST *JSBALONE_CTL* To control the land surface scheme JSBACH.
!
!      This namelist is read only in the stand-alone jsbach experiments
!                                                                 
!      
!      NAME     TYPE      PURPOSE                                       
!      ----     ----      -------                                       
!                                                                       
!     *OUT_EXPNAME*      CHARACTER   Experiment name 
!     *LRESUME*          LOGICAL     Restarted run?
!     *DT_START*         INTEGER     Start date of experiment
!                                      Vector: (yr, mo, dy, hr, min, sec)
!     *DT_RESUME*        INTEGER     Restart date of experiment
!                                      Vector: (yr, mo, dy, hr, min, sec)
!     *DT_STOP*          INTEGER     Stop date of experiment
!                                      Vector: (yr, mo, dy, hr, min, sec)
!     *NO_DAYS*          INTEGER     Number of days to run the experiment for
!     *NO_STEPS*         INTEGER     Number of time steps to run the experiment for
!     *NO_CYCLES*        INTEGER     Number of restart cycles
!     *DELTA_TIME*       INTEGER     Model time step
!     *LTIMER*           LOGICAL     Check model performance using the timer?
!     *PUTRERUN*      IO_TIME_EVENT  Periode of restart writing
!     *PUTDATA*       IO_TIME_EVENT  Periode of writing output
!     *TRIGFILES*     IO_TIME_EVENT  Periode of output file generation 
!
NAMELIST /JSBALONE_CTL/ &
      OUT_EXPNAME, LRESUME, DT_START, DT_RESUME, DT_STOP, NO_DAYS, NO_STEPS, &
      NO_CYCLES, DELTA_TIME, LTIMER, PUTRERUN, PUTDATA, TRIGFILES
!                                                                       
!     ------------------------------------------------------------------
