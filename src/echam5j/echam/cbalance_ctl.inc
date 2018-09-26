!*    NAMELIST *CBALANCE_CTL* To control the carbon pools.
!                                                                 
!      
!   NAME            TYPE      PURPOSE                                       
!   ----            ----      -------                                       
!                                                                       
!  *READ_CPOOLS*     LOGICAL    Inintialize carbon pools with data from an 
!                               external file 
!  *CPOOL_FILE_NAME* CHARACTER  Name of the file containing initial data for
!                               carbon pools (only used if read_cpools=true)

NAMELIST /CBALANCE_CTL/ &
      READ_CPOOLS, CPOOL_FILE_NAME
!                                                                       
!     ------------------------------------------------------------------
