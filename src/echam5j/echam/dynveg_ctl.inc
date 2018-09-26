!*    NAMELIST *DYNVEG_CTL* To control JSBACH dynamic vegetation.
!                                                                 
!      
!      NAME                TYPE       PURPOSE                                       
!      ----                ----       -------                                       
!                                                                       
!     *LCTLIB_DYNVEG*      CHARACTER  land cover type library file for the 
!                                     dynamic vegetation
!     *READ_FPC*           LOGICAL    read initial cover fraction from an 
!                                     external file
!     *FPC_FILE_NAME*      CHARACTER  name of the external vegetation file
!     *DYNVEG_ALL*         LOGICAL    activate competition between woody types
!                                     and grass
!     *DYNVEG_FEEDBACK*    LOGICAL    activate feedback of dynveg on jsbach
!     *INIT_RUNNING_MEANS* LOGICAL    initialize running means of climate in
!                                     the climate buffer
!     *READ_CLIMBUF*       LOGICAL    read initial climate variables from an
!                                     external file
!     *CLIMBUF_FILE_NAME*  CHARACTER  name of the external climate file

NAMELIST /DYNVEG_CTL/ &
      LCTLIB_DYNVEG, READ_FPC, FPC_FILE_NAME, DYNVEG_ALL, &
      DYNVEG_FEEDBACK, INIT_RUNNING_MEANS, READ_CLIMBUF, CLIMBUF_FILE_NAME
!                                                                       
!     ------------------------------------------------------------------
