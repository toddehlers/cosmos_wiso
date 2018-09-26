   MODULE mo_control_add

!***********************************************************************
!
!**** *MODULE mo_control_add* - control variables for add modules.
!
!    
!     
!     Purpose
!     -------
!     - declaration
!
!
!**********************************************************************
      implicit none

      
! Logical unit number for I/O.

      INTEGER :: io_stdo_add = 16       !  standard out.
      INTEGER :: io_stdi_add = 15       !  standard in.

      INTEGER :: io_rsti_add = 37      !  restart in. 
      INTEGER :: io_rsto_add = 37      !  restart out. 
      
      INTEGER :: io_in_gpre16 = 410
      INTEGER :: io_in_gpre18 = 411
      INTEGER :: io_in_gpreD = 414
      INTEGER :: io_in_rval16 = 412
      INTEGER :: io_in_rval18 = 413
      INTEGER :: io_in_rvalD = 415
      

! Control variables

      REAL    :: dt_add            !  time step length [sec].
      REAL    :: dtb_add           !  time step length [days].
      INTEGER :: ndtday_add        !  time steps per day.

      INTEGER :: ldtadd           !  time step number from add restart file
      INTEGER :: ldtrunadd        !  actual time steps of run.


      INTEGER :: icycliadd        !  switch for cyclicity.
      INTEGER :: ndtrunadd        !  total no. of time steps of run.

      INTEGER :: addstartyear            !  year of ocean restart file
      INTEGER :: addstartmonth           !  month of ocean restart file
      INTEGER :: addstartday             !  day of ocean restart file

! MPIOM is using variable lyear already !



   

      REAL    :: rmasko = -9.e33   !  value at wet/?dry cells in ocean. 
                                   ! will overwrite value from namelist!

      INTEGER :: kchck_add = 0          !  switch for extended print control (0/1). 
      
      END MODULE mo_control_add
