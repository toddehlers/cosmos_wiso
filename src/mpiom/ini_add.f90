
      SUBROUTINE INI_ADD(kpaufr,kpicycli,pdt,kpndtrun,kpie,kpje,kpke   &
     &           ,pddpo,ptho,psao,pdlxp,pdlyp,ptiestu,ptiestw          &
     &           ,kplyear,kplmonth,kplday,kpldtoce,pyears,pmonts       &
     &           ,pgila,pgiph)

!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/ini_bgc.f90,v $\\
!$Revision: 1.2.2.1.20.2 $\\
!$Date: 2006/04/03 11:27:49 $\\

!****************************************************************
!
!**** *INI_ADD* - initialize marine add module.
!
!     
!     
!     Purpose
!     -------
!    
!     - initialize add variables
!     
!     - read restart fields
!     - calculate budgets
!
!     Method
!     -------
!     - Noch zu tun : Biharmonic tracer diffusion (mit AULAPTS > ALMZER)
!                     Convective adjustment
!                     sea ice growth influence on tracer concentration
!
!**   Interface.
!     ----------
!
!     *CALL*       *INI_ADD(kpaufr,kpicycli,pdt,kpndtrun,kpie,kpje,kpke
!                          ,pddpo,ptho,psao,pdlxp,pdlyp,ptiestu
!                          ,kplyear,kplmonth,kplday,kpldtoce)*
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpaufr*   - 1/0 for read / do not read restart file
!     *INTEGER* *kpicycli* - flag for cyclicity.
!     *REAL*    *pdt*      - ocean model time step [sec].
!     *INTEGER* *kpndtrun* - total no. of time steps of run.
!     *INTEGER* *kpie*     - zonal dimension of model grid.
!     *INTEGER* *kpje*     - meridional dimension of model grid.
!     *INTEGER* *kpke*     - vertical dimension of model grid.
!     *REAL*    *pddpo*    - size of scalar grid cell (3rd REAL) [m].
!     *REAL*    *ptho*     - potential temperature [deg C].
!     *REAL*    *psao*     - salinity [psu.].
!     *REAL*    *pdlxp*    - size of scalar grid cell (zonal) [m].
!     *REAL*    *pdlyp*    - size of scalar grid cell (meridional) [m].
!     *REAL*    *pgila*   - geographical longitude of grid points [degree E].
!     *REAL*    *pgiph*   - geographical latitude  of grid points [degree N].
!     *REAL*    *ptiestu*  - depth of level [m].
!     *INTEGER* *kplyear*  - year  in ocean restart date
!     *INTEGER* *kplmonth* - month in ocean restart date
!     *INTEGER* *kplday*   - day   in ocean restart date
!     *INTEGER* *kpldtoce* - step  in ocean restart date
!
!     Externals
!     ---------
!     READ_NAMELIST_ADD CALL INI_TIMESER_ADD BELEG_ADD AUFR_ADD 
!     CHCK_ADD 
!
!**********************************************************************

      USE mo_contra     
      USE mo_addmean
      USE mo_control_add
      use mo_param1_add 
      use mo_mpi
      USE mo_parallel, ONLY: global_sum
      USE mo_commo1, ONLY: weto



      implicit none
      INTEGER :: pyears,pmonts,kpie,kpje,kpke
      INTEGER :: kplyear,kplmonth,kplday,kpldtoce
      INTEGER :: kpaufr,kpicycli,kpndtrun,k
      INTEGER :: ii,jj,kk
      
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: ptho (kpie,kpje,kpke)
      REAL :: psao (kpie,kpje,kpke)
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL :: pgila(kpie*2,kpje*2)
      REAL :: pgiph(kpie*2,kpje*2)
      REAL :: ptiestu(kpke+1),ptiestw(kpke+1)

      REAL :: pdt

#ifdef DMSASSIM
      REAL :: dmsrms_ref
      REAL :: assim_fac(2)
      INTEGER :: ndmsdr,idmsdr,i,dmssuc(5),dmsdir(5)
#endif      

          
      !                    
! Set control constants ( mo_control_add )
!
      dt_add = pdt                   !  time step length [sec]
      ndtday_add=NINT(86400./dt_add)  !  time steps per day [no.]
      dtb_add=1./ndtday_add          !  time step length [days]
      
      icycliadd = kpicycli
      ndtrunadd = kpndtrun

!
! Initialize some namelist parameters
!
      

     
      mean_3D_freq_add = 3
!
! Initialize time step counter of run.
!
      ldtrunadd = 0

!                        
! Set namelist defaults and read namelist
!
     
     
      CALL READ_NAMELIST_ADD

!                        
! Initialize add time series.
!
     

      CALL INI_TIMESER_ADD(kpke,ptiestw)

!
! Compute total ocean area
!
      totarea_add = 0.0
      DO ii=2,kpie-1
         DO jj=2,kpje-1
            totarea_add = totarea_add + pdlxp(ii,jj)*pdlyp(ii,jj)*weto(ii,jj,1)
         END DO
      END DO
      CALL global_sum(totarea_add)
!                        
! Initialize addmean
!
      meantime_3d_add = 0
      meantime_2d_add = 0


     
      if (mean_3D_freq_add.eq.1) nmeantime_3d_add = ndtrunadd*dtb_add
     
      if (mean_3D_freq_add.eq.2) nmeantime_3d_add = 12*pyears + pmonts
     
      if (mean_3D_freq_add.eq.3) nmeantime_3d_add = 1

     
      if (mean_3D_freq_add.eq.4) nmeantime_3d_add = 1

      
      CALL ALLOC_MEM_ADDMEAN(kpie,kpje,kpke)

! Allocate memory for conservative tracer flied
     
      CALL ALLOC_MEM_CONTRA(kpie,kpje,kpke)
      CALL ALLOC_MEM_FORC(kpie,kpje)

!                        
! Initialize ocean tracers.
! 
     WRITE(0,*) 'Call of beleg_add'

      CALL BELEG_ADD(kpie,kpje,kpke,pddpo)

!     WRITE(0,*) 'End of beleg_add'
     
!                        
! Read restart fields from restart file
!
      IF(kpaufr.eq.1) THEN
         CALL AUFR_ADD(kpie,kpje,kpke,pddpo,kplyear,kplmonth,   &
     &                 kplday,kpldtoce)
      ENDIF



! Open add output files
    
    
      CALL OPEN_ADDMEAN_3D(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,pgila,pgiph     &
     &                    ,ptiestu,kplyear,kplmonth,kplday,kpldtoce)
     



#ifdef DMSASSIM
       


     IF(p_pe == p_io) THEN

       OPEN(321,FILE='SUCHWURZ',ACCESS='SEQUENTIAL',FORM='FORMATTED')
       do i=1,100000
       read(321,3210,err=1466,end=1466) ndmsdr,idmsdr,dmssuc,dmsdir,      &
     &                                  dmsrms_ref,assim_fac,dmspar   
       write(io_stdo_add,*)'INI_ADD :',i,ndmsdr,idmsdr,dmssuc,dmsdir,      &
     &                                  dmsrms_ref,assim_fac,dmspar   
       enddo
  
1466   continue
3210  format(12i4,8e15.9)
      CLOSE(321)

     ENDIF
!      

!      
#endif  /*DMSASSIM*/
   
      


!
! Global inventory of all tracers
!

      CALL INVENTORY_ADD(kpie,kpje,kpke)



      RETURN
      END
