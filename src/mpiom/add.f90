      SUBROUTINE ADD(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,psictho,psicomo,&
     &ptiestu,pdpio,kplyear,kplmon,kplday,kmonlen,kldtmon,kldtday)
     
!**********************************************************************
!
!**** *ADD* - .
!
! 
!     Purpose
!     -------
!     - time step addtional conservative tracers
!
!     Method:
!     ------
!     kchck=1 is used to check max/min of add arrays on wet/dry cells.
!     Note: prosil is used only for k=1,2. It is adressed, however, for
!           k=1,4 in loop 100. To save memory, 
!
!     *CALL*  *ADD(kpie,kpje,kpke,pfswr,ptho,pddpo,pdlxp,pdlyp,pdpio)*
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *PARAM1_add.h*        - .
!     
!     *COMMON*     *UNITS_ADD.h*        - .
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st REAL of model grid.
!     *INTEGER* *kpje*    - 2nd REAL of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) REAL of model grid.
!     *REAL*    *pfswr*   - solar radiation [W/m**2].
!     *REAL*    *psictho*  - sea ice thickness [m]
!     *REAL*    *psicomo*  - sea ice compactness [fra]
!     *REAL*    *ptho*    - potential temperature [deg C].
!     *REAL*    *psao*    - salinity [psu.].
!     *REAL*    *ppao*    - sea level pressure [Pascal].
!     *REAL*    *pddpo*   - size of scalar grid cell (depth) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (longitudinal) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (latitudinal) [m].
!     *REAL*    *pdpio*   - inverse size of grid cell (depth)[m].
!     *INTEGER* *kldtmon* - time step within current month (from mpi-om) (1=1st time step of month, 2=2nd etc)
!     *INTEGER* *kldtday* - time step within current day (from mpi-om)
!
!     Externals
!     ---------
!     OCPROD 
!**********************************************************************

      USE mo_contra
      
      
      USE mo_addmean
      USE mo_control_add
      USE mo_timeser_add
      use mo_param1_add
      USE mo_parallel
     
      use mo_commo1, only : dt
 
#ifdef ADDCONTRA
#if defined (__coupled)
      use mo_fluxes1, only: aoflfrwo16, aoflfrwo18, aoflfrwhdo, &
                            aoflfrio16, aoflfrio18, aoflfrihdo
#endif
      
implicit none
      
      INTEGER :: kpie,kpje,kpke
      INTEGER :: i, j
      

      
      REAL,intent(in) :: pddpo(kpie,kpje,kpke)
      REAL,intent(in) :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL,intent(in) :: pdpio(kpie,kpje,kpke)          
      REAL,intent(in) :: ptiestu(kpke+1)
      REAL,intent(in) :: psictho(kpie,kpje),psicomo(kpie,kpje)
      REAL :: h2o18_tep(kpie,kpje),h2o16_tep(kpie,kpje)
      REAL :: hDo16_tep(kpie,kpje)
      REAL :: io18_tep(kpie,kpje),io16_tep(kpie,kpje)
      REAL :: iD_tep(kpie,kpje)
      REAL :: o18_tep(kpie,kpje),o16_tep(kpie,kpje)
      REAL :: D_tep(kpie,kpje)
      REAL :: h2o18_pre(kpie,kpje),h2o16_pre(kpie,kpje)
      REAL :: hDo16_pre(kpie,kpje)
      REAL :: h2o18_rff(kpie,kpje),h2o16_rff(kpie,kpje)
      REAL :: hDo16_rff(kpie,kpje)
      REAL :: Ri_o18(kpie,kpje),Rw_o18(kpie,kpje)
      REAL :: Ri_D(kpie,kpje),Rw_D(kpie,kpje)
      REAL :: mice(kpie,kpje)
      REAL :: o18mice(kpie,kpje),o16mice(kpie,kpje),Dmice(kpie,kpje)
      REAL :: mice_tmp(kpie,kpje), oce_tmp(kpie,kpje)
      INTEGER ,intent(in) :: kplyear,kplmon,kplday,kmonlen,kldtmon,kldtday

      REAL :: sum_o16, sum_o18, sum_hdo, sum_vol, vol
      REAL :: zdelta, znumer, zdenom
      REAL :: zwisomin, zwisosec
      
      integer :: ndtday 
      NDTDAY=NINT(86400./DT)
     
      zwisomin= 1.e-6_dp
      zwisosec= 1.e-8_dp

!
! Increment add time step counter of run (initialized in INI_ADD).
!
      ldtrunadd = ldtrunadd + 1
!
! Increment bgc time step counter of experiment (initialized if IAUFR=0).
!
      ldtadd = ldtadd + 1
!
! Increment timeseries-1 sample counter (initialized in INI_TIMSER_ADD).
!

!
!--------------------------------------------------------------------

#ifdef PBGC_CK_TIMESTEP     ! ADD share other cppflags with PBGC
      WRITE(io_stdo_add,*)' '
      WRITE(io_stdo_add,*)'before ADD: call INVENTORY'
      call INVENTORY_ADD(kpie,kpje,kpke)
#endif
       call maschk_add(kpie,kpje,kpke,-1)
!
!-------------------------------------------------------------------!

!-------------------------------------------------------------------      
! First timestep of the day --> kldtday eq 1     
     
     meancnt_add_2D = meancnt_add_2D + 1     
     meancnt_add_3D = meancnt_add_3D + 1 

!------------------------------------------------------------------
!
! Add the boundary forcing (precipitation = P + E ; runoff)

!     WRITE(0,*) 'Max SO18 BFF ', MAXVAL(ocectra(:,:,1,ih2o18)), 'Min SO18 BFF ', MINVAL(ocectra(:,:,1,ih2o18)) 
!     WRITE(0,*) 'Max SO16 BFF ', MAXVAL(ocectra(:,:,1,ih2o16)), 'Min SO16 BFF ', MINVAL(ocectra(:,:,1,ih2o16))     

     sum_vol = 0.
     sum_o16 = 0.
     sum_o18 = 0.
     sum_hdo = 0.
          
     DO j=1,kpje
      DO i=1,kpie
           h2o18_pre(i,j) = 0.
           hDo16_pre(i,j) = 0.
           h2o16_pre(i,j) = 0.
           h2o18_rff(i,j) = 0.
           hDo16_rff(i,j) = 0.
           h2o16_rff(i,j) = 0. 
           o18_tep(i,j) = 0.
           D_tep(i,j) = 0.
           o16_tep(i,j) = 0.
           IF(pddpo(i,j,1).GT.0.5) THEN
              vol = pddpo(i,j,1)*pdlxp(i,j)*pdlyp(i,j)
              sum_vol = sum_vol + vol

#if defined (__coupled)
!MW for coupled model: convert fresh water flux from m/s to kg/m**3, river runoff already included in fresh water flux 
! - sum up both fresh water fluxes here as storage of snow isotope values on sea ice is ignored in GROWTH routine, so far
! - take into consideration the different definition of SMOW in ECHAM5 as compared to MPIOM (see ECHAM5 routine setwiso.f90)

              h2o18_pre(i,j) = ((aoflfrwo18(i,j)+aoflfrio18(i,j))*1000.*DT/pddpo(i,j,1))/((20./18.)*100.)
              hDo16_pre(i,j) = ((aoflfrwhdo(i,j)+aoflfrihdo(i,j))*1000.*DT/pddpo(i,j,1))/((19./18.)*2.*1000.) 
              h2o16_pre(i,j) = (aoflfrwo16(i,j)+aoflfrio16(i,j))*1000.*DT/pddpo(i,j,1)

              sum_o16 = sum_o16 + h2o16_pre(i,j)*vol
              sum_o18 = sum_o18 + h2o18_pre(i,j)*vol
              sum_hdo = sum_hdo + hDo16_pre(i,j)*vol
         
#else
!MW for uncoupled model: convert fresh water flux from m/s to kg/m**3, and convert river runoff from m**3/s to kg/m**3 
              h2o18_pre(i,j) = pre18(i,j)*1000.*DT/pddpo(i,j,1)
              hDo16_pre(i,j) = preD(i,j)*1000.*DT/pddpo(i,j,1)
              h2o16_pre(i,j) = pre16(i,j)*1000.*DT/pddpo(i,j,1)
              h2o18_rff(i,j) = (riv18(i,j)*1000.*DT)/vol
              hDo16_rff(i,j) = (rivD(i,j)*1000.*DT)/vol
              h2o16_rff(i,j) = (riv16(i,j)*1000.*DT)/vol
              
              sum_o16 = sum_o16 + (h2o16_pre(i,j)+h2o16_rff(i,j))*vol
              sum_o18 = sum_o18 + (h2o18_pre(i,j)+h2o18_rff(i,j))*vol
              sum_hdo = sum_hdo + (hDo16_pre(i,j)+hDo16_rff(i,j))*vol
         
#endif

          ENDIF
       ENDDO
     ENDDO

     CALL global_sum(sum_vol)
     CALL global_sum(sum_o16)
     CALL global_sum(sum_o18)
     CALL global_sum(sum_hdo)

     DO j=1,kpje
      DO i=1,kpie
         IF(pddpo(i,j,1).GT.0.5) THEN
#if defined (__coupled)

            o16_tep(i,j) = h2o16_pre(i,j)-sum_o16/sum_vol        

            zdelta=2.0052/1000.0
            znumer=h2o18_pre(i,j)
            zdenom=h2o16_pre(i,j)
            IF (abs(zdenom).gt.zwisomin) zdelta = znumer/zdenom
            IF (abs(1.-zdelta).lt.zwisosec) zdelta=1.  ! cut off rounding errors
            o18_tep(i,j) = zdelta * o16_tep(i,j)

            zdelta=0.15576/1000.0
            znumer=hDo16_pre(i,j)
            zdenom=h2o16_pre(i,j)
            IF (abs(zdenom).gt.zwisomin) zdelta = znumer/zdenom
            IF (abs(1.-zdelta).lt.zwisosec) zdelta=1.  ! cut off rounding errors
            D_tep(i,j) = zdelta * o16_tep(i,j)

#else

            o16_tep(i,j) = (h2o16_pre(i,j)+h2o16_rff(i,j))-sum_o16/sum_vol        

            zdelta=2.0052/1000.0
            znumer=(h2o18_pre(i,j)+h2o18_rff(i,j))
            zdenom=(h2o16_pre(i,j)+h2o16_rff(i,j))
            IF (abs(zdenom).gt.zwisomin) zdelta = znumer/zdenom
            IF (abs(1.-zdelta).lt.zwisosec) zdelta=1.  ! cut off rounding errors
            o18_tep(i,j) = zdelta * o16_tep(i,j)

            zdelta=0.15576/1000.0
            znumer=(hDo16_pre(i,j)+hDo16_rff(i,j))
            zdenom=(h2o16_pre(i,j)+h2o16_rff(i,j))
            IF (abs(zdenom).gt.zwisomin) zdelta = znumer/zdenom
            IF (abs(1.-zdelta).lt.zwisosec) zdelta=1.  ! cut off rounding errors
            D_tep(i,j) = zdelta * o16_tep(i,j)

#endif         
            ocectra(i,j,1,ih2o18) = ocectra(i,j,1,ih2o18)+o18_tep(i,j)
            ocectra(i,j,1,ihDo16) = ocectra(i,j,1,ihDo16)+D_tep(i,j)
            ocectra(i,j,1,ih2o16) = ocectra(i,j,1,ih2o16)+o16_tep(i,j)
          
          ENDIF
       ENDDO
     ENDDO
     
!     WRITE(0,*) 'Max Pre18 ', MAXVAL(pre18), 'Min Pre18 ', MINVAL(pre18) 
!     WRITE(0,*) 'Max Pre16 ', MAXVAL(pre16), 'Min Pre16 ', MINVAL(pre16)
!     WRITE(0,*) 'Max Riv18 ', MAXVAL(riv18), 'Min Riv18 ', MINVAL(riv18) 
!     WRITE(0,*) 'Max Riv16 ', MAXVAL(riv16), 'Min Riv16 ', MINVAL(riv16) 
!     WRITE(0,*) 'Max T_Pre18 ', MAXVAL(h2o18_pre), 'Min T_Pre18 ', MINVAL(h2o18_pre) 
!     WRITE(0,*) 'Max T_Pre16 ', MAXVAL(h2o16_pre), 'Min T_Pre16 ', MINVAL(h2o16_pre)     
!     WRITE(0,*) 'Max T_Riv18 ', MAXVAL(h2o18_rff), 'Min T_Riv18 ', MINVAL(h2o18_rff) 
!     WRITE(0,*) 'Max T_Riv16 ', MAXVAL(h2o16_rff), 'Min T_Riv16 ', MINVAL(h2o16_rff)


!     WRITE(0,*) 'Max SO18 BFICE ', MAXVAL(o18_tep), 'Min SO18 BFICE ', MINVAL(o18_tep) 
!     WRITE(0,*) 'Max SO16 BFICE ', MAXVAL(o16_tep), 'Min SO16 BFICE ', MINVAL(o16_tep)
!     WRITE(0,*) 'Max O18 BFICE ', MAXVAL(ocectra(:,:,1,ih2o18)), 'Min SO18 BFICE ', MINVAL(ocectra(:,:,1,ih2o18)) 
!     WRITE(0,*) 'Max O16 BFICE ', MAXVAL(ocectra(:,:,1,ih2o16)), 'Min SO16 BFICE ', MINVAL(ocectra(:,:,1,ih2o16))
      
!-------------------------------------------------------------------
!
! Calculate ice formation fractionation
!  
      DO j=1,kpje
      DO i=1,kpie
         Rw_o18(i,j) = 0.
         Rw_D(i,j) = 0.
         Ri_o18(i,j) = 0.
         Ri_D(i,j) = 0. 
         mice(i,j) = 0.       
         o18mice(i,j) = 0.
         o16mice(i,j) = 0.
         Dmice(i,j) = 0.
         h2o16_tep(i,j) = 0.
         h2o18_tep(i,j) = 0.
         hDo16_tep(i,j) = 0.
         io16_tep(i,j) = 0.
         io18_tep(i,j) = 0.
         iD_tep(i,j) = 0.
     
         IF(pddpo(i,j,1).GT.0.5) THEN
	       Rw_o18(i,j) = ocectra(i,j,1,ih2o18)/ocectra(i,j,1,ih2o16)
           Rw_D(i,j) = ocectra(i,j,1,ihDo16)/ocectra(i,j,1,ih2o16)   
           Ri_o18(i,j) = 1.003*Rw_o18(i,j)
           Ri_D(i,j) = 1.01865*Rw_D(i,j)
           mice(i,j) = 910.*psictho(i,j)*psicomo(i,j)*pdlxp(i,j)*pdlyp(i,j)
           o18mice(i,j) = (mice(i,j)*(Ri_o18(i,j)/(Ri_o18(i,j)+Ri_D(i,j)+1)))/(pddpo(i,j,1)*pdlxp(i,j)*pdlyp(i,j))
           Dmice(i,j) = (mice(i,j)*(Ri_D(i,j)/(Ri_o18(i,j)+Ri_D(i,j)+1)))/(pddpo(i,j,1)*pdlxp(i,j)*pdlyp(i,j)) 
           o16mice(i,j) = (mice(i,j)/(Ri_o18(i,j)+Ri_D(i,j)+1))/(pddpo(i,j,1)*pdlxp(i,j)*pdlyp(i,j))
 
           mice_tmp(i,j) = mice(i,j)/(pddpo(i,j,1)*pdlxp(i,j)*pdlyp(i,j))
    
          

           IF((mice_tmp(i,j).LT.910.).AND.(o18mice(i,j).LT.1.82624).AND.(Dmice(i,j).LT.0.14407)) THEN             

              h2o16_tep(i,j) = o16mice(i,j) - icectra(i,j,1,ih2o16_ice) 
              io16_tep(i,j) = icectra(i,j,1,ih2o16_ice) + h2o16_tep(i,j)
              o16_tep(i,j) = ocectra(i,j,1,ih2o16) - h2o16_tep(i,j)
              ocectra(i,j,1,ih2o16) = o16_tep(i,j)
              icectra(i,j,1,ih2o16_ice) = io16_tep(i,j)
          
              h2o18_tep(i,j) = o18mice(i,j) - icectra(i,j,1,ih2o18_ice) 
              io18_tep(i,j) = icectra(i,j,1,ih2o18_ice) + h2o18_tep(i,j)
              o18_tep(i,j) = ocectra(i,j,1,ih2o18) - h2o18_tep(i,j)
              ocectra(i,j,1,ih2o18) = o18_tep(i,j)
              icectra(i,j,1,ih2o18_ice) = io18_tep(i,j)

              hDo16_tep(i,j) = Dmice(i,j) - icectra(i,j,1,ihDo16_ice) 
              iD_tep(i,j) = icectra(i,j,1,ihDo16_ice) + hDo16_tep(i,j)
              D_tep(i,j) = ocectra(i,j,1,ihDo16) - hDo16_tep(i,j)
              ocectra(i,j,1,ihDo16) = D_tep(i,j)
              icectra(i,j,1,ihDo16_ice) = iD_tep(i,j)
            ENDIF

         ENDIF
                      
      ENDDO
      ENDDO

#endif

!-----------------------------------------------------------------------
!       WRITE(0,*) 'Max of ice fra ', MAXVAL(psicomo)
      
!       WRITE(0,*) 'Max O18 ICE', MAXVAL(o18mice), 'Min O18 ICE ', MINVAL(o18mice)
!       WRITE(0,*) 'Max O16 ICE', MAXVAL(o16mice), 'Min O16 ICE ', MINVAL(o16mice)
!       WRITE(0,*) 'Max O18 AFICE', MAXVAL(h2o18_tep), 'Min O18 AFICE ', MINVAL(h2o18_tep)
!       WRITE(0,*) 'Max O16 AFICE', MAXVAL(h2o16_tep), 'Min O16 AFICE ', MINVAL(h2o16_tep)
!
!     Calculate the 3d mean field

!     WRITE(0,*) 'check pddpo max ', MAXVAL(pddpo(:,:,:))
      
      CALL OCPROD_ADD(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,pdpio,kplmon)
 
      call maschk_add(kpie,kpje,kpke,21)

    
     

      IF (MOD(ldtrunadd,nfreqts1).EQ.0) THEN
      
        CALL AVRG_TIMESER_ADD
      ENDIF




#ifdef PBGC_CK_TIMESTEP 
      WRITE(io_stdo_add,*)' '
      WRITE(io_stdo_add,*)'after ADD: call INVENTORY'
      call INVENTORY_ADD(kpie,kpje,kpke)
#endif	 
      call maschk_add(kpie,kpje,kpke,101)

      RETURN
      END
