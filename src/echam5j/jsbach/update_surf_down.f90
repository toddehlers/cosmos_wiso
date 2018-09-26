SUBROUTINE update_surf_down(                    &    
       kbdim                                    &
     , pqm1, ptsl, pspeed10, ptm1               &
     , pwl, pcvw, pwlmx                         &
     , psn, pcvs, psnc, pgld                    &
     , pgrndcapc                                &
!!$  , paphm1, paphms1, ptte                       ! ATTENTION!
     , pevapl, pevapot, pevwsd                  &
     , prain, psnow                             &
     , lmask, lpglac                            &
     , palac                                    &
     , psnacl, psnmel, progl                    &
     , pmlres, tte_corr                         &
!---wiso-code
     , lwisofracl                               &
     , pwisosn, pwisosnc, pwisogld              &
     , psnglac,    pwisosnglac                  &
     , pwisosnmel, pwisowl                      &
     , pwisoevapl, pwisoevapot, pwisoevwsd      &
     , pwisorain, pwisosnow, pwisorogl          &
     , pwisoalac, pwisosnacl                    &
     , pwisomlres                               &
!---wiso-code-end
     ) 
  !
  !     ------------------------------------------------------------------
  ! ptsl          Sfc temp [K] at timestep t+dt (unfiltered)
  ! pspeed10      Wind speed at 10m height [m/s] ... from 'vdiff'
  ! ptm1          Air temp [K] at timestep t-dt (filtered)
  ! pwl           Water content [m] in skin reservoir
  !               (vegetation and bare soil)
  ! pcvw          Skin reservoir fraction (= pwl/pwlmx, see *vdiff*)
  ! pwlmx         Skin reservoir [m] (calculated in *vdiff* as a function
  !               of vegetation index and leaf area index)
  ! psn           Snow depth [m water equivalent] at the ground
  ! pcvs          Fractional snow cover (function of psn in *physc*)
  ! psnc          Snow depth [m water equivalent] at the canopy
  ! pgld          Glacier depth (including snow) [m water equivalent]
  ! pgrndcapc     Heat capacity of the uppermost soil layer [j/m**2/K]
  ! pevapl        Total evapotranspiration, including sublimation [kg/m**2/s]
  ! pevapot       Potential evaporation/sublimation [kg/m**2/s]
  ! pevwsd        Evapotranspiration without sublimation and evaporation from interception reservoir [kg/m**2]
  ! prain         Total rainfall [kg/m**2/s]
  ! psnow         Total snowfall [kg/m**2/s]
  ! lpglac        Logical glacier mask
  ! palac         Precipitation minus sublimation at glacier points [m water equivalent]
  ! psnacl        Snow budget at non-glacier points [kg/m**2] (accumul.)
  ! psnmel        Snow melt [kg/m**2] (accumulated for diagnostics)
  ! progl         Glacier runoff (rain+snow/ice melt) [kg/m**2] (accumul.)
  !
!---wiso-code
  !
  ! psnglac       Snow depth [m water equivalent] on glaciers
  !
  ! pwisowl	tracer for Water content [m] in skin reservoir
  !               (vegetation and bare soil)
  ! pwisosn	    tracer for Snow depth [m water equivalent] at the ground
  ! pwisosnglac tracer for Snow depth [m water equivalent] on glaciers
  ! pwisorain	tracer for total rainfall [kg/m**2/s]
  ! pwisosnow	tracer for total snowfall [kg/m**2/s]
  ! pwisoevapl	tracer for Total evaporation, including sublimation [kg/m**2/s]
  ! pwisoevapot	tracer for Potential evaporation/sublimation [kg/m**2/s]
  ! pwisoevwsd  tracer for Evapotranspiration without sublimation and evaporation from interception reservoir [kg/m**2]
  ! pwisosnc    tracer for snow depth at the cancopy
  ! pwisogld    tracer for Glacier depth (including snow) [m water equivalent]
  ! pwisoalac   tracer for Precipitation minus sublimation at glacier points
  ! pwisosnmel  tracer for Snow melt
  ! pwisorogl   tracer for Glacier runoff (rain+snow/ice melt)
  ! pwisosnacl  tracer for snow budget at non-glacier points
  !
!---wiso-code-end
  !
  ! The following local variables represent the respective fluxes
  ! integrated over one timestep (delta_time) and divided by the density of
  ! water (rhoh2o). Units are m water equivalent.
  !
  ! zraind        Total rain
  ! zsnowd        Total snowfall
  ! zevttd        Total evaporation
  ! zevsnd        Sublimation from snow
  ! zevwld        Evaporation from the skin reservoir
  ! zrogl         Runoff at glacier points (rain and melt, but no calving)
  ! zsnmel        Snow/ice melt at land and glacier points
  ! zsncmelt      Snow melt in the canopy
  ! zsn           Snow budget at non-glacier points (snowfall-subl-melt)
  ! pmlres        Residual melt water available for infiltration into the
  !               non-frozen soil after being intercepted by the
  !               skin reservoir
  !
!---wiso-code
  !
  ! zwisoraind	tracer for Total rain
  ! zwisosnowd	tracer for Total snow
  ! zwisoevttd	tracer for Total evaporation
  ! zwisoevsnd	tracer for Sublimation from snow
  ! zwisoevwld	tracer for Evaporation from the skin reservoir
  ! zwisosnmel	tracer for Snow/ice melt at land and glacier points
  ! zwisosncmelt	tracer for Snow melt in the canopy
  ! zwisosn     tracer for snow budget at non-glacier points
  ! pwisomlres  tracer for Residual melt water available for infiltration into the
  !             non-frozen soil after being intercepted by the skin reservoir
  ! zwisorogl   tracer for Runoff at glacier points (rain and melt, but no calving)
  ! zgld_tmp    temporary file for Glacier depth
  ! zsn_tmp     temporary file for snow depth
  !
!---wiso-code-end
  !
  !     *SURF* - UPDATES LAND VALUES OF TEMPERATURE, MOISTURE AND SNOW.
  !              CALCULATE FLUXES OF TOTAL RAIN, TOTAL SNOW AND EVAPO-
  !              RATION FROM THE THREE RESERVOIRS (SN, WS, WL)
  !              CONVERT FLUXES (KG/M**2*S) TO CHANGES OF WATER LEVELS (M)
  !              DURING TIMESTEP DELTA_TIME.
  !
  !     J.F.GELEYN     E.C.M.W.F.     08/06/82.
  !     MODIFIED BY
  !     C.A.BLONDIN    E.C.M.W.F.    18/12/86.
  !     MODIFIED BY L.DUMENIL      MET.INST.HH     20/05/88
  !     J.-P. SCHULZ   MPI - 1997 : IMPLEMENTATION OF IMPLICIT
  !                                 COUPLING BETWEEN LAND SURFACE
  !                                 AND ATMOSPHERE.
  !     MODIFIED BY E. ROECKNER    MPI - SEPT 1998
  !     MODIFIED BY M. ESCH        MPI - APR  1999
  !     MODIFIED BY E. ROECKNER    MPI - JAN  2001
  !     MODIFIED BY I. Kirchner    MPI - MARCH 2001 date/time control
  !     MODIFIED BY E. ROECKNER    MPI - SEP  2002 interception reservoir 
  !                                                for snow changed
  !     MODIFIED BY L. KORNBLUEH   MPI - JAN  2003 removed MERGE
  !
  !     MODIFICATION
  !
  !     PURPOSE
  !
  !     INTERFACE.
  !
  !          *SURF* IS CALLED FROM *PHYSC*.
  !
  !     METHOD.
  !
  !     EXTERNALS.
  !
  !          NONE.
  !
  !     REFERENCE.
  !
  !          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
  !     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.
  !
  USE mo_parameters
  USE mo_control
  USE mo_param_switches
  USE mo_physc2
  USE mo_constants
  USE mo_vegetation
  USE mo_radiation
  USE mo_time_control,   ONLY : lstart, delta_time
  USE mo_kind,           ONLY : dp

!---wiso-code
  USE mo_wiso,           ONLY : lwiso, nwiso, twisorhoh2o, tnat, snglacmn, snglacmx
!---wiso-code-end

  IMPLICIT NONE

  INTEGER :: kbdim, jl
  REAL(dp) ::                                                               &
       ptsl(kbdim)                                                          &
       , pwl(kbdim)                                                         &
       , psn(kbdim),          pgld(kbdim)                                   &
       , psnc(kbdim),         pspeed10(kbdim)                               &
       , pgrndcapc(kbdim),    pmlres(kbdim)
  REAL(dp) ::                                                               &
       ptm1(kbdim)                                                          &
       , pqm1(kbdim),    tte_corr(kbdim)
       !!$    paphm1(kbdim), paphms1(kbdim) , ptte(kbdim)                   &
  REAL(dp) ::                                                               &
       pcvs(kbdim),      pcvw(kbdim),     pwlmx(kbdim)                      &
       , pevapl(kbdim),  pevapot(kbdim),  pevwsd(kbdim)                     &
       , prain(kbdim),   psnow(kbdim)                                       &
       , palac(kbdim),   psnacl(kbdim),   psnmel(kbdim),  progl(kbdim)
  LOGICAL ::                                                                &
       lmask(kbdim),     lpglac(kbdim)
  
  !  local arrays
  REAL(dp) ::                                                               &
       zraind (kbdim),   zsnowd(kbdim),   zevttd(kbdim)                     &
       , zevsnd(kbdim),  zevwld(kbdim)                                      &
       , zrogl(kbdim)                                                       &
       , zsnmel(kbdim),  zsn(kbdim)                                         &
       , zsncmelt(kbdim)
!!$    , zdp(kbdim),          
  REAL(dp) :: zlfdcp(kbdim)
  
  !  local scalars
  REAL(dp) ::                                                               &
       zsmelt, zsnmlt,                                                      &
       zmprcp,                                                              &
       zc2, zc3, zsncp, zexpt, zexpw, zsncmax, zsncwind
  REAL(dp) :: zdtime, zrcp, zsncfac

!---wiso-code
  ! declaration parameters water isotopes
  INTEGER :: jt
  INTEGER, OPTIONAL :: lwisofracl

  REAL(dp), OPTIONAL :: pwisorain(kbdim,nwiso),   pwisosnow(kbdim,nwiso)           &
     , pwisowl(kbdim,nwiso),    pwisosn(kbdim,nwiso)                               &
     , psnglac(kbdim),          pwisosnglac(kbdim,nwiso)                           &
     , pwisosnc(kbdim,nwiso),   pwisosnmel(kbdim,nwiso)

  REAL(dp), OPTIONAL :: pwisorogl(kbdim,nwiso)                                     &
     , pwisoevapl(kbdim,nwiso), pwisoevapot(kbdim,nwiso), pwisoevwsd(kbdim,nwiso)  &
     , pwisogld(kbdim,nwiso),   pwisosnacl(kbdim,nwiso)                            &
     , pwisoalac(kbdim,nwiso),  pwisomlres(kbdim,nwiso)

  ! water isotopes: local arrays
  REAL(dp) ::                                                                      &
       zwisoraind(kbdim,nwiso), zwisosnowd(kbdim,nwiso), zwisoevttd(kbdim,nwiso)   &
     , zwisoevsnd(kbdim,nwiso), zwisoevwld(kbdim,nwiso), zwisosncmelt(kbdim,nwiso) &
     , zwisorogl(kbdim,nwiso) , zwisosnmel(kbdim,nwiso), zwisosn(kbdim,nwiso)      &
     , zsnowd_tmp(kbdim)      , zmprcp_tmp(kbdim)      , zsncp_tmp(kbdim)          &
     , zevsnd_tmp(kbdim)      , zsubl_tmp(kbdim)       , zsnc_tmp(kbdim)           &

     , zsnc_tmp2(kbdim)       , zsncwind_tmp(kbdim)    , zwl_tmp(kbdim)            &
     , zsnmlt_tmp(kbdim)      , zsnglac_tmp(kbdim)                                 &
     , zgld_tmp(kbdim)        , zsn_tmp(kbdim)         , ztsl_tmp(kbdim)      

  LOGICAL :: lo_sub(kbdim), lo_zros(kbdim), lo_snglacmn(kbdim), lo_snglacmx(kbdim)

  !  water isotopes: local scalars
  REAL(dp) ::                                                                      &
       zwisosnmlt, zwisomprcp,                                                     &
       zwisosncp, zwisosubl,                                                       &
       zwisosncmax, zwisosncwind, zdelta,                                          &
       zwisomin, zwisosec
!---wiso-code-end

  ! Parameters
  zdtime = delta_time
  zc2=1.87E5_dp
  zc3=1.56E5_dp
  zsncfac=rhoh2o*g/zdtime
  
!---wiso-code
  IF (lwiso) THEN

  zwisomin = 1.e-12_dp
  zwisosec = 1.e-10_dp
  
  END IF
!---wiso-code-end

  !     ------------------------------------------------------------------
  !
  !*    1.     Convert water fluxes to [m water equivalent * timestep]
  !
  DO 110 jl=1,kbdim
     palac(jl)=0._dp
     zrogl(jl)=0._dp                    ! Runoff at glacier points (rain and melt, but no calving)
     pmlres(jl)=0._dp                   ! Residual melt water available for infiltration into the
                                        ! non-frozen soil after being intercepted by the
                                        ! skin reservoir
     zsn(jl)=0._dp                      ! Snow budget at non-glacier points (snowfall-subl-melt)
     zsnmel(jl)=0._dp                   ! Snow/ice melt at land and glacier points
     zsncmelt(jl)=0._dp                 ! Snow melt in the canopy
     zraind(jl)=prain(jl)*zdtime/rhoh2o ! Total rain
     zsnowd(jl)=psnow(jl)*zdtime/rhoh2o ! Total snowfall
     IF(lmask(jl)) THEN
        zrcp=1._dp/(cpd*(1._dp+vtmpc2*MAX(0.0_dp,pqm1(jl))))
        zlfdcp(jl)=alf*zrcp
        zevttd(jl)=pevapl(jl)*zdtime/rhoh2o                              ! Total evaporation
        zevsnd(jl)=pcvs(jl)*pevapot(jl)*zdtime/rhoh2o                    ! Sublimation from snow
        zevwld(jl)=(1._dp-pcvs(jl))*pcvw(jl)*pevapot(jl)*zdtime/rhoh2o   ! Evaporation from the skin reservoir
        pevwsd(jl)=zevttd(jl)-zevsnd(jl)                                 ! Evapotranspiration without sublimation and
                                                                         ! evaporation from interception reservoir [kg/m**2]
     END IF
110 END DO

!---wiso-code
 IF (lwiso) THEN

!   Convert water fluxes to [m water equivalent * timestep] - water isotopes
 DO jt=1,nwiso
   DO jl=1,kbdim
     zwisorogl(jl,jt)    =0._dp
     zwisosn(jl,jt)      =0._dp
     zwisosnmel(jl,jt)   =0._dp
     zwisosncmelt(jl,jt) =0._dp
     zwisoraind(jl,jt)   = pwisorain(jl,jt) * zdtime/twisorhoh2o(jt)
     zwisosnowd(jl,jt)   = pwisosnow(jl,jt) * zdtime/twisorhoh2o(jt)

     IF(lmask(jl)) THEN
!        IF (lwisofracl==0) THEN
!          calculate isotope ratio of snow reservoir
           zdelta=tnat(jt)
           IF (psn(jl).gt.zwisomin .and. pwisosn(jl,jt).gt.zwisomin) zdelta = pwisosn(jl,jt)/psn(jl)
           IF (ABS(1._dp-zdelta).lt.zwisosec) zdelta=1._dp  ! cut off rounding errors
           zwisoevsnd(jl,jt)=pcvs(jl)*zdelta*pwisoevapot(jl,jt)*zdtime/twisorhoh2o(jt)
!          calculate isotope ratio of skin reservoir
           zdelta=tnat(jt)
           IF (pwl(jl).gt.zwisomin .and. pwisowl(jl,jt).gt.zwisomin) zdelta = pwisowl(jl,jt)/pwl(jl)
           IF (ABS(1._dp-zdelta).lt.zwisosec) zdelta=1._dp  ! cut off rounding errors
           zwisoevwld(jl,jt)=(1._dp-pcvs(jl))*pcvw(jl)*zdelta*pwisoevapot(jl,jt)*zdtime/twisorhoh2o(jt)
!        ELSE
!!          calculate sublimation of snow reservoir
!           zwisoevsnd(jl,jt)=pcvs(jl)*pwisoevapot(jl,jt)*zdtime/twisorhoh2o(jt)
!!          calculate evaporation of skin reservoir
!           zwisoevwld(jl,jt)=(1._dp-pcvs(jl))*pcvw(jl)*pwisoevapot(jl,jt)*zdtime/twisorhoh2o(jt)
!        END IF
        zwisoevttd(jl,jt)=pwisoevapl(jl,jt)*zdtime/twisorhoh2o(jt)
        pwisoevwsd(jl,jt)=zwisoevttd(jl,jt)-zwisoevsnd(jl,jt)
     END IF
   END DO
 END DO
 
 END IF
!---wiso-code-end

  IF(lsurf) THEN

      !     ------------------------------------------------------------------
      !
      !*    2.     Budgets of snow (canopy and ground) and glaciers
      !
      !*    2.1    Snow changes in the canopy (interception of snowfall,
      !            sublimation, melting, unloading due to wind)
      !
      DO 210 jl=1,kbdim
        IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
           zsn(jl)=zsnowd(jl)+zevsnd(jl)
           zsncmax=MAX(0.0_dp,pwlmx(jl)-cwlmax)
           zsnowd_tmp(jl)=zsnowd(jl)                                          !---wiso-code: remember total sn
           zmprcp=MIN(zsnowd(jl)*cvinter,zsncmax-psnc(jl))
           zmprcp_tmp(jl)=zmprcp                                              !---wiso-code: remember interception
           zsncp=psnc(jl)+zmprcp
           zsnowd(jl)=zsnowd(jl)-zmprcp
           zsncp_tmp(jl)=zsncp                                                !---wiso-code: remember snow canopy
           zevsnd_tmp(jl)=zevsnd(jl)                                          !---wiso-code: remember atmospheric flux
           lo_sub(jl) = (zevsnd(jl).GT.0.0_dp)                                !---wiso-code: flag for sublimation vs. depositio
           psnc(jl)=MIN(MAX(0._dp,zsncp+zevsnd(jl)),zsncmax)
           zsubl_tmp(jl)=psnc(jl)-zsncp                                       !---wiso-code: remember sublimation
           zevsnd(jl)=zevsnd(jl)-(psnc(jl)-zsncp)
           zexpt=MAX(0._dp,ptm1(jl)+3._dp-tmelt)*zdtime/zc2
           zexpw=pspeed10(jl)*zdtime/zc3
           zsnc_tmp(jl)=psnc(jl)                                              !---wiso-code: remember new snow canopy
           zsncmelt(jl)=psnc(jl)*(1._dp-EXP(-zexpt))
           psnc(jl)=psnc(jl)-zsncmelt(jl)
           zsnc_tmp2(jl)=psnc(jl)                                             !---wiso-code: remember new snow canopy
           zsncwind=psnc(jl)*(1._dp-EXP(-zexpw))
           zsncwind_tmp(jl)=zsncwind                                          !---wiso-code: remember unloading due to wind
           psnc(jl)=psnc(jl)-zsncwind
           zsnowd(jl)=zsnowd(jl)+zsncwind
           tte_corr(jl)=zsncmelt(jl)*zsncfac*zlfdcp(jl)
           !
           !   pwl(jl)=pwl(jl)+zsncmelt(jl) see section 2.5
           !
        ELSE
           psnc(jl)=0._dp
        END IF
210  END DO


!---wiso-code
    IF (lwiso) THEN

    ! Snow changes in the canopy (interception of snowfall,
    ! sublimation, melting, unloading due to wind)
     DO jt=1,nwiso
         DO jl=1,kbdim
             IF (lmask(jl).AND..NOT.lpglac(jl)) THEN
                 zwisosn(jl,jt)=zwisosnowd(jl,jt)+zwisoevsnd(jl,jt)
                 zdelta=tnat(jt)
                 IF (zsnowd_tmp(jl).gt.zwisomin .and. zwisosnowd(jl,jt).gt.zwisomin) zdelta=zwisosnowd(jl,jt)/zsnowd_tmp(jl)
                 IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp   ! cut off rounding errors
                 zwisomprcp=zmprcp_tmp(jl) * zdelta
     	         zwisosncp=pwisosnc(jl,jt)+zwisomprcp
                 zwisosnowd(jl,jt)=zwisosnowd(jl,jt)-zwisomprcp
                 zdelta=tnat(jt)
                 IF (lo_sub(jl)) THEN                           ! deposition, assume identical delta value as the atmospheric flux
                    IF (zevsnd_tmp(jl).gt.zwisomin .and. zwisoevsnd(jl,jt).gt.zwisomin) zdelta=zwisoevsnd(jl,jt)/zevsnd_tmp(jl)
                 ELSE                                           ! sublimation, assume identical delta value as the canopy layer
                    IF (zsncp_tmp(jl).gt.zwisomin .and. zwisosncp.gt.zwisomin) zdelta=zwisosncp/zsncp_tmp(jl)
                 ENDIF  
                 IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp   ! cut off rounding errors
                 zwisosubl=zsubl_tmp(jl) * zdelta
                 pwisosnc(jl,jt)=zwisosncp+zwisosubl
                 zwisoevsnd(jl,jt)=zwisoevsnd(jl,jt)-zwisosubl
                 zdelta=tnat(jt)
                 IF (zsnc_tmp(jl).gt.zwisomin .and. pwisosnc(jl,jt).gt.zwisomin) zdelta=pwisosnc(jl,jt)/zsnc_tmp(jl)
                 IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp   ! cut off rounding errors
                 zwisosncmelt(jl,jt)=zsncmelt(jl) * zdelta
                 pwisosnc(jl,jt)=pwisosnc(jl,jt)-zwisosncmelt(jl,jt)
                 zdelta=tnat(jt)
                 IF (zsnc_tmp2(jl).gt.zwisomin .and. pwisosnc(jl,jt).gt.zwisomin) zdelta=pwisosnc(jl,jt)/zsnc_tmp2(jl)
                 IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp   ! cut off rounding errors
                 zwisosncwind=zsncwind_tmp(jl) * zdelta
                 pwisosnc(jl,jt)=pwisosnc(jl,jt)-zwisosncwind
                 zwisosnowd(jl,jt)=zwisosnowd(jl,jt)+zwisosncwind
             ELSE
                 pwisosnc(jl,jt)=0._dp
             END IF
         END DO
     END DO
     
    END IF
!---wiso-code-end

    !*    2.2    Snowfall and sublimation on land (excluding glaciers)
    !
     DO 220 jl=1,kbdim
        IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
           psn(jl)=psn(jl)+zsnowd(jl)+zevsnd(jl)
           zsn_tmp(jl)=psn(jl)                                       !---wiso-code: remember snow fall
           IF (psn(jl).LT.0._dp) THEN
              pevwsd(jl)=pevwsd(jl)+psn(jl)
              psn(jl)=0._dp
           END IF
        ELSE
           psn(jl)=0._dp
        END IF
220  END DO

!---wiso-code
    IF (lwiso) THEN

    ! Snowfall and sublimation on land (excluding glaciers) - water isotopes
     DO jt=1,nwiso
         DO jl=1,kbdim
             IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
                 pwisosn(jl,jt)=pwisosn(jl,jt)+zwisosnowd(jl,jt)+zwisoevsnd(jl,jt)
                 IF (zsn_tmp(jl).LT.0._dp) THEN
                     pwisoevwsd(jl,jt)=pwisoevwsd(jl,jt)+pwisosn(jl,jt)
                     pwisosn(jl,jt)=0._dp
                 END IF
             ELSE
                 pwisosn(jl,jt)=0._dp
             END IF
         END DO
     END DO
     
    END IF
!---wiso-code-end

    !*    2.3    Snowfall and sublimation on glaciers and diagnostics
    !
     DO 230 jl=1,kbdim
        IF (lpglac(jl)) THEN
           pgld(jl)=pgld(jl)+zsnowd(jl)+zevsnd(jl)
           palac(jl)=zraind(jl)+zsnowd(jl)+zevttd(jl)
           zrogl(jl)=zraind(jl)
        END IF
230  END DO

!---wiso-code
    IF (lwiso) THEN

     ! Snowfall and sublimation on glaciers and diagnostics - water isotopes
     DO jt=1,nwiso
         DO jl=1,kbdim
             IF (lpglac(jl)) THEN
                 pwisogld(jl,jt)=pwisogld(jl,jt)+zwisosnowd(jl,jt)+zwisoevsnd(jl,jt)
                 pwisoalac(jl,jt)=zwisoraind(jl,jt)+zwisosnowd(jl,jt)+zwisoevttd(jl,jt)
                 zwisorogl(jl,jt)=zwisoraind(jl,jt)
             END IF
         END DO
     END DO
     
     ! snow changes on glaciers - add precipitation and evaporation to reservoir snglac
     DO jl=1,kbdim
       lo_snglacmx(jl)=.FALSE.
       lo_snglacmn(jl)=.FALSE.
       IF (lpglac(jl)) THEN
         psnglac(jl)=psnglac(jl)+zsnowd(jl)+zraind(jl)
         IF (psnglac(jl).LT.snglacmn) THEN
            zsnglac_tmp(jl)=psnglac(jl)-zsnowd(jl)-zraind(jl)
            lo_snglacmn(jl)=.TRUE.
            psnglac(jl)=snglacmn
         END IF
         IF (psnglac(jl).GT.snglacmx) THEN
            zsnglac_tmp(jl)=psnglac(jl)
            lo_snglacmx(jl)=.TRUE.
            psnglac(jl)=snglacmx
         END IF
       ELSE
         psnglac(jl)=0._dp
       ENDIF
     END DO

     ! snow changes on glaciers - add snow fall to reservoir snglac - water isotopes
     DO jt=1,nwiso
       DO jl=1,kbdim
         IF (lpglac(jl)) THEN
           pwisosnglac(jl,jt)=pwisosnglac(jl,jt)+zwisosnowd(jl,jt)+zwisoraind(jl,jt)
           IF (lo_snglacmn(jl)) THEN
             pwisosnglac(jl,jt)=pwisosnglac(jl,jt)-zwisosnowd(jl,jt)-zwisoraind(jl,jt)
             zdelta=tnat(jt)
             IF (zsnglac_tmp(jl).gt.snglacmn) zdelta=pwisosnglac(jl,jt)/zsnglac_tmp(jl)
             IF (ABS(1._dp-zdelta).lt.zwisosec) zdelta=1._dp   ! cut off rounding errors
             pwisosnglac(jl,jt)=zdelta*snglacmn
           END IF
           IF (lo_snglacmx(jl)) THEN
             zdelta=pwisosnglac(jl,jt)/zsnglac_tmp(jl)
             IF (ABS(1._dp-zdelta).lt.zwisosec) zdelta=1._dp   ! cut off rounding errors
             pwisosnglac(jl,jt)=zdelta*snglacmx
           END IF
         ELSE
           pwisosnglac(jl,jt)=0._dp
         ENDIF
       END DO
     END DO

    END IF
!---wiso-code-end

    !
    !*    2.4    Snow and glacier melt
    !
     IF (.NOT. lstart) THEN
        DO 240 jl=1,kbdim
           zgld_tmp(jl)=pgld(jl)                    !---wiso-code: remember glacier depth
           ztsl_tmp(jl)=ptsl(jl)                    !---wiso-code: remember land surface temperatur
           zsn_tmp(jl)=psn(jl)                      !---wiso-code: remember snow depth
           IF (lmask(jl)) THEN
              IF (ptsl(jl).GT.tmelt) THEN
                 IF (lpglac(jl)) THEN
                    zsnmel(jl)=pgrndcapc(jl)*(ptsl(jl)-tmelt)/(alf*rhoh2o)
                    pgld(jl)=pgld(jl)-zsnmel(jl)
                    zrogl(jl)=zrogl(jl)+zsnmel(jl)
                    ptsl(jl)=tmelt
                 ELSE IF (psn(jl).GT.0._dp) THEN
                    zsmelt=pgrndcapc(jl)*(ptsl(jl)-tmelt)/(alf*rhoh2o)
                    zsnmel(jl)=MIN(psn(jl),zsmelt)
                    ptsl(jl)=ptsl(jl)-zsnmel(jl)*alf*rhoh2o/pgrndcapc(jl)
                    psn(jl)=MAX(psn(jl)-zsnmel(jl),0._dp)
                 END IF
              END IF
           END IF
240     END DO
     END IF

!---wiso-code
    IF (lwiso) THEN

    ! Snow and glacier melt - water isotopes
     IF (.NOT. lstart) THEN
         DO jt=1,nwiso
             DO jl=1,kbdim
                 IF (lmask(jl)) THEN
                     IF (ztsl_tmp(jl).GT.tmelt) THEN
                         IF (lpglac(jl)) THEN
                             ! melt part of snow (no additional fractionation)
                             zdelta=tnat(jt)
                             IF (zgld_tmp(jl).gt.zwisomin .and. pwisogld(jl,jt).gt.zwisomin) zdelta=pwisogld(jl,jt)/zgld_tmp(jl)
                             IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp ! cut off rounding errors
                             zwisosnmel(jl,jt)=zsnmel(jl) * zdelta
                             pwisogld(jl,jt)=pwisogld(jl,jt)-zwisosnmel(jl,jt)
                             zwisorogl(jl,jt)=zwisorogl(jl,jt)+zwisosnmel(jl,jt)
                         ELSE IF (zsn_tmp(jl).GT.0._dp) THEN
                             ! melt part of snow (no additional fractionation)
                             zdelta=tnat(jt)
                             IF (zsn_tmp(jl).gt.zwisomin .and. pwisosn(jl,jt).gt.zwisomin) zdelta=pwisosn(jl,jt)/zsn_tmp(jl)
                             IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp ! cut off rounding errors
                             zwisosnmel(jl,jt)=zsnmel(jl) * zdelta
                             pwisosn(jl,jt)=MAX(pwisosn(jl,jt)-zwisosnmel(jl,jt),0._dp)
                         END IF
                     END IF
                 END IF
             END DO
         END DO
     END IF
    
    END IF
!---wiso-code-end
    
    !*    2.5    Snow budget and meltwater (glacier-free land only)
    !
     DO 250 jl=1,kbdim
        IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
           pwl(jl)=pwl(jl)+zsncmelt(jl)                     ! Add melt water on canopy to skin reservoir
           zwl_tmp(jl)=pwl(jl)                              !---wiso-code: remember skin reservoir
           zsnmlt_tmp(jl)=MAX(0._dp,pwl(jl)-pwlmx(jl))      !---wiso-code: remember additional meltwater
           zsnmlt=zsnmel(jl)+MAX(0._dp,pwl(jl)-pwlmx(jl))   ! Excess water on canopy drips to ground and
           pwl(jl)=MIN(pwlmx(jl),pwl(jl))                   ! skin reservoir on canopy is set to maximum
           pmlres(jl)=zsnmlt
           zsnmel(jl)=zsnmel(jl)+zsncmelt(jl)
           zsn(jl)=zsn(jl)-zsnmel(jl)
        END IF
250  END DO

!---wiso-code
    IF (lwiso) THEN

    ! Snow budget and meltwater (glacier-free land only) - water isotopes
     DO jt=1,nwiso
         DO jl=1,kbdim
             IF (lmask(jl).AND. .NOT.lpglac(jl)) THEN
                  pwisowl(jl,jt)=pwisowl(jl,jt)+zwisosncmelt(jl,jt)
                 zdelta=tnat(jt)
                 IF (zwl_tmp(jl).gt.zwisomin .and. pwisowl(jl,jt).gt.zwisomin) zdelta=pwisowl(jl,jt)/zwl_tmp(jl)
                 IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp ! cut off rounding errors
                 zwisosnmlt=zwisosnmel(jl,jt) + (zsnmlt_tmp(jl) * zdelta)
                 pwisowl(jl,jt)=pwisowl(jl,jt) - (zsnmlt_tmp(jl) * zdelta)
                 pwisomlres(jl,jt)=zwisosnmlt
                 zwisosnmel(jl,jt)=zwisosnmel(jl,jt)+zwisosncmelt(jl,jt)
                 zwisosn(jl,jt)=zwisosn(jl,jt)-zwisosnmel(jl,jt)
             END IF
         END DO
     END DO
     
    END IF
!---wiso-code-end

    ! Accumulate water fluxes
    ! Note: conversion from m water equivalent to kg/m^2s - multiply by rhoh2o / zdtime
    !       accumulation - multiply by zdtime
      DO 601 jl=1,kbdim
         IF (lmask(jl)) THEN
            IF (.NOT. lpglac(jl)) psnacl(jl) = psnacl(jl)  +zsn(jl) * rhoh2o
            psnmel(jl) = psnmel(jl)  +zsnmel(jl) * rhoh2o
            IF (lpglac(jl)) progl(jl)  = progl(jl)   +zrogl(jl)  *rhoh2o
         END IF
 601  END DO

!---wiso-code
    IF (lwiso) THEN

    ! Accumulate water fluxes - Water Isotopes
    ! Note: conversion from m water equivalent to kg/m^2s - multiply by rhoh2o / zdtime
    !       accumulation - multiply by zdtime
      DO jt=1,nwiso
          DO jl=1,kbdim
              IF (lmask(jl)) THEN
                  IF (.NOT.lpglac(jl)) pwisosnacl(jl,jt) = pwisosnacl(jl,jt)  + zwisosn(jl,jt) * twisorhoh2o(jt)
                  pwisosnmel(jl,jt) = pwisosnmel(jl,jt)  + zwisosnmel(jl,jt) * twisorhoh2o(jt)
                  IF (lpglac(jl)) pwisorogl(jl,jt)  = pwisorogl(jl,jt)   + zwisorogl(jl,jt)  * twisorhoh2o(jt)
              END IF
          END DO
      END DO
      
    !Corrections of minor water pools
      DO jt = 1,nwiso
          DO jl = 1,kbdim
              IF (psn(jl).LT.zwisomin)  pwisosn(jl,jt)=psn(jl)*tnat(jt)
              IF (psnc(jl).LT.zwisomin) pwisosnc(jl,jt)=psnc(jl)*tnat(jt)
              IF (pgld(jl).LT.zwisomin) pwisogld(jl,jt)=pgld(jl)*tnat(jt)
          END DO
      END DO
    END IF
!---wiso-code-end

  END IF
  
  RETURN

END SUBROUTINE update_surf_down
