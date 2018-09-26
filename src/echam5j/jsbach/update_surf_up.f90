SUBROUTINE update_surf_up ( kbdim, itile          &
     , pwl,             pcvw,           pwlmx     &
     , pws,             pwsmx                     &
     , ptsoil,          pcvs                      &
     , porostd                                    &
     , pevapl,          pevapot,        pevwsd    &
     , prain                                      &
     , lmask,           lpglac                    &
     , prunoff,         pdrain                    &
     , pros_hd,         pdrain_hd                 &
     , pmlres                                     &
!---wiso-code
     , lwisofracl                                 &
! - 1D from mo_memory_g3b - water isotopes
     , pwisows,         pwisowl                   &
     , pwisorain			          &
     , pwisorunoff,     pwisodrain                &
! - variables internal to physics - water isotopes
     , pwisoevapot,     pwisoevwsd                &
     , pwisoros_hd,     pwisodrain_hd             &
     , pwisomlres                                 &
!---wiso-code-end
     )

!     ------------------------------------------------------------------
! pwl           Water content [m] in skin reservoir
!               (vegetation and bare soil)
! pwlmx         Skin reservoir [m] (calculated in *vdiff* as a function
!               of vegetation index and leaf area index)
! pcvw          Skin reservoir fraction (= pwl/pwlmx, see *vdiff*)
! pws           Soil water content [m]
! pwsmx         Water holding capacity [m] of the soil
! ptsoil        Temperature [K] in upper soil layer
! pcvs          Fractional snow cover (function of psn in *physc*)
! porostd       Subgrid standard diviation [m] used in runoff scheme
! pevapl        Total evapotranspiration, including sublimation [kg/m**2/s]
! pevapot       Potential evaporation/sublimation [kg/m**2/s]
! pevwsd        Evapotranspiration without sublimation and evaporation from interception reservoir [kg]/m**2]
! prain         Totalrainfall [kg/m**2/s]
! prunoff       Total runoff [kg/m**2] at non-glacier points (accumul.)
! pdrain        Drainage at non-glacier points [kg/m**2] (accumul.)
! pros_hd       Runoff for HD-Model (does NOT include drainage) [m]
! pdrain_hd     Drainage for HD-Model [m]
! lpglac        Logical glacier mask
!
!---wiso-code
!
! pwisows	tracer for Soil water content [m]
! pwisowl	tracer for Water content [m] in skin reservoir
!               (vegetation and bare soil)
! pwisorunoff	tracer for Total runoff [kg/m**2] at non-glacier points (accumul.)
! pwisorain     tracer for Totalrainfall [kg/m**2/s]
! pwisodrain	tracer for Drainage at non-glacier points [kg/m**2] (accumul.)
! pwisoevapot	tracer for Potential evaporation/sublimation [kg/m**2/s]
! pwisoevwsd    Evapotranspiration without sublimation and evaporation from interception reservoir [kg]/m**2]
! pwisoros_hd	tracer for Runoff for HD-Model (does NOT include drainage) [m]
! pwisodrain_hd	tracer for Drainage for HD-Model [m]
!
!---wiso-code-end
!
! The following local variables represent the respective fluxes
! integrated over one timestep (delta_time) and divided by the density of
! water (rhoh2o). Units are m water equivalent.
!
! zraind        Total rain
! zevwld        Evaporation from the skin reservoir
! zros          Total runoff (including drainage) at non-glacier points
! zdrain        Drainage at non-glacier points
! zsn           Snow budget at non-glacier points (snowfall-subl-melt)
! pmlres        Residual melt water available for infiltration into the
!               non-frozen soil after being intercepted by the
!               skin reservoir
!
!---wiso-code
!
! zwisoraind	tracer for Total rain
! zwisoevwld	tracer for Evaporation from the skin reservoir
! zwisoros      tracer for Total runoff (including drainage) at non-glacier points
! zwisodrain	tracer for Drainage at non-glacier points
! zwisomlres	tracer for Residual melt water available for infiltration into the
!
!---wiso-code-end
!
!       Rest folgt spaeter ....
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
  USE mo_control,                   ONLY: ngl
  USE mo_param_switches
  USE mo_physc2
  USE mo_constants
  USE mo_vegetation
  USE mo_radiation
  USE mo_time_control,              ONLY : lstart, delta_time
#ifdef STANDALONE
  USE mo_jsbach_comm_to_echam5mods, ONLY: nlat
#endif
  USE mo_kind,                      ONLY: dp

!---wiso-code
  USE mo_wiso,                      ONLY: lwiso, nwiso, twisorhoh2o, tnat
!---wiso-code-end

  IMPLICIT NONE
 
  INTEGER :: kbdim, jl, itile
  REAL(dp) ::                                                       &
       pws(kbdim),           pwl(kbdim),           pwsmx(kbdim)     &
     , prunoff(kbdim),       pdrain(kbdim)                          &
     , porostd(kbdim),       pmlres(kbdim)     
  REAL(dp) ::                                                       &
       ptsoil(kbdim) 
  REAL(dp) ::                                                       &
       pcvs(kbdim),          pcvw(kbdim),          pwlmx(kbdim)     &
     , pevapl(kbdim),        pevapot(kbdim), pevwsd(kbdim)          &
     , prain(kbdim)                                                 &
     , pros_hd(kbdim),       pdrain_hd(kbdim)
  LOGICAL ::                                                        &
       lmask(kbdim), lpglac(kbdim)
!
!  local arrays
!
  REAL(dp) ::                                                       &
       zraind (kbdim)                                               &
     , zevwld(kbdim)                                                &
     , zros(kbdim),          zdrain(kbdim)
!
!  local scalars
!
  REAL(dp) ::                                                         &
       zorvari, zorvars, zdrmin, zdrmax, zdrexp, zsmelt,              &
       zmprcp, zwlp, zwptr, zwdtr, zwslim, zconw2, zconw3, zconw4,    &
       zroeff, zbws, zb1, zbm, zconw1, zlyeps, zvol, zinfil, zprfl,   &
       zlysic, zwsup, zsncp, zexpt, zexpw, zsncmax, zsncwind
  REAL(dp) :: zdtime, zrcp, zsncfac

!---wiso-code
! declaration parameters water isotopes
  INTEGER  :: jt
  INTEGER, OPTIONAL  :: lwisofracl

  REAL(dp), OPTIONAL ::                                             &
       pwisows(kbdim,nwiso),        pwisowl(kbdim,nwiso)            &
       , pwisorunoff(kbdim,nwiso),  pwisodrain(kbdim,nwiso)         &
       , pwisomlres(kbdim,nwiso)                

  REAL(dp), OPTIONAL ::                                             &
         pwisoevapot(kbdim,nwiso),  pwisoevwsd(kbdim,nwiso)         &
       , pwisorain(kbdim,nwiso)                                     &
       , pwisoros_hd(kbdim,nwiso),  pwisodrain_hd(kbdim,nwiso)

! water isotopes: local arrays
      
  REAL(dp) ::                                                       &
       zwisoraind(kbdim,nwiso),     zwisoevwld(kbdim,nwiso)         &
       , zwisoros(kbdim,nwiso),     zwisodrain(kbdim,nwiso)         &
       , zws_tmp(kbdim),	    zb1_tmp(kbdim)                  &
       , zros_tmp(kbdim),           zwsup_tmp(kbdim)		    &
       , zprfl_tmp(kbdim),          zraind_tmp(kbdim)               &
       , zmprcp_tmp(kbdim),         zwlp_tmp(kbdim)                 &
       , zevwld_tmp(kbdim),         zevap_tmp(kbdim)             

  LOGICAL :: lo_evap(kbdim), lo_zros(kbdim)

!  water isotopes: local scalars

  REAL(dp) :: zprfltmp(kbdim)

  REAL(dp) ::                                                       &
       zwisomprcp,  zwisowlp,   zwisoinfil,                         &
       zwisoprfl,   zwisowsup,  zdelta,                             &
       zwisomin,    zwisosec,   zwisoevap
!---wiso-code-end

!      Parameters

  zdtime = delta_time
  zorvari=100._dp
#ifndef STANDALONE
  zorvars=1000._dp*64._dp/ngl
#else
  zorvars=1000._dp*64._dp/nlat
#endif
  zdrmin=0.001_dp/(3600._dp*1000._dp)
  zdrmax=0.1_dp/(3600._dp*1000._dp)
  zdrexp=1.5_dp
  zsncfac=rhoh2o*g/zdtime

!---wiso-code
    IF (lwiso) THEN

  zwisomin = 1.e-12_dp
  zwisosec = 1.e-10_dp
  
  END IF
!---wiso-code-end

!     -----------------------------------------------------------------

!
!*    1.     Convert water fluxes to [m water equivalent * timestep]
!
  zros = 0._dp
  zdrain = 0._dp
  zraind = prain * zdtime / rhoh2o
  zevwld = 0._dp
  DO 110 jl=1,kbdim
      IF (lmask(jl)) THEN
!!$        zros(jl)=0._dp
!!$        zdrain(jl)=0._dp
!!$        zraind(jl)=prain(jl) * zdtime/rhoh2o
          zevwld(jl)=(1._dp-pcvs(jl))*pcvw(jl)*pevapot(jl)*zdtime/rhoh2o
      END IF
 110 END DO

!---wiso-code
  IF (lwiso) THEN

  ! Convert water fluxes to [m water equivalent * timestep] - water isotopes
  DO jt=1,nwiso
      DO jl=1,kbdim
          zwisoros(jl,jt)=0._dp
          zwisodrain(jl,jt)=0._dp
          zwisoraind(jl,jt)= pwisorain(jl,jt)*zdtime/twisorhoh2o(jt)
          zwisoevwld(jl,jt)=0._dp
          IF (lmask(jl)) THEN
             !IF (lwisofracl==0) THEN
                ! calculate isotope ratio of skin reservoir
                zdelta=tnat(jt)
                IF (pwl(jl).gt.zwisomin .and. pwisowl(jl,jt).gt.zwisomin) zdelta = pwisowl(jl,jt)/pwl(jl)
                IF (abs(1._dp-zdelta).lt.zwisosec) zdelta=1._dp  ! cut off rounding errors
                zwisoevwld(jl,jt)=(1._dp-pcvs(jl))*pcvw(jl)*zdelta*pwisoevapot(jl,jt)*zdtime/twisorhoh2o(jt)
             !ELSE
             !   zwisoevwld(jl,jt)=(1._dp-pcvs(jl))*pcvw(jl)*pwisoevapot(jl,jt)*zdtime/twisorhoh2o(jt)
             !END IF
          END IF
      END DO
  END DO
  
  END IF
!---wiso-code-end

  IF(lsurf) THEN

!
!     ------------------------------------------------------------------
!
!*    4.     Water budget
!
!*    4.1    Skin reservoir (vegetation and bare soil)
!
      DO 410 jl=1,kbdim
         IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
!
!*    4.1.1  Interception of rain
!
            zraind_tmp(jl)=zraind(jl)                                          !---wiso-code: remember total rain
            zmprcp=MIN(zraind(jl)*cvinter,pwlmx(jl)-pwl(jl))
            zmprcp_tmp(jl)=zmprcp                                              !---wiso-code: remember interception
            zwlp=pwl(jl)+zmprcp
            zraind(jl)=zraind(jl)-zmprcp
!
!*    4.1.2  Evaporation or dew collection
!
            zwlp_tmp(jl)=zwlp                                                  !---wiso-code: remember skin reservoir

            zevwld_tmp(jl)=zevwld(jl)                                          !---wiso-code: remember atmospheric flux

            lo_evap(jl) = (zevwld(jl).GT.0.0_dp)                               !---wiso-code: flag for evaporation vs dew collection
            pwl(jl)=MIN(MAX(0._dp,zwlp+zevwld(jl)),pwlmx(jl))
            zevap_tmp(jl)=pwl(jl)-zwlp                                         !---wiso-code: remember evaporation/dew collection
            pevwsd(jl)=pevwsd(jl)-(pwl(jl)-zwlp)
         ELSE
            pwl(jl)=0._dp
         END IF
 410  END DO

!---wiso-code
    IF (lwiso) THEN

    ! Water budget - water isotopes
    ! Skin reservoir (vegetation and bare soil) - water isotopes
     DO jt=1,nwiso
         DO jl=1,kbdim
             IF (lmask(jl).AND..NOT.lpglac(jl)) THEN

                 ! Interception of rain - water isotopes
                 zdelta=tnat(jt)

                 IF (zraind_tmp(jl).gt.zwisomin .and. zwisoraind(jl,jt).gt.zwisomin) zdelta=zwisoraind(jl,jt)/zraind_tmp(jl)

                 IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp   ! cut off rounding errors

                 zwisomprcp=zmprcp_tmp(jl) * zdelta

                 zwisowlp=pwisowl(jl,jt)+zwisomprcp

                 zwisoraind(jl,jt)=zwisoraind(jl,jt)-zwisomprcp


                 ! Evaporation or dew collection - water isotopes
                 zdelta=tnat(jt)

                 IF (lo_evap(jl)) THEN                     ! dew collection, assume identical delta value as the atmospheric flux

                   IF (zevwld_tmp(jl).gt.zwisomin .and. zwisoevwld(jl,jt).gt.zwisomin) zdelta=zwisoevwld(jl,jt)/zevwld_tmp(jl)

                 ELSE                                      ! evaporation, assume identical delta value as the skin reservoir

                   IF (zwlp_tmp(jl).gt.zwisomin .and. zwisowlp.gt.zwisomin) zdelta=zwisowlp/zwlp_tmp(jl)

                 ENDIF  

                 IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp   ! cut off rounding errors

                 zwisoevap=zevap_tmp(jl) * zdelta

                 pwisowl(jl,jt)=zwisowlp+zwisoevap

                 pwisoevwsd(jl,jt)=pwisoevwsd(jl,jt)-zwisoevap
             ELSE
                 pwisowl(jl,jt)=0._dp
             END IF
         END DO
     END DO
     
     END IF
!---wiso-code-end

!
!*    4.2    Soil reservoir
!
      DO 420 jl=1,kbdim
          IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
             zwptr=0.90_dp*pwsmx(jl)
             zwdtr=0.90_dp*pwsmx(jl)
             zwslim=0.05_dp*pwsmx(jl)
             zconw2=pwsmx(jl)-zwdtr
             zconw3=zdrmax-zdrmin
             zconw4=pwsmx(jl)-zwptr
             zroeff=MAX(0._dp, porostd(jl)-zorvari)   &
                          /(porostd(jl)+zorvars)
             zbws=MAX(MIN(zroeff,0.5_dp),0.01_dp)
             zb1=1._dp+zbws
             zb1_tmp(jl)=zb1                                         !---wiso-code: remember value
             zbm=1._dp/zb1
             zconw1=pwsmx(jl)*zb1
             zlyeps=0._dp
             zvol=0._dp
             zinfil=0._dp
!
!*    4.2.1  Surface runoff, infiltration and evaporation from soil
!
             IF (pevwsd(jl) >= 0.0_dp) THEN
               zprfl=pmlres(jl)+zraind(jl)+pevwsd(jl)
             ELSE
               pws(jl)=pws(jl)+pevwsd(jl)
               zprfl=pmlres(jl)+zraind(jl)
             END IF
             zprfl_tmp(jl)=zprfl                                     !---wiso-code: remember total precipitation value
             IF (ptsoil(jl).LT.tmelt) THEN
                zros(jl)=zprfl
             ELSE
                lo_zros(jl)= (zprfl.GT.0._dp.AND.pws(jl).GT.zwslim)  !---wiso-code: remember if runoff occurs
                IF (zprfl.GT.0._dp.AND.pws(jl).GT.zwslim) THEN
                   IF (pws(jl).GT.pwsmx(jl)) THEN
                      zlyeps=pws(jl)-pwsmx(jl)
                   ELSE
                      zlyeps=0._dp
                   END IF
                   zlysic=(pws(jl)-zlyeps)/pwsmx(jl)
                   zlysic=MIN(zlysic,1._dp)
                   zvol=(1._dp-zlysic)**zbm-zprfl/zconw1
                   zros(jl)=zprfl-(pwsmx(jl)-pws(jl))
                   IF (zvol.GT.0._dp) THEN
                      zros(jl)=zros(jl)+pwsmx(jl)*zvol**zb1
                   END IF
                   zros(jl)=MAX(zros(jl),0._dp)
                   zros_tmp(jl)=zros(jl)                             !---wiso-code: remember runoff value
                   zinfil=zprfl-zros(jl)
                ELSE
                   zros(jl)=0._dp
                   zinfil=zprfl
                END IF
                pws(jl)=pws(jl)+zinfil
             END IF

!*    4.2.2  Drainage and total runoff
             zws_tmp(jl)=pws(jl)                                     !---wiso-code: remember value of pws before drainage
             IF (pws(jl).LE.zwslim) THEN
                zdrain(jl)=0._dp
             ELSE
                IF (ptsoil(jl).GT.tmelt) THEN
                   zdrain(jl)=zdrmin*pws(jl)/pwsmx(jl)
                   IF (pws(jl).GT.zwdtr) THEN
                      zdrain(jl)=zdrain(jl)+zconw3*                  &
                                ((pws(jl)-zwdtr)/zconw2)**zdrexp
                   END IF
                   zdrain(jl)=zdrain(jl)*zdtime
                   zdrain(jl)=MIN(zdrain(jl),pws(jl)-zwslim)
                   pws(jl)=pws(jl)-zdrain(jl)
                ELSE
                   zdrain(jl)=0._dp
                END IF
             END IF
             zwsup=MAX(pws(jl)-pwsmx(jl),0._dp)
             zwsup_tmp(jl)=zwsup                                     !---wiso-code: remember soil water surplus
             pws(jl)=pws(jl)-zwsup
             zros(jl)=zros(jl)+zdrain(jl)+zwsup
          ELSE
             pws(jl)=0._dp
          END IF
 420   END DO
 
!---wiso-code
       IF (lwiso) THEN

       ! Soil reservoir - water isotopes
       DO jt=1,nwiso
           DO jl=1,kbdim
               IF (lmask(jl).AND..NOT.lpglac(jl)) THEN
                   zwslim=0.05_dp*pwsmx(jl)
                   zwisoinfil=0._dp
                   
                   ! Surface runoff, infiltration and evaporation from soil - water isotopes
                   IF (pevwsd(jl) >= 0.0_dp) THEN
                       zwisoprfl=pwisomlres(jl,jt)+zwisoraind(jl,jt)+pwisoevwsd(jl,jt)
                   ELSE
                       pwisows(jl,jt)=pwisows(jl,jt)+pwisoevwsd(jl,jt)
                       zwisoprfl=pwisomlres(jl,jt)+zwisoraind(jl,jt)
                   END IF
      
                   ! assume that surface runoff has same delta value as precipitation (no mixing with soil water pool ws)
                   zdelta=tnat(jt)
                   IF (zprfl_tmp(jl).gt.zwisomin .and. zwisoprfl.gt.zwisomin) zdelta=zwisoprfl/zprfl_tmp(jl)
                   IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp ! cut off rounding errors
                   IF (ptsoil(jl).LT.tmelt) THEN
                       zwisoros(jl,jt)=zwisoprfl
                   ELSE 
                       IF (lo_zros(jl)) THEN
                           zwisoros(jl,jt)=zros_tmp(jl)*zdelta
                           zwisoinfil=zwisoprfl-zwisoros(jl,jt)
                       ELSE
                           zwisoros(jl,jt)=0._dp
                           zwisoinfil=zwisoprfl
                       END IF
                       pwisows(jl,jt)=pwisows(jl,jt)+zwisoinfil
                   END IF
      
                   ! Drainage and total runoff - water isotopes
                   ! drainage water and soil water surplus have same delta value as soil water ws
                   zdelta=tnat(jt)
                   IF (zws_tmp(jl).gt.zwisomin .and. pwisows(jl,jt).gt.zwisomin) zdelta=pwisows(jl,jt)/zws_tmp(jl)
                   IF (ABS(1._dp-zdelta).lt.zwisosec) zdelta=1._dp ! cut off rounding errors
      
                   IF (zws_tmp(jl).LE.zwslim) THEN
                       zwisodrain(jl,jt)=0._dp
                   ELSE
                       IF (ptsoil(jl).GT.tmelt) THEN
                           zwisodrain(jl,jt)=zdrain(jl)*zdelta
                           pwisows(jl,jt)=pwisows(jl,jt)-zwisodrain(jl,jt)
                       ELSE
                           zwisodrain(jl,jt)=0._dp
                       END IF
                   END IF
                   zwisowsup=zwsup_tmp(jl)*zdelta
                   pwisows(jl,jt)=pwisows(jl,jt)-zwisowsup
                   zwisoros(jl,jt)=zwisoros(jl,jt)+zwisodrain(jl,jt)+zwisowsup
               ELSE
                   pwisows(jl,jt)=0._dp
               END IF
           END DO
       END DO
       
       END IF
!---wiso-code-end

!*     4.2.3  Runoff and drainage for the HD-Model
!
      DO 423 jl=1,kbdim
         IF (lmask(jl)) THEN
            pros_hd(jl)=zros(jl)-zdrain(jl)
            pdrain_hd(jl)=zdrain(jl)
         END IF
 423  END DO

!---wiso-code
    IF (lwiso) THEN

    ! Runoff and drainage for the HD-Model - water isotopes
    DO jt=1,nwiso
        DO jl=1,kbdim
            IF (lmask(jl)) THEN
                pwisoros_hd(jl,jt)=zwisoros(jl,jt)-zwisodrain(jl,jt)
                pwisodrain_hd(jl,jt)=zwisodrain(jl,jt)
            END IF
        END DO
    END DO
    
    END IF
!---wiso-code-end
!
!     ------------------------------------------------------------------
!
!*    6.     Accumulate fluxes for diagnostics
!
!     6.1    Water fluxes
!
      ! Note: conversion from m water equivalent to kg/m^2s - multiply by rhoh2o / zdtime
      !       accumulation - multiply by zdtime
      DO 601 jl=1,kbdim
         IF (lmask(jl)) THEN
            prunoff(jl)= prunoff(jl) +zros(jl)   *rhoh2o
            pdrain(jl) = pdrain(jl)  +zdrain(jl) *rhoh2o
         END IF
 601  END DO

!---wiso-code
    IF (lwiso) THEN

!   Water fluxes - water isotopes
    DO jt=1,nwiso
        DO jl=1,kbdim
            IF (lmask(jl)) THEN
                pwisorunoff(jl,jt)= pwisorunoff(jl,jt) +zwisoros(jl,jt)   *twisorhoh2o(jt)
                pwisodrain(jl,jt) = pwisodrain(jl,jt)  +zwisodrain(jl,jt) *twisorhoh2o(jt)
            END IF
        END DO
    END DO
!    Corrections of minor water pools
     DO jt = 1,nwiso
         DO jl = 1,kbdim
             IF (pwl(jl).LT.zwisomin) pwisowl(jl,jt)=pwl(jl)*tnat(jt)
             IF (pws(jl).LT.zwisomin) pwisows(jl,jt)=pws(jl)*tnat(jt)
         END DO
     END DO
     
     END IF
!---wiso-code-end

!     ------------------------------------------------------------------
  END IF
!
  RETURN

END SUBROUTINE update_surf_up
