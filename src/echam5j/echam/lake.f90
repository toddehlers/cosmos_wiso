  SUBROUTINE lake ( kproma                                             &
                  , pseaice,    psiced,    palake                      &
                  , ptsi,       ptsw                                   &
                  , pahflw,     pahfsw,    pfluxres                    &
                  , ptrflw,     psoflw                                 &
                  , pevapi,     psni,      pcvsi                       &
                  , pahfres,    pfri                      )  

!
!  ---------------------------------------------------------------------
!
  USE mo_kind,           ONLY: dp
  USE mo_constants,      ONLY: tmelt, rhoh2o, alf
  USE mo_time_control,   ONLY: delta_time
!
  IMPLICIT NONE
!
  INTEGER, INTENT (IN):: kproma
!
! Arguments
!
  REAL(dp):: pseaice(kproma),   psiced(kproma),   palake(kproma)       &
            ,ptsi(kproma),      ptsw(kproma)                           &
            ,pahflw(kproma),    pahfsw(kproma),   pfluxres(kproma)     &
            ,ptrflw(kproma),    psoflw(kproma)                         &
            ,pevapi(kproma),    psni(kproma),     pcvsi(kproma)        &
            ,pahfres(kproma),   pfri(kproma)
  INTEGER :: jl
  REAL(dp):: zalpha, zalphas, zrho_sn, ziscond, zcpice, zrhoice        &
           , zdice, zrhoilf, zicecap, zdtrilf, zfreez, zdmix           &
           , zcpwater, zmixcap, zmcapdt, zmcaprilf, zfluxw, zts        &
           , zfres, zconhflx, zsubice, zhi
  REAL(dp):: zdtime

! Executable statements
!
! 1. Set up constants
!
  zdtime = delta_time
  zalpha=2.1656_dp
  zalphas=0.31_dp
  zrho_sn=330._dp
  ziscond=zalpha/zalphas*rhoh2o/zrho_sn
  zcpice=2106._dp
  zrhoice=910._dp
  zdice=0.10_dp
  zrhoilf=zrhoice*alf
  zicecap=zrhoice*zcpice*zdice
  zdtrilf=zdtime/zrhoilf
  zfreez=-zdice/zdtrilf
  zdmix=10._dp
  zcpwater=4218._dp
  zmixcap=rhoh2o*zcpwater*zdmix
  zmcapdt=zdtime/zmixcap
  zmcaprilf=zmixcap/zrhoilf
!
! 2. Lake temperature and ice thickness
!
  DO jl=1,kproma
!
  IF (palake(jl).GE.0.5_dp) THEN                           ! lake points
!
     IF (pseaice(jl).LT.0.5_dp) THEN                       ! open water
!
        zfluxw=pahflw(jl)+pahfsw(jl)+ptrflw(jl)+psoflw(jl)
!
!       Lake temperature (ptsw)
!
        zts=ptsw(jl)+zmcapdt*(zfluxw+pfluxres(jl))
        ptsi(jl)=tmelt
        pfluxres(jl)=0._dp
        psiced(jl)=0._dp
        IF (zts.GE.tmelt) THEN                  ! open water (unchanged)
           ptsw(jl)=zts
        ELSE                                    ! check ice formation
           ptsw(jl)=tmelt
           zfres=(zts-tmelt)/zmcapdt            ! < 0.
           IF (zfres.LE.zfreez) THEN            ! ice formation
              psiced(jl)=zmcaprilf*(tmelt-zts)  ! > zdice
              pseaice(jl)=1._dp
           ELSE
              pfluxres(jl)=zfres
           END IF
        END IF
!  ---------------------------------------------------------------------
     ELSE IF (psiced(jl).GE.zdice) THEN
!
!       Ice thickness (psiced)
!
        zconhflx=zalpha*(ptsi(jl)-tmelt)/(psiced(jl)+ziscond*psni(jl))
        zsubice=(1._dp-pcvsi(jl))*pevapi(jl)*zdtime/zrhoice
        zhi=psiced(jl)-zdtrilf*(zconhflx+pfluxres(jl))+zsubice
        ptsw(jl)=tmelt
        IF (zhi.GE.zdice) THEN
           psiced(jl)=zhi
           pseaice(jl)=1._dp
           pfluxres(jl)=0._dp
        ELSE IF (zhi.LE.0._dp) THEN               ! complete melting
           ptsw(jl)=tmelt-zhi/zmcaprilf        ! ptsw > tmelt
           psiced(jl)=0._dp
           pseaice(jl)=0._dp
           pfluxres(jl)=0._dp
        ELSE                                   ! incomplete melting
           psiced(jl)=zdice
           pseaice(jl)=1._dp
           pfluxres(jl)=(zdice-zhi)/zdtrilf
           pahfres(jl)=pahfres(jl)-zdtime*pfri(jl)*pfluxres(jl)
        END IF
     END IF
   END IF
  END DO
!  ---------------------------------------------------------------------
     RETURN
  END SUBROUTINE lake
