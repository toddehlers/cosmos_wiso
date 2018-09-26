MODULE mo_ocice  
 !     hiblers ice model as described in
 !     w.d.hibler iii, a dynamic thermodynamic sea ice model. j. phys. oceanogr.
 !     9, 815 - 846, 1979

 !
 !     modifications:
 !     uwe 2.2.00
 !       include advection of ice velocities
 !       correct tau-w term using mixture of old and new velocities
 !        new fileds included:
 !         speedu    difference in speed between water and ice (old)
 !         speedv
 !     uwe 8.3.00
 !        include tauwat, proper determination of water ice stress
 !        to be used in ocwind
 !     uwe 28.6.00
 !        make sw-penetration applicable under sea ice too
 !     uwe 10/00
 !        include implicite treatment of mass advection in iteration 
 !         of ice velocities
  USE mo_param1
  USE mo_parallel
  USE mo_commo1
  USE mo_commoau1
  USE mo_commoau2
  USE mo_units

#ifdef PBGC
  USE mo_param1_bgc , only: nocetra
  USE mo_carbch , only: ocetra
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef ADDCONTRA
  USE mo_param1_add , only: nocectra
  USE mo_contra , only: ocectra
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifdef __coupled
  USE mo_fluxes1
#endif

IMPLICIT NONE




CONTAINS

SUBROUTINE ocice



INTEGER i,j
REAL sicomo1
REAL siotho(ie,je),sioomo(ie,je),siosno(ie,je)
  !     taf atmospheric temperature
  !     hth00: arbitrary constant for discrimination between thin and thick ice
  !          set to 1.5 m   not in use
  !     sicth  ice thickness
  !     hicce,hiccp the coefficients e and p from (7), (8)
  !     hiccp is divided by mean density of 1000 kg / m**3
  !     sicom  ice compactness
  !     sicdi  ice flow divergence
  !     sicsh  ice flow shear
  !     hibzet,hibdel,hibet are the fields of (7) - (9)
  !     entmel:melting enthalpy of ice = 320 millions ws/m**3
  !      stebol: stefan boltzmann constant =5.67 10**(-8)w/m**2/k**4
  !     sichec: heat conductivity of ice = 2 w/m/k   cuwe not in use
  !
  !      part 4: increase of existing ice
  ! in theory a distinction should be made bewteen the ice-covered part of a 
  ! grid-cell for which this do-loop applies and the ice free part, which
  ! should be treated as do loops 20 and 25. such a distinction enhances the
  ! overall ice-growth which results in irrealistically high values in some
  ! grid points and subsequent numerical instabilities. inclusion of diffusion
  ! on icethickness and compactness also gave rise to instabilities. 
  !
  !      cutoff value for prognostic calculation of ice velocities
  !
!$OMP PARALLEL  private(i,j,sicomo1)
!$OMP DO
  DO j=2,je1
    DO i=2,ie1
      sicomo1 = sicomo(i,j)
      sictho(i,j)=MAX(0.,sictho(i,j)*weto(i,j,1))
      IF(sictho(i,j) <= 0.) sicomo1=0.
      sicomo1=MAX(0.,sicomo1*weto(i,j,1))
      sicsno(i,j)=MAX(0.,sicsno(i,j)*weto(i,j,1))
      sicomo(i,j)=MIN(1.,sicomo1)
    END DO
  END DO
!$OMP END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('p',sicomo,'mo_ocice 100')
  CALL bounds_exch('p',sicsno,'mo_ocice 101')
  CALL bounds_exch('p',sictho,'mo_ocice 102')
!$OMP END SINGLE
!#endif

! store old values for sr growth

!$OMP DO
  DO j=1,je
    DO i=1,ie
      siotho(i,j)=sictho(i,j)
      sioomo(i,j)=sicomo(i,j)
      siosno(i,j)=sicsno(i,j)
    END DO
  END DO
!$OMP END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('u',sicuo,'mo_ocice 103')
  CALL bounds_exch('v',sicve,'mo_ocice 104')
!$OMP END SINGLE
!#endif

!$OMP END PARALLEL

CALL ice_dynamics
CALL ice_advection
CALL ice_thermodynamics(siotho, sioomo, siosno)
  
END SUBROUTINE ocice

SUBROUTINE ice_dynamics

  INTEGER i,j,iter

  REAL rhoicsn,rhowaic,uepsi
  REAL eps11(ie,je),eps22(ie,je),eps12(ie,je)
  REAL uh(ie,je),vh(ie,je)
  REAL ux,uy,vx,vy
  REAL e12,e11,e22,hibcc,pst,rst,sicm,sicom,dxdx,dydy,dxdy
  REAL alpalt,alpneu
  REAL siouo(ie,je),siove(ie,je)
  REAL effico(ie,je),effice(ie,je)                  
  REAL cweffo(ie,je),cweffe(ie,je)                
  REAL zsic(ie,je)
  REAL zschaltd
  REAL zschalt
  REAL ucor(ie,je),vcor(ie,je)
  REAL zdiff(ie,je)
  REAL zzl(ie,je),zzr(ie,je),zzo(ie,je),zzu(ie,je)
  REAL rs,rv,speedmi,reval
  REAL speedu(ie,je),speedv(ie,je)  

  reval=0.01
  speedmi=0.01
  rhowaic=rhowat/rhoice
  rhoicsn=rhoice/rhosno
  uepsi=1.e-8

!$OMP PARALLEL  private(i,j,ux,uy,vx,vy,e12,e11,e22,hibcc, &
!$OMP   pst,rst,sicm,sicom,dxdx,dydy,dxdy,rs,rv)

! DN: Commenting out this block gives binary different results. Why??

!$OMP DO
  DO j=1,je
     DO i=1,ie
        sicomp(i,j)=sicomo(i,j)
    END DO
  END DO
!$OMP END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('p',sicomp,'mo_ocice 105')
!$OMP END SINGLE
!#endif

!$OMP DO
      DO J=1,JE
       DO I=1,IE
         SICUDO(I,J)=SICUO(I,J)
         SICVDE(I,J)=SICVE(I,J)
       ENDDO
      ENDDO
!$OMP END DO



!      ice dynamics

!$OMP DO
  DO j=1,je
    DO i=1,ie
      eps11(i,j)=0.
      eps12(i,j)=0.
      eps22(i,j)=0.
      uh(i,j)=0.
      vh(i,j)=0.
      zdiff(i,j)=0.
      effico(i,j)=0.
      effice(i,j)=0.
      speedu(i,j)=0.
      speedv(i,j)=0.
      zzl(i,j)=0.
      zzr(i,j)=0.
      zzo(i,j)=0.
      zzu(i,j)=0.
      ucor(i,j)=0.
      vcor(i,j)=0.
      siouo(i,j)=0.
      siove(i,j)=0.
   END DO
  END DO
!$OMP END DO
  
!$OMP DO
  DO j=2,je1
    DO i=2,ie1
      ux=(sicuo(i,j)-sicuo(i-1,j))/dlxp(i,j)
      vy=(sicve(i,j-1)-sicve(i,j))/dlyp(i,j)
      vx=(sicve(i+1,j)-sicve(i,j))/dlxp(i,j)
      uy=(sicuo(i,j)-sicuo(i,j+1))/dlyp(i,j)
      eps11(i,j)=ux
      eps22(i,j)=vy
      eps12(i,j)=0.5*(vx+uy)
      IF(amsuo(i,j,1).GT.0.5)                                          &
           speedu(i,j)=MAX(SQRT((uepsi+sicuo(i,j)-uko(i,j,1))**2           &
           +(0.25*(sicve(i,j)+sicve(i+1,j)+sicve(i,j-1)+sicve(i+1,j-1)    &
           -vke(i,j,1)-vke(i+1,j,1)-vke(i,j-1,1)-vke(i+1,j-1,1)))**2)   &
           ,speedmi)
      IF(amsue(i,j,1).GT.0.5)                                          &
           speedv(i,j)=MAX(SQRT((uepsi+sicve(i,j)-vke(i,j,1))**2           &
           +(0.25*(sicuo(i,j)+sicuo(i-1,j)+sicuo(i,j+1)+sicuo(i-1,j+1)    &
           -uko(i,j,1)-uko(i-1,j,1)-uko(i,j+1,1)-uko(i-1,j+1,1)))**2)   &
           ,speedmi)
    END DO
  END DO
!$OMP END DO


!$OMP SINGLE
!#ifdef bounds_exch_save
  CALL bounds_exch('p',eps11,'mo_ocice 106')
  CALL bounds_exch('p',eps22,'mo_ocice 107')
  CALL bounds_exch('s',eps12,'mo_ocice 108')
!#endif
  CALL bounds_exch('u',speedu,'mo_ocice 109')
  CALL bounds_exch('v',speedv,'mo_ocice 110')
!$OMP END SINGLE
  
!$OMP DO
  DO j=2,je1
    DO i=2,ie1
      e12=0.25*(eps12(i-1,j-1)+eps12(i,j-1)+eps12(i-1,j)+eps12(i,j))
      !      argu=((eps11(i,j)-eps22(i,j))**2+4.*e12**2)/hicce**2             &
      !    & +(eps11(i,j)+eps22(i,j))**2
      hibdelo(i,j)=SQRT(((eps11(i,j)-eps22(i,j))**2+4.*e12**2)/hicce**2       &
           +(eps11(i,j)+eps22(i,j))**2)
    END DO
  END DO
!$OMP END DO
 
!$OMP DO
  DO j=2,je1
    DO i=2,ie1
      e11=0.25*(eps11(i,j)+eps11(i+1,j)+eps11(i,j+1)+eps11(i+1,j+1))
      e22=0.25*(eps22(i,j)+eps22(i+1,j)+eps22(i,j+1)+eps22(i+1,j+1))
      hibdele(i,j)=SQRT(  ((e11-e22)**2+4.*eps12(i,j)**2)/hicce**2     &
           +(e11+e22)**2)
    END DO
  END DO
!$OMP END DO
!#ifdef bounds_exch_save
!$OMP SINGLE
   CALL bounds_exch('p',hibdelo,'mo_ocice 111')
   CALL bounds_exch('s',hibdele,'mo_ocice 112a')
!$OMP END SINGLE
!#endif
!uwe introduce hibcc, corresponds to c IN hibler, jpo 9,825

  hibcc=20.
 
  alpneu=1.
  
  alpalt=1.-alpneu

!$OMP DO
  DO j=2,je1
    DO i=2,ie1
      
!uwe    hibler (17), 
      pst=hiccp*sictho(i,j)*EXP(-hibcc*(1.-sicomo(i,j)))*weto(i,j,1)
      sicm=(sictho(i,j)+sictho(i+1,j)+sictho(i,j+1)+sictho(i+1,j+1))            &
           /(weto(i,j,1)+weto(i+1,j,1)+weto(i,j+1,1)+weto(i+1,j+1,1)+1.e-20)
      sicom=(sicomo(i,j)+sicomo(i+1,j)+sicomo(i,j+1)+sicomo(i+1,j+1))           &
           /(weto(i,j,1)+weto(i+1,j,1)+weto(i,j+1,1)+weto(i+1,j+1,1)+1.e-20)
      sicom=MIN(sicom,1.)
      sicom=MAX(sicom,0.)
      sicm=MAX(sicm,0.)
      rst=hiccp*sicm*EXP(-hibcc*(1.-sicom))
      uh(i,j)=sicuo(i,j)
      vh(i,j)=sicve(i,j)
      siouo(i,j)=sicuo(i,j)
      siove(i,j)=sicve(i,j)
      IF(hibdelo(i,j).GT.almzer)THEN
        hibzeto(i,j)=alpneu*pst*0.5/MAX(hibdelo(i,j),almzer)                    &
             +alpalt*hibzeto(i,j)
        hibeto(i,j)=alpneu*hibzeto(i,j)/hicce**2+alpalt*hibeto(i,j)
      ELSE
        hibzeto(i,j)=0.
        hibeto(i,j)=0.
      ENDIF
      IF(hibdele(i,j).GT.almzer) THEN
        hibzete(i,j)=alpneu*rst*0.5/MAX(hibdele(i,j),almzer)           &
             +alpalt*hibzete(i,j)
        hibete(i,j)=alpneu*hibzete(i,j)/hicce**2+alpalt*hibete(i,j)
      ELSE
        hibzete(i,j)=0.
        hibete(i,j)=0.
      ENDIF
      !c       spur=eps11(i,j)+eps22(i,j)
      !
      effico(i,j)=0.5*(sictho(i,j)+sictho(i+1,j)                       &
           +rhoicsn*(sicsno(i,j)+sicsno(i+1,j)))*amsuo(i,j,1)
      effice(i,j)=0.5*(sictho(i,j)+sictho(i,j+1)                       &
           +rhoicsn*(sicsno(i,j)+sicsno(i,j+1)))*amsue(i,j,1)
    ENDDO
  ENDDO
!$OMP END DO

!$OMP SINGLE
!#ifdef bounds_exch_save
  CALL bounds_exch('u',uh,'mo_ocice 112')
  CALL bounds_exch('v',vh,'mo_ocice 113')
  CALL bounds_exch('u',siouo,'mo_ocice 114')
  CALL bounds_exch('v',siove,'mo_ocice 115')
  CALL bounds_exch('p',hibzeto,'mo_ocice 116')
  CALL bounds_exch('p',hibeto,'mo_ocice 117')
  CALL bounds_exch('s',hibzete,'mo_ocice 118')
  CALL bounds_exch('s',hibete,'mo_ocice 119')
  CALL bounds_exch('u',effico,'mo_ocice 120')
  CALL bounds_exch('v',effice,'mo_ocice 121')
!#endif
!$OMP END SINGLE


!      enhance friction coefficient for very thin ice
!      ==> smoothe transition towards water velocities for effic ==> 0.

!$OMP DO
  DO j=1,je
    DO i=1,ie
      !c         cweffo(i,j)=cw*MAX(1.,1./(25.*(effico(i,j)+1.e-6)**2))
      !c         cweffe(i,j)=cw*MAX(1.,1./(25.*(effice(i,j)+1.e-6)**2))
      cweffo(i,j)=cw*MAX(1.,1./(10.*effico(i,j)+1.e-6))**1
      cweffe(i,j)=cw*MAX(1.,1./(10.*effice(i,j)+1.e-6))**1
      sicudo(i,j)=(0.8*sicudo(i,j)+0.2*sicuo(i,j))*amsuo(i,j,1)
      sicvde(i,j)=(0.8*sicvde(i,j)+0.2*sicve(i,j))*amsue(i,j,1)
      zsic(i,j)=sictho(i,j)+rhosnic*sicsno(i,j)
    END DO
  END DO
!$OMP END DO
  
  !      old ice velocities, u on v-point and v on u-point
  !
  !      zschaltd: upwind type diffusion
 
  zschaltd=1.
!$OMP DO
  DO j=2,je1
    DO i=2,ie1
      ucor(i,j)=0.25*(uh(i,j)+uh(i-1,j)+uh(i,j+1)+uh(i-1,j+1))
      vcor(i,j)=0.25*(vh(i,j)+vh(i+1,j)+vh(i,j-1)+vh(i+1,j-1))
      !
      !     sea level change due to diffusion of sea ice
      !
      zdiff(i,j)=zschaltd*0.5*dt*rhoicwa*                              &
           (ABS(sicudo(i-1,j))*dlyu(i-1,j)*(zsic(i-1,j)-zsic(i,j))       &
           +ABS(sicudo(i,j))*dlyu(i,j)*(zsic(i+1,j)-zsic(i,j))           &
           +ABS(sicvde(i,j))*dlxv(i,j)*(zsic(i,j+1)-zsic(i,j))           &
           +ABS(sicvde(i,j-1))*dlxv(i,j-1)*(zsic(i,j-1)-zsic(i,j)))      &
           /(dlxp(i,j)*dlyp(i,j))
      !
      !     thickness protection, avoid negative layer thickness
      !       induce divergent ice transports
      !
      !c     x  +  MAX(0.,zsic(i,j)-0.6667*(dzw(1)+zo(i,j)))
      !c     x     /MAX(0.01,dzw(1)+zo(i,j)-zsic(i,j))
      !
    END DO
  END DO
!$OMP END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('u',ucor,'mo_ocice 122')
  CALL bounds_exch('v',vcor,'mo_ocice 123')
  CALL bounds_exch('p',zdiff,'mo_ocice 124')
!$OMP END SINGLE
!#endif

  zschalt=1.
  !      main iteration

  !      iteration_loop: 

DO iter = 1, 40
   
!$OMP DO
    DO j=1,je
      DO i=1,ie
        uh(i,j)=0.
        vh(i,j)=0.
      END DO
    END DO
!$OMP END DO
    
!$OMP DO
    DO j=2,je1
      DO i=2,ie1
        zzl(i,j)=dt*rhoicwa*(                                            &
             dlyu(i-1,j)*sicuo(i-1,j)*effico(i-1,j)                       &
             +dlxv(i,j)*sicve(i,j)*effice(i,j)                             &
             -dlxv(i,j-1)*sicve(i,j-1)*effice(i,j-1))                      &
             /(dlxp(i,j)*dlyp(i,j))
        zzr(i,j)=dt*rhoicwa*(                                            &
             -dlyu(i+1,j)*sicuo(i+1,j)*effico(i+1,j)                       &
             +dlxv(i+1,j)*sicve(i+1,j)*effice(i+1,j)                       &
             -dlxv(i+1,j-1)*sicve(i+1,j-1)*effice(i+1,j-1))                &
             /(dlxp(i+1,j)*dlyp(i+1,j))
        speedu(i,j)=0.5*(MAX(SQRT((uepsi+sicuo(i,j)-uko(i,j,1))**2       &
             +(0.25*(sicve(i,j)+sicve(i+1,j)+sicve(i,j-1)+sicve(i+1,j-1)    &
             -vke(i,j,1)-vke(i+1,j,1)-vke(i,j-1,1)-vke(i+1,j-1,1)))**2)   &
             ,speedmi)+speedu(i,j))
        speedv(i,j)=0.5*(MAX(SQRT((uepsi+sicve(i,j)-vke(i,j,1))**2       &
             +(0.25*(sicuo(i,j)+sicuo(i-1,j)+sicuo(i,j+1)+sicuo(i-1,j+1)    &
             -uko(i,j,1)-uko(i-1,j,1)-uko(i,j+1,1)-uko(i-1,j+1,1)))**2)   &
             ,speedmi)+speedv(i,j))
      END DO
    END DO
!$OMP END DO

!$OMP SINGLE
!#ifdef bounds_exch_save
    CALL bounds_exch('p',zzl,'mo_ocice 124a')
    CALL bounds_exch('p',zzr,'mo_ocice 125')
!#endif
    CALL bounds_exch('u',speedu,'mo_ocice 126')
    CALL bounds_exch('v',speedv,'mo_ocice 127')
!$OMP END SINGLE

!$OMP DO
    DO j=2,je1
      DO i=2,ie1
         
        IF (effico(i,j).GT.reval)THEN
          dxdx=dlxp(i,j)**2
          dydy=dlyp(i,j)**2
          dxdy=dlxp(i,j)*dlyp(i,j)
          rs=effico(i,j)*(siouo(i,j)+0.5*dt*ftwou(i,j)*(vcor(i,j)          &
               +0.25*(sicve(i,j)+sicve(i+1,j)+sicve(i,j-1)+sicve(i+1,j-1))))     &
               +cweffo(i,j)*dt*uko(i,j,1)*speedu(i,j)                           &
               + dt*hiccp*(sictho(i,j)*EXP(-hibcc*(1.-sicomo(i,j)))              &
               -sictho(i+1,j)*EXP(-hibcc*(1.-sicomo(i+1,j))))/dlxu(i,j)       &
               + dt*effico(i,j)*g*(zo(i,j)-zo(i+1,j)+zdiff(i,j)-zdiff(i+1,j)     &
               +zschalt*(zzl(i,j)-zzr(i,j)))/dlxu(i,j)                       &
#ifdef __coupled
               !sv 31.08.99 included fluxes option
               !   (i.e. distinguish bewteen tau over water and ice)
               !sv here only tau over ice is considered (i.e. txo is replaced by aofltxio)
               !sv tau over water is used IN sbr ocwind
               +dt*rhowaic*aofltxio(i,j)*(sicomo(i,j)+sicomo(i+1,j))*0.5
#else
               +dt*rhowaic*txo(i,j)*(sicomo(i,j)+sicomo(i+1,j))*0.5
#endif
               !svx 31.08.99

          rv=                                                              &
               +(hibzeto(i,j)+hibeto(i,j))    *sicuo(i-1,j)/dxdx                 &
               +(hibzeto(i+1,j)+hibeto(i+1,j))*sicuo(i+1,j)/dxdx                 &
               +(hibzeto(i+1,j)-hibeto(i+1,j))*(sicve(i+1,j-1)-sicve(i+1,j))     &
               /dxdy                                                            &
               -(hibzeto(i,j)-hibeto(i,j))*(sicve(i,j-1)-sicve(i,j))/dxdy        &
               +(hibete(i,j-1)*sicuo(i,j-1)+hibete(i,j)*sicuo(i,j+1))/dydy       &
               +hibete(i,j-1)*(sicve(i+1,j-1)-sicve(i,j-1))/dxdy                 &
               -hibete(i,j)*(sicve(i+1,j)-sicve(i,j))/dxdy

          uh(i,j)=(rs+rv*dt)/(effico(i,j)                                  &
               +cweffo(i,j)*dt*speedu(i,j)                                 &
               +dt*effico(i,j)*g*dt*rhoicwa*effico(i,j)*zschalt*              &
               (1./(dlxp(i,j)*dlyp(i,j))+1./(dlxp(i+1,j)*dlyp(i+1,j)))         &
               *dlyu(i,j)/dlxu(i,j)                                              &
               +dt*(hibzeto(i,j)+hibeto(i,j)+hibzeto(i+1,j)+hibeto(i+1,j))/dxdx &
               +dt*(hibete(i,j-1)+hibete(i,j))/dydy)
        ENDIF

      END DO
    END DO
!$OMP END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch('u',uh,'mo_ocice 128')
!$OMP END SINGLE
!#endif

!$OMP DO
    DO j=2,je1
      DO i=2,ie1
        zzo(i,j)=dt*rhoicwa*(                                           &
             dlyu(i-1,j)*sicuo(i-1,j)*effico(i-1,j)                      &
             -dlyu(i,j)*sicuo(i,j)*effico(i,j)                            &
             -dlxv(i,j-1)*sicve(i,j-1)*effice(i,j-1))                     &
             /(dlxp(i,j)*dlyp(i,j))
        zzu(i,j)=dt*rhoicwa*(                                           &
             dlyu(i-1,j+1)*sicuo(i-1,j+1)*effico(i-1,j+1)                &
             -dlyu(i,j+1)*sicuo(i,j+1)*effico(i,j+1)                      &
             +dlxv(i,j+1)*sicve(i,j+1)*effice(i,j+1))                     &
             /(dlxp(i,j+1)*dlyp(i,j+1))
      END DO
    END DO
!$OMP END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch('p',zzo,'mo_ocice 129')
    CALL bounds_exch('p',zzu,'mo_ocice 130')
!$OMP END SINGLE
!#endif

!$OMP DO
    DO j=2,je1      
      DO i=2,ie1
        IF(effice(i,j).GT.reval) THEN
          dxdx=dlxp(i,j)**2
          dydy=dlyp(i,j)**2
          dxdy=dlxp(i,j)*dlyp(i,j)
          rs=effice(i,j)*(siove(i,j)-0.5*dt*ftwov(i,j)*(ucor(i,j)           &
               +0.25*(sicuo(i,j)+sicuo(i-1,j)+sicuo(i,j+1)+sicuo(i-1,j+1))))     &
               +cweffe(i,j)*dt*vke(i,j,1)*speedv(i,j)                           &
               +dt*hiccp*(sictho(i,j+1)*EXP(-hibcc*(1.-sicomo(i,j+1)))          &
               -sictho(i,j)*EXP(-hibcc*(1.-sicomo(i,j))))/dlyv(i,j)            &
               +effice(i,j)*g*dt*(zo(i,j+1)-zo(i,j)+zdiff(i,j+1)-zdiff(i,j)     &
               +zschalt*(zzu(i,j)-zzo(i,j)))/dlyv(i,j)                      &
#ifdef __coupled
               !sv 31.08.99 included fluxes option
               !    (i.e. distinguish between tau over water and ice)
               !sv  here only tau over ice is considered (i.e. tye is replaced by aofltyie)
               !sv         tau over water is used IN sbr ocwind
               +rhowaic*aofltyie(i,j)*dt*(sicomo(i,j)+sicomo(i,j+1))*0.5
          !svx 31.08.99
#else   
               +rhowaic*tye(i,j)*dt*(sicomo(i,j)+sicomo(i,j+1))*0.5 
#endif
          rv=(hibzeto(i,j)+hibeto(i,j))*sicve(i,j-1)/dydy                      &
               +   (hibzeto(i,j+1)+hibeto(i,j+1))*sicve(i,j+1)/dydy            &
               +(hibete(i,j)*sicve(i+1,j)+hibete(i-1,j)*sicve(i-1,j))/dxdx     &
               +(hibzeto(i,j)-hibeto(i,j))*(sicuo(i,j)-sicuo(i-1,j))/dxdy      &
               -(hibzeto(i,j+1)-hibeto(i,j+1))*(sicuo(i,j+1)-sicuo(i-1,j+1))   &
               /dxdy                                                           &
               +hibete(i,j)*(sicuo(i,j)-sicuo(i,j+1))/dxdy                     &
               -hibete(i-1,j)*(sicuo(i-1,j)-sicuo(i-1,j+1))/dxdy    
          !
          vh(i,j)=(rs+rv*dt)/(effice(i,j)                                      &
               +cweffe(i,j)*dt*speedv(i,j)                                     &
               +dt*effice(i,j)*g*dt*zschalt*effice(i,j)*rhoicwa*               &
               (1./(dlxp(i,j)*dlyp(i,j))+1./(dlxp(i,j+1)*dlyp(i,j+1)))         &
               *dlxv(i,j)/dlyv(i,j)                                            &
               +dt*(hibzeto(i,j)+hibeto(i,j)+hibzeto(i,j+1)+hibeto(i,j+1))/dydy&
               +dt*(hibete(i,j)+hibete(i-1,j))/dxdx)
          !
        ENDIF
      END DO
    END DO
!$OMP END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch('v',vh,'mo_ocice 131')
!$OMP END SINGLE
!#endif

!$OMP DO
    DO j=2,je1
      DO i=2,ie1
        IF(effico(i,j).GT.reval)THEN
          sicuo(i,j)=0.5*(uh(i,j)+sicuo(i,j))*amsuo(i,j,1)
        ELSE
          sicuo(i,j)=uko(i,j,1)
        ENDIF
        IF(effice(i,j).GT.reval)THEN
          sicve(i,j)=0.5*(vh(i,j)+sicve(i,j))*amsue(i,j,1)
        ELSE
          sicve(i,j)=vke(i,j,1)
        ENDIF
      END DO
    END DO
!$OMP END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch('u',sicuo,'mo_ocice 132')
    CALL bounds_exch('v',sicve,'mo_ocice 133')
!$OMP END SINGLE
!#endif

    !      write(io_stdout,*) ' iteration   u ',iter
    !      write(io_stdout,681)((amsuo(i,j,1)*sicuo(i,j),j=1,25),i=1,ie)
    !      write(io_stdout,*) ' iteration   v ',iter
    !      write(io_stdout,681)((amsue(i,j,1)*sicve(i,j),j=1,25),i=1,ie)

    
  END DO ! iteration_loop


!$OMP SINGLE
  CALL contro(-99999)
!$OMP END SINGLE


!$OMP DO
  DO j=2,je1
    DO i=2,ie1
      tauwatu(i,j)=cw*rhoicwa*speedu(i,j)*(sicuo(i,j)-uko(i,j,1))
      tauwatv(i,j)=cw*rhoicwa*speedv(i,j)*(sicve(i,j)-vke(i,j,1))
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL


#ifndef bounds_exch_tp
  IF(have_g_js) THEN
     tauwatu(:,1)=0.
     tauwatv(:,1)=0.
  ENDIF
  IF(have_g_je) THEN
     tauwatu(:,je)=0.
     tauwatv(:,je)=0.
  ENDIF
#endif

!#ifdef bounds_exch_save
  CALL bounds_exch('u',tauwatu,'mo_ocice 133')
  CALL bounds_exch('v',tauwatv,'mo_ocice 134')
!#endif

END SUBROUTINE ice_dynamics

SUBROUTINE ice_advection

  REAL uh(ie,je),vh(ie,je),wh(ie,je),zschalt



  !     advection of sea ice and snow
  !
  !cccc      IF (zschaltd.eq.0)THEN
  !
  !     upwind
  !

  INTEGER i,j

  REAL uwin,uwou,uein,ueou,vsin,vsou,vnin,vnou
  REAL area,areain

  zschalt=1.
  wh = 0

!$OMP PARALLEL private(i,j,uwin,uwou,uein,ueou,vsin,vsou,vnin,vnou,area,areain)

!$OMP DO
  DO j=2,je1
    DO i=2,ie1
      area=dlxp(i,j)*dlyp(i,j)
      areain=1./area
      uwin=0.5*(sicuo(i-1,j)+ABS(sicuo(i-1,j)))*dt*dlyu(i-1,j)
      uwou=0.5*(ABS(sicuo(i-1,j))-sicuo(i-1,j))*dt*dlyu(i-1,j)
      uein=0.5*(ABS(sicuo(i,j))-sicuo(i,j))*dt*dlyu(i,j)
      ueou=0.5*(ABS(sicuo(i,j))+sicuo(i,j))*dt*dlyu(i,j)
      vsin=0.5*(ABS(sicve(i,j))+sicve(i,j))*dt*dlxv(i,j)
      vsou=0.5*(ABS(sicve(i,j))-sicve(i,j))*dt*dlxv(i,j)
      vnin=0.5*(ABS(sicve(i,j-1))-sicve(i,j-1))*dt*dlxv(i,j-1)
      vnou=0.5*(ABS(sicve(i,j-1))+sicve(i,j-1))*dt*dlxv(i,j-1)

      uh(i,j)=(sictho(i,j)*(area-uwou-ueou-vsou-vnou)               &
             +sictho(i-1,j)*uwin+sictho(i+1,j)*uein+                  &
              sictho(i,j+1)*vsin+sictho(i,j-1)*vnin)*areain

        vh(i,j)=(sicomo(i,j)*(area-uwou-ueou-vsou-vnou)               &
             +sicomo(i-1,j)*uwin+sicomo(i+1,j)*uein+                  &
             sicomo(i,j+1)*vsin+sicomo(i,j-1)*vnin)*areain

        wh(i,j)=(sicsno(i,j)*(area-uwou-ueou-vsou-vnou)               &
             +sicsno(i-1,j)*uwin+sicsno(i+1,j)*uein+                  &
             sicsno(i,j+1)*vsin+sicsno(i,j-1)*vnin)*areain

    ENDDO
  ENDDO
!$OMP END DO

!!$!!!$  DO j=2,je1
!!$!!!$    DO i=2,ie1
!!$!!!$      uwin=0.5*(sicuo(i-1,j)+ABS(sicuo(i-1,j)))
!!$!!!$      uwou=0.5*(ABS(sicuo(i-1,j))-sicuo(i-1,j))
!!$!!!$      uein=0.5*(ABS(sicuo(i,j))-sicuo(i,j))
!!$!!!$      ueou=0.5*(ABS(sicuo(i,j))+sicuo(i,j))
!!$!!!$      vsin=0.5*(ABS(sicve(i,j))+sicve(i,j))
!!$!!!$      vsou=0.5*(ABS(sicve(i,j))-sicve(i,j))
!!$!!!$      vnin=0.5*(ABS(sicve(i,j-1))-sicve(i,j-1))
!!$!!!$      vnou=0.5*(ABS(sicve(i,j-1))+sicve(i,j-1))
!!$!!!$      uh(i,j)=sictho(i,j)*(1.-(uwou+ueou)*dt/dlxp(i,j)                  &
!!$!!!$           -dt*(vsou*dlxv(i,j)+vnou*dlxv(i,j-1))/(dlyp(i,j)*dlxp(i,j)))     &
!!$!!!$           + dt*(uwin*sictho(i-1,j)+uein*sictho(i+1,j))/dlxp(i,j)           &
!!$!!!$           +dt*( vsin*dlxv(i,j)*sictho(i,j+1)                               &
!!$!!!$           +    vnin*dlxv(i,j-1)*sictho(i,j-1))/(dlyp(i,j)*dlxp(i,j))
!!$!!!$      vh(i,j)=sicomo(i,j)*(1.-(uwou+ueou)*dt/dlxp(i,j)                  &
!!$!!!$           -dt*(vsou*dlxv(i,j)+vnou*dlxv(i,j-1))/(dlyp(i,j)*dlxp(i,j)))     &
!!$!!!$           + dt*(uwin*sicomo(i-1,j)+uein*sicomo(i+1,j))/dlxp(i,j)           &
!!$!!!$           +dt*( vsin*dlxv(i,j)*sicomo(i,j+1)                               &
!!$!!!$           +    vnin*dlxv(i,j-1)*sicomo(i,j-1))/(dlyp(i,j)*dlxp(i,j))
!!$!!!$      !as***add continuity of snow:
!!$!!!$      wh(i,j)=sicsno(i,j)*(1.-(uwou+ueou)*dt/dlxp(i,j)                  &
!!$!!!$           -dt*(vsou*dlxv(i,j)+vnou*dlxv(i,j-1))/(dlyp(i,j)*dlxp(i,j)))     &
!!$!!!$           + dt*(uwin*sicsno(i-1,j)+uein*sicsno(i+1,j))/dlxp(i,j)           &
!!$!!!$           +dt*( vsin*dlxv(i,j)*sicsno(i,j+1)                               &
!!$!!!$           +    vnin*dlxv(i,j-1)*sicsno(i,j-1))/(dlyp(i,j)*dlxp(i,j))
!!$!!!$    END DO
!!$!!!$  END DO

!!$  ! no bounds exchange necessary on uh, vh, wh since they are used only
!!$  ! IN the inner domain below
!$OMP DO
  DO j=2,je1
    DO i=2,ie1
!uwe  include effect on sealevel zo
      zo(i,j)=zo(i,j)+(uh(i,j)-sictho(i,j))*rhoicwa*zschalt             &
                     +(wh(i,j)-sicsno(i,j))*rhosnwa*zschalt
      sictho(i,j)=uh(i,j)
      sicomo(i,j)=vh(i,j)

      !as***add snow:
      sicsno(i,j)=wh(i,j)
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL


  CALL bounds_exch('p',zo,'mo_ocice 135')
  CALL bounds_exch('p',sictho,'mo_ocice 136')
  CALL bounds_exch('p',sicomo,'mo_ocice 137')
!#ifdef bounds_exch_save
  CALL bounds_exch('p',sicsno,'mo_ocice 138')
!#endif
  CALL contro(-9998)




END SUBROUTINE ice_advection

SUBROUTINE ice_thermodynamics(siotho, sioomo, siosno)

  USE MO_TRO

#ifdef FB_BGC_OCE
#ifdef PBGC  
  USE mo_biomod, ONLY: abs_oce
#else 
  USE mo_attenmap
#endif   
#endif

  REAL, INTENT (INOUT) :: siotho(ie,je),sioomo(ie,je),siosno(ie,je)

#ifdef FB_BGC_OCE 
  REAL swsum(ie,je),swsumi(ie,je)
  REAL swrab(ie,je,ke)
  REAL heatabs(ie,je),heatabb(ie,je)
#else
  REAL swsum, swsumi
  REAL swrab(ke),heatabs(ie,je),heatabb(ie,je)
#endif      

  INTEGER i,j,k,l
  REAL opendep
  heatabs(:,:)=0.
  heatabb(:,:)=0.


!$OMP PARALLEL private(i,j,k)

!     length scale for penetration of sw-radiation [m]
#ifdef FB_BGC_OCE   
  call bounds_exch('p',abs_oce,'in ice_thermodynamics')
!$OMP DO
  DO j=1,je
    DO i=1,ie
      swsum(i,j)=abs_oce(i,j,2)
      swsumi(i,j)=1./abs_oce(i,j,2)
    END DO
  END DO
!$OMP END DO
#else	       	  
  opendep=11.
#ifdef OPEND55
  opendep=5.5
#endif
!uwe    new sw-penetration
  swsum=EXP(-tiestw(2)/opendep)
  swsumi=1./swsum
#endif /* FB_BGC_OCE */  

#ifdef FB_BGC_OCE  
  DO k=1,ke-1
!$OMP DO
    DO j=1,je
      DO i=1,ie
        swrab(i,j,k)=swsumi(i,j)*(abs_oce(i,j,k)-abs_oce(i,j,k+1))
      END DO
    END DO
!$OMP END DO
  END DO
  k=ke
!$OMP DO
  DO j=1,je
    DO i=1,ie
      swrab(i,j,k)=swsumi(i,j)*abs_oce(i,j,k) 
    END DO
  END DO
!$OMP END DO
#else	       	 
  DO k=1,ke
    swrab(k)=swsumi*(EXP(-tiestw(k)/opendep)                         &
         -EXP(-tiestw(k+1)/opendep))
  ENDDO
#endif      

!$OMP END PARALLEL

#ifdef SAOCLOSE
  CALL dilcor_gtrf
#endif 


#ifdef PBGC
  CALL dilcor_gtrf2
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef ADDCONTRA
  CALL dilcor_gtrf2
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifdef __coupled
  !sv 31.08.99 included option for forcing with fluxes 
  !sv 6.9.99 set heat flux correction to zero for the time being

  aofldhwo(:,:) = 0.


  CALL growth(siotho,sioomo,siosno,sictho,sicomo,sicsno,            &
       aoflfrwo,aoflfrio,aofldhwo,aoflnhwo,aoflchio,aoflrhio,       &
       sao,tho,dzw(1),zo,dt,fwo,weto,                               &
       preco,prech,                                                 &
       aoflshwo,qswo,qlwo,qseo,qlao,heatabs,swsum)


#else /*__coupled*/
!#ifdef  bounds_exch_save
CALL bounds_exch('p',alat,'mo_ocice 139')
CALL bounds_exch('p',siotho,'mo_ocice 140')
CALL bounds_exch('p',sioomo,'mo_ocice 141')
CALL bounds_exch('p',siosno,'mo_ocice 142')
CALL bounds_exch('p',sictho,'mo_ocice 143')
CALL bounds_exch('p',sicomo,'mo_ocice 144')
CALL bounds_exch('p',sicsno,'mo_ocice 145')
CALL bounds_exch('p',rpreco,'mo_ocice 146')
CALL bounds_exch('p',tairo,'mo_ocice 147')
CALL bounds_exch('p',tdo,'mo_ocice 148')
CALL bounds_exch('p',aclo,'mo_ocice 149')
CALL bounds_exch('p',pao,'mo_ocice 150')
CALL bounds_exch('p',fu10,'mo_ocice 151')
CALL bounds_exch('p',fswr,'mo_ocice 152')
CALL bounds_exch('p',sao,'mo_ocice 153')
CALL bounds_exch('p',tho,'mo_ocice 154')
CALL bounds_exch('p',zo,'mo_ocice 155')
CALL bounds_exch('p',fwo,'mo_ocice 156')
CALL bounds_exch('p',sho,'mo_ocice 157')
CALL bounds_exch('p',weto,'mo_ocice 158')
CALL bounds_exch('p',qswo,'mo_ocice 159')
CALL bounds_exch('p',qlwo,'mo_ocice 160')
CALL bounds_exch('p',qseo,'mo_ocice 161')
CALL bounds_exch('p',qlao,'mo_ocice 162')
CALL bounds_exch('p',heatabs,'mo_ocice 165')
!#endif
CALL bounds_exch('p',preco,'mo_ocice 163')
CALL bounds_exch('p',prech,'mo_ocice 164')



  CALL growth(alat,siotho,sioomo,siosno,sictho,sicomo,sicsno,       &
       rpreco,tairo,tdo,aclo,pao,fu10,fswr,                         &
       sao,tho,dzw(1),zo,dt,fwo,sho,weto,                           &
       qswo,qlwo,qseo,qlao,preco,prech,heatabs,swsum)


#endif /*__coupled*/
!#ifdef bounds_exch_save
CALL bounds_exch('p',siotho,'mo_ocice 166')
CALL bounds_exch('p',sioomo,'mo_ocice 167')
CALL bounds_exch('p',siosno,'mo_ocice 168')
CALL bounds_exch('p',sictho,'mo_ocice 169')
CALL bounds_exch('p',sicomo,'mo_ocice 170')
CALL bounds_exch('p',sicsno,'mo_ocice 171')
CALL bounds_exch('p',rpreco,'mo_ocice 172')
CALL bounds_exch('p',tairo,'mo_ocice 173')
CALL bounds_exch('p',tdo,'mo_ocice 174')
CALL bounds_exch('p',aclo,'mo_ocice 175')
CALL bounds_exch('p',pao,'mo_ocice 176')
CALL bounds_exch('p',fu10,'mo_ocice 177')
CALL bounds_exch('p',fswr,'mo_ocice 178')
CALL bounds_exch('p',sao,'mo_ocice 179')
CALL bounds_exch('p',tho,'mo_ocice 180')
CALL bounds_exch('p',zo,'mo_ocice 181')
CALL bounds_exch('p',fwo,'mo_ocice 182')
CALL bounds_exch('p',sho,'mo_ocice 183')
CALL bounds_exch('p',weto,'mo_ocice 184')
CALL bounds_exch('p',qswo,'mo_ocice 185')
CALL bounds_exch('p',qlwo,'mo_ocice 186')
CALL bounds_exch('p',qseo,'mo_ocice 187')
CALL bounds_exch('p',qlao,'mo_ocice 188')
CALL bounds_exch('p',preco,'mo_ocice 189')
CALL bounds_exch('p',prech,'mo_ocice 190')
CALL bounds_exch('p',heatabs,'mo_ocice 191')
!#endif

#ifdef PBGC
   do l=1,nocetra
      CALL dilcor_ptrf2(ocetra(1,1,1,l))
   enddo
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef ADDCONTRA
   do l=1,nocectra
      CALL dilcor_ptrf2(ocectra(1,1,1,l))
   enddo
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifdef SAOCLOSE
     CALL dilcor_ptrf
#endif

!$OMP PARALLEL PRIVATE(I,J,K)

!$OMP DO
  DO j=1,je
    DO i=1,ie
      heatabs(i,j)=weto(i,j,1)*heatabs(i,j)*dt/cc
      heatabb(i,j)=0.
    END DO
  END DO   
!$OMP END DO


#ifdef FB_BGC_OCE  
  DO k=ke,2,-1
!$OMP DO
     DO j=1,je
        DO i=1,ie
           heatabb(i,j)=heatabb(i,j)+swrab(i,j,k)*heatabs(i,j)
           tho(i,j,k)=tho(i,j,k)+heatabb(i,j)*dpio(i,j,k)
           heatabb(i,j)=heatabb(i,j)*(1.-weto(i,j,k))
        ENDDO
     END DO
!$OMP END DO
  END DO  ! k-loop

#else /*FB_BGC_OCE*/

  DO k=ke,2,-1
!$OMP DO
     DO j=1,je
        DO i=1,ie
           heatabb(i,j)=heatabb(i,j)+swrab(k)*heatabs(i,j)
           tho(i,j,k)=tho(i,j,k)+heatabb(i,j)*dpio(i,j,k)
           heatabb(i,j)=heatabb(i,j)*(1.-weto(i,j,k))
        ENDDO
     END DO
!$OMP END DO
  END DO  ! k-loop

#endif /*FB_BGC_OCE*/      

!$OMP END PARALLEL


  CALL contro(-9997)


END SUBROUTINE ice_thermodynamics


END MODULE mo_ocice

