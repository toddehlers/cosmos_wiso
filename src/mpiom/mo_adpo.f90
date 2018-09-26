MODULE MO_ADPO

  USE MO_PARAM1

#ifdef bounds_exch_tp
  REAL ,ALLOCATABLE :: TRPHELP(:,:),TRMHELP(:,:)
#endif 
  real ,allocatable :: woh1(:,:,:),woh2(:,:,:)
  real, allocatable :: svo(:,:,:),dlxypa(:,:)

CONTAINS

SUBROUTINE ALLOC_MEM_ADPO

#ifdef bounds_exch_tp
   ALLOCATE(TRPHELP(IE,KE),TRMHELP(IE,KE))
#endif 
  allocate(woh1(ie,je,kep),woh2(ie,je,kep))
  allocate(dlxypa(ie,je),svo(ie,je,ke))

END SUBROUTINE ALLOC_MEM_ADPO



SUBROUTINE ocadpo_base
  !
  !     Computes bottom boundary layer velocities (iocad.eq.4)
  !     for each timestep to be used in sbrt ocadpo_trf
  !
  !     By Stephan Lorenz 08/2007
  !

  USE mo_param1
  USE mo_parallel
  USE mo_commo1
  USE mo_commoau1
  USE mo_commobbl

!$OMP PARALLEL PRIVATE(i,j,k,suminf)

!$OMP WORKSHARE
   dlxypa = dlxp * dlyp
!$OMP END WORKSHARE

!$OMP WORKSHARE
   woh1 = wo
   woh2 = wo
!  WOBACK(I,J,K)=WO(I,J,K)
!$OMP END WORKSHARE

!$OMP DO
  DO k=1,ke
    DO j=1,je
      DO i=1,ie
        svo(i,j,k)=weto(i,j,k)* dlxypa(i,j) *ddpo(i,j,k)
      ENDDO
    ENDDO
  ENDDO

!     COMPUTE NEW VERTICAL VELOCITIES INCLUDING BBL TRANSPORTS

  if (iocad==4) then

!$OMP DO
    DO J=2,JE1
      DO K=KE,1,-1
        DO I=2,IE1
          SUMINF=0.
        
          IF (KDWUBBL(I-1,J).EQ.K) SUMINF=SUMINF+UBBL(I-1,J)
          IF (KUPUBBL(I-1,J).EQ.K) SUMINF=SUMINF-UBBL(I-1,J)
          IF (KDWUBBL(I,J).EQ.K) SUMINF=SUMINF-UBBL(I,J)
          IF (KUPUBBL(I,J).EQ.K) SUMINF=SUMINF+UBBL(I,J)
                               
          IF (KDWVBBL(I,J).EQ.K) SUMINF=SUMINF+VBBL(I,J)
          IF (KUPVBBL(I,J).EQ.K) SUMINF=SUMINF-VBBL(I,J)
          IF (KDWVBBL(I,J-1).EQ.K) SUMINF=SUMINF-VBBL(I,J-1)
          IF (KUPVBBL(I,J-1).EQ.K) SUMINF=SUMINF+VBBL(I,J-1)
                              
          woh1(I,J,K)=woh1(I,J,K+1)+SUMINF/dlxypa(I,J)
          
        ENDDO
      ENDDO
    ENDDO
 
!$OMP DO
    DO K=1,KE
      DO J=2,JE1
        DO I=2,IE1
!         WO(I,J,K)=WOBACK(I,J,K)+WO(I,J,K)
          woh2(I,J,K)=woh1(I,J,K)+woh2(I,J,K)
        ENDDO
      ENDDO
    ENDDO

!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch('p',woh2,'ocadpo_base 1')
!$OMP END SINGLE
!#endif

  endif ! iocad==4

!     Compute old volume of surface layer according to old sealevel 
!       before update of ZO in update_zo
  
!$OMP DO
  DO J=1,JE
    DO I=1,IE
      svo(I,J,1)=svo(I,J,1)+ dlxypa(i,j) *WETO(I,J,1)*              &
           (ZO(I,J)-woh2(I,J,1)*DT-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)
    ENDDO
  ENDDO

!$OMP END PARALLEL  
  
END SUBROUTINE ocadpo_base

SUBROUTINE ocadpo_trf(trf)

  !     ROUTINE OCADPO
  !
  !     COMPUTES ADVECTION OF ONE TRACERFIELD
  !
  !     BY ERNST MAIER-REIMER 1/1999
  !     MODIFIED 1/2000 UWE MIKOLAJEWICZ
  !     MAKE MASS CONSERVING WITH FREE SURFACE LAYER
  !     MODIFIED 3/2000 UWE MIKOLAJEWICZ
  !     INCLUDE BOTTOM BOUNDARY LAYER TRANSPORTS IN
  !     ADVECTION SCHEME (IOCAD.EQ.4 :: ADPO + SLOPECON_ADPO)
  !     NEEDS TO BE RUN TOGETHER WITH SR SLOPETRANS
  !
  !     INPUT/OUTPUT  
  !     TRF(IE,JE,KE)     TRACER FIELD 
  !
  !     USES VELOCITIES (UKO,VKE,WO)
  !
  !     Changes by S. Lorenz, J.-O. Biesmann, 2007-08 optimised
  !       vertical velocities updated in ocadpo_base

  USE mo_param1
  USE mo_parallel
  USE mo_commo1
  USE mo_commoau1
  USE mo_commobbl

  REAL :: trf(ie,je,ke)
  REAL :: s2o(ie,je,ke)
  REAL :: trp(ie,je,kep), trm(ie,je,0:ke)
  REAL :: wtp(ie,je,kep), wtm(ie,je,0:ke)
  REAL :: b1u, b1um, b1v, b1vm, zzsurf


#ifndef TRACER_OMP
!$OMP PARALLEL PRIVATE(i,j,k,klo,kup,suminf,abl,zwabl,up,um,scha,sch,schi,rn,   &
!$OMP                  zzsurf,b1u,b1um,b1v,b1vm)
#endif

!!#ifdef bounds_exch_save
!!$OMP SINGLE
!! CALL bounds_exch('p',trf,'ocadpo_trf 1') 
!!$OMP END SINGLE
!!#endif

! #slo# - initialise local variables
#ifndef TRACER_OMP
!$OMP WORKSHARE
#endif
   s2o(:,:,:)=0.
   trp(:,:,:)=0.
   trm(:,:,:)=0.
   wtp(:,:,:)=0.
   wtm(:,:,:)=0.
#ifndef TRACER_OMP
!$OMP END WORKSHARE
#endif

if ( iocad == 1 ) then

!
! Vertical, Zonal, and Meridional Transports for iocad == 1 - Pure Upwind
!

!      VERTICAL TRANSPORTS

#ifndef TRACER_OMP
!$OMP DO
#endif
  DO J=2,JE1
    DO k=1,ke
      zzsurf = 1.
      if (k == 1) zzsurf=0.
      DO I=2,IE1
        UP=0.5*DT * dlxypa(i,j) * (woh2(I,J,K)+ABS(woh2(I,J,K)))*zzsurf
        UM=0.5*DT * dlxypa(i,j) * (ABS(woh2(I,J,K+1))-woh2(I,J,K+1))
        TRP(I,J,K)=UP*TRF(I,J,K)
        TRM(I,J,K)=UM*TRF(I,J,K)
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      END DO
    END DO
  END DO

#ifndef TRACER_OMP
!$OMP DO
#endif
  DO J=2,JE1
    DO K=1,KE
      DO I=2,IE1
        IF(WETO(I,J,K).GT.0.5) THEN
          RN=svo(I,J,K)+WTP(I,J,K+1)-WTP(I,J,K)-WTM(I,J,K)+WTM(I,J,K-1)
          TRF(I,J,K)=(TRF(I,J,K)*svo(I,J,K)+TRP(I,J,K+1)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I,J,K-1))/RN
          s2o(I,J,K)=RN
        ENDIF
      END DO
    END DO
  END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('p',TRF,'ocadpo_trf 3')   ! #slo#: not necessary?
  CALL bounds_exch('p',s2o,'ocadpo_trf 4')   ! #slo#: not necessary?
!$OMP END SINGLE
!#endif

!      TRANSPORTS IN X-DIRECTION
    
#ifndef TRACER_OMP
!$OMP DO
#endif
  DO J=2,JE1
    DO K=1,KE
      DO I=2,IE1

        UP=0.5*DT*DDUO(I,J,K)*DLYU(I,J)*(UKO(I,J,K)+ABS(UKO(I,J,K)))

        UM=0.5*DT*DDUO(I-1,J,K)*(ABS(UKO(I-1,J,K))-UKO(I-1,J,K))         &
                *DLYU(I-1,J)

        TRP(I,J,K)=UP*TRF(I,J,K)
        TRM(I,J,K)=UM*TRF(I,J,K)
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      END DO
    END DO
  END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('u+',TRP,'ocadpo_trf 5')
  CALL bounds_exch('u+',TRM,'ocadpo_trf 6')
  CALL bounds_exch('u+',WTP,'ocadpo_trf 7')
  CALL bounds_exch('u+',WTM,'ocadpo_trf 8')
!$OMP END SINGLE
!#endif  

#ifndef TRACER_OMP
!$OMP DO
#endif
  DO J=2,JE1
    DO K=1,KE
      DO I=2,IE1
        IF(WETO(I,J,K).GT.0.5) THEN   
          RN=s2o(I,J,K)+WTP(I-1,J,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I+1,J,K)
          TRF(I,J,K)=(TRF(I,J,K)*s2o(I,J,K)+TRP(I-1,J,K)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I+1,J,K))/RN
          s2o(I,J,K)=RN
        ENDIF
      ENDDO
    ENDDO
  END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('p',TRF,'ocadpo_trf 9')    ! #slo#: not necessary?
  CALL bounds_exch('p',s2o,'ocadpo_trf 10')   ! #slo#: not necessary?
!$OMP END SINGLE
!#endif 

!      TRANSPORTS IN Y-DIRECTION

#ifndef TRACER_OMP
!$OMP DO
#endif
  DO J=2,JE1
    DO K=1,KE
      DO I=2,IE1
        UP=0.5*DT*DDUE(I,J-1,K)*DLXV(I,J-1)                              &
             *(VKE(I,J-1,K)+ABS(VKE(I,J-1,K)))
        UM=0.5*DT*DDUE(I,J,K)*DLXV(I,J)*(ABS(VKE(I,J,K))-VKE(I,J,K))

        TRP(I,J,K)=UP*TRF(I,J,K)
        TRM(I,J,K)=UM*TRF(I,J,K)
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      END DO
    END DO
  END DO

!
! Vertical, Zonal, and Meridional Transports for iocad == 3, 4
!

ELSE     ! i.e. iocad==3 or iocad==4

!      VERTICAL TRANSPORTS

#ifndef TRACER_OMP
!$OMP DO
#endif
    DO J=2,JE1
! k=1
      DO I=2,IE1
        UM=0.5*DT * dlxypa(i,j) * (ABS(woh2(I,J,2))-woh2(I,J,2))
        TRP(I,J,1)=0.
        TRM(I,J,1)=UM*TRF(I,J,1)
        WTP(I,J,1)=0.
        WTM(I,J,1)=UM
      END DO
! k=ke
      DO I=2,IE1
        UP=0.5*DT * dlxypa(i,j) * (woh2(I,J,KE)+ABS(woh2(I,J,KE)))
        UM=0.5*DT * dlxypa(i,j) * (ABS(woh2(I,J,KE+1))-woh2(I,J,KE+1))
        TRP(I,J,KE)=UP*TRF(I,J,KE)
        TRM(I,J,KE)=UM*TRF(I,J,KE)
        WTP(I,J,KE)=UP
        WTM(I,J,KE)=UM
      END DO

    DO K=2,KE1
      KLO=MIN(K+1,KE)
      KUP=MAX(K-1,1)

      DO I=2,IE1
        ABL=ABS(TRF(I,J,KUP)-TRF(I,J,KLO))
        ZWABL=ABS(TRF(I,J,KUP)+TRF(I,J,KLO)-2.*TRF(I,J,K))
        UP=0.5*DT * dlxypa(i,j) * (woh2(I,J,K)+ABS(woh2(I,J,K)))
        UM=0.5*DT * dlxypa(i,j) * (ABS(woh2(I,J,K+1))-woh2(I,J,K+1))
        SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))

#ifdef SMOADV
        SCH=MIN(WETO(I,J,KLO),SCHA)
#else
        SCH=MIN(WETO(I,J,KLO),SCHA*svo(I,J,K)/(UP+UM+1.E-20))
#endif

        SCHI=1.-SCH  

        TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J,KUP)))
        TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J,KLO)))
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      END DO
    END DO

  END DO


#ifndef TRACER_OMP
!$OMP DO
#endif
  DO J=2,JE1
    DO K=1,KE
      DO I=2,IE1
        IF(WETO(I,J,K).GT.0.5) THEN
          RN=svo(I,J,K)+WTP(I,J,K+1)-WTP(I,J,K)-WTM(I,J,K)+WTM(I,J,K-1)
          TRF(I,J,K)=(TRF(I,J,K)*svo(I,J,K)+TRP(I,J,K+1)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I,J,K-1))/RN
          s2o(I,J,K)=RN
        ENDIF
      END DO
    END DO
  END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('p',TRF,'ocadpo_trf 3')   ! #slo#: not necessary?
  CALL bounds_exch('p',s2o,'ocadpo_trf 4')   ! #slo#: not necessary?
!$OMP END SINGLE
!#endif

!      TRANSPORTS IN X-DIRECTION
    
#ifndef TRACER_OMP
!$OMP DO
#endif
  DO J=2,JE1
    DO K=1,KE

      DO I=2,IE1
        b1u  = 0.
        b1um = 0.
         if ( iocad == 4 ) then
           IF (KDWUBBL(I,J).EQ.K) b1u=UBBL(I,J)
           IF (KUPUBBL(I,J).EQ.K) b1u=-UBBL(I,J)
           IF (KDWUBBL(I-1,J).EQ.K) b1um=UBBL(I-1,J)
           IF (KUPUBBL(I-1,J).EQ.K) b1um=-UBBL(I-1,J)
         endif
        ABL=ABS(TRF(I+1,J,K)-TRF(I-1,J,K))
        ZWABL=ABS(TRF(I+1,J,K)+TRF(I-1,J,K)-2.*TRF(I,J,K))

        UP=0.5*DT*DDUO(I,J,K)*DLYU(I,J)*(UKO(I,J,K)+ABS(UKO(I,J,K)))     &
             + 0.5*DT*(b1u+ABS(b1u))

        UM=0.5*DT*DDUO(I-1,J,K)*(ABS(UKO(I-1,J,K))-UKO(I-1,J,K))         &
             *DLYU(I-1,J)                                                &
             + 0.5*DT*(ABS(b1um)-b1um)

        SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))*WETO(I,J,K)                &
             *WETO(I+1,J,K)*WETO(I-1,J,K)
#ifdef SMOADH
        SCH=MIN(1.,SCHA)
#else
        SCH=MIN(1.,SCHA*s2o(I,J,K)/(UP+UM+1.E-20))
#endif

        SCHI=1.-SCH

        TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I+1,J,K)))
        TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I-1,J,K)))
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      END DO
    END DO

  END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('u+',TRP,'ocadpo_trf 5')
  CALL bounds_exch('u+',TRM,'ocadpo_trf 6')
  CALL bounds_exch('u+',WTP,'ocadpo_trf 7')
  CALL bounds_exch('u+',WTM,'ocadpo_trf 8')
!$OMP END SINGLE
!#endif  

#ifndef TRACER_OMP
!$OMP DO
#endif
  DO J=2,JE1
    DO K=1,KE
      DO I=2,IE1
        IF(WETO(I,J,K).GT.0.5) THEN   
          RN=s2o(I,J,K)+WTP(I-1,J,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I+1,J,K)
          TRF(I,J,K)=(TRF(I,J,K)*s2o(I,J,K)+TRP(I-1,J,K)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I+1,J,K))/RN
          s2o(I,J,K)=RN
        ENDIF
      ENDDO
    ENDDO
  END DO

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('p',TRF,'ocadpo_trf 9')    ! #slo#: not necessary?
  CALL bounds_exch('p',s2o,'ocadpo_trf 10')   ! #slo#: not necessary?
!$OMP END SINGLE
!#endif 

!      TRANSPORTS IN Y-DIRECTION

#ifndef TRACER_OMP
!$OMP DO
#endif
  DO J=2,JE1
    DO K=1,KE

#ifdef bounds_exch_tp
       if(p_joff.eq.0)then
          do i=1,ie
             trf(i,1,k)=2.*trf(i,2,k)-trf(i,3,k)
          enddo
       endif
#endif

      DO I=2,IE1
              b1v = 0.
              b1vm = 0.
              if ( iocad == 4 ) then
                IF (K.EQ.KDWVBBL(I,J)) b1v=VBBL(I,J)
                IF (K.EQ.KUPVBBL(I,J)) b1v=-VBBL(I,J)
                IF (K.EQ.KDWVBBL(I,J-1)) b1vm=VBBL(I,J-1)
                IF (K.EQ.KUPVBBL(I,J-1)) b1vm=-VBBL(I,J-1)
              endif

        ABL=ABS(TRF(I,J-1,K)-TRF(I,J+1,K))
        ZWABL=ABS(TRF(I,J-1,K)+TRF(I,J+1,K)-2.*TRF(I,J,K))
        UP=0.5*DT*DDUE(I,J-1,K)*DLXV(I,J-1)                              &
             *(VKE(I,J-1,K)+ABS(VKE(I,J-1,K)))                           &
                + 0.5*DT*(b1vm+ABS(b1vm))
        UM=0.5*DT*DDUE(I,J,K)*DLXV(I,J)*(ABS(VKE(I,J,K))-VKE(I,J,K))     &
                + 0.5*DT*(ABS(b1v)-b1v)
        SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))*WETO(I,J,K)                &
             *WETO(I,J-1,K)*WETO(I,J+1,K)
#ifdef SMOADH
        SCH=MIN(1.,SCHA)
#else
        SCH=MIN(1.,SCHA*s2o(I,J,K)/(UP+UM+1.E-20))
#endif

        SCHI=1.-SCH
        TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J-1,K)))
        TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J+1,K)))
        WTP(I,J,K)=UP
        WTM(I,J,K)=UM
      END DO
    END DO
  END DO

ENDIF    ! iocad

!#ifdef bounds_exch_save  
!$OMP SINGLE
  CALL bounds_exch('vv',TRP,'ocadpo_trf 11')
  CALL bounds_exch('vv',TRM,'ocadpo_trf 12')
  CALL bounds_exch('vv',WTP,'ocadpo_trf 13')
  CALL bounds_exch('vv',WTM,'ocadpo_trf 14')
!$OMP END SINGLE
!#endif


#ifndef TRACER_OMP
!$OMP DO
#endif
  DO J=2,JE1
    DO K=1,KE
      DO I=2,IE1
        IF(WETO(I,J,K).GT.0.5) THEN
          
          RN=s2o(I,J,K)+WTP(I,J+1,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I,J-1,K)
          TRF(I,J,K)=(TRF(I,J,K)*s2o(I,J,K)+TRP(I,J+1,K)-TRP(I,J,K)      &
               -TRM(I,J,K)+TRM(I,J-1,K))/RN
          s2o(I,J,K)=RN
        ENDIF
      ENDDO
    ENDDO
  END DO

!#ifdef bounds_exch_save  
!$OMP SINGLE
  CALL bounds_exch('p',TRF,'ocadpo_trf 15')
!$OMP END SINGLE
!#endif 

#ifndef TRACER_OMP
!$OMP END PARALLEL  
#endif

END SUBROUTINE ocadpo_trf

END MODULE MO_ADPO
