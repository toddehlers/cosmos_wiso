SUBROUTINE octdiff_trf(trf)


  USE mo_param1
  USE mo_parallel
  USE mo_commo1
  USE mo_levitus
  USE mo_commoau1
  USE mo_units
  USE mo_octdiff
  USE mo_mpi


  !     COMPUTES DIFFUSION OF TEMPERATURE trf AND SALINITY SAO
  !
  !     UWE MIKOLAJEWICZ 12/99
  !
  !     VERTICAL DIFFUSION IMPLICITE
  !     HORIZONTAL DIFFUSION BOTH HARMONIC AND BIHARMONIC, BOTH EXPLICITE
  !                                AH00         AULAPTS
  !     ACTUAL COEFFICIENTS SPATIALLY VARIABLE WITH RESOLUTION
  !                                AHX,AHY      AULX,AULY
  !
  !     VARIABLE THICKNESS OF SURFACE LAYER NOW INCLUDED
  !
  !     Changes R. Johanni, 2003-11-10:
  !     DSLOPX, DSLOPY, DVSLOP were calculated but nowhere used -> removed
  !
  !     Changes R. Smith, 2004-09-27
  !     tracer-independent matrices calculated separately first for use
  !     with multiple tracers
  !
  !     Changes S. J. Lorenz, J.-O. Biesmann, 08/2007, optimisation

  IMPLICIT NONE

  INTEGER i, j, k

  REAL tv_diff, tfluz, zzsur0
  REAL dlxyp(ie,je), dlxyp1(ie,je), dlxyp1j(ie,je)
  REAL dlxui(ie,je), dlyvi(ie,je)
  REAL tflux(ie,je,ke), tfluy(ie,je,ke)
  REAL trf(ie,je,ke), ten(ie,je,ke), t2o(ie,je,ke)


!#ifdef bounds_exch_save
!#slo#: Here - no bounds_exch necessary
!! !$OMP SINGLE
!!   CALL bounds_exch('p',trf,'octdiff_trf 1')
!! !$OMP END SINGLE
!#endif

#ifndef TRACER_OMP
!$OMP PARALLEL PRIVATE(i,j,k,tv_diff,tfluz,zzsur0)

!$OMP WORKSHARE
#endif
  t2o(:,:,:)=trf(:,:,:)*vol_term(:,:,:)
#ifndef TRACER_OMP
!$OMP END WORKSHARE
#endif

  IF (ah00.GT.almzer) THEN

#ifndef TRACER_OMP
!$OMP WORKSHARE
#endif
     tflux(:,:,:)=0.0
     tfluy(:,:,:)=0.0
#ifndef TRACER_OMP
!$OMP END WORKSHARE
#endif


     !X DIRECTION - tracer dependent

#ifdef ISOPYK

!  prepare for multiplication using inverse geometry

#ifndef TRACER_OMP
!$OMP DO
#endif           
        DO j=1,je
           DO i=1,ie
              dlxyp(i,j)=dlxp(i,j)*dlyp(i,j)
              dlxui(i,j) = 1. / dlxu(i,j)
              dlyvi(i,j) = 1. / dlyv(i,j)
           ENDDO
           DO i=1,ie1
              dlxyp1(i,j)=dlxp(i+1,j)*dlyp(i+1,j)
           ENDDO
        ENDDO

#ifndef TRACER_OMP
!$OMP DO
#endif           
        DO j=1,je1
           DO i=1,ie
              dlxyp1j(i,j)=dlxp(i,j+1)*dlyp(i,j+1)
           ENDDO
        ENDDO


#ifndef TRACER_OMP
!$OMP DO
#endif           
     DO j=1,je
        DO k=1,ke
           DO i=1,ie1
              IF(amsuo(i,j,k).GT.0.)THEN
                 tv_diff=(trf(i+1,j,k)-trf(i,j,k))*dlxui(i,j)
                 !
                 !       triangle left,upw
                 !
                 tflux(i,j,k)=tflux(i,j,k)+xcoeff_lo(i,j,k)*(tv_diff  &
                      +sloplox(i,j,k)*                                &
                      (trf(i,j,ko(j,k))-trf(i,j,k))*di(k))
                 tfluz=zcoeff_lox(i,j,k)*(sloplox(i,j,k)*             &
                      (trf(i,j,ko(j,k))-trf(i,j,k))*                  &
                      dlxyp(i,j)*di(k)+1.*dlxyp(i,j)*tv_diff)

                 t2o(i,j,ko(j,k))=t2o(i,j,ko(j,k))-tfluz
                 t2o(i,j,k)=t2o(i,j,k)+tfluz
                 !
                 !       triangle left down
                 !
                 tflux(i,j,k)=tflux(i,j,k)+xcoeff_lu(i,j,k)*(tv_diff  &
                      +sloplux(i,j,k)*                                &
                      (trf(i,j,k)-trf(i,j,ku(j,k)))*di(ku(j,k)))            
                 tfluz=zcoeff_lux(i,j,k)*(sloplux(i,j,k)*             &
                      (trf(i,j,k)-trf(i,j,ku(j,k)))*                  &
                      dlxyp(i,j)*di(ku(j,k))+1.*dlxyp(i,j)*tv_diff)

                 t2o(i,j,k)=t2o(i,j,k)-tfluz
                 t2o(i,j,ku(j,k))=t2o(i,j,ku(j,k))+tfluz

              ENDIF
           END DO
        END DO
     END DO


#ifndef TRACER_OMP
!$OMP DO
#endif           
     DO j=1,je
        DO k=1,ke
           DO i=1,ie1
              IF(amsuo(i,j,k).GT.0.)THEN
                 tv_diff=(trf(i+1,j,k)-trf(i,j,k))*dlxui(i,j)
                 !
                 !       triangle right,upw
                 !
                 tflux(i,j,k)=tflux(i,j,k)+xcoeff_ro(i,j,k)*(tv_diff &
                      +sloprox(i,j,k)*                               &
                      (trf(i+1,j,ko(j,k))-trf(i+1,j,k))*di(k))
                 tfluz=zcoeff_rox(i,j,k)*(sloprox(i,j,k)*            &
                      (trf(i+1,j,ko(j,k))-trf(i+1,j,k))*             &
                      dlxyp1(i,j)*di(k)+1.*dlxyp1(i,j)*tv_diff)

                 t2o(i+1,j,ko(j,k))=t2o(i+1,j,ko(j,k))-tfluz
                 t2o(i+1,j,k)=t2o(i+1,j,k)+tfluz
                 !
                 !       triangle right down
                 !
                 tflux(i,j,k)=tflux(i,j,k)+xcoeff_ru(i,j,k)*(tv_diff &
                      +sloprux(i,j,k)*                               &
                      (trf(i+1,j,k)-trf(i+1,j,ku(j,k)))*di(ku(j,k)))
                 tfluz=zcoeff_rux(i,j,k)*(sloprux(i,j,k)*            &
                      (trf(i+1,j,k)-trf(i+1,j,ku(j,k)))*             &
                      dlxyp1(i,j)*di(ku(j,k))+1.*dlxyp1(i,j)*tv_diff)

                 t2o(i+1,j,k)=t2o(i+1,j,k)-tfluz
                 t2o(i+1,j,ku(j,k))=t2o(i+1,j,ku(j,k))+tfluz

              ENDIF
           END DO
        END DO
     END DO

#else /*ISOPYK*/

#ifndef TRACER_OMP
!$OMP DO
#endif           
     DO j=1,je
        DO k=1,ke
           DO i=1,ie1
              IF(amsuo(i,j,k).GT.0.)THEN
                 tv_diff=(trf(i+1,j,k)-trf(i,j,k))*dlxui(i,j)
                 tflux(i,j,k)=tv_diff*xflux(i,j,k)
              ENDIF
           END DO
        END DO
     END DO

#endif /*ISOPYK*/

!$OMP SINGLE
     call bounds_exch('u',tflux,'octdiff_trf 2 ')
!$OMP END SINGLE

#ifndef TRACER_OMP
!$OMP DO
#endif           
     DO k = 1, ke
        DO j = 1, je
!hh not needed due to following bounds_exch
!              t2o(1,j,k)   = t2o(1,j,k)   + tflux(1,j,k)
!              t2o(ie,j,k)  = t2o(ie,j,k)  - tflux(ie1,j,k)
           DO i = 2, ie1
              t2o(i,j,k)   = t2o(i,j,k)   + tflux(i,j,k) - tflux(i-1,j,k)
           END DO
        END DO
     END DO 

!$OMP SINGLE
     call bounds_exch('p',t2o,'octdiff_trf 3 ')
!$OMP END SINGLE


     !Y DIRECTION - tracer dependent values

#ifdef ISOPYK

#ifndef TRACER_OMP
!$OMP DO
#endif           
     DO j=1,je1
        DO k=1,ke
           DO i=1,ie
              IF(amsue(i,j,k).GT.0.)THEN
                 tv_diff=(trf(i,j+1,k)-trf(i,j,k))*dlyvi(i,j)
                 !
                 !       triangle left,upw
                 !
                 tfluy(i,j,k)=tfluy(i,j,k)+ycoeff_lo(i,j,k)*(tv_diff &
                      +sloploy(i,j,k)*                               &
                      (trf(i,j,ko(j,k))-trf(i,j,k))*di(k))
                 tfluz=zcoeff_loy(i,j,k)*(sloploy(i,j,k)*            &
                      (trf(i,j,ko(j,k))-trf(i,j,k))*                 &
                      dlxyp(i,j)*di(k)+1.*dlxyp(i,j)*tv_diff)

                 t2o(i,j,ko(j,k))=t2o(i,j,ko(j,k))-tfluz
                 t2o(i,j,k)=t2o(i,j,k)+tfluz
                 !
                 !       triangle left down
                 !
                 tfluy(i,j,k)=tfluy(i,j,k)+ycoeff_lu(i,j,k)*(tv_diff &
                      +slopluy(i,j,k)*                               &
                      (trf(i,j,k)-trf(i,j,ku(j,k)))*di(ku(j,k)))
                 tfluz=zcoeff_luy(i,j,k)*(slopluy(i,j,k)*            &
                      (trf(i,j,k)-trf(i,j,ku(j,k)))*                 &
                      dlxyp(i,j)*di(ku(j,k))+1.*dlxyp(i,j)*tv_diff)

                 t2o(i,j,k)=t2o(i,j,k)-tfluz
                 t2o(i,j,ku(j,k))=t2o(i,j,ku(j,k))+tfluz

              END IF
           END DO
        END DO
     END DO


#ifndef TRACER_OMP
!$OMP DO
#endif           
    DO j=1,je1
       DO k=1,ke
          DO i=1,ie
             IF(amsue(i,j,k).GT.0.)THEN
                 tv_diff=(trf(i,j+1,k)-trf(i,j,k))*dlyvi(i,j)
                 !
                 !       triangle right,upw
                 !
                 tfluy(i,j,k)=tfluy(i,j,k)+ycoeff_ro(i,j,k)*(tv_diff &
                      +sloproy(i,j,k)*                               &
                      (trf(i,j+1,ko(j,k))-trf(i,j+1,k))*di(k))
                 tfluz=zcoeff_roy(i,j,k)*(sloproy(i,j,k)*            &
                      (trf(i,j+1,ko(j,k))-trf(i,j+1,k))*             &
                      dlxyp1j(i,j)*di(k)+1.*dlxyp1j(i,j)*tv_diff)

                 t2o(i,j+1,ko(j,k))=t2o(i,j+1,ko(j,k))-tfluz
                 t2o(i,j+1,k)=t2o(i,j+1,k)+tfluz
                 !
                 !       triangle right down
                 !
                 tfluy(i,j,k)=tfluy(i,j,k)+ycoeff_ru(i,j,k)*(tv_diff &
                      +slopruy(i,j,k)*                               &
                      (trf(i,j+1,k)-trf(i,j+1,ku(j,k)))*di(ku(j,k)))
                 tfluz=zcoeff_ruy(i,j,k)*(slopruy(i,j,k)*            &
                      (trf(i,j+1,k)-trf(i,j+1,ku(j,k)))*             &
                      dlxyp1j(i,j)*di(ku(j,k))+1.*dlxyp1j(i,j)*tv_diff)

                 t2o(i,j+1,k)=t2o(i,j+1,k)-tfluz
                 t2o(i,j+1,ku(j,k))=t2o(i,j+1,ku(j,k))+tfluz

              END IF
           END DO
        END DO
     END DO

#else /*ISOPYK*/

#ifndef TRACER_OMP
!$OMP DO
#endif           
     DO j=1,je1
        DO k=1,ke
           DO i=1,ie                
              IF(amsue(i,j,k).GT.0.)THEN        
                 tv_diff=(trf(i,j+1,k)-trf(i,j,k))*dlyvi(i,j)
                 tfluy(i,j,k)=yflux(i,j,k)*tv_diff
              ENDIF
           END DO
        END DO
     END DO

#endif

!        SUMMATION
!
!$OMP SINGLE
     CALL bounds_exch('v',tfluy,'octdiff_trf 4')
!$OMP END SINGLE

#ifndef TRACER_OMP
! #slo# this loop cannot be parallelised in k - due to outerunroll ??
! !!CDIR OUTERUNROLL=20
!$OMP DO
#endif           
     DO k=1,ke
        DO j=2,je1
           DO i=1,ie 
              t2o(i,j,k)=t2o(i,j,k)+tfluy(i,j,k)-tfluy(i,j-1,k)
!              t2o(i,j+1,k)=t2o(i,j+1,k)-tfluy(i,j,k)
           END DO
        END DO
     END DO

!$OMP SINGLE
     CALL bounds_exch('p',t2o,'octdiff trf 5')
!$OMP END SINGLE

#ifndef TRACER_OMP
!$OMP DO
#endif           
     DO j=2,je1
        DO k=1,ke
           DO i=1,ie
              ten(i,j,k)=t2o(i,j,k)/vol_term(i,j,k)-trf(i,j,k)
           END DO
        END DO
     END DO

  ELSE ! ah00.gt.almzer

#ifndef TRACER_OMP
!$OMP WORKSHARE
#endif           
    ten(:,:,:)=0.
#ifndef TRACER_OMP
!$OMP END WORKSHARE
#endif

  END IF ! ah00.gt.almzer

  ! VERTICAL DIFFUSION -tracer dependent

#ifndef TRACER_OMP
!$OMP DO
#endif           
   DO k=1,ke
     DO j=1,je 
       DO i=1,ie 
         t2o(i,j,k)=trf(i,j,k)
       END DO
     END DO
   END DO

  !     IMPLICIT VERTICAL DIFFUSION

#ifndef TRACER_OMP
!$OMP DO
#endif           
  DO j=2,je-1 
    DO k=2,ke
      DO i=1,ie 
        t2o(i,j,k)=t2o(i,j,k)-tridsy(i,j,k-1,1)*t2o(i,j,k-1)
      END DO
    END DO  
  END DO

#ifndef TRACER_OMP
!$OMP DO
#endif           
  DO j = 2,je-1
    DO i = 1,ie
      trf(i,j,ke)=t2o(i,j,ke)/tridsy(i,j,ke,2)
    END DO
  END DO

#ifndef TRACER_OMP
!$OMP DO
#endif           
  DO j = 2, je-1
    DO k = ke-1, 1, -1      
      DO i = 1,ie 
        trf(i,j,k)=(t2o(i,j,k)-tridsy(i,j,k,3)*trf(i,j,k+1))              &
                  /tridsy(i,j,k,2)
      END DO
    END DO
  END DO


!$OMP SINGLE
  CALL bounds_exch('p',trf,'octdiff_trf 6')
!$OMP END SINGLE


#ifndef TRACER_OMP
!$OMP DO
#endif           
  DO j = 2,je1
    DO k = 1,ke
      DO i = 1,ie
        trf(i,j,k)=trf(i,j,k)+ten(i,j,k)
      END DO
    END DO
  END DO

!$OMP SINGLE
  CALL bounds_exch('p',trf,'octdiff_trf 7')
!$OMP END SINGLE

  !
  ! BIHARMONIC tracer DIFFUSION - tracer dependent
  !
  IF (aulapts.GT.almzer) THEN

#ifndef TRACER_OMP
!$OMP DO
#endif           
    DO k=1,ke
      DO j=2,je1
        DO i=2,ie1
          t2o(i,j,k)=weto(i,j,k)*(                                           &
              (weto(i-1,j,k)*(trf(i-1,j,k)-trf(i,j,k))*dlxui(i-1,j)          &
              +weto(i+1,j,k)*(trf(i+1,j,k)-trf(i,j,k))*dlxui(i,j))/dlxp(i,j) &
              +(weto(i,j-1,k)*(trf(i,j-1,k)-trf(i,j,k))*dlyvi(i,j-1)         &
              +weto(i,j+1,k)*(trf(i,j+1,k)-trf(i,j,k))*dlyvi(i,j))/dlyp(i,j))
        END DO
      END DO
    END DO

!$OMP SINGLE
    CALL bounds_exch('p',t2o,'octdiff_trf 8')
!$OMP END SINGLE

#ifndef TRACER_OMP
!$OMP DO
#endif           
    DO k=1,ke
      zzsur0=0.
      if (k == 1) zzsur0=1.
      DO j=2,je1
        DO i=2,ie1
          trf(i,j,k)=trf(i,j,k)                                              &
               - (weto(i,j,k)/((ddpo(i,j,k)+almzer                           &
               +zzsur0*(zo(i,j)-rhoicwa*sictho(i,j)-rhosnwa*sicsno(i,j)))    &
               *dlxyp(i,j)))*dzw(k)*(                                        &
               weto(i-1,j,k)*aulx(i-1,j)*dlyu(i-1,j)*                        &
               (t2o(i-1,j,k)-t2o(i,j,k))*dlxui(i-1,j)                        &
               +weto(i+1,j,k)*aulx(i,j)*dlyu(i,j)*                           &
               (t2o(i+1,j,k)-t2o(i,j,k))*dlxui(i,j)                          &
               +weto(i,j-1,k)*auly(i,j-1)*dlxv(i,j-1)*                       &
               (t2o(i,j-1,k)-t2o(i,j,k))*dlyvi(i,j-1)                        &
               +weto(i,j+1,k)*auly(i,j)*dlxv(i,j)*                           &
               (t2o(i,j+1,k)-t2o(i,j,k))*dlyvi(i,j))
        END DO
      END DO
    END DO

!$OMP SINGLE
    CALL bounds_exch('p',trf,'octdiff_trf 9')
!$OMP END SINGLE

  ENDIF ! aulapts > almzer

#ifndef TRACER_OMP
!$OMP END PARALLEL
#endif           

END SUBROUTINE octdiff_trf

