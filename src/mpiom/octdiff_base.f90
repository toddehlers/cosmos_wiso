SUBROUTINE octdiff_base

  USE mo_kind, ONLY: dp
  USE mo_param1
  USE mo_parallel
  USE mo_commo1
  USE mo_commoau1
  USE mo_units
  USE mo_octdiff

  !     COMPUTES DIFFUSION OF TEMPERATURE THO AND SALINITY SAO
  !
  !     UWE MIKOLAJEWICZ 12/99
  !
  !     VERTICAL DIFFUSION IMPLICITE
  !     HORIZONTAL DIFFUSION BOTH HARMONIC AND BIHARMONIC, BOTH EXPLICITE
  !                                AH00         AULAPTS
  !     ACTUAL COEFFICIENTS SPATIALLY VARIABLE WITH RESOLUTION
  !                               AHX,AHY      AULX,AULY
  !
  !     VARIABLE THICKNESS OF SURFACE LAYER NOW INCLUDED
  !
  !     Changes R. Johanni, 2003-11-10:
  !     DSLOPX, DSLOPY, DVSLOP were calculated but nowhere used -> removed
  !
  !     Changes R. Smith, 2004-09-27
  !     tracer-independent matrices calculated separately first for use
  !     with multiple tracers

  IMPLICIT NONE

  REAL(dp) :: ahx(ie,je), ahy(ie,je), ahxd(ie,je), ahyd(ie,je), zzsurf(ke)

  INTEGER :: i, j, k
  REAL(dp) :: otmp
  REAL(dp) :: scale, slcut, stabmin, cutdel
  REAL(dp) ::rho_diff, xterm, zterm, yterm
  !
  DO k = 1, ke
    zzsurf(k) = 0.0_dp
  ENDDO
  zzsurf(1) = 1.0_dp
  !
  !     MINIMUM FOR STABILITY
  !
  stabmin = 1.e-7_dp
  cutdel  = 0.02_dp
  !
!$OMP PARALLEL PRIVATE(i,j,k,scale,slcut,rho_diff,xterm,zterm,yterm,otmp)
  !
  IF(ah00 > almzer)THEN
    !
!$OMP DO
    DO j = 1, je
      DO i = 1, ie
        ahx(i,j)  = ah00*MAX(dlxu(i,j),dlyu(i,j))
        ahy(i,j)  = ah00*MAX(dlxv(i,j),dlyv(i,j))      !emr vpoints
        ahxd(i,j) = ahx(i,j)*dt*dlyu(i,j)
        ahyd(i,j) = ahy(i,j)*dt*dlxv(i,j)
      ENDDO
    ENDDO
    
!!$!$OMP SINGLE
!!$  call bounds_exch('p',zo)
!!$  call bounds_exch('p',sictho)
!!$  call bounds_exch('p',sicsno)
!!$  call bounds_exch('p',weto)
!!$!$OMP END SINGLE

!$OMP DO
    DO k = 1, ke                  
      DO j = 1, je
        DO i = 1, ie
          vol_term(i,j,k) = dlxp(i,j)*dlyp(i,j)*(zzsurf(k)*(zo(i,j) &
               -rhoicwa*sictho(i,j)-rhosnwa*sicsno(i,j))            &
               +ddpo(i,j,k)+(1.0_dp-weto(i,j,k)))
        ENDDO
      ENDDO
    ENDDO

!$OMP DO
    DO k = 1, ke                  
      DO j = 1, je 
        ko(j,k) = MAX(k-1,1)
        ku(j,k) = MIN(k+1,ke)
      ENDDO
    ENDDO
!$OMP DO
    DO k=1,ke                  
      DO j=1,je 
        zsurf(j,k) = REAL(k-ko(j,k),dp)
        zbott(j,k) = REAL(ku(j,k)-k,dp)
      ENDDO
    ENDDO
    !
    !X DIRECTION- tracer independent values
!$OMP DO
    DO j = 2, je-1
      DO k = 1, ke 
        DO i = 1, ie1
          IF (amsuo(i,j,k) > 0.0_dp) THEN
#ifdef ISOPYK
            slcut = cutdel*dzw(k)**2/(ahx(i,j)*dt)
            rho_diff = (rhoo(i+1,j,k)-rhoo(i,j,k))/dlxu(i,j)
            xterm = 0.25_dp*ahxd(i,j)*dduo(i,j,k)
            zterm = 0.25_dp*dt*ahx(i,j)
            !
            !       TRIANGLE LEFT,UPW
            !       
            sloplox(i,j,k) = rho_diff/MAX(stabio(i,j,k),stabmin)
            scale = 1.0_dp/MAX(1.0_dp,(sloplox(i,j,k)/slcut)**2)
            xcoeff_lo(i,j,k)  = xterm*zsurf(j,k)*scale
            zcoeff_lox(i,j,k) = zterm*zsurf(j,k)*scale*sloplox(i,j,k)
            !       
            !       TRIANGLE LEFT DOWN
            !        
            sloplux(i,j,k) = rho_diff/MAX(stabio(i,j,ku(j,k)),stabmin)
            scale = 1.0_dp/MAX(1.0_dp,(sloplux(i,j,k)/slcut)**2)
            xcoeff_lu(i,j,k)  = xterm*zbott(j,k)*weto(i,j,ku(j,k))*scale
            zcoeff_lux(i,j,k) = zterm*zbott(j,k)*scale*sloplux(i,j,k)*      &
                 weto(i,j,ku(j,k))
            !       
            !       TRIANGLE RIGHT,UPW
            !       
            sloprox(i,j,k) = rho_diff/MAX(stabio(i+1,j,k),stabmin)
            scale = 1.0_dp/MAX(1.0_dp,(sloprox(i,j,k)/slcut)**2)
            xcoeff_ro(i,j,k)  = xterm*zsurf(j,k)*scale
            zcoeff_rox(i,j,k) = zterm*zsurf(j,k)*scale*sloprox(i,j,k)
            !
            !       TRIANGLE RIGHT DOWN
            !
            sloprux(i,j,k) = rho_diff/MAX(stabio(i+1,j,ku(j,k)),stabmin)
            scale = 1.0_dp/MAX(1.0_dp,(sloprux(i,j,k)/slcut)**2)
            xcoeff_ru(i,j,k)  = xterm*zbott(j,k)*weto(i+1,j,ku(j,k))*scale
            zcoeff_rux(i,j,k) = zterm*zbott(j,k)*scale*sloprux(i,j,k)*     &
                 weto(i+1,j,ku(j,k))
#else
            xflux(i,j,k) = ahxd(i,j)*dduo(i,j,k)
#endif
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    !Y DIRECTION - Tracer independent values
    !
!$OMP DO
    DO j = 1, je1
      DO k = 1, ke
        DO i = 2, ie1
          IF (amsue(i,j,k) > 0.0_dp) THEN
#ifdef ISOPYK
            slcut = cutdel*dzw(k)**2/(ahy(i,j)*dt)
            rho_diff = (rhoo(i,j+1,k)-rhoo(i,j,k))/dlyv(i,j)
            yterm = 0.25_dp*ahyd(i,j)*ddue(i,j,k)
            zterm = 0.25_dp*dt*ahy(i,j)
            !
            !       TRIANGLE LEFT,UPW
            !
            sloploy(i,j,k) = rho_diff/MAX(stabio(i,j,k),stabmin)
            scale = 1.0_dp/MAX(1.0_dp,(sloploy(i,j,k)/slcut)**2)
            ycoeff_lo(i,j,k)  = yterm*zsurf(j,k)*scale
            zcoeff_loy(i,j,k) = zterm*zsurf(j,k)*scale*sloploy(i,j,k)
            !
            !       TRIANGLE LEFT DOWN
            !
            slopluy(i,j,k) = rho_diff/MAX(stabio(i,j,ku(j,k)),stabmin)
            scale = 1.0_dp/MAX(1.0_dp,(slopluy(i,j,k)/slcut)**2)
            ycoeff_lu(i,j,k)  = yterm*zbott(j,k)*weto(i,j,ku(j,k))*scale
            zcoeff_luy(i,j,k) = zterm*zbott(j,k)*scale*slopluy(i,j,k)*     &
                 weto(i,j,ku(j,k))
            !
            !       TRIANGLE RIGHT,UPW
            !
            sloproy(i,j,k) = rho_diff/MAX(stabio(i,j+1,k),stabmin)
            scale = 1.0_dp/MAX(1.0_dp,(sloproy(i,j,k)/slcut)**2)
            ycoeff_ro(i,j,k)  = yterm*zsurf(j,k)*scale
            zcoeff_roy(i,j,k) = zterm*zsurf(j,k)*scale*sloproy(i,j,k)

            !
            !       TRIANGLE RIGHT DOWN
            !
            slopruy(i,j,k) = rho_diff/MAX(stabio(i,j+1,ku(j,k)),stabmin)
            scale = 1.0_dp/MAX(1.0_dp,(slopruy(i,j,k)/slcut)**2)
            ycoeff_ru(i,j,k)  = yterm*zbott(j,k)*weto(i,j+1,ku(j,k))*scale
            zcoeff_ruy(i,j,k) = zterm*zbott(j,k)*scale*slopruy(i,j,k)*     &
                 weto(i,j+1,ku(j,k))
#else
            yflux(i,j,k) = ahyd(i,j)*ddue(i,j,k)
#endif
          ENDIF
        ENDDO
      ENDDO
    ENDDO

!#ifndef ISOPYK    
!!$!$OMP SINGLE
!  call bounds_exch('u+',xflux)
!  call bounds_exch('v+',yflux)
!!$!$OMP END SINGLE
!#endif

  ENDIF ! ah00.gt.almzer
  !
  ! VERTICAL DIFFUSION - tracer independent values
  !     IMPLICIT VERTICAL DIFFUSION
  !
!CDIR NOLOOPCHG
  DO k = 1,ke
    ! INCLUDE ACTUAL LEVEL THICKNESS
!$OMP DO
    DO j = 2,je1
      tridsy(:,j,k,1) = - dt*dvo(:,j,k)*weto(:,j,k)*di(k)         &
           /(ddpo(:,j,k)+zzsurf(k)*(zo(:,j)-rhoicwa*sictho(:,j)   &
           -rhosnwa*sicsno(:,j))+almzer)
      tridsy(:,j,k,3) = - dt*dvo(:,j,k+1) * di(k+1)               &
           /(ddpo(:,j,k)+zzsurf(k)*(zo(:,j)-rhoicwa*sictho(:,j)   &
           -rhosnwa*sicsno(:,j))+almzer)
      tridsy(:,j,k,2) = 1.0_dp - tridsy(:,j,k,1) - tridsy(:,j,k,3)
    END DO
  ENDDO
  !

!!$!CDIR NOLOOPCHG
!!$  DO k = 2,ke
!!$!$OMP DO
!!$    DO j = 2,je1
!!$      tridsy(:,j,k-1,1) = tridsy(:,j,k,1) / tridsy(:,j,k-1,2)
!!$      tridsy(:,j,k,2)   = tridsy(:,j,k,2) - tridsy(:,j,k-1,3)     &
!!$                         *tridsy(:,j,k,1) / tridsy(:,j,k-1,2)
!!$    END DO
!!$  ENDDO


!HH This workaround is suggested by NEC
!   it is needed to avoid different results
!   for openmp and non-openmp optimisation


!CDIR NOLOOPCHG
   DO k = 2,ke
!$OMP DO
     DO j = 2,je-1
       DO i=1,ie
       otmp = 1.0/tridsy(i,j,k-1,2)
       tridsy(i,j,k-1,1) = tridsy(i,j,k,1) * otmp
       tridsy(i,j,k,2)   = tridsy(i,j,k,2) - tridsy(i,j,k-1,3)     &
                         * tridsy(i,j,k,1) * otmp
       END DO
     END DO
   END DO




  !
  ! BIHARMONIC DIFFUSION - tracer independent values
  !
  IF (aulapts > almzer) THEN
!$OMP DO
#ifndef AULREDSC
!CDIR NOUNROLL
    DO j = 1, je
      DO i = 1, ie
        aulx(i,j) = aulapts*dlxu(i,j)**4
        auly(i,j) = aulapts*dlyv(i,j)**4
      ENDDO
    ENDDO
#else
!CDIR NOUNROLL
    DO j = 1, je
      DO i = 1, ie
        aulx(i,j) = 1.e5_dp*aulapts*dlxu(i,j)**3
        auly(i,j) = 1.e5_dp*aulapts*dlyv(i,j)**3
      ENDDO
    ENDDO
#endif
  ENDIF
  !
!$OMP END PARALLEL
  !
END SUBROUTINE octdiff_base












