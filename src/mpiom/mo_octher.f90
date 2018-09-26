 MODULE mo_octher

  USE mo_param1
  USE mo_param3      ! This is not actually used here but only in subroutine adisitj which is called from octher 
                     ! and is supposed to be expanded inline. The purpose of this dependency is to force compilation
                     ! of mo_param3.f90 before compilation of mo_octher.f90, because mo_param3.mod is needed for
                     ! the inline expansion of adisitj.
  USE mo_mpi
  USE mo_parallel
  USE mo_commo1
  USE mo_commoau1
  USE mo_commoau2
  USE mo_units



#ifdef PBGC
  USE mo_tro ,only: dilcor_gtrf2 , dilcor_ptrf2
  USE mo_param1_bgc, only: nocetra
  USE mo_carbch, only: ocetra
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef ADDCONTRA
  USE mo_tro, only: dilcor_gtrf2, dilcor_ptrf2
  USE mo_param1_add, only: nocectra
  USE mo_contra, only: ocectra
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

 CONTAINS

  SUBROUTINE octher

  USE mo_diagnosis, only: calc_potential_energy_release 

  integer :: l
  
  !-----------------------------------------------------------------------
  !
  !     sbr octher computes
  !
  !          relax_surf     -  boundary forcing on salt, temperature and zeta
  !          river_runoff -  only in the uncoupled model
  !          convection   -  baroclinic pressure in each layer +  convective adjustment
  !          calc_rinum   -  richardson-number depending coefficients for
  !                          vertical diffusion of momentum (avo) and
  !                          temperature and salinity (dvo)
  !          calc_den     -  update of density field 
  !          
  !
  !
  !-----------------------------------------------------------------------
  
#ifndef __coupled

#ifdef PBGC
    call dilcor_gtrf2    
#endif

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
#ifdef ADDCONTRA
    call dilcor_gtrf2    
#endif
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

    CALL relax_surf      ! only in the uncoupled model     

    CALL river_runoff

#ifdef PBGC
    do l=1,nocetra
       call dilcor_ptrf2(ocetra(:,:,1,l))
    enddo
#endif

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
#ifdef ADDCONTRA
    do l=1,nocectra
       call dilcor_ptrf2(ocectra(:,:,1,l))
    enddo
#endif
!------------------------------------------------------------------------
!------------------------------------------------------------------------

#endif

#ifndef NURDIF
    if (LCONVDIAG) call calc_potential_energy_release(1)
#endif
    CALL convection
#ifndef NURDIF
    if (LCONVDIAG) call calc_potential_energy_release(2)
#endif
    CALL calc_rinum

    CALL calc_dens

  END SUBROUTINE octher


  SUBROUTINE relax_surf
        
    INTEGER    i,j
    REAL       reiscon,zsurf(ke),oldso

!     boundary forcing on temperature, salt and zeta
!     relaxation time on salinity   : 1./ relsal

!$omp parallel private (i,j,reiscon,oldso)
!$omp do        

    DO j=1,je
       DO i=1,ie
          
          eminpo(i,j)=0.
          
          reiscon=0.
          IF(sicomo(i,j).LE.0.01) THEN
             reiscon=1.
          ENDIF
          
#ifdef EISREST
          reiscon=1.-sicomo(i,j)
#endif
          oldso=sao(i,j,1)
          sao(i,j,1)=sao(i,j,1)+dt*relsal*reiscon*(relsao(i,j)-sao(i,j,1))
          tho(i,j,1)=tho(i,j,1)+dt*reltem*(reltho(i,j)-tho(i,j,1))

#ifdef ANOMALY_FORCING
          reltem=3.0e-6
          IF (ibek(i,j) .EQ. 1) THEN
             tho(i,j,1)=tho(i,j,1)+dt*reltem*(reltho(i,j)-3.0-tho(i,j,1))
          ENDIF
#endif
          eminpo(i,j)=(ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa             &
                     -sicsno(i,j)*rhosnwa)                                 &
                     *(MAX(oldso,1.e-3)/MAX(sao(i,j,1),1.e-3)-1.)
       ENDDO
    ENDDO
!$omp end do

    dti=1./dt

!$omp do  
    DO j=1,je
       DO i=1,ie
          zo(i,j)=(zo(i,j)+eminpo(i,j))*weto(i,j,1)
          eminpo(i,j)=eminpo(i,j)*dti                       !uwe eminpo in m/s
       ENDDO
    ENDDO
!$omp end do
!$omp end parallel

  END SUBROUTINE relax_surf



  SUBROUTINE river_runoff_ini

#ifndef __coupled
!  initialisation of river input locations

#ifndef RIVER_GIRIV
      INTEGER :: i,n,m
      REAL :: dist, pibog

      pibog=180./(2.*asin(1.))

!hh-------------------------------------------------------------------      
!hh   river runoff data from ieee extra file: observations 12mon.clim.

      IF(p_pe==p_io) THEN
         OPEN(io_in_rval,file='runoff_obs',form='unformatted')
         OPEN(io_in_rpos,file='runoff_pos',form='unformatted')
         DO i=1,numriv
            READ(io_in_rval)ibla
            READ(io_in_rval)(rival(i,m),m=1,12)
            READ(io_in_rpos)ibla
            READ(io_in_rpos)rivlat(i),rivlon(i)
         ENDDO
         CLOSE(io_in_rval)
         CLOSE(io_in_rpos)
      ENDIF
      CALL p_bcast(rival,p_io)
      CALL p_bcast(rivlat,p_io)
      CALL p_bcast(rivlon,p_io)


      DO n=1,numriv
      CALL suchij(rivlat(n),rivlon(n),1,irivi(n),irivj(n),dist,1.)
      WRITE(io_stdout,*)'river nr: ',rivlat(n),rivlon(n),n              &
     &                       ,irivi(n),irivj(n)                         &
     &                       ,dist,                                     &
     &                       weto_g(irivi(n),irivj(n),1),               &
     &                       gila_g(2*irivi(n),2*irivj(n))*pibog,       &
     &                       giph_g(2*irivi(n),2*irivj(n))*pibog
      ENDDO

#endif /*RIVER_GIRIV*/


#ifdef GLACCALV
      IF(p_pe==p_io) THEN
         OPEN(io_in_glac,file='gletscher_5653',form='unformatted')
         DO i=1,numglac
          READ(io_in_glac)glaclat(i),glaclon(i),glacval(i)
         ENDDO
         CLOSE (io_in_glac)
      ENDIF
      CALL p_bcast(glaclat,p_io)
      CALL p_bcast(glaclon,p_io)
      CALL p_bcast(glacval,p_io)

      DO n=1,numglac
      CALL suchij(glaclat(n),glaclon(n),1,iglac(n),jglac(n),dist,1.)

      WRITE(io_stdout,*)'glacier nr: ',glaclat(n),glaclon(n),n          &
     &                       ,iglac(n),jglac(n)                         &
     &                       ,dist,                                     &
     &                       weto_g(iglac(n),jglac(n),1),               &
     &                       gila_g(2*iglac(n),2*jglac(n))*pibog,       &
     &                       giph_g(2*iglac(n),2*jglac(n))*pibog
      ENDDO
#endif /*GLACCALV*/

#endif /*__coupled*/

  END SUBROUTINE river_runoff_ini


  SUBROUTINE river_runoff

    INTEGER i,j,n,monmon
    REAL driv,zzzdz,awert,ewert


    rivrun(:,:)=0.                   !  runoff diagnostic

#ifndef RIVER_GIRIV
!
!  UPDATE FORCING ONCE PER DAY

       IF(LDAYS.LE.(MONLEN(LMONTS)+1)/2)THEN
        MONMON=LMONTS-1
        IF(MONMON.LT.1)MONMON=12
        AWERT=0.5+(FLOAT(LDAYS)-0.5)/FLOAT(MONLEN(LMONTS))
       ELSE
        MONMON=LMONTS+1
        IF(MONMON.GT.1)MONMON=1
        AWERT=0.5+(FLOAT(MONLEN(LMONTS)-LDAYS)+0.5)                     &
     &      /FLOAT(MONLEN(LMONTS))
       ENDIF
        EWERT=1.-AWERT

       DO N=1,NUMRIV
         FRIV(N)=(AWERT*RIVAL(N,LMONTS)+EWERT*RIVAL(N,MONMON))
       ENDDO

!$omp parallel private(i,j,n,zzzdz)
!$omp do
    DO n=1,numriv
       i=irivi(n)-p_ioff
       j=irivj(n)-p_joff
       IF(i>=1 .AND. i<=ie .AND. j>=1 .AND. j<=je) THEN
          zzzdz=MAX(almzer,ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa      &
               -sicsno(i,j)*rhosnwa)
!         if(weto(i,j,1).lt.0.5)write(io_stdout,*)'alarm! river ',n,i,j
          ddriv(n)=friv(n)*dt/(dlxp(i,j)*dlyp(i,j))
!uwe      use actual layerthickness for mass/salt conservation
          sao(i,j,1)=sao(i,j,1)*zzzdz/(zzzdz+ddriv(n))
          zo(i,j)=zo(i,j)+ddriv(n)
          prech(i,j)=prech(i,j)+ddriv(n)/dt
          preco(i,j)=preco(i,j)+ddriv(n)/dt
          rivrun(i,j)=ddriv(n)/dt
!        write(io_stdout,*)' river ',n,ddriv(n),zo(i,j),sao(i,j,1)
       ENDIF
    ENDDO
!$omp end do
!$omp end parallel
#else

!$omp parallel private(i,j,driv,zzzdz)
!$omp do
    DO j=1,je
       DO i=1,ie   
          zzzdz=MAX(almzer,ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa        &
               -sicsno(i,j)*rhosnwa)
          !       if(weto(i,j,1).lt.0.5.and.giriv(i,j).gt.0.)then
          !           write(io_stdout,*)'alarm! river ',n,i,j
          !       endif
          driv=giriv(i,j)*dt/(dlxp(i,j)*dlyp(i,j))
          !uwe      use actual layerthickness for mass/salt conservation
          sao(i,j,1)=sao(i,j,1)*zzzdz/(zzzdz+driv)
          zo(i,j)=zo(i,j)+driv
          preco(i,j)=preco(i,j)+driv/dt
          prech(i,j)=prech(i,j)+driv/dt
          rivrun(i,j)=driv/dt
       ENDDO
    ENDDO
!$omp end do
!$omp end parallel

#endif/*RIVER_GIRIV*/


#ifdef GLACCALV

!     glacier calving
!
      tempfac=tfreez-entmel/rocp

!     assumption glacier calving as very cold water to conserve heat

      DO n=1,numglac
       i=iglac(n)-p_ioff
       j=jglac(n)-p_joff
       IF(i>=1 .AND. i<=ie .AND. j>=1 .AND. j<=je) THEN
         wert=glacval(n)/(dlxp(i,j)*dlyp(i,j))

!        ice compactness minimum glacinput*0.5

         thick=ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa                  &
     &       +sicsno(i,j)*rhosnwa
         sao(i,j,1)=sao(i,j,1)*thick/(thick+wert*dt)
         tho(i,j,1)=(tho(i,j,1)*thick+wert*dt*tempfac)                  &
     &              /(thick+wert*dt)
         zo(i,j)=zo(i,j)+wert*dt
       ENDIF
      ENDDO

!      update p-e field for diagnostics
       DO n=1,numglac
       i=iglac(n)-p_ioff
       j=jglac(n)-p_joff
       IF(i>=1 .AND. i<=ie .AND. j>=1 .AND. j<=je) THEN
         preco(i,j)=preco(i,j)+glacval(n)/(dlxp(i,j)*dlyp(i,j))      
         prech(i,j)=prech(i,j)+glacval(n)/(dlxp(i,j)*dlyp(i,j))      
       ENDIF
       ENDDO
#endif



  END SUBROUTINE river_runoff

  SUBROUTINE convection

    USE mo_mean

    INTEGER i,j,k
    
    REAL rhuppo(ie,je)
    REAL zsurf(ke)
    REAL stabio1,disti,dtts
    REAL tupper,tlower,supper,slower,ddhelp,sq,tq,helpswi,switc


 
!$omp parallel private(i,j,k, disti,stabio1,tq,sq,switc,    &
!$omp                  helpswi,tupper, tlower,                    &
!$omp                  supper, slower, ddhelp)

    !
    !=====================================================================
    !
    !     b)


    !
    !     baroclinic pressure and stability
    !
    !
    !---------------------------------------------------------------------
    !
    !     b.1) upper layer
    !

    zsurf(1)=1.
    DO k=2,ke
       zsurf(k)=0.
    ENDDO

!$omp do
    DO j=1,je

       DO i=1,ie
          vk1e(i,j,1)=0.
          shelp(i,j)=sao(i,j,1)
          thelp(i,j)=tho(i,j,1)
       ENDDO
       
       CALL adisitj(thelp,shelp,preff(1),j)
       CALL rho1j(thelp,shelp,preff(1),rhelp,j)
       
       DO i=1,ie
          stabio(i,j,1) = 0.
          po(i,j,1)     = g*tiestu(1)*0.00097*rhelp(i,j)
          s1o(i,j,1)    = rhelp(i,j)
       ENDDO
       
       DO k=2,ke
          
          disti=1./dz(k)
          DO i=1,ie
             shelp(i,j)=sao(i,j,k)
             thelp(i,j)=tho(i,j,k)
          ENDDO
          
          CALL adisitj(thelp,shelp,preff(k),j)
          CALL rho1j(thelp,shelp,preff(k),rhelp,j)
          
          DO i=1,ie
             shelp(i,j)=sao(i,j,k-1)
             thelp(i,j)=tho(i,j,k-1)
          ENDDO
          
          CALL adisitj(thelp,shelp,preff(k),j)
          CALL rho1j(thelp,shelp,preff(k),rhuppo,j)
          
          DO i=1,ie
             s1o(i,j,k)=rhelp(i,j)
             stabio1 = disti * ( rhelp(i,j) - rhuppo(i,j) )
             po(i,j,k) = po(i,j,k-1) + g*dz(k)*0.00049*(rhelp(i,j)+rhuppo(i,j))
             !
#ifndef PLUME
#ifndef UMKLAP
             !uwe     salt conservation 12/99
             tq=((ddpo(i,j,k-1)+zsurf(k-1)*(zo(i,j)-sictho(i,j)*rhoicwa          &
                  -sicsno(i,j)*rhosnwa))        &
                  *tho(i,j,k-1)+ddpo(i,j,k)*tho(i,j,k))                         &
                  /(ddpo(i,j,k)+ddpo(i,j,k-1)                                   &
                  +zsurf(k-1)*(zo(i,j)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa) &
                  +(1.-weto(i,j,k)))
             sq=((ddpo(i,j,k-1)+zsurf(k-1)*(zo(i,j)-sictho(i,j)*rhoicwa          &
                  -sicsno(i,j)*rhosnwa))        &
                  *sao(i,j,k-1)+ddpo(i,j,k)*sao(i,j,k))                         &
                  /(ddpo(i,j,k)+ddpo(i,j,k-1)                                   &
                  +zsurf(k-1)*(zo(i,j)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa) &
                  +(1.-weto(i,j,k)))
#endif /*UMKLAP*/
             !
             !lk not used
             !lk        switc=(half-sign(half,stabio(i,j,k)))*weto(i,j,k)
             switc=MAX(0.,-stabio1/(1.e-11+ABS(stabio1)))*weto(i,j,k)
             !-et      converei=converei+switc
             !
#ifdef NURDIF
             helpswi=switc
             switc=0.
#endif /*NURDIF*/
             !
             vk1e(i,j,k)=switc
             !
#ifndef UMKLAP
             tho(i,j,k-1) = tq * switc + (one-switc) * tho(i,j,k-1)
             tho(i,j,k)   = tq * switc + (one-switc) * tho(i,j,k)
             sao(i,j,k-1) = sq * switc + (one-switc) * sao(i,j,k-1)
             sao(i,j,k)   = sq * switc + (one-switc) * sao(i,j,k)
#else /*UMKLAP*/
             IF (switc.GE.0.5) THEN
                tupper=tho(i,j,k-1)
                tlower=tho(i,j,k)
                supper=sao(i,j,k-1)
                slower=sao(i,j,k)
                ddhelp=ddpo(i,j,k-1)+zsurf(k-1)*(zo(i,j)-sictho(i,j)*rhoicwa &
                     -sicsno(i,j)*rhosnwa)
                !
                IF(ddpo(i,j,k).GT.ddhelp) THEN
                   tho(i,j,k-1)=tlower
                   sao(i,j,k-1)=slower
                   tho(i,j,k)=tlower+(tupper-tlower)*(ddhelp/ddpo(i,j,k))
                   sao(i,j,k)=slower+(supper-slower)*(ddhelp/ddpo(i,j,k))
                ELSE
                   tho(i,j,k)=tupper
                   sao(i,j,k)=supper
                   tho(i,j,k-1)=tupper+(tlower-tupper)*(ddpo(i,j,k)/ddhelp)
                   sao(i,j,k-1)=supper+(slower-supper)*(ddpo(i,j,k)/ddhelp)
                ENDIF
                !
             ENDIF
#endif /*UMKLAP*/
             !
             stabio1=(one-switc)*stabio1
             stabio(i,j,k)=MAX(stabio1,0.)
             !
#ifdef NURDIF
             switc=helpswi
#endif /*NURDIF*/
             !
             IF(kcondep(i,j).EQ.k-1) kcondep(i,j) = kcondep(i,j)+NINT(switc)
#endif /*PLUME*/
          ENDDO
       ENDDO
       !
    ENDDO ! j-loop
!$omp end do
    

#ifdef PLUME
!$omp single
!sjm plume convection
    dtts=dt
    CALL nlopps(tho,sao,dzw,ddpo,preff,kcondep,dtts)
!$omp end single
#endif /*PLUME*/
 
!$omp end parallel
    
  END SUBROUTINE convection

  SUBROUTINE calc_rinum

    USE mo_mean

    INTEGER i,j,k,ko,ku
    REAL rinumo,stabeps,switc,avo1,dvo1,dudo,gor,hho
    REAL topnah,wpendep,wtdecay,relne,relax,dudz,drdz0,cra,crd,zsurf(ke)

!$omp parallel private(i,j,k, rinumo,ku,ko,stabeps,wtdecay,dudo,hho    &
!$omp                 ,switc,avo1,dvo1)

    !======================================================================
    !
    !     d)
    !
    !                   calculation of richardson number dependent
    !                      vertical eddy viscosity   av     and
    !                      vertical eddy diffusivity dv
    !
    !
    !         rinum : ( (g/rho)*d(rho)/dz ) / ( (d(velocity)/dz)**2 )
    !                  richardson number (even,odd ==> rinume,rinumo)
    !         gor   : g/rho
    !
    !         av0   : numerical value of vertical eddy viscosity in case
    !                 of neutral stability, i.e. free turbulence
    !
    !---------------------------------------------------------------------
    zsurf(1)=1.
    DO k=2,ke
       zsurf(k)=0.
    ENDDO

    !
    relne=0.4
    relax=1.-relne
    !
    dudz=1.e4
    drdz0=1.e-3
    !
    gor=g/1025.
    !
    !--------------------------------------------------------------------
    !
    !     c.1)
    !
    !     vertical eddy viscosity  (for momentum equation)
    !
    !--------------------------------------------------------------------
    !
    !     d.1)    mixed-layer turbulence
    !
    !  amplitude of turbulence decays by factor wtdecay every model level.
    !  turbulence stops once density difference reaches the equivalent of
    !   wtdt
    !  temperature difference.
    !  turbulence under ice is / is not  enhanced.
    !                             ==
    cra=5.
    crd=5.

    wpendep=40.
    relne=0.4
    relax=1.-relne
    wtdecay=EXP(-dzw(1)/wpendep)

!$omp do
    DO j=1,je
       DO i=1,ie
          !hh       t1o(i,j,1)=wt*2.*(1.-sicomo(i,j))*fu10(i,j)**3*wtdecay
          !hh       s1o(i,j,1)=wt*(1.-sicomo(i,j))*fu10(i,j)**3*wtdecay
#ifdef REDWMICE
          t1o(i,j,1)=wa*(1.-sicomo(i,j))**2*fu10(i,j)**3
          s1o(i,j,1)=wt*(1.-sicomo(i,j))**2*fu10(i,j)**3
#else
          t1o(i,j,1)=wa*(1.-sicomo(i,j))*fu10(i,j)**3
          s1o(i,j,1)=wt*(1.-sicomo(i,j))*fu10(i,j)**3
#endif /*REDWMICE*/

       ENDDO
    ENDDO
!$omp end do
    
    if (iMEAN.ne.0)then
       if (LDIFFDIAG) then
          wtmix(:,:,1)=s1o(:,:,1)
          rinu(:,:,1)=0.
       endif
    endif


!$omp do
    DO j=2,je1

       DO k=2,ke
          ku=MIN(k+1,ke)
          ko=MAX(k-1,1)
          stabeps=cstabeps/dzw(ko)
          wtdecay=EXP(-dzw(ko)/wpendep)

          topnah=0.0
          DO i=2,ie1

             t1o(i,j,1)=t1o(i,j,1)*wtdecay*stabeps                            &
                  /(stabeps+0.5*(stabio(i,j,k)+(1.-zsurf(ko))*stabio(i,j,ko)))
             
             s1o(i,j,1)=s1o(i,j,1)*wtdecay*stabeps                            &
                  /(stabeps+0.5*(stabio(i,j,k)+(1.-zsurf(ko))*stabio(i,j,ko)))
   
             if (iMEAN.ne.0)then           
                if (LDIFFDIAG) then
                   wtmix(i,j,k)=s1o(i,j,1)
                endif
             endif
             

             dudo=almzer + weto(i,j,k) * di(k)**2                                  &
                  * (   ( uko(i-1,j,k) - uko(i-1,j,ko) )**2                        &
                  + ( uko(i,j,k)   - uko(i,j,ko)   )**2                        &
                  + ( vke(i,j-1,k) - vke(i,j-1,ko) )**2                        &
                  + ( vke(i,j,k)   - vke(i,j,ko)   )**2   ) *0.5
             
             hho=weto(i,j,k)*(amsuo(i,j,k)+amsuo(i-1,j,k)+amsue(i,j-1,k)       &
                  +amsue(i,j,k))*0.25
             
             rinumo=hho*MAX(gor*stabio(i,j,k)/dudo,0.)

             if (iMEAN.ne.0)then
                if (LDIFFDIAG) then
                   rinu(i,j,k)=rinumo
                endif
             endif
                          
             switc=(half-SIGN(half,stabio(i,j,k)))*weto(i,j,k)

             avo1 = avo(i,j,k)
             avo1 = (relax*MIN(avo1,av0+abackv(k))+relne*(t1o(i,j,1)     &
                  +av0/((1.+ cra*rinumo)**2)+abackv(k)+topnah))*weto(i,j,k)   
             
             avo(i,j,k)=weto(i,j,k)*MAX(cavocon*(1.e-11-stabio(i,j,k))      &
                  /(1.e-11+ABS(stabio(i,j,k))),avo1)
             !
             dvo1 = dvo(i,j,k)

             dvo1 = (relax*MIN(dvo1,dv0+s1o(i,j,1))+relne*(s1o(i,j,1)    &
                  +dv0/((1.+ crd*rinumo)**3)+dbackv(k)+topnah))*weto(i,j,k)   
             
#ifndef UMKLAP
#ifdef NURMISCH
             dvo(i,j,k) = dvo1
#else /*NURMISCH*/
             dvo(i,j,k)=weto(i,j,k)*MAX(cdvocon*(1.e-11-stabio(i,j,k))        &
                  /(1.e-11+ABS(stabio(i,j,k))),dvo1)
#endif /*NURMISCH*/
#else
             dvo(i,j,k) = dvo1
#endif /*UMKLAP*/

             stabio(i,j,k)=MAX(stabio(i,j,k),0.)

          ENDDO
       ENDDO
    ENDDO ! j-loop
!$omp end do
    

!$omp single
    CALL bounds_exch('p',s1o(:,:,1),'mo_octher 10')
    CALL bounds_exch('p',t1o(:,:,1),'mo_octher 11')
!#ifdef bounds_exch_save
    CALL bounds_exch('p',dvo,'mo_octher 12')
    CALL bounds_exch('p',avo,'mo_octher 13')
    CALL bounds_exch('p',stabio,'mo_octher 14')
!#endif    
!$omp end single


    avo(:,:,kep)=0.
    dvo(:,:,kep)=0.

!$omp end parallel

  END SUBROUTINE calc_rinum



  SUBROUTINE calc_dens

    INTEGER i,j,k
!$omp parallel private(i,j,k)
!$omp do
    DO j=1,je

       DO i=1,ie
          shelp(i,j)=0.0
          thelp(i,j)=0.0
       ENDDO

       !uwe  include downward propagation of tho and sao in land

       DO k=2,ke
          DO i=1,ie
             IF(weto(i,j,k).LT.0.5)THEN
                tho(i,j,k)=tho(i,j,k-1)
                sao(i,j,k)=sao(i,j,k-1)
             ENDIF
          ENDDO
       ENDDO

       !uwe  compute pressure and density new

       DO k=1,ke
          DO i=1,ie
             thelp(i,j)=tho(i,j,k)
             shelp(i,j)=sao(i,j,k)
          ENDDO
          CALL adisitj(thelp,shelp,preff(k),j)
          CALL rho1j(thelp,shelp,preff(k),rhelp,j)

          !uwe   po noch sauber bestimmen, aus rhoo noch machen!!!

          DO i=1,ie
             rhoo(i,j,k)=rhelp(i,j)
          ENDDO
       ENDDO
    ENDDO 
!$omp end do
!$omp end parallel

  END SUBROUTINE calc_dens

 END MODULE mo_octher
