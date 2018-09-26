MODULE mo_ncar_ocean_fluxes

  USE mo_kind
  USE mo_parallel
  USE mo_param1, ONLY: ie,je,ie_g,je_g
  USE mo_units
  USE mo_commo1, ONLY: lmonts,lmont1,monlen,ldays,lday
  USE mo_commoau1, ONLY: tfrez,tmelt,clb,con


  IMPLICIT NONE

  REAL,ALLOCATABLE,DIMENSION(:,:) :: coru10,corv10,cort10,       &
                                     corq10,corqdlw,corqdsw,     &
                                     corprec,corrunoff
CONTAINS

  SUBROUTINE alloc_mem_core

    USE mo_param1, ONLY: ie,je
 
    IMPLICIT NONE

    ALLOCATE(coru10(ie,je),corv10(ie,je),cort10(ie,je),corq10(ie,je),    &
         corqdlw(ie,je),corqdsw(ie,je),corprec(ie,je),corrunoff(ie,je))   

  END SUBROUTINE alloc_mem_core

  SUBROUTINE open_core 

    USE mo_param1, ONLY: ie,je
    IMPLICIT NONE

    INTEGER :: i,j  

    OPEN(io_in_coru10,file='CORU10',form='unformatted')
    OPEN(io_in_corv10,file='CORV10',form='unformatted')
    OPEN(io_in_cort10,file='CORT10',form='unformatted')
    OPEN(io_in_corq10,file='CORQ10',form='unformatted')
    OPEN(io_in_corlw,file='CORQDLW',form='unformatted')
    OPEN(io_in_corsw,file='CORQDSW',form='unformatted')
    OPEN(io_in_corpr,file='CORPREC',form='unformatted')
    OPEN(io_in_corro,file='CORRIV',form='unformatted')
 
    DO j=1,je
       DO i=1,ie
          coru10(i,j)=0.      ! 10m u-wind speed [m s-1]
          corv10(i,j)=0.      ! 10m v-wind speed [m s-1]
          cort10(i,j)=0.      ! 10m air temperature [C]
          corq10(i,j)=0.      ! 10m humidity [??]
          corqdlw(i,j)=0.     ! downward longwave[W m-2]
          corqdsw(i,j)=0.     ! downward shortwave[W m-2]
          corprec(i,j)=0.     ! precipitation [kg m-2]
          corrunoff(i,j)=0.   ! runoff [kg m-2]
       ENDDO
    ENDDO

  END SUBROUTINE open_core

  SUBROUTINE spool_core



    IMPLICIT NONE

    INTEGER(i8) :: ii1,ii2,ii3,ii4
    INTEGER     :: lmont  

    !spool the fields to the actual month

    IF ( lmont1 .GT. 1 ) THEN
       DO lmont=1,lmont1-1
          DO lday=1,monlen(lmont)

             WRITE(IO_STDOUT,*)'in spool'
 
             IF(p_pe==p_io) THEN
                READ(io_in_coru10)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_coru10,coru10)

             IF(p_pe==p_io) THEN
                READ(io_in_corv10)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corv10,corv10)

             IF(p_pe==p_io) THEN
                READ(io_in_cort10)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_cort10,cort10)

             IF(p_pe==p_io) THEN
                READ(io_in_corq10)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corq10,corq10)
 
             IF(p_pe==p_io) THEN
                READ(io_in_corlw)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corlw,corqdlw)
   
             IF(p_pe==p_io) THEN
                READ(io_in_corsw)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corsw,corqdsw)
   
             IF(p_pe==p_io) THEN
                READ(io_in_corpr)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corpr,corprec)

             IF(p_pe==p_io) THEN
                READ(io_in_corro)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corro,corrunoff)
   
             WRITE(io_stdout,*)'spool: ',lmont,lday,ii1
          ENDDO
       ENDDO
       WRITE(io_stdout,*) 'forcing data is spooled to start of month ',lmont1
    ENDIF

  END SUBROUTINE spool_core


  SUBROUTINE read_core 

  IMPLICIT NONE

  INTEGER(i8) :: ii1,ii2,ii3,ii4  

             WRITE(io_stdout,*)'in read'

             IF(p_pe==p_io) THEN
                READ(io_in_coru10)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_coru10,coru10)

             IF(p_pe==p_io) THEN
                READ(io_in_corv10)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corv10,corv10)

             IF(p_pe==p_io) THEN
                READ(io_in_cort10)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_cort10,cort10)

             IF(p_pe==p_io) THEN
                READ(io_in_corq10)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corq10,corq10)
 
             IF(p_pe==p_io) THEN
                READ(io_in_corlw)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corlw,corqdlw)
   
             IF(p_pe==p_io) THEN
                READ(io_in_corsw)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corsw,corqdsw)
   
             IF(p_pe==p_io) THEN
                READ(io_in_corpr)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corpr,corprec)

         
             IF(p_pe==p_io) THEN
                READ(io_in_corro)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corro,corrunoff)

             WRITE(io_stdout,*)'after read'

  END SUBROUTINE read_core

  SUBROUTINE buget_ocean_core(uo,vo,to,qla,qse,qnsw,qnlw,qpre,qnet,fw,taux,tauy)

    USE mo_param1, ONLY: ie,je
    IMPLICIT NONE
    
    REAL,INTENT(in),DIMENSION(ie,je) :: uo,vo,to

    REAL,INTENT(out),DIMENSION(ie,je) :: qpre,qla,qse,qnsw,qnlw,qnet
    REAL,INTENT(out),DIMENSION(ie,je) :: fw,taux,tauy

    REAL,DIMENSION(ie,je) :: evap,du,dv,uvdel,qo,z,ce,cd,ch,ustar,bstar
    REAL :: ra,q1,q2,f1,vlamda,cp,flamda,albw,tv,emiss,stebo
    INTEGER :: i,j,n_itts


    n_itts=2                    ! number of iterations in ncar_ocean_fluxes  
                                ! ocean: n_itts=2 ; sea ice: n_itts=5 
    ra=1.22                     ! near surface air density [kg/m3]
    q1=640380.                   ! coefficiant of q saturation function [kg/m3]
    q2=-5107.4                  ! coefficiant of q saturation function [k]
    f1=0.98                     ! "saturation effect" applied over sea water
    flamda=3.337e5              ! latent heat of fusion [J kg-1] 
    vlamda=2.5e6                ! latent heat of vaporisation [J kg-1]  
    cp=1005                     ! specific heat of air [J kg-1 K-1]
    albw=0.066                  ! albedo of sea water
    stebo=5.67e-8               ! Stefan-Boltzmann Constant [W m-2 K-4]
    emiss=1.0                   ! emissivity of sea water
    
    qpre(:,:)=0.
    qla(:,:)=0.
    qse(:,:)=0.
    qnsw(:,:)=0.
    qnlw(:,:)=0.
    qnet(:,:)=0.
    fw(:,:)=0.
    taux(:,:)=0.
    tauy(:,:)=0.
    evap(:,:)=0.
    du(:,:)=0.
    dv(:,:)=0.
    uvdel(:,:)=0.
    qo(:,:)=0.
    z(:,:)=0.
    ce(:,:)=0.
    cd(:,:)=0.
    ch(:,:)=0.
    ustar(:,:)=0.
    bstar(:,:)=0.

    DO j=2,je-1
       DO i=2,ie-1
          du(i,j)=(coru10(i,j)-uo(I,j))
          dv(i,j)=(corv10(i,j)-vo(i,j))
          uvdel(i,j)=SQRT(du(i,j)**2+dv(i,j)**2)
          qo(i,j)=(1./ra)*f1*q1*EXP(q2/to(i,j))          ! L-Y eqn. 5
          z(i,j)=10.
       ENDDO
    ENDDO

    CALL bounds_exch('u',du,'mo_ncar_ocean_fluxes 1')
    CALL bounds_exch('v',dv,'mo_ncar_ocean_fluxes 2')
    CALL bounds_exch('p',uvdel,'mo_ncar_ocean_fluxes 3')
    CALL bounds_exch('p',qo,'mo_ncar_ocean_fluxes 4')
    CALL bounds_exch('p',z,'mo_ncar_ocean_fluxes 5')


    CALL ncar_ocean_fluxes(n_itts,uvdel,cort10,to,corq10,qo,z,cd,ch,ce,ustar,bstar)

    call bounds_exch('p',ce,'mo_ncar_ocean_fluxes 6')
    call bounds_exch('p',cd,'mo_ncar_ocean_fluxes 7')
    call bounds_exch('p',ch,'mo_ncar_ocean_fluxes 8')

    DO j=2,je-1    
       DO i=2,ie-1
 
!          ce(i,j)=1.75e-3

          ! calculate evaporation         
          evap(i,j)=ra*ce(i,j)*(corq10(i,j)-qo(i,j))*uvdel(i,j)    ! L-Y eqn. 4a 

          ! calculate precipitaion heat flux  
          IF (cort10(i,j)>273.16) THEN                                
             qpre(i,j)=flamda*corprec(i,j)                        ! L-Y eqn. 14
          ENDIF

          ! calculate latent heat flux
          qla(i,j)=vlamda*evap(i,j)                           ! L-Y eqn. 4b 

          ! calculate sensible heat flux
          tv = cort10(i,j)

          qse(i,j)=ra*cp*ch(i,j)*(tv-to(i,j))*uvdel(i,j)  ! L-Y eqn. 4c 

          ! calculate net short wave heat flux

          qnsw(i,j)=corqdsw(i,j)*(1.-albw)                       ! L-Y eqn. 11
          

          ! calculate net long wave heat flux


          qnlw(i,j)=corqdlw(i,j)-emiss*stebo*to(i,j)**4             ! L-Y eqn. 12  

          ! calculate net heat flux 
          
          qnet(i,j)=-(qnsw(i,j)+qnlw(i,j) + qla(i,j)                      &
	  	               +qse(i,j)  + qpre(i,j)) /clb          ! L-Y eqn. 3a

          ! calculate net fw flux
          fw(i,j)=corprec(i,j)
!                +evap(i,j)+corrunoff(i,j)

          ! calculate wind stress
          taux(i,j)=ra*cd(i,j)*uvdel(i,j)*du(i,j)             ! L-Y eqn. 4d
          tauy(i,j)=ra*cd(i,j)*uvdel(i,j)*dv(i,j)             ! L-Y eqn. 4d
       ENDDO
    ENDDO

    call  bounds_exch('u',taux,'mo_ncar_ocean_fluxes 9')
    call  bounds_exch('v',tauy,'mo_ncar_ocean_fluxes 10')
    call  bounds_exch('p',fw,'mo_ncar_ocean_fluxes 11')
    call  bounds_exch('p',qnet,'mo_ncar_ocean_fluxes 12')
    call  bounds_exch('p',qnlw,'mo_ncar_ocean_fluxes 13')
    call  bounds_exch('p',qnsw,'mo_ncar_ocean_fluxes 14')
    call  bounds_exch('p',qse,'mo_ncar_ocean_fluxes 15')
    call  bounds_exch('p',qla,'mo_ncar_ocean_fluxes 16')
    call  bounds_exch('p',qpre,'mo_ncar_ocean_fluxes 17')
    call  bounds_exch('p',evap,'mo_ncar_ocean_fluxes 18')
    
  END SUBROUTINE buget_ocean_core

 
  SUBROUTINE buget_ice_core(ui,vi,ti,sictho,sicsno,          &
                            qla,qse,qnsw,qnlw,qres,qcon,fw,taux,tauy)

    USE mo_param1, ONLY: ie,je
    IMPLICIT NONE

    REAL,INTENT(in),DIMENSION(ie,je) :: ui,vi,ti,sictho,sicsno

    REAL,INTENT(out),DIMENSION(ie,je) :: qla,qse,qnsw,qnlw,qres,qcon
    REAL,INTENT(out),DIMENSION(ie,je) :: fw,taux,tauy

    REAL,DIMENSION(ie,je) :: evap,qnet,du,dv,uvdel,qo,z,ce,cd,ch,ustar,bstar
    REAL :: ra,q1,q2,f1,vlamda,cp,flamda,albw,tv,emiss,stebo
    REAL :: albvdr, albndr, albvdf, albndf
    INTEGER :: i,j,n_itts

    REAL :: tb

    n_itts=5                    ! number of iterations in ncar_ocean_fluxes  
    ! ocean: n_itts=2 ; sea ice: n_itts=5 see L&Y page 10


    ra=1.22                     ! near surface air density [kg/m3]
    q1=640380                   ! coefficiant of q saturation function [kg/m3]
    q2=-5107.4                  ! coefficiant of q saturation function [k]
    f1=0.98                     ! "saturation effect" applied over sea water
    flamda=3.337e5              ! latent heat of fusion [J kg-1] 
    vlamda=2.5e6                ! latent heat of vaporisation [J kg-1]  
    cp=1005                     ! specific heat of air [J kg-1 K-1]

    albw=0.066                  ! albedo of sea water
    albvdr=0.95                 ! Fraction of dsw 0.29 ; visible direct albedo of snow/sea ice
    albvdf=0.85                 ! Fraction of dsw 0.31 ; visible difuse albedo of snow/sea ice
    albndr=0.5                  ! Fraction of dsw 0.24 ; near infrared direct albedo of snow/sea ice   
    albndf=0.4                  ! Fraction of dsw 0.16 ; near infrared difuse albedo of snow/sea ice

    stebo=5.67e-8               ! Stefan-Boltzmann Constant [W m-2 K-4]
    emiss=0.95                  ! emissivity of sea ice 

    tb=tfrez+tmelt

!    call flx_blk_albedo(ti,sictho,sicsno,palb,palcn,palbp,palcnp)
!           palb         , & !  albedo of ice under overcast sky
!           palcn        , & !  albedo of ocean under overcast sky
!           palbp        , & !  albedo of ice under clear sky
!           palcnp           !  albedo of ocean under clear sky
    

    DO j=2,je-1
       DO i=2,ie-1
          du(i,j)=(coru10(i,j)-ui(I,j))
          dv(i,j)=(corv10(i,j)-vi(i,j))
          uvdel(i,j)=SQRT(du(i,j)**2+dv(i,j)**2)
          qo(i,j)=(1./ra)*f1*q1*EXP(q2/ti(i,j))          ! L-Y eqn. 5
          z(i,j)=10.
       ENDDO
    ENDDO

!    CALL ncar_ocean_fluxes(n_itts,uvdel,cort10,ti,corq10,qo,z,cd,ch,ce,ustar,bstar)

    DO j=2,je-1
       DO i=2,ie-1

          ch(i,j) = 1.63e-3                                   ! L-Y eqn. 21 
          cd(i,j) = 1.63e-3                                   ! L-Y eqn. 21 
          ce(i,j) = 1.63e-3                                   ! L-Y eqn. 21 

          ! calculate evaporation         
          evap(i,j)=ra*ce(i,j)*(corq10(i,j)-qo(i,j))*uvdel(i,j)   ! L-Y eqn. 20a 

          ! calculate latent heat flux
          qla(i,j)=vlamda*evap(i,j)                           ! L-Y eqn. 20b 

          ! calculate sensible heat flux
          tv = cort10(i,j)

          qse(i,j)=ra*cp*ch(i,j)*(tv-ti(i,j))*uvdel(i,j) ! L-Y eqn. 20c 

          ! calculate net short wave heat flux                ! L-Y eqn. 22
          qnsw(i,j)=corqdsw(i,j)*(0.29*(1.-albvdr)             &  
                              +0.31*(1.-albvdf)             &  
                              +0.24*(1.-albndr)             &  
                              +0.16*(1.-albndf))        


          ! calculate net long wave heat flux

          qnlw(i,j)=emiss*corqdlw(i,j)-emiss*stebo*ti(i,j)**4  ! L-Y eqn. 23  

          ! calculate net heat flux 
          qnet(i,j)=qnsw(i,j)+qnlw(i,j)+qla(i,j)       &
               +qse(i,j)                                          ! L-Y eqn. 19a

          ! calculate net fw flux
          fw(i,j)=corprec(i,j)
!                             +evap(i,j)                         ! L-Y eqn. 19b

          ! calculate wind stress
          taux(i,j)=ra*cd(i,j)*uvdel(i,j)*du(i,j)                ! L-Y eqn. 20d
          tauy(i,j)=ra*cd(i,j)*uvdel(i,j)*dv(i,j)                ! L-Y eqn. 20d

          qcon(i,j)=((tb-ti(i,j))/sictho(i,j)*con)/clb           ! conductive hflx

          qres(i,j)=(-qnet(i,j)/clb)-qcon(i,j)                   ! residual hflx  
                     

       ENDDO
    ENDDO

    call  bounds_exch('u',taux,'mo_ncar_ocean_fluxes 17')
    call  bounds_exch('v',tauy,'mo_ncar_ocean_fluxes 18')
    call  bounds_exch('p',fw,'mo_ncar_ocean_fluxes 19')
    call  bounds_exch('p',qnet,'mo_ncar_ocean_fluxes 20')
    call  bounds_exch('p',qnlw,'mo_ncar_ocean_fluxes 21')
    call  bounds_exch('p',qnsw,'mo_ncar_ocean_fluxes 22')
    call  bounds_exch('p',qse,'mo_ncar_ocean_fluxes 23')
    call  bounds_exch('p',qla,'mo_ncar_ocean_fluxes 24')
    call  bounds_exch('p',evap,'mo_ncar_ocean_fluxes 25')
    call  bounds_exch('p',qres,'mo_ncar_ocean_fluxes 26')
    call  bounds_exch('p',qcon,'mo_ncar_ocean_fluxes 27')

  END SUBROUTINE buget_ice_core


  SUBROUTINE ncar_ocean_fluxes(n_itts,udel,t,ts,q,qs,z,cd,ch,ce,ustar,bstar)

    USE mo_param1, ONLY: ie,je
    IMPLICIT NONE

    REAL, PARAMETER :: GRAV   = 9.80
    REAL, PARAMETER :: VONKARM = 0.40

    REAL,INTENT(in) :: udel(ie,je), t(ie,je), ts(ie,je),        &
         q(ie,je), qs(ie,je), z(ie,je)
    REAL,INTENT(inout) :: cd(ie,je), ch(ie,je), ce(ie,je),       &
         ustar(ie,je), bstar(ie,je)

    REAL :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
    REAL :: cd_rt                                ! full drag coefficients @ z
    REAL :: zeta, x2, x, psi_m, psi_h            ! stability parameters
    REAL :: u, u10, tv, tstar, qstar, z0, xx, stab

!    INTEGER, PARAMETER :: n_itts = 2
    INTEGER               i, j, jj,n_itts

    DO i=2,ie-1
       DO j=2,je-1
          !          if (avail(i,j)) then
          tv = t(i,j)*(1.+0.608*q(i,j));
          u = MAX(udel(i,j), 0.5);                                 ! 0.5 m/s floor on wind (undocumented NCAR)
          u10 = u;                                                ! first guess 10m wind

          cd_n10 = (2.7/u10+0.142+0.0764*u10)/1.e3;                ! L-Y eqn. 6a
          cd_n10_rt = SQRT(cd_n10);
          ce_n10 =                     34.6 *cd_n10_rt/1.e3;       ! L-Y eqn. 6b
          stab = 0.5 + SIGN(0.5,t(i,j)-ts(i,j))
          ch_n10 = (18.0*stab+32.7*(1.-stab))*cd_n10_rt/1.e3;       ! L-Y eqn. 6c

          cd(i,j) = cd_n10;                                         ! first guess for exchange coeff's at z
          ch(i,j) = ch_n10;
          ce(i,j) = ce_n10;
          DO jj=1,n_itts                                           ! Monin-Obukhov iteration
             cd_rt = SQRT(cd(i,j));
             ustar(i,j) = cd_rt*u;                                   ! L-Y eqn. 7a
             tstar    = (ch(i,j)/cd_rt)*(t(i,j)-ts(i,j));                ! L-Y eqn. 7b
             qstar    = (ce(i,j)/cd_rt)*(q(i,j)-qs(i,j));                ! L-Y eqn. 7c
             bstar(i,j) = grav*(tstar/tv+qstar/(q(i,j)+1/0.608));
             zeta     = vonkarm*bstar(i,j)*z(i,j)/(ustar(i,j)*ustar(i,j)); ! L-Y eqn. 8a
             zeta     = SIGN( MIN(ABS(zeta),10.0), zeta );         ! undocumented NCAR
             x2 = SQRT(ABS(1.-16.*zeta));                            ! L-Y eqn. 8b
             x2 = MAX(x2, 1.0);                                    ! undocumented NCAR
             x = SQRT(x2);

             IF (zeta > 0) THEN
                psi_m = -5.*zeta;                                    ! L-Y eqn. 8c
                psi_h = -5.*zeta;                                    ! L-Y eqn. 8c
             ELSE
                psi_m = LOG((1.+2.*x+x2)*(1.+x2)/8.)-2.*(ATAN(x)-ATAN(1.0)); ! L-Y eqn. 8d
                psi_h = 2*LOG((1.+x2)/2.);                                ! L-Y eqn. 8e
             END IF

             u10 = u/(1.+cd_n10_rt*(LOG(z(i,j)/10.)-psi_m)/vonkarm);       ! L-Y eqn. 9
             cd_n10 = (2.7/u10+0.142+0.0764*u10)/1.e3;                  ! L-Y eqn. 6a again
             cd_n10_rt = SQRT(cd_n10);
             ce_n10 = 34.6*cd_n10_rt/1.e3;                              ! L-Y eqn. 6b again
             stab = 0.5 + SIGN(0.5,zeta)
             ch_n10 = (18.0*stab+32.7*(1.-stab))*cd_n10_rt/1.e3;         ! L-Y eqn. 6c again
             z0 = 10.*EXP(-vonkarm/cd_n10_rt);                          ! diagnostic

             xx = (LOG(z(i,j)/10.)-psi_m)/vonkarm;
             cd(i,j) = cd_n10/(1.+cd_n10_rt*xx)**2;                       ! L-Y 10a
             xx = (LOG(z(i,j)/10.)-psi_h)/vonkarm;
             ch(i,j) = ch_n10/(1.+ch_n10*xx/cd_n10_rt)**2;                !       b
             ce(i,j) = ce_n10/(1.+ce_n10*xx/cd_n10_rt)**2;                !       c
          END DO
          !          end if
       END DO
    END DO
 
  


  END SUBROUTINE ncar_ocean_fluxes

  SUBROUTINE flx_blk_albedo( tice, sictho, sicsno, palb , palcn , palbp , palcnp )
  
     USE MO_PARAM1

    !!----------------------------------------------------------------------
    !!               ***  ROUTINE flx_blk_albedo  ***
    !!
    !! ** Purpose :   Computation of the albedo of the snow/ice system
    !!      as well as the ocean one
    !!
    !! ** Method  : - Computation of the albedo of snow or ice (choose the
    !!      right one by a large number of tests
    !!              - Computation of the albedo of the ocean
    !!
    !! References :
    !!      Shine and Hendersson-Sellers 1985, JGR, 90(D1), 2243-2250.
    !!
    !! History :
    !!  8.0   !  01-04  (LIM 1.0)
    !!  8.5   !  03-07  (C. Ethe, G. Madec)  Optimization (old name:shine)
    !!----------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL, DIMENSION(ie,je), INTENT(out) ::  &
           palb         , & !  albedo of ice under overcast sky
           palcn        , & !  albedo of ocean under overcast sky
           palbp        , & !  albedo of ice under clear sky
           palcnp                  !  albedo of ocean under clear sky
    
    !! * Local variables
    INTEGER ::   i,j        ! dummy loop indices
    REAL ::   &
         c1     = 0.05  , & ! constants values
         c2     = 0.1   , &
         albice = 0.50  , & !  albedo of melting ice in the arctic and antarctic (Shine & Hendersson-Sellers)
         cgren  = 0.06  , & !  correction of the snow or ice albedo to take into account
         !  effects of cloudiness (Grenfell & Perovich, 1984)
         alphd  = 0.80  , & !  coefficients for linear interpolation used to compute
         alphdi = 0.72  , & !  albedo between two extremes values (Pyane, 1972)
         alphc  = 0.65  , &
         zmue   = 0.4   , &     !  cosine of local solar altitude
         zzero   = 0.0  ,  &
         zone    = 1.0  ,  &
         rt0_snow=0.0   ,  &
         rt0_ice=0.0


    REAL ::   &
         zmue14         , & !  zmue**1.4
         zalbpsnm       , & !  albedo of ice under clear sky when snow is melting
         zalbpsnf       , & !  albedo of ice under clear sky when snow is freezing
         zalbpsn        , & !  albedo of snow/ice system when ice is coverd by snow
         zalbpic        , & !  albedo of snow/ice system when ice is free of snow
         zithsn         , & !  = 1 for hsn >= 0 ( ice is cov. by snow ) ; = 0 otherwise (ice is free of snow)
         zitmlsn        , & !  = 1 freezinz snow (sist >=rt0_snow) ; = 0 melting snow (sist<rt0_snow)
         zihsc1         , & !  = 1 hsn <= c1 ; = 0 hsn > c1
         zihsc2                   !  = 1 hsn >= c2 ; = 0 hsn < c2

    REAL, DIMENSION(ie,je) ::  &
         zalbfz         , & !  ( = alphdi for freezing ice ; = albice for melting ice )
         zficeth                  !  function of ice thickness
    LOGICAL , DIMENSION(ie,je) ::  &
         llmask
  
    REAL, DIMENSION(ie,je) ::      &  
          tice   ,                 &   !: Sea-Ice Surface Temperature (Kelvin )
          sicsno  ,                &   !: Snow thickness
          sictho                       !: Ice thickness
    
    
    !-------------------------
    !  Computation of  zficeth
    !--------------------------
    
    llmask = (sicsno == 0.0) .AND. ( tice >= rt0_ice )
    WHERE ( llmask )   !  ice free of snow and melts
       zalbfz = albice
    ELSEWHERE
       zalbfz = alphdi
    END WHERE
    
    DO i = 1, ie
       DO j = 1, je
          IF( sictho(i,j) > 1.5 ) THEN
             zficeth(i,j) = zalbfz(i,j)
          ELSEIF( sictho(i,j) > 1.0  .AND. sictho(i,j) <= 1.5 ) THEN
             zficeth(i,j) = 0.472 + 2.0 * ( zalbfz(i,j) - 0.472 ) * ( sictho(i,j) - 1.0 )
          ELSEIF( sictho(i,j) > 0.05 .AND. sictho(i,j) <= 1.0 ) THEN
             zficeth(i,j) = 0.2467 + 0.7049 * sictho(i,j)                                  &
                          - 0.8608 * sictho(i,j) * sictho(i,j)                 &
                          + 0.3812 * sictho(i,j) * sictho(i,j) * sictho (i,j)
          ELSE
             zficeth(i,j) = 0.1 + 3.6 * sictho(i,j)
          ENDIF
       END DO
    END DO
    
    !-----------------------------------------------
    !    Computation of the snow/ice albedo system
    !-------------------------- ---------------------
    
    !    Albedo of snow-ice for clear sky.
    !-----------------------------------------------
    DO i = 1, ie
       DO j = 1, je
          !  Case of ice covered by snow.
          !  melting snow
          zihsc1       = 1.0 - MAX ( zzero , SIGN ( zone , - ( sicsno(i,j) - c1 ) ) )
          zalbpsnm     = ( 1.0 - zihsc1 ) & 
               * ( zficeth(i,j) + sicsno(i,j) * ( alphd - zficeth(i,j) ) / c1 ) &
               &                 + zihsc1   * alphd
          !  freezing snow
          zihsc2       = MAX ( zzero , SIGN ( zone , sicsno(i,j) - c2 ) )
          zalbpsnf     = ( 1.0 - zihsc2 ) * & 
               ( albice + sicsno(i,j) * ( alphc - albice ) / c2 ) &
               &                 + zihsc2   * alphc
          
          zitmlsn      =  MAX ( zzero , SIGN ( zone , tice(i,j) - rt0_snow ) )
          zalbpsn      =  zitmlsn * zalbpsnf + ( 1.0 - zitmlsn ) * zalbpsnm
          
          !  Case of ice free of snow.
          zalbpic      = zficeth(i,j)
          
          ! albedo of the system
          zithsn       = 1.0 - MAX ( zzero , SIGN ( zone , - sicsno(i,j) ) )
          palbp(i,j) =  zithsn * zalbpsn + ( 1.0 - zithsn ) *  zalbpic
       END DO
    END DO
    
    !    Albedo of snow-ice for overcast sky.
    !----------------------------------------------
    palb(:,:)   = palbp(:,:) + cgren
    
    !--------------------------------------------
    !    Computation of the albedo of the ocean
    !-------------------------- -----------------
    
    
    !  Parameterization of Briegled and Ramanathan, 1982
    zmue14      = zmue**1.4
    palcnp(:,:) = 0.05 / ( 1.1 * zmue14 + 0.15 )
    
    !  Parameterization of Kondratyev, 1969 and Payne, 1972
    palcn(:,:)  = 0.06
    
  END SUBROUTINE flx_blk_albedo

  SUBROUTINE normpem

! close fw budget for each timestep as suggested in the core release notes
  
  USE mo_kind
  USE mo_parallel
  USE mo_param1, ONLY: ie,je,ie_g,je_g
  USE mo_units
  USE mo_commo1, ONLY: zo,sao,eminpo,dlxp,dlyp,sictho &
       ,sicsno,ddpo,weto,dt
  USE mo_commoau2, ONLY: prech
  USE mo_commoau1, ONLY: rhoicwa,rhosnwa

  REAL :: rimbl,rarea, re, zold, znew, delz, draft
  INTEGER :: i,j

  rimbl=0.
  rarea=0.
  DO j=2,je-1
     DO i=2,ie-1        
        rimbl=rimbl+(prech(i,j)+eminpo(i,j))          &
             *dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
        rarea=rarea+dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
     ENDDO
  ENDDO

  CALL global_sum(rimbl)
  CALL global_sum(rarea)

  re=rimbl/rarea

  WRITE(IO_STDOUT,*)'core:  before normalisation P-E+R+r [Sv]: ',re*rarea*1.e-6
 
    DO j=2,je-1
     DO i=2,ie-1       
	if (weto(i,j,1).ge.0.5) then	
    	    draft=ddpo(i,j,1)+zo(i,j)-((sictho(i,j)*rhoicwa) &
        	-(sicsno(i,j)*rhosnwa))
    	    zold=zo(i,j)
    	    znew=zold-(re*dt)
    	    zo(i,j)=znew
    	    eminpo(i,j)=eminpo(i,j)-re
    	    delz=znew-zold       
    	    sao(i,j,1)=sao(i,j,1)*(draft/(delz+draft))
	endif     
     ENDDO
    ENDDO     

  rimbl=0.
  rarea=0.
  DO j=2,je-1
     DO i=2,ie-1        
        rimbl=rimbl+(prech(i,j)+eminpo(i,j))          &
             *dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
        rarea=rarea+dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
     ENDDO
  ENDDO

  CALL global_sum(rimbl)
  CALL global_sum(rarea)

  re=rimbl/rarea

  WRITE(IO_STDOUT,*)'core: after normalisation P-E+R+r [Sv]: ',re*rarea*1.e-6

  CALL bounds_exch('p',zo,'mo_ncar_ocean_fluxes 28')
  CALL bounds_exch('p',sao,'mo_ncar_ocean_fluxes 29')
  CALL bounds_exch('p',eminpo,'mo_ncar_ocean_fluxes 39')

  END SUBROUTINE normpem

END MODULE mo_ncar_ocean_fluxes















