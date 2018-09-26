MODULE mo_omip

  USE mo_kind
  USE mo_parallel
  USE mo_param1, ONLY: ie,je,ie_g,je_g 
  USE mo_units, ONLY: io_in_gwix,io_in_gwiy,io_in_gtem,     &
                      io_in_gpre,io_in_gswr,io_in_gtde,     &
                      io_in_gu10,io_in_rval,io_in_gslp,     &
                      io_in_glwr,io_in_gclo,io_stdout

  USE mo_commo1, ONLY: lmonts,lmont1,monlen,ldays,lday,     &
                       tafo,txo,tye,fprec,fswr,lyear,       &
                       ftdew,fu10,fclou,giriv,fslp

  USE mo_commoau1

  IMPLICIT NONE

  
CONTAINS

  SUBROUTINE open_omip

    IMPLICIT NONE

    OPEN(io_in_gwix,file='GIWIX',status='unknown',                  &
         access='sequential',form='unformatted')

    OPEN(io_in_gwiy,file='GIWIY',status='unknown',                  &
         access='sequential',form='unformatted')

    OPEN(io_in_gtem,file='GITEM',status='unknown',                  &
         access='sequential',form='unformatted')

    OPEN(io_in_gpre,file='GIPREC',status='unknown',                 &
         access='sequential',form='unformatted')

    OPEN(io_in_gswr,file='GISWRAD',status='unknown',                &
         access='sequential',form='unformatted')

    OPEN(io_in_gtde,file='GITDEW',status='unknown',                 &
         access='sequential',form='unformatted')

    OPEN(io_in_gu10,file='GIU10',status='unknown',                  &
         access='sequential',form='unformatted')

#ifdef RIVER_GIRIV
    OPEN(io_in_rval,file='GIRIV',status='unknown',                  &
         access='sequential',form='unformatted')
#endif /*RIVER_GIRIV*/

#ifdef FORCE_SLP
    OPEN(io_in_gslp,file='GIPRESS',status='unknown',                &
         access='sequential',form='unformatted')
#endif

#ifdef FORCE_DWLW
    OPEN(io_in_glwr,file='GIDWLW',status='unknown',                 &
         access='sequential',form='unformatted')
#else
    OPEN(io_in_gclo,file='GICLOUD',status='unknown',                &
         access='sequential',form='unformatted')
#endif /*FORCE_DWLW*/


  END SUBROUTINE open_omip



  SUBROUTINE spool_omip

    IMPLICIT NONE

    INTEGER(i8) :: ii1,ii2,ii3,ii4
    INTEGER     :: lmont  

    !spool the fields to the actual month

    IF ( lmont1 >= 1 ) THEN
       DO lmont=1,lmont1-1
          DO lday=1,monlen(lmont)

             WRITE(IO_STDOUT,*)'in spool'
             IF(p_pe==p_io) THEN
                READ(io_in_gtem)ii1,ii2,ii3,ii4  
             ENDIF
             CALL read_slice(io_in_gtem,tafo)

             IF(p_pe==p_io) THEN
                READ(io_in_gwix)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_gwix,txo)

             IF(p_pe==p_io) THEN
                READ(io_in_gwiy)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_gwiy,tye)

             IF(p_pe==p_io) THEN
                READ(io_in_gpre)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_gpre,fprec)

             IF(p_pe==p_io) THEN
                READ(io_in_gswr)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_gswr,fswr)

             IF(p_pe==p_io) THEN
                READ(io_in_gtde)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_gtde,ftdew)

             IF(p_pe==p_io) THEN
                READ(io_in_gu10)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_gu10,fu10)

#ifdef RIVER_GIRIV
             IF(p_pe==p_io) THEN
                READ(io_in_rval)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_rval,giriv)
#endif /*RIVER_GIRIV*/

#ifdef FORCE_SLP
             IF(p_pe==p_io) THEN
                READ(io_in_gslp)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_gslp,fslp)
#endif /*FORCE_SLP*/

#ifdef FORCE_DWLW
             IF(p_pe==p_io) THEN
                READ(io_in_glwr)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_glwr,fclou)
#else
             IF(p_pe==p_io) THEN
                READ(io_in_gclo)ii1,ii2,ii3,ii4  
             ENDIF
             CALL read_slice(io_in_gclo,fclou)
#endif /*FORCE_DWLW*/

             WRITE(io_stdout,*)'spool: ',lmont,lday,ii1
          ENDDO
       ENDDO
       WRITE(io_stdout,*) 'forcing data is spooled to start of month ',lmont1
    ENDIF

  END SUBROUTINE spool_omip



  SUBROUTINE read_omip

    IMPLICIT NONE

    INTEGER(i8) :: ii1,ii2,ii3,ii4  

       IF(p_pe==p_io) THEN
          READ(io_in_gtem)ii1,ii2,ii3,ii4  
       ENDIF
       CALL read_slice(io_in_gtem,tafo)
       IF (p_pe==p_io) THEN
          WRITE(0,'(a,i4.4,''-'',i2.2,''-'',i2.2,a,i8)') &
               'OMIP: read forcing: ', lyear,lmonts,ldays, &
               ' from file:', ii1
       ENDIF

       IF(p_pe==p_io) THEN
          READ(io_in_gwix)ii1,ii2,ii3,ii4  
       ENDIF
       CALL read_slice(io_in_gwix,txo)

       IF(p_pe==p_io) THEN
          READ(io_in_gwiy)ii1,ii2,ii3,ii4  
       ENDIF
       CALL read_slice(io_in_gwiy,tye)

       IF(p_pe==p_io) THEN 
          READ(io_in_gpre)ii1,ii2,ii3,ii4  
       ENDIF
       CALL read_slice(io_in_gpre,fprec)

       IF(p_pe==p_io) THEN
          READ(io_in_gswr)ii1,ii2,ii3,ii4  
       ENDIF
       CALL read_slice(io_in_gswr,fswr)

       IF(p_pe==p_io) THEN 
          READ(io_in_gtde)ii1,ii2,ii3,ii4
       ENDIF
       CALL read_slice(io_in_gtde,ftdew)

       IF(p_pe==p_io) THEN
          READ(io_in_gu10)ii1,ii2,ii3,ii4
       ENDIF
       CALL read_slice(io_in_gu10,fu10)

#ifdef RIVER_GIRIV
       IF(p_pe==p_io)THEN
          READ(io_in_rval)ii1,ii2,ii3,ii4
       ENDIF
       CALL read_slice(io_in_rval,giriv)
#endif /*RIVER_GIRIV*/

#ifdef FORCE_SLP
       IF(p_pe==p_io)THEN
          READ(io_in_gslp)ii1,ii2,ii3,ii4
       ENDIF
       CALL read_slice(io_in_gslp,fslp)
#endif /*FORCE_SLP*/

#ifdef FORCE_DWLW
       IF(p_pe==p_io) THEN
          READ(io_in_glwr)ii1,ii2,ii3,ii4
       ENDIF
       CALL read_slice(io_in_glwr,fclou)
       WRITE(io_stdout,*) 'attn: option force_dwlw: fclou cont. dwrd lw'
#else
       IF(p_pe==p_io) THEN 
          READ(io_in_gclo)ii1,ii2,ii3,ii4
       ENDIF
       CALL read_slice(io_in_gclo,fclou)
#endif /*FORCE_DWLW*/

       !      For periodic boundaries

       CALL bounds_exch('u',TXO,'mo_omip 1')
       CALL bounds_exch('v',TYE,'mo_omip 2') ! option "vf" syscronises also line 2
       CALL bounds_exch('p',TAFO,'mo_omip 3')
       CALL bounds_exch('p',FCLOU,'mo_omip 4')
       CALL bounds_exch('p',FSWR,'mo_omip 5')
       CALL bounds_exch('p',FU10,'mo_omip 6')
       CALL bounds_exch('p',FPREC,'mo_omip 7')
       CALL bounds_exch('p',FTDEW,'mo_omip 8')


!call p_barrier
!write(0,*)'test 0',p_pe,maxloc(tye(:,:)),maxval(tye(:,:))

  END SUBROUTINE read_omip

!!$     DO j=1,je
!!$        DO i=1,ie
!!$           aclo(i,j)=fclou(i,j)
!!$           rpreco(i,j)=fprec(i,j)
!!$        ENDDO
!!$     ENDDO

!!$  !
!!$  !  update forcing once per day
!!$  !
!!$  IF(ldays.LE.(monlen(lmonts)+1)/2)THEN
!!$     monmon=lmonts-1
!!$     IF(monmon.LT.1)monmon=12
!!$     awert=0.5+(float(ldays)-0.5)/float(monlen(lmonts))
!!$  ELSE
!!$     monmon=lmonts+1
!!$     IF(monmon>=1)monmon=1
!!$     awert=0.5+(float(monlen(lmonts)-ldays)+0.5)                     &
!!$          &      /float(monlen(lmonts))
!!$  ENDIF
!!$  ewert=1.-awert
!!$  !
!!$#ifndef RIVER_GIRIV
!!$  DO n=1,numriv
!!$     friv(n)=(awert*rival(n,lmonts)+ewert*rival(n,monmon))
!!$  ENDDO
!!$#endif /* RIVER_GIRIV*/
!!$  !
!!$  !as**forcing for growth:
!!$  DO j=1,je
!!$     DO i=1,ie
!!$        tairo(i,j)=tafo(i,j)
!!$        tdo(i,j)=ftdew(i,j)
!!$     ENDDO
!!$  ENDDO
!!$  !


  SUBROUTINE vardrag(dragl,drags,k,sst,airt,humida,humido,uvabsol)

!===============================================================================
!
!     purpose:  evaluation of drag coefficient for momentum,sensible and
!               latent heat
!
!     method:   large and pond  (1981,1982)
!
!     author:   josef maximilian oberhuber
!
!     date:     20. july 1991
!
!     modified for da-silva-intercomparison: 19. june 1996 by jmo
!
!     modified for including in hope:       02. july 1996 by fr
!                                               (frank roeske)
!
!
!===============================================================================

    USE mo_param1
    IMPLICIT NONE 

    REAL, INTENT(IN), DIMENSION(ie*je) :: airt,uvabsol,sst,humido,humida
    INTEGER, INTENT(IN) :: k
                               
    REAL, INTENT(INOUT), DIMENSION(ie*je) :: dragl,drags

    !----> define local variables

    REAL,  DIMENSION(ie*je) :: help,psim,psils,uvarian,drag
    REAL :: deltav,deltat,t0,flag,zoverl,x,        &
            xh,fvonx,dfdx,cm,csn,cmn,cln,deltaq !,sum1,sum2,sum3

    INTEGER :: iter,i


    !----> initialize parameters

    REAL, PARAMETER ::     &
        consts  = 0.0327 , &
        constl  = 0.0346 , &
        refhgt  = 10.     ,&
        charnck = 0.0187  ,&
        rkarman = 0.4     ,&
        grav    = 9.81    ,&
        rhalf   = 0.5

!   charnck = 0.0064


    !----> determine stability parameter

    DO i=1,k
       uvarian(i)=0.2*uvabsol(i)
       deltav=uvabsol(i)*(uvabsol(i)+2.*uvarian(i))
       deltat=sst(i)-airt(i)
       deltaq=humido(i)-humida(i)
       help(i)=grav*refhgt/(charnck*deltav)
       t0=airt(i)*(1.+1.7e-6*airt(i)*humida(i))
       flag=(rhalf+SIGN(rhalf,deltat))
       zoverl=-70.*refhgt*(deltat+2.5e-6*t0**2*deltaq)*(1.-flag)       &
              -100.*refhgt*(deltat+1.7e-6*t0**2*deltaq)*flag
       zoverl=zoverl/(t0*deltav)
       x=SQRT(SQRT(1.+16.*ABS(zoverl)))
       psim(i)=-(1.-flag)*7.*zoverl                                    &
              +flag*(2.*LOG((1.+x)*.5)+LOG((1.+x**2)*.5)               &
              +2.*(ATAN(1.)-ATAN(x)))
       psils(i)=(1.-flag)*psim(i)+flag*2.*LOG((1.+x**2)*.5)
       drag(i)=1.4e-3
    ENDDO

    !---->solve for drag coefficient

    DO iter=1,5
       !sum3=0.
       DO i=1,k
          xh = LOG(help(i)/drag(i))-psim(i)
          fvonx=SQRT(drag(i))*xh-rkarman
          dfdx=(2.*SQRT(drag(i)))/(xh-2./help(i))
          cm=drag(i)-fvonx*dfdx
          cm=MIN(MAX(cm,5.e-4),3.e-3)
        !  sum3=sum3+ABS(cm-drag(i))
          drag(i)=cm
       ENDDO
!       sum3=sum3/float(k)
    ENDDO


    !---->determine transfer coefficients for latent heat
    !                                 and for sensible heat

    DO i=1,k
       cmn=(rkarman/LOG(help(i)/drag(i)))**2
       cln=constl*rkarman/LOG(help(i)/drag(i))
       dragl(i)=cln*SQRT(drag(i)/cmn)                                 &
            /(1.-cln*psils(i)/(SQRT(cmn)*rkarman))
       csn=consts*rkarman/LOG(help(i)/drag(i))
       drags(i)=csn*SQRT(drag(i)/cmn)                                 &
            /(1.-csn*psils(i)/(SQRT(cmn)*rkarman))
    ENDDO
       
    dragl(1:k) = MIN( MAX(dragl(1:k),5.e-4),3.e-3 )
    drags(1:k) = MIN( MAX(drags(1:k),5.e-4),3.e-3 )
!    sum1       = SUM( drags(1:k) )
!    sum2       = SUM( dragl(1:k) )

  END SUBROUTINE vardrag




 SUBROUTINE budget(alat, fhs, fhb, t, lrhs, a, tair, td, acl, pa, ug, &
                   sln, om, hice, alb, a2, qsw, qlw, qse, qla)
!
!=======================================================================
!
!  !DESCRIPTION:
!   Calculation of growth rates for the ice-covered part of a grid
!   cell with standard bulk formulas.
!
!  !REVISION HISTORY
!   Developed by A. Stoessel, MPI, Hamburg, 1989      
!
!   Modified by S.Legutke, DKRZ, 16.04.95:
!     - separate calculation of surface/bottom melt/growth.
!
!   Modified by Uwe, MPI, 15.1.2000
!     - remove error in sw calculation (diagnostic)
!     - include update of temperature in sw-fluxes in iteration
!     - make albedo continuously dependent of temeprature
!
!   Modified by Helmuth Haak, MPI, 19.04.2000
!     - longwave-flux from berliand(1952)
!
!   Modified by Dirk Notz, MPI, 27.04.2006
!     - changed to style-guide style (incl. transition to f90)
!
!  !METHOD:
!   Ice or snow surface temperatures, respectively, are calculated by
!   iteration (regula falsi)
!
!  ! RELEVANT COMPILE-TIME SWITCHES
!    BULK_KARA: use bulk equations according to Kara, B. A., P. A. Rochford, 
!      and H. E. Hurlburt, Efficient and accurate bulk parameterizations
!      of air-sea fluxes for use in general circulation models, J. Atmos.
!      Ocean. Tech., 17, 14211438, 2000.
!=======================================================================
!
  USE mo_param1, ONLY: ie,je
  USE mo_commoau1
  USE mo_kind
  IMPLICIT NONE
  REAL, DIMENSION(ie,je), INTENT(in)   ::  &
        sln,     & ! incoming surface solar rad.     [W/m^2] 
        tair,    & ! air temperatures                [deg C]
        td,      & ! dew point temperatures          [K]
        acl,     & ! fractional cloud cover          [frac.]
        pa,      & ! atmospheric surface presure     [Pa]
        ug,      & ! 2 m wind speed                  [m/s]
        om,      & ! land/sea mask                   [0/1]
        hice,    & ! ice thickness                   [m]
        alb,     & ! surface albedo
        a2,      & ! snow flag
        alat       ! Latitude
  
  REAL, DIMENSION(ie, je), INTENT(inout)   :: &
        t          ! surface snow/ice temperature      [C]
  
  REAL, DIMENSION(ie, je), INTENT(out)   ::   &
        fhs,     & ! surface melt thick ice            [m/s]
        fhb,     & ! bottom melt/growth of thick ice   [m/s]
        qsw,     & ! absorbed solar radiation          [W/m^2]
        qlw,     & ! outgoing longwave heat flux       [W/m^2]
        qse,     & ! sensible heat flux                [W/m^2]
        qla        ! latent heat flux                  [W/m^2]
  
  REAL,  DIMENSION(ie * je) ::  &
        sln1, tair1, td1, acl1, pa1, ug1,  hice1, alb1, t1, fhs1, a21,&
        fhb1, qse1, qla1,  &   ! One dimensional arrays of the above
        esta,    & ! water vapor pressure at 2 m height  [Pa]
        esti,    & ! water vapor pressure at ice surface [Pa]
        fakts,   & !
        flai,    & !
        qlw_out1,& ! outgoing longwave radiation [W/m^2]
        qsw1,    & ! net absorbed shortwave radiation [W/m^2]
        qlw_in1, & ! incoming longwave radiation [W/m^2]
        xlat,    & !
        stp,     & ! surface temperature previous iteration
        stpp,    & ! surface temperature pre-previous iteration
        qnet,    & ! net heat flux at ice surface [W/m^2]
        qnetp,   & ! net heat flux at ice surface [W/m^2] at prev. time step
        dragl,   & !
        drags,   & !
        ugg,     & ! 
        feu,     & ! 
        fdiff,   & ! 
        flag,    & !
        sphumida,& ! saturation specific humidity
        rhoair_k,& !
        qlwin1  ,& ! Incoming longwave radiation dummy field
        sphumido   !
        
  REAL, DIMENSION(ie, je, 2), INTENT(in) :: &
        a          ! ice compactness [frac.]
  
  
  REAL :: &
        tb                        ! ice temperature at ice-ocean interface [K]
        
  REAL, PARAMETER ::           &
        rgas     = 287.1_wp,   & ! specific gas constant of dry air  [J/kg/K]
        cpa      = 1004.67_wp, & ! specific heat of dry air          [J/kg/K] 
        fr_fac   = 1.1925_wp,  & ! Frank Roeske's budget closing factor
        albtrans = 0.5_wp        ! Width of albedo transition melt to freez:
                                 ! large vale: sharp transition
                                 ! low value:  smoothe transition
  
  INTEGER, PARAMETER ::        &
        imax   = 10              ! maximum number of iteration steps
  
  INTEGER ::               &
        k,                     & ! number of grid cells with sea ice
        i,                     & ! counter in DO-loops (x index cells)
        j,                     & ! counter in DO-loops (y index cells)
        iter                     ! counter in iteration DO-loop

  INTEGER, INTENT(IN) ::   &
        lrhs                     !
        
#ifdef FORCE_DWLW
  d3=0.97*5.67e-8
#endif
  
  tb = tfrez+tmelt               ! Set ice-ocean interf. temp. (in K)
  
  
!-----------------------------------------------------------------------
!  store external variables into one-dimensional array
!-----------------------------------------------------------------------
!
#if defined QLOBERL && defined FORCE_DWLW
  WRITE(io_stdout,*) 'warning: options qloberl and force_dwlw make no sense'
#endif
#if defined BULK_KARA && defined DASILVA
  WRITE(io_stdout,*) 'warning: options BULK_KARA .AND. DASILVA make no sense'
#endif
     
!  initialize with zero

    hice1(:)=0.0
    alb1 (:)=0.0
    a21  (:)=0.0
    t1   (:)=0.0
    tair1(:)=0.0
    td1  (:)=0.0
    acl1 (:)=0.0
    pa1  (:)=0.0
    ug1  (:)=0.0
    sln1 (:)=0.0
    xlat (:)=0.0

    fhs(:,:)=0.0
    fhb(:,:)=0.0
    qsw(:,:)=0.0
    qlw(:,:)=0.0
    qse(:,:)=0.0
    qla(:,:)=0.0
  
    k=0
    DO i=1,ie
      DO j=1,je
        IF((om(i,j) >= 0.5) .AND. (a(i,j,lrhs) >= 1.e-6)) THEN
    
          k=k+1
          hice1(k) = hice(i,j)
          alb1 (k) = alb (i,j)
          a21  (k) = a2  (i,j)
          t1   (k) = t   (i,j) + tmelt   ! Change surface temp units from C to K
          tair1(k) = tair(i,j) + tmelt   ! Change air temp units from C to K
          td1  (k) = td  (i,j) - tmelt   ! Change dew point units from K to C
          acl1 (k) = acl (i,j)
          pa1  (k) = pa  (i,j)
          ug1  (k) = MAX(ug(i,j),1.e-6)
          sln1 (k) = sln (i,j)
          xlat (k) = MIN(ABS(alat(i,j)),60.)

        ENDIF
      ENDDO
    ENDDO
  
    IF (k == 0) RETURN
  
  
!-----------------------------------------------------------------------
!  Compute water vapor pressure according to "Buck Research Manual
!  (1996); updated from Buck, A. L., New equations for computing vapor pressure 
!  and enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981"
!  in 2m (esta) and at ice surface(esti)

    esta(1:k)  = 611.21 * EXP( (18.729-td1(1:k)/227.3) &
               * td1(1:k) / (td1(1:k)+257.87) )
    esti(1:k)  = 611.15 * EXP( (23.036-(MIN(t1(1:k),tmelt)-tmelt)/333.7) &
               * (MIN(t1(1:k),tmelt)-tmelt) / (MAX(t1(1:k),200.)-tmelt+279.82))
  
!-----------------------------------------------------------------------
  
  
#if defined BULK_KARA

    sphumida(1:k) = 0.62197 * esta(1:k) / (pa1(1:k)-0.378*esta(1:k))
    sphumido(1:k) = 0.62197 * esti(1:k) / (pa1(1:k)-0.378*esti(1:k))
    rhoair_k(1:k) = pa1(1:k) / ( rgas*tair1(1:k) * (1.0+0.61*sphumida(1:k)) )
    ugg(1:k)   = MAX(2.5, MIN(32.5,ug1(1:k)) )
    dragl(1:k) = MAX(0.5e-3,&
                 1.e-3*(0.8195+0.0506*ugg(1:k)-0.0009*ugg(1:k)*ugg(1:k)&
                 + (-0.0154+0.5698/ugg(1:k)-0.6743/(ugg(1:k)*ugg(1:k)))*&
                 (t1(1:k)-tair1(1:k)) ) )
    dragl(1:k) = MIN(dragl(1:k),3.0E-3)
    drags(1:k) = 0.96*dragl(1:k)

#elif defined DASILVA
    
    sphumida(1:k) = 0.622*esta(1:k)/(1.e5-0.378*esta(1:k))
    sphumido(1:k) = 0.622*esti(1:k)/(1.e5-0.378*esti(1:k))
    CALL vardrag(dragl, drags, k, td1+tmelt, t1, sphumida, sphumido, ug1)

#else
    
    dragl(1:k) = 1. 
    drags(1:k) = 1. 
    sphumida(1:k) = 0. 
    sphumido(1:k) = 0. 

#endif

!-----------------------------------------------------------------------
!  Compute cloudiness factor for downgoing longwave radiation
!               (Koch 1988; Beitr.Phys.Atmosph.,61(4),344-354);
!                absorbed surface radiation;
!          or    Berliand, 1952
!-----------------------------------------------------------------------
  
  
#ifndef FORCE_DWLW
  
#ifndef QLOBERL
    fakts(1:k) = 1. + 0.3*acl1(1:k)**2
#else
    fakts(1:k) = 1. - (0.5+0.4/90.*xlat(1:k)) * acl1(1:k)**2
#endif // QLOBERL

#endif // FORCE_DWLW


 !-----------------------------------------------------------------------
 !  Make first guess for surface temperature and heat balance
 !-----------------------------------------------------------------------
    stp(1:k)   =  t1(1:k)
!uwe   make albedo temperature dependent
!         flag=(0.5+sign(0.5,t1(n)-tmelt))
    flag(1:k)  =  1./ ( 1.+albtrans * (stp(1:k)-tmelt)**2 )
    alb1(1:k)  =  flag(1:k) * (a21(1:k)*albsnm + (1.-a21(1:k)) * albm)  &
                  + (1.-flag(1:k)) * (a21(1:k)*albsn+(1.-a21(1:k)) *albi)
    flai (1:k) = (1.-alb1(1:k))*sln1(1:k)

#ifndef QLOBERL
  
#ifdef FORCE_DWLW
    qlwin1(1:k) = acl1(1:k)
#else
    feu(1:k) = 0.601 + 5.95*1.e-7*esta(1:k)*exp(1500./tair1(1:k))
    qlwin1(1:k) = fakts(1:k)*feu(1:k)*d3*tair1(1:k)**4
#endif // FORCE_DWLW
  

  
#ifdef BULK_KARA
  qnet(1:k) = d3*stp(1:k)**4-flai(1:k)-qlwin1(1:k)                                  &
       - rhoair_k(1:k)*cpa*ug1(1:k)*(tair1(1:k)-stp(1:k))*drags(1:k)*fr_fac         &
       - rhoair_k(1:k)*subl*ug1(1:k)*(sphumida(1:k)-sphumido(1:k))*fr_fac           &
       * dragl(1:k)                                                                 &
       + (stp(1:k)-tb)/hice1(1:k)*con                                            
#else
  qnet(1:k) = d3*stp(1:k)**4-flai(1:k)-qlwin1(1:k)                                  &
       - d1*ug1(1:k)*(tair1(1:k)-stp(1:k))*drags(1:k)                               &
       - d2i*ug1(1:k)*(esta(1:k)-esti(1:k))*0.623/pa1(1:k)*dragl(1:k)               &
       + (stp(1:k)-tb)/hice1(1:k)*con
#endif // BULK_KARA
      
           
#else // QLOBERL

    feu(1:k)= 0.39-0.05*sqrt(esta(1:k)/100.)

#ifdef BULK_KARA
    qnet(1:k) = -flai(1:k)+(fakts(1:k)*feu(1:k)*d3*tair1(1:k)**4                 &
         +4.*d3*(tair1(1:k)**3)*(stp(1:k)-tair1(1:k)))                           &
         -rhoair_k(1:k)*cpa*ug1(1:k)*(tair1(1:k)-stp(1:k))*drags(1:k)*fr_fac     &
         -rhoair_k(1:k)*subl*ug1(1:k)*(sphumida(1:k)-sphumido(1:k))*fr_fac       &
         * dragl(1:k)                                                            &
         + (stp(1:k)-tb)/hice1(1:k)*con
#else
    qnet(1:k) = -flai(1:k)+(fakts(1:k)*feu(1:k)*d3*tair1(1:k)**4                 &
              +4.*d3*(tair1(1:k)**3)*(stp(1:k)-tair1(1:k)))                      &
            - d1*ug1(1:k)*(tair1(1:k)-stp(1:k))*drags(1:k)                       &
            - d2i*ug1(1:k)*(esta(1:k)-esti(1:k))*0.623/pa1(1:k)*dragl(1:k)       &
            + (stp(1:k)-tb)/hice1(1:k)*con
#endif // BULK_KARA

            
#endif // QLOBERL

    t1(1:k) = t1(1:k)+1.

!-----------------------------------------------------------------------
!  Calculate the surface temperature (start of iteration procedure)
!-----------------------------------------------------------------------

    DO iter=1,imax
       esti(1:k)     =  611.15*exp((23.036-(min(t1(1:k),tmelt)-tmelt)/333.7)*  &
                        (min(t1(1:k),tmelt)-tmelt)/(max(t1(1:k),200.)-tmelt    &
                        +279.82))
#ifdef BULK_KARA
     sphumido(1:k) = 0.62197*esti(1:k)/(1.e5-0.378*esti(1:k))
#endif

     stpp (1:k) = stp(1:k)
     qnetp  (1:k) = qnet(1:k)
     stp  (1:k) = t1(1:k)
     flag (1:k) = 1./(1.+albtrans*(stp(1:k)-tmelt)**2)
     flai (1:k) = (1.-alb1(1:k))*sln1(1:k)


#ifndef QLOBERL
 
#ifdef FORCE_DWLW
     qlwin1  (1:k) = acl1(1:k)
#else
     feu  (1:k) = 0.601+5.95*1.e-7*esta(1:k)*exp(1500./tair1(1:k))
     qlwin1  (1:k) = fakts(1:k)*feu(1:k)*d3*tair1(1:k)**4
#endif // FORCE_DWLW
  

  
#ifdef BULK_KARA

     qnet(1:k) = d3*stp(1:k)**4-flai(1:k)-qlwin1                                &
          - rhoair_k(1:k)*cpa*ug1(1:k)*(tair1(1:k) -stp(1:k) )                  &
          * drags(1:k)*fr_fac-rhoair_k(1:k)*subl*ug1(1:k)                       &
          * (sphumida(1:k)-sphumido(1:k))*fr_fac*dragl(1:k)                     &
          + (stp(1:k)-tb)/hice1(1:k)*con

#else 

     qnet(1:k) = d3*stp(1:k)**4-flai(1:k)-qlwin1                                &
          - d1*ug1(1:k)*(tair1(1:k)-stp(1:k) )*drags(1:k)                       &
          - d2i*ug1(1:k)*(esta(1:k)-esti(1:k))*0.623/pa1(1:k)                   &
          * dragl(1:k)                                                          &
          + (stp(1:k)-tb)/hice1(1:k)*con

#endif // BULK_KARA
  


#else // QLOBERL

     feu(1:k) = 0.39-0.05*sqrt(esta(1:k)/100.)

#ifdef BULK_KARA

     qnet(1:k) = -flai(1:k)+(fakts(1:k)*feu(1:k)*d3*tair1(1:k)**4              &
          + 4.*d3*(tair1(1:k)**3)*(stp(1:k)-tair1(1:k)))                       &
          - rhoair_k(1:k)*cpa*ug1(1:k)*(tair1(1:k)-stp(1:k))*drags(1:k)        &
          * fr_fac-rhoair_k(1:k)*subl*ug1(1:k)*(sphumida(1:k)                  &
          - sphumido(1:k))*fr_fac*dragl(1:k)                                   &
          + (stp(1:k)-tb)/hice1(1:k)*con
#else 

     qnet(1:k) = -flai(1:k)+(fakts(1:k)*feu(1:k)*d3*tair1(1:k)**4              &
          + 4.*d3*(tair1(1:k)**3)*(stp(1:k)-tair1(1:k)))                       &
          - d1*ug1(1:k)*(tair1(1:k)-stp(1:k))*drags(1:k)                       &
          - d2i*ug1(1:k)*(esta(1:k)-esti(1:k))* 0.623/pa1(1:k)                 &
          * dragl(1:k)                                                         &
          + (stp(1:k)-tb)/hice1(1:k)*con

#endif // BULK_KARA


#endif // QLOBERL
  
     fdiff (1:k)   = qnet(1:k)-qnetp(1:k)
     t1    (1:k)   = max(stp(1:k)-(stp(1:k)-stpp(1:k))*qnet(1:k)/              &
                     max(abs(fdiff(1:k)),1.e-10)*sign(1.,fdiff(1:k)),100.)
  
  END DO
  
!-----------------------------------------------------------------------
!  Calculate growth rates with updated heat balance equation:
!            fhs is melting rate at the snow/ice surface
!            fhb is melting/growth rate at the bottom of the ice.
!-----------------------------------------------------------------------

 
    t1  (1:k) = min(t1(1:k),tmelt)

#ifdef QLOBERL 
    feu (1:k) = 0.39-0.05*sqrt(esta(1:k)/100.)
#else
    feu (1:k) = 0.601+5.95*1.e-7*esta(1:k)*exp(1500./tair1(1:k))
#endif

!DEL         A1       = (0.5+SIGN(0.5,T1(1:k)-TMELT))
  !

#ifdef BULK_KARA
    qse1(1:k) = rhoair_k(1:k)*cpa*ug1(1:k)*(tair1(1:k)-t1(1:k))            &
                * drags(1:k)*fr_fac             
    qla1(1:k) = rhoair_k(1:k)*subl*ug1(1:k)*(sphumida(1:k)-sphumido(1:k))  &
                * fr_fac*dragl(1:k)
#else
    qse1(1:k) = d1*ug1(1:k)*(tair1(1:k)-t1(1:k))*drags(1:k)
    qla1(1:k) = d2i*ug1(1:k)*(esta(1:k)-esti(1:k))*0.623/pa1(1:k)*dragl(1:k)
#endif /*BULK_KARA*/
  
#ifdef QLOBERL 
    qlw_out1   (1:k) = 0.
#else
    qlw_out1   (1:k) = d3*t1(1:k)**4
#endif // QLOBERL

    qsw1   (1:k) =  (1.-alb1(1:k))*sln1(1:k)

#ifdef FORCE_DWLW
    qlw_in1   (1:k) =  acl1(1:k)
#else // FORCE_DWLW
#ifdef QLOBERL 
    qlw_in1   (1:k) = - (fakts(1:k)*feu(1:k)*d3*tair1(1:k)**4               &  
                 + 4.*d3*(tair1(1:k)**3)*(t1(1:k)-tair1(1:k)))
#else
    qlw_in1   (1:k) = fakts(1:k)*feu(1:k)*d3*tair1(1:k)**4
#endif // QLOBERL
#endif // FORCE_DWLW
  
    fhs1 (1:k) = (qlw_out1(1:k)-qsw1(1:k)-qlw_in1(1:k)-qse1(1:k)-qla1(1:k)           &
                 - (tb-t1(1:k))/hice1(1:k)*con)/clb
    fhb1 (1:k) =  (tb-t1(1:k))/hice1(1:k)*con /clb
  
!-----------------------------------------------------------------------
!  Unscramble for two-dimensional field
!-----------------------------------------------------------------------

    k=0
    DO  i=1,IE
      DO  j=1,JE
        IF((om(i,j) >= 0.5) .AND. (a (i,j,lrhs) >= 1.e-6)) THEN
          k        = k+1
          fhs(i,j) = fhs1 (k)
          fhb(i,j) = fhb1 (k)
          t  (i,j) = t1   (k) - tmelt
          qsw(i,j) = qsw1   (k)
          qlw(i,j) = qlw_in1   (k) - qlw_out1(k)
          qse(i,j) = qse1(k)
          qla(i,j) = qla1(k)
        ENDIF
      ENDDO
    ENDDO
    RETURN
 
 END SUBROUTINE budget
 


  
 SUBROUTINE obudget(alat,qt,tair,td,acl,pa,ug,sln,fo,om,           &
       qsw,qlw,qse,qla)

!=======================================================================
!  programmed by:
!  --------------
!     a.stoessel                 mpi, hamburg                      1989
!
!  modified
!  --------
!  s.legutke                     *dkrz*          13.8.95
!    - net atmospheric heat flux used for surface temp. > tfrez
!
!  helmuth 19.04.2000
!  - longwave-flux from berliand(1952)
!
!
!  purpose:
!  --------
!     calculates growth rates of new ice in the ice free part of a grid
!     cell with bulk formulas or net atm. heat flux
!
!  method:
!  -------
!     heat budget equation for open water
!
!  interface:
!  ----------
!     input:
!     *qt*       sea surface=oml temperature [deg c]
!     *tair*     2m air temperatures         [deg c]
!uwe     *sair*     surface temperatures        [k]
!     *td*       dew point temperatures      [k]
!     *acl*      fractional cloud cover      [frac.]
!     for cpp option force_dwlw acl contains downward longwave radiation!!
!     *pa*       atmospheric surface presure [pa]
!     *ug*       2 m wind speed              [m/sec]
!     *sln*      incoming surface solar rad. [w/m**2]
!     *om*       land/sea mask               [0/1]
!     *a*        ice compactness             [frac.]
!uwe     *aohflx*   atmospheric heat flux       [w/m**2]
!uwe     *aodflx*   heat flux correction        [w/m**2]
!
!     output:
!     *fo*       growth rate of thin ice     [m/sec]
!     *qsw*      absorbed solar radiation    [w/m**2]
!     *qlw*      outgoing longwave heat flux [w/m**2]
!     *qse*      sensible heat flux          [w/m**2]
!     *qla*      latent heat flux            [w/m**2]
!
!
!  externals:
!  ----------
!     vapor: calculates vapor pressure
!
!=======================================================================

    USE mo_param1, ONLY: ie,je
    USE mo_commoau1
    IMPLICIT NONE

    REAL,intent(in),DIMENSION(ie,je) :: alat,qt,sln,tair,td,acl,pa,ug,om
    REAL,intent(out),DIMENSION(ie,je) :: qsw,qlw,qse,qla,fo

    REAL, DIMENSION(ie*je) ::  tair1,td1,acl1,pa1,ug1,sln1,qse1,      &
         qla1,fh1,qt1,esta,estw,fakts,                              &
         qlw_out1,qsw1,qlw_in1,xlat,drags,dragl,sphumido,                         &
         sphumida, feu, rhoair_k, ugg

    REAL, PARAMETER ::         &
        rgas     = 287.1_wp,   & ! specific gas constant of dry air  [J/kg/K]
        cpa      = 1004.67_wp, & ! specific heat of dry air          [J/kg/K] 
        fr_fac   = 1.1925_wp     ! Frank Roeske's budget closing factor

    INTEGER :: i,j,k

 !-----------------------------------------------------------------------
 !  Select the grid cells and store them in selective 1-dimensional array
 !-----------------------------------------------------------------------
    
#if defined FORCE_DWLW && defined QLOBERL
    WRITE(io_stdout,*) 'warning: cpp options qloberl and force_dwlw together'
#endif /*FORCE_DWLW*/

#ifdef FORCE_DWLW
    d3=0.97*5.67e-8
#endif

    qlw_out1(:) = 0.0
    qlw_in1(:)  = 0.0
    qsw1(:)     = 0.0
    sln1(:)     = 0.0

    k=0
    DO i=1,ie
       DO j=1,je
          IF(om(i,j) >= 0.5) THEN

             k = k+1
             qt1   (k) = qt(i,j)+tmelt
             tair1 (k) = tair(i,j)+tmelt
             td1   (k) = td(i,j)
             acl1  (k) = acl(i,j)
             pa1   (k) = pa(i,j)
             ug1   (k) = MAX(ug(i,j),1.e-6)
             sln1  (k) = sln(i,j)
             xlat  (k) = MIN(ABS(alat(i,j)),60.)

          ENDIF
       ENDDO
    ENDDO

!-----------------------------------------------------------------------
!  Compute water vapor pressure according to "Buck Research Manual
!  (1996); updated from Buck, A. L., New equations for computing vapor pressure 
!  and enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981"
!  in 2m (esta) and at sea surface(estw)
!
       esta(1:k) = 611.21*EXP((18.729-(td1(1:k)-tmelt)/227.3)   &
                   * (td1(1:k)-tmelt)/(td1(1:k)-tmelt+257.87))
       estw(1:k) = 0.9815*611.21*EXP((18.729-(qt1(1:k)-tmelt)/227.3) &
                   * (qt1(1:k)-tmelt)/(qt1(1:k)-tmelt+257.87))


!-----------------------------------------------------------------------
!  compute cloudiness factor for downgoing longwave radiation
!               (koch 1988; beitr.phys.atmosph.,61(4),344-354);
!-----------------------------------------------------------------------


#if defined DASILVA

     sphumida(1:k) = 0.622*esta(1:k)/(1.e5-0.378*esta(1:k))
     sphumido(1:k) = 0.622*estw(1:k)/(1.e5-0.378*estw(1:k))
    CALL vardrag(dragl, drags, k, qt1, tair1, sphumida, sphumido, ug1)

#elif defined BULK_KARA

    sphumida(1:k) = 0.62197 * esta(1:k) / (pa1(1:k)-0.378*esta(1:k))
    sphumido(1:k) = 0.62197 * estw(1:k) / (pa1(1:k)-0.378*estw(1:k))
    rhoair_k(1:k) = pa1(1:k) / ( rgas*tair1(1:k) * (1.0+0.61*sphumida(1:k)) )
    ugg(1:k)   = MAX(2.5, MIN(32.5,ug1(1:k)) )
    dragl(1:k) = MAX(0.5e-3, 1.e-3*(0.8195+0.0506*ugg(1:k)-0.0009*ugg(1:k) &
                 * ugg(1:k) + (-0.0154+0.5698/ugg(1:k)-0.6743/ugg(1:k)**2) &
                 * (qt1(1:k)-tair1(1:k)) ) )
    dragl(1:k) = MIN(dragl(1:k),3.0E-3)
    drags(1:k) = 0.96*dragl(1:k)

#else

    dragl(1:k) = 1. 
    drags(1:k) = 1. 
    sphumida(1:k) = 0. 
    sphumido(1:k) = 0. 

#endif


#ifndef FORCE_DWLW

#ifndef QLOBERL
    fakts(1:k) = 1. + 0.3*acl1(1:k)**2
#else
    fakts(1:k) = 1. - (0.5+0.4/90.*xlat(1:k)) * acl1(1:k)**2
#endif // QLOBERL

#endif // FORCE_DWLW

!-----------------------------------------------------------------------
!  CALCULATE HEAT FLUXES AND GROWTH RATES:
!            USE HEAT BALANCE EQUATION WHERE THE ATMOSPHERE HAD ICE
!            USE NET ATMOSPHERIC HEAT FLUX + FLUX CORRECTION ELSEWHERE
!-----------------------------------------------------------------------

#ifndef FORCE_DWLW

#ifndef QLOBERL
     feu(1:k)      = 0.605+5.95*1.e-7*esta(1:k)*EXP(1500./tair1(1:k))
#else
     feu(1:k)      = 0.39-0.05*SQRT(esta(1:k)/100.)
#endif //QLOBERL

#endif //FORCE_DWLW


#ifdef BULK_KARA
     rhoair_k(1:k) = pa1(1:k)/(rgas*tair1(1:k)*(1.0+0.61*sphumida(1:k)) )
     qse1(1:k)    = rhoair_k(1:k)*cpa*ug1(1:k)*(tair1 (1:k)-qt1(1:k) )       &
                     * drags(1:k)*fr_fac
     qla1(1:k)    = (2.501-0.00237*(qt1(1:k)-tmelt))*1.e+6 *rhoair_k(1:k)     &
                     * ug1(1:k)*(sphumida(1:k)-sphumido(1:k))*dragl(1:k)*fr_fac
#else
     qse1(1:k)   = d1 *ug1(1:k)*(tair1 (1:k)-qt1(1:k) ) *drags(1:k)
     qla1(1:k)   = d2w*ug1(1:k)*(esta(1:k)-estw(1:k))*.623/pa1(1:k)*dragl(1:k)
#endif //BULK_KARA


#ifndef QLOBERL
     qlw_out1(1:k)      = d3*qt1(1:k)**4
#else
     qlw_out1(1:k)      = 0.
#endif


     qsw1(1:k)      = (1.-albw)*sln1(1:k)


#ifdef FORCE_DWLW

     qlw_in1(1:k)      = acl1(1:k)

#else
  
#ifndef QLOBERL
     qlw_in1(1:k)      = fakts(1:k)*feu(1:k)*d3*tair1(1:k)**4
#else
     qlw_in1(1:k)      = -(fakts(1:k)*feu(1:k)*d3*tair1(1:k)**4   &
                    +4.*d3*(tair1(1:k)**3)*(qt1(1:k)-tair1(1:k))) 
                    
#endif /*QLOBERL*/

#endif /*FORCE_DWLW*/

!             sw = qsw1
!             lw = qlw_in1-qlw_out1
!             se = qse1
!             la = qla1

    fh1(1:k)    = (qlw_out1(1:k)-qsw1(1:k)-qlw_in1(1:k)-qse1(1:k)-qla1(1:k))/clb


!-----------------------------------------------------------------------
!  unscramble into two-dimensional field
!-----------------------------------------------------------------------

    k=0
    qsw(:,:) = 0.0
    qlw(:,:) = 0.0
    qse(:,:) = 0.0
    qla(:,:) = 0.0
    fo (:,:) = 0.0

    DO i=1,ie
       DO j=1,je
          IF(om(i,j) >= 0.5) THEN

             k=k+1
             fo (i,j) = fh1(k)
             qsw(i,j) = qsw1(k)
             qlw(i,j) = qlw_in1(k)-qlw_out1(k)
             qse(i,j) = qse1(k)
             qla(i,j) = qla1(k)
       
       ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE obudget

END MODULE mo_omip
