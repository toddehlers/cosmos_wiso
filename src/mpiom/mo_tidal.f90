      MODULE mo_tidal

      USE mo_param1

      IMPLICIT NONE

      integer :: mmccdt

      REAL, ALLOCATABLE :: acc(:,:),ass(:,:),acs(:,:)
      REAL, ALLOCATABLE :: alatr(:,:),alonr(:,:),tipoto(:,:)
      REAL, ALLOCATABLE :: alatru(:,:),alonru(:,:)
      REAL, ALLOCATABLE :: alatrv(:,:),alonrv(:,:)
      REAL, ALLOCATABLE :: colato(:,:),silato(:,:)
      REAL, ALLOCATABLE :: colatou(:,:),silatou(:,:)
      REAL, ALLOCATABLE :: colatov(:,:),silatov(:,:)
      REAL, ALLOCATABLE :: colono(:,:),silono(:,:)
      REAL, ALLOCATABLE :: colonou(:,:),silonou(:,:)
      REAL, ALLOCATABLE :: colonov(:,:),silonov(:,:)
      REAL, ALLOCATABLE :: boxare(:,:),boxaro(:,:)

      CONTAINS

      SUBROUTINE alloc_mem_tidal

      ALLOCATE(acc(ie,je),ass(ie,je),acs(ie,je))
      ALLOCATE(alatr(ie,je),alonr(ie,je),tipoto(ie,je))
      ALLOCATE(alatru(ie,je),alonru(ie,je))
      ALLOCATE(alatrv(ie,je),alonrv(ie,je))
      ALLOCATE(colato(ie,je),silato(ie,je))
      ALLOCATE(colatou(ie,je),silatou(ie,je))
      ALLOCATE(colatov(ie,je),silatov(ie,je))
      ALLOCATE(colono(ie,je),silono(ie,je))
      ALLOCATE(colonou(ie,je),silonou(ie,je))
      ALLOCATE(colonov(ie,je),silonov(ie,je))
      ALLOCATE(boxare(ie,je),boxaro(ie,je))  

      END SUBROUTINE alloc_mem_tidal

      Subroutine foreph_ini

      Use MO_PARAM1
      Use MO_PARALLEL
      Use MO_COMMO1
      Use MO_UNITS


      IMPLICIT NONE

      integer :: ioff,joff,i,j,jcc,moph

!     calc time index for tidal model 
      mmccdt = 0
      jcc=0.
      moph=0.

      CALL eph(lyear1,lmont1,jcc,moph)

      mmccdt=(jcc+moph)*NINT(86400./DT)
        
      if (icontro.ne.0) then
      WRITE(0,*)'tidal: phase relative to 2000 :'    &
      ,'year= ',lyear1                                &
      ,'month= ',lmont1                              &
      ,'yearoff= ',jcc                               &
      ,' monoff= ',moph                              &
      ,'mmccdt= ',mmccdt
      endif

!      ioff = 2*p_ioff
!      joff = 2*p_joff

!      GILA(:,:) = GILA_G(ioff+1:ioff+2*ie,joff+1:joff+2*je)
!      GIPH(:,:) = GIPH_G(ioff+1:ioff+2*ie,joff+1:joff+2*je)

      DO J=1,JE
         DO I=1,IE
            alatr(i,j)=giph(2*i,2*j)
            colato(i,j)=cos(giph(2*i,2*j))
            silato(i,j)=sin(giph(2*i,2*j))

            alonr(i,j)=gila(2*i,2*j)
            colono(I,J)=COS(GILA(2*I,2*J))!:: cxe new lines begin
            silono(I,J)=SIN(GILA(2*I,2*J))

            acc(i,j) = colono(I,J)*colono(I,J)
            ass(i,j) = silono(I,J)*silono(I,J)
            acs(i,j) = colono(I,J)*silono(I,J)

            boxaro(i,j) = dlxu(i,j)*dlyu(i,j)
            boxare(i,j) = dlxv(i,j)*dlyv(i,j)

         end do
      end do

      DO J=2,JE-1
         DO I=2,IE-1
            
            alatru(i,j)=giph(2*i+1,2*j)
            colatou(i,j)=cos(giph(2*i+1,2*j))
            silatou(i,j)=sin(giph(2*i+1,2*j))

            alatrv(i,j)=giph(2*i,2*j+1)
            colatov(i,j)=cos(giph(2*i,2*j+1))
            silatov(i,j)=sin(giph(2*i,2*j+1))

            alonru(i,j)=gila(2*i+1,2*j)
            colonou(I,J)=COS(GILA(2*I+1,2*J))
            silonou(I,J)=SIN(GILA(2*I+1,2*J))

            alonrv(i,j)=gila(2*i,2*j+1)
            colonov(I,J)=COS(GILA(2*I,2*J+1))
            silonov(I,J)=SIN(GILA(2*I,2*J+1))

         end do
      end do

      call bounds_exch('u+',alatru,'foreph 1') 
      call bounds_exch('u+',colatou,'foreph 2') 
      call bounds_exch('u+',silatou,'foreph 3') 
      call bounds_exch('v+',alatrv,'foreph 4') 
      call bounds_exch('v+',colatov,'foreph 5') 
      call bounds_exch('v+',silatov,'foreph 6') 


      end Subroutine foreph_ini

      Subroutine foreph
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     This program calculates the tidal potential due to sun/moon.
!     Xueen Chen. Jun 22th, 2005.
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Use MO_PARAM1
      Use MO_PARALLEL
      Use MO_COMMO1
      Use MO_UNITS


      IMPLICIT NONE


      Real :: dres(3,2),crim3,rkomp,erdrad,rekts,dekls
      real :: cris3, rektm,deklm,deklm2,dekls2,sidm,sidmq
      real :: rkosp,codm,codmq,sids,sidsq,cods,codsq,sidm2
      real :: sids2,sip2,cop2,argp,hamp,hasp
      integer :: i,j


      mmccdt = mmccdt + 1

      rkomp = -4.113e-07
      rkosp = 0.46051*rkomp
      erdrad=6371000.

      Call ephvsop87(dres)
      rekts=dres(1,1)
      dekls=dres(2,1)
      cris3=dres(3,1)

      rektm=dres(1,2)
      deklm=dres(2,2)
      crim3=dres(3,2)

      if (icontro.ne.0) then
         Print*, 'ephvsop87(dres) => dres'
         Print*, dres
         Print*,'in potential tipoto: rkomp*erdrad,rkosp*erdrad'
         Print*, rkomp*erdrad,rkosp*erdrad
         !      Print*, 'alatr(20,54) & alonr(20,54): ',alatr(20,54),alonr(20,54)
      endif

      Do j=1,je
!!      Do j=2,je1

      deklm2 = deklm*2.
      dekls2 = dekls*2.
      sidm   = Sin(deklm)
      sidmq  = sidm*sidm
      codm   = Cos(deklm)
      codmq  = codm*codm
      sids   = Sin(dekls)
      sidsq  = sids*sids
      cods   = Cos(dekls)
      codsq  = cods*cods
      sidm2=Sin(deklm2)
      sids2=Sin(dekls2)
 
      DO i=1,ie

      sip2=2.*silato(i,j)*colato(i,j)
      cop2=colato(i,j)**2-silato(i,j)**2
      argp=alonr(i,j)

      hamp = rektm+argp
      hasp = rekts+argp

      tipoto(i,j)=erdrad*rkomp*crim3*(3.*(silato(i,j)**2-1./3.)         &
     &          *(sidmq-1./3.)+SIN(2.*alatr(i,j))*sidm2*COS(hamp)+      &
     & colato(i,j)**2*codmq*COS(2.*hamp))                               &
     &+erdrad*rkosp*cris3*(3.*(silato(i,j)**2-1./3.)*(sidsq-1./3.)      &
     &+SIN(2.*alatr(i,j))*sids2*COS(hasp)+                              &
     & colato(i,j)**2*codsq*COS(2.*hasp))          

      END DO
      END DO

       Call bounds_exch('p',tipoto,'ocpheme 1')
       Call bounds_exch('u+',dlxu,'foreph 7')
       Call bounds_exch('v+',dlyv,'foreph 8')
       Call bounds_exch('u+',amsuo,'foreph 9')
       Call bounds_exch('v+',amsue,'foreph 10')

      End Subroutine foreph

      Subroutine tipouv

      Use MO_PARAM1
      Use MO_PARALLEL
      Use MO_COMMO1
      Use MO_UNITS


      IMPLICIT NONE

      real, parameter :: sal=0.69   !Earth Tides ( maik thomas, emr  pers. comm )  

      integer :: i,j,k


      Do k=1,ke
      do j=2,je1
      do i=2,ie1
!      if(dlxu(i,j)*dlyv(i,j).lt.1.)              &
!     &  write(0,*)i,j,dlxu(i,j),dlyv(i,j)                &
!     &    ,nint(amsuo(i,j,1)),nint(amsue(i,j,1))
         uoo(i,j,k) = uoo(i,j,k)+sal*(amsuo(i,j,k)*dt*                 &
     &   (tipoto(i,j)-tipoto(i+1,j))/dlxu(i,j))
         voe(i,j,k) = voe(i,j,k)+sal*(amsue(i,j,k)*dt*                 &
     &   (tipoto(i,j+1)-tipoto(i,j))/dlyv(i,j))
      end do
      end do
      end do

      Call bounds_exch('u',uoo,'ocpheme 2')
      Call bounds_exch('v',voe,'ocpheme 3')

!!$      Print*, 'Tidal potential induced zonal and meridional velocities'
!!$
!!$      Write(6,6118)(1000*(tipoto(i,54)                                  &
!!$     &       -tipoto(i+1,54))/dlxu(i,54)*amsuo(i,54,1)*dt,i = 10,120,10)
!!$ 6118 Format('Tpotent  u?',12f9.3)
!!$      Write(6,6119)(1000*(tipoto(i,55)                                  &
!!$     &       -tipoto(i  ,54))/dlyv(i,54)*amsue(i,54,1)*dt,i = 10,120,10)
!!$ 6119 Format('Tpotent  v?',12f9.3)
!!$      Write(6,6120)(uoo(i,54,1),i=10,120,10)
!!$ 6120 Format('TotalVel u?',12f9.3)
!!$      Write(6,6121)(voe(i,54,1),i=10,120,10)
!!$ 6121 Format('TotalVel v?',12f9.3)
!!$      Write(6,6123)(weto(i,54,1),i=10,120,10)
!!$ 6123 Format('Mask  weto?',12f9.3)
!!$      Write(6,6122)(tipoto(i,54),i=10,120,10)
!!$ 6122 Format('TPotential?',12f9.3)
!!$      Write(6,6124)(ZO(i,54),i=10,120,10)
!!$ 6124 Format('Sea levels?',12f9.3)


      Return
      End Subroutine tipouv


      Subroutine ephvsop87(res2)
!     This subroutine calcuates the ephemeredes of sun/moon
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Use MO_PARAM1
      Use MO_MPI
      Use MO_PARALLEL
      Use MO_COMMO1
      Use MO_UNITS
!     IMPLICIT NONE

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Ann-Margrit Hellmich: 07/14/2003
! Xueen Chen: 02.08.2005
!  calculation of
!  right ascension, declination and geocentric distance
!  for Sun, Moon and major planets
!     +
!  optional preparation for calculation of potential 
!  according to ephaufb.f by Maik Thomas
! 
      Integer fnut
      Real pi2,pic,T,sidt,ecl,nutob,nutl,res(3,2)
      Real res2(3,2)

! C pi, ...
      pi2=dacos(-1.d0)*2
      pic=dacos(-1.d0)/Dble(180.d0)

! D cxe Transformation of Time into fractional julian centuries T
! cxe i.e. julian century T is of the beginning to each timestep
! we set it as 99999.99

! cxe Julian date 2451545 means 01.01.2000 at 12:00 UT
! cxe Why should 01.01.2000 12:00 UT be refer time ???
      T=(9999.99+(mmccdt-1)*dt/86400.-2451545.d0)/36525.
!      tidtim=t

!      Print*,'Subroutine ephvsop: julian century T at start of mmccdt => ',T,mmccdt
! C corresponding siderical time Greenwich
      Call sidt2(pic,pi2,T,sidt)

! C set flags 
!  set fnut (perform nutation -> 1; don't -> 0)
      fnut=0

! C obliquity of the ecliptic
      Call obliq(fnut,pic,pi2,T,ecl,nutob,nutl)

! C Sun
      Call sun_n(fnut,pic,pi2,T,ecl,nutl,res)

! C Moon      
      Call moon(fnut,pic,pi2,T,ecl,nutl,res)

! C modifications in preparaion of calculation of potentials
      Call aufb2(sidt,res,res2)
! 
      End Subroutine ephvsop87 

      Subroutine time2(iiyr,iimo,iida,iihr,iimin,iisec,T)
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
! author        : Ann-Margrit Hellmich
! 
! last modified : 03/17/2003
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
!  convertion of date and UTC to Julian Centuries since J2000 T
!  according to Duffett, 1990
! 
       real :: daf,t,ajd
       integer :: iida, iimo, iiyr,iisec,iihr,iimin,ia,ib
! 
      daf=iida+iihr/Dble(24.d0)+iimin/Dble(3600.d0)+iisec/Dble(86400.d0)
      If(iimo.Le.2) Then
        iimo=iimo+12
        iiyr=iiyr-1
      Endif
      IA=Int(iiyr/100)
      IB=2-IA+IA/4
!  Julian Date
      AJD=Int(365.25*(iiyr+4716))+Int(30.6001*(iimo+1))+daf+IB-1524.5
!  Julian Centuries since J2000
      T=(AJD-2451545.0)/Dble(36525.0)
! 
      End Subroutine time2

      Subroutine sidt2(pic,pi2,T,sidt)
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
! author        : Ann-Margrit Hellmich
! 
! last modified : 03/17/2003
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
!  convertion of Julian Centuries since J2000 T to siderical time sidt
!  according to Duffett, 1990
! 
      Real pic,pi2,T,sidt,JD,T2,T3
! 
      T2=T*T
      T3=T2*T
!      
!  Julian days
      JD=T*36525.0+2451545.0
! 
!  mean siderial time of Greenwich in rad
      sidt=(280.46061837+360.98564736629*(JD-2451545.0)                 &
     &+0.000387933*T2-T3/Dble(38710000))*pic
!  between 0 and 2pi :
      If(sidt.Lt.0) Then
        Call negangle2(pi2,sidt)
      Endif
      If(sidt.Ge.pi2) Then
        Call langle2(pi2,sidt)
      Endif
! 
      End Subroutine sidt2

      Subroutine obliq(fnut,pic,pi2,T,ecl,nutob,nutl)
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
! author        : Ann-Margrit Hellmich
! 
! last modified : 03/17/2003
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
!  calculation of obliquity of ecliptic
!  according to Duffett, 1990
! 
      Integer fnut
      Real pic,pi2,T,ecl,nutob,nutl,A,B,C,T1,T2,T3
      Real L1,L2,D1,D2,M1,M2,N1,N2
! 
!  see page 57
! correction terms if requested
      If(fnut.Eq.1) Then
        T1=T+1
        T2=T1*T1
! 
        A=100.0021358*T1
        B=360.0*(A-Int(A))
        L1=279.6967+0.000303*T2+B
        L2=2.0*L1*pic
        A=1336.855231*T1
        B=360.0*(A-Int(A))
        D1=270.4342-0.001133*T2+B
        D2=2.0*D1*pic
        A=99.99736056*T1
        B=360.0*(A-Int(A))
        M1=(358.4758-0.00015*T2+B)*pic
        A=1325.552359*T1
        B=360.0*(A-Int(A))
        M2=(296.1046+0.009192*T2+B)*pic
        A=5.372616667*T1
        B=360.0*(A-Int(A))
        N1=(259.1833+0.002078*T2-B)*pic
        N2=2*N1
! correction term for nutation in longitude
        nutl=((-17.2327-0.01737*T1)*Sin(N1)                             &
     &+(-1.2729-0.00013*T1)*Sin(L2)+0.2088*Sin(N2)                      &
     &-0.2037*Sin(D2)+(0.1261-0.00031*T1)*Sin(M1)                       &
     &+0.0675*Sin(M2)-(0.0497-0.00012*T1)*Sin(L2+M1)                    &
     &-0.0342*Sin(D2-N1)-0.0261*Sin(D2+M2)                              &
     &+0.0214*Sin(L2-M1)-0.0149*Sin(L2-D2+M2)                           &
     &+0.0124*Sin(L2-N1)+0.0114*Sin(D2-M2))/Dble(3600.0)*pic
! correction term for nutation in obliquity of the ecliptic
        nutob=((9.21+0.00091*T1)*Cos(N1)                                &
     &+(0.5522-0.00029*T1)*Cos(L2)-0.0904*Cos(N2)                       &
     &+0.0884*Cos(D2)+0.0216*Cos(L2+M1)                                 &
     &+0.0183*Cos(D2-N1)+0.0113*Cos(D2+M2)                              &
     &-0.0093*Cos(L2-M1)-0.0066*Cos(L2-N1))/Dble(3600.0)*pic
! 
      Else
        nutob=0.0
        nutl=0.0
      Endif
! obliquity of the ecliptic
      T2=T*T
      T3=T2*T
      C=46.815*T+0.0006*T2-0.00181*T3
      ecl=(23.43929167-C/Dble(3600.0))*pic+nutob
! 
      End Subroutine obliq

      Subroutine Sun_n(fnut,pic,pi2,T,ecl,nutl,res)
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
! author        : Ann-Margrit Hellmich
! 
! last modified : 03/17/2003
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
!  calculation of position of the Sun
!  according to Duffett, 1990
! 
      Integer fnut
      Real pic,pi2,T,T1,T2,T3,ecl,A,B,nutl
      Real L,M1,EC,AT,AE
      Real A1,B1,C1,D1,E1,H1,D2,D3,S1,S2,S3,SW,X1,X2,res(3,2)
! 
!  see page 116
      T1=T+1
      T2=T1*T1
      T3=T2*T1
! 
      A=100.0021359*T1
      B=360.0*(A-Int(A))
      L=(279.69668+0.0003025*T2+B)*pic
      A=99.99736042*T1
      B=360.0*(A-Int(A))
      M1=(358.47583-0.00015*T2+0.0000033*T3+B)*pic
      EC=0.01675104-0.0000418*T1-0.000000126*T2
! 
!  true and eccentric anomaly in rad
      Call anomaly(pi2,M1,EC,AT,AE)
! 
!  various arguments in rad
      A=62.55209472*T1
      B=360.0*(A-Int(A))
      A1=(153.23+B)*pic
      A=125.1041894*T1
      B=360.0*(A-Int(A))
      B1=(216.57+B)*pic
      A=91.56766028*T1
      B=360.0*(A-Int(A))
      C1=(312.69+B)*pic
      A=1236.853095*T1
      B=360.0*(A-Int(A))
      D1=(350.74-0.00144*T2+B)*pic
      E1=(231.19+20.2*T1)*pic
      A=183.1353208*T1
      B=360.0*(A-Int(A))
      H1=(353.4+B)*pic
! 
      D2=(0.00134*Cos(A1)+0.00154*Cos(B1)+0.002*Cos(C1)                 &
     &+0.00179*Sin(D1)+0.00178*Sin(E1))*pic
      D3=0.00000543*Sin(A1)+0.00001575*Sin(B1)                          &
     &+0.00001627*Sin(C1)+0.00003076*Cos(D1)                            &
     &+0.00000927*Sin(H1)
! 
!  geocentric ecliptic coordinates of the Sun
      S1=AT+L-M1+D2
      If(fnut.Eq.1) Then
        S1=S1+nutl
      Endif
      If(S1.Lt.0) Call negangle2(pi2,S1)
      If(S1.Ge.pi2) Call langle2(pi2,S1)
      S2=0.0
      S3=1.0000002*(1.0-EC*Cos(AE))+D3
! 
!  geocentric equatorial coordinates of the Sun
      SW=-1
      Call eqecl(pic,pi2,S1,S2,X1,X2,ecl,SW)
      res(1,1)=X1
      res(2,1)=X2
      res(3,1)=S3
! 
      End Subroutine sun_n

      Subroutine moon(fnut,pic,pi2,T,ecl,nutl,res)
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
! author        : Ann-Margrit Hellmich
! 
! last modified : 03/18/2003
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
!  calculation of position of the Moon
!  according to Duffett, 1990
! 
      Integer fnut
      Real pic,pi2,T,T1,T2,T3,ecl,nutl,A,B,C,SW
      Real Q,M1,M2,M3,M4,M5,M6,ML,MS,MD,ME,MF,NA,S1,S2,S3,S4,E,E2
      Real L,G,W1,W2,PM,MO1,MO2,MO3,X1,X2,res(3,2)
! 
!  see page 157
      T1=T+1
      T2=T1*T1
      T3=T2*T1
! 
      Q=T1*36525
      M1=Q/27.32158213
      M1=360*(M1-Int(M1))
      M2=Q/365.2596407
      M2=360*(M2-Int(M2))
      M3=Q/27.55455094
      M3=360*(M3-Int(M3))
      M4=Q/29.53058868
      M4=360*(M4-Int(M4))
      M5=Q/27.21222039
      M5=360*(M5-Int(M5))
      M6=Q/6798.363307
      M6=360*(M6-Int(M6))
      ML=270.434164+M1-0.001133*T2+0.0000019*T3
      MS=358.475833+M2-0.00015*T2+0.0000033*T3
      MD=296.104608+M3+0.009192*T2+0.0000144*T3
      ME=350.737486+M4-0.001436*T2+0.0000019*T3
      MF=11.250889+M5-0.003211*T2-0.0000003*T3
      NA=(259.183275-M6+0.002078*T2+0.0000022*T3)*pic
      S2=Sin(NA)
      A=(51.2+20.2*T1)*pic
      S1=Sin(A)
      B=(346.56+132.87*T1-0.0091731*T2)*pic
      S3=0.003964*Sin(B)
      C=NA+(275.05-2.3*T1)*pic
      S4=Sin(C)
      ML=(ML+0.000233*S1+S3+0.001964*S2)*pic
      MS=(MS-0.001778*S1)*pic
      MD=(MD+0.000817*S1+S3+0.002541*S2)*pic
      MF=(MF+S3-0.024691*S2-0.004328*S4)*pic
      ME=(ME+0.002011*S1+S3+0.001964*S2)*pic
      E=1-0.002495*T1+0.00000752*T2
      E2=E*E
      
!  ecliptic longitude MO1
      L=6.28875*Sin(MD)+1.274018*Sin(2*ME-MD)                           &
     &+0.658309*Sin(2*ME)+0.213616*Sin(2*MD)                            &
     &-E*0.185596*Sin(MS)-0.114336*Sin(2*MF)                            &
     &+0.058793*Sin(2*(ME-MD))                                          &
     &+0.057212*E*Sin(2*ME-MS-MD)+0.05332*Sin(2*ME+MD)                  &
     &+0.045874*E*Sin(2*ME-MS)+0.041024*E*Sin(MD-MS)                    &
     &-0.034718*Sin(ME)-E*0.030465*Sin(MD+MS)                           &
     &+0.015326*Sin(2*(ME-MF))-0.012528*Sin(2*MF+MD)                    &
     &-0.01098*Sin(2*MF-MD)+0.010674*Sin(4*ME-MD)                       &
     &+0.010034*Sin(3*MD)+0.008548*Sin(4*ME-2*MD)                       &
     &-E*0.00791*Sin(MS-MD+2*ME)-E*0.006783*Sin(2*ME+MS)                &
     &+0.005162*Sin(MD-ME)+E*0.005*Sin(ME+MS)                           &
     &+0.003862*Sin(4*ME)+E*0.004049*Sin(MD-MS+2*ME)                    &
     &+0.003996*Sin(2*(MD+ME))+0.003665*Sin(2*ME-3*MD)                  &
     &+E*0.002695*Sin(2*MD-MS)+0.002602*Sin(MD-2*(MF+ME))               &
     &+E*0.002396*Sin(2*(ME-MD)-MS)-0.002349*Sin(ME+MD)                 &
     &+E2*0.002249*Sin(2*(ME-MS))-E*0.002125*Sin(MS+2*MD)               &
     &-E2*0.002079*Sin(2*MS)+E2*0.002059*Sin(2*(ME-MS)-MD)              &
     &-0.001773*Sin(2*(ME-MF)+MD)-0.001595*Sin(2*(ME+MF))               &
     &+E*0.00122*Sin(4*ME-MS-MD)-0.00111*Sin(2*(MD+MF))                 &
     &+0.000892*Sin(MD-3*ME)-E*0.000811*Sin(MS+MD+2*ME)                 &
     &+E*0.000761*Sin(4*ME-MS-2*MD)                                     &
     &+E2*0.000704*Sin(MD-2*(MS+ME))                                    &
     &+E*0.000693*Sin(MS-2*(MD-ME))                                     &
     &+E*0.000598*Sin(2*(ME-MF)-MS)                                     &
     &+0.00055*Sin(MD+4*ME)+0.000538*Sin(4*MD)                          &
     &+E*0.000521*Sin(4*ME-MS)+0.000486*Sin(2*MD-ME)                    &
     &+E2*0.000717*Sin(MD-2*MS)
      MO1=ML+L*pic
      If(fnut.Eq.1) Then
        MO1=MO1+nutl
      Endif
      If(MO1.Lt.0) Call negangle2(pi2,MO1)
      If(MO1.Ge.pi2) Call langle2(pi2,MO1)

!  ecliptic latitude MO2
      G=5.128189*Sin(MF)+0.280606*Sin(MD+MF)                            &
     &+0.277693*Sin(MD-MF)+0.173238*Sin(2*ME-MF)                        &
     &+0.055413*Sin(2*ME+MF-MD)+0.046272*Sin(2*ME-MF-MD)                &
     &+0.032573*Sin(2*ME+MF)+0.017198*Sin(2*MD+MF)                      &
     &+0.009267*Sin(2*ME-MF+MD)+0.008823*Sin(2*MD-MF)                   &
     &+E*0.008247*Sin(2*ME-MS-MF)+0.004323*Sin(2*(ME+MD)-MF)            &
     &+0.0042*Sin(2*ME+MD+MF)+E*0.003372*Sin(MF-MS-2*ME)                &
     &+E*0.002472*Sin(2*ME-MD+MF-MS)                                    &
     &+E*0.002222*Sin(2*ME+MF-MS)                                       &
     &+E*0.002072*Sin(2*ME-MD-MF-MS)                                    &
     &+E*0.001877*Sin(MF-MS+MD)+0.001828*Sin(4*ME-MD-MF)                &
     &-E*0.001803*Sin(MS+MF)-0.00175*Sin(3*MF)                          &
     &+E*0.00157*Sin(MD-MF-MS)-0.001487*Sin(ME+MF)                      &
     &-E*0.001481*Sin(MF+MS+MD)+E*0.001417*Sin(MF-MS-MD)                &
     &+E*0.00135*Sin(MF-MS)+0.00133*Sin(MF-ME)                          &
     &+0.001106*Sin(MF+3*MD)+0.00102*Sin(4*ME-MF)                       &
     &+0.000833*Sin(MF+4*ME-MD)+0.000781*Sin(MD-3*MF)                   &
     &+0.00067*Sin(MF+3*ME-2*MD)+0.000606*Sin(2*ME-3*MF)                &
     &+0.000597*Sin(2*(ME+MD)-MF)                                       &
     &+E*0.000492*Sin(2*ME+MD-MS-MF)+0.00045*Sin(2*(MD-ME)-MF)          &
     &+0.000439*Sin(3*ME-MF)+0.000423*Sin(MF+2*(ME+MD))                 &
     &+0.000422*Sin(2*ME-3*MD-MF)-E*0.000367*Sin(MF+MS+2*ME-MD)         &
     &-E*0.000353*Sin(MF+MS+2*ME)+0.000331*Sin(MF+4*ME)                 &
     &+E*0.000317*Sin(2*ME+MD-MS+MF)                                    &
     &+E2*0.000306*Sin(2*(ME-MS)-MF)-0.000283*Sin(MD+3*MF)
      W1=0.0004664*Cos(NA)
      W2=0.0000754*Cos(C)
      MO2=G*pic*(1.0-W1-W2)

!  horizontal parallax PM
      PM=0.950724+0.051818*Cos(MD)+0.009531*Cos(2*ME-MD)                &
     &+0.007843*Cos(2*ME)+0.002824*Cos(2*MD)                            &
     &+0.000857*Cos(2*ME+MD)+E*0.000533*Cos(2*ME-MS)                    &
     &+E*0.000401*Cos(2*ME-MD-MS)                                       &
     &+E*0.00032*Cos(MD-MS)-0.000271*Cos(ME)                            &
     &-E*0.000264*Cos(MD+MS)-0.000198*Cos(2*MF-MD)                      &
     &+0.000173*Cos(3*MD)+0.000167*Cos(4*ME-MD)                         &
     &-E*0.000111*Cos(MS)+0.000103*Cos(4*ME-2*MD)                       &
     &-0.000084*Cos(2*MD-2*ME)-E*0.000083*Cos(2*ME+MS)                  &
     &+0.000079*Cos(2*ME+2*MD)+0.000072*Cos(4*ME)                       &
     &+E*0.000064*Cos(2*ME-MS+MD)-E*0.000063*Cos(2*ME+MS-MD)            &
     &+E*0.000041*Cos(MS+ME)+E*0.000035*Cos(2*MD-MS)                    &
     &-0.000033*Cos(3*MD-2*ME)-0.00003*Cos(MD+ME)                       &
     &-0.000029*Cos(2*(MF-ME))-E*0.000029*Cos(2*MD+MS)                  &
     &+E2*0.000026*Cos(2*(ME-MS))-0.000023*Cos(2*(MF-ME)+MD)            &
     &+E*0.000019*Cos(4*ME-MD-MS)
      PM=PM*pic

!  geocentric distance MO3 in km
      MO3=6378.14/Sin(PM)

!  geocentric equatorial coordinates of the Moon
      SW=-1
      Call eqecl(pic,pi2,MO1,MO2,X1,X2,ecl,SW)
      res(1,2)=X1
      res(2,2)=X2
      res(3,2)=MO3
! 
      End Subroutine moon

      Subroutine anomaly(pi2,AM,EC,AT,AE)
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
! author        : Ann-Margrit Hellmich
! 
! last modified : 03/18/2003
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
!  calculation of true anomaly AT and eccentric anomaly AE
!  given mean anomaly AM and eccentricity EC
!  according to Duffett, 1990
! 
      Real pi2,AM,EC,AT,AE,M,D,A
! 
!  see page 113
      M=AM-pi2*Int(AM/pi2)
      AE=M
  1   D=AE-(EC*Sin(AE))-M
      If(Abs(D).Ge.0.000006) Then
        D=D/(1.0-EC*Cos(AE))
        AE=AE-D
        Goto 1
      Else
        A=Sqrt((1.0+EC)/(1.0-EC))*Tan(AE/2.0)
        AT=2.0*Atan(A)
      Endif
! 
      End Subroutine anomaly

      Subroutine eqecl(pic,pi2,X,Y,P,Q,ecl,SW)
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
! author        : Ann-Margrit Hellmich
! 
! last modified : 03/30/2003
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
!  conversion of ecliptic into equatorial coordinates
!  according to Duffett, 1990
! 
! if SW=+1: equatorial (X,Y..alpha,delta) to ecliptic (P,Q..lambda,beta)
! if SW=-1: equatorial (X,Y..lambda,beta) to ecliptic (P,Q..alpha,delta)
! 
      Real pic,pi2,ecl,P,Q,X,Y,SW
! 
!  see page 62
      P=Atan2((Sin(X)*Cos(ecl)+Tan(Y)*Sin(ecl)*SW),Cos(X))
      If(P.Lt.0) Call negangle2(pi2,P)
      If(P.Ge.pi2) Call langle2(pi2,P)
      Q=Asin(Sin(Y)*Cos(ecl)-Cos(Y)*Sin(ecl)*Sin(X)*SW)     
! 
      End Subroutine eqecl

      Subroutine aufb2(sidt,res,res2)
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
! author        : Ann-Margrit Hellmich
! 
! last modified : 01/07/2003
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
! modifications according to "ephaufb.f" by Maik Thomas
! (rekt(rad)->sid.time.green.-r.asc.; dekl(rad); cri3->(a/r)^3
! for Sun and Moon)
! 
      Real sidt,h(3)
      Real res(3,2),res2(3,2)
! 
      res2(1,1)=sidt-res(1,1)
      res2(1,2)=sidt-res(1,2)
      res2(2,1)=res(2,1)
      res2(2,2)=res(2,2)
      h(1)=1./res(3,1)
      h(2)=384400./res(3,2)
      res2(3,1)=h(1)*h(1)*h(1)
      res2(3,2)=h(2)*h(2)*h(2)
! 
      Return
      End Subroutine aufb2

      Subroutine negangle2(pi2,x)
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
! author        : Ann-Margrit Hellmich
! 
! last modified : 08/22/2002
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
!  transformation of negative angles
!  to angles in interval [0°;360°)
! 
      Logical endvar
      Real pi2,x
! 
      endvar=.False.
    1 If(.Not.endvar) Then
        x=x+pi2
        endvar=(x.Ge.0)
        Goto 1
      Endif
! 
      Return
      End  Subroutine  negangle2 

      Subroutine langle2(pi2,x)
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
! author        : Ann-Margrit Hellmich
! 
! last modified : 08/22/2002
! 
!! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
!  transformation of large angles
!  to angles in interval [0°;360°)
! 
      Logical endvar
      Real pi2,x
! 
      endvar=.False.
    1 If(.Not.endvar) Then
        x=x-pi2
        endvar=(x.Lt.pi2)
        Goto 1
      Endif
! 
      Return
      End Subroutine    langle2

SUBROUTINE eph(jul,mon,jahrph,moph)
  !     Phase bezueglich 1.1.2000 0 h
  !     NB 2000 war Schaltjahr
  INTEGER ::  mole(12),jul,mon,jahrph,moph,jpl,jj,m
  DATA mole/31,28,31,30,31,30,31,31,30,31,30,31/

  mole(2)=28
#ifdef RYEAR
  IF (   (MOD(jUL,4).EQ.0                                       &
       .AND. MOD(jul,100).NE.0)                                 &
       .OR.  MOD(jul,400).EQ.0  )       mole(2)=29
#endif
  IF(jul.LT.2000) THEN
     jahrph=0
     DO jj=jul,1999
        jpl=365
#ifdef RYEAR
        IF (   (MOD(jj,4).EQ.0                                   &
             .AND. MOD(jj,100).NE.0)                                  &
             .OR.  MOD(jj,400).EQ.0  ) jpl = 366
#endif
        jahrph=jahrph-jpl
     ENDDO
  ENDIF

  IF(jul.EQ.2000) jahrph=0

  IF(jul.GT.2000) THEN
     jahrph=0
     DO jj=2000,jul-1
        jpl=365
#ifdef RYEAR
        IF (   (MOD(jj,4).EQ.0                                   &
             .AND. MOD(jj,100).NE.0)                                  &
             .OR.  MOD(jj,400).EQ.0  ) jpl = 366
#endif
        jahrph=jahrph+jpl
     ENDDO
  ENDIF

  moph=0

  IF (mon.GT.1)THEN
     DO m=1,mon-1
        moph=moph+mole(m)
     ENDDO
  ENDIF
  !     jahrph=phase zu beginn des Jahres
  !     gesamtphase=jahrph+moph

END SUBROUTINE eph

      END MODULE mo_tidal


