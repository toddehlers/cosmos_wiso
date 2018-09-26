SUBROUTINE radheat (kproma, kbdim, klev,  klevp1               &
                  , krow                                       &
                  , ldblrad                                    &
                  , pi0                                        &
                  , ptm1,         pqm1                         &
                  , ptrsof,       ptrsol                       &
                  , pemtef,       pemter                       &
                  , pemtef0,      ptrsof0                      &
                  , pemter1,      ptrsol1                      &
                  , pemtef01,     ptrsof01                     &
                  , pnetht_sw,    pnetht_lw                    &
                  , psrad0,       psrads                       &
                  , psradl,       psrafl                       &
                  , psrad0u,      psradsu                      &
                  , psraf0,       psrafs                       &
                  , psrad0d                                    &
                  , ptrad0,       ptrads                       &
                  , ptradl,       ptrafl                       &
                  , ptraf0,       ptrafs                       &
                  , ptradsu                                    &
                  , palbedo                                    &
                  , paphm1,       papm1                        &
                  , ptslnew,      ptte                         &
                  , pradtemp ,    pztsnew     )
!
!
!
!**** *RADHEAT* - COMPUTES TEMPERATURE CHANGES DUE TO RADIATION.
!
!
!     SUBJECT.
!     --------
!
!          THIS ROUTINE COMPUTES THE TENDENCIES OF THE ATMOSPHERE'S
!     TEMPERATURE DUE TO THE EFFECTS OF LONG WAVE AND SHORT WAVE
!     RADIATION. THE COMPUTATION IS DONE ON THE T-1 TIME LEVEL USING
!     VALUES OF ATMOSPHERIC TRANSMISIVITIES AND EMISSIVITIES THAT HAVE
!     BEEN STORED AT THE LAST FULL RADIATION TIME STEP. THE SURFACE
!     SOLAR FLUX LATER TO BE USED IN THE SOIL PROCESS CALCULATIONS IS
!     ALSO STORED.
!
!**   INTERFACE.
!     ----------
!
!          *RADHEAT* IS CALLED FROM *PHYSC*.
!
!     INPUT ARGUMENTS.
!     ----- ---------
!
!
!     OUTPUT ARGUMENTS.
!     ------ ---------
!
!
!     METHOD.
!     -------
!
!     PRODUCT OF SOLAR
!     INFLUX BY TRANSMISSIVITIES LEADS TO SOLAR FLUXES. THEN THE
!     TEMPERATURES ARE INTERPOLATED/EXTRAPOLATED TO THE LAYER BOUNDARIES
!     (AT THE BOTTOM ONE TAKES THE SURFACE TEMPERATURE) AND A PRODUCT BY
!     EMISSIVITIES OF SIGMA*T**4 GIVES THERMAL FLUXES. THE TWO FLUXES
!     ARE ADDED AND DIVERGENCES COMPUTED TO GIVE HEATING RATES.
!
!     EXTERNALS.
!     ----------
!
!          *SOLANG*.
!
!     AUTHOR.
!     ------
!
!     U. SCHLESE    DKRZ-HAMBURG    JUNE 1995

!     Modifications
!     U. Schlese, December 1999:  version for coupling
!     U. Schlese, July 2000, *solang* removed, weighted surface fluxes
!     I. Kirchner, May 2002, tendency diagnose bugfix for surface fluxes
!     S.J. Lorenz, Jan.2008, volcanic forcing (M.A.Thomas)
!
USE mo_kind,              ONLY: dp
USE mo_control,           ONLY: ltdiag
USE mo_constants,         ONLY: g, cpd, vtmpc2, stbo
USE mo_radiation,         ONLY: cemiss
USE mo_diag_tendency,     ONLY: pdiga
USE mo_time_control,      ONLY: delta_time
USE mo_accuflx_mem,       ONLY: d_aflx_sw, d_aflx_lw, d_aflx_swc,      &
                                d_aflx_lwc
!
IMPLICIT NONE
!
INTEGER :: kproma,kbdim,klev,klevp1,krow

! logicals for volcanic forcing and double radiation:
LOGICAL, INTENT(in) :: ldblrad

REAL(dp) ::                                                            &
        pi0(kbdim),                                                    &
        ptm1(kbdim,klev),      pqm1(kbdim,klev),                       &
        ptrsol(kbdim,klevp1),  ptrsof(kbdim,2),                        &
        pemter(kbdim,klevp1),  pemtef(kbdim,2),                        &
        ptrsof0(kbdim,klevp1), pemtef0(kbdim,klevp1),                  &
        pemter1(kbdim,klevp1), ptrsol1(kbdim,klevp1),                  &
        ptrsof01(kbdim,klevp1), pemtef01(kbdim,klevp1),                &
        pnetht_sw(kbdim,klev), pnetht_lw(kbdim,klev),                  &
        psrad0(kbdim),         psrads(kbdim),                          &
        psradl(kbdim),         psrafl(kbdim),                          &
        psrad0u(kbdim),        psradsu(kbdim),                         &
        psraf0(kbdim),         psrafs(kbdim),                          &
        psrad0d(kbdim),                                                &
        ptrad0(kbdim),         ptrads(kbdim),                          &
        ptradl(kbdim),         ptrafl(kbdim),                          &
        ptraf0(kbdim),         ptrafs(kbdim),                          &
        ptradsu(kbdim),                                                &
        palbedo(kbdim),                                                &
        paphm1(kbdim,klevp1),  papm1(kbdim,klev),                      &
        ptslnew(kbdim),        ptte(kbdim,klev),                       &
        ztrdown(kbdim),        pradtemp(kbdim),                        &
        ztrdown0(kbdim)
!
!
!    Local arrays
!
REAL(dp) :: zti(kbdim,klevp1),                                         &
            zflxs(kbdim,klevp1),  zflxt(kbdim,klevp1),                 &
            zflxs0(kbdim,klevp1), zflxt0(kbdim,klevp1),                &
            ztsnew(kbdim),        zteffl4(kbdim),                      &
            pztsnew(kbdim)
REAL(dp) :: ztrps(kbdim),  ztrpt(kbdim),  ztrpss(kbdim), ztrpts(kbdim)
REAL(dp) :: zdtime, zcons3, zdtdt, zfltop, zflbot,  zflts,             &
            zfltt, zffact, zsr0u, zsrsu, zdp1, zdp2, ztrsu

!    Extra local variables for volcanic calculations:
REAL(dp) :: zdtdt_sw,  &
            zdtdt_lw,  &
            zdtdt1_sw, zdtdt1_lw, ztrdown1, ztrdown01

REAL(dp) :: zflxs1(kbdim,klevp1), zflxt1(kbdim,klevp1),               &
            zflxs01(kbdim,klevp1), zflxt01(kbdim,klevp1),             &
            zflx_sw(kbdim,klevp1), zflx_lw(kbdim,klevp1),             &
            zflx_swc(kbdim,klevp1), zflx_lwc(kbdim,klevp1)

REAL(dp) :: zheat_sw(kbdim,klev), zheat_lw(kbdim,klev)

REAL(dp) :: zflxtoa_sw(kbdim), zflxtoa_lw(kbdim)

INTEGER :: jrow, jk, jl
!
!
! ----------------------------------------------------------------------
!
!*     1.   COMPUTATIONAL CONSTANTS.
!           ------------- ----------
!
100 CONTINUE
!
!
  zcons3=g/cpd
  zdtime = delta_time
!
!
  jrow = krow
!
!     ------------------------------------------------------------------
!
!*         3.     TEMPERATURES AT LAYERS' BOUDARIES.
!                 ------------ -- ------- ----------
!
300  CONTINUE
!
!*         3.1     INTERPOLATION PROPER.
!
310  CONTINUE
     DO 312 jk=2,klev
        DO 311 jl=1,kproma
           zti(jl,jk)=(ptm1(jl,jk-1)*papm1(jl,jk-1)                    &
                     *(papm1(jl,jk)-paphm1(jl,jk))                     &
                     +ptm1(jl,jk)*papm1(jl,jk)                         &
                     *(paphm1(jl,jk)-papm1(jl,jk-1)))                  &
                     /(paphm1(jl,jk)*(papm1(jl,jk)-papm1(jl,jk-1)))
311     END DO
312  END DO
!
!*        3.2     SURFACE AND TOP OF ATMOSPHERE TEMPERATURE.
!
320  CONTINUE
     DO 321 jl=1,kproma
        zti(jl,klevp1) = pradtemp(jl)
        zti(jl,1)=ptm1(jl,1)-papm1(jl,1)*(ptm1(jl,1)-zti(jl,2))     &
                  /(papm1(jl,1)-paphm1(jl,2))
321  END DO
!
!     ------------------------------------------------------------------
!
!*         4.    UPDATE FLUXES AND COMPUTE HEATING RATES.
!                ------ ------ --- ------- ------- -----
!
!    4.1 Fluxes at top of the atmosphere
!
     DO 401 jl=1,kproma
       zflxs(jl,1)=pi0(jl)*ptrsol(jl,1)
       zflxt(jl,1)=pemter(jl,1)
       zflxs0(jl,1)=pi0(jl)*ptrsof0(jl,1)
       zflxt0(jl,1)=pemtef0(jl,1)
401  END DO
!
!  for double radiation: - all/clear sky fluxes
   IF(ldblrad) THEN
     DO 4011 jl=1,kproma
       zflxs1(jl,1)=pi0(jl)*ptrsol1(jl,1)
       zflxt1(jl,1)=pemter1(jl,1)
       zflxs01(jl,1)=pi0(jl)*ptrsof01(jl,1)
       zflxt01(jl,1)=pemtef01(jl,1)
4011  END DO
   ENDIF
!
!
!     4.2  Fluxes and heating rates except for lowest layer
!
      DO 403 jk=1,klev-1
        DO 402 jl=1,kproma
          zfltop=zflxs(jl,jk)+zflxt(jl,jk)
          zflxs(jl,jk+1)=pi0(jl)*ptrsol(jl,jk+1)
          zflxt(jl,jk+1)=pemter(jl,jk+1)
          zflbot=zflxs(jl,jk+1)+zflxt(jl,jk+1)
          zdtdt=-zcons3*(zflbot-zfltop)/((paphm1(jl,jk+1)              &
                -paphm1(jl,jk))*(1._dp+vtmpc2*pqm1(jl,jk)))
          ptte(jl,jk)=ptte(jl,jk)+zdtdt
!
          zflxs0(jl,jk+1)=pi0(jl)*ptrsof0(jl,jk+1)
          zflxt0(jl,jk+1)=pemtef0(jl,jk+1)
402     END DO
403  END DO
!
!  for double radiation: 
   IF(ldblrad) THEN
      DO 4013 jk=1,klev-1
        DO 4012 jl=1,kproma
          zflxs1(jl,jk+1)=pi0(jl)*ptrsol1(jl,jk+1)
          zflxt1(jl,jk+1)=pemter1(jl,jk+1)
          zflxs01(jl,jk+1)=pi0(jl)*ptrsof01(jl,jk+1)
          zflxt01(jl,jk+1)=pemtef01(jl,jk+1)
4012     END DO
4013  END DO
   ENDIF
!
!
!     4.3  Lowest layer
!
     DO 404 jl=1,kproma      
       ztrdown(jl)=pemter(jl,klevp1)+cemiss*stbo*zti(jl,klevp1)**4
       zflxt(jl,klevp1)=ztrdown(jl)-cemiss*stbo*pztsnew(jl)**4
       ztrdown0(jl)=pemtef0(jl,klevp1)+cemiss*stbo*zti(jl,klevp1)**4
       zflxt0(jl,klevp1)=ztrdown0(jl)-cemiss*stbo*pztsnew(jl)**4
       zflxs(jl,klevp1)=pi0(jl)*ptrsol(jl,klevp1)
       zflxs0(jl,klevp1)=pi0(jl)*ptrsof0(jl,klevp1)
       zfltop=zflxs(jl,klev)+zflxt(jl,klev)
       zflbot=zflxs(jl,klevp1)+zflxt(jl,klevp1)
       zdtdt=-zcons3*(zflbot-zfltop)/((paphm1(jl,klevp1)               &
                 -paphm1(jl,klev))*(1._dp+vtmpc2*pqm1(jl,klev)))
       ptte(jl,klev)=ptte(jl,klev)+zdtdt
404  END DO
!
!
!  for double radiation: 
   IF(ldblrad) THEN
     DO 4014 jl=1,kproma      
!   all sky fluxes
       ztrdown1=pemter1(jl,klevp1)+cemiss*stbo*zti(jl,klevp1)**4
       zflxt1(jl,klevp1)=ztrdown1-cemiss*stbo*pztsnew(jl)**4
       zflxs1(jl,klevp1)=pi0(jl)*ptrsol1(jl,klevp1)
!   clear sky fluxes
       ztrdown01=pemtef01(jl,klevp1)+cemiss*stbo*zti(jl,klevp1)**4
       zflxt01(jl,klevp1)=ztrdown01-cemiss*stbo*pztsnew(jl)**4
       zflxs01(jl,klevp1)=pi0(jl)*ptrsof01(jl,klevp1)
4014  END DO
   ENDIF

IF (ldblrad) THEN
!  calculations for volcanic forcing:
    DO jk = 1, klev
      DO jl = 1, kproma

        zdtdt_sw=-zcons3*(zflxs(jl,jk+1)-zflxs(jl,jk))/((paphm1(jl,jk+1)    &
                 -paphm1(jl,jk))*(1.+vtmpc2*pqm1(jl,jk)))
        zdtdt_lw=-zcons3*(zflxt(jl,jk+1)-zflxt(jl,jk))/((paphm1(jl,jk+1)    &
                 -paphm1(jl,jk))*(1.+vtmpc2*pqm1(jl,jk)))

        zdtdt1_sw=-zcons3*(zflxs1(jl,jk+1)-zflxs1(jl,jk))/((paphm1(jl,jk+1) &
                 -paphm1(jl,jk))*(1.+vtmpc2*pqm1(jl,jk)))
        zdtdt1_lw=-zcons3*(zflxt1(jl,jk+1)-zflxt1(jl,jk))/((paphm1(jl,jk+1) &
                 -paphm1(jl,jk))*(1.+vtmpc2*pqm1(jl,jk)))

        zheat_sw(jl,jk) = zdtdt_sw - zdtdt1_sw
        zheat_lw(jl,jk) = zdtdt_lw - zdtdt1_lw

       pnetht_sw(jl,jk) = pnetht_sw(jl,jk)+zheat_sw(jl,jk)*zdtime
       pnetht_lw(jl,jk) = pnetht_lw(jl,jk)+zheat_lw(jl,jk)*zdtime


      ENDDO
     ENDDO

    DO jk = 1, klevp1
     DO jl = 1, kproma

         zflx_sw(jl,jk) = zflxs(jl,jk) - zflxs1(jl,jk)
         zflx_lw(jl,jk) = zflxt(jl,jk) - zflxt1(jl,jk)
         zflx_swc(jl,jk) = zflxs0(jl,jk) - zflxs01(jl,jk)
         zflx_lwc(jl,jk) = zflxt0(jl,jk) - zflxt01(jl,jk)


      d_aflx_sw(jl,jk,jrow) = d_aflx_sw(jl,jk,jrow) + zflx_sw(jl,jk)*zdtime
      d_aflx_lw(jl,jk,jrow) = d_aflx_lw(jl,jk,jrow) + zflx_lw(jl,jk)*zdtime
      d_aflx_swc(jl,jk,jrow) = d_aflx_swc(jl,jk,jrow) + zflx_swc(jl,jk)*zdtime
      d_aflx_lwc(jl,jk,jrow) = d_aflx_lwc(jl,jk,jrow) + zflx_lwc(jl,jk)*zdtime

     ENDDO
    ENDDO
ENDIF
!
     IF (ltdiag) THEN
       ! tendency diagnostics
       DO jk = 1, klev
         DO jl = 1,kproma
           zflts = zflxs(jl,jk+1)-zflxs(jl,jk)
           zfltt = zflxt(jl,jk+1)-zflxt(jl,jk)
           zffact = - zcons3/((paphm1(jl,jk+1)-paphm1(jl,jk)) &
                * (1._dp+vtmpc2*pqm1(jl,jk)) )
           pdiga(jl,jk,24,jrow) = pdiga(jl,jk,24,jrow) + zfltt*zffact
           pdiga(jl,jk,25,jrow) = pdiga(jl,jk,25,jrow) + zflts*zffact
         END DO
       END DO
     END IF
 
!     ------------------------------------------------------------------
!
!*         5.     Diagnostics of top and surface fluxes
!
!
!
     DO 510 jl = 1, kproma
       psrad0(jl) = psrad0(jl) + zdtime*zflxs(jl,1)
       ptrad0(jl) = ptrad0(jl) + zdtime*zflxt(jl,1)
       zsr0u = zflxs(jl,1) - pi0(jl)
       psrad0u(jl) = psrad0u(jl) + zdtime*zsr0u
       psrads(jl) = psrads(jl) + zdtime*zflxs(jl,klevp1)
       ptrads(jl) = ptrads(jl) + zdtime*zflxt(jl,klevp1)
       zsrsu = -zflxs(jl,klevp1)*(1._dp/(1._dp-palbedo(jl))-1._dp)
       psradsu(jl) = psradsu(jl) + zdtime*zsrsu
       psrad0d(jl) = psrad0d(jl) + zdtime*pi0(jl)
       psraf0(jl) = psraf0(jl) + zdtime*pi0(jl)*ptrsof(jl,1)
       psrafs(jl) = psrafs(jl) + zdtime*pi0(jl)*ptrsof(jl,2)
510  END DO
!
!
! Diagnostics of fluxes at 200mb
!
  DO 524 jk=1,klev
     DO 523 jl=1,kproma
        IF (paphm1(jl,jk) .LE. 20000._dp .AND.                         &
            paphm1(jl,jk+1) .GE. 20000._dp) THEN
            zdp1=paphm1(jl,jk)-paphm1(jl,jk+1)
            zdp2=paphm1(jl,jk)-20000._dp
            ztrps(jl)=zflxs(jl,jk)-(zflxs(jl,jk)-zflxs(jl,jk+1))       &
                        *(zdp2/zdp1)
            ztrpt(jl)=zflxt(jl,jk)-(zflxt(jl,jk)-zflxt(jl,jk+1))       &
                        *(zdp2/zdp1)
            ztrpss(jl)=zflxs0(jl,jk)-(zflxs0(jl,jk)-zflxs0(jl,jk+1))   &
                        *(zdp2/zdp1)
            ztrpts(jl)=zflxt0(jl,jk)-(zflxt0(jl,jk)-zflxt0(jl,jk+1))   &
                        *(zdp2/zdp1)
        END IF
523 END DO
524 END DO
  DO 528 jl=1,kproma
     psradl(jl)  =psradl(jl) + ztrps(jl)*zdtime
     ptradl(jl)  =ptradl(jl) + ztrpt(jl)*zdtime
     psrafl(jl)  =psrafl(jl) + ztrpss(jl)*zdtime
     ptrafl(jl)  =ptrafl(jl) + ztrpts(jl)*zdtime
528 END DO
!
!
     DO 520 jl=1,kproma
       ztrsu = -cemiss*stbo*pztsnew(jl)**4+(cemiss-1._dp)              &
               *(stbo*zti(jl,klevp1)**4+pemter(jl,klevp1)/cemiss)
       ptradsu(jl)=ptradsu(jl)+zdtime*ztrsu
       ptraf0(jl) =ptraf0(jl) +zdtime*pemtef(jl,1)
       ptrafs(jl) =ptrafs(jl) +zdtime*(pemtef(jl,2)+cemiss*stbo        &
                                  *(zti(jl,klevp1)**4-pztsnew(jl)**4))
520  END DO
!
  RETURN
END SUBROUTINE radheat
