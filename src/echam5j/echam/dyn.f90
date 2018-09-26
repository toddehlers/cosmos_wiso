SUBROUTINE dyn(petadot,ulz,vmaxz)

  ! Description:
  !
  ! Computes adiabatic tendencies and auxilliary hybrid variables.
  !
  ! Method:
  !
  ! The primary purpose is to compute adiabatic tendencies.
  ! Auxilliary variables connected with the vertical difference
  ! scheme are saved as they are subsequently required to compute
  ! input to physical parameterizations.
  !
  ! *dyn* is called from *gpc*. 
  ! Input is from long-term storage and modules *mo_hyb* and *mo_gaussgrid*
  ! Output is to long-term storage.
  !
  ! Externals:
  ! *pres* and *auxhyb* are called to calculate auxilliary variables.
  ! *geopot* is called to calculate full-level geopotentials.
  ! *locate*, *alloc* and *unloc* are called to manage storage.
  !
  ! *External documentation of the model equations, of
  ! the organization of the vertical calculation, and of the
  ! organization of the spectral horizontal calculation.
  !
  ! Authors:
  !
  ! A. J. Simmons, ECMWF, January 1982, original source
  ! U. Schlese, DKRZ, June 1991, changes for new advection scheme
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, August 1998, tendency diagnostics
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
  USE mo_scan_buffer,   ONLY: alnpr, alpha, alps, alpste, d, qte, rh,  &
                              t, tte, u, ul, v, vo, vervel, vol, vom,  &
                              xite, xlte, dalpsl, dalpsm,              &
                              dtl, dtm, xtte
  USE mo_memory_gl,     ONLY: q, xl, xi
!---wiso-code
  USE mo_memory_wiso,   ONLY: wisoq,   wisoxl,   wisoxi,               & 
                              wisoqte, wisoxlte, wisoxite 
  USE mo_wiso,          ONLY: lwiso
!---wiso-code-end
  USE mo_memory_g3a,    ONLY: apsm, geospm
  USE mo_memory_g3b,    ONLY: aps
  USE mo_control,       ONLY: ltdiag, nlev, nlevp1, nvclev, vct
  USE mo_gaussgrid,     ONLY: gl_rsqcst
  USE mo_geoloc,        ONLY: coriol_2d, cst_2d, rcsth_2d
  USE mo_hyb,           ONLY: alpham, alrrdic, ardprc, cetah, cpg,     &
                              delb, delpr, nlevm1, nlmsgl, nlmslp,     &
                              nplev, nplvp1, ralpha, rddelb, rdlnp0i,  &
                              rdt0ral, rlnpr, t0icao
  USE mo_constants,     ONLY: cpd, rcpd, rd, vtmpc1, vtmpc2
  USE mo_diag_tendency, ONLY: pdiga, pdsga
  USE mo_global_op,     ONLY: maxval_zonal, minval_zonal
  USE mo_decomposition, ONLY: ldc => local_decomposition
  USE mo_memory_gl,     ONLY: xt
  USE mo_transpose,     ONLY: reorder
  USE mo_column,        ONLY: get_col_pt, &! prescribe ver.vel. in SCM run
                              lcotra,     &! column model or 'traject' run
                              int_vtr      ! flag: vertical transport in SCM
  USE mo_advection

  IMPLICIT NONE

  !  Array arguments (now indexed by continuous latitude index)

  REAL(dp) ,INTENT(OUT) :: petadot(ldc% nproma ,nlevp1, ldc% ngpblks)
  REAL(dp) ,INTENT(OUT) :: ulz    (             nlev  , ldc% nglat)
  REAL(dp) ,INTENT(OUT) :: vmaxz  (             nlev  , ldc% nglat)

  !  Local array bounds
  INTEGER :: nglat, nglon, nproma, ngpblks, nbdim

  !  Local scalars: 
  REAL(dp) :: dpde, sdoth, zalpha, zalphm, zardrc, zcorio, zcpdrrd,    &
              zcpg, zcst, zdelb, zdlp, zdt, zetdpde, zln, zlnpr,       &
              zrcsth, zrddlb, zrdwrp, zrgr, zrrd, zrtgr, ztq, zvcikp
  INTEGER :: ikp, jrow, jk, jl, jglat

  !  Local arrays: 
  REAL(dp) :: zgeos  (ldc%nproma)         !
  REAL(dp) :: ztbar  (ldc%nproma, nlev)   !
  REAL(dp) :: ztv    (ldc%nproma, nlev)   !
  REAL(dp) :: aph    (ldc%nproma ,nlevp1) ! pressure
  REAL(dp) :: zdpsl  (ldc%nproma        ) ! surface pressure gradient
  REAL(dp) :: zdpsm  (ldc%nproma        ) ! surface pressure gradient
  REAL(dp) :: delp   (ldc%nproma ,nlev  ) ! pressure difference across layers
  REAL(dp) :: zrdelp (ldc%nproma ,nlev  ) ! reciprocal of *pdelp.*
  REAL(dp) :: zvgrps (ldc%nproma ,nlev  ) ! v * grad(ps)
  REAL(dp) :: zsdiv  (ldc%nproma ,nlevp1) ! surface pressure tendency
  REAL(dp) :: zru    (ldc%nglon  ,nlev, ldc% nglat)
  REAL(dp) :: zrv    (ldc%nglon  ,nlev, ldc% nglat)
  REAL(dp) :: zrcst  (ldc%nglat)

  !  External subroutines 
  EXTERNAL auxhyb, geopot, pres

  !  Intrinsic functions 
  INTRINSIC EXP, SQRT


  !  Executable statements 

  !-- 1. Preliminary calculations

  !-- 1.1 Set local values

  zcpdrrd = cpd/rd
  zrrd    = 1.0_dp/rd

  !  Local array bounds

  nglon   = ldc% nglon   ! number of longitudes
  nglat   = ldc% nglat   ! number of latitudes
  ngpblks = ldc% ngpblks ! number of rows
  nbdim   = ldc% nproma

!CSD$ PARALLEL DO PRIVATE(nproma, dpde, sdoth, zalpha, zalphm, zardrc,     &
!CSD$&  zcorio, zcpg, zcst, zdelb, zdlp, zdt, zetdpde, zln, zlnpr, zrcsth, &
!CSD$&  zrddlb, zrdwrp, zrgr, zrtgr, ztq, zvcikp, ikp, jrow, jk, jl,       &
!CSD$&  zgeos, ztbar, ztv, aph, zdpsl, zdpsm, delp, zrdelp, zvgrps, zsdiv, &
!CSD$&  zru, zrv, zrcst)
!$OMP PARALLEL PRIVATE(nproma, dpde, sdoth, zalpha, zalphm, zardrc,       &
!$OMP  zcorio, zcpg, zcst, zdelb, zdlp, zdt, zetdpde, zln, zlnpr, zrcsth, &
!$OMP  zrddlb, zrdwrp, zrgr, zrtgr, ztq, zvcikp, ikp, jrow, jk, jl,       &
!$OMP  zgeos, ztbar, ztv, aph, zdpsl, zdpsm, delp, zrdelp, zvgrps, zsdiv, &
!$OMP  zru, zrv, zrcst)
!$OMP DO
  DO jrow = 1, ngpblks

    IF ( jrow == ngpblks ) THEN
      nproma = ldc% npromz
    ELSE
      nproma = ldc% nproma
    END IF

!-- 1.2 Compute surface pressure and its gradient

    aph(1:nproma,nlevp1) = EXP(alps(1:nproma,jrow))
    zdpsl(1:nproma) = aph(1:nproma,nlevp1)*dalpsl(1:nproma,jrow)
    zdpsm(1:nproma) = aph(1:nproma,nlevp1)*dalpsm(1:nproma,jrow)

!-- 1.3 Compute half-level pressures and auxilliary variables.

    CALL pres(aph,nbdim,aph(1,nlevp1),nproma)
    CALL auxhyb(delp,zrdelp,alnpr(:,:,jrow),alpha(:,:,jrow),   &
              aph,nbdim,nproma)

!-- 1.4 Compute v.grad(ps)

    DO jk = nplvp1, nlev
      DO jl = 1, nproma
        zvgrps(jl,jk) = u(jl,jk,jrow)*zdpsl(jl)                    &
                                        +(v(jl,jk,jrow)*zdpsm(jl))
      END DO
    END DO

!-- 2. Sum divergence and compute surface pressure tendency

!-- 2.1 Compute pressure-level sums

    zsdiv (1:nproma,1) = 0._dp

    DO jk = 1, nplev
      ikp = jk + 1
      zdlp = delpr(jk)

      DO jl = 1, nproma
        zsdiv(jl,ikp) = d(jl,jk,jrow)*zdlp + zsdiv(jl,jk)
      END DO
    END DO

!-- 2.2 Compute hybrid-level sums

    DO jk = nplvp1, nlev
      ikp = jk + 1
      zdelb = delb(jk)

      DO jl = 1, nproma
        zsdiv(jl,ikp) = d(jl,jk,jrow)*delp(jl,jk) +                &
                                  zdelb*zvgrps(jl,jk) + zsdiv(jl,jk)
      END DO
    END DO

!-- prescribe vertical velocity in SCM run

    IF (lcotra) CALL get_col_pt (zsdiv, d, zvgrps, delp, aph, jrow)

    IF (ltdiag) THEN
       ! store pressure tendency term for each layer
       DO jk = 1,nlev
          pdiga(1:nglon,jk,20,jrow) = (zsdiv(:,jk+1) -                 &
                                            zsdiv(:,jk))/aph(:,nlevp1)
       ENDDO
    END IF

!-- 2.3 Tendency of logarithm of surface pressure

    DO jl = 1, nproma
      alpste(jl,jrow) = -zsdiv(jl,nlevp1)/aph(jl,nlevp1)
    END DO

    IF (ltdiag) pdsga(1:nglon,1,jrow) = pdsga(1:nglon,1,jrow) +        &
                                                    alpste(:,jrow)

!-- 3. Compute reference temperature and deviation
!      of virtual temprature

    DO jl = 1, nproma
      zgeos(jl) = rd*alps(jl,jrow)
    END DO

    DO jk = nlev, 1, -1

      DO jl = 1, nproma
        ztbar(jl,jk) = t0icao*EXP(alrrdic*(zgeos(jl)-                  &
                                       alpha(jl,jk,jrow)-rdlnp0i))
        zgeos(jl) = zgeos(jl) - alnpr(jl,jk,jrow)
        ztv(jl,jk) = t(jl,jk,jrow)*                                &
          (1._dp+vtmpc1*q(jl,jk,jrow)-(xl(jl,jk,jrow)+xi(jl,jk,jrow))) &
           - ztbar(jl,jk)
      END DO
    END DO

!-- 4. Compute vertical advection
!      and do zonal mean and box diagnostics.

!-- 4.1 Compute vertical advection

    vom    (1:nproma,1,jrow) = 0._dp
    vol    (1:nproma,1,jrow) = 0._dp
    tte    (1:nproma,1,jrow) = 0._dp
    qte    (1:nproma,1,jrow) = 0._dp
    xlte   (1:nproma,1,jrow) = 0._dp
    xite   (1:nproma,1,jrow) = 0._dp
    vervel (1:nproma,1,jrow) = 0._dp

!---wiso-code
    IF (lwiso) THEN
      wisoqte  (1:nproma,1,:,jrow) = 0._dp
      wisoxlte (1:nproma,1,:,jrow) = 0._dp
      wisoxite (1:nproma,1,:,jrow) = 0._dp
    ENDIF
!---wiso-code-end

    DO jk = 1, nlevm1
      ikp = jk + 1
      zvcikp = vct(nvclev+ikp)

      DO jl = 1, nproma
        sdoth = 0.5_dp*(zvcikp*zsdiv(jl,nlevp1)-zsdiv(jl,ikp))
        vom(jl,jk ,jrow) = (u(jl,jk,jrow)-u(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,jk)) + vom(jl,jk,jrow)
        vom(jl,ikp,jrow) = (u(jl,jk,jrow)-u(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,ikp))
        vol(jl,jk ,jrow) = (v(jl,jk,jrow)-v(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,jk)) + vol(jl,jk,jrow)
        vol(jl,ikp,jrow) = (v(jl,jk,jrow)-v(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,ikp))
        tte (jl,jk,jrow) = (t(jl,jk,jrow)-t(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,jk)) + tte(jl,jk,jrow)
        tte (jl,ikp,jrow) =(t(jl,jk,jrow)-t(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,ikp))

        ! compute eta-dot for semi Lagrangian scheme

        zetdpde = 2._dp*sdoth
        dpde = (aph(jl,jk+2)-aph(jl,jk))/(cetah(jk+2)-cetah(jk))
        IF (iadvec == semi_lagrangian .AND. .NOT. ldc% col_1d)         &
          petadot(jl,jk+1,jrow) = zetdpde/dpde
        ! compute eta-dot for Spitfire scheme
        IF (iadvec == spitfire .AND. .NOT. ldc% col_1d)                &
          petadot(jl,jk+1,jrow) = zetdpde/dpde   
      END DO
    END DO

    IF (ldc% col_1d .AND. int_vtr==1) THEN
      DO jk = 1, nlevm1
        ikp = jk + 1
        zvcikp = vct(nvclev+ikp)

        DO jl = 1, nproma
          sdoth = 0.5_dp*(zvcikp*zsdiv(jl,nlevp1)-zsdiv(jl,ikp))

          ! vertical transport in single column model
          ! this is a poor man's replacement for the spitfire transport

            qte(jl,jk ,jrow)     = (q(jl,jk,jrow)-q(jl,ikp,jrow))* &
              (sdoth*zrdelp(jl,jk )) + qte(jl,jk,jrow)
            qte(jl,ikp,jrow)     = (q(jl,jk,jrow)-q(jl,ikp,jrow))* &
              (sdoth*zrdelp(jl,ikp))
            xlte(jl,jk ,jrow)  = (xl(jl,jk,jrow)-xl(jl,ikp,jrow))* &
              (sdoth*zrdelp(jl,jk )) + xlte(jl,jk,jrow)
            xlte(jl,ikp,jrow)  = (xl(jl,jk,jrow)-xl(jl,ikp,jrow))* &
              (sdoth*zrdelp(jl,ikp))
            xite(jl,jk ,jrow)  = (xi(jl,jk,jrow)-xi(jl,ikp,jrow))* &
              (sdoth*zrdelp(jl,jk )) + xite(jl,jk,jrow)
            xite(jl,ikp,jrow)  = (xi(jl,jk,jrow)-xi(jl,ikp,jrow))* &
              (sdoth*zrdelp(jl,ikp))
            xtte(jl,jk ,:,jrow)=                                   &
                                 (xt(jl,jk,:,jrow)-xt(jl,ikp,:,jrow))* &
                       (sdoth*zrdelp(jl,jk )) + xtte(jl,jk,:,jrow)
            xtte(jl,ikp,:,jrow)=                                   &
                                 (xt(jl,jk,:,jrow)-xt(jl,ikp,:,jrow))* &
                        (sdoth*zrdelp(jl,ikp))
        END DO
!---wiso-code
        IF (lwiso) THEN
         DO jl = 1, nproma
          sdoth = 0.5_dp*(zvcikp*zsdiv(jl,nlevp1)-zsdiv(jl,ikp))
            wisoqte(jl,jk ,:,jrow)  = (wisoq(jl,jk,:,jrow)-wisoq(jl,ikp,:,jrow))* &
              (sdoth*zrdelp(jl,jk )) + wisoqte(jl,jk,:,jrow)
            wisoqte(jl,ikp,:,jrow)  = (wisoq(jl,jk,:,jrow)-wisoq(jl,ikp,:,jrow))* &
              (sdoth*zrdelp(jl,ikp))
            wisoxlte(jl,jk ,:,jrow) = (wisoxl(jl,jk,:,jrow)-wisoxl(jl,ikp,:,jrow))* &
              (sdoth*zrdelp(jl,jk )) + wisoxlte(jl,jk,:,jrow)
            wisoxlte(jl,ikp,:,jrow) = (wisoxl(jl,jk,:,jrow)-wisoxl(jl,ikp,:,jrow))* &
              (sdoth*zrdelp(jl,ikp))
            wisoxite(jl,jk ,:,jrow) = (wisoxi(jl,jk,:,jrow)-wisoxi(jl,ikp,:,jrow))* &
              (sdoth*zrdelp(jl,jk )) + wisoxite(jl,jk,:,jrow)
            wisoxite(jl,ikp,:,jrow) = (wisoxi(jl,jk,:,jrow)-wisoxi(jl,ikp,:,jrow))* &
              (sdoth*zrdelp(jl,ikp))
         END DO
        END IF
!---wiso-code-end
      END DO
    ENDIF
 
    IF (ltdiag) THEN
       ! store vertical advection increment
       pdiga(1:nglon,:, 2,jrow) = pdiga(1:nglon,:, 2,jrow) + vom(:,:,jrow)
       pdiga(1:nglon,:, 7,jrow) = pdiga(1:nglon,:, 7,jrow) + vol(:,:,jrow)
       pdiga(1:nglon,:,13,jrow) = pdiga(1:nglon,:,13,jrow) + tte(:,:,jrow)
       ! prepare next parts
       pdiga(1:nglon,:, 1,jrow) = pdiga(1:nglon,:, 1,jrow) - vom(:,:,jrow)
       pdiga(1:nglon,:, 6,jrow) = pdiga(1:nglon,:, 6,jrow) - vol(:,:,jrow)
       pdiga(1:nglon,:,14,jrow) = pdiga(1:nglon,:,14,jrow) - tte(:,:,jrow)
    END IF

    DO jl = 1, nproma
      petadot(jl,     1,jrow) = 0._dp
      petadot(jl,nlevp1,jrow) = 0._dp
    END DO

!-- 5. Compute energy conversion term for pressure levels
!      and do zonal mean and box diagnostics.

!-- 5.1 Compute energy conversion term for pressure levels

    DO jk = 1, nplev
      zardrc = ardprc(jk)
      zalphm = alpham(jk)

!OCL NOALIAS
      DO jl = 1, nproma
        ztq = (ztbar(jl,jk)+ztv(jl,jk))/(1._dp+vtmpc2*q(jl,jk,jrow))
        zdt = -ztq*(zsdiv(jl,jk)*zardrc+zalphm*d(jl,jk,jrow))
        tte(jl,jk,jrow) = tte(jl,jk,jrow) + zdt
      END DO
    END DO

!-- 6. Compute pressure-gradient terms,complete calculation
!      of energy conversion term

!-- 6.1 Hybrid levels

    DO jk = nplvp1, nlmsgl
      zdelb = delb(jk)
      zrddlb = rddelb(jk)
      zcpg = cpg(jk)

!OCL NOVREC,NOALIAS
      DO jl = 1, nproma
        zrgr = (zrddlb+zcpg*alnpr(jl,jk,jrow)*zrdelp(jl,jk))*      &
                                                          zrdelp(jl,jk)
        zrtgr = zrgr*ztv(jl,jk)
        zcst = cst_2d(jl,jrow)
        vom(jl,jk,jrow) = vom(jl,jk,jrow) - zcst*zrtgr*zdpsl(jl)
        vol(jl,jk,jrow) = vol(jl,jk,jrow) - zcst*zrtgr*zdpsm(jl)
        zrdwrp = (zrgr*zvgrps(jl,jk)-(zrdelp(jl,jk)*(zsdiv(jl,jk)*     &
                 alnpr(jl,jk,jrow)+alpha(jl,jk,jrow)*zdelb*    &
                 zvgrps(jl,jk))+                                       &
                 alpha(jl,jk,jrow)*d(jl,jk,jrow)))*rcpd
        ztq = (ztbar(jl,jk)+ztv(jl,jk))/(1._dp+vtmpc2*q(jl,jk,jrow))
        zdt = ztq*zrdwrp
        tte(jl,jk,jrow) = tte(jl,jk,jrow) + zdt

      END DO
    END DO

!-- 6.2 Sigma levels

    DO jk = nlmslp, nlev
      zalphm = alpham(jk)
      zlnpr = rlnpr(jk)
      zalpha = ralpha(jk)

      zrddlb = rddelb(jk)
!OCL NOALIAS
      DO jl = 1, nproma
        zrgr = zrddlb*zrdelp(jl,jk)
        zrtgr = zrgr*ztv(jl,jk)
        zcst = cst_2d(jl,jrow)
        vom(jl,jk,jrow) = vom(jl,jk,jrow) - zcst*zrtgr*zdpsl(jl)
        vol(jl,jk,jrow) = vol(jl,jk,jrow) - zcst*zrtgr*zdpsm(jl)
        zrdwrp = (zrgr*zvgrps(jl,jk)-(zrdelp(jl,jk)*(zsdiv(jl,jk)*     &
                 zlnpr+zalphm*                                         &
                 zvgrps(jl,jk))+zalpha*d(jl,jk,jrow)))*rcpd
        ztq = (ztbar(jl,jk)+ztv(jl,jk))/(1._dp+vtmpc2*q(jl,jk,jrow))
        zdt = ztq*zrdwrp
        tte(jl,jk,jrow) = tte(jl,jk,jrow) + zdt

      END DO
    END DO
 
    IF (ltdiag) THEN
       ! store energy conversion increment
       pdiga(1:nglon,:,14,jrow) = pdiga(1:nglon,:,14,jrow) + tte (:,:,jrow)
       ! prepare next part
       pdiga(1:nglon,:,12,jrow) = pdiga(1:nglon,:,12,jrow) - tte (:,:,jrow)
    END IF

!-- 6.3 Compute vertical velocity for mass-flux scheme

    DO jk = 1, nplev
      zardrc = ardprc(jk)
      zalphm = alpham(jk)
      DO jl = 1, nproma
        vervel(jl,jk,jrow) = &
          -(zsdiv(jl,jk)*zardrc+zalphm*d(jl,jk,jrow))*zcpdrrd
      END DO
    END DO

    DO jk = nplvp1, nlmsgl
      zdelb = delb(jk)
      zrddlb = rddelb(jk)
      zcpg = cpg(jk)
      DO jl = 1, nproma
        zrgr = (zrddlb+zcpg*alnpr(jl,jk,jrow)*zrdelp(jl,jk))*      &
                                                         zrdelp(jl,jk)
        vervel(jl,jk,jrow) = (zrgr*zvgrps(jl,jk)-zrdelp(jl,jk)*    &
                (zsdiv(jl,jk)*alnpr(jl,jk,jrow)+                   &
                alpha(jl,jk,jrow)*zdelb*zvgrps(jl,jk))-            &
                alpha(jl,jk,jrow)*d(jl,jk,jrow))*zrrd
      END DO
    END DO

    DO jk = nlmslp, nlev
      zalphm = alpham(jk)
      zlnpr = rlnpr(jk)
      zalpha = ralpha(jk)
      zrddlb = rddelb(jk)
      DO jl = 1, nproma
        zrgr = zrddlb*zrdelp(jl,jk)
        vervel(jl,jk,jrow) = (zrgr*zvgrps(jl,jk)-zrdelp(jl,jk)*    &
                (zsdiv(jl,jk)*zlnpr+zalphm*zvgrps(jl,jk))-             &
                 zalpha*d(jl,jk,jrow))*zrrd
      END DO
    END DO

    DO jk = 1, nlev
      DO jl = 1, nproma
        vervel(jl,jk,jrow) = vervel(jl,jk,jrow)                &
                                    * 0.5_dp*(aph(jl,jk)+aph(jl,jk+1))
      END DO
    END DO

!-- 7. Compute geopotential

!-- 7.1 Compute deviation of geopotential height at surface

    DO jl = 1, nproma
      zln = rd*alps(jl,jrow) - rdlnp0i
      zgeos(jl) = geospm(jl,jrow) + rdt0ral*EXP(alrrdic*zln)
    END DO

!-- 7.2 Compute deviation of geopotential height

    CALL geopot(rh(:,:,jrow),ztv,alnpr(:,:,jrow),              &
                 alpha(:,:,jrow), zgeos, nbdim, nproma)

!-- 8. Compute horizontal advection terms

    DO jk = 1, nlev
!OCL NOVREC,NOALIAS
      DO jl = 1, nproma
        zcorio = coriol_2d(jl,jrow)
        zrcsth = rcsth_2d(jl,jrow)
        rh(jl,jk,jrow)=zrcsth*(u(jl,jk,jrow)*u(jl,jk,jrow)+&
          v(jl,jk,jrow)*v(jl,jk,jrow)) + rh(jl,jk,jrow)
        vom(jl,jk,jrow)=(vo(jl,jk,jrow)+zcorio)*               &
          v(jl,jk,jrow)+vom(jl,jk,jrow)
        vol(jl,jk,jrow) = -(vo(jl,jk,jrow)+zcorio)*            &
          u(jl,jk,jrow)+vol(jl,jk,jrow)
        zdt = -u(jl,jk,jrow)*dtl(jl,jk,jrow)                   &
              - v(jl,jk,jrow)*dtm(jl,jk,jrow)
        tte (jl,jk,jrow) = tte(jl,jk,jrow) + zdt

      END DO
    END DO
 
    IF (ltdiag) THEN
       ! pressure gradient term, horizontal advection and coriolisterm
       pdiga(1:nglon,:, 1,jrow) = pdiga(1:nglon,:, 1,jrow) + vom(:,:,jrow)
       pdiga(1:nglon,:, 6,jrow) = pdiga(1:nglon,:, 6,jrow) + vol(:,:,jrow)
       ! horizontal advection
       pdiga(1:nglon,:,12,jrow) = pdiga(1:nglon,:,12,jrow) + tte(:,:,jrow)
       ! G-term, potential energy and kinetic energy
       pdiga(1:nglon,:,11,jrow) = pdiga(1:nglon,:,11,jrow) - rh (:,:,jrow)
    END IF

!-- 10. Duplicate ps

    aps(1:nproma,jrow) = aph(1:nproma,nlevp1)
    apsm(1:nproma,jrow) = aph(1:nproma,nlevp1)

  END DO
!$OMP END DO
!$OMP END PARALLEL
!CSD$ END PARALLEL DO

!-- Compute maximum !u!+!v!

  DO jrow = 1, nglat
    jglat = ldc% glat(jrow)          ! global continuous north -> south
    zrcst(jrow) = gl_rsqcst(jglat)
  END DO

  CALL reorder(zru, u)
  CALL reorder(zrv, v)

  ul(:,:) = (maxval_zonal(zru(:,:,:))  + &
             minval_zonal(zru(:,:,:))) * 0.5_dp

  ulz(:,:) = ul(:,:)

  vmaxz(:,:) = SQRT(maxval_zonal(zru(:,:,:)*zru(:,:,:)                 &
                               + zrv(:,:,:)*zrv(:,:,:)))               &
                 * SPREAD(zrcst,dim=1,ncopies=nlev)
 
END SUBROUTINE dyn
