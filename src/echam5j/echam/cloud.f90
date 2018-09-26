SUBROUTINE cloud (         kproma, kbdim, ktdia, klev, klevp1          &
!-----------------------------------------------------------------------
! - INPUT  2D .
                         , paphm1,   pvervel                           &
                         , papm1,    papp1,    pacdnc                  &
                         , pqm1,     ptm1,     ptvm1                   &
                         , pxlm1,    pxim1,    pxtec                   &
                         , pxvar,    pxskew,   pqtec                   &
                         , pbetaa,   pbetab                            &
                         , pvdiffp,  phmixtau, pvmixtau                &
                         , pgeo,     pbetass                           &
! - INPUT  1D .
                         , knvb                                        &
! - OUTPUT 2D .
                         , paclc,    paclcac,  prelhum                 &
! - INPUT/OUTPUT 1D .
                         , paclcov,  paprl,    pqvi                    &
                         , pxlvi,    pxivi                             &
! - OUTPUT 1D .
                         , pssfl,    prsfl                             &
! - INPUT/OUTPUT 2D .
                         , pqte,     ptte                              &
                         , pxlte,    pxite                             &
! - INPUT/OUTPUT 1D .
                         , paprs                                       &
!---wiso-code
                         , lwiso, kwiso                                &
                         , pwisoqm1                                    &
                         , pwisoxlm1, pwisoxim1                        &
                         , pwisoxtec                                   &
                         , pwisoqtec                                   &
                         , pwisoaprl, pwisoqvi                         &
                         , pwisoxlvi, pwisoxivi                        &
                         , pwisossfl, pwisorsfl                        &               
                         , pwisoqte                                    &
                         , pwisoxlte, pwisoxite                        &
                         , pwisoaprs)
!---wiso-code-end
!
!     *Cloud* computes large-scale water phase changes, precipitation,
!             cloud cover, and vertical integrals of specific humidity,
!             cloud liquid water content and cloud ice (diagnostics).
!
!     Subject.
!     --------
!
!          This routine computes the tendencies of the four prognostic
!          variables (temperature t, specific humidity q, cloud liquid
!          water xl, cloud ice xi) due to phase changes (condensation/
!          deposition, evaporation/sublimation of rain/snow falling
!          into the unsaturated part of the grid box, melting of snow,
!          melting/freezing of cloud ice/cloud water, sedimentation of
!          cloud ice, and precipitation formation in warm, cold and
!          mixed phase clouds.
!          The precipitation at the surface (rain and snow) is used in
!          later for computing the land surface hydrology in *surf*.
!          The cloud parameters (cloud cover, cloud liquid water and
!          cloud ice are used for the calculation of radiation at the
!          next timestep.
!          Attention: 
!          In the current version the advective tendencies of skewness 
!          and variance are set to zero.
!
!     INTERFACE.
!     ----------
!
!     *Call cloud*
!
!     Input arguments.
!     ----- ----------
!  - 2D
!  paphm1   : pressure at half levels                              (n-1)
!  papm1    : pressure at full levels                              (n-1)
!  papp1    : pressure at full levels                              (n-1)
!  pacdnc   : cloud droplet number concentration (specified)
!  pqm1     : specific humidity                                    (n-1)
!  ptm1     : temperature                                          (n-1)
!  pxlm1    : cloud liquid water                                   (n-1)
!  pxim1    : cloud ice                                            (n-1)
!  pxtec    : detrained convective cloud liquid water or cloud ice (n)
!  pxvar    : distribution width (b-a)                             (n-1)
!  pxskew   : beta shape parameter "q"                             (n-1)
!  pbetaa   : the beta distribution minimum a                      (n-1)
!  pbetab   : the beta distribution maximum b                      (n-1)
!  pvdiffp  : the rate of change of q due to vdiff scheme          (n-1)
!  phmixtau : mixing timescale**-1 for horizontal turbulence       (n)
!  pvmixtau : mixing timescale**-1 for horizontal turbulence       (n)
!
!  - 1D
!  knvb     :
!
!     Output arguments.
!     ------ ----------
!  - 1D
!  prsfl    : surface rain flux
!  pssfl    : surface snow flux
!
!     Input, Output arguments.
!     ------------ ----------
!  - 2D
!  paclc    : cloud cover  (now diagnosed in cover)
!  paclcac  : cloud cover, accumulated
!  paclcov  : total cloud cover
!  paprl    : total stratiform precipitation (rain+snow), accumulated
!  pqvi     : vertically integrated spec. humidity, accumulated
!  pxlvi    : vertically integrated cloud liquid water, accumulated
!  pxivi    : vertically integrated cloud ice, accumulated
!  ptte     : tendency of temperature
!  pqte     : tendency of specific humidity
!  pxlte    : tendency of cloud liquid water
!  pxite    : tendency of cloud ice
!  - 1D
!  paprs    : Snowfall, accumulated
!
!     Externals.
!     ----------
!
!     Method.
!     -------
!     see References
!
!     References.
!     ----------
!
!     Lohmann and Roeckner, 1996: Clim. Dyn. 557-572
!     Levkov et al., 1992: Beitr. Phys. Atm. 35-58.          (ice phase)
!     Beheng, 1994: Atmos. Res. 193-206.                    (warm phase)
!     Lenderink et al., 1998; KNMI-REPORT NO. 98-13       (condensation)
!     Tompkins 2002, J. Atmos. Sci.                        (cloud cover)
!
!     Authors.
!     -------
!     M.Esch        MPI-Hamburg  1999
!     G.Lenderink   KNMI, de Bilt 1998
!     U.Lohmann     MPI-Hamburg  1995
!
!     Modifications.
!     --------------
!     E.Roeckner    MPI-Hamburg  2000
!     A.Tompkins    MPI-Hamburg  2000
!     U.Schlese     MPI-Hamburg  2003
!
USE mo_kind,           ONLY : dp
USE mo_constants,      ONLY : cpd, vtmpc2, g, rd, alv, als, rv, api    &
                            , vtmpc1, rhoh2o
USE mo_convect_tables, ONLY : lookuperror, lookupoverflow, jptlucu1    &
                            , jptlucu2, tlucua, tlucub, tlucuaw
USE mo_param_switches, ONLY : lcover
USE mo_cloud,          ONLY : cqtmin, tmelt, cvtfall, crhosno, cn0s    &
                            , cthomi, csecfrl, ncctop, cvarmin         &
                            , cbeta_pq, cbeta_pq_max, nbetaq, cbetaqs  &
                            , rbetak, nbetax, tbetai0, tbetai1, cauloc &
                            , clmax, clmin, jbmin, jbmax, lonacc       &
                            , ccraut, ceffmin, ceffmax, crhoi, ccsaut  &
                            , ccsacl, cbeta_cs, ccwmin 
USE mo_time_control,   ONLY : delta_time, time_step_len

!---wiso-code

USE mo_wiso,           ONLY : wiso_frac_liq, wiso_frac_liq_ice         &
                            , tnat, nwisotyp                           &
                            , talphal1, talphal2, talphal3, talphal4   &
                            , talphas1, talphas2, talphas3, talphas4   &
                            , tsatbase, tsatfac, tdifrel               &
                            , cwisomin, cwisosec
!---wiso-code-end

!
  IMPLICIT NONE
!
  INTEGER :: kbdim, klevp1, klev, kproma, ktdia

!---wiso-code

  LOGICAL :: lwiso
  INTEGER :: kwiso

!---wiso-code-end

  REAL(dp):: paphm1(kbdim,klevp1) ,pvervel(kbdim,klev)                 &
           , papm1(kbdim,klev)    ,pqm1(kbdim,klev)                    &
           , papp1(kbdim,klev)    ,ptm1(kbdim,klev)                    &
           , ptvm1(kbdim,klev)    ,pxlm1(kbdim,klev)                   &
           , pxim1(kbdim,klev)    ,pxtec(kbdim,klev)                   &
           , pqtec(kbdim,klev)                                         &
           , pxvar(kbdim,klev)    ,pxskew(kbdim,klev)                  &
           , pbetaa(kbdim,klev)   ,pbetab(kbdim,klev)                  &
           , pvdiffp(kbdim,klev)  ,phmixtau(kbdim,klev)                &
           , pvmixtau(kbdim,klev) ,pgeo(kbdim,klev)                    &
           , pbetass(kbdim,klev)
  REAL(dp):: pxlvi(kbdim)         ,pxivi(kbdim)
  REAL(dp):: paclc(kbdim,klev)    ,paclcac(kbdim,klev)
  REAL(dp):: pacdnc(kbdim,klev)   ,prelhum(kbdim,klev)
  REAL(dp):: paclcov(kbdim)       ,paprl(kbdim)                        &
           , pqvi(kbdim)          ,pssfl(kbdim)
  REAL(dp):: ptte(kbdim,klev)     ,pqte(kbdim,klev)
  REAL(dp):: pxlte(kbdim,klev)    ,pxite(kbdim,klev)
  REAL(dp):: paprs(kbdim)         ,prsfl(kbdim)

!---wiso-code

  REAL(dp), OPTIONAL :: pwisoqm1(kbdim,klev,kwiso)                                &
                      , pwisoxlm1(kbdim,klev,kwiso)                               &
                      , pwisoxim1(kbdim,klev,kwiso), pwisoxtec(kbdim,klev,kwiso)  &
                      , pwisoqtec(kbdim,klev,kwiso)
  REAL(dp), OPTIONAL :: pwisoxlvi(kbdim,kwiso),      pwisoxivi(kbdim,kwiso)
  REAL(dp), OPTIONAL :: pwisoaprl(kbdim,kwiso)                                    &
                      , pwisoqvi(kbdim,kwiso),       pwisossfl(kbdim,kwiso)
  REAL(dp), OPTIONAL :: pwisoqte(kbdim,klev,kwiso)
  REAL(dp), OPTIONAL :: pwisoxlte(kbdim,klev,kwiso), pwisoxite(kbdim,klev,kwiso)
  REAL(dp), OPTIONAL :: pwisoaprs(kbdim,kwiso),      pwisorsfl(kbdim,kwiso)
  
!---wiso-code-end

!
!   Temporary arrays
!
  REAL(dp):: zclcpre(kbdim)                                            &
           , zcnd(kbdim)          ,zdep(kbdim)        ,zdp(kbdim)      &
           , zevp(kbdim)          ,zxievap(kbdim)     ,zxlevap(kbdim)  &
           , zfrl(kbdim)          ,zimlt(kbdim)       ,zsmlt(kbdim)    &
           , zrpr(kbdim)          ,zspr(kbdim)        ,zsub(kbdim)     &
           , zxlte(kbdim)         ,zxite(kbdim)       ,zxiflux(kbdim)  &
           , zsacl(kbdim)         ,zdz(kbdim)                          &
           , zlsdcp(kbdim)        ,zlvdcp(kbdim)      ,zximlt(kbdim)   &
           , ztp1tmp(kbdim)       ,zqp1tmp(kbdim)     ,zxisub(kbdim)   &
           , zrfl(kbdim)          ,zsfl(kbdim)                         &
           , zxlb(kbdim)          ,zxib(kbdim)                         &
           , zrho(kbdim,klev)     ,zclcov(kbdim)      ,zclcaux(kbdim)  &
           , zqvi(kbdim)          ,zxlvi(kbdim)       ,zxivi(kbdim)    &
           , zbetaqt(kbdim)       ,zwide(kbdim)                        &
           , zbetacl(kbdim)       ,zturbvar(kbdim)                     &
           , zturbskew(kbdim)                                          &
           , zconvvar(kbdim)      ,zconvskew(kbdim)   ,zvartg(kbdim)   &
           , zmicroskew(kbdim)    ,zgenti(kbdim)      ,zgentl(kbdim)   &
           , zxvarte(kbdim)       ,zxskewte(kbdim)                     &
           , zgeoh(kbdim,klevp1)
!
  REAL(dp):: zbqp1, zbbap1, zbap1, ztt, zgent, zqcdif, zqp1, ztp1      &
           , zqp1b, zskew, zbetai0, zbetai1, zskewp1, zvarp1, zifrac   &
           , zvarmx, zdtime, zauloc
  INTEGER:: iqidx, ixidx, jb
!  LOGICAL   lo, lo1, lo2, locc    !-changed for wiso_code
  LOGICAL   lo, lo1, lo2(kbdim), locc(kbdim)
  INTEGER   knvb(kbdim)
  INTEGER:: jl, jk, it, it1
  REAL(dp):: zepsec, zxsec, zqsec, ztmst, zcons1, zcons2, zrcp, zcons  &
           , ztdif, zsnmlt, zximelt, zclcstar, zdpg, zqrho, zesi, zqsi &
           , zsusati, zb1, zb2, zcoeff, zxsp1, zclambs, zcfac4c, zzeps &
           , zsubi, zxrp1, zesw, zesat, zqsw, zsusatw, zdv, zast, zbst &
!           , zzepr, zxip1, zxifall, zal1, zal2, zxised, zxim1evp       &  !---wiso-code: arrry dimension changed
           , zzepr, zxip1, zxifall, zal1, zal2, zxised(kbdim), zxim1evp(kbdim)       &
!           , zxlm1evp, zxidt, zxldt, zxidtstar, zxldtstar, zlc, zqsm1  &  !---wiso-code: arrry dimension changed
           , zxlm1evp(kbdim), zxidt, zxldt, zxidtstar(kbdim), zxldtstar(kbdim), zlc, zqsm1  &
           , zqst1, zdqsdt, zlcdqsdt, zqvdt, zdtdt, zdtdtstar, zdqsat  &
           , zxilb, zrelhum, zqtau, zpp, zqq, zeta, zprod, zaa, zes    &
           , zcor, zqsp1tmp, zoversat, zrhtest, zqcon, zdepcor         &
           , zdepos, zcond, zradl, zf1, zraut, zexm1, zexp, zrac1      &
           , zrac2, zrieff, zrih, zri, zcolleffi, zc1, zdt2, zsaut     &
           , zsaci1, zsaci2, zsacl1, zsacl2, zlamsm, zxsp2, zzdrr      &
           , zzdrs, zpretot, zpredel, zpresum, zmdelb, zmqp1, zxlp1    &
           , zxlold, zxiold, zdxicor, zdxlcor, zcndcor

!---wiso-code

  REAL(dp)::                                                           &
             zwisocnd(kbdim,kwiso),  zwisodep(kbdim,kwiso)             &
           , zwisoevp(kbdim, kwiso), zwisoxievap(kbdim,kwiso)          &
           , zwisoxlevap(kbdim,kwiso)                                  &
           , zwisofrl(kbdim,kwiso),  zwisoimlt(kbdim,kwiso)            &
           , zwisorpr(kbdim,kwiso),  zwisospr(kbdim,kwiso)             &
           , zwisosub(kbdim,kwiso)                                     &
           , zwisoxlte(kbdim,kwiso), zwisoxite(kbdim,kwiso)            &
           , zwisoxiflux(kbdim, kwiso)                                 &
           , zwisosacl(kbdim,kwiso)                                    &
           , zwisoximlt(kbdim,kwiso)                                   &
           , zwisoqp1tmp(kbdim,kwiso)                                  &
           , zwisoxisub(kbdim,kwiso)                                   &
           , zwisorfl(kbdim,kwiso),  zwisosfl(kbdim,kwiso)             &
           , zwisoxlb(kbdim,kwiso),  zwisoxib(kbdim,kwiso)             &
           , zwisoxilb(kbdim,kwiso)                                    &
           , zwisoqvi(kbdim,kwiso),  zwisoxlvi(kbdim,kwiso)            &
           , zwisoxivi(kbdim,kwiso)                                    &
           , zwisogenti(kbdim,kwiso), zwisogentl(kbdim,kwiso)          &
           , zwisoxim1evp(kbdim,kwiso), zwisoxlm1evp(kbdim,kwiso)      &
           , zwisoxidtstar(kbdim,kwiso), zwisoxldtstar(kbdim,kwiso)
           
  REAL(dp):: zwisofracliq(kbdim,kwiso), zwisofracice(kbdim,kwiso), zsatval

  REAL(dp):: zsfl_tmp(kbdim), zsnmlt_tmp(kbdim)                        &
           , zxip1_tmp(kbdim)                                          &
           , zxiflux_tmp(kbdim), zximelt_tmp(kbdim)                    &
           , zxiflux_tmp2(kbdim)                                       &
           , zzepr_tmp(kbdim), zqsw_tmp(kbdim)                         &
           , zdep_tmp(kbdim), zcnd_tmp(kbdim)                          &
           , zxlb_tmp(kbdim), zxib_tmp(kbdim)                          &
           , zxilb_tmp(kbdim)                                          &
           , zxlb_tmp2(kbdim),zxib_tmp2(kbdim)                         &
           , zxlb_tmp3(kbdim)                                          &
           , zqp1_tmp(kbdim), zqp1_tmp2(kbdim)                         &
           , zqp1_tmp3(kbdim)                                          &
           , ztp1_tmp(kbdim), ztp1_tmp2(kbdim)                         &
           , zqcdif_tmp(kbdim)                                         &
           , zqvdt_tmp(kbdim),ztp1tmp_tmp(kbdim)                       &
           , zzdrs_tmp(kbdim)                                          &
           , zdepcor_tmp(kbdim), zcndcor_tmp(kbdim)                    &
           , zfrl_tmp(kbdim), zfrl_tmp2(kbdim)                         &
           , zrpr_tmp(kbdim), zspr_tmp(kbdim)                          &
           , zxlb_difw(kbdim), zxib_difc(kbdim)                        &
           , zqsp1tmp_tmp(kbdim), zclcaux_tmp(kbdim)                   &
           , ztpone(kbdim), zqpone(kbdim), zqspone(kbdim)              &
           , zxlpone(kbdim)
  LOGICAL    lo_wiso(kbdim),  lo2_wiso(kbdim)                          &
           , lo3_wiso(kbdim), lo4_wiso(kbdim)                          &
           , lo_zqcdif(kbdim), lo2_tmp(kbdim)                          &
           , lo_dep(kbdim), lo_cnd(kbdim)                              &
           , lo_xlcor(kbdim), lo_xicor(kbdim)                          &
           , lo2_new(kbdim)                                            &
           , lo_zxldt(kbdim), lo_zxidt(kbdim)                          &
           , lo_depcond(kbdim)
  INTEGER::  jt
  REAL(dp):: zwisosnmlt, zwisoximelt                                   &
           , zwisoxip1,  zwisoxised                                    &
           , zwisoxidt, zwisoxldt                                      &
           , zwisozepr, zwisoqsw                                       &
           , zwisocndcor, zwisodepcor                                  &
           , zwisodepos, zwisocond                                     &
           , zwisozdrr, zwisozdrs                                      &
           , zwisoqp1, zwisoxlp1, zwisoxlold, zwisoxiold               &
           , zwisodxicor, zwisodxlcor                                  &
           , zwisoqvdt, zwisoqcdif                                     &
           , zdelta, zwisodenom, zwisonumer

!---wiso-code-end

!
! Executable statements
!
  lookupoverflow = .FALSE.
!
!   Security parameters
!
  zepsec = 1.0e-12_dp
  zxsec  = 1.0_dp-zepsec
  zqsec  = 1.0_dp-cqtmin
!
!   Computational constants
!
  zdtime = delta_time
  ztmst  = time_step_len
  zcons1 = cpd*vtmpc2
  zcons2 = 1._dp/(ztmst*g)
!
!     ------------------------------------------------------------------
!
!       1.   Top boundary conditions, air density and geopotential
!            height at half levels
!
!       1.1   Set to zero precipitation fluxes etc.
!
  DO 111 jl = 1,kproma
     zclcpre(jl)   = 0.0_dp
     zxiflux(jl)   = 0.0_dp
     zrfl(jl)      = 0.0_dp
     zsfl(jl)      = 0.0_dp
111 END DO

!---wiso-code
  IF (lwiso) THEN
  
  DO jt = 1,kwiso
    DO jl = 1,kproma
      zwisoxiflux(jl,jt)   = 0.0_dp
      zwisorfl(jl,jt)      = 0.0_dp
      zwisosfl(jl,jt)      = 0.0_dp
    END DO
  END DO

  END IF
!---wiso-code-end
!
!       1.2   Air density
!
  DO 122 jk = ktdia,klev
     DO 121 jl = 1,kproma
        zrho(jl,jk)   = papm1(jl,jk)/(rd*ptvm1(jl,jk))
        pxtec(jl,jk)  = MAX(pxtec(jl,jk),0.0_dp)
        pqtec(jl,jk)  = MAX(pqtec(jl,jk),0.0_dp)
121  END DO
122 END DO

!---wiso-code
  IF (lwiso) THEN

  DO jt = 1,kwiso
    DO jk = ktdia,klev
      DO jl = 1,kproma
        pwisoxtec(jl,jk,jt)  = MAX(pwisoxtec(jl,jk,jt),0.0_dp)
        pwisoqtec(jl,jk,jt)  = MAX(pwisoqtec(jl,jk,jt),0.0_dp)
      END DO
    END DO
  END DO

  END IF
!---wiso-code-end
!
!       1.3   Geopotential at half levels
!
  DO 132 jk = 2,klev
     DO 131 jl = 1,kproma
        zgeoh(jl,jk)   = 0.5_dp*(pgeo(jl,jk)+pgeo(jl,jk-1))
131  END DO
132 END DO
   DO 133 jl = 1,kproma
      zgeoh(jl,1)      = pgeo(jl,1)+(pgeo(jl,1)-zgeoh(jl,2))
      zgeoh(jl,klevp1) = 0.0_dp
133 END DO
!
  DO 831 jk=ktdia,klev

!
!     ------------------------------------------------------------------
!
!       2.    Set to zero local tendencies (increments)
!
     DO 201 jl = 1,kproma
        zcnd(jl)       = 0.0_dp
        zdep(jl)       = 0.0_dp
        zfrl(jl)       = 0.0_dp
        zrpr(jl)       = 0.0_dp
        zspr(jl)       = 0.0_dp
        zimlt(jl)      = 0.0_dp
        zevp(jl)       = 0.0_dp
        zsub(jl)       = 0.0_dp
        zximlt(jl)     = 0.0_dp
        zxisub(jl)     = 0.0_dp
        zsmlt(jl)      = 0.0_dp
        zsacl(jl)      = 0.0_dp
        zgenti(jl)     = 0.0_dp
        zgentl(jl)     = 0.0_dp
        zxievap(jl)    = 0.0_dp
        zxlevap(jl)    = 0.0_dp
        zvartg(jl)     = 0.0_dp
        zconvvar(jl)   = 0.0_dp
        zconvskew(jl)  = 0.0_dp
        zturbvar(jl)   = 0.0_dp
        zturbskew(jl)  = 0.0_dp
        zmicroskew(jl) = 0.0_dp
        zxvarte(jl)    = 0.0_dp
        zxskewte(jl)   = 0.0_dp   
        zdp(jl)        = paphm1(jl,jk+1)-paphm1(jl,jk)
        zdz(jl)        = (zgeoh(jl,jk)-zgeoh(jl,jk+1))/g
        zrcp           = 1._dp/(cpd+zcons1*MAX(pqm1(jl,jk),0.0_dp))
        zlvdcp(jl)     = alv*zrcp
        zlsdcp(jl)     = als*zrcp
201  END DO

!---wiso-code
  IF (lwiso) THEN

    DO jt = 1,kwiso
      DO jl = 1,kproma
        zwisocnd(jl,jt)    = 0.0_dp
        zwisodep(jl,jt)    = 0.0_dp
        zwisofrl(jl,jt)    = 0.0_dp
        zwisorpr(jl,jt)    = 0.0_dp
        zwisospr(jl,jt)    = 0.0_dp
        zwisoimlt(jl,jt)   = 0.0_dp
        zwisoevp(jl,jt)    = 0.0_dp
        zwisosub(jl,jt)    = 0.0_dp
        zwisoximlt(jl,jt)  = 0.0_dp
        zwisoxisub(jl,jt)  = 0.0_dp
        zwisosacl(jl,jt)   = 0.0_dp
        zwisogenti(jl,jt)  = 0.0_dp
        zwisogentl(jl,jt)  = 0.0_dp
        zwisoxievap(jl,jt) = 0.0_dp
        zwisoxlevap(jl,jt) = 0.0_dp
      END DO
     END DO

  END IF
!---wiso-code-end
!
!     ------------------------------------------------------------------
!
!       3.   Modification of incoming precipitation fluxes by
!            melting, sublimation and evaporation
!
     IF (jk .GT. 1) THEN
!
!DIR$ CONCURRENT
        DO 331 jl = 1,kproma
!
!       3.1   Melting of snow and ice
!
           zcons     = zcons2*zdp(jl)/(zlsdcp(jl)-zlvdcp(jl))
           ztdif     = MAX(0.0_dp,ptm1(jl,jk)-tmelt)
           zsnmlt    = MIN(zxsec*zsfl(jl),zcons*ztdif)
           zsfl_tmp(jl)   = zsfl(jl)                                 !---wiso-code: store old value of snowfall
           zsnmlt_tmp(jl) = zsnmlt                                   !---wiso-code: store amount of snowmelt
           zrfl(jl)  = zrfl(jl)+zsnmlt
           zsfl(jl)  = zsfl(jl)-zsnmlt
           zsmlt(jl) = zsnmlt/(zcons2*zdp(jl))
           zximelt   = MIN(zxsec*zxiflux(jl),zcons*ztdif)
           zxiflux_tmp(jl) = zxiflux(jl)                             !---wiso-code: store old value of ice flux           
           zximelt_tmp(jl) = zximelt                                 !---wiso-code: store amount of ice melt           
           zxiflux(jl)=zxiflux(jl)-zximelt
           zximlt(jl) =zximelt/(zcons2*zdp(jl))
           IF (ztdif.GT.0.0_dp) THEN
            zimlt(jl) = MAX(0.0_dp,pxim1(jl,jk)+pxite(jl,jk)*ztmst)
           ELSE
            zimlt(jl) = 0.0_dp
           END IF
!
!       3.2   Sublimation of snow and ice (Lin et al., 1983)
!
           IF (zclcpre(jl) .GT. 0.0_dp) THEN
              zclcstar = zclcpre(jl)
              zdpg     = zdp(jl)/g
              zqrho    = 1.3_dp/zrho(jl,jk)
              it       = NINT(ptm1(jl,jk)*1000._dp)
              IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
              it = MAX(MIN(it,jptlucu2),jptlucu1)
              zesi     = tlucua(it)/papm1(jl,jk)
              zesi     = MIN(zesi,0.5_dp)
              zqsi     = zesi/(1._dp-vtmpc1*zesi)
              zsusati  = MIN(pqm1(jl,jk)/zqsi-1.0_dp,0.0_dp)
              zb1      = zlsdcp(jl)**2/(2.43e-2_dp*rv*(ptm1(jl,jk)**2))
              zb2      = 1._dp/(zrho(jl,jk)*zqsi*0.211e-4_dp)
              zcoeff   = 3.e6_dp*2._dp*api*zsusati/                    &
                                          (zrho(jl,jk)*(zb1+zb2))
!
              IF (zsfl(jl) .GT. cqtmin) THEN
               zxsp1    = (zsfl(jl)/                                   &
                               (zclcpre(jl)*cvtfall))**(1._dp/1.16_dp)
               zclambs  = (zxsp1/(api*crhosno*cn0s))**0.25_dp
               zcfac4c  = 0.78_dp*zclambs**2+232.19_dp*zqrho**0.25_dp  &
                                                    *zclambs**2.625_dp
               zzeps    = MAX(-zxsec*zsfl(jl)/zclcpre(jl),             &
                                                  zcoeff*zcfac4c*zdpg)
               zsub(jl) = -zzeps/zdpg*ztmst*zclcstar
               lo_wiso(jl) = zsub(jl).LT.MAX(zxsec*(zqsi-pqm1(jl,jk)),0.0_dp)  !---wiso-code: store water source of zsub
               zsub(jl) = MIN(zsub(jl),                                &
                                  MAX(zxsec*(zqsi-pqm1(jl,jk)),0.0_dp))
               zsub(jl) = MAX(zsub(jl),0.0_dp)
              END IF
!
              IF (zxiflux(jl) .GT. cqtmin) THEN
               zxsp1    = (zxiflux(jl)/                                &
                               (zclcpre(jl)*cvtfall))**(1._dp/1.16_dp)
               zclambs  = (zxsp1/(api*crhosno*cn0s))**0.25_dp
               zcfac4c  = 0.78_dp*zclambs**2+232.19_dp*zqrho**0.25_dp  &
                                                    *zclambs**2.625_dp
               zzeps    = MAX(-zxsec*zxiflux(jl)/zclcpre(jl),          &
                                                  zcoeff*zcfac4c*zdpg)
               zsubi    = -zzeps/zdpg*ztmst*zclcstar
               lo2_wiso(jl) = zsubi.LT.MAX(zxsec*(zqsi-pqm1(jl,jk)),0.0_dp)    !---wiso-code: store water source of zxisub
               zsubi    = MIN(zsubi,                                   &
                                  MAX(zxsec*(zqsi-pqm1(jl,jk)),0.0_dp))
               zsubi    = MAX(zsubi,0.0_dp)
               zxiflux_tmp2(jl) = zxiflux(jl)                                  !---wiso-code: store old value of ice flux
               zxiflux(jl)=zxiflux(jl)-zsubi*zcons2*zdp(jl)
               zxisub(jl) =zsubi
              END IF
           END IF
!
!       3.3   Evaporation of rain (Rotstayn, 1997)
!
           IF (zclcpre(jl) .GT. 0.0_dp .AND. zrfl(jl) .GT. cqtmin) THEN
              zclcstar = zclcpre(jl)
              zdpg     = zdp(jl)/g
              zqrho    = 1.3_dp/zrho(jl,jk)
              zxrp1    = (zrfl(jl)/(zclcpre(jl)*12.45_dp*SQRT(zqrho))) &
                                                       **(8._dp/9._dp)
              it       = NINT(ptm1(jl,jk)*1000._dp)
              IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
              it = MAX(MIN(it,jptlucu2),jptlucu1)
              zesw     = tlucuaw(it)/papm1(jl,jk)
              zesat    = zesw*papm1(jl,jk)*rv/rd
              zesw     = MIN(zesw,0.5_dp)
              zqsw     = zesw/(1._dp-vtmpc1*zesw)
              zsusatw  = MIN(pqm1(jl,jk)/zqsw-1.0_dp,0.0_dp)
              zdv      = 2.21_dp/papm1(jl,jk)
              zast     = alv*(alv/(rv*ptm1(jl,jk))-1.0_dp)/            &
                                            (0.024_dp*ptm1(jl,jk))
              zbst     = rv*ptm1(jl,jk)/(zdv*zesat)
              zzepr    = 870._dp*zsusatw*(zrfl(jl)/                    &
                                             zclcpre(jl))**0.61_dp     &
                                     /(SQRT(zrho(jl,jk))*(zast+zbst))
              zzepr    = MAX(-zxsec*zrfl(jl)/zclcpre(jl),zzepr*zdpg)
              zzepr_tmp(jl) = -zzepr/zdpg*ztmst                      !---wiso-code: store rainfall amount available for evaporation
              zevp(jl) = -zzepr/zdpg*ztmst*zclcstar
              zqsw_tmp(jl) = MAX(zxsec*(zqsw-pqm1(jl,jk)),0.0_dp)    !---wiso-code: store rainfall amount available for evaporation
              lo3_wiso(jl) = zevp(jl).LT.MAX(zxsec*(zqsw-pqm1(jl,jk)),0.0_dp)   !---wiso-code: store water source of zevp
              zevp(jl) = MIN(zevp(jl),                                 &
                                  MAX(zxsec*(zqsw-pqm1(jl,jk)),0.0_dp))
              zevp(jl) = MAX(zevp(jl),0.0_dp)
           END IF
331     END DO
!

!---wiso-code
  IF (lwiso) THEN

!       Caculate fractionation factors for isotopic changes during evaporation of rain
        CALL wiso_frac_liq(kproma,kbdim,kwiso,ptm1(:,jk),zwisofracliq)

        DO jt = 1,kwiso
          DO jl = 1,kproma
          
!       Melting of snow and ice for water isotopes
!       - assume no fractionation during melting due to low diffusivity of water isotopes in ice
 
            ztdif     = MAX(0.0_dp,ptm1(jl,jk)-tmelt)
            zdelta=tnat(jt)
            IF (zsfl_tmp(jl).GT.cwisomin) zdelta = zwisosfl(jl,jt)/zsfl_tmp(jl)
            IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta = 1._dp
            zwisosnmlt       = zsnmlt_tmp(jl)*zdelta
            zwisorfl(jl,jt)  = zwisorfl(jl,jt)+zwisosnmlt
            zwisosfl(jl,jt)  = zwisosfl(jl,jt)-zwisosnmlt
            zdelta=tnat(jt)
            IF (zxiflux_tmp(jl).GT.cwisomin) zdelta = zwisoxiflux(jl,jt)/zxiflux_tmp(jl)
            IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta = 1._dp
            zwisoximelt        = zximelt_tmp(jl)*zdelta
            zwisoxiflux(jl,jt) = zwisoxiflux(jl,jt)-zwisoximelt
            zwisoximlt(jl,jt)  = zwisoximelt/(zcons2*zdp(jl))
            IF (ztdif.GT.0.0_dp) THEN
              zwisoimlt(jl,jt) = MAX(0.0_dp,pwisoxim1(jl,jk,jt)+pwisoxite(jl,jk,jt)*ztmst)
            ELSE
             zwisoimlt(jl,jt)  = 0.0_dp
            END IF

!       Sublimation of snow and ice - water isotopes 
!       - assume no fractionation during sublimation due to low diffusivity of water isotopes in ice

            IF (zclcpre(jl) .GT. 0.0_dp) THEN
              IF (zsfl(jl) .GT. cqtmin) THEN
                zdelta=tnat(jt)
                IF (lo_wiso(jl)) THEN
                  IF (zsfl(jl).GT.cwisomin) zdelta = zwisosfl(jl,jt)/zsfl(jl)
                ELSE
                  IF (pqm1(jl,jk).GT.cwisomin) zdelta = pwisoqm1(jl,jk,jt)/pqm1(jl,jk)
                ENDIF
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisosub(jl,jt) = zsub(jl)*zdelta
              END IF
              IF (zxiflux_tmp2(jl) .GT. cqtmin) THEN
                zdelta=tnat(jt)
                IF (lo2_wiso(jl)) THEN
                  IF (zxiflux_tmp2(jl).GT.cwisomin) zdelta = zwisoxiflux(jl,jt)/zxiflux_tmp2(jl)
                ELSE
                  IF (pqm1(jl,jk).GT.cwisomin) zdelta = pwisoqm1(jl,jk,jt)/pqm1(jl,jk)
                ENDIF
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisoxisub(jl,jt)  = zxisub(jl)*zdelta
                zwisoxiflux(jl,jt) = zwisoxiflux(jl,jt)-zwisoxisub(jl,jt)*zcons2*zdp(jl)
              END IF
            END IF

!       Evaporation of rain - water isotopes
!       - assume fractionation during liquid-vapour phase change occurs within a closed system

            IF (zclcpre(jl) .GT. 0.0_dp .AND. zrfl(jl) .GT. cqtmin) THEN
              If (lo3_wiso(jl)) THEN
                ! calculate isotope amount of rainfall available for evaporation (assume same delta as in rainfall rfl)
                zdelta=tnat(jt)
                if (zrfl(jl).GT.cwisomin) zdelta = zwisorfl(jl,jt)/zrfl(jl)
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisozepr = zzepr_tmp(jl) * zdelta
                zwisodenom=zwisofracliq(jl,jt)*zzepr_tmp(jl)+(1._dp-zwisofracliq(jl,jt))*zevp(jl)
                zdelta=tnat(jt)
                IF (ABS(zwisodenom).GT.cwisomin) zdelta=zwisozepr/zwisodenom
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisoevp(jl,jt) = zevp(jl)*zdelta
              ELSE
                ! calculate isotope amount of rainfall available for evaporation (assume same delta as in pqm1)
                zdelta=tnat(jt)
                if (pqm1(jl,jk).GT.cwisomin) zdelta = pwisoqm1(jl,jk,jt)/pqm1(jl,jk)
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisoqsw = zqsw_tmp(jl) * zdelta
                zwisodenom=zwisofracliq(jl,jt)*zqsw_tmp(jl)+(1._dp-zwisofracliq(jl,jt))*zevp(jl)
                zdelta=tnat(jt)
                IF (ABS(zwisodenom).GT.cwisomin) zdelta=zwisoqsw/zwisodenom
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisoevp(jl,jt) = zevp(jl)*zdelta
              ENDIF
              
            END IF

          END DO
        END DO

  END IF
!---wiso-code-end

        IF (lookupoverflow) CALL lookuperror ('cloud (1)    ')
!
     END IF
!
!DIR$ CONCURRENT
     DO 610 jl=1,kproma
!
!     ------------------------------------------------------------------
!       4.    Sedimentation of cloud ice from grid-mean values.
!             Updating the tendency 'pxite' to include sedimentation.
!             At jk=klev, the sedimentation sink is balanced by
!             precipitation at the surface (through 'zzdrs', see 7.3).
!             Finally: In-cloud cloud water/ice.
!
        zxip1         = pxim1(jl,jk)+pxite(jl,jk)*ztmst-zimlt(jl)
        zxip1         = MAX(zxip1,EPSILON(1._dp))
        zxip1_tmp(jl) = zxip1                                        !---wiso-code: store value of zxip1
        zxifall       = cvtfall*(zrho(jl,jk)*zxip1)**0.16_dp
        zal1          = zxifall*g*zrho(jl,jk)*ztmst/zdp(jl)
        zal2          = zxiflux(jl)/(zrho(jl,jk)*zxifall)
        zxised(jl)    = zxip1*EXP(-zal1)+zal2*(1._dp-EXP(-zal1))     !---wiso-code: arrry dimension changed
        zxiflux(jl)   = zxiflux(jl)+(zxip1-zxised(jl))*zcons2*zdp(jl)!---wiso-code: arrry dimension changed
        pxite(jl,jk)  = (zxised(jl)-pxim1(jl,jk))/ztmst              !---wiso-code: arrry dimension changed
!
!             In-cloud water/ice calculated from respective grid-means,
!             partial cloud cover, advective/diffusive tendencies,
!             detrained cloud water/ice and ice sedimentation.
!             In-cloud values are required for cloud microphysics.
!
!---wiso-code: change array dimensions of several variables in the following IF-block
        zclcaux(jl) = paclc(jl,jk)
        zclcaux_tmp(jl) = zclcaux(jl)                                !---wiso-code: store old value of zclcaux(jl)
        locc(jl)    = zclcaux(jl) .GT. 0.0_dp
        lo2(jl)     = (ptm1(jl,jk) .LT. cthomi) .OR.                   &
                      (ptm1(jl,jk) .LT. tmelt .AND. zxised(jl) .GT. csecfrl)
        IF (lo2(jl)) THEN                 !     ice cloud
           zxite(jl)          = pxtec(jl,jk)
           zxlte(jl)          = 0.0_dp
           IF (locc(jl)) THEN
              zxib(jl)        = pxim1(jl,jk)/zclcaux(jl)
              zxlb(jl)        = pxlm1(jl,jk)/zclcaux(jl)
              zxim1evp(jl)    = 0.0_dp
              zxlm1evp(jl)    = 0.0_dp
              zxidt           = (pxite(jl,jk)+zxite(jl))*ztmst
              zxldt           =  pxlte(jl,jk)*ztmst+zximlt(jl)+zimlt(jl)
              IF (zxidt .GT. 0.0_dp) THEN
                 lo_zxidt(jl) = .TRUE.                               !---wiso-code: remember if zxidt is positive/negative 
                 zxidtstar(jl)    = zxidt
                 zxib(jl)     = zxib(jl)+zxidt
              ELSE
                 lo_zxidt(jl) = .FALSE.                              !---wiso-code: remember if zxidt is positive/negative 
                 zxidtstar(jl)    = 0.0_dp
                 zxib(jl)     = zxib(jl)+MAX(zxidt/zclcaux(jl),        &
                                                       -zxib(jl))
                 pxite(jl,jk) = MAX(pxite(jl,jk),                      &
                                 -(pxim1(jl,jk)/ztmst+zxite(jl)))
              END IF
              IF (zxldt .GT. 0.0_dp) THEN
                 lo_zxldt(jl) = .TRUE.                               !---wiso-code: remember if zxldt is positive/negative 
                 zxldtstar(jl)    = zxldt
                 zxlb(jl)     = zxlb(jl)+zxldt
              ELSE
                 lo_zxldt(jl) = .FALSE.                              !---wiso-code: remember if zxldt is positive/negative 
                 zxldtstar(jl)    = 0.0_dp
                 zxlb(jl)     = zxlb(jl)+MAX(zxldt/zclcaux(jl),        &
                                                       -zxlb(jl))
                 pxlte(jl,jk) = MAX(pxlte(jl,jk),-pxlm1(jl,jk)/ztmst)
              END IF
           ELSE                      !    cloud cover = 0.
              zxib(jl)        = 0.0_dp
              zxlb(jl)        = 0.0_dp
              zxidt           = (pxite(jl,jk)+zxite(jl))*ztmst
              zxldt           =  pxlte(jl,jk)*ztmst+zximlt(jl)+zimlt(jl)
              IF (zxidt .GT. 0.0_dp) THEN
                 lo_zxidt(jl) = .TRUE.                               !---wiso-code: remember if zxidt is positive/negative 
                 zxidtstar(jl)    = zxidt
                 zxim1evp(jl)     = pxim1(jl,jk)
              ELSE
                 lo_zxidt(jl) = .FALSE.                              !---wiso-code: remember if zxidt is positive/negative 
                 zxidtstar(jl)    = 0.0_dp
                 pxite(jl,jk) = MAX(pxite(jl,jk),                      &
                                  -(pxim1(jl,jk)/ztmst+zxite(jl)))
                 zxim1evp(jl)     = pxim1(jl,jk)+(pxite(jl,jk)+zxite(jl))  &
                                                                 *ztmst
              END IF
              IF (zxldt .GT. 0.0_dp) THEN
                 lo_zxldt(jl) = .TRUE.                               !---wiso-code: remember if zxldt is positive/negative 
                 zxldtstar(jl)    = zxldt
                 zxlm1evp(jl)     = pxlm1(jl,jk)
              ELSE
                 lo_zxldt(jl) = .FALSE.                              !---wiso-code: remember if zxldt is positive/negative 
                 zxldtstar(jl)    = 0.0_dp
                 pxlte(jl,jk) = MAX(pxlte(jl,jk),-pxlm1(jl,jk)/ztmst)
                 zxlm1evp(jl)     = pxlm1(jl,jk)+pxlte(jl,jk)*ztmst
              END IF
           END IF
        ELSE                           !    water cloud
           zxlte(jl)          = pxtec(jl,jk)
           zxite(jl)          = 0.0_dp
           IF (locc(jl)) THEN
              zxlb(jl)        = pxlm1(jl,jk)/zclcaux(jl)
              zxib(jl)        = pxim1(jl,jk)/zclcaux(jl)
              zxlm1evp(jl)        = 0.0_dp
              zxim1evp(jl)        = 0.0_dp
              zxldt           = (pxlte(jl,jk)+zxlte(jl))*ztmst         &
                                             +zximlt(jl)+zimlt(jl)
              zxidt           =  pxite(jl,jk)*ztmst
              IF (zxldt .GT. 0.0_dp) THEN
                 lo_zxldt(jl) = .TRUE.                               !---wiso-code: remember if zxldt is positive/negative 
                 zxldtstar(jl)    = zxldt
                 zxlb(jl)     = zxlb(jl)+zxldt
              ELSE
                 lo_zxldt(jl) = .FALSE.                              !---wiso-code: remember if zxldt is positive/negative 
                 zxldtstar(jl)    = 0.0_dp
                 zxlb(jl)     = zxlb(jl)+MAX(zxldt/zclcaux(jl),        &
                                                          -zxlb(jl))
                 pxlte(jl,jk) = MAX(pxlte(jl,jk),                      &
                                 -(pxlm1(jl,jk)/ztmst+zxlte(jl)))
              END IF
              IF (zxidt .GT. 0.0_dp) THEN
                 lo_zxidt(jl) = .TRUE.                               !---wiso-code: remember if zxidt is positive/negative 
                 zxidtstar(jl)    = zxidt
                 zxib(jl)     = zxib(jl)+zxidt
              ELSE
                 lo_zxidt(jl) = .FALSE.                              !---wiso-code: remember if zxidt is positive/negative 
                 zxidtstar(jl)    = 0.0_dp
                 zxib(jl)     = zxib(jl)+MAX(zxidt/zclcaux(jl),        &
                                                          -zxib(jl))
                 pxite(jl,jk) = MAX(pxite(jl,jk),-pxim1(jl,jk)/ztmst)
              END IF
           ELSE                          !    cloud cover = 0.
              zxlb(jl)        = 0.0_dp
              zxib(jl)        = 0.0_dp
              zxldt           = (pxlte(jl,jk)+zxlte(jl))*ztmst         &
                                             +zximlt(jl)+zimlt(jl)
              zxidt           =  pxite(jl,jk)*ztmst
              IF (zxldt .GT. 0.0_dp) THEN
                 lo_zxldt(jl) = .TRUE.                               !---wiso-code: remember if zxldt is positive/negative 
                 zxldtstar(jl)    = zxldt
                 zxlm1evp(jl)     = pxlm1(jl,jk)
              ELSE
                 lo_zxldt(jl) = .FALSE.                              !---wiso-code: remember if zxldt is positive/negative 
                 zxldtstar(jl)    = 0.0_dp
                 pxlte(jl,jk) = MAX(pxlte(jl,jk),                      &
                                  -(pxlm1(jl,jk)/ztmst+zxlte(jl)))
                 zxlm1evp(jl)     = pxlm1(jl,jk)+(pxlte(jl,jk)+zxlte(jl))  &
                                                                *ztmst
              END IF
              IF (zxidt .GT. 0.0_dp) THEN
                 lo_zxidt(jl) = .TRUE.                               !---wiso-code: remember if zxidt is positive/negative 
                 zxidtstar(jl)    = zxidt
                 zxim1evp(jl)     = pxim1(jl,jk)
              ELSE
                 lo_zxidt(jl) = .FALSE.                              !---wiso-code: remember if zxidt is positive/negative 
                 zxidtstar(jl)    = 0.0_dp
                 pxite(jl,jk) = MAX(pxite(jl,jk),-pxim1(jl,jk)/ztmst)
                 zxim1evp(jl)     = pxim1(jl,jk)+pxite(jl,jk)*ztmst
              END IF
           END IF
        END IF
!
!     ------------------------------------------------------------------
!       5.    Condensation/deposition and evaporation/sublimation
!
!             zlc       =  L_{v/s} / c_p
!             zlcdqsdt  = L dq_sat / c_p dT
!             zdqsdt    = dq_sat / dT
!
        zrcp        = 1._dp/(cpd+zcons1*MAX(pqm1(jl,jk),0.0_dp))
        zlvdcp(jl)  = alv*zrcp
        zlsdcp(jl)  = als*zrcp
        zlc         = MERGE(zlsdcp(jl),zlvdcp(jl),lo2(jl))
        it          = NINT(ptm1(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsm1       = MERGE(tlucua(it),tlucuaw(it),lo2(jl))/papm1(jl,jk)
        zqsm1       = MIN(zqsm1,0.5_dp)
        zqsm1       = zqsm1/(1._dp-vtmpc1*zqsm1)
        it1         = it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1       = MERGE(tlucua(it1),tlucuaw(it1),lo2(jl))/papm1(jl,jk)
        zqst1       = MIN(zqst1,0.5_dp)
        zqst1       = zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt      = (zqst1-zqsm1)*1000._dp
        zlcdqsdt    = zlc*zdqsdt
        zxievap(jl) = (1.0_dp-zclcaux(jl))*zxidtstar(jl)+zxim1evp(jl)
        zxlevap(jl) = (1.0_dp-zclcaux(jl))*zxldtstar(jl)+zxlm1evp(jl)
        zqvdt       = pqte(jl,jk)*ztmst+zevp(jl)+zsub(jl)              &
                             +zxievap(jl)+zxlevap(jl)+zxisub(jl)
        zdtdt       = ptte(jl,jk)*ztmst-zlvdcp(jl)*(zevp(jl)           &
                       +zxlevap(jl))                                   &
                       -zlsdcp(jl)*(zsub(jl)+zxievap(jl)+zxisub(jl))   &
                       -(zlsdcp(jl)-zlvdcp(jl))                        &
                       *(zsmlt(jl)+zximlt(jl)+zimlt(jl))
        zqp1        = MAX(pqm1(jl,jk)+zqvdt,0.0_dp)
        ztp1        = ptm1(jl,jk)+zdtdt
        zdtdtstar   = zdtdt+zclcaux(jl)*(zlc*pqte(jl,jk)*ztmst         &
                           +zlvdcp(jl)*(zevp(jl)+zxlevap(jl))          &
                         +zlsdcp(jl)*(zsub(jl)+zxievap(jl)+zxisub(jl)))
        zdqsat      = zdtdtstar*zdqsdt/(1._dp+zclcaux(jl)*zlcdqsdt)
        zxib(jl)    = MAX(zxib(jl),0.0_dp)
        zxlb(jl)    = MAX(zxlb(jl),0.0_dp)
        zxilb       = zxib(jl)+zxlb(jl)

!---wiso-code: store temporary values
        zqvdt_tmp(jl) = zqvdt                                                  !---wiso-code: store value of zqvdt   
        zqp1_tmp(jl) = zqp1                                                    !---wiso-code: store value of zqp1   
        ztp1_tmp(jl) = ztp1                                                    !---wiso-code: store value of ztp1   
        zxlb_tmp(jl) = zxlb(jl)                                                !---wiso-code: store old value of zxlb(jl)   
        zxib_tmp(jl) = zxib(jl)                                                !---wiso-code: store old value of zxib(jl)   
        zxilb_tmp(jl) = zxilb                                                  !---wiso-code: store old value of zxilb(jl)   
!---wiso-code-end

!
!       Diagnostics: relative humidity
!
        zrelhum        = pqm1(jl,jk)/zqsm1
        zrelhum        = MAX(MIN(zrelhum,1._dp),0._dp)
        prelhum(jl,jk) = zrelhum
!
        IF (lcover .AND. jk.GE.ncctop) THEN
!
!       define variables needed for cover scheme
!
!       zbetaqt = total water
!       zbetass = saturation mixing ratio adjusted to match qv
!       zwide   = current diagnosed distribution width
!
           zbetacl(jl) = MAX(0.0_dp,pxlm1(jl,jk))+                     &
                                               MAX(0.0_dp,pxim1(jl,jk))
           zbetaqt(jl) = MAX(cqtmin,pqm1(jl,jk))+zbetacl(jl)
           zvartg(jl)  = MAX(cqtmin,cvarmin*pqm1(jl,jk))
           zwide(jl)   = MAX(zvartg(jl),pbetab(jl,jk)-pbetaa(jl,jk))
           zskew       = MAX(MIN(pxskew(jl,jk),cbeta_pq_max),cbeta_pq)
           iqidx       = INT((nbetaq/cbetaqs)*                         &
                             LOG((zskew-cbeta_pq)/rbetak+1._dp)+0.5_dp)
!
!
!       5.1 Turbulence: Skewness - equation solved implicitly
!           This solver only works if phmixtau has non-zero timescale
!
           zqtau         = phmixtau(jl,jk)+pvmixtau(jl,jk)
           zbqp1         = cbeta_pq-(cbeta_pq-pxskew(jl,jk))           &
                                                    *EXP(-zqtau*zdtime)
           zbqp1         = MAX(MIN(zbqp1,cbeta_pq_max),cbeta_pq)
           zturbskew(jl) = (zbqp1-pxskew(jl,jk))/zdtime
!
!       5.2 Turbulence: variance - equation solved implicitly
!
           zpp          = cbeta_pq
           zqq          = pxskew(jl,jk)
           zeta         = (zpp+zqq)**2*(zpp+zqq+1._dp)/(zpp*zqq)
           zprod        = zeta*pvdiffp(jl,jk)/zwide(jl)
           zbbap1       = zprod/zqtau+zvartg(jl)-(zprod/zqtau          &
                            +zvartg(jl)-zwide(jl))*EXP(-zqtau*zdtime)
           zbbap1       = MAX(zbbap1,zvartg(jl))
           zbbap1       = MIN(zbbap1,zbetaqt(jl)*(cbeta_pq+zbqp1)      &
                                                            /cbeta_pq)
           zturbvar(jl) = (zbbap1-zwide(jl))/zdtime
           zbap1        = zbetaqt(jl)-zbbap1*cbeta_pq/(cbeta_pq+zbqp1)
!
!              translated into apparent xl,xi,q and heat sources
!              first order effect only, effect of evaporation of
!              cloud on qsat taken into account in thermodynamic budget
!              but does not change the mixing term here since that
!              would require iteration and is therefore neglected
!
!              calculate values after one timestep
!
           iqidx   = INT((nbetaq/cbetaqs)*                             &
                     LOG((zbqp1-cbeta_pq)/rbetak+1._dp)+0.5_dp)
           ztt     = (pbetass(jl,jk)-zbap1)*cbeta_pq/                  &
                     ((zbetaqt(jl)-zbap1)*(cbeta_pq+zbqp1))
           ztt     = nbetax*MAX(MIN(ztt,1.0_dp),0.0_dp)
           ixidx   = INT(ztt)
           IF (ixidx == nbetax) THEN
              zbetai0 = 1.0_dp
              zbetai1 = 1.0_dp
           ELSE
              zbetai0 = (ztt-ixidx)*tbetai0(iqidx,ixidx+1)             &
                             +(ixidx+1-ztt)*tbetai0(iqidx,ixidx)
              zbetai1 = (ztt-ixidx)*tbetai1(iqidx,ixidx+1)             &
                             +(ixidx+1-ztt)*tbetai1(iqidx,ixidx)
           ENDIF
           zqp1b      = (zbetaqt(jl)-zbap1)*zbetai1 -                  &
                        (pbetass(jl,jk)-zbap1)*zbetai0 + pbetass(jl,jk)
           lo_wiso(jl)= (pqm1(jl,jk)-zqp1b).GT.(-zxilb*zclcaux(jl))            !---wiso-code: remember water source of zgent
           zgent      = MAX(pqm1(jl,jk)-zqp1b,-zxilb*zclcaux(jl))              !              (true=pqm1-zqp1b, false=zxilb)
           lo2_wiso(jl)= (zgent.GT.zqsec*zqp1)                                 !---wiso-code: remember water source of zgent
           zgent      = MIN(zgent,zqsec*zqp1)              ! limit to qv       !              (true=zqp1)
           zifrac     = zxib(jl)/MAX(zepsec,zxilb)
           zifrac     = MAX(MIN(zifrac,1.0_dp),0.0_dp)
           zgenti(jl) = zgent*zifrac
           zgentl(jl) = zgent*(1.0_dp-zifrac)
           IF (locc(jl)) THEN
              zxib(jl) = MAX(zxib(jl)+zgenti(jl)/zclcaux(jl),0.0_dp)
              zxlb(jl) = MAX(zxlb(jl)+zgentl(jl)/zclcaux(jl),0.0_dp)
           END IF
           zxilb       = zxib(jl)+zxlb(jl)
!
!       5.3 Deposition/sublimation of cloud ice and condensation/
!           evaporation of liquid water due to changes in water vapour
!           and temperature (advection, convection, turbulent mixing,
!           evaporation of rain, sublimation and melting of snow).
!           Translate PDF laterally to calculate cloud
!           after one timestep
!
           zqvdt     = zqvdt-zgent
           zdtdt     = zdtdt+zlvdcp(jl)*zgentl(jl)+zlsdcp(jl)*zgenti(jl)
           zqp1      = MAX(pqm1(jl,jk)+zqvdt,0.0_dp)
           ztp1      = ptm1(jl,jk)+zdtdt
           zdtdtstar = zdtdt+zclcaux(jl)*(zlc*pqte(jl,jk)*ztmst        &
              +zlvdcp(jl)*(zevp(jl)+zxlevap(jl)-zgentl(jl))            &
              +zlsdcp(jl)*(zsub(jl)+zxievap(jl)+zxisub(jl)-zgenti(jl)))
           zdqsat    = zdtdtstar*zdqsdt/(1._dp+zclcaux(jl)*zlcdqsdt)
           ztt       = (pbetass(jl,jk)-zqvdt+zdqsat-zbap1)/zbbap1
           ztt       = nbetax*MAX(MIN(ztt,1.0_dp),0.0_dp)
           ixidx     = INT(ztt)
           IF (ixidx == nbetax) THEN
              zbetai0 = 1.0_dp
              zbetai1 = 1.0_dp
           ELSE
              zbetai0 = (ztt-ixidx)*tbetai0(iqidx,ixidx+1)             &
                       +(ixidx+1-ztt)*tbetai0(iqidx,ixidx)
              zbetai1 = (ztt-ixidx)*tbetai1(iqidx,ixidx+1)             &
                       +(ixidx+1-ztt)*tbetai1(iqidx,ixidx)
           ENDIF
           zaa        = pbetaa(jl,jk)

!---wiso-code: store temporary values
           zxlb_tmp2(jl) = zxlb(jl)                                            !---wiso-code: store changed value of zxlb
           zxib_tmp2(jl) = zxib(jl)                                            !---wiso-code: store changed value of zxib
           zqp1_tmp2(jl) = zqp1                                                !---wiso-code: store changed value of zqp1   
           ztp1_tmp2(jl) = ztp1                                                !---wiso-code: store changed value of ztp1   
!---wiso-code-end

           zqcdif     = (zbetaqt(jl)-zaa)*(1._dp-zbetai1)              &
                         +(zaa+zqvdt-pbetass(jl,jk)-zdqsat)            &
                         *(1._dp-zbetai0)
           zqcdif_tmp(jl) = MAX(0.0_dp,zqcdif)+zbetacl(jl)                     !---wiso-code: store value of zqcdif
                                                                               !              add cloud water zbetacl to model the
                                                                               !              very low delta values of Antarctica
           zqcdif     = MAX(0.0_dp,zqcdif)-zbetacl(jl)
           lo3_wiso(jl)= (zqcdif.GT.(-zxilb*zclcaux(jl)))                      !---wiso-code: remember water source of zqcdif
           zqcdif     = MAX(zqcdif,-zxilb*zclcaux(jl))                         !              (true=zbetaq-terms, false=zxilb)
           lo4_wiso(jl)= (zqcdif.LT.zqsec*zqp1)                                !---wiso-code: remember water source for zqcdif
           zqcdif     = MIN(zqcdif,zqsec*zqp1)             ! limit to qv       !              (true=zqcdif, false=zqp1)
!
           lo_zqcdif(jl) = (zqcdif.LT.0.0_dp)                                  !---wiso-code: remember if zqcdif is negative
           IF (zqcdif .LT. 0.0_dp) THEN                 ! cloud dissipation
              zifrac   = zxib(jl)/MAX(zepsec,zxilb)
              zifrac   = MAX(MIN(zifrac,1.0_dp),0.0_dp)
              zdep(jl) = zqcdif*zifrac
              zcnd(jl) = zqcdif*(1.0_dp-zifrac)
           ELSE                                      ! cloud generation

              IF (lo2(jl)) THEN                      ! deposition
                 zdep(jl) = zqcdif
                 zcnd(jl) = 0.0_dp
              ELSE                                   ! condensation
                 zcnd(jl) = zqcdif
                 zdep(jl) = 0.0_dp
              END IF
           END IF
!
        ELSE !lcover=.false. or jk < ncctop
           zqcdif         = (zqvdt-zdqsat)*zclcaux(jl)
           zqcdif_tmp(jl) = zqsec*zqp1                                         !---wiso-code: store value of zqcdif
                                                                               !              use total vapor zqp1 here to model the
                                                                               !              very low delta values of Antarctica
           lo3_wiso(jl)   = (zqcdif.GT.(-zxilb*zclcaux(jl)))                   !---wiso-code: remember water source of zqcdif
           zqcdif         = MAX(zqcdif,-zxilb*zclcaux(jl))                     !              (true=zqcdif, false=zxilb)
           lo4_wiso(jl)   = (zqcdif.LT.zqsec*zqp1)                             !---wiso-code: remember water source for zqcdif
           zqcdif         = MIN(zqcdif,zqsec*zqp1)                             !              (true=zqcdif, false=zqp1)

           lo_zqcdif(jl)  = (zqcdif.LT.0.0_dp)                                 !---wiso-code: remember if zqcdif is negative
           IF (zqcdif .LT. 0.0_dp) THEN                 ! cloud dissipation
              zifrac      = zxib(jl)/MAX(zepsec,zxilb)
              zifrac      = MAX(MIN(zifrac,1.0_dp),0.0_dp)
              zdep(jl)    = zqcdif*zifrac
              zcnd(jl)    = zqcdif*(1.0_dp-zifrac)
           ELSE                                      ! cloud generation
              IF (lo2(jl)) THEN                      ! deposition
                 zdep(jl) = zqcdif
                 zcnd(jl) = 0.0_dp
              ELSE                                   ! condensation
                 zcnd(jl) = zqcdif
                 zdep(jl) = 0.0_dp
              END IF
           END IF
        END IF !lcover
!
!       5.4 Accounting for cloud evaporation in clear air and
!           checking for supersaturation
!
        zdep_tmp(jl) = zdep(jl)                            !---wiso-code: store old value of zdep(jl)   
        zcnd_tmp(jl) = zcnd(jl)                            !---wiso-code: store old value of zcnd(jl)   
        zqp1_tmp3(jl) = zqp1                               !---wiso-code: store actual value of zqp1 at begin of section 5.4   

        ztp1tmp(jl) = ztp1+zlvdcp(jl)*zcnd(jl)+zlsdcp(jl)*zdep(jl)
        ztp1tmp_tmp(jl) = ztp1tmp(jl)                      !---wiso-code: store value of ztp1tmp   
        zqp1tmp(jl) = zqp1-zcnd(jl)-zdep(jl)
        zxip1       = MAX(zxised(jl)+zxite(jl)*ztmst-zxievap(jl)           &
                               +zgenti(jl)+zdep(jl),0.0_dp)
        lo2(jl)     = (ztp1tmp(jl) .LT. cthomi) .OR.                       &
                      (ztp1tmp(jl) .LT. tmelt .AND. zxip1 .GT. csecfrl)
        lo2_tmp(jl) = lo2(jl)                              !---wiso-code: store value of lo2(jl)
        it          = NINT(ztp1tmp(jl)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes         = MERGE(tlucua(it),tlucuaw(it),lo2(jl))/papp1(jl,jk)
        zes         = MIN(zes,0.5_dp)
        LO          = zes<0.4_dp
        zcor        = 1._dp/(1._dp-vtmpc1*zes)
        zqsp1tmp    = zes*zcor
        zqsp1tmp_tmp(jl) = zqsp1tmp                        !---wiso-code: store value of zqsp1tmp
        zoversat    = zqsp1tmp*0.01_dp
        zrhtest     = MIN(pqm1(jl,jk)/zqsm1,1._dp)*zqsp1tmp
        it1         = it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1       = MERGE(tlucua(it1),tlucuaw(it1),lo2(jl))/papp1(jl,jk)
        zqst1       = MIN(zqst1,0.5_dp)
        zqst1       = zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt      = (zqst1-zqsp1tmp)*1000._dp
        zlc         = MERGE(zlsdcp(jl),zlvdcp(jl),lo2(jl))
        zlcdqsdt    = MERGE(zlc*zdqsdt,zqsp1tmp*zcor*tlucub(it),LO)
        zqcon       = 1._dp/(1._dp+zlcdqsdt)
!
        IF (lo2(jl)) THEN                                    ! ice cloud
           IF (zqp1tmp(jl) .GT. zqsp1tmp+zoversat) THEN
              zdepcor     = (zqp1tmp(jl)-zqsp1tmp-zoversat)*zqcon
              zdepcor_tmp(jl) = zdepcor                  !---wiso-code: store value of zdepcor 
              zdep(jl)    = zdep(jl)+zdepcor
           END IF
           lo_dep(jl)=(zdep(jl).GT.0.0_dp .AND. zqp1tmp(jl).LT.zrhtest .AND. zqsp1tmp.LE.zqsm1) !---wiso-code: store dep. condition
           IF (zdep(jl) .GT. 0.0_dp .AND. zqp1tmp(jl) .LT. zrhtest     &
                                 .AND. zqsp1tmp .LE. zqsm1) THEN
              zdep(jl)    = zqp1-zrhtest
              zdep(jl)    = MAX(zdep(jl),0.0_dp)
           END IF
        ELSE                                             ! water cloud
           IF (zqp1tmp(jl) .GT. zqsp1tmp+zoversat) THEN
              zcndcor     = (zqp1tmp(jl)-zqsp1tmp-zoversat)*zqcon
              zcndcor_tmp(jl) = zcndcor                  !---wiso-code: store value of zcndcor 
              zcnd(jl)    = zcnd(jl)+zcndcor
           END IF
           lo_cnd(jl)=(zcnd(jl).GT.0.0_dp .AND. zqp1tmp(jl).LT.zrhtest .AND. zqsp1tmp.LE.zqsm1) !---wiso-code: store cond. condition
           IF (zcnd(jl) .GT. 0.0_dp .AND. zqp1tmp(jl) .LT. zrhtest     &
                                 .AND. zqsp1tmp .LE. zqsm1) THEN
              zcnd(jl)    = zqp1-zrhtest
              zcnd(jl)    = MAX(zcnd(jl),0.0_dp)
           END IF
        END IF
!
!       5.5 Change of in-cloud water due to deposition/sublimation and
!           condensation/evaporation (input for cloud microphysics)
!
        lo_depcond(jl)=.FALSE.                             !---wiso-code: store depos/cond condition

        zrelhum=zqp1tmp(jl)/zqsp1tmp
        zdepos =MAX(zdep(jl)+zgenti(jl),0.0_dp)
        zcond  =MAX(zcnd(jl)+zgentl(jl),0.0_dp)
        IF (locc(jl)) THEN
           zxib(jl) = MAX(zxib(jl)+zdep(jl)/zclcaux(jl),0.0_dp)
           zxlb(jl) = MAX(zxlb(jl)+zcnd(jl)/zclcaux(jl),0.0_dp)
        ELSEIF (zdepos>0.0_dp .OR. zcond>0.0_dp) THEN
           lo_depcond(jl)=.TRUE.                           !---wiso-code: store depos/cond condition
           zclcaux(jl)=MAX(MIN(zrelhum,1.0_dp),0.01_dp)
           zxib(jl)   = zdepos/zclcaux(jl)
           zxlb(jl)   = zcond /zclcaux(jl)
        END IF
        ztp1tmp(jl) = ztp1+zlvdcp(jl)*zcnd(jl)+zlsdcp(jl)*zdep(jl)
!
!     ------------------------------------------------------------------
!       6.    Freezing of cloud water
!
!       6.1   Freezing of cloud water for T < 238 K
!
        IF (ztp1tmp(jl) .LE. cthomi) THEN
           zfrl(jl)  = zxlb(jl)*zclcaux(jl)
           zxib(jl)  = zxib(jl)+zxlb(jl)
           zxlb(jl)  = 0.0_dp
        END IF
!
        zfrl_tmp(jl) = zfrl(jl)                            !---wiso-code: store old value of zfrl(jl)   
        zxlb_tmp3(jl) = zxlb(jl)                           !---wiso-code: store old value of zxlb(jl)   

!       6.2   Freezing of cloud water between 238 and 273 K
!
        lo           = zxlb(jl) .GT. 0.0_dp                            &
                                         .AND. ztp1tmp(jl) .LT. tmelt  &
                                         .AND. ztp1tmp(jl) .GT. cthomi
        IF (lo) THEN
           zfrl(jl)  = 100._dp*(EXP(0.66_dp*(tmelt-ztp1tmp(jl)))       &
                          -1._dp)*zrho(jl,jk)/(rhoh2o*pacdnc(jl,jk))
           zfrl(jl)  = zxlb(jl)*(1._dp-1._dp/                          &
                                      (1._dp+zfrl(jl)*ztmst*zxlb(jl)))
           zradl     = (0.75_dp*zxlb(jl)*zrho(jl,jk)                   &
                            /(api*rhoh2o*pacdnc(jl,jk)))**(1._dp/3._dp)
           zf1       = 4._dp*api*zradl*pacdnc(jl,jk)*2.e5_dp           &
                                 *(tmelt-3._dp-ztp1tmp(jl))/zrho(jl,jk)
           zf1       = MAX(0.0_dp,zf1)
           zfrl(jl)  = zfrl(jl)+ztmst*1.4e-20_dp*zf1
           zfrl(jl)  = MAX(0.0_dp,MIN(zfrl(jl),zxlb(jl)))
           zfrl_tmp2(jl) = zfrl(jl)                          !---wiso-code: store value of zfrl(jl)   
           zxlb(jl)  = zxlb(jl)-zfrl(jl)
           zxib(jl)  = zxib(jl)+zfrl(jl)
           zfrl(jl)  = zfrl(jl)*zclcaux(jl)
        END IF
!
610  END DO


!---wiso-code
  IF (lwiso) THEN

!       Caculate fractionation factors for condensation/deposition and evaporation/sublimation
        CALL wiso_frac_liq_ice(kproma,kbdim,kwiso,ztp1_tmp(:),zwisofracliq,zwisofracice)

        DO jt = 1,kwiso
          DO jl = 1,kproma

!           ------------------------------------------------------------------
!           Sedimentation of cloud ice from grid-mean values and
!           Updating the tendency 'pxite' to include sedimentation - water isotopes

!           zwisoxised is determined as a mixture from zwisoxiflux and zwisoxip1

            zwisoxip1         = pwisoxim1(jl,jk,jt)+pwisoxite(jl,jk,jt)*ztmst-zwisoimlt(jl,jt)
            zwisoxip1         = MAX(zwisoxip1,EPSILON(1._dp)*tnat(jt))
            zxifall           = cvtfall*(zrho(jl,jk)*zxip1_tmp(jl))**0.16_dp
            zal1              = zxifall*g*zrho(jl,jk)*ztmst/zdp(jl)
            zal2              = zwisoxiflux(jl,jt)/(zrho(jl,jk)*zxifall)
            zwisoxised        = zwisoxip1*EXP(-zal1)+zal2*(1._dp-EXP(-zal1))
            zwisoxiflux(jl,jt)= zwisoxiflux(jl,jt)+(zwisoxip1-zwisoxised)*zcons2*zdp(jl)
            pwisoxite(jl,jk,jt) = (zwisoxised-pwisoxim1(jl,jk,jt))/ztmst

!           In-cloud water/ice calculated from respective grid-means,
!           partial cloud cover, advective/diffusive tendencies,
!           detrained cloud water/ice and ice sedimentation - water isotopes

            lo2_new(jl)       = (ptm1(jl,jk) .LT. cthomi) .OR.                   &
                              (ptm1(jl,jk) .LT. tmelt .AND. zxised(jl) .GT. csecfrl)
                          
            IF (lo2_new(jl)) THEN                                         ! ice cloud
              zwisoxite(jl,jt)      = pwisoxtec(jl,jk,jt)
              zwisoxlte(jl,jt)      = 0.0_dp
              IF (locc(jl)) THEN                                   ! cloud cover > 0.
                zwisoxib(jl,jt)     = pwisoxim1(jl,jk,jt)/zclcaux_tmp(jl)
                zwisoxlb(jl,jt)     = pwisoxlm1(jl,jk,jt)/zclcaux_tmp(jl)
                zwisoxim1evp(jl,jt) = 0.0_dp
                zwisoxlm1evp(jl,jt) = 0.0_dp
                zwisoxidt           = (pwisoxite(jl,jk,jt)+zwisoxite(jl,jt))*ztmst
                zwisoxldt           =  pwisoxlte(jl,jk,jt)*ztmst+zwisoximlt(jl,jt)+zwisoimlt(jl,jt)
                IF (lo_zxidt(jl)) THEN
                  zwisoxidtstar(jl,jt) = zwisoxidt
                  zwisoxib(jl,jt)   = zwisoxib(jl,jt)+zwisoxidt
                ELSE
                  zwisoxidtstar(jl,jt) = 0.0_dp
                  zwisoxib(jl,jt)   = zwisoxib(jl,jt)+MAX(zwisoxidt/zclcaux_tmp(jl),           &
                                      -zwisoxib(jl,jt))
                  pwisoxite(jl,jk,jt) = MAX(pwisoxite(jl,jk,jt),                           &
                                       -(pwisoxim1(jl,jk,jt)/ztmst+zwisoxite(jl,jt)))
                END IF
                IF (lo_zxldt(jl)) THEN
                  zwisoxldtstar(jl,jt) = zwisoxldt
                  zwisoxlb(jl,jt)   = zwisoxlb(jl,jt)+zwisoxldt
                ELSE
                  zwisoxldtstar(jl,jt) = 0.0_dp
                  zwisoxlb(jl,jt)   = zwisoxlb(jl,jt)+MAX(zwisoxldt/zclcaux_tmp(jl),           &
                                      -zwisoxlb(jl,jt))
                  pwisoxlte(jl,jk,jt) = MAX(pwisoxlte(jl,jk,jt),-pwisoxlm1(jl,jk,jt)/ztmst)
                END IF
              ELSE                                             ! cloud cover = 0.
                zwisoxib(jl,jt)     = 0.0_dp
                zwisoxlb(jl,jt)     = 0.0_dp
                zwisoxidt           = (pwisoxite(jl,jk,jt)+zwisoxite(jl,jt))*ztmst
                zwisoxldt           =  pwisoxlte(jl,jk,jt)*ztmst+zwisoximlt(jl,jt)+zwisoimlt(jl,jt)
                IF (lo_zxidt(jl)) THEN
                  zwisoxidtstar(jl,jt) = zwisoxidt
                  zwisoxim1evp(jl,jt) = pwisoxim1(jl,jk,jt)
                ELSE
                  zwisoxidtstar(jl,jt) = 0.0_dp
                  pwisoxite(jl,jk,jt) = MAX(pwisoxite(jl,jk,jt),                           &
                                       -(pwisoxim1(jl,jk,jt)/ztmst+zwisoxite(jl,jt)))
                  zwisoxim1evp(jl,jt)= pwisoxim1(jl,jk,jt)+(pwisoxite(jl,jk,jt)+zwisoxite(jl,jt))*ztmst
                END IF
                IF (lo_zxldt(jl)) THEN
                  zwisoxldtstar(jl,jt) = zwisoxldt
                  zwisoxlm1evp(jl,jt)= pwisoxlm1(jl,jk,jt)
                ELSE
                  zwisoxldtstar(jl,jt) = 0.0_dp
                  pwisoxlte(jl,jk,jt) = MAX(pwisoxlte(jl,jk,jt),-pwisoxlm1(jl,jk,jt)/ztmst)
                  zwisoxlm1evp(jl,jt) = pwisoxlm1(jl,jk,jt)+pwisoxlte(jl,jk,jt)*ztmst
                END IF
              END IF
            ELSE                                               ! water cloud
              zwisoxlte(jl,jt)      = pwisoxtec(jl,jk,jt)
              zwisoxite(jl,jt)      = 0.0_dp
              IF (locc(jl)) THEN                                   ! cloud cover > 0.
                zwisoxlb(jl,jt)     = pwisoxlm1(jl,jk,jt)/zclcaux_tmp(jl)
                zwisoxib(jl,jt)     = pwisoxim1(jl,jk,jt)/zclcaux_tmp(jl)
                zwisoxlm1evp(jl,jt) = 0.0_dp
                zwisoxim1evp(jl,jt) = 0.0_dp
                zwisoxldt           = (pwisoxlte(jl,jk,jt)+zwisoxlte(jl,jt))*ztmst         &
                                      +zwisoximlt(jl,jt)+zwisoimlt(jl,jt)
                zwisoxidt           = pwisoxite(jl,jk,jt)*ztmst
                IF (lo_zxldt(jl)) THEN
                  zwisoxldtstar(jl,jt) = zwisoxldt
                  zwisoxlb(jl,jt)     = zwisoxlb(jl,jt)+zwisoxldt
                ELSE
                  zwisoxldtstar(jl,jt) = 0.0_dp
                  zwisoxlb(jl,jt)   = zwisoxlb(jl,jt)+MAX(zwisoxldt/zclcaux_tmp(jl),           &
                                      -zwisoxlb(jl,jt))
                  pwisoxlte(jl,jk,jt) = MAX(pwisoxlte(jl,jk,jt),                           &
                                       -(pwisoxlm1(jl,jk,jt)/ztmst+zwisoxlte(jl,jt)))
                END IF
                IF (lo_zxidt(jl)) THEN
                  zwisoxidtstar(jl,jt) = zwisoxidt
                  zwisoxib(jl,jt)   = zwisoxib(jl,jt)+zwisoxidt
                ELSE
                  zwisoxidtstar(jl,jt) = 0.0_dp
                  zwisoxib(jl,jt)   = zwisoxib(jl,jt)+MAX(zwisoxidt/zclcaux_tmp(jl),           &
                                      -zwisoxib(jl,jt))
                  pwisoxite(jl,jk,jt) = MAX(pwisoxite(jl,jk,jt),-pwisoxim1(jl,jk,jt)/ztmst)
                END IF
              ELSE                                             ! cloud cover = 0.
                zwisoxlb(jl,jt)     = 0.0_dp
                zwisoxib(jl,jt)     = 0.0_dp
                zwisoxldt           = (pwisoxlte(jl,jk,jt)+zwisoxlte(jl,jt))*ztmst         &
                                      +zwisoximlt(jl,jt)+zwisoimlt(jl,jt)
                zwisoxidt           = pwisoxite(jl,jk,jt)*ztmst
                IF (lo_zxldt(jl)) THEN
                  zwisoxldtstar(jl,jt) = zwisoxldt
                  zwisoxlm1evp(jl,jt)= pwisoxlm1(jl,jk,jt)
                ELSE
                  zwisoxldtstar(jl,jt) = 0.0_dp
                  pwisoxlte(jl,jk,jt) = MAX(pwisoxlte(jl,jk,jt),                           &
                                        -(pwisoxlm1(jl,jk,jt)/ztmst+zwisoxlte(jl,jt)))
                  zwisoxlm1evp(jl,jt) = pwisoxlm1(jl,jk,jt)+(pwisoxlte(jl,jk,jt)+zwisoxlte(jl,jt))*ztmst
                END IF
                IF (lo_zxidt(jl)) THEN
                  zwisoxidtstar(jl,jt) = zwisoxidt
                  zwisoxim1evp(jl,jt)= pwisoxim1(jl,jk,jt)
                ELSE
                  zwisoxidtstar(jl,jt) = 0.0_dp
                  pwisoxite(jl,jk,jt) = MAX(pwisoxite(jl,jk,jt),-pwisoxim1(jl,jk,jt)/ztmst)
                  zwisoxim1evp(jl,jt)= pwisoxim1(jl,jk,jt)+pwisoxite(jl,jk,jt)*ztmst
                END IF
              END IF
            END IF
            
!           ------------------------------------------------------------------
!           Condensation/deposition and evaporation/sublimation - water isotopes
!
            zwisoxievap(jl,jt) = (1.0_dp-zclcaux_tmp(jl))*zwisoxidtstar(jl,jt)+zwisoxim1evp(jl,jt)
            zwisoxlevap(jl,jt) = (1.0_dp-zclcaux_tmp(jl))*zwisoxldtstar(jl,jt)+zwisoxlm1evp(jl,jt)
            zwisoqvdt          = pwisoqte(jl,jk,jt)*ztmst+zwisoevp(jl,jt)+zwisosub(jl,jt)          &
                                +zwisoxievap(jl,jt)+zwisoxlevap(jl,jt)+zwisoxisub(jl,jt)
            zwisoqp1           = MAX(pwisoqm1(jl,jk,jt)+zwisoqvdt,0.0_dp)
            zwisoxib(jl,jt)    = MAX(zwisoxib(jl,jt),0.0_dp)
            zwisoxlb(jl,jt)    = MAX(zwisoxlb(jl,jt),0.0_dp)
            zwisoxilb(jl,jt)   = zwisoxib(jl,jt)+zwisoxlb(jl,jt)

            IF (lcover .AND. jk.GE.ncctop) THEN

!           translated into apparent xl,xi,q and heat sources
!           first order effect only, effect of evaporation of
!           cloud on qsat taken into account in thermodynamic budget
!           but does not change the mixing term here since that
!           would require iteration and is therefore neglected - water isotopes
!
!           calculate values after one timestep

              zwisodenom=zxlb_tmp(jl)+zxib_tmp(jl)                                       !source: zxlb+zxib
              zwisonumer=zwisoxlb(jl,jt)+zwisoxib(jl,jt)
              zwisodenom=MERGE(pqm1(jl,jk),zwisodenom,lo_wiso(jl))                       !source: pqm1
              zwisonumer=MERGE(pwisoqm1(jl,jk,jt),zwisonumer,lo_wiso(jl))
              zwisodenom=MERGE(zqp1_tmp(jl),zwisodenom,lo2_wiso(jl))                     !source: zqp1
              zwisonumer=MERGE(zwisoqp1,zwisonumer,lo2_wiso(jl))

              IF (zgenti(jl).GT.0._dp) THEN ! deposition of vapour to frozen cloud ice (fractionation as open system)
                IF (ABS(zwisodenom).GT.cwisomin) THEN
                  zdelta = 1._dp - (zgenti(jl)/zwisodenom)
                  zwisogenti(jl,jt) = zwisonumer*(1._dp-zdelta**zwisofracice(jl,jt))
                ELSE
                  zdelta=tnat(jt)
                  zwisogenti(jl,jt) = zgenti(jl)*zdelta
                ENDIF
              ELSE                          ! sublimation of ice in cloud (no fractionation due to low diffusivity of isotopes)
                zdelta=tnat(jt)
                IF (ABS(zwisodenom).GT.cwisomin) zdelta = zwisonumer/zwisodenom
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisogenti(jl,jt) = zgenti(jl)*zdelta
              ENDIF

              IF (zgentl(jl).GT.0._dp) THEN ! condensation of vapour to liquid cloud water (fractionation as closed system)
                zwisodenom=zwisodenom+(zwisofracliq(jl,jt)-1._dp)*zgentl(jl)
                zdelta=tnat(jt)
                IF (ABS(zwisodenom).GT.cwisomin) zdelta = zwisofracliq(jl,jt)*zwisonumer/zwisodenom
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisogentl(jl,jt) = zgentl(jl)*zdelta
              ELSE                          ! evaporation of liquid water (fractionation as closed system) 
                zwisodenom=zwisofracliq(jl,jt)*zwisodenom+(1._dp-zwisofracliq(jl,jt))*zgentl(jl)
                zdelta=tnat(jt)
                IF (ABS(zwisodenom).GT.cwisomin) zdelta = zwisonumer/zwisodenom
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisogentl(jl,jt) = zgentl(jl)*zdelta
              ENDIF

              IF (locc(jl)) THEN
                zwisoxib(jl,jt) = MAX(zwisoxib(jl,jt)+zwisogenti(jl,jt)/zclcaux_tmp(jl),0.0_dp)
                zwisoxlb(jl,jt) = MAX(zwisoxlb(jl,jt)+zwisogentl(jl,jt)/zclcaux_tmp(jl),0.0_dp)
              END IF

!           Deposition/sublimation of cloud ice and condensation/
!           evaporation of liquid water due to changes in water vapour
!           and temperature (advection, convection, turbulent mixing,
!           evaporation of rain, sublimation and melting of snow).
!
              zwisoqvdt = zwisoqvdt-(zwisogenti(jl,jt)+zwisogentl(jl,jt))      ! here: zwisogent = zwisogenti+zwisogentl 
              zwisoqp1  = MAX(pwisoqm1(jl,jk,jt)+zwisoqvdt,0.0_dp)

!             Calculate new fractionation coefficients depending on new temperature value ztp1_tmp2
!             (explicit in-line calculation to enable code vectorization)
              
              zwisofracliq(jl,jt) =    &
                (exp(talphal1(jt)/(ztp1_tmp2(jl)**2._dp)+talphal2(jt)/ztp1_tmp2(jl)+talphal3(jt)))**talphal4(jt) ! frac. liquid
              zwisofracice(jl,jt) =    &
                (exp(talphas1(jt)/(ztp1_tmp2(jl)**2._dp)+talphas2(jt)/ztp1_tmp2(jl)+talphas3(jt)))**talphas4(jt) ! frac. ice

              IF (nwisotyp(jt).ne.1.and.ztp1_tmp2(jl).lt.tmelt) THEN        ! effective fractionation over ice if necessary
                 zsatval=tsatbase-tsatfac*(ztp1_tmp2(jl)-tmelt)
                 zwisofracice(jl,jt) = zwisofracice(jl,jt)*(zsatval/(1._dp+zwisofracice(jl,jt)*(zsatval-1._dp)*tdifrel(jt)))
              ENDIF

              ! calculate isotope equivalent of zqcdif (assume same isotope ratio as for zqp1)
              zdelta=tnat(jt)
              IF (zqp1_tmp2(jl).GT.cwisomin) zdelta = zwisoqp1/zqp1_tmp2(jl)
              IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
              zwisoqcdif = zqcdif_tmp(jl) * zdelta

              zwisodenom=zxlb_tmp2(jl)+zxib_tmp2(jl)                         !source: zxlb+zxib
              zwisonumer=zwisoxlb(jl,jt)+zwisoxib(jl,jt)
              zwisodenom=MERGE(zqcdif_tmp(jl),zwisodenom,lo3_wiso(jl))       !source: zqcdif
              zwisonumer=MERGE(zwisoqcdif,zwisonumer,lo3_wiso(jl))
              zwisodenom=MERGE(zwisodenom,zqp1_tmp2(jl),lo4_wiso(jl))        !source: zqp1
              zwisonumer=MERGE(zwisonumer,zwisoqp1,lo4_wiso(jl))
                
              IF (lo_zqcdif(jl)) THEN                                          ! cloud dissipation
              
                ! sublimation of ice in cloud (no fractionation due to low diffusivity of isotopes in ice)
                zdelta=tnat(jt)
                IF (ABS(zwisodenom).GT.cwisomin) zdelta = zwisonumer/zwisodenom
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisodep(jl,jt) = zdep_tmp(jl)*zdelta

                ! evaporation of liquid water (fractionation as closed system)                  
                zwisodenom=zwisofracliq(jl,jt)*zwisodenom+(1._dp-zwisofracliq(jl,jt))*zcnd_tmp(jl)
                zdelta=tnat(jt)
                IF (ABS(zwisodenom).GT.cwisomin) zdelta = zwisonumer/zwisodenom
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisocnd(jl,jt) = zcnd_tmp(jl)*zdelta

              ELSE                                                             ! cloud generation
              
                IF (lo2_new(jl)) THEN            ! deposition to frozen cloud ice (fractionation as open system)

                  IF (ABS(zwisodenom).GT.cwisomin) THEN
                    zdelta = 1._dp - (zdep_tmp(jl)/zwisodenom)
                    zwisodep(jl,jt) = zwisonumer*(1._dp-zdelta**zwisofracice(jl,jt))
                  ELSE
                    zdelta=tnat(jt)
                    zwisodep(jl,jt) = zdep_tmp(jl)*zdelta
                  ENDIF
                  zwisocnd(jl,jt) = 0.0_dp

                ELSE                             ! condensation to liquid cloud water (fractionation as closed system)
                
                  zwisodenom=zwisodenom+(zwisofracliq(jl,jt)-1._dp)*zcnd_tmp(jl)
                  zdelta=tnat(jt)
                  IF (ABS(zwisodenom).GT.cwisomin) zdelta = zwisofracliq(jl,jt)*zwisonumer/zwisodenom
                  IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                  zwisocnd(jl,jt) = zcnd_tmp(jl)*zdelta
                  zwisodep(jl,jt) = 0.0_dp

                END IF

              END IF

            ELSE                                                              !lcover=.false. or jk < ncctop
            
              ! calculate isotope equivalent of zqcdif (assume same isotope ratio as for zqp1)
              zdelta=tnat(jt)
              IF (zqp1_tmp(jl).GT.cwisomin) zdelta = zwisoqp1/zqp1_tmp(jl)
              IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
              zwisoqcdif = zqcdif_tmp(jl) * zdelta

              zwisodenom=zxlb_tmp(jl)+zxib_tmp(jl)                           !source: zxlb+zxib
              zwisonumer=zwisoxlb(jl,jt)+zwisoxib(jl,jt)
              zwisodenom=MERGE(zqcdif_tmp(jl),zwisodenom,lo3_wiso(jl))       !source: zqcdif
              zwisonumer=MERGE(zwisoqcdif,zwisonumer,lo3_wiso(jl))
              zwisodenom=MERGE(zwisodenom,zqp1_tmp(jl),lo4_wiso(jl))         !source: zqp1
              zwisonumer=MERGE(zwisonumer,zwisoqp1,lo4_wiso(jl))
                
              IF (lo_zqcdif(jl)) THEN                                          ! cloud dissipation (zqcdif .LT. 0.)

                ! sublimation of ice in cloud (no fractionation due to low diffusivity of isotopes in ice)
                zdelta=tnat(jt)
                IF (ABS(zwisodenom).GT.cwisomin) zdelta = zwisonumer/zwisodenom
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisodep(jl,jt) = zdep_tmp(jl)*zdelta

                ! evaporation of liquid water (fractionation as closed system)                  
                zwisodenom=zwisofracliq(jl,jt)*zwisodenom+(1._dp-zwisofracliq(jl,jt))*zcnd_tmp(jl)
                zdelta=tnat(jt)
                IF (ABS(zwisodenom).GT.cwisomin) zdelta = zwisonumer/zwisodenom
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisocnd(jl,jt) = zcnd_tmp(jl)*zdelta

              ELSE                                                             ! cloud generation
              
                IF (lo2_new(jl)) THEN            ! deposition to frozen cloud ice (fractionation as open system)

                  IF (ABS(zwisodenom).GT.cwisomin) THEN
                    zdelta = 1._dp - (zdep_tmp(jl)/zwisodenom)
                    zwisodep(jl,jt) = zwisonumer*(1._dp-zdelta**zwisofracice(jl,jt))
                  ELSE
                    zdelta=tnat(jt)
                    zwisodep(jl,jt) = zdep_tmp(jl)*zdelta
                  ENDIF
                  zwisocnd(jl,jt) = 0.0_dp

                ELSE                             ! condensation to liquid cloud water (fractionation as closed system)
                                   
                  zwisodenom=zwisodenom+(zwisofracliq(jl,jt)-1._dp)*zcnd_tmp(jl)
                  zdelta=tnat(jt)
                  IF (ABS(zwisodenom).GT.cwisomin) zdelta = zwisofracliq(jl,jt)*zwisonumer/zwisodenom
                  IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                  zwisocnd(jl,jt) = zcnd_tmp(jl)*zdelta
                  zwisodep(jl,jt) = 0.0_dp

                END IF
                
              END IF
            
            ENDIF                                              !lcover

!           Accounting for cloud evaporation in clear air and
!           checking for supersaturation                       - water isotopes

!           Calculate new fractionation coefficients depending on new temperature value ztp1tmp_tmp
!           (explicit in-line calculation to enable code vectorization)
              
            zwisofracliq(jl,jt) =      &
              (exp(talphal1(jt)/(ztp1tmp_tmp(jl)**2._dp)+talphal2(jt)/ztp1tmp_tmp(jl)+talphal3(jt)))**talphal4(jt) ! frac. liquid
            zwisofracice(jl,jt) =      &
              (exp(talphas1(jt)/(ztp1tmp_tmp(jl)**2._dp)+talphas2(jt)/ztp1tmp_tmp(jl)+talphas3(jt)))**talphas4(jt) ! frac. ice

            IF (nwisotyp(jt).ne.1.and.ztp1tmp_tmp(jl).lt.tmelt) THEN        ! effective fractionation over ice if necessary
              zsatval=tsatbase-tsatfac*(ztp1tmp_tmp(jl)-tmelt)
              zwisofracice(jl,jt) = zwisofracice(jl,jt)*(zsatval/(1._dp+zwisofracice(jl,jt)*(zsatval-1._dp)*tdifrel(jt)))
            ENDIF

            zwisoqp1tmp(jl,jt) = zwisoqp1-zwisocnd(jl,jt)-zwisodep(jl,jt)

            zoversat = zqsp1tmp_tmp(jl)*0.01_dp
           
            IF (lo2_tmp(jl)) THEN                                              ! ice cloud
              IF (zqp1tmp(jl) .GT. zqsp1tmp_tmp(jl)+zoversat) THEN
              ! deposition of vapour to frozen cloud ice (fractionation as open system), source: zqp1tmp
                zwisodenom=zqp1tmp(jl)
                IF (ABS(zwisodenom).GT.cwisomin) THEN
                  zdelta = 1._dp-(zdepcor_tmp(jl)/zwisodenom)
                  zwisodepcor = zwisoqp1tmp(jl,jt)*(1._dp-zdelta**zwisofracice(jl,jt))
                ELSE
                  zdelta=tnat(jt)
                  zwisodepcor = zdepcor_tmp(jl)*zdelta
                ENDIF
                zwisodep(jl,jt)= zwisodep(jl,jt)+zwisodepcor
              END IF
              IF (lo_dep(jl)) THEN
              ! deposition of vapour to frozen cloud ice (fractionation as open system), source: zqp1
                zwisodenom=zqp1_tmp3(jl)
                IF (ABS(zwisodenom).GT.cwisomin) THEN
                  zdelta = 1._dp-(zdep(jl)/zwisodenom)
                  zwisodep(jl,jt) = zwisoqp1*(1._dp-zdelta**zwisofracice(jl,jt))
                ELSE
                  zdelta=tnat(jt)
                  zwisodep(jl,jt) = zdep(jl)*zdelta
                ENDIF
              ENDIF
            ELSE                                                               ! water cloud
              IF (zqp1tmp(jl) .GT. zqsp1tmp_tmp(jl)+zoversat) THEN
              ! condensation of vapour to liquid cloud water (fractionation as closed system), source: zqp1tmp
                zwisodenom=zqp1tmp(jl)+(zwisofracliq(jl,jt)-1._dp)*zcndcor_tmp(jl)
                zdelta=tnat(jt)
                IF (ABS(zwisodenom).GT.cwisomin) zdelta = zwisofracliq(jl,jt)*zwisoqp1tmp(jl,jt)/zwisodenom
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisocndcor = zcndcor_tmp(jl)*zdelta
                zwisocnd(jl,jt)= zwisocnd(jl,jt)+zwisocndcor
              END IF
              IF (lo_cnd(jl)) THEN
              ! condensation of vapour to liquid cloud water (fractionation as closed system), source: zqp1
                zwisodenom=zqp1_tmp3(jl)+(zwisofracliq(jl,jt)-1._dp)*zcnd(jl)
                zdelta=tnat(jt)
                IF (ABS(zwisodenom).GT.cwisomin) zdelta = zwisofracliq(jl,jt)*zwisoqp1/zwisodenom
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisocnd(jl,jt) = zcnd(jl)*zdelta
              ENDIF
            ENDIF

!           Change of in-cloud water due to deposition/sublimation and
!           condensation/evaporation (input for cloud microphysics) - water isotopes

            zwisodepos =MAX(zwisodep(jl,jt)+zwisogenti(jl,jt),0.0_dp)
            zwisocond  =MAX(zwisocnd(jl,jt)+zwisogentl(jl,jt),0.0_dp)
            IF (locc(jl)) THEN
              zwisoxib(jl,jt) = MAX(zwisoxib(jl,jt)+zwisodep(jl,jt)/zclcaux_tmp(jl),0.0_dp)
              zwisoxlb(jl,jt) = MAX(zwisoxlb(jl,jt)+zwisocnd(jl,jt)/zclcaux_tmp(jl),0.0_dp)
            ELSEIF (lo_depcond(jl)) THEN              ! here: use new value of zclcaux(jl)
              zwisoxib(jl,jt) = zwisodepos/zclcaux(jl)
              zwisoxlb(jl,jt) = zwisocond /zclcaux(jl)
            END IF

!           ------------------------------------------------------------------
!           Freezing of cloud water
!           - assume no fractionation for phase change between liquid and frozen cloud water
!             (this assumption might not be correct in reality)
!
!           Freezing of cloud water for T < 238 K
!
            IF (ztp1tmp(jl) .LE. cthomi) THEN
              zwisofrl(jl,jt) = zwisoxlb(jl,jt)*zclcaux(jl)
              zwisoxib(jl,jt) = zwisoxib(jl,jt)+zwisoxlb(jl,jt)
              zwisoxlb(jl,jt) = 0.0_dp
            END IF

!           Freezing of cloud water between 238 and 273 K
!
            lo                = zxlb_tmp3(jl) .GT. 0.0_dp                                  &
                                         .AND. ztp1tmp(jl) .LT. tmelt                      &
                                         .AND. ztp1tmp(jl) .GT. cthomi
            IF (lo) THEN
              IF (zfrl_tmp2(jl).NE.zfrl_tmp(jl)) THEN
                zdelta=tnat(jt)
                IF (zxlb_tmp3(jl).GT.cwisomin) zdelta = zwisoxlb(jl,jt)/zxlb_tmp3(jl)
                IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
                zwisofrl(jl,jt) = zfrl_tmp2(jl)*zdelta
                zwisoxlb(jl,jt) = zwisoxlb(jl,jt)-zwisofrl(jl,jt)
                zwisoxib(jl,jt) = zwisoxib(jl,jt)+zwisofrl(jl,jt)
                zwisofrl(jl,jt) = zwisofrl(jl,jt)*zclcaux(jl)
              END IF
            END IF

          END DO ! end of jl=1,kproma loop
        END DO ! end of jt=1,kwiso loop

  END IF
!---wiso-code-end

!
     IF (lookupoverflow) CALL lookuperror ('cloud (2)    ')
!
!     ------------------------------------------------------------------
!       7.  Cloud physics and precipitation fluxes at the surface
!
     DO 701 jl = 1,kproma
        locc(jl) = zclcaux(jl) .GT. 0.0_dp
        zclcstar = MIN(zclcaux(jl),zclcpre(jl))
        zauloc   = cauloc*zdz(jl)/5000._dp
        zauloc   = MAX(MIN(zauloc,clmax),clmin)
!
        jb=knvb(jl)
        lo=(jb.GE.jbmin .AND. jb.LE.jbmax .AND. pvervel(jl,jk).GT.0._dp)
        lo1=(jk.EQ.jb .OR. jk.EQ.jb+1)
        IF(lo .AND. lo1 .AND. lonacc) zauloc= 0.0_dp
!
        zqrho    = 1.3_dp/zrho(jl,jk)
        zxlb(jl) = MAX(zxlb(jl),1.e-20_dp)
        zxib(jl) = MAX(zxib(jl),1.e-20_dp)
        IF (zclcpre(jl) .GT. 0.0_dp) THEN
           zxrp1 = (zrfl(jl)/(zclcpre(jl)*12.45_dp*SQRT(zqrho)))       &
                                                        **(8._dp/9._dp)
           zxsp1 = (zsfl(jl)/(zclcpre(jl)*cvtfall))**(1._dp/1.16_dp)
        ELSE
           zxrp1 = 0.0_dp
           zxsp1 = 0.0_dp
        END IF

        zxlb_tmp(jl) = zxlb(jl)                          !---wiso-code: store old value of zxlb(jl)   
        zxib_tmp(jl) = zxib(jl)                          !---wiso-code: store old value of zxib(jl)   
        zrpr_tmp(jl) = zrpr(jl)                          !---wiso-code: store old value of zrpr(jl)   
        zspr_tmp(jl) = zspr(jl)                          !---wiso-code: store old value of zspr(jl)   
!
!       7.1   Warm clouds: Coalescence processes after Beheng (1994):
!             Autoconversion of cloud droplets and collection of cloud
!             droplets by falling rain. Accretion of cloud droplets by
!             falling snow (zsacl) is calculated under 7.2
!
        IF (locc(jl) .AND. (zxlb(jl) > cqtmin .OR. zxib(jl) > cqtmin)) THEN
           zraut    = ccraut*1.2e27_dp/zrho(jl,jk)                     &
                                            *(pacdnc(jl,jk)*1.e-6_dp)  &
                             **(-3.3_dp)*(zrho(jl,jk)*1.e-3_dp)**4.7_dp
           zexm1    = 4.7_dp-1.0_dp
           zexp     = -1._dp/zexm1
           zraut    = zxlb(jl)*(1._dp-(1._dp+zraut*ztmst*zexm1*zxlb(jl)&
                                                       **zexm1)**zexp)
           zxlb(jl) = zxlb(jl)-zraut
           zrac1    = 6._dp*zxrp1*ztmst
           zrac1    = zxlb(jl)*(1._dp-EXP(-zrac1))
           zxlb(jl) = zxlb(jl)-zrac1
           zrac2    = 6._dp*zauloc*zrho(jl,jk)*zraut*ztmst
           zrac2    = zxlb(jl)*(1._dp-EXP(-zrac2))
           zxlb(jl) = zxlb(jl)-zrac2
           zxlb_difw(jl) = zraut+zrac1+zrac2                         !---wiso-code: store sum of changes to zxlb for warm clouds 
           zrpr(jl) = zrpr(jl)+zclcaux(jl)*(zraut+zrac2)+zclcstar*zrac1
!
!       7.2  Cold clouds:
!            Conversion of cloud ice to snow after Levkov et al. 1992:
!            Aggregation of ice crystals to snow and accretion of ice
!            by falling snow.
!            Accrection of cloud droplets by falling snow.
!            Effective radius of ice crystals after Moss (1995)
!
           zrieff    = 83.8_dp*(zxib(jl)*zrho(jl,jk)*1000._dp)**0.216_dp
           zrieff    = MIN(MAX(zrieff,ceffmin),ceffmax)
           zrih      = -2261._dp                                       &
                       +SQRT(5113188._dp+2809._dp*zrieff*zrieff*zrieff)
           zri       = 1.e-6_dp*zrih**(1._dp/3._dp)
           zcolleffi = EXP(0.025_dp*(ztp1tmp(jl)-tmelt))
           zc1       = 17.5_dp*zrho(jl,jk)/crhoi*zqrho**0.33_dp
           zdt2      = -6._dp/zc1*LOG10(zri*1.e4_dp)
           zsaut     = ccsaut/zdt2
           zsaut     = zxib(jl)*(1._dp-1._dp/                          &
                                          (1._dp+zsaut*ztmst*zxib(jl)))
           zxib(jl)  = zxib(jl)-zsaut
           zsaci1    = 0.0_dp
           zsaci2    = 0.0_dp
           zsacl1    = 0.0_dp
           zsacl2    = 0.0_dp
           IF (zxsp1 .GT. cqtmin) THEN
              zlamsm    = (zxsp1/(api*crhosno*cn0s))**0.8125_dp
              zsaci1    = api*cn0s*3.078_dp*zlamsm*zqrho**0.5_dp
              zsacl1    = zxlb(jl)*(1._dp-EXP(-zsaci1*ccsacl*ztmst))
              zxlb(jl)  = zxlb(jl)-zsacl1
              zsacl1    = zclcstar*zsacl1
              zsaci1    = zsaci1*zcolleffi*ztmst
              zsaci1    = zxib(jl)*(1._dp-EXP(-zsaci1))
              zxib(jl)  = zxib(jl)-zsaci1
           END IF
           zxsp2        = zauloc*zrho(jl,jk)*zsaut
           IF (zxsp2 .GT. cqtmin) THEN
              zlamsm    = (zxsp2/(api*crhosno*cn0s))**0.8125_dp
              zsaci2    = api*cn0s*3.078_dp*zlamsm*zqrho**0.5_dp
              zsacl2    = zxlb(jl)*(1._dp-EXP(-zsaci2*ccsacl*ztmst))
              zxlb(jl)  = zxlb(jl)-zsacl2
              zsacl2    = zclcaux(jl)*zsacl2
              zsaci2    = zsaci2*zcolleffi*ztmst
              zsaci2    = zxib(jl)*(1._dp-EXP(-zsaci2))
              zxib(jl)  = zxib(jl)-zsaci2
           END IF
           zxib_difc(jl) = zsaut+zsaci1+zsaci2                         !---wiso-code: store sum of changes to zxib for cold clouds 
           zsacl(jl)    = zsacl1+zsacl2
           zspr(jl)     = zspr(jl)+zclcaux(jl)*(zsaut+zsaci2)          &
                                  +zclcstar*zsaci1
        END IF
!
!       7.3 Updating precipitation fluxes. In the lowest layer (klev),
!           the sedimentation sink of cloud ice is balanced
!           by precipitation at the surface (through 'zzdrs').
!           Fraction of precipitating clouds (zclcpre) used for the
!           calculation of evaporation/sublimation of rain/snow in
!           the next layer
!
        zzdrr          = zcons2*zdp(jl)*zrpr(jl)
        zzdrs          = zcons2*zdp(jl)*(zspr(jl)+zsacl(jl))
        IF (jk .EQ. klev) THEN
           zzdrs       = zzdrs+zxiflux(jl)
           zcons       = zcons2*zdp(jl)/(zlsdcp(jl)-zlvdcp(jl))
           zsnmlt      = MIN(zxsec*zzdrs,zcons                         &
                                       *MAX(0._dp,(ztp1tmp(jl)-tmelt)))
           zzdrs_tmp(jl)   = zzdrs                                   !---wiso-code: store old value of sedimentation
           zsnmlt_tmp(jl) = zsnmlt                                   !---wiso-code: store amount of snowmelt
           zzdrr       = zzdrr+zsnmlt
           zzdrs       = zzdrs-zsnmlt
           zsmlt(jl)   = zsmlt(jl)+zsnmlt/(zcons2*zdp(jl))
        END IF
        zpretot        = zrfl(jl)+zsfl(jl)
        zpredel        = zzdrr+zzdrs
        lo=(zpretot .GT. zpredel)
        zclcpre(jl)    = MERGE(zclcpre(jl),zclcaux(jl),lo)
        zpresum        = zpretot+zpredel
        IF (zpresum .GT. cqtmin) THEN
           zclcpre(jl) = MAX(zclcpre(jl),(zclcaux(jl)*zpredel          &
                                         +zclcpre(jl)*zpretot)/zpresum)
           zclcpre(jl) = MIN(zclcpre(jl),1.0_dp)
           zclcpre(jl) = MAX(zclcpre(jl),0.0_dp)
        ELSE
           zclcpre(jl) = 0.0_dp
        END IF
        zrfl(jl)       = zrfl(jl)+zzdrr-zcons2*zdp(jl)*zevp(jl)
        zsfl(jl)       = zsfl(jl)+zzdrs-zcons2*zdp(jl)*zsub(jl)
!
701  END DO

!---wiso-code
  IF (lwiso) THEN

!           ------------------------------------------------------------------
!           Cloud physics and precipitation fluxes at the surface - water isotopes

        DO jt = 1,kwiso
          DO jl = 1,kproma
            locc(jl)     = zclcaux(jl) .GT. 0.0_dp
            zwisoxlb(jl,jt) = MAX(zwisoxlb(jl,jt),1.e-20_dp)
            zwisoxib(jl,jt) = MAX(zwisoxib(jl,jt),1.e-20_dp)

            IF (locc(jl) .AND. (zxlb_tmp(jl) > cqtmin .OR. zxib_tmp(jl) > cqtmin)) THEN
            
!           For both warm and cold clouds: 
!           - assume no further fractionation for different cloud coalescence and conversion processes
!

!             warm clouds
              zdelta=tnat(jt)
              IF (zxlb_tmp(jl).GT.cwisomin) zdelta = zwisoxlb(jl,jt)/zxlb_tmp(jl)
              IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
              zwisoxlb(jl,jt) = zwisoxlb(jl,jt)-zxlb_difw(jl)*zdelta
              zwisorpr(jl,jt) = zwisorpr(jl,jt)+(zrpr(jl)-zrpr_tmp(jl))*zdelta
!             cold clouds
              zwisosacl(jl,jt)= zsacl(jl)*zdelta
              zwisoxlb(jl,jt) = zwisoxlb(jl,jt)-zwisosacl(jl,jt)
              zdelta=tnat(jt)
              IF (zxib_tmp(jl).GT.cwisomin) zdelta = zwisoxib(jl,jt)/zxib_tmp(jl)
              IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
              zwisoxib(jl,jt) = zwisoxib(jl,jt)-zxib_difc(jl)*zdelta
              zwisospr(jl,jt) = zwisospr(jl,jt)+(zspr(jl)-zspr_tmp(jl))*zdelta
              
            ENDIF

!           Updating precipitation fluxes - water isotopes
!           In the lowest layer (klev), the sedimentation sink of cloud ice 
!           is balanced by precipitation at the surface (through 'zwisozdrs').
!           - assume no fractionation during melting due to low diffusivity of water isotopes in ice

            zwisozdrr         = zcons2*zdp(jl)*zwisorpr(jl,jt)
            zwisozdrs         = zcons2*zdp(jl)*(zwisospr(jl,jt)+zwisosacl(jl,jt))
            IF (jk .EQ. klev) THEN
              zwisozdrs       = zwisozdrs+zwisoxiflux(jl,jt)
              zdelta=tnat(jt)
              IF (zzdrs_tmp(jl).GT.cwisomin) zdelta = zwisozdrs/zzdrs_tmp(jl)
              IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta = 1._dp
              zwisosnmlt      = zsnmlt_tmp(jl)*zdelta
              zwisozdrr       = zwisozdrr+zwisosnmlt
              zwisozdrs       = zwisozdrs-zwisosnmlt
            END IF

            zwisorfl(jl,jt)   = zwisorfl(jl,jt)+zwisozdrr-zcons2*zdp(jl)*zwisoevp(jl,jt)
            zwisosfl(jl,jt)   = zwisosfl(jl,jt)+zwisozdrs-zcons2*zdp(jl)*zwisosub(jl,jt)

          END DO
        END DO

  END IF
!---wiso-code-end

!
!     ------------------------------------------------------------------
!       8.    Updating tendencies of t, q, xl, xi and final cloud cover
!
     DO 811 jl = 1,kproma
!
!       8.10   Cloud cover scheme tendencies
!
        IF (lcover .AND. jk.GE.ncctop) THEN
           locc(jl)        = zclcaux(jl) .GT. 0.0_dp
!
!          Source terms from convection
!          Skewness:
!
           zconvskew(jl)   = cbeta_cs * (pxtec(jl,jk)+pqtec(jl,jk))    &
                                  /pbetass(jl,jk)
           zconvskew(jl)   = MIN(zconvskew(jl),                        &
                                (cbeta_pq_max-pxskew(jl,jk))/zdtime)
!
!          Convective width now diagnosed, assuming 'a' unchanged:
!
           IF (pqm1(jl,jk) >= pbetass(jl,jk)) THEN
              zskewp1      = pxskew(jl,jk)+zconvskew(jl)*zdtime
              zbbap1       = zwide(jl)*(cbeta_pq+zskewp1)/             &
                                       (cbeta_pq+pxskew(jl,jk))
              zconvvar(jl) = (zbbap1-zwide(jl))/zdtime
           ELSE
              zconvvar(jl) = 0.0_dp
           ENDIF
!
!       8.11 Simple linearized effect of microphysics on skewness
!
           IF (pbetaa(jl,jk) < pbetass(jl,jk) .AND.                    &
               pbetab(jl,jk) > pbetass(jl,jk)) THEN
              zmdelb = (zxlte(jl)+zxite(jl))*ztmst                     &
                       -zrpr(jl)-zsacl(jl)-zspr(jl)+zcnd(jl)+zdep(jl)  &
                       +zgenti(jl)+zgentl(jl)
              zmdelb = MAX(0.0_dp,MIN(1.0_dp,-zmdelb/                  &
                                              MAX(zepsec,zbetacl(jl))))
              zmdelb = (pbetass(jl,jk)-pbetab(jl,jk))*zmdelb
              zmqp1  = (pbetab(jl,jk)+zmdelb-pbetaa(jl,jk))            &
                        *cbeta_pq/(zbetaqt(jl)-pbetaa(jl,jk))          &
                                                          - cbeta_pq
              zmqp1  = MAX(MIN(zmqp1,cbeta_pq_max),cbeta_pq)
              zmicroskew(jl) = MIN(0.0_dp,(zmqp1-pxskew(jl,jk))/zdtime)
           ENDIF
!
!       8.2   New skewness and variance
!
           zxskewte(jl)    = zconvskew(jl)                             &
                             +zmicroskew(jl)+zturbskew(jl)
           zxvarte(jl)     = zconvvar(jl)+zturbvar(jl)
!
           zvarp1          = pxvar(jl,jk)+zxvarte(jl)*zdtime
           zskewp1         = pxskew(jl,jk)+zxskewte(jl)*zdtime
!
           pxskew(jl,jk)   = MAX(MIN(zskewp1,cbeta_pq_max),cbeta_pq)
           zvarmx          = zbetaqt(jl)*(1._dp+pxskew(jl,jk)/cbeta_pq)
           pxvar(jl,jk)    = MAX(MIN(zvarp1,zvarmx),zvartg(jl))
!
        ENDIF !lcover and jk >= ncctop
!
!       8.3   Tendencies of thermodynamic variables
!             Attn: The terms zxisub and zximlt do not appear in
!                   pxite because these processes have already been
!                   included in pxite via changes in cloud ice
!                   sedimentation (see 3.1, 3.2 and 4)
!
        pqte(jl,jk)  = pqte(jl,jk)                                     &
                        +(-zcnd(jl)-zgentl(jl)+zevp(jl)+zxlevap(jl)    &
                          -zdep(jl)-zgenti(jl)+zsub(jl)+zxievap(jl)    &
                                          +zxisub(jl))/ztmst
        ptte(jl,jk)  = ptte(jl,jk)+(zlvdcp(jl)                         &
                        *(zcnd(jl)+zgentl(jl)-zevp(jl)-zxlevap(jl))    &
                                  +zlsdcp(jl)                          &
                        *(zdep(jl)+zgenti(jl)-zsub(jl)-zxievap(jl)     &
                        -zxisub(jl))+(zlsdcp(jl)-zlvdcp(jl))           &
                        *(-zsmlt(jl)-zimlt(jl)-zximlt(jl)+zfrl(jl)     &
                                           +zsacl(jl)))/ztmst
        pxlte(jl,jk) = pxlte(jl,jk)+zxlte(jl)                          &
                        +(zimlt(jl)+zximlt(jl)-zfrl(jl)-zrpr(jl)       &
                        -zsacl(jl)+zcnd(jl)+zgentl(jl)-zxlevap(jl))    &
                                                       /ztmst
        pxite(jl,jk) = pxite(jl,jk)+zxite(jl)+(zfrl(jl)-zspr(jl)       &
                          +zdep(jl)+zgenti(jl)-zxievap(jl))/ztmst
        ztp1         = ptm1(jl,jk)+ptte(jl,jk)*ztmst
        zqp1         = pqm1(jl,jk)+pqte(jl,jk)*ztmst
        zxlp1        = pxlm1(jl,jk)+pxlte(jl,jk)*ztmst
        zxip1        = pxim1(jl,jk)+pxite(jl,jk)*ztmst
!
!       8.4   Corrections: Avoid negative cloud water/ice
!
        zxlold = zxlp1
        lo             = (zxlp1 .LT. ccwmin)
        lo_xlcor(jl)   = lo                                !---wiso-code: store value of lo for isotope correction
        zxlp1          = MERGE(0.0_dp,zxlp1,lo)
        zdxlcor        = (zxlp1-zxlold)/ztmst
        zxiold         = zxip1
        lo1            = (zxip1 .LT. ccwmin)
        lo_xicor(jl)   = lo1                               !---wiso-code: store value of lo1 for isotope correction
        zxip1          = MERGE(0.0_dp,zxip1,lo1)
        zdxicor        = (zxip1-zxiold)/ztmst
        paclc(jl,jk)   = MERGE(0.0_dp,paclc(jl,jk),lo.AND.lo1)
        paclcac(jl,jk) = paclcac(jl,jk)+paclc(jl,jk)*zdtime
        pxlte(jl,jk)   = pxlte(jl,jk)+zdxlcor
        pxite(jl,jk)   = pxite(jl,jk)+zdxicor
        pqte(jl,jk)    = pqte(jl,jk)-zdxlcor-zdxicor
        ptte(jl,jk)    = ptte(jl,jk)+zlvdcp(jl)*zdxlcor                &
                                    +zlsdcp(jl)*zdxicor
!
811  END DO

!---wiso-code
  IF (lwiso) THEN

!           ------------------------------------------------------------------
!           Updating tendencies of q, xl, xi - water isotopes

        DO jt = 1,kwiso
          DO jl = 1,kproma

!           Tendencies of thermodynamic variables - water isotopes
!           Attn: The terms zwisoxisub and zwisoximlt do not appear in
!                   pwisoxite because these processes have already been
!                   included in pwisoxite via changes in cloud ice
!                   sedimentation (see 3.1, 3.2 and 4)

            pwisoqte(jl,jk,jt)  = pwisoqte(jl,jk,jt)                                     &
                                 +(-zwisocnd(jl,jt)-zwisogentl(jl,jt)+zwisoevp(jl,jt)+zwisoxlevap(jl,jt)    &
                                   -zwisodep(jl,jt)-zwisogenti(jl,jt)+zwisosub(jl,jt)+zwisoxievap(jl,jt)    &
                                   +zwisoxisub(jl,jt))/ztmst
            pwisoxlte(jl,jk,jt) = pwisoxlte(jl,jk,jt)+zwisoxlte(jl,jt)                          &
                                 +(zwisoimlt(jl,jt)+zwisoximlt(jl,jt)-zwisofrl(jl,jt)-zwisorpr(jl,jt)       &
                                 -zwisosacl(jl,jt)+zwisocnd(jl,jt)+zwisogentl(jl,jt)-zwisoxlevap(jl,jt))/ztmst
            pwisoxite(jl,jk,jt) = pwisoxite(jl,jk,jt)+zwisoxite(jl,jt)                                      &
                                 +(zwisofrl(jl,jt)-zwisospr(jl,jt)                                          &
                                 +zwisodep(jl,jt)+zwisogenti(jl,jt)-zwisoxievap(jl,jt))/ztmst
                                 

            zwisoqp1            = pwisoqm1(jl,jk,jt) +pwisoqte(jl,jk,jt) *ztmst
            zwisoxlp1           = pwisoxlm1(jl,jk,jt)+pwisoxlte(jl,jk,jt)*ztmst
            zwisoxip1           = pwisoxim1(jl,jk,jt)+pwisoxite(jl,jk,jt)*ztmst
            
!           Corrections: Avoid negative cloud water/ice - water isotopes
!
            zwisoxlold = zwisoxlp1
            zwisoxlp1          = MERGE(0.0_dp,zwisoxlp1,lo_xlcor(jl))
            zwisodxlcor        = (zwisoxlp1-zwisoxlold)/ztmst
            zwisoxiold         = zwisoxip1
            zwisoxip1          = MERGE(0.0_dp,zwisoxip1,lo_xicor(jl))
            zwisodxicor        = (zwisoxip1-zwisoxiold)/ztmst
            pwisoxlte(jl,jk,jt)= pwisoxlte(jl,jk,jt)+zwisodxlcor
            pwisoxite(jl,jk,jt)= pwisoxite(jl,jk,jt)+zwisodxicor
            pwisoqte(jl,jk,jt) = pwisoqte(jl,jk,jt)-zwisodxlcor-zwisodxicor
            
          END DO
        END DO

!           Isotopic equilibration of vapour, liquid water and rainfall
!
        DO jl=1,kproma
          ztpone(jl)=ptm1(jl,jk)+ptte(jl,jk)*ztmst 
          zqpone(jl)=pqm1(jl,jk)+pqte(jl,jk)*ztmst
          zqspone(jl)=zqsp1tmp_tmp(jl)
          zxlpone(jl)=pxlm1(jl,jk)+pxlte(jl,jk)*ztmst
        END DO

        CALL cloudadjwiso(kproma,kbdim,klev,kwiso,jk,             &
                          ztmst,zcons2,zdp,                       &
                          ztpone,zqpone,zqspone,zxlpone,          &
                          paclc,                                  &
                          pwisoqm1,  pwisoqte,                    &
                          pwisoxlm1, pwisoxlte,                   &
                          zrfl, zwisorfl)   

  END IF
!---wiso-code-end

831 END DO    ! Vertical loop

!
!     ------------------------------------------------------------------
!       9.    Diagnostics
!
!       9.1   Accumulated precipitation at the surface
!
  DO 911 jl    = 1,kproma
     prsfl(jl) = zrfl(jl)
     pssfl(jl) = zsfl(jl)
     paprl(jl) = paprl(jl)+zdtime*(prsfl(jl)+pssfl(jl))
     paprs(jl) = paprs(jl)+zdtime*pssfl(jl)
911 END DO
!
!       9.2   Total cloud cover
!
    DO 921 jl    = 1,kproma
      zclcov(jl) = 1.0_dp-paclc(jl,1)
921 END DO
    DO 923 jk      = 2,klev
      DO 922 jl    = 1,kproma
        zclcov(jl) = zclcov(jl)                                        &
                            *(1._dp-MAX(paclc(jl,jk),paclc(jl,jk-1)))  &
                            /(1._dp-MIN(paclc(jl,jk-1),zxsec))
922   END DO
923 END DO
    DO 924 jl     = 1,kproma
      zclcov(jl)  = 1.0_dp-zclcov(jl)
      paclcov(jl) = paclcov(jl)+zdtime*zclcov(jl)
924 END DO
!
!       9.3   Vertical integrals of humidity, cloud water and cloud ice
!
    DO 931 jl   = 1,kproma
      zqvi(jl)  = 0.0_dp
      zxlvi(jl) = 0.0_dp
      zxivi(jl) = 0.0_dp
931 END DO
!
    DO 933 jk     = ktdia,klev
      DO 932 jl   = 1,kproma
        zdpg      = (paphm1(jl,jk+1)-paphm1(jl,jk))/g
        zqvi(jl)  = zqvi(jl)+pqm1(jl,jk)*zdpg
        zxlvi(jl) = zxlvi(jl)+pxlm1(jl,jk)*zdpg
        zxivi(jl) = zxivi(jl)+pxim1(jl,jk)*zdpg
932   END DO
933 END DO
!
    DO 934 jl   = 1,kproma
      pqvi(jl)  = pqvi(jl)+zdtime*zqvi(jl)
      pxlvi(jl) = pxlvi(jl)+zdtime*zxlvi(jl)
      pxivi(jl) = pxivi(jl)+zdtime*zxivi(jl)
934 END DO
!

!---wiso-code
  IF (lwiso) THEN

!       Diagnostics - Water Isotopes
!
!       Accumulated precipitation at the surface
!
    DO jt = 1,kwiso
      DO jl = 1,kproma
        pwisorsfl(jl,jt) = zwisorfl(jl,jt)
        pwisossfl(jl,jt) = zwisosfl(jl,jt)
        pwisoaprl(jl,jt) = pwisoaprl(jl,jt)+zdtime*(pwisorsfl(jl,jt)+pwisossfl(jl,jt))
        pwisoaprs(jl,jt) = pwisoaprs(jl,jt)+zdtime*pwisossfl(jl,jt)
      END DO
    END DO

!       Vertical integrals of humidity, cloud water and cloud ice
!
    DO jt = 1,kwiso
      DO jl = 1,kproma
        zwisoqvi(jl,jt)  = 0.0_dp
        zwisoxlvi(jl,jt) = 0.0_dp
        zwisoxivi(jl,jt) = 0.0_dp
      END DO
    END DO

    DO jt = 1,kwiso
      DO jk = ktdia,klev
        DO jl = 1,kproma
          zdpg             = (paphm1(jl,jk+1)-paphm1(jl,jk))/g
          zwisoqvi(jl,jt)  = zwisoqvi(jl,jt) +pwisoqm1(jl,jk,jt) *zdpg
          zwisoxlvi(jl,jt) = zwisoxlvi(jl,jt)+pwisoxlm1(jl,jk,jt)*zdpg
          zwisoxivi(jl,jt) = zwisoxivi(jl,jt)+pwisoxim1(jl,jk,jt)*zdpg
        END DO
      END DO
    END DO

    DO jt = 1,kwiso
      DO jl = 1,kproma
        pwisoqvi(jl,jt)  = pwisoqvi(jl,jt)+zdtime*zwisoqvi(jl,jt)
        pwisoxlvi(jl,jt) = pwisoxlvi(jl,jt)+zdtime*zwisoxlvi(jl,jt)
        pwisoxivi(jl,jt) = pwisoxivi(jl,jt)+zdtime*zwisoxivi(jl,jt)
      END DO
    END DO

  END IF
!---wiso-code-end

  RETURN
END SUBROUTINE cloud
