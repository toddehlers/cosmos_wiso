SUBROUTINE cover (         kproma,   kbdim, ktdia, klev, klevp1        &
!
! - INPUT  1D .
                         , ktype,    pfrw                              &
! - OUTPUT 1D .
                         , knvb,     printop                           &
! - INPUT  2D .
                         , paphm1,   papm1                             &
                         , pqm1,     ptm1,     pxlm1                   &
                         , pxim1,    pvervel                           &
                         , pxvar,    pxskew                            &
! - OUTPUT 2D .
                         , paclc                                       &
                         , pbetaa,   pbetab                            &
                         , pbetass                                     &
       )
!-----------------------------------------------------------------------
!     *Cover*
!          Diagnoses cloud cover for current timestep
!
!     Subject.
!     --------
!     In ECHAM4 all cloud cover calculations were performed
!     in routine CLOUD.  The routine diagnosed the cloud cover
!     from the previous timestep using m1 variables, but then
!     used a "first guess" calculation of T and Q to give a
!     cloud cover estimate for the next timestep.  This meant
!     that the radiation scheme used cloud cover values that were
!     from a different timestep to the temperature and water vapour
!     values, and also that the T and Q values were anyway
!     preliminary.  Finally, the cover calculation was performed
!     twice each timestep, when one calculation suffices.
!
!     This scheme calculates cover diagnostically and is called
!     once at the beginning of each timestep.  It uses the
!     standard relative humidity calculation from Lohmann and
!     Roeckner (96), or the method from the new prognostic
!     scheme of Tompkins.  The choice of which scheme to use is
!     controlled by the parameter switch LCOVER, which is set in
!     namelist PHYSCTL along with lsurf etc... Note that even if
!     lcover=.false. (RH scheme) you can't restart this model version
!     from restart files saved from a different model version, since
!     the two extra prognostic equations, pxvar and pxskew are still
!     stored even though they are not actively used.  However, this means
!     that once you have restart files from this version, you are able
!     change lcover at will.
!
!     In the new scheme the variable xskew is provided
!     as outlined in the reference, this variable represents
!     directly the Beta distribution shape parameter "q"
!     The shape parameter "p" (zbetap) a tunable parameter and it is
!     recommended that this be set to a low constant 1.5<p<2.0
!     (This may be changed later to prognostic to allow negative skewness
!     from downdraft detrainment, see ref. for details).
!
!     from xi,xl,q,xskew and zbetap, the Beta distribution is defined
!     and cloud cover is diagnosable.  For the iteration, Ridders' method
!     is used (see Numerical Recipes).
!
!     Attention: 
!     In the current version the advective tendencies of skewness 
!     and variance are set to zero.
!
!     INTERFACE.
!     ----------
!
!     *Call cover*
!
!     Input arguments.
!     ----- ----------
!  - 1D
!  ktype    : type of convection
!  - 2D
!  paphm1   : pressure at half levels                              (n-1)
!  papm1    : pressure at full levels                              (n-1)
!  pqm1     : specific humidity                                    (n-1)
!  ptm1     : temperature                                          (n-1)
!  pxlm1    : cloud liquid water                                   (n-1)
!  pxim1    : cloud ice                                            (n-1)
!  pxvar    : the beta distribution width "b-a"                    (n-1)
!  pxskew   : the beta distribution shape parameter "q"            (n-1)
!  pbetaa   : the beta distribution minimum a                      (n-1)
!  pbetab   : the beta distribution maximum b                      (n-1)
!
!     Output arguments.
!     ------ ----------
!  - 2D
!  paclc    : cloud cover
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
!     Diagnostic CC scheme: Lohmann and Roeckner 96, Clim. Dyn.
!     Prognostic CC scheme: Tompkins 2002, J. Atmos. Sci. 
!
!     Authors.
!     -------
!     A. Tompkins    MPI-Hamburg  2000
!     K. Ketelesen   NEC, April 2002
!     L. Kornblueh   MPI, April 2002
!
!     Modifications.
!     --------------
!     view cover
!       v2: first working version
!       v5: lookup table added
!       v8: zriddr and functions replaced for vectorization
!       v9: optimizations, longer vector loop and less indirect addressing
!          - introduction of additional arrays 
!          - change structure if "beta function scheme" IF block
!          - scattered loops "ictit" and "ictdg" are collected over "kproma" and "klev"
!          - intoducing additional arrays to hold data in "ictit" and "ictdg"
!            addressing scheme
!
  USE mo_kind,            ONLY : dp
  USE mo_constants,       ONLY : rd, cpd, vtmpc1
  USE mo_convect_tables,  ONLY : lookuperror, lookupoverflow           &
                               , jptlucu1, jptlucu2                    &
                               , tlucua, tlucuaw
  USE mo_cloud,           ONLY : jbmin1, ncctop, cqtmin, cbeta_pq      &
                               , cthomi, tmelt, csecfrl, jbmin, jbmax  &
                               , csatsc, ccwmin, cbeta_pq_max, nbetaq  &
                               , cbetaqs, rbetak, nbetax, tbetai0      &
                               , tbetai1, cvarmin, cmmrmax, crt, crs   &
                               , nex, crhsc 
  USE mo_param_switches,  ONLY : lcover
  
  IMPLICIT NONE

  INTEGER :: kbdim, klevp1, klev, kproma, jl, jk, ktdia, it, jb, jt
  REAL(dp):: rdcpd, zdthdp, zcor, zbetai0, zbetai1, zrhc, zsat, zqr
!   Input arrays
  INTEGER::  ktype(kbdim)
  REAL(dp) ::    pfrw(kbdim)
  REAL(dp) ::    paphm1(kbdim,klevp1)                                  &
            ,papm1(kbdim,klev)    ,pqm1(kbdim,klev)                    &
            ,ptm1(kbdim,klev)     ,pxlm1(kbdim,klev)                   &
            ,pxim1(kbdim,klev)    ,pvervel(kbdim,klev)                 &
            ,pxvar(kbdim,klev)    ,pxskew(kbdim,klev)
  REAL(dp) ::    pbetass(kbdim,klev)

!   Output arrays
  REAL(dp) ::    paclc(kbdim,klev)                                     &
            ,pbetaa(kbdim,klev) ,pbetab(kbdim,klev)
  REAL(dp) ::    printop(kbdim)
  INTEGER :: knvb(kbdim)

!
!   Temporary arrays
!
  REAL(dp)   ::  zdthmin(kbdim)      ,ztheta(kbdim,klev)
!
!   Pointers and counters for iteration and diagnostic loop:
!
  INTEGER :: iptit(kproma*klev), iptdg(kproma*klev)                    &
           , iqidx_l1(kproma*klev)
  INTEGER :: ictdg, ixidx, ictit, icount
  INTEGER :: iqidx(kproma,klev)
!
!   variables required for zriddr iteration scheme:
!
  REAL(dp) :: fh(kproma*klev), fl(kproma*klev), fm(kproma*klev)        &
            , fnew(kproma*klev), xh(kproma*klev), xl(kproma*klev)      &
            , xm(kproma*klev), xnew(kproma*klev)

  REAL(dp) :: unused = -1.11e30_dp
  REAL(dp) :: ztt, zx1, zx2, zvar, zvartarget, zss, zqt, zpp, zqq      &
            , zaa, zbb, zsqarg
  REAL(dp) :: zqsm1(kproma*klev), zskew1(kproma*klev)                  &
            , zskew2(kproma*klev)
  REAL(dp) :: zbetaqt(kproma,klev), zbetacl(kproma,klev)

  LOGICAL :: lo, lo1, lo2
  LOGICAL :: liter(kproma*klev)

  REAL(dp) :: pqm1_l1(kproma*klev), pqm1_l2(kproma*klev)
  REAL(dp) :: zbetaqt_l1(kproma*klev), zbetacl_l1(kproma*klev)         &
            , zbetass_l1(kproma*klev)
  REAL(dp) :: zbetaqt_l2(kproma*klev), zbetacl_l2(kproma*klev)         &
            , zbetass_l2(kproma*klev)

  INTEGER :: jk1(kproma*klev), jk2(kproma*klev)

!   number of iteration loops.  In a small test, most gridpoints
!   converged in 2 to 4 iterations, with occassionally as many as
!   8 needed for an accuracy of 0.01*(ql+qi). Don't set too high since
!   the loops are all performed, albeit often as "empty loops" (i.e.
!   LITER is false because of convergence), due to vectorisation needs.
!
  INTEGER, PARAMETER :: niter=60            ! max loop number
  REAL(dp), PARAMETER :: ziter_acc=0.01_dp  ! iteration accuracy*(ql+qi)
!
!  Executable statements
!
  lookupoverflow = .FALSE.
!
!   Computational constants
!
  rdcpd  = rd/cpd
!
!   Initialize variables
!
  DO jl=1,kproma
    zdthmin(jl)=0.0_dp
    knvb(jl)=1
    printop(jl)=1.0_dp
  END DO
!
  DO jk = ktdia,klev
     DO jl = 1,kproma
        ztheta(jl,jk) = ptm1(jl,jk)*(1.0e5_dp/papm1(jl,jk))**rdcpd
     END DO
  END DO
!
!       1.3   Checking occurrence of low-level inversion
!             (below 1000 m, sea points only, no convection)
!
  DO jk = klev,jbmin1,-1
     DO jl = 1,kproma
        IF (pfrw(jl) .GT. 0.0_dp .AND. ktype(jl) .EQ. 0) THEN
           zdthdp      = (ztheta(jl,jk)-ztheta(jl,jk-1))               &
                                      /(papm1(jl,jk)-papm1(jl,jk-1))
           lo          = zdthdp .LT. zdthmin(jl)
           zdthmin(jl) = MERGE(zdthdp,zdthmin(jl),lo)
           knvb(jl)    = MERGE(jk,knvb(jl),lo)
        END IF
     END DO
  END DO

!
!   Tunable parameters now in mo_cloud.f90
!
  IF (ncctop > 1) THEN
    DO jk=1,ncctop-1
      DO jl=1,kproma
        pbetaa(jl,jk)=0.0_dp
        pbetab(jl,jk)=0.0_dp
        pxvar(jl,jk)=cqtmin
        pxskew(jl,jk)=cbeta_pq
      END DO
    END DO
  END IF
!
  ictit=0
  ictdg=0
  DO jk=ktdia,klev
    IF (lcover .AND. jk >= ncctop) THEN  ! beta function scheme
!
!       1.   Calculate the saturation mixing ratio
!
      DO jl=1,kproma
        lo2=(ptm1(jl,jk) < cthomi).OR.(ptm1(jl,jk) < tmelt             &
             .AND.pxim1(jl,jk) > csecfrl)
        it=NINT(ptm1(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsm1(jl)=MERGE(tlucua(it),tlucuaw(it),lo2)/papm1(jl,jk)
        zqsm1(jl)=MIN(zqsm1(jl),0.5_dp)
        zcor=1._dp/(1._dp-vtmpc1*zqsm1(jl))
        zqsm1(jl)=zqsm1(jl)*zcor
        jb=knvb(jl)
        lo=(jb.GE.jbmin .AND. jb.LE.jbmax .AND. pvervel(jl,jk).GT.0._dp)
        lo1=(jk.EQ.jb .OR. jk.EQ.jb+1)
        IF(lo .AND. lo1) THEN
          zqsm1(jl)=zqsm1(jl)*csatsc
          printop(jl)=REAL(jb,dp)
        ENDIF
      ENDDO
!
      IF (lookupoverflow) CALL lookuperror ('cover       ')
!
!       2.    calculate cloud cover
!
!   Don't need to iterate at every gridpoint, thus make pointers
!
!DIR$ CONCURRENT
      DO jl=1,kproma
        zbetacl(jl,jk)=MAX(0._dp,pxlm1(jl,jk))+MAX(0._dp,pxim1(jl,jk))
        zbetaqt(jl,jk)=MAX(cqtmin,pqm1(jl,jk))+zbetacl(jl,jk)
        pbetass(jl,jk)=MAX(pqm1(jl,jk),zqsm1(jl)) !safety
        lo2=((pxim1(jl,jk)>ccwmin .OR. pxlm1(jl,jk)>ccwmin) .AND.      &
             pqm1(jl,jk)<pbetass(jl,jk))
        IF (lo2) THEN
          ictit=ictit+1
          iptit(ictit)=jl
          liter(ictit)=.true.
          pqm1_l1(ictit)    = pqm1(jl,jk)
          jk1(ictit)        = jk
          zbetaqt_l1(ictit) = zbetaqt(jl,jk)
          zbetacl_l1(ictit) = zbetacl(jl,jk)
          zbetass_l1(ictit) = pbetass(jl,jk)
          zskew1(ictit)=MAX(MIN(pxskew(jl,jk),cbeta_pq_max),cbeta_pq)
          iqidx(jl,jk)=INT((nbetaq/cbetaqs)*                           &
               LOG((zskew1(ictit)-cbeta_pq)/rbetak+1.)+0.5_dp)
          iqidx_l1(ictit)   = iqidx(jl,jk)
        ELSE
          ictdg=ictdg+1
          iptdg(ictdg)=jl
          pqm1_l2(ictdg)    = pqm1(jl,jk)
          jk2(ictdg)        = jk
          zbetaqt_l2(ictdg) = zbetaqt(jl,jk)
          zbetacl_l2(ictdg) = zbetacl(jl,jk)
          zbetass_l2(ictdg) = pbetass(jl,jk)
          zskew2(ictdg)=MAX(MIN(pxskew(jl,jk),cbeta_pq_max),cbeta_pq)
          iqidx(jl,jk)=INT((nbetaq/cbetaqs)*                           &
               LOG((zskew2(ictdg)-cbeta_pq)/rbetak+1.)+0.5_dp)
        ENDIF
!
!   setup index for skewness, q, in lookup table (doesn't change)
!
      ENDDO
    ENDIF !lcover
  END DO
!
!
!   Partially cloudy: Iterative gridpoints
!   uses func using ridders' method, return the root of a function
!   Method: See Numerical Recipes
!
!
  DO jl=1,ictit
!
!   Lower bound < 0 to catch occassional overflow
!
    zx1=-0.1_dp
    ztt=(zbetass_l1(jl)-zx1)*cbeta_pq/                                 &
         ((zbetaqt_l1(jl)-zx1)*                                        &
         (cbeta_pq+zskew1(jl)))
    ztt=nbetax*MAX(MIN(ztt,1.0_dp),0.0_dp)
    ixidx=INT(ztt)
    IF (ixidx==nbetax) THEN
      zbetai0=1.0_dp
      zbetai1=1.0_dp
    ELSE
      zbetai0=(ztt-ixidx)*tbetai0(iqidx_l1(jl),ixidx+1)                &
           +(ixidx+1._dp-ztt)*tbetai0(iqidx_l1(jl),ixidx)
      zbetai1=(ztt-ixidx)*tbetai1(iqidx_l1(jl),ixidx+1)                &
           +(ixidx+1._dp-ztt)*tbetai1(iqidx_l1(jl),ixidx)
    ENDIF
    fl(jl)=(zbetaqt_l1(jl)-zx1)*zbetai1 -                              &
         (zbetass_l1(jl)-zx1)*zbetai0 +                                &
         zbetass_l1(jl)-MAX(cqtmin,pqm1_l1(jl))
!
!   3 conditions for the maximum iteration bound: a<qs,a<qt,b>qs
!
    zx2=MIN(((cbeta_pq+zskew1(jl))*                                    &
         zbetaqt_l1(jl)-cbeta_pq*zbetass_l1(jl))                       &
         /zskew1(jl),                                                  &
         zbetaqt_l1(jl),zbetass_l1(jl))
    ztt=(zbetass_l1(jl)-zx2)*cbeta_pq/                                 &
         ((zbetaqt_l1(jl)-zx2)*                                        &
         (cbeta_pq+zskew1(jl)))
    ztt=nbetax*MAX(MIN(ztt,1.0_dp),0.0_dp)
    ixidx=INT(ztt)
    IF (ixidx==nbetax) THEN
      zbetai0=1.0_dp
      zbetai1=1.0_dp
    ELSE
      zbetai0=(ztt-ixidx)*tbetai0(iqidx_l1(jl),ixidx+1)                &
           +(ixidx+1._dp-ztt)*tbetai0(iqidx_l1(jl),ixidx)
      zbetai1=(ztt-ixidx)*tbetai1(iqidx_l1(jl),ixidx+1)                &
           +(ixidx+1._dp-ztt)*tbetai1(iqidx_l1(jl),ixidx)
    ENDIF
    fh(jl)=(zbetaqt_l1(jl)-zx2)*zbetai1 -                              &
         (zbetass_l1(jl)-zx2)*zbetai0 +                                &
         zbetass_l1(jl)-MAX(cqtmin,pqm1_l1(jl))
    xl(jl)=zx1
    xh(jl)=zx2
    if (xl(jl)==xh(jl)) liter(jl)=.false.
    xnew(jl)=unused
  ENDDO
!
  DO jt=1,niter   ! short iteration loop
    DO jl=1,ictit   ! longer loop over indices
      IF (liter(jl)) THEN
        xm(jl)=0.5_dp*(xl(jl)+xh(jl))
        ztt=(zbetass_l1(jl)-xm(jl))*cbeta_pq/                          &
             ((zbetaqt_l1(jl)-xm(jl))*                                 &
             (cbeta_pq+zskew1(jl)))
        ztt=nbetax*MAX(MIN(ztt,1.0_dp),0.0_dp)
        ixidx=INT(ztt)
        IF (ixidx==nbetax) THEN
          zbetai0=1.0_dp
          zbetai1=1.0_dp
        ELSE
          zbetai0=(ztt-ixidx)*tbetai0(iqidx_l1(jl),ixidx+1)            &
               +(ixidx+1._dp-ztt)*tbetai0(iqidx_l1(jl),ixidx)
          zbetai1=(ztt-ixidx)*tbetai1(iqidx_l1(jl),ixidx+1)            &
               +(ixidx+1._dp-ztt)*tbetai1(iqidx_l1(jl),ixidx)
        ENDIF
        fm(jl)=(zbetaqt_l1(jl)-xm(jl))*zbetai1-                        &
             (zbetass_l1(jl)-xm(jl))*zbetai0+                          &
             zbetass_l1(jl)-MAX(cqtmin,pqm1_l1(jl))
        zsqarg=MAX(fm(jl)**2-fl(jl)*fh(jl),0._dp)
        ztt=SQRT(zsqarg)
        IF (ztt.eq.0._dp) THEN
          ztt=xnew(jl)
        ELSE
          ztt=xm(jl)+(xm(jl)-xl(jl))                                   &
               *(sign(1._dp,fl(jl)-fh(jl))*fm(jl)/ztt) !update formula
        ENDIF
        IF (abs(xnew(jl)-ztt).le.ziter_acc*zbetacl_l1(jl))             &
             liter(jl)=.false.
        xnew(jl)=ztt
        ztt=(zbetass_l1(jl)-xnew(jl))*cbeta_pq/                        &
             ((zbetaqt_l1(jl)-xnew(jl))*                               &
             (cbeta_pq+zskew1(jl)))
        ztt=nbetax*MAX(MIN(ztt,1.0_dp),0.0_dp)
        ixidx=INT(ztt)
        IF (ixidx==nbetax) THEN
          zbetai0=1.0_dp
          zbetai1=1.0_dp
        ELSE
          zbetai0=(ztt-ixidx)*tbetai0(iqidx_l1(jl),ixidx+1)            &
               +(ixidx+1._dp-ztt)*tbetai0(iqidx_l1(jl),ixidx)
          zbetai1=(ztt-ixidx)*tbetai1(iqidx_l1(jl),ixidx+1)            &
               +(ixidx+1._dp-ztt)*tbetai1(iqidx_l1(jl),ixidx)
        ENDIF
        fnew(jl)=(zbetaqt_l1(jl)-xnew(jl))*zbetai1-                    &
             (zbetass_l1(jl)-xnew(jl))*zbetai0+                        &
             zbetass_l1(jl)-MAX(cqtmin,pqm1_l1(jl))
!
!   bookkeeping to keep the root bracketed on next iteration.
!
        IF (sign(fm(jl),fnew(jl)).ne.fm(jl)) THEN
          xl(jl)=xm(jl)
          fl(jl)=fm(jl)
          xh(jl)=xnew(jl)
          fh(jl)=fnew(jl)
        ELSE IF(sign(fl(jl),fnew(jl)).ne.fl(jl)) THEN
          xh(jl)=xnew(jl)
          fh(jl)=fnew(jl)
        ELSE IF (sign(fh(jl),fnew(jl)).ne.fh(jl)) THEN
          xl(jl)=xnew(jl)
          fl(jl)=fnew(jl)
        ENDIF
      ENDIF ! liter
    ENDDO !jl
    icount = count(liter(1:ictit))
    if( icount == 0 ) exit
  ENDDO !niter
!
!   Set a and diagnose b
!
!DIR$ CONCURRENT
  DO jl=1,ictit
    zvartarget=MAX(cqtmin,cvarmin*pqm1_l1(jl))
    zss=zbetass_l1(jl)
    zqt=zbetaqt_l1(jl)
    zpp=cbeta_pq
    zqq=zskew1(jl)
    zaa=MAX(xnew(jl),0.0_dp)
    zaa=MIN(zaa,zqt-zvartarget*zpp/(zpp+zqq))
    zbb=(zqt-zaa)*(zpp+zqq)/zpp+zaa
    zbb=MAX(zbb,zaa+zvartarget)
    pbetaa(iptit(jl),jk1(jl))=zaa
    pbetab(iptit(jl),jk1(jl))=zbb
  ENDDO
!
!   Overcast or clear sky: Diagnostic gridpoints
!
!DIR$ CONCURRENT
  DO jl=1,ictdg
    zvartarget=MAX(cqtmin,cvarmin*pqm1_l2(jl))
    zvar=MAX(pxvar(iptdg(jl),jk2(jl)),zvartarget) ! b-a width
!
!   Defined to make equations easier to understand:
!
    zss=zbetass_l2(jl)  ! qsat
    zqt=zbetaqt_l2(jl)  ! qtotal
    zpp=cbeta_pq        ! p shape factor
    zqq=zskew2(jl)      ! q shape factor
!
!   Limit a>0 and ensure (a,b)<qs or (a,b)>qs
!
    zaa=MAX(0._dp,zqt-zvar*zpp/(zpp+zqq)) !a from eqn(12) limited
    zbb=(zqt-zaa)*(zpp+zqq)/zpp+zaa    !b from eqn(12) again
    IF (zaa<zss .AND. zbb>zss) THEN
      IF (zqt<zss) THEN ! set b=qsat
        zbb = zss
        zaa = MAX(0.0_dp,(zqt*(zpp+zqq)-zbb*zpp)/zqq)
      ELSE ! set a=qsat
        zaa = zss
        zbb=(zqt-zaa)*(zpp+zqq)/zpp+zaa ! b from (12) again
      ENDIF
    ENDIF
    zbb=MAX(zbb,zaa+zvartarget)
    pbetaa(iptdg(jl),jk2(jl))=zaa
    pbetab(iptdg(jl),jk2(jl))=zbb
  ENDDO
!
  DO jk=ktdia,klev
    IF (lcover .AND. jk >= ncctop) THEN  ! beta function scheme
      DO jl=1,kproma
!
!   Define variance
!
        pxvar(jl,jk) = pbetab(jl,jk)-pbetaa(jl,jk)
!
!       set cloud fraction
!
        IF (pxim1(jl,jk)<=ccwmin .AND. pxlm1(jl,jk)<=ccwmin) THEN
          paclc(jl,jk)=0.0_dp
        ELSEIF( pqm1(jl,jk) >= pbetass(jl,jk)) THEN
          paclc(jl,jk)=1.0_dp
        ELSE
          ztt=(pbetass(jl,jk)-pbetaa(jl,jk))/                          &
               (pbetab(jl,jk)-pbetaa(jl,jk))
          ztt=nbetax*MAX(MIN(ztt,1.0_dp),0.0_dp)
          ixidx=INT(ztt)
          IF (ixidx==nbetax) THEN
            zbetai0=1.0_dp
          ELSE
            zbetai0=(ztt-ixidx)*tbetai0(iqidx(jl,jk),ixidx+1)          &
                 +(ixidx+1._dp-ztt)*tbetai0(iqidx(jl,jk),ixidx)
          ENDIF
          paclc(jl,jk) = 1.0_dp - zbetai0
          paclc(jl,jk) = MAX(zbetacl(jl,jk)/cmmrmax,paclc(jl,jk))
!    Fractional cloud cover > 0.01 required for middle atmosphere
          paclc(jl,jk) = MAX(MIN(paclc(jl,jk),1.0_dp),0.01_dp)
        ENDIF
      ENDDO !jl loop
    ELSE  ! lcover=.false.
!
!       1.   Calculate the saturation mixing ratio
!
      DO jl=1,kproma
        lo2=(ptm1(jl,jk).LT.cthomi).OR.(ptm1(jl,jk).LT.tmelt           &
             .AND.pxim1(jl,jk).GT.csecfrl)
        it=NINT(ptm1(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsm1(jl)=MERGE(tlucua(it),tlucuaw(it),lo2)/papm1(jl,jk)
        zqsm1(jl)=MIN(zqsm1(jl),0.5_dp)
        zcor=1._dp/(1._dp-vtmpc1*zqsm1(jl))
        zqsm1(jl)=zqsm1(jl)*zcor
      ENDDO

      IF (lookupoverflow) CALL lookuperror ('cover       ')
!
!       Threshold relative humidity, qsat and cloud cover
!       This is from cloud, and is the original calculation for
!       cloud cover, based on relative humidity
!       (Lohmann and Roeckner, Clim. Dyn.  96)
!
      DO jl=1,kproma
!
        zrhc=crt+(crs-crt)*EXP(1._dp-(paphm1(jl,klevp1)                &
             /papm1(jl,jk))**nex)
        zsat=1._dp
        jb=knvb(jl)
        lo=(jb.GE.jbmin .AND. jb.LE.jbmax .AND. pvervel(jl,jk).GT.0._dp)
        lo1=(jk.EQ.jb .OR. jk.EQ.jb+1)
        IF (lo .AND. lo1) THEN
          printop(jl)=REAL(jb,dp)
          zsat=csatsc*zrhc
          zrhc=crhsc*zrhc
        END IF
        zqr=zqsm1(jl)*zsat*zrhc
        paclc(jl,jk)=(pqm1(jl,jk)-zqr)/(zqsm1(jl)*zsat-zqr)
        paclc(jl,jk)=MAX(paclc(jl,jk),0.0_dp)
        paclc(jl,jk)=MIN(paclc(jl,jk),1.0_dp)
        paclc(jl,jk)=1._dp-SQRT(1._dp-paclc(jl,jk))
        IF (pxim1(jl,jk)<=ccwmin .AND. pxlm1(jl,jk)<=ccwmin) THEN
          paclc(jl,jk)=0.0_dp
        ENDIF
        paclc(jl,jk) = MAX(MIN(paclc(jl,jk),1.0_dp),0.0_dp)
      END DO !jl
    ENDIF !lcover
  END DO  !jk

!
END SUBROUTINE cover
