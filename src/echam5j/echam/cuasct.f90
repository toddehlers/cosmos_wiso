SUBROUTINE cuasct(   kproma, kbdim, klev, klevp1, klevm1,              &
           ptenh,    pqenh,    puen,     pven,                         &
           ktrac,                                                      &
           pxtenh,   pxten,    pxtu,     pmfuxt,                       &
           pten,     pqen,     pqsen,                                  &
           pgeo,     pgeoh,    paphp1,                                 &
           pqte,     pverv,    klwmin,                                 &
           ldcum,    ldland,   ktype,    klab,                         &
           ptu,      pqu,      plu,      puu,      pvu,                &
           pmfu,     pmfub,    pentr,                                  &
           pmfus,    pmfuq,                                            &
           pmful,    plude,    pqude,    pdmfup,                       &
           pcpen,    pcpcu,                                            &
           kcbot,    kctop,    kctop0,                                 &
!---wiso-code
           lwiso, kwiso,                                               &
           pwisoqenh,                                                  &
           pwisoqen,                                                   &
           pwisoqu,  pwisolu,                                          &
           pwisomfuq,                                                  &
           pwisomful,pwisolude,pwisoqude,pwisodmfup                    &
           )
!---wiso-code-end
!
!          THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
!          FOR CUMULUS PARAMETERIZATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
!          (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
!           FLUXES AS WELL AS PRECIPITATION RATES)
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUMASTRT*.
!
!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!          AND THEN CALCULATE MOIST ASCENT FOR
!          ENTRAINING/DETRAINING PLUME.
!          ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
!          SHALLOW AND DEEP CUMULUS CONVECTION.
!          IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
!          CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
!          (CLOUD BASE VALUES CALCULATED IN *CUBASMC*)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
!          *CUENTR*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
!          *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!
!          REFERENCE
!          ---------
!          (TIEDTKE,1989)
!
USE mo_kind,         ONLY : dp
USE mo_control,      ONLY : nn
USE mo_constants,    ONLY : g, tmelt, vtmpc1
USE mo_cumulus_flux, ONLY : lmfdudv, lmfmid, nmctop, cmfcmin, cprcon   &
                          , cmfctop
USE mo_time_control, ONLY : time_step_len

!---wiso-code

USE mo_wiso,         ONLY : wiso_frac_liq_ice, twisoice, cwisomin, cwisosec

!---wiso-code-end

IMPLICIT NONE

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1, ktrac

!---wiso-code

LOGICAL, INTENT (IN) :: lwiso
INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

INTEGER :: jl, jk, jt, ik, icall 
!
REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pten(kbdim,klev),        pqen(kbdim,klev),                 &
            pgeo(kbdim,klev),        pgeoh(kbdim,klev),                &
            paphp1(kbdim,klevp1),                                      &
            pqsen(kbdim,klev),       pqte(kbdim,klev),                 &
            pverv(kbdim,klev)
!
REAL(dp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                  &
            puu(kbdim,klev),         pvu(kbdim,klev),                  &
            pmfu(kbdim,klev),                                          &
            pmfub(kbdim),            pentr(kbdim),                     &
            pmfus(kbdim,klev),       pmfuq(kbdim,klev),                &
            plu(kbdim,klev),         plude(kbdim,klev),                &
            pqude(kbdim,klev),                                         &
            pmful(kbdim,klev),       pdmfup(kbdim,klev)
REAL(dp) :: pcpen(kbdim,klev),       pcpcu(kbdim,klev)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                     &
            klab(kbdim,klev),        kcbot(kbdim),                     &
            kctop(kbdim),            kctop0(kbdim)
LOGICAL  :: ldcum(kbdim),            ldland(kbdim)
!
REAL(dp) :: zdmfen(kbdim),           zdmfde(kbdim),                    &
            zmfuu(kbdim),            zmfuv(kbdim),                     &
            zpbase(kbdim),           zqold(kbdim)
REAL(dp) :: zph(kbdim)
LOGICAL  :: loflag(kbdim)
REAL(dp) :: pxtenh(kbdim,klev,ktrac),pxten(kbdim,klev,ktrac),          &
            pxtu(kbdim,klev,ktrac),  pmfuxt(kbdim,klev,ktrac)
!
!---wiso-code

REAL(dp), OPTIONAL ::  pwisoqenh(kbdim,klev,kwiso),                               &
                       pwisoqen(kbdim,klev,kwiso)

REAL(dp), OPTIONAL  :: pwisoqu(kbdim,klev,kwiso),                                 &
                       pwisomfuq(kbdim,klev,kwiso),                               &
                       pwisolu(kbdim,klev,kwiso),pwisolude(kbdim,klev,kwiso),     &
                       pwisoqude(kbdim,klev,kwiso),                               &
                       pwisomful(kbdim,klev,kwiso),pwisodmfup(kbdim,klev,kwiso)

!---wiso-code-end

REAL(dp) :: zcons2, ztglace, zmfmax, zfac, zmftest, zqeen, zseen       &
          , zscde, zqude, zmfusk, zmfuqk, zmfulk, zxteen, zxtude       &
          , zmfuxtk, zbuo, zdnoprc, zprcon, zlnew, zz, zdmfeu, zdmfdu  &
          , zzdmf, zdlev
!
!---wiso-code

REAL(dp) :: zprec_tmp(kbdim,klev)

REAL(dp) :: zwisofracliq(kbdim,kwiso), zwisofracice(kbdim,kwiso)

REAL(dp) :: zwisoqeen,                                                 &
            zwisoqude, zwisomfuqk, zwisomfulk

REAL(dp) :: zqliq, zqice,                                              &
            zql,  zqv,   zwisoql, zwisoqv,                             &
            zquot,zquot2,zqvo,    zdelta

LOGICAL  :: lo
            
!---wiso-code-end

!      INTRINSIC FUNCTIONS
INTRINSIC MAX, MIN
!
!
!----------------------------------------------------------------------
!
!*    1.           SPECIFY PARAMETERS
!                  ------------------
!
100 CONTINUE
  zcons2=1._dp/(g*time_step_len)
  ztglace=tmelt-13._dp
  IF(klev == 11) THEN
    IF(nn == 21) THEN
      zdlev=1.5E4_dp
    ELSE IF(nn == 31) THEN
      zdlev=2.0E4_dp
    ELSE
      zdlev=3.0E4_dp
    ENDIF
  ELSE
   zdlev=3.0E4_dp
  ENDIF
!
!
!----------------------------------------------------------------------
!
!     2.           SET DEFAULT VALUES
!                  ------------------
!
200 CONTINUE
  DO 210 jl=1,kproma
     zmfuu(jl)=0._dp
     zmfuv(jl)=0._dp
     IF(.NOT.ldcum(jl)) ktype(jl)=0
210 END DO
  DO 230 jk=1,klev
     DO 220 jl=1,kproma
        plu(jl,jk)=0._dp
        pmfu(jl,jk)=0._dp
        pmfus(jl,jk)=0._dp
        pmfuq(jl,jk)=0._dp
        pmful(jl,jk)=0._dp
        plude(jl,jk)=0._dp
        pqude(jl,jk)=0._dp
        pdmfup(jl,jk)=0._dp
        IF(.NOT.ldcum(jl).OR.ktype(jl).EQ.3) klab(jl,jk)=0
        IF(.NOT.ldcum(jl).AND.paphp1(jl,jk).LT.4.e4_dp) kctop0(jl)=jk
        IF(jk.LT.kcbot(jl)) klab(jl,jk)=0
220  END DO

!---wiso-code
  IF (lwiso) THEN
  
     DO jt=1,kwiso
       DO jl=1,kproma
         pwisolu(jl,jk,jt)=0._dp
         pwisomfuq(jl,jk,jt)=0._dp
         pwisomful(jl,jk,jt)=0._dp
         pwisolude(jl,jk,jt)=0._dp
         pwisoqude(jl,jk,jt)=0._dp
         pwisodmfup(jl,jk,jt)=0._dp
       END DO
     END DO

  END IF
!---wiso-code-end

     DO 2204 jt=1,ktrac
        DO 2202 jl=1,kproma
           pmfuxt(jl,jk,jt)=0._dp
2202    END DO
2204 END DO
!
230 END DO
!
!
!----------------------------------------------------------------------
!
!     3.0          INITIALIZE VALUES AT LIFTING LEVEL
!                  ----------------------------------
!
300 CONTINUE
  DO 310 jl=1,kproma
     kctop(jl)=klevm1
     IF(.NOT.ldcum(jl)) THEN
        kcbot(jl)=klevm1
        pmfub(jl)=0._dp
        pqu(jl,klev)=0._dp
     END IF
     pmfu(jl,klev)=pmfub(jl)
     pmfus(jl,klev)=pmfub(jl)*(pcpcu(jl,klev)*ptu(jl,klev)             &
                                       +pgeoh(jl,klev))
     pmfuq(jl,klev)=pmfub(jl)*pqu(jl,klev)
     IF(lmfdudv) THEN
        zmfuu(jl)=pmfub(jl)*puu(jl,klev)
        zmfuv(jl)=pmfub(jl)*pvu(jl,klev)
     END IF
310 END DO
!
!---wiso-code
  IF (lwiso) THEN

  DO jt=1,kwiso
    DO jl=1,kproma
      IF(.NOT.ldcum(jl)) THEN
        pwisoqu(jl,klev,jt)=0._dp
      ENDIF
     pwisomfuq(jl,klev,jt)=pmfub(jl)*pwisoqu(jl,klev,jt)
    END DO
  END DO

  END IF
!---wiso-code-end

  DO 3112 jt=1,ktrac
     DO 3110 jl=1,kproma
        IF(.NOT.ldcum(jl)) THEN
           pxtu(jl,klev,jt)=0._dp
        ENDIF
        pmfuxt(jl,klev,jt)=pmfub(jl)*pxtu(jl,klev,jt)
3110 END DO
3112 END DO
!
  DO 320 jl=1,kproma
     ldcum(jl)=.FALSE.
320 END DO
!
!
!----------------------------------------------------------------------
!
!     4.           DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
!                  BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
!                  BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!                  THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
!                  -------------------------------------------------
!
400 CONTINUE
  DO 480 jk=klevm1,2,-1
!
!                  SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!                  IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
!                  ----------------------------------------------------
!
     ik=jk
     IF(lmfmid.AND.ik.LT.klevm1.AND.ik.GT.nmctop) THEN

       IF (lwiso) THEN

        CALL cubasmc(kproma,   kbdim,    klev,     ik,       klab,     &
                     pten,     pqen,     pqsen,    puen,     pven,     &
                     ktrac,                                            &
                     pxten,    pxtu,     pmfuxt,                       &
                     pverv,    pgeo,     pgeoh,    ldcum,    ktype,    &
                     pmfu,     pmfub,    pentr,    kcbot,              &
                     ptu,      pqu,      plu,      puu,      pvu,      &
                     pmfus,    pmfuq,    pmful,    pdmfup,   zmfuu,    &
                     pcpen,                                            &
                     zmfuv,                                            &
!---wiso-code
                     lwiso, kwiso,                                     &
                     pwisoqen,                                         &
                     pwisoqu,  pwisolu,                                &
                     pwisomfuq,pwisomful,pwisodmfup)
!---wiso-code-end

       ELSE
       
        CALL cubasmc(kproma,   kbdim,    klev,     ik,       klab,     &
                     pten,     pqen,     pqsen,    puen,     pven,     &
                     ktrac,                                            &
                     pxten,    pxtu,     pmfuxt,                       &
                     pverv,    pgeo,     pgeoh,    ldcum,    ktype,    &
                     pmfu,     pmfub,    pentr,    kcbot,              &
                     ptu,      pqu,      plu,      puu,      pvu,      &
                     pmfus,    pmfuq,    pmful,    pdmfup,   zmfuu,    &
                     pcpen,                                            &
                     zmfuv,                                            &
!---wiso-code
                     lwiso, kwiso)
!---wiso-code-end
       
       END IF
       
     ENDIF
!
     DO 410 jl=1,kproma
        IF(klab(jl,jk+1).EQ.0) klab(jl,jk)=0
        loflag(jl)=klab(jl,jk+1).GT.0
        zph(jl)=paphp1(jl,jk)
        IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
           zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
           IF(pmfub(jl).GT.zmfmax) THEN
              zfac=zmfmax/pmfub(jl)
              pmfu(jl,jk+1)=pmfu(jl,jk+1)*zfac
              pmfus(jl,jk+1)=pmfus(jl,jk+1)*zfac
              pmfuq(jl,jk+1)=pmfuq(jl,jk+1)*zfac
              zmfuu(jl)=zmfuu(jl)*zfac
              zmfuv(jl)=zmfuv(jl)*zfac
           END IF
        END IF
410  END DO

!---wiso-code
  IF (lwiso) THEN

     DO jt=1,kwiso
       DO jl=1,kproma
         IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
           zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
           IF(pmfub(jl).GT.zmfmax) THEN
             zfac=zmfmax/pmfub(jl)
             pwisomfuq(jl,jk+1,jt)=pwisomfuq(jl,jk+1,jt)*zfac
           END IF
         END IF
       END DO
     END DO

  END IF
!---wiso-code-end

     DO 4102 jt=1,ktrac
        DO 4101 jl=1,kproma
           IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
              zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
              IF(pmfub(jl).GT.zmfmax) THEN
                 zfac=zmfmax/pmfub(jl)
                 pmfuxt(jl,jk+1,jt)=pmfuxt(jl,jk+1,jt)*zfac
              END IF
           END IF
4101    END DO
4102 END DO
!
! RESET PMFUB IF NECESSARY
!
     DO 4103 jl=1,kproma
        IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
           zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
           pmfub(jl)=MIN(pmfub(jl),zmfmax)
        END IF
4103 END DO
!
!
!*                 SPECIFY ENTRAINMENT RATES IN *CUENTRT*
!                  --------------------------------------
!
     ik=jk
     CALL cuentrt(kproma,  kbdim,    klev,     klevp1,   ik,           &
                 ptenh,    pqenh,    pqte,     paphp1,                 &
                 klwmin,   ldcum,    ktype,    kcbot,    kctop0,       &
                 zpbase,   pmfu,     pentr,                            &
                 zdmfen,   zdmfde)
!
!
!
!                  DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
!                  ---------------------------------------------------
!
     DO 420 jl=1,kproma
        IF(loflag(jl)) THEN
           IF(jk.LT.kcbot(jl)) THEN
              zmftest=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)
              zmfmax=MIN(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2)
              zdmfen(jl)=MAX(zdmfen(jl)-MAX(zmftest-zmfmax,0._dp),0._dp)
           END IF
           zdmfde(jl)=MIN(zdmfde(jl),0.75_dp*pmfu(jl,jk+1))
           pmfu(jl,jk)=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)
           zqeen=pqenh(jl,jk+1)*zdmfen(jl)
           zseen=(pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))        &
                                               *zdmfen(jl)
           zscde=(pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1))          &
                                               *zdmfde(jl)
           zqude=pqu(jl,jk+1)*zdmfde(jl)
           pqude(jl,jk)=zqude
           plude(jl,jk)=plu(jl,jk+1)*zdmfde(jl)
           zmfusk=pmfus(jl,jk+1)+zseen-zscde
           zmfuqk=pmfuq(jl,jk+1)+zqeen-zqude
           zmfulk=pmful(jl,jk+1)    -plude(jl,jk)
           plu(jl,jk)=zmfulk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
           pqu(jl,jk)=zmfuqk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
           ptu(jl,jk)=(zmfusk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))-        &
                                pgeoh(jl,jk))/pcpcu(jl,jk)
           ptu(jl,jk)=MAX(100._dp,ptu(jl,jk))
           ptu(jl,jk)=MIN(400._dp,ptu(jl,jk))
           zqold(jl)=pqu(jl,jk)
        END IF
420  END DO
!
!
!---wiso-code
  IF (lwiso) THEN

     DO jt=1,kwiso
       DO jl=1,kproma
         IF(loflag(jl)) THEN
           zwisoqeen=pwisoqenh(jl,jk+1,jt)*zdmfen(jl)
           zwisoqude=pwisoqu(jl,jk+1,jt)*zdmfde(jl)
           pwisoqude(jl,jk,jt)=zwisoqude
           pwisolude(jl,jk,jt)=pwisolu(jl,jk+1,jt)*zdmfde(jl)
           zwisomfuqk=pwisomfuq(jl,jk+1,jt)+zwisoqeen-zwisoqude
           zwisomfulk=pwisomful(jl,jk+1,jt)-pwisolude(jl,jk,jt)
           pwisolu(jl,jk,jt)=zwisomfulk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
           pwisoqu(jl,jk,jt)=zwisomfuqk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
         ENDIF
       END DO
     END DO

  END IF
!---wiso-code-end

     DO 4204 jt=1,ktrac
        DO 4202 jl=1,kproma
           IF(loflag(jl)) THEN
              zxteen=pxtenh(jl,jk+1,jt)*zdmfen(jl)
              zxtude=pxtu(jl,jk+1,jt)*zdmfde(jl)
              zmfuxtk=pmfuxt(jl,jk+1,jt)+zxteen-zxtude
              pxtu(jl,jk,jt)=zmfuxtk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
           ENDIF
4202    END DO
4204 END DO
!
!
!                  DO CORRECTIONS FOR MOIST ASCENT
!                  BY ADJUSTING T,Q AND L IN *CUADJTQ*
!                  -----------------------------------
!
     ik=jk
     icall=1
     CALL cuadjtq(kproma,   kbdim,    klev,     ik,                    &
                  zph,      ptu,      pqu,      loflag,   icall)
!

!---wiso-code
  IF (lwiso) THEN

! calculate fractionation factors for liquid-vapour and solid-vapour phase change

     CALL wiso_frac_liq_ice(kproma,kbdim,kwiso,ptu(:,ik),zwisofracliq,zwisofracice)

!    adjusting isotope vapour *pwisoqu* depending on fract.
     DO jt=1,kwiso
       DO jl=1,kproma
         IF (loflag(jl)) THEN
! divide into ice and liquid
           IF (ptu(jl,jk).gt.tmelt) THEN
             zqliq=1._dp
             zqice=0._dp
           ELSEIF (ptu(jl,jk).lt.twisoice) THEN
             zqice=1._dp
             zqliq=0._dp
           ELSE
             zqice=(tmelt-ptu(jl,jk))/(tmelt-twisoice)
             zqliq=1._dp-zqice
           ENDIF
! fractionation of ice and liquid
           zql=zqliq*(plu(jl,jk)+zqold(jl)-pqu(jl,jk))
           zqv=zqold(jl)+zqliq*(pqu(jl,jk)-zqold(jl))
           zwisoql=zqliq*pwisolu(jl,jk,jt)
           zwisoqv=pwisoqu(jl,jk,jt)
           zquot=zqv+zwisofracliq(jl,jt)*zql
           zqvo=zqv
           IF (zquot.gt.cwisomin.and.zqvo.gt.cwisomin) THEN
             zdelta=(zwisofracliq(jl,jt)*zql*zwisoqv-zwisoql*zqv)/zquot 
             zqv=pqu(jl,jk)
             zwisoqv=zwisoqv-zdelta
             zquot2=zqv/zqvo
             lo=abs(1._dp-zquot2).lt.cwisosec
             zquot2=MERGE(1._dp,zquot2,lo)
             zdelta=zdelta+zwisoqv*(1._dp-zquot2**zwisofracice(jl,jt))
             pwisoqu(jl,jk,jt)=pwisoqu(jl,jk,jt)-zdelta
             IF(loflag(jl).AND.pqu(jl,jk).LT.zqold(jl)) THEN  ! add condensate to liquid as it is done
               pwisolu(jl,jk,jt)=pwisolu(jl,jk,jt)+zdelta     ! for normal water in next loop 440
             END IF
           ENDIF
         ENDIF
       END DO
     END DO

  END IF
!---wiso-code-end

!DIR$ IVDEP
!OCL NOVREC
     DO 440 jl=1,kproma
        IF(loflag(jl).AND.pqu(jl,jk).LT.zqold(jl)) THEN
           klab(jl,jk)=2
           plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
           zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk)-plu(jl,jk))-       &
                    ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))
           IF(klab(jl,jk+1).EQ.1) zbuo=zbuo+0.5_dp
           IF(zbuo.GT.0._dp.AND.pmfu(jl,jk).GE.0.1_dp*pmfub(jl)) THEN
              kctop(jl)=jk
              ldcum(jl)=.TRUE.
              zdnoprc=MERGE(zdlev,1.5e4_dp,ldland(jl))
              zprcon=MERGE(0._dp,cprcon,                               &
                               zpbase(jl)-paphp1(jl,jk).LT.zdnoprc)
              zlnew=plu(jl,jk)/                                        &
                           (1._dp+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
              pdmfup(jl,jk)=MAX(0._dp,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
!---wiso-code
! Store fraction of precipitation formed from cloud water
              IF(plu(jl,jk).gt.cwisomin) THEN
                zprec_tmp(jl,jk)=(plu(jl,jk)-zlnew)/plu(jl,jk)
              ENDIF
              IF (abs(1._dp-zprec_tmp(jl,jk)).lt.cwisosec) zprec_tmp(jl,jk)=1._dp
!---wiso-code-end
              plu(jl,jk)=zlnew
           ELSE
              klab(jl,jk)=0
              pmfu(jl,jk)=0._dp
!---wiso-code
! No precipitation formed from cloud water
              zprec_tmp(jl,jk)=0._dp
!---wiso-code-end
           END IF
        END IF
440  END DO

!---wiso-code
  IF (lwiso) THEN

     DO jt=1,kwiso
       DO jl=1,kproma
        IF(loflag(jl).AND.pqu(jl,jk).LT.zqold(jl)) THEN
! Calculate water isotope precipitation and new mass fluxes for the water isotopes in updrafts
           pwisodmfup(jl,jk,jt)=MAX(0._dp,pwisolu(jl,jk,jt)*pmfu(jl,jk)*zprec_tmp(jl,jk))
           pwisolu(jl,jk,jt)=(1._dp-zprec_tmp(jl,jk))*pwisolu(jl,jk,jt)
         END IF
       END DO
     END DO

  END IF
!---wiso-code-end

     DO 455 jl=1,kproma
        IF(loflag(jl)) THEN
           pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
           pmfus(jl,jk)=(pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))         &
                                 *pmfu(jl,jk)
           pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
        END IF
455  END DO

!---wiso-code
  IF (lwiso) THEN

     DO jt=1,kwiso
       DO jl=1,kproma
         IF(loflag(jl)) THEN
           pwisomful(jl,jk,jt)=pwisolu(jl,jk,jt)*pmfu(jl,jk)
           pwisomfuq(jl,jk,jt)=pwisoqu(jl,jk,jt)*pmfu(jl,jk)
         END IF
       END DO
     END DO

  END IF
!---wiso-code-end

     DO 4554 jt=1,ktrac
        DO 4552 jl=1,kproma
           IF(loflag(jl)) THEN
              pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)
           ENDIF
4552    END DO
4554 END DO
!
     IF(lmfdudv) THEN
        DO 460 jl=1,kproma
           IF(loflag(jl)) THEN
              IF(ktype(jl).EQ.1.OR.ktype(jl).EQ.3) THEN
                 zz=MERGE(3._dp,2._dp,zdmfen(jl).EQ.0._dp)
              ELSE
                 zz=MERGE(1._dp,0._dp,zdmfen(jl).EQ.0._dp)
              END IF
              zdmfeu=zdmfen(jl)+zz*zdmfde(jl)
              zdmfdu=zdmfde(jl)+zz*zdmfde(jl)
              zdmfdu=MIN(zdmfdu,0.75_dp*pmfu(jl,jk+1))
              zmfuu(jl)=zmfuu(jl)+                                     &
                             zdmfeu*puen(jl,jk)-zdmfdu*puu(jl,jk+1)
              zmfuv(jl)=zmfuv(jl)+                                     &
                             zdmfeu*pven(jl,jk)-zdmfdu*pvu(jl,jk+1)
              IF(pmfu(jl,jk).GT.0._dp) THEN
                 puu(jl,jk)=zmfuu(jl)*(1._dp/pmfu(jl,jk))
                 pvu(jl,jk)=zmfuv(jl)*(1._dp/pmfu(jl,jk))
              END IF
           END IF
460     END DO
     END IF
!
480 END DO
!
!
!----------------------------------------------------------------------
!
!     5.           DETERMINE CONVECTIVE FLUXES ABOVE NON-BUOYANCY LEVEL
!                  ----------------------------------------------------
!                  (NOTE: CLOUD VARIABLES LIKE T,Q AND L ARE NOT
!                         AFFECTED BY DETRAINMENT AND ARE ALREADY KNOWN
!                         FROM PREVIOUS CALCULATIONS ABOVE)
!
500 CONTINUE
  DO 510 jl=1,kproma
     IF(kctop(jl).EQ.klevm1) ldcum(jl)=.FALSE.
     kcbot(jl)=MAX(kcbot(jl),kctop(jl))
510 END DO
!DIR$ IVDEP
  DO 530 jl=1,kproma
     IF(ldcum(jl)) THEN
        jk=kctop(jl)-1
        zzdmf=cmfctop
        zdmfde(jl)=(1._dp-zzdmf)*pmfu(jl,jk+1)
        plude(jl,jk)=zdmfde(jl)*plu(jl,jk+1)
        pqude(jl,jk)=zdmfde(jl)*pqu(jl,jk+1)
        pmfu(jl,jk)=pmfu(jl,jk+1)-zdmfde(jl)
        pdmfup(jl,jk)=0._dp
        pmfus(jl,jk)=(pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
        pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
        pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
        plude(jl,jk-1)=pmful(jl,jk)
        pqude(jl,jk-1)=pmfuq(jl,jk)
     END IF
530 END DO

!---wiso-code
  IF (lwiso) THEN

! Calculate water isotope fluxes above the non-buoyancy level
     DO jt=1,kwiso
!DIR$ IVDEP
       DO jl=1,kproma
         IF (ldcum(jl)) THEN
           jk=kctop(jl)-1
           pwisolude(jl,jk,jt)=zdmfde(jl)*pwisolu(jl,jk+1,jt)
           pwisoqude(jl,jk,jt)=zdmfde(jl)*pwisoqu(jl,jk+1,jt)
           pwisodmfup(jl,jk,jt)=0._dp
           pwisomfuq(jl,jk,jt)=pwisoqu(jl,jk,jt)*pmfu(jl,jk)
           pwisomful(jl,jk,jt)=pwisolu(jl,jk,jt)*pmfu(jl,jk)
           pwisolude(jl,jk-1,jt)=pwisomful(jl,jk,jt)
           pwisoqude(jl,jk-1,jt)=pwisomfuq(jl,jk,jt)
         END IF
       END DO
     END DO

  END IF
!---wiso-code-end

  DO 5312 jt=1,ktrac
     DO 5310 jl=1,kproma
        IF(ldcum(jl)) THEN
           jk=kctop(jl)-1
           pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)
        ENDIF
5310 END DO
5312 END DO
!
  IF(lmfdudv) THEN
!DIR$      IVDEP
     DO 540 jl=1,kproma
        IF(ldcum(jl)) THEN
           jk=kctop(jl)-1
           puu(jl,jk)=puu(jl,jk+1)
           pvu(jl,jk)=pvu(jl,jk+1)
        END IF
540  END DO
  END IF
!
  RETURN
END SUBROUTINE cuasct
