SUBROUTINE cudlfs(   kproma, kbdim, klev, klevp1,                      &
           ptenh,    pqenh,    puen,     pven,                         &
           ktrac,                                                      &
           pxtenh,   pxtu,     pxtd,     pmfdxt,                       &
           pgeoh,    paphp1,                                           &
           ptu,      pqu,      puu,      pvu,                          &
           ldcum,    kcbot,    kctop,    pmfub,    prfl,               &
           ptd,      pqd,      pud,      pvd,                          &
           pmfd,     pmfds,    pmfdq,    pdmfdp,                       &
           pcpcu,                                                      &
           kdtop,    lddraf,                                           &
!---wiso-code
           lwiso, kwiso,                                               &
           pwisoqenh,                                                  &
           pwisoqu,                                                    &
           pwisoqd,                                                    &
           pwisomfdq,pwisodmfdp,                                       &
           pdmfup,   pwisodmfup                                        &
           )
!---wiso-code-end

!          THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
!          CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES
!
!          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
!          FOR MASSFLUX CUMULUS PARAMETERIZATION
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
!          AND UPDRAFT VALUES T,Q,U AND V AND ALSO
!          CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
!          IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.
!
!          METHOD.
!          --------
!          CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
!          MOIST ENVIRONMENTAL AIR AND CLOUD AIR.
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR CALCULATING WET BULB T AND Q AT LFS
!
USE mo_kind,          ONLY : dp
USE mo_constants,     ONLY : vtmpc1
USE mo_cumulus_flux,  ONLY : lmfdudv, lmfdd, cmfdeps
!
!---wiso-code

USE mo_constants,     ONLY : tmelt
USE mo_wiso,          ONLY : talphal1, talphal2, talphal3, cwisomin

!---wiso-code-end

IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kbdim, klev, ktrac, kproma, klevp1

!---wiso-code

LOGICAL, INTENT (IN) :: lwiso
INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pgeoh(kbdim,klev),       paphp1(kbdim,klevp1),             &
            ptu(kbdim,klev),         pqu(kbdim,klev),                  &
            puu(kbdim,klev),         pvu(kbdim,klev),                  &
            pmfub(kbdim),            prfl(kbdim)
!
REAL(dp) :: ptd(kbdim,klev),         pqd(kbdim,klev),                  &
            pud(kbdim,klev),         pvd(kbdim,klev),                  &
            pmfd(kbdim,klev),        pmfds(kbdim,klev),                &
            pmfdq(kbdim,klev),       pdmfdp(kbdim,klev)
REAL(dp) :: pcpcu(kbdim,klev)
INTEGER  :: kcbot(kbdim),            kctop(kbdim),                     &
            kdtop(kbdim)
LOGICAL  :: ldcum(kbdim),            lddraf(kbdim)
!
!---wiso-code

REAL(dp), OPTIONAL :: pwisoqenh(kbdim,klev,kwiso),                               &
                      pwisoqu(kbdim,klev,kwiso)

REAL(dp), OPTIONAL :: pwisoqd(kbdim,klev,kwiso),                                 &
                      pwisomfdq(kbdim,klev,kwiso),pwisodmfdp(kbdim,klev,kwiso)

REAL(dp), OPTIONAL :: pdmfup(kbdim,klev),      pwisodmfup(kbdim,klev,kwiso)

!---wiso-code-end

REAL(dp) :: ztenwb(kbdim,klev),      zqenwb(kbdim,klev),               &
            zcond(kbdim)
REAL(dp) :: zph(kbdim)
LOGICAL  :: llo2(kbdim)
REAL(dp) :: pxtenh(kbdim,klev,ktrac),pxtu(kbdim,klev,ktrac),           &
            pxtd(kbdim,klev,ktrac),  pmfdxt(kbdim,klev,ktrac)
LOGICAL  :: llo3(kbdim)
INTEGER  :: jl, jk, ke, is, ik, icall, jt
REAL(DP) :: zttest, zqtest, zbuo, zmftop
!
!---wiso-code

REAL(dp) :: zrain(kbdim),            zwisorain(kbdim,kwiso)

REAL(DP) :: zqv, zql, zwisoqv, zwisoql,                                &
            zt, zwisofracliq, zquot, zdelta

!---wiso-code-end

!
!----------------------------------------------------------------------
!
!     1.           SET DEFAULT VALUES FOR DOWNDRAFTS
!                  ---------------------------------
!
100 CONTINUE
  DO 110 jl=1,kproma
     lddraf(jl)=.FALSE.
     kdtop(jl)=klevp1
110 END DO
!
!---wiso-code
  IF (lwiso) THEN
  
    zrain(:) = 0._dp
    zwisorain(:,:) = 0._dp

  END IF
!---wiso-code-end

  IF(.NOT.lmfdd) go to 300
!
!
!----------------------------------------------------------------------
!
!     2.           DETERMINE LEVEL OF FREE SINKING BY
!                  DOING A SCAN FROM TOP TO BASE OF CUMULUS CLOUDS
!
!                  FOR EVERY POINT AND PROCEED AS FOLLOWS:
!
!                    (1) DETEMINE WET BULB ENVIRONMENTAL T AND Q
!                    (2) DO MIXING WITH CUMULUS CLOUD AIR
!                    (3) CHECK FOR NEGATIVE BUOYANCY
!
!                  THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
!                  OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
!                  TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
!                  EVAPORATION OF RAIN AND CLOUD WATER)
!                  ----------------------------------------------------
!
200 CONTINUE
!
  ke=klev-3
  DO 290 jk=3,ke
!
!
!     2.1          CALCULATE WET-BULB TEMPERATURE AND MOISTURE
!                  FOR ENVIRONMENTAL AIR IN *CUADJTQ*
!                  -------------------------------------------
!
!---wiso-code
  IF (lwiso) THEN

! Calculate Precipitation on this level
     DO jl=1,kproma
       zrain(jl)=zrain(jl)+pdmfup(jl,jk)
     END DO
     DO jt=1,kwiso
       DO jl=1,kproma
         zwisorain(jl,jt)=zwisorain(jl,jt)+pwisodmfup(jl,jk,jt)
       END DO
     END DO
     
  END IF
!---wiso-code-end

210  CONTINUE
     is=0
     DO 212 jl=1,kproma
        ztenwb(jl,jk)=ptenh(jl,jk)
        zqenwb(jl,jk)=pqenh(jl,jk)
        zph(jl)=paphp1(jl,jk)
        llo2(jl)=ldcum(jl).AND.prfl(jl).GT.0._dp.AND..NOT.lddraf(jl)   &
                          .AND.(jk.LT.kcbot(jl).AND.jk.GT.kctop(jl))
        is=is+MERGE(1,0,llo2(jl))
212  END DO
     IF(is.EQ.0) go to 290
!
     ik=jk
     icall=2
     CALL cuadjtq(kproma, kbdim, klev, ik,                             &
                  zph,      ztenwb,   zqenwb,   llo2,     icall)
!
!
!     2.2          DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
!                  AND CHECK FOR NEGATIVE BUOYANCY.
!                  THEN SET VALUES FOR DOWNDRAFT AT LFS.
!                  ----------------------------------------
!
220  CONTINUE
!DIR$ IVDEP
!OCL NOVREC
     DO 222 jl=1,kproma
        llo3(jl)=.FALSE.
        IF(llo2(jl)) THEN
           zttest=0.5_dp*(ptu(jl,jk)+ztenwb(jl,jk))
           zqtest=0.5_dp*(pqu(jl,jk)+zqenwb(jl,jk))
           zbuo=zttest*(1._dp+vtmpc1*zqtest)-                          &
                         ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))
           zcond(jl)=pqenh(jl,jk)-zqenwb(jl,jk)
           zmftop=-cmfdeps*pmfub(jl)
           IF(zbuo.LT.0._dp.AND.prfl(jl).GT.10._dp*zmftop*zcond(jl))   &
                                                                  THEN
              llo3(jl)=.TRUE.
              kdtop(jl)=jk
              lddraf(jl)=.TRUE.
              ptd(jl,jk)=zttest
              pqd(jl,jk)=zqtest
              pmfd(jl,jk)=zmftop
              pmfds(jl,jk)=pmfd(jl,jk)*(pcpcu(jl,jk)*ptd(jl,jk)        &
                                               +pgeoh(jl,jk))
              pmfdq(jl,jk)=pmfd(jl,jk)*pqd(jl,jk)
              pdmfdp(jl,jk-1)=-0.5_dp*pmfd(jl,jk)*zcond(jl)
              prfl(jl)=prfl(jl)+pdmfdp(jl,jk-1)
           END IF
        END IF
222  END DO
!
!---wiso-code
  IF (lwiso) THEN

     DO jl=1,kproma
       IF(llo3(jl)) THEN
         zrain(jl)=zrain(jl)+pdmfdp(jl,jk-1)
       ENDIF
     END DO

     DO jt=1,kwiso
       DO jl=1,kproma
         IF(llo3(jl)) THEN
          zwisoql=zwisorain(jl,jt)
          zwisoqv=pwisoqu(jl,jk,jt)+pwisoqenh(jl,jk,jt)
          zql=zrain(jl)
          zqv=pqd(jl,jk)
          zt=ptd(jl,jk)

! fractionation over water
          zwisofracliq=exp(talphal1(jt)/(zt**2._dp)+talphal2(jt)/zt+talphal3(jt))

          zquot=0.5_dp*(pmfd(jl,jk)*zqv-zwisofracliq*zql)
          IF (abs(zquot).lt.cwisomin) GOTO 2220

          zdelta=(zwisoql*zqv-0.5_dp*zwisofracliq*zql*zwisoqv)/zquot

          pwisoqd(jl,jk,jt)=0.5_dp*(zwisoqv-zdelta)
          pwisomfdq(jl,jk,jt)=pmfd(jl,jk)*pwisoqd(jl,jk,jt)
          pwisodmfdp(jl,jk-1,jt)=-0.5_dp*pmfd(jl,jk)*zdelta

        ENDIF
 2220   CONTINUE
 
      END DO
    END DO

  END IF
!---wiso-code-end

     DO 2224 jt=1,ktrac
        DO 2222 jl=1,kproma
           IF(llo3(jl)) THEN
              pxtd(jl,jk,jt)=0.5_dp*(pxtu(jl,jk,jt)+pxtenh(jl,jk,jt))
              pmfdxt(jl,jk,jt)=pmfd(jl,jk)*pxtd(jl,jk,jt)
           ENDIF
2222    END DO
2224 END DO
!
     IF(lmfdudv) THEN
        DO 224 jl=1,kproma
           IF(pmfd(jl,jk).LT.0._dp) THEN
              pud(jl,jk)=0.5_dp*(puu(jl,jk)+puen(jl,jk-1))
              pvd(jl,jk)=0.5_dp*(pvu(jl,jk)+pven(jl,jk-1))
           END IF
224     END DO
     END IF
!
290 END DO
!
300 CONTINUE
  RETURN
END SUBROUTINE cudlfs
