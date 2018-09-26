SUBROUTINE cudtdq(kproma, kbdim, klev, klevp1, ktopm2, ldcum, ktrac,   &
                  paphp1,   pten,     ptte,     pqte,                  &
                  pxtte,    pxtec,    pmfuxt,   pmfdxt,                &
                  pmfus,    pmfds,    pmfuq,    pmfdq,                 &
                  pmful,    pdmfup,   pdmfdp,   plude,                 &
                  pdpmel,   prfl,     psfl,                            &
                  pcpen,    pqtec,    pqude,                           &
                  prsfc,    pssfc,    paprc,    paprs)
!
!
!**** *CUDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
!                DOES GLOBAL DIAGNOSTICS
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!**   INTERFACE.
!     ----------
!
!          *CUDTDQ* IS CALLED FROM *CUMASTR*
!
USE mo_kind,         ONLY : dp
USE mo_constants,    ONLY : alv, als, alf, tmelt, g
USE mo_tracer,       ONLY : trlist
USE mo_time_control, ONLY : delta_time
!
IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktopm2, ktrac
!
LOGICAL  llo1
!
REAL(dp) :: ptte(kbdim,klev),        pqte(kbdim,klev),                 &
            pten(kbdim,klev),        paphp1(kbdim,klevp1),             &
            paprc(kbdim),            paprs(kbdim),                     &
            prsfc(kbdim),            pssfc(kbdim)
REAL(dp) :: pmfus(kbdim,klev),       pmfds(kbdim,klev),                &
            pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                &
            pmful(kbdim,klev),       plude(kbdim,klev),                &
            pdmfup(kbdim,klev),      pdmfdp(kbdim,klev),               &
            pqtec(kbdim,klev),       pqude(kbdim,klev),                &
            pxtec(kbdim,klev),       prfl(kbdim)
REAL(dp) :: pdpmel(kbdim,klev),      psfl(kbdim)
REAL(dp) :: pcpen(kbdim,klev)
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: zmelt(kbdim)
REAL(dp) :: zsheat(kbdim)
REAL(dp) :: pxtte(kbdim,klev,ktrac), pmfuxt(kbdim,klev,ktrac),         &
            pmfdxt(kbdim,klev,ktrac)
!
REAL(dp) :: zrcpm ! reciprocal value of specific heat of moist air
INTEGER  :: jl, jk, jt
REAL(dp) :: zdiagt, zalv, zdtdt, zdqdt, zdxtdt
!
!----------------------------------------------------------------------
!
!*    1.0          SPECIFY PARAMETERS
!                  ------------------
!
100 CONTINUE
  zdiagt=delta_time
!
!
!----------------------------------------------------------------------
!
!*    2.0          INCREMENTATION OF T AND Q TENDENCIES
!                  ------------------------------------
!
200 CONTINUE
  DO 210 jl=1,kproma
     zmelt(jl)=0._dp
     zsheat(jl)=0._dp
210 END DO
!
  DO 250 jk=ktopm2,klev
!
     IF(jk.LT.klev) THEN
        DO 220 jl=1,kproma
           IF(ldcum(jl)) THEN
              llo1=(pten(jl,jk)-tmelt).GT.0._dp
              zalv=MERGE(alv,als,llo1)
              zrcpm=1._dp/pcpen(jl,jk)
              zdtdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*         &
                                  (pmfus(jl,jk+1)-pmfus(jl,jk)+        &
                                   pmfds(jl,jk+1)-pmfds(jl,jk)-        &
                                   alf*pdpmel(jl,jk)-                  &
                                   zalv*(pmful(jl,jk+1)-pmful(jl,jk)-  &
                                   plude(jl,jk)-                       &
                                  (pdmfup(jl,jk)+pdmfdp(jl,jk))))
              ptte(jl,jk)=ptte(jl,jk)+zdtdt
              zdqdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*               &
                                  (pmfuq(jl,jk+1)-pmfuq(jl,jk)+        &
                                   pmfdq(jl,jk+1)-pmfdq(jl,jk)+        &
                                   pmful(jl,jk+1)-pmful(jl,jk)-        &
                                   plude(jl,jk)-                       &
                                  (pdmfup(jl,jk)+pdmfdp(jl,jk)))
              pqte(jl,jk)=pqte(jl,jk)+zdqdt
              pxtec(jl,jk)=(g/(paphp1(jl,jk+1)-                        &
                            paphp1(jl,jk)))*plude(jl,jk)
              pqtec(jl,jk)=(g/(paphp1(jl,jk+1)-                        &
                            paphp1(jl,jk)))*pqude(jl,jk)
              zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
           END IF
220     END DO
!
        IF (trlist% anyconv /= 0) THEN
           DO 2204 jt=1,ktrac
              IF (trlist% ti(jt)% nconv == 1) THEN
                DO 2202 jl=1,kproma
                   IF(ldcum(jl)) THEN
                     zdxtdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))        &
                                 *(pmfuxt(jl,jk+1,jt)-pmfuxt(jl,jk,jt) &
                                  +pmfdxt(jl,jk+1,jt)-pmfdxt(jl,jk,jt))
                     pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
                   ENDIF
2202            END DO
              ENDIF
2204       END DO
        ENDIF
!
!
     ELSE
        DO 230 jl=1,kproma
           IF(ldcum(jl)) THEN
              llo1=(pten(jl,jk)-tmelt).GT.0._dp
              zalv=MERGE(alv,als,llo1)
              zrcpm=1._dp/pcpen(jl,jk)
              zdtdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*        &
                     (pmfus(jl,jk)+pmfds(jl,jk)+alf*pdpmel(jl,jk)-     &
                                 zalv*(pmful(jl,jk)+pdmfup(jl,jk)      &
                                +pdmfdp(jl,jk)+plude(jl,jk)))
              ptte(jl,jk)=ptte(jl,jk)+zdtdt
              zdqdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
                        (pmfuq(jl,jk)+pmfdq(jl,jk)+plude(jl,jk)+       &
                        (pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
              pqte(jl,jk)=pqte(jl,jk)+zdqdt
              pxtec(jl,jk)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))         &
                           *plude(jl,jk)
              pqtec(jl,jk)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))         &
                           *pqude(jl,jk)
              zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
           END IF
230     END DO
!
        IF (trlist% anyconv /= 0) THEN
           DO 2304 jt=1,ktrac
              IF (trlist% ti(jt)% nconv == 1) THEN
                DO 2302 jl=1,kproma
                   IF(ldcum(jl)) THEN
                      zdxtdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))      &
                             *(pmfuxt(jl,jk,jt)+pmfdxt(jl,jk,jt))
                      pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
                   ENDIF
2302            END DO
              END IF
2304       END DO
        ENDIF
!
     END IF
!
250 END DO
!
!
!---------------------------------------------------------------------
!
!      3.          UPDATE SURFACE FIELDS
!                  ---------------------
!
300 CONTINUE
  DO 310 jl=1,kproma
     prsfc(jl)=prfl(jl)
     pssfc(jl)=psfl(jl)
     paprc(jl)=paprc(jl)+zdiagt*(prfl(jl)+psfl(jl))
     paprs(jl)=paprs(jl)+zdiagt*psfl(jl)
310 END DO
!
  RETURN
END SUBROUTINE cudtdq
