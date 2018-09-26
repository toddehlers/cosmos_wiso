SUBROUTINE cududv(   kproma,   kbdim,    klev,     klevp1,             &
           ktopm2,   ktype,    kcbot,    paphp1,   ldcum,              &
           puen,     pven,     pvom,     pvol,                         &
           puu,      pud,      pvu,      pvd,                          &
           pmfu,     pmfd)
!
!
!**** *CUDUDV* - UPDATES U AND V TENDENCIES,
!                DOES GLOBAL DIAGNOSTIC OF DISSIPATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!**   INTERFACE.
!     ----------
!
!          *CUDUDV* IS CALLED FROM *CUMASTR*
!
USE mo_kind,       ONLY: dp
USE mo_constants,  ONLY: g
!
IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktopm2
!
REAL(dp):: puen(kbdim,klev),        pven(kbdim,klev),                  &
           pvol(kbdim,klev),        pvom(kbdim,klev),                  &
           paphp1(kbdim,klevp1)
REAL(dp):: puu(kbdim,klev),         pud(kbdim,klev),                   &
           pvu(kbdim,klev),         pvd(kbdim,klev),                   &
           pmfu(kbdim,klev),        pmfd(kbdim,klev)
INTEGER :: ktype(kbdim),            kcbot(kbdim)
LOGICAL :: ldcum(kbdim)
!
REAL(dp):: zmfuu(kbdim,klev),       zmfdu(kbdim,klev),                 &
           zmfuv(kbdim,klev),       zmfdv(kbdim,klev)
!
INTEGER :: jl, jk, ik, ikb
REAL(dp):: zzp, zdudt, zdvdt
!
!----------------------------------------------------------------------
!
!*    1.0          CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
!                  ----------------------------------------------
!
100 CONTINUE
  IF(ktopm2.EQ.1) THEN
    DO 107 jk=2,klev
       ik=jk-1
       DO 106 jl=1,kproma
          IF(ldcum(jl)) THEN
             zmfuu(jl,jk)=pmfu(jl,jk)*(puu(jl,jk)-puen(jl,ik))
             zmfuv(jl,jk)=pmfu(jl,jk)*(pvu(jl,jk)-pven(jl,ik))
             zmfdu(jl,jk)=pmfd(jl,jk)*(pud(jl,jk)-puen(jl,ik))
             zmfdv(jl,jk)=pmfd(jl,jk)*(pvd(jl,jk)-pven(jl,ik))
          END IF
106    END DO
107  END DO
    DO 105 jl=1,kproma
      IF(ldcum(jl)) THEN
        zmfuu(jl,1)=zmfuu(jl,2)
        zmfuv(jl,1)=zmfuv(jl,2)
        zmfdu(jl,1)=zmfdu(jl,2)
        zmfdv(jl,1)=zmfdv(jl,2)
      END IF
105 END DO
  ELSE
    DO 120 jk=ktopm2,klev
       ik=jk-1
       DO 110 jl=1,kproma
          IF(ldcum(jl)) THEN
             zmfuu(jl,jk)=pmfu(jl,jk)*(puu(jl,jk)-puen(jl,ik))
             zmfuv(jl,jk)=pmfu(jl,jk)*(pvu(jl,jk)-pven(jl,ik))
             zmfdu(jl,jk)=pmfd(jl,jk)*(pud(jl,jk)-puen(jl,ik))
             zmfdv(jl,jk)=pmfd(jl,jk)*(pvd(jl,jk)-pven(jl,ik))
          END IF
110    END DO
120  END DO
  END IF
  DO 140 jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
     DO 130 jl=1,kproma
        IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
           ikb=kcbot(jl)
           zzp=((paphp1(jl,klevp1)-paphp1(jl,jk))/                     &
                         (paphp1(jl,klevp1)-paphp1(jl,ikb)))
           zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
           zmfuu(jl,jk)=zmfuu(jl,ikb)*zzp
           zmfuv(jl,jk)=zmfuv(jl,ikb)*zzp
           zmfdu(jl,jk)=zmfdu(jl,ikb)*zzp
           zmfdv(jl,jk)=zmfdv(jl,ikb)*zzp
        END IF
130  END DO
140 END DO
!
  DO 190 jk=ktopm2,klev
!
     IF(jk.LT.klev) THEN
        DO 160 jl=1,kproma
           IF(ldcum(jl)) THEN
              zdudt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*               &
                          (zmfuu(jl,jk+1)-zmfuu(jl,jk)+                &
                           zmfdu(jl,jk+1)-zmfdu(jl,jk))
              zdvdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*               &
                          (zmfuv(jl,jk+1)-zmfuv(jl,jk)+                &
                           zmfdv(jl,jk+1)-zmfdv(jl,jk))
              pvom(jl,jk)=pvom(jl,jk)+zdudt
              pvol(jl,jk)=pvol(jl,jk)+zdvdt
           END IF
160     END DO
!
     ELSE
        DO 170 jl=1,kproma
           IF(ldcum(jl)) THEN
              zdudt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
                           (zmfuu(jl,jk)+zmfdu(jl,jk))
              zdvdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
                           (zmfuv(jl,jk)+zmfdv(jl,jk))
              pvom(jl,jk)=pvom(jl,jk)+zdudt
              pvol(jl,jk)=pvol(jl,jk)+zdvdt
           END IF
170     END DO
     END IF
!
190 END DO
!
!
  RETURN
END SUBROUTINE cududv
