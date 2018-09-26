SUBROUTINE cubasmc(  kproma, kbdim, klev, kk, klab,                    &
           pten,     pqen,     pqsen,    puen,     pven,               &
           ktrac,                                                      &
           pxten,    pxtu,     pmfuxt,                                 &
           pverv,    pgeo,     pgeoh,    ldcum,    ktype,              &
           pmfu,     pmfub,    pentr,    kcbot,                        &
           ptu,      pqu,      plu,      puu,      pvu,                &
           pmfus,    pmfuq,    pmful,    pdmfup,   pmfuu,              &
           pcpen,                                                      &
           pmfuv,                                                      &
!---wiso-code
           lwiso, kwiso,                                               &
           pwisoqen,                                                   &
           pwisoqu,  pwisolu,                                          &
           pwisomfuq,pwisomful,pwisodmfup                              &
           )
!---wiso-code-end
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES CLOUD BASE VALUES
!          FOR MIDLEVEL CONVECTION
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          IT RETURNS CLOUDBASE VALUES FOR MIDLEVEL CONVECTION
!
!          METHOD.
!          --------
!          S. TIEDTKE (1989)
!
!          EXTERNALS
!          ---------
!          NONE
!
USE mo_kind,          ONLY : dp
USE mo_constants,     ONLY : g
USE mo_cumulus_flux,  ONLY : lmfdudv, entrmid, cmfcmin, cmfcmax
!
IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kbdim, klev, ktrac, kproma, kk

!---wiso-code

LOGICAL, INTENT (IN) :: lwiso
INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

REAL(dp) :: pten(kbdim,klev),        pqen(kbdim,klev),                 &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pqsen(kbdim,klev),       pverv(kbdim,klev),                &
            pgeo(kbdim,klev),        pgeoh(kbdim,klev)
!
REAL(dp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                  &
            puu(kbdim,klev),         pvu(kbdim,klev),                  &
            plu(kbdim,klev),         pmfu(kbdim,klev),                 &
            pmfub(kbdim),            pentr(kbdim),                     &
            pmfus(kbdim,klev),       pmfuq(kbdim,klev),                &
            pmful(kbdim,klev),       pdmfup(kbdim,klev),               &
            pmfuu(kbdim),            pmfuv(kbdim)
REAL(dp) :: pcpen(kbdim,klev)
INTEGER  :: ktype(kbdim),            kcbot(kbdim),                     &
            klab(kbdim,klev)
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: pxten(kbdim,klev,ktrac), pxtu(kbdim,klev,ktrac),           &
            pmfuxt(kbdim,klev,ktrac)

!---wiso-code

REAL(dp), OPTIONAL :: pwisoqen(kbdim,klev,kwiso)

REAL(dp), OPTIONAL :: pwisoqu(kbdim,klev,kwiso),                                 &
                      pwisolu(kbdim,klev,kwiso),                                 &
                      pwisomfuq(kbdim,klev,kwiso),                               &
                      pwisomful(kbdim,klev,kwiso),pwisodmfup(kbdim,klev,kwiso)

!---wiso-code-end

LOGICAL  :: llo3(kbdim)
INTEGER  :: jl, jt
REAL(dp) :: zzzmb
!
!
!----------------------------------------------------------------------
!
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------
!
100 CONTINUE
!DIR$ IVDEP
!OCL NOVREC
  DO 150 jl=1,kproma
     llo3(jl)=.FALSE.
     IF(.NOT.ldcum(jl).AND.klab(jl,kk+1).EQ.0                          &
              .AND.pqen(jl,kk).GT.0.90_dp*pqsen(jl,kk)) THEN
        llo3(jl)=.TRUE.
        ptu(jl,kk+1)=(pcpen(jl,kk)*pten(jl,kk)                         &
                        +pgeo(jl,kk)-pgeoh(jl,kk+1))/pcpen(jl,kk)
        pqu(jl,kk+1)=pqen(jl,kk)
        plu(jl,kk+1)=0._dp
        zzzmb=MAX(cmfcmin,-pverv(jl,kk)/g)
        zzzmb=MIN(zzzmb,cmfcmax)
        pmfub(jl)=zzzmb
        pmfu(jl,kk+1)=pmfub(jl)
        pmfus(jl,kk+1)=pmfub(jl)*(pcpen(jl,kk)*ptu(jl,kk+1)            &
                                        +pgeoh(jl,kk+1))
        pmfuq(jl,kk+1)=pmfub(jl)*pqu(jl,kk+1)
        pmful(jl,kk+1)=0._dp
        pdmfup(jl,kk+1)=0._dp
        kcbot(jl)=kk
        klab(jl,kk+1)=1
        ktype(jl)=3
        pentr(jl)=entrmid
        IF(lmfdudv) THEN
           puu(jl,kk+1)=puen(jl,kk)
           pvu(jl,kk+1)=pven(jl,kk)
           pmfuu(jl)=pmfub(jl)*puu(jl,kk+1)
           pmfuv(jl)=pmfub(jl)*pvu(jl,kk+1)
        END IF
     END IF
150 END DO

!---wiso-code
  IF (lwiso) THEN
  
!DIR$ IVDEP
!OCL NOVREC
  DO jt=1,kwiso
    DO jl=1,kproma
      IF (llo3(jl)) THEN
        pwisoqu(jl,kk+1,jt)=pwisoqen(jl,kk,jt)
        pwisolu(jl,kk+1,jt)=0._dp
        pwisomfuq(jl,kk+1,jt)=pmfub(jl)*pwisoqu(jl,kk+1,jt)
        pwisomful(jl,kk+1,jt)=0._dp
        pwisodmfup(jl,kk+1,jt)=0._dp
      ENDIF
    END DO
  END DO

  END IF
!---wiso-code-end

!DIR$ IVDEP
!OCL NOVREC
  DO 1504 jt=1,ktrac
     DO 1502 jl=1,kproma
        IF (llo3(jl)) THEN
           pxtu(jl,kk+1,jt)=pxten(jl,kk,jt)
           pmfuxt(jl,kk+1,jt)=pmfub(jl)*pxtu(jl,kk+1,jt)
        ENDIF
1502 END DO
1504 END DO
!
!
  RETURN
END SUBROUTINE cubasmc
