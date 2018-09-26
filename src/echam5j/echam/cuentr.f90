SUBROUTINE cuentr(   kproma, kbdim, klev, klevp1, kk,                  &
           ptenh,    pqenh,    pqte,     paphp1,                       &
           klwmin,   ldcum,    ktype,    kcbot,    kctop0,             &
           ppbase,   pmfu,     pentr,    podetr,                       &
           khmin,    pgeoh,                                            &
           pdmfen,   pdmfde)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
!          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          AND UPDRAFT VALUES T,Q ETC
!          IT RETURNS ENTRAINMENT/DETRAINMENT RATES
!
!          METHOD.
!          --------
!          S. TIEDTKE (1989)
!
!          EXTERNALS
!          ---------
!          NONE
!
USE mo_kind,           ONLY: dp
USE mo_constants,      ONLY: g, rd, vtmpc1
USE mo_cumulus_flux,   ONLY: centrmax, cmfcmin
!
IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kbdim, klev, klevp1, kproma, kk
!
REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            paphp1(kbdim,klevp1),                                      &
            pmfu(kbdim,klev),        pqte(kbdim,klev),                 &
            pentr(kbdim),            ppbase(kbdim)
REAL(dp) :: podetr(kbdim,klev)
REAL(dp) :: pgeoh (kbdim,klev)
INTEGER  :: khmin (kbdim)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                     &
            kcbot(kbdim),            kctop0(kbdim)
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: pdmfen(kbdim),           pdmfde(kbdim)
!
LOGICAL  :: llo1,llo2
!
INTEGER  :: jl, ikb, ikt, ikh, iklwmin
REAL(dp) :: zrg, zrrho, zdprho, zpmid, zentr, zentest, zzmzk, ztmzk    &
          , zorgde, zarg
!
!----------------------------------------------------------------------
!
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------
!
100 CONTINUE
!
!
!*    1.1          SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
!                  --------------------------------------------
!
110 CONTINUE
!
!
!*    1.2          SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
!                  -----------------------------------------
!
120 CONTINUE
  zrg=1._dp/g
  DO 125 jl=1,kproma
     ppbase(jl)=paphp1(jl,kcbot(jl))
     zrrho=(rd*ptenh(jl,kk+1)*(1._dp+vtmpc1*pqenh(jl,kk+1)))           &
                    /paphp1(jl,kk+1)
     zdprho=(paphp1(jl,kk+1)-paphp1(jl,kk))*zrg
     zpmid=0.5_dp*(ppbase(jl)+paphp1(jl,kctop0(jl)))
     zentr=pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho
     llo1=kk.LT.kcbot(jl).AND.ldcum(jl)
     pdmfde(jl)=MERGE(zentr,0._dp,llo1)
     llo2=llo1.AND.ktype(jl).EQ.2.AND.                                 &
                   (ppbase(jl)-paphp1(jl,kk).LT.0.2e5_dp.OR.           &
                    paphp1(jl,kk).GT.zpmid)
     pdmfen(jl)=MERGE(zentr,0._dp,llo2)
     iklwmin=MAX(klwmin(jl),kctop0(jl)+2)
     llo2=llo1 .AND. ktype(jl).EQ.3 .AND. kk .GE. iklwmin
     IF(llo2) pdmfen(jl)=zentr
     IF(llo2 .AND. pqenh(jl,kk+1).GT.1.E-5_dp) THEN
        pmfu(jl,kk+1) = MAX(pmfu(jl,kk+1),cmfcmin)
        zentest = MAX(pqte(jl,kk),0._dp)/pqenh(jl,kk+1)
        zentest = MIN(centrmax,zentest/(pmfu(jl,kk+1)*zrrho))
        pdmfen(jl) = zentr+zentest*pmfu(jl,kk+1)*zrrho*zdprho
     ENDIF
     llo2=llo1 .AND. ktype(jl).EQ.1 .AND.                              &
                   (kk.GE.iklwmin.OR.paphp1(jl,kk).GT.zpmid)
     IF(llo2) pdmfen(jl)=zentr
!
!    organized detrainment, detrainment starts at khmin
!
     llo2=llo1 .AND. ktype(jl).EQ.1
     ikb=kcbot(jl)
     podetr(jl,kk)=0._dp
     IF(llo2.AND.kk.LE.khmin(jl).AND.kk.GE.kctop0(jl)) THEN
        ikt=kctop0(jl)
        ikh=khmin(jl)
        IF(ikh.GT.ikt) THEN
           zzmzk  =-(pgeoh(jl,ikh)-pgeoh(jl,kk))*zrg
           ztmzk  =-(pgeoh(jl,ikh)-pgeoh(jl,ikt))*zrg
           zarg  =3.1415_dp*(zzmzk/ztmzk)*0.5_dp
           zorgde=TAN(zarg)*3.1415_dp*0.5_dp/ztmzk
           zdprho=(paphp1(jl,kk+1)-paphp1(jl,kk))*(zrg*zrrho)
           podetr(jl,kk)=MIN(zorgde,centrmax)*pmfu(jl,kk+1)*zdprho
        ENDIF
     ENDIF
125 END DO
!
  RETURN
END SUBROUTINE cuentr
