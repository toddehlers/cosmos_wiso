SUBROUTINE cuini(kproma, kbdim, klev, klevp1, klevm1,                  &
           pten,     pqen,     pqsen,    pxen,     puen,     pven,     &
           ptven,    ktrac,                                            &
           pxten,    pxtenh,   pxtu,     pxtd,     pmfuxt,   pmfdxt,   &
           pverv,    pgeo,     paphp1,   pgeoh,                        &
           ptenh,    pqenh,    pqsenh,   pxenh,    klwmin,             &
           ptu,      pqu,      ptd,      pqd,                          &
           puu,      pvu,      pud,      pvd,                          &
           pmfu,     pmfd,     pmfus,    pmfds,                        &
           pmfuq,    pmfdq,    pdmfup,   pdmfdp,                       &
           pcpen,    pcpcu,                                            &
           pdpmel,   plu,      plude,    pqude,    klab,               &
!---wiso-code
           lwiso, kwiso,                                               &
           pwisoqen, pwisoxen,                                         &
           pwisoqenh,pwisoxenh,                                        &
           pwisoqu,  pwisoqd,                                          &
           pwisomfuq,pwisomfdq,pwisodmfup,pwisodmfdp,                  &
           pwisolu,  pwisolude,pwisoqude                               &
           )         
!---wiso-code-end
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE
!          -------
!
!          THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
!          TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
!          DETERMINES LEVEL OF MAXIMUM VERTICAL VELOCITY
!          AND INITIALIZES VALUES FOR UPDRAFTS AND DOWNDRAFTS
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!
!          METHOD.
!          --------
!          FOR EXTRAPOLATION TO HALF LEVELS SEE TIEDTKE(1989)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* TO SPECIFY QS AT HALF LEVELS
!
USE mo_kind,      ONLY: dp
USE mo_constants, ONLY: rd, cpd

!---wiso-code

USE mo_wiso,      ONLY: tnat, cwisomin

!---wiso-code-end

IMPLICIT NONE

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1, ktrac

!---wiso-code

LOGICAL, INTENT (IN) :: lwiso
INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

REAL(dp):: pten(kbdim,klev),          pqen(kbdim,klev),                &
           puen(kbdim,klev),          pven(kbdim,klev),                &
           pqsen(kbdim,klev),         pverv(kbdim,klev),               &
           pgeo(kbdim,klev),          pgeoh(kbdim,klev),               &
           paphp1(kbdim,klevp1),      ptenh(kbdim,klev),               &
           pxenh(kbdim,klev),         pxen(kbdim,klev),                &
           ptven(kbdim,klev),                                          &
           pqenh(kbdim,klev),         pqsenh(kbdim,klev)
REAL(dp):: pcpen(kbdim,klev),         pcpcu(kbdim,klev)
!
REAL(dp):: ptu(kbdim,klev),           pqu(kbdim,klev),                 &
           ptd(kbdim,klev),           pqd(kbdim,klev),                 &
           puu(kbdim,klev),           pud(kbdim,klev),                 &
           pvu(kbdim,klev),           pvd(kbdim,klev),                 &
           pmfu(kbdim,klev),          pmfd(kbdim,klev),                &
           pmfus(kbdim,klev),         pmfds(kbdim,klev),               &
           pmfuq(kbdim,klev),         pmfdq(kbdim,klev),               &
           pdmfup(kbdim,klev),        pdmfdp(kbdim,klev),              &
           plu(kbdim,klev),           plude(kbdim,klev),               &
           pqude(kbdim,klev)
REAL(dp):: pdpmel(kbdim,klev)
INTEGER :: klab(kbdim,klev),          klwmin(kbdim)
!
REAL(dp):: zwmax(kbdim)
REAL(dp):: zph(kbdim)
LOGICAL :: loflag(kbdim)
REAL(dp):: pxten(kbdim,klev,ktrac),   pxtenh(kbdim,klev,ktrac),        &
           pxtu(kbdim,klev,ktrac),    pxtd(kbdim,klev,ktrac),          &
           pmfuxt(kbdim,klev,ktrac),  pmfdxt(kbdim,klev,ktrac)

!---wiso-code

REAL(dp), OPTIONAL :: pwisoqen(kbdim,klev,kwiso) ,                                &
                      pwisoxenh(kbdim,klev,kwiso), pwisoxen(kbdim,klev,kwiso),    &
                      pwisoqenh(kbdim,klev,kwiso)
!
REAL(dp), OPTIONAL :: pwisoqu(kbdim,klev,kwiso),                                  &
                      pwisoqd(kbdim,klev,kwiso),                                  &
                      pwisomfuq(kbdim,klev,kwiso), pwisomfdq(kbdim,klev,kwiso),   &
                      pwisodmfup(kbdim,klev,kwiso),pwisodmfdp(kbdim,klev,kwiso),  &
                      pwisolu(kbdim,klev,kwiso),   pwisolude(kbdim,klev,kwiso),   &
                      pwisoqude(kbdim,klev,kwiso)

!---wiso-code-end

INTEGER :: jk, jl, jt, ik, icall
REAL(dp):: zarg, zcpm, zzs
!
!  INTRINSIC FUNCTIONS
INTRINSIC MAX, MIN

!---wiso-code

  REAL(dp)            :: zdelta

!---wiso-code-end

!
!----------------------------------------------------------------------
!
!*    1.           SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
!*                 ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
!*                 FIND LEVEL OF MAXIMUM VERTICAL VELOCITY
!                  ----------------------------------------------
!
100 CONTINUE
    DO 101 jk=1,klev
       DO 102 jl=1,kproma
!          pcpen(jl,jk)=cpd*(1.+vtmpc2*MAX(pqen(jl,jk),0.0_dp))
          pcpen(jl,jk)=cpd
102    END DO
101 END DO
  DO 105 jl=1,kproma
     zarg=paphp1(jl,klevp1)/paphp1(jl,klev)
     pgeoh(jl,klev)=rd*ptven(jl,klev)*LOG(zarg)
105 END DO
  DO 107 jk=klevm1,2,-1
     DO 106 jl=1,kproma
        zarg=paphp1(jl,jk+1)/paphp1(jl,jk)
        pgeoh(jl,jk)=pgeoh(jl,jk+1)+rd*ptven(jl,jk)*LOG(zarg)
106  END DO
107 END DO
  DO 130 jk=2,klev
     DO 110 jl=1,kproma
        zcpm=(pcpen(jl,jk)+pcpen(jl,jk-1))*0.5_dp
        ptenh(jl,jk)=(MAX(pcpen(jl,jk-1)*pten(jl,jk-1)+pgeo(jl,jk-1),  &
                   pcpen(jl,jk)*pten(jl,jk)+pgeo(jl,jk))-pgeoh(jl,jk)) &
                    /zcpm
        pqsenh(jl,jk)=pqsen(jl,jk-1)
        zph(jl)=paphp1(jl,jk)
        loflag(jl)=.TRUE.
110  END DO
!
     DO 1104 jt=1,ktrac
        DO 1102 jl=1,kproma
           pxtenh(jl,jk,jt)=(pxten(jl,jk,jt)+pxten(jl,jk-1,jt))*0.5_dp
1102    END DO
1104 END DO
!
!
     ik=jk
     icall=0
     CALL cuadjtq(kproma, kbdim, klev, ik,                             &
                  zph,      ptenh,    pqsenh,   loflag,   icall)
!
     DO 120 jl=1,kproma
        pxenh(jl,jk)=(pxen(jl,jk)+pxen(jl,jk-1))*0.5_dp
        pqenh(jl,jk)=MIN(pqen(jl,jk-1),pqsen(jl,jk-1))                 &
                          +(pqsenh(jl,jk)-pqsen(jl,jk-1))
        pqenh(jl,jk)=MAX(pqenh(jl,jk),0._dp)
!        pcpcu(jl,jk)=cpd*(1.+vtmpc2*pqenh(jl,jk))
        pcpcu(jl,jk)=cpd
120  END DO

!---wiso-code
  IF (lwiso) THEN
  
! Half level values of water isotopes in updrafts.
! Arith. Weighted Mean of the relations (Water isotope/normal Water) for
! and arith. Mean of the absolute values for cloud water

! do not use negative "normal water" values for delta calculation unless
! water on both levels is negative

     DO jt=1,kwiso
       DO jl=1,kproma
         pwisoxenh(jl,jk,jt)=(pwisoxen(jl,jk,jt)+pwisoxen(jl,jk-1,jt))*0.5_dp
         IF ((abs(pqen(jl,jk-1)).lt.cwisomin).and.(abs(pqen(jl,jk)).lt.cwisomin)) THEN
           zdelta=tnat(jt)    
         ELSEIF ((pqen(jl,jk-1).le.0).and.(pqen(jl,jk).gt.0)) THEN
           zdelta=pwisoqen(jl,jk,jt)/pqen(jl,jk)
         ELSEIF ((pqen(jl,jk-1).gt.0).and.(pqen(jl,jk).le.0)) THEN
           zdelta=pwisoqen(jl,jk-1,jt)/pqen(jl,jk-1)
         ELSE
           zdelta=(pwisoqen(jl,jk,jt)+pwisoqen(jl,jk-1,jt))/(pqen(jl,jk)+pqen(jl,jk-1))
         ENDIF
         pwisoqenh(jl,jk,jt)=zdelta*pqenh(jl,jk)
       END DO
     END DO

  END IF
!---wiso-code-end

130 END DO
!
  DO 140 jl=1,kproma
     ptenh(jl,klev)=(pcpen(jl,klev)*pten(jl,klev)+pgeo(jl,klev)-       &
                           pgeoh(jl,klev))/pcpen(jl,klev)
     pxenh(jl,klev)=pxen(jl,klev)
     pqenh(jl,klev)=pqen(jl,klev)
     pcpcu(jl,1)=pcpen(jl,1)
     ptenh(jl,1)=pten(jl,1)
     pxenh(jl,1)=pxen(jl,1)
     pqenh(jl,1)=pqen(jl,1)
     pgeoh(jl,1)=pgeo(jl,1)
     klwmin(jl)=klev
     zwmax(jl)=0._dp
140 END DO
!
!---wiso-code
  IF (lwiso) THEN

     DO jt=1,kwiso
       DO jl=1,kproma
         pwisoxenh(jl,klev,jt)=pwisoxen(jl,klev,jt)
         pwisoqenh(jl,klev,jt)=pwisoqen(jl,klev,jt)
         pwisoxenh(jl,1,jt)=pwisoxen(jl,1,jt)
         pwisoqenh(jl,1,jt)=pwisoqen(jl,1,jt)   
       END DO
     END DO

  END IF
!---wiso-code-end

  DO 1404 jt=1,ktrac
     DO 1402 jl=1,kproma
        pxtenh(jl,klev,jt)=pxten(jl,klev,jt)
        pxtenh(jl,1,jt)=pxten(jl,1,jt)
1402 END DO
1404 END DO
!
!
  DO 160 jk=klevm1,2,-1
     DO 150 jl=1,kproma
        zzs=MAX(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk),                &
                      pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))
        ptenh(jl,jk)=(zzs-pgeoh(jl,jk))/pcpcu(jl,jk)
150  END DO
160 END DO
!
  DO 190 jk=klev,3,-1
!DIR$ IVDEP
!OCL NOVREC
     DO 180 jl=1,kproma
        IF(pverv(jl,jk).LT.zwmax(jl)) THEN
           zwmax(jl)=pverv(jl,jk)
           klwmin(jl)=jk
        END IF
180  END DO
190 END DO
!
!
!-----------------------------------------------------------------------
!*    2.0          INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
!*                 ---------------------------------------------
!
200 CONTINUE
  DO 230 jk=1,klev
     ik=jk-1
     IF(jk.EQ.1) ik=1
     DO 220 jl=1,kproma
        ptu(jl,jk)=ptenh(jl,jk)
        ptd(jl,jk)=ptenh(jl,jk)
        pqu(jl,jk)=pqenh(jl,jk)
        pqd(jl,jk)=pqenh(jl,jk)
        plu(jl,jk)=0._dp
        puu(jl,jk)=puen(jl,ik)
        pud(jl,jk)=puen(jl,ik)
        pvu(jl,jk)=pven(jl,ik)
        pvd(jl,jk)=pven(jl,ik)
        pmfu(jl,jk)=0._dp
        pmfd(jl,jk)=0._dp
        pmfus(jl,jk)=0._dp
        pmfds(jl,jk)=0._dp
        pmfuq(jl,jk)=0._dp
        pmfdq(jl,jk)=0._dp
        pdmfup(jl,jk)=0._dp
        pdmfdp(jl,jk)=0._dp
        pdpmel(jl,jk)=0._dp
        plude(jl,jk)=0._dp
        pqude(jl,jk)=0._dp
        klab(jl,jk)=0
220  END DO
!
!---wiso-code
  IF (lwiso) THEN

     DO jt=1,kwiso
       DO jl=1,kproma
         pwisoqu(jl,jk,jt)=pwisoqenh(jl,jk,jt)
         pwisoqd(jl,jk,jt)=pwisoqenh(jl,jk,jt)
         pwisolu(jl,jk,jt)=0._dp
         pwisomfuq(jl,jk,jt)=0._dp
         pwisomfdq(jl,jk,jt)=0._dp
         pwisodmfup(jl,jk,jt)=0._dp
         pwisodmfdp(jl,jk,jt)=0._dp
         pwisolude(jl,jk,jt)=0._dp
         pwisoqude(jl,jk,jt)=0._dp
       END DO
     END DO

  END IF
!---wiso-code-end

     DO 2204 jt=1,ktrac
        DO 2202 jl=1,kproma
           pxtu(jl,jk,jt)=pxtenh(jl,jk,jt)
           pxtd(jl,jk,jt)=pxtenh(jl,jk,jt)
           pmfuxt(jl,jk,jt)=0._dp
           pmfdxt(jl,jk,jt)=0._dp
2202    END DO
2204 END DO
!
230 END DO
!
  RETURN
END SUBROUTINE cuini
