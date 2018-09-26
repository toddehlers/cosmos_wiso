SUBROUTINE cubase(   kproma, kbdim, klev, klevp1, klevm1,              &
           ptenh,    pqenh,    pgeoh,    paph,                         &
           ptu,      pqu,      plu,                                    &
           puen,     pven,     puu,      pvu,                          &
           pcpcu,                                                      &
           ldcum,    kcbot,    klab,                                   &
!---wiso-code
           lwiso, kwiso,                                               &
           pwisoqu,  pwisolu                                           &
           )
!---wiso-code-end
!
!          THIS ROUTINE CALCULATES CLOUD BASE VALUES (T AND Q)
!          FOR CUMULUS PARAMETERIZATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD BASE VALUES FOR CU-PARAMETRIZATION
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD BASE VALUES AND FLAGS AS FOLLOWS;
!                 KLAB=1 FOR SUBCLOUD LEVELS
!                 KLAB=2 FOR CONDENSATION LEVEL
!
!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!          (NON ENTRAINING PLUME,I.E.CONSTANT MASSFLUX)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT
!
USE mo_kind,         ONLY: dp
USE mo_constants,    ONLY: vtmpc1
USE mo_cumulus_flux, ONLY: lmfdudv

IMPLICIT NONE

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1

!---wiso-code

LOGICAL, INTENT (IN) :: lwiso
INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

REAL(dp):: ptenh(kbdim,klev),       pqenh(kbdim,klev),                 &
           pgeoh(kbdim,klev),       paph(kbdim,klevp1)
!
REAL(dp):: ptu(kbdim,klev),         pqu(kbdim,klev),                   &
           plu(kbdim,klev)
REAL(dp):: puen(kbdim,klev),        pven(kbdim,klev),                  &
           puu(kbdim,klev),         pvu(kbdim,klev)
REAL(dp):: pcpcu(kbdim,klev)
INTEGER :: klab(kbdim,klev),        kcbot(kbdim)
LOGICAL :: ldcum(kbdim)
!
!---wiso-code

REAL(dp), OPTIONAL :: pwisoqu(kbdim,klev,kwiso),                                  &
                      pwisolu(kbdim,klev,kwiso)

!---wiso-code-end

REAL(dp):: zqold(kbdim)
REAL(dp):: zph(kbdim)
LOGICAL :: loflag(kbdim)

INTEGER :: jl, jk, is, ik, ikb, icall
REAL(dp):: zbuo, zz

!---wiso-code

REAL(dp):: zqcond(kbdim,klev),zwisoqcond(kbdim,klev,kwiso),zqu_tmp(kbdim)
INTEGER :: jt

!---wiso-code-end

!
!
!----------------------------------------------------------------------
!
!     1.           INITIALIZE VALUES AT LIFTING LEVEL
!                  ----------------------------------
!
100 CONTINUE
  DO 110 jl=1,kproma
     klab(jl,klev)=1
     kcbot(jl)=klevm1
     ldcum(jl)=.FALSE.
     puu(jl,klev)=puen(jl,klev)*(paph(jl,klevp1)-paph(jl,klev))
     pvu(jl,klev)=pven(jl,klev)*(paph(jl,klevp1)-paph(jl,klev))
110 END DO
!
!
!----------------------------------------------------------------------
!
!     2.0          DO ASCENT IN SUBCLOUD LAYER,
!                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
!                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!                  CHECK FOR BUOYANCY AND SET FLAGS
!                  -------------------------------------
!

!---wiso-code
  IF (lwiso) THEN
  
  zqcond(:,:) = 0._dp                ! dummy condensate of normal water at cloud base 
  zwisoqcond(:,:,:) = 0._dp          ! dummy condensate of water isotope at cloud base

  END IF
!---wiso-code-end

200 CONTINUE
  DO 290 jk=klevm1,2,-1
     is=0
     DO 210 jl=1,kproma
        is=is+MERGE(1,0,klab(jl,jk+1).EQ.1)
        loflag(jl)=klab(jl,jk+1).EQ.1
        zph(jl)=paph(jl,jk)
210  END DO
     IF(is.EQ.0) go to 290
     DO 220 jl=1,kproma
        IF(loflag(jl)) THEN
           pqu(jl,jk)=pqu(jl,jk+1)
           ptu(jl,jk)=(pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1)      &
                           -pgeoh(jl,jk))/pcpcu(jl,jk)
           zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk))-ptenh(jl,jk)      &
                                   *(1._dp+vtmpc1*pqenh(jl,jk))+0.5_dp
           IF(zbuo.GT.0._dp) klab(jl,jk)=1
           zqold(jl)=pqu(jl,jk)
        END IF
220  END DO
!
!---wiso-code
  IF (lwiso) THEN

     DO jt=1,kwiso
       DO jl=1,kproma
         IF(loflag(jl)) THEN
           pwisoqu(jl,jk,jt)=pwisoqu(jl,jk+1,jt)
         END IF
       END DO
     END DO

  END IF
!---wiso-code-end

     ik=jk
     icall=1
     zqu_tmp(:)=pqu(:,ik)                      !--- wiso-code: Store old value of pqu
     CALL cuadjtq(kproma, kbdim, klev, ik,                             &
                  zph,      ptu,      pqu,      loflag,   icall)
!
!DIR$ IVDEP
!OCL NOVREC
     DO 240 jl=1,kproma
        zqcond(jl,jk)=zqu_tmp(jl)-pqu(jl,jk)   !--- wiso-code: Storing the condensate of normal water
        IF(loflag(jl).AND.pqu(jl,jk).LT.zqold(jl)) THEN
           klab(jl,jk)=2
           plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
           zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk)-plu(jl,jk))-       &
                        ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))+0.5_dp
           IF(zbuo.GT.0.) THEN
              kcbot(jl)=jk
              ldcum(jl)=.TRUE.
           END IF
        END IF
240  END DO

!---wiso-code
  IF (lwiso) THEN

! Calculate the Condensate of water isotope (including fractionation)
! -> Subroutine cuadjwisoq returns new values of pwisoqu
    ik = jk
    CALL cuadjwisoq(kproma,   kbdim,      klev,      kwiso,    ik,     &
                    ptu,      pqu,        pwisoqu,   pwisolu,          &
                    zqcond,   zwisoqcond, loflag)

  END IF
!---wiso-code-end
!
!             CALCULATE AVERAGES OF U AND V FOR SUBCLOUD ARA,.
!             THE VALUES WILL BE USED TO DEFINE CLOUD BASE VALUES.
!
     IF(lmfdudv) THEN
        DO 250 jl=1,kproma
           IF(jk.GE.kcbot(jl)) THEN
              puu(jl,klev)=puu(jl,klev)+                               &
                             puen(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
              pvu(jl,klev)=pvu(jl,klev)+                               &
                             pven(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
           END IF
250     END DO
     END IF
!
290 END DO
!
!
  IF(lmfdudv) THEN
     DO 310 jl=1,kproma
        IF(ldcum(jl)) THEN
           ikb=kcbot(jl)
           zz=1._dp/(paph(jl,klevp1)-paph(jl,ikb))
           puu(jl,klev)=puu(jl,klev)*zz
           pvu(jl,klev)=pvu(jl,klev)*zz
        ELSE
           puu(jl,klev)=puen(jl,klevm1)
           pvu(jl,klev)=pven(jl,klevm1)
        END IF
310  END DO
  END IF
!
  RETURN
END SUBROUTINE cubase
