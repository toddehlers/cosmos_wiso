SUBROUTINE cucall(   kproma, kbdim, klev, klevp1, klevm1, ilab,        &
                     ktrac,                                            &
                     pxtm1,    pxtte,                                  &
                     ptm1,     pqm1,     pum1,     pvm1,               &
                     pxlm1,    pxim1,                                  &
                     ptte,     pqte,     pvom,     pvol,               &
                     pxlte,    pxite,                                  &
                     pverv,    pxtec,    pqtec,    pqhfla,             &
                     papp1,    paphp1,   pgeo,                         &
                     prsfc,    pssfc,    paprc,    paprs,              &
                     ktype,    ldland,                                 &
                     ptopmax,                                          &
!---wiso-code
                     lwiso, kwiso,                                     &
                     pwisoqm1,                                         &
                     pwisoxlm1,pwisoxim1,                              &
                     pwisoqte,                                         &
                     pwisoxlte,pwisoxite,                              &
                     pwisoxtec,pwisoqtec,                              &
                     pwisorsfc,pwisossfc,pwisoaprc, pwisoaprs)
!---wiso-code-end
                     
!
!
!          *CUCALL* - MASTER ROUTINE - PROVIDES INTERFACE FOR:
!                     *CUMASTR* (CUMULUS PARAMETERIZATION)
!
!           M.TIEDTKE      E.C.M.W.F.     12/1989
!
!**   PURPOSE.
!     --------
!
!          *CUCALL* - INTERFACE FOR *CUMASTR*:
!                     PROVIDES INPUT FOR CUMASTR
!                     RECEIVES UPDATED TENDENCIES, PRECIPITATION.
!
!**   INTERFACE.
!     ----------
!
!          *CUCALL* IS CALLED FROM *PHYSC*
!
!     EXTERNALS.
!     ----------
!
!          CUMASTR, CUMASTRT OR CUMASTRH
!
  USE mo_kind,           ONLY: dp
  USE mo_time_control,   ONLY: time_step_len
  USE mo_constants,      ONLY: vtmpc1    ! vtmpc1=rv/rd-1
  USE mo_param_switches, ONLY: iconv
  USE mo_convect_tables, ONLY: tlucua, jptlucu1, jptlucu2,             &
                               lookuperror, lookupoverflow

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: klev, klevm1, klevp1, kproma, kbdim, ktrac

  REAL(dp)::  ptm1(kbdim,klev),         pqm1(kbdim,klev),              &
              pum1(kbdim,klev),         pvm1(kbdim,klev),              &
              ptte(kbdim,klev),         pqte(kbdim,klev),              &
              pvom(kbdim,klev),         pvol(kbdim,klev),              &
              pverv(kbdim,klev),        pgeo(kbdim,klev),              &
              papp1(kbdim,klev),        paphp1(kbdim,klevp1)
  REAL(dp)::  paprc(kbdim),             paprs(kbdim),                  &
              prsfc(kbdim),             pssfc(kbdim)
  INTEGER ::  ktype(kbdim)
  REAL(dp)::  pqhfla(kbdim)
  REAL(dp)::  ptopmax(kbdim)
  INTEGER ::  ilab(kbdim,klev)
  REAL(dp)::  pxtec(kbdim,klev),        pqtec(kbdim,klev)
  REAL(dp)::  pxlm1(kbdim,klev),        pxim1(kbdim,klev),             &
              pxlte(kbdim,klev),        pxite(kbdim,klev)

!---wiso-code

  LOGICAL, INTENT (IN) :: lwiso
  INTEGER, INTENT (IN) :: kwiso

  REAL(dp), OPTIONAL ::  pwisoqm1(kbdim,klev,kwiso), pwisoqte(kbdim,klev,kwiso)
  REAL(dp), OPTIONAL ::  pwisoaprc(kbdim,kwiso),     pwisoaprs(kbdim,kwiso),      &
                         pwisorsfc(kbdim,kwiso),     pwisossfc(kbdim,kwiso)
  REAL(dp), OPTIONAL ::  pwisoxtec(kbdim,klev,kwiso),pwisoqtec(kbdim,klev,kwiso)
  REAL(dp), OPTIONAL ::  pwisoxlm1(kbdim,klev,kwiso),pwisoxim1(kbdim,klev,kwiso), &
                         pwisoxlte(kbdim,klev,kwiso),pwisoxite(kbdim,klev,kwiso)
              
!---wiso-code-end

  REAL(dp)::  ztp1(kbdim,klev),         zqp1(kbdim,klev),              &
              zxp1(kbdim,klev),         ztvp1(kbdim,klev),             &
              zup1(kbdim,klev),         zvp1(kbdim,klev),              &
              ztu(kbdim,klev),          zqu(kbdim,klev),               &
              zlu(kbdim,klev),          zlude(kbdim,klev),             &
              zqude(kbdim,klev),                                       &
              zmfu(kbdim,klev),         zmfd(kbdim,klev),              &
              zqsat(kbdim,klev),        zrain(kbdim)
  INTEGER ::  itopec2(kbdim)
  INTEGER ::  icbot(kbdim),             ictop(kbdim)
  REAL(dp)::  zxtp1(kbdim,klev,ktrac),  zxtu(kbdim,klev,ktrac),        &
              pxtm1(kbdim,klev,ktrac),  pxtte(kbdim,klev,ktrac)
  REAL(dp)::  ztopmax(kbdim)
  LOGICAL ::  locum(kbdim),             ldland(kbdim)

  !  Local scalars: 
  REAL(dp):: ztmst, zxlp1, zxip1
  INTEGER :: ilevmin, jk, jl, jt, it

!---wiso-code

  REAL(dp)::  zwisoqp1(kbdim,klev,kwiso),                              &
              zwisoxp1(kbdim,klev,kwiso),                              &
              zwisoqu(kbdim,klev,kwiso),                               &
              zwisolu(kbdim,klev,kwiso),  zwisolude(kbdim,klev,kwiso), &
              zwisoqude(kbdim,klev,kwiso),                             &
              zwisorain(kbdim,kwiso)

  !  Local scalars: 
  REAL(dp)::  zwisoxlp1, zwisoxip1

!---wiso-code-end

  !  Executable statements 

  lookupoverflow = .FALSE.

!
!-----------------------------------------------------------------------
!*    1.           CALCULATE T,Q AND QS AT MAIN LEVELS
!*                 -----------------------------------
!
!
100 CONTINUE
  ztmst=time_step_len
  DO 120 jk=1,klev
     DO 110 jl=1,kproma
        ztp1(jl,jk)=ptm1(jl,jk)+ptte(jl,jk)*ztmst
        zqp1(jl,jk)=MAX(0._dp,pqm1(jl,jk)+pqte(jl,jk)*ztmst)
        zxlp1=pxlm1(jl,jk)+pxlte(jl,jk)*ztmst
        zxip1=pxim1(jl,jk)+pxite(jl,jk)*ztmst
        zxp1(jl,jk)=MAX(0._dp,zxlp1+zxip1)
        ztvp1(jl,jk)=ztp1(jl,jk)*(1._dp+vtmpc1*zqp1(jl,jk)-zxp1(jl,jk))
        zup1(jl,jk)=pum1(jl,jk)+pvom(jl,jk)*ztmst
        zvp1(jl,jk)=pvm1(jl,jk)+pvol(jl,jk)*ztmst
        it = INT(ztp1(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsat(jl,jk)=tlucua(it)/papp1(jl,jk)
        zqsat(jl,jk)=MIN(0.5_dp,zqsat(jl,jk))
        zqsat(jl,jk)=zqsat(jl,jk)/(1._dp-vtmpc1*zqsat(jl,jk))
110  END DO

     IF (lookupoverflow) CALL lookuperror ('cucall      ')

     DO 1104 jt=1,ktrac
        DO 1102 jl=1,kproma
           zxtp1(jl,jk,jt)=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst
1102    END DO
1104 END DO

120 END DO

!---wiso-code
  IF (lwiso) THEN
  
  DO jt=1,kwiso
    DO jk=1,klev
      DO jl=1,kproma
        zwisoqp1(jl,jk,jt)=MAX(0._dp,pwisoqm1(jl,jk,jt)+pwisoqte(jl,jk,jt)*ztmst)
        zwisoxlp1=pwisoxlm1(jl,jk,jt)+pwisoxlte(jl,jk,jt)*ztmst
        zwisoxip1=pwisoxim1(jl,jk,jt)+pwisoxite(jl,jk,jt)*ztmst
        zwisoxp1(jl,jk,jt)=MAX(0._dp,zwisoxlp1+zwisoxip1)
      END DO
    END DO
  END DO

  END IF
!---wiso-code-end

  DO 130 jl=1,kproma
     zrain(jl)=0._dp
     locum(jl)=.FALSE.
130 END DO

!---wiso-code
  IF (lwiso) THEN

  DO jt=1,kwiso
    DO jl=1,kproma
      zwisorain(jl,jt)=0._dp
    END DO
  END DO

  END IF
!---wiso-code-end
!
!
!-----------------------------------------------------------------------
!
!*    2.     CALL 'CUMASTR'(MASTER-ROUTINE FOR CUMULUS PARAMETERIZATION)
!*           -----------------------------------------------------------
!
!
200 CONTINUE

  IF (lwiso) THEN

  SELECT CASE (iconv)
  CASE(1)
     CALL cumastr(kproma, kbdim, klev, klevp1, klevm1, ilab,           &
                  ztp1,     zqp1,     zxp1,     zup1,   zvp1,          &
                  ztvp1,    ktrac,    ldland,                          &
                  zxtp1,    zxtu,     pxtte,                           &
                  pverv,    zqsat,    pqhfla,                          &
                  paphp1,   pgeo,                                      &
                  ptte,     pqte,     pvom,     pvol,                  &
                  prsfc,    pssfc,    paprc,    paprs,  pxtec,         &
                  pqtec,    zqude,                                     &
                  locum,    ktype,    icbot,    ictop,                 &
                  ztu,      zqu,      zlu,      zlude,                 &
                  zmfu,     zmfd,     zrain,                           &
!---wiso-code
                  lwiso, kwiso,                                        &
                  zwisoqp1, zwisoxp1,                                  &
                  pwisoqte,                                            &
                  pwisorsfc,pwisossfc,pwisoaprc,pwisoaprs,pwisoxtec,   &
                  pwisoqtec,zwisoqude,                                 &
                  zwisoqu,  zwisolu,  zwisolude,                       &
                  zwisorain  )
!---wiso-code-end
  CASE(2)
     CALL cumastrt(kproma, kbdim, klev, klevp1, klevm1, ilab,          &
                  ztp1,     zqp1,     zxp1,     zup1,   zvp1,          &
                  ztvp1,    ktrac,    ldland,                          &
                  zxtp1,    zxtu,     pxtte,                           &
                  pverv,    zqsat,    pqhfla,                          &
                  paphp1,   pgeo,                                      &
                  ptte,     pqte,     pvom,     pvol,                  &
                  prsfc,    pssfc,    paprc,    paprs,  pxtec,         &
                  pqtec,    zqude,                                     &
                  locum,    ktype,    icbot,    ictop,                 &
                  ztu,      zqu,      zlu,      zlude,                 &
                  zmfu,     zmfd,     zrain,                           &
!---wiso-code
                  lwiso, kwiso,                                        &
                  zwisoqp1, zwisoxp1,                                  &
                  pwisoqte,                                            &
                  pwisorsfc,pwisossfc,pwisoaprc,pwisoaprs,pwisoxtec,   &
                  pwisoqtec,zwisoqude,                                 &
                  zwisoqu,  zwisolu,  zwisolude,                       &
                  zwisorain  )
!---wiso-code-end
  CASE(3)
     CALL cumastrh(kproma, kbdim, klev, klevp1, klevm1, ilab,          &
                  ztp1,     zqp1,     zxp1,     zup1,   zvp1,          &
                  ztvp1,    ktrac,    ldland,                          &
                  zxtp1,    zxtu,     pxtte,                           &
                  pverv,    zqsat,    pqhfla,                          &
                  paphp1,   pgeo,                                      &
                  ptte,     pqte,     pvom,     pvol,                  &
                  prsfc,    pssfc,    paprc,    paprs,  pxtec,         &
                  pqtec,    zqude,                                     &
                  locum,    ktype,    icbot,    ictop,                 &
                  ztu,      zqu,      zlu,      zlude,                 &
                  zmfu,     zmfd,     zrain,                           &
!---wiso-code
                  lwiso, kwiso,                                        &
                  zwisoqp1, zwisoxp1,                                  &
                  pwisoqte,                                            &
                  pwisorsfc,pwisossfc,pwisoaprc,pwisoaprs,pwisoxtec,   &
                  pwisoqtec,zwisoqude,                                 &
                  zwisoqu,  zwisolu,  zwisolude,                       &
                  zwisorain  )
!---wiso-code-end

  END SELECT
  
  ELSE

  SELECT CASE (iconv)
  CASE(1)
     CALL cumastr(kproma, kbdim, klev, klevp1, klevm1, ilab,           &
                  ztp1,     zqp1,     zxp1,     zup1,   zvp1,          &
                  ztvp1,    ktrac,    ldland,                          &
                  zxtp1,    zxtu,     pxtte,                           &
                  pverv,    zqsat,    pqhfla,                          &
                  paphp1,   pgeo,                                      &
                  ptte,     pqte,     pvom,     pvol,                  &
                  prsfc,    pssfc,    paprc,    paprs,  pxtec,         &
                  pqtec,    zqude,                                     &
                  locum,    ktype,    icbot,    ictop,                 &
                  ztu,      zqu,      zlu,      zlude,                 &
                  zmfu,     zmfd,     zrain,                           &
!---wiso-code
                  lwiso, kwiso)
!---wiso-code-end
  CASE(2)
     CALL cumastrt(kproma, kbdim, klev, klevp1, klevm1, ilab,          &
                  ztp1,     zqp1,     zxp1,     zup1,   zvp1,          &
                  ztvp1,    ktrac,    ldland,                          &
                  zxtp1,    zxtu,     pxtte,                           &
                  pverv,    zqsat,    pqhfla,                          &
                  paphp1,   pgeo,                                      &
                  ptte,     pqte,     pvom,     pvol,                  &
                  prsfc,    pssfc,    paprc,    paprs,  pxtec,         &
                  pqtec,    zqude,                                     &
                  locum,    ktype,    icbot,    ictop,                 &
                  ztu,      zqu,      zlu,      zlude,                 &
                  zmfu,     zmfd,     zrain,                           &
!---wiso-code
                  lwiso, kwiso)
!---wiso-code-end
  CASE(3)
     CALL cumastrh(kproma, kbdim, klev, klevp1, klevm1, ilab,          &
                  ztp1,     zqp1,     zxp1,     zup1,   zvp1,          &
                  ztvp1,    ktrac,    ldland,                          &
                  zxtp1,    zxtu,     pxtte,                           &
                  pverv,    zqsat,    pqhfla,                          &
                  paphp1,   pgeo,                                      &
                  ptte,     pqte,     pvom,     pvol,                  &
                  prsfc,    pssfc,    paprc,    paprs,  pxtec,         &
                  pqtec,    zqude,                                     &
                  locum,    ktype,    icbot,    ictop,                 &
                  ztu,      zqu,      zlu,      zlude,                 &
                  zmfu,     zmfd,     zrain,                           &
!---wiso-code
                  lwiso, kwiso)
!---wiso-code-end

  END SELECT

  END IF
!
!
! ------------------------------------------------------------------
!
!*     3.     PRESSURE ALTITUDE OF CONVECTIVE CLOUD TOPS.
!             -------- -------- -- ---------- ----- -----
!
300 CONTINUE
!
  ilevmin=klev-4
!
  DO 301 jl=1,kproma
     itopec2(jl)=klevp1
301 END DO
!
  DO 303 jk=1,ilevmin
     DO 302 jl=1,kproma
        IF(ilab(jl,jk).EQ.2 .AND. itopec2(jl).EQ.klevp1) THEN
           itopec2(jl)=jk
        END IF
302  END DO
303 END DO
!
  ztopmax(1:kproma) = ptopmax(1:kproma)

  DO 304 jl=1,kproma
     IF(itopec2(jl).EQ.1) THEN
        ptopmax(jl)=papp1(jl,1)
     ELSE IF(itopec2(jl).NE.klevp1) THEN
        ptopmax(jl)=paphp1(jl,itopec2(jl))
     ELSE
        ptopmax(jl)=99999._dp
     END IF
     ptopmax(jl)=MIN(ptopmax(jl),ztopmax(jl))
304 END DO
!
!
!---------------------------------------------------------------------
!
  RETURN
END SUBROUTINE cucall
