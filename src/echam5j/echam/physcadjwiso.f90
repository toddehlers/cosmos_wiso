SUBROUTINE physcadjwiso (kproma,kbdim,klev,kwiso,           &
                         pqm1,     pxlm1,     pxim1,        &
                         pqte,     pxlte,     pxite,        &
                         pwisoqm1, pwisoxlm1, pwisoxim1,    &
                         pwisoqte, pwisoxlte, pwisoxite)
!
!      M. WERNER           AWI, BREMERHAVEN      2009
!
!      PURPOSE
!      -------
!      ADJUSTMENT OF WATER ISOTOPE TRACERS TO DEFAULT MODEL WATER VARIABLES Q, XL AND XI
!      (ATTENTION: IT IS ASSUMEND, THAT ISOTOPE TRACER #1 IS SET TO H2-16O)
!
!      INTERFACE
!      ---------
!      THIS SUBROUTINE IS CALLED FROM
!        *PHYSC*
!
!      INPUT  NORMAL AND ISOTOPE VALUES AND TENDENCIES OF
!             VAPOUR, CLOUD WATER AND CLOUD ICE, 
!      OUTPUT NEW ISOTOPE VALUES AND TENDENCIES 
!
!      EXTERNALS
!      ---------
!      NONE

  USE mo_kind,              ONLY: dp
  USE mo_time_control,      ONLY: time_step_len
  USE mo_wiso,              ONLY: nwiso, tnat, cwisomin

  IMPLICIT NONE

! input arguments  
  INTEGER, INTENT(IN) :: kproma,kbdim,klev,kwiso

  REAL(dp), INTENT(IN) :: pqm1(kbdim,klev),                       &
                          pxlm1(kbdim,klev),                      &
                          pxim1(kbdim,klev),                      &
                          pqte(kbdim,klev),                       &
                          pxlte(kbdim,klev),                      &
                          pxite(kbdim,klev)

! input/output arguments  
  REAL(dp), INTENT(INOUT) :: pwisoqm1(kbdim,klev,kwiso),                       &
                             pwisoxlm1(kbdim,klev,kwiso),                      &
                             pwisoxim1(kbdim,klev,kwiso),                      &
                             pwisoqte(kbdim,klev,kwiso),                       &
                             pwisoxlte(kbdim,klev,kwiso),                      &
                             pwisoxite(kbdim,klev,kwiso)
  
! local variables
  REAL(dp) :: zcorr(kbdim)
  
  LOGICAL  :: lcorr(kbdim)

  REAL(dp) :: zqpone,     zxlpone,     zxipone,           &
              zwisoqpone, zwisoxlpone, zwisoxipone
              
  REAL(dp) :: ztwodt

  INTEGER :: jl,jk,jt
    

  ztwodt=time_step_len

!
! time step m1
  DO jk=1,klev
    zcorr(:) = 1._dp
    lcorr(:) = .FALSE.
    jt=1
    DO jl = 1,kproma
      IF (pqm1(jl,jk) .GT.cwisomin .AND. pwisoqm1(jl,jk,jt) .GT.cwisomin) THEN
        lcorr(jl)=.TRUE.
        zcorr(jl) =pqm1(jl,jk) /pwisoqm1(jl,jk,1)
      ENDIF
    ENDDO
    DO jt = 1,kwiso
      DO jl = 1,kproma
        IF (lcorr(jl)) THEN 
          pwisoqm1 (jl,jk,jt)=zcorr(jl) *pwisoqm1 (jl,jk,jt)
        ELSE
          pwisoqm1(jl,jk,jt) = pqm1(jl,jk) * tnat(jt)  
        ENDIF
      ENDDO
    ENDDO
    zcorr(:) = 1._dp
    lcorr(:) = .FALSE.
    jt=1
    DO jl = 1,kproma
      IF (pxlm1(jl,jk) .GT.cwisomin .AND. pwisoxlm1(jl,jk,jt) .GT.cwisomin) THEN
        lcorr(jl)=.TRUE.
        zcorr(jl) =pxlm1(jl,jk) /pwisoxlm1(jl,jk,1)
      ENDIF
    ENDDO
    DO jt = 1,kwiso
      DO jl = 1,kproma
        IF (lcorr(jl)) THEN 
          pwisoxlm1 (jl,jk,jt)=zcorr(jl) *pwisoxlm1 (jl,jk,jt)
        ELSE
          pwisoxlm1(jl,jk,jt) = pxlm1(jl,jk) * tnat(jt)  
        ENDIF
      ENDDO
    ENDDO
    zcorr(:) = 1._dp
    lcorr(:) = .FALSE.
    jt=1
    DO jl = 1,kproma
      IF (pxim1(jl,jk) .GT.cwisomin .AND. pwisoxim1(jl,jk,jt) .GT.cwisomin) THEN
        lcorr(jl)=.TRUE.
        zcorr(jl) =pxim1(jl,jk) /pwisoxim1(jl,jk,1)
      ENDIF
    ENDDO
    DO jt = 1,kwiso
      DO jl = 1,kproma
        IF (lcorr(jl)) THEN 
          pwisoxim1 (jl,jk,jt)=zcorr(jl) *pwisoxim1 (jl,jk,jt)
        ELSE
          pwisoxim1(jl,jk,jt) = pxim1(jl,jk) * tnat(jt)  
        ENDIF
      ENDDO
    ENDDO    
  ENDDO

! time step p1
  DO jk=1,klev
    zcorr(:) = 1._dp
    lcorr(:) = .FALSE.
    DO jl = 1,kproma
      jt=1
      zqpone=pqm1(jl,jk)+pqte(jl,jk)*ztwodt
      zwisoqpone=pwisoqm1(jl,jk,jt)+pwisoqte(jl,jk,jt)*ztwodt
      IF (zqpone .GT.cwisomin .AND. zwisoqpone .GT.cwisomin) THEN
        lcorr(jl)=.TRUE.
        zcorr(jl) =zqpone /zwisoqpone
      ENDIF
    ENDDO
    DO jt = 1,kwiso
      DO jl = 1,kproma
        zqpone=pqm1(jl,jk)+pqte(jl,jk)*ztwodt
        zwisoqpone=pwisoqm1(jl,jk,jt)+pwisoqte(jl,jk,jt)*ztwodt
        IF (lcorr(jl)) THEN 
          zwisoqpone=zcorr(jl) *zwisoqpone
        ELSE
          zwisoqpone = zqpone * tnat(jt)  
        ENDIF
        pwisoqte(jl,jk,jt)=(zwisoqpone-pwisoqm1(jl,jk,jt))/ztwodt
      ENDDO
    ENDDO
    zcorr(:) = 1._dp
    lcorr(:) = .FALSE.
    DO jl = 1,kproma
      jt=1
      zxlpone=pxlm1(jl,jk)+pxlte(jl,jk)*ztwodt
      zwisoxlpone=pwisoxlm1(jl,jk,jt)+pwisoxlte(jl,jk,jt)*ztwodt
      IF (zxlpone .GT.cwisomin .AND. zwisoxlpone .GT.cwisomin) THEN
        lcorr(jl)=.TRUE.
        zcorr(jl) =zxlpone /zwisoxlpone
      ENDIF
    ENDDO
    DO jt = 1,kwiso
      DO jl = 1,kproma
        zxlpone=pxlm1(jl,jk)+pxlte(jl,jk)*ztwodt
        zwisoxlpone=pwisoxlm1(jl,jk,jt)+pwisoxlte(jl,jk,jt)*ztwodt
        IF (lcorr(jl)) THEN 
          zwisoxlpone=zcorr(jl) *zwisoxlpone
        ELSE
          zwisoxlpone = zxlpone * tnat(jt)  
        ENDIF
        pwisoxlte(jl,jk,jt)=(zwisoxlpone-pwisoxlm1(jl,jk,jt))/ztwodt
      ENDDO
    ENDDO
    zcorr(:) = 1._dp
    lcorr(:) = .FALSE.
    DO jl = 1,kproma
      jt=1
      zxipone=pxim1(jl,jk)+pxite(jl,jk)*ztwodt
      zwisoxipone=pwisoxim1(jl,jk,jt)+pwisoxite(jl,jk,jt)*ztwodt
      IF (zxipone .GT.cwisomin .AND. zwisoxipone .GT.cwisomin) THEN
        lcorr(jl)=.TRUE.
        zcorr(jl) =zxipone /zwisoxipone
      ENDIF
    ENDDO
    DO jt = 1,kwiso
      DO jl = 1,kproma
        zxipone=pxim1(jl,jk)+pxite(jl,jk)*ztwodt
        zwisoxipone=pwisoxim1(jl,jk,jt)+pwisoxite(jl,jk,jt)*ztwodt
        IF (lcorr(jl)) THEN 
          zwisoxipone=zcorr(jl) *zwisoxipone
        ELSE
          zwisoxipone = zxipone * tnat(jt)  
        ENDIF
        pwisoxite(jl,jk,jt)=(zwisoxipone-pwisoxim1(jl,jk,jt))/ztwodt
      ENDDO
    ENDDO
  ENDDO
  
  
!      Adjustment of water isotope tracers for negative default model water variables q, xl and xi
!      - for negative default water values the tracers are set to negative values, too, with a delta value equal to SMOW
!
! time step m1
  DO jt = 1, kwiso
    DO jk=1,klev
      DO jl = 1,kproma
        IF (pqm1(jl,jk).LT.0.0_dp) THEN
          pwisoqm1(jl,jk,jt)=tnat(jt)*pqm1(jl,jk)
        ENDIF
        IF (pxlm1(jl,jk).LT.0.0_dp) THEN
          pwisoxlm1(jl,jk,jt)=tnat(jt)*pxlm1(jl,jk)
        ENDIF
        IF (pxim1(jl,jk).LT.0.0_dp) THEN
          pwisoxim1(jl,jk,jt)=tnat(jt)*pxim1(jl,jk)
        ENDIF
      END DO
    END DO
  END DO
! time step p1
  DO jt = 1, kwiso
    DO jk=1,klev
      DO jl = 1,kproma
        zqpone=pqm1(jl,jk)+pqte(jl,jk)*ztwodt
        IF (zqpone.LT.0.0_dp) THEN
          zwisoqpone=tnat(jt)*zqpone
          pwisoqte(jl,jk,jt)=(zwisoqpone-pwisoqm1(jl,jk,jt))/ztwodt
        ENDIF
        zxlpone=pxlm1(jl,jk)+pxlte(jl,jk)*ztwodt
        IF (zxlpone.LT.0.0_dp) THEN
          zwisoxlpone=tnat(jt)*zxlpone
          pwisoxlte(jl,jk,jt)=(zwisoxlpone-pwisoxlm1(jl,jk,jt))/ztwodt
        ENDIF
        zxipone=pxim1(jl,jk)+pxite(jl,jk)*ztwodt
        IF (zxipone.LT.0.0_dp) THEN
          zwisoxipone=tnat(jt)*zxipone
          pwisoxite(jl,jk,jt)=(zwisoxipone-pwisoxim1(jl,jk,jt))/ztwodt
        ENDIF
      END DO
    END DO
  END DO
  

  RETURN
END SUBROUTINE physcadjwiso
