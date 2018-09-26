SUBROUTINE cuadjtq(  kproma, kbdim, klev, kk,             &
           pp,       pt,       pq,       ldflag,   kcall)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!          D.SALMOND         CRAY(UK))      12/8/91
!
!          PURPOSE.
!          --------
!          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM SUBROUTINES:
!              *CUBASE*   (T AND Q AT CONDENSTION LEVEL)
!              *CUASC*    (T AND Q AT CLOUD LEVELS)
!              *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
!          INPUT ARE UNADJUSTED T AND Q VALUES,
!          IT RETURNS ADJUSTED VALUES OF T AND Q
!          NOTE: INPUT PARAMETER KCALL DEFINES CALCULATION AS
!               KCALL=0    ENV. T AND QS IN*CUINI*
!               KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
!               KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
!
!          EXTERNALS
!          ---------
!          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
!          FOR CONDENSATION CALCULATIONS.
!          THE TABLES ARE INITIALISED IN *SETPHYS*.
!

  USE mo_kind,           ONLY: dp
  USE mo_constants,      ONLY: vtmpc1
  USE mo_convect_tables, ONLY: tlucua,   & ! table a
                               tlucub,   & ! table b
                               tlucuc,   & ! table c
                               jptlucu1, jptlucu2, &
                               lookuperror, lookupoverflow

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: kcall, kk, klev, kproma, kbdim

  !  Array arguments with intent(In):
  REAL(dp), INTENT (IN) :: pp(kbdim)
  LOGICAL, INTENT (IN) :: ldflag(kbdim)

  !  Array arguments with intent(InOut):
  REAL(dp), INTENT (INOUT) :: pq(kbdim,klev), pt(kbdim,klev)

  !  Local scalars: 
  REAL(dp):: zcond1, zqst1, zdqsdt, zqsat, zes, zcor, zlcdqsdt
  INTEGER :: isum, jl, it, it1
  LOGICAL :: LO

  !  Local arrays: 
  REAL(dp):: zcond(kbdim)


  !  Executable statements 

  lookupoverflow = .FALSE.

  zcond = 0.0_dp
!
!----------------------------------------------------------------------
!
!     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
!                  -----------------------------------------------------
!
200 CONTINUE
  IF (kcall.EQ.1 ) THEN

     isum=0
!DIR$ IVDEP
!OCL NOVREC
     DO 210 jl=1,kproma
        IF(ldflag(jl)) THEN
           it = NINT(pt(jl,kk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
           zes=tlucua(it)/pp(jl)
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           zqst1=tlucua(it1)/pp(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=zqst1/(1._dp-vtmpc1*zqst1)
           zdqsdt=(zqst1-zqsat)*1000._dp
           zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
           zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
           zcond(jl)=MAX(zcond(jl),0._dp)
           pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
           pq(jl,kk)=pq(jl,kk)-zcond(jl)
           IF(ABS(zcond(jl)) > 0._dp) isum=isum+1
        END IF
210  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (1) ')

     IF(isum.EQ.0) go to 230

!DIR$ IVDEP
!OCL NOVREC
     DO 220 jl=1,kproma
        IF(ldflag(jl)) THEN
           IF(ABS(zcond(jl)) > 0._dp) THEN
           it = NINT(pt(jl,kk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
           zes=tlucua(it)/pp(jl)
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           zqst1=tlucua(it1)/pp(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=zqst1/(1._dp-vtmpc1*zqst1)
           zdqsdt=(zqst1-zqsat)*1000._dp
           zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
           zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
           pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
           pq(jl,kk)=pq(jl,kk)-zcond1
           END IF
        END IF
220  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (2) ')

230  CONTINUE

  END IF

  IF(kcall.EQ.2) THEN

     isum=0
!DIR$ IVDEP
!OCL NOVREC
     DO 310 jl=1,kproma
        IF(ldflag(jl)) THEN
           it = NINT(pt(jl,kk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
           zes=tlucua(it)/pp(jl)
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           zqst1=tlucua(it1)/pp(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=zqst1/(1._dp-vtmpc1*zqst1)
           zdqsdt=(zqst1-zqsat)*1000._dp
           zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
           zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
           zcond(jl)=MIN(zcond(jl),0._dp)
           pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
           pq(jl,kk)=pq(jl,kk)-zcond(jl)
           IF(ABS(zcond(jl)) > 0._dp) isum=isum+1
        END IF
310  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (3) ')

     IF(isum.EQ.0) go to 330

!DIR$ IVDEP
!OCL NOVREC
     DO 320 jl=1,kproma
        IF(ldflag(jl) .AND. ABS(zcond(jl)) > 0._dp) THEN
           it = NINT(pt(jl,kk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
           zes=tlucua(it)/pp(jl)
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           zqst1=tlucua(it1)/pp(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=zqst1/(1._dp-vtmpc1*zqst1)
           zdqsdt=(zqst1-zqsat)*1000._dp
           zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
           zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
           pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
           pq(jl,kk)=pq(jl,kk)-zcond1
        END IF
320  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (4) ')

330  CONTINUE

  END IF

  IF(kcall.EQ.0) THEN

     isum=0
!DIR$ IVDEP
!OCL NOVREC
     DO 410 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/pp(jl)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/pp(jl)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
        zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
        pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
        pq(jl,kk)=pq(jl,kk)-zcond(jl)
        IF(ABS(zcond(jl)) > 0._dp) isum=isum+1
410  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (5) ')

     IF(isum.EQ.0) go to 430

!DIR$ IVDEP
!OCL NOVREC
     DO 420 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/pp(jl)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/pp(jl)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
        zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
        pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
        pq(jl,kk)=pq(jl,kk)-zcond1
420  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (6) ')

430  CONTINUE

  END IF

  IF(kcall.EQ.4) THEN

!DIR$ IVDEP
!OCL NOVREC
     DO 510 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/pp(jl)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/pp(jl)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
        zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
        pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
        pq(jl,kk)=pq(jl,kk)-zcond(jl)
510  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (7) ')

!DIR$ IVDEP
!OCL NOVREC
     DO 520 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/pp(jl)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/pp(jl)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
        zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
        pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
        pq(jl,kk)=pq(jl,kk)-zcond1
520  END DO

     IF (lookupoverflow) CALL lookuperror ('cuadjtq (8) ')

  END IF

  RETURN
END SUBROUTINE cuadjtq
