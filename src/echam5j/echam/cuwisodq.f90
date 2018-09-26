SUBROUTINE cuwisodq(kproma, kbdim, klev, klevp1, ktopm2, ldcum, kwiso, &
                  paphp1,   pwisoqte,                                  &
                  pwisoxtec,                                           &
                  pwisomfuq,pwisomfdq,                                 &
                  pwisomful,pwisodmfup,pwisodmfdp,pwisolude,           &
                  pwisoqtec,pwisoqude)
!
! DESCRIPTION:
!
! UPDATES WATER ISOTOPE TENDENCIES, PRECIPITATION RATES DOES GLOBAL DIAGNOSTICS
!
! METHOD:
!
! *CUWISODQ* IS CALLED FROM *CUMASTR*
!
! AUTHORS:
!
! G.HOFFMANN, MPI MET, HAMBURG, 1992 
! ADAPTED TO F90: M. WERNER, MPI BGC, JENA, 2004 
! ADAPTED TO ECHAM5: M. WERNER, AWI, BREMERHAVEN, 2009
!

USE mo_kind,         ONLY : dp
USE mo_constants,    ONLY : g
USE mo_time_control, ONLY : delta_time

IMPLICIT NONE

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktopm2, kwiso

REAL(dp) :: pwisoqte(kbdim,klev,kwiso),                                &
            paphp1(kbdim,klevp1)

REAL(dp) :: pwisomfuq(kbdim,klev,kwiso), pwisomfdq(kbdim,klev,kwiso),  &
            pwisomful(kbdim,klev,kwiso), pwisolude(kbdim,klev,kwiso),  &
            pwisodmfup(kbdim,klev,kwiso),pwisodmfdp(kbdim,klev,kwiso), &
            pwisoqtec(kbdim,klev,kwiso), pwisoqude(kbdim,klev,kwiso),  &
            pwisoxtec(kbdim,klev,kwiso)

LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: zwisodqdt

INTEGER  :: jl, jk, jt

!  Executable statements 

  DO jk=ktopm2,klev
!
    IF(jk.LT.klev) THEN
      DO jt=1,kwiso
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            zwisodqdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*             &
                      (pwisomfuq(jl,jk+1,jt)-pwisomfuq(jl,jk,jt)+      &
                       pwisomfdq(jl,jk+1,jt)-pwisomfdq(jl,jk,jt)+      &
                       pwisomful(jl,jk+1,jt)-pwisomful(jl,jk,jt)-      &
                       pwisolude(jl,jk,jt)-                            &
                      (pwisodmfup(jl,jk,jt)+pwisodmfdp(jl,jk,jt)))
            pwisoqte(jl,jk,jt)=pwisoqte(jl,jk,jt)+zwisodqdt
            pwisoxtec(jl,jk,jt)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))    &
                                *pwisolude(jl,jk,jt)
            pwisoqtec(jl,jk,jt)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))    &
                                *pwisoqude(jl,jk,jt)
          ENDIF
        END DO
      END DO
    ELSE
      DO jt=1,kwiso
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            zwisodqdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*            &
                       (pwisomfuq(jl,jk,jt)+pwisomfdq(jl,jk,jt)+       &
                        pwisolude(jl,jk,jt)+                           &
                       (pwisomful(jl,jk,jt)+pwisodmfup(jl,jk,jt)+      &
                        pwisodmfdp(jl,jk,jt)))
            pwisoqte(jl,jk,jt)=pwisoqte(jl,jk,jt)+zwisodqdt
            pwisoxtec(jl,jk,jt)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))    &
                                *pwisolude(jl,jk,jt)
            pwisoqtec(jl,jk,jt)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))    &
                                *pwisoqude(jl,jk,jt)
          ENDIF
        END DO
      END DO
    END IF
!
  END DO

  RETURN
END SUBROUTINE cuwisodq
