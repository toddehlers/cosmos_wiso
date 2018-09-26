SUBROUTINE helmo(kdt)

  ! Description:
  !
  ! Computes matrix needed to invert helmoltz equation.
  !
  ! Method:
  !
  ! This subroutine computes the matrix *cn* needed in *sccd*
  ! to invert the *helmoltz equation.
  !
  ! *helmo* is called first from *initialize*: call *helmo*(1).
  ! But, after the first time step *cn* needs to be recomputed. Thus
  ! *helmo* is called a second time, from *stepon*: call *helmo*(2).
  !
  ! Externals:
  ! *getrf* *getri*    called to invert matrix.
  !
  ! 1-appendix b1:organisation of the spectral model.
  ! m.j       5/10/81.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! T. Diehl, DKRZ, July 1999, parallel version
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: finish
  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_tmp_buffer,    ONLY: cn
  USE mo_control,       ONLY: nlev
  USE mo_hyb,           ONLY: bb
  USE mo_constants,     ONLY: a
  USE mo_doctor,        ONLY: nout

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kdt

  !  Local scalars: 
  REAL(dp):: zdtsa2, zk
  INTEGER :: info, jk, jl, jll, i, nns

  !  Local arrays: 
  REAL(dp):: zb(nlev,nlev), zw(2*nlev)
  INTEGER :: ipiv(nlev), nindex(lc%nns)

  !  External subroutines 

  EXTERNAL dgetrf, dgetri

  !  Executable statements 

  nns = lc%nns
  nindex = lc%nindex

!-- 1. Compute matrix *cn*

  zdtsa2 = (kdt/a)**2

  do i=1,nns
     jk=nindex(i)
     zk = (jk-1)*jk*zdtsa2

     DO jl = 1, nlev
        DO jll = 1, nlev
           zb(jl,jll) = zk*bb(jl,jll)
           IF (jl==jll) zb(jl,jll) = zb(jl,jll) + 1._dp
        END DO
     END DO

     CALL dgetrf(nlev,nlev,zb,nlev,ipiv,info)
     IF (info/=0) THEN
        WRITE (nout,*) 'S/DGETRF returns:', info
        CALL finish('helmo','Run terminated.')
     END IF

     CALL dgetri(nlev,zb,nlev,ipiv,zw,2*nlev,info)
     IF (info/=0) THEN
        WRITE (nout,*) 'S/DGETRI returns:', info
        CALL finish('helmo','Run terminated.')
     END IF
     
     DO jl = 1, nlev
        DO jll = 1, nlev
           cn(jl,jll,jk) = zb(jl,jll)
        END DO
        
     END DO

  END DO

  RETURN
END SUBROUTINE helmo
