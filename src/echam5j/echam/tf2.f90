!OCL NOALIAS

SUBROUTINE tf2

  ! Description:
  !
  ! This subroutine completes the second part of the
  ! time filtering:xm1=xm1+eps*x.
  !
  ! Method:
  !
  ! *tf2* is called from *scans1l*.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, May 1991, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 1998, tendency diagnostics
  ! U. Schlese, DKRZ, November 1999, cloud ice added
  ! I. Kirchner, MPI, December 2000, time control
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
  USE mo_control,       ONLY: ltdiag, nlev
  USE mo_semi_impl,     ONLY: eps
  USE mo_time_control,  ONLY: get_time_step, INIT_STEP
  USE mo_tracer,        ONLY: trlist
  USE mo_scan_buffer,   ONLY: alps, d, t, u, v, vo, dudl, dvdl 
  USE mo_memory_gl,     ONLY: q, xl, xi, xt
  USE mo_memory_g1a,    ONLY: alpsm1, dm1, qm1, tm1, vom1, xlm1, xim1, &
                              xtm1
  USE mo_memory_g2a,    ONLY: um1, vm1, dudlm1, dvdlm1
  USE mo_diag_tendency, ONLY: ptfh1, ptfh2, pdiga, pdsga, lset_fh1
  USE mo_decomposition, ONLY: ldc => local_decomposition
!---wiso-code
  USE mo_memory_wiso,   ONLY: wisoq,   wisoxl,   wisoxi,               &
                              wisoqm1, wisoxlm1, wisoxim1
  USE mo_wiso,          ONLY: lwiso, nwiso
!---wiso-code-end

  IMPLICIT NONE

  !  Local scalars: 
  REAL(dp):: zeps
  INTEGER :: jk, jl, jt, jrow
  INTEGER :: nglon, nproma, ngpblks, istep

  !  Executable statements 

  ! local array bounds

  nglon   = ldc% nglon   ! number of longitudes
  ngpblks = ldc% ngpblks ! number of rows

  ! local latitude index, continuous north -> south

  !-- 1. Second part of time filtering

  istep = get_time_step()
  IF ((istep - INIT_STEP) <= 1) THEN
     zeps = 0._dp
  ELSE
     zeps = eps
  END IF

!CSD$ PARALLEL DO PRIVATE(nproma,jk,jl,jt)
!$OMP PARALLEL PRIVATE(nproma,jk,jl,jt,jrow)
!$OMP DO
  DO jrow = 1, ngpblks

    IF ( jrow == ngpblks ) THEN
      nproma = ldc% npromz
    ELSE
      nproma = ldc% nproma
    END IF

    DO jk = 1, nlev
!DIR$ IVDEP
      DO jl = 1, nproma
        vom1   (jl,jk,jrow) = vom1(jl,jk,jrow)    + zeps*vo (jl,jk,jrow)
        dm1    (jl,jk,jrow) = dm1 (jl,jk,jrow)    + zeps* d (jl,jk,jrow)
        qm1    (jl,jk,jrow) = qm1 (jl,jk,jrow)    + zeps* q (jl,jk,jrow)
        xlm1   (jl,jk,jrow) = xlm1(jl,jk,jrow)    + zeps* xl(jl,jk,jrow)
        xim1   (jl,jk,jrow) = xim1(jl,jk,jrow)    + zeps* xi(jl,jk,jrow)
        tm1    (jl,jk,jrow) = tm1 (jl,jk,jrow)    + zeps* t (jl,jk,jrow)
        um1    (jl,jk,jrow) = um1 (jl,jk,jrow)    + zeps* u (jl,jk,jrow)
        vm1    (jl,jk,jrow) = vm1 (jl,jk,jrow)    + zeps* v (jl,jk,jrow) 

        dudlm1 (jl,jk,jrow) = dudlm1(jl,jk,jrow)  + zeps* dudl(jl,jk,jrow)
        dvdlm1 (jl,jk,jrow) = dvdlm1(jl,jk,jrow)  + zeps* dvdl(jl,jk,jrow)
     END DO
    END DO

!---wiso-code

    IF (lwiso) THEN
       DO jt = 1, nwiso
          DO jk = 1, nlev
!DIR$ IVDEP
             DO jl = 1, nproma
                wisoqm1 (jl,jk,jt,jrow) = wisoqm1 (jl,jk,jt,jrow) + zeps*wisoq (jl,jk,jt,jrow)
                wisoxlm1(jl,jk,jt,jrow) = wisoxlm1(jl,jk,jt,jrow) + zeps*wisoxl(jl,jk,jt,jrow)
                wisoxim1(jl,jk,jt,jrow) = wisoxim1(jl,jk,jt,jrow) + zeps*wisoxi(jl,jk,jt,jrow)
             END DO
          END DO
       END DO
    ENDIF

!---wiso-code-end

    DO jt = 1, trlist% ntrac
      IF (trlist% ti(jt)% nint /= 1) CYCLE
      DO jk = 1, nlev
!DIR$ IVDEP
        DO jl = 1, nproma
          xtm1(jl,jk,jt,jrow) = xtm1(jl,jk,jt,jrow) + zeps*xt(jl,jk,jt,jrow)
        END DO
      END DO
    END DO

    DO jl = 1, nproma
      alpsm1(jl,jrow) = alpsm1(jl,jrow) + zeps*alps(jl,jrow)
    END DO

  END DO
!$OMP END DO
!$OMP END PARALLEL
!CSD$ END PARALLEL DO

  IF (ltdiag) THEN

    DO jrow = 1, ngpblks
      IF (.NOT. lset_fh1) THEN
        lset_fh1 = DOT_PRODUCT(ptfh2(:,jrow),ptfh2(:,jrow)) > EPSILON(1.0_dp)
      END IF
      IF (lset_fh1) THEN
        ptfh1(1:nglon,:,1,jrow) = vom1  (:,:,jrow) - ptfh1(1:nglon,:,1,jrow)
        ptfh1(1:nglon,:,2,jrow) = dm1   (:,:,jrow) - ptfh1(1:nglon,:,2,jrow)
        ptfh1(1:nglon,:,3,jrow) = tm1   (:,:,jrow) - ptfh1(1:nglon,:,3,jrow)
        ptfh2(1:nglon    ,jrow) = alpsm1(:,jrow)   - ptfh2(1:nglon    ,jrow)
      ELSE
        ptfh1(:,:,:,jrow) = 0.0_dp
        ptfh2(:    ,jrow) = 0.0_dp
      ENDIF
      ! accumulate time filter part
      pdiga(1:nglon,:,21:23,jrow) = pdiga(1:nglon,:,21:23,jrow) + ptfh1(1:nglon,:,1:3,jrow)
      pdsga(1:nglon,      2,jrow) = pdsga(1:nglon,      2,jrow) + ptfh2(1:nglon      ,jrow)
    ENDDO

  ENDIF

END SUBROUTINE tf2
