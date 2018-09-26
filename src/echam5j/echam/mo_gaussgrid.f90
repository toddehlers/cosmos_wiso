MODULE mo_gaussgrid

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  ! ---------------------------------------------------------------
  !
  ! module *mo_gaussgrid* - quantities related to the gaussian grid.
  !
  ! ---------------------------------------------------------------

  REAL(dp), ALLOCATABLE :: gl_gw(:)       ! Gaussian weights
  REAL(dp), ALLOCATABLE :: gl_gmu(:)      ! mu = sin(Gaussian latitudes)
  REAL(dp), ALLOCATABLE :: gl_coriol(:)   ! coriolis parameter, 2*omega*mu
  REAL(dp), ALLOCATABLE :: gl_twomu(:)    ! 2*mu
  REAL(dp), ALLOCATABLE :: gl_cst(:)      ! square of cos(latitude):(1-mu**2)
  REAL(dp), ALLOCATABLE :: gl_sqcst(:)    ! sqrt(1-mu**2)
  REAL(dp), ALLOCATABLE :: gl_rcsth(:)    ! half reciprocal of *cst*:1/(2*cst)
  REAL(dp), ALLOCATABLE :: gl_racst(:)    ! 1/(a*(1-mu**2)
  REAL(dp), ALLOCATABLE :: gl_rsqcst(:)   ! 1/sqrt(1-mu**2)
  REAL(dp), ALLOCATABLE :: gl_budw(:)     ! weights for global budgets,
                                          ! budw = gw/nlon
  REAL(dp), ALLOCATABLE :: gridarea(:)    ! area of a grid cell [m**2]
  REAL(dp), ALLOCATABLE :: philat(:)      ! gaussian latitude   [degrees]
  REAL(dp), ALLOCATABLE :: philon(:)      ! longitudes          [degrees]
  REAL(dp), ALLOCATABLE :: sinlon(:)      ! sin(longitude).
  REAL(dp), ALLOCATABLE :: coslon(:)      ! cos(longitude).

CONTAINS

  SUBROUTINE inigau

    ! Description:
    !
    ! Preset constants in *mo_gaussgrid*.
    !
    ! Method:
    !
    ! *inigau* is called from *setdyn*.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, December 1982, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A. Rhodin, MPI, January 1999, subroutine inigau -> module mo_gaussgrid
    ! U. Schlese, MPI, January 2001, grid area added
    ! A. Rhodin, MPI, June 2001, philon, philat added
    ! L. Kornblueh, MPI, October 2001, provide non 'ping-pong' properties
    ! U. Schulzweida, MPI, May 2002, change 'ping-pong' to N->S
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_control,   ONLY: ngl, nhgl, nlon
    USE mo_constants, ONLY: a, api, omega

    IMPLICIT NONE

    !  Local scalars: 
    REAL(dp) :: zcst, zl, zsqcst
    INTEGER :: jgl, jlon

    !  Local arrays: 
    REAL(dp) :: zgmu(ngl), zgw(ngl)

    !  Intrinsic functions 
    INTRINSIC COS, SIN, SQRT


    !  Executable statements 

    !-- 0. Allocate module provided fields

    IF (.NOT. ALLOCATED(gl_gw)) THEN
      ALLOCATE(gl_gw(ngl))      
      ALLOCATE(gl_gmu(ngl))     
      ALLOCATE(gl_coriol(ngl))  
      ALLOCATE(gl_twomu(ngl))   
      ALLOCATE(gl_cst(ngl))     
      ALLOCATE(gl_sqcst(ngl))   
      ALLOCATE(gl_rcsth(ngl))   
      ALLOCATE(gl_racst(ngl))   
      ALLOCATE(gl_rsqcst(ngl))  
      ALLOCATE(gl_budw(ngl))    
      ALLOCATE(gridarea(ngl))    
      ALLOCATE(philat(ngl))    
      ALLOCATE(philon(nlon))    
      ALLOCATE(sinlon(2*nlon))    
      ALLOCATE(coslon(2*nlon))    
    END IF

    !-- 1. Compute Gaussian latitudes and weights

    CALL gauaw(zgmu,zgw,ngl)

    DO jgl = 1, nhgl
      gl_gw(jgl)        = zgw(jgl)*0.5_dp
      gl_gw(ngl-jgl+1)  = zgw(jgl)*0.5_dp
      gl_gmu(jgl)       = zgmu(jgl)
      gl_gmu(ngl-jgl+1) = zgmu(jgl)
      
      !-- 2. Derive some other constants
      
      gl_coriol(jgl)       =  2*omega*zgmu(jgl)
      gl_coriol(ngl-jgl+1) = -2*omega*zgmu(jgl)
      gl_twomu(jgl)        =  2*zgmu(jgl)
      gl_twomu(ngl-jgl+1)  = -2*zgmu(jgl)
      gl_budw(jgl)         = gl_gw(jgl)/nlon
      gl_budw(ngl-jgl+1)   = gl_gw(jgl)/nlon
      zcst                 = 1.0_dp-zgmu(jgl)**2
      zsqcst               = SQRT(zcst)
      gl_cst(jgl)          = zcst
      gl_cst(ngl-jgl+1)    = zcst
      gl_sqcst(jgl)        = zsqcst
      gl_sqcst(ngl-jgl+1)  = zsqcst
      gl_rsqcst(jgl)       = 1.0_dp/zsqcst
      gl_rsqcst(ngl-jgl+1) = 1.0_dp/zsqcst
      gl_rcsth(jgl)        = 0.5_dp/zcst
      gl_rcsth(ngl-jgl+1)  = 0.5_dp/zcst
      gl_racst(jgl)        = 1.0_dp/(a*zcst)
      gl_racst(ngl-jgl+1)  = 1.0_dp/(a*zcst)
    END DO


    DO jlon = 1, nlon               ! double size for rotated domains 
      zl = 2*api*(jlon-1.0_dp)/nlon    ! on decomposed grid
      sinlon(jlon) = SIN(zl)
      sinlon(jlon+nlon) = sinlon(jlon)
      coslon(jlon) = COS(zl)
      coslon(jlon+nlon) = coslon(jlon)
      philon(jlon) = 360._dp*(jlon-1.0_dp)/nlon
    END DO

    !  Grid area stored from N - > S !

    DO jgl=1,nhgl      
      gridarea(jgl)       = gl_budw(jgl)*4*api*a**2
      gridarea(ngl+1-jgl) = gridarea(jgl)
      philat  (jgl)       = 180._dp/api*ASIN(zgmu(jgl))
      philat  (ngl-jgl+1) = -philat(jgl)
    END DO

  END SUBROUTINE inigau
!------------------------------------------------------------------------------
  SUBROUTINE gauaw (pa, pw, nlat)

    ! Description:
    !
    ! Compute abscissas and weights for gaussian integration.
    !
    ! Method:
    !

    USE mo_constants, ONLY: api

    IMPLICIT NONE

    !  Scalar arguments 
    INTEGER :: nlat

    !  Array arguments 
    REAL(dp) :: pa(nlat), pw(nlat)
    ! *pa*  - array, length at least *k,* to receive abscis abscissas.
    ! *pw*  - array, length at least *k,* to receive weights.


    !  Local scalars: 
    REAL(dp), PARAMETER :: epsil = EPSILON(0.0_dp)
    INTEGER, PARAMETER :: itemax = 20

    INTEGER :: iter, ins2, isym, jn, jgl
    REAL(dp):: za, zw, z, zan
    REAL(dp):: zk, zkm1, zkm2, zx, zxn, zldn, zmod

    !  Intrinsic functions 
    INTRINSIC ABS, COS, MOD, TAN

    !  Executable statements 

    ins2 = nlat/2+MOD(nlat,2)

    ! Find first approximation of the roots of the
    ! Legendre polynomial of degree nlat
    
    DO jgl = 1, ins2
       z = REAL(4*jgl-1,dp)*api/REAL(4*nlat+2,dp)
       pa(jgl) = COS(z+1.0_dp/(TAN(z)*REAL(8*nlat**2,dp)))
    END DO

    ! Computes roots and weights
    ! Perform the Newton loop
    ! Find 0 of Legendre polynomial with Newton loop

    DO jgl = 1, ins2

       za = pa(jgl)
    
       DO iter = 1, itemax+1
          zk = 0.0_dp

          ! Newton iteration step
    
          zkm2 = 1.0_dp
          zkm1 = za
          zx = za
          DO jn = 2, nlat
             zk = (REAL(2*jn-1,dp)*zx*zkm1-REAL(jn-1,dp)*zkm2)/REAL(jn,dp)
             zkm2 = zkm1
             zkm1 = zk
          END DO
          zkm1 = zkm2
          zldn = (REAL(nlat)*(zkm1-zx*zk))/(1.0_dp-zx*zx)
          zmod = -zk/zldn
          zxn = zx+zmod
          zan = zxn
    
          ! computes weight
    
          zkm2 = 1.0_dp
          zkm1 = zxn
          zx = zxn
          DO jn = 2,nlat
             zk = (REAL(2*jn-1,dp)*zx*zkm1-REAL(jn-1,dp)*zkm2)/REAL(jn,dp)
             zkm2 = zkm1
             zkm1 = zk
          END DO
          zkm1 = zkm2
          zw = (1.0_dp-zx*zx)/(REAL(nlat*nlat,dp)*zkm1*zkm1)
          za = zan
          IF (ABS(zmod) <= epsil) EXIT
       END DO

       pa(jgl) = zan
       pw(jgl) = 2*zw
    
    ENDDO

!DIR$ IVDEP
!OCL NOVREC

    DO jgl = 1, nlat/2
       isym = nlat-jgl+1
       pa(isym) = -pa(jgl)
       pw(isym) = pw(jgl)
    ENDDO

  END SUBROUTINE gauaw
!------------------------------------------------------------------------------
  SUBROUTINE cleanup_gaussgrid
    IF (ALLOCATED(gl_gw)) THEN
      DEALLOCATE (gl_gw)
      DEALLOCATE (gl_gmu)
      DEALLOCATE (gl_coriol)
      DEALLOCATE (gl_twomu)
      DEALLOCATE (gl_cst)
      DEALLOCATE (gl_sqcst)
      DEALLOCATE (gl_rcsth)
      DEALLOCATE (gl_racst)
      DEALLOCATE (gl_rsqcst)
      DEALLOCATE (gl_budw)
      DEALLOCATE (gridarea)
      DEALLOCATE (philat)
      DEALLOCATE (philon)
      DEALLOCATE (sinlon)
      DEALLOCATE (coslon)
    END IF
  END SUBROUTINE cleanup_gaussgrid
!------------------------------------------------------------------------------
END MODULE mo_gaussgrid
