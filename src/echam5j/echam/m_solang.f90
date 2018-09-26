MODULE m_solang

CONTAINS

SUBROUTINE solang_extended_const(p1,p2,p3,pamu0,prdayl)

  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: ldc => local_decomposition
  USE mo_gaussgrid,     ONLY: coslon, sinlon, gl_twomu, gl_sqcst
  USE mo_transpose,     ONLY: reorder

  !  Scalar arguments
  REAL(dp), INTENT(in) :: p1, p2, p3
  !  Array arguments
  REAL(dp), INTENT(out) :: pamu0(:,:), prdayl(:,:)

  INTEGER :: jlat, lnlat, lnlon, jglat, nlof
  REAL(dp)    :: zsin, ztim1, ztim2, ztim3
  REAL(dp)    :: zpamu0(ldc% nglon, ldc% nglat)
  REAL(dp)    :: zprdayl(ldc% nglon, ldc% nglat)


  lnlat = ldc%nglat
  lnlon = ldc%nglon        ! local number of longitudes

  DO jlat = 1, lnlat

    jglat = ldc%glat(jlat)
    nlof  = ldc%glon(jlat)  ! local longitude offset to global fields

    zsin  = 0.5_dp*gl_twomu(jglat)

    ztim1 =  p1*zsin
    ztim2 = -p2*gl_sqcst(jglat)
    ztim3 =  p3*gl_sqcst(jglat)

    CALL solang(lnlon, coslon(1+nlof:lnlon+nlof),                      &
                sinlon(1+nlof:lnlon+nlof),                             &
                ztim1, ztim2, ztim3,                                   &
                zpamu0(:,jlat), zprdayl(:,jlat))

  END DO

  CALL reorder (pamu0, zpamu0)
  CALL reorder (prdayl, zprdayl)

END SUBROUTINE solang_extended_const

SUBROUTINE solang(klon, coslon, sinlon, ptim1, ptim2, ptim3,           &
                  pamu0, prdayl)

  ! Description:
  !
  ! For solar zenith angle and relative daylength.
  !
  ! Method:
  !
  ! This routine gives different results depending on a logical
  ! switch. If ldiur is true one obtains actual solar zenith angles
  ! and values of one or zero depending on the sign of the former. If
  ! ldiur is false one gets the same answers at all points, i.e. mean
  ! value of the daytime solar zenith angle and relative length of
  ! the day.
  !
  ! *solang* is called from *radmod* and from *radheat*.
  ! There are three dummy arguments: *ptim1*, *ptim2* and *ptim3*
  ! are latitude dependent parameter about the sun's position.
  ! The routine returns solar zenith angles and relative day
  ! lengths to the long term storage.
  !
  ! Staightforward in the case "on". For the case "off" the
  ! type of "on" computation is repeated  with 128 points and the
  ! relevant mean values are computed and stored.
  !
  ! Authors:
  !
  ! J. F. Geleyn, ECMWF, June 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! M.A. Giorgetta, MPI, May 2000, modified for ECHAM5
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! L. Kornblueh, MPI, Januray 2003, removed MERGE
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,         ONLY: dp
  USE mo_radiation,    ONLY: ldiur
  USE mo_physc1,       ONLY: cosrad, sinrad

  IMPLICIT NONE

  ! INPUT
  ! -----

  INTEGER ,INTENT(in)               :: klon   ! number of local longitudes
  REAL(dp), INTENT(in), DIMENSION(klon) :: coslon ! cos(longitude)
  REAL(dp), INTENT(in), DIMENSION(klon) :: sinlon ! sin(longitude)
  REAL(dp), INTENT(in)                  :: ptim1  ! 
  REAL(dp), INTENT(in)                  :: ptim2  ! 
  REAL(dp), INTENT(in)                  :: ptim3  ! 

  ! OUTPUT
  ! ------

  REAL(dp) ,INTENT(out),DIMENSION(klon) :: pamu0  ! 
  REAL(dp) ,INTENT(out),DIMENSION(klon) :: prdayl ! 

  !  Local scalars: 
  REAL(dp) :: zs1, zs2, ztim1, ztim2, ztim3
  INTEGER :: jl
  LOGICAL :: lo

  !  Local arrays: 
  REAL(dp) :: zmu0(128), zrdayl(128)

  !  Intrinsic functions 
  INTRINSIC SUM

  !  Executable statements 

  ztim1 = ptim1
  ztim2 = ptim2
  ztim3 = ptim3

!-- 1. Computations if diurnal cycle "on"

  IF (ldiur) THEN
    DO jl = 1, klon
      pamu0(jl) = ztim1 + ztim2*coslon(jl) + ztim3*sinlon(jl)
      IF (pamu0(jl) >= 0.0_dp) THEN
        pamu0(jl)  = pamu0(jl)
        prdayl(jl) = 1.0_dp
      ELSE  
        pamu0(jl)  = 0.0_dp
        prdayl(jl) = 0.0_dp
      END IF

    END DO

!-- 2. Computations if diurnal cycle "off"

  ELSE
    DO jl = 1, 128
      zmu0(jl) = ztim1 + ztim2*cosrad(jl) + ztim3*sinrad(jl)
      IF (zmu0(jl) >= 0.0_dp) THEN
        zmu0(jl)   = zmu0(jl)
        zrdayl(jl) = 1.0_dp
      ELSE
        zmu0(jl)   = 0.0_dp
        zrdayl(jl) = 0.0_dp
      END IF
    END DO
    zs1 = SUM(zmu0(1:128))
    zs2 = SUM(zrdayl(1:128))
    IF (ABS(zs2) > 0._dp) THEN
      zs1 = zs1/zs2
      zs2 = zs2/128._dp
    END IF
    DO jl = 1, klon
      pamu0(jl) = zs1
      prdayl(jl) = zs2
    END DO
  END IF

END SUBROUTINE solang

END MODULE m_solang
