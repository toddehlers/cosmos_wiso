SUBROUTINE iniphy

  ! Description:
  !
  ! Initialises physical constants of uncertain value.
  !
  ! Method:
  !
  ! This routine sets the values for the physical constants used
  ! in the parameterization routines (except for the radiation
  ! black-box) whenever these values are not well enough known to
  ! forbid any tuning or whenever they are subject to an arbitrary
  ! choice of the modeller. These constants will be in *comph2*.
  !
  ! *iniphy* is called from *physc*.
  !
  ! Authors:
  !
  ! J. F. Geleyn, ECMWF, December 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! M. Esch, MPI, June 1999, ECHAM5-modifications
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,         ONLY: dp
  USE mo_control,      ONLY: nlev
  USE mo_constants,    ONLY: g
  USE mo_physc2,       ONLY: clam, ckap, cb, cc, cd, cchar, cfreec     &
                           , cgam, cvdifts, ctfreez, cz0ice, cqsncr    &
                           , csncri, cevapcu, cdel, cmid, clice        &
                           , cqcon, cqdif, cgh2o, cwlmax, corvari      &
                           , corvars
  USE mo_hyb,          ONLY: ceta
  USE mo_vegetation,   ONLY: cvroots, cvrootd, cvrootc, cvrad          &
                           , cvinter, cva, cvb, cvc, cvk, cvbc         &
                           , cvkc, cvabc
  USE mo_cumulus_flux, ONLY: cuparam
  USE mo_cloud,        ONLY: sucloud

  IMPLICIT NONE

  !  Local scalars: 
  INTEGER :: jk

  !  Intrinsic functions 
  INTRINSIC SQRT


  !  Executable statements 

!-- 1. Setting of constants

!-- 1.1 Constants for *vdiff*

  clam = 150._dp
  ckap = 0.4_dp
  cb = 5._dp
  cc = 5._dp
  cd = 5._dp
  cchar = 0.018_dp
  cfreec = 0.001_dp
  cgam = 1.25_dp
  cvdifts = 1.5_dp

!-- 1.2 Constants for *vdiff*, *clsst* and *atmice*

  ctfreez = 271.38_dp
  cz0ice = 0.001_dp

!-- 1.3 Constants for *vdiff* and *surf*

  cqsncr = 0.95_dp
  csncri = 5.85036E-3_dp

!-- 1.4 Constant for massflux convection scheme

  CALL cuparam

!-- 1.5 Highest level ncctop where condensation is allowed

  CALL sucloud
!
  DO jk = 1, nlev
   cevapcu(jk) = 1.93E-6_dp*261._dp*SQRT(1.E3_dp/(38.3_dp*0.293_dp)*   &
                                              SQRT(ceta(jk)))*0.5_dp/g
  END DO

!-- 1.6 Constants for *surf*
!
!  THICKNESS OF SOIL LAYERS
!
      CDEL(1)=0.065_dp
      CDEL(2)=0.254_dp
      CDEL(3)=0.913_dp
      CDEL(4)=2.902_dp
      CDEL(5)=5.700_dp
!
!  DEPTH OF MIDS OF SOIL LAYERS
!
      CMID(1)=CDEL(1)*0.5_dp
      CMID(2)=CDEL(1)+CDEL(2)*0.5_dp
      CMID(3)=CDEL(1)+CDEL(2)+CDEL(3)*0.5_dp
      CMID(4)=CDEL(1)+CDEL(2)+CDEL(3)+CDEL(4)*0.5_dp
      CMID(5)=CDEL(1)+CDEL(2)+CDEL(3)+CDEL(4)+CDEL(5)*0.5_dp
!
  clice = 0.3_dp
  cqcon = 1.E-10_dp
  cqdif = 1.E-7_dp


!-- 1.7 Constants for *vdiff* and *surf*

  cgh2o = 4.18E6_dp

  cwlmax = 2.E-4_dp

  cvroots = .50_dp
  cvrootd = .50_dp
  cvrootc = 0._dp
  cvrad = 0.55_dp
  cvinter = 0.25_dp
  corvari = 50._dp**2
  corvars = 750._dp**2

  cva = 5000._dp
  cvb = 10._dp
  cvc = 100._dp
  cvk = .9_dp
  cvbc = cvb*cvc
  cvkc = cvk*cvc
  cvabc = (cva+cvbc)/cvc

  RETURN
END SUBROUTINE iniphy
