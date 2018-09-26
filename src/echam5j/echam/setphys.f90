SUBROUTINE setphys

  !
  ! preset physical constants and process according namelist 
  !
  ! M. Jarraud, ECMWF, December 1982, original version
  ! A. Tompkins, MPI, April 2000, added LUTs for cover
  ! L. Kornblueh, MPI, April 2003, moved calculation of convection parameterization LUTs to mo_convect_tables
  ! U. Schlese, M.Esch June 2007, added volcanic aerosol
  !
  ! This subroutine preset ,modifies and derive
  ! constants in modules used in the parameterization
  ! of diabatic processes except the constants internal
  ! to the radiation "black box".
  !

  USE mo_kind,           ONLY: dp 
  USE mo_mpi,            ONLY: p_io, p_parallel, p_parallel_io, p_bcast  
  USE mo_doctor,         ONLY: nerr, nout
  USE mo_control,        ONLY: nn, l_volc, ngl
  USE mo_constants,      ONLY: api
  USE mo_param_switches, ONLY: lphys, lrad, lvdiff, lcond, lsurf, lcover, lconv, lgwdrag, lice, iconv
  USE mo_physc1,         ONLY: sinrad, cosrad
  USE mo_cumulus_flux,   ONLY: lmfpen, lmfscv, lmfmid, lmfdd, lmfdudv
  USE mo_convect_tables, ONLY: init_convect_tables  
  USE mo_cloud,          ONLY: rbetak, cbeta_pq_max, cbeta_pq, nbetax, nbetaq, cbetaqs, tbetai0, tbetai1
  USE mo_namelist,       ONLY: position_nml, nnml, POSITIONED
  USE mo_exception,      ONLY: finish
  USE mo_volc_data,      ONLY: jpd, aod, reff, extrat, ssa, asym, init_volc_tables, read_volc_data

  IMPLICIT NONE

  INCLUDE 'physctl.inc'
  
  INTEGER, PARAMETER ::   ilrad = 128     ! this variable is independent from nlon of T42 !!!    

  REAL(dp) :: zl, zq, zx

  INTEGER :: jlon, it, iq
  INTEGER :: ierr                         ! error return value from position_nml

  REAL(dp), EXTERNAL :: betai

  !     ------------------------------------------------------------
  !
  !*        1.       preset constants.
  !                  ------ ----------
  !
100 CONTINUE
  !
  !*        1.1      preset constants in MO_PHYSC1
  !
110 CONTINUE
  !
  DO jlon = 1, ilrad
    zl = 2.0_dp*api*(jlon-1.0_dp)/ilrad
    sinrad(jlon) = SIN(zl)
    cosrad(jlon) = COS(zl)
  END DO
  !
  !*        1.2      preset constants in MO_PARAM_SWITCHES.
  !
120 CONTINUE
  !
  lphys   = .TRUE.
  lrad    = .TRUE.
  lvdiff  = .TRUE.
  lconv   = .TRUE.
  iconv   = 1
  lcond   = .TRUE.
  lsurf   = .TRUE.
  lcover  = .TRUE.
  lmfpen  = .TRUE.
  lmfscv  = .TRUE.
  lmfmid  = .TRUE.
  lmfdd   = .TRUE.
  lmfdudv = .TRUE.
  IF (nn == 21) THEN
     lgwdrag = .FALSE.
  ELSE
     lgwdrag = .TRUE.
  ENDIF
  lice    = .TRUE.
  !
  !
  !*         1.3     Initialise lookup tables for CUADJTQ
  !
130 CONTINUE
  !
  CALL init_convect_tables
  !
  !*         1.4     Initialise lookup tables for COVER
  !
140 CONTINUE
  !
  rbetak = (cbeta_pq_max-cbeta_pq)/(EXP(cbetaqs)-1.0_dp)
  DO it = 0, nbetax 
    DO iq = 0, nbetaq
      zx = REAL(it,dp)/REAL(nbetax,dp)
      zq = rbetak*(EXP(cbetaqs*REAL(iq,dp)/REAL(nbetaq,dp))-1.0_dp)+cbeta_pq
      tbetai0(iq,it) = betai(cbeta_pq       ,zq,zx)
      tbetai1(iq,it) = betai(cbeta_pq+1.0_dp,zq,zx)
    ENDDO
  ENDDO
  !
  !*         1.5     Initialise lookup tables for aerosol radiation parameters
  !
150 CONTINUE
  !
  IF(l_volc) THEN
    ALLOCATE ( aod(ngl,0:jpd+1))
    ALLOCATE (reff(ngl,0:jpd+1))
    IF (p_parallel_io) THEN
       WRITE(nout,*) 'l_volc =.TRUE.  --> volcanic forcing on'
       CALL init_volc_tables
    ENDIF
    IF (p_parallel) THEN
       CALL p_bcast (extrat, p_io)
       CALL p_bcast (ssa, p_io)
       CALL p_bcast (asym, p_io)
    ENDIF
    IF (p_parallel_io) THEN
       CALL read_volc_data
    ENDIF
    IF (p_parallel) THEN
       CALL p_bcast (aod, p_io)
       CALL p_bcast (reff, p_io)
    ENDIF
  ELSE
       IF (p_parallel_io) &
           WRITE(nout,*) 'l_volc =.FALSE.  --> volcanic forcing off'
  ENDIF
  !
  !     ------------------------------------------------------------
  !
  !*        2.       READ NAMELIST.
  !                  ---- ---------
  !
200 CONTINUE
  !
  IF (p_parallel_io) THEN
     CALL position_nml ('PHYSCTL', status=ierr)
     SELECT CASE (ierr)
     CASE (POSITIONED)
       READ (nnml, physctl)
     END SELECT
  ENDIF
  IF (p_parallel) THEN
     CALL p_bcast (lphys, p_io)
     CALL p_bcast (lrad, p_io)
     CALL p_bcast (lvdiff, p_io)
     CALL p_bcast (lcond, p_io)
     CALL p_bcast (lsurf, p_io)
     CALL p_bcast (lcover, p_io)
     CALL p_bcast (lconv, p_io)
     CALL p_bcast (iconv, p_io)
     CALL p_bcast (lmfpen, p_io)
     CALL p_bcast (lgwdrag, p_io)
     CALL p_bcast (lice, p_io)
  ENDIF
!
!     ------------------------------------------------------------
!
!*        3.       MODIFY CONSTANTS.
!                  ------ ----------
!
300 CONTINUE
!
!
!*        3.1      MODIFY CONSTANTS IN *mo_param_switches*.
!
310 CONTINUE
  IF(.NOT.lphys) THEN
     lrad=.FALSE.
     lvdiff=.FALSE.
     lgwdrag=.FALSE.
     lconv=.FALSE.
     lcond=.FALSE.
     lsurf=.FALSE.
     lcover=.FALSE.
     lice=.FALSE.
  END IF
!
  IF(.NOT.lconv) THEN
     lmfpen=.FALSE.
     lmfscv=.FALSE.
     lmfmid=.FALSE.
     lmfdd=.FALSE.
     lmfdudv=.FALSE.
  ELSE
     ! Check iconv
     ! -----------
     SELECT CASE (iconv)
     CASE(1)
       IF (p_parallel_io) &
            WRITE(nout,*) 'iconv = 1 --> Convection: Nordeng (default)'
     CASE(2)
       IF (p_parallel_io) &
            WRITE(nout,*) 'iconv = 2 --> Convection: Tiedtke'
     CASE(3)
       IF (p_parallel_io) &
            WRITE(nout,*) 'iconv = 3 --> Convection: Hybrid'
     CASE default
       IF (p_parallel_io) &
            WRITE(nerr,*) 'iconv =',iconv,' in physctl namelist not supported'
       CALL finish('setphys','Run terminated')
     END SELECT
  ENDIF
!
  IF(.NOT.lmfpen) THEN
     lconv=.FALSE.
     lmfscv=.FALSE.
     lmfmid=.FALSE.
     lmfdd=.FALSE.
     lmfdudv=.FALSE.
  END IF
!
!
!*        3.2      SET UP CONSTANTS IN *mo_physc2*.
!
320 CONTINUE
!
  CALL iniphy
!
!     ------------------------------------------------------------
!
!*         5.     WRITE NAMELIST.
!                 ----- ---------
!
500 CONTINUE
  IF (.NOT. p_parallel) THEN
    WRITE(nout,physctl)
  END IF
!
!     ------------------------------------------------------------
!
  RETURN
END SUBROUTINE setphys

!-----------------------------------------------------------------------

FUNCTION betai(p,q,x) 

  ! Description:
  !
  ! Uses betacf, gammln; returns the incomplete beta function I x (a; b).  
  !
  ! Method:
  !   See Numerical Recipes (Fortran)
  !
  ! Author:
  !   A. Tompkins, MPI, 2000
  !

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  REAL(dp)             :: betai
  REAL(dp), INTENT(in) :: p,q, & ! beta shape parameters
                          x      ! integration limit
  !  local scalars: 
  REAL(dp) :: bt,betacf,gammln 

  IF (x > 0.0_dp .AND. x < 1.0_dp ) THEN  ! factors in front of the continued fraction. 
    bt = EXP(gammln(p+q)-gammln(p)-gammln(q)+p*LOG(x)+q*LOG(1.0_dp-x)) 
  ELSE  
    bt = 0.0_dp 
  ENDIF 
  IF (x < (p+1.0_dp)/(p+q+2.0_dp)) THEN ! use continued fraction directly.
    betai = bt*betacf(p,q,x)/p 
  ELSE     ! use continued fraction after making the symmetry transformation. 
    betai = 1.0_dp-bt*betacf(q,p,1.0_dp-x)/q !
  ENDIF 

END FUNCTION betai

!-----------------------------------------------------------------------

FUNCTION betacf(p,q,x) 

  ! Description:
  !
  !  used by betai: evaluates continued fraction for incomplete 
  !  beta function by modi ed lentz's method ( x 5.2). 
  !  first step of lentz's method. 
  !
  ! Method:
  !   See Numerical Recipes (Fortran)
  !
  ! Author:
  !   A. Tompkins, MPI, 2000
  !

  USE mo_kind,      ONLY: dp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE

  REAL(dp)             :: betacf
  REAL(dp), INTENT(in) :: p,q, & ! beta shape parameters
                          x      ! integration limit

  INTEGER :: maxit = 100, m, m2 
  REAL(dp) :: zeps = 3.e-7_dp, fpmin = 1.e-30_dp, &
              aa, c, d, del, h, qab, qam, qap 

  qab = p+q 

  ! these q's will be used in factors that occur in the coe cients (6.4.6). 

  qap = p+1.0_dp 
  qam = p-1.0_dp 
  c = 1.0_dp 
  d = 1.0_dp-qab*x/qap 
  IF (ABS(d) < fpmin) d = fpmin 
  d = 1.0_dp/d 
  h = d 
  m = 1
  del = 2.0_dp 
  DO WHILE (ABS(del-1.0_dp) > zeps)
    m2 = 2*m 
    aa = m*(q-m)*x/((qam+m2)*(p+m2))
    d = 1.0_dp+aa*d  ! one step (the even one) of the recurrence. 
    IF (ABS(d) < fpmin) d = fpmin 
    c = 1.0_dp+aa/c 
    IF (ABS(c) < fpmin) c = fpmin 
    d = 1.0_dp/d 
    h = h*d*c 
    aa = -(p+m)*(qab+m)*x/((p+m2)*(qap+m2)) 
    d = 1.0_dp+aa*d ! next step of the recurrence (the odd one). 
    IF (ABS(d) < fpmin) d = fpmin 
    c = 1.0_dp+aa/c 
    IF (ABS(c) < fpmin) c = fpmin 
    d = 1.0_dp/d 
    del = d*c 
    h = h*del 
    m = m+1            ! AMT
    IF (m > maxit) THEN 
      CALL finish ('betacf','a or b too big, or maxit too small in betacf') 
      del = 1.0_dp
    ENDIF  
  ENDDO 
  betacf = h 

END FUNCTION betacf

!-----------------------------------------------------------------------

FUNCTION gammln(xx) 

  ! Description:
  !
  ! Gamma function calculation
  ! returns the value ln[g(xx)] for xx > 0.
  !
  ! Method:
  !   See Numerical Recipes
  !
  ! Author:
  !   A. Tompkins, MPI, 2000
  !
  
  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  REAL(dp)             :: gammln
  REAL(dp), INTENT(in) :: xx            

  INTEGER :: j 

  REAL(dp) :: ser, tmp, x, y
  REAL(dp), PARAMETER :: cof(6) = (/ &
       76.18009172947146_dp, -86.50532032941677_dp, &
       24.01409824083091_dp, -1.231739572450155_dp, &
       0.1208650973866179e-2_dp, -0.5395239384953e-5_dp /)
  REAL(dp), PARAMETER :: stp = 2.5066282746310005_dp 

  x = xx 
  y = x 
  tmp = x+5.5_dp 
  tmp = (x+0.5_dp)*LOG(tmp)-tmp 
  ser = 1.000000000190015_dp 
  DO j =1, 6
    y = y+1.0_dp 
    ser = ser+cof(j)/y 
  ENDDO 
  gammln = tmp+LOG(stp*ser/x) 

END FUNCTION gammln
