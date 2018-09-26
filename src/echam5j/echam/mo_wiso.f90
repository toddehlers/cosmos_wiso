MODULE mo_wiso
  !--------------------------------------------------------------
  ! Stable Water Isotopes (18-O, HDO) in the Hydrological Cycle
  !
  ! Echam5 Code based on previous Echam4_iso Model
  !
  ! Martin Werner, Marcus Herold, Petra Langebroek
  ! AWI Bremerhaven, 2009
  !-------------------------------------------------------------- 

  USE mo_kind,            ONLY: dp
  USE mo_constants,       ONLY: tmelt
  USE mo_cloud,           ONLY: cthomi

  IMPLICIT NONE

  LOGICAL :: lwiso                 ! switch water isotope diagnostics on/off
  LOGICAL :: lwiso_rerun           ! switch for initial/resumed water isotope values
  INTEGER, PARAMETER :: nwisomax=4 ! maximum number of water isotopes  
  INTEGER :: nwiso                 ! number of water isotope tracers
  INTEGER :: nwisotyp(nwisomax)    ! number type of the water isotope tracer: 1=h218o 2=hdo 3=h2o  

  CHARACTER(LEN=7), DIMENSION(nwisomax) :: wisoq_names,   wisoxl_names,   wisoxi_names,   &
                                           wisoqm1_names, wisoxlm1_names, wisoxim1_names

  ! water isotope properties
  REAL(dp) :: tnat(nwisomax)            ! natural isotope-ratio
  REAL(dp) :: twisorhoh2o(nwisomax)     ! density of isotopic water
  REAL(dp) :: toce(nwisomax)            ! mean ocean concentration of the tracer
  REAL(dp) :: twisoatm(nwisomax)        ! initial isotope deviation from SMOW in the atmosphere
  REAL(dp) :: tdifrel(nwisomax)         ! relation of the diffusivities rel. water
  REAL(dp) :: twisosur1(nwisomax)       ! factors for initial deviation from SMOW in the surf. layers
  REAL(dp) :: twisosur2(nwisomax)       ! factors for initial deviation from SMOW in the surf. layers
  REAL(dp) :: tkinsl(nwisomax)          ! factors for the kinetic fractionation over ocean
  REAL(dp) :: tkinfa1(nwisomax)         ! factors for the kinetic fractionation over ocean
  REAL(dp) :: tkinfa2(nwisomax)         ! factors for the kinetic fractionation over ocean
  REAL(dp) :: tkinfa3(nwisomax)         ! factors for the kinetic fractionation over ocean
  REAL(dp) :: talphal1(nwisomax)        ! factors for the computation of the fractionation over liquid
  REAL(dp) :: talphal2(nwisomax)        ! factors for the computation of the fractionation over liquid
  REAL(dp) :: talphal3(nwisomax)        ! factors for the computation of the fractionation over liquid
  REAL(dp) :: talphal4(nwisomax)        ! factors for the computation of the fractionation over liquid
  REAL(dp) :: talphas1(nwisomax)        ! factors for the computation of the fractionation over ice
  REAL(dp) :: talphas2(nwisomax)        ! factors for the computation of the fractionation over ice 
  REAL(dp) :: talphas3(nwisomax)        ! factors for the computation of the fractionation over ice 
  REAL(dp) :: talphas4(nwisomax)        ! factors for the computation of the fractionation over liquid

! default saturation values
  REAL(dp), PARAMETER :: tsatbase=1.01_dp  ! base value for temperatur-dependency of supersaturation during kinetic effects
  REAL(dp), PARAMETER :: tsatfac=0.0045_dp ! temperatur-dependency of supersaturation during kinetic effects
! saturation values for new fractionation factors by M. Ellehoj
!  REAL(dp), PARAMETER :: tsatbase=1.0_dp   ! base value for temperatur-dependency of supersaturation during kinetic effects
!  REAL(dp), PARAMETER :: tsatfac=0.002_dp  ! temperatur-dependency of supersaturation during kinetic effects

  REAL(dp), PARAMETER :: tdifexp=0.58_dp   ! a constant for alpha_eff for equilibrium below cloud base
  REAL(dp), PARAMETER :: thumwiso1=0.75_dp ! coefficient for effective fractionation below cumulus clouds
  REAL(dp), PARAMETER :: thumwiso2=0.25_dp ! coefficient for effective fractionation below cumulus clouds
  REAL(dp), PARAMETER :: twisoeqcu=0.5_dp  ! coefficient for raindrop equilibrium below cumulus clouds
  REAL(dp), PARAMETER :: twisoeqls=0.9_dp  ! coefficient for raindrop equilibrium below large scale clouds
  
  REAL(dp), PARAMETER :: twisoice=cthomi   ! temperature limit for homogenous ice cumulus clouds
                                           ! (set to same default Echam5 value (-35C) as for large-scale clouds)
                                           
  REAL(dp), PARAMETER :: snglacmn=0.001_dp ! minimum threshold value for the snow layer depth on glaciers 
  REAL(dp), PARAMETER :: snglacmx=0.02_dp  ! maximum threshold value for the snow layer depth on glaciers 

  INTEGER,  PARAMETER :: lwisofracl=0      ! switch for evapotranspiration over land:
                                           ! 0: no fractionation over land surface
                                           ! 1: fractionation for soil evaporation, only
                                           ! 2: fractionation for soil evaporation and plant transpiration
                                                 
  INTEGER,  PARAMETER :: lwisokinl=0       ! switch for kinetic fractionation over land:
                                           ! 0: no kinetic fractionation over land surface
                                           ! 1: kinetic fractionation according to Brutsaert (1975) as for frac. over water
                                           ! 2: kinetic fractionation according to Mathieu & Bariac (1996)
                                                 
  REAL(dp), PARAMETER :: cwisomin= 1.e-15_dp  ! mininum value for calculating different water isotopic ratios
  REAL(dp), PARAMETER :: cwisosec= 1.e-12_dp  ! security parameter for calculation of delta vaues
  
  ! prescribed ocean surface water isotopes (delta values)
  INTEGER :: nisw  = 25   ! *nisw*      logical file unit for surface seawater isotopes file


  ! Define namelist
  ! ---------------
  NAMELIST /wisoctl/ lwiso, lwiso_rerun, nwiso

CONTAINS


  SUBROUTINE wiso_frac_liq(kproma,kbdim,kwiso,pt,pwisofracliq)
  ! ---------------------------------------------------
  !
  ! wiso_frac_liq calculates fractionation coefficient for vapour-liquid
  ! for a band of longitudinal grid points, simultaneously 
  !
  ! ---------------------------------------------------
     
    IMPLICIT NONE

  ! input arguments
    INTEGER, INTENT(IN)     :: kproma, kbdim, kwiso
    REAL(dp), INTENT(IN)    :: pt(kbdim)
  
  ! input/output arguments
    REAL(dp), INTENT(INOUT) :: pwisofracliq(kbdim,kwiso)
  
  ! local variables
    INTEGER     :: jl,jt
    
  ! fractionation over liquid water
    DO jt=1,kwiso
       DO jl=1,kproma
          pwisofracliq(jl,jt) = (exp(talphal1(jt)/(pt(jl)**2._dp)+talphal2(jt)/pt(jl)+talphal3(jt)))**talphal4(jt)
       END DO
    END DO
  
    RETURN
  END SUBROUTINE wiso_frac_liq


  SUBROUTINE wiso_frac_liq_ice(kproma,kbdim,kwiso,pt,pwisofracliq,pwisofracice)
  ! ---------------------------------------------------
  !
  ! wiso_frac_ice calculates fractionation coefficients for vapour-liquid and vapour-ice
  ! for a band of longitudinal grid points, simultaneously 
  !
  ! ---------------------------------------------------
     
    USE mo_constants,     ONLY: tmelt
  
    IMPLICIT NONE

  ! input arguments
    INTEGER, INTENT(IN)     :: kproma, kbdim, kwiso
    REAL(dp), INTENT(IN)    :: pt(kbdim)
  
  ! input/output arguments
    REAL(dp), INTENT(INOUT) :: pwisofracliq(kbdim,kwiso),pwisofracice(kbdim,kwiso)
  
  ! local variables
    INTEGER     :: jl,jt
    REAL(dp)    :: zsatval
    
  ! fractionation over liquid water
    DO jt=1,kwiso
       DO jl=1,kproma
         pwisofracliq(jl,jt) = (exp(talphal1(jt)/(pt(jl)**2._dp)+talphal2(jt)/pt(jl)+talphal3(jt)))**talphal4(jt)
       END DO
    END DO
  
  ! fractionation over ice
    DO jt=1,kwiso
       DO jl=1,kproma
         pwisofracice(jl,jt) = (exp(talphas1(jt)/(pt(jl)**2._dp)+talphas2(jt)/pt(jl)+talphas3(jt)))**talphas4(jt)
         IF (nwisotyp(jt).ne.1.and.pt(jl).lt.tmelt) THEN ! effective fractionation over ice if necessary
            zsatval=tsatbase-tsatfac*(pt(jl)-tmelt)
            pwisofracice(jl,jt)=pwisofracice(jl,jt)*(zsatval/(1._dp+pwisofracice(jl,jt)*(zsatval-1._dp)*tdifrel(jt)))
         ENDIF
       END DO
    END DO
   
    RETURN
  END SUBROUTINE wiso_frac_liq_ice


  FUNCTION wiso_calc_delta(jt,pwiso,pdefault)
  
  IMPLICIT NONE
  
  REAL(dp), INTENT(IN) :: pwiso, pdefault
  INTEGER,  INTENT(IN) :: jt
  
  REAL(dp) :: wiso_calc_delta
  
  IF (ABS(pdefault).gt.cwisomin) THEN
     wiso_calc_delta = ((pwiso/pdefault)/tnat(jt)-1.0_dp)*1000.0_dp
  ELSE
     wiso_calc_delta = -9999.9_dp
  ENDIF
    
  END FUNCTION wiso_calc_delta
  
  
END MODULE mo_wiso
