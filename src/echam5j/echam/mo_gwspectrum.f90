MODULE mo_gwspectrum

  !========================================================================
  !
  ! module *mo_gwspectrum* contains:
  ! A) gwsctl, the namelist for the set up of the Hines parameterization 
  !    for a gravity wave spectrum and its source spectrum. 
  ! B) Some internal switches and constants specific for the 
  !    Hines gw parameterization.
  !
  !
  ! ----------------------------------------------------------------------
  !
  ! A) Parameters and switches in namelist gwsctl
  ! =============================================
  ! see subroutine setgws for setup of namelist and modifications
  ! read from job control file 
  !
  ! Switches:
  ! ---------
  !
  ! lextro     : true for  hines' doppler spreading extrowave 
  !              parameterization (Hines, 1997a,b)
  !
  ! lfront     : true for gw emerging from fronts and background
  !              (Charron and Manzini, 2002)  
  !
  ! lozpr      : true for background enhancement associated with 
  !              precipitation (Manzini et al., 1997)
  !
  ! iheatcal   : switch for activating upper atmosphere processes
  !              iheatcal = 1 to calculate heating rates and diffusion 
  !                           coefficient.
  !              iheatcal = 0  only momentum flux deposition
  !
  ! Parameters:
  ! -----------
  !
  !
  ! rmscon        : root mean square gravity wave wind at lowest level (m/s)
  !
  ! emiss_lev     : number of levels above the ground at which gw
  !                 are emitted (attention: this is vertical resolution
  !                 dependent. Must be generalized)
  !
  ! kstar         : typical gravity wave horizontal wavenumber (1/m)
  !
  ! m_min         : minimum bound in  vertical wavenumber (1/m)
  !               
  ! rms_front     : rms frontal gw wind at source level  (m/s)
  !
  ! front_thres   : minimum value of the frontogenesis function
  !                 for which gw are emitted from fronts [(K/m)^2/hour]
  !
  !
  !  pcrit        : critical precipitation value (mm/d) above which 
  !                 rms gw wind enhancement is applied
  !
  !  pcons        : adimensional facto rfor  background enhancement 
  !                 associated with precipitation 
  ! 
  ! B) Internal swithches and constants 
  ! ===================================
  ! 
  !  naz        : actual number of horizontal azimuths used
  !               (code set up presently for only naz = 4 or 8)
  !
  !  slope      : slope of incident vertical wavenumber spectrum
  !               (code set up presently for slope 1., 1.5 or 2.)
  !               (if m_min not zero, code set up only for ??)
  !
  !  f1         : "fudge factor" used in calculation of trial value of
  !                azimuthal cutoff wavenumber m_alpha (1.2 <= f1 <= 1.9)
  !
  !  f2         : "fudge factor" used in calculation of maximum
  !                permissible instabiliy-induced cutoff wavenumber 
  !                (0.1 <= f2 <= 1.4)
  !
  !  f3         : "fudge factor" used in calculation of maximum 
  !               permissible molecular viscosity-induced cutoff wavenumber 
  !                m_sub_m_mol 
  !
  !  f5         : "fudge factor" used in calculation of heating rate
  !                (1 <= f5 <= 3)
  !
  !  f6         : "fudge factor" used in calculation of turbulent 
  !                diffusivity coefficient 
  !
  !  ksmin      : additional parameter to define a latitudinally varying kstar
  !                (ksmin = 1.e-5 [1/m] used maecham4)        
  !  ksmax      : additional parameter to define a latitudinally varying kstar
  !                (ksmax = 1.e-4  [1/m] used maecham4) 
  !                kstar = ksmin/( coslat+(ksmin/ksmax) ) as for maecham4
  !
  !  icutoff    : switch  for activating exponenetial damp of gwd, 
  !               heating and diffusion arrays above alt_cutoff
  !               icutoff = 1 : The exponentially damp is on
  !               icutoff = 0 : The arrays are not modified (default)
  !
  !  alt_cutoff : altitude (meters) above which exponential decay applied
  !
  !  smco       : smoother used to smooth cutoff vertical wavenumbers
  !               and total rms winds before calculating drag or heating.
  !                 (==> a 1:smco:1 stencil used; smco >= 1.)
  !
  !  nsmax      : number of times smoother applied ( >= 1)
  !
  !===========================================================================

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PUBLIC

  ! Switches in namelist gwsctl
  ! ---------------------------

  LOGICAL :: lextro   = .TRUE.
  LOGICAL :: lfront   = .FALSE.
  LOGICAL :: lozpr    = .FALSE.
  INTEGER :: iheatcal =   0

  ! Parameters in namelist gwsctl
  ! ----------------------------

  REAL(dp) :: rmscon     = 1.0_dp

  INTEGER  :: emiss_lev  = 7              

  REAL(dp) :: kstar      = 5.0e-5_dp        ! = 2*pi/(126 km)
  REAL(dp) :: m_min      = 0.0_dp

  REAL(dp) :: rms_front   = 2.0_dp
  REAL(dp) :: front_thres = 0.12_dp              ! default for T42

  REAL(dp) :: pcrit       = 5.0_dp
  REAL(dp) :: pcons       = 4.75_dp

  !---------------------------------------------------------------

  ! Define namelist
  !----------------
  NAMELIST/gwsctl/                                           &
       lextro, lfront, lozpr, iheatcal, rmscon, emiss_lev,   &
       kstar, m_min, rms_front, front_thres,  pcrit,  pcons  

  !---------------------------------------------------------------

  ! Internal switches and constants 
  ! --------------------------------

  INTEGER :: naz   = 8

  REAL(dp) :: slope = 1.0_dp

  REAL(dp) :: f1    = 1.5_dp 
  REAL(dp) :: f2    = 0.3_dp 
  REAL(dp) :: f3    = 1.0_dp 
  REAL(dp) :: f5    = 1.0_dp 
  REAL(dp) :: f6    = 0.5_dp   

  REAL(dp) :: ksmin = 1.e-5_dp       
  REAL(dp) :: ksmax = 1.e-4_dp       

  INTEGER  :: icutoff    = 0   
  REAL(dp) :: alt_cutoff = 105.e3_dp

  REAL(dp) :: smco       = 2.0_dp      !  (test value: smco = 1.0)
  INTEGER  :: nsmax      = 5           !  (test value: nsmax = 2)

  !-----------------------------------------------------------------

END MODULE mo_gwspectrum
