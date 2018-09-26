!+ <A one line description of this module> 
! 
MODULE mo_zenith
  ! 
  ! Description: 
  !   <Say what this module is for> 
  ! 
  ! Current Code Owner: <Name of person responsible for this code> 
  ! 
  ! History: 
  !  
  ! Version   Date     Comment 
  ! -------   ----     ------- 
  ! <version> <date>   Original code. <Your name> 
  ! 
  ! Code Description: 
  !   Language:           Fortran 90. 
  !   Software Standards: "European Standards for Writing and  
  !     Documenting Exchangeable Fortran 90 Code". 
  ! 
  ! Modules used: 
  ! 
  USE mo_mpi,         ONLY: p_parallel_io, p_bcast, p_io
  USE mo_kind,        ONLY: dp
  USE mo_doctor,      ONLY: nout
  USE mo_exception,   ONLY: finish
  USE mo_radiation,   ONLY: cecc, cobld, clonp

  IMPLICIT NONE 

  LOGICAL :: module_configured  = .FALSE.
  LOGICAL :: module_initialized = .FALSE.

CONTAINS 

  !=================================================================================================
  SUBROUTINE compute_declination(jsec_add, declination)

    ! Compute solar declination angle

    ! Called from *jsbach_start_timestep* in *mo_jsbach* at the start of each time step.

    USE mo_time_control,    ONLY: l_orbvsop87, julian_date, get_date_components, get_year_len, &
         current_date, start_date, next_date
    USE mo_time_conversion, ONLY: TC_get
    USE mo_time_base,       ONLY: Set_JulianDay
    USE mo_constants,       ONLY: api
    USE mo_orbit,           ONLY: orbit

    INTEGER,  INTENT(in)  :: jsec_add     ! Number of seconds to be added to the Julian date
    REAL(dp), INTENT(out) :: declination  ! Solar declination

    TYPE(julian_date)  :: jul_date, jul_date_pal
    REAL(dp) :: zdoy                    !! Julian date    
    REAL(dp)    :: zvemar0, zyearve, zdoyve
    INTEGER :: yr, mo, dy, hr, mn, se, jday, jsec
    REAL(dp)     :: zdummy1, zdummy2, zdummy3

    CALL get_date_components(next_date, yr, mo, dy, hr, mn, se)
    CALL TC_get             (next_date, jday, jsec)

!!$    CALL get_date_components(current_date, yr, mo, dy, hr, mn, se)
!!$    CALL TC_get             (current_date, jday, jsec)
!!$    jsec =jsec + jsec_add
    CALL Set_JulianDay(yr, mo, dy, jsec, jul_date)

    IF (.NOT. l_orbvsop87) THEN           ! This used to be lpaleo=TRUE
       ! CALL orbit(zdoyve, cdisse, zdec, zra)
       CALL get_date_components(start_date, yr ,mo ,dy ,hr ,mn ,se)
       CALL Set_JulianDay(1900, 1, 1, 0, jul_date_pal)
       zdoy =  (jul_date    %day + jul_date    %fraction) &
            - (jul_date_pal%day + jul_date_pal%fraction)

       ! Determine day of vernal equinox in March 1900
       zvemar0 = 78.41_dp - 0.0078_dp*(1900-1987) + 0.25_dp*MOD(1900,4)
       zyearve = zdoy + get_year_len() - zvemar0
       zdoyve = MOD(zyearve/get_year_len(),1.0_dp)*2._dp*api

       CALL orbit(zdoyve, zdummy1, declination, zdummy2)

    ELSE
       !CALL orbit(zdoy, zra,zdec,zdis,zgha)
       zdoy = jul_date%day + jul_date%fraction

       CALL orbit(zdoy, zdummy1, declination, zdummy2, zdummy3)

    ENDIF

  END SUBROUTINE compute_declination
  !
  !=================================================================================================
  SUBROUTINE init_zenith

    ! Pre-compute some work quantities for domain

    module_initialized = .TRUE.

  END SUBROUTINE init_zenith
  !
  !=================================================================================================
  FUNCTION get_zenith_angle(kdim, domain, time_add, declination) RESULT(cangle)

    ! Compute cosine of zenith angle for current time and all land points of current CALL to LSS.
    ! This returns the value corrected for relative daylength
    
    ! *** This assumes currently that the model simulates the diurnal cycle ***

    USE mo_jsbach_grid , ONLY: domain_type, kindex, nidx
    USE mo_time_control, ONLY: get_clock, current_date

    INTEGER,           INTENT(in) :: kdim        !! Dimension of result vector (must be allocated)
    TYPE(domain_type), INTENT(in) :: domain
    REAL(dp),          INTENT(in) :: time_add    !! 
    REAL(dp),          INTENT(in) :: declination !! solar declination angle

    REAL(dp) :: cangle(kdim)           !! Cosine of solar zenith angle

    REAL(dp) :: ztim1(kdim), ztim2(kdim), ztim3(kdim)
    REAL(dp) :: zclock                 !! Daytime fraction in radians

    ! Note that components nidx, kidx0, kidx1 and kindex are only set at the beginning of
    ! jsbach_interface, while component nland is set in jsbach_init by the decomposition. This means that, from
    ! the jsbach driver, get_zenith_angle must be called with kdim = domain%nland, while kdim = nidx is
    ! only possible if called from jsbach_interface or any of the subroutines called from there.

    IF (kdim /= domain%nland .AND. kdim /= nidx) &
         CALL finish('get_zenith_angle','Wrong dimension!')

    IF (.NOT. module_initialized) CALL init_zenith

    ! Cf. Eq. 10.53 in Jacobsen, Fundamental Atmospheric Modeling
    !      -cos(zclock)*cos(lon) + sin(zclock)*sin(lon) = 
    !    = -cos(zclock+lon) = cos(zclock+lon-pi) = cos(H_a) where H_a local hour angle
    ! => ztim2*cos(lon) + ztim3*sin(lon) = cos(declination)*cos(lat)*cos(H_a)

    zclock = get_clock(current_date)
    zclock = zclock + time_add
    IF (kdim == domain%nland) THEN
       ztim1(:) =   SIN(declination)               * domain%sinlat(:)
       ztim2(:) = - COS(declination) * COS(zclock) * domain%coslat(:)
       ztim3(:) =   COS(declination) * SIN(zclock) * domain%coslat(:)
       cangle(:) = ztim1(:) + ztim2(:) * domain%coslon(:) + ztim3(:) * domain%sinlon(:)
    ELSE
       ztim1(:) = SIN(declination)               * domain%sinlat(kindex)
       ztim2(:) = - COS(declination) * COS(zclock) * domain%coslat(kindex)
       ztim3(:) = COS(declination) * SIN(zclock) * domain%coslat(kindex)
       cangle(:) = ztim1(:) + ztim2(:) * domain%coslon(kindex) + ztim3(:) * domain%sinlon(kindex)
    ENDIF


    ! For diurnal cycle = "on" relative daylength is either 0 or 1 depending on 
    ! whether cangle is negative or positive.
    ! Correction for relative daylength:
    cangle = MAX(0._dp, cangle)

  END FUNCTION get_zenith_angle

END MODULE mo_zenith
