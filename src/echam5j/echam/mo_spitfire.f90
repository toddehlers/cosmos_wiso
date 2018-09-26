 MODULE mo_spitfire

  !
  !     SPITFIRE advection
  !
  !     (SPlit Implementation of Transport
  !     using Flux Integral REpresentations)
  !     by Phil Rasch, NCAR, 1998
  !  
  ! P. Rasch, NCAR, 1998, Original source
  ! U. Schlese, DKRZ, February 1998, modified for Echam4
  ! U. Schlese, DKRZ, November 1999, modified for Echam5
  ! U. Schlese, DKRZ, November 1999, modified for greater precision
  ! L. Kornblueh, MPI, July 2001, module packing and encapsulating,
  !                               precision standardization
  ! T. Diehl, DKRZ, September 2000, modified for parallel version
  ! T. Diehl, DKRZ, June 2001, modified for round-off problems
  !

  USE mo_kind,          ONLY: wp
  USE mo_doctor,        ONLY: nerr
  USE mo_exception,     ONLY: finish
  USE mo_decomposition, ONLY: gdc => global_decomposition, &
                              ldc => local_decomposition

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_spitfire
  PUBLIC :: setup_spitfire
  PUBLIC :: spitfire_transport
  PUBLIC :: spitfire_tendencies
!---wiso-code
!  PUBLIC :: pole_filter
  PUBLIC :: pole_filter_wiso
!---wiso-code-end

  ! function declarations

#ifdef __uxp__
  REAL(wp), EXTERNAL :: minmod, medan
#endif

  ! basic grid point resolution parameters
  !
  INTEGER ::   nplon    ! global number of longitudes per latitude 
  INTEGER ::   nplond   ! nplon + extra point(s) for interpolation or
                        ! extra point for wall centered grid x-direction
  INTEGER ::   nplat    ! global number of gaussian latitudes
  INTEGER ::   nplatext ! 2*nplat+1: number of points on meridian 
                        ! for rotated grid + extra point(s) for interpolation
  !
  INTEGER ::   plon    ! number of longitudes
  INTEGER ::   plev    ! number of vertical levels
  INTEGER ::   plat    ! number of latitudes
  INTEGER ::   pcnst   ! number of constituents (including water vapor)
  INTEGER ::   plevp   ! plev + 1
  !
  REAL(wp) :: dymin, dzmin, dx, arad
  !
  INTEGER :: njeven, nkeven
  !
  REAL(wp), ALLOCATABLE :: ub(:,:,:)  ! u-component of wind
  REAL(wp), ALLOCATABLE :: vb(:,:,:)  ! v-component of wind
  REAL(wp), ALLOCATABLE :: fb(:,:,:,:)  ! constituent fields !! indexing!!
  !
  REAL(wp), ALLOCATABLE :: ub_g(:,:)  ! global level of u-component
  REAL(wp), ALLOCATABLE :: vb_g(:,:)  ! global level of v-comp.
  REAL(wp), ALLOCATABLE :: fb_g(:,:)  ! global level of tracer
  !
  REAL(wp), ALLOCATABLE :: uwdp(:,:), vwdp(:,:), wwdp(:,:,:)
  REAL(wp), ALLOCATABLE :: uw(:,:), vw(:,:)
  REAL(wp), ALLOCATABLE :: zw(:,:) 
  REAL(wp), ALLOCATABLE :: xc(:), xw(:), yc(:), yw(:), ylat(:), ywext(:)
  REAL(wp), ALLOCATABLE :: dmu(:) 
  REAL(wp), ALLOCATABLE :: ylat_l(:)
  !
  INTEGER, ALLOCATABLE :: idep(:,:), jdep(:,:), kdep(:,:,:)
  INTEGER, ALLOCATABLE :: idom(:,:), jdom(:,:)
  INTEGER, ALLOCATABLE :: bini(:,:), bini_l(:,:)
  INTEGER, ALLOCATABLE :: nptslat(:), nptslat_l(:)
  INTEGER, ALLOCATABLE :: jeven(:), keven(:)

  ! Integer for dimensioning certain work arrays
  INTEGER :: nwork

  ! Buffers for pole filter communication
!!!  REAL,    ALLOCATABLE :: bufs_x2p (:,:)
!!!  REAL,    ALLOCATABLE :: bufr_x2p (:,:)

CONTAINS

!!============================================================================
!!============================================================================
!!
!! Section 1: init_spitfire is public interface 
!!
  SUBROUTINE init_spitfire

    USE mo_advection,     ONLY: nalatd, nalond, nalev, nacnst, nalat, nalon 
    USE mo_tracer,        ONLY: jps, trlist
!---wiso-code
    USE mo_wiso,          ONLY: lwiso, nwiso
!---wiso-code-end

    plon    = ldc%nglon
    plev    = ldc%nlev
    plat    = ldc%nglat

    pcnst   = jps + trlist% nadvec + 1

!---wiso-code

!   add number of wiso-fields
    IF (lwiso) THEN
      pcnst   = pcnst + nwiso*3
    END IF

!---wiso-code-end

    plevp   = ldc%nlev+1

    nplon    = ldc%nlon
    nplond   = ldc%nlon + 1
    nplat    = ldc%nlat
    nplatext = 2*ldc%nlat + 1
    nwork    = MAX(nplond,nplatext)

    nalev  = plev
    nacnst = pcnst
    nalat  = plat
    nalon  = plon
    nalatd = ldc%nglat/2
    nalond = ldc%nglon+1

    IF (.NOT. ALLOCATED(uwdp)) THEN
      ALLOCATE (                                      &
        ub        (plon,     plat,    plev        ),  &
        vb        (plon,     plat,    plev        ),  &
        fb        (plon,     plat,    plev, pcnst ),  &
        ub_g      (nplon,    nplat                ),  &
        vb_g      (nplon,    nplat                ),  &
        fb_g      (nplon,    nplat                ),  &
        uwdp      (nplond,   nplat                ),  &
        vwdp      (nplatext, nplon/2              ),  &
        wwdp      (plon,     plev,    plat        ),  &
        uw        (nplond,   nplat                ),  &
        vw        (nplond,   nplat+1              ),  &
        zw        (plon,     plevp                ),  &
        xc        (nplond                         ),  &
        xw        (nplond                         ),  &
        yc        (nplat                          ),  &
        yw        (nplat+1                        ),  &
        ylat      (nplat                          ),  &
        ylat_l    (plat                           ),  &
        dmu       (nplat                          ),  &
        ywext     (nplatext                       ),  &
        idep      (nplond,   nplat                ),  &
        idom      (nplond,   nplat                ),  &
        jdep      (nplatext, nplon/2              ),  &
        jdom      (nplatext, nplon/2              ),  &
        kdep      (plon,     plev,    plat        ),  &
        bini      (nplon,    nplat                ),  &
        bini_l    (nplon,    plat                 ),  &
        nptslat   (nplat                          ),  &
        nptslat_l (plat                           ))

    END IF

    ub(:,:,:)   = 0.0_wp
    vb(:,:,:)   = 0.0_wp
    fb(:,:,:,:) = 0.0_wp

    CALL initrgrd1

    IF (.NOT. ALLOCATED (jeven)) THEN
      ALLOCATE (jeven(njeven), keven(nkeven))
    END IF

    CALL initrgrd2

    CALL inifilter
  
  END SUBROUTINE init_spitfire

!!============================================================================

  SUBROUTINE initrgrd1

    USE mo_constants,     ONLY: api, a
    USE mo_gaussgrid,     ONLY: gauaw
    USE mo_hyb,           ONLY: cetah

    ! Local arrays
    REAL(wp) ::  zgauw(nplat)  ! Gaussian weights for both hemispheres
                               ! summing up to 2
    REAL(wp) :: zgmu(nplat)    ! Gaussian latitudes (dummy)

    ! Local variables
    INTEGER :: i, j, k

    REAL(wp) :: pi         
    
    pi = api
    arad = a

    ! Call gauaw to get Gaussian weights tailored for spitfire,
    ! i.e. for all latitudes and summing up to 2.

    CALL gauaw(zgmu,zgauw,nplat)

    !     grid for vertical transforms

!    WRITE (nerr,*) ' number of interfaces ', plevp
    DO k = 1,plevp
!      WRITE (nerr,*) ' initrgrd: cetah ', k, cetah(k)
      DO i = 1,plon
        zw(i,k) = cetah(k)
      END DO
    END DO

    dzmin = 1.0e36_wp

    DO k = 1,plev
      !        write (nerr,*) ' dzmin dz ', dzmin, abs(cetah(k+1)-cetah(k))
      dzmin = MIN(dzmin,ABS(cetah(k+1)-cetah(k)))
    END DO

!    WRITE (nerr,*) ' smallest z grid increment is ', dzmin

    nkeven = ABS(cetah(plevp)-cetah(1))/dzmin + 2

!    WRITE (nerr,*) 'regular z grid requires nypts = ', nkeven, dzmin

    !     grid for horizontal transforms

    yw(1) = -1.0_wp
    DO j = 1,nplat
      yw(j+1) = yw(j) + zgauw(j)
      dmu(j) = zgauw(j)
    END DO
    !     write (nerr,*) ' np sin lat ', yw(plat+1)
    yw(1)       = -0.5_wp*pi
    yw(nplat+1) =  0.5_wp*pi

    ywext(1)       = -1.0_wp
    ywext(nplat+1) =  1.0_wp

    DO j = 2,nplat
      ywext(j) = yw(j)
      yw(j) = ASIN(yw(j))
    END DO
    DO j = 1,nplat
      ylat(j) = (yw(j+1)+yw(j))*0.5_wp
      yc(j) = ylat(j)
    END DO

    ! local version of ylat for pole_filter
    DO j = 1,plat
       ylat_l(j) = (yw(ldc%glat(j)+1) + yw(ldc%glat(j)))*0.5_wp
    END DO

    dx = 2*pi/nplon
    DO i = 1,nplon+1
      xw(i) = (i-1.0_wp)*dx
      xw(i) = (i-1.5_wp)*dx
    END DO

    !     cell centered quantities
    DO i = 1,nplon
      xc(i) = (xw(i) + xw(i+1))*0.5_wp
    END DO
    xc(nplon+1) = xc(nplon) + dx

    DO j = nplat+2,nplatext
      ywext(j) = ywext(j-1) + (ywext(j-nplat)-ywext(j-nplat-1))
    END DO


    dymin = 1.0e36_wp
    DO j = 1,nplatext-1
      dymin = MIN(dymin,ABS(ywext(j+1)-ywext(j)))
    END DO

    njeven = ABS(ywext(nplatext)-ywext(1))/dymin + 2
!    WRITE (nerr,*) 'even y grid requires nypts = ', njeven, dymin

    IF (njeven < nplatext) THEN
!      WRITE (nerr,*) ' initrgrd1: something is terribly wrong ',           &
!           njeven, nplatext
      call finish ('initrgrd1','something is terribly wrong')
    ENDIF

!    WRITE (nerr,*) ' initrgrd1 complete '

  END SUBROUTINE initrgrd1
!!============================================================================
  SUBROUTINE initrgrd2

    USE mo_hyb,       ONLY: cetah
    USE mo_exception, ONLY: finish

    ! Local variables
    INTEGER :: i, j, k

    REAL(wp) :: zeven, yeven

!    WRITE (nerr,*) ' entering initrgrd2 '

    dzmin = (cetah(plevp)-cetah(1))/(nkeven-1)
!    WRITE (nerr,*) ' dzmin reset to ', dzmin
    DO k = 1,nkeven
      zeven = (k-1)*dzmin + cetah(1)
      IF (k == nkeven) THEN
      ! due to round-off errors, zeven can be larger than cetah(plevp) for
      ! k=nkeven ; since zeven should be exactly equal to cetah(plevp) for
      ! j=nkeven by construction, we set it to this value.
         zeven = cetah(plevp)
      END IF
      keven(k) = 0
      DO i = 1,plev
        IF (zeven >= cetah(i) .AND. zeven <= cetah(i+1)) THEN
          keven(k) = i
        ENDIF
      END DO

      IF (keven(k) == 0) THEN
        WRITE (nerr,*) 'couldnt map zeven ', k, zeven
        CALL finish('initgrd2', 'couldnt map zeven')
      ELSE
        i = keven(k)
        !           write (nerr,*) 'even grid ', k, ' at z ', zeven, &
        !                       ' is on uneven grid interval ',   &
        !                       i, cetah(i), cetah(i+1)
      ENDIF
    END DO

    dymin = (ywext(nplatext)-ywext(1))/(njeven-1)
!    WRITE (nerr,*) ' dymin reset to ', dymin
    DO j = 1,njeven
      yeven = (j-1)*dymin + ywext(1)
      IF (j == njeven) then 
        ! due to round-off errors, yeven can be larger than ywext(nplatext) for
        ! j=njeven ; since yeven should be exactly equal to ywext(nplatext) for 
        ! j=njeven by construction, we set it to this value.
        yeven = ywext(nplatext)
      ENDIF
      jeven(j) = 0
      DO i = 1,nplatext-1
        IF (yeven >= ywext(i) .AND. yeven <= ywext(i+1)) THEN
          jeven(j) = i
        ENDIF
      END DO
      IF (jeven(j) == 0) THEN
        WRITE (nerr,*) 'couldnt map yeven ', j, yeven
        CALL finish('initgrd2', 'couldnt map yeven')
      ELSE
        i = jeven(j)
        !           write (nerr,*) 'even grid ', j, ' at y ', yeven, &
        !                       ' is on uneven grid interval ',   &
        !                       jeven(j), ywext(i), ywext(i+1)
      ENDIF
    END DO

!    WRITE (nerr,*) ' initrgrd2 complete '

  END SUBROUTINE initrgrd2
!!============================================================================
  SUBROUTINE inifilter
    !
    !     Phil Rasch, NCAR, 1998
    !

    INTEGER :: i, j, n

    REAL(wp) :: reduce
    INTEGER :: npole, nfilt, nmin
    INTEGER, PARAMETER :: pow2(8) = (/1, 2, 4, 8, 16, 32, 64, 128 /)

    reduce = 0.6_wp  ! factor by which to reduce the stencil
                     ! (0 (max reduction) to 1. (standard))

    ! compute nptslat_l and bini_l for routine filter
    DO j = 1,nplat
      nptslat(j) = MAX(1,MIN(nplon,INT(reduce/(COS(ylat(j))+1.e-10_wp))))
      IF (nptslat(j) > 1) THEN
        npole = ABS(nptslat(j)-pow2(1))
        nmin = 1
        DO n = 1,8
          IF (ABS(nptslat(j)-pow2(n)) <= npole) THEN
            npole = ABS(nptslat(j)-pow2(n))
            nmin = n
          ENDIF
        END DO
        nptslat(j) = pow2(nmin)
      ENDIF
    END DO

    DO j = 1,nplat
      nfilt = nptslat(j)
      DO i = 1,nplon
        bini(i,j) = (i+nfilt-1)/nfilt
      END DO
!      WRITE (nerr,*) j, ' # in new filter ', nptslat(j)
    END DO

!   compute nptslat_l and bini_l for routine pole_filter
    DO j = 1,plat
      nptslat_l(j) = MAX(1,MIN(nplon,INT(reduce/(COS(ylat_l(j))+1.e-10_wp))))
      IF (nptslat_l(j) > 1) THEN
        npole = ABS(nptslat_l(j)-pow2(1))
        nmin = 1
        DO n = 1,8
          IF (ABS(nptslat_l(j)-pow2(n)) <= npole) THEN
            npole = ABS(nptslat_l(j)-pow2(n))
            nmin = n
          ENDIF
        END DO
        nptslat_l(j) = pow2(nmin)
      ENDIF
    END DO

    DO j = 1,plat
      nfilt = nptslat_l(j)
      DO i = 1,nplon
        bini_l(i,j) = (i+nfilt-1)/nfilt
      END DO
!      WRITE (nerr,*) j, ' # in new filter ', nptslat_l(j)
    END DO

  END SUBROUTINE inifilter

!!============================================================================
!!============================================================================

!OCL NOALIAS

  SUBROUTINE setup_spitfire

    !  preliminary for parallel version NOT running parallel !!!!
    !  i.e. plon = all longitudes of one latitude
    
    ! Fills 3dim arrays for the SPITFIRE scheme
    
    USE mo_tracer,        ONLY: jps, trlist
    USE mo_gaussgrid,     ONLY: gl_sqcst
    USE mo_hyb,           ONLY: cetah
    USE mo_scan_buffer,   ONLY: u, v
    USE mo_memory_g1a,    ONLY: qm1, xlm1, xim1, xtm1, alpsm1
!---wiso-code
    USE mo_memory_wiso,   ONLY: wisoqm1, wisoxlm1, wisoxim1
    USE mo_wiso,          ONLY: lwiso, nwiso
!---wiso-code-end
    
    REAL(wp) :: zlimit, zrcst
    INTEGER :: i, jk, jl, jt, l, jglat, ita
    
    REAL(wp) :: zphm1(plon,plevp,plat), zpsm1(plon,plat)
    REAL(wp) :: zpde (plon,plev, plat)        ! density in eta-coordin.
    
    !  External subroutines 
    EXTERNAL pres
    
    !  Intrinsic functions 
    INTRINSIC ABS
    
    !  Executable statements 
    
    ! Set scalar fields to zero if the absolute value
    ! is smaller than a hardware dependent value
    ! to prevent slt scheme from underflow.
    ! do this even if spitfire is skipped to ensure same results in SCM
    
    zlimit = 1.E-200_wp

    DO i = 1, plat

      l = plat-i+1  ! S->N latitude index

      WHERE (ABS(qm1(:,:,l)) < zlimit) 
        qm1(:,:,l) = 0.0_wp
      END WHERE
    
      WHERE (ABS(xlm1(:,:,l)) < zlimit)
        xlm1(:,:,l) = 0.0_wp
      END WHERE
    
      WHERE (ABS(xim1(:,:,l)) < zlimit)
        xim1(:,:,l) = 0.0_wp
      END WHERE
    
      WHERE (ABS(xtm1(:,:,:,l)) < zlimit) 
        xtm1(:,:,:,l) = 0.0_wp
      END WHERE

      zpsm1(:,l) = EXP(alpsm1(1:plon,l))
    
      ! Compute half level pressure at t-dt for zpdel.
    
      CALL pres(zphm1(:,:,l), plon, zpsm1(:,l), plon)

      ! Fill up 3d arrays for Spitfire from south to north
      ! and divide winds by cos(latitude)
      
      jglat  = ldc%glat(l)   ! global index (continuous from north to south)
      zrcst  = 1.0_wp/gl_sqcst(jglat)                 ! 1./cos(latitude)
      
      ! arrays for SPITFIRE must have the format south -> north
      DO jk = 1, plev
        DO jl = 1, plon
          zpde(jl,jk,i) = (zphm1(jl,jk+1,l)-zphm1(jl,jk,l))  &
               /(cetah(jk+1)-cetah(jk))
          ub(jl,i,jk)   = u(jl,jk,l)*zrcst
          vb(jl,i,jk)   = v(jl,jk,l)*zrcst
          fb(jl,i,jk,1) = qm1(jl,jk,l)*zpde(jl,jk,i)
          fb(jl,i,jk,2) = xlm1(jl,jk,l)*zpde(jl,jk,i)
          fb(jl,i,jk,3) = xim1(jl,jk,l)*zpde(jl,jk,i)
        END DO
      END DO
      
      ita = 0
!---wiso-code

      IF (lwiso) THEN
        DO jt = 1, nwiso

          ita = ita + 1
          DO jk = 1, plev
            DO jl = 1, plon
              fb(jl,i,jk,jps+ita) = wisoqm1(jl,jk,jt,l)*zpde(jl,jk,i)
            END DO
          END DO
          ita = ita + 1
          DO jk = 1, plev
            DO jl = 1, plon
              fb(jl,i,jk,jps+ita) = wisoxlm1(jl,jk,jt,l)*zpde(jl,jk,i)
            END DO
          END DO
          ita = ita + 1
          DO jk = 1, plev
            DO jl = 1, plon
              fb(jl,i,jk,jps+ita) = wisoxim1(jl,jk,jt,l)*zpde(jl,jk,i)
            END DO
          END DO
          
        END DO
      END IF 

!---wiso-code-end

      DO jt = 1, trlist% ntrac
        IF (trlist% ti(jt)% ntran == 0 ) CYCLE
        ita = ita + 1
        DO jk = 1, plev
          DO jl = 1, plon
            fb(jl,i,jk,jps+ita) = xtm1(jl,jk,jt,l)*zpde(jl,jk,i)
          END DO
        END DO
      END DO
      
      ! put density of air on last tracer
      DO jk = 1,plev
        DO jl = 1, plon
          fb(jl,i,jk,pcnst) = zpde(jl,jk,i)
        END DO
      END DO
    END DO

  END SUBROUTINE setup_spitfire

!!============================================================================
!!============================================================================

!OCL NOALIAS

  SUBROUTINE spitfire_tendencies(psm1cor, pscor)
    !
    !**** *SPITTEN* - TENDENCIES OF *SPITFIRE* ADVECTION, FIXING OF AIR MASS
    !
    !        U. SCHLESE    DKRZ - HAMBURG    FEB-98
    !        U. SCHLESE    DKRZ - OCT 99  preliminary version for ECHAM5
    !
    !  I.Kirchner, Dec/2000, date/time control
    !
    USE mo_time_control,   ONLY: time_step_len
    USE mo_tracer,         ONLY: jps, trlist
    USE mo_memory_g1a,     ONLY: alpsm1, qm1, xlm1, xim1, xtm1
    USE mo_scan_buffer,    ONLY: alps, qte, xlte, xite, xtte
!---wiso-code
    USE mo_memory_wiso,    ONLY: wisoqm1, wisoxlm1, wisoxim1,   &
                                 wisoqte, wisoxlte, wisoxite
    USE mo_wiso,           ONLY: lwiso, nwiso
!---wiso-code-end

    REAL(wp), INTENT(in) :: pscor
    REAL(wp), INTENT(in) :: psm1cor
    
    INTEGER :: ilat
    INTEGER :: jk
    INTEGER :: jl
    INTEGER :: jlat
    INTEGER :: jt
    INTEGER :: ita
    REAL(wp) :: zlpsc
    REAL(wp) :: zlpsm1c
    REAL(wp) :: ztmst
    
    ztmst = time_step_len

    DO jlat  = 1, plat                  ! local index north -> south
      ilat  = plat-jlat+1               ! local latitude from south to north
      !
      !
      ! ------------------------------------------------------------------
      !
      !*         COMPUTE TENDENCIES 
      !          ------- ---------- 
      !
      !
      DO  jk=1,plev
        DO  jl=1,plon
          qte    (jl,jk,jlat) = qte(jl,jk,jlat)                    &
                                  + (fb(jl,ilat,jk,1)- qm1(jl,jk,jlat))/ztmst
          xlte   (jl,jk,jlat) = xlte(jl,jk,jlat)                    &
                                  + (fb(jl,ilat,jk,2)-xlm1(jl,jk,jlat))/ztmst
          xite   (jl,jk,jlat) = xite(jl,jk,jlat)                    &
                                  + (fb(jl,ilat,jk,3)-xim1(jl,jk,jlat))/ztmst
        END DO
      END DO
      !
      ita = 0

!---wiso-code

      IF (lwiso) THEN
        DO jt = 1, nwiso

          ita = ita + 1
          DO jk=1,plev
            DO jl=1,plon
              wisoqte(jl,jk,jt,jlat)=wisoqte(jl,jk,jt,jlat)+                 &
                   (fb(jl,ilat,jk,ita+jps)-wisoqm1(jl,jk,jt,jlat))/ztmst
            END DO
          END DO
          ita = ita + 1
          DO jk=1,plev
            DO jl=1,plon
              wisoxlte(jl,jk,jt,jlat)=wisoxlte(jl,jk,jt,jlat)+                 &
                   (fb(jl,ilat,jk,ita+jps)-wisoxlm1(jl,jk,jt,jlat))/ztmst
            END DO
          END DO
          ita = ita + 1
          DO jk=1,plev
            DO jl=1,plon
              wisoxite(jl,jk,jt,jlat)=wisoxite(jl,jk,jt,jlat)+                 &
                   (fb(jl,ilat,jk,ita+jps)-wisoxim1(jl,jk,jt,jlat))/ztmst
            END DO
          END DO

        END DO
      END IF 

!---wiso-code-end

      DO jt=1,trlist% ntrac
        IF (trlist% ti(jt)% ntran == 0) CYCLE
        ita = ita + 1
        DO jk=1,plev
          DO jl=1,plon
            xtte(jl,jk,jt,jlat)=xtte(jl,jk,jt,jlat)+                 &
                 (fb(jl,ilat,jk,ita+jps)-xtm1(jl,jk,jt,jlat))/ztmst
          END DO
        END DO
      END DO
      !
      !    ------------------------------------------------------------------
      !
      !*          FIX MASS OF AIR
      !           --- ---- -- ---
      !
      zlpsm1c=LOG(psm1cor)
      zlpsc=LOG(pscor)
      !
      DO jl=1,plon
        alpsm1(jl,jlat)=alpsm1(jl,jlat)+zlpsm1c
        alps(jl,jlat)=alps(jl,jlat)+zlpsc
      END DO
      !
      !     -----------------------------------------------------------------
      !
    END DO

  END SUBROUTINE spitfire_tendencies

!!============================================================================
!!============================================================================

  SUBROUTINE spitfire_transport (etadot, ztodt, istep)

    USE mo_mpi,           ONLY: p_isend, p_recv, p_barrier, p_wait, p_nprocs, &
                                p_pe, p_communicator_d
    USE mo_transpose,     ONLY: gather_gp_level, scatter_gp_level, indx, reorder

    REAL(wp), INTENT(in) :: etadot(ldc%nproma,plevp,ldc%ngpblks)
    REAL(wp), INTENT(in) :: ztodt

    INTEGER, INTENT(in)  :: istep
    
    REAL(wp) :: zetadot(plon,plevp,plat)
    INTEGER :: i, ilat, k, m, nproc, pe, jlat
    INTEGER :: kidx(ldc%spe-1:ldc%d_nprocs-1), midx(ldc%spe-1:ldc%d_nprocs-1)

    CALL reorder(zetadot, etadot)

    ! do the vertical advection 

    DO ilat = 1,plat
      jlat  = plat-ilat+1               ! local latitude from south to north
      CALL setv( zetadot(:,:,jlat), ztodt ,ilat) 
      CALL vadvn (fb,ilat)
    ENDDO

    ! do the horizontal transport

    IF (p_nprocs == 1) THEN
       DO k = 1,plev
          CALL seth (ub(:,:,k), vb(:,:,k), ztodt)
          DO m = 1,pcnst
             CALL hadvn (fb(:,:,k,m), ztodt, istep)
          END DO
       END DO
    ELSE
       nproc = ldc%spe - 1
       DO k = 1,plev
          DO m = 1,pcnst
             pe = gdc(nproc+1)%pe
             CALL p_isend( ub(1,1,k),   pe, 1111, plon*plat)
             CALL p_isend( vb(1,1,k),   pe, 2222, plon*plat)
             CALL p_isend( fb(1,1,k,m), pe, 3333, plon*plat)
             kidx(nproc) = k
             midx(nproc) = m
             nproc = nproc+1
             IF (nproc == ldc%d_nprocs .OR. (k == plev .AND. m == pcnst)) THEN
                IF (indx(p_pe,gdc) < nproc+1) THEN
                   CALL gather_gp_level(ub_g,gdc,1111)
                   CALL gather_gp_level(vb_g,gdc,2222)
                   CALL gather_gp_level(fb_g,gdc,3333)

                   CALL seth ( ub_g, vb_g,  ztodt)
                   CALL hadvn( fb_g, ztodt, istep)
                   
                   CALL scatter_gp_level(fb_g,gdc,4444)
                   
                END IF

                CALL p_barrier(p_communicator_d)

                DO i = ldc%spe-1, nproc-1
                   pe = gdc(i+1)%pe
                   CALL p_recv(fb(:,:,kidx(i),midx(i)),pe,4444)
                END DO

                CALL p_barrier(p_communicator_d)
                CALL p_wait

                nproc = ldc%spe - 1

             END IF
          END DO
       END DO
    END IF

    ! make mixing ratios from updated tracer and updated densities

    DO ilat = 1,plat
      DO m = 1, pcnst-1  ! <--  new density is stored on last constituent !
        DO k = 1,plev
          DO i = 1,plon
            fb(i,ilat,k,m) = fb(i,ilat,k,m)/fb(i,ilat,k,pcnst)           
          END DO
        ENDDO
      END DO
    END DO

  END SUBROUTINE spitfire_transport

!!============================================================================
!!============================================================================
  SUBROUTINE cfdotmc (x, f, fdot)

    !     Phil Rasch, NCAR, 1998

    !     calculate the derivative for the interpolating polynomial
    !     multi column version

    REAL(wp), INTENT(in) :: x(plon, plevp), f(plon, plevp)
    REAL(wp), INTENT(out):: fdot(plon, plevp)      ! derivative at nodes

    !     assumed variable distribution
    !  x1.......x2.......x3.......x4.......x5.......x6     1,plevp points
    !  f1.......f2.......f3.......f4.......f5.......f6     1,plevp points
    !  ...sh1.......sh2......sh3......sh4......sh5....     1,plev points
    !  .........d2.......d3.......d4.......d5.........     2,plev points
    !  .........s2.......s3.......s4.......s5.........     2,plev points
    !  .............dh2......dh3......dh4.............     2,plev-1 points
    !  .............eh2......eh3......eh4.............     2,plev-1 points
    !  ..................e3.......e4..................     3,plev-1 points
    !  .................ppl3......ppl4................     3,plev-1 points
    !  .................ppr3......ppr4................     3,plev-1 points
    !  .................t3........t4..................     3,plev-1 points
    !  ................fdot3.....fdot4................     3,plev-1 points

    INTEGER :: i, k

    REAL(wp) :: s(plon,plev)     ! first divided differences at nodes
    REAL(wp) :: sh(plon,plev)    ! first divided differences between nodes
    REAL(wp) :: d(plon,plev)     ! second divided differences at nodes
    REAL(wp) :: dh(plon,plev)    ! second divided differences between nodes
    REAL(wp) :: e(plon,plev)     ! third divided differences at nodes
    REAL(wp) :: eh(plon,plev)    ! third divided differences between nodes
    REAL(wp) :: pp               ! p prime
    REAL(wp) :: ppl(plon,plev)   ! p prime on left
    REAL(wp) :: ppr(plon,plev)   ! p prime on right

    REAL(wp) :: qpl, qpr, ttt, t, tmin, tmax
    REAL(wp) :: delxh(plon,plev)

    DO k = 1,plev
      !        first divided differences between nodes
      DO i = 1, plon
        delxh(i,k) = (x(i,k+1)-x(i,k))
        sh(i,k) = (f(i,k+1)-f(i,k))/delxh(i,k)
      END DO

      ! first and second divided differences at nodes

      IF (k >= 2) THEN
        DO i = 1,plon
          d(i,k) = (sh(i,k)-sh(i,k-1))/(x(i,k+1)-x(i,k-1))
          s(i,k) = minmod(sh(i,k),sh(i,k-1))
        END DO
      ENDIF
    END DO

    IF (plev == 1) THEN
      DO i = 1,plon
        fdot(i,1) = 0.0_wp
        fdot(i,2) = 0.0_wp
      END DO
      RETURN
    ENDIF

    !     second and third divided diffs between nodes

    DO k = 2,plev-1
      DO i = 1, plon
        eh(i,k) = (d(i,k+1)-d(i,k))/(x(i,k+2)-x(i,k-1))
        dh(i,k) = minmod(d(i,k),d(i,k+1))
      END DO
    END DO

    !     treat the boundaries

    DO i = 1,plon
      e(i,2) = eh(i,2)
      e(i,plev) = eh(i,plev-1)
      !        outside level
      fdot(i,1) = sh(i,1) - d(i,2)*delxh(i,1)                        &
           - eh(i,2)*delxh(i,1)*(x(i,1)-x(i,3))
      fdot(i,1) = minmod(fdot(i,1),3*sh(i,1))
      fdot(i,plevp) = sh(i,plev) + d(i,plev)*delxh(i,plev)           &
           + eh(i,plev-1)*delxh(i,plev)*(x(i,plevp)-x(i,plev-1))
      fdot(i,plevp) = minmod(fdot(i,plevp),3*sh(i,plev))
      !        one in from boundary
      fdot(i,2) = sh(i,1) + d(i,2)*delxh(i,1)                        &
           - eh(i,2)*delxh(i,1)*delxh(i,2)
      fdot(i,2) = minmod(fdot(i,2),3*s(i,2))
      fdot(i,plev) = sh(i,plev) - d(i,plev)*delxh(i,plev)            &
           - eh(i,plev-1)*delxh(i,plev)*delxh(i,plev-1)
      fdot(i,plev) = minmod(fdot(i,plev),3*s(i,plev))
    END DO

    DO k = 3,plev-1
      DO i = 1,plon
        e(i,k) = minmod(eh(i,k),eh(i,k-1))
      END DO
    END DO

    DO k = 3,plev-1
      DO i = 1,plon
        !           p prime at k-0.5
        ppl(i,k)=sh(i,k-1) + dh(i,k-1)*delxh(i,k-1)
        !           p prime at k+0.5
        ppr(i,k)=sh(i,k)   - dh(i,k)  *delxh(i,k)
        t = minmod(ppl(i,k),ppr(i,k))
        !           derivate from parabola thru f(i,k-1), f(i,k), and f(i,k+1)
        pp = sh(i,k-1) + d(i,k)*delxh(i,k-1)
        !           quartic estimate of fdot
        fdot(i,k) = pp                                  &
             - delxh(i,k-1)*delxh(i,k)                  &
             *(  eh(i,k-1)*(x(i,k+2)-x(i,k  ))          &
             + eh(i,k  )*(x(i,k  )-x(i,k-2))            &
             )/(x(i,k+2)-x(i,k-2))
        !           now limit it
        qpl = sh(i,k-1)                                 &
             + delxh(i,k-1)*minmod(d(i,k-1)+e(i,k-1)    &
             *(x(i,k)-x(i,k-2)),                        &
             d(i,k)  -e(i,k)*delxh(i,k))
        qpr = sh(i,k)                                   &
             + delxh(i,k  )*minmod(d(i,k)  +e(i,k)      &
             *delxh(i,k-1),                             &
             d(i,k+1)+e(i,k+1)*(x(i,k)-x(i,k+2)))

        fdot(i,k) = medan(fdot(i,k), qpl, qpr)

        ttt = minmod(qpl, qpr)
        tmin = MIN(0.0_wp,3*s(i,k),1.5_wp*t,ttt)
        tmax = MAX(0.0_wp,3*s(i,k),1.5_wp*t,ttt)

        fdot(i,k) = fdot(i,k) + minmod(tmin-fdot(i,k), tmax-fdot(i,k))

      END DO

    END DO

  END SUBROUTINE cfdotmc

!!============================================================================

  SUBROUTINE cfintv (plevp, plon, x, f, fdot, xin, psistar, kdep2)

    !     Phil Rasch, NCAR, 1998

    !     do plon simultaneous interpolations over vectors of length plevp
    !     in arrays dimensioned plon,plev

    !     input
    INTEGER, INTENT(in) :: plevp, plon
    REAL(wp), INTENT(in) :: x(plon, plevp), f(plon, plevp)
    REAL(wp), INTENT(in) :: fdot(plon, plevp), xin(plon)

    INTEGER, INTENT(in) :: kdep2(plon)

    !     real fxdot(plond)
    !     real fxdd(plond)
    REAL(wp), INTENT(out) :: psistar(plon)

    INTEGER :: i, k

    REAL(wp) :: dx, s, c2, c3, xx
    REAL(wp) :: psi1, psi2, psi3, psim
    REAL(wp) :: cfnew, cfint
    REAL(wp) :: xins

    DO i = 1,plon
      xins = MAX(x(i,1), MIN(xin(i), x(i,plevp)))
      k = kdep2(i)
      dx = (x(i,k+1)-x(i,k))
      s = (f(i,k+1)-f(i,k))/dx
      c2 = (3*s-2*fdot(i,k)-fdot(i,k+1))/dx
      c3 = (fdot(i,k)+fdot(i,k+1)-2*s)/dx**2
      xx = (xins-x(i,k))
      !        fxdot(i) =  (3*c3*xx + 2*c2)*xx + fdot(i,k)
      !        fxdd(i) = 6*c3*xx + 2*c2
      cfint = ((c3*xx + c2)*xx + fdot(i,k))*xx + f(i,k)
      psistar(i) = cfint
      !
      !        leonards limiter
      !        limit the interpolant
      psi1 = f(i,k)+(f(i,k+1)-f(i,k))*xx/dx
      IF (k == 1) THEN
        psi2 = f(i,1)
      ELSE
        psi2 = f(i,k) + (f(i,k)-f(i,k-1))*xx/(x(i,k)-x(i,k-1))
      ENDIF
      IF (k+1 == plevp) THEN
        psi3 = f(i,plevp)
      ELSE
        psi3 = f(i,k+1)                                             &
             - (f(i,k+2)-f(i,k+1))*(dx-xx)/(x(i,k+2)-x(i,k+1))
      ENDIF
      psim = medan(psi1, psi2, psi3)
      cfnew = medan(cfint, psi1, psim)
      !        if (abs(cfnew-cfint)/(abs(cfnew)+abs(cfint)+1.e-36)
      !        $        .gt..03) then
      !        write (nerr,*) ' cfintv limiting important ', cfint, cfnew
      !        endif
      psistar(i) = cfnew
      !
    END DO

  END SUBROUTINE cfintv

!!============================================================================

  SUBROUTINE getf2 (npts, xw, vel, deltat, depflag, idep, idom, xwd)
    !
    !     Phil Rasch, NCAR, 1998
    !

    !     calculate the wall departure point and intervals

    !     xw are coords of cell walls
    !     vel are velocites at cell walls
    !     deltat is time step

    !     assume number of cell centers is npts-1
    !     number of cell walls is npts


    !     .....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
    !     ....velw1.....velw2.....velw3.....velw4.....velw5.....velw6

    INTEGER, INTENT(in) :: npts         ! number of cell walls
    REAL(wp), INTENT(in) :: vel(npts)   ! velocity at cell walls
    REAL(wp), INTENT(in) :: xw(npts)    ! coords of cell walls
    INTEGER, INTENT(in) :: depflag      ! even or uneven grid flag
    REAL(wp), INTENT(in) :: deltat      ! time step

    INTEGER, INTENT(out)  :: idep(npts)  ! departure point interval
    INTEGER, INTENT(out)  :: idom(npts)  ! flag the wrapping of dep points
    REAL(wp), INTENT(out) :: xwd(npts)   ! departure point

    INTEGER :: i, ii, j, k, l

    REAL(wp) :: xint, tmp, dxm

    !     xwd is the departure point for a
    !     an arrival point originating at the wall
    !
    !     idom flags whether the departure point is outside and
    !     greater than the max value of the interval (1), or
    !     less than the min value of the interval (-1)

    xint = xw(npts)-xw(1)
    DO i = 1,npts
      xwd(i) = xw(i)-vel(i)*deltat
      idom(i) = 0
      IF (xwd(i) < xw(1)) THEN
        idom(i) = -1
        tmp = xwd(i)+xint
        xwd(i) = tmp
      ENDIF
      IF (xwd(i) >= xw(npts)) THEN
        idom(i) = 1
        tmp = xwd(i)-xint
        xwd(i) = tmp
        ! take care of potential round off error
        IF (xwd(i) < xw(1)) xwd(i) = xw(1) 
      ENDIF
    END DO

    IF (depflag == 1) THEN
      ! guess the interval assuming x is distributed uniformly
      dxm = xint/(npts-1)
      DO ii = 1,npts
        i = INT((xwd(ii)-xw(1))/dxm)+1
        IF ((xwd(ii)-xw(i))*(xw(i+1)-xwd(ii)) >= 0) THEN
          idep(ii) = i
        ELSE IF (xwd(ii) < xw(i)) THEN
          idep(ii) = i-1
        ELSE IF (xwd(ii) > xw(i+1)) THEN
          idep(ii) = i+1
        ENDIF
      END DO
    ELSE IF (depflag == 2) THEN
      ! guess the interval using a precalculated mapping
      ! from an even distribution to
      ! the actual grid
      dxm = xint/(njeven-1)
      ! set xwd(npts) manually for y-advection:
      ! meridional velocity at poles=0 and: 
      ! IF (xwd(i) >= xw(npts)) THEN ... (see above)
      xwd(npts) = xw(1)
      DO ii = 1,npts
        j = INT((xwd(ii)-xw(1))/dxm)+1
        i = jeven(j)
        IF ((xwd(ii)-xw(i))*(xw(i+1)-xwd(ii)) >= 0) THEN
          idep(ii) = i
        ELSE IF (xwd(ii) < xw(i)) THEN
          idep(ii) = i-1
        ELSE IF (xwd(ii) > xw(i+1)) THEN
          idep(ii) = i+1
        ENDIF
      END DO
    ELSE
      WRITE (nerr,*) ' cfint1x: illegal value for depflag ', depflag
      CALL finish('getf2','illegal value for depflag')
    ENDIF

    ! now check to see whether we actually found it
    j = 0
    DO ii = 1,npts-1
      !         write (nerr,*) ' looking for interval for point ', xwd(ii)
      i = idep(ii)
      !         write (nerr,*) ' we think it is in interval ', i,
      !     $        ' defined by endpoints ', xw(i), xw(i+1)
      IF ((xwd(ii)-xw(i))*(xw(i+1)-xwd(ii)) < 0) THEN
        j = ii
        !            write (nerr,*) ' problems ', xw(i), xwd(ii), xw(i+1)
        !            call abort
      ENDIF
    END DO
    IF (j /= 0) THEN
      ii = j
      i = idep(ii)
      WRITE (nerr,*) ' getf2: somethings wrong at dep calc ', ii, i
      WRITE (nerr,*) ' depflag is ', depflag
      WRITE (nerr,*) ' dep pt xwd(ii) ', xwd(ii)
      WRITE (nerr,*) ' we thought it should be between ',xw(i),' and ', xw(i+1)
      DO k = 1,npts
        WRITE (nerr,*) ' wall points are ', k, xw(k)
      END DO
      WRITE (nerr,*) ' dxm, xint ', dxm, xint, xw(1), xw(npts)
      k = INT((xwd(ii)-xw(1))/dxm)+1
      WRITE (nerr,*) ' xwd(ii) corresponds to even index ', k
      l = jeven(k)
      WRITE (nerr,*) ' even index ', k, ' corresponds to wall int ', l
      CALL finish ('getf2')
    ENDIF

  END SUBROUTINE getf2

!!============================================================================

  SUBROUTINE getys2 (vw, vwext, i)

    !     move the velocities from normal lon/lat grid to the rotated grid
    
    REAL(wp), INTENT(in) :: vw(nplond, nplat+1)
    INTEGER, INTENT(in)  :: i

    REAL(wp), INTENT(out) :: vwext(nplatext)

    INTEGER :: ii, j

    ii = MOD(i+nplon/2-1,nplon) + 1

!dir$ ivdep
    DO j = 1,nplat
      vwext(j) = vw(i,j)
      vwext(2*nplat-j+2) = -vw(ii,j)
    END DO

    vwext(nplat+1) = vw(i,nplat+1)
    vwext(2*nplat+1) = vwext(1)

  END SUBROUTINE getys2

!!============================================================================

  SUBROUTINE getflxa2 (npts,  nj,    xw,           &
                       psi,   idimpsi,             &
                       flux,  idimflux,            &
                       jdep2, jdom2, vwdp2, idim)

    !
    !     Phil Rasch, NCAR, 1998
    !

    !     calculate the horizontal mass fluxes at cell walls

    !     xw are coords of cell walls
    !     psi is integral of mean tracer densities at cell walls

    !     assume number of cell centers is npts-1
    !     number of cell walls is npts

    !     tracers are updated from the output from this code as
    !     do i = 1,npts-1
    !     q(i) = q(i) - (flux(i+1)-flux(i))/(xw(i+1)-xw(i))
    !     end do

    !.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
    !.... psiw1.....psiw2.....psiw3.....psiw4.....psiw5.....psiw6
    !.... velw1.....velw2.....velw3.....velw4.....velw5.....velw6
    !.........phi1......phi2.......phi3.....phi4.......phi5.......

    INTEGER, INTENT(in)  :: npts            ! number of points
    INTEGER, INTENT(in)  :: idimpsi, idimflux, idim
    INTEGER, INTENT(in)  :: nj
    REAL(wp), INTENT(in) :: xw(npts)        ! cell walls
    REAL(wp), INTENT(in) :: psi(idimpsi,nj) ! integral of densities
    INTEGER, INTENT(in)  :: jdep2(idim,nj)  ! departure point interval
    INTEGER, INTENT(in)  :: jdom2(idim,nj)  ! whether the departure 
                                            ! point wrapped the domain
    REAL(wp), INTENT(in) :: vwdp2(idim,nj)  ! departure point

    REAL(wp), INTENT(out) :: flux(idimflux,nj) ! flux over the time step

    INTEGER :: i, ii

    REAL(wp) :: fdot(nwork,nj)     ! derivative of psi
    REAL(wp) :: psistar(nwork,nj)  ! value of interpolated psi at dep points
    REAL(wp) :: fext(-3:npts+4,nj)
    REAL(wp) :: dx(npts,nj)

    REAL(wp) :: dxh(npts,nj)
    REAL(wp) :: dx2(npts,nj)
    REAL(wp) :: dx3h(npts,nj)

    REAL(wp) :: s(npts,nj)     ! first divided differences at nodes
    REAL(wp) :: sh(npts,nj)    ! first divided differences between nodes
    REAL(wp) :: d(npts,nj)     ! second divided differences at nodes
    REAL(wp) :: dh(npts,nj)    ! second divided differences between nodes
    REAL(wp) :: e(npts,nj)     ! third divided differences at nodes
    REAL(wp) :: eh(npts,nj)    ! third divided differences between nodes

    !  calculate derivatives

    DO ii = 1,nj
!CDIR EXPAND
      CALL cfdot1dp2 (xw, psi(:,ii), npts, fdot(:,ii),  dxh, dx2,  dx3h, &
                      s, sh, d, dh, e, eh, ii, nj)
    END DO

    !  interpolate to departure points

    DO ii = 1,nj 
!CDIR EXPAND
      CALL cfint1x2 (xw,   psi(:,ii), fdot(:,ii), npts,            &
                     npts, psistar(:,ii),                          &
                     jdep2(:,ii), vwdp2(:,ii), fext, dx, nj, ii  )
    END DO

    !     fluxes

    DO ii = 1,nj
      DO i = 1,npts                            ! periodicity correction
        flux(i,ii) = (psi(i,ii)-psistar(i,ii)) - jdom2(i,ii)*psi(npts,ii)
      END DO
    END DO

  CONTAINS

    SUBROUTINE cfdot1dp2 (x, f, npts, fdot, dxh, dx2, dx3h, &
                          s, sh, d, dh, e, eh, j, nj)
      !
      !     Phil Rasch, NCAR, 1998
      !

      !     calculate the derivative for the interpolating polynomial
      !     periodic BCs
      !     for extending the values,
      !     we assume x(1)-x(0) = x(npts)-x(npts-1),
      !     f(1)-f(0) = f(npts)-f(npts-1)

      INTEGER, INTENT(in)  :: j, nj
      INTEGER, INTENT(in)  :: npts
      REAL(wp), INTENT(in) :: x(npts)        ! nodes where data are
      REAL(wp), INTENT(in) :: f(npts)        ! values of data

      REAL(wp), INTENT(out) :: fdot(npts)    ! derivative of f at x

      !     assumed variable distribution
      !     x1.......x2.......x3.......x4.......x5.......x6     1,npts points
      !     f1.......f2.......f3.......f4.......f5.......f6     1,npts points
      !     ...sh1.......sh2......sh3......sh4......sh5....     1,npm1 points
      !     .........d2.......d3.......d4.......d5.........     2,npm1 points
      !     .........s2.......s3.......s4.......s5.........     2,npm1 points
      !     .............dh2......dh3......dh4.............     2,npm2 points
      !     .............eh2......eh3......eh4.............     2,npm2 points
      !     ..................e3.......e4..................     3,npm2 points
      !     .................ppl3......ppl4................     3,npm2 points
      !     .................ppr3......ppr4................     3,npm2 points
      !     .................t3........t4..................     3,npm2 points
      !     ................fdot3.....fdot4................     3,npm2 points


      !     work variables

      INTEGER :: i, npm1, npm2

      REAL(wp) :: s(npts,nj)     ! first divided differences at nodes
      REAL(wp) :: sh(npts,nj)    ! first divided differences between nodes
      REAL(wp) :: d(npts,nj)     ! second divided differences at nodes
      REAL(wp) :: dh(npts,nj)    ! second divided differences between nodes
      REAL(wp) :: e(npts,nj)     ! third divided differences at nodes
      REAL(wp) :: eh(npts,nj)    ! third divided differences between nodes

      REAL(wp) :: pp           ! p prime
      REAL(wp) :: ppl          ! p prime on left
      REAL(wp) :: ppr          ! p prime on right

      REAL(wp) :: qpl, qpr, ttt, t, tmin, tmax

      INTEGER :: ii, im1, ip1

      REAL(wp) :: dxh(npts,nj), dx2(npts,nj), dx3h(npts,nj)

      REAL(wp) :: fd1, fd2

      npm1 = npts-1
      npm2 = npts-2

      !     a bunch of divided differences are related functions
      DO i = 1,npm1
        dxh(i,j) = x(i+1)-x(i)
        sh(i,j) = (f(i+1)-f(i))/dxh(i,j)
      END DO

      DO i = 2,npm1
        dx2(i,j) = (dxh(i,j)+dxh(i-1,j))     ! equals x(i+1)-x(i-1)
        d(i,j) = (sh(i,j)-sh(i-1,j))/dx2(i,j)
        s(i,j) = minmod(sh(i,j),sh(i-1,j))
      END DO

      DO i = 2,npm2
        dx3h(i,j) = (x(i+2)-x(i-1))
        eh(i,j) = (d(i+1,j)-d(i,j))/dx3h(i,j)
        dh(i,j) = minmod(d(i,j),d(i+1,j))
      END DO

      DO i = 3,npm2
        e(i,j) = minmod(eh(i,j),eh(i-1,j))
      END DO

      sh(npts,j) = sh(1,j)
      dxh(npts,j) = dxh(1,j)
      dx2(1,j) = dxh(1,j) + dxh(npm1,j)
      dx2(npts,j) = dx2(1,j)
      dx3h(npts,j) = dxh(npm1,j) + dxh(npts,j) + dxh(2,j)
      dx3h(1,j) = dx3h(npts,j)
      dx3h(npm1,j) = dxh(npts-2,j) + dxh(npm1,j) + dxh(npts,j)
      d(npts,j) = (sh(npts,j)-sh(npm1,j))/dx2(npts,j)
      d(1,j) = d(npts,j)
      s(npts,j) = minmod(sh(npts,j),sh(npm1,j))
      s(1,j) = s(npts,j)
      eh(npm1,j) = (d(npts,j)-d(npm1,j))/dx3h(npm1,j)
      dh(npm1,j) = minmod(d(npm1,j),d(npts,j))
      eh(npts,j) = (d(2,j)-d(npts,j))/dx3h(npts,j)
      eh(1,j) = eh(npts,j)
      dh(npts,j) = minmod(d(npts,j),d(2,j))
      dh(1,j) = dh(npts,j)
      e(2,j) = minmod(eh(2,j), eh(1,j))
      e(1,j) = minmod(eh(1,j), eh(npm1,j))
      e(npts,j) = e(1,j)
      e(npm1,j) = minmod(eh(npm1,j),eh(npts-2,j))

      !     now calculate the derivatives in interior
      DO i = 3,npm2

        im1 = i-1
        ip1 = i+1

        !        p prime on left side (what huynh calls i-0.5)
        ppl=sh(im1,j) + dh(im1,j)*dxh(im1,j)
        !        p prime on right side (what huynh calls i+0.5)
        ppr=sh(i,j)   - dh(i,j)  *dxh(i,j)
        t = minmod(ppl,ppr)
        !        tmax = sign(1.,t)*max(3.*abs(s(i,j)),1.5*abs(t))
        !        fdot(i) = 0.5*(ppl+ppr) !  derivative M3-A
        !        fdot(i) = minmod(fdot(i),tmax)


        !        derivative from parabola thru f(im1), f(i), and f(ip1)
        pp = sh(im1,j) + d(i,j)*dxh(im1,j)
        fd1 = pp - dxh(im1,j)*dxh(i,j)                       &
             *(eh(im1,j)*dx2(ip1,j) + eh(i,j)*dx2(im1,j))    &
             /(dx2(ip1,j)+dx2(im1,j))      !(x(i+2)-x(i-2))
        !                         /(x(i+2)-x(i-2))

        qpl = sh(im1,j)                                              &
             + (dxh(im1,j))*minmod(d(im1,j) + e(im1,j)*dx2(im1,j),   &
             d(i,j)   - e(i,j)  *dxh(i,j))
        qpr = sh(i,j)                                                &
             + (dxh(i,j))  *minmod(d(i,j)   + e(i,j)  *dxh(i,j),     &
             d(ip1,j) - e(ip1,j)*dx2(ip1,j))
        ttt = minmod(qpl, qpr)
        tmin = MIN(0.0_wp,3*s(i,j),1.5_wp*t,ttt)
        tmax = MAX(0.0_wp,3*s(i,j),1.5_wp*t,ttt)
        fd2 = medan(fd1, qpl, qpr)
        fdot(i) = fd2 + minmod(tmin-fd2, tmax-fd2)

      END DO

!CDIR UNROLL=8
      DO ii = -1, 2
        i = MOD(npts+ii-1,npts) + 1
        im1 = MOD(npts+i-3,npts-1) + 1
        ip1 = MOD(npts+i-1,npts-1) + 1

        !        p prime on left side (what huynh calls i-0.5)
        ppl=sh(im1,j) + dh(im1,j)*dxh(im1,j)
        !        p prime on right side (what huynh calls i+0.5)
        ppr=sh(i,j)   - dh(i,j)  *dxh(i,j)
        t = minmod(ppl,ppr)
        !        tmax = sign(1.,t)*max(3.*abs(S(i,j)),1.5*abs(t))
        !        fdot(i) = 0.5*(ppl+ppr) !  derivative M3-A
        !        fdot(i) = minmod(fdot(i),tmax)


        !        derivative from parabola thru f(im1), f(i), and f(ip1)
        pp = sh(im1,j) + d(i,j)*dxh(im1,j)
        fd1 = pp - dxh(im1,j)*dxh(i,j)                              &
             *(eh(im1,j)*dx2(ip1,j) + eh(i,j)*dx2(im1,j))           &
             /(dx2(ip1,j)+dx2(im1,j))      !(x(i+2)-x(i-2))

        qpl = sh(im1,j)                                              &
             + (dxh(im1,j))*minmod(d(im1,j) + e(im1,j)*dx2(im1,j),   &
             d(i,j)   - e(i,j)  *dxh(i,j))
        qpr = sh(i,j)                                                &
             + (dxh(i,j))  *minmod(d(i,j)   + e(i,j)  *dxh(i,j),     &
             d(ip1,j) - e(ip1,j)*dx2(ip1,j))
        ttt = minmod(qpl, qpr)
        tmin = MIN(0.0_wp,3*s(i,j),1.5_wp*t,ttt)
        tmax = MAX(0.0_wp,3*s(i,j),1.5_wp*t,ttt)
        fd2 = medan(fd1, qpl, qpr)
        fdot(i) = fd2 + minmod(tmin-fd2, tmax-fd2)

      END DO

    END SUBROUTINE cfdot1dp2

    SUBROUTINE cfint1x2 (x, f, fdot, ndata, npts, psistar,   &
                         jdep2, vwdp2, fext, dx, nj ,j)
      !
      !     Phil Rasch, NCAR, 1998
      !

      !     interpolate a function specified by (x, f, fdot)
      !     at points vwdp2
      !     return values of psistar

      !     assume periodic in the first dimension

      !     input
      INTEGER, INTENT(in)  :: ndata
      INTEGER, INTENT(in)  :: j, nj
      REAL(wp), INTENT(in) :: x(ndata)    ! coords of function
      REAL(wp), INTENT(in) :: f(ndata)    ! value of function at x
      REAL(wp), INTENT(in) :: fdot(ndata) ! deriv of function at x
      INTEGER, INTENT(in)  :: npts        ! number of points to interpolate to
      INTEGER, INTENT(in)  :: jdep2(npts) ! interval where each interpolation 
                                          ! is done
      REAL(wp), INTENT(in) :: vwdp2(npts) ! coord to construct the interpolant
                                          ! for
      REAL(wp), INTENT(out) :: fext(-3:npts+4,nj)
      REAL(wp), INTENT(out) :: dx(npts,nj)

      REAL(wp), INTENT(out) :: psistar(npts)  ! interpolated value

      INTEGER :: i, ii
      REAL(wp) :: s, c2, c3, xx

      REAL(wp) :: psi1, psi2, psi3, psim
      REAL(wp) :: cfnew, cfint
      INTEGER :: ip1, im1

      !     construct an extended version of the grid
      DO ii = 1,ndata-1
        fext(ii,j) = f(ii)
        dx(ii,j) = (x(ii+1)-x(ii))
      END DO
      fext(ndata,j) = f(ndata)
      dx(ndata,j) = dx(1,j)

!  next = 3
!dir$ ivdep
!CDIR UNROLL=8
!OCL NOVREC
      DO ii = 0,3
        fext(ii+ndata+1,j) = f(ii+2)  + (f(ndata)-f(1))
        fext(-ii,j)     = f(ndata-ii-1) - (f(ndata)-f(1))
      END DO

      !     now do the interpolation at each point
      DO ii = 1,npts

        !        i is the departure point interval
        i = jdep2(ii)
        xx = (vwdp2(ii)-x(i))
        ip1 = MOD(ndata+i-2,ndata-1) + 2
        !        ip2 = mod(ndata+i-1,ndata-1) + 2
        im1 = MOD(ndata+i-4,ndata-1) + 2
        !        now interpolate
        s = (f(ip1)-f(i))/dx(i,j)
        c2 = (3*s-2*fdot(i)-fdot(ip1))/dx(i,j)
        c3 = (fdot(i)+fdot(ip1)-2*s)/dx(i,j)**2
        !        fxdot =  (3*c3*xx + 2*c2)*xx + fdot(i)
        !        fxdd = 6*c3*xx + 2*c2
        cfint = ((c3*xx + c2)*xx + fdot(i))*xx + f(i)
        !
        !        leonards limiter
        psi1 = fext(i,j)   + (fext(i+1,j)-fext(i,j))*xx/dx(i,j)
        psi2 = fext(i,j)   + (fext(i,j)-fext(i-1,j))*xx/dx(im1,j)
        psi3 = fext(i+1,j) - (fext(i+2,j)-fext(i+1,j))*(dx(i,j)-xx)/dx(ip1,j)
        psim = medan(psi1, psi2, psi3)
        cfnew = medan(cfint, psi1, psim)
        cfint = cfnew
        !
        psistar(ii) = cfint

      END DO

    END SUBROUTINE cfint1x2

  END SUBROUTINE getflxa2

!!============================================================================

  SUBROUTINE getflxv ( phi, flux, dp, kdep2)
    !
    !     Phil Rasch, NCAR, 1998
    !
    !
    !     calculate the vertical mass fluxes
    !
    !     zw are coords of cell walls
    !     phi are mean tracer densities in cell (nominally at cell centers)
    !
    !     tracers are updated from the output from this code as
    !     do k = 1,plev
    !     do i = 1,plon
    !     q(i,k) = q(i,k) - (flux(i,k+1)-flux(i,k))/(zw(i,k+1)-zw(i,k))
    !     end do
    !     end do
    !
    !     .....zw1........zw2........zw3........zw4........zw5........zw6
    !     .... psiw1......psiw2......psiw3......psiw4......psiw5......psiw6
    !     .... velw1......velw2......velw3......velw4......velw5......velw6
    !     ...........phi1.......phi2........phi3......phi4........phi5.....

    REAL(wp), INTENT(in) :: phi(plon,plev) ! densities
    REAL(wp), INTENT(in) :: dp(plon,plev)    ! departure point of cell walls
    INTEGER, INTENT(in)  :: kdep2(plon,plev) ! interval of departure point

    REAL(wp), INTENT(out) :: flux(plon,plevp)

    REAL(wp) :: psi(plon,plevp)
    REAL(wp) :: fdot(plon,plevp)
    !     real fxdot(plond)
    !     real fxdd(plond)

    REAL(wp) :: psistar(nplond)
    INTEGER :: i, k

    DO i = 1,plon
      !        integral of phi
      psi(i,1) = 0.0_wp
      !        fluxes at boundaries
      flux(i,1) = 0.0_wp
      flux(i,plevp) = 0.0_wp
    END DO

    !     integral function
    DO k = 2,plevp
      DO i = 1,plon
        psi(i,k) = phi(i,k-1)*(zw(i,k)-zw(i,k-1)) + psi(i,k-1)
      END DO
    END DO

    !     calculate the derivatives for the interpolating polynomial
    CALL cfdotmc (zw, psi, fdot)

    DO k = 2,plev
      CALL cfintv(plevp, plon, zw, psi, fdot, dp(:,k),        &
           psistar, kdep2(:,k))
      DO i = 1,plon
        flux(i,k) = (psi(i,k)-psistar(i))
      END DO
    END DO

  END SUBROUTINE getflxv

!!============================================================================

  SUBROUTINE hadvn (den, deltat, n)

    ! horizontal advection algorithm

    REAL(wp), INTENT(in) :: deltat
    INTEGER, INTENT(in)  :: n    ! time step index used for toggling operators

    REAL(wp), INTENT(inout) :: den(nplon,nplat) ! density on input/output

    INTEGER :: i, j
    REAL(wp) :: denx(nplond,nplat), deny(nplond,nplat), dent(nplon,nplat)
    REAL(wp) :: xtend(nplond,nplat), ytend(nplond,nplat)

    REAL(wp) :: yflux(nplond,nplat+1)        ! fluxes on y cell walls
    REAL(wp) :: xflux(nplond,nplat)          ! fluxes on x cell walls

    ! .....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
    ! .... psiw1.....psiw2.....psiw3.....psiw4.....psiw5.....psiw6
    ! .... velw1.....velw2.....velw3.....velw4.....velw5.....velw6
    ! ..........xc1.......xc2........xc3......xc4........xc5.......
    ! ..........dx1.......dx2........dx3......dx4........dx5.......
    ! .........phi1......phi2.......phi3.....phi4.......phi5.......


    !        toggle order of x and y operations each time step
    IF (MOD(n,2) == 1) THEN

      !           y advection
      CALL yupdate(den, deny, ytend, yflux)


      DO j = 1,nplat
!!argchk        DO i = 1,nplon+1
        DO i = 1,nplon
          !                 very bad advective approximation
          dent(i,j) = deny(i,j) + den(i,j)                    &
               *deltat*(vw(i,j+1)-vw(i,j))/dmu(j)
        END DO
      END DO

      !           x advection
      CALL xupdate(dent, denx, xtend, xflux)

      DO j = 1,nplat
!!argchk        DO i = 1,nplon+1
        DO i = 1,nplon
          den(i,j) = deny(i,j) + xtend(i,j)
        END DO
      END DO


    ELSE

      !           x advection
      CALL xupdate(den, denx, xtend, xflux)

      DO j = 1,nplat
        DO i = 1,nplon
          !                 very bad advective approximation
          dent(i,j) = denx(i,j) + den(i,j)                    &
               *deltat*(uw(i+1,j)-uw(i,j))/dx
        END DO
!!argchk        dent(nplon+1,j) = dent(1,j)
      END DO


      !           y advection
      CALL yupdate(dent, deny, ytend, yflux)

      DO j = 1,nplat
!!argchk        DO i = 1,nplon+1
        DO i = 1,nplon
         den(i,j) = denx(i,j) + ytend(i,j)
        END DO
      END DO

    ENDIF

    !        filter the field near the poles
    CALL filter (den)

  END SUBROUTINE hadvn

!!============================================================================

  SUBROUTINE seth (u3, v3, deltat)
    !
    !     Phil Rasch, NCAR, 1998
    !
    !     calculate velocity information for all tracers

    !     calculate the departuture points (uwdp, vwdp)
    !     their intervals (idep, jdep)
    !     and
    !     whether they wrap outside the domain

    !.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
    !.... velw1.....velw2.....velw3.....velw4.....velw5.....velw6
    !..........xc1.......xc2........xc3......xc4........xc5.......
    !..........dx1.......dx2........dx3......dx4........dx5.......

    REAL(wp), INTENT(in) :: u3(nplon,nplat)
    REAL(wp), INTENT(in) :: v3(nplon,nplat)
    REAL(wp), INTENT(in) :: deltat

    INTEGER :: i, j

    REAL(wp) :: vwext(nplatext)              ! work variable
    REAL(wp) :: c1, c2, c3                   ! work variable

    c1 = 0.5_wp/arad

    !        map from cell centered velocities on ccm grid
    !        to cell wall velocities
    !        convert from traditional velocities
    !        to lamdadot and mudot (mu=sin(phi))
    !        or      from traditional velocities to lamdadot and phidot

    !------------------------------

    ! wall centered meridional velocities

    ! a) set meridional velocity at pole walls to zero
    DO i = 1,nplon+1
      vw(i,1) = 0.0_wp
      vw(i,nplat+1) = 0.0_wp
    END DO
    ! b) compute meridional velocity at walls between cells
    DO j = 2,nplat
      c2 = COS(yw(j))*c1
      DO i = 1,nplon
        vw(i,j) = (v3(i,j-1)+v3(i,j))*c2
      END DO
      vw(nplon+1,j) = vw(1,j)
    END DO

    !------------------------------

    ! wall centered zonal velocities

    DO j = 1,nplat
      c3 = c1/COS((yw(j+1)+yw(j))*0.5_wp)
      !    a) compute zonal velocity into first cell
      uw(1,j) =  (u3(nplon,j)+u3(1,j))*c3
      !    b) compute zonal velocity into following cells
      DO i = 2,nplon
        uw(i,j) = (u3(i-1,j)+u3(i,j))*c3
      END DO
      !    c) close zonal cycle
      uw(nplon+1,j) = uw(1,j)
    END DO

    !------------------------------

    !        take the wall centered winds,
    !        and calculate their departure interval and increment

    !        first do the y cell walls
    DO i = 1,nplon/2

      !           put v winds into the reshaped arrays
      CALL getys2 (vw, vwext, i)

      !           calculate the departure point for particles
      !            located at each cell and the interval it occurs in
      !
      CALL getf2 (nplatext, ywext, vwext, deltat,                 &
           2,                                                     &  
           jdep(:,i), jdom(:,i), vwdp(:,i))

    END DO

    !        now do the x cell walls
    DO j = 1,nplat

      !           calculate the departure point for particles
      !           located at each cell and the interval it occurs in
      !
      CALL getf2 (nplon+1, xw, uw(:,j), deltat,                   &
           1,                                                     &  
           idep(:,j), idom(:,j), uwdp(:,j))

    END DO

  END SUBROUTINE seth

!!============================================================================

  SUBROUTINE setv (etadot, deltat, j)
    !
    !     Phil Rasch, NCAR, 1998
    !

    !     calculate the departure points of vertical cell walls
    !     and their intervals given vertical velocities and time step

    !     input

    !***NOTE THAT THIS IS NOT IN STANDARD ORDER, OR DIMENSIONED PLOND!!

    REAL(wp), INTENT(in) :: etadot(plon,plevp)           ! vertical velocities
    REAL(wp), INTENT(in) :: deltat                       ! time step
    INTEGER,  INTENT(in) :: j                            ! latitude

    INTEGER :: i, jj, k, kk
    REAL(wp) :: xins


    DO i = 1,plon
      wwdp(i,1,j)     = 0.0_wp
!!      wwdp(i,plevp,j) = 0.0_wp
    END DO

    !        calculate departure points of cell walls
    DO k = 2,plev
      DO i = 1,plon
        !
        !              the departure point
        wwdp(i,k,j) = zw(i,k)-etadot(i,k)*deltat

        !              guess the interval
        !              using a precalculated mapping from an even distribution
        !              to the actual grid, then refine it
        !              xins = medan(zw(i,1), wwdp(i,k,j), zw(i,plevp))
        xins = MAX(zw(i,1), MIN(wwdp(i,k,j),zw(i,plevp)))
        kdep(i,k,j) = 0
        jj = INT((xins-zw(i,1))/dzmin)+1
        kk = keven(jj)
        IF (xins < zw(i,kk)) THEN
          kdep(i,k,j) = kk-1
        ELSE IF (xins > zw(i,kk+1)) THEN
          kdep(i,k,j) = kk+1
        ELSE
          kdep(i,k,j) = kk
        ENDIF

      END DO
    END DO

  END SUBROUTINE setv

!!============================================================================

  SUBROUTINE vadvn (dn, ilat)

    ! vertical advection algorithm

    INTEGER, INTENT(in) :: ilat

    REAL(wp), INTENT(inout) :: dn (plon, plat, plev, pcnst)

    REAL(wp) :: dnxz(plon,plev), flux(plon,plevp)
    INTEGER :: i, k, m

    !        move 3d fields to standard lon/height slice
    DO m = 1,pcnst
      DO k = 1,plev
        DO i = 1,plon
          dnxz(i,k) = dn(i,ilat,k,m)
        END DO
!!        dnxz(plon+1,k) = dnxz(1,k)
      END DO

      !        calculate the mass moving past each face
      CALL getflxv( dnxz, flux, wwdp(:,:,ilat), kdep(:,:,ilat))

      !        update the species array for vertical transport
      DO k = 1,plev
        DO i = 1,plon
          dn(i,ilat,k,m) = dn(i,ilat,k,m)           &
               - (flux(i,k+1)-flux(i,k))            &
               /(zw(i,k+1)-zw(i,k))
        END DO
!        dn(plon+1,ilat,k,m) = dn(1,ilat,k,m)
      END DO                         ! do k = 1,plev
    ENDDO                     ! do m=1,pcnst

  END SUBROUTINE vadvn

!!============================================================================

  SUBROUTINE xupdate(den, denx, xtend, xflux)
    !
    !     Phil Rasch, NCAR, 1998
    !

    !     calculate the new densities after advection in y direction

    REAL(wp), INTENT(in) :: den(nplon,nplat)

    REAL(wp), INTENT(out) :: denx(nplond,nplat)
    REAL(wp), INTENT(out) :: xtend(nplond,nplat)
    REAL(wp), INTENT(out) :: xflux(nplond,nplat)

    INTEGER :: i, j

    REAL(wp) :: deni(nplond,nplat)
    REAL(wp) :: tmp


    !     construct x update and calculate tendencies

    !     calculate the x integral of the density
    DO j = 1,nplat
      deni(1,j) = 0.0_wp
    END DO

    DO i = 2,nplon+1
      DO j = 1,nplat
        deni(i,j) = deni(i-1,j) + den(i-1,j)*(xw(i)-xw(i-1))
      END DO
    END DO

    !     calculate the x fluxes

    CALL getflxa2(nplond, nplat, xw,          &
                  deni,   nplond,             &
                  xflux,  nplond,             &
                  idep, idom, uwdp, nplond )

    !     add wrap point
    DO j = 1,nplat
      xflux(nplon+1,j) = xflux(1,j)
    END DO

    !     update the fields
    DO j = 1,nplat
      DO i = 1,nplon
        tmp = - (xflux(i+1,j) - xflux(i,j))/dx
        xtend(i,j) = tmp
        denx(i,j) = den(i,j) + tmp
      END DO
    END DO
    DO j = 1,nplat
      denx(nplon+1,j) = denx(1,j)
      xtend(nplon+1,j) = xtend(1,j)
    END DO

  END SUBROUTINE xupdate

!!============================================================================

  SUBROUTINE yupdate(den, deny, ytend, yflux)

    !     input
    REAL(wp), INTENT(in) :: den(nplon,nplat)                ! input density

    !     output
    REAL(wp), INTENT(out) :: deny(nplond,nplat)             ! output density
    REAL(wp), INTENT(out) :: ytend(nplond,nplat)
    REAL(wp), INTENT(out) :: yflux(nplond,nplat+1)

    INTEGER :: i, j 

    REAL(wp) :: yfluxext(nplatext,nplon/2)
    REAL(wp) :: den2(nplatext,nplon/2)
    REAL(wp) :: deni(nplatext,nplon/2)
    REAL(wp) :: tmp

!CDIR NODEP
!OCL NOVREC
    DO i = 1,nplon/2
!CDIR EXPAND
      !     move densities to rotated grid
      CALL getyslice (den, den2(:,i), i)
    END DO

    !     compute the y integral of the densities
    DO i = 1,nplon/2
      deni(1,i) = 0.0_wp
    END DO
    DO j = 2,nplatext
      DO i = 1,nplon/2
        deni(j,i) = deni(j-1,i) + den2(j-1,i)*(ywext(j)-ywext(j-1))
      END DO
    END DO

    !     calculate fluxes
    CALL getflxa2(nplatext, nplon/2, ywext, deni, nplatext,     &
                  yfluxext, nplatext,                           &
                  jdep, jdom, vwdp, nplatext)

    !        move fluxes back to regular lon/lat grid
!CDIR NODEP
!OCL NOVREC
    DO i = 1,nplon/2
!CDIR EXPAND
      CALL putyslice(yfluxext(:,i), yflux, i)
    END DO

    !     add wrap point
!dir$ ivdep
    DO j = 1,nplat+1
      yflux(nplon+1,j) = yflux(1,j)
    END DO


    !     calculate tendencies and update densities
    DO j = 1,nplat
      !        *works            yc = (yw(j+1)+yw(j))/2.
      !        *works            dy = cos(yc)*(yw(j+1)-yw(j))

      DO i = 1,nplon
        tmp = - (yflux(i,j+1)- yflux(i,j))/dmu(j)
        !           $                                      )/dy
        !           tmp = - (  yflux(i,j+1)*cos(yw(j+1))
        !           $              - yflux(i,j)*cos(yw(j))
        !           $             )/dy
        ytend(i,j) = tmp
        deny(i,j) = den(i,j) + tmp
      END DO
    END DO

    DO j = 1,nplat
      deny(nplon+1,j)  = deny(1,j)
      ytend(nplon+1,j) = ytend(1,j)
    END DO

  CONTAINS

    SUBROUTINE getyslice (den, denext, i)
      !
      !     Phil Rasch, NCAR, 1998
      !

      !     move the densities from a normal lon/lat grid to the rotated grid

      !     input
      REAL(wp), INTENT(in) :: den(nplon, nplat)         ! input densities
      INTEGER, INTENT(in)  :: i                         ! longitude to extrace

      !     output
      REAL(wp), INTENT(out) :: denext(nplatext)         ! rotated densities

      INTEGER :: ii, j

      ii = MOD(i+nplon/2-1,nplon) + 1     ! x index on other side of globe
!cdir$ ivdep
      DO j = 1,nplat
        denext(j) = den(i,j)
        denext(2*nplat-j+1) = den(ii,j)
      END DO

      denext(2*nplat+1) = denext(1)

    END SUBROUTINE getyslice

    SUBROUTINE putyslice(yfluxext, yflux, i)
      !
      !     Phil Rasch, NCAR, 1998
      !

      !     map the fluxes from the rotated grid back to the lon/lat grid

      INTEGER, INTENT(in)  :: i
      REAL(wp), INTENT(in) :: yfluxext(nplatext)

      REAL(wp), INTENT(inout) :: yflux(nplond,nplat+1)

      INTEGER :: ii, j

      ii = MOD(i+nplon/2-1,nplon) + 1

      DO j = 2,nplat
        yflux(i,j) = yfluxext(j)
        yflux(ii,j) = -yfluxext(2*nplat-j+2)
      END DO
      yflux(i,nplat+1) =  yfluxext(nplat+1)
      yflux(ii,nplat+1) = -yfluxext(nplat+1)
      yflux(i,1) =  yfluxext(1)
      yflux(ii,1) =  -yfluxext(1)

    END SUBROUTINE putyslice

  END SUBROUTINE yupdate

!!============================================================================
!!============================================================================

  SUBROUTINE filter (field)

    !     Phil Rasch, NCAR, 1998

    REAL(wp), INTENT(inout) :: field(nplon,nplat)

    INTEGER :: i, j

!!    REAL(wp) :: temp(nplond)
    REAL(wp) :: temp(nplon)

    !     update solution

    DO j = 1,nplat
      IF (nptslat(j) > 1) THEN
        DO i = 1,nplon
          temp(i) = 0.0_wp
        END DO
        DO i = 1,nplon
          temp(bini(i,j)) = temp(bini(i,j)) + field(i,j)
        END DO
        DO i = 1,nplon
          field(i,j) = temp(bini(i,j))/nptslat(j)
        END DO
!!        field(nplon+1,j) = field(1,j)
      ENDIF
    END DO

  END SUBROUTINE filter

!!============================================================================

  SUBROUTINE pole_filter(ptte, pqte, klat)
    !
    !      FILTER SETTINGS FROM *INIFILTER* ARE USED
    !
    !
    !      INPUT:  PTTE, PQTE  -  UNFILTERED TENDENCIES
    !
    !      OUTPUT: PTTE, PQTE  -  FILTERED TENDENCIES
    !

    USE mo_mpi,           ONLY: p_isend, p_irecv, p_wait

    INTEGER, INTENT(in) :: klat
    REAL(wp), INTENT(inout) :: ptte(plon,plev), pqte(plon,plev)
    
    REAL(wp) :: zont(nplon,plev), zonq(nplon,plev)
    REAL(wp) :: zt(nplon,plev), zq(nplon,plev)
    INTEGER  :: jl, jk, jlg, i, k, nsrc

    REAL(wp) :: bufs (2* plon*plev) ! send buffer
    REAL(wp) :: bufr (2*nplon*plev) ! receive buffer

    IF(nptslat_l(klat) == 1) RETURN  ! no filtering required

    if (ldc% nprocb > 1) then
      !---------------------------------
      ! copy ptte, pqte into send buffer
      !---------------------------------
      bufs (:plon*plev   ) = reshape (ptte, (/plon*plev/))
      bufs ( plon*plev+1:) = reshape (pqte, (/plon*plev/))
      !------------------------------
      ! loop over all PEs in this row
      !------------------------------
      k = 1
      do i = ldc%spe, ldc%epe
        IF (gdc(i)% set_a /= ldc% set_a) CYCLE
        if (gdc(i)% pe /= ldc% pe) then
          !-----------------------------------------------
          ! initiate unblocked send receive for other pe's
          !-----------------------------------------------
          nsrc = 2 * gdc(i)% nglon * plev          
          CALL p_isend(bufs    ,gdc(i)% pe ,p_tag=1)
          CALL p_irecv(bufr(k:),gdc(i)% pe ,p_tag=1, p_count=nsrc)
          k = k + nsrc
        endif
      end do
      !-------------------------------------------------
      ! copy values from receive buffer to zonal segment
      !-------------------------------------------------
      CALL p_wait
      k = 0
      do i = ldc%spe, ldc%epe
        IF (gdc(i)% set_a /= ldc% set_a) CYCLE
        if (gdc(i)% pe == ldc% pe) then
          !----------------------
          ! direct copy for my pe
          !----------------------
          zont (ldc% glon(klat)+1:ldc% glon(klat)+plon ,:) = ptte(:,:)
          zonq (ldc% glon(klat)+1:ldc% glon(klat)+plon ,:) = pqte(:,:)
        else
          !--------------------------------
          ! copy from buffer for other pe's
          !-------------------------------- 
          nsrc = gdc(i)% nglon * plev
          zont (gdc(i)% glon(klat)+1:gdc(i)% glon(klat)+plon ,:) = &
            reshape (bufr(k+1:k+nsrc), (/gdc(i)% nglon, plev/))
          k = k + nsrc
          zonq (gdc(i)% glon(klat)+1:gdc(i)% glon(klat)+plon ,:) = &
            reshape (bufr(k+1:k+nsrc), (/gdc(i)% nglon, plev/))
          k = k + nsrc
        endif
      end do
    else
      zont(:,:) = ptte(:,:)
      zonq(:,:) = pqte(:,:)
    endif

    zt(:,:) = 0.0_wp
    zq(:,:) = 0.0_wp

    DO jl = 1, nplon      ! redundantly performed by all PEs in a row
      DO jk = 1, plev
        zt(bini_l(jl,klat),jk) = zt(bini_l(jl,klat),jk)+zont(jl,jk)
        zq(bini_l(jl,klat),jk) = zq(bini_l(jl,klat),jk)+zonq(jl,jk)
      END DO
    END DO

    DO jk = 1, plev
      DO jl = 1, plon
        jlg = ldc%glon(klat) + jl
        ptte(jl,jk) = zt(bini_l(jlg,klat),jk)/nptslat_l(klat)
        pqte(jl,jk) = zq(bini_l(jlg,klat),jk)/nptslat_l(klat)
     END DO
   END DO

  END SUBROUTINE pole_filter

!!============================================================================

!---wiso-code

!! This subroutine includes water isotopes and replaces the original pole_filter routine

  SUBROUTINE pole_filter_wiso(ptte, pqte, pwisoqte, klat)
    !
    !      FILTER SETTINGS FROM *INIFILTER* ARE USED
    !
    !
    !      INPUT:  PTTE, PQTE  -  UNFILTERED TENDENCIES
    !
    !      OUTPUT: PTTE, PQTE  -  FILTERED TENDENCIES
    !

    USE mo_mpi,           ONLY: p_isend, p_irecv, p_wait
    USE mo_wiso,          ONLY: lwiso, nwiso

    INTEGER, INTENT(in) :: klat
    REAL(wp), INTENT(inout) :: ptte(plon,plev), pqte(plon,plev), pwisoqte(plon,plev,nwiso)
    
    REAL(wp) :: zont(nplon,plev), zonq(nplon,plev), zonwisoq(nplon,plev,nwiso)
    REAL(wp) :: zt(nplon,plev), zq(nplon,plev), zwisoq(nplon,plev,nwiso)
    INTEGER  :: jl, jk, jlg, i, k, nsrc, jt

    REAL(wp) :: bufs ((2+nwiso)* plon*plev) ! send buffer
    REAL(wp) :: bufr ((2+nwiso)*nplon*plev) ! receive buffer

    IF(nptslat_l(klat) == 1) RETURN  ! no filtering required

    if (ldc% nprocb > 1) then
      !---------------------------------
      ! copy ptte, pqte into send buffer
      !---------------------------------
      bufs (1                   :plon*plev   )      = reshape (ptte, (/plon*plev/))
      bufs (1+plon*plev         :2*plon*plev )      = reshape (pqte, (/plon*plev/))
      IF (lwiso) THEN
        do jt=0,nwiso-1
         bufs (1+(2+jt)*plon*plev:(2+jt+1)*plon*plev) = reshape (pwisoqte(:,:,jt), (/plon*plev/)) 
        end do
      END IF
      !------------------------------
      ! loop over all PEs in this row
      !------------------------------
      k = 1
      do i = ldc%spe, ldc%epe
        IF (gdc(i)% set_a /= ldc% set_a) CYCLE
        if (gdc(i)% pe /= ldc% pe) then
          !-----------------------------------------------
          ! initiate unblocked send receive for other pe's
          !-----------------------------------------------
          nsrc = 2 * gdc(i)% nglon * plev          
          CALL p_isend(bufs    ,gdc(i)% pe ,p_tag=1)
          CALL p_irecv(bufr(k:),gdc(i)% pe ,p_tag=1, p_count=nsrc)
          k = k + nsrc
        endif
      end do
      !-------------------------------------------------
      ! copy values from receive buffer to zonal segment
      !-------------------------------------------------
      CALL p_wait
      k = 0
      do i = ldc%spe, ldc%epe
        IF (gdc(i)% set_a /= ldc% set_a) CYCLE
        if (gdc(i)% pe == ldc% pe) then
          !----------------------
          ! direct copy for my pe
          !----------------------
          zont (ldc% glon(klat)+1:ldc% glon(klat)+plon ,:) = ptte(:,:)
          zonq (ldc% glon(klat)+1:ldc% glon(klat)+plon ,:) = pqte(:,:)
          IF (lwiso) THEN
            do jt=1,nwiso
              zonwisoq (ldc% glon(klat)+1:ldc% glon(klat)+plon ,:,jt) = pwisoqte(:,:,jt)
            end do
          END IF
        else
          !--------------------------------
          ! copy from buffer for other pe's
          !-------------------------------- 
          nsrc = gdc(i)% nglon * plev
          zont (gdc(i)% glon(klat)+1:gdc(i)% glon(klat)+plon ,:) = &
            reshape (bufr(k+1:k+nsrc), (/gdc(i)% nglon, plev/))
          k = k + nsrc
          zonq (gdc(i)% glon(klat)+1:gdc(i)% glon(klat)+plon ,:) = &
            reshape (bufr(k+1:k+nsrc), (/gdc(i)% nglon, plev/))
          k = k + nsrc
          IF (lwiso) THEN
            do jt=1,nwiso
              zonwisoq (gdc(i)% glon(klat)+1:gdc(i)% glon(klat)+plon ,:,jt) = &
                reshape (bufr(k+1:k+nsrc), (/gdc(i)% nglon, plev/))
              k = k + nsrc
            end do
          END IF
        endif
      end do
    else
      zont(:,:) = ptte(:,:)
      zonq(:,:) = pqte(:,:)
      IF (lwiso) THEN
        zonwisoq(:,:,:) = pwisoqte(:,:,:)
      END IF
    endif

    zt(:,:) = 0.0_wp
    zq(:,:) = 0.0_wp
    zwisoq(:,:,:) = 0.0_wp

    DO jl = 1, nplon      ! redundantly performed by all PEs in a row
      DO jk = 1, plev
        zt(bini_l(jl,klat),jk) = zt(bini_l(jl,klat),jk)+zont(jl,jk)
        zq(bini_l(jl,klat),jk) = zq(bini_l(jl,klat),jk)+zonq(jl,jk)
      END DO
    END DO

    IF (lwiso) THEN
      DO jt = 1, nwiso      ! redundantly performed by all PEs in a row - water isotopes
        DO jl = 1, nplon
          DO jk = 1, plev
              zwisoq(bini_l(jl,klat),jk,jt) = zwisoq(bini_l(jl,klat),jk,jt)+zonwisoq(jl,jk,jt)
          END DO
        END DO
      END DO
    END IF

    DO jk = 1, plev
      DO jl = 1, plon
        jlg = ldc%glon(klat) + jl
        ptte(jl,jk) = zt(bini_l(jlg,klat),jk)/nptslat_l(klat)
        pqte(jl,jk) = zq(bini_l(jlg,klat),jk)/nptslat_l(klat)
      END DO
    END DO

    IF (lwiso) THEN
      DO jt = 1, nwiso
        DO jk = 1, plev
          DO jl = 1, plon
            jlg = ldc%glon(klat) + jl
            pwisoqte(jl,jk,jt) = zwisoq(bini_l(jlg,klat),jk,jt)/nptslat_l(klat)
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE pole_filter_wiso

!---wiso-code-end

!!============================================================================

#ifndef __uxp__
  FUNCTION minmod(a, b) RESULT(c)
    REAL(wp), INTENT(in) :: a, b
    REAL(wp) :: c
    
    c = 0.5_wp*(SIGN(1.0_wp,a) + SIGN(1.0_wp,b))*MIN(ABS(a),ABS(b))
    
  END FUNCTION minmod
  
  FUNCTION medan(a, b, c) RESULT(d)
    REAL(wp), INTENT(in) :: a, b, c
    REAL(wp)  :: d
    REAL(wp) :: z1, z2
  
    z1 = b-a
    z2 = c-a
    d = a + 0.5_wp*(SIGN(1.0_wp,z1) + SIGN(1.0_wp,z2))*MIN(ABS(z1),ABS(z2))
  
    ! original version requires 2 level nested inlining ;-)
    ! d = a + minmod(b-a,c-a)
  
  END FUNCTION medan
#endif 

END MODULE mo_spitfire

#ifdef __uxp__

  FUNCTION minmod(a, b) RESULT(c)

    USE mo_kind, ONLY: wp	

    REAL(wp), INTENT(in) :: a, b
    REAL(wp) :: c
    
    c = 0.5_wp*(SIGN(1.0_wp,a) + SIGN(1.0_wp,b))*MIN(ABS(a),ABS(b))
    
  END FUNCTION minmod
  
  FUNCTION medan(a, b, c) RESULT(d)

    USE mo_kind, ONLY: wp	

    REAL(wp), INTENT(in) :: a, b, c
    REAL(wp) :: d
    REAL(wp) :: z1, z2
  
    z1 = b-a
    z2 = c-a
    d = a + 0.5_wp*(SIGN(1.0_wp,z1) + SIGN(1.0_wp,z2))*MIN(ABS(z1),ABS(z2))
  
    ! original version requires 2 level nested inlining ;-)
    ! d = a + minmod(b-a,c-a)
  
  END FUNCTION medan

#endif
