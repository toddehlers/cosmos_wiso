MODULE mo_geoloc

  USE mo_kind,        ONLY: dp
  USE mo_radiation,   ONLY: io3, iaero

  !
  ! This module holds quantities related to the geolocation to be used
  ! by the physical parameterizations.
  !
  ! Author: 
  ! A. Rhodin,       MPI, June 2001, original source
  ! U. Schulzweida,  MPI, May  2002, blocking (nproma)
  ! A. Rhodin,       DWD, Aug  2002, lon/lat index arrays
  ! L. Kornblueh,    MPI, Oct  2003, added budw_2d and budw (for scan loop)      
  !                                  changed dp to wp 
  !
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  !
  ! Public entities
  !

  PRIVATE

  !
  ! module procedures
  !

  PUBLIC :: init_geoloc    ! allocate memory, set time independent parameters
  PUBLIC :: set_geoloc     ! set time dependent parameters
  PUBLIC :: cleanup_geoloc ! deallocate module variables

  !
  ! Parameters (2D)
  !

  REAL(dp) ,POINTER, PUBLIC :: gboxarea_2d (:,:) ! grid box surface area
  REAL(dp) ,POINTER, PUBLIC :: philat_2d   (:,:) ! latitude  (degree)
  REAL(dp) ,POINTER, PUBLIC :: philon_2d   (:,:) ! longitude (degree)
  REAL(dp) ,POINTER, PUBLIC :: sqcst_2d    (:,:) !     cos(Gaussian latitudes)
  REAL(dp) ,POINTER, PUBLIC :: twomu_2d    (:,:) ! 2 * sin(Gaussian latitudes)
  REAL(dp) ,POINTER, PUBLIC :: budw_2d     (:,:)
  REAL(dp) ,POINTER, PUBLIC :: cst_2d      (:,:)
  REAL(dp) ,POINTER, PUBLIC :: rcsth_2d    (:,:)
  REAL(dp) ,POINTER, PUBLIC :: coriol_2d   (:,:)
  REAL(dp) ,POINTER, PUBLIC :: racst_2d    (:,:)

  REAL(dp), POINTER, PUBLIC :: amu0_x(:,:)
  REAL(dp), POINTER, PUBLIC :: rdayl_x(:,:)
  REAL(dp), POINTER, PUBLIC :: amu0m_x(:,:)
  REAL(dp), POINTER, PUBLIC :: rdaylm_x(:,:)
  REAL(dp), POINTER, PUBLIC :: zaes_x(:,:)
  REAL(dp), POINTER, PUBLIC :: zael_x(:,:)
  REAL(dp), POINTER, PUBLIC :: zaeu_x(:,:)
  REAL(dp), POINTER, PUBLIC :: zaed_x(:,:)
  REAL(dp), POINTER, PUBLIC :: zozq_x(:,:)
  REAL(dp), POINTER, PUBLIC :: zozh_x(:,:)

  INTEGER, ALLOCATABLE, PUBLIC :: ilon (:,:) ! global longitude index array
  INTEGER, ALLOCATABLE, PUBLIC :: ilat (:,:) ! global latitude index array
  
  !
  ! Parameters (1D, used within latitude loop below scan1)
  !


  REAL(dp), POINTER, PUBLIC :: rz1u0(:,:,:)

  !------------------------------------------------------------------------------
CONTAINS
  !------------------------------------------------------------------------------
  SUBROUTINE init_geoloc
    !
    ! allocate memory, set time independent parameters
    ! to be called during startup
    !
    USE mo_memory_base,   ONLY: t_stream, new_stream, add_stream_element, &
                                default_stream_setting, AUTO
    USE mo_decomposition, ONLY: ldc => local_decomposition
    USE mo_transpose,     ONLY: reorder
    USE mo_gaussgrid,     ONLY: gridarea, philon, philat, gl_sqcst, gl_coriol, &
                                gl_racst, gl_cst, gl_rcsth, gl_twomu, gl_budw

    TYPE (t_stream) ,POINTER :: geoloc

    INTEGER  :: jlat, jlon, jg, nlon, nglat, nglon
    REAL(dp) :: tmp (ldc% nglon, ldc% nglat)
    REAL(dp) :: tmp2(ldc% nproma, ldc% ngpblks)
    !
    ! request a new memory buffer
    !
    CALL new_stream (geoloc, 'geoloc',                         &
         lpost=.FALSE., lrerun=.FALSE., linit=.FALSE.)
    !
    ! request 2D fields
    !
    CALL default_stream_setting (geoloc, lpost=.FALSE., code=AUTO)

    CALL add_stream_element (geoloc, 'gboxarea', gboxarea_2d, units='m^2')
    CALL add_stream_element (geoloc, 'philat',   philat_2d,   units='degree')
    CALL add_stream_element (geoloc, 'philon',   philon_2d,   units='degree')

    CALL add_stream_element (geoloc, 'sqcst',  sqcst_2d,  units='')
    CALL add_stream_element (geoloc, 'twomu',  twomu_2d,  units='')
    CALL add_stream_element (geoloc, 'budw',   budw_2d,   units='')
    CALL add_stream_element (geoloc, 'cst',    cst_2d,    units='')
    CALL add_stream_element (geoloc, 'rcsth',  rcsth_2d,  units='')
    CALL add_stream_element (geoloc, 'coriol', coriol_2d, units='')
    CALL add_stream_element (geoloc, 'racst',  racst_2d,  units='')

    CALL add_stream_element (geoloc, 'rz1u0',  rz1u0,    units='')

    CALL add_stream_element (geoloc, 'amu0',   amu0_x,   units='')
    CALL add_stream_element (geoloc, 'rdayl',  rdayl_x,  units='')
    CALL add_stream_element (geoloc, 'amu0m',  amu0m_x,  units='')
    CALL add_stream_element (geoloc, 'rdaylm', rdaylm_x, units='')

    IF ( iaero == 2 .OR. iaero == 3 .OR. iaero == 4 ) THEN
      CALL add_stream_element (geoloc, 'zaes', zaes_x, units='')
      CALL add_stream_element (geoloc, 'zael', zael_x, units='')
      CALL add_stream_element (geoloc, 'zaeu', zaeu_x, units='')
      CALL add_stream_element (geoloc, 'zaed', zaed_x, units='')
    END IF

    IF ( io3 == 2 ) THEN
      CALL add_stream_element (geoloc, 'zozq', zozq_x, units='')
      CALL add_stream_element (geoloc, 'zozh', zozh_x, units='')
    END IF

    ALLOCATE (ilon (ldc% nproma, ldc% ngpblks))
    ALLOCATE (ilat (ldc% nproma, ldc% ngpblks))

    !
    ! set time independent parameters
    !

    nlon  = ldc% nlon
    nglat = ldc% nglat
    nglon = ldc% nglon
    CALL reorder (gboxarea_2d ,SPREAD( gridarea (ldc%glat(:)),1 ,nglon))
    CALL reorder (philat_2d   ,SPREAD( philat   (ldc%glat(:)),1 ,nglon))
    CALL reorder (sqcst_2d    ,SPREAD( gl_sqcst (ldc%glat(:)),1 ,nglon))
    CALL reorder (twomu_2d    ,SPREAD( gl_twomu (ldc%glat(:)),1 ,nglon))
    CALL reorder (budw_2d     ,SPREAD( gl_budw  (ldc%glat(:)),1 ,nglon))
    CALL reorder (cst_2d      ,SPREAD( gl_cst   (ldc%glat(:)),1 ,nglon))
    CALL reorder (rcsth_2d    ,SPREAD( gl_rcsth (ldc%glat(:)),1 ,nglon))
    CALL reorder (coriol_2d   ,SPREAD( gl_coriol(ldc%glat(:)),1 ,nglon))
    CALL reorder (racst_2d    ,SPREAD( gl_racst (ldc%glat(:)),1 ,nglon))

    DO jlat = 1, ldc% nglat     ! local  latitude index N -> S
      jg   =    ldc% glat(jlat) ! global latitude index N -> S
      DO jlon = 1, ldc% nglon   ! local  longitude index
        tmp(jlon,jlat) = philon (jlon + ldc% glon(jlat))
      END DO
    END DO
    CALL reorder (philon_2d, tmp)

    tmp = SPREAD (ldc% glat(:) ,1 ,nglon)
    CALL reorder (tmp2 ,tmp)
    ilat = tmp2

    DO jlat = 1, ldc% nglat     ! local  latitude index N -> S
      jg   =    ldc% glat(jlat) ! global latitude index N -> S
      DO jlon = 1, ldc% nglon   ! local  longitude index
        tmp(jlon,jlat) = jlon + ldc% glon(jlat)
      END DO
    END DO
    CALL reorder (tmp2 ,tmp)
    ilon = tmp2

  END SUBROUTINE init_geoloc
  !------------------------------------------------------------------------------
  SUBROUTINE set_geoloc
    !
    ! set time dependent parameters
    ! to be called during the time stepping
    !  
  END SUBROUTINE set_geoloc
  !----------------------------------------------------------------------------
  SUBROUTINE cleanup_geoloc
    DEALLOCATE (ilon)
    DEALLOCATE (ilat)
  END SUBROUTINE cleanup_geoloc
  !----------------------------------------------------------------------------
END MODULE mo_geoloc
