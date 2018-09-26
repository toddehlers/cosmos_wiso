MODULE mo_memory_wiso
  !
  ! This file contains the routines to allocate memory for water isotope-related variables
  ! and to define the output stream 'wiso' for water isotopes
  !
  ! Authors:
  ! -------
  ! M. Werner, AWI, 2012

  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_linked_list,   ONLY: t_stream
  USE mo_memory_base,   ONLY: delete_stream, add => add_stream_element, &
                              default_stream_setting,                   &
                              ABOVESUR2, ABOVESUR10, BELOWSUR, HYBRID_H,&
                              HYBRID
  USE mo_netCDF,        ONLY: max_dim_name
  USE mo_wiso,          ONLY: nwiso, lwiso,                                 &
                              wisoq_names,   wisoxl_names,   wisoxi_names,  &
                              wisoqm1_names, wisoxlm1_names, wisoxim1_names
  
  IMPLICIT NONE

  REAL(dp), ALLOCATABLE, PUBLIC :: wisoqte   (:,:,:,:)
  REAL(dp), ALLOCATABLE, PUBLIC :: wisoxlte  (:,:,:,:)
  REAL(dp), ALLOCATABLE, PUBLIC :: wisoxite  (:,:,:,:)

  LOGICAL, SAVE, PRIVATE :: lnot_used   = .TRUE.

  PRIVATE

  !--- Service routines ----------------------------------------------------------------

  PUBLIC :: construct_wiso            ! construct the wiso table
  PUBLIC :: destruct_wiso             ! destruct  the wiso_g1b table
  PUBLIC :: m_bufscan_wiso
  PUBLIC :: cleanup_scanbuffer_wiso

  PUBLIC :: wiso                      ! the wiso output stream


  !--- Declarations for stream wiso ----------------------------------------------------

  TYPE (t_stream), POINTER :: wiso

  !--- Declaration of predefined fields within this table ----------------------------- 

  ! water isotope variables analog to stream g1b

  REAL(dp), POINTER, PUBLIC :: wisoqf(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxlf(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxif(:,:,:,:)


  ! water isotope variables analog to stream g3b

  REAL(dp), POINTER, PUBLIC :: wisows(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisowl(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisosn(:,:,:)

  REAL(dp), POINTER, PUBLIC :: wisoaprl(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoaprc(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoaprs(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoevap(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisorunoff(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisogld(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisosnmel(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoruntoc(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoapmegl(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoqvi(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxlvi(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxivi(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisodrain(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisosnacl(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisorogl(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxtec(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisosnc(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoapmeb(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoapmebco(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisorain(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoqtnew(:,:,:)
  !
  !  variables for fractional surface coverage
  !
  REAL(dp), POINTER, PUBLIC :: wisoevapiac(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoevapwac(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoevaplac(:,:,:)
  !
  !  variables for ocean coupling only
  !
  REAL(dp), POINTER, PUBLIC :: wisoawfre(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoafre_residual(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoaifre(:,:,:)
  !
  !  variables for coupling with HD-model and calving model only
  !
  REAL(dp), POINTER, PUBLIC :: wisoaros(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoadrain(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisodisch(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoapmecal(:,:,:)
  !
  !  variable for water isotopes in ocean surface water
  !
  REAL(dp), POINTER, PUBLIC :: wisosw_d(:,:,:)
  !
  !  variable for water isotopes in snow on glaciers
  !
  REAL(dp), POINTER, PUBLIC :: snglac(:,:)       ! new ECHAM5 variable: snow depth on glaciers
  REAL(dp), POINTER, PUBLIC :: wisosnglac(:,:,:) ! (required to simulate delta of snow on glaciers used for sublimation in vdiff)


  ! water isotope variables analog to stream g3a

  REAL(dp), POINTER, PUBLIC :: wisosnm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisowlm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisowsm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoaprlm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoaprcm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoaprsm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoevapm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisorunoffm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisosnmelm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoruntocm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoapmeglm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoqvim(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxlvim(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxivim(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisodrainm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisosnaclm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoroglm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisosw_dm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: snglacm(:,:)
  REAL(dp), POINTER, PUBLIC :: wisosnglacm(:,:,:)


  ! water isotope variables analog to stream gl

  REAL(dp), POINTER, PUBLIC :: wisoq  (:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxl (:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxi (:,:,:,:)

  REAL(dp), POINTER :: pwisoq  (:,:,:,:,:)
  REAL(dp), POINTER :: pwisoxl (:,:,:,:,:)
  REAL(dp), POINTER :: pwisoxi (:,:,:,:,:)
  
  ! water isotope variables analog to stream g1a

  REAL(dp), POINTER, PUBLIC :: wisoqm1  (:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxlm1 (:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxim1 (:,:,:,:)

  REAL(dp), POINTER :: pwisoqm1 (:,:,:,:,:)
  REAL(dp), POINTER :: pwisoxlm1 (:,:,:,:,:)
  REAL(dp), POINTER :: pwisoxim1 (:,:,:,:,:)

!----------------------------------------------------------------------------------------

CONTAINS

!---------------------------------------------------------------------------------------

  SUBROUTINE construct_wiso (lnlon, lnlev, lnwiso, lngl, &
                             nlon,  nlev,   nwiso, ngl)
    
    USE mo_control,   ONLY: lcouple, lhd

    INTEGER, INTENT (in) :: lnlon, lnlev, lnwiso, lngl
    INTEGER, INTENT (in) ::  nlon,  nlev,  nwiso, ngl

    INTEGER                      :: dim1(4), dim1p(4)
    INTEGER                      :: dim2(3), dim2p(3)
    INTEGER                      :: dim3(2), dim3p(2)
    CHARACTER (len=max_dim_name) :: dim1n(4),dim2n(3),dim3n(2)
    INTEGER                      :: i, zi
    REAL(dp), POINTER            :: p3(:,:,:), p4(:,:,:,:)


    CALL default_stream_setting (wiso              &
                                ,lrerun=.TRUE.     &
                                ,lpost=.TRUE.      &
                                ,table=128 ,bits=16)


    ! construct the water isotope equivalent of the g1b table

    dim1p = (/ lnlon, lnlev, lnwiso, lngl /)
    dim1  = (/  nlon,  nlev,  nwiso,  ngl /)
    dim1n = (/  "lon   ","lev   ","nwiso ","lat   "/)

    dim2p = (/  lnlon,  lnlev, lngl  /)
    dim2  = (/   nlon,   nlev,  ngl  /)
    dim2n = (/   "lon ","lev ","lat "/)

    dim3p = (/  lnlon,  lngl  /)
    dim3  = (/   nlon,   ngl  /)
    dim3n = (/   "lon ","lat "/)

    CALL add (wiso, 'wisoqf', wisoqf, dim1p, dim1, lrerun=.FALSE., lpost=.FALSE.)
    CALL add (wiso, 'wisoxlf',wisoxlf,dim1p, dim1, lrerun=.FALSE., lpost=.FALSE.)
    CALL add (wiso, 'wisoxif',wisoxif,dim1p, dim1, lrerun=.FALSE., lpost=.FALSE.)


    ! construct the water isotope equivalent of the g3b table

    CALL add (wiso,'wisoqtnew',  wisoqtnew   ,lpost=.FALSE.,contnorest=.true.,klev=nwiso)

    CALL add (wiso,'wisows',     wisows,     code=51,             longname='soil wetness - water isotopes',          &
                                                                                   units='m',        klev=nwiso)
    CALL add (wiso,'wisosn',     wisosn,     code=52,             longname='snow depth - water isotopes',            &
                                                                                   units='m',        klev=nwiso)
    CALL add (wiso,'wisoaprl',   wisoaprl,   code=53,laccu=.TRUE.,longname='l.scale precipitation - water isotopes', &
                                                                          bits=24, units='kg/m**2s', klev=nwiso)
    CALL add (wiso,'wisoaprc',   wisoaprc,   code=54,laccu=.TRUE.,longname='conv. precipitation - water isotopes',   &
                                                                          bits=24, units='kg/m**2s', klev=nwiso)
    CALL add (wiso,'wisoaprs',   wisoaprs,   code=55,laccu=.TRUE.,longname='snow fall - water isotopes',             &
                                                                          bits=24, units='kg/m**2s', klev=nwiso)
    CALL add (wiso,'wisoxivi',   wisoxivi,   code=56,laccu=.TRUE.,longname='vert. int. cloud ice - water isotopes',  &
                                                                                   units='kg/m**2',  klev=nwiso)
    CALL add (wiso,'wisorunoff', wisorunoff, code=57,laccu=.TRUE.,longname='runoff and drainage - water isotopes',   &
                                                                          bits=24, units='kg/m**2s', klev=nwiso)
    CALL add (wiso,'wisodrain',  wisodrain,  code=58,laccu=.TRUE.,longname='drainage - water isotopes',              &
                                                                          bits=24, units='kg/m**2s', klev=nwiso)
    CALL add (wiso,'wisoevap',   wisoevap,   code=59,laccu=.TRUE.,longname='evaporation - water isotopes',           &
                                                                          bits=24, units='kg/m**2s', klev=nwiso)
    CALL add (wiso,'wisoevapiac',wisoevapiac,code=60,laccu=.TRUE.,longname='evaporation over ice - water isotopes',  &
                                                                                   units='kg/m**2s', klev=nwiso)
    CALL add (wiso,'wisoevapwac',wisoevapwac,code=61,laccu=.TRUE.,longname='evaporation over water - water isotopes',&
                                                                                   units='kg/m**2s', klev=nwiso)
    CALL add (wiso,'wisoevaplac',wisoevaplac,code=62,laccu=.TRUE.,longname='evaporation over land - water isotopes', &
                                                                                   units='kg/m**2s', klev=nwiso)
    CALL add (wiso,'wisowl',     wisowl,     code=63,             longname='skin reservoir content - water isotopes',&
                                                                                   units='m',        klev=nwiso)
    CALL add (wiso,'wisogld',    wisogld,    code=64,             longname='glacier depth - water isotopes',         & 
                                                                                   units='m',        klev=nwiso)
    CALL add (wiso,'wisorogl',   wisorogl,   code=65,laccu=.TRUE.,longname='glacier runoff - water isotopes',        & 
                                                                                   units='kg/m**2s', klev=nwiso)
    CALL add (wiso,'wisosnmel',  wisosnmel,  code=66,laccu=.TRUE.,longname='snow melt - water isotopes',             & 
                                                                                   units='kg/m**2s', klev=nwiso)
    CALL add (wiso,'wisoruntoc', wisoruntoc, code=74,laccu=.TRUE.,lpost=.FALSE.,longname='surf.runoff into ocean - water isotopes',&
                                                                          bits=24 ,units='kg/m**2s', klev=nwiso)
    CALL add (wiso,'wisoapmegl', wisoapmegl, code=67,laccu=.TRUE.,longname='P-E over land ice - water isotopes',     & 
                                                                          bits=24, units='kg/m**2s', klev=nwiso)
    CALL add (wiso,'wisosnacl',  wisosnacl,  code=68,laccu=.TRUE.,longname='snow accu. over land - water isotopes',  & 
                                                                                   units='kg/m**2s', klev=nwiso)
    CALL add (wiso,'wisoqvi',    wisoqvi,    code=69,laccu=.TRUE.,longname='vert. int. water vapor - water isotopes',&
                                                                                   units='kg/m**2',  klev=nwiso)
    CALL add (wiso,'wisoxlvi',   wisoxlvi,   code=70,laccu=.TRUE.,longname='vert. int. cloud water - water isotopes',&
                                                                                   units='kg/m**2',  klev=nwiso)
    CALL add (wiso,'wisosnc',    wisosnc,    code=71,             longname='snow depth at canopy - water isotopes',  & 
                                                                                   units='m',        klev=nwiso)
    CALL add (wiso,'wisosw_d',   wisosw_d,   code=72,             longname='delta of ocean surface - water isotopes',& 
                                                                                   units='permil',   klev=nwiso)
    CALL add (wiso,'wisoapmeb',  wisoapmeb,  code=73,laccu=.TRUE.,longname='vert.integr.tendencies of water - water isotopes', &
                                                                          bits=24, units='kg/m**2s',klev=nwiso)
    CALL add (wiso,'wisosnglac', wisosnglac, code=76,             longname='snow depth on glaciers - water isotopes',&
                                                                                   units='m',        klev=nwiso)
    CALL add (wiso,'snglac',     snglac, ldims=dim3p, gdims=dim3, dimnames=dim3n, code=75, longname='snow depth on glaciers',  &
                                                                                   units='m')

    ! add fields not written to the output stream
    !
    CALL add (wiso,'wisoxtec',wisoxtec,lpost=.FALSE. ,ktrac=nwiso)

    !  variables for ocean coupling only
    !
    CALL add (wiso,'wisoapmebco',  wisoapmebco,   lpost=.FALSE., contnorest=.true., klev=nwiso)
    CALL add (wiso,'wisorain',     wisorain,      lpost=.FALSE., contnorest=.true., klev=nwiso)
    !
    IF (lcouple) THEN
      CALL add (wiso,'wisoawfre',  wisoawfre,     lpost=.FALSE., klev=nwiso)
      CALL add (wiso,'wisoafre_residual',wisoafre_residual,lpost=.FALSE.,contnorest=.true.,klev=nwiso)
      CALL add (wiso,'wisoaifre',  wisoaifre,     lpost=.FALSE., klev=nwiso)
    END IF

    !  variables for coupling with HD-model only
    !
    IF (lhd) THEN
      CALL add (wiso,'wisoaros',     wisoaros,    lpost=.FALSE., klev=nwiso)
      CALL add (wiso,'wisoadrain',   wisoadrain,  lpost=.FALSE., klev=nwiso)
      CALL add (wiso,'wisodisch',    wisodisch,   lpost=.FALSE., klev=nwiso)
      CALL add (wiso,'wisoapmecal',  wisoapmecal, lpost=.FALSE., klev=nwiso)
    END IF

    ! construct the water isotope equivalent of the g3a table

    wisosnm       => wisosn
    wisowlm       => wisowl
    wisowsm       => wisows
    wisoaprlm     => wisoaprl
    wisoaprcm     => wisoaprc
    wisoaprsm     => wisoaprs
    wisoevapm     => wisoevap
    wisorunoffm   => wisorunoff
    wisosnmelm    => wisosnmel
    wisoruntocm   => wisoruntoc
    wisoapmeglm   => wisoapmegl
    wisoqvim      => wisoqvi
    wisoxlvim     => wisoxlvi
    wisoxivim     => wisoxivi
    wisodrainm    => wisodrain
    wisosnaclm    => wisosnacl
    wisoroglm     => wisorogl
    wisosw_dm     => wisosw_d
    wisosnglacm   => wisosnglac
    snglacm       => snglac


    ! construct the water isotope equivalent of the gl table

    !
    ! Allocate three 5d-array (with dummy index 5) for water isotopes.
    ! This array is referenced by the 4-d arrays 'wisoq', 'wisoxl', 'wsioxi',
    ! and by the 3-d arrays of individual water isotope tracers.
    !
    ALLOCATE (pwisoq  (lnlon, lnlev, nwiso, lngl, 1))
    ALLOCATE (pwisoxl (lnlon, lnlev, nwiso, lngl, 1))
    ALLOCATE (pwisoxi (lnlon, lnlev, nwiso, lngl, 1))
    !
    ! Set meta information on water isotope arrays.
    ! Set restart flag, obtain reference to memory info entry.
    !
    CALL default_stream_setting (wiso              &
                                ,lrerun=.TRUE.     &
                                ,lpost=.FALSE.)
    !
    p4 => pwisoq(:,:,:,:,1)
    CALL add (wiso, 'wisoq', wisoq, ldims=dim1p, gdims=dim1,dimnames=dim1n, &
                             longname='specific humidity - water isotopes', p4=p4,table=128)
    !
    p4 => pwisoxl(:,:,:,:,1)
    CALL add (wiso, 'wisoxl', wisoxl, ldims=dim1p, gdims=dim1,dimnames=dim1n,               &
                             longname='cloud water - water isotopes',       p4=p4,table=128)
    !
    p4 => pwisoxi(:,:,:,:,1)
    CALL add (wiso, 'wisoxi', wisoxi, ldims=dim1p, gdims=dim1,dimnames=dim1n,               &
                             longname='cloud ice - water isotopes',         p4=p4,table=128)
    !
    ! provide additional meta-information for individual tracers.
    !
    zi = 240     ! Water isotope tracer fields start at code #241
    !
    DO i = 1, nwiso
      p4 => pwisoq(:,:,i,:,:)
      zi = zi+1
      CALL add (wiso, wisoq_names(i), p3, ldims=dim2p, gdims=dim2, dimnames=dim2n,           &
                               longname='specific humidity - water isotopes',units='kg/kg', &
                               lpost=.TRUE., code = zi,p4=p4,table=128)
      !
      p4 => pwisoxl(:,:,i,:,:)
      zi = zi+1
      CALL add (wiso, wisoxl_names(i), p3, ldims=dim2p, gdims=dim2, dimnames=dim2n,          &
                               longname='cloud water - water isotopes',      units='kg/kg', &
                               lpost=.TRUE., code = zi,p4=p4,table=128)
      !
      p4 => pwisoxi(:,:,i,:,:)
      zi = zi+1
      CALL add (wiso, wisoxi_names(i), p3, ldims=dim2p, gdims=dim2, dimnames=dim2n,          &
                               longname='cloud ice - water isotopes',        units='kg/kg', &
                               lpost=.TRUE., code = zi,p4=p4,table=128)
    END DO


    ! construct the water isotope equivalent of the gla table

    !
    ! Allocate three 5d-arrays (with dummy index 5) for water isotope tracers.
    ! These arrays are referenced by the 4-d arrays 'WISOQM1','WISOXLM1','WISOXIM1'
    ! and by the 3-d arrays of individual tracers.
    !
    ALLOCATE (pwisoqm1 (lnlon, lnlev, nwiso, lngl, 1))
    ALLOCATE (pwisoxlm1(lnlon, lnlev, nwiso, lngl, 1))
    ALLOCATE (pwisoxim1(lnlon, lnlev, nwiso, lngl, 1))
    !
    ! Set meta information on water isotope arrays.
    ! Set restart flag, obtain reference to memory info entry.
    !
    CALL default_stream_setting (wiso              &
                                ,lrerun=.TRUE.     &
                                ,lpost=.FALSE.)
    !
    p4 => pwisoqm1(:,:,:,:,1)
    CALL add (wiso, 'wisoqm1', wisoqm1, dim1p, dim1, &
           dimnames   = dim1n,                               &
           lrerun     = .TRUE.,                              &
           contnorest = .TRUE.,                              &
           p4         = p4)
    !
    ! provide additional meta-information for individual WISOQM1 tracers.
    !
    DO i = 1, nwiso
      p4 => pwisoqm1(:,:,i,:,:)
      CALL add (wiso, wisoqm1_names(i),       &
                                   p3, dim2p, dim2,          &
                      dimnames   = dim2n,                    &
                      units      = 'kg/kg',                  &
                      lrerun     = .TRUE.,                   &
                      contnorest = .TRUE.,                   &
                      lpost      = .FALSE.,                  &
                      p4         = p4)
    END DO
    p4 => pwisoxlm1(:,:,:,:,1)
    CALL add (wiso, 'wisoxlm1', wisoxlm1, dim1p, dim1, &
           dimnames   = dim1n,                               &
           lrerun     = .TRUE.,                              &
           contnorest = .TRUE.,                              &
           p4         = p4)
    !
    ! provide additional meta-information for individual WISOXLM1 tracers.
    !
    DO i = 1, nwiso
      p4 => pwisoxlm1(:,:,i,:,:)
      CALL add (wiso, wisoxlm1_names(i),      &
                                   p3, dim2p, dim2,          &
                      dimnames   = dim2n,                    &
                      units      = 'kg/kg',                  &
                      lrerun     = .TRUE.,                   &
                      contnorest = .TRUE.,                   &
                      lpost      = .FALSE.,                  &
                      p4         = p4)
    END DO
    p4 => pwisoxim1(:,:,:,:,1)
    CALL add (wiso, 'wisoxim1', wisoxim1, dim1p, dim1, &
           dimnames   = dim1n,                               &
           lrerun     = .TRUE.,                              &
           contnorest = .TRUE.,                              &
           p4         = p4)
    !
    ! provide additional meta-information for individual WISOXIM1 tracers.
    !
    DO i = 1, nwiso
      p4 => pwisoxim1(:,:,i,:,:)
      CALL add (wiso, wisoxim1_names(i),      &
                                   p3, dim2p, dim2,          &
                      dimnames   = dim2n,                    &
                      units      = 'kg/kg',                  &
                      lrerun     = .TRUE.,                   &
                      contnorest = .TRUE.,                   &
                      lpost      = .FALSE.,                  &
                      p4         = p4)
    END DO



  END SUBROUTINE construct_wiso

!---------------------------------------------------------------------------------------

  SUBROUTINE destruct_wiso

    CALL delete_stream (wiso)

    IF(ASSOCIATED (pwisoq))    DEALLOCATE (pwisoq)
    IF(ASSOCIATED (pwisoxl))   DEALLOCATE (pwisoxl)
    IF(ASSOCIATED (pwisoxi))   DEALLOCATE (pwisoxi)

    IF(ASSOCIATED (pwisoqm1))  DEALLOCATE (pwisoqm1)
    IF(ASSOCIATED (pwisoxlm1)) DEALLOCATE (pwisoxlm1)
    IF(ASSOCIATED (pwisoxim1)) DEALLOCATE (pwisoxim1)

  END SUBROUTINE destruct_wiso

!---------------------------------------------------------------------------------------

  SUBROUTINE m_bufscan_wiso

    USE mo_wiso,          ONLY: nwiso

    INTEGER :: ngpblks, nlev, nproma, ngl

    IF (lnot_used) THEN

       ngl     = ldc% nglat
       ngpblks = ldc% ngpblks
       nlev    = ldc% nlev
       nproma  = ldc% nproma

       ALLOCATE (wisoqte  (nproma,nlev,nwiso,ngpblks)) ;wisoqte  = 0.0_dp
       ALLOCATE (wisoxlte (nproma,nlev,nwiso,ngpblks)) ;wisoxlte = 0.0_dp
       ALLOCATE (wisoxite (nproma,nlev,nwiso,ngpblks)) ;wisoxite = 0.0_dp

       lnot_used = .FALSE.

    ENDIF

  END SUBROUTINE m_bufscan_wiso

!---------------------------------------------------------------------------------------

  SUBROUTINE cleanup_scanbuffer_wiso
    !------------------------------------
    ! deallocate variables in this module
    !------------------------------------

    IF (.NOT. lnot_used) THEN

       DEALLOCATE (wisoqte  )
       DEALLOCATE (wisoxlte )
       DEALLOCATE (wisoxite )

       lnot_used   = .TRUE.

    END IF

  END SUBROUTINE cleanup_scanbuffer_wiso

!---------------------------------------------------------------------------------------

END MODULE mo_memory_wiso
