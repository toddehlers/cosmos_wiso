SUBROUTINE inipost

  ! Description:
  !
  ! Preset constants in mo_control.
  !
  ! Method:
  !
  ! Preset module *mo_post*.
  !
  ! *inipost* is called from subroutine *initialize*
  !
  ! all variables of module *mo_post* are defined
  !
  ! Authors:
  !
  ! E. Kirk, MI, March 1989, original source
  ! U. Schlese, MPI, January 1990, g3x-fields added
  ! U. Schlese, DKRZ, May 1992, statistics of horizontal diffusion 
  ! U. Schlese, DKRZ, April 1993 g4x-fields added
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! L. Kornblueh, MPI, April 1998, added NWP forecast mode
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! M. Esch, MPI, June 1999, ECHAM5-modifications
  ! U. Schlese, DKRZ, March 2000, obsolete/new codes
  !                               8 Bit packing increased to 16 Bit
  ! U. Schlese, DKRZ, new codes for HD discharge and Glacier calving,
  !                   packing of water balance codes increased to 24 bits
  ! I. Kirchner, MPI, November 2000, date/time control
  ! A. Rhodin, MPI, June 2001, changes for output stream interface
  ! for more details see file AUTHORS
  !

  USE mo_doctor,        ONLY: nout
  USE mo_post,          ONLY: lppspe, lppd, lppvo, lppt, lppp, lppq    &
                            , lppxl, lppxi
  USE mo_mpi,           ONLY: p_io, p_parallel, p_parallel_io, p_bcast
  USE mo_namelist,      ONLY: position_nml, nnml, POSITIONED
  USE mo_memory_base,   ONLY: set_stream_element_info
  USE mo_memory_sp,     ONLY: sp  !spectral                output stream buffer
  USE mo_memory_gl,     ONLY: gl  !gridpoint (transported) output stream buffer
!---wiso-code
  USE mo_wiso,          ONLY: lwiso, nwiso, wisoq_names, wisoxl_names, wisoxi_names
!---wiso-code-end

  IMPLICIT NONE

  !  Local scalars:
  INTEGER       :: ierr

!---wiso-code
  INTEGER       :: jt
!---wiso-code-end

  INCLUDE 'postctl.inc'

  !  Executable statements 

  !-- 1. Define variables

  lppspe = .TRUE.
  lppt = .TRUE.
  lppd = .TRUE.
  lppvo = .TRUE.
  lppq = .TRUE.
  lppp = .TRUE.
  lppxl = .TRUE.
  lppxi = .TRUE.

  !
  !  variables for fractional surface coverage (91 -126)
  !

  !-- 2. Read namelist postctl

  IF (p_parallel_io) THEN
    CALL position_nml ('POSTCTL',status=ierr)
    IF (ierr == POSITIONED) READ (nnml, postctl)
  ENDIF
  IF (p_parallel) THEN
    CALL p_bcast (lppspe, p_io)
    CALL p_bcast (lppd, p_io)
    CALL p_bcast (lppvo, p_io)
    CALL p_bcast (lppt, p_io)
    CALL p_bcast (lppp, p_io)
    CALL p_bcast (lppq, p_io)
    CALL p_bcast (lppxl, p_io)   
    CALL p_bcast (lppxi, p_io)   
  ENDIF

  IF ( .NOT. lppspe) THEN
    lppd     = .FALSE.
    lppvo    = .FALSE.
    lppt     = .FALSE.
    lppp     = .FALSE.
    lppq     = .FALSE.
    lppxl    = .FALSE.
    lppxi    = .FALSE.
  END IF

  ! set print flags in spectral output stream buffer

  CALL set_stream_element_info (sp ,'sd'   ,lpost=lppd)    ! divergence
  CALL set_stream_element_info (sp ,'svo'  ,lpost=lppvo)   ! vorticity
  CALL set_stream_element_info (sp ,'st'   ,lpost=lppt)    ! temperature
  CALL set_stream_element_info (sp ,'sp'   ,lpost=lppp)    ! surface pressure

  ! set print flags in gridpoint (spitfire transport) output stream buffer

  CALL set_stream_element_info (gl ,'q'    ,lpost=lppq)    ! humidity
  CALL set_stream_element_info (gl ,'xl'   ,lpost=lppxl)   ! cloud water
  CALL set_stream_element_info (gl ,'xi'   ,lpost=lppxi)   ! cloud ice

!---wiso-code

  ! set print flags of water isotopes in gridpoint (spitfire transport) output stream buffer to .FALSE.,
  ! as they are written out to extra water isotope stream (in mo_memory_wiso)
  
  IF (lwiso) THEN
    DO jt=1,nwiso
      CALL set_stream_element_info (gl ,wisoq_names(jt)    ,lpost=.FALSE.)   ! humidity
      CALL set_stream_element_info (gl ,wisoxl_names(jt)   ,lpost=.FALSE.)   ! cloud water
      CALL set_stream_element_info (gl ,wisoxi_names(jt)   ,lpost=.FALSE.)   ! cloud ice
    END DO
  END IF

!---wiso-code-end


  !-- 3. Write namelist postctl

  IF (.NOT. p_parallel) THEN
    WRITE (nout,postctl)
  END IF

END SUBROUTINE inipost
