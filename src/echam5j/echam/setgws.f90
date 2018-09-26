SUBROUTINE setgws

  !- Description:
  !
  !  Modify preset variables of module mo_gwspectrum which control
  !  the configuration of the Hines gravity waves parameterization 
  !
  !-  Method:
  !
  !  Read the gwsclt namelist and modify constants.
  !
  !  *setgws* is called from *initialize*.
  !
  !- Authors:
  !
  !   E. Manzini, MPI,  January 2002
  !   (re-write from M. Charron original)
  !

  USE mo_kind           ,ONLY: dp
  USE mo_doctor         ,ONLY: nout,nerr
  USE mo_mpi            ,ONLY: p_parallel,p_parallel_io,p_bcast,p_io
  USE mo_exception      ,ONLY: finish
  USE mo_gwspectrum     ,ONLY: lextro,lfront,lozpr,iheatcal,       &
                               rmscon,emiss_lev,kstar,m_min,       &
                               rms_front,front_thres,pcrit,pcons,  &
                               naz,slope,smco,gwsctl
  USE mo_control        ,ONLY: nn
  USE mo_namelist       ,ONLY: position_nml, nnml, POSITIONED

  IMPLICIT NONE
 
  ! Local variable
  INTEGER :: ierr

  ! Executable statements 

  ! Read gwsctl namelist to modify mo_gwspectrum
  ! ============================================

  IF (p_parallel_io) THEN
     CALL position_nml ('GWSCTL', status=ierr)
     SELECT CASE (ierr)
     CASE (POSITIONED)
       READ (nnml, gwsctl)
     END SELECT
  ENDIF
  IF (p_parallel) THEN
     CALL p_bcast (lextro, p_io)
     CALL p_bcast (lfront, p_io)
     CALL p_bcast (lozpr, p_io)
     CALL p_bcast (iheatcal, p_io)
     CALL p_bcast (rmscon, p_io)
     CALL p_bcast (emiss_lev, p_io)
     CALL p_bcast (kstar, p_io)
     CALL p_bcast (m_min, p_io)
     CALL p_bcast (rms_front, p_io)
     CALL p_bcast (front_thres, p_io)
     CALL p_bcast (pcrit, p_io)
     CALL p_bcast (pcons, p_io)
  ENDIF

  ! Write modified namelist
  ! =======================
  IF (.NOT. p_parallel) THEN
     WRITE(nout,*) 'Namelist gwsctl modified by setgws: '
     WRITE(nout,gwsctl)
  END IF


  ! Check lextro
  ! ============
  IF (.NOT.lextro) THEN
     IF (p_parallel_io) &
       WRITE(nout,*)    'lextro = .FALSE. --> no hines gw spectrum'
  ENDIF

  IF (lextro) THEN

     ! Initialization for Hines gravity wave parameterization
     ! =======================================================

     !  check that things are setup correctly and log error if not


     IF (lfront) THEN
        IF (p_parallel_io) &
           WRITE(nout,*)    'lfront = .TRUE.  --> gravity waves from fronts' 
        IF (nn==31)  THEN
           front_thres  = 0.10_dp
        ENDIF
        IF ( (nn /= 31) .AND. (nn /= 42) ) THEN
           IF (p_parallel_io) &
              WRITE(nerr,*) 'LFRONT=.true. not implemented for nn=', nn 
           CALL finish('setgws','Run terminated')
        ENDIF
     ENDIF

     IF (lozpr) THEN
        IF (p_parallel_io) &
           WRITE(nerr,*) 'LOZPR=.TRUE. not implemented for maecham5'
        CALL finish('setgws','Run terminated') 
     ENDIF

     SELECT CASE(iheatcal)
     CASE(0)
        IF (p_parallel_io) &
           WRITE(nout,*) 'standard maecham5: momentum flux deposition ONLY'
     CASE(1)
        IF (p_parallel_io) &
           WRITE(nout,*) ' iheatcal = 1 : momentum deposion, heating and diff coeff'
        IF ( ABS(slope-1._dp) > EPSILON(1._dp) ) THEN
           IF (p_parallel_io) &
              WRITE(nerr,*) 'Error: Heat and coeff for slope not 1 not implemented'
           CALL finish('setgws','Run terminated')
        ENDIF
     END SELECT

     IF ( naz /= 8 )  then
        IF (p_parallel_io) &
           WRITE(nerr,*) 'Error: naz not 8 not implemented' 
        CALL finish('setgws','Run terminated')
     endif

     IF (ABS(slope-1._dp) > EPSILON(1._dp)  ) then
        IF (p_parallel_io) &
           WRITE(nerr,*) 'Error:  slope not 1. not implemented' 
        CALL finish('setgws','Run terminated')
     endif

     IF ( (m_min > EPSILON(1._dp)) .AND. (ABS(slope-1._dp) > EPSILON(1._dp)) ) then
        IF (p_parallel_io) &
           WRITE(nerr,*) 'Error:  m_min not zero for slope not 1. not implemented' 
        CALL finish('setgws','Run terminated')
     endif

     IF ( smco < 1._dp ) then
        IF (p_parallel_io) &
           WRITE(nerr,*) 'Error:  smco less than 1. not implemented' 
        CALL finish('setgws','Run terminated')
     endif
     !

     !------------------------------------------------------------
     ! set up for Hammonia:
     !  rmscon       = 1.0
     !  m_minimum    = 3.141593e-4     ! = 2*pi/(20 km)
     !  iheatcal    = 1
     !------------------------------------------------------------

  ENDIF
  !
  RETURN
END SUBROUTINE setgws
