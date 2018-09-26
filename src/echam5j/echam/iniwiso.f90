SUBROUTINE iniwiso

  ! Description:
  !
  ! Initialize water isotope fields
  !
  ! Method:
  !
  ! The water isotope fields are set to analogue values
  ! as the normal water values 
  ! (as defined for a cold model start in *ioinitial* 
  !  or read in for a model re-start in *iorestart*)
  !
  ! Authors:
  !
  ! M. Werner, AWI, 2014
  !
  !

  USE mo_kind,          ONLY: dp
  USE mo_mpi,           ONLY: p_io, p_pe
  USE mo_doctor,        ONLY: nout

  USE mo_time_control,  ONLY: lresume


  USE mo_memory_gl,     ONLY: q, xl, xi
  USE mo_memory_g1a,    ONLY: qm1, xlm1, xim1

  USE mo_memory_g3b,         ONLY: sn, ws, wsmx
  USE mo_decomposition,      ONLY: ldc => local_decomposition
  USE mo_clim,               ONLY: tslclim
  USE mo_constants,          ONLY: tmelt

  USE mo_wiso,               ONLY: lwiso, lwiso_rerun,                 &
                                   nwiso, nwisotyp, twisoatm,          &
                                   twisosur1, twisosur2, tnat, snglacmx
  USE mo_memory_wiso,        ONLY: wisoq, wisoxl, wisoxi,              &
                                   wisoqm1, wisoxlm1, wisoxim1,        &
                                   wisosn, wisows, snglac, wisosnglac

  IMPLICIT NONE

  REAL(dp) :: tsurface(ldc%nproma,ldc%ngpblks),zwisosurf(ldc%nproma,nwiso,ldc%ngpblks)
  INTEGER  :: jt, jmonth

  ! Executable Statements

  ! Set atmospheric isotope values for vapour, liquid and solid water

  ! set isotopes values for time step t:
  DO jt=1,nwiso      
    wisoq(:,:,jt,:)  = q(:,:,:)*tnat(jt)*(1._dp+twisoatm(jt))        ! vapour isotopes = SMOW * q * (1 + initial atm. delta value)
    wisoxl(:,:,jt,:) = xl(:,:,:)*tnat(jt)                            ! liquid water isotopes = SMOW * xl
    wisoxi(:,:,jt,:) = xi(:,:,:)*tnat(jt)                            ! solid water isotopes = SMOW * xl
  END DO

  IF (lresume .AND. (.NOT. lwiso_rerun)) THEN
    ! set isotopes values for time step t-1:
    DO jt=1,nwiso      
      wisoqm1(:,:,jt,:)  = qm1(:,:,:)*tnat(jt)*(1._dp+twisoatm(jt))  ! vapour isotopes = SMOW * q * (1 + initial atm. delta value)
      wisoxlm1(:,:,jt,:) = xlm1(:,:,:)*tnat(jt)                      ! liquid water isotopes = SMOW * xl
      wisoxim1(:,:,jt,:) = xim1(:,:,:)*tnat(jt)                      ! solid water isotopes = SMOW * xl
    END DO  
  END IF
    
  ! Initialize water isotope land surface fields

  ! set snow depth on glaciers to constant value *snglacmx*

  snglac(:,:) = snglacmx
     
 
  ! calculate annual mean temperature from climate land surface temperatures

  tsurface(:,:) = 0._dp
  DO jmonth = 1, 12
    tsurface(:,:) = tsurface(:,:) + tslclim(:,:,jmonth) 
  END DO
  tsurface(:,:) = tsurface(:,:)/12._dp

  ! calculate temperature-dependent initial isotope values of land water reservoirs (soil, snow, skin layer)
  DO jt=1,nwiso
    ! zwisosurf = inital_deviation_1 from SMOW at surface * surf.temp. - initial_deviation_2
    zwisosurf(:,jt,:)=twisosur1(jt)*(tsurface(:,:) - tmelt) - twisosur2(jt)
    IF (nwisotyp(jt).eq.2) THEN  ! if tracer = 18o: zwisosurf = min (-4./1000., zwisosurf)
      zwisosurf(:,jt,:)=min( -4._dp/1000._dp,zwisosurf(:,jt,:))
    ELSE IF (nwisotyp(jt).eq.3) THEN  ! if tracer = hdo: zwisosurf = min (-22./1000., zwisosurf)
      zwisosurf(:,jt,:)=min(-22._dp/1000._dp,zwisosurf(:,jt,:))
    ELSE IF (nwisotyp(jt).eq.4) THEN  ! if tracer = o17: (zwisosurf_o17+1) = (zwisosurf_o18+1)**0.528
      zwisosurf(:,jt,:)=(zwisosurf(:,2,:)+1._dp)**0.528_dp - 1._dp
    ENDIF
  END DO

  DO jt=1,nwiso
    ! set tracer reservoirs = nat.isotope-ratio * corresponding wetness reservoir * (1.+zwisosurf)
    wisosn(:,jt,:)=tnat(jt)*sn(:,:)*(1._dp+zwisosurf(:,jt,:))
    wisows(:,jt,:)=tnat(jt)*ws(:,:)*(1._dp+zwisosurf(:,jt,:))
    wisosnglac(:,jt,:)=tnat(jt)*snglac(:,:)*(1._dp+zwisosurf(:,jt,:))
    
    wisows(:,jt,:) = MIN(wisows(:,jt,:),wsmx(:,:)) ! maximum ws value for any tracer: wsmx (isotope-independent)
  END DO

  ! Prepare further grid point surface fields

  IF (lresume .AND. (.NOT. lwiso_rerun)) THEN
    CALL resume_g3_wiso()
    IF (p_pe == p_io) &
       WRITE (nout,'(a,/,/)') 'Water isotopes resumed from default ECHAM fields.'
  ELSE
    CALL init_g3_wiso()
    IF (p_pe == p_io) &
       WRITE (nout,'(a,/,/)') 'Water isotopes set to initial values.'
  ENDIF
 
  RETURN


CONTAINS

  SUBROUTINE init_g3_wiso ()

    !
    ! init_g3_wiso - initialize some further water isotope surface fields
    !
    ! set *wiso* variables to zero analogue to the default ECHAM fields

    USE mo_kind,          ONLY: dp
    USE mo_control,       ONLY: lhd
    Use mo_memory_wiso

    IMPLICIT NONE

    !  Executable Statements

    ! Initialize *wiso* variables not read

    wisoaprlm(:,:,:)  = 0.0_dp
    wisoaprcm(:,:,:)  = 0.0_dp
    wisoaprsm(:,:,:)  = 0.0_dp
    wisoevapm(:,:,:)  = 0.0_dp
    wisorunoffm(:,:,:)= 0.0_dp
    wisosnmelm(:,:,:) = 0.0_dp
    wisoruntocm(:,:,:)= 0.0_dp
    wisoapmeglm(:,:,:)= 0.0_dp
    wisoqvim(:,:,:)   = 0.0_dp
    wisoxlvim(:,:,:)  = 0.0_dp
    wisoxivim(:,:,:)  = 0.0_dp
    wisowl(:,:,:)     = 0.0_dp
    wisogld(:,:,:)    = 0.0_dp
    wisosnc(:,:,:)    = 0.0_dp
    wisoxtec(:,:,:,:) = 0.0_dp

    wisodrainm(:,:,:) = 0.0_dp
    wisosnaclm(:,:,:) = 0.0_dp

    wisoevapwac(:,:,:)   = 0._dp
    wisoevapiac(:,:,:)   = 0._dp
    wisoevaplac(:,:,:)   = 0._dp
    wisoapmeb(:,:,:)     = 0.0_dp
    wisoapmebco(:,:,:)   = 0.0_dp
    wisoqtnew(:,:,:)     = 0.0_dp
    wisorain(:,:,:)   = 0.0_dp

! variables for exchange with HD model and glacier calving model

    IF(lhd) THEN
      wisoaros(:,:,:)      = 0._dp
      wisoadrain(:,:,:)    = 0._dp
      wisodisch(:,:,:)     = 0._dp
      wisoapmecal(:,:,:)   = 0._dp
    END IF 
    
  END SUBROUTINE init_g3_wiso


  SUBROUTINE resume_g3_wiso ()

    !
    ! resume_g3_wiso - set rerun values for some further water isotope fields
    !
    ! Initialize *wiso* variables from default ECHAM restart fields
    ! (assume at start a delta value of zero permill) 
    !

    USE mo_control,       ONLY: lhd
    USE mo_memory_g3a 
    USE mo_memory_g3b 
    USE mo_surface_memory
    USE mo_memory_wiso

    IMPLICIT NONE

    !  Executable Statements

    ! set variables of *echam_wiso* stream
    
    DO jt=1,nwiso
      wisoaprl(:,jt,:)     = tnat(jt)*aprl(:,:)
      wisoaprc(:,jt,:)     = tnat(jt)*aprc(:,:)
      wisoaprs(:,jt,:)     = tnat(jt)*aprs(:,:)
      wisoevap(:,jt,:)     = tnat(jt)*evap(:,:)
      wisorunoff(:,jt,:)   = tnat(jt)*runoff(:,:)
      wisosnmel(:,jt,:)    = tnat(jt)*snmel(:,:)
      wisoruntoc(:,jt,:)   = tnat(jt)*runtoc(:,:)
      wisoapmegl(:,jt,:)   = tnat(jt)*apmegl(:,:)
      wisoqvi(:,jt,:)      = tnat(jt)*qvi(:,:)
      wisoxlvi(:,jt,:)     = tnat(jt)*xlvi(:,:)
      wisoxivi(:,jt,:)     = tnat(jt)*xivi(:,:)
      wisowl(:,jt,:)       = tnat(jt)*wl(:,:)
      wisogld(:,jt,:)      = tnat(jt)*gld(:,:)
      wisosnc(:,jt,:)      = tnat(jt)*snc(:,:)
      wisoxtec(:,:,jt,:)   = tnat(jt)*xtec(:,:,:)
    
      wisodrain(:,jt,:)    = tnat(jt)*drain(:,:)
      wisosnacl(:,jt,:)    = tnat(jt)*snacl(:,:)

      wisoevapwac(:,jt,:)  = tnat(jt)*evapwac(:,:)
      wisoevapiac(:,jt,:)  = tnat(jt)*evapiac(:,:)
      wisoevaplac(:,jt,:)  = tnat(jt)*evaplac(:,:)
      wisoapmeb(:,jt,:)    = tnat(jt)*apmeb(:,:)
      wisoapmebco(:,jt,:)  = tnat(jt)*apmebco(:,:)
      wisoqtnew(:,jt,:)    = tnat(jt)*qtnew(:,:)
      wisorain(:,jt,:)     = tnat(jt)*rain(:,:)     
    END DO

    ! variables for exchange with HD model and glacier calving model
    IF(lhd) THEN
      DO jt=1,nwiso
        wisoaros(:,jt,:)    = tnat(jt)*aros(:,:)
        wisoadrain(:,jt,:)  = tnat(jt)*adrain(:,:)
        wisodisch(:,jt,:)   = tnat(jt)*disch(:,:)
        wisoapmecal(:,jt,:) = tnat(jt)*apmecal(:,:)
      END DO
    END IF 


    ! set variables of *surf_wiso* stream

    DO jt=1,nwiso
      box%wiso_evaporation_inst(:,jt,:)   = tnat(jt)*box%evaporation_inst(:,:)
      box%wiso_evaporation_acc(:,jt,:)    = tnat(jt)*box%evaporation_acc(:,:)

      land%wiso_evaporation_inst(:,jt,:)  = tnat(jt)*land%evaporation_inst(:,:)
      land%wiso_evaporation_acc(:,jt,:)   = tnat(jt)*land%evaporation_acc(:,:)
      ocean%wiso_evaporation_inst(:,jt,:) = tnat(jt)*ocean%evaporation_inst(:,:)
      ocean%wiso_evaporation_acc(:,jt,:)  = tnat(jt)*ocean%evaporation_acc(:,:)
      ice%wiso_evaporation_inst(:,jt,:)   = tnat(jt)*ice%evaporation_inst(:,:)
      ice%wiso_evaporation_acc(:,jt,:)    = tnat(jt)*ice%evaporation_acc(:,:)

      land%wiso_evaporation_pot(:,jt,:)   = tnat(jt)*land%evaporation_pot(:,:) 
      land%wiso_surface_qsat_new(:,jt,:)  = tnat(jt)*land%surface_qsat_new(:,:)
      land%wiso_surface_qsat(:,jt,:)      = tnat(jt)*land%surface_qsat(:,:)

      land%zwisoeqnl(:,jt,:)              = tnat(jt)*land%zeqnl(:,:)
      land%zwisofqnl(:,jt,:)              = tnat(jt)*land%zfqnl(:,:)
      ocean%zwisoeqnw(:,jt,:)             = tnat(jt)*ocean%zeqnw(:,:)
      ocean%zwisofqnw(:,jt,:)             = tnat(jt)*ocean%zfqnw(:,:)
      ice%zwisoeqni(:,jt,:)               = tnat(jt)*ice%zeqni(:,:)
      ice%zwisofqni(:,jt,:)               = tnat(jt)*ice%zfqni(:,:)

      land%zwisoqsl(:,jt,:)               = tnat(jt)*land%zqsl(:,:)
      land%zwisoqklevl(:,jt,:)            = tnat(jt)*land%zqklevl(:,:)
      ocean%zwisoqsw(:,jt,:)              = tnat(jt)*ocean%zqsw(:,:)
      ocean%zwisoqklevw(:,jt,:)           = tnat(jt)*ocean%zqklevw(:,:)
      ice%zwisoqsi(:,jt,:)                = tnat(jt)*ice%zqsi(:,:)
      ice%zwisoqklevi(:,jt,:)             = tnat(jt)*ice%zqklevi(:,:)

      land%zwisocair(:,jt,:)              = tnat(jt)*land%zcair(:,:)
      land%zwisocsat(:,jt,:)              = tnat(jt)*land%zcsat(:,:)

      jwisorsfl(:,jt,:)                   = tnat(jt)*jrsfl(:,:)
      jwisorsfc(:,jt,:)                   = tnat(jt)*jrsfc(:,:)
      jwisossfl(:,jt,:)                   = tnat(jt)*jssfl(:,:)
      jwisossfc(:,jt,:)                   = tnat(jt)*jssfc(:,:)
    END DO
    
  END SUBROUTINE resume_g3_wiso


END SUBROUTINE iniwiso
