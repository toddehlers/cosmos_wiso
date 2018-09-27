     PROGRAM MPIOM
!
!    MPIOM  
!
!     Version:
!     ---------
!     $URL: http://svn.zmaw.de/svn/cosmos/branches/cosmos-landveg/src/mod/mpiom/src/mpiom.f90 $
!     $Rev: 1334 $
!
!    VERSION OF THE HOPE OGCM ON C-GRID
!
!    DEVELOPED BY MAIER-REIMER TILL 1997
!
!    VERSION OF UWE MIKOLAJEWICZ,JOHANN JUNGCLAUS AND HELMUTH HAAK, 12/99
!    INCLUDES CONFORMAL MAPPING     
!
!    MODIFIED :
!    ---------
!   HOPS65 : CREATED JULI 25, 2000  H. HAAK, J. JUNGCLAUS
!   HOPS67 : CREATED  NOV 09, 2000  H. HAAK, J. JUNGCLAUS, U.MIKOLAJEWICZ
!   HOPS68 : CREATED JUNE 20, 2001  H. HAAK, J. JUNGCLAUS, U.MIKOLAJEWICZ
!                        Nov. 2001  O. Boeringer  - TRANSFER TO FORTRAN90
!                                   S. Legutke    - interface to HAMOCC5
!   HOPS69 : CREATED   JAN 5, 2002  H. HAAK, J. JUNGCLAUS, U.MIKOLAJEWICZ
!   MPI-OM :           JAN 14 2003  S. Legutke    - created mo_couple.F90
!   MPI-OM :           JUN 03 2003  J. Jungclaus  - update GIRIV, RELTEM, ZOCORR
!   MPI-OM :           AUG 15 2007  S. Lorenz     - omp-parallel tracer loops
!   MPI-OM :           JUN 23 2010  X. Xu         - isotope tracer module
!*************************************************************************
!
!
!                       MP MPIOM
!                      SBR BELEG
!                          BODEN
!                          CORIOL
!                          ITPREP
!                            !
!                 !  --->  THERMODYNAMIC FORCING
!                 !        WIND FORCING
!      TIME       !        DECOMPOSITION INTO BAROTROPIC AND BAROCLINIC FIELD
!      STEPPING   !        BAROTROPIC SYSTEM
!                 !        BAROCLINIC SYSTEM
!                 !        MOMENTUM ADVECTION
!                 !        TRACER ADVECTION
!                 !        TRACER EDDY DIFFUSION
!                 !        TRACER DIFFUSION
!                 !        MOMENTUM DIFFUSION
!                 ! <---     !
!                            !
!                          OUTPUT-ROUTINES

!          ------------------------------
! NUMBER OF VERTICAL LAYERS IS KE.
!
!************************************************************
! PARAMETER :  IE   NUMBER OF GRID POINTS IN X 
!              JE                            Y
!              KE                            Z
!              KBB=NMAX,ILL=MATR MUST BE SET BY USER!!!
!
! SOME IMPORTANT VARIABLES AND FIELDS :
!             DT        TIME STEP
!             TIESTU    DEPTH OF HORIZONTAL VELOCITY POINTS
!             TIESTW    DEPTH OF VERTICAL VELOCITY POINTS
!             ZO        SEA SURFACE ELEVATION
!             AMSUE/O   LAND/SEA-MASK FOR VECTORFIELD (LAND=0/SEA=1)
!             WETO                    FOR SCALARFIELD       " 
!             UKO       ZONAL VELOCITY COMPONENT
!             VKE       MERIDIONAL "       "
!             WO        VERTICAL VELOCITY COMPONENT
!             THO       TEMPERATURE
!             SAO       SALINITY
!             PO        PRESSURE
!             TXO       ZONAL WIND STRESS
!             TYE       MERIDIONAL WIND STRESS
!             AVO       VARIABLE VERTICAL EDDY VISCOSITY
!                       (DEFINED ON SCALAR POINTS!!)
!             DVO       VARIABLE VERTICAL EDDY DIFFUSIVITY
!
!
! E=MERIDIONAL VELOCITY POINTS, O=ZONAL VELOCITY POINTS/SCALAR POINTS
!***********************************************************************
!                                                                      *
!UWE  AKTUALISIEREN!!!!!!!!!!!!!!!!
!
!     DATA SETS AND FORTRAN CHANNEL NUMBER CODING                      *
!                                                                      *
!     UNIT NO.   NAME     DESCRIPTION                      LOCATION    *
!
!HH    IO_IN_Z370 Z37000  RESTART-INFORMATION   MPIOM.F90,SBR AUFR,SBR AUFW
!HH    IO_IN_Z380 Z38000  RESTART-INFORMATION   MPIOM.F90,SBR AUFR,SBR AUFW
!HH    IO_IN_INIT INITEM  TEMPERATURE LEVITUS          MPIOM.F90,SBR LEVIRE
!HH    IO_IN_INIS INISAL  SALINITY LEVITUS             MPIOM.F90,SBR LEVIRE
!HH    IO_IN_SURS SURSAL  SEA SURFACE SALINITY LEVITUS MPIOM.F90,SBR LEVIRE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_LEVITUS
      USE MO_COMMOAU1
      USE MO_COMMOAU2
      USE MO_COMMOAU3
      USE MO_DIAGNOSIS
      USE MO_OCTHER
      USE MO_OCICE  
      USE MO_TRO


#ifdef CORE
     USE MO_NCAR_OCEAN_FLUXES
#else
     USE MO_OMIP
#endif

#ifdef NUDGE_ECCO
     USE MO_NUDGE_TS
#endif

#ifdef TIDAL
      USE MO_TIDAL
#endif
 
#ifdef FLUXCORRECT
      USE MO_COMMO_FLUXCORR
#endif /*FLUXCORRECT*/

      USE MO_ADPO
      USE MO_COMMOBBL

      USE MO_MEAN

      USE MO_ELICOM
      USE MO_PARA2

      USE MO_UNITS
      USE MO_OCECTLCOM
      USE MO_OCTDIFF

#if defined (__coupled)
      USE mo_fluxes1, ONLY : alloc_mem_fluxes1
      USE mo_couple

!, ONLY: couple_prep, couple_correct_ini, couple_init, &
!                           couple_get_a2o, couple_put_o2a, couple_calendar , &
!                           couple_end
#endif

#ifdef FB_BGC_OCE
#ifndef PBGC
      USE mo_attenmap
#endif
#endif

#ifdef PBGC

!#ifndef __coupled
!      USE mo_fluxes1
!#endif
      USE mo_carbch
      USE mo_control_bgc
      USE mo_biomod
      USE mo_sedmnt
      USE mo_param1_bgc 

!      REAL, POINTER :: layer1_bgc(:,:)
!      REAL, POINTER :: layer1_new(:,:)
      REAL, POINTER :: bgcddpo(:,:,:)
      REAL, POINTER :: bgcdpio(:,:,:)
      REAL vol_old,vol_new
      INTEGER ndtrun
#endif /*PBGC*/

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef ADDCONTRA


      USE mo_contra
      USE mo_control_add     
      USE mo_param1_add
      USE mo_omip_add 


      REAL, POINTER :: add_ddpo(:,:,:)
      REAL, POINTER :: add_dpio(:,:,:) 
      REAL add_vol_old, add_vol_new
      INTEGER add_ndtrun
#endif /*ADDCONTRA*/ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
!     USE MO_OCECTL        module removed, source included below
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     DECLARATIONS
      INTEGER*8 IDATE,III
      INTEGER NACTYEAR
      REAL ABACK,CAULAPTS,CAULAPUV,CRELSAL,CRELTEM,CWT,CWA,DBACK,GFDL_DIFF
      REAL PIBOG,RRELSAL,RRELTEM,DIST
      REAL AWERT, DDPOMAX, DDPOMIN, EWERT


      REAL CDZW(500)    
 
      INTEGER I,ICONVA,IWIRB,J,K,M,N,NANF,NDTDAY,ENDDAY
      INTEGER ICOU,IJJU,JMANF,JMEND,JMM,JMMM,LDTRUN
      INTEGER IERR, LDTDAY, LEN_FEB, LMON1, LMON2, LREAD
      INTEGER LREADMAX, MONMON, NACYEAR, NDTYEAR, IMAL
      INTEGER L,IREADC
      REAL*8 ttts, tttr, ttt

#ifdef SORTEST
      REAL,allocatable :: z1(:,:),z2(:,:)
#endif

#ifdef CLOCK
  REAL :: zwtime,zutime,zstime, zrtime

!  External functions 
  REAL, EXTERNAL :: util_walltime
  INTEGER, EXTERNAL :: util_cputime
#endif



! OtB  cannot use variables in more than one module
! ---------------------------------------------------------------------
!
!*    *NAMELIST* *OCECTL*   - Defines namelist parameters.
!                           - included in *MPIOM*.
!
!*     VARIABLE  TYPE        PURPOSE.
!      --------  ----        --------
!      *DT*      *REAL*      Ocean time step in sec
!
! ---------------------------------------------------------------------
!
      NAMELIST /OCEDIM/ IE_G,JE_G,KE     !declared in mo_param1

      NAMELIST /OCECTL/                                                &
              DT                                                       &
             ,CAULAPTS, CAULAPUV, CAH00, AUS                           &
             ,AV0,DV0,CWT,CWA,CSTABEPS,DBACK,ABACK,CRELSAL,CRELTEM     &
             ,ICONVA, IWIRB                                            &
             ,RRELTEM, RRELSAL                                         &
             ,CDVOCON,CAVOCON,IBOLK                                    &
             ,IOCAD,isnflg                                             &
             ,H0, HMIN, ARMIN, ARMAX, HSNTOICE, SICTHMIN, SICE         &
             ,D3                                                       &
             ,IAUFR, IAUFW                                             &
             ,ISTART,I3DREST                                           &
             ,NYEARS, NMONTS, NDAYS, IMEAN, LY_END, LY_START,LM_START  &
             ,EXPTID,ICONTRO,IOASISFLUX,imocdiag                       &
             ,LFORCEDIAG,LHFLDIAG,LCONVDIAG,LDIFFDIAG,LGMDIAG          &
             ,LGRIDINFO,LCALCDIFI                                      &
             ,ITSDIAG,LTSTRANSPOSE

      NAMELIST /OCEDZW/ CDZW

      NAMELIST /NPROCS/ nprocx, nprocy ! declared in mo_parallel

#ifdef __IFC /* Intel Fortran Compiler */
      INTEGER ieee_handler
      EXTERNAL handler
!
      n = ieee_handler('set','division',handler)
      n = ieee_handler('set','overflow',handler)
      n = ieee_handler('set','invalid',handler)
#endif


#ifdef CLOCK
! Initialize wallclock timer
  zwtime = util_walltime()
#endif

!     Initialize MPI

      CALL p_start

!     SPECIFY LOGICAL I/O UNITS

      CALL SETUNITS
!  OPEN STD OUT FILE
!  This is needed for the MPI version since we do not want
!  the output of all processors intermixed (and basically
!  we want only the ouput of processor 0)

      IF (io_stdout > 0) THEN
        CALL OPEN_STDOUT(io_stdout,'oceout')
      ENDIF

!     Open input file and read the number of processors along
!     x- and y-direction

      IF(p_pe==p_io) THEN
        OPEN(IO_IN_OCTL,FILE='OCECTL',STATUS='UNKNOWN',                 &
     &          ACCESS='SEQUENTIAL',FORM='FORMATTED')

        READ(IO_IN_OCTL,OCEDIM)   ! read dimensions
!        write(0,*)'read dimensions'
        READ(IO_IN_OCTL,NPROCS)
!	 write(0,*)'read decomposition'
      ENDIF

      CALL p_bcast(ie_g,p_io)
      CALL p_bcast(je_g,p_io)
      CALL p_bcast(ke,p_io)

!     Domain decomposition, setting local IE and JE
!       write(0,*)'vor deco'

      CALL p_deco
!       write(0,*)'nach deco',ie_g,je_g,ke,ie,je


!     Set some dependent parameters and allocate arrays

      CALL set_param1
      CALL alloc_mem_commo1
      CALL alloc_mem_octdiff
      CALL alloc_mem_commoau2
      CALL alloc_mem_commoau3
      CALL alloc_mem_diag

      CALL alloc_mem_dilcor



#ifdef __coupled
      CALL alloc_mem_fluxes1
#endif /* __coupled */

!#ifdef CONVDIAG
!      CALL alloc_mem_commconv
!#endif /*CONVDIAG*/

#ifdef PBGC
!      ALLOCATE( layer1_bgc(ie,je) )
!      ALLOCATE( layer1_new(ie,je) )
      ALLOCATE( bgcddpo(ie,je,ke) )
      ALLOCATE( bgcdpio(ie,je,ke) )
#endif /*PBGC*/      
!    write(0,*)'vor elicom'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef ADDCONTRA

      ALLOCATE( add_ddpo(ie,je,ke) )
      ALLOCATE( add_dpio(ie,je,ke) )
#endif /*ADDCONTRA*/      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CALL alloc_mem_elicom
      CALL alloc_mem_para2

#ifdef bounds_exch_put
      CALL set_put_window
#endif
#ifdef CORE
       WRITE(IO_STDOUT,*) 'allocate core'
       CALL alloc_mem_core
#endif
!      write(0,*)'vor beleg'

      CALL BELEG_ZERO


!  DEFAULT LENGTH OF INTEGRATION

!      write(0,*)'vor namelist'


      NYEARS=0
      NMONTS=1


      NANF=0
      NNNDT = 0
      NDAYS = 0
      DH = 0.
      AH = 0.

      ALMZER=1.E-19

      monlen(:) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

#ifdef YEAR360
      monlen(:) = 30
#endif /*YEAR360*/

#ifdef DEBUG_ONEDAY
      monlen(:) = 1
#endif /*DEBUG_ONEDAY*/

      WINTUR(1)=0.
      WINTUR(2)=1.E-3
      DO K=3,KEP
         WINTUR(K)=0.3*WINTUR(K-1)
      ENDDO

!     DZwW(K) DECOFTES THE THICKNESS OF LAYER K IN METERS
!     SUM OF DZW(1-KE) IS THE BOTTOM OF THE MODEL, I.E. SOLID GROUND !!!
!     DZ(K) DECOFTES THE DISTANCE BETWEEN VECTOR POINTS BUT WILL BE
!     COMPUTED IN SBR BODEN FOR CONSISTENCY WITH W-POINTS
      DO K=1,KE
         SAF(K)=34.8
         TAF(K)=1.
         DZW(K)=0.
      ENDDO  
!
!----------------------------------------------------------------------
! DEFAULT PARAMETER SETTINGS - CAN BE OVERWRITEN BY THE NAMELIST
! ISTART :: START OPTIONS
! ISTART == 0: COMPLETELY NEW SETUP, 
!              topography read from anta and written to topo (!!! WARNING !!!)
!              start from climatology 
! ISTART == 1: new run, topography read from topo
!              start from horizontally uniform ts-profile (taf,saf)
! ISTART == 2: new run, topography read from topo
!              start from climatology
! ISTART == 3: continuing run (default)
         ISTART =3
!
! I3DREST OPTIONS FOR 3-D RESTORING
! I3DREST == 0: NO RESTORINg (DEFAULT)
! I3DREST == 1: RESTORING to annual climatology
! I3DREST == 2: RESTORING to monthly climatology
         I3DREST=0
!--------------------------------------------------------------------
!  TRACER ADVECTION ROUTINES : SEVERAL OPTIONS 
!  IOCAD ==  1: UPWIND  
!  IOCAD ==  2: not used
!  IOCAD ==  3: ADPO 
!  IOCAD ==  4: ADPO + SLOPECON_ADPO
!  IOCAD ==  5: ADFS
!  IOCAD ==  6: not yet --> QUICK
!  IOCAD ==  7: mot yet --> QUICK2

         IOCAD=4

        ICONTRO=0

!HH    DEFAULT TIME STEPS
         DT=1800.

!UWE   CONSTANTS FOR BIHARMONIC DIFFUSION
         CAULAPTS=0.0002

!UWE   CONSTANTS FOR HARMONIC DIFFUSION
         CAH00=0.

#ifdef ISOPYK
         CAULAPTS=0.
#endif

!UWE   CONSTANTS FOR BIHARMONIC FRICTION
         CAULAPUV=0.0045

!UWE   CONSTANTS FOR HARMONIC FRICTION
         AUS=3.E-6

#ifdef ISOPYK
         CAH00=1000.
#endif

!HH    CONSTANT FOR DIFFUSION IN OCTHER  
        DV0=0.5E-2
        AV0=0.5E-2
        CDVOCON=20.
        CAVOCON=0.

        IBOLK=1000

!HH    CONSTANT FOR WINDMIXING IN OCTHER 
        CWT=0.5E-3
        CWA=-9e33   ! cwa is set to cwt if the special value 
                    !  is not overwritten in the namelist   
        CSTABEPS=0.05

!HH    BACKGROUNDDIFFUSION
        ABACK=5.E-5
        DBACK=5.E-5

!HH    DEFAULT MEAN OUTPUT
       IMEAN=2

!HH    RELAXATION TIME SALINITY 
        CRELSAL = 3.E-7

!HH    RELAXATION TIME TEMERATURE
        CRELTEM = 0.

!JJ    DEFAULT END OF RUN IN YEAR
       LY_END=9999

!JJ    DEFAULT OFFSET FOR YEAR COUNTER (NEG. VALUE == NO ACTION TAKEN)
       LY_START=-999
       LM_START=-999
!----------------------------------------------------------------------
!     ICONVA : 1   WITH CONVECTIVE ADJUSTMENT
!              0    NO      "          "
      ICONVA=1
!-------------------------------------------------------
!     IWIRB  : 1   COMPUTE VARIABLE EDDY VISCOSITY
!              0   CONSTANT EDDY VISCOSITY
      IWIRB=1
!----------------------------------------------------------------------
! SWITCH FOR RESPECTIVE RESTART FILES
!
      IFLAG=1

!----------------------------------------------------------------------

!  SWITCH FOR WRITING OASIS COUPLED FLUXES

    IOASISFLUX=0
    IMOCDIAG=0

    LGMDIAG=.FALSE.
    LFORCEDIAG=.FALSE.
    LHFLDIAG=.FALSE.
    LCONVDIAG=.FALSE.
    LDIFFDIAG=.FALSE.
    LGRIDINFO=.FALSE.
    LCALCDIFI=.FALSE.

!----------------------------------------------------------------------
!FR   Initialize timeseries writing
! itsdiag=0     :  No Output
! itsdiag=1     :  one snapshot per day
! itsdiag=2     :  monthly averaged snapshots
! itsdiag=3     :  yearly  averaged snapshots
! itsdiag=4     :  output every timestep
! itsdiag=5     :  daily   average
! itsdiag=7     :  monthly average of daily means
! itsdiag=7     :  yearly  average of daily means 
! ltstranspose=T:  one value for each code in diagnosis (extra-format)
! ltstranspose=F:  all codes in one array  in diagnosis (pseudo-extra)

      itsdiag=1
      ltstranspose=.TRUE.
      ltswrite=.TRUE.     !  set to .FALSE. after first call of diagnosis


!  READ OCEAN NAMELIST

      IF(p_pe==p_io) THEN

!        write(0,*)'vor ocectl'
        READ(IO_IN_OCTL,OCECTL)
!        Write(0,*)'after ocectl'
        READ(IO_IN_OCTL,OCEDZW)
!	Write(0,*)'after ocedzw'
        CLOSE(IO_IN_OCTL)

      ENDIF


      if (imean.eq.0) then

	IF ( LGMDIAG ) Write(0,*)'imean = 0 => LGMDIAG disabled'
        LGMDIAG=.FALSE.
	IF ( LFORCEDIAG ) Write(0,*)'imean = 0 => LFORCEDIAG disabled'
        LFORCEDIAG=.FALSE.
	IF ( LHFLDIAG ) Write(0,*)'imean = 0 => LHFLDIAG disabled'
	LHFLDIAG=.FALSE.
#ifdef __coupled
	IF ( LHFLDIAG ) Write(0,*)'coupled => LHFLDIAG disabled'
	LHFLDIAG=.FALSE.
#endif   
	IF ( LCONVDIAG ) Write(0,*)'imean = 0 => LCONVDIAG disabled'
	LCONVDIAG=.FALSE.
	IF ( LDIFFDIAG ) Write(0,*)'imean = 0 => LDIFFDIAG disabled'
	LDIFFDIAG=.FALSE.
        IF ( LGRIDINFO ) Write(0,*)'imean = 0 => LGRIDINFO disabled'
	LGRIDINFO=.FALSE.
        IF ( LCALCDIFI ) Write(0,*)'imean = 0 => LCALCDIFI disabled'
	LCALCDIFI=.FALSE.
      endif

      IF ( IBOLK .EQ. 0 )  LGMDIAG=.FALSE.
      IF (p_pe==p_io) THEN
        IF ( IBOLK .EQ. 0 ) Write(0,*)'ibolk = 0 => GM scheme disabled'
        IF ( IBOLK .GT. 0 ) Write(0,*)'ibolk > 0 => GM scheme enabled'
#ifndef ISOPYK
        IF ( IBOLK .GT. 0 ) Write(0,*)                                  &
           'ATTN => GM scheme enabled with horizontal ',                &
           ' instead of isopycnal diffusion'
#endif
        IF ( IBOLK .LT. 0 ) Write(0,*)                                  &
           'ibolk < 0 => GM + Visbek scheme enabled'
#ifndef ISOPYK
        IF ( IBOLK .LT. 0 ) Write(0,*)                                  &
           'ATTN => GM+Visbek scheme enabled with horizontal ',         &
           ' instead of isopycnal diffusion'
#endif
      ENDIF


      CALL p_bcast(DT,p_io)
      CALL p_bcast(CAULAPTS,p_io)
      CALL p_bcast(CAULAPUV,p_io)
      CALL p_bcast(CAH00,p_io)
      CALL p_bcast(AUS,p_io)
      CALL p_bcast(AV0,p_io)
      CALL p_bcast(DV0,p_io)
      CALL p_bcast(CWT,p_io)
      if (cwa.lt.0.) cwa=cwt
      CALL p_bcast(CWA,p_io)
      CALL p_bcast(CSTABEPS,p_io)
      CALL p_bcast(DBACK,p_io)
      CALL p_bcast(ABACK,p_io)
      CALL p_bcast(CRELSAL,p_io)
      CALL p_bcast(CRELTEM,p_io)
      CALL p_bcast(ICONVA,p_io)
      CALL p_bcast(IWIRB,p_io)
      CALL p_bcast(RRELTEM,p_io)
      CALL p_bcast(RRELSAL,p_io)
      CALL p_bcast(CDVOCON,p_io)
      CALL p_bcast(CAVOCON,p_io)
      CALL p_bcast(IBOLK,p_io)
      CALL p_bcast(ISNFLG,p_io)
      CALL p_bcast(IOCAD,p_io)
      CALL p_bcast(ICONTRO,p_io)
      CALL p_bcast(H0,p_io)
      CALL p_bcast(HMIN,p_io)
      CALL p_bcast(ARMIN,p_io)
      CALL p_bcast(ARMAX,p_io)
      CALL p_bcast(HSNTOICE,p_io)
      CALL p_bcast(SICTHMIN,p_io)
      CALL p_bcast(SICE,p_io)
      CALL p_bcast(D3,p_io)
      CALL p_bcast(IAUFR,p_io)
      CALL p_bcast(IAUFW,p_io)
      CALL p_bcast(NYEARS,p_io)
      CALL p_bcast(NMONTS,p_io)
      CALL p_bcast(NDAYS,p_io)
      CALL p_bcast(IMEAN,p_io)
      CALL p_bcast(LY_END,p_io)
      CALL p_bcast(LY_START,p_io)
      CALL p_bcast(LM_START,p_io)
      CALL p_bcast(EXPTID,p_io)
      CALL p_bcast(ISTART,p_io)
      CALL p_bcast(I3DREST,p_io)
      CALL p_bcast(IOASISFLUX,p_io)
      CALL p_bcast(IMOCDIAG,p_io)
      CALL p_bcast(LFORCEDIAG,p_io)
      CALL p_bcast(LHFLDIAG,p_io)
      CALL p_bcast(LCONVDIAG,p_io)
      CALL p_bcast(LDIFFDIAG,p_io)
      CALL p_bcast(LGMDIAG,p_io)
      CALL p_bcast(LGRIDINFO,p_io)
      CALL p_bcast(LCALCDIFI,p_io)
      CALL p_bcast(ITSDIAG,p_io)
      CALL p_bcast(LTSTRANSPOSE,p_io)
      CALL p_bcast(LTSWRITE,p_io)

      if ( p_pe==p_io ) then
      do k=1,ke
              DZW(k)=CDZW(k)
!         write(0,*) k,cdzw(k),dzw(k)
 
      enddo
      endif

      CALL p_bcast(dzw,p_io)

      DTI=1./DT
      NDTDAY=NINT(86400./DT)

! Write final namelist parameters
!
      WRITE(IO_STDOUT,*)' SECONDS PER TIMESTEP   (DT): ',DT
      WRITE(IO_STDOUT,*)' TIME STEPS PER DAY (NDTDAY): ',NDTDAY
      WRITE(IO_STDOUT,*)'                     (IOCAD): ',IOCAD
      WRITE(IO_STDOUT,*)'                  (CAULAPUV): ',CAULAPUV
      WRITE(IO_STDOUT,*)'                  (CAULAPTS): ',CAULAPTS
      WRITE(IO_STDOUT,*)'                       (AUS): ',AUS
      WRITE(IO_STDOUT,*)'                     (CAH00): ',CAH00
      WRITE(IO_STDOUT,*)'                       (AV0): ',AV0
      WRITE(IO_STDOUT,*)'                       (DV0): ',DV0
      WRITE(IO_STDOUT,*)'                       (CWT): ',CWT
      WRITE(IO_STDOUT,*)'                       (CWA): ',CWA
      WRITE(IO_STDOUT,*)'                  (CSTABEPS): ',CSTABEPS
      WRITE(IO_STDOUT,*)'                     (DBACK): ',DBACK
      WRITE(IO_STDOUT,*)'                     (ABACK): ',ABACK
      WRITE(IO_STDOUT,*)'                   (CRELSAL): ',CRELSAL
      WRITE(IO_STDOUT,*)'                   (CRELTEM): ',CRELTEM
      WRITE(IO_STDOUT,*)'                   (CDVOCON): ',CDVOCON
      WRITE(IO_STDOUT,*)'                   (CAVOCON): ',CAVOCON
      WRITE(IO_STDOUT,*)'                     (IBOLK): ',IBOLK
      WRITE(IO_STDOUT,*)'                     (IMEAN): ',IMEAN
      WRITE(IO_STDOUT,*)'                  (LY_START): ',LY_START
      WRITE(IO_STDOUT,*)'                  (LM_START): ',LM_START
      WRITE(IO_STDOUT,*)'                    (LY_END): ',LY_END
      WRITE(IO_STDOUT,*)'                    (ISTART): ',ISTART
      WRITE(IO_STDOUT,*)'                   (I3DREST): ',I3DREST
      WRITE(IO_STDOUT,*)'                   (ICONTRO): ',ICONTRO
      WRITE(IO_STDOUT,*)'                (IOASISFLUX): ',IOASISFLUX
      WRITE(IO_STDOUT,*)'                  (IMOCDIAG): ',IMOCDIAG
      WRITE(IO_STDOUT,*)'                   (LGMDIAG): ',LGMDIAG
      WRITE(IO_STDOUT,*)'                (LFORCEDIAG): ',LFORCEDIAG
      WRITE(IO_STDOUT,*)'                  (LHFLDIAG): ',LHFLDIAG
      WRITE(IO_STDOUT,*)'                 (LDIFFDIAG): ',LDIFFDIAG
      WRITE(IO_STDOUT,*)'                 (LCONVDIAG): ',LCONVDIAG
      WRITE(IO_STDOUT,*)'                 (LGRIDINFO): ',LGRIDINFO
      WRITE(IO_STDOUT,*)'                 (LCALCDIFI): ',LCALCDIFI
      WRITE(IO_STDOUT,*)'                   (ITSDIAG): ',ITSDIAG
      WRITE(IO_STDOUT,*)'              (LTSTRANSPOSE): ',LTSTRANSPOSE

      IF (ISTART .LT. 0 .OR. ISTART .GT. 3)                              &
          WRITE(IO_STDOUT,*)'ISTART NOT SUPPORTED!!!!!'
      IF (I3DREST .LT. 0 .OR. I3DREST .GT. 2)                            &
          WRITE(IO_STDOUT,*)'I3DREST NOT SUPPORTED!!!!!'     

      WRITE(IO_STDOUT,*)'THIS JOB WILL TRY TO INTEGRATE '               &
                        ,NYEARS,' YEARS AND '                           &
                        ,NMONTS,' MONTHS AND'                           &
                        ,NDAYS,'DAYS'




      IF (I3DREST .GT. 0 .OR. ISTART .LT. 3) THEN
       CALL init_levitus(ie,je,ke)
      ENDIF


!  DEFAULT ...
       IAUFR=1
       IAUFW=1
      IF (ISTART .LT. 3) IAUFR=0

      AULAPTS=CAULAPTS*DT/3600.
      AULAPUV=CAULAPUV*DT/3600.
      AH00=CAH00/4.E5
      WT=CWT/(6.**3)
      WA=CWA/(6.**3)
   
!HH   CHECK LAYER THICKNESS
      DO K=1,KE
        IF(DZW(K).EQ.0.) THEN
          WRITE(IO_STDOUT,*)' LAYER: ',K,' THICKNESS IS ZERO !!!'
          CALL ABSTURZ 
        ELSE
          WRITE(IO_STDOUT,*) K, ' LAYERTHICKNESS    (DZW): ',DZW(K)
        ENDIF
      ENDDO


!HH   BACKGROUNDDIFFUSION
      DO K=1,KEP
        ABACKV(K) = ABACK
        DBACKV(K) = DBACK 
      ENDDO

      TIESTW(1) = 0.
      DO K=1,KE
      TIESTW(K+1)  = TIESTW(K) + DZW(K)
      ENDDO
      PI=4.*ATAN(1.)
      PIBOG=180./PI

#ifdef DBACKPROFIL
      IDBACK=0
      DO K=1,KEP
         ABACKV(K) = ABACK
         DBACKV(K) = DBACK + (1-IDBACK) * 1.E-4 *                       &
     &               (0.8 + 1.05/PI*ATAN(4.5*1.E-3*(TIESTW(K)-2500.)))* &
     &               (0.5+SIGN(0.5,TIESTW(K)-500.)) *                   &
     &               SQRT(ABS((TIESTW(K)-500.)/(3500.-500.)))
      ENDDO
#endif

      DO K=1,KEP
         GFDL_DIFF = 1.E-4 *                                            &
     &            (0.8 + 1.05/PI*ATAN(4.5*1.E-3*(TIESTW(K)-2500.)))
#ifdef DBACKGFDL
         DBACKV(K) = GFDL_DIFF
#endif
       ENDDO

#ifdef DBACKGFDL2
         DBACKV(K) = 0.5*(DBACK+GFDL_DIFF)
#endif
       DO K=1,KEP
         WRITE(IO_STDOUT,6002)'BACKGROUND DIFFUSIVITY AT '              &
     &                        ,INT(TIESTW(K))                           &
     &         ,'M : HOPE : ',DBACKV(K),' GFDL : ',GFDL_DIFF
         WRITE(IO_STDOUT,6002)'BACKGROUND VISCOSITY AT ',INT(TIESTW(K)) &
     &         ,'M : HOPE : ',ABACKV(K)
       ENDDO

 6002  FORMAT(1X,A27,I5,A10,E12.3,A8,E12.3)

!HH    RELAXATION TIME SALINITY 
!      RELSAL = CRELSAL*20./DZW(1)
      RELSAL = CRELSAL
      IF (RELSAL.GT. ALMZER) THEN 
!         WRITE(IO_STDOUT,26668)1./(RELSAL*24.*3600.)
         if (p_pe==p_io) then
          WRITE(0,*) 'RELAXATION DZW(1) [m]  =',DZW(1)
          WRITE(0,*)' RELAXATION TIME [DAYS] =',1./(RELSAL*24.*3600.)
          WRITE(0,*)' PISTON VELOCITY [m/s]  =',DZW(1)*(RELSAL)
          WRITE(0,*)' RELAXATION TIME (relative to 20m) [DAYS] =',1./(RELSAL*24.*3600.*DZW(1)/20.)
       endif
      ELSE
         if (p_pe==p_io) then
            WRITE(IO_STDOUT,*) 'SSS relaxation switched off !!'
         endif
      ENDIF

!26668 FORMAT('  RELAXATION TIME SALINITY COUPLING : ',F10.2,' DAYS')
!HH    RELAXATION TIME TEMPERATURE
!      RELTEM = CRELTEM*20./DZW(1)
      RELTEM = CRELTEM
      IF (reltem.GT. ALMZER) THEN 
!         WRITE(IO_STDOUT,26669)1./(RELTEM*24.*3600.)
       if (p_pe==p_io) then
          WRITE(IO_STDOUT,*)'  RELAXATION TIME TEMPEARTURE COUPLING : '&
              ,DZW(1),' m /',(RELTEM*24.*3600.),' DAYS'
       endif
     ELSE
        if (p_pe==p_io) then
           WRITE(IO_STDOUT,*) 'SST relaxation switched off !!'
        endif
      ENDIF
!26669 FORMAT('  RELAXATION TIME SST COUPLING : ',F10.2,' DAYS')


!-----------------------------------------------------------------------

      IF (iocad .NE. 5) THEN
         CALL alloc_mem_adpo
      ENDIF
      IF (iocad .EQ. 4) THEN
         CALL alloc_mem_commobbl
      ENDIF
      
      IF (IBOLK .NE. 0) THEN
         CALL alloc_mem_gmbolus
      ENDIF

!    write(0,*)'vor mean'

      if(imean.ne.0) then
         CALL alloc_mem_mean
      endif


   
#if defined (__coupled)
!----------------------------------------------------------------------
! Prepare coupling
!
     CALL couple_prep
#endif

!----------------------------------------------------------------------
! OPEN FILES
      IF(p_pe==p_io) THEN
        OPEN(IO_IN_ARCG,FILE='arcgri'                                   &
     &                 ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

        OPEN(IO_IN_INIT,FILE='INITEM',STATUS='UNKNOWN',                 &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        OPEN(IO_IN_INIS,FILE='INISAL',STATUS='UNKNOWN',                 &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        OPEN(IO_IN_SURS,FILE='SURSAL',STATUS='UNKNOWN',                 &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        IF (RELTEM .GT. ALMZER) THEN
          OPEN(IO_IN_SURT,FILE='SURTEM',STATUS='UNKNOWN',               &
     &            ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        ENDIF

#if defined (__coupled)
#ifdef FLUXCORRECT
        CALL couple_correct_ini 
#endif /*FLUXCORRECT*/
#else /*(__coupled)*/


#ifdef CORE

       WRITE(IO_STDOUT,*) 'open core'
       CALL open_core 
#else 
      WRITE(IO_STDOUT,*) 'open omip'
      CALL open_omip
#endif

!------------------------------------------------------
!------------------------------------------------------
#ifdef ADDCONTRA
  
       WRITE(IO_STDOUT,*) 'open omip_add'
       CALL open_omip_add                
#endif
!------------------------------------------------------
!------------------------------------------------------

#endif /*(__coupled)*/

#ifdef AMOCEMR
!----------
        OPEN(IO_OU_SCHA,FILE='SCHALL',STATUS='UNKNOWN',                 &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
!---------
#endif /*AMOCEMR*/
      ENDIF ! p_pe == p_io

!-------------------------------------------------------------------

!     DU/DT - F*V = -G*( STABN* (DZ/DX)(NEW) + STABO* (DZ/DX)(OLD) )

!     DZ/DT = H*( CONN* (DU/DX)(NEW) + CONO* (DU/DX)(OLD) + ... )

      STABN=0.6
      STABO=1.-STABN

      CONN=0.50
      CONO=1.-CONN

      WRITE(IO_STDOUT,*)' STABN = ',STABN,' CONN = ',CONN
!
!-------------------------------------------------------------------
!  OMEGA  :  ANGULAR VELOCITY OF OUR NICE BLUE PLANET
!  RADIUS :  RADIUS OF THE ABOVE MENTIONED PLANET
!  G      :  GRAVITATIONAL ACCELERATION ( CONSTANT 9.81 M/S**2
!                TO ACCOUNT FOR THE EXISTENCE OF MINI BLACK HOLES )
!  ROCP   :  HEAT CAPACITY PER CUBICMETER
!  ROCD   :  WIND STRESS DRAG COEFFICIENT
!------------------------------------------------------------------
      PI=4.*ATAN(1.)
      OMEGA=7.292E-5
      RADIUS=6371.E+3
      G=9.81
      GHN=G*STABN
      GHO=G*STABO
      ROCP=4.E06
      ROCD=1.2*1.5E-3
!-------------------------------------------------------------------
! PARAMETER ICE MODEL 

      ACLO(:,:)=0.7
      PAO(:,:)=101300.
      RPRECO(:,:)=0.0
      FRSE(:,:)=0.0

      ISNFLG=1

      ALBI=0.75
      ALBM=0.66
      ALBW=0.10
      ALBSN=0.85
      ALBSNM=0.75
#ifdef ALBMSN07
      ALBSN=0.82
      ALBSNM=0.7
      ALBM=0.63
#endif /*ALBMSN07*/

#ifdef ALBOMIP
!SJM TUNE ALBEDOS FOR OMIP WITH BERYLIAND
      ALBI=0.75    !ICE < 0
      ALBM=0.70    !ICE > 0
      ALBW=0.10    !WATER
      ALBSN=0.85   !SNOW < 0
      ALBSNM=0.70  !SNOW > 0
#endif /*ALBOMIP*/
#ifdef ALBNCEP
!HH   TUNE ALBEDOS FOR NCEP WITH KOCH
      ALBI=0.76
      ALBM=0.71
      ALBW=0.1
      ALBSN=0.86
      ALBSNM=0.72
#endif
#ifdef ALBMELTHI
      ALBSNM=0.8
      ALBSN=0.87
      ALBI=0.78
      ALBM=0.7
#endif /*ALBMELTHI*/
      TMELT=273.16
      TFREZ=-1.9
      CC=4.2E6
      CW=0.0045
      CLO=3.02E8
      CLB=2.70E8
      RHOAIR=1.3E+00
      RHOWAT=1.025E+03
      RHOICE=0.91E+03
      RHOSNO=0.33E+03
      RHOICWA=RHOICE/RHOWAT
      RHOSNWA=RHOSNO/RHOWAT
      RHOSNIC=RHOSNO/RHOICE
      CON=2.1656
      CONSN=0.31
      H0=0.5
      ARMIN=0.15
      ARMAX=1.0
      HMIN=0.05
!UWE  MAXIMUM ICE THICKNESS FOR SNOW-->ICE CONVERSION
!      HSNTOICE = 17.
       HSNTOICE = 0.45 * DZW(1)
!UWE  MINIMUM ICE THICKNESS IN NEW ICE GROWTH
      SICTHMIN=0.5
!
      VAPL=2.5E6
      SUBL=2.834E6
      SICE=5.0
      D3=5.5E-08
!
!JJ OLD VALUES UP TO HOPS 62!
      D1=RHOAIR*1004.*1.75E-3
      D2I=RHOAIR*SUBL*1.75E-3
      D2W=RHOAIR*VAPL*1.75E-3

#ifdef DRAGGILL
!HH   VALUES FROM GILL ATMOSPHERE OCEAN DYNAMICS
      D1=RHOAIR*1004.*1.1E-3
      D2I=RHOAIR*SUBL*1.5E-3
      D2W=RHOAIR*VAPL*1.5E-3 
#endif
#ifdef DASILVA
      D1=RHOAIR*1004.
      D2I=RHOAIR*SUBL
      D2W=RHOAIR*VAPL
#endif
#ifdef BULK_KARA
      D1=RHOAIR*1004.67
      D2I=RHOAIR*SUBL
      D2W=RHOAIR*VAPL
#endif

      TFREEZ=TFREZ
      ENTMEL=320.E6
      SICHEC=2.
      HICCE=2.
      HICCP=20.
      HTH00=0.5
      STEBOL=5.67E-8

      N=JTO
      N1=N-1
      N2=N-2
      N3=N-3
      N4=N-4

      M=IE
      M1=M-1
      M2=M-2
      M3=M-3
      M4=M-4

      KB=KBB
      KM=KB+1
      KBM=KB+KM

      DEUTO(:,:)=0.
      DEUTE(:,:)=0.
      DEPTO(:,:)=0.

      IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
      CALL read_slice(IO_IN_ARCG,DEUTO)



      IF (ISTART .EQ. 0) THEN

      CALL bounds_exch('p',DEUTO,'mpiom 10')


         DO J=1,JE
            DO I=1,IE
               DEPTO(I,J)=DEUTO(I,J)
            ENDDO
         ENDDO


         DO I=1,IE
#ifndef bounds_exch_tp
            IF(have_g_js) THEN
               DEPTO(I,1)=0.
               DEPTO(I,2)=0.
            ENDIF
#endif
            IF(have_g_je) THEN
               DEPTO(I,JE)=0.
               DEPTO(I,JE1)=0.
            ENDIF
         ENDDO


       DO 8261 J=2,JE1
       DO 8261 I=2,IE1
          DEUTO(I,J)=DEPTO(I,J)
          IF(DEPTO(I,J).LT.1.) GO TO 8261
          IF(DEPTO(I-1,J).LT.1..AND.DEPTO(I+1,J).LT.1.                  &
     &        .AND.DEPTO(I,J-1).LT.1..AND.DEPTO(I,J+1).LT.1.)           &
     &        DEUTO(I,J)=0.
8261   CONTINUE

       CALL bounds_exch('u+',DEUTO,'mpiom 11')


       DO 8262 J=1,JE
       DO 8262 I=1,IE
          DEPTO(I,J)=DEUTO(I,J)
8262   CONTINUE

       ENDIF ! ISTART

       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, DLXP)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, DLXU)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, DLXV)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, DLYP)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, DLYU)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, DLYV)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, FTWOU)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, FTWOV)
       IF(p_pe==p_io) CLOSE(IO_IN_ARCG)




       DLXP(:,:)=MAX(1.,DLXP(:,:))
       DLXU(:,:)=MAX(1.,DLXU(:,:))
       DLXV(:,:)=MAX(1.,DLXV(:,:))
       DLYP(:,:)=MAX(1.,DLYP(:,:))
       DLYU(:,:)=MAX(1.,DLYU(:,:))
       DLYV(:,:)=MAX(1.,DLYV(:,:))


      DO J=1,JE
       DO I=2,IE1
        DLXPSI(I,J)=0.5*(DLXV(I,J)+DLXV(I+1,J))
        DLYPSI(I,J)=0.5*(DLYV(I,J)+DLYV(I+1,J))
       enddo
      enddo

       CALL bounds_exch('v+',DLXV,'mpiom 12')
       CALL bounds_exch('p',DLXP,'mpiom 13')
       CALL bounds_exch('u+',DLXU,'mpiom 14')
       CALL bounds_exch('p',DLYP,'mpiom 15')
       CALL bounds_exch('u+',DLYU,'mpiom 16')
       CALL bounds_exch('v+',DLYV,'mpiom 17')
       CALL bounds_exch('s',DLYPSI,'mpiom 17a')
       CALL bounds_exch('s',DLXPSI,'mpiom 17b')

!SL 
!SL GRID DEFORMATION
!SL 
       DO J=2,JE1
         DO I=2,IE1
             CURVAV(I,J)=(DLXP(I,J+1)-DLXP(I,J))/(DLYV(I,J)*DLXV(I,J))
         ENDDO
       ENDDO

       CALL bounds_exch('v+',CURVAV,'mpiom 18')


       DO 13788 I=1,IE
          IF(have_g_je) THEN
             CURVAV(I,JE)=CURVAV(I,JE-1)
             DLYP(I,JE)=DLYP(I,JE-1)
             DLYV(I,JE)=DLYV(I,JE-1)
             DLXV(I,JE)=DLXV(I,JE-1)
             DLXP(I,JE)=DLXP(I,JE-1)
             DLXU(I,JE)=DLXU(I,JE-1)
             DLYU(I,JE)=DLYU(I,JE-1)
          ENDIF
13788  CONTINUE

      CALL GATHER(DLXP,DLXP_G,p_io)
      CALL GATHER(DLYP,DLYP_G,p_io)
      CALL GATHER(DLXU,DLXU_G,p_io)
      CALL GATHER(DLYU,DLYU_G,p_io)
      CALL GATHER(DLXV,DLXV_G,p_io)
      CALL GATHER(DLYV,DLYV_G,p_io)

      CALL BELEG
!      PRINT*,'nach beleg'
!
      CALL BODEN
!      PRINT*,'nach boden'
!
!-------------------------------------------------------------------
!OtB LMONTS used in LEVIRE but not defined yet
      LMONTS=0
       IF (I3DREST .GT. 0) THEN
! READ 3D LEVITUS DATA FOR USE IN 3D RESTORING
!
!!$      IF(p_pe==p_io) THEN
!!$        REWIND(IO_IN_INIT)
!!$        REWIND(IO_IN_INIS)
!!$      ENDIF
      CALL LEVIRE(-1)
      ENDIF ! I3DREST
!
!-------------------------------------------------------------------
! SPECIFY CORIOLIS PARAMETER ETC.
!
      CALL CORIOL
!
!
! INITIALIZE DIAGNOSTICS

      CALL DIAG_INI   

!   write(0,*)'nach diag_ini'
!
!-------------------------------------------------------------------
! DIAGNOSTICS : MASKS OF BASINS
!               9 : GLOBAL
!
      IF(p_pe==p_io) OPEN(IO_IN_BGIN,FILE='BEK',FORM='FORMATTED')
! 
      DO J=1,JE_G
      DO I=1,IE_G
         IBEK_G(I,J)=9*NINT(WETO_G(I,J,1))
      ENDDO
      ENDDO
      JMMM=(JE_G-1)/120
      IF (ISTART .GT. 0 ) THEN
      IF(p_pe==p_io) THEN
        DO JMM=0,JMMM
        JMANF=1+JMM*120
        JMEND=MIN((JMM+1)*120,JE_G)
        DO I=2,IE_G-1
#ifdef VERSIONGR03
          READ(IO_IN_BGIN,'(I4,2X,120I1)')IJJU,                         &
     &                    (IBEK_G(I,J),J=JMEND,JMANF,-1)
#else
          READ(IO_IN_BGIN,'(I3,2X,120I1)')IJJU,                         &
     &                    (IBEK_G(I,J),J=JMEND,JMANF,-1)
#endif
        ENDDO
        ENDDO
      ENDIF
      CALL p_bcast(IBEK_G,p_io)
      ENDIF !ISTART
!
      DO J=1,JE_G
       IBEK_G(1,J)=IBEK_G(IE_G-1,J)
       IBEK_G(IE_G,J)=IBEK_G(2,J)
      ENDDO

      DO J=2,JE_G-1
       DO I=2,IE_G-1
        IBEK_G(I,J)=IBEK_G(I,J)*NINT(WETO_G(I,J,1))
        IF(WETO_G(I,J,1).GT.0.5.AND.IBEK_G(I,J).EQ.0)THEN
        IBEK_G(I,J)=MAX(IBEK_G(I+1,J),IBEK_G(I-1,J),IBEK_G(I,J+1),      &
     &                  IBEK_G(I,J-1))
        ENDIF
       ENDDO
      ENDDO

      IF(p_pe==p_io) THEN
        REWIND(IO_IN_BGIN)
        DO JMM=0,JMMM
        JMANF=1+JMM*120
        JMEND=MIN((JMM+1)*120,JE_G)
        WRITE(IO_STDOUT,*)'IBEK, JM ',JMM,JMEND,JMANF
        DO I=2,IE_G-1
#ifdef VERSIONGR03
          WRITE(IO_IN_BGIN,'(I4,2X,120I1)')I,                           &
     &                     (IBEK_G(I,J),J=JMEND,JMANF,-1)
#else
          WRITE(IO_IN_BGIN,'(I3,2X,120I1)')I,                           &
     &                     (IBEK_G(I,J),J=JMEND,JMANF,-1)
#endif

        ENDDO
        ENDDO
        CLOSE(IO_IN_BGIN)
      ENDIF

      IBEK(:,:) = IBEK_G(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)

!-------------------------------------------------------------------

      WRITE(IO_STDOUT,*) 'DZ ', DZ,DI,DZW,TIESTU
      WRITE(IO_STDOUT,*)'WETO:'
      DO JMM=0,JMMM
         JMANF=1+JMM*120
         JMEND=MIN((JMM+1)*120,JE_G)
         WRITE(IO_STDOUT,*)'JMM ',JMM,JMANF,JMEND
         WRITE(IO_STDOUT,6061)0,(MOD(J,10),J=JMEND,JMANF,-1)
         WRITE(IO_STDOUT,*)'        J  <=== '
!         DO I=1,IE_G
!            WRITE(IO_STDOUT,6061)I,(NINT(WETO_G(I,J,1))*MOD(J,10)       &
!                 -10*NINT(WETO_G(I,J,1)-1.),J=JMEND,JMANF,-1)
!         ENDDO
      ENDDO

      DO J=1,JE_G
         ICOU=0
         DO I=2,IE_G-1
            IF(WETO_G(I,J,1).GT.0.5)ICOU=ICOU+1
         ENDDO
            WRITE(IO_STDOUT,*)'ZAEHL : ',J,ICOU
         ENDDO
 6061       FORMAT(I4,1X,120I1)
!            Print*,'nach Zaehl'
!-------------------------------------------------------------------

      IF (ISTART .LT. 3) THEN
         DO K=1,KE
            DO J=1,JE
               DO I=1,IE
                  THO(I,J,K)=TLEVI(I,J,K)
                  SAO(I,J,K)=SLEVI(I,J,K)
               ENDDO
            ENDDO
         ENDDO
     ENDIF !ISTART








     IF(IOCAD.EQ.4)THEN
! PREPARE FOR BOTTOM BOUNDARY PARAMETRIZATIONS
        CALL FINDBOT(KBOT,WETO,IE,JE,KE)
        CALL FINDALFA
!      PRINT*,'nach SLOPECON_ADPO'
     ENDIF

!-------------------------------------------------------------------
!  provide matrix for barotropic mode
!     uwe  add iteration of matrix
!
      IF(ielimi.GE.1) THEN
         CALL trian
      ELSE
         CALL itprep
         CALL trotest2
      ENDIF
!      print*,'nach trotest'
!-------------------------------------------------------------------


    call river_runoff_ini
!      write(0,*)'no river_runoff_ini'

! INCLUDE WRITEOUT OF CPP OPTIONS AND PARAMETER SETTINGS
      CALL CPPOUTPUT
!-----------------------------------------------------------------------
!     
!     TIMESTEPPING IS BASED ON 3 LOOPS
!          OUTER LOOP : DO 1000  LYEAR = LYEAR1, LYEAR2
!        1.INNER LOOP : DO 1100  LMONT = 1 OR LMONT1, LMONT2 OR 12
!        2.INNER LOOP : DO 1010  LDAY  = 1,MONLEN(LMONTS)
!        3.INNER LOOP : DO 1001  LDTDAY  = 1,NDTDAY
! ADDITIONAL COUNTERS
!        LYEARS = ACTUAL YEAR
!        LMONTS = ACTUAL MONTH
!        LDAYS  = ACTUAL DAY
!     
!     TIME COUNTER ARE :
!               LDT       : TIME STEP COUNTER FOR EXPERIMENT
!               LDTRUN    : TIME STEP COUNTER FOR RUN
!               LDTYEAR   : TIME STEP COUNTER FOR YEAR
!               LDTMONTH  : TIME STEP COUNTER FOR MONTH
!               NDTMONTH  : TIME STEPS IN ACTUAL MONTH
!               NDTDAY  : TIME STEPS PER DAY
      WRITE(IO_STDOUT,*)'START FROM PREVIOUS CALCULATION YES=1/NO=0 : ' &
     &                  ,IAUFR
      WRITE(IO_STDOUT,*)'WRITE RESTART FILE              YES=1/NO=0 : ' &
     &                  ,IAUFW




!
!
!-----------------------------------------------------------------------
!
!     START FROM STATUS LEVITUS OR HORIZONTAL STRATIFICATION
! 
!  CALL LEVIRE WITH ARGUMENT =  0 : HORIZONTAL STRATIFICATION
!  CALL LEVIRE WITH ARGUMENT = 13 : 3D STRATIFICATION
!  CALL LEVIRE WITH ARGUMENT = -2 : SURFACE SALINITY



      CALL LEVIRE(-2)
      IF (RELTEM .GT. ALMZER) THEN
      CALL LEVIRE(-3)
      ENDIF


      IF(IAUFR.EQ.0) THEN
         IF (ISTART .EQ. 0) CALL LEVIRE(13)
         IF (ISTART .EQ. 1) CALL LEVIRE(0)
         IF (ISTART .EQ. 2) THEN
         
	 CALL LEVIRE(13)
         DO J=1,JE
          DO I=1,IE
           IF (THO(I,J,1)*WETO(I,J,1).LT.-0.5) THEN
             SICTHO(I,J)=MAX(0.,3.*(-0.5-THO(I,J,1)))
             SICOMO(I,J)=1.
           ENDIF
          ENDDO
         ENDDO
       ENDIF ! ISTART

!  CHECK GLOBAL SALT CONTENT
 
       IF ( ICONTRO .ne. 0 ) THEN
          CALL CONTRO(24)
       ENDIF
 
!  SET AGE OF LEVITUS DATA
 
         LYEARS=0
         LMONTS=0 
         LDAYS=0 
         LDT = 0
         WRITE(IO_STDOUT,*)'TIME COUNTER SET TO ZERO ' 
 
!  SET LOGICAL UNIT FOR RESTART FILE; FOR FIRST RUN WRITE TO IO_IN_Z370
 
         IUNIT = IO_IN_Z380
         IFLAG = -1
!     WRITE RESTART FILES WITH INITIAL STRATIFICATION                   
!     MAKE SURE THAT BOTH RESTART FILES EXIST ON EXIT
         
	CALL AUFW
         CALL AUFW

      ENDIF

!
!-----------------------------------------------------------------------

!     START FROM RESTART FILES Z37000 OR Z38000

      IF (IAUFR .EQ. 1) THEN

         CALL AUFR

 
!  CHECK GLOBAL SALT CONTENT
        IF ( ICONTRO .ne. 0 ) THEN
         CALL CONTRO(13)
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
! SET LIMITS OF YEAR/MONTH TIME STEPPING LOOP

      LDTRUN = 0

      LYEAR1 = LYEARS

      LMONT1 = LMONTS + 1

      IF (LMONT1 .GT. 12) THEN
         LMONT1=1
         LYEAR1=LYEAR1+1
      ENDIF

      LMONT2 = LMONT1 + MAX(0,12*NYEARS + NMONTS - 1)

      LYEAR2 = LYEAR1 + (LMONT2-1)/12

      LMONT2 = MOD(LMONT2-1,12)+1

#ifdef AMOCEMR
!-------------
       CALL AMOCPR
!--------------
#endif /*AMOCEMR*/


#ifdef PBGC
!     
!--------------------------------------------------------------------
! Initialize bgc modules
!  
!
!  Length of run in model time steps:
      DO  lyear=lyear1,lyear2
#ifdef RYEAR
#ifndef __coupled
        nacyear=lyear
#endif /*__coupled*/
        CALL MONLEN_FEB(nacyear,len_feb)
        monlen(2) = len_feb
#endif /*RYEAR*/
       lmon1=1
       lmon2=12
       IF(lyear1.EQ.lyear2) THEN
          lmon1=lmont1
          lmon2=lmont2
       ELSE
          IF(lyears.EQ.lyear2) THEN
             lmon1=1
             lmon2=lmont2
          ENDIF          
          IF(lyears.EQ.lyear1) THEN
             lmon1=lmont1
             lmon2=12
          ENDIF
       ENDIF
      ndtrun = 0
      DO  lmont=lmon1,lmon2
      DO  lday=1,monlen(lmont)
         ndtrun = ndtrun + NINT(86400./dt)
      ENDDO
      ENDDO
      ENDDO
      WRITE(IO_STDOUT,*)'Total no. of time steps this run: ',ndtrun

      CALL INI_BGC(iaufr,icycli,dt,ndtrun,ie,je,ke                    &
     &                  ,ddpo,tho,sao,dlxp,dlyp,tiestu,tiestw         &
     &                  ,lyears,lmonts,ldays,ldt,nyears,nmonts        &
     &                  ,gila,giph)

      IF (iocad == 5) THEN
        WRITE(io_stdout,*) ' '
        WRITE(io_stdout,*) ' WARNING:'
        WRITE(io_stdout,*) ' iocad == 5 i.e. advection routine OCADFS', &
                           ' is chosen but BGC Model HAMOCC is active!'
        WRITE(io_stdout,*) ' There is NO ADVECTION of BGC tracers'
      ENDIF
        


#endif /*PBGC*/

WRITE(0,*) 'Open Addcontra'
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifdef ADDCONTRA
    
!--------------------------------------------------------------------
! Initialize add modules
!  
!
!  Length of run in model time steps:
      DO  lyear=lyear1,lyear2

#ifdef RYEAR
#ifndef __coupled
        nacyear=lyear
#endif /*__coupled*/
        CALL MONLEN_FEB(nacyear,len_feb)
        monlen(2) = len_feb
#endif /*RYEAR*/

       lmon1=1
       lmon2=12
       IF(lyear1.EQ.lyear2) THEN
          lmon1=lmont1
          lmon2=lmont2
       ELSE
          IF(lyears.EQ.lyear2) THEN
             lmon1=1
             lmon2=lmont2
          ENDIF          
          IF(lyears.EQ.lyear1) THEN
             lmon1=lmont1
             lmon2=12
          ENDIF
       ENDIF
      add_ndtrun = 0
      DO  lmont=lmon1,lmon2
      DO  lday=1,monlen(lmont)
         add_ndtrun = add_ndtrun + NINT(86400./dt)
      ENDDO
      ENDDO
      ENDDO
      WRITE(IO_STDOUT,*)'Total no. of time steps this run: ',add_ndtrun

      CALL INI_ADD(iaufr,icycli,dt,add_ndtrun,ie,je,ke                    &
     &                  ,ddpo,tho,sao,dlxp,dlyp,tiestu,tiestw         &
     &                  ,lyears,lmonts,ldays,ldt,nyears,nmonts        &
     &                  ,gila,giph)

      IF (iocad == 5) THEN
        WRITE(io_stdout,*) ' '
        WRITE(io_stdout,*) ' WARNING:'
        WRITE(io_stdout,*) ' iocad == 5 i.e. advection routine OCADFS', &
                           ' is chosen but ADDCONTRA Model is active!'
        WRITE(io_stdout,*) ' There is NO ADVECTION of ADD tracers'
      ENDIF
        


#endif /*ADDCONTRA*/

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifdef FB_BGC_OCE
#ifndef PBGC
      WRITE(IO_STDOUT,*)''
      WRITE(IO_STDOUT,*)'SW radiation via attenuation map activ!'
      WRITE(IO_STDOUT,*)''

      CALL ALLOC_MEM_ATTEN(ie,je,ke)
#endif
#endif
!--------------------------------------------------------------------
#if defined (__coupled)

!  Initialise the coupled model : 
!     nactyear = first ocean year of run (from restart data) at entry
!     nactyear = first year of actual date (from KONTROL.h)  at exit
!  
      nactyear = lyear1
      CALL couple_init(nactyear)
#endif /*(__coupled)*/
!--------------------------------------------------------------------

#ifdef RYEAR
      WRITE(IO_STDOUT,*)                                                &
     &'THIS RUN ASSUMES REALISTIC YEAR WITH LEAP FEBRUARIES'
#else

#ifdef YEAR360
      WRITE(IO_STDOUT,*)                                                &
     &'THIS RUN ASSUMES IDEALIZED YEARS (12*30 DAYS EACH)'
#else
      WRITE(IO_STDOUT,*)                                                &
     &'THIS RUN ASSUMES IDEALIZED YEARS WITH 365 DAYS'
#endif /*YEAR360*/

#endif /*RYEAR*/

      WRITE(IO_STDOUT,*)                                                &
     &' INTEGRATION PERIOD IS YEAR ',LYEAR1,' TO ',LYEAR2

#ifdef TIDAL
!    allocate memory for tidal model 
      CALL alloc_mem_tidal
      CALL foreph_ini
#endif


      DO 1000 LYEAR=LYEAR1,LYEAR2
         LYEARS=LYEAR
         LDTYEAR = 0 

         WRITE(IO_STDOUT,*) ' TO BE CALCULATED NOW : YEAR = ',LYEARS

#ifdef RYEAR
#ifndef __coupled
         NACTYEAR=LYEAR
#endif /*__coupled*/
         CALL MONLEN_FEB(NACTYEAR,LEN_FEB)
         MONLEN(2) = LEN_FEB
         WRITE(IO_STDOUT,*)                                             &
     & 'FEBRUARY OF YEAR ',NACTYEAR,' HAS ',MONLEN(2),' DAYS'
#endif /*RYEAR*/

#ifndef __coupled


#ifdef CORE
       WRITE(IO_STDOUT,*) 'call spool core'
       CALL spool_core
#else
       WRITE(IO_STDOUT,*) 'call spool omip'
       CALL spool_omip
#endif

!-------------------------------------------------------------------
!-------------------------------------------------------------------
#ifdef ADDCONTRA
 
       WRITE(IO_STDOUT,*) 'call spool omip_add'
       CALL spool_omip_add
#endif 
!-------------------------------------------------------------------
!-------------------------------------------------------------------

#endif /*(__coupled)*/


!--------------------------------------------------------------------
!  SET NUMBER OF TIME STEPS PER YEAR

       ndtyear=0
       DO i=1,12
         ndtyear=ndtyear+monlen(i)*ndtday
       ENDDO


!  SET MONTH LOOP LIMITS FOR EACH YEAR

       LMON1=1
       LMON2=12
       IF(LYEAR1.EQ.LYEAR2) THEN
          LMON1=LMONT1
          LMON2=LMONT2
       ELSE
          IF(LYEARS.EQ.LYEAR2) THEN
             LMON1=1
             LMON2=LMONT2
          ENDIF          
          IF(LYEARS.EQ.LYEAR1) THEN
             LMON1=LMONT1
             LMON2=12
          ENDIF
       ENDIF
          WRITE(IO_STDOUT,*) 'MONTHLY LOOP FROM MONTH ',LMON1,          &
     &                        ' TO',LMON2

!:: MONTH LOOP
       DO 1100 LMONT=LMON1,LMON2
          LMONTS=LMONT
          LDTMONTH = 0

          NDTMONTH = MONLEN(LMONTS) * NDTDAY

          WRITE(IO_STDOUT,18823) LYEARS,LMONTS
18823     FORMAT(' TO BE CALCULATED NOW : YEAR ',I4,' MONTH ',I4)


#ifdef FB_BGC_OCE
#ifndef PBGC
      WRITE(IO_STDOUT,*)''
      WRITE(IO_STDOUT,*)'Read attenuation map at month: ',lmonts
      WRITE(IO_STDOUT,*)''

      CALL READ_ATTEN(ie,je,ke,ddpo,lmonts)
#endif
#endif


#ifdef RESTORE_MON
!:: REWIND SURSAL TO BEGIN OF YEAR
       IF(p_pe==p_io) THEN
         REWIND(IO_IN_SURS)
         IF (RELTEM .GT. ALMZER) REWIND(IO_IN_SURT)
       ENDIF
!:: READ APPROPRIATE SURSAL FOR MONTHLY RESTORING
       WRITE(IO_STDOUT,*)                                               &
     & 'READING MONTHLY SURFACE SALINITY IN MONTH ',LMONT
       DO IREADC=1,LMONT
        CALL LEVIRE(-2)
        IF (RELTEM .GT. ALMZER)                                         &
     &  CALL LEVIRE(-3)
       ENDDO
#endif /*RESTORE_MON*/
!::
       IF (I3DREST .GT. 1) THEN
!:: REWIND INISAL/INITEM TO BEGIN OF YEAR
       IF(p_pe==p_io) THEN
         REWIND(IO_IN_INIT)
         REWIND(IO_IN_INIS)
       ENDIF
!:: READ APPROPRIATE INISAL/INITEM FOR MONTHLY RESTORING
       WRITE(IO_STDOUT,*)                                               &
     & 'READING MONTHLY LEVITUS FIELDS IN MONTH ',LMONT
       DO IREADC=1,LMONT
        CALL LEVIRE(-1)
       ENDDO
       ENDIF ! I3DREST

!::================================================================

      if (NDAYS.gt.0.and.NDAYS.lt.MONLEN(LMONTS)) then
         ENDDAY=NDAYS
      else
         if (NDAYS.gt.MONLEN(LMONTS)) then
            write(0,*)'ATTN: ndays greater than monlen(lmonts)'
            write(0,*)'=> ndays is set to ',monlen(lmonts)
         endif
         ENDDAY=MONLEN(LMONTS)
      endif

!::==================================================================

      DO 1010 LDAY=1,ENDDAY
      LDAYS=LDAY
      DO 1001 LDTDAY = 1,NDTDAY
      LDTDAYC = LDTDAY

#if defined (TIMECHECK) || defined (__coupled)
      ttts = p_time()
#endif

!*       6.  TIMESTEP THE OCEAN.
!            -------- --- ------
      LDT = LDT + 1
      LDTRUN = LDTRUN + 1
      LDTYEAR = LDTYEAR + 1
      LDTMONTH =LDTMONTH + 1

!  current date: yyyy-mm-dd (W3C XML Schema)
      IF (p_pe==p_io) THEN
        WRITE(0,'(a,i4.4,''-'',i2.2,''-'',i2.2,a,i3,a,i7)') &
             'MPIOM: current date: ', lyear,lmonts,lday, &
             ' begin of timestep (day):', ldtday,'  (run):', ldtrun
      ENDIF


#if defined (__coupled)

! Get data from coupler

          CALL couple_get_a2o(ldtrun)
#else /*(__coupled)*/



          IF (MOD(LDTDAY-1,NDTDAY) .EQ. 0) THEN
          WRITE(IO_STDOUT,*)'READ SURFACE FORCING AT TIMESTEP ',LDTDAY  &
                            ,LDAYS,LMONTS,LYEARS

#ifdef CORE
       WRITE(IO_STDOUT,*) 'call read core'
       CALL read_core 
#else
       WRITE(IO_STDOUT,*) 'call read omip'
       CALL read_omip
#endif

!--------------------------------------------------------------------
!--------------------------------------------------------------------
#ifdef ADDCONTRA
       WRITE(IO_STDOUT,*) 'call read omip_add'
       CALL read_omip_add
#endif
!-------------------------------------------------------------------
!-------------------------------------------------------------------


    ENDIF




#endif /*(__coupled)*/

#ifdef bounds_exch_tp
      call vcheck(51,voe)
      if(mod(ldtrun,24).eq.1) then
      s1o=voe
      CALL bounds_exch('vd',s1o,'mpiom 10')
      voe(:,2,:)=0.5*(voe(:,2,:)+s1o(:,1,:))
     endif
#endif


#ifdef PBGC
#ifdef TIMECHECK
      tttr = p_time()
#endif
!--------------------------------------------------------------------
! Call biogeochemistry modules.
      ddpomin=100.
      ddpomax=1.
      DO J=1,JE
      DO I=1,IE
        bgcdpio(i,j,1)=0.
        bgcddpo(i,j,1)=0.
!        layer1_bgc(i,j)=0.
        IF (ddpo(i,j,1).GT.0.)THEN
          bgcddpo(i,j,1)= DDPO(I,J,1)                                  &
     &    +ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA        
          bgcdpio(i,j,1)=1/bgcddpo(i,j,1)
!          layer1_bgc(i,j)=bgcddpo(i,j,1)
        ENDIF
      ENDDO
      ENDDO    
      DO K=2,KE
      DO J=1,JE
      DO I=1,IE 
        bgcddpo(i,j,k)= DDPO(I,J,K) 
        bgcdpio(i,j,k)= DPIO(i,j,k)
      ENDDO
      ENDDO
      ENDDO
#ifdef __coupled
      CALL BGC(ie,je,ke,                                               &
     &         AOFLSHWO,sicomo,tho,sao,bgcddpo,dlxp,dlyp,              &
     &         tiestu,bgcdpio,AOFLWSVO,wo,pao,lyears,lmonts,ldays,     &
     &         monlen(lmonts),ldtmonth,ldtday)
#else      

      call bounds_exch('p',fswr,'mpiom 19')
      call bounds_exch('p',sicomo,'mpiom 20')
      call bounds_exch('p',tho,'mpiom 20')
      call bounds_exch('p',sao,'mpiom 21')
      call bounds_exch('p',bgcddpo,'mpiom 22')
      call bounds_exch('p',dlxp,'mpiom 23')      
      call bounds_exch('p',dlyp,'mpiom 24')
      call bounds_exch('p',fu10,'mpiom 25')
      call bounds_exch('p',wo,'mpiom 26')
      call bounds_exch('p',pao,'mpiom 27')
      call bounds_exch('p',bgcdpio,'mpiom 28')




      CALL BGC(ie,je,ke,                                               &
     &         fswr,sicomo,tho,sao,bgcddpo,dlxp,dlyp,                  &
     &         tiestu,bgcdpio,fu10,wo,pao,lyears,lmonts,ldays,         &
     &         monlen(lmonts),ldtmonth,ldtday)
#endif

#endif /*PBGC*/

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef ADDCONTRA

! WRITE(0,*) 'Start ADD'
      ddpomin=100.
      ddpomax=1.
      DO J=1,JE
      DO I=1,IE
        add_dpio(i,j,1)=0.
        add_ddpo(i,j,1)=0.
        IF (ddpo(i,j,1).GT.0.)THEN
          add_ddpo(i,j,1)= DDPO(I,J,1)                                  &
     &    +ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA        
          add_dpio(i,j,1)=1/add_ddpo(i,j,1)

        ENDIF
      ENDDO
      ENDDO    
      DO K=2,KE
      DO J=1,JE
      DO I=1,IE 
        add_ddpo(i,j,k)= DDPO(I,J,K) 
        add_dpio(i,j,k)= DPIO(i,j,k)
      ENDDO
      ENDDO
      ENDDO  

     
      call bounds_exch('p',add_ddpo,'mpiom 22')
      call bounds_exch('p',dlxp,'mpiom 23')      
      call bounds_exch('p',dlyp,'mpiom 24')
      
      call bounds_exch('p',add_dpio,'mpiom 28')
!     WRITE(0,*) 'Max of add_ddpo', MAXVAL(add_ddpo(:,:,:))

 CALL ADD(ie,je,ke,add_ddpo,dlxp,dlyp,SICTHO,SICOMO,        &
     &         tiestu,add_dpio,lyears,lmonts,ldays,         &
     &         monlen(lmonts),ldtmonth,ldtday)
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#ifdef TIMECHECK
      ttt = p_time()
      WRITE(io_stdout,123)'Time BGC',ttt-tttr
#endif

!--------------------------------------------------------------------
! OCTHER : THERMODYNAMIC FORCING

#ifdef TIMECHECK
      tttr = p_time()
#endif

#ifdef TIDAL
      CALL foreph
#endif /*tidal*/

      CALL OCICE

!      goto 999
      CALL OCTHER

!      write(0,*) 'without octher !!!'
       

#ifdef TIMECHECK
      ttt = p_time()
      WRITE(io_stdout,123)'Time OCTHER',ttt-tttr
#endif
  123 FORMAT(a,f10.3)
!--------------------------------------------------------------------


#ifdef CORE
      WRITE(IO_STDOUT,*) 'call norm pem'
      CALL NORMPEM 
#endif



      IF (ibolk .LT. 0) THEN
! update gm coefficients once/day !
         IF (ldtday .EQ. 1) CALL calcgmvis
      ENDIF


       IF ( ICONTRO .ne. 0 ) THEN
          CALL CONTRO(1)
       ENDIF
!--------------------------------------------------------------------
! OCWIND : 

      CALL OCWIND
!      PRINT*,'nach ocwind'

!--------------------------------------------------------------------
#ifdef TIDAL
!  TIDES
      CALL tipouv
#endif /*tidal*/

!--------------------------------------------------------------------
!  OCTIMF : UPDATE VELOCITY
!     
      CALL OCTIMF
!--------------------------------------------------------------------
!  OCMODMOM : DECOMPOSITION INTO BAROTROPIC AND BAROCLINIC FIELD   

#ifdef TIMECHECK
      tttr = p_time()
#endif
      CALL OCMODMOM

#ifdef TIMECHECK
      ttt = p_time()
      WRITE(io_stdout,123)'Time OCMODMOM',ttt-tttr
#endif

!--------------------------------------------------------------------

      CALL OCBARP

!      PRINT*,'nach ocbarp'
!
!--------------------------------------------------------------------
   
#ifdef TIMECHECK
      tttr = p_time()
#endif

      CALL OCCLIT

!      PRINT*,'nach occlit'
#ifdef TIMECHECK
      ttt = p_time()
      WRITE(io_stdout,123)'Time OCCLIT',ttt-tttr
#endif

!--------------------------------------------------------------------
     
#ifdef TIMECHECK
      tttr = p_time()
#endif
       IF ( ICONTRO .ne. 0 ) THEN
          call contro(38)
       endif

      IF(IELIMI.GE.1) THEN
       CALL BARTIM
      ELSE
#ifdef SORTEST
        z1o(:,:)=0
#endif
        CALL TRONEU
#ifdef SORTEST
        allocate( z1(ie,je),z2(ie,je))
        z1(:,:)=z1o(:,:)
        z1o(:,:)=0.

        CALL TRONEU2

        z2(:,:)=z1o(:,:)

        z1o(:,:)= z1(:,:)   !use z1o from troneu

        z1(:,:)=z1(:,:)-z2(:,:)
        write(0,*)'z1',p_pe,maxval(z1),maxloc(z1)

        call p_barrier

        deallocate(z1,z2)
!        stop
#endif

      ENDIF
       IF ( ICONTRO .ne. 0 ) THEN
       call contro(39)
    endif
!      PRINT*,'nach troneu'
#ifdef TIMECHECK
      ttt = p_time()
      WRITE(io_stdout,123)'Time Solver',ttt-tttr
#endif
60023 FORMAT(4E20.12)
!
!--------------------------------------------------------------------
      CALL OCVTRO
!      PRINT*,'nach ocvtro'
      CALL update_zo

!      call maschk(ie,je,ke,51)
!
!--------------------------------------------------------------------
!     
#ifdef TIMECHECK
      tttr = p_time()
#endif
       IF ( ICONTRO .ne. 0 ) THEN
       call contro(40)
    endif
      CALL OCVTOT

!      PRINT*,'nach ocvtot'
       IF ( ICONTRO .ne. 0 ) THEN
       call contro(41)
       endif
#ifdef TIMECHECK
      ttt = p_time()
      WRITE(io_stdout,123)'Time OCVTOT',ttt-tttr
#endif
      if(imean.ne.0)then
         CALL WRTE_MFL(LDTDAY,LDAYS,LMONTS,LYEARS,IMEAN,NANF)
      endif
! 
!UWE   RESET STABIO TO DRHODZ
       DO K=1,KE
        DO J=1,JE
         DO I=1,IE 
          STABIO(I,J,K)=(1000./DZ(K))*STABIO(I,J,K)
         ENDDO
        ENDDO
       ENDDO
       IF ( ICONTRO .ne. 0 ) THEN
       call contro(42)
    endif


#ifndef bounds_exch_tp
       CALL calc_psi
       CALL calc_mixedlayerdepth
#endif
!      CALL maschk(ie,je,ke,42)       


       IF(itsdiag.GE.4) THEN
          CALL diagnosis
       ELSE
         IF (mod(ldtday,ndtday).EQ.ndtday/2) THEN
            CALL diagnosis
         ENDIF
       ENDIF

!       call maschk(ie,je,ke,43)

       call calc_moc

#ifdef TIMECHECK
      tttr = p_time()
#endif
       IF ( ICONTRO .ne. 0 ) THEN
          call contro(43)
       endif
       CALL OCUAD(UOO)
       CALL OCVAD(VOE)
       IF ( ICONTRO .ne. 0 ) THEN
          call contro(44)
       endif
!!        write(0,*)'without ocuvad'

#ifdef TIMECHECK
      ttt = p_time()
      WRITE(io_stdout,123)'Time OCUAD/OCVAD',ttt-tttr
#endif

      IF(IOCAD.EQ.4)THEN               ! SLOPECON_ADPO
         CALL SLOPETRANS
         IF ( ICONTRO .ne. 0 ) THEN
            CALL CONTRO(2)
         endif
      ENDIF

#ifdef TIMECHECK
      tttr = p_time()
#endif


       IF(IOCAD.EQ.5) CALL OCADFS
#ifdef TIMECHECK
      ttt = p_time()
      WRITE(io_stdout,123)'Time OCADFS',ttt-tttr
#endif
 

!--------------------------------------------------------------------
!  ADVECTION | EDDY-DIFFUSION | DIFFUSION 
!--------------------------------------------------------------------

#ifdef TIMECHECK
      tttr = p_time()
#endif

#ifdef PBGC

!  BGC tracer advection 
!
!!$! Dilution of bgc tracer due to fluxes of water.
!!$! BRINE contains ice volume change  (see SBR GROWTH).
!!$! PRECO contains atmospheric fluxes (see SBR GROWTH).
!!$! Tracer concentration is assumed to be zero in water added.
!!$! SICDIO contains old mixed layer depths before update of 
!!$! ZO, SICTHO, SICSNO in SBR GROWTH.
!!$
!!$        DO J=1,JE
!!$        DO I=1,IE
!!$          layer1_new(i,j)=0.
!!$          IF (ddpo(i,j,1).GT.0.)THEN
!!$            layer1_new(i,j)=DDPO(I,J,1)+ZO(I,J)                    &
!!$     &       -WO(I,J,1)*DT-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA
!!$          ENDIF
!!$        ENDDO
!!$        ENDDO        
!!$        
!!$      CALL DILUTE_BGC(IE,JE,layer1_bgc,layer1_new)


!     IF(IOCAD.EQ.4)THEN
!!$       call maschk(ie,je,ke,44)
!!$ !$OMP PARALLEL
!!$ !$OMP DO
!!$         DO l=1,nocetra
!!$            DO K=1,KE
!!$               DO J=1,JE
!!$                  DO I=1,IE
!!$                     OCETRA(I,J,K,L)=MAX(0.,OCETRA(I,J,K,L))
!!$                  ENDDO
!!$               ENDDO
!!$            ENDDO
!!$ !$OMP ENDDO
!!$ !$OMP END PARALLEL  
!!$         ENDDO
!     ENDIF !  iocad.eq.4
#endif /*PBGC*/             

       CALL bounds_exch('p',rhoo,'mpiom 29')

!  Preparation for advection
       IF (iocad .NE. 5) CALL ocadpo_base

!  Preparation for diffusion
       CALL octdiff_base

!  Preparation for eddy-diffusion
       IF ( IBOLK .NE. 0 ) CALL ocjitr_base 

#ifdef PBGC
!  Loop over tracers (S, T and advected BGC tracers)
       ntracerloop=ntraad+2
#else
!  Loop over tracers (S and T)
        ntracerloop=2
#endif /*PBGC*/ 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef ADDCONTRA
!  Loop over tracers (S, T and advected ADD tracers)
       ntracerloop=nctraad+2
#else
!  Loop over tracers (S and T)
        ntracerloop=2
#endif /*ADDCONTRA*/  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~             

#ifdef TRACER_OMP

!  OpenMP routine-level parallel (outside routines) tracer calculations
!   - advection, eddy-diffusion, and diffusion are calculated completely
!     parallel for each tracer
!   - use a maximum of ntracerloop OpenMP threads
!   - without   HAMOCC: maximum of 2 OpenMP threads (for S and T)
!   - including HAMOCC: maximum of ntraad+2 threads (COSMOS: 16)

!$OMP PARALLEL PRIVATE(L,LL)
!$OMP DO
#endif
        DO LL=1,ntracerloop
          L=LL-2
          IF (LL==1) THEN

!  Advection of Temperature
            IF (iocad .NE. 5) CALL ocadpo_trf(tho)
            CALL bounds_exch('p',tho,'mpiom 30')

!  GM eddy-diffusion of Temperature
            IF ( IBOLK .NE. 0 ) CALL ocjitr_trf(tho)

          ELSE IF (LL==2) THEN

!  Advection of Salt
            IF (iocad .NE. 5) CALL ocadpo_trf(sao)
            CALL bounds_exch('p',sao,'mpiom 31')

!  GM eddy-diffusion of Salt
            IF ( IBOLK .NE. 0 ) CALL ocjitr_trf(sao)
!           call maschk(ie,je,ke,47)

#ifdef PBGC
          ELSE IF (LL>2) THEN

! Set tracer values at dry cells to bottom cell value as done
! for temperature and salinity in SBR OCTHER.

            DO K=2,KE
               DO J=1,JE
                  DO I=1,IE
                     IF(WETO(I,J,K).LT.0.5) THEN
                        OCETRA(I,J,K,L)=OCETRA(I,J,K-1,L)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO

!  BGC tracer advection
            IF (iocad .NE. 5) CALL ocadpo_trf(ocetra(1,1,1,L))

! Reset tracer values at dry cells to undefined.

#ifndef TRACER_OMP
!$OMP PARALLEL
!$OMP DO
#endif           
            DO K=1,KE
               DO J=1,JE
                  DO I=1,IE
                     IF(WETO(I,J,K).LT.0.5) THEN
                        OCETRA(I,J,K,L)=RMASKO
!                     ELSE      
!                        OCETRA(I,J,K,L)=MAX(0.,OCETRA(I,J,K,L)) 
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
#ifndef TRACER_OMP
!$OMP END DO
!$OMP END PARALLEL
#endif           

!  GM eddy-diffusion of BGC-tracer
            IF ( IBOLK .NE. 0 ) CALL ocjitr_trf(ocetra(1,1,1,L))

#endif /*PBGC*/             
!-----------------------------------------------------------------
#ifdef ADDCONTRA
          ELSE IF (LL>2) THEN

! Set tracer values at dry cells to bottom cell value as done
! for temperature and salinity in SBR OCTHER.

            DO K=2,KE
               DO J=1,JE
                  DO I=1,IE
                     IF(WETO(I,J,K).LT.0.5) THEN
                        OCECTRA(I,J,K,L)=OCECTRA(I,J,K-1,L)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO

!  ADD tracer advection
            IF (iocad .NE. 5) CALL ocadpo_trf(ocectra(1,1,1,L))

! Reset tracer values at dry cells to undefined.

#ifndef TRACER_OMP
!$OMP PARALLEL
!$OMP DO
#endif           
            DO K=1,KE
               DO J=1,JE
                  DO I=1,IE
                     IF(WETO(I,J,K).LT.0.5) THEN
                        OCECTRA(I,J,K,L)=RMASKO
 !                     ELSE      
 !                       OCECTRA(I,J,K,L)=MAX(0.,OCECTRA(I,J,K,L)) 
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
#ifndef TRACER_OMP
!$OMP END DO
!$OMP END PARALLEL
#endif           

!  GM eddy-diffusion of ADD-tracer
            IF ( IBOLK .NE. 0 ) CALL ocjitr_trf(ocectra(1,1,1,L))

#endif /*ADDCONTRA*/                        
!---------------------------------------------------------------------------
          ENDIF !  LL
        ENDDO ! LL
#ifdef TRACER_OMP
!$OMP ENDDO
!$OMP END PARALLEL  
#endif           

       if (LCONVDIAG) call calc_potential_energy_release(1)

#ifdef TRACER_OMP
!$OMP PARALLEL PRIVATE(L,LL)
!$OMP DO
#endif
        DO LL=1,ntracerloop
          L=LL-2
          IF (LL==1) THEN

!  Diffusion of Temperature
            CALL octdiff_trf(tho)

          ELSE IF (LL==2) THEN

!  Diffusion of Salt
            CALL octdiff_trf(sao)

#ifdef PBGC
          ELSE IF (LL>2) THEN

!  BGC tracer diffusion
            CALL octdiff_trf(ocetra(1,1,1,L))
#endif /*PBGC*/             
!-----------------------------------------------------------
#ifdef ADDCONTRA
          ELSE IF (LL>2) THEN

!  ADD tracer diffusion
            CALL octdiff_trf(ocectra(1,1,1,L))
#endif /*ADDCONTRA*/
!-----------------------------------------------------------  
          ENDIF !  LL
        ENDDO ! LL
#ifdef TRACER_OMP
!$OMP ENDDO
!$OMP END PARALLEL  
#endif           

       if (LCONVDIAG) call calc_potential_energy_release(3)

!           call maschk(ie,je,ke,48)
!  
       CALL RELAX_TS
#ifdef NUDGE_ECCO

       CALL READ_ECCO
       
       CALL NUDGE_T(tho,rreltem)

       CALL NUDGE_S(sao,rrelsal)
#endif

       IF ( ICONTRO .ne. 0 ) THEN
          CALL CONTRO(4)
       endif

       CALL OCTIMF
       IF ( ICONTRO .ne. 0 ) THEN
          CALL CONTRO(5)
       endif

#ifdef TIMECHECK
          ttt = p_time()
          WRITE(io_stdout,123)'Time ADVECTION/DIFFUSION/GM',ttt-tttr
#endif

!--------------------------------------------------------------------
!  MOMENTUM DIFFUSION
!
#ifdef TIMECHECK
       tttr = p_time()
#endif
       IF (AUS .GT. 1.E-12) CALL OCSCHEP
#ifdef TIMECHECK
       ttt = p_time()
       WRITE(io_stdout,123)'Time OCSCHEP',ttt-tttr
#endif

!--------------------------------------------------------------------
!  BIHARMONIC MOMENTUM DIFFUSION, BOTTOM FRICTION, 
!  SHEAR DEPENDENT DIFFUSION

#ifdef TIMECHECK
       tttr = p_time()
#endif

       CALL OCVISC

#ifdef TIMECHECK
       ttt = p_time()
       WRITE(io_stdout,123)'Time OCVISC',ttt-tttr
#endif
       IF ( ICONTRO .ne. 0 ) THEN
          CALL CONTRO(6)
       endif

       CALL OCTIMF

#ifdef FLUXCORRECT
!JJ     MONTHLY AVERAGE AND YEARLY AVERAGE FOR FLUX CORRECTION
        IF (LDAYS.EQ.1) THEN
          NANF=LDTDAY
        ELSEIF (LDAYS.GT.1) THEN
          NANF=LDTDAY+((LDAYS-1)*NDTDAY)
        ENDIF
        NEND=NDTDAY*MONLEN(LMONTS)
! CALCULATE END OF YEAR FOR 12 MONTH
        MYEND=0
        DO I=1,12
          MYEND=MYEND+MONLEN(I)
        ENDDO
!:: CALCULATE ACTUAL DAY IN YEAR
        MYACT=0
        DO I=1,LMONTS-1
          MYACT=MYACT+MONLEN(I)
        ENDDO
        MYACT=MYACT+LDAYS-1
        MYACT=MYACT*NDTDAY+LDTDAY
        MYEND=MYEND*NDTDAY
#ifdef FLUXCORRECT
        CALL FLUXC_MMEAN2D(LDAYS,LMONTS,LYEARS,                         &
     &             NANF,NEND,MYACT,MYEND)
#endif /*FLUXCORRECT*/
#endif /*FLUXCORRECT*/

       if (iMEAN.ne.0)then
          CALL WRTE_MEAN(LDTDAY,LDAYS,LMONTS,LYEARS,IMEAN,NANF)
       endif
#ifdef PBGC
       CALL WRTE_MEANBGC
#endif

!-----------------------------------------------------------
#ifdef ADDCONTRA
       CALL WRTE_MEANADD
      
#endif
!-----------------------------------------------------------

 999   CONTINUE

#if defined (__coupled)

! Put data to coupler

          CALL couple_put_o2a(ldtrun)
#endif /*(__coupled)*/
!
#if defined (__coupled)

! Update coupled run calendar
 
     CALL couple_calendar(dt,io_stdout)
#endif

#if defined (TIMECHECK)
      ttt = p_time()
      IF (p_pe==p_io) THEN
        WRITE(io_stdout,'(a,f7.3,a)') &
     &       'Time for 1 timestep: ', ttt-ttts, ' s'
      ENDIF
      WRITE(io_stdout,123) 'Time for 1 timestep',ttt-ttts
      WRITE(0,123) 'Time for 1 timestep',ttt-ttts
#endif


 1001 CONTINUE

! END OF ONE DAY
!----------------------------------------------------------------------

#ifndef ZO_NOCORRECT

!     correct ZO to a global mean of zero once per day
      call correct_zo

#endif /*ZO_NOCORRECT*/


 1010 CONTINUE


!#ifdef PBGC
!      CALL AVRG_BGCMEAN(ie,je,ke)
!#endif /*PBGC*/  


#ifndef RESYEAR

!  WRITE RESTART FILE AT END OF MONTH

      IF(IAUFW.EQ.1)THEN
        CALL AUFW

#ifdef PBGC
        CALL AUFW_BGC(ie,je,ke,ddpo,gila,giph,tiestu                    &
     &                ,lyears,lmonts,ldays,ldt)
#endif /*PBGC*/

!--------------------------------------------------------------------------------
#ifdef ADDCONTRA
        CALL AUFW_ADD(ie,je,ke,ddpo,gila,giph,tiestu                    &
     &                ,lyears,lmonts,ldays,ldt)
#endif /*ADDCONTRA*/
!--------------------------------------------------------------------------------
        
      ELSE
        WRITE(IO_STDOUT,*)'STOPPED WITHOUT WRITING RESTART FILE,        &
     &                     IAUFW= ',IAUFW
      ENDIF
#endif

#ifdef AMLDDIAG
          CALL WRTE_AMLDDIAG(LYEARS,LMONTS,MONLEN(LMONTS))
#endif /*AMLDDIAG*/

          IF ( LCALCDIFI) THEN
             CALL CALC_DIFI((LYEARS*10000)+(LMONTS*100)+LDAYS)
          ENDIF

! END OF ONE MONTH
!----------------------------------------------------------------------
#ifdef FLUXCORRECT
      CALL FLUX_CORRECT
#endif /*FLUXCORRECT*/

 1100 CONTINUE

#ifdef RESYEAR
!  WRITE RESTART FILE AT END OF YEAR
      IF(IAUFW.EQ.1)THEN
        CALL AUFW
!        CALL AUFW_CDI
#ifdef PBGC
        CALL AUFW_BGC(ie,je,ke,ddpo,gila,giph,tiestu                    &
     &                ,lyears,lmonts,ldays,ldt)
#endif /* PBGC*/

!----------------------------------------------------------
#ifdef ADDCONTRA
        CALL AUFW_ADD(ie,je,ke,ddpo,gila,giph,tiestu                    &
     &                ,lyears,lmonts,ldays,ldt)
#endif /* ADDCONTRA*/
!----------------------------------------------------------

      ELSE
        WRITE(IO_STDOUT,*)'STOPPED WITHOUT WRITING RESTART FILE         &
     &        AT END OF YEAR,                                                 &
     &                     IAUFW= ',IAUFW
      ENDIF
#endif /*RESYEAR*/

! END OF ONE YEAR
!----------------------------------------------------------------------

1000  CONTINUE

! END OF TIME STEPPING
!----------------------------------------------------------------------

#ifdef KONVDIAG
       CALL WRTE_KONVDIAG(LYEAR2,LMONT2,MONLEN(LMONT2))
#endif /*KONVDIAG*/

if (LGRIDINFO) then
       CALL WRTE_GRIDINFO(LYEAR2,LMONT2,MONLEN(LMONT2))
endif

#ifdef FLUXCORRECT
        CALL FLUXC_WRTE
#endif /*FLUXCORRECT*/

#ifdef PBGC
!----------------------------------------------------------------------
! Finish cleanly with marine bgc
!    
      CALL END_BGC(ie,je,ke,ddpo,dlxp,dlyp,gila,giph,tiestu     &
     &                ,lyears,lmonts,ldays,ldt)
#endif /*PBGC*/

!--------------------------------------------------------------------------
#ifdef ADDCONTRA
!----------------------------------------------------------------------
! Finish cleanly with marine add
!    
      CALL END_ADD(ie,je,ke,ddpo,dlxp,dlyp,gila,giph,tiestu     &
     &                ,lyears,lmonts,ldays,ldt)
#endif /*ADDCONTRA*/
!---------------------------------------------------------------------------

#if defined (__coupled)

! Finish cleanly
 
      CALL couple_end
#endif

!     Branch target for finishing cleanly
99999 CONTINUE
      CALL print_stats

      WRITE(nerr,*) 'NORMAL END OF MPIOM'

#ifdef CLOCK
  IF (util_cputime(zutime, zstime) == -1) THEN
     WRITE(nerr,*)'Cannot determine used CPU time'
  ELSE
     zwtime = util_walltime()
     zrtime = (zutime+zstime)/zwtime

     IF(p_pe==p_io) THEN
       WRITE (nerr,'(a,f10.2,a)') ' Wallclock        : ', zwtime, ' s'
       WRITE (nerr,'(a,f10.2,a)') ' CPU-time (user)  : ', zutime, ' s'
       WRITE (nerr,'(a,f10.2,a)') ' CPU-time (system): ', zstime, ' s'
       WRITE (nerr,'(a,f10.2,a)') ' Ratio            : ', 100*zrtime, ' %'
     END IF
  END IF
#endif

#ifdef bounds_exch_put
      CALL close_put_window
#endif
      CALL p_stop
      END

#ifdef __IFC /* Intel Fortran Compiler */
      SUBROUTINE handler(isig,icode)

      WRITE(io_stdout,*) 'Handler: ',isig,icode

      END
#endif










