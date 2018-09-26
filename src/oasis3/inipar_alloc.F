      SUBROUTINE inipar_alloc
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 0 *
C               * -------------     ------- *
C               *****************************
C
C**** *inipar_alloc*  - Get main run parameters to allocate arrays 
C
C     Purpose:
C     -------
C     Reads out run parameters.
C
C**   Interface:
C     ---------
C       *CALL*  *inipar_alloc*
C
C     Input:
C     -----
C     None
C
C     Output:
C     ------
C     None
C
C     Workspace:
C     ---------
C     None
C
C     Externals:
C     ---------
C     parse
C
C     Reference:
C     ---------
C     See OASIS manual (1995) 
C
C     History:
C     -------
C     Version Programmer  Date      Description
C     ------- ----------  ----      -----------  
C       2.5   A. Caubel   02/04/02  created from inipar.f
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* -----------------Include files and use of modules---------------------------
C
      USE mod_kinds_oasis
      USE mod_parameter
      USE mod_parallel
      USE mod_string
      USE mod_analysis
      USE mod_anais
      USE mod_rainbow
      USE mod_extrapol
      USE mod_unitncdf
      USE mod_experiment
      USE mod_timestep
      USE mod_coast
#if defined use_comm_MPI1 || defined use_comm_MPI2 || defined use_comm_GSIP || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
      USE mod_clim
#endif
      USE mod_calendar
      USE mod_gauss
      USE mod_hardware
      USE mod_unit
      USE mod_label
      USE mod_printing
#ifdef use_netCDF
#include <netcdf.inc>
#endif
C
C* ---------------------------- Local declarations --------------------
C
      CHARACTER*80 clline, clline_aux, clvari
      CHARACTER*9 clword, clfield, clstring, clmod, clchan
      CHARACTER*3 clind
      CHARACTER*2 cldeb
      CHARACTER*1 clequa
      CHARACTER*8 clwork
      CHARACTER*8 clstrg
      CHARACTER*7 cl_bsend

      CHARACTER(len=8), DIMENSION(:), ALLOCATABLE :: cl_aux
      INTEGER (kind=ip_intwp_p) il_varid, il_len, il_err, il_maxanal 
      INTEGER (kind=ip_intwp_p) nlonbf_notnc, nlatbf_notnc, 
     $    nlonaf_notnc, nlataf_notnc
      INTEGER (kind=ip_intwp_p) iind, il_redu, ib, il_aux, il_auxbf,
     $     il_auxaf, istatus, il_id
C
C* ---------------------------- Poema verses --------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Get basic info for the simulation 
C        ---------------------------------
C
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) 
     $    '           ROUTINE inipar_alloc  -  Level 0'
      WRITE (UNIT = nulou,FMT = *) 
     $    '           *********************    *******'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' Initialization of main run 
     $    parameters'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' Reading input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '

C
C* Initialization
      cchan = '    '
      ig_direct_nfield = 0
      ig_nfield = 0
      lg_oasis_field = .true.
      lg_vector = .FALSE.
C* Initialize character keywords to locate appropriate input
C
      clfield  = ' $NFIELDS'
      clchan   = ' $CHANNEL'
      clstring = ' $STRINGS'
      clmod    = ' $NBMODEL'
C
C* Get number of models involved in this simulation
C
      REWIND nulin
 100  CONTINUE
      READ (UNIT = nulin,FMT = 1001,END = 110) clword
      IF (clword .NE. clmod) GO TO 100
      READ (UNIT = nulin,FMT = 1002) clline
      CALL parse (clline, clvari, 1, jpeighty, ilen)
      IF (ilen .LE. 0) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Nothing on input for $NBMODEL '
          WRITE (UNIT = nulou,FMT = *) ' Default value will be used '
          WRITE (UNIT = nulou,FMT = *) ' '
      ELSE IF (ilen .GT. 1) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Input variable length is incorrect'
          WRITE (UNIT = nulou,FMT = *) ' There are too many models '
          WRITE (UNIT = nulou,FMT = *) ' ilen = ', ilen  
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Check $NBMODEL variable spelling '
          WRITE (UNIT = nulou,FMT = *) ' Default value will be used '
      ELSE
          READ (clvari,FMT = 1003) ig_nmodel
      ENDIF
C
C* Print out the number of models
C
      CALL prtout
     $    ('The number of models for this run is nmodel =', ig_nmodel,
     $    1)
      CALL alloc_experiment
C
C --> Get the message passing technique we are using
C
      REWIND nulin
 120  CONTINUE
      READ (UNIT = nulin,FMT = 1001,END = 130) clword
      IF (clword .NE. clchan) GO TO 120
      READ (UNIT = nulin,FMT = 1002) clline
      CALL parse (clline, clvari, 1, jpeighty, ilen)
      IF (ilen .LE. 0) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Nothing on input for $CHANNEL '
          CALL HALTE('STOP in inipar')
      ELSE IF (ilen .GT. 0 .AND. ilen .NE. 4) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Input variable length is incorrect'
          WRITE (UNIT = nulou,FMT = *) ' ilen = ', ilen  
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Check $CHANNEL variable spelling '
          CALL HALTE('STOP in inipar')
      ELSE
          cchan = clvari(1:4)
C
          IF (cchan .EQ.'MPI1' .or. cchan .EQ. 'MPI2'
     $        .or. cchan .EQ. 'GSIP') THEN
              CALL parse (clline, clvari, 2, jpeighty, ilen)
#if defined use_comm_MPI1 || defined use_comm_MPI2 || defined use_comm_GSIP || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_GSIP && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
              IF (ilen .eq. 0) THEN
                  lg_bsend = .true.
              ELSE
                  cl_bsend = clvari(1:7)
                  IF (cl_bsend .eq. 'NOBSEND') THEN
                      lg_bsend = .false.
                  ELSE
                      lg_bsend = .true.
                  ENDIF
              ENDIF
#endif
C* Get one additional line for each model
              DO 186 jm = 1, ig_nmodel
                READ (UNIT = nulin,FMT = 1002) clline
C*              Get the total number of processors for the model
                CALL parse (clline, clvari, 1, jpeighty, ilen)
                READ (clvari,FMT = 1004) nbtotproc(jm)
                WRITE (UNIT = nulou,FMT = *) ' '
                WRITE (UNIT=nulou,FMT='(''The total number of proc'',
     $          '' for model'', I2, '' is'', I2)')
     $          jm, nbtotproc(jm)
                WRITE (UNIT = nulou,FMT = *) ' '
C
C*              Get the nbr of processors involved in the coupling for the model
                CALL parse (clline, clvari, 2, jpeighty, ilen)
                IF (ilen .LE. 0) THEN
                    WRITE (UNIT = nulou,FMT = *) '      ***WARNING***'
                    WRITE (UNIT = nulou,FMT = *) 
     $     'No input for number of processors involved in the coupling'
                    WRITE (UNIT = nulou,FMT = *) 'for model', jm
                    WRITE (UNIT = nulou,FMT = *)
     $                  'Total number of processors will be used'
                    nbcplproc(jm)=nbtotproc(jm)
                ELSE
                    READ (clvari,FMT = 1004) nbcplproc(jm)
                ENDIF
                WRITE (UNIT = nulou,FMT = *) ' '
                WRITE (UNIT = nulou,FMT ='(''The number of processors'',
     $          '' involved in the coupling for model'', I2, '' is'', 
     $          I2)') jm, nbcplproc(jm)
                WRITE (UNIT = nulou,FMT = *) ' '
C
 186          CONTINUE
     
          ELSE IF(cchan .EQ. 'PVM3') THEN
              WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
              WRITE (UNIT = nulou,FMT = *) 
     $            ' PVM3 no longer supported'
              WRITE (UNIT = nulou,FMT = *) 
     $            ' Please keep on using Oasis 2.4'
              CALL HALTE('STOP in inipar')
          ENDIF
C
C* Print out the message passing technique
C
          CALL prcout
     $    (' The message passing used in OASIS is cchan =', cchan, 1)
      ENDIF
C
C* Formats
C
 1001 FORMAT(A9)
 1002 FORMAT(A80)
 1003 FORMAT(I1)
 1004 FORMAT(I3)

C
C*    2. Get field information
C        ---------------------
C
C* Read total number of fields exchanged
C
      REWIND nulin
 200  CONTINUE
      READ (UNIT = nulin,FMT = 2001,END = 210) clword
      IF (clword .NE. clfield) GO TO 200
      READ (UNIT = nulin,FMT = 2002) clline
      CALL parse(clline, clvari, 1, jpeighty, ilen)
      IF (ilen .LE. 0) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Nothing on input for $NFIELDS '
          WRITE (UNIT = nulou,FMT = *) ' Default value will be used '
          WRITE (UNIT = nulou,FMT = *) ' '
      ELSE
          READ (clvari,FMT = 2003) ig_total_nfield
      ENDIF
C
C* Print out the total number of fields exchanged
C
      CALL prtout
     $    ('The number of exchanged fields is nfield =', 
     $     ig_total_nfield, 1)
C
C* Alloc field number array 
C
      ALLOCATE (ig_number_field(ig_total_nfield),stat=il_err)
      IF (il_err.NE.0) CALL prtout 
     $    ('Error: ig_number_field allocation of inipar_alloc',il_err,1)
      ig_number_field(:)=0
C
C* Alloc field status array (logical indicating if the field goes through 
C* Oasis or not)
C
      ALLOCATE (lg_state(ig_total_nfield), stat=il_err)
      IF (il_err.NE.0) CALL prtout 
     $    ('Error: lg_state allocation of inipar_alloc',il_err,1)
      lg_state(:)=.false.
C
C* Alloc status of all the fields
C 
      ALLOCATE (ig_total_state(ig_total_nfield), stat=il_err)
      IF (il_err.NE.0) CALL prtout 
     $    ('Error: ig_total_state allocation of inipar_alloc',il_err,1)
      ig_total_state(:)=0
C
C* Alloc input field name array 
C
      ALLOCATE (cg_output_field(ig_total_nfield), stat=il_err)
      IF (il_err.NE.0) CALL prtout 
     $    ('Error: cg_output_field allocation of inipar_alloc',il_err,1)
      cg_output_field(:)=' ' 
C
C* Alloc number of analyses array 
C
      ALLOCATE (ig_total_ntrans(ig_total_nfield),stat=il_err)
      IF (il_err.NE.0) CALL prtout 
     $    ('Error: ig_total_ntrans"allocation of inipar_alloc',il_err,1)
      ig_total_ntrans (:) = 0
C 
C* Alloc array of restart file names, input and output file names
C
      ALLOCATE (cg_restart_file(ig_total_nfield),stat=il_err)
      IF (il_err.NE.0) CALL prtout 
     $  ('Error: cg_restart_FILE allocation of inipar_alloc',il_err,1)
      cg_restart_file(:)=' '
      ALLOCATE (cg_input_file(ig_total_nfield), stat=il_err)
      IF (il_err.NE.0) CALL prtout 
     $ ('Error in "cg_input_file"allocation of inipar_alloc',il_err,1)
      cg_input_file(:)=' '
C 
C* Alloc array of source and target locator prefix
C
      ALLOCATE (cga_locatorbf(ig_total_nfield),stat=il_err)
      IF (il_err.NE.0) CALL prtout 
     $  ('Error: cga_locatorbf allocation of inipar_alloc',il_err,1)
      cga_locatorbf(:)=' '

      ALLOCATE (cga_locatoraf(ig_total_nfield),stat=il_err)
      IF (il_err.NE.0) CALL prtout 
     $  ('Error: cga_locatoraf allocation of inipar_alloc',il_err,1)
      cga_locatoraf(:)=' '  
C
C* Get the SSCS for all fields
C
      REWIND nulin
 220  CONTINUE
      READ (UNIT = nulin,FMT = 2001,END = 230) clword
      IF (clword .NE. clstring) GO TO 220
C
C* Loop on total number of fields (NoF)
C
      DO 240 jf = 1, ig_total_nfield
C
C* First line
C
        READ (UNIT = nulin,FMT = 2002) clline
C* Get output field symbolic name
        CALL parse(clline, clvari, 2, jpeighty, ilen)
        cg_output_field(jf) = clvari
C* Get total number of analysis
        CALL parse(clline, clvari, 5, jpeighty, ilen)
        READ (clvari,FMT = 2003) ig_total_ntrans(jf)
C* Get field STATUS for OUTPUT fields
        CALL parse(clline, clvari, 6, jpeighty, ilen)
        IF (clvari(1:6) .eq. 'OUTPUT') THEN
           ig_direct_nfield = ig_direct_nfield + 1
           lg_state(jf) = .false.
           ig_total_state(jf) = ip_output
        ELSE
C* Get field status (direct or through oasis) and the number  
C* of direct and indirect fields if not PIPE nor NONE
            IF (cchan .ne. 'PIPE' .and. cchan .ne. 'NONE') THEN
                CALL parse(clline, clvari, 7, jpeighty, ilen)
                IF (clvari(1:8).eq.'EXPORTED') THEN
                    ig_nfield = ig_nfield + 1
                    lg_state(jf) = .true.
                    ig_number_field(jf) = ig_nfield
                    ig_total_state(jf) = ip_exported
C* Get restart file name
                    CALL parse(clline, clvari, 6, jpeighty, ilen)
                    cg_restart_file(jf) = clvari
                ELSEIF (clvari(1:7) .eq. 'IGNORED' ) THEN
                    IF (cchan .eq. 'GSIP') THEN
                        WRITE (UNIT = nulou,FMT = *)
     $                      'Direct exchange of fields between models'
                        WRITE (UNIT = nulou,FMT = *)
     $                      'not supported with GSIP'
                        CALL HALTE('STOP in inipar_alloc (GSIP)')
                    ENDIF
                    ig_direct_nfield = ig_direct_nfield + 1
                    lg_state(jf) = .false.
                    ig_total_state(jf) = ip_ignored
                    CALL parse(clline, clvari, 6, jpeighty, ilen)
C* Get restart file name
                    cg_restart_file(jf) = clvari
                ELSEIF (clvari(1:6) .eq. 'EXPOUT') THEN
                    ig_nfield = ig_nfield + 1
                    lg_state(jf) = .true.
                    ig_number_field(jf) = ig_nfield
                    ig_total_state(jf) = ip_expout
                    CALL parse(clline, clvari, 6, jpeighty, ilen)
C* Get restart file name               
                    cg_restart_file(jf) = clvari
                ELSEIF (clvari(1:6) .eq. 'IGNOUT' ) THEN
                    ig_direct_nfield = ig_direct_nfield + 1
                    lg_state(jf) = .false.
                    ig_total_state(jf) = ip_ignout
                    CALL parse(clline, clvari, 6, jpeighty, ilen)
C* Get restart file name 
                    cg_restart_file(jf) = clvari
                ELSEIF (clvari(1:9).eq. 'AUXILARY') THEN  
                    ig_nfield = ig_nfield + 1
                    lg_state(jf) = .true.
                    ig_number_field(jf) = ig_nfield
                    ig_total_state(jf) = ip_auxilary
                    CALL parse(clline, clvari, 6, jpeighty, ilen)
C* Get restart file name
                    cg_restart_file(jf) = clvari
                ELSEIF (clvari(1:5) .eq. 'INPUT') THEN
                    ig_direct_nfield = ig_direct_nfield + 1
                    lg_state(jf) = .false.
                    ig_total_state(jf) = ip_input
                    CALL parse(clline, clvari, 6, jpeighty, ilen)
C* Get input file name
                    cg_input_file(jf) = clvari
                ENDIF
            ELSE
C*          Get field status if PIPE or NONE
              CALL parse(clline, clvari, 8, jpeighty, ilen)
              IF (clvari .ne. 'EXPORTED' .and. clvari .ne. 'AUXILARY') 
     $             THEN
                 CALL prtout 
     $                ('Error in namcouple for status of field',jf,1)
                 WRITE (UNIT = nulou,FMT = *) 
     $                '==> Must be EXPORTED or AUXILARY'
                 IF (clvari(1:7) .eq. 'IGNORED') 
     $                WRITE (UNIT = nulou,FMT = *)
     $           'Direct communication is only for CLIM/MPI1 or MPI2'
                 WRITE (UNIT = nulou,FMT = *) 
     $               'Maybe you forgot the output FILE name which'
                 WRITE (UNIT = nulou,FMT = *) 
     $               'is mandatory for PIPE or NONE techniques'
                 CALL HALTE('STOP in inipar') 
              ELSE IF (clvari .eq. 'EXPORTED') THEN
                 ig_nfield = ig_nfield + 1
                 lg_state(jf) = .true.
                 ig_number_field(jf) = ig_nfield
                 ig_total_state(jf) = ip_exported
C* Get restart file name
                 CALL parse(clline, clvari, 6, jpeighty, ilen)
                 cg_restart_file(jf) = clvari
              ELSE IF (clvari .eq. 'AUXILARY') THEN
                 ig_nfield = ig_nfield + 1
                 lg_state(jf) = .true.
                 ig_number_field(jf) = ig_nfield
                 ig_total_state(jf) = ip_auxilary
C* Get restart file name
                 CALL parse(clline, clvari, 6, jpeighty, ilen)
                 cg_restart_file(jf) = clvari
              ENDIF
          ENDIF
      ENDIF
      IF (lg_state(jf)) THEN
           IF (ig_total_ntrans(jf) .eq. 0) THEN
              WRITE (UNIT = nulou,FMT = *)
     $             'If there is no analysis for the field',jf,
     $             'then the status must not be "EXPORTED"' 
              WRITE (UNIT = nulou,FMT = *)' "AUXILARY" or "EXPOUT" '
              CALL HALTE('STOP in inipar_alloc') 
           ENDIF
           READ (UNIT = nulin,FMT = 2002) clline
           CALL skip(clline, jpeighty)
           READ (UNIT = nulin,FMT = 2002) clline
           CALL skip(clline, jpeighty)
           READ (UNIT = nulin,FMT = 2002)clline_aux
           DO ja=1,ig_total_ntrans(jf)
              CALL parse(clline_aux, clvari, ja, jpeighty, ilen)
              IF (clvari.eq.'CORRECT'.or.clvari.eq.'BLASOLD'.
     $             or.clvari.eq.'BLASNEW') THEN
                 READ (UNIT = nulin,FMT = 2002) clline
                 CALL parse(clline, clvari, 2, jpeighty, ilen)
                 READ(clvari,FMT = 2003) il_aux
                 DO ib = 1, il_aux
                    READ (UNIT = nulin,FMT = 2002) clline
                    CALL skip(clline, jpeighty)
                  ENDDO
              ELSE IF (clvari.eq.'NOINTERP') THEN
                  CONTINUE
              ELSE
                  READ (UNIT = nulin,FMT = 2002) clline
                  CALL skip(clline, jpeighty)
              ENDIF
            ENDDO
        ELSE
           IF (ig_total_state(jf) .ne. ip_input) THEN
              READ (UNIT = nulin,FMT = 2002) clline
              CALL skip(clline, jpeighty)
           ENDIF
           IF (ig_total_state(jf) .ne. ip_input .and. 
     $          ig_total_ntrans(jf) .gt. 0 ) THEN
              READ (UNIT = nulin,FMT = 2002) clline
              CALL parse(clline, clvari, 1, jpeighty, ilen)
              IF (clvari(1:8) .ne. 'LOCTRANS') THEN
                 WRITE (UNIT = nulou,FMT = *)
     $              'You want a transformation which is not available !'
                 WRITE (UNIT = nulou,FMT = *)
     $              'Only local transformations are available for '
                 WRITE (UNIT = nulou,FMT = *)
     $              'fields exchanged directly or output fields '
                 CALL HALTE('STOP in inipar_alloc') 
              ENDIF
              DO ja=1,ig_total_ntrans(jf)
                 READ (UNIT = nulin,FMT = 2002) clline
                 CALL skip(clline, jpeighty)
              ENDDO
           ENDIF 
        ENDIF          
        
 240    CONTINUE
        IF (ig_nfield.eq.0) THEN
            lg_oasis_field = .false.
            WRITE (nulou,*)'==> All the fields are exchanged directly'
        ENDIF
        
C
C* Number of different restart files
C
        allocate (cl_aux(ig_total_nfield))
        cl_aux(:)=' '
        DO jf = 1,ig_total_nfield
          IF (jf.eq.1) THEN
              cl_aux(1) = cg_restart_file(1)
              il_aux = 1
          ELSEIF (jf.gt.1) THEN
              IF (ALL(cl_aux.ne.cg_restart_file(jf))) THEN
                  il_aux = il_aux + 1  
                  cl_aux(il_aux) = cg_restart_file(jf)
              ENDIF
          ENDIF
        ENDDO
        deallocate(cl_aux)
        ig_nbr_rstfile = il_aux

        IF (lg_oasis_field) THEN
C 
C*      Alloc array needed for INTERP and initialize them
C
            ALLOCATE (cintmet(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: cintmet allocation of inipar_alloc',il_err,1)
            ALLOCATE (naismfl(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: naismfl allocation of inipar_alloc',il_err,1)
            ALLOCATE (naismvoi(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: naismvoi allocation of inipar_alloc',il_err,1)
            ALLOCATE (naisgfl(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: naisgfl allocation of inipar_alloc',il_err,1)
            ALLOCATE (naisgvoi(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: naisgvoi allocation of inipar_alloc',il_err,1)
            cintmet(:)=' '
            naismfl(:) = 1
            naismvoi(:) = 1
            naisgfl(:) = 1
            naisgvoi(:) = 1
C     
C*          Alloc arrays needed for EXTRAP and initialize them
C     
            ALLOCATE (cextmet(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: cextmet allocation of inipar_alloc',il_err,1)
            ALLOCATE (nninnfl(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout
     $          ('Error: nninnfl allocation of inipar_alloc',il_err,1)
            ALLOCATE (nninnflg(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout
     $          ('Error: nninnflg allocation of inipar_alloc',il_err,1)
            ALLOCATE (neighbor(ig_nfield), stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: neighbor allocation of inipar_alloc',il_err,1)
            ALLOCATE (nextfl(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: nextfl allocation of inipar_alloc',il_err,1)
            cextmet(:)=' '
            nninnfl(:) = 1
            nninnflg(:) = 1
            neighbor(:) = 1
            nextfl(:) = 1
C     
C*          Alloc arrays needed for BLAS... analyses and initialize them 
C     
            ALLOCATE (nbofld(ig_nfield), stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: nbofld allocation of inipar_alloc',il_err,1)
            ALLOCATE (nbnfld(ig_nfield), stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: nbnfld allocation of inipar_alloc',il_err,1)
            nbofld(:) = 1
            nbnfld(:) = 1
C     
C*          Alloc arrays needed for MOZAIC and initialize them
C     
            ALLOCATE (nmapvoi(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: nmapvoi allocation of  inipar_alloc',il_err,1)
            ALLOCATE (nmapfl(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: nmapfl allocation of inipar_alloc',il_err,1)
            nmapvoi(:) = 1
            nmapfl(:) = 1
C     
C*          Alloc arrays needed for SUBGRID and initialize them
C     
            ALLOCATE (nsubfl(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: nsubfl allocation of inipar_alloc',il_err,1)
            ALLOCATE (nsubvoi(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: nsubvoi allocation of inipar_alloc',il_err,1)
            nsubfl(:) = 1
            nsubvoi(:) = 1
C     
C*          Alloc arrays needed for GLORED and REDGLO and initialize them 
C     
            ALLOCATE (ntronca(ig_nfield), stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: ntronca allocation of inipar_alloc',il_err,1)
            ntronca(:) = 0

C     
C*          Alloc array needed for analyses parameters 
C     
            ALLOCATE (cficbf(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: cficbf allocation of inipar_alloc',il_err,1)
            cficbf(:)=' '
            ALLOCATE (cficaf(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: cficaf allocation of inipar_alloc',il_err,1)
            cficaf(:)=' '
C     
C*         Alloc arrays needed for grid dimensions of direct fields and 
C*         indirect fields
C     
            ALLOCATE (nlonbf(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: nlonbf allocation of inipar_alloc',il_err,1)
            nlonbf(:)=0
            ALLOCATE (nlatbf(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: nlatbf allocation of inipar_alloc',il_err,1)
            nlatbf(:)=0
            ALLOCATE (nlonaf(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: nlonaf allocation of inipar_alloc',il_err,1)
            nlonaf(:)=0
            ALLOCATE (nlataf(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: nlataf allocation of inipar_alloc',il_err,1)
            nlataf(:)=0
C     
C*         Alloc arrays needed for grid number associated to each field
C 
            ALLOCATE (ig_grid_nbrbf(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $      ('Error: ig_grid_nbrbf allocation of inipar_alloc',il_err,1)
            ig_grid_nbrbf(:)=0
            ALLOCATE (ig_grid_nbraf(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $      ('Error: ig_grid_nbraf allocation of inipar_alloc',il_err,1)
            ig_grid_nbraf(:)=0
            
C     
C*          Alloc number of analyses array 
C     
            ALLOCATE (ig_ntrans(ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: ig_ntrans allocation of inipar_alloc',il_err,1)
            ig_ntrans(:)=0
            DO ib = 1, ig_total_nfield
              IF (lg_state(ib))
     $            ig_ntrans(ig_number_field(ib))=ig_total_ntrans(ib)
            ENDDO
C     
C*          Maximum number of analyses 
C     
            il_maxanal = imaxim(ig_ntrans,ig_nfield)
C     
C*          Alloc array of restart file names
C     
            ALLOCATE (cficinp(ig_nfield), stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: cficinp allocation of inipar_alloc',il_err,1)
            cficinp(:)=' '
            DO ib = 1, ig_total_nfield
              IF (lg_state(ib)) 
     $            cficinp(ig_number_field(ib))=cg_restart_file(ib)
            END DO
#ifdef use_netCDF
            istatus=NF_OPEN(cg_restart_file(1), NF_NOWRITE, il_id)
            IF (istatus .eq. NF_NOERR) THEN
                lncdfrst = .true.
            ELSE
#endif
                lncdfrst = .false.
#ifdef use_netCDF
            ENDIF
#endif
            WRITE(nulou, *) 'lncdfrst =', lncdfrst
C     
C*          Alloc array needed to get analysis names
C 
            ALLOCATE (canal(il_maxanal,ig_nfield),stat=il_err)
            IF (il_err.NE.0) CALL prtout 
     $          ('Error: canal allocation of inipar_alloc',il_err,1)
            canal(:,:)=' '
        ENDIF
C
C*      Get analysis parameters 
C
        REWIND nulin
 221    CONTINUE
        READ (UNIT = nulin,FMT = 2001,END = 230) clword
        IF (clword .NE. clstring) GO TO 221
C
C*      Loop on total number of fields (NoF)
C      
        DO 250 jf=1,ig_total_nfield
C
C*        Initialization
C
          nlonbf_notnc = 0
          nlatbf_notnc = 0
          nlonaf_notnc = 0
          nlataf_notnc = 0
C
C*        Skip first line read before
C
          READ (UNIT = nulin,FMT = 2002) clline
          CALL skip(clline, jpeighty)
C     
C* Second line
C
C* In the indirect case, reading of second, third, fourth line and analyses 
C* lines
C
          IF (ig_total_state(jf) .ne. ip_input) THEN
              READ (UNIT = nulin,FMT = 2002) clline
C*            First determine what information is on the line
              CALL parse(clline, clvari, 3, jpeighty, ilen)
              IF (ilen .lt. 0) THEN
C*                  
C*                IF only two words on the line, then they are the locator 
C*                prefixes and the grids file must be in NetCDF format       
                  CALL parse(clline, clvari, 1, jpeighty, ilen)
                  IF (lg_state(jf)) 
     $                cficbf(ig_number_field(jf)) = clvari
                  cga_locatorbf(jf) = clvari(1:4)
                  CALL parse(clline, clvari, 2, jpeighty, ilen)
                  IF (lg_state(jf))
     $                cficaf(ig_number_field(jf)) = clvari
                  cga_locatoraf(jf) = clvari(1:4)
                  lncdfgrd = .true.
              ELSE 
                  READ(clvari,FMT = 2010) clind, clequa, iind
                  IF (clind .EQ. 'SEQ' .or. clind .EQ. 'DEL' .or. 
     $                clind .eq. 'XTS' .or. clind .eq. 'LAG' .and.
     $                clequa .EQ. '=') THEN
                      
C*                    If 3rd word is an index, then first two words are 
C*                    locator prefixes and grids file must be NetCDF format
                      CALL parse(clline, clvari, 1, jpeighty, ilen) 
                      IF (lg_state(jf)) 
     $                    cficbf(ig_number_field(jf)) = clvari 
                      cga_locatorbf(jf) = clvari(1:4)
                      CALL parse(clline, clvari, 2, jpeighty, ilen)
                      IF (lg_state(jf)) 
     $                    cficaf(ig_number_field(jf)) = clvari
                      cga_locatoraf(jf) = clvari(1:4)
                      lncdfgrd = .true.
                  ELSE
C*              If not, the first 4 words are grid dimensions and next
C*              2 words are locator prefixes, and grids file may be or
C*              not in NetCDF format 
                      CALL parse(clline, clvari, 1, jpeighty, ilen)
C*                    Get number of longitudes for initial field
                      READ(clvari,FMT = 2004) nlonbf_notnc
                      CALL parse(clline, clvari, 2, jpeighty, ilen)
C*                    Get number of latitudes for initial field
                      READ(clvari,FMT = 2004) nlatbf_notnc
                      CALL parse(clline, clvari, 3, jpeighty, ilen)
C*                    Get number of longitudes for final field
                      READ(clvari,FMT = 2004) nlonaf_notnc
                      CALL parse(clline, clvari, 4, jpeighty, ilen)
C*                    Get number of latitudes for final field
                      READ(clvari,FMT = 2004) nlataf_notnc
                      CALL parse(clline, clvari, 5, jpeighty, ilen)
C*                    Get root name grid-related files (initial field)
                      IF (lg_state(jf)) 
     $                    cficbf(ig_number_field(jf)) = clvari
                      cga_locatorbf(jf) = clvari(1:4)
                      CALL parse(clline, clvari, 6, jpeighty, ilen)
C*                    Get root name for grid-related files (final field)
                      IF (lg_state(jf)) 
     $                    cficaf(ig_number_field(jf)) = clvari
                      cga_locatoraf(jf) = clvari(1:4)
                      nlonbf(ig_number_field(jf)) = nlonbf_notnc
                      nlatbf(ig_number_field(jf)) = nlatbf_notnc
                      nlonaf(ig_number_field(jf)) = nlonaf_notnc
                      nlataf(ig_number_field(jf)) = nlataf_notnc
C
                   ENDIF
               ENDIF
C
C*           Read the P 2 P 0 line for exported, expout or auxilary
C
             IF (lg_state(jf)) THEN
              READ (UNIT = nulin,FMT = 2002) clline
              CALL skip(clline, jpeighty)
             ENDIF
C     
C*            Read next line of strings
C             --->>> Stuff related to field transformation
C
             IF (ig_total_ntrans(jf) .gt. 0) then 
              READ (UNIT = nulin,FMT = 2002) clline
              CALL skip(clline, jpeighty)
              DO 260 ja = 1, ig_total_ntrans(jf)
                CALL parse(clline, clvari, ja, jpeighty, ilen)
C*              Get the whole set of analysis to be performed
                IF (lg_state(jf)) 
     $              canal(ja,ig_number_field(jf)) = clvari
 260          CONTINUE
C     
C*            Now read specifics for each transformation
C
           DO 270 ja = 1, ig_total_ntrans(jf)
C     
C*            Read next line unless if analysis is NOINTERP (no input)
C
             IF (lg_state(jf)) THEN
C            For EXPORTED, AUXILARY, or EXPOUT:
              IF(canal(ja,ig_number_field(jf)) .NE. 'NOINTERP') THEN 
                  READ (UNIT = nulin,FMT = 2002) clline
                  CALL skip(clline, jpeighty)
              ENDIF
              IF (canal(ja,ig_number_field(jf)) .EQ. 'MOZAIC') THEN
                  CALL parse(clline, clvari, 3, jpeighty, ilen)
C*                Get dataset number
                  READ(clvari,FMT = 2005) nmapfl(ig_number_field(jf))
                  CALL parse(clline, clvari, 4, jpeighty, ilen)
C*                Get max nbr of neighbors for the grids associated to current field
                  READ(clvari,FMT = 2003) nmapvoi(ig_number_field(jf))
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'EXTRAP')THEN
                  CALL parse(clline, clvari, 1, jpeighty, ilen)
C*                Get extrapolation method
                  cextmet(ig_number_field(jf)) = clvari
C*                If choice is NINENN, read one more data
                  IF (cextmet(ig_number_field(jf)) .EQ. 'NINENN') THEN
                      CALL parse(clline, clvari, 4, jpeighty, ilen)
C*                    Get NINENN dataset identificator
                      READ(clvari,FMT = 2005)
     $                    nninnfl(ig_number_field(jf))
                      IF (nninnfl(ig_number_field(jf)) .EQ. 0) THEN
                          WRITE(UNIT = nulou,FMT = *)'  ***WARNING***'
                          WRITE(UNIT = nulou,FMT = *) 
     $   ' **WARNING** The EXTRAP/NINENN dataset id cannot be 0' 
                          CALL HALTE('STOP in inipar')
                      ENDIF
C*                    If choice is WEIGHT, read more data
                  ELSE IF (cextmet(ig_number_field(jf)) .EQ. 'WEIGHT') 
     $                    THEN
                      CALL parse(clline, clvari, 2, jpeighty, ilen)
C  Get number of neighbors used in extrapolation
C  If extrapolation method is NINENN, next variable is the MINIMUM
C  number of neighbors required (among the 8 closest) to perform
C  the extrapolation (cannot be greater than 4 for convergence). 
C  In case it is WEIGHT, it is the MAXIMUM number
C  of neighbors required by the extrapolation operation.
C     
                      READ(clvari,FMT = 2003) 
     $                    neighbor(ig_number_field(jf))
                      CALL parse(clline, clvari, 5, jpeighty, ilen)
C*                    Get dataset number
                      READ(clvari,FMT = 2005) 
     $                    nextfl(ig_number_field(jf))
                  ENDIF
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'INTERP')THEN
                  CALL parse(clline, clvari, 1, jpeighty, ilen)
C*                Get interpolation method
                  cintmet(ig_number_field(jf)) = clvari
C*                If interpolation uses ANAIS(G-M), read in more data
                  IF (cintmet(ig_number_field(jf)) .EQ. 'SURFMESH') THEN
                      CALL parse(clline, clvari, 4, jpeighty, ilen)
C*                    Get Anaism dataset identificator
                      READ(clvari,FMT = 2005)
     $                    naismfl(ig_number_field(jf))
                      CALL parse(clline, clvari, 5, jpeighty, ilen)
C*                    Get max number of neighbors for the grids related to current field
                      READ(clvari,FMT = 2003) 
     $                    naismvoi(ig_number_field(jf))
                  ENDIF 
                  IF (cintmet(ig_number_field(jf)) .EQ. 'GAUSSIAN') THEN
                      CALL parse(clline, clvari, 4, jpeighty, ilen)
C*                    Get Anaisg dataset identificator
                      READ(clvari,FMT = 2005)
     $                    naisgfl(ig_number_field(jf))
                      CALL parse(clline, clvari, 5, jpeighty, ilen)
C*                    Get max number of neighbors for the grids related to current field
                      READ(clvari,FMT = 2003) 
     $                    naisgvoi(ig_number_field(jf))
                  ENDIF
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'SCRIPR')THEN
C* Get field type (scalar/vector)
              CALL parse(clline, clvari, 3, jpeighty, ilen)
              READ(clvari,FMT = 2009) clstrg
              IF(clstrg .EQ. 'VECTOR_I' .OR. clstrg .EQ. 'VECTOR_J')THEN
                  lg_vector  = .TRUE.
              ENDIF
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'REDGLO')THEN
                  CALL parse(clline, clvari, 1, jpeighty, ilen)
C*                Get half the nbr of lats for reduced<->global gaussian grid switch
                  READ(clvari,FMT = 2008) cldeb, 
     $                ntronca(ig_number_field(jf))
                  IF (cldeb .NE. 'NO') THEN
                      CALL prcout
     $                    ('ERROR in namcouple for analysis', 
     $                    canal(ja,ig_number_field(jf)), 1) 
                      WRITE (UNIT = nulou,FMT = *) 
     $             'Since version 2.3, the information on the reduced'
                      WRITE (UNIT = nulou,FMT = *) 
     $             'grid in namcouple has to be NOxx WHERE xx is half'
                      WRITE (UNIT = nulou,FMT = *) 
     $             'the number of latitude lines.'
                      CALL HALTE ('STOP in inipar')
                  ENDIF
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'GLORED')THEN
                  CALL parse(clline, clvari, 1, jpeighty, ilen)
C*                Get gaussian troncature for reduced <-> global gaussian grid switch
                  READ(clvari,FMT = 2008) cldeb, 
     $                ntronca(ig_number_field(jf))
                  IF (cldeb .NE. 'NO') THEN
                      CALL prcout
     $                    ('ERROR in namcouple for analysis', 
     $                    canal(ja,ig_number_field(jf)), 1) 
                      WRITE (UNIT = nulou,FMT = *) 
     $              'Since version 2.3, the information on the reduced'
                      WRITE (UNIT = nulou,FMT = *) 
     $              'grid in namcouple has to be NOxx WHERE xx is half'
                      WRITE (UNIT = nulou,FMT = *) 
     $              'the number of latitude lines.'
                      CALL HALTE ('STOP in inipar')
                  ENDIF
                  CALL parse(clline, clvari, 4, jpeighty, ilen)
C*                Get NINENN dataset identificator
                  READ(clvari,FMT = 2005) nninnflg(ig_number_field(jf))
                  IF (nninnflg(ig_number_field(jf)) .EQ. 0) THEN
                     WRITE(UNIT = nulou,FMT = *) 
     $ '**WARNING** The EXTRAP/NINENN dataset identificator in GLORED' 
                     WRITE(UNIT = nulou,FMT = *) 
     $                    'cannot be 0'
                     CALL HALTE('STOP in inipar')
                  ENDIF
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'BLASOLD') 
     $                THEN
                  CALL parse(clline, clvari, 2, jpeighty, ilen)
C*                Get number of additional fields in linear formula
                  READ(clvari,FMT = 2003) nbofld (ig_number_field(jf))
                  DO ib = 1,nbofld (ig_number_field(jf))
                    READ (UNIT = nulin,FMT = 2002) clline
                    CALL skip(clline, jpeighty)
                  ENDDO
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'BLASNEW') 
     $                THEN
                  CALL parse(clline, clvari, 2, jpeighty, ilen)
C*                Get number of additional fields in linear formula
                  READ(clvari,FMT = 2003) nbnfld (ig_number_field(jf))
                  DO ib = 1,nbnfld (ig_number_field(jf))
                    READ (UNIT = nulin,FMT = 2002) clline
                    CALL skip(clline, jpeighty)
                  ENDDO
C*                Get fields to restore subgrid variability 
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'SUBGRID') 
     $                THEN
                  CALL parse(clline, clvari, 3, jpeighty, ilen)
C*                Get dataset number
                  READ(clvari,FMT = 2005) nsubfl(ig_number_field(jf))
                  CALL parse(clline, clvari, 4, jpeighty, ilen)
C*                Get max number of neighbors for the grids related to current field
                  READ(clvari,FMT = 2003) nsubvoi(ig_number_field(jf))
              ENDIF
             ELSE
C             For IGNORED, IGNOUT and OUTPUT, only one line for LOCTRANS
              READ (UNIT = nulin,FMT = 2002) clline
              CALL skip(clline, jpeighty)
             ENDIF             
 270        CONTINUE
        ENDIF
      ENDIF
C
C*      End of loop on NoF
C     
 250  CONTINUE
C
      IF (lg_oasis_field) THEN 
C     
C*       Search maximum number of fields to be combined in the BLASxxx analyses
C     
          ig_maxcomb = imaxim(nbofld,ig_nfield)
          IF (imaxim(nbnfld,ig_nfield).gt.ig_maxcomb)
     $        ig_maxcomb = imaxim(nbnfld,ig_nfield)
C     
C*          Search maximum number of underlying neighbors for SURFMESH interpolation
C     
          ig_maxwoa = imaxim(naismvoi,ig_nfield)
          WRITE(nulou,*)
     $        'Max number of underlying neighbors for SURFMESH : ', 
     $        ig_maxwoa
          WRITE(nulou,*)' '
C     
C*          Search maximum number of neighbors for GAUSSIAN interpolation
C     
          ig_maxnoa = imaxim(naisgvoi,ig_nfield)
          WRITE(nulou,*)
     $        'Max number of neighbors for GAUSSIAN interp : ',
     $        ig_maxnoa
          WRITE(nulou,*)' '
C     
C*          Search maximum number of underlying neighbors for MOZAIC interpolation
C     
          ig_maxmoa = imaxim(nmapvoi,ig_nfield)
          WRITE(nulou,*)
     $    'Maximum number of underlying neighbors for MOZAIC interp: ',
     $        ig_maxmoa
          WRITE(nulou,*)' '
C     
C*          Search maximum number of overlaying neighbors for SUBGRID interpolation
C     
          ig_maxsoa = imaxim(nsubvoi,ig_nfield)
          WRITE(nulou,*)
     $    'Maximum number of overlaying neighbors for SUBGRID interp :',
     $        ig_maxsoa
          WRITE(nulou,*)' '
C     
C*          Search maximum number of different SURFMESH interpolations
C     
          ig_maxnfm = imaxim(naismfl,ig_nfield)
          WRITE(nulou,*)
     $        'Maximum number of different SURFMESH interpolations : ',
     $        ig_maxnfm
          WRITE(nulou,*)' '
C     
C*          Search maximum number of different GAUSSIAN interpolations
C     
          ig_maxnfg = imaxim(naisgfl,ig_nfield)
          WRITE(nulou,*)
     $        'Maximum number of different GAUSSIAN interpolations : ',
     $        ig_maxnfg
          WRITE(nulou,*)' '
C     
C*          Search maximum number of different MOZAIC interpolations
C     
          ig_maxnfp = imaxim(nmapfl,ig_nfield)
          WRITE(nulou,*)
     $        'Maximum number of different MOZAIC interpolations : ',
     $        ig_maxnfp
          WRITE(nulou,*)' '
C     
C*          Search maximum number of different SUBGRID interpolations
C     
          ig_maxnfs = imaxim(nsubfl,ig_nfield)
          WRITE(nulou,*)
     $        'Maximum number of different SUBGRID interpolations : ',
     $        ig_maxnfs
          WRITE(nulou,*)' '
C     
C*          Search maximum number of different NINENN extrapolations
C     
          ig_maxnfn = imaxim(nninnfl,ig_nfield)
          IF (imaxim(nninnflg,ig_nfield).gt.ig_maxnfn)
     $        ig_maxnfn = imaxim(nninnflg,ig_nfield)
          WRITE(nulou,*)
     $        'Maximum number of different NINENN extrapolations : ',
     $        ig_maxnfn
          WRITE(nulou,*)' '
C     
C*          Search maximum number of neighbors for extrapolation  
C     
          ig_maxext = imaxim(neighbor,ig_nfield)
          WRITE(nulou,*)
     $        'Maximum number of neighbors for extrapolation : ',
     $        ig_maxext
          WRITE(nulou,*)' '
C     
C*          Search maximum number of different extrapolation
C     
          ig_maxnbn = imaxim(nextfl,ig_nfield)
C  
      ENDIF
C*    Formats
C
 2001 FORMAT(A9)
 2002 FORMAT(A80)
 2003 FORMAT(I4)
 2004 FORMAT(I8)
 2005 FORMAT(I2)
 2008 FORMAT(A2,I4)
 2009 FORMAT(A8)
 2010 FORMAT(A3,A1,I2)
C
C*    3. End of routine
C        --------------
C
      WRITE(UNIT = nulou,FMT = *) ' '
      WRITE(UNIT = nulou,FMT = *) 
     $      '          ---------- End of routine inipar_alloc ---------'
      CALL FLUSH (nulou)
      RETURN
C
C*    Error branch output
C
 110  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *)
     $    ' No active $NBMODEL data found in input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) 
     $    ' We STOP!!! Check the file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL HALTE ('STOP in inipar_alloc')
 130  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *) 
     $    ' No active $MACHINE data found in input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) 
     $    ' We STOP!!! Check the file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL HALTE ('STOP in inipar')
 210  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *) 
     $    ' No active $FIELDS data found in input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) 
     $    ' We STOP!!! Check the file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL HALTE ('STOP in inipar_alloc')
 230  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *) 
     $    ' No active $STRING data found in input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) 
     $    ' We STOP!!! Check the file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL HALTE ('STOP in inipar_alloc')

      END
















