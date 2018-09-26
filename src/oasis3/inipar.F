      SUBROUTINE inipar
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 0 *
C               * -------------     ------- *
C               *****************************
C
C**** *inipar*  - Get run parameters
C
C     Purpose:
C     -------
C     Reads and prints out run parameters.
C
C**   Interface:
C     ---------
C       *CALL*  *inipar*
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
C       1.0   L. Terray   94/01/01  created
C       1.1   L. Terray   94/08/01  modified: change in namelist
C                                   nice flag + new case for nmode
C       1.1   L. Terray   94/10/01  modified: change printing
C       2.0b  L. Terray   95/07/24  modified: new structure
C       2.0   L. Terray   96/02/01  modified: lecture of cdqdt for
C                                   subgrid and add mozaic analysis
C                                   Lecture of a unit for filling
C       2.1   L. Terray   96/09/25  Changes to mozaic and subgrid
C                                   analysis, addition of nfend and
C                                   nintflx, check[in-out] analysis
C                                   addition of nointerp case.
C       2.2   L. Terray   97/02/12  Printing bug on analysis sub-
C                                   grid (SOLAR) corrected
C       2.2   L. Terray   97/02/20  Printing bug on analysis ANAIS
C                                   corrected
C       2.2   L. Terray   97/12/14  Add new input: MODINFO and new
C                                   extrapolation technique
C       2.3   S. Valcke   99/03/14  cjobnam with 3 or 4 characters
C       2.3   S. Valcke   99/03/25  troncature as NOxxxx in namcouple
C       2.3   S. Valcke   99/03/30  READ/WRITE flag and dataset index
C                                   for NINENN weights
C       2.3   S. Valcke   99/04/30  NLOGPRT for printing levels
C       2.3   L. Terray   99/09/15  changed periodicity variables
C                                   and input them as field parameters
C       2.4   J. Latour   99/10/28  Add new input: CHATYPE for type of
C	    	     	   message passing : MPI2 or PVM3
C       2.4   S. Valcke   2K/02/04  Additional readings for CLIM/MPI2
C       2.5   S. Valcke   2K/09/04  Remove $MACHINE, clmach, cmach
C       2.5   S. Valcke   2K/09/04  $CHATYPE in $CHANNEL
C       2.5   S. Valcke   2K/09/05  Add input line with integral flag
C                                   for  check[in-out], remove nintflx 
C       2.5   S. Valcke   2K/09/05  Remove fld 3rd input line for CLIM
C       2.5   S. Valcke   2K/09/08  Changed input lines for PVM3&MPI2 
C       2.5   J. Latour   01/11/28  Add MPI1 startup
C       2.5   A. Caubel   02/05/15  Mods for dynamic allocation
C       2.5   S. Valcke   02/06/12  PVM3 no longer supported
C       2.5   V. Gayler   01/09/20  Scrip-Remapping
C       3.0   S. Legutke  03/04/24  proposal of CF compliant cfldlab entries
C                                   added labels for echam5/mpi-om
C                                   grouping into classes
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------------------- Include files ---------------------------
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
      CHARACTER*80 clline, clvari
      CHARACTER*9 clword, clstring, clprint, clcal, clchan
      CHARACTER*9 cljob, clmod, cltime, clseq, cldate, clhead
      CHARACTER*8 cl_print_trans, cl_print_state
      CHARACTER*3 clinfo, clind
      CHARACTER*1 clequa
      CHARACTER*64 cl_cfname,cl_cfunit
      INTEGER (kind=ip_intwp_p) iind, il_aux
      INTEGER (kind=ip_intwp_p) il_file_unit, id_error
      INTEGER (kind=ip_intwp_p) il_max_entry_id, il_no_of_entries
      INTEGER (kind=ip_intwp_p) il_i, il_pos
      LOGICAL llseq, lldel, llxts, lllag, ll_exist
      INTEGER lastplace
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
     $    '           ROUTINE inipar  -  Level 0'
      WRITE (UNIT = nulou,FMT = *) 
     $    '           **************     *******'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' Initialization of run parameters'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' Reading input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
C
C* Initialize character keywords to locate appropriate input
C
      clstring = ' $STRINGS'
      cljob    = ' $JOBNAME'
      clchan   = ' $CHANNEL'
      clmod    = ' $NBMODEL'
      cltime   = ' $RUNTIME'
      clseq    = ' $SEQMODE'
      cldate   = ' $INIDATE'
      clhead   = ' $MODINFO'
      clprint  = ' $NLOGPRT'
      clcal    = ' $CALTYPE'
C
C* Initialize some variables 
      ndate = 0 ; nmseq = 1 ; ntime = 432000 ; niter = 5 
      nstep = 86400 ; nitfn=4
      cjobnam = 'DEF'
      lmodinf = .TRUE. 

C
C* CF long names for exchange fields
      INQUIRE (file='cf_name_table.txt', exist=ll_exist)

      IF (ll_exist) THEN
          WRITE (nulou,*) 'inipar: Reading CF name table!'
          il_file_unit = 99
          OPEN (file='cf_name_table.txt', unit=il_file_unit, 
     $        form='formatted', status='old')

          READ (unit=il_file_unit,fmt=*,iostat=id_error)
          READ (unit=il_file_unit,fmt=*,iostat=id_error) 
     $        il_max_entry_id, il_no_of_entries

          IF (id_error.ne.0) THEN 
              WRITE (nulou,*) 'inipar :cf_name_table.txt:' 
     $            ,' Reading of first record failed!'
              CALL halte('STOP in inipar')
          ENDIF

          IF (il_max_entry_id.gt.0) THEN 
              allocate (cfldlab(1:il_max_entry_id),STAT=id_error)
              IF (id_error.ne.0) THEN 
                  write(nulou,*) 'inipar: Allocation of cfldlab failed!'
                  CALL halte('STOP in inipar')
              ENDIF
          ELSE
              WRITE (nulou,*) 'inipar: cf_name_table.txt:', 
     $            'The number of entries is less than 0 !'
              CALL halte('STOP in inipar')                
          ENDIF

          READ (unit=il_file_unit,fmt=*,iostat=id_error)
          DO il_i=1,il_no_of_entries
            READ (unit=il_file_unit,fmt=*,iostat=id_error) 
     $          il_pos,cl_cfname,cl_cfunit

            IF (id_error.eq.0) THEN 
                IF (il_pos .le. il_max_entry_id) THEN 
                    cfldlab(il_pos)=trim(cl_cfname)
                ELSE
                    WRITE (nulou,*) 'inipar: cf_name_table.txt:',
     $               'Record ',il_i,': numlab =',il_pos,' out of range!'
                    CALL halte('STOP in inipar')  
                ENDIF
            ELSE
                WRITE (nulou,*) 'inipar: cf_name_table.txt:',
     $              'Reading record ',il_i,' failed!'
                CALL halte('STOP in inipar') 
            ENDIF
          END DO
      ELSE
          WRITE (nulou,*) 'inipar: cf_name_table.txt missing'
          CALL halte('STOP in inipar') 
      ENDIF
      CLOSE(il_file_unit)
C
C* First get experiment name 
C
      REWIND nulin
 100  CONTINUE
      READ (UNIT = nulin,FMT = 1001,END = 110) clword
      IF (clword .NE. cljob) GO TO 100
      READ (UNIT = nulin,FMT = 1002) clline
      CALL parse (clline, clvari, 1, jpeighty, ilen)
      IF (ilen .LE. 0) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Nothing on input for $JOBNAME '
          WRITE (UNIT = nulou,FMT = *) ' Default value will be used '
          WRITE (UNIT = nulou,FMT = *) ' '
        ELSE IF (ilen .GT. 0 .AND. ilen .NE. 3 .AND. ilen .NE .4 ) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Input variable length is incorrect'
          WRITE (UNIT = nulou,FMT = *) ' ilen = ', ilen  
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Check $JOBNAME variable spelling '
          WRITE (UNIT = nulou,FMT = *) ' Default value will be used '
          WRITE (UNIT = nulou,FMT = *) ' '
        ELSE
          IF (ilen .EQ. 3) THEN
              WRITE (cjobnam,FMT='(A1,A3)') ' ',clvari
          ELSE IF (ilen .EQ. 4) THEN
              WRITE (cjobnam,FMT='(A4)') clvari
          ENDIF
      ENDIF
C 
C* Print out experiment name
C
      CALL prcout
     $    ('The experiment name for this run is cjobnam =', cjobnam,1)
C
C* Get number of models involved in this simulation
C
      REWIND nulin
 120  CONTINUE
      READ (UNIT = nulin,FMT = 1001,END = 130) clword
      IF (clword .NE. clmod) GO TO 120
      READ (UNIT = nulin,FMT = 1002) clline
C
C* Get model names
C
      DO 140 jm = 1, ig_nmodel
        imodel = jm + 1
        CALL parse (clline, clvari, imodel, jpeighty, ilen)
        cmodnam(jm) = clvari
C
C* Print out model names
C
        WRITE (UNIT = nulou,FMT ='
     $      (''   Name for model '',I1,'' is '',A6,/)') 
     $      jm, cmodnam(jm)
 140  CONTINUE
C
C* Get model maximum unit number used if they appear on the line
C
      DO 142 jm = 1, ig_nmodel
        imodel = jm + 1 + ig_nmodel
        CALL parse (clline, clvari, imodel, jpeighty, ilen)
        IF (ilen .gt. 0) THEN
            READ (clvari,FMT = 1004) iga_unitmod(jm)
C
C* Print out model minimum logfile unit number
C
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT=nulou,FMT='(''The maximum Fortran unit number'',
     $          '' used in model'', I2, '' is '', I2)')
     $          jm, iga_unitmod(jm)
            WRITE (UNIT = nulou,FMT = *) ' '
C
C* Verify that maximum unit number is larger than 1024; 
C* if not, use 1024.
            IF (iga_unitmod(jm) .lt. 1024) iga_unitmod(jm)=1024
        ELSE
            WRITE (UNIT = nulou, FMT = *)
     $      ' WARNING: You did not give in the namcouple the maximum',
     $      ' Fortran unit numbers used in your models.',
     $      ' Oasis will suppose that units above 1024 are free !'
            iga_unitmod(jm)=1024
        ENDIF
 142      CONTINUE
C
C* Get hardware info for this OASIS simulation
C
      REWIND nulin
 160  CONTINUE
      READ (UNIT = nulin,FMT = 1001,END = 170) clword
      IF (clword .NE. clchan) GO TO 160
      READ (UNIT = nulin,FMT = 1002) clline
      CALL skip(clline, jpeighty)
      IF(cchan .EQ. 'MPI2' .OR. cchan .EQ. 'MPI1' 
     $    .OR. cchan .EQ. 'GSIP') THEN
C* Get one additional line for each model
          DO 186 jm = 1, ig_nmodel
            READ (UNIT = nulin,FMT = 1002) clline
C
C*    Get the launching arguments for the model
C
            CALL parseblk (clline, clvari, 3, jpeighty, ilen)
            IF (ilen .LE. 0) THEN
                WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
                WRITE (UNIT = nulou,FMT = *) 
     $     'No launching argument for model', jm
                WRITE (UNIT = nulou,FMT = *) ' '
                cmpiarg(jm)=' '
            ELSE
                cmpiarg(jm)=clvari
                WRITE (UNIT = nulou,FMT = *) ' '
                WRITE (UNIT =nulou,FMT='
     $ (''The launching argument for model '', I2, '' is'')') jm
                WRITE (UNIT = nulou,FMT = *) cmpiarg(jm)
                WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) 'ilen ',ilen
            ENDIF
            
C
 186      CONTINUE
C
         ENDIF
C
C* Get total time for this simulation
C
      REWIND nulin
 190  CONTINUE
      READ (UNIT = nulin,FMT = 1001,END = 191) clword
      IF (clword .NE. cltime) GO TO 190
      READ (UNIT = nulin,FMT = 1002) clline
      CALL parse (clline, clvari, 1, jpeighty, ilen)
      IF (ilen .LE. 0) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Nothing on input for $RUNTIME '
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Default value of 5 days will be used '
          WRITE (UNIT = nulou,FMT = *) ' '
        ELSE
          READ (clvari,FMT = 1004) ntime
      ENDIF
C
C* Print out total time
C
      CALL prtout
     $    ('The total time for this run is ntime =', ntime, 1)
C
C* Get initial date for this simulation
C
      REWIND nulin
 192  CONTINUE
      READ (UNIT = nulin,FMT = 1001,END = 193) clword
      IF (clword .NE. cldate) GO TO 192
      READ (UNIT = nulin,FMT = 1002) clline
      CALL parse (clline, clvari, 1, jpeighty, ilen)
      IF (ilen .LE. 0) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Nothing on input for $INIDATE '
          WRITE (UNIT = nulou,FMT = *) ' Default value will be used '
          WRITE (UNIT = nulou,FMT = *) ' '
        ELSE
          READ (clvari,FMT = 1004) ndate
      ENDIF
C
C* Print out initial date
C
      CALL prtout
     $    ('The initial date for this run is ndate = ', ndate, 1)
C
C* Get number of sequential models involved in this simulation
C
      REWIND nulin
 194  CONTINUE
      READ (UNIT = nulin,FMT = 1001,END = 195) clword
      IF (clword .NE. clseq) GO TO 194
      READ (UNIT = nulin,FMT = 1002) clline
      CALL parse (clline, clvari, 1, jpeighty, ilen)
      IF (ilen .LE. 0) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Nothing on input for $SEQMODE '
          WRITE (UNIT = nulou,FMT = *) ' Default value will be used '
          WRITE (UNIT = nulou,FMT = *) ' '
        ELSE IF (ilen .GT. 1) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Input variable length is incorrect'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Sequential models are too many'
          WRITE (UNIT = nulou,FMT = *) ' ilen = ', ilen  
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Check $SEQMODE variable spelling '
          WRITE (UNIT = nulou,FMT = *) ' Default value will be used '
        ELSE
          READ (clvari,FMT = 1003) nmseq
      ENDIF
C
C* Print out the number of sequential models
C
      CALL prtout
     $    ('The number of sequential fields is nmseq =', nmseq, 1)
C
C* Get the information mode for this simulation
C
      REWIND nulin
 196  CONTINUE
      READ (UNIT = nulin,FMT = 1001,END = 197) clword
      IF (clword .NE. clhead) GO TO 196
      READ (UNIT = nulin,FMT = 1002) clline
      CALL parse (clline, clvari, 1, jpeighty, ilen)
      IF (ilen .LE. 0) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Nothing on input for $MODINFO '
          WRITE (UNIT = nulou,FMT = *) ' Default value will be used '
          WRITE (UNIT = nulou,FMT = *) ' '
        ELSE IF (ilen .GT. 0 .AND. ilen .NE. 3) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Input variable length is incorrect'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Info mode uncorrectly specified'
          WRITE (UNIT = nulou,FMT = *) ' ilen = ', ilen  
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Check $MODINFO variable spelling '
          WRITE (UNIT = nulou,FMT = *) ' Default value will be used '
        ELSE
          clinfo = clvari
          IF (clinfo .EQ. 'YES') THEN 
              lmodinf = .TRUE. 
            ELSE 
              lmodinf = .FALSE. 
          ENDIF 
      ENDIF
C
C* Print out the information mode
C
      CALL prcout
     $    ('The information mode is activated ? ==>', clinfo, 1)
C
C* Get the printing level for this simulation
C
      REWIND nulin
 198  CONTINUE
      READ (UNIT = nulin,FMT = 1001,END = 199) clword
      IF (clword .NE. clprint) GO TO 198
      READ (UNIT = nulin,FMT = 1002) clline
      CALL parse (clline, clvari, 1, jpeighty, ilen)
      IF (ilen .LE. 0) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Nothing on input for $NLOGPRT '
          WRITE (UNIT = nulou,FMT = *) ' Default value 2 will be used '
          WRITE (UNIT = nulou,FMT = *) ' '
      ELSE IF (ilen .NE. 1) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Input variable length is incorrect'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Printing level uncorrectly specified'
          WRITE (UNIT = nulou,FMT = *) ' ilen = ', ilen  
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Check $NLOGPRT variable spelling '
          WRITE (UNIT = nulou,FMT = *) ' Default value will be used '
      ELSE
          READ (clvari,FMT = 1003) nlogprt
      ENDIF
C
C* Print out the printing level
C
      CALL prtout
     $    ('The printing level is nlogprt =', nlogprt, 1)
C
C* Get the calendar type for this simulation
C
      REWIND nulin
 200  CONTINUE
      READ (UNIT = nulin,FMT = 1001,END = 201) clword
      IF (clword .NE. clcal) GO TO 200
      READ (UNIT = nulin,FMT = 1002) clline
      CALL parse (clline, clvari, 1, jpeighty, ilen)
      IF (ilen .LE. 0) THEN
          WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
          WRITE (UNIT = nulou,FMT = *) 
     $        ' Nothing on input for $CALTYPE '
          WRITE (UNIT = nulou,FMT = *) ' Default value 1 will be used '
          WRITE (UNIT = nulou,FMT = *) ' '
          ncaltype = 1
      ELSE
          READ (clvari,FMT = 1003) ncaltype
      ENDIF
C
C* Print out the calendar type
C
      CALL prtout
     $   ('The calendar type is ncaltype =', ncaltype, 1)
      IF (ncaltype .EQ. 1) THEN
          CALL prcout
     $        ('Gregorian calendar', ' ', 1)
      ELSE IF (ncaltype .EQ. 0) THEN
          CALL prcout
     $        ('365 day calendar (no leap years)', ' ', 1)
      ELSE
          CALL prtout
     $        ('The number of days per month =', ncaltype, 1)
      ENDIF
C
C* Formats
C
 1001 FORMAT(A9)
 1002 FORMAT(A80)
 1003 FORMAT(I3)
 1004 FORMAT(I8)
C
C*    2. Get field information
C        ---------------------
C
C* Init. array needed for local transformation  
C
      ig_local_trans(:) = ip_instant
C
C* Init. arrays needed for ANAIS(G-M),mapping and subgrid interpolation
C
      IF (lg_oasis_field) THEN
         lcoast = .TRUE.
         DO 215 jz = 1, ig_nfield
            linit(jz) = .TRUE.
            lmapp(jz) = .TRUE.
            lsubg(jz) = .TRUE.
            lextra(jz) = .TRUE.
            varmul(jz) = 1.
            lsurf(jz) = .FALSE.
 215     CONTINUE 
C     
C* Initialize flag indicating IF EXTRAP/NINENN parameter sets have 
C* already been calculated or read (.TRUE.) or not (.FALSE.)
C     
         DO 217 jfn = 1, ig_maxnfn
            lweight(jfn) = .FALSE.
 217     CONTINUE
      ENDIF
C
C* Get the SSCS for all fields
C
      REWIND nulin
 220  CONTINUE
      READ (UNIT = nulin,FMT = 2001,END = 230) clword
      IF (clword .NE. clstring) GO TO 220
C
C  Initialize restart name index
C
      il_aux = 0
C
C* Loop on total number of fields (NoF)
C
      DO 240 jf = 1, ig_total_nfield
C
C* Read first two lines of strings for field n = 1,2...,ig_total_nfield
C      --->>> Main characteristics of fields
C
C* First line
C
         READ (UNIT = nulin,FMT = 2002) clline
         CALL parse(clline, clvari, 1, jpeighty, ilen)
C* Get output field symbolic name
         cg_input_field(jf) = clvari
        IF (lg_state(jf)) cnaminp(ig_number_field(jf)) = 
     $        cg_input_field(jf)
         IF (lg_state(jf)) cnamout(ig_number_field(jf)) = 
     $        cg_output_field(jf)
         CALL parse(clline, clvari, 3, jpeighty, ilen)
C* Get field label number
         READ (clvari,FMT = 2003) ig_numlab(jf)
         IF (lg_state(jf)) numlab(ig_number_field(jf)) = ig_numlab(jf)
         CALL parse(clline, clvari, 4, jpeighty, ilen)
C* Get field exchange frequency
         IF (clvari(1:4) .eq. 'ONCE') THEN
C
C* The case 'ONCE' means that the coupling period will be equal to the 
C* time of the simulation
C
            ig_freq(jf) = ntime
         ELSE
         READ (clvari,FMT = 2004) ig_freq(jf)
         IF (ig_freq(jf) .eq. 0) THEN
            GOTO 236
         ELSEIF (ig_freq(jf) .gt. ntime) THEN
           WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
           WRITE (UNIT = nulou,FMT = *) 
     $          'The coupling period of the field ',jf
           WRITE (UNIT = nulou,FMT = *) 
     $          'is greater than the time of the simulation '
           WRITE (UNIT = nulou,FMT = *) 
     $          'This field will not be exchanged !'
         ENDIF
         ENDIF
         IF (lg_state(jf)) nfexch(ig_number_field(jf)) = ig_freq(jf)
C* Fill up restart file number and restart file name arrays
         IF (cg_restart_file(jf).ne.' ') THEN
             IF (jf.eq.1) THEN
                 il_aux = il_aux + 1
                 ig_no_rstfile(jf) = il_aux
                 cg_name_rstfile (ig_no_rstfile(jf)) = 
     $               cg_restart_file(jf)
             ELSEIF (jf.gt.1) THEN
                 IF (ALL(cg_name_rstfile.ne.cg_restart_file(jf))) THEN
                     il_aux = il_aux + 1
                     ig_no_rstfile(jf) = il_aux
                     cg_name_rstfile (ig_no_rstfile(jf))= 
     $                   cg_restart_file(jf)
                 ELSE 
                     DO ib = 1, jf - 1 
                       IF(cg_name_rstfile(ig_no_rstfile(ib)).eq.
     $                     cg_restart_file(jf)) THEN
                           ig_no_rstfile(jf) = ig_no_rstfile(ib)
                       ENDIF
                     ENDDO
                 ENDIF
             ENDIF
         ENDIF
         CALL parse(clline, clvari, 7, jpeighty, ilen)
C* For all techniques beside PIPE and NONE technique, get eventually
C* the field STATUS
         IF (cchan .ne. 'PIPE' .and. cchan .ne. 'NONE') THEN
             IF (clvari(1:8).eq.'EXPORTED' .or. 
     $           clvari(1:8).eq.'AUXILARY') THEN
                 cstate(ig_number_field(jf)) = clvari
             ELSEIF (clvari(1:6) .eq. 'EXPOUT') THEN
                 cstate(ig_number_field(jf)) = 'EXPORTED'
             ENDIF
C*
         ELSE
            IF (lg_state(jf)) cficout(ig_number_field(jf)) = clvari
C*          Get field status
            CALL parse(clline, clvari, 8, jpeighty, ilen)
            IF (lg_state(jf)) cstate(ig_number_field(jf)) = clvari
            IF (lg_state(jf)) then
            IF (cstate(ig_number_field(jf)) .ne. 'EXPORTED'
     $          .and. cstate(ig_number_field(jf)) .ne. 'AUXILARY') THEN
                CALL prtout 
     $              ('Error in namcouple for status of field',jf,1)
                WRITE (UNIT = nulou,FMT = *) 
     $              '==> Must be EXPORTED or AUXILARY'
                WRITE (UNIT = nulou,FMT = *) 
     $              'Maybe you forgot the output FILE name which'
                WRITE (UNIT = nulou,FMT = *) 
     $              'is mandatory for PIPE or NONE techniques'
                CALL HALTE('STOP in inipar') 
            ENDIF
            ENDIF
        ENDIF
C
C* Second line
C
        IF (ig_total_state(jf) .ne. ip_input) THEN
           READ (UNIT = nulin,FMT = 2002) clline
C     *      First determine what information is on the line
           CALL parse(clline, clvari, 3, jpeighty, ilen)
           IF (ilen .lt. 0) THEN
C     *          IF only two words on the line, then they are the locator 
C     *          prefixes and the grids file must be in NetCDF format       
              ig_lag(jf)=0
              ig_total_nseqn(jf)=1
              IF (lg_state(jf)) then
                  nseqn(ig_number_field(jf)) = 1
                  nlagn(ig_number_field(jf)) = 0
              ENDIF
              llseq=.FALSE.
              lldel=.FALSE.
              llxts=.FALSE.
              lllag=.FALSE.
              WRITE (UNIT=nulou,FMT=3043) jf
              IF(nmseq .gt. 1 .and. .not. llseq) GO TO 231
           ELSE 
              READ(clvari,FMT = 2011) clind, clequa, iind
              IF (clind .EQ. 'SEQ' .or. clind .EQ. 'DEL' .or. 
     $             clind .eq. 'XTS' .or. clind .EQ. 'LAG' .and.
     $             clequa .EQ. '=') THEN
C     *              If 3rd word is an index, then first two words are 
C     *              locator prefixes and grids file must be NetCDF format
                 ilind1=3
                 ilind2=6
              ELSE
C     *              If not, the first 4 words are grid dimensions and next
C     *              2 words are locator prefixes, and grids file may be or
C     *              not in NetCDF FORMAT.
                 ilind1=7
                 ilind2=10
              ENDIF
C     *          Get possibly additional indices
              ig_lag(jf)=0
              ig_total_nseqn(jf)=1
              IF (lg_state(jf)) then
                  nseqn(ig_number_field(jf)) = 1
                  nlagn(ig_number_field(jf)) = 0
              ENDIF
              llseq=.FALSE.
              lldel=.FALSE.
              llxts=.FALSE.
              lllag=.FALSE.
C     
              DO 245 ilind=ilind1, ilind2
                 CALL parse(clline, clvari, ilind, jpeighty, ilen)
                 IF(ilen .eq. -1) THEN
                    IF (nlogprt .EQ. 2) THEN 
C                       IF(.not. lldel) WRITE (UNIT=nulou,FMT=3039) jf
C                       IF(.not. llxts) WRITE (UNIT=nulou,FMT=3040) jf
                       IF(.not. lllag) WRITE (UNIT=nulou,FMT=3043) jf
                    ENDIF
                    IF(nmseq .gt. 1 .and. .not. llseq) GO TO 231
                    GO TO 247
                 ELSE
                    READ(clvari,FMT = 2011) clind, clequa, iind
                    IF (clind .EQ. 'SEQ') THEN
                       IF (iind .gt. nmseq) THEN
                          GO TO 232
                       ELSE IF (iind .eq. 0) THEN
                          GO TO 234
                       ELSE
                          ig_total_nseqn(jf)=iind
                          IF (lg_state(jf))
     $                        nseqn(ig_number_field(jf)) = iind
                          llseq=.TRUE.
                       ENDIF
                    ELSE IF (clind .eq. 'LAG') THEN
                       ig_lag(jf)=iind
                       IF (lg_state(jf))
     $                     nlagn(ig_number_field(jf)) = iind
                       lllag=.TRUE.
                       WRITE (UNIT = nulou,FMT = 3044)jf,ig_lag(jf)
                    ENDIF              
                 ENDIF
 245          CONTINUE
          ENDIF
       ENDIF

C
 247    CONTINUE
C
C* Third line
C
        IF (lg_state(jf)) THEN
           READ (UNIT = nulin,FMT = 2002) clline
           CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get source grid periodicity type
           csper(ig_number_field(jf)) = clvari
           IF(csper(ig_number_field(jf)) .NE. 'P' .AND. 
     $          csper(ig_number_field(jf)) .NE. 'R') THEN
              CALL prtout
     $      ('ERROR in namcouple for source grid type of field', jf, 1)
              WRITE (UNIT = nulou,FMT = *) '==> must be P or R'
              CALL HALTE('STOP in inipar')
           ENDIF
C     
           CALL parse(clline, clvari, 2, jpeighty, ilen)
C     * Get nbr of overlapped longitudes for the Periodic type source grid
           READ(clvari,FMT = 2005) nosper(ig_number_field(jf))
           CALL parse(clline, clvari, 3, jpeighty, ilen)
C     * Get target grid periodicity type
           ctper(ig_number_field(jf)) = clvari
           IF(ctper(ig_number_field(jf)) .NE. 'P' .AND. 
     $          ctper(ig_number_field(jf)) .NE. 'R') THEN
              CALL prtout
     $      ('ERROR in namcouple for target grid type of field', jf, 1)
              WRITE (UNIT = nulou,FMT = *) '==> must be P or R'
              CALL HALTE('STOP in inipar')
           ENDIF
C     
           CALL parse(clline, clvari, 4, jpeighty, ilen)
C     * Get nbr of overlapped longitudes for the Periodic type target grid
           READ(clvari,FMT = 2005) notper(ig_number_field(jf))
C     
C     Define stuff related to parallel decomposition. For now, as oasis
C     is always monoproc, cparal(ig_number_field(jf))='SERIAL'. 
C     
           IF (cchan .EQ. 'MPI2' .OR. cchan .EQ. 'MPI1'
     $         .OR. cchan .EQ. 'GSIP') THEN
              cparal(ig_number_field(jf)) = 'SERIAL'
           ENDIF
       ENDIF
C
C* Get the local transformation
C
        IF (.not. lg_state(jf)) THEN
           IF (ig_total_state(jf) .ne. ip_input .and. 
     $          ig_total_ntrans(jf) .gt. 0 ) THEN
              READ (UNIT = nulin,FMT = 2002) clline
              CALL skip(clline, jpeighty)
              DO ja=1,ig_total_ntrans(jf)
                 READ (UNIT = nulin,FMT = 2002) clline 
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
                 IF (clvari(1:7) .eq. 'INSTANT') THEN 
                    ig_local_trans(jf) = ip_instant
                 ELSEIF (clvari(1:7) .eq. 'AVERAGE') THEN
                    ig_local_trans(jf) = ip_average
                 ELSEIF (clvari(1:7) .eq. 'ACCUMUL') THEN
                    ig_local_trans(jf) = ip_accumul
                 ELSEIF (clvari(1:5) .eq. 'T_MIN') THEN
                    ig_local_trans(jf) = ip_min
                 ELSEIF (clvari(1:5) .eq. 'T_MAX') THEN
                    ig_local_trans(jf) = ip_max   
                 ELSE
                    CALL prtout
     $ ('ERROR in namcouple for local transformations of field', jf, 1)
                    WRITE (UNIT = nulou,FMT = *) 
     $    '==> Must be INSTANT, AVERAGE, ACCUMUL, T_MIN or T_MAX'
                    CALL HALTE('STOP in inipar')  
                 ENDIF
              ENDDO
           ENDIF
       ELSE
         READ (UNIT = nulin,FMT = 2002) clline
              CALL skip(clline, jpeighty)
C     
C     * Now read specifics for each transformation
C 
           DO 270 ja = 1, ig_ntrans(ig_number_field(jf))
C     
C     * Read next line unless if analysis is NOINTERP (no input)
C     
              IF(canal(ja,ig_number_field(jf)) .NE. 'NOINTERP') THEN 
                 READ (UNIT = nulin,FMT = 2002) clline
                 CALL skip(clline, jpeighty)
              ENDIF
              IF (canal(ja,ig_number_field(jf)) .EQ. 'LOCTRANS') THEN
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
                 IF (clvari(1:7) .eq. 'INSTANT') THEN 
                    ig_local_trans(jf) = ip_instant
                 ELSEIF (clvari(1:7) .eq. 'AVERAGE') THEN
                    ig_local_trans(jf) = ip_average
                 ELSEIF (clvari(1:7) .eq. 'ACCUMUL') THEN
                    ig_local_trans(jf) = ip_accumul
                 ELSEIF (clvari(1:5) .eq. 'T_MIN') THEN
                    ig_local_trans(jf) = ip_min
                 ELSEIF (clvari(1:5) .eq. 'T_MAX') THEN
                    ig_local_trans(jf) = ip_max   
                 ELSE
                    CALL prtout
     $ ('ERROR in namcouple for local transformations of field', jf, 1)
                    WRITE (UNIT = nulou,FMT = *) 
     $    '==> Must be INSTANT, AVERAGE, ACCUMUL, T_MIN or T_MAX'
                    CALL HALTE('STOP in inipar')  
                 ENDIF
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'MASK') THEN
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get mask value
                 READ(clvari,FMT = 2006) amskval(ig_number_field(jf))
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'MASKP') THEN
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get the Mask value for the output field
                 READ(clvari,FMT = 2006)amskvalnew(ig_number_field(jf))
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'MOZAIC')THEN
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get file name for grid mapping
                 cgrdmap(ig_number_field(jf)) = clvari
                 CALL parse(clline, clvari, 2, jpeighty, ilen)
C     * Get related logical unit 
                 READ(clvari,FMT = 2005) nlumap(ig_number_field(jf))
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'INVERT')THEN
                 ig_invert(jf) = 1
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get lat-lon ordering for initial fields
                 cxordbf(ig_number_field(jf)) = clvari
                 IF(trim(adjustl(clvari)).eq.'NORSUD')
     $                    ig_invert(jf)=ig_invert(jf)+1
                 CALL parse(clline, clvari, 2, jpeighty, ilen)
                 cyordbf(ig_number_field(jf)) = clvari
                 IF(trim(adjustl(clvari)).eq.'ESTWST')
     $                    ig_invert(jf)=ig_invert(jf)+2
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'CHECKIN')THEN
C     * Get field integral flag
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
                 READ(clvari,FMT = 2010) clind, clequa, 
     $                ntinpflx(ig_number_field(jf))
                 IF(clind .NE. 'INT') GO TO 235
                 IF (ntinpflx(ig_number_field(jf)) .eq. 1) 
     $                lsurf(ig_number_field(jf))= .TRUE. 
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'CHECKOUT') 
     $                THEN
C     *Get field integral flag
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
                 READ(clvari,FMT = 2010) clind, clequa, 
     $                ntoutflx(ig_number_field(jf))
                 IF(clind .NE. 'INT') GO TO 235
                 IF (ntoutflx(ig_number_field(jf)) .eq. 1) 
     $                lsurf(ig_number_field(jf))= .TRUE.
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'NOINTERP') 
     $                THEN
C     * No interpolation case
                 CONTINUE
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'REVERSE') 
     $                THEN
                 ig_reverse(jf) = 1
C     * Get lat-lon ordering for final fields
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
                 cxordaf(ig_number_field(jf)) = clvari
                 IF(trim(adjustl(clvari)).eq.'NORSUD')
     $                    ig_reverse(jf)=ig_reverse(jf)+1
                 CALL parse(clline, clvari, 2, jpeighty, ilen)
                 cyordaf(ig_number_field(jf)) = clvari
                 IF(trim(adjustl(clvari)).eq.'ESTWST')
     $                    ig_reverse(jf)=ig_reverse(jf)+2
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'EXTRAP')THEN
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get extrapolation method
                 cextmet(ig_number_field(jf)) = clvari
                 CALL parse(clline, clvari, 2, jpeighty, ilen)
C     * Get number of neighbors used in extrapolation
C     If extrapolation method is NINENN, next variable is the MINIMUM
C     number of neighbors required (among the 8 closest) to perform
C     the extrapolation (cannot be greater than 4 for convergence). 
C     In case it is WEIGHT, it is the MAXIMUM number
C     of neighbors required by the extrapolation operation.
C     
                 READ(clvari,FMT = 2003) neighbor(ig_number_field(jf))
                 IF (cextmet(ig_number_field(jf)) .EQ. 'NINENN' .AND. 
     $                neighbor(ig_number_field(jf)) .GT. 4) THEN
                    neighbor(ig_number_field(jf))=4
                    WRITE(UNIT = nulou,FMT = *) '        ***WARNING***'
                    WRITE(UNIT = nulou,FMT = *) 
     $                   'For EXTRAP/NINENN extrapolation' 
                    WRITE(UNIT = nulou,FMT = *) 
     $                   'the number of neighbors has been set to 4'
                 ENDIF
C     * If choice is NINENN, read one more data
                 IF (cextmet(ig_number_field(jf)) .EQ. 'NINENN') THEN
                    CALL parse(clline, clvari, 3, jpeighty, ilen)
C     * Get NINENN weights read/write flag
                    READ(clvari,FMT = 2005) niwtn(ig_number_field(jf))
                 ENDIF
C     * If choice is WEIGHT, read more data
                 IF (cextmet(ig_number_field(jf)) .EQ. 'WEIGHT') THEN
                    CALL parse(clline, clvari, 3, jpeighty, ilen)
C     * Get file name for grid mapping
                    cgrdext(ig_number_field(jf)) = clvari
                    CALL parse(clline, clvari, 4, jpeighty, ilen)
C     * Get related logical unit 
                    READ(clvari,FMT = 2005) nluext(ig_number_field(jf))
                 ENDIF 
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'INTERP')THEN
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get interpolation method
                 cintmet(ig_number_field(jf)) = clvari
                 CALL parse(clline, clvari, 2, jpeighty, ilen)
C     * Get source grid type
                 cgrdtyp(ig_number_field(jf)) = clvari
                 CALL parse(clline, clvari, 3, jpeighty, ilen)
C     * Get field type (scalar or vector)
                 cfldtyp(ig_number_field(jf)) = clvari
C     * If interpolation uses ANAIS(G-M), read in more data
                 IF (cintmet(ig_number_field(jf)) .EQ. 'SURFMESH') THEN
                    CALL parse(clline, clvari, 6, jpeighty, ilen)
C     * Get Anaism weights read/write flag
                    READ(clvari,FMT = 2005) niwtm(ig_number_field(jf))
                 ENDIF 
                 IF (cintmet(ig_number_field(jf)) .EQ. 'GAUSSIAN') THEN
                    CALL parse(clline, clvari, 6, jpeighty, ilen)
C     * Read variance multiplicator for gaussian weights
                    READ(clvari,FMT = 2006) varmul(ig_number_field(jf))
                    CALL parse(clline, clvari, 7, jpeighty, ilen)
C     * Get Anaisg weights read/write flag
                    READ(clvari,FMT = 2005) niwtg(ig_number_field(jf))
                 ENDIF
          ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'SCRIPR') THEN
C* Get Scrip remapping method
              CALL parse(clline, clvari, 1, jpeighty, ilen)
              READ(clvari,FMT = 2009) cmap_method(ig_number_field(jf))
C* Get source grid type
              CALL parse(clline, clvari, 2, jpeighty, ilen)
              READ(clvari,FMT = 2009) cgrdtyp(ig_number_field(jf))
              IF (cmap_method(ig_number_field(jf)) .eq. 'BICUBIC' 
     $            .and. cgrdtyp(ig_number_field(jf)) .ne. 'LR'
     $            .and. cgrdtyp(ig_number_field(jf)) .ne. 'D') THEN
                  WRITE (UNIT = nulou,FMT = *) '    '
                  CALL prtout
     $                ('ERROR in namcouple for type of field', jf, 1)
                  WRITE (UNIT = nulou,FMT = *) 
     $    'BICUBIC interpolation cannot be used if grid is not LR or D'
                  CALL HALTE('STOP in inipar') 
              ENDIF
              IF (cmap_method(ig_number_field(jf)) .eq. 'BILINEAR' 
     $            .and. cgrdtyp(ig_number_field(jf)) .ne. 'LR'
     $            .and. cgrdtyp(ig_number_field(jf)) .ne. 'D') THEN
                  WRITE (UNIT = nulou,FMT = *) '    '
                  CALL prtout
     $                ('ERROR in namcouple for type of field', jf, 1)
                  WRITE (UNIT = nulou,FMT = *) 
     $    'BILINEAR interpolation cannot be used if grid is not LR or D'
                  CALL HALTE('STOP in inipar') 
              ENDIF
C* Get field type (scalar/vector)
              CALL parse(clline, clvari, 3, jpeighty, ilen)
              READ(clvari,FMT = 2009) cfldtype(ig_number_field(jf))
              IF(cfldtype(ig_number_field(jf)) .NE. 'SCALAR' .AND. 
     $             cfldtype(ig_number_field(jf)) .NE. 'VECTOR' .AND.
     $             cfldtype(ig_number_field(jf)) .NE. 'VECTOR_I' .AND.
     $             cfldtype(ig_number_field(jf)) .NE. 'VECTOR_J') THEN
                  WRITE (UNIT = nulou,FMT = *) '    '
                  CALL prtout
     $                ('ERROR in namcouple for type of field', jf, 1)
                  WRITE (UNIT = nulou,FMT = *) 
     $                '==> must be SCALAR, VECTOR_I or VECTOR_J'
                  CALL HALTE('STOP in inipar')
              ENDIF
C* Get restriction type for SCRIP search
              CALL parse(clline, clvari, 4, jpeighty, ilen)
              READ(clvari,FMT = 2009) crsttype(ig_number_field(jf))
              IF (cgrdtyp(ig_number_field(jf)) .EQ. 'D') THEN
                  IF (cmap_method(ig_number_field(jf)) .EQ. 'BILINEAR'
     $                    .or.
     $                 cmap_method(ig_number_field(jf)) .EQ. 'BICUBIC')
     $                    THEN
                      IF (crsttype(ig_number_field(jf)) .NE. 'LATITUDE')
     $                    THEN
                          WRITE (UNIT = nulou,FMT = *) '    '
                          CALL prtout
     $             ('ERROR in namcouple for restriction of field',jf,1)
                          WRITE (UNIT = nulou,FMT = *) 
     $             '==> LATITUDE must be chosen for reduced grids (D)'
                          CALL HALTE('STOP in inipar')
                      ELSE  
                          crsttype(ig_number_field(jf)) = 'REDUCED'
                      ENDIF
                  ENDIF
              ENDIF
C
              IF(crsttype(ig_number_field(jf)) .NE. 'LATITUDE' .AND. 
     $            crsttype(ig_number_field(jf)) .NE. 'LATLON' .AND.
     $            crsttype(ig_number_field(jf)) .NE. 'REDUCED') THEN
                  WRITE (UNIT = nulou,FMT = *) '    '
                  CALL prtout
     $            ('ERROR in namcouple for restriction of field',jf,1)
                  WRITE (UNIT = nulou,FMT = *) 
     $                '==> must be LATITUDE or LATLON'
                  CALL HALTE('STOP in inipar')
              ENDIF
C*
C* Get number of search bins for SCRIP search
              CALL parse(clline, clvari, 5, jpeighty, ilen)
              READ(clvari,FMT = 2003) nbins(ig_number_field(jf))
C* Get normalize option for CONSERV
              IF (cmap_method(ig_number_field(jf)) .EQ. 'CONSERV') THEN
                  CALL parse(clline, clvari, 6, jpeighty, ilen)
                  READ(clvari,FMT = 2009)cnorm_opt(ig_number_field(jf))
                  IF (cnorm_opt(ig_number_field(jf)) .NE. 'FRACAREA'
     $                .AND. 
     $		      cnorm_opt(ig_number_field(jf)) .NE. 'DESTAREA' 
     $                .AND. 
     $                cnorm_opt(ig_number_field(jf)) .NE. 'FRACNNEI') 
     $                THEN
                      WRITE (UNIT = nulou,FMT = *) '    '
                      CALL prtout
     $        ('ERROR in namcouple for normalize option of field',jf,1)
                      WRITE (UNIT = nulou, FMT = *) 
     $                '==> must be FRACAREA, DESTAREA, or FRACNNEI'
                      CALL HALTE('STOP in inipar')
                  ENDIF
C* Get order of remapping for CONSERV
                  CALL parse(clline, clvari, 7, jpeighty, ilen)
                  IF (ilen .LE. 0) THEN
                      WRITE (UNIT = nulou,FMT = *) '    '
                      CALL prtout
     $        ('ERROR in namcouple for CONSERV for field',jf,1)
                      WRITE (UNIT = nulou,FMT = *) 
     $        '==> SECOND or FIRST must be indicated at end of line'
                      CALL HALTE('STOP in inipar')
                  ENDIF
                  READ(clvari,FMT = 2009) corder(ig_number_field(jf))
              ELSE
                  cnorm_opt(ig_number_field(jf))='NONORM'
              ENDIF
C* Get number of neighbours for DISTWGT and GAUSWGT
              IF (cmap_method(ig_number_field(jf)) .EQ. 'DISTWGT' .or.
     $            cmap_method(ig_number_field(jf)) .EQ. 'GAUSWGT') THEN
                  CALL parse(clline, clvari, 6, jpeighty, ilen)
                  IF (ilen .LE. 0) THEN
                      WRITE (UNIT = nulou,FMT = *) '    '
                      CALL prtout
     $        ('ERROR in namcouple for field',jf,1)
                      WRITE (UNIT = nulou,FMT = *) 
     $        '==> Number of neighbours must be indicated on the line'
                      CALL HALTE('STOP in inipar')
                  ELSE
                    READ(clvari,FMT=2003)nscripvoi(ig_number_field(jf))
                  ENDIF 
              ENDIF
C* Get gaussian variance for GAUSWGT
              IF (cmap_method(ig_number_field(jf)) .EQ. 'GAUSWGT') THEN
                  CALL parse(clline, clvari, 7, jpeighty, ilen)
                  IF (ilen .LE. 0) THEN
                      WRITE (UNIT = nulou,FMT = *) '    '
                      CALL prtout
     $        ('ERROR in namcouple for GAUSWGT for field',jf,1)
                      WRITE (UNIT = nulou,FMT = *) 
     $        '==> Variance must be indicated at end of line'
                      CALL HALTE('STOP in inipar')
                  ELSE
                      READ(clvari,FMT=2006) varmul(ig_number_field(jf))
                  ENDIF
              ENDIF
C*Get associated file name and information about rotation to cartesien
              IF (cfldtype(ig_number_field(jf))=='VECTOR_I' .or.
     $             cfldtype(ig_number_field(jf))=='VECTOR_J') THEN
                 IF(cmap_method(ig_number_field(jf)) .EQ. 'DISTWGT')
     $                lastplace=7
                 IF(cmap_method(ig_number_field(jf)) .EQ. 'GAUSWGT')
     $                lastplace=8
                 IF(cmap_method(ig_number_field(jf)) .EQ. 'BILINEAR')
     $                lastplace=6
                 IF(cmap_method(ig_number_field(jf)) .EQ. 'BICUBIC')
     $                lastplace=6
                 IF(cmap_method(ig_number_field(jf)) .EQ. 'CONSERV')
     $                lastplace=8
                 CALL parse(clline, clvari, lastplace, jpeighty, ilen)
                 IF (ilen .le. 0) THEN
                      WRITE (UNIT = nulou,FMT = *) ' '
                      WRITE (UNIT = nulou,FMT = *) 
     $                    '==> A field associated must be indicated'
                      CALL HALTE('STOP in inipar')
                 ENDIF
                 cg_assoc_input_field(ig_number_field(jf))=clvari
C*Rotation?
                 CALL parse(clline, clvari, lastplace+1, jpeighty, 
     $                ilen) 
                 IF (ilen .le. 0) THEN
                    lrotate(ig_number_field(jf)) = .false.
                 ELSEIF(clvari .le. 'PROJCART') THEN
                    lrotate(ig_number_field(jf)) = .true.
                    WRITE (UNIT = nulou,FMT = *)
     $                   'rotation to cartesian for field : ', jf
                 ELSE
                    WRITE (UNIT = nulou,FMT = *) ' '
                    CALL prtout
     $                   ('ERROR in namcouple for vector in SCRIPR
     $                   for field',jf,1)     
                    WRITE (UNIT = nulou,FMT = *)
     $                   'must be PROJCART or nothing' 
                    CALL HALTE('STOP in inipar') 
                 ENDIF
              END IF
C
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'FILLING') 
     $                THEN
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get data file name (used to complete the initial field array)
                 cfilfic(ig_number_field(jf)) = clvari
                 CALL parse(clline, clvari, 2, jpeighty, ilen)
C     * Get logical unit connected to previous file
                 READ(clvari,FMT = 2005) nlufil(ig_number_field(jf))
                 CALL parse(clline, clvari, 3, jpeighty, ilen)
C     * Get filling method
                 cfilmet(ig_number_field(jf)) = clvari
C     * If current field is SST
                 IF(cfilmet(ig_number_field(jf))(4:6) .EQ. 'SST') THEN
                    CALL parse(clline, clvari, 4, jpeighty, ilen)
C     * Get flag for coast mismatch correction
                    READ(clvari,FMT = 2005) nfcoast
                    IF (cfilmet(ig_number_field(jf))(1:3) .EQ. 'SMO')
     $                  THEN
                        CALL parse(clline, clvari, 5, jpeighty, ilen)
C     * Get field name for flux corrective term 
                        cfldcor = clvari
                        CALL parse(clline, clvari, 6, jpeighty, ilen)
C     * Get logical unit used to write flux corrective term
                        READ(clvari,FMT = 2005) nlucor
                    ENDIF
                 ENDIF 
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'CONSERV') 
     $                THEN            
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get conservation method
                 cconmet(ig_number_field(jf)) = clvari
                 lsurf(ig_number_field(jf)) = .TRUE.
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'REDGLO') THEN
C     * Get extrapolation flag to go from reduced to global gaussian grid
                 CALL parse(clline, clvari, 2, jpeighty, ilen)
                 cmskrd(ig_number_field(jf)) = clvari
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'GLORED') THEN
                 CALL parse(clline, clvari, 2, jpeighty, ilen)
C     * Get number of neighbors used in EXTRAP/NINENN extrapolation always
C     performed within GLORED (cannot be greater than 4 for convergence).
                 READ(clvari,FMT = 2003) neighborg(ig_number_field(jf))
                 CALL parse(clline, clvari, 3, jpeighty, ilen)
                 IF (neighborg(ig_number_field(jf)) .GT. 4) THEN
                    neighborg(ig_number_field(jf))=4
                    WRITE(UNIT = nulou,FMT = *) '        ***WARNING***'
                    WRITE(UNIT = nulou,FMT = *) 
     $                   'For EXTRAP/NINENN extrapolation in GLORED' 
                    WRITE(UNIT = nulou,FMT = *) 
     $                   'the number of neighbors has been set to 4'
                 ENDIF
C     * Get EXTRAP/NINENN weights read/write flag
                 READ(clvari,FMT = 2005) niwtng(ig_number_field(jf))
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'CORRECT') 
     $                THEN
C     * Get flux correction parameters
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get main field multiplicative coefficient
                 READ(clvari,FMT = 2006) afldcoef(ig_number_field(jf))
                 CALL parse(clline, clvari, 2, jpeighty, ilen)
C     * Get number of auxilary fields in correction formula
                 READ(clvari,FMT = 2003) ncofld (ig_number_field(jf))
C     * Read auxilary field parameters
                 icofld = ncofld(ig_number_field(jf))
                 DO 280 jc = 1, icofld
                    READ (UNIT = nulin,FMT = 2002) clline   
                    CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get symbolic names for additional fields
                    ccofld(jc,ig_number_field(jf)) = clvari
                    CALL parse(clline, clvari, 2, jpeighty, ilen)
C     * Get multiplicative coefficients for  additional fields
                    READ(clvari,FMT = 2006) 
     $                   acocoef (jc,ig_number_field(jf))
                    CALL parse(clline, clvari, 3, jpeighty, ilen)
C     * Get file names for external data files 
                    ccofic(jc,ig_number_field(jf)) = clvari
C     * Get related logical units 
                    CALL parse(clline, clvari, 4, jpeighty, ilen)
                    READ(clvari,FMT = 2005) 
     $                   nludat(jc,ig_number_field(jf))
 280             CONTINUE
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'BLASOLD')THEN
C     * Get linear combination parameters for initial fields
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get main field multiplicative coefficient
                 READ(clvari,FMT = 2006) afldcobo(ig_number_field(jf))
                 DO 290 jc = 1, nbofld(ig_number_field(jf))
                    READ (UNIT = nulin,FMT = 2002) clline   
                    CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get symbolic names for additional fields
                    cbofld(jc,ig_number_field(jf)) = clvari
                    CALL parse(clline, clvari, 2, jpeighty, ilen)
C     * Get multiplicative coefficients for  additional fields
                    READ(clvari,FMT = 2006) 
     $                   abocoef (jc,ig_number_field(jf))
 290             CONTINUE
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'BLASNEW')THEN
C     * Get linear combination parameters for final fields
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get main field multiplicative coefficient
                 READ(clvari,FMT = 2006) afldcobn(ig_number_field(jf))
                 DO 291 jc = 1, nbnfld(ig_number_field(jf))
                    READ (UNIT = nulin,FMT = 2002) clline   
                    CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get symbolic names for additional fields
                    cbnfld(jc,ig_number_field(jf)) = clvari
                    CALL parse(clline, clvari, 2, jpeighty, ilen)
C     * Get multiplicative coefficients for  additional fields
                    READ(clvari,FMT = 2006) 
     $                   abncoef (jc,ig_number_field(jf))
 291             CONTINUE
C     * Get fields to restore subgrid variability 
              ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'SUBGRID')THEN
                 CALL parse(clline, clvari, 1, jpeighty, ilen)
C     * Get file name for subgrid interpolation
                 cgrdsub(ig_number_field(jf)) = clvari
                 CALL parse(clline, clvari, 2, jpeighty, ilen)
C     * Get related logical unit 
                 READ(clvari,FMT = 2005) nlusub(ig_number_field(jf))
                 CALL parse(clline, clvari, 5, jpeighty, ilen)
C     * Get type of subgrid interpolation (solar or non solar flux)
                 ctypsub(ig_number_field(jf)) = clvari
                 CALL parse(clline, clvari, 6, jpeighty, ilen)
C     * Get additional field name on coarse grid
                 cfldcoa(ig_number_field(jf)) = clvari
                 CALL parse(clline, clvari, 7, jpeighty, ilen)
C     * Get additional field name on fine grid
                 cfldfin(ig_number_field(jf)) = clvari
                 IF (ctypsub(ig_number_field(jf)) .EQ. 'NONSOLAR') THEN
                    CALL parse(clline, clvari, 8, jpeighty, ilen)
C     * Get coupling ratio on coarse grid
                    cdqdt(ig_number_field(jf)) = clvari
                 ENDIF 
              ELSE 
                 WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
                 WRITE (UNIT = nulou,FMT = *)
     $                ' Type of analysis not implemented yet '
                 WRITE (UNIT = nulou,FMT = *) 
     $                ' The analysis required in OASIS is :'
                 WRITE (UNIT = nulou,FMT = *) ' canal = ', 
     $                canal(ja,ig_number_field(jf))
                 WRITE (UNIT = nulou,FMT = *) 
     $                ' with ja = ', ja, ' jf = ', jf
                 WRITE (UNIT = nulou,FMT = *) ' '
                 CALL HALTE ('STOP in inipar')
             ENDIF
 270       CONTINUE
       ENDIF
C
C* End of loop on NoF
C 
 240  CONTINUE
C
C*Get the associated number file for case vector
C
      DO jf = 1, ig_total_nfield
         IF (cfldtype(ig_number_field(jf))=='VECTOR_I' .or.
     $        cfldtype(ig_number_field(jf))=='VECTOR_J') THEN
            DO jff = 1, ig_total_nfield
               IF(cnaminp(ig_number_field(jff)) .eq.
     $              cg_assoc_input_field(ig_number_field(jf)))THEN
                  ig_assoc_input_field(ig_number_field(jf))=
     $                 ig_number_field(jff)
                  exit
               ENDIF
            ENDDO
C
C*Verify if interpolations are the same for the 2 components of the vector field
C
            IF(cmap_method(ig_number_field(jf)) .ne. 
     $           cmap_method(ig_assoc_input_field(
     $           ig_number_field(jf)))) THEN
               WRITE (UNIT = nulou,FMT = *) 
     $              'Interpolations must be the same for the 2'
               WRITE (UNIT = nulou,FMT = *)
     $              'components in vector case'
               CALL HALTE('STOP in inipar')
            END IF
         ENDIF
      ENDDO
C
C* Minimum coupling period
C
      ig_total_frqmin = iminim(ig_freq, ig_total_nfield)
C
C* Formats
C
 2001 FORMAT(A9)
 2002 FORMAT(A80)
 2003 FORMAT(I4)
 2004 FORMAT(I8)
 2005 FORMAT(I2)
 2006 FORMAT(E15.6)
 2008 FORMAT(A2,I4)
 2009 FORMAT(A8)
 2010 FORMAT(A3,A1,I2)
 2011 FORMAT(A3,A1,I8)
C
C*    3. Printing
C        --------
C* Warning: no indentation for the next if (nightmare ...)
      IF (nlogprt .GE. 1) THEN 
      DO 310 jf = 1, ig_total_nfield
         IF (ig_total_state(jf) .eq. ip_exported ) THEN
            cl_print_state = 'EXPORTED'
         ELSEIF (ig_total_state(jf) .eq. ip_ignored ) THEN
            cl_print_state = 'IGNORED'
         ELSEIF (ig_total_state(jf) .eq. ip_ignout ) THEN
            cl_print_state = 'IGNOUT'
         ELSEIF (ig_total_state(jf) .eq. ip_expout ) THEN 
            cl_print_state = 'EXPOUT'
         ELSEIF (ig_total_state(jf) .eq. ip_input ) THEN 
            cl_print_state = 'INPUT'
         ELSEIF (ig_total_state(jf) .eq. ip_output ) THEN 
            cl_print_state = 'OUTPUT'
         ELSEIF (ig_total_state(jf) .eq. ip_auxilary ) THEN 
            cl_print_state = 'AUXILARY'
         ENDIF
         IF (ig_local_trans(jf) .eq. ip_instant) THEN
            cl_print_trans = 'INSTANT'
         ELSEIF (ig_local_trans(jf) .eq. ip_average) THEN
             cl_print_trans = 'AVERAGE'
         ELSEIF (ig_local_trans(jf) .eq. ip_accumul) THEN
            cl_print_trans = 'ACCUMUL'
         ELSEIF (ig_local_trans(jf) .eq. ip_min) THEN
            cl_print_trans = 'T_MIN'
         ELSEIF (ig_local_trans(jf) .eq. ip_max) THEN
            cl_print_trans = 'T_MAX'   
         ENDIF
C* Local indexes
      IF (.not. lg_state(jf)) THEN
         ilab = ig_numlab(jf)
         WRITE (UNIT = nulou,FMT = 3001) jf
         WRITE (UNIT = nulou,FMT = 3002)
         WRITE (UNIT = nulou,FMT = 3003)
         WRITE (UNIT = nulou,FMT = 3004)
         IF (ig_total_state(jf) .eq. ip_input .or. 
     $        ig_total_state(jf) .eq. ip_output) THEN
            WRITE (UNIT = nulou,FMT = 3121)
     $           cg_input_field(jf), cg_output_field(jf), cfldlab(ilab), 
     $           ig_freq(jf), cl_print_trans,
     $           cl_print_state, ig_total_ntrans(jf)
         ELSE  
            WRITE (UNIT = nulou,FMT = 3116)
     $           cg_input_field(jf), cg_output_field(jf), cfldlab(ilab), 
     $           ig_freq(jf), cl_print_trans, ig_total_nseqn(jf), 
     $           ig_lag(jf), cl_print_state, ig_total_ntrans(jf)
         ENDIF
      ELSE
         ilab = numlab(ig_number_field(jf))
         ifcb = ilenstr(cficbf(ig_number_field(jf)),jpeight)
         ifca = ilenstr(cficaf(ig_number_field(jf)),jpeight)
         WRITE (UNIT = nulou,FMT = 3001) jf
         WRITE (UNIT = nulou,FMT = 3002)
         WRITE (UNIT = nulou,FMT = 3003)
         WRITE (UNIT = nulou,FMT = 3004) 
         IF (cchan.EQ.'MPI2' .OR. cchan.EQ.'MPI1' 
     $       .OR. cchan.EQ.'GSIP') THEN
            WRITE (UNIT = nulou,FMT = 3005)
     $           cnaminp(ig_number_field(jf)), 
     $           cnamout(ig_number_field(jf)), cfldlab(ilab), 
     $           nfexch(ig_number_field(jf)),
     $           nseqn(ig_number_field(jf)),
     $           ig_lag(jf),
     $           cl_print_state,
     $           ig_ntrans(ig_number_field(jf)), 
     $           cparal(ig_number_field(jf))
         ELSE
               WRITE (UNIT = nulou,FMT = 3115)
     $           cnaminp(ig_number_field(jf)), 
     $           cnamout(ig_number_field(jf)), cfldlab(ilab), 
     $           nfexch(ig_number_field(jf)),
     $           nseqn(ig_number_field(jf)), 
     $           cstate(ig_number_field(jf)), 
     $           ig_ntrans(ig_number_field(jf))
         ENDIF
      ENDIF
C* Warning: no indentation for the next if (nightmare ...)
        IF (nlogprt .EQ. 2) THEN
C* Warning: no indentation for the next if (nightmare ...)            
        IF (.not. lg_state(jf)) THEN
           IF (ig_total_state(jf) .eq. ip_ignored .or. 
     $         ig_total_state(jf) .eq. ip_ignout ) THEN
              WRITE (UNIT = nulou,FMT = 3117) cg_restart_file(jf)
           ELSEIF (ig_total_state(jf) .eq. ip_input) THEN
              WRITE (UNIT = nulou,FMT = 3118) cg_input_file(jf)
           ENDIF
        ELSE
           IF (ig_total_state(jf) .eq. ip_exported .or. 
     $          ig_total_state(jf) .eq. ip_expout .or. 
     $          ig_total_state(jf) .eq. ip_auxilary )
     $          WRITE (UNIT = nulou,FMT = 3117) cg_restart_file(jf)
C* Warning: no indentation for the next if (nightmare ...)           
        WRITE (UNIT = nulou,FMT = 3007)
     $      csper(ig_number_field(jf)), nosper(ig_number_field(jf)), 
     $      ctper(ig_number_field(jf)), notper(ig_number_field(jf))
        WRITE (UNIT = nulou,FMT = 3008)
     $      cficbf(ig_number_field(jf))(1:ifcb)//cglonsuf, 
     $       cficbf(ig_number_field(jf))(1:ifcb)//cglatsuf,
     $      cficbf(ig_number_field(jf))(1:ifcb)//cmsksuf, 
     $       cficbf(ig_number_field(jf))(1:ifcb)//csursuf,
     $      cficaf(ig_number_field(jf))(1:ifca)//cglonsuf, 
     $       cficaf(ig_number_field(jf))(1:ifca)//cglatsuf,
     $      cficaf(ig_number_field(jf))(1:ifca)//cmsksuf, 
     $       cficaf(ig_number_field(jf))(1:ifca)//csursuf
        WRITE (UNIT = nulou,FMT = 3009) 
        WRITE (UNIT = nulou,FMT = 3010)
        DO 320 ja = 1, ig_ntrans(ig_number_field(jf))
          WRITE (UNIT = nulou,FMT = 3011) ja, 
     $          canal(ja,ig_number_field(jf))
          IF (canal(ja,ig_number_field(jf)) .EQ. 'MASK') THEN
              WRITE(UNIT = nulou,FMT = 3012) 
     $            amskval(ig_number_field(jf))
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'MASKP') THEN
              WRITE(UNIT = nulou,FMT = 3042) 
     $              amskvalnew(ig_number_field(jf))
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'MOZAIC') THEN
              WRITE(UNIT = nulou,FMT = 3013) 
     $              cgrdmap(ig_number_field(jf)), 
     $              nlumap(ig_number_field(jf)),
     $              nmapfl(ig_number_field(jf)), 
     $              nmapvoi(ig_number_field(jf))
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'INVERT') THEN
              WRITE(UNIT = nulou,FMT = 3014) 
     $              cxordbf(ig_number_field(jf))
              WRITE(UNIT = nulou,FMT = 3015) 
     $             cyordbf(ig_number_field(jf))
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'REVERSE') THEN
              WRITE(UNIT = nulou,FMT = 3016) 
     $              cxordaf(ig_number_field(jf))
              WRITE(UNIT = nulou,FMT = 3017) 
     $             cyordaf(ig_number_field(jf))
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'EXTRAP') THEN
              WRITE(UNIT = nulou,FMT = 3018) 
     $              cextmet(ig_number_field(jf)), 
     $              neighbor(ig_number_field(jf))
              IF (cextmet(ig_number_field(jf)) .EQ. 'WEIGHT') THEN 
                  WRITE(UNIT = nulou,FMT = 3019) 
     $                cgrdext(ig_number_field(jf)), 
     $                nluext(ig_number_field(jf)), 
     $                nextfl(ig_number_field(jf))
              ELSE IF (cextmet(ig_number_field(jf)) .EQ. 'NINENN') THEN 
                  WRITE(UNIT = nulou,FMT = 3038) 
     $                niwtn(ig_number_field(jf)), 
     $                nninnfl(ig_number_field(jf))
              ENDIF
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'INTERP') THEN
              WRITE(UNIT = nulou,FMT = 3020) 
     $              cintmet(ig_number_field(jf)), 
     $              cgrdtyp(ig_number_field(jf)),
     $              cfldtyp(ig_number_field(jf))
              IF (cintmet(ig_number_field(jf)) .EQ. 'SURFMESH') THEN 
                  WRITE(UNIT = nulou,FMT = 3021) 
     $                naismfl(ig_number_field(jf)), 
     $                naismvoi(ig_number_field(jf)), 
     $                niwtm(ig_number_field(jf))
              ENDIF 
              IF (cintmet(ig_number_field(jf)) .EQ. 'GAUSSIAN') THEN 
                  WRITE(UNIT = nulou,FMT = 3021) 
     $                naisgfl(ig_number_field(jf)), 
     $                naisgvoi(ig_number_field(jf)), 
     $                niwtg(ig_number_field(jf))
                  WRITE(UNIT = nulou,FMT = 3022) 
     $                 varmul(ig_number_field(jf))
              ENDIF
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'SCRIPR') THEN
              WRITE(UNIT = nulou,FMT = 3045) 
     $              cmap_method(ig_number_field(jf)), 
     $              cfldtype(ig_number_field(jf)), 
     $              cnorm_opt(ig_number_field(jf)),
     $              crsttype(ig_number_field(jf)), 
     $              nbins(ig_number_field(jf))
              IF (cmap_method(ig_number_field(jf)) .EQ. 'CONSERV') THEN 
                  WRITE(UNIT = nulou,FMT = 3046) 
     $                corder(ig_number_field(jf))
              ENDIF  
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'FILLING') THEN
              WRITE(UNIT = nulou,FMT = 3023) 
     $              cfilfic(ig_number_field(jf)), 
     $              nlufil(ig_number_field(jf)),
     $              cfilmet(ig_number_field(jf))
              IF(cfilmet(ig_number_field(jf))(1:6) .EQ. 'SMOSST')
     $            WRITE(UNIT = nulou,FMT = 3024) 
     $            nfcoast, cfldcor, nlucor
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'CONSERV') THEN            
              WRITE(UNIT = nulou,FMT = 3025) 
     $              cconmet(ig_number_field(jf))
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'REDGLO') THEN
              WRITE(UNIT = nulou,FMT = 3026) 
     $              ntronca(ig_number_field(jf)), 
     $              cmskrd(ig_number_field(jf))
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'CORRECT') THEN
              WRITE(UNIT = nulou,FMT = 3027) 
     $              cnamout(ig_number_field(jf)), 
     $              afldcoef(ig_number_field(jf))
              WRITE(UNIT = nulou,FMT=3028) ncofld(ig_number_field(jf))
              icofld = ncofld(ig_number_field(jf))
              DO 330 jc = 1, icofld
                WRITE(UNIT = nulou,FMT = 3029) 
     $              ccofic(jc,ig_number_field(jf)),
     $                nludat(jc,ig_number_field(jf))
                WRITE (UNIT = nulou,FMT = 3030) 
     $              ccofld(jc,ig_number_field(jf)), 
     $               acocoef(jc,ig_number_field(jf))
 330          CONTINUE
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'BLASOLD') THEN
              WRITE(UNIT = nulou,FMT = 3027) 
     $              cnaminp(ig_number_field(jf)), 
     $              afldcobo(ig_number_field(jf))
              WRITE(UNIT = nulou,FMT=3028) nbofld(ig_number_field(jf))
              DO 340 jc = 1, nbofld(ig_number_field(jf))
                WRITE (UNIT = nulou,FMT = 3030) 
     $              cbofld(jc,ig_number_field(jf)), 
     $                abocoef (jc,ig_number_field(jf))
 340          CONTINUE
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'BLASNEW') THEN
              WRITE(UNIT = nulou,FMT = 3027) 
     $              cnamout(ig_number_field(jf)), 
     $              afldcobn(ig_number_field(jf))
              WRITE(UNIT = nulou,FMT=3028) nbnfld(ig_number_field(jf))
              DO 350 jc = 1, nbnfld(ig_number_field(jf))
                WRITE (UNIT = nulou,FMT = 3030) 
     $              cbnfld(jc,ig_number_field(jf)), 
     $                abncoef (jc,ig_number_field(jf))
 350          CONTINUE
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'SUBGRID') THEN
              WRITE(UNIT = nulou,FMT = 3031) 
     $              cgrdsub(ig_number_field(jf)), 
     $              nlusub(ig_number_field(jf)),
     $              nsubfl(ig_number_field(jf)), 
     $              nsubvoi(ig_number_field(jf)), 
     $              ctypsub(ig_number_field(jf))
              IF (ctypsub(ig_number_field(jf)) .EQ. 'NONSOLAR') THEN 
                  WRITE(UNIT = nulou,FMT = 3032) 
     $                cdqdt(ig_number_field(jf)),
     $                cfldcoa(ig_number_field(jf)), 
     $                cfldfin(ig_number_field(jf))
                ELSE IF (ctypsub(ig_number_field(jf)) .EQ. 'SOLAR') THEN
                  WRITE(UNIT = nulou,FMT = 3033)
     $                cfldfin(ig_number_field(jf)), 
     $                  cfldcoa(ig_number_field(jf))
              ENDIF
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'CHECKIN') THEN
                WRITE(UNIT = nulou,FMT = 3034) 
     $              ntinpflx(ig_number_field(jf))
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'CHECKOUT') THEN
                WRITE(UNIT = nulou,FMT = 3035) 
     $              ntoutflx(ig_number_field(jf)) 
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'GLORED') THEN
              WRITE(UNIT = nulou,FMT = 3036) 
     $              ntronca(ig_number_field(jf)) 
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'NOINTERP') THEN
                WRITE(UNIT = nulou,FMT = 3037)
            ELSE IF (canal(ja,ig_number_field(jf)) .EQ. 'LOCTRANS') THEN
               WRITE(UNIT = nulou,FMT = 3047) cl_print_trans
            ELSE 
              WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
              WRITE (UNIT = nulou,FMT = *)
     $            ' Type of analysis not implemented yet '
              WRITE (UNIT = nulou,FMT = *) 
     $            ' The analysis required in OASIS is :'
              WRITE (UNIT = nulou,FMT = *) ' canal = ', 
     $             canal(ja,ig_number_field(jf))
              WRITE (UNIT = nulou,FMT = *) 
     $            ' with ja = ', ja, ' jf = ', jf
              WRITE (UNIT = nulou,FMT = *) ' '
              CALL HALTE ('STOP in inipar')
          ENDIF
 320    CONTINUE
      ENDIF
      ENDIF
 310  CONTINUE
      ENDIF
C
C* Formats
C
 3001 FORMAT(//,15X,'  FIELD NUMBER ',I3)
 3002 FORMAT(15X,'  ************  ')
 3003 FORMAT(/,10X,'  Field parameters ')
 3004 FORMAT(10X,'  ****************  ',/)
 3005 FORMAT(/,10X,'  Input field symbolic name       = ',A8,
     $       /,10X,'  Output field symbolic name      = ',A8,
     $       /,10X,'  Field long name                 = ',
     $       /,18X,A53,
     $       /,10X,'  Field exchange frequency        = ',I8,
     $       /,10X,'  Model sequential index          = ',I2,
     $       /,10X,'  Field Lag                       = ',I8,
C     $       /,10X,'  Model delay flag                = ',I2,
C     $       /,10X,'  Extra time step flag            = ',I2,
     $       /,10X,'  Field I/O status                = ',A8,
     $       /,10X,'  Number of basic operations      = ',I4,
     $       /,10X,'  Parallel decomposition strategy = ',A8,/)
 3115 FORMAT(/,10X,'  Input field symbolic name       = ',A8,
     $       /,10X,'  Output field symbolic name      = ',A8,
     $       /,10X,'  Field long name                 = ',
     $       /,18X,A53,
     $       /,10X,'  Field exchange frequency        = ',I8,
     $       /,10X,'  Model sequential index          = ',I2,
C     $       /,10X,'  Model delay flag                = ',I2,
C     $       /,10X,'  Extra time step flag            = ',I2,
     $       /,10X,'  Field I/O status                = ',A8,
     $       /,10X,'  Number of basic operations      = ',I4,/)
 3116 FORMAT(/,10X,'  Input field symbolic name       = ',A8,
     $       /,10X,'  Output field symbolic name      = ',A8,
     $       /,10X,'  Field long name                 = ',
     $       /,18X,A53,
     $       /,10X,'  Field exchange frequency        = ',I8,
     $       /,10X,'  Local transformation            = ',A8,
     $       /,10X,'  Model sequential index          = ',I2,
     $       /,10X,'  Field Lag                       = ',I8, 
     $       /,10X,'  Field I/O status                = ',A8,
     $       /,10X,'  Number of basic operations      = ',I4,/)
 3117 FORMAT(/,10X,'  Restart file name               = ',A8,/)
 3118 FORMAT(/,10X,'  Input file name                 = ',A32,/)
 3121 FORMAT(/,10X,'  Input field symbolic name       = ',A8,
     $       /,10X,'  Output field symbolic name      = ',A8,
     $       /,10X,'  Field long name                 = ',
     $       /,18X,A53,
     $       /,10X,'  Field exchange frequency        = ',I8,
     $       /,10X,'  Local transformation            = ',A8,
     $       /,10X,'  Field I/O status                = ',A8,
     $       /,10X,'  Number of basic operations      = ',I4,/)
 3006 FORMAT(/,10X,'  Input file name                 = ',A8,
     $       /,10X,'  Output file name                = ',A8,/)
 3007 FORMAT(
     $       /,10X,'  Source grid periodicity type is      = ',A8,
     $       /,10X,'  Number of overlapped grid points is  = ',I2,
     $       /,10X,'  Target grid periodicity type is      = ',A8,
     $       /,10X,'  Number of overlapped grid points is  = ',I2,/)
 3008 FORMAT(/,10X,'  Source longitude file string    = ',A8,
     $       /,10X,'  Source latitude file string     = ',A8,
     $       /,10X,'  Source mask file string         = ',A8,
     $       /,10X,'  Source surface file string      = ',A8,
     $       /,10X,'  Target longitude file string    = ',A8,
     $       /,10X,'  Target latitude file string     = ',A8,
     $       /,10X,'  Target mask file string         = ',A8,
     $       /,10X,'  Target surface file string      = ',A8,/)
 3009 FORMAT(/,10X,'  ANALYSIS PARAMETERS ')
 3010 FORMAT(10X,'  ******************* ',/)
 3011 FORMAT(/,5X,'  ANALYSIS number ',I2,' is ',A8,
     $       /,5X,'  ***************  ',/)
 3012 FORMAT(5X,' Value for masked points is        = ',E15.6)
 3013 FORMAT(5X,' Grid mapping file = ',A8,' linked to unit = ',I2,
     $     /,5X,' Dataset identificator number      = ',I2,
     $     /,5X,' Maximum number of neighbors is    = ',I4)
 3014 FORMAT(5X,' Source grid latitude order is    = ',A8)
 3015 FORMAT(5X,' Source grid longitude order is     = ',A8)
 3016 FORMAT(5X,' Target grid latitude order is    = ',A8)
 3017 FORMAT(5X,' Target grid longitude order is     = ',A8)
 3018 FORMAT(5X,' Extrapolation method is           = ',A8,
     $     /,5X,' Number of neighbors used is       = ',I2)
 3019 FORMAT(5X,' Extrapolation file = ',A8,' linked to unit = ',I2,
     $     /,5X,' Dataset identificator number      = ',I2)
 3020 FORMAT(5X,' Interpolation method is           = ',A8,
     $     /,5X,' Source grid type is               = ',A8,
     $     /,5X,' Field type is                     = ',A8)
 3021 FORMAT(5X,' Pointer for ANAIS storage is      = ',I2,
     $     /,5X,' Maximum number of neighbors is    = ',I4,
     $     /,5X,' Write/Read flag for weights is    = ',I2)
 3022 FORMAT(5X,' Variance multiplicator for ANAISG = ',E15.6)
 3023 FORMAT(5X,' Data to fill up field is in file  = ',A8,
     $     /,5X,' Connected to logical unit number  = ',I2,
     $     /,5X,' Filling method to blend field is  = ',A8)
 3024 FORMAT(5X,' Flag for coasts mismatch is       = ',I2, 
     $     /,5X,' Name for flux correction field is = ',A8,
     $     /,5X,' It is written on logical unit     = ',I2)
 3025 FORMAT(5X,' Conservation method for field is  = ',A8)
 3026 FORMAT(5X,' Half number of latitudes for gaussian grid is = ',I3,
     $     /,5X,' Extrapolation flag is             = ',A8)
 3027 FORMAT(5X,' Field ',A8,' is multiplied by Cst = ',E15.6)
 3028 FORMAT(5X,' It is combined with N fields    N = ',I2)
 3029 FORMAT(5X,' Data file = ',A8,' linked to unit = ',I2)
 3030 FORMAT(5X,'   With field ',A8,'   coefficient = ',E15.6)
 3031 FORMAT(5X,' Subgrid data file = ',A8,' linked to unit = ',I2,
     $     /,5X,' Dataset identificator number      = ',I2,
     $     /,5X,' Maximum number of neighbors is    = ',I4,
     $     /,5X,' Type of subgrid interpolation is  = ',A8)
 3032 FORMAT(5X,' Subgrid variability is restored with addition of',
     $       /,5X,A8,' x (',A8,' - ',A8,')')
 3033 FORMAT(5X,' Subgrid variability is restored multiplying by',
     $       /,5X,'( 1 - ',A8,') / ( 1 - ',A8,')')
 3034 FORMAT(5X,' Integral calculation flag is =', I2)
 3035 FORMAT(5X,' Integral calculation flag is =', I2)
 3036 FORMAT(5X,' Half number of latitudes for gaussian grid is = ',I3)
 3037 FORMAT(5X,' No interpolation for this field ')
 3038 FORMAT(5X,' Write/Read flag for weights is    = ',I2,
     $     /,5X,' Dataset identificator number      = ',I2)
C 3039 FORMAT(/,5X,' No delay flag in namcouple for field', I3,
C     $     /,5X,' Default value DEL=0 will be used ')
C 3040 FORMAT(/,5X,' No extra timestep flag in namcouple for field', I3,
C     $     /,5X,' Default value XTS=0 will be used ') 
 3041 FORMAT(/,5X,' WARNING: Extra timestep flag > 1 for field', I3,
     $     /,5X,' XTS=1 will be used instead') 
 3042 FORMAT(5X,' Value for exported masked points is = ',E15.6)
 3043 FORMAT(/,5X,'No lag in namcouple for the field', I3,
     $    /,5X,' Default value LAG=0 will be used ')
 3044 FORMAT(/,5X,'The lag for the field ',I3,3X,'is : ',I8)
 3045 FORMAT(5X,' Remapping method is               = ',A8,
     $     /,5X,' Field type is                     = ',A8,
     $     /,5X,' Normalization option is           = ',A8,
     $     /,5X,' Seach restriction type is         = ',A8,
     $     /,5X,' Number of search bins is          = ',I4)
 3046 FORMAT(5X,' Order of remapping is             = ',A8)
 3047 FORMAT(5X,' Local transformation  = ',A8) 


C
C
C*    4. End of routine
C        --------------
C
      WRITE(UNIT = nulou,FMT = *) ' '
      WRITE(UNIT = nulou,FMT = *) 
     $    '          ---------- End of routine inipar ---------'
      CALL FLUSH (nulou)
      RETURN
C
C* Error branch output
C
 110  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *) 
     $    ' No active $JOBNAME data found in input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) 
     $    ' We STOP!!! Check the file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL HALTE ('STOP in inipar')
 130  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *)
     $    ' No active $NBMODEL data found in input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) 
     $    ' We STOP!!! Check the file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL HALTE ('STOP in inipar')
 170  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *) 
     $    ' No active $MACHINE data found in input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) 
     $    ' We STOP!!! Check the file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL HALTE ('STOP in inipar')
      
 181  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *)
     $    ' No active $CHATYPE data found in input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *)
     $    ' We STOP!!! Check the file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL HALTE ('STOP inipar')
 191  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *)
     $    ' No active $RUNTIME data found in input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) 
     $    ' We STOP!!! Check the file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL HALTE ('STOP in inipar')
 193  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *)
     $    ' No active $INIDATE data found in input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) 
     $    ' We STOP!!! Check the file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL HALTE ('STOP in inipar')
 195  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *)
     $    ' No active $SEQMODE data found in input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) 
     $    ' We STOP!!! Check the file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL HALTE ('STOP in inipar')
 197  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *)
     $    ' No active $MODINFO data found in input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) 
     $    ' We STOP!!! Check the file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL HALTE ('STOP in inipar')
 199  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *)
     $    ' No active $NLOGPRT found in input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) 
     $    ' We STOP!!! Check the file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL HALTE ('STOP in inipar')
 201  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *)
     $    ' No active $CALTYPE found in input file namcouple'
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
      CALL HALTE ('STOP in inipar')
 230  CONTINUE
      WRITE (UNIT = nulou,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulou,FMT = *) 
     $    ' No active $STRING data found in input file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) ' '
      WRITE (UNIT = nulou,FMT = *) 
     $    ' We STOP!!! Check the file namcouple'
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL HALTE ('STOP in inipar')
 231  CONTINUE
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL prtout ('ERROR in namcouple for field', jf, 1)
      WRITE (UNIT = nulou,FMT = *) 
     $             'NO index of sequential position and $SEQMODE > 1'
      CALL halte('STOP in inipar.f')
 232  CONTINUE
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL prtout ('ERROR in namcouple for field', jf, 1)
      WRITE (UNIT = nulou,FMT = *) 
     $             'Index of sequential position greater than $SEQMODE'
      CALL halte('STOP in inipar.f') 
 233  CONTINUE
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL prtout ('ERROR in namcouple for field', jf, 1)
      WRITE (UNIT = nulou,FMT = *) 
     $'Check the 2nd line for either the index of sequential position, 
     $the delay flag, or the extra timestep flag.'
      CALL halte('STOP in inipar.f')
 234  CONTINUE
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL prtout ('ERROR in namcouple for field', jf, 1)
      WRITE (UNIT = nulou,FMT = *) 
     $             'Index of sequential position equals 0'
      WRITE (UNIT = nulou,FMT = *) 
     $             '(Should be 1 -default value- IF $SEQMODE=1)'
      CALL halte('STOP in inipar.f') 
 235  CONTINUE
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL prtout ('ERROR in namcouple for field', jf, 1)
      WRITE (UNIT = nulou,FMT = *) 
     $      'An input line with integral calculation flag' 
      WRITE (UNIT = nulou,FMT = *) 
     $      '("INT=0" or "INT=1")'
      WRITE (UNIT = nulou,FMT = *) 
     $      'is now required for analysis CHECKIN or CHECKOUT'
      CALL halte('STOP in inipar.f') 
 236  CONTINUE
      WRITE (UNIT = nulou,FMT = *) ' '
      CALL prtout ('ERROR in namcouple for field', jf, 1)
      WRITE (UNIT = nulou,FMT = *) 
     $     'The coupling period must not be 0 !'
      WRITE (UNIT = nulou,FMT = *) 
     $     'If you do not want to exchange this field at all'
       WRITE (UNIT = nulou,FMT = *) 
     $     'give a coupling period longer than the total run time.'
       CALL halte('STOP in inipar.f') 
      END
