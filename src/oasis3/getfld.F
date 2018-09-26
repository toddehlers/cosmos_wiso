      SUBROUTINE getfld (kindex, kfield, kiter)
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 1 *
C               * -------------     ------- *
C               *****************************
C
C**** *getfld* - reading routine
C
C
C     Purpose:
C     -------
C     Read coupling fields for iteration kiter
C
C**   Interface:
C     ---------
C       *CALL*  *getfld (kindex, kfield, kiter)*
C
C     Input:
C     -----
C                kindex : current active fields index array
C                kfield : current active fields total number
C                kiter  : iteration number
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
C     PIPE_Recv, CLIM_Import, SVIPC_read
C
C     Reference:
C     ---------
C     See OASIS manual (1995)
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       2.0       L. Terray      95/09/01  created
C       2.1       L. Terray      96/08/05  modified: correct printing of
C                                          diagnostic (variable clname)
C       2.2       S. Valcke      97/08/22  added: introduction of SVIPC
C       2.2       L.Terray       97/09/20  added: extra test in delay mode
C                                          + general cleaning + test on
C                                          info mode
C       2.3       L. Terray      99/03/01  modified: reading in delay mode
C                                          (negative values mean one must
C                                           read files, 0 means channel)
C       2.3       S. Valcke      99/04/30  added: printing levels
C       2.3       L. Terray      99/09/15  added: GMEM branch
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* -----------------Include files and USE of modules---------------------------
C
      USE mod_kinds_oasis
#if defined use_comm_MPI1 || defined use_comm_MPI2 || defined use_comm_GSIP || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_GSIP && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
      USE mod_clim
#endif
      USE mod_parameter
      USE mod_string
      USE mod_analysis
      USE mod_memory
      USE mod_sipc
      USE mod_unitncdf
      USE mod_experiment
      USE mod_timestep
      USE mod_unit
      USE mod_hardware
      USE mod_label
      USE mod_printing
#ifdef use_netCDF
#include <netcdf.inc>
#endif
C
C* ---------------------------- Argument declarations -------------------
C
      INTEGER (kind=ip_intwp_p) kindex(kfield)
C
C* ---------------------------- Local declarations ----------------------
C
      INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: info, 
     $    iflag
      CHARACTER*8 clname, clfic
      CHARACTER*53 clabel
      INTEGER (kind=ip_intwp_p) itime(3)
      CHARACTER*3 cljobnam
      INTEGER (kind=ip_intwp_p) il_varid, il_dimid, il_dimlen
      INTEGER (kind=ip_intwp_p) il_nvar, il_natt, il_att
      INTEGER (kind=ip_intwp_p) ist(3), icnt(3)
      INTEGER (kind=ip_intwp_p) il_indice_since, il_cf_units
      LOGICAL l_unit_found
      CHARACTER(len=60) cl_attname
      CHARACTER(len=17) cl_units, cl_units_cf
      CHARACTER(len=260) cl_atttext
      INTEGER (kind=ip_intwp_p) , parameter :: ip_nb_cf_units = 63
      CHARACTER(len=17), dimension(ip_nb_cf_units) :: 
     $         cl_cf_time_units = 
     $         (/'second           ', 'minute           ',
     $           'hour             ', 'day              ',
     $           's                ', 'sec              ',
     $           'shake            ', 'sidereal_day     ',
     $           'sidereal_hour    ', 'sidereal_minute  ',
     $           'sidereal_second  ', 'sidereal_year    ',
     $           'tropical_year    ', 'lunar_month      ',
     $           'common_year      ', 'leap_year        ',
     $           'Julian_year      ', 'Gregorian_year   ',
     $           'sidereal_month   ', 'tropical_month   ',
     $           'd                ', 'min              ',
     $           'hr               ', 'h                ',
     $           'fortnight        ', 'week             ',
     $           'jiffy            ', 'jiffies          ',
     $           'year             ', 'yr               ',
     $           'a                ', 'eon              ',
     $           'month            ',
     $           'seconds          ', 'minutes          ',
     $           'hours            ', 'days             ',
     $           'secs             ',
     $           'shakes           ', 'sidereal_days    ',
     $           'sidereal_hours   ', 'sidereal_minutes ',
     $           'sidereal_seconds ', 'sidereal_years   ',
     $           'tropical_years   ', 'lunar_months     ',
     $           'common_years     ', 'leap_years       ',
     $           'Julian_years     ', 'Gregorian_years  ',
     $           'sidereal_months  ', 'tropical_months  ',
     $           'ds               ', 'mins             ',
     $           'hrs              ', 'hs               ',
     $           'fortnights       ', 'weeks            ',
     $           'years            ', 'yrs              ',
     $           'as               ', 'eons             ',
     $           'months           '/)
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Allocation and initialization
C        -----------------------------
C
      IF (nlogprt .GE. 2) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '           ROUTINE getfld  -  Level 1'
          WRITE (UNIT = nulou,FMT = *) 
     $    '           **************     *******'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' Get coupling fields'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
      ENDIF
      imrca = 0
      imrcb = 0
      itime(:)=0
      ist(:)=0
      icnt(:)=0
      ALLOCATE(iflag(ig_nfield))
      ALLOCATE(info(ig_nfield))
      CALL izero (iflag, ig_nfield)
      CALL izero (info, ig_nfield)
C
C*    2. Loop on active fields for iteration kiter
C        -----------------------------------------
!$omp parallel do default (shared)
!$omp+ private(jf,iloc,iadrold,isizold,ilabel)
!$omp+ private(clabel,clname,clfic,iunit)
!$omp+ private(iseqn,il_varid,n_reaty)

C
      DO 210 jf = 1, kfield
C
C* Assign local variables
C
         iloc = kindex(jf)
         iadrold = nadrold(iloc)
         isizold = nsizold(iloc)
         ilabel = numlab(iloc)
         clabel = cfldlab(ilabel)
         clname = cnaminp(iloc)
         clfic = cficinp(iloc)
         iunit = nluinp(iloc)
         iseqn = nseqn(iloc)
         ilagn = nlagn(iloc)
C
C* Print field name
C
         IF (nlogprt .GE. 1) THEN
            CALL prcout('Reading of field : ', clname, 2)
         ENDIF
         IF (nlogprt .GE. 2) THEN
            CALL prcout('Field long name : ', clabel, 2)
         ENDIF
C
C* - Get coupling fields
C
C* Specific treatment for first iteration 
C
#if defined use_comm_NONE || defined use_comm_SIPC
         ilagn = 1
#endif
         IF (kiter .EQ. 0 .AND. ilagn .GT. 0) THEN
C
C* Reconnect data file to appropriate unit
C
            IF (lncdfrst) THEN
               ist(1)=1 
               ist(2)=1
               ist(3)=kiter+1
               icnt(1)=nlonbf(iloc) 
               icnt(2)=nlatbf(iloc)
               icnt(3)=1
C  
#ifdef use_netCDF   
               CALL hdlerr
     $          (NF_INQ_VARID(nc_inpid(iloc),clname,il_varid),'getfld')
#if defined use_comm_NONE
               nc_invartimeid=0
               istat = NF_INQ_VARID(nc_inpid(iloc),'time_counter',
     $              nc_invartimeid ) 
               IF (istat .ne. NF_NOERR) then
                   istat = NF_INQ_VARID(nc_inpid(iloc),'time',
     $                                  nc_invartimeid)
               ENDIF
               IF (istat .ne. NF_NOERR) then
                   istat = NF_INQ_VARID(nc_inpid(iloc),'TIME_COUNTER',
     $                                  nc_invartimeid)
               ENDIF
               IF (istat .ne. NF_NOERR) then
                   istat = NF_INQ_VARID(nc_inpid(iloc),'TIME',
     $                                  nc_invartimeid)
               ENDIF
               IF (istat .ne. NF_NOERR) then
                   istat = NF_INQ_UNLIMDIM(nc_inpid(iloc),
     $                     il_dimid)
                   IF (istat .eq. NF_NOERR) then
                     istat = NF_INQ_DIMNAME(nc_inpid(iloc),il_dimid,
     $                       nc_invartime_name)
                       IF (istat .eq. NF_NOERR) then
                           istat = NF_INQ_VARID(nc_inpid(iloc),
     $                             nc_invartime_name, nc_invartimeid)
                       ENDIF
                   ENDIF
               ENDIF
               l_unit_found=.FALSE.
               IF (istat .ne. NF_NOERR) then
                   CALL hdlerr(NF_INQ_NVARS(nc_inpid(iloc),il_nvar),
     $                         'getfld')
                  DO nc_invartimeid=1,il_nvar
                    CALL hdlerr(NF_INQ_VARNATTS(nc_inpid(iloc),
     $                             nc_invartimeid, il_natt),'getfld')
                    DO il_att=1,il_natt
                       CALL hdlerr(NF_INQ_ATTNAME(nc_inpid(iloc),
     $                          nc_invartimeid, il_att, cl_attname),
     $                          'getfld')
                       IF ( cl_attname .eq. 'units' ) THEN
                         CALL hdlerr(NF_GET_ATT_TEXT(
     $                          nc_inpid(iloc),nc_invartimeid ,'units',
     $                          cl_atttext),'getfld')
                         il_indice_since=INDEX(cl_atttext,'since')
                         IF ( il_indice_since .ne. 0 ) THEN
                            cl_units=cl_atttext(1:il_indice_since-1)
                            cl_units=ADJUSTL(cl_units)
                            cl_units=TRIM(cl_units)
                            DO il_cf_units =1, ip_nb_cf_units
                               cl_units_cf=cl_cf_time_units(il_cf_units)
                              IF ( cl_units .EQ. cl_units_cf )
     $                           l_unit_found=.TRUE.
                              CALL case_change(cl_units_cf,1)
                              IF ( cl_units .EQ. cl_units_cf)
     $                           l_unit_found=.TRUE.
                            ENDDO
                            IF ( l_unit_found )  THEN
                              istat = NF_NOERR
                              GO TO 230
                            ENDIF
                         ENDIF
                       ENDIF
                    ENDDO
                  ENDDO
 230              CONTINUE
               ENDIF
C
               IF (istat .ne. NF_NOERR) then
                  WRITE(nulou,*) 'WARNING: Time variable not found'
                  WRITE(nulou,*) ''
                  IF ( nitfn .gt. 0) THEN
                     WRITE(nulou,*) 'Check coherence between'
                     WRITE(nulou,*) 'namcouple timestep number'
                     WRITE(nulou,*) 'and input file dimension'
                     WRITE (nulou, *)
                     CALL HALTE('STOP in getfld')
                  END IF
             WRITE(nulou,*) 'Interpolation only for one time iteration'
                  nc_invartime_name = 'time_counter'
               ELSE
                  CALL hdlerr
     $              (NF_INQ_VARNAME(nc_inpid(iloc),nc_invartimeid,
     $               nc_invartime_name),'getfld')
C
C                 Check coherence 
C                 between namcouple timestep number
C                 and input file dimension
C    
                  CALL hdlerr
     $              (NF_INQ_VARDIMID(nc_inpid(iloc),nc_invartimeid,
     $               il_dimid),'getfld')
                  CALL hdlerr
     $              (NF_INQ_DIMLEN(nc_inpid(iloc),il_dimid,
     $               il_dimlen),'getfld')
                  IF ( il_dimlen .le. nitfn ) then
                      WRITE(nulou,*) 'Incoherence between'
                      WRITE(nulou,*) 'namcouple timestep number'
                      WRITE(nulou,*) 'and input file dimension'
                      WRITE(nulou,*) ''
                      CALL HALTE('STOP in getfld')
                  ENDIF
              ENDIF
#else
               nc_invartime_name = 'time_counter'
#endif
               CALL hdlerr
     $              (NF_INQ_VARTYPE(nc_inpid(iloc),il_varid,n_reaty),
     $              'getfld')
               IF (n_reaty .eq. NF_FLOAT) THEN
c                  IF (kind(rl_testvar). ne. 4) THEN
                  IF ( ip_realwp_p .ne. ip_single_p ) THEN
                     WRITE (nulou, *)
     $         'Incoherence: Oasis is compiled in double precision'
                     WRITE (nulou, *)
     $      'but Restart auxiliary file contains single precision REAL'
                     CALL HALTE('STOP in getfld')
                  ENDIF
                  CALL hdlerr(NF_GET_VARA_REAL(nc_inpid(iloc), 
     $                 il_varid, ist, icnt, 
     $                fldold(iadrold:iadrold+isizold-1)), 'getfld')
C
#if defined use_comm_NONE
                  IF (nc_invartimeid .ne. 0) then
                      istat = NF_GET_VARA_REAL(nc_inpid(iloc),
     $                     nc_invartimeid, ist(3), icnt(3), rtime_val)
                      IF (istat .eq. NF_NOERR) then
                         istatus=NF_GET_VARA_REAL(nc_inpid(iloc),
     $                    nc_invartimeid, ist(3), icnt(3), rtime_val)
                      ENDIF
                  ENDIF
C
#endif
               ELSE IF (n_reaty .eq. NF_DOUBLE) THEN
c                  IF (kind(rl_testvar). ne. 8) THEN
                  IF ( ip_realwp_p .ne. ip_double_p ) THEN
                     WRITE (nulou, *)
     $          'Incoherence: Oasis is compiled in single precision'
                     WRITE (nulou, *)
     $    'but Restart auxiliary file contains double precision REAL'
                     CALL HALTE('STOP in getfld')
                  ENDIF               
                  CALL hdlerr(NF_GET_VARA_DOUBLE(nc_inpid(iloc), 
     $                 il_varid, ist, icnt, 
     $                fldold(iadrold:iadrold+isizold-1) ), 'getfld')
C
#if defined use_comm_NONE
                  IF (nc_invartimeid .ne. 0) then
                      istat = NF_GET_VARA_DOUBLE(nc_inpid(iloc),
     $                     nc_invartimeid, ist(3), icnt(3), rtime_val)
                      IF (istat .eq. NF_NOERR) then
                         istatus=NF_GET_VARA_DOUBLE(nc_inpid(iloc),
     $                    nc_invartimeid, ist(3), icnt(3), rtime_val)
                      ENDIF
                  ENDIF
C
#endif
               ENDIF
#endif
           ELSE
cSL
C
C* Reconnect data file to appropriate unit; this first close is necessary !
C
                    CLOSE(iunit)
                    OPEN (UNIT=iunit,FILE=clfic,STATUS ='UNKNOWN',
     $                      FORM ='UNFORMATTED',IOSTAT = iost)
                    IF (iost .eq. 0) then
                        IF (nlogprt .GE. 2) THEN
                            WRITE(UNIT = nulou, FMT = 2100)iunit,clfic
                        ENDIF
                    ELSE
                        CALL prtout
     $                  ('No binary transfert file for field ',iloc,1)
                        CALL HALTE('STOP in getfld') 
                    ENDIF
cSL
               IF (lmodinf) THEN 
                  CALL locreadh (clname, cljobnam, itime,  
     $                 fldold(iadrold), isizold, iunit, iflag(iloc)) 
               ELSE 
                  CALL locread (clname, fldold(iadrold), isizold,
     $                 iunit, iflag(iloc))
               ENDIF
               IF (iflag(iloc) .NE. 0) THEN
                  CALL prcout ('WARNING: problem in reading field',
     $                 clname, 1)
                  CALL prtout ('Error reading logical unit', iunit, 2)
                  CALL HALTE('STOP in getfld')
               ENDIF
C
           ENDIF
C
C* Next iterations
C
         ELSE
C
C* PIPE and NONE case
C
#ifdef use_comm_PIPE
C
C* Wait for message
C
            CALL PIPE_Recv (clname, kiter)
#endif
#if defined use_comm_PIPE || defined use_comm_NONE
C
C* Reconnect data file to appropriate unit IF not already done
C
            IF (lncdfrst) THEN
               ist(1)=1 
               ist(2)=1
               ist(3)=1+kiter
               icnt(1)=nlonbf(iloc)
               icnt(2)=nlatbf(iloc)
               icnt(3)=1
#ifdef use_netCDF
               CALL hdlerr(NF_OPEN(clfic,NF_NOWRITE,nc_inpid(iloc)),
     $             'getfld')
               istatus=NF_INQ_VARID(nc_inpid(iloc),clname,il_varid)
               
               istatus=NF_INQ_VARID(nc_inpid(iloc),nc_invartime_name,
     $                        nc_invartimeid)
               IF (n_reaty .eq. NF_FLOAT) THEN
                  CALL hdlerr(NF_GET_VARA_REAL(nc_inpid(iloc), 
     $                 il_varid,ist,icnt,
     $                 fldold(iadrold:iadrold+isizold-1)),'getfld') 
                  IF (nc_invartimeid .ne. 0) 
     $            CALL hdlerr(NF_GET_VARA_REAL(nc_inpid(iloc),
     $                 nc_invartimeid, ist(3), icnt(3), rtime_val)
     $                 , 'getfld') 
               ELSE IF (n_reaty .eq. NF_DOUBLE) THEN
                  CALL hdlerr(NF_GET_VARA_DOUBLE(nc_inpid(iloc), 
     $                 il_varid,ist,icnt,
     $                 fldold(iadrold:iadrold+isizold-1)),'getfld')
                  IF (nc_invartimeid .ne. 0) 
     $               CALL hdlerr(NF_GET_VARA_DOUBLE(nc_inpid(iloc),
     $                 nc_invartimeid, ist(3), icnt(3), rtime_val)
     $                 , 'getfld')
               ENDIF
               CALL hdlerr(NF_CLOSE(nc_inpid(iloc)), 'getfld')
#endif
            ELSE
               OPEN (UNIT=iunit,FILE=clfic,STATUS ='UNKNOWN',
     $              FORM ='UNFORMATTED',IOSTAT = iost)
               IF (iost .eq. 0) then
                  IF (nlogprt .GE. 2) THEN
                     WRITE(UNIT = nulou, FMT = 2100)iunit,clfic
                  ENDIF
               ELSE
                  CALL prtout
     $                 ('No binary transfert file for field ',iloc,1)
                  CALL HALTE('STOP in getfld') 
               ENDIF
               IF (lmodinf) THEN 
                  CALL locreadh (clname, cljobnam, itime,  
     $                 fldold(iadrold),isizold,iunit,iflag(iloc)) 
               ELSE
                  CALL locread (clname, fldold(iadrold), 
     $                 isizold, iunit, iflag(iloc))
               ENDIF 
               IF (iflag(iloc) .NE. 0) THEN 
                  clname = cnaminp(iloc)
                  CALL prcout('Problem reading field',clname,1)
                  CALL prtout ('Error reading lunit', iunit, 2)
                  CALL HALTE ('STOP in getfld')
               ENDIF
               CLOSE (UNIT = iunit, IOSTAT = ios)
               IF (ios .EQ. 0) THEN
                  IF (nlogprt .GE. 2) THEN
                     WRITE(UNIT=nulou,FMT=2200) iunit, clfic
                  ENDIF
               ELSE
                  CALL prtout('Problem in closing unit',iunit,2)
                  CALL prtout('Error message nbrr is= ', ios, 2)
                  CALL HALTE('STOP in getfld')
               ENDIF 
            ENDIF
C
C* SIPC or GMEM case
C
#elif defined use_comm_SIPC || defined use_comm_GMEM
C
            ipbytecha=kind('A')
            ipbyteint=kind(itime)
            ipbyterea=kind(fldold)
C
C* Read encapsulated infos in field-specific shared memory pool
C  (experiment name, initial date, time since start, iteration number)
C
            IF (lmodinf) THEN 
               isizeinp = 3*ipbytecha
               CALL SVIPC_read(mpoolidin(iloc), cljobnam 
     $              , isizeinp, imrca)
               isizeinp = 3*ipbyteint
               CALL SVIPC_read(mpoolidin(iloc), itime 
     $              , isizeinp, imrcb)
C
C* Find error if any
C
               IF (imrca .LT. 0 .OR. imrcb .LT. 0) THEN
                  CALL prcout
     $   ('Problem in reading encapsulated infos for field', clname, 1)
                  CALL HALTE ('STOP in getfld')
               ELSE
                  CALL prcout
     $   ('Read encapsulated infos in SHM pool for field', clname, 1)
               ENDIF 
            ENDIF
C
C* Read part of macro array in field-specific shared memory pool
C
            isizeinp = isizold * ipbyterea
            WRITE(UNIT = nulou, FMT = *) 
     $           'Reading field data from pool = ',mpoolidin(iloc) 
            CALL SVIPC_read(mpoolidin(iloc),
     $           fldold(iadrold), isizeinp, imrc)
C
C* Find error if any
C
            IF (imrc .LT. 0) THEN
               CALL prcout
     $           ('Problem in reading the field in SHM pool:',clname,1)
               CALL HALTE ('STOP in getfld')
            ELSE
               CALL prcout
     $            ('Read field in SHM pool:', clname, 1)
            ENDIF
C
C* CLIM case
C
#elif defined use_comm_MPI1 || defined use_comm_MPI2 || defined use_comm_GSIP || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_GSIP && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE) 
!$omp critical
            CALL CLIM_Import 
     $           (ig_portin_id(iloc), kiter*nstep, fldold(iadrold), 
     $           info(iloc))
!$omp end critical
            IF (info(iloc) .NE. CLIM_Ok) THEN 
               CALL prcout('Problem in reading field on port',
     $              clname, 1)
               CALL prtout('error code number', info(iloc), 2)
               CALL HALTE ('STOP in getfld')
            ENDIF
#else
            CALL prcout 
     $           ('Wrong channel option for field', clname, 1)
#endif
         IF (lmodinf .AND. nlogprt .GE. 2) THEN 
            WRITE (UNIT = nulou,FMT = *) 
     $           ' Encapsulated data for current field is :' 
            CALL prcout ('Experiment name', cljobnam, 1)
            CALL prtout ('Initial date', itime(1), 2)
            CALL prtout ('Iteration number', itime(2), 2)
            CALL prtout ('Time since start', itime(3), 2)
         ENDIF 
      ENDIF 
 210  CONTINUE
C
C* Close units or netcdf restart files if first iteration
C
      IF (kiter .EQ. 0) THEN
          DO 220 jf = 1, kfield
            iloc = kindex(jf)
            ilagn = ig_lag(iloc)
            IF (ilagn .GT. 0) THEN
            isamefic=0
            DO 225 jj = 1, jf-1
              IF (nluinp(jf) .eq. nluinp(jj)) isamefic=isamefic+1
 225        END DO
            IF (isamefic .lt. 1) THEN
#ifdef use_netCDF
                IF (lncdfrst) THEN
                    CALL hdlerr(NF_CLOSE(nc_inpid(jf)), 'getfld')
                ELSE
#endif
                    CLOSE (UNIT = nluinp(jf), IOSTAT = ios)
                    IF (ios .NE. 0) THEN
                        CALL prtout('Problem in closing unit',iunit,2)
                        CALL prtout('Error message number is = ',ios,2)
                        CALL HALTE('STOP in getfld')
                    ELSE
                        IF (nlogprt .GE. 2) THEN
                            WRITE(UNIT=nulou, FMT=2200) iunit, clfic
                        ENDIF
                    ENDIF
#ifdef use_netCDF
                ENDIF
#endif
            ENDIF
            ENDIF
 220      CONTINUE
      ENDIF
C
C* Formats
C
 2100 FORMAT(/,5X,' Unit ',I2,' has been connected to file ',A8)
 2200 FORMAT(/,5X,' Unit ',I2,' has been disconnected from file ',A8)
C
C
C*    3. End of routine
C        --------------
C
      IF (nlogprt .GE. 2) THEN 
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '          --------- End of routine getfld ---------'
          CALL FLUSH (nulou)
      ENDIF
      DEALLOCATE(iflag)
      DEALLOCATE(info)
      RETURN
      END
C
C====================================================================
C
          subroutine case_change( string, type )
C
C Purpose
C _______
C
C   This converts each lower case alphabetic letter in STRING to upper
C   case, or vice versa.
C
C
C Arguments
C _________
C
C   string  (INPUT/OUTPUT) character(len=*)
C           The string to convert
C
C   TYPE    (INPUT/OUTPUT) integer(i1b)
C           Define the conversion. Specifically,
C
C               If TYPE = 1 = TOUPPER,    conversion is lower to upper
C               If TYPE = 2 = TOLOWER,    conversion is upper to lower
C               If TYPE = 3 = CAPITALIZE, use upper for first letter; lower for rest
C
C           Definitions of TOUPPER, TOLOWER and CAPITALIZE may be obtained from
C           the host module Strings.
C
C
C Further Details
C _______________
C
C   All non-alphabetic characters are left unchanged.
C   It uses the underlying machine collating sequence.
C
C     History:
C     -------
C       Version   Programmer         Date      Description
C       -------   ----------         ----      -----------
C       3.0       P. Terray (LODyC)  03/12/19  created
C
C _________________________________________________________________________________
C
C
C
C SPECIFICATIONS FOR ARGUMENTS
C ____________________________
C
          character(len=*), intent(inout) :: string
C
          integer,     intent(in)    :: type
C
C
C SPECIFICATIONS FOR LOCAL VARIABLES
C __________________________________
C
         integer :: i
         integer, parameter :: shift = ichar( 'a' ) - ichar( 'A' )
C
C
C EXECUTABLE STATEMENTS
C _____________________
C
        if ( type==1 ) then
            do i = 1, len( string )
                select case( string(i:i) )
                  case( 'a':'i' , 'j':'r' , 's':'z' )
                    string(i:i) = char( ichar(string(i:i)) - shift  )
                end select
            end do
        else
            do i = 1, len( string )
                select case( string(i:i) )
                  case( 'A':'I' , 'J':'R' , 'S':'Z' )
                    string(i:i) = char( ichar(string(i:i)) + shift  )
                end select
            end do
            if ( type==3 ) then
                select case( string(1:1) )
                  case( 'a':'i' , 'j':'r' , 's':'z' )
                    string(1:1) = char( ichar(string(1:1)) - shift  )
                end select
            end if
        end if
        return
C
C
C END OF SUBROUTINE case_change
C _____________________________
C
        end
C

