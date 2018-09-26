      SUBROUTINE givfld (kindex, kfield, kiter)
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 1 *
C               * -------------     ------- *
C               *****************************
C
C**** *givfld* - writing routine
C
C
C     Purpose:
C     -------
C     Write out coupling fields for iteration kiter
C
C**   Interface:
C     ---------
C       *CALL*  *givfld (kindex, kfield, kiter)*
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
C     PIPE_Send, CLIM_Export, SVIPC_write
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
C       2.1       L. Terray      96/08/07  modified: addition of cstate
C                                          to prevent field transfer
C       2.2       S. Valcke      97/08/22  added: introduction of SVIPC
C       2.2       L. Terray      97/12/14  added: test on mode info +
C                                          general cleaning
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
      USE mod_calendar
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
      INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: idimlon, 
     $    idimlat, ivarid
      INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: ilongrd, 
     $    ilatgrd
      CHARACTER*8 clname, clfic, clstat, clwork, clstrg
      CHARACTER*53 clabel
      CHARACTER*32 cl_att_name
      INTEGER (kind=ip_intwp_p) itime(3)
      INTEGER (kind=ip_intwp_p) ist(3), icnt(3), ivardim(3)
      INTEGER (kind=ip_intwp_p) idimtime, il_timevarid
      INTEGER (kind=ip_intwp_p) il_nb_att, il_ind_att
      CHARACTER*4 clon, clat
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
     $    '           ROUTINE givfld  -  Level 1'
          WRITE (UNIT = nulou,FMT = *) 
     $    '           **************     *******'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' Give coupling fields'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
      ENDIF
      infos = CLIM_Ok
      itime(:)=0
      ist(:)=0
      icnt(:)=0
      iflags = 0
      istop = 0
      imrca = 0
      imrcb = 0
      ALLOCATE (iflag(ig_nfield))
      ALLOCATE (info(ig_nfield))
      ALLOCATE (idimlon(ig_nfield))
      ALLOCATE (idimlat(ig_nfield))
      ALLOCATE (ivarid(ig_nfield))
      ALLOCATE (ilongrd(ig_nfield))
      ALLOCATE (ilatgrd(ig_nfield))      
      CALL izero (iflag, ig_nfield)
      CALL izero (info, ig_nfield)
      idimlon(:)=0
      idimlat(:)=0
      ivarid(:)=0
      ilongrd(:)=0
      ilatgrd(:)=0
C
C
C*    1.5 Create binary or netcdf output files for PIPE or NONE technique
C        ----------------------------------------------------------------
C
      IF (kiter .eq. 0 .or. cchan .eq. 'PIPE') THEN
      IF (cchan .eq. 'PIPE' .or. cchan .eq. 'NONE' ) THEN
         DO 120 jf = 1, kfield
           iloc = kindex(jf)
           isamefic=0
           DO 125 jj = 1, jf-1
             ilocp=kindex(jj)
             IF (nluout(iloc) .eq. nluout(ilocp)) THEN
                 isamefic=isamefic+1
                 nc_outid(jf) =  nc_outid(jj)
             ENDIF
 125       END DO
           IF (isamefic .lt. 1) THEN
               iunit = nluout(iloc)
               clfic = cficout(iloc)
#ifdef use_netCDF
               IF (lncdfrst) THEN
                   istatus=NF_CREATE(clfic, NF_CLOBBER, nc_outid(jf))
                   IF (istatus .ne. NF_NOERR) THEN
                      CALL prtout
     $           ('Cannot create netcdf output file for field ',iloc,1)
                      CALL HALTE('STOP in givfld') 
                   ENDIF
C
C                  For each new file created define all possible 
C                  dimensions and coordinate variables 
C
                   DO 130 jfint=1, kfield
                     ilocint = kindex(jfint)
                     isamedimlon = 0
                     ilon = 0
                     isamedimlat = 0
                     ilat = 0
                     isamegrd = 0
                     igrd = 0
C
C                    Test if same output file 
                     IF (nluout(ilocint) .eq. iunit) THEN
C
C                        Loop on all previous fields
                         DO 135 jfintp=1, jfint-1
                           ilocintp=kindex(jfintp)
C
C                          Test if same x dimension
                           IF(nluout(ilocint) .eq. nluout(ilocintp)
     $                .and. nlonaf(ilocint) .eq. nlonaf(ilocintp)) THEN
                               isamedimlon = isamedimlon + 1
                               idimlon(jfint) = idimlon(jfintp)
                           ENDIF
C
C                          Test if same y dimension 
                           IF(nluout(ilocint) .eq. nluout(ilocintp)
     $                .and. nlataf(ilocint) .eq. nlataf(ilocintp)) THEN
                               isamedimlat = isamedimlat + 1
                               idimlat(jfint) = idimlat(jfintp)
                           ENDIF
C
C                          Test IF same grid
                           IF(nluout(ilocint).eq.nluout(ilocintp) .and.
     $                     cficaf(ilocint) .eq. cficaf(ilocintp)) THEN
                               isamegrd = isamegrd + 1
                               ilongrd(jfint) = ilongrd(jfintp)
                               ilatgrd(jfint) = ilatgrd(jfintp)
                           ENDIF
 135                     CONTINUE
C
C                        Define x dimension if not already defined
                         IF (isamedimlon .lt. 1) THEN
                             ilon=ilon+1
                             clon='lon'//char(ilon+48)
                             CALL hdlerr(NF_DEF_DIM(nc_outid(jf), clon,
     $                       nlonaf(ilocint), idimlon(jfint)),'givfld')
                         ENDIF
C
C                        Define y dimension if not already defined
                         IF (isamedimlat .lt. 1) THEN
                             ilat=ilat+1
                             clat='lat'//char(ilon+48)
                             CALL hdlerr(NF_DEF_DIM(nc_outid(jf), clat,
     $                       nlataf(ilocint), idimlat(jfint)),'givfld')
                         ENDIF
C
C                        Define t dimension and variable at first 
C                        time step if necessary
                            CALL hdlerr(NF_DEF_DIM(nc_outid(jf),
     $                           nc_invartime_name,NF_UNLIMITED,
     $                           idimtime),'givfld')
                            CALL hdlerr(NF_DEF_VAR(nc_outid(jf),
     $                           nc_invartime_name,n_reaty,1,idimtime,
     $                           il_timevarid),'givfld')
                            icount = 80
C
C                        Copy time variable attributes
C
                         IF ( nitfn .gt. 0) THEN
                            CALL hdlerr(NF_OPEN(cficinp(1),NF_NOWRITE,
     $                           nc_inpid(1)),'givfld')
                            CALL hdlerr(NF_INQ_VARNATTS(nc_inpid(1),
     $                           nc_invartimeid, il_nb_att), 'givfld')
                            DO il_ind_att=1, il_nb_att
                                CALL hdlerr(NF_INQ_ATTNAME(nc_inpid(1),
     $                               nc_invartimeid, il_ind_att,
     $                               cl_att_name), 'givfld')
                                CALL hdlerr(NF_COPY_ATT(nc_inpid(1),
     $                               nc_invartimeid, cl_att_name,
     $                               nc_outid(jf), il_timevarid)
     $                               , 'givfld')
                            END DO
                            CALL hdlerr(NF_CLOSE(nc_inpid(1)))
                         ENDIF
C
C                        Create coord. variable IF not already created
                         IF (isamegrd .lt. 1) THEN
                             ivardim(1)=idimlon(jfint)
                             ivardim(2)=idimlat(jfint)
                             clwork = cficaf(jfint)
                             icount = ilenstr(clwork,jpeight)
                             clstrg = clwork(1:icount)//cglonsuf
                             CALL hdlerr(NF_DEF_VAR(nc_outid(jf),
     $               clstrg,n_reaty,2,ivardim,ilongrd(jfint)),'givfld')
                             clstrg = clwork(1:icount)//cglatsuf
                             CALL hdlerr(NF_DEF_VAR(nc_outid(jf), 
     $               clstrg,n_reaty,2,ivardim,ilatgrd(jfint)),'givfld')
                         ENDIF
                     ENDIF
 130               CONTINUE
               ELSE
#endif
                   OPEN (UNIT = iunit,FILE = clfic,STATUS = 'UNKNOWN',
     $             FORM = 'UNFORMATTED',IOSTAT = iost)
                   IF (iost .ne. 0) THEN
                       CALL prtout
     $           ('Cannot create binary output file for field ',iloc,1)
                       CALL HALTE('STOP in givfld') 
                   ENDIF
#ifdef use_netCDF
               ENDIF
#endif
           ENDIF

C
C* Create the jf netcdf variables
C
#ifdef use_netCDF
           IF (lncdfrst) THEN
               ivardim(1)=idimlon(jf)
               ivardim(2)=idimlat(jf)
C                 IF ( nitfn .gt. 0) THEN
                     ivardim(3)=idimtime
                     CALL hdlerr(NF_DEF_VAR(nc_outid(jf), cnamout(iloc), 
     $                    n_reaty, 3, ivardim, ivarid(jf)), 'givfld')
C                 ELSE
C                    CALL hdlerr(NF_DEF_VAR(nc_outid(jf), cnamout(iloc), 
C    $                    n_reaty, 2, ivardim, ivarid(jf)), 'givfld')
C                 ENDIF
           ENDIF
#endif
 120     CONTINUE
C
C* Define header if info mode is on :
C  experiment name, initial date, iteration number, time since start
C
         IF (lmodinf) THEN 
             itime(1) = ndate
             itime(2) = kiter
             itime(3) = kiter * nstep
             IF (nlogprt .GE. 2) THEN
                 WRITE (UNIT = nulou,FMT = *) 
     $              ' Encapsulated data for current field is :'
                 CALL prcout ('Experiment name', cjobnam, 1)
                 CALL prtout ('Initial date', itime(1), 2)
                 CALL prtout ('Iteration number', itime(2), 2)
                 CALL prtout ('Time since start', itime(3), 2)
             ENDIF
         ENDIF
C
#ifdef use_netCDF
         IF (lncdfrst) THEN
C     
C            Loop on all fields for this iteration
             DO 140 jf = 1, kfield
               iloc=kindex(jf)
               isamefic=0
C
C              Test IF define mode has been already left. 
               DO 145 jj = 1, jf-1
                 IF (nc_outid(jf) .eq. nc_outid(jj)) THEN
                     isamefic=isamefic+1
                 ENDIF
 145           CONTINUE
C
               IF (isamefic .lt. 1) THEN
C
C*                 Put header as global attribute IF lmodinf
                   WRITE(nulou,*)'lmodinf=', lmodinf
                   IF (lmodinf) THEN
                      CALL hdlerr(NF_PUT_ATT_TEXT(nc_outid(jf),
     $                NF_GLOBAL,'Experiment_name',4,cjobnam),'givfld')
                      CALL hdlerr(NF_PUT_ATT_INT(nc_outid(jf),NF_GLOBAL
     $                ,'Initial_date',NF_INT,1,itime(1)),'givfld') 
                      CALL hdlerr(NF_PUT_ATT_INT(nc_outid(jf),NF_GLOBAL
     $                ,'Iteration_number',NF_INT,1,itime(2)),'givfld')
                      CALL hdlerr(NF_PUT_ATT_INT(nc_outid(jf),NF_GLOBAL
     $                ,'Time_since_start',NF_INT,1,itime(3)),'givfld')
                   ENDIF
C
C*                 Leave define mode for netcdf files
                   istatus=NF_ENDDEF(nc_outid(jf))
                   IF(istatus .ne. NF_NOERR) THEN
                       CALL prtout
     $   ('Cannot leave define mode for output file for field ',iloc,1)
                       CALL HALTE('STOP in givfld') 
                   ENDIF    
               ENDIF
C
C              Test IF coordinate variables have already been put. 
               isamegrd = 0
               iadrnew = nadrnew_grid(jf)
               isiznew = nsiznew(jf)
               DO 155 jj = 1, jf-1
                 ilocp=kindex(jj)
                 IF (nluout(iloc) .eq. nluout(ilocp) .and. cficaf(iloc)
     $               .eq. cficaf(ilocp)) isamegrd=isamegrd+1
 155           CONTINUE
C
C              If not, put coordinate variable to file jf
               IF (isamegrd .lt. 1) THEN
                   ist(1)=1 ; ist(2)=1
                   icnt(1)=nlonaf(iloc) ; icnt(2)=nlataf(iloc)
                   IF (n_reaty .eq. NF_FLOAT) THEN
                       istatus=NF_PUT_VARA_REAL 
     $                    (nc_outid(jf),ilongrd(jf),ist, icnt,
     $                     xgrnew(iadrnew:iadrnew+isiznew-1))

                   ELSE IF (n_reaty .eq. NF_DOUBLE) THEN
                       istatus=NF_PUT_VARA_DOUBLE 
     $                    (nc_outid(jf),ilongrd(jf),ist, icnt,
     $                     xgrnew(iadrnew:iadrnew+isiznew-1))
                   ENDIF
                   IF (n_reaty .eq. NF_FLOAT) THEN
                       istatus=NF_PUT_VARA_REAL 
     $                    (nc_outid(jf),ilatgrd(jf),ist, icnt,
     $                     ygrnew(iadrnew:iadrnew+isiznew-1))

                   ELSE IF (n_reaty .eq. NF_DOUBLE) THEN
                       istatus=NF_PUT_VARA_DOUBLE 
     $                    (nc_outid(jf),ilatgrd(jf),ist, icnt,
     $                     ygrnew(iadrnew:iadrnew+isiznew-1))
                   ENDIF
               ENDIF
 140         CONTINUE
         ENDIF
#endif
C
      ENDIF
C
      ENDIF
C
C*    2. Loop on active fields for iteration kiter
C        -----------------------------------------
C
!$omp parallel do default (shared)
!$omp+ private (jf,iloc,iadrnew,isiznew,clname,ilabel)
!$omp+ private (clabel,iunit,clstat)
!$omp+ reduction (+:infos)

      DO 210 jf = 1, kfield
C
C* Assign local variables
C
        iloc = kindex(jf)
        iadrnew = nadrnew(iloc)
        isiznew = nsiznew(iloc)
        clname = cnamout(iloc)
        ilabel = numlab(iloc)
        clabel = cfldlab(ilabel)
        iunit = nluout(iloc)
        clstat = cstate(iloc)
C
C* Test if field must be exported
C
        IF (clstat .EQ. 'EXPORTED') THEN 
C
C* Print field name
C
           IF (nlogprt .GE. 1) THEN
              CALL prcout('Writing of field : ', clname, 2)
              CALL prcout('Field definition : ', clabel, 2)
           ENDIF

C
C* - Give coupling fields
C
C* PIPE or NONE case
C
           IF (cchan .EQ. 'PIPE' .OR. cchan .EQ. 'NONE') THEN
C 
#ifdef use_netCDF    
              IF (lncdfrst) THEN
CEM
C  Get variables ID
CEM
                 istatus=NF_INQ_VARID(nc_outid(jf),clname,
     $                   ivarid(jf))
                 istatus=NF_INQ_VARID(nc_outid(jf),nc_invartime_name,
     $                   il_timevarid)
                 ist(1)=1 ; ist(2)=1 ; ist(3)=1+kiter
                 icnt(1)=nlonaf(iloc) ; icnt(2)=nlataf(iloc)
                 icnt(3)=1
                 IF (n_reaty .eq. NF_FLOAT) THEN
                    CALL hdlerr(NF_PUT_VARA_REAL
     $                   (nc_outid(jf),il_timevarid,ist(3),
     $                   icnt(3),r_time_val), 'givfld')
                    istatus=NF_PUT_VARA_REAL
     $                   (nc_outid(jf),ivarid(jf),ist, icnt,
     $                   fldnew(iadrnew:iadrnew+isiznew-1))
                    IF (istatus .ne. 0)  THEN
                       iflag(jf) = istatus
                       iflags= iflags + 1
                    ENDIF
                    
                 ELSE IF (n_reaty .eq. NF_DOUBLE) THEN
                    CALL hdlerr(NF_PUT_VARA_DOUBLE
     $                   (nc_outid(jf),il_timevarid,ist(3),
     $                    icnt(3),rtime_val), 'givfld')
                    istatus=NF_PUT_VARA_DOUBLE 
     $                   (nc_outid(jf),ivarid(jf),ist, icnt,
     $                   fldnew(iadrnew:iadrnew+isiznew-1))
                    IF (istatus .ne. 0)  THEN
                       iflag(jf) = istatus
                       iflags= iflags + 1
                    ENDIF
                 ENDIF
              ELSE
#endif
C* Write new field on unit iunit with header if required
C 
                 IF (lmodinf) THEN 
                    CALL locwrith (clname, cjobnam, itime, 
     $                   fldnew(iadrnew), isiznew, iunit, iflag(jf))
                 ELSE 
                    CALL locwrite (clname, fldnew(iadrnew), isiznew, 
     $                   iunit, iflag(jf))
                 ENDIF 
                 iflags = iflags + iflag(jf)
#ifdef use_netCDF
             ENDIF
#endif
          ELSE IF (cchan .EQ. 'SIPC' .OR. cchan .EQ. 'GMEM') THEN
C
C* SIPC case
C
#if defined use_comm_SIPC || defined use_comm_GMEM
C
              ipbytecha=kind('A')
              ipbyteint=kind(itime)
              ipbyterea=kind(fldnew)
C
C* Write encapsulated infos in field-specific shared memory pool
C
              IF (lmodinf) THEN 
                 isizeout = 3*ipbytecha
                 CALL SVIPC_write(mpoolidou(iloc), cjobnam 
     $                , isizeout, imrca)
                 isizeout = 3*ipbyteint
                 CALL SVIPC_write(mpoolidou(iloc), itime 
     $                ,isizeout,imrcb)
C
C* Find error if any
C
                 IF (imrca .LT. 0 .OR. imrcb .LT. 0) THEN
                    CALL prcout 
     $  ('Problem in writing encapsulated infos for field', clname, 1)
                    istop = 1
                 ENDIF 
                 IF (nlogprt .GE. 2) THEN
                    CALL prcout
     $    ('Wrote encapsulated infos in pool for field',clname,1)
                 ENDIF 
              ENDIF 
C
C* Write part of macro array in field-specific shared memory pool
C
              isizeout = isiznew * ipbyterea
              IF (nlogprt .GE. 2) THEN
                 WRITE(UNIT = nulou, FMT = *) 
     $                'Writing field data to pool = ',mpoolidou(iloc)
              ENDIF
              CALL SVIPC_write(mpoolidou(iloc),
     $             fldnew(iadrnew), isizeout, imrc)
C
C* Find error and stop if any
C
              IF (imrc .LT. 0) THEN
                 CALL prcout
     $               ('Problem in writing field in SHM pool:',clname,1) 
                 istop = 1
              ELSE IF (nlogprt .GE. 2) THEN
                 CALL prcout
     $                ('Wrote field in SHM pool:', clname, 1)
              ENDIF
C
C* CLIM case
C
#endif
          ELSE IF (cchan .EQ. 'MPI2' .or. cchan .EQ. 'MPI1' 
     $             .or. cchan .EQ. 'GSIP') THEN
#if defined use_comm_MPI1 || defined use_comm_MPI2  || defined use_comm_GSIP || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_GSIP && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE) 
C
C* Write new field on port clname
C
!$omp critical
              CALL CLIM_Export 
     $            (ig_portout_id(iloc), kiter*nstep, 
     $            fldnew(iadrnew), info(jf))
              infos = infos + info(jf)
!$omp end critical
#endif
          ELSE
              CALL prcout 
     $            ('Wrong CHANNEL option for field', clname, 1)
          ENDIF
      ELSE IF (clstat .EQ. 'AUXILARY') THEN 
          iflag(jf) = 0
          info(jf) = 0
      ENDIF
 210  CONTINUE
C
      IF (cchan .EQ. 'MPI2' .or. cchan .EQ. 'MPI1' 
     $    .or. cchan .EQ. 'GSIP' ) THEN
C*
#if defined use_comm_MPI1 || defined use_comm_MPI2  || defined use_comm_GSIP || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_GSIP && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE) 
C*    Stop if problem in writing in CLIM case
C
          IF (infos .NE. CLIM_Ok) THEN 
              DO 220 jf = 1, kfield
                IF (info(jf) .NE. CLIM_Ok) THEN
                    CALL prcout
     $                  ('WARNING: problem in writing field on port',
     $                  clname, 1)
                    CALL prtout
     $                  ('error code number', info(jf), 2)
                ENDIF
 220          CONTINUE
              CALL HALTE ('STOP in givfld') 
          ENDIF
C
#endif
      ELSE IF (cchan .EQ. 'PIPE' .or. cchan .eq. 'NONE') THEN
#if defined use_comm_PIPE || defined use_comm_NONE
C*    Stop if problem in writing in PIPE or NONE case
C
          IF (iflags .NE. 0) THEN 
              DO 230 jf = 1, kfield
                IF (iflag(jf) .NE. 0) THEN
                    iloc = kindex(jf)
                    iunit = nluout(iloc)
                    clname = cnamout(iloc)
                    CALL prcout
     $                  ('WARNING: problem in writing field',
     $                  clname, 1)
                    CALL prtout
     $                  ('Error writing on logical unit', iunit, 2)
                ENDIF
 230          CONTINUE
              CALL HALTE ('STOP in givfld')
          ENDIF
C
#endif
      ELSE IF (cchan .EQ. 'SIPC' .or. cchan .eq. 'GMEM') THEN
#if defined use_comm_SIPC || defined use_comm_GMEM
C*    Stop if problem in writing in SIPC or GMEM case
C
          IF (istop .NE. 0) CALL HALTE ('STOP in givfld')
#endif
C
      ENDIF
C
C*    3. PIPE Case: flush data files and send message
C        --------------------------------------------
C
      IF (cchan .EQ. 'PIPE' .or. cchan .eq. 'NONE') THEN
C 
          DO 310 jf = 1, kfield
C
C* Assign local variables
C
            iloc = kindex(jf)
            clname = cnamout(iloc)
            iunit = nluout(iloc)
            clfic = cficout(iloc)
            clstat = cstate(iloc)
C
C* Test if field must be exported
C
            IF (clstat .EQ. 'EXPORTED') THEN 
C
C* Close data file
C
                iloc=kindex(jf)
                isamefic=0
                DO 3005 jj = 1, jf-1
                  IF (nc_outid(jf) .eq. nc_outid(jj)) THEN
                      isamefic=isamefic+1
                  ENDIF
 3005           CONTINUE
C
                IF (isamefic .lt. 1) THEN
#ifdef use_netCDF
                    IF (lncdfrst) THEN
                      IF (kiter.eq.nitfn) then
                        CALL hdlerr(NF_CLOSE(nc_outid(jf)), 'givfld')
                     ENDIF
                    ELSE
#endif
                        CLOSE(UNIT = iunit, ERR = 3010, IOSTAT = ios)
                        IF (nlogprt .GE. 2) THEN
                            WRITE(UNIT = nulou,FMT = 3100) iunit,clfic
                        ENDIF
 3010                   CONTINUE
                        IF (ios .NE. 0) THEN
                            CALL prtout('Problem in closing unit',
     $                          iunit, 2)
                            CALL prtout('Error message nbr is= ',ios,2)
                            CALL HALTE('STOP in givfld')
                        ENDIF 
#ifdef use_netCDF
                    ENDIF
#endif
                ENDIF
C
C* Send message on pipe clname for PIPE technique
C
#ifdef use_comm_PIPE
                CALL PIPE_Send (clname, kiter)
#endif
            ELSE IF (clstat .EQ. 'AUXILARY' .AND. nlogprt .GE. 2) THEN 
                WRITE(UNIT = nulou, FMT = 3300) clname, clfic
            ENDIF 
 310      CONTINUE
      ENDIF
 3100 FORMAT(/,5X,' Unit ',I2,' has been disconnected from file ',A8)
 3200 FORMAT(/,5X,' Unit ',I2,' has been reconnected to file ',A8)
 3300 FORMAT(/,5X,' Auxilary field ',A8,' is not written on file ',A8)
C
C
C*    4. End of routine
C        --------------
C
      IF (nlogprt .GE. 2) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '          --------- End of routine givfld ---------'
          CALL FLUSH (nulou)
      ENDIF
      DEALLOCATE (info)
      DEALLOCATE (iflag)
      DEALLOCATE (idimlon)
      DEALLOCATE (idimlat)
      DEALLOCATE (ivarid)
      DEALLOCATE (ilongrd)
      DEALLOCATE (ilatgrd)
      RETURN
      END

