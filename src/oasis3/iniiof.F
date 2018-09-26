      SUBROUTINE iniiof
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 0 *
C               * -------------     ------- *
C               *****************************
C
C**** *iniiof*  - Open files
C
C     Purpose:
C     -------
C     Open files for grid related parameters, data transfer and 
C     auxilary outputs
C
C**   Interface:
C     ---------
C       *CALL*  *iniiof*
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
C     None
C
C     Reference:
C     ---------
C     See OASIS manual (1995) 
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       1.0       L. Terray      94/01/01  created
C       2.0beta   L. Terray      95/08/23  modified: new structure
C       2.0       L. Terray      96/02/01  modified: no opening of
C                                          CLIM trace file within oasis
C       2.1       L. Terray      96/08/26  modified: open data file for
C                                          SUBGRID Analysis
C       2.3       S. Valcke      99/03/30  modified: open data file for
C                                          NINENN analysis
C       2.3       S. Valcke      99/04/30  added: printing levels
C       2.5       S. Valcke      01/03/23  removed: nulrd opening
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------- Include files and USE of modules---------------------------
C
      USE mod_parameter
      USE mod_unitncdf
      USE mod_string
      USE mod_anais
      USE mod_extrapol
      USE mod_analysis
      USE mod_unit
      USE mod_hardware
      USE mod_label
      USE mod_printing
#ifdef use_netCDF
#include <netcdf.inc>
#endif
C
C* ---------------------------- Local declarations ----------------------
C
      CHARACTER*8 clfic
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C
C* This routine will be called only if one field (at least) goes through Oasis
C
      IF (lg_oasis_field) THEN
C     
C*    1. Open files
C        ----------
C     
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) 
     $           '           ROUTINE iniiof  -  Level 0'
            WRITE (UNIT = nulou,FMT = *) 
     $           '           **************     *******'
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) ' opening grid related files'
            WRITE (UNIT = nulou,FMT = *) ' '
            iost = 0
            WRITE (UNIT = nulou,FMT = *) ' opening gcms grid file '
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
         IF (.not. lncdfgrd) THEN
            OPEN (UNIT = nulgr,FILE = cgrdnam,STATUS='OLD',
     $           FORM ='UNFORMATTED',ERR = 110,IOSTAT = iost)
            IF (nlogprt .GE. 1) THEN
               WRITE (UNIT = nulou,FMT = 1001) nulgr, cgrdnam
               WRITE (UNIT = nulou,FMT = *) ' '
            ENDIF
 110        CONTINUE
c           Try whether a new netcdf grids file was created at runtime
            IF (iost .NE. 0) THEN
               lncdfgrd = .true.
            ENDIF
         ENDIF
#ifdef use_netCDF
         IF (lncdfgrd) THEN 
            istatus=NF_OPEN(cgrdnam//'.nc', NF_NOWRITE, nc_grdid)
            IF (istatus .ne. NF_NOERR) THEN
               WRITE (UNIT = nulou,FMT = *) 
     $              ' ===>>>> NO grids neither grids.nc files.'
               CALL halte('STOP in iniiof')
            ENDIF
         ENDIF
#endif
C     
C* Masks file
C     
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = *) '   opening gcm masks file '
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
#ifdef use_netCDF
         IF (lncdfgrd) THEN
            istatus=NF_OPEN(cmsknam//'.nc', NF_NOWRITE, nc_mskid)
            IF (istatus .ne. NF_NOERR) THEN
               WRITE (UNIT = nulou,FMT = *) 
     $              ' ===>>>> : error in opening masks.nc. We STOP!'
               WRITE (UNIT = nulou,FMT = *) ' '
               CALL HALTE('STOP in iniiof')   
            ENDIF
         ELSE
#endif
            OPEN (UNIT = nulma,FILE = cmsknam,STATUS='OLD',
     $           FORM ='UNFORMATTED',ERR = 120,IOSTAT = iost)
            IF (nlogprt .GE. 1) THEN
               WRITE (UNIT = nulou,FMT = 1001) nulma, cmsknam
               WRITE (UNIT = nulou,FMT = *) ' '
            ENDIF
 120        CONTINUE
            IF (iost .ne. 0) THEN
               WRITE (UNIT = nulou,FMT = *) 
     $              ' ===>>>> : error opening masks file'
               WRITE (UNIT = nulou,FMT = *) 
     $              ' =======   =====               ===='
               WRITE (UNIT = nulou,FMT = *) 
     $              ' logical unit ',nulma,' error number = ',
     $              iost
               WRITE (UNIT = nulou,FMT = *) 
     $              ' We STOP. Verify the file ', cmsknam
               WRITE (UNIT = nulou,FMT = *) ' '
               CALL HALTE('STOP in iniiof')   
            ENDIF
#ifdef use_netCDF
         ENDIF
#endif
C     
C* Surfaces file
C     
         DO 125 jf = 1, ig_nfield
            IF (lsurf(jf)) lg_areas = .TRUE.
 125     CONTINUE
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = *) 'Trying open gcm surface file'
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
#ifdef use_netCDF
         IF (lncdfgrd) THEN
            istatus=NF_OPEN(csurnam//'.nc', NF_NOWRITE, nc_surid)
            IF (istatus .ne. NF_NOERR .and. lg_areas) THEN
               WRITE (UNIT = nulou,FMT = *) 
     $              ' ===>>>> : error in opening areas.nc. We STOP!'
               WRITE (UNIT = nulou,FMT = *) ' '
               CALL HALTE('STOP in iniiof')   
            ENDIF
         ELSE
#endif
            OPEN (UNIT = nulsu,FILE = csurnam,
     $           FORM ='UNFORMATTED',ERR = 130,IOSTAT = iost)
            IF (nlogprt .GE. 1) THEN
               WRITE (UNIT = nulou,FMT = 1001) nulsu, csurnam
               WRITE (UNIT = nulou,FMT = *) ' '
            ENDIF
 130        CONTINUE
            IF (iost .ne. 0 .and. lg_areas) THEN
               WRITE (UNIT = nulou,FMT = *) 
     $              ' ===>>>> : error opening surfaces file'
               WRITE (UNIT = nulou,FMT = *) 
     $              ' =======   =====                  ===='
               WRITE (UNIT = nulou,FMT = *) 
     $              ' logical unit ',nulsu,' error number = ', iost
               WRITE (UNIT = nulou,FMT = *) 
     $              ' We STOP. Verify the file ', csurnam
               WRITE (UNIT = nulou,FMT = *) ' '
               CALL HALTE ('STOP in iniiof')   
            ENDIF
#ifdef use_netCDF
         ENDIF
#endif
C     
C* Trace file for CLIM 
C     
         IF (nlogprt .GE. 1) THEN
            IF (cchan .EQ. 'MPI2' .or. cchan .EQ. 'MPI1' 
     $           .or. cchan .EQ. 'GSIP') THEN
               WRITE (UNIT = nulou,FMT = *) 
     $              '      The CLIM trace file is opened'
               WRITE (UNIT = nulou,FMT = *) 
     $              '      within the CLIM program '
               WRITE (UNIT = nulou,FMT = *) 
     $              ' The name of the file is oasis.prt '
               WRITE (UNIT = nulou,FMT = *) 
     $              ' linked to logical unit = ', nultr
            ENDIF
         ENDIF
C     
C* Formats
C     
 1001    FORMAT(10X,' open unit = ',I3,'    file ',A8,' ok')
C     
C     
C*    2. Deal with restart/transfer files
C     --------------------------------
C     
C* Files for input fields
C     
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) 
     $           ' opening restart/transfer files '
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
#ifdef use_netCDF
         lncdfrst = .true.
#endif
         DO jf = 1, ig_nfield
            isamefic=0
            DO 215 jj = 1, jf-1
               IF (nluinp(jf) .eq. nluinp(jj)) THEN
                  isamefic=isamefic+1
                  nc_inpid(jf)=nc_inpid(jj)
               ENDIF
 215        END DO
            IF (isamefic .lt. 1) THEN
               iunit = nluinp(jf)
               clfic = cficinp(jf)
c  Find the index of jf in total number of field
               DO jj = 1, ig_total_nfield
	          if ( jf .eq. ig_number_field(jj)) ig_ind = jj
               END DO
#if !defined use_comm_NONE
               IF (ig_lag(ig_ind) .gt. 0) THEN
#endif
#ifdef use_netCDF
                  IF (lncdfrst) THEN
                     istatus=NF_OPEN(clfic, NF_NOWRITE, nc_inpid(jf))
                     IF (istatus .ne. NF_NOERR) THEN
#endif
                        OPEN (UNIT = iunit,FILE = clfic,STATUS='OLD',
     $                       FORM = 'UNFORMATTED',IOSTAT = iost)
                        IF (iost .ne. 0) THEN
                           CALL prtout
     $                          ('No restart input file for field ',
     $                          jf, 1)
                           CALL HALTE('STOP in iniiof')
                        ELSE
                           IF (nlogprt .GE. 1) THEN
                              WRITE (UNIT = nulou,FMT = 1001)iunit,clfic 
                           ENDIF
                           lncdfrst = .false.
                        ENDIF
#ifdef use_netCDF
                     ELSE
                        IF (nlogprt .GE. 1) THEN
                           WRITE (UNIT = nulou,FMT = *)
     $                          '   opened netcdf restart file ', clfic
                        ENDIF
                     ENDIF
                  ELSE
                     OPEN (UNIT = iunit,FILE = clfic,STATUS='OLD',
     $                    FORM = 'UNFORMATTED',IOSTAT = iost)
                     IF (iost .ne. 0) THEN
                        CALL prtout
     $                       ('No binary restart input file for field ',
     $                       jf, 1)
                        CALL HALTE('STOP in iniiof')
                     ELSE
                        IF (nlogprt .GE. 1) THEN
                           WRITE (UNIT = nulou,FMT = 1001) iunit, clfic
                        ENDIF
                     ENDIF
                  ENDIF
#endif
#if !defined use_comm_NONE
              ENDIF
#endif
            ENDIF
         ENDDO
C     
C     *    3. Anais weights and output files and NINENN weight file
C     -----------------------------------------------------
C     
C     * Weights file for ANAISM
C     
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) 
     $           ' open interpolation related files '
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) 
     $           '      need file with surface mesh weights'
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
         OPEN (UNIT = nulcc,FILE = cwanaism,
     $        FORM ='UNFORMATTED',ERR = 310,IOSTAT = iost)
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = 1001) nulcc, cwanaism
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
 310     CONTINUE
         IF (iost .ne. 0) THEN
            WRITE (UNIT = nulou,FMT = *) 
     $           ' ===>>>> : error opening surfmesh weights file'
            WRITE (UNIT = nulou,FMT = *) 
     $           ' =======   =====                          ===='
            WRITE (UNIT = nulou,FMT = *) 
     $           ' logical unit ',nulcc,' error number = ',
     $           iost
            WRITE (UNIT = nulou,FMT = *) 
     $           ' We STOP. Verify the file', cwanaism
            WRITE (UNIT = nulou,FMT = *) ' '
            CALL HALTE ('STOP in iniiof')   
         ENDIF
C     * Weights file for ANAISG
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = *) 
     $           '      need file with curvilinear grid weights'
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
         OPEN (UNIT = nulgg,FILE = cwanaisg,
     $        FORM ='UNFORMATTED',ERR = 320,IOSTAT = iost)
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = 1001) nulgg, cwanaisg
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
 320     CONTINUE 
         IF (iost .ne. 0) THEN
            WRITE (UNIT = nulou,FMT = *) 
     $           ' ===>>>> : error opening curvilinear weights file'
            WRITE (UNIT = nulou,FMT = *) 
     $           ' =======   =====                             ===='
            WRITE (UNIT = nulou,FMT = *) 
     $           ' logical unit ',nulgg,' error number = ',
     $           iost
            WRITE (UNIT = nulou,FMT = *) 
     $           ' We STOP. Verify the file', cwanaisg
            WRITE (UNIT = nulou,FMT = *) ' '
            CALL HALTE ('STOP in iniiof')   
         ENDIF
C     * ANAIS output FILE
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = *) '      open ANAIS output file'
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
         OPEN (UNIT = nulan,FILE = cnaisout, STATUS='NEW',
     $        FORM ='FORMATTED',ERR = 330,IOSTAT = iost)
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = 1001) nulan, cnaisout
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
 330     CONTINUE
         IF (iost .ne. 0) THEN
            WRITE (UNIT = nulou,FMT = *) 
     $           ' ===>>>> : error opening ANAIS output file'
            WRITE (UNIT = nulou,FMT = *) 
     $           ' =======   =====                      ===='
            WRITE (UNIT = nulou,FMT = *) 
     $           ' logical unit ',nulan,' error number = ',
     $           iost
            WRITE (UNIT = nulou,FMT = *) 
     $           ' We STOP. Verify the file', cnaisout
            WRITE (UNIT = nulou,FMT = *) ' '
            CALL HALTE ('STOP in iniiof')   
         ENDIF
C     
C     * Weights file for NINENN
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = *) 
     $           'need file with weights, address and iteration number'
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
#if !defined key_openmp
         OPEN (UNIT = nulgn,FILE = cwninenn,
     $        FORM ='UNFORMATTED',ERR = 340,IOSTAT = iost)
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = 1001) nulgn, cwninenn
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
 340     CONTINUE 
         IF (iost .ne. 0) THEN
            WRITE (UNIT = nulou,FMT = *) 
     $           ' ===>>>> : error opening weight file'
            WRITE (UNIT = nulou,FMT = *) 
     $           ' =======   =====                ===='
            WRITE (UNIT = nulou,FMT = *) 
     $           ' logical unit ',nulgn,' error number = ',
     $           iost
            WRITE (UNIT = nulou,FMT = *) 
     $           ' We STOP. Verify the file', cwninenn
            WRITE (UNIT = nulou,FMT = *) ' '
            CALL HALTE ('STOP in iniiof')   
         ENDIF
#endif
C
C*    4. Deal with data files used in analysis routines
C        ----------------------------------------------
C
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) ' open analysis related files '
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
         DO 410 jf = 1, ig_nfield
            DO 420 ja = 1, ig_ntrans(jf)
               IF(canal(ja,jf) .EQ. 'CORRECT') THEN
                  DO 430 jc = 1, ncofld(jf)
                     iunit = nludat(jc,jf)
                     clfic = ccofic(jc,jf)
                     OPEN(UNIT = iunit,FILE = clfic,
     $                    FORM ='UNFORMATTED',ERR = 440,IOSTAT = iost)
                     IF (nlogprt .GE. 1) THEN
                        WRITE (UNIT = nulou,FMT = 1001) iunit, clfic
                        WRITE (UNIT = nulou,FMT = *) ' '
                     ENDIF
 440                 IF (iost .ne. 0) THEN
                        WRITE (UNIT = nulou,FMT = *) 
     $                       ' ===>>>> : error opening data file'
                        WRITE (UNIT = nulou,FMT = *) 
     $                       ' =======   =====              ===='
                        WRITE (UNIT = nulou,FMT = *) 
     $                       ' logical unit ',iunit,' error number = ',
     $                       iost
                        WRITE (UNIT = nulou,FMT = *) 
     $                       ' We STOP. Verify the file ', clfic
                        WRITE (UNIT = nulou,FMT = *) ' '
                        CALL HALTE('STOP in iniiof')
                     ENDIF
 430              CONTINUE 
               ELSE IF(canal(ja,jf) .EQ. 'MOZAIC') THEN
                  iunit = nlumap(jf)
                  clfic = cgrdmap(jf)
                  OPEN(UNIT = iunit,FILE = clfic, 
     $                 FORM ='UNFORMATTED',ERR = 450,IOSTAT = iost)
                  IF (nlogprt .GE. 1) THEN
                     WRITE (UNIT = nulou,FMT = 1001) iunit, clfic
                     WRITE (UNIT = nulou,FMT = *) ' '
                  ENDIF
 450              IF (iost .ne. 0) THEN
                     WRITE (UNIT = nulou,FMT = *) 
     $                    ' ===>>>> : error opening mapping file'
                     WRITE (UNIT = nulou,FMT = *) 
     $                    ' =======   =====                 ===='
                     WRITE (UNIT = nulou,FMT = *) 
     $                    ' logical unit ',iunit,' error number = ',
     $                    iost
                     WRITE (UNIT = nulou,FMT = *) 
     $                    ' We STOP. Verify the file ', clfic
                     WRITE (UNIT = nulou,FMT = *) ' '
                     CALL HALTE('STOP in iniiof')
                  ENDIF
               ELSE IF(canal(ja,jf) .EQ. 'FILLING') THEN
                  iunit = nlufil(jf)
                  clfic = cfilfic(jf)
                  OPEN(UNIT = iunit,FILE = clfic,
     $                 FORM ='UNFORMATTED',ERR = 460,IOSTAT = iost)
                  IF (nlogprt .GE. 1) THEN
                     WRITE (UNIT = nulou,FMT = 1001) iunit, clfic
                     WRITE (UNIT = nulou,FMT = *) ' '
                  ENDIF
 460              IF (iost .ne. 0) THEN
                     WRITE (UNIT = nulou,FMT = *) 
     $                    ' ===>>>> : error opening filling data file'
                     WRITE (UNIT = nulou,FMT = *) 
     $                    ' =======   =====         ======= ===='
                     WRITE (UNIT = nulou,FMT = *) 
     $                    ' logical unit ',iunit,' error number = ',
     $                    iost
                     WRITE (UNIT = nulou,FMT = *) 
     $                    ' We STOP. Verify the file ', clfic
                     WRITE (UNIT = nulou,FMT = *) ' '
                     CALL HALTE('STOP in iniiof')
                  ENDIF
               ELSE IF(canal(ja,jf) .EQ. 'SUBGRID') THEN
                  iunit = nlusub(jf)
                  clfic = cgrdsub(jf)
                  OPEN(UNIT = iunit,FILE = clfic,
     $                 FORM ='UNFORMATTED',ERR = 470,IOSTAT = iost)
                  IF (nlogprt .GE. 1) THEN
                     WRITE (UNIT = nulou,FMT = 1001) iunit, clfic
                     WRITE (UNIT = nulou,FMT = *) ' '
                  ENDIF
 470              IF (iost .ne. 0) THEN
                     WRITE (UNIT = nulou,FMT = *) 
     $                    ' ===>>>> : error opening subgrid data file'
                     WRITE (UNIT = nulou,FMT = *) 
     $                    ' =======   =====         ======= ===='
                     WRITE (UNIT = nulou,FMT = *) 
     $                    ' logical unit ',iunit,' error number = ',
     $                    iost
                     WRITE (UNIT = nulou,FMT = *) 
     $                    ' We STOP. Verify the file ', clfic
                     WRITE (UNIT = nulou,FMT = *) ' '
                     CALL HALTE('STOP in iniiof')
                  ENDIF
               ENDIF
 420        CONTINUE 
 410     CONTINUE 
C
C
C*    5. End of routine
C        --------------
C
         IF (nlogprt .GE. 1) THEN
            WRITE(UNIT = nulou,FMT = *) ' '
            WRITE(UNIT = nulou,FMT = *) 
     $           '          ---------- End of routine iniiof ---------'
            CALL FLUSH (nulou)
         ENDIF
      ENDIF
      RETURN
      END
