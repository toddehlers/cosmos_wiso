      SUBROUTINE inigrd
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 0 *
C               * -------------     ------- *
C               *****************************
C
C**** *inigrd*  - Initialize gcm's grids
C
C     Purpose:
C     -------
C     Get atmospheric and oceanic grids and masks as well as scale factors
C     for interpolation purposes.
C
C     ***NOTE***
C     ----------
C     Grids, masks and scale factors must be ordered with the OASIS convention.
C     They go from South to North and from Greenwhich to the east.
C     EXCEPTION: array maskr (mask for reduced gaussian grid) must
C     be ordered as in the atmospheric code (North to South in 
C     Arpege case).
C
C**   Interface:
C     ---------
C       *CALL*  *inigrd*
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
C     locread, locrint
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
C       2.0       L. Terray      95/08/31  modified : new structure
C       2.3       S. Valcke      99/04/30  added: printing levels
C       2.5       S. Valcke      01/03/23  added: netCDF auxilary files
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------------------- Include files ---------------------------
C
      USE mod_kinds_oasis
      USE mod_parameter
      USE mod_memory
      USE mod_string
      USE mod_unitncdf
      USE mod_analysis
      USE mod_label
      USE mod_unit
      USE mod_printing
#ifdef use_netCDF
#include <netcdf.inc>
#endif
C
C* ---------------------------- Local declarations ----------------------
C
      CHARACTER*8 clwork
      CHARACTER*8 clstrg,clstrg1,clstrg2,clstrg3,clstrg4
      CHARACTER*8 clstrg5,clstrg6
      INTEGER (kind=ip_intwp_p) il_varid
      INTEGER (kind=ip_intwp_p) ist(2), icnt(2)
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C
C* This routine will be called only if one field (at least) goes through Oasis
C
      IF (ig_nfield .ne.0) THEN
C
C*    1. Initialization
C        --------------
C
      IF (nlogprt .GE. 1) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '           ROUTINE inigrd  -  Level 0'
          WRITE (UNIT = nulou,FMT = *) 
     $    '           **************     *******'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    ' get gcms grids, masks and surfaces'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
      ENDIF
      ist(:)=0
      icnt(:)=0
C
C* Initializa error codes
C
      infos = 0
      iflag = 0
C
C
C*    2. Get gcm's grids, masks and mesh surfaces
C        ----------------------------------------
C
!$omp parallel do default (shared)
!$omp+ private (jf,iadrold,isizold,iadrnew,isiznew)
!$omp+ private (clwork,icount,clstrg,clstrg1)
!$omp+ private (clstrg2,clstrg3,clstrg4,clstrg5)
!$omp+ private (clstrg6)
!$omp+ private (istatus,iunit1,inunit2,iunit3)
!$omp+ private (iunit4,iunit5,iunit6,iunit7)
!$omp+ private (iflag,ist,icnt,nc_file,n_intty)
!$omp+ private (n_reaty,il_varid)
!$omp+ reduction (+:infos)

      DO 200 jf = 1, ig_nfield
C
C* Print
C
        IF (nlogprt .GE. 2) THEN
            CALL prtout
     $          ('Start reading source grid data for field nbr=',jf,2)
        ENDIF
C
C* Assign local variables to array sizes and pointers
C
        iadrold = nadrold_grid(jf)
        isizold = nsizold(jf)
        iadrnew = nadrnew_grid(jf)
        isiznew = nsiznew(jf)
C
C* Get the right locator string for initial grids files
C
        clwork = cficbf(jf)
        icount = ilenstr(clwork,jpeight)
        clstrg1 = clwork(1:icount)//cglonsuf
C
C* Read longitudes of initial grid associated to field jf
C
        iunit1 = nulgr
#ifdef use_netCDF
        IF (lncdfgrd) THEN
C Next BLOCK done only once per field
            ist(1)=1 ; ist(2)=1
            icnt(1)=nlonbf(jf) ; icnt(2)=nlatbf(jf)
C
            nc_file = nc_grdid
            CALL hdlerr
     $          (NF_INQ_VARID(nc_file, clstrg1, il_varid), 'inigrd')
            CALL hdlerr
     $          (NF_INQ_VARTYPE(nc_file, il_varid, n_reaty), 'inigrd')
            IF (n_reaty .eq. NF_FLOAT) THEN
                 IF ( ip_realwp_p .ne. ip_single_p ) THEN
                    WRITE (nulou, *)
     $             'Incoherence: Oasis is compiled in double precision'
                    WRITE (nulou, *)
     $   'but Grid auxiliary file contains single precision REAL'
                    CALL HALTE('STOP in getfld')
                 ENDIF
               CALL hdlerr(NF_GET_VARA_REAL
     $              (nc_file, il_varid, ist, icnt, 
     $               xgrold(iadrold:iadrold+isizold-1)),'inigrd')
            ELSE IF (n_reaty .eq. NF_DOUBLE) THEN
                 IF ( ip_realwp_p .ne. ip_double_p ) THEN
                    WRITE (nulou, *)
     $             'Incoherence: Oasis is compiled in single precision'
                    WRITE (nulou, *)
     $   'but Grid auxiliary file contains double precision REAL'
                    CALL HALTE('STOP in getfld')
                 ENDIF               
               CALL hdlerr(NF_GET_VARA_DOUBLE 
     $              (nc_file, il_varid, ist, icnt, 
     $               xgrold(iadrold:iadrold+isizold-1)),'inigrd')
            ELSE
                CALL prcout ('Problem with type of array =',clstrg1,1) 
                CALL HALTE('STOP in inigrd')
            ENDIF
        ELSE
#endif
            CALL locread (clstrg1,xgrold(iadrold),isizold,iunit1,iflag)
            IF (iflag .NE. 0) THEN
                CALL prtout('Problem in reading unit =',iunit1,2)
                CALL prcout
     $          ('Problem with array linked to string =',clstrg1,1)
                infos = infos + iflag
            ENDIF
#ifdef use_netCDF
        ENDIF
#endif
        IF (nlogprt .GE. 2) THEN
            CALL prtout('Read source grid longitudes for field #',jf,1)
        ENDIF
C
C* Read latitudes of initial grid associated to field jf
C
        clstrg2 = clwork(1:icount)//cglatsuf
#ifdef use_netCDF
        IF (lncdfgrd) THEN
            CALL hdlerr
     $          (NF_INQ_VARID(nc_file, clstrg2, il_varid), 'inigrd')
            IF (n_reaty .eq. NF_FLOAT) THEN
                CALL hdlerr(NF_GET_VARA_REAL 
     $              (nc_file, il_varid, ist, icnt, 
     $              ygrold(iadrold:iadrold+isizold-1)), 'inigrd')
            ELSE IF (n_reaty .eq. NF_DOUBLE) THEN
                CALL hdlerr(NF_GET_VARA_DOUBLE 
     $              (nc_file, il_varid, ist, icnt, 
     $              ygrold(iadrold:iadrold+isizold-1)), 'inigrd')
            ENDIF
        ELSE
#endif
            CALL locread(clstrg2,ygrold(iadrold),isizold,iunit1,iflag)
            IF (iflag .NE. 0) THEN
                CALL prtout('Problem in reading unit =',iunit1,2)
                CALL prcout
     $          ('Problem with array linked to string =',clstrg2,1)
                infos = infos + iflag
            ENDIF
#ifdef use_netCDF
        ENDIF
#endif
        IF (nlogprt .GE. 2) THEN
            CALL prtout('Read source grid latitudes for field # ',jf,1)
        ENDIF
C
C* Read initial mask associated to field jf
C
        clstrg3 = clwork(1:icount)//cmsksuf
        iunit2 = nulma
#ifdef use_netCDF
        IF (lncdfgrd) THEN
            CALL hdlerr
     $          (NF_INQ_VARID(nc_mskid, clstrg3, il_varid), 'inigrd')
            CALL hdlerr
     $          (NF_INQ_VARTYPE(nc_mskid, il_varid, n_intty), 'inigrd')
            IF (n_intty .eq. NF_INT) THEN
                CALL hdlerr(NF_GET_VARA_INT 
     $              (nc_mskid, il_varid, ist, icnt, 
     $              mskold(iadrold:iadrold+isizold-1)), 'inigrd')
            ELSE
                CALL prcout ('Problem with type of array =',clstrg3,1) 
                CALL HALTE('STOP in inigrd')
            ENDIF
        ELSE
#endif
            CALL locrint (clstrg3,mskold(iadrold),isizold,iunit2,iflag)
            IF (iflag .NE. 0) THEN
                CALL prtout('Problem in reading unit =',iunit2,2)
                CALL prcout
     $          ('Problem with array linked to string =',clstrg3,1)
                infos = infos + iflag
            ENDIF
#ifdef use_netCDF
        ENDIF
#endif
        IF (nlogprt .GE. 2) THEN
            CALL prtout
     $          ('Read source grid mask for field number =',jf,1)
        ENDIF
C
C* Read initial mesh surface associated to field jf
C
        clstrg4 = clwork(1:icount)//csursuf
        iunit3 = nulsu
#ifdef use_netCDF
        IF (lncdfgrd) THEN
            istatus = 
     $          NF_INQ_VARID(nc_surid, clstrg4, il_varid)
            IF (istatus .ne. NF_NOERR .and. lsurf(jf)) THEN
                CALL prtout('Cannot read surface information for fld '
     $              , jf, 1)
                CALL HALTE ('STOP in inigrd')
            ELSE IF (istatus .eq. NF_NOERR) THEN
                IF (n_reaty .eq. NF_FLOAT) THEN
                    istatus = NF_GET_VARA_REAL(nc_surid, 
     $                   il_varid, ist, icnt, 
     $                  surold(iadrold:iadrold+isizold-1))
                    IF (istatus .ne. NF_NOERR .and. lsurf(jf)) THEN
                        CALL prtout('Cannot read source grid surface 
     $ information for field ' , jf, 1)
                        CALL HALTE ('STOP in inigrd')
                    ENDIF
                ELSE IF (n_reaty .eq. NF_DOUBLE) THEN
                    istatus = NF_GET_VARA_DOUBLE(nc_surid, 
     $                   il_varid, ist, icnt, 
     $                  surold(iadrold:iadrold+isizold-1))
                    IF (istatus .ne. NF_NOERR .and. lsurf(jf)) THEN
                        CALL prtout('Cannot read source grid surface 
     $ information for field ' , jf, 1)
                        CALL HALTE ('STOP in inigrd')
                    ENDIF
                ENDIF
            ENDIF
        ELSE
#endif
            CALL locread (clstrg4, surold(iadrold), isizold,
     $          iunit3,iflag)
            IF (iflag .NE. 0 .and. lsurf(jf)) THEN
                CALL prtout('Problem in reading unit =',iunit3,2)
                CALL prcout('Problem with array linked to string ',
     $              clstrg4,1)
                infos = infos + iflag
            ENDIF
#ifdef use_netCDF
        ENDIF
#endif
        IF (nlogprt .GE. 2) THEN
            CALL prtout
     $          ('Read source grid surface for field number ',jf,1)
        ENDIF
C
C* Get the right locator string for final grids files
C
!$omp critical
        clwork = cficaf(jf)
        icount = ilenstr(clwork,jpeight)
        clstrg5 = clwork(1:icount)//cglonsuf
        IF (nlogprt .GE. 2) THEN
            CALL prtout
     $          ('Start reading target grid data for field nbr =',jf,2)
        ENDIF
C
C* Read longitudes of final grid associated to field jf
C
        iunit4 = nulgr
#ifdef use_netCDF
        IF (lncdfgrd) THEN
C Next BLOCK done only once per field
            ist(1)=1 ; ist(2)=1
            icnt(1)=nlonaf(jf) ; icnt(2)=nlataf(jf)
C
            nc_file = nc_grdid
            CALL hdlerr
     $          (NF_INQ_VARID(nc_file, clstrg5, il_varid), 'inigrd')
            IF (n_reaty .eq. NF_FLOAT) THEN
                CALL hdlerr(NF_GET_VARA_REAL
     $              (nc_file, il_varid, ist, icnt, 
     $              xgrnew(iadrnew:iadrnew+isiznew-1)), 'inigrd')
            ELSE IF (n_reaty .eq. NF_DOUBLE) THEN
                CALL hdlerr(NF_GET_VARA_DOUBLE 
     $              (nc_file, il_varid, ist, icnt, 
     $              xgrnew(iadrnew:iadrnew+isiznew-1)), 'inigrd')
            ENDIF
        ELSE
#endif
            CALL locread (clstrg5,xgrnew(iadrnew),isiznew,iunit4,iflag)
            IF (iflag .NE. 0) THEN
                CALL prtout('Problem in reading unit =',iunit4,2)
                CALL prcout
     $          ('Problem with array linked to string =',clstrg5,1)
                infos = infos + iflag
            ENDIF
#ifdef use_netCDF
        ENDIF
#endif
!$omp end critical
        IF (nlogprt .GE. 2) THEN
            CALL prtout
     $          ('Read target grid longitudes for field number =',jf,1)
        ENDIF
C
C* Read latitudes of final grid associated to field jf
C
!$omp critical
        clstrg6 = clwork(1:icount)//cglatsuf
        iunit5 = nulgr
#ifdef use_netCDF
        IF (lncdfgrd) THEN
            CALL hdlerr
     $          (NF_INQ_VARID(nc_file, clstrg6, il_varid), 'inigrd')
            IF (n_reaty .eq. NF_FLOAT) THEN
                CALL hdlerr(NF_GET_VARA_REAL 
     $              (nc_file, il_varid, ist, icnt, 
     $              ygrnew(iadrnew:iadrnew+isiznew-1)), 'inigrd')
            ELSE IF (n_reaty .eq. NF_DOUBLE) THEN
                CALL hdlerr(NF_GET_VARA_DOUBLE 
     $              (nc_file, il_varid, ist, icnt, 
     $              ygrnew(iadrnew:iadrnew+isiznew-1)), 'inigrd')
            ENDIF
        ELSE
#endif
            CALL locread (clstrg6,ygrnew(iadrnew),isiznew,iunit5,iflag)
            IF (iflag .NE. 0) THEN
                CALL prtout('Problem in reading unit =',nulgr,2)
                CALL prcout
     $          ('Problem with array linked to string =',clstrg6,1)
                infos = infos + iflag
            ENDIF
#ifdef use_netCDF
        ENDIF
#endif
!$omp end critical
        IF (nlogprt .GE. 2) THEN
            CALL prtout
     $          ('Read target grid latitudes for field number =',jf,1)
        ENDIF
C
C* Read final mask associated to field jf
C
!$omp critical
        clstrg = clwork(1:icount)//cmsksuf
        iunit6 = nulma
#ifdef use_netCDF
        IF (lncdfgrd) THEN
            CALL hdlerr
     $          (NF_INQ_VARID(nc_mskid, clstrg, il_varid), 'inigrd')
            IF (n_intty .eq. NF_INT) THEN
                CALL hdlerr(NF_GET_VARA_INT 
     $              (nc_mskid, il_varid, ist, icnt, 
     $              msknew(iadrnew:iadrnew+isiznew-1)), 'inigrd')
            ENDIF
        ELSE
#endif
            CALL locrint (clstrg,msknew(iadrnew),isiznew,iunit6,iflag)
            IF (iflag .NE. 0) THEN
                CALL prtout('Problem in reading unit =',iunit6,2)
                CALL prcout
     $          ('Problem with array linked to string =',clstrg,1)
                infos = infos + iflag
            ENDIF
#ifdef use_netCDF
        ENDIF
#endif
        IF (nlogprt .GE. 2) THEN
            CALL prtout
     $          ('Read target grid mask for field number =',jf,1)
        ENDIF
C
!$omp end critical
C* Read final mesh surface associated to field jf
C
!$omp critical
        clstrg = clwork(1:icount)//csursuf
        iunit7 = nulsu
#ifdef use_netCDF
        IF (lncdfgrd) THEN
            istatus = 
     $          NF_INQ_VARID(nc_surid, clstrg, il_varid)
            IF (istatus .ne. NF_NOERR .and. lsurf(jf)) THEN
                CALL prtout('Cannot read surface information for fld '
     $              , jf, 1)
                CALL HALTE ('STOP in inigrd')
            ELSE IF (istatus .eq. NF_NOERR) THEN
                IF (n_reaty .eq. NF_FLOAT) THEN
                    istatus = NF_GET_VARA_REAL(nc_surid, 
     $              il_varid, ist, icnt, 
     $                  surnew(iadrnew:iadrnew+isiznew-1))
                    IF (istatus .ne. NF_NOERR .and. lsurf(jf)) THEN
                        CALL prtout('Cannot read target grid surface 
     $ information for field ' , jf, 1)
                        CALL HALTE ('STOP in inigrd')
                    ENDIF
                ELSE IF (n_reaty .eq. NF_DOUBLE) THEN
                    istatus = NF_GET_VARA_DOUBLE(nc_surid, 
     $               il_varid, ist, icnt, 
     $                  surnew(iadrnew:iadrnew+isiznew-1))
                    IF (istatus .ne. NF_NOERR .and. lsurf(jf)) THEN
                        CALL prtout('Cannot read source grid surface 
     $ information for field ' , jf, 1)
                        CALL HALTE ('STOP in inigrd')
                    ENDIF
                ENDIF
            ENDIF
        ELSE
#endif
            CALL locread (clstrg, surnew(iadrnew), isiznew,
     $          iunit7,iflag)
            IF (iflag .NE. 0 .and. lsurf(jf)) THEN
                CALL prtout('Problem in reading unit =',iunit7,2)
                CALL prcout('Problem with array linked to string ',
     $              clstrg,1)
                infos = infos + iflag
            ENDIF
#ifdef use_netCDF
        ENDIF
#endif
!$omp end critical
        IF (nlogprt .GE. 2) THEN
            CALL prtout
     $          ('Read target grid surface for field number ',jf,1)
        ENDIF
 200  CONTINUE
#ifdef use_netCDF
      IF (lncdfgrd) THEN
          CALL hdlerr(NF_CLOSE(nc_grdid), 'inigrd')
          CALL hdlerr(NF_CLOSE(nc_mskid), 'inigrd')
          IF (lg_areas) CALL hdlerr(NF_CLOSE(nc_surid), 'inigrd')
      ENDIF
#endif
C
C* Finish up if problem in reading
C
      IF (infos .NE. 0) CALL HALTE('STOP in inigrd')
C
C
C*    3. End of routine
C        --------------
C
      IF (nlogprt .GE. 1) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '          --------- End of routine inigrd ---------'
          CALL FLUSH (nulou)
      ENDIF
      ENDIF
      RETURN
      END


