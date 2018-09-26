      SUBROUTINE preproc (kindex, kfield)
      USE mod_kinds_oasis
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 1 *
C               * -------------     ------- *
C               *****************************
C
C**** *preproc* - preprocess routine
C
C
C     Purpose:
C     -------
C     Do the field preprocessing
C
C**   Interface:
C     ---------
C       *CALL*  *preproc (kindex, kfield)*
C
C     Input:
C     -----
C                kindex : field identificator array (integer 1D)
C                kfield : number of fields for current iteration (integer)
C
C     Output:
C     ------
C     None
C
C     Workspace:
C     ---------
C                zcocoef : additional field coefficients for correct (real 1D)
C                inip    : array for reduced grid data for redglo (integer 1D)
C                clcofld : additional field names for correct (character 1D)
C                clcofic : array to handle data file names (character 1D)
C                iunit   : array to handle I/O units of data files (integer 1D)
C
C     Externals:
C     ---------
C     correct, extrap, extraw, invert, masq, redglo, chkfld
C
C     Reference:
C     ---------
C     See OASIS manual (1998)
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       2.0       L. Terray      95/09/01  created
C       2.1       L. Terray      96/09/25  modified: Call to chkfld
C       2.2       L. Terray      97/12/16  Added: new extrapolation
C                                          and change call to extrap
C       2.3       L. Terray      99/03/01  modified: call to extrap
C       2.3       S. Valcke      99/03/16  modified for T213 and T319
C       2.3       S. Valcke      99/03/16  modified for T213 and T319
C       2.3       S. Valcke      99/03/26  changed troncature for number of 
C                                          latitude between equator and pole
C       2.3       S. Valcke      99/03/30  changed arguments in CALL to extrap
C       2.3       S. Valcke      99/04/30  added: printing levels
C       2.3       L. Terray      99/09/15  changed periodicity variables
C       2.5       S. Valcke      00/09/05  Changed iintflx for itinpflx 
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* -----------------Include files and USE of modules---------------------------
C
      USE mod_parameter
      USE mod_string
      USE mod_analysis
      USE mod_memory
      USE mod_extrapol
      USE mod_unit
      USE mod_gauss
      USE mod_label
      USE mod_printing
C
C* ---------------------------- Argument declarations -------------------
C
      INTEGER (kind=ip_intwp_p) kindex(kfield)
C
C* ---------------------------- Local declarations ----------------------
C
      INTEGER (kind=ip_intwp_p) inip(320)
      INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: iunit
      REAL (kind=ip_realwp_p),DIMENSION(:), ALLOCATABLE :: zcocoef
      CHARACTER*8 clxordbf, clyordbf, clextmet, clname, clmsk, clper
      CHARACTER*8 clfic, clstrg
      CHARACTER(len=8), DIMENSION(:), ALLOCATABLE :: clcofic, clcofld
      CHARACTER*32 clabel
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Initialization and allocation of local arrays
C        ---------------------------------------------
C
      IF (nlogprt .GE. 1) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '           ROUTINE preproc  -  Level 1'
          WRITE (UNIT = nulou,FMT = *) 
     $    '           ***************     *******'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    ' Preprocessing of coupling fields'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
      ENDIF
C
      ALLOCATE (iunit(ig_maxcomb),stat=il_err)
      IF (il_err.NE.0) CALL prtout ('Error in "iunit" allocation of 
     $    cookart ',il_err,1)
      ALLOCATE (zcocoef(ig_maxcomb),stat=il_err)
      IF (il_err.NE.0) CALL prtout ('Error in "zcocoef" allocation of 
     $    cookart ',il_err,1)
      ALLOCATE (clcofic(ig_maxcomb),stat=il_err)
      IF (il_err.NE.0) CALL prtout ('Error in "clcofic" allocation of 
     $    cookart ',il_err,1)
      ALLOCATE (clcofld(ig_maxcomb),stat=il_err)
      IF (il_err.NE.0) CALL prtout ('Error in "clcofld" allocation of 
     $    cookart ',il_err,1)
      iunit(:)=0
      zcocoef(:)=0
      clcofic(:)=' '
      clcofld(:)=' '
C
C* Zeroes work array
C
      CALL szero (work,ig_work)
      CALL izero (nwork,ig_nwork)
C
C
C*    2. Do the job
C        ----------
!$omp parallel do default (shared)
!$omp+ private (ja,jf,ji,jc)
!$omp+ private (ifield,clname,ilataf)
!$omp+ private (clabel,clextmet,ineibor)
!$omp+ private (ipdeb,ilatbf,ilonbf)
!$omp+ private (iadrold,iadrold_grid,isizold)
!$omp+ private (clfic,ilun,clper)
!$omp+ private (clstrg,ilabel)
!$omp+ private (iloc,clintmet,clgrdtyp,clfldtyp)
!$omp+ private (icofld)
!$omp+ private (clfilfic,itinpflx,zfldcoef)
!$omp+ private (iper,itronca,clmsk,iredu)
!$omp+ private (zmskval,clxordbf,clyordbf)

C
      DO 210 jf = 1, kfield
C
C* Assign local variables
C
        ifield = kindex(jf)
        ilabel = numlab(ifield)
        clname = cnaminp(ifield)
        clabel = cfldlab(ilabel)
        iadrold = nadrold(ifield)
        iadrold_grid = nadrold_grid(ifield)
        isizold = nsizold(ifield)
        ilonbf = nlonbf(ifield)
        ilatbf = nlatbf(ifield)
        itinpflx = ntinpflx(ifield)
C
C* Print field name
C
        IF (nlogprt .GE. 1) THEN
            CALL prcout('Treatment of field :', clname, 2)
        ENDIF
C
C* - Do preprocessing analysis
C
        DO 220 ja = 1, ig_ntrans(ifield)
          IF (canal(ja,ifield) .EQ. 'MASK') THEN
C
C* --->>> Mask
C
              zmskval = amskval(ifield)
!$omp critical
              CALL masq (fldold(iadrold), isizold, zmskval,
     $                   mskold(iadrold_grid))
!$omp end critical
C
C* --->>> Invert
C
            ELSE IF (canal(ja,ifield) .EQ. 'INVERT') THEN 
              clxordbf = cxordbf(ifield)
              clyordbf = cyordbf(ifield)
              CALL invert (fldold(iadrold), ilonbf,
     $                     ilatbf, clxordbf, clyordbf)
C
C* --->>> Checkin: perform basic checks on input field
C
            ELSE IF (canal(ja,ifield) .EQ. 'CHECKIN') THEN 
              CALL chkfld(clname, clabel, 
     $            fldold(iadrold), mskold(iadrold_grid), 
     $              surold(iadrold_grid),
     $            isizold, ilonbf, itinpflx)
C
C* --->>> Flux correction
C
            ELSE IF (canal(ja,ifield) .EQ. 'CORRECT') THEN
C
C* Assign local variables to main field coefficient
C
              zfldcoef = afldcoef(ifield)
C
C* Get loop index to read additional fields, coefficients, filenames
C  and related logical units
C
              icofld = ncofld(ifield)
              DO 230 jc = 1, icofld
                clcofld(jc) = ccofld(jc,ifield)
                zcocoef(jc) = acocoef(jc,ifield)
                clcofic(jc) = ccofic(jc,ifield)
                iunit(jc)   = nludat(jc,ifield)
 230          CONTINUE
C
C* Zero work array
C
              CALL szero (work, ig_work)
C
C* Do the job
C 
              CALL correct (fldold(iadrold), isizold,
     $                      zfldcoef, icofld, zcocoef(1),
     $                      work(1), iunit(1), clcofic(1), 
     $                      clcofld(1))
C
C* --->>> Extrap
C
            ELSE IF (canal(ja,ifield) .EQ. 'EXTRAP') THEN
              clextmet = cextmet(ifield)
              ineibor = neighbor(ifield)
C
C* 8-nearest neighbors extrapolation
C
              IF (clextmet .EQ. 'NINENN') THEN 
                  zmskval = amskval(ifield)
                  clper = csper(ifield)
                  iper = nosper(ifield)
C
C* Zero work array
C
                  CALL szero (work, ig_work)
C
C* Do it now
C
#ifdef key_openmp
                  niwtn=0
#endif
!$omp critical
                  CALL extrap (fldold(iadrold), zmskval, work(1),
     $                         mskold(iadrold_grid), ilonbf, ilatbf, 
     $                         ineibor, clextmet, clper, iper,  
     $                         niwtn(ifield), nninnfl(ifield))
C
!$omp end critical
C* N-weighted neighbors extrapolation
C
              ELSE IF (clextmet .EQ. 'WEIGHT') THEN 
                  clfic = cgrdext(ifield)
                  ilun  = nluext(ifield)
                  iloc  = nextfl(ifield)
                  clstrg = cficbf(ifield)//cficbf(ifield)
                  ipdeb = (nextfl(ifield)-1)*ig_maxext*ig_maxgrd+1
C
C* Do it now
C
                  CALL extraw (fldold(iadrold), mskold(iadrold_grid),
     $                isizold, clfic, ilun, clstrg, iloc,  
     $                aextra(ipdeb), nextra(ipdeb), ineibor,
     $                lextra(ifield))
              ENDIF 
              lextrapdone(ifield) = .true.
C
C* --->>> Redglo
C
            ELSE IF (canal(ja,ifield) .EQ. 'REDGLO') THEN
              itronca = ntronca(ifield)
              clmsk = cmskrd(ifield)
              zmskval = amskred
C
C* get number of longitudes by latitude circle and total number of points
C  reduced grid
C
C* Zero work array
C
              CALL szero (work, ig_work)
              CALL izero (nwork, ig_nwork)
              CALL izero (inip, 320)
              IF (itronca .EQ. 16) THEN
                  DO 240 ji = 1, itronca
                    inip(ji) = ninip16(ji)
 240              CONTINUE
                  iredu = nredu16
                ELSE IF (itronca .EQ. 24)  THEN 
                  DO 250 ji = 1, itronca
                    inip(ji) = ninip24(ji)
 250              CONTINUE 
                  iredu = nredu24
                ELSE IF (itronca .EQ. 32)  THEN 
                  DO 260 ji = 1, itronca
                    inip(ji) = ninip32(ji)
 260              CONTINUE 
                  iredu = nredu32
                ELSE IF (itronca .EQ. 48)  THEN 
                  DO 270 ji = 1, itronca
                    inip(ji) = ninip48(ji)
 270              CONTINUE
                   iredu = nredu48   
                ELSE IF (itronca .EQ. 80)  THEN 
                  DO 280 ji = 1, itronca 
                    inip(ji) = ninip80(ji)
 280              CONTINUE
                  iredu = nredu80
                ELSE IF (itronca .EQ. 160)  THEN 
                  DO 285 ji = 1, itronca
                    inip(ji) = ninip160(ji)
 285              CONTINUE
                  iredu = nredu160
                ELSE
                  CALL prtout
     $          ('WARNING!!! Oasis cannot treat this grid with 2*NO 
     $          latitude lines with NO = ', itronca, 2)
                  CALL prtout
     $                ('Implement data for NO =', itronca, 2)
                  CALL HALTE('STOP in preproc')
              ENDIF
              DO 290 ji = ilatbf/2 + 1, ilatbf
                inip(ji) = inip(ilatbf - ji + 1)
 290          CONTINUE
              CALL redglo (fldold(iadrold), work(1), iredu, inip,
     $                     ilonbf, ilatbf, nwork(1), nwork(1+iredu),
     $                     itronca, zmskval, clmsk)
            ELSE
              CONTINUE 
          END IF
 220      CONTINUE
 210    CONTINUE 
C
C
C*    3. Deallocation and end of routine
C        -------------------------------
C
        DEALLOCATE(iunit)
        DEALLOCATE(zcocoef)
        DEALLOCATE(clcofic)
        DEALLOCATE(clcofld)
C
      IF (nlogprt .GE. 1) THEN 
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '          --------- End of routine preproc ---------'
          CALL FLUSH (nulou)
      ENDIF
      RETURN
      END
