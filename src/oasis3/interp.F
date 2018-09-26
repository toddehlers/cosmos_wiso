      SUBROUTINE interp (kindex, kfield)
C****
C               *****************************
C       
C        * OASIS ROUTINE  -  LEVEL 1 *
C               * -------------     ------- *
C               *****************************
C
C**** *interp* - Control routine for interpolation
C
C
C     Purpose:
C     -------
C     Monitor the field interpolation and auxilary analysis
C
C**   Interface:
C     ---------
C       *CALL*  *interp (kindex, kfield)*
C
C     Input:
C     -----
C                kindex : current active fields index array
C                kfield : current active fields total number
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
C     fiasco, blasold, filling, mozaic
C
C     Reference:
C     ---------
C     See OASIS manual (1995)
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       2.0beta   L. Terray      95/09/01  created
C       2.0       L. Terray      96/02/01  modified: mozaic interpolation
C       2.1       L. Terray      96/08/05  modified: Add new arrays for
C                                          mapping data(weight, adresses)
C                                          addition of no interpolation
C       2.3       S. Valcke      99/04/30  added: printing levels
C       2.3       L. Terray      99/09/15  changed: periodicity variables
C       2.5       Gayler/Declat  01/11/09  scrip remapping
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ----------------Include files and USE of modules ---------------------------
C
      USE mod_kinds_oasis
      USE mod_parameter
      USE mod_extrapol
      USE mod_string
      USE mod_analysis
      USE mod_memory
      USE mod_anais
      USE mod_rainbow
      USE mod_unit
      USE mod_printing
C
C* ---------------------------- Argument declarations -------------------
C
      INTEGER (kind=ip_intwp_p) kindex(kfield)
C
C* ---------------------------- Local declarations ----------------------
C
      REAL (kind=ip_realwp_p), DIMENSION(:), ALLOCATABLE :: zbocoef
      INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: iaddr, 
     $    isize, iflag
      CHARACTER(len=8),DIMENSION(:), ALLOCATABLE :: clbofld 
      CHARACTER*8 clintmet, clgrdtyp, clfldtyp, clfilfic, clfilmet
      CHARACTER*8 clfic, clstrg, clname
      CHARACTER*8 clsper, cltper
      LOGICAL llchk
      INTEGER (kind=ip_intwp_p) :: il_sper, il_tper
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
     $    '           ROUTINE interp  -  Level 1'
          WRITE (UNIT = nulou,FMT = *) 
     $    '           **************     *******'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' Monitor field interpolation'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
      ENDIF
C
      ALLOCATE (zbocoef(ig_maxcomb),stat=il_err)
      IF (il_err.NE.0) CALL prtout ('Error in "zbocoef" allocation of 
     $    cookart ',il_err,1)
      ALLOCATE (iaddr(ig_maxcomb),stat=il_err)
      IF (il_err.NE.0) CALL prtout ('Error in "iaddr" allocation of 
     $    cookart ',il_err,1)
      ALLOCATE (isize(ig_maxcomb),stat=il_err)
      IF (il_err.NE.0) CALL prtout ('Error in "isize" allocation of 
     $    cookart ',il_err,1)
      ALLOCATE (iflag(ig_maxcomb),stat=il_err)
      IF (il_err.NE.0) CALL prtout ('Error in "iflag" allocation of 
     $    cookart ',il_err,1)
      ALLOCATE (clbofld(ig_maxcomb),stat=il_err)
      IF (il_err.NE.0) CALL prtout ('Error in "clbofld" allocation of 
     $    cookart ',il_err,1)
      zbocoef(:)=0
      iaddr(:)=0
      isize(:)=0
      iflag(:)=0
      clbofld(:)=' '
C
C
C*    2. Do the job
C        ----------
C
!$omp parallel do default (shared)
!$omp+ private (ja,jf,jj,ji,jd,jb,jc,jk,jl)
!$omp+ private (ifield,clsper,clname,ilataf,ilonaf)
!$omp+ private (il_tper,il_sper,cltper)
!$omp+ private (ipdeb,ilatbf,ilonbf,isiznew)
!$omp+ private (ibofld,zfldcobo,iadrnew,iadrnew_grid)
!$omp+ private (iadrold,iadrold_grid,isizold)
!$omp+ private (isize,iflag)
!$omp+ private (ipointer,iaddr,isiztot,clfic)
!$omp+ private (clfilmet,iunit,ivoisin,clstrg)
!$omp+ private (iloc,clintmet,clgrdtyp,clfldtyp)
!$omp+ private (clfilfic,llchk)

      DO 210 jf = 1, kfield
C
C* Assign local variables
C
        ifield = kindex(jf)
        iadrold = nadrold(ifield)
        iadrold_grid = nadrold_grid(ifield)
        isizold = nsizold(ifield)
        iadrnew = nadrnew(ifield)
        iadrnew_grid = nadrnew_grid(ifield)
        isiznew = nsiznew(ifield)
        ilonbf = nlonbf(ifield)
        ilatbf = nlatbf(ifield)
        ilonaf = nlonaf(ifield)
        ilataf = nlataf(ifield)
        clname = cnaminp(ifield)
        clsper = csper(ifield)
        cltper = ctper(ifield)
        il_sper = nosper(ifield)
        il_tper= notper(ifield)
C
C* Print field name
C
        IF (nlogprt .GE. 1) THEN
            CALL prcout('Treatment of field : ', clname, 2)
        ENDIF
C
C* - Do interpolation
C
        DO 220 ja = 1, ig_ntrans(ifield)
          IF (canal(ja,ifield) .EQ. 'INTERP') THEN
C
C* --->>> Interp
C
C* Assign local variables
              clintmet = cintmet(ifield)
              clgrdtyp = cgrdtyp(ifield)
              clfldtyp = cfldtyp(ifield)
C
C* Zero work array
C
              CALL szero (work, ig_work)
!$omp critical
              CALL fiasco
C* Data about fields, grids, masks and surfaces
     $             (fldnew(iadrnew),
     $             xgrnew(iadrnew_grid), ygrnew(iadrnew_grid), 
     $             surnew(iadrnew_grid),
     $             msknew(iadrnew_grid),ilonaf, ilataf, cltper,il_tper,
     $             fldold(iadrold),
     $             xgrold(iadrold_grid), ygrold(iadrold_grid), 
     $             surold(iadrold_grid),
     $             mskold(iadrold_grid),ilonbf, ilatbf, clsper,il_sper,
C* Work arrays for the different interpolators
     $             work(1), 
     $             work(1+ilonaf), 
     $             work(1+ilonaf+ilataf),
     $             work(1+ilonaf+ilataf+isizold),
     $             work(1+ilonaf+ilataf+isizold+isiznew),
     $             work(1+ilonaf+ilataf+2*isizold+isiznew),
C* Define field number and type, grid and interpolation 
     $             ifield, clintmet, clgrdtyp, clfldtyp)
!$omp end critical
C
C* --->>> Nointerp
C
          ELSE IF (canal(ja,ifield) .EQ. 'NOINTERP') THEN
C
C* output field is equal to input field
C
C* First, check dimension
C
                llchk = ilataf .EQ. ilatbf .AND.
     $              ilonaf-il_tper .EQ. ilonbf-il_sper
                IF (.NOT. llchk) CALL prcout('WARNING: size mismatch
     $              in NOINTERP between old and new field ',clname,2)
                IF (.NOT. llchk) CALL HALTE('STOP in interp')
C* Do the assign
                DO 230 jj = 1, ilatbf
                  DO 233 ji = 1, ilonbf-il_sper
                    jk = (jj-1)*ilonaf+ji
                    jl = (jj-1)*ilonbf+ji
                    fldnew(iadrnew - 1 + jk) = fldold(iadrold - 1 + jl)
 233              CONTINUE 
                  DO 236 ji = 1, il_tper
                    jk = (jj-1)*ilonaf+ilonaf-il_tper+ji
                    jl = (jj-1)*ilonbf+ji
                    fldnew(iadrnew - 1 + jk) = fldold(iadrold - 1 + jl)
 236              CONTINUE 
 230            CONTINUE 
C
C* --->>> Mozaic
C
            ELSE IF (canal(ja,ifield) .EQ. 'MOZAIC') THEN
C
C* assign local variables and get pointer for mapping interpolation
C
                clfic = cgrdmap(ifield)
                iunit = nlumap(ifield)
                iloc = nmapfl(ifield)
                ipdeb = (nmapfl(ifield)-1)*ig_maxmoa*ig_maxgrd+1
                ivoisin = nmapvoi(ifield)
                clstrg = cficbf(ifield)//cficaf(ifield)
                CALL mozaic (fldnew(iadrnew), isiznew,
     $                       fldold(iadrold), isizold,
     $                       clfic, iunit, clstrg, iloc,
     $                       amapp(ipdeb), nmapp(ipdeb), 
     $                       ivoisin, lmapp(ifield))
C* --->>> SCRIP remapping
C
#ifndef key_noSCRIP
            ELSE IF (canal(ja,ifield) .EQ. 'SCRIPR') THEN
C
C* perform the scrip remapping
C
               clgrdtyp = cgrdtyp(ifield)
C
C* Vector case : VECTOR_I 
C
                IF (cfldtype(ifield) .eq. 'VECTOR_I') THEN

                   ifield_J = ig_assoc_input_field(ifield)
                   iadrold_J = nadrold(ifield_J)
                   iadrold_grid_J = nadrold_grid(ifield_J)
                   iadrnew_J = nadrnew(ifield_J)
                   isizold_J = nsizold(ifield_J)
                   isiznew_J = nsiznew(ifield_J)
                   ilonbf_J = nlonbf(ifield_J)
                   ilatbf_J = nlatbf(ifield_J)
                   iadrnew_grid_J = nadrnew_grid(ifield_J)

                   CALL scriprmp_vector (
     $                  fldnew(iadrnew), fldnew(iadrnew_J), 
     $                  fldold(iadrold), fldold(iadrold_J),  
     $                  cficbf(ifield),  cficbf(ifield_J), 
     $                  cficaf(ifield),  cficaf(ifield_J), 
     $                  isizold, isiznew, 
     $                  mskold(iadrold_grid), mskold(iadrold_grid_J), 
     $                  msknew(iadrnew_grid), msknew(iadrnew_grid_J),
     $                  xgrold(iadrold_grid),  ygrold(iadrold_grid),
     $                  xgrold(iadrold_grid_J),ygrold(iadrold_grid_J),
     $                  xgrnew(iadrnew_grid),  ygrnew(iadrnew_grid), 
     $                  xgrnew(iadrnew_grid_J),ygrnew(iadrnew_grid_J), 
     $                  ilonbf, ilatbf, 
     $                  ilonaf, ilataf, 
     $                  cmap_method(ifield), clgrdtyp, 
     $                  il_sper, il_tper, 
     $                  clsper, cltper,
     $                  cnorm_opt(ifield), corder(ifield), 
     $                  crsttype(ifield), nbins(ifield),
     $                  lextrapdone(ifield), lrotate(ifield),
     $                  varmul(ifield), nscripvoi(ifield))
                   
C
C* Vector case : VECTOR_J 
C
                ELSE IF (cfldtype(ifield) .eq. 'VECTOR_J') THEN

                   WRITE (UNIT = nulou,FMT = *) 
     $    ' This is a VECTOR_J field. It is calculated at the same '
                   WRITE (UNIT = nulou,FMT = *) 
     $    ' time as the corresponding VECTOR_I field'

                ELSE

C
C* Scalar case
C
                     CALL scriprmp (
     $                   fldnew(iadrnew), fldold(iadrold), 
     $                   isizold, isiznew, 
     $                   mskold(iadrold_grid), msknew(iadrnew_grid),
     $                   xgrold(iadrold_grid), ygrold(iadrold_grid), 
     $                   ilonbf, ilatbf, 
     $                   xgrnew(iadrnew_grid), ygrnew(iadrnew_grid), 
     $                   ilonaf, ilataf, 
     $                   cmap_method(ifield), clgrdtyp, 
     $                   il_sper, il_tper, clsper, cltper,
     $                   cficbf(ifield), cficaf(ifield), 
     $                   cnorm_opt(ifield), 
     $                   corder(ifield), crsttype(ifield), 
     $                   nbins(ifield),
     $                   lextrapdone(ifield), varmul(ifield), 
     $                   nscripvoi(ifield)) 
                 END IF
#endif
C
C* --->>> Blasold
C
            ELSE IF (canal(ja,ifield) .EQ. 'BLASOLD') THEN
C
C* Assign local variables
C
              zfldcobo = afldcobo(ifield)
              ibofld = nbofld(ifield)
              DO 240 jc = 1, ibofld
                clbofld(jc) = cbofld(jc,ifield)
                zbocoef(jc) = abocoef(jc,ifield)
 240          CONTINUE
C
C* Get the additional fields (pointers and sizes)
C
              CALL szero( work, ig_work)
              DO 250 jc = 1, ibofld
                IF (clbofld(jc) .EQ. 'CONSTANT') THEN
                    isize(jc) = isizold
                  ELSE 
                    DO 260 jb = 1, ig_nfield
C
C* Check field names input list
C
                      IF (clbofld(jc) .EQ. cnaminp(jb)) THEN
                          iflag(jc) = jb
                      ENDIF 
 260                CONTINUE
                    ipointer  = nadrold(iflag(jc))
                    isize(jc) = nsizold(iflag(jc))
                ENDIF 
                IF (jc .EQ. 1) THEN
                    iaddr(jc) = 1
                  ELSE
                    iaddr(jc) = 1 + isize(jc-1)
                ENDIF
C
C* Assign values to temporary array work
C
                IF (clbofld(jc) .EQ. 'CONSTANT') THEN
                    DO 270 jd = 1, isize(jc)
                      work(iaddr(jc)+jd-1) = 1.0
 270                CONTINUE 
                  ELSE 
                    DO 280 jd = 1, isize(jc)
                      work(iaddr(jc)+jd-1) = fldold(ipointer+jd-1)
 280                CONTINUE 
                ENDIF 
 250          CONTINUE
C
C* Get total size for array work ( sum of additional fields sizes)
C 
              isiztot = iaddr(ibofld) + isize(ibofld) - 1
              CALL blasold (fldold(iadrold), isizold, ifield,
     $                      zfldcobo, ibofld, iaddr, isize,
     $                      zbocoef, isiztot, work)
C
C* --->>> Filling
C
            ELSE IF (canal(ja,ifield) .EQ. 'FILLING') THEN
C
C* Assign local variables
C
              clfilfic = cfilfic(ifield)
              iunit = nlufil(ifield)
              clfilmet = cfilmet(ifield)
C
C* Zero work array
C
              CALL szero (work, ig_work)
C
C* Address of overlapping grids array used in Anaism 
C
              ipointer = (naismfl(ifield)-1)*ig_maxgrd + 1
              CALL filling (fldnew(iadrnew), work(1), work(isiznew+1),
     $                      work(2*isiznew+1),
     $                      xgrnew(iadrnew_grid), ygrnew(iadrnew_grid),
     $                      ilonaf, ilataf, 
     $                      msknew(iadrnew_grid), nmesh(ipointer), iunit,
     $                      clfilfic, clfilmet)
            ELSE
              CONTINUE
          END IF
 220    CONTINUE
 210  CONTINUE 
C
C
C*    3. Deallocation and end of routine
C        -------------------------------
C
      DEALLOCATE(zbocoef)
      DEALLOCATE(iaddr)
      DEALLOCATE(isize)
      DEALLOCATE(iflag)
      DEALLOCATE(clbofld)

      IF (nlogprt .GE. 1) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '          --------- End of routine interp ---------'
          CALL FLUSH (nulou)
      ENDIF
      RETURN
      END
