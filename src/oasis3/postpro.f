      SUBROUTINE postpro (kindex, kfield)
      USE mod_kinds_oasis
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 1 *
C               * -------------     ------- *
C               *****************************
C
C**** *postpro* - postprocessing routine
C
C
C     Purpose:
C     -------
C     Do the field postprocessing
C
C**   Interface:
C     ---------
C       *CALL*  *postpro (kindex, kfield)*
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
C     reverse, revmsk, masq, extrap, glored, chkfld
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
C       2.1       L. Terray      96/09/25  modified: call to chkfld and
C                                          addition of amskred
C       2.2       L. Terray      97/12/31  modified: call to extrap
C       2.3       L. Terray      99/03/01  modified: call to extrap
C       2.3       S. Valcke      99/04/15  modified: CALL to extrap
C       2.3       S. Valcke      99/04/30  added: printing levels
C       2.3       L. Terray      99/09/15  changed periodicity variables
C       2.3       S. Valcke      99/10/14  CALL to extrap corrected
C       2.5       S. Valcke      00/09/05  Changed iintflx for itoutflx 
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------- Include files and USE of modules---------------------------
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
      CHARACTER*8 clxordaf, clyordaf, clextmet, clname, clper
      CHARACTER*32 clabel
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Initialization
C        --------------
C
      IF (nlogprt .GE. 1) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '           ROUTINE postpro  -  Level 1'
          WRITE (UNIT = nulou,FMT = *) 
     $    '           ***************     *******'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    ' Postprocessing of coupling fields'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
      ENDIF
C
C
C*    2. Do the job
C        ----------
C
!$omp parallel do default (shared)
!$omp+ private (jf,ifield,iadrnew,iadrnew_grid)
!$omp+ private (isiznew,clname,clyordaf)
!$omp+ private (ilabel,clabel,ilonaf,ilataf,clxordaf)
!$omp+ private (itoutflx,ja,ji,clextmet,clper)
!$omp+ private (itronca,ineibor,iper)

      DO 210 jf = 1, kfield
C
C* Assign local variables
C
        ifield = kindex(jf)
        iadrnew = nadrnew(ifield)
        iadrnew_grid = nadrnew_grid(ifield)
        isiznew = nsiznew(ifield)
        clname = cnamout(ifield)
        ilabel = numlab(ifield)
        clabel = cfldlab(ilabel)
        ilonaf = nlonaf(ifield)
        ilataf = nlataf(ifield)
        clxordaf = cxordaf(ifield)
        clyordaf = cyordaf(ifield)
        itoutflx = ntoutflx(ifield)
C
C* Print field name
C
        IF (nlogprt .GE. 1) THEN
            CALL prcout('Treatment of field : ', clname, 2)
        ENDIF
C
C* - Do postprocessing analysis
C
        DO 220 ja = 1, ig_ntrans(ifield)
C
C* --->>> MaskP
C
          IF (canal(ja,ifield) .EQ. 'MASKP') THEN 
            
              CALL masq(fldnew(iadrnew), isiznew, amskvalnew(ifield),
     $                      msknew(iadrnew_grid))
C
C* --->>> Reverse
C
          ELSE IF (canal(ja,ifield) .EQ. 'REVERSE') THEN 
              CALL reverse (fldnew(iadrnew), ilonaf,
     $                      ilataf, clxordaf, clyordaf)
C
C* --->>> Checkout: perform basic checks on the field to be exported
C
          ELSE IF (canal(ja,ifield) .EQ. 'CHECKOUT') THEN
              CALL chkfld(clname, clabel, 
     $            fldnew(iadrnew), msknew(iadrnew_grid), 
     $            surnew(iadrnew_grid),
     $            isiznew, ilonaf, itoutflx)
C
C* --->>> Glored
C
          ELSE IF (canal(ja,ifield) .EQ. 'GLORED') THEN
C
C* Do extrapolation on full grid to assign sea values to land points
C
C* - First we mask the field
C    We use a predefined value for the mask
C    If necessary, we reorder the mask as the field might have been already
C    reversed ( while array msknew is ordered along OASIS conventions). 
C
C* Put mask in work array
C
              CALL izero (nwork, ig_nwork)
              DO 230 ji = 1, isiznew
                nwork(ji) = msknew(iadrnew_grid + ji -1)
 230          CONTINUE 
C
C* Reverse mask if necessary
C
              CALL revmsk (nwork(1), ilonaf, ilataf,
     $            clxordaf, clyordaf)
              zmskval = amskred
              CALL masq (fldnew(iadrnew), isiznew, zmskval,
     $            nwork(1))
C
C* - Then we extrapolate
C    We use predefined values for extrapolation parameters
C
              clextmet = 'NINENN'
              ineibor = neighborg(ifield)
C
C* Grid periodicity
C
              clper = ctper(ifield)
              iper = notper(ifield)
C
C* Zero work array
C
              CALL szero (work, ig_work)
C
C* Do it now
C
              CALL extrap (fldnew(iadrnew), zmskval, work(1), 
     $            nwork(1), ilonaf, ilataf,
     $            ineibor, clextmet, clper, iper,
     $            niwtng(ifield), nninnflg(ifield))
C
C* Do the interpolation full to reduced gaussian grid
C
              itronca = ntronca(ifield)
              CALL szero (work, ig_work)
              CALL glored (fldnew(iadrnew), work(1),
     $                     ilonaf, ilataf, itronca)
          ELSE
              CONTINUE
          END IF
 220    CONTINUE
 210  CONTINUE 
C
C
C*    3. End of routine
C        --------------
C
      IF (nlogprt .GE. 1) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '          --------- End of routine postpro ---------'
          CALL FLUSH (nulou)
      ENDIF
      RETURN
      END
