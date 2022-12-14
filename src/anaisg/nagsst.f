      SUBROUTINE nagsst (pflda, kmska, kvmska, kngxa, kngya,
     $                   prbtoa, kbtoa, kwbtoa,
     $                   pfldb, kngxb, kngyb)
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 3 *
C               * -------------     ------- *
C               *****************************
C
C**** *nagsst* -  Interpolation  Anais-global  without constraints 
C
C     Purpose:
C     -------
C     Interpolate with a nearest neighbor method
C
C**   Interface:
C     ---------
C       *CALL*  *nagsst(pflda, kmska, kvmska, kngxa, kngya,
C                       prbtoa, kbtoa, kwbtoa,
C                       pfldb, kngxb, kngyb)*
C
C     Input:
C     -----
C                kmska  : mask for target grid (integer 2D)
C                kvmska : the value of the mask for target grid
C                kngxa  : number of longitudes for target grid
C                kngya  : number of latitudes for target grid
C                prbtoa : weights for Anaisg interpolation (real 2D)
C                kbtoa  : source grid neighbors adresses (integer 2D)
C                kwbtoa : maximum number of nearest neighbors
C                pfldb  : field on source grid (real 2D)
C                kngxb  : number of longitudes for source grid
C                kngyb  : number of latitudes for source grid
C
C     Output:
C     ------
C                pflda: field on target grid (real 2D)
C
C     Workspace:
C     ---------
C     None
C
C     External:
C     --------
C     qlsst
C
C     References:
C     ----------
C     O. Thual, Simple ocean-atmosphere interpolation. 
C               Part A: The method, EPICOA 0629 (1992)
C               Part B: Software implementation, EPICOA 0630 (1992)
C     See also OASIS manual (1995)
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      ----------- 
C       1.1       O. Thual       93/04/15  created 
C       2.0       L. Terray      95/10/01  modified: new structure
C       2.3       S. Valcke      99/04/30  added: printing levels
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------------------- Include files ---------------------------
C
      USE mod_kinds_oasis
      USE mod_unit
      USE mod_printing
C
C* ---------------------------- Argument declarations -------------------
C     
      REAL (kind=ip_realwp_p) pflda(kngxa,kngya), pfldb(kngxb,kngyb)
      REAL (kind=ip_realwp_p) prbtoa(kwbtoa,kngxb*kngyb)
      INTEGER (kind=ip_intwp_p) kmska(kngxa,kngya), 
     $    kbtoa(kwbtoa,kngxb*kngyb)
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Initializations and checkings
C        -----------------------------
C
      IF (nlogprt .GE. 2) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '           ROUTINE nagsst  -  Level 3'
          WRITE (UNIT = nulou,FMT = *) 
     $    '           **************     *******'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    ' Does Anais-global interpolation'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
      ENDIF
C
C* Define global dimensions
C
      inga = kngxa * kngya
      ingb = kngxb * kngyb
C
C
C*    2. Interpolation
C        -------------
C
      CALL qlsst (pflda, prbtoa, kbtoa, kwbtoa, inga, pfldb,
     $            ingb, kmska, kvmska)
C
C
C*    3. End of routine
C        --------------
C
      IF (nlogprt .GE. 2) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '          --------- End of routine nagsst ---------'
          CALL FLUSH (nulou)
      ENDIF
      RETURN
      END
