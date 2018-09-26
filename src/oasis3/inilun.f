      SUBROUTINE inilun
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 0 *
C               * -------------     ------- *
C               *****************************
C
C**** *inilun*  - Initialize logical unit numbers
C
C     Purpose:
C     -------
C     Creates and prints logical unit numbers used to deal with
C     grids, masks and surfaces files as well as anais-related files
C
C**   Interface:
C     ---------
C       *CALL*  *inilun*
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
C       2.0       L. Terray      95/08/23  modified: new structure
C       2.2       L. Terray      97/10/10  added: unit nudum for SVIPC
C       2.3       S. Valcke      99/03/30  added: unit nulgn for NINENN
C       2.3       S. Valcke      99/04/30  added: printing levels
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------- Include files and USE of modules---------------------------
C
      USE mod_parameter
      USE mod_string
      USE mod_unit
      USE mod_printing
      USE mod_hardware
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Assign unit numbers and initialize comlun
C        -----------------------------------------
C
C* First we open output file for coupler
C
      IF (nlogprt .GE. 1) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '           ROUTINE inilun  -  Level 0'
          WRITE (UNIT = nulou,FMT = *) 
     $    '           **************     *******'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    ' Set up logical unit numbers'
          WRITE (UNIT = nulou,FMT = *) ' '
      ENDIF
C
C* If there is no field going through Oasis then assign only the trace file 
C  unit for CLIM
      IF (.not. lg_oasis_field) THEN
         nultr = 7
         IF (nlogprt .GE. 1)
     $        WRITE (UNIT = nulou,FMT = *)'Trace file unit for CLIM :',
     $        nultr
      ELSE
C* Grids file
         nulgr = 11
C* Masks file
         nulma = 12
C* Surfaces file
         nulsu = 13
C* File for reduced grid masks
         nulrd = 14
C* Trace file for CLIM and PVM
         nultr = 7
C* Output file for ANAIS interpolation
         nulan = 8
C* Dummy file for SVIPC library
         nudum = 9
C* Anaism weights file
         nulcc = 16
C* Anaisg weights file
         nulgg = 17
C* NINENN weight and address file
         nulgn = 18
C     
C     
C*    2. Print comlun
C        ------------
C     
         IF (nlogprt .GE. 1) THEN 
            WRITE (UNIT = nulou,FMT ='(
     $           '' nulin ='',i3,'' nulou ='',i3,
     $           '' nulgr ='',i3,'' nulma ='',i3,
     $           '' nulsu ='',i3,'' nultr ='',i3,   
     $           '' nulcc ='',i3,'' nulgg ='',i3,'' nulgn ='',i3,
     $           '' nulan ='',i3,'' nulrd ='',i3,'' nudum ='',i3,/)')
     $           nulin, nulou, nulgr, nulma, nulsu, 
     $           nultr, nulcc, nulgg, nulgn, nulan, nulrd, nudum
         ENDIF 
C     
C     
C     
C*    3. Assign unit numbers to input and output binary files
C        ----------------------------------------------------
C     
         nluinp(1)=21
         DO 310 jf = 2, ig_nfield
            isamefic=0
            DO 320 jj = 1, jf-1
               IF (cficinp(jf) .eq. cficinp(jj)) THEN
                  isamefic=1
                  nluinp(jf) = nluinp(jj)
               ENDIF
 320        CONTINUE
            IF (isamefic .lt. 1) nluinp(jf) = 20 + jf
 310     CONTINUE
C     
C*    For PIPE technique only
C     
         IF (cchan. eq. 'PIPE' .or. cchan. eq. 'NONE') THEN
            nluout(1)=21 + ig_nfield
            DO 330 jf = 2, ig_nfield
               DO 340 jj = 1, jf-1
                  IF(cficout(jf) .eq. cficout(jj)) THEN
                     nluout(jf) = nluout(jj)
                  ELSE
                     nluout(jf) = 20 + ig_nfield + jf
                  ENDIF
 340           CONTINUE
 330        CONTINUE
         ENDIF
      ENDIF
C
C*     4. End of routine
C         --------------
C
      IF (nlogprt .GE. 1) THEN 
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '          --------- End of routine inilun ---------'
          CALL FLUSH (nulou)
      ENDIF 
      RETURN
      END






