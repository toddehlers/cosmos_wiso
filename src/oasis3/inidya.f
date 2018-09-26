      SUBROUTINE inidya
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 0 *
C               * -------------     ------- *
C               *****************************
C
C**** *inidya*  - Initialize simulated dynamic allocation
C
C     Purpose:
C     -------
C     Set up pseudo-dynamic allocation of arrays.
C     Define adresses and lengths of small arrays in big arrays
C
C**   Interface:
C     ---------
C       *CALL*  *inidya*
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
C     ineed
C
C     Externals:
C     ---------
C     chkmem
C
C     Reference:
C     ---------
C     See OASIS manual (1995) 
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       2.0       L. Terray      95/08/23  created 
C       2.3       S. Valcke      99/04/30  added: printing levels
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------- Include files and USE of modules---------------------------
C
      USE mod_kinds_oasis
      USE mod_parameter
      USE mod_string
      USE mod_memory
      USE mod_unit
      USE mod_printing
C
C* ---------------------------- Local declarations ----------------------
C
      INTEGER (kind=ip_intwp_p) il_memtot(2)
      INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: 
     $     nsizold_grid_aux, nsiznew_grid_aux, nadrold_grid_aux, 
     $     nadrnew_grid_aux
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
C*    1. Set up memory dynamic allocation
C     --------------------------------
C     
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) 
     $           '           ROUTINE inidya  -  Level 0'
            WRITE (UNIT = nulou,FMT = *) 
     $           '           **************     *******'
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) 
     $           ' Set up memory dynamic allocation'
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) ' '
         ENDIF
C
C* Allocation of local arrays
C
         ALLOCATE (nsizold_grid_aux(maxval(ig_grid_nbrbf)))
         nsizold_grid_aux(:)=0
         ALLOCATE (nsiznew_grid_aux(maxval(ig_grid_nbraf)))
         nsiznew_grid_aux(:)=0
         ALLOCATE (nadrold_grid_aux(maxval(ig_grid_nbrbf)))
         nadrold_grid_aux(:)=0
         ALLOCATE (nadrnew_grid_aux(maxval(ig_grid_nbraf)))
         nadrnew_grid_aux(:)=0
C     
C* Zero size and adress arrays
C     
         CALL izero (nsizold, ig_nfield)
         CALL izero (nsiznew, ig_nfield)
         CALL izero (nadrold, ig_nfield)
         CALL izero (nadrnew, ig_nfield)
         CALL izero (nadrold_grid, ig_nfield)
         CALL izero (nadrnew_grid, ig_nfield)
C     
C* Zero main storage area
C     
C* - Integer arrays: masks
C     
         CALL izero (mskold, ig_maxold_grid)
         CALL izero (msknew, ig_maxnew_grid)
C     
C* - Real arrays: fields, grids, surfaces
C     
         CALL szero (fldold, ig_maxold)
         CALL szero (fldnew, ig_maxnew)
         CALL szero (xgrold, ig_maxold_grid)
         CALL szero (ygrold, ig_maxold_grid)
         CALL szero (xgrnew, ig_maxnew_grid)
         CALL szero (ygrnew, ig_maxnew_grid)
         CALL szero (surold, ig_maxold_grid)
         CALL szero (surnew, ig_maxnew_grid)

C     
C* Get the sizes for each small array
C     
         DO 110 jf = 1, ig_nfield
            nsizold(jf) = nlonbf(jf) * nlatbf(jf)
            nsiznew(jf) = nlonaf(jf) * nlataf(jf)
            nsizold_grid_aux (ig_grid_nbrbf(jf)) = 
     $           nlonbf(jf) * nlatbf(jf)
            nsiznew_grid_aux (ig_grid_nbraf(jf)) = 
     $           nlonaf(jf) * nlataf(jf)
 110     CONTINUE
C     
C* Get the pointers for each small array
C     
         nadrold(1) = 1
         nadrnew(1) = 1
         nadrold_grid_aux(1) = 1
         nadrnew_grid_aux(1) = 1
         nadrold_grid(1) = 1
         nadrnew_grid(1) = 1
         
         DO ib = 2, maxval(ig_grid_nbrbf)
            nadrold_grid_aux(ib) = nadrold_grid_aux(ib-1) + 
     $           nsizold_grid_aux(ib-1)
         ENDDO
         DO ib = 2, maxval(ig_grid_nbraf)
            nadrnew_grid_aux(ib) = nadrnew_grid_aux(ib-1) + 
     $           nsiznew_grid_aux(ib-1)
         ENDDO

         DO 120 jf = 2, ig_nfield
            nadrold(jf) = nadrold(jf-1) + nsizold(jf-1)
            nadrnew(jf) = nadrnew(jf-1) + nsiznew(jf-1)
            nadrold_grid(jf) = nadrold_grid_aux(ig_grid_nbrbf(jf))
            nadrnew_grid(jf) = nadrnew_grid_aux(ig_grid_nbraf(jf))
 120     CONTINUE
C     
C* Print memory required for field and grid arrays
C
         il_memtot(1) = ((ig_maxold + 3*ig_maxold_grid)*ip_realwp_p
     $   + ig_maxold_grid*ip_intwp_p)/1000
         il_memtot(2) = ((ig_maxnew + 3*ig_maxnew_grid)*ip_realwp_p
     $   + ig_maxnew_grid*ip_intwp_p)/1000
C
         IF (nlogprt .GE. 1) THEN
             WRITE (UNIT = nulou,FMT = 1001) il_memtot(1), il_memtot(2)
         ENDIF
 1001    FORMAT(1H ,5X,
     $   'Memory (kB) requested for source field and grid arrays',
     $    I10,/, 1H ,5X,
     $   'Memory (kB) requested for target field and grid arrays', I10)
C     
C*    2. End of routine
C     --------------
C     
         DEALLOCATE (nsizold_grid_aux)
         DEALLOCATE (nsiznew_grid_aux)
         DEALLOCATE (nadrold_grid_aux)
         DEALLOCATE (nadrnew_grid_aux)
C
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou,FMT = *) ' '
            WRITE (UNIT = nulou,FMT = *) 
     $           '          --------- End of routine inidya ---------'
            CALL FLUSH (nulou)
         ENDIF
      ENDIF
      RETURN
      END


