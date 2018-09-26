      SUBROUTINE driver
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL C *
C               * -------------     ------- *
C               *****************************
C
C**** *driver* - Main OASIS routine
C
C
C     Purpose:
C     -------
C     Drive and control the simulation between GCMs and coupler.
C     Contain the time loop. A coupled simulation with OASIS 2.0
C     starts with the interpolation of the boundary conditions
C     from their original grid to the target grid in contrast with
C     the previous versions. Consequently, the GCMs pause initially
C     until the coupling variables have been interpolated.
C
C     N.B: Note that the time loop goes from 0 to niter-1 in contrast
C          with previous versions. The iteration 0 of oasis DOES NOT
C          increment the simulation time.
C
C**   Interface:
C     ---------
C       *CALL*  *driver*
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
C               iindex : index array for field identificators
C
C     Externals:
C     ---------
C                        Initialization
C                        **************
C     inipar, inilun, iniiof, inidya, initim, inicmc, chkpar, inigrd,
C
C                        Temporal loop
C                        ************* 
C     getfld, preproc, interp, cookart, postpro, givfld, reset, updtim
C
C                        Synchronization
C                        ***************
C     modsgc, waitpc
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
C       2.0beta   L. Terray      95/07/24  modified: new structure
C       2.0       L. Terray      96/02/01  modified: change in time loop
C       2.1       O. Marti, L.T  96/09/25  added: extra time step
C       2.2       S. Valcke, L.T 97/11/13  added: SIPC call to modsgc
C                                                 mode no message passing
C       2.3       L. Terray      99/09/15  added: GMEM branch
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* -----------------Include files and USE of modules---------------------------
C
      USE mod_kinds_oasis
      USE mod_string
      USE mod_analysis
      USE mod_memory
      USE mod_parameter
      USE mod_experiment
      USE mod_timestep
      USE mod_unit
      USE mod_hardware
C
C* ---------------------------- Local declarations ----------------------
C
      INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: iindex
      LOGICAL lltime, llseqn, llend
      INTEGER (kind=ip_intwp_p) :: il_flag
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Initialization
C        --------------
C
C* - Initialize main run parameters
C
      CALL inipar_alloc
C
C* - "iindex" allocation
C
      IF (lg_oasis_field) THEN 
         ALLOCATE(iindex(ig_nfield))
         iindex(:)=0
      ENDIF
C
C* - Allocate arrays (first round)
C
      il_flag=1
      CALL alloc(il_flag)
C
C* - Initialize other run parameters
C
      CALL inipar  
C
C* - Initialize logical unit numbers
C
      CALL inilun
cvg C
cvg C* - Open necessary files
cvg C
cvg       CALL iniiof

C
C* - Initialize time information
C
      CALL initim
C
C* - Initialize communication between processes
C
      CALL inicmc
C
C* - Allocate arrays (second round)
C
      il_flag=2
      CALL alloc(il_flag)
C
cvg>>>
C
C* - Set up dynamic allocation for all grid-related fields
C
      CALL inidya
C
C* - Open necessary files
C
      CALL iniiof
cvg<<<

C
C* - Check run parameters compatibility between GCM's and coupler
C
      CALL chkpar
C
C* - Initialize GCM's grids
C
      CALL inigrd
C
C* - We will go through the next time loop only if one field (at least) goes
C    through Oasis

      IF (lg_oasis_field) THEN
C
C*    2. Time loop
C        ---------
C
C* Loop on number of iterations
C* First iteration takes place at ndate (INIDATE in namcouple) and
C* last iteration at one timestep before the end of the simulation
C* (INIDATE + RUNTIME in namcouple), 
C
         DO 210 jt = 0, nitfn
C
C* Assign local variable for iteration number
C
            iter = jt
C
C* Get time counter. 
C
            icount = iter * nstep
C
C* Update calendar date
C
            CALL updtim (iter)
C
C* Loop on number of sequential models
C
            DO 220 jm = 1, nmseq
C
C* Loop on number of fields to find active fields for current iteration
C
               ifield = 0
               DO 230 jf = 1, ig_total_nfield
C
C* Treat the field only IF it has to go through Oasis
C
                  IF (lg_state(jf)) THEN
C     
C* Get conditional logical flags for doing analysis set
C
C
C* Treat the field only if time smaller than end of simulation
                     llend  = icount .LT. ntime
C
C* Treat the field only if iteration corresponds to one 
C* of its coupling timesteps
                     ifnow = nfexch(ig_number_field(jf))
                     lltime = mod(icount,ifnow) .EQ. 0
C
C* If sequential order, treat the field only if it is consumed
C* in present loop
                     llseqn = nseqn(ig_number_field(jf)) .EQ. jm
C
C* Conditional test to fill up iindex array
C
                     IF (llseqn .AND. lltime .AND. llend) 
     $                    THEN
                        ifield = ifield + 1
C
C* Fill up iindex array with active fields at iteration jt
C
                        iindex(ifield) = ig_number_field(jf)
                     ENDIF
                  ENDIF
 230           CONTINUE
C
C* There are ifield fields to be exchanged for iteration jt
C
               IF (ifield .GT. 0) THEN
C
C* Get fields
C
                  CALL getfld (iindex, ifield, iter)
C
C* Do preprocessing
C
                  CALL preproc (iindex, ifield)
C
C* Do the interpolation
C
                  CALL interp (iindex, ifield)
C
C* Do the nitty gritty stuff

C
                  CALL cookart (iindex, ifield)
C
C* Do postprocessing
C
                  CALL postpro (iindex, ifield)
C
C* If last iteration in PIPE or SIPC case, switch sigcld handler
                  IF (iter .EQ. nitfn) THEN
                     IF (cchan .EQ. 'PIPE' .OR. cchan .EQ. 'SIPC')
     $                    CALL modsgc
                  ENDIF
C
C* Give back fields
C     
                  CALL givfld (iindex, ifield, iter)
               ENDIF
C
C* End of loop over the sequential models
C
 220        CONTINUE
C
C* Reset macro arrays
C
        CALL reset
C
C* End of iterative loop
C
 210  CONTINUE
      ENDIF
C
C*    3. Wait until end of child processes
C        ---------------------------------
C
      CALL waitpc
C
C
C*    4. End of routine
C        --------------
C
C* "iindex" deallocation
C
      IF (lg_oasis_field) DEALLOCATE(iindex)
C
C* Deallocation of arrays allocated in "inialloc" routine      
C
      CALL dealloc
C
      RETURN
      END
