      SUBROUTINE inicmc
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL C *
C               * -------------     ------- *
C               *****************************
C
C**** *inicmc*  - Initialize coupled mode communication for coupler
C
C     Purpose:
C     -------
C     Use either PIPE, CLIM, SIPC or GMEM library to start communication with
C     the models being coupled. The PIPE library uses named pipes(fifo)
C     while CLIM uses message passing library MPI-1 or MPI-2.
C     The SIPC library uses shared memory segments based on system V IPC.
C     The GMEM library uses the global memory concept of NEC machines.
C     In cases SIPC, GMEM or PIPE, signal handling is implemented to trap 
C     oasis or child status changes.
C
C**   Interface:
C     ---------
C       *CALL*  *inicmc*
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
C     (CLIM-PIPE-SIPC)_Init, (CLIM-PIPE-SIPC)_Define, (CLIM-PIPE-SIPC)_Stepi, 
C     CLIM_Start, chksgc
C
C     Reference:
C     ---------
C     See OASIS manual (2000)
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       1.0       L. Terray      94/01/01  created
C       2.0       L. Terray      95/09/01  modified : new structure
C       2.1       L. Terray      96/09/03  mofified : norm DOCTOR (loop
C                                          index jp --> jl)
C       2.2       S. Valcke      97/06/20  added: introduction of SVIPC
C       2.2       L. Terray      97/09/20  general cleaning + call chksgc
C       2.3       S. Valcke      99/04/30  added: printing levels
C       2.3       L. Terray      99/09/15  added: GMEM branch
C       2.4       S. Valcke      2K/05/01  added: llmodel for CLIM/MPI2
C       2.5       S. Valcke      2K/09/04  Remove cmach
C       2.5       S. Valcke      2K/09/08  -Remove llmodel
C                                          -Modified CALL to CLIM_Init
C       2.5       J. Latour      01/11/28  added : MPI1 use of mpirun
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------- Include files and USE of modules---------------------------
C
      USE mod_kinds_oasis
#if defined use_comm_MPI1 || defined use_comm_MPI2 || defined use_comm_GSIP || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_GSIP && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
      USE mod_clim
      USE mod_clim_def
      USE mod_comclim
#endif
      USE mod_parameter
      USE mod_parallel
      USE mod_string
      USE mod_unitncdf
      USE mod_experiment
      USE mod_timestep
      USE mod_unit
      USE mod_hardware
      USE mod_printing
      USE mod_analysis
      USE mod_label
      USE mod_gauss
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_GSIP && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
#include <mpif.h>
#endif
C
C* ---------------------------- Local declarations ----------------------
C
      INTEGER (kind=ip_intwp_p) il_dimid(2)
      INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: info 
      INTEGER (kind=ip_intwp_p) iparal(3)
      INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE :: il_sizold,
     $    il_siznew, il_maxaux
      CHARACTER*8 clpinp, clpout
      CHARACTER*8 clwork, clstrg
      CHARACTER(len=8), DIMENSION(:), ALLOCATABLE :: cl_auxaf, cl_auxbf
      INTEGER (kind=ip_intwp_p) :: il_ibyt, il_bufsendsize, 
     $    il_bufsendsizebyt, il_varid, ii, io_size, integer_io_size
      INTEGER (kind=ip_intwp_p) :: integer_byte_size
      REAL(kind=ip_double_p) rl_testvar
      INTEGER (kind=ip_intwp_p) :: il_mpierr
c
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Allocation and initializations 
C        ------------------------------
C
      IF (nlogprt .GE. 1) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '           ROUTINE inicmc  -  Level C'
          WRITE (UNIT = nulou,FMT = *) 
     $    '           **************     *******'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' Process stuff initialization'
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) ' '
      ENDIF
      iparal(:)=0
C
C* Set up signal handling
C
      CALL chksgc
C
C* Allocate and initialize error codes (only if all fields aren't exchanged 
C  directly)
C
      IF (lg_oasis_field) THEN
         ALLOCATE(info(ig_nfield))
         CALL izero (info, ig_nfield)
      ENDIF
   
C
c
C* PIPE initialization
c     
#ifdef use_comm_PIPE
C     
C* Open pipes between models to exchange checking data
C
      infos = 0
      CALL PIPE_Init (cjobnam, cmodnam, ig_nmodel, infos)
      IF (infos .NE. 0) THEN 
         WRITE (UNIT = nulou,FMT = *) 
     $        ' WARNING : Problem with model pipe initialization'
         WRITE (UNIT = nulou,FMT = *)
     $        ' =======   Error code number = ',infos
         CALL HALTE('STOP in inicmc')
      ENDIF
C     
C* SVIPC or GMEM initialization
C     
#elif defined use_comm_SIPC || defined use_comm_GMEM
C     
C* Open shared memory pools used by models to exchange initial infos 
C     
      infos = 0
      CALL SIPC_Init (cjobnam, cmodnam, ig_nmodel, infos)
      IF (infos .NE. 0) THEN 
         WRITE (UNIT = nulou,FMT = *)  
     $ 'WARNING : Problem with model shared memory pool initialization'
         WRITE (UNIT = nulou,FMT = *)
     $        ' =======   Error code number = ',infos
         CALL HALTE('STOP in inicmc')
      ENDIF
C
C* CLIM/MPI2 or MPI1 initialization
C
#elif defined use_comm_MPI1 || defined use_comm_MPI2 || defined use_comm_GSIP || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_GSIP && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
      ione = 1
      IF (cchan .EQ. 'MPI1' .or. cchan .EQ. 'GSIP') cmpiarg(:)=' '
      CALL CLIM_Init_Oasis(cjobnam, 'Oasis', 
     $     ig_nmodel, cmodnam, nbcplproc, nbtotproc, cmpiarg,
     $     iga_unitmod, nultr, infos)
      IF (infos .NE. CLIM_Ok) THEN 
         WRITE (UNIT = nulou,FMT = *) 
     $        ' WARNING : Problem with MPI initialization '
         WRITE (UNIT = nulou,FMT = *)
     $        ' =======  CLIM error code number = ',infos
         CALL HALTE('STOP in inicmc')
      ENDIF
      IF (.not. lg_oasis_field) THEN 
         WRITE (UNIT = nulou,FMT = *) 
     $        '===>>> All the fields are exchanged directly ' 
         WRITE (UNIT = nulou,FMT = *) 
     $        ' Oasis calls only his MPI initialization and 
     $        sleeps until his MPI finalization'
      ENDIF
#else
C     
C* No message passing case
C 
      IF (cchan .EQ. 'NONE') THEN 
         IF (nlogprt .GE. 1) THEN
            WRITE (UNIT = nulou, FMT = *) 
     $           ' NOTE : No message passing used in this run'
            WRITE (UNIT = nulou,FMT = *)
     $           ' ======= '
         ENDIF
      ELSE
         CALL prcout('WARNING: your $CHANNEL option is wrong:',
     $        cchan,1)
         CALL HALTE('Wrong $CHANNEL option: STOP in inicmc')
      ENDIF
C      
#endif
C
C* Put here stuff moved from inipar_alloc because grids.nc (and therefore
C  grid dimensions) is not available yet.
C*      If grids.nc file exists, get grid dimensions there.
C
        il_dimid(:)=0
        il_auxbf=0
        il_auxaf=0
        allocate(cl_auxbf(ig_nfield))
        cl_auxbf(:)=' '
        allocate(cl_auxaf(ig_nfield))
        cl_auxaf(:)=' '
        ig_maxold_grid = 0
        ig_maxnew_grid = 0
C
#ifdef use_netCDF

        IF (lncdfgrd) CALL hdlerr
     $      (NF_OPEN(cgrdnam//'.nc', NF_NOWRITE, nc_grdid), 'inicmc0')
#endif
C       
        DO 252 jf=1,ig_total_nfield
          IF (ig_total_state(jf) .ne. ip_input) THEN
#ifdef use_netCDF
              IF (lncdfgrd) THEN
C
C*                Initial dimensions
                  clwork = cga_locatorbf(jf)
                  icount = ilenstr(clwork,jpeight)
                  clstrg = clwork(1:icount)//cglonsuf
                  CALL hdlerr
     $                (NF_INQ_VARID(nc_grdid, clstrg, il_varid), 
     $                'inicmc1')
                  CALL hdlerr
     $                (NF_INQ_VARDIMID(nc_grdid,il_varid, il_dimid),
     $                'inicmc2')
                  CALL hdlerr
     $                (NF_INQ_DIMLEN(nc_grdid, il_dimid(1), il_len),
     $                'inicmc3')
                  IF (lg_state(jf))
     $                nlonbf(ig_number_field(jf)) = il_len
                  CALL hdlerr
     $                (NF_INQ_DIMLEN(nc_grdid, il_dimid(2), il_len),
     $                'inicmc4')
                  IF (lg_state(jf))
     $                nlatbf(ig_number_field(jf)) = il_len
C     
C*                Final dimensions
                  clwork = cga_locatoraf(jf)
                  icount = ilenstr(clwork,jpeight)
                  clstrg = clwork(1:icount)//cglonsuf
                  CALL hdlerr
     $                (NF_INQ_VARID(nc_grdid, clstrg, il_varid), 
     $                'inicmc5')
                  CALL hdlerr
     $                (NF_INQ_VARDIMID(nc_grdid,il_varid, il_dimid),
     $                'inicmc6')
                  CALL hdlerr
     $                (NF_INQ_DIMLEN(nc_grdid, il_dimid(1), il_len),
     $                'inicmc7')
                  IF (lg_state(jf))
     $                nlonaf(ig_number_field(jf)) = il_len
                  CALL hdlerr
     $                (NF_INQ_DIMLEN(nc_grdid, il_dimid(2), il_len),
     $                'inicmc8')
                  IF (lg_state(jf))
     $                nlataf(ig_number_field(jf)) = il_len
C
C
              ENDIF
#endif
           ENDIF
C
C* Get the number of different grids used 
C
           IF (lg_state(jf)) then
              IF (ig_number_field(jf).eq.1) THEN
                 cl_auxbf(1) = cficbf(1)
                 cl_auxaf(1) = cficaf(1)
                 ig_grid_nbrbf(1)=1
                 ig_grid_nbraf(1)=1
                 il_auxbf = 1
                 il_auxaf = 1
                 ig_maxold_grid = nlonbf(ig_number_field(jf)) * 
     $                nlatbf(ig_number_field(jf))
                 ig_maxnew_grid = nlonaf(ig_number_field(jf)) * 
     $                nlataf(ig_number_field(jf))
              ELSEIF (ig_number_field(jf).gt.1) THEN
                 IF (ALL(cl_auxbf.ne.
     $                cficbf(ig_number_field(jf)))) THEN
                    il_auxbf=il_auxbf + 1
                    cl_auxbf(il_auxbf)=cficbf(ig_number_field(jf))
                    ig_grid_nbrbf(ig_number_field(jf))=il_auxbf
                    ig_maxold_grid = ig_maxold_grid + 
     $                   (nlonbf(ig_number_field(jf)) *
     $                   nlatbf(ig_number_field(jf)))
                 ELSE
                    DO ib = 1, il_auxbf
                       IF (cficbf(ig_number_field(jf)).eq.
     $                      cl_auxbf(ib))
     $                      ig_grid_nbrbf(ig_number_field(jf))= ib
                    ENDDO
                 ENDIF
                 IF (ALL(cl_auxaf.ne.
     $                cficaf(ig_number_field(jf)))) THEN
                    il_auxaf=il_auxaf + 1
                    cl_auxaf(il_auxaf)=cficaf(ig_number_field(jf))
                    ig_grid_nbraf(ig_number_field(jf))=il_auxaf
                    ig_maxnew_grid = ig_maxnew_grid + 
     $                   (nlonaf(ig_number_field(jf)) *
     $                   nlataf(ig_number_field(jf)))
                 ELSE
                    DO ib = 1, il_auxaf
                       IF (cficaf(ig_number_field(jf)).eq.
     $                      cl_auxaf(ib))
     $                      ig_grid_nbraf(ig_number_field(jf))= ib
                    ENDDO
                 ENDIF  
             ENDIF
         ENDIF
 252   END DO
C
C*     Close netcdf FILE
#ifdef use_netCDF
       IF (lncdfgrd) call hdlerr (NF_CLOSE(nc_grdid),'inicmc9')
#endif
C 
       IF (lg_oasis_field) THEN 
          ALLOCATE (il_sizold(ig_nfield))
          ALLOCATE (il_siznew(ig_nfield))
          il_sizold(:) = 0
          il_siznew(:) = 0
          DO jf=1,ig_nfield
            il_sizold (jf) = nlonbf(jf)*nlatbf(jf)
            il_siznew (jf) = nlonaf(jf)*nlataf(jf)
            ig_maxgrd = imaxim(il_sizold, ig_nfield)
            IF (imaxim(il_siznew, ig_nfield).gt.ig_maxgrd)
     $          ig_maxgrd = imaxim(il_siznew, ig_nfield)
          ENDDO
          WRITE(nulou,*)
     $        'Maximum size of the different grids of indirect fields:',
     $        ig_maxgrd
          WRITE(nulou,*)' '
C
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_GSIP && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
          ig_CLIMmax = 3 + CLIM_Clength + 2*ig_total_nfield * 
     $    (CLIM_Clength + 5 + CLIM_ParSize)
          IF (ig_CLIMmax.lt.ig_maxgrd ) 
     $        ig_CLIMmax = ig_maxgrd
          ALLOCATE (pkwork(ig_CLIMmax), stat = il_ERR)
          IF (il_ERR.NE.0) WRITE(nulou,*)
     $        'Error in pkwork allocation in inicmc '
          pkwork(:)=0 
#elif defined use_comm_GSIP
          ig_CLIMmax = ig_maxgrd
          ALLOCATE (pkwork(ig_CLIMmax), stat = il_ERR)
          IF (il_ERR.NE.0) WRITE(nulou,*)
     $        'Error in pkwork allocation in inicmc '
          pkwork(:)=0 
          
#endif          
C     
C*          Search sum of grids sizes of fields before and after interpolation
C     
          ig_maxold = isumi (il_sizold, ig_nfield)
          ig_maxnew = isumi (il_siznew, ig_nfield)
          WRITE(nulou,*)
     $        'Sum of grid sizes of fields before interpolation : ',
     $        ig_maxold
          WRITE(nulou,*)' '
          WRITE(nulou,*)
     $        'Sum of grid sizes of fields after interpolation : ',
     $        ig_maxnew
          WRITE(nulou,*)' '
C     
C*          Search dimension of "nwork" array
C     
          ALLOCATE (il_maxaux(ig_nfield))
          il_maxaux(:)=0
          DO jf = 1,ig_nfield
            il_redu = 0
            IF (ntronca(jf) .NE. 0) THEN
                IF (ntronca(jf) .EQ. 16) THEN
                    il_redu = nredu16
                ELSE IF (ntronca(jf) .EQ. 24)  THEN 
                    il_redu = nredu24
                ELSE IF (ntronca(jf) .EQ. 32)  THEN 
                    il_redu = nredu32
                ELSE IF (ntronca(jf) .EQ. 48)  THEN 
                    il_redu = nredu48   
                ELSE IF (ntronca(jf) .EQ. 80)  THEN 
                    il_redu = nredu80
                ELSE IF (ntronca(jf) .EQ. 160)  THEN 
                    il_redu = nredu160
                ELSE 
                    CALL prtout
     $                  ('WARNING!!! Oasis cannot treat this grid with 
     $                  2*NO latitude lines with NO = ',
     $                  ntronca(jf), 2)
                ENDIF
            ENDIF
            IF (il_siznew(jf).ge.(il_redu+nlatbf(jf))) THEN
                il_maxaux(jf) = il_siznew(jf)
            ELSE
                il_maxaux(jf) = il_redu+nlatbf(jf)
            ENDIF
          END DO
          ig_nwork = imaxim (il_maxaux, ig_nfield)
C     
C*          Search dimension of "work" array
C     
          CALL izero(il_maxaux, ig_nfield)
          DO jf = 1, ig_nfield
            il_maxaux(jf) = nbnfld(jf)*il_siznew(jf)
            IF (il_maxaux(jf).lt.(nbofld(jf)*il_sizold(jf)))
     $          il_maxaux(jf) = nbofld(jf)*il_sizold(jf)
            IF (il_maxaux(jf).lt.(nlonaf(jf)+nlataf(jf)
     $          +2*il_sizold(jf)+2*il_siznew(jf))) 
     $          il_maxaux(jf) = nlonaf(jf) + nlataf(jf)
     $          +2*il_sizold(jf) + 2*il_siznew(jf)
          END DO
          ig_work = imaxim(il_maxaux, ig_nfield)
          IF (ig_work.lt.(3*ig_maxgrd)) 
     $        ig_work = 3*ig_maxgrd
C     
          DEALLOCATE (il_maxaux)
          DEALLOCATE (il_sizold)
          DEALLOCATE (il_siznew)
      ENDIF
      DEALLOCATE (cl_auxbf)
      DEALLOCATE (cl_auxaf)
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_GSIP && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
      IF (lg_bsend) THEN
          ii=0
          rl_testvar = 0.0_ip_double_p
          integer_byte_size = BIT_SIZE(ii)/8
          INQUIRE (iolength=io_size) ii
          integer_io_size = io_size
          INQUIRE (iolength=io_size) rl_testvar
          il_ibyt = io_size/integer_io_size*integer_byte_size
          il_bufsendsize = 2*ig_clim_nfield*
     $        (ig_maxgrd + MPI_BSEND_OVERHEAD/il_ibyt+1)
          ALLOCATE (dg_bufsend(il_bufsendsize), stat = il_err)
          IF (il_ERR.ne.0) WRITE(nulou,*)'Error in dg_bufsend
     $        allocation in CLIM_Init_Oasis'
          dg_bufsend(:)=0
          il_bufsendsizebyt = il_bufsendsize * il_ibyt 
C
C*   Attach a buffer able to contain the maximum possible number of
C    coupling fields with maximum possible
C    size (ig_maxgrd) supposed to be declared as DOUBLE.  
          CALL MPI_Buffer_Attach(dg_bufsend,il_bufsendsizebyt,il_mpierr)
          WRITE(nulou,*)'Attached buffer of size=', il_bufsendsizebyt
      ENDIF
#endif
C
C* END of stuff moved from inipar_alloc         
C
C* Define ports and start the communication only if one field (at least) 
C  goes through Oasis
C
      IF (lg_oasis_field) THEN
C     
C*    2. Define ports or files for data exchange
C     ---------------------------------------
C     
C     * PIPE case: define and open pipe files (FIFO)
C     
#ifdef use_comm_PIPE
         infos = 0
         DO 250 jf = 1, ig_nfield
C     
C     * Define pipes actual names for each field
C     
C     - Reading
C     
            clpinp = cnaminp(jf)
C     - Writing
            clpout = cnamout(jf)
C     
C     * Create pipes
C     
            CALL PIPE_Define (cnaminp(jf), cnamout(jf), 
     $           clpinp, clpout , info(jf))
            infos = infos + info(jf)
 250     CONTINUE
         IF (infos .NE. 0) THEN 
            WRITE (UNIT = nulou,FMT = *) 
     $           ' WARNING : Problem with field pipe initialization'
            WRITE (UNIT = nulou,FMT = *)
     $           ' =======   Error code number = ',infos
            CALL HALTE('STOP in inicmc')
         ENDIF
c         ENDIF
C     
C     * SIPC case
C     
#elif defined use_comm_SIPC || defined use_comm_GMEM
         infos = 0
         DO 260 jf = 1, ig_nfield
C     
C     * Create two shared memory pools for each field
C     
            CALL SIPC_Define (jf, info(jf))  
            infos = infos + info(jf)
 260     CONTINUE
         IF (infos .NE. 0) THEN 
            WRITE (UNIT = nulou,FMT = *) 
     $           'WARNING : Problem with field shared memory pool 
     $           initialization'
            WRITE (UNIT = nulou,FMT = *)
     $           ' =======   Error code number = ',infos
            CALL HALTE('STOP in inicmc')
         ELSE
C     
C     * If everything went alright, open one dummy file signaling  
C     that pools for exchange are opened
C     
            OPEN (UNIT = nudum, FILE = 'DUMMY_SIPC', STATUS 
     $           = 'NEW')
            CLOSE (UNIT = nudum)
            IF (nlogprt .GE. 1) THEN
               WRITE (UNIT = nulou,FMT = *)
     $              '* ===>>> : file DUMMY_SIPC created'
               WRITE (UNIT = nulou,FMT = *) 
     $              '  ------   -----------------------'
            ENDIF
        ENDIF
C     
C     * CLIM case: define ports
C
#elif defined use_comm_MPI1 || defined use_comm_MPI2 || defined use_comm_GSIP  || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_GSIP && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE) 
C     
C     * Input ports for incoming fields
C
         infos = CLIM_Ok
         iparal(CLIM_Strategy) = CLIM_Serial
         iparal(CLIM_Offset) = 0
         DO 210 jf = 1, ig_nfield
            iparal(CLIM_Length) = nlonbf(jf) * nlatbf(jf)
            CALL CLIM_Define (ig_portin_id(jf),cnaminp(jf),CLIM_In, 
     $           CLIM_Real,iparal, info(jf))
            infos = infos + info(jf)
 210     CONTINUE
         IF (infos .NE. CLIM_Ok) THEN
            DO 220 jf = 1, ig_nfield
               IF (info(jf) .NE. CLIM_Ok) THEN
                  WRITE (UNIT = nulou,FMT = *) 
     $                 ' WARNING : Problem with port ', cnaminp(jf)
                  WRITE (UNIT = nulou,FMT = *)
     $                 ' =======   Error code number = ',info(jf)
               ENDIF 
 220        CONTINUE 
            CALL HALTE ('STOP in inicmc')
         ENDIF 
C     
C     * Output ports for outgoing fields : serial decomposition
C     
         infos = CLIM_Ok
         iparal(CLIM_Strategy) = CLIM_Serial
         iparal(CLIM_Offset) = 0
         DO 230 jf = 1, ig_nfield
            iparal(CLIM_Length) = nlonaf(jf) * nlataf(jf)
            CALL CLIM_Define (ig_portout_id(jf), cnamout(jf), 
     $           CLIM_Out, CLIM_Real, iparal, info(jf))
            infos = infos + info(jf)
 230     CONTINUE 
         IF (infos .NE. CLIM_Ok) THEN
            DO 240 jf = 1, ig_nfield
               IF (info(jf) .NE. CLIM_Ok) THEN
                  WRITE (UNIT = nulou,FMT = *) 
     $                 ' WARNING : Problem with port ', cnamout(jf)
                  WRITE (UNIT = nulou,FMT = *)
     $                 ' =======   Error code number = ',info(jf)
               ENDIF 
 240        CONTINUE 
            CALL HALTE ('STOP in inicmc')
         ENDIF
#endif
C     
C*    3. Start the communication and get timestep information
C     ----------------------------------------------------
#ifdef use_comm_PIPE
         infos = 0
         DO 320 jm = 1, ig_nmodel
            CALL PIPE_Stepi (cmodnam(jm), jm, 
     $           istep, ifcpl, idt, infos)
            IF (infos .NE. 0) THEN
               WRITE (UNIT = nulou,FMT = *) 
     $              ' WARNING : Problem in getting step info 
     $              from model ', cmodnam(jm)
               WRITE (UNIT = nulou,FMT = *)
     $              ' =======   Error code number = ',infos
               CALL HALTE('STOP in inicmc')
            ELSE
               IF (nlogprt .GE. 1) THEN
                  WRITE (UNIT = nulou,FMT = *)
     $                 ' Got step information from model ', 
     $                 cmodnam(jm)
               ENDIF
            ENDIF
            mstep(jm) = istep
            mfcpl(jm) = ifcpl
            mdt(jm) = idt
 320     CONTINUE 
C     
C     * SIPC or GMEM Case
C     
#elif defined use_comm_SIPC || defined use_comm_GMEM
         infos = 0
         DO 330 jm = 1, ig_nmodel
            CALL SIPC_Stepi (cmodnam(jm), jm, 
     $           istep, ifcpl, idt, infos)
            IF (infos .NE. 0) THEN
               CALL prcout
     $              ('WARNING: Problem in getting step info from 
     $              model', cmodnam(jm), 1)
               CALL prtout('Error code number = ',infos,1)
               CALL HALTE('STOP in inicmc')
            ELSE
               CALL prcout
     $              ('Got step information from model ',
     $              cmodnam(jm), 1)
            ENDIF
            mstep(jm) = istep
            mfcpl(jm) = ifcpl
            mdt(jm) = idt
 330     CONTINUE 
c
C     * CLIM CASE
C
#elif defined use_comm_MPI1 || defined use_comm_MPI2 || defined use_comm_GSIP || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_GSIP && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
c 
         CALL CLIM_Start (infos)
         IF (infos .NE. CLIM_Ok) THEN
            WRITE (UNIT = nulou,FMT = *) 
     $           ' WARNING : Problem in starting CLIM '
            WRITE (UNIT = nulou,FMT = *)
     $           ' =======   Error code number = ',infos
         ENDIF
         DO 310 jm = 1, ig_nmodel
            CALL CLIM_Stepi (cmodnam(jm), infos)
            IF (infos .lt. nbcplproc(jm)) THEN
               WRITE (UNIT = nulou,FMT = *) 
     $              'PROBLEM: Got initial informations only from',
     $              infos
               WRITE (UNIT = nulou,FMT = *)
     $              'processes of model',  cmodnam(jm) 
               CALL halte ('STOP in inicmc.f')
            ELSE
               IF (nlogprt .GE. 1) THEN  
                  WRITE (UNIT = nulou,FMT = *)
     $                 ' Got step informations from model ', 
     $                 cmodnam(jm)
               ENDIF
            ENDIF
 310     CONTINUE
#endif
      ENDIF
C
C
C*    4. End of routine
C        --------------
C
      IF (lg_oasis_field) DEALLOCATE(info)
C
      IF (nlogprt .GE. 1) THEN
          WRITE (UNIT = nulou,FMT = *) ' '
          WRITE (UNIT = nulou,FMT = *) 
     $    '          --------- End of routine inicmc ---------'
          CALL FLUSH (nulou)
      ENDIF
      RETURN
      END
