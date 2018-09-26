       SUBROUTINE CLIM_Init_Oasis(cdexp,cdnam,knmods,
     &                     cdmodnam, kbcplproc, kbtotproc,cdmpiarg, 
     &                     ida_unitmod, ktrout,kinfo)
c
c*    *** Init ***   CLIM 3.0
c
c     purpose:
c     --------
c        start models
c
c     interface:
c     ----------
c        cdexp  : experience descriptor
c        cdnam  : name of the calling model
c	 knmods	: number of models to be launched
c        cdmodnam : names of models to be launched 
c        kbcplproc : number of processes implied in the coupling
c                    for each model to be launched
c        kbtotproc : total number of procs for each model to be launched
c        cdmpiarg: nothing significant for MPI1 and launching argument 
c                  for each model to be launched for MPI2
c        ida_unitmod: minimum logfile unit number for each model
c        ktrout : unit number for trace file
c	 kinfo	: output status
c
c     lib mpi:
c     -------
c        mpi-1 (or mpi-2)
c
c     author:
c     -------
c       Sophie Valcke  - CERFACS (08/09/00 -created from CLIM_Init)
c       Jean Latour - F.S.E. - Version MPMD start with mpi-1
c                   - implies the use of mpiexec server process on VPPs
c                     or the command mpirun on most platfomrs
c       Arnaud Caubel - FECIT (08/02) - removed some arguments and added 
c                                       dynamic allocation
c
c     ----------------------------------------------------------------
c
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
      USE mod_kinds_oasis
      USE mod_clim
      USE mod_comclim
      USE mod_parameter
      USE mod_timestep
      USE mod_string
      USE mod_unitncdf
      USE mod_hardware
      USE mod_calendar
cvg>>
      USE mod_label
cvg<<
      IMPLICIT NONE
#include <mpif.h>
c     ----------------------------------------------------------------
      CHARACTER*(*) cdexp, cdnam
      INTEGER (kind=ip_intwp_p)	    knmods, ktrout,kinfo
      CHARACTER*6   cdmodnam(knmods)
      INTEGER (kind=ip_intwp_p)     kbtotproc(knmods), 
     $    kbcplproc(knmods), ida_unitmod(knmods)
      CHARACTER*80  cdmpiarg(knmods)
c     ----------------------------------------------------------------
      INTEGER (kind=ip_intwp_p)     iost, ip, iprcou, imodst, iprmod
      INTEGER (kind=ip_intwp_p)     imyrank, inumproc, mynummod
      INTEGER (kind=ip_intwp_p)     il_maxcplproc, iargs, il_mpisize, 
     $    kcomm_local
      INTEGER (kind=ip_intwp_p)      il_CLIM_Maxport, il_CLIM_MaxLink
      INTEGER (kind=ip_intwp_p)     ji, jj, jl, jn, ibuff, il_err, ib
      INTEGER (kind=ip_intwp_p)     icolor, icolorl, ikey, il_size, 
     $    iposbuf, info
      INTEGER (kind=ip_intwp_p)       impi_intercomp(knmods), 
     $    impi_newcomm(knmods+1)
      INTEGER (kind=ip_intwp_p)       iarrb(2)
      INTEGER (kind=ip_intwp_p)     mpi_status(MPI_STATUS_SIZE)
      INTEGER (kind=ip_intwp_p)     icolmods(knmods)
      INTEGER (kind=ip_intwp_p), PARAMETER :: iargmax=20
      INTEGER (kind=ip_intwp_p), PARAMETER :: itagcol=9876
      REAL (kind=ip_realwp_p), DIMENSION(:), ALLOCATABLE ::  rl_work
      REAL (kind=ip_realwp_p) rl_testvar
      CHARACTER*32  cdmodnam_aux(knmods)
      CHARACTER*24  cdspawnarg(iargmax)
      CHARACTER*(CLIM_Clength),DIMENSION(:), ALLOCATABLE :: cunames
      INTEGER(kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: iarr_err
      INTEGER(kind=ip_intwp_p) :: integer_byte_size, ii, io_size, 
     $    integer_io_size, il_ibyt, il_int, il_char, il_log, il_rl_work
      LOGICAL :: ll_log
      CHARACTER*1 :: cl_char 
cvg>>
      INTEGER (kind=ip_intwp_p)  grids_start, grids_done
      INTEGER (kind=ip_intwp_p)  itype, ilen, itag 
      INTEGER (kind=ip_intwp_p)  idproc
      CHARACTER*8                cgrdnamnc 
      LOGICAL                    existent
cvg<<
c     ----------------------------------------------------------------
c
c*    0. get arguments and some Initilizations
c     ----------------------------------------
c
      kinfo  = CLIM_Ok
      cgroup=' '
      cgroup=cdexp
      cmynam=' '
      cmynam=cdnam
      impi_intercomp(:)=0
      impi_newcomm(:)=0
      iarrb(:)=0
      icolmods(:)=0
      cdmodnam_aux(:)=' '
      cdspawnarg(:)=' '
      mpi_status(:)=0
      rl_testvar=0.0_ip_realwp_p
      ibuff=0_ip_intwp_p
c
      nports = 0
      nlinks = 0
      nbsend = 0
      nbrecv = 0
c
c*    1. open trace file 
c     ------------------
c
      iost = 0
      nulprt = ktrout
      WRITE (cnaprt,FMT='(A,''.prt'')') cdnam
      OPEN (UNIT=nulprt, FILE=cnaprt, STATUS='UNKNOWN',
     *      FORM='FORMATTED', ERR=110, IOSTAT=iost)
c
  110 CONTINUE
      IF (iost.ne.0) THEN
         WRITE(6,*) 'ABORT in Init_Oasis - unable to open trace file ', 
     $        iost
         WRITE(6,*) nulprt, ' ', cnaprt
         CALL FLUSH(nulprt)
         CALL ABORT
      ENDIF
c
      WRITE (nulprt, FMT='(A)')  'Init_Oasis - CLIM/MPI'
      WRITE (nulprt, FMT='(A,A)') 'Init_Oasis - experiment name: ', 
     $    cgroup
      WRITE (nulprt, FMT='(A,A)') 'Init_Oasis - name of the model: ',
     $    cdnam
c
c*    2. Start global MPI environment for OASIS and models
c     ----------------------------------------------------
c
      CALL MPI_INIT(mpi_err)
      IF(mpi_err .NE. MPI_SUCCESS) GOTO 211
      WRITE(UNIT = nulprt,FMT = *)' MPI_INIT done '
c
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,inumproc,mpi_err)
      IF(mpi_err .NE. MPI_SUCCESS) GOTO 212
      WRITE(UNIT = nulprt,FMT = *)' MPI_COMM_SIZE done ',inumproc
c
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,imyrank,mpi_err)
      IF(mpi_err .NE. MPI_SUCCESS) GOTO 213
      WRITE(UNIT = nulprt,FMT = *)' MPI_COMM_RANK done ',imyrank
c
c
c*    3. Test which MPI version we want to use
c*    ----------------------------------------
c
#ifdef use_comm_MPI1
c
c*    MPI_COMM_WORLD is the global communicator for all processes
c*    it includes Oasis and all model processes
c
         CALL MPI_COMM_DUP(MPI_COMM_WORLD,mpi_comm,mpi_err)
         WRITE(nulprt,*)'Init_Oasis - - comm_dup done= ',mpi_comm
c
c*    3.1 If MPI1 context then split global communicator into disjoint 
c*        communicators local for each model
c     ----------------------------------------------------------------    
c
c*    3.1.1 generates a "color" from the model name (4 characters)
c*    ----------------------------------------------------------
c
         ALLOCATE (cunames(inumproc), stat = il_err)
         IF (il_err.NE.0) WRITE (nulprt,*) 'Init_Oasis - Error in 
     $        "cunames" allocation '
         cunames(:)=' '
         call MPI_Allgather(cmynam,CLIM_Clength,MPI_CHARACTER,
     &        cunames,CLIM_Clength,MPI_CHARACTER,
     &        mpi_comm,mpi_err)
         icolor=1
         DO WHILE ((trim(cmynam).ne.trim(cunames(icolor))).and.
     &        (icolor.le.inumproc))
            icolor=icolor+1
         END DO 
         IF (icolor.le.inumproc) THEN 
            WRITE(nulprt,*)'Init_OASIS - - table lookup:',
     $           cunames(icolor)
            icolor=icolor*100
            WRITE(nulprt,*)'Init_OASIS - - color after table lookup:',
     &           icolor
         ELSE 
            Write(nulprt,*) 'Init_OASIS - - Could not find myself in '
            Write(nulprt,*) '               the model namespace!'
            Write(nulprt,*) 'Init_OASIS - - Check namcouple and model'
            Write(nulprt,*) '               name tags for consistency!'
            call MPI_ABORT (mpi_comm, 0, mpi_err)
         ENDIF

c     
c*       3.1.2 computes color of all models
c        --------------------------------
c     
         DO ji=1,knmods
            cdmodnam_aux(ji) = cdmodnam(ji)
         END DO
         DO jl = 1,knmods
            icolorl=1
            do while ((trim(cdmodnam_aux(jl)).ne.
     $           trim(cunames(icolorl))).and.(icolorl.le.inumproc))
               icolorl=icolorl+1
            END DO
            IF (icolorl.le.inumproc) THEN 
               icolmods(jl)=100*icolorl
               WRITE(nulprt,*) 'MPI color of model ',jl,' ',
     *              cdmodnam_aux(jl),' is : ',icolmods(jl) 
            ELSE 
               Write(nulprt,*) 'Init_OASIS - - No match in the'
               Write(nulprt,*) '               model namespace!'
               Write(nulprt,*) 'Init_OASIS - - Check namcouple '
               Write(nulprt,*) 'model name tags for consistency!'
               call flush(nulprt)
               call MPI_ABORT (mpi_comm, 0, mpi_err)
            ENDIF 
         END DO
         DEALLOCATE (cunames)
c     
c*       3.1.3 split MPI_COMM_WORLD in local, disjoints, communicators
c*       -----------------------------------------------------------
c     
         ikey = 1
         call MPI_COMM_SPLIT(MPI_COMM_WORLD, icolor, ikey,
     *        kcomm_local, mpi_err)
         IF(mpi_err .NE. MPI_SUCCESS) GOTO 216
c     
c*       3.1.4 Send to each process (involved or not in coupling) its
c*             model number "mynummod". This is obtained by comparing the
c*             "color" of this process with the list of colors icolmods
c        ----------------------------------------------------------------

         DO jn = 1,inumproc-1
            CALL MPI_Recv(ibuff,1,MPI_INTEGER,jn,itagcol,mpi_comm,
     *           mpi_status, mpi_err)
            mynummod = -1
            DO jl = 1,knmods
               IF ( ibuff .EQ. icolmods(jl) ) mynummod = jl 
            ENDDO
c     
            IF ( mynummod .EQ. -1 ) GOTO 215
c     
            CALL MPI_Send(mynummod,1,MPI_INTEGER,jn,itagcol,mpi_comm,
     *           mpi_err)
         ENDDO
         CALL MPI_COMM_SIZE(mpi_comm,il_mpisize,mpi_err)
#else
c
c*    3.2 If MPI2 context then spawn models
c*    --------------------------------------
c 
         CALL MPI_COMM_DUP(MPI_COMM_WORLD,impi_newcomm(1),mpi_err)
         WRITE(nulprt,*)'Init - - comm_dup done= ',impi_newcomm(1)

         DO jl=1,knmods
c
c       Parse optional arguments for the model
c
c*      print optional arguments cdmpiarg(jl)
c
            WRITE(nulprt,*)'Init_Oasis - - model number : ',jl
            WRITE(nulprt,*)'Init_Oasis - - model name = ',cdmodnam(jl)
            WRITE(nulprt,*)'Init_Oasis - - argument string = ',
     $           cdmpiarg(jl)
c
            call CLIM_Parse(cdmpiarg(jl), cdspawnarg, iargs, iargmax)
c
            WRITE(nulprt,*)'Init_Oasis - - number of args = ',iargs
c     
c       Individual spawn of each model
c 
            ALLOCATE (iarr_err(kbtotproc(jl)), stat = il_err)
            IF (il_err.NE.0) THEN
                WRITE (nulprt,*) 'Error iarr_err allocation'
                CALL HALTE ('STOP in CLIM_Init_Oasis')
            ENDIF

            CALL MPI_COMM_SPAWN(cdmodnam(jl),cdspawnarg,
     &           kbtotproc(jl), MPI_INFO_NULL, 0,
     &           impi_newcomm(jl), impi_intercomp(jl), 
     &           iarr_err, mpi_err)
            DO ib = 1, kbtotproc(jl) 
              IF(iarr_err(ib) .NE. MPI_SUCCESS) GOTO 214
            END DO
            DEALLOCATE(iarr_err)
            IF(mpi_err .NE. MPI_SUCCESS) GOTO 214
            WRITE(UNIT = nulprt,FMT = *)' MPI_COMM_SPAWN done '
c     
            CALL MPI_Intercomm_merge(impi_intercomp(jl), .false., 
     &           impi_newcomm(jl+1), il_err)
            CALL MPI_Comm_size(impi_newcomm(jl+1), mpi_size, mpi_err)
            CALL MPI_Comm_rank(impi_newcomm(jl+1), mpi_rank, mpi_err)
c     
            WRITE(nulprt,*)'Init_Oasis - - oasis new size = ',mpi_size
            WRITE(nulprt,*)'Init_Oasis - - oasis new rank = ',mpi_rank
            iarrb(1) = jl+1
            iarrb(2) = knmods
            CALL MPI_BCAST(iarrb, 2, MPI_INTEGER, 0,
     &           impi_newcomm(jl+1), mpi_err)
c
         ENDDO
c
         WRITE(nulprt,*)'Init_Oasis - - oasis all models started '
         CALL MPI_COMM_DUP(impi_newcomm(knmods+1), mpi_comm, mpi_err)
         DO jl=1,knmods+1
            call MPI_COMM_FREE(impi_newcomm(jl), mpi_err)
         ENDDO
         CALL MPI_Comm_size(mpi_comm, il_mpisize, mpi_err)
c

#endif
         
C*    4.0 Initialize global variables, send information to 
c         model processes and allocate arrays.
c      --------------------------------------------------------

      ig_clim_nfield = ig_total_nfield
      integer_byte_size = BIT_SIZE(ii)/8
      INQUIRE (iolength=io_size) ii
      integer_io_size = io_size
      il_int = io_size/integer_io_size*integer_byte_size
      INQUIRE (iolength=io_size) rl_testvar
      il_ibyt = io_size/integer_io_size*integer_byte_size
      INQUIRE (iolength=io_size) cl_char
      il_char = io_size/integer_io_size*integer_byte_size
      INQUIRE (iolength=io_size) ll_log
      il_log = io_size/integer_io_size*integer_byte_size
      il_rl_work = (29 + 64*ig_clim_nfield) * (il_char/il_ibyt + 1) + 
     $    (9 + 3*knmods + 9*ig_clim_nfield) * (il_int/il_ibyt + 1) +
     $    3 * (il_log/il_ibyt + 1)
      ALLOCATE (rl_work(il_rl_work), stat=il_err)
      IF (il_ERR.ne.0) WRITE(nulprt,*)'Error in "rl_work" 
     $    allocation in CLIM_Init_Oasis routine ! '
      rl_work(:)=0
      il_size = il_rl_work * il_ibyt
      il_CLIM_Maxport = 2 * ig_clim_nfield
      ig_frqmin = ig_total_frqmin
      ig_ntime = ntime
c
      DO jn = 1, il_mpisize-1
        CALL MPI_Send(knmods, 1, MPI_INTEGER, jn, itagcol,
     $      mpi_comm, mpi_err)
        IF (mpi_ERR.ne.MPI_SUCCESS) THEN
            WRITE(UNIT = nulprt,FMT = *)'Init_Oasis - Problem with 
     $          sending of information to models !'
            CALL HALTE ('STOP in CLIM_Init_Oasis')
        ENDIF
        CALL MPI_Send(ig_clim_nfield, 1, MPI_INTEGER, jn, itagcol+1,
     $      mpi_comm, mpi_err)
        IF (mpi_ERR.ne.MPI_SUCCESS) THEN
            WRITE(UNIT = nulprt,FMT = *)'Init_Oasis - Problem with 
     $          sending of information to models !'
            CALL HALTE ('STOP in CLIM_Init_Oasis')
        ENDIF
      END DO
      iposbuf = 0
      call MPI_Pack (cgroup, 8, MPI_CHARACTER, rl_work, il_size,
     $     iposbuf, mpi_comm, info)
      call MPI_Pack (ntime, 1, MPI_INTEGER, rl_work, il_size,
     $     iposbuf, mpi_comm, info)
      call MPI_Pack (kbcplproc, knmods, MPI_INTEGER, rl_work, il_size,
     $     iposbuf, mpi_comm, info)
      call MPI_Pack (kbtotproc, knmods, MPI_INTEGER, rl_work, il_size,
     $     iposbuf, mpi_comm, info)
      call MPI_Pack (ida_unitmod, knmods, MPI_INTEGER, rl_work, il_size,
     $     iposbuf, mpi_comm, info)
      call MPI_Pack (ig_frqmin, 1, MPI_INTEGER, rl_work, il_size,
     $     iposbuf, mpi_comm, info)
      call MPI_Pack (cg_input_field, 8*ig_clim_nfield, MPI_CHARACTER, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call MPI_Pack (cg_output_field, 8*ig_clim_nfield, MPI_CHARACTER, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call MPI_Pack (ig_lag, ig_clim_nfield, MPI_INTEGER, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call MPI_Pack (ig_reverse, ig_clim_nfield, MPI_INTEGER, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call MPI_Pack (ig_invert, ig_clim_nfield, MPI_INTEGER, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call MPI_Pack (ig_freq, ig_clim_nfield, MPI_INTEGER, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call MPI_Pack (ig_total_nseqn, ig_clim_nfield, MPI_INTEGER, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call MPI_Pack (cg_restart_file, 8*ig_clim_nfield, MPI_CHARACTER, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call MPI_Pack (ig_no_rstfile, ig_clim_nfield, MPI_INTEGER, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call MPI_Pack (ig_nbr_rstfile, 1, MPI_INTEGER, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call MPI_Pack (lncdfrst, 1, MPI_LOGICAL, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call MPI_Pack (lg_oasis_field, 1, MPI_LOGICAL, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call MPI_Pack (ig_total_state, ig_clim_nfield, MPI_INTEGER, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call MPI_Pack (ig_local_trans, ig_clim_nfield, MPI_INTEGER, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call  MPI_Pack (cg_input_file, 32*ig_clim_nfield, MPI_CHARACTER,
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call  MPI_Pack (ig_date, 6, MPI_INTEGER,
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call  MPI_Pack (ig_numlab, ig_clim_nfield, MPI_INTEGER,
     $     rl_work, il_size, iposbuf, mpi_comm, info)
CSV>>
      call  MPI_Pack (cgrdnam, 5, MPI_CHARACTER,
     $     rl_work, il_size, iposbuf, mpi_comm, info) 
      call  MPI_Pack (cglonsuf, 4, MPI_CHARACTER,
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call  MPI_Pack (cglatsuf, 4, MPI_CHARACTER,
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call  MPI_Pack (crnlonsuf, 4, MPI_CHARACTER,
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call  MPI_Pack (crnlatsuf, 4, MPI_CHARACTER,
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call  MPI_Pack (cga_locatorbf, 4*ig_clim_nfield, MPI_CHARACTER,
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call  MPI_Pack (cga_locatoraf, 4*ig_clim_nfield, MPI_CHARACTER,
     $     rl_work, il_size, iposbuf, mpi_comm, info)
      call MPI_Pack (lg_bsend, 1, MPI_LOGICAL, 
     $     rl_work, il_size, iposbuf, mpi_comm, info)
CSV<< 
      DO jn = 1, il_mpisize-1
        CALL MPI_Send(rl_work, iposbuf, MPI_PACKED, jn, itagcol+2,
     $      mpi_comm, mpi_err)
        IF (mpi_ERR.ne.MPI_SUCCESS) THEN
            WRITE(UNIT = nulprt,FMT = *)'Init_Oasis - Problem with 
     $          sending of information to models !'
            CALL HALTE ('STOP in CLIM_Init_Oasis')
        ENDIF
      END DO
      il_maxcplproc = 1
      DO ji = 1, knmods
        IF (il_maxcplproc.lt.kbcplproc(ji)) 
     $      il_maxcplproc = kbcplproc(ji)
      END DO

      DEALLOCATE(rl_work)
c      
C*    4.1 Define ncplprocs, the total number of processes involved
C*    in the coupling, including Oasis and initialize il_CLIM_MaxLink
C     ------------------------------------------------------------
c
      ncplprocs=1
      DO 3 ji = 1, knmods
        ncplprocs = ncplprocs + kbcplproc(ji)
 3    CONTINUE
c
      il_CLIM_MaxLink = ncplprocs * il_CLIM_Maxport

C     4.2 Allocate and initialize arrays defined in comclim MODULE 
c     ------------------------------------------------------------
c
      ALLOCATE (modtid(0:ncplprocs-1), stat=il_err)
      IF (il_ERR.ne.0) WRITE(nulprt,*)'Error in "modtid" allocation in 
     $    CLIM_Init_Oasis routine ! '
      modtid(:)=0
       ALLOCATE (ncode(0:ncplprocs-1), stat=il_err)
      IF (il_ERR.ne.0) WRITE(nulprt,*)'Error in "ncode" allocation in 
     $    CLIM_Init_Oasis routine ! '
      ncode(:)=0
       ALLOCATE (cnames(0:ncplprocs-1), stat=il_err)
      IF (il_ERR.ne.0) WRITE(nulprt,*)'Error in "cnames" allocation in 
     $    CLIM_Init_Oasis routine ! '
      cnames(:)=' '
c
      ALLOCATE (myport(5+il_maxcplproc,il_CLIM_Maxport), 
     $    stat = il_ERR)
      IF (il_ERR.ne.CLIM_Ok) WRITE(UNIT = nulprt,FMT = *)' Problem 
     $    in myport allocation in CLIM_Init_Oasis !'
      myport(:,:)=0
      ALLOCATE (mydist(CLIM_ParSize, il_CLIM_Maxport), 
     $    stat = il_ERR)
      IF (il_ERR.ne.CLIM_Ok) WRITE(UNIT = nulprt,FMT = *)' Problem 
     $    in mydist allocation in CLIM_Init_Oasis !'
      mydist(:,:)=0
      ALLOCATE (cports(il_CLIM_Maxport), stat = il_ERR)
      IF (il_ERR.ne.CLIM_Ok) WRITE(UNIT = nulprt,FMT = *)' Problem 
     $    in mydist allocation in CLIM_Init_Oasis !'
      cports(:)=' '
      ALLOCATE(clrport(il_CLIM_Maxport), stat = il_ERR)
      IF (il_ERR.ne.CLIM_Ok) WRITE(nulprt,*)'Error in clrport 
     $    allocation in CLIM_Init_Oasis !'
      clrport(:)=' '
      ALLOCATE(irdist(CLIM_ParSize, il_CLIM_Maxport), stat = il_ERR)
      IF (il_ERR.ne.0) WRITE(nulprt,*)'Error in irdist allocation 
     $    in CLIM_Init_Oasis'
      irdist(:,:)=0
      ALLOCATE(irport(5,il_CLIM_Maxport ), stat = il_ERR)
      IF (il_ERR.ne.0) WRITE(nulprt,*)'Error in irport allocation 
     $    in CLIM_Init_Oasis'
      irport(:,:)=0
      ALLOCATE (mylink(4+CLIM_ParSize,il_CLIM_MaxLink), stat = il_ERR)
      IF (il_ERR.ne.CLIM_Ok) WRITE(UNIT = nulprt,FMT = *)' Problem 
     $    in mylink allocation in CLIM_Init_Oasis !'
      mylink(:,:)=0
      DO 10 ip = 0, ncplprocs-1
        modtid(ip) = -1
        cnames(ip) = ' '
 10   CONTINUE

c
C*    4.3 Define modtid, the vector giving, for each process involved
C*        in the coupling,  its number in mpi_comm (i.e in all model 
C*        processes involved OR NOT in the coupling) 
C         -----------------------------------------------------------
c
c     For coupler
      iprcou = 0
      imodst = 0
      modtid(0) = 0
c     For models
      DO 5 ji = 1, knmods
        IF (ji .eq. 1) THEN
           imodst = 1
        ELSE
           imodst = imodst + kbtotproc(ji-1)
        ENDIF
        iprmod = 0
        DO 7 jj = 1, kbcplproc(ji)
          iprcou = iprcou + 1
          iprmod = iprmod + 1
          modtid(iprcou) = imodst + iprmod - 1
 7      CONTINUE
 5    CONTINUE

C*    4.4 Define mynum, the number of the actual processor 
c         in all processors involved in the coupling 
C         -------------------------------------------------
c
c*    For the coupler (if mono-processor)
      mynummod = 0
      mynum = 0
c
      WRITE(nulprt,*)'Init_Oasis - - Oasis mynummod = ', mynummod
      WRITE(nulprt,*)'Init_Oasis - - Oasis mynum = ', mynum
      WRITE(nulprt,*)'Init_Oasis - - modtid() = ', modtid

C*    4.5 Initialize writing of the grids file by the models
C         --------------------------------------------------
c
      grids_start=0
      WRITE(nulprt,*)'Init_Oasis - - check grids file:'
      INQUIRE(FILE = cgrdnam, EXIST = existent)
      IF (.NOT. existent) THEN
         cgrdnamnc=cgrdnam//'.nc'
         INQUIRE(FILE = cgrdnamnc, EXIST = existent)
         IF (.NOT. existent) THEN
            grids_start=1
            WRITE(nulprt,*)
     $           '             - no grids file found'
            lncdfgrd = .true.
         ELSE
            WRITE(nulprt,*)
     $           'Init_Oasis - - grids file "', cgrdnamnc, '" found'
            lncdfgrd = .true.
         ENDIF
      ELSE
         WRITE(nulprt,*)
     $        'Init_Oasis - - grids file "', cgrdnam, '" found'
            lncdfgrd = .false.
      ENDIF
      
      idproc = 1
      DO jl = 1, knmods
c
c*    Send starting flag (grids_start) to all model processors
c
         DO ip = idproc, idproc + kbtotproc(jl) - 1
            WRITE(nulprt,FMT='(A)')'Send - grids_start'
            ilen = 1
            itype = MPI_INTEGER
            itag = itagcol+3
            CALL MPI_Send(grids_start, ilen, itype, ip, itag,
     $           mpi_comm, MPI_ERR)
            IF (MPI_ERR .EQ. MPI_SUCCESS) THEN
               WRITE(nulprt,FMT='(A,I2,A,I6,A,I2,A,I3,A,I5,A,I1)')
     $              'Send - <dest:',ip,'> <comm:',mpi_comm,'> <len:',
     $              ilen,'> <type:',itype,'> <tag:',itag,'> ::  ',
     $              grids_start
               CALL FLUSH(nulprt)
            ELSE
               WRITE(nulprt,*)'Init_Oasis - - Problem sending
     $              flag to initialize grids writing !'
               WRITE(nulprt,*)'MPI_ERR: ', MPI_ERR
               CALL HALTE ('STOP in CLIM_Init_Oasis')
            ENDIF
         ENDDO

         IF (.NOT. existent) THEN
c
c*          Send grid file names and suffixes (mod_label) to master
c*          proc of each model
c
c*          grids file name
c 
            WRITE(nulprt,FMT='(A)')'Send - cgrdnam'
            ilen = 5
            itype = MPI_CHARACTER
            itag = itagcol+4
            CALL MPI_Send(cgrdnam, ilen, itype, idproc, itag, mpi_comm,
     $                    MPI_ERR)
            IF (MPI_ERR .EQ. MPI_SUCCESS) THEN
               WRITE(nulprt,FMT='(A,I2,A,I6,A,I2,A,I3,A,I5,A,A5)')
     $              'Send - <dest:',idproc,'> <comm:',mpi_comm,
     $              '> <len:', ilen,'> <type:',itype,'> <tag:',itag,
     $              '> ::  ', cgrdnam
               CALL FLUSH(nulprt)
            ELSE
               WRITE(nulprt,*) 
     $              'Init_Oasis - - Problem sending grids file name !'
               CALL HALTE ('STOP in CLIM_Init_Oasis')
            ENDIF
c     
c*          masks file name
c
            WRITE(nulprt,FMT='(A)')'Send - cmsknam'
            ilen = 5
            itype = MPI_CHARACTER
            itag = itagcol+5
            CALL MPI_Send(cmsknam, ilen, itype, idproc, itag, mpi_comm,
     $                    MPI_ERR)
            IF (MPI_ERR .EQ. MPI_SUCCESS) THEN
               WRITE(nulprt,FMT='(A,I2,A,I6,A,I2,A,I3,A,I5,A,A5)')
     $              'Send - <dest:',idproc,'> <comm:',mpi_comm,
     $              '> <len:', ilen,'> <type:',itype,'> <tag:',itag,
     $              '> ::  ', cmsknam
               CALL FLUSH(nulprt)
            ELSE
               WRITE(nulprt,*) 
     $              'Init_Oasis - - Problem sending masks file name !'
               CALL HALTE ('STOP in CLIM_Init_Oasis')
            ENDIF
c
c*          areas file name
c
            WRITE(nulprt,FMT='(A)')'Send - csurnam'
            ilen = 5
            itype = MPI_CHARACTER
            itag = itagcol+6
            CALL MPI_Send(csurnam, ilen, itype, idproc, itag, mpi_comm,
     $                    MPI_ERR)
            IF (MPI_ERR .EQ. MPI_SUCCESS) THEN
               WRITE(nulprt,FMT='(A,I2,A,I6,A,I2,A,I3,A,I5,A,A5)')
     $              'Send - <dest:',idproc,'> <comm:',mpi_comm,
     $              '> <len:', ilen,'> <type:',itype,'> <tag:',itag,
     $              '> ::  ', csurnam
               CALL FLUSH(nulprt)
            ELSE
               WRITE(nulprt,*) 
     $              'Init_Oasis - - Problem sending areas file name !'
               CALL HALTE ('STOP in CLIM_Init_Oasis')
            ENDIF
c
c*          longitudes suffix
c
            WRITE(nulprt,FMT='(A)')'Send - cglonsuf'
            ilen = 4
            itype = MPI_CHARACTER
            itag = itagcol+7
            CALL MPI_Send(cglonsuf, ilen, itype, idproc, itag, mpi_comm,
     $                    MPI_ERR)
            IF (MPI_ERR .EQ. MPI_SUCCESS) THEN
               WRITE(nulprt,FMT='(A,I2,A,I6,A,I2,A,I3,A,I5,A,A4)')
     $              'Send - <dest:',idproc,'> <comm:',mpi_comm,
     $              '> <len:', ilen,'> <type:',itype,'> <tag:',itag,
     $              '> ::  ', cglonsuf
               CALL FLUSH(nulprt)
            ELSE
               WRITE(nulprt,*)
     $          'Init_Oasis - - Problem sending suffix for longitudes !'
               CALL HALTE ('STOP in CLIM_Init_Oasis')
            ENDIF
c
c*          latitudes suffix
c
            WRITE(nulprt,FMT='(A)')'Send - cglatsuf'
            ilen = 4
            itype = MPI_CHARACTER
            itag = itagcol+8
            CALL MPI_Send(cglatsuf, ilen, itype, idproc, itag, mpi_comm,
     $                    MPI_ERR)
            IF (MPI_ERR .EQ. MPI_SUCCESS) THEN
               WRITE(nulprt,FMT='(A,I2,A,I6,A,I2,A,I3,A,I5,A,A4)')
     $              'Send - <dest:',idproc,'> <comm:',mpi_comm,
     $              '> <len:', ilen,'> <type:',itype,'> <tag:',
     $              itag,'> ::  ', cglatsuf
               CALL FLUSH(nulprt)
            ELSE
               WRITE(nulprt,*)
     $           'Init_Oasis - - Problem sending suffix for latitudes !'
               CALL HALTE ('STOP in CLIM_Init_Oasis')
            ENDIF
c
c*          corner longitudes suffix
c
            WRITE(nulprt,FMT='(A)')'Send - crnlonsuf'
            ilen = 4
            itype = MPI_CHARACTER
            itag = itagcol+9
            CALL MPI_Send(crnlonsuf, ilen, itype, idproc, itag,mpi_comm,
     $                    MPI_ERR)
            IF (MPI_ERR .EQ. MPI_SUCCESS) THEN
               WRITE(nulprt,FMT='(A,I2,A,I6,A,I2,A,I3,A,I5,A,A4)')
     $              'Send - <dest:',idproc,'> <comm:',mpi_comm,
     $              '> <len:', ilen,'> <type:',itype,'> <tag:',
     $              itag,'> ::  ', crnlonsuf
               CALL FLUSH(nulprt)
            ELSE
               WRITE(nulprt,*)'Init_Oasis - - '
               WRITE(nulprt,*)
     $              '    Problem sending suffix for corner longitudes !'
               CALL HALTE ('STOP in CLIM_Init_Oasis')
            ENDIF
c
c*          corner latitudes suffix
c
            WRITE(nulprt,FMT='(A)')'Send - crnlatsuf'
            ilen = 4
            itype = MPI_CHARACTER
            itag = itagcol+10
            CALL MPI_Send(crnlatsuf, ilen, itype, idproc, itag,mpi_comm,
     $                    MPI_ERR)
            IF (MPI_ERR .EQ. MPI_SUCCESS) THEN
               WRITE(nulprt,FMT='(A,I2,A,I6,A,I2,A,I3,A,I5,A,A4)')
     $              'Send - <dest:',idproc,'> <comm:',mpi_comm,
     $              '> <len:', ilen,'> <type:',itype,'> <tag:',
     $              itag,'> ::  ', crnlatsuf
               CALL FLUSH(nulprt)
            ELSE
               WRITE(nulprt,*)'Init_Oasis - - '
               WRITE(nulprt,*)
     $              '    Problem sending suffix for corner latitudes !'
               CALL HALTE ('STOP in CLIM_Init_Oasis')
            ENDIF
c
c*          masks suffix
c
            WRITE(nulprt,FMT='(A)')'Send - cmsksuf'
            ilen = 4
            itype = MPI_CHARACTER
            itag = itagcol+11
            CALL MPI_Send(cmsksuf, ilen, itype, idproc, itag, mpi_comm,
     $                    MPI_ERR)
            IF (MPI_ERR .EQ. MPI_SUCCESS) THEN
               WRITE(nulprt,FMT='(A,I2,A,I6,A,I2,A,I3,A,I5,A,A4)')
     $              'Send - <dest:',idproc,'> <comm:',mpi_comm,
     $              '> <len:', ilen,'> <type:',itype,'> <tag:',itag,
     $              '> ::  ', cmsksuf
               CALL FLUSH(nulprt)
            ELSE
               WRITE(nulprt,*)
     $              'Init_Oasis - - Problem sending suffix for masks !'
               CALL HALTE ('STOP in CLIM_Init_Oasis')
            ENDIF
c
c*          areas suffix
c     
            WRITE(nulprt,FMT='(A)')'Send - csursuf'
            ilen = 4
            itype = MPI_CHARACTER
            itag = itagcol+12
            CALL MPI_Send(csursuf, ilen, itype, idproc, itag, mpi_comm,
     $                    MPI_ERR)
            IF (MPI_ERR .EQ. MPI_SUCCESS) THEN
               WRITE(nulprt,FMT='(A,I2,A,I6,A,I2,A,I3,A,I5,A,A4)')
     $              'Send - <dest:',idproc,'> <comm:',
     $              mpi_comm,'> <len:', ilen,'> <type:',itype,
     $              '> <tag:',itag,'> ::  ', csursuf
               CALL FLUSH(nulprt)
            ELSE
               WRITE(nulprt,*)
     $              'Init_Oasis - - Problem sending suffix for areas !'
               CALL HALTE ('STOP in CLIM_Init_Oasis')
            ENDIF
c
c*          Wait until writing is finished
c
            grids_done=0
            WRITE(nulprt,FMT='(A)')'Recv - grids_done'
            ilen = 1
            itype = MPI_INTEGER
            itag = itagcol+13
            CALL MPI_Recv(grids_done, ilen, itype, idproc, itag,
     $           mpi_comm, mpi_status, MPI_ERR)
            IF (MPI_ERR .EQ. MPI_SUCCESS) THEN
               WRITE(nulprt,FMT='(A,I2,A,I6,A,I2,A,I3,A,I5,A,I1)')
     $              'Recv - <from:',idproc,'> <comm:',
     $              mpi_comm,'> <len:', ilen,'> <type:',itype,
     $              '> <tag:',itag,'> ::  ', grids_done
               CALL FLUSH(nulprt)
            ELSE
               WRITE(nulprt,*)'Init_Oasis - - Problem receiving
     $              flag that grids file is written!'
               WRITE(nulprt,*)'MPI_ERR: ', MPI_ERR
               CALL HALTE ('STOP in CLIM_Init_Oasis')
            ENDIF
         ENDIF
         idproc = idproc + kbtotproc(jl)
      ENDDO
c 
C*    5. Normal EXIT
C     --------------
C
      WRITE (nulprt,FMT='(A)') 'Init_Oasis - -'
      CALL FLUSH(nulprt)
      RETURN
C
C*    6. Error : stop
C     ---------------

 211  CONTINUE
      WRITE (UNIT = nulprt,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulprt,FMT = *)
     $    ' Problem when starting MPI_Init_Oasis !!! ',
     $    ' Mpi error code = ',mpi_err
      WRITE (UNIT = nulprt,FMT = *) ' '
      WRITE (UNIT = nulprt,FMT = *) ' '
      WRITE (UNIT = nulprt,FMT = *)
     $    ' We STOP!!! Check the job parameters '
      WRITE (UNIT = nulprt,FMT = *) ' '
      CALL HALTE ('STOP in CLIM_Init_Oasis')
 212  CONTINUE
      WRITE (UNIT = nulprt,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulprt,FMT = *)
     $    ' Problem with MPI_Comm_Size function !!! ',
     $    ' Mpi error code = ',mpi_err
      WRITE (UNIT = nulprt,FMT = *)
     $    ' Number of processors = ',inumproc
      WRITE (UNIT = nulprt,FMT = *) ' '
      WRITE (UNIT = nulprt,FMT = *) ' '
      WRITE (UNIT = nulprt,FMT = *)
     $    ' We STOP!!! Check the job parameters '
      WRITE (UNIT = nulprt,FMT = *) ' '
      CALL HALTE ('STOP in CLIM_Init_Oasis')
 213  CONTINUE
      WRITE (UNIT = nulprt,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulprt,FMT = *)
     $    ' Problem with MPI_Comm_Rank function !!! ',
     $    ' Mpi error code = ',mpi_err
      WRITE (UNIT = nulprt,FMT = *)
     $    ' Rank of processor is : ',imyrank
      WRITE (UNIT = nulprt,FMT = *) ' '
      WRITE (UNIT = nulprt,FMT = *) ' '
      WRITE (UNIT = nulprt,FMT = *)
     $    ' We STOP!!! Check the job parameters '
      WRITE (UNIT = nulprt,FMT = *) ' '
      CALL HALTE ('STOP in CLIM_Init_Oasis')
  214  CONTINUE
      WRITE (UNIT = nulprt,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulprt,FMT = *)
     $    ' The spawn of the models has a problem '
      WRITE (UNIT = nulprt,FMT = *) ' '
      WRITE (UNIT = nulprt,FMT = *)
     $    ' In spawn mpi_err ',mpi_err
      WRITE (UNIT = nulprt,FMT = *) ' '
      WRITE (UNIT = nulprt,FMT = *)
     $    ' We STOP!!! Check the job parameters '
      WRITE (UNIT = nulprt,FMT = *) ' '
      CALL HALTE ('STOP in CLIM_Init_Oasis')     
 215  CONTINUE
      WRITE (UNIT = nulprt,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulprt,FMT = *)
     $    ' The color of one process is incorrect '
      WRITE (UNIT = nulprt,FMT = *) ' '
      WRITE (UNIT = nulprt,FMT = *)
     $    ' proc number =  ',jn, ' color = ',ibuff
      WRITE (UNIT = nulprt,FMT = *) ' '
      WRITE (UNIT = nulprt,FMT = *)
     $    ' We STOP!!! Check the job parameters '
      WRITE (UNIT = nulprt,FMT = *) ' '
      CALL HALTE ('STOP in CLIM_Init_Oasis')                 
 216  CONTINUE
      WRITE (UNIT = nulprt,FMT = *) '        ***WARNING***'
      WRITE (UNIT = nulprt,FMT = *)
     $    ' Problem with MPI_Comm_Split function !!! ',
     $    ' Mpi error code = ',mpi_err
      WRITE (UNIT = nulprt,FMT = *) ' '
      WRITE (UNIT = nulprt,FMT = *) ' '
      CALL HALTE ('STOP in CLIM_Init_Oasis')      
c
#endif
      RETURN
      END
