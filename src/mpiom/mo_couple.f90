!-----------------------------------------------------------------------
! BOP
!
! !MODULE:  mo_couple.F90
!
MODULE mo_couple
!
! !DESCRIPTION:
!
!   contains all additional code used for coupling with OASIS3
!
! !USES:
!
#if defined __coupled || __prism

#ifdef NAG
         USE f90_unix_io, ONLY: flush
#endif

  USE mo_param1,          ONLY: ie,je
  USE mo_units,           ONLY: io_in_octl, io_stdout, io_stderr
  USE mo_commo1,          ONLY: iaufr,lyear1,lmont1,dt,ldt,lmonts, &
                                v1e,z1e,u1e,b1e,fu10, &
                                weto, sicomo, sictho, sicsno, tho &
                                ,uko,vke,sicuo,sicve    &
                                ,weto_g,amsuo,amsue, ioasisflux  

#ifdef ADDCONTRA
  USE mo_param1_add,      ONLY: ih2o16,ih2o18,ihDo16
#endif
#ifdef __cpl_dust
  USE mo_carbch,          ONLY: dustdep
#endif
#ifdef __cpl_dms
  USE mo_carbch,          ONLY: dmsflux, dms
#endif
#ifdef __cpl_co2
  USE mo_carbch,          ONLY: co2conc, co2flux_cpl, co2trans, suppco2, atm
  USE mo_param1_bgc,      ONLY: iatmco2
#endif

  USE mo_commoau1,        ONLY: rhowat, tmelt, rhosno, con, consn

  USE mo_fluxes1,         ONLY: aofltxwo, aofltxio, aofltywe, aofltyie, &
                                aoflnhwo, aoflchio, aoflrhio, &
                                aoflfrwo, aoflfrio, aoflwsvo, &
                                aoflshwo
#ifdef ADDCONTRA
  USE mo_fluxes1,         ONLY: aoflfrwo16, aoflfrwo18, aoflfrwhdo, &
                                aoflfrio16, aoflfrio18, aoflfrihdo
  USE mo_contra,          ONLY: ocectra
#endif
  USE mo_thacc

  USE mo_mpi
  USE mo_parallel
  USE mo_kind

  USE mod_kinds_model,    ONLY: ip_realwp_p
  USE mod_prism_proto
  USE mod_comprism_proto
  USE mod_prism_get_proto
  USE mod_prism_put_proto
  USE mod_prism_def_partition_proto, ONLY: prism_def_partition_proto

  IMPLICIT NONE

  PUBLIC

  INTEGER, PARAMETER :: pwp = ip_realwp_p   ! working precision (psmile)

  INTEGER, ALLOCATABLE :: paral(:) ! parallel strategy
  INTEGER :: nsegments             ! no. of segments
  INTEGER :: parsize               ! 
  INTEGER :: nbtotproc             ! total no. of procs
  INTEGER :: nbcplproc = 1         ! no. of procs involved in data exchange
  INTEGER :: commlocal             ! local communicator
  INTEGER :: comp_id               ! model id
  INTEGER :: rank                  ! rank of processor in local communicator
  INTEGER :: info, ierror             ! message, error codes
  INTEGER :: var_shape(2)          !
  INTEGER :: part_id               ! 
  INTEGER :: var_nodims(2)         !

  CHARACTER (len=6) :: modnam = 'mpi-om' ! model name
#ifndef __cpl_co2
#ifdef __cpl_dust
  INTEGER, PARAMETER :: nflda2o = 16 ! no. of fields passed from atmosphere to ocean
#else
#ifdef ADDCONTRA
  INTEGER, PARAMETER :: nflda2o = 21 ! no. of fields passed from atmosphere to ocean
#else
  INTEGER, PARAMETER :: nflda2o = 15 ! no. of fields passed from atmosphere to ocean
#endif
#endif
#ifdef __cpl_dms
  INTEGER, PARAMETER :: nfldo2a = 8  ! no. of fields passed from ocean to atmosphere
#else
#ifdef ADDCONTRA
  INTEGER, PARAMETER :: nfldo2a = 9  ! no. of fields passed from ocean to atmosphere
#else
  INTEGER, PARAMETER :: nfldo2a = 6  ! no. of fields passed from ocean to atmosphere
#endif
#endif
#else /*__cpl_co2*/
#ifdef __cpl_dust
  INTEGER, PARAMETER :: nflda2o = 18 ! no. of fields passed from atmosphere to ocean
#else
#ifdef ADDCONTRA
  INTEGER, PARAMETER :: nflda2o = 20 ! no. of fields passed from atmosphere to ocean
#else
  INTEGER, PARAMETER :: nflda2o = 17 ! no. of fields passed from atmosphere to ocean
#endif
#endif
#ifdef __cpl_dms
  INTEGER, PARAMETER :: nfldo2a = 10 ! no. of fields passed from ocean to atmosphere
#else
#ifdef ADDCONTRA
  INTEGER, PARAMETER :: nfldo2a = 11  ! no. of fields passed from ocean to atmosphere
#else
  INTEGER, PARAMETER :: nfldo2a = 8  ! no. of fields passed from ocean to atmosphere
#endif
#endif
#endif /*__cpl_co2*/

  CHARACTER (len=8) :: clstr8a2o(nflda2o)   ! port names of exchange fields received
  CHARACTER (len=8) :: clstr8o2a(nfldo2a)   ! port names of exchange fields sent

  INTEGER    :: a2o_freq             ! exchange frequency in secs for receiving fields
  REAL(pwp)  :: o2a_time = 0.0_pwp   ! time passed since last couple_put_o2a
  INTEGER    :: o2a_freq             ! exchange frequency in secs for sending fields
 
  INTEGER :: run_date_secs = 0     ! time (sec) passed since start of run

  INTEGER*4 portin_id (nflda2o) ! Port IDs of exchange fields received
  INTEGER*4 portout_id(nfldo2a) ! Port IDs of exchange fields sent

  INTEGER :: nmseq                 ! run mode : sequential/concurrent = 2/1
  INTEGER :: num_fld_recvd = 0     ! counter for received fields
  INTEGER :: num_fld_sent  = 0     ! counter for fields sent

  INTEGER :: jfld

  REAL(wp), POINTER ::         gl_txwo(:,:), gl_txio(:,:), gl_wsvo(:,:), & 
                               gl_u1e(:,:), gl_v1e(:,:), gl_b1e(:,:), &
                               gl_z1e(:,:), gl_tywe(:,:), gl_tyie(:,:), &
                               gl_frio(:,:), gl_frwo(:,:), gl_rhio(:,:), &
                               gl_chio(:,:), gl_nhwo(:,:), gl_shwo(:,:)

#ifdef ADDCONTRA
  REAL(wp), POINTER :: gl_frwo16(:,:), gl_frwo18(:,:), gl_frwhdo(:,:), &
                       gl_frio16(:,:), gl_frio18(:,:), gl_frihdo(:,:)
#endif

#ifdef __cpl_dms
  REAL(wp), POINTER :: gl_dmsflux((:,:), gl_dms(:,:)
#endif
#ifdef __cpl_co2
  REAL(wp), POINTER :: gl_co2trans(:,:), gl_suppco2(:,:)
#endif


#ifdef FLUXCORRECT
  USE mo_commo_fluxcorr,  ONLY: fluko_hfl,fluko_frw
  REAL    ::  fhflmax, fhflmin
#endif

  CHARACTER (len=8) :: cnfileow    ! file to write ocean surface conditions.

!
! !DESCRIPTION:
!
! - contains all code related to coupling with oasis3
! - Note: this interface to the psmile library uses variables from
!   the psmile modules mod_comprism_proto.
!   In the next version (oasis4) these variables will all be received
!   by a call of a prism_get_... routine. Then an update of the code will
!   be needed. These variables are recognized by their prefix 'ig_'.
!
! !REVISION HISTORY:
! 03.05.22 Stephanie Legutke, MPI-Met, M&D
!          - created
! 03.08.07 Stephanie Legutke, MPI-Met, M&D
!          - update psmile
! 04.10.01 Noel Keenlyside, IFM-GEOMAR
!          - message passing version (accumulation done on global arrays)
!
! EOP
!-----------------------------------------------------------------------
! $Id: mo_couple.f90,v 1.4.2.1.10.1.4.2.2.3.4.1 2006/03/31 08:20:35 m211054 Exp $
!-----------------------------------------------------------------------
! 

CHARACTER(len=*), PARAMETER, PRIVATE :: cl_version_string = &
      '$Id: mo_couple.f90,v 1.4.2.1.10.1.4.2.2.3.4.1 2006/03/31 08:20:35 m211054 Exp $'

CONTAINS


!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  couple_prep
!
! !INTERFACE:
!
  SUBROUTINE couple_prep
!
! !USES:
!
! !RETURN VALUE:
!
  IMPLICIT NONE
!
! !PARAMETERS:
!
! !DESCRIPTION:
!
! - prepare coupling.
!   
! !REVISION HISTORY:
! 03.05.22  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Standard output file oceout is connected to unit ',io_stdout
      WRITE(io_stdout,*) ' '
      CALL FLUSH(io_stdout)
    ENDIF
#endif

  END SUBROUTINE couple_prep



!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  couple_correct_ini
!
! !INTERFACE:
!
  SUBROUTINE couple_correct_ini
!
! !USES:
!
! !RETURN VALUE:
!
    IMPLICIT NONE
!
! !PARAMETERS:
!
! !DESCRIPTION:
!
!- flux correction initialization.
!
! !REVISION HISTORY:
! 03.05.22  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

#ifdef FLUXCORRECT
    CALL FLUXC_INI
#endif /*FLUXCORRECT*/

#ifdef __synout
    IF (p_parallel_io) THEN
       WRITE(io_stdout,*) ' '
       WRITE(io_stdout,*) ' Flux correction initialized'
       WRITE(io_stdout,*) ' '
       CALL FLUSH(io_stdout)
    ENDIF
#endif

  END SUBROUTINE couple_correct_ini


!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  couple_init
!
! !INTERFACE:
!
  SUBROUTINE couple_init(actyear)
!
! !USES:
!
! !RETURN VALUE:
!
    IMPLICIT NONE
!
    INTEGER, INTENT(inout) :: actyear    ! year at start of run from restart file (IN)
                                       ! year at start of run from oasis (OUT)!
    ! !PARAMETERS:
    !
    ! !DESCRIPTION:
    !
    !- Initialize communication.
    !
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    !
    ! EOP
    !-----------------------------------------------------------------------
    !WRITE(0,*) ' couple_init 1'
#ifdef __synout
    IF (p_parallel_io) THEN
       WRITE(IO_STDOUT,*) ' '
       WRITE(IO_STDOUT,*) ' Start couple_init'
       WRITE(IO_STDOUT,*) ' *****************'
       WRITE(IO_STDOUT,*) ' '
       CALL FLUSH(IO_STDOUT)           
    ENDIF
#endif
    !    Allocate arrays for accumulation
    CALL alloc_mem_thacc

    CALL gather(AMSUO(:,:,1),AMSUO_G_L1,p_io)

    CALL gather(AMSUE(:,:,1),AMSUE_G_L1,p_io)

    IF (p_parallel_io) THEN
      ALLOCATE (gl_txwo(ie_g,je_g))
      ALLOCATE (gl_u1e (ie_g,je_g))
      ALLOCATE (gl_v1e (ie_g,je_g))
      ALLOCATE (gl_tywe(ie_g,je_g))
      ALLOCATE (gl_txio(ie_g,je_g))
      ALLOCATE (gl_b1e (ie_g,je_g))
      ALLOCATE (gl_z1e (ie_g,je_g))
      ALLOCATE (gl_tyie(ie_g,je_g))
      ALLOCATE (gl_frio(ie_g,je_g))
      ALLOCATE (gl_frwo(ie_g,je_g))
      ALLOCATE (gl_rhio(ie_g,je_g))
      ALLOCATE (gl_chio(ie_g,je_g))
      ALLOCATE (gl_nhwo(ie_g,je_g))
      ALLOCATE (gl_shwo(ie_g,je_g))
      ALLOCATE (gl_wsvo(ie_g,je_g))
#ifdef ADDCONTRA
      ALLOCATE (gl_frwo16(ie_g,je_g))
      ALLOCATE (gl_frwo18(ie_g,je_g))
      ALLOCATE (gl_frwhdo(ie_g,je_g))
      ALLOCATE (gl_frio16(ie_g,je_g))
      ALLOCATE (gl_frio18(ie_g,je_g))
      ALLOCATE (gl_frihdo(ie_g,je_g))
#endif
    else
      gl_txwo => NULL()
      gl_u1e  => NULL()
      gl_v1e  => NULL()
      gl_tywe => NULL()
      gl_txio => NULL()
      gl_b1e  => NULL()
      gl_z1e  => NULL()
      gl_tyie => NULL()
      gl_frio => NULL()
      gl_frwo => NULL()
      gl_rhio => NULL()
      gl_chio => NULL()
      gl_nhwo => NULL()
      gl_shwo => NULL()
      gl_wsvo => NULL()
#ifdef ADDCONTRA
      gl_frwo16 => NULL()
      gl_frwo18 => NULL()
      gl_frwhdo => NULL()
      gl_frio16 => NULL()
      gl_frio18 => NULL()
      gl_frihdo => NULL()
#endif
    endif
#ifdef __cpl_dms
    IF (p_parallel_io) THEN
       ALLOCATE (gl_dmsflux(ie_g,je_g))
       ALLOCATE (gl_dms(ie_g,je_g))
    ELSE
       gl_dmsflux => NULL()
       gl_dms     => NULL()
    END IF
#endif
#ifdef __cpl_co2
    IF (p_parallel_io) THEN
       ALLOCATE (gl_co2trans(ie_g,je_g))
       ALLOCATE (gl_suppco2(ie_g,je_g))
    ELSE
       gl_co2trans => NULL()
       gl_suppco2 => NULL()
    END IF
#endif

!WRITE(0,*) ' couple_init 2'

    CALL gather(aofltxwo, gl_txwo, p_io)
    CALL gather(u1e,      gl_u1e,  p_io)
    CALL gather(v1e,      gl_v1e,  p_io)
    CALL gather(aofltywe, gl_tywe, p_io)
    CALL gather(aofltxio, gl_txio, p_io)
    CALL gather(b1e,      gl_b1e,  p_io)
    CALL gather(z1e,      gl_z1e,  p_io)
    CALL gather(aofltyie, gl_tyie, p_io)
    CALL gather(aoflfrio, gl_frio, p_io)
    CALL gather(aoflfrwo, gl_frwo, p_io)
    CALL gather(aoflrhio, gl_rhio, p_io)
    CALL gather(aoflchio, gl_chio, p_io)
    CALL gather(aoflnhwo, gl_nhwo, p_io)
    CALL gather(aoflshwo, gl_shwo, p_io)
    CALL gather(aoflwsvo, gl_wsvo, p_io)
#ifdef ADDCONTRA
    CALL gather(aoflfrwo16, gl_frwo16, p_io)
    CALL gather(aoflfrwo18, gl_frwo18, p_io)
    CALL gather(aoflfrwhdo, gl_frwhdo, p_io)
    CALL gather(aoflfrio16, gl_frio16, p_io)
    CALL gather(aoflfrio18, gl_frio18, p_io)
    CALL gather(aoflfrihdo, gl_frihdo, p_io)
#endif
#ifdef __cpl_dms
    CALL gather(dmsflux,  gl_dmsflux, p_io)
    CALL gather(dms,      gl_dms,     p_io)
#endif
#ifdef __cpl_co2
    CALL gather(co2trans,  gl_co2trans, p_io)
    CALL gather(suppco2,  gl_suppco2, p_io)
#endif

#if ! (defined (use_comm_MPI1) || defined (use_comm_MPI2))
    CALL message('','In coupled mode specification of either')
    CALL message('','-Duse_comm_MP1 or -Duse_comm_MPI2 is required.')
    CALL finish('p_start','Aborting run, please recompile.')
#endif
  !
!WRITE(0,*) ' couple_init 3'
#ifdef use_comm_MPI2
  !-- Join the communicator group.
  !  
    ierror = PRISM_Ok
    CALL prism_init_comp_proto (comp_id, modnam, ierror)
    IF (ierror /= PRISM_Ok) THEN
       CALL couple_abort (modnam,'couple_init','pb init_comp',io_stderr)
    ENDIF
#endif
  !
!WRITE(0,*) ' couple_init 4'
#ifdef __synout
    IF (p_parallel_io) THEN
       WRITE(IO_STDOUT,*) ' '
       WRITE(IO_STDOUT,*) 'after prism_init_comp_proto ...'
       WRITE(IO_STDOUT,*) ' *****************'
       WRITE(IO_STDOUT,*) ' '
       CALL FLUSH(IO_STDOUT)           
    ENDIF
#endif

    !
    !-- Write grids file for oasis (only one thread)
    !
    CALL grids_writing(io_stdout)
    !WRITE(0,*) ' couple_init 5'
    !
    !-- receive local communicator
    !  
#ifdef use_comm_MPI1
    commlocal = p_all_comm
#endif
#ifdef use_comm_MPI2
    ierror = PRISM_Ok
    CALL prism_get_localcomm_proto(commlocal, ierror)
    IF (ierror /= PRISM_Ok) THEN
       CALL couple_abort (modnam,'couple_init','pb get_localcomm',io_stderr)
    ENDIF
#endif
    nbtotproc = p_nprocs
    rank      = p_pe
!lk      CALL MPI_Comm_Size(commlocal, nbtotproc, ierror)
!lk      CALL MPI_Comm_Rank(commlocal, rank, ierror)

    ! Fields passed from the ocean to the atmosphere
    !
    clstr8o2a(1:nfldo2a) = cg_cnaminp(1:nfldo2a)
    !
    ! Fields passed from the atmosphere to the ocean
    !
    clstr8a2o(1:nflda2o) = cg_cnamout(nfldo2a+1:nfldo2a+nflda2o)
    !
    !-- Calendar year
    ! 
    actyear = ig_inidate(1)
    !
    !-- Create OASIS initial file (ocean surface conditions).
    ! 
    cnfileow=cg_clim_rstfile(1)
    nmseq = MAXVAL (ig_clim_seq)
    IF ( nmseq > 1 ) THEN
       WRITE( io_stdout,* ) 'Calling ini_wrte to create ', cnfileow
       CALL ini_wrte(io_in_octl,io_stdout)
    ENDIF
    !WRITE(0,*) ' couple_init 5'
    !
    !-- Decomposition strategy
    !  
    nsegments = 1
    parsize = 3
    ALLOCATE(paral(parsize))
    paral ( clim_strategy ) = clim_serial 
    paral ( clim_length   ) = jpdim_oce
    paral ( clim_offset   ) = 0

    IF (p_parallel_io) THEN
       ierror   = PRISM_Ok
       CALL prism_def_partition_proto (part_id, paral, ierror)
       IF (ierror /= PRISM_Ok) THEN
          CALL couple_abort (modnam,'couple_init','pb def_partition',io_stderr)
       ENDIF

       !
       !-- Definitions for ports of incoming fields
       !  
       var_nodims(1) = 1
       var_nodims(2) = 0
       var_shape(1)= 1
       var_shape(2)= paral (clim_length)

       DO jfld = 1,nflda2o

          ierror = PRISM_Ok
          CALL prism_def_var_proto &
               (portin_id(jfld), clstr8a2o(jfld), &
               part_id, var_nodims, PRISM_in,  &
               var_shape, PRISM_Real, ierror )
          IF (ierror /= PRISM_Ok) THEN
             WRITE (io_stdout, *) &
                  'WARNING : Problem with import port ',clstr8a2o(jfld)
             WRITE (io_stdout, *) &
                  'WARNING : Port ID returned is : ',portin_id(jfld)
             WRITE (io_stdout, *) '=======   Error code number = ',ierror
          ELSE
             WRITE (io_stdout,*)'   '
             WRITE (io_stdout,*) &
                  ' couple_init: Import port ',clstr8a2o(jfld),' defined'
             WRITE (io_stdout,*) &
                  ' couple_init: Port ID returned is : ',portin_id(jfld)
             WRITE (io_stdout,*)' With exchange freq.    ', &
                  ig_def_freq(portin_id(jfld))
          ENDIF

       END DO

       !
       !-- Definitions for ports of outgoing fields
       !  
       DO jfld = 1,nfldo2a

          ierror = PRISM_Ok
          CALL prism_def_var_proto (portout_id(jfld), clstr8o2a(jfld), &
               part_id, var_nodims, PRISM_Out, &
               var_shape, PRISM_Real, ierror )
          IF (ierror /= PRISM_Ok) THEN
             WRITE (io_stdout, *) &
                  'WARNING : Problem with export port ',clstr8o2a(jfld)
             WRITE (io_stdout, *) &
                  'WARNING : Port ID returned is : ',portout_id(jfld)
             WRITE (io_stdout, *) '=======   Error code number = ',ierror
          ELSE
             WRITE (io_stdout,*)'   '
             WRITE (io_stdout,*) &
                  ' couple_init: Export port ',clstr8o2a(jfld),' defined'
             WRITE (io_stdout,*) &
                  ' couple_init: Port ID returned is : ',portout_id(jfld)
             WRITE (io_stdout,*)' With exchange freq.    ', &
                  ig_def_freq(portout_id(jfld))
          ENDIF
       END DO
       WRITE (io_stdout, *) ' '
      
       ierror = PRISM_Ok
       CALL prism_enddef_proto(ierror)
       IF (ierror /= PRISM_Ok) THEN
          WRITE (io_stdout, *) 'WARNING : Problem with prism_enddef_proto'
          WRITE (io_stdout, *) '=======   Error code number = ',ierror
       ENDIF
       !
       !-- Check 
       !   ...coupling-control parameters against those received from oasis.
       !
       CALL chck_par
       !     
       !-- Display control parameters received from oasis3.
       !   ------------------------------------------------
       WRITE(io_stdout,*)' Couple_init:'
       WRITE(io_stdout,*)' -----------------------------------------------------'
       WRITE(io_stdout,*)' model name                = ',modnam
       WRITE(io_stdout,*)' model time step           = ',dt
       WRITE(io_stdout,*)' model restart-file year   = ',lyear1
       WRITE(io_stdout,*)' -----------------------------------------------------'
       WRITE(io_stdout,*)' '
       WRITE(io_stdout,*)' Parameters received from oasis3 :'
       WRITE(io_stdout,*)' mynum                         ', mynum
       WRITE(io_stdout,*)' mytid                         ', mytid
       WRITE(io_stdout,*)' model number                  ', comp_id
       WRITE(io_stdout,*)' ig_mynummod                   ', ig_mynummod
       WRITE(io_stdout,*)' Total no. of processors       ', nbtotproc
       WRITE(io_stdout,*)' No. of procs for data exchange', nbcplproc
       WRITE(io_stdout,*)' depomposition ID              ', part_id
       WRITE(io_stdout,*)' total time of the simulation  ', ig_ntime
       WRITE(io_stdout,*)' number of fields exchanged    ', ig_clim_nfield
       WRITE(io_stdout,*)' initial date                  ', ig_inidate
       WRITE(io_stdout,*)' Lag of exported fields        ', ig_def_lag
       WRITE(io_stdout,*)' coupling period of fields     ', ig_def_freq
       WRITE(io_stdout,*)' sequential index of fields    ', ig_def_seq
       WRITE(io_stdout,*)' ig_def_norstfile              ', ig_def_norstfile
       WRITE(io_stdout,*)' number of restart files       ', ig_nbr_rstfile
       WRITE(io_stdout,*)' name of restart files         ', cg_def_rstfile
       WRITE(io_stdout,*)' restart files in netcdf       ', lg_ncdfrst
!!lk       WRITE(io_stdout,*)' ig_aux                        ', ig_aux
       WRITE(io_stdout,*)' name of exchanged fields (in) ', cg_cnaminp
       WRITE(io_stdout,*)' name of exchanged fields (out)', cg_cnamout
       WRITE(io_stdout,*)' cg_ignout_field               ', cg_ignout_field
       WRITE(io_stdout,*)' any field going through Oasis?', lg_oasis_field
       WRITE(io_stdout,*)' seconds between 2 sent        ', o2a_freq
       WRITE(io_stdout,*)' seconds between 2 receive     ', a2o_freq
       WRITE(io_stdout,*)' calender run year             ', actyear
       IF(nmseq /= 1) THEN
          WRITE(io_stdout,*)' the models run sequentially !'
          WRITE(io_stdout,*)' no. of sequential fields: ',nmseq
       ELSE
          WRITE(io_stdout,*)' all fields have sequential no. 1 !'
          WRITE(io_stdout,*)' the models run concurrently !'
       ENDIF

#ifndef __accu_by_psmile
       ! (Transfer instantaneous field => ig_clim_trans=1
       IF (SUM(ig_clim_trans(1:nfldo2a)) .NE. nfldo2a ) THEN
         WRITE (io_stdout,*) 'ERROR: Wrong time transformation for o2a fields!'
         WRITE (io_stdout,*) 'Accumulation by mpi-om requires to specify'
         WRITE (io_stdout,*) 'timtranso2a=INSTANT for all fields in the run script.'
         CALL flush (io_stdout)
         CALL p_abort
       ENDIF
#else
       ! (Transfer averaged field => ig_clim_trans=2
       IF (SUM(ig_clim_trans(1:nfldo2a)) .NE. 2*nfldo2a ) THEN
         WRITE (io_stdout,*) 'ERROR: Wrong time transformation for o2a fields!'
         WRITE (io_stdout,*) 'Accumulation by psmile requires to specify'
         WRITE (io_stdout,*) 'timtranso2a=AVERAGE for all fields in the run script.'
         CALL flush (io_stdout)
         CALL p_abort
       ENDIF
#endif

    ENDIF  ! p_parallel_io
    !
    !-- Initilize exchange fields
    !
    CALL ini_o2a(io_stdout)
    !
    CALL p_bcast(o2a_freq,p_io)
    CALL p_bcast(o2a_time,p_io)

#ifdef __synout
    IF (p_parallel_io) THEN
       WRITE(io_stdout,*) ' '
       WRITE(io_stdout,*) ' End of couple_init'
       WRITE(io_stdout,*) ' ******************'
       WRITE(io_stdout,*) ' '
       CALL FLUSH(io_stdout)
    ENDIF
#endif

  END SUBROUTINE couple_init


!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  couple_get_a2o
!
! !INTERFACE:
!
  SUBROUTINE couple_get_a2o(ldtrun)
!
! !USES:
!
! !RETURN VALUE:
!
  IMPLICIT NONE
  INTEGER, INTENT(IN)   ::  ldtrun ! model run step
!
! !PARAMETERS:
!
   REAL(pwp)  :: exfld(jpdim_ocei,jpdim_ocej) ! buffer for receiving exchange fields
!
! !DESCRIPTION:
!
!- Get data from coupler.
!
! !REVISION HISTORY:
! 03.05.22  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

#ifdef __synout
     IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Start of couple_get_a2o'
      WRITE(io_stdout,*) ' ***********************'
      WRITE(io_stdout,*) ' '
      CALL FLUSH(io_stdout)           
     ENDIF
#endif

      !
      !-- Import exchange fields; try at all dates
      ! 

      DO jfld = 1,nflda2o

        IF (p_parallel_io) THEN
#ifdef __synout
         WRITE(io_stdout,*)'    '
         WRITE(io_stdout,*)' ==> prism_get for ',clstr8a2o(jfld)
         WRITE(io_stdout,*) &
         ' ==> at run-step ',ldtrun,' (seconds passed=',run_date_secs,')'
         CALL FLUSH(io_stdout)
#endif
         info = PRISM_Ok
         CALL prism_get_proto(portin_id(jfld),run_date_secs,exfld,info)

         IF (info == PRISM_Recvd        .OR. &
             info == PRISM_FromRest     .OR. &
             info == PRISM_RecvOut      .OR. &
             info == PRISM_FromRestOut         ) num_fld_recvd = num_fld_recvd + 1

         CALL digest_get_Id &
            (io_stdout,io_stderr,info,clstr8a2o(jfld),run_date_secs,num_fld_recvd)

        ENDIF
        CALL p_bcast(num_fld_recvd,p_io)
        CALL p_bcast(info, p_io)

         IF (info == PRISM_Recvd        .OR. &
             info == PRISM_FromRest     .OR. &
             info == PRISM_RecvOut      .OR. &
             info == PRISM_FromRestOut         ) &
            CALL PUT_A2O(clstr8a2o(jfld),exfld,jpdim_ocei,jpdim_ocej,io_stdout)
      ENDDO

      !
      !-- If all fields have been received ...postprocessing
      !

      IF ( num_fld_recvd == nflda2o )      &
         CALL post_a2o(num_fld_recvd, io_stdout, ldtrun)

#ifdef __synout
     IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' End of couple_get_a2o'
      WRITE(io_stdout,*) ' *********************'
      WRITE(io_stdout,*) ' '
      CALL FLUSH(io_stdout)
     ENDIF
#endif

  END SUBROUTINE couple_get_a2o


!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  couple_put_o2a
!
! !INTERFACE:
!
  SUBROUTINE couple_put_o2a(ldtrun)
!
! !USES:
!
! !RETURN VALUE:
!
    IMPLICIT NONE
    INTEGER, INTENT(IN)   ::  ldtrun ! model run step
    !
    ! !PARAMETERS:
    !
    REAL(pwp)  :: exfld(jpdim_ocei,jpdim_ocej) ! buffer array for sending.
    !
    ! !DESCRIPTION:
    !
    !- Send data to the coupler. All o2a fields are treated here.
    !- The exchange fields is first copied to the sending buffer.
    !- Then prism_put is called with the buffer.
    !- Two options:
    !    with cpp accu_by_psmile :
    !             accumulation and normalisation is done by psmile.
    !    else:
    !             accumulation and normalisation is done by the model.
    !
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    !
    ! EOP
    !-----------------------------------------------------------------------

#ifdef __synout
    IF (p_parallel_io) THEN
       WRITE(io_stdout,*) ' '
       WRITE(io_stdout,*) ' Start of couple_put_o2a'
       WRITE(io_stdout,*) ' ***********************'
       WRITE(io_stdout,*) ' '
       CALL FLUSH(io_stdout)           
    END IF
#endif

#ifndef __accu_by_psmile
    !
    !-- Accumulate exchange data
    !

    CALL acc_o2a(io_stdout)

    !
    !-- At the end of every coupled time step : normalize
    !
    IF ( MOD(INT(o2a_time), o2a_freq) == 0) THEN
       IF (p_parallel_io) CALL avg_o2a(io_stdout)
    ENDIF
#else
    !
    !-- Prepare fields for transfer to psmile
    !
    CALL prep_o2a(io_stdout)
#endif

    !
    !-- Export exchange fields; call prism_put every model step
    !  

!   WRITE(io_stdout,*) 'WISO-TEST MPIOM PUT_O2A 830 GL_O16ACC', SUM(gl_o16acc)
!   WRITE(io_stdout,*) 'WISO-TEST MPIOM PUT_O2A 830 GL_O18ACC', SUM(gl_o18acc)
!   WRITE(io_stdout,*) 'WISO-TEST MPIOM PUT_O2A 830 GL_HDOACC', SUM(gl_hdoacc)

    DO jfld = 1,nfldo2a

#ifndef __accu_by_psmile
       IF (MOD(INT(o2a_time),o2a_freq) == 0) &
#endif
            CALL GET_O2A(clstr8o2a(jfld),exfld,jpdim_ocei,jpdim_ocej,io_stdout)

       IF (p_parallel_io) THEN
#ifdef __synout
          WRITE(io_stdout,*)'    '
          WRITE(io_stdout,*)' ==> prism_put for ',clstr8o2a(jfld)
          WRITE(io_stdout,*) &
               ' ==> at run-step ',ldtrun,' (seconds passed=',run_date_secs,')'
          CALL FLUSH(io_stdout)
#endif
          info = PRISM_Ok
          CALL prism_put_proto(portout_id(jfld),run_date_secs,exfld,info)

          IF (info == PRISM_Sent     .OR. &
               info == PRISM_ToRest   .OR. &
               info == PRISM_SentOut  .OR. &
               info == PRISM_ToRestOut       ) num_fld_sent = num_fld_sent + 1

          CALL digest_put_Id &
               (io_stdout,io_stderr,info,clstr8o2a(jfld),run_date_secs,num_fld_sent)
       END IF
    ENDDO

#ifndef __accu_by_psmile
    !
    !--  Reset exchange fields / accumulation interval
    !  

    IF ( MOD(INT(o2a_time),o2a_freq) == 0    .AND. &
         ldtrun*INT(dt) /= ig_ntime) THEN
       CALL ini_o2a(io_stdout)
    ENDIF
#endif
    !
    !-- Reset counter when fields are complete.
    !  
    IF ( num_fld_sent == nfldo2a ) THEN

       IF (p_parallel_io) WRITE(io_stdout,*) ' All ',nfldo2a,' fields sent.'
       num_fld_sent = 0

    ENDIF


#ifdef __synout
    IF (p_parallel_io) THEN
       WRITE(io_stdout,*) ' '
       WRITE(io_stdout,*) ' End of couple_put_o2a'
       WRITE(io_stdout,*) ' *********************'
       WRITE(io_stdout,*) ' '
       CALL FLUSH(io_stdout)
    END IF !p_parallel_io
#endif

    ! Update o2a_time
    CALL p_bcast(o2a_time,p_io)

  END SUBROUTINE couple_put_o2a

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  couple_end
!
! !INTERFACE:
!
  SUBROUTINE couple_end
!
! !USES:
!
! !RETURN VALUE:
!
  IMPLICIT NONE
!
! !PARAMETERS:
!
!
! !DESCRIPTION:
!
!- Stop communication with oasis.
!
! !REVISION HISTORY:
! 03.05.22  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

     IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Start of couple_end '
      WRITE(io_stdout,*) ' *******************'
      WRITE(io_stdout,*) ' '
      CALL FLUSH(io_stdout)           
     END IF

  !
  !-- Deallocation of memory
  !  
      ierror = PRISM_Ok
      CALL prism_terminate_proto(ierror)
      IF (ierror /= PRISM_Ok) THEN
         WRITE (io_stdout, *) 'An error occured couple_end: Error = ',ierror
      ENDIF

#ifdef __synout
     IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' End of couple_end'
      WRITE(io_stdout,*) ' *****************'
      WRITE(io_stdout,*) ' '
      CALL FLUSH(io_stdout)           
     END IF
#endif

  END SUBROUTINE couple_end

!-----------------------------------------------------------------------
! BOP
!
! !ROUTINE:  couple_abort
!
! !INTERFACE:

    SUBROUTINE couple_abort(str1,str2,str3,kerr)

    IMPLICIT NONE

    INTEGER,          INTENT(in) :: kerr      ! std error out unit
    CHARACTER(len=*), INTENT(in) :: str1      ! name of calling model
    CHARACTER(len=*), INTENT(in) :: str2      ! name of calling routine
    CHARACTER(len=*), INTENT(in) :: str3      ! error message

! !DESCRIPTION:
!
! - stop the MPI application with messages.
!
! !REVISION HISTORY:
! 03.08.07  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------
! $Id: mo_couple.f90,v 1.4.2.1.10.1.4.2.2.3.4.1 2006/03/31 08:20:35 m211054 Exp $
!-----------------------------------------------------------------------
! 

!lk      INCLUDE 'mpif.h'      

!lk      INTEGER :: ierr           ! error code

    CALL prism_abort_proto (comp_id, TRIM(str2), TRIM(str3))

!lk      WRITE (kerr,'(a,a,a,a,a)') 'Calling MPI_ABORT in ',str2,' by model ',str1,':'
!lk      WRITE (kerr,'(a,a)')   '        ',str3

!lk      CALL MPI_ABORT (MPI_COMM_WORLD, 0, ierr)

!lk      IF (ierr /= MPI_SUCCESS) THEN
!lk         WRITE (kerr,'(a)') ' MPI_ABORT failed'
!lk         WRITE (kerr,'(a,i4)') ' Error =  ', ierr
!lk         STOP
!lk      END IF
!lk      CALL FLUSH(kerr)

      END SUBROUTINE couple_abort

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  grids_writing
!
! !INTERFACE:

  SUBROUTINE grids_writing(kout)

! !RETURN VALUE:

! !USES:
  USE MO_PARAM1,            ONLY: ie_g, je_g 
  USE MO_COMMO1,            ONLY: gila_g, giph_g, weto_g, dlxp_g    & 
                             ,dlyp_g, dlxu_g, dlxv_g, dlyu_g, dlyv_g
  USE mod_prism_grids_writing

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: kout        ! unit number for standard output
                                     
  INTEGER       ::   gwrite           ! flag to state whether grids writing is
                                      ! needed or not (1 / 0)
  INTEGER       ::   i,j              ! looping indicees
  INTEGER       ::   ip1,im1          ! i+1, i-1
  INTEGER       ::   mask(ie_g-2,je_g)    ! inverse land sea mask
  INTEGER       ::   msks(ie_g-2,je_g)    ! land sea mask of scalar grid
  INTEGER       ::   msku(ie_g-2,je_g)    ! land sea mask of vector-u grid
  INTEGER       ::   mskv(ie_g-2,je_g)    ! land sea mask of vector-v grid
  REAL          ::   lon(2*ie_g,2*je_g)   ! longitudes of doubled array
  REAL          ::   lat(2*ie_g,2*je_g)   ! latitudes of doubled array
  REAL          ::   lons(ie_g-2,je_g)    ! longitudes of scalars
  REAL          ::   lats(ie_g-2,je_g)    ! latitudes of scalars
  REAL          ::   lonu(ie_g-2,je_g)    ! longitudes of vector u-components
  REAL          ::   latu(ie_g-2,je_g)    ! latitudes of vector u-components
  REAL          ::   lonv(ie_g-2,je_g)    ! longitudes of vector v-components
  REAL          ::   latv(ie_g-2,je_g)    ! latitudes of vector v-components
  REAL          ::   lonc(ie_g-2,je_g)    ! corner grid points of scalar grid
  REAL          ::   latc(ie_g-2,je_g)    ! corner grid points of scalar grid
  REAL          ::   clons(ie_g-2,je_g,4) ! corner longitudes of scalars
  REAL          ::   clats(ie_g-2,je_g,4) ! corner latitudes of scalars
  REAL          ::   clonu(ie_g-2,je_g,4) ! corner longitudes of vector u-components
  REAL          ::   clatu(ie_g-2,je_g,4) ! corner latitudes of vector u-components
  REAL          ::   clonv(ie_g-2,je_g,4) ! corner longitudes of vector v-components
  REAL          ::   clatv(ie_g-2,je_g,4) ! corner latitudes of vector v-components
  REAL          ::   areas(ie_g-2,je_g)   ! grid cell area of scalar
  REAL          ::   areau(ie_g-2,je_g)   ! grid cell area of vector u-components
  REAL          ::   areav(ie_g-2,je_g)   ! grid cell area of vector v-components
  REAL          ::   pi               ! PI
  REAL          ::   pi180            ! 180/PI
  CHARACTER*4   ::   grdacr           ! grid acronym (as used in namcouple)

! !DESCRIPTION:
!
! - Write grids and masks file for oasis.
!
!  WRITE(0,*)'grids_writing 1 '
!
! !REVISION HISTORY:
! July 7, 2003  V. Gayler - created
!
! EOP
!-----------------------------------------------------------------------
#ifdef __synout
   WRITE(kout,*)' '
   WRITE(kout,*)' Start of grids_writing'
   WRITE(kout,*)' **********************'
   WRITE(kout,*)' '
#endif

!  WRITE(0,*)'grids_writing 3a '

  IF ( p_parallel_io ) THEN
! Write grids out only on p_parallel_io
!
!-- start writing the arrays
!
#ifdef __synout
  WRITE(kout,*)'prism_start_grids_writing'
#endif

!  WRITE(0,*)'grids_writing 4 '
  call prism_start_grids_writing(gwrite)
!  WRITE(0,*)'grids_writing 5 '

  IF (gwrite == 1) THEN
!
!-- create arrays of longitudes and latitudes
!
#if defined (__noinvert)
     DO j = 1, 2*je_g
        lon(:,j)=gila_g(:,j)
        lat(:,j)=giph_g(:,j)
     ENDDO
#else
!
!--  follow the OASIS conventions: S -> N 
!
     DO j = 1, 2*je_g
        lon(:,j)=gila_g(:,2*je_g+1-j)
        lat(:,j)=giph_g(:,2*je_g+1-j)
     ENDDO
#endif
!
!--  convert from radiant to degree
!
     pi=ATAN(1.)*4.0
     pi180=180./pi
     lon(:,:)=lon(:,:) * pi180
     lat(:,:)=lat(:,:) * pi180

     WHERE (lon(:,:) < 0.)
        lon(:,:) = lon(:,:) + 360.
     END WHERE
!
!--  extract scalar/vector grid points
!
!     2*ij                            
!      :                              s: scalar
!      :                              u: vector-u
!      5   s  u  s  u  s  u           v: vector-v
!      4   v  c  v  c  v  c           c: grid cell corners of scalars
!      3   s  u  s  u  s  u
!      2   v  c  v  c  v  c
!      1   s  u  s  u  s  u           Line 0 and 1 are identical
!          v  c  v  c  v  c
!         
!          1  2  3  4  5  6 ... 2*ie
!
!
     DO i = 1, ie_g-2
        DO j = 1, je_g
#if defined (__noinvert)
!--     scalar
           lats(i,j) = lat(i*2+2,j*2)
           lons(i,j) = lon(i*2+2,j*2)
!--     vector - u                  
           latu(i,j) = lat(i*2+3,j*2)
           lonu(i,j) = lon(i*2+3,j*2)
        ENDDO
!--     vector - v                  
        DO j = 1, je_g-1
           latv(i,j) = lat(i*2+2,j*2+1)
           lonv(i,j) = lon(i*2+2,j*2+1)
        ENDDO
        latv(i,je_g) = lat(i*2+2,j*2)
        lonv(i,je_g) = lon(i*2+2,j*2)
!--     corners of scalar grid cells                  
        DO j = 1, je_g-1
           latc(i,j) = lat(i*2+3,j*2+1)
           lonc(i,j) = lon(i*2+3,j*2+1)
        ENDDO
        latc(i,je_g) = lat(i*2+3,j*2)
        lonc(i,je_g) = lon(i*2+3,j*2)
#else
!--     scalar
           lats(i,j) = lat(i*2+2,j*2-1)
           lons(i,j) = lon(i*2+2,j*2-1)
!--     vector - u                  
           latu(i,j) = lat(i*2+3,j*2-1)
           lonu(i,j) = lon(i*2+3,j*2-1)
        ENDDO
!--     vector - v                  
        latv(i,1) = lat(i*2+2,1)
        lonv(i,1) = lon(i*2+2,1)
        DO j = 2, je_g
           latv(i,j) = lat(i*2+2,j*2-2)
           lonv(i,j) = lon(i*2+2,j*2-2)
        ENDDO
!--     corners of scalar grid cells                  
        latc(i,1) = lat(i*2+3,1)
        lonc(i,1) = lon(i*2+3,1)
        DO j = 2, je_g
           latc(i,j) = lat(i*2+3,j*2-2)
           lonc(i,j) = lon(i*2+3,j*2-2)
        ENDDO
#endif
     ENDDO
!
!--  create corner arrays for SCRIP interpolation
!
     DO i = 1, ie_g-2
#if defined (__noinvert)
!
!--     scalar
!
        DO j = 1, je_g-1
           clons(i,j,1) = lon(2*i+3,2*j-1)
           clons(i,j,2) = lon(2*i+1,2*j-1)
           clons(i,j,3) = lon(2*i+1,2*j+1)
           clons(i,j,4) = lon(2*i+3,2*j+1)
           clats(i,j,1) = lat(2*i+3,2*j-1)
           clats(i,j,2) = lat(2*i+1,2*j-1)
           clats(i,j,3) = lat(2*i+1,2*j+1)
           clats(i,j,4) = lat(2*i+3,2*j+1)
        ENDDO
        j = je_g
        clons(i,j,1) = lon(2*i+3,2*j-1)
        clons(i,j,2) = lon(2*i+1,2*j-1)
        clons(i,j,3) = lon(2*i+1,2*j  )
        clons(i,j,4) = lon(2*i+3,2*j  )
        clats(i,j,1) = lat(2*i+3,2*j-1)
        clats(i,j,2) = lat(2*i+1,2*j-1)
        clats(i,j,3) = lat(2*i+1,2*j  )
        clats(i,j,4) = lat(2*i+3,2*j  )
!
!--     vector - u
!
        DO j = 1, je_g-1
           clonu(i,j,1) = lon(2*i+4,2*j-1)
           clonu(i,j,2) = lon(2*i+2,2*j-1)
           clonu(i,j,3) = lon(2*i+2,2*j+1)
           clonu(i,j,4) = lon(2*i+4,2*j+1)
           clatu(i,j,1) = lat(2*i+4,2*j-1)
           clatu(i,j,2) = lat(2*i+2,2*j-1)
           clatu(i,j,3) = lat(2*i+2,2*j+1)
           clatu(i,j,4) = lat(2*i+4,2*j+1)
        ENDDO
        j = je_g
        clonu(i,j,1) = lon(2*i+4,2*j-1)
        clonu(i,j,2) = lon(2*i+2,2*j-1)
        clonu(i,j,3) = lon(2*i+2,2*j  )
        clonu(i,j,4) = lon(2*i+4,2*j  )
        clatu(i,j,1) = lat(2*i+4,2*j-1)
        clatu(i,j,2) = lat(2*i+2,2*j-1)
        clatu(i,j,3) = lat(2*i+2,2*j  )
        clatu(i,j,4) = lat(2*i+4,2*j  )
!
!--     vector - v
!
        DO j = 1, je_g-2
           clonv(i,j,1) = lon(2*i+3,2*j  )
           clonv(i,j,2) = lon(2*i+1,2*j  )
           clonv(i,j,3) = lon(2*i+1,2*j+2)
           clonv(i,j,4) = lon(2*i+3,2*j+2)
           clatv(i,j,1) = lat(2*i+3,2*j  )
           clatv(i,j,2) = lat(2*i+1,2*j  )
           clatv(i,j,3) = lat(2*i+1,2*j+2)
           clatv(i,j,4) = lat(2*i+3,2*j+2)
        ENDDO
        j = je_g-1
        clonv(i,j,1) = lon(2*i+3,2*j  )
        clonv(i,j,2) = lon(2*i+1,2*j  )
        clonv(i,j,3) = lon(2*i+1,2*j+1)
        clonv(i,j,4) = lon(2*i+3,2*j+1)
        clatv(i,j,1) = lat(2*i+3,2*j  )
        clatv(i,j,2) = lat(2*i+1,2*j  )
        clatv(i,j,3) = lat(2*i+1,2*j+1)
        clatv(i,j,4) = lat(2*i+3,2*j+1)
        j = je_g
        clonv(i,j,1) = lon(2*i+3,2*j-1)
        clonv(i,j,2) = lon(2*i+1,2*j-1)
        clonv(i,j,3) = lon(2*i+1,2*j  )
        clonv(i,j,4) = lon(2*i+3,2*j  )
        clatv(i,j,1) = lat(2*i+3,2*j-1)
        clatv(i,j,2) = lat(2*i+1,2*j-1)
        clatv(i,j,3) = lat(2*i+1,2*j  )
        clatv(i,j,4) = lat(2*i+3,2*j  )
#else
        im1 = i-1
        ip1 = i+1
        IF (im1 == 0) im1 = ie_g-2
        IF (ip1 == ie_g-1) ip1 = 1
!
!--     scalar
!
        DO j = 1, je_g-1
           clons(i,j,1) = lonc(i  ,j+1)
           clons(i,j,2) = lonc(im1,j+1)
           clons(i,j,3) = lonc(im1,j  )
           clons(i,j,4) = lonc(i  ,j  )
           clats(i,j,1) = latc(i  ,j+1)
           clats(i,j,2) = latc(im1,j+1)
           clats(i,j,3) = latc(im1,j  )
           clats(i,j,4) = latc(i  ,j  )
        ENDDO
        clons(i,je_g,1) = lonu(i  ,je_g)
        clons(i,je_g,2) = lonu(im1,je_g)
        clons(i,je_g,3) = lonc(im1,je_g)
        clons(i,je_g,4) = lonc(i  ,je_g)
        clats(i,je_g,1) = latu(i  ,je_g)
        clats(i,je_g,2) = latu(im1,je_g)
        clats(i,je_g,3) = latc(im1,je_g)
        clats(i,je_g,4) = latc(i  ,je_g)
!
!--     vector - u
!
        DO j = 1, je_g-1
           clonu(i,j,1) = lonv(ip1,j+1)
           clonu(i,j,2) = lonv(i  ,j+1)
           clonu(i,j,3) = lonv(i  ,j  )
           clonu(i,j,4) = lonv(ip1,j  )
           clatu(i,j,1) = latv(ip1,j+1)
           clatu(i,j,2) = latv(i  ,j+1)
           clatu(i,j,3) = latv(i  ,j  )
           clatu(i,j,4) = latv(ip1,j  )
        ENDDO
        clonu(i,je_g,1) = lons(ip1,je_g)
        clonu(i,je_g,2) = lons(i  ,je_g)
        clonu(i,je_g,3) = lonv(i  ,je_g)
        clonu(i,je_g,4) = lonv(ip1,je_g)
        clatu(i,je_g,1) = lats(ip1,je_g)
        clatu(i,je_g,2) = lats(i  ,je_g)
        clatu(i,je_g,3) = latv(i  ,je_g)
        clatu(i,je_g,4) = latv(ip1,je_g)
!
!--     vector - v
!
        DO j = 3, je_g
           clonv(i,j,1) = lonu(i  ,j  )
           clonv(i,j,2) = lonu(im1,j  )
           clonv(i,j,3) = lonu(im1,j-1)
           clonv(i,j,4) = lonu(i  ,j-1)
           clatv(i,j,1) = latu(i  ,j  )
           clatv(i,j,2) = latu(im1,j  )
           clatv(i,j,3) = latu(im1,j-1)
           clatv(i,j,4) = latu(i  ,j-1)
        ENDDO
        clonv(i,1,1) = lonc(i  ,2)
        clonv(i,1,2) = lonc(im1,2)
        clonv(i,1,3) = lonu(im1,1)
        clonv(i,1,4) = lonu(i  ,1)
        clatv(i,1,1) = latc(i  ,2)
        clatv(i,1,2) = latc(im1,2)
        clatv(i,1,3) = latu(im1,1)
        clatv(i,1,4) = latu(i  ,1)
        clonv(i,2,1) = lonu(i  ,2)
        clonv(i,2,2) = lonu(im1,2)
        clonv(i,2,3) = lonc(im1,2)
        clonv(i,2,4) = lonc(i  ,2)
        clatv(i,2,1) = latu(i  ,2)
        clatv(i,2,2) = latu(im1,2)
        clatv(i,2,3) = latc(im1,2)
        clatv(i,2,4) = latc(i  ,2)
#endif
     ENDDO

!
!-- create mask arrays, following the OASIS convention
!
     msks(:,:) = 0
     DO j = 1, je_g
#if defined (__noinvert)
        mask(:,j)=weto_g(2:ie_g-1, j, 1)
#else
        mask(:,j)=weto_g(2:ie_g-1, je_g+1-j, 1)
#endif
     ENDDO
     WHERE (mask(:,:) == 0)
        msks = 1
     END WHERE

     msku(:,:) = 0
     DO j = 1, je_g
#if defined (__noinvert)
        mask(:,j)=AMSUO_G_L1(2:ie_g-1, j)
#else
        mask(:,j)=AMSUO_G_L1(2:ie_g-1, je_g+1-j)
#endif
     ENDDO
     WHERE (mask(:,:) == 0)
        msku = 1
     END WHERE

     mskv(:,:) = 0
     DO j = 1, je_g
#if defined (__noinvert)
        mask(:,j)=AMSUE_G_L1(2:ie_g-1, j)
#else
        mask(:,j)=AMSUE_G_L1(2:ie_g-1, je_g+1-j)
#endif
     ENDDO
     WHERE (mask(:,:) == 0)
        mskv = 1
     END WHERE
  
!
!-- create area arrays, following the OASIS convention
!
     DO j = 1, je_g
#if defined (__noinvert)
        areas(:,j)=dlxp_g(2:ie_g-1,j)*dlyp_g(2:ie_g-1,j)
        areau(:,j)=dlxu_g(2:ie_g-1,j)*dlyu_g(2:ie_g-1,j)
        areav(:,j)=dlxv_g(2:ie_g-1,j)*dlyv_g(2:ie_g-1,j)
#else
        areas(:,j)=dlxp_g(2:ie_g-1,je_g+1-j)*dlyp_g(2:ie_g-1,je_g+1-j)
        areau(:,j)=dlxu_g(2:ie_g-1,je_g+1-j)*dlyu_g(2:ie_g-1,je_g+1-j)
        areav(:,j)=dlxv_g(2:ie_g-1,je_g+1-j)*dlyv_g(2:ie_g-1,je_g+1-j)
#endif
     ENDDO
     WHERE (msks(:,:) == 1)
        areas(:,:) = 0
     END WHERE
     WHERE (msku(:,:) == 1)
        areau(:,:) = 0
     END WHERE
     WHERE (mskv(:,:) == 1)
        areav(:,:) = 0
     END WHERE

!
!-- write grids
!
     grdacr='oces'
#ifdef __synout
     WRITE(kout,*)'prism_write_grid: ', grdacr
#endif
     call prism_write_grid (grdacr, ie_g-2, je_g, lons(:,:), lats(:,:))
     grdacr='oceu'
#ifdef __synout
     WRITE(kout,*)'prism_write_grid: ', grdacr
#endif
     call prism_write_grid (grdacr, ie_g-2, je_g, lonu(:,:), latu(:,:))
     grdacr='ocev'
#ifdef __synout
     WRITE(kout,*)'prism_write_grid: ', grdacr
#endif
     call prism_write_grid (grdacr, ie_g-2, je_g, lonv(:,:), latv(:,:))
!
!-- write corners
!
!   writing corners is optional. If they are missing in the grids file, the
!   corners are calculated by scrip (with good results).
     grdacr='oces'
#ifdef __synout
     WRITE(kout,*)'prism_write_corner: ', grdacr
#endif
     call prism_write_corner (grdacr, ie_g-2, je_g, 4, clons(:,:,:), clats(:,:,:))
     grdacr='oceu'
#ifdef __synout
     WRITE(kout,*)'prism_write_corner: ', grdacr
#endif
     call prism_write_corner (grdacr, ie_g-2, je_g, 4, clonu(:,:,:), clatu(:,:,:))
     grdacr='ocev'
#ifdef __synout
     WRITE(kout,*)'prism_write_corner: ', grdacr
#endif
     call prism_write_corner (grdacr, ie_g-2, je_g, 4, clonv(:,:,:), clatv(:,:,:))
!
!-- write masks
!
     grdacr='oces'
#ifdef __synout
     WRITE(kout,*)'prism_write_mask: ', grdacr
#endif
     call prism_write_mask (grdacr, ie_g-2, je_g, msks(:,:))
     grdacr='oceu'
#ifdef __synout
     WRITE(kout,*)'prism_write_mask: ', grdacr
#endif
     call prism_write_mask (grdacr, ie_g-2, je_g, msku(:,:))
     grdacr='ocev'
#ifdef __synout
     WRITE(kout,*)'prism_write_mask: ', grdacr
#endif
     call prism_write_mask (grdacr, ie_g-2, je_g, mskv(:,:))
!
!-- write areas
!
     grdacr='oces'
#ifdef __synout
     WRITE(kout,*)'prism_write_area: ', grdacr
#endif
     call prism_write_area (grdacr, ie_g-2, je_g, areas(:,:))
     grdacr='oceu'
#ifdef __synout
     WRITE(kout,*)'prism_write_area: ', grdacr
#endif
     call prism_write_area (grdacr, ie_g-2, je_g, areau(:,:))
     grdacr='ocev'
#ifdef __synout
     WRITE(kout,*)'prism_write_area: ', grdacr
#endif
     call prism_write_area (grdacr, ie_g-2, je_g, areav(:,:))
!
!-- all arrays are written
!
#ifdef __synout
     WRITE(kout,*)'prism_terminate_grids_writing'
#endif
     call prism_terminate_grids_writing

#ifdef __synout
  ELSE
     WRITE(kout,*)'  grids files are existing'
     WRITE(kout,*)'  no writing needed'
#endif

  ENDIF

#ifdef __synout
   WRITE(kout,*)' '
   WRITE(kout,*)' End of grids_writing'
   WRITE(kout,*)' ********************'
   WRITE(kout,*)' '
#endif
  ENDIF ! p_parallel_io

  RETURN

  END SUBROUTINE grids_writing

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  avg_o2a
!
! !INTERFACE:
!
  SUBROUTINE avg_o2a(kout)
!
! !USES:
!
! !RETURN VALUE:
!
  IMPLICIT NONE
  INTEGER, INTENT(in) :: kout
!
! !PARAMETERS:
!
  REAL(pwp)                :: fakt  ! normalizing factor
  INTEGER                  :: jx,jy ! loop indices
!
! !DESCRIPTION:
!
!- average exchange fields at the end of a coupled time step.
!- overlapped cells are not exchanged.
!- transform into appropriate units if needed
!- land points should be set to zero
!- truncate ice thickness to max allowed value 
!  (might be needed due to truncation erros).
!!
! !REVISION HISTORY:
! 03.05.22  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

!     1.   Transform and average the field.
!             -------------------------------

      fakt = 1.0_pwp/o2a_time
#ifdef __synout
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Averaging of exchange fields.'
      WRITE(io_stdout,*) ' With factor ',fakt
      CALL FLUSH(io_stdout)
#endif

      DO jx=1,jpdim_ocei
         DO jy=1,jpdim_ocej
            gl_sitoacc(jx,jy) = gl_sitoacc(jx,jy)*fakt*weto_g(jx+1,jy,1)
            gl_sicoacc(jx,jy) = gl_sicoacc(jx,jy)*fakt*weto_g(jx+1,jy,1)
            gl_sntoacc(jx,jy) = gl_sntoacc(jx,jy)*fakt*weto_g(jx+1,jy,1)
            gl_sstacc (jx,jy) = gl_sstacc(jx,jy)*fakt*weto_g(jx+1,jy,1)+tmelt
            gl_socuacc (jx,jy) = gl_socuacc(jx,jy)*fakt*weto_g(jx+1,jy,1)
            gl_socvacc (jx,jy) = gl_socvacc(jx,jy)*fakt*weto_g(jx+1,jy,1)
#ifdef ADDCONTRA
            gl_o16acc (jx,jy) = gl_o16acc(jx,jy)*fakt*weto_g(jx+1,jy,1)
            gl_o18acc (jx,jy) = gl_o18acc(jx,jy)*fakt*weto_g(jx+1,jy,1)
            gl_hdoacc (jx,jy) = gl_hdoacc(jx,jy)*fakt*weto_g(jx+1,jy,1)

! convert to delta values (H216O is set to zero permill)           
            IF (gl_o16acc (jx,jy) > 0.0) THEN
              gl_o18acc (jx,jy) = (gl_o18acc(jx,jy)/gl_o16acc(jx,jy)/2005.2e-6 - 1.0) * 1000.
              gl_hdoacc (jx,jy) = (gl_hdoacc(jx,jy)/gl_o16acc(jx,jy)/155.76e-6 - 1.0) * 1000.
            ELSE
              gl_o18acc (jx,jy) = 0.0           
              gl_hdoacc (jx,jy) = 0.0           
            ENDIF
            gl_o16acc (jx,jy) = 0.0
#endif

         ENDDO
      ENDDO


!     3.    No computation of effectice ice thickness.
!           Exchange sea ice in ice-meter
!                    snow in fresh water equivalent.
!           ----------------------------------------

      DO jx=1,jpdim_ocei
         DO JY=1,jpdim_ocej
            IF (gl_sicoacc(jx,jy) > 0.0 ) THEN
                gl_sitoacc(jx,jy) = MIN(1000.0,(gl_sitoacc(jx,jy)/gl_sicoacc(jx,jy)))
                gl_sntoacc(jx,jy) = MIN(1000.0,(gl_sntoacc(jx,jy)/gl_sicoacc(jx,jy)))
                gl_sntoacc(jx,jy) = gl_sntoacc(jx,jy)*rhosno/rhowat
            ELSE
               gl_sitoacc(jx,jy) = 0.0
               gl_sntoacc(jx,jy) = 0.0
               gl_sicoacc(jx,jy) = 0.0
            ENDIF
         ENDDO
      ENDDO


  END SUBROUTINE avg_o2a

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  acc_o2a
!
! !INTERFACE:
!
  SUBROUTINE acc_o2a(kout)
!
! !USES:
!
! !RETURN VALUE:
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: kout
!
! !PARAMETERS:
!
  INTEGER             :: jx,jy,jym1 ! loop indices
!
! !DESCRIPTION:
!
!- accumulate ocean SST, ice thickness and snow depth .
!                   tho, sictho, sicomo, sicsno
!!
! !REVISION HISTORY:
! 03.05.22  S. Legutke - created
!
! EOP
!----------------------------------------------------------------------

! Gather fields for accumulation
  call gather(tho(:,:,1),gl_sst,p_io)
  call gather(sictho,gl_sictho,p_io)
  call gather(sicomo,gl_sicomo,p_io)
  call gather(sicsno,gl_sicsno,p_io)
  call gather(uko(:,:,1),gl_socu,p_io)
  call gather(vke(:,:,1),gl_socv,p_io)
  call gather(sicuo(:,:),gl_sicu,p_io)
  call gather(sicve(:,:),gl_sicv,p_io)
#ifdef ADDCONTRA
  call gather(ocectra(:,:,1,ih2o16),gl_o16,p_io)
  call gather(ocectra(:,:,1,ih2o18),gl_o18,p_io)
  call gather(ocectra(:,:,1,ihDo16),gl_hdo,p_io)
#endif

  if (p_parallel_io) then
     DO jy = 1,jpdim_ocej
        DO jx = 1,jpdim_ocei
           gl_sitoacc(jx,jy) =  gl_sitoacc(jx,jy) + gl_sictho(jx+1,jy)*dt
           gl_sicoacc(jx,jy) =  gl_sicoacc(jx,jy) + gl_sicomo(jx+1,jy)*dt
           gl_sntoacc(jx,jy) =  gl_sntoacc(jx,jy) + gl_sicsno(jx+1,jy)*dt
           gl_sstacc (jx,jy) =  gl_sstacc (jx,jy) + gl_sst   (jx+1,jy)*dt
           jym1=jy-1
           IF (jym1 == 0) jym1=1
           gl_socuacc (jx,jy) =  gl_socuacc (jx,jy) + weto_g(jx+1,jy,1)*           &
              ((1.-gl_sicomo(jx+1,jy))*(gl_socu(jx+1,jy)+gl_socu(jx,jy))           &
            + gl_sicomo(jx+1,jy)*(gl_sicu(jx+1,jy)+gl_sicu(jx,jy)))*dt*0.5
           gl_socvacc (jx,jy) =  gl_socvacc (jx,jy) + weto_g(jx+1,jy,1)*           &
             ((1.-gl_sicomo(jx+1,jy))*(gl_socv(jx+1,jy)+gl_socv(jx+1,jym1))        &
            + gl_sicomo(jx+1,jy)*(gl_sicv(jx+1,jy)+gl_sicv(jx+1,jym1)))*dt*0.5
#ifdef ADDCONTRA
           gl_o16acc (jx,jy) =  gl_o16acc (jx,jy) + gl_o16(jx+1,jy)*dt
           gl_o18acc (jx,jy) =  gl_o18acc (jx,jy) + gl_o18(jx+1,jy)*dt
           gl_hdoacc (jx,jy) =  gl_hdoacc (jx,jy) + gl_hdo(jx+1,jy)*dt
#endif
            ENDDO
         ENDDO
      o2a_time = o2a_time +  dt
#ifdef __synout
      WRITE(kout,*) ' Accumulating since ',o2a_time
      WRITE(kout,*)' gl_sitoacc: '                   &
     &             ,(gl_sitoacc(jpdim_ocei/2+1,jy),jy=3,jpdim_ocej,15)
      WRITE(kout,*)' gl_sicoacc: '                   &
     &             ,(gl_sicoacc(jpdim_ocei/2+1,jy),jy=3,jpdim_ocej,15)
      WRITE(kout,*)' gl_sntoacc: '                   &
     &             ,(gl_sntoacc(jpdim_ocei/2+1,jy),jy=3,jpdim_ocej,15)
      WRITE(kout,*)' gl_sstacc: '                   &
     &             ,(gl_sstacc(jpdim_ocei/2+1,jy),jy=3,jpdim_ocej,15)
#ifdef ADDCONTRA
      WRITE(kout,*)' gl_o16acc: '                   &
     &             ,(gl_o16acc(jpdim_ocei/2+1,jy),jy=3,jpdim_ocej,15)
      WRITE(kout,*)' gl_o18acc: '                   &
     &             ,(gl_o18acc(jpdim_ocei/2+1,jy),jy=3,jpdim_ocej,15)
      WRITE(kout,*)' gl_hdoacc: '                   &
     &             ,(gl_hdoacc(jpdim_ocei/2+1,jy),jy=3,jpdim_ocej,15)
#endif
      CALL FLUSH(kout)
#endif
   endif ! p_parallel_io

! Brodcast ocean time to all proccessors
   call p_bcast(o2a_time,p_io)

 END SUBROUTINE acc_o2a

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  ini_o2a
!
! !INTERFACE:
!
SUBROUTINE ini_o2a(kout)
!
! !USES:
!
! !RETURN VALUE:
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: kout
!
! !PARAMETERS:
!
  INTEGER             :: jx,jy ! loop indices
!
! !DESCRIPTION:
!
!- initialise arrays for  accumulating ocean SST (and ice/snow 
!  thickness) at the beginning of each coupled time step.
!!
! !REVISION HISTORY:
! 03.05.22  S. Legutke - created
!
! EOP
!----------------------------------------------------------------------

  IF (p_parallel_io) THEN
     DO jy = 1,jpdim_ocej
        DO jx = 1,jpdim_ocei
           gl_sitoacc(jx,jy) =  0.0
           gl_sicoacc(jx,jy) =  0.0
           gl_sntoacc(jx,jy) =  0.0
           gl_sstacc (jx,jy) =  0.0
           gl_socuacc (jx,jy) = 0.0
           gl_socvacc (jx,jy) = 0.0
#ifdef ADDCONTRA
           gl_o16acc (jx,jy) =  0.0
           gl_o18acc (jx,jy) =  0.0
           gl_hdoacc (jx,jy) =  0.0
#endif
        ENDDO
     ENDDO
  ENDIF



  o2a_time = 0.0_pwp

#ifdef __synout
  IF (p_parallel_io) THEN
     WRITE(kout,*) ' Accumulation time and exchange fields initialized.'
     CALL FLUSH(kout)
  END IF
#endif

END SUBROUTINE ini_o2a

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  ini_wrte
!
! !INTERFACE:
!
  SUBROUTINE ini_wrte(kunit,kout)
!
! !USES:
!
    USE MO_PARAM1
    USE MO_COMMO1
    USE MO_COMMOAU1
    USE MO_THACC
    !
    ! !RETURN VALUE:
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(INOUT) :: kout
    INTEGER, INTENT(IN)    :: kunit
    !
    ! !PARAMETERS:
    !
    INTEGER      :: jx,jy,jym1                    ! loop indices
    REAL(pwp)    :: exfld(jpdim_ocei,jpdim_ocej)
    
    !
    ! !DESCRIPTION:
    ! - If no restart files for the interface OASIS are
    ! available from a preceeding run, i.e. the run is the first of the 
    ! experiment (nmseq>1), at least one model needs data to start
    ! the run and provide the data for the other models.
    ! Here the ocean surface conditions file  is created using the actual
    ! SST and ice conditions of the ocean restart file or of the initial 
    ! ocean conditions.
    !- change unit
    !
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    !
    ! EOP
    !----------------------------------------------------------------------

    !
    !--   Celsius --> Kelvin
    !     Compute effective ice thickness and initialize
    !
    ! Gather fields for accumulation
    CALL gather(tho(:,:,1),gl_sst,p_io)
    CALL gather(sictho,gl_sictho,p_io)
    CALL gather(sicomo,gl_sicomo,p_io)
    CALL gather(sicsno,gl_sicsno,p_io)
    CALL gather(uko(:,:,1),gl_socu,p_io)
    CALL gather(vke(:,:,1),gl_socv,p_io)
    CALL gather(sicuo(:,:),gl_sicu,p_io)
    CALL gather(sicve(:,:),gl_sicv,p_io)
#ifdef ADDCONTRA
    CALL gather(ocectra(:,:,1,ih2o16),gl_o16,p_io)
    CALL gather(ocectra(:,:,1,ih2o18),gl_o18,p_io)
    CALL gather(ocectra(:,:,1,ihDo16),gl_hdo,p_io)
#endif

    IF (p_parallel_io) THEN
       WRITE(kout,*)'processing global fields'         
       DO jx=1,jpdim_ocei
          DO jy=1,jpdim_ocej
             IF (gl_sicomo(jx+1,jy) > 0.0) THEN
                gl_sicoacc(jx,jy)=gl_sicomo(jx+1,jy)
                gl_sntoacc(jx,jy)=MIN(1000.0,gl_sicsno(jx+1,jy)/gl_sicomo(jx+1,jy))
                gl_sitoacc(jx,jy)=MIN(1000.0,gl_sictho(jx+1,jy)/gl_sicomo(jx+1,jy))
             ELSE
                gl_sitoacc(jx,jy) =  0.0
                gl_sicoacc(jx,jy) =  0.0
                gl_sntoacc(jx,jy) =  0.0
             ENDIF
             gl_sstacc (jx,jy) = gl_sst(jx+1,jy)*weto_g(jx+1,jy,1)+tmelt
             gl_sicoacc(jx,jy) = gl_sicoacc(jx,jy)*weto_g(jx+1,jy,1)
             gl_sntoacc(jx,jy) = gl_sntoacc(jx,jy)*weto_g(jx+1,jy,1)*rhosno/rhowat
             gl_sitoacc(jx,jy) = gl_sitoacc(jx,jy)*weto_g(jx+1,jy,1)
             gl_socuacc(jx,jy) = weto_g(jx+1,jy,1)*                                 &
                  ((1.-gl_sicomo(jx+1,jy))*(gl_socu(jx+1,jy)+gl_socu(jx,jy))           &
                  + gl_sicomo(jx+1,jy)*(gl_sicu(jx+1,jy)+gl_sicu(jx,jy)))*0.5
             jym1=jy-1
             IF ( jym1 == 0) jym1 = 1
             gl_socvacc (jx,jy) =  weto_g(jx+1,jy,1)*                                &
                  ((1.-gl_sicomo(jx+1,jy))*(gl_socv(jx+1,jy)+gl_socv(jx+1,jym1))        &
                  + gl_sicomo(jx+1,jy)*(gl_sicv(jx+1,jy)+gl_sicv(jx+1,jym1)))*0.5
#ifdef ADDCONTRA
             gl_o16acc (jx,jy) = gl_o16(jx+1,jy)*weto_g(jx+1,jy,1)
             gl_o18acc (jx,jy) = gl_o18(jx+1,jy)*weto_g(jx+1,jy,1)
             gl_hdoacc (jx,jy) = gl_hdo(jx+1,jy)*weto_g(jx+1,jy,1)

! convert to delta values (H216O is set to zero permill)           
            IF (gl_o16acc (jx,jy) > 0.0) THEN
              gl_o18acc (jx,jy) = (gl_o18acc(jx,jy)/gl_o16acc(jx,jy)/2005.2e-6 - 1.0) * 1000.
              gl_hdoacc (jx,jy) = (gl_hdoacc(jx,jy)/gl_o16acc(jx,jy)/155.76e-6 - 1.0) * 1000.
            ELSE
              gl_o18acc (jx,jy) = 0.0           
              gl_hdoacc (jx,jy) = 0.0           
            ENDIF
            gl_o16acc (jx,jy) = 0.0
#endif
          ENDDO
       ENDDO

       CALL rotate_2_north_east(gl_socuacc,gl_socvacc,jpdim_ocei,jpdim_ocej)

!   WRITE(io_stdout,*) 'WISO-TEST MPIOM PUT_O2A 1861 GL_O16ACC', SUM(gl_o16acc)
!   WRITE(io_stdout,*) 'WISO-TEST MPIOM PUT_O2A 1861 GL_O18ACC', SUM(gl_o18acc)
!   WRITE(io_stdout,*) 'WISO-TEST MPIOM PUT_O2A 1861 GL_HDOACC', SUM(gl_hdoacc)

       !     
       !*      2.     Write exchange fields.
       !              ----------------------
       !     
       !*    1.0     Open OASIS data file.
       !             ---------------------
#ifdef __synout
       WRITE(kout,*)'Opening file ',cnfileow,' UNIT=',kunit
#endif
       OPEN (kunit,STATUS='UNKNOWN',FILE=cnfileow,FORM='UNFORMATTED')
    END IF
    !     
    !*    2.0    Write exchange fields to OASIS data file.
    !            -----------------------------------------
    DO jfld = 1,nfldo2a
       CALL GET_O2A(clstr8o2a(jfld),exfld,jpdim_ocei,jpdim_ocej,kout)
       IF (p_parallel_io) THEN
          WRITE(kunit) clstr8o2a(jfld)
          WRITE(kunit) exfld
#ifdef __synout
          WRITE(kout,*)'ini_wrte : ',clstr8o2a(jfld),             &
               &   (exfld(jpdim_ocei/2+1,jy),jy=3,jpdim_ocej,15)
#endif
       END IF
    ENDDO

    IF (p_parallel_io) THEN
       CALL FLUSH(kunit)
       CLOSE(kunit)

#ifdef __synout
       WRITE(kout,*)
       WRITE(kout,*)' End of ini_wrte'
       WRITE(kout,*)' ***************'
       CALL FLUSH(kout)
#endif

    ENDIF ! p_parallel_io

  END SUBROUTINE ini_wrte

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  get_o2a
!
! !INTERFACE:
!
  SUBROUTINE get_o2a(clfield,pfield,jpdimi,jpdimj,kout)
!
! !USES:
!
  USE MO_PARAM1
  USE MO_THACC
  USE MO_COMMO1


! !RETURN VALUE:
!
  IMPLICIT NONE
!

  INTEGER, INTENT(IN) :: jpdimi,jpdimj         ! field dimensions
  INTEGER, INTENT(IN) :: kout                  ! unit for std output

  CHARACTER(len=8), INTENT(IN) :: clfield      ! fields symbolic name

  REAL, INTENT(out)   :: pfield(jpdimi,jpdimj) ! target field


 
!
! !PARAMETERS:
!
  INTEGER :: index                 ! field search index
  INTEGER :: jx,jy                 ! loop indices
!
! !DESCRIPTION:
! - moves the model field to be exchange to the array sent to OASIS. 
!
! !REVISION HISTORY:
! 03.05.22  S. Legutke - created
!
! EOP
!----------------------------------------------------------------------
!       OPEN(497,FILE='THOOC.ext',FORM='UNFORMATTED',POSITION='APPEND')
!       OPEN(498,FILE='POCUO.ext',FORM='UNFORMATTED',POSITION='APPEND')
!       OPEN(499,FILE='POCVO.ext',FORM='UNFORMATTED',POSITION='APPEND')

      index = 0

      index = index + 1
      IF ( clfield == clstr8o2a(index) ) THEN
       IF (p_parallel_io) THEN
         DO jx=1,jpdimi
            DO jy=1,jpdimj
               pfield (jx,jy) = gl_sstacc (jx,jy)
            ENDDO
         ENDDO
#ifdef __synout
       WRITE(kout,*)' get_o2a: ',clstr8o2a(index)         &
     &              ,(pfield(jpdimi/2+1,jy),jy=3,jpdimj,15)
#endif /*__synout*/
       END IF
      GOTO 1000
      ENDIF

      index = index + 1
      IF ( clfield == clstr8o2a(index) ) THEN
       IF (p_parallel_io) THEN
      DO jx=1,jpdimi
         DO jy=1,jpdimj
             pfield (jx,jy) = gl_sitoacc(jx,jy)
         ENDDO
      ENDDO
#ifdef __synout
      WRITE(kout,*)' get_o2a: ',clstr8o2a(index)                    &
     &             ,(pfield(jpdimi/2+1,jy),jy=3,jpdimj,15)
#endif /*__synout*/
      ENDIF
      GOTO 1000
      ENDIF

      index = index + 1
      IF ( clfield == clstr8o2a(index) ) THEN
       IF (p_parallel_io) THEN
      DO jx=1,jpdimi
         DO jy=1,jpdimj
             pfield (jx,jy) = gl_sicoacc(jx,jy)
         ENDDO
      ENDDO
#ifdef __synout
      WRITE(kout,*)' get_o2a: ',clstr8o2a(index)                    &
     &             ,(pfield(jpdimi/2+1,jy),jy=3,jpdimj,15)
#endif /*__synout*/
      ENDIF
      GOTO 1000
      ENDIF

      index = index + 1
      IF ( clfield == clstr8o2a(index) ) THEN
       IF (p_parallel_io) THEN
      DO jx=1,jpdimi
         DO jy=1,jpdimj
             pfield (jx,jy) = gl_sntoacc(jx,jy)
         ENDDO
      ENDDO
#ifdef __synout
      WRITE(kout,*)' get_o2a: ',clstr8o2a(index)                    &
     &             ,(pfield(jpdimi/2+1,jy),jy=3,jpdimj,15)
#endif /*__synout*/
      ENDIF
      GOTO 1000
      ENDIF

      index = index + 1
      IF ( clfield == clstr8o2a(index) ) THEN
       IF (p_parallel_io) THEN
      DO jx=1,jpdimi
         DO jy=1,jpdimj
             pfield (jx,jy) = gl_socuacc(jx,jy)
         ENDDO
      ENDDO
#ifdef __synout
      WRITE(kout,*)' get_o2a: ',clstr8o2a(index)                    &
     &             ,(pfield(jpdimi/2+1,jy),jy=3,jpdimj,15)
#endif /*__synout*/
      ENDIF
      GOTO 1000
      ENDIF

      index = index + 1
      IF ( clfield == clstr8o2a(index) ) THEN
       IF (p_parallel_io) THEN
      DO jx=1,jpdimi
         DO jy=1,jpdimj
             pfield (jx,jy) = gl_socvacc(jx,jy)
         ENDDO
      ENDDO
#ifdef __synout
      WRITE(kout,*)' get_o2a: ',clstr8o2a(index)                    &
     &             ,(pfield(jpdimi/2+1,jy),jy=3,jpdimj,15)
#endif /*__synout*/
      ENDIF
      GOTO 1000
      ENDIF

#ifdef ADDCONTRA

! here:  pass isotope concentration of first ocean layer (H216O, H218O, HDO) to atmosphere model
!       (fix tracer order is 1=H216O, 2=H218O, 3=HDO)
 

      index = index + 1
      IF ( clfield == clstr8o2a(index) ) THEN
       IF (p_parallel_io) THEN
      DO jx=1,jpdimi
         DO jy=1,jpdimj
             pfield (jx,jy) = gl_o16acc(jx,jy)
         ENDDO
      ENDDO
#ifdef __synout
      WRITE(kout,*)' get_o2a: ',clstr8o2a(index)                    &
     &             ,(pfield(jpdimi/2+1,jy),jy=3,jpdimj,15)
#endif /*__synout*/
      ENDIF
      GOTO 1000
      ENDIF

      index = index + 1
      IF ( clfield == clstr8o2a(index) ) THEN
       IF (p_parallel_io) THEN
      DO jx=1,jpdimi
         DO jy=1,jpdimj
             pfield (jx,jy) = gl_o18acc(jx,jy)
         ENDDO
      ENDDO
#ifdef __synout
      WRITE(kout,*)' get_o2a: ',clstr8o2a(index)                    &
     &             ,(pfield(jpdimi/2+1,jy),jy=3,jpdimj,15)
#endif /*__synout*/
      ENDIF
      GOTO 1000
      ENDIF

      index = index + 1
      IF ( clfield == clstr8o2a(index) ) THEN
       IF (p_parallel_io) THEN
      DO jx=1,jpdimi
         DO jy=1,jpdimj
             pfield (jx,jy) = gl_hdoacc(jx,jy)
         ENDDO
      ENDDO
#ifdef __synout
      WRITE(kout,*)' get_o2a: ',clstr8o2a(index)                    &
     &             ,(pfield(jpdimi/2+1,jy),jy=3,jpdimj,15)
#endif /*__synout*/
      ENDIF
      GOTO 1000
      ENDIF

#endif

#ifdef __cpl_dms

! ATTN: dmsflux is only needed for diagnostics; for performace reasons it should be removed from the coupling interface
  
      index = index + 1
      IF ( clfield == clstr8o2a(index) ) THEN

         CALL gather(dmsflux,gl_dmsflux,p_io)

         IF (p_parallel_io) THEN
            DO jx=1,jpdimi
               DO jy=1,jpdimj
                  pfield (jx,jy) = gl_dmsflux(jx+1,jy)*32./dt         ! [kmol/[m**2 timestep]] -->  [kg(S)/[m**2 s]]
               ENDDO
            ENDDO

#ifdef __synout
            WRITE(kout,*)' get_o2a: ',clstr8o2a(index)                    &
                 &             ,(pfield(jpdimi/2+1,jy),jy=3,jpdimj,15)
#endif /*__synout*/
         ENDIF
         GOTO 1000
      ENDIF

      index = index + 1
      IF ( clfield == clstr8o2a(index) ) THEN

         CALL gather(dms,gl_dms,p_io)

         IF (p_parallel_io) THEN
            DO jx=1,jpdimi
               DO jy=1,jpdimj
                  pfield (jx,jy) = gl_dms(jx+1,jy)*1.e9              ![kmol/m**3 --> nanomol/l]
               ENDDO
            ENDDO
#ifdef __synout
            WRITE(kout,*)' get_o2a: ',clstr8o2a(index)         &
                 &                     ,(pfield(jpdimi/2+1,jy),jy=3,jpdimj,15)
#endif /*__synout*/
         ENDIF
         GOTO 1000
      ENDIF
#endif /*__cpl_dms*/

#ifdef __cpl_co2
      index = index + 1
      IF ( clfield == clstr8o2a(index) ) THEN

         CALL gather(co2trans,gl_co2trans,p_io)

         IF (p_parallel_io) THEN
            DO jx=1,jpdimi
               DO jy=1,jpdimj
                  pfield (jx,jy) = gl_co2trans(jx+1,jy)
               ENDDO
            ENDDO
#ifdef __synout
            WRITE(kout,*)' get_o2a: ',clstr8o2a(index)         &
                 &                     ,(pfield(jpdimi/2+1,jy),jy=3,jpdimj,15)
#endif /*__synout*/
         END IF
         GOTO 1000
      ENDIF

      index = index + 1
      IF ( clfield == clstr8o2a(index) ) THEN

         CALL gather(suppco2,gl_suppco2,p_io)

         IF (p_parallel_io) THEN
            DO jx=1,jpdimi
               DO jy=1,jpdimj
                  pfield (jx,jy) = gl_suppco2(jx+1,jy)             ![ppm CO2]
               ENDDO
            ENDDO
#ifdef __synout
            WRITE(kout,*)' get_o2a: ',clstr8o2a(index)         &
                 &                     ,(pfield(jpdimi/2+1,jy),jy=3,jpdimj,15)
#endif /*__synout*/
         END IF
         GOTO 1000
      ENDIF
#endif /*__cpl_co2*/

      CALL couple_abort (modnam,'get_o2a','pb Invalid locator string: '//clfield,io_stderr)

 1000 CONTINUE

#ifdef __synout
      CALL FLUSH(kout)
#endif /*__synout*/

      END SUBROUTINE get_o2a


!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE: put_a2o
!
! !INTERFACE:
!
  SUBROUTINE put_a2o(clfield,pfield,jpdimi,jpdimj,kout)
!
! !USES:
!
  USE MO_PARAM1
  USE MO_COMMO1
  USE MO_FLUXES1
! !RETURN VALUE:
!
  IMPLICIT NONE
!
  INTEGER, INTENT(IN) :: jpdimi,jpdimj         ! field dimensions
  INTEGER, INTENT(IN) :: kout                  ! unit for std output

  CHARACTER(len=8), INTENT(IN) :: clfield      ! fields symbolic name

  REAL, INTENT(IN) :: pfield(jpdimi,jpdimj) ! source field

  REAL :: zfield(ie_g,je_g)

!#ifdef __cpl_co2 
!  REAL :: gl_co2conc(ie_g,je_g)
!  REAL :: gl_co2flux_cpl(ie_g,je_g)
!#endif
!#ifdef __cpl_dust 
!  REAL :: gl_dustdep(ie_g,je_g)
!#endif

!
! !PARAMETERS:
!
  INTEGER :: jx,jy                 ! loop indeces
!
! !DESCRIPTION:
! - only the 'inner' oceanic grid region excluding the cyclic boundary 
! columns are exchanged by OASIS.
! put_a2o distributes them to the respective array defined
! for the complete ocean grid.
!
! !REVISION HISTORY:
! 03.05.22  S. Legutke - created
!
! EOP
!----------------------------------------------------------------------

!
!JJ 01.10.99
!JJ ATTN: FOR READING AND WIND STRESS ROTATION TEMP 2D FIELDS
!  ARE USED
!::      V1E===wind stress y water on u-point
!::      Z1E===wind stress y ice on u-point
!::      U1E===wind stress x water on v-point
!::      B1E===wind stress x ice on v-point
!*       2.    Forcing with fluxes only.
!              -------------------------
!
!JJ 27.9.99 INTRODUCE WIND STRESS COMPONENTS ON
!    VELOCITY POINTS, USE SURFACE FIELDS TX,TY AS TEMPORARY
!    VARABLES BEFORE ROTATION INTO NEW GRID ORIENTATION
!*       2.1   Zonal wind stress over water.
!*       2.1a  Zonal wind stress over water on u-point
!              -----------------------------

      IF ( CLFIELD == 'TXWOCEAU' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*AMSUO_G_L1(JX,JY)
            ENDDO
         ENDDO
         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO

       ENDIF
       call scatter_arr(zfield, aofltxwo, p_io)
         GOTO 1000
      ENDIF

!*       2.1b  Zonal wind stress over water on v-point
!              -----------------------------

      IF ( CLFIELD == 'TXWOCEAV' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*AMSUO_G_L1(JX,JY)
            ENDDO
         ENDDO

         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(2,JY)
         ENDDO
       ENDIF
      
       call scatter_arr(zfield, u1e, p_io)
         GOTO 1000
      ENDIF
!
!*       2.2   Meridional wind stress water.
!*       2.2a  Meridional wind stress water on u-point
!              -----------------------------

      IF ( CLFIELD == 'TYWOCEAU' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*AMSUE_G_L1(JX,JY) 
            ENDDO
         ENDDO
         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, v1e, p_io)
         GOTO 1000
      ENDIF

!*       2.2b  Meridional wind stress water on v-point
!              -----------------------------

      IF ( CLFIELD == 'TYWOCEAV' ) THEN

       IF(p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*AMSUE_G_L1(JX,JY) 
            ENDDO
         ENDDO
         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF
 
       call scatter_arr(zfield, AOFLTYWE, p_io)
         GOTO 1000
      ENDIF

!
!*       2.3   Zonal wind stress ice.
!*       2.3a  Zonal wind stress ice on u-point
!              ----------------------

      IF ( CLFIELD == 'TXIOCEAU' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*AMSUO_G_L1(JX,JY)
            ENDDO
         ENDDO

         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, AOFLTXIO, p_io)
         GOTO 1000
      ENDIF

!*       2.3b  Zonal wind stress ice on v-point
!              ----------------------

      IF ( CLFIELD == 'TXIOCEAV' ) THEN

       IF (p_parallel_io) then
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*AMSUO_G_L1(JX,JY)
            ENDDO
         ENDDO

         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, b1e, p_io)

         GOTO 1000
      ENDIF


!*       2.4   Meridional wind stress ice.
!*       2.4a  Meridional wind stress ice on u-point
!              ---------------------------

      IF ( CLFIELD == 'TYIOCEAU' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*AMSUE_G_L1(JX,JY)
            ENDDO
         ENDDO
         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, z1e, p_io)

         GOTO 1000
      ENDIF

!*       2.4b  Meridional wind stress ice on v-point
!              ---------------------------

      IF ( CLFIELD == 'TYIOCEAV' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*AMSUE_G_L1(JX,JY)
            ENDDO
         ENDDO
         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, AOFLTYIE, p_io)
         GOTO 1000

      ENDIF

!*       2.5   Freshwater flux over ice.
!              -------------------------

      IF ( CLFIELD == 'FRIOCEAN' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
            ENDDO
         ENDDO
         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, aoflfrio, p_io)
         GOTO 1000
      ENDIF

!*       2.6   Freshwater flux over water.
!              --------------------------

      IF ( CLFIELD == 'FRWOCEAN' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
            ENDDO
         ENDDO
         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, aoflfrwo, p_io)
         GOTO 1000
      ENDIF

!*       2.7   Residual heat flux over ice.
!              ----------------------------

      IF ( CLFIELD == 'RHIOCEAN' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
            ENDDO
         ENDDO

         DO JY = 1,JE_G   
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, aoflrhio, p_io)

         GOTO 1000
      ENDIF

!*       2.8   Conductive heat flux through ice.
!              ---------------------------------

      IF ( CLFIELD == 'CHIOCEAN' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1,JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
            ENDDO
         ENDDO
         DO JY = 1,JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, aoflchio, p_io)
         GOTO 1000
      ENDIF

!*       2.9   Net heat flux over water.
!              -------------------------

      IF ( CLFIELD == 'NHWOCEAN' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1,JE_G
            DO  JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
            ENDDO   
         ENDDO   
         DO JY = 1,JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO 
       ENDIF
       
        call scatter_arr(zfield, AOFLNHWO, p_io)
         GOTO 1000
      ENDIF

!*       2.10   Solar heat flux.
!               ----------------

      IF ( CLFIELD == 'SHWOCEAN' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1,JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
            ENDDO
         ENDDO
         DO JY = 1,JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO 
       ENDIF
       
       call scatter_arr(zfield, AOFLSHWO, p_io)
         GOTO 1000
      ENDIF

!*       2.11   Wind Stress Velocity.
!               --------------------

      IF ( CLFIELD == 'WSVOCEAN' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1,JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
            ENDDO
         ENDDO
         DO JY = 1,JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, AOFLWSVO, p_io)
         GOTO 1000
      ENDIF

#ifdef ADDCONTRA

!*             Freshwater flux over water - water isotopes
!              -------------------------------------------

      IF ( CLFIELD == 'FRWOCO16' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
            ENDDO
         ENDDO
         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, aoflfrwo16, p_io)
         GOTO 1000
      ENDIF

      IF ( CLFIELD == 'FRWOCO18' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
            ENDDO
         ENDDO
         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, aoflfrwo18, p_io)
         GOTO 1000
      ENDIF

      IF ( CLFIELD == 'FRWOCHDO' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
            ENDDO
         ENDDO
         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, aoflfrwhdo, p_io)
         GOTO 1000
      ENDIF

!*             Freshwater flux over ice - water isotopes
!              -----------------------------------------

      IF ( CLFIELD == 'FRIOCO16' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
            ENDDO
         ENDDO
         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, aoflfrio16, p_io)
         GOTO 1000
      ENDIF

      IF ( CLFIELD == 'FRIOCO18' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
            ENDDO
         ENDDO
         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, aoflfrio18, p_io)
         GOTO 1000
      ENDIF

      IF ( CLFIELD == 'FRIOCHDO' ) THEN

       IF (p_parallel_io) THEN
         DO JY = 1, JE_G
            DO JX = 2,IE_G-1
               zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
            ENDDO
         ENDDO
         DO JY = 1, JE_G
            zfield( 1,JY) = zfield(IE_G-1,JY)
            zfield(IE_G,JY) = zfield(  2,JY)
         ENDDO
       ENDIF

       call scatter_arr(zfield, aoflfrihdo, p_io)
         GOTO 1000
      ENDIF

#endif

#ifdef __cpl_dust

!*       2.12   Dust dep. flux
!               --------------------

      IF ( CLFIELD == 'DEPOCEAN' ) THEN

         IF (p_parallel_io) THEN

            DO JY = 1,JE_G
               DO JX = 2,IE_G-1
                  zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
               ENDDO
            ENDDO

            DO JY = 1,JE_G
               zfield( 1,JY) = zfield(IE_G-1,JY)
               zfield(IE_G,JY) = zfield(  2,JY)
            ENDDO

         ENDIF

         call scatter_arr(zfield, dustdep, p_io)

         GOTO 1000

      ENDIF
#endif

#ifdef __cpl_co2

!*       2.12   CO2 concentration
!               --------------------

      IF ( CLFIELD == 'CO2CONOC' ) THEN

         IF (p_parallel_io) THEN

            DO JY = 1,JE_G
               DO JX = 2,IE_G-1
                  zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
               ENDDO
            ENDDO

            DO JY = 1,JE_G
               zfield( 1,JY) = zfield(IE_G-1,JY)
               zfield(IE_G,JY) = zfield(  2,JY)
            ENDDO

         ENDIF

         call scatter_arr(zfield, co2conc, p_io)

         atm(:,:,iatmco2) = co2conc(:,:)*(28.970/44.011)*1.e6 

         GOTO 1000
         
      ENDIF

!*       2.12   CO2 flux
!               ----------

      IF ( CLFIELD == 'CO2FLXOC' ) THEN

         IF (p_parallel_io) THEN

            DO JY = 1,JE_G
               DO JX = 2,IE_G-1
                  zfield(JX,JY) = pfield(JX-1,JY)*WETO_G(JX,JY,1)
               ENDDO
            ENDDO

            DO JY = 1,JE_G
               zfield( 1,JY) = zfield(IE_G-1,JY)
               zfield(IE_G,JY) = zfield(  2,JY)
            ENDDO

         ENDIF

         call scatter_arr(zfield, co2flux_cpl, p_io)

         GOTO 1000
         
      ENDIF

#endif

      CALL couple_abort (modnam,'put_a2o','pb Invalid locator string: '//CLFIELD,io_stderr)

 1000 CONTINUE

!   if (CLFIELD == 'FRWOCEAN') WRITE(io_stdout,*) 'WISO-TEST MPIOM PUT_A2O 2747 FRWOCEAN', SUM(aoflfrwo)
!   if (CLFIELD == 'FRWOCO16') WRITE(io_stdout,*) 'WISO-TEST MPIOM PUT_A2O 2747 FRWOCO16', SUM(aoflfrwo16)
!   if (CLFIELD == 'FRWOCO18') WRITE(io_stdout,*) 'WISO-TEST MPIOM PUT_A2O 2747 FRWOCO18', SUM(aoflfrwo18)
!   if (CLFIELD == 'FRWOCHDO') WRITE(io_stdout,*) 'WISO-TEST MPIOM PUT_A2O 2747 FRWOCHDO', SUM(aoflfrwhdo)


#ifdef __synout
      IF (p_parallel_io) THEN
      WRITE(KOUT,*) ' put_a2o : ',CLFIELD                           &
     &                            ,(pfield((ie_g-2)/2,jy),jy=2,je_g,14)
      ENDIF

#endif

 6001 FORMAT(1x,4x,6(f6.3,3x))
 6002 FORMAT(1x,6(f6.3,3x))

      END SUBROUTINE put_a2o

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE: post_a2o
!
! !INTERFACE:
!
  SUBROUTINE post_a2o(k_fld_recvd, kout, kdtrun)
!
! !USES:
  USE mo_fluxes1, only : wrte_flux_extra 
!
! !RETURN VALUE:
!
  IMPLICIT NONE
!
  INTEGER, INTENT(inout) :: k_fld_recvd    ! counter of received fields
  INTEGER, INTENT(in)    :: kout
  INTEGER, INTENT(in)   ::  kdtrun         ! model time step 
!
! !PARAMETERS:
!
  INTEGER               ::  ncorrect
  INTEGER               ::  jx,jy          ! loop indices
!
! !DESCRIPTION:
! - 
!
! !REVISION HISTORY:
! 03.05.22  S. Legutke - created
!
! EOP
!----------------------------------------------------------------------

#ifdef __synout
     IF (p_parallel_io) THEN
      WRITE(kout,*)
      WRITE(kout,*)' Start of  post_a2o'
      WRITE(kout,*)' ******************'
      WRITE(kout,*)
      CALL FLUSH(kout)
     ENDIF
#endif
      !
      !-- Rotate received wind stress vector (JJ 01.10.99)
      !

         CALL rotate_u(aofltxwo,v1e,ie,je)
         CALL rotate_u(aofltxio,z1e,ie,je)
         CALL rotate_v(u1e,aofltywe,ie,je)
         CALL rotate_v(b1e,aofltyie,ie,je)

         IF(icycli >= 1) THEN
            CALL bounds_exch('u',aofltxwo)
            CALL bounds_exch('u',aofltxio)
            CALL bounds_exch('v',aofltywe)
            CALL bounds_exch('v',aofltyie)
         ENDIF

         !-- Divide windstress by ref. density (JJ 01.10.99)

#ifdef FLUXCORRECT
         WRITE(kout,*) 'Adding flux correction! '
         fhflmax=-999.
         fhflmin=999.                  
#endif /*FLUXCORRECT*/


         DO jx=1,ie
            DO jy=1,je
               aofltxwo(jx,jy)=aofltxwo(jx,jy)/rhowat
               aofltywe(jx,jy)=aofltywe(jx,jy)/rhowat
               aofltyie(jx,jy)=aofltyie(jx,jy)/rhowat
               aofltxio(jx,jy)=aofltxio(jx,jy)/rhowat
#ifdef FLUXCORRECT

#ifdef TEMPCORRECT
               aoflnhwo(jx,jy)=aoflnhwo(jx,jy)+fluko_hfl(jx,jy)
               aoflchio(jx,jy)=aoflchio(jx,jy)+fluko_hfl(jx,jy)
               IF (fhflmax < fluko_hfl(jx,jy)) fhflmax=fluko_hfl(jx,jy)
               IF (fhflmin > fluko_hfl(jx,jy)) fhflmin=fluko_hfl(jx,jy)
#endif /*TEMPCORRECT*/

#ifdef SALTCORRECT
               aoflfrwo(jx,jy)=aoflfrwo(jx,jy)+fluko_frw(jx,jy)
#endif /*SALTCORRECT*/

#endif /*FLUXCORRECT*/

            ENDDO
         ENDDO

#ifdef FLUXCORRECT
         WRITE(kout,*) 'max hfl fluxco ',fhflmax
         WRITE(kout,*) 'min hfl fluxco ',fhflmin
#endif /*FLUXCORRECT*/

         !-- Modify residual/conductive heat flux if negative

         ncorrect = 0
         DO jy = 1,je
            DO jx = 1,ie
               IF (aoflrhio(jx,jy) < 0.0 .AND. weto(jx,jy,1) > 0.5) THEN
                  IF ( ncorrect <= 10 ) THEN
                     IF ( ncorrect == 0 ) WRITE(kout,*)          &
     &                   'Model time step of run : ',kdtrun
                     WRITE(kout,*) 'Modify residual heat flux at :'&
     &                   ,jx,jy,aoflrhio(jx,jy),aoflchio(jx,jy)
                     ncorrect = ncorrect + 1
                  ENDIF
                  aoflchio(jx,jy) = aoflchio(jx,jy) + aoflrhio(jx,jy)
                  aoflrhio(jx,jy) = 0.0
               ENDIF
            ENDDO
         ENDDO

         !
         !-- Move wind stress velocity array used in OCTHER
         !
         DO jy=1,je
            DO jx=1,ie
               fu10(jx,jy)=aoflwsvo(jx,jy)
            ENDDO
         ENDDO

         k_fld_recvd = 0
!
#ifdef CMIP_READ_FLUX
! OVERREAD FLUXES BY STORED FIELDS
         CALL READ_FLUX_EXTRA
#endif
         if ( IOASISFLUX.eq.99 ) then
            CALL WRTE_FLUX_EXTRA
         endif
 
#ifdef __synout
      WRITE(kout,*)
      WRITE(kout,*)' End of  post_a2o'
      WRITE(kout,*)' ****************'
      CALL FLUSH(kout)
#endif
      END SUBROUTINE post_a2o

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE: prep_o2a
!
! !INTERFACE:
!
  SUBROUTINE prep_o2a(kout)
!
! !USES:
!
!
! !RETURN VALUE:
!
  IMPLICIT NONE
!
  INTEGER, INTENT(in)    :: kout
!
! !PARAMETERS:
!
  INTEGER                :: jx,jy,jym1 ! loop indices
!
! !DESCRIPTION:
! - prepare fields before sending to psmile.
!
! !REVISION HISTORY:
! 03.05.22  S. Legutke - created
!
! EOP
!----------------------------------------------------------------------
      !
      !-- Masking and unit transformation.
      !

! Gather arrays for accumulation
  call gather(tho(:,:,1),gl_sst,p_io)
  call gather(sictho,gl_sictho,p_io)
  call gather(sicomo,gl_sicomo,p_io)
  call gather(sicsno,gl_sicsno,p_io)
  call gather(uko(:,:,1),gl_socu,p_io)
  call gather(vke(:,:,1),gl_socv,p_io)
  call gather(sicuo(:,:),gl_sicu,p_io)
  call gather(sicve(:,:),gl_sicv,p_io)
#ifdef ADDCONTRA
  call gather(ocectra(:,:,1,ih2o16),gl_o16,p_io)
  call gather(ocectra(:,:,1,ih2o18),gl_o18,p_io)
  call gather(ocectra(:,:,1,ihDo16),gl_hdo,p_io)
#endif  
      if (p_parallel_io) then
      DO jx=1,jpdim_ocei
         DO jy=1,jpdim_ocej

            gl_sitoacc(jx,jy) = gl_sictho(jx+1,jy)*weto_g(jx+1,jy,1)
            gl_sicoacc(jx,jy) = gl_sicomo(jx+1,jy)*weto_g(jx+1,jy,1)
            gl_sntoacc(jx,jy) = gl_sicsno(jx+1,jy)*weto_g(jx+1,jy,1)

            gl_sstacc (jx,jy) = gl_sst(jx+1,jy)*weto_g(jx+1,jy,1)+tmelt

            gl_socuacc(jx,jy) = weto_g(jx+1,jy,1)*                                 &
              ((1.-gl_sicomo(jx+1,jy))*(gl_socu(jx+1,jy)+gl_socu(jx,jy))           &
              + gl_sicomo(jx+1,jy)*(gl_sicu(jx+1,jy)+gl_sicu(jx,jy)))*0.5
           jym1=jy-1
           IF (jym1 == 0) jym1=1
           gl_socvacc (jx,jy) =  weto_g(jx+1,jy,1)*                                &
             ((1.-gl_sicomo(jx+1,jy))*(gl_socv(jx+1,jy)+gl_socv(jx+1,jym1))        &
             + gl_sicomo(jx+1,jy)*(gl_sicv(jx+1,jy)+gl_sicv(jx+1,jym1)))*0.5
#ifdef ADDCONTRA
            gl_o16acc (jx,jy) = gl_o16(jx+1,jy)*weto_g(jx+1,jy,1)
            gl_o18acc (jx,jy) = gl_o18(jx+1,jy)*weto_g(jx+1,jy,1)
            gl_hdoacc (jx,jy) = gl_hdo(jx+1,jy)*weto_g(jx+1,jy,1)
            
! convert to delta values (H216O is set to zero permill)           
            IF (gl_o16acc (jx,jy) > 0.0) THEN
              gl_o18acc (jx,jy) = (gl_o18acc(jx,jy)/gl_o16acc(jx,jy)/2005.2e-6 - 1.0) * 1000.
              gl_hdoacc (jx,jy) = (gl_hdoacc(jx,jy)/gl_o16acc(jx,jy)/155.76e-6 - 1.0) * 1000.
            ELSE
              gl_o18acc (jx,jy) = 0.0           
              gl_hdoacc (jx,jy) = 0.0           
            ENDIF
            gl_o16acc (jx,jy) = 0.0
#endif
         ENDDO
      ENDDO

      CALL rotate_2_north_east(gl_socuacc,gl_socvacc,jpdim_ocei,jpdim_ocej)

      DO jx=1,jpdim_ocei
         DO JY=1,jpdim_ocej
            IF (gl_sicoacc(jx,jy) > 0.0 ) THEN
                gl_sitoacc(jx,jy) = MIN(1000.0,(gl_sitoacc(jx,jy)/gl_sicoacc(jx,jy)))
                gl_sntoacc(jx,jy) = MIN(1000.0,(gl_sntoacc(jx,jy)/gl_sicoacc(jx,jy)))
                gl_sntoacc(jx,jy) = gl_sntoacc(jx,jy)*rhosno/rhowat
            ELSE
               gl_sitoacc(jx,jy) = 0.0
               gl_sntoacc(jx,jy) = 0.0
               gl_sicoacc(jx,jy) = 0.0
            ENDIF
         ENDDO
      ENDDO
 
#ifdef __synout
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Preparation of fields for sending done.'
      CALL FLUSH(io_stdout)
#endif

      endif ! p_parallel_io

      END SUBROUTINE prep_o2a


!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE: couple_calendar
!
! !INTERFACE:
!
  SUBROUTINE couple_calendar(pdt,kout)
!
! !USES:
!
!
! !RETURN VALUE:
!
  IMPLICIT NONE
!
  INTEGER, INTENT(in)   :: kout  ! unit std output
  REAL,    INTENT(in)   :: pdt   ! model time step length (secs)
!
! !PARAMETERS:
!
!
! !DESCRIPTION:
! - time control of coupling aspects
!
! !REVISION HISTORY:
! 03.05.22  S. Legutke - created
!
! EOP
!----------------------------------------------------------------------

      run_date_secs = run_date_secs + INT(pdt)

#ifdef __synout
      WRITE(kout,*) ' '
      WRITE(kout,*) ' Calendar is updated.'
      CALL FLUSH(kout)
#endif

      END SUBROUTINE couple_calendar
!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  chck_par
!
! !INTERFACE:

  SUBROUTINE chck_par

! !USES:
!
#ifdef NAG
    USE f90_unix,           ONLY: flush
#endif  

! !RETURN VALUE:

  IMPLICIT NONE

! !DESCRIPTION:
!
! - Checks coupled model control parameter against
!   those of the calling model mpi-om.
!
! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

    !     
    !-- Check whether month of restart file is ok.
    !   ------------------------------------------
      IF(lmont1 /= ig_inidate(2)) THEN
         WRITE(io_stderr,*)'  Restart file has not the right month:'
         WRITE(io_stderr,*)'      from restart file = ',lmont1
         WRITE(io_stderr,*)'      from coupler      = ',ig_inidate(2)
         CALL couple_abort (modnam,'chck_par', &
             'Restart file has not the right month',io_stderr)
      ENDIF
      
      a2o_freq = 0
      DO jfld = 1,nflda2o
         a2o_freq = MAX(a2o_freq,ig_def_freq(portin_id(jfld)))
      END DO
      !
      !-- ... whether model allows for exchange algorithm.
      !  
      DO jfld = 1,nflda2o
         IF(a2o_freq /= ig_def_freq(portin_id(jfld))) THEN
              CALL couple_abort (modnam,'chck_par', &
              'The algorithsm is not allowed ',io_stderr)
         ENDIF
      END DO

      o2a_freq = 0
      DO jfld = 1,nfldo2a
         o2a_freq = MAX(o2a_freq,ig_def_freq(portout_id(jfld)))
      END DO

      !
      !-- ...  whether model allows for exchange algorithm.
      !  
      DO jfld = 1,nfldo2a
         IF(o2a_freq /= ig_def_freq(portout_id(jfld))) THEN
              CALL couple_abort (modnam,'chck_par', &
                  'Algorithm not allowed when averaging in model',io_stderr)
         ENDIF
      END DO


  END SUBROUTINE chck_par

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  digest_get_Id
!
! !INTERFACE:

  SUBROUTINE digest_get_Id(kout,kerr,kinfo,clfield,kdate,kcount)

! !USES:
!

! !RETURN VALUE:

  IMPLICIT NONE

  INTEGER,            INTENT(in)  :: kout    ! unit std output
  INTEGER,            INTENT(in)  :: kerr    ! unit error output
  INTEGER,            INTENT(in)  :: kinfo   ! info ID passed from psmile
  INTEGER,            INTENT(in)  :: kdate   ! date (seconds)
  INTEGER,OPTIONAL,   INTENT(in)  :: kcount  ! exchanges field count

  CHARACTER(len=8),   INTENT(in)  :: clfield ! fields symbolic name

!
! !DESCRIPTION:
!
! - Print info after prism_get. Abort if error.
!
! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

      INTEGER              :: icount = -1

#ifdef __synout
      IF (PRESENT(kcount)) THEN
         icount = kcount
      ENDIF

      WRITE(kout,*)' At date (seconds) ',kdate,clfield
      IF (icount /= -1 ) THEN
         WRITE(kout,*)' (field no. ',kcount,')'
      ENDIF
#endif

      IF ( kinfo /= PRISM_Ok          .AND. &
        kinfo /= PRISM_FromRest    .AND. &
        kinfo /= PRISM_Input       .AND. &
        kinfo /= PRISM_RecvOut     .AND. &
        kinfo /= PRISM_FromRestOut .AND. &
        kinfo /= PRISM_Recvd            )  THEN
        WRITE (kout, *)' Problem with port = ',clfield
        WRITE (kout, *)' Error code number = ',kinfo
        WRITE (kout, *)' Seconds passed    = ',kdate
        CALL couple_abort (modnam,'couple_put','pb prism_get',kerr)
#ifdef __synout
      ELSEIF (kinfo == PRISM_Recvd) THEN
         WRITE(kout,*) ' was received from another model'
         CALL FLUSH(kout)
      ELSEIF (kinfo == PRISM_Ok) THEN
         WRITE(kout,*)' was not received; no error.'
         CALL FLUSH(kout)
      ELSEIF (kinfo == PRISM_FromRest) THEN
         WRITE(kout,*) ' was received from restart file.'
         CALL FLUSH(kout)
      ELSEIF (kinfo == PRISM_Input) THEN
         WRITE(kout,*) ' was received from input file.'
         CALL FLUSH(kout)
      ELSEIF (kinfo == PRISM_RecvOut) THEN
         WRITE(kout,*) ' was received from input file and other model.'
         CALL FLUSH(kout)
      ELSEIF (kinfo == PRISM_FromRestOut) THEN
         WRITE(kout,*) ' was received from input file ' &
                      ,'and written to an output file.'
         CALL FLUSH(kout)
#endif
      ENDIF 

  END SUBROUTINE digest_get_Id

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  digest_put_Id
!
! !INTERFACE:

  SUBROUTINE digest_put_Id(kout,kerr,kinfo,clfield,kdate,kcount)

! !USES:
!

! !RETURN VALUE:

  IMPLICIT NONE

  INTEGER,            INTENT(in)  :: kout    ! unit std output
  INTEGER,            INTENT(in)  :: kerr    ! unit error output
  INTEGER,            INTENT(in)  :: kinfo   ! info ID passed from psmile
  INTEGER,            INTENT(in)  :: kdate   ! date (seconds)
  INTEGER,OPTIONAL,   INTENT(in)  :: kcount  ! exchanges field count

  CHARACTER(len=8), INTENT(in) :: clfield      ! fields symbolic name

!
! !DESCRIPTION:
!
! - Print info after prism_put. Abort if error.
!
! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

      INTEGER              :: icount = -1

#ifdef __synout
      IF (PRESENT(kcount)) THEN
         icount = kcount
      ENDIF

      WRITE(kout,*)' At date (seconds) ',kdate,clfield
      IF (icount /= -1 ) THEN
         WRITE(kout,*)' (field no. ',kcount,')'
      ENDIF
#endif

      IF ( kinfo /= PRISM_Ok        .AND. &
         kinfo /= PRISM_LocTrans  .AND. &
         kinfo /= PRISM_ToRest    .AND. &
         kinfo /= PRISM_Output    .AND. &
         kinfo /= PRISM_SentOut   .AND. &
         kinfo /= PRISM_ToRestOut .AND. &
         kinfo /= PRISM_Sent            )  THEN
         WRITE (kout, *)' Problem with port = ',clfield
         WRITE (kout, *)' Error code number = ',kinfo
         WRITE (kout, *)' Seconds passed    = ',kdate
         CALL couple_abort (modnam,'couple_put','pb prism_put',kerr)
#ifdef __synout
      ELSEIF (kinfo == PRISM_Sent) THEN
         WRITE(kout,*) ' was sent to another model'
         CALL FLUSH(kout)
      ELSEIF (kinfo == PRISM_Ok) THEN
         WRITE(kout,*)' was not sent; no error.'
         CALL FLUSH(kout)
      ELSEIF (kinfo == PRISM_LocTrans) THEN
         WRITE(kout,*) ' was used in local transformation.'
         CALL FLUSH(kout)
      ELSEIF (kinfo == PRISM_ToRest) THEN
         WRITE(kout,*) ' was written to a restart file.'
         CALL FLUSH(kout)
      ELSEIF (kinfo == PRISM_Output) THEN
         WRITE(kout,*) ' was output to a file.'
         CALL FLUSH(kout)
      ELSEIF (kinfo == PRISM_SentOut) THEN
         WRITE(kout,*) ' was sent to another model and output to a file.'
         CALL FLUSH(kout)
      ELSEIF (kinfo == PRISM_ToRestOut) THEN
         WRITE(kout,*) ' was sent to another model and written to a restart file.'
         CALL FLUSH(kout)
#endif
      ENDIF 

  END SUBROUTINE digest_put_Id

#endif /* __coupled || __prism */
END MODULE mo_couple
