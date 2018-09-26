MODULE mo_couple
!
! !DESCRIPTION:
!
!   contains additional code used for coupling with OASIS3
!   use with MPI2 only
!  
!   coupling strategies : 
!     compile flag cpl_mpiom  if used with MPI-OM;
!     compile flag cpl_co2    if used with MPI-OM, HAMOCC and interactive CO2
!     compile flag cpl_hope    if used with HOPE; see DKRZ rep. 18
!     compile flag cpl_fluxes4 if only 'open water' fluxes of heat, 
!                               freshwater, and momentum is provided;
!
!   exchanged fields (see SBR collect):
!     over water:   pawsol  downwelling solar radiation (used for solar 
!                           penetration)
!                   pawsta  ustar**3 (C-grid HOPE only)
!                   pawhea  net heat flux
!                   pawfre  liquid freshwater flux
!                   pawust  zonal wind stress
!                   pawvst  meridional wind stress
!
!     over sea ice: paicon  conductive heat flux
!                   paiqre  residual heat flux (used for surface melt)
!                   paifre  solid freshwater flux
!                   paiust  zonal wind stress
!                   paivst  meridional wind stress
!
!     cpl_mpiom   defined: 
!                   aicon,aiqre not multiplied with sea ice concentration
!                   awhea       not multiplied with open water fraction
!                   awust,awvst not multiplied with open water fraction
!     cpl_hope     defined: 
!                   aicon,aiqre     multiplied with sea ice concentration
!                   awhea           multiplied with open water fraction
!                   awust,awvst     multiplied with open water fraction
!     always     multiplied with sea ice concentration: aifre
!     always not  "       "    "   "       "          : awsol,awsta,aiust,
!                                                       aivst
!                  weighted average of both           : awfre
!
! !USES:
!

#ifdef NAG
         USE f90_unix_io, ONLY: flush
#endif

  USE mo_kind,            ONLY: dp ! working precision (echam5)
  USE mo_control,         ONLY: nlon,ngl
  USE mo_time_control,    ONLY: delta_time,get_time_step,lresume,l_putocean, &
                                lbreak, stop_date, next_date, &
                                write_date, l_putrerun
  USE mo_memory_g3b,      ONLY: awsol, &
                                awhea,awfre,awust,awvst,awsta, &
                                aicon,aiqre,aifre,aiust,aivst, &
                                seaice,tsw,siced,sni,alake,slf,slm, &
                                ocu,ocv,afre_residual
#ifdef __cpl_co2
  USE mo_co2,             ONLY: co2trans, co2ocean, co2atmos, co2flux_cpl
#endif
  USE mo_constants,       ONLY: rhoh2o
  USE mo_mpi,             ONLY: p_bcast, p_io, p_parallel_io, p_parallel, p_all_comm
  USE mo_decomposition,   ONLY: gd => global_decomposition
  USE mo_transpose,       ONLY: gather_gp, scatter_gp     
  USE mo_exception,       ONLY: message, finish

!---wiso-code
  USE mo_memory_wiso,     ONLY: wisoawfre, wisoaifre, wisoafre_residual, wisosw_d
  USE mo_wiso,            ONLY: lwiso
!---wiso-code-end

#if defined (__prism)
  USE mod_kinds_model,    ONLY: ip_realwp_p
  USE mod_prism_proto
  USE mod_comprism_proto
  USE mod_prism_get_proto
  USE mod_prism_put_proto
  USE mod_prism_def_partition_proto, ONLY: prism_def_partition_proto
#endif

  USE mo_io,              ONLY: slm_glob

  IMPLICIT NONE    

  PUBLIC

#if defined (__prism)
  INTEGER, PARAMETER :: pwp = ip_realwp_p   ! working precision (psmile)

  INTEGER, ALLOCATABLE :: paral(:) ! parallel strategy
  INTEGER :: nsegments             ! no. of segments
  INTEGER :: parsize               ! 
  INTEGER :: nbtotproc             ! total no. of procs
  INTEGER :: nbcplproc = 1         ! no. of procs involved in data exchange
  INTEGER :: commlocal             ! local communicator
  INTEGER :: comp_id               ! model id
!!$  INTEGER :: size                  ! size of local communicator
  INTEGER :: rank                  ! rank of processor in local communicator
  INTEGER :: info, ierror          ! message, error codes
  INTEGER :: var_shape(2)          !
  INTEGER :: part_id               ! 
  INTEGER :: var_nodims(2)         !

  CHARACTER (len=6) :: modnam = 'echam5' ! Model name

#if defined __cpl_mpiom
#if defined __cpl_co2
!  INTEGER, PARAMETER :: nflda2o = 17 ! no. of fields passed from atmosphere to ocean
  INTEGER, PARAMETER :: nflda2o = 23 ! no. of fields passed from atmosphere to ocean - water isotopes 
                                     ! 6 additional fields: fresh water o16, o18, hdo, snow over sea ice o16, o18, hdo
!  INTEGER, PARAMETER :: nfldo2a =  8 ! no. of fields passed from ocean to atmosphere
  INTEGER, PARAMETER :: nfldo2a = 11 ! no. of fields passed from ocean to atmosphere - water isotopes 
                                     ! 3 additional fields: surface water o16, o18, hdo
#else
!---wiso-code
!  INTEGER, PARAMETER :: nflda2o = 15 ! no. of fields passed from atmosphere to ocean
  INTEGER, PARAMETER :: nflda2o = 21 ! no. of fields passed from atmosphere to ocean - water isotopes 
                                     ! 6 additional fields: fresh water o16, o18, hdo, snow over sea ice o16, o18, hdo
!  INTEGER, PARAMETER :: nfldo2a =  6 ! no. of fields passed from ocean to atmosphere
  INTEGER, PARAMETER :: nfldo2a =  9 ! no. of fields passed from ocean to atmosphere - water isotopes 
                                     ! 3 additional fields: surface water o16, o18, hdo
!---wiso-code-end
#endif
#endif

  CHARACTER (len=8) :: clstr8a2o(nflda2o)   ! Port names of exchange fields sent
  CHARACTER (len=8) :: clstr8o2a(nfldo2a)   ! Port names of exchange fields received

  INTEGER*4 :: portout_id(nflda2o)    ! Port IDs of exchange fields sent
  INTEGER*4 :: portin_id(nfldo2a)     ! Port IDs of exchange fields received

  REAL(pwp) :: couple_a2o_time = 0.0_pwp ! time passed since last couple_put_a2o
#else
  REAL(dp) :: couple_a2o_time = 0.0_dp ! time passed since last couple_put_a2o
#endif


  INTEGER :: nmseq                       ! Maximum sequential number of exchange fields

  INTEGER :: run_date_secs  = 0          ! time (sec) passed since start of run

  INTEGER :: couple_resum = 0            ! model step when coupled run is resumed

  INTEGER :: jfld                        ! loop index

  REAL(dp), POINTER ::         gl_slf   (:,:), gl_alake(:,:), gl_tsw(:,:), &
                               gl_seaice(:,:), gl_siced(:,:), gl_sni(:,:), &
                               gl_ocu(:,:), gl_ocv(:,:), gl_slm(:,:)

!wiso-code

  REAL(dp), POINTER ::         gl_o16(:,:), gl_o18(:,:), gl_hdo(:,:)

!wiso-code-end

#ifdef __cpl_co2
  REAL(dp), POINTER ::         gl_co2trans(:,:)
  REAL(dp), POINTER ::         gl_co2ocean(:,:)
#endif

!
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
! 03.05.14 Stephanie Legutke, MPI/M&D
!       - created
! 03.05.14 Veronika Gayler,   MPI/M&D
!       -
! 03.07.08 Luis Kornblueh,    MPI
!       - minor modifications
! 03.08.04 Johann Jungclaus MPI
!       - include call for water budget correction
!
! EOP
!-----------------------------------------------------------------------
! 

CHARACTER(len=*), PARAMETER, PRIVATE :: cl_version_string = '$Id: mo_couple.f90,v 1.5 2003/10/07 12:31:22 k204059 Exp $'

CONTAINS

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  couple_init
!
! !INTERFACE:

  SUBROUTINE couple_init(kout,kerr)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kout      ! unit number for standard output
    INTEGER, INTENT(in) :: kerr      ! error out unit

! !DESCRIPTION:
!
! - Initialize coupling.
! - This routine should not be called in uncoupled mode !
!
! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

    INTEGER :: k    ! loop index

#if defined (__prism)

#ifdef __synout
    IF (p_parallel_io) THEN
       WRITE(kout,*)' '
       WRITE(kout,*)' Start couple_init'
       WRITE(kout,*)' *****************'
       WRITE(kout,*)' '
       CALL FLUSH(kout)           
    ENDIF
#endif

    ! get control variables for selecting features on global arrays

      IF (p_parallel_io) THEN

        ALLOCATE (gl_slf  (nlon,ngl))
        ALLOCATE (gl_slm  (nlon,ngl))
        ALLOCATE (gl_alake(nlon,ngl))

        ALLOCATE (gl_tsw   (nlon,ngl))
        ALLOCATE (gl_seaice(nlon,ngl))
        ALLOCATE (gl_siced (nlon,ngl))
        ALLOCATE (gl_sni   (nlon,ngl))
        ALLOCATE (gl_ocu   (nlon,ngl))
        ALLOCATE (gl_ocv   (nlon,ngl))
#ifdef __cpl_co2
        ALLOCATE (gl_co2trans(nlon,ngl))
        ALLOCATE (gl_co2ocean(nlon,ngl))
#endif

!wiso-code
        IF (lwiso) THEN
          ALLOCATE (gl_o16   (nlon,ngl))
          ALLOCATE (gl_o18   (nlon,ngl))
          ALLOCATE (gl_hdo   (nlon,ngl))
        END IF
!wiso-code-end

      ELSE

        gl_slf   => NULL()
        gl_slm   => NULL()
        gl_alake => NULL()

        gl_tsw    => NULL()
        gl_seaice => NULL()
        gl_siced  => NULL()
        gl_sni    => NULL()
        gl_ocu    => NULL()
        gl_ocv    => NULL()
#ifdef __cpl_co2
        gl_co2trans => NULL()
        gl_co2ocean => NULL()
#endif

!wiso-code
        IF (lwiso) THEN
          gl_o16    => NULL()
          gl_o18    => NULL()
          gl_hdo    => NULL()
        END IF
!wiso-code-end

      END IF

      CALL gather_gp(gl_slf, slf, gd)
      CALL gather_gp(gl_slm, slm, gd)
      CALL gather_gp(gl_alake, alake, gd)

    ! land values must be saved

      CALL gather_gp(gl_tsw, tsw, gd)
      CALL gather_gp(gl_seaice, seaice, gd)
      CALL gather_gp(gl_siced, siced, gd)
      CALL gather_gp(gl_sni, sni, gd)
      CALL gather_gp(gl_ocu, ocu, gd)
      CALL gather_gp(gl_ocv, ocv, gd)
#ifdef __cpl_co2
      CALL gather_gp(gl_co2trans, co2trans, gd)
      CALL gather_gp(gl_co2ocean, co2ocean, gd)
#endif

!wiso-code
      IF (lwiso) THEN
        CALL gather_gp(gl_o16, wisosw_d(:,1,:), gd)
        CALL gather_gp(gl_o18, wisosw_d(:,2,:), gd)
        CALL gather_gp(gl_hdo, wisosw_d(:,3,:), gd)
      END IF
!wiso-code-end

    !
    !-- Get model step when coupled run is started.
    !  
    couple_resum = get_time_step()
    IF (p_parallel_io) THEN
      WRITE(kout,*)' Coupled run is started at model step  ',couple_resum
      WRITE(kout,*)' '
    END IF

#ifdef __synout
    !
    !-- Check land/lake/sea masks.
    !  
    IF (p_parallel_io) CALL check_alake(kout)
    !
    !-- Print land/lake/sea masks.
    !  
    IF (p_parallel_io) CALL print_masks(kout,kerr)
#endif
    !
#if ! (defined (use_comm_MPI1) || defined (use_comm_MPI2))

     CALL message('','In coupled mode specification of either')
     CALL message('','-Duse_comm_MPI1 or -Duse_comm_MPI2 is required.')
     CALL finish('p_start','Aborting run, please recompile.')

#endif
#ifdef use_comm_MPI2
    !
    !--   Join the communicator group.
    ! 
    ierror = PRISM_Ok
    CALL prism_init_comp_proto (comp_id, modnam, ierror)
    IF (ierror /= PRISM_Ok) THEN
       CALL couple_abort (modnam,'couple_init','pb prism_init_comp',kerr)
    ENDIF
    !
#endif
    !
    !-- Write grids file for oasis
    !
    IF (p_parallel_io) THEN
       CALL grids_writing(kout)
    ENDIF
    !
#ifdef use_comm_MPI1
    commlocal=p_all_comm
#endif
#ifdef use_comm_MPI2
    !
    !-- Extract context for local communicator commlocal.
    !
    ierror   = PRISM_Ok
    CALL prism_get_localcomm_proto(commlocal, ierror)
    IF (ierror /= PRISM_Ok) THEN
          CALL couple_abort (modnam,'couple_init','pb get localcomm',kerr)
    ENDIF
    !
#endif
    !
    CALL MPI_Comm_Size(commlocal, nbtotproc, ierror)
    CALL MPI_Comm_Rank(commlocal, rank, ierror)

    !
    !-- Set exchange-field locators.
    !
#if defined __cpl_mpiom
    !
    ! Set locator for fields passed from the ocean to the atmosphere.
    !
    clstr8o2a(1:nfldo2a) = cg_cnamout(1:nfldo2a)
    !
    ! Set locator for fields passed from the atmosphere to the ocean.
    !
    clstr8a2o(1:nflda2o) = cg_cnaminp(nfldo2a+1:nfldo2a+nflda2o)

#endif

    !
    !-- Define decomposition strategy for field exchange: none.
    !
    nsegments = 1
    parsize = 3
    ALLOCATE(paral(parsize))
    paral ( clim_strategy ) = clim_serial 
    paral ( clim_length   ) = nlon*ngl
    paral ( clim_offset   ) = 0

    !
    !-- Associate a part_id to the decomposition of the process.
    !  
    IF (p_parallel_io) THEN
       ierror   = PRISM_Ok
       call prism_def_partition_proto (part_id, paral, ierror)
       IF (ierror /= PRISM_Ok) THEN
          CALL couple_abort (modnam,'couple_init','pb def_partition',kerr)
       ENDIF

    !
    !-- Define ports of incoming fields.
    !  
       var_nodims(1) = 1
       var_nodims(2) = 0
       var_shape(1)  = 1
       var_shape(2)  = paral(clim_length)

       DO jfld = 1,nfldo2a

          ierror = PRISM_Ok
          CALL prism_def_var_proto &
          (portin_id(jfld), clstr8o2a(jfld), &
          part_id, var_nodims, PRISM_In,  &
          var_shape, PRISM_Real, ierror )
          IF (ierror /= PRISM_Ok) THEN
             WRITE (kout, *) &
             'WARNING : Problem with import port ',clstr8o2a(jfld)
             WRITE (kout, *) &
             'WARNING : Port ID returned is : ',portin_id(jfld)
             WRITE (kout, *)' =======   Error code number = ',ierror
          ELSE
             WRITE (kout,*)'   '
             WRITE (kout,*)' couple_init: Import port ',clstr8o2a(jfld),' defined'
             WRITE (kout,*) &
             ' couple_init: Port ID returned is : ',portin_id(jfld)
             WRITE (kout,*)' With port status    ',PRISM_In
             WRITE (kout,*)' With port data type ',PRISM_Double
             WRITE (kout,*)' With decomposition  ',(paral(k),k=1,3)
          ENDIF 

       END DO

    !
    !--   Defined ports of outgoing fields.
    !  

       DO jfld = 1,nflda2o

          ierror = PRISM_Ok
          CALL prism_def_var_proto (portout_id(jfld), clstr8a2o(jfld), &
          part_id, var_nodims, PRISM_Out, &
          var_shape, PRISM_Real, ierror )
          IF (ierror /= PRISM_Ok) THEN
             WRITE (kout, *) &
             ' WARNING : Problem with export port ',clstr8a2o(jfld)
             WRITE (kout, *) &
             ' WARNING : Port ID returned is : ',portout_id(jfld)
             WRITE (kout, *)' =======   Error code number = ',ierror
          ELSE
             WRITE (kout,*)'   '
             WRITE (kout,*) &
             ' couple_init: Export port ',clstr8a2o(jfld),' defined'
             WRITE (kout,*) &
             ' couple_init: Port ID returned is : ',portout_id(jfld)
             WRITE (kout,*)' with port status      ',PRISM_Out
             WRITE (kout,*)' with port data type   ',PRISM_Double
             WRITE (kout,*)' with decomposition    ',(paral(k),k=1,3)
          ENDIF 

       END DO

       ierror = PRISM_Ok
       CALL prism_enddef_proto(ierror)
       IF (ierror /= PRISM_Ok) THEN
          WRITE (kout, *)' WARNING : Problem with prism_enddef_proto'
          WRITE (kout, *)' =======   Error code number = ',ierror
       ENDIF 

    ENDIF

    !
    !-- Initilize exchange fields.
    !

    CALL ini_a2o(kout)

    !
    !-- Check coupling-control parameters
    !   against those received from oasis.
    !  

    CALL chck_par(lresume,kout,kerr,nlon,ngl)
       
    IF (p_parallel_io) THEN
       nmseq = maxval (ig_def_seq)
    ENDIF
    CALL p_bcast(nmseq, p_io)

    !     
    !-- Print coupling parameters received from oasis3.
    !  
    IF (p_parallel_io) THEN 
       WRITE (kout,*)' Couple_init:'
       WRITE (kout,*)' -----------------------------------------------------'
       WRITE (kout,*)' Model name      = ',modnam
       WRITE (kout,*)' Model time step = ',delta_time
       WRITE (kout,*)' -----------------------------------------------------'
       WRITE (kout,*)' '       
       WRITE (kout,*)' Parameters received from oasis3 (mod_comprism_proto):'
       WRITE (kout,*)' mynum                         ', mynum
       WRITE (kout,*)' mytid                         ', mytid
       WRITE (kout,*)' Model number                  ', comp_id
       WRITE (kout,*)' ig_mynummod                   ', ig_mynummod
       WRITE (kout,*)' Total no. of processors       ', nbtotproc
       WRITE (kout,*)' No. of procs for data exchange', nbcplproc
       WRITE (kout,*)' Depomposition ID              ', part_id
       WRITE (kout,*)' minimum exchange frequency    ', ig_frqmin
       WRITE (kout,*)' total time of the simulation  ', ig_ntime
       WRITE (kout,*)' number of fields exchanged    ', ig_clim_nfield
       WRITE (kout,*)' initial date                  ', ig_inidate
       WRITE (kout,*)' Lag of exported fields        ', ig_def_lag
       WRITE (kout,*)' coupling period of fields     ', ig_def_freq
       WRITE (kout,*)' sequential index of fields    ', ig_def_seq
       WRITE (kout,*)' ig_def_norstfile              ', ig_def_norstfile
       WRITE (kout,*)' number of restart files       ', ig_nbr_rstfile
       WRITE (kout,*)' name of restart files         ', cg_def_rstfile
       WRITE (kout,*)' restart files in netcdf       ', lg_ncdfrst
!!lk       WRITE (kout,*)' ig_aux                        ', ig_aux
       WRITE (kout,*)' name of exchanged fields (in) ', cg_cnaminp
       WRITE (kout,*)' name of exchanged fields (out)', cg_cnamout
       WRITE (kout,*)' cg_ignout_field               ', cg_ignout_field
       WRITE (kout,*)' any field going through Oasis?', lg_oasis_field
       WRITE (kout,*)' '
       WRITE(kout,*)' '
       IF(nmseq /= 1) THEN
          WRITE(kout,*)' The models run sequentially !'
          WRITE(kout,*)' No. of sequential exchange fields: ',nmseq
       ELSE
          WRITE(kout,*)' All fields have sequential no. 1 !'
          WRITE(kout,*)' The models run concurrently!'
       ENDIF
       WRITE(kout,*)' '
    END IF

#ifdef __synout
    IF (p_parallel_io) THEN
       WRITE(kout,*)' '
       WRITE(kout,*)' End of couple_init'
       WRITE(kout,*)' ******************'
       WRITE(kout,*)' '
       CALL FLUSH(kout)
    ENDIF
#endif

#endif

  END SUBROUTINE couple_init

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  couple_get_o2a
!
! !INTERFACE:

  SUBROUTINE couple_get_o2a(kout,kerr)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kout      ! unit number for standard output
    INTEGER, INTENT(in) :: kerr      ! error out unit

! !DESCRIPTION:
!  - Receives ocean data from the coupling system; called in stepon.
!
! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

    INTEGER   :: istep ! model step since initialisation

#if defined (__prism)

    REAL(pwp) :: exfld(nlon,ngl) ! buffer for receiving exchange fields

#ifdef __synout
      IF (p_parallel_io) THEN
         WRITE(kout,*)' '
         WRITE(kout,*)' Start of couple_get_o2a'
         WRITE(kout,*)' ***********************'
         WRITE(kout,*)' '
         CALL FLUSH(kout)
      END IF
#endif
      istep = get_time_step()

#ifdef __synout
      IF (p_parallel_io) THEN
         WRITE(kout,*)                 ' model step since init. : ',istep
         WRITE(kout,*)                 ' model run step         : ',istep-couple_resum
         WRITE(kout,*)                 ' run time (sec) passed  : ',run_date_secs
         CALL write_date(next_date,    ' next_date              : ')
         WRITE(kout,*)' '
      END IF
#endif

      DO jfld = 1,nfldo2a

         info = PRISM_Ok
         IF ( p_parallel_io ) THEN
#ifdef __synout
            WRITE(kout,*)'    '
            WRITE(kout,*)' prism_get field with ',clstr8o2a(jfld),' ...'
#endif
            CALL prism_get_proto(portin_id(jfld), run_date_secs, exfld, info)

            CALL digest_get_Id &
                 (kout,kerr,info,clstr8o2a(jfld),run_date_secs)
 
         ENDIF 
         CALL p_bcast(info, p_io)

         IF (info == PRISM_Recvd        .OR. &
             info == PRISM_FromRest     .OR. &
             info == PRISM_RecvOut      .OR. &
             info == PRISM_FromRestOut         ) THEN
            CALL put_o2a(clstr8o2a(jfld),exfld,nlon,ngl,kout,kerr)
         ENDIF 

      ENDDO

      !
      !-- Accumulate time passed since last averaging.
      !
       couple_a2o_time = couple_a2o_time + delta_time


#ifdef __synout
      IF (p_parallel_io) THEN
         WRITE(kout,*)' '
         WRITE(kout,*)' End of couple_get_o2a'
         WRITE(kout,*)' *********************'
         WRITE(kout,*)' '
         CALL FLUSH(kout)
      END IF
#endif

#endif

  END SUBROUTINE couple_get_o2a

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  couple_put_a2o
!
! !INTERFACE:

  SUBROUTINE couple_put_a2o(kout,kerr)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kout      ! unit number for standard output
    INTEGER, INTENT(in) :: kerr      ! error out unit
! !DESCRIPTION:
!
! - Sends exchange fields (for ocean) to oasis.
!
! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------
 
#if defined (__prism)

    REAL(pwp) :: exfld(nlon,ngl) ! array sent by prism_put

    INTEGER :: istep ! model step since initialisation

#ifdef __synout
      IF (p_parallel_io) THEN
         WRITE(kout,*)' '
         WRITE(kout,*)' Start of couple_put_a2o'
         WRITE(kout,*)' ***********************'
         WRITE(kout,*)' '
         CALL FLUSH(kout)
      END IF
#endif
      istep = get_time_step()

#ifdef __synout
      IF (p_parallel_io) THEN
         WRITE(kout,*)' Run time (sec) passed                   : ',run_date_secs
         WRITE(kout,*)' Accum. time (sec) since last couple_put : ',couple_a2o_time
         WRITE(kout,*)                 ' model run step         : ',istep-couple_resum
         WRITE(kout,*)                 ' model step since init. : ',istep
         CALL FLUSH(kout)
      END IF
#endif
      !
      !-- Export exchange fields
      !  
      IF (l_putocean) THEN

        CALL water_budget_corr

        DO jfld = 1,nflda2o

           CALL get_a2o(clstr8a2o(jfld),exfld,nlon,ngl,kout,kerr)

           IF (p_parallel_io) THEN
#ifdef __synout
              WRITE(kout,*)' '
              WRITE(kout,*)' Exporting field no. ',jfld,':',clstr8a2o(jfld)
              CALL FLUSH(kout)
#endif
              info = PRISM_Ok
              CALL prism_put_proto(portout_id(jfld), run_date_secs, exfld, info)

              CALL digest_put_Id &
              (kout,kerr,info,clstr8a2o(jfld),run_date_secs)
           END IF

        ENDDO

        !
        !-- Reset exchange fields and accumulation interval.
        !  

        IF ( .NOT. lbreak .OR. .NOT. l_putrerun ) THEN
           CALL ini_a2o(kout)
           couple_a2o_time = 0.0_pwp
        END IF

      ELSE

#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*)' No call of prism_put at date (sec)', run_date_secs
            CALL FLUSH(kout)
         END IF
#endif
      END IF

#ifdef __synout
      IF (p_parallel_io) THEN
         WRITE(kout,*)' '
         WRITE(kout,*)' End of couple_put_a2o'
         WRITE(kout,*)' *********************'
         WRITE(kout,*)' '
         CALL FLUSH(kout)
      END IF
#endif

#endif

  END SUBROUTINE couple_put_a2o
  
!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  couple_end
!
! !INTERFACE:

  SUBROUTINE couple_end(kout,kerr)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kout      ! unit number for standard output
    INTEGER, INTENT(in) :: kerr      ! unit number for error output

! !DESCRIPTION:
!
! - Terminates the coupled run.
!

! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------
 
#if defined (__prism)
 
    REAL(pwp) :: exfld(nlon,ngl) ! array sent by prism_put

    REAL(pwp) ::  zrcouple

#ifdef __synout
    IF (p_parallel_io) THEN
       WRITE(kout,*)' '
       WRITE(kout,*)' Start of couple_end '
       WRITE(kout,*)' *******************'
       WRITE(kout,*)' '
       CALL FLUSH(kout)           
    END IF
#endif

    IF (.NOT.l_putocean) THEN

       IF (p_parallel_io) THEN
          WRITE(kout,*)' Warning : this run does not stop ', &
                       ' at the end of coupled time step.'
          WRITE(kout,*)' Time (sec) passed since last put : ',couple_a2o_time
          WRITE(kout,*)' Date in prism_put                : ',ig_ntime
          CALL FLUSH(kout)
       END IF
       zrcouple = 1.0_pwp/couple_a2o_time
       awust(:,:) = awust(:,:)*zrcouple
       awvst(:,:) = awvst(:,:)*zrcouple
       aiust(:,:) = aiust(:,:)*zrcouple
       aivst(:,:) = aivst(:,:)*zrcouple
       aifre(:,:) = aifre(:,:)*zrcouple/rhoh2o
       awfre(:,:) = awfre(:,:)*zrcouple/rhoh2o
       aiqre(:,:) = aiqre(:,:)*zrcouple
       aicon(:,:) = aicon(:,:)*zrcouple
       awhea(:,:) = awhea(:,:)*zrcouple
       awsol(:,:) = awsol(:,:)*zrcouple
       awsta(:,:) = awsta(:,:)*zrcouple
!---wiso-code
       IF (lwiso) THEN
         wisoaifre(:,:,:) = wisoaifre(:,:,:)*zrcouple/rhoh2o
         wisoawfre(:,:,:) = wisoawfre(:,:,:)*zrcouple/rhoh2o
       END IF
!---wiso-code
       CALL couple_put_a2o(kout,kerr)
    ENDIF

!-- Write restart file.
!
    IF ( nmseq > 1) THEN
       CALL couple_restart_a2o(kout,kerr)
    END IF

!CSL  DO jfld = 1,nflda2o
!CSL
!CSLCSL     IF ( nmseq > 1) THEN
!CSLCSL
!CSLCSL       CALL get_a2o(clstr8a2o(jfld),exfld,nlon,ngl,kout,kerr)
!CSL
!CSL       IF (p_parallel_io) THEN
!CSL#ifdef __synout
!CSL          WRITE(kout,*)' Writing to restart : ',clstr8a2o(jfld)
!CSL          WRITE(kout,*)' sequential index   : ',ig_def_seq(portout_id(jfld))
!CSL          CALL FLUSH(kout)
!CSL#endif
!CSL          info = PRISM_Ok
!CSLCSL          CALL prism_put_proto(portout_id(jfld), ig_ntime, exfld, ierror)
!CSL
!CSL          CALL prism_put_restart_proto(portout_id(jfld), run_date_secs, info)
!CSL          CALL digest_put_Id &
!CSL              (kout,kerr,info,clstr8a2o(jfld),run_date_secs)
!CSL
!CSLCSL          IF (ierror /= PRISM_Ok) THEN
!CSLCSL             WRITE (kerr, *)'couple_put.prism_put:'
!CSLCSL             WRITE (kerr, *)'Problem with port/no. ',clstr8a2o(jfld) &
!CSLCSL                             ,portout_id(jfld)
!CSLCSL             WRITE (kerr, *)'Error code number = ',ierror
!CSLCSL             WRITE (kerr, *)'Restart writing.'
!CSLCSL             CALL couple_abort (modnam,'couple_put','pb prism_put',kerr)
!CSLCSL#ifdef __synout
!CSLCSL          ELSE
!CSLCSL            WRITE(kout,*)' '
!CSLCSL            WRITE(kout,*)' Field written to restart : ',clstr8a2o(jfld)
!CSLCSL            CALL FLUSH(kout)
!CSLCSL#endif
!CSLCSL          ENDIF 
!CSL       END IF
!CSL     END IF
!CSL  ENDDO

  !     
  !-- Deallocate memory.
  !
  ierror = PRISM_Ok
  CALL prism_terminate_proto( ierror )
  IF (ierror /= PRISM_Ok) THEN
     WRITE (kerr, *)'couple_end: rank=',rank
     WRITE (kerr, *)'pb prism_terminate '
     CALL couple_abort (modnam,'couple_end','pb prism_terminate',kerr)
  END IF

#ifdef __synout
  IF (p_parallel_io) THEN
     WRITE(kout,*)' '
     WRITE(kout,*)' End of couple_end '
     WRITE(kout,*)' *****************'
     WRITE(kout,*)' '
     CALL FLUSH(kout)           
  END IF
#endif

#endif

  END SUBROUTINE couple_end

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  couple_restart_a2o
!
! !INTERFACE:

  SUBROUTINE couple_restart_a2o(kout,kerr)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kout      ! unit number for standard output
    INTEGER, INTENT(in) :: kerr      ! unit number for error output

! !DESCRIPTION:
!
! - Writes restart file with exchange fields at end of coupled run.
!

! !REVISION HISTORY:
! 03.08.19  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------
 
#if defined (__prism)
 
#ifdef __synout
    IF (p_parallel_io) THEN
       WRITE(kout,*)' '
       WRITE(kout,*)' Start of couple_restart_a2o '
       WRITE(kout,*)' ****************************'
       WRITE(kout,*)' '
       CALL FLUSH(kout)           
    END IF
#endif

      IF (p_parallel_io) THEN
         DO jfld = 1,nflda2o

            info = PRISM_Ok
            CALL prism_put_restart_proto &
                 (portout_id(jfld), run_date_secs, info)
            CALL digest_put_Id &
                 (kout,kerr,info,clstr8a2o(jfld),run_date_secs)
         ENDDO
      END IF
#ifdef __synout
  IF (p_parallel_io) THEN
     WRITE(kout,*)' '
     WRITE(kout,*)' End of couple_restart_a2o '
     WRITE(kout,*)' **************************'
     WRITE(kout,*)' '
     CALL FLUSH(kout)           
  END IF
#endif

#endif

  END SUBROUTINE couple_restart_a2o

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  check_alake
!
! !INTERFACE:

  SUBROUTINE check_alake(kout)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kout      ! unit number for standard output

! !DESCRIPTION:
!
! - Checks cell partitioning into lake, land, ocean surfaces.
!
! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------
 
    INTEGER :: jl, jrow

#if defined (__prism)

!!CSL    REAL(dp) :: cellarea

!CSL    IF (p_parallel) RETURN ! check not implemented

    ! lakes are only allowed on cells without ocean.
    ! therefore check whether this is fullfilled [else
    ! enforce by setting lake part of cell to zero if ocean part > 0]
    ! [Also: check whether slf+alake=1 in non-ocean cells.]

    DO jl=1,nlon
      DO jrow=1,ngl

        !  Print cells with lake
!CSL        IF ( gl_alake(jl,jrow) > 0.0_dp ) THEN
!CSL          WRITE(kout,*)'lake in cell ',jl,jrow,' : ',gl_alake(jl,jrow),gl_slf(jl,jrow)
!CSL          WRITE(kout,*)'lake in cell ',jl,jrow,' : ',gl_alake(jl,jrow)
!CSL          !  Modify if inconsistent
!CSL          IF ( gl_slf(jl,jrow)+gl_alake(jl,jrow)>1.) THEN
!CSL            WRITE(kout,*)'Inconsistent land/lake partioning at lon/lat !'
!CSL            cellarea=gl_slf(jl,jrow)+gl_alake(jl,jrow)
!CSL            gl_slf(jl,jrow)=gl_slf(jl,jrow)/cellarea
!CSL            gl_alake(jl,jrow)=gl_alake(jl,jrow)/cellarea
!CSL            WRITE(kout,*)'Changed to :',gl_alake(jl,jrow),gl_slf(jl,jrow)
!CSL          END IF
!CSL#ifdef __cpl_hope
!CSL          !  No lakes allowed: move lake surface to land surface
!CSL          WRITE(kout,*)'Add alake surface to slf at lon/lat index : ',jl,jrow
!CSL          gl_slf(jl,jrow) = gl_slf(jl,jrow)+gl_alake(jl,jrow)
!CSL          gl_alake(jl,jrow)=0.0_dp
!CSL          IF ( gl_slf(jl,jrow) < 1.) THEN
!CSL            WRITE(kout,*)'Unexpected sea/land : ',gl_alake(jl,jrow),gl_slf(jl,jrow)
!CSL          END IF
!CSL#endif
!CSL        END IF
      END DO
    END DO
      
#endif
 
  END SUBROUTINE check_alake

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  print_sst
!
! !INTERFACE:

  SUBROUTINE print_sst(kout)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kout      ! unit number for standard output

! !DESCRIPTION:
!
! - Prints window of SST (array tsw) around min/max.
!
! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------
    
#if defined (__prism)

    REAL(dp) :: tswgmax, tswgmin, tswmax, tswmin
    INTEGER :: igmax, igmin, jgmax, jgmin, imax, imin
    INTEGER :: jx, jy, jpdx, jpdy

    jpdy =ngl
    jpdx =nlon
    tswgmax = -999.0_dp
    tswgmin = +999.0_dp
    igmax = 1
    igmin = 1
    jgmax = 1
    jgmin = 1
    DO JY = 1, jpdy
      tswmax = -999.0_dp
      tswmin = +999.0_dp
      imax = 1
      imin = 1
      DO JX = 1,jpdx
        IF(gl_tsw(jx,jy) > tswmax) THEN 
          IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
            !                  IF(gl_tsw(jx,jy) > 1.) THEN
            tswmax = gl_tsw(jx,jy)
            imax = jx
          ENDIF
        ENDIF
        IF(gl_tsw(jx,jy)  <  tswmin) THEN 
          IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
            !                  IF(gl_tsw(jx,jy) > 1.) THEN
            tswmin = gl_tsw(jx,jy)
            imin = jx
          ENDIF
        ENDIF
      ENDDO
      IF(tswmax > tswgmax) THEN 
        tswgmax = tswmax
        igmax = imax
        jgmax = jy
      ENDIF
      IF(tswmin < tswgmin) THEN 
        tswgmin = tswmin
        igmin = imin
        jgmin = jy
      ENDIF
      WRITE(kout,'(a,i5,a)', advance = 'no') 'Max/min of SST at j= ',jy,' : '
      WRITE(kout,'(2i5,2f12.7)')  imax,imin,gl_tsw(imax,jy),gl_tsw(imin,jy)
    ENDDO

    WRITE(kout,*)' '
    WRITE(kout,*)' SST near global max :',igmax,jgmax,gl_tsw(igmax,jgmax)
    WRITE(kout,*)' ---------------------------------------------------'
    DO JY = max(jgmax-5,1),min(jgmax+5,jpdy)
      WRITE(kout,6001)(gl_tsw(jx,jy),jx=max(igmax-5,1),min(igmax+5,jpdx))
    ENDDO

    WRITE(kout,*)' '
    WRITE(kout,*)' slf there (1/0=land/sea only):',gl_slf(igmax,jgmax)
    WRITE(kout,*)' ------------------------------------------------'
    DO JY = max(jgmax-5,1),min(jgmax+5,jpdy)
      WRITE(kout,6002)(gl_slf(jx,jy),jx=max(igmax-5,1),min(igmax+5,jpdx))
    ENDDO

    WRITE(kout,*)'  '
    WRITE(kout,*)' SST near global min :',igmin,jgmin,gl_tsw(igmin,jgmin)
    WRITE(kout,*)' ---------------------------------------------------'
    DO JY = max(jgmin-5,1),min(jgmin+5,jpdy)
      WRITE(kout,6001)(gl_tsw(jx,jy),jx=max(igmin-5,1),min(igmin+5,jpdx))
    ENDDO

    WRITE(kout,*)'  '
    WRITE(kout,*)' slf there (1/0=land/sea only):',gl_slf(igmin,jgmin)
    WRITE(kout,*)' ------------------------------------------------'
    DO JY = max(jgmin-5,1),min(jgmin+5,jpdy)
      WRITE(kout,6002)(gl_slf(jx,jy),jx=max(igmin-5,1),min(igmin+5,jpdx))
    ENDDO

6001 FORMAT (11(1x,f7.2,1x))
6002 FORMAT (11(1x,f7.4,1x))

    WRITE (kout, '(/)')

    tswgmax = -999.0_dp
    tswgmin = +999.0_dp
    igmax = 1
    igmin = 1
    jgmax = 1
    jgmin = 1
    DO JY = 1, jpdy
      tswmax = -999.0_dp
      tswmin = +999.0_dp
      imax = 1
      imin = 1
      DO JX = 1,jpdx
        IF(gl_tsw(jx,jy) > tswmax) THEN 
          IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy)==0.0_dp) THEN
            tswmax = gl_tsw(jx,jy)
            imax = jx
          ENDIF
        ENDIF
        IF(gl_tsw(jx,jy) < tswmin) THEN 
          IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
            tswmin = gl_tsw(jx,jy)
            imin = jx
          ENDIF
        ENDIF
      ENDDO
      IF(tswmax > tswgmax) THEN 
        tswgmax = tswmax
        igmax = imax
        jgmax = jy
      ENDIF
      IF(tswmin < tswgmin) THEN 
        tswgmin = tswmin
        igmin = imin
        jgmin = jy
      ENDIF
      WRITE(kout,'(a,i5,a)', advance = 'no') 'Max/min of SST at j= ',jy,' : '
      WRITE(kout,'(2i5,2f12.7)')  imax,imin,gl_tsw(imax,jy),gl_tsw(imin,jy)
    ENDDO

    WRITE (kout, '(/)')

#endif

  END SUBROUTINE print_sst

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  print_masks
!
! !INTERFACE:

  SUBROUTINE print_masks(kout,kerr)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kout      ! unit number for standard output
    INTEGER, INTENT(in) :: kerr      ! error out unit

! !DESCRIPTION:
!
! - Prints land/sea mask characteristics.
!
! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

    REAL(dp)    :: slfgmin, slfgmax, slfmin, slfmax, slfsum, alakesum

    INTEGER :: igmin, igmax, imin, imax
    INTEGER :: jgmin, jgmax
    INTEGER :: jx, jy

    CHARACTER(len=3) :: cmask(nlon,ngl)

#if defined (__prism)

    !-- Land mask characteristics

    slfgmin = 2.0_dp
    slfgmax =-1.0_dp
    DO JY = 1, ngl
      imax = 0
      imin = 0
      slfsum = 0.0_dp
      alakesum = 0.0_dp
      slfmin = 2.0_dp
      slfmax =-1.0_dp
      DO JX = 1,nlon
        IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
          IF(gl_slf(jx,jy) > slfmax) THEN
            slfmax = gl_slf(jx,jy)
            imax = jx
          ENDIF
        ENDIF
        IF(gl_slf(jx,jy) > 0.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
          IF(gl_slf(jx,jy) < slfmin) THEN
            slfmin = gl_slf(jx,jy)
            imin = jx
          ENDIF
        ENDIF
        slfsum   = slfsum   + gl_slf(jx,jy)
        alakesum = alakesum + gl_alake(jx,jy)
        IF(slfmax > slfgmax) THEN 
          slfgmax = slfmax
          igmax = imax
          jgmax = jy
        ENDIF
        IF(slfmin  <  slfgmin) THEN 
          slfgmin = slfmin
          igmin = imin
          jgmin = jy
        ENDIF
        !            WRITE(kout,*)' Sum of slf at ',jy,' : ', slfsum
        !            IF(alakesum /= 0.0_dp) &
        !                 WRITE(kout,*)' Warning: Sum of alake  : ', alakesum
!CSL        IF(gl_alake(jx,jy) > 0.0_dp) THEN
!CSL          IF(gl_slf(jx,jy)+gl_alake(jx,jy) > 1.0_dp) THEN
!CSL            WRITE(kout,*)' Warning: Correcting slf at ', &
!CSL                 jx,jy,gl_slf(jx,jy),gl_alake(jx,jy)
!CSL            gl_slf(jx,jy)=1.-gl_alake(jx,jy)
!CSL          ENDIF
!CSL        ENDIF
      ENDDO
    ENDDO
    WRITE(kout,*)' Global min of fractional land area on cells ' &
         ,'with dry part but without lakes : ' &
         ,slfgmin,igmin,jgmin
    WRITE(kout,*)' Global max of fractional land area on cells ' &
         ,'with wet part but without lakes : ' &
         ,slfgmax,igmax,jgmax

    !       Print land mask

    DO JY = 1, ngl
      DO JX = 1,nlon
        IF     (gl_slf(jx,jy) == 0.0_dp) THEN
          cmask(jx,jy) = '  .'
        ELSEIF (gl_slf(jx,jy) == 1.0_dp) THEN
          cmask(jx,jy) = ' **'
        ELSEIF (gl_slf(jx,jy)<1.0_dp .AND. gl_slf(jx,jy) > 0.0_dp) THEN
          WRITE(cmask(jx,jy),'(1X,I2.2)') INT(gl_slf(jx,jy)*100.0_dp)
        ELSE
          WRITE(kout,*)' Invalid fractional land area on ' &
               ,jx,jy,gl_slf(jx,jy)
          CALL couple_abort (modnam,'print_masks','pb fractional land mask',kerr)
        ENDIF
      ENDDO
    ENDDO
    WRITE(kout,*)' Partial land mask (%, jx=1,64) :'
    WRITE(kout,*)'   .: no land / **: no water / >%'
    WRITE(kout,*)' '
    WRITE(kout,6002) (mod(jx,10),jx=1,nlon/2)
    DO JY = 1, ngl
      WRITE(kout,6003) (cmask(jx,jy),jx=1,nlon/2)
    ENDDO
    WRITE(kout,*)' '
    WRITE(kout,*)' Partial sea/land mask (%, jx=65,128) :'
    WRITE(kout,*)' '
    WRITE(kout,6002) (mod(jx,10),jx=nlon/2+1,nlon)
    DO JY = 1, ngl
      WRITE(kout,6003) (cmask(jx,jy),jx=nlon/2+1,nlon)
    ENDDO

    !       Print lake mask

    DO JY = 1, ngl
      DO JX = 1,nlon
        cmask(jx,jy) = '  .'
        IF(gl_slf(jx,jy) > 0.0_dp) cmask(jx,jy) = ' **'
        IF(gl_alake(jx,jy) > 0.0_dp) THEN
          WRITE(cmask(jx,jy),'(1X,I2)') INT(gl_alake(jx,jy)*100.0_dp+1.0_dp)
        ENDIF
      ENDDO
    ENDDO
    WRITE(kout,*)' Lake mask (%, jx=1,64) :'
    WRITE(kout,*)'   .: no land / **: no lake / <%'
    WRITE(kout,*)' '
    WRITE(kout,6002) (mod(jx,10),jx=1,nlon/2)
    DO JY = 1, ngl
      WRITE(kout,6003) (cmask(jx,jy),jx=1,nlon/2)
    ENDDO
    WRITE(kout,*)' '
    WRITE(kout,*)' Lake mask (%, jx=65,128) :'
    WRITE(kout,*)' '
    WRITE(kout,6002) (mod(jx,10),jx=nlon/2+1,nlon)
    DO JY = 1, ngl
      WRITE(kout,6003) (cmask(jx,jy),jx=nlon/2+1,nlon)
    ENDDO

    !       Print ocean mask

    DO JY = 1, ngl
      DO JX = 1,nlon
        cmask(jx,jy) = '  .'
        IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
          WRITE(cmask(jx,jy),'(1X,I2)') INT((1.-gl_slf(jx,jy))*100.0_dp)
        ENDIF
      ENDDO
    ENDDO
    WRITE(kout,*)' Ocean mask (%, jx=1,64) :'
    WRITE(kout,*)'   .: no ocean / **: only / >%'
    WRITE(kout,*)' '
    WRITE(kout,6002) (mod(jx,10),jx=1,nlon/2)
    DO JY = 1, ngl
      WRITE(kout,6003) (cmask(jx,jy),jx=1,nlon/2)
    ENDDO
    WRITE(kout,*)' '
    WRITE(kout,*)' Ocean mask (%, jx=65,128) :'
    WRITE(kout,*)' '
    WRITE(kout,6002) (mod(jx,10),jx=nlon/2+1,nlon)
    DO JY = 1, ngl
      WRITE(kout,6003) (cmask(jx,jy),jx=nlon/2+1,nlon)
    ENDDO
6003 FORMAT (64(A3))
6002 FORMAT (64(1X,I2))

#endif

  END SUBROUTINE print_masks

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  chck_par
!
! !INTERFACE:

  SUBROUTINE chck_par(lmres,kout,kerr,klon,kgl)

! !USES:
!
#ifdef NAG
    USE f90_unix,           ONLY: flush
#endif  

    USE mo_time_conversion, ONLY: print_date, OPERATOR(==), OPERATOR(<)

! !RETURN VALUE:

    IMPLICIT NONE

    INTEGER, INTENT(in) :: klon, kgl ! dimensions of model grid
    INTEGER, INTENT(in) :: kout      ! unit number for standard output
    INTEGER, INTENT(in) :: kerr      ! error out unit

    LOGICAL, INTENT(in) :: lmres     ! logical switch for restart mode

! !DESCRIPTION:
!
! - Checks coupled model control parameter against
!   those of the calling model ECHAM5.
!
! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

#if defined (__prism)

    INTEGER             :: nfldsize  ! fields size of exchange fields

    CHARACTER(len=32)   :: date_text

    nfldsize = klon*kgl
    
    !     
    !-- Check whether month of restart fiel is ok.
    !   ------------------------------------------
    
    IF (stop_date < next_date .OR. stop_date == next_date) THEN
      WRITE(kerr,*) 'stop date <= next date!'
      CALL print_date(next_date, mess=date_text)
      WRITE(kerr,*) 'next date = '//TRIM(date_text)
      CALL print_date(stop_date, mess=date_text)
      WRITE(kerr,*) 'stop date = '//TRIM(date_text)
      WRITE(kerr,*) 'STOP in chck_par'
      CALL couple_abort (modnam,'chck_par','pb with dates',kerr)
    ENDIF

    !     
    !-- Check restart mode.
    !   -------------------

    IF (lmres) THEN
      WRITE(kout,*) ' Model is restarted!'
    ELSE
      WRITE(kout,*) ' Model is initialized!'
    ENDIF
    WRITE(kout,*) ' '

#endif

  END SUBROUTINE chck_par

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  grids_writing
!
! !INTERFACE:

  SUBROUTINE grids_writing(kout)

! !RETURN VALUE:

! !USES:
  USE mo_gaussgrid,           ONLY: philon, philat, gridarea 

#if defined (__prism)
  USE mod_prism_grids_writing
#endif

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: kout            ! unit number for standard output
  INTEGER, PARAMETER   :: nc = 4          ! number of corners per grid cell
  INTEGER              :: gwrite          ! flag to state whether grids writing
                                          ! is needed or not. (1 / 0)
  INTEGER              :: i, j            ! looping indices
  CHARACTER*4          :: grdacr          ! grid acronym (as used in namcouple)
  CHARACTER*4          :: grdacr_2        ! grid acronym (as used in namcouple)
  REAL(dp)             :: lon(nlon,ngl)     ! 2dim array of longitudes
  REAL(dp)             :: lat(nlon,ngl)     ! 2dim array of latitudes 
  REAL(dp)             :: clon(nlon,ngl,nc) ! 3dim array of corner longitudes
  REAL(dp)             :: clat(nlon,ngl,nc) ! 3dim array of corner latitudes 
  INTEGER              :: mask(nlon,ngl)    ! land see mask (0 for all cells
                                            !    with wet area fraction > 0.) 
  INTEGER              :: mask_slm(nlon,ngl)! land see mask (corresponding to
                                            !    the echam land sea mask slm) 
  REAL(dp)             :: area(nlon,ngl)    ! grid cell areas corresponding to slf
  REAL(dp)             :: area_slm(nlon,ngl)!           corresponding to slm
  REAL(dp)             :: slf_sn(nlon,ngl)! land see mask in oasis order (S->N)
!
! !DESCRIPTION:
!
! - Write grids, masks and areas for oasis.
!
! !REVISION HISTORY:
! July 7, 2003  V. Gayler - created
!
! EOP
!-----------------------------------------------------------------------

#if defined (__prism)

#ifdef __synout
   WRITE(kout,*)' '
   WRITE(kout,*)' Start of grids_writing'
   WRITE(kout,*)' **************************'
   WRITE(kout,*)' '
#endif

#ifdef __synout
  WRITE(kout,*)'prism_start_grids_writing'
#endif
  call prism_start_grids_writing(gwrite)

  IF (gwrite == 1) THEN
!
!-- create 2d arrays of longitudes and latitudes following OASIS convention
!   (S -> N, E -> W)

#if defined (__noinvert)
     DO i = 1, nlon
        lat(i,:) = philat(:)
     ENDDO
     DO j = 1, ngl
        lon(:,j) = philon(:)
     ENDDO
#else
     DO i = 1, nlon
        DO j = 1, ngl
           lat(i,j) = philat(ngl+1-j)
        ENDDO
     ENDDO
     DO j = 1, ngl
        lon(:,ngl+1-j) = philon(:)
     ENDDO
#endif
     WHERE (lon(:,:) < 0._dp)
        lon = lon + 360._dp
     END WHERE
     WHERE (lon(:,:) >= 360._dp)
        lon = lon - 360._dp
     END WHERE
!
!-- create 3d arrays of grid cell corner longitudes and latitudes
!   (grid cell corners must be written in counterclockwise sense)
!      The writing of corners is optional. If they are missing in the grids
!      file they will be calculated by scrip. In this case it is better to
!      calculate the corners by the model. (Scrip-corner will not reach the
!      poles.) 

     DO i = 1, nlon-1
        clon(i,:,1) = 0.5_dp * (lon(i+1,1) + lon(i,1))
        clon(i,:,4) = 0.5_dp * (lon(i+1,1) + lon(i,1))
     ENDDO
     DO i = 2, nlon
        clon(i,:,2) = 0.5_dp * (lon(i-1,1) + lon(i,1))
        clon(i,:,3) = 0.5_dp * (lon(i-1,1) + lon(i,1))
     ENDDO
     clon(nlon,:,1) = 0.5_dp * (lon(1,1) + lon(nlon,1) + 360._dp)
     clon(nlon,:,4) = 0.5_dp * (lon(1,1) + lon(nlon,1) + 360._dp)
     clon(   1,:,2) = 0.5_dp * (lon(nlon,1) + lon(1,1) - 360._dp)
     clon(   1,:,3) = 0.5_dp * (lon(nlon,1) + lon(1,1) - 360._dp)
     IF (clon(nlon,1,1) >= 360._dp) clon(nlon,:,1) = clon(nlon,:,1) - 360._dp
     IF (clon(nlon,1,4) >= 360._dp) clon(nlon,:,4) = clon(nlon,:,4) - 360._dp
     IF (clon(nlon,1,2) < 0._dp) clon(nlon,:,2) = clon(nlon,:,2) + 360._dp
     IF (clon(nlon,1,3) < 0._dp) clon(nlon,:,3) = clon(nlon,:,3) + 360._dp

#if defined (__noinvert)
     DO j = 1, ngl-1
        clat(:,j,3) = 0.5_dp * (lat(1,j+1) + lat(1,j))
        clat(:,j,4) = 0.5_dp * (lat(1,j+1) + lat(1,j))
     ENDDO
     DO j = 2, ngl
        clat(:,j,1) = 0.5_dp * (lat(1,j-1) + lat(1,j))
        clat(:,j,2) = 0.5_dp * (lat(1,j-1) + lat(1,j))
     ENDDO
     clat(:,1,1) = 90._dp
     clat(:,1,2) = 90._dp
     clat(:,ngl,3) = -90._dp
     clat(:,ngl,4) = -90._dp
#else
     DO j = 2, ngl
        clat(:,j,1) = 0.5_dp * (lat(1,j-1) + lat(1,j))
        clat(:,j,2) = 0.5_dp * (lat(1,j-1) + lat(1,j))
     ENDDO
     DO j = 1, ngl-1
        clat(:,j,3) = 0.5_dp * (lat(1,j+1) + lat(1,j))
        clat(:,j,4) = 0.5_dp * (lat(1,j+1) + lat(1,j))
     ENDDO
     clat(:,1,1) = -90._dp
     clat(:,1,2) = -90._dp
     clat(:,ngl,3) = 90._dp
     clat(:,ngl,4) = 90._dp
#endif

!
!-- create 2d land ocean mask following OASIS convention (S -> N, E -> W)
!
#if defined (__noinvert)
     slf_sn(:,:) = gl_slf(:,:) + gl_alake(:,:)
#else
     slf_sn = 0._dp
     DO j = 1, ngl
        slf_sn(:,j) = gl_slf(:,ngl+1-j) + gl_alake(:,ngl+1-j)
     ENDDO
#endif
!
!-- create 2d integer sea land mask
!
     mask(:,:) = 0
     WHERE (slf_sn(:,:) > 0.9999999_dp)
        mask(:,:) = 1
     END WHERE

#if defined (__noinvert)
     mask_slm(:,:) = INT(gl_slm(:,:)) + INT(gl_alake(:,:))
#else
     DO j = 1, ngl
        mask_slm(:,j) = INT(gl_slm(:,ngl+1-j)) + INT(gl_alake(:,ngl+1-j))
     ENDDO
#endif
!
!-- create 2d array of of grid cell area
!     
!     The area array is used for CONSERV analysis of OASIS. 
!     In ECHAM5, fluxes over water and fluxes over land are calculated
!     separately. But the final flux going into the lowest grid cell is either
!     the flux calculated over water or the flux calculated over land, 
!     depending on the 1/0-land-sea-mask (slm); (compare ioinitial).
!     In this configuration CONSERV analysis is used for arrays going
!     FROM THE ATMOSPHERE TO THE OCEAN. Here all grid cells with a wet 
!     area fraction greater 0 are considered.

#if defined (__noinvert)
     DO j = 1, ngl
        area(:,j) = gridarea(j)
     ENDDO
#else
     DO j = 1, ngl
        area(:,j) = gridarea(ngl+1-j)
     ENDDO
#endif
     area_slm(:,:) = area(:,:)

     WHERE (slf_sn(:,:) > 0.0_dp .AND. mask(:,:) == 0)
        area(:,:) = (1.0_dp-slf_sn(:,:)) * area(:,:)
     END WHERE

     WHERE (mask(:,:) == 1)
        area(:,:) = 0.0_dp
     END WHERE

     WHERE (mask_slm(:,:) == 1)
        area_slm(:,:) = 0.0_dp
     END WHERE

     grdacr='atmo'
     grdacr_2='atml'   ! the mask on this grid corresponds to slm (0/1 mask)
#ifdef __synout
     WRITE(kout,*)'prism_write_grid: ', grdacr, ', ', grdacr_2
#endif     
     CALL prism_write_grid (grdacr,   nlon, ngl, lon(:,:), lat(:,:))
     CALL prism_write_grid (grdacr_2, nlon, ngl, lon(:,:), lat(:,:))

#ifdef __synout
     WRITE(kout,*)'prism_write_corner: ', grdacr
#endif     
     CALL prism_write_corner (grdacr,   nlon, ngl, nc, clon(:,:,:), clat(:,:,:))
     CALL prism_write_corner (grdacr_2, nlon, ngl, nc, clon(:,:,:), clat(:,:,:))
 
#ifdef __synout
     WRITE(kout,*)'prism_write_mask: ', grdacr
#endif     
     CALL prism_write_mask (grdacr,   nlon, ngl, mask(:,:))
     CALL prism_write_mask (grdacr_2, nlon, ngl, mask_slm(:,:))

#ifdef __synout
     WRITE(kout,*)'prism_write_area: ', grdacr
#endif     
     CALL prism_write_area (grdacr,   nlon, ngl, area(:,:))
     CALL prism_write_area (grdacr_2, nlon, ngl, area_slm(:,:))

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
   WRITE(kout,*)' ************************'
   WRITE(kout,*)' '
#endif

#endif

  RETURN

  END SUBROUTINE grids_writing

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  ini_a2o
!
! !INTERFACE:

  SUBROUTINE ini_a2o(kout)

! !RETURN VALUE:

  IMPLICIT NONE

  INTEGER, INTENT(in) :: kout      ! unit number for standard output

! !DESCRIPTION:
!
! - Initializes outgoing exchange fields.
!
! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

#if defined (__prism)

#ifdef __synout
   IF (p_parallel_io) THEN
      WRITE(kout,*)' '
      WRITE(kout,*)' Start of INI_A2O'
      WRITE(kout,*)' ****************'
      WRITE(kout,*)' '
    END IF
#endif
    awhea(:,:) = 0.0_dp
    awsol(:,:) = 0.0_dp
    awfre(:,:) = 0.0_dp
    awust(:,:) = 0.0_dp
    awvst(:,:) = 0.0_dp
    awsta(:,:) = 0.0_dp
    aicon(:,:) = 0.0_dp
    aiqre(:,:) = 0.0_dp
    aifre(:,:) = 0.0_dp
    aiust(:,:) = 0.0_dp
    aivst(:,:) = 0.0_dp
    afre_residual(:,:) = 0.0_dp
!---wiso-code
    IF (lwiso) THEN
      wisoawfre(:,:,:) = 0.0_dp
      wisoaifre(:,:,:) = 0.0_dp
      wisoafre_residual(:,:,:) = 0.0_dp
    END IF
!---wiso-code-end
    
#ifdef __synout
    IF (p_parallel_io) THEN
       WRITE(kout,*)' '
       WRITE(kout,*)' End of INI_A2O'
       WRITE(kout,*)' **************'
       WRITE(kout,*)' '
    END IF
#endif
    
#endif
    
  END SUBROUTINE ini_a2o


!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  get_a2o
!
! !INTERFACE:

  SUBROUTINE get_a2o(clfield,pfield,jpdim_atmi,jpdim_atmj,kout,kerr)

! !USES:
!---wiso-code
  USE mo_wiso,           ONLY: nwiso
!---wiso-code-end

! !RETURN VALUE:
  IMPLICIT NONE

  CHARACTER(len=8), INTENT(in)    :: clfield

  INTEGER,      INTENT(in) :: jpdim_atmi, jpdim_atmj
#if defined (__prism)
  REAL(pwp), TARGET        :: pfield(jpdim_atmi,jpdim_atmj)
#else
  REAL(dp), TARGET         :: pfield(jpdim_atmi,jpdim_atmj)
#endif
  INTEGER,      INTENT(in) :: kout
  INTEGER,      INTENT(in) :: kerr

! !DESCRIPTION:
!
! - Move model field to exchange array.
!
! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

      INTEGER :: j, index
!---wiso-code
      INTEGER :: jt
!---wiso-code

      REAL(dp), POINTER :: gl(:,:)

#if defined (__prism)

#ifdef __synout
      IF (p_parallel_io) THEN
         WRITE(kout,*) ' '
         WRITE(kout,*) ' Start of get_a2o'
         WRITE(kout,*) ' ****************'
         WRITE(kout,*) ' '
      END IF
#endif

      index = 0

      index = index + 1
!
!- zonal wind stress over water
!
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield

         CALL gather_gp(gl, awust, gd)

#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF

#if defined __cpl_mpiom
      index = index + 1
!
!- zonal wind stress over water
!
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, awust, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF
#endif

      index = index + 1
!
!- meridional wind stress over water
!
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, awvst, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF

#if defined __cpl_mpiom
      index = index + 1
!
!- meridional wind stress over water
!
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, awvst, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF
#endif

#ifndef __cpl_fluxes4
      index = index + 1
!
!- zonal wind stress over ice
!
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, aiust, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF
#endif

#if defined __cpl_mpiom
      index = index + 1
!
!- zonal wind stress over ice
!
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, aiust, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF
#endif

#ifndef __cpl_fluxes4
      index = index + 1
!
!- meridional wind stress over ice
!
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, aivst, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF
#endif

#if defined __cpl_mpiom
      index = index + 1
!
!- meridional wind stress over ice
!
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, aivst, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF
#endif

#ifndef __cpl_fluxes4
      index = index + 1
!
!- snow flux on ice
!
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, aifre, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF
#endif

      index = index + 1
!
!- water flux into the ocean
!
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, awfre, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF

#ifndef __cpl_fluxes4
!
!- residual heat flux (sea-ice topmelt heat flux)
!
      index = index + 1
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, aiqre, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF
#endif

#ifndef __cpl_fluxes4
      index = index + 1
!
!- heat flux in sea ice
!
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, aicon, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF
#endif

      index = index + 1
!
!- heat flux into the ocean
!
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, awhea, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF

      index = index + 1
!
!- solar radiation into the ocean
!
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, awsol, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF

#if defined __cpl_mpiom
!
!- 10m wind speed
!
      index = index + 1
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, awsta, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF
#endif

!---wiso-code
    IF (lwiso) THEN

! here:  pass fresh water flux (H216O, H218O, HD16O) to ocean model
!        (assume tracer oder in *wisoawfre*: 1=H216O, 2=H218O, 3=HD16O)
 
      index = index + 1
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, wisoawfre(:,1,:), gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF

      index = index + 1
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, wisoawfre(:,2,:), gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF

      index = index + 1
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, wisoawfre(:,3,:), gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF


! here:  pass snow over sea ice flux (H216O, H218O, HD16O) to ocean model
!        (assume tracer oder in *wisoaifre*: 1=H216O, 2=H218O, 3=HD16O)
 
      index = index + 1
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, wisoaifre(:,1,:), gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF

      index = index + 1
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, wisoaifre(:,2,:), gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF

      index = index + 1
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, wisoaifre(:,3,:), gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF

    END IF
!---wiso-code-end

#if defined __cpl_co2
!
!- CO2 concentration (used for diagnostics in hamocc)
!
      index = index + 1
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, co2atmos, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF
!
!- CO2 flux
!
      index = index + 1
      IF ( clfield == clstr8a2o(index) ) THEN
         gl => pfield
         CALL gather_gp(gl, co2flux_cpl, gd)
#ifdef __synout
         IF (p_parallel_io) THEN
            WRITE(kout,*) &
            ' ',clstr8a2o(index),' : ',(pfield(jpdim_atmi/2,j),j=1,jpdim_atmj,15)
         END IF
#endif
         GOTO 1000
      ENDIF

#endif

#ifdef __synout
      CALL FLUSH(kout)
#endif
      CALL couple_abort (modnam,'get_a2o','Invalid locator string in get_a2o',kerr)

 1000 CONTINUE

#ifdef __synout
      IF (p_parallel_io) THEN
         WRITE(kout,*) ' '
         WRITE(kout,*) ' End of get_a2o'
         WRITE(kout,*) ' **************'
         WRITE(kout,*) ' '
      END IF
#endif

#endif

  END SUBROUTINE get_a2o
  
!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE:  put_o2a
!
! !INTERFACE:

  SUBROUTINE put_o2a(clfield,pfield,jpdx,jpdy,kout,kerr)

! !USES:

! !RETURN VALUE:

    IMPLICIT NONE

    CHARACTER(len=8), INTENT(in) :: clfield           ! exchange field symbolic name

    INTEGER, INTENT(in)          :: jpdx,jpdy         ! exchange field dimensions
#if defined (__prism)
    REAL(pwp), INTENT(in)        :: pfield(jpdx,jpdy) ! array with received field
#else
    REAL(dp), INTENT(in)         :: pfield(jpdx,jpdy) ! array with received field
#endif
    INTEGER, INTENT(in)          :: kout              ! unit standard output
    INTEGER, INTENT(in)          :: kerr              ! error out unit

! !DESCRIPTION:
!
! - Moves the data received from OASIS to the respective 
!   fields defined on the model grid.
!   Since seaice/siced/sni are also used for lake surfaces,
!   use only those values where slf<1 and alake=0 (criterion for
!   ocean surface in grid cell).
! 
! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

    INTEGER :: j, jx, jy, index
    INTEGER :: ii1,ii2,ii3,ii4

!CSL    REAL(dp), POINTER, SAVE :: gl_slf(:,:), gl_alake(:,:), gl_tsw(:,:), &
!CSL                               gl_seaice(:,:), gl_siced(:,:), gl_sni(:,:)
!CSL
    LOGICAL, SAVE :: linit_glexflds = .TRUE. 

#if defined (__prism)

#ifdef __synout
    IF (p_parallel_io) THEN
       WRITE(kout,*) ' '
       WRITE(kout,*) ' Start of put_o2a'
       WRITE(kout,*) ' ****************'
       WRITE(kout,*) ' '
    END IF
#endif

    ! start receiving and processing

    index = 0
      
    index = index + 1

    IF ( clfield == clstr8o2a(index) ) THEN

      CALL gather_gp(gl_tsw, tsw, gd)
      IF (p_parallel_io) THEN
  
        do jy = 1, jpdy
          do jx = 1,jpdx
            if (gl_slf(jx,jy) < 1.0_dp .and. gl_alake(jx,jy) == 0.0_dp) THEN
              gl_tsw(jx,jy) = pfield(jx,jy)
            endif
          enddo
        ENDDO
        ! Do not allow unrealistic sst's on ocean cells.
        
        do jy = 1,jpdy
          do jx = 1,jpdx
            IF( gl_tsw(jx,jy) < 270.0_dp .AND. gl_slf(jx,jy) < 1.0_dp ) THEN
              WRITE(kout,*) 'WARNING !!!','slf,alake,tsw:', &
                   gl_slf(jx,jy),gl_alake(jx,jy),gl_tsw(jx,jy),jx,jy
!vg              gl_tsw(jx,jy)= 271.15_dp
!vg!!              CALL couple_abort (modnam,'put_o2a','pb with partial land/sea mask',kerr)
           ENDIF
          ENDDO
        ENDDO

#ifdef __synout
        WRITE(kout,*) ' ',clfield,' : ',(pfield(jpdx/2,j),j=1,jpdy,14)
#endif
      END IF

      CALL scatter_gp (gl_tsw, tsw, gd)
#ifdef __synout
      IF (p_parallel_io) CALL print_sst(kout)
#endif
      GOTO 1000
    ENDIF

    
#ifndef __cpl_fluxes4
    index = index + 1
    IF ( clfield == clstr8o2a(index) ) THEN

      CALL gather_gp(gl_siced, siced, gd)
      IF (p_parallel_io) THEN

        do jy = 1, jpdy
          do jx = 1,jpdx
            IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
              gl_siced(jx,jy) = pfield(jx,jy)
            END IF
          ENDDO
        ENDDO

#ifdef __synout
        WRITE(kout,*) ' ',clfield,' : ',(pfield(jpdx/2,j),j=1,jpdy,14)
#endif
      END IF

      CALL scatter_gp (gl_siced, siced, gd)
      GOTO 1000
    ENDIF



    index = index + 1
    IF ( clfield == clstr8o2a(index) ) THEN

      CALL gather_gp(gl_seaice, seaice, gd)
      IF (p_parallel_io) THEN

        do jy = 1, jpdy
          do jx = 1,jpdx
            IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
              gl_seaice(jx,jy) = pfield(jx,jy)
            END IF
          ENDDO
        ENDDO

#ifdef __synout
        WRITE(kout,*) ' ',clfield,' : ',(pfield(jpdx/2,j),j=1,jpdy,14)
#endif
      END IF

      CALL scatter_gp (gl_seaice, seaice, gd)
      GOTO 1000
    ENDIF

    index = index + 1
    IF ( clfield == clstr8o2a(index) ) THEN

      CALL gather_gp(gl_sni, sni, gd)
      IF (p_parallel_io) THEN

        do jy = 1, jpdy
          do jx = 1,jpdx
            IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
              gl_sni(jx,jy) = pfield(jx,jy)
            END IF
          ENDDO
        ENDDO
        
#ifdef __synout
        WRITE(kout,*) ' ',clfield,' : ',(pfield(jpdx/2,j),j=1,jpdy,14)
#endif
      END IF

      CALL scatter_gp(gl_sni, sni, gd)
      GOTO 1000
    ENDIF
!jj test add. fields
    index = index + 1
    IF ( clfield == clstr8o2a(index) ) THEN

      CALL gather_gp(gl_ocu, ocu, gd)
      IF (p_parallel_io) THEN

        do jy = 1, jpdy
          do jx = 1,jpdx
            IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
              gl_ocu(jx,jy) = pfield(jx,jy)
            END IF
          ENDDO
        ENDDO
        
#ifdef __synout
        WRITE(kout,*) ' ',clfield,' : ',(pfield(jpdx/2,j),j=1,jpdy,14)
#endif
      END IF

      CALL scatter_gp(gl_ocu, ocu, gd)
!       ii1=1
!       ii2=3
!       ii3=-100
!       ii4=jpdx*jpdy
!       WRITE(398) ii1,ii2,ii3,ii4
!       WRITE(398) pfield
!       WRITE(398,*) 'u comp'
!       DO jx=1,jpdx
!        WRITE(398,*) (NINT(pfield(jx,jy)),jy=1,jpdy)
!       ENDDO
      GOTO 1000
    ENDIF
    index = index + 1
    IF ( clfield == clstr8o2a(index) ) THEN

      CALL gather_gp(gl_ocv, ocv, gd)
      IF (p_parallel_io) THEN

        do jy = 1, jpdy
          do jx = 1,jpdx
            IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
              gl_ocv(jx,jy) = pfield(jx,jy)
            END IF
          ENDDO
        ENDDO
        
#ifdef __synout
        WRITE(kout,*) ' ',clfield,' : ',(pfield(jpdx/2,j),j=1,jpdy,14)
#endif
      END IF

!       ii1=1
!       ii2=4
!       ii3=-100
!       ii4=jpdx*jpdy
!       WRITE(399) ii1,ii2,ii3,ii4
!       WRITE(399) pfield
!       WRITE(398,*) 'v comp'
!       DO jx=1,jpdx
!        WRITE(398,*) (NINT(pfield(jx,jy)),jy=1,jpdy)
!       ENDDO
!       CLOSE(398)
      CALL scatter_gp(gl_ocv, ocv, gd)
      GOTO 1000
    ENDIF

!wiso-code
    IF (lwiso) THEN

! here:  pass delta value of ocean surface water (H216O, H218O, HD16O) to atmosphere model
!        (assume tracer oder in *wisosw_d*: 1=H216O, 2=H218O, 3=HD16O)
 
    index = index + 1
    IF ( clfield == clstr8o2a(index) ) THEN

      CALL gather_gp(gl_o16, wisosw_d(:,1,:), gd)
      IF (p_parallel_io) THEN

        do jy = 1, jpdy
          do jx = 1,jpdx
            IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
              gl_o16(jx,jy) = pfield(jx,jy)
            END IF
          ENDDO
        ENDDO

#ifdef __synout
        WRITE(kout,*) ' ',clfield,' : ',(pfield(jpdx/2,j),j=1,jpdy,14)
#endif
      END IF

      CALL scatter_gp (gl_o16, wisosw_d(:,1,:), gd)
      GOTO 1000
    ENDIF

    index = index + 1
    IF ( clfield == clstr8o2a(index) ) THEN

      CALL gather_gp(gl_o18, wisosw_d(:,2,:), gd)
      IF (p_parallel_io) THEN

        do jy = 1, jpdy
          do jx = 1,jpdx
            IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
              gl_o18(jx,jy) = pfield(jx,jy)
            END IF
          ENDDO
        ENDDO

#ifdef __synout
        WRITE(kout,*) ' ',clfield,' : ',(pfield(jpdx/2,j),j=1,jpdy,14)
#endif
      END IF

      CALL scatter_gp (gl_o18, wisosw_d(:,2,:), gd)
      GOTO 1000
    ENDIF

    index = index + 1
    IF ( clfield == clstr8o2a(index) ) THEN

      CALL gather_gp(gl_hdo, wisosw_d(:,3,:), gd)
      IF (p_parallel_io) THEN

        do jy = 1, jpdy
          do jx = 1,jpdx
            IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
              gl_hdo(jx,jy) = pfield(jx,jy)
            END IF
          ENDDO
        ENDDO

#ifdef __synout
        WRITE(kout,*) ' ',clfield,' : ',(pfield(jpdx/2,j),j=1,jpdy,14)
#endif
      END IF

      CALL scatter_gp (gl_hdo, wisosw_d(:,3,:), gd)
      GOTO 1000
    ENDIF

    END IF
!wiso-code-end

#ifdef __cpl_co2
    index = index + 1
    IF ( clfield == clstr8o2a(index) ) THEN

      CALL gather_gp(gl_co2trans, co2trans, gd)
      IF (p_parallel_io) THEN

        do jy = 1, jpdy
          do jx = 1,jpdx
            IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
              gl_co2trans(jx,jy) = pfield(jx,jy) * 1.e6
            END IF
          ENDDO
        ENDDO
        
#ifdef __synout
        WRITE(kout,*) ' ',clfield,' : ',(pfield(jpdx/2,j),j=1,jpdy,14)
#endif
      END IF

      CALL scatter_gp(gl_co2trans, co2trans, gd)
      GOTO 1000
    ENDIF
    index = index + 1
    IF ( clfield == clstr8o2a(index) ) THEN

      CALL gather_gp(gl_co2ocean, co2ocean, gd)
      IF (p_parallel_io) THEN

        do jy = 1, jpdy
          do jx = 1,jpdx
            IF(gl_slf(jx,jy) < 1.0_dp .AND. gl_alake(jx,jy) == 0.0_dp) THEN
              gl_co2ocean(jx,jy) = pfield(jx,jy)
            END IF
          ENDDO
        ENDDO
        
#ifdef __synout
        WRITE(kout,*) ' ',clfield,' : ',(pfield(jpdx/2,j),j=1,jpdy,14)
#endif
      END IF

      CALL scatter_gp(gl_co2ocean, co2ocean, gd)
      GOTO 1000
    ENDIF

#endif

#endif
         CALL couple_abort (modnam,'put_o2a','Invalid locator string in put_o2a',kerr)

1000 CONTINUE

#ifdef __synout
    IF (p_parallel_io) THEN    
         WRITE(kout,*) ' '
         WRITE(kout,*) ' End of put_o2a'
         WRITE(kout,*) ' **************'
         WRITE(kout,*) ' '
         CALL flush(kout)
      ENDIF
#endif

#endif

  END SUBROUTINE put_o2a


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
  INTEGER,  INTENT(in)   :: kout  ! unit std output
  REAL(dp), INTENT(in)   :: pdt   ! model time step length (secs)
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

    !
    !-- Accumulation of time passed in this run.
    !
    !  
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

#ifdef DEBUG

#ifndef NOMPI
    INCLUDE 'mpif.h'

    INTEGER :: ierr           ! error code

    WRITE (kerr,'(a,a,a,a,a)') 'Calling MPI_ABORT in ',str2,' by model ',str1,':'
    WRITE (kerr,'(a,a)')   '        ',str3
    
    CALL MPI_ABORT (MPI_COMM_WORLD, 0, ierr)

    IF (ierr /= MPI_SUCCESS) THEN
      WRITE (kerr,'(a)') ' MPI_ABORT failed'
      WRITE (kerr,'(a,i4)') ' Error =  ', ierr
      STOP
    END IF
    CALL FLUSH(kerr)
#endif

#else

#if defined (__prism)
    EXTERNAL :: prism_abort_proto

    CALL prism_abort_proto (comp_id, TRIM(str2), TRIM(str3))

#endif

#endif

      END SUBROUTINE couple_abort

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

#if defined (__prism)

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
        CALL couple_abort (modnam,'couple_get','pb prism_get',kerr)
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

#endif

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

#if defined (__prism)

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

#endif

  END SUBROUTINE digest_put_Id
!----------------------------------------------------------------------

  SUBROUTINE water_budget_corr

    !  freshwater (p-e and residual flux) correction applied to ice-free ocean points

    !  Author:
    !  U. Schlese, MPI, July 2003
    !  E. Roeckner, U. Schlese, MPI, February 2007
    !  T. Raddatz, MPI, November 2007

    USE mo_control,       ONLY: ngl, nlon
    USE mo_memory_g3b,    ONLY: apmebco, rain, awfre, afre_residual, slf, slm, alake, seaice
    USE mo_decomposition, ONLY: gl_dc => global_decomposition
    USE mo_transpose,     ONLY: gather_gp, scatter_gp
    USE mo_gaussgrid,     ONLY: gl_budw
    USE mo_constants,     ONLY: rhoh2o
    USE mo_mpi,           ONLY: p_pe, p_io
!---wiso-code
    USE mo_memory_wiso,   ONLY: wisoapmebco, wisorain, wisoawfre, wisoafre_residual
    USE mo_wiso,          ONLY: nwiso
!---wiso-code-end

    IMPLICIT NONE

#if defined (__prism)

    ! local arrays

    REAL(dp), POINTER  :: zpmeb(:,:)   ! p-e correction, accumulated [kg/m**2]
    REAL(dp), POINTER  :: zrain(:,:)   ! total rain, accumulated     [kg/m**2]
    REAL(dp), POINTER  :: zawfre(:,:)  ! liquid freshwater flux [m/s]
    REAL(dp), POINTER  :: zafre_residual(:,:)  ! liquid + solid freshwater flux before discharge [m/s]
    REAL(dp), POINTER  :: zslf(:,:)    ! fractional land cover
    REAL(dp), POINTER  :: zslm(:,:)    ! echam land-sea mask (1. for land and 0. for ocean)
    REAL(dp), POINTER  :: zalake(:,:)  ! lake mask [0,1]
    REAL(dp), POINTER  :: zseaice(:,:) ! fraction of water covered by sea ice
    REAL(dp), POINTER  :: zwater(:,:)  ! fraction of grid area covered by ice-free ocean

    REAL(dp)           :: zpmebz(ngl), zrainz(ngl), zresiz(ngl)

!---wiso-code
    REAL(dp), POINTER :: zdummy(:,:)
    REAL(dp), TARGET  :: zwisopmeb(nlon,nwiso,ngl)   ! p-e correction, accumulated [kg/m**2] - water isotopes
    REAL(dp), TARGET  :: zwisorain(nlon,nwiso,ngl)   ! total rain, accumulated     [kg/m**2] - water isotopes
    REAL(dp), TARGET  :: zwisoawfre(nlon,nwiso,ngl)  ! liquid freshwater flux [m/s] - water isotopes
    REAL(dp), TARGET  :: zwisoafre_residual(nlon,nwiso,ngl)  ! liquid + solid freshwater fluxbefore discharge [m/s] - water isotopes

    REAL(dp)           :: zwisopmebz(nwiso,ngl), zwisorainz(nwiso,ngl), zwisoresiz(nwiso,ngl)
!---wiso-code-end

    ! Local scalars

    INTEGER :: jlat
    REAL(dp)    :: zpmeb_glob, zrain_glob, zresi_glob

!---wiso-code
    INTEGER :: jt
    REAL(dp)    :: zwisopmeb_glob(nwiso), zwisorain_glob(nwiso), zwisoresi_glob(nwiso), zdelta
!---wiso-code-end

    IF (p_pe == p_io) THEN
      ALLOCATE (zpmeb(nlon,ngl))
      ALLOCATE (zrain(nlon,ngl))
      ALLOCATE (zawfre(nlon,ngl))
      ALLOCATE (zafre_residual(nlon,ngl))
      ALLOCATE (zslf(nlon,ngl))
      ALLOCATE (zslm(nlon,ngl))
      ALLOCATE (zalake(nlon,ngl))
      ALLOCATE (zseaice(nlon,ngl))
      ALLOCATE (zwater(nlon,ngl))
    END IF

!---wiso-code
    IF (lwiso) THEN
      IF (p_pe == p_io) THEN
        ALLOCATE (zdummy(nlon,ngl))
      END IF
    END IF
!---wiso-code-end

    CALL gather_gp (zpmeb, apmebco,gl_dc)
    CALL gather_gp (zrain, rain,   gl_dc)
    CALL gather_gp (zawfre,awfre,  gl_dc)
    CALL gather_gp (zafre_residual,afre_residual,gl_dc)
    CALL gather_gp (zslf,  slf,    gl_dc)
    CALL gather_gp (zslm,  slm,    gl_dc)
    CALL gather_gp (zalake,alake,  gl_dc)
    CALL gather_gp (zseaice,seaice,gl_dc)

!---wiso-code
    IF (lwiso) THEN
      DO jt=1,nwiso
        zdummy => zwisopmeb(:,jt,:)
        CALL gather_gp (zdummy, wisoapmebco(:,jt,:),gl_dc)
        zdummy => zwisorain(:,jt,:)
        CALL gather_gp (zdummy, wisorain(:,jt,:),   gl_dc)
        zdummy => zwisoawfre(:,jt,:)
        CALL gather_gp (zdummy,wisoawfre(:,jt,:),  gl_dc)
        zdummy => zwisoafre_residual(:,jt,:)
        CALL gather_gp (zdummy,wisoafre_residual(:,jt,:),gl_dc)
      END DO
    END IF
!---wiso-code-end

    IF (p_pe == p_io) THEN

!  ice-free ocean fractional area

      zwater(:,:) = (1._dp - zslf(:,:)) * (1._dp - zseaice(:,:)) * (1._dp - zalake(:,:))

!  total rain and p-e correction [kg/m2] averaged over coupling time step [m/s]

      zrain(:,:) = zrain(:,:)/(rhoh2o*couple_a2o_time)
      zpmeb(:,:) = zpmeb(:,:)/(rhoh2o*couple_a2o_time)

!---wiso-code
      IF (lwiso) THEN
        zwisorain(:,:,:) = zwisorain(:,:,:)/(rhoh2o*couple_a2o_time)
        zwisopmeb(:,:,:) = zwisopmeb(:,:,:)/(rhoh2o*couple_a2o_time)
      END IF
!---wiso-code-end

!  latitudinal sum of rain, p-e correction and residual flux correction [m/s over the whole globe]

      DO jlat = 1,ngl
        zrainz(jlat) = SUM(zrain(:,jlat) * zwater(:,jlat)) * gl_budw(jlat)
        zpmebz(jlat) = SUM(zpmeb(:,jlat)) * gl_budw(jlat)
        zresiz(jlat) = SUM(zafre_residual(:,jlat) * (zslf(:,jlat) - zslm(:,jlat)) * &
                           (1._dp - zalake(:,jlat))) * gl_budw(jlat)
      END DO

!---wiso-code
      IF (lwiso) THEN
        DO jt = 1,nwiso
          DO jlat = 1,ngl
            zwisorainz(jt,jlat) = SUM(zwisorain(:,jt,jlat) * zwater(:,jlat)) * gl_budw(jlat)
            zwisopmebz(jt,jlat) = SUM(zwisopmeb(:,jt,jlat)) * gl_budw(jlat)
            zwisoresiz(jt,jlat) = SUM(zwisoafre_residual(:,jt,jlat) * (zslf(:,jlat) - zslm(:,jlat)) * &
                                     (1._dp - zalake(:,jlat))) * gl_budw(jlat)
          END DO
        END DO
      END IF
!---wiso-code-end

!  global sum of rain over ice-free ocean, p-e correction and residual flux correction [m/s over the whole globe]

      zrain_glob = SUM(zrainz(:))
      zpmeb_glob = SUM(zpmebz(:))
      zresi_glob = SUM(zresiz(:))

!---wiso-code
      IF (lwiso) THEN
        DO jt = 1,nwiso
          zwisorain_glob(jt) = SUM(zwisorainz(jt,:))
          zwisopmeb_glob(jt) = SUM(zwisopmebz(jt,:))
          zwisoresi_glob(jt) = SUM(zwisoresiz(jt,:))
        END DO
      END IF
!---wiso-code-end

!  p-e correction and residual flux correction applied to the ocean freshwater flux [m/s]

      zwater(:,:) = (1._dp - zseaice(:,:)) * (1._dp - zalake(:,:))
      zawfre(:,:) = zawfre(:,:) + zrain(:,:) * zwater(:,:) * (zpmeb_glob + zresi_glob) / zrain_glob

!---wiso-code
      IF (lwiso) THEN
        DO jt = 1,nwiso
          zwisoawfre(:,jt,:) = zwisoawfre(:,jt,:) + zwisorain(:,jt,:) * zwater(:,:) * (zwisopmeb_glob(jt) + zwisoresi_glob(jt)) &
                                                                                      / zwisorain_glob(jt)
        END DO
      END IF
!---wiso-code-end

    END IF

    CALL scatter_gp (zawfre, awfre, gl_dc)

!---wiso-code
    IF (lwiso) THEN
      DO jt=1,nwiso
        zdummy(:,:) = zwisoawfre(:,jt,:)
        CALL scatter_gp (zdummy, wisoawfre(:,jt,:), gl_dc)
      END DO
    END IF
!---wiso-code-end

    IF (p_pe == p_io) THEN
      DEALLOCATE (zpmeb)
      DEALLOCATE (zawfre)
      DEALLOCATE (zafre_residual)
      DEALLOCATE (zrain)
      DEALLOCATE (zslf)
      DEALLOCATE (zslm)
      DEALLOCATE (zalake)
      DEALLOCATE (zseaice)
      DEALLOCATE (zwater)
    END IF
#endif

    RETURN
  END SUBROUTINE water_budget_corr
END MODULE mo_couple

