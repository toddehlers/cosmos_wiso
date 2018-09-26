MODULE mo_sub_nml
!------------------------------------------------------------------------------
! namelist interface to the tracer and output stream facilities
!
! Authors:
! Andreas Rhodin MPI/DWD April 2002
!------------------------------------------------------------------------------
  !-------------
  ! Modules used
  !-------------
  USE mo_kind,        ONLY: dp               ! working precition kind parameter
  USE mo_tracer,      ONLY: trlist,         &! tracer list
                            jptrac,         &! index of last (unused) tracer
                            AUTO,           &! take arbitrary GRIB code number
                            n_t =>new_tracer ! tracer request routine
  USE mo_tracdef,     ONLY: ln               ! length of name-string
  USE mo_mpi,         ONLY: p_pe,           &! processor index of actual PE
                            p_io,           &! I/O processor index
                            p_parallel_io,  &! flag for p_pe==p_io
                            p_parallel,     &! flag for parallel processing
                            p_bcast          ! broadcast routine
  USE mo_namelist,    ONLY: position_nml,   &! function: search namelist group
                            POSITIONED,     &! return value: OK
                            nnml             ! unit number of namelist file
  USE mo_exception,   ONLY: finish
  USE mo_time_event,  ONLY: io_time_event
  USE mo_time_control,ONLY: get_new_ev_putdata
  USE mo_memory_base, ONLY: memory_info, t_stream, nstreams,ostreams,&
                            set_stream_element_info, get_stream_element_info
  IMPLICIT NONE

  !----------------
  ! Public entities
  !----------------
  PRIVATE
  PUBLIC :: request_tracer_nml      ! request tracers via namelist
  PUBLIC :: set_tracer_nml          ! set tracer meta data via namelist
  PUBLIC :: set_stream_element_nml  ! set meta data from namelist
  PUBLIC :: set_stream_nml          ! set stream parameters
  
  CONTAINS
!==============================================================================
  SUBROUTINE request_tracer_nml
  !============================================================
  ! request tracers as specified by namelist group /NEW_TRACER/
  !============================================================
    !--------------------------
    ! namelist group definition
    !--------------------------
    CHARACTER (len=ln) :: name   ! name of tracer
    CHARACTER (len=ln) :: units  ! units
    INTEGER            :: ninit  ! initialisation flag
    INTEGER            :: nrerun ! restart flag
    REAL(dp)           :: vini   ! initialization value
    REAL(dp)           :: tdecay ! exponential decay time
    INTEGER            :: nwrite ! write to postprocessing file (0 or 1)
    INTEGER            :: code   ! GRIB code  number
    INTEGER            :: table  ! GRIB table number
    INTEGER            :: bits   ! number of bits used for GRIB encoding
    INTEGER            :: ntran  ! transport flag
    INTEGER            :: nvdiff ! vertical diffusion flag
    INTEGER            :: nconv  ! convection flag
    INTEGER            :: nint   ! integration flag

    NAMELIST /NEW_TRACER/ name, units, ninit, vini, nrerun, tdecay, nwrite, &
                          code, table, bits, ntran, nvdiff, nconv, nint
    !----------------
    ! local variables
    !----------------
    LOGICAL :: first  ! true for first attempt to position namelist group
    INTEGER :: ierr   ! error return value

    !----------------------------------------------------
    ! loop over occurences of namelist group /NEW_TRACER/
    !----------------------------------------------------
    first = .TRUE.
    DO
      IF (p_pe==p_io) THEN
        CALL position_nml ('NEW_TRACER', REWIND=first, status=ierr)
        first = .FALSE.
        !------------------------------
        ! if namelist group is present:
        !------------------------------
        SELECT CASE (ierr)
        CASE (POSITIONED)
          !----------------------------------------------------
          ! set default values from unused entry in tracer list
          !----------------------------------------------------
          name     = ''
          units    = trlist% ti(jptrac)% units
          ninit    = trlist% ti(jptrac)% ninit
          nrerun   = trlist% ti(jptrac)% nrerun
          vini     = trlist% ti(jptrac)% vini
          tdecay   = trlist% ti(jptrac)% tdecay
          nwrite   = trlist% ti(jptrac)% nwrite
          code     = AUTO
          table    = trlist% ti(jptrac)% table
          bits     = trlist% ti(jptrac)% gribbits
          ntran    = trlist% ti(jptrac)% ntran
          nvdiff   = trlist% ti(jptrac)% nvdiff
          nconv    = trlist% ti(jptrac)% nconv
          nint     = trlist% ti(jptrac)% nint
          !--------------
          ! read namelist
          !--------------
          READ (nnml, new_tracer)
        END SELECT
      ENDIF
      !-----------------------------------------
      ! broadcast parameters to other processors
      !-----------------------------------------
      IF (p_parallel) THEN
        CALL    p_bcast (ierr,  p_io)
        IF (ierr==POSITIONED) THEN
          CALL p_bcast (name,   p_io)
          CALL p_bcast (units,  p_io)
          CALL p_bcast (ninit,  p_io)
          CALL p_bcast (nrerun, p_io)
          CALL p_bcast (vini,   p_io)
          CALL p_bcast (tdecay, p_io)
          CALL p_bcast (nwrite, p_io)
          CALL p_bcast (code,   p_io)
          CALL p_bcast (table,  p_io)
          CALL p_bcast (bits,   p_io)
          CALL p_bcast (ntran,  p_io)
          CALL p_bcast (nvdiff, p_io)
          CALL p_bcast (nconv,  p_io)
          CALL p_bcast (nint,   p_io)
        ENDIF
      ENDIF
      IF (ierr/=POSITIONED) EXIT
      !-------------------
      ! request new tracer
      !-------------------
      CALL n_t (name       = name,       &
                modulename = 'namelist', &
                units      = units,      &
                ninit      = ninit,      &
                vini       = vini,       &
                nrerun     = nrerun,     &
                tdecay     = tdecay,     &
                nwrite     = nwrite,     &
                code       = code,       &
                table      = table,      &
                bits       = bits,       &
                ntran      = ntran,      &
                nvdiff     = nvdiff,     &
                nconv      = nconv,      &
                nint       = nint        )
    END DO

  END SUBROUTINE request_tracer_nml
!------------------------------------------------------------------------------
  SUBROUTINE set_tracer_nml
  !============================================================
  ! request tracers as specified by namelist group /SET_TRACER/
  !============================================================
    !--------------------------
    ! namelist group definition
    !--------------------------
    CHARACTER (len=ln) :: name   ! name of tracer
    CHARACTER (len=ln) :: units  ! units
    INTEGER            :: ninit  ! initialisation flag
    INTEGER            :: nrerun ! restart flag
    REAL(dp)           :: vini   ! initialization value
    REAL(dp)           :: tdecay ! exponential decay time
    INTEGER            :: nwrite ! write to postprocessing file (0 or 1)
    INTEGER            :: code   ! GRIB code  number
    INTEGER            :: table  ! GRIB table number
    INTEGER            :: bits   ! number of bits used for GRIB encoding
    INTEGER            :: ntran  ! transport flag
    INTEGER            :: nvdiff ! vertical diffusion flag
    INTEGER            :: nconv  ! convection flag
    INTEGER            :: nint   ! integration flag

    NAMELIST /SET_TRACER/ name, units, ninit, vini, nrerun, tdecay, nwrite, &
                          code, table, bits, ntran, nvdiff, nconv, nint
    !----------------
    ! local variables
    !----------------
    LOGICAL :: first  ! true for first attempt to position namelist group
    INTEGER :: ierr   ! error return value
    INTEGER :: it     ! tracer index
    !----------------------------------------------------
    ! loop over occurences of namelist group /SET_TRACER/
    !----------------------------------------------------
    first = .TRUE.
l1: DO
      IF (p_pe==p_io) THEN
        CALL position_nml ('SET_TRACER', REWIND=first, status=ierr)
        first = .FALSE.
        !------------------------------
        ! if namelist group is present:
        !------------------------------
        SELECT CASE (ierr)
        CASE (POSITIONED)
          !----------------------------------------------------
          ! set default values from unused entry in tracer list
          !----------------------------------------------------
          name     = ''
          units    = '999'
          ninit    = -999
          nrerun   = -999
          vini     = -999._dp
          tdecay   = -999._dp
          nwrite   = -999
          code     = -999
          table    = -999
          bits     = -999
          ntran    = -999
          nvdiff   = -999
          nconv    = -999
          nint     = -999
          !--------------
          ! read namelist
          !--------------
          READ (nnml, set_tracer)
        END SELECT
      ENDIF
      !-----------------------------------------
      ! broadcast parameters to other processors
      !-----------------------------------------
      IF (p_parallel) THEN
        CALL    p_bcast (ierr,  p_io)
        IF (ierr==POSITIONED) THEN
          CALL p_bcast (name,   p_io)
          CALL p_bcast (units,  p_io)
          CALL p_bcast (ninit,  p_io)
          CALL p_bcast (nrerun, p_io)
          CALL p_bcast (vini,   p_io)
          CALL p_bcast (tdecay, p_io)
          CALL p_bcast (nwrite, p_io)
          CALL p_bcast (code,   p_io)
          CALL p_bcast (table,  p_io)
          CALL p_bcast (bits,   p_io)
          CALL p_bcast (ntran,  p_io)
          CALL p_bcast (nvdiff, p_io)
          CALL p_bcast (nconv,  p_io)
          CALL p_bcast (nint,   p_io)
        ENDIF
      ENDIF
      IF (ierr/=POSITIONED) EXIT
      !--------------------
      ! find tracer by name
      !--------------------
      do it = 1, trlist% ntrac
        if (name == trlist% ti(it)% fullname) then
          if (units  /= '999')    trlist% ti(jptrac)% units    = units 
          if (ninit  /= -999 )    trlist% ti(jptrac)% ninit    = ninit
          if (nrerun /= -999 )    trlist% ti(jptrac)% nrerun   = nrerun
          if (vini   /= -999._dp) trlist% ti(jptrac)% vini     = vini
          if (tdecay /= -999._dp) trlist% ti(jptrac)% tdecay   = tdecay
          if (nwrite /= -999 )    trlist% ti(jptrac)% nwrite   = nwrite
          if (code   /= -999 )    trlist% ti(jptrac)% code     = code
          if (table  /= -999 )    trlist% ti(jptrac)% table    = table
          if (bits   /= -999 )    trlist% ti(jptrac)% gribbits = bits
          if (ntran  /= -999 )    trlist% ti(jptrac)% ntran    = ntran
          if (nvdiff /= -999 )    trlist% ti(jptrac)% nvdiff   = nvdiff
          if (nconv  /= -999 )    trlist% ti(jptrac)% nconv    = nconv
          if (nint   /= -999 )    trlist% ti(jptrac)% nint     = nint
          cycle l1
        end if
      end do
      call finish('set_tracer_nml','unknown tracer name: '//name)
    END DO l1

  END SUBROUTINE set_tracer_nml
!==============================================================================
  SUBROUTINE set_stream_element_nml
  USE mo_namelist,     ONLY: position_nml, nnml, POSITIONED

    CHARACTER(len=64) :: stream, name, longname ,units
    INTEGER           :: table, code, bits
    INTEGER           :: lpost, lrerun

    NAMELIST /SET_STREAM_ELEMENT/ stream, name, longname ,units, &
                                  table, code, bits,             &
                                  lpost ,lrerun

    INTEGER                    :: i
    INTEGER                    :: ierr
    LOGICAL                    :: first
    TYPE(memory_info)          :: info
    TYPE(t_stream)    ,POINTER :: strm

    first = .TRUE.
    DO
      IF (p_parallel_io) THEN
        CALL position_nml ('SET_STREAM_ELEMENT', REWIND=first, status=ierr)
        first = .FALSE.
        SELECT CASE (ierr)
        CASE (POSITIONED)
          stream   = ''
          name     = ''
          longname = ''
          units    = ''
          table    = -HUGE(table)
          code     = -HUGE(code)
          bits     = -HUGE(bits)
          lpost    = -HUGE(lpost)
          lrerun   = -HUGE(lrerun)
          READ (nnml,set_stream_element)
        END SELECT
      ENDIF
      IF (p_parallel) THEN
        CALL p_bcast (ierr,     p_io)
        CALL p_bcast (stream,   p_io)
        CALL p_bcast (name,     p_io)
        CALL p_bcast (longname, p_io)
        CALL p_bcast (units,    p_io)
        CALL p_bcast (table,    p_io)
        CALL p_bcast (code,     p_io)
        CALL p_bcast (bits,     p_io)
        CALL p_bcast (lpost,    p_io)
        CALL p_bcast (lrerun,   p_io)
      ENDIF
      IF (ierr /= POSITIONED) EXIT

      DO i=1,nstreams
        IF (stream == '' .OR. stream == ostreams(i)% name) THEN
          strm => ostreams(i)
          CALL get_stream_element_info (strm, name, info)
          IF (info% name /= '') EXIT
        ENDIF
      END DO

      IF (info% name == '') &
        CALL finish ('set_stream_element_nml','cannot find '//name)

      IF (lpost  >= 0) info% lpost  = (lpost  == 1)
      IF (lrerun >= 0) info% lrerun = (lrerun == 1)

      CALL set_stream_element_info (strm, name, longname = longname,   &
                                                units    = units,      &
                                                table    = table,      &
                                                code     = code,       &
                                                bits     = bits,       &
                                                lpost    = info% lpost,&
                                                lrerun   = info% lrerun)
    END DO

  END SUBROUTINE set_stream_element_nml
!------------------------------------------------------------------------------
  SUBROUTINE set_stream_nml
  USE mo_namelist,     ONLY: position_nml, nnml, POSITIONED
  
    CHARACTER(len=16)   :: stream   ! output stream to change
    INTEGER             :: filetype ! 'GRIB' or 'NetCDF'
    CHARACTER(len=8)    :: post_suf ! suffix of output  file
    CHARACTER(len=8)    :: rest_suf ! suffix of restart file
    CHARACTER(len=8)    :: init_suf ! suffix of initial file
    INTEGER             :: lpost    ! in standard output file
    INTEGER             :: lrerun   ! in standard restartfile
    INTEGER             :: linit    ! in standard initialfile
    TYPE(io_time_event) :: interval ! output interval

    NAMELIST /SET_STREAM/ stream, filetype, post_suf, rest_suf, init_suf, &
                          lpost, lrerun, linit, interval

    INTEGER                    :: i
    INTEGER                    :: ierr
    LOGICAL                    :: first
    TYPE(t_stream)    ,POINTER :: strm

    first = .TRUE.
    DO
      IF (p_parallel_io) THEN
        CALL position_nml ('SET_STREAM', REWIND=first, status=ierr)
        first = .FALSE.
        SELECT CASE (ierr)
        CASE (POSITIONED)
          stream   = ''
          filetype = 0
          post_suf = ''
          rest_suf = ''
          init_suf = ''
          lpost    = -1
          lrerun   = -1
          linit    = -1
          interval = io_time_event (0,'','',0)
          READ (nnml,set_stream)
        END SELECT
      ENDIF
      IF (p_parallel) THEN
        CALL   p_bcast (ierr,                 p_io)
        IF (ierr==POSITIONED) THEN
          CALL p_bcast (stream,               p_io)
          CALL p_bcast (filetype,             p_io)
          CALL p_bcast (post_suf,             p_io)
          CALL p_bcast (rest_suf,             p_io)
          CALL p_bcast (init_suf,             p_io)
          CALL p_bcast (lpost,                p_io)
          CALL p_bcast (lrerun,               p_io)
          CALL p_bcast (linit,                p_io)
          CALL p_bcast (interval% counter,    p_io)
          CALL p_bcast (interval% unit,       p_io)
          CALL p_bcast (interval% adjustment, p_io)
          CALL p_bcast (interval% offset,     p_io)
        ENDIF
      ENDIF
      IF (ierr /= POSITIONED) EXIT

      DO i=1,nstreams
        IF (stream == ostreams(i)% name) THEN
          strm => ostreams(i)
          IF (filetype          /= 0) strm% filetype = filetype
          IF (post_suf          /='') strm% post_suf = post_suf
          IF (rest_suf          /='') strm% rest_suf = rest_suf          
          IF (init_suf          /='') strm% init_suf = init_suf
          IF (lpost             /=-1) strm% lpost    = lpost ==1
          IF (lrerun            /=-1) strm% lrerun   = lrerun==1
          IF (linit             /=-1) strm% linit    = linit ==1
          IF (interval% counter /= 0) strm% post_idx = &
                        get_new_ev_putdata (interval)
          EXIT
        ENDIF
      END DO

    END DO

  END SUBROUTINE set_stream_nml
!==============================================================================
END MODULE mo_sub_nml
