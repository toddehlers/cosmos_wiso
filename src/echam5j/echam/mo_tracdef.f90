MODULE mo_tracdef

!
! definition of tracer information data type and variable
!

!------------------------------------------------------------------------------
  !
  ! Modules used
  !
  USE mo_kind,           ONLY: dp
  USE mo_parameters,     ONLY: jptrac      ! maximum number of progn. tracers
  USE mo_linked_list,    ONLY: memory_info ! I/O meta information data type
  USE mo_time_conversion,ONLY: time_days   ! date (days,seconds) data type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: trlist                         ! tracer info list variable
  PUBLIC :: t_trlist, t_trinfo, t_p_mi     ! "" type definitions
  PUBLIC :: t_flag, time_days              !    type definition
  PUBLIC :: memory_info                    !    type definition
  PUBLIC :: ln, ll, lf, nf                 ! length of names in trlist
  PUBLIC :: submlist                       ! submodel list
  PUBLIC :: t_submlist                     ! submodel list element data type
  PUBLIC :: nsubm                          ! number of submodel list entries

!------------------------------------------------------------------------------
  !
  ! Type declarations
  !

  !
  ! Individual settings for each tracer
  !
  INTEGER, PARAMETER :: ln = 24 ! length of name (char) components
  INTEGER, PARAMETER :: ll = 64 ! length of longname    component
  INTEGER, PARAMETER :: lf =  8 ! length of flag character string
  INTEGER, PARAMETER :: nf = 10 ! number of user defined flags
  INTEGER, PARAMETER :: ns = 20 ! max number of submodels

  TYPE t_flag
    CHARACTER(len=lf) :: c      ! character string
    REAL(dp)          :: v      ! value
  END TYPE t_flag

  TYPE t_trinfo
    !
    ! identification of transported quantity
    !
    CHARACTER(len=ln) :: basename   ! name (instead of xt..)
    CHARACTER(len=ln) :: subname    ! optional for 'colored' tracer
    CHARACTER(len=ln) :: fullname   ! name_subname
    CHARACTER(len=ln) :: modulename ! name of requesting sub-model
    CHARACTER(len=ln) :: units      ! units
    CHARACTER(len=ll) :: longname   ! long name
    INTEGER           :: tag        ! tag of requesting routine
    !
    ! Requested resources ...
    !
    INTEGER           :: nbudget    ! calculate budgets        (default 0)
    INTEGER           :: ntran      ! perform transport        (default 1)
    INTEGER           :: nfixtyp    ! type of mass fixer       (default 1)
    INTEGER           :: nvdiff     ! vertical diffusion flag  (default 1)
    INTEGER           :: nconv      ! convection flag          (default 1)
    INTEGER           :: nwetdep    ! wet deposition flag      (default 0)
    INTEGER           :: nsedi      ! sedimentation flag       (default 0)
    REAL(dp)          :: tdecay     ! decay time (exponential) (default 0.sec)
    INTEGER           :: nsemis     ! surface emission flag    (default 0)
    INTEGER           :: nint       ! integration flag         (default 1) 
    LOGICAL           :: lexpdry    ! explicit coupling       (default .false.)
    !
    ! initialization and restart
    !
    INTEGER           :: ninit      ! initialization request flag
    INTEGER           :: nrerun     ! rerun flag
    REAL(dp)          :: vini       ! initialisation value     (default 0.)
    !
    ! Flags used for postprocessing
    !
    INTEGER           :: nwrite     ! write flag            (default 1)
    INTEGER           :: code       ! tracer code,          (default 235...)
    INTEGER           :: table      ! tracer code table     (default 0)
    INTEGER           :: gribbits   ! bits for encoding     (default 16)
    !
    ! Flags to be used by chemistry or tracer modules
    !
    INTEGER           :: ndrydep    ! dry deposition flag   (default 0)
    REAL(dp)          :: henry      ! Henry coefficient [?] (default 1.e-10)
    REAL(dp)          :: dryreac    ! dry reactivity coeff. (default 0.)
    INTEGER           :: nsoluble   ! soluble flag          (default 0)
    INTEGER           :: nphase     ! phase indicator       (default 0)
    INTEGER           :: mode       ! mode                  (default 0)
    REAL(dp)          :: moleweight ! molecular weight      (default 0.)
    TYPE(t_flag)      :: myflag (nf)! user defined flag
    TYPE(time_days)   :: tupdatel   ! last update time
    TYPE(time_days)   :: tupdaten   ! next update time
    !
    ! Indicate actions actually performed by ECHAM
    !
    INTEGER           :: init       ! initialisation method used
  END TYPE t_trinfo

  !
  ! Reference to memory buffer information for each tracer
  ! used to access the 'restart' flags
  !
  TYPE t_p_mi                           ! pointers to memory info type
    TYPE (memory_info), POINTER :: xt   ! tracers         ,meta information
    TYPE (memory_info), POINTER :: xtm1 ! tracers at t-dt ,meta information
  END TYPE t_p_mi

  !
  ! Basic data type definition for tracer info list
  !
  TYPE t_trlist
    !
    ! global tracer list information
    !
    INTEGER         :: ntrac        ! number of tracers specified
    INTEGER         :: anyfixtyp    ! mass fixer types used
    INTEGER         :: anywetdep    ! wet deposition requested for any tracer
    INTEGER         :: anydrydep    ! wet deposition requested for any tracer
    INTEGER         :: anysedi      ! sedimentation  requested for any tracer
    INTEGER         :: anysemis     ! surface emission flag for any tracer
    INTEGER         :: anyconv      ! convection flag
    INTEGER         :: anyvdiff     ! vertical diffusion flag
    INTEGER         :: nadvec       ! number of advected tracers
    LOGICAL         :: oldrestart   ! true to read old restart format
    !
    ! individual information for each tracer
    !
    TYPE (t_trinfo) :: ti  (jptrac) ! Individual settings for each tracer
    !
    ! reference to memory buffer info
    !
    TYPE (t_p_mi)   :: mi  (jptrac) ! memory buffer information for each tracer

    TYPE (memory_info), POINTER :: mixt   ! memory buffer information for XT
    TYPE (memory_info), POINTER :: mixtm1 ! memory buffer information for XTM1
  END TYPE t_trlist

  !
  ! submodule list entry data type
  !
  TYPE t_submlist
    CHARACTER(len=ln) :: modulename ! name of sub-model
    INTEGER           :: id         ! id   of sub-model
  END TYPE t_submlist

!------------------------------------------------------------------------------
  !
  ! module variables
  !

  TYPE(t_trlist)   ,SAVE ,TARGET :: trlist        ! tracer list

  TYPE(t_submlist) ,SAVE         :: submlist (ns) ! submodel list
  INTEGER          ,SAVE         :: nsubm = 0     ! number of entries in list

END MODULE mo_tracdef
