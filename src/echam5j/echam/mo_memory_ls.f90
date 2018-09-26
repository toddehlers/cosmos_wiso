MODULE mo_memory_ls

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: delete_stream, add_stream_element, &
                            default_stream_setting, SPECTRAL

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_ls ! construct the ls table
  PUBLIC :: destruct_ls  ! destruct  the ls table

  PUBLIC :: ls           ! the ls table

  PUBLIC :: add_stream_element

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC  :: lvo (:,:,:)
  REAL(dp), POINTER, PUBLIC  :: ld  (:,:,:) 
  REAL(dp), POINTER, PUBLIC  :: ltp (:,:,:)
  REAL(dp), POINTER, PUBLIC  :: lu0 (:,:) 

  ! declaration of table with 3d-field entries

  TYPE (t_stream), POINTER :: ls

CONTAINS

  SUBROUTINE construct_ls (lnlev, lnlevp1, lnnp1, lnsp, nlev, nnp1, nsp)

    INTEGER, INTENT (in) :: lnlev, lnlevp1, lnnp1, lnsp ! local array bounds
    INTEGER, INTENT (in) :: nlev, nnp1, nsp             ! global array bounds

    ! construct the ls table
    !
    ! all information specific to this table is set in this subroutine

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    ! assign pointers

    CALL default_stream_setting (ls, repr   = SPECTRAL, &
                                     lpost  = .FALSE., &
                                     lrerun = .FALSE.)

    CALL add_stream_element (ls,'svo',lvo, (/ lnlev,   2, lnsp /),&
                                           (/ nlev,    2, nsp  /))
    CALL add_stream_element (ls,'sd', ld,  (/ lnlev,   2, lnsp /),&
                                           (/ nlev,    2, nsp  /))
    CALL add_stream_element (ls,'stp',ltp, (/ lnlevp1, 2, lnsp /),&
                                           (/ nlev+1,  2, nsp  /))
    CALL add_stream_element (ls,'su0',lu0, (/ lnlev,      lnnp1/),&
                                           (/ nlev,       nnp1 /))

  END SUBROUTINE construct_ls

  SUBROUTINE destruct_ls

    CALL delete_stream (ls)

  END SUBROUTINE destruct_ls

END MODULE mo_memory_ls

