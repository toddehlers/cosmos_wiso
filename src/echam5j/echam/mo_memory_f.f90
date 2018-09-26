#if defined(__uxp__) || defined(__SX__) || defined (ES)
#define FAST_AND_DIRTY 1
#endif

MODULE mo_memory_f
  !
  ! declaration of predefined fields within this module 
  !
  ! used in inverse Legendre transformation
  !

  ! Modules used

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: add_stream_element, delete_stream, &
                            default_stream_setting, FOURIER
  USE mo_netCDF,      ONLY: max_dim_name
  IMPLICIT NONE

  ! public entities

  PRIVATE
  PUBLIC :: construct_f ! construct the f table
  PUBLIC :: destruct_f  ! destruct  the f table
  PUBLIC :: f           ! the f table
#ifdef FAST_AND_DIRTY
  PUBLIC :: resort_restart_write_memory_f
  PUBLIC :: resort_restart_read_memory_f
#endif

  ! pointers for *f1* space.

  REAL(dp), POINTER, PUBLIC  :: fsvo(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: favo(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fsu(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fau(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fsv(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fav(:,:,:,:)

  ! pointers for *f3* space.

  REAL(dp), POINTER, PUBLIC  :: fsd(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fad(:,:,:,:)

  ! pointers for *f4* space.

  REAL(dp), POINTER, PUBLIC  :: fstp(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fatp(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fstpm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fatpm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fsu0(:,:)
  REAL(dp), POINTER, PUBLIC  :: fau0(:,:)
  REAL(dp), POINTER, PUBLIC  :: fsdu0(:,:)
  REAL(dp), POINTER, PUBLIC  :: fadu0(:,:)

  ! used in direct Legendre transformation
  !
  ! pointers for *f5* space.

  REAL(dp), POINTER, PUBLIC :: fszl(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fazl(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fszm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fazm(:,:,:,:)

  ! pointers for *f7* space.

  REAL(dp), POINTER, PUBLIC :: fsdl(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fadl(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fsdm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fadm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fsr(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: far(:,:,:,:)

  ! pointers for *f8* space.

  REAL(dp), POINTER, PUBLIC :: fstp1(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fatp1(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fsul(:,:)
  REAL(dp), POINTER, PUBLIC :: faul(:,:)

#ifdef FAST_AND_DIRTY
  ! pointers for *f1* space IO.

  REAL(dp), POINTER, PUBLIC  :: ofsvo(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: ofavo(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: ofsu(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: ofau(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: ofsv(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: ofav(:,:,:,:)

  ! pointers for *f3* space IO.

  REAL(dp), POINTER, PUBLIC  :: ofsd(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: ofad(:,:,:,:)

  ! pointers for *f4* space IO.

  REAL(dp), POINTER, PUBLIC  :: ofstp(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: ofatp(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: ofstpm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: ofatpm(:,:,:,:)
#endif

  ! declaration of table with 3d-field entries

  TYPE (t_stream), POINTER :: f

CONTAINS

  SUBROUTINE construct_f (lnlev, lnlevp1, lnmp1, lnhgl, nlev, nmp1, nhgl)

    INTEGER, INTENT (in) :: lnlev, lnlevp1, lnmp1, lnhgl
    INTEGER, INTENT (in) ::  nlev,           nmp1,  nhgl
    INTEGER :: dim1(4), dim1p(4)
    INTEGER :: dim2(4), dim2p(4)
    INTEGER :: dim3(2), dim3p(2)
    CHARACTER (max_dim_name) :: dim1n(4), dim2n(4), dim3n(2)

#ifdef FAST_AND_DIRTY
    INTEGER :: dim1x(4), dim1xp(4), dim2x(4), dim2xp(4)
#endif


    ! construct the f table
    !
    ! all information specific to this table is set in this subroutine
    !
    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields
    !
    ! assign pointers

    dim1p = (/ lnlev,     2,       lnmp1,    lnhgl    /)
    dim1  = (/  nlev,     2,        nmp1,     nhgl    /)
    dim1n = (/ "lev    ","complex","nmp1   ","nhgl   "/)

    dim2p = (/ lnlevp1,   2,       lnmp1,    lnhgl    /)
    dim2  = (/  nlev+1,   2,        nmp1,     nhgl    /)
    dim2n = (/ "ilev   ","complex","nmp1   ","nhgl   "/)

#ifdef FAST_AND_DIRTY
    dim1xp = (/ 2*lnlev+1,     1,       lnmp1+1,    lnhgl    /)
    dim1x  = (/  2*nlev+1,     1,        nmp1+1,     nhgl    /)

    dim2xp = (/ 2*lnlevp1+1,   1,       lnmp1+1,    lnhgl    /)
    dim2x  = (/  2*(nlev+1)+1,   1,        nmp1+1,     nhgl    /)
#endif

    dim3p = (/ lnlev, lnhgl /)
    dim3  = (/  nlev,  nhgl /)
    dim3n = (/ "lev ","nhgl"/)

    ! Arrays used by inverse transform

    CALL default_stream_setting (f, repr   = FOURIER, &
                                    lpost  = .FALSE., &
                                    lrerun = .TRUE.)

#ifdef FAST_AND_DIRTY
    CALL add_stream_element (f, 'fsvo',  fsvo,  dim1xp, dim1x, dimnames=dim1n, lrerun=.FALSE.)
    CALL add_stream_element (f, 'favo',  favo,  dim1xp, dim1x, dimnames=dim1n, lrerun=.FALSE.)
    CALL add_stream_element (f, 'fsu',   fsu,   dim1xp, dim1x, dimnames=dim1n, lrerun=.FALSE.)
    CALL add_stream_element (f, 'fau',   fau,   dim1xp, dim1x, dimnames=dim1n, lrerun=.FALSE.)
    CALL add_stream_element (f, 'fsv',   fsv,   dim1xp, dim1x, dimnames=dim1n, lrerun=.FALSE.)
    CALL add_stream_element (f, 'fav',   fav,   dim1xp, dim1x, dimnames=dim1n, lrerun=.FALSE.)
    CALL add_stream_element (f, 'fsd',   fsd,   dim1xp, dim1x, dimnames=dim1n, lrerun=.FALSE.)
    CALL add_stream_element (f, 'fad',   fad,   dim1xp, dim1x, dimnames=dim1n, lrerun=.FALSE.)
    CALL add_stream_element (f, 'fstp',  fstp,  dim2xp, dim2x, dimnames=dim2n, lrerun=.FALSE.)
    CALL add_stream_element (f, 'fatp',  fatp,  dim2xp, dim2x, dimnames=dim2n, lrerun=.FALSE.)
    CALL add_stream_element (f, 'fstpm', fstpm, dim2xp, dim2x, dimnames=dim2n, lrerun=.FALSE.)
    CALL add_stream_element (f, 'fatpm', fatpm, dim2xp, dim2x, dimnames=dim2n, lrerun=.FALSE.)
    ! required for output ----------------------------------------------------
    CALL add_stream_element (f, 'fsvo',  ofsvo,  dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'favo',  ofavo,  dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fsu',   ofsu,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fau',   ofau,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fsv',   ofsv,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fav',   ofav,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fsd',   ofsd,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fad',   ofad,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fstp',  ofstp,  dim2p, dim2, dimnames=dim2n)
    CALL add_stream_element (f, 'fatp',  ofatp,  dim2p, dim2, dimnames=dim2n)
    CALL add_stream_element (f, 'fstpm', ofstpm, dim2p, dim2, dimnames=dim2n)
    CALL add_stream_element (f, 'fatpm', ofatpm, dim2p, dim2, dimnames=dim2n)
#else
    CALL add_stream_element (f, 'fsvo',  fsvo,  dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'favo',  favo,  dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fsu',   fsu,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fau',   fau,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fsv',   fsv,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fav',   fav,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fsd',   fsd,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fad',   fad,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fstp',  fstp,  dim2p, dim2, dimnames=dim2n)
    CALL add_stream_element (f, 'fatp',  fatp,  dim2p, dim2, dimnames=dim2n)
    CALL add_stream_element (f, 'fstpm', fstpm, dim2p, dim2, dimnames=dim2n)
    CALL add_stream_element (f, 'fatpm', fatpm, dim2p, dim2, dimnames=dim2n)
#endif
    CALL add_stream_element (f, 'fsu0',  fsu0,  dim3p, dim3, dimnames=dim3n)
    CALL add_stream_element (f, 'fau0',  fau0,  dim3p, dim3, dimnames=dim3n)
    CALL add_stream_element (f, 'fsdu0', fsdu0, dim3p, dim3, dimnames=dim3n)
    CALL add_stream_element (f, 'fadu0', fadu0, dim3p, dim3, dimnames=dim3n)

    ! Arrays used by direct transform (not yet in memory buffer)

    CALL default_stream_setting (f, lrerun= .FALSE.)

    CALL add_stream_element (f, 'fsdl',  fsdl, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fadm',  fadm, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fsr',   fsr,  dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fszl',  fszl, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fazm',  fazm, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fstp1', fstp1,dim2p, dim2, dimnames=dim2n)
    CALL add_stream_element (f, 'fsul',  fsul, dim3p, dim3, dimnames=dim3n)


    CALL add_stream_element (f, 'fadl', fadl, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fsdm', fsdm, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'far',  far,  dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fazl', fazl, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fszm', fszm, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fatp1',fatp1,dim2p, dim2, dimnames=dim2n)
    CALL add_stream_element (f, 'faul', faul, dim3p, dim3, dimnames=dim3n)

  END SUBROUTINE construct_f

  SUBROUTINE destruct_f

    CALL delete_stream (f)

  END SUBROUTINE destruct_f

#ifdef FAST_AND_DIRTY

  SUBROUTINE resort_restart_write_memory_f
 
    INTEGER :: lnlev, lnlevp1, lnmp1, lnhgl 
    INTEGER :: i, j, k

    lnlev   = SIZE(ofsvo,1)
    lnlevp1 = SIZE(ofstp,1)
    lnmp1   = SIZE(ofsvo,3) 
    lnhgl   = SIZE(ofsvo,4)
   
    DO k = 1, lnhgl
      DO j = 1, lnmp1
        DO i = 1, lnlev
          ofsvo(i, 1, j, k) = fsvo(i,       1, j, k) 
          ofsvo(i, 2, j, k) = fsvo(lnlev+i, 1, j, k) 
          ofavo(i, 1, j, k) = favo(i,       1, j, k) 
          ofavo(i, 2, j, k) = favo(lnlev+i, 1, j, k) 
          ofsu (i, 1, j, k) = fsu (i,       1, j, k) 
          ofsu (i, 2, j, k) = fsu (lnlev+i, 1, j, k) 
          ofau (i, 1, j, k) = fau (i,       1, j, k) 
          ofau (i, 2, j, k) = fau (lnlev+i, 1, j, k) 
          ofsv (i, 1, j, k) = fsv (i,       1, j, k) 
          ofsv (i, 2, j, k) = fsv (lnlev+i, 1, j, k) 
          ofav (i, 1, j, k) = fav (i,       1, j, k) 
          ofav (i, 2, j, k) = fav (lnlev+i, 1, j, k) 
          ofsd (i, 1, j, k) = fsd (i,       1, j, k) 
          ofsd (i, 2, j, k) = fsd (lnlev+i, 1, j, k) 
          ofad (i, 1, j, k) = fad (i,       1, j, k) 
          ofad (i, 2, j, k) = fad (lnlev+i, 1, j, k) 
        ENDDO
      ENDDO
    ENDDO
                  
    DO k = 1, lnhgl
      DO j = 1, lnmp1
        DO i = 1, lnlevp1
          ofstp (i, 1, j, k) = fstp (i,         1, j, k) 
          ofstp (i, 2, j, k) = fstp (lnlevp1+i, 1, j, k) 
          ofatp (i, 1, j, k) = fatp (i,         1, j, k)
          ofatp (i, 2, j, k) = fatp (lnlevp1+i, 1, j, k)
          ofstpm(i, 1, j, k) = fstpm(i,         1, j, k)
          ofstpm(i, 2, j, k) = fstpm(lnlevp1+i, 1, j, k)
          ofatpm(i, 1, j, k) = fatpm(i,         1, j, k)
          ofatpm(i, 2, j, k) = fatpm(lnlevp1+i, 1, j, k)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE resort_restart_write_memory_f

  SUBROUTINE resort_restart_read_memory_f
 
    INTEGER :: lnlev, lnlevp1, lnmp1, lnhgl 
    INTEGER :: i, j, k

    lnlev   = SIZE(ofsvo,1)
    lnlevp1 = SIZE(ofstp,1)
    lnmp1   = SIZE(ofsvo,3) 
    lnhgl   = SIZE(ofsvo,4)
   
    DO k = 1, lnhgl
      DO j = 1, lnmp1
        DO i = 1, lnlev
          fsvo(i,       1, j, k) = ofsvo(i, 1, j, k) 
          fsvo(lnlev+i, 1, j, k) = ofsvo(i, 2, j, k) 
          favo(i,       1, j, k) = ofavo(i, 1, j, k) 
          favo(lnlev+i, 1, j, k) = ofavo(i, 2, j, k) 
          fsu (i,       1, j, k) = ofsu (i, 1, j, k) 
          fsu (lnlev+i, 1, j, k) = ofsu (i, 2, j, k) 
          fau (i,       1, j, k) = ofau (i, 1, j, k) 
          fau (lnlev+i, 1, j, k) = ofau (i, 2, j, k) 
          fsv (i,       1, j, k) = ofsv (i, 1, j, k) 
          fsv (lnlev+i, 1, j, k) = ofsv (i, 2, j, k) 
          fav (i,       1, j, k) = ofav (i, 1, j, k) 
          fav (lnlev+i, 1, j, k) = ofav (i, 2, j, k) 
          fsd (i,       1, j, k) = ofsd (i, 1, j, k) 
          fsd (lnlev+i, 1, j, k) = ofsd (i, 2, j, k) 
          fad (i,       1, j, k) = ofad (i, 1, j, k) 
          fad (lnlev+i, 1, j, k) = ofad (i, 2, j, k) 
        ENDDO
      ENDDO
    ENDDO
                  
    DO k = 1, lnhgl
      DO j = 1, lnmp1
        DO i = 1, lnlevp1
          fstp (i,         1, j, k) = ofstp (i, 1, j, k)  
          fstp (lnlevp1+i, 1, j, k) = ofstp (i, 2, j, k) 
          fatp (i,         1, j, k) = ofatp (i, 1, j, k)
          fatp (lnlevp1+i, 1, j, k) = ofatp (i, 2, j, k)
          fstpm(i,         1, j, k) = ofstpm(i, 1, j, k)
          fstpm(lnlevp1+i, 1, j, k) = ofstpm(i, 2, j, k)
          fatpm(i,         1, j, k) = ofatpm(i, 1, j, k)
          fatpm(lnlevp1+i, 1, j, k) = ofatpm(i, 2, j, k)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE resort_restart_read_memory_f

#endif

END MODULE mo_memory_f
