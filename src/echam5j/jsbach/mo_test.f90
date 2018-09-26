MODULE mo_test

  USE mo_kind,        ONLY: dp
  USE mo_jsbach_grid, ONLY: grid_type, domain_type
  USE mo_jsbach, ONLY: debug, missing_value
  USE mo_linked_list, ONLY: t_stream, GAUSSIAN, NETCDF, GRIB
  USE mo_exception, ONLY: message, int2string, real2string
  USE mo_memory_g3b

#if defined (__SX__) && defined (_OPENMP)
  USE omp_lib,          ONLY: omp_get_thread_num, &
                              omp_get_num_threads
#endif

  IMPLICIT NONE

  REAL(dp), POINTER, DIMENSION(:) :: test_hack_var1, test_hack_var2, test_hack_var3

  REAL(dp), POINTER, DIMENSION(:,:) :: &
       test0, test1, test2, test3, test4, test5, test6, test7, test8, test9, &
       test10, test11, test12, test13, test14, test15, test16, test17, test18, test19

  REAL(dp), POINTER, DIMENSION(:,:) :: &
       etest0, etest1, etest2, etest3, etest4, etest5, etest6, etest7, etest8, etest9, &
       etest10, etest11, etest12, etest13, etest14, etest15, etest16, etest17, etest18, etest19

  REAL(dp), POINTER, DIMENSION(:,:) :: &
       diff0, diff1, diff2, diff3, diff4, diff5, diff6, diff7, diff8, diff9, &
       diff10, diff11, diff12, diff13, diff14, diff15, diff16, diff17, diff18, diff19

  TYPE(t_stream), SAVE, POINTER :: IO_test
  LOGICAL, SAVE, POINTER :: domain_mask(:,:)
  INTEGER, SAVE :: domain_ndim, domain_nblocks, domain_nland
  REAL(dp), POINTER :: cover_fract(:,:)
  LOGICAL, POINTER :: cover_mask(:,:)
  INTEGER :: ntiles
  INTEGER, SAVE :: ilon0, ilon1, ilat0, ilat1
  INTEGER :: kidx0, kidx1, kblock
!$OMP THREADPRIVATE(kidx0, kidx1, kblock)

  INTERFACE test
     MODULE PROCEDURE test_0d_name
     MODULE PROCEDURE test_1d
     MODULE PROCEDURE test_1d_name
     MODULE PROCEDURE test_2d
     MODULE PROCEDURE test_2d_name
  END INTERFACE

  INTERFACE put_test_jvar
     MODULE PROCEDURE put_test_jvar_0d
  END INTERFACE

  INTERFACE put_test_var
     MODULE PROCEDURE put_test_var_1d
     MODULE PROCEDURE put_test_var_0d
  END INTERFACE

  INTERFACE put_test_var_nomask
     MODULE PROCEDURE put_test_var_nomask_1d
     MODULE PROCEDURE put_test_var_nomask_0d
  END INTERFACE

CONTAINS

  SUBROUTINE test_init_memory(grid, domain, stream)

    USE mo_netCDF,      ONLY : max_dim_name
    USE mo_memory_base, ONLY: new_stream, default_stream_setting, add =>add_stream_element
    USE mo_time_event,   ONLY: io_time_event, TIME_INC_MINUTES, TRIG_FIRST
    USE mo_time_control, ONLY: delta_time

    TYPE(grid_type), INTENT(in) :: grid
    TYPE(domain_type), INTENT(in) :: domain
    TYPE(t_stream), POINTER, OPTIONAL :: stream

    INTEGER  :: dim2p(2), dim2(2)
    CHARACTER(LEN=max_dim_name) :: dim2n(2)

    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'test', filetype=NETCDF, lpost=.TRUE., lrerun=.FALSE., &
               interval = io_time_event(INT(delta_time/60),TIME_INC_MINUTES, TRIG_FIRST, 0))
          ! Set default stream options
          CALL default_stream_setting(stream, table  = 198, repr=GAUSSIAN, lpost=.TRUE., lrerun=.FALSE.)
       ENDIF
       IO_test => stream
    ELSE
       ! Add new stream
       CALL new_stream(IO_test, 'test', filetype=NETCDF, lpost=.TRUE., lrerun=.FALSE., &
            interval = io_time_event(INT(delta_time/60),TIME_INC_MINUTES, TRIG_FIRST, 0))
       ! Set default stream options
       CALL default_stream_setting(IO_test, table  = 198, repr=GAUSSIAN, lpost=.TRUE., lrerun=.FALSE.)
    ENDIF

    dim2p = (/ domain%ndim, domain%nblocks /)
    dim2  = (/ grid%nlon, grid%nlat /)
    dim2n(1) = 'lon'
    dim2n(2) = 'lat'

    CALL add(IO_test, 'test0', test0, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=1)
    CALL add(IO_test, 'test1', test1, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=2)
    CALL add(IO_test, 'test2', test2, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=3)
    CALL add(IO_test, 'test3', test3, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=4)
    CALL add(IO_test, 'test4', test4, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=5)
    CALL add(IO_test, 'test5', test5, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=6)
    CALL add(IO_test, 'test6', test6, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=7)
    CALL add(IO_test, 'test7', test7, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=8)
    CALL add(IO_test, 'test8', test8, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=9)
    CALL add(IO_test, 'test9', test9, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=10)
    CALL add(IO_test, 'test10', test10, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=11)
    CALL add(IO_test, 'test11', test11, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=12)
    CALL add(IO_test, 'test12', test12, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=13)
    CALL add(IO_test, 'test13', test13, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=14)
    CALL add(IO_test, 'test14', test14, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=15)
    CALL add(IO_test, 'test15', test15, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=16)
    CALL add(IO_test, 'test16', test16, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=17)
    CALL add(IO_test, 'test17', test17, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=18)
    CALL add(IO_test, 'test18', test18, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=19)
    CALL add(IO_test, 'test19', test19, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=20)

    CALL add(IO_test, 'etest0', etest0, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=51)
    CALL add(IO_test, 'etest1', etest1, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=52)
    CALL add(IO_test, 'etest2', etest2, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=53)
    CALL add(IO_test, 'etest3', etest3, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=54)
    CALL add(IO_test, 'etest4', etest4, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=55)
    CALL add(IO_test, 'etest5', etest5, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=56)
    CALL add(IO_test, 'etest6', etest6, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=57)
    CALL add(IO_test, 'etest7', etest7, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=58)
    CALL add(IO_test, 'etest8', etest8, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=59)
    CALL add(IO_test, 'etest9', etest9, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=60)
    CALL add(IO_test, 'etest10', etest10, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=61)
    CALL add(IO_test, 'etest11', etest11, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=62)
    CALL add(IO_test, 'etest12', etest12, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=63)
    CALL add(IO_test, 'etest13', etest13, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=64)
    CALL add(IO_test, 'etest14', etest14, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=65)
    CALL add(IO_test, 'etest15', etest15, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=66)
    CALL add(IO_test, 'etest16', etest16, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=67)
    CALL add(IO_test, 'etest17', etest17, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=68)
    CALL add(IO_test, 'etest18', etest18, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=69)
    CALL add(IO_test, 'etest19', etest19, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=70)

    CALL add(IO_test, 'diff0', diff0, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=101)
    CALL add(IO_test, 'diff1', diff1, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=102)
    CALL add(IO_test, 'diff2', diff2, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=103)
    CALL add(IO_test, 'diff3', diff3, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=104)
    CALL add(IO_test, 'diff4', diff4, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=105)
    CALL add(IO_test, 'diff5', diff5, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=106)
    CALL add(IO_test, 'diff6', diff6, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=107)
    CALL add(IO_test, 'diff7', diff7, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=108)
    CALL add(IO_test, 'diff8', diff8, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=109)
    CALL add(IO_test, 'diff9', diff9, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=110)
    CALL add(IO_test, 'diff10', diff10, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=111)
    CALL add(IO_test, 'diff11', diff11, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=112)
    CALL add(IO_test, 'diff12', diff12, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=113)
    CALL add(IO_test, 'diff13', diff13, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=114)
    CALL add(IO_test, 'diff14', diff14, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=115)
    CALL add(IO_test, 'diff15', diff15, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=116)
    CALL add(IO_test, 'diff16', diff16, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=117)
    CALL add(IO_test, 'diff17', diff17, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=118)
    CALL add(IO_test, 'diff18', diff18, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=119)
    CALL add(IO_test, 'diff19', diff19, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=120)

    test0   = missing_value
    test1   = missing_value
    test2   = missing_value
    test3   = missing_value
    test4   = missing_value
    test5   = missing_value
    test6   = missing_value
    test7   = missing_value
    test8   = missing_value
    test9   = missing_value
    test10  = missing_value
    test11  = missing_value
    test12  = missing_value
    test13  = missing_value
    test14  = missing_value
    test15  = missing_value
    test16  = missing_value
    test17  = missing_value
    test18  = missing_value
    test19  = missing_value
    etest0  = missing_value
    etest1  = missing_value
    etest2  = missing_value
    etest3  = missing_value
    etest4  = missing_value
    etest5  = missing_value
    etest6  = missing_value
    etest7  = missing_value
    etest8  = missing_value
    etest9  = missing_value
    etest10 = missing_value
    etest11 = missing_value
    etest12 = missing_value
    etest13 = missing_value
    etest14 = missing_value
    etest15 = missing_value
    etest16 = missing_value
    etest17 = missing_value
    etest18 = missing_value
    etest19 = missing_value
    diff0   = missing_value
    diff1   = missing_value
    diff2   = missing_value
    diff3   = missing_value
    diff4   = missing_value
    diff5   = missing_value
    diff6   = missing_value
    diff7   = missing_value
    diff8   = missing_value
    diff9   = missing_value
    diff10  = missing_value
    diff11  = missing_value
    diff12  = missing_value
    diff13  = missing_value
    diff14  = missing_value
    diff15  = missing_value
    diff16  = missing_value
    diff17  = missing_value
    diff18  = missing_value
    diff19  = missing_value

    ALLOCATE(test_hack_var1(domain%ndim))
    ALLOCATE(test_hack_var2(domain%ndim))
    ALLOCATE(test_hack_var3(domain%ndim))

  END SUBROUTINE test_init_memory

  SUBROUTINE init_test(grid, domain, itiles, IO_stream)

    USE mo_jsbach_grid, ONLY: grid_type, domain_type

    TYPE(grid_type), INTENT(in)          :: grid
    TYPE(domain_type), INTENT(in)        :: domain
    INTEGER, INTENT(in)                  :: itiles
    TYPE(t_stream), POINTER, OPTIONAL    :: IO_stream

    INTEGER :: i

    CALL test_init_memory(grid, domain, IO_stream)

    ALLOCATE(domain_mask(domain%ndim,domain%nblocks))
    domain_mask = domain%mask
    domain_ndim = domain%ndim
    domain_nblocks = domain%nblocks
    domain_nland = domain%nland
    ntiles = itiles
    ALLOCATE(cover_fract(grid%nland,ntiles))
    ALLOCATE(cover_mask(grid%nland,ntiles))

    DO i=1,grid%nlon
       IF (ABS(domain%lon(1)-grid%lon(i))<EPSILON(1._dp)) ilon0 = i
       IF (ABS(domain%lon(domain%nland)-grid%lon(i))<EPSILON(1._dp)) ilon1 = i
    END DO
    DO i=1,grid%nlat
       IF (ABS(domain%lat(1)-grid%lat(i))<EPSILON(1._dp)) ilat0 = i
       IF (ABS(domain%lat(domain%nland)-grid%lat(i))<EPSILON(1._dp)) ilat1 = i
    END DO

  END SUBROUTINE init_test

  SUBROUTINE init_test_2(i1,i2, iblock, fract, mask, domain)

    INTEGER, INTENT(in) :: i1, i2, iblock
    REAL(dp), INTENT(in)                     :: fract(:,:)
    LOGICAL, INTENT(in)                  :: mask(:,:)
    TYPE(domain_type), INTENT(in) :: domain

    kidx0 = i1
    kidx1 = i2
    kblock = iblock

    cover_fract(kidx0:kidx1,:) = fract
    cover_mask(kidx0:kidx1,:) = mask

    domain_mask(:,iblock) = domain%mask(:,iblock)
    domain_ndim = domain%ndim
    domain_nblocks = domain%nblocks
    domain_nland = domain%nland

  END SUBROUTINE init_test_2

  SUBROUTINE test_2d_name(jvar, index, ecode, text, mask)

    USE mo_utils, ONLY: average_tiles
    USE mo_memory_base, ONLY: get_stream_element, set_stream_element_info
    USE mo_exception, ONLY: int2string

    REAL(dp), INTENT(in) :: jvar(:,:)
    INTEGER, INTENT(in) :: index
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: text
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ecode
    LOGICAL, INTENT(in), OPTIONAL :: mask(:,:)
    
    REAL(dp) :: jvar_avg(SIZE(jvar,1))
    REAL(dp), POINTER :: evar(:,:), evar2(:,:)
    REAL(dp), POINTER :: testvar(:,:), diffvar(:,:)
    CHARACTER(LEN=8) :: cname

    IF (PRESENT(mask)) THEN
       CALL average_tiles(jvar, cover_mask(kidx0:kidx1,:).AND.mask, cover_fract(kidx0:kidx1,:), jvar_avg)
    ELSE
       CALL average_tiles(jvar, cover_mask(kidx0:kidx1,:), cover_fract(kidx0:kidx1,:), jvar_avg)
    END IF
    cname = 'test'//int2string(index)
    IF (PRESENT(text)) THEN
       CALL set_stream_element_info(IO_test, cname, longname=text)
    END IF
    CALL get_stream_element(IO_test, cname, testvar)
    testvar(:,kblock) = UNPACK(jvar_avg, domain_mask(:,kblock), missing_value)

    IF (PRESENT(ecode)) THEN
       IF (LEN_TRIM(ecode)>5 .AND. ecode(1:MIN(5,LEN(ecode))) == 'etest') THEN
          CALL get_stream_element(IO_test, ecode, evar)
       ELSE
          CALL get_stream_element(g3b, ecode, evar)
          cname = 'etest'//int2string(index)
          CALL get_stream_element(IO_test, cname, evar2)
          evar2(:,kblock) = MERGE(evar(:,kblock), missing_value, domain_mask(:,kblock))
       END IF
       cname = 'diff'//int2string(index)
       CALL get_stream_element(IO_test, cname, diffvar)
       diffvar(:,kblock) = MERGE(testvar(:,kblock) - evar(:,kblock), missing_value, domain_mask(:,kblock))
    END IF

    IF (debug .AND. kblock==domain_nblocks) THEN
       IF (PRESENT(ecode)) THEN
          PRINT*, text, ' (difference) ', MINVAL(diffvar), &
               MAXVAL(MERGE(diffvar,-missing_value,diffvar<0.95_dp*missing_value))
       ELSE
          PRINT*, text, MINVAL(testvar), MAXVAL(MERGE(testvar,-missing_value,testvar<0.95_dp*missing_value))
       END IF
    END IF

  END SUBROUTINE test_2d_name

  SUBROUTINE test_2d(jvar, evar, index, text, mask)

    USE mo_utils, ONLY: average_tiles
    USE mo_memory_base, ONLY: get_stream_element, set_stream_element_info
    USE mo_exception, ONLY: int2string

    REAL(dp), INTENT(in) :: jvar(:,:)
    REAL(dp), INTENT(in) :: evar(:)
    INTEGER, INTENT(in) :: index
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: text
    LOGICAL, INTENT(in), OPTIONAL :: mask(:,:)
    
    REAL(dp) :: jvar_avg(SIZE(jvar,1))
    REAL(dp), POINTER :: testvar(:,:), evar2(:,:)
    CHARACTER(LEN=8) :: cname

    IF (PRESENT(mask)) THEN
       CALL average_tiles(jvar, cover_mask(kidx0:kidx1,:).AND.mask, cover_fract(kidx0:kidx1,:), jvar_avg)
    ELSE
       CALL average_tiles(jvar, cover_mask(kidx0:kidx1,:), cover_fract(kidx0:kidx1,:), jvar_avg)
    END IF
    cname = 'test'//int2string(index)
    CALL get_stream_element(IO_test, cname, testvar)
    IF (PRESENT(text)) THEN
       CALL set_stream_element_info(IO_test, cname, longname=text)
    END IF
    testvar(:,kblock) = UNPACK(jvar_avg, domain_mask(:,kblock), missing_value)
    cname = 'etest'//int2string(index)
    CALL get_stream_element(IO_test, cname, evar2)
    evar2(:,kblock) = MERGE(evar, missing_value, domain_mask(:,kblock))

    testvar(:,kblock) = MERGE(testvar(:,kblock) - evar(:), missing_value, domain_mask(:,kblock))

    IF (debug .AND. kblock==domain_nblocks) &
         PRINT*, text, MINVAL(testvar), MAXVAL(MERGE(testvar,-missing_value,testvar<0.95_dp*missing_value))

  END SUBROUTINE test_2d

  SUBROUTINE test_1d_name(jvar, index, ecode, text)

    USE mo_utils, ONLY: average_tiles
    USE mo_memory_base, ONLY: get_stream_element, set_stream_element_info
    USE mo_exception, ONLY: int2string

    REAL(dp), INTENT(in) :: jvar(:)
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ecode
    INTEGER, INTENT(in) :: index
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: text
    
    REAL(dp) :: jvar_avg(SIZE(jvar))
    REAL(dp), POINTER :: evar(:,:), evar2(:,:)
    REAL(dp), POINTER :: testvar(:,:), diffvar(:,:)
    CHARACTER(LEN=8) :: cname

#if defined (__SX__) && defined (_OPENMP)
    INTEGER :: tid
#endif

#if defined (__SX__) && defined (_OPENMP)
!!$    IF (debug) THEN
!!$       tid = omp_get_thread_num()
!!$       CALL message('test_1d_name', 'OpenMP thread #'//int2string(tid)//' '//int2string(kblock))
!!$    END IF
#endif

    cname = 'test'//int2string(index)
    CALL get_stream_element(IO_test, cname, testvar)
    IF (PRESENT(text)) THEN
       CALL set_stream_element_info(IO_test, cname, longname=text)
    END IF
    testvar(:,kblock) = UNPACK(jvar, domain_mask(:,kblock), missing_value)


    IF (PRESENT(ecode)) THEN
       IF (LEN_TRIM(ecode)>5 .AND. ecode(1:MIN(5,LEN(ecode))) == 'etest') THEN
          CALL get_stream_element(IO_test, ecode, evar)
       ELSE
          CALL get_stream_element(g3b, ecode, evar)
          cname = 'etest'//int2string(index)
          CALL get_stream_element(IO_test, cname, evar2)
          evar2(:,kblock) = MERGE(evar(:,kblock), missing_value, domain_mask(:,kblock))
       END IF
       cname = 'diff'//int2string(index)
       CALL get_stream_element(IO_test, cname, diffvar)
       diffvar(:,kblock) = MERGE(testvar(:,kblock) - evar(:,kblock), missing_value, domain_mask(:,kblock))
    END IF

    IF (debug .AND. kblock==domain_nblocks) THEN
       IF (PRESENT(ecode)) THEN
          PRINT*, text, ' (difference) ', MINVAL(diffvar), &
               MAXVAL(MERGE(diffvar,-missing_value,testvar<0.95_dp*missing_value))
       ELSE
          PRINT*, text, MINVAL(testvar), MAXVAL(MERGE(testvar,-missing_value,testvar<0.95_dp*missing_value))
       END IF
    END IF

  END SUBROUTINE test_1d_name

  SUBROUTINE test_0d_name(index, ecode, text)

    USE mo_utils, ONLY: average_tiles
    USE mo_memory_base, ONLY: get_stream_element, set_stream_element_info
    USE mo_exception, ONLY: int2string

    INTEGER, INTENT(in) :: index
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ecode
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: text
    
    REAL(dp), POINTER :: evar(:,:), evar2(:,:)
    REAL(dp), POINTER :: testvar(:,:), diffvar(:,:)
    CHARACTER(LEN=8) :: cname

!    cname = 'test'//int2string(index)
    CALL get_stream_element(IO_test, cname, testvar)
    IF (PRESENT(text)) THEN
       CALL set_stream_element_info(IO_test, cname, longname=text)
    END IF

    IF (PRESENT(ecode)) THEN
       IF (LEN_TRIM(ecode)>5 .AND. ecode(1:MIN(5,LEN(ecode))) == 'etest') THEN
          CALL get_stream_element(IO_test, ecode, evar)
       ELSE
          CALL get_stream_element(g3b, ecode, evar)
          cname = 'etest'//int2string(index)
          CALL get_stream_element(IO_test, cname, evar2)
          evar2(:,kblock) = MERGE(evar(:,kblock), missing_value, domain_mask(:,kblock))
       END IF
       cname = 'diff'//int2string(index)
       CALL get_stream_element(IO_test, cname, diffvar)
       diffvar(:,kblock) = MERGE(testvar(:,kblock) - evar(:,kblock), missing_value, domain_mask(:,kblock))
    END IF

    IF (debug .AND. kblock==domain_nblocks) THEN
       IF (PRESENT(ecode)) THEN
          PRINT*, text, ' (difference) ', MINVAL(diffvar), &
               MAXVAL(MERGE(diffvar,-missing_value,testvar<0.95_dp*missing_value))
       ELSE
          PRINT*, text, MINVAL(testvar), MAXVAL(MERGE(testvar,-missing_value,testvar<0.95_dp*missing_value))
       END IF
    END IF

  END SUBROUTINE test_0d_name

  SUBROUTINE put_test_jvar_0d(jvar, index, iland)

    USE mo_memory_base, ONLY: get_stream_element
    USE mo_exception, ONLY: int2string
 
    REAL(dp), INTENT(in) :: jvar
    INTEGER, INTENT(in) :: index, iland

    REAL(dp) :: var(domain_nland)
    REAL(dp), POINTER :: testvar(:,:)
    CHARACTER(LEN=8) :: cname

    cname = 'test'//int2string(index)
    CALL get_stream_element(IO_test, cname, testvar)
    var = PACK(testvar(:,kblock), domain_mask(:,kblock))
    var(iland) = jvar
    testvar(:,kblock) = UNPACK(var, domain_mask(:,kblock), missing_value)

  END SUBROUTINE put_test_jvar_0d

  SUBROUTINE test_1d(jvar, evar, index, text)

    USE mo_utils, ONLY: average_tiles
    USE mo_memory_base, ONLY: get_stream_element, set_stream_element_info
    USE mo_exception, ONLY: int2string

    REAL(dp), INTENT(in) :: jvar(:)
    REAL(dp), INTENT(in) :: evar(:)
    INTEGER, INTENT(in) :: index
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: text
    
    REAL(dp) :: jvar_avg(SIZE(jvar,1))
    REAL(dp), POINTER :: testvar(:,:), evar2(:,:)
    CHARACTER(LEN=8) :: cname

    cname = 'test'//int2string(index)
    CALL get_stream_element(IO_test, cname, testvar)
    IF (PRESENT(text)) THEN
       CALL set_stream_element_info(IO_test, cname, longname=text)
    END IF
    testvar(:,kblock) = UNPACK(jvar, domain_mask(:,kblock), missing_value)
    cname = 'etest'//int2string(index)
    CALL get_stream_element(IO_test, cname, evar2)
    evar2(:,kblock) = MERGE(evar, missing_value, domain_mask(:,kblock))

    testvar(:,kblock) = MERGE(testvar(:,kblock) - evar(:), missing_value, domain_mask(:,kblock))

    IF (debug .AND. kblock==domain_nblocks) THEN
       IF (PRESENT(text)) THEN
          PRINT*, text, MINVAL(testvar), MAXVAL(MERGE(testvar,-missing_value,testvar<0.95_dp*missing_value))
       ELSE
          PRINT*, MINVAL(testvar), MAXVAL(MERGE(testvar,-missing_value,testvar<0.95_dp*missing_value))
       END IF
    END IF

   
  END SUBROUTINE test_1d


  SUBROUTINE put_test_var_1d(jrow, evar, index)

    USE mo_exception, ONLY: int2string
    USE mo_memory_base, ONLY: get_stream_element, set_stream_element_info

    REAL(dp), INTENT(in) :: evar(:)
    INTEGER, INTENT(in) :: index, jrow
!    INTEGER, INTENT(in) :: iblock

    REAL(dp), POINTER :: testvar(:,:)

    CHARACTER(LEN=9) :: cname

    cname = 'etest'//int2string(index)

    CALL get_stream_element(IO_test, cname, testvar)

    testvar(:,jrow) = MERGE(evar, missing_value, domain_mask(:,jrow))
    

  END SUBROUTINE put_test_var_1d

  SUBROUTINE put_test_var_0d(jrow, evar, index, iproma)

    USE mo_exception, ONLY: int2string
    USE mo_memory_base, ONLY: get_stream_element, set_stream_element_info

    REAL(dp), INTENT(in) :: evar
    INTEGER, INTENT(in) :: index
    INTEGER, INTENT(in) :: iproma, jrow

    REAL(dp), POINTER :: testvar(:,:)

    CHARACTER(LEN=9) :: cname

    cname = 'etest'//int2string(index)

    CALL get_stream_element(IO_test, cname, testvar)

    IF (domain_mask(iproma,jrow)) THEN
       testvar(iproma,jrow) = evar
    ELSE
       testvar(iproma,jrow) = missing_value
    END IF
    
  END SUBROUTINE put_test_var_0d

  SUBROUTINE put_test_var_nomask_1d(jrow, evar, index)

    USE mo_exception, ONLY: int2string
    USE mo_memory_base, ONLY: get_stream_element, set_stream_element_info

    REAL(dp), INTENT(in) :: evar(:)
    INTEGER, INTENT(in) :: index, jrow
!    INTEGER, INTENT(in) :: iblock

    REAL(dp), POINTER :: testvar(:,:)
    INTEGER :: nproma

    CHARACTER(LEN=9) :: cname

    nproma = SIZE(evar)

    cname = 'etest'//int2string(index)

!!$    CALL message('put_test_var_nomask_1d',TRIM(cname)//' - '//TRIM(int2string(nproma))//', '//&
!!$         TRIM(int2string(jrow))//' - '//TRIM(real2string(MINVAL(evar(1:nproma))))//', '//&
!!$         TRIM(real2string(MAXVAL(evar(1:nproma)))))

    CALL get_stream_element(IO_test, cname, testvar)

    testvar(1:nproma,jrow) = evar(:)

  END SUBROUTINE put_test_var_nomask_1d

  SUBROUTINE put_test_var_nomask_0d(jrow, evar, index, iproma)

    USE mo_exception, ONLY: int2string
    USE mo_memory_base, ONLY: get_stream_element, set_stream_element_info

    INTEGER, INTENT(in) :: jrow
    REAL(dp), INTENT(in) :: evar
    INTEGER, INTENT(in) :: index
    INTEGER, INTENT(in) :: iproma

    REAL(dp), POINTER :: testvar(:,:)

    CHARACTER(LEN=9) :: cname

    cname = 'etest'//int2string(index)

    CALL get_stream_element(IO_test, cname, testvar)

    testvar(iproma,jrow) = evar
    
  END SUBROUTINE put_test_var_nomask_0d

END MODULE mo_test
