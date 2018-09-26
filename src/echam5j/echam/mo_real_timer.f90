MODULE mo_real_timer

  !
  ! utility for real time measurements
  ! on AIX a wrapper for the fast read_real_time function is used
  ! on SX a CPU counter is used
  ! on Linux/Intel a CPU register is read
  !
  ! Authors:
  !
  ! Initial version by:
  !
  ! J. Behrend, GWDG, March 2002, initial version for AIX
  !
  ! History:
  !
  ! L. Kornblueh, MPI, April 2002, extended for NEC SX and Linux/Intel
  ! S. Shingu. NEC, March 2004, bugfix
  ! Luis Kornblueh, MPI fuer Meteorologie, Hamburg, May 2004, adapted 
  !          for use in ICON
  !

  USE mo_kind,        ONLY: dp
  USE mo_doctor,      ONLY: nerr
  USE mo_exception,   ONLY: finish, message, message_text
  USE mo_util_String, ONLY: separator

#ifdef _OPENMP
#if (! defined __PGI)
 USE omp_lib,       ONLY: omp_get_thread_num, omp_get_max_threads, &
                          omp_in_parallel, omp_get_num_threads
!                          omp_get_dynamic, omp_set_dynamic
#endif
#endif

#ifndef NOMPI
  USE mo_mpi,       ONLY: MPI_STATUS_SIZE,                              &
                          p_recv, p_send, p_barrier, p_abort, p_gather, &
                          p_pe, p_nprocs, p_all_comm 
#endif

  IMPLICIT NONE

  PRIVATE

#if (defined _OPENMP) && (defined __PGI)
  INTEGER, EXTERNAL :: omp_get_max_threads
  INTEGER, EXTERNAL :: omp_get_thread_num
  INTEGER, EXTERNAL :: omp_get_num_threads
  LOGICAL, EXTERNAL :: omp_in_parallel
#endif

  ! simple access - no statistics 

  PUBLIC :: time_mark_type, set_time_mark, get_time_val

  ! more informative, threadsafe statistics:

  PUBLIC :: new_timer, del_timer, timer_start, timer_stop
  PUBLIC :: timer_val, timer_last, timer_count, timer_average
  PUBLIC :: exercise_timer, timer_reset, timer_reset_all
  PUBLIC :: timer_report

#ifndef NOMPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif

  REAL(dp), EXTERNAL :: util_walltime

  INTEGER, PARAMETER :: timer_max=128

  INTEGER :: top_timer=0

  INTEGER, PARAMETER :: rt_undef_stat=0
  INTEGER, PARAMETER :: rt_on_stat=1
  INTEGER, PARAMETER :: rt_off_stat=2

  ! minimal internal time needed to do one measuremenet

  REAL(dp) :: tm_shift = 0.0_dp

  ! average total overhead for one timer_start & timer_stop pair

  REAL(dp) :: tm_overhead = 0.0_dp 

  ! shared part of timer

  TYPE srt_type
     LOGICAL           :: reserved             ! usage flag
     CHARACTER(len=80) :: text                 ! description of timer
  END TYPE srt_type

#ifdef __ibm__ 
  TYPE time_mark_type
    INTEGER :: t(4) ! 'raw' timer mark
  END TYPE time_mark_type
#else
  TYPE time_mark_type
    REAL(dp) :: t
  END TYPE time_mark_type
#endif

  ! thread private part of timer

  TYPE rt_type
     SEQUENCE
#ifdef __ibm__ 
     INTEGER  :: mark1(4)
#else
     REAL(dp) :: mark1               ! last start time
#endif
     REAL(dp) :: tot                 ! total sum of active time 
     REAL(dp) :: min                 ! min. active time
     REAL(dp) :: max                 ! max. ..
     REAL(dp) :: last                ! last measurement
     INTEGER  :: call_n              ! number of calls
     INTEGER  :: stat                ! status
  END TYPE rt_type

  TYPE(srt_type), PARAMETER :: srt_init = srt_type(.FALSE., 'noname')
  TYPE(srt_type) :: srt(timer_max)

  TYPE(rt_type), PARAMETER :: rt_init = rt_type( &
#ifdef __ibm__ 
        (/ 0, 0, 0, 0/),  & ! mark1
#else
        0.0_dp, &           ! mark1
#endif
        0.0_dp, &           ! tot
        HUGE(0.0_dp), &     ! min
        0.0_dp, &           ! max
        0.0_dp, &           ! last
        0, &                ! call_n
        rt_undef_stat)      ! stat

  TYPE(rt_type) :: rt(timer_max)

#if defined (__SX__) || defined (ES)
! necessary due to a bug in NECs OMP implementation
LOCAL COMMON /th_real_timer/ rt
#else
#ifdef _OPENMP
  COMMON /th_real_timer/ rt
!$omp threadprivate(/th_real_timer/)
#endif
#endif

  LOGICAL :: need_init = .TRUE.

  INTERFACE 
    SUBROUTINE util_init_real_time()
      IMPLICIT NONE
    END SUBROUTINE util_init_real_time

    SUBROUTINE util_get_real_time_size(sz)
      IMPLICIT NONE
      INTEGER, INTENT(out) :: sz           
    END SUBROUTINE util_get_real_time_size

#if defined (__ibm__)
    SUBROUTINE util_read_real_time(t)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: t(*)
    END SUBROUTINE util_read_real_time

    SUBROUTINE util_diff_real_time(t1,t2,dt)
      USE mo_kind, ONLY: dp
      IMPLICIT NONE
      INTEGER, INTENT(in)  :: t1(*), t2(*)
      REAL(dp), INTENT(out) :: dt
    END SUBROUTINE util_diff_real_time
#else
    SUBROUTINE util_read_real_time(t)
      USE mo_kind, ONLY: dp
      IMPLICIT NONE
      REAL(dp), INTENT(in) :: t
    END SUBROUTINE util_read_real_time

    SUBROUTINE util_diff_real_time(t1,t2,dt)
      USE mo_kind, ONLY: dp
      IMPLICIT NONE
      REAL(dp), INTENT(in)  :: t1, t2
      REAL(dp), INTENT(out) :: dt
    END SUBROUTINE util_diff_real_time
#endif

  END INTERFACE

CONTAINS

  SUBROUTINE init
    INTEGER :: sz

    INTEGER   :: ii   = 0
    REAL (dp) :: dd = 0.0_dp
    INTEGER :: io_size, integer_byte_size , integer_io_size, realdp_byte_size

!#ifdef _OPENMP
!    IF (omp_get_dynamic()) THEN
!      CALL message ('mo_real_timer (init)', &
!           'OMP dynamic is true - does not work with this program')
!      CALL message ('mo_real_timer (init)', &
!           'Set OMP dynamic to false')
!      CALL omp_set_dynamic(.FALSE.)
!    END IF
!#endif

    CALL util_init_real_time()
    CALL util_get_real_time_size(sz)

#if defined (__ibm__)
    IF (BIT_SIZE(rt(1)%mark1)*SIZE(rt(1)%mark1) < sz*8) &
         CALL real_timer_abort(0,'buffer size for time stamps too small')
#else
    integer_byte_size = BIT_SIZE(ii)/8
    INQUIRE (iolength=io_size) ii
    integer_io_size = io_size
    INQUIRE (iolength=io_size) dd
    realdp_byte_size = io_size/integer_io_size*integer_byte_size
    IF (realdp_byte_size < sz) &
         CALL real_timer_abort(0,'buffer size for time stamps too small')
#endif
    
    need_init = .FALSE.

    CALL estimate_overhead
  END SUBROUTINE init

!----

  SUBROUTINE estimate_overhead    
    INTEGER, PARAMETER :: n = 100 ! tests need about n microsecs on pwr4
    REAL(dp) :: dt_min

    tm_shift = 0.0_dp
    CALL m1(dt_min)
    tm_shift = dt_min

    CALL m2(dt_min)
    tm_overhead = dt_min

  CONTAINS

    SUBROUTINE m1(dt0)
      REAL(dp), INTENT(out) :: dt0
      TYPE(time_mark_type) :: mark
      REAL(dp) :: dt
      INTEGER :: i

      dt0 = 1.0_dp
      DO i = 1, n
        CALL set_time_mark(mark)
        dt = get_time_val(mark)
        IF (dt < dt0) dt0 = dt
      ENDDO      

    END SUBROUTINE m1

    SUBROUTINE m2(dt0)
      REAL(dp), INTENT(out) :: dt0
      TYPE(time_mark_type) :: mark1, mark2
      REAL(dp) :: dt1, dt2
      INTEGER :: i

      dt0 = 1.0_dp
      DO i = 1, n
        CALL set_time_mark(mark2)
        CALL set_time_mark(mark1)
        dt1 = get_time_val(mark1)
        dt2 = get_time_val(mark2)
        IF (dt2 < dt0) dt0 = dt2
        IF (dt2 < dt1) CALL real_timer_abort(reason='estimate_overhead:internal error')
      ENDDO

    END SUBROUTINE m2

  END SUBROUTINE estimate_overhead

!----

  SUBROUTINE timer_reset_all

    INTEGER :: it

    DO it = 1, top_timer
      rt(it)%tot    = 0.0_dp
      rt(it)%min    = 0.0_dp
      rt(it)%max    = 0.0_dp
      rt(it)%last   = 0.0_dp
      rt(it)%call_n = 0
    ENDDO

  END SUBROUTINE timer_reset_all

  SUBROUTINE timer_reset(it)
    INTEGER, INTENT(in) :: it

    IF (it < 1 .OR. it > top_timer) &
         CALL real_timer_abort(it,'timer_reset: timer out of bounds')

    rt(it)%tot    = 0.0_dp
    rt(it)%min    = 0.0_dp
    rt(it)%max    = 0.0_dp
    rt(it)%last   = 0.0_dp
    rt(it)%call_n = 0

  END SUBROUTINE timer_reset


  SUBROUTINE timer_reset_field(it_field)
    INTEGER, INTENT(in) :: it_field(:)

    INTEGER :: iit, it

    DO iit = LBOUND(it_field,1), UBOUND(it_field,1)
      it = it_field(iit)
      IF (it < 1 .OR. it > top_timer) &
           CALL real_timer_abort(it,'timer_reset_field: timer out of bounds')

      rt(it)%tot    = 0.0_dp
      rt(it)%min    = 0.0_dp
      rt(it)%max    = 0.0_dp
      rt(it)%last   = 0.0_dp
      rt(it)%call_n = 0

    ENDDO

  END SUBROUTINE timer_reset_field

!---

  INTEGER FUNCTION new_timer(text)
    CHARACTER(len=*), INTENT(in), OPTIONAL :: text

    INTEGER::jt

#ifdef _OPENMP
    IF ( omp_in_parallel() ) &
         CALL real_timer_abort(0,'new_timer called in parallel region')
#endif

    IF (need_init) CALL init
    top_timer = top_timer+1

    IF (top_timer > timer_max) THEN
       CALL message('new_timer','list of timers:')
       DO jt = 1, timer_max
          WRITE (nerr,*)  jt, srt(jt)%text
       ENDDO
       CALL message ('new_timer','timer_max is too small')
       CALL real_timer_abort(jt,'new_timer failed')
    ENDIF
    srt(top_timer) = srt_init
    srt(top_timer)%reserved = .TRUE.
    IF (PRESENT(text)) srt(top_timer)%text = text

!$omp parallel
    rt(top_timer) = rt_init
    new_timer = top_timer
!$omp end parallel

  END FUNCTION new_timer

!---

  SUBROUTINE del_timer(it)
    INTEGER, INTENT(in) :: it
    
    INTEGER :: jt

#ifdef _OPENMP
    IF ( omp_in_parallel() ) &
         CALL real_timer_abort(0,'del_timer called in parallel region')
#endif

    srt(it)%reserved = .FALSE.
    DO jt = top_timer, 1, -1
       IF (srt(jt)%reserved) EXIT
       top_timer = jt-1
    ENDDO

  END SUBROUTINE del_timer

!---

  SUBROUTINE timer_start(it)
    INTEGER, INTENT(in) :: it

    IF (it < 1 .OR. it > top_timer) &
         CALL real_timer_abort(it,'timer_start: timer out of bounds')
    IF (rt(it)%stat == rt_on_stat) &
         CALL real_timer_abort(it,'timer_start: timer_stop call missing')

    call util_read_real_time(rt(it)%mark1)

    rt(it)%stat = rt_on_stat

  END SUBROUTINE timer_start

!---

  SUBROUTINE set_time_mark(mark0)
    TYPE(time_mark_type), INTENT(out) :: mark0

    ! simple timer - no statistics

    CALL util_read_real_time(mark0%t)

  END SUBROUTINE set_time_mark

  REAL(dp) FUNCTION get_time_val(mark0)
    TYPE(time_mark_type), INTENT(in) :: mark0

    ! simple timer - no statistics

    REAL(dp) :: dt

    TYPE(time_mark_type)::mark

    CALL util_read_real_time(mark%t)
    CALL util_diff_real_time(mark0%t,mark%t,dt)

    get_time_val = dt-tm_shift
    
  END FUNCTION get_time_val

!---

  REAL(dp) FUNCTION timer_val(it)
    INTEGER, INTENT(in) :: it
#if defined (__ibm__)
    INTEGER :: mark2(4)
#else
    REAL(dp) :: mark2
#endif
    REAL(dp) :: dt

    IF (it < 1 .OR. it > top_timer) &
         CALL real_timer_abort(it,'timer_val: invalid timer id')

    timer_val = rt(it)%tot
    IF (rt(it)%stat == rt_on_stat) THEN
      CALL util_read_real_time(mark2)
      CALL util_diff_real_time(rt(it)%mark1,mark2,dt)
      timer_val = timer_val+dt-tm_shift
    ENDIF
    
  END FUNCTION timer_val


  REAL(dp) FUNCTION timer_average(it)
    INTEGER, INTENT(in) :: it

    IF (it < 1 .OR. it > top_timer) &
         CALL real_timer_abort(it,'timer_average: invalid timer id')
    IF (rt(it)%stat == rt_on_stat) &
         CALL real_timer_abort(it,'timer_average: timer still active')
    IF (rt(it)%call_n == 0) &
         CALL real_timer_abort(it,'timer_average: timer never called')

    timer_average = rt(it)%tot/rt(it)%call_n

  END FUNCTION timer_average

!---

  REAL(dp) FUNCTION timer_last(it)
    INTEGER, INTENT(in) :: it

    IF (it < 1 .OR. it > top_timer) &
         CALL real_timer_abort(it,'timer_last: invalid timer id')

    timer_last = rt(it)%last

  END FUNCTION timer_last

!---

  INTEGER FUNCTION timer_count(it)
    INTEGER, INTENT(in) :: it

    IF (it < 1 .OR. it > top_timer) &
         CALL real_timer_abort(it,'timer_count: invalid timer id')

    timer_count = rt(it)%call_n

  END FUNCTION timer_count

!---

  SUBROUTINE timer_stop(it)
    INTEGER, INTENT(in) :: it
#if defined (__ibm__)
    INTEGER :: mark2(4)
#else
    REAL(dp) :: mark2
#endif
    REAL(dp) :: dt

    IF (it < 1 .OR. it > top_timer) &
         CALL real_timer_abort(it,'timer_stop: invalid timer id')

    IF (rt(it)%stat /= rt_on_stat) THEN
      IF (rt(it)%stat == rt_off_stat) THEN
        CALL real_timer_abort(it,'timer_stop: timer_start call missing')
      ELSE
        CALL real_timer_abort(it,'timer_stop: undefined timer')
      ENDIF
    ENDIF

    CALL util_read_real_time(mark2)
    CALL util_diff_real_time(rt(it)%mark1,mark2,dt)

    dt = dt-tm_shift
    rt(it)%last = dt
    rt(it)%tot = rt(it)%tot + dt
    rt(it)%call_n = rt(it)%call_n+1
    IF (dt < rt(it)%min) rt(it)%min = dt
    IF (dt > rt(it)%max) rt(it)%max = dt
    rt(it)%stat = rt_off_stat

  END SUBROUTINE timer_stop

!---

  SUBROUTINE timer_report(ithread,itimer)
    INTEGER, INTENT(in), OPTIONAL :: ithread     ! show this thread if present
    INTEGER, INTENT(in), OPTIONAL :: itimer      ! show this timer if present

    INTEGER :: tid, it1
#ifdef _OPENMP
    INTEGER :: itid
#endif

#ifndef NOMPI
    INTEGER :: ibuf(2)
#endif 

#ifdef _OPENMP
    IF ( omp_in_parallel() ) &
         CALL real_timer_abort(0,'timer_report called in parallel region')
#endif

    IF (PRESENT(itimer)) THEN
      it1 = itimer
    ELSE
      it1 = -1 ! report of all timers
    ENDIF

#ifndef NOMPI
    ! order mpi:
    IF (p_pe > 0) THEN
      CALL p_recv(ibuf(1), p_pe-1, 12345)
    ENDIF
#endif

    ! order omp:
!$omp parallel private(itid,tid)

!$omp do ordered
#ifdef _OPENMP
    DO itid = 1, omp_get_num_threads()
      tid = omp_get_thread_num()
#else
      tid = 1
#endif
!$omp ordered
      IF (PRESENT(ithread)) THEN
        IF ( tid == ithread  ) CALL my_report0(it1)
      ELSE
        CALL my_report0(it1)
      ENDIF
!$omp flush
!$omp end ordered
#ifdef _OPENMP
    ENDDO
#endif
!$omp end do

!$omp end parallel

#ifndef NOMPI
    IF (p_pe < p_nprocs-1) THEN
      CALL p_send(ibuf(1), p_pe+1, 12345)
    ENDIF
    CALL p_barrier(p_all_comm) 
#endif

  END SUBROUTINE timer_report

  SUBROUTINE my_report0(it1)
    INTEGER, INTENT(in) :: it1
    INTEGER :: it
    LOGICAL :: any_active

    IF (it1 > 0) THEN 
      IF (rt(it1)%stat /= rt_undef_stat) CALL my_report1(it1)
    ELSE
      
!$omp master
      WRITE (nerr,*)
      WRITE (nerr,separator)
      WRITE (nerr,separator)
      WRITE (nerr,'(a)') 'Timer report: '
      !     0        1         2         3         4         5         6         7         8
      !    '12345678901234567890123456789012345678901234567890123456789012345678901234567890'
#ifndef NOMPI
      IF (p_nprocs > 1) THEN
        WRITE (nerr,'(a,i4,a)') &
           'PE: ', p_pe, '                 calls  t_min       t_average   t_max       t_total    '
      ELSE
        WRITE (nerr,'(a)') &
           '                         calls  t_min       t_average   t_max       t_total    '
      END IF
#else
        WRITE (nerr,'(a)') &
           '                         calls  t_min       t_average   t_max       t_total    '
#endif
      WRITE (nerr,separator)
!$omp end master
        
      any_active = .FALSE.
      DO it = 1, top_timer
        IF (rt(it)%stat /= rt_undef_stat) THEN
          any_active = .TRUE.
          CALL my_report1(it)
        ENDIF
      ENDDO
      IF (any_active) &
           WRITE (nerr,separator)
    ENDIF

  END SUBROUTINE my_report0

  SUBROUTINE my_report1(it)
    INTEGER, INTENT(in) :: it
    
    REAL(dp) :: avg, total
    CHARACTER(len=12) :: min_str, avg_str, max_str, tot_str
    INTEGER :: tid
    
#ifdef _OPENMP
    tid = omp_get_thread_num()
#else
    tid = -1
#endif    

    total = timer_val(it)
    
    avg = rt(it)%tot/MAX(1,rt(it)%call_n)
    
    IF ( rt(it)%call_n > 0 ) THEN
      min_str=time_str(rt(it)%min)
      avg_str=time_str(avg)
      max_str=time_str(rt(it)%max)
    ELSE
      min_str=''
      avg_str=''
      max_str=''
    ENDIF
    tot_str=time_str(total)
    IF (tid == -1) THEN
      WRITE (nerr,'(a22,i6,4a12)') '    '//srt(it)%text, &
           rt(it)%call_n,min_str,avg_str,max_str,tot_str
    ELSE
      WRITE (nerr,'(i2,a20,i6,4a12)') tid,': '//srt(it)%text, &
           rt(it)%call_n,min_str,avg_str,max_str,tot_str
    ENDIF

  END SUBROUTINE my_report1

  CHARACTER(len=12) FUNCTION time_str(ts)
    REAL(dp), INTENT(in) :: ts
    
    REAL(dp) :: rest

    INTEGER :: d, h, m, s

    CHARACTER(len=2) :: d_str, h_str, m_str, s_str
    CHARACTER(len=12) :: x
    
    rest = ts
    
    d = INT(rest/(3600*24))
    rest = rest-d*(3600*24)
    IF (d > 99) THEN
      x = '>99d'
      time_str = ADJUSTR(x)
      RETURN
    ENDIF
    WRITE(d_str,'(I2.2)') d
    
    h = INT(rest/3600)
    rest = rest-h*3600
    WRITE(h_str,'(I2.2)') h
    
    m = INT(rest/60)
    rest = rest-m*60
    WRITE(m_str,'(I2.2)') m
    
    s = INT(rest)
    rest = rest-s
    WRITE(s_str,'(I2.2)') s
    
    IF (d > 0) THEN        
      x = TRIM(d_str)//'d'//TRIM(h_str)//'h'
    ELSEIF (h > 0) THEN
      x = TRIM(h_str)//'h'//TRIM(m_str)//'m'
    ELSEIF (m > 0) THEN
      x = TRIM(m_str)//'m'//TRIM(s_str)//'s'
    ELSEIF (ts >= 1.0_dp) THEN
      WRITE(x,'(F6.3,A)') ts, 's'
    ELSE
      WRITE(x,'(F6.5,A)') ts, 's'
    ENDIF
    time_str = ADJUSTR(x)
    
  END FUNCTION time_str
  
!---

  SUBROUTINE exercise_timer

    ! Try to detect possible interferences with other tasks running
    ! on the same node. In this case, real-time measurement makes no sense.

    INTEGER ::  tn, v
    REAL(dp) :: t1, t2
    REAL(dp) :: dt
    REAL(dp), ALLOCATABLE :: buf(:)
#ifdef _OPENMP
    INTEGER :: tid
    REAL(dp) :: rmin, rmax, thread_sum, tavg
    REAL(dp), ALLOCATABLE :: thread_dt(:),hybrid_dt(:,:)
#endif
#ifdef NOMPI
    INTEGER :: p_nprocs = 1
#endif

    IF (need_init) CALL init

#ifdef _OPENMP
    IF (omp_in_parallel()) &
         CALL real_timer_abort(0,'exercise_timer called in parallel region')

    tn = omp_get_max_threads()
#else
    tn = 1
#endif

    ALLOCATE(buf(0:p_nprocs-1))
#ifdef _OPENMP
    ALLOCATE(hybrid_dt(0:tn-1,0:p_nprocs-1),thread_dt(0:tn-1))
    thread_dt(:) = 0.0_dp
    hybrid_dt(:,:) = 0.0_dp

    thread_sum = 0
#endif
!$omp parallel private(tid,t1,t2,dt,v),shared(thread_sum,thread_dt)
#ifdef _OPENMP
    tid = omp_get_thread_num()
#endif
    t1 = util_walltime()
    v = 1
    CALL hom_cpu_work(v,p_nprocs*tn)
    t2 = util_walltime()
    dt = t2-t1
!    IF (v < 0) PRINT*,'v=',v ! branch never reached
#ifdef _OPENMP
    thread_dt(tid) = dt
#endif
!$omp end parallel

#if (! defined NOMPI) && (defined _OPENMP)
    CALL p_gather (thread_dt, hybrid_dt, 0)

    IF (p_pe == 0) THEN
       tavg = SUM(hybrid_dt)/(p_nprocs*tn)
       rmin = MINVAL(hybrid_dt)/tavg
       rmax = MAXVAL(hybrid_dt)/tavg
       IF (rmax-rmin > 0.2_dp) THEN
          WRITE (nerr,'(a,2f10.4)') &
               'exercise_timer: WARNING - unexpected inbalance of CPU share - rmin, rmax: ', &
               rmin, rmax
          CALL util_system("ps -ef >&2")
       ELSE
          WRITE  (nerr,'(a,2f10.4)') 'exercise_timer: CPU balance ok - rmin, rmax: ', &
               rmin, rmax
       ENDIF
    ENDIF
    
    DEALLOCATE(buf,thread_dt,hybrid_dt)
#endif

  CONTAINS

    SUBROUTINE hom_cpu_work(val, ttn)
      INTEGER, INTENT(inout) :: val
      INTEGER, INTENT(in) :: ttn
      INTEGER :: i, j, k

      ! some cpu-load
      DO i=1,1024
         DO j=1,100
            DO k=1,100
               val=MOD(val+i+j+k,19)+ttn
            ENDDO
         ENDDO
      ENDDO
      
    END SUBROUTINE hom_cpu_work
    

  END SUBROUTINE exercise_timer

!---

  SUBROUTINE real_timer_abort(it,reason)
    INTEGER, OPTIONAL, INTENT(in) :: it
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: reason    

    WRITE (message_text,*)  'error in module mo_real_timer:'
    CALL message ('', TRIM(message_text))
    IF (PRESENT(it)) THEN
      WRITE (message_text,*) 'timer handle: ', it
      CALL message ('', TRIM(message_text))
      IF (it < 1 .OR. it > top_timer) THEN
        WRITE (message_text,*) 'timer name: unspecified'
      CALL message ('', TRIM(message_text))
      ELSE
        WRITE (message_text,*) 'timer name: ', TRIM(srt(it)%text)      
      CALL message ('', TRIM(message_text))
      ENDIF
    ENDIF
    IF (PRESENT(reason)) THEN
      WRITE (message_text,*) '            ', reason
      CALL message ('', TRIM(message_text))
    ENDIF

    CALL finish('real_timer_abort','Abort')

  END SUBROUTINE real_timer_abort

END MODULE mo_real_timer
