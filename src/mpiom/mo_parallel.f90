MODULE mo_parallel

  USE mo_kind,   ONLY: dp 
  USE mo_param1, ONLY: ke, ie, je, ie_g, je_g, icycli
  USE mo_mpi

  IMPLICIT NONE

  INTEGER, PARAMETER :: maxproc = 1024

  ! Number of subdivisions in x/y direction

  INTEGER :: nprocx, nprocy
  INTEGER :: nprocxy

  ! Our own offset

  INTEGER :: p_ioff, p_joff

  ! Flag if we have the boundaries

  LOGICAL :: have_g_is, have_g_ie, have_g_js, have_g_je

  ! Start of the inner domains in x/y direction
  ! p_lim_x(0) = 2
  ! p_lim_x(i) = start of inner domain i
  ! p_lim_x(nprocx) = ie_g

  INTEGER, PRIVATE :: p_lim_x(0:maxproc), p_lim_y(0:maxproc)

  ! For every processor: number in x/y direction (0-based)

  INTEGER, PRIVATE :: p_num_x(0:maxproc), p_num_y(0:maxproc)

  ! Global offsets and sizes for each processor (both for outer domains)

  INTEGER, PRIVATE :: p_ioff_g(0:maxproc), p_joff_g(0:maxproc)
  INTEGER, PRIVATE :: p_size_x(0:maxproc), p_size_y(0:maxproc)

  REAL, PRIVATE :: t2d=0.0_dp, t3d=0.0_dp, ts, te
  INTEGER,  PRIVATE :: n2d=0, n3d=0, count=0

  INTERFACE bounds_exch
    MODULE PROCEDURE bounds_exch_2d
    MODULE PROCEDURE bounds_exch_3d
    MODULE PROCEDURE bounds_exch_hh_2d
    MODULE PROCEDURE bounds_exch_hh_3d
  END INTERFACE

  INTERFACE para_check
    MODULE PROCEDURE para_check_2d
    MODULE PROCEDURE para_check_3d
  END INTERFACE

  INTERFACE global_sum
    MODULE PROCEDURE global_sum_i
    MODULE PROCEDURE global_sum_r
    MODULE PROCEDURE global_sum_1d
    MODULE PROCEDURE global_sum_2d
  END INTERFACE

CONTAINS

  !-----------------------------------------------------------------------

  SUBROUTINE p_deco

    !   domain decomposition

    INTEGER :: i, nx, ny
#ifdef bounds_exch_tp
    INTEGER :: ii
#endif

    IF(p_nprocs > 1) THEN

#ifdef DEBUG
      WRITE(nerr,*) 'Process ',p_pe,' of ',p_nprocs,' is alive'
#endif

      ! set some variables

      IF(p_pe==0) THEN

        ! nprocx and nprocy must be set by the calling process

        IF(nprocx==0 .OR. nprocy==0) THEN
          WRITE(nerr,*) 'ERROR: nprocx or nprocy not set'
          CALL p_abort
        ENDIF

#ifdef bounds_exch_tp
        IF (mod(nprocx,2) /= 0 .AND. nprocx  /= 1)  THEN
           WRITE(nerr,*) 'ERROR: for the tripolar setup  nprocx should be even (or one) '
          CALL p_abort
        ENDIF
         
#endif 

        nprocxy = nprocx*nprocy

        IF(nprocxy /= p_nprocs .AND. nprocxy /= p_nprocs-1) THEN
          WRITE(nerr,*)'Number of processors = ',p_nprocs
          WRITE(nerr,*)'nprocx = ',nprocx,' nprocy = ',nprocy
          WRITE(nerr,*)'Number of processors doesnt fit!'
          CALL p_abort
        ENDIF

        IF(((ie_g-2)/nprocx)<3 .OR. ((je_g-2)/nprocy)<3) THEN
          WRITE(nerr,*)'Decomposed domain gets too small'
          WRITE(nerr,*)'We need at least 3 rows in every direction'
          CALL p_abort
        ENDIF

      ENDIF

      ! broadcast nprocx and nprocy

      CALL p_bcast(nprocx,0)
      CALL p_bcast(nprocy,0)
      nprocxy = nprocx*nprocy

      ! Decomposition - compute domain limits
#ifdef bounds_exch_tp
      DO i=0,(nprocx/2)
        p_lim_x(i) = 2 + i*(ie_g-2)/nprocx
        ii=nprocx-i
        p_lim_x(ii) = (ie_g - (p_lim_x(i)-2))
      ENDDO
#else
      DO i=0,nprocx
        p_lim_x(i) = 2 + i*(ie_g-2)/nprocx
      ENDDO
#endif
      DO i=0,nprocy
        p_lim_y(i) = 2 + i*(je_g-2)/nprocy
      ENDDO

      ! Set number of processors in x and y direction

      DO i=0,nprocx-1
        p_num_x(i:nprocxy-1:nprocx) = i
      ENDDO
      DO i=0,nprocy-1
        p_num_y(i*nprocx:(i+1)*nprocx-1) = i
      ENDDO

      ! Offsets and sizes

      DO i=0,nprocxy-1
        nx = p_num_x(i)
        ny = p_num_y(i)

        p_ioff_g(i) = p_lim_x(nx)-2
        p_joff_g(i) = p_lim_y(ny)-2
        p_size_x(i) = p_lim_x(nx+1) - p_lim_x(nx) + 2
        p_size_y(i) = p_lim_y(ny+1) - p_lim_y(ny) + 2
      ENDDO

      ! Get our own values

      IF(p_pe<nprocxy) THEN
        ie = p_size_x(p_pe)
        je = p_size_y(p_pe)
        p_ioff = p_ioff_g(p_pe)
        p_joff = p_joff_g(p_pe)
        have_g_is = p_num_x(p_pe) == 0
        have_g_ie = p_num_x(p_pe) == nprocx-1
        have_g_js = p_num_y(p_pe) == 0
        have_g_je = p_num_y(p_pe) == nprocy-1
      ELSE
        ie = ie_g
        je = je_g
        p_ioff = 0
        p_joff = 0
        have_g_is = .TRUE.
        have_g_ie = .TRUE.
        have_g_js = .TRUE.
        have_g_je = .TRUE.
      ENDIF

#ifdef DEBUG
      WRITE(nerr,'(a,i2,a,2i4,a,2i4)')'Proc ',p_pe,' offset: ',p_ioff,p_joff, &
           ' Size: ',ie,je
#endif

    ELSE

      nprocx = 1
      nprocy = 1
      nprocxy = 1
      p_lim_x(0) = 2
      p_lim_x(1) = ie_g
      p_lim_y(0) = 2
      p_lim_y(1) = je_g
      p_num_x(0) = 0
      p_num_y(0) = 0
      p_ioff_g(0) = 0
      p_joff_g(0) = 0
      p_size_x(0) = ie_g
      p_size_y(0) = je_g
      ie = ie_g
      je = je_g
      p_ioff = 0
      p_joff = 0
      have_g_is = .TRUE.
      have_g_ie = .TRUE.
      have_g_js = .TRUE.
      have_g_je = .TRUE.

#ifdef DEBUG
      WRITE(nerr,*) 'Running on a single processor'
#endif

    ENDIF

  END SUBROUTINE p_deco

  !-----------------------------------------------------------------------

  SUBROUTINE gather_arr(arrl, arrg, pe)

    ! Gathers all local array parts into a global array on PE pe

    REAL,    INTENT(IN)  :: arrl(:,:)
    REAL,    INTENT(OUT) :: arrg(:,:)
    INTEGER, INTENT(IN)  :: pe

    INTEGER :: n, iis, iie, jjs, jje, ioff, joff
    REAL, ALLOCATABLE :: aux(:,:)

    IF(p_pe/=pe) THEN
      IF(p_pe<nprocxy) CALL p_send(arrl,pe,1111)
    ELSE
      DO n=0,nprocxy-1
        iis = 1
        iie = p_size_x(n)
        jjs = 1
        jje = p_size_y(n)
        ALLOCATE(aux(iie,jje))
        IF (n==pe) THEN
          aux = arrl
        ELSE
          CALL p_recv(aux,n,1111)
        ENDIF
        ! Copy only outer limits into arrg
        IF(p_num_x(n) /= 0        ) iis = iis+1
        IF(p_num_x(n) /= nprocx-1 ) iie = iie-1
        IF(p_num_y(n) /= 0        ) jjs = jjs+1
        IF(p_num_y(n) /= nprocy-1 ) jje = jje-1
        ioff = p_ioff_g(n)
        joff = p_joff_g(n)
        arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje) = aux(iis:iie,jjs:jje)
        DEALLOCATE(aux)
      ENDDO
    ENDIF

  END SUBROUTINE gather_arr

  SUBROUTINE gather3_arr(arrl, arrg, pe)

    ! Gathers all local array parts into a global array on PE pe

    REAL,    INTENT(IN)  :: arrl(:,:,:)
    REAL,    INTENT(OUT) :: arrg(:,:,:)
    INTEGER, INTENT(IN)  :: pe

    INTEGER :: n, iis, iie, jjs, jje, ioff, joff
    REAL, ALLOCATABLE :: aux(:,:,:)

    IF(p_pe/=pe) THEN
      IF(p_pe<nprocxy) CALL p_send(arrl,pe,1111)
    ELSE
      DO n=0,nprocxy-1
        iis = 1
        iie = p_size_x(n)
        jjs = 1
        jje = p_size_y(n)
        ALLOCATE(aux(iie,jje,ke))
        IF (n==pe) THEN
          aux = arrl
        ELSE
          CALL p_recv(aux,n,1111)
        ENDIF
        ! Copy only outer limits into arrg
        IF(p_num_x(n) /= 0        ) iis = iis+1
        IF(p_num_x(n) /= nprocx-1 ) iie = iie-1
        IF(p_num_y(n) /= 0        ) jjs = jjs+1
        IF(p_num_y(n) /= nprocy-1 ) jje = jje-1
        ioff = p_ioff_g(n)
        joff = p_joff_g(n)
        arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje,1:ke) = aux(iis:iie,jjs:jje,1:ke)
        DEALLOCATE(aux)
      ENDDO
    ENDIF

  END SUBROUTINE gather3_arr

  SUBROUTINE gather(arrl, arrg, pe)

    ! Gathers all local array parts into a global array on PE pe

    REAL,    INTENT(IN)  :: arrl(:,:)
    REAL,    POINTER     :: arrg(:,:)
    INTEGER, INTENT(IN)  :: pe

    INTEGER :: n, iis, iie, jjs, jje, ioff, joff
    REAL, ALLOCATABLE :: aux(:,:)

    IF(p_pe/=pe) THEN
      IF(p_pe<nprocxy) CALL p_send(arrl,pe,1111)
    ELSE
      DO n=0,nprocxy-1
        iis = 1
        iie = p_size_x(n)
        jjs = 1
        jje = p_size_y(n)
        ALLOCATE(aux(iie,jje))
        IF (n==pe) THEN
          aux = arrl
        ELSE
          CALL p_recv(aux,n,1111)
        ENDIF
        ! Copy only outer limits into arrg
        IF(p_num_x(n) /= 0        ) iis = iis+1
        IF(p_num_x(n) /= nprocx-1 ) iie = iie-1
        IF(p_num_y(n) /= 0        ) jjs = jjs+1
        IF(p_num_y(n) /= nprocy-1 ) jje = jje-1
        ioff = p_ioff_g(n)
        joff = p_joff_g(n)
        arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje) = aux(iis:iie,jjs:jje)
        DEALLOCATE(aux)
      ENDDO
    ENDIF

  END SUBROUTINE gather
  !-----------------------------------------------------------------------

  SUBROUTINE scatter_arr(arrg, arrl, pe)

    ! Scatters global array (arrg) on PE pe to the local array (arrl) on all PEs
    ! Since this routine is not used in perfomance critical parts
    ! we use the simple broadcast version

    REAL,    INTENT(INOUT) :: arrg(:,:)
    REAL,    INTENT(OUT)   :: arrl(:,:)
    INTEGER, INTENT(IN)    :: pe

    CALL p_bcast(arrg,pe)

    arrl(:,:) = arrg(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)

  END SUBROUTINE scatter_arr

  SUBROUTINE scatter3_arr(arrg, arrl, pe)

    ! Scatters global array (arrg) on PE pe to the local array (arrl) on all PEs
    ! Since this routine is not used in perfomance critical parts
    ! we use the simple broadcast version

    REAL,    INTENT(INOUT) :: arrg(:,:,:)
    REAL,    INTENT(OUT)   :: arrl(:,:,:)
    INTEGER, INTENT(IN)    :: pe

    CALL p_bcast(arrg,pe)

    arrl(:,:,:) = arrg(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je,:)

  END SUBROUTINE scatter3_arr


  SUBROUTINE scatter(arrg, arrl, pe)

    ! Scatters global array (arrg) on PE pe to the local array (arrl) on all PEs
    ! Since this routine is not used in perfomance critical parts
    ! we use the simple broadcast version

    REAL,    INTENT(INOUT) :: arrg(:,:)
    REAL,    INTENT(OUT)   :: arrl(:,:)
    INTEGER, INTENT(IN)    :: pe

    CALL p_bcast(arrg,pe)

    arrl(:,:) = arrg(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)

  END SUBROUTINE scatter

  !-----------------------------------------------------------------------

  SUBROUTINE bounds_exch_2d(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,text)

    ! Exchanges boundaries of 2D arrays
   ! Since 2D boundary exchange has an impact on performance,
    ! we allow that up to 10 2D arrays are exchanged at the same time

#ifdef bounds_exch_put
    USE mo_boundsexch
#endif
    REAL, INTENT(INOUT) :: a0(:,:)
    REAL, INTENT(INOUT), OPTIONAL :: a1(:,:),a2(:,:),a3(:,:),a4(:,:),a5(:,:)
    REAL, INTENT(INOUT), OPTIONAL :: a6(:,:),a7(:,:),a8(:,:),a9(:,:)
    CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text

    INTEGER nm, np, n

#ifndef bounds_exch_put
    REAL xr1(je,10),xr2(je,10),yr1(ie,10),yr2(ie,10)
#endif
    REAL xs1(je,10),xs2(je,10),ys1(ie,10),ys2(ie,10)

    ts = p_time()

    ! x-direction

    xs1(:,1) = a0(2,:)
    xs2(:,1) = a0(ie-1,:)
    n = 1
    IF(PRESENT(a1)) THEN; n=n+1; xs1(:,n) = a1(2,:); xs2(:,n) = a1(ie-1,:); ENDIF
    IF(PRESENT(a2)) THEN; n=n+1; xs1(:,n) = a2(2,:); xs2(:,n) = a2(ie-1,:); ENDIF
    IF(PRESENT(a3)) THEN; n=n+1; xs1(:,n) = a3(2,:); xs2(:,n) = a3(ie-1,:); ENDIF
    IF(PRESENT(a4)) THEN; n=n+1; xs1(:,n) = a4(2,:); xs2(:,n) = a4(ie-1,:); ENDIF
    IF(PRESENT(a5)) THEN; n=n+1; xs1(:,n) = a5(2,:); xs2(:,n) = a5(ie-1,:); ENDIF
    IF(PRESENT(a6)) THEN; n=n+1; xs1(:,n) = a6(2,:); xs2(:,n) = a6(ie-1,:); ENDIF
    IF(PRESENT(a7)) THEN; n=n+1; xs1(:,n) = a7(2,:); xs2(:,n) = a7(ie-1,:); ENDIF
    IF(PRESENT(a8)) THEN; n=n+1; xs1(:,n) = a8(2,:); xs2(:,n) = a8(ie-1,:); ENDIF
    IF(PRESENT(a9)) THEN; n=n+1; xs1(:,n) = a9(2,:); xs2(:,n) = a9(ie-1,:); ENDIF

    IF (nprocx>1 .AND. p_pe<nprocxy ) THEN
      ! Get processor numbers of neighbors
      nm = MERGE(p_pe-1,p_pe+nprocx-1,p_num_x(p_pe)/=0)
      np = MERGE(p_pe+1,p_pe-nprocx+1,p_num_x(p_pe)/=nprocx-1)
      
#ifdef bounds_exch_isend
      count = je*n
      CALL p_isend(xs1(1,1),nm,1,p_count=count)
      CALL p_isend(xs2(1,1),np,2,p_count=count)
      CALL p_irecv(xr1(1,1),np,1,p_count=count)
      CALL p_irecv(xr2(1,1),nm,2,p_count=count)
      CALL p_wait
#else
#ifdef bounds_exch_put
      CALL p_win_fence(p_win2dx)
      CALL p_put(xs1(:,1:n),nm,p_win2dx(1))
      CALL p_put(xs2(:,1:n),np,p_win2dx(2))
      CALL p_win_fence(p_win2dx)
#else
      CALL p_sendrecv(xs1(:,1:n),nm,xr1(:,1:n),np,1)
      CALL p_sendrecv(xs2(:,1:n),np,xr2(:,1:n),nm,2)
#endif
#endif
    ELSE
      xr1(:,1:n) = xs1(:,1:n)
      xr2(:,1:n) = xs2(:,1:n)
    ENDIF
    
    IF(icycli/=0 .OR. .NOT. have_g_is) THEN
      a0( 1,:) = xr2(:,1)
      n = 1
      IF(PRESENT(a1)) THEN; n=n+1; a1( 1,:) = xr2(:,n); ENDIF
      IF(PRESENT(a2)) THEN; n=n+1; a2( 1,:) = xr2(:,n); ENDIF
      IF(PRESENT(a3)) THEN; n=n+1; a3( 1,:) = xr2(:,n); ENDIF
      IF(PRESENT(a4)) THEN; n=n+1; a4( 1,:) = xr2(:,n); ENDIF
      IF(PRESENT(a5)) THEN; n=n+1; a5( 1,:) = xr2(:,n); ENDIF
      IF(PRESENT(a6)) THEN; n=n+1; a6( 1,:) = xr2(:,n); ENDIF
      IF(PRESENT(a7)) THEN; n=n+1; a7( 1,:) = xr2(:,n); ENDIF
      IF(PRESENT(a8)) THEN; n=n+1; a8( 1,:) = xr2(:,n); ENDIF
      IF(PRESENT(a9)) THEN; n=n+1; a9( 1,:) = xr2(:,n); ENDIF
    ENDIF
    IF(icycli/=0 .OR. .NOT. have_g_ie) THEN
      a0(ie,:) = xr1(:,1)
      n = 1
      IF(PRESENT(a1)) THEN; n=n+1; a1(ie,:) = xr1(:,n); ENDIF
      IF(PRESENT(a2)) THEN; n=n+1; a2(ie,:) = xr1(:,n); ENDIF
      IF(PRESENT(a3)) THEN; n=n+1; a3(ie,:) = xr1(:,n); ENDIF
      IF(PRESENT(a4)) THEN; n=n+1; a4(ie,:) = xr1(:,n); ENDIF
      IF(PRESENT(a5)) THEN; n=n+1; a5(ie,:) = xr1(:,n); ENDIF
      IF(PRESENT(a6)) THEN; n=n+1; a6(ie,:) = xr1(:,n); ENDIF
      IF(PRESENT(a7)) THEN; n=n+1; a7(ie,:) = xr1(:,n); ENDIF
      IF(PRESENT(a8)) THEN; n=n+1; a8(ie,:) = xr1(:,n); ENDIF
      IF(PRESENT(a9)) THEN; n=n+1; a9(ie,:) = xr1(:,n); ENDIF
    ENDIF
      
    ! y-direction
    ! Note that there is no action required if nprocy==1
      
    IF (nprocy>1 .AND. p_pe<nprocxy) THEN
      
      ys1(:,1) = a0(:,2)
      ys2(:,1) = a0(:,je-1)
      n = 1
        
      IF(PRESENT(a1)) THEN; n=n+1; ys1(:,n) = a1(:,2); ys2(:,n) = a1(:,je-1); ENDIF
      IF(PRESENT(a2)) THEN; n=n+1; ys1(:,n) = a2(:,2); ys2(:,n) = a2(:,je-1); ENDIF
      IF(PRESENT(a3)) THEN; n=n+1; ys1(:,n) = a3(:,2); ys2(:,n) = a3(:,je-1); ENDIF
      IF(PRESENT(a4)) THEN; n=n+1; ys1(:,n) = a4(:,2); ys2(:,n) = a4(:,je-1); ENDIF
      IF(PRESENT(a5)) THEN; n=n+1; ys1(:,n) = a5(:,2); ys2(:,n) = a5(:,je-1); ENDIF
      IF(PRESENT(a6)) THEN; n=n+1; ys1(:,n) = a6(:,2); ys2(:,n) = a6(:,je-1); ENDIF
      IF(PRESENT(a7)) THEN; n=n+1; ys1(:,n) = a7(:,2); ys2(:,n) = a7(:,je-1); ENDIF
      IF(PRESENT(a8)) THEN; n=n+1; ys1(:,n) = a8(:,2); ys2(:,n) = a8(:,je-1); ENDIF
      IF(PRESENT(a9)) THEN; n=n+1; ys1(:,n) = a9(:,2); ys2(:,n) = a9(:,je-1); ENDIF
          
      ! Get processor numbers of neighbors
      ! For the sake of simplicity  the p_sendrecv call is as
      ! for periodic boundary conditions
      nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
      np = MOD(p_pe+nprocxy+nprocx,nprocxy)
#ifdef bounds_exch_isend
      count = ie*n
      CALL p_isend(ys1(1,1),nm,1,p_count = count)
      CALL p_isend(ys2(1,1),np,2,p_count = count)
      CALL p_irecv(yr1(1,1),np,1,p_count = count)
      CALL p_irecv(yr2(1,1),nm,2,p_count = count)
      CALL p_wait
#else
#ifdef bounds_exch_put
      CALL p_win_fence(p_win2dy)
      CALL p_put(ys1(:,1:n),nm,p_win2dy(1))
      CALL p_put(ys2(:,1:n),np,p_win2dy(2))
      CALL p_win_fence(p_win2dy)
#else
      CALL p_sendrecv(ys1(:,1:n),nm,yr1(:,1:n),np,1)
      CALL p_sendrecv(ys2(:,1:n),np,yr2(:,1:n),nm,2)
#endif
#endif
      IF(.NOT. have_g_js) THEN
        a0(:, 1) = yr2(:,1)
        n = 1
        IF(PRESENT(a1)) THEN; n=n+1; a1(:, 1) = yr2(:,n); ENDIF
        IF(PRESENT(a2)) THEN; n=n+1; a2(:, 1) = yr2(:,n); ENDIF
        IF(PRESENT(a3)) THEN; n=n+1; a3(:, 1) = yr2(:,n); ENDIF
        IF(PRESENT(a4)) THEN; n=n+1; a4(:, 1) = yr2(:,n); ENDIF
        IF(PRESENT(a5)) THEN; n=n+1; a5(:, 1) = yr2(:,n); ENDIF
        IF(PRESENT(a6)) THEN; n=n+1; a6(:, 1) = yr2(:,n); ENDIF
        IF(PRESENT(a7)) THEN; n=n+1; a7(:, 1) = yr2(:,n); ENDIF
        IF(PRESENT(a8)) THEN; n=n+1; a8(:, 1) = yr2(:,n); ENDIF
        IF(PRESENT(a9)) THEN; n=n+1; a9(:, 1) = yr2(:,n); ENDIF
      ENDIF
      IF(.NOT. have_g_je) THEN
        a0(:,je) = yr1(:,1)
        n = 1
        IF(PRESENT(a1)) THEN; n=n+1; a1(:,je) = yr1(:,n); ENDIF
        IF(PRESENT(a2)) THEN; n=n+1; a2(:,je) = yr1(:,n); ENDIF
        IF(PRESENT(a3)) THEN; n=n+1; a3(:,je) = yr1(:,n); ENDIF
        IF(PRESENT(a4)) THEN; n=n+1; a4(:,je) = yr1(:,n); ENDIF
        IF(PRESENT(a5)) THEN; n=n+1; a5(:,je) = yr1(:,n); ENDIF
        IF(PRESENT(a6)) THEN; n=n+1; a6(:,je) = yr1(:,n); ENDIF
        IF(PRESENT(a7)) THEN; n=n+1; a7(:,je) = yr1(:,n); ENDIF
        IF(PRESENT(a8)) THEN; n=n+1; a8(:,je) = yr1(:,n); ENDIF
        IF(PRESENT(a9)) THEN; n=n+1; a9(:,je) = yr1(:,n); ENDIF
      ENDIF
    ENDIF

    te = p_time()
    t2d = t2d + te-ts
    n2d = n2d + 1
    
    IF(p_nprocs > nprocxy) THEN
      ! Test mode
      IF(PRESENT(text)) THEN
        CALL para_check_2d(a0,text)
        IF(PRESENT(a1)) CALL para_check_2d(a1,text)
        IF(PRESENT(a2)) CALL para_check_2d(a2,text)
        IF(PRESENT(a3)) CALL para_check_2d(a3,text)
        IF(PRESENT(a4)) CALL para_check_2d(a4,text)
        IF(PRESENT(a5)) CALL para_check_2d(a5,text)
        IF(PRESENT(a6)) CALL para_check_2d(a6,text)
        IF(PRESENT(a7)) CALL para_check_2d(a7,text)
        IF(PRESENT(a8)) CALL para_check_2d(a8,text)
        IF(PRESENT(a9)) CALL para_check_2d(a9,text)
      ELSE
        CALL para_check_2d(a0,'bounds_exch_2d')
        IF(PRESENT(a1)) CALL para_check_2d(a1,'bounds_exch_2d')
        IF(PRESENT(a2)) CALL para_check_2d(a2,'bounds_exch_2d')
        IF(PRESENT(a3)) CALL para_check_2d(a3,'bounds_exch_2d')
        IF(PRESENT(a4)) CALL para_check_2d(a4,'bounds_exch_2d')
        IF(PRESENT(a5)) CALL para_check_2d(a5,'bounds_exch_2d')
        IF(PRESENT(a6)) CALL para_check_2d(a6,'bounds_exch_2d')
        IF(PRESENT(a7)) CALL para_check_2d(a7,'bounds_exch_2d')
        IF(PRESENT(a8)) CALL para_check_2d(a8,'bounds_exch_2d')
        IF(PRESENT(a9)) CALL para_check_2d(a9,'bounds_exch_2d')
      ENDIF
    ENDIF
        
  END SUBROUTINE bounds_exch_2d
  !-----------------------------------------------------------------------

  SUBROUTINE bounds_exch_hh_2d(ttt,a0,text)

    ! Exchanges boundaries of 2D arrays
    ! adopted to support orca type grid (hh, 01/2006)

#ifdef bounds_exch_check
    use mo_param1, only : ie,je
    INTEGER NN,jj,ii
    REAL,ALLOCATABLE :: aa0(:,:),aa1(:,:)
#endif

    REAL, INTENT(INOUT) :: a0(:,:)
    CHARACTER (LEN=*), INTENT(IN) :: ttt
    CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
    INTEGER nm, np
    REAL xr1(je),xr2(je),yr1(ie),yr2(ie)
    REAL xs1(je),xs2(je),ys1(ie),ys2(ie)

#ifdef bounds_exch_tp
    REAL ysl2(ie),yrl2(ie),ysl3(ie),yrl3(ie),ysl4(ie),yrl4(ie)
    REAL pq
    INTEGER i, ir, il, ne
#endif
#ifdef bounds_exch_check
    ALLOCATE(aa0(ie,je),aa1(ie,je))
    aa0(:,:) = a0(:,:)
#endif


    ! x-direction

    xs1(:) = a0(2,:)
    xs2(:) = a0(ie-1,:)

    IF (nprocx>1 .AND. p_pe<nprocxy ) THEN
      ! Get processor numbers of neighbors
      nm = MERGE(p_pe-1,p_pe+nprocx-1,p_num_x(p_pe)/=0)
      np = MERGE(p_pe+1,p_pe-nprocx+1,p_num_x(p_pe)/=nprocx-1)
      
      CALL p_sendrecv(xs1,nm,xr1,np,1)
      CALL p_sendrecv(xs2,np,xr2,nm,2)

    ELSE
      xr1(:)=xs1(:)
      xr2(:)=xs2(:)
    ENDIF
    
    IF(icycli/=0 .OR. .NOT. have_g_is) THEN
      a0(1,:) = xr2(:)
    ENDIF

    IF(icycli/=0 .OR. .NOT. have_g_ie) THEN
      a0(ie,:) = xr1(:)
    ENDIF
      
    ! y-direction
    ! Note that there is no action required if nprocy==1
      
    IF (nprocy>1 .AND. p_pe<nprocxy) THEN
      
      ys1(:) = a0(:,2)
      ys2(:) = a0(:,je-1)

      ! Get processor numbers of neighbors
      ! For the sake of simplicity  the p_sendrecv call is as
      ! for periodic boundary conditions
      nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
      np = MOD(p_pe+nprocxy+nprocx,nprocxy)

      CALL p_sendrecv(ys1,nm,yr1,np,1)
      CALL p_sendrecv(ys2,np,yr2,nm,2)

      IF(.NOT. have_g_js) THEN
        a0(:, 1) = yr2(:)
      ENDIF
      IF(.NOT. have_g_je) THEN
        a0(:,je) = yr1(:)
      ENDIF
      
    ENDIF

#ifdef bounds_exch_tp
                                                           ! do the top margin first 
    IF(have_g_js) THEN

    ne=(nprocx-1)-p_pe                                      ! find the neighbour cpu 
       
       ysl3(:) = a0(:,3)                                    ! define send/recive buffer
       ysl2(:) = a0(:,2)
       ysl4(:) = a0(:,4)


       if ( p_pe /= ne ) then                      
          CALL p_sendrecv(ysl3,ne,yrl3,ne,3)                ! exchange with neighbour cpu
          CALL p_sendrecv(ysl2,ne,yrl2,ne,4)                ! exchange with neighbour cpu
          CALL p_sendrecv(ysl4,ne,yrl4,ne,5)                ! exchange with neighbour cpu
       else                                                 ! dont do send/recive on the same cpu
          yrl3(:)=ysl3(:)
          yrl2(:)=ysl2(:)
          yrl4(:)=ysl4(:)
       endif

       IF(ttt=='p') THEN                                    ! p point without sign change 
!          DO i=2,ie-1                              
          DO i=1,ie
             il=i                             
             ir=ie+1-i 
             a0(il,2) = yrl3(ir)                            ! syncronise line 2 with line 3
             a0(il,1) = yrl4(ir)                            ! syncronise line 1 with line 4
          END DO

       ENDIF

       IF (ttt=='p-') THEN                                  ! p point  with sign change 
!          DO i=2,ie-1                              
          DO i=1,ie                              
             il=i                             
             ir=ie+1-i 
             a0(il,2) = -yrl3(ir)                            ! syncronise line 2 with line 3
             a0(il,1) = -yrl4(ir)                            ! syncronise line 1 with line 4
          END DO
       ENDIF

       IF (ttt=='v') THEN                                   ! v point with sign change
!          DO i=2,ie-1
          DO i=1,ie
             il=i
             ir=ie+1-i     
             a0(il,1) = -yrl3(ir)                           ! syncronise line 1 with line 3
          END DO
       ENDIF


       IF (ttt=='vf') THEN                                   ! v point with sign change
!          DO i=2,ie-1
          DO i=1,ie
             il=i
             ir=ie+1-i     
             a0(il,1) = -yrl3(ir)                           ! syncronise line 1 with line 3
          END DO

!          DO i=2,ie/2-1
          DO i=1,ie/2
             il=i
             ir=ie+1-i     
             pq=0.5*(a0(il,2)+yrl2(ir))
             a0(il,2) =  pq
             a0(ir,2) = -pq
          enddo
       ENDIF

       IF (ttt=='v+') THEN                                  ! v point without sign change
!          DO i=2,ie-1
          DO i=1,ie
             il=i
             ir=ie+1-i     
             a0(il,1) = yrl3(ir)                            ! syncronise line 1 with line 3
          END DO

!          DO i=2,ie-1
          DO i=1,ie
             il=i
             ir=ie+1-i     
             pq=0.5*(a0(il,2)+yrl2(ir))
             a0(il,2) = pq                                  ! syncronise line 2 with line 2
          END DO


       ENDIF

       IF (ttt=='vv') THEN  
          ! nothing to do
       ENDIF

       IF (ttt=='vd') THEN                      ! v point with sign change
!          DO i=2,ie-1
          DO i=1,ie
             il=i
             ir=ie+1-i     
             a0(il,1)=-yrl2(ir)
          END DO
       ENDIF

       IF (ttt=='u') THEN                                   ! u point with sign change
!          DO i=2,ie-1
          DO i=1,ie-1
             il=i
             ir=ie-i     
             a0(il,2) = -yrl3(ir)                          ! syncronise line 2 with line 3 
          END DO
       ENDIF

       IF (ttt=='u+') THEN                                  ! u point without sign change
!          DO i=2,ie-1
          DO i=1,ie-1
             il=i
             ir=ie-i     
             a0(il,2) = yrl3(ir)                           ! syncronise line 2 with line 3 
          END DO
       ENDIF

       IF (ttt=='uu') THEN  
          ! nothing to do
       ENDIF

       IF (ttt=='s') THEN                                   ! psi point without sign change
!           DO i=2,ie-1
           DO i=1,ie-1
             il=i
             ir=ie-i     
             a0(il,1) = yrl2(ir)          ! syncronise line 1 with line 2
          END DO
       ENDIF

       IF (ttt=='s-') THEN                                  ! psi point with sign change
!           DO i=2,ie-1
           DO i=1,ie-1
             il=i
             ir=ie-i     
             a0(il,1) = -yrl2(ir)          ! syncronise line 1 with line 2
          END DO
       ENDIF
    ENDIF
#endif

!    te = p_time()
!    t2d = t2d + te-ts
!    n2d = n2d + 1

#ifdef bounds_exch_check

    aa1(:,:) = a0(:,:)
    nn=0
    DO jj=1,je
       DO ii=1,ie
          IF (aa1(ii,jj).ne.aa0(ii,jj)) THEN
            nn=nn+1
          ENDIF
       ENDDO
    ENDDO
    call global_sum(nn)
    if ( p_pe .eq. p_io ) then
       if ( nn .eq. 0 )   WRITE(0,*) 'attn: pe ',p_pe,' of ',nprocxy,' needless bounds_exch! ',text
    endif 
    DEALLOCATE(aa0,aa1)

#endif

    
    IF(p_nprocs > nprocxy) THEN
      ! Test mode
      IF(PRESENT(text)) THEN
        CALL para_check_2d(a0,text)
      ELSE
        CALL para_check_2d(a0,'bounds_exch_2d')
      ENDIF
    ENDIF
        
  END SUBROUTINE bounds_exch_hh_2d
                     
  SUBROUTINE bounds_exch_3d(arr,text)
        
    ! Exchanges boundaries of 3D arrays
        
#ifdef bounds_exch_put
    USE mo_boundsexch
        
    REAL xs1(je,ke+1), xs2(je,ke+1)
    REAL ys1(ie,ke+1), ys2(ie,ke+1)
#endif
    REAL, INTENT(INOUT) :: arr(:,:,:)
    CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
        
    INTEGER nm, np, kk
#ifndef bounds_exch_put
    REAL x3r1(je,UBOUND(arr,3)), x3r2(je,UBOUND(arr,3))
    REAL y3r1(ie,UBOUND(arr,3)), y3r2(ie,UBOUND(arr,3))
    REAL xs1(je,UBOUND(arr,3)), xs2(je,UBOUND(arr,3))
    REAL ys1(ie,UBOUND(arr,3)), ys2(ie,UBOUND(arr,3))
#endif

    INTEGER :: i, j, k
        
    ts = p_time()

    kk = UBOUND(arr,3)
        
    ! x-direction
        
    DO k = 1, kk
      DO j = 1, je 
        xs1(j,k) = arr(2,j,k)
        xs2(j,k) = arr(ie-1,j,k)
      ENDDO
    ENDDO

    IF (nprocx>1 .AND. p_pe<nprocxy) THEN
      ! Get processor numbers of neighbors
      nm = MERGE(p_pe-1,p_pe+nprocx-1,p_num_x(p_pe)/=0)
      np = MERGE(p_pe+1,p_pe-nprocx+1,p_num_x(p_pe)/=nprocx-1)
#ifdef bounds_exch_isend
      count = size(xs1)
      CALL p_isend(xs1(1,1),nm,1,p_count = count)
      CALL p_isend(xs2(1,1),np,2,p_count = count)
      CALL p_irecv(x3r1(1,1),np,1,p_count = count)
      CALL p_irecv(x3r2(1,1),nm,2,p_count = count)
      CALL p_wait
#else
#ifdef bounds_exch_put
      CALL p_win_fence(p_win3dx)
      CALL p_put(xs1,nm,p_win3dx(1))
      CALL p_put(xs2,np,p_win3dx(2))
      CALL p_win_fence(p_win3dx)
#else
      CALL p_sendrecv(xs1,nm,x3r1,np,1)
      CALL p_sendrecv(xs2,np,x3r2,nm,2)
#endif
#endif
    ELSE
      x3r1(:,:) = xs1(:,:)
      x3r2(:,:) = xs2(:,:)
    ENDIF
      
    IF(icycli/=0 .OR. p_num_x(p_pe)/=0)        arr( 1,:,:) = x3r2(:,1:kk)
    IF(icycli/=0 .OR. p_num_x(p_pe)/=nprocx-1) arr(ie,:,:) = x3r1(:,1:kk)
        
    ! y-direction
    ! Note that there is no action required if nprocy==1
        
    IF (nprocy>1 .AND. p_pe<nprocxy) THEN
      DO k = 1, kk
        DO i = 1, ie 
          ys1(i,k) = arr(i,2,k)
          ys2(i,k) = arr(i,je-1,k)
        ENDDO
      ENDDO
      ! Get processor numbers of neighbors
      ! For the sake of simplicity  the p_sendrecv call is as
      ! for periodic boundary conditions
      nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
      np = MOD(p_pe+nprocxy+nprocx,nprocxy)
#ifdef bounds_exch_isend
      count = size(ys1)
      CALL p_isend(ys1(1,1),nm,1,p_count = count)
      CALL p_isend(ys2(1,1),np,2,p_count = count)
      CALL p_irecv(y3r1(1,1),np,1,p_count = count)
      CALL p_irecv(y3r2(1,1),nm,2,p_count = count)
      CALL p_wait
#else
#ifdef bounds_exch_put
      CALL p_win_fence(p_win3dy)
      CALL p_put(ys1,nm,p_win3dy(1))
      CALL p_put(ys2,np,p_win3dy(2))
      CALL p_win_fence(p_win3dy)
#else
      CALL p_sendrecv(ys1,nm,y3r1,np,1)
      CALL p_sendrecv(ys2,np,y3r2,nm,2)
#endif
#endif
      IF(p_num_y(p_pe)/=0)        arr(:, 1,:) = y3r2(:,1:kk)
      IF(p_num_y(p_pe)/=nprocy-1) arr(:,je,:) = y3r1(:,1:kk)
    ENDIF
        
!    te = p_time()
!    t3d = t3d + te-ts
!    n3d = n3d + 1
        
    IF(p_nprocs > nprocxy) THEN
      ! Test mode
      IF(PRESENT(text)) THEN
        CALL para_check_3d(arr,text)
      ELSE
        CALL para_check_3d(arr,'bounds_exch_3d')
      ENDIF
    ENDIF
        
  END SUBROUTINE bounds_exch_3d
      
  !-----------------------------------------------------------------------
    SUBROUTINE bounds_exch_hh_3d(ttt,arr,text)
        
    ! Exchanges boundaries of 3D arrays

#ifdef bounds_exch_check
    use mo_param1, only : ie,je
    INTEGER NN,jj,ii
    REAL,ALLOCATABLE :: aa0(:,:,:),aa1(:,:,:)
#endif
        
    REAL, INTENT(INOUT) :: arr(:,:,:)
    CHARACTER (LEN=*), INTENT(IN) :: ttt
    CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
    INTEGER nm, np, kk 

    REAL x3r1(je,UBOUND(arr,3)), x3r2(je,UBOUND(arr,3))
    REAL y3r1(ie,UBOUND(arr,3)), y3r2(ie,UBOUND(arr,3))
    REAL xs1(je,UBOUND(arr,3)), xs2(je,UBOUND(arr,3))
    REAL ys1(ie,UBOUND(arr,3)), ys2(ie,UBOUND(arr,3))

    INTEGER :: i, j, k

#ifdef bounds_exch_tp
    REAL ysl2(ie,ke+1), ysl3(ie,ke+1), yrl2(ie,ke+1), yrl3(ie,ke+1) ,ysl4(ie,ke+1),yrl4(ie,ke+1)
    REAL pq(ke+1)
    INTEGER ne, il, ir
#endif


        
    ts = p_time()

    kk = UBOUND(arr,3)
        
#ifdef bounds_exch_check
    ALLOCATE(aa0(ie,je,kk),aa1(ie,je,kk))
    do k=1,kk
       aa0(:,:,k) = arr(:,:,k)
    enddo
#endif

    ! x-direction
        
    DO k = 1, kk
      DO j = 1, je 
        xs1(j,k) = arr(2,j,k)
        xs2(j,k) = arr(ie-1,j,k)
      ENDDO
    ENDDO

    IF (nprocx>1 .AND. p_pe<nprocxy) THEN
      ! Get processor numbers of neighbors
      nm = MERGE(p_pe-1,p_pe+nprocx-1,p_num_x(p_pe)/=0)
      np = MERGE(p_pe+1,p_pe-nprocx+1,p_num_x(p_pe)/=nprocx-1)

      CALL p_sendrecv(xs1,nm,x3r1,np,1)
      CALL p_sendrecv(xs2,np,x3r2,nm,2)
    ELSE
      x3r1(:,:) = xs1(:,:)
      x3r2(:,:) = xs2(:,:)
    ENDIF
      
    IF(icycli/=0 .OR. p_num_x(p_pe)/=0)        arr( 1,:,:) = x3r2(:,1:kk)
    IF(icycli/=0 .OR. p_num_x(p_pe)/=nprocx-1) arr(ie,:,:) = x3r1(:,1:kk)
        
    ! y-direction
    ! Note that there is no action required if nprocy==1
        
    IF (nprocy>1 .AND. p_pe<nprocxy) THEN
      DO k = 1, kk
        DO i = 1, ie 
          ys1(i,k) = arr(i,2,k)
          ys2(i,k) = arr(i,je-1,k)
        ENDDO
      ENDDO
      ! Get processor numbers of neighbors
      ! For the sake of simplicity  the p_sendrecv call is as
      ! for periodic boundary conditions
      nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
      np = MOD(p_pe+nprocxy+nprocx,nprocxy)

      CALL p_sendrecv(ys1,nm,y3r1,np,1)
      CALL p_sendrecv(ys2,np,y3r2,nm,2)

      IF(p_num_y(p_pe)/=0)        arr(:, 1,:) = y3r2(:,1:kk)
      IF(p_num_y(p_pe)/=nprocy-1) arr(:,je,:) = y3r1(:,1:kk)
    ENDIF

#ifdef bounds_exch_tp


    IF(have_g_js) THEN

       ne=(nprocx-1)-p_pe                    ! find the neighbour cpu 

       ysl3(:,1:kk) = arr(:,3,:)          ! define send/recive buffer
       ysl2(:,1:kk) = arr(:,2,:)
       ysl4(:,1:kk) = arr(:,4,:)


       IF(p_pe /= ne) then 
          CALL p_sendrecv(ysl4,ne,yrl4,ne,22)   ! exchange with neighbour cpu
          CALL p_sendrecv(ysl3,ne,yrl3,ne,33)   ! exchange with neighbour cpu
          CALL p_sendrecv(ysl2,ne,yrl2,ne,44)   ! exchange with neighbour cpu
       else                                     ! dont do send/recive on the same cpu
          yrl3(:,:)=ysl3(:,:)
          yrl2(:,:)=ysl2(:,:)
          yrl4(:,:)=ysl4(:,:)
       endif

        IF(ttt=='p') THEN                                    ! p point without sign change 
!          DO i=2,ie-1                              
          DO i=1,ie
             il=i                             
             ir=ie+1-i 
             arr(il,2,1:kk) = yrl3(ir,1:kk)                            ! syncronise line 2 with line 3
             arr(il,1,1:kk) = yrl4(ir,1:kk)                            ! syncronise line 1 with line 4
          END DO

       ENDIF

       IF (ttt=='p-') THEN                                  ! p point  with sign change 
!          DO i=2,ie-1                              
          DO i=1,ie
             il=i                             
             ir=ie+1-i 
             arr(il,2,1:kk) = -yrl3(ir,1:kk)                            ! syncronise line 2 with line 3
             arr(il,1,1:kk) = -yrl4(ir,1:kk)                            ! syncronise line 1 with line 4
          END DO
       ENDIF

       IF (ttt=='v') THEN                                   ! v point with sign change
!          DO i=2,ie-1
          DO i=1,ie
             il=i
             ir=ie+1-i     
             arr(il,1,1:kk) = -yrl3(ir,1:kk)                           ! syncronise line 1 with line 3
          END DO
       ENDIF

       IF (ttt=='v+') THEN                                  ! v point without sign change
!          DO i=2,ie-1
          DO i=1,ie
             il=i
             ir=ie+1-i     
             arr(il,1,1:kk) = yrl3(ir,1:kk)                            ! syncronise line 1 with line 3
          END DO

!          DO i=2,ie-1
          DO i=1,ie
             il=i
             ir=ie+1-i     
             pq(1:kk)=0.5*(arr(il,2,1:kk)+yrl2(ir,1:kk))
             arr(il,2,1:kk) = pq(1:kk)                                  ! syncronise line 2 with line 2
          END DO

       ENDIF


       IF (ttt=='vf') THEN                                   ! v point with sign change
!          DO i=2,ie-1
          DO i=1,ie
             il=i
             ir=ie+1-i     
             arr(il,1,1:kk) = -yrl3(ir,1:kk)                           ! syncronise line 1 with line 3
          END DO

!          DO i=2,ie/2-1
          DO i=1,ie/2
             il=i
             ir=ie+1-i     
             pq=0.5*(arr(il,2,1:kk)+yrl2(ir,1:kk))
             arr(il,2,1:kk) =  pq
             arr(ir,2,1:kk) = -pq
          enddo
       ENDIF



       IF (ttt=='vd') THEN                      ! v point with sign change
!          DO i=2,ie-1
          DO i=1,ie
             il=i
             ir=ie+1-i     
             arr(il,1,1:kk)=-yrl2(ir,1:kk)
          END DO
       ENDIF


       IF (ttt=='u') THEN                                   ! u point with sign change
!          DO i=2,ie-1
          DO i=1,ie-1
             il=i
             ir=ie-i     
             arr(il,2,1:kk) = -yrl3(ir,1:kk)                          ! syncronise line 2 with line 3 
          END DO
       ENDIF

       IF (ttt=='u+') THEN                                  ! u point without sign change
!          DO i=2,ie-1
          DO i=1,ie-1
             il=i
             ir=ie-i     
             arr(il,2,1:kk) = yrl3(ir,1:kk)                           ! syncronise line 2 with line 3 
          END DO
       ENDIF

       IF (ttt=='s') THEN                                   ! psi point without sign change
!           DO i=2,ie-1
           DO i=1,ie-1
             il=i
             ir=ie-i     
             arr(il,1,1:kk) = yrl2(ir,1:kk)          ! syncronise line 1 with line 2
          END DO
       ENDIF

       IF (ttt=='s-') THEN                                  ! psi point with sign change
!           DO i=2,ie-1
           DO i=1,ie-1
             il=i
             ir=ie-i     
             arr(il,1,1:kk) = -yrl2(ir,1:kk)          ! syncronise line 1 with line 2
          END DO
       ENDIF
    ENDIF
#endif

#ifdef bounds_exch_check
    do k=1,kk
       aa1(:,:,k) = arr(:,:,k)
    enddo
    nn=0
    DO k=1,kk
       DO jj=1,je
          DO ii=1,ie
            IF (aa1(ii,jj,k).ne.aa0(ii,jj,k)) THEN
               nn=nn+1
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   call global_sum(nn)
   if ( p_pe .eq. p_io ) then
   if ( nn .eq. 0 )   WRITE(0,*) 'attn: pe ',p_pe,' of ',nprocxy,' needless bounds_exch! ',text 
   endif
   DEALLOCATE(aa0,aa1)



#endif


    te = p_time()
    t3d = t3d + te-ts
    n3d = n3d + 1
        
    IF(p_nprocs > nprocxy) THEN
      ! Test mode
      IF(PRESENT(text)) THEN
        CALL para_check_3d(arr,text)
      ELSE
        CALL para_check_3d(arr,'bounds_exch_3d')
      ENDIF
    ENDIF
        

  END SUBROUTINE bounds_exch_hh_3d
      
  !-----------------------------------------------------------------------
  
  SUBROUTINE read_slice(iunit,arr)
    
    ! Reads a 2D array and scatters it to all
    
    INTEGER,INTENT(IN) :: iunit
    REAL, INTENT(OUT) :: arr(:,:)
    
    REAL arr_g(ie_g,je_g)
    
    IF(p_pe==p_io) READ(iunit) arr_g
    CALL p_bcast(arr_g,p_io)
    
    arr(:,:) = arr_g(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)
    
  END SUBROUTINE read_slice
      
  !-----------------------------------------------------------------------
      
  SUBROUTINE write_slice(iunit,arr)
        
    ! Gathers a 2D array and writes it to iunit
    
    INTEGER,INTENT(IN) :: iunit
    REAL, INTENT(IN) :: arr(:,:)
    
    REAL arr_g(ie_g,je_g)
    
    CALL gather_arr(arr,arr_g,p_io)
    IF(p_pe==p_io) WRITE(iunit) arr_g
    
  END SUBROUTINE write_slice

  SUBROUTINE write_slice_sp(iunit,arr)
        
    ! Gathers a 2D array and writes it to iunit
    use mo_kind

    INTEGER,INTENT(IN) :: iunit
    REAL, INTENT(IN) :: arr(:,:)
    
    REAL arr_g(ie_g,je_g)
    
    CALL gather_arr(arr,arr_g,p_io)
    IF(p_pe==p_io) WRITE(iunit) real(arr_g,sp)
    
  END SUBROUTINE write_slice_sp
          
  !-----------------------------------------------------------------------
      
  SUBROUTINE global_sum_r(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9)
    
    ! Build global sum for real scalars
    ! For performance reasons we permit up to 10 arguments in a single call
    
    REAL, INTENT(INOUT):: s0
    REAL, INTENT(INOUT), OPTIONAL:: s1,s2,s3,s4,s5,s6,s7,s8,s9
    
    REAL s(10)
    INTEGER n
    
    s(1) = s0
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s(n) = s1; ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s(n) = s2; ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s(n) = s3; ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s(n) = s4; ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s(n) = s5; ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s(n) = s6; ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s(n) = s7; ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s(n) = s8; ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s(n) = s9; ENDIF
         
    CALL global_sum_1d(s(1:n))
      
    s0 = s(1)
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s1 = s(n); ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s2 = s(n); ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s3 = s(n); ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s4 = s(n); ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s5 = s(n); ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s6 = s(n); ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s7 = s(n); ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s8 = s(n); ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s9 = s(n); ENDIF
          
  END SUBROUTINE global_sum_r
        
  !-----------------------------------------------------------------------
      
  SUBROUTINE global_sum_i(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9)
        
    ! Build global sum for integer scalars
    ! For performance reasons we permit up to 10 arguments in a single call
    
    INTEGER, INTENT(INOUT):: s0
    INTEGER, INTENT(INOUT), OPTIONAL:: s1,s2,s3,s4,s5,s6,s7,s8,s9
    
    REAL s(10)
    INTEGER n
    
    s(1) = s0
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s(n) = s1; ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s(n) = s2; ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s(n) = s3; ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s(n) = s4; ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s(n) = s5; ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s(n) = s6; ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s(n) = s7; ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s(n) = s8; ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s(n) = s9; ENDIF
          
    CALL global_sum_1d(s(1:n))
          
    s0 = NINT(s(1))
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s1 = NINT(s(n)); ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s2 = NINT(s(n)); ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s3 = NINT(s(n)); ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s4 = NINT(s(n)); ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s5 = NINT(s(n)); ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s6 = NINT(s(n)); ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s7 = NINT(s(n)); ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s8 = NINT(s(n)); ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s9 = NINT(s(n)); ENDIF
          
  END SUBROUTINE global_sum_i
        
  !-----------------------------------------------------------------------
        
  SUBROUTINE global_sum_1d(s)
    
    ! Build global sum for real 1D array
    
    REAL, INTENT(INOUT):: s(:)
    
    REAL r(SIZE(s)), q(SIZE(s)), errmax, err
    INTEGER n, i, ierr
    
    ! If we are running in test mode, use result from last PE and
    ! check if the others are approximatly right, else do real summation
    
    n = SIZE(s)
    
    IF(p_nprocs>nprocxy) THEN
      
      ! Test mode
      IF(p_pe>=nprocxy) THEN
        q(:) = s(:) ! Save s
        s(:) = 0    ! Remove our contribution to global sum
      ENDIF
      
      ! Get sum on working PEs
      r = p_sum(s)
      
      ! Check against saved value on test PE
      IF(p_pe==nprocxy) THEN
        errmax = 0;
        DO i=1,n
          IF(q(i)/=r(i)) THEN
            err = ABS(q(i)-r(i))/MAX(ABS(q(i)),ABS(r(i)))
            errmax = MAX(errmax,err)
          ENDIF
        ENDDO
        WRITE(nerr,*) 'Global Sum Max Err = ',errmax
        IF(errmax>1.e-5) CALL p_abort
        s(:) = q(:) ! Restore s
      ENDIF
      
      ! We use the value from test PE in order to get always
      ! identical results during test
      CALL p_bcast(s,nprocxy)
      
    ELSE
      
      r = p_sum(s)
      s(:) = r(:)
      
    ENDIF
    
  END SUBROUTINE global_sum_1d
        
  !-----------------------------------------------------------------------
        
  SUBROUTINE global_sum_2d(s)
    
    ! Build global sum for real 2D array
    
    REAL, INTENT(INOUT):: s(:,:)
    
    REAL r(SIZE(s))
    INTEGER ierr
    
    IF(p_nprocs>nprocxy) THEN
      ! Test mode - use global_sum_1d
      r = RESHAPE(s,(/SIZE(s)/))
      CALL global_sum_1d(r)
      s = RESHAPE(r,SHAPE(s))
    ELSE
      r = p_sum(RESHAPE(s,(/SIZE(s)/)))
      s = RESHAPE(r,SHAPE(s))
    ENDIF
    
  END SUBROUTINE global_sum_2d
  !-----------------------------------------------------------------------
        
  SUBROUTINE global_mean_2d(arr,su)

    USE MO_KIND 
    USE MO_PARAM1, only : ie,je,ie_g,je_g
    USE MO_COMMO1, only : dlxp, dlyp, weto    

    ! Build global mean for real 2D array

    INTEGER             ::  i,j
    REAL, INTENT(IN)    ::  arr(:,:)
    REAL, INTENT(INOUT) ::  su

    REAL                ::  as
    REAL                ::  arr2(ie,je), sarr(ie,je)


    REAL, POINTER       ::  arr_g(:,:)
    REAL, POINTER       ::  sarr_g(:,:)

    if (p_pe==p_io) then
       ALLOCATE(arr_g(ie_g,je_g),sarr_g(ie_g,je_g))
    else
       arr_g => NULL()
       sarr_g => NULL()
    endif

    do i=1,ie
       do j=1,je 
          arr2(i,j)=arr(i,j)*dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
          sarr(i,j)=dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
       enddo
    enddo

    call gather(arr2,arr_g,p_io)
    call gather(sarr,sarr_g,p_io)

    if ( p_pe == p_io ) then
       su=0.
       as=0. 
       do i=2,ie_g-1
          do j=2,je_g-1
             su=su+arr_g(i,j)
             as=as+sarr_g(i,j)
          enddo
       enddo
       if (as==0.) then
         su=0.	
       else	
         su=su/as	
       endif
    endif

    call p_bcast(su,p_io)

  END SUBROUTINE global_mean_2d        
  !-----------------------------------------------------------------------
        
  SUBROUTINE global_sum_2d_pio(arr,su)

    USE MO_KIND 
    USE MO_PARAM1, only : ie_g,je_g

    INTEGER             ::  i,j,jb
    REAL, INTENT(IN)    ::  arr(:,:)
    REAL, INTENT(INOUT) ::  su

    REAL, ALLOCATABLE       ::  arr_g(:,:)

    if (p_pe==p_io) then
       ALLOCATE(arr_g(ie_g,je_g))
    else
       ALLOCATE(arr_g(0,0))
    endif

    call gather_arr(arr,arr_g,p_io)

    jb=2
#ifdef bounds_exch_tp    
    jb=3
#endif
    if ( p_pe == p_io ) then
       su=0.
       do i=2,ie_g-1
          do j=jb,je_g-1
             su=su+arr_g(i,j)
          enddo
       enddo
    endif
    call p_bcast(su,p_io)

    DEALLOCATE(arr_g)

  END SUBROUTINE global_sum_2d_pio        
  !-----------------------------------------------------------------------
        
  SUBROUTINE global_max(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9)
          
    ! Build global max for real scalars
    ! For performance reasons we permit up to 10 arguments in a single call
    
    REAL, INTENT(INOUT):: s0
    REAL, INTENT(INOUT), OPTIONAL:: s1,s2,s3,s4,s5,s6,s7,s8,s9
    
    REAL s(10), r(10)
    INTEGER n, ierr
    s(:)=0.
    r(:)=0.
    s(1) = s0
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s(n) = s1; ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s(n) = s2; ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s(n) = s3; ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s(n) = s4; ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s(n) = s5; ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s(n) = s6; ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s(n) = s7; ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s(n) = s8; ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s(n) = s9; ENDIF
      
    r = p_max(s)
            
    s0 = r(1)
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s1 = r(n); ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s2 = r(n); ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s3 = r(n); ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s4 = r(n); ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s5 = r(n); ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s6 = r(n); ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s7 = r(n); ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s8 = r(n); ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s9 = r(n); ENDIF
     
  END SUBROUTINE global_max
                            
  !-----------------------------------------------------------------------
                            
  SUBROUTINE global_min(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9)
    
    ! Build global min for real scalars
    ! For performance reasons we permit up to 10 arguments in a single call
    
    REAL, INTENT(INOUT):: s0
    REAL, INTENT(INOUT), OPTIONAL:: s1,s2,s3,s4,s5,s6,s7,s8,s9
    
    REAL s(10), r(10)
    INTEGER n, ierr
    s(:)=0.
    r(:)=0.  
    
    s(1) = s0
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s(n) = s1; ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s(n) = s2; ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s(n) = s3; ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s(n) = s4; ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s(n) = s5; ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s(n) = s6; ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s(n) = s7; ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s(n) = s8; ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s(n) = s9; ENDIF
     
    r = p_min(s)
    
    s0 = r(1)
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s1 = r(n); ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s2 = r(n); ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s3 = r(n); ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s4 = r(n); ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s5 = r(n); ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s6 = r(n); ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s7 = r(n); ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s8 = r(n); ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s9 = r(n); ENDIF
    
  END SUBROUTINE global_min
  
  !-----------------------------------------------------------------------
  
  REAL FUNCTION global_array_sum(arr)
    
    ! Builds the global sum of a 2D array by gathering it on one PE,
    ! calculating th sum on this PE and broadcasting the result
    ! As opposed to calculating the local sum and the using a global_sum call,
    ! this method should give identical results indepentently of the number
    ! of processors, at least when optimization is switched off.
    
    REAL, INTENT(IN) :: arr(:,:)
    
    INTEGER i,j
    REAL arr_g(ie_g,je_g), s, sum
    
    CALL gather_arr(arr,arr_g,0)
    
    IF(p_pe==0) THEN
      sum = 0
      DO j=2,je_g-1
        DO i=2,ie_g-1
          sum = sum + arr_g(i,j)
        ENDDO
      ENDDO
    ENDIF
    
    CALL p_bcast(sum,0)
    
    global_array_sum = sum
    
    IF(p_pe==nprocxy) THEN
      s = 0
      DO j=2,je-1
        DO i=2,ie-1
          s = s + arr(i,j)
        ENDDO
      ENDDO
      
      IF(s/=sum) THEN
        WRITE(nerr,*) 'global_array_sum: ',s,sum,ABS(s-sum)
        CALL p_abort
      ENDIF
    ENDIF
    
  END FUNCTION global_array_sum
  
  !-----------------------------------------------------------------------
  
  SUBROUTINE para_check_2d(arr,text,lev)
    
    ! If running in test mode, checks if the results on the PEs running
    ! in parallel are identical with the results on the test PE
    ! This works only if optimization is switched off
    
    REAL, INTENT(IN) :: arr(:,:)
    CHARACTER (LEN=*), INTENT(IN) :: text
    INTEGER, INTENT(IN), OPTIONAL :: lev
    
    REAL, ALLOCATABLE :: aux(:,:)
    INTEGER iie, jje, ioff, joff, i, j, n, num
    
    IF(p_nprocs<=nprocxy) RETURN
    
    
    IF(p_pe<nprocxy) THEN
      CALL p_send(arr,nprocxy,1111)
    ELSE
      num = 0
      DO n=0,nprocxy-1
        iie = p_size_x(n)
        jje = p_size_y(n)
        ioff = p_ioff_g(n)
        joff = p_joff_g(n)
        ALLOCATE(aux(iie,jje))
        CALL p_recv(aux,n,1111)
        DO j=1,jje
          DO i=1,iie
            IF(aux(i,j) /= arr(i+ioff,j+joff)) THEN
              IF(num==0) THEN
                IF(PRESENT(lev)) THEN
                  WRITE(nerr,*) 'Consistency Check Error! K=',lev,text
                ELSE
                  WRITE(nerr,*) 'Consistency Check Error! ',text
                ENDIF
              ENDIF
              WRITE(nerr,'(3i4,2e25.16)') n,i,j,aux(i,j),arr(i+ioff,j+joff)
              num = num+1
            ENDIF
          ENDDO
        ENDDO
        DEALLOCATE(aux)
      ENDDO
      IF(num>0) THEN
        WRITE(nerr,*) num,' errors'
        CALL p_abort
      ENDIF
    ENDIF
    
    CALL p_barrier(p_all_comm)
    
  END SUBROUTINE para_check_2d
  
  !-----------------------------------------------------------------------
  
  SUBROUTINE para_check_3d(arr,text)
    
    ! interface for 3D arrays
    
    REAL, INTENT(IN) :: arr(:,:,:)
    CHARACTER (LEN=*), INTENT(IN) :: text
    
    INTEGER :: k
    
    DO k=1,UBOUND(arr,3)
      CALL para_check_2d(arr(:,:,k),text,k)
    ENDDO
    
  END SUBROUTINE para_check_3d
  
  !-----------------------------------------------------------------------
  
  SUBROUTINE stop_all(text)
    
    ! Does an emergency stop by calling p_abort
    
    CHARACTER (LEN=*), INTENT(IN) :: text
    
    WRITE(nerr,*) text
    
    CALL p_abort
    
  END SUBROUTINE stop_all
  
  !-----------------------------------------------------------------------
  
  SUBROUTINE print_stats
    
    REAL :: t2, t3
    INTEGER :: n
    
    IF(p_pe==p_io) THEN
      WRITE(nerr,*) 'Number of 2D boundary exchanges: ',n2d
      WRITE(nerr,*) 'Number of 3D boundary exchanges: ',n3d
    ENDIF
    
    DO n=0,nprocxy-1
      t2 = t2d
      t3 = t3d
      CALL p_bcast(t2,n)
      CALL p_bcast(t3,n)
      IF(p_pe==p_io) THEN
        WRITE(nerr,'(a,i4,a,2f10.3)') 'PE: ',n, &
             ' Times for 2D/3D boundary exchanges: ',t2,t3
      ENDIF
    ENDDO
    
  END SUBROUTINE print_stats
  
#ifdef bounds_exch_put
  SUBROUTINE set_put_window
    
    USE mo_mpi
    USE mo_boundsexch
    
    CALL alloc_mem_boundsexch
    
    CALL p_win_create(xr1,p_win2dx(1))
    CALL p_win_create(xr2,p_win2dx(2))
    CALL p_win_fence(p_win2dx)
    CALL p_win_create(yr1,p_win2dy(1))
    CALL p_win_create(yr2,p_win2dy(2))
    CALL p_win_fence(p_win2dy)
    CALL p_win_create(x3r1,p_win3dx(1))
    CALL p_win_create(x3r2,p_win3dx(2))
    CALL p_win_fence(p_win3dx)
    CALL p_win_create(y3r1,p_win3dy(1))
    CALL p_win_create(y3r2,p_win3dy(2))
    CALL p_win_fence(p_win3dy)
    
  END SUBROUTINE set_put_window
  
  SUBROUTINE close_put_window

    USE mo_mpi
    USE mo_boundsexch
    
    CALL p_win_free(p_win2dx)
    CALL p_win_free(p_win2dy)
    CALL p_win_free(p_win3dx)
    CALL p_win_free(p_win3dy)
    
  END SUBROUTINE close_put_window
#endif /*bounds_exch_put*/
  
  
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE sethaloN(a0,a1,ihalo)

! ao(ie,je), a1(ie+2,je+2) 


    REAL, INTENT(INOUT) :: a0(:,:),a1(:,:)
    INTEGER nm, np, ihalo,i,j,ii,jj,iii,jjj,nul,nll,nur,nlr 
    REAL xr1(je,ihalo-1),xr2(je,ihalo-1),yr1(ie,ihalo-1),yr2(ie,ihalo-1)
    REAL xs1(je,ihalo-1),xs2(je,ihalo-1),ys1(ie,ihalo-1),ys2(ie,ihalo-1)
    REAL xsul(ihalo-1,ihalo-1),xsur(ihalo-1,ihalo-1)    &
        ,xsll(ihalo-1,ihalo-1),xslr(ihalo-1,ihalo-1)
    REAL xrul(ihalo-1,ihalo-1),xrur(ihalo-1,ihalo-1)    &
        ,xrll(ihalo-1,ihalo-1),xrlr(ihalo-1,ihalo-1)

    a1(:,:)=0. 

    do i=1,ie
       do j=1,je
          a1(i+ihalo-1,j+ihalo-1)=a0(i,j)
       enddo
    enddo

    ! x-direction  

!    write(0,*)'xdir'

    do ii=1,ihalo-1   
       xs1(:,ii) = a0(ii+2,:)
       xs2(:,ii) = a0(ie-(ii+1),:)
    enddo
    

!    do ii=1,ihalo-1
!     do j=1,je
!       write(0,*)"xs1(",j,":",ii,") =", xs1(j,ii)
!     enddo
!    enddo
!    do ii=1,ihalo-1
!     do j=1,je
!       write(0,*)"xs2(",j,":",ii,") =", xs2(j,ii)
!     enddo
!    enddo


    IF (nprocx>1 .AND. p_pe<nprocxy ) THEN
      ! Get processor numbers of neighbors
      nm = MERGE(p_pe-1,p_pe+nprocx-1,p_num_x(p_pe)/=0)
      np = MERGE(p_pe+1,p_pe-nprocx+1,p_num_x(p_pe)/=nprocx-1)


      CALL p_sendrecv(xs1,nm,xr1,np,201)
      CALL p_sendrecv(xs2,np,xr2,nm,202)






    ELSE
      xr1(:,:)=xs1(:,:)
      xr2(:,:)=xs2(:,:)

! ecken kopieren

    ENDIF
    
!      do ii=1,ihalo-1
!     do j=1,je
!       write(0,*)"xr1(",j,":",ii,") =", xr1(j,ii)
!     enddo
!    enddo
!    do ii=1,ihalo-1
!     do j=1,je
!       write(0,*)"xr2(",j,":",ii,") =", xr2(j,ii)
!     enddo
!    enddo

!------ update left halo -----------------

    IF(icycli/=0 .OR. .NOT. have_g_is) THEN
      do ii=1,ihalo-1
       do j=1,je
          a1(ii,j+1) = xr2(j,ihalo-ii)
       enddo
      enddo
!
    ENDIF

!------- update right halo -----------------

    IF(icycli/=0 .OR. .NOT. have_g_ie) THEN

       do ii=1,ihalo-1
        do j=1,je
         a1(ie+(ii+ihalo-1),j+1) = xr1(j,ii)
       enddo
      enddo
!
    ENDIF




!    return

    ! y-direction
      
!    write(0,*)'ydir'
      
      do ii=1,ihalo-1  
        ys1(:,ii) = a0(:,1+ii+1)
        ys2(:,ii) = a0(:,je-(ii+1))
      enddo

    IF (nprocy>1 .AND. p_pe<nprocxy) THEN

      ! Get processor numbers of neighbors
      nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
      np = MOD(p_pe+nprocxy+nprocx,nprocxy)

      CALL p_sendrecv(ys1,nm,yr1,np,301)
      CALL p_sendrecv(ys2,np,yr2,nm,302)

   ELSE
      yr1(:,:)=ys1(:,:)
      yr2(:,:)=ys2(:,:)
   ENDIF
!-------- update north halo-----------------

   IF(icycli/=0 .OR. .NOT. have_g_js) THEN
    do ii=1,ihalo-1
       do i=1,ie
          a1(i+1,ii) = yr2(i,ihalo-ii)
       enddo
    enddo
!
   ENDIF

!-------- update south halo-----------------

   IF(icycli/=0 .OR. .NOT. have_g_je) THEN
      do ii=1,ihalo-1 
         do i=1,ie
           a1(i+1,je+(ii+1)) = yr1(i,ii)
         enddo
      enddo
   ENDIF
   


! -------diagonal------------------

   if (ihalo > 2 ) then

      do ii=1,ihalo-1
         do jj=1,ihalo-1
            
            iii=ihalo-ii
            jjj=ihalo-jj
            
            xsul(ii,jj) = a0(ii+1,jj+1)

            xsur(ii,jj) = a0(ie-1-iii,jj+1)
            
            xsll(ii,jj) = a0(ii+1,je-1-jjj)
            
            xslr(ii,jj) = a0(ie-1-iii,je-1-jjj)
            
         enddo
      enddo
      
      if (nprocx>1 .and. nprocy>1) then
         
         nul = MOD(p_pe+nprocxy-nprocx,nprocxy)-1
         if (have_g_is) nul=p_pe-1
         
         nll = MOD(p_pe+nprocxy+nprocx,nprocxy)-1
         if (have_g_is) nll=nprocxy-1 
         
         nur = MOD(p_pe+nprocxy-nprocx,nprocxy)+1
         if (have_g_ie) nur=nur-nprocx-1

         nlr = MOD(p_pe+nprocxy+nprocx,nprocxy)+1
         if (have_g_ie) nur=p_pe+1

         CALL p_sendrecv(xsul,nul,xrul,nul,203)
         CALL p_sendrecv(xsur,nur,xrur,nur,204)
         CALL p_sendrecv(xsll,nll,xrll,nll,205)
         CALL p_sendrecv(xslr,nlr,xrlr,nlr,206)

         

         do ii=1,ihalo-1
            do jj=1,ihalo-1

               iii=ie+(ihalo-1)
               jjj=je+(ihalo-1)
               
               a1(ii,jj) = xrul(ii,jj)            !------ update upper left -----------------
               
               a1(iii+ii,jj) = xrur(ii,jj)        !------- update upper right -----------------
               
               a1(ii,jjj+jj) = xrll(ii,jj)        !------ update lower left -----------------

               a1(iii+ii,jjj+jj) = xrlr(ii,jj)    !------ update lower left -----------------

            enddo
         enddo
         

      endif  ! nprocx>1 .and. nprocy>1

   endif  ! ihalo > 2
  
  END SUBROUTINE sethaloN


  SUBROUTINE sethalo2(a0,a1,ihalo)

! ao(ie,je), a1(ie+2,je+2) 


    REAL, INTENT(INOUT) :: a0(:,:),a1(:,:)
    INTEGER nm, np, ihalo,i,j
    REAL xr1(je),xr2(je),yr1(ie),yr2(ie)
    REAL xs1(je),xs2(je),ys1(ie),ys2(ie)

    a1(:,:)=0. 

    do i=1,ie
       do j=1,je
          a1(i+1,j+1)=a0(i,j)
       enddo
    enddo

    ! x-direction halo=2

!    write(0,*)'xdir'

    xs1(:) = a0(3,:)
    xs2(:) = a0(ie-2,:)

    IF (nprocx>1 .AND. p_pe<nprocxy ) THEN
      ! Get processor numbers of neighbors
      nm = MERGE(p_pe-1,p_pe+nprocx-1,p_num_x(p_pe)/=0)
      np = MERGE(p_pe+1,p_pe-nprocx+1,p_num_x(p_pe)/=nprocx-1)
      
      CALL p_sendrecv(xs1,nm,xr1,np,201)
      CALL p_sendrecv(xs2,np,xr2,nm,202)

    ELSE
      xr1(:)=xs1(:)
      xr2(:)=xs2(:)
    ENDIF
    
    IF(icycli/=0 .OR. .NOT. have_g_is) THEN
      do j=1,je
         a1(1,j+1) = xr2(j)         
      enddo
    ENDIF

    IF(icycli/=0 .OR. .NOT. have_g_ie) THEN
      do j=1,je
         a1(ie+2,j+1) = xr1(j)
      enddo
    ENDIF
      
    ! y-direction
      
!    write(0,*)'ydir'
      
      ys1(:) = a0(:,3)
      ys2(:) = a0(:,je-2)

    IF (nprocy>1 .AND. p_pe<nprocxy) THEN

      ! Get processor numbers of neighbors
      nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
      np = MOD(p_pe+nprocxy+nprocx,nprocxy)

      CALL p_sendrecv(ys1,nm,yr1,np,301)
      CALL p_sendrecv(ys2,np,yr2,nm,302)

   ELSE
      yr1(:)=ys1(:)
      yr2(:)=ys2(:)
   ENDIF

   IF(icycli/=0 .OR. .NOT. have_g_js) THEN
      do i=1,ie
         a1(i+1, 1) = yr2(i)
      enddo
   ENDIF
   IF(icycli/=0 .OR. .NOT. have_g_je) THEN
      do i=1,ie
         a1(i+1,je+2) = yr1(i)
      enddo
   ENDIF
   
  
  END SUBROUTINE sethalo2


END MODULE mo_parallel
