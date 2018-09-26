SUBROUTINE init_decomposition

  USE mo_doctor,        ONLY: nout
  USE mo_control,       ONLY: nproca, nprocb, ngl, nlon, nlev, nm      &
                            , nn, nk, lcolumn, nproma
  USE mo_exception,     ONLY: finish, message
  USE mo_mpi,           ONLY: p_nprocs, p_pe, p_io, p_parallel         &
                            , p_parallel_io
  USE mo_decomposition, ONLY: local_decomposition,                     & 
                              global_decomposition,                    &
                              print_decomposition,                     &
                              decompose, debug_seriell
  USE mo_transpose,     ONLY: indx
  USE mo_column,        ONLY: lat_1d, lon_1d !, inicolumn
  USE mo_advection,     ONLY: iadvec
#ifdef DEBUG
  USE mo_util_string,   ONLY: separator
#endif

  IMPLICIT NONE

  INTEGER :: p,i
  
  LOGICAL :: lrot, lfull_m
  INTEGER :: debug

  lrot          = .TRUE.  ! true: no rotation of longitudes
  lfull_m       = .FALSE. ! true: full m-columns per PE
  debug_seriell = .FALSE. ! true: same results as ser.version if nprocb == 1

  ! debug = 0,1 : PE0 takes full domain (always no rotation)
  !         -1  : no special treatment of PE 0
  !          0  : gather from PE 0
  !          1  : gather from PEs > 0

#ifdef NOMPI
  IF (nproca*nprocb > 1) THEN
     CALL message ('', &
          'Parallel execution selected with wrong compiler options')
     CALL message ('', 'Please, recompile without -DNOMPI') 
     CALL finish ('init_decomposition', 'Program aborted')
  END IF
#endif

  IF (p_nprocs == nproca*nprocb) THEN
     debug = -1
  ELSE IF ( p_nprocs == nproca*nprocb+1) THEN
     debug = 0
  ELSE
     CALL finish ('init_decomposition',                                &
          'Number of runtime PEs doesn''t fit nproca*nprocb(+1)')
  END IF

  IF (p_parallel .AND. p_parallel_io) THEN 
     WRITE (nout,'(/,a,i4,a,i3,a,i3)')                                 &
          ' Total number of PEs: ', p_nprocs,                          &
          ' set A: ', nproca, ' set B: ', nprocb
  ENDIF

  ALLOCATE (global_decomposition(1:p_nprocs))

  ! derive decomposition

  IF (lcolumn) THEN
    CALL decompose (global_decomposition, 0, nproca, nprocb, ngl,      &
         nlon, nlev, nm, nn, nk, iadvec, norot=lrot, debug=debug,      &
         lfull_m=lfull_m, lats_1d=lat_1d(1), lons_1d=lon_1d(1))
  ELSE
    CALL decompose (global_decomposition, nproma, nproca, nprocb,      &
         ngl, nlon, nlev, nm, nn, nk, iadvec, norot=lrot, debug=debug, &
         lfull_m=lfull_m)
  ENDIF

  ! keep index values, not id values

  DO p = 1, p_nprocs
     global_decomposition(p)%mapmesh(:,:) =                            &
        indx(global_decomposition(p)%mapmesh(:,:),global_decomposition)
  END DO

  ! copy global decomposition table entry to local decomposition

  DO i = 1, p_nprocs
     IF (global_decomposition(i)% pe == p_pe) THEN
        local_decomposition = global_decomposition(i)
     END IF
  END DO

  IF (p_parallel_io) THEN
    WRITE(nout,*) '---------------------------------------------'
    WRITE(nout,*) ' Blocking information:'
    WRITE(nout,*) '   lreg  nproma ngpblks  npromz'
    WRITE(nout,'(L5,3I8)') local_decomposition% lreg,                  &
         local_decomposition% nproma,                                  &
         local_decomposition% ngpblks,                                 &
         local_decomposition% npromz
    WRITE(nout,*) '---------------------------------------------'
  END IF

#ifdef DEBUG
  ! decomposition printout

  IF (p_pe == p_io .AND. p_parallel) THEN
     DO i = 1, p_nprocs
        WRITE (nout,separator)
        CALL print_decomposition(global_decomposition(i))
     END DO
     WRITE (nout,separator)
     WRITE (nout,*)
  END IF
#endif

END SUBROUTINE init_decomposition
