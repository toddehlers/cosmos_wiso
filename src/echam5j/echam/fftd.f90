#if defined (__SX__) && defined (_OPENMP)
#define VECTOR 1
#endif
SUBROUTINE fftd

#ifdef FFT991
  USE mo_kind,          ONLY: dp
  USE mo_fft991,        ONLY: fft991cy
#else
  USE mo_fft992,        ONLY: fft992
#endif
  USE mo_buffer_fft,    ONLY: fftz, nvar
  USE mo_decomposition, ONLY: dc => local_decomposition
#if defined (VECTOR) && defined (_OPENMP)
  USE omp_lib,          ONLY: omp_get_thread_num, &
                              omp_get_max_threads
#endif

  IMPLICIT NONE

  !  Local scalars:
  INTEGER :: inc, isign
  INTEGER :: nlon, nlp2, nlev, nlat
#if (! defined FFT991) && (defined _OPENMP)
#ifdef VECTOR
  INTEGER :: ivar, tid, nthreads, chunk, rest
  INTEGER, ALLOCATABLE, SAVE :: istart(:), icount(:)
#else
  INTEGER :: ilat, ivar
#endif
#endif
  LOGICAL :: col_1d

  !  Local arrays:
#ifdef FFT991
  REAL(dp) :: zwork((dc%nlon+2) * dc%nflevp1 * dc%nflat * nvar)
#endif

!-- 2. direct *Fourier transforms

!-- 2.1 Set constants

  inc    = 1
  isign  = -1
  nlon   = dc% nlon
  nlp2   = nlon + 2
  nlev   = dc %nflevp1
  nlat   = dc% nflat
  col_1d = dc% col_1d

!-- 2.2 fft(*vo*,*d*,*t*,*alps*,*u*,*v*,*dtl*,*dtm*,*dalpsl*,*dalpsm*)

  IF (.NOT.col_1d) THEN
#ifdef FFT991
    CALL fft991cy(fftz,zwork,inc,nlp2,nlon,nvar*nlev*nlat,isign)
#else
#ifdef _OPENMP
#ifdef VECTOR
    IF (.NOT. ALLOCATED(istart)) THEN
      nthreads = omp_get_max_threads()
      ALLOCATE(istart(0:nthreads), icount(0:nthreads))
      istart(0) = 1
      DO tid = 0, nthreads-1
        chunk = nlat/nthreads 
        rest  = MOD(nlat, nthreads)
        if (tid < rest) chunk = chunk+1
        icount(tid) = chunk
        istart(tid+1) = istart(tid)+chunk
      ENDDO
    ENDIF
!$OMP PARALLEL PRIVATE(tid)
    tid = omp_get_thread_num()
    DO ivar = 1, nvar
        CALL fft992(fftz(1,1,istart(tid),ivar),inc,nlp2,nlon,nlev*icount(tid),isign)
    END DO
!$OMP END PARALLEL
#else
!$OMP PARALLEL PRIVATE(ilat)
    DO ivar = 1, nvar
!$OMP DO
      DO ilat = 1, nlat
        CALL fft992(fftz(1,1,ilat,ivar),inc,nlp2,nlon,nlev,isign)
      ENDDO
!$OMP END DO
    ENDDO
!$OMP END PARALLEL
#endif
#else
    CALL fft992(fftz,inc,nlp2,nlon,nvar*nlev*nlat,isign)
#endif
#endif
  ENDIF

END SUBROUTINE fftd
