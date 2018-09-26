      SUBROUTINE CLIM_Start(kinfo)
c
c*    *** Start ***   CLIM 3.0
c
c     purpose:
c     --------
c        beginning of the coupled run
c
c     interface:
c     ----------
c	 kinfo	: output status
c
c     lib mp:
c     -------
c        mpi-2
c
c     author:
c     -------
c        Eric Sevault   - METEO FRANCE
c        Laurent Terray - CERFACS
c        Jean Latour    - F.S.E.   (mpi-2)
c        Arnaud Caubel  - FECIT - removed kmxtag as argument
c
c     ----------------------------------------------------------------
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
      USE mod_kinds_oasis
      USE mod_clim
      USE mod_comclim
c     ----------------------------------------------------------------
      INTEGER (kind=ip_intwp_p) kinfo
c     ----------------------------------------------------------------
c
      kinfo  = CLIM_Ok
      CALL CLIM_Start_MPI(kinfo)
      WRITE (nulprt,FMT='(A)') 'Returning from Start -- '
      CALL FLUSH(nulprt)
c
#endif
      RETURN
      END
