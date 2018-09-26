      SUBROUTINE CLIM_Stepi(cdnam,kinfo)
c
c*    *** Stepi ***   CLIM 2.0
c
c     purpose:
c     --------
c        retrieve info on cdnam models time steps
c
c     interface:
c     ----------
c        cdnam   : character string name of the model
c        kinfo  : output status or nbr of procs implied 
c                 in the coupling for cdnam
c
c     lib mp:
c     -------
c        mpi-2
c
c     author:
c     -------
c        Eric Sevault   - METEO FRANCE
c        Laurent Terray - CERFACS
c        Arnaud Caubel - FECIT (08/02 - removed some arguments)
c     ----------------------------------------------------------------
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
      USE mod_kinds_oasis
      USE mod_clim
      USE mod_comclim
c     ----------------------------------------------------------------
      CHARACTER*(*) cdnam
      INTEGER (kind=ip_intwp_p)	    kinfo
c     ----------------------------------------------------------------
      CHARACTER*32  cltest
      INTEGER (kind=ip_intwp_p)       ip
c     ----------------------------------------------------------------
c
c*    0. Entering
c     -----------
c
      WRITE (nulprt, FMT='(A)')'Stepi called'
      kinfo = CLIM_Ok
c
c*    1. check for the cdnam model in my list
c     --------------------------------------
c
      kinfo = CLIM_BadName
c
      cltest=cdnam
c
      kinfo=0
      DO 100 ip=0,ncplprocs-1
         IF ( cnames(ip).eq.cltest ) THEN
C IF more than one process for the model, 
C info coming from last proc is kept as output
            kinfo = kinfo + 1
         ENDIF
  100 CONTINUE
c
c     ----------------------------------------------------------------
c
      CALL FLUSH(nulprt)
#endif
      RETURN
      END
