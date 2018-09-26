      SUBROUTINE CLIM_Quit(kinfo)
c
c*    *** Quit ***   CLIM 3.0
c
c     purpose:
c     --------
c        leave mpi
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
c        Arnaud Caubel - Fecit - Added deallocation of CLIM arrays
c     ----------------------------------------------------------------
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
      USE mod_kinds_oasis
      USE mod_clim
      USE mod_comclim
#include <mpif.h>
c     ----------------------------------------------------------------
      INTEGER (kind=ip_intwp_p)	kinfo
c     ----------------------------------------------------------------
      INTEGER (kind=ip_intwp_p) info, il_rank, ji
c     ----------------------------------------------------------------
      INTEGER (kind=ip_intwp_p)       il_bufsizebyt
      REAL (kind=ip_double_p)  dl_bufaddr
c
c*    0. Entering
c     -----------
c
      WRITE(nulprt,*)'entering clim_quit'
      kinfo = CLIM_Ok
c
      CALL MPI_Comm_Rank(mpi_comm,il_rank,mpi_err)
c
      DO ji = 0, ncplprocs-1 
        IF (il_rank.eq.modtid(ji)) THEN
            DEALLOCATE (ncode, stat=il_err)
            IF (il_ERR.ne.0) WRITE(nulprt,*)'Error in "ncode" 
     $          deallocation in CLIM_Quit routine ! '
            DEALLOCATE (cnames, stat=il_err)
            IF (il_ERR.ne.0) WRITE(nulprt,*)'Error in "cnames" 
     $          deallocation in CLIM_Quit routine ! '
            DEALLOCATE (myport, stat = il_ERR)
            IF (il_ERR.ne.CLIM_Ok) WRITE(nulprt,*)' Problem 
     $          in "myport" deallocation in CLIM_Quit !'
            DEALLOCATE (mydist,stat = il_ERR)
            IF (il_ERR.ne.CLIM_Ok) WRITE(nulprt,*)' Problem 
     $          in "mydist" deallocation in CLIM_Quit !'
            DEALLOCATE (cports, stat = il_ERR)
            IF (il_ERR.ne.CLIM_Ok) WRITE(nulprt,*)' Problem 
     $          in "cports" allocation in CLIM_Quit !'
            DEALLOCATE(clrport, stat = il_ERR)
            IF (il_ERR.ne.CLIM_Ok) WRITE(nulprt,*)'Error in clrport 
     $          deallocation in CLIM_Quit !'
            DEALLOCATE(irdist, stat = il_ERR)
            IF (il_ERR.ne.0) WRITE(nulprt,*)'Error in "irdist" 
     $          deallocation in CLIM_Quit'
            DEALLOCATE(irport, stat = il_ERR)
            IF (il_ERR.ne.0) WRITE(nulprt,*)'Error in "irport" 
     $          deallocation in CLIM_Quit'
            DEALLOCATE (mylink, stat = il_ERR)
            IF (il_ERR.ne.CLIM_Ok) WRITE(nulprt,*)' Problem 
     $          in "mylink" allocation in CLIM_Quit !'
            DEALLOCATE (pkwork, stat = il_ERR)
            IF (il_ERR.ne.CLIM_Ok) WRITE(nulprt,*)' Problem 
     $          in "pkwork" allocation in CLIM_Quit !'
            IF (lg_bsend) then
               CALL MPI_Buffer_Detach(dl_bufaddr,il_bufsizebyt,mpi_err)
               DEALLOCATE (dg_bufsend, stat = il_ERR)
               IF (il_ERR.ne.CLIM_Ok) WRITE(nulprt,*)' Problem 
     $             in "dg_bufsend" allocation in CLIM_Quit !'
            ENDIF
        ENDIF
      END DO
      DEALLOCATE (modtid, stat=il_err)
      IF (il_ERR.ne.0) WRITE(nulprt,*)'Error in "modtid" 
     $    deallocation in CLIM_Quit routine ! '
c
c
      CALL MPI_Finalize ( info )
c
      WRITE(nulprt,FMT='(A,I3,A)')
     *     'Quit - exit status <mpi ',info,'>'
c
c     ----------------------------------------------------------------
c
      CLOSE (nulprt)
#endif
      RETURN
      END
