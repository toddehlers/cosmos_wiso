      SUBROUTINE CLIM_Export(id_port_id,kstep,pfield,kinfo)
c
c*    *** Export ***   
c
c     purpose:
c     --------
c        give pfield to models connected to port id_port_id at time kstep
c
c     interface:
c     ----------
c        id_port_id : port number of the field
c	 kstep	: current time in seconds
c	 pfield	: buffer of reals
c	 kinfo	: output status
c
c     lib mp:
c     -------
c        mpi-1
c
c     author:
c     -------
c        Eric Sevault   - METEO FRANCE
c        Laurent Terray - CERFACS
c        Jean Latour    - F.S.E.     
c        Arnaud Caubel  - Fecit
c
c     ----------------------------------------------------------------
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
      USE mod_kinds_oasis
      USE mod_clim
      USE mod_comclim
#include <mpif.h>
c     ----------------------------------------------------------------
      INTEGER (kind=ip_intwp_p)       kstep, kinfo
      CHARACTER(len=1), DIMENSION(*) :: pfield
c     ----------------------------------------------------------------
      INTEGER (kind=ip_intwp_p)    info
      INTEGER (kind=ip_intwp_p)	   isend, ip, iport, ilk, iseg, is, 
     $    ilgb, imod, itid, itag, ilen, ioff, ityp, ibyt
      INTEGER (kind=ip_intwp_p)    iposbuf
c     ----------------------------------------------------------------
c
c*    0. Entering
c     -----------
c
C #slo - Output reduced
#ifndef __SILENT
      WRITE(nulprt, *)'entering CLIM_Export'
#endif
      kinfo = CLIM_Ok
c
c*    1. check for this port in my list
c     ---------------------------------
c
      isend = 0
      iport = -1
c
      IF (myport(1,id_port_id).eq.CLIM_Out) iport=id_port_id

      IF (iport.lt.0) THEN
         kinfo = CLIM_BadPort
         WRITE(nulprt,FMT='(A,A)')
     *         'Export - WARNING - Invalid port out: ',
     $        cports(id_port_id)
         GO TO 1010
      ENDIF
c
c*    2. check for connected ports (in)
c     ---------------------------------
c
#ifndef __SILENT
      WRITE(nulprt,FMT='(A,A)') 'Export - ', cports(iport)
#endif
c
      ityp = myport(2,iport)
      ibyt = myport(3,iport)
c
      DO 290 ip=1,myport(5,iport)
c
         ilk  = myport(5+ip,iport)
         imod = mylink(1,ilk)
         itid = mylink(2,ilk)
         itag = mylink(3,ilk) - kstep / ig_frqmin
         iseg = mylink(4,ilk)
c
             ilgb = 0
             iposbuf = 0
             DO 240 is=1,iseg
               ioff = mylink(4+2*is-1,ilk) * ibyt + 1
               ilen = mylink(4+2*is,ilk)
c
               IF ( ityp .EQ. CLIM_Real ) THEN
                   IF ( ip_realwp_p .eq. ip_single_p ) THEN
                       CALL MPI_Pack ( pfield(ioff), ilen, 
     *                     MPI_REAL,pkwork, ig_maxtype, iposbuf, 
     *                     mpi_comm, info )
                   ELSE IF ( ip_realwp_p .eq. ip_double_p ) THEN
                       CALL MPI_Pack ( pfield(ioff), ilen, 
     *                     MPI_DOUBLE_PRECISION,pkwork, ig_maxtype, 
     *                     iposbuf, mpi_comm, info )
                   ENDIF
               ELSE
                   WRITE(nulprt,*)'Export - pb type incorrect ',ityp
		   kinfo = CLIM_BadType
                   GO TO 1010
               ENDIF
               ilgb = ilgb + ilen
 240        CONTINUE
            IF (info.ne.0 .or. ilgb*ibyt .gt. ig_maxtype) THEN
               kinfo = CLIM_Pack
               WRITE(nulprt,FMT='(A,I3,I8,A)')
     *              'Export - pb pack<mpi ',info,ilgb*ibyt,'>'
            ELSE
c
                IF (lg_bsend) then
c*
c*   Buffered send
c*   -> if fields are not sent and received in the same order, and
c*   and on architectures on which MPI_Send is not implemented with a 
c*   mailbox (e.g. NEC SX5)
c*
                    CALL MPI_BSend ( pkwork, iposbuf, MPI_PACKED, itid,
     *                  itag, mpi_comm, info )
                ELSE
c
c*   Standard blocking send: To be used
c*   -> if fields are necessarily sent and received in the same order, 
c*   -> or on architectures on which MPI_Send is implemented with a 
c*      mailbox (e.g. VPPs); in this case, make sure that your mailbox
c*      size is large enough.
c
                    CALL MPI_Send ( pkwork, iposbuf, MPI_PACKED, itid,
     *                  itag, mpi_comm, info )
c
                ENDIF
c
                 IF (info.eq.CLIM_ok) THEN
                     isend = isend + 1
                     nbsend = nbsend + ilgb * ibyt
C #slo - Output reduced
#ifndef __SILENT
                     WRITE(nulprt,FMT='(A,I2,A,I9,A,I7,A,I2,A,I10,A)')
     *                   'Export - <dest:',imod,
     *                   '> <step:',kstep,
     *                   '> <len:',ilgb,
     *                   '> <type:',ibyt,
     *                   '> <tag:',itag,'>'
#endif
                 ELSE
                     kinfo = CLIM_Pvm
                     WRITE(nulprt,FMT='(A,I3,A)')
     *                   'Export - pb send <mpi ',info,'>'
                 ENDIF
             ENDIF
c
c
 290    CONTINUE
c
#ifndef __SILENT
        WRITE(nulprt,FMT='(A,I3,A)') 
     *     'Export - ',isend,' fields exported'
#endif
c
c     ----------------------------------------------------------------
c
 1010 CONTINUE
      CALL FLUSH(nulprt)
#endif
      RETURN
      END
