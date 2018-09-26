      SUBROUTINE CLIM_Import(id_port_id,kstep,pfield,kinfo)
c
c*    *** Import ***   
c
c     purpose:
c     --------
c        recv pfield from models connected to port id_port_id
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
      INTEGER (kind=ip_intwp_p)     i4kstep
      INTEGER (kind=ip_intwp_p)     info
      INTEGER (kind=ip_intwp_p)	    irecv, imod, ilk, iseg, is, ilgb,
     $    itid, itag, ilen, ioff, ityp, ibyt
      INTEGER (kind=ip_intwp_p)     iposbuf, istatus(MPI_STATUS_SIZE), 
     $    imaxbyt
c
#ifdef __DEBUG
      INTEGER(kind=ip_intwp_p)     icount
      INTEGER(kind=ip_intwp_p), parameter :: icountmax=600
      LOGICAL                      iflag
#endif
c     ----------------------------------------------------------------
c
c*    0. Entering
c     -----------
c
      i4kstep = kstep
      istatus(:)=0
      kinfo = CLIM_Ok
C #slo - Output reduced
#ifndef __SILENT
      WRITE(nulprt, *)'entering CLIM_Import'
#endif
c
c*    1. check for this port in my list
c     ---------------------------------
c
      irecv = 0
      iport = -1
c
      IF (myport(1,id_port_id).eq.CLIM_In) iport=id_port_id
      IF (iport.lt.0) THEN
         kinfo = CLIM_BadPort
         WRITE(nulprt,FMT='(A,A)')
     *         'Import - WARNING - Invalid port out: ', 
     $        cports(id_port_id)
         GO TO 1010
      ENDIF
c
c*    2. check for connected ports (in)
c     ---------------------------------
c
#ifndef __SILENT
      WRITE(nulprt,FMT='(A,A)') 'Import - ', cports(iport)
      call flush(nulprt)
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
         itag = mylink(3,ilk) - i4kstep / ig_frqmin
         iseg = mylink(4,ilk)
c
c*   Implementation with "blocking" receives : the program will wait
c*   indefinitely until a message is received (this may generate a
c*   deadlock if all models are waiting on a receive).
c*   However this method will be more efficient in most cases than the
c*   receives with a time-out loop. 
c
#ifdef __DEBUG
c
c added to detect deadlocks (J. Latour, 2004-04-28)
c
         CALL MPI_Iprobe ( itid, itag, mpi_comm, iflag, istatus, info )
         WRITE(nulprt,*) 'probing for tid = ',itid,' tag = ',itag,
     *       ' comm = ',mpi_comm,' result is : ',iflag
         CALL flush(nulprt) 

         IF (.NOT.iflag) THEN
             icount = 0
             WAITLOOP:  DO
             CALL  MPI_Iprobe (itid,itag,mpi_comm,iflag,istatus,info )
             icount = icount + 1
             IF ( iflag ) EXIT WAITLOOP
             IF ( icount .GE. icountmax ) THEN
                 WRITE(nulprt,*) 'probing for tid =',itid,'tag= ',itag,
     *               ' still negative after ',icountmax,
     *               ' seconds : Abort the job'
                 CALL flush(nulprt)
                 CALL MPI_ABORT (mpi_comm, 0, mpi_err)
             ENDIF
             call sleep(1)
             END DO WAITLOOP
             WRITE(nulprt,*)'probing for tid= ',itid,'icount= ', icount
             call flush(nulprt)
         ENDIF
#endif
            CALL MPI_Recv ( pkwork, ig_maxtype, MPI_PACKED, itid,
     *                       itag, mpi_comm, istatus, info ) 
            CALL MPI_Get_count ( istatus, MPI_PACKED, imaxbyt, info )
#ifndef __SILENT
            WRITE(nulprt,*)'imaxbyt ',imaxbyt
            WRITE(nulprt,*)'info ',info 
            call flush(nulprt)
#endif
c
            IF ( info .EQ. CLIM_ok  .AND.  imaxbyt .GT. 0) THEN
               ilgb = 0
               iposbuf = 0
               DO 260 is=1,iseg
                  ioff = mylink(4+2*is-1,ilk) * ibyt + 1
                  ilen = mylink(4+2*is,ilk)
c
                  IF ( ityp .EQ. CLIM_Real ) THEN
                     IF ( ip_realwp_p .eq. ip_single_p ) THEN
                         CALL MPI_Unpack ( pkwork, ig_maxtype, iposbuf,
     *                       pfield(ioff), ilen, MPI_REAL, 
     *                       mpi_comm, info)
                     ELSE IF ( ip_realwp_p .eq. ip_double_p ) THEN
                         CALL MPI_Unpack ( pkwork, ig_maxtype, iposbuf,
     *                       pfield(ioff), ilen, MPI_DOUBLE_PRECISION, 
     *                       mpi_comm, info)
                     ENDIF
                 ELSE
                     WRITE(nulprt,*)'Import - pb type incorrect ',ityp
                     kinfo = CLIM_BadType
                     GO TO 1010
                 ENDIF
                 ilgb = ilgb + ilen
 260           CONTINUE
               IF (ilgb*ibyt .le. imaxbyt) THEN
                   irecv  = irecv + 1
                   nbrecv = nbrecv + ilgb * ibyt
#ifndef __SILENT
                   WRITE(nulprt,FMT='(A,I2,A,I9,A,I7,A,I2,A,I10,A)')
     *                 'Import - <from:',imod,
     *                 '> <step:',kstep,
     *                 '> <len:',ilgb,
     *                 '> <type:',ibyt,
     *                 '> <tag:',itag,'>'
#endif
               ELSE
                   kinfo = CLIM_Unpack
                   WRITE(nulprt,FMT='(A,I3,A)')
     *                 'Import - pb unpack <mpi ',info,'>'
               ENDIF
           ELSE
               kinfo = CLIM_TimeOut
               WRITE(nulprt,FMT='(A,I3,A)')
     *              'Import - abnormal exit from trecv <mpi ',info,'>'
           ENDIF
c
 290  CONTINUE
c
#ifndef __SILENT
      WRITE(nulprt,FMT='(A,I3,A)')
     *     'Import - ',irecv,' fields imported'
#endif
c
c     ----------------------------------------------------------------
c
 1010 CONTINUE
      CALL FLUSH(nulprt)
#endif
      RETURN
      END
