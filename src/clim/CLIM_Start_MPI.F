      SUBROUTINE CLIM_Start_MPI(kinfo)
c
c*    *** Start-mpi ***   CLIM 3.0
c
c     purpose:
c     --------
c        beginning of the coupled run (MPI-2 only)
c
c     interface:
c     ----------
c	 kinfo	: output status
c
c     method:
c     -------
c        With MPI2, use the new intra-communicator created by the
c                   MPI_Intercomm_merge in CLIM_Init (models)
c                   or in inicmc via CLIM_Init (coupler Oasis)
c
c     lib mp:
c     -------
c        mpi-1 (or mpi-2)
c
c     author:
c     -------
c        Eric Sevault   - METEO FRANCE
c        Laurent Terray - CERFACS
c        Jean Latour    - F.S.E.   (mpi-2)
c        Sophie Valcke  - CERFACS (08/09/00 -modified)
c        Arnaud Caubel  - FECIT (08/02 - removed kmxtag as argument)
c
c     ----------------------------------------------------------------
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
      USE mod_kinds_oasis
      USE mod_clim 
      USE mod_comclim
#include <mpif.h>
c     ----------------------------------------------------------------
      INTEGER (kind=ip_intwp_p) kinfo
c     ----------------------------------------------------------------
      INTEGER (kind=ip_intwp_p)    ip, info1, info, ilgdt, ipos, is,
     $    itag1, isdmod, incp, ipotag,nsegid
      INTEGER (kind=ip_intwp_p)    irempo, il_maxtag
c
      INTEGER (kind=ip_intwp_p), PARAMETER ::    ibuflen=3
      INTEGER (kind=ip_intwp_p)    iposbuf, ibuff_mpi(ibuflen) 
      INTEGER (kind=ip_intwp_p)    istatus(MPI_STATUS_SIZE), imaxbyt
      REAL (kind=ip_realwp_p) rl_testvar
      LOGICAL      ll_flag
      INTEGER (kind=ip_intwp_p), allocatable :: ireq(:)
      INTEGER (kind=ip_intwp_p)              :: ireqrecv
      INTEGER (kind=ip_intwp_p)              :: irc
      INTEGER (kind=ip_intwp_p)              :: rpkwork(ig_CLIMmax)
      INTEGER (kind=ip_intwp_p) :: ii, il_bytesize, il_iosize
      INTEGER (kind=ip_intwp_p) :: il_int_iosize, il_rbyt
c     ----------------------------------------------------------------
c
c
c     0. define some variables 
c     ------------------------
c
      WRITE(nulprt,*)'Start - -  '
c           
      istatus(:)=0
      ibuff_mpi(:)=0
      rl_testvar=0.0_ip_realwp_p
c
      call MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, CLIM_MaxTag, 
     $    ll_flag,info)
      if(info.ne.MPI_SUCCESS) then
        write(nulprt, * )'CLIM_Start_MPI -- MPI_Attr_get failed !',info
        call MPI_ABORT (mpi_comm, 0, info)
      endif
      if ( .not. ll_flag ) then
         write(nulprt, * ) 'Warning:  MPI_Attr_get did not return '
         write(nulprt, * ) ' a valid CLIM_MaxTag!'
         write(nulprt, * ) ' CLIM_MaxTag is set to 32767!'
         CLIM_MaxTag = 32767
      else
         write(nulprt, * ) 'CLIM_MaxTag is now ', CLIM_MaxTag
!
!     limit CLIM_MaxTag to 2^^30 - 1 (J. Latour)
!
         if ( CLIM_MaxTag .gt. 1073741823 ) CLIM_MaxTag =  1073741823
         write(nulprt, * ) 'CLIM_MaxTag is now ', CLIM_MaxTag
!
      endif
      il_maxtag = CLIM_MaxTag - 1
      itag1 = CLIM_MaxTag
      ilgdt = CLIM_ParSize
C
      il_bytesize = BIT_SIZE(ii)/8
      INQUIRE (iolength=il_iosize) ii
      il_int_iosize = il_iosize
      INQUIRE (iolength=il_iosize) rl_testvar
      il_rbyt = il_iosize/il_int_iosize*il_bytesize
      ig_maxtype = ig_CLIMmax*il_rbyt
C
c*    1. broadcast usefull informations
c     ---------------------------------
c
      ibuff_mpi(1) = mynum 
      ibuff_mpi(2) = modtid(mynum) 
      ibuff_mpi(3) = nports
c
      WRITE(nulprt,FMT='(A,A)') 'Start- send- MODEL ',cmynam
      WRITE(nulprt,FMT='(A,I9)')'Start- send-  num  :',mynum
      WRITE(nulprt,FMT='(A,I9)')'Start- send-  tid  :',modtid(mynum)
      WRITE(nulprt,FMT='(A,I9)')'Start- send-  nport:',nports
c
      iposbuf = 0
      call MPI_Pack ( ibuff_mpi, ibuflen, MPI_INTEGER, pkwork,
     *                ig_maxtype, iposbuf, mpi_comm, info )
      call MPI_Pack ( cmynam, CLIM_Clength, MPI_CHARACTER, pkwork,
     *                ig_maxtype, iposbuf, mpi_comm, info )
      DO 210 ji=1,nports
         call MPI_Pack ( cports(ji), CLIM_Clength, MPI_CHARACTER,
     *           pkwork, ig_maxtype, iposbuf, mpi_comm, info )
         call MPI_Pack ( myport(1,ji), 5, MPI_INTEGER, pkwork,
     *                   ig_maxtype, iposbuf, mpi_comm, info )
         call MPI_Pack ( mydist(1,ji), ilgdt, MPI_INTEGER, pkwork,
     *                   ig_maxtype, iposbuf, mpi_comm, info )
 210   CONTINUE 
c
       ALLOCATE(ireq(0:ncplprocs))
       irc = 0
c
      DO 220 ji = 0 , ncplprocs-1
c
c       Send to all processors involved in the coupling except itself
        IF ( ji .NE. mynum ) THEN
            CALL MPI_Isend ( pkwork, iposbuf, MPI_PACKED, modtid(ji),
     *                      itag1, mpi_comm, ireq(irc), info )
          irc = irc+1
c
        ENDIF 
 220  CONTINUE

      WRITE (nulprt,FMT='(A,I3,A,I8,A)') 
     * 'Start - broadcast from mynum = ',mynum,' <MPI ',info,'>'
c
c*    3. get these infos from other models and check ports
c     ----------------------------------------------------
c
      imaxbyt=0
      ireqrecv = 0
      DO 380 ip=1,ncplprocs-1
c
        CALL MPI_Irecv(rpkwork, ig_maxtype, MPI_PACKED, MPI_ANY_SOURCE,
     *                 itag1, mpi_comm, ireqrecv, info1 )
        CALL MPI_Wait ( ireqrecv, istatus, info )
        CALL MPI_Get_count ( istatus, MPI_PACKED, imaxbyt, info )
        IF ( info .NE. 0 ) THEN
            kinfo = CLIM_Unpack
            WRITE(nulprt,FMT='(A,I3,A)')
     *                 'Import - pb unpack <mpi ',info,'>'
            GO TO 1010
        ENDIF
c END of receive
c
c       
        iposbuf = 0
        CALL MPI_Unpack ( rpkwork, ig_maxtype, iposbuf, ibuff_mpi,
     *                      ibuflen, MPI_INTEGER, mpi_comm, info )
        IF ( info .NE. 0 ) THEN
            kinfo = CLIM_Unpack
            WRITE(nulprt,FMT='(A,I3,A)')
     *                 'Import - pb unpack <mpi ',info,'>'
            GO TO 1010
        ENDIF
c
        isdmod         = ibuff_mpi(1)
        modtid(isdmod) = ibuff_mpi(2)
        irempo         = ibuff_mpi(3)
c
        call MPI_Unpack ( rpkwork, ig_maxtype, iposbuf, cnames(isdmod),
     *                   CLIM_Clength, MPI_CHARACTER, mpi_comm, info )
        DO 310 ji=1,irempo
          call MPI_Unpack ( rpkwork, ig_maxtype, iposbuf, clrport(ji),
     *                   CLIM_Clength, MPI_CHARACTER, mpi_comm, info )
          call MPI_Unpack( rpkwork, ig_maxtype, iposbuf, irport(1,ji), 
     *                           5, MPI_INTEGER, mpi_comm, info )
          call MPI_Unpack( rpkwork, ig_maxtype, iposbuf, irdist(1,ji),
     *                           ilgdt, MPI_INTEGER, mpi_comm, info )
 310    CONTINUE
c
        WRITE (nulprt,FMT='(A,A)')  'Start - MODEL ',cnames(isdmod)
        WRITE (nulprt,FMT='(A,I9)') 'Start -   num  :',isdmod
        WRITE (nulprt,FMT='(A,I9)') 'Start -   tid  :',modtid(isdmod)
        WRITE (nulprt,FMT='(A,I9)') 'Start -   nport:',irempo
c
        ncode(isdmod) = 0
c
        DO 350 ji=1,nports
          DO 340 jj=1,irempo
            IF (cports(ji).eq.clrport(jj).and.
     *          myport(1,ji)+irport(1,jj).eq.1) THEN
                IF ((mydist(1,ji).eq.irdist(1,jj).and.
     *              myport(4,ji).ne.irport(4,jj)).or.
     *              (mydist(CLIM_Strategy,ji).eq.CLIM_Serial.and.
     *              myport(4,ji).lt.irport(4,jj)).or.
     *              (irdist(CLIM_Strategy,jj).eq.CLIM_Serial.and.
     *              irport(4,jj).lt.myport(4,ji))) THEN
                    kinfo = CLIM_IncSize
                    WRITE(nulprt,FMT='(A,A,A,I2,A,I2,A,I2,I2)')
     *                  'Start - WARNING Incompatible sizes - field',
     *                  cports(ji),'model ',mynum,' and model ',
     *                  isdmod,': ',myport(4,ji),irport(4,jj)
                ELSEIF (myport(2,ji).ne.irport(2,jj)) THEN
                    kinfo = CLIM_BadType
                    WRITE(nulprt,FMT='(A,A,A,I2,A,I2,A,I2,I2)')
     *                  'Start - WARNING Incompatible types - field',
     *                  cports(ji),'model ',mynum,' and model ',
     *                  isdmod,': ',myport(2,ji),irport(2,jj)
                ELSE
                    IF (myport(1,ji).eq.CLIM_Out) THEN
                        incp   = ig_ntime / ig_frqmin + 1
                        ipotag = CLIM_MaxTag - (ji-1)*incp - 1
                    ELSE
                        incp   = ig_ntime / ig_frqmin + 1
                        ipotag = CLIM_MaxTag - (jj-1)*incp - 1
                    ENDIF
                    IF ((ipotag-incp).lt.il_maxtag) THEN
                        il_maxtag = ipotag - incp
                    ENDIF
                    IF (mydist(CLIM_Strategy,ji).ne.CLIM_Serial.and
     *               .irdist(CLIM_Strategy,jj).ne.CLIM_Serial) THEN
                        IF (mydist(CLIM_Strategy,ji).ne.
     *                      irdist(CLIM_Strategy,jj)) THEN
                            CALL halte('STOP in Clim_Start_MPI')
                        ELSE
                            IF (mydist(CLIM_Segments,ji).eq.
     *                      irdist(CLIM_Segments,jj)) THEN
                                nsegid = 0                            
                                DO  is=1,2*mydist(CLIM_Segments,ji)
                                  IF (mydist(CLIM_Segments+is,ji).eq.
     *                               irdist(CLIM_Segments+is,jj)) THEN 
                                      nsegid = nsegid + 1
                                  ENDIF
                                END DO
                                IF (nsegid.eq.
     *                              (2*mydist(CLIM_Segments,ji))) THEN
                                    nlinks = nlinks + 1
                                    ipos   = 5
                                    myport(ipos,ji)=myport(ipos,ji)+ 1
                                    myport(ipos+myport(ipos,ji),ji)
     *                                     = nlinks
                                    mylink(1,nlinks) = isdmod
                                    mylink(2,nlinks) = modtid(isdmod)
                                    mylink(3,nlinks) = ipotag
                                    mylink(4,nlinks) = 1
                                    mylink(5,nlinks) = 0
                                    mylink(6,nlinks) = myport(4,ji)
                                    WRITE(nulprt,FMT='(A,A)')
     *                                  'Start - LINK ',cports(ji)
                                    WRITE(nulprt,FMT='(A,I2,A,I1,A,I2,
     *                                 A,I1,A,I10,A,I10)')
     *                                  'Start - [model ',
     *                                  mynum,'/io ',myport(1,ji),'] - 
     *                                    [model ',
     *                                  isdmod,'/io ',irport(1,jj),'] 
     *                                   using tags ',
     *                                  ipotag,' to ',ipotag-incp+1
                                ENDIF    
                            ENDIF
                        ENDIF
                    ELSE                    
                        nlinks = nlinks + 1
                        ipos   = 5
                        myport(ipos,ji) = myport(ipos,ji) + 1
                        myport(ipos+myport(ipos,ji),ji) = nlinks
                        mylink(1,nlinks) = isdmod
                        mylink(2,nlinks) = modtid(isdmod)
                        mylink(3,nlinks) = ipotag
                        mylink(4,nlinks) = 1
                        mylink(5,nlinks) = 0
                        mylink(6,nlinks) = myport(4,ji)
                        IF (mydist(CLIM_Strategy,ji).eq.CLIM_Serial.and
     *                   .irdist(CLIM_Strategy,jj).ne.CLIM_Serial) THEN
                            mylink(4,nlinks) = irdist(CLIM_Segments,jj)
                            DO 330 is=1,2*irdist(CLIM_Segments,jj)
                              mylink(4+is,nlinks) = irdist
     *                                      (CLIM_Segments+is,jj)
 330                        CONTINUE
                        ENDIF
                        WRITE(nulprt,FMT='(A,A)')
     *                      'Start - LINK ',cports(ji)
                        WRITE(nulprt,FMT='(A,I2,A,I1,A,I2,A,
     *                       I1,A,I10,A,I10)')
     *                      'Start - [model ',
     *                      mynum,'/io ',myport(1,ji),'] - [model ',
     *                      isdmod,'/io ',irport(1,jj),'] using tags ',
     *                      ipotag,' to ',ipotag-incp+1
                    ENDIF                    
                ENDIF
            ENDIF
 340      CONTINUE
 350    CONTINUE
c
        WRITE (nulprt,FMT='(A)') 'Start - -'
 380  CONTINUE
c
c     MPI_wait on the above MPI_Isend so to not change the content of pkwork
c     before the sending is complete.
      irc = 0
      DO 420 ji = 0 , ncplprocs-1
        IF ( ji .NE. mynum ) THEN
            CALL MPI_Wait ( ireq(irc), istatus, info )
            WRITE (nulprt, FMT='(A)') 'After MPI_wait'
            irc = irc+1
        ENDIF 
 420  CONTINUE
      DEALLOCATE(ireq)
c
c     ----------------------------------------------------------------
c
 1010 CONTINUE
      WRITE (nulprt,FMT='(A)') 'Returning from Start-mpi  -- '
      CALL FLUSH(nulprt)
#endif
      RETURN
      END
