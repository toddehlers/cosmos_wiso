MODULE mod_clim_def
!
!
!**** DEFINE
!
!     Purpose:
!      Routines CLIM_Define and CLIM_Defport are now in a module.
!     
!     Interface:
!       none
!    
!     Method:
!       Uses assumed shape array method to dimension local arrays.        
!
!     External:
!       none
!
!     Files:
!       none
!   
!     References:
!
!     History:
!     --------
!       Version   Programmer     Date        Description
!       ------------------------------------------------
!       2.5       A.Caubel       2002/06     created
!
!*-----------------------------------------------------------------------

CONTAINS
!
 SUBROUTINE CLIM_Define(id_nports, cdport,kinout,ktype,kparal,kinfo)
!
!*    *** Define ***   CLIM 2.0
!
!     purpose:
!     --------
!        define a port
!
!     interface:
!     ----------
!        id_nports : port number of the field
!        cdport : symbolic name of the field
!        kinout : port status in/out
!        ktype  : type of data
!        kparal : type of parallel decomposition
!	 kinfo	: output status
!
!     lib mp:
!     -------
!        mpi-1
!
!     author:
!     -------
!        Eric Sevault   - METEO FRANCE
!        Laurent Terray - CERFACS
!        Arnaud Caubel - Fecit
!     ----------------------------------------------------------------
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
   USE mod_kinds_oasis
   USE mod_clim
   USE mod_comclim
!     ----------------------------------------------------------------
   INTEGER(kind=ip_intwp_p)       kinout, ktype, kinfo, id_nports
   INTEGER(kind=ip_intwp_p), DIMENSION(:) :: kparal
   CHARACTER*(*) cdport
!     ----------------------------------------------------------------
!
!*    0. Entering
!     -----------
!
   kinfo = CLIM_Ok
!
!*    1. define the port
!     ------------------
!
   WRITE(nulprt,*) 'Define - port name, status : ',cdport,kinout
   WRITE(nulprt,*) 'Define - data type: ',ktype
!
   IF ( kinout.EQ.CLIM_InOut ) THEN
       CALL CLIM_Defport (id_nports,cdport, CLIM_In, ktype, kparal, kinfo )
       IF (kinfo.EQ.CLIM_Ok) THEN
           CALL CLIM_Defport (id_nports,cdport,CLIM_Out,ktype,kparal,kinfo )
       ENDIF
   ELSE
       CALL CLIM_Defport ( id_nports,cdport, kinout,  ktype, kparal, kinfo )
   ENDIF
!
!     ----------------------------------------------------------------
!
   CALL FLUSH(nulprt)
#endif
   RETURN
 END SUBROUTINE CLIM_Define
!
! ====================================================================

 SUBROUTINE CLIM_Defport(id_nports,cdport,kinout,ktype,kparal,kinfo)
!
!*    *** Define ***   CLIM 2.2
!
!     purpose:
!     --------
!        define a port
!
!     interface:
!     ----------
!        id_nports : port number of the field
!        cdport : symbolic name of the field
!        kinout : port status in/out
!        ktype  : type of data
!        kparal : type of parallel decomposition
!        kinfo  : output status
!
!     lib mp:
!     -------
!        mpi-2
!
!     author:
!     -------
!        Eric Sevault   - METEO FRANCE
!        Laurent Terray - CERFACS
!        Jean Latour    - F.S.E.   (mpi-2)
!
!     ----------------------------------------------------------------
#if defined use_comm_MPI1 || defined use_comm_MPI2 || (!defined use_comm_MPI1 && !defined use_comm_MPI2 && !defined use_comm_SIPC && !defined use_comm_GMEM && !defined use_comm_PIPE && !defined use_comm_NONE)
   USE mod_kinds_oasis
   USE mod_clim
   USE mod_comclim
!     ----------------------------------------------------------------
   INTEGER(kind=ip_intwp_p)       kinout, ktype, kinfo
   INTEGER(kind=ip_intwp_p), DIMENSION(:) :: kparal
   CHARACTER*(*) cdport
!     ----------------------------------------------------------------
   INTEGER(kind=ip_intwp_p)     ip, is, id_nports
   REAL(kind=ip_realwp_p)          rl_testvar
   REAL(kind=ip_double_p)  dl_testvar
   CHARACTER*32  cltest
   INTEGER (kind=ip_intwp_p) :: il_testvar, il_bytesize, il_iosize
   INTEGER (kind=ip_intwp_p) :: il_int_iosize, il_rbyt
!     ----------------------------------------------------------------
!
!*    1. check if this port already exist
!     -----------------------------------
!
   cltest= cdport
   il_testvar = 0
   rl_testvar = 0.0
   dl_testvar = 0.0
!
   DO ip=1,nports
     IF (cltest.EQ.cports(ip).AND.kinout.EQ.myport(1,ip)) THEN
         kinfo = CLIM_DoubleDef
         WRITE(nulprt,FMT='(A,A)') &
            'Define - WARNING - duplicate definition of port ', & 
            cdport
         GO TO 1010
     ENDIF
   ENDDO
!
!*    2. save arguments as half a link
!     --------------------------------
!
   nports = nports + 1
   id_nports = nports
   cports(nports) = cdport
!
   myport(1,nports) = kinout
   myport(5,nports) = 0
!
   IF ( ktype.EQ.CLIM_Real ) THEN
       myport(2,nports) = CLIM_Real
       il_bytesize = BIT_SIZE(il_testvar)/8
       INQUIRE (iolength=il_iosize) il_testvar
       il_int_iosize = il_iosize
       INQUIRE (iolength=il_iosize) rl_testvar
       il_rbyt = il_iosize/il_int_iosize*il_bytesize
       myport(3,nports) = il_rbyt
   ELSE
       kinfo = CLIM_BadType
       WRITE(nulprt,FMT='(A,I4)') &
          'Define - WARNING - Bad data type:',ktype
   ENDIF
!
   IF (kparal(CLIM_Strategy).EQ.CLIM_Serial) THEN
!
       mydist(CLIM_Strategy,nports)   = CLIM_Serial
       mydist(CLIM_Segments,nports)   = 1
       mydist(CLIM_Segments+1,nports) = 0
       mydist(CLIM_Segments+2,nports) = kparal(CLIM_Length)
       myport(4,nports) = kparal(CLIM_Length)
!
   ELSEIF (kparal(CLIM_Strategy).EQ.CLIM_Apple) THEN
!
       mydist(CLIM_Strategy,nports)   = CLIM_Apple
       mydist(CLIM_Segments,nports)   = 1
       mydist(CLIM_Segments+1,nports) = kparal(CLIM_Offset)
       mydist(CLIM_Segments+2,nports) = kparal(CLIM_Length)
       myport(4,nports) = kparal(CLIM_Length)
!
   ELSEIF (kparal(CLIM_strategy).EQ.CLIM_Box) THEN
!
       mydist(CLIM_Strategy,nports)   = CLIM_Box
       mydist(CLIM_Segments,nports)   = kparal(CLIM_SizeY)
       DO is=1,kparal(CLIM_SizeY)
         mydist(CLIM_Segments+2*is-1,nports) = &
            kparal(CLIM_Offset) + (is-1) * kparal(CLIM_LdX)
         mydist(CLIM_Segments+2*is,nports) = kparal(CLIM_SizeX)
       ENDDO
       myport(4,nports) = kparal(CLIM_SizeX) * kparal(CLIM_SizeY)
!
   ELSEIF (kparal(CLIM_strategy).EQ.CLIM_Orange) THEN
!
       mydist(CLIM_Strategy,nports)   = CLIM_Orange
       mydist(CLIM_Segments,nports)   = kparal(CLIM_Segments)
       myport(4,nports) = 0
       DO is=1,2*kparal(CLIM_Segments)
         mydist(CLIM_Segments+is,nports) = kparal(CLIM_Segments+is)
         IF (MOD(is,2).EQ.0) THEN
             myport(4,nports) = myport(4,nports) + &
                kparal(CLIM_Segments+is)
         ENDIF
       ENDDO
   ENDIF
!
!     ----------------------------------------------------------------
!
1010 CONTINUE
   CALL FLUSH(nulprt)
#endif
   RETURN
 END SUBROUTINE CLIM_Defport
END MODULE mod_clim_def
