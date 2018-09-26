      SUBROUTINE BODEN
!*************************************************************
!
!     BBBBB    OOO   DDDDD   EEEEEE  NN   N
!     B    B  O   O  D    D  E       NNN  N
!     BBBBB   O   O  D    D  EEEEE   N NN N
!     B    B  O   O  D    D  E       N  NNN
!     BBBBB    OOO   DDDDD   EEEEEE  N   NN
!
!
!--------------------------------------------------------------
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_UNITS
      IMPLICIT NONE
!
      CHARACTER*6 BLORD
     
!      CHARACTER*1 CDEUTO(IE_G,JE_G)
!      REAL DEUTO_G(IE_G,JE_G), DEPTO_G(IE_G,JE_G)
      REAL AMSUE_G(IE_G,JE_G), AMSUO_G(IE_G,JE_G)
      REAL CMUE,CMUO,CMWE,CMWO,CUE,CUO,CWE,CWO,DEPALT,DEPMAX,DEPMIN
      REAL SUGGO,SUM,SUPPO,TMAX, dpp,dpp13
      INTEGER I,I1,I2,IANF,IEND,III,IPR,IPRANZ,J,JEV,JJ,JJJ,K
      INTEGER LAUF,NDIF1,ir,il,iu,ju
#ifdef bounds_exch_tp
      integer, pointer :: mus(:,:)
#endif

!------------------------------------------------------------------------
!
!   HERE   :  READ TOPOGRAPHY FROM EXTERNAL FILE ON ARRAY DEPTH
!
!      (BE CAREFUL WITH THE SIGN OF THE TOPOGR. DATA ! , THIS VERSION
!       USES POSITIVE DEPTH DATA BECAUSE SYNBAPS COMES WITH POSITIVE
!       VALUES)
!
#ifdef bounds_exch_tp
      if ( p_pe == p_io ) then
        allocate(mus(ie_g,je_g))      
        mus=0
      else
        mus => NULL()
      endif



#endif
      SUM=0.
      DO K=1,KE
         SUM=SUM+DZW(K)
      ENDDO
!
      DEPMAX=SUM
      DEPMIN=DZW(1)+DZW(2)+1.
!
      WRITE(IO_STDOUT,*)' --------------------------------------------- &
     &                    --------'
      WRITE(IO_STDOUT,*)' MINIMUM WATER DEPTH IS SET TO ',DEPMIN
      WRITE(IO_STDOUT,*)' MAXIMUM WATER DEPTH IS SET TO ',DEPMAX
      WRITE(IO_STDOUT,*)' --------------------------------------------- &
     &                    --------'
!

      depto_g(:,:)=0.0
!      deuto_g(:,:)=0.0


      IF (ISTART .EQ. 0) THEN
      DO 1296 J=1,JE
      DO 1296 I=1,IE
!
      IF(DEPTO(I,J).LE.1.) DEPTO(I,J)=0.
      IF(DEPTO(I,J).GT.1.AND.DEPTO(I,J).LT.DEPMIN)                      &
     & DEPTO(I,J)=DEPMIN
1296  CONTINUE
      ENDIF
!
      TMAX=0.
      DO J=1,JE
       DO I=1,IE
        TMAX=MAX(TMAX,DEPTO(I,J))
       ENDDO
      ENDDO
      CALL global_max(TMAX)
      WRITE(IO_STDOUT,*)'TOPOGRAPHIE CHECK: ',TMAX
!
      IF (ISTART .EQ. 0) CALL gather_arr(DEPTO,DEPTO_G,p_io)
!
      IF (p_pe==p_io) THEN
         OPEN(IO_IN_GIGI,FILE='topo'                                    &
     &               ,ACCESS='SEQUENTIAL',FORM='FORMATTED')
         DO I1=2,IE_G-1,20
            I2=MIN(I1+19,IE_G-1)
            WRITE(IO_STDOUT,*) 'LESSTREIFEN ',I1,I2
            IF (ISTART .EQ. 0) THEN
            WRITE(IO_IN_GIGI,*)'STREIFEN ',I1,I2
            DO J=1,JE_G
               WRITE(IO_IN_GIGI,63885)J,(NINT(DEPTO_G(I,J)),I=I1,I2)
            ENDDO
            ELSE

            READ(IO_IN_GIGI,*) BLORD
            DO J=1,JE_G
#ifdef bounds_exch_tp 
              READ(IO_IN_GIGI,63886)JJJ,(DEPTO_G(I,J),I=I1,I2)
#else
              READ(IO_IN_GIGI,63885)JJJ,(DEPTO_G(I,J),I=I1,I2)
#endif
            ENDDO
         ENDIF

 63886       FORMAT(I5,20F6.0)    !format for tp04
 63885       FORMAT(I5,20F5.0)         

      ENDDO
         
         CLOSE(IO_IN_GIGI)

#ifdef bounds_exch_tp     
    
!         do i=2,ie_g/2-1
!           ir=i
!           il=ie_g+1-i
!           dpp=0.5*(depto_g(ir,2)+depto_g(il,2))
!           dpp13=0.5*(depto_g(ir,3)+depto_g(il,3))
!           depto_g(ir,2)=dpp13
!           depto_g(il,2)=dpp13
!           depto_g(ir,3)=dpp13
!           depto_g(il,3)=dpp13
!           depto_g(ir,4)=dpp13
!           depto_g(il,4)=dpp13
!           write(0,*)'global',i,ir,il,dpp
!         enddo
         where ( depto_g > 0.5 ) mus=13 
         do i=2,ie_g/2
         ir=ie_g+1-i
         write(0,6637)i,ir,(depto_g(i,j),j=12,3,-1),(depto_g(ir,j),j=3,12)
6637     format(2i4,10f5.0,'|',10f5.0)
         enddo
         do i=2,ie_g/2
         ir=ie_g+1-i
         write(0,6638)i,ir,(mus(i,j),j=62,3,-1),(mus(ir,j),j=3,62)
6638     format(2i3,60i1,'|',60i1)
         enddo    

         deallocate(mus)
#endif
      ENDIF

      IF (ISTART .NE. 0) CALL scatter_arr(DEPTO_G,DEPTO,p_io)

!     FLACHE TEILE DER TOPOGRAPHIE KORRIGIEREN
!

      DO J=1,JE
       DO I=1,IE
        IF(DEPTO(I,J).GT.0.5) THEN
           DEPALT=DEPTO(I,J)
           DEPTO(I,J)=MAX(DEPTO(I,J),DZW(1)+DZW(2))
           depto(i,j)=MIN(depto(i,j),depmax)

#ifndef bounds_exch_tp
           if ( icontro > 0 ) then
              if (depalt.ne.depto(i,j)) then    
                 WRITE(IO_STDOUT,*)'TOPOGRAPHY CHANGED AT: ' &
                     ,I+p_ioff,J+p_joff,DEPALT,'==>',DEPTO(I,J)
              endif
           endif
#endif
        ELSE
           DEPTO(I,J)=0.	   
        ENDIF
       ENDDO
      ENDDO

      DO J=1,JE
       DO I=1,IE
        DO K=1,KE
         IF(ABS(DEPTO(I,J)-TIESTU(K)).LT.0.5) DEPTO(I,J)=TIESTU(K)+1.
        ENDDO
       ENDDO
      ENDDO

      CALL bounds_exch('p',DEPTO,'boden 10')
!
!     Now the calculations on DEPTO are finished, we gather and 
!     distribute global DEPTO_G (which is also needed for WETO_G)
!
      CALL gather_arr(DEPTO,DEPTO_G,0)
      CALL p_bcast(DEPTO_G,0)

      DO 296 J=1,JE1
      DO 296 I=1,IE1
!
      DEUTE(I,J)=MIN(DEPTO(I,J),DEPTO(I,J+1))
      DEUTO(I,J)=MIN(DEPTO(I,J),DEPTO(I+1,J))
296   CONTINUE
      CALL bounds_exch('v+',DEUTE,'boden 11')
      CALL bounds_exch('u+',DEUTO,'boden 12')
!
!----------------------------------------------------------------------
!     COMPUTE BOUNDARY TOPOGRAPHY , I.E. 3 ROWS
!     FROM THE NORTHERN AND FROM THE SOUTHERN BOUNDARY ARE SET
!     TO ZERO DEPTH
!
!     CALCULATION OF DEPTH AT SCALAR POINTS (HP : DEPTH AT PRESSURE PT.)
!        HP IS THE MAXIMUM DEPTH OF THE 4 SURROUNDING U-POINT DEPTHS (H)

      LAUF=0
      JJ=0
      NUM_G(:,:) = 0
      DO III=2,IE_G-1
         JJ=JJ+1
         IF(MOD(III,2).EQ.0)THEN
            I=III/2+1
         ELSE
            I=IE_G-(III/2)
         ENDIF
         DO J=1,JE_G
            IF(DEPTO_G(I,J).GE.1.) THEN
               LAUF=LAUF+1
               NUM_G(I,J)=LAUF
            ENDIF
         ENDDO
      ENDDO
      IF(ICYCLI.GE.1)THEN
         DO J=1,JE_G
            NUM_G(1,J)=NUM_G(IE_G-1,J)
            NUM_G(IE_G,J)=NUM_G(2,J)
         ENDDO
      ENDIF

!     Set local NUM

      NUM(1:ie,1:je) = NUM_G(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)
!
      WRITE(IO_STDOUT,*)' FEUCHTPUNKTE',LAUF
!
!============WERTE FUER TRIANGULARISIERUNG
      MATR=LAUF
      NMAX=0
!
      DO J=2,JE_G-1
      DO I=2,IE_G-1
        IF(NUM_G(I,J).EQ.0 .OR. NUM_G(I-1,J).EQ.0) CYCLE
        NDIF1=NUM_G(I,J)-NUM_G(I-1,J)
        IF(NDIF1.GT.NMAX) NMAX=NDIF1
      ENDDO
      ENDDO
!
      WRITE(IO_STDOUT,*)' MAXDIF ', NMAX
!======LAND/OZEAN-STEUERFELDER UND SCHICHTDICKEN
!
      DO 444 K=1,KE
      DO 444 J=1,JE
         JEV=2*J
      DO 444 I=1,IE
!------------------------------  EVEN LINES 2,4,6,...,2*JE ---------
!
      IF(TIESTU(K+1).LT.DEUTE(I,J)) DDUE(I,J,K)=DZW(K)
      IF(TIESTU(K).LT.DEUTE(I,J) .AND. TIESTU(K+1).GE. DEUTE(I,J))      &
     &    DDUE(I,J,K)=DEUTE(I,J)-TIESTW(K)
!
      IF(DDUE(I,J,K).GT.ZERO) AMSUE(I,J,K)=1.
!====================================================================
!
!------------------------------   ODD LINES 1,3,5,...,2*JE-1 -------
!
      IF(TIESTU(K+1).LT.DEUTO(I,J)) DDUO(I,J,K)=DZW(K)
      IF(TIESTU(K).LT.DEUTO(I,J) .AND. TIESTU(K+1).GE. DEUTO(I,J))      &
     &    DDUO(I,J,K)=DEUTO(I,J)-TIESTW(K)
      IF(DDUO(I,J,K).GT.ZERO) AMSUO(I,J,K)=1.
      IF(TIESTU(K).LT.DEPTO(I,J)) WETO(I,J,K)=1.
      IF(TIESTU(K+1).LT.DEPTO(I,J)) DDPO(I,J,K)=DZW(K)
      IF(TIESTU(K).LT.DEPTO(I,J) .AND. TIESTU(K+1).GE. DEPTO(I,J))      &
     &    DDPO(I,J,K)=DEPTO(I,J)-TIESTW(K)
!
       IF(DDPO(I,J,K) .EQ. ZERO) GOTO 444
       DPIO(I,J,K)=1./DDPO(I,J,K)
  444 CONTINUE

      WETO_G(:,:,:) = 0.
      DO K=1,KE
      DO J=1,JE_G
      DO I=1,IE_G
        IF(TIESTU(K).LT.DEPTO_G(I,J)) WETO_G(I,J,K)=1.
      ENDDO
      ENDDO
      ENDDO


      DDPSIO(:,:,:)=0.


      DO K=1,KE
         DO J=1,JE1
            DO I=1,IE1

               if ( icontro > 0 ) then
               IF(DDUE(I,J,K).GT.0..AND.DDPO(I,J,K)*DDPO(I,J+1,K).LT.1.) THEN
               WRITE(IO_STDOUT,*) ' DDUE ',I+p_ioff,J+p_joff,K   &
                          ,DDUE(I,J,K),DDPO(I,J,K),DDPO(I,J+1,K)
               ENDIF
         
               IF(DDUO(I,J,K).GT.0..AND.DDPO(I,J,K)*DDPO(I+1,J,K).LT.1.) THEN
                     WRITE(IO_STDOUT,*) ' DDUO ',I+p_ioff,J+p_joff,k &
                          ,DDUO(I,J,K),DDPO(I,J,K),DDPO(I+1,J,K)
               ENDIF
               endif

               SUPPO=0.
               SUGGO=0.
               IF(DDPO(I,J,K).GT.1.) THEN
                  SUPPO=SUPPO+DDPO(I,J,K)
                  SUGGO=SUGGO+1.
               ENDIF
               IF(DDPO(I+1,J,K).GT.1.) THEN
                  SUPPO=SUPPO+DDPO(I+1,J,K)
                  SUGGO=SUGGO+1.
               ENDIF
               IF(DDPO(I,J+1,K).GT.1.) THEN
                  SUPPO=SUPPO+DDPO(I,J+1,K)
                  SUGGO=SUGGO+1.
               ENDIF
               IF(DDPO(I+1,J+1,K).GT.1.) THEN
                  SUPPO=SUPPO+DDPO(I+1,J+1,K)
                  SUGGO=SUGGO+1.
               ENDIF
               IF(SUGGO.GE.1.) DDPSIO(I,J,K)=SUPPO/SUGGO

            enddo
         enddo
      enddo


      CALL bounds_exch('s',DDPSIO,'boden 13')
!
      DO 456 J=1,JE
      DO 456 I=1,IE
       KCONDEP(I,J)=NINT(WETO(I,J,1))-99*NINT(1.-WETO(I,J,1))
!
       DEUTIE(I,J)=0.
       DEUTIO(I,J)=0.
!
       IF(DEUTE(I,J).GT.ONE) DEUTIE(I,J)=1./DEUTE(I,J)
       IF(DEUTO(I,J).GT.ONE) DEUTIO(I,J)=1./DEUTO(I,J)
!
  456 CONTINUE
!
!====================================================================
!
!     COMPUTE NUMBER OF WET GRIDPOINTS FOR ODD/EVEN AND SCALAR/VECTOR-
!             POINTS AND STORE ON CNUMWE(KE) : WETE   SCALAR EVEN
!                                 CNUMWO(KE) : WETO   SCALAR ODD
!                                 CNUMUE(KE) : AMSUE  VECTOR EVEN
!                                 CNUMUO(KE) : AMSUO  VECTOR ODD
!             USED FOR CALCULATION OF LAYER MEAN VALUES
!
!             FOR ZONAL INTEGRATION COMPUTE CMERWE(JE,KE)
!                                           CMERWO(..,..)
!                                           ....UE   .
!                                           ....UO   .
!
      DO 19442 K=1,KE
      CWE=0.0
      CWO=0.0
      CUE=0.0
      CUO=0.0
      CALL gather_arr(AMSUE(:,:,k),AMSUE_G,0)
      CALL p_bcast(AMSUE_G,0)
      CALL gather_arr(AMSUO(:,:,k),AMSUO_G,0)
      CALL p_bcast(AMSUO_G,0)
      DO 19443 J=1,JE_G
      CMWE=0.0
      CMWO=0.0
      CMUE=0.0
      CMUO=0.0
      DO 19444 I=1,IE_G
!
      IF(WETO_G(I,J,K).GE.0.5) CMWO=CMWO+1.
      IF(AMSUE_G(I,J).GE.0.5)CMUE=CMUE+1.
      IF(AMSUO_G(I,J).GE.0.5)CMUO=CMUO+1.
19444 CONTINUE
      JJ = J-p_joff
      IF(JJ>=1 .AND. JJ<=JE) THEN
        CMERWE(JJ,K)=CMWE
        IF(CMWE.LT.0.1)CMERWE(JJ,K)=1.
        CMERWO(JJ,K)=CMWO
        IF(CMWO.LT.0.1)CMERWO(JJ,K)=1.
        CMERUE(JJ,K)=CMUE
        IF(CMUE.LT.0.1)CMERUE(JJ,K)=1.
        CMERUO(JJ,K)=CMUO
        IF(CMUO.LT.0.1)CMERUO(JJ,K)=1.
      ENDIF
      CWE=CWE+CMWE
      CWO=CWO+CMWO
      CUE=CUE+CMUE
      CUO=CUO+CMUO
19443 CONTINUE
      CNUMWE(K)=CWE
      CNUMWO(K)=CWO
      CNUMUE(K)=CUE
      CNUMUO(K)=CUO
19442 CONTINUE
      WRITE(IO_STDOUT,19662)
19662 FORMAT(' WET POINTS IN LAYERS 1-13 FOR SCALAR EVEN/ODD AND VECTOR &
     &EVEN/ODD (TOTAL NUMBER) ')
      WRITE(IO_STDOUT,19663)(CNUMWE(K),K=1,KE)
      WRITE(IO_STDOUT,19663)(CNUMWO(K),K=1,KE)
      WRITE(IO_STDOUT,19663)(CNUMUE(K),K=1,KE)
      WRITE(IO_STDOUT,19663)(CNUMUO(K),K=1,KE)
19663 FORMAT(' ',20F6.0)


      CALL bounds_exch('v+',DEUTIE,'boden 14')
      CALL bounds_exch('u+',DEUTIO,'boden 15')
      CALL bounds_exch('v+',DEUTE,'boden 16')
      CALL bounds_exch('u+',DEUTO,'boden 17')
      CALL bounds_exch('u+',DDUO,'boden 18')
      CALL bounds_exch('v+',DDUE,'boden 19')
      CALL bounds_exch('v+',AMSUE,'boden 20')
      CALL bounds_exch('u+',AMSUO,'boden 21')
!#ifdef bounds_exch_save
      CALL bounds_exch('p',WETO,'boden 22')
      CALL bounds_exch('p',DDPO,'boden 23')
!#endif
!
!
!----------------------------------------------------------------------
!  PRINT DEPTHS AT VECTOR CELLS
!
!      CALL gather_arr(DEUTO,DEUTO_G,0)
!      CALL p_bcast(DEUTO_G,0)
!!$      DO I=1,IE_G
!!$         DO J=1,JE_G
!!$            CDEUTO(I,J)='7'
!!$            IF(DEUTO_G(I,J).LE.6000.) CDEUTO(I,J)='6'
!!$            IF(DEUTO_G(I,J).LE.5000.) CDEUTO(I,J)='5'
!!$            IF(DEUTO_G(I,J).LE.4000.) CDEUTO(I,J)='4'
!!$            IF(DEUTO_G(I,J).LE.3000.) CDEUTO(I,J)='3'
!!$            IF(DEUTO_G(I,J).LE.2000.) CDEUTO(I,J)='2'
!!$            IF(DEUTO_G(I,J).LE.1000.) CDEUTO(I,J)='1'
!!$            IF(DEUTO_G(I,J).EQ.   0.) CDEUTO(I,J)='*'
!!$         ENDDO
!!$      ENDDO
!!$      IPRANZ=(2*IE_G-1)/120 + 1
!!$      DO IPR=1,IPRANZ
!!$       IANF =    (IPR-1)*60    + 1
!!$       IEND =     IPR   *60
!!$       WRITE(IO_STDOUT,*) '    1M - 1000M : 1'
!!$       WRITE(IO_STDOUT,*) ' 1001M - 2000M : 2'
!!$       WRITE(IO_STDOUT,*) ' 2001M - 3000M : 3'
!!$       WRITE(IO_STDOUT,*) ' 3001M - 4000M : 4'
!!$       WRITE(IO_STDOUT,*) ' 4001M - 5000M : 5'
!!$       WRITE(IO_STDOUT,6337) (I,I=IANF+1,IEND,2) 
!!$       DO J=1,JE_G
!!$        IF(MOD(J,10).NE.0) THEN
!!$         WRITE(IO_STDOUT,6334) J,(CDEUTO(MOD(I-1,IE_G)+1,J),I=IANF,IEND)
!!$        ELSE
!!$         WRITE(IO_STDOUT,6335) J,(CDEUTO(MOD(I-1,IE_G)+1,J),I=IANF,IEND)
!!$        ENDIF
!!$       ENDDO
!!$      ENDDO
!!$!
!!$ 6334 FORMAT(6X,I3,3X,6(9(A1,1X),(A1,".")))
!!$ 6335 FORMAT(6X,I3,3X,6(9(A1,"-"),(A1,".")))
!!$ 6337 FORMAT(6X,'DEPTHS AT VECTOR CELLS:',/,11X,30I4,/)
!!$!
!=====================================================================
      WRITE(IO_STDOUT,*)'SR BODEN UEBERLEBT!!!'
!#ifdef bounds_exch_save
      CALL bounds_exch('p',DLXP,'boden 24')
      CALL bounds_exch('p',DLYP,'boden 25')
!#endif
!
      RETURN
      END
