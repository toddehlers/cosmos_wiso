      SUBROUTINE BELEG_ZERO
!
!     BELEG_ZERO : INITIALIZES ALL FIELDS DEFINED in MO_COMMO1 by ZERO  
!---------------------------------------------------------------
      USE mo_param1
      USE mo_commo1
      USE mo_commoau2
      USE mo_commoau3
      USE mo_commobbl


!#ifdef MEAN         
!      USE MO_MEAN
!#endif

      USE mo_elicom

      USE mo_para2

!:: COMMON BLOCK GRID     
!:: THREE- DIM FIELDS
      DO K=1,KE
       DO J=1,JE
        DO I=1,IE
         DDUE(I,J,K)=ZERO
         DDUO(I,J,K)=ZERO
         DDPO(I,J,K)=ZERO
         DPIO(I,J,K)=ZERO
         DDPSIO(I,J,K)=ZERO
        ENDDO
       ENDDO
      ENDDO
    

!:: TWO- DIM FIELDS
      DO J=1,JE
       DO I=1,IE
        DLXP(I,J)=ZERO
        DLYP(I,J)=ZERO
        DLXU(I,J)=ZERO
        DLYV(I,J)=ZERO
        DLXV(I,J)=ZERO
        DLYU(I,J)=ZERO
        DTDXUO(I,J)=ZERO
        DTDXUE(I,J)=ZERO
        DTDYO(I,J)=ZERO
        DTDXPO(I,J)=ZERO
        DPYO(I,J)=ZERO
        DPYE(I,J)=ZERO
        DTDXPE(I,J)=ZERO
        DEUTE(I,J)=ZERO
        DEUTO(I,J)=ZERO
        DEUTIE(I,J)=ZERO
        DEUTIO(I,J)=ZERO
        DEPTO(I,J)=ZERO
       ENDDO
      ENDDO


!:: ONE- DIM FIELDS
      DO K=1,KE     
       DWI(K)=ZERO
       DZW(K)=ZERO
      ENDDO
      DO K=1,KEP    
       DZ(K)=ZERO
       DI(K)=ZERO
       TIESTU(K)=ZERO
       TIESTW(K)=ZERO
      ENDDO


!:: COMMON BLOCK GEO
!:: TWO-DIMENSIONAL FIELDS
      DO J=1,JE
      DO I=1,IE
       ALAT(I,J)=ZERO
      ENDDO
      ENDDO



      GILA(:,:)=ZERO
      GIPH(:,:)=ZERO



!:: COMMON BLOCK MASK
      DO K=1,KE
      DO J=1,JE
       DO I=1,IE
        WETO(I,J,K)=ZERO
        AMSUE(I,J,K)=ZERO
        AMSUO(I,J,K)=ZERO
       ENDDO
      ENDDO
      ENDDO

!:: COMMON BLOCK FIELDS3D
      DO K=1,KE
      DO J=1,JE
       DO I=1,IE
        RHOO(I,J,K)=ZERO
        UKO(I,J,K)=ZERO
        VKE(I,J,K)=ZERO
        THO(I,J,K)=ZERO
        SAO(I,J,K)=ZERO
        PO(I,J,K)=ZERO
        VK1E(I,J,K)=ZERO
        UK1O(I,J,K)=ZERO
        VOE(I,J,K)=ZERO
        UOO(I,J,K)=ZERO
        T1O(I,J,K)=ZERO
        S1O(I,J,K)=ZERO
        STABIO(I,J,K)=ZERO
       ENDDO
      ENDDO
      ENDDO

      DO K=1,KEP
      DO J=1,JE
       DO I=1,IE
        DVO(I,J,K)=ZERO
        AVO(I,J,K)=ZERO
        WO(I,J,K)=ZERO
       ENDDO
      ENDDO
      ENDDO

!:: COMMON BLOCK FIELDS2D
      DO J=1,JE
       DO I=1,IE
       ZO(I,J)=ZERO
       Z1E(I,J)=ZERO
       Z1O(I,J)=ZERO
       U1E(I,J)=ZERO
       V1E(I,J)=ZERO
       U1O(I,J)=ZERO
       V1O(I,J)=ZERO
       UROT(I,J)=ZERO
       UDIV(I,J)=ZERO
       USNO(I,J)=ZERO
       VSNE(I,J)=ZERO
       USO(I,J)=ZERO
       VSE(I,J)=ZERO
       UCOS(I,J)=ZERO
       VCOS(I,J)=ZERO
       SHELP(I,J)=ZERO
       THELP(I,J)=ZERO
       RHELP(I,J)=ZERO
       TLOW(I,J)=ZERO
       SLOW(I,J)=ZERO
       TSUP(I,J)=ZERO
       SSUP(I,J)=ZERO
       PXOIN(I,J)=ZERO
       PYEIN(I,J)=ZERO
       UZO(I,J)=ZERO
       VZE(I,J)=ZERO
       FTWOU(I,J)=ZERO
       FTWOV(I,J)=ZERO
       CURVAV(I,J)=ZERO
       B1O(I,J)=ZERO
       B1E(I,J)=ZERO
       NUM(I,J)=0   
       BRINE(I,J)=ZERO
       KCONDEP(I,J)=0   
      ENDDO
      ENDDO
!:: COMMON BLOCK FIELDS1D
      DO K=1,KE
      DOTRA(K)=ZERO
      DOTRT(K)=ZERO
      UPTRA(K)=ZERO
      UPTRT(K)=ZERO
      PREFF(K)=ZERO
      ROREF(K)=ZERO
      ENDDO
      DO K=1,KEP
      WINTUR(K)=ZERO
      DBACKV(K)=ZERO
      ABACKV(K)=ZERO
      ENDDO
      DO i=1,ill
      SKAL(I)=ZERO
      ENDDO
      DO i=1,ilt
      B(I)=ZERO
      X(I)=ZERO
      ENDDO
!:: COMMON BLOCK HOPEICF
      DO J=1,JE
       DO I=1,IE
       TICE(I,J)=ZERO
       SICTHO(I,J)=ZERO
       SICOMO(I,J)=ZERO
       SICDIO(I,J)=ZERO
       SICSHO(I,J)=ZERO
       SICPO(I,J)=ZERO
       SICUO(I,J)=ZERO
       SICVE(I,J)=ZERO
       HIBDELO(I,J)=ZERO
       HIBZETO(I,J)=ZERO
       HIBETO(I,J)=ZERO
       HIBDELE(I,J)=ZERO
       HIBZETE(I,J)=ZERO
       HIBETE(I,J)=ZERO
       TAUWATU(I,J)=ZERO
       TAUWATV(I,J)=ZERO
       SICSNO(I,J)=ZERO
       FWO(I,J)=ZERO
       SHO(I,J)=ZERO
       PSIUWE(I,J)=ZERO
      ENDDO
      ENDDO
!:: COMMON BLOCK DIAGCOMN
      DO J=1,JE
       DO I=1,IE
       HFLM(I,J)=ZERO
       PMEM(I,J)=ZERO
       EISTRX(I,J)=ZERO
       EISTRY(I,J)=ZERO
       DO MM=1,12
        AMLD(I,J,MM) = ZERO
       ENDDO
       ENDDO
      ENDDO
!::
      DO J=1,JE
       DO I=1,IE
        KBOT(I,J)=0   
        TXO(I,J)=ZERO
        TYE(I,J)=ZERO
        TAFO(I,J)=ZERO
        HEATO(I,J)=ZERO
        EMINPO(I,J)=ZERO
        FCLOU(I,J)=ZERO
        FSWR(I,J)=ZERO
        FU10(I,J)=ZERO
        FPREC(I,J)=ZERO
        FTDEW(I,J)=ZERO
        RELSAO(I,J)=ZERO
        RELTHO(I,J)=ZERO
       ENDDO
       ENDDO
#ifndef RIVER_GIRIV
       DO NR=1,NUMRIV
        FRIV(NR)=ZERO
        RIVLON(NR)=ZERO
        RIVLAT(NR)=ZERO
        DDRIV(NR)=ZERO
        IRIVI(NR)=0   
        IRIVJ(NR)=0    
        DO NM=1,12
         RIVAL(NR,NM)=ZERO
        ENDDO
       ENDDO
#else
       DO I=1,IE
       DO j=1,JE
         GIRIV(i,j)=ZERO
       ENDDO
       ENDDO
#endif /*RIVER_GIRIV*/

!:: COMMON BLOCK COM11
       DO K=1,KE
        CNUMWE(K)=ZERO
        CNUMWO(K)=ZERO
        CNUMUE(K)=ZERO
        CNUMUO(K)=ZERO
        DO J=1,JE
         CMERWE(J,K)=ZERO
         CMERWO(J,K)=ZERO
         CMERUE(J,K)=ZERO
         CMERUO(J,K)=ZERO
        ENDDO
       ENDDO

!:: USE FILE MO_COMMOAU2
!:: COMMON bLOCK SIOBOD
       DO J=1,JE
        DO I=1,IE
         ACLO(I,J)=ZERO
         PAO(I,J)=ZERO
         RPRECO(I,J)=ZERO
         FRSE(I,J)=ZERO
         PRECO(I,J)=ZERO
         PRECH(I,J)=ZERO
         QSWO(I,J)=ZERO
         QLWO(I,J)=ZERO
         QSEO(I,J)=ZERO
         QLAO(I,J)=ZERO
         TAIRO(I,J)=ZERO
         TDO(I,J)=ZERO
         SICUDO(I,J)=ZERO
         SICVDE(I,J)=ZERO
        ENDDO
       ENDDO
!:: INCLUDE FILE PARA2.h
!:: COMMON BLOCK ITHILF
       DO J=1,JE
        DO I=1,IE
         UF(I,J)=ZERO
         VF(I,J)=ZERO
         FF(I,J)=ZERO
         XX(I,J)=ZERO
        ENDDO
       ENDDO

!:: INCLUDE FILE MO_ELICOM.h
!:: COMMON BLOCK DUMDUM
       DO J=1,ILL
        DO I =1,IMM
         PGL(I,J)=ZERO
        ENDDO
       ENDDO

      END
