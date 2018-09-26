SUBROUTINE TRIAN
  !**********************************************************
  ! UNTERPROGRAMM ZUR TRIANGULARISIERUNG DER MATRIX
  ! LOESUNG LINEARER GL., GAUSS-SCHER ALGOR.
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  USE MO_PARAM1
  USE MO_MPI
  USE MO_PARALLEL
  USE MO_ELICOM
  USE MO_COMMO1
  USE MO_UNITS

  IMPLICIT NONE

  INTEGER I, IN, IO, IS, IW, J, JA, JJ, JO, K, K1, KITER, L, LL, lmax
  REAL FAKE, FAKN, FAKS, FAKW, RE, RR1, WILKIN, ZZZ
  REAL FW(IE,JE), FE(IE,JE), FN(IE,JE), FS(IE,JE)
  REAL FW_G(IE_G,JE_G), FE_G(IE_G,JE_G)
  REAL FN_G(IE_G,JE_G), FS_G(IE_G,JE_G)

  
  DO j = 1, ill
    DO i = 1, kbm
      pgl(i,j) = 0.0
    ENDDO
  ENDDO

  zzz = -g*conn*stabn*dt*dt
  DO j = 2, je-1
    DO i = 2, ie-1
      fakw = deuto(i-1,j)*dlyu(i-1,j)/dlxu(i-1,j)
      fake = deuto(i,j)  *dlyu(i,j)  /dlxu(i,j)
      fakn = deute(i,j-1)*dlxv(i,j-1)/dlyv(i,j-1)
      faks = deute(i,j)  *dlxv(i,j)  /dlyv(i,j)

      fw(i,j) = fakw*zzz
      fe(i,j) = fake*zzz
      fn(i,j) = fakn*zzz
      fs(i,j) = faks*zzz
    ENDDO
  ENDDO

  CALL gather_arr(FW,FW_G,p_io)
  CALL gather_arr(FE,FE_G,p_io)
  CALL gather_arr(FN,FN_G,p_io)
  CALL gather_arr(FS,FS_G,p_io)

!open(1,file="stencil.dat",access='sequential',form='unformatted')
!write (1) fw
!write (1) fe
!write (1) fn
!write (1) fs
!close (1)

  if (p_pe == p_io) then
    lmax = 0
    do j = 2, je_g-1
      do i = 2, ie_g-1
        if (num_g(i,j) == 0) cycle
        l  = num_g(i,j)
        in = num_g(i,j-1)-l
        is = num_g(i,j+1)-l
        iw = num_g(i-1,j)-l
        io = num_g(i+1,j)-l
        if (num_g(i,j-1) <= 0) in = 0
        if (num_g(i,j+1) <= 0) is = 0
        if (num_g(i-1,j) <= 0) iw = 0
        if (num_g(i+1,j) <= 0) io = 0
        pgl(km+iw,l) = fw_g(i,j)
        pgl(km+io,l) = fe_g(i,j)
        pgl(km+in,l) = fn_g(i,j)
        pgl(km+is,l) = fs_g(i,j)
        pgl(km,l)    = dlxp_g(i,j)*dlyp_g(i,j) &
                      -fw_g(i,j)-fe_g(i,j)-fn_g(i,j)-fs_g(i,j)
        lmax = MAX(lmax,l)
      enddo
    enddo


!
!open(1,file='band_matrix.dat',form='unformatted',access='sequential')
!write (1) SIZE(pgl,1), SIZE(pgl,2)
!write (1) pgl
!write (1) ie_g, je_g
!write (1) num_g
!close(1)
!     
    DO L = 1, lmax
      
      SKAL(L)=1./PGL(KM,L)
      DO LL=1,KBM
        PGL(LL,L)=PGL(LL,L)*SKAL(L)
      ENDDO
    ENDDO

    WILKIN=1.
    WRITE(IO_STDOUT,*)' LAUF ',L,MATR
!
    DO KITER=2,MATR
      K1=KITER-1
      RE=PGL(KM,K1)
      JA=KM
      DO K=1,KB
        PGL(K,K1)=0.
        JA=JA-1
        IF(K+K1.GT.MATR) CYCLE
        JO=JA+KB
        RR1=PGL(JA,K+K1)/RE
        PGL(K,K1)=RR1
        IF(ABS(RR1).GT.WILKIN) WILKIN=ABS(RR1)
        DO JJ=JA,JO
          PGL(JJ,K+K1)=PGL(JJ,K+K1)-RR1*PGL(JJ+K,K1)
        ENDDO
      ENDDO
    ENDDO
602 FORMAT(' MAX.ELIMI',E14.5)
    WRITE(IO_STDOUT,602) WILKIN
  ENDIF
      
  CALL p_bcast(PGL,p_io)
  !
END SUBROUTINE TRIAN
