SUBROUTINE OCVTOT

  USE MO_PARAM1
  USE MO_MPI
  USE MO_PARALLEL
  USE MO_COMMO1
  USE MO_UNITS
  
  IMPLICIT NONE

!  REAL :: psiz(kep)
  REAL :: surdif(ie,je)


  INTEGER :: i, j, k

  REAL :: vake, uako, uiny,wel,zel,bel

!$OMP PARALLEL PRIVATE(i,j,k,vake,uako)

!$OMP DO
      DO  K=1,KE
        DO  J=1,JE
          DO  I=1,IE
            T1O(I,J,K)=UOO(I,J,K)
            S1O(I,J,K)=VOE(I,J,K)
            VKE(I,J,K)=AMSUE(I,J,K)*(VSE(I,J)+VKE(I,J,K))
            UKO(I,J,K)=AMSUO(I,J,K)*(USO(I,J)+UKO(I,J,K))
!            VKE(I,J,K)=AMSUE(I,J,K)*VSE(I,J)
!            UKO(I,J,K)=AMSUO(I,J,K)*USO(I,J)
          ENDDO
        ENDDO
      ENDDO
!      do i=17,26
!      WRITE(IO_STDOUT,*)'trovel',i,uko(i,2,1),uko(i,2,1)+uko(83-i,2,1),&
!     & vke(i,1,1),vke(i,1,1)+  vke(84-i,2,1)
!      enddo
!
!RJ 6991  FORMAT(6E20.12)
!RJ       EPSM=1.E-20
!RJ       USUA=0.
!
!$OMP DO
      DO J=1,JE
        DO I=1,IE
          WO(I,J,KEP) = ZERO
        ENDDO
      ENDDO
!
!RJ       USUM=0.
!$OMP DO
      DO K=1,KE
        DO J=1,JE
          DO I=1,IE
            VAKE=VKE(I,J,K)
            UAKO=UKO(I,J,K)
            VKE(I,J,K)=CONN*VKE(I,J,K)+VOE(I,J,K)*CONO
            UKO(I,J,K)=CONN*UKO(I,J,K)+UOO(I,J,K)*CONO
            UOO(I,J,K)=UAKO
            VOE(I,J,K)=VAKE
          ENDDO
        ENDDO
      ENDDO
!
!
!======================================================================
!
!     B)
!
!     VERTICAL VELOCITY = VERTICAL INTEGRAL OF DIVERGENCE OF
!                             HORIZONTAL VELOCITY FIELD
!
!$OMP DO
      DO J=1,JE
        DO I=1,IE
          WO(I,J,KEP) = ZERO
        ENDDO
      ENDDO
!
!
!
!RJ       UDIMA=0.
!RJ       K=1
!RJ       UPTRA(K)=0.
!RJ       DOTRA(K)=0.

          UPTRT(:)=0.

          DOTRT(:)=0.
!
!
!$OMP DO
      DO K=KEp,1,-1
        DO J=1,JE
          DO I=1,IE
            WO(I,J,K) = ZERO
          ENDDO
        ENDDO
      ENDDO
!
!$OMP DO
      DO J=2,JE1
        DO K=KE,1,-1
!
          DO I=2,IE1
!
            WO(I,J,K) = WO(I,J,K+1)+  WETO(I,J,K) * (                 &
                  UKO(I-1,J,K) * DDUO(I-1,J,K)*DLYU(I-1,J) &
                 - UKO(I,J,K)   * DDUO(I,J,K)*DLYU(I,J)              &
        +      VKE(I,J,K)     * DDUE(I,J,K)*DLXV(I,J)                  &
        &-    VKE(I,J-1,K)   * DDUE(I,J-1,K)*DLXV(I,J-1))/(dlyp(i,j)*DLXP(I,J))
!
            if (icontro.ne.0) then
               UPTRT(K)=UPTRT(K)+DLXP(I,J)*DLYP(I,J)*(WO(I,J,K)+ABS(WO(I,J,K)))
               DOTRT(K)=DOTRT(K)-DLXP(I,J)*DLYP(I,J)*(WO(I,J,K)-ABS(WO(I,J,K)))
            endif

         ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL

      if (icontro.ne.0) then
         uptrt=uptrt*1.e-6
         dotrt=dotrt*1.e-6

         call global_sum(uptrt)
         call global_sum(dotrt)
         if (p_pe == p_io) then
         write(0,*)'TOTAL UPWELLING [sv]'
        do k=1,ke
           write(0,621)k,uptrt(k),dotrt(k)
        enddo
621      format(i2,2f12.2)
      endif
      endif

      CALL bounds_exch('p',WO,'ocvtot 1')

!!$      IF(NNNDT.EQ.360) THEN
!!$        DO  I=2,IE-1
!!$          PSIZ(KEP)=0.
!!$          DO K=KE,1,-1
!!$            UINY=0.
!!$            DO J=2,JE-1
!!$              UINY=UINY+UKO(I,J,K)*DDUO(I,J,K)*DLYU(I,J)
!!$            ENDDO
!!$            PSIZ(K)=PSIZ(K+1)+1.E-6*UINY
!!$          ENDDO
!!$
!!$!     RJ: Don't know if this is needed for something
!!$!         If it is needed it should be summed up over rows (not over all!!!!)
!!$          WRITE(IO_STDOUT,*) 'Sum on local PE:'
!!$          WRITE(IO_STDOUT,6636) PSIZ
!!$        ENDDO
!!$6636    FORMAT(21F6.0)
!!$      ENDIF

      if (icontro.ne.0) then      
      surdif(:,:)=wo(:,:,1)*dt-z1o(:,:)

      WRITE(0,*) 'surdif check (values < 1.0e-12): ', &
           MINVAL(surdif), MAXVAL(surdif)
      endif

!!$      do i=17,26
!!$      WRITE(IO_STDOUT,723)(surdif(i,j),j=1,6),(z1o(i,j),j=1,6)
!!$      enddo
!!$      wel=0.
!!$      zel=0.
!!$      bel=0.
!!$      j=2
!!$      do i=2,42
!!$      wel=wel+wo(i,j,1)*dt*weto(i,j,1)*dlxp(i,j)*dlyp(i,j)
!!$      zel=zel+z1o(i,j)*weto(i,j,1)*dlxp(i,j)*dlyp(i,j)
!!$      bel=bel+b1o(i,j)*weto(i,j,1)
!!$      enddo
!!$!      j=3
!!$!      do i=2,ie1
!!$!      wel=wel+wo(i,j,1)*dt*weto(i,j,1)*dlxp(i,j)*dlyp(i,j)
!!$!      zel=zel+z1o(i,j)*weto(i,j,1)*dlxp(i,j)*dlyp(i,j)
!!$!      bel=bel+b1o(i,j)*weto(i,j,1)
!!$!      enddo
!!$!      write(6,*)'flaeche',wel,zel,bel
!!$!      wel=0.
!!$!      zel=0.
!!$!      bel=0.
!!$      do j=3,je
!!$      do i=2,ie1
!!$      wel=wel+wo(i,j,1)*dt*weto(i,j,1)*dlxp(i,j)*dlyp(i,j)
!!$      zel=zel+z1o(i,j)*weto(i,j,1)*dlxp(i,j)*dlyp(i,j)
!!$      bel=bel+b1o(i,j)*weto(i,j,1)
!!$      enddo
!!$      enddo
!!$
!!$      write(6,*)'flaeche',wel,zel,bel
!!$723   format(12e10.3)      

    END SUBROUTINE OCVTOT
