
SUBROUTINE OCCLIT

  USE MO_PARAM1
  USE MO_PARALLEL
  USE MO_COMMO1
  USE MO_UNITS
  
  IMPLICIT NONE
  
  INTEGER :: i,j,k,iter,itermax, jb
  REAL :: brems, contra, ove, speed, under, vorw, uko1, vke1
  !
  REAL :: bruva(ie,je),contrij(ie,je),contrij_g(ie_g,je_g)
  !
  REAL :: wpo(ie,je,ke) ! changed from kep to ke by LK
  !
  !  ITERATION OF BAROCLINIC SYSTEM DIRECTLY IN BAROCLINIC VELOCITIES
  !
!$OMP PARALLEL PRIVATE(i,j,k,iter,itermax,under,uko1,vke1,ove,speed,brems,vorw)
  !
!$OMP DO
  DO J=1,JE
    DO I=1,IE
      PXOIN(I,J)    = 0.
      BRUVA(I,J)    = 0.
      PYEIN(I,J)    = 0.
      STABIO(I,J,1) = 0.
      SLOW(I,J)     = 0.
    ENDDO
  ENDDO
  !
!$OMP DO
  DO J=2,JE1
    DO K=1,KE
      DO I=2,IE1
        PXOIN(I,J)=PXOIN(I,J)+DDUO(I,J,K)*DTDXUO(I,J)*(PO(I+1,J,K)-PO(I,J,K))
        PYEIN(I,J)=PYEIN(I,J)+DDUE(I,J,K)*DPYE(I,J)*(PO(I,J,K)-PO(I,J+1,K))
        BRUVA(I,J)=BRUVA(I,J)+TIESTW(K+1)*STABIO(I,J,K)
      ENDDO
    ENDDO
  ENDDO
  !
!$OMP DO
  DO J=2,JE1
    DO I=2,IE1
      PXOIN(I,J)=PXOIN(I,J)*DEUTIO(I,J)
      PYEIN(I,J)=PYEIN(I,J)*DEUTIE(I,J)
    ENDDO
  ENDDO
  !
#ifdef bounds_exch_tp
   call vcheck(11,vke)
   call vcheck(12,vk1e)
#endif

!$OMP SINGLE
  CALL bounds_exch('u',PXOIN,'occlit 1')
  call bounds_exch('v',PYEIN,'occlit 2')
  call bounds_exch('p',bruva,'occlit 3')
!$OMP END SINGLE

!$OMP DO
  DO J=1,JE
    DO K=1,KE
      DO I=1,IE
        UKO(I,J,K)=UKO(I,J,K)*DDUO(I,J,K)
        VKE(I,J,K)=VKE(I,J,K)*DDUE(I,J,K)
        UK1O(I,J,K)=UK1O(I,J,K)*DDUO(I,J,K)
        VK1E(I,J,K)=VK1E(I,J,K)*DDUE(I,J,K)
      ENDDO
    ENDDO
  ENDDO
  !

#ifdef bounds_exch_tp
   call vcheck(13,vke)
   call vcheck(14,vk1e)
#endif


!$OMP DO
  DO J=2,JE1
    DO K=KE,1,-1  
      DO I=2,IE1
        WO(I,J,K) = WO(I,J,K+1) + DTI * WETO(I,J,K) * (                   &
               DTDXPO(I,J) * (   UKO(I-1,J,K)*DLYU(I-1,J)                 &
                               - UKO(I,J,K)*DLYU(I,J)    )/DLYP(I,J)      &
             + DTDYO(I,J)  * (   VKE(I,J,K)*DLXV(I,J)                     &
                               - VKE(I,J-1,K)*DLXV(I,J-1))/DLXP(I,J)  )
      ENDDO
    ENDDO
  ENDDO
  !
!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('p',WO,'occlit 4')
!$OMP END SINGLE
!#endif
  !
!$OMP DO
  DO J=1,JE
    DO K=1,KE
      DO I=1,IE
        T1O(I,J,K)=WO(I,J,K)
        SLOW(I,J)=SLOW(I,J)+STABIO(I,J,K)*DT*WO(I,J,K)
        WPO(I,J,K)=SLOW(I,J)
      ENDDO
    ENDDO
  ENDDO
  !
  !OtB Needs more iterations than 6
  !     ITERMAX=6
  !     ITERMAX=16
  !Uwe says 8 is cheaper
  ITERMAX=8 

#ifdef bounds_exch_tp
  ITERMAX=12
#endif
  !
  iteration_loop: DO ITER=1,ITERMAX
    !
!$OMP DO
    DO J=1,JE
      DO I=1,IE
        UCOS(I,J)  = 0.
        VCOS(I,J)  = 0.
        PXOIN(I,J) = 0.
        PYEIN(I,J) = 0.
      ENDDO
    ENDDO
    !
!$OMP DO
    DO J=2,JE1
      DO K=1,KE
        DO I=2,IE1
          PXOIN(I,J)=PXOIN(I,J)+DDUO(I,J,K)*DTDXUO(I,J)*                   &
               (WPO(I+1,J,K)-WPO(I,J,K))*AMSUO(I,J,K)
          PYEIN(I,J)=PYEIN(I,J)+DDUE(I,J,K)*DPYE(I,J)*                     &
               (WPO(I,J,K)-WPO(I,J+1,K))*AMSUE(I,J,K)
        ENDDO
      ENDDO
    ENDDO
    !
!$OMP DO
    DO J=2,JE1
      DO I=2,IE1
        PXOIN(I,J)=PXOIN(I,J)*DEUTIO(I,J)
        PYEIN(I,J)=PYEIN(I,J)*DEUTIE(I,J)
      ENDDO
    ENDDO
    !
!$OMP SINGLE
  CALL bounds_exch('u',PXOIN,'occlit 5')
  CALL bounds_exch('v',PYEIN,'occlit 6')
!$OMP END SINGLE
  !
!$OMP DO
    DO J=2,JE1
      DO K=1,KE
        DO I=2,IE1
          uko1       = uk1o(i,j,k)+stabn*(pxoin(i,j)                    &
                      -dtdxuo(i,j)*(wpo(i+1,j,k)-wpo(i,j,k)))*dduo(i,j,k)
          vke1       = vk1e(i,j,k)+stabn*(pyein(i,j)                    &
                      -dpye(i,j)*(wpo(i,j,k)-wpo(i,j+1,k)))*ddue(i,j,k)
          uko(i,j,k) = uko1*amsuo(i,j,k)
          vke(i,j,k) = vke1*amsue(i,j,k)

          ucos(i,j)  = ucos(i,j)+uko(i,j,k)
          vcos(i,j)  = vcos(i,j)+vke(i,j,k)

        ENDDO
      ENDDO
    ENDDO
    !
!$OMP DO
    DO J=2,JE1
      DO K=1,KE
        DO I=2,IE1
          UKO(I,J,K)=UKO(I,J,K)-UCOS(I,J)*DEUTIO(I,J)*DDUO(I,J,K)
          VKE(I,J,K)=VKE(I,J,K)-VCOS(I,J)*DEUTIE(I,J)*DDUE(I,J,K)
        ENDDO
      ENDDO
    ENDDO
    !
!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch('u',UKO,'occlit 7')
    CALL bounds_exch('v',VKE,'occlit 8')
!$OMP END SINGLE
!#endif
    !
    IF (ITER == ITERMAX) CYCLE
    !
!$OMP DO
    DO J=1,JE
      DO I=1,IE
        TLOW(I,J)    = 0.
        SLOW(I,J)    = 0.
        CONTRIJ(I,J) = 0
      ENDDO
    ENDDO
    !
!$OMP DO
    DO J=2,JE1
      DO K=KE,1,-1
        DO I=2,IE1
          !
!          T1O(I,J,K)= SLOW(I,J) + DTI * WETO(I,J,K) * (                    &
!               DTDXPO(I,J) * ( UKO(I-1,J,K)*DLYU(I-1,J)                    &
!                              -UKO(I,J,K)  *DLYU(I,J)  )/DLYP(I,J)         &
!             + DTDYO(I,J)  * ( VKE(I,J,K)  *DLXV(I,J)                      &
!                              -VKE(I,J-1,K)*DLXV(I,J-1))/DLXP(I,J)  )

!emr new formulation 
          T1O(I,J,K)= SLOW(I,J) + weto(i,j,k)* (                    &
       &          UKO(I-1,J,K)*DLYU(I-1,J)                    &
                 -UKO(I,J,K)  *DLYU(I,J)         &
             +   VKE(I,J,K)  *DLXV(I,J)                      &
       &      -VKE(I,J-1,K)*DLXV(I,J-1))/(dlyp(i,j)*DLXP(I,J)  )


          SLOW(I,J)=T1O(I,J,K)
        ENDDO
      ENDDO
    ENDDO
    !
!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch('p',T1O,'occlit 9')
!$OMP END SINGLE
!#endif
    !
    IF (iter == 1) THEN
      under = 0.0
    ELSE
      under = 1.0
    ENDIF

      contrij_g=0.

    jb=2
#ifdef bounds_exch_tp
     if(p_joff.eq.0) jb=3
#endif
    !
!$OMP DO
    DO J=jb,JE1
      DO K=2,KE
        DO I=2,IE1
          OVE=WPO(I,J,K)
          TLOW(I,J)=TLOW(I,J)+G*STABIO(I,J,K)*DT*(CONN*T1O(I,J,K)+WO(I,J,K))
!          SPEED=2.*G*BRUVA(I,J)*(DTDXUO(I,J)**2+DTDYO(I,J)**2)
!          BREMS=UNDER*SPEED*CONN*STABN/(1.+SPEED)
          SPEED=G*BRUVA(I,J)*(DTDXUO(I,J)**2+DTDYO(I,J)**2)
          BREMS=2.*UNDER*SPEED*CONN*STABN/(1.+SPEED)
          VORW=1.-BREMS
          WPO(I,J,K)=VORW*TLOW(I,J)+BREMS*WPO(I,J,K)
          CONTRIJ(I,J)=CONTRIJ(I,J)+(OVE-WPO(I,J,K))**2
        ENDDO
      ENDDO
    ENDDO
    
!$OMP SINGLE
!#ifdef bounds_exch_save
    CALL bounds_exch('p',WPO,'occlit 10')
!#endif
    CALL gather_arr(contrij,contrij_g,p_io)
!$OMP END SINGLE  

  if (icontro.ne.0) then 
  IF(p_pe==p_io) THEN
    CONTRA=0.
    DO  J=2,JE_G-1
      DO  I=2,IE_G-1
        CONTRA=CONTRA+CONTRIJ_G(I,J)
      ENDDO
    ENDDO
    WRITE(0,*)'ITER CONTRA: ',iter,CONTRA,maxloc(CONTRIJ_G)
  ENDIF
  endif

  END DO iteration_loop
  
!$OMP DO
  DO K=1,KE
    DO J=1,JE
      DO I=1,IE
        uko(i,j,k)=uko(i,j,k)*amsuo(i,j,k)/(almzer+dduo(i,j,k))
        vke(i,j,k)=vke(i,j,k)*amsue(i,j,k)/(almzer+ddue(i,j,k))
        uk1o(i,j,k)=uk1o(i,j,k)*amsuo(i,j,k)/(almzer+dduo(i,j,k))
        vk1e(i,j,k)=vk1e(i,j,k)*amsue(i,j,k)/(almzer+ddue(i,j,k))
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL
  !
  CALL gather_arr(contrij,contrij_g,p_io)
  
  IF(p_pe==p_io) THEN
    CONTRA=0.
    DO  J=2,JE_G-1
      DO  I=2,IE_G-1
        CONTRA=CONTRA+CONTRIJ_G(I,J)
      ENDDO
    ENDDO
    !
    WRITE(IO_STDOUT,*)'CONTRA: ',CONTRA
    !
    IF(CONTRA.GT.1.E-3) THEN
      WRITE(IO_STDOUT,*) 'CONTRI : ',SUM(CONTRIJ_G(2:IE_G-1,2:JE_G-1),1)
      WRITE(IO_STDOUT,*) ' '
      WRITE(IO_STDOUT,*) 'CONTRJ : ',SUM(CONTRIJ_G(2:IE_G-1,2:JE_G-1),2)
      CALL ABSTURZ
    ENDIF
  ENDIF

END SUBROUTINE OCCLIT



