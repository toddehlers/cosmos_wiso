SUBROUTINE OCMODMOM

  USE MO_PARAM1
  USE MO_PARALLEL
  USE MO_COMMO1
  USE MO_UNITS

  IMPLICIT NONE
  !
  !-=====================================================================
  !
  !    DECOMPOSITION INTO BAROTROPIC AND BAROCLINIC FIELD
  !
  REAL UCOR(IE,JE,KE),VCOR(IE,JE,KE)
  INTEGER i,j,k,iter
  !
  !---------------------------------------------------------------------
  !
!$OMP PARALLEL PRIVATE(i,j,k,iter)
  !  
!$OMP DO
  DO J=1,JE
    DO I=1,IE
      U1O(I,J)=ZERO
      U1E(I,J)=ZERO
      V1O(I,J)=ZERO
      V1E(I,J)=ZERO
      USO(I,J)=ZERO
      VSE(I,J)=ZERO
#ifdef bounds_exch_tp
       CURVAV(I,J)=0.
#endif
    ENDDO
  ENDDO
  !


!$OMP SINGLE

!#ifdef bounds_exch_save
    CALL bounds_exch('u',UKO,'ocmodmom 1')
    CALL bounds_exch('v',VKE,'ocmodmom 2')
    CALL bounds_exch('u',UOO,'ocmodmom 3')
    CALL bounds_exch('v',VOE,'ocmodmom 4')
    CALL bounds_exch('u+',dduo,'ocmodmom 7')
    CALL bounds_exch('v+',ddue,'ocmodmom 8')
!#endif

    CALL bounds_exch('p',po,'ocmodmom 5')
    CALL bounds_exch('p',zo,'ocmodmom 6')

!$OMP END SINGLE
#ifdef bounds_exch_tp
    call vcheck(1,voe)
#endif
!$OMP DO
  DO K=1,KE
    DO J=1,JE
      DO I=1,IE
        UKO(I,J,K)=UKO(I,J,K)*DDUO(I,J,K)
        VKE(I,J,K)=VKE(I,J,K)*DDUE(I,J,K)
        UK1O(I,J,K)=UKO(I,J,K)
        VK1E(I,J,K)=VKE(I,J,K)
        STABIO(I,J,K)=0.001*DZ(K)*STABIO(I,J,K)*WETO(I,J,K)
        PO(I,J,K)=PO(I,J,K)+G*ZO(I,J)
      ENDDO
    ENDDO
  ENDDO
  !
  iteration_loop: DO ITER=1,12
    !
!$OMP DO
    DO J=1,JE
      DO K=1,KE
        DO I=1,IE
          UCOR(I,J,K)=STABN*UK1O(I,J,K)+STABO*UKO(I,J,K)
          VCOR(I,J,K)=STABN*VK1E(I,J,K)+STABO*VKE(I,J,K)
        ENDDO
      ENDDO
    ENDDO
    !



!$OMP DO
    DO J=2,JE1
      DO K=1,KE
        DO I=2,IE1


          !CUWENEU   INCLUDE EARTH CURVATURE
          !C
          UK1O(I,J,K)=UKO(I,J,K)+DTDXUO(I,J)*(PO(I,J,K)-PO(I+1,J,K))   &
               *DDUO(I,J,K)                                            &
               +0.25*DT*(                                              & 
               FTWOU(I,J)                                              &
               +0.25*(CURVAV(I,J)+CURVAV(I+1,J)+CURVAV(I,J-1)          &
               +CURVAV(I+1,J))*UCOR(I,J,K)                             &
               )*(VCOR(I,J,K)+VCOR(I+1,J,K)                            &
               +VCOR(I,J-1,K)+VCOR(I+1,J-1,K))*AMSUO(I,J,K)

          VK1E(I,J,K)=VKE(I,J,K)+DPYE(I,J)*(PO(I,J+1,K)-PO(I,J,K))     &
               *DDUE(I,J,K)                                            &
               -0.25*DT*(FTWOV(I,J)                                    &
               +CURVAV(I,J)*0.25                                       &
               *(UCOR(I,J,K)+UCOR(I-1,J,K)+UCOR(I,J+1,K)               &
               +UCOR(I-1,J+1,K))                                       &
               )*(UCOR(I,J,K)+UCOR(I-1,J,K)                            &
               +UCOR(I,J+1,K)+UCOR(I-1,J+1,K))*AMSUE(I,J,K)
        ENDDO
      ENDDO
    ENDDO
    !
!#ifdef bounds_exch_save
!$OMP SINGLE
    CALL bounds_exch('u',UK1O,'ocmodmom 9')
    CALL bounds_exch('v',VK1E,'ocmodmom 10')
!$OMP END SINGLE
!#endif
    !
  ENDDO iteration_loop
  !
  !     CALCULATION OF BAROTROPIC VELOCITIES ON FIELDS U1 AND V1
  !
!$OMP DO
  DO J=2,JE1
    DO K=1,KE
      DO I=1,IE
        po(i,j,k) = po(i,j,k)-g*zo(i,j)
        v1e(i,j)  = v1e(i,j)+ddue(i,j,k)*voe(i,j,k)
        u1o(i,j)  = u1o(i,j)+dduo(i,j,k)*uoo(i,j,k)
        vse(i,j)  = vse(i,j)+vk1e(i,j,k)*amsue(i,j,k)
        uso(i,j)  = uso(i,j)+uk1o(i,j,k)*amsuo(i,j,k)
      ENDDO
    ENDDO
  ENDDO
  !
!$OMP DO
  DO K=1,KE
    DO J=2,JE1
      DO I=1,IE
        uko(i,j,k)=uoo(i,j,k)-amsuo(i,j,k)*deutio(i,j)*u1o(i,j)
        vke(i,j,k)=voe(i,j,k)-amsue(i,j,k)*deutie(i,j)*v1e(i,j)
        uk1o(i,j,k)=uk1o(i,j,k)*(amsuo(i,j,k)/(almzer+dduo(i,j,k))) &
              -amsuo(i,j,k)*deutio(i,j)*uso(i,j)
        vk1e(i,j,k)=vk1e(i,j,k)*(amsue(i,j,k)/(almzer+ddue(i,j,k))) &
              -amsue(i,j,k)*deutie(i,j)*vse(i,j)
      ENDDO
    ENDDO
  ENDDO
  !


#ifdef bounds_exch_tp
  call vcheck(3,vk1e)
#endif

!$OMP SINGLE
    CALL bounds_exch('u',UK1O,'ocmodmom 9')
    CALL bounds_exch('v',VK1E,'ocmodmom 11')
    CALL bounds_exch('u',UKO,'ocmodmom 9')
    CALL bounds_exch('v',VKE,'ocmodmom 12')
    CALL bounds_exch('u',U1O,'ocmodmom 9')
    CALL bounds_exch('v',V1E,'ocmodmom 10')
    CALL bounds_exch('u',UsO,'ocmodmom 9')
    CALL bounds_exch('v',VsE,'ocmodmom 10')
    CALL bounds_exch('p',po,'ocmodmom 10')
!$OMP END SINGLE



!$OMP END PARALLEL



  !
END SUBROUTINE OCMODMOM
