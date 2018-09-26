MODULE mo_tro

! module file containing the sbrs for iterative sea level solver (SOR)

IMPLICIT NONE

  REAL ,ALLOCATABLE :: dilcor(:,:),susalo(:,:),suwath(:,:)

CONTAINS

  SUBROUTINE itprep

    ! rj: rewritten for mpi-parallelization nov 2003
    ! calculates some arrays needed for the iterative solver

    USE mo_param1, ONLY: ie,je
    USE mo_para2, ONLY: xx, uf, vf, ff
    USE mo_commo1, ONLY: dlxp, dlyp,dlxu,dlxv,dlyu,dlyv,deute,deuto &
                        ,g,conn,stabn,dt
    USE mo_parallel, ONLY: bounds_exch

    IMPLICIT NONE

    INTEGER i,j

    DO i=2,ie-1
       DO j=2,je-1
          xx(i,j) = dlxp(i,j)*dlyp(i,j)
          uf(i,j) = deuto(i,j)*(dlyu(i,j)/dlxu(i,j))*g*conn*stabn*dt**2
          vf(i,j) = deute(i,j)*(dlxv(i,j)/dlyv(i,j))*g*conn*stabn*dt**2
       ENDDO
    ENDDO

    CALL bounds_exch('p',xx,'mo_tro 1')
    CALL bounds_exch('u+',uf,'mo_tro 2')
    CALL bounds_exch('v+',vf,'mo_tro 3')

    ff(:,:) = 0
    DO j=2,je-1
       DO i=2,ie-1
          ff(i,j) = 1./(xx(i,j)+uf(i,j)+uf(i-1,j)+vf(i,j)+vf(i,j-1))
       ENDDO
    ENDDO

    CALL bounds_exch('p',ff,'mo_tro 4')

  END SUBROUTINE itprep

  SUBROUTINE trotest_lk

    USE mo_param1,   ONLY: ie_g, je_g
    USE mo_para2,    ONLY: xx, uf, vf, ff, sorpar
    USE mo_commo1,   ONLY: weto_g,z1o,b1o
    USE mo_units,    ONLY: io_stdout
    USE mo_parallel, ONLY: p_pe, p_io, nprocxy, p_recv, p_send, p_bcast, &
         gather_arr,scatter_arr
    !
    !=======================================================================
    !     purpose :
    !
    !     
    !uwe  preparation for iterative solution of zeta field
    !     optimises sorpar
    !
    ! rj: rewritten for mpi-parallelization nov 2003
    !

    IMPLICIT NONE

    INTEGER, PARAMETER :: isormax = 250

    REAL :: sor_omega(isormax), convergence(isormax)

    REAL :: b1o_g(ie_g,je_g), z1o_g(ie_g,je_g)
    REAL :: xx_g(ie_g,je_g) 
    REAL :: uf_g(ie_g,je_g), vf_g(ie_g,je_g), ff_g(ie_g,je_g)

    INTEGER :: i, j, is, n, iss, ies, isor, kiter, isor_start, isor_end
    INTEGER :: i0(1),jb
    REAL :: chabas, sorpai, sumcha, zalt
    REAL :: sqrtrnd, rnd

    ! build global fields on io pe

    CALL gather_arr(xx,xx_g,p_io)
    CALL gather_arr(uf,uf_g,p_io)
    CALL gather_arr(vf,vf_g,p_io)
    CALL gather_arr(ff,ff_g,p_io)

!    ! distribute global fields to all pe 
!
!    CALL p_bcast(xx_g,p_io)
!    CALL p_bcast(uf_g,p_io)
!    CALL p_bcast(vf_g,p_io)
!    CALL p_bcast(ff_g,p_io)

    !uwe  initialize with random numbers
    !hh random_number doesn't work with fujitsu f90!
    !hh    call random_number(xhh)
    !hh    b1o(i,j)=(xhh-0.5)*weto(i,j,1)
    !otb  using a simple one

    ! rj: we do that in the following way so that the result
    !     does not depend on the decomposition
    !     why is that done for every step in the original version ????
    !     we do it only once here

    IF (p_pe==p_io) THEN 
    b1o_g(:,:) = 0.0
    DO j = 1, je_g
       DO i = 2, ie_g-1
          rnd = sqrtrnd()
          IF (i >=1 .AND. i <= ie_g .AND. j >= 1 .AND. j <= je_g) THEN
             b1o_g(i,j) = (rnd-0.5)*weto_g(i,j,1)*xx_g(i,j)
          ENDIF
       ENDDO
    ENDDO
    ENDIF

    CALL scatter_arr(b1o_g,b1o,p_io)

!    b1o_g(ie_g,:) = b1o_g(2,:)
!    b1o_g(1,:)    = b1o_g(ie_g-1,:)

    jb=2
    DO isor = 1, isormax
       sor_omega(isor) = 1.70+isor*0.001
    ENDDO

    convergence(:) = 9.e99

    CALL decompose (isormax, nprocxy, p_pe, isor_start, isor_end)

    !$omp parallel private (isor, sorpar, sorpai, z1o_g, kiter, sumcha, i, j, is, &
    !$omp                   zalt, chabas)

    !$omp do
    DO isor = isor_start, isor_end
       !
       sorpar = sor_omega(isor)
       sorpai = 1.0-sorpar
       !
       z1o_g(:,:) = 0.
       !
       DO kiter = 1, 250

         IF (p_pe==p_io) THEN

          sumcha = 0.0

          ! odd/odd and even/even z1o elements

          DO j=jb,je_g-1
             is = MOD(j,2)+2
             DO i=is,ie_g-1,2
                zalt = z1o_g(i,j)
                z1o_g(i,j) = sorpar*ff_g(i,j)*(b1o_g(i,j)                      &
                     +  uf_g(i,j)*z1o_g(i+1,j) + uf_g(i-1,j)*z1o_g(i-1,j)      &
                     +  vf_g(i,j)*z1o_g(i,j+1) + vf_g(i,j-1)*z1o_g(i,j-1))     &
                     + sorpai*z1o_g(i,j)
                sumcha=sumcha+(zalt-z1o_g(i,j))**2
             ENDDO
          ENDDO

          ENDIF

          CALL scatter_arr(z1o_g,z1o,p_io)
          CALL gather_arr(z1o,z1o_g,p_io)


!          z1o_g(ie_g,:) = z1o_g(2,:)
!          z1o_g(1,:)    = z1o_g(ie_g-1,:)

          ! odd/even and even/odd z1o elements

	  IF (p_pe==p_io) THEN

          DO j=jb,je_g-1
             is = MOD(j+1,2)+2
             DO i=is,ie_g-1,2
                zalt = z1o_g(i,j)
                z1o_g(i,j) = sorpar*ff_g(i,j)*(b1o_g(i,j)                      &
                     +  uf_g(i,j)*z1o_g(i+1,j) + uf_g(i-1,j)*z1o_g(i-1,j)      &
                     +  vf_g(i,j)*z1o_g(i,j+1) + vf_g(i,j-1)*z1o_g(i,j-1))     &
                     + sorpai*z1o_g(i,j)
                sumcha=sumcha+(zalt-z1o_g(i,j))**2
             ENDDO
          ENDDO

	  ENDIF

          CALL scatter_arr(z1o_g,z1o,p_io)
          CALL gather_arr(z1o,z1o_g,p_io)

!          z1o_g(ie_g,:) = z1o_g(2,:)
!          z1o_g(1,:)    = z1o_g(ie_g-1,:)

          IF (kiter == 1) chabas = sumcha

       ENDDO ! kiter loop

       convergence(isor) = sumcha/chabas

    ENDDO ! isor loop

    !$omp end parallel

    IF (p_pe == p_io) THEN
       DO n = 0, nprocxy
          IF  (p_pe == p_io) CYCLE
          CALL decompose (isormax, nprocxy, n, iss, ies)
          CALL p_recv(convergence(iss), n, 200, &
               p_count=ies-iss+1) 
       ENDDO
    ELSE
       CALL p_send(convergence(isor_start),p_io, 200, &
            p_count=isor_end-isor_start+1)
    ENDIF

    IF (p_pe == p_io) THEN
       i0 = MINLOC(convergence)
       sorpar = sor_omega(i0(1))
    ENDIF

    CALL p_bcast(sorpar,p_io)

    WRITE(io_stdout,*) 'trotest: sor parameter = ',sorpar

  CONTAINS 

    SUBROUTINE decompose (ns, npes, pe, istart, iend)

      INTEGER, INTENT(in) :: ns       ! domain size
      INTEGER, INTENT(in) :: npes     ! pes availbale for that decomposition
      INTEGER, INTENT(in) :: pe       ! this pe

      INTEGER, INTENT(out) :: istart  ! start index of sub domain
      INTEGER, INTENT(out) :: iend    ! end index of sub domain

      INTEGER :: nlocal, nremain

      nlocal = ns/npes
      istart = pe*nlocal+1
      nremain = MOD(ns,npes)
      istart = istart+MIN(pe,nremain)
      IF (pe < nremain) THEN
         nlocal = nlocal+1
      ENDIF
      iend = istart+nlocal-1
      IF (iend > ns .OR. pe == npes-1) iend = ns

    END SUBROUTINE decompose

  END SUBROUTINE trotest_lk

  SUBROUTINE TROTEST
      USE MO_PARAM1
      USE MO_PARA2
      USE MO_COMMO1
      USE MO_UNITS

      USE MO_PARALLEL
!
!=======================================================================
!     PURPOSE :
!
!
!UWE  PREPARATION FOR ITERATIVE SOLUTION OF ZETA FIELD
!     OPTIMISES SORPAR
!
! RJ: Rewritten for MPI-Parallelization Nov 2003
!

      IMPLICIT NONE
      INTEGER I, J, II, JJ, IS, ISOR, KITER
      REAL CHABAS, CONTRIC, CONVMIN, SORPAI, SORPMIN, SUMCHA, ZALT
      REAL SQRTRND, RND

      SORPAR=1.75

      CONVMIN=9.E99
      SORPMIN=SORPAR
!
      DO ISOR=1,200
      SORPAR=SORPAR+0.001
!
!UWE  INITIALIZE WITH RANDOM NUMBERS
!hh random_number doesn't work with fujitsu f90!
!hh    CALL RANDOM_NUMBER(XHH)
!hh    B1O(I,J)=(XHH-0.5)*WETO(I,J,1)
!OtB  Using a simple one

! RJ: We do that in the following way so that the result
!     does not depend on the decomposition
!     Why is that done for every step in the original version ????
!     We do it only once here

      IF (ISOR==1) THEN
        DO JJ=1,JE_G
        DO II=2,IE_G-1
          I = II - p_ioff
          J = JJ - p_joff
          RND = SQRTRND()
          IF(I>=1 .AND. I<=IE .AND. J>=1 .AND. J<=JE) THEN
            B1O(I,J)=(RND-0.5)*WETO(I,J,1)*XX(I,J)
          ENDIF
        ENDDO
        ENDDO

        CALL bounds_exch('p',B1O,'mo_tro 5')
      ENDIF

      Z1O(:,:) = 0.

      SORPAI=1.-SORPAR

      DO KITER=1,200

        SUMCHA=0.

        ! odd/odd and even/even Z1O elements

        DO J=2,JE-1
          IS = MOD(p_ioff+p_joff+j,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch('p',Z1O,'mo_tro 6')

        ! odd/even and even/odd Z1O elements

        DO J=2,JE-1
          IS = MOD(p_ioff+p_joff+j+1,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch('p',Z1O,'mo_tro 7')

        CALL global_sum(SUMCHA)

        IF(KITER.EQ.1) CHABAS=SUMCHA

      ENDDO ! KITER Loop

      CONTRIC=SUMCHA/CHABAS

      IF (icontro.NE.0)THEN
         WRITE(6,*)'sorpa convergence ',sorpar,contric
      ENDIF

      IF(CONTRIC.LT.CONVMIN)THEN
        SORPMIN=SORPAR
        CONVMIN=CONTRIC
      ENDIF

      ENDDO ! ISOR Loop

      SORPAR=SORPMIN
      WRITE(IO_STDOUT,*) 'TROTEST SORPAR=',SORPAR
    
    END SUBROUTINE TROTEST

    SUBROUTINE TROTEST2
      USE MO_PARAM1
      USE MO_PARA2
      USE MO_COMMO1
      USE MO_UNITS

      USE MO_PARALLEL
!
!=======================================================================
!     PURPOSE :
!
!
!UWE  PREPARATION FOR ITERATIVE SOLUTION OF ZETA FIELD
!     OPTIMISES SORPAR
!
! RJ: Rewritten for MPI-Parallelization Nov 2003
!

      IMPLICIT NONE
      INTEGER I, J, II, JJ, IS, ISOR, KITER,jb
      REAL CHABAS, CONTRIC, CONVMIN, SORPAI, SORPMIN, SUMCHA, ZALT
      REAL SQRTRND, RND

      SORPAR=0.9
!
      CONVMIN=9.E99
      SORPMIN=SORPAR

      jb=2
#ifdef bounds_exch_tp
     if(p_joff.eq.0)jb=3
#endif


!     1) find SORPMIN between 1. and 2.0
      DO ISOR=1,10
      SORPAR=SORPAR+0.1

      IF (ISOR==1) THEN
        DO JJ=1,JE_G
        DO II=2,IE_G-1
          I = II - p_ioff
          J = JJ - p_joff
          RND = SQRTRND()
          IF(I>=1 .AND. I<=IE .AND. J>=1 .AND. J<=JE) THEN
            B1O(I,J)=(RND-0.5)*WETO(I,J,1)*XX(I,J)
          ENDIF
        ENDDO
        ENDDO

        CALL bounds_exch('p',B1O,'mo_tro 5')
      ENDIF

      Z1O(:,:) = 0.

      SORPAI=1.-SORPAR

      DO KITER=1,200

        SUMCHA=0.

        ! odd/odd and even/even Z1O elements

        DO J=jb,JE-1
          IS = MOD(p_ioff+p_joff+j,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch('p',Z1O,'mo_tro 6')

        ! odd/even and even/odd Z1O elements

        DO J=jb,JE-1
          IS = MOD(p_ioff+p_joff+j+1,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch('p',Z1O,'mo_tro 7')

        CALL global_sum(SUMCHA)

        IF(KITER.EQ.1) CHABAS=SUMCHA

      ENDDO ! KITER Loop

      CONTRIC=SUMCHA/CHABAS

      IF (icontro.NE.0)THEN
         WRITE(6,*)'sorpa convergence 1: ',sorpar,contric
      ENDIF


      IF(CONTRIC.LT.CONVMIN)THEN
        SORPMIN=SORPAR
        CONVMIN=CONTRIC
      ENDIF

      ENDDO ! ISOR Loop


!     2) find SORPMIN between -0.1 and +0.1

      convmin=9.e99
      SORPAR=SORPMIN-0.1
      SORPMIN=SORPAR

      DO ISOR=1,20
      SORPAR=SORPAR+0.01

!!$      IF (ISOR==1) THEN
!!$        DO JJ=1,JE_G
!!$        DO II=2,IE_G-1
!!$          I = II - p_ioff
!!$          J = JJ - p_joff
!!$          RND = SQRTRND()
!!$          IF(I>=1 .AND. I<=IE .AND. J>=1 .AND. J<=JE) THEN
!!$            B1O(I,J)=(RND-0.5)*WETO(I,J,1)*XX(I,J)
!!$          ENDIF
!!$        ENDDO
!!$        ENDDO
!!$
!!$        CALL bounds_exch('p',B1O,'mo_tro 5')
!!$      ENDIF

      Z1O(:,:) = 0.

      SORPAI=1.-SORPAR

      DO KITER=1,200

        SUMCHA=0.

        ! odd/odd and even/even Z1O elements

        DO J=jb,JE-1
          IS = MOD(p_ioff+p_joff+j,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch('p',Z1O,'mo_tro 6')

        ! odd/even and even/odd Z1O elements

        DO J=jb,JE-1
          IS = MOD(p_ioff+p_joff+j+1,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch('p',Z1O,'mo_tro 7')

        CALL global_sum(SUMCHA)

        IF(KITER.EQ.1) CHABAS=SUMCHA

      ENDDO ! KITER Loop

      CONTRIC=SUMCHA/CHABAS

      IF (icontro.NE.0)THEN
         WRITE(6,*)'sorpa convergence 2: ',sorpar,contric
      ENDIF


      IF(CONTRIC.LT.CONVMIN)THEN
        SORPMIN=SORPAR
        CONVMIN=CONTRIC
      ENDIF

      ENDDO ! ISOR Loop

!     3) find SORPMIN between -0.01 and +0.01
      
      convmin=9.e99
      SORPAR=SORPMIN-0.01
      SORPMIN=SORPAR

      DO ISOR=1,20
      SORPAR=SORPAR+0.001

!!$      IF (ISOR==1) THEN
!!$        DO JJ=1,JE_G
!!$        DO II=2,IE_G-1
!!$          I = II - p_ioff
!!$          J = JJ - p_joff
!!$          RND = SQRTRND()
!!$          IF(I>=1 .AND. I<=IE .AND. J>=1 .AND. J<=JE) THEN
!!$            B1O(I,J)=(RND-0.5)*WETO(I,J,1)*XX(I,J)
!!$          ENDIF
!!$        ENDDO
!!$        ENDDO
!!$
!!$        CALL bounds_exch('p',B1O,'mo_tro 5')
!!$      ENDIF

      Z1O(:,:) = 0.

      SORPAI=1.-SORPAR

      DO KITER=1,200

        SUMCHA=0.

        ! odd/odd and even/even Z1O elements

        DO J=jb,JE-1
          IS = MOD(p_ioff+p_joff+j,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch('p',Z1O,'mo_tro 6')

        ! odd/even and even/odd Z1O elements

        DO J=jb,JE-1
          IS = MOD(p_ioff+p_joff+j+1,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch('p',Z1O,'mo_tro 7')

        CALL global_sum(SUMCHA)

        IF(KITER.EQ.1) CHABAS=SUMCHA

      ENDDO ! KITER Loop

      CONTRIC=SUMCHA/CHABAS

      IF (icontro.GT.0)THEN
         WRITE(6,*)'sorpa convergence 3: ',sorpar,contric
      ENDIF


      IF(CONTRIC.LT.CONVMIN)THEN
        SORPMIN=SORPAR
        CONVMIN=CONTRIC
      ENDIF

      ENDDO ! ISOR Loop




      SORPAR=SORPMIN
      WRITE(IO_STDOUT,*) 'TROTEST SORPAR2=',SORPAR
    
    END SUBROUTINE TROTEST2


  SUBROUTINE troneu
    !**********************************************************************
    !
    !
    !     tttttt  rrrrr    ooo    nn   n  eeeeee  u    u
    !       tt    r    r  o   o   n n  n  e       u    u
    !       tt    rrrrr   o   o   n  n n  eeeee   u    u
    !       tt    r rr    o   o   n   nn  e       u    u
    !       tt    r   rr   ooo    n    n  eeeeee  uuuuuu
    !
    !
    !**********************************************************************
    !
    USE mo_param1, ONLY: ie,je
    USE mo_para2, ONLY: ff,uf,vf,sorpar
    USE mo_commo1, ONLY: weto,dt,cono,u1o,uzo,conn,dlyu,v1e,vze   &
                        ,dlxv,b1o,z1o,icontro
    USE mo_units
    USE mo_parallel, ONLY: p_ioff,p_joff,bounds_exch

    IMPLICIT NONE
    INTEGER i, j, is, kiter,jb
    REAL sorpai,zurr,zalt
    !
    !=======================================================================
    !     sbr tropit
    !
    !     purpose :
    !
    !     a) iterative solution of zeta field
    !
    !
    !   iterative solution of the system
    !   z*(dx*dy+g*dt**2(hw*dyw/dxw+ho*dyo/dxo+hs*dxs/dys+hn*dxn/dyn) )
    !
    ! =            g*dt**2*(hw*zww*dyw/dxw+ho*zoo*dyo/dxo
    !                      +hs*zss*dxs/dys+hn*znn*dxn/dyn  )
    !+g*dt**3*(f*hw*(zsw-znw)+f*ho*(zno-zso)+fs*hs*(zso-zsw)+fn*hn*(znw-zno)
    !                 +b
    !
    !   where b contains the winstress, divergence of old flow and old z
    !
    !
    ! rj: rewritten for mpi-parallelization nov 2003


    DO j=2,je-1
       DO i=2,ie-1

          b1o(i,j) = weto(i,j,1) *dt* (                                       &
                     (cono * u1o(i-1,j) +  conn * uzo(i-1,j) ) * dlyu(i-1,j)  &
                   - (cono * u1o(i,j)   +  conn * uzo(i,j)   ) * dlyu(i,j)    &
                   + (cono * v1e(i,j)   +  conn * vze(i,j)   ) * dlxv(i,j)    &
                   - (cono * v1e(i,j-1) +  conn * vze(i,j-1) ) * dlxv(i,j-1))

       ENDDO
    ENDDO

    CALL bounds_exch('p',b1o,'mo_tro 8')

    DO j=1,je
       DO i=1,ie
          z1o(i,j)=0.
       ENDDO
    ENDDO


     jb=2
#ifdef bounds_exch_tp
      if(p_joff.eq.0)jb=3
#endif



    sorpai=1.-sorpar
    
    DO kiter=1,300
       zurr=0.

       ! SOR with red/black splitting
       ! 1step: update of the "red" points 

             IF (icontro.NE.0) THEN

       DO j=jb,je-1
          is = MOD(p_ioff+p_joff+j,2)+2
          DO i=is,ie-1,2
             zalt=z1o(i,j)
             z1o(i,j) = sorpar*ff(i,j)*(b1o(i,j)                         &
                      +  uf(i,j)*z1o(i+1,j) + uf(i-1,j)*z1o(i-1,j)       &
                      +  vf(i,j)*z1o(i,j+1) + vf(i,j-1)*z1o(i,j-1))      &
                      + sorpai*z1o(i,j)
                zurr=zurr+(zalt-z1o(i,j))**2
          ENDDO
       ENDDO

             ELSE

       DO j=jb,je-1
          is = MOD(p_ioff+p_joff+j,2)+2
          DO i=is,ie-1,2
             zalt=z1o(i,j)
             z1o(i,j) = sorpar*ff(i,j)*(b1o(i,j)                         &
                      +  uf(i,j)*z1o(i+1,j) + uf(i-1,j)*z1o(i-1,j)       &
                      +  vf(i,j)*z1o(i,j+1) + vf(i,j-1)*z1o(i,j-1))      &
                      + sorpai*z1o(i,j)
          ENDDO
       ENDDO

             ENDIF
    
       CALL bounds_exch('p',z1o,'mo_tro 9')

       ! 2step: update of the "black" points 

             IF (icontro.NE.0) THEN

       DO j=jb,je-1
          is = MOD(p_ioff+p_joff+j+1,2)+2
          DO i=is,ie-1,2
             zalt=z1o(i,j)
             z1o(i,j) = sorpar*ff(i,j)*(b1o(i,j)                        &
                      +  uf(i,j)*z1o(i+1,j) + uf(i-1,j)*z1o(i-1,j)      &
                      +  vf(i,j)*z1o(i,j+1) + vf(i,j-1)*z1o(i,j-1))     &
                      + sorpai*z1o(i,j)
                zurr=zurr+(zalt-z1o(i,j))**2
          ENDDO
       ENDDO

             ELSE

       DO j=jb,je-1
          is = MOD(p_ioff+p_joff+j+1,2)+2
          DO i=is,ie-1,2
             zalt=z1o(i,j)
             z1o(i,j) = sorpar*ff(i,j)*(b1o(i,j)                        &
                      +  uf(i,j)*z1o(i+1,j) + uf(i-1,j)*z1o(i-1,j)      &
                      +  vf(i,j)*z1o(i,j+1) + vf(i,j-1)*z1o(i,j-1))     &
                      + sorpai*z1o(i,j)
          ENDDO
       ENDDO

             ENDIF

       CALL bounds_exch('p',z1o,'mo_tro 10')

    ENDDO
    IF (icontro.NE.0) THEN
    WRITE(0,*)'zurr',zurr,'sorpar',sorpar
    ENDIF
  END SUBROUTINE troneu

  SUBROUTINE troneu2
    !**********************************************************************
    !
    !
    !     tttttt  rrrrr    ooo    nn   n  eeeeee  u    u
    !       tt    r    r  o   o   n n  n  e       u    u
    !       tt    rrrrr   o   o   n  n n  eeeee   u    u
    !       tt    r rr    o   o   n   nn  e       u    u
    !       tt    r   rr   ooo    n    n  eeeeee  uuuuuu
    !
    !
    !**********************************************************************
    !
    USE mo_param1, ONLY: ie,je
    USE mo_para2, ONLY: ff,uf,vf,sorpar
    USE mo_commo1, ONLY: weto,dt,cono,u1o,uzo,conn,dlyu,v1e,vze   &
                        ,dlxv,b1o,z1o,icontro
    USE mo_units
    USE mo_parallel, ONLY: p_ioff,p_joff,bounds_exch,sethalo2,sethaloN,p_pe

    IMPLICIT NONE
    INTEGER i, j, is, kiter,jb,ii,jj
    Integer, parameter :: ihalo=2,ihm=2*(ihalo-1) 
    REAL sorpai,zurr,zalt(ie,je),zz1o(ie+ihm,je+ihm)
    REAL zff(ie+ihm,je+ihm),zb1o(ie+ihm,je+ihm),zuf(ie+ihm,je+ihm),zvf(ie+ihm,je+ihm)
    !
    !=======================================================================
    !     sbr tropit
    !
    !     purpose :
    !
    !     a) iterative solution of zeta field
    !
    !
    !   iterative solution of the system
    !   z*(dx*dy+g*dt**2(hw*dyw/dxw+ho*dyo/dxo+hs*dxs/dys+hn*dxn/dyn) )
    !
    ! =            g*dt**2*(hw*zww*dyw/dxw+ho*zoo*dyo/dxo
    !                      +hs*zss*dxs/dys+hn*znn*dxn/dyn  )
    !+g*dt**3*(f*hw*(zsw-znw)+f*ho*(zno-zso)+fs*hs*(zso-zsw)+fn*hn*(znw-zno)
    !                 +b
    !
    !   where b contains the winstress, divergence of old flow and old z
    !
    !
    ! rj: rewritten for mpi-parallelization nov 2003


    DO j=2,je-1
       DO i=2,ie-1

          b1o(i,j) = weto(i,j,1) *dt* (                                       &
                     (cono * u1o(i-1,j) +  conn * uzo(i-1,j) ) * dlyu(i-1,j)  &
                   - (cono * u1o(i,j)   +  conn * uzo(i,j)   ) * dlyu(i,j)    &
                   + (cono * v1e(i,j)   +  conn * vze(i,j)   ) * dlxv(i,j)    &
                   - (cono * v1e(i,j-1) +  conn * vze(i,j-1) ) * dlxv(i,j-1))

       ENDDO
    ENDDO

    CALL bounds_exch('p',b1o,'mo_tro 8')

    DO j=1,je
       DO i=1,ie
          z1o(i,j)=0.
       ENDDO
    ENDDO

    DO j=1,je+ihm
       DO i=1,ie+ihm
          zz1o(i,j)=0.
       ENDDO
    ENDDO


     jb=2
#ifdef bounds_exch_tp
      if(p_joff.eq.0)jb=3
#endif


    sorpai=1.-sorpar
 

!    do i=1,ie
!    do j=1,je
!       z1o(i,j)=float(i)
!    enddo
!    enddo



!    call sethalo2(z1o,zz1o,2)

!    do j=1,je+2
!    write(0,*) '2',j,zz1o(1,j),zz1o(2,j),zz1o(3,j),zz1o(4,j)
!    enddo

   call sethaloN(z1o,zz1o,ihalo)

!    do j=1,je+2
!    write(0,*) 'links ',p_pe,j,zz1o(1,j),zz1o(2,j),zz1o(3,j),zz1o(4,j)
!    enddo

!    do j=1,je+2
!    write(0,*) 'rechts',p_pe,j,zz1o(ie-1,j),zz1o(ie,j),zz1o(ie+1,j),zz1o(ie+2,j)
!    enddo




!   stop


    call sethaloN(b1o,zb1o,ihalo)
    call sethaloN(ff,zff,ihalo)
    call sethaloN(uf,zuf,ihalo)
    call sethaloN(vf,zvf,ihalo)


   
    DO kiter=1,300
       zurr=0.

       ! SOR with red/black splitting
       ! 1step: update of the "red" points 

       if ( icontro.NE.0 ) THEN      
          DO j=jb,je-1
             is = MOD(p_ioff+p_joff+j,2)+2
             DO i=is,ie-1,2
                zalt(i,j)=z1o(i,j) 
             ENDDO
          ENDDO
       endif
          
       DO j=jb,je+ihm-1
          is = MOD(p_ioff+p_joff+j,2)+2
          DO i=is,ie+ihm-1,2
             zz1o(i,j) = sorpar*zff(i,j)*(zb1o(i,j)                      &
                  +  zuf(i,j)*zz1o(i+1,j) + zuf(i-1,j)*zz1o(i-1,j)       &
                  +  zvf(i,j)*zz1o(i,j+1) + zvf(i,j-1)*zz1o(i,j-1))      &
                  + sorpai*zz1o(i,j)
          ENDDO
       ENDDO

       if ( icontro.NE.0 ) THEN      
          do i=1,ie
             do j=1,je
                z1o(i,j)=zz1o(i+ihalo-1,j+ihalo-1)
             enddo
          enddo
          
          DO j=jb,je-1
             is = MOD(p_ioff+p_joff+j,2)+2
             DO i=is,ie-1,2
                zurr=zurr+(zalt(i,j)-z1o(i,j))**2
             enddo
          enddo
        endif
          
       ! 2step: update of the "black" points 
       
       if ( icontro.NE.0 ) THEN      
          DO j=jb,je-1
             is = MOD(p_ioff+p_joff+j+1,2)+2
             DO i=is,ie-1,2
                zalt(i,j)=z1o(i,j) 
             ENDDO
          ENDDO
       endif
       
       
       DO j=jb,je+ihm-1
          is = MOD(p_ioff+p_joff+j+1,2)+2
          DO i=is,ie+ihm-1,2
             zz1o(i,j) = sorpar*zff(i,j)*(zb1o(i,j)                      &
                  +  zuf(i,j)*zz1o(i+1,j) + zuf(i-1,j)*zz1o(i-1,j)       &
                  +  zvf(i,j)*zz1o(i,j+1) + zvf(i,j-1)*zz1o(i,j-1))      &
                  + sorpai*zz1o(i,j)
          ENDDO
       ENDDO

       if ( icontro.NE.0)  THEN      
          do i=1,ie
             do j=1,je
                z1o(i,j)=zz1o(i+ihalo-1,j+ihalo-1)
             enddo
          enddo
          DO j=jb,je-1
             is = MOD(p_ioff+p_joff+j+1,2)+2
             DO i=is,ie-1,2
                zurr=zurr+(zalt(i,j)-z1o(i,j))**2
             enddo
          enddo
       endif


! ATTN : these can be reduced for larger halos than 2    
       if ( mod(kiter,1).eq.0) then 

          do i=1,ie
             do j=1,je
                z1o(i,j)=zz1o(i+ihalo-1,j+ihalo-1)
             enddo
          enddo

          CALL bounds_exch('p',z1o,'mo_tro 10')
          call sethaloN(z1o,zz1o,ihalo)

       endif


    ENDDO   !kiter
    IF (icontro.NE.0) THEN
    WRITE(0,*)'zurr',zurr,'sorpar',sorpar
    ENDIF

  END SUBROUTINE troneu2



  SUBROUTINE OCVTRO

      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_UNITS

      REAL :: DTHI,DTH

      INTEGER :: I,J,K

      DTHI = 0.5 * DTI
      DTH  = 0.5 * DT
!
      DO 1 J=1,JE
      DO 1 I=1,IE
      UZO(I,J)=UZO(I,J)*DEUTIO(I,J)
      VZE(I,J)=VZE(I,J)*DEUTIE(I,J)
!
1     CONTINUE
!      WRITE(IO_STDOUT,*) 'INE OCVTR ',UZO(5,5),VZE(5,5),UCOR(5,5),VCOR(5,5)
!

       DO 11 J=2,JE1
       DO 11 I=2,IE1
       USO(I,J)=UZO(I,J)                                                &
     & +G*STABN*DT*(Z1O(I,J)-Z1O(I+1,J))                       &
     & *AMSUO(I,J,1)/dlxu(i,j)
       VSE(I,J)=VZE(I,J)                                                &
     & +G*STABN*Dt*(Z1O(I,J+1)-Z1O(I,J))                         &
     & *AMSUE(I,J,1)/dlyv(i,j)
11     CONTINUE
!
!#ifdef bounds_exch_save
       CALL bounds_exch('u',USO,'ocvtro 1')
       CALL bounds_exch('v',VSE,'ocvtro 2')
!#endif
!
3001   CONTINUE
!
!      WRITE(IO_STDOUT,*) 'TROPVEL ',USO(5,5),VSE(5,5)
6881  FORMAT(6E20.12)  
!
      K=1
      DO 7715 J=2,JE1
      DO 7715 I=2,IE1
!      ZO(I,J)=ZO(I,J)+Z1O(I,J)  ! ATTN: moved to sbr update_zo
!
!
7715  CONTINUE
!
!#ifdef bounds_exch_save
      CALL bounds_exch('p',Z1O,'ocvtro 3')
!#endif
      
      END SUBROUTINE OCVTRO



  SUBROUTINE update_zo

   USE MO_PARAM1 ,ONLY: ie,je 
   USE MO_COMMO1, ONLY: ZO, Z1O, surdis, dt
#ifdef PBGC
   USE mo_param1_bgc, ONLY: nocetra, ntraad
   USE mo_carbch, ONLY: ocetra
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef ADDCONTRA
   USE mo_param1_add, ONLY: nocectra, nctraad
   USE mo_contra, ONLY: ocectra
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   INTEGER :: i,j,l

#ifdef PBGC
      CALL dilcor_gtrf2
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef ADDCONTRA
      CALL dilcor_gtrf2
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      DO i=1,ie
         DO j=1,je 
            ZO(i,j)=Zo(i,j)+Z1O(I,J)
!             ZO(i,j)=Zo(i,j)+surdis(I,J)*dt
         ENDDO
      ENDDO

#ifdef PBGC
      DO l=ntraad+1,nocetra
      CALL dilcor_ptrf2(ocetra(1,1,1,l))
      ENDDO

#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef ADDCONTRA
      DO l=nctraad+1,nocectra
      CALL dilcor_ptrf2(ocectra(1,1,1,l))
      ENDDO

#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  END SUBROUTINE update_zo

SUBROUTINE alloc_mem_dilcor

      USE MO_PARAM1, ONLY: ie,je

      ALLOCATE(susalo(ie,je),suwath(ie,je))
      ALLOCATE(dilcor(ie,je))


END SUBROUTINE alloc_mem_dilcor

SUBROUTINE dilcor_gtrf

      USE mo_param1,ONLY: ie,je
      USE mo_commo1,ONLY: sao,zo,ddpo,sictho        &
      ,sicsno
      USE MO_COMMOAU1,only:rhoicwa,rhosnwa,sice

      INTEGER :: i,j

!  WRITE(0,*)'in dilcor_pre',size(susalo(:,1)),size(susalo(1,:))

      DO J=1,JE
        DO I=1,IE
           susalo(i,j)=(SAO(I,J,1)*(ZO(I,J)+ddpo(i,j,1)             &
                -SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)           &
                +SICE*SICTHO(I,J))
           suwath(i,j)=ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa      &
                -sicsno(i,j)*rhosnwa
           
        ENDDO
     ENDDO

END SUBROUTINE dilcor_gtrf

SUBROUTINE dilcor_gtrf2

      USE mo_param1,ONLY: ie,je
      USE mo_commo1,ONLY: zo,ddpo,sictho        &
      ,sicsno
      USE MO_COMMOAU1,only:rhoicwa,rhosnwa
      INTEGER :: i,j

!  WRITE(0,*)'in dilcor_pre',size(susalo(:,1)),size(susalo(1,:))

      DO J=1,JE
        DO I=1,IE
           suwath(i,j)=ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa      &
                -sicsno(i,j)*rhosnwa
        ENDDO
     ENDDO

END SUBROUTINE dilcor_gtrf2

SUBROUTINE dilcor_ptrf

  USE mo_param1,ONLY: ie,je
  USE mo_commo1,ONLY: sao,zo,ddpo,sictho        &
       ,sicsno,weto
  USE MO_COMMOAU1,only:rhoicwa,rhosnwa,sice

  INTEGER :: i,j,l
  REAL :: wathne

!  WRITE(0,*)'in dilcor_post',size(dilcor(:,1)),size(dilcor(1,:))

  dilcor(:,:)=1.
  DO j=1,je
     DO i=1,ie
        IF(weto(i,j,1).GT.0.5)THEN
           
           wathne=zo(i,j)+ddpo(i,j,1)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa

           sao(i,j,1)=(susalo(i,j)-sice*sictho(i,j))/wathne
           
        ENDIF        
     ENDDO
  ENDDO

END SUBROUTINE dilcor_ptrf

SUBROUTINE dilcor_ptrf2(TRF)

  USE mo_param1,ONLY: ie,je
  USE mo_commo1,ONLY: zo,ddpo,sictho        &
       ,sicsno,weto
  USE MO_COMMOAU1,only:rhoicwa,rhosnwa
  INTEGER :: i,j
  REAL :: wathne
  REAL :: trf(ie,je,1)

!  WRITE(0,*)'in dilcor_post',size(dilcor(:,1)),size(dilcor(1,:))

  dilcor(:,:)=1.
  DO j=1,je
     DO i=1,ie
        IF(weto(i,j,1).GT.0.5)THEN
           wathne=zo(i,j)+ddpo(i,j,1)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa
           dilcor(i,j)=suwath(i,j)/wathne
        ENDIF
     ENDDO
  ENDDO

  DO j=1,je
     DO i=1,ie
        trf(i,j,1)=trf(i,j,1)*dilcor(i,j)
     ENDDO
  ENDDO

END SUBROUTINE dilcor_ptrf2

SUBROUTINE correct_zo

      USE MO_PARAM1
      USE MO_COMMO1 
      USE MO_PARALLEL
      USE MO_UNITS
#ifdef PBGC
      USE MO_PARAM1_BGC,only: nocetra
      USE MO_CARBCH,only: ocetra      
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef ADDCONTRA
      USE MO_PARAM1_ADD, only: nocectra
      USE MO_CONTRA, only: ocectra
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      INTEGER :: I,J,L

      REAL,ALLOCATABLE  :: ZO_g(:,:), Z1O_g(:,:)

      REAL SUGG, SUZZ, SUZZ1, ZQQ, ZQQ1

      ALLOCATE(ZO_g(ie_g,je_g), Z1O_g(ie_g,je_g))

#ifdef PBGC
      call dilcor_gtrf2
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef ADDCONTRA
      call dilcor_gtrf2
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! SET SEA LEVEL (ZO and Z1O) BACK TO GLOBAL ZERO. Done on one processor
! to ensure reproducibility.


      CALL gather_arr(ZO,ZO_g,p_io)
      CALL gather_arr(Z1O,Z1O_g,p_io)

      IF (p_pe==p_io) THEN
         SUGG=0.0
         SUZZ=0.0
         SUZZ1=0.0
         DO  J=1,je_g
            DO  I=2,ie_g-1
               IF(WETO_g(I,J,1).GT. 0.5) THEN
                  IF(ZO_g(I,J).LT.-3.) THEN
                     WRITE(IO_STDOUT,*) 'ZO LT -3! at i=',i,'j=',j
                  ENDIF
                  SUGG=SUGG+DLXP_g(I,J)*DLYP_g(I,J)
                  SUZZ=SUZZ+DLXP_g(I,J)*DLYP_g(I,J)*ZO_g(I,J)
                  SUZZ1=SUZZ1+DLXP_g(I,J)*DLYP_g(I,J)*Z1O_g(I,J)
               ENDIF
            ENDDO
         ENDDO
         ZQQ=SUZZ/SUGG
         ZQQ1=SUZZ1/SUGG
         WRITE(IO_STDOUT,*) ' MEAN ZETA ', ZQQ,ZQQ1
         DO J=1,je_g
            DO I=1,ie_g
               ZO_g(I,J)=(ZO_g(I,J)-ZQQ)*WETO_g(I,J,1)
               Z1O_g(I,J)=(Z1O_g(I,J)-ZQQ1)*WETO_g(I,J,1)
            ENDDO
         ENDDO
      END IF
      CALL scatter_arr(ZO_g,ZO,p_io)
      CALL scatter_arr(Z1O_g,Z1O,p_io)
      CALL p_bcast(ZQQ,p_io)
      CALL bounds_exch('p',ZO,'mpiom 32')
      CALL bounds_exch('p',Z1O,'mpiom 33')

#ifdef PBGC
      do l=1,nocetra
         call dilcor_ptrf2(ocetra(1,1,1,l))
      enddo
#endif/*PBGC*/             

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef ADDCONTRA
      do l=1,nocectra
         call dilcor_ptrf2(ocectra(1,1,1,l))
      enddo
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE correct_zo

END MODULE mo_tro






