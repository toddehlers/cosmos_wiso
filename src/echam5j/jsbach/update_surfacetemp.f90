SUBROUTINE update_surfacetemp(klon, pcp,                             &
  &            pescoe, pfscoe, peqcoe, pfqcoe,                       &
  &            psold, pqsold, pdqsold,                               &
  &            pnetrad, pgrdfl,                                      &
  &            pcfh, pcair, pcsat, pfracsu, pgrdcap,                 &
  &            psnew, pqsnew)


  USE mo_time_control,     ONLY: time_step_len
#ifndef STANDALONE
    USE mo_physc2, ONLY: cvdifts
#endif
  USE mo_radiation,        ONLY: cemiss
  USE mo_constants,        ONLY: stbo, cpd, vtmpc2, alv, als
  USE mo_kind,             ONLY: dp

  IMPLICIT NONE

  INTEGER,  INTENT(in)    :: klon
  REAL(dp),     INTENT(in)    :: pcp(klon), pfscoe(klon), pescoe(klon), pfqcoe(klon), peqcoe(klon)
  REAL(dp),     INTENT(in)    :: psold(klon), pqsold(klon), pdqsold(klon)
  REAL(dp),     INTENT(in)    :: pnetrad(klon), pgrdfl(klon)
  REAL(dp),     INTENT(in)    :: pcfh(klon), pcair(klon), pcsat(klon), pfracsu(klon)
  REAL(dp),     INTENT(in)    :: pgrdcap(klon)
  REAL(dp),     INTENT(out)   :: psnew(klon), pqsnew(klon)
!
  INTEGER :: jl
  REAL(dp) :: zcolin(klon), zcohfl(klon), zcoind(klon), zicp(klon), zca(klon), zcs(klon)
  REAL(dp) :: pdt, pemi, pboltz
  REAL(dp) :: pc16, platev, platsu
  REAL(dp) :: ztpfac1, ztmst

#ifdef STANDALONE
    REAL(dp), PARAMETER :: cvdifts = 1.0_dp
#endif


!-------------------------------------------------------------------------------------
! Constants

  ztpfac1 = cvdifts
  ztmst   = time_step_len
  pdt     = ztpfac1*ztmst                  ! zcons29 in 'old' vdiff
  pemi    = cemiss                         ! emissivity
  pboltz  = stbo
  pc16    = cpd * vtmpc2
  platev  = alv
  platsu  = als

!************************************************************************************
!
     zicp(:) = 1._dp / pcp(:)
!
     zca(:)    = platsu * pfracsu(:) +  platev * (pcair(:) - pfracsu(:))
     zcs(:)    = platsu * pfracsu(:) +  platev * (pcsat(:) - pfracsu(:))
!
     zcolin(:) = pgrdcap(:)*zicp(:) +                                             &
          &      pdt * (zicp(:) * 4._dp * pemi * pboltz * ((zicp(:) * psold(:))**3) -     &
          &      pcfh(:) * (zca(:) * peqcoe(:) - zcs(:) -                      &
          &      zicp(:) * pc16 * psold(:) *                        &
          &      (pcair(:) * peqcoe(:) - pcsat(:)))*        &
          &      zicp(:) * pdqsold(:))
!
     zcohfl(:) = -pdt * pcfh(:) * (pescoe(:) - 1._dp)
!
     zcoind(:) = pdt * (pnetrad(:) + pcfh(:) * pfscoe(:) +  pcfh(:)*          &
          &      ((zca(:) * peqcoe(:) - zcs(:)) * pqsold(:) + zca(:) * pfqcoe(:) -   &
          &      zicp(:) * pc16 * psold(:) *                                &
          &      ((pcair(:) * peqcoe(:) - pcsat(:)) * pqsold(:) +    &
          &      pcair(:) * pfqcoe(:))) + pgrdfl(:))
!
    psnew(:)  = (zcolin(:) * psold(:) + zcoind(:)) / (zcolin(:) + zcohfl(:))
    pqsnew(:) = pqsold(:) + zicp(:) * pdqsold(:) * (psnew(:) - psold(:))

END SUBROUTINE update_surfacetemp
