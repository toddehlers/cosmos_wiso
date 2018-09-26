SUBROUTINE hdiff

  ! Description:
  !
  ! Horizontal diffusion.
  !
  ! Method:
  !
  ! This subroutine performs horizontal diffusion.
  !
  ! *hdiff* is called from *stepon*
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, June 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, August 1998, tendency diagnostics
  ! T. Diehl, DKRZ, July 1999, parallel version 
  ! I. Kirchner, MPI, December 2000, time control
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_memory_sp,     ONLY: sd, stp, svo
  USE mo_control,       ONLY: lmidatm, ltdiag, nk, nkp1, nlev, nlevp1, &
                              ndiahdf
  USE mo_truncation,    ONLY: ntrn
  USE mo_semi_impl,     ONLY: hdamp, vcrit, vmax
  USE mo_hdiff,         ONLY: damhih, difd, dift, diftcor, difvo,      &
                              enstdif, ldiahdf, nlvstd1, nlvstd2
  USE mo_constants,     ONLY: a
  USE mo_diff,          ONLY: iq, ncdif
  USE mo_diag_tendency, ONLY: pdvor, pdtem, pddiv
  USE mo_time_control,  ONLY: time_step_len, delta_time, get_time_step

  IMPLICIT NONE

  !  Local scalars: 
  REAL(dp):: zaa, zcons, zcor, zdamp, zdifd, zdift, zdifvo, znn        &
           , ztwodt, zzcor, zztemp, zzzn
  INTEGER :: ilev2, in, is, jk, jlev, jn, jndif, jr, ic, i, snsp, nns

  !  Local arrays: 
  REAL(dp):: zcor1(nkp1,nlev), zcor2(nkp1,nlev), znd(nkp1,2*nlev)      &
           , zndifd(nkp1,nlev), zndift(nkp1,nlev), zndifvo(nkp1,nlev)  &
           , znt(nkp1,nlev), znvo(nkp1,2*nlev), zzus(nkp1,nlev)        &
           , zhigh(nkp1,nlev)

  INTEGER :: np1(lc%snsp), nindex(lc%nns)
  INTEGER :: istep
  !  Executable statements 

  istep = get_time_step()

  snsp   = lc%snsp
  nns    = lc%nns
  np1    = lc%np1
  nindex = lc%nindex

  ! Each PE loops only over n's owned by itself 

  do i=1,nns
     jn=nindex(i)
     zzus(jn,:) = 0._dp
  enddo

!-- 1. Compute diffusion correction factor

  ztwodt = time_step_len
  ilev2  = 2*nlev
  zaa    = 1._dp/(a*a)
  zcons  = 2.5_dp*delta_time/a

  do i=1,nns
     jn=nindex(i)
     if (jn >= 2) then
        in = jn - 1
        zdamp = 1._dp
        IF (lmidatm) THEN
           DO jk = nlev, 1, -1
              IF (in>ntrn(jk)) zdamp = zdamp*hdamp
              IF (jk<=nlvstd2 .AND. jk>=nlvstd1)    zdamp=zdamp*enstdif
              zzus(jn,jk) = zdamp
           END DO
        ELSE
           IF (jn<=nk/3) THEN
              DO jk = nlev, 1, -1
                 IF (in>ntrn(jk)) zdamp = zdamp*hdamp
                 IF (jk<=nlvstd2 .AND. jk>=nlvstd1) zdamp=zdamp*enstdif
                 zzus(jn,jk) = zdamp
              END DO
           ELSE
              DO jk = 1, nlev
                 zzus(jn,jk) = 1._dp
              END DO
           END IF
        END IF
     endif
  END DO

  IF (lmidatm) THEN
     do i=1,nns
        jn=nindex(i)
        DO jk = 1,nlev
           zhigh(jn,jk)=1._dp
        END DO
     END DO

     IF (.NOT.ldiahdf) THEN
        do i=1,nns
           jn=nindex(i)
           if (jn >= 3) then
              znn = jn-1
              DO jk = 1,nlev
                 zcor = (1._dp+(znn*vmax(jk)-vcrit)) 
                 IF (zcor .GT. 1._dp) zhigh(jn,jk) = damhih
              END DO
           endif
        END DO
     ELSE
      ! same loop but WITH statistics
        do i=1,nns
           jn=nindex(i)
           if (jn >= 3) then
              znn = jn-1
              DO jk = 1,nlev
                 zcor = (1._dp+(znn*vmax(jk)-vcrit)) 
                 IF (zcor .GT. 1._dp) THEN
                    zhigh(jn,jk) = damhih
                    WRITE(ndiahdf,*) istep, jk, jn-1, vmax(jk)   
                 ENDIF
              END DO
           endif
        END DO
     END IF
  END IF

  ! Vertical loop

  DO jk = 1, nlev

     jndif = ncdif(jk)
     zztemp = (nkp1*nk*zaa)**(1-iq(jk))*ztwodt
     zdifd = difd*zztemp
     zdifvo = difvo*zztemp
     zdift = dift*zztemp
     
     do i=1,nns
        jn=nindex(i)
        zndifd(jn,jk) = 1._dp
        zndifvo(jn,jk) = 1._dp
        zndift(jn,jk) = 1._dp
     END DO
     
     do i=1,nns
        jn=nindex(i)
        if (jn >= jndif+1) then
           zzzn = jn - 1 - jndif
           zztemp = (zaa*zzzn*(zzzn+1._dp))**iq(jk)*zzus(jn,jk)
           IF (lmidatm) THEN
              zndifd(jn,jk)  = 1._dp + zhigh(jn,jk)*zdifd*zztemp
              zndifvo(jn,jk) = 1._dp + zhigh(jn,jk)*zdifvo*zztemp
              zndift(jn,jk)  = 1._dp + zhigh(jn,jk)*zdift*zztemp
           ELSE
              zndifd(jn,jk)  = 1._dp + zdifd*zztemp
              zndifvo(jn,jk) = 1._dp + zdifvo*zztemp
              zndift(jn,jk)  = 1._dp + zdift*zztemp
           ENDIF
        endif
     END DO
     
     IF ( .NOT. lmidatm) THEN
        IF ( .NOT. ldiahdf) THEN
           do i=1,nns
              jn=nindex(i)
              if (jn >= 3) then
                 znn = jn - 1
                 zcor = (1._dp+zcons*(znn*vmax(jk)-vcrit))
                 IF (zcor>1._dp) THEN
                    zndifd(jn,jk) = zcor*zndifd(jn,jk)
                    zndifvo(jn,jk) = zcor*zndifvo(jn,jk)
                 END IF
              endif
           END DO
        ELSE
           
           ! Same loop but with statistics
           
           do i=1,nns
              jn=nindex(i)
              if (jn >= 3) then
                 znn = jn - 1
                 zcor = (1._dp+zcons*(znn*vmax(jk)-vcrit))
                 IF (zcor>1._dp) THEN
                    WRITE (ndiahdf,*) istep, jk, jn - 1, zcor, zndifd(jn,jk)
                    zndifd(jn,jk) = zcor*zndifd(jn,jk)
                    zndifvo(jn,jk) = zcor*zndifvo(jn,jk)
                 END IF
              endif
           END DO
        END IF
     END IF
     
     do i=1,nns
        jn=nindex(i)
        if (jn >= 3) then
           znd(jn,jk) = 1._dp/zndifd(jn,jk)
           znd(jn,jk+nlev) = znd(jn,jk)
           znvo(jn,jk) = 1._dp/zndifvo(jn,jk)
           znvo(jn,jk+nlev) = znvo(jn,jk)
        else
           znd(jn,jk) = 1._dp
           znd(jn,jk+nlev) = 1._dp
           znvo(jn,jk) = 1._dp
           znvo(jn,jk+nlev) = 1._dp
        endif
     END DO
     
     ! Temperature
     
     do i=1,nns
        jn=nindex(i)
        if (jn >= 2) then
           zzcor = 1._dp
           IF (jn-1>ntrn(jk)) zzcor = 0._dp
           znn = jn - 1._dp
           IF (lmidatm) THEN
              zcor = 1._dp
           ELSE
              zcor = (1._dp+zcons*(znn*vmax(jk)-vcrit))
              IF (zcor>1._dp) zndift(jn,jk) = zcor*zndift(jn,jk)
              IF (zcor<=1._dp) zcor = 1._dp
           END IF
           zcor2(jn,jk) = zzcor*diftcor(jk)
           zcor1(jn,jk) = zcor2(jn,jk)/zcor
        else
           zcor1(jn,jk) = 0._dp
           zcor2(jn,jk) = 0._dp
        endif
     END DO
     
     do i=1,nns
        jn=nindex(i)
        znt(jn,jk) = 1._dp/zndift(jn,jk)
     END DO
  END DO
  
!-- 2. Modify fields

  DO is = 1, snsp

     ic = np1(is)
 
      IF (ltdiag) THEN
         ! remove vorticity and divergence without diffusion part
         pddiv(:,:,is,8) = pddiv(:,:,is,8) - sd (:,:,is)
         pdvor(:,:,is,8) = pdvor(:,:,is,8) - svo(:,:,is)
      ENDIF

      sd (:,1,is) = sd (:,1,is)*znd (ic,:nlev)
      sd (:,2,is) = sd (:,2,is)*znd (ic, nlev+1:)
      svo(:,1,is) = svo(:,1,is)*znvo(ic,:nlev)
      svo(:,2,is) = svo(:,2,is)*znvo(ic, nlev+1:)
 
      IF (ltdiag) THEN
         ! store diffusion part of vorticity and divergence
         pddiv(:,:,is,8) = pddiv(:,:,is,8) + sd (:,:,is)
         pdvor(:,:,is,8) = pdvor(:,:,is,8) + svo(:,:,is)
         ! remove temperature without diffusion part
         pdtem(:,:,is,11) = pdtem(:,:,is,11) - stp(1:nlev,:,is)
      ENDIF

      DO jr = 1, 2

!DIR$ IVDEP
!CDIR NODEP
!OCL NOVREC,NOALIAS

        DO jlev = 1, nlev
          stp(jlev,jr,is) = stp(nlevp1,jr,is)*zcor1(ic,jlev) +         &
               (stp(jlev,jr,is)-zcor2(ic,jlev)*stp(nlevp1,jr,is))*     &
               znt(ic,jlev)
        END DO

      END DO
      ! add temperature with diffusion part
      IF (ltdiag) pdtem(:,:,is,11) = pdtem(:,:,is,11) + stp(1:nlev,:,is)

   END DO

  RETURN
END SUBROUTINE hdiff
