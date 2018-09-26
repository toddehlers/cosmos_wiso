      SUBROUTINE WRTE_MEANADD

!****************************************************************
!     save add output
!**********************************************************************

      USE mo_param1, ONLY: ie,je,ke
      USE MO_COMMO1, ONLY: ldays,lmonts,ldtdayc,dt,ddpo,monlen
      USE mo_addmean, ONLY: mean_3d_freq_add,  meantime_2d_add, meantime_3d_add

      IMPLICIT NONE

      INTEGER :: NDTDAY,nanf3d,nend3d,i

      NDTDAY=NINT(86400./DT)
     
      nanf3d=0
     
      nend3d=-1








! for 3d output
  ! daily average
      IF (mean_3d_freq_add.EQ.1) THEN
        nanf3d=ldtdayc
        nend3d=ndtday
! monthly average
     ELSEIF (mean_3d_freq_add.EQ.2) THEN
        IF (ldays.EQ.1) THEN
          nanf3d=ldtdayc
        ELSEIF (ldays.GT.1) THEN
          nanf3d=ldtdayc+((ldays-1)*ndtday)
        ENDIF
        nend3d=ndtday*monlen(lmonts)
! yearly average
      ELSEIF (mean_3d_freq_add.EQ.3) THEN
        IF ((ldays.EQ.1).AND.(lmonts.EQ.1)) THEN
          nanf3d=ldtdayc
         ELSEIF ((ldays.GT.1).AND.(lmonts.EQ.1)) THEN
          nanf3d=ldtdayc+((ldays-1)*ndtday)
         ELSEIF (lmonts.GT.1) THEN
            nanf3d=0
            DO i=1,lmonts-1
               nanf3d=nanf3d+monlen(i)
            ENDDO
            nanf3d=nanf3d*ndtday
            nanf3d=nanf3d+(ldtdayc+((ldays-1)*ndtday))
        ENDIF
        nend3d=0
        DO i=1,12
          nend3d=nend3d+monlen(i)
        ENDDO
        nend3d=nend3d*ndtday
      ELSEIF (mean_3d_freq_add.EQ.4) THEN

! every timestep   
        nanf3d=1
        nend3d=1

! no diagnostic output
       ELSEIF (mean_3d_freq_add.EQ.0) THEN
         CONTINUE
      nanf3d=0
      nend3d=-1

       ELSE
         STOP 'stop in wrte_bgcmean due to wrong parameter.'
      ENDIF



 
    

      IF (nanf3d.EQ.nend3d) THEN
         meantime_2d_add = meantime_2d_add + 1
         meantime_3d_add = meantime_3d_add + 1
         CALL avrg_addmean_3d(ie,je,ke) 
         CALL write_addmean_3d(ie,je,ke,ddpo)


      ENDIF


      END






