SUBROUTINE relax_ts

  USE mo_param1
  USE mo_parallel
  USE mo_commo1
  USE mo_levitus
  USE mo_units

  IMPLICIT NONE

  INTEGER i, j, k

  REAL fakres

  IF (I3DREST .GT. 0) THEN  

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('p',tlevi,'relax_ts 1')
  CALL bounds_exch('p',slevi,'relax_ts 2')
!$OMP END SINGLE
!#endif

    fakres=dt/(86400.*30.*4.)
!$OMP SINGLE
    WRITE(io_stdout,*)'3D-RESTORING!!!'
!$OMP END SINGLE
    DO k=4,ke
!$OMP DO
      DO j=2,je-1
        DO i=2,ie-1
          IF(ltlev(i,j,k))                                                 &
               tho(i,j,k)=tho(i,j,k)+fakres*(tlevi(i,j,k)-tho(i,j,k))
          IF(lslev(i,j,k))                                                 &
               sao(i,j,k)=sao(i,j,k)+fakres*(slevi(i,j,k)-sao(i,j,k))
        END DO
      END DO
!$OMP END DO
    END DO
  ENDIF

!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('p',tho,'relax_ts 3')
  CALL bounds_exch('p',sao,'relax_ts 4')
!$OMP END SINGLE
!#endif
  
END SUBROUTINE relax_ts

