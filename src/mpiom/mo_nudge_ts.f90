
MODULE mo_nudge_ts

  USE mo_param1, ONLY : ie,je,ke
  USE mo_parallel 
  USE mo_commo1, ONLY : dt, ldtdayc,ldays,lmonts
  USE mo_units

  REAL,ALLOCATABLE :: tecco(:,:,:),tecco_ano(:,:,:)
  REAL,ALLOCATABLE :: secco(:,:,:),secco_ano(:,:,:)
  REAL :: rtim_t
  REAL :: rtim_s


CONTAINS 

SUBROUTINE read_ecco

  INTEGER*8 :: IDATE

  if ( (ldtdayc.eq.1).and.(lmonts.eq.1).and.(ldays.eq.1) ) then

     ALLOCATE(tecco(ie,je,ke),tecco_ano(ie,je,ke)           &
                   ,secco(ie,je,ke),secco_ano(ie,je,ke))

     OPEN(68,FILE='CLIMTEM',STATUS='UNKNOWN',                 &
                    ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
     OPEN(69,FILE='ECCOTEMANO',STATUS='UNKNOWN',                 &
                       ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
     
     OPEN(66,FILE='CLIMSAL',STATUS='UNKNOWN',                 &
                    ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
     OPEN(67,FILE='ECCOSALANO',STATUS='UNKNOWN',                 &
                    ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  endif

  if ( ldtdayc.eq.1 .and.ldays.eq.1 ) then
     do k=1,ke
        IF(p_pe==p_io) READ(68) idate
        CALL read_slice(68,tecco(:,:,k))

        IF(p_pe==p_io) READ(69) idate
        if (k.eq.1) write(0,*)'read_ecco T',idate
        CALL read_slice(69,tecco_ano(:,:,k))

        IF(p_pe==p_io) READ(66) idate
        CALL read_slice(66,secco(:,:,k))

        IF(p_pe==p_io) READ(67) idate
        if (k.eq.1) write(0,*)'read_ecco S',idate
        CALL read_slice(67,secco_ano(:,:,k))
     enddo
  endif 

END SUBROUTINE read_ecco


SUBROUTINE nudge_t(trf,rtim_t)

  IMPLICIT NONE

  INTEGER :: i,j,k
  REAL :: rtim_t
  REAL    :: tt,trf(ie,je,ke)

      IF (rtim_t.GT. 0) THEN 
         if (p_pe==p_io) then
          WRITE(0,*)' TEMP NUDGE TIME [DAYS] =',1./(rtim_t*24.*3600.)
       endif
      ELSE
         if (p_pe==p_io) then
            WRITE(IO_STDOUT,*) 'TEMP NUDGING switched off !!'
         endif
      ENDIF


!$omp do
    DO k=2,ke
      DO j=2,je-1
        DO i=2,ie-1
               tt=tecco(i,j,k)+tecco_ano(i,j,k)
               trf(i,j,k)=trf(i,j,k)+dt*rtim_t*(tt-trf(i,j,k))
        END DO
      END DO
    END DO
!$omp end do


!$omp single
  CALL bounds_exch('p',trf,'nudge_t 3')
!$omp end single

END SUBROUTINE nudge_t


SUBROUTINE nudge_s(trf,rtim_s)

  IMPLICIT NONE

  INTEGER :: i,j,k
  REAL :: rtim_s
  REAL    :: tt,trf(ie,je,ke)

      IF (rtim_s.GT. 0) THEN 
         if (p_pe==p_io) then
          WRITE(0,*)' SAL NUDGE TIME [DAYS] =',1./(rtim_s*24.*3600.)
       endif
      ELSE
         if (p_pe==p_io) then
            WRITE(IO_STDOUT,*) 'SAL NUDGING switched off !!'
         endif
      ENDIF


!$omp do
    DO k=2,ke
      DO j=2,je-1
        DO i=2,ie-1
               tt=secco(i,j,k)+secco_ano(i,j,k)
               trf(i,j,k)=trf(i,j,k)+dt*rtim_s*(tt-trf(i,j,k))
        END DO
      END DO
    END DO
!$omp end do


!$omp single
  CALL bounds_exch('p',trf,'nudge_s 3')
!$omp end single

END SUBROUTINE nudge_s



END MODULE mo_nudge_ts
