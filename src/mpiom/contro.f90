      SUBROUTINE CONTRO(L)
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMOAU1
      USE MO_UNITS
      IMPLICIT NONE
      INTEGER l,i,j,k
      REAL sum1, sums(IE,JE)

      if (icontro.ne.0) then
!     Summing up in the following way should give identical results
!     independent of the number of processors
!     (at least if optimization is turned off)

      sums(:,:) = 0
      DO K=1,KE
       DO J=1,JE
        DO I=2,IE-1
         sums(i,j)=sums(i,j)+DLXP(I,J)*DLYP(I,J)*DDPO(I,J,K)*SAO(I,J,K)
        ENDDO
       ENDDO
      ENDDO
      DO J=1,JE
       DO I=2,IE-1
        sums(i,j)=sums(i,j)+DLXP(I,J)*DLYP(I,J)*WETO(I,J,1)*(SAO(I,J,1) &
     &    *(ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)            &
     &    +SICE*SICTHO(I,J))
       ENDDO
      ENDDO

      SUM1 = global_array_sum(sums)

      if (p_pe==p_io) then
      WRITE(IO_STDOUT,*)'SALZCHECK: ',SUM1,L
!      WRITE(IO_STDOUT,*)'SALZCHECK: ',SUM1,L
      endif
      if (icontro.eq.999) then
         CALL WRTE_DEBUG(L)
      endif

      endif

      RETURN
      END
