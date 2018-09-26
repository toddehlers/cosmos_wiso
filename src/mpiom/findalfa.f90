      SUBROUTINE FINDALFA
!:: SBR TO CALCULATE SITES OF POSSIBLE
!:: SLOPE CONVECTION PRIOR TO TIME-STEPPING
!:: JHJ SEPT. 2, 1999
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMOBBL
      USE MO_UNITS    

      ALPX(:,:) = 0
      ALPY(:,:) = 0

!:: SWEEP IN X DIR
      ISZX=0
      DO J=2,JE1
       DO I=2,IE1
        ALP= WETO(I,J,1)*WETO(I+1,J,1)*                                 &
     &        (DEPTO(I+1,J)-DEPTO(I,J))/DLXU(I,J)
        IF (ABS(ALP) .GT. 1.E-8) THEN
        ISZX=ISZX+1
        ALPX(I,J)=ALP
        ENDIF
       ENDDO
      ENDDO
      CALL global_sum(ISZX)
      WRITE(IO_STDOUT,*)'MAX NUMBER OF SLOPE SITES X = ',ISZX
!:: SWEEP IN Y DIR
      ISZY=0
      DO J=2,JE1
       DO I=2,IE1
        ALP= WETO(I,J,1)*WETO(I,J+1,1)*                                 &
     &        (DEPTO(I,J)-DEPTO(I,J+1))/DLYV(I,J)
        IF (ABS(ALP) .GT. 1.E-8) THEN
          ISZY=ISZY+1
          ALPY(I,J)=ALP
        ENDIF
       ENDDO
      ENDDO
      CALL global_sum(ISZY)
      WRITE(IO_STDOUT,*)'MAX NUMBER OF SLOPE SITES Y = ',ISZY

      CALL bounds_exch('u',ALPX,'findalpha 1')
      CALL bounds_exch('v',ALPY,'findalpha 2')

!JJ INITIALIZE BBL THICKNESS ANS FLUX FIELDS (DIAGNOSTIC)
      DO J=1,JE
       DO I=1,IE
        BBLFLX(I,J)=0.0
       ENDDO
      ENDDO
      RETURN
      END
