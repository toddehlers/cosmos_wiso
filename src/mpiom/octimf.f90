      SUBROUTINE OCTIMF
      USE MO_PARAM1
      USE MO_COMMO1
      DO  K=1,KE
      DO  J=1,JE
      DO  I=1,IE
       VKE(I,J,K)=VOE(I,J,K)
       UKO(I,J,K)=UOO(I,J,K)
      ENDDO
      ENDDO
      ENDDO
      RETURN
      END
