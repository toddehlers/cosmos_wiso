      SUBROUTINE OCBARP
      USE MO_PARAM1
      USE MO_COMMO1
      USE MO_UNITS
      IMPLICIT NONE
      INTEGER I,J
!
!
       DO J=1,JE
         DO I=1,IE
           UZO(I,J)=USO(I,J)
           VZE(I,J)=VSE(I,J)
         ENDDO
       ENDDO
!
!--------------------------------------------------------------------
!
!     B)
!
!     WINDSTRESS INPUT ON BAROCLINIC VELOCITIES
!
!=====================================================================
!
!     C)
!
       DO 31 J=1,JE
       DO 31 I=1,IE
       VSE(I,J)=VSE(I,J)*AMSUE(I,J,1)
       VZE(I,J)=VZE(I,J)*AMSUE(I,J,1)
       USO(I,J)=USO(I,J)*AMSUO(I,J,1)
       UZO(I,J)=UZO(I,J)*AMSUO(I,J,1)
31     CONTINUE
!
       RETURN
       END
