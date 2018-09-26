SUBROUTINE SWUVO3 &
  &( KIDIA,KFDIA,KLON,KNU,KABS &
  &, PU, PTR &
  & )
  
!**** *SWUVO3* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS

!     PURPOSE.
!     --------
!           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR OZONE
!     IN THE UV and VISIBLE SPECTRAL INTERVALS.

!**   INTERFACE.
!     ----------
!          *SWUVO3* IS CALLED FROM *SW1S*.


!        EXPLICIT ARGUMENTS :
!        --------------------
! KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
! KABS   :                     ; NUMBER OF ABSORBERS
! PU     : (KLON,KABS)         ; ABSORBER AMOUNT
!     ==== OUTPUTS ===
! PTR    : (KLON,KABS)         ; TRANSMISSION FUNCTION

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          TRANSMISSION FUNCTION ARE COMPUTED USING SUMS OF EXPONENTIALS

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------
!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        Ph.DUBUISSON/B.BONNEL L.O.A.

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 00-12-18
   
!-----------------------------------------------------------------------

USE MO_KIND  , ONLY : DP

USE MO_SW    , ONLY : NEXPO3, REXPO3


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER :: KABS
INTEGER :: KFDIA
INTEGER :: KIDIA
INTEGER :: KLON
INTEGER :: KNU

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL(DP) :: PU(KLON,KABS)
REAL(DP) :: PTR(KLON,KABS)

!-----------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

!     LOCAL INTEGER SCALARS
INTEGER :: JA, JL, IEXP, JX


IEXP=NEXPO3(KNU)
!print *,'IEXP(',KNU,')=',IEXP
!print *,(REXPO3(KNU,1,JX),JX=1,IEXP)
!print *,(REXPO3(KNU,2,JX),JX=1,IEXP)

DO JA = 1,KABS
  DO JL=KIDIA,KFDIA
    PTR(JL,JA)=0._DP
  END DO
    
  DO JX=1,IEXP
    DO JL = KIDIA,KFDIA
      PTR(JL,JA) = PTR(JL,JA) &
        &+REXPO3(KNU,1,JX)*EXP(-REXPO3(KNU,2,JX)*PU(JL,JA))
    END DO
  END DO    
ENDDO

RETURN
END SUBROUTINE SWUVO3
