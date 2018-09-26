SUBROUTINE SWTT1 ( KIDIA,KFDIA,KBDIM,KNU,KABS,KKIND, PU, PTR )

!**** *SWTT1* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS

!     PURPOSE.
!     --------
!           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
!     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
!     INTERVALS.

!**   INTERFACE.
!     ----------
!          *SWTT1* IS CALLED FROM *SW1S*.


!        EXPLICIT ARGUMENTS :
!        --------------------
! KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
! KABS   :                     ; NUMBER OF ABSORBERS
! KKIND  : (KABS)              ; INDICES OF THE ABSORBERS
! PU     : (KBDIM,KABS)         ; ABSORBER AMOUNT
!     ==== OUTPUTS ===
! PTR    : (KBDIM,KABS)         ; TRANSMISSION FUNCTION

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
!     AND HORNER'S ALGORITHM.

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

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 95-01-20
   
!-----------------------------------------------------------------------

USE MO_KIND  , ONLY : DP

USE MO_SW    , ONLY : APAD     ,BPAD     ,D


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER :: KABS
INTEGER :: KFDIA
INTEGER :: KIDIA
INTEGER :: KBDIM
INTEGER :: KNU



!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER :: KKIND(KABS)
REAL(DP):: PU(KBDIM,KABS)
REAL(DP):: PTR(KBDIM,KABS)

!-----------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL(DP):: ZR1(KBDIM), ZR2(KBDIM), ZU(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IA, JA, JL


!-----------------------------------------------------------------------

!*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION


DO JA = 1,KABS
  IA=KKIND(JA)
  DO JL = KIDIA,KFDIA
    ZU(JL) = PU(JL,JA)
    ZR1(JL) = APAD(KNU,IA,1) + ZU(JL) * (APAD(KNU,IA,2) + ZU(JL)&
     &* ( APAD(KNU,IA,3) + ZU(JL) * (APAD(KNU,IA,4) + ZU(JL)&
     &* ( APAD(KNU,IA,5) + ZU(JL) * (APAD(KNU,IA,6) + ZU(JL)&
     &* ( APAD(KNU,IA,7) ))))))

    ZR2(JL) = BPAD(KNU,IA,1) + ZU(JL) * (BPAD(KNU,IA,2) + ZU(JL)&
     &* ( BPAD(KNU,IA,3) + ZU(JL) * (BPAD(KNU,IA,4) + ZU(JL)&
     &* ( BPAD(KNU,IA,5) + ZU(JL) * (BPAD(KNU,IA,6) + ZU(JL)&
     &* ( BPAD(KNU,IA,7) ))))))


!*         2.      ADD THE BACKGROUND TRANSMISSION


    PTR(JL,JA) = (ZR1(JL)/ZR2(JL)) * (1._DP-D(KNU,IA)) + D(KNU,IA)
  ENDDO
ENDDO

RETURN
END SUBROUTINE SWTT1
