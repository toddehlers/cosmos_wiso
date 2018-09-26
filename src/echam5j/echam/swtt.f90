SUBROUTINE SWTT ( KIDIA, KFDIA, KBDIM, KNU, KA , PU, PTR)

!**** *SWTT* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS

!     PURPOSE.
!     --------
!           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
!     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
!     INTERVALS.

!**   INTERFACE.
!     ----------
!          *SWTT* IS CALLED FROM *SW1S*, *SWNI*.


!        EXPLICIT ARGUMENTS :
!        --------------------
! KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
! KA     :                     ; INDEX OF THE ABSORBER
! PU     : (KBDIM)             ; ABSORBER AMOUNT
!     ==== OUTPUTS ===
! PTR    : (KBDIM)             ; TRANSMISSION FUNCTION

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
!        ORIGINAL : 88-12-15
   
!-----------------------------------------------------------------------

USE MO_KIND  , ONLY : DP

USE MO_SW    , ONLY : APAD     ,BPAD     ,D


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER :: KA
INTEGER :: KFDIA
INTEGER :: KIDIA
INTEGER :: KBDIM
INTEGER :: KNU



!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL(DP):: PU(KBDIM), PTR(KBDIM)

!-----------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL(DP):: ZR1(KBDIM), ZR2(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: JL


!-----------------------------------------------------------------------

!*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION


DO JL = KIDIA,KFDIA
  ZR1(JL) = APAD(KNU,KA,1) + PU(JL) * (APAD(KNU,KA,2) + PU(JL)&
   &* ( APAD(KNU,KA,3) + PU(JL) * (APAD(KNU,KA,4) + PU(JL)&
   &* ( APAD(KNU,KA,5) + PU(JL) * (APAD(KNU,KA,6) + PU(JL)&
   &* ( APAD(KNU,KA,7) ))))))

  ZR2(JL) = BPAD(KNU,KA,1) + PU(JL) * (BPAD(KNU,KA,2) + PU(JL)&
   &* ( BPAD(KNU,KA,3) + PU(JL) * (BPAD(KNU,KA,4) + PU(JL)&
   &* ( BPAD(KNU,KA,5) + PU(JL) * (BPAD(KNU,KA,6) + PU(JL)&
   &* ( BPAD(KNU,KA,7) ))))))


!*         2.      ADD THE BACKGROUND TRANSMISSION



  PTR(JL) = (ZR1(JL) / ZR2(JL)) * (1._DP - D(KNU,KA)) + D(KNU,KA)
ENDDO

RETURN
END SUBROUTINE SWTT
