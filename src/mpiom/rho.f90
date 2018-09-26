FUNCTION RHO(S,T,P)
  USE mo_kind, ONLY: dp

!*********************************************************************
!
!
!     RRRRR   H    H   OOO
!     R    R  H    H  O   O
!     RRRRR   HHHHHH  O   O
!     R  RR   H    H  O   O
!     R   RR  H    H   OOO
!
!*****************************************************************
! WIRD FUER DEN REFERENZ-ZUSTAND VERWENDET
!++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL(dp), PARAMETER :: a0 =  999.842594_dp,     &
                         a1 =  6.793952e-2_dp,  &
                         a2 = -9.095290e-3_dp,  &
                         a3 =  1.001685e-4_dp,  &
                         a4 = -1.120083e-6_dp,  &
                         a5 =  6.536332e-9_dp
  REAL(dp), PARAMETER :: b0 =  8.24493e-1_dp,   &
                         b1 = -4.0899e-3_dp,    &
                         b2 =  7.6438e-5_dp,    &
                         b3 = -8.2467e-7_dp,    &
                         b4 =  5.3875e-9_dp
  REAL(dp), PARAMETER :: c0 = -5.72466e-3_dp,   &
                         c1 =  1.0227e-4_dp,    &
                         c2 = -1.6546e-6_dp
  REAL(dp), PARAMETER :: d0 =  4.8314e-4_dp
  REAL(dp), PARAMETER :: e0 =  19652.21_dp,     &
                         e1 =  148.4206_dp,     &
                         e2 = -2.327105_dp,     &
                         e3 =  1.360477e-2_dp,  &
                         e4 = -5.155288e-5_dp
  REAL(dp), PARAMETER :: f0 =  54.6746_dp,      &
                         f1 = -0.603459_dp,     &
                         f2 =  1.09987e-2_dp,   &
                         f3 = -6.1670e-5_dp
  REAL(dp), PARAMETER :: g0 =  7.944e-2_dp,     &
                         g1 =  1.6483e-2_dp,    &
                         g2 = -5.3009e-4_dp
  REAL(dp), PARAMETER :: h0 =  3.239908_dp,     &
                         h1 =  1.43713e-3_dp,   &
                         h2 =  1.16092e-4_dp,   &
                         h3 = -5.77905e-7_dp
  REAL(dp), PARAMETER :: ai0 =  2.2838e-3_dp,   &
                         ai1 = -1.0981e-5_dp,   &
                         ai2 = -1.6078e-6_dp
  REAL(dp), PARAMETER :: aj0 =  1.91075e-4_dp
  REAL(dp), PARAMETER :: ak0 =  8.50935e-5_dp,  &
                         ak1 = -6.12293e-6_dp,  &
                         ak2 =  5.2787e-8_dp
  REAL(dp), PARAMETER :: am0 = -9.9348e-7_dp,   &
                         am1 =  2.0816e-8_dp,   &
                         am2 =  9.1697e-10_dp


      S3H=SQRT(S**3)
      RHOW=A0+T*(A1+T*(A2+T*(A3+T*(A4+T*A5))))
      AKW=E0+T*(E1+T*(E2+T*(E3+T*E4)))
      AW=H0+T*(H1+T*(H2+T*H3))
      BW=AK0+T*(AK1+T*AK2)
      B=BW+S*(AM0+T*(AM1+T*AM2))
      A=AW+S*(AI0+T*(AI1+AI2*T))+AJ0*S3H
      AKST0=AKW+S*(F0+T*(F1+T*(F2+T*F3)))+S3H*(G0+T*(G1+G2*T))
      AKSTP=AKST0+P*(A+B*P)
      RHST0=RHOW+S*(B0+T*(B1+T*(B2+T*(B3+T*B4))))+D0*S**2               &
     &+S3H*(C0+T*(C1+C2*T))
      RHO=RHST0/(1.-P/AKSTP)

      RETURN
END FUNCTION RHO
