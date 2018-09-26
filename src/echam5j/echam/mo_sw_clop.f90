MODULE mo_sw_clop

!- Description:
!
!  This module contains:
!
!  Coefficients for the parameterization of the extinction (tau),
!  single scattering albedo (omega) and asymmetry factor (gamma)
!  of cloud particles. Coefficients are given for cloud liquid droplets
!  and ice crystals.
!
!  Coefficients are given for 2 or 4 band resolution of the SW spectrum.
!  Reference : Rockel et al. (1991).
!
!  2 bands: as in ECHAM4
!           visible : band 1 : 0.25 - 0.68 micrometer
!           near IR : band 2 : 0.68 - 4.0  micrometer
!
!  4 bands: as in ECMWF model
!           visible : band 1 : 0.25 - 0.69 micrometer
!           near IR : band 2 : 0.69 - 1.19 micrometer
!                     band 3 : 1.19 - 2.38 micrometer
!                     band 4 : 2.38 - 4.00 micrometer
!
!  Coefficients: xNBPI, N = 2 or 4 for SW resolution
!                       B = 1,..N band index
!                       P = 1 for liquid, 2 for ice phase
!                       I = coefficient index
!
!  extinction  : x=a , single scattering : x=b , asymmetry factor  : x=c
!
!  Author:
!
!  Marco Giorgetta, MPI, December 1998

   USE mo_kind, ONLY: dp
   IMPLICIT NONE

!  2 bands
!   extinction
!    band 1
!     liquid
      REAL(dp), PARAMETER :: a21l0= 1.8706058417_dp
      REAL(dp), PARAMETER :: a21l1=-1.0756364457_dp
!     ice
      REAL(dp), PARAMETER :: a21i0= 1.9056067246_dp
      REAL(dp), PARAMETER :: a21i1=-1.0318784654_dp
!    band 2
!     liquid
      REAL(dp), PARAMETER :: a22l0= 1.9655460426_dp
      REAL(dp), PARAMETER :: a22l1=-1.0778999732_dp
!     ice
      REAL(dp), PARAMETER :: a22i0= 2.1666771102_dp
      REAL(dp), PARAMETER :: a22i1=-1.0634702711_dp
!   single scattering
!    band 1
!     liquid
      REAL(dp), PARAMETER :: b21l0= 0.99999999_dp
!     ice
      REAL(dp), PARAMETER :: b21i0= 0.99999999_dp
!    band 2
!     liquid
      REAL(dp), PARAMETER :: b22l0= 0.9854369057_dp
      REAL(dp), PARAMETER :: b22l1= 0.013584242533_dp
      REAL(dp), PARAMETER :: b22l2=-0.024856960461_dp
      REAL(dp), PARAMETER :: b22l3= 0.0055918314369_dp
!     ice
      REAL(dp), PARAMETER :: b22i0= 0.98475089485_dp
      REAL(dp), PARAMETER :: b22i1= 0.0053152066002_dp
      REAL(dp), PARAMETER :: b22i2=-0.0061150583857_dp
      REAL(dp), PARAMETER :: b22i3=-0.0032775655896_dp
!   asymmetry factor
!    band 1
!     liquid
      REAL(dp), PARAMETER :: c21l0= 0.78756640717_dp
      REAL(dp), PARAMETER :: c21l1= 0.10660598895_dp
      REAL(dp), PARAMETER :: c21l2=-0.031012468401_dp
!     ice
      REAL(dp), PARAMETER :: c21i0= 0.7700034985_dp
      REAL(dp), PARAMETER :: c21i1= 0.19598466851_dp
      REAL(dp), PARAMETER :: c21i2=-0.11836420885_dp
      REAL(dp), PARAMETER :: c21i3= 0.025209205131_dp
!    band 2
!     liquid
      REAL(dp), PARAMETER :: c22l0= 0.79208639802_dp
      REAL(dp), PARAMETER :: c22l1=-0.044930076174_dp
      REAL(dp), PARAMETER :: c22l2= 0.18980672305_dp
      REAL(dp), PARAMETER :: c22l3=-0.082590933352_dp
!     ice
      REAL(dp), PARAMETER :: c22i0= 0.83631171237_dp
      REAL(dp), PARAMETER :: c22i1=-0.19965998649_dp
      REAL(dp), PARAMETER :: c22i2= 0.46130320487_dp
      REAL(dp), PARAMETER :: c22i3=-0.29719270332_dp
      REAL(dp), PARAMETER :: c22i4= 0.062554483594_dp

!  4 bands
!   extinction
!    band 1
!     liquid
      REAL(dp), PARAMETER :: a41l0= 1.8362_dp
      REAL(dp), PARAMETER :: a41l1=-1.0665_dp
!     ice
      REAL(dp), PARAMETER :: a41i0= 1.9787_dp
      REAL(dp), PARAMETER :: a41i1=-1.0365_dp
!    band 2
!     liquid
      REAL(dp), PARAMETER :: a42l0= 2.0731_dp
      REAL(dp), PARAMETER :: a42l1=-1.1079_dp
!     ice
      REAL(dp), PARAMETER :: a42i0= 2.1818_dp
      REAL(dp), PARAMETER :: a42i1=-1.0611_dp
!    band 3
!     liquid
      REAL(dp), PARAMETER :: a43l0= 1.8672_dp
      REAL(dp), PARAMETER :: a43l1=-1.0420_dp
!     ice
      REAL(dp), PARAMETER :: a43i0= 1.9608_dp
      REAL(dp), PARAMETER :: a43i1=-1.0212_dp
!    band 4
!     liquid
      REAL(dp), PARAMETER :: a44l0= 1.0787_dp
      REAL(dp), PARAMETER :: a44l1=-0.79772_dp
!     ice
      REAL(dp), PARAMETER :: a44i0= 1.2558_dp
      REAL(dp), PARAMETER :: a44i1=-0.88622_dp
!   single scattering
!    band 1
!     liquid
      REAL(dp), PARAMETER :: b41l0= 1._dp
      REAL(dp), PARAMETER :: b41l1=-2.2217e-7_dp
!     ice
      REAL(dp), PARAMETER :: b41i0= 1._dp
      REAL(dp), PARAMETER :: b41i1=-1.143e-7_dp
!    band 2
!     liquid
      REAL(dp), PARAMETER :: b42l0= 1._dp
      REAL(dp), PARAMETER :: b42l1=-1.6712e-5_dp
!     ice
      REAL(dp), PARAMETER :: b42i0= 0.99999_dp
      REAL(dp), PARAMETER :: b42i1=-7.9238e-6_dp
!    band 3
!     liquid
      REAL(dp), PARAMETER :: b43l0= 0.99936_dp
      REAL(dp), PARAMETER :: b43l1=-0.0013632_dp
!     ice
      REAL(dp), PARAMETER :: b43i0= 0.99975_dp
      REAL(dp), PARAMETER :: b43i1=-0.001662_dp
      REAL(dp), PARAMETER :: b43i2= 6.9726e-6_dp
!    band 4
!     liquid
      REAL(dp), PARAMETER :: b44l0= 0.90032_dp
      REAL(dp), PARAMETER :: b44l1=-0.091955_dp
!     ice
      REAL(dp), PARAMETER :: b44i0= 0.89779_dp
      REAL(dp), PARAMETER :: b44i1=-0.0802_dp
!   asymmetry factor
!    band 1
!     liquid
      REAL(dp), PARAMETER :: c41l0= 0.78063_dp
      REAL(dp), PARAMETER :: c41l1= 0.12600_dp
      REAL(dp), PARAMETER :: c41l2=-0.042412_dp
!     ice
      REAL(dp), PARAMETER :: c41i0= 0.79602_dp
      REAL(dp), PARAMETER :: c41i1= 0.10183_dp
      REAL(dp), PARAMETER :: c41i2=-0.028648_dp
!    band 2
!     liquid
      REAL(dp), PARAMETER :: c42l0= 0.74102_dp
      REAL(dp), PARAMETER :: c42l1= 0.16315_dp
      REAL(dp), PARAMETER :: c42l2=-0.050268_dp
!     ice
      REAL(dp), PARAMETER :: c42i0= 0.77176_dp
      REAL(dp), PARAMETER :: c42i1= 0.11995_dp
      REAL(dp), PARAMETER :: c42i2=-0.030557_dp
!    band 3
!     liquid
      REAL(dp), PARAMETER :: c43l0= 0.70730_dp
      REAL(dp), PARAMETER :: c43l1= 0.18299_dp
      REAL(dp), PARAMETER :: c43l2=-0.045693_dp
!     ice
      REAL(dp), PARAMETER :: c43i0= 0.74691_dp
      REAL(dp), PARAMETER :: c43i1= 0.13514_dp
      REAL(dp), PARAMETER :: c43i2=-0.027140_dp
!    band 4
!     liquid
      REAL(dp), PARAMETER :: c44l0= 0.70554_dp
      REAL(dp), PARAMETER :: c44l1= 0.88798_dp
      REAL(dp), PARAMETER :: c44l2=-1.8214_dp
      REAL(dp), PARAMETER :: c44l3= 1.5775_dp
      REAL(dp), PARAMETER :: c44l4=-0.46293_dp
!     ice
      REAL(dp), PARAMETER :: c44i0= 0.77614_dp
      REAL(dp), PARAMETER :: c44i1= 0.15186_dp
      REAL(dp), PARAMETER :: c44i2=-0.031452_dp
 
END MODULE mo_sw_clop
