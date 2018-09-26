MODULE mo_lw_clop

!- Description:
!
!  This module contains:
!
!  Coefficients for the parameterization of the emissivity
!  of cloud particles. Coefficients are given for cloud liquid droplets
!  and ice crystals.
!
!  Coefficients are given for 16 band resolution of the LW spectrum
!  as used by the RRTM scheme.
!
!  Reference : follows ECMWF-CY23R1,
!              Smith and Shi (1992), Ebert and Curry (1992)
!
!  Author:
!
!  Marco Giorgetta, MPI, December 2000

   USE mo_kind    ,ONLY: dp
   USE mo_parrrtm ,ONLY: jpband  ! number of bands in rrtm

   IMPLICIT NONE

   ! repcug = C7 in mass absorption coefficient formula
   REAL(dp), PARAMETER, DIMENSION(jpband) :: & 
        rebcug = (/ 0.718_dp, 0.726_dp, 1.136_dp, 1.320_dp, &
                    1.505_dp, 1.290_dp, 0.911_dp, 0.949_dp, &
                    1.021_dp, 1.193_dp, 1.279_dp, 0.626_dp, &
                    0.647_dp, 0.668_dp, 0.690_dp, 0.690_dp /)

   ! repcuh = C8 in mass absorption coefficient formula
   REAL(dp), PARAMETER, DIMENSION(jpband) :: & 
        rebcuh = (/ 0.0069_dp, 0.0060_dp, 0.0024_dp, 0.0004_dp, &
                   -0.0016_dp, 0.0003_dp, 0.0043_dp, 0.0038_dp, &
                    0.0030_dp, 0.0013_dp, 0.0005_dp, 0.0054_dp, &
                    0.0052_dp, 0.0050_dp, 0.0048_dp, 0.0048_dp /)
 
END MODULE mo_lw_clop
