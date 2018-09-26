MODULE mo_lookup_tables

  !
  ! This module initializes and holds the lookup tables of JSBACH
  ! (adapted from the ECHAM5 routine "setphys")
  !
  USE mo_kind,         ONLY: dp
  
  IMPLICIT NONE
  private

  INTEGER,PARAMETER :: jptlucu1 =  50000  ! lookup table lower bound (in milli-Kelvin)
  INTEGER,PARAMETER :: jptlucu2 = 400000  ! lookup table upper bound (in milli-Kelvin)

  ! === BEGIN OF PUBLIC PART =====================================================================================================

  public :: init_lookup_table !! subroutine that initializes the lookup tables --- to be called when initializing JSBACH

  REAL(dp),public, save :: sat_pressure_scaled(jptlucu1:jptlucu2) ! lookup table for e_s*Rd/Rv, where e_s is the saturation vapour 
                                                              ! pressure in Pascal and Rd and Rv are the gas constants of dry 
                                                              ! air (Rd=287.05 J/(K*kg)) and vapour (Rv=461.51 J/(K*kg)). The 
                                                              ! index denotes the temperature in  milli-Kelvin, e.g. 
                                                              ! sat_pressure_scaled(300000) is e_s*Rd/Rv at 300 Kelvin.


  ! === END OF PUBLIC PART =======================================================================================================

  REAL(dp), PARAMETER :: tmelt = 273.15_dp      ! melting temperature of ice/snow
  REAL(dp), PARAMETER :: rd    = 287.05_dp      ! gas constant for dry air in J/(K*kg)
  REAL(dp), PARAMETER :: rv    = 461.51_dp      ! gas constant for water vapour in J/(K*kg)

  contains

  ! --- init_lookup_table ---------------------------------------------------------------------------------------------------------
  ! This routine initializes the lookup table(s)
  subroutine init_lookup_table
    REAL(dp) :: zavl1,zavl2,zavl3,zavl4,zavl5
    REAL(dp) :: zavi1,zavi2,zavi3,zavi4,zavi5
    REAL(dp) :: zavm1,zavm2,zavm3,zavm4,zavm5
    integer :: i,it
    REAL(dp) :: tt,zt

    ! --- initialize lookup table for saturation pressure ---

    zavl1=-6096.9385_dp
    zavl2=21.2409642_dp
    zavl3=-2.711193_dp
    zavl4=1.673952_dp
    zavl5=2.433502_dp

    zavi1=-6024.5282_dp
    zavi2=29.32707_dp
    zavi3=1.0613868_dp
    zavi4=-1.3198825_dp
    zavi5=-0.49382577_dp      

    tt=50._dp
    DO i=jptlucu1,jptlucu2
       zt=tt*1000._dp
       it=NINT(zt)
       IF(tt-tmelt.GT.0._dp) THEN
          zavm1=zavl1
          zavm2=zavl2
          zavm3=zavl3
          zavm4=zavl4
          zavm5=zavl5
       ELSE
          zavm1=zavi1
          zavm2=zavi2
          zavm3=zavi3
          zavm4=zavi4
          zavm5=zavi5
       END IF
       sat_pressure_scaled(it)=EXP((zavm1/tt+zavm2+zavm3*tt/100._dp+zavm4*tt*tt/1.e5_dp+zavm5*LOG(tt)))*rd/rv
       tt=tt+0.001_dp
    end DO
  end subroutine init_lookup_table
end MODULE mo_lookup_tables
