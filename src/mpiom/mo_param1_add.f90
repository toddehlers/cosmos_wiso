      MODULE mo_param1_add


!***********************************************************************
!
!**** *MODULE mo_param1_add* - add(additional conservative) tracer parameters.
!
!     
!  
!     Purpose
!     -------
!     - declaration and memory allocation
!
!**********************************************************************
      implicit none
      
     
      INTEGER, PARAMETER :: kwrbioz_add=8                            ! euphotic layers


! advected tracers
      INTEGER, PARAMETER :: i_base_adv_add=3,                              &
     &                      ih2o18  =1,                                    &
     &                      ihDo16  =3,                                    &
     &                      ih2o16  =2                             



! total number of advected tracers
      INTEGER, PARAMETER :: nctraad=i_base_adv_add

      INTEGER, PARAMETER :: nocectra  = nctraad

                           
 ! non-advected tracers in ice
       INTEGER, PARAMETER ::                                               &
      &                     ih2o18_ice = 1,                                &
      &                     ih2o16_ice = 2,                                &
      &                     ihDo16_ice = 3,                                &
      &                     i_base_ice = 3
     
      INTEGER, PARAMETER :: nicectra  = i_base_ice


      END MODULE mo_param1_add
