      MODULE mo_timeser_add

!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/mo_timeser_bgc.f90,v $\\
!$Revision: 1.3.10.1.2.2.4.1.2.2.2.3.2.1 $\\
!$Date: 2006/04/03 11:27:49 $\\
!$Name: mpiom_1_2_0 $\\
!
!***********************************************************************
!
!**** *MODULE mo_timeser_add* - Parameter and memory for time series.
!
!    
!     Modified
!     --------
!     
!     Purpose
!     -------
!     - declaration
!
!     *nfreqts1*     *INTEGER*  - sampling frequency of time series 1.
!     *nts1*         *INTEGER*  - no. of elements in time series 1.
!     *nvarts1_add*      *INTEGER*  - no. of variables sampled for each element.
!     *its1*         *INTEGER*  - 1st index of grid cell for samples in time series 1.
!     *jts1*         *INTEGER*  - 2nd index of grid cell for samples in time series 1.
!     *lts1*         *INTEGER*  - actual sample counter of time series 1.
!
!**********************************************************************
      implicit none

      INTEGER, PARAMETER :: nts=8

      REAL, DIMENSION (:,:,:), ALLOCATABLE :: ts1

      INTEGER :: lents1,nelets1

      INTEGER :: io_timeser_add = 17  ! logical unit number for time series output

      INTEGER :: its1(nts),jts1(nts),nts1,nfreqts1,lts1
      INTEGER :: k1ts1(nts),k2ts1(nts),k3ts1(nts)

      REAL :: rlonts1(nts),rlatts1(nts)
      REAL :: rdep1ts1(nts),rdep2ts1(nts),rdep3ts1(nts)



      INTEGER, PARAMETER ::                                    &  
     &          i_add_ts   =6,                                 &                               
     &          itsh2o18  =1,                                  &
     &          itshDo16  =3,                                  &
     &          itsh2o16  =2,                                  & 
     &          itsh2o18_ice =4,                               &
     &          itsh2o16_ice =5,                               &
     &          itshDo16_ice =6                                
    
      INTEGER, PARAMETER :: nvarts1_add =i_add_ts


      END MODULE mo_timeser_add
