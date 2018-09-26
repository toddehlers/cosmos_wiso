      MODULE mo_addmean


!***********************************************************************
!
!**** *MODULE mo_addmean* - Variables for addmean.
!
!     Patrick Wetzel    *MPI-Met, HH*    09.12.02
!  
!     Purpose
!     -------
!     - declaration and memory allocation for output fields
!
!**********************************************************************
      implicit none

!----------------------------------------------------------------
! ice non adv fields

      INTEGER, PARAMETER ::                                   &
     &          jh2o18_ice =1,                                  &
     &          jhDo16_ice =3,                                  &
     &          jh2o16_ice =2,                                  & 
     &          i_bsc_ice_add =3

      INTEGER, PARAMETER :: naddice =i_bsc_ice_add

!----------------------------------------------------------------
!  3d 'total' fields
     
      INTEGER, PARAMETER ::                                   &
     &          jh2o18_t =1,                                  &
     &          jhDo16_t =3,                                  &
     &          jh2o16_t =2,                                  & 
     &          i_bsc_t3d_add=3

      
      INTEGER, PARAMETER :: naddt3d =i_bsc_t3d_add
     
      
      REAL, DIMENSION (:),       ALLOCATABLE :: stepspm    
       
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: addice               ! ice summed-up field
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: addt3d               ! 3d summed-up field
      
      
      INTEGER :: meancnt_add_2D, meancnt_add_3D, meancnt_add_ice

      INTEGER :: mean_2D_freq_add, mean_3D_freq_add, mean_ice_freq_add
      INTEGER :: meantime_2d_add, nmeantime_2d_add
      INTEGER :: meantime_3d_add, nmeantime_3d_add
      INTEGER :: meantime_ice_add, nmeantime_ice_add

      INTEGER :: nc_2d_id_add, nc_3d_id_add
      
      
      CONTAINS

      SUBROUTINE ALLOC_MEM_ADDMEAN(kpie,kpje,kpke)

      use mo_control_add
      use mo_param1_add 	
      
      INTEGER :: kpie,kpje,kpke
  
        WRITE(io_stdo_add,*)'Memory allocation for variable stepspm ...'
	    WRITE(io_stdo_add,*)'First dimension    : ',nmeantime_2d_add
	
        ALLOCATE (stepspm(nmeantime_2d_add))  

        WRITE(io_stdo_add,*)'Memory allocation for variable addt2d ...'
        WRITE(io_stdo_add,*)'First dimension    : ',kpie
        WRITE(io_stdo_add,*)'Second dimension   : ',kpje
        WRITE(io_stdo_add,*)'Third dimension    : ',kpke
        WRITE(io_stdo_add,*)'Forth dimension    : ',naddice
	
        ALLOCATE (addice(kpie,kpje,kpke,naddice))	  
	
        WRITE(io_stdo_add,*)'Memory allocation for variable addt3d ...'
        WRITE(io_stdo_add,*)'First dimension    : ',kpie
        WRITE(io_stdo_add,*)'Second dimension   : ',kpje
	    WRITE(io_stdo_add,*)'Third dimension    : ',kpke 
	    WRITE(io_stdo_add,*)'Forth dimension    : ',naddt3d
	
        ALLOCATE (addt3d(kpie,kpje,kpke,naddt3d))		

       

      END SUBROUTINE ALLOC_MEM_ADDMEAN

      END MODULE mo_addmean
