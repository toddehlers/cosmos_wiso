      MODULE mo_contra
!***********************************************************************
!
!**** *MODULE mo_contra* - Variables for additional conservative tracers.
!
!    
!     - declaration and memory allocation
!
!     *ocectra*       *REAL*  - .
!     *h2o18*         *REAL*  - .
!     *h2o16*         *REAL*   
!     
!     - forcing field
!     *pre16*         *REAL*
!     *pre18*
!     *riv16*
!     *riv18*
!
!**********************************************************************
      
      implicit none
      
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: ocectra
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: icectra
     
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: h2o18
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: hDo16
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: h2o16
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: h2o18_ice
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: hDo16_ice
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: h2o16_ice    
 
      REAL, DIMENSION (:,:), ALLOCATABLE :: pre16
      REAL, DIMENSION (:,:), ALLOCATABLE :: pre18
      REAL, DIMENSION (:,:), ALLOCATABLE :: preD
      REAL, DIMENSION (:,:), ALLOCATABLE :: riv16 
      REAL, DIMENSION (:,:), ALLOCATABLE :: riv18 
      REAL, DIMENSION (:,:), ALLOCATABLE :: rivD

      REAL :: totarea_add
      

      CONTAINS

      SUBROUTINE ALLOC_MEM_CONTRA(kpie,kpje,kpke)

      use mo_control_add
      use mo_param1_add 

      INTEGER :: kpie,kpje,kpke
      
        WRITE(io_stdo_add,*)'Memory allocation for variable oceanic passive tracer'
        WRITE(io_stdo_add,*)'First dimension    : ',kpie
        WRITE(io_stdo_add,*)'Second dimension   : ',kpje
        WRITE(io_stdo_add,*)'Third dimension    : ',kpke
        WRITE(io_stdo_add,*)'Forth dimension    : ',nocectra

        ALLOCATE (ocectra(kpie,kpje,kpke,nocectra))

        WRITE(io_stdo_add,*)'Memory allocation for variable ice passive tracer'
        WRITE(io_stdo_add,*)'First dimension    : ',kpie
        WRITE(io_stdo_add,*)'Second dimension   : ',kpje
        WRITE(io_stdo_add,*)'Third dimension    : ',kpke
        WRITE(io_stdo_add,*)'Forth dimension    : ',nicectra

        ALLOCATE (icectra(kpie,kpje,kpke,nicectra)) 
        
      END SUBROUTINE ALLOC_MEM_CONTRA


      SUBROUTINE ALLOC_MEM_FORC(kpie,kpje)

      use mo_control_add
      use mo_param1_add 

      INTEGER :: kpie,kpje
      
        WRITE(io_stdo_add,*)'Memory allocation for forcing field of variable oceanic passive tracer'
        WRITE(io_stdo_add,*)'First dimension    : ',kpie
        WRITE(io_stdo_add,*)'Second dimension   : ',kpje
       
        ALLOCATE (pre16(kpie,kpje))
        ALLOCATE (pre18(kpie,kpje))
        ALLOCATE (preD(kpie,kpje))
        ALLOCATE (riv16(kpie,kpje))
        ALLOCATE (riv18(kpie,kpje))
        ALLOCATE (rivD(kpie,kpje))
        
      END SUBROUTINE ALLOC_MEM_FORC

      END MODULE mo_contra
