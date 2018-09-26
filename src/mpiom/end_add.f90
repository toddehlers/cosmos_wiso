      SUBROUTINE END_ADD(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,      &
     &        pgila,pgiph,ptiestu,kplyear,kplmon,kplday,kpldtoce)
!****************************************************************
!
!**** *END_ADD* - finish with marine addtional conservetive module.
!
!    
!     
!     Purpose
!     -------
!     - call inventory
!     - save time series
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!     called by mpiom.f90
!
!     *CALL*       *END_ADD(list....)*
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st REAL :: of model grid.
!     *INTEGER* *kpje*    - 2nd REAL :: of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) REAL :: of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd REAL ::) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st REAL ::) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd REAL ::) [m].
!     *REAL*    *pgila*   - geographical longitude of grid points [degree E].
!     *REAL*    *pgiph*   - geographical latitude  of grid points [degree N].
!     *REAL*    *ptiestu* - depth of layers [m].
!
!
!     Externals
!     ---------
!     INVENTORY_ADD SAVE_TIMESER_ADD
!
!**********************************************************************
      USE mo_contra

      USE mo_addmean
      USE mo_control_add
      use mo_param1_add

      use mo_parallel
implicit none

      INTEGER :: kpie,kpje,kpke,i,j,k
      INTEGER :: kplyear,kplmon,kplday,kpldtoce
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL :: pgila(kpie*2,kpje*2)
      REAL :: pgiph(kpie*2,kpje*2)
      REAL :: ptiestu(kpke)




!       write(0,*) 'calling close_addmean_3d from end_add'
      CALL CLOSE_ADDMEAN_3D



!
! Global inventory of all tracers
!
      CALL INVENTORY_ADD(kpie,kpje,kpke)

!
! save the time series of add
!
      CALL SAVE_TIMESER_ADD

     

      RETURN
      END
