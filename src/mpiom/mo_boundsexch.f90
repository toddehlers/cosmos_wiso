      MODULE MO_BOUNDSEXCH

      USE mo_param1

      IMPLICIT NONE

      INTEGER,POINTER ::  p_win2dx(:),p_win2dy(:),p_win3dx(:),p_win3dy(:)

!:: TWO-DIMENSIONAL FIELDS
      REAL,POINTER ::  xr1(:,:),xr2(:,:),yr1(:,:),yr2(:,:), &
     &                 x3r1(:,:),x3r2(:,:),y3r1(:,:),y3r2(:,:) 

!
      CONTAINS

      SUBROUTINE alloc_mem_boundsexch

       ALLOCATE( p_win2dx(2),p_win2dy(2),p_win3dx(2),p_win3dy(2) )

!:: TWO-DIMENSIONAL FIELDS
       ALLOCATE( xr1(je,10),xr2(je,10),yr1(ie,10),yr2(ie,10),&
     &           x3r1(je,ke+1),x3r2(je,ke+1),                        &
     &           y3r1(ie,ke+1),y3r2(ie,ke+1) )

      END SUBROUTINE alloc_mem_boundsexch

      END MODULE MO_BOUNDSEXCH
