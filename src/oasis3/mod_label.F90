MODULE mod_label
!
! -- label.h   18-08-95   Version 2.0beta   Author: Laurent Terray
!    *******   01-02-96   Version 2.0 : disparition of ctranam
!              07-08-96   Version 2.1 : rise of cfldlab dimension
!              25-09-96   Version 2.1 : rise of cfldlab length
!@
!@  Contents : variables related to labelling files and character variable
!@  --------
!@
!@ -- cfldlab : field label definition (1D)
!@
!@ -- cgrdnam : name for grid file
!@
!@ -- cmsknam : name for mask file
!@
!@ -- csurnam : name for surface file
!@
!@ -- cglonsuf : suffix name for longitude grid file locator
!@
!@ -- cglatsuf : suffix name for longitude latitude grid file locator
!@
!@ -- cmsksuf : suffix name for mask file locator
!@
!@ -- crnlonsuf : suffix name for corner longitude grid file locator
!@
!@ -- crnlatsuf : suffix name for corner latitude grid file locator
!@
!@ -- csursuf : suffix name for surface file locator
!@
!@ -- crednam : name for reduced gaussian grid masks file
!@
! -------------------------------------------------------------------
!
  CHARACTER(len=64),ALLOCATABLE::cfldlab(:)
  CHARACTER(len=5), PARAMETER :: cgrdnam = 'grids'  
  CHARACTER(len=5), PARAMETER :: cmsknam = 'masks' 
  CHARACTER(len=5), PARAMETER :: csurnam = 'areas' 
  CHARACTER(len=5), PARAMETER :: crednam = 'maskr'
  CHARACTER(len=4), PARAMETER :: cglonsuf = '.lon'
  CHARACTER(len=4), PARAMETER :: cglatsuf = '.lat'
  CHARACTER(len=4), PARAMETER :: crnlonsuf = '.clo'
  CHARACTER(len=4), PARAMETER :: crnlatsuf = '.cla'
  CHARACTER(len=4), PARAMETER :: cmsksuf = '.msk'
  CHARACTER(len=4), PARAMETER :: csursuf = '.srf'
!
END MODULE mod_label
! -------------------------------------------------------------------

