      SUBROUTINE CLOSE_ADDMEAN_3D
      
!****************************************************************
!
!**** *WRITE_ADDMEAN_3D* - calculate and write 3-dimensional add mean data.
!

!     --------
!
!     Purpose
!     -------
!     Write add mean data.
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *CLOSE_ADDMEAN_3D*
!
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************

      USE mo_addmean
      USE mo_contra

      use mo_param1_add 

      USE mo_control_add
      USE mo_parallel
 
      implicit none


      INCLUDE 'netcdf.inc'
      INTEGER :: ncstat


      
!-----------------------------------------------------------------------

      if (p_pe==p_io) then

 
      WRITE(io_stdo_add,*) 'Close 3D addmean data at step:', meantime_3d_add


!
! Close File
!
!-----------------------------------------------------------------------
!
      ncstat = NF_CLOSE(nc_3d_id_add)
      IF ( ncstat .NE. NF_NOERR ) call STOP_all('CLOSE_ADDMEAN_3D: Problem with netCDF200')

      end if


      RETURN
      END
