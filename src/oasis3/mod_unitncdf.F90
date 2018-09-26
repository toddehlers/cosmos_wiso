MODULE mod_unitncdf
!
! -- unitncdf.h: 23-03-01   Version 2.5   Author: S. Valcke 
!    ********** 
!@
!@  Contents : netcdf LOGICAL flags and file ids 
!@  --------
!@
!@ -- nc_grdid: netCDF id for grids.nc file
!@ -- nc_mskid: netCDF id for masks.nc file
!@ -- nc_surid: netCDF id for areas.nc file
!@ -- nc_inpid: netCDF id for input field file
!@ -- nc_scpid: netCDF id for SCRIP matrix file
!@ -- nc_outid: netCDF id for output field file
!@ -- nc_invartimeid: netCDF id for time variable (input file)
!@
!@ -- n_reaty: Type of REAL variables
!@ -- n_intty: Type of INTEGER variables
!@ -- lncdfgrd : LOGICAL true IF all grid auxilary files are netCDF
!@ -- lncdfrst : LOGICAL true IF restart auxilary files are netCDF
!@ -- lncdfana : LOGICAL true IF analysis auxilary files are netCDF  
!@
!@ -- rtime_val : time counter as read in interpolation input file
!     ---------------------------------------------------------------   
!
  USE mod_kinds_oasis
  INTEGER (kind=ip_intwp_p) :: nc_grdid, nc_mskid, nc_surid, n_reaty, n_intty
  INTEGER (kind=ip_intwp_p) :: nc_scpid
  INTEGER (kind=ip_intwp_p),DIMENSION(:),ALLOCATABLE ::  nc_inpid, nc_outid 
  INTEGER (kind=ip_intwp_p) :: nc_invartimeid
!
  LOGICAL :: lncdfgrd, lncdfrst, lncdfana
!
  REAL(kind=ip_realwp_p) :: rtime_val
!
  CHARACTER(len=20) :: nc_invartime_name
!
!* ------------------------------------------------------------------   
!
END MODULE mod_unitncdf
!
!*====================================================================









