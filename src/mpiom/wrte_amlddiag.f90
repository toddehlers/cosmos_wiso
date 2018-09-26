      SUBROUTINE WRTE_AMLDDIAG(KYEAR,KMONTH,KDAY)
!C****************************************************************
!C
!C**** *WRTE_AMLDDIAG* - save information on mixed layer depth.
!C
!     CH,    *MPI-Met, HH*    10.04.01
!     S.Legutke,        *MPI-MaD, HH*    01.10.01
!     - separate routine extracted from OLLIE (MAIN)
!     - netCDF version possible (with cond.comp. PNETCDFO)
!     - kind sp HH 12.10.04 
!     Purpose
!     -------
!     
!
!     Method
!     -------
!     
!
!**   Interface.
!     ----------
  
!     *CALL*       *WRTE_AMLDDIAG(KYEAR,KMONTH,KDAY)*
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *COMMO1.h*     - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS.h*      - std I/O logical units.
!
!**   Interface to calling routine (parameter list):
!     ----------------------------------------------
!
!     *INTEGER* *KYEAR*   - actual year.
!     *INTEGER* *KMONT*   - actual month.
!     *INTEGER* *KDAY*    - actual day.
!
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************

      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_UNITS
      USE MO_KIND

      INTEGER(KIND=i4) I4I1,I4I2,I4I3,I4I4
      REAL(KIND=sp) FF_G(IE_G,JE_G)

#ifdef PNETCDFO
      INCLUDE 'netcdf.inc'

      INTEGER ncid,ncstat,ncvarid,                                      &
     &        nclonid,nclatid,nclevid,nclontoid,nclattoid
#endif

      REAL AMLD_G(IE_G,JE_G)

      CALL gather_arr(AMLD(:,:,LMONTS),AMLD_G,p_io)

      IF(p_pe/=p_io) RETURN ! Only I/O pe does the write

#ifdef PNETCDFO
!
! Append to NetCDF file : NOTE: the global attribute 'date' is set when 
!                               the file is opened

!
! Open NetCDF file
!
      ncstat = NF_OPEN('hopc.nc',NF_NOCLOBBER, ncid)
      IF ( ncstat .NE. NF_NOERR )                                       &
     & CALL STOP_ALL('WRTE_AMLDDIAG: Problem with netCDF1')
!
! Set DEFINE mode
!
      ncstat = NF_REDEF(ncid)
      IF ( ncstat .NE. NF_NOERR )                                       &
     & CALL STOP_ALL('WRTE_AMLDDIAG: Problem with define mode')
!
! Get dimension IDs
!
      ncstat = NF_INQ_DIMID(ncid,'lon',nclonid)
      ncstat = NF_INQ_DIMID(ncid,'lat',nclatid)
      ncdims(1) = nclonid
      ncdims(2) = nclatid

!
! Define variables
!
      ncstat = NF_DEF_VAR(ncid,'amld',NF_DOUBLE,2,ncdims,ncvarid)
      IF (ncstat.NE.NF_NOERR) CALL STOP_ALL('WRTE_AMLDDIAG: STOP at 1')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',5, 'meter')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,17, 'mixed layer depth')
!
! Write variables
!
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,amld_g)

      ncstat = NF_CLOSE(ncid)
      IF ( ncstat .NE. NF_NOERR )                                       &
     &   CALL STOP_ALL('AUFW: Problem with netCDF151')

#else
!
! Write to disk (EXTRA format)
!
!          OPEN(IO_OU_AMLD,      STATUS='UNKNOWN',                       &
!     &                  ACCESS='SEQUENTIAL',                            &
!     &                POSITION='APPEND',                                &
!     &                    FORM='UNFORMATTED')
          I4I1=(KYEAR*10000)+(KMONTH*100)+KDAY
          I4I2=83
          I4I3=0
          I4I4=(IE_G*JE_G)

          do J=1,JE_G
          do i=1,ie_G
	   FF_G(i,j)=real(AMLD_G(I,J),sp)
          enddo
          enddo



          WRITE(IO_OU_AMLD)I4I1,I4I2,I4I3,I4I4
          WRITE(IO_OU_AMLD)FF_G
!          CLOSE(IO_OU_AMLD)
#endif
       RETURN
       END
