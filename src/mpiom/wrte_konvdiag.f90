      SUBROUTINE WRTE_KONVDIAG(KYEAR2,KMONT2,KDAYS)
!****************************************************************
!
!**** *WRTE_KONVDIAG* - save information on convective overturning.
!
!     CH,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    01.10.01
!     - separate routine extracted from OLLIE (MAIN)
!     - netCDF version possible (with cond.comp. PNETCDFO)
!
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
!
!     *CALL*       *WRTE_KONVDIAG(KYEAR2,KMONT2,KDAYS)*
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *COMMO1.h*     - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS.h*      - std I/O logical units.
!
!**   Interface to calling routine (parameter list):
!     ----------------------------------------------
!
!     *INTEGER* *KYEAR2*   - actual year.
!     *INTEGER* *KMONT2*   - actual month.
!     *INTEGER* *KDAYS*    - actual day.
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
      REAL AUX_G(IE_G,JE_G)
      INTEGER kcondep_g(IE_G,JE_G)

      CALL gather_arr(REAL(kcondep),AUX_G,p_io)

      IF(p_pe/=p_io) RETURN ! Only I/O pe does the write

      kcondep_g = aux_g

#ifdef PNETCDFO
!
! Append to NetCDF file : NOTE: the global attribute 'date' is set when 
!                               the file is opened

!
! Open NetCDF file
!
      ncstat = NF_OPEN('hopc.nc',NF_NOCLOBBER, ncid)
      IF ( ncstat .NE. NF_NOERR )                                       &
     & CALL STOP_ALL('WRTE_KONVDIAG: Problem with netCDF1')
!
! Set DEFINE mode
!
      ncstat = NF_REDEF(ncid)
      IF ( ncstat .NE. NF_NOERR )                                       &
     & CALL STOP_ALL('WRTE_KONVDIAG: Problem with define mode')
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
      ncstat = NF_DEF_VAR(ncid,'kcondep',NF_INT,2,ncdims,ncvarid)
      IF (ncstat.NE.NF_NOERR)CALL STOP_ALL('WRTE_KONVDIAG: STOP at 1')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',5, 'layer number')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,48, 'maximum layer number to which convection reached')
!
! Write variables
!
      ncstat = NF_PUT_VAR_INT(ncid,ncvarid,kcondep_g)

      ncstat = NF_CLOSE(ncid)
      IF ( ncstat .NE. NF_NOERR )                                       &
     & CALL STOP_ALL('AUFW: Problem with netCDF151')
#else
!
! Write to disk (EXTRA format)
!
      IO_OU_F090=90
!      OPEN(IO_OU_F090,STATUS='UNKNOWN',                               &
!     &    ACCESS='SEQUENTIAL',                                        &
!     &    POSITION='APPEND',                                          &
!     &    FORM='UNFORMATTED')
      I4I1=((KYEAR2*10000)+(KMONT2*100)+KDAYS)
      I4I2= 69
      I4I3=0
      I4I4=(IE_G*JE_G)
      do j=1,je_g
      do i=1,ie_g
      ff_g(i,j)=real(FLOAT(KCONDEP_g(i,j)),sp)
      enddo
      enddo


     WRITE(IO_OU_F090)I4I1,I4I2,I4I3,I4I4
      WRITE(IO_OU_F090)FF_G
!      CLOSE(IO_OU_F090)
#endif

      RETURN
      END
