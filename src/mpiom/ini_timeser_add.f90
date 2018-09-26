      SUBROUTINE INI_TIMESER_ADD(kpke,ptiestw)

!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/ini_timeser_bgc.f90,v $\\
!$Revision: 1.2.10.1.2.2.4.1.2.2.2.3.2.1 $\\
!$Date: 2006/04/03 11:27:49 $\\
!$Name: mpiom_1_2_0 $\\
!
!****************************************************************
!
!**** *INI_TIMER_ADD* - initialize add time series.
!
!     S.Legutke,        *MPI-MaD, HH*    13.08.01
!
!     Modified
!     --------
!     
!     Purpose
!     -------
!     Calculation of :
!     - grid indices of element's lat/lon positions.
!     - number of samples in run.
!     - allocate memory.
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *INI_TIMSER_ADD
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *PARAM1_ADD.h* - declaration of ocean/sediment tracer.
!     *COMMON*     *COMMO1_ADD.h* - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS_ADD.h*  - std I/O logical units.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* **   - .
!     *REAL*    **   - .
!
!     Externals
!     ---------
!     ALLOCATE
!
!**********************************************************************
!ik the time series include traps at three different depths now; these
!ik will be checked vs. the model depths

      USE mo_timeser_add
      USE mo_control_add
implicit none
      INTEGER :: kpke,k,l,i1,i2,i3
      REAL :: dist
      REAL :: ptiestw(kpke+1)
!
! Initialize time series sample counter.
!
      lts1 = 1
!                        
! Search grid cell index of time series 1.
!
      DO l=1,nts
         CALL SUCHIJ(rlatts1(l),rlonts1(l),1,its1(l),jts1(l),dist,1.)
         WRITE(io_stdo_add,*)                                 &
     &       'Element ',l,' of time series 1 is sampled at '  &
     &       ,rlatts1(l),' N ',rlonts1(l),                    &
     &       ' E. Grid index is ',its1(l),' and ',jts1(l)
      ENDDO

!ik assign depth indices for time series 1
      DO l=1,nts
        DO k = 1, kpke
           if (ptiestw(k).lt.rdep1ts1(l)) k1ts1(l) = k      
           if (ptiestw(k).lt.rdep2ts1(l)) k2ts1(l) = k   
           if (ptiestw(k).lt.rdep3ts1(l)) k3ts1(l) = k   
        ENDDO
      ENDDO        
              
!                        
! Allocate memory for ...
!

! ... time series 1
      IF(nfreqts1 .NE. 0) THEN
	 IF (MOD(ndtrunadd,nfreqts1).EQ.0) THEN
           lents1  = ndtrunadd / nfreqts1
	 ELSE
	   lents1  = ndtrunadd / nfreqts1 + 1
	 ENDIF
         nelets1 = nts + 1 ! No. of stations + global average
	 
         WRITE(io_stdo_add,*)'Memory allocated for timeseries 1 ...'
         WRITE(io_stdo_add,*)'No. of variables   : ',nvarts1_add
         WRITE(io_stdo_add,*)'No. of stations    : ',nelets1
         WRITE(io_stdo_add,*)'No. of step in run : ',ndtrunadd
         WRITE(io_stdo_add,*)'No. of timesteps   : ',lents1
         ALLOCATE (ts1(nvarts1_add,nelets1,lents1))
         DO i3=1, lents1
         DO i2=1,nelets1
         DO i1=1,nvarts1_add
            ts1(i1,i2,i3) = 0.0
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      

      RETURN
      END
