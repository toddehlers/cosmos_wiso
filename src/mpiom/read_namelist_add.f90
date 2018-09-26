    SUBROUTINE READ_NAMELIST_ADD

!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/read_namelist.f90,v $\\
!$Revision: 1.2.2.1.4.1.2.2.4.1.2.2.2.3.2.1 $\\
!$Date: 2006/04/03 11:27:49 $\\
!$Name: mpiom_1_2_0 $\\
!
!****************************************************************
!
!**** *READ_NAMELIST_ADD* - set namelist variables.
!
!    
!     --------
!     
!     Purpose
!     -------
!     - set default values of namelist variables
!     - read from namelist
!     - print values of namelist variables
!
!     Method
!     -------
!     -
!
!**   Interface.
!     ----------
!     called by ini_add
!
!     *CALL*       *READ_NAMELIST_ADD
!
!     *COMMON*     *PARAM1_ADD.h*   - declaration of ocean tracer.
!     *COMMON*     *COMMO1_ADD.h*   - ocean tracer arrays.
!     *COMMON*     *NAMELIST_ADD.h* - acc modules namelist.
!     *COMMON*     *UNITS_ADD.h*    - std I/O logical units.
!     *MODULE*     *mo_timeser_ADD* - parameter and memory for time series.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************
!ik added trap depths for time series

      USE mo_control_add      
      USE mo_timeser_add      
      use mo_param1_add 
      use mo_addmean 

      use mo_mpi



implicit none      

      INTEGER :: L

      NAMELIST /ADDCTL/ rmasko, kchck_add                     &
     &                 ,mean_3D_freq_add             &      
     &                 ,io_stdo_add,io_stdi_add                       &
     &                 ,io_rsti_add,io_rsto_add,io_timeser_add        &
     &                 ,rlonts1,rlatts1,nfreqts1                      &
     &                 ,rdep1ts1,rdep2ts1,rdep3ts1                    


! Coordinates etc. of time series 1 
      nfreqts1=0
      DO l=1,nts
         rlonts1(l)=-1.0
         rlatts1(l)=-1.0
         rdep1ts1(l) = 100.
         rdep2ts1(l) = 200.
         rdep3ts1(l) = 300.
      ENDDO

!                        
! Read ADDCTL namelist
!
      IF(p_pe==p_io) THEN
        OPEN(io_stdi_add,FILE='NAMELIST_ADD',STATUS='UNKNOWN',           &
     &        ACCESS='SEQUENTIAL',FORM='FORMATTED')

        READ(io_stdi_add,ADDCTL)

        CLOSE(io_stdi_add)
      ENDIF

      CALL p_bcast(rmasko,p_io)
      CALL p_bcast(kchck_add,p_io)
      CALL p_bcast(mean_3D_freq_add,p_io)
      CALL p_bcast(io_stdo_add,p_io)
      CALL p_bcast(io_stdi_add,p_io)
      CALL p_bcast(io_rsti_add,p_io)
      CALL p_bcast(io_rsto_add,p_io)
      CALL p_bcast(io_timeser_add,p_io)
      CALL p_bcast(rlonts1,p_io)
      CALL p_bcast(rlatts1,p_io)
      CALL p_bcast(nfreqts1,p_io)
      CALL p_bcast(rdep1ts1,p_io)
      CALL p_bcast(rdep2ts1,p_io)
      CALL p_bcast(rdep3ts1,p_io)




!
!  Open std out file
!

      CALL OPEN_STDOUT(io_stdo_add,'addout')

      WRITE(io_stdo_add,*)                                             &
     &'****************************************************************'
      WRITE(io_stdo_add,*)                                             &
     &'* '
      WRITE(io_stdo_add,*)                                             &
     &'* Values of ADDCTL namelist variables : '
      WRITE(io_stdo_add,*)                                             &
     &'* Time series :'
      WRITE(io_stdo_add,*)                                             &
     &'*                                nts      = ',nts,'  elements.'
      WRITE(io_stdo_add,*)                                             &
     &'*     sampled at / averaged over nfreqts1 = ',nfreqts1          &
     &                                             ,' model time steps.'
      DO l=1,nts
        WRITE(io_stdo_add,*)                                           &
     &'*         Position of element ',l,' at    = '                   &
     &                                          ,rlonts1(l),rlatts1(l)
      ENDDO

      WRITE(io_stdo_add,*)                                             &
     &'* Standard output unit        io_stdo_add = ',io_stdo_add
      WRITE(io_stdo_add,*)                                             &
     &'* '
     
      WRITE(io_stdo_add,*) '* mean_3D_freq is set to: ', mean_3D_freq_add





      WRITE(io_stdo_add,*)                                             &
     &'* '
      WRITE(io_stdo_add,*)                                             &
     &'****************************************************************'

      RETURN
      END
