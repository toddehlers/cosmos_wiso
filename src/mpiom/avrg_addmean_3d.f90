      SUBROUTINE AVRG_ADDMEAN_3D(kpie,kpje,kpke)
!*******************************************************************
!
!**** *AVRG_BGCMEAN* - average addmean data.
!
!   
!     
!     Purpose
!     -------
!
!     Method
!     -------
!     -
!
!**   Interface.
!     ----------
!
!     *CALL*       *AVRG_ADDMEAN_3D(kpie,kpje,kpke)*
!
!     *MODULES*     *mo_contra* - ocean tracer arrays.
!     *MODULES*     *mo_control_add*  - std I/O logical units.
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


      USE mo_contra
      USE mo_control_add
      USE mo_addmean
      use mo_param1_add 

implicit none    

      INTEGER :: kpie,kpje,kpke,i,j,k,l




      WRITE(io_stdo_add,*) 'AVRG_ADDMEAN_3D using:   mean_3D_freq_add='  &
                                                      ,mean_3D_freq_add
      WRITE(io_stdo_add,*) 'AVRG_ADDMEAN_3D using:    meantime_3d_add='  &
                                                      ,meantime_3d_add   
      WRITE(io_stdo_add,*) 'AVRG_ADDMEAN_3D using: meancnt_add_3D='&
                                                     ,meancnt_add_3D




      WRITE(io_stdo_add,*) 'AVRG_ADDMEAN_3D at step: ',meantime_3d_add
      WRITE(io_stdo_add,*) 'avrg. total 3D fields with ', meancnt_add_3D
        
	DO l=1,naddt3d
        DO k=1,kpke
        DO j=1,kpje
        DO i=1,kpie	
           addt3d(i,j,k,l)=addt3d(i,j,k,l)/meancnt_add_3D
        ENDDO
        ENDDO
        ENDDO
        ENDDO

	DO l=1,naddice
        DO k=1,kpke
        DO j=1,kpje
        DO i=1,kpie	
           addice(i,j,k,l)=addice(i,j,k,l)/meancnt_add_3D
        ENDDO
        ENDDO
        ENDDO
        ENDDO

! Reset counter  
    
      meancnt_add_3D=0


      RETURN
      END
