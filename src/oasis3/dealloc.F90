SUBROUTINE dealloc
!
!**** DEALLOC
!
!     Purpose:
!       Deallocate arrays defined in the modules
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "inialloc".       
!
!     External:
!       none
!
!     Files:
!       none
!   
!     References:
!
!     History:
!     --------
!       Version   Programmer     Date        Description
!       ------------------------------------------------
!       2.5       A.Caubel       2002/03/18  created
!
!*-----------------------------------------------------------------------
!
!** + DECLARATIONS
!
!** ++ Use of modules
!
  USE mod_parameter
!
!*-----------------------------------------------------------------------
  IF (lg_oasis_field) THEN 
!
!** + Module anais
!
   CALL dealloc_anais
!
!** + Module analysis
!
   CALL dealloc_analysis
!
!** + Module coast
!
   CALL dealloc_coast
!
!** + Module extrapol
!
   CALL dealloc_extrapol 
!
!** + Module memory
!
   CALL dealloc_memory
!
!** + Module nproc
!
   CALL dealloc_nproc
!
!** + Module parallel
!
   CALL dealloc_parallel
!
!** + Module pipe
!
#ifdef use_comm_PIPE
   CALL dealloc_pipe
#endif
!
!** + Module rainbow
!
   CALL dealloc_rainbow
!
!** + Module sipc
!
#ifdef use_comm_SIPC
   CALL dealloc_sipc
#endif
!
!** + Module mod_gsip
!
#ifdef use_comm_GSIP
     CALL dealloc_gsip
#endif
!
!** + Module timestep
!
   CALL dealloc_timestep
!
!** + Module unitncdf
!
   CALL dealloc_unitncdf

ENDIF
!
!** + Module experiment
!
CALL dealloc_experiment
!
!** + Module string
!
CALL dealloc_string
!
!*------------------------------------------------------------------------
!
END SUBROUTINE dealloc
!
!*========================================================================
