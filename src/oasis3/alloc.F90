SUBROUTINE alloc (id_flag)
!
!
!**** INIALLOC
!
!     Purpose:
!       Allocate arrays defined in the modules
!
!     Interface:
!       none
!    
!     Method:
!       Uses run parameters read in "inipar_alloc" routine to 
!       allocate arrays.       
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
!* ---------------------------- Local declarations ----------------------
!
     INTEGER (kind=ip_intwp_p) :: id_flag
!
!*-----------------------------------------------------------------------
!
IF (id_flag .EQ. 1) THEN
!** + Module mod_experiment
!
!  CALL alloc_experiment
!
!** + Module mod_string
!
  CALL alloc_string
!
!** + Following routines will be called only if one field (at least) goes
!     through Oasis  
!
  IF (lg_oasis_field) THEN 
!
!** + Module mod_anais
!
     CALL alloc_anais1
!
!** + Module mod_analysis
!
     CALL alloc_analysis
!
!** + Module mod_coast
!
     CALL alloc_coast
!
!** + Module mod_extrapol
!
     CALL alloc_extrapol1 
!
!** + Module mod_memory
!
     CALL alloc_memory1
!
!** + Module mod_nproc
!
     CALL alloc_nproc
!
!** + Module mod_parallel
!
     CALL alloc_parallel
!
!** + Module mod_pipe
!
#ifdef use_comm_PIPE
     CALL alloc_pipe
#endif
!
!** + Module mod_rainbow
!
     CALL alloc_rainbow1
!
!** + Module mod_sipc
!
#ifdef use_comm_SIPC
     CALL alloc_sipc
#endif
!
!** + Module mod_gsip
!
#ifdef use_comm_GSIP
     CALL alloc_gsip
#endif
!
!!** + Module mod_timestep
!
     CALL alloc_timestep
!
!** + Module mod_unitncdf
!
     CALL alloc_unitncdf
!
  ENDIF
ELSE
  IF (lg_oasis_field) THEN 
!
!** + Module mod_anais
!
     CALL alloc_anais2
!
!** + Module mod_extrapol
!
     CALL alloc_extrapol2 
!
!** + Module mod_memory
!
     CALL alloc_memory2
!
!** + Module mod_rainbow
!
     CALL alloc_rainbow2
!
   ENDIF
ENDIF
!*------------------------------------------------------------------------
!
END SUBROUTINE alloc
!
!*========================================================================

