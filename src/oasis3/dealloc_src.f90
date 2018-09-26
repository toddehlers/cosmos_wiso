!                       ****************************
!                       * DEALLOCATION SUBROUTINES *
!                       ****************************
SUBROUTINE dealloc_anais
!
!**** DEALLOC_ANAIS
!
!     Purpose:
!       Deallocate arrays defined in the "anais" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_anais" or "inipar_alloc"
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
  USE mod_kinds_oasis
  USE mod_anais
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  DEALLOCATE (varmul, stat=il_err)
  DEALLOCATE (niwtm)
  DEALLOCATE (niwtg)
  DEALLOCATE (linit)
  DEALLOCATE (ngint)
  DEALLOCATE (nmint)
  DEALLOCATE (nmesh)
  DEALLOCATE (agint)
  DEALLOCATE (amint)
  IF (allocated(naismfl)) DEALLOCATE (naismfl)
  IF (allocated(naismvoi)) DEALLOCATE (naismvoi)
  IF (allocated(naisgfl)) DEALLOCATE (naisgfl)
  IF (allocated(naisgvoi)) DEALLOCATE (naisgvoi)
  
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_anais
!
!*========================================================================
SUBROUTINE dealloc_analysis
!
!**** DEALLOC_ANALYSIS
!
!     Purpose:
!       Deallocate arrays defined in the "analysis" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_analysis".     
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
  USE mod_analysis
!
!*-----------------------------------------------------------------------
!
  IF (allocated(ntronca)) DEALLOCATE (ntronca)
  DEALLOCATE (ncofld)
  DEALLOCATE (neighborg)
  DEALLOCATE (nludat)
  DEALLOCATE (nlufil)
  DEALLOCATE (nlumap)
  DEALLOCATE (nlusub)
  DEALLOCATE (nluext)
  DEALLOCATE (nosper)
  DEALLOCATE (notper)
  DEALLOCATE (ntinpflx)
  DEALLOCATE (ntoutflx)
  DEALLOCATE (amskval)
  DEALLOCATE (amskvalnew)
  DEALLOCATE (acocoef)
  DEALLOCATE (abocoef)
  DEALLOCATE (abncoef)
  DEALLOCATE (afldcoef)
  DEALLOCATE (afldcobo)
  DEALLOCATE (afldcobn)
  DEALLOCATE (cxordbf)
  DEALLOCATE (cyordbf)
  DEALLOCATE (cxordaf)
  DEALLOCATE (cyordaf)
  DEALLOCATE (cgrdtyp)
  DEALLOCATE (cfldtyp)
  DEALLOCATE (cfilfic)
  DEALLOCATE (cfilmet)
  DEALLOCATE (cconmet)
  DEALLOCATE (cfldcoa)
  DEALLOCATE (cfldfin)
  DEALLOCATE (ccofld)
  DEALLOCATE (cbofld)
  DEALLOCATE (cbnfld)
  DEALLOCATE (ccofic)
  DEALLOCATE (cdqdt)
  DEALLOCATE (cgrdmap)
  DEALLOCATE (cmskrd)
  DEALLOCATE (cgrdsub)
  DEALLOCATE (ctypsub)
  DEALLOCATE (cgrdext)
  DEALLOCATE (csper)
  DEALLOCATE (ctper)
  DEALLOCATE (lsurf)
  DEALLOCATE (nscripvoi)
  DEALLOCATE (cmap_method)
  DEALLOCATE (cfldtype)
  DEALLOCATE (crsttype)
  DEALLOCATE (nbins)
  DEALLOCATE (cnorm_opt)
  DEALLOCATE (corder)
!* Vector case
  IF (lg_vector) THEN
      DEALLOCATE (cg_assoc_input_field)
      DEALLOCATE (ig_assoc_input_field)
      DEALLOCATE (lrotate)
  ENDIF
!
  IF (allocated(cintmet)) DEALLOCATE(cintmet)
  IF (allocated(cextmet)) DEALLOCATE(cextmet)
  IF (allocated(neighbor)) DEALLOCATE(neighbor)
  IF (allocated(nextfl)) DEALLOCATE(nextfl)
  IF (allocated(nbofld)) DEALLOCATE(nbofld)
  IF (allocated(nbnfld)) DEALLOCATE(nbnfld)
  IF (allocated(nmapvoi)) DEALLOCATE(nmapvoi)
  IF (allocated(nmapfl)) DEALLOCATE(nmapfl)
  IF (allocated(nsubfl)) DEALLOCATE(nsubfl)
  IF (allocated(nsubvoi)) DEALLOCATE(nsubvoi)
  
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_analysis
!
!*========================================================================
SUBROUTINE dealloc_coast
!
!**** DEALLOC_COAST
!
!     Purpose:
!       Deallocate arrays defined in the "coast" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_coast".     
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


  USE mod_coast
!
!*-----------------------------------------------------------------------
!
  DEALLOCATE (npcoast)
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_coast
!
!*========================================================================
SUBROUTINE dealloc_experiment
!
!**** DEALLOC_EXPERIMENT
!
!     Purpose:
!       Deallocate arrays defined in the "experiment" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_experiment".       
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
  USE mod_experiment
!
!*-----------------------------------------------------------------------
!
  DEALLOCATE (nbcplproc)
  DEALLOCATE (nbtotproc)
  DEALLOCATE (cmodnam)
  DEALLOCATE (cmpiarg)
  DEALLOCATE (iga_unitmod)
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_experiment
!
!*========================================================================
SUBROUTINE dealloc_extrapol
!
!**** DEALLOC_EXTRAPOL
!
!     Purpose:
!       Deallocate arrays defined in the "extrapol" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_extrapol".       
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
  USE mod_extrapol
!
!*-----------------------------------------------------------------------
!
  DEALLOCATE (niwtn)
  DEALLOCATE (niwtng)
  IF (allocated(nninnfl)) DEALLOCATE(nninnfl)
  IF (allocated(nninnflg)) DEALLOCATE (nninnflg)
  DEALLOCATE (lextra)
  DEALLOCATE (lweight)
  DEALLOCATE (aextra)
  DEALLOCATE (nextra)
  DEALLOCATE (lextrapdone)
  IF (allocated(iaddress)) DEALLOCATE (iaddress)
  IF (allocated(iincre)) DEALLOCATE (iincre)
  IF (allocated(zweights)) DEALLOCATE (zweights)
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_extrapol
!
!*========================================================================
SUBROUTINE dealloc_memory
!
!**** DEALLOC_MEMORY
!
!     Purpose:
!       Deallocate arrays defined in the "memory" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_memory".       
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
  USE mod_memory
!
!*-----------------------------------------------------------------------
!
  DEALLOCATE (nsizold)
  DEALLOCATE (nsiznew)
  DEALLOCATE (nadrold)
  DEALLOCATE (nadrold_grid)
  DEALLOCATE (nadrnew)
  DEALLOCATE (nadrnew_grid)
  DEALLOCATE (mskold)
  DEALLOCATE (msknew)
  DEALLOCATE (fldold)
  DEALLOCATE (xgrold)
  DEALLOCATE (ygrold)
  DEALLOCATE (surold)
  DEALLOCATE (fldnew)
  DEALLOCATE (xgrnew)
  DEALLOCATE (ygrnew)
  DEALLOCATE (surnew)
  DEALLOCATE (nwork)
  DEALLOCATE (work)
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_memory
!
!*========================================================================
SUBROUTINE dealloc_nproc
!
!**** DEALLOC_NPROC
!
!     Purpose:
!       Deallocate arrays defined in the "nproc" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_nproc".       
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
  USE mod_nproc
!
!*-----------------------------------------------------------------------
!
  DEALLOCATE (nproc)
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_nproc
!
!*========================================================================
SUBROUTINE dealloc_parallel
!
!**** DEALLOC_ANALYSIS
!
!     Purpose:
!       Deallocate arrays defined in the "parallel" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_parallel".       
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
  USE mod_parallel
!
!*-----------------------------------------------------------------------
!
  DEALLOCATE (nparal)
  DEALLOCATE (cparal)
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_parallel
!
!*========================================================================
SUBROUTINE dealloc_pipe
!
!**** DEALLOC_PIPE
!
!     Purpose:
!       Deallocate arrays defined in the "pipe" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_pipe".       
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
  USE mod_pipe
!
!*-----------------------------------------------------------------------
!
  DEALLOCATE (cprnam)
  DEALLOCATE (cpwnam)
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_pipe
!
!*========================================================================
SUBROUTINE dealloc_rainbow
!
!**** DEALLOC_RAINBOW
!
!     Purpose:
!       Deallocate arrays defined in the "rainbow" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_rainbow".       
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
  USE mod_rainbow
!
!*-----------------------------------------------------------------------
!
  DEALLOCATE (lmapp)
  DEALLOCATE (lsubg)
  DEALLOCATE (amapp)
  DEALLOCATE (asubg)
  DEALLOCATE (nmapp)
  DEALLOCATE (nsubg)
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_rainbow
!
!*========================================================================
SUBROUTINE dealloc_sipc
!
!**** DEALLOC_SIPC
!
!     Purpose:
!       Deallocate arrays defined in the "sipc" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_sipc".       
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
  USE mod_sipc
!
!*-----------------------------------------------------------------------
!
  DEALLOCATE (mpoolidin)
  DEALLOCATE (mpoolidou)
  DEALLOCATE (mpoolinitr)
  DEALLOCATE (mpoolinitw)
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_sipc
!*========================================================================
SUBROUTINE dealloc_gsip
!
!**** DEALLOC_GSIP
!
!     Purpose:
!       Deallocate arrays defined in the "gsip" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_gsip".       
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
!       3_2-3     S. Valcke      2004/10/28   created
!
!*-----------------------------------------------------------------------
!
!** + DECLARATIONS
!
!** ++ Use of modules
!
  USE mod_gsip
!
!*-----------------------------------------------------------------------
!
  DEALLOCATE (iga_gsipw)
  DEALLOCATE (iga_gsipr)
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_gsip
!
!*========================================================================
SUBROUTINE dealloc_string
!
!**** DEALLOC_STRING
!
!     Purpose:
!       Deallocate arrays defined in the "string" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_string".       
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
  USE mod_string
  USE mod_parameter
!
!*-----------------------------------------------------------------------
!
  IF (lg_oasis_field) THEN 
     DEALLOCATE (numlab)
     DEALLOCATE (nfexch)
     DEALLOCATE (nluinp)
     DEALLOCATE (nluout)
     DEALLOCATE (nseqn)
     DEALLOCATE (nlagn)
     DEALLOCATE (cnaminp)
     DEALLOCATE (cnamout)
     DEALLOCATE (cficinp)
     DEALLOCATE (cficout)
     DEALLOCATE (cstate)
     DEALLOCATE (ig_portin_id)
     DEALLOCATE (ig_portout_id)
     DEALLOCATE (cficbf)
     DEALLOCATE (cficaf)
     DEALLOCATE (nlonbf)
     DEALLOCATE (nlatbf)
     DEALLOCATE (nlonaf)
     DEALLOCATE (nlataf)
     DEALLOCATE (ig_ntrans)
     DEALLOCATE (canal)
     DEALLOCATE (ig_grid_nbrbf)
     DEALLOCATE (ig_grid_nbraf)
  ENDIF

  DEALLOCATE (ig_lag)
  DEALLOCATE (lg_state)
  DEALLOCATE (ig_no_rstfile)
  DEALLOCATE (cg_name_rstfile)
  DEALLOCATE (ig_numlab)
  DEALLOCATE (ig_freq)
  DEALLOCATE (ig_total_nseqn)
  DEALLOCATE (cg_input_field)
  DEALLOCATE (cg_output_field)
  DEALLOCATE (ig_total_state)
  DEALLOCATE (ig_local_trans)
  DEALLOCATE (cg_input_file)
  DEALLOCATE (ig_number_field)
  DEALLOCATE (ig_total_ntrans)
  DEALLOCATE (cg_restart_file)
  DEALLOCATE (cga_locatorbf)
  DEALLOCATE (cga_locatoraf)
  DEALLOCATE (ig_invert)
  DEALLOCATE (ig_reverse)
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_string
!
!*========================================================================
SUBROUTINE dealloc_timestep
!
!**** DEALLOC_TIMESTEP
!
!     Purpose:
!       Deallocate arrays defined in the "timestep" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_timestep".       
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
  USE mod_timestep
!
!*-----------------------------------------------------------------------
!
  DEALLOCATE (mstep)
  DEALLOCATE (mfcpl)
  DEALLOCATE (mdt)
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_timestep
!
!*========================================================================
SUBROUTINE dealloc_unitncdf
!
!**** DEALLOC_UNITNCDF
!
!     Purpose:
!       Deallocate arrays defined in the "unitncdf" module
!
!     Interface:
!       none
!    
!     Method:
!       Deallocation of arrays allocated in "alloc_unitncdf".       
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
  USE mod_unitncdf
!
!*-----------------------------------------------------------------------
!
  DEALLOCATE (nc_inpid)
  DEALLOCATE (nc_outid)
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE dealloc_unitncdf
!
!*========================================================================
