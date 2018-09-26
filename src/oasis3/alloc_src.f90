!                       ***************************
!                       *  ALLOCATION SUBROUTINES *
!                       ***************************
SUBROUTINE alloc_anais1
!
!**** ALLOC_ANAIS
!
!     Purpose:
!       Allocate arrays defined in the "anais" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_anais
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (varmul(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout &
     ('Error in "varmul"allocation of anais module',il_err,1)
  varmul(:)=0
  ALLOCATE (niwtm(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "niwtm"allocation of anais module',il_err,1)
  niwtm(:)=0
  ALLOCATE (niwtg(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "niwtg"allocation of anais module',il_err,1)
  niwtg(:)=0
  ALLOCATE (linit(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "linit"allocation of anais module',il_err,1)
  linit(:)=.false.
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_anais1
!
!*========================================================================
!
SUBROUTINE alloc_anais2
!
!**** ALLOC_ANAIS
!
!     Purpose:
!       Allocate arrays defined in the "anais" module
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
!       3.0       S. Valcke       2004/01/05  created
!
!*-----------------------------------------------------------------------
!
!** + DECLARATIONS
!
!** ++ Use of modules
!
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_anais
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (ngint(ig_maxnoa*ig_maxnfg*ig_maxgrd), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ngint"allocation of anais module',il_err,1)
  ngint(:)=0
  ALLOCATE (nmint(ig_maxwoa*ig_maxnfm*ig_maxgrd), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nmint"allocation of anais module',il_err,1)
  nmint(:)=0
  ALLOCATE (nmesh(ig_maxnfm*ig_maxgrd), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nmesh"allocation of anais module',il_err,1)
  nmesh(:)=0
  ALLOCATE (agint(ig_maxnoa*ig_maxnfg*ig_maxgrd), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "agint"allocation of anais module',il_err,1)
  agint(:)=0
  ALLOCATE (amint(ig_maxwoa*ig_maxnfm*ig_maxgrd), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "amint"allocation of anais module',il_err,1)
  amint(:)=0
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_anais2
!
!*========================================================================
SUBROUTINE alloc_analysis
!
!**** ALLOC_ANALYSIS
!
!     Purpose:
!       Allocate arrays defined in the "analysis" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_analysis
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (ncofld(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ncofld"allocation of analysis module',il_err,1)
  ncofld(:)=0
  ALLOCATE (neighborg(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "neighborg"allocation of analysis module',il_err,1)
  neighborg(:)=0
  ALLOCATE (nludat(ig_maxcomb,ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nludat"allocation of analysis module',il_err,1)
  nludat(:,:)=0
  ALLOCATE (nlufil(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nlufil"allocation of analysis module',il_err,1)
  nlufil(:)=0
  ALLOCATE (nlumap(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nlumap"allocation of analysis module',il_err,1)
  nlumap(:)=0
  ALLOCATE (nlusub(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nlusub"allocation of analysis module',il_err,1)
  nlusub(:)=0
  ALLOCATE (nluext(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nluext"allocation of analysis module',il_err,1)
  nluext(:)=0
  ALLOCATE (nosper(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nosper"allocation of analysis module',il_err,1)
  nosper(:)=0
  ALLOCATE (notper(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "notper"allocation of analysis module',il_err,1)
  notper(:)=0
  ALLOCATE (ntinpflx(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ntinpflx"allocation of analysis module',il_err,1)
  ntinpflx(:)=0
  ALLOCATE (ntoutflx(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ntoutflx"allocation of analysis module',il_err,1)
  ntoutflx(:)=0
  ALLOCATE (amskval(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "amskval"allocation of analysis module',il_err,1)
  amskval(:)=0
  ALLOCATE (amskvalnew(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "amskvalnew"allocation of analysis module',il_err,1)
  amskvalnew(:)=0
  ALLOCATE (acocoef(ig_maxcomb,ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "acocoef"allocation of analysis module',il_err,1)
  acocoef(:,:)=0
  ALLOCATE (abocoef(ig_maxcomb,ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "abocoef"allocation of analysis module',il_err,1)
  abocoef(:,:)=0
  ALLOCATE (abncoef(ig_maxcomb,ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "abncoef"allocation of analysis module',il_err,1)
  abncoef(:,:)=0
  ALLOCATE (afldcoef(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "afldcoef"allocation of analysis module',il_err,1)
  afldcoef(:)=0
  ALLOCATE (afldcobo(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "afldcobo"allocation of analysis module',il_err,1)
  afldcobo(:)=0
  ALLOCATE (afldcobn(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "afldcobn"allocation of analysis module',il_err,1)
  afldcobn(:)=0
  ALLOCATE (cxordbf(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cxordbf"allocation of analysis module',il_err,1)
  cxordbf(:)=' '
  ALLOCATE (cyordbf(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cyordbf"allocation of analysis module',il_err,1)
  cyordbf(:)=' '
  ALLOCATE (cxordaf(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cxordaf"allocation of analysis module',il_err,1)
  cxordaf(:)=' '
  ALLOCATE (cyordaf(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cyordaf"allocation of analysis module',il_err,1)
  cyordaf(:)=' '
  ALLOCATE (cgrdtyp(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cgrdtyp"allocation of analysis module',il_err,1)
  cgrdtyp(:)=' '
  ALLOCATE (cfldtyp(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cfldtyp"allocation of analysis module',il_err,1)
  cfldtyp(:)=' '
  ALLOCATE (cfilfic(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cfilfic"allocation of analysis module',il_err,1)
  cfilfic(:)=' '
  ALLOCATE (cfilmet(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cfilmet"allocation of analysis module',il_err,1)
  cfilmet(:)=' '
  ALLOCATE (cconmet(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cconmet"allocation of analysis module',il_err,1)
  cconmet(:)=' '
  ALLOCATE (cfldcoa(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cfldcoa"allocation of analysis module',il_err,1)
  cfldcoa(:)=' '
  ALLOCATE (cfldfin(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cfldfin"allocation of analysis module',il_err,1)
  cfldfin(:)=' '
  ALLOCATE (ccofld(ig_maxcomb,ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ccofld"allocation of analysis module',il_err,1)
  ccofld(:,:)=' '
  ALLOCATE (cbofld(ig_maxcomb,ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cbofld"allocation of analysis module',il_err,1)
  cbofld(:,:)=' '
  ALLOCATE (cbnfld(ig_maxcomb,ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout &
     ('Error in "cbnfld"allocation of analysis module',il_err,1)
  cbnfld(:,:)=' '
  ALLOCATE (ccofic(ig_maxcomb,ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ccofic"allocation of analysis module',il_err,1)
  ccofic(:,:)=' '
  ALLOCATE (cdqdt(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cdqdt"allocation of analysis module',il_err,1)
  cdqdt(:)=' '
  ALLOCATE (cgrdmap(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cgrdmap"allocation of analysis module',il_err,1)
  cgrdmap(:)=' '
  ALLOCATE (cmskrd(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cmskrd"allocation of analysis module',il_err,1)
  cmskrd(:)=' '
  ALLOCATE (cgrdsub(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cgrdsub"allocation of analysis module',il_err,1)
  cgrdsub(:)=' '
  ALLOCATE (ctypsub(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ctypsub"allocation of analysis module',il_err,1)
  ctypsub(:)=' '
  ALLOCATE (cgrdext(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout &
     ('Error in "cgrdext"allocation of analysis module',il_err,1)
  cgrdext(:)=' '
  ALLOCATE (csper(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "csper"allocation of analysis module',il_err,1)
  csper(:)=' '
  ALLOCATE (ctper(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ctper"allocation of analysis module',il_err,1)
  ctper(:)=' '
  ALLOCATE (lsurf(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "lsurf"allocation of analysis module',il_err,1)
  lsurf(:)=.false.
  ALLOCATE (nscripvoi(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in nscripvoi allocation of analysis module',il_err,1)
  nscripvoi(:)=0
! 
!* Alloc array needed for SCRIP 
!
  ALLOCATE (cmap_method(ig_nfield),stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cmap_method" allocation of inipar_alloc',il_err,1)
  cmap_method(:)=' '
  ALLOCATE (cfldtype(ig_nfield),stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cfldtype"allocation of inipar_alloc',il_err,1)
  cfldtype(:)=' '
  ALLOCATE (crsttype(ig_nfield),stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "crsttype"allocation of inipar_alloc',il_err,1)
  crsttype(:)=' '
  ALLOCATE (nbins(ig_nfield),stat=il_err)
  IF (il_err.NE.0) CALL prtout &
     ('Error in "nbins"allocation of inipar_alloc',il_err,1)
  nbins(:)=0
  ALLOCATE (cnorm_opt(ig_nfield),stat=il_err)
  IF (il_err.NE.0) CALL prtout &
     ('Error in "cnorm_opt"allocation of inipar_alloc',il_err,1)
  cnorm_opt(:)=' '
  ALLOCATE (corder(ig_nfield),stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "corder"allocation of inipar_alloc',il_err,1)
  corder(:)=' '
!
!Vector case:
!
  IF (lg_vector) THEN
      ALLOCATE (cg_assoc_input_field(ig_total_nfield),stat=il_err)
      IF (il_err.NE.0) CALL prtout & 
     ('Error in "cg_assoc_input_field"allocation of inipar_alloc',il_err,1)
      cg_assoc_input_field(:)=' '
      ALLOCATE (ig_assoc_input_field(ig_total_nfield),stat=il_err)
      IF (il_err.NE.0) CALL prtout & 
     ('Error in "ig_assoc_input_field"allocation of inipar_alloc',il_err,1)
      ig_assoc_input_field(:)=0
      ALLOCATE (lrotate(ig_total_nfield),stat=il_err)
      IF (il_err.NE.0) CALL prtout & 
     ('Error in "lrotate"allocation of inipar_alloc',il_err,1)
  ENDIF
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_analysis
!
!*========================================================================
SUBROUTINE alloc_coast
!
!**** ALLOC_COAST
!
!     Purpose:
!       Allocate arrays defined in the "coast" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_coast
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (npcoast(ig_maxwoa,6), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "npcoast"allocation of coast module',il_err,1)
  npcoast(:,:)=0
  
!
!* ----------------------------------------------------------------------
!
END SUBROUTINE alloc_coast
!
!*=======================================================================
SUBROUTINE alloc_experiment
!
!**** ALLOC_EXPERIMENT
!
!     Purpose:
!       Allocate arrays defined in the "experiment" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_experiment
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (nbcplproc(ig_nmodel), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nbcplproc"allocation of experiment module',il_err,1)
  nbcplproc(:)=0
  ALLOCATE (nbtotproc(ig_nmodel), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nbtotproc"allocation of experiment module',il_err,1)
  nbtotproc(:)=0
  ALLOCATE (cmodnam(ig_nmodel), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cmodnam"allocation of experiment module',il_err,1)
  cmodnam(:)=' '
  ALLOCATE (cmpiarg(ig_nmodel), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cmpiarg"allocation of experiment module',il_err,1)
  cmpiarg(:)=' '
  ALLOCATE (iga_unitmod(ig_nmodel), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in iga_unitmod allocation of experiment module',il_err,1)
  iga_unitmod(:)=0
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_experiment
!
!*========================================================================
SUBROUTINE alloc_extrapol1
!
!**** ALLOC_EXTRAPOL
!
!     Purpose:
!       Allocate arrays defined in the "extrapol" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_extrapol
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (niwtn(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "niwtn"allocation of extrapol module',il_err,1)
  niwtn(:)=0
  ALLOCATE (niwtng(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "niwtng"allocation of extrapol module',il_err,1)
  niwtng(:)=0
  ALLOCATE (lextra(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "lextra"allocation of extrapol module',il_err,1)
  lextra(:)=.false.
  ALLOCATE (lweight(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "lweight"allocation of extrapol module',il_err,1)
  lweight(:)=.false.
  ALLOCATE (lextrapdone(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout &
     ('Error in "lextrapdone" allocation',il_err,1)
  lextrapdone(:) = .FALSE.
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_extrapol1
!
!*========================================================================
!
SUBROUTINE alloc_extrapol2
!
!**** ALLOC_EXTRAPOL
!
!     Purpose:
!       Allocate arrays defined in the "extrapol" module
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
!       3.0       S. Valcke      2004/01/05  created
!
!*-----------------------------------------------------------------------
!
!** + DECLARATIONS
!
!** ++ Use of modules
!
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_extrapol
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (aextra(ig_maxext*ig_maxnbn*ig_maxgrd), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "aextra"allocation of extrapol module',il_err,1)
  aextra(:)=0
  ALLOCATE (nextra(ig_maxext*ig_maxnbn*ig_maxgrd), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nextra"allocation of extrapol module',il_err,1)
  nextra(:)=0
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_extrapol2
!
!*========================================================================
SUBROUTINE alloc_memory1
!
!**** ALLOC_MEMORY
!
!     Purpose:
!       Allocate arrays defined in the "memory" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_memory
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (nsizold(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout &  
     ('Error in "nsizold"allocation of memory module',il_err,1)
  nsizold(:)=0
  ALLOCATE (nsiznew(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nsiznew"allocation of memory module',il_err,1)
  nsiznew(:)=0
  ALLOCATE (nadrold(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nadrold"allocation of memory module',il_err,1)
  nadrold(:)=0
  ALLOCATE (nadrold_grid(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nadrold_grid"allocation of memory module',il_err,1)
  nadrold_grid(:)=0
  ALLOCATE (nadrnew(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nadrnew"allocation of memory module',il_err,1)
  nadrnew(:)=0
  ALLOCATE (nadrnew_grid(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nadrnew_grid"allocation of memory module',il_err,1)
  nadrnew_grid(:)=0
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_memory1
!
!*========================================================================
SUBROUTINE alloc_memory2
!
!**** ALLOC_MEMORY
!
!     Purpose:
!       Allocate arrays defined in the "memory" module
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
!       3.0       S. Valcke       2004/01/05  created
!
!*-----------------------------------------------------------------------
!
!** + DECLARATIONS
!
!** ++ Use of modules
!
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_memory
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (mskold(ig_maxold_grid), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "mskold"allocation of memory module',il_err,1)
  mskold(:)=0
  ALLOCATE (msknew(ig_maxnew_grid), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "msknew"allocation of memory module',il_err,1)
  msknew(:)=0
  ALLOCATE (fldold(ig_maxold), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "fldold"allocation of memory module',il_err,1)
  fldold(:)=0
  ALLOCATE (xgrold(ig_maxold_grid), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "xgrold"allocation of memory module',il_err,1)
  xgrold(:)=0
  ALLOCATE (ygrold(ig_maxold_grid), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ygrold"allocation of memory module',il_err,1)
  ygrold(:)=0
  ALLOCATE (surold(ig_maxold_grid), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "surold"allocation of memory module',il_err,1)
  surold(:)=0
  ALLOCATE (fldnew(ig_maxnew), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "fldnew"allocation of memory module',il_err,1)
  fldnew(:)=0
  ALLOCATE (xgrnew(ig_maxnew_grid), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "xgrnew"allocation of memory module',il_err,1)
  xgrnew(:)=0
  ALLOCATE (ygrnew(ig_maxnew_grid), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ygrnew"allocation of memory module',il_err,1)
  ygrnew(:)=0
  ALLOCATE (surnew(ig_maxnew_grid), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "surnew"allocation of memory module',il_err,1)
  surnew(:)=0
  ALLOCATE (nwork(ig_nwork), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nwork"allocation of memory module',il_err,1)
  nwork(:)=0
  ALLOCATE (work(ig_work), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "work"allocation of memory module',il_err,1)
  work(:)=0
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_memory2
!
!*========================================================================
SUBROUTINE alloc_nproc
!
!**** ALLOC_NPROC
!
!     Purpose:
!       Allocate arrays defined in the "nproc" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_nproc
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (nproc(ig_nmodel), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nproc"allocation of nproc module',il_err,1)
  nproc(:)=0
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_nproc
!
!*========================================================================
SUBROUTINE alloc_parallel
!
!**** ALLOC_ANALYSIS
!
!     Purpose:
!       Allocate arrays defined in the "parallel" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_parallel
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (nparal(3,ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nparal"allocation of parallel module',il_err,1)
  nparal(:,:)=0
  ALLOCATE (cparal(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cparal"allocation of parallel module',il_err,1)
  cparal(:)=' '  
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_parallel
!
!*========================================================================
SUBROUTINE alloc_pipe
!
!**** ALLOC_PIPE
!
!     Purpose:
!       Allocate arrays defined in the "pipe" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_pipe
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (cprnam(ig_nmodel), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cprnam"allocation of pipe module',il_err,1)
  cprnam(:)=' '
  ALLOCATE (cpwnam(ig_nmodel), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cpwnam"allocation of pipe module',il_err,1)
  cpwnam(:)=' '
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_pipe
!
!*========================================================================
SUBROUTINE alloc_rainbow1
!
!**** ALLOC_RAINBOW
!
!     Purpose:
!       Allocate arrays defined in the "rainbow" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_rainbow
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (lmapp(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "lmapp"allocation of rainbow module',il_err,1)
  lmapp(:)=.false.
  ALLOCATE (lsubg(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "lsubg"allocation of rainbow module',il_err,1)
  lsubg(:)=.false.
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_rainbow1
!
!*========================================================================
SUBROUTINE alloc_rainbow2
!
!**** ALLOC_RAINBOW
!
!     Purpose:
!       Allocate arrays defined in the "rainbow" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_rainbow
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (amapp(ig_maxmoa*ig_maxnfp*ig_maxgrd), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "amapp"allocation of rainbow module',il_err,1)
  amapp(:)=0
  ALLOCATE (asubg(ig_maxsoa*ig_maxnfs*ig_maxgrd), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "asubg"allocation of rainbow module',il_err,1)
  asubg(:)=0
  ALLOCATE (nmapp(ig_maxmoa*ig_maxnfp*ig_maxgrd), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nmapp"allocation of rainbow module',il_err,1)
  nmapp(:)=0
  ALLOCATE (nsubg(ig_maxsoa*ig_maxnfs*ig_maxgrd), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "asubg"allocation of rainbow module',il_err,1)
  nsubg(:)=0
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_rainbow2
!
!*========================================================================
SUBROUTINE alloc_sipc
!
!**** ALLOC_SIPC
!
!     Purpose:
!       Allocate arrays defined in the "sipc" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_sipc
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (mpoolidin(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "mpoolidin"allocation of sipc module',il_err,1)
  mpoolidin(:)=0
  ALLOCATE (mpoolidou(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "mpoolidou"allocation of sipc module',il_err,1)
  mpoolidou(:)=0
  ALLOCATE (mpoolinitr(ig_nmodel), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "mpoolinitr"allocation of sipc module',il_err,1)
  mpoolinitr(:)=0
  ALLOCATE (mpoolinitw(ig_nmodel), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "mpoolinitw"allocation of sipc module',il_err,1)
  mpoolinitw(:)=0
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_sipc
!
!*========================================================================
SUBROUTINE alloc_gsip
!
!**** ALLOC_GSIP
!
!     Purpose:
!       Allocate arrays defined in the "gsip" module
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
!       3_2-4     S. Valcke      2004/10/15  created
!
!*-----------------------------------------------------------------------
!
!** + DECLARATIONS
!
!** ++ Use of modules
!
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_experiment
  USE mod_gsip
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err, il_totproc
!
!*-----------------------------------------------------------------------
!
  il_totproc = sum ( nbtotproc(:) )
  ALLOCATE (iga_gsipw(il_totproc), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "iga_gsipw" allocation of gsip module',il_err,1)
  iga_gsipw(:)=0
  ALLOCATE (iga_gsipr(il_totproc), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "iga_gsipr" allocation of gsip module',il_err,1)
  iga_gsipr(:)=0
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_gsip
!
!*========================================================================
SUBROUTINE alloc_string
!
!**** ALLOC_STRING
!
!     Purpose:
!       Allocate arrays defined in the "string" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_string
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (cg_name_rstfile(ig_nbr_rstfile), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cg_name_rstfile"allocation of string module',il_err,1)
  cg_name_rstfile(:)=' '
  ALLOCATE (ig_lag(ig_total_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ig_lag"allocation of string module',il_err,1) 
  ig_lag(:)=0
  ALLOCATE (ig_no_rstfile(ig_total_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ig_no_rstfile"allocation of string module',il_err,1)
  ig_no_rstfile(:)=1
  ALLOCATE (cg_input_field(ig_total_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cg_input_field"allocation of string module',il_err,1)
  cg_input_field(:)=' '
  ALLOCATE (ig_numlab(ig_total_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ig_numlab"allocation of string module',il_err,1)
  ig_numlab(:)=0
  ALLOCATE (ig_freq(ig_total_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ig_freq"allocation of string module',il_err,1)
  ig_freq(:)=0
  ALLOCATE (ig_total_nseqn(ig_total_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ig_total_nseqn"allocation of string module',il_err,1)
  ig_total_nseqn(:)=0
  ALLOCATE (ig_local_trans(ig_total_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ig_local_trans"allocation of string module',il_err,1)
  ig_local_trans(:)=0
  ALLOCATE (ig_invert(ig_total_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ig_invert" allocation of string module',il_err,1) 
  ig_invert(:)=0
  ALLOCATE (ig_reverse(ig_total_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ig_reverse" allocation of string module',il_err,1) 
  ig_reverse(:)=0
!
!** + Allocate following arrays only if one field (at least) goes
!     through Oasis
!
  IF (lg_oasis_field) THEN
  ALLOCATE (numlab(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "numlab"allocation of string module',il_err,1)
  numlab(:)=0
  ALLOCATE (nfexch(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nfexch"allocation of string module',il_err,1)
  nfexch(:)=0
  ALLOCATE (nluinp(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nluinp"allocation of string module',il_err,1)
  nluinp(:)=0
  ALLOCATE (nluout(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nluout"allocation of string module',il_err,1)
  nluout(:)=0
  ALLOCATE (nseqn(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nseqn"allocation of string module',il_err,1)
  nseqn(:)=0
  ALLOCATE (nlagn(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nlagn" allocation of string module',il_err,1)
  nlagn(:)=0
  ALLOCATE (cnaminp(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cnaminp"allocation of string module',il_err,1)
  cnaminp(:)=' '
  ALLOCATE (cnamout(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cnamout"allocation of string module',il_err,1)
  cnamout(:)=' '
  ALLOCATE (cficout(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cficout"allocation of string module',il_err,1)
  cficout(:)=' '
  ALLOCATE (cstate(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "cstate"allocation of string module',il_err,1)
  cstate(:)=' '
  ALLOCATE (ig_portin_id(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ig_portin_id"allocation of string module',il_err,1)
  ig_portin_id(:)=0
  ALLOCATE (ig_portout_id(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "ig_portout_id"allocation of string module',il_err,1)
  ig_portout_id(:)=0
  ENDIF
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_string
!
!*========================================================================
SUBROUTINE alloc_timestep
!
!**** ALLOC_TIMESTEP
!
!     Purpose:
!       Allocate arrays defined in the "timestep" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_timestep
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (mstep(ig_nmodel), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "mstep"allocation of timestep module',il_err,1)
  mstep(:)=0
  ALLOCATE (mfcpl(ig_nmodel), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "mfcpl"allocation of timestep module',il_err,1)
  mfcpl(:)=0
  ALLOCATE (mdt(ig_nmodel), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "mdt"allocation of timestep module',il_err,1)
  mdt(:)=0
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_timestep
!
!*========================================================================
SUBROUTINE alloc_unitncdf
!
!**** ALLOC_UNITNCDF
!
!     Purpose:
!       Allocate arrays defined in the "unitncdf" module
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
  USE mod_kinds_oasis
  USE mod_parameter
  USE mod_unitncdf
!
!** ++ Local declarations
!
  INTEGER (kind=ip_intwp_p) :: il_err
!
!*-----------------------------------------------------------------------
!
  ALLOCATE (nc_inpid(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nc_inpid"allocation of unitncdf module',il_err,1)
  nc_inpid(:)=0
  ALLOCATE (nc_outid(ig_nfield), stat=il_err)
  IF (il_err.NE.0) CALL prtout & 
     ('Error in "nc_outid"allocation of unitncdf module',il_err,1)
  nc_outid(:)=0
!
!*-----------------------------------------------------------------------
!
END SUBROUTINE alloc_unitncdf
!
!*========================================================================




