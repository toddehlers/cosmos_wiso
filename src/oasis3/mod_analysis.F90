MODULE mod_analysis
!
! -- analysis.h   18-08-95   Version 2.0beta   Author: Laurent Terray
!    **********   01-02-96   Version 2.0 : addition of cdqdt,cgrdmap,nlumap
!                                          and nmapfl
!                 05-08-96   Version 2.1 : addition of nlusub, cgrdsub 
!                                          ctypsub and nsubfl (subgrid)
!                 29-08-96   Version 2.1 : new input for subgrid, mozaic
!                                          and anais, maximum number of
!                                          neighbors. New input for redglo
!                 14-04-99   Version 2.3 : new input for extrap called by
!                                          glored. Change periodicity 
!                                          variables.
!                 05-09-00   Version 2.5 : integral flags for check[in-out]
!                 18-03-02   Version 2.5 : changed in module
!                 27-06-02   Version 2.5 : scrip added
!@
!@  Contents : variables related to the set of analysis for each field
!@  --------
!@
!@ Mask :
!@
!@ -- amskval : mask value for the imported field
!@
!@ MaskP :
!@
!@ -- amskvalnew : mask value for the exported field
!@
!@
!@ Mozaic :
!@
!@ -- cgrdmap : file name for grid mapping used in mozaic (1D)
!@
!@ -- nlumap  : logical units associated to previous files (1D)
!@
!@ -- nmapfl  : mapping dataset identificator number (1D)
!@
!@ -- nmapvoi : maximum number of neighbors (1D)
!@
!@ Invert :
!@
!@ -- cxordbf : field ordering (longitude) before interpolation (1D)
!@
!@ -- cyordbf : field ordering (latitude) before interpolation (1D)
!@
!@ Reverse :
!@
!@ -- cxordaf : field ordering (longitude) after interpolation (1D)
!@
!@ -- cyordaf : field ordering (latitude) after interpolation (1D)
!@
!@ Extrap :
!@
!@ -- cextmet : extrapolation method (1D)
!@
!@ -- neighbor : number of neighbors used in extrapolation (1D)
!@
!@ -- neighborg : number of neighbors used in extrapolation when extrap
!@                is called by GLORED (1D)
!@
!@ -- cgrdext : file names for data used in extrapolation WEIGHT (1D)
!@
!@ -- nluext  : logical units associated to previous files (1D)
!@
!@ -- nextfl  : extrapolation dataset identificator number (1D)
!@
!@ Interp :
!@
!@ -- cintmet : interpolation method (1D)
!@
!@ -- cgrdtyp : source grid type (1D)
!@
!@ -- csper   : source grid periodicity type (P=Periodic, R=Regional)
!@
!@ -- ctper   : target grid periodicity type (P=Periodic, R=Regional)
!@
!@ -- cfldtyp : field type (scalar or vector) (1D)
!@
!@ Scrip :
!@
!@ -- cmap_method : scrip remapping method ('CONSERV', 'BILINEAR', 'BICUBIC',
!@                                          'GAUSWGT'or 'DISTWGT')
!@ -- cfldtype    : field type ('SCALAR' or 'VECTOR_I' or 'VECTOR_J')
!@
!@ -- cnorm_opt   : normalization option ('FRACAREA', 'DESTAREA' or 'NONE')
!@
!@ -- crsttype    : restriction type for SCRIP search ('LATITUDE', 'LATLON',
!@                                                     or 'REDUCED')
!@
!@ -- nbins       : number of search bins
!@
!@ -- corder      : remapping order ('FIRST' or 'SECOND')
!@
!@ --  nscripvoi  : number of neighbour for 'GAUSWGT'or 'DISTWGT'
!@
!@ -- ig_assoc_input_field : number id of the associated field for vector case
!@
!@ -- cg_assoc_input_field : name of the associated field for vector case
!@
!@ -- lrotate :  logical for the rotation to cartesian referential in vector case
!@
!@ Filling :
!@
!@ -- cfilfic : file name for climatological field to complete model field (1D)
!@
!@ -- nlufil  : logical units connected to previous files (1D)
!@
!@ -- cficmet : filling method (1D)
!@
!@ -- cfldcor : field name for flux correction term due to SST filling
!@
!@ -- nlucor  : logical unit used to write flux correction term
!@
!@ !onserv :
!@
!@ -- cconmet : conservation method (1D)
!@
!@ Glored, redglo :
!@
!@ -- ntronca : gaussian troncature for reduced <-> global gaussian grid (1D)
!@
!@ -- cmskrd : extrapolation flag to handle processing reduced grid values (1D)
!@
!@ Correct :
!@
!@ -- afldcoef : main field multiplicative coefficient (1D)
!@
!@ -- ncofld : number of additional fields in correction formula (1D)
!@
!@ -- ccofic : file names for additional data files (2D)
!@
!@ -- nludat : logical units associated to previous files (2D)
!@
!@ -- ccofld : symbolic names for additional fields (2D)
!@
!@ -- acocoef : multiplicative coefficients for additional fields (2D)
!@
!@ Blasold :
!@
!@ -- afldcobo : main field multiplicative coefficient (1D)
!@
!@ -- nbofld : number of additional fields in linear combination formula (1D)
!@
!@ -- cbofld : symbolic names for additional fields (2D)
!@
!@ -- abocoef : multiplicative coefficients for additional fields (2D)
!@
!@ Blasnew :
!@
!@ -- afldcobn : main field multiplicative coefficient (1D)
!@
!@ -- nbnfld : number of additional fields in linear combination formula (1D)
!@
!@ -- cbnfld : symbolic names for additional fields (2D)
!@
!@ -- abncoef : multiplicative coefficients for additional fields (2D)
!@
!@ Subgrid :
!@
!@ -- cgrdsub : file name for subgrid data used in subgrid (1D)
!@
!@ -- nlusub  : logical units associated to previous files (1D)
!@
!@ -- ctypsub : type of subgrid interpolation (solar or non solar) (1D)
!@
!@ -- nsubfl  : subgrid dataset identificator number (1D)
!@
!@ -- nsubvoi : maximum number of neighbors (1D)
!@
!@ -- cfldcoa : coarse grid field name (1D)
!@
!@ -- cfldfin : fine grid field name (1D)
!@
!@ -- cdqdt   : first order taylor coefficient (1D)
!@
!@ Checkin :
!@
!@ -- ntinpflx : input field integral flag
!@
!@ Checkout :
!@
!@ -- ntoutflx : output field integral flag
!@
!@ Checkin, Checkout, Conserv, Gaussian
!@
!@ -- lsurf: Logical indicating if grid surface information has to be given 
!@            in auxiliary file areas or areas.nc for each field
!@ -- lg_areas: Logical indicating if grid surface file has to exist.
!@
!@ SCRIPR
!@
!@ -- lg_vector: Logical indicating if vector fields are treated
!@
!     -------------------------------------------------------------------
!
  USE mod_kinds_oasis
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: neighbor, ntronca, ncofld
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: neighborg, nbofld, nbnfld
  INTEGER (kind=ip_intwp_p), DIMENSION(:,:), ALLOCATABLE :: nludat
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: nlufil, nlumap 
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: nmapfl, nmapvoi, nlusub
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: nsubfl, nsubvoi, nluext
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE :: nextfl, nosper, notper
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE ::  ntinpflx, ntoutflx, nbins
  INTEGER (kind=ip_intwp_p) :: nlucor
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE ::  nscripvoi
  INTEGER (kind=ip_intwp_p), DIMENSION(:), ALLOCATABLE ::  ig_assoc_input_field
  !
  REAL (kind=ip_realwp_p), DIMENSION(:), ALLOCATABLE :: amskval, amskvalnew
  REAL (kind=ip_realwp_p), DIMENSION(:,:), ALLOCATABLE :: acocoef, abocoef, abncoef
  REAL (kind=ip_realwp_p), DIMENSION(:), ALLOCATABLE :: afldcoef, afldcobo, afldcobn
  !
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: cxordbf, cyordbf
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: cxordaf, cyordaf
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: cextmet, cintmet
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: cgrdtyp, cfldtyp
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: cfilfic, cfilmet 
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: cconmet, cfldcoa
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: cfldfin
  CHARACTER(len=8), DIMENSION(:,:),ALLOCATABLE :: ccofld, cbofld
  CHARACTER(len=8), DIMENSION(:,:),ALLOCATABLE :: cbnfld, ccofic
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: cdqdt, cgrdmap
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: cmskrd, cgrdsub
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: ctypsub, cgrdext
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: csper,ctper
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: cmap_method, corder
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: cnorm_opt, cfldtype
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: crsttype
  CHARACTER(len=8) :: cfldcor
  CHARACTER(len=8), DIMENSION(:),ALLOCATABLE :: cg_assoc_input_field
  !      
  LOGICAL, DIMENSION(:),ALLOCATABLE :: lsurf, lrotate
  LOGICAL :: lg_areas 
  LOGICAL :: lg_vector  

  !    -------------------------------------------------------------------

END MODULE mod_analysis




