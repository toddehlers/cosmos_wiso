MODULE mod_parameter
!
! -- parameter.h   01-09-95   Version 2.0   Author: Laurent Terray
!    ***********   05-08-96   Version 2.1 
!                             Mods: Add parameters for mapping 
!                                   and subgrid (L. Terray)
!                  20-06-97   Version 2.2 Mods: Add jpbyteint (S. Valcke)
!                  16-12-97                     Add jpext and jpnbn (L. Terray)
!                  31-12-97                     Suppress jpwrk (L. Terray)
!                  12-10-99   Version 2.3 Mods: Add jpnfn (S. Valcke)
!                  28-01-00   Version 2.4 Mods :jpparal=102, jpfield=10 (J.Latour)
!                  18-03-02   Version 2.5 changed in module (A.Caubel)
!@
!@  Contents : parameter file for OASIS
!@  --------
!@
!@ -- ig_total_nfield : number of coupling fields
!@
!@ -- ig_direct_nfield : number of fields not going through Oasis
!@
!@ -- ig_nfield : number of fields going through Oasis
!@
!@ -- lg_oasis_field : logical indicating if all fields go through Oasis or not
!@                 (.false. if all fields are exchanged directly)
!@
!@ -- ig_nmodel : maximum number of models
!@
!@ -- ig_maxcomb : maximum number of fields to be combined in the BLASxxx analyses
!@
!@ -- ig_maxold : Memory size of the macro arrays handling fields values
!@
!@ -- ig_maxnew : Memory size of the macro arrays handling fields values 
!@
!@ -- ig_maxold_grid : Memory size of the arrays handling field grid-related 
!@                     data before interpolation
!@
!@ -- ig_maxnew_grid : Memory size of the arrays handling field grid-related 
!@                     data after interpolation
!@
!@ -- ig_nwork : Size of "nwork" array
!@
!@ -- ig_work : Size of "work" array
!@
!@ -- ig_maxgrd : maximum grid size of indirect fields
!@
!@ -- ig_total_maxgrd : maximum grid size of all fields defined in the 
!@                      namcouple (direct and indirect fields)
!@
!@ -- ig_maxwoa : maximum number of underlying neighbors for SURFMESH interpolation
!@
!@ -- ig_maxnoa : number of neighbors for GAUSSIAN interpolation
!@
!@ -- ig_maxmoa : maximum number of underlying neighbors for MOZAIC interpolation
!@
!@ -- ig_maxsoa : maximum number of overlaying neighbors for SUBGRID interpolation
!@
!@ -- ig_maxnfm : maximum number of different SURFMESH interpolations
!@
!@ -- ig_maxnfg : maximum number of different GAUSSIAN interpolations
!@
!@ -- ig_maxnfp : maximum number of different MOZAIC interpolations
!@
!@ -- ig_maxnfs : maximum number of different SUBGRID interpolations
!@
!@ -- ig_maxnfn : maximum number of different NINENN extrapolations
!@
!@ -- ig_maxext : maximum number of neighbors for extrapolation
!@
!@ -- ig_maxnbn : maximum number of different extrapolation
!@
!     -------------------------------------------------------------------
!
!* Useful numerical values
  USE mod_kinds_oasis
  INTEGER (kind=ip_intwp_p),PARAMETER :: jpeight = 8, jpfour = 4, jpeighty = 80
!
!-----Field status
!
  INTEGER(kind=ip_intwp_p), PARAMETER :: ip_exported  = 1
  INTEGER(kind=ip_intwp_p), PARAMETER :: ip_ignored   = 2
  INTEGER(kind=ip_intwp_p), PARAMETER :: ip_input     = 3
  INTEGER(kind=ip_intwp_p), PARAMETER :: ip_output    = 4
  INTEGER(kind=ip_intwp_p), PARAMETER :: ip_expout    = 5
  INTEGER(kind=ip_intwp_p), PARAMETER :: ip_ignout    = 6
  INTEGER(kind=ip_intwp_p), PARAMETER :: ip_auxilary  = 7
!
!* Essential parameters to dimension the simulation
  INTEGER (kind=ip_intwp_p) :: ig_nmodel, ig_nfield, ig_direct_nfield, ig_total_nfield
  LOGICAL :: lg_oasis_field
  INTEGER (kind=ip_intwp_p) :: ig_maxcomb
!* Parameters to dimension the main big arrays
  INTEGER (kind=ip_intwp_p) :: ig_nwork, ig_work
  INTEGER (kind=ip_intwp_p) :: ig_maxold, ig_maxnew
  INTEGER (kind=ip_intwp_p) :: ig_maxold_grid, ig_maxnew_grid
!* Parameters related to maxima of grid-dimension
  INTEGER (kind=ip_intwp_p) :: ig_maxgrd, ig_total_maxgrd
!* Parameters related to ANAIS(M-G), MOZAIC and SUBGRID interpolation
  INTEGER (kind=ip_intwp_p) :: ig_maxwoa, ig_maxnoa, ig_maxmoa, ig_maxsoa
  INTEGER (kind=ip_intwp_p) :: ig_maxnfm, ig_maxnfg, ig_maxnfp, ig_maxnfs, ig_maxnfn
!* Parameters related to extrapolation method WEIGHT
  INTEGER (kind=ip_intwp_p) :: ig_maxext, ig_maxnbn
!
!     -------------------------------------------------------------------
!
END MODULE MOD_PARAMETER
