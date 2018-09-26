!***************************************************************************
!                                                                          *
!                RRTM :  RAPID RADIATIVE TRANSFER MODEL                    *
!                                                                          *
!             ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                 *
!                        840 MEMORIAL DRIVE                                *
!                        CAMBRIDGE, MA 02139                               *
!                                                                          *
!                           ELI J. MLAWER                                  *
!                         STEVEN J. TAUBMAN~                               *
!                         SHEPARD A. CLOUGH                                *
!                        ~currently at GFDL                                *
!                       email:  mlawer@aer.com                             *
!                                                                          *
!        The authors wish to acknowledge the contributions of the          *
!        following people:  Patrick D. Brown, Michael J. Iacono,           *
!        Ronald E. Farren, Luke Chen, Robert Bergstrom.                    *
!                                                                          *
!***************************************************************************
!     Reformatted for F90 by JJMorcrette, ECMWF, 980714                    *
!***************************************************************************
! *** mji ***
! *** This version of RRTM has been altered to interface with either
!     the ECMWF numerical weather prediction model or the ECMWF column 
!     radiation model (ECRT) package. 

!     Revised, April, 1997;  Michael J. Iacono, AER, Inc.
!          - initial implementation of RRTM in ECRT code
!     Revised, June, 1999;  Michael J. Iacono and Eli J. Mlawer, AER, Inc.
!          - to implement generalized maximum/random cloud overlap
!     Bug-Fix, November 1999: Michael J. Iacono, AER, Inc.
!         and  December 1999: JJMorcrette  initializations within cloud 
!                             overlap
!     Modified for ECHAM5, April 2000, Marco A. Giorgetta, MPI
!     Update to ECMWF-Cy23R1, Dec2000, Marco A. Giorgetta, MPI
!     Shift longitude loop to rrtm subroutines, Jan2001, Marco A. Giorgetta, MPI

SUBROUTINE RRTM_RRTM_140GP        &
  !-- in
  ( kproma, kbdim, klev,                   &
    ppave, ptave, ptl,            &
    ptbound, psemiss,             &
    pcldfrac,kcldlyr,             &
    pcoldry, pwkl, pwx,           &
    ptaucld, ptauaerl,            &
  !-- out
    ptotuflux, ptotufluc,         &
    ptotdflux, ptotdfluc,         &
    psemit )

  ! *** This program is the driver for RRTM, the AER rapid model.  
  !     For each atmosphere the user wishes to analyze, this routine
  !     a) sets the atmospheric profile 
  !     b) calls SETCOEF to calculate various quantities needed for 
  !        the radiative transfer algorithm
  !     c) calls RTRN to do the radiative transfer calculation for
  !        clear or cloudy sky
  !     d) writes out the upward, downward, and net flux for each
  !        level and the heating rate for each layer


  !  Input
  USE mo_kind   ,ONLY : dp
  USE MO_PARRRTM,ONLY : JPBAND   ,JPXSEC  ,JPGPT    ,JPINPX

  !------------------------------Arguments--------------------------------

  IMPLICIT NONE

  ! Input arguments
  !
  INTEGER, INTENT(in)                             :: kproma    ! number of longitudes
  INTEGER, INTENT(in)                             :: kbdim     ! first dimension of 2-d arrays
  INTEGER, INTENT(in)                             :: klev      ! number of levels
  REAL(dp),INTENT(in),DIMENSION(kbdim,klev)       :: ppave     ! full level pressure [mb]
  REAL(dp),INTENT(in),DIMENSION(kbdim,klev)       :: ptave     ! full level temperature [K]
  REAL(dp),INTENT(in),DIMENSION(kbdim,klev+1)     :: ptl       ! half level temperature [K]
  REAL(dp),INTENT(in),DIMENSION(kbdim)            :: ptbound   ! surface temperature [K]
  REAL(dp),INTENT(in),DIMENSION(kbdim,jpband)     :: psemiss   ! surface emissivity in each band []
  INTEGER, INTENT(in),DIMENSION(kbdim,klev)       :: kcldlyr   ! cloud indicator, 0:clear, 1: cloudy
  REAL(dp),INTENT(in),DIMENSION(kbdim,klev)       :: pcldfrac  ! layer cloud fraction with respect to
  !                                                             cloud fraction of total column []
  REAL(dp),INTENT(in),DIMENSION(kbdim,klev)       :: pcoldry   ! number of molecules/cm2 of dry air and
  !                                                             water vapor [#/cm2]
  REAL(dp),INTENT(in),DIMENSION(kbdim,jpinpx,klev):: pwkl      ! number of molecules/cm2 of N species
  !                                                             in [#/cm2], N=JPINPX
  !                                                             1: H2O
  !                                                             2: CO2
  !                                                             3: O3
  !                                                             4: N2O
  !                                                             5: ------ empty ------
  !                                                             6: CH4
  !                                                             7... : -- empty ------
  REAL(dp),INTENT(in),DIMENSION(kbdim,jpxsec,klev):: pwx       ! number of molecules/cm2 of N species
  !                                                             in [1e20/cm2], N=JPXSEC
  !                                                             1: ------ empty ------
  !                                                             2: CFC11
  !                                                             3: CFC12
  !                                                             4... : -- empty ------
  REAL(dp),INTENT(in),DIMENSION(kbdim,klev,jpband):: ptaucld   ! optical thickness of clouds
  !                                                             in each band []
  REAL(dp),INTENT(in),DIMENSION(kbdim,klev,jpband):: ptauaerl  ! optical thickness of aerosols
  !                                                             in each band []
  !
  ! Output arguments
  !
  REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1)    :: ptotuflux ! upward flux, total sky
  REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1)    :: ptotufluc ! upward flux, clear sky
  REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1)    :: ptotdflux ! downward flux, total sky
  REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1)    :: ptotdfluc ! downward flux, clear sky
  REAL(dp),INTENT(out),DIMENSION(kbdim)           :: psemit    ! surface emissivity

  !------------------------------RRTM variables---------------------------
  !
  ! These local variables have been defines in older versions in modules mo_rrtXYZ.
  ! Modules mo_rrtXYZ have been composed from fortran77 common blocks.

  !- from mo_rrtatm, COMMON INTFAC
  REAL(dp),DIMENSION(kbdim,klev)        :: FAC00
  REAL(dp),DIMENSION(kbdim,klev)        :: FAC01
  REAL(dp),DIMENSION(kbdim,klev)        :: FAC10
  REAL(dp),DIMENSION(kbdim,klev)        :: FAC11
  REAL(dp),DIMENSION(kbdim,klev)        :: FORFAC
  !
  !- from mo_rrtatm, COMMON INTIND
  INTEGER, DIMENSION(kbdim,klev)        :: JP
  INTEGER, DIMENSION(kbdim,klev)        :: JT
  INTEGER, DIMENSION(kbdim,klev)        :: JT1
  !
  !- from mo_rrtatm, COMMON PRECISE
  ! ONEMINUS is a real number just smaller than 1
  ! used in rrtm_taugbN some subroutines
  REAL(dp),PARAMETER                   :: ONEMINUS =1._dp-1.E-06_dp

  !
  !- from mo_rrtatm, COMMON PROFDATA             
  REAL(dp),DIMENSION(kbdim,klev)        :: COLH2O
  REAL(dp),DIMENSION(kbdim,klev)        :: COLCO2
  REAL(dp),DIMENSION(kbdim,klev)        :: COLO3
  REAL(dp),DIMENSION(kbdim,klev)        :: COLN2O
  REAL(dp),DIMENSION(kbdim,klev)        :: COLCH4
  REAL(dp),DIMENSION(kbdim,klev)        :: CO2MULT
  INTEGER, DIMENSION(kbdim)             :: LAYTROP
  INTEGER, DIMENSION(kbdim)             :: LAYSWTCH
  INTEGER, DIMENSION(kbdim)             :: LAYLOW
  !
  !- from mo_rrtatm, COMMON SELF             
  REAL(dp),DIMENSION(kbdim,klev)        :: SELFFAC
  REAL(dp),DIMENSION(kbdim,klev)        :: SELFFRAC
  INTEGER, DIMENSION(kbdim,klev)        :: INDSELF
  !
  !- from mo_rrtatm, COMMON SP             
  REAL(dp),DIMENSION(kbdim,jpgpt,klev)  :: PFRAC
  !
  ! former module mo_rrtctl
  ! -----------------------
  INTEGER, PARAMETER                   :: ISTART  = 1
  INTEGER, PARAMETER                   :: IEND    = JPBAND

  ! former module mo_rrtgd
  ! ----------------------
  ! - except for tau which is defined in rrtm_gasabs1a_140gp
  ! - equivalence is avoided
  REAL(dp), DIMENSION(kbdim,jpgpt*klev) :: ABSS1
  REAL(dp), DIMENSION(kbdim,jpgpt,klev) :: OD
  REAL(dp), DIMENSION(kbdim,jpgpt*klev) :: TAUSF1
  

  !  Calculate information needed by the radiative transfer routine
  !  that is specific to this atmosphere, especially some of the 
  !  coefficients and indices needed to compute the optical depths
  !  by interpolating data from stored reference atmospheres. 
  
  CALL RRTM_SETCOEF_140GP (KPROMA,KBDIM,KLEV,pcoldry,pwkl &
       &, FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1 &
       &, COLH2O,COLCO2,COLO3,COLN2O,COLCH4,CO2MULT &
       &, LAYTROP,LAYSWTCH,LAYLOW,ppave,ptave,SELFFAC,SELFFRAC,INDSELF)

  CALL RRTM_GASABS1A_140GP (KPROMA,KBDIM,KLEV,ABSS1,OD,TAUSF1,pcoldry,pwx &
       &, ptauaerl,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,ONEMINUS &
       &, COLH2O,COLCO2,COLO3,COLN2O,COLCH4,CO2MULT &
       &, LAYTROP,LAYSWTCH,LAYLOW,SELFFAC,SELFFRAC,INDSELF,PFRAC)

  !- Call the radiative transfer routine.

  !  Clear and cloudy parts of column are treated together in RTRN.
  !  Clear radiative transfer is done for clear layers and cloudy radiative
  !  transfer is done for cloudy layers as identified by icldlyr.

  CALL RRTM_RTRN1A_140GP (KPROMA,KBDIM,KLEV,ISTART,IEND,kcldlyr,pcldfrac,ptaucld,ABSS1 &
       &, OD,TAUSF1,ptotdfluc,ptotdflux,ptotufluc,ptotuflux &
       &, ptave,ptl,ptbound,PFRAC,psemiss,psemit)

END SUBROUTINE RRTM_RRTM_140GP
