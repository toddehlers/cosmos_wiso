MODULE MO_PARAM1

  USE mo_kind, ONLY: i8

  IMPLICIT NONE

  ! Global model dimensions
   
  INTEGER :: IE_G, JE_G, KE 

  ! Local dimensions (set by domain decompostion)
  
  INTEGER IE, JE
  
  ! Derived Parameters (set in set_param1 below)
  
  INTEGER ITO, JTO, IE1, IE2, IT3, KEP, KE1, JE1, IT4, JE2, JT2,    &
       IT1, JEH, IEJTO, IEJEKE, IEJE
  
  ! SOLVER
  
#ifdef SOR
  INTEGER, PARAMETER :: IELIMI=0,ICYCLI=1
#else
  INTEGER, PARAMETER :: IELIMI=1,ICYCLI=1
#endif
  
  INTEGER KBB, ILL, IMM, ILT, ILLT, KBBT
  
  INTEGER, PARAMETER :: NBOX=9
  
  !HH   NUMBER OF RIVERS
#ifdef OMIP_RIV
  INTEGER, PARAMETER :: NUMRIV=6366
#else
  INTEGER, PARAMETER :: NUMRIV=52
#endif
#ifdef GLACCALV
  INTEGER, PARAMETER :: NUMGLAC=5653
#endif
  INTEGER(i8) :: ibla, ii1, ii2, ii3, ii4, ii5
  
CONTAINS

  SUBROUTINE set_param1

    ! Sets all derived parameters dpeneding on local settings for IE and JE
    ! Since the times of FORTRAN 66 are long ago, these parameters
    ! shouldn't be used any more

!   VERSIONT43 
    if (ie_g.eq.130.and.je_g.eq.211) then
	KBBT=340
	ILLT=15900
    endif

!   VERSIONGR30
    if (ie_g.eq.122.and.je_g.eq.101) then
	KBBT=218
	ILLT=8430
    endif

! VERSIONGR60
    if (ie_g.eq.60.and.je_g.eq.50) then
	KBBT=218
	ILLT=8430
    endif

! VERSIONTP40
    if (ie_g.eq.82.and.je_g.eq.41) then
	KBBT=340
	ILLT=15900
    endif

! VERSIONTP04
    if (ie_g.eq.802.and.je_g.eq.401) then
	KBBT=800
	ILLT=200000
    endif

! VERSIONGR15
    if (ie_g.eq.256.and.je_g.eq.220) then
	KBBT=400
	ILLT=36504
    endif

! VERSIONGR09
    if (ie_g.eq.400.and.je_g.eq.338) then
	KBBT=700
	ILLT=95000
    endif

! VERSIONGR03
    if (ie_g.eq.1080.and.je_g.eq.680) then
	KBBT=2000
	ILLT=295000
    endif

! VERSIONGIN
    if (ie_g.eq.182.and.je_g.eq.84) then
	KBBT=118
	ILLT=7430
    endif



    ito    = 2*ie
    jto    = 2*je
    ie1    = ie-1
    ie2    = ie-2
    it3    = ito+3
    kep    = ke+1
    ke1    = ke-1
    je1    = je-1
    it4    = ito-4
    je2    = je-2
    jt2    = jto-2
    it1    = ito-1
    jeh    = je/2
    iejto  = ie*jto
    iejeke = ie*je*ke
    ieje   = ie*je
    
#ifdef SOR
    kbb=4
    ill=9
#else
    kbb=kbbt
    ill=illt
#endif
    imm=2*kbb+1
    ilt=ill+kbb
    
  END SUBROUTINE set_param1

END MODULE mo_param1
