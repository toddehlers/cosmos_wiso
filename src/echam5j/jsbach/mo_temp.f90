!+ <A one line description of this module> 
! 
Module mo_temp

  ! 
  ! Description: 
  !   <Say what this module is for> 
  ! 
  ! Current Code Owner: <Name of person responsible for this code> 
  ! 
  ! History: 
  !  
  ! Version   Date        Comment 
  ! -------   ----        ------- 
  ! 0.1       2001/10/02  Original code. Reiner Schnur
  ! 
  ! Code Description: 
  !   Language:           Fortran 90. 
  !   Software Standards: "European Standards for Writing and  
  !     Documenting Exchangeable Fortran 90 Code". 
  ! 
  ! Modules used: 
  ! 
  Use mo_kind,   Only : dp

       ! Imported Type Definitions: 
       
       ! Imported Parameters: 
       
       ! Imported Scalar Variables with intent (in): 
       
       ! Imported Scalar Variables with intent (out): 
       
       ! Imported Array Variables with intent (in): 
       
       ! Imported Array Variables with intent (out): 
       
       ! Imported Routines: 
       
       ! <Repeat from Use for each module...> 
       
       ! Declarations must be of the form: 
       ! <type>   <VariableName>      ! Description/ purpose of variable 
       
  IMPLICIT NONE
  
  ! Global (i.e. public) Declarations: 
  ! Global Type Definitions: 
  
  ! Global Parameters: 
  
  ! Global Scalars: 
  
  ! Global Arrays: 
  
  ! Local (i.e. private) Declarations: 
  INTEGER,  ALLOCATABLE, TARGET  :: zint1d(:)
  INTEGER,  ALLOCATABLE, TARGET  :: zzint1d(:)
  REAL(dp),     ALLOCATABLE, TARGET     :: zreal1d(:)
  REAL(dp), ALLOCATABLE, TARGET :: zdouble1d(:)
  LOGICAL,  ALLOCATABLE, TARGET  :: zlogical1d(:)
  INTEGER,  ALLOCATABLE, TARGET  :: zint2d(:,:)
  REAL(dp),     ALLOCATABLE, TARGET     :: zreal2d(:,:)
  REAL(dp),     ALLOCATABLE, TARGET     :: zzreal2d(:,:)
  REAL(dp), ALLOCATABLE, TARGET :: zdouble2d(:,:)
  LOGICAL,  ALLOCATABLE, TARGET  :: zlogical2d(:,:)
  INTEGER,  ALLOCATABLE, TARGET  :: zint3d(:,:,:)
  REAL(dp),     ALLOCATABLE, TARGET     :: zreal3d(:,:,:)
  REAL(dp),    ALLOCATABLE, TARGET     :: zzreal3d(:,:,:)
  REAL(dp), ALLOCATABLE, TARGET :: zdouble3d(:,:,:)
  INTEGER,  ALLOCATABLE, TARGET  :: zint4d(:,:,:,:)
  REAL(dp),    ALLOCATABLE, TARGET     :: zreal4d(:,:,:,:)
  REAL(dp), ALLOCATABLE, TARGET :: zdouble4d(:,:,:,:)

  REAL(dp), POINTER         :: zreal2d_ptr(:,:) => NULL()
  REAL(dp), POINTER         :: zreal3d_ptr(:,:,:) => NULL()
  REAL(dp), POINTER         :: zreal4d_ptr(:,:,:,:) => NULL()
  

End module mo_temp
