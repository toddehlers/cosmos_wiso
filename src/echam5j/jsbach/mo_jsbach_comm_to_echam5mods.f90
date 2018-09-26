!+ <A one line description of this module> 
! 
MODULE mo_jsbach_comm_to_echam5mods
  ! 
  ! Description: 
  !   <Say what this module is for> 
  ! 
  ! Current Code Owner: <Name of person responsible for this code> 
  ! 
  ! History: 
  !  
  ! Version   Date     Comment 
  ! -------   ----     ------- 
  ! <version> <date>   Original code. <Your name> 
  ! 
  ! Code Description: 
  !   Language:           Fortran 90. 
  !   Software Standards: "European Standards for Writing and  
  !     Documenting Exchangeable Fortran 90 Code". 
  ! 

  USE mo_kind, ONLY: dp

  IMPLICIT NONE 

  INTEGER           :: nlon, nlat, nland
  INTEGER,  POINTER :: kpoints(:)
  REAL(dp), POINTER :: lon(:), lat(:)
  LOGICAL,  POINTER :: mask(:,:)

  INTEGER           :: domain_nlon, domain_nlat, domain_nland
  LOGICAL,  POINTER :: domain_mask(:,:)

END MODULE mo_jsbach_comm_to_echam5mods

!- End of module header
