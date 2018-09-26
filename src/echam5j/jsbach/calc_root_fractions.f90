!+ Computes the fraction of roots in each soil layer
! 
PURE FUNCTION calc_root_fractions ( & 
     root_depth, root_fract, soil_depth ) RESULT(veg_root)

  ! Description: 
  !   This routine computes the fraction of roots in each soil layer based on the
  !   root zone distribution (depth and fraction). Roots are assumed to be
  !   linearly distributed within each root zone.
  ! 
  !   This is a pure function so it can be called in a FORALL statement. It
  !   requires to an explicit interface in the calling procedure.
  !
  ! Method: 
  !   <Say how it does it: include references to external documentation> 
  ! 
  ! Current Code Owner: <Name of person responsible for this code> 
  ! 
  ! History: 
  ! Version   Date        Comment 
  ! -------   ----        ------- 
  ! 0.1       2001/10/14  Original code. Reiner Schnur
  !                       Based on VIC version 4.0.3 (C code)
  ! 
  ! Code Description: 
  !   Language:           Fortran 90. 
  !   Software Standards: "European Standards for Writing and  
  !     Documenting Exchangeable Fortran 90 Code". 
  ! 
  ! Modules used: 

  USE mo_kind, ONLY: dp

  IMPLICIT NONE 
  
  ! Function arguments 
  REAL(dp), INTENT(in) :: root_depth(:),  & !! Depth and fraction of roots for each
                      root_fract(:)     !! root zone
  REAL(dp), INTENT(in) :: soil_depth(:)     !! Depth of each soil layer

  ! Function result
  REAL(dp) :: veg_root(SIZE(soil_depth))    !! Fraction of roots in each soil layer
  
  ! Local scalars (must not be initialized in declaration since this is a 
  ! pure function!)
  INTEGER :: nRootZones          !! Number of root zones
  INTEGER :: nSoilLayers         !! Number of soil layers

  INTEGER :: iroot
  INTEGER :: ilayer
  REAL(dp)    :: zsum
  REAL(dp)    :: zstep
  REAL(dp)    :: lsum
  REAL(dp)    :: lstep
  REAL(dp)    :: sum_depth
  REAL(dp)    :: sum_fract
  REAL(dp)    :: zmin_depth
  REAL(dp)    :: zmin_fract
  REAL(dp)    :: zmax
  REAL(dp)    :: dum
  INTEGER :: i


  !- End of header --------------------------------------------------------------- 

  nRootZones  = SIZE(root_depth)
  nSoilLayers = SIZE(root_depth)

  iroot  = 1
  ilayer = 1
  zsum   = 0.0_dp
  lstep  = soil_depth(1)
  lsum   = lstep
  sum_depth = 0.0_dp
  sum_fract = 0.0_dp

  DO WHILE(iroot <= nRootZones)

     zstep = root_depth(iroot)
     IF (zsum+zstep <= lsum .AND. zsum >= lsum-lstep) THEN
        ! Case 1: Root zone completely in soil layer
        sum_fract = sum_fract + root_fract(iroot)
     ELSE
        ! Case 2: Root zone partially soil layer
        IF (zsum < lsum-lstep) THEN
           ! Root zone starts in previous soil layer             
           zmin_depth = lsum-lstep
           zmin_fract = root_fract(iroot) * &
                (zmin_depth - zsum) / (zsum+zstep - zsum)
        ELSE
           ! Root zone starts in current soil layer
           zmin_depth = zsum
           zmin_fract = 0.0_dp
        END IF
        IF (zsum+zstep <= lsum) THEN
           ! Root zone ends in current soil layer
           zmax = zsum + zstep
        ELSE
           ! Root zone extends beyond bottom of current soil layer
           zmax = lsum
        END IF
        sum_fract = sum_fract + (root_fract(iroot)-zmin_fract) * &
             (zmax - zsum) / (zsum+zstep - zsum)
     END IF
      
     ! Update current root zone and soil layer
     IF (zsum+zstep < lsum) THEN
        zsum = zsum + zstep
        iroot = iroot + 1
     ELSE IF (zsum+zstep == lsum) THEN
        zsum = zsum + zstep
        iroot = iroot + 1
        IF (ilayer <= nSoilLayers) THEN
           veg_root(ilayer) = sum_fract
           sum_fract = 0.0_dp
        ELSE IF (ilayer == nSoilLayers+1) THEN
           lstep = zsum + zstep - lsum
           IF (iroot < nRootZones) THEN
              DO i=iroot+1,nRootZones
                 lstep = lstep + root_depth(i)
              END DO
           END IF
           lsum = lsum + lstep
        END IF
     ELSE IF (zsum+zstep > lsum) THEN
        IF (ilayer <= nSoilLayers) THEN
           veg_root(ilayer) = sum_fract
           sum_fract = 0.0_dp
        END IF
        ilayer = ilayer + 1
        IF (ilayer <= nSoilLayers) THEN
           lstep = soil_depth(ilayer)
           lsum = lsum + lstep
        ELSE IF (ilayer == nSoilLayers+1) THEN
           lstep = zsum + zstep - lsum
           IF (iroot < nRootZones) THEN
              DO i=iroot+1,nRootZones
                 lstep = lstep + root_depth(i)
              END DO
           END IF
           lsum = lsum + lstep
        END IF
      END IF
      
   END DO
   
   IF (sum_fract > 0 .AND. ilayer > nSoilLayers) THEN
      veg_root(nSoilLayers) = veg_root(nSoilLayers) + sum_fract
   ELSE IF (sum_fract > 0) THEN
      veg_root(ilayer) = veg_root(ilayer) + sum_fract
   END IF

   dum = 0.0_dp
   DO i=1,nSoilLayers
      IF (veg_root(i) < 1.E-4_dp) veg_root(i) = 0.0_dp
      dum = dum + veg_root(i)
   END DO

   IF (dum == 0.0_dp) THEN
      veg_root = 0.0_dp        !! Error .. test for ALL(veg_root==0.0) in caller
   ELSE
      DO i=1,nSoilLayers
         veg_root(i) = veg_root(i) / dum
      END DO
   END IF

END FUNCTION calc_root_fractions
