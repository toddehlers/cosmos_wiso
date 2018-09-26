SUBROUTINE sudif

  ! Description:
  !
  ! Initializes module *mo_diff* for horizontal diffusion subroutines
  !
  ! Authors:
  !
  ! M. Esch, MPI, September 1993, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! M. Esch, MPI, August 1999, modifications for ECHAM5
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_control,     ONLY: nlev, nn, lmidatm
  USE mo_diff,        ONLY: ncdif, iq
  USE mo_exception,   ONLY: finish, message

  IMPLICIT NONE

  !  Local scalars: 
  INTEGER :: jn

  !  Local arrays: 
  INTEGER :: ihq319(nlev), ihq511(nlev)
  INTEGER :: ihq213(nlev), ihq255(nlev)
  INTEGER :: ihq159(nlev)
  INTEGER :: ihq106(nlev), ihq85(nlev)
  INTEGER :: ihq63(nlev),  ihq42(nlev)
  INTEGER :: ihq31(nlev),  ihq21(nlev)

  !  Executable statements 

!-- 1. Set parameter
  
  IF (lmidatm) THEN
    ihq511(:)  = 2
    ihq319(:)  = 2
    ihq255(:)  = 2
    ihq159(:)  = 3
    ihq106(:)  = 4
    ihq85(:)   = 4
    ihq63(:)   = 4
    ihq42(:)   = 5
    ihq31(:)   = 5
    ihq21(:)   = 2
  ELSE
    ihq511(1:3) = 1
    ihq511(4:)  = 2

    ihq319(1:3) = 1
    ihq319(4:)  = 2

    ihq255(1:3) = 1
    ihq255(4:)  = 2

    ihq213(1:3) = 1
    ihq213(4)   = 2
    ihq213(5:)  = 3

    ihq159(1:3) = 1
    ihq159(4)   = 2
    ihq159(5:)  = 3

    ihq106(1:3) = 1
    ihq106(4)   = 2
    ihq106(5)   = 3
    ihq106(6:)  = 4

    ihq85(1:3)  = 1
    ihq85(4)    = 2
    ihq85(5)    = 3
    ihq85(6:)   = 4

    ihq63(1:3)  = 1
    ihq63(4)    = 2
    ihq63(5)    = 3
    ihq63(6:)   = 4

    ihq42(1:3)  = 1
    ihq42(4)    = 2
    ihq42(5)    = 3
    ihq42(6)    = 4
    ihq42(7:)   = 5

    ihq31(1:3)  = 1
    ihq31(4)    = 2
    ihq31(5)    = 3
    ihq31(6)    = 4
    ihq31(7:)   = 5

    ihq21(1)    = 1
    ihq21(2)    = 2
    ihq21(3)    = 3
    ihq21(4)    = 4
    ihq21(5)    = 6
    ihq21(6)    = 8
    ihq21(7:)   = 10
  ENDIF

  ncdif(1:nlev) = 0

!-- 2. Copy to iq

  DO jn = 1, nlev
    IF (nn==511) THEN
      iq(jn) = ihq511(jn)      
    ELSE IF (nn==319) THEN
      iq(jn) = ihq319(jn)      
    ELSE IF (nn==255) THEN
      iq(jn) = ihq255(jn)      
    ELSE IF (nn==213) THEN
      iq(jn) = ihq213(jn)
    ELSE IF (nn==159) THEN
      iq(jn) = ihq159(jn)      
    ELSE IF (nn==106) THEN
      iq(jn) = ihq106(jn)
    ELSE IF (nn==85) THEN
      iq(jn) = ihq85(jn)
    ELSE IF (nn==63) THEN
      iq(jn) = ihq63(jn)
    ELSE IF (nn==42) THEN
      iq(jn) = ihq42(jn)
    ELSE IF (nn==31) THEN
      iq(jn) = ihq31(jn)
    ELSE IF (nn==21) THEN
      iq(jn) = ihq21(jn)
    ELSE
      CALL message('sudif',' This model resolution is not supported')
      CALL finish('sudif','Run terminated.')
    END IF
  END DO

  RETURN
END SUBROUTINE sudif
