      SUBROUTINE CPPOUTPUT
      USE MO_PARAM1
      USE MO_COMMO1
      USE MO_UNITS
!:: SUMMARIZE COMPILER OPTION AND PARAMTER SETTING FOR OUTPUT
!
       write(io_stdout,*)' '
       write(io_stdout,*)'LIST OF COMPILER OPTIONS AND PARAMETERS'
       write(io_stdout,*)' '
!#ifdef MEAN   
!       write(io_stdout,*)'WRITING MEAN FIELDS'
!#endif
!
#ifdef EISREST
       write(io_stdout,*)'RESTORING UNDER ICE!'
#endif
!
#ifdef AULREDSC
       write(io_stdout,*)'DIFFUSION PROPORTIONAL TO DX,DY**3!'
#else
       write(io_stdout,*)'DIFFUSION PROPORTIONAL TO DX,DY**4!'
#endif
!
#ifdef REDWMICE
       write(io_stdout,*)'REDUCED WIND MIXING UNDER ICE'
#endif 
!
       if (iocad.eq.4)then
          write(io_stdout,*)'SLOPE CONVECTION THROUGH ADVECTION in ADPO!'
       endif
!
       IF (I3DREST .GT. 0) write(io_stdout,*)'3D-RESTORING!!!!'
!
!OtB_OLDMOMADV ifdef OLDMOMADV 
!OtB_OLDMOMADV       write(io_stdout,*)'USING OLD MOMENTUM ADVECTION!'
!OtB_OLDMOMADV endif
!
 if (iocad.eq.3.or.iocad.eq.5)then
       write(io_stdout,*)'USING ADPO TRACER ADVECTION!'
 endif

 if (iocad.eq.5)then
       write(io_stdout,*)'USING ADFS TRACER ADVECTION!'
 endif

#ifdef QUICK2    
       write(io_stdout,*)'USING QUICK2 SCHEME!'
#endif
!
#ifdef FREESLIP
       write(io_stdout,*)'BIHARMONIC FRICTION WITH FREE SLIP'
#endif
!
#ifdef ISOPYK  
       write(io_stdout,*)'USING ISOPYCNIC DIFFUSION         '
#endif
!
       IF (ibolk .NE. 0) THEN
          WRITE(io_stdout,*)'USING GENT MCWILLIAMS EDDY PARAMETERIZATION'
          IF (ibolk .LT. 0) THEN
             WRITE(io_stdout,*)'USING VISBECK ET AL. EDDY COEFFICIENT'
          ENDIF
       ENDIF
!
#ifdef HARM     
       write(io_stdout,*)'HARMONIC MOMENTUM DIFFUSUIN'
       write(io_stdout,*)'aus= ',aus
#endif
!
      WRITE(io_stdout,*)' AH00= ',ah00
      WRITE(io_stdout,*)' AH= ',ah
!
#ifdef DBACKGFDL
      WRITE(io_stdout,*)' OPTION BACKGROUND VERT. DIFF AFTER GFDL'
#endif
!
#ifdef DBACK3E5 
      WRITE(io_stdout,*)' OPTION BACKGROUND VERT. DIFF 3E5'
#endif
!
#ifdef DBACK0
      WRITE(io_stdout,*)' OPTION BACKGROUND VERT. DIFF 0.0'
#endif
!
#ifdef PLUME
      WRITE(io_stdout,*)' OPTION PLUME CONVECTION'
#endif
       write(io_stdout,*)' '
      RETURN
      END
