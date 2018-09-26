SUBROUTINE inidoc

  ! Description:
  !
  ! Preset constants in mo_doctor.
  !
  ! Method:
  !
  ! *inidoc* is called from *initialize*.
  !
  ! Authors:
  !
  ! J. K. Gibson, ECMWF, April 1983, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,   ONLY: dp
  USE mo_doctor, ONLY: nvers                                           &
                     , ylabel1, ylabel2, ylabel3, ylabel4              &
                     , ylabel5, ylabel6, ylabel7, ylabel8
  USE mo_mpi,    ONLY: p_pe, p_io, p_bcast
  USE mo_netcdf, ONLY: global_att, put_att

  IMPLICIT NONE

  !  Local scalars: 
  INTEGER :: nlena, nlenb, nlenc
  REAL(dp):: zvers
  CHARACTER (3)   :: yovers
  CHARACTER (11)   :: yoversion  !---wiso-code change
  CHARACTER (256) :: os_name, user_name, host_name
  CHARACTER (8)   :: ydate
  CHARACTER (10)  :: ytime

  !  External subroutines
  EXTERNAL util_os_system, util_user_name, util_node_name

  !  Executable statements 

  yovers = '5.4'
  yoversion = '5.4.01-wiso'    !---wiso-code change
  READ (yovers,'(f8.1)') zvers
  nvers = NINT(zvers*10.0_dp) + 10
  IF (p_pe == p_io) THEN
     CALL util_os_system (os_name, nlena)
     CALL util_user_name (user_name, nlenb)
     CALL util_node_name (host_name, nlenc)
  END IF
  CALL p_bcast (os_name, p_io)
  CALL p_Bcast (nlena, p_io)
  CALL p_bcast (user_name, p_io)
  CALL p_Bcast (nlenb, p_io)
  CALL p_bcast (host_name, p_io)
  CALL p_Bcast (nlenc, p_io)

  ylabel1(:) = ' '
  ylabel2(:) = ' '
  ylabel3(:) = ' '
  ylabel4(:) = ' '
  ylabel5(:) = ' '
  ylabel6(:) = ' '
  ylabel7(:) = ' '
  ylabel8(:) = ' '

  CALL DATE_AND_TIME(ydate, ytime)

  ylabel1(1:) = ' Atmospheric model version ' // yoversion
  ylabel2(1:) = ' Library 17-Oct-2007'
  ylabel3(1:) = ' Lin & Rood ADVECTION is default'
  ylabel4(1:) = ' Modified ECMWF physics'
  ylabel5(1:) = ' Modified ECMWF radiation'
  ylabel6(1:) = ' Date - ' // ydate(1:8) // ' Time - ' // ytime(1:6)
  ylabel7(2:) = user_name(1:nlenb) // ' on ' // host_name(1:nlenc) 
  ylabel8(2:) = os_name(1:nlena)
  !---------------------------------------------------
  ! set global attributes to be written to NetCDF file
  !---------------------------------------------------
  CALL put_att (global_att,'source','ECHAM'//yovers//&
                           ' Max-Planck-Institute for Meteorology, Hamburg')
  CALL put_att (global_att,'echam_version',yoversion)
  CALL put_att (global_att,'institution',&
                           'Max-Planck-Institute for Meteorology, Hamburg')
  CALL put_att (global_att,'advection','Lin & Rood')
  call put_att (global_att,'physics','Modified ECMWF physics')
  call put_att (global_att,'radiation','Modified ECMWF radiation')
  call put_att (global_att,'date_time',ydate(1:8)//' '//ytime(1:6))
  call put_att (global_att,'operating_system',os_name(1:nlena))
  call put_att (global_att,'user_name',user_name(1:nlenb))
  call put_att (global_att,'host_name',host_name(1:nlenc))
  RETURN
END SUBROUTINE inidoc
