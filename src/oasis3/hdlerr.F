      subroutine hdlerr(istatus, cd_routine)
#ifdef use_netCDF
      USE mod_kinds_oasis
      USE mod_unit
      integer (kind=ip_intwp_p) :: istatus
      CHARACTER*(*) cd_routine
#include <netcdf.inc>
      IF (istatus .ne. NF_NOERR) THEN
          WRITE(nulou,*) NF_STRERROR(istatus)
          CALL prcout('STOP in ', cd_routine, 2)
          CALL halte('STOP')
      ENDIF
#endif
      return
      end
