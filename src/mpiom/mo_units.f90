MODULE mo_units
  ! ---------------------------------------------------------------------
  !
  !*    *common* *units*      - logical units for i/o of the ocean model.
  !
  !                           - included in all subroutines.
  !                               
  !
  !     s. legutke          *dkrz*           07.01.00
  !
  !*    variable     type        purpose.
  !     --------     ----        --------
  !     *io_stderr*  *integer*   logical unit for ocean std eror
  !     *io_stdout*  *integer*   logical unit for ocean stdout
  !     *io_in_octl* *integer*   logical unit for ocean namelist (ocectl)
  !     ** *integer*   
  !
  ! ---------------------------------------------------------------------
  !

  INTEGER :: io_stdout , io_stderr , io_in_octl,                          &
       io_in_arcg, io_in_gigi, io_in_anta,                          &
       io_in_z370, io_in_z380, io_in_surt,                          &
       io_in_rval, io_in_rpos, io_in_surs, io_in_glac,              &
       io_in_gtem, io_in_gwix, io_in_gwiy, io_in_gclo,io_in_glwr,   &
       io_in_gpre, io_in_gswr, io_in_gtde, io_in_gu10,io_in_gslp,   &
       io_in_inis, io_in_init, io_in_inii,                          &
       io_in_bgin 

  INTEGER :: io_ou_f125, io_ou_f126, io_ou_scha, io_ou_somi, io_ou_ray,   &
       io_ou_mtho, io_ou_msao,              &
       io_ou_emip, io_ou_mmzo, io_ou_flum, io_ou_mpem,              &
       io_ou_sict, io_ou_sico, io_ou_sicu, io_ou_sicv,              &
       io_ou_sics, io_ou_tafo, io_ou_fclo, io_ou_fpre,              &
       io_ou_fswr, io_ou_ftde, io_ou_qswo, io_ou_qlwo,              &
       io_ou_qseo, io_ou_qlao, io_ou_prec, io_ou_amld, io_ou_psiu,  &
       io_ou_weto, io_ou_gila, io_ou_giph, io_ou_dept, io_ou_zmld,  &
       io_ou_sictru, io_ou_sictrv,io_ou_txo,io_ou_tye,io_ou_mwgo,   &
       io_ou_mavo,io_ou_mdvo,io_ou_mwo,io_ou_dlxp,io_ou_dlyp,       & 
       io_ou_wu10,io_ou_deut,io_ou_dlxu,io_ou_dlyu,io_ou_wtmi,io_ou_rinu,&
       io_ou_dqswo,io_ou_dqlwo,io_ou_dqseo,io_ou_dqlao,io_ou_dqtho, &
       io_ou_dqswi,io_ou_dqlwi,io_ou_dqsei,io_ou_dqlai,io_ou_dqthi, &
       io_ou_dticeo,io_ou_mbolx,io_ou_mboly,io_ou_f090,io_ou_tmceo, &
       io_ou_tmcdo,io_ou_ukomfl,io_ou_vkemfl,io_ou_rivrun,          &
       io_ou_aonhw,io_ou_aoshw,io_ou_aorhi,io_ou_aochi,io_ou_aofrw, &
       io_ou_aofri,io_ou_aotxw,io_ou_aotyw,io_ou_aotxi,io_ou_aotyi, &
       io_ou_aowsv,io_ou_mocg,io_ou_moca,io_ou_difi,io_ou_tsvar,    &
       io_ou_tsdesc,io_ou_tsunit

#ifdef ADDCONTRA
  INTEGER :: io_ou_aofrwo16, io_ou_aofrio16,                        &
             io_ou_aofrwo18, io_ou_aofrio18,                        &
             io_ou_aofrwhdo, io_ou_aofrihdo
#endif

  INTEGER :: io_in_coru10,io_in_corv10,io_in_cort10,io_in_corq10,         &
       io_in_corlw,io_in_corsw,io_in_corpr,io_in_corro

CONTAINS

  SUBROUTINE setunits

!  USE mo_commo1, only :lgmdiag,ldiffdiag                           &
!                        ,lhfldiag,lgridinfo,lconvdiag,imean


#ifdef __coupled
    io_stdout=8
    io_stderr=0
#else
    io_stdout=6
#ifdef __stdout0
    io_stdout=0
#endif

#endif
    io_in_octl=10

    !     modelgrid and topofile
    io_in_arcg=21
    io_in_gigi=91
    io_in_anta=81

    !     restart files
    io_in_z370=27
    io_in_z380=28

    !     forcing
    io_in_rval=20
    io_in_rpos=30
    io_in_surs=38
    io_in_surt=338
    io_in_glac=39
    io_in_gtem=37
    io_in_gwix=40
    io_in_gwiy=41
    io_in_gclo=51
    io_in_gpre=52
    io_in_gswr=53
    io_in_gtde=54
    io_in_gu10=55

#ifdef FORCE_DWLW
    io_in_glwr=56
#endif /*FORCE_DWLW*/
#ifdef FORCE_SLP
    io_in_gslp=57
#endif /*FORCE_SLP*/

#ifdef CORE
    io_in_coru10=51
    io_in_corv10=52
    io_in_cort10=53
    io_in_corq10=54
    io_in_corlw=55
    io_in_corsw=56
    io_in_corpr=57
    io_in_corro=58
#endif

    !     newstart/3drestoring
    io_in_inis=33
    io_in_init=32

    !     diagnostic
    io_in_bgin=115

    !     set the output units

    io_ou_f125=125
    io_ou_f126=126

#ifdef KONVDIAG
    io_ou_f090=73
#endif
    io_ou_amld=73

    !     amoc
    io_ou_scha=47
    io_ou_somi=83
    io_ou_ray=92

    io_ou_mocg=75
    io_ou_moca=75


!    if (IMEAN.ne.0)then
       io_ou_mtho=71
       io_ou_msao=71
       io_ou_ukomfl=71
       io_ou_vkemfl=71
       io_ou_emip=71
       io_ou_mmzo=71
       io_ou_flum=71
       io_ou_mpem=71
       io_ou_sict=71
       io_ou_sico=71
       io_ou_sicu=71
       io_ou_sicv=71
       io_ou_sics=71
       io_ou_tafo=71
       io_ou_wu10=71
       io_ou_fclo=71
       io_ou_fpre=71
       io_ou_fswr=71
       io_ou_ftde=71
       io_ou_qswo=71
       io_ou_qlwo=71
       io_ou_qseo=71
       io_ou_qlao=71
       io_ou_prec=71
       io_ou_zmld=71

       io_ou_psiu=71
       io_ou_sictru=71
       io_ou_sictrv=71
       io_ou_txo=71
       io_ou_mwo=71
       io_ou_tye=71
       io_ou_rivrun=71
       
!       if (LGMDIAG) then
          io_ou_mwgo=71
          io_ou_mbolx=71
          io_ou_mboly=71
!       endif

!       if (LDIFFDIAG) then
          io_ou_mavo=71
          io_ou_mdvo=71
          io_ou_wtmi=71
          io_ou_rinu=71
!       endif

!       if (LHFLDIAG) then
          io_ou_dqswo=71
          io_ou_dqlwo=71
          io_ou_dqseo=71
          io_ou_dqlao=71
          io_ou_dqtho=71
          io_ou_dqswi=71
          io_ou_dqlwi=71
          io_ou_dqsei=71
          io_ou_dqlai=71
          io_ou_dqthi=71
          io_ou_dticeo=71
!       endif
    
!    endif
 

#ifdef __coupled
    io_ou_aonhw=270
    io_ou_aoshw=270
    io_ou_aorhi=270
    io_ou_aochi=270
    io_ou_aofrw=270
    io_ou_aofri=270
    io_ou_aotxw=270
    io_ou_aotyw=270
    io_ou_aotxi=270
    io_ou_aotyi=270
    io_ou_aowsv=270
#ifdef ADDCONTRA
    io_ou_aofrwo16=270
    io_ou_aofrio16=270
    io_ou_aofrwo18=270
    io_ou_aofrio18=270
    io_ou_aofrwhdo=270
    io_ou_aofrihdo=270
#endif
#endif /*__coupled*/


    
!    IF (LGRIDINFO) then
       io_ou_weto=72
       io_ou_gila=73
       io_ou_giph=73
       io_ou_dept=72
       io_ou_dlxp=72
       io_ou_dlyp=72
       io_ou_deut=72
       io_ou_dlxu=72
       io_ou_dlyu=72
!    endif

!    if (LCONVDIAG) then
       io_ou_tmceo=71
       io_ou_tmcdo=71
!    endif

      io_ou_difi=74

      io_ou_tsvar=100
      io_ou_tsdesc=101
      io_ou_tsunit=102

  END SUBROUTINE setunits

END MODULE mo_units




