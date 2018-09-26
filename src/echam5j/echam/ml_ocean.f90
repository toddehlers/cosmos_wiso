  SUBROUTINE ml_ocean ( kproma                                         &
                  , pslm                                               &
                  , lonorth                                            &
                  , pseaice,    psiced,    palake                      &
                  , ptsi,       ptsw                                   &
                  , pahflw,     pahfsw,    pfluxres                    &
                  , ptrflw,     psoflw                                 &
                  , pamlcorr,   pamlcorac, pamlheatac                  &
                  , pevapi,     psni,      pcvsi                       &
                  , pahfres,    pfri                      )  

    !
    !  ---------------------------------------------------------------------
    !
    USE mo_kind,      ONLY: dp
    USE mo_exception, ONLY: finish
    !
    ! Arguments
    !
    INTEGER, INTENT(in) :: kproma
    REAL(dp) ::                                                              &
          pseaice(kproma),     psiced(kproma),     palake(kproma)            &
         ,ptsi(kproma),        ptsw(kproma)                                  &
         ,pahflw(kproma),      pahfsw(kproma),     pfluxres(kproma)          &
         ,ptrflw(kproma),      psoflw(kproma)                                &
         ,pamlcorr(kproma),    pamlcorac(kproma),  pamlheatac(kproma)        &
         ,pevapi(kproma),      psni(kproma),       pcvsi(kproma)             &
         ,pahfres(kproma),     pfri(kproma),       pslm(kproma)
  LOGICAL :: lonorth(kproma)      ! .true. for northern latitude

  !
  ! do not use this subroutine
  !
  CALL finish('ml_ocean','MLO: Option not yet available')
  !  ---------------------------------------------------------------------
END SUBROUTINE ml_ocean
