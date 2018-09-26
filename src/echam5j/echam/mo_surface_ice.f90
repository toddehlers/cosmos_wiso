MODULE mo_surface_ice

  USE mo_convect_tables,   ONLY: tlucua, jptlucu1, jptlucu2, &
       lookuperror, lookupoverflow
  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  REAL(dp), SAVE :: calbmni     ! minimum (bare sea ice)
  REAL(dp), SAVE :: calbmxi     ! maximum (bare sea ice)
  REAL(dp), SAVE :: calbmns     ! minimum (snow on ice)
  REAL(dp), SAVE :: calbmxs     ! maximum (snow on ice)
  
CONTAINS
  !
  !-------------------------------------------------------------------------------------------------
  SUBROUTINE init_albedo_ice

    USE mo_control,   ONLY: nn, lcouple, lipcc
    USE mo_exception, ONLY: finish

    ! Define albedo parameters depending on coupled/uncoupled runs and depending on resolution (in coupled case)
    calbmxi = 0.75_dp
    IF(lcouple .OR. lipcc) THEN
       IF (nn == 31) THEN
          calbmns = 0.65_dp
          calbmxs = 0.8_dp
          calbmni = 0.55_dp
       ELSE IF (nn == 42) THEN
          calbmns = 0.70_dp
          calbmxs = 0.85_dp
          calbmni = 0.60_dp
       ELSE IF (nn == 63 .OR. nn == 106 .OR. nn==213) THEN
          calbmns = 0.70_dp
          calbmxs = 0.85_dp
          calbmni = 0.60_dp
       ELSE
          CALL finish ('init_albedo_ice', 'Truncation not supported in coupled runs.')
       ENDIF
    ELSE IF (nn == 319 .OR. nn == 511) THEN
      calbmns = 0.75_dp
      calbmxs = 0.85_dp
      calbmni = 0.65_dp
    ELSE
       calbmns = 0.6_dp
       calbmxs = 0.8_dp
       calbmni = 0.5_dp
    ENDIF

  END SUBROUTINE init_albedo_ice
  !
  !---------------------------------------------------------------------------------------------------------
  SUBROUTINE update_albedo_ice(mask, surf_temp, snow, p_albedo)
    
    USE mo_constants, ONLY:  &
         tmelt                 ! Melting temperature of ice/snow
    
    LOGICAL, INTENT(in)  :: mask(:)
    REAL(dp),    INTENT(in)  :: surf_temp(:), snow(:)
    REAL(dp),    INTENT(out) :: p_albedo(:)
    
    ! Local variables
    REAL(dp), DIMENSION(SIZE(mask)) :: alb_min, alb_max              & ! Minimum and maximum snow albedo
         , temp_upper_limit              & ! Upper temperature limit for cold snow albedo
         , alb_sens                        ! Change of snow albedo per deg C

    temp_upper_limit = tmelt-1.0_dp

    p_albedo = MERGE(0.55_dp,0.0_dp,mask)

    WHERE (mask)
       
       ! Minimum and maximum albedo
       WHERE (snow > 0.01_dp)               ! sea ice covered with snow
          alb_min = calbmns
          alb_max = calbmxs
       ELSEWHERE                         ! bare sea ice
          alb_min = calbmni
          alb_max = calbmxi
       END WHERE
       
       ! Temperature-dependent snow albedo
       WHERE (surf_temp >= tmelt)
          p_albedo = alb_min
       ELSEWHERE (surf_temp < temp_upper_limit)
          p_albedo = alb_max
       ELSEWHERE
          p_albedo = alb_min + &
               (((alb_max-alb_min) / (tmelt-temp_upper_limit)) *    & ! Change of snow albedo per deg C
               (tmelt-surf_temp))
       END WHERE
       
    END WHERE
    
  END SUBROUTINE update_albedo_ice
  !
  !---------------------------------------------------------------------------------------------------------
  SUBROUTINE update_z0_ice(kdim, mask, paz0i)

    USE mo_physc2,         ONLY: cz0ice
    
    INTEGER, INTENT(in)     :: kdim
    LOGICAL, INTENT(in)     :: mask(kdim)
    REAL(dp), INTENT(out)       :: paz0i(kdim)

    paz0i = cz0ice

    WHERE(mask)
       paz0i(:) = cz0ice
    END WHERE

  END SUBROUTINE UPDATE_Z0_ICE
  !---------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE update_stress_ice(     &
       kdim          , mask         &
     , zcfmi       , zudif        &
     , zvdif       , pustri       &
     , pvstri                     &
     )

    USE mo_time_control, ONLY: time_step_len
    USE mo_constants,    ONLY: g
    
    INTEGER, INTENT(in)    :: kdim
    LOGICAL, INTENT(in)    :: mask(kdim)
    REAL(dp), INTENT(in)       :: zcfmi(kdim)
    REAL(dp), INTENT(in)       :: zudif(kdim), zvdif(kdim)
    REAL(dp), INTENT(out)      :: pustri(kdim), pvstri(kdim)

    REAL(dp) :: ztmst, zcons15

    pustri = 0._dp
    pvstri = 0._dp

    ztmst   = time_step_len
    zcons15 = 1._dp / (g * ztmst)

    WHERE(mask)
       pustri(:) = zcons15 * zcfmi(:) * zudif(:)
       pvstri(:) = zcons15 * zcfmi(:) * zvdif(:)
    END WHERE

  END SUBROUTINE update_stress_ice

  !---------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE precalc_ice( &
         kdim &
       , mask &
       , ptsi, paphm1 &
       , pqm1, zx &
       , ptm1, zqss &
       , zteta1, ztvir1 &
       , zfaxe, paclc &
       , zlteta1, pgeom1 &
       , pum1, pvm1 &
       , pocu, pocv &
       , paz0i, zghabl &
       , zqsi, zcpti &
       , zrii, zcfmi &
       , zchi, zcfhi &
       , zbni, zbmi &
       , zbhi, zustari &
       , ztkevi , zusti &
!---wiso-code
       , lwiso, kwiso &
       , pwisosw_d &
       , zwisoqsi &
!---wiso-code-end
       )

    USE mo_constants,    ONLY: rd, cpd, vtmpc1, g, vtmpc2,          &
!---wiso-code
                               tmelt
!---wiso-code-end
    USE mo_physc2,       ONLY: ckap, cc, cb, cvdifts, cfreec, cgam
    USE mo_time_control, ONLY: time_step_len
!---wiso-code
    USE mo_wiso,         ONLY: nwisotyp,                         &
                               talphas1, talphas2, talphas3, talphas4, &
                               tsatbase,                         &
                               toce, tdifrel,                    &
                               tsatfac
!---wiso-code-end

    INTEGER, INTENT(in)    :: kdim
!---wiso-code
    LOGICAL, INTENT(in)    :: lwiso
    INTEGER, INTENT(in)    :: kwiso
!---wiso-code-end
    LOGICAL, INTENT(in)    :: mask(kdim)
    ! ACHTUNG klevp1:  paphm1
    ! ACHTUNG klev:    pqm1, ptm1, zqss, zx, ztvir1, zteta1,zface, paclc, zlteta1, pgeom1,
    REAL(dp),    INTENT(IN)    :: ptsi(kdim), paphm1(kdim), pqm1(kdim), zx(kdim)
    REAL(dp),    INTENT(in)    :: ptm1(kdim), zqss(kdim), zteta1(kdim)
    REAL(dp),    INTENT(IN)    :: ztvir1(kdim), zfaxe(kdim), paclc(kdim)
    REAL(dp),    INTENT(in)    :: zlteta1(kdim), pum1(kdim), pvm1(kdim), paz0i(kdim)
    REAL(dp),    INTENT(in)    :: pocu(kdim), pocv(kdim)
    REAL(dp),    INTENT(IN)    :: pgeom1(kdim), zghabl(kdim)
     
    REAL(dp),    INTENT(out)   :: zqsi(kdim), zcpti(kdim), zrii(kdim)
    REAL(dp),    INTENT(out)   :: zcfmi(kdim), zchi(kdim), zcfhi(kdim)
    REAL(dp),    INTENT(out)   :: zbni(kdim), zbmi(kdim), zbhi(kdim)
    REAL(dp),    INTENT(out)   :: zustari(kdim), ztkevi(kdim), zusti(kdim)
!---wiso-code
    REAL(dp),    INTENT(in), OPTIONAL    :: pwisosw_d(kdim,kwiso)
    REAL(dp),    INTENT(out), OPTIONAL   :: zwisoqsi(kdim,kwiso)
!---wiso-code-end
    
    !  SAVE VARIABLE FOR DIAGNOSE AFTER VERTICAL COLUMN IST CALCULATED WITHIN VDIFF
        
    ! Local Variables
    INTEGER  :: it(kdim)
!---wiso-code
    INTEGER  :: jt, jl
!---wiso-code-end
    
    REAL(dp)     :: zes(kdim), zqmitte(kdim), zqtmit(kdim), ztmit(kdim), zqsmit(kdim), zvirmitte(kdim)
    REAL(dp)     :: ztemitte(kdim), zfux(kdim), zfox(kdim), zmult1(kdim), zqddif(kdim), zbuoy(kdim)
    REAL(dp)     :: zmult2(kdim), zmult3(kdim), zmult4(kdim), zmult5(kdim), zdus1(kdim), zdus2(kdim)
    REAL(dp)     :: zteldif(kdim), zalo(kdim), zaloh(kdim), zcons(kdim), zdthv(kdim) 
    REAL(dp)     :: zucfi(kdim), zscfi(kdim), zcfnchw(kdim), zchnw(kdim), zcr(kdim)
    REAL(dp)     :: zcdn2m(kdim), zcdni(kdim), zcfm2m(kdim), ztesw(kdim), ztvi(kdim), zcdnw(kdim)
    REAL(dp)     :: zcfnci(kdim), ztesi(kdim), zdu2oc(kdim), zcdnr(kdim), zconvs(kdim), zwsti(kdim)
    REAL(dp)     :: zmonob(kdim), zstabf(kdim), zero(kdim)
  
    REAL(dp)     :: zkappa, zkap, zrvrd, zcons9, zcons8, zcons11, zcons12, zepsec, zcons17, zepdu2
    REAL(dp)     :: zrdrv, zepz0o, zepsr, zcons6, zustf
    REAL(dp)     :: zsmn, zshn, zm1, zm2, zm4, zh1, zh2, zwstf

!---wiso-code
    REAL(dp)     :: ztsi_tmp(kdim), zwisofraci(kdim,kwiso), zsatval, zice(kdim,kwiso)
!---wiso-code-end
    
    ! Constants
    zkappa          = rd / cpd
    zkap            = ckap
    zrvrd           = vtmpc1 + 1._dp
    zepdu2          = 1.0_dp
    zcons8          = 2._dp * cb
    zcons9          = 3._dp * cb
    zepsec          = 1.e-2_dp
    zrdrv           = 1._dp / zrvrd
    zepz0o          = 2._dp
    zepsr           = 1.e-10_dp
    
    zcons11         = 3._dp * cb * cc
    zcons12         = cvdifts * time_step_len * g / rd
    zcons17         = 1._dp / zkap**2
    zcons6          = 1._dp / 3._dp
    zh1             = 2.22_dp
    zh2             = 0.22_dp
    zm1             = 1.24_dp
    zm2             = 2.37_dp
    zm4             = 3.69_dp
    zshn            = zh1 * zh2 * SQRT(2._dp)
    zsmn            = zshn * zm1 * zm2 / zm4
    zustf           = 1._dp / zsmn**2
    zwstf           = 0.2_dp
    zero(:)         = 0._dp

    zqsi = 0._dp
    zcpti = 0._dp
    zrii = 0._dp
    zcfmi = 0._dp
    zchi = 0._dp
    zcfhi = 0._dp
    zbni = 0._dp
    zbmi = 0._dp
    zbhi = 0._dp
    zustari = 0._dp
    ztkevi = 0._dp
    zusti = 0._dp

!---wiso-code
    IF (lwiso) THEN

    DO jt=1,kwiso
       zwisoqsi(:,jt) = 0._dp
    END DO
    
    END IF
!---wiso-code-end

    !------------------------------------------------------
    !      2.2   surface humidity and virtual temperature
    WHERE(mask)
       it(:)          = NINT(ptsi(:) * 1000._dp)
!       IF (it(:) < jptlucu1 .OR. it(:) > jptlucu2) lookupoverflow = .TRUE.
       it(:)          = MAX(MIN(it(:), jptlucu2), jptlucu1)
       zes(:)         = tlucua(it(:)) / paphm1(:)
       zqsi(:)        = zes(:) / (1._dp- vtmpc1 * zes(:))
       zcpti(:)       = ptsi(:) * cpd * (1._dp+ vtmpc2 * zqsi(:))
       ztesi(:)       = ptsi(:) * (1.e5_dp / paphm1(:))**zkappa
       ztvi(:)        = ztesi(:) * (1._dp + vtmpc1 * zqsi(:))
    ENDWHERE

!---wiso-code
    IF (lwiso) THEN

    ! surface humidity for ice - water isotopes
    DO jt=1,kwiso
       DO jl=1,kdim
            ztsi_tmp(jl)=ptsi(jl)
            zwisofraci(jl,jt) = (exp(talphas1(jt)/(ztsi_tmp(jl)**2._dp)+talphas2(jt)/ztsi_tmp(jl)+talphas3(jt)))**talphas4(jt)
            IF (nwisotyp(jt).ne.1.and.ztsi_tmp(jl).lt.tmelt) THEN ! effective fractionation over ice if necessary
                zsatval=tsatbase-tsatfac*(ztsi_tmp(jl)-tmelt)
                zwisofraci(jl,jt)=zwisofraci(jl,jt)*(zsatval/(1._dp+zwisofraci(jl,jt)*(zsatval-1._dp)*tdifrel(jt)))
            END IF
       END DO
    END DO

    DO jt=1,kwiso
        WHERE(mask)
            ! "corrections" of seaice surface isotope values be prescibed seasurface delta values
            zice(:,jt)=toce(jt)*(1._dp+(pwisosw_d(:,jt)/1000._dp))
            zwisofraci(:,jt)=zice(:,jt)/zwisofraci(:,jt)
            !calulation of tracer humidity saturation
            zwisoqsi(:,jt)=zwisofraci(:,jt)*zqsi(:)
        END WHERE
    END DO
    
    END IF
!---wiso-code-end

    !     ------------------------------------------------------------------
    !        3.     COMPUTATION OF THE EXCHANGE COEFFICIENTS.
    !        3.1       COMPUTATION OF BASIC QUANTITIES: WIND SHEAR,
    !                  RICHARDSON NUMBER,SQUARED MIXING LENGTHS, UNSTABLE
    !                  AND STABLE CASE COMMON FACTORS AND NEUTRAL CASE
    !                  COMMON PART OF THE DRAG COEFFICIENTS.
    !
    WHERE(mask)
       zdu2oc(:)       = MAX(zepdu2 , pum1(:)**2 + pvm1(:)**2)
    !
    ! correction for water and ice points
    !
       zdu2oc(:)    = MAX(zepdu2,(pum1(:)-pocu(:))**2                  &
                                +(pvm1(:)-pocv(:))**2)
       zcons(:)      = zcons12 * paphm1(:) / (ptm1(:) * (1._dp + vtmpc1 * pqm1(:) - zx(:)))       
       zqmitte(:)    = (pqm1(:) + zqsi(:)) / 2._dp
       zqtmit(:)     = zx(:) * 0.5_dp + zqmitte(:)
       ztmit(:)      = (ptm1(:) + ptsi(:)) / 2._dp
       zqsmit(:)     = (zqss(:) + zqsi(:)) / 2._dp
       ztemitte(:)   = (zteta1(:) + ztesi(:)) / 2._dp
       zvirmitte(:)  = (ztvir1(:) + ztvi(:)) / 2._dp
       zfux(:)       = zfaxe(:) / (cpd * ztmit(:))
       zfox(:)       = zfaxe(:) / (rd * ztmit(:))
       zmult1(:)     = 1._dp + vtmpc1 * zqtmit(:)
       zmult2(:)     = zfux(:) * zmult1(:) - zrvrd
       zmult3(:)     = zrdrv * zfox(:) * zqsmit(:) / (1._dp + zrdrv * zfox(:) * zfux(:) * zqsmit(:))
       zmult5(:)     = zmult1(:) - zmult2(:) * zmult3(:)
       zmult4(:)     = zfux(:) * zmult5(:) - 1._dp
       zdus1(:)      = paclc(:) * zmult5(:) + (1._dp - paclc(:)) * zmult1(:)
       zdus2(:)      = paclc(:) * zmult4(:) + (1._dp - paclc(:)) * vtmpc1
       zteldif(:)    = zlteta1(:) - ztesi(:)
       zqddif(:)     = pqm1(:) + zx(:) - zqsi(:)
       zbuoy(:)      = zdus1(:) * zteldif(:) + zdus2(:) * ztemitte(:) * zqddif(:)
       zrii(:)       = pgeom1(:) * zbuoy(:) / (zvirmitte(:) * zdu2oc(:))
    END WHERE
    WHERE(mask)
       zalo(:)       = LOG(1._dp+ pgeom1(:) / (g * paz0i(:)))
    END WHERE
    WHERE(mask)
       zcdni(:)      = (ckap / zalo(:))**2
    END WHERE
    WHERE(mask)
       zucfi(:)      = 1._dp / (1._dp + zcons11 * zcdni(:) * SQRT(ABS(zrii(:)) &
            * (1._dp + pgeom1(:) / (g * paz0i(:)))))
    END WHERE
    WHERE(mask)
       zscfi(:)      = SQRT(1._dp+ ABS(zrii(:)))
       zcfnci(:)     = zcons(:) * SQRT(zdu2oc(:)) * zcdni(:)
    END WHERE
    WHERE(mask)
       zdthv(:)      = MAX(0._dp,(ztvi(:) - ztvir1(:)))
       zwsti(:)      = zdthv(:) * SQRT(zdu2oc(:)) / zvirmitte(:)
    ENDWHERE
    !----------------------------------------------------------------------------------------------
    !     3.2  DIMENSIONLESS HEAT TRANSFER COEFFICIENTS MULTIPLIED
    !          BY PRESSURE THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE
    !
    WHERE(mask) 
       WHERE(zrii(:).GT.0._dp)
          zcfmi(:)  = zcfnci(:) / (1._dp + zcons8 * zrii(:) / zscfi(:))
          zcfhi(:)  = zcfnci(:) / (1._dp + zcons8 * zrii(:) * zscfi(:))
          zchi(:)   = zcfhi(:) / zcfnci(:) * zcdni(:)
       ELSEWHERE
          zcfmi(:)  = zcfnci(:) * (1._dp - zcons8 * zrii(:) * zucfi(:))
          zcfhi(:)  = zcfnci(:) * (1._dp - zcons9 * zrii(:) * zucfi(:))
          zchi(:)   = zcfhi(:) / zcfnci(:) * zcdni(:)
       ENDWHERE
    ENDWHERE
    !-----------------------------------------------------------------------------------------------
    !     interpolation functions for diagnostics
    !
    WHERE(mask)
       zbni(:)      = zkap / SQRT(zcdni(:))
       zbmi(:)      = MAX(zepsec, SQRT(zcfmi(:) * zcdni(:) * zcons17 / zcfnci(:)))
       zbhi(:)      = MAX(zepsec, zchi(:) / zbmi(:) * zcons17)
       zbmi(:)      = 1._dp / zbmi(:)
       zbhi(:)      = 1._dp / zbhi(:)
    ENDWHERE
    
    !-----------------------------------------------------------------------------------------------
    !*       3.4       COMPUTATION OF THE PBL EXTENSION.
    !
    WHERE(mask)
       WHERE(paz0i(:).GT.zepz0o)
          zcdn2m(:)  = (zkap / LOG(1._dp+ pgeom1(:) / (g * zepz0o)))**2
       ELSEWHERE
          zcdn2m(:)   = zcdni(:)
       ENDWHERE
       zcdnr(:)       = zcdn2m(:) / zcdni(:)
       WHERE( paz0i(:).GT.zepz0o.AND.zrii(:).LT.0._dp)
          zcfm2m(:)   = zcfnci(:) * zcdnr(:) * (1._dp - zcons8 *zrii(:)) &
               / (1._dp + zcons11 * zcdn2m(:) * SQRT(ABS(zrii(:)) &
               * (1._dp + pgeom1(:) / (g * zepz0o))))
       ELSEWHERE
          zcfm2m(:)   = zcfmi(:) * zcdnr(:)
       ENDWHERE
       zusti(:)       = zcfm2m(:) * SQRT(zdu2oc(:))
       zustari(:)     = SQRT(zusti(:) * ptm1(:) &
            * (1._dp + vtmpc1 * pqm1(:) - zx(:)) &
            / (zcons12 * paphm1(:)))
    ENDWHERE
    
    !----------------------------------------------------------------------------------------------
    !      CONVECTIVE VELOCITY SCALE, MONIN-OBUKHOV LENGTH AND
    !      TKE BOUNDARY CONDITION (MAILHOT/BENOIT, 1982)
    WHERE(mask)
       WHERE(zwsti(:) .GT. zepsr)
          zconvs(:)    = (zwsti(:) * zchi(:) * zghabl(:))**zcons6
          zmonob(:)    = (zustari(:)**3) / (zkap * g * zwsti(:) * zchi(:))
          zstabf(:)    = (pgeom1(:) / (g * zmonob(:)))**(zcons6 * 2._dp)
          zstabf(:)    = MIN(zustf * 3._dp, zstabf(:))
       ELSEWHERE
          zconvs=0._dp
          zstabf=0._dp
       ENDWHERE
       ztkevi(:)       = (zustf + zstabf(:)) * (zustari(:)**2) + zwstf * (zconvs(:)**2)
       
    ENDWHERE
  END SUBROUTINE PRECALC_ICE
  !------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE richtmeyr_ice(kdim, mask, &
       klev, klevp1, klevm1, & 
       paphm1, zcfh, &
       zebsh, zqdif, &
       ztdif, zcfhi, &
       zetni, zftni, &
       zeqni, zfqni, &
!---wiso-code
       lwiso, kwiso, &
       ! Inputvariables: diffusion coefficient of humidity
       zwisoqdif, &
       ! Outputvariables: EN and FN coefficient of the Richtmyer-Morton-Scheme, water isotopes
       zwisoeqni, zwisofqni &
!---wiso-code-end
       )

    USE mo_physc2,       ONLY: cvdifts

    INTEGER,   INTENT(IN)  :: kdim, klev, klevp1, klevm1
!---wiso-code
    LOGICAL, INTENT(in)   :: lwiso
    INTEGER, INTENT(in)   :: kwiso
!---wiso-code-end
    LOGICAL,   INTENT(IN)  :: mask(kdim)
    REAL(dp),      INTENT(IN)  :: zebsh(kdim, klev), zcfhi(kdim)
    REAL(dp),      INTENT(IN)  :: paphm1(kdim,klevp1), zcfh(kdim,klev), ztdif(kdim,klev)
    REAL(dp),      INTENT(IN)  :: zqdif(kdim,klev)
!---wiso-code
    REAL(dp), INTENT(in), OPTIONAL  :: zwisoqdif(kdim,klev,kwiso)
!---wiso-code-end
    REAL(dp),      INTENT(OUT) :: zetni(kdim), zftni(kdim), zeqni(kdim), zfqni(kdim)
!---wiso-code
    REAL(dp), INTENT(OUT), OPTIONAL :: zwisoeqni(kdim,kwiso), zwisofqni(kdim,kwiso)
!---wiso-code-end
    
    REAL(dp)   :: zdisci(kdim), zdisqi(kdim)
    REAL(dp)   :: zqdp(kdim), zfac(kdim)
    REAL(dp)   :: ztpfac1
!---wiso-code
    INTEGER    :: jt, jl
    REAL(dp)   :: zwisodisqi(kdim,kwiso)
!---wiso-code-end
    
    !------------------------
    ! CONSTANTS

    ztpfac1   = cvdifts
    zetni = 0._dp
    zftni = 0._dp
    zeqni = 0._dp
    zfqni = 0._dp

!---wiso-code
    IF (lwiso) THEN

    zwisoeqni = 0._dp
    zwisofqni = 0._dp
    
    END IF
!---wiso-code-end

    WHERE(mask)
       zqdp(:)       = 1._dp / (paphm1(:,klevp1) - paphm1(:,klev))
       zfac(:)       = zcfh(:,klevm1) * zqdp(:)
       !-------------------------------------------------------------------
       !*  CALCULATION OF THE EN AND FN COEFFICIENTS OF THE RICHTMYER-
       !*  MORTON-SCHEME CONCERNING THE EQUATION:
       !
       !*  XN = EN * XS + FN
       !
       !*  WITH XN = S_ATM  OR  XN = QATM : ATM. VALUE OF S OR Q
       !*  AND  XS = SSURF  OR  XS = QSAT : SURFACE VALUE OF S OR SAT. SPEC.
       !*                                   HUM. AT THE SURFACE
       !
       zdisci(:)     = 1._dp / (1._dp + zfac(:) * (1._dp - zebsh(:,klevm1)) + zcfhi(:) * zqdp(:))
       zdisqi(:)     = zdisci(:)
       zetni(:)      = zdisci(:) * zcfhi(:) * zqdp(:)
       zftni(:)      = zdisci(:) * (ztdif(:,klev) + zfac(:) * ztdif(:,klevm1)) * ztpfac1
       zeqni(:)      = zdisqi(:) * zcfhi(:) * zqdp(:)
       zfqni(:)      = zdisqi(:) * (zqdif(:,klev) + zfac(:) * zqdif(:,klevm1)) * ztpfac1
    ENDWHERE

!---wiso-code
    IF (lwiso) THEN

    DO jt=1,kwiso
        WHERE(mask)
            zqdp(:)       = 1._dp / (paphm1(:,klevp1) - paphm1(:,klev))
            zfac(:)       = zcfh(:,klevm1) * zqdp(:)
            zwisodisqi(:,jt)    = 1._dp / (1._dp + zfac(:) * (1._dp - zebsh(:,klevm1)) + zcfhi(:) * zqdp(:))
            zwisoeqni(:,jt)     = zwisodisqi(:,jt) * zcfhi(:) * zqdp(:)
            zwisofqni(:,jt)     = zwisodisqi(:,jt) * (zwisoqdif(:,klev,jt) + zfac(:) * zwisoqdif(:,klevm1,jt)) * ztpfac1
        ENDWHERE
    END DO
    
    END IF
!---wiso-code-end

  END SUBROUTINE richtmeyr_ice
  !-----------------------------------------------------------------------------------
  !  
  
  SUBROUTINE update_ice (kdim, mask                     &
       , zetni, zftni                                   &
       , zeqni, zfqni                                   &
       , zcpti, zqsi                                    &
       , ztklevi, zqklevi                               &
!---wiso-code
       , lwiso, kwiso &
       ! Input - calculated in subroutine "Richtmeyer-ice"
       , zwisoeqni, zwisofqni                           &
       ! Input - calculated in subroutine "precalc_ice"
       , zwisoqsi                                       &
       ! Output
       , zwisoqklevi &
!---wiso-code-end
       )

    INTEGER, INTENT(IN) :: kdim
!---wiso-code
    LOGICAL, INTENT(IN) :: lwiso
    INTEGER, INTENT(IN) :: kwiso
!---wiso-code-end
    LOGICAL, INTENT(IN) :: mask(kdim)
    REAL(dp), INTENT(IN)    :: zetni(kdim), zftni(kdim), zeqni(kdim), zfqni(kdim)
    REAL(dp), INTENT(IN)    :: zcpti(kdim), zqsi(kdim)
!---wiso-code
    REAL(dp), INTENT(IN), OPTIONAL    :: zwisoeqni(kdim,kwiso), zwisofqni(kdim,kwiso)
    REAL(dp), INTENT(IN), OPTIONAL    :: zwisoqsi(kdim,kwiso)
!---wiso-code-end
    REAL(dp), INTENT(OUT)   :: ztklevi(kdim), zqklevi(kdim)
!---wiso-code
    REAL(dp), INTENT(OUT), OPTIONAL   :: zwisoqklevi(kdim,kwiso)

    ! local variable
    INTEGER      :: jt, jl
!---wiso-code-end
    
    !*  CALCULATION OF SKLEV AND QKLEV USING THE NEW SURFACE VALUES
    !*  ZSNEW AND ZQSNEW WHICH WERE CALCULATED IN SUBROUTINE SURFTEMP
    
    ztklevi = 0._dp
    zqklevi = 0._dp
!---wiso-code
    IF (lwiso) THEN
    
    zwisoqklevi = 0._dp
    
    END IF
!---wiso-code-end

    WHERE(MASK)
       ztklevi(:)   = zetni(:) * zcpti(:) + zftni(:)
       zqklevi(:)   = zeqni(:) * zqsi(:) + zfqni(:)
    ENDWHERE

!---wiso-code
    IF (lwiso) THEN

    DO jt=1,kwiso
        WHERE(MASK)
            zwisoqklevi(:,jt) = zwisoeqni(:,jt) * zwisoqsi(:,jt) + zwisofqni(:,jt)
        ENDWHERE
    END DO
    
    END IF
!---wiso-code-end

  END SUBROUTINE update_ice
  !--------------------------------------------------------------------------------------
  !
  SUBROUTINE postproc_ice(kdim, klev &
       , klevp1, mask                &
       , zcfhi, zqsi                 &
       , zqdif, pqm1                 &
       , pgeom1, ztdif               &
       , zcptgz, zcpti               &
       , ptsi , zbni                 &
       , zbhi, zrii                  &
       , ptm1, papm1                 &
       , zx, pum1                    &
       , pvm1, pevapi                &
       , pahfsi, pahfli              &
       , zspeedi, zt2i               &
       , paphm1, zbmi                &
       , zdew2i, zu10i               &
       , zv10i                       &
!---wiso-code
       ! Input
       , lwiso, kwiso                  &
       , zwisoqsi, zwisoqdif, pwisoqm1 &
       ! Output
       , pwisoevapi                    &
!---wiso-code-end
       )

    USE mo_constants,    ONLY: vtmpc1, vtmpc2, g, cpd, alv, als, tmelt, c3les, c4les, c3ies, c4ies, rd &
         ,c2es
    USE mo_physc2,       ONLY: cvdifts
    USE mo_time_control, ONLY: time_step_len

    INTEGER, INTENT(in)  :: kdim, klev, klevp1
!---wiso-code
    LOGICAL, INTENT(in)  :: lwiso
    INTEGER, INTENT(in)  :: kwiso
!---wiso-code-end
    LOGICAL, INTENT(in)  :: mask(kdim)
    REAL(dp), INTENT(in)     :: zcfhi(kdim), zqsi(kdim)
    REAL(dp), INTENT(in)     :: zqdif(kdim), pqm1(kdim), pgeom1(kdim) 
    REAL(dp), INTENT(in)     :: ztdif(kdim), zcptgz(kdim)
    REAL(dp), INTENT(in)     :: zcpti(kdim), ptsi(kdim), zbni(kdim), zbhi(kdim), zrii(kdim)
    REAL(dp), INTENT(in)     :: ptm1(kdim), papm1(kdim), zx(kdim)
    REAL(dp), INTENT(in)     :: pum1(kdim), pvm1(kdim), paphm1(kdim,klevp1)
    REAL(dp), INTENT(in)     :: zbmi(kdim)
!---wiso-code
    REAL(dp), INTENT(in), OPTIONAL     :: zwisoqsi(kdim,kwiso), zwisoqdif(kdim,kwiso)
    REAL(dp), INTENT(in), OPTIONAL     :: pwisoqm1(kdim,klev,kwiso)
!---wiso-code-end

    REAL(dp), INTENT(out)    :: pevapi(kdim), pahfsi(kdim), pahfli(kdim), zspeedi(kdim), zt2i(kdim)
    REAL(dp), INTENT(out)    :: zdew2i(kdim), zu10i(kdim), zv10i(kdim)
!---wiso-code
    REAL(dp), INTENT(out), OPTIONAL    :: pwisoevapi(kdim,kwiso)
!---wiso-code-end

    INTEGER  :: it(kdim)
    REAL(dp)     :: ztmst, zcons15, ztpfac1, ztpfac2, ztpfac3,zcons16, zhuv, zhtq, zephum
    REAL(dp)     :: zcoefi(kdim), zqnlev(kdim), zqhfli(kdim)
    REAL(dp)     :: ztnlev(kdim), zthfli(kdim)
    REAL(dp)     :: zrat(kdim), zcbn(kdim), zcbs(kdim), zcbu(kdim), zmerge(kdim), zred(kdim)
    REAL(dp)     :: zh2m(kdim), zqs1(kdim), zrh2m(kdim), zcvm3(kdim), zcvm4(kdim)
    REAL(dp)     :: zaph2m(kdim), zqs2(kdim), zq2m(kdim), zfrac(kdim)
    REAL(dp)     :: zmerge1(kdim)
!---wiso-code
    INTEGER      :: jt, jl
    REAL(dp)     :: zwisoqnlev(kdim,kwiso)
    REAL(dp)     :: zwisoqhfli(kdim,kwiso)
!---wiso-code-end

    !CONSTANTS
    ztmst         = time_step_len
    ztpfac1       = cvdifts
    ztpfac2       = 1._dp / ztpfac1
    ztpfac3       = 1._dp - ztpfac2
    zcons15       = 1._dp / (g * ztmst)
    zcons16       = cpd * vtmpc2
    zhuv          = 10._dp * g
    zhtq          = 2._dp * g
    zephum        = 5.e-2_dp

    pevapi = 0._dp
    pahfsi = 0._dp
    pahfli = 0._dp
    zspeedi = 0._dp
    zt2i = 0._dp
    zdew2i = 0._dp
    zu10i =  0._dp
    zv10i = 0._dp

!---wiso-code
    IF (lwiso) THEN

    pwisoevapi     = 0._dp
    
    END IF
!---wiso-code-end

    WHERE(mask)

       !*         5.8     Surface fluxes of heat and moisture
       !*    Moisture fluxes
       !
       zcoefi(:)      = zcons15 * zcfhi(:)
       zqnlev(:)      = zqdif(:) ! - ztpfac3 * pqm1(:)
       zqhfli(:)      = zcoefi(:) * (zqnlev(:) - ztpfac2 * zqsi(:))
       !*    Sensible heat fluxes
       ztnlev(:)      = ztdif(:) ! - ztpfac3 * zcptgz(:)
       zthfli(:)      = zcoefi(:) * (ztnlev(:) - ztpfac2 * zcpti(:))
       zthfli(:)      = zthfli(:) - ptsi(:) * zcons16 * zqhfli(:)
       !    Accumulated sensible heat flux and evaporation
       pevapi(:)      = zqhfli(:)
       pahfsi(:)      = zthfli(:)
       !     Latent heat fluxes
       pahfli(:)      = als * zqhfli(:)

       !
       !     Compute new t2m, t2m_max t2m_min
       !
       zrat(:)        = zhtq / pgeom1(:)
       zcbn(:)        = LOG(1._dp + (EXP (zbni(:)) - 1._dp) * zrat(:))
       zcbs(:)        = -(zbni(:) - zbhi(:)) * zrat(:)
       zcbu(:)        = -LOG(1.+ (EXP (zbni(:) - zbhi(:)) - 1._dp) * zrat(:))
       WHERE(zrii(:).GT.0._dp)
          zmerge(:)   = zcbs(:)
       ELSEWHERE
          zmerge(:)   = zcbu(:)
       ENDWHERE
       zred(:)        = (zcbn(:) + zmerge(:)) / zbhi(:)
       zh2m(:)        = zcpti(:) + zred(:) * (zcptgz(:) - zcpti(:))
       zt2i(:)        = (zh2m(:) - zhtq) / (cpd * (1._dp + vtmpc2 * pqm1(:)))  
       !
       !           5.96   2M DEW POINT
       !
       it(:)          = NINT(ptm1(:) * 1000._dp)
       !    WHERE (it(:) < jptlucu1 .OR. it(:) > jptlucu2) 
       !       lookupoverflow = .TRUE.
       !    ENDWHERE
       it(:)          = MAX(MIN(it,jptlucu2),jptlucu1)
       zqs1(:)        = tlucua(it(:)) / papm1(:)
       zqs1(:)        = zqs1(:) / (1._dp - vtmpc1 * zqs1(:))
       zrh2m(:)       = MAX(zephum, pqm1(:) / zqs1(:))
       WHERE(zt2i(:) .GT. tmelt)
          zcvm3(:)    = c3les
          zcvm4(:)    = c4les
       ELSEWHERE
          zcvm3(:)    = c3ies
          zcvm4(:)    = c4ies
       ENDWHERE
       zaph2m(:)      = paphm1(:,klevp1) * &
            (1._dp - zhtq / (rd * zt2i(:) * (1._dp + vtmpc1 * pqm1(:) - zx(:))))
       it(:)          = NINT(zt2i(:) * 1000._dp)
       !    WHERE (it(:) < jptlucu1 .OR. it(:) > jptlucu2) lookupoverflow = .TRUE.
       it(:)          = MAX( MIN( it(:), jptlucu2), jptlucu1)
       zqs2(:)        = tlucua(it(:)) / zaph2m(:)
       zqs2(:)        = zqs2(:) / (1._dp - vtmpc1 * zqs2(:))
       zq2m(:)        = zrh2m(:) * zqs2(:)
       zfrac(:)       = LOG(zaph2m(:) * zq2m(:) / (c2es * (1._dp + vtmpc1 * zq2m(:)))) / zcvm3(:)
       zdew2i(:)      = MIN(zt2i(:), (tmelt - zfrac(:) * zcvm4(:)) / (1._dp - zfrac(:)))
       !
       !*          5.97   10M WIND COMPONENTS, MAX 10M WINDSPEED
       !
       zrat(:)        = zhuv / pgeom1(:)
       zcbn(:)        = LOG(1._dp + (EXP (zbni(:)) - 1._dp) * zrat(:))
       zcbs(:)        = -(zbni(:) - zbmi(:)) * zrat(:)
       zcbu(:)        = -LOG(1._dp + (EXP (zbni(:) - zbmi(:)) - 1._dp) * zrat(:))
       WHERE(zrii(:).GT.0._dp)
          zmerge1(:)  = zcbs(:)
       ELSEWHERE
          zmerge1(:)  = zcbu(:)
       ENDWHERE
       zred(:)        = (zcbn(:) + zmerge1(:)) / zbmi(:)
       zu10i(:)       = zred(:) * pum1(:)
       zv10i(:)       = zred(:) * pvm1(:)
       zspeedi(:)     = SQRT(zu10i(:)**2 + zv10i(:)**2)
      
    ENDWHERE

!---wiso-code
    IF (lwiso) THEN

    DO jt=1,kwiso
        WHERE(mask)
            ! Surface fluxes of moisture over the ice
            zwisoqnlev(:,jt) = zwisoqdif(:,jt) ! - ztpfac3*pwisoqm1(:,klev,jt)
            zwisoqhfli(:,jt) = zcoefi(:) * (zwisoqnlev(:,jt) - ztpfac2 * zwisoqsi(:,jt))
            ! Evaporation water isotopes
            pwisoevapi(:,jt) = zwisoqhfli(:,jt)
        ENDWHERE
    END DO
    
    END IF
!---wiso-code-end

  END SUBROUTINE postproc_ice

  SUBROUTINE s_licetemp(                &
       kdim          , psiced         &
     , psni        , palake         &
     , ptsi        , ptrfli         &
     , psofli      , pahfice        &
     , pfluxres    , pahfcon        &
     , pahfres     , pevapi         &
     , psnowd      , pahfsi         &
     , pahfli      , pcvsi          &
     , pfri                         &
     )

    USE mo_param_switches, ONLY: lice
    USE mo_constants,      ONLY: tmelt, rhoh2o, alf
    USE mo_time_control,   ONLY: delta_time

    INTEGER, INTENT(in) ::  kdim
    REAL(dp), INTENT(in)    ::  psiced(kdim), palake(kdim)
    REAL(dp), INTENT(inout) ::  psni(kdim),  ptsi(kdim)          
    REAL(dp), INTENT(in)    ::  ptrfli(kdim), psofli(kdim)          
    REAL(dp), INTENT(inout) ::  pahfice(kdim),     pfluxres(kdim)                          
    REAL(dp), INTENT(inout) ::  pahfcon(kdim),     pahfres(kdim)
    REAL(dp), INTENT(in)    ::  pahfsi(kdim),      pahfli(kdim),     pcvsi(kdim)           
    REAL(dp), INTENT(in)    ::  pfri(kdim), psnowd(kdim),    pevapi(kdim)          
  

!!$  REAL(dp) :: zdtime, zalpha, zalphas, zrho_sn, ziscond, zcpice, zrhoice, zcpcon, zcpdt, zdice
!!$  REAL(dp) :: zevsnd(kdim), zicefl(kdim), zsflx(kdim), zmelfac(kdim), zsmelt(kdim)
!!$  REAL(dp) :: zsniced(kdim)
    REAL(dp) :: zdtime, zalpha, zalphas, zrho_sn, ziscond, zcpice, zrhoice, zcpcon, zcpdt, zdice
    REAL(dp) :: zevsnd, zicefl, zsflx, zmelfac, zsmelt
    REAL(dp) :: zsniced, zsnowd
    INTEGER  :: jl

    ! CONSTANTS

    zdtime                    = delta_time
    zalpha                    = 2.1656_dp
    zalphas                   = 0.31_dp
    zrho_sn                   = 330._dp
    ziscond                   = zalpha / zalphas * rhoh2o / zrho_sn
    zcpice                    = 2106._dp
    zrhoice                   = 910._dp
    zdice                     = 0.10_dp
    zcpcon                    = zrhoice * zcpice * zdice
    zcpdt                     = zcpcon / zdtime
    IF (lice) THEN
       !
       DO jl=1,kdim
          IF (palake(jl).GE.0.5_dp) THEN
             IF (psiced(jl).GE.zdice) THEN                         ! ice
                zsnowd=psnowd(jl)*zdtime/rhoh2o
                zevsnd=pcvsi(jl)*pevapi(jl)*zdtime/rhoh2o
                psni(jl)=MAX(psni(jl)+zsnowd+zevsnd,0._dp)
                zsniced=psiced(jl)+ziscond*psni(jl)
                zicefl=zalpha*tmelt/zsniced
                zsflx=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)+pfluxres(jl)
                pfluxres(jl)=0._dp
                ptsi(jl)=(zcpdt*ptsi(jl)+zsflx+zicefl)/                  &
                                                (zcpdt+zalpha/zsniced)
                IF (ptsi(jl).GT.tmelt) THEN
                   zmelfac=(zcpdt+zalpha/zsniced)*zdtime/(alf*rhoh2o)
                   zsmelt=MIN(zmelfac*(ptsi(jl)-tmelt),psni(jl))
                   psni(jl)=psni(jl)-zsmelt
                   zsniced=psiced(jl)+ziscond*psni(jl)
                   ptsi(jl)=ptsi(jl)-zsmelt/zmelfac
                   pahfres(jl)=pahfres(jl)+zsmelt*alf*rhoh2o*pfri(jl)
                END IF
                IF (ptsi(jl).GT.tmelt) THEN
                   pfluxres(jl)=(zcpdt+zalpha/zsniced)*(ptsi(jl)-tmelt)
                   ptsi(jl)=tmelt
                END IF
                pahfres(jl)=pahfres(jl)+zdtime*pfri(jl)*pfluxres(jl)
                pahfice(jl)=zalpha*(ptsi(jl)-tmelt)/zsniced
             ELSE                                                 ! water
                pahfice(jl)=0._dp
                ptsi(jl)=tmelt
                psni(jl)=0._dp
             END IF
             pahfcon(jl)=pahfcon(jl)+zdtime*pfri(jl)*pahfice(jl)
          END IF
       END DO
       !
       !        Necessary computations if subroutine is bypassed
       !
    ELSE
       DO jl = 1, kdim
          ptsi(jl)=tmelt
          psni(jl)=0._dp
       END DO
    END IF
!!$  IF (lice) THEN
!!$        WHERE (palake(:) .GE. 0.5)
!!$           WHERE (psiced(:) .GE. zdice)           !---------------------------------------------------ice
!!$!! convective and large scale direct: psnowd(:)     = (pssfl(jl)+pssfc(jl))*zdtime/rhoh2o
!!$              zevsnd(:)     = pcvsi(:) * pevapi(:) * zdtime / rhoh2o
!!$              psni(:)       = MAX(psni(:) + psnowd(:) + zevsnd(:) , 0.)
!!$              zsniced(:)    = psiced(:) + ziscond * psni(:)
!!$              zicefl(:)     = zalpha * tmelt / zsniced(:)
!!$              zsflx(:)      = ptrfli(:) + psofli(:) + pahfsi(:) + pahfli(:) + pfluxres(:)
!!$              pfluxres(:)   = 0.
!!$              ptsi(:)       = (zcpdt * ptsi(:) + zsflx(:) + zicefl(:)) / (zcpdt + zalpha / zsniced(:))
!!$              WHERE (ptsi(:) .GT. tmelt)
!!$                 zmelfac(:) = (zcpdt + zalpha / zsniced(:)) * zdtime / (alf * rhoh2o)
!!$                 zsmelt(:)  = MIN(zmelfac(:) * (ptsi(:) - tmelt), psni(:))
!!$                 psni(:)    = psni(:) - zsmelt(:)
!!$                 zsniced(:) = psiced(:) + ziscond * psni(:)
!!$                 ptsi(:)    = ptsi(:) - zsmelt(:) / zmelfac
!!$                 pahfres(:) = pahfres(:) + zsmelt(:) * alf * rhoh2o * pfri(:)
!!$              END WHERE
!!$              WHERE (ptsi(:) .GT. tmelt) 
!!$                 pfluxres(:)= (zcpdt + zalpha / zsniced(:)) * (ptsi(:) - tmelt)
!!$                 ptsi(:)    = tmelt
!!$              END WHERE
!!$              pahfres(:)    = pahfres(:) + zdtime * pfri(:) * pfluxres(:)
!!$              pahfice(:)    = zalpha * (ptsi(:) - tmelt) / zsniced(:)
!!$           ELSEWHERE                               ! ------------------------------------------------water
!!$              pahfice(:)    = 0.
!!$              ptsi(:)       = tmelt
!!$              psni(:)       = 0.
!!$           END WHERE
!!$           pahfcon(:)       = pahfcon(:) + zdtime * pfri(:) * pahfice(:)
!!$        END WHERE
!!$     
!!$     !----------------------------------------------------Necessary computations if subroutine is bypassed
!!$     
!!$  ELSE
!!$        ptsi(:)             = tmelt
!!$        psni(:)             = 0.
!!$  END IF

  END SUBROUTINE s_licetemp

  SUBROUTINE s_sicetemp(            &
       kdim          , psiced       &
     , psni        , palake       &
     , pslf                       &
     , ptsi        , ptrfli       &
     , psofli      , pahfice      &
     , pfluxres    , pqres        &
     , pahfcon     , pahfres      &
     , pahfsi      , pahfli       &
     , pfri        , mask         &
     )

    USE mo_param_switches, ONLY: lice
    USE mo_constants,      ONLY: tmelt, rhoh2o, alf
    USE mo_time_control,   ONLY: delta_time
    USE mo_physc2,         ONLY: ctfreez
    USE mo_control,        ONLY: lmlo, lcouple

    INTEGER, INTENT(in ) :: kdim
    REAL(dp), INTENT(in)     :: psiced(kdim), palake(kdim)          
    REAL(dp), INTENT(in)     :: pslf(kdim)
    REAL(dp), INTENT(inout)  :: ptsi(kdim), pqres(kdim)  
    REAL(dp), INTENT(in)     :: ptrfli(kdim), psofli(kdim)          
    REAL(dp), INTENT(inout)  :: pahfice(kdim), pfluxres(kdim), psni(kdim)        
    REAL(dp), INTENT(inout)  :: pahfcon(kdim), pahfres(kdim)   
    REAL(dp), INTENT(in)     :: pahfsi(kdim), pahfli(kdim)        
    REAL(dp), INTENT(in)     :: pfri(kdim)
    LOGICAL, INTENT(in)  :: mask(kdim)

    REAL(dp) :: zdtime, zalpha, zalphas, zrho_sea, zrho_sn, ziscond
    REAL(dp) :: zcpice, zrhoice, zdice, zcpcon, zcpdt

!!$  REAL(dp) :: zsniced(kdim), zicefl(kdim)
    REAL(dp) :: zsniced, zicefl, zsflx(kdim)
    INTEGER :: jl

    ! CONSTANTS

    zdtime                   = delta_time
    zalpha                   = 2.1656_dp
    zalphas                  = 0.31_dp
    zrho_sea                 = 1025._dp
    zrho_sn                  = 330._dp
    ziscond                  = zalpha / zalphas * zrho_sea / zrho_sn
    zcpice                   = 2106._dp
    zrhoice                  = 910._dp
    zdice                    = 0.10_dp
    zcpcon                   = zrhoice * zcpice * zdice
    zcpdt                    = zcpcon / zdtime
  
    zdtime = delta_time
    zalpha = 2.1656_dp
    zalphas=0.31_dp
    zrho_sea=1025._dp
    zrho_sn=330._dp
    ziscond=zalpha/zalphas*zrho_sea/zrho_sn
    zcpice = 2106._dp
    zrhoice = 910._dp
    zdice = 0.10_dp
    zcpcon = zrhoice*zcpice*zdice
    zcpdt = zcpcon/zdtime

!-- 2. Compute new skin-temperature
    
    IF (lice) THEN
       !
       IF (.NOT.lmlo .AND. .NOT.lcouple) THEN     ! ECHAM5
          DO jl=1,kdim
             IF (palake(jl).EQ.0._dp) THEN
                IF (psiced(jl).GE.zdice) THEN
                   zsniced=psiced(jl)+ziscond*psni(jl)
                   zicefl=zalpha*ctfreez/zsniced
                   zsflx(jl)=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)
                   ptsi(jl)=(zcpdt*ptsi(jl)+zsflx(jl)+zicefl)/                    &
                        (zcpdt+zalpha/zsniced)
                   IF (ptsi(jl).GT.tmelt) THEN
                      pqres(jl)=(zalpha/zsniced+zcpdt)*(ptsi(jl)-tmelt)
                      ptsi(jl)=tmelt
                   ELSE
                      pqres(jl)=0._dp
                   END IF
                   pahfice(jl)=zalpha*(ptsi(jl)-ctfreez)/zsniced
                ELSE
                   pqres(jl)=0._dp
                   pahfice(jl)=0._dp
                   ptsi(jl)=tmelt
                END IF
                pahfcon(jl)=pahfcon(jl)+zdtime*pfri(jl)*pahfice(jl)
                pahfres(jl)=pahfres(jl)+zdtime*pfri(jl)*pqres(jl)
             END IF
          END DO

       ELSE IF (lcouple) THEN                      ! ECHAM5-IPCC
          DO jl=1,kdim
             IF (palake(jl).EQ.0._dp .AND. pslf(jl).LT.1.0_dp) THEN
                zsniced=MAX(psiced(jl),zdice)+ziscond*psni(jl)
                zicefl=zalpha*ctfreez/zsniced
                zsflx(jl)=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)
                ptsi(jl)=(zcpdt*ptsi(jl)+zsflx(jl)+zicefl)/                &
                                        (zcpdt+zalpha/zsniced)
                IF (ptsi(jl).GT.tmelt) THEN
                   pqres(jl)=(zalpha/zsniced+zcpdt)*(ptsi(jl)-tmelt)
                   ptsi(jl)=tmelt
                ELSE
                   pqres(jl)=0._dp
                END IF
                pahfice(jl)=zalpha*(ptsi(jl)-ctfreez)/zsniced
                pahfcon(jl)=pahfcon(jl)+zdtime*pfri(jl)*pahfice(jl)
                pahfres(jl)=pahfres(jl)+zdtime*pfri(jl)*pqres(jl)
             ELSE
                pqres(jl)=0._dp
                pahfice(jl)=0._dp
                ptsi(jl)=tmelt
             END IF
          END DO

       ELSE         ! lmlo

          DO jl=1,kdim
             IF (palake(jl).EQ.0._dp) THEN
                IF (psiced(jl).GE.zdice) THEN                         ! ice
                   zsniced=psiced(jl)+ziscond*psni(jl)
                   zicefl=zalpha*ctfreez/zsniced
                   zsflx(jl)=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)        &
                        +pfluxres(jl)
                   pfluxres(jl)=0._dp
                   ptsi(jl)=(zcpdt*ptsi(jl)+zsflx(jl)+zicefl)/                  &
                        (zcpdt+zalpha/zsniced)
                   IF (ptsi(jl).GT.tmelt) THEN
                      pfluxres(jl)=(zcpdt+zalpha/zsniced)*(ptsi(jl)-tmelt)
                      ptsi(jl)=tmelt
                   END IF
                   pahfres(jl)=pahfres(jl)+zdtime*pfri(jl)*pfluxres(jl)
                   pahfice(jl)=zalpha*(ptsi(jl)-ctfreez)/zsniced
                ELSE                                                 ! water
                   pahfice(jl)=0._dp
                   ptsi(jl)=tmelt
                   psni(jl)=0._dp
                END IF
                pahfcon(jl)=pahfcon(jl)+zdtime*pfri(jl)*pahfice(jl)
             END IF
          END DO
       END IF
!
!       Necessary computations if subroutine is bypassed
!
    ELSE
       DO jl = 1, kdim
          ptsi(jl)=tmelt
          psni(jl)=0._dp
       END DO
    END IF

  END SUBROUTINE s_sicetemp

  SUBROUTINE ice_rad(               &
       kdim          , mask         &
      ,ztrdown      , ptsi         &
      ,pi0          , ptrsol       &
      ,palsoi       , palbedo      &
      ,psofli       , ptrfli       &
      )

    USE mo_radiation,         ONLY: cemiss
    USE mo_constants,         ONLY: stbo

    INTEGER, INTENT(in)    :: kdim
    LOGICAL, INTENT(in)    :: mask(kdim)
    REAL(dp),    INTENT(in)    :: ptsi(kdim), pi0(kdim), ptrsol(kdim)
    REAL(dp),    INTENT(in)    :: palsoi(kdim), ztrdown(kdim), palbedo(kdim)
    REAL(dp),    INTENT(out)   :: psofli(kdim), ptrfli(kdim)

    REAL(dp) :: zflxs(kdim)

    psofli = 0._dp
    ptrfli = 0._dp

    WHERE(mask)    
       ptrfli(:)  = ztrdown(:) - cemiss * stbo * ptsi(:)**4
       zflxs(:)   = pi0(:) * ptrsol(:)     
       psofli(:)  = (1._dp - palsoi(:)) * zflxs(:) / (1._dp - palbedo(:))
    END WHERE

  END SUBROUTINE ice_rad

END MODULE mo_surface_ice
