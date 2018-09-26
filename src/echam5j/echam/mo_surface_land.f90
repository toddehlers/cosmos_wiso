MODULE MO_SURFACE_LAND
  
  USE mo_convect_tables,   ONLY: tlucua, jptlucu1, jptlucu2, &
       lookuperror, lookupoverflow
!   USE mo_surface_memory,       ONLY: stest0, stest1, stest2, stest3, stest4, stest5,&
!       stest6, stest7, stest8, stest9
   USE mo_kind, ONLY: dp

  IMPLICIT NONE
  
CONTAINS
  !
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_stress_land(  &
       kdim          , mask       &
       , zcfml       , zudif      &
       , zvdif       , pustrl     &
       , pvstrl                   &
       )

    USE mo_time_control, ONLY: time_step_len
    USE mo_constants,    ONLY: g
    
    INTEGER, INTENT(in)    :: kdim
    LOGICAL, INTENT(in)    :: mask(kdim)
    REAL(dp), INTENT(in)       :: zcfml(kdim)
    REAL(dp), INTENT(in)       :: zudif(kdim), zvdif(kdim)
    REAL(dp), INTENT(out)      :: pustrl(kdim), pvstrl(kdim)

    REAL(dp) :: ztmst, zcons15

    pustrl = 0._dp
    pvstrl = 0._dp

    ztmst   = time_step_len
    zcons15 = 1._dp / (g * ztmst)

    WHERE(mask)
       pustrl(:) = zcons15 * zcfml(:) * zudif(:)
       pvstrl(:) = zcons15 * zcfml(:) * zvdif(:)
    END WHERE

  END SUBROUTINE UPDATE_STRESS_LAND
  !
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE precalc_land(  &
         kdim		    &
       , mask    	    &
       , jrow               &
       , ptslm1   , paphm1  &
       , pum1     , pvm1    &
       , pqm1     , zx      &
       , zqss     , zteta1  &
       , ztvir1   , zfaxe   &
       , paclc    , zlteta1 &
       , pgeom1   , paz0lh, paz0lm &
       , ptm1     , zghabl  &              
       , zdqsl    , zril    &
       , zqsl     , zcfncl  &
       , zchl     , zcfhl   &
       , zbnl     , zbhnl   &
       , zbml     , zbhl    &
       , zustarl  , ztkevl  &
       , zhsoil   , zcfml   &
       , zustl    , zcptl   &
       , zcpq     , zcair   &
       , zcsat              &
       )
    
    USE mo_constants,    ONLY: rd, cpd, vtmpc1, g, vtmpc2
    USE mo_physc2,       ONLY: ckap, cc, cb, cvdifts, cfreec, cgam
    USE mo_time_control, ONLY: time_step_len

    INTEGER, INTENT(in)    :: kdim, jrow
    LOGICAL, INTENT(in)    :: mask(kdim)
    REAL(dp),    INTENT(IN)    :: ptslm1(kdim), paphm1(kdim), pum1(kdim), pvm1(kdim), zhsoil(kdim)
    REAL(dp),    INTENT(in)    :: pqm1(kdim), zx(kdim), zqss(kdim)
    REAL(dp),    INTENT(in)    :: zteta1(kdim), ztvir1(kdim)
    REAL(dp),    INTENT(in)    :: zfaxe(kdim), zlteta1(kdim), zghabl(kdim)
    REAL(dp),    INTENT(in)    :: pgeom1(kdim), paz0lh(kdim), paz0lm(kdim), paclc(kdim), ptm1(kdim)
    REAL(dp),    INTENT(in)    :: zcsat(kdim), zcair(kdim)
    
    REAL(dp),    INTENT(out)   :: zdqsl(kdim), zril(kdim)
    REAL(dp),    INTENT(out)   :: zqsl(kdim) 
    REAL(dp),    INTENT(out)   :: zcfncl(kdim), zchl(kdim), zcfhl(kdim)
    REAL(dp),    INTENT(out)   :: zbnl(kdim), zbhnl(kdim), zbml(kdim), zbhl(kdim), zustarl(kdim)
    REAL(dp),    INTENT(out)   :: ztkevl(kdim), zustl(kdim), zcptl(kdim), zcpq(kdim)
    REAL(dp),    INTENT(out)   :: zcfml(kdim)
    
    INTEGER   :: it(kdim), it1(kdim), i
    REAL(dp)      :: zqls(kdim), zqst1(kdim), ztesl(kdim), ztvl(kdim) 
    REAL(dp)      :: zdu2(kdim), zqmitte(kdim), zqmit(kdim), ztmit(kdim), zqsmit(kdim)
    REAL(dp)      :: ztemitte(kdim)
    REAL(dp)      :: zvirmitte(kdim), zfux(kdim), zfox(kdim), zmult1(kdim),zmult2(kdim), zmult3(kdim)
    REAL(dp)      :: zmult4(kdim), zmult5(kdim), zdus1(kdim), zdus2(kdim), zteldif(kdim), zqddif(kdim)
    REAL(dp)      :: zbuoy(kdim), zchsnow(kdim), zchland(kdim), zchmean(kdim), zcdnl(kdim)
    REAL(dp)      :: zchnl(kdim), zucfl(kdim), zucfhl(kdim), zscfl(kdim), zcons(kdim), zcfnchl(kdim)
    REAL(dp)      :: zdthv(kdim), zcdn2m(kdim), zcdnr(kdim), zcfm2m(kdim), zes(kdim), zqtmit(kdim)
    REAL(dp)      :: zwstl(kdim), zconvs(kdim), zmonob(kdim), zstabf(kdim)
    
    REAL(dp)      :: zcpd, zrd, zkappa, zepdu2, zrvrd, zrdrv, zblend, zkap, zcons11, zcons12, ztmst
    REAL(dp)      :: ztpfac1, zcons8, zepsec, zcons17, zepz0o, zcons9, zwslim, zepsr, zcons6, zustf
    REAL(dp)      :: zsmn, zshn, zm1, zm2, zm4, zh1, zh2, zwstf, zplmin

    ! Constants
    ztmst         = time_step_len
    ztpfac1       = cvdifts
    zepsec        = 1.e-2_dp
    zplmin        = 0.35_dp
    zwslim        = zplmin
    zepsr         = 1.e-10_dp
    
    zcpd          = cpd
    zrd           = rd
    zkappa        = zrd/zcpd
    zepdu2        = 1.0_dp
    zrvrd         = vtmpc1+1._dp
    zrdrv         = 1._dp/ zrvrd
    zblend        = 100._dp
    zkap          = ckap
    zcons11       = 3._dp * cb * cc
    zcons12       = ztpfac1 * ztmst * g / rd
    zcons8        = 2._dp * cb
    zcons9        = 3._dp * cb
    zcons17       = 1._dp / zkap**2
    zepz0o        = 2._dp
    zcons6        = 1._dp / 3._dp
    zh1           = 2.22_dp
    zh2           = 0.22_dp
    zm1           = 1.24_dp
    zm2           = 2.37_dp
    zm4           = 3.69_dp
    zshn          = zh1 * zh2 * SQRT(2._dp)
    zsmn          = zshn * zm1 * zm2 / zm4
    zustf         = 1._dp / zsmn**2
    zwstf         = 0.2_dp

    zdqsl = 0._dp
    zril = 0._dp
    zqsl = 0._dp
    zcfncl = 0._dp
    zchl = 0._dp
    zcfhl = 0._dp
    zbnl = 0._dp
    zbhnl = 0._dp
    zbml = 0._dp
    zbhl = 0._dp
    zustarl = 0._dp    
    ztkevl = 0._dp
    zustl = 0._dp
    zcptl = 0._dp
    zcpq = 0._dp
    zcfml = 0._dp

    !-----------------------------------------------------------------------
    !      2.2   surface humidity and virtual temperature
    !
    WHERE(mask)
       it(:)        = NINT(ptslm1(:) * 1000._dp)
!       IF (it(:) < jptlucu1 .OR. it(:) > jptlucu2) lookupoverflow = .TRUE.
       it(:)        = MAX(MIN(it(:), jptlucu2), jptlucu1)
       zes(:)       = tlucua(it(:)) / paphm1(:)
       zqsl(:)      = zes(:) / (1._dp - vtmpc1 * zes(:))


       it1(:)       = it(:) + 1
       it1(:)       = MAX(MIN(it1(:), jptlucu2), jptlucu1)
       zqst1(:)     = tlucua(it1(:)) / paphm1(:)
       zqst1(:)     = zqst1(:) / (1._dp - vtmpc1 * zqst1(:))
       zdqsl(:)     = (zqst1(:) - zqsl(:)) * 1000._dp
!!       pws(:)       = MIN(pws(:),pwsmx(:))
!!       zwstop(:)    = MIN(0.1, pwsmx(:))
!!       zwslev(:)    = pwsmx(:) - zwstop(:)
!!       WHERE (pws(:) .GT. zwslev(:) .AND. pws(:) .GT. zwslim(:) * pwsmx(:))
!!          zhum(:)   = 0.5_dp * (1.- COS((pws(:) - zwslev(:)) * api / zwstop(:)))
!!       ELSEWHERE
!!          zhum(:)   = 0._dp
!!       ENDWHERE
!!       zhsoil(:)    = pcvs(:) + (1._dp - pcvs(:)) * (pcvw(:) + (1._dp - pcvw(:)) * zhum(:))
!!       WHERE (pqm1(:) .GT. zqsl(:))
!!          zhsoil(:) = 1._dp
!!       ELSEWHERE
!!          zhsoil(:) = zhsoil(:)
!!       ENDWHERE
       ztesl(:)     = ptslm1(:) * (1.e5_dp/ paphm1(:))**zkappa
       ztvl(:)      = ztesl(:) * (1._dp + vtmpc1 * zhsoil(:) * zqsl(:))
    END WHERE

    !     ------------------------------------------------------------------
    !        3.     COMPUTATION OF THE EXCHANGE COEFFICIENTS.
    !        3.1       COMPUTATION OF BASIC QUANTITIES: WIND SHEAR,
    !                  RICHARDSON NUMBER,SQUARED MIXING LENGTHS, UNSTABLE
    !                  AND STABLE CASE COMMON FACTORS AND NEUTRAL CASE
    !                  COMMON PART OF THE DRAG COEFFICIENTS.
    !

    where(mask)
       zdu2(:)     = MAX(zepdu2, pum1(:)**2 + pvm1(:)**2)
       zqmitte(:)  = (pqm1(:) + zhsoil(:)*zqsl(:)) / 2._dp
       zqtmit(:)   = zx(:) * 0.5_dp + zqmitte(:)
       ztmit(:)    = (ptm1(:) + ptslm1(:)) / 2._dp


       zqsmit(:)   = (zqss(:) + zqsl(:)) / 2._dp
       ztemitte(:) = (zteta1(:) + ztesl(:)) / 2._dp
       zvirmitte(:)= (ztvir1(:) + ztvl(:)) / 2._dp
       zfux(:)     = zfaxe(:) / (zcpd * ztmit(:))
       zfox(:)     = zfaxe(:) / (zrd * ztmit(:))
       zmult1(:)   = 1._dp + vtmpc1 * zqtmit(:)
       zmult2(:)   = zfux(:) * zmult1(:) - zrvrd
       zmult3(:)   = zrdrv * zfox(:) * zqsmit(:) / (1._dp + zrdrv * zfox(:) * zfux(:) * zqsmit(:))
       zmult5(:)   = zmult1(:) - zmult2(:) * zmult3(:)
       zmult4(:)   = zfux(:) * zmult5(:) - 1._dp
       zdus1(:)    = paclc(:) * zmult5(:) + (1._dp - paclc(:)) * zmult1(:)
       zdus2(:)    = paclc(:) * zmult4(:) + (1._dp - paclc(:)) * vtmpc1
       zteldif(:)  = zlteta1(:) - ztesl(:)
       zqddif(:)   = (pqm1(:) + zx(:)) - (zqsl(:) * zhsoil(:))


       zbuoy(:)    = zdus1(:) * zteldif(:) + zdus2(:) * ztemitte(:) * zqddif(:)
       zril(:)     = pgeom1(:) * zbuoy(:) / (zvirmitte(:) * zdu2(:))

       zcdnl(:)    = (zkap / LOG(1._dp + pgeom1(:) / (g * paz0lm(:))))**2

  
       zchnl(:)    = (zkap / LOG(1._dp + pgeom1(:) / (g * paz0lh(:))))**2


       zucfl(:)    = 1._dp / (1._dp + zcons11 * zcdnl(:) * SQRT(ABS(zril(:)) * ( 1._dp &
            + pgeom1(:) / (g * paz0lm(:)))))
       zucfhl(:)   = 1._dp / (1._dp + zcons11 * zchnl(:) * SQRT(ABS(zril(:)) * (1._dp  &
            + pgeom1(:) / (g * paz0lh))))
       zscfl(:)    = SQRT(1._dp + ABS(zril(:)))
       zcons(:)    = zcons12 * paphm1(:) / &
            (ptm1(:) * (1._dp + vtmpc1 * pqm1(:) - zx(:)))
       zcfncl(:)   = zcons(:) * SQRT(zdu2(:)) * zcdnl(:)
       zcfnchl(:)  = zcons(:) * SQRT(zdu2(:)) * zchnl(:)
       zdthv(:)    = MAX(0._dp,(ztvl(:) - ztvir1(:)))
       zwstl(:)    = zdthv(:) * SQRT(zdu2(:)) / zvirmitte(:)

    END WHERE
    !----------------------------------------------------------------------------------------------
    !     3.2  DIMENSIONLESS HEAT TRANSFER COEFFICIENTS MULTIPLIED
    !          BY PRESSURE THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE
    !

    WHERE(mask)
       WHERE(zril(:).GT.0._dp)
          zcfml(:) = zcfncl(:) / (1._dp + zcons8 * zril(:) / zscfl(:))
          zcfhl(:) = zcfnchl(:) / (1._dp + zcons8 * zril(:) * zscfl(:))
          zchl(:)  = zcfhl(:) / zcfnchl(:) * zchnl(:)
       ELSEWHERE
          zcfml(:) = zcfncl(:) * (1._dp - zcons8 * zril(:) * zucfl(:))
          zcfhl(:) = zcfnchl(:) * (1._dp - zcons9 * zril(:) * zucfhl(:))
          zchl(:)  = zcfhl(:) / zcfnchl(:) * zchnl(:)
       END WHERE
    END WHERE

    !-----------------------------------------------------------------------------------------------
    !     interpolation functions for diagnostics
    !

   WHERE(mask)
       zbnl(:)     = zkap / SQRT(zcdnl(:))
       zbhnl(:)    = zkap / SQRT(zchnl(:))
       zbml(:)     = MAX(zepsec, SQRT(zcfml(:) * zcdnl(:) * zcons17 / zcfncl(:)))
       zbhl(:)     = MAX(zepsec, zchl(:) / zbml(:) * zcons17)
       zbml(:)     = 1._dp / zbml(:)
       zbhl(:)     = 1._dp / zbhl(:)
    END WHERE

    !-----------------------------------------------------------------------------------------------
    !*       3.4       COMPUTATION OF THE PBL EXTENSION.
    !

  WHERE(mask)
       WHERE(paz0lm(:) .GT. zepz0o)
          zcdn2m(:)= (zkap / LOG(1._dp+ pgeom1(:) / (g * zepz0o)))**2
       ELSEWHERE
          zcdn2m(:)= zcdnl(:)
       END WHERE
       zcdnr(:)    = zcdn2m(:) / zcdnl(:)
       WHERE(paz0lm(:) .GT. zepz0o .AND. zril(:) .LT. 0._dp)
          zcfm2m(:)= zcfncl(:) * zcdnr(:) * (1._dp - zcons8 * zril(:) &
               / (1._dp + zcons11 * zcdn2m(:) * SQRT(ABS(zril(:)) &
               * (1._dp + pgeom1(:) / (g * zepz0o)))))
       ELSEWHERE
          zcfm2m(:)= zcfml(:) * zcdnr(:)
       END WHERE
       zustl(:)    = zcfm2m(:) * SQRT(zdu2(:))
       zustarl(:)  = SQRT(zustl(:) * ptm1(:) &
            * (1._dp + vtmpc1 * pqm1(:) - zx(:)) &
            / (zcons12 * paphm1(:)))
    END WHERE

    !----------------------------------------------------------------------------------------------
    !      CONVECTIVE VELOCITY SCALE, MONIN-OBUKHOV LENGTH AND
    !      TKE BOUNDARY CONDITION (MAILHOT/BENOIT, 1982)

    WHERE(mask)
       WHERE(zwstl(:) .GT. zepsr)
          zconvs(:)= (zwstl(:) * zchl(:) * zghabl(:))**zcons6
          zmonob(:)= (zustarl(:)**3) / (zkap * g * zwstl(:) * zchl(:))
          zstabf(:)= (pgeom1(:) / (g * zmonob(:)))**(zcons6 * 2._dp)
          zstabf(:)= MIN(zustf * 3._dp, zstabf(:))
       ELSEWHERE
          zconvs(:)= 0._dp
          zstabf(:)= 0._dp
       END WHERE
       ztkevl(:)   = (zustf + zstabf(:)) * (zustarl(:)**2) + zwstf * (zconvs(:)**2)
    END WHERE

    WHERE(mask)
       zcptl(:)    = ptslm1(:) * cpd * (1._dp + vtmpc2 * &
            ( zcsat(:) * zqsl(:) + (1._dp-zcair(:)) * pqm1(:)))
       zcpq(:)     = zcptl(:) / ptslm1(:)
    END WHERE

  END SUBROUTINE PRECALC_LAND
  !------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE RICHTMEYR_land( kdim, mask, &
       jrow,                             &
       klev, klevp1, klevm1,             &
       paphm1, zcfh,                     &
       zebsh, zqdif,                     &
       ztdif,                            &
       zcfhl, zcair,                     &
       zcsat,                            &
       zetnl, zftnl,                     &
       zeqnl, zfqnl,                     &
!---wiso-code
       lwiso, kwiso,                     &
       ! Hilfsvariablen fuer qdif in mo_soil
       pfrl, pfrw, pfri,                 &
       zcfhw, zcfhi,                     &
       ocean_mask, ice_mask,             &
       ! Inputvariables: kinetic fractionaiton factor, diffusion coefficient of ???humidity???
       zwisoqdif,                        &
       ! Outputvariables: EN and FN coefficient of the Richtmyer-Morton-Scheme, water isotopes
       zwisoeqnl, zwisofqnl,             & 
       zwiso_nenner,                     &
       ! Koeffizient fuer Berechnung qdif (land only) analog zu mo_surface_boundary line 107
       zwiso_helpqdif                    &
!---wiso-code-end
       )

    USE mo_physc2,       ONLY: cvdifts
    USE mo_time_control, ONLY: lstart

    INTEGER,   INTENT(IN)  :: kdim, klev, klevp1, klevm1, jrow
!---wiso-code
    LOGICAL,   INTENT(in)  :: lwiso
    INTEGER,   INTENT(in)  :: kwiso
!---wiso-code-end
    LOGICAL,   INTENT(IN)  :: mask(kdim)
    REAL(dp),      INTENT(IN)  :: zcfhl(kdim), zcair(kdim), zcsat(kdim)
    REAL(dp),      INTENT(IN)  :: paphm1(kdim,klevp1), zcfh(kdim,klev), ztdif(kdim,klev)
    REAL(dp),      INTENT(IN)  :: zqdif(kdim,klev), zebsh(kdim,klev)
!---wiso-code
    REAL(dp),      INTENT(IN), OPTIONAL  :: pfrl(kdim), pfrw(kdim), pfri(kdim)
    REAL(dp),      INTENT(IN), OPTIONAL  :: zcfhw(kdim), zcfhi(kdim)
    LOGICAL,       INTENT(IN), OPTIONAL  :: ocean_mask(kdim), ice_mask(kdim)
    REAL(dp),      INTENT(in), OPTIONAL  :: zwisoqdif(kdim,klev,kwiso)
!---wiso-code-end
    REAL(dp),      INTENT(OUT) :: zetnl(kdim), zftnl(kdim), zeqnl(kdim), zfqnl(kdim)
!---wiso-code
    REAL(dp),      INTENT(OUT), OPTIONAL :: zwisoeqnl(kdim,kwiso), zwisofqnl(kdim,kwiso) 
    REAL(dp),      INTENT(OUT), OPTIONAL :: zwiso_nenner(kdim,kwiso)
    REAL(dp),      INTENT(OUT), OPTIONAL :: zwiso_helpqdif(kdim)
!---wiso-code-end

    ! Locals
    REAL(dp)   :: zdiscl(kdim), zdisql(kdim)
    REAL(dp)   :: zqdp(kdim), zfac(kdim)
    REAL(dp)   :: ztpfac1
!---wiso-code
    REAL(dp)   :: zgql(kdim), zgqw(kdim), zgqi(kdim)
    REAL(dp)   :: zero(kdim)
    REAL(dp)   :: ztpfac2
    INTEGER    :: jt, jl
    REAL(dp)   :: zdelta
!---wiso-code-end

    INTEGER :: i
    !------------------------
    ! CONSTANTS
    ztpfac1   = cvdifts
    zetnl = 0._dp
    zftnl = 0._dp
    zeqnl = 0._dp
    zfqnl = 0._dp

!---wiso-code
    IF (lwiso) THEN

    ztpfac2    = 1._dp / cvdifts
    zero(:)    = 0._dp
    zwiso_helpqdif(:)  = 0._dp
    DO jt=1,kwiso
       zwisoeqnl(:,jt)    = 0._dp
       zwisofqnl(:,jt)    = 0._dp
       zwiso_nenner(:,jt) = 0._dp
    END DO
    
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
       
       zdiscl(:)    = 1._dp / (1._dp + zfac(:) * (1._dp - zebsh(:,klevm1)) + zcfhl(:) * zqdp(:))
       zdisql(:)    = 1._dp / (1._dp + zfac(:) * (1._dp - zebsh(:,klevm1)) + zcair(:) * zcfhl(:) * zqdp(:))
       zetnl(:)     = zdiscl(:) * zcfhl(:) * zqdp(:)
       zftnl(:)     = zdiscl(:) * (ztdif(:,klev) + zfac(:) * ztdif(:,klevm1)) * ztpfac1
       zeqnl(:)     = zdisql(:) * zcsat(:) * zcfhl(:) * zqdp(:)
       zfqnl(:)     = zdisql(:) * (zqdif(:,klev) + zfac(:) * zqdif(:,klevm1)) * ztpfac1

    END WHERE

!---wiso-code
    IF (lwiso) THEN

   DO jt=1,kwiso
      DO jl=1,kdim
        IF (mask(jl)) THEN
        zwiso_nenner(jl,jt) = (1._dp + zfac(jl) * (1._dp - zebsh(jl,klevm1)))
        zwisoeqnl(jl,jt) = zcfhl(jl) * zqdp(jl)
        zwisofqnl(jl,jt) = (zwisoqdif(jl,klev,jt) + zfac(jl) * zwisoqdif(jl,klevm1,jt)) * ztpfac1
        END IF
      END DO
   IF (lstart) THEN
      zwisoeqnl(:,jt) = 0._dp
   END IF
   END DO

   ! Berechnung von zwiso_helpqdif
    zgql(:)    = pfrl(:) * MERGE(zcfhl(:),zero,mask)
    zgqw(:)    = pfrw(:) * MERGE(zcfhw(:),zero,ocean_mask)
    zgqi(:)    = pfri(:) * MERGE(zcfhi(:),zero,ice_mask)

    WHERE (pfrl(:) .LT. 1._dp)
       zgql(:) = MERGE(zgql(:),zero,mask) * MERGE(zcair(:),zero,mask)
    ELSE WHERE
       zgql(:) = MERGE(zgql(:),zero,mask)
    END WHERE

    zwiso_helpqdif(:)  = (zgql(:) + zgqw(:) + zgqi(:) ) / ztpfac2
    zwiso_helpqdif(:)  = zgql(:)/zwiso_helpqdif(:)
    
    END IF
!---wiso-code-end

  END SUBROUTINE RICHTMEYR_LAND
  !-----------------------------------------------------------------------------------
  !  
  SUBROUTINE UPDATE_LAND (                   &
         kdim , mask                         &
       , zetnl, zftnl                        &
       , zeqnl, zfqnl                        &
       , zcptlnew, zqslnew                   &
       , ztklevl, zqklevl                    &
       )

    INTEGER, INTENT(IN) :: kdim
    LOGICAL, INTENT(IN) :: mask(kdim)
    REAL(dp), INTENT(IN)    :: zetnl(kdim), zftnl(kdim), zeqnl(kdim), zfqnl(kdim)
    REAL(dp), INTENT(IN)    :: zcptlnew(kdim), zqslnew(kdim)
    REAL(dp), INTENT(OUT)   :: ztklevl(kdim), zqklevl(kdim)
    
    !*  CALCULATION OF SKLEV AND QKLEV USING THE NEW SURFACE VALUES
    !*  ZSNEW AND ZQSNEW WHICH WERE CALCULATED IN SUBROUTINE SURFTEMP

    ztklevl = 0._dp
    zqklevl = 0._dp

    WHERE(MASK)
       
       ztklevl(:) = zetnl(:) * zcptlnew(:) + zftnl(:)
       zqklevl(:) = zeqnl(:) * zqslnew(:)  + zfqnl(:)
       
    END WHERE

  END SUBROUTINE UPDATE_LAND

  !-----------------------------------------------------------------------------------
  !  
!---wiso-code
  SUBROUTINE UPDATE_LAND_WISO(                         &
       kwiso                                           &
       , kdim , mask                                   &
       , zwisoeqnl, zwisofqnl                          &
       , zqslnew, zwisoqslnew                          & 
       , zwisocsat, zwisocair                          &
       , zwisocsat_fra, zwisocair_fra                  &
       , zqdif                                         &
       , zwisoqklevl                                   &
       )

    USE mo_physc2,       ONLY: cvdifts
    USE mo_time_control, ONLY: lstart

    INTEGER, INTENT(IN) :: kdim
    INTEGER, INTENT(IN) :: kwiso
    LOGICAL, INTENT(IN) :: mask(kdim)
    REAL(dp), INTENT(IN)    :: zqslnew(kdim)
    REAL(dp), INTENT(IN)    :: zwisoeqnl(kdim,kwiso), zwisofqnl(kdim,kwiso)
    REAL(dp), INTENT(IN)    :: zwisocsat(kdim,kwiso), zwisocair(kdim,kwiso)
    REAL(dp), INTENT(IN)    :: zwisocsat_fra(kdim,kwiso), zwisocair_fra(kdim,kwiso)
    REAL(dp), INTENT(IN)    :: zwisoqslnew(kdim,kwiso)
    REAL(dp), INTENT(IN)    :: zqdif(kdim)
    REAL(dp), INTENT(OUT)   :: zwisoqklevl(kdim,kwiso)

    ! Locals:
    INTEGER    :: jt, jl
    REAL(dp)   :: ztpfac1, ztpfac2

    ztpfac1   = cvdifts
    
    !*  CALCULATION OF SKLEV AND QKLEV USING THE NEW SURFACE VALUES
    !*  ZSNEW AND ZQSNEW WHICH WERE CALCULATED IN SUBROUTINE SURFTEMP

    zwisoqklevl = 0._dp
    ztpfac1 = cvdifts
    ztpfac2 = 1._dp/ztpfac1

    IF (lstart) THEN
       DO jt=1,kwiso
          WHERE(MASK)
             zwisoqklevl(:,jt) = zwisofqnl(:,jt)            
          END WHERE
       END DO
    ELSE
       DO jt=1,kwiso
          WHERE(MASK)
             zwisoqklevl(:,jt) = ( zwisoeqnl(:,jt) * ( zwisocsat(:,jt) * zqslnew(:) &
                                                      - zwisocair(:,jt) * zqdif(:) * ztpfac1 ) ) &
                                 + ( zwisoeqnl(:,jt) * zwisocsat_fra(:,jt) * zwisoqslnew(:,jt) ) &
                                 + zwisofqnl(:,jt)          
          END WHERE
       END DO
    END IF

  END SUBROUTINE UPDATE_LAND_WISO
!---wiso-code-end

  !-----------------------------------------------------------------------------------
  !
  SUBROUTINE postproc_land( &
    kdim      , mask        &
    , pgeom1  , zbnl        &
    , zbml    , pum1        &
    , pvm1    , zril        &
    , zspeedl , zbhnl       &
    , zbhl    , zcptl       &
    , zcptgz  , pqm1        &
    , zt2l    , ptm1        &
    , papm1   , paphm1      &
    , klevp1  , zx          &
    , zdew2l  , zu10l       &
    , zv10l                 &
    )

    USE mo_constants,    ONLY: g, vtmpc2, cpd, vtmpc1, tmelt, c3les, c4les, c3ies, c4ies &
         , rd, c2es

    INTEGER,  INTENT(in)   :: kdim, klevp1
    LOGICAL,  INTENT(in)   :: mask(kdim)
    REAL(dp),     INTENT(in)   :: pgeom1(kdim), zbnl(kdim), zbml(kdim), zril(kdim)
    REAL(dp),     INTENT(in)   :: pum1(kdim), pvm1(kdim)
    REAL(dp),     INTENT(out)  :: zspeedl(kdim)
    REAL(dp),     INTENT(in)   :: zbhnl(kdim), zbhl(kdim), zcptl(kdim)
    REAL(dp),     INTENT(in)   :: zcptgz(kdim), pqm1(kdim)
    REAL(dp),     INTENT(out)  :: zt2l(kdim)
    REAL(dp),     INTENT(in)   :: ptm1(kdim), papm1(kdim), paphm1(kdim,klevp1), zx(kdim)
    REAL(dp),     INTENT(out)  :: zdew2l(kdim), zu10l(kdim), zv10l(kdim)

    REAL(dp)     :: zhtq, zephum, zhuv
    REAL(dp)     :: zrat(kdim), zcbn(kdim), zcbs(kdim), zmerge1(kdim), zred(kdim)
    REAL(dp)     :: zmerge(kdim), zcbu(kdim)
    INTEGER  :: it(kdim)
    REAL(dp)     :: zqs1(kdim), zrh2m(kdim), zcvm3(kdim), zcvm4(kdim), zqs2(kdim)
    REAL(dp)     :: zaph2m(kdim), zq2m(kdim), zfrac(kdim), zh2m(kdim)

    integer :: i

    zspeedl = 0._dp
    zt2l = 0._dp
    zdew2l = 0._dp
    zu10l = 0._dp
    zv10l = 0._dp

    zhtq          = 2._dp * g
    zephum        = 5.e-2_dp
    zhuv          = 10._dp * g
    !--------------------------------------------------------------------------------
    ! 10M WIND COMPONENTS, MAX 10M WINDSPEED
    WHERE(mask)
       zrat(:)        = zhuv / pgeom1(:)
       zcbn(:)        = LOG(1._dp+ (EXP (zbnl(:)) - 1._dp) * zrat(:))
       zcbs(:)        = -(zbnl(:) - zbml(:)) * zrat(:)
       zcbu(:)        = -LOG(1._dp +(EXP (zbnl(:) - zbml(:)) - 1._dp) * zrat(:))
       WHERE(zril(:).GT.0._dp)
          zmerge1(:)  = zcbs(:)
       ELSEWHERE
          zmerge1(:)  = zcbu(:)
       END WHERE
       zred(:)        = (zcbn(:) + zmerge1(:)) / zbml(:)
       zu10l(:)       = zred(:) * pum1(:)
       zv10l(:)       = zred(:) * pvm1(:)
       zspeedl(:)     = SQRT(zu10l(:)**2 + zv10l(:)**2)
    END WHERE
    !--------------------------------------------------------------------------------
    ! 2M temperature
    WHERE(mask)    
       zrat(:)        = zhtq / pgeom1(:)
       zcbn(:)        = LOG(1._dp+ (EXP (zbnl(:)) - 1._dp) * zrat(:))

       zcbn(:) = LOG( 1._dp + (EXP ( zbhnl(:)) - 1._dp) * zrat(:))
       zcbs(:) = -(zbhnl(:) - zbhl(:)) * zrat(:)
       zcbu(:) = -LOG(1._dp + (EXP ( zbhnl(:) - zbhl(:)) -1._dp) * zrat(:))
       WHERE(zril(:) .GT. 0._dp)
          zmerge(:)   = zcbs(:)
       ELSEWHERE
          zmerge(:)   = zcbu(:)
       END WHERE
       zred(:)    = (zcbn(:) + zmerge(:)) / zbhl(:)
       zh2m(:)    = zcptl(:) + zred(:) * ( zcptgz(:) - zcptl(:))
       zt2l(:)    = (zh2m(:) - zhtq) / (cpd * (1._dp + vtmpc2 * pqm1(:)))
    END WHERE
    !--------------------------------------------------------------------------------
    ! 2M dew point
    WHERE(mask)    
       it(:)          = NINT(ptm1(:) * 1000._dp)
       it(:)          = MAX(MIN(it,jptlucu2),jptlucu1)
       zqs1(:)        = tlucua(it(:)) / papm1(:)
       zqs1(:)        = zqs1(:) / (1._dp - vtmpc1 * zqs1(:))
       zrh2m(:)       = MAX(zephum, pqm1(:) / zqs1(:))
       WHERE(zt2l(:) .GT. tmelt)
          zcvm3(:)    = c3les
          zcvm4(:)    = c4les
       ELSEWHERE
          zcvm3(:)    = c3ies
          zcvm4(:)    = c4ies
       END WHERE
       zaph2m(:)      = paphm1(:,klevp1) * &
            (1._dp - zhtq / (rd * zt2l(:) * (1._dp + vtmpc1 * pqm1(:) - zx(:))))
    END WHERE
!!$    do i=1,SIZE(ptm1)
!!$       if(mask(i)) then
!!$          it(i) = NINT(zt2l(i) * 1000._dp)
!!$          it(i) = MAX(MIN(it(i),jptlucu2),jptlucu1)
!!$          zqs2(i)        = tlucua(it(i)) / zaph2m(i)
!!$          zqs2(i)        = zqs2(i) / (1._dp - vtmpc1 * zqs2(i))
!!$          zq2m(i)        = zrh2m(i) * zqs2(i)
!!$          print*,zt2l(i), it(i), zaph2m(i), zrh2m(i), zq2m(i),c2es,vtmpc1,zcvm3(i)
!!$          zfrac(i)       = LOG(zaph2m(i) * zq2m(i) / (c2es * (1._dp + vtmpc1 * zq2m(i)))) / zcvm3(i)
!!$          zdew2l(i)      = MIN(zt2l(i), (tmelt - zfrac(i) * zcvm4(i)) / (1._dp - zfrac(i)))
!!$       end if
!!$    end do
    WHERE(mask)
       it(:)          = NINT(zt2l(:) * 1000._dp)
       it(:)          = MAX( MIN( it(:), jptlucu2), jptlucu1)
       zqs2(:)        = tlucua(it(:)) / zaph2m(:)
       zqs2(:)        = zqs2(:) / (1._dp - vtmpc1 * zqs2(:))
       zq2m(:)        = zrh2m(:) * zqs2(:)
       zfrac(:)       = LOG(zaph2m(:) * zq2m(:) / (c2es * (1._dp+ vtmpc1 * zq2m(:)))) / zcvm3(:)
       zdew2l(:)      = MIN(zt2l(:), (tmelt - zfrac(:) * zcvm4(:)) / (1._dp - zfrac(:)))
    END WHERE

  END SUBROUTINE POSTPROC_LAND

  SUBROUTINE land_rad(        &
     kdim          , mask         &
     ,ztrdown      , ptslm1       &
     ,pi0          , ptrsol       &
     ,palsoi       , palbedo      &
     ,psofli       , ptrfli       &
     ,ptslnew      , zteffl4      &
     )

  USE mo_radiation,         ONLY: cemiss
  USE mo_constants,         ONLY: stbo

  INTEGER, INTENT(in)    :: kdim
  LOGICAL, INTENT(in)    :: mask(kdim)
  REAL(dp),    INTENT(in)    :: ptslm1(kdim), pi0(kdim), ptrsol(kdim), ptslnew(kdim)
  REAL(dp),    INTENT(in)    :: palsoi(kdim), ztrdown(kdim), palbedo(kdim)
  REAL(dp),    INTENT(out)   :: psofli(kdim), ptrfli(kdim), zteffl4(kdim)


  REAL(dp) :: zflxs(kdim)

  psofli  = 0._dp
  ptrfli  = 0._dp
  zteffl4 = 0._dp
  WHERE(mask)    
     zteffl4(:) = ptslm1(:)**3*(4._dp*ptslnew(:)-3._dp*ptslm1(:))
   
     ptrfli(:)  = ztrdown(:) - cemiss * stbo * zteffl4(:)
     zflxs(:)   = pi0(:) * ptrsol(:)     
     psofli(:)  = (1._dp - palsoi(:)) * zflxs(:) / (1._dp - palbedo(:))
      
  END WHERE

END SUBROUTINE land_rad

END MODULE MO_SURFACE_LAND
