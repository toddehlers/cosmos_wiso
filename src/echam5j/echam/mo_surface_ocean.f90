MODULE mo_surface_ocean
  
  USE mo_convect_tables,   ONLY: tlucua, jptlucu1, jptlucu2, &
       lookuperror, lookupoverflow
  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  REAL(dp), PARAMETER :: calbsea = 0.07_dp    ! sea albedo
  
CONTAINS
 
  !-------------------------------------------------------------------------------------------------
  ! Set the sea albedo to 0.07 
  SUBROUTINE update_albedo_ocean(mask, p_albedo)
      
    LOGICAL, INTENT(in)  :: mask(:)
    REAL(dp),    INTENT(out) :: p_albedo(:)
    
    p_albedo = MERGE(calbsea, 0.0_dp, mask)

  END SUBROUTINE update_albedo_ocean
  
  !-------------------------------------------------------------------------------------------------
  
  SUBROUTINE update_z0_ocean(kdim, klevp1, &
       mask, zcfmw, zudif, zvdif, ptm1, pqm1, zx,&
       paphm1, paz0w&
       )
    
    USE mo_constants,    ONLY: g, rd, vtmpc1
    USE mo_physc2,       ONLY: cchar
    USE mo_time_control, ONLY: time_step_len

    INTEGER, INTENT(in)     :: kdim, klevp1
    LOGICAL, INTENT(in)     :: mask(kdim)
    REAL(dp), INTENT(in)        :: zcfmw(kdim), zudif(kdim), zvdif(kdim), ptm1(kdim)
    REAL(dp), INTENT(in)        :: pqm1(kdim), zx(kdim), paphm1(kdim, klevp1)
    REAL(dp), INTENT(out)       :: paz0w(kdim)

    REAL(dp) :: zepzzo, ztmst, zcons14

    ztmst  = time_step_len
    zcons14 = cchar * rd / (g**2 * ztmst)
    zepzzo=1.5e-05_dp

    paz0w = 0.0005_dp

    WHERE(mask)
       
       paz0w(:)  = MAX(zepzzo, zcons14 * zcfmw(:) &
            * SQRT(zudif(:)**2 + zvdif(:)**2) * ptm1(:) &
            * (1._dp + vtmpc1 * pqm1(:) - zx(:)) / paphm1(:,klevp1))
    END WHERE

  END SUBROUTINE UPDATE_Z0_OCEAN
  
  !-------------------------------------------------------------------------------------------------
  
  SUBROUTINE update_stress_ocean(kdim, mask, zcfmw, zudif, zvdif &
       , pustrw, pvstrw)

    USE mo_time_control, ONLY: time_step_len
    USE mo_constants,    ONLY: g
    
    INTEGER, INTENT(in)    :: kdim
    LOGICAL, INTENT(in)    :: mask(kdim)
    REAL(dp), INTENT(in)       :: zcfmw(kdim)
    REAL(dp), INTENT(in)       :: zudif(kdim), zvdif(kdim)
    REAL(dp), INTENT(out)      :: pustrw(kdim), pvstrw(kdim)
    
    REAL(dp) :: ztmst, zcons15
    
    ztmst   = time_step_len
    zcons15 = 1._dp / (g * ztmst)

    pustrw = 0._dp
    pvstrw = 0._dp

    Where (mask)
       pustrw(:) = zcons15 * zcfmw(:) * zudif(:)
       pvstrw(:) = zcons15 * zcfmw(:) * zvdif(:)
    END Where
    
  END SUBROUTINE UPDATE_STRESS_OCEAN
  
  !-----------------------------------------------------------------------------------------------------
  
  SUBROUTINE precalc_ocean( kdim   &
       , mask &
       , ptsw, paphm1 & 
       , pqm1, zx &
       , ptm1, zqss &
       , zteta1, ztvir1 &
       , zfaxe ,paclc &
       , zlteta1, pum1 &
       , pvm1                                                          &
       , pocu, pocv                                                    &
       , paz0w &
       , pgeom1, zghabl &
       , zqsw, zcptw &
       , zriw, zcfhw &
       , zchw, zbnw &
       , zbmw, zbhw &
       , zustarw, ztkevw &
       , zcfmw , zustw &
!---wiso-code
       , lwiso, kwiso &
       , pwisosw_d &
       , zwisoqsw &
!---wiso-code-end
       )

    USE mo_constants,    ONLY: rd, cpd, vtmpc1, g, vtmpc2
    USE mo_physc2,       ONLY: ckap, cc, cb, cvdifts, cfreec, cgam
    USE mo_time_control, ONLY: time_step_len
!---wiso-code
    USE mo_wiso,         ONLY: talphal1, talphal2, talphal3, talphal4, &
                               tsatbase,                               &
                               toce, nwiso
!---wiso-code-end

    INTEGER, INTENT(in)    :: kdim
!---wiso-code
    LOGICAL, INTENT(in)    :: lwiso
    INTEGER, INTENT(in)    :: kwiso
!---wiso-code-end
    LOGICAL, INTENT(in)    :: mask(kdim)
    REAL(dp),    INTENT(IN)    :: ptsw(kdim), paphm1(kdim), pqm1(kdim), zx(kdim)
    REAL(dp),    INTENT(in)    :: ptm1(kdim), zqss(kdim), zteta1(kdim)
    REAL(dp),    INTENT(IN)    :: ztvir1(kdim), zfaxe(kdim), paclc(kdim), zlteta1(kdim)
    REAL(dp),    INTENT(in)    :: pum1(kdim), pvm1(kdim), paz0w(kdim)
    REAL(dp),    INTENT(in)    :: pocu(kdim), pocv(kdim)
    REAL(dp),    INTENT(IN)    :: pgeom1(kdim), zghabl(kdim)
    
    REAL(dp),    INTENT(out)   :: zqsw(kdim), zcptw(kdim), zriw(kdim)
    REAL(dp),    INTENT(out)   :: zcfhw(kdim), zchw(kdim)
    REAL(dp),    INTENT(out)   :: zbnw(kdim), zbmw(kdim), zbhw(kdim), zustarw(kdim)
    REAL(dp),    INTENT(out)   :: ztkevw(kdim), zcfmw(kdim), zustw(kdim)
!---wiso-code
    REAL(dp),    INTENT(in), OPTIONAL    :: pwisosw_d(kdim,kwiso)
    REAL(dp),    INTENT(out), OPTIONAL   :: zwisoqsw(kdim,kwiso)
!---wiso-code-end

    ! Achtung, (klevp1): paphm1 
    ! ACHTUNG, (klev): zteta1, zvir1, pqm1, zfaxe, paclc, zlteta1, pgeom1, pum1, pvm1
    
    ! Local Variables
    INTEGER  :: it(kdim)
    REAL(dp)     :: zes(kdim), zqmitte(kdim), zqtmit(kdim), ztmit(kdim), zqsmit(kdim), zvirmitte(kdim)
    REAL(dp)     :: ztemitte(kdim), zfux(kdim), zfox(kdim), zmult1(kdim), zbuoy(kdim), zdu2oc(kdim)
    REAL(dp)     :: zmult2(kdim), zmult3(kdim), zmult4(kdim), zmult5(kdim), zdus1(kdim), zdus2(kdim)
    REAL(dp)     :: zteldif(kdim), z0h(kdim), zalo(kdim), zaloh(kdim), zcons(kdim), zdthv(kdim) 
    REAL(dp)     :: zucfw(kdim), zscfw(kdim), zcfnchw(kdim), zchnw(kdim), zcr(kdim)
    REAL(dp)     :: zcdn2m(kdim), zcdnr(kdim), zcfm2m(kdim), ztesw(kdim), ztvw(kdim), zcdnw(kdim)
    REAL(dp)     :: zcfncw(kdim), zwstw(kdim), zconvs(kdim), zmonob(kdim), zstabf(kdim)
    
    REAL(dp)     :: zkappa, zrvrd, zepdu2, zkap, zcons8, zcons11, zcons12, zepsec, zcons17, zepz0o
    REAL(dp)     :: zrdrv, zepsr, zcons6, zustf
    REAL(dp)     :: zsmn, zshn, zm1, zm2, zm4, zh1, zh2, zwstf

!---wiso-code
    INTEGER      :: jt, jl
    REAL(dp)     :: ztsw_tmp(kdim), zwisofracw(kdim,kwiso), zocean(kdim,kwiso)
!---wiso-code-end
    
    ! Constants
    zepsec      = 1.e-2_dp
    zkappa      = rd / cpd
    zrvrd       = vtmpc1 + 1._dp
    zepdu2      = 1.0_dp
    zkap        = ckap
    zcons8      = 2._dp * cb
    zepz0o      = 2._dp
    zcons11     = 3._dp * cb * cc
    zcons12     = cvdifts * time_step_len * g / rd
    zcons17     = 1._dp / zkap**2
    zrdrv       = 1._dp / zrvrd
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
    zepsr         = 1.e-10_dp
 
    zqsw =  0._dp
    zcptw =  0._dp
    zriw =  0._dp
    zcfhw =  0._dp
    zchw =  0._dp
    zbnw =  0._dp
    zbmw =  0._dp
    zbhw =  0._dp
    zustarw =  0._dp
    ztkevw =  0._dp
    zcfmw =  0._dp
    zustw =  0._dp

!---wiso-code
    IF (lwiso) THEN

    zwisoqsw =  0._dp
    
    END IF
!---wiso-code-end

    !-----------------------------------------------------------------------
    !      2.2   surface humidity and virtual temperature
           
    WHERE(mask)
       it(:)    = NINT(ptsw(:) * 1000._dp)
!       WHERE (it(:) .LT. jptlucu1 .OR. it (:) .GT. jptlucu2) 
!          lookupoverflow = .TRUE.
!       ENDWHERE
       it(:)    = MAX(MIN(it(:), jptlucu2), jptlucu1)
       zes(:)   = tlucua(it(:)) / paphm1(:)
       zqsw(:)  = zes(:) / (1._dp-vtmpc1 * zes(:))
       zcptw(:) = ptsw(:) * cpd * (1._dp + vtmpc2 * zqsw(:))
       ztesw(:) = ptsw(:) * (1.e5_dp / paphm1(:))**zkappa
       ztvw(:)  = ztesw(:) * (1._dp + vtmpc1 * zqsw(:))
    ENDWHERE
    
!---wiso-code
    IF (lwiso) THEN

    ! fractionation over water
    DO jt=1,kwiso        ! loop over types of water isotopes
        WHERE(mask)      ! loop over longitudes includes mask
            ztsw_tmp(:)=ptsw(:)
            zwisofracw(:,jt)=(exp(talphal1(jt)/(ztsw_tmp(:)**2)+talphal2(jt)/ztsw_tmp(:)+talphal3(jt)))**talphal4(jt)
            ! "corrections" of ocean surface isotope values by prescriped delta values
            zocean(:,jt)=toce(jt)*(1._dp+(pwisosw_d(:,jt)/1000._dp))
            zwisofracw(:,jt)=zocean(:,jt)/zwisofracw(:,jt)
            !calculation of tracer humidity saturation
            zwisoqsw(:,jt)=zwisofracw(:,jt)*zqsw(:)
        ENDWHERE
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
    !
    ! correction for water and ice points
    !
       zdu2oc(:)    = MAX(zepdu2,(pum1(:)-pocu(:))**2                  &
                                +(pvm1(:)-pocv(:))**2)
       zcons(:)     = zcons12 * paphm1(:) / (ptm1(:) * (1._dp + vtmpc1 * pqm1(:) - zx(:)))
       zqmitte(:)   = (pqm1(:) + zqsw(:)) / 2._dp
       zqtmit(:)    = zx(:) * 0.5_dp + zqmitte(:)
       ztmit(:)     = (ptm1(:) + ptsw(:)) / 2._dp
       zqsmit(:)    = (zqss(:) + zqsw(:)) / 2._dp
       ztemitte(:)  = (zteta1(:) + ztesw(:)) / 2._dp
       zvirmitte(:) = (ztvir1(:) + ztvw(:)) / 2._dp
       zfux(:)      = zfaxe(:) / (cpd * ztmit(:))
       zfox(:)      = zfaxe(:) / (rd * ztmit(:))
       zmult1(:)    = 1._dp + vtmpc1 * zqtmit(:)
       zmult2(:)    = zfux(:) * zmult1 - zrvrd
       zmult3(:)    = zrdrv * zfox(:) * zqsmit(:) / (1._dp + zrdrv * zfox(:) * zfux(:) * zqsmit(:))
       zmult5(:)    = zmult1(:) - zmult2(:) * zmult3(:)
       zmult4(:)    = zfux(:) * zmult5(:) - 1._dp
       zdus1(:)     = paclc(:) * zmult5(:) + ( 1._dp - paclc(:)) * zmult1(:)
       zdus2(:)     = paclc(:) * zmult4(:) + ( 1._dp - paclc(:)) * vtmpc1
       zteldif(:)   = zlteta1(:) - ztesw(:)
       zbuoy(:)     = zdus1(:) * zteldif(:) + zdus2(:) * ztemitte(:) * ((pqm1(:) + zx(:)) - zqsw(:))
       zriw(:)      = pgeom1(:) * zbuoy(:) / (zvirmitte(:) * zdu2oc(:))
       z0h(:)       = paz0w(:) * EXP(2._dp-86.276_dp * paz0w(:)**0.375_dp)
       zalo(:)      = LOG(1._dp + pgeom1(:) / (g * paz0w(:)))
       zaloh(:)     = LOG(1._dp + pgeom1(:) / (g * z0h(:)))
       zcdnw(:)     = (zkap / zalo(:))**2
       zchnw(:)     = zkap**2 / (zalo(:) * zaloh(:))
       zucfw(:)     = 1._dp / (1._dp + zcons11 * zcdnw(:) * SQRT(ABS(zriw(:)) &
            * (1._dp + pgeom1(:) / (g * paz0w(:)))))
       zscfw(:)     = SQRT(1._dp + ABS(zriw(:)))
       zcfncw(:)    = zcons(:) * SQRT(zdu2oc(:)) * zcdnw(:)
       zcfnchw(:)   = zcons(:) * SQRT(zdu2oc(:)) * zchnw(:)
       zdthv(:)     = MAX(0._dp,(ztvw(:) - ztvir1(:)))
       zwstw(:)     = zdthv(:) * SQRT(zdu2oc(:)) / zvirmitte(:)
       zcr(:)       = (cfreec / (zchnw(:) * SQRT(zdu2oc(:)))) * ABS(zbuoy(:))**(1._dp/3._dp)
    ENDWHERE
    
    !----------------------------------------------------------------------------------------------
    !     3.2  DIMENSIONLESS HEAT TRANSFER COEFFICIENTS MULTIPLIED
    !          BY PRESSURE THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE
    !
    WHERE(mask) 
       WHERE(zriw(:).GT.0._dp)
          zcfmw(:)  = zcfncw(:) / (1._dp + zcons8 * zriw(:) / zscfw(:))
          zcfhw(:)  = zcfnchw(:) / (1._dp + zcons8 * zriw(:) * zscfw(:))
          zchw(:)   = zcfhw(:) / zcfnchw(:) * zchnw(:)
       ELSEWHERE
          zcfmw(:)  = zcfncw(:) *(1._dp - zcons8 * zriw(:) * zucfw(:))
          zcfhw(:)  = zcfnchw(:) * (1._dp + zcr(:)**cgam)**(1._dp/cgam)
          zchw(:)   = zcfhw(:) / zcfnchw(:) * zchnw(:)
       ENDWHERE
    ENDWHERE
    
    !-----------------------------------------------------------------------------------------------
    !     interpolation functions for diagnostics
    !
    WHERE(mask)
       zbnw(:)       = zkap / SQRT(zcdnw(:))
       zbmw(:)       = MAX(zepsec, SQRT(zcfmw(:) * zcdnw(:) * zcons17 / zcfncw(:)))
       zbhw(:)       = MAX(zepsec, zchw(:) / zbmw(:) * zcons17)
       zbmw(:)       = 1._dp / zbmw(:)
       zbhw(:)       = 1._dp / zbhw(:)
    ENDWHERE
    
    !-----------------------------------------------------------------------------------------------
    !*       3.4       COMPUTATION OF THE PBL EXTENSION.
    !
    WHERE(mask)
       WHERE(paz0w(:).GT.zepz0o)
          zcdn2m(:)   = (zkap / LOG(1._dp + pgeom1(:) / (g * zepz0o)))**2
       ELSEWHERE
          zcdn2m(:)   = zcdnw(:)
       ENDWHERE
       
       zcdnr(:)       = zcdn2m(:) / zcdnw(:)
       
       WHERE(paz0w(:).GT.zepz0o.AND.zriw(:).LT.0._dp)
          zcfm2m(:)   = zcfncw(:) * zcdnr(:) * (1._dp - zcons8 *zriw(:)) &
               / (1._dp+ zcons11 * zcdn2m(:) * SQRT(ABS(zriw(:)) &
               * (1._dp+ pgeom1(:) / (g * zepz0o))))
       ELSEWHERE
          zcfm2m(:)   = zcfmw(:)*zcdnr(:)
       ENDWHERE
       zustw(:)    = zcfm2m(:) * SQRT(zdu2oc(:))
       zustarw(:)  = SQRT(zustw(:) * ptm1(:) &
            * (1._dp + vtmpc1 * pqm1(:) - zx(:)) &
            / (zcons12 * paphm1(:)))
    ENDWHERE
    !----------------------------------------------------------------------------------------------
    !      CONVECTIVE VELOCITY SCALE, MONIN-OBUKHOV LENGTH AND
    !      TKE BOUNDARY CONDITION (MAILHOT/BENOIT, 1982)
    WHERE(mask)
       WHERE(zwstw(:) .GT. zepsr)
          zconvs(:) = (zwstw(:) * zchw(:) * zghabl(:))**zcons6
          zmonob(:) = (zustarw(:)**3) / (zkap * g * zwstw(:) * zchw(:))
          zstabf(:) = (pgeom1(:) / (g * zmonob(:)))**(zcons6 * 2._dp)
          zstabf(:) = MIN(zustf*3._dp, zstabf(:))
       ELSEWHERE
          zconvs=0._dp
          zstabf=0._dp
       ENDWHERE
       ztkevw(:)    = (zustf + zstabf(:)) * (zustarw(:)**2) + zwstf * (zconvs(:)**2)
    ENDWHERE
    
    
  END SUBROUTINE precalc_ocean
  
  !------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE RICHTMEYR_ocean(kdim, mask, &
       klev, klevp1, klevm1, &
       paphm1, zcfh, &
       zebsh, zqdif, &
       ztdif, zcfhw, &
       zetnw, zftnw, &
       zeqnw, zfqnw, &
!---wiso-code
       lwiso, kwiso, &
       ! Inputvariables: kinetic fractionaiton factor, diffusion coefficient of humidity
       zwisokinw, zwisoqdif, &
       ! Outputvariables: EN and FN coefficient of the Richtmyer-Morton-Scheme, water isotopes
       zwisoeqnw, zwisofqnw &
!---wiso-code-end
       )

    USE mo_physc2,       ONLY: cvdifts

    INTEGER, INTENT(IN)   :: kdim, klev, klevp1, klevm1
!---wiso-code
    LOGICAL, INTENT(in)   :: lwiso
    INTEGER, INTENT(in)   :: kwiso
!---wiso-code-end
    LOGICAL, INTENT(IN)   :: mask(kdim)
    REAL(dp), INTENT(IN)  :: zcfhw(kdim)
    REAL(dp), INTENT(IN)  :: paphm1(kdim,klevp1), zcfh(kdim,klev), ztdif(kdim,klev)
    REAL(dp), INTENT(IN)  :: zqdif(kdim,klev), zebsh(kdim,klev)
!---wiso-code
    REAL(dp), INTENT(in), OPTIONAL  :: zwisokinw(kdim,kwiso), zwisoqdif(kdim,klev,kwiso)
!---wiso-code-end
    REAL(dp), INTENT(OUT) :: zetnw(kdim), zftnw(kdim), zeqnw(kdim), zfqnw(kdim)
    ! klevm1   : zebsh
    ! klevm1   : zcfh
!---wiso-code
    REAL(dp), INTENT(OUT), OPTIONAL :: zwisoeqnw(kdim,kwiso), zwisofqnw(kdim,kwiso)
!---wiso-code-end

    ! Local variables
    REAL(dp)   :: zdiscw(kdim), zdisqw(kdim)
    REAL(dp)   :: zqdp(kdim), zfac(kdim)
    REAL(dp)   :: ztpfac1
!---wiso-code
    INTEGER                :: jt, jl
    REAL(dp), DIMENSION(2) :: zwisodisqw(kdim,kwiso)
!---wiso-code-end

    !------------------------
    ! CONSTANTS
    ztpfac1   = cvdifts
    
    zetnw = 0._dp
    zftnw = 0._dp
    zeqnw = 0._dp
    zfqnw = 0._dp

!---wiso-code
    IF (lwiso) THEN

    zwisoeqnw = 0._dp
    zwisofqnw = 0._dp
    
    END IF
!---wiso-code-end

    WHERE(mask)
       zqdp(:)       = 1._dp / (paphm1(:,klevp1) - paphm1(:,klev))
       zfac(:)       = zcfh(:,klevm1) * zqdp(:)

       !*  CALCULATION OF THE EN AND FN COEFFICIENTS OF THE RICHTMYER-
       !*  MORTON-SCHEME CONCERNING THE EQUATION:
       !
       !*  XN = EN * XS + FN
       !
       !*  WITH XN = S_ATM  OR  XN = QATM : ATM. VALUE OF S OR Q
       !*  AND  XS = SSURF  OR  XS = QSAT : SURFACE VALUE OF S OR SAT. SPEC.
       !*                                   HUM. AT THE SURFACE
       !
       zdiscw(:)    = 1._dp / (1._dp + zfac(:) * (1._dp - zebsh(:,klevm1)) + zcfhw(:) * zqdp(:))
       zdisqw(:)    = zdiscw(:)
       zetnw(:)     = zdiscw(:) * zcfhw(:) * zqdp(:)
       zftnw(:)     = zdiscw(:) * (ztdif(:,klev) + zfac(:) * ztdif(:,klevm1)) * ztpfac1
       zeqnw(:)     = zdisqw(:) * zcfhw(:) * zqdp(:)
       zfqnw(:)     = zdisqw(:) * (zqdif(:,klev) + zfac(:) * zqdif(:,klevm1)) * ztpfac1
       
    ENDWHERE

!---wiso-code
    IF (lwiso) THEN

    DO jt=1,kwiso
        WHERE(mask)
            zwisodisqw(:,jt)    = 1._dp / (1._dp + zfac(:) * (1._dp - zebsh(:,klevm1)) + zwisokinw(:,jt) * zcfhw(:) * zqdp(:))
            zwisoeqnw(:,jt)     = zwisodisqw(:,jt) * zcfhw(:) * zqdp(:)
            zwisofqnw(:,jt)     = zwisodisqw(:,jt) * (zwisoqdif(:,klev,jt) + zfac(:) * zwisoqdif(:,klevm1,jt)) * ztpfac1
        END WHERE
    END DO
    
    END IF
!---wiso-code-end

  END SUBROUTINE RICHTMEYR_OCEAN
  !---------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE UPDATE_OCEAN (kdim, mask &
       , zetnw, zftnw &
       , zeqnw, zfqnw &
       , zcptw, zqsw &
       , ztklevw, zqklevw &
!---wiso-code
       , lwiso, kwiso &
       ! Input - calculated in subroutine "Richtmeyer-Ocean"
       , zwisoeqnw, zwisofqnw &
       ! Input - calculated in subroutine "precalc_ocean"
       , zwisokinw, zwisoqsw &
       ! Output
       , zwisoqklevw &
!---wiso-code-end
       )

!---wiso-code
    USE mo_wiso,         ONLY: tnat
!---wiso-code-end

    INTEGER, INTENT(IN) :: kdim
!---wiso-code
    LOGICAL, INTENT(IN) :: lwiso
    INTEGER, INTENT(IN) :: kwiso
!---wiso-code-end
    LOGICAL, INTENT(IN) :: mask(kdim)
    REAL(dp), INTENT(IN)    :: zetnw(kdim), zftnw(kdim), zeqnw(kdim), zfqnw(kdim)
    REAL(dp), INTENT(IN)    :: zcptw(kdim), zqsw(kdim)
!---wiso-code
    REAL(dp), INTENT(IN), OPTIONAL    :: zwisoeqnw(kdim,kwiso), zwisofqnw(kdim,kwiso)
    REAL(dp), INTENT(IN), OPTIONAL    :: zwisokinw(kdim,kwiso)
    REAL(dp), INTENT(IN), OPTIONAL    :: zwisoqsw(kdim,kwiso)
!---wiso-code-end
    REAL(dp), INTENT(OUT)   :: ztklevw(kdim), zqklevw(kdim)
!---wiso-code
    REAL(dp), INTENT(OUT), OPTIONAL   :: zwisoqklevw(kdim,kwiso)

    ! Locals
    INTEGER  :: jt, jl
!---wiso-code-end
    
    !*  CALCULATION OF SKLEV AND QKLEV USING THE NEW SURFACE VALUES
    !*  ZSNEW AND ZQSNEW WHICH WERE CALCULATED IN SUBROUTINE SURFTEMP
    !
    
    ztklevw = 0._dp
    zqklevw = 0._dp
!---wiso-code
    IF (lwiso) THEN

    zwisoqklevw = 0._dp
    
    END IF
!---wiso-code-end

    WHERE(MASK)
       ztklevw(:)   = zetnw(:) * zcptw(:) + zftnw(:)
       zqklevw(:)   = zeqnw(:) * zqsw(:) + zfqnw(:)
    END WHERE

!---wiso-code
    IF (lwiso) THEN

    DO jt=1,kwiso
        WHERE(MASK)
            zwisoqklevw(:,jt) = zwisoeqnw(:,jt) * zwisokinw(:,jt) * zwisoqsw(:,jt) + zwisofqnw(:,jt)
        END WHERE
    END DO
    
    END IF
!---wiso-code-end

  END SUBROUTINE UPDATE_OCEAN
  !---------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE POSTPROC_OCEAN( &
       kdim,      klev,      &
       klevp1,    mask,      &
       zcfhw,     zqsw,      & 
       zqdif,     pqm1,      &
       pgeom1,    ztdif,     &
       zcptgz,    zcptw,     &
       ptsw,      zbnw,      &
       zbhw,      zriw,      &
       ptm1,      papm1,     &
       zx,        pum1,      &
       pvm1,                 &
       pocu,      pocv,      & 
       pevapw,               &
       pahfsw,    pahflw,    &
       pwind10w,  zt2w,      &
       paphm1,    zbmw,      &
       zdew2w,    zu10w,     &
       zv10w,                &
!---wiso-code
       lwiso, kwiso,         &
       ! Input
       zwisoqsw, zwisoqdif,  &
       pwisoqm1, zwisokinw,  &
       ! Output
       pwisoevapw            &
!---wiso-code-end
       )

    USE mo_constants,    ONLY: vtmpc1, vtmpc2, g, cpd, alv, tmelt, c3les, c4les, c3ies, c4ies, rd &
         ,c2es
    USE mo_physc2,       ONLY: cvdifts
    USE mo_time_control, ONLY: time_step_len

!---wiso-code
    USE mo_wiso, ONLY: nwiso, tnat
!---wiso-code-end

    INTEGER, INTENT(in)  :: kdim, klev, klevp1
!---wiso-code
    LOGICAL, INTENT(in)  :: lwiso
    INTEGER, INTENT(in)  :: kwiso
!---wiso-code-end
    LOGICAL, INTENT(in)  :: mask(kdim)
    REAL(dp), INTENT(in)     :: zcfhw(kdim), zqsw(kdim)
    REAL(dp), INTENT(in)     :: zqdif(kdim), pqm1(kdim), pgeom1(kdim) 
    REAL(dp), INTENT(in)     :: ztdif(kdim), zcptgz(kdim)
    REAL(dp), INTENT(in)     :: zcptw(kdim), ptsw(kdim), zbnw(kdim), zbhw(kdim), zriw(kdim)
    REAL(dp), INTENT(in)     :: ptm1(kdim), papm1(kdim), zx(kdim)
    REAL(dp), INTENT(in)     :: pum1(kdim), pvm1(kdim), paphm1(kdim,klevp1)
    REAL(dp), INTENT(in)     :: pocu(kdim), pocv(kdim)
    REAL(dp), INTENT(in)     :: zbmw(kdim)
!---wiso-code
    REAL(dp), INTENT(in), OPTIONAL :: zwisoqsw(kdim,kwiso), zwisoqdif(kdim,kwiso), zwisokinw(kdim,kwiso)
    REAL(dp), INTENT(in), OPTIONAL :: pwisoqm1(kdim,klev,kwiso)
!---wiso-code-end
    REAL(dp), INTENT(out)    :: pevapw(kdim), pahfsw(kdim), pahflw(kdim), pwind10w(kdim)
    REAL(dp), INTENT(out)    :: zt2w(kdim), zu10w(kdim), zv10w(kdim), zdew2w(kdim)
!---wiso-code
    REAL(dp), INTENT(out), OPTIONAL :: pwisoevapw(kdim,kwiso)
!---wiso-code-end

    ! Locals
    INTEGER  :: it(kdim)
    REAL(dp)     :: ztmst, zcons15, ztpfac1, ztpfac2, ztpfac3,zcons16, zhuv, zhtq, zephum
    REAL(dp)     :: zcoefw(kdim), zqnlev(kdim), zqhflw(kdim)
    REAL(dp)     :: ztnlev(kdim), zthflw(kdim)
    REAL(dp)     :: zrat(kdim), zcbn(kdim), zcbs(kdim), zcbu(kdim), zmerge(kdim), zred(kdim)
    REAL(dp)     :: zh2m(kdim), zqs1(kdim), zrh2m(kdim), zcvm3(kdim), zcvm4(kdim)
    REAL(dp)     :: zaph2m(kdim), zqs2(kdim), zq2m(kdim), zfrac(kdim)
    REAL(dp)     :: zmerge1(kdim)
!---wiso-code
    INTEGER 	 :: jt, jl
    REAL(dp)     :: zwisoqnlev(kdim,kwiso), zwisoqhflw(kdim,kwiso)
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
   
    pevapw = 0._dp
    pahfsw = 0._dp
    pahflw = 0._dp
    pwind10w = 0._dp
    zt2w = 0._dp
    zdew2w = 0._dp
    zu10w = 0._dp
    zv10w = 0._dp

!---wiso-code
    IF (lwiso) THEN

    pwisoevapw = 0._dp
    
    END IF
!---wiso-code-end
    WHERE(mask)
       
       !*         5.8     Surface fluxes of heat and moisture
       !
       !*    Moisture fluxes
       zcoefw(:)     = zcons15 * zcfhw(:)
       zqnlev(:)     = zqdif(:) ! - ztpfac3 * pqm1(:)
       zqhflw(:)     = zcoefw(:) * (zqnlev(:) - ztpfac2 * zqsw(:))
       !*    Sensible heat fluxes
       ztnlev(:)     = ztdif(:) ! - ztpfac3 * zcptgz(:)
       zthflw(:)     = zcoefw(:) * (ztnlev(:) - ztpfac2 * zcptw(:))
       zthflw(:)     = zthflw(:) - ptsw(:) * zcons16 * zqhflw(:)
      
       pevapw(:)     = zqhflw(:)
       pahfsw(:)     = zthflw(:)
       !     Latent heat fluxes
       pahflw(:)     = alv * zqhflw(:)
       !
       !     Compute new t2m, t2m_max t2m_min
       !
       zrat(:)       = zhtq / pgeom1(:)
       zcbn(:)       = LOG(1._dp + (EXP (zbnw(:)) - 1._dp) * zrat(:) )
       zcbs(:)       = -(zbnw(:) - zbhw(:)) * zrat(:)
       zcbu(:)       = -LOG(1._dp + (EXP (zbnw(:) - zbhw(:)) - 1._dp) * zrat(:))
       WHERE(zriw(:) .GT. 0._dp)
          zmerge(:)  = zcbs(:)
       ELSEWHERE
          zmerge(:)  = zcbu(:)
       ENDWHERE
       zred(:)       = (zcbn(:) + zmerge(:)) / zbhw(:)
       zh2m(:)       = zcptw(:) + zred(:) * (zcptgz(:) - zcptw(:))
       zt2w(:)       = (zh2m(:) - zhtq ) / (cpd * (1._dp + vtmpc2 * pqm1(:)))   
       !
       !           5.96   2M DEW POINT
       !
       it(:)         = NINT(ptm1(:) * 1000._dp)
!       WHERE (it(:) < jptlucu1 .OR. it(:) > jptlucu2) 
!          lookupoverflow = .TRUE.
!       ENDWHERE
       it(:)         = MAX( MIN ( it(:) , jptlucu2) , jptlucu1)
       zqs1(:)       = tlucua(it(:)) / papm1(:)
       zqs1(:)       = zqs1(:) / (1._dp - vtmpc1 * zqs1(:))
       zrh2m(:)      = MAX(zephum , pqm1(:) / zqs1(:))  
       WHERE(zt2w(:) .GT. tmelt)
          zcvm3(:)   = c3les
          zcvm4(:)   = c4les
       ELSEWHERE
          zcvm3(:)   = c3ies
          zcvm4(:)   = c4ies
       ENDWHERE
       zaph2m(:)     = paphm1(:,klevp1) * &
            (1._dp - zhtq / ( rd * zt2w(:) * (1._dp + vtmpc1 * pqm1(:) - zx(:))))
       it(:)         = NINT(zt2w(:) * 1000._dp)
!       WHERE (it(:) < jptlucu1 .OR. it(:) > jptlucu2) lookupoverflow = .TRUE.
       it(:)         = MAX( MIN(it(:), jptlucu2), jptlucu1)
       zqs2(:)       = tlucua(it(:)) / zaph2m(:)
       zqs2(:)       = zqs2(:) / (1._dp - vtmpc1 * zqs2(:))
       zq2m(:)       = zrh2m(:) * zqs2(:)
       zfrac(:)      = LOG(zaph2m(:) * zq2m(:) / (c2es * (1._dp+ vtmpc1 * zq2m(:)))) / zcvm3(:)
       zdew2w(:)     = MIN(zt2w(:), (tmelt - zfrac(:) * zcvm4(:)) / (1._dp - zfrac(:)))
       !
       !*          5.97   10M WIND COMPONENTS, MAX 10M WINDSPEED
       !
       zrat(:)       = zhuv / pgeom1(:)
       zcbn(:)       = LOG(1._dp + (EXP (zbnw(:)) - 1._dp) * zrat(:) )

       zcbs(:)       = -(zbnw(:) - zbmw(:)) * zrat(:)
       zcbu(:)       = -LOG(1._dp + (EXP (zbnw(:) - zbmw(:)) - 1._dp) * zrat(:))
       WHERE(zriw(:) .GT. 0._dp)
          zmerge1(:) =zcbs(:)
       ELSEWHERE
          zmerge1(:) =zcbu(:)
       ENDWHERE
       zred(:)       = (zcbn(:) + zmerge1(:)) / zbmw(:)
       zu10w(:)      = zred(:) * pum1(:)
       zv10w(:)      = zred(:) * pvm1(:)
       pwind10w(:)   = zred(:)*SQRT((pum1(:)-pocu(:))**2               &
                                   +(pvm1(:)-pocv(:))**2)
    ENDWHERE

!---wiso-code
    IF (lwiso) THEN

    DO jt=1,kwiso
        WHERE(mask)
            ! Surface fluxes of moisture over the ocean
            zcoefw(:)        = zcons15 * zcfhw(:)
            zwisoqnlev(:,jt) = zwisoqdif(:,jt) 
            zwisoqhflw(:,jt) = zcoefw(:) * zwisokinw(:,jt) * (zwisoqnlev(:,jt) - ztpfac2 * zwisoqsw(:,jt))
            ! Evaporation water isotopes
            pwisoevapw(:,jt) = zwisoqhflw(:,jt)
        ENDWHERE
    END DO
    
    END IF
!---wiso-code-end
    
  END SUBROUTINE POSTPROC_OCEAN
  !---------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE s_lake(          &
       kdim       , pseaice &
       , psiced   , palake  &
       , ptsi     , ptsw    &
       , pahflw   , pahfsw  &
       , pfluxres , ptrflw  &
       , psoflw   , pevapi  &
       , psni     , pcvsi   &
       , pahfres  , pfri    &
       )
    
    USE mo_constants,      ONLY: tmelt, rhoh2o, alf
    USE mo_time_control,   ONLY: delta_time

    INTEGER,   INTENT(in)   :: kdim
    REAL(dp),      INTENT(in)   :: palake(kdim)
    REAL(dp),      INTENT(inout):: psiced(kdim), pseaice(kdim)
    REAL(dp),      INTENT(inout):: ptsi(kdim), ptsw(kdim)  
    REAL(dp),      INTENT(in)   :: pahflw(kdim), pahfsw(kdim)
    REAL(dp),      INTENT(inout):: pfluxres(kdim) 
    REAL(dp),      INTENT(in)   :: ptrflw(kdim), psoflw(kdim)  
    REAL(dp),      INTENT(in)   :: pevapi(kdim), psni(kdim)    
    REAL(dp),      INTENT(inout):: pahfres(kdim)
    REAL(dp),      INTENT(in)   :: pfri(kdim), pcvsi(kdim)

    REAL(dp)   :: zfluxw(kdim), zts(kdim), zfres(kdim), zconhflx(kdim)
    REAL(dp)   :: zsubice(kdim), zhi(kdim) 
    REAL(dp)   :: zmcapdt, zdtime, zmixcap, zcpwater, zdmix, zfreez, zdice
    REAL(dp)   :: zdtrilf, zrhoilf, zrhoice, zmcaprilf, zalpha, zalphas
    REAL(dp)   :: zrho_sn, ziscond

    INTEGER :: jl

    ! CONSTANTS
    zdtime            = delta_time
    zalpha            = 2.1656_dp
    zalphas           = 0.31_dp
    zrho_sn           = 330._dp
    ziscond           = zalpha/zalphas*rhoh2o/zrho_sn
    zcpwater          = 4218._dp
    zdmix             = 10._dp
    zdice             = 0.10_dp
    zrhoice           = 910._dp
    zrhoilf           = zrhoice * alf
    zdtrilf           = zdtime/zrhoilf

    zfreez            = -zdice / zdtrilf

    zmixcap           = rhoh2o * zcpwater * zdmix
    zmcapdt           = zdtime / zmixcap
    zmcaprilf         = zmixcap / zrhoilf

!!$    do jl=1,kdim
!!$       if(palake(jl) .GE. 0.5_dp) THEN 
!!$          if(pseaice(jl) .LT. 0.5_dp) write(*,*) 'hallo seaice',pseaice(jl), jl
!!$          if(psiced(jl) .GE. zdice) write(*,*) 'hallo siced',psiced(jl), jl
!!$       END if
!!$    end do
!!$

    WHERE (palake(:).GE.0.5_dp)   
       ! lake points
       WHERE (pseaice(:) .LT. 0.5_dp)                                    ! open water
          
          zfluxw(:)   = pahflw(:) + pahfsw(:) + ptrflw(:) + psoflw(:)
          
          !--------       Lake temperature (ptsw)         ------------------------------------
          
          zts(:)            = ptsw(:) + zmcapdt * (zfluxw(:) + pfluxres(:))
          ptsi(:)           = tmelt
          pfluxres(:)       = 0._dp
          psiced(:)         = 0._dp
          WHERE (zts(:) .GE. tmelt)                                   ! open water (unchanged)
             ptsw(:)        = zts
          ELSEWHERE                                                   ! check ice formation
             ptsw(:)        = tmelt
             zfres(:)       = (zts(:) - tmelt) / zmcapdt                    ! < 0.
             WHERE (zfres(:) .LE. zfreez)                             ! ice formation
                psiced(:)   = zmcaprilf * (tmelt - zts(:))            ! > zdice
                pseaice(:)  = 1._dp
             ELSEWHERE
                pfluxres(:) = zfres(:)
             END WHERE
          END WHERE
       ELSEWHERE (psiced(:) .GE. zdice) 
         
          !--------       Ice thickness (psiced)        --------------------------------------
          
          zconhflx(:)       = zalpha * (ptsi(:) - tmelt) / (psiced(:) + ziscond * psni(:))
          zsubice(:)        = (1._dp - pcvsi(:)) * pevapi(:) * zdtime / zrhoice
          zhi(:)            = psiced(:) - zdtrilf * (zconhflx(:) + pfluxres(:)) + zsubice(:)
          ptsw(:)           = tmelt
          WHERE (zhi(:) .GE. zdice)
             psiced(:)      = zhi(:)
             pseaice(:)     = 1._dp
             pfluxres(:)    = 0._dp
          ELSEWHERE (zhi(:) .LE. 0._dp)                                   ! complete melting
             ptsw(:)        = tmelt - zhi(:) / zmcaprilf               ! ptsw > tmelt
             psiced(:)      = 0._dp
             pseaice(:)     = 0._dp
             pfluxres(:)    = 0._dp
          ELSEWHERE                                   ! incomplete melting
             psiced(:)      = zdice
             pseaice(:)     = 1._dp
             pfluxres(:)    = (zdice - zhi(:)) / zdtrilf
             pahfres(:)     = pahfres(:) - zdtime * pfri(:) * pfluxres(:)
          END WHERE
       END WHERE
    END WHERE
    
    
  END SUBROUTINE s_lake

  SUBROUTINE ocean_rad(               &
     kdim          , mask         &
     ,ztrdown      , ptsi         &
     ,pi0          , ptrsol       &
     ,palsoi       , p_albedo      &
     ,psofli       , ptrfli       &
     )

  USE mo_radiation,         ONLY: cemiss
  USE mo_constants,         ONLY: stbo

  INTEGER, INTENT(in)    :: kdim
  LOGICAL, INTENT(in)    :: mask(kdim)
  REAL(dp),    INTENT(in)    :: ptsi(kdim), pi0(kdim), ptrsol(kdim)
  REAL(dp),    INTENT(in)    :: palsoi(kdim), ztrdown(kdim), p_albedo(kdim)
  REAL(dp),    INTENT(out)   :: psofli(kdim), ptrfli(kdim)


  REAL(dp) :: zflxs(kdim)

  psofli = 0._dp
  ptrfli = 0._dp

  WHERE(mask)    
     
   
     ptrfli(:)  = ztrdown(:) - cemiss * stbo * ptsi(:)**4
     zflxs(:)   = pi0(:) * ptrsol(:)     
     psofli(:)  = (1._dp - palsoi(:)) * zflxs(:) / (1._dp - p_albedo(:))
      
  END WHERE

END SUBROUTINE ocean_rad

  
END MODULE MO_SURFACE_OCEAN
