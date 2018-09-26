MODULE MO_SURFACE_BOUNDARY

  USE mo_kind, ONLY: dp

  IMPLICIT NONE
  
CONTAINS
  !
  !------------------------------------------------------------------------------------------------- 
  ! Some pre calculations for derived parameters for the lowest atmosphere level
  !
  SUBROUTINE atm_conditions ( kdim, &
       pxlm1   , pxim1,   &
       pgeom1  , ptm1,    &
       pqm1    , papm1,   &
       zx      , zcptgz,  &
       zteta1  , ztvir1,  &
       zfaxe   , zlteta1, &
       zqss               &
       )
    
    USE mo_constants, ONLY: cpd, vtmpc2, rd, vtmpc1, alv, als, tmelt
    USE mo_convect_tables ,   ONLY: tlucua, jptlucu1, jptlucu2
    
    INTEGER, INTENT(in)    :: kdim
    REAL(dp), INTENT(in)       :: pxlm1(kdim), pxim1(kdim), pgeom1(kdim)
    REAL(dp), INTENT(in)       :: ptm1(kdim), pqm1(kdim), papm1(kdim)
    REAL(dp), INTENT(out)      :: zx(kdim), zcptgz(kdim), zteta1(kdim), ztvir1(kdim), zfaxe(kdim)
    REAL(dp), INTENT(out)      :: zlteta1(kdim), zqss(kdim)
    
    INTEGER :: it(kdim)
    REAL(dp)    :: zkappa, zes(kdim)
    
    !--CONSTANTS
    !
    zkappa       = rd / cpd
    
    !-- MAIN CALCULATIONS
    !
    zx(:)        = pxlm1(:) + pxim1(:) 
    zcptgz(:)    = pgeom1(:) + ptm1(:) * cpd * (1._dp+ vtmpc2 * pqm1(:))
    zteta1(:)    = ptm1(:) * (100000._dp / papm1(:))**zkappa
    ztvir1(:)    = zteta1(:) * (1._dp+ vtmpc1 * pqm1(:) - zx(:))
    WHERE( ptm1(:) .GE. tmelt)
       zfaxe(:)  = alv
    ELSEWHERE
       zfaxe(:)  = als
    ENDWHERE
    zlteta1(:)   = zteta1(:) - zfaxe(:) / cpd * zteta1(:) / ptm1(:) * zx(:)
    it(:)        = NINT(ptm1(:) * 1000._dp)
    it(:)        = MAX (MIN(it(:), jptlucu2), jptlucu1)
    zes(:)       = MIN(tlucua(it(:)) / papm1(:) ,0.5_dp)
    zqss(:)      = zes(:) / (1._dp- vtmpc1 * zes(:))
    
  END SUBROUTINE atm_conditions
  !----------------------------------------------------------------------------------------------
  !
  SUBROUTINE blend_zq_zt(     &
       kdim                   &
       , pfrl                 &
       , pfrw     , pfri      &
       , zcfhl    , zcfhw     & 
       , zcfhi    , zcair     &
       , ztklevl  , ztklevw   &
       , ztklevi  , zqklevl   &
       , zqklevw  , zqklevi   &    
       , ztdif    , zqdif     &
       , land     , ocean     &
       , ice                  &
       )

    USE mo_physc2,         ONLY: cvdifts

    INTEGER,   INTENT(in)  :: kdim
    REAL(dp),      INTENT(in)  :: pfrl(kdim), pfrw(kdim), pfri(kdim)
    REAL(dp),      INTENT(in)  :: zcfhl(kdim), zcfhw(kdim), zcfhi(kdim)
    REAL(dp),      INTENT(in)  :: zcair(kdim)
    REAL(dp),      INTENT(in)  :: ztklevl(kdim), ztklevw(kdim), ztklevi(kdim)
    REAL(dp),      INTENT(in)  :: zqklevl(kdim), zqklevw(kdim), zqklevi(kdim)
    REAL(dp),      INTENT(out) :: ztdif(kdim), zqdif(kdim)

    LOGICAL,   INTENT(in)  :: land(kdim), ocean(kdim), ice(kdim)

    REAL(dp)       :: ztpfac2
    REAL(dp)       :: zgtsum(kdim), zgqsum(kdim)
    REAL(dp)       :: zgtl(kdim), zgtw(kdim), zgti(kdim)
    REAL(dp)       :: zgql(kdim), zgqw(kdim), zgqi(kdim)
    REAL(dp)       :: zero(kdim)

    INTEGER    :: jl

    ztpfac2    = 1._dp / cvdifts
    zero(:)    = 0._dp

    zgtl(:)    = pfrl(:) * MERGE(zcfhl(:),zero,land)
    zgtw(:)    = pfrw(:) * MERGE(zcfhw(:),zero,ocean)
    zgti(:)    = pfri(:) * MERGE(zcfhi(:),zero,ice)
    zgtsum(:)  = (zgtl(:) + zgtw(:) + zgti(:) ) / ztpfac2

    WHERE (pfrl(:) .LT. 1._dp)
       zgql(:) = MERGE(zgtl(:),zero,land) * MERGE(zcair(:),zero,land)
    ELSEWHERE
       zgql(:) = MERGE(zgtl(:),zero,land)
    ENDWHERE

    zgqw(:)    = MERGE(zgtw(:),zero,ocean)
    zgqi(:)    = MERGE(zgti(:),zero,ice)

    zgqsum(:)  = (zgql(:) + zgqw(:) + zgqi(:) ) / ztpfac2
    ztdif(:)   = (zgtl(:) * MERGE(ztklevl(:),zero,land) & 
         + zgtw(:) * MERGE(ztklevw(:),zero,ocean) &
         + zgti(:) * MERGE(ztklevi(:),zero,ice))  / zgtsum(:)
    zqdif(:)   = (zgql(:) * MERGE(zqklevl(:),zero,land) &
         + zgqw(:) * MERGE(zqklevw(:),zero,ocean) &
         + zgqi(:) * MERGE(zqklevi(:),zero,ice)) / zgqsum(:)
  
  END SUBROUTINE blend_zq_zt
  !----------------------------------------------------------------------------------------------

!---wiso-code
  SUBROUTINE blend_zq_zt_wiso(kdim    &
       , pfrl, pfrw, pfri             &
       , zcfhl, zcfhw, zcfhi          &
       , land, ocean, ice             &
       , zcair, zwisocair             &
       , zwisocair_fra                &
       , zwisoqklevl                  &
       , zwisoqklevw  , zwisoqklevi   &
       , zwisoqdif                    &
       )

    USE mo_physc2,         ONLY: cvdifts
    USE mo_wiso,           ONLY: nwiso

    INTEGER,   INTENT(in)  :: kdim
    REAL(dp),      INTENT(in)  :: pfrl(kdim), pfrw(kdim), pfri(kdim)
    REAL(dp),      INTENT(in)  :: zcfhl(kdim), zcfhw(kdim), zcfhi(kdim)
    REAL(dp),      INTENT(in)  :: zcair(kdim)
    REAL(dp),      INTENT(in)  :: zwisocair(kdim,nwiso)
    REAL(dp),      INTENT(in)  :: zwisocair_fra(kdim,nwiso)
    REAL(dp),      INTENT(in)  :: zwisoqklevl(kdim,nwiso), zwisoqklevw(kdim,nwiso), zwisoqklevi(kdim,nwiso)
    REAL(dp),      INTENT(out) :: zwisoqdif(kdim,nwiso)

    LOGICAL,   INTENT(in)  :: land(kdim), ocean(kdim), ice(kdim)

    REAL(dp)       :: ztpfac2
    REAL(dp)       :: zgqw(kdim), zgqi(kdim)
    REAL(dp)       :: zero(kdim)
    REAL(dp)       :: zgql_wiso(kdim,nwiso), zgqsum_wiso(kdim,nwiso)
    INTEGER        :: jt, jl

    ztpfac2    = 1._dp / cvdifts
    zero(:)    = 0._dp

    zgqw(:)    = pfrw(:) * MERGE(zcfhw(:),zero,ocean)
    zgqi(:)    = pfri(:) * MERGE(zcfhi(:),zero,ice)
    
    DO jt=1,nwiso 
    WHERE (pfrl(:) .LT. 1._dp)
          zgql_wiso(:,jt) = (pfrl(:) * MERGE(zcfhl(:),zero,land)) * MERGE(zcair(:),zero,land)
    ELSE WHERE
          zgql_wiso(:,jt) = (pfrl(:) * MERGE(zcfhl(:),zero,land))
    END WHERE

    zgqsum_wiso(:,jt)     = (zgql_wiso(:,jt) + zgqw(:) + zgqi(:)) / ztpfac2
    
       zwisoqdif(:,jt)   = (zgql_wiso(:,jt) * MERGE(zwisoqklevl(:,jt),zero,land)                     &
                         +  zgqw(:) * MERGE(zwisoqklevw(:,jt),zero,ocean)                            &
                         +  zgqi(:) * MERGE(zwisoqklevi(:,jt),zero,ice)) / zgqsum_wiso(:,jt)
    END DO

  END SUBROUTINE blend_zq_zt_wiso
!---wiso-code-end
  !----------------------------------------------------------------------------------------------
  !
  SUBROUTINE average_ustar( &
       kdim      , klevp1   &
       , pfrl    , pfrw     &
       , pfri    , zustl    &
       , zustw   , zusti    &
       , ptm1    , pqm1     &
       , zx      , paphm1   &
       , zustarm , land     &
       , ocean   , ice      &
       )
        
    USE mo_constants,      ONLY: vtmpc1, g, rd
    USE mo_physc2,         ONLY: cvdifts
    USE mo_time_control,   ONLY: time_step_len
    
    INTEGER, INTENT(IN)    :: kdim, klevp1
    REAL(dp),    INTENT(IN)    :: pfrl(kdim), pfrw(kdim), pfri(kdim)
    REAL(dp),    INTENT(IN)    :: zustl(kdim), zustw(kdim), zusti(kdim)
    REAL(dp),    INTENT(in)    :: ptm1(kdim), pqm1(kdim), zx(kdim)
    REAL(dp),    INTENT(in)    :: paphm1(kdim,klevp1)
    LOGICAL, INTENT(in)    :: land(kdim), ocean(kdim), ice(kdim)
    
    REAL(dp),    INTENT(OUT) :: zustarm(kdim)
    
    REAL(dp)  :: zcons12, ztmst, ztpfac1, ztkemin
    REAL(dp)  :: zust(kdim), zero(kdim)
    
    !   CONSTANTS
    ztmst   = time_step_len
    ztpfac1 = cvdifts 
    zcons12 = ztpfac1*ztmst*g/rd   
    zero(:) = 0._dp
 
    !----------------------------------------------------------     
    !      TKE BOUNDARY CONDITION (MAILHOT/BENOIT, 1982)
    !      For pbl-height calculation
    zust(:)        = pfrl(:) * MERGE(zustl(:),zero,land) + pfrw(:) * MERGE(zustw(:),zero,ocean) &
         + pfri(:) * MERGE(zusti(:),zero,ice)    
    zustarm(:)     = SQRT(zust(:) * ptm1(:) *  &
         ( 1._dp + vtmpc1 * pqm1(:) - zx(:)) / (zcons12 * paphm1(:,klevp1)))

  END SUBROUTINE average_ustar
  !
  !-------------------------------------------------------------------------------------------------
  SUBROUTINE average_tvh_qsurf( &
       kdim       , ptslm1      & 
       , zhsoil   , ptsw        &
       , ptsi     , pfrl        &
       , pfrw     , pfri        &
       , zcsat    , zqsw        &
       , zqsi     , ztvh        &
       , zqsurf   , zqsl        &
       , land     , ocean       &
       , ice                    &
       )

    USE mo_constants,      ONLY: vtmpc1

    INTEGER, INTENT(in)   :: kdim
    REAL(dp), INTENT(in)      :: ptslm1(kdim), zhsoil(kdim), ptsw(kdim), ptsi(kdim)
    REAL(dp), INTENT(in)      :: pfrl(kdim), pfrw(kdim), pfri(kdim)
    REAL(dp), INTENT(in)      :: zcsat(kdim), zqsw(kdim), zqsi(kdim), zqsl(kdim)
    REAL(dp), INTENT(out)     :: ztvh(kdim), zqsurf(kdim)
    LOGICAL, INTENT(in)   :: land(kdim), ocean(kdim), ice(kdim)

    REAL(dp) :: ztvlan(kdim), ztvsea(kdim), ztvice(kdim), zero(kdim)

    zero(:)        = 0._dp

!    ztvlan(:)      = ptslm1(:) * (1. + vtmpc1 * zhsoil(:)*zqsl(:))
    ztvlan(:)      = MERGE(ptslm1(:),zero,land)  * (1._dp + vtmpc1 * MERGE(zhsoil(:)*zqsl(:),zero,land))
    ztvsea(:)      = MERGE(ptsw(:)  ,zero,ocean) * (1._dp + vtmpc1 * MERGE(zqsw(:) ,zero,ocean))
    ztvice(:)      = MERGE(ptsi(:)  ,zero,ice)   * (1._dp + vtmpc1 * MERGE(zqsi(:) ,zero,ice))
    ztvh(:)        = pfrl(:) * ztvlan(:) + pfrw(:) * ztvsea(:) + pfri(:) * ztvice(:)
    zqsurf(:)      = pfrl(:) * MERGE(zqsl(:),zero,land) * MERGE(zcsat(:),zero,land) &
         + pfrw(:) * MERGE(zqsw(:),zero,ocean) + pfri(:) * MERGE(zqsi(:),zero,ice)

  END SUBROUTINE average_tvh_qsurf
  !
  !-------------------------------------------------------------------------------------------------
  SUBROUTINE surface_box_average(        &
       kdim      , land     &
       , ice     , ocean    &
       , box_avg , jrow     &
       , m_land  , m_ocean  &
       , m_ice   , f_land   &
       , f_ocean , f_ice    &
       )

    INTEGER,     INTENT(in)  :: kdim, jrow
    REAL(dp),    INTENT(in)  :: land(kdim), ice(kdim), ocean(kdim)
    LOGICAL,     INTENT(in)  :: m_land(kdim), m_ocean(kdim), m_ice(kdim)
    REAL(dp),    INTENT(in)  :: f_land(kdim), f_ocean(kdim), f_ice(kdim)
    REAL(dp),    INTENT(out) :: box_avg(kdim)
    REAL(dp)                 :: zero(kdim)

    zero(:)         = 0._dp
    box_avg(1:kdim) = f_land * MERGE(land,zero,m_land) + &
         f_ice * MERGE(ice,zero,m_ice)  + f_ocean * MERGE(ocean,zero,m_ocean)
    
  END SUBROUTINE surface_box_average

!---wiso-code
  SUBROUTINE surface_box_avg_wiso(&
       kdim      , land           &
       , ice     , ocean          &
       , box_avg , jrow , kwiso   &
       , m_land  , m_ocean        &
       , m_ice   , f_land         &
       , f_ocean , f_ice          &
       )

    INTEGER,     INTENT(in)  :: kdim, jrow, kwiso
    REAL(dp),    INTENT(in)  :: land(kdim,kwiso), ice(kdim,kwiso), ocean(kdim,kwiso)
    LOGICAL,     INTENT(in)  :: m_land(kdim), m_ocean(kdim), m_ice(kdim)
    REAL(dp),    INTENT(in)  :: f_land(kdim), f_ocean(kdim), f_ice(kdim)
    REAL(dp),    INTENT(out) :: box_avg(kdim,kwiso)
    REAL(dp)                 :: zero(kdim)
    INTEGER                  :: jt, jl

    zero(:)         = 0._dp
    DO jt=1,kwiso
       box_avg(1:kdim,jt) = f_land * MERGE(land(1:kdim,jt),zero,m_land) + &
                         f_ice * MERGE(ice(1:kdim,jt),zero,m_ice)  + f_ocean * MERGE(ocean(1:kdim,jt),zero,m_ocean)
    END DO
    
  END SUBROUTINE surface_box_avg_wiso
!---wiso-code-end
  !
  !-------------------------------------------------------------------------------------------------
  !
  SUBROUTINE longwave_down_rad(     &
       kdim          , pemter &
       ,ptslm1       , ptsw    &
       ,ptsi         , pfrl    &
       ,pfrw         , pfri    &
       ,ztrdown         &
       )

    USE mo_radiation,        ONLY: cemiss
    USE mo_constants,        ONLY: stbo

    INTEGER, INTENT(in)    :: kdim
    REAL(dp), INTENT(in)       :: pemter(kdim), ptslm1(kdim), ptsw(kdim)
    REAL(dp), INTENT(in)       :: ptsi(kdim), pfrl(kdim), pfrw(kdim), pfri(kdim)
    REAL(dp), INTENT(out)      :: ztrdown(kdim)

    REAL(dp)    :: zteff4
    INTEGER :: jl
    
     DO jl = 1,kdim
        zteff4=pfrl(jl)*ptslm1(jl)**4                                  &
              +pfri(jl)*ptsi(jl)**4                                    &
              +pfrw(jl)*ptsw(jl)**4
        ztrdown(jl)=pemter(jl)+cemiss*stbo*zteff4
     END DO

   END SUBROUTINE longwave_down_rad
  !
  !-------------------------------------------------------------------------------------------------
  !
END MODULE MO_SURFACE_BOUNDARY
