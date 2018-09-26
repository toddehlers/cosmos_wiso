SUBROUTINE update_soiltemp (klon,nsoil                       &
                         , pts, psn                          &
                         , psodif, prgcgn                    &
                         , pgrndc, pgrndd                    &
                         , ptsoil                            &
                         , pgrndcapc, pgrndhflx              &
                         , lmask, ldglac)
!
!   AUTHOR:  FREDERIC HOURDIN     30/01/92
!
!            ADAPTED TO THE LMD-GCM BY JAN POLCHER  26/02/92
!            ADAPTED TO THE ECHAM-GCM BY JAN-PETER SCHULZ, MPI  03/02/96
!
!            J.-P. SCHULZ   MPI - OCTOBER 1997 :
!               ROUTINE USED FOR IMPLEMENTATION OF AN IMPLICIT
!               COUPLING BETWEEN LAND SURFACE AND ATMOSPHERE IN THE
!               ECHAM4 GCM.
!            U.SCHLESE DKRZ - NOVEMBER 1999  MODIFIED FOR ECHAM5
!            U.Schlese DKRZ - February 2000  new soil temperatures
!            L Kornblueh, MPI, January 2003, removed MERGE
!            
!            Adapted to JSBACH by Thomas Raddatz, Mai 2004
!
!   OBJECTIVE:  COMPUTATION OF:
!               THE GROUND TEMPERATURE EVOLUTION
!               THE GROUND SPECIFIC HEAT "CAPCAL"
!               THE SURFACE DIFFUSIVE FLUX FROM GROUND "F0"
!
!
!   METHOD:  IMPLICIT TIME INTEGRATION
!
!   CONSECUTIVES GROUND TEMPERATURES ARE RELATED BY:
!           T(K+1) = C(K) + D(K)*T(K)  (1)
!   THE COEFFICIENTS C (=pgrndc) AND D (=pgrndd) ARE COMPUTED AT THE
!   T-DT TIME-STEP.
!   ROUTINE STRUCTURE:
!   1)NEW TEMPERATURES ARE COMPUTED  USING (1)
!   2)C AND D COEFFICIENTS ARE COMPUTED FROM THE NEW TEMPERATURE
!     PROFILE FOR THE T+DT TIME-STEP
!   3)THE COEFFICIENTS A AND B ARE COMPUTED WHERE THE DIFFUSIVE
!     FLUXES AT THE T+DT TIME-STEP IS GIVEN BY
!            FDIFF = A + B TS(T+DT)
!     OR     FDIFF = F0 + CAPCAL (TS(T+DT)-TS(T))/DT
!            WITH F0 = A + B (TS(T))
!                 CAPCAL = B*DT
!
!     ------------------------------------------------------------------
!
!   DECLARATIONS:
!
USE mo_jsbach_constants   , ONLY: RhoH2O
USE mo_time_control       , ONLY: lstart, delta_time
USE mo_kind               , ONLY: dp

IMPLICIT NONE
!
!-----------------------------------------------------------------------
!  ARGUMENTS
!
  INTEGER, Intent(in)  ::  klon                 !! length of the vector
  INTEGER, Intent(in)  ::  nsoil                !! number of soil layers (fixed to 5) 
  REAL(dp), Intent(in)     ::  pts(klon)            !! surface temperature at top of soil [K]
  REAL(dp), Intent(in)     ::  psn(klon)            !! equivalent snow depth [m water]
  REAL(dp), Intent(in)     ::  psodif(klon)         !! soil temperature diffusivity [m^2/s]
  REAL(dp), Intent(in)     ::  prgcgn(klon)         !! soil heat capacity [J/m^3K]
  REAL(dp), Intent(inout)  ::  pgrndc(klon,nsoil)   !!
  REAL(dp), Intent(inout)  ::  pgrndd(klon,nsoil)   !!
  REAL(dp), Intent(inout)  ::  ptsoil(klon,nsoil)   !! soil temperature [K]
  REAL(dp), Intent(out)    ::  pgrndcapc(klon)      !!
  REAL(dp), Intent(out)    ::  pgrndhflx(klon)      !! ground heat flux
  LOGICAL, Intent(in)  ::  lmask(klon)
  LOGICAL, Intent(in)  ::  ldglac(klon)         !! glacier mask
!
!     ------------------------------------------------------------------
!
!  local Variables
!
  INTEGER :: jk
  REAL(dp) :: zso_cond(klon), zso_capa(klon)
  REAL(dp) :: z1(klon)
  REAL(dp) :: zd1(nsoil)
  REAL(dp) :: zdz1(klon,nsoil),   zdz2(klon,nsoil)
  REAL(dp) :: zkappa(klon,nsoil), zcapa(klon,nsoil)
  REAL(dp) :: zsnow_h(klon), zx1(klon), zx2(klon)
  REAL(dp) :: cdel(nsoil),cmid(nsoil)
  REAL(dp) :: zrici, zdifiz, zsn_cond, zsn_dens, zsn_capa
!
!     ------------------------------------------------------------------
!
!*    1.  SPECIFYING THE DEPTHS OF THE TEMPERATURE LEVELS.
!
!*    1.1 SOME CONSTANTS.
!
  zrici = 2.09e+06_dp                                 !! volumetric heat capacity of ice [j/m**3/k]
  zdifiz = 12.e-07_dp                                 !! temperature diffusivity of ice  [m**2/s]
  zsn_cond = 0.31_dp                                  !! snow thermal conductivity [j/s/m/k]
  zsn_dens = 330.0_dp                                 !! snow density              [kg/m**3]
  zsn_capa = 634500.0_dp                              !! snow  heat capacity   [j/m**3/k]
  cdel = (/0.065_dp,0.254_dp,0.913_dp,2.902_dp,5.700_dp/)         !! thicknesses of soil layers [m]
  cmid = (/0.0325_dp,0.192_dp,0.7755_dp,2.683_dp,6.984_dp/)       !! depth of mids of soil layers [m]
!
!*    1.2 COMPUTING SOME USEFUL CONSTANTS.
!
  DO jk = 1,nsoil-1
     zd1(jk) = 1._dp / (cmid(jk+1) - cmid(jk))
  END DO
!
!*    1.3 COMPUTE OF THE SOIL THERMAL CONDUCTIVITY [J/S/M/K] FROM
!*        THE SOIL TEMPERATURE DIFFUSIVITY [M**2/S].
!
  WHERE (lmask(:) .AND. ldglac(:))
    zso_capa(:) = zrici
    zso_cond(:) = zso_capa(:) * zdifiz
  ELSEWHERE (lmask(:))
    zso_capa(:) = prgcgn(:)
    zso_cond(:) = zso_capa(:) * psodif(:)
  END WHERE
!
!*    1.4 PRE-SET THERMAL CONDUCTIVITY AT ALL LEVELS.
!
  DO jk = 1,nsoil
     WHERE (lmask)
        zkappa(:,jk) = zso_cond(:)
        zcapa(:,jk)  = zso_capa(:)
     END WHERE
  END DO
!
!   --------------------------------------------------------------
!   COMPUTATION OF THE GROUND TEMPERATURES USING THE CGRD AND DGRD
!   COEFFICIENTS COMPUTED AT THE PREVIOUS TIME-STEP
!   --------------------------------------------------------------
!
  IF(.NOT.lstart) THEN
!
!   Upper layer
!
    ptsoil(:,1) = pts(:)
!
!   Deeper layers
!
    DO jk = 1,nsoil-1
      WHERE(lmask) ptsoil(:,jk+1) = pgrndc(:,jk) + pgrndd(:,jk) * ptsoil(:,jk)
    END DO
  END IF
!
!   ---------------------------------------------------------------
!   COMPUTATION OF THE CGRD AND DGRD COEFFICIENTS FOR THE NEXT STEP
!   ---------------------------------------------------------------
!
  WHERE (lmask) 
     zsnow_h(:) = psn(:) * RhoH2O / zsn_dens
!
!*       Special treatment for first layer
!
     WHERE ( zsnow_h(:) > cmid(2) )
        zcapa(:,1) = zsn_capa
        zkappa(:,1) = zsn_cond
     ELSEWHERE( zsnow_h(:) > 0.0_dp .AND. zsnow_h(:) <= cmid(2) )
        zx1 = zsnow_h(:) / cmid(2)
        zx2 = ( cmid(2) - zsnow_h(:)) / cmid(2)
        zcapa(:,1) = zx1 * zsn_capa + zx2 * zso_capa(:)
        zkappa(:,1) = 1._dp / ( zx1 / zsn_cond + zx2 / zso_cond(:) )
     ELSEWHERE
        zcapa(:,1) = zso_capa(:)
        zkappa(:,1) = zso_cond(:)
     ENDWHERE
  END WHERE
!
  DO jk = 2, nsoil - 2
    WHERE (lmask)
       WHERE ( zsnow_h(:) > cmid(jk+1) )
          zcapa(:,jk) = zsn_capa
          zkappa(:,jk) = zsn_cond
       ELSEWHERE ( zsnow_h(:) > cmid(jk) .AND. zsnow_h(:) <= cmid(jk+1) )
          zx1 = (zsnow_h(:) - cmid(jk)) * zd1(jk)
          zx2 = ( cmid(jk+1) - zsnow_h(:)) * zd1(jk)
          zcapa(:,jk) = zx1 * zsn_capa + zx2 * zso_capa(:)
          zkappa(:,jk) = 1._dp / ( zx1 / zsn_cond + zx2 / zso_cond(:) )
       ELSEWHERE
          zcapa(:,jk) = zso_capa(:)
          zkappa(:,jk) = zso_cond(:)
       ENDWHERE
    END WHERE
  END DO
!
  DO jk=1,nsoil
    WHERE (lmask) zdz2(:,jk) = zcapa(:,jk) * cdel(jk) / delta_time
  END DO
!
  DO jk=1,nsoil-1
    WHERE (lmask) zdz1(:,jk) = zd1(jk) * zkappa(:,jk)
  END DO
!
  WHERE (lmask)
     z1(:) = zdz2(:,nsoil) + zdz1(:,nsoil-1)
     pgrndc(:,nsoil-1) = zdz2(:,nsoil) * ptsoil(:,nsoil) / z1(:)
     pgrndd(:,nsoil-1) = zdz1(:,nsoil-1) / z1(:)
  END WHERE
!
  DO jk=nsoil-1,2,-1
     WHERE (lmask)
        z1(:) = 1._dp / (zdz2(:,jk) + zdz1(:,jk-1) + zdz1(:,jk) * (1._dp - pgrndd(:,jk)))
        pgrndc(:,jk-1) = (ptsoil(:,jk) * zdz2(:,jk) + zdz1(:,jk) * pgrndc(:,jk)) * z1(:)
        pgrndd(:,jk-1) = zdz1(:,jk-1) * z1(:)
     END WHERE
  END DO
!
!   ---------------------------------------------------------
!   COMPUTATION OF THE SURFACE DIFFUSIVE FLUX FROM GROUND AND
!   CALORIFIC CAPACITY OF THE GROUND:
!   ---------------------------------------------------------
!
  WHERE (lmask)
     pgrndhflx(:) = zdz1(:,1) * (pgrndc(:,1) + (pgrndd(:,1) - 1._dp) * ptsoil(:,1))
     pgrndcapc(:) = (zdz2(:,1) * delta_time + delta_time * (1._dp - pgrndd(:,1)) * zdz1(:,1))
  ELSEWHERE
     pgrndhflx = 0._dp
     pgrndcapc = 0._dp
  END WHERE

END SUBROUTINE update_soiltemp
