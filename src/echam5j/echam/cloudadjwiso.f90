SUBROUTINE cloudadjwiso(kproma,kbdim,klev,kwiso,kk,              &
                        ptmst,pcons2,pdp,                        &
                        ptpone,pqpone,pqspone,pxlpone,           &
                        paclc,                                   &
                        pwisoqm1,  pwisoqte,                     &
                        pwisoxlm1, pwisoxlte,                    &
                        prfl, pwisorfl)   
!
!      G.HOFFMANN          MPI MET, HAMBURG      1992 
!      M. WERNER,          MPI BGC, JENA,        2004
!      M. WERNER           AWI, BREMERHAVEN      2009
!
!      PURPOSE
!      -------
!      CALCULATES THE FRACTIONATION FACTOR (LIQUID/VAPOUR),
!      GETS CLOUD LIQUID ISOTOPE AND VAPOUR AND THEN
!      LIQUID PRECIPITATION AND VAPOUR IN EQUILIBRIUM
!
!      INTERFACE
!      ---------
!      THIS SUBROUTINE IS CALLED FROM
!        *CLOUD*
!
!      INPUT  TEMPERATURE , 
!             CLOUD LIQUID WATER,
!             VAPOUR AND SATURATION MIXING RATIO OF VAPOUR,
!             PRECIPITATION AND FRACTIONAL CLOUD COVER,
!             OLD ISOTOPE MIXING RATIO AND TENDENCY,
!             ISOTOPE PRECIPITATION
!      OUTPUT NEW ISOTOPE TENDENCIES 
!
!      EXTERNALS
!      ---------
!      NONE

  USE mo_kind,           ONLY : dp
  USE mo_cloud,          ONLY : tmelt  
  USE mo_wiso,           ONLY : twisoeqls, wiso_frac_liq, tnat, cwisomin, cwisosec

  IMPLICIT NONE

! input arguments  
  INTEGER, INTENT(IN) :: kproma,kbdim,klev,kwiso,kk

  REAL(dp), INTENT(IN) :: ptmst,pcons2

  REAL(dp), INTENT(IN) :: pdp(kbdim),ptpone(kbdim),pqpone(kbdim),              &
                          pqspone(kbdim),pxlpone(kbdim),                       &
                          pwisoqm1(kbdim,klev,kwiso),                          &
                          pwisoxlm1(kbdim,klev,kwiso),                         &
                          prfl(kbdim),paclc(kbdim,klev)

! input/output arguments  
  REAL(dp), INTENT(INOUT) :: pwisoqte(kbdim,klev,kwiso),                       &
                             pwisoxlte(kbdim,klev,kwiso),                      &
                             pwisorfl(kbdim,kwiso)
  
! local variables
  REAL(dp) :: zwisofracliq(kbdim,kwiso)

  REAL(dp) :: zdelta, zwisoqpone,zwisoxlpone,                                    &
              zrafac,zqrain,zwisoqrain,                                        &
              zqinc,zqoutc,zwisoqinc,zwisoqoutc,zdeltab,zdeltamb,zdenom

  INTEGER :: jl,jt
    

  CALL wiso_frac_liq(kproma,kbdim,kwiso,ptpone,zwisofracliq)

  DO jt=1,kwiso        
      DO jl=1,kproma
          
! Equilibration Part I: bring liquid cloud water in isotopic equilibrium with surrounding vapour

       zwisoqpone =pwisoqm1(jl,kk,jt) +pwisoqte(jl,kk,jt) *ptmst 
       zwisoxlpone=pwisoxlm1(jl,kk,jt)+pwisoxlte(jl,kk,jt)*ptmst
          
       IF (pqpone(jl).ge.0.0_dp .and. pxlpone(jl).ge.0.0_dp .and. zwisoqpone.ge.0.0_dp .and. zwisoxlpone.ge.0.0_dp) THEN

        IF (pxlpone(jl).ne.0.0_dp) THEN ! if some liquid cloud water will still exist in next time step: 
                                        ! assume fractionation between liquid and vapour (closed system)

         IF (pqpone(jl).gt.cwisomin) THEN  ! apply fractionation for positive vapour values, only

          IF ((pqpone(jl)+zwisofracliq(jl,jt)*pxlpone(jl)).lt.cwisomin) THEN
            zdelta=0.0_dp
          ELSE
            zdelta=(zwisofracliq(jl,jt)*pxlpone(jl)*zwisoqpone-zwisoxlpone*pqpone(jl))         &
                  /(pqpone(jl)+zwisofracliq(jl,jt)*pxlpone(jl))
          ENDIF

          ! limit changes to absolut amount of available vapour and liquid water
          IF (zwisoqpone.GE.0._dp) THEN
            zdelta=MIN(zdelta,zwisoqpone)
            zdelta=MAX(zdelta,-zwisoxlpone)
          ELSE
            zdelta=MAX(zdelta,zwisoqpone)
            zdelta=MIN(zdelta,-zwisoxlpone)
          ENDIF
          
          pwisoqte(jl,kk,jt) =pwisoqte(jl,kk,jt) -zdelta/ptmst
          pwisoxlte(jl,kk,jt)=pwisoxlte(jl,kk,jt)+zdelta/ptmst
         
         ENDIF
        
        ELSE  ! if no default cloud water will exist in next time step: adjust isotope tendencies, accordingly

         zdelta=pwisoxlm1(jl,kk,jt)+pwisoxlte(jl,kk,jt)*ptmst
         pwisoqte(jl,kk,jt) =pwisoqte(jl,kk,jt) +zdelta/ptmst
         pwisoxlte(jl,kk,jt)=pwisoxlte(jl,kk,jt)-zdelta/ptmst
 
        ENDIF
        
       ENDIF

! Equilibration Part II: bring precipitation in isotopic equilibrium with surrounding vapour,
! distinguish here between rain still within a cloud versus rain already outside a cloud

       ! convert rain to same units as water vapour
       zrafac=pcons2*pdp(jl)
       zqrain=prfl(jl)/zrafac
       zwisoqrain=pwisorfl(jl,jt)/zrafac
       zwisoqpone=pwisoqm1(jl,kk,jt)+pwisoqte(jl,kk,jt)*ptmst

       IF (pqpone(jl).ge.0.0_dp .and. zqrain.ge.0.0_dp .and. zwisoqpone.ge.0.0_dp .and. zqrain.ge.0.0_dp) THEN

! if some rain will exist in next time step:
! assume fractionation between rain and vapour (closed system)

        IF (zqrain.ne.0.0_dp) THEN
           
         IF (pqpone(jl).gt.cwisomin) THEN   ! apply fractionation for positive vapour values, only

          ! if a cloud exist: assume that vapour in cloud (zqinc,zwisoinc) is at saturation value with
          ! a delta value equal to the value of the whole grid cell
          IF (paclc(jl,kk).gt.1.e-10_dp) THEN
            zqinc=pqspone(jl)
            zdelta=tnat(jt)
            IF (pqpone(jl).GT.cwisomin) zdelta=zwisoqpone/pqpone(jl)
            IF (ABS(1._dp-zdelta).LT.cwisosec) zdelta=1._dp
            zwisoqinc=pqspone(jl)*zdelta
          ELSE
            zqinc=pqpone(jl) 
            zwisoqinc=zwisoqpone
          ENDIF 

          ! if some cloud-free areas exist: calculate normal vapour and isotope value outside cloud (zqoutc, zwisoqoutc)
          ! as the residual of the whole grid box - cloud values
          IF (paclc(jl,kk).lt.(1.0_dp-1.0e-10_dp)) THEN
            zqoutc=(pqpone(jl)-paclc(jl,kk)*zqinc)/(1._dp-paclc(jl,kk))
            zwisoqoutc=(zwisoqpone-paclc(jl,kk)*zwisoqinc)/(1._dp-paclc(jl,kk))
          ELSE
            zqoutc=pqpone(jl)
            zwisoqoutc=zwisoqpone
          ENDIF

          ! correct for negative values
          IF ((zwisoqoutc).lt.0.0_dp) THEN
            zwisoqoutc=zwisoqpone
            zwisoqinc=zwisoqpone
            zqinc=pqpone(jl)
            zqoutc=pqpone(jl)
          ENDIF

          ! calculate isotope equilibration of precipitation in the cloudy part of the grid box
          zdenom=zwisofracliq(jl,jt)*twisoeqls*zqrain+zqinc
          IF (zdenom.lt.cwisomin) THEN
            zdeltab=0.0_dp
          ELSE
            zdeltab=twisoeqls*(zwisoqinc*zwisofracliq(jl,jt)*zqrain-zqinc*zwisoqrain)/zdenom
          ENDIF


          ! calculate isotope equilibration of precipitation in the cloud-free part of the grid box
          zdenom=zwisofracliq(jl,jt)*twisoeqls*zqrain+zqoutc
          IF (zdenom.lt.cwisomin) THEN
            zdeltamb=0.0_dp
          ELSE
            zdeltamb=twisoeqls*(zwisoqoutc*zwisofracliq(jl,jt)*zqrain-zqoutc*zwisoqrain)/zdenom
          ENDIF


          ! calculate isotope changes of the whole grid box (cloudy and cloud-free part)
          zdelta=paclc(jl,kk)*zdeltab+(1.-paclc(jl,kk))*zdeltamb
          
         ELSE ! nothing happens if no positive vapour exist in next time step (pqpone(jl).le.cwisomin) 
         
          zdelta=0.0_dp

         ENDIF

! if no rain will exist in next time step (=all rain evaporates):
! - add isotope rain tendency to vapour without any fractionation  
        ELSE

          zdelta=-pwisorfl(jl,jt)/zrafac

        ENDIF

! add up isotopic changes caused by equilibration of rain with vapour

        pwisorfl(jl,jt)=pwisorfl(jl,jt)+zdelta*zrafac
        pwisoqte(jl,kk,jt)=pwisoqte(jl,kk,jt)-zdelta/ptmst
        
       ENDIF

      END DO
   END DO

  RETURN
END SUBROUTINE cloudadjwiso
