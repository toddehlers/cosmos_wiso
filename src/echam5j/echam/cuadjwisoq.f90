SUBROUTINE cuadjwisoq(kproma,   kbdim,      klev,      kwiso,    kk, &
                      ptu,      pqu,        pwisoqu,   pwisolu,      &
                      pqcond,   pwisoqcond, ldflag)
!
!      G.HOFFMANN          MPI MET, HAMBURG      1992 
!      M. WERNER,          MPI BGC, JENA,        2004
!      M. WERNER           AWI, BREMERHAVEN      2009
!
!      PURPOSE
!      -------
!      CALCULATES THE FRACTIONATION FACTOR (LIQUID/VAPOUR SOLID/VAPOUR)
!      AND PUTS THE LIQUID (OR SOLID) ISOTOPE CONDENSATE IN (EFFECTIVE) 
!      EQUILIBRIUM WITH THE SURROUNDING VAPOUR
!
!      INTERFACE
!      ---------
!      THIS SUBROUTINE IS CALLED FROM
!        *CUBASE* (CONDENSATION OF CLOUD LIQUID ISOTOPE WATER)

!      INPUT T, 
!            TWO WATER VALUES (E.G.,OLD VAPOUR AND CONDENSATE)
!            TWO ISOTOPE VALUES (E.G.,OLD ISOTOPE VAPOUR AND LIQUID)
!      OUTPUT ONE OR TWO ISOTOPE VALUES (E.G.,NEW ISOTOPE CONDENSATE)
!
!      EXTERNALS
!      ---------
!      NONE

  USE mo_kind,           ONLY: dp
  USE mo_constants,      ONLY: tmelt
  USE mo_wiso,           ONLY: wiso_frac_liq_ice, cwisomin, cwisosec, twisoice

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT(IN)      :: kproma, kbdim, klev, kwiso, kk
  
  !  Array arguments with intent(In):
  REAL(dp), INTENT(IN)     :: ptu(kbdim,klev), pqu(kbdim,klev), pqcond(kbdim,klev)
  LOGICAL, INTENT(IN)      :: ldflag(kbdim)

  !  Array arguments with intent(InOut):
  REAL(dp), INTENT(INOUT)  :: pwisoqu(kbdim,klev,kwiso), pwisolu(kbdim,klev,kwiso)
  REAL(dp), INTENT(INOUT)  :: pwisoqcond(kbdim,klev,kwiso)
  
  !  Local arrays: 
  REAL(dp) :: zwisofracliq(kbdim,kwiso), zwisofracice(kbdim,kwiso)

  !  Local scalars: 
  REAL(dp) :: zqliq,  zqice
  REAL(dp) :: zqvapo, zqvapl,    zqconl, zqconi 
  REAL(dp) :: zwisoqvapo
  REAL(dp) :: zdenom, zevrel,    zquot
  INTEGER  :: jl, jt
  LOGICAL  :: LO

!  Executable statements 

! calculate fractionation factors for liq.-vap. and sol.-vap.
  CALL wiso_frac_liq_ice(kproma,kbdim,kwiso,ptu(:,kk),zwisofracliq,zwisofracice)

! adjusting isotope vapour and liquid depending on fract.
  DO jt=1,kwiso
    DO jl=1,kproma

! divide into ice and liquid
      IF (ptu(jl,kk).gt.tmelt) THEN
        zqliq=1._dp
        zqice=0._dp
      ELSEIF (ptu(jl,kk).lt.twisoice) THEN
        zqice=1._dp
        zqliq=0._dp
      ELSE
        zqice=(tmelt-ptu(jl,kk))/(tmelt-twisoice)
        zqliq=1._dp-zqice
      ENDIF
     
! fractionation of ice and liquid
      IF (ldflag(jl).and.pqcond(jl,kk).gt.cwisomin) THEN
      
        zqvapo=pqu(jl,kk)+pqcond(jl,kk)
        zqconl=zqliq*pqcond(jl,kk)
        zqvapl=zqvapo-zqconl

        zdenom=zqvapl+zwisofracliq(jl,jt)*zqconl
        lo=abs(zdenom).lt.cwisomin
        IF (lo) THEN
          pwisoqcond(jl,kk,jt)=0._dp
        ELSE
          pwisoqcond(jl,kk,jt)=zwisofracliq(jl,jt)*zqconl*pwisoqu(jl,kk,jt)/zdenom
        ENDIF
 
        zqconi=zqice*pqcond(jl,kk)
        zwisoqvapo=pwisoqu(jl,kk,jt)-pwisoqcond(jl,kk,jt)
          
        IF (zqvapl.gt.cwisomin) THEN
          zevrel=zqconi/zqvapl
          IF (abs(1.-zevrel).lt.cwisosec) zevrel=1._dp
        ELSE
          zevrel=0._dp
        ENDIF
           
        zquot=1._dp-zevrel
        zdenom=(1._dp-zquot**zwisofracice(jl,jt))

        pwisoqcond(jl,kk,jt)=pwisoqcond(jl,kk,jt)+zwisoqvapo*zdenom

        pwisoqu(jl,kk,jt)=pwisoqu(jl,kk,jt)-pwisoqcond(jl,kk,jt)
        pwisolu(jl,kk,jt)=pwisolu(jl,kk,jt)+pwisoqcond(jl,kk,jt)

      ELSE
      
        zqvapo=pqu(jl,kk)+pqcond(jl,kk)
        zwisoqvapo=pwisoqu(jl,kk,jt)
        
        zevrel=1._dp
        IF ((ABS(pqcond(jl,kk)).GT.cwisomin).AND.(ABS(zqvapo).GT.cwisomin)) zevrel=pqcond(jl,kk)/zqvapo
        IF (ABS(1._dp-zevrel).LT.cwisosec) zevrel=1._dp
        
        pwisoqu(jl,kk,jt)=pwisoqu(jl,kk,jt)-zwisoqvapo*zevrel
        
      ENDIF
    
    END DO
  END DO

  RETURN

END SUBROUTINE cuadjwisoq
