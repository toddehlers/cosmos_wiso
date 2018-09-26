SUBROUTINE cuwisoequ(kproma, kbdim, klev, klevp1, ktopm2, ldcum, kwiso, &
                     kctop,          kcbot,       paphp1,              &
                     ptp1,           ptte,                             & 
                     pqp1,           pqte,                             &
                     pwisoqp1,       pwisoqte,                         &
                     ptu,            pqu,         pwisoqu,             &
                     pten,                                             &
                     prfl_tmp,       psfl_tmp,                         &
                     pwisorfl,       pwisosfl,                         &
                     pwisodmfup,     pwisodmfdp,                       &
                     pmelt_tmp,                                        &
                     pprec_frac,     prain_tmp,                        &
                     pwisorsfc,      pwisossfc,                        &
                     pwisoaprc,      pwisoaprs)
                     
!
! PURPOSE:
!
! GET ISOTOPE CONCENTRATION IN RAIN INTO EQUILIBRIUM WITH 
! THE CONCENTRATION IN VAPOUR.
!
! INTERFACE:
!
! THIS SUBROUTINE IS CALLED FROM
! *CUMASTR*
!
! AUTHORS:
!
! G.HOFFMANN, MPI MET, HAMBURG, 1992 
! ADAPTED TO F90: M. WERNER, MPI BGC, JENA, 2004 
! ADAPTED TO ECHAM5: M. WERNER, AWI, BREMERHAVEN, 2009
!


USE mo_kind,           ONLY: dp
USE mo_constants,      ONLY: tmelt, g, vtmpc1
USE mo_convect_tables, ONLY: tlucua,                                   &
                             tlucub,                                   &
                             jptlucu1, jptlucu2,                       &
                             lookuperror, lookupoverflow
USE mo_time_control,   ONLY: time_step_len, delta_time
USE mo_wiso,           ONLY: talphal1, talphal2, talphal3,             &
                             thumwiso1, thumwiso2, twisoeqcu,          &
                             tdifrel,   tdifexp,                       &
                             nwisotyp,  tnat, twisoice,                &
                             cwisomin,  cwisosec

IMPLICIT NONE
!  
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktopm2, kwiso

REAL(dp):: paphp1(kbdim,klevp1),                                       &
           ptp1(kbdim,klev),           ptte(kbdim,klev),               &             
           pqp1(kbdim,klev),           pqte(kbdim,klev),               &
           pwisoqp1(kbdim,klev,kwiso), pwisoqte(kbdim,klev,kwiso),     &
           ptu(kbdim,klev),            pqu(kbdim,klev),                &
           pwisoqu(kbdim,klev,kwiso),                                  &
           pten(kbdim,klev),                                           &
           prfl_tmp(kbdim),            psfl_tmp(kbdim),                &
           pwisorfl(kbdim,kwiso),      pwisosfl(kbdim,kwiso),          &
           pwisodmfup(kbdim,klev,kwiso), pwisodmfdp(kbdim,klev,kwiso), &
           pmelt_tmp(kbdim,klev),                                      &
           pprec_frac(kbdim,klev),    prain_tmp(kbdim,klev),           &
           pwisorsfc(kbdim,kwiso),    pwisossfc(kbdim,kwiso),          &
           pwisoaprc(kbdim,kwiso),    pwisoaprs(kbdim,kwiso)

INTEGER :: kctop(kbdim), kcbot(kbdim)

LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: zrelq(kbdim,klev),                                         &
            zwisopsubcl(kbdim,kwiso)

REAL(dp) :: zdiagt,                                                    &
            zes,      zcor,       zqsat,                               &
            zwisorfl, zwisorfln,  zwisodrfl,                           &
            zt,       zqv,        zwisoqv,    zql,      zwisoql,       &
            zhumass,  zwisofracliq,  zqliq,      zqice, zquot,  zdelt, &
            zpromill, zdmaxo18,   zdmaxdeu,   zdabsmax,                &
            zqob,     zwisoqob,   zdeltaob,                            &
            zqun,     zwisoqun,   zdeltaun,                            &
            zqtop,    zwisoqtop,  zdeltatop,                           &
            zqmiddle, zwisoqmiddle,zdeltamiddle,                       &
            zqbottom, zwisoqbottom,zdeltabottom,                       &
            zqsum,    zwisoqsum,                                       &
            zdeltasum,                                                 &
            zwisosum_tmp,                                              &
            zwisorsum,zwisodpevap

INTEGER  :: jl, jk, jt, it, jkoben, jtf1, jtf2

INTEGER  :: kmix(kbdim), knull(kbdim) 

LOGICAL  :: lomix1, lomix2, lomix(kbdim), lonot

! SPECIFY PARAMETERS

  zdiagt=delta_time

  pwisorfl(:,:) = 0._dp
  pwisosfl(:,:) = 0._dp
  zwisopsubcl(:,:) = 0._dp

! calculation of the relative humidity below the cloud base

#ifndef NOLOOKUP
  lookupoverflow = .FALSE.
#endif

  DO jk=1,klev
    DO jl=1,kproma

      ptp1(jl,jk)=ptp1(jl,jk)+ptte(jl,jk)*time_step_len
      pqp1(jl,jk)=pqp1(jl,jk)+pqte(jl,jk)*time_step_len
       
      IF (jk.ge.kcbot(jl)) THEN

        it = NINT(ptp1(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/paphp1(jl,jk)
        zes=MIN(0.5_dp,zes)
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        IF (zqsat.lt.cwisomin) THEN
          zrelq(jl,jk)=0._dp
        ELSE IF (zqsat .LT. pqp1(jl,jk)) THEN
          zrelq(jl,jk)=1._dp
        ELSE
          zrelq(jl,jk)=pqp1(jl,jk)/zqsat
          IF ((1._dp-zrelq(jl,jk)).LT.cwisosec) zrelq(jl,jk)=1._dp
        ENDIF

      ELSE
        zrelq(jl,jk)=1._dp
      ENDIF
      
    END DO
  END DO

  IF (lookupoverflow) CALL lookuperror ('cuwisoequ')

! getting into equilibrium rain and vapour in and below the cloud

  DO jt=1,kwiso

! first part: calculate isotope precipitation like normal precipitation (done in *cuflx*)
 
    DO jk=ktopm2,klev
      DO jl=1,kproma
        IF (ldcum(jl)) THEN
          IF(pten(jl,jk).gt.tmelt) THEN
            pwisorfl(jl,jt)=pwisorfl(jl,jt)+pwisodmfup(jl,jk,jt)+pwisodmfdp(jl,jk,jt)
!           corrected rainfall for amount of melted water isotopes            
            pwisorfl(jl,jt)=pwisorfl(jl,jt)+pmelt_tmp(jl,jk)*pwisosfl(jl,jt)
            pwisosfl(jl,jt)=pwisosfl(jl,jt)-pmelt_tmp(jl,jk)*pwisosfl(jl,jt)
          ELSE
            pwisosfl(jl,jt)=pwisosfl(jl,jt)+pwisodmfup(jl,jk,jt)+pwisodmfdp(jl,jk,jt)
          END IF
        ENDIF
      END DO
    END DO

! check for negative rainfall and snowfall values
    DO jl=1,kproma
      if(prfl_tmp(jl).le.0.) pwisorfl(jl,jt)=0.
      if(psfl_tmp(jl).le.0.) pwisosfl(jl,jt)=0.
      zwisopsubcl(jl,jt)=pwisorfl(jl,jt)+pwisosfl(jl,jt)
    END DO
            
    DO jk=ktopm2,klev
      DO jl=1,kproma
        IF (ldcum(jl).AND.jk.ge.kcbot(jl).AND.ABS(zwisopsubcl(jl,jt)).gt.1.e-20_dp) THEN 

          zwisorfl=zwisopsubcl(jl,jt)
          zwisorfln=zwisorfl*pprec_frac(jl,jk)
          zwisodrfl=MIN(0._dp,zwisorfln-zwisorfl)
          pwisoqte(jl,jk,jt)=pwisoqte(jl,jk,jt)-zwisodrfl*g/(paphp1(jl,jk+1)-paphp1(jl,jk))
          zwisosum_tmp=pwisorfl(jl,jt)+pwisosfl(jl,jt)
          IF (ABS(zwisosum_tmp).LT.1.e-20_dp) THEN
            IF (zwisosum_tmp.GE.0._dp) THEN
              zwisosum_tmp=+1.e-20_dp
            ELSE 
              zwisosum_tmp=-1.e-20_dp
            ENDIF
          ENDIF
          pwisorfl(jl,jt)=zwisorfln*pwisorfl(jl,jt)/zwisosum_tmp
          pwisosfl(jl,jt)=zwisorfln*pwisosfl(jl,jt)/zwisosum_tmp

! second part: bring isotopes in precipitation in equilibrium with the surrounding water vapour  
! (only isotopes in rain water go into equilibrium)    
          IF (jk.lt.kcbot(jl).and.jk.ge.kctop(jl)) THEN ! if-condition = cloud layer
            zt=ptu(jl,jk)
            zqv=pqu(jl,jk)
            zwisoqv=pwisoqu(jl,jk,jt)
          ELSE                                          ! subcloud layer
            zt=ptp1(jl,jk)
            zqv=pqp1(jl,jk)
            zwisoqv=pwisoqp1(jl,jk,jt)+pwisoqte(jl,jk,jt)*time_step_len
          ENDIF
               
	      zql=prain_tmp(jl,jk)/((paphp1(jl,jk+1)-paphp1(jl,jk))/g/time_step_len)     
	      zwisoql=pwisorfl(jl,jt)/((paphp1(jl,jk+1)-paphp1(jl,jk))/g/time_step_len)
	        
!         calculate fractionation coefficient for liquid-vapour phase change
          zwisofracliq=exp(talphal1(jt)/(zt**2._dp)+talphal2(jt)/zt+talphal3(jt))

          IF (zqv.gt.cwisomin.and.zwisoqv.gt.cwisomin) THEN
	
! effective fractionation  below the cloud base
	        zhumass=thumwiso1+thumwiso2*zrelq(jl,jk)
	        zwisofracliq=zwisofracliq*zhumass/(zwisofracliq*(zhumass-1._dp)*(tdifrel(jt)**tdifexp)+1._dp)
	            
	        IF (zt.GT.tmelt) THEN   ! liquid and ice fraction
	          zqliq=1._dp
	        ELSEIF (zt.LT.twisoice) THEN
	          zqliq=0._dp
	        ELSE
	          zqice=(tmelt-zt)/(tmelt-twisoice)
	          zqliq=1._dp-zqice
	        ENDIF
	
	        zquot=zqv+zwisofracliq*twisoeqcu*zqliq*zql

            zdelt=0._dp
	        IF (zquot.gt.cwisomin) zdelt=twisoeqcu*zqliq*(zwisofracliq*zql*zwisoqv-zwisoql*zqv)/zquot

	        pwisorfl(jl,jt)=pwisorfl(jl,jt)+zdelt*((paphp1(jl,jk+1)-paphp1(jl,jk))/g/time_step_len)
	        pwisoqte(jl,jk,jt)=pwisoqte(jl,jk,jt)-zdelt/time_step_len

! for negative humidity values: set vapour in simple equilibrium to convective rain      
	      ELSE
	        
	        zdelt=tnat(jt)
		    IF (zql.gt.cwisomin) zdelt=zwisoql/zql
		        
		    zwisoqv=(1._dp/zwisofracliq)*zdelt*zqv
		    pwisoqte(jl,jk,jt)=(zwisoqv-pwisoqp1(jl,jk,jt))/time_step_len

	      ENDIF
            
          zwisopsubcl(jl,jt)=pwisorfl(jl,jt)+pwisosfl(jl,jt)

        ENDIF
        
      END DO
    END DO
  
  END DO

  DO jt=1,kwiso
    DO jl=1,kproma
      zwisorsum=pwisorfl(jl,jt)+pwisosfl(jl,jt)
      zwisodpevap=zwisopsubcl(jl,jt)-zwisorsum
      pwisorfl(jl,jt)=pwisorfl(jl,jt)+zwisodpevap*pwisorfl(jl,jt)*(1._dp/MAX(1.e-20_dp,zwisorsum))
      pwisosfl(jl,jt)=pwisosfl(jl,jt)+zwisodpevap*pwisosfl(jl,jt)*(1._dp/MAX(1.e-20_dp,zwisorsum))
    END DO
  END DO

! do additional vertical mixing of water isotopes
! if a certain vertical gradient of delta values is exceeded
! do mixing only if flag *lonot* is set
  lonot=.TRUE.
  IF (lonot) THEN
      
! set some constants 
    zpromill=1000._dp
!    zdmaxo18=1000._dp
!    zdmaxdeu=8000._dp  
    zdmaxo18=150._dp
    zdmaxdeu=1200._dp  
    zdabsmax=1000._dp
    jkoben=1
    kmix(:) = 0
    knull(:) = 0
    jtf1 = 0
    jtf2 = 0
      
! jtf1=tracerno. of h2-18o, jtf2=tracerno. of hdo 
    DO jt=1,kwiso
      IF (nwisotyp(jt).eq.2) jtf1=jt   
      IF (nwisotyp(jt).eq.3) jtf2=jt
    END DO

! loop over all levels, isotopes and longitudes
    DO jk=jkoben+1,klev-1
      lomix(:) = .FALSE.
      DO jt=1,kwiso
!CDIR$ IVDEP
        DO jl=1,kproma
        
! check, if delta-o18-gradient between two levels is gt. zdmaxo18
          IF (jtf1.gt.0) THEN
            zwisoqob=pwisoqp1(jl,jk-1,jtf1)+pwisoqte(jl,jk-1,jtf1)*time_step_len
            zqob=pqp1(jl,jk-1)
            zdeltaob=0._dp
            IF ((zwisoqob.gt.cwisomin).and.(zqob.gt.cwisomin)) zdeltaob=(zwisoqob/(zqob*tnat(jtf1))-1._dp)*zpromill
            zwisoqun=pwisoqp1(jl,jk,jtf1)+pwisoqte(jl,jk,jtf1)*time_step_len
            zqun=pqp1(jl,jk)
            zdeltaun=0._dp
            IF ((zwisoqun.gt.cwisomin).and.(zqun.gt.cwisomin)) zdeltaun=(zwisoqun/(zqun*tnat(jtf1))-1._dp)*zpromill
          ELSE
            zdeltaob=0._dp
            zdeltaun=0._dp
          ENDIF
          lomix1=abs(zdeltaob-zdeltaun).gt.zdmaxo18

! check, if delta-hdo-gradient between two levels is gt. zdmaxdeu
          IF (jtf2.gt.0) THEN
            zwisoqob=pwisoqp1(jl,jk-1,jtf2)+pwisoqte(jl,jk-1,jtf2)*time_step_len
            zqob=pqp1(jl,jk-1)
            zdeltaob=0._dp
            IF ((zwisoqob.gt.cwisomin).and.(zqob.gt.cwisomin)) zdeltaob=(zwisoqob/(zqob*tnat(jtf2))-1._dp)*zpromill
            zwisoqun=pwisoqp1(jl,jk,jtf2)+pwisoqte(jl,jk,jtf2)*time_step_len
            zqun=pqp1(jl,jk)
            zdeltaun=0._dp
            IF ((zwisoqun.gt.cwisomin).and.(zqun.gt.cwisomin)) zdeltaun=(zwisoqun/(zqun*tnat(jtf2))-1._dp)*zpromill
          ELSE
            zdeltaob=0._dp
            zdeltaun=0._dp
          ENDIF
          lomix2=abs(zdeltaob-zdeltaun).gt.zdmaxdeu
          
! check exactly once if one of the two gradients is too big
! if yes: mix three levels (mix only water isotopes, not normal water)    

          IF ((lomix1.or.lomix2).and.(jt.eq.1)) lomix(jl)=.TRUE.

          IF (lomix(jl).and.nwisotyp(jt).ne.1) THEN
           
! get old delta values of three boxes (lev-1,lev,lev+1)
! (neglect any box with negative water values...)           
            zwisoqtop=pwisoqp1(jl,jk-1,jt)+pwisoqte(jl,jk-1,jt)*time_step_len
            zqtop=pqp1(jl,jk-1)
            zdeltatop=0._dp
            IF ((zwisoqtop.gt.cwisomin).and.(zqtop.gt.cwisomin))               & 
               zdeltatop=(zwisoqtop/(zqtop*tnat(jt))-1._dp)*zpromill
          
            zwisoqmiddle=pwisoqp1(jl,jk,jt)+pwisoqte(jl,jk,jt)*time_step_len
            zqmiddle=pqp1(jl,jk)
            zdeltamiddle=0._dp
            IF ((zwisoqmiddle.gt.cwisomin).and.(zqmiddle.gt.cwisomin))         &
               zdeltamiddle=(zwisoqmiddle/(zqmiddle*tnat(jt))-1._dp)*zpromill

            zwisoqbottom=pwisoqp1(jl,jk+1,jt)+pwisoqte(jl,jk+1,jt)*time_step_len
            zqbottom=pqp1(jl,jk+1)
            zdeltabottom=0._dp
            IF ((zwisoqbottom.gt.cwisomin).and.(zqbottom.gt.cwisomin))         &
               zdeltabottom=(zwisoqbottom/(zqbottom*tnat(jt))-1._dp)*zpromill

! Calculate mean delta value of the three boxes (exclude boxes with too high delta-values)
            zqsum=0._dp
            zwisoqsum=0._dp
              
            IF (abs(zdeltatop).lt.zdabsmax) THEN 
              zqsum=zqsum+zqtop
              zwisoqsum=zwisoqsum+zwisoqtop
            ENDIF

            IF (abs(zdeltamiddle).lt.zdabsmax) THEN               
              zqsum=zqsum+zqmiddle
              zwisoqsum=zwisoqsum+zwisoqmiddle
            ENDIF
         
            IF (abs(zdeltabottom).lt.zdabsmax) THEN
              zqsum=zqsum+zqbottom
              zwisoqsum=zwisoqsum+zwisoqbottom
            ENDIF
        
            zdeltasum=0._dp
            IF ((zwisoqsum.gt.cwisomin).and.(zqsum.gt.cwisomin)) zdeltasum=(zwisoqsum/(zqsum*tnat(jt))-1._dp)*zpromill

! Calculate new isotope values of the three boxes
            zwisoqtop=(zdeltasum/zpromill+1)*zqtop*tnat(jt)
            zwisoqmiddle=(zdeltasum/zpromill+1)*zqmiddle*tnat(jt)
            zwisoqbottom=(zdeltasum/zpromill+1)*zqbottom*tnat(jt)

! Calculate new isotope tendencies
            pwisoqte(jl,jk-1,jt)=(zwisoqtop-pwisoqp1(jl,jk-1,jt))/time_step_len
            pwisoqte(jl,jk,jt)=(zwisoqmiddle-pwisoqp1(jl,jk,jt))/time_step_len
            pwisoqte(jl,jk+1,jt)=(zwisoqbottom-pwisoqp1(jl,jk+1,jt))/time_step_len

          ENDIF
          
        END DO
      END DO
    END DO
      
  ENDIF

! UPDATE SURFACE FIELDS
  DO jt=1,kwiso
    DO jl=1,kproma
      pwisorsfc(jl,jt)=pwisorfl(jl,jt)
      pwisossfc(jl,jt)=pwisosfl(jl,jt)
      pwisoaprc(jl,jt)=pwisoaprc(jl,jt)+zdiagt*(pwisorfl(jl,jt)+pwisosfl(jl,jt))
      pwisoaprs(jl,jt)=pwisoaprs(jl,jt)+zdiagt*pwisosfl(jl,jt)
    END DO
  END DO

  RETURN

END SUBROUTINE cuwisoequ
