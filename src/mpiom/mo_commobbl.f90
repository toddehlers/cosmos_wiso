MODULE mo_commobbl

  USE mo_param1
  IMPLICIT NONE

  REAL, ALLOCATABLE:: alpx(:,:),alpy(:,:)                           &
       , bblflx(:,:),ubbl(:,:),vbbl(:,:)

  REAL, ALLOCATABLE :: woback(:,:,:),bbu(:,:),bbv(:,:)

  INTEGER, ALLOCATABLE :: kupubbl(:,:),kdwubbl(:,:)                 &
       ,kupvbbl(:,:),kdwvbbl(:,:)


CONTAINS

  SUBROUTINE alloc_mem_commobbl

    ALLOCATE(woback(ie,je,kep),bbu(ie,je),bbv(ie,je))

    ALLOCATE(alpx(ie,je),alpy(ie,je),bblflx(ie,je),                   &
         ubbl(ie,je),vbbl(ie,je),                             &
         kupubbl(ie,je),kdwubbl(ie,je),                       &
         kupvbbl(ie,je),kdwvbbl(ie,je) )

                alpx(:,:)=0.0
                alpy(:,:)=0.0
                bblflx(:,:)=0.0
                ubbl(:,:)=0.0
                vbbl(:,:)=0.0
                kupubbl(:,:)=0   
                kdwubbl(:,:)=0   
                kupvbbl(:,:)=0   
                kdwvbbl(:,:)=0
   
  END SUBROUTINE alloc_mem_commobbl



  SUBROUTINE SLOPETRANS

    !:: SBR TO CALCULATE BBL TRANSPORT VOR USE
    !:: IN TRACER ADVECTION                       
    !:: JHJ FEB.24. 2000    

    USE MO_PARAM1
    USE MO_PARALLEL
    USE MO_COMMO1



    INTEGER :: i,j,k,ivz,IDO,ISH,KUP,KDW,IMIST,JVZ,JDO,JSH 

    REAL :: GEE, REFDENS, YMUE,FLXFAC,BBLMAX,BBLMIN,ALP       &
            ,AALP,DDRO1,DZWBBL,UBOT,VBOT

    REAL :: aux(ie,je)



    GEE=9.81
    REFDENS=1025.0

! LINEAR FRICTION COEFF FOR CAMPIN&GOOSSE FLUX

    YMUE=1.E-6 
    FLXFAC=GEE/(REFDENS*YMUE)
    ! MAX BBL LAYER THICKNESS
    BBLMAX=500.
    BBLMIN=100.

    DO J=1,JE
       DO I=1,IE
          UBBL(I,J)=0.0       
          VBBL(I,J)=0.0        
          KUPUBBL(I,J)=0
          KUPVBBL(I,J)=0
          KDWUBBL(I,J)=0
          KDWVBBL(I,J)=0
          BBLFLX(I,J)=0.0
       ENDDO
    ENDDO

!:: SWEEP X

    DO J=2,JE-1
       DO I=2,IE-1

          ALP=ALPX(I,J)
          AALP=ABS(ALP)         
!:: SET ZONAL INDEX FOR SHELF (SH) AND DEEP OCEAN (DO)
          IVZ=NINT(ALP/(1.E-11+AALP))

!:: IF SIGN(SLOPE) IS POSITIVE (I.E. DEP(I+1)> DEP(I))
!:: THEN THE FLOW IS FROM ISH=I TO IDO=I+1 AND VICE VERSA

          IDO=I+MAX(0,IVZ)
          ISH=I+MAX(0,-IVZ)
          KUP=KBOT(ISH,J)
          IF (KUP.LT.1) THEN
             DDRO1=0.
          ELSE
             DDRO1=RHOO(ISH,J,KUP)-RHOO(IDO,J,KUP)
          ENDIF

!:: CHECK CONDITION FOR DOWNSLOPE CONV
 
          IF (DDRO1 .GT. 1.E-8) THEN

!:: SEARCH LEVEL OF EQUAL IN-SITU DENSITY
             KDW=KUP
             IMIST=0
             DO K=KUP+1,KBOT(IDO,J)
                IF ( (RHOO(ISH,J,K) .GT. RHOO(IDO,J,K))      &
                         .AND. (IMIST .EQ. 0) ) THEN
                   KDW=K
                ELSE
                   IMIST=1
                ENDIF
             ENDDO
!:: DO SLOPE CONV ONLY IF THERE IS A LESS DENSE DEEPER
!:: LAYER IN THE NEIGHBORING BOX, OR, IF BOTH CELLS ARE
!:: IN THE BOTTOMMOST LAYER BUT HAVE DIFFERENT DEPTHS:

             IF (KDW.GT.KUP) THEN

!:: SAVE UP AND DOWNLOCATION
                KUPUBBL(I,J)=KUP
                KDWUBBL(I,J)=KDW

!:SET THICKNESS OF DOWNFLOW TO BML THICKNESS OF BBL 
!:: OR TO THE TOTAL THICKNESS
                DZWBBL=MIN(BBLMAX,DDUO(I,J,KUP))

                
            
!:: USE BOTTOM VELOCITY AT POINT I,J ONLY IF
!   IT IS DIRECTED DOWNSLOPE

                UBOT=UKO(I,J,KUP)*DZWBBL
                IF (ALP*UBOT .GT. 0.0) THEN
                   UBBL(I,J)=UBOT*DLYU(I,J)
                ENDIF


             ENDIF
        
          ENDIF
        
       ENDDO   
    ENDDO
    
!!$OMP SINGLE 
    CALL bounds_exch('u',UBBL,'mo_commobbl 1')
!!$OMP END SINGLE 

    AUX(:,:) = REAL(KUPUBBL(:,:))
!!$OMP SINGLE 
   CALL bounds_exch('u',aux,'mo_commobbl 2')
!!$OMP END SINGLE 
    KUPUBBL(:,:) = INT(AUX(:,:))

    AUX(:,:) = REAL(KDWUBBL(:,:))
!!$OMP SINGLE 
    CALL bounds_exch('u',aux,'mo_commobbl 3')
!!$OMP END SINGLE 
    KDWUBBL(:,:) = INT(AUX(:,:))


!:: SWEEP Y

    DO J=2,JE-1
       DO I=2,IE-1

          ALP=ALPY(I,J)
          AALP=ABS(ALP) 

!:: SET MERIDIONAL INDEX FOR SHELF (SH) AND DEEP OCEAN (DO)
!:: IF SIGN(SLOPE) IS POSITIVE (I.E. DEP(J)> DEP(J+1))
!:: THEN THE FLOW IS FROM JSH=J+1 TO JDO=J AND VICE VERSA

          JVZ=NINT(ALP/(1.E-11+AALP)) 
          JDO=J+MAX(0,-JVZ)
          JSH=J+MAX(0,JVZ)
          KUP=KBOT(I,JSH)
          IF (KUP.LT.1) THEN
             DDRO1=0.
          ELSE
             DDRO1=RHOO(I,JSH,KUP)-RHOO(I,JDO,KUP)
          ENDIF

!:: CHECK CONDITION FOR DOWNSLOPE CONV

          IF (DDRO1 .GT. 1.E-8) THEN

!:: SEARCH LEVEL OF EQUAL IN-SITU DENSITY

             KDW=KUP
             IMIST=0
             DO K=KUP+1,KBOT(I,JDO)
                IF ((RHOO(I,JSH,K) .GT. RHOO(I,JDO,K))                           &
                     &      .AND. (IMIST .EQ. 0)) THEN
                   KDW=K
                ELSE
                   IMIST=1
                ENDIF
             ENDDO

!:: DO SLOPE CONV ONLY IF THERE IS A LESS DENSE DEEPER
!:: LAYER IN THE NEIGHBORING BOX, OR, IF BOTH CELLS ARE
!:: IN THE BOTTOMMOST LAYER BUT HAVE DIFFERENT DEPTHS:

             IF (KDW.GT.KUP) THEN

!:: SAVE UP AND DOWNLOCATION
                KUPVBBL(I,J)=KUP
                KDWVBBL(I,J)=KDW
!:SET THICKNESS OF DOWNFLOW TO BML THICKNESS OF BBL M
!:: OR TO THE TOTAL THICKNESS
                DZWBBL=MIN(BBLMAX,DDUE(I,J,KUP))

!:: USE BOTTOM VELOCITY AT POINT I,J
! IT IS DIRECTED DOWNSLOPE
                VBOT=VKE(I,J,KUP)*DZWBBL
                IF (ALP*VBOT .GT. 0.0) THEN
                   VBBL(I,J)=VBOT*DLXV(I,J)
                ENDIF
             ENDIF
          ENDIF
       ENDDO
    ENDDO

!!$OMP SINGLE 
    CALL bounds_exch('v',VBBL,'mo_commobbl 4')
!!$OMP END SINGLE 
    
    AUX(:,:) = REAL(KUPVBBL(:,:))

!!$OMP SINGLE
    CALL bounds_exch('u',aux,'mo_commobbl 5')
!!$OMP END SINGLE 

    KUPVBBL(:,:) = INT(AUX(:,:))

    AUX(:,:) = REAL(KDWVBBL(:,:))

!!$OMP SINGLE 
    CALL bounds_exch('u',aux,'mo_commobbl 6')
!!$OMP END SINGLE 

    KDWVBBL(:,:) = INT(AUX(:,:))



    !:: SAVE ABS OF TRANSPORT FOR DIAGNOSTIC

    DO J=1,JE
       DO I=1,IE
          BBLFLX(I,J)=SQRT(UBBL(I,J)**2+VBBL(I,J)**2)*1.E-6
       ENDDO
    ENDDO


  END SUBROUTINE SLOPETRANS


END MODULE mo_commobbl
