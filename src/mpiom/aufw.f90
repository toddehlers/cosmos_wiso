      SUBROUTINE AUFW
!****************************************************************
!
!     AAAAAA  U    U  FFFFFF  W    W
!     A    A  U    U  F       W    W
!     AAAAAA  U    U  FFFFF   W WW W
!     A    A  U    U  F       WWWWWW
!     A    A  UUUUUU  F       WW  WW
!
!
!*****************************************************************
! SUBROUTINE AUFW
!
!     THIS SBR WRITES ROUTINELY AT CERTAIN TIME INTERVALS
!     A RESTART FILE ON Z37000 OR Z38000
!
!-----------------------------------------------------------------
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMOAU2
      USE MO_UNITS

#ifdef __coupled
      USE MO_FLUXES1
#endif

      IMPLICIT NONE
      INTEGER*8 IDATE,ICODE,KCODE,IEDIM,IDATE2
      INTEGER k, l
      REAL rldt, rldays
!
      IF(p_pe==p_io) THEN
        OPEN(IO_IN_Z370,FILE='Z37000',STATUS='UNKNOWN'                  &
     &                 ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

        OPEN(IO_IN_Z380,FILE='Z38000',STATUS='UNKNOWN'                  &
     &                 ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

        REWIND(IO_IN_Z370)
        REWIND(IO_IN_Z380)
      ENDIF
!
      IUNIT=IUNIT+IFLAG
      WRITE(IO_STDOUT,*) 'AUFW: IUNIT IFLAG',IUNIT,IFLAG
      IF(p_pe==p_io) REWIND(IUNIT)
      WRITE(IO_STDOUT,*)' ++++++ WRITING RESTART FILE ON ',IUNIT        &
     &,' AFTER ',LDT                                                    &
     &,' TIME STEPS . CALCULATED YEAR : ',LYEARS,' MONTH : ',LMONTS
!:: WRITE ALL DATA IN EXTRA FORMAT
      IDATE= (LYEARS*10000)+(LMONTS*100)+LDAYS
      RLDT=REAL(LDT)
      RLDAYS=REAL(LDAYS)
!:: WRITE TIME STEP INFORMATION
      KCODE=1
      IEDIM=2
      ICODE =999
      IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
      IF(p_pe==p_io) WRITE(IUNIT)RLDT,RLDAYS                   
!  WRITE ZONAL VELOCITY COMPONENT
      ICODE=3
      IEDIM=IE_G*JE_G
      DO K=1,KE
       KCODE=TIESTU(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,UKO(:,:,K))
      ENDDO
!  WRITE MERIDIONAL  VELOCITY COMPONENT
      ICODE=4
      DO K=1,KE
       KCODE=TIESTU(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,VKE(:,:,K))
      ENDDO
!  WRITE TEMPERATURE                     
      ICODE=2
      DO K=1,KE
       KCODE=TIESTU(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,THO(:,:,K))
      ENDDO
!  WRITE SALINITY                        
      ICODE=5
      DO K=1,KE
       KCODE=TIESTU(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,SAO(:,:,K))
      ENDDO
!:: WRITE 2-D FIELDS
      KCODE=0
!  WRITE SEA SURFACE ELEVATION           
      ICODE=1
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,ZO)
!  WRITE SEA SURFACE ELEVATION CHANGE    
      ICODE=82
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,Z1O)
!  WRITE SEA ICE THICKNESS               
      ICODE=13
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,SICTHO)
!  WRITE SEA ICE CONCENTRATION           
      ICODE=15
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,SICOMO)
!  WRITE ZONAL SEA ICE VELOCITY                
      ICODE=35
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,SICUO)
!  WRITE MERIDIONAL SEA ICE VELOCITY                
      ICODE=36
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,SICVE)
!  WRITE SNOW                                       
      ICODE=141
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,SICSNO)
!  WRITE HIBLER ETA/ZETA FIELDS                     
      ICODE=501 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,HIBETE)
      ICODE=502 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,HIBETO)
      ICODE=503 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,HIBZETE)
      ICODE=504 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,HIBZETO)
!  WRITE VERTICAL DIFFUSIVITY DVO                            
      ICODE=111
      DO K=1,KEP
       KCODE=TIESTW(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,DVO(:,:,K))
      ENDDO
!  WRITE VERTICAL FRICTION AVO
      ICODE=110
      DO K=1,KEP
       KCODE=TIESTW(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,AVO(:,:,K))
      ENDDO
!  WRITE VERTICAL VELOCITY (DIAGNOSTIC)
      ICODE=7
      DO K=1,KEP
       KCODE=TIESTW(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,WO(:,:,K))
      ENDDO
      ICODE=99
      KCODE=0
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,TICE)
#ifdef MIST
!:: WRITE 2-D DIAGNOSTIC FIELDS
      KCODE=0
!:: WRITE ZONAL ICE TRANSPORTS
      ICODE=142
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,EISTRX)
!:: WRITE MERIDIONAL ICE TRANSPORTS
      ICODE=143
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,EISTRY)
!:: WRITE MAX. CONVECTION DEPTH    
      ICODE=69 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,FLOAT(KCONDEP(:,:)))
!:: WRITE HEAT FLUX                
      ICODE=70 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,HFLM)
!:: WRITE P-E                      
      ICODE=79 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,PMEM)
!:: WRITE ZONAL WIND STRESS              
      ICODE=52 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
#ifdef __coupled
       CALL write_slice(IUNIT,AOFLTXWO)
#else
       CALL write_slice(IUNIT,TXO)
#endif /*__coupled*/
!:: WRITE MERIDIONAL WIND STRESS              
      ICODE=53 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
#ifdef __coupled
       CALL write_slice(IUNIT,AOFLTYWE)
#else
       CALL write_slice(IUNIT,TYE)
#endif /*__coupled*/
!:: WRITE SHORT WAVE RAD FLUX      
      ICODE=176
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,QSWO)
!:: WRITE LONG WAVE RAD FLUX      
      ICODE=177
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,QLWO)
!:: WRITE SENS HEAT FLUX          
      ICODE=146
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,QSEO)
!:: WRITE LATENT HEAT FLUX          
      ICODE=147
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,QLAO)
!:: WRITE STREAM FUNCTION           
      ICODE=27 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,PSIUWE)
!:: WRITE MONTHLY MIXED LAYER DEPTH         
      ICODE=83
      DO L=1,12
       IDATE2= (LYEARS*10000)+(L*100)+LDAYS
       IF(p_pe==p_io) WRITE(IUNIT) IDATE2,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,AMLD(:,:,L))
      ENDDO
#endif
!:: WRITE DEPTO
      ICODE=84
      kcode=0
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,DEPTO)



#ifdef MIST
!:: WRITE DLXP
      ICODE=85
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,DLXP)
!:: WRITE DLYP
      ICODE=86
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,DLYP)
!:: WRITE WETO
      ICODE=506
      DO K=1,KE
       KCODE=TIESTU(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,WETO(:,:,K))
      ENDDO
#endif 
      IFLAG=-IFLAG
      IF(p_pe==p_io) THEN
        REWIND(IUNIT)
        CLOSE(IO_IN_Z370)
        CLOSE(IO_IN_Z380)
      ENDIF
      RETURN
      END
