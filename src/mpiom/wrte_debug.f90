      SUBROUTINE wrte_debug(III)

      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMOAU2
      USE MO_UNITS


      IMPLICIT NONE

      INTEGER*8 IDATE,ICODE,KCODE,IEDIM,IDATE2
      INTEGER k, l, III
!
      IUNIT=777
      IF(p_pe==p_io) THEN
        OPEN(IUNIT,FILE='D37000',FORM='UNFORMATTED')
      ENDIF

      IDATE=III
     

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
!:: WRITE DEPTO
      ICODE=84
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,DEPTO)

     END SUBROUTINE wrte_debug
