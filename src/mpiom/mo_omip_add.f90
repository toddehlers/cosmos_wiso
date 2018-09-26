MODULE mo_omip_add

  USE mo_kind
  USE mo_parallel
  USE mo_param1, ONLY: ie,je,ie_g,je_g 
  USE mo_control_add, ONLY: io_stdo_add,io_in_gpre16,io_in_gpre18,      &
                      io_in_gpreD,                                      & 
                      io_in_rval16,io_in_rval18,                        &
                      io_in_rvalD
  USE mo_contra, ONLY: pre16,pre18,preD,riv16,riv18,rivD

  USE mo_commo1, ONLY: lmonts,lmont1,monlen,ldays,lday,     &
                       lyear

  USE mo_commoau1

  IMPLICIT NONE

  
CONTAINS

  SUBROUTINE open_omip_add

    IMPLICIT NONE

    OPEN(io_in_gpre16,file='GIPREC16',status='unknown',                  &
         access='sequential',form='unformatted')

    OPEN(io_in_gpre18,file='GIPREC18',status='unknown',                  &
         access='sequential',form='unformatted')

    OPEN(io_in_gpreD,file='GIPRECD',status='unknown',                  &
         access='sequential',form='unformatted')

    OPEN(io_in_rval16,file='GIRIV16',status='unknown',                  &
         access='sequential',form='unformatted')

    OPEN(io_in_rval18,file='GIRIV18',status='unknown',                 &
         access='sequential',form='unformatted')

    OPEN(io_in_rvalD,file='GIRIVD',status='unknown',                 &
         access='sequential',form='unformatted')

  END SUBROUTINE open_omip_add



  SUBROUTINE spool_omip_add

    IMPLICIT NONE

    INTEGER(i8) :: ii1,ii2,ii3,ii4
    INTEGER     :: lmont  

    !spool the fields to the actual month

    IF ( lmont1 >= 1 ) THEN
       DO lmont=1,lmont1-1
          DO lday=1,monlen(lmont)

             WRITE(IO_STDO_ADD,*)'in spool'
             IF(p_pe==p_io) THEN
                READ(io_in_gpre16)ii1,ii2,ii3,ii4  
             ENDIF
             CALL read_slice(io_in_gpre16,pre16)

             IF(p_pe==p_io) THEN
                READ(io_in_gpre18)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_gpre18,pre18)

             IF(p_pe==p_io) THEN
                READ(io_in_gpreD)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_gpreD,preD)             

             IF(p_pe==p_io) THEN
                READ(io_in_rval16)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_rval16,riv16)

             IF(p_pe==p_io) THEN
                READ(io_in_rval18)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_rval18,riv18)

             IF(p_pe==p_io) THEN
                READ(io_in_rvalD)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_rvalD,rivD)

             WRITE(io_stdo_add,*)'spool: ',lmont,lday,ii1
          ENDDO
       ENDDO
       WRITE(io_stdo_add,*) 'forcing data is spooled to start of month ',lmont1
    ENDIF

  END SUBROUTINE spool_omip_add



  SUBROUTINE read_omip_add

    IMPLICIT NONE

    INTEGER(i8) :: ii1,ii2,ii3,ii4  

       IF(p_pe==p_io) THEN
          READ(io_in_gpre16)ii1,ii2,ii3,ii4  
       ENDIF
       CALL read_slice(io_in_gpre16,pre16)
       IF (p_pe==p_io) THEN
          WRITE(0,'(a,i4.4,''-'',i2.2,''-'',i2.2,a,i8)') &
               'ECHM: read forcing: ', lyear,lmonts,ldays, &
               ' from file:', ii1
       ENDIF

       IF(p_pe==p_io) THEN
          READ(io_in_gpre18)ii1,ii2,ii3,ii4  
       ENDIF
       CALL read_slice(io_in_gpre18,pre18)

       IF(p_pe==p_io) THEN
          READ(io_in_gpreD)ii1,ii2,ii3,ii4  
       ENDIF
       CALL read_slice(io_in_gpreD,preD)

       IF(p_pe==p_io) THEN
          READ(io_in_rval16)ii1,ii2,ii3,ii4  
       ENDIF
       CALL read_slice(io_in_rval16,riv16)

       IF(p_pe==p_io) THEN 
          READ(io_in_rval18)ii1,ii2,ii3,ii4  
       ENDIF
       CALL read_slice(io_in_rval18,riv18)

       IF(p_pe==p_io) THEN 
          READ(io_in_rvalD)ii1,ii2,ii3,ii4  
       ENDIF
       CALL read_slice(io_in_rvalD,rivD)

       !      For periodic boundaries

       
       CALL bounds_exch('p',PRE18,'mo_omip_add 1')
       CALL bounds_exch('p',PRE16,'mo_omip_add 2')
       CALL bounds_exch('p',PRED,'mo_omip_add 3')



  END SUBROUTINE read_omip_add


END MODULE mo_omip_add
