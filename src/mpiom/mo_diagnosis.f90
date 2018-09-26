MODULE mo_diagnosis

  USE mo_param1
  USE mo_commo1
  USE mo_mpi
  USE mo_parallel

#ifdef PBGC
  USE mo_carbch, ONLY: atm_co2
#endif

  IMPLICIT NONE

  INTEGER, POINTER :: ibek(:,:), ibek_g(:,:)
  REAL, POINTER :: sigh(:,:),zmld(:,:)

  REAL,ALLOCATABLE :: psiuwe_g(:,:)

  INTEGER :: grid_family

  
  INTEGER :: nmecen
  PARAMETER (nmecen=5)

  REAL,ALLOCATABLE ::                                               &
       tberi(:),tfaroe(:),tdenmar(:)                                &
       ,tdavis(:),tmedi(:),hflb(:)                                   &
       ,wflb(:),eiscb(:),eisab(:)                                    &
       ,arem(:,:),tquer(:,:),squer(:,:)                              &
!---wiso-code
       ,o16quer(:,:),o18quer(:,:),hdoquer(:,:)                       &
!---wiso-code-end       
       ,dvquer(:,:),avquer(:,:)                                      &
       ,tmerci(:,:)                                                  &
       ,aabw(:),rnadw(:)                                             &
       ,tvquer(:,:),svquer(:,:)                                      &
       ,tdquer(:,:)


  INTEGER iban1,jban1,iban2,jban2                                   &
       ,idra1,jdra1,idra2,jdra2                                     &
       ,iber1,jber1,iber2,jber2                                     &
       ,idberi,jdberi,iddenm,jddene,jddenw                          &
       ,idmedi,jdmedi                                               &
       ,idfarw,idfare,jdfaro,kdfaro                                 &
       ,iddavi,jddavw,jddave                                        &
       ,idfram,jdfraw,jdfrae                                        &
       ,iddav1,jddav1,iddav2,jddav2                                 &
       ,idguln,idguls,jdgulw,jdgule                                 &
       ,idkurn,idkurs,jdkurw,jdkure                                 &
       ,jdgss,isued,inord                                           &
       ,kaabw,knadw,jlink(nmecen),jrech(nmecen),imerci(nmecen)      &
       ,iddene,iddenw,jddenm                                        &
       ,jddenn,jddens                                               &
       ,jdfarn,jdfars,idfaro                                        &
       ,idgulw,idgule,jdguln,jdguls                                 &
       ,idkurw,idkure,jdkurn,jdkurs                                 &
       ,ilink(nmecen),irech(nmecen),jmerci(nmecen)


    REAL, ALLOCATABLE :: tmerc(:,:,:),sum_tmerc(:,:,:)
    real, allocatable :: avg_timeser(:)


!FR   Must be global variables:
    REAL    :: psigolf,psikuro,psiaspg,cflux,sfram
    INTEGER :: ltq1,ltq2,ltq3,ltq4
    INTEGER :: ldenmar,lfaroe


CONTAINS

  SUBROUTINE alloc_mem_diag

    ! define her the grid families, i.e. take into account grid orientation
    ! for throughflows are different in gin and grob-type grids

    grid_family=1

    ! versiongin
    IF (ie_g.EQ.182.AND.je_g.EQ.84) THEN
       grid_family=2
    ENDIF


    ALLOCATE(tberi(ke),tfaroe(ke),tdenmar(ke),tdavis(ke),tmedi(ke)    &
         ,arem(ke,nbox),tquer(ke,nbox),squer(ke,nbox),dvquer(ke,nbox) &
!---wiso-code
         ,o16quer(ke,nbox),o18quer(ke,nbox),hdoquer(ke,nbox)          &
!---wiso-code-end       
         ,avquer(ke,nbox),tmerci(ke,nmecen),tvquer(ke,nmecen)         &
         ,svquer(ke,nmecen),tdquer(ke,nmecen),aabw(nmecen),rnadw(nmecen) &
         ,wflb(nbox),hflb(nbox),eiscb(nbox),eisab(nbox))

    ALLOCATE(ibek(ie,je),sigh(ie,je),zmld(ie,je))
    ALLOCATE(ibek_g(ie_g,je_g),psiuwe_g(ie_g,je_g))
    
    ALLOCATE(tmerc(180,2,kep),sum_tmerc(180,2,kep))
!---wiso-code
!    ALLOCATE(avg_timeser(200))
    ALLOCATE(avg_timeser(300))
!---wiso-code-end       

  END SUBROUTINE alloc_mem_diag

  SUBROUTINE diag_ini  

    USE mo_param1
    USE mo_mpi
    USE mo_commo1
    USE mo_commoau1
    USE mo_commoau2
    USE mo_units


    INTEGER i,j,k,n
    REAL dist

    ! initialisations
    !:: common diagnostics
    DO j=1,je
       DO i=1,ie
          ibek(i,j)=0
       ENDDO
    ENDDO


#ifdef DIAG
    ! open file for timeseries
    IF(p_pe==p_io) then
       OPEN(io_ou_f125,file='TIMESER.asc')
       OPEN(io_ou_f126,file='TIMESER.ext',form='unformatted')
    endif
    !:: initialize diag fields
    DO k=1,ke
       tberi(k)=zero
       tfaroe(k)=zero
       tdenmar(k)=zero
       tdavis(k)=zero
       tmedi(k)=zero
    ENDDO
    DO n=1,nbox
       hflb(n)=zero
       wflb(n)=zero
       eiscb(n)=zero
       eisab(n)=zero
       DO k=1,ke
          arem(k,n)=zero
          tquer(k,n)=zero
!---wiso-code
          o16quer(k,n)=zero
          o18quer(k,n)=zero
          hdoquer(k,n)=zero
!---wiso-code-end       
          squer(k,n)=zero
          dvquer(k,n)=zero
          avquer(k,n)=zero
       ENDDO
    ENDDO



    DO n=1,nmecen

       IF (grid_family.EQ. 2) THEN
          jlink(n)=0
          jrech(n)=0
          imerci(n)=0
       ELSE
          ilink(n)=0
          irech(n)=0
          jmerci(n)=0
       ENDIF

       aabw(n)=zero
       rnadw(n)=zero

       DO k=1,ke
          tmerci(k,n)=zero
          tvquer(k,n)=zero
          svquer(k,n)=zero
          tdquer(k,n)=zero
       ENDDO

    ENDDO


    IF(ie_g.EQ.182.AND.je_g.EQ.84)THEN           ! versiongin

       imerci(1)=50
       imerci(2)=123
       imerci(3)=134
       imerci(4)=148
       imerci(5)=172
       jlink(1)=64
       jlink(2)=22
       jlink(3)=20
       jlink(4)=35
       jlink(5)=35
       jrech(1)=64
       jrech(2)=55
       jrech(3)=65
       jrech(4)=72
       jrech(5)=72     

    ENDIF

    IF(ie_g.EQ.130.AND.je_g.EQ.211)THEN         ! versiont43

       jmerci(1)=55
       jmerci(2)=48 
       jmerci(3)=55 
       jmerci(4)=69 
       jmerci(5)=145
       ilink(1)=7 
       ilink(2)=42
       ilink(3)=50
       ilink(4)=42
       ilink(5)=50
       irech(1)=13
       irech(2)=88
       irech(3)=76
       irech(4)=72
       irech(5)=85

    ENDIF

    IF(ie_g.EQ.60.AND.je_g.EQ.50)THEN         ! versiongr60
       jmerci(1)=19
       jmerci(2)=14 
       jmerci(3)=17 
       jmerci(4)=21 
       jmerci(5)=34 
       ilink(1)=1 
       ilink(2)=21
       ilink(3)=19
       ilink(4)=17
       ilink(5)=22
       irech(1)=3
       irech(2)=36
       irech(3)=36
       irech(4)=37
       irech(5)=36
    ENDIF


    IF(ie_g.EQ.122.AND.je_g.EQ.101)THEN         ! versiongr30

       jmerci(1)=41
       jmerci(2)=12 
       jmerci(3)=33 
       jmerci(4)=41 
       jmerci(5)=71 
       ilink(1)=120 
       ilink(2)=42
       ilink(3)=44
       ilink(4)=32
       ilink(5)=43
       irech(1)=122
       irech(2)=72
       irech(3)=74
       irech(4)=70
       irech(5)=75
    ENDIF


    IF(ie_g.EQ.256.AND.je_g.EQ.220)THEN         ! versiongr15          
       jmerci(1)=85
       jmerci(2)=34
       jmerci(3)=70
       jmerci(4)=94
       jmerci(5)=155
       ilink(1)=1
       ilink(2)=107
       ilink(3)=78
       ilink(4)=60
       ilink(5)=87
       irech(1)=15
       irech(2)=153
       irech(3)=157
       irech(4)=150
       irech(5)=160
    ENDIF

    IF(ie_g.EQ.400.AND.je_g.EQ.338)THEN         ! versiongr09
       jmerci(1)=135
       jmerci(2)=40
       jmerci(3)=109
       jmerci(4)=123
       jmerci(5)=235
       ilink(1)=398
       ilink(2)=139
       ilink(3)=145
       ilink(4)=105
       ilink(5)=142
       irech(1)=400
       irech(2)=237
       irech(3)=244
       irech(4)=231
       irech(5)=248
    ENDIF


    IF(ie_g.EQ.182.AND.je_g.EQ.84)THEN           ! versiongin

       !        bering
       idberi=50
       jdberi=63

       !        denmark
       iddenm=115
       jddene=40
       jddenw=50

       !        mediterranean outflow
       idmedi=146
       jdmedi=33

       !        faroer-bank
       jdfaro=33
       idfarw=121
       idfare=124
       kdfaro=10

       IF (ke.EQ.40)   kdfaro=15

       IF (ke.EQ.30)   kdfaro=17

       !        davis
       iddavi=92
       jddavw=67
       jddave=69

       !        ice fram
       idfram=92
       jdfraw=33
       jdfrae=45

       !        ice davis
       iddav1=15
       jddav1=64
       iddav2=36
       jddav2=70

       !        gulfstream: max of streamfunction
       idguln=130
       idguls=158
       jdgulw=62
       jdgule=84

       !        kuroshio: max of streamfunction
       idkurn=22
       idkurs=50
       jdkurw=36
       jdkure=53

       !         find gulfstream separation   
       jdgss=61
       isued=158
       inord=135

       !        aabw = max below level kaabw
       kaabw=11
       IF (ke.EQ.40)  kaabw=18
       IF (ke.EQ.30)  kaabw=16

       !        nadw = min above level knadw 
       knadw=9
       IF (ke.EQ.40) knadw=14
       IF (ke.EQ.30) knadw=13
       IF (ke.EQ.23) knadw=12

    ENDIF           ! versiongin



    IF (ie_g .EQ. 60 .AND. je_g .EQ. 50 ) THEN   ! versiongr60
       !        bering
       idberi=2
       jdberi=19
       !
       !        denmark
       iddenm=34
       jddenn=1
       jddens=6
       !
       jddenm=1
       iddenw=34
       iddene=34
       !:: mediterranean outflow
       idmedi=34    
       jdmedi=21
       !
       !        faroer-bank
       idfaro=36
       jdfarn=8
       jdfars=14
       kdfaro=12

       IF ( ke .EQ. 40 ) kdfaro=15
       IF ( ke .EQ. 30 ) kdfaro=17
       IF ( ke .EQ. 23 ) kdfaro=17

       !        davis
       iddavi=16
       jddavw=1
       jddave=9
       !
       !        ice fram
       idfram=52
       jdfraw=3
       jdfrae=11
       !
       !        ice davis
       iddav1=16
       jddav1=1
       iddav2=16
       jddav2=9
       !      gulfstream: max of streamfunction
       idgulw=35
       idgule=45
       jdguln=33
       jdguls=50

       !
       !      kuroshio: max of streamfunction
       idkurw=105
       idkure=120
       jdkurn=50
       jdkurs=60
       !
       !      find gulfstream separation   
       jdgss=45
       isued=45
       inord=35
       !
       !      aabw = max below level kaabw
       kaabw=11

       IF ( ke .EQ. 40 ) kaabw=18
       IF ( ke .EQ. 30 ) kaabw=16
       IF ( ke .EQ. 23 ) kaabw=18



       !      nadw = min above level knadw 
       knadw=9

       IF ( ke .EQ. 40 ) knadw=14
       IF ( ke .EQ. 30 ) knadw=13
       IF ( ke .EQ. 23 ) knadw=12


    ENDIF

    IF (ie_g .EQ. 122 .AND. je_g .EQ. 101 ) THEN   ! versiongr30

       !        bering
       idberi=2
       jdberi=41
       !
       !        denmark
       iddenm=69
       jddenn=2
       jddens=16
       !
       jddenm=2
       iddenw=48
       iddene=72
       !:: mediterranean outflow
       idmedi=69
       jdmedi=43
       !
       !        faroer-bank
       idfaro=73
       jdfarn=15
       jdfars=30
       kdfaro=12

       IF ( ke .EQ. 40 ) kdfaro=15
       IF ( ke .EQ. 30 ) kdfaro=17
       IF ( ke .EQ. 23 ) kdfaro=17


       !        davis
       iddavi=11
       jddavw=20
       jddave=29
       !
       !        ice fram
       idfram=101
       jdfraw=3
       jdfrae=23
       !
       !        ice davis
       iddav1=13
       jddav1=31
       iddav2=15
       jddav2=33
       !      gulfstream: max of streamfunction
       idgulw=35
       idgule=45
       jdguln=33
       jdguls=50

       !
       !      kuroshio: max of streamfunction
       idkurw=105
       idkure=120
       jdkurn=50
       jdkurs=60
       !
       !      find gulfstream separation   
       jdgss=45
       isued=45
       inord=35
       !
       !      aabw = max below level kaabw
       kaabw=11


       IF ( ke .EQ. 40 ) kaabw=18
       IF ( ke .EQ. 30 ) kaabw=16
       IF ( ke .EQ. 23 ) kaabw=18


       !      nadw = min above level knadw 
       knadw=9

       IF ( ke .EQ. 40 ) knadw=14
       IF ( ke .EQ. 30 ) knadw=13
       IF ( ke .EQ. 23 ) knadw=12

    ENDIF


    IF(ie_g.EQ.256.AND.je_g.EQ.220)THEN         ! versiongr15          

       ! bering
       idberi=8
       jdberi=85
       ! denmark
       iddenm=152
       jddenn=2
       jddens=34

       jddenm=2
       iddenw=48
       iddene=72

       ! mediterranean outflow
       idmedi=147
       jdmedi=92

       ! faroer-bank
       idfaro=157
       jdfarn=35
       jdfars=65
       kdfaro=12

       IF (ke.EQ.40) kdfaro=14
       IF (ke.EQ.30) kdfaro=17
       IF (ke.EQ.23) kdfaro=17

       ! davis
       iddavi=22
       jddavw=49
       jddave=80

       ! ice fram
       idfram=214
       jdfraw=3
       jdfrae=50

       ! ice davis
       iddav1=28
       jddav1=42
       iddav2=29
       jddav2=424

       ! gulfstream: max of streamfunction
       idgulw=74
       idgule=96
       jdguln=72
       jdguln=72
       jdguls=98

       ! kuroshio: max of streamfunction
       idkurw=105
       idkure=120
       jdkurn=230
       jdkurs=250

       ! find gulfstream separation
       jdgss=90
       isued=74
       inord=84

       ! aabw = max below level kaabw
       kaabw=11

       IF (ke.EQ.40) kaabw=18
       IF (ke.EQ.30) kaabw=16
       IF (ke.EQ.23) kaabw=18

       ! nadw = min above level knadw
       knadw=9

       IF (ke.EQ.40) knadw=14
       IF (ke.EQ.30) knadw=13
       IF (ke.EQ.23) knadw=12

    ENDIF

    IF(ie_g.EQ.400.AND.je_g.EQ.338)THEN         ! versiongr09          

       !        bering
       idberi=6
       jdberi=135
       !
       !        denmark
       iddenm=228
       jddenn=7
       jddens=52
       !
       jddenm=6
       iddenw=159
       iddene=238
       !:: mediterranean outflow
       idmedi=228
       jdmedi=142
       !
       !        faroer-bank
       idfaro=241
       jdfarn=50
       jdfars=99
       kdfaro=40

       IF (ke.EQ.40) kdfaro=15
       IF (ke.EQ.30) kdfaro=17
       IF (ke.EQ.23) kdfaro=17

       !        davis
       iddavi=37
       jddavw=66
       jddave=96
       !        ice fram
       idfram=333
       jdfraw=10
       jdfrae=76
       !
       !        ice davis
       iddav1=43
       jddav1=102
       iddav2=50
       jddav2=109
       !      gulfstream: max of streamfunction
       idgulw=35
       idgule=45
       jdguln=33
       jdguls=50
       !
       !      kuroshio: max of streamfunction
       idkurw=105
       idkure=120
       jdkurn=50
       jdkurs=60
       !
       !      find gulfstream separation
       jdgss=45
       isued=45
       inord=35
       !
       !      aabw = max below level kaabw
       kaabw=11

       IF (ke.EQ.40) kaabw=18
       IF (ke.EQ.30) kaabw=16
       IF (ke.EQ.23) kaabw=18

       !      nadw = min above level knadw
       knadw=9
       IF (ke.EQ.40) knadw=14
       IF (ke.EQ.30) knadw=13
       IF (ke.EQ.23) knadw=12

    ENDIF
    !

    IF(ie_g.EQ.130.AND.je_g.EQ.211)THEN         ! versiont43

       !         bering
       idberi=10
       jdberi=55

       !         denmark strait
       jddenm=40
       iddene=67
       iddenw=68

       iddenm=67
       jddenn=34
       jddens=40

       !         mediterranean
       idmedi=75
       jdmedi=65

       !         faroe bank
       idfaro=76
       jdfarn=45
       jdfars=50

       !         fram strait     
       idfram=101
       jdfraw=7
       jdfrae=25

       !         davis strait
       iddavi=28 
       jddavw=39
       jddave=39

       !         gulfstream: max of streamfunction
       idgulw=40 
       idgule=65 
       jdguln=50
       jdguls=80

       !         kuroshio: max of streamfunction
       idkurw=115
       idkure=125
       jdkurn=60
       jdkurs=80

       !         find gulfstream separation   
       jdgss=61
       isued=158
       inord=135

       !         aabw = max below level kaabw
       kaabw=11
       !         nadw = min above level knadw 
       knadw=9

       IF (ke.EQ.40) THEN
          kaabw=33
          knadw=23
          kdfaro=27
       ENDIF

       IF (ke.EQ.23) THEN
          kaabw=81
          knadw=12
       ENDIF

    ENDIF


#endif /*diag*/


    !:: search for streamfunction points
    CALL suchij(-20.,-60.,1,idra2,jdra2,dist,0.)
    CALL suchij(-89.,-99.,1,idra1,jdra1,dist,0.)
    CALL suchij(30.,100.,1,iban1,jban1,dist,0.)
    CALL suchij(-25.,135.,1,iban2,jban2,dist,0.)
    CALL suchij(67.,160.,1,iber1,jber1,dist,0.)
    CALL suchij(65.,-150.,1,iber2,jber2,dist,0.)
    !
    WRITE(io_stdout,*)'drake: ',idra1,jdra1,weto_g(idra1,jdra1,1)    &
         &                    ,idra2,jdra2,weto_g(idra2,jdra2,1)
    WRITE(io_stdout,*)'banda: ',iban1,jban1,weto_g(iban1,jban1,1)    &
         &                    ,iban2,jban2,weto_g(iban2,jban2,1)
    WRITE(io_stdout,*)'bering: ',iber1,jber1,weto_g(iber1,jber1,1)   &
         &                    ,iber2,jber2,weto_g(iber2,jber2,1)
    !::

  END SUBROUTINE diag_ini

  SUBROUTINE calc_mixedlayerdepth

    USE mo_param1
    USE mo_mpi
    USE mo_parallel
    USE mo_commo1
    USE mo_commoau1
    USE mo_commoau2
    USE mo_units

    REAL zzzuwe,zzzcou,zzz,sigcrit

    INTEGER i,j,k,m

    !   compute mixed-layer depth

!$OMP PARALLEL PRIVATE(i,j,k,m)

    sigcrit=0.125
    m=lmonts


    sigh(:,:)=zero
    zmld(:,:)=zero

    !$OMP DO
    DO j=1,je       
       DO i=1,ie
          sigh(i,j)=weto(i,j,1)*sigcrit
          zmld(i,j)=weto(i,j,1)*tiestu(1)
       ENDDO
    ENDDO
    !$OMP END DO

    DO k=2,ke
       !$OMP DO
       DO j=1,je
          DO i=1,ie
             IF(weto(i,j,k)*sigh(i,j).GT.1.e-6)THEN
                zzz=MIN(sigh(i,j)/(ABS(stabio(i,j,k))+almzer),dz(k))
                sigh(i,j)=MAX(0.,sigh(i,j)-zzz*stabio(i,j,k))
                zmld(i,j)=zmld(i,j)+zzz
             ELSE
                sigh(i,j)=0.
             ENDIF
          ENDDO
       ENDDO
       !$OMP END DO
    ENDDO
    !
    DO j=2,je-1
       DO i=2,ie-1
          amld(i,j,m)=MAX(amld(i,j,m),zmld(i,j))
       ENDDO
    ENDDO

    zzzuwe = SUM(zmld(2:ie-1,2:je-1)*weto(2:ie-1,2:je-1,1))
    zzzcou = SUM(weto(2:ie-1,2:je-1,1))
!$OMP SINGLE
    CALL global_sum(zzzuwe,zzzcou)
!$OMP END SINGLE
    zzzcou = zzzcou + almzer

    WRITE(io_stdout,*)'mean mld: ',zzzuwe/zzzcou
    !hh      end mixed layer depths

!$OMP END PARALLEL

  END SUBROUTINE calc_mixedlayerdepth



  SUBROUTINE calc_icecutoff
  

    USE mo_param1
    USE mo_mpi
    USE mo_parallel
    USE mo_commo1
    USE mo_commoau1
    USE mo_commoau2
    USE mo_units
#ifdef PBGC
    USE mo_carbch, ONLY: ocetra
    USE mo_param1_bgc, ONLY: nocetra
    INTEGER l
#endif

    REAL swmin,eismax,volice,areice,volsno,arges,zmi,zma  & 
         ,schwell,ddice,dlxpdlyp

    REAL svold,svnew,draftold,draftnew

    INTEGER i,j,icou,ierr
    
    swmin=0.
    icou=0
    eismax=0.
    volice=0.
    areice=0.
    volsno=0.
    arges=0.
    zmi=0.
    zma=0.

    ierr=0


    !uwe     ensure that effective thickness of uppermost waterlayer
    !        does not fall below a critical value
    !        in that case ice is converted to water
    !        not heat conserving!!!!
    !         

    DO j=1,je
       DO i=1,ie
          IF(weto(i,j,1).GT.0.5)THEN
             schwell=0.7*(dzw(1)+zo(i,j))

             IF(rhoicwa*sictho(i,j)+rhosnwa*sicsno(i,j).GT.schwell)THEN
                !-et
                WRITE(io_stdout,*)'ice cut-off at ',i,j,sictho(i,j),sao(i,j,1)
                ierr=1

                draftold=ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa
                svold=sao(i,j,1)*draftold+sice*sictho(i,j)

                ddice=MIN((sictho(i,j)*rhoicwa+sicsno(i,j)*rhosnwa     &
                     -schwell)/rhoicwa,sictho(i,j))
                
                sictho(i,j)=sictho(i,j)-ddice
                
                draftnew=ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa

                !uwe       salt conservation!
!                sao(i,j,1)=(ddice*sice+sao(i,j,1)                      &
!                     *(schwell-ddice*rhoicwa))/schwell

                sao(i,j,1) =(ddice*sice+sao(i,j,1)*draftold)/draftnew

                draftnew=ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa
                svnew=sao(i,j,1)*draftnew+sice*sictho(i,j)

                WRITE(io_stdout,*)'ice cut-off Salt Content Change : ', svnew-svold

         
#ifdef PBGC
                DO l=1,nocetra
                   ocetra(i,j,1,l)=ocetra(i,j,1,l)*draftold/draftnew
                ENDDO
#endif
                WRITE(io_stdout,*)'eisneu,sneu: ',i,j,sictho(i,j),sao(i,j,1)
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    !
    DO j=2,je-1
       DO i=2,ie-1
          IF(weto(i,j,1).GT.0.5.AND.giph(2*i,2*j).GT.0.)THEN
             icou=icou+1
             !        wetie=wetie+dlxp(1,j)*dlyp(1,j)
             arges=arges+dlxp(i,j)*dlyp(i,j)
             swmin=MIN(swmin,fswr(i,j))
             eismax=MAX(eismax,sictho(i,j))
             zmi=MIN(zmi,zo(i,j))
             zma=MAX(zma,zo(i,j))
             volice=volice+dlxp(i,j)*dlyp(i,j)*sictho(i,j)
             volsno=volsno+dlxp(i,j)*dlyp(i,j)*sicsno(i,j)
             !-et
             dlxpdlyp=dlxp(i,j)*dlyp(i,j)
             IF(sicomo(i,j).LT.0.2)THEN
                dlxpdlyp=dlxpdlyp*MAX(0.,(sicomo(i,j)-0.1)/0.1)
             ENDIF
             areice=areice+dlxpdlyp
          ENDIF
       ENDDO
    ENDDO

    CALL global_sum(icou)
    CALL global_sum(arges,volice,volsno,areice)
    CALL global_max(eismax,zma)
    CALL global_min(swmin,zmi)

    WRITE(io_stdout,*)'eis: ',areice*1.e-12,volice*1.e-12          &
         ,volsno*1.e-12,eismax,' zeta: ',zmi,zma


  END SUBROUTINE calc_icecutoff

  SUBROUTINE calc_psi

    USE mo_commo1
    USE mo_parallel
    USE mo_param1
    
    INTEGER I,J,K

    REAL zhilf_g(ie_g,je_g)

    zhilf_g(:,:)=0.0
    psiuwe_g(:,:)=0.0

    psiuwe(:,:)=0.

    DO j=2,je
       DO k=1,ke
          DO i=1,ie
             psiuwe(i,j)=psiuwe(i,j)+uko(i,j,k)*dlyu(i,j)*dduo(i,j,k)
          ENDDO
       ENDDO
    ENDDO

  !  now gather and broadcast global psiuwe array, the following loop
  !  sums up globally

    !        gather global psiuwe array

    CALL gather_arr(psiuwe,psiuwe_g,p_io)
    CALL p_bcast(psiuwe_g,p_io)

    IF(p_pe==p_io) THEN
       DO j=2,je_g
          DO i=1,ie_g
             psiuwe_g(i,j)=psiuwe_g(i,j-1)+psiuwe_g(i,j)
          ENDDO
       ENDDO
       DO j=2,je_g-1
          DO i=2,ie_g-1
             zhilf_g(i,j)=0.25*(psiuwe_g(i,j)+psiuwe_g(i,j-1) &
                  +psiuwe_g(i-1,j)+psiuwe_g(i-1,j-1))
          ENDDO
       ENDDO
       DO j=1,je_g
          DO i=1,ie_g
             psiuwe_g(i,j)=zhilf_g(i,j)
          ENDDO
       ENDDO
    ENDIF


    ! for later local use:
!!    psiuwe(1:ie,1:je) = psiuwe_g(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)

    CALL scatter_arr(psiuwe_g,psiuwe,p_io)


  END SUBROUTINE calc_psi

  SUBROUTINE calc_potential_energy_release(is)

    USE mo_commo1, only : sao, tho, zo, preff, weto, g, ddpo, dt 
    USE mo_param1, only : ie,je,ke
    USE mo_mean, only : tmepo,tmcdo,tmceo
    
    INTEGER :: I,J,K,IS
    real :: swi,thick,sh(ie,je),th(ie,je), rh(ie,je)

    if (is.eq.1) then
       
       ! mixing diagnostics
       ! potential energy

       tmepo(:,:)=0.

       DO j=1,je
          DO k=1,ke
             DO i=1,ie
                sh(i,j)=sao(i,j,k)
                th(i,j)=tho(i,j,k)
             ENDDO
          
             CALL adisitj(th,sh,preff(k),j)
             CALL rho1j(th,sh,preff(k),rh,j)
          
             swi=0.
             if (k.eq.1) swi=1.

             DO i=1,ie
                thick = ddpo(i,j,k)+zo(i,j)*swi
                tmepo(i,j)=tmepo(i,j)-weto(i,j,k)*rh(i,j)*thick*g*thick
             ENDDO
          ENDDO
       ENDDO

    endif
    
    if (is.eq.2) then

       ! mixing diagnostics
       ! potential energy after convection

       tmcdo(:,:)=0.
       
       DO j=1,je
           DO k=1,ke
             DO i=1,ie
                sh(i,j)=sao(i,j,k)
                th(i,j)=tho(i,j,k)
             ENDDO
          
             CALL adisitj(th,sh,preff(k),j)
             CALL rho1j(th,sh,preff(k),rh,j)
          
             swi=0.
             if (k.eq.1) swi=1.

             DO i=1,ie
                thick = ddpo(i,j,k)+zo(i,j)*swi
                tmcdo(i,j)=tmcdo(i,j)-weto(i,j,k)*rh(i,j)*thick*g*thick
             ENDDO
          ENDDO
       ENDDO

       tmcdo(:,:)=(tmepo(:,:)-tmcdo(:,:))/dt

    endif

    if (is.eq.3) then

       ! mixing diagnostics
       ! potential energy after vertical mixing

       tmceo(:,:)=0. 

       DO j=1,je
           DO k=1,ke
             DO i=1,ie
                sh(i,j)=sao(i,j,k)
                th(i,j)=tho(i,j,k)
             ENDDO
          
             CALL adisitj(th,sh,preff(k),j)
             CALL rho1j(th,sh,preff(k),rh,j)
          
             swi=0.
             if (k.eq.1) swi=1.

             DO i=1,ie
                thick = ddpo(i,j,k)+zo(i,j)*swi
                tmceo(i,j)=tmceo(i,j)-weto(i,j,k)*rh(i,j)*thick*g*thick
             ENDDO
          ENDDO
       ENDDO
       tmceo(:,:)=(tmepo(:,:)-tmceo(:,:))/dt
    endif   


  END SUBROUTINE calc_potential_energy_release

  SUBROUTINE diagnosis

    USE mo_param1
    USE mo_mpi
    USE mo_parallel
    USE mo_commo1
    USE mo_commoau1
    USE mo_commoau2
    USE mo_units
!---wiso-code
    USE mo_contra,          ONLY: ocectra
    USE mo_param1_add,      ONLY: ih2o16, ih2o18, ihDo16
!---wiso-code-end

#if defined (PBGC) && defined (__cpl_co2)
    USE mo_carbch, ONLY: co2flux
#endif

    INTEGER i,j,k,m,n
    INTEGER jj,jjj,ii,iii

    REAL rkmaxlon, rkmaxlat,grarad                              &
               ,rgmaxlat,rgmaxlon,rsmaxlat,rsmaxlon,faks          &
              ,rberi,RINM,RUTM,SMEDI,RIND,RUTD,RINF,RUTF,RINA &
              ,RUTA,SDAVIS,SDAV1,SDAV2,CPRHO



    CALL calc_icecutoff

!    CALL calc_psi

    

#ifdef DIAG
    !******************************************************************

    IF(p_pe == p_io) THEN
       pi=4.*ATAN(1.)
       grarad=180./pi

       ! gulfstreamindex
       psigolf=0.
       rgmaxlat=0.
       rgmaxLON=0.

       !!cdir novector
       DO j=1,je_g 
          DO i=1,ie_g
             IF (ibek_g(i,j).EQ.4.OR.ibek_g(i,j).EQ.5) THEN
                psigolf=MAX(psigolf,psiuwe_g(i,j))
                IF (psigolf .EQ. psiuwe_g(i,j)) THEN
                   rgmaxlat=grarad*giph_g((2*i)+1,(2*j)+1)
                   rgmaxlon=grarad*gila_g((2*i)+1,(2*j)+1)
                ENDIF
             ENDIF
          ENDDO
       ENDDO

       ! atl spgyre index
       psiaspg=0.
       rsmaxlat=0.
       rsmaxlon=0.
       !!cdir novector
       DO j=1,je_g
          DO i=1,ie_g
             IF (ibek_g(i,j).EQ.4) THEN
                psiaspg=MIN(psiaspg,psiuwe_g(i,j))
                IF (psiaspg .EQ. psiuwe_g(i,j)) THEN
                   rsmaxlat=grarad*giph_g((2*i)+1,(2*j)+1)
                   rsmaxlon=grarad*gila_g((2*i)+1,(2*j)+1)
                ENDIF
             ENDIF
          ENDDO
       ENDDO

       ! kuroshioindex
       psikuro=0.
       rkmaxlat=0.
       rkmaxlon=0.  
       !!cdir novector
       DO j=1,je_g-1
          DO i=1,ie_g-1
             IF (ibek_g(i,j).EQ.7) THEN
                psikuro=MAX(psikuro,psiuwe_g(i,j))
                IF (psikuro .EQ. psiuwe_g(i,j)) THEN
                   rkmaxlat=grarad*giph_g((2*i)+1,(2*j)+1)
                   rkmaxlon=grarad*gila_g((2*i)+1,(2*j)+1)
                ENDIF
             ENDIF
          ENDDO
       ENDDO

       WRITE(io_stdout,*)'golfstream: ', NINT(psigolf*1.e-6),                  &
            'at lat/lon ',rgmaxlat,rgmaxlon
       WRITE(io_stdout,*)'a-subpolar gyer: ', NINT(psiaspg*1.e-6),             &
            'at lat/lon ',rsmaxlat,rsmaxlon
       WRITE(io_stdout,*)'kuroshio: ', NINT(psikuro*1.e-6),                    &
            'at lat/lon ',rkmaxlat,rkmaxlon
       WRITE(io_stdout,*)'banda: '                                             &
            ,NINT((psiuwe_g(iban1,jban1)-psiuwe_g(iban2,jban2))*1.e-6)         &
            ,'  drake: '                                                       &
            ,NINT((psiuwe_g(idra1,jdra1)-psiuwe_g(idra2,jdra2))*1.e-6)         &
            ,'  bering*10: '                                                   &
            ,NINT((psiuwe_g(iber1,jber1)-psiuwe_g(iber2,jber2))*1.e-5)

    ENDIF
    !****************************************************************************
#endif

    faks=1./float(30*(lmont2+1-lmont1))

    DO j=1,je
       DO i=1,ie
          hflm(i,j)=hflm(i,j)+(qswo(i,j)+qlwo(i,j)+qlao(i,j)+qseo(i,j))*faks
          pmem(i,j)=pmem(i,j)+(prech(i,j)+eminpo(i,j))*faks
       ENDDO
    ENDDO

    DO j=1,je1 
       DO i=1,ie1
          eistrx(i,j)=eistrx(i,j)                                    &
               +sicuo(i,j)*0.5*(sictho(i,j)+sictho(i+1,j)               &
               +(sicsno(i,j)+sicsno(i+1,j))*rhosnic)*faks
          eistry(i,j)=eistry(i,j)                                    &
               +sicve(i,j)*0.5*(sictho(i,j)+sictho(i,j+1)               &
               +(sicsno(i,j)+sicsno(i,j+1))*rhosnic)*faks
       ENDDO
    ENDDO

    CALL bounds_exch('u',eistrx,'mo_diagnosis 1')
    CALL bounds_exch('v',eistry,'mo_diagnosis 2')

!====================================================================

    !
#ifdef DIAG
    !        diagnostic i
    !          transporte
    WRITE(io_stdout,*)'vor 714'
    !
    !cdir novector
    DO k=ke,1,-1
       tmedi(k)=0.0
       tberi(k)=0.
       tdavis(k)=0.
       tfaroe(k)=0.
       tdenmar(k)=0.
       IF(k.LT.ke)tmedi(k)=tmedi(k+1)
       IF(k.LT.ke)tberi(k)=tberi(k+1)
       IF(k.LT.ke)tfaroe(k)=tfaroe(k+1)
       IF(k.LT.ke)tdenmar(k)=tdenmar(k+1)
       IF(k.LT.ke)tdavis(k)=tdavis(k+1)
       !
       !hh          bering 
       IF (grid_family.EQ.2) THEN 
          i=idberi-p_ioff
          j=jdberi-p_joff
          rberi=0.
          IF(i>1 .AND. i<ie .AND. j>1 .AND. j<je) THEN
             rberi = uko(i,j,k)*dlyu(i,j)*dduo(i,j,k)
          ENDIF

          !hh          mediterranean
          i=idmedi-p_ioff
          j=jdmedi-p_joff
          rinm=0.
          rutm=0.
          smedi=0.
          IF(i>1 .AND. i<ie .AND. j>1 .AND. j<je) THEN
             rinm=rinm+MAX(0.,vke(i,j,k))*dlxv(i,j)*ddue(i,j,k)
             rutm=rutm+MIN(0.,vke(i,j,k))*dlxv(i,j)*ddue(i,j,k)
             smedi=sao(i,j,k)*weto(i,j,k)
          ENDIF

          !hh          denmark
          i=iddenm-p_ioff
          rind=0.
          rutd=0.
          IF(i>1 .AND. i<ie) THEN
             DO j=MAX(2,jddene-p_joff),MIN(je-1,jddenw-p_joff)
                rind=rind+MAX(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
                rutd=rutd+MIN(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
             ENDDO
          ENDIF

          !hh          faroe k+=> 700m
          !hh          if (k.ge.kdfaro) then
          j=jdfaro-p_joff
          rinf=0.
          rutf=0.
          IF(j>1 .AND. j<je) THEN
             DO i=MAX(2,idfarw-p_ioff),MIN(ie-1,idfare-p_ioff)
                rinf=rinf+MIN(0.,vke(i,j,k))*dlxv(i,j)*ddue(i,j,k)
                rutf=rutf+MAX(0.,vke(i,j,k))*dlxv(i,j)*ddue(i,j,k)
             ENDDO
          ENDIF

          !hh          davis
          rina=0.
          ruta=0.   
          i=iddavi-p_ioff
          IF(i>1 .AND. i<ie) THEN
             DO  j=MAX(2,jddavw-p_joff),MIN(je-1,jddave-p_joff)
                rina=rina+MAX(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
                ruta=ruta+MIN(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
             ENDDO
          ENDIF
       ELSE
          i=idberi-p_ioff
          j=jdberi-p_joff
          rberi=0.
          IF(i>1 .AND. i<ie .AND. j>1 .AND. j<je) THEN
             rberi=vke(i,j,k)*dlxv(i,j)*ddue(i,j,k)
          ENDIF
          !
          !hh          mediterranean
          i=idmedi-p_ioff
          j=jdmedi-p_joff
          rinm=0.
          rutm=0.
          smedi=0.
          IF(i>1 .AND. i<ie .AND. j>1 .AND. j<je) THEN
             rinm=rinm+MAX(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
             rutm=rutm+MIN(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
             smedi=sao(i,j,k)*weto(i,j,k)
          ENDIF
          ! 
          !hh          denmark
          j=jddenm-p_joff
          rind=0.
          rutd=0.
          IF(j>1 .AND. j<je) THEN
             DO i=MAX(2,iddene-p_ioff),MIN(ie-1,iddenw-p_ioff)
                rind=rind+MAX(0.,vke(i,j,k))*dlxv(i,j)*ddue(i,j,k)
                rutd=rutd+MIN(0.,vke(i,j,k))*dlxv(i,j)*ddue(i,j,k)
             ENDDO
          ENDIF
          i=iddenm-p_ioff
          IF(i>1 .AND. i<ie) THEN
             DO j=MAX(2,jddenn-p_joff),MIN(je-1,jddens-p_joff)
                rind=rind+MIN(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
                rutd=rutd+MAX(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
             ENDDO
          ENDIF
          !
          !hh          faroe k+=> 700m
          !hh          if (k.ge.kdfaro) then
          i=idfaro-p_ioff
          rinf=0.
          rutf=0.
          IF(i>1 .AND. i<ie) THEN
             DO j=MAX(2,jdfarn-p_joff),MIN(je-1,jdfars-p_joff)
                rinf=rinf+MIN(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
                rutf=rutf+MAX(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
             ENDDO
          ENDIF


          !hh          davis
          rina=0.
          ruta=0.   
          i=iddavi-p_ioff
          IF(i>1 .AND. i<ie) THEN
             DO  j=MAX(2,jddavw-p_joff),MIN(je-1,jddave-p_joff)
                rina=rina+MAX(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
                ruta=ruta+MIN(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
             ENDDO
          ENDIF
       ENDIF

       CALL global_sum(rberi,smedi,rinm,rutm,rind,rutd,rinf,rutf,rina,ruta)

       tberi(k)=rberi+tberi(k)
       tdavis(k)=rina+ruta+tdavis(k)
       tmedi(k)=rinm+rutm+tmedi(k)
       tfaroe(k)=rinf+rutf+tfaroe(k)
       tdenmar(k)=rind+rutd+tdenmar(k)

       !jj test
       IF (MOD(ldays,30).EQ.5) THEN
          WRITE(io_stdout,'(f5.0,a,2f6.3,a,2f6.3,a,f6.3,a,f6.3,a,2f6.3,f6.2)') &
               tiestu(k),                           &
               ' denm:  ', (rind+rutd)*1.e-6, tdenmar(k)*1.e-6, &
               ' faroe: ', (rinf+rutf)*1.e-6, tfaroe(k) *1.e-6, &
               ' ber:   ',                    tberi(k)  *1.e-6, &
               ' dav:   ',                    tdavis(k) *1.e-6, &
               ' med:   ', (rinm+rutm)*1.e-6, tmedi(k)  *1.e-6, smedi
       ENDIF

    ENDDO

    !          ice transport fram strait
    sfram=0.

    i=idfram-p_ioff
    IF(i>1 .AND. i<ie) THEN
       !cdir novector
       DO j=MAX(2,jdfraw-p_joff),MIN(je-1,jdfrae-p_joff)
          !hh             sfram=sfram+0.5*(sicuo(i,j)+sicuo(i-1,j))
          !hh     x            *(sictho(i,j)+sicsno(i,j)*rhosno/rhoice)*dlyp(i,j)
          sfram=sfram+((ABS(sicuo(i,j))+sicuo(i,j))                    &
               *(sictho(i,j)+sicsno(i,j)*rhosnic)                 &
               +((sicuo(i,j)-ABS(sicuo(i,j))))               &
               *(sictho(i+1,j)+sicsno(i+1,j)*rhosnic))            &
               *0.5*dlyu(i,j)
          !
       ENDDO
    ENDIF
    CALL global_sum(sfram)
    sdavis=0.
    i=iddavi-p_ioff
    IF(i>1 .AND. i<ie) THEN
       !cdir novector
       DO j=MAX(2,jddavw-p_joff),MIN(je-1,jddave-p_joff)
          !hh          sdavis=sdavis+0.5*(sicuo(i,j)+sicuo(i-1,j))
          !hh  x             *(sictho(i,j)+sicsno(i,j)*rhosnic)*dlyp(i,j)
          sdavis=sdavis+((ABS(sicuo(i,j))+sicuo(i,j))               &
               *(sictho(i,j)+sicsno(i,j)*rhosnic)                 &
               +((sicuo(i,j)-ABS(sicuo(i,j))))               &
               *(sictho(i+1,j)+sicsno(i+1,j)*rhosnic))            &
               *0.5*dlyu(i,j)
          !
       ENDDO
    ENDIF
    CALL global_sum(sdavis)

    IF (grid_family.EQ.2) THEN
       i=iddav1-p_ioff
       j=jddav1-p_joff
       sdav1=0
       IF(i>1 .AND. i<ie .AND. j>1 .AND. j<je) THEN
          sdav1=((ABS(sicuo(i,j))+sicuo(i,j))                       &
               *(sictho(i,j)+sicsno(i,j)*rhosnic)                 &
               +((sicuo(i,j)-ABS(sicuo(i,j))))               &
               *(sictho(i+1,j)+sicsno(i+1,j)*rhosnic))            &
               *0.5*dlyu(i,j)
       ENDIF
       !
       i=iddav2-p_ioff
       j=jddav2-p_joff
       sdav2=0
       IF(i>1 .AND. i<ie .AND. j>1 .AND. j<je) THEN
          sdav2=((ABS(sicuo(i,j))+sicuo(i,j))                       &
               *(sictho(i,j)+sicsno(i,j)*rhosnic)                 &
               +((sicuo(i,j)-ABS(sicuo(i,j))))               &
               *(sictho(i+1,j)+sicsno(i+1,j)*rhosnic))            &
               *0.5*dlyu(i,j)
       ENDIF
       CALL global_sum(sdav1,sdav2)

       IF(MOD(ldays,10).EQ.5)                                       &
            WRITE(io_stdout,*)'eis fram ', NINT(sfram), ' davis ',NINT(sdavis)
    ENDIF

    DO n=1,nbox
       hflb(n)=0.
       wflb(n)=0.
       eiscb(n)=0.
       eisab(n)=0.
       DO k=1,ke
          arem(k,n)=0.
          tquer(k,n)=0.
!---wiso-code
          o16quer(k,n)=0.
          o18quer(k,n)=0.
          hdoquer(k,n)=0.
!---wiso-code-end       
          squer(k,n)=0.
          !
          avquer(k,n)=0.
          dvquer(k,n)=0.
       ENDDO
    ENDDO
    !hh test
    !cdir novector
    DO k=1,ke
       DO j=2,je-1
          DO i=2,ie-1
             n=ibek(i,j)
             IF(weto(i,j,k).GT.0.5.AND.n.GT.0)THEN
                arem(k,n)=arem(k,n)+dlxp(i,j)*dlyp(i,j)
                tquer(k,n)=tquer(k,n)+dlxp(i,j)*dlyp(i,j)*tho(i,j,k)
!---wiso-code
                o16quer(k,n)=o16quer(k,n)+dlxp(i,j)*dlyp(i,j)*ocectra(i,j,k,ih2o16)
                o18quer(k,n)=o18quer(k,n)+dlxp(i,j)*dlyp(i,j)*ocectra(i,j,k,ih2o18)
                hdoquer(k,n)=hdoquer(k,n)+dlxp(i,j)*dlyp(i,j)*ocectra(i,j,k,ihDo16)
!---wiso-code-end       
                squer(k,n)=squer(k,n)+dlxp(i,j)*dlyp(i,j)*sao(i,j,k)
                avquer(k,n)=avquer(k,n)+dlxp(i,j)*dlyp(i,j)*avo(i,j,k)
                dvquer(k,n)=dvquer(k,n)+dlxp(i,j)*dlyp(i,j)*dvo(i,j,k)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    WRITE(io_stdout,*)'check -5'
    !cdir novector
    DO j=2,je-1
       DO i=2,ie-1
          n=ibek(i,j)
          IF(weto(i,j,1).GT.0.5.AND.n.GT.0)THEN
             hflb(n)=hflb(n)+dlxp(i,j)*dlyp(i,j)*                   &
                  (qswo(i,j)+qlwo(i,j)+qlao(i,j)+qseo(i,j))
             wflb(n)=wflb(n)+dlxp(i,j)*dlyp(i,j)                    &
                  *(prech(i,j)+eminpo(i,j))
             eiscb(n)=eiscb(n)+dlxp(i,j)*dlyp(i,j)                  &
                  *(sictho(i,j)+(sicsno(i,j)*rhosnic))
             IF(sicomo(i,j).GE.0.15)                                &
                  eisab(n)=eisab(n)+dlxp(i,j)*dlyp(i,j)
          ENDIF
       ENDDO
    ENDDO
    CALL global_sum(hflb)
    CALL global_sum(wflb)
    CALL global_sum(eiscb)
    CALL global_sum(eisab)
    CALL global_sum(arem)
    arem(:,:) = arem(:,:) + 1.e-10
    CALL global_sum(tquer)
!---wiso-code
    CALL global_sum(o16quer)
    CALL global_sum(o18quer)
    CALL global_sum(hdoquer)
!---wiso-code-end       
    CALL global_sum(squer)
    CALL global_sum(avquer)
    CALL global_sum(dvquer)




    DO k=1,ke
       DO n=1,nmecen
          tvquer(k,n)=0.
          svquer(k,n)=0.
       ENDDO
    ENDDO

    DO k=ke,1,-1
       DO n=1,nmecen
          IF(k.EQ.ke)THEN
             tvquer(k,n)=0.
             svquer(k,n)=0.
          ELSE
             tvquer(k,n)=tvquer(k+1,n)
             svquer(k,n)=svquer(k+1,n)
          ENDIF

          cprho=rocp
          IF (grid_family.EQ.2) THEN

             i=imerci(n)-p_ioff
             IF(i>1 .AND. i<ie) THEN
                !cdir novector
                DO j=MAX(2,jlink(n)-p_joff),MIN(je-1,jrech(n)-p_joff)
                   tvquer(k,n)=tvquer(k,n)+uko(i,j,k)*dlyu(i,j)*dduo(i,j,k)&
                        *(tho(i,j,k)+tho(i+1,j,k))*0.5*cprho
                   svquer(k,n)=svquer(k,n)+uko(i,j,k)*dlyu(i,j)*dduo(i,j,k)&
                        *(sao(i,j,k)+sao(i+1,j,k))*0.5
                ENDDO
             ENDIF

          ELSE 
             j=jmerci(n)-p_joff
             IF(j>1 .AND. j<je) THEN
                !cdir novector
                DO i=MAX(2,ilink(n)-p_ioff),MIN(ie-1,irech(n)-p_ioff)
                   tvquer(k,n)=tvquer(k,n)+vke(i,j,k)*dlxv(i,j)*ddue(i,j,k)&
                        *(tho(i,j,k)+tho(i,j+1,k))*0.5*cprho
                   svquer(k,n)=svquer(k,n)+vke(i,j,k)*dlxv(i,j)*ddue(i,j,k)&
                        *(sao(i,j,k)+sao(i,j+1,k))*0.5
                ENDDO
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    CALL global_sum(tvquer)
    CALL global_sum(svquer)

    WRITE(io_stdout,*)'atlantic overturning: '
    DO k=ke,1,-1
       DO n=1,nmecen
          IF(k.EQ.ke)THEN
             tmerci(k,n)=0.
          ELSE
             tmerci(k,n)=tmerci(k+1,n)
          ENDIF
          IF (grid_family.EQ.2) THEN
             i=imerci(n)
             IF(weto_g(i,jlink(n),k)+weto_g(i,jrech(n),k).GT.0.5)THEN
                WRITE(io_stdout,*)'alarm! jlink,jrech ',n,i,jlink(n)   &
                     ,weto_g(i,jlink(n),k),jrech(n),weto_g(i,jrech(n),k)
             ENDIF
             i=imerci(n)-p_ioff
             IF(i>1 .AND. i<ie) THEN
                !cdir novector
                DO j=MAX(2,jlink(n)-p_joff),MIN(je-1,jrech(n)-p_joff)
                   tmerci(k,n)=tmerci(k,n)+uko(i,j,k)*dlyu(i,j)*dduo(i,j,k)
                ENDDO
             ENDIF
          ELSE
             j=jmerci(n)
             IF(weto_g(ilink(n),j,k)+weto_g(irech(n),j,k).GT.0.5)THEN
                WRITE(io_stdout,*)'alarm! ilink,irech ',n,j,ilink(n)   &
                     ,irech(n),weto_g(ilink(n),j,k),weto_g(irech(n),j,k)
             ENDIF
             j=jmerci(n)-p_joff
             IF(j>1 .AND. j<je) THEN
                !cdir novector
                DO i=MAX(2,ilink(n)-p_ioff),MIN(ie-1,irech(n)-p_ioff)
                   tmerci(k,n)=tmerci(k,n)+vke(i,j,k)*dlxv(i,j)*ddue(i,j,k)
                ENDDO
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    CALL global_sum(tmerci)

    DO k=1,ke
       WRITE(io_stdout,'(f8.0,10f10.3)')tiestw(k),(tmerci(k,n)*1.e-6,n=1,nmecen)
    ENDDO

    IF (grid_family.EQ.2) THEN
       WRITE(io_stdout,'(f8.0,10f10.3)')tiestw(1)                   &      
            ,(float(imerci(n)),n=1,nmecen)
    ELSE  
       WRITE(io_stdout,'(f8.0,10f10.3)')tiestw(1)                   &
            ,(float(jmerci(n)),n=1,nmecen)
    ENDIF

    !hh        find aabw(max) and nadw(min) at 
    !hh        aabw below level kaabw
    !hh        nadw above level knadw 
    DO n=1,nmecen
       aabw(n)=0.
       rnadw(n)=0.
    ENDDO

    DO n=2,nmecen
       IF (grid_family.EQ.2) THEN
          !cdir novector
          DO k=ke,kaabw,-1
             aabw(n)=MIN(tmerci(k,n),aabw(n))
          ENDDO
       ELSE  
          !cdir novector
          DO k=ke,kaabw,-1
             aabw(n)=MAX(tmerci(k,n),aabw(n))
          ENDDO
       ENDIF

       IF (grid_family.EQ.2) THEN
          !cdir novector
          DO k=ke,knadw,-1
             rnadw(n)=MAX(tmerci(k,n),rnadw(n))
          ENDDO
       ELSE
          !cdir novector 
          DO k=ke,knadw,-1
             rnadw(n)=MIN(tmerci(k,n),rnadw(n))
          ENDDO
       ENDIF
    ENDDO
    !hh           write(io_stdout,*)'check 4'
    !cdir novector
    DO n=1,nbox-1
       eisab(nbox)=eisab(nbox)+eisab(n)
       eiscb(nbox)=eiscb(nbox)+eiscb(n)
       hflb(nbox)=hflb(nbox)+hflb(n)
       wflb(nbox)=wflb(nbox)+wflb(n)
       !cdir novector
       DO k=1,ke
          arem(k,nbox)=arem(k,nbox)+arem(k,n)
          squer(k,nbox)=squer(k,nbox)+squer(k,n)
          tquer(k,nbox)=tquer(k,nbox)+tquer(k,n)
!---wiso-code
          o16quer(k,nbox)=o16quer(k,nbox)+o16quer(k,n)
          o18quer(k,nbox)=o18quer(k,nbox)+o18quer(k,n)
          hdoquer(k,nbox)=hdoquer(k,nbox)+hdoquer(k,n)
!---wiso-code-end       
       ENDDO
    ENDDO
    !hh           write(io_stdout,*)'check 5'
    !
    !hh        default 20 levels
    ldenmar=9
    lfaroe=11
    ltq1=1
    ltq2=7
    ltq3=11
    ltq4=16

    IF (ke.EQ.40) THEN
       ldenmar=16
       lfaroe=19
       ltq1=1
       ltq2=13
       ltq3=21
       ltq4=31
    ENDIF


    IF (ke.EQ.23) THEN
       ldenmar=11
       lfaroe=12
       ltq1=1
       ltq2=9
       ltq3=14
       ltq4=19
    ENDIF

    IF (ke.EQ.30) THEN

       ldenmar=12
       lfaroe=16
       ltq1=1
       ltq2=9
       ltq3=16
       ltq4=22
    ENDIF

#if defined (PBGC) && defined (__cpl_co2)

     !sum up the global ocean mean co2flux GtC/yr 

     cflux=0.

     ii=2
     iii=ie-1
     jj=2
     jjj=je-1

!     IF (have_g_is) ii=2
!     IF (have_g_ie) iii=ie-1     
!     IF (have_g_js) jj=3
!     IF (have_g_je) jjj=je-1

     DO i=ii,iii
        DO j=jj,jjj
           cflux=cflux+co2flux(i,j)*dlxp(i,j)*dlyp(i,j)        &
                        *(12./44.011)*86400.*365.*weto(i,j,1)
        ENDDO
     ENDDO

     CALL global_sum(cflux)
#endif


    IF (p_pe == p_io) THEN

!FR  ITSDIAG is integer and included in NAMELIST 'OCECTL':
    if ( itsdiag.ne.0 ) then
      call write_timeseries
    endif

    WRITE(io_ou_f125,'(10e13.6)')float(lyears)                   &
         +(float(lmonts-1)+float(ldays)/float(monlen(lmonts)))/12.    &
         ,psigolf, psikuro                                            &
         ,((psiuwe_g(iban1,jban1)-psiuwe_g(iban2,jban2))*1.e0)        &
         ,((psiuwe_g(idra1,jdra1)-psiuwe_g(idra2,jdra2))*1.e0)        &
         ,((psiuwe_g(iber1,jber1)-psiuwe_g(iber2,jber2))*1.e0)        &
         ,psiaspg                                                     &
#if defined (PBGC) && defined (__cpl_co2)
         ,atm_co2,cflux,(aabw(m),rnadw(m),m=2,nmecen)                 &
#else
         ,(aabw(m),rnadw(m),m=1,nmecen)                               &
#endif
         ,(tvquer(1,m),svquer(1,m),tmerci(1,m),m=1,nmecen)            &
         ,(tvquer(1,m)-tvquer(1,1),svquer(1,m)-svquer(1,1),m=2,nmecen)&
         ,tberi(1),tdenmar(ldenmar),tfaroe(lfaroe),sfram              &
         ,(eisab(n),eiscb(n),hflb(n),wflb(n)                          &
         ,tquer(ltq1,n)/arem(ltq1,n),squer(ltq1,n)/arem(ltq1,n)       &
         ,tquer(ltq2,n)/arem(ltq2,n),squer(ltq2,n)/arem(ltq2,n)       &
         ,tquer(ltq3,n)/arem(ltq3,n),squer(ltq3,n)/arem(ltq3,n)       &
         ,tquer(ltq4,n)/arem(ltq4,n),squer(ltq4,n)/arem(ltq4,n)       &
!---wiso-code
         ,(o18quer(ltq1,n)/o16quer(ltq1,n)/2005.2e-6 - 1.)*1000.,(hdoquer(ltq1,n)/o16quer(ltq1,n)/155.76e-6 - 1.)*1000.  &
         ,(o18quer(ltq2,n)/o16quer(ltq2,n)/2005.2e-6 - 1.)*1000.,(hdoquer(ltq2,n)/o16quer(ltq2,n)/155.76e-6 - 1.)*1000.  &
         ,(o18quer(ltq3,n)/o16quer(ltq3,n)/2005.2e-6 - 1.)*1000.,(hdoquer(ltq3,n)/o16quer(ltq3,n)/155.76e-6 - 1.)*1000.  &
         ,(o18quer(ltq4,n)/o16quer(ltq4,n)/2005.2e-6 - 1.)*1000.,(hdoquer(ltq4,n)/o16quer(ltq4,n)/155.76e-6 - 1.)*1000.  &
!---wiso-code-end       
         ,n=1,nbox)
    
    ENDIF

#endif /*diag*/


    !      **************************************************************
    !      *******************end of common diagnostic*******************
    !      **************************************************************
  END SUBROUTINE diagnosis


subroutine write_timeseries
! Subroutine to write timeseries data in the PSEUDO EXTRA Format.


USE mo_kind
USE mo_param1
USE mo_commo1
USE mo_units, only: io_ou_tsvar, io_ou_tsdesc, io_ou_tsunit, io_ou_f126


!FR  Some more variables:
INTEGER          :: i,m,n
INTEGER(KIND=i4) :: lentimser, noffset, idate, ione
INTEGER          :: nanf,nend,ndtday,lmonts1
REAL, DIMENSION(:), allocatable  :: timeser
REAL(KIND=sp), DIMENSION(:), allocatable  :: write_timeser

NDTDAY=NINT(86400./DT)


!FR  The length of time series is determined and then memory is allocated.
    lentimser = 0+2+3+1
#if defined (PBGC) && defined (__cpl_co2)
    lentimser = lentimser + 2 + 2*(nmecen-1)
#else
    lentimser = lentimser + 2*nmecen
#endif
    lentimser = lentimser + 3*nmecen + 2*(nmecen-1) + 4 + nbox*(4+4*2)
!---wiso-code
    lentimser = lentimser + nbox*(4*2)
!---wiso-code-end

    allocate(      timeser( lentimser ) )          
    allocate(write_timeser( lentimser ) )          
        
!FR   All output parameters are gathered in one vector.
    timeser(1) = psigolf
    timeser(2) = psikuro
    timeser(3) = ((psiuwe_g(iban1,jban1)-psiuwe_g(iban2,jban2))*1.e0)
    timeser(4) = ((psiuwe_g(idra1,jdra1)-psiuwe_g(idra2,jdra2))*1.e0) 
    timeser(5) = ((psiuwe_g(iber1,jber1)-psiuwe_g(iber2,jber2))*1.e0) 
    timeser(6) = psiaspg
!FR
! Write three files with information for the time series plots:
! 1. Variables as used here in the Fortran code.
! 2. Description as used in scientific language.
! 3. Units as recommended by the System International (SI).
    if ( ltswrite ) then
990 format(i3,1x,a9,i3,a1,i3,a11,i3,a1,i3,a1)
      write(io_ou_tsvar,'(i3,1x,a7)') 1,'psigolf'
      write(io_ou_tsvar,'(i3,1x,a7)') 2,'psikuro'
      write(io_ou_tsvar,990) 3,  &
           'psiuwe_g(',iban1,',',jban1,')-psiuwe_g(',iban2,',',jban2,')'
      write(io_ou_tsvar,990) 4,  &
           'psiuwe_g(',idra1,',',jdra1,')-psiuwe_g(',idra2,',',jdra2,')'
      write(io_ou_tsvar,990) 5,  &
           'psiuwe_g(',iber1,',',jber1,')-psiuwe_g(',iber2,',',jber2,')'
      write(io_ou_tsvar,'(i3,1x,a7)') 6,'psiaspg'
      write(io_ou_tsdesc,'(i3,1x,a27)') 1,'Stream function Gulf stream'
      write(io_ou_tsdesc,'(i3,1x,a24)') 2,'Stream function Kuroshio'
      write(io_ou_tsdesc,'(i3,1x,a28)') 3,'Stream function Banda Strait'
      write(io_ou_tsdesc,'(i3,1x,a29)') 4,'Stream function Drake Passage'
      write(io_ou_tsdesc,'(i3,1x,a29)') 5,'Stream function Bering Strait'
      write(io_ou_tsdesc,'(i3,1x,a44)') 6,  &
                         'Stream function North Atlantic subpolar gyre'
      write(io_ou_tsunit,'(i3,1x,a5)') 1,'m^3/s'
      write(io_ou_tsunit,'(i3,1x,a5)') 2,'m^3/s'
      write(io_ou_tsunit,'(i3,1x,a5)') 3,'m^3/s'
      write(io_ou_tsunit,'(i3,1x,a5)') 4,'m^3/s'
      write(io_ou_tsunit,'(i3,1x,a5)') 5,'m^3/s'
      write(io_ou_tsunit,'(i3,1x,a5)') 6,'m^3/s'
    endif
    noffset = 6
#if defined (PBGC) && defined (__cpl_co2)
    timeser(noffset + 1) = atm_co2
    timeser(noffset + 2) = cflux
!
    if ( ltswrite ) then
      write(io_ou_tsvar,'(i3,1x,a7)') noffset + 1,'atm_co2'
      write(io_ou_tsvar,'(i3,1x,a5)') noffset + 2,'cflux'
      write(io_ou_tsdesc,'(i3,1x,a15)') noffset + 1,'Atmospheric CO2'
      write(io_ou_tsdesc,'(i3,1x,a11)') noffset + 2,'Carbon flux'
      write(io_ou_tsunit,'(i3,1x,a3)') noffset + 1,'ppm'
      write(io_ou_tsunit,'(i3,1x,a7)') noffset + 2,'mol/m^2'
    endif
    noffset = noffset + 2

    do m=2,nmecen
      timeser(noffset + (2*(m-1)-1) ) = aabw(m)
      timeser(noffset + (2*(m-1)  ) ) = rnadw(m)
    enddo
!
    if ( ltswrite ) then
991 format(i3,1x,a5,i1,a1)
992 format(i3,1x,a6,i1,a1)
      do m=2,nmecen
        write(io_ou_tsvar,991) noffset + (2*(m-1)-1),'aabw(',m,')'
        write(io_ou_tsvar,992) noffset + (2*(m-1)  ),'rnadw(',m,')'
        write(io_ou_tsdesc,'(i3,1x,a22)') noffset + (2*(m-1)-1),  &
                           'Antarctic bottom water'
        write(io_ou_tsdesc,'(i3,1x,a25)') noffset + (2*(m-1)  ),  &
                           'North Atlantic deep water'
        write(io_ou_tsunit,'(i3,1x,a5)' ) noffset + (2*(m-1)-1),'m^3/s'
        write(io_ou_tsunit,'(i3,1x,a5)' ) noffset + (2*(m-1)  ),'m^3/s'
      enddo
    endif
    noffset = noffset + (2*(nmecen-1))
#else
    do m=1,nmecen
      timeser(noffset + (2*(m  )-1) ) = aabw(m)
      timeser(noffset + (2*(m  )  ) ) = rnadw(m)
    enddo
!
    if ( ltswrite ) then
991 format(i3,1x,a5,i1,a1)
992 format(i3,1x,a6,i1,a1)
      do m=1,nmecen
        write(io_ou_tsvar,991) noffset + (2*(m  )-1),'aabw(',m,')'
        write(io_ou_tsvar,992) noffset + (2*(m  )  ),'rnadw(',m,')'
        write(io_ou_tsdesc,'(i3,1x,a22)') noffset + (2*(m  )-1),  &
                           'Antarctic bottom water'
        write(io_ou_tsdesc,'(i3,1x,a25)') noffset + (2*(m  )  ),  &
                           'North Atlantic deep water'
        write(io_ou_tsunit,'(i3,1x,a5)' ) noffset + (2*(m  )-1),'m^3/s'
        write(io_ou_tsunit,'(i3,1x,a5)' ) noffset + (2*(m  )  ),'m^3/s'
      enddo
    endif
    noffset = noffset + (2*nmecen)
#endif

    do m=1,nmecen
      timeser(noffset + (3*(m  )-2) ) = tvquer(1,m)
      timeser(noffset + (3*(m  )-1) ) = svquer(1,m)
      timeser(noffset + (3*(m  )  ) ) = tmerci(1,m)
    enddo
!
    if ( ltswrite ) then
993 format(i3,1x,a9,i1,a10,i1,a2,i3)
      do m=1,nmecen
        write(io_ou_tsvar,993) noffset + (3*(m  )-2),  &
             'tvquer(1,',m,'), jmerci(',m,')=',jmerci(m)
        write(io_ou_tsvar,993) noffset + (3*(m  )-1),  &
             'svquer(1,',m,'), jmerci(',m,')=',jmerci(m)
        write(io_ou_tsvar,993) noffset + (3*(m  )  ),  &
             'tmerci(1,',m,'), jmerci(',m,')=',jmerci(m)
        write(io_ou_tsdesc,'(i3,1x,a32)') noffset + (3*(m  )-2), &
                           'Meridional temperature transport'
        write(io_ou_tsdesc,'(i3,1x,a29)') noffset + (3*(m  )-1), &
                           'Meridional salinity transport'
        write(io_ou_tsdesc,'(i3,1x,a20)') noffset + (3*(m  )  ),  &
                           'Atlantic overturning'
        write(io_ou_tsunit,'(i3,1x,a8)')noffset + (3*(m  )-2),'W*m^3/kg'
        write(io_ou_tsunit,'(i3,1x,a9)')noffset + (3*(m  )-1),'m^3/s*psu'
        write(io_ou_tsunit,'(i3,1x,a5)')noffset + (3*(m  )  ),'m^3/s'
      enddo
    endif
    noffset = noffset + (3*(nmecen  ))

    do m=2,nmecen
      timeser(noffset + (2*(m-1)-1) ) = tvquer(1,m)-tvquer(1,1)
      timeser(noffset + (2*(m-1)  ) ) = svquer(1,m)-svquer(1,1)
    enddo
!
    if ( ltswrite ) then
994 format(i3,1x,a9,i1,a22,i1,a2,i3)
      do m=2,nmecen
        write(io_ou_tsvar,994) noffset + (2*(m-1)-1), &
             'tvquer(1,',m,')-tvquer(1,1), jmerci(',m,')=',jmerci(m)
        write(io_ou_tsvar,994) noffset + (2*(m-1)  ), &
             'svquer(1,',m,')-svquer(1,1), jmerci(',m,')=',jmerci(m)
        write(io_ou_tsdesc,'(i3,1x,a43)') noffset + (2*(m-1)-1),  &
                           'Meridional temperature transport difference'
        write(io_ou_tsdesc,'(i3,1x,a40)') noffset + (2*(m-1)  ),  &
                           'Meridional salinity transport difference'
        write(io_ou_tsunit,'(i3,1x,a8)')noffset + (2*(m-1)-1),'W*m^3/kg' 
        write(io_ou_tsunit,'(i3,1x,a9)')noffset + (2*(m-1)  ),'m^3/s*psu'
      enddo
    endif
    noffset = noffset + (2*(nmecen-1))

    timeser(noffset + 1) = tberi(1)
    timeser(noffset + 2) = tdenmar(ldenmar)
    timeser(noffset + 3) = tfaroe(lfaroe)
    timeser(noffset + 4) = sfram 
!
    if ( ltswrite ) then
995 format(i3,1x,a8,i2,a1)
996 format(i3,1x,a7,i2,a1)
      write(io_ou_tsvar,'(i3,1x,a8)') noffset + 1,'tberi(1)'
      write(io_ou_tsvar,995) noffset + 2,'tdenmar(',ldenmar,')'
      write(io_ou_tsvar,996) noffset + 3,'tfaroe(',lfaroe,')'
      write(io_ou_tsvar,'(i3,1x,a5)') noffset + 4,'sfram'
      write(io_ou_tsdesc,'(i3,1x,a23)') noffset + 1,  &
                                       'Transport Bering Strait'
      write(io_ou_tsdesc,'(i3,1x,a24)') noffset + 2,  &
                                       'Transport Denmark Strait'
      write(io_ou_tsdesc,'(i3,1x,a20)') noffset + 3,  &
                                       'Transport Faroe Bank'
      write(io_ou_tsdesc,'(i3,1x,a25)') noffset + 4,  &
                                       'Ice transport Fram Strait'
      write(io_ou_tsunit,'(i3,1x,a5)') noffset + 1, 'm^3/s'
      write(io_ou_tsunit,'(i3,1x,a5)') noffset + 2, 'm^3/s'
      write(io_ou_tsunit,'(i3,1x,a5)') noffset + 3, 'm^3/s'
      write(io_ou_tsunit,'(i3,1x,a5)') noffset + 4, 'm^3/s'
    endif
    noffset = noffset + 4

!---wiso-code

!---wiso-code: comment out the following code block to change offset for including o16 and hdo tracer
!    do n=1,nbox
!      timeser(noffset + (12*(n  )-11) ) = eisab(n)
!      timeser(noffset + (12*(n  )-10) ) = eiscb(n)
!      timeser(noffset + (12*(n  )- 9) ) = hflb (n)
!      timeser(noffset + (12*(n  )- 8) ) = wflb (n)
!      timeser(noffset + (12*(n  )- 7) ) = tquer(ltq1,n)/arem(ltq1,n)
!      timeser(noffset + (12*(n  )- 6) ) = squer(ltq1,n)/arem(ltq1,n)
!      timeser(noffset + (12*(n  )- 5) ) = tquer(ltq2,n)/arem(ltq2,n)
!      timeser(noffset + (12*(n  )- 4) ) = squer(ltq2,n)/arem(ltq2,n)
!      timeser(noffset + (12*(n  )- 3) ) = tquer(ltq3,n)/arem(ltq3,n)
!      timeser(noffset + (12*(n  )- 2) ) = squer(ltq3,n)/arem(ltq3,n)
!      timeser(noffset + (12*(n  )- 1) ) = tquer(ltq4,n)/arem(ltq4,n)
!      timeser(noffset + (12*(n  )   ) ) = squer(ltq4,n)/arem(ltq4,n)
!    enddo
!!
!    if ( ltswrite ) then
!997 format(i3,1x,a6,i1,a1)
!998 format(i3,1x,a5,i1,a1)
!999 format(i3,1x,a6,i2,a1,i1,a7,i2,a1,i1,a1)
!      do n=1,nbox
!        write(io_ou_tsvar,997) noffset + (12*(n  )-11),'eisab(',n,')'
!        write(io_ou_tsvar,997) noffset + (12*(n  )-10),'eiscb(',n,')'
!        write(io_ou_tsvar,998) noffset + (12*(n  )- 9),'hflb(',n,')'
!        write(io_ou_tsvar,998) noffset + (12*(n  )- 8),'wflb(',n,')'
!        write(io_ou_tsvar,999) noffset + (12*(n  )- 7),  &
!                          'tquer(',ltq1,',',n,')/arem(',ltq1,',',n,')'
!        write(io_ou_tsvar,999) noffset + (12*(n  )- 6),  &
!                          'squer(',ltq1,',',n,')/arem(',ltq1,',',n,')'
!        write(io_ou_tsvar,999) noffset + (12*(n  )- 5),  &
!                          'tquer(',ltq2,',',n,')/arem(',ltq2,',',n,')'
!        write(io_ou_tsvar,999) noffset + (12*(n  )- 4),  &
!                          'squer(',ltq2,',',n,')/arem(',ltq2,',',n,')'
!        write(io_ou_tsvar,999) noffset + (12*(n  )- 3),  &
!                          'tquer(',ltq3,',',n,')/arem(',ltq3,',',n,')'
!        write(io_ou_tsvar,999) noffset + (12*(n  )- 2),  &
!                          'squer(',ltq3,',',n,')/arem(',ltq3,',',n,')'
!        write(io_ou_tsvar,999) noffset + (12*(n  )- 1),  &
!                          'tquer(',ltq4,',',n,')/arem(',ltq4,',',n,')'
!        write(io_ou_tsvar,999) noffset + (12*(n  )   ),  &
!                          'squer(',ltq4,',',n,')/arem(',ltq4,',',n,')'
!        write(io_ou_tsdesc,'(i3,1x, a8)') noffset + (12*(n  )-11),  &
!                           'Ice area'
!        write(io_ou_tsdesc,'(i3,1x,a10)') noffset + (12*(n  )-10),  &
!                           'Ice volume'
!        write(io_ou_tsdesc,'(i3,1x,a11)') noffset + (12*(n  )- 9),  &
!                           'Global heat'
!        write(io_ou_tsdesc,'(i3,1x,a18)') noffset + (12*(n  )- 8),  &
!                           'Global fresh water'
!        write(io_ou_tsdesc,'(i3,1x,a29)') noffset + (12*(n  )- 7),  &
!                           'Globally averaged temperature'
!        write(io_ou_tsdesc,'(i3,1x,a26)') noffset + (12*(n  )- 6),  &
!                           'Globally averaged salinity'
!        write(io_ou_tsdesc,'(i3,1x,a29)') noffset + (12*(n  )- 5),  &
!                           'Globally averaged temperature'
!        write(io_ou_tsdesc,'(i3,1x,a26)') noffset + (12*(n  )- 4),  &
!                           'Globally averaged salinity'
!        write(io_ou_tsdesc,'(i3,1x,a29)') noffset + (12*(n  )- 3),  &
!                           'Globally averaged temperature'
!        write(io_ou_tsdesc,'(i3,1x,a26)') noffset + (12*(n  )- 2),  &
!                           'Globally averaged salinity'
!        write(io_ou_tsdesc,'(i3,1x,a29)') noffset + (12*(n  )- 1),  &
!                           'Globally averaged temperature'
!        write(io_ou_tsdesc,'(i3,1x,a26)') noffset + (12*(n  )   ),  &
!                           'Globally averaged salinity'
!        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (12*(n  )-11),'m^2'
!        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (12*(n  )-10),'m^3'
!        write(io_ou_tsunit,'(i3,1x, a1)') noffset + (12*(n  )- 9),'W'
!        write(io_ou_tsunit,'(i3,1x, a5)') noffset + (12*(n  )- 8),'m^3/s'
!        write(io_ou_tsunit,'(i3,1x, a8)') noffset + (12*(n  )- 7),'degree C'
!        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (12*(n  )- 6),'psu'
!        write(io_ou_tsunit,'(i3,1x, a8)') noffset + (12*(n  )- 5),'degree C'
!        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (12*(n  )- 4),'psu'
!        write(io_ou_tsunit,'(i3,1x, a8)') noffset + (12*(n  )- 3),'degree C'
!        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (12*(n  )- 2),'psu'
!        write(io_ou_tsunit,'(i3,1x, a8)') noffset + (12*(n  )- 1),'degree C'
!        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (12*(n  )   ),'psu'
!      enddo
!    endif
!    noffset = noffset + (12*(nbox  ))

    do n=1,nbox
      timeser(noffset + (20*(n  )-19) ) = eisab(n)
      timeser(noffset + (20*(n  )-18) ) = eiscb(n)
      timeser(noffset + (20*(n  )-17) ) = hflb (n)
      timeser(noffset + (20*(n  )-16) ) = wflb (n)
      timeser(noffset + (20*(n  )-15) ) = tquer(ltq1,n)/arem(ltq1,n)
      timeser(noffset + (20*(n  )-14) ) = squer(ltq1,n)/arem(ltq1,n)
      timeser(noffset + (20*(n  )-13) ) = tquer(ltq2,n)/arem(ltq2,n)
      timeser(noffset + (20*(n  )-12) ) = squer(ltq2,n)/arem(ltq2,n)
      timeser(noffset + (20*(n  )-11) ) = tquer(ltq3,n)/arem(ltq3,n)
      timeser(noffset + (20*(n  )-10) ) = squer(ltq3,n)/arem(ltq3,n)
      timeser(noffset + (20*(n  )- 9) ) = tquer(ltq4,n)/arem(ltq4,n)
      timeser(noffset + (20*(n  )- 8) ) = squer(ltq4,n)/arem(ltq4,n)

      timeser(noffset + (20*(n  )- 7) ) = (o18quer(ltq1,n)/o16quer(ltq1,n)/2005.2e-6 - 1.)*1000.
      timeser(noffset + (20*(n  )- 6) ) = (hdoquer(ltq1,n)/o16quer(ltq1,n)/155.76e-6 - 1.)*1000.
      timeser(noffset + (20*(n  )- 5) ) = (o18quer(ltq2,n)/o16quer(ltq2,n)/2005.2e-6 - 1.)*1000.
      timeser(noffset + (20*(n  )- 4) ) = (hdoquer(ltq2,n)/o16quer(ltq2,n)/155.76e-6 - 1.)*1000.
      timeser(noffset + (20*(n  )- 3) ) = (o18quer(ltq3,n)/o16quer(ltq3,n)/2005.2e-6 - 1.)*1000.
      timeser(noffset + (20*(n  )- 2) ) = (hdoquer(ltq3,n)/o16quer(ltq3,n)/155.76e-6 - 1.)*1000.
      timeser(noffset + (20*(n  )- 1) ) = (o18quer(ltq4,n)/o16quer(ltq4,n)/2005.2e-6 - 1.)*1000.
      timeser(noffset + (20*(n  )   ) ) = (hdoquer(ltq4,n)/o16quer(ltq4,n)/155.76e-6 - 1.)*1000.
    enddo
!
    if ( ltswrite ) then
997 format(i3,1x,a6,i1,a1)
998 format(i3,1x,a5,i1,a1)
999 format(i3,1x,a6,i2,a1,i1,a7,i2,a1,i1,a1)
      do n=1,nbox
        write(io_ou_tsvar,997) noffset + (20*(n  )-19),'eisab(',n,')'
        write(io_ou_tsvar,997) noffset + (20*(n  )-18),'eiscb(',n,')'
        write(io_ou_tsvar,998) noffset + (20*(n  )-17),'hflb(',n,')'
        write(io_ou_tsvar,998) noffset + (20*(n  )-16),'wflb(',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )-15),  &
                          'tquer(',ltq1,',',n,')/arem(',ltq1,',',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )-14),  &
                          'squer(',ltq1,',',n,')/arem(',ltq1,',',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )-13),  &
                          'tquer(',ltq2,',',n,')/arem(',ltq2,',',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )-12),  &
                          'squer(',ltq2,',',n,')/arem(',ltq2,',',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )-11),  &
                          'tquer(',ltq3,',',n,')/arem(',ltq3,',',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )-10),  &
                          'squer(',ltq3,',',n,')/arem(',ltq3,',',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )- 9),  &
                          'tquer(',ltq4,',',n,')/arem(',ltq4,',',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )- 8),  &
                          'squer(',ltq4,',',n,')/arem(',ltq4,',',n,')'

        write(io_ou_tsvar,999) noffset + (20*(n  )- 7),  &
                          'delta O18(',ltq1,',',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )- 6),  &
                          'delta HDO(',ltq1,',',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )- 5),  &
                          'delta O18(',ltq2,',',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )- 4),  &
                          'delta HDO(',ltq2,',',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )- 3),  &
                          'delta O18(',ltq3,',',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )- 2),  &
                          'delta HDO(',ltq3,',',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )- 1),  &
                          'delta O18(',ltq4,',',n,')'
        write(io_ou_tsvar,999) noffset + (20*(n  )   ),  &
                          'delta HDO(',ltq4,',',n,')'

        write(io_ou_tsdesc,'(i3,1x, a8)') noffset + (20*(n  )-19),  &
                           'Ice area'
        write(io_ou_tsdesc,'(i3,1x,a10)') noffset + (20*(n  )-18),  &
                           'Ice volume'
        write(io_ou_tsdesc,'(i3,1x,a11)') noffset + (20*(n  )-17),  &
                           'Global heat'
        write(io_ou_tsdesc,'(i3,1x,a18)') noffset + (20*(n  )-16),  &
                           'Global fresh water'
        write(io_ou_tsdesc,'(i3,1x,a29)') noffset + (20*(n  )-15),  &
                           'Globally averaged temperature'
        write(io_ou_tsdesc,'(i3,1x,a26)') noffset + (20*(n  )-14),  &
                           'Globally averaged salinity'
        write(io_ou_tsdesc,'(i3,1x,a29)') noffset + (20*(n  )-13),  &
                           'Globally averaged temperature'
        write(io_ou_tsdesc,'(i3,1x,a26)') noffset + (20*(n  )-12),  &
                           'Globally averaged salinity'
        write(io_ou_tsdesc,'(i3,1x,a29)') noffset + (20*(n  )-11),  &
                           'Globally averaged temperature'
        write(io_ou_tsdesc,'(i3,1x,a26)') noffset + (20*(n  )-10),  &
                           'Globally averaged salinity'
        write(io_ou_tsdesc,'(i3,1x,a29)') noffset + (20*(n  )- 9),  &
                           'Globally averaged temperature'
        write(io_ou_tsdesc,'(i3,1x,a26)') noffset + (20*(n  )- 8),  &
                           'Globally averaged salinity'

        write(io_ou_tsdesc,'(i3,1x,a27)') noffset + (20*(n  )- 7),  &
                           'Globally averaged delta O18'
        write(io_ou_tsdesc,'(i3,1x,a27)') noffset + (20*(n  )- 6),  &
                           'Globally averaged delta HDO'
        write(io_ou_tsdesc,'(i3,1x,a27)') noffset + (20*(n  )- 5),  &
                           'Globally averaged delta O18'
        write(io_ou_tsdesc,'(i3,1x,a27)') noffset + (20*(n  )- 4),  &
                           'Globally averaged delta HDO'
        write(io_ou_tsdesc,'(i3,1x,a27)') noffset + (20*(n  )- 3),  &
                           'Globally averaged delta O18'
        write(io_ou_tsdesc,'(i3,1x,a27)') noffset + (20*(n  )- 2),  &
                           'Globally averaged delta HDO'
        write(io_ou_tsdesc,'(i3,1x,a27)') noffset + (20*(n  )- 1),  &
                           'Globally averaged delta O18'
        write(io_ou_tsdesc,'(i3,1x,a27)') noffset + (20*(n  )   ),  &
                           'Globally averaged delta HDO'

        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (20*(n  )-19),'m^2'
        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (20*(n  )-18),'m^3'
        write(io_ou_tsunit,'(i3,1x, a1)') noffset + (20*(n  )-17),'W'
        write(io_ou_tsunit,'(i3,1x, a5)') noffset + (20*(n  )-16),'m^3/s'
        write(io_ou_tsunit,'(i3,1x, a8)') noffset + (20*(n  )-15),'degree C'
        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (20*(n  )-14),'psu'
        write(io_ou_tsunit,'(i3,1x, a8)') noffset + (20*(n  )-13),'degree C'
        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (20*(n  )-12),'psu'
        write(io_ou_tsunit,'(i3,1x, a8)') noffset + (20*(n  )-11),'degree C'
        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (20*(n  )-10),'psu'
        write(io_ou_tsunit,'(i3,1x, a8)') noffset + (20*(n  )- 9),'degree C'
        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (20*(n  )- 8),'psu'

        write(io_ou_tsunit,'(i3,1x, a8)') noffset + (20*(n  )- 7),'permill'
        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (20*(n  )- 6),'permill'
        write(io_ou_tsunit,'(i3,1x, a8)') noffset + (20*(n  )- 5),'permill'
        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (20*(n  )- 4),'permill'
        write(io_ou_tsunit,'(i3,1x, a8)') noffset + (20*(n  )- 3),'permill'
        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (20*(n  )- 2),'permill'
        write(io_ou_tsunit,'(i3,1x, a8)') noffset + (20*(n  )- 1),'permill'
        write(io_ou_tsunit,'(i3,1x, a3)') noffset + (20*(n  )   ),'permill'
      enddo
    endif
    noffset = noffset + (20*(nbox  ))

!---wiso-code-end

    ltswrite=.false.

!FR    Control: 'noffset' should equal 'lentimser':
      IF ( noffset .NE. lentimser ) THEN
         write(0,*) 'noffset: ',noffset,' lentimser: ',lentimser
      ELSE
!FR  Namelist parameter ITSDIAG switch to compute timeseries data:
! itsdiag=0:  No Output
! itsdiag=1:  one snapshot per day
! itsdiag=2:  monthly averaged snapshots
! itsdiag=3:  yearly  averaged snapshots
! itsdiag=4:  output every timestep
! itsdiag=5:  daily   average
! itsdiag=6:  monthly average of daily means
! itsdiag=7:  yearly  average of daily means 

! snapshots: one per day or each timestep depending on calling frequency
        IF (itsdiag.EQ.1.OR.itsdiag.EQ.4) THEN
           nanf=1
           nend=1

! Monthly averaged snapshots
        ELSEIF (itsdiag.EQ.2) THEN
          nanf=ldays
          nend=monlen(lmonts)

! Yearly  averaged snapshots
        ELSEIF (itsdiag.EQ.3) THEN
          IF (lmonts.EQ.1) THEN
            nanf=ldays
          ELSE
            nanf=0
            lmonts1=lmonts-1
            DO I=1,lmonts1
              nanf=nanf+monlen(i)
            ENDDO
            nanf=nanf+ldays
          ENDIF

          nend=0
          DO I=1,12
             nend=nend+monlen(i)
          ENDDO

! Daily average
        ELSEIF (itsdiag.EQ.5) THEN
          nanf=ldtdayc
          nend=ndtday

! Monthly average
        ELSEIF (itsdiag.EQ.6) THEN
          IF (ldays.EQ.1) THEN
             nanf=ldtdayc
          ELSEIF (ldays.GT.1) THEN
             nanf=ldtdayc+((ldays-1)*ndtday)
          ENDIF

          nend=ndtday*monlen(lmonts)

! Yearly average
        ELSEIF (ITSDIAG.EQ.7) THEN
          IF (LMONTS.EQ.1) THEN
            IF (LDAYS.EQ.1) THEN
              nanf=ldtdayc
            ELSE
              nanf=ldtdayc+((ldays-1)*ndtday)
            ENDIF
          ELSE
            nanf=0
            lmonts1=lmonts-1
            DO I=1,LMONTS1
              nanf=nanf+monlen(i)
            ENDDO
            nanf=nanf*ndtday
            IF (ldays.EQ.1) THEN
              nanf=nanf+ldtdayc
            ELSE
              nanf=nanf+ldtdayc+((ldays-1)*ndtday)
            ENDIF
          ENDIF

          nend=0
          DO I=1,12
             nend=nend+monlen(i)
          ENDDO
          nend=nend*ndtday

        ELSE
          STOP 'Stop in diagnosis due to wrong parameter.'
        ENDIF

        IF (nanf.EQ.1) avg_timeser(:)=0.

        do n=1,lentimser    
          avg_timeser(n)=avg_timeser(n)+(timeser(n)-avg_timeser(n))/nanf
        enddo

!FR    Write in EXTRA or PSEUDO EXTRA Format:
!      EXTRA format allows to transpose the whole time series data.
        IF (nanf.EQ.nend) THEN
!
! #slo: if write-frequency higher one day, time specification is missing
! #slo: a time specification should be appended to the idate header
! #slo: needs update of FR's tsplot script
!
          idate=(lyears*10000) + (lmonts*100) + ldays 
          ione=1

          DO n=1,lentimser
            write_timeser(n) = avg_timeser(n)
          ENDDO
 
          IF ( ltstranspose ) then
!FR    EXTRA Format - this is used for plot-script "tsplot":
! #slo: write n arrays of length one (only) for all n codes
            DO n=1,lentimser
              WRITE(io_ou_f126) idate,n,ione,ione
              WRITE(io_ou_f126) write_timeser(n)
            ENDDO
          ELSE
!FR    PSEUDO-EXTRA Format:
! #slo: write one array of code-values of length lentimser (the number of codes)
            WRITE(io_ou_f126) idate,ione,ione,lentimser
            WRITE(io_ou_f126) write_timeser
          ENDIF
        ENDIF

      ENDIF

      DEALLOCATE(      timeser)
      DEALLOCATE(write_timeser)

end subroutine write_timeseries


  SUBROUTINE calc_moc

    USE mo_param1, ONLY :ie_g,je_g,ke,kep
 
    USE mo_parallel, ONLY :gather_arr
    
    USE mo_commo1, ONLY : wo,depto,weto,      &
         tiestw,dzw,dlxp_g,dlyp_g,alat_g,weto_g,    &
         lyears,ldays,lmonts,ldtdayc,dt,imocdiag

    USE mo_units, ONLY :  io_ou_mocg,io_ou_moca

    USE mo_mean, ONLY: calc_avgint

    IMPLICIT NONE

    REAL :: wo_g(ie_g,Je_g,kep)


    REAL :: zlat
    INTEGER :: i,j,k,jbrei,lbrei,l,lb,i4                   &
              ,nanf,nend
    !change type of fort.75 header variables to 8 Byte integer
    !(necessary for processing fort.75 with CDO versions > 1.4.1)
    INTEGER*8 :: i1,i2,i3,param180=180

    tmerc(:,:,:)=0.

    IF (IMOCDIAG.GE.1 .AND. IMOCDIAG.LE.4) THEN
       call calc_avgint(imocdiag,nanf,nend)
    ELSEIF (IMOCDIAG.EQ.0) THEN
       RETURN
    ELSE
       STOP 'Stop in calc_moc due to wrong parameter.'
    ENDIF

!      WRITE(0,*)'in calc moc :',nanf,nend

    IF (IMOCDIAG.NE.0) THEN

    DO k=1,kep
       CALL gather_arr(wo(:,:,k),wo_g(:,:,k),p_io )
    ENDDO


!          compute moc 

    IF (p_pe==p_io) THEN
    jbrei=3
    DO i=2,ie_g-1
       DO j=1,je_g
          IF(weto_g(i,j,1).EQ.1.)THEN
             !     1   suedpol
             !     180 nordpol
             lbrei=NINT(90.+alat_g(i,j))
             lbrei=MAX(lbrei,1)
             lbrei=MIN(lbrei,180)

             DO k=1,ke
                zlat=MIN(dlxp_g(i,j),dlyp_g(i,j))/(float(2*jbrei)*111111.)
                DO l=-jbrei,jbrei
                   lbrei=NINT(90.+alat_g(i,j)+float(l)*zlat)
                   lbrei=MAX(lbrei,1)
                   lbrei=MIN(lbrei,180)

                   tmerc(lbrei,1,k)=tmerc(lbrei,1,k)-dlxp_g(i,j)*dlyp_g(i,j) &
                        *wo_g(i,j,k)/float(2*jbrei+1)

                   IF((ibek_g(i,j).LE.5))THEN
                      tmerc(lbrei,2,k)=tmerc(lbrei,2,k)-dlxp_g(i,j)*dlyp_g(i,j)  &
                           *wo_g(i,j,k)/float(2*jbrei+1)

                   ENDIF
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO

    DO lb=179,1,-1
       DO k=1,kep
          tmerc(lb,1,k)=tmerc(lb+1,1,k)+tmerc(lb,1,k)
          tmerc(lb,2,k)=tmerc(lb+1,2,k)+tmerc(lb,2,k)
       ENDDO
    ENDDO

    IF (nanf.EQ.1) sum_tmerc(:,:,:)=0.

    sum_tmerc(:,:,:)=sum_tmerc(:,:,:)+(tmerc(:,:,:)-sum_tmerc(:,:,:))/nanf

    IF (nanf.EQ.nend) THEN
       i1=(lyears*10000)+(lmonts*100)+ldays 
       i2=100
       i4=180*kep
       DO k=1,kep
          i3=NINT(tiestw(k))   
          WRITE(io_ou_mocg)i1,i2,i3,param180
          WRITE(io_ou_mocg) sum_tmerc(:,1,k)
       ENDDO
    ENDIF
    IF (nanf.EQ.nend) THEN
       i1=(lyears*10000)+(lmonts*100)+ldays 
       i2=101
       i4=180*kep
       DO k=1,kep
          i3=NINT(tiestw(k))   
          WRITE(io_ou_moca)i1,i2,i3,param180
          WRITE(io_ou_moca) sum_tmerc(:,2,k)
       ENDDO
    ENDIF

    ENDIF ! p_pe==p_pio
    ENDIF	


  END SUBROUTINE calc_moc

SUBROUTINE calc_difi(idate)

  USE mo_kind
  USE mo_commoau1
  USE mo_units


!
!  calculates 2d-fields helpful for diagnostics of ocean fluxes
!  e.g. sea ice melting from below etc.
!  snapshots at the end of months
!  interface:
!  idate parameter for header of output
!

  REAL hice(ie,je),tdz(ie,je),sdz(ie,je)
  INTEGER idate 
  INTEGER (kind=sp) i41,i42,i43,i44
  REAL (kind=sp) ff_4(ie_g,je_g)
  INTEGER i,j,k
  REAL       ff_g(ie_g,je_g)
!
!  initialize fields

  hice(:,:)=0.
  tdz(:,:)=0.
  sdz(:,:)=0.

  DO k=1,ke
   DO j=1,je
   DO i=1,ie
   tdz(i,j)=tdz(i,j)+ddpo(i,j,k)*tho(i,j,k)*weto(i,j,k)
   sdz(i,j)=sdz(i,j)+ddpo(i,j,k)*sao(i,j,k)*weto(i,j,k)
   ENDDO
  ENDDO
  ENDDO

  DO j=1,je
   DO i=1,ie
    hice(i,j)=(rhosnwa*sicsno(i,j)+rhoicwa*sictho(i,j))*weto(i,j,1)
    sdz(i,j)=sdz(i,j)+(sao(i,j,1)*(zo(i,j)-rhosnwa*sicsno(i,j))       &
            -(sao(i,j,1)-5.)*rhoicwa*sictho(i,j))*weto(i,j,1)
    tdz(i,j)=tdz(i,j)+tho(i,j,1)*(zo(i,j)-rhosnwa*sicsno(i,j)         &
                                -rhoicwa*sictho(i,j))*weto(i,j,1)
   ENDDO
  ENDDO


      CALL gather_arr(hice(:,:),ff_g,p_io)
      i41=idate
      i42=201
      i43=-100
      i44=ie_g*je_g


      IF(p_pe.EQ.p_io) THEN
       ff_4=ff_g
       WRITE(io_ou_difi)i41,i42,i43,i44
       WRITE(io_ou_difi)ff_4
      ENDIF

      CALL gather_arr(tdz,ff_g,p_io)
      i41=idate
      i42=202
      i43=-100
      i44=ie_g*je_g


      IF(p_pe.EQ.p_io) THEN
       ff_4=ff_g
       WRITE(io_ou_difi)i41,i42,i43,i44
       WRITE(io_ou_difi)ff_4
      ENDIF

      CALL gather_arr(sdz,ff_g,p_io)
      i41=idate
      i42=203
      i43=-100
      i44=ie_g*je_g


      IF(p_pe.EQ.p_io) THEN
       ff_4=ff_g
       WRITE(io_ou_difi)i41,i42,i43,i44
       WRITE(io_ou_difi)ff_4
      ENDIF



  END  SUBROUTINE calc_difi

  


END MODULE mo_diagnosis
