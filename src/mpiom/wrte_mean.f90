      SUBROUTINE WRTE_MEAN(kdtday,kdays,kmonts,kyears,kmean,nanf)
!C****************************************************************
!C
!C**** *WRTE_MEAN* - save mean diagnostic output.
!
!C     CHH,    *MPI-Met, HH*   14.01.99
!C
!C     Modified
!C     --------
!C     S.Legutke,        *MPI-MaD, HH*    01.10.01
!C     - separate routine extracted from OLLIE (MAIN)
!C
!C     Purpose
!C     -------
!C     Accumulate fields, average, and save.
!C
!C     Method
!C     -------
!CHH   NO OUTPUT     : KMEAN=0
!!CHH   MONTLY AVERAGE: KMEAN=2
!CHH   YEARLY AVERAGE: KMEAN=3 
!C     
!C
!C**   Interface.
!C     ----------
!C!
!C     *CALL*       *WRTE_MEAN(kdtday,kdays,kmonts,kmean)*
!C
!C     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!C     *COMMON*     *COMMO1.h*     - ocean/sediment tracer arrays.
!C     *COMMON*     *UNITS.h*      - std I/O logical units.
!C
!C**   Interface to calling routine (parameter list):
!C     ----------------------------------------------
!C
!C     *INTEGER* *KYEARS*   - actual year.
!C     *INTEGER* *KMONT*   - actual month.
!C     *INTEGER* *KDAYS*    - actual day.
!C
!C
!C     Externals
!C     ---------
!C     none.
!C
!C**************************************************************************
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMOAU1
      USE MO_COMMOAU2
      USE MO_UNITS
      USE MO_MEAN
      USE MO_DIAGNOSIS
#ifdef __coupled
      USE MO_FLUXES1
      USE mo_couple,   ONLY: o2a_time, o2a_freq
#endif /*__coupled*/
      NDTDAY=NINT(86400./DT)

      DO I=1,IE1
        DO J=1,JE1
! Heat flux
          FLUM(I,J)=QSWO(I,J)+QLWO(I,J)+QLAO(I,J)+QSEO(I,J)
          PEM(I,J)=PRECH(I,J)+EMINPO(I,J)
! Sea ice transport (x-direction)
          SICTRU(I,J)=( (ABS(SICUO(I,J))+SICUO(I,J))                    &
     &                 *(SICTHO(I,J)+SICSNO(I,J)*RHOSNIC)               &
     &               +( SICUO(I,J)-ABS(SICUO(I,J)))                     &
     &                 *(SICTHO(I+1,J)+SICSNO(I+1,J)*RHOSNIC) )*0.5
! Sea ice transport (y-direction)
          SICTRV(I,J)=( (ABS(SICVE(I,J))+SICVE(I,J))                    &
     &                 *(SICTHO(I,J+1)+SICSNO(I,J+1)*RHOSNIC)           &
     &               +( SICVE(I,J)-ABS(SICVE(I,J)))                     &
     &                 *(SICTHO(I,J)+SICSNO(I,J)*RHOSNIC) )*0.5         
        ENDDO
      ENDDO
!

      CALL CALC_AVGINT(kmean,nanf,nend)

      if (kmean.ne.0) then


!HH 3D temporal means
       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MTHO            &
            ,THO,SUM_THO,2,0)
       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MSAO            &
            ,SAO,SUM_SAO,5,0)
!       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MUKO            &
!       ,UKO,SUM_UKO,3,1)
!       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MVKE            &
!       ,VKE,SUM_VKE,4,2)
       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MWO             &
            ,WO,SUM_WO,7,0)

       if (LGMDIAG) then
          CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MWGO            &
               ,WGO,SUM_WGO,207,0)
       endif


       if (LDIFFDIAG) then
          CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_WTMI            &
               ,WTMIX,SUM_WTMIX,112,0)
          CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_WTMI            &
               ,RINU,SUM_RINU,214,0)
          CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MAVO            &
               ,AVO,SUM_AVO,110,0)
          CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MDVO            &
               ,DVO,SUM_DVO,111,0)
       endif


!HH     2D ZEITMITTEL
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_EMIP           &
             ,EMINPO,SUM_EMINPO,67,0)
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MMZO           &
             ,ZO,SUM_ZO,1,0)
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_FLUM           &
             ,FLUM,SUM_FLUM,70,0)
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MPEM           &
             ,PEM,SUM_PEM,79,0)
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_SICT           &
             ,SICTHO,SUM_SICTHO,13,0)
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_SICO           &
             ,SICOMO,SUM_SICOMO,15,0)
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_SICU           &
             ,SICUO,SUM_SICUO,35,1)
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_SICV           &
             ,SICVE,SUM_SICVE,36,2)
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_SICS           &
             ,SICSNO,SUM_SICSNO,141,0)
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_SICTRU         &
             ,SICTRU,SUM_SICTRU,142,1)
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_SICTRV         &
             ,SICTRV,SUM_SICTRV,143,2)
#ifndef __coupled 
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_RIVRUN           &
     &              ,RIVRUN,SUM_RIVRUN,105,0)
#endif

        if (LGMDIAG) then
           CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_MBOLX          &
                ,BOLX,SUM_BOLX,159,2)
           CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_MBOLY          &
                ,BOLY,SUM_BOLY,160,2)
        endif

        if (LHFLDIAG) then
           CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_DQSWO            &
                ,DQSWO,SUM_DQSWO,247,0)
           CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_DQLWO            &
                ,DQLWO,SUM_DQLWO,248,0)
           CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_DQSEO            &
                ,DQSEO,SUM_DQSEO,249,0)
           CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_DQLAO            &
                ,DQLAO,SUM_DQLAO,250,0)
           CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_DQTHO            &
                ,DQTHO,SUM_DQTHO,251,0)
           CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_DQSWI            &
                ,DQSWI,SUM_DQSWI,252,0)
           CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_DQLWI            &
                ,DQLWI,SUM_DQLWI,253,0)
           CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_DQSEI            &
                ,DQSEI,SUM_DQSEI,254,0)
           CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_DQLAI            &
                ,DQLAI,SUM_DQLAI,246,0)
           CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_DQTHI            &
                ,DQTHI,SUM_DQTHI,245,0)
           CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_DTICEO           &
                ,DTICEO,SUM_DTICEO,212,0)
        endif

        if (LFORCEDIAG) then
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_TAFO            &
                ,TAFO,SUM_TAFO,92,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_FCLO            &
                ,FCLOU,SUM_FCLOU,164,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_FPRE            &
                ,FPREC,SUM_FPREC,213,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_FSWR            &
                ,FSWR,SUM_FSWR,80,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_FTDE            &
                ,FTDEW,SUM_FTDEW,81,0)


#ifdef __coupled 
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_TXO             &
                ,AOFLTXWO,SUM_TXO,52,1)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_TYE             &
                ,AOFLTYWE,SUM_TYE,53,2)
#else
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_TXO             &
                ,TXO,SUM_TXO,52,1)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_TYE             &
                ,TYE,SUM_TYE,53,2)
#endif
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_WU10            &
                ,FU10,SUM_WU10,171,0)
        endif


        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_QSWO            &
             ,QSWO,SUM_QSWO,176,0)
#ifndef __coupled
!       in __coupled setup qlwo = qlao = 0  
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_QLWO            &
             ,QLWO,SUM_QLWO,177,0)
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_QLAO            &
             ,QLAO,SUM_QLAO,147,0)
#endif /*__coupled*/

#ifndef __coupled
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_QSEO            &
             ,QSEO,SUM_QSEO,146,0)
#else
!       in __coupled setup qseo = net - sw hetaflux       
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_QSEO            &
             ,QSEO,SUM_QSEO,148,0)
#endif /*__coupled*/

        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_PREC            &
             ,PRECH,SUM_PRECH,65,0)	     
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_PSIU            &
             ,PSIUWE,SUM_PSIUWE,27,0)
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_ZMLD            &
             ,ZMLD,SUM_ZMLD,183,0)

        if (LCONVDIAG) then
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_TMCEO           &
                ,TMCEO,SUM_TMCEO,157,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_TMCDO           &
                ,TMCDO,SUM_TMCDO,158,0)
        endif

     endif


#ifdef __coupled
        CALL CALC_AVGINT(ioasisflux,nanf1,nend1)        

        if (IOASISFLUX.ne.0 .and. IOASISFLUX.ne.99) then
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF1,NEND1,IO_OU_AONHW           &
                ,AOFLNHWO,SUM_AOFLNHW,130,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF1,NEND1,IO_OU_AOSHW           &
                ,AOFLSHWO,SUM_AOFLSHW,131,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF1,NEND1,IO_OU_AORHI           &
                ,AOFLRHIO,SUM_AOFLRHI,132,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF1,NEND1,IO_OU_AOCHI           &
                ,AOFLCHIO,SUM_AOFLCHI,133,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF1,NEND1,IO_OU_AOFRW           &
                ,AOFLFRWO,SUM_AOFLFRW,134,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF1,NEND1,IO_OU_AOFRI           &
                ,AOFLFRIO,SUM_AOFLFRI,135,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF1,NEND1,IO_OU_AOTXW           &
                ,AOFLTXWO,SUM_AOFLTXW,136,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF1,NEND1,IO_OU_AOTYW           &
                ,AOFLTYWE,SUM_AOFLTYW,137,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF1,NEND1,IO_OU_AOTXI           &
                ,AOFLTXIO,SUM_AOFLTXI,138,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF1,NEND1,IO_OU_AOTYI           &
                ,AOFLTYIE,SUM_AOFLTYI,139,0)
           CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF1,NEND1,IO_OU_AOWSV           &
                ,AOFLWSVO,SUM_AOFLWSV,140,0)
        endif
#endif /*__coupled*/


      

      RETURN
      END
