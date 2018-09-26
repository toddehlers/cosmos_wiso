      SUBROUTINE WRTE_MFL(kdtday,kdays,kmonts,kyears,kmean,nanf)
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

      CALL CALC_AVGINT(kmean,nanf,nend)

      if (kmean.ne.0) then

       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_UKOMFL            &
                   ,UKO,SUM_UKOMFL,3,1)
       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_VKEMFL            &
                   ,VKE,SUM_VKEMFL,4,2)

      ENDIF

      RETURN
      END
