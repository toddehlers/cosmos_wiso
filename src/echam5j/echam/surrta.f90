SUBROUTINE surrta

  USE mo_mpi, ONLY: p_parallel_io, p_bcast, p_io

  !============================================================================
  !
  ! Description:
  !
  !   Reads data for mo_rrtaN (N=1:16) from formatted fortran file rrtadata
  !
  !   M.A. Giorgetta, MPI, June 2000
  !   L. Kornblueh, MPI, November 2001, fix read in parallel case
  !
  !============================================================================


  IMPLICIT NONE

  CHARACTER(23), PARAMETER :: readform ='(/,  /,(tr1,5(es16.9)))'
  INTEGER                  :: i

  ! read data from formatted file 'rrtadata'
  ! ----------------------------------------

  IF (p_parallel_io) THEN
    OPEN(11,file='rrtadata',status='old',form='formatted',action='read')

    DO i=1,12
      READ(11,*)
    END DO
  END IF

  CALL read_rrta1
  CALL read_rrta2
  CALL read_rrta3
  CALL read_rrta4
  CALL read_rrta5
  CALL read_rrta6
  CALL read_rrta7
  CALL read_rrta8
  CALL read_rrta9
  CALL read_rrta10
  CALL read_rrta11
  CALL read_rrta12
  CALL read_rrta13
  CALL read_rrta14
  CALL read_rrta15
  CALL read_rrta16

  CLOSE(11)

CONTAINS

  SUBROUTINE read_rrta1
    USE mo_rrta1 ,ONLY:absa,absb,fracrefa,fracrefb,forref,selfref
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) absb(:,:)
    CALL p_bcast(absb, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefb(:)
    CALL p_bcast(fracrefb, p_io)
    IF (p_parallel_io) READ(11,readform) forref(:)
    CALL p_bcast(forref, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)

!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(ABSB)
!CDIR DU_UPDATE(SELFREF)

  END SUBROUTINE read_rrta1

  SUBROUTINE read_rrta2
    USE mo_rrta2 ,ONLY:absa,absb,fracrefa,fracrefb,forref,selfref,refparam
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) absb(:,:)
    CALL p_bcast(absb, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:,:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefb(:)
    CALL p_bcast(fracrefb, p_io)
    IF (p_parallel_io) READ(11,readform) forref(:)
    CALL p_bcast(forref, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)
    IF (p_parallel_io) READ(11,readform) refparam(:)
    CALL p_bcast(refparam, p_io)
!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(ABSB)
!CDIR DU_UPDATE(SELFREF)
!CDIR DU_UPDATE(FRACREFA)

  END SUBROUTINE read_rrta2

  SUBROUTINE read_rrta3
    USE mo_rrta3 ,ONLY:absa,absb,fracrefa,fracrefb,forref,selfref,&
         &absn2oa,absn2ob,etaref,h2oref,n2oref,co2ref,strrat
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) absb(:,:)
    CALL p_bcast(absb, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:,:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefb(:,:)
    CALL p_bcast(fracrefb, p_io)
    IF (p_parallel_io) READ(11,readform) forref(:)
    CALL p_bcast(forref, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)
    IF (p_parallel_io) READ(11,readform) absn2oa(:)
    CALL p_bcast(absn2oa, p_io)
    IF (p_parallel_io) READ(11,readform) absn2ob(:)
    CALL p_bcast(absn2ob, p_io)
    IF (p_parallel_io) READ(11,readform) etaref(:)
    CALL p_bcast(etaref, p_io)
    IF (p_parallel_io) READ(11,readform) h2oref(:)
    CALL p_bcast(h2oref, p_io)
    IF (p_parallel_io) READ(11,readform) n2oref(:)
    CALL p_bcast(n2oref, p_io)
    IF (p_parallel_io) READ(11,readform) co2ref(:)
    CALL p_bcast(co2ref, p_io)
    IF (p_parallel_io) READ(11,readform) strrat
    CALL p_bcast(strrat, p_io)
!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(ABSB)
!CDIR DU_UPDATE(SELFREF)
!CDIR DU_UPDATE(FRACREFA)
!CDIR DU_UPDATE(FRACREFB)
!CDIR DU_UPDATE(H2OREF)
!CDIR DU_UPDATE(N2OREF)
!CDIR DU_UPDATE(CO2REF)
!CDIR DU_UPDATE(ETAREF)
  END SUBROUTINE read_rrta3

  SUBROUTINE read_rrta4
    USE mo_rrta4 ,ONLY:absa,absb,fracrefa,fracrefb,selfref,strrat1,strrat2
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) absb(:,:)
    CALL p_bcast(absb, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:,:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefb(:,:)
    CALL p_bcast(fracrefb, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)
    IF (p_parallel_io) READ(11,readform) strrat1
    CALL p_bcast(strrat1, p_io)
    IF (p_parallel_io) READ(11,readform) strrat2
    CALL p_bcast(strrat2, p_io)
!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(ABSB)
!CDIR DU_UPDATE(SELFREF)
!CDIR DU_UPDATE(FRACREFA)
!CDIR DU_UPDATE(FRACREFB)

  END SUBROUTINE read_rrta4

  SUBROUTINE read_rrta5
    USE mo_rrta5 ,ONLY:absa,absb,ccl4,fracrefa,fracrefb,selfref,strrat1,strrat2
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) absb(:,:)
    CALL p_bcast(absb, p_io)
    IF (p_parallel_io) READ(11,readform) ccl4(:)
    CALL p_bcast(ccl4, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:,:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefb(:,:)
    CALL p_bcast(fracrefb, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)
    IF (p_parallel_io) READ(11,readform) strrat1
    CALL p_bcast(strrat1, p_io)
    IF (p_parallel_io) READ(11,readform) strrat2
    CALL p_bcast(strrat2, p_io)
!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(ABSB)
!CDIR DU_UPDATE(SELFREF)
!CDIR DU_UPDATE(FRACREFA)
!CDIR DU_UPDATE(FRACREFB)

  END SUBROUTINE read_rrta5

  SUBROUTINE read_rrta6
    USE mo_rrta6 ,ONLY:absa,absco2,cfc11adj,cfc12,fracrefa,selfref
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) absco2(:)
    CALL p_bcast(absco2, p_io)
    IF (p_parallel_io) READ(11,readform) cfc11adj(:)
    CALL p_bcast(cfc11adj, p_io)
    IF (p_parallel_io) READ(11,readform) cfc12(:)
    CALL p_bcast(cfc12, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)

!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(SELFREF)

  END SUBROUTINE read_rrta6

  SUBROUTINE read_rrta7
    USE mo_rrta7 ,ONLY:absa,absb,absco2,fracrefa,fracrefb,selfref,strrat
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) absb(:,:)
    CALL p_bcast(absb, p_io)
    IF (p_parallel_io) READ(11,readform) absco2(:)
    CALL p_bcast(absco2, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:,:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefb(:)
    CALL p_bcast(fracrefb, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)
    IF (p_parallel_io) READ(11,readform) strrat
    CALL p_bcast(strrat, p_io)
!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(ABSB)
!CDIR DU_UPDATE(SELFREF)
!CDIR DU_UPDATE(FRACREFA)
!CDIR DU_UPDATE(FRACREFB)

  END SUBROUTINE read_rrta7

  SUBROUTINE read_rrta8
    USE mo_rrta8 ,ONLY:absa,absb,fracrefa,fracrefb,selfref,absco2a,absco2b,&
         &absn2oa,absn2ob,cfc12,cfc22adj,h2oref,n2oref,o3ref
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) absb(:,:)
    CALL p_bcast(absb, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefb(:)
    CALL p_bcast(fracrefb, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)
    IF (p_parallel_io) READ(11,readform) absco2a(:)
    CALL p_bcast(absco2a, p_io)
    IF (p_parallel_io) READ(11,readform) absco2b(:)
    CALL p_bcast(absco2b, p_io)
    IF (p_parallel_io) READ(11,readform) absn2oa(:)
    CALL p_bcast(absn2oa, p_io)
    IF (p_parallel_io) READ(11,readform) absn2ob(:)
    CALL p_bcast(absn2ob, p_io)
    IF (p_parallel_io) READ(11,readform) cfc12(:)
    CALL p_bcast(cfc12, p_io)
    IF (p_parallel_io) READ(11,readform) cfc22adj(:)
    CALL p_bcast(cfc22adj, p_io)
    IF (p_parallel_io) READ(11,readform) h2oref(:)
    CALL p_bcast(h2oref, p_io)
    IF (p_parallel_io) READ(11,readform) n2oref(:)
    CALL p_bcast(n2oref, p_io)
    IF (p_parallel_io) READ(11,readform) o3ref(:)
    CALL p_bcast(o3ref, p_io)

!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(ABSB)
!CDIR DU_UPDATE(SELFREF)
!CDIR DU_UPDATE(H2OREF)
!CDIR DU_UPDATE(N2OREF)
!CDIR DU_UPDATE(O3REF)

  END SUBROUTINE read_rrta8

  SUBROUTINE read_rrta9
    USE mo_rrta9 ,ONLY:absa,absb,fracrefa,fracrefb,selfref,absn2o,ch4ref,&
         &etaref,h2oref,n2oref,strrat
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) absb(:,:)
    CALL p_bcast(absb, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:,:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefb(:)
    CALL p_bcast(fracrefb, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)
    IF (p_parallel_io) READ(11,readform) absn2o(:)
    CALL p_bcast(absn2o, p_io)
    IF (p_parallel_io) READ(11,readform) ch4ref(:)
    CALL p_bcast(ch4ref, p_io)
    IF (p_parallel_io) READ(11,readform) etaref(:)
    CALL p_bcast(etaref, p_io)
    IF (p_parallel_io) READ(11,readform) h2oref(:)
    CALL p_bcast(h2oref, p_io)
    IF (p_parallel_io) READ(11,readform) n2oref(:)
    CALL p_bcast(n2oref, p_io)
    IF (p_parallel_io) READ(11,readform) strrat
    CALL p_bcast(strrat, p_io)
!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(ABSB)
!CDIR DU_UPDATE(SELFREF)
!CDIR DU_UPDATE(ABSN2O)
!CDIR DU_UPDATE(FRACREFA)
!CDIR DU_UPDATE(N2OREF)
!CDIR DU_UPDATE(H2OREF)
!CDIR DU_UPDATE(CH4REF)
!CDIR DU_UPDATE(ETAREF)

  END SUBROUTINE read_rrta9

  SUBROUTINE read_rrta10
    USE mo_rrta10,ONLY:absa,absb,fracrefa,fracrefb
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) absb(:,:)
    CALL p_bcast(absb, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefb(:)
    CALL p_bcast(fracrefb, p_io)

!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(ABSB)

  END SUBROUTINE read_rrta10

  SUBROUTINE read_rrta11
    USE mo_rrta11,ONLY:absa,absb,fracrefa,fracrefb,selfref
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) absb(:,:)
    CALL p_bcast(absb, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefb(:)
    CALL p_bcast(fracrefb, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)

!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(ABSB)
!CDIR DU_UPDATE(SELFREF)

  END SUBROUTINE read_rrta11

  SUBROUTINE read_rrta12
    USE mo_rrta12,ONLY:absa,fracrefa,selfref,strrat
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:,:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)
    IF (p_parallel_io) READ(11,readform) strrat
    CALL p_bcast(strrat, p_io)

!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(SELFREF)
!CDIR DU_UPDATE(FRACREFA)

  END SUBROUTINE read_rrta12

  SUBROUTINE read_rrta13
    USE mo_rrta13,ONLY:absa,fracrefa,selfref,strrat
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:,:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)
    IF (p_parallel_io) READ(11,readform) strrat
    CALL p_bcast(strrat, p_io)

!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(SELFREF)
!CDIR DU_UPDATE(FRACREFA)

  END SUBROUTINE read_rrta13

  SUBROUTINE read_rrta14
    USE mo_rrta14,ONLY:absa,absb,fracrefa,fracrefb,selfref
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) absb(:,:)
    CALL p_bcast(absb, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefb(:)
    CALL p_bcast(fracrefb, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)

!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(ABSB)
!CDIR DU_UPDATE(SELFREF)

  END SUBROUTINE read_rrta14

  SUBROUTINE read_rrta15
    USE mo_rrta15,ONLY:absa,fracrefa,selfref,strrat
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:,:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)
    IF (p_parallel_io) READ(11,readform) strrat
    CALL p_bcast(strrat, p_io)

!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(SELFREF)
!CDIR DU_UPDATE(FRACREFA)

  END SUBROUTINE read_rrta15

  SUBROUTINE read_rrta16
    USE mo_rrta16,ONLY:absa,fracrefa,selfref,strrat
    IF (p_parallel_io) READ(11,'(//)')
    IF (p_parallel_io) READ(11,readform) absa(:,:)
    CALL p_bcast(absa, p_io)
    IF (p_parallel_io) READ(11,readform) fracrefa(:,:)
    CALL p_bcast(fracrefa, p_io)
    IF (p_parallel_io) READ(11,readform) selfref(:,:)
    CALL p_bcast(selfref, p_io)
    IF (p_parallel_io) READ(11,readform) strrat
    CALL p_bcast(strrat, p_io)

!CDIR DU_UPDATE(ABSA)
!CDIR DU_UPDATE(SELFREF)
!CDIR DU_UPDATE(FRACREFA)

  END SUBROUTINE read_rrta16

END SUBROUTINE surrta
