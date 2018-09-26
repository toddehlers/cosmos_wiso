      MODULE MO_FLUXES1
      USE mo_kind, ONLY: dp
      USE MO_PARAM1
      IMPLICIT NONE
! ---------------------------------------------------------------------
!
!*    *COMMON* *FLUXES1* - Atmospheric fluxes on the ocean grid.
!
!     S. Legutke          *DKRZ*           30.05.97
!sv   S. Venzke           *MPI*            21.07.99
!sv   restricted to odd fields for c-Grid and added wind stress velocity
!
!*    VARIABLE        TYPE       PURPOSE.
!     --------        ----       --------
!     *AOFLTXW*       *REAL*    Zonal wind stress on water
!     *AOFLTXI*       *REAL*    Zonal wind stress on snow/ice
!     *AOFLTYW*       *REAL*    Meridional wind stress on water
!     *AOFLTYI*       *REAL*    Meridional wind stress on snow/ice
!     *AOFLFRW*       *REAL*    liquid freshwater flux (over water and ice)
!     *AOFLFRI*       *REAL*    solid freshwater flux (over ice only)
!     *AOFLRHI*       *REAL*    residual heat flux used to melt snow/ice
!     *AOFLNHW*       *REAL*    net heat flux over water
!     *AOFLSHW*       *REAL*    downwelling solar radiation
!     *AOFLCHI*       *REAL*    conductive heat flux through ice
!sv   *AOFLWSV        *REAL*    wind stress velocity
!
!*    PURPOSE
!     -------
!     The ocean can be driven either by fluxes (gpp option FLUXES, used
!     with the coupled model) or by surface variables and fluxes. 
!     In the first case the forcing fields reside in COMMON block FLUXES, 
!     in the second case, they are in COMMON block OCEVAL.
!    Note: if (Prec.-evap.(over ice))*comp. is downward and the atmos.
!     temperature at the blending height is above 0 (i.e. it is raining), 
!     then this is added to AOFLFRW before being passed to the ocean 
!     This saves transfer of 1 field. AOFLFRW also includes runoff.
!    Note: downwelling solar radiation is multiplied by surface albedo in
!     the atmosphere model already. This is not the case with the other version.
!    Note: all fluxes besides freshwater are per m^2, and not per grid cell, i.e.
!     they are not multiplied with ice compactness.
!    Note: the relaxation heat flux is declared in COMMON FLUXES2.h
!
! ----------------------------------------------------------------------
!
#ifdef __coupled
      REAL, POINTER :: &
     &AOFLTXWO(:,:),AOFLTYWE(:,:),                                  &
     &AOFLTXIO(:,:),AOFLTYIE(:,:),                                  &
     &AOFLFRIO(:,:),AOFLFRWO(:,:),                                  &
     &AOFLRHIO(:,:),AOFLCHIO(:,:),                                  &
     &AOFLNHWO(:,:),AOFLSHWO(:,:),                                  &
     &AOFLWSVO(:,:),AOFLDHWO(:,:)
#ifdef ADDCONTRA
      REAL, POINTER :: &
     &AOFLFRWO16(:,:),AOFLFRWO18(:,:),AOFLFRWHDO(:,:),              &
     &AOFLFRIO16(:,:),AOFLFRIO18(:,:),AOFLFRIHDO(:,:)
#endif
#endif
      CONTAINS

      SUBROUTINE alloc_mem_fluxes1
!
! Allocate memory for arrays used in coupling, which are assigned in mo_couple
! nsk, 03.09.2004

#ifdef __coupled
        ALLOCATE(AOFLTXWO(IE,JE),AOFLTYWE(IE,JE),                       &
     &AOFLTXIO(IE,JE),AOFLTYIE(IE,JE),                                  &
     &AOFLFRIO(IE,JE),AOFLFRWO(IE,JE),                                  &
     &AOFLRHIO(IE,JE),AOFLCHIO(IE,JE),                                  &
     &AOFLNHWO(IE,JE),AOFLSHWO(IE,JE),                                  &
     &AOFLWSVO(IE,JE),AOFLDHWO(IE,JE))

        AOFLTXWO(:,:) = 0.0_dp
        AOFLTYWE(:,:) = 0.0_dp
        AOFLTXIO(:,:) = 0.0_dp
        AOFLTYIE(:,:) = 0.0_dp
        AOFLFRIO(:,:) = 0.0_dp
        AOFLFRWO(:,:) = 0.0_dp
        AOFLRHIO(:,:) = 0.0_dp
        AOFLCHIO(:,:) = 0.0_dp
        AOFLNHWO(:,:) = 0.0_dp
        AOFLSHWO(:,:) = 0.0_dp
        AOFLWSVO(:,:) = 0.0_dp
        AOFLDHWO(:,:) = 0.0_dp

#ifdef ADDCONTRA
        ALLOCATE(                                                       &
     &AOFLFRWO16(IE,JE),AOFLFRWO18(IE,JE),AOFLFRWHDO(IE,JE),            &
     &AOFLFRIO16(IE,JE),AOFLFRIO18(IE,JE),AOFLFRIHDO(IE,JE))

        AOFLFRWO16(:,:) = 0.0_dp
        AOFLFRWO18(:,:) = 0.0_dp
        AOFLFRWHDO(:,:) = 0.0_dp
        AOFLFRIO16(:,:) = 0.0_dp
        AOFLFRIO18(:,:) = 0.0_dp
        AOFLFRIHDO(:,:) = 0.0_dp
#endif

#endif

      END SUBROUTINE alloc_mem_fluxes1

      SUBROUTINE WRTE_FLUX_EXTRA
#ifdef __coupled
!C****************************************************************
!C
!C**** *WRTE_FLUX_EXTRA* - save atmosph. fluxes at coupled timestep
!
!C     JJ,    *MPI-Met, HH*    03.08.2003
!C     SL,    *MPI-Met, M&D*   08.08.2003
!C      - field name correction
!C     NSK,   *IFM-GEOMAR*     28.09.2004  
!C      - MPI version
!C
!C     Modified
!C     --------
!C
!C     Purpose
!C     -------
!C     sace atmosph flux fileds at each coupled timestep
!C
!C     Method
!C     -------
!C**   Interface.
!C     ----------
!C!
!C     *CALL*       *WRTE_MEAN(kdtday,kdays,kmonts,kmean)*
!C
!C     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!C     *COMMON*     *COMMO1.h*     - ocean/sediment tracer arrays.
!C     *COMMON*     *UNITS.h*      - std I/O logical units.
!C
!C     Externals
!C     ---------
!C     none.
!C
!C**************************************************************************
      USE MO_PARAM1
      USE MO_COMMO1

      USE MO_UNITS

      USE MO_PARALLEL
      USE MO_KIND

      INTEGER(KIND=i4) I4I1,I4I2,I4I3,I4I4,IDATE

!      REAL(wp), POINTER ::         glfld(:,:)
      REAL(wp), ALLOCATABLE ::         glfld(:,:)

      IDATE= (LYEARS*10000)+(LMONTS*100)+LDAYS
      i4i1=idate
      i4i2=270
      i4i3=0
      i4i4=ie_g*je_g
!     WRITE(IO_STDOUT,*) 'in sbr flux_extra' 
!     WRITE(IO_STDOUT,*) 'i4i1= ',i4i1
!     WRITE(IO_STDOUT,*) 'i4i2= ',i4i2
!     WRITE(IO_STDOUT,*) 'i4i3= ',i4i3
!     WRITE(IO_STDOUT,*) 'i4i4= ',i4i4
!

      IF ( p_parallel_io ) THEN
         ALLOCATE (glfld(ie_g,je_g))
!!$      ELSE
!!$         glfld => NULL()
      ENDIF

      call gather_arr(aoflnhwo,glfld,p_io)
      IF (p_parallel_io) THEN
         OPEN(IO_OU_AONHW, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AONHW) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AONHW) real(glfld,sp)
      CLOSE(IO_OU_AONHW)
      ENDIF
!
      call gather_arr(aoflshwo,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=271
         OPEN(IO_OU_AOSHW, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOSHW) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOSHW) real(glfld,sp)
      CLOSE(IO_OU_AOSHW)
      ENDIF
!
      call gather_arr(aoflrhio,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=272
         OPEN(IO_OU_AORHI, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AORHI) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AORHI) real(glfld,sp)
      CLOSE(IO_OU_AORHI)
      ENDIF
!
      call gather_arr(aoflchio,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=273
         OPEN(IO_OU_AOCHI, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOCHI) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOCHI) real(glfld,sp)
      CLOSE(IO_OU_AOCHI)
      ENDIF
!
      call gather_arr(aoflfrwo,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=274
         OPEN(IO_OU_AOFRW, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOFRW) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOFRW) real(glfld,sp)
      CLOSE(IO_OU_AOFRW)
      ENDIF
!
      call gather_arr(aoflfrio,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=275
         OPEN(IO_OU_AOFRI, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOFRI) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOFRI) real(glfld,sp)
      CLOSE(IO_OU_AOFRI)
      ENDIF
!
      call gather_arr(aofltxwo,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=276
         OPEN(IO_OU_AOTXW, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOTXW) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOTXW) real(glfld,sp)
      CLOSE(IO_OU_AOTXW)
      ENDIF
!
      call gather_arr(aofltywe,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=277
         OPEN(IO_OU_AOTYW, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOTYW) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOTYW) real(glfld,sp)
      CLOSE(IO_OU_AOTYW)
      ENDIF
!
      call gather_arr(aofltxio,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=278
         OPEN(IO_OU_AOTXI, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOTXI) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOTXI) real(glfld,sp)
      CLOSE(IO_OU_AOTXI)
      ENDIF
!
      call gather_arr(aofltyie,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=279
         OPEN(IO_OU_AOTYI, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOTYI) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOTYI) real(glfld,sp)
      CLOSE(IO_OU_AOTYI)
      ENDIF
!
      call gather_arr(aoflwsvo,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=280
         OPEN(IO_OU_AOWSV, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOWSV) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOWSV) real(glfld,sp)
      CLOSE(IO_OU_AOWSV)
      ENDIF
!      
#ifdef ADDCONTRA
!
      call gather_arr(aoflfrwo16,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=281
         OPEN(IO_OU_AOFRWO16, STATUS='UNKNOWN',                         &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOFRWO16) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOFRWO16) real(glfld,sp)
      CLOSE(IO_OU_AOFRWO16)
      ENDIF
!
      call gather_arr(aoflfrwo18,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=282
         OPEN(IO_OU_AOFRWO18, STATUS='UNKNOWN',                         &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOFRWO18) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOFRWO18) real(glfld,sp)
      CLOSE(IO_OU_AOFRWO18)
      ENDIF
!
      call gather_arr(aoflfrwhdo,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=283
         OPEN(IO_OU_AOFRWHDO, STATUS='UNKNOWN',                         &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOFRWHDO) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOFRWHDO) real(glfld,sp)
      CLOSE(IO_OU_AOFRWHDO)
      ENDIF
!
      call gather_arr(aoflfrio16,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=281
         OPEN(IO_OU_AOFRIO16, STATUS='UNKNOWN',                         &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOFRIO16) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOFRIO16) real(glfld,sp)
      CLOSE(IO_OU_AOFRIO16)
      ENDIF
!
      call gather_arr(aoflfrio18,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=282
         OPEN(IO_OU_AOFRIO18, STATUS='UNKNOWN',                         &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOFRIO18) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOFRIO18) real(glfld,sp)
      CLOSE(IO_OU_AOFRIO18)
      ENDIF
!
      call gather_arr(aoflfrihdo,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=283
         OPEN(IO_OU_AOFRIHDO, STATUS='UNKNOWN',                         &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOFRIHDO) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOFRIHDO) real(glfld,sp)
      CLOSE(IO_OU_AOFRIHDO)
      ENDIF
!
#endif

      IF (p_parallel_io) THEN      
         DEALLOCATE(glfld)
      ENDIF
!
      RETURN
#endif
      END SUBROUTINE WRTE_FLUX_EXTRA

      END MODULE MO_FLUXES1
