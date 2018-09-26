      SUBROUTINE EXTR(kpie,kpje,field,pfldmsk,pmask,kunit)
!     ------------------------------------------------------------------------------
!     
!**** *EXTR* - computes extreme field entries at wet/dry cells.
!     
!     Purpose.
!     --------
!     Debugging and control.
!
!     Method
!     --------
!     pfldmsk is used as land/sea mask.
!     
!**   Interface.
!     ----------
!     *CALL* *EXTR(field,pfldmsk,pmask)
!     
!     
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *UNITS_BGC.h*  - std I/O logical units.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *REAL*   *field*   - field to be checked.
!     *REAL*   *pmask*   - test value at dry points.
!     *REAL*   *pfldmsk* - land sea mask [>0 for wet cells].
!     
!     ------------------------------------------------------------------------------

      USE mo_param1, ONLY: ie_g, je_g
      USE mo_mpi
      USE mo_parallel

      DIMENSION field  (kpie,kpje)
      DIMENSION pfldmsk(kpie,kpje)

      REAL field_g  (ie_g,je_g)
      REAL pfldmsk_g(ie_g,je_g)

      anzd  = 0.

      tmaxw = -1.e50
      tmaxd = -1.e50
      tminw = +1.e50
      tmind = +1.e50

      imaxw = 0
      imaxd = 0
      iminw = 0
      imind = 0

      jmaxw = 0
      jmaxd = 0
      jminw = 0
      jmind = 0

      CALL gather_arr(field,field_g,0)
      CALL p_bcast(field_g,0)
      CALL gather_arr(pfldmsk,pfldmsk_g,0)
      CALL p_bcast(pfldmsk_g,0)

      DO  j=1,je_g
         DO  i=1,ie_g

            IF(pfldmsk_g(i,j).LT.1.e-10) THEN

               IF ( tmaxd .LE. field_g(i,j) ) THEN
                  tmaxd = field_g(i,j)
                  imaxd = I
                  jmaxd = J
               ENDIF
               IF ( tmind .GE. field_g(i,j) ) THEN
                  tmind = field_g(i,j)
                  imind = I
                  jmind = J
               ENDIF
               IF( field_g(i,j).NE. pmask .AND. anzd.LT.1.) THEN
                  WRITE(kunit,*)                                        &
     &                'field .ne. pmask at dry cell i,j=('              &
     &                ,i,',',j,') : ',pfldmsk_g(i,j),field_g(i,j)
                  anzd = anzd + 1.
               ENDIF

            ELSE

               IF ( tmaxw .LE. field_g(i,j) ) THEN
                  tmaxw = field_g(i,j)
                  imaxw = I
                  jmaxw = J
               ENDIF
               IF ( tminw .GE. field_g(i,j) ) THEN
                  tminw = field_g(i,j)
                  iminw = I
                  jminw = J
               ENDIF

            ENDIF

         ENDDO
      ENDDO
      
      WRITE(kunit,1)'        Max. in wet  cells : ',tmaxw,            &
     &          ' at i,j= ',imaxw,jmaxw
      WRITE(kunit,1)'        Min. in wet  cells : ',tminw,            &
     &          ' at i,j= ',iminw,jminw                               

      IF( tmaxd.NE.pmask .OR. tmind.NE.pmask) THEN
         WRITE(kunit,1)'        Max. in dry  cells : ',tmaxd,           &
     &             ' at i,j= ',imaxd,jmaxd
         WRITE(kunit,1)'        Min. in dry  cells : ',tmind,           &
     &             ' at i,j= ',imind,jmind
      ENDIF
    1 FORMAT(a,e24.16,a,2i5)


      RETURN
      END
