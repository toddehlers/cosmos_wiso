      MODULE MO_THACC
      USE MO_PARAM1
      USE MO_MPI,             ONLY:  p_parallel_io
!-----------------------------------------------------------------------
!
!*   COMMON *THACC*  - mean (coupled time step) ocean/ice surface fields.
!
!
!    S.Legutke           *DKRZ*           25.06.97
!
!        - fields needed by the atmosphere to calculate fluxes.
!
!*   *THOACC*    -  mean SST of coupled time step on scalar grid.
!*   *SITOACC*   -  mean sea ice thickness     "  "
!*   *SNTOACC*   -  mean snow depth            "  "
!*   *SICOACC*   -  mean sea ice concentration "  "
!
!    Modified by
!    -----------
!    S.VENZKE                MPI          5.5.99
!    -  modfied for C-grid
!    S.Legutke               DKRZ           12.1.00
!    -  exchange fields declared only for cells which are exchanged
!    N. Keenlyside           IFM-GEOMAR   2.10.04
!    -  message passing version  (field accumulation on global arrays)
!  
!-----------------------------------------------------------------------
!
        INTEGER :: jpdim_ocei,jpdim_ocej,jpdim_oce
!

#ifdef __coupled
      REAL, POINTER :: AMSUE_G_L1(:,:),AMSUO_G_L1(:,:)

      REAL, POINTER :: &
     &             GL_SSTACC(:,:),GL_SITOACC(:,:)                       &
     &            ,GL_SICOACC(:,:),GL_SNTOACC (:,:)                     &
     &            ,gl_sst(:,:),gl_sicomo(:,:)                           &
     &            ,gl_sicsno(:,:),gl_sictho(:,:) 

      REAL, POINTER :: GL_SOCUACC(:,:),GL_SOCVACC(:,:)                  &
     &            ,gl_socu(:,:),gl_socv(:,:)                            &
     &            ,gl_sicu(:,:),gl_sicv(:,:)
#ifdef ADDCONTRA
      REAL, POINTER :: &
     &             gl_o16acc(:,:),gl_o18acc(:,:),gl_hdoacc(:,:)         &
     &            ,gl_o16(:,:),   gl_o18(:,:),   gl_hdo(:,:)
#endif      
#endif

      CONTAINS

      SUBROUTINE alloc_mem_thacc

        jpdim_ocei = ie_g-2 
        jpdim_ocej = je_g    
        jpdim_oce  = jpdim_ocej * jpdim_ocei 

!
! Allocate memory for arrays used in coupling, which are assigned in mo_couple
! nsk, 03.09.2004

#ifdef __coupled
      IF (p_parallel_io) THEN
        ALLOCATE(AMSUE_G_L1(IE_G,JE_G),AMSUO_G_L1(IE_G,JE_G))
        ALLOCATE(GL_SSTACC(jpdim_ocei,jpdim_ocej)                       &
     &          ,GL_SITOACC(jpdim_ocei,jpdim_ocej)                      &
     &          ,GL_SICOACC(jpdim_ocei,jpdim_ocej)                      &
     &          ,GL_SNTOACC(jpdim_ocei,jpdim_ocej)                      &
     &          ,gl_sst(ie_g,je_g),gl_sicomo(ie_g,je_g)                 &
     &          ,gl_sicsno(ie_g,je_g),gl_sictho(ie_g,je_g)) 
#ifdef ADDCONTRA
        ALLOCATE(gl_o16acc(jpdim_ocei,jpdim_ocej)                       &
     &          ,gl_o18acc(jpdim_ocei,jpdim_ocej)                       &
     &          ,gl_hdoacc(jpdim_ocei,jpdim_ocej)                       &
     &          ,gl_o16(ie_g,je_g)                                      &
     &          ,gl_o18(ie_g,je_g)                                      &
     &          ,gl_hdo(ie_g,je_g))
#endif
      ELSE
        AMSUE_G_L1 =>  NULL()
        AMSUO_G_L1 =>  NULL()
        GL_SSTACC =>  NULL() 
        GL_SITOACC => NULL()
        GL_SICOACC =>  NULL()
        GL_SNTOACC =>  NULL()
        gl_sst =>  NULL()
        gl_sicomo =>  NULL()
        gl_sicsno =>  NULL()
        gl_sictho =>  NULL()
#ifdef ADDCONTRA
        gl_o16acc =>  NULL() 
        gl_o18acc =>  NULL() 
        gl_hdoacc =>  NULL() 
        gl_o16    =>  NULL() 
        gl_o18    =>  NULL() 
        gl_hdo    =>  NULL() 
#endif
      ENDIF
       IF (p_parallel_io) THEN
          ALLOCATE(GL_SOCUACC(jpdim_ocei,jpdim_ocej)                    &
      &           ,GL_SOCVACC(jpdim_ocei,jpdim_ocej))
          ALLOCATE(gl_socu(IE_G,JE_G),gl_socv(IE_G,JE_G)                &
      &           ,gl_sicu(IE_G,JE_G),gl_sicv(IE_G,JE_G))
       ELSE
          GL_SOCUACC =>  NULL()
          GL_SOCVACC =>  NULL()
          gl_socu    =>  NULL()
          gl_socv    =>  NULL()
          gl_sicu    =>  NULL()
          gl_sicv    =>  NULL()
       ENDIF
#endif /* __coupled */
      END SUBROUTINE alloc_mem_thacc

      END MODULE MO_THACC




