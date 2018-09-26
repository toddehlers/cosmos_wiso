MODULE mo_cbalone_dynveg

  USE mo_kind,          ONLY: dp

  IMPLICIT none

  PUBLIC :: init_offline_dynveg
  PUBLIC :: update_offline_dynveg

  PRIVATE



CONTAINS

  SUBROUTINE init_offline_dynveg(grid, vegetation, dynveg, &
       dynveg_clim, rock_fract, read_fpc, filename, &
       dynveg_params, dynveg_options)

    USE mo_dynveg,         ONLY: cover_fract_to_fpc, dynveg_params_type, &
                                 dynveg_options_type, fpc_to_cover_fract, &
                                 scale_fpc
    USE mo_cbalone_memory, ONLY: grid_offline_type, dynveg_offline_type, &
                                 dynveg_clim_type, vegetation_offline_type
    USE mo_land_surface,   ONLY: scale_cover_fract

    INCLUDE 'netcdf.inc'
   
    TYPE(grid_offline_type),       INTENT(in)    :: grid
    TYPE(vegetation_offline_type), INTENT(inout) :: vegetation
    TYPE(dynveg_offline_type),     INTENT(inout) :: dynveg
    TYPE(dynveg_clim_type),        INTENT(inout) :: dynveg_clim
    TYPE(dynveg_params_type),      INTENT(in)    :: dynveg_params
    TYPE(dynveg_options_type),     INTENT(in)    :: dynveg_options
    REAL(dp),                      INTENT(inout) :: rock_fract(:)  
    LOGICAL,                       INTENT(in)    :: read_fpc
    CHARACTER(LEN=*),              INTENT(in)    :: filename

    REAL(dp), ALLOCATABLE :: array_read(:,:,:,:)
    REAL(dp), ALLOCATABLE :: sum_fpc(:)
    INTEGER,  ALLOCATABLE :: ndimlen(:), nvartyp(:), nvardim(:)
    INTEGER,  ALLOCATABLE :: nvardimid(:,:), nvaratt(:)
    CHARACTER(LEN=50), ALLOCATABLE :: fdimnam(:), fvarnam(:)
    INTEGER :: iret, ndim, nvar, i, noTimeSteps, idim, itile, ncid
    INTEGER :: nattr, nunlimdim
    INTEGER :: nin(100)

    ALLOCATE(dynveg%act_fpc(grid%nland,vegetation%ntiles))
    ALLOCATE(dynveg%pot_fpc(grid%nland,vegetation%ntiles))
    ALLOCATE(dynveg%burned_fpc(grid%nland,vegetation%ntiles))
    ALLOCATE(dynveg%damaged_fpc(grid%nland,vegetation%ntiles))
    ALLOCATE(dynveg%bare_fpc(grid%nland))
    ALLOCATE(dynveg%desert_fpc(grid%nland))
    ALLOCATE(dynveg%max_green_bio(grid%nland,vegetation%ntiles))
    ALLOCATE(dynveg%sum_green_bio_memory(grid%nland))
    dynveg%pot_fpc(:,:) = 0._dp
    dynveg%burned_fpc(:,:) = 0._dp
    dynveg%damaged_fpc(:,:) = 0._dp
    dynveg%max_green_bio(:,:) = 0._dp
    dynveg%sum_green_bio_memory(:) = 0._dp

    IF (dynveg_options%dynveg_feedback) THEN
       ALLOCATE(dynveg%carbon_2_LeafLitterPool(grid%nland))
       ALLOCATE(dynveg%carbon_2_WoodLitterPool(grid%nland))
       ALLOCATE(dynveg%carbon_2_slowSoilPool_damage(grid%nland))
       ALLOCATE(dynveg%carbon_2_slowSoilPool_fire(grid%nland))
       ALLOCATE(dynveg%carbon_2_atmos(grid%nland))
       dynveg%carbon_2_LeafLitterPool(:)=0._dp
       dynveg%carbon_2_WoodLitterPool(:)=0._dp
       dynveg%carbon_2_slowSoilPool_damage(:)=0._dp
       dynveg%carbon_2_slowSoilPool_fire(:)=0._dp
       dynveg%carbon_2_atmos(:)=0._dp
    END IF

    ALLOCATE(dynveg_clim%prev_year_npp(grid%nland,vegetation%ntiles,31))
    ALLOCATE(dynveg_clim%bio_exist(grid%nland,vegetation%ntiles,31))
    ALLOCATE(dynveg_clim%rel_hum_air(grid%nland,31))
    ALLOCATE(dynveg_clim%max_wind10(grid%nland,31))
    ALLOCATE(dynveg_clim%prev_day_max_wind10(grid%nland,31))
    dynveg_clim%bio_exist(:,:,:) = 1._dp

    IF (read_fpc) THEN    !! Read FPCs from file

       WRITE (*,*) 'Initial FPCs read from file ',filename
       iret = nf_open(filename,nf_nowrite,ncid)
       call check_err(iret)

       ! check what is in the input-file
       iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
       call check_err(iret)

       ! get the dimension name and length
       allocate (ndimlen(ndim))
       allocate (fdimnam(ndim))
       do i=1,ndim
          iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
          call check_err(iret)
       end do

       ! get variable names, types and shapes
       allocate (fvarnam(nvar))
       allocate (nvartyp(nvar))
       allocate (nvardim(nvar))
       allocate (nvardimid(nvar,100))
       allocate (nvaratt(nvar))
       do i=1,nvar
          iret = nf_inq_var(ncid,i,fvarnam(i),nvartyp(i),nvardim(i),nin,nvaratt(i))
          call check_err(iret)
          do idim=1,nvardim(i)
             nvardimid(i,idim)=nin(idim)
          end do
       end do

       ! get the data
       DO i=1,nvar
          IF (fvarnam(i) == 'act_fpc') THEN
             ALLOCATE(array_read(ndimlen(nvardimid(i,1)), ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)), ndimlen(nvardimid(i,4))))
             IF (ndimlen(nvardimid(i,1)) /= grid%nlon .OR. ndimlen(nvardimid(i,2)) /= grid%nlat) &
                  STOP 'init_offline_dynveg: dimensions of act_fpc are wrong'
             iret = NF_GET_VAR_DOUBLE(ncid,i,array_read(:,:,:,:))
             CALL check_err(iret)
             noTimeSteps = ndimlen(nvardimid(i,4))
             DO itile = 1,vegetation%ntiles
                dynveg%act_fpc(:,itile) = PACK(array_read(:,:,itile,noTimeSteps),MASK=grid%mask)
             END DO
             DEALLOCATE(array_read)
          END IF
          IF (fvarnam(i) == 'bare_fpc') THEN
             ALLOCATE(array_read(ndimlen(nvardimid(i,1)), ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)), ndimlen(nvardimid(i,4))))
             IF (ndimlen(nvardimid(i,1)) /= grid%nlon .OR. ndimlen(nvardimid(i,2)) /= grid%nlat) &
                  STOP 'init_offline_dynveg: dimensions of bare_fpc are wrong'
             iret = NF_GET_VAR_DOUBLE(ncid,i,array_read(:,:,:,:))
             CALL check_err(iret)
             noTimeSteps = ndimlen(nvardimid(i,4))
             dynveg%bare_fpc(:) = PACK(array_read(:,:,1,noTimeSteps),MASK=grid%mask)
             DEALLOCATE(array_read)
          END IF
          IF (fvarnam(i) == 'desert_fpc') THEN
             ALLOCATE(array_read(ndimlen(nvardimid(i,1)), ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)), ndimlen(nvardimid(i,4))))
             IF (ndimlen(nvardimid(i,1)) /= grid%nlon .OR. ndimlen(nvardimid(i,2)) /= grid%nlat) &
                  STOP 'init_offline_dynveg: dimensions of desert_fpc are wrong'
             iret = nf_get_var_double(ncid,i,array_read(:,:,:,:))
             CALL check_err(iret)
             noTimeSteps = ndimlen(nvardimid(i,4))
             dynveg%desert_fpc(:) = PACK(array_read(:,:,1,noTimeSteps),MASK=grid%mask)
             DEALLOCATE(array_read)
          END IF
       END DO
       WRITE(*,*) "init_offline_dynveg: FPCs read from file"

       ! close the file
       iret = nf_close(ncid)
       CALL check_err(iret)

       ! deallocate
       DEALLOCATE(fdimnam,ndimlen)
       DEALLOCATE(fvarnam,nvartyp,nvardim,nvardimid,nvaratt)

    ELSE                  !! Initialize FPCs from echam cover fractctions

       dynveg%act_fpc(:,:) = 0._dp
       dynveg%bare_fpc(:) = 0._dp
       dynveg%desert_fpc(:) = 0._dp
       CALL scale_cover_fract(grid%nland, vegetation%ntiles, vegetation%glacier(:), vegetation%coverFract(:,:))
       CALL cover_fract_to_fpc(grid%nland, vegetation%ntiles, vegetation%coverFract(:,:), &
            vegetation%vegRatioMax(:), rock_fract(:), vegetation%glacier(:), &
            dynveg_params%dynamic_PFT(:), dynveg%act_fpc(:,:), &
            dynveg%bare_fpc(:), dynveg%desert_fpc(:))

    END IF

    WHERE (vegetation%glacier(:))
       rock_fract(:) = 1._dp
    END WHERE
    DO itile=1,vegetation%ntiles
       WHERE (rock_fract(:) >= 1._dp - EPSILON(1._dp))
          dynveg%act_fpc(:,itile) = 0._dp
       END WHERE
    END DO

    CALL scale_fpc(grid%nland, vegetation%ntiles, vegetation%glacier(:), dynveg_params%dynamic_pft(:), &
            dynveg%act_fpc(:,:), dynveg%bare_fpc(:))

    IF (dynveg_options%dynveg_feedback) THEN

!   -- call fpc_to_cover_fract to get consistent veg_ratio_max
       CALL fpc_to_cover_fract(grid%nland, vegetation%ntiles, dynveg%act_fpc(:,:), dynveg%bare_fpc(:),               &
                               dynveg%desert_fpc(:), rock_fract(:), vegetation%glacier(:),                           &
                               dynveg_params%woody_pft(:), dynveg_params%dynamic_pft(:), vegetation%coverFract(:,:), &
                               vegetation%vegRatioMax(:))
    ENDIF

  END SUBROUTINE init_offline_dynveg

!------------------------------------------------------------------------------
  SUBROUTINE update_offline_dynveg(day, run_year, run_year_first, new_year, &
       grid, vegetation, dynveg, dynveg_clim, cbalance, dynveg_params, &
       dynveg_options, specificLeafArea_C, rock_fract, veg_fract_correction)
!------------------------------------------------------------------------------

    USE mo_dynveg,         ONLY: potential_tree_fpc, desert_fraction, &
                                 frac_wood_2_litter_wind, frac_wood_2_litter_fire, &
                                 frac_wood_2_atmos, dynveg_params_type, &
                                 dynveg_options_type, & 
                                 fpc_daily, fpc_to_cover_fract, scale_fpc
    USE mo_cbalone_memory, ONLY: grid_offline_type, dynveg_offline_type, &
                                 dynveg_clim_type, vegetation_offline_type, &
                                 cbal_offline_type
    USE mo_cbal_cpools,    ONLY: relocate_carbon, relocate_carbon_desert, &
                                 relocate_carbon_fire, relocate_carbon_damage
    USE mo_utils,          ONLY: average_tiles
    USE mo_land_surface,   ONLY: fract_small

    INTEGER,                       INTENT(in)    :: day
    INTEGER,                       INTENT(in)    :: run_year
    INTEGER,                       INTENT(in)    :: run_year_first
    LOGICAL,                       INTENT(in)    :: new_year
    TYPE(grid_offline_type),       INTENT(in)    :: grid
    TYPE(vegetation_offline_type), INTENT(inout) :: vegetation
    TYPE(dynveg_offline_type),     INTENT(inout) :: dynveg
    TYPE(dynveg_clim_type),        INTENT(in)    :: dynveg_clim
    TYPE(cbal_offline_type),       INTENT(inout) :: cbalance
    TYPE(dynveg_params_type),      INTENT(in)    :: dynveg_params
    TYPE(dynveg_options_type),     INTENT(in)    :: dynveg_options
    REAL(dp),                      INTENT(in)    :: specificLeafArea_C(:,:)
    REAL(dp),                      INTENT(in)    :: rock_fract(:)
    REAL(dp),                      INTENT(in)    :: veg_fract_correction(:,:)

    LOGICAL               :: first_year, second_year
    REAL(dp), ALLOCATABLE :: veg_ratio_max_old(:)
    REAL(dp), ALLOCATABLE :: cover_fract_old(:,:)
    REAL(dp)              :: act_fpc_old(grid%nland,vegetation%ntiles)
    REAL(dp)              :: cpool_litter(grid%nland) !! amount of litter [mol(C)/m^2(grid box)]
    INTEGER               :: itile

    first_year = (run_year == run_year_first)
    second_year = (run_year == run_year_first + 1)
!
!-- annual calculations
!
    IF (new_year) THEN
       IF (.NOT. first_year) THEN
          CALL potential_tree_fpc(grid%nland, vegetation%ntiles, dynveg_options%dynveg_all, & 
                                  dynveg_params%woody_pft(:), dynveg_params%dynamic_pft(:), &
                                  dynveg_clim%prev_year_npp(:,:,day), dynveg_clim%bio_exist(:,:,day), &
                                  dynveg%act_fpc(:,:), dynveg%pot_fpc(:,:))

          CALL desert_fraction(second_year, grid%nland, vegetation%ntiles,                       &
                               dynveg_params%woody_pft(:), dynveg_params%dynamic_pft(:),         &
                               dynveg%act_fpc(:,:), dynveg%bare_fpc(:), specificLeafArea_C(:,:), &
                               dynveg%max_green_bio(:,:), dynveg%sum_green_bio_memory(:),        &
                               dynveg%desert_fpc(:))

          IF (dynveg_options%dynveg_feedback) THEN

!         -- conversion of FPCs to JSBACH cover fractions and new vegetated fraction (veg_ratio_max)
             ALLOCATE(veg_ratio_max_old(grid%nland))
             ALLOCATE(cover_fract_old(grid%nland,vegetation%ntiles))
             veg_ratio_max_old(:) = vegetation%vegRatioMax(:)
             cover_fract_old(:,:) = vegetation%coverFract(:,:)
             CALL fpc_to_cover_fract(grid%nland, vegetation%ntiles, dynveg%act_fpc(:,:), dynveg%bare_fpc(:),               &
                                     dynveg%desert_fpc(:), rock_fract(:), vegetation%glacier(:),                           &
                                     dynveg_params%woody_pft(:), dynveg_params%dynamic_pft(:), vegetation%coverFract(:,:), &
                                     vegetation%vegRatioMax(:))

!         -- relocate carbon accouting for the change in the cover fraction of each PFT
             CALL relocate_carbon(cover_fract_old(:,:), vegetation%coverFract(:,:), veg_fract_correction(:,:), &
                                  fract_small, cbalance%Cpool_green(:,:), cbalance%Cpool_woods(:,:), &
                                  cbalance%Cpool_reserve(:,:), cbalance%Cpool_litter_leaf(:,:), &
                                  cbalance%Cpool_litter_wood(:,:), cbalance%Cpool_fast(:,:), &
                                  cbalance%Cpool_slow(:,:))

!         -- scaling cpools to account for change in vegetated fraction in order to conserve the total mass of carbon
             CALL relocate_carbon_desert(grid%nland, vegetation%ntiles, vegetation%vegRatioMax(:),                          &
                                         veg_ratio_max_old(:), fract_small,                                                 &
                                         cbalance%LAI_yDayMean(:,:,day), specificLeafArea_C(:,:),                           &
                                         cbalance%Cpool_green(:,:), cbalance%Cpool_woods(:,:), cbalance%Cpool_reserve(:,:), &
                                         cbalance%Cpool_litter_leaf(:,:), cbalance%Cpool_litter_wood(:,:),                  &
                                         cbalance%Cpool_fast(:,:), cbalance%Cpool_slow(:,:))

             DEALLOCATE(veg_ratio_max_old)
             DEALLOCATE(cover_fract_old)
          END IF
       END IF
    END IF
!
!-- daily calculations
!
!-- determine maximum of green biomass for the current year
    DO itile = 1,vegetation%ntiles
       IF (dynveg_params%dynamic_pft(itile)) THEN
          dynveg%max_green_bio(:,itile) = MAX(dynveg%max_green_bio(:,itile),cbalance%cpool_green(:,itile))
       END IF
    END DO

!-- determine grid box average litter
    IF (dynveg_options%dynveg_feedback) THEN
       CALL average_tiles((cbalance%cpool_litter_leaf(:,:) + cbalance%cpool_litter_wood(:,:)) * &
            veg_fract_correction(:,:) * SPREAD(vegetation%vegRatioMax(:),NCOPIES=vegetation%ntiles,DIM=2), &
            vegetation%is_present(:,:), vegetation%coverFract(:,:), cpool_litter(:))
       act_fpc_old(:,:) = dynveg%act_fpc(:,:)
    ELSE
       cpool_litter(:) = 0._dp
    ENDIF

!-- calculation of FPC (competition of vegetation types and disturbances (fire, wind break))
    CALL fpc_daily(grid%nland, vegetation%ntiles,  dynveg_options%dynveg_all, dynveg_params%woody_pft(:), &
                   dynveg_params%dynamic_pft(:), dynveg_params%tau_pft, &
                   dynveg_clim%prev_year_npp(:,:,day), dynveg_clim%bio_exist(:,:,day), &
                   dynveg%pot_fpc(:,:), dynveg_clim%rel_hum_air(:,day), &
                   dynveg_clim%prev_day_max_wind10(:,day), dynveg_clim%max_wind10(:,day), cpool_litter(:), &
                   dynveg%act_fpc(:,:), dynveg%burned_fpc(:,:), dynveg%damaged_fpc(:,:), dynveg%bare_fpc(:))

!-- rescaling of actual fpc and bare fpc
 
    CALL scale_fpc(grid%nland, vegetation%ntiles, vegetation%glacier(:), dynveg_params%dynamic_pft(:), &
                   dynveg%act_fpc(:,:), dynveg%bare_fpc(:))

    IF (dynveg_options%dynveg_feedback) THEN

       dynveg%carbon_2_WoodLitterPool(:) = 0._dp

!   -- changes in carbon pools caused by wind break
       CALL relocate_carbon_damage(grid%nland, vegetation%ntiles, act_fpc_old(:,:), dynveg%damaged_fpc(:,:), &
            vegetation%coverFract(:,:), fract_small, veg_fract_correction(:,:), frac_wood_2_litter_wind, &
            cbalance%Cpool_green(:,:), cbalance%Cpool_woods(:,:), cbalance%Cpool_reserve(:,:), &
            cbalance%Cpool_litter_leaf(:,:), cbalance%Cpool_litter_wood(:,:), cbalance%Cpool_slow(:,:), &
            dynveg%carbon_2_slowSoilPool_damage(:), dynveg%carbon_2_LeafLitterPool(:), &
            dynveg%carbon_2_WoodLitterPool(:))

!   -- changes in carbon pools caused by fire
       CALL relocate_carbon_fire(grid%nland, vegetation%ntiles, act_fpc_old(:,:), dynveg%burned_fpc(:,:), &
            vegetation%coverFract(:,:), fract_small, veg_fract_correction(:,:),  frac_wood_2_atmos, frac_wood_2_litter_fire, &
            cbalance%Cpool_green(:,:), cbalance%Cpool_woods(:,:), cbalance%Cpool_reserve(:,:), &
            cbalance%Cpool_litter_leaf(:,:), cbalance%Cpool_litter_wood(:,:), cbalance%Cpool_slow(:,:), &
            dynveg%carbon_2_WoodLitterPool(:), dynveg%carbon_2_slowSoilPool_fire(:), dynveg%carbon_2_atmos(:))

!   -- scale fluxes from vegetated area to whole grid box
       dynveg%carbon_2_slowSoilPool_damage(:) = dynveg%carbon_2_slowSoilPool_damage(:) * vegetation%vegRatioMax(:)
       dynveg%carbon_2_LeafLitterPool(:) = dynveg%carbon_2_LeafLitterPool(:) * vegetation%vegRatioMax(:)
       dynveg%carbon_2_WoodLitterPool(:) = dynveg%carbon_2_WoodLitterPool(:) * vegetation%vegRatioMax(:)
       dynveg%carbon_2_slowSoilPool_fire(:) = dynveg%carbon_2_slowSoilPool_fire(:) * vegetation%vegRatioMax(:)
       dynveg%carbon_2_atmos(:) = dynveg%carbon_2_atmos(:) * vegetation%vegRatioMax(:)

    END IF

  END SUBROUTINE update_offline_dynveg

!------------------------------------------------------------------------------
  SUBROUTINE check_err(iret,text)
!------------------------------------------------------------------------------

    INCLUDE 'netcdf.inc'

    INTEGER,                  INTENT(in) :: iret
    CHARACTER(len=*),OPTIONAL,INTENT(in) :: text

    IF (iret /= NF_NOERR) then
       IF (PRESENT(text)) WRITE(*,*) 'check_err(): '//TRIM(text)
       WRITE(*,*) 'check_err(): netcdf says: '//TRIM(nf_strerror(iret))
       STOP
    END IF
  END SUBROUTINE check_err

END MODULE mo_cbalone_dynveg
