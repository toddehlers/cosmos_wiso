MODULE mo_cbalone_memory
  USE mo_kind,             ONLY : dp

  implicit none

  private

  public :: cbal_offline_type,grid_offline_type,vegetation_offline_type,landcover_change_type
  public :: dynveg_clim_type, dynveg_offline_type
  public :: initializePools !! reads initial carbon pools from filepack
  public :: initVegetation  !! initializes the structure "vegetation" from file
  public :: initGrid        !! initializes the structure "grid" from file
  public :: initCbalance    !! initializes the structure "cbalance" (except carbon pools)
  public :: initLandCoverChange !! initializes the structure for landcover change emissions

  TYPE cbal_offline_type
     real(dp),pointer :: Cpool_green(:,:)    !! C-pool for leaves, fine roots, vegetative organs and other green (living) parts 
                                             !! .. of vegetation [mol(C)/m^2(canopy)]
     real(dp),pointer :: Cpool_reserve(:,:)  !! C-pool for carbohydrate reserve (sugars, starches) that allows plants to survive 
                                             !! .. bad times[mol(C)/m^2(canopy)]
     real(dp),pointer :: Cpool_woods(:,:)    !! C-pool for stems, thick roots and other (dead) structural material of living 
                                             !! .. plants [mol(C)/m^2(canopy)]
     real(dp),pointer :: Cpool_litter_leaf(:,:) !! C-pool for litter originating from leaves [mol(C)/m^2(canopy)]
     real(dp),pointer :: Cpool_litter_wood(:,:) !! C-pool for litter originating from woody parts of the plants [mol(C)/m^2(canopy)]
     real(dp),pointer :: Cpool_fast(:,:)     !! C-pool for below ground organic material that quickly decomposes partly[mol(C)/m^2(canopy)]
     real(dp),pointer :: Cpool_slow(:,:)     !! C-pool for below ground organic material (coming from the fast pool) [mol(C)/m^2(canopy)]

     real(dp),pointer :: soil_respiration(:,:)    !! mean daily rate of heterotrophic (soil) respiration [mol(CO2)/m^2(ground)] 
     real(dp),pointer :: NPP_flux_correction(:,:) !! Daily updated flux correction from yesterdays carbon balance [mol(CO2)/m^2(canopy) s]
     real(dp),pointer :: litter_flux(:,:)         !! That part of NPP that could not be stored in the carbon pools of the living plant
                                                  !! .. but had to be dropped into the fast pool

     real(dp),pointer :: LAI_yDayMean(:,:,:)            !! mean value of LAI yesterday (from  LAI_sum())
     real(dp),pointer :: LAI_previousDayMean(:,:)       !! mean value of LAI the day before yesterday
     real(dp),pointer :: NPP_yDayMean(:,:,:)            !! mean value of NPP-Rate yesterday (from NPP_sum()) [mol(CO2)/(m^2(canopy) s)]
     real(dp),pointer :: topSoilTemp_yDayMean(:,:,:)    !! mean value of upper layer soil temperature yesterday (from topSoilTemp_sum()) [K]
     real(dp),pointer :: alpha_yDayMean(:,:,:)          !! mean value of water stress coefficient alpha yesterday (from alpha_sum())

     real(dp),pointer :: box_Cpools_total(:)      !! Sum of all carbon pools [mol(C)/m^2(grid box)]
     real(dp),pointer :: box_NEP_wholeRun(:)      !! Net ecosystem productivity followed over the whole run [mol(C)/m^2(grid box)]

     real(dp),pointer :: avg_Cpool_green(:,:)     !! As Cpool_green() but time averaged and relative to grid box area [mol(C)/m^2(grid box)] 
     real(dp),pointer :: avg_Cpool_reserve(:,:)   !! As Cpool_reserve() but time averaged and relative to grid box area [mol(C)/m^2(grid box)]    
     real(dp),pointer :: avg_Cpool_woods(:,:)     !! As Cpool_woods() but time averaged and relative to grid box area [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_Cpool_litter_leaf(:,:) !! As Cpool_litter_leaf() but time averaged and relative to grid box area [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_Cpool_litter_wood(:,:)   !! As Cpool_litter_wood() but time averaged and relative to grid box area [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_Cpool_fast(:,:)          !! As Cpool_fast() but time averaged and relative to grid box area [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_Cpool_slow(:,:)          !! As Cpool_slow() but time averaged and relative to grid box area [mol(C)/m^2(grid box)

     real(dp),pointer :: avg_soil_respiration(:,:)    !! Average soil respiration rate [mol(CO2)/m^2(grid box) s]
     real(dp),pointer :: avg_NPP_yDayMean(:,:)        !! Average NPP rate [mol(CO2)/m^2(grid box) s]
     real(dp),pointer :: avg_NPP_flux_correction(:,:) !! Average flux correction for NPP [mol(CO2)/m^2(grid box) s]
     real(dp),pointer :: avg_litter_flux(:,:)         !! That part of NPP that could not be stored in the carbon pools of the living plant
                                                      !! .. but had to be dropped into the fast pool [mol(CO2)/m^2(grid box) s]
     real(dp),pointer :: avg_box_NEP(:)               !! Net ecosystem productivity [kg(CO2)/m^2(grid box) s]

  end TYPE cbal_offline_type

  TYPE vegetation_offline_type
     integer           :: ntiles         !! Number of tiles
     integer,pointer   :: coverType(:,:) !! Cover type for each grid point and tile
     real(dp),pointer  :: coverFract(:,:)!! Cover fraction of vegetated area for each grid cell and tile
     real(dp),pointer  :: vegRatioMax(:) !! Fraction of vegetated area in each grid cell
     logical,pointer   :: glacier(:)     !! Logical glacier mask
     logical,pointer   :: is_present(:,:) ! 
  end TYPE vegetation_offline_type

  TYPE grid_offline_type
     INTEGER           :: nlon, nlat                !! Number of global longitudes/latitudes
     INTEGER           :: nland                     !! Number of global land points
     REAL(dp), POINTER :: lon(:), lat(:)  !! Global longitudes and latitudes
     LOGICAL,  POINTER :: mask(:,:)       !! Global land mask
     REAL(dp), POINTER :: mask_real(:,:)  !! Real representation of mask (for netcdf output)
     INTEGER,  POINTER :: kpoints(:)      !! Index of land boxes into global grid

  END TYPE grid_offline_type

  TYPE landcover_change_type
     logical      :: do_landcover_change = .false.  !! Indicates whether landcover change shall be performed
     real(dp),pointer :: LCC_sum_box_C2atmos(:)        !! carbon released to atmosphere by landcover change [mol(C)/m^2(grid box)]
     real(dp),pointer :: LCC_sum_box_C2fastSoilPool(:) !! carbon relocated from green and reserve pool to fast soil pool by landcover 
                                                       !!    change [mol(C)/m^2(grid box)]
     real(dp),pointer :: LCC_sum_box_C2slowSoilPool(:) !! carbon relocated from wood pool to slow soil pool by landcover change 
                                                       !!    [mol(C)/m^2(grid box)]
     real(dp),pointer :: LCC_emissions_wholeRun(:)     !! Landcover emissions summed through whole run [kg(CO2)/m^2(grid box)]
  END TYPE landcover_change_type

  TYPE dynveg_clim_type                            !! input for dynamic vegetation (dimensions including time)
     REAL(dp), POINTER :: prev_year_npp(:,:,:) 
     REAL(dp), POINTER :: bio_exist(:,:,:)
     REAL(dp), POINTER :: rel_hum_air(:,:)
     REAL(dp), POINTER :: max_wind10(:,:)
     REAL(dp), POINTER :: prev_day_max_wind10(:,:)
  END TYPE dynveg_clim_type
  TYPE dynveg_offline_type                          !! output of dynamic vegetation
     REAL(dp), POINTER :: act_fpc(:,:)
     REAL(dp), POINTER :: pot_fpc(:,:)
     REAL(dp), POINTER :: burned_fpc(:,:)
     REAL(dp), POINTER :: damaged_fpc(:,:)
     REAL(dp), POINTER :: bare_fpc(:)
     REAL(dp), POINTER :: desert_fpc(:)
     REAL(dp), POINTER :: max_green_bio(:,:)
     REAL(dp), POINTER :: sum_green_bio_memory(:)
     REAL(dp), POINTER :: carbon_2_LeafLitterPool(:)
     REAL(dp), POINTER :: carbon_2_WoodLitterPool(:)
     REAL(dp), POINTER :: carbon_2_slowSoilPool_damage(:)
     REAL(dp), POINTER :: carbon_2_slowSoilPool_fire(:)
     REAL(dp), POINTER :: carbon_2_atmos(:)
  END TYPE dynveg_offline_type

  !! Parameters

  logical,parameter  ::  debug = .false.  !! If this is set to .true. debug information is printed out

contains

  ! --- initializePools() --------------------------------------------------------------------------------
  !
  ! Initializes 4 carbon pools from a file (any output file from Cbalone can be taken as initialization file)
  ! If carbon pools are read in from a file, they must be in LatLon format, packed format is currently not
  ! supported.
  !
  subroutine initializePools(grid,cbalance,vegetation,filename)
    USE mo_cbal_cpools,  ONLY: printCpoolParameters
    include 'netcdf.inc'
    type(grid_offline_type),intent(in)   :: grid
    type(cbal_offline_type),intent(inout):: cbalance
    type(vegetation_offline_type),intent(inout):: vegetation
    character(len=*),OPTIONAL,intent(in) :: filename

    !! --- internal variables

    real(dp),allocatable,dimension(:,:,:,:)  :: read_Cpools
    ! error status return
    integer  ::  iret
    ! netCDF id
    integer  ::  ncid
    ! dimension
    integer  ::  ndim,nvar,nattr,nunlimdim
    integer,allocatable,dimension(:)  ::  ndimlen
    character(len=50),allocatable,dimension(:)  ::  fdimnam
    ! variables
    integer ::  nin(100)
    integer,allocatable ,dimension(:)  ::  nvartyp
    integer,allocatable,dimension(:)   ::  nvardim
    integer,allocatable,dimension(:,:) ::  nvardimid
    integer,allocatable,dimension(:)   ::  nvaratt
    character(len=128),allocatable,dimension(:)  ::  fvarnam
    integer :: noTimeSteps
    ! indices
    integer  ::  i,ii,itile

    !! Printout parameters 

    call printCpoolParameters(6)

    !! --- initialize the pools
    
    if( .not. present(filename) ) then            !! Initialize from scratch
       cbalance%Cpool_green(:,:)   = 0.0_dp
       cbalance%Cpool_woods(:,:)   = 0.0_dp
       cbalance%Cpool_reserve(:,:) = 0.0_dp
       cbalance%Cpool_litter_leaf(:,:)   = 0.0_dp
       cbalance%Cpool_litter_wood(:,:)   = 0.0_dp
       cbalance%Cpool_fast(:,:)    = 0.0_dp
       cbalance%Cpool_slow(:,:)    = 0.0_dp
    else                                          !! Initialize from file
       !! open the file
       write (*,*) 'initializePools(): read Cpools from: ',filename
       iret = nf_open(filename,nf_nowrite,ncid)
       call check_err(iret)
       ! check what is in the input-file
       iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
       call check_err(iret)
       if (debug) write (*,*) 'initializePools(): ndim,nvar,nattr,nunlimdim',ndim,nvar,nattr,nunlimdim
       ! get the dimension name and length
       allocate (ndimlen(ndim))
       allocate (fdimnam(ndim))
       do i=1,ndim
          iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
          call check_err(iret)
          if (debug) write (*,*) 'initializePools(): i,fdimnam,ndimlen',i,trim(fdimnam(i)),ndimlen(i)
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
          do ii=1,nvardim(i)
             nvardimid(i,ii)=nin(ii)
          end do
          if (debug) write (*,*) 'initializePools(): i,fvarnam,nvartyp,nvardim,nvardimid,nvaratt',&
               i,trim(fvarnam(i)),nvartyp(i),nvardim(i),nvardimid(i,1),nvaratt(i)
       end do
       ! get data
       noTimeSteps = 0
       do i=1,nvar
          if (fvarnam(i) == 'Cpool_green') then
             allocate (read_Cpools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                write (0,*) 'initializePools(): dimensions Cpool_green wrong'
                stop
             endif
             iret = nf_get_var_real(ncid,i,read_Cpools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   write(*,*) "initializePools(): inconsistent length of time series for carbon pools in "//trim(filename)
                   stop
                end if
             end if
             do itile = 1,vegetation%ntiles
                cbalance%Cpool_green(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Cpools)
          else if(fvarnam(i) == 'Cpool_reserve') then
             allocate (read_Cpools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                write (0,*) 'initializePools(): dimensions Cpool_reserve wrong'
                stop
             endif
             iret = nf_get_var_real(ncid,i,read_Cpools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   write(*,*) "initializePools(): inconsistent length of time series for carbon pools in "//trim(filename)
                   stop
                end if
             end if
             do itile = 1,vegetation%ntiles
                cbalance%Cpool_reserve(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Cpools)
          else if(fvarnam(i) == 'Cpool_woods') then
             allocate (read_Cpools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                write (0,*) 'initializePools(): dimensions Cpool_woods wrong'
                stop
             endif
             iret = nf_get_var_real(ncid,i,read_Cpools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   write(*,*) "initializePools(): inconsistent length of time series for carbon pools in "//trim(filename)
                   stop
                end if
             end if
             do itile = 1,vegetation%ntiles
                cbalance%Cpool_woods(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Cpools)
          else if(fvarnam(i) == 'Cpool_litter_leaf') then
             allocate (read_Cpools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                write (0,*) 'initializePools(): dimensions Cpool_litter_leaf wrong'
                stop
             endif
             iret = nf_get_var_real(ncid,i,read_Cpools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   write(*,*) "initializePools(): inconsistent length of time series for carbon pools in "//trim(filename)
                   stop
                end if
             end if
             do itile = 1,vegetation%ntiles
                cbalance%Cpool_litter_leaf(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Cpools)
          else if(fvarnam(i) == 'Cpool_litter_wood') then
             allocate (read_Cpools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                write (0,*) 'initializePools(): dimensions Cpool_litter_wood wrong'
                stop
             endif
             iret = nf_get_var_real(ncid,i,read_Cpools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   write(*,*) "initializePools(): inconsistent length of time series for carbon pools in "//trim(filename)
                   stop
                end if
             end if
             do itile = 1,vegetation%ntiles
                cbalance%Cpool_litter_wood(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Cpools)
          else if(fvarnam(i) == 'Cpool_fast') then
             allocate (read_Cpools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                write (0,*) 'initializePools(): dimensions Cpool_fast wrong'
                stop
             endif
             iret = nf_get_var_real(ncid,i,read_Cpools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   write(*,*) "initializePools(): inconsistent length of time series for carbon pools in "//trim(filename)
                   stop
                end if
             end if
             do itile = 1,vegetation%ntiles
                cbalance%Cpool_fast(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Cpools)
          else if(fvarnam(i) == 'Cpool_slow') then
             allocate (read_Cpools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                write (0,*) 'initializePools(): dimensions Cpool_slow wrong'
                stop
             endif
             iret = nf_get_var_real(ncid,i,read_Cpools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   write(*,*) "initializePools(): inconsistent length of time series for carbon pools in "//trim(filename)
                   stop
                end if
             end if
             do itile = 1,vegetation%ntiles
                cbalance%Cpool_slow(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Cpools)
          end if
       end do
       write(*,*) "initializePools(): Carbon pools read in. Took values from last time step ",noTimeSteps," found in file."
       
       ! close the file
       iret = nf_close(NCID)
       call check_err(iret)
       ! deallocate
       deallocate(fdimnam,ndimlen)
       deallocate(fvarnam,nvartyp,nvardim,nvardimid,nvaratt)
    end if
    
  end subroutine initializePools

  ! --- initVegetation() -----------------------------------------------------------------------------------

  subroutine initVegetation(filename,grid,vegetation)
    include 'netcdf.inc'
    character(len=*),             intent(in)    :: filename
    TYPE(grid_offline_type),      intent(in)    :: grid
    TYPE(vegetation_offline_type),intent(inout) :: vegetation

    !! locals

    integer :: iret,ncid,ndim,nvar,nattr,nunlimdim
    integer,allocatable            :: ndimlen(:)
    character(len=128),allocatable :: fdimnam(:),fvarnam(:)
    integer, allocatable           :: nvartyp(:),nvardim(:),nvaratt(:)
    integer,allocatable            :: nvardimid(:,:)
    integer,allocatable            :: coverType(:,:,:)
    real(dp),allocatable           :: coverFract(:,:,:)
    real(dp),allocatable           :: vegRatioMax(:,:)
    integer,allocatable            :: glac(:,:)
    integer,allocatable            :: glacier01(:)
    integer :: i,ii,itile
    integer ::  nin(100)
    
    !! open the file
    write (*,*) 'initVegetation(): read vegetation data from: ',trim(filename)
    iret = nf_open(filename,nf_nowrite,ncid)
    call check_err(iret)
    ! check what is in the input-file
    iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
    call check_err(iret)
    if (debug) write (*,*) 'initVegetation(): ndim,nvar,nattr,nunlimdim',ndim,nvar,nattr,nunlimdim
    ! get the dimension name and length
    allocate (ndimlen(ndim))
    allocate (fdimnam(ndim))
    do i=1,ndim
       iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
       call check_err(iret)
       if (debug) write (*,*) 'initVegetation(): i,fdimnam,ndimlen',i,trim(fdimnam(i)),ndimlen(i)
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
       do ii=1,nvardim(i)
          nvardimid(i,ii)=nin(ii)
       end do
       if (debug) write (*,*) 'initVegetation(): i,fvarnam,nvartyp,nvardim,nvardimid,nvaratt',&
            i,trim(fvarnam(i)),nvartyp(i),nvardim(i),nvardimid(i,1),nvaratt(i)
    end do
    ! get ntiles
    vegetation%ntiles = 0
    do i=1,ndim
       if (fdimnam(i) == 'ntiles') vegetation%ntiles = ndimlen(i)
    end do
    if (vegetation%ntiles == 0) then
       write (0,*) 'initVegetation(): ntiles not found' 
       stop
    endif
    if(debug) write(*,*) "initVegetation(): ntiles=",vegetation%ntiles

    ! get cover type data
    do i=1,nvar
       if (fvarnam(i) == 'cover_type') then
          allocate (coverType(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
               ndimlen(nvardimid(i,3))))
          if (debug) write (*,*) 'initVegetation(): cover_type dimensions: ', &
               ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)),ndimlen(nvardimid(i,3))
          iret = nf_get_var_int(ncid,i,coverType)
          call check_err(iret)
       elseif (fvarnam(i) == 'cover_fract') then
          allocate (coverFract(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
               ndimlen(nvardimid(i,3))))
          if (debug) write (*,*) 'initVegetation(): cover_fract dimensions: ', &
               ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)),ndimlen(nvardimid(i,3))
          iret = nf_get_var_double(ncid,i,coverFract)
          call check_err(iret)
       elseif (fvarnam(i) == 'veg_ratio_max') then
          allocate (vegRatioMax(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2))))
          if (debug) write (*,*) 'initVegetation(): veg_ratio_max dimensions: ', &
               ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2))
          iret = nf_get_var_double(ncid,i,vegRatioMax)
          call check_err(iret)
       elseif (fvarnam(i) == 'glac') then
          allocate (glac(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2))))
          if (debug) write (*,*) 'initVegetation(): glac dimensions: ', &
               ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2))
          iret = nf_get_var_int(ncid,i,glac)
          call check_err(iret)
       end if
    end do

    ! pack cover type

    if(.not. allocated(coverType)) then
       write (0,*) 'Variable cover_type not found in input file'
       stop
    else
       allocate (vegetation%coverType(grid%nland,vegetation%ntiles))
       do itile=1,vegetation%ntiles
          vegetation%coverType(:,itile) = PACK(coverType(:,:,itile),mask=grid%mask)
       end do
       deallocate(coverType)
    end if
    
    ! pack cover fractions

    if(.not. allocated(coverFract)) then
       write (0,*) 'Variable cover_fract not found in input file'
       stop
    else
       allocate (vegetation%coverFract(grid%nland,vegetation%ntiles))
       do itile=1,vegetation%ntiles
          vegetation%coverFract(:,itile) = PACK(coverFract(:,:,itile),mask=grid%mask)
       end do
       deallocate(coverFract)
    end if

    ! pack vegetation fractions

    if(.not. allocated(vegRatioMax)) then
       write (0,*) 'Variable veg_ratio_max not found in input file'
       stop
    else
       allocate (vegetation%vegRatioMax(grid%nland))
       vegetation%vegRatioMax(:) = PACK(vegRatioMax(:,:),mask=grid%mask)
       deallocate(vegRatioMax)
    end if

    ! pack glacier mask

    if(.not. allocated(glac)) then
       write (0,*) 'Variable glac not found in input file'
       stop
    else
       allocate(glacier01(grid%nland))
       allocate(vegetation%glacier(grid%nland))
       glacier01(:) = PACK(glac(:,:),mask=grid%mask)
       WHERE (glacier01(:) == 1)
          vegetation%glacier = .TRUE.
       ELSEWHERE
          vegetation%glacier = .FALSE.
       END WHERE
       deallocate(glac)
       deallocate(glacier01)
    end if

    ! allocate is_present

    allocate (vegetation%is_present(grid%nland,vegetation%ntiles))
    vegetation%is_present = .TRUE.

    ! final deallocations

    deallocate (ndimlen)
    deallocate (fdimnam)
    deallocate (fvarnam)
    deallocate (nvartyp)
    deallocate (nvardim)
    deallocate (nvardimid)
    deallocate (nvaratt)

  end subroutine initVegetation

  ! --- initGrid() -----------------------------------------------------------------------------------

  subroutine initGrid(filename,grid)
    include 'netcdf.inc'
    character(len=*),intent(in)           :: filename
    type(grid_offline_type),intent(inout) :: grid

    double precision,allocatable,dimension(:,:) :: slm
    ! local variables
    ! error status return
    integer  ::  iret
    ! netCDF id
    integer  ::  ncid
    ! dimensions
    integer  ::  ndim,nvar,nattr,nunlimdim
    integer,allocatable,dimension(:)  ::  ndimlen
    character(len=50),allocatable,dimension(:)  ::  fdimnam
    ! variables
    integer ::  nin(100)
    integer,allocatable ,dimension(:)     ::  nvartyp
    integer,allocatable,dimension(:)      ::  nvardim
    integer,allocatable,dimension(:,:)    ::  nvardimid
    integer,allocatable,dimension(:)      ::  nvaratt
    character(len=128),allocatable,dimension(:)  ::  fvarnam
    ! indices
    integer  ::  i,ii,j,k,icount
    ! open the file
    if (debug) write (*,*) 'initGrid(): read grid information from: ',filename
    iret = nf_open(filename,nf_nowrite,ncid)
    call check_err(iret,"initGrid(): File "//trim(filename)//" not found.")
    ! check what is in the file
    iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
    call check_err(iret)
    if (debug) write (*,*) 'initGrid(): ndim,nvar,nattr,nunlimdim',ndim,nvar,nattr,nunlimdim
    ! get the dimension name and length
    allocate (ndimlen(ndim))
    allocate (fdimnam(ndim))
    do i=1,ndim
       iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
       call check_err(iret)
       if (debug) write (*,*) 'initGrid(): i,fdimnam,ndimlen',i,trim(fdimnam(i)),ndimlen(i)
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
       do ii=1,nvardim(i)
          nvardimid(i,ii)=nin(ii)
       end do
       if (debug) write (*,*) 'initGrid(): i,fvarnam,nvartyp,nvardim,nvardimid,nvaratt',&
            i,trim(fvarnam(i)),nvartyp(i),nvardim(i),nvardimid(i,1),nvaratt(i)
    end do
    ! get nlon, nlat, nland and mask (slm > 0.5)
    do i=1,nvar
       if (fvarnam(i) == 'slm') then
          allocate (slm(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2))))
          allocate (grid%mask(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2))))
          if (debug) write (*,*) 'initGrid(): reading land-sea-mask'
          iret = nf_get_var_double(ncid,i,slm)
          call check_err(iret)
          grid%nlon = ndimlen(nvardimid(i,1))
          grid%nlat = ndimlen(nvardimid(i,2))
          ! set nland
          icount = 0
          grid%mask(:,:) = .false.
          do j=1,grid%nlon
             do k=1,grid%nlat
                if (slm(j,k) > 0.5_dp) then
                   icount = icount + 1
                   grid%mask(j,k) = .true.
                end if
             end do
          end do
          grid%nland = icount
       endif
    end do
    if (.not. allocated(slm)) then
       write (0,*) 'initGrid(): slm not found - nland unknown'
       stop
    endif
    nullify(grid%lon)
    nullify(grid%lat)
    do i=1,nvar
       if (fvarnam(i) == 'lon') then
          allocate (grid%lon(grid%nlon))
          if(debug) write (*,*) 'initGrid(): reading longitudes'
          iret = nf_get_var_double(ncid,i,grid%lon(:))
          call check_err(iret)
       else if (fvarnam(i) == 'lat') then
          allocate (grid%lat(grid%nlat))    
          if(debug)  write (*,*) 'initGrid(): reading latitudes'
          iret = nf_get_var_double(ncid,i,grid%lat(:))
          call check_err(iret)
       end if
    end do
    if (.not. associated(grid%lon)) then
       write (0,*) 'initGrid(): lon not found'
       stop
    endif
    if (.not. associated(grid%lat)) then
       write (0,*) 'initGrid(): lat not found'
       stop
    endif
    if (debug) write (*,*) 'initGrid(): nland,nlon,nlat',grid%nland,grid%nlon,grid%nlat
    ! close the file
    iret = nf_close(NCID)
    call check_err(iret)
    ! deallocate local arrays
    deallocate(fdimnam,ndimlen)
    deallocate(fvarnam,nvartyp,nvardim,nvardimid,nvaratt)
    deallocate(slm)

  end subroutine initGrid

  ! --- initCbalance()  -----------------------------------------------------------------------------------

  subroutine initCbalance(grid,ntiles,cbalance)
    TYPE(grid_offline_type),intent(in) :: grid
    integer,intent(in)                 :: ntiles
    TYPE(cbal_offline_type),intent(out):: cbalance
    
    !! --- Allocate cpools

    allocate (cbalance%Cpool_green(grid%nland,ntiles))
    allocate (cbalance%Cpool_reserve(grid%nland,ntiles))
    allocate (cbalance%Cpool_woods(grid%nland,ntiles))
    allocate (cbalance%Cpool_litter_leaf(grid%nland,ntiles))
    allocate (cbalance%Cpool_litter_wood(grid%nland,ntiles))
    allocate (cbalance%Cpool_fast(grid%nland,ntiles))
    allocate (cbalance%Cpool_slow(grid%nland,ntiles))

    allocate (cbalance%soil_respiration(grid%nland,ntiles))
    allocate (cbalance%NPP_flux_correction(grid%nland,ntiles))
    allocate (cbalance%litter_flux(grid%nland,ntiles))

    allocate (cbalance%LAI_yDayMean(grid%nland,ntiles,31))
    allocate (cbalance%LAI_previousDayMean(grid%nland,ntiles))
    cbalance%LAI_previousDayMean(:,:) = 0.0_dp
    allocate (cbalance%NPP_yDayMean(grid%nland,ntiles,31))
    allocate (cbalance%topSoilTemp_yDayMean(grid%nland,ntiles,31))
    allocate (cbalance%alpha_yDayMean(grid%nland,ntiles,31))

    allocate (cbalance%box_Cpools_total(grid%nland))
    allocate (cbalance%box_NEP_wholeRun(grid%nland))

    allocate (cbalance%avg_Cpool_green(grid%nland,ntiles))
    allocate (cbalance%avg_Cpool_reserve(grid%nland,ntiles))
    allocate (cbalance%avg_Cpool_woods(grid%nland,ntiles))
    allocate (cbalance%avg_Cpool_litter_leaf(grid%nland,ntiles))
    allocate (cbalance%avg_Cpool_litter_wood(grid%nland,ntiles))
    allocate (cbalance%avg_Cpool_fast(grid%nland,ntiles))
    allocate (cbalance%avg_Cpool_slow(grid%nland,ntiles))

    allocate (cbalance%avg_soil_respiration(grid%nland,ntiles))
    allocate (cbalance%avg_NPP_yDayMean(grid%nland,ntiles))
    allocate (cbalance%avg_NPP_flux_correction(grid%nland,ntiles))
    allocate (cbalance%avg_litter_flux(grid%nland,ntiles))
    allocate (cbalance%avg_box_NEP(grid%nland))

  end subroutine initCbalance

  subroutine initLandCoverChange(grid,LC_change)
    TYPE(grid_offline_type),    intent(in)    :: grid
    TYPE(landcover_change_type),intent(inout) :: LC_change
    allocate (LC_change%LCC_sum_box_C2atmos(grid%nland))
    LC_change%LCC_sum_box_C2atmos(:) = 0.0_dp
    allocate (LC_change%LCC_sum_box_C2fastSoilPool(grid%nland))
    LC_change%LCC_sum_box_C2fastSoilPool(:) = 0.0_dp
    allocate (LC_change%LCC_sum_box_C2slowSoilPool(grid%nland))
    LC_change%LCC_sum_box_C2slowSoilPool(:) = 0.0_dp
    allocate (LC_change%LCC_emissions_wholeRun(grid%nland))
    LC_change%LCC_emissions_wholeRun(:)=0.0_dp
  end subroutine initLandCoverChange

  ! --- check_err --------------------------------------------------------------------------------------------

  subroutine check_err(IRET,text)
    include 'netcdf.inc'
    integer,intent(in)          :: iret
    character(len=*),optional,intent(in) :: text
    if (iret /= NF_NOERR) then
       if(present(text)) write(*,*) 'check_err(): '//trim(text)
       write(*,*) 'check_err(): netcdf says: '//trim(nf_strerror(iret))
       stop
    end if
  end subroutine check_err

END MODULE mo_cbalone_memory
