module mo_cbalone_io

  USE mo_kind,          ONLY: dp
  implicit none

  public :: getDailyData
  public :: writeSingleTimeStep
  public :: read_landcover_fractions

  private

  logical,parameter  ::  debug = .false.

  TYPE netcdfInfo_type
     character(len=64) :: varName  = ""
     integer           :: varType  = -1   !! E.g. NF_FLOAT or NF_DOUBLE
     integer           :: dimType  = -1   !! BOX_TYPE or  TILES_TYPE (see below)  
     integer           :: varId    = -1   !! Is set when calling netcdf definition routine
     real(dp),pointer  :: array_box(:)     => NULL() !! For access to the data (only if NF_FLOAT  and BOX_TYPE)
     real(dp),pointer  :: array_tiles(:,:) => NULL() !! For access to box data (only if NF_FLOAT  and TILES_TYPE)
  end TYPE netcdfInfo_type

  !! For the definition of Netcdf dimensions: 
  integer,parameter :: BOX_TYPE = 3          !! lat x lon x time
  integer,parameter :: TILES_TYPE = 4        !! lat x lon x tiles x time

  !! Module variables

  integer :: numberOfOutputVariables
  TYPE(netcdfInfo_type),pointer :: netcdfInfo(:)

contains
  
  ! --- getDailyData()  -----------------------------------------------------------------------------------
  !
  ! This routine reads in a time series of daily forcing data, called LAI_yDayMean, NPP_yDayMean, topSoilTemp_yDayMean
  ! and alpha_yDayMean. The routine automatically detects whether data are packed or in LonLat format.
  !
  subroutine getDailyData(filename, grid, run_dynveg, cbalance, dynveg_clim, ntiles, nday)

    USE mo_cbalone_memory, ONLY: cbal_offline_type, grid_offline_type, dynveg_clim_type
    include 'netcdf.inc'
    character(len=*),intent(in)            :: filename
    type(grid_offline_type),intent(in)     :: grid
    LOGICAL,               INTENT(in)      :: run_dynveg
    type(cbal_offline_type),intent(inout)  :: cbalance 
    TYPE(dynveg_clim_type),INTENT(inout)   :: dynveg_clim 
    integer,intent(in)                     :: ntiles
    integer,intent(out)                    :: nday               !! number of days in the considered month
    ! local variables

    logical :: lai_found, npp_found, temp_found, alpha_found
    logical :: prev_year_npp_found, bio_exist_found, rel_hum_found
    logical :: max_wind10_found , prev_day_max_wind10_found

    real,     allocatable :: hlp2D_sp(:,:),hlp3D_sp(:,:,:),hlp4D_sp(:,:,:,:)
    real(dp), allocatable :: hlp2D(:,:),hlp3D(:,:,:),hlp4D(:,:,:,:)

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
    integer,allocatable ,dimension(:)  ::  nvartyp
    integer,allocatable,dimension(:)   ::  nvardim
    integer,allocatable,dimension(:,:) ::  nvardimid
    integer,allocatable,dimension(:)   ::  nvaratt
    character(len=128),allocatable,dimension(:)  ::  fvarnam
    ! indices
    integer  ::  i,ii,itile,day

    ! preparations

    lai_found = .false.
    npp_found = .false.
    temp_found = .false.
    alpha_found = .false.
    prev_year_npp_found = .false.
    bio_exist_found = .false.
    prev_day_max_wind10_found = .false.
    max_wind10_found = .false.
    rel_hum_found = .false.
    
    ! open the input-file

    if(debug) write (*,*) 'getDailyData(): read from: ',filename
    iret = nf_open(filename,nf_nowrite,ncid)
    call check_err(iret,"Error in getDailyData() when opening file "//trim(filename))

    ! check what is in the input-file

    iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
    call check_err(iret,"Error in getDailyData() when inquiring contents of "//trim(filename))
    if (debug) write (*,*) 'getDailyData(): ndim,nvar,nattr,nunlimdim',ndim,nvar,nattr,nunlimdim

    ! get the dimension name and length

    allocate (ndimlen(ndim))
    allocate (fdimnam(ndim))
    do i=1,ndim
       iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
       call check_err(iret,"Error in getDailyData() when inquiring dimension of "//trim(fdimnam(i))//" in "//trim(filename))
       if (debug) write (*,*) 'getDailyData(): i,fdimnam,ndimlen',i,trim(fdimnam(i)),ndimlen(i)
    end do

    ! set ndate (number of output intervals of the monthly file) and noutput_per_day

    nday = 0
    do i=1,ndim
       if (fdimnam(i) == 'time') nday = ndimlen(i)
    end do
    if (debug) write (*,*) 'getDailyData(): nday: ',nday
    if (nday == 0) then
       write(*,*) 'getDailyData(): variable time not found in '//trim(filename)
       stop
    end if

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
       if (debug) write (*,*) 'getDailyData(): i,fvarnam,nvartyp,nvardim,nvardimid,nvaratt',&
            i,trim(fvarnam(i)),nvartyp(i),nvardim(i),nvardimid(i,1),nvaratt(i)
    end do

    ! get data
    
    do i=1,nvar
       select case (trim(fvarnam(i)))
       case('LAI_yDayMean','NPP_yDayMean','alpha_yDayMean','prev_year_npp','bio_exist')
          select case (nvardim(i))
          case(3) !! The array has 3 dimensions => We find packed data
             if (ndimlen(nvardimid(i,1)) /= grid%nland .or. ndimlen(nvardimid(i,2)) /= ntiles) then
                write(*,*) "getDailyData(): Dimension error for variable "//trim(fvarnam(i))
                write(*,*) "getDailyData(): nland=",grid%nland," but ndimlen(nvardimid(i,1)=",ndimlen(nvardimid(i,1))
                write(*,*) "getDailyData(): ntiles=",ntiles," but ndimlen(nvardimid(i,2)=",ndimlen(nvardimid(i,2))
                stop
             end if
             if (ndimlen(nvardimid(i,3)) /= nday) then
                write(*,*) "getDailyData(): Number of time steps wrong for input variable "//trim(fvarnam(i))
                stop
             end if

             allocate (hlp3D(grid%nland,ntiles,nday))
             select case (nvartyp(i))
             case(NF_REAL)
                allocate (hlp3D_sp(grid%nland,ntiles,nday))
                iret = nf_get_var_real(ncid,i,hlp3D_sp)
                call check_err(iret,"Error(1) in getDailyData() reading variable "//trim(fvarnam(i))//" from "//trim(filename))
                hlp3D = REAL(hlp3D_sp)
                deallocate(hlp3D_sp)
             case(NF_DOUBLE)
                allocate (hlp3D(grid%nland,ntiles,nday))
                iret = nf_get_var_real(ncid,i,hlp3D)
                call check_err(iret,"Error(1) in getDailyData() reading variable "//trim(fvarnam(i))//" from "//trim(filename))
             case default
                write(*,*) "ERROR: unsupported type of variable ", trim(fvarnam(i)), " in file ", trim(filename)
                stop
             end select
          case(4)!! The array has 4 dimensions => We find latlon data
             if(ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat &
                  .or. ndimlen(nvardimid(i,3)) /= ntiles) then
                write(*,*) "getDailyData(): Dimension error for variable "//trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getDailyData(): nlon=",grid%nlon," but found nndimlen(nvardimid(i,1))=",ndimlen(nvardimid(i,1))
                write(*,*) "getDailyData(): nlat=",grid%nlat," but found nndimlen(nvardimid(i,2))=",ndimlen(nvardimid(i,2))
                write(*,*) "getDailyData(): ntiles=",ntiles," but ndimlen(nvardimid(i,3)=",ndimlen(nvardimid(i,3))
                stop
             end if
             if (ndimlen(nvardimid(i,4)) /= nday) then
                write(*,*) "getDailyData(): Number of time steps wrong for input variable "//&
                     trim(fvarnam(i))//" in file "//trim(filename)//":"
                stop
             end if
             allocate (hlp4D(grid%nlon,grid%nlat,ntiles,nday))
             select case (nvartyp(i))
             case(NF_REAL)
                allocate (hlp4D_sp(grid%nlon,grid%nlat,ntiles,nday))
                iret = nf_get_var_real(ncid,i,hlp4D_sp)
                call check_err(iret,"Error(2)in getDailyData() reading variable "//trim(fvarnam(i))//" from "//trim(filename))
                hlp4D=real(hlp4d_sp)
                deallocate(hlp4D_sp)
             case(NF_DOUBLE)
                iret = nf_get_var_real(ncid,i,hlp4D)
                call check_err(iret,"Error(2)in getDailyData() reading variable "//trim(fvarnam(i))//" from "//trim(filename))
             case default
                write(*,*) "ERROR: unsupported type of variable ", trim(fvarnam(i)), " in file ", trim(filename)
                stop
             end select

             ! pack the array
             allocate (hlp3D(grid%nland,ntiles,nday))
             do itile=1,ntiles
                do day=1,nday
                   hlp3D(:,itile,day) = pack(hlp4D(:,:,itile,day),mask=grid%mask)
                end do
             end do
             deallocate(hlp4D)
          case default
             write(*,*) "getDailyData(): ERROR: Variable "//trim(fvarnam(i))//&
                          " from file "//trim(filename)//" is neither of type 'packed' nor 'LonLat'."
             stop
          end select
       case('topSoilTemp_yDayMean','rel_hum_air','prev_day_max_wind10','max_wind10')
          select case (nvardim(i))
          case(2) !! The array has 2 dimensions => We find packed data
             if (ndimlen(nvardimid(i,1)) /= grid%nland) then
                write(*,*) "getDailyData(): Dimension error for variable "//trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getDailyData(): nland=",grid%nland," but ndimlen(nvardimid(i,1)=",ndimlen(nvardimid(i,1))
                stop
             end if
             if (ndimlen(nvardimid(i,2)) /= nday) then
                write(*,*) "getDailyData(): Number of time steps wrong for input variable "//trim(fvarnam(i))//&
                     " in file "//trim(filename)//":"
                write(*,*) "getDailyData(): nday=",nday," but ndimlen(nvardimid(i,2))=",ndimlen(nvardimid(i,2))
                stop
             end if
             allocate (hlp2D(grid%nland,nday))
             select case (nvartyp(i))
             case(NF_REAL)
               allocate (hlp2D_sp(grid%nland,nday))
               iret = nf_get_var_real(ncid,i,hlp2D_sp)
               call check_err(iret,"Error(3) in getDailyData() reading variable "//trim(fvarnam(i))//" from "//trim(filename))
               hlp2D=real(hlp2d_sp)
               deallocate(hlp2D_sp)
             case(NF_DOUBLE)
               iret = nf_get_var_real(ncid,i,hlp2D)
               call check_err(iret,"Error(3) in getDailyData() reading variable "//trim(fvarnam(i))//" from "//trim(filename))
             case default
                write(*,*) "ERROR: unsupported type of variable ", trim(fvarnam(i)), " in file ", trim(filename)
                stop
             end select

          case(3)!! The array has 3 dimensions => We find latlon data
             if(ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                write(*,*) "getDailyData(): Dimension error for variable "//trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getDailyData(): nlon=",grid%nlon," but found nndimlen(nvardimid(i,1))=",ndimlen(nvardimid(i,1))
                write(*,*) "getDailyData(): nlat=",grid%nlat," but found nndimlen(nvardimid(i,2))=",ndimlen(nvardimid(i,2))
                stop
             end if
             if (ndimlen(nvardimid(i,3)) /= nday) then
                write(*,*) "getDailyData(): Number of time steps wrong for input variable "//&
                     trim(fvarnam(i))//" in file "//trim(filename)//":"
                write(*,*) "getDailyData(): nday=",nday," but ndimlen(nvardimid(i,3))=",ndimlen(nvardimid(i,3))
                stop
             end if
             allocate (hlp3D(grid%nlon,grid%nlat,nday))
             select case (nvartyp(i))
             case(NF_REAL)
                allocate (hlp3D_sp(grid%nlon,grid%nlat,nday))
                iret = nf_get_var_real(ncid,i,hlp3D_sp)
                call check_err(iret,"Error(4) in getDailyData() when reading variable "//trim(fvarnam(i))//" from "//trim(filename))
                hlp3D=real(hlp3D_sp)
                deallocate(hlp3D_sp)
             case(NF_DOUBLE)
                iret = nf_get_var_real(ncid,i,hlp3D)
                call check_err(iret,"Error(4) in getDailyData() when reading variable "//trim(fvarnam(i))//" from "//trim(filename))
             case default
                write(*,*) "ERROR: unsupported type of variable ", trim(fvarnam(i)), " in file ", trim(filename)
                stop
             end select

             ! pack the array
             allocate (hlp2D(grid%nland,nday))
             do day=1,nday
                hlp2D(:,day) = pack(hlp3D(:,:,day),mask=grid%mask)
             end do
             deallocate(hlp3D)
          case default
             write(*,*) "getDailyData(): ERROR: Variable "//trim(fvarnam(i))//&
                  " from file "//trim(filename)//" is neither of type 'packed' nor 'LonLat'."
             stop
          end select
       end select

       !! copy data to correct variable

       select case(fvarnam(i))
       case('LAI_yDayMean')
          cbalance%LAI_yDayMean(:,:,1:nday) = hlp3D(:,:,1:nday)
          deallocate(hlp3D)
          lai_found=.true.
       case('NPP_yDayMean')
          cbalance%NPP_yDayMean(:,:,1:nday) = hlp3D(:,:,1:nday)
          deallocate(hlp3D)
          npp_found=.true.
       case('topSoilTemp_yDayMean')
          do itile=1,ntiles
             cbalance%topSoilTemp_yDayMean(:,itile,1:nday) = hlp2D(:,1:nday)
          end do
          deallocate(hlp2D)
          temp_found=.true.
       case('alpha_yDayMean')
          cbalance%alpha_yDayMean(:,:,1:nday) = hlp3D(:,:,1:nday)
          deallocate(hlp3D)
          alpha_found=.true.
       case('prev_year_npp')
          IF (run_dynveg) dynveg_clim%prev_year_npp(:,:,1:nday) = hlp3D(:,:,1:nday)
          deallocate(hlp3D)
          prev_year_npp_found=.true.
       case('bio_exist')
          IF (run_dynveg) dynveg_clim%bio_exist(:,:,1:nday) = hlp3D(:,:,1:nday)
          deallocate(hlp3D)
          bio_exist_found=.true.
       case('rel_hum_air')
          IF (run_dynveg) dynveg_clim%rel_hum_air(:,1:nday) = hlp2D(:,1:nday)
          deallocate(hlp2D)
          rel_hum_found=.true.
       case('prev_day_max_wind10')
          IF (run_dynveg) dynveg_clim%prev_day_max_wind10(:,1:nday) = hlp2D(:,1:nday)
          deallocate(hlp2D)
          prev_day_max_wind10_found=.true.
       case('max_wind10')
          IF (run_dynveg) dynveg_clim%max_wind10(:,1:nday) = hlp2D(:,1:nday)
          deallocate(hlp2D)
          max_wind10_found=.true.
       end select

    end do

    ! close the file

    iret = nf_close(NCID)
    call check_err(iret,"Error(4) in getDailyData() when closing "//trim(filename))

    ! check whether all data have been found

    if(.not. LAI_found) then
       write(*,*) "ERROR: Did not find variable LAI_yDayMean in "//trim(filename)//"."
       stop
    end if
    if(.not. NPP_found) then
       write(*,*) "ERROR: Did not find variable NPP_yDayMean in "//trim(filename)//"."
       stop
    end if
    if(.not. temp_found) then
       write(*,*) "ERROR: Did not find variable topSoilTemp_yDayMean in "//trim(filename)//"."
       stop
    end if
    if(.not. alpha_found) then
       write(*,*) "ERROR: Did not find variable alpha_yDayMean in "//trim(filename)//"."
       stop
    end if
    IF (run_dynveg) THEN
       if(.not. prev_year_npp_found) then
          write(*,*) "ERROR: Did not find variable prev_year_npp in "//trim(filename)//"."
          stop
       end if
       if(.not. bio_exist_found) then
          write(*,*) "ERROR: Did not find variable bio_exist in "//trim(filename)//"."
          stop
       end if
       if(.not. rel_hum_found) then
          write(*,*) "ERROR: Did not find variable rel_hum_air in "//trim(filename)//"."
          stop
       end if
       if(.not. prev_day_max_wind10_found) then
          write(*,*) "ERROR: Did not find variable prev_day_max_wind10 in "//trim(filename)//"."
          stop
       end if
       if(.not. max_wind10_found) then
          write(*,*) "ERROR: Did not find variable max_wind10 in "//trim(filename)//"."
          stop
       end if
    ENDIF

    ! deallocate local arrays

    deallocate(fdimnam,ndimlen)
    deallocate(fvarnam,nvartyp,nvardim,nvardimid,nvaratt)
  end subroutine getDailyData

  ! --- initNetcdfInfo() -------------------------------------------------------------------------------------------

  subroutine initNetcdfInfo(cbalance, LC_change, run_dynveg, dynveg_feedback, dynveg)

    USE mo_cbalone_memory, ONLY: cbal_offline_type,landcover_change_type, dynveg_offline_type
    type(cbal_offline_type)    ,intent(in) :: cbalance
    type(landcover_change_type),intent(in) :: LC_change
    LOGICAL                    ,INTENT(in) :: run_dynveg
    LOGICAL                    ,INTENT(in) :: dynveg_feedback
    TYPE(dynveg_offline_type)  ,INTENT(in) :: dynveg

    include 'netcdf.inc'

    integer :: n

    numberOfOutputVariables = 21
    IF (LC_change%do_landcover_change) numberOfOutputVariables = numberOfOutputVariables + 4
    IF (run_dynveg)                    numberOfOutputVariables = numberOfOutputVariables + 8
    IF (run_dynveg .AND. dynveg_feedback) numberOfOutputVariables = numberOfOutputVariables + 5
    allocate(netcdfInfo(1:numberOfOutputVariables))
    n=0

    n=n+1
    netcdfInfo(n)%varName = "Cpool_green"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles => cbalance%Cpool_green

    n=n+1
    netcdfInfo(n)%varName = "Cpool_reserve"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles => cbalance%Cpool_reserve

    n=n+1
    netcdfInfo(n)%varName = "Cpool_woods"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles => cbalance%Cpool_woods

    n=n+1
    netcdfInfo(n)%varName = "Cpool_litter_leaf"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles => cbalance%Cpool_litter_leaf

    n=n+1
    netcdfInfo(n)%varName = "Cpool_litter_wood"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles => cbalance%Cpool_litter_wood

    n=n+1
    netcdfInfo(n)%varName = "Cpool_fast"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles => cbalance%Cpool_fast

    n=n+1
    netcdfInfo(n)%varName = "Cpool_slow"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_DOUBLE
    netcdfInfo(n)%array_tiles => cbalance%Cpool_slow

    n=n+1
    netcdfInfo(n)%varName = "box_Cpools_total"
    netcdfInfo(n)%dimType = BOX_TYPE
    netcdfInfo(n)%varType = NF_FLOAT
    netcdfInfo(n)%array_box => cbalance%box_Cpools_total

    n=n+1
    netcdfInfo(n)%varName = "box_NEP_wholeRun"
    netcdfInfo(n)%dimType = BOX_TYPE
    netcdfInfo(n)%varType = NF_FLOAT
    netcdfInfo(n)%array_box => cbalance%box_NEP_wholeRun

    n=n+1
    netcdfInfo(n)%varName = "avg_Cpool_green"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_FLOAT
    netcdfInfo(n)%array_tiles => cbalance%avg_Cpool_green

    n=n+1
    netcdfInfo(n)%varName = "avg_Cpool_reserve"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_FLOAT 
    netcdfInfo(n)%array_tiles => cbalance%avg_Cpool_reserve

    n=n+1
    netcdfInfo(n)%varName = "avg_Cpool_woods"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_FLOAT
    netcdfInfo(n)%array_tiles => cbalance%avg_Cpool_woods 

    n=n+1
    netcdfInfo(n)%varName = "avg_Cpool_litter_leaf"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_FLOAT
    netcdfInfo(n)%array_tiles => cbalance%avg_Cpool_litter_leaf

    n=n+1
    netcdfInfo(n)%varName = "avg_Cpool_litter_wood"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_FLOAT
    netcdfInfo(n)%array_tiles => cbalance%avg_Cpool_litter_wood 

    n=n+1
    netcdfInfo(n)%varName = "avg_Cpool_fast"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_FLOAT
    netcdfInfo(n)%array_tiles => cbalance%avg_Cpool_fast

    n=n+1
    netcdfInfo(n)%varName = "avg_Cpool_slow"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_FLOAT
    netcdfInfo(n)%array_tiles => cbalance%avg_Cpool_slow

    n=n+1
    netcdfInfo(n)%varName = "avg_soil_respiration"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_FLOAT
    netcdfInfo(n)%array_tiles =>  cbalance%avg_soil_respiration

    n=n+1
    netcdfInfo(n)%varName = "avg_NPP_yDayMean"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_FLOAT
    netcdfInfo(n)%array_tiles =>  cbalance%avg_NPP_yDayMean

    n=n+1
    netcdfInfo(n)%varName = "avg_NPP_flux_correction"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_FLOAT
    netcdfInfo(n)%array_tiles =>  cbalance%avg_NPP_flux_correction

    n=n+1
    netcdfInfo(n)%varName = "avg_litter_flux"
    netcdfInfo(n)%dimType = TILES_TYPE
    netcdfInfo(n)%varType = NF_FLOAT
    netcdfInfo(n)%array_tiles =>  cbalance%avg_litter_flux

    n=n+1
    netcdfInfo(n)%varName = "avg_box_NEP"
    netcdfInfo(n)%dimType = BOX_TYPE
    netcdfInfo(n)%varType = NF_FLOAT
    netcdfInfo(n)%array_box =>  cbalance%avg_box_NEP

    if(LC_change%do_landcover_change) then

       n=n+1
       netcdfInfo(n)%varName = "LCC_sum_box_C2atmos"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_FLOAT
       netcdfInfo(n)%array_box =>  LC_change%LCC_sum_box_C2atmos

       n=n+1
       netcdfInfo(n)%varName = "LCC_sum_box_C2fastSoilPool"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_FLOAT
       netcdfInfo(n)%array_box =>  LC_change%LCC_sum_box_C2fastSoilPool

       n=n+1
       netcdfInfo(n)%varName = "LCC_sum_box_C2slowSoilPool"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_FLOAT
       netcdfInfo(n)%array_box =>  LC_change%LCC_sum_box_C2slowSoilPool

       n=n+1
       netcdfInfo(n)%varName = "LCC_emissions_wholeRun"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_FLOAT
       netcdfInfo(n)%array_box =>  LC_change%LCC_emissions_wholeRun

    end if

    IF (run_dynveg) THEN
       
       n=n+1
       netcdfInfo(n)%varName = "act_fpc"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_FLOAT
       netcdfInfo(n)%array_tiles =>  dynveg%act_fpc

       n=n+1
       netcdfInfo(n)%varName = "pot_fpc"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_FLOAT
       netcdfInfo(n)%array_tiles =>  dynveg%pot_fpc

       n=n+1
       netcdfInfo(n)%varName = "burned_fpc"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_FLOAT
       netcdfInfo(n)%array_tiles =>  dynveg%burned_fpc

       n=n+1
       netcdfInfo(n)%varName = "damaged_fpc"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_FLOAT
       netcdfInfo(n)%array_tiles =>  dynveg%damaged_fpc

       n=n+1
       netcdfInfo(n)%varName = "bare_fpc"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_FLOAT
       netcdfInfo(n)%array_box =>  dynveg%bare_fpc

       n=n+1
       netcdfInfo(n)%varName = "desert_fpc"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_FLOAT
       netcdfInfo(n)%array_box =>  dynveg%desert_fpc

       n=n+1
       netcdfInfo(n)%varName = "max_green_bio"
       netcdfInfo(n)%dimType = TILES_TYPE
       netcdfInfo(n)%varType = NF_FLOAT
       netcdfInfo(n)%array_tiles =>  dynveg%max_green_bio

       n=n+1
       netcdfInfo(n)%varName = "sum_green_bio_memory"
       netcdfInfo(n)%dimType = BOX_TYPE
       netcdfInfo(n)%varType = NF_FLOAT
       netcdfInfo(n)%array_box =>  dynveg%sum_green_bio_memory

       IF (dynveg_feedback) THEN
          n=n+1
          netcdfInfo(n)%varName = "carbon_2_LeafLitterPool"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_FLOAT
          netcdfInfo(n)%array_box => dynveg%carbon_2_LeafLitterPool

          n=n+1
          netcdfInfo(n)%varName = "carbon_2_WoodLitterPool"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_FLOAT
          netcdfInfo(n)%array_box => dynveg%carbon_2_WoodLitterPool

          n=n+1
          netcdfInfo(n)%varName = "carbon_2_slowSoilPool_damage"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_FLOAT
          netcdfInfo(n)%array_box => dynveg%carbon_2_slowSoilPool_damage

          n=n+1
          netcdfInfo(n)%varName = "carbon_2_slowSoilPool_fire"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_FLOAT
          netcdfInfo(n)%array_box => dynveg%carbon_2_slowSoilPool_fire

          n=n+1
          netcdfInfo(n)%varName = "carbon_2_atmos"
          netcdfInfo(n)%dimType = BOX_TYPE
          netcdfInfo(n)%varType = NF_FLOAT
          netcdfInfo(n)%array_box => dynveg%carbon_2_atmos
       END IF

    END IF

    if(n /= numberOfOutputVariables) stop 'initNetcdfInfo(): programming error!!'

  end subroutine initNetcdfInfo

  ! --- writeSingleTimeStep() --------------------------------------------------------------------------------------------

  subroutine writeSingleTimeStep(filename, grid, cbalance, run_dynveg, dynveg_feedback, dynveg, ntiles, LC_change,&
                                 timeStep, finish)

    USE mo_cbalone_memory, ONLY: grid_offline_type, cbal_offline_type, vegetation_offline_type, landcover_change_type, &
                                 dynveg_offline_type

    include 'netcdf.inc'
    character(len=*),intent(in)            :: filename
    type(grid_offline_type),intent(in)     :: grid
    type(cbal_offline_type),intent(in)     :: cbalance
    LOGICAL, INTENT(in)                    :: run_dynveg
    LOGICAL, INTENT(in)                    :: dynveg_feedback
    TYPE(dynveg_offline_type), INTENT(in)  :: dynveg
    integer,intent(in)                     :: ntiles
    type(landcover_change_type),intent(in) :: LC_change
    real(dp),intent(in)                    :: timeStep !! indicator for time (in whatever units you like)
    logical,intent(in),optional            :: finish   !! If "true": routine closes the output file, deallocates memory and exits.

    !! state variables

    logical,save :: initialized = .false.
    integer,save :: fileId   ! netCDF file ID
    integer,save :: box_type_DimIds(1:3)
    integer,save :: tiles_type_DimIds(1:4)
    integer,save :: dim,dimIds(1:4)
    integer,save :: varId_lat,varId_lon,varId_tiles,varId_time
    real(dp),pointer,save ::  array_4D(:,:,:,:)   !! dummy array: lat x lon x tiles x 1 (1 time step)
    real(dp),pointer,save ::  array_3D(:,:,:)     !! dummy array: lat x lon x 1 (1 time step)
    real    ,pointer,save ::  array_sp(:,:,:)     !! dummy array in single precission
    integer,save :: timeStepCounter


    !! local variables

    integer :: iret                                                   ! error status return
    integer :: dimid_lat(1),dimid_lon(1),dimid_tiles(1),dimid_time(1) ! netCDF dimension IDs
    integer :: ivar,itile,dimId
    integer :: tilesArray(1:ntiles)
    real(dp):: timeArray(1:1)
    integer :: start(1:4),count(1:4)


    !! Check for finishing

    if(present(finish)) then
        if(initialized) then 

          !! close the output-file
          iret = nf_close(fileId)
          call check_err(iret)

          !! deallocate memory
          deallocate(array_4D,array_3D)

          !! turn into non-initialized status
          initialized = .false.

          !! exit
          return

       end if
    end if

    !! Initialization

    if(.not. initialized) then

       timeStepCounter = 0

       !! initialize info structure for netcdf output fields

       call initNetcdfInfo(cbalance, LC_change, run_dynveg, dynveg_feedback, dynveg)

       ! control output
       if (debug) write(*,*) "writeSingleTimeStep(): open output file ",trim(filename)

       !! allocate memory for future use
       allocate (array_4D(grid%nlon,grid%nlat,ntiles,1))
       allocate (array_3D(grid%nlon,grid%nlat,1))

       ! open the output-file
       if (debug) write (*,*) 'writeSingleTimeStep(): write Cpools: ',filename
       iret = nf_create(filename,nf_clobber,fileId)
       call check_err(iret)
       ! define dimensions
       iret = nf_def_dim(fileId,'lon',grid%nlon,dimId)
       call check_err(iret)
       dimid_lon(1) = dimId
       tiles_type_DimIds(1) = dimId
       box_type_DimIds(1)   = dimId

       iret = nf_def_dim(fileId,'lat',grid%nlat,dimId)
       call check_err(iret)
       dimid_lat(1) = dimId
       tiles_type_DimIds(2) = dimId
       box_type_DimIds(2)   = dimid

       iret = nf_def_dim(fileId,'tiles',ntiles,dimId)
       call check_err(iret)
       dimid_tiles(1) = dimId
       tiles_type_DimIds(3) = dimId

       iret = nf_def_dim(fileId,'time',NF_UNLIMITED,dimId)
       call check_err(iret)
       dimid_time(1) = dimId
       tiles_type_DimIds(4) = dimId
       box_type_DimIds(3)   = dimId

       !! set dimensions for arrays

       ! define longitudes, latitudes, tiles and time as variables
       iret = nf_def_var(fileId,'lat',NF_FLOAT,1,dimid_lat(1:1),varId_lat)
       call check_err(iret)
       iret = nf_def_var(fileId,'lon',NF_FLOAT,1,dimid_lon(1:1),varId_lon)
       call check_err(iret)
       iret = nf_def_var(fileId,'tiles',NF_INT,1,dimid_tiles(1:1),varId_tiles)
       call check_err(iret)
       iret = nf_def_var(fileId,'time',NF_DOUBLE,1,dimid_time(1:1),varId_time)
       call check_err(iret)

       ! define all other output variables

       do ivar =1,numberOfOutputVariables
          select case(netcdfInfo(ivar)%dimType)
          case(BOX_TYPE)
             dim=3
             dimIds(1:3) = box_type_DimIds(1:3)
          case(TILES_TYPE)
             dim=4
             dimIds(1:4) = tiles_type_DimIds(1:4)
          case default
             stop "writeSingleTimeStep(): non-existent dimension type"
          end select
          iret = nf_def_var(fileId,                                &
                            trim(netcdfInfo(ivar)%varName),        &
                            netcdfInfo(ivar)%varType,              &
                            dim,                                   &
                            dimIds(1:dim),                         &
                            netcdfInfo(ivar)%varId                 )
          call check_err(iret)
          IF (netcdfInfo(ivar)%varType == nf_double) THEN
             iret = nf_put_att_double(fileId, netcdfInfo(ivar)%varId, '_FillValue', &
                 nf_double,1,nf_fill_double)
          END IF
          call check_err(iret)
       end do

       !! end definition mode

       iret = nf_enddef(FILEID)
       call check_err(IRET)

       !! put variables (except time) that derive from dimensions

       iret = nf_put_var_double(fileId,varid_lat,grid%lat(1:grid%nlat))
       call check_err(iret)

       iret = nf_put_var_double(fileId,varid_lon,grid%lon(1:grid%nlon))
       call check_err(iret)

       do itile=1,ntiles
          tilesArray(itile) = itile
       end do
       iret = nf_put_var_int(fileId,varId_tiles,tilesArray(1:ntiles))
       call check_err(iret)

       !! remember initialization status

       initialized = .true.

    end if

    !! increase counter for time steps

    timeStepCounter = timeStepCounter + 1

    !! put data

    timeArray(1)=timeStep
    start(1)=timeStepCounter
    count(1)=1
    iret = nf_put_vara_double(fileId,varId_time,start(1:1),count(1:1),timeArray(1:1))
    call check_err(iret)

    do ivar =1,numberOfOutputVariables
       select case(netcdfInfo(ivar)%dimType)
       case(BOX_TYPE)
          start(1:3)=(/ 1,1,timeStepCounter /)
          count(1:3)=(/ grid%nlon,grid%nlat,1 /)
          if(debug) write(*,*) "Box-case: Try to write variable "//trim(netcdfInfo(ivar)%varName)
          array_3D(:,:,1) = UNPACK(netcdfInfo(ivar)%array_box(:),mask=grid%mask,field=nf_fill_double)

          select case(netcdfInfo(ivar)%varType)
          case(NF_FLOAT)
             allocate (array_sp(grid%nlon,grid%nlat,1))
             array_sp(:,:,1)=array_3D(1:grid%nlon,1:grid%nlat,1)
             iret = nf_put_vara_real(fileId,netcdfInfo(ivar)%varId,start(1:3),count(1:3),array_sp(:,:,:))
             call check_err(iret)
             deallocate (array_sp)
          case(NF_DOUBLE)
             iret = nf_put_vara_double(fileId,netcdfInfo(ivar)%varId,start(1:3),count(1:3),array_3D(1:grid%nlon,1:grid%nlat,1:1))
             call check_err(iret)
          case default
             stop "Non-existent case lhv�.jv"
          end select

       case(TILES_TYPE)
          start(1:4) = (/ 1,1,1,timeStepCounter /)
          count(1:4) = (/ grid%nlon,grid%nlat,ntiles,1 /)
          if(debug) write(*,*) "Tiles-case: Try to write variable "//trim(netcdfInfo(ivar)%varName)
          do itile = 1,ntiles
             array_4D(:,:,itile,1) = UNPACK(netcdfInfo(ivar)%array_tiles(:,itile),mask=grid%mask,field=nf_fill_double)
          enddo

          select case(netcdfInfo(ivar)%varType)
          case(NF_FLOAT)
             allocate(array_sp(grid%nlon,grid%nlat,ntiles))
             array_sp = array_4D(1:grid%nlon,1:grid%nlat,1:ntiles,1)
             iret = nf_put_vara_real(fileId, netcdfInfo(ivar)%varId, start(1:4), count(1:4), array_sp(:,:,:))
             call check_err(iret,"Tiles-case: Error happened when trying to write "//trim(netcdfInfo(ivar)%varName))
             deallocate(array_sp)
          case(NF_DOUBLE)
             iret = nf_put_vara_double(fileId, netcdfInfo(ivar)%varId, start(1:4), count(1:4), &
                                       array_4D(1:grid%nlon,1:grid%nlat,1:ntiles,1))
             call check_err(iret,"Tiles-case: Error happened when trying to write "//trim(netcdfInfo(ivar)%varName))
          case default
             stop "Non-existent case .kugb-.lqjb"
          end select

       case default
          stop "writeSingleTimeStep(): non-existent dimension type"
       end select
    end do

  end subroutine writeSingleTimeStep


  !! --- read_landcover_fractions() ----------------------------------------------------
  !!
  !! This routine reads in new landcover fractions
  !!
  subroutine read_landcover_fractions(filename,grid,ntiles,cover_fract)
    USE mo_cbalone_memory, ONLY: grid_offline_type
    USE mo_kind,           ONLY: dp
    include 'netcdf.inc'
    character(len=*),        intent(in)    :: filename
    TYPE(grid_offline_type), intent(in)    :: grid
    integer,                 intent(in)    :: ntiles
    real(dp),                intent(out)   :: cover_fract(:,:)

    !! locals

    integer :: iret,ncid,ndim,nvar,nattr,nunlimdim
    integer,allocatable            :: ndimlen(:)
    character(len=128),allocatable :: fdimnam(:),fvarnam(:)
    integer, allocatable           :: nvartyp(:),nvardim(:),nvaratt(:)
    integer,allocatable            :: nvardimid(:,:)
    real(dp),allocatable           :: coverFract_2D(:,:,:)
    integer :: i,ii,itile
    integer ::  nin(100)
    
    !! open the file
    write (*,*) 'read_landcover_fractions(): read vegetation data from: ',trim(filename)
    iret = nf_open(filename,nf_nowrite,ncid)
    call check_err(iret)
    ! check what is in the input-file
    iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
    call check_err(iret)
    if (debug) write (*,*) 'read_landcover_fractions(): ndim,nvar,nattr,nunlimdim',ndim,nvar,nattr,nunlimdim
    ! get the dimension name and length
    allocate (ndimlen(ndim))
    allocate (fdimnam(ndim))
    do i=1,ndim
       iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
       call check_err(iret)
       if (debug) write (*,*) 'read_landcover_fractions(): i,fdimnam,ndimlen',i,trim(fdimnam(i)),ndimlen(i)
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
       if (debug) write (*,*) 'read_landcover_fractions(): i,fvarnam,nvartyp,nvardim,nvardimid,nvaratt',&
            i,trim(fvarnam(i)),nvartyp(i),nvardim(i),nvardimid(i,1),nvaratt(i)
    end do

    ! get cover type data
    do i=1,nvar
       if (fvarnam(i) == 'cover_fract') then

          !! Check field dimensions

          if(ndimlen(nvardimid(i,1)) /= grid%nlon) then
             write(*,*) "read_landcover_fractions(): ERROR: When reading field cover_fract from "//trim(filename)//" :"
             write(*,*) "read_landcover_fractions(): Expected longitudinal dimension ",grid%nlon, &
                                         "but found ",ndimlen(nvardimid(i,1)),"."
             stop
          end if
          if(ndimlen(nvardimid(i,2)) /= grid%nlat) then
             write(*,*) "read_landcover_fractions(): ERROR: When reading field cover_fract from "//trim(filename)//" :"
             write(*,*) "read_landcover_fractions(): Expected latitudinal dimension ",grid%nlat, &
                                         "but found ",ndimlen(nvardimid(i,2)),"."
             stop
          end if
          if(ndimlen(nvardimid(i,3)) /= ntiles) then
             write(*,*) "read_landcover_fractions(): ERROR: When reading field cover_fract from "//trim(filename)//" :"
             write(*,*) "read_landcover_fractions(): Expected tile dimension ",ntiles, &
                                         "but found ",ndimlen(nvardimid(i,3)),"."
          end if
          if (debug) write (*,*) 'read_landcover_fractions(): cover_fract dimensions: ', &
               ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)),ndimlen(nvardimid(i,3))

          !! Allocate temporary memory

          allocate (coverFract_2D(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
               ndimlen(nvardimid(i,3))))
          if (debug) write (*,*) 'read_landcover_fractions(): cover_type dimensions: ', &
               ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)),ndimlen(nvardimid(i,3))
          iret = nf_get_var_double(ncid,i,coverFract_2D)
          call check_err(iret)

          ! pack cover fractions

          do itile=1,ntiles
             cover_fract(:,itile) = PACK(coverFract_2D(:,:,itile),mask=grid%mask)
          end do
          deallocate(coverFract_2D)

       end if
    end do

    ! final deallocations

    deallocate (ndimlen)
    deallocate (fdimnam)
    deallocate (fvarnam)
    deallocate (nvartyp)
    deallocate (nvardim)
    deallocate (nvardimid)
    deallocate (nvaratt)

    !! close the output-file

    iret = nf_close(ncid)
    call check_err(iret,"Error in read_landcover_fractions() when closing "//trim(filename))

  end subroutine read_landcover_fractions

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

end module mo_cbalone_io
