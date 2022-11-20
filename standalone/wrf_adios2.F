!    files and reads input files that are somewhat compatible with wrf files
!

module wrf_adios2
use adios2
use module_fr_sfire_util , only : crash, interpolate_2d, continue_at_boundary
use module_domain, only: domain
use wrf_netcdf, only : grid_info, nf90_max_name, nf90_float



implicit none
include 'mpif.h'

! default file names
character(len=*),parameter::inputfile='fire_input.nc'
character(len=*),parameter::outputfile='fire_output.nc'

! control whether the dimensions of the fire grid variables should have 
! the same (incorrect) sizes that wrf outputs and same dimension names 
! 
logical::compat_fire_grid=.true.
logical::debug_print=.false.
logical::debug_print_f=.true.
logical::read_check=.false.

! output variable type
integer,parameter::vartype=nf90_float,field_type=104

! max number of dimensions
integer, parameter:: mdims=4

! variable names
character(len=nf90_max_name)::                        &
                              var_xtime='XTIME',        &
                              unit_xtime='min',                &
                              desc_xtime='minutes since simulation start',                &
                              var_itimestep='ITIMESTEP',        &
                              unit_itimestep='',                &
                              desc_itimestep='',                &
                              var_nfuel_cat='NFUEL_CAT',        &
                              unit_nfuel_cat='',                &
                              desc_nfuel_cat='',                &
                              var_dzdxf='DZDXF',                &
                              unit_dzdxf='',                    &
                              desc_dzdxf='',                    &
                              var_dzdyf='DZDYF',                &
                              unit_dzdyf='',                    &
                              desc_dzdyf='',                    &
                              var_zsf='ZSF',                    &
                              unit_zsf='',                      &
                              desc_zsf='',                      &
                              var_lfn='LFN',                    &
                              unit_lfn='',                      &
                              desc_lfn='',                      &
                              var_tign_g='TIGN_G',              &
                              unit_tign_g='',                   &
                              desc_tign_g='',                   &
                              var_tign_in='TIGN_IN',            &
                              unit_tign_in='',                  &
                              desc_tign_in='',                  &
                              var_fmc_g='FMC_G',                &
                              unit_fmc_g='',                    &
                              desc_fmc_g='',                    &
                              var_fmc_gc='FMC_GC',              &
                              unit_fmc_gc='',                   &
                              desc_fmc_gc='',                   &
                              var_fxlong='FXLONG',              &
                              unit_fxlong='',                   &
                              desc_fxlong='',                   &
                              var_fxlat='FXLAT',                &
                              unit_fxlat='',                    &
                              desc_fxlat='',                    &
                              var_unit_fxlong='UNIT_FXLONG',    &
                              unit_unit_fxlong='',              &
                              desc_unit_fxlong='',              &
                              var_unit_fxlat='UNIT_FXLAT',      &
                              unit_unit_fxlat='',               &
                              desc_unit_fxlat='',               &
                              var_uf='UF',                      &
                              unit_uf='',                       &
                              desc_uf='',                       &
                              var_vf='VF',                      &
                              unit_vf='',                       &
                              desc_vf='',                       &
                              var_fuel_frac='FUEL_FRAC',        &
                              unit_fuel_frac='',                &
                              desc_fuel_frac='',                &
                              var_fire_area='FIRE_AREA',        &
                              unit_fire_area='',                &
                              desc_fire_area='',                &
                              var_fgrnhfx='FGRNHFX',            &
                              unit_fgrnhfx='',                  &
                              desc_fgrnhfx='',                  &
                              var_fgrnqfx='FGRNQFX',            &
                              unit_fgrnqfx='J/m^2',             &
                              desc_fgrnqfx='heat flux',         &
                              var_ros='ROS',                    &
                              unit_ros='m/s',                   &
                              desc_ros='ROS',                   &
                              var_r_0='R_0',                    &
                              unit_r_0='m/s',                   &
                              desc_r_0='R_0',                   &
                              var_flineint='FLINEINT',          &
                              unit_flineint='J/m/s',            &
                              desc_flineint='Byram fireline intensity',         &
                              var_flineint2='FLINEINT2',        & 
                              unit_flineint2='J/m/s^2',         & 
                              desc_flineint2='New fireline intensity' , &
                              var_f_ros0='F_ROS0',          &
                              unit_f_ros0='m/s',            &
                              desc_f_ros0='base rate of spread in all directions',  &
                              var_f_rosx='F_ROSX',          &
                              unit_f_rosx='m/s',            &
                              desc_f_rosx='potential rate of spread in direction X',  &
                              var_f_rosy='F_ROSY',          &
                              unit_f_rosy='m/s',            &
                              desc_f_rosy='potential rate of spread in direction Y',  &
                              var_f_ros='F_ROS',          &
                              unit_f_ros='m/s',            &
                              desc_f_ros='potential fire max spread rate in any direction',  &
                              var_f_int='F_INT',          &
                              unit_f_int='J/m^2/s',            &
                              desc_f_int='potential fire reaction intensity for risk rating',  &
                              var_f_lineint='F_LINEINT',          &
                              unit_f_lineint='J/m/s',            &
                              desc_f_lineint='potential Byram fireline intensity for risk rating',  &
                              var_f_lineint2='F_LINEINT2',          &
                              unit_f_lineint2='J/m/s^2',            &
                              desc_f_lineint2='potential alternative fireline intensity for risk rating'

character(len=nf90_max_name),parameter::var_times='Times' 

! grid information structure
integer, parameter::max_times_length=19


type(adios2_adios):: adios ! ADIOS
type(adios2_io) :: ioGet ! IO


contains

!
!*** READING 
!

subroutine read_global_info(adios2_id,info, ierr)
  ! get fire grid sizes from input file
  implicit none
  !*** arguments

  type(adios2_engine), intent(in)::adios2_id 
  type(grid_info),intent(inout)::info
  integer, intent(out):: ierr

  !*** local
  type(adios2_attribute):: adios2_attr

  integer::it,ndims,idim,ratio,stagid,i
  integer,dimension(mdims)::dims,sr,dimids,stag
  real::dx,dy
  character(len=128)::msg
  character(len=NF90_MAX_NAME)::dimnames(mdims),dimname,stagname,stagnames(mdims)
  character(len=*), parameter:: subgrid='_subgrid'

  integer:: status = adios2_step_status_not_ready

  !*** executable
  

  ! Feng:
  ! First read all dims information, look at VAR_NAME/Dims

  ! get time string information
  ! We don't have this dim, instead read
  !call read_var_info(filename,var_times,ndims,dims,dimnames)
  !if(ndims.ne.2)call crash(trim(var_times)//' must have 2 dimensions')
  !info%dim_time_string = dimnames(1)
  !info%len_time_string = dims(1)
  !if(dims(1).gt.max_times_length)call crash('time string too long')
  info%dim_time_string = "DateStrLen"
  info%len_time_string = LEN_TRIM(info%dim_time_string)
    
  ! get dimensions
  ! NFUEL_CAT
  

  ! store in info what you can now


  ! attributes
  do while (status == adios2_step_status_not_ready)
    write (*,*) "Stream not ready, sleep 1 sec..."
    call sleep(1)
    call adios2_begin_step(adios2_id, adios2_step_mode_read, 120.0, status, ierr)
  enddo

  if (status /= adios2_step_status_ok) then
    ierr = -1
    call check(ierr)
    return
  endif
  !call adios2_begin_step(adios2_id, ierr)
  ! if (status == adios2_step_status_end_of_stream) then
    
  !   ierr = adios2_step_status_end_of_stream
  !   return
  ! end if
  ! use engine as parmeter, and call begin step ahead of it
  call read_var_info(var_nfuel_cat,ndims,dims,dimnames)

  info%io_nfirex=dims(1)
  info%io_nfirey=dims(2)
  !info%ntimes=dims(3)
  info%ntimes=100
  info%dim_fire_x=dimnames(1)
  info%dim_fire_y=dimnames(2)
  !info%dim_time=dimnames(3)
  info%dim_time="Time"

  call att_read(adios2_id,'DX',info%fdx)
  call att_read(adios2_id,'DY',info%fdy)
  call att_read(adios2_id,'DT',info%dt)

  sr=0
  stag=0
  stagnames=' '
  if(compat_fire_grid)then
    do idim=1,2
      dimname=dimnames(idim)
      i=index(dimname,subgrid,.true.)
      write (*,*) "len of subgrid,", len(subgrid), "len of dimname", len_trim(dimname)
      if ( i + len(subgrid) .eq. len_trim(dimname)+1) then ! subgrid, need to fix
          stagname=dimname(1:i-1) // '_stag'
          ! TO read dim information from attr

          !call check(nf90_get_att(adios2_id,nf90_global,att_name,att_val))
          call adios2_inquire_attribute(adios2_attr, ioGet, '_DIM_'//stagname, ierr)
          call check(ierr)
          call adios2_attribute_data(stag(idim), adios2_attr, ierr)
          call check(ierr)

          !call att_read(adios2_id, '__DIM__'//stagname, stag(idim))
          !call check(nf90_inq_dimid(adios2_id, stagname, stagid))
          !call check(nf90_inquire_dimension(adios2_id, stagid, len=stag(idim)))
          sr(idim) = dims(idim)/stag(idim)
          if(debug_print)write(*,'(3a,i5,a,i5)')'dimension ',trim(stagname), &
            ' length ',stag(idim),' ratio ',sr(idim)
      endif
      stagnames(idim)=stagname
      dims(idim) = dims(idim) - sr(idim)
    enddo
    if((sr(1).ne.0.and.sr(2).eq.0).or.sr(1).lt.0.or.sr(2).lt.0)then
      write(msg,'(a,2i5)')'bad subgrid refinement ratios',sr(1),sr(2)
      call crash(msg)
    endif

  endif

  call adios2_end_step(adios2_id, ierr)

  

  ! store the rest in info
  info%nfirex=dims(1)
  info%nfirey=dims(2)
  info%sr_x=sr(1)
  info%sr_y=sr(2)
  info%nstagx=stag(1)
  info%nstagy=stag(2)
  info%dim_atm_x_s=stagnames(1)
  info%dim_atm_y_s=stagnames(2)
  if(sr(1).ne.0)then
      info%fdx=info%fdx/sr(1)
      info%fdy=info%fdy/sr(2)
  endif
    
  if(debug_print)then
    write(*,'(4(a,1x))')'dimension names:',trim(info%dim_fire_x),trim(info%dim_fire_y), &
      trim(info%dim_time)
    write(*,'(a,2i6)')'fire grid dimensions:',info%nfirex,info%nfirey
    write(*,'(a,i6)')'number of time frames',info%ntimes
    write(*,'(3(a,f8.4,1x))')'stepsizes fdx=',info%fdx,'fdy=',info%fdy,'dt=',info%dt
  endif

end subroutine read_global_info

subroutine att_read(adios2_id,att_name,att_val)
  type(adios2_engine), intent(in)::adios2_id 
  character(len=*),intent(in)::att_name
  real, intent(out)::att_val

  type(adios2_attribute) :: adios2_attr

  integer::ierr
  !call check(nf90_get_att(adios2_id,nf90_global,att_name,att_val))
  call adios2_inquire_attribute(adios2_attr, ioGet, att_name, ierr)
  call check(ierr)
  call adios2_attribute_data(att_val, adios2_attr, ierr)
  call check(ierr)
  if(debug_print)write(*,'(a,1x,a,g20.5)')'attribute',trim(att_name),att_val
end subroutine att_read

!
!***
!

subroutine print_var_info(filename,varname)

  !*** arguments
  character(len=*), intent(in)::filename,varname

  !*** local
  integer::adios2_id,ndims,dimlengths(mdims),dimids(mdims),i,type
  character(len=NF90_MAX_NAME)::dimnames(mdims)

  !*** executable
  if(debug_print)write(*,'(4a)')'reading from adios2 io',' dimensions of variable ',trim(varname)
  call read_var_info(varname,ndims,dimlengths,dimnames)
  if(debug_print)write(*,'(3a,4(2a,i5,1x))')'variable ',trim(varname),' dimensions ', &
    (trim(dimnames(i)),'=',dimlengths(i),i=1,ndims)

end subroutine print_var_info




subroutine read_var_info(varname,ndims,dimlengths,dimnames,type)

  ! get variable dimensions from a file
  implicit none

  ! arguments
  character(len=*), intent(in):: varname ! variable name
  integer, intent(out)::ndims            ! number of dimensions of this variable
  integer, intent(out)::dimlengths(mdims)      ! the dimensions
  character(len=NF90_MAX_NAME),intent(out)::dimnames(mdims) ! dimension names
  integer, intent(out), optional :: type

  ! local
 
  type(adios2_variable) :: adios2_var
  type(adios2_attribute) :: adios2_attr
  integer::ierr

  integer::varid,dimid,idimids(mdims),idim,i,dimlen,xtype
  integer ,allocatable :: Buffer(:)
  character(len=NF90_MAX_NAME)::dimname
  character(len=NF90_MAX_NAME)::Buffer_dimnames(mdims)

  allocate(Buffer(1), STAT=ierr)

  ! executable
  if(debug_print)write(*,'(4a)')'reading from ADIOS2',' variable ',trim(varname)
  !call adios2_open(adios2_id, ioGet, filename, adios2_mode_read, ierr )
  !call check(ierr)

  !call check(nf90_inq_varid(adios2_id, varname, varid))
  !call check(nf90_inquire_variable(adios2_id,varid,ndims=ndims,xtype=xtype))
  call adios2_inquire_variable(adios2_var, ioGet, varname, ierr)
  call check(ierr)
  !if (adios2_var % valid == .false.)
   ! call crash('variable not found')

  write (*,*) "Query dims!"
  call adios2_variable_ndims(ndims, adios2_var, ierr)
  call adios2_variable_type(xtype, adios2_var, ierr)
  write(*, *) "ndims: ", ndims, "xtype", xtype




  if(ndims.gt.mdims)then
      write(*,1)'variable ',trim(varname),' has ',ndims,' dimensions >',mdims
      call crash('variable has too many dimensions')
  1 format(3a,i5,a,(4i6))
  endif

  call adios2_inquire_attribute(adios2_attr, ioGet, trim(varname)//'/Dims', ierr)
  call check(ierr)
  call adios2_attribute_data(Buffer_dimnames, adios2_attr, ierr)
  call check(ierr)
  
  !call check(nf90_inquire_variable(adios2_id, varid, dimids=idimids))
  do idim=1,ndims
      dimnames(idim) = ctrim(Buffer_dimnames(idim))
      write (*,*) LEN(Buffer_dimnames(idim)), LEN(dimnames(idim))
      dimname = dimnames(idim)

      if(debug_print)write(*,*) 'inquiring dimension ',idim,' name', TRIM(dimname), ', dim_name length=', LEN_TRIM(dimname)
      
      call adios2_inquire_attribute(adios2_attr, ioGet, '_DIM_'//dimnames(idim), ierr)
      call check(ierr, msg="Inquire Attribute: "//'_DIM_'//dimnames(idim))
      
      call adios2_attribute_data(Buffer, adios2_attr, ierr)
      call check(ierr)
      dimlen = Buffer(1)

      if(debug_print)write(*,*) 'got dimension length',dimlen
      dimlengths(idim) = dimlen  
  enddo
  !call adios2_close(adios2_id, ierr)
  if(present(type)) type=xtype
  2 format(3a,i4,a,(4i6))
  if(debug_print)write(*,2)'variable ',trim(varname), &
    ' type',xtype,' dimensions ',(dimlengths(idim),idim=1,ndims)
  
    deallocate(Buffer, STAT=ierr)
end subroutine read_var_info

! https://community.intel.com/t5/Intel-Fortran-Compiler/How-to-remove-NUL-character-from-end-of-variable-in-output-file/td-p/1179485
function ctrim(instring) result(outstring)
  implicit none
  character(len=*), intent(in) :: instring
  character(len=len(instring)) :: outstring

  integer :: pos

  pos = index( instring, achar(0) )
  if ( pos > 0 ) then
      outstring = instring(1:pos-1)
  else 
      outstring = instring
  endif
end function ctrim

!
!***
!

! get a real value (real number)
subroutine read_real(adios2_id,iframe,varname,v)
  implicit none
  !*** arguments
  type(adios2_engine), intent(in)::adios2_id             ! id of netcdf file open in data mode
  integer, intent(inout):: iframe          ! number of frame in the file
  character(len=*),intent(in)::varname  ! the variable name
  real :: v                ! value
  !*** local
  integer::varid
  real:: val(1)

  
  type(adios2_variable) :: adios_var
  integer(kind=8), dimension(1) :: istart, icount
  integer:: ierr

  !*** executable

  if(debug_print)write(*,'(3a,i5)')'reading real ',trim(varname),' timestep ',iframe
  
  
 
  !call check(nf90_inq_varid(adios2_id,varname,varid),'cannot find '//trim(varname))
  !call check(nf90_get_var(adios2_id,varid,val,start=(/iframe/),count=(/1/)), &
  !    'error reading '//trim(varname))
  istart = (/0/)
  icount = (/1/)
  call adios2_inquire_variable(adios_var, ioGet, varname, ierr)
  call adios2_set_selection(adios_var, 1, istart, icount, ierr)
  call adios2_get(adios2_id, adios_var, val, ierr)

  v=val(1)
  if(debug_print)write(*,'(2a,i5,a,g18.8)')trim(varname),'(',iframe,')=',v

end subroutine read_real

!
!***
!

subroutine read_fire_var(adios2_id,info,iframe,varname,v)
  implicit none
  !*** arguments
  type(adios2_engine), intent(in)::adios2_id              ! id of netcdf file open in data mode
  type(grid_info),intent(in)::info      ! dimensions
  integer, intent(inout):: iframe          ! number of frame in the file
  character(len=*),intent(in)::varname  ! the variable name
  real, pointer :: v(:,:)            ! values
  !*** local
  integer::varid,nx,ny,start(3),count(3)
  type(adios2_variable) :: adios_var
  integer(kind=8), dimension(2) :: istart, icount
  integer:: ierr, i, j

  real, pointer:: v2(:,:)

  !*** executable
  nx=info%nfirex
  ny=info%nfirey
  
  if(debug_print)write(*,'(2a)')'reading variable ',trim(varname)
  if(debug_print)write(*,'(a,4i10)')'lower bounds',lbound(v)
  if(debug_print)write(*,'(a,4i10)')'upper bounds',ubound(v)
  !call check(nf90_inq_varid(adios2_id,varname,varid),'cannot find '//trim(varname))
  !call check(nf90_get_var(adios2_id,varid,v(1:nx,1:ny),start=(/1,1,iframe/),count=(/nx,ny,1/)), &
  !    'error reading '//trim(varname))
  !if(debug_print)write(*,'(a,2i5,a,1x,e15.5,3(1x,a,e15.5))')'dimensions',nx,ny, &
  !     'min',minval(v(1:nx,1:ny)),'max',maxval(v(1:nx,1:ny)),'(1,1)=',v(1,1),'end=',v(nx,ny)

  istart= (/0,0/)
  icount = (/nx, ny/)

  allocate(v2(nx, ny))
  
  call adios2_inquire_variable(adios_var, ioGet, varname, ierr)
  call check(ierr)
  call adios2_set_selection(adios_var, 2, istart ,icount, ierr)
  call check(ierr)
  call adios2_get(adios2_id, adios_var, v2, adios2_mode_sync,ierr)

  do j=1,ny
    do i=1,nx
      v(i,j)=v2(i,j)
    enddo
  enddo
  call check(ierr)
  
  end subroutine read_fire_var

subroutine read_step_vars(adios2_id,info, iframe, grid, ierr)
  ! read all variables from input file
  ! compared with the netcdf version, this doesn't not have a iframe (cannot specify which step to read)
  ! but calling this function once will advance a step (iframe will +=1)
  ! when ierr != 0, 
  implicit none

  !*** arguments
  type(adios2_engine), intent(in) :: adios2_id ! engine
  type(grid_info),intent(inout)::info   ! dimensions
  integer, intent(inout):: iframe          ! number of frame in the file (not used by ADIOS2)
  type(domain),intent(inout)::grid      ! the mother of all arrays
  integer, intent(out)::ierr ! need to check whether it's the end

  !*** local
  integer::varid
  integer,dimension(4)::s,c

  type(adios2_variable) :: adios_var
  type(adios2_attribute) :: adios2_attr

  integer:: tmp, ierr2, status

  !*** executable
  if(debug_print_f)write(*,*)'read_step_vars: reading from ADIOS2' 
  !call adios2_open(adios2_id, ioGet, filename, adios2_mode_read, ierr)

  !call check(nf90_inq_varid(adios2_id,var_times,varid),'cannot find '//trim(var_times))
  !call check(nf90_get_var(adios2_id,varid,info%times,start=(/1,iframe/),count=(/info%len_time_string,1/)), &
  !    'error reading '//trim(var_times))
#if 0
  call adios2_inquire_attribute(adios2_attr, ioGet, var_times , ierr)
  call check(ierr)
  call adios2_attribute_data(info%times, adios2_attr, ierr)
  call check(ierr)
#else
  info%times="TimeADIOS2Steps"

  ! TODO: I need to print this time
  if(debug_print_f)write(*,'(2a)')'Time ',info%times
  write (*, *) "## Now begin Step", iframe
#endif
  call adios2_begin_step(adios2_id, adios2_step_mode_read, 120.0, status, ierr)
  !call adios2_begin_step(adios2_id, ierr)
  if (status == adios2_step_status_end_of_stream) then
    
    ierr = adios2_step_status_end_of_stream
    return
  end if
  
  call read_real(adios2_id,iframe,var_unit_fxlong,grid%unit_fxlong)
  call read_real(adios2_id,iframe,var_unit_fxlat,grid%unit_fxlat)

  call read_fire_var(adios2_id,info,iframe,var_nfuel_cat,grid%nfuel_cat) 
  call read_fire_var(adios2_id,info,iframe,var_dzdxf,grid%dzdxf) 
  call read_fire_var(adios2_id,info,iframe,var_dzdyf,grid%dzdyf) 
  call read_fire_var(adios2_id,info,iframe,var_zsf,grid%zsf) 
  call read_fire_var(adios2_id,info,iframe,var_fxlong,grid%fxlong) 
  call read_fire_var(adios2_id,info,iframe,var_fxlat,grid%fxlat) 
  call read_fire_var(adios2_id,info,iframe,var_fmc_g,grid%fmc_g) 
  !call read_fire_var(adios2_id,info,iframe,var_unit_fxlong,grid%unit_fxlong) 
  !call read_fire_var(adios2_id,info,iframe,var_unit_fxlat,grid%unit_fxlat) 
  call read_fire_var(adios2_id,info,iframe,var_uf,grid%uf) 
  call read_fire_var(adios2_id,info,iframe,var_vf,grid%vf) 

  iframe = iframe +1

  call adios2_end_step(adios2_id, ierr)

  !call adios2_close(adios2_id, ierr)
end subroutine read_step_vars


subroutine check(ierr,msg,cont)
  implicit none
  integer,intent(in)::ierr
  character(len=*), optional, intent(in)::msg
  logical, intent(in), optional::cont
  character(len=128)::message
  if(ierr.ne.(adios2_error_none))then
    write(*,*) "ADIOS2 error, type=", ierr
    !write(6,'(2I)')"ADIOSs error: ",ierr
    if(present(msg))then
      message=msg
    else
      message="ADIOS2 ERROR"
    endif
    if(present(cont))then
      if(cont)return
    endif
    call crash(message)
  endif
  end subroutine check

subroutine init_adios2(filename, adios2_id, ierr)
  implicit none
  character(len=*), intent(in)::filename
  type(adios2_engine), intent(out) :: adios2_id ! engine
  integer, intent(out) :: ierr
  integer comm_size, comm_rank

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, comm_rank, ierr)
  call adios2_init(adios, 'adios2.xml', MPI_COMM_WORLD, ierr)
  call adios2_declare_io(ioGet, adios, 'fire_input', ierr)

  call adios2_open(adios2_id, ioGet, filename, adios2_mode_read, ierr )
  call check(ierr)
end subroutine init_adios2

subroutine finalize_adios2(adios2_id, ierr)
  implicit none
  type(adios2_engine), intent(inout) :: adios2_id ! engine
  integer, intent(out) :: ierr

  call adios2_close(adios2_id, ierr)
  call check(ierr)
  call adios2_finalize(adios, ierr)
end subroutine finalize_adios2


end module wrf_adios2
