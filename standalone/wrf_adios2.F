!    files and reads input files that are somewhat compatible with wrf files
!

module wrf_adios2
use adios2
use module_fr_sfire_util , only : crash, interpolate_2d, continue_at_boundary
use module_domain, only: domain
use wrf_netcdf, only : grid_info, nf90_max_name, nf90_float
implicit none

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


type(adios2_adios):: adios
type(adios2_io) :: ioGet



contains

!
!*** READING 
!

subroutine read_info(filename,info)
  ! get fire grid sizes from input file
  implicit none


  !*** arguments
  character(len=*),intent(in)::filename
  type(grid_info),intent(inout)::info

  !*** local

  type(adios2_engine) :: adios2_id

  integer::ierr,it,ndims,idim,ratio,stagid,i
  integer,dimension(mdims)::dims,sr,dimids,stag
  real::dx,dy
  character(len=128)::msg
  character(len=NF90_MAX_NAME)::dimnames(mdims),dimname,stagname,stagnames(mdims)
  character(len=*), parameter:: subgrid='_subgrid'

  !*** executable
  call adios2_declare_io(ioGet, adios, 'WRFOUT', ierr)

  ! Feng:
  ! First read all dims information, look at VAR_NAME/Dims

  ! get time string information
  ! We don't have this dim, instead read
  call read_var_info(filename,var_times,ndims,dims,dimnames)
  if(ndims.ne.2)call crash(trim(var_times)//' must have 2 dimensions')
  info%dim_time_string = dimnames(1)
  info%len_time_string = dims(1)
  if(dims(1).gt.max_times_length)call crash('time string too long')
    
  ! get dimensions
  ! NFUEL_CAT
  call read_var_info(filename,var_nfuel_cat,ndims,dims,dimnames)

  ! store in info what you can now
  info%io_nfirex=dims(1)
  info%io_nfirey=dims(2)
  info%ntimes=dims(3)
  info%dim_fire_x=dimnames(1)
  info%dim_fire_y=dimnames(2)
  info%dim_time=dimnames(3)

  ! attributes

  call adios2_open(adios2_id, ioGet, filename, adios2_mode_read, MPI_COMM_SELF, ierr )

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
      if ( i + len(subgrid) .eq. len_trim(dimname)+1) then ! subgrid, need to fix
          stagname=dimname(1:i-1) // '_stag'
          ! TO read dim information from attr
          call att_read(adios2_id, '__DIM__'//stagname, stag(idim))
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

  call check(adios2_close(adios2_id, ierr))

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

end subroutine read_info

subroutine dim_read(adios2_id,dim_name,dim_len)
character(len=*),intent(in)::dim_name
integer, intent(in)::adios2_id
integer,intent(out):: dim_len
integer:: dim_id

call att_read(adios2_id, '__DIM__'//stagname, stag(idim))
!call check(nf90_inq_dimid(adios2_id,dim_name,dim_id))
!call check(nf90_inquire_dimension(adios2_id,dim_id,len=dim_len))
if(debug_print)write(*,'(a,1x,a,i6)')'dimension',trim(dim_name),dim_len
end subroutine dim_read

subroutine att_read(adios2_id,att_name,att_val)
  integer, intent(in)::adios2_id
  character(len=*),intent(in)::att_name
  real, intent(out)::att_val

  integer::ierr
  !call check(nf90_get_att(adios2_id,nf90_global,att_name,att_val))
  call check(adios2_inquire_attribute(attr_val, adios2_io, att_name, ierr))
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
  if(debug_print)write(*,'(4a)')'reading file ',filename,' dimensions of variable ',trim(varname)
  call read_var_info(filename,varname,ndims,dimlengths,dimnames)
  if(debug_print)write(*,'(3a,4(2a,i5,1x))')'variable ',trim(varname),' dimensions ', &
    (trim(dimnames(i)),'=',dimlengths(i),i=1,ndims)

end subroutine print_var_info


subroutine read_var_info(filename,varname,ndims,dimlengths,dimnames,type)

  ! get variable dimensions from a file
  implicit none

  ! arguments
  character(len=*), intent(in):: filename ! variable name
  character(len=*), intent(in):: varname ! variable name
  integer, intent(out)::ndims            ! number of dimensions of this variable
  integer, intent(out)::dimlengths(mdims)      ! the dimensions
  character(len=NF90_MAX_NAME),intent(out)::dimnames(mdims) ! dimension names
  integer, intent(out), optional :: type

  ! local
  type(adios2_engine) :: adios2_id
  type(adios2_variable) :: adios_var
  integer::ierr

  integer::varid,dimid,idimids(mdims),idim,i,dimlen,xtype
  character(len=NF90_MAX_NAME)::dimname

  ! executable
  if(debug_print)write(*,'(4a)')'reading file ',trim(filename),' variable ',trim(varname)
  call adios2_open(adios2_id, ioGet, filename, adios2_mode_read, MPI_COMM_SELF, ierr )

  !call check(nf90_inq_varid(adios2_id, varname, varid))
  !call check(nf90_inquire_variable(adios2_id,varid,ndims=ndims,xtype=xtype))
  call adios2_inquire_variable(adios2_var, ioGet, varname, ierr)
  if (adios2_var % valid == .false.)
    call crash('variable not found')
  call adios2_variable_ndims(ndims, adios2_var, ierr)
  call adios2_variable_type(xtype, adios2_var, ierr)


  if(ndims.gt.mdims)then
      write(*,1)'variable ',trim(varname),' has ',ndims,' dimensions >',mdims
      call crash('variable has too many dimensions')
  1 format(3a,i5,a,(4i6))
  endif

  call check(adios2_inquire_attribute(dimnames, adios2_io, trim(varname)//'/Dims', ierr))
  !call check(nf90_inquire_variable(adios2_id, varid, dimids=idimids))
  do idim=1,ndims
      dimname = dimnames(idim)
      if(debug_print)write(*,'(a,i3,a,i6)')'inquiring dimension ',idim,' name', dimname
      
      adios2_inquire_attribute(dimlen, adios2_io, '_DIMS_'//dimnames(idim), ierr

      if(debug_print)write(*,'(3a,i5)')'got dimension length',dimlen
      dimlengths(idim) = dimlen  
  enddo
  call check(adios2_close(adios2_id, ierr))
  if(present(type)) type=xtype
  2 format(3a,i4,a,(4i6))
  if(debug_print)write(*,2)'variable ',trim(varname), &
    ' type',xtype,' dimensions ',(dimlengths(idim),idim=1,ndims)
end subroutine read_var_info

!
!***
!

! get a real value (real number)
subroutine read_real(adios2_id,iframe,varname,v)
  implicit none
  !*** arguments
  type(adios2_engine), intent(in)::adios2_id             ! id of netcdf file open in data mode
  integer, intent(in):: iframe          ! number of frame in the file
  character(len=*),intent(in)::varname  ! the variable name
  real :: v                ! value
  !*** local
  integer::varid
  real:: val(1)

  
  type(adios2_variable) :: adios_var 

  !*** executable

  if(debug_print)write(*,'(3a,i5)')'reading real ',trim(varname),' timestep ',iframe
  
  
 
  !call check(nf90_inq_varid(adios2_id,varname,varid),'cannot find '//trim(varname))
  !call check(nf90_get_var(adios2_id,varid,val,start=(/iframe/),count=(/1/)), &
  !    'error reading '//trim(varname))
  call check(adios2_inquire_variable(adios_var, ioGet, varname, ierr))
  call adios2_set_selection(adios_var, 1, /1/, /1/, ierr)
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
  integer, intent(in):: iframe          ! number of frame in the file
  character(len=*),intent(in)::varname  ! the variable name
  real, pointer :: v(:,:)            ! values
  !*** local
  integer::varid,nx,ny,start(3),count(3)
  type(adios2_variable) :: adios_var

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
  
  call check(adios2_inquire_variable(adios_var, adios2_id, varname, ierr))
  call adios2_set_selection(adios_var, 2, /1,1/, /nx, ny/, ierr)
  call adios2_get(adios2_id, adios_var, val, ierr)
  
  end subroutine read_fire_var

subroutine read_vars(filename,info,iframe,grid)
  ! read all variables from input file
  implicit none

  !*** arguments
  character(len=*),intent(in)::filename ! the input file
  type(grid_info),intent(inout)::info   ! dimensions
  integer, intent(in):: iframe          ! number of frame in the file (not used by ADIOS2)
  type(domain),intent(inout)::grid      ! the mother of all arrays

  !*** local
  integer::varid,ierr
  integer,dimension(4)::s,c

  type(adios2_engine) :: adios2_id
  type(adios2_variable) :: adios_var

  !*** executable
  if(debug_print_f)write(*,'(3a,i4)')'read_vars: reading file ',trim(filename),' frame ',iframe 
  call adios2_open(adios2_id, ioGet, filename, adios2_mode_read, MPI_COMM_SELF, ierr)

  !call check(nf90_inq_varid(adios2_id,var_times,varid),'cannot find '//trim(var_times))
  !call check(nf90_get_var(adios2_id,varid,info%times,start=(/1,iframe/),count=(/info%len_time_string,1/)), &
  !    'error reading '//trim(var_times))
  call adios2_inquire_attribute(info%times, adios2_io, var_times , ierr
  if(debug_print_f)write(*,'(2a)')'Time ',info%times

  call adios2_begin_step()

  
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

  call adios2_end_step(adios2_id, ierr)

  call check(adios2_close(adios2_id, ierr))
end subroutine read_vars

!
!***
!

subroutine check(adios2_err,msg,cont)
  implicit none
  integer,intent(in)::adios2_err
  character(len=*), optional, intent(in)::msg
  logical, intent(in), optional::cont
  character(len=128)::message
  if(adios2_err.ne.adios2_error_none)then
    ! write(6,'(2a)')"NetCDF error: ",trim(nf90_strerror(ncerr))
    write(6,'(2a)')"ADIOS2 error: ",adios2_err)
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

!
!***
!

subroutine init_adios2(adios)
  implicit none
  type(adios2_adios), intent(inout) :: adios
  integer, intent(out) :: ierr
  call adios2_init(adios, MPI_COMM_WORLD, ierr)

subroutine finalize_adios2(adios2)
  implicit none
  type(adios2_adios), intent(in) :: adios
  integer, intent(out) :: ierr

  call adios2_finalize(adios, ierr)


end module wrf_adios2
