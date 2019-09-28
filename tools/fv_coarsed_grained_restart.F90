module fv_coarse_grained_restart_mod

  use fms_io_mod,      only: register_restart_field, set_filename_appendix, restart_file_type, save_restart
  use fms_io_mod,          only: set_domain
  use field_manager_mod,       only: MODEL_ATMOS
  use fv_arrays_mod,   only: fv_atmos_type
  use mpp_mod,         only: mpp_error, mpp_npes, FATAL
  use mpp_domains_mod, only : mpp_get_ntile_count, domain2d, mpp_get_compute_domain, mpp_define_io_domain
  use mpp_domains_mod, only: domain2d, EAST, WEST, NORTH, CENTER, SOUTH, CORNER
  use time_manager_mod, only: time_type
  use tracer_manager_mod, only: tr_get_tracer_names=>get_tracer_names, &
                                     get_tracer_names, &
                                     set_tracer_profile
  use external_sst_mod,        only: use_ncep_sst
  
  use fv_coarse_grained_diagnostics_mod, only: define_cubic_mosaic
  use fv_coarse_grained_diagnostics_mod, only: coarse_grain_variable, area_and_delp_coarse_grain_variable, zonal_block_edge_coarse_grain_variable, meridional_block_edge_coarse_grain_variable

  implicit none
  private

  public :: fv_io_init_coarse, fv_io_register_restart_coarse, fv_io_write_restart_coarse
  
  integer :: is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse, npz  
  logical :: initialized = .false.
  
  contains
  
  subroutine fv_io_init_coarse(Atm)

    type(fv_atmos_type), intent(in) :: Atm(:)

    npz = Atm(1)%npz
    call mpp_get_compute_domain(Atm(1)%domain, is, ie, js, je)
    call mpp_get_compute_domain(Atm(1)%coarse_restart%coarse_domain, is_coarse, ie_coarse, &
         js_coarse, je_coarse)

    initialized = .true.
    
  end subroutine fv_io_init_coarse

  subroutine  fv_io_register_restart_coarse(Atm)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    
    character(len=64) :: fname, tracer_name
    character(len=6)  :: gn, stile_name
    integer           :: id_restart
    integer           :: n, nt, ntracers, ntprog, ntdiag, ntileMe, ntiles

    real, allocatable :: u_coarse(:,:,:), v_coarse(:,:,:), w_coarse(:,:,:),&
         delz_coarse(:,:,:), ze0_coarse(:,:,:), pt_coarse(:,:,:), &
         delp_coarse(:,:,:), phis_coarse(:,:), ua_coarse(:,:,:), &
         va_coarse(:,:,:), u_srf_coarse(:,:), v_srf_coarse(:,:), &
         sgh_coarse(:,:), oro_coarse(:,:), q_coarse(:,:,:), qdiag_coarse(:,:,:)

    write(*,*) 'coarse restart initialized?', initialized
    ! call set_domain(coarse_domain)
    
    ntileMe = size(Atm(:)) 
    ntprog = size(Atm(1)%q,4) 
    ntdiag = size(Atm(1)%qdiag,4) 
    ntracers = ntprog+ntdiag

    !--- set the 'nestXX' appendix for all files using fms_io
    if (Atm(1)%grid_number > 1) then
       write(gn,'(A4, I2.2)') "nest", Atm(1)%grid_number
    else
       gn = ''
    end if
    call set_filename_appendix(gn)

    !--- fix for single tile runs where you need fv_core.res.nc and fv_core.res.tile1.nc
    ntiles = mpp_get_ntile_count(Atm(1)%coarse_restart%coarse_domain)
    if(ntiles == 1 .and. .not. Atm(1)%neststruct%nested) then
       stile_name = '.tile1'
    else
       stile_name = ''
    endif

    fname = 'fv_core_coarse.res.nc'
    id_restart = register_restart_field(Atm(1)%fv_restart_coarse, fname, 'ak', Atm(1)%ak(:), no_domain=.true.)
    id_restart = register_restart_field(Atm(1)%fv_restart_coarse, fname, 'bk', Atm(1)%bk(:), no_domain=.true.) 

    do n = 1, ntileMe
       fname = 'fv_core_coarse.res'//trim(stile_name)//'.nc'

       id_restart =  register_restart_field(Atm(n)%fv_tile_restart_coarse, fname, 'u', Atm(n)%coarse_restart%u, &
                     domain=Atm(n)%coarse_restart%coarse_domain, position=NORTH,tile_count=n)
       id_restart =  register_restart_field(Atm(n)%fv_tile_restart_coarse, fname, 'v', Atm(n)%coarse_restart%v, &
                     domain=Atm(n)%coarse_restart%coarse_domain, position=EAST,tile_count=n)
       
       if (.not.Atm(n)%flagstruct%hydrostatic) then
          id_restart =  register_restart_field(Atm(n)%fv_tile_restart_coarse, fname, 'W', Atm(n)%coarse_restart%w, &
                        domain=Atm(n)%coarse_restart%coarse_domain, mandatory=.false., tile_count=n)
          id_restart =  register_restart_field(Atm(n)%fv_tile_restart_coarse, fname, 'DZ', Atm(n)%coarse_restart%delz, &
                        domain=Atm(n)%coarse_restart%coarse_domain, mandatory=.false., tile_count=n)
          
          ! if ( Atm(n)%flagstruct%hybrid_z ) then

          !    ! allocate(ze0_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
             
          !    ! TODO: is this the right method of coarse graining?
          !    ! call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
          !    !     Atm(n)%ze0(is:ie,js:je,1:npz), Atm(n)%coarse_restart%ze0, coarsening_factor)
      
          !    id_restart =  register_restart_field(Atm(n)%fv_tile_restart_coarse, fname, 'ZE0', Atm(n)%coarse_restart%ze0, &
          !         domain=coarse_domain, mandatory=.false., tile_count=n)

          !    ! deallocate(ze0_coarse)
          ! endif
       endif
       
       id_restart =  register_restart_field(Atm(n)%fv_tile_restart_coarse, fname, 'T', Atm(n)%coarse_restart%pt, &
                     domain=Atm(n)%coarse_restart%coarse_domain, tile_count=n)
       id_restart =  register_restart_field(Atm(n)%fv_tile_restart_coarse, fname, 'delp', Atm(n)%coarse_restart%delp, &
                     domain=Atm(n)%coarse_restart%coarse_domain, tile_count=n)
       id_restart =  register_restart_field(Atm(n)%fv_tile_restart_coarse, fname, 'phis', Atm(n)%coarse_restart%phis, &
                     domain=Atm(n)%coarse_restart%coarse_domain, tile_count=n)
       
       !--- include agrid winds in restarts for use in data assimilation 
       if (Atm(n)%flagstruct%agrid_vel_rst) then
          
         id_restart =  register_restart_field(Atm(n)%fv_tile_restart_coarse, fname, 'ua', Atm(n)%coarse_restart%ua, &
                       domain=Atm(n)%coarse_restart%coarse_domain, tile_count=n, mandatory=.false.)
         id_restart =  register_restart_field(Atm(n)%fv_tile_restart_coarse, fname, 'va', Atm(n)%coarse_restart%va, &
                       domain=Atm(n)%coarse_restart%coarse_domain, tile_count=n, mandatory=.false.)

      endif
      
       fname = 'fv_srf_wnd_coarse.res'//trim(stile_name)//'.nc'
       id_restart =  register_restart_field(Atm(n)%rsf_restart_coarse, fname, 'u_srf', Atm(n)%coarse_restart%u_srf, &
                     domain=Atm(n)%coarse_restart%coarse_domain, tile_count=n)
       id_restart =  register_restart_field(Atm(n)%rsf_restart_coarse, fname, 'v_srf', Atm(n)%coarse_restart%v_srf, &
                     domain=Atm(n)%coarse_restart%coarse_domain, tile_count=n)
       
       if ( Atm(n)%flagstruct%fv_land ) then          
          !-------------------------------------------------------------------------------------------------
          ! Optional terrain deviation (sgh) and land fraction (oro)
          ! allocate(sgh_coarse(is_coarse:ie_coarse,js_coarse:je_coarse))
          ! allocate(oro_coarse(is_coarse:ie_coarse,js_coarse:je_coarse))

          fname = 'mg_drag_coarse.res'//trim(stile_name)//'.nc'
          id_restart =  register_restart_field(Atm(n)%mg_restart_coarse, fname, 'ghprime', Atm(n)%coarse_restart%sgh, &
                        domain=Atm(n)%coarse_restart%coarse_domain, tile_count=n)  

          fname = 'fv_land_coarse.res'//trim(stile_name)//'.nc'
          id_restart = register_restart_field(Atm(n)%lnd_restart_coarse, fname, 'oro', Atm(n)%coarse_restart%oro, &
                        domain=Atm(n)%coarse_restart%coarse_domain, tile_count=n)
       endif
       
       fname = 'fv_tracer_coarse.res'//trim(stile_name)//'.nc'

       do nt = 1, ntprog
          call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
          call set_tracer_profile (MODEL_ATMOS, nt, Atm(n)%coarse_restart%q(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,nt))
          id_restart = register_restart_field(Atm(n)%tra_restart_coarse, fname,tracer_name, Atm(n)%coarse_restart%q(:,:,:,nt), &
               domain=Atm(n)%coarse_restart%coarse_domain, mandatory=.false., tile_count=n)
       enddo
       
       do nt = ntprog+1, ntracers
          call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
          call set_tracer_profile (MODEL_ATMOS, nt, Atm(n)%coarse_restart%qdiag(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,nt))
          id_restart = register_restart_field(Atm(n)%tra_restart_coarse, fname, tracer_name, Atm(n)%coarse_restart%qdiag(:,:,:,nt), &
               domain=Atm(n)%coarse_restart%coarse_domain, mandatory=.false., tile_count=n)
       enddo

    enddo

    ! call set_domain(Atm(1)%domain)
    
  end subroutine  fv_io_register_restart_coarse

  subroutine fv_io_write_restart_coarse(Atm, timestamp)
    type(fv_atmos_type),        intent(inout) :: Atm
    character(len=*), optional, intent(in) :: timestamp

    ! TODO: do we need to support this?
    ! if ( (use_ncep_sst .or. Atm%flagstruct%nudge) .and. .not. Atm%gridstruct%nested ) then
    !    call save_restart(Atm%SST_restart, timestamp)
    ! endif

    call coarse_grain_restart_data(Atm)
    
    call save_restart(Atm%fv_restart_coarse, timestamp)
    call save_restart(Atm%fv_tile_restart_coarse, timestamp)
    call save_restart(Atm%rsf_restart_coarse, timestamp)

    if ( Atm%flagstruct%fv_land ) then
       call save_restart(Atm%mg_restart_coarse, timestamp)
       call save_restart(Atm%lnd_restart_coarse, timestamp)
    endif

    call save_restart(Atm%tra_restart_coarse, timestamp)
    
  end subroutine fv_io_write_restart_coarse

  subroutine coarse_grain_restart_data(Atm)
    type(fv_atmos_type) :: Atm
    integer :: ntprog, ntdiag, ntracers, nt
    
    ntprog = size(Atm%q,4) 
    ntdiag = size(Atm%qdiag,4) 
    ntracers = ntprog+ntdiag
    
    call zonal_block_edge_coarse_grain_variable(Atm%gridstruct%dx(is:ie,js:je+1), &
          Atm%u(is:ie,js:je+1,1:npz), Atm%coarse_restart%u, Atm%coarse_restart%coarsening_factor)
    call meridional_block_edge_coarse_grain_variable(Atm%gridstruct%dy(is:ie+1,js:je), &
          Atm%v(is:ie+1,js:je,1:npz), Atm%coarse_restart%v, Atm%coarse_restart%coarsening_factor)
    call area_and_delp_coarse_grain_variable(Atm%gridstruct%area(is:ie,js:je), &
         Atm%delp(is:ie,js:je,1:npz), Atm%pt(is:ie,js:je,1:npz), &
         Atm%coarse_restart%pt, Atm%coarse_restart%coarsening_factor)

    if (.not.Atm%flagstruct%hydrostatic) then
       call area_and_delp_coarse_grain_variable(Atm%gridstruct%area(is:ie,js:je), &
            Atm%delp(is:ie,js:je,1:npz), Atm%w(is:ie,js:je,1:npz), &
            Atm%coarse_restart%w, Atm%coarse_restart%coarsening_factor)
       call coarse_grain_variable(Atm%gridstruct%area(is:ie,js:je), &
            Atm%delz(is:ie,js:je,1:npz), Atm%coarse_restart%delz, Atm%coarse_restart%coarsening_factor)
    endif

    call area_and_delp_coarse_grain_variable(Atm%gridstruct%area(is:ie,js:je), &
         Atm%delp(is:ie,js:je,1:npz), Atm%pt(is:ie,js:je,1:npz), &
         Atm%coarse_restart%pt, Atm%coarse_restart%coarsening_factor)
    call coarse_grain_variable(Atm%gridstruct%area(is:ie,js:je), &
         Atm%delp(is:ie,js:je,1:npz), Atm%coarse_restart%delp, Atm%coarse_restart%coarsening_factor)
    call coarse_grain_variable(Atm%gridstruct%area(is:ie,js:je), &
         Atm%phis(is:ie,js:je), Atm%coarse_restart%phis, Atm%coarse_restart%coarsening_factor)

    if (Atm%flagstruct%agrid_vel_rst) then      
       call coarse_grain_variable(Atm%gridstruct%area(is:ie,js:je), &
            Atm%ua(is:ie,js:je,1:npz), Atm%coarse_restart%ua, Atm%coarse_restart%coarsening_factor)
       call coarse_grain_variable(Atm%gridstruct%area(is:ie,js:je), &
            Atm%va(is:ie,js:je,1:npz), Atm%coarse_restart%va, Atm%coarse_restart%coarsening_factor)
    endif

    call coarse_grain_variable(Atm%gridstruct%area(is:ie,js:je), &
         Atm%u_srf(is:ie,js:je), Atm%coarse_restart%u_srf, Atm%coarse_restart%coarsening_factor)
    call coarse_grain_variable(Atm%gridstruct%area(is:ie,js:je), &
         Atm%v_srf(is:ie,js:je), Atm%coarse_restart%v_srf, Atm%coarse_restart%coarsening_factor)

    if ( Atm%flagstruct%fv_land ) then
       call coarse_grain_variable(Atm%gridstruct%area(is:ie,js:je), &
            Atm%sgh(is:ie,js:je), Atm%coarse_restart%sgh, Atm%coarse_restart%coarsening_factor)
       call coarse_grain_variable(Atm%gridstruct%area(is:ie,js:je), &
            Atm%oro(is:ie,js:je), Atm%coarse_restart%oro, Atm%coarse_restart%coarsening_factor)
    endif

    do nt = 1, ntprog
       call area_and_delp_coarse_grain_variable(Atm%gridstruct%area(is:ie,js:je), &
            Atm%delp(is:ie,js:je,1:npz), Atm%q(is:ie,js:je,1:npz,nt), &
            Atm%coarse_restart%q(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,nt), Atm%coarse_restart%coarsening_factor)
    enddo
    
    do nt = ntprog+1, ntracers
       call area_and_delp_coarse_grain_variable(Atm%gridstruct%area(is:ie,js:je), &
            Atm%delp(is:ie,js:je,1:npz), Atm%qdiag(is:ie,js:je,1:npz,nt), &
            Atm%coarse_restart%qdiag(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,nt), Atm%coarse_restart%coarsening_factor)
    enddo
    
  end subroutine coarse_grain_restart_data
  
end module fv_coarse_grained_restart_mod
