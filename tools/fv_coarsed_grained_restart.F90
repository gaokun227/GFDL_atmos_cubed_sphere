module fv_coarse_grained_restart_mod

  use fms_io_mod,      only: register_restart_field, set_filename_appendix, restart_file_type, save_restart
  use field_manager_mod,       only: MODEL_ATMOS
  use fv_arrays_mod,   only: fv_atmos_type
  use mpp_mod,         only: mpp_error, mpp_npes, FATAL
  use mpp_domains_mod, only : mpp_get_ntile_count, domain2d, mpp_get_compute_domain
  use mpp_domains_mod, only: domain2d, EAST, WEST, NORTH, CENTER, SOUTH, CORNER
  use time_manager_mod, only: time_type
  use tracer_manager_mod, only: tr_get_tracer_names=>get_tracer_names, &
                                     get_tracer_names, &
                                     set_tracer_profile
  use external_sst_mod,        only: use_ncep_sst
  
  use fv_coarse_grained_diagnostics_mod, only: define_cubic_mosaic
  use fv_coarse_grained_diagnostics_mod, only: coarse_grain_variable, zonal_block_edge_coarse_grain_variable, meridional_block_edge_coarse_grain_variable

  implicit none
  private

  public :: fv_io_init_coarse, fv_io_register_restart_coarse, fv_io_write_restart_coarse
  
  type(restart_file_type) :: fv_restart_coarse, fv_tile_restart_coarse, fv_mg_restart_coarse
  type(restart_file_type) :: fv_tra_restart_coarse, fv_lnd_restart_coarse, fv_srf_restart_coarse
  type(domain2d) :: coarse_domain
  integer :: coarsening_factor, target_resolution
  integer :: is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse, npz  

  contains
  
  subroutine fv_io_init_coarse(Atm, cf)

    type(fv_atmos_type), intent(in) :: Atm(:)
    integer, intent(in) :: cf

    integer :: layout(2)
    integer :: n, npx, native_resolution

    n = 1
    npx = Atm(n)%gridstruct%npx_g
    coarsening_factor = cf
    native_resolution = npx - 1
    target_resolution = native_resolution / coarsening_factor
    layout = Atm(n)%layout
    call define_cubic_mosaic(coarse_domain, target_resolution, target_resolution, layout)

    npz = Atm(n)%npz
    call mpp_get_compute_domain(Atm(n)%domain, is, ie, js, je)
    call mpp_get_compute_domain(coarse_domain, is_coarse, ie_coarse, &
         js_coarse, je_coarse)
    
  end subroutine fv_io_init_coarse

  subroutine  fv_io_register_restart_coarse(Atm, coarsening_factor)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    integer, intent(in) :: coarsening_factor
    
    character(len=64) :: fname, tracer_name
    character(len=6)  :: gn, stile_name
    integer           :: id_restart
    integer           :: n, nt, ntracers, ntprog, ntdiag, ntileMe, ntiles

    real, allocatable :: u_coarse(:,:,:), v_coarse(:,:,:), w_coarse(:,:,:),&
         delz_coarse(:,:,:), ze0_coarse(:,:,:), pt_coarse(:,:,:), &
         delp_coarse(:,:,:), phis_coarse(:,:), ua_coarse(:,:,:), &
         va_coarse(:,:,:), u_srf_coarse(:,:), v_srf_coarse(:,:), &
         sgh_coarse(:,:), oro_coarse(:,:), q_coarse(:,:,:), qdiag_coarse(:,:,:)
    
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
    ntiles = mpp_get_ntile_count(coarse_domain)
    if(ntiles == 1 .and. .not. Atm(1)%neststruct%nested) then
       stile_name = '.tile1'
    else
       stile_name = ''
    endif

    fname = 'fv_core_coarse.res.nc'
    id_restart = register_restart_field(fv_restart_coarse, fname, 'ak', Atm(1)%ak(:), no_domain=.true.)
    id_restart = register_restart_field(fv_restart_coarse, fname, 'bk', Atm(1)%bk(:), no_domain=.true.) 

    do n = 1, ntileMe
       fname = 'fv_core_coarse.res'//trim(stile_name)//'.nc'

       allocate(u_coarse(is_coarse:ie_coarse,js_coarse:je_coarse+1,1:npz))
       allocate(v_coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse,1:npz))

       call zonal_block_edge_coarse_grain_variable(Atm(n)%gridstruct%dx(is:ie,js:je), &
            Atm(n)%u(is:ie,js:je+1,1:npz), u_coarse, coarsening_factor)
       call meridional_block_edge_coarse_grain_variable(Atm(n)%gridstruct%dy(is:ie,js:je), &
            Atm(n)%v(is:ie+1,js:je,1:npz), v_coarse, coarsening_factor)
       
       id_restart =  register_restart_field(fv_tile_restart_coarse, fname, 'u', u_coarse, &
                     domain=coarse_domain, position=NORTH,tile_count=n)
       id_restart =  register_restart_field(fv_tile_restart_coarse, fname, 'v', v_coarse, &
                     domain=coarse_domain, position=EAST,tile_count=n)

       deallocate(u_coarse, v_coarse)
       
       if (.not.Atm(n)%flagstruct%hydrostatic) then

          ! TODO: is this the right method of coarse graining for these variables?
          allocate(w_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
          allocate(delz_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))

          call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
               Atm(n)%w(is:ie,js:je,1:npz), w_coarse, coarsening_factor)
          call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
               Atm(n)%delz(is:ie,js:je,1:npz), delz_coarse, coarsening_factor)
          
          id_restart =  register_restart_field(fv_tile_restart_coarse, fname, 'W', w_coarse, &
                        domain=coarse_domain, mandatory=.false., tile_count=n)
          id_restart =  register_restart_field(fv_tile_restart_coarse, fname, 'DZ', delz_coarse, &
                        domain=coarse_domain, mandatory=.false., tile_count=n)

          deallocate(w_coarse, delz_coarse)
          
          if ( Atm(n)%flagstruct%hybrid_z ) then

             allocate(ze0_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
             
             ! TODO: is this the right method of coarse graining?
             call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
                  Atm(n)%ze0(is:ie,js:je,1:npz), ze0_coarse, coarsening_factor)
      
             id_restart =  register_restart_field(fv_tile_restart_coarse, fname, 'ZE0', ze0_coarse, &
                  domain=coarse_domain, mandatory=.false., tile_count=n)

             deallocate(ze0_coarse)
          endif
       endif
       
       allocate(pt_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
       allocate(delp_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
       allocate(phis_coarse(is_coarse:ie_coarse,js_coarse:je_coarse))
       
       call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
            Atm(n)%pt(is:ie,js:je,1:npz), pt_coarse, coarsening_factor)
       call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
            Atm(n)%delp(is:ie,js:je,1:npz), delp_coarse, coarsening_factor)
       call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
            Atm(n)%phis(is:ie,js:je), phis_coarse, coarsening_factor)
       
       id_restart =  register_restart_field(fv_tile_restart_coarse, fname, 'T', pt_coarse, &
                     domain=coarse_domain, tile_count=n)
       id_restart =  register_restart_field(fv_tile_restart_coarse, fname, 'delp', delp_coarse, &
                     domain=coarse_domain, tile_count=n)
       id_restart =  register_restart_field(fv_tile_restart_coarse, fname, 'phis', phis_coarse, &
                     domain=coarse_domain, tile_count=n)

       deallocate(pt_coarse, delp_coarse, phis_coarse)
       
       !--- include agrid winds in restarts for use in data assimilation 
       if (Atm(n)%flagstruct%agrid_vel_rst) then

          allocate(ua_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
          allocate(va_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
          
          call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
               Atm(n)%ua(is:ie,js:je,1:npz), ua_coarse, coarsening_factor)
          call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
               Atm(n)%va(is:ie,js:je,1:npz), va_coarse, coarsening_factor)
          
         id_restart =  register_restart_field(fv_tile_restart_coarse, fname, 'ua', ua_coarse, &
                       domain=coarse_domain, tile_count=n, mandatory=.false.)
         id_restart =  register_restart_field(fv_tile_restart_coarse, fname, 'va', va_coarse, &
                       domain=coarse_domain, tile_count=n, mandatory=.false.)

         deallocate(ua_coarse, va_coarse)
      endif

      allocate(u_srf_coarse(is_coarse:ie_coarse,js_coarse:je_coarse))
      allocate(v_srf_coarse(is_coarse:ie_coarse,js_coarse:je_coarse))

      call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
           Atm(n)%u_srf(is:ie,js:je), u_srf_coarse, coarsening_factor)
      call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
           Atm(n)%v_srf(is:ie,js:je), v_srf_coarse, coarsening_factor)
      
       fname = 'fv_srf_wnd_coarse.res'//trim(stile_name)//'.nc'
       id_restart =  register_restart_field(fv_srf_restart_coarse, fname, 'u_srf', u_srf_coarse, &
                     domain=coarse_domain, tile_count=n)
       id_restart =  register_restart_field(fv_srf_restart_coarse, fname, 'v_srf', v_srf_coarse, &
                     domain=coarse_domain, tile_count=n)

       deallocate(u_srf_coarse, v_srf_coarse)
       
       if ( Atm(n)%flagstruct%fv_land ) then          
          !-------------------------------------------------------------------------------------------------
          ! Optional terrain deviation (sgh) and land fraction (oro)
          allocate(sgh_coarse(is_coarse:ie_coarse,js_coarse:je_coarse))
          allocate(oro_coarse(is_coarse:ie_coarse,js_coarse:je_coarse))

          call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
               Atm(n)%sgh(is:ie,js:je), sgh_coarse, coarsening_factor)
          call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
               Atm(n)%oro(is:ie,js:je), oro_coarse, coarsening_factor)
          
          fname = 'mg_drag_coarse.res'//trim(stile_name)//'.nc'
          id_restart =  register_restart_field(fv_mg_restart_coarse, fname, 'ghprime', sgh_coarse, &
                        domain=coarse_domain, tile_count=n)  

          fname = 'fv_land_coarse.res'//trim(stile_name)//'.nc'
          id_restart = register_restart_field(fv_lnd_restart_coarse, fname, 'oro', oro_coarse, &
                        domain=coarse_domain, tile_count=n)

          deallocate(sgh_coarse, oro_coarse)
       endif
       
       fname = 'fv_tracer_coarse.res'//trim(stile_name)//'.nc'

       do nt = 1, ntprog
          allocate(q_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
          call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
               Atm(n)%q(is:ie,js:je,1:npz,nt), q_coarse, coarsening_factor)        
          call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
          ! set all tracers to an initial profile value
          call set_tracer_profile (MODEL_ATMOS, nt, q_coarse(:,:,:))
          id_restart = register_restart_field(fv_tra_restart_coarse, fname, tracer_name, q_coarse(:,:,:), &
               domain=coarse_domain, mandatory=.false., tile_count=n)
          deallocate(q_coarse)
       enddo
       
       do nt = ntprog+1, ntracers
          allocate(qdiag_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
          call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), &
               Atm(n)%qdiag(is:ie,js:je,1:npz,nt), qdiag_coarse, coarsening_factor)
          call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
          ! set all tracers to an initial profile value
          call set_tracer_profile (MODEL_ATMOS, nt, qdiag_coarse(:,:,:))
          id_restart = register_restart_field(fv_tra_restart_coarse, fname, tracer_name, qdiag_coarse(:,:,:), &
               domain=coarse_domain, mandatory=.false., tile_count=n)
          deallocate(qdiag_coarse)
       enddo

    enddo

  end subroutine  fv_io_register_restart_coarse

  subroutine fv_io_write_restart_coarse(Atm, timestamp)
    type(fv_atmos_type),        intent(inout) :: Atm
    character(len=*), optional, intent(in) :: timestamp

    ! TODO: do we need to support this?
    ! if ( (use_ncep_sst .or. Atm%flagstruct%nudge) .and. .not. Atm%gridstruct%nested ) then
    !    call save_restart(Atm%SST_restart, timestamp)
    ! endif

    call save_restart(fv_restart_coarse, timestamp)
    call save_restart(fv_tile_restart_coarse, timestamp)
    call save_restart(fv_srf_restart_coarse, timestamp)

    if ( Atm%flagstruct%fv_land ) then
       call save_restart(fv_mg_restart_coarse, timestamp)
       call save_restart(fv_lnd_restart_coarse, timestamp)
    endif

    call save_restart(fv_tra_restart_coarse, timestamp)
    
  end subroutine fv_io_write_restart_coarse
  
end module fv_coarse_grained_restart_mod
