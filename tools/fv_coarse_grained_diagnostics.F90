module fv_coarse_grained_diagnostics_mod

  use constants_mod,   only: pi=>pi_8
  use diag_manager_mod, only: diag_axis_init, register_diag_field, register_static_field, send_data
  use fv_arrays_mod,   only: fv_atmos_type, fv_grid_bounds_type, R_GRID
  use fv_grid_utils_mod, only: gnomonic_grids, cell_center2
  use fv_mp_mod,       only: domain_decomp
  use mpp_mod,         only: mpp_error, mpp_npes, FATAL
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_define_mosaic, domain2d, mpp_define_io_domain
  use time_manager_mod, only: time_type

  implicit none
  private
  
  public :: fv_coarse_grained_diagnostics_init, fv_coarse_grained_diagnostics

  interface coarse_grain_variable
     module procedure coarse_grain_variable_3d
     module procedure coarse_grain_variable_2d
  end interface coarse_grain_variable
  
  real, parameter:: missing_value = -1.e10
  real, parameter    ::     rad2deg = 180./pi
  integer :: id_ps_coarse, id_vort_coarse, id_coarse_lon, id_coarse_lat, id_coarse_lont, id_coarse_latt
  integer :: coarsening_factor
  integer :: is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse, npz

  real(kind=R_GRID), allocatable :: grid_coarse(:,:,:), agrid_coarse(:,:,:)
  contains
  
  subroutine fv_coarse_grained_diagnostics_init(Atm, Time, cf, id_pfull)
    type(time_type), intent(in) :: Time
    type(fv_atmos_type), intent(in) :: Atm(:)
    integer, intent(in) :: cf
    integer, intent(in) :: id_pfull
    type(domain2d)     :: coarse_domain
    integer :: coarse_axes(3), coarse_axes_t(3)
    integer :: npx, npy, target_resolution, native_resolution
    integer :: i, j, n, id_coarse_x, id_coarse_y, id_coarse_xt, id_coarse_yt, npx_target, npy_target
    integer :: offset
    integer :: layout(2), io_layout(2)
    character(len=256) :: errmsg

    real, allocatable :: grid_x_coarse(:), grid_y_coarse(:), grid_xt_coarse(:), grid_yt_coarse(:)
    
    n = 1
    call mpp_get_compute_domain(Atm(n)%domain, is, ie, js, je)
    coarsening_factor = cf
    
    ! TODO: add checks to make sure the coarsening factor is compatible with
    ! both the original resolution and the domain decomposition.
    npx = Atm(n)%gridstruct%npx_g
    npy = Atm(n)%gridstruct%npy_g
    npz = Atm(n)%npz
    layout = Atm(n)%layout
    io_layout = Atm(n)%io_layout
    native_resolution = npx - 1
    
    ! Check that the coarsening factor evenly divides the native domain
    if (mod(native_resolution, coarsening_factor) .ne. 0) then
       write(errmsg, *) 'fv_coarse_grained_diagnostics_init: coarsening factor ',&
            coarsening_factor, ' does not evenly divide native resolution',&
            native_resolution, '.'
       call mpp_error(FATAL, errmsg)
    endif
    
    target_resolution = native_resolution / coarsening_factor
    npx_target = target_resolution + 1
    npy_target = target_resolution + 1

    ! Check that the coarsening factor is appropriate for the domain
    ! decomposition
    if (mod(native_resolution, layout(1)) > 0 .or. mod(target_resolution, layout(1)) > 0) then
       write(errmsg, *) 'fv_coarse_grained_diagnostics_init: domain decomposition layout ',&
            layout(1), ' does not evenly divide both the native and coarse grid.'
       call mpp_error(FATAL, errmsg)
    endif
    
    allocate(grid_x_coarse(target_resolution + 1), grid_y_coarse(target_resolution + 1))
    grid_x_coarse = (/ (i, i=1, target_resolution + 1) /)
    grid_y_coarse = (/ (j, j=1, target_resolution + 1) /)
    
    allocate(grid_xt_coarse(target_resolution), grid_yt_coarse(target_resolution))
    grid_xt_coarse = (/ (i, i=1, target_resolution) /)
    grid_yt_coarse = (/ (j, j=1, target_resolution) /)
    
    call define_cubic_mosaic(coarse_domain, target_resolution, target_resolution, layout)
    call mpp_define_io_domain(coarse_domain, io_layout)
    call mpp_get_compute_domain(coarse_domain, is_coarse, ie_coarse, &
         js_coarse, je_coarse)
    
    id_coarse_x = diag_axis_init('grid_x_coarse', grid_x_coarse, &
         'degrees_E', 'x', 'Corner longitude', set_name='coarse_grid', &
         Domain2=coarse_domain, tile_count=n)
    id_coarse_y = diag_axis_init('grid_y_coarse', grid_y_coarse, &
         'degrees_N', 'y', 'Corner latitude', set_name='coarse_grid', &
         Domain2=coarse_domain, tile_count=n)
    
    id_coarse_xt = diag_axis_init('grid_xt_coarse', grid_xt_coarse, &
         'degrees_E', 'x', 'T-cell longitude', set_name='coarse_grid', &
         Domain2=coarse_domain, tile_count=n)
    id_coarse_yt = diag_axis_init('grid_yt_coarse', grid_yt_coarse, &
         'degrees_N', 'y', 'T-cell latitude', set_name='coarse_grid', &
         Domain2=coarse_domain, tile_count=n)

    coarse_axes(1) = id_coarse_x
    coarse_axes(2) = id_coarse_y
    coarse_axes(3) = id_pfull
    
    coarse_axes_t(1) = id_coarse_xt
    coarse_axes_t(2) = id_coarse_yt
    coarse_axes_t(3) = id_pfull
    
    id_ps_coarse = register_diag_field('dynamics', 'ps_coarse', &
         coarse_axes_t(1:2), Time, 'surface pressure', 'Pa', &
         missing_value=missing_value)

    id_vort_coarse = register_diag_field('dynamics', 'vort_coarse', &
         coarse_axes_t(1:3), Time, 'vorticity', 'per s', &
         missing_value=missing_value)

    ! Compute the longitude and latitude of the coarse grid based on values
    ! from fine grid.
    
    allocate(grid_coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse+1, 1:2))
    allocate(agrid_coarse(is_coarse:ie_coarse,js_coarse:je_coarse, 1:2))

    grid_coarse = Atm(n)%gridstruct%grid(is:ie+1:coarsening_factor,js:je+1:coarsening_factor,:)

    ! agrid_coarse is constructed differently depending on whether coarsening
    ! factor is even or odd (in the even case the values are taken from the
    ! native grid corners; in the odd case the values are taken from the native
    ! grid centers).
    offset = coarsening_factor / 2
    if (mod(coarsening_factor, 2) .eq. 0) then
       agrid_coarse = Atm(n)%gridstruct%grid(is+offset:ie+1:coarsening_factor,js+offset:je+1:coarsening_factor,:)
    else
       agrid_coarse = Atm(n)%gridstruct%agrid(is+offset:ie:coarsening_factor,js+offset:je:coarsening_factor,:)
    endif
    
    id_coarse_lon = register_static_field('dynamics', 'grid_lon_coarse', coarse_axes(1:2),  &
         'longitude', 'degrees_E')
    id_coarse_lat = register_static_field('dynamics', 'grid_lat_coarse', coarse_axes(1:2),  &
         'latitude', 'degrees_N')
    
    id_coarse_lont = register_static_field('dynamics', 'grid_lont_coarse', coarse_axes_t(1:2),  &
         'longitude', 'degrees_E')
    id_coarse_latt = register_static_field('dynamics', 'grid_latt_coarse', coarse_axes_t(1:2),  &
         'latitude', 'degrees_N')
    
    deallocate(grid_xt_coarse, grid_yt_coarse)
    
  end subroutine fv_coarse_grained_diagnostics_init
  
  subroutine coarse_grain_variable_2d(area, variable, coarse_grained_variable)

    real, intent(in) :: area(is:ie,js:je)
    real, intent(in) :: variable(is:ie,js:je)
    real, intent(out) :: coarse_grained_variable(is_coarse:ie_coarse,js_coarse:je_coarse)
 
    real, allocatable :: area_weighted_variable(:,:)
    integer :: i, j, n, i_coarse, j_coarse, a

    allocate(area_weighted_variable(is:ie,js:je))
    area_weighted_variable = area * variable
    
    n = 1
    a = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse_grained_variable(i_coarse,j_coarse) = sum(area_weighted_variable(i:i + a,j:j + a)) / sum(area(i:i + a,j:j + a))
       enddo
    enddo

    deallocate(area_weighted_variable)
    
  end subroutine coarse_grain_variable_2d

  subroutine coarse_grain_variable_3d(area, variable, coarse_grained_variable)

    real, intent(in) :: area(is:ie,js:je)
    real, intent(in) :: variable(is:ie,js:je,1:npz)
    real, intent(out) :: coarse_grained_variable(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)
 
    real, allocatable :: area_weighted_variable(:,:)
    integer :: k
    
    ! a = coarsening_factor - 1
    ! do k = 1, npz
    !    area_weighted_variable = area * variable(is:ie,js:je,k)
    !    do i = is, ie, coarsening_factor
    !       i_coarse = (i - 1) / coarsening_factor + 1
    !       do j = js, je, coarsening_factor
    !          j_coarse = (j - 1) / coarsening_factor + 1
    !          coarse_grained_variable(i_coarse,j_coarse,k) = sum(area_weighted_variable(i:i + a,j:j + a)) / sum(area(i:i + a,j:j + a))
    !       enddo
    !    enddo
    ! enddo

    do k = 1, npz
       call coarse_grain_variable_2d(area, variable(is:ie,js:je,k),&
            coarse_grained_variable(is_coarse:ie_coarse,js_coarse:je_coarse,k))
    enddo
    
  end subroutine coarse_grain_variable_3d
  
  subroutine fv_coarse_grained_diagnostics(Atm, Time, vort)

    type(fv_atmos_type), intent(in), target :: Atm(:)
    type(time_type), intent(in) :: Time
    real, intent(in) :: vort(:,:,:)
    logical :: used
    real, allocatable :: ps_coarse(:,:)
    real, allocatable :: vort_coarse(:,:,:)
    integer :: n = 1

    allocate(ps_coarse(is_coarse:ie_coarse, js_coarse:je_coarse),&
         vort_coarse(is_coarse:ie_coarse, js_coarse:je_coarse, 1:npz))
    
    call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), Atm(n)%ps(is:ie,js:je), ps_coarse)
    if( id_ps_coarse > 0 ) used = send_data(id_ps_coarse, ps_coarse, Time)
    deallocate(ps_coarse)

    call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), vort(is:ie,js:je,1:npz), vort_coarse)
    if( id_vort_coarse > 0 ) used = send_data(id_vort_coarse, vort_coarse, Time)
    deallocate(vort_coarse)
    
    if (id_coarse_lon > 0) used = send_data(id_coarse_lon, rad2deg *grid_coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse+1,1), Time)
    if (id_coarse_lat > 0) used = send_data(id_coarse_lat, rad2deg *grid_coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse+1,2), Time)
    
    if (id_coarse_lont > 0) used = send_data(id_coarse_lont, rad2deg * agrid_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1), Time)
    if (id_coarse_latt > 0) used = send_data(id_coarse_latt, rad2deg * agrid_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,2), Time)    

  end subroutine fv_coarse_grained_diagnostics

  ! This subroutine is copied from fms/horiz_interp/horiz_interp_test.F90;
  ! domain_decomp in fv_mp_mod.F90 does something similar, but it does a
  ! few other unnecessary things (and requires more arguments).
  subroutine define_cubic_mosaic(domain, ni, nj, layout)
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: layout(:)
    integer,        intent(in)    :: ni, nj
    integer   :: pe_start(6), pe_end(6)
    integer   :: global_indices(4,6), layout2d(2,6)
    integer, dimension(12)        :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(12)        :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact
    integer :: p, npes_per_tile, i

    ntiles = 6
    num_contact = 12
    p = 0
    npes_per_tile = mpp_npes()/ntiles
    do i = 1, 6
       layout2d(:,i) = layout(:)
       global_indices(1,i) = 1
       global_indices(2,i) = ni
       global_indices(3,i) = 1
       global_indices(4,i) = nj
       pe_start(i) = p
       p = p + npes_per_tile
       pe_end(i) = p-1
    enddo

    !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1;     tile2(1) = 2
    istart1(1) = ni;  iend1(1) = ni;  jstart1(1) = 1;      jend1(1) = nj
    istart2(1) = 1;   iend2(1) = 1;   jstart2(1) = 1;      jend2(1) = nj
    !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
    tile1(2) = 1; tile2(2) = 3
    istart1(2) = 1;      iend1(2) = ni;  jstart1(2) = nj;  jend1(2) = nj
    istart2(2) = 1;      iend2(2) = 1;   jstart2(2) = nj;  jend2(2) = 1
    !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
    tile1(3) = 1;     tile2(3) = 5
    istart1(3) = 1;   iend1(3) = 1;      jstart1(3) = 1;   jend1(3) = nj
    istart2(3) = ni;  iend2(3) = 1;      jstart2(3) = nj;  jend2(3) = nj
    !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
    tile1(4) = 1; tile2(4) = 6
    istart1(4) = 1;      iend1(4) = ni;  jstart1(4) = 1;   jend1(4) = 1
    istart2(4) = 1;      iend2(4) = ni;  jstart2(4) = nj;  jend2(4) = nj
    !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
    tile1(5) = 2;        tile2(5) = 3
    istart1(5) = 1;      iend1(5) = ni;  jstart1(5) = nj;  jend1(5) = nj
    istart2(5) = 1;      iend2(5) = ni;  jstart2(5) = 1;   jend2(5) = 1
    !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
    tile1(6) = 2; tile2(6) = 4
    istart1(6) = ni;  iend1(6) = ni;  jstart1(6) = 1;      jend1(6) = nj
    istart2(6) = ni;  iend2(6) = 1;   jstart2(6) = 1;      jend2(6) = 1
    !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
    tile1(7) = 2; tile2(7) = 6
    istart1(7) = 1;   iend1(7) = ni;  jstart1(7) = 1;   jend1(7) = 1
    istart2(7) = ni;  iend2(7) = ni;  jstart2(7) = nj;  jend2(7) = 1
    !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
    tile1(8) = 3; tile2(8) = 4
    istart1(8) = ni;  iend1(8) = ni;  jstart1(8) = 1;      jend1(8) = nj
    istart2(8) = 1;   iend2(8) = 1;   jstart2(8) = 1;      jend2(8) = nj
    !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
    tile1(9) = 3; tile2(9) = 5
    istart1(9) = 1;      iend1(9) = ni;  jstart1(9) = nj;  jend1(9) = nj
    istart2(9) = 1;      iend2(9) = 1;   jstart2(9) = nj;  jend2(9) = 1
    !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
    tile1(10) = 4; tile2(10) = 5
    istart1(10) = 1;     iend1(10) = ni; jstart1(10) = nj; jend1(10) = nj
    istart2(10) = 1;     iend2(10) = ni; jstart2(10) = 1;  jend2(10) = 1
    !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
    tile1(11) = 4; tile2(11) = 6
    istart1(11) = ni; iend1(11) = ni; jstart1(11) = 1;     jend1(11) = nj
    istart2(11) = ni; iend2(11) = 1;  jstart2(11) = 1;     jend2(11) = 1
    !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
    tile1(12) = 5; tile2(12) = 6
    istart1(12) = ni; iend1(12) = ni; jstart1(12) = 1;     jend1(12) = nj
    istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;     jend2(12) = nj
    call mpp_define_mosaic(global_indices, layout2d, domain, ntiles, num_contact, tile1, tile2, &
         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
         pe_start, pe_end, symmetry = .true., name = 'coarse cubic mosaic'  )

    return

  end subroutine define_cubic_mosaic
  
end module fv_coarse_grained_diagnostics_mod
