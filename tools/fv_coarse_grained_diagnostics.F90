module fv_coarse_grained_diagnostics_mod

  use diag_manager_mod, only: diag_axis_init, register_diag_field, send_data
  use fv_arrays_mod,   only: fv_atmos_type, fv_grid_bounds_type
  use fv_grid_tools_mod, only: init_grid
  use fv_mp_mod,       only: domain_decomp
  ! use fv_control_mod,   only: layout, npx
  use mpp_mod,         only: mpp_npes
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_define_mosaic, domain2d, mpp_define_io_domain
  use time_manager_mod, only: time_type

  implicit none
  private
  
  public :: fv_coarse_grained_diagnostics_init, fv_coarse_grained_diagnostics
  
  real, parameter:: missing_value = -1.e10
  integer :: id_ps_coarse
  integer :: coarsening_factor = 2  ! TODO: convert to namelist parameter
  integer :: layout(2) = (/2, 2/)  ! TODO: obtain from model
  integer :: io_layout(2) = (/1, 1/)  ! TODO: obtain from model
  ! integer :: npx = 49  ! TODO: obtain from model (use Atm(n)%gridstruct%npx_g?)
  integer :: is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse

  contains
  
  subroutine fv_coarse_grained_diagnostics_init(Atm, Time)
    type(time_type), intent(in) :: Time
    type(fv_atmos_type), intent(in) :: Atm(:)
    type(domain2d)     :: coarse_domain, coarse_coupler_domain
    type(fv_grid_bounds_type) :: bd
    integer :: coarse_axes(2)
    integer :: npx, target_resolution
    integer :: i, j, n, id_coarse_xt, id_coarse_yt, npx_target, npy_target
    integer :: npes_per_tile, num_contact, tile
    logical :: square_domain
    integer :: nregions = 6
    integer :: grid_type = 6
    logical :: nested = .false.
    
    integer, allocatable :: pelist(:)
    real, allocatable :: grid_xt(:), grid_yt(:)
    ! type(fv_atmos_type), allocatable, target :: Atm_coarse(:)
    
    n = 1
    
    ! TODO: add checks to make sure the coarsening factor is compatible with
    ! both the original resolution and the domain decomposition.
    ! TODO: check that grid_type is equal to 6 (this is not implemented for
    ! other kinds of grids).
    npx = Atm(n)%gridstruct%npx_g
    target_resolution = (npx - 1) / coarsening_factor
    npx_target = target_resolution + 1
    npy_target = target_resolution + 1

    allocate(grid_xt(target_resolution), grid_yt(target_resolution))
    grid_xt = (/ (i, i=1, target_resolution) /)
    grid_yt = (/ (j, j=1, target_resolution) /)
    
    ! Initialize values of is, ie, js, and je
    bd = Atm(n)%bd
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    
    ! Initialize value of coarse_domain and *[se]_coarse
    ! allocate(pelist(lbound(Atm(n)%pelist,1):ubound(Atm(n)%pelist,1)))
    ! pelist = Atm(n)%pelist
    ! tile = Atm(n)%tile_of_mosaic

    ! This might be the more canonical way of creating the domain; it
    ! automatically handles setting the io_domain appropriately
    ! call domain_decomp(npx_target, npy_target, nregions, grid_type, nested,&
    !     layout, io_layout, bd, tile, square_domain, npes_per_tile,&
    !     coarse_domain, coarse_coupler_domain, num_contact, pelist)
    ! deallocate(pelist)

    ! -------------------- Alternate code path for creating the cubic domain
    call define_cubic_mosaic(coarse_domain, target_resolution, target_resolution, layout)
    call mpp_define_io_domain(coarse_domain, io_layout)
    ! --------------------
    
    call mpp_get_compute_domain(coarse_domain, is_coarse, ie_coarse, &
         js_coarse, je_coarse)
    
    ! Create diagnostic axes
    ! TODO: Add static fields to contain the actual longitude and latitude
    ! coordinates.
    id_coarse_xt = diag_axis_init('grid_xt_coarse', grid_xt, &
         'degrees_E', 'x', 'T-cell longitude', set_name='coarse_grid', &
         Domain2=coarse_domain, tile_count=n)
    id_coarse_yt = diag_axis_init('grid_yt_coarse', grid_yt, &
         'degrees_N', 'y', 'T-cell latitude', set_name='coarse_grid', &
         Domain2=coarse_domain, tile_count=n)

    coarse_axes(1) = id_coarse_xt
    coarse_axes(2) = id_coarse_yt
    
    ! Register coarse diagnostics
    id_ps_coarse = register_diag_field('dynamics', 'ps_coarse', &
         coarse_axes(1:2), Time, 'surface pressure', 'Pa', &
         missing_value=missing_value)
    
    deallocate(grid_xt, grid_yt)
    
  end subroutine fv_coarse_grained_diagnostics_init

  subroutine coarse_grain_variable(area, variable, coarse_grained_variable)

    real, intent(in) :: area(is:ie,js:je)
    real, intent(in) :: variable(is:ie,js:je)
    real, intent(out) :: coarse_grained_variable(is_coarse:ie_coarse,js_coarse:je_coarse)
 
    real, allocatable :: area_weighted_variable(:,:)
    integer :: i, j, n, i_coarse, j_coarse, a

    ! write(*,*) 'FINES', lbound(variable,1), ubound(variable,1), 'COARSES',&
    !      & lbound(coarse_grained_variable,1), ubound(coarse_grained_variable&
    !      &,1), 'BOUNDS', is, ie
    ! write(*,*) 'FINES', lbound(variable,2), ubound(variable,2), 'COARSES',&
    !      & lbound(coarse_grained_variable,2), ubound(coarse_grained_variable&
    !      &,2), 'BOUNDS', js, je

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
    
  end subroutine coarse_grain_variable

  subroutine fv_coarse_grained_diagnostics(Atm, Time)

    type(fv_atmos_type), intent(in), target :: Atm(:)
    type(time_type), intent(in) :: Time
    logical :: used
    real, allocatable :: ps_coarse(:,:)
    integer :: n = 1

    allocate(ps_coarse(is_coarse:ie_coarse, js_coarse:je_coarse))
    ! write(*,*) 'FINE', lbound(Atm(n)%ps,1), ubound(Atm(n)%ps,1), 'COARSE', lbound(ps_coarse,1), ubound(ps_coarse,1)
    ! write(*,*) 'FINE', lbound(Atm(n)%ps,2), ubound(Atm(n)%ps,2), 'COARSE', lbound(ps_coarse,2), ubound(ps_coarse,2)
    call coarse_grain_variable(Atm(n)%gridstruct%area(is:ie,js:je), Atm(n)%ps(is:ie,js:je), ps_coarse)
    if( id_ps_coarse > 0 ) used = send_data(id_ps_coarse, ps_coarse, Time)
    deallocate(ps_coarse)
    
  end subroutine fv_coarse_grained_diagnostics

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
