!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of fvGFS.                                       *
!*                                                                     *
!* fvGFS is free software; you can redistribute it and/or modify it    *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!*                                                                     *
!* fvGFS is distributed in the hope that it will be useful, but        *
!* WITHOUT ANY WARRANTY; without even the implied warranty of          *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   *
!* General Public License for more details.                            *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
module fv_nesting_mod

   use mpp_domains_mod,     only: mpp_update_domains
   use mpp_domains_mod,     only: mpp_global_field
   use field_manager_mod,   only: MODEL_ATMOS
   use tracer_manager_mod,  only: get_tracer_index
   use fv_sg_mod,           only: neg_adj3
   use mpp_domains_mod,     only: mpp_get_data_domain, mpp_get_compute_domain, mpp_get_global_domain
   use mpp_domains_mod,     only: DGRID_NE, mpp_update_domains, domain2D
   use fv_restart_mod,      only: d2a_setup, d2c_setup
   use mpp_mod,             only: mpp_sync_self, mpp_sync, mpp_send, mpp_recv, mpp_error, FATAL, mpp_pe, WARNING, NOTE
   use mpp_domains_mod,     only: mpp_global_sum, BITWISE_EFP_SUM, BITWISE_EXACT_SUM
   use boundary_mod,        only: update_coarse_grid
   use boundary_mod,        only: nested_grid_BC_send, nested_grid_BC_recv, nested_grid_BC_save_proc
   use boundary_mod,        only: nested_grid_BC, nested_grid_BC_apply_intT
   use fv_mp_mod,           only: is, ie, js, je, isd, ied, jsd, jed, isc, iec, jsc, jec
   use fv_arrays_mod,       only: fv_grid_type, fv_flags_type, fv_atmos_type, fv_nest_type, fv_diag_type, fv_nest_BC_type_3D
   use fv_arrays_mod,       only: allocate_fv_nest_BC_type, fv_atmos_type, fv_grid_bounds_type, deallocate_fv_nest_BC_type
   use fv_grid_utils_mod,   only: ptop_min, g_sum, cubed_to_latlon, f_p
   use init_hydro_mod,      only: p_var
   use constants_mod,       only: grav, pi=>pi_8, radius, hlv, rdgas, cp_air, rvgas, cp_vapor, kappa
   use fv_mapz_mod,         only: mappm, remap_2d
   use fv_timing_mod,       only: timing_on, timing_off
   use fv_mp_mod,           only: is_master
   use fv_mp_mod,           only: mp_reduce_sum
   use fv_diagnostics_mod,  only: sphum_ll_fix, range_check
   use sw_core_mod,         only: divergence_corner, divergence_corner_nest

implicit none
   logical :: RF_initialized = .false.
   logical :: bad_range
   real, allocatable ::  rf(:), rw(:)
   integer :: kmax=1
   !Arrays for global grid total energy, used for grid nesting
   real, allocatable :: te_2d_coarse(:,:)
   real, allocatable :: dp1_coarse(:,:,:)

   !For nested grid buffers
   !Individual structures are allocated by nested_grid_BC_recv
   type(fv_nest_BC_type_3d) :: u_buf, v_buf, uc_buf, vc_buf, delp_buf, delz_buf, pt_buf, w_buf, divg_buf, pe_u_buf,pe_v_buf,pe_b_buf
   type(fv_nest_BC_type_3d), allocatable:: q_buf(:)
   !#ifdef USE_COND
   real, dimension(:,:,:), allocatable, target :: dum_West, dum_East, dum_North, dum_South
   !#endif

private
public :: twoway_nesting, setup_nested_grid_BCs, set_physics_BCs

contains

!!!!NOTE: Later we can add a flag to see if remap BCs are needed
!!!  if not we can save some code complexity and cycles by skipping it

 subroutine setup_nested_grid_BCs(npx, npy, npz, zvir, ncnst,     &
                        u, v, w, pt, delp, delz,q, uc, vc, &
#ifdef USE_COND
                        q_con, &
#ifdef MOIST_CAPPA
                        cappa, &
#endif
#endif
                        nested, inline_q, make_nh, ng, &
                        gridstruct, flagstruct, neststruct, &
                        nest_timestep, tracer_nest_timestep, &
                        domain, parent_grid, bd, nwat, ak, bk)

   
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(IN) :: zvir

    integer, intent(IN) :: npx, npy, npz
    integer, intent(IN) :: ncnst, ng, nwat
    logical, intent(IN) :: inline_q, make_nh,nested
    real, intent(IN), dimension(npz) :: ak, bk

    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u ! D grid zonal wind (m/s)
    real, intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   bd%isd:        ,bd%jsd:        ,1:)  !  W (m/s)
    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: delz(bd%isd:        ,bd%jsd:        ,1:)  ! height thickness (m)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) ! specific humidity and constituents
    real, intent(inout) :: uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) ! (uc,vc) mostly used as the C grid winds
    real, intent(inout) :: vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)
#ifdef USE_COND
    real, intent(inout) :: q_con(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
#ifdef MOIST_CAPPA
    real, intent(inout) :: cappa(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
#endif
#endif
    integer, intent(INOUT) :: nest_timestep, tracer_nest_timestep
    type(fv_atmos_type), intent(INOUT) :: parent_grid

    type(fv_grid_type), intent(INOUT) :: gridstruct
    type(fv_flags_type), intent(INOUT) :: flagstruct
    type(fv_nest_type), intent(INOUT), target :: neststruct
    type(domain2d), intent(INOUT) :: domain
    real :: divg(bd%isd:bd%ied+1,bd%jsd:bd%jed+1, npz)
    real :: ua(bd%isd:bd%ied,bd%jsd:bd%jed)
    real :: va(bd%isd:bd%ied,bd%jsd:bd%jed)
    real :: pe_ustag(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz+1)
    real :: pe_vstag(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz+1)
    real :: pe_bstag(bd%isd:bd%ied+1,bd%jsd:bd%jed+1,npz+1)
    real, parameter :: a13 = 1./3.

    integer :: i,j,k,n,p, sphum, npz_coarse
    logical :: do_pd

   type(fv_nest_BC_type_3d) :: delp_lag_BC, lag_BC, pe_lag_BC, pe_eul_BC
   type(fv_nest_BC_type_3d) :: lag_u_BC, pe_u_lag_BC, pe_u_eul_BC
   type(fv_nest_BC_type_3d) :: lag_v_BC, pe_v_lag_BC, pe_v_eul_BC
   type(fv_nest_BC_type_3d) :: lag_b_BC, pe_b_lag_BC, pe_b_eul_BC

    !local pointers
    logical, pointer :: child_grids(:)

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

    child_grids => neststruct%child_grids
    

    !IF nested, set up nested grid BCs for time-interpolation
    !(actually applying the BCs is done in dyn_core

    nest_timestep = 0
    if (.not. inline_q) tracer_nest_timestep = 0


    if (neststruct%nested .and. (.not. (neststruct%first_step) .or. make_nh) ) then
       do_pd = .true.
       call set_BCs_t0(ncnst, flagstruct%hydrostatic, neststruct) 
    else
       !On first timestep the t0 BCs are not initialized and may contain garbage
       do_pd = .false.
    end if

    !compute uc/vc for nested-grid BCs
    !!! CLEANUP: if we compute uc/vc here we don't need to do on the first call of c_sw, right?
    if (ANY(neststruct%child_grids)) then
       call timing_on('COMM_TOTAL')
       !!! CLEANUP: could we make this a non-blocking operation?
       !!! Is this needed? it is on the initialization step.
       call mpp_update_domains(delp, domain) !This is needed to make sure delp is updated for pe calculations
       call mpp_update_domains(u, v, &       
            domain, gridtype=DGRID_NE, complete=.true.)
       call timing_off('COMM_TOTAL')
!$OMP parallel do default(none) shared(isd,jsd,ied,jed,is,ie,js,je,npx,npy,npz, &
!$OMP       gridstruct,flagstruct,bd,u,v,uc,vc,nested,divg) &
!$OMP       private(ua,va)
       do k=1,npz
          call d2c_setup(u(isd,jsd,k),  v(isd,jsd,k),   &
               ua, va, &
               uc(isd,jsd,k), vc(isd,jsd,k), flagstruct%nord>0, &
               isd,ied,jsd,jed, is,ie,js,je, npx,npy, &
               gridstruct%grid_type, gridstruct%nested, &
               gridstruct%se_corner, gridstruct%sw_corner, &
               gridstruct%ne_corner, gridstruct%nw_corner, &
               gridstruct%rsin_u, gridstruct%rsin_v, &
               gridstruct%cosa_s, gridstruct%rsin2 )
          if (nested) then
             call divergence_corner_nest(u(isd,jsd,k), v(isd,jsd,k), ua, va, divg(isd,jsd,k), gridstruct, flagstruct, bd)
          else
             call divergence_corner(u(isd,jsd,k), v(isd,jsd,k), ua, va, divg(isd,jsd,k), gridstruct, flagstruct, bd)
          endif
       end do       
    endif

!! Nested grid: receive from parent grid (Lagrangian coordinate, npz_coarse)
    if (neststruct%nested) then

       npz_coarse = neststruct%parent_grid%npz

       if (.not. allocated(q_buf)) then
          allocate(q_buf(ncnst))
       endif

       call nested_grid_BC_recv(neststruct%nest_domain, 0, 0,  npz_coarse, bd, &
            delp_buf)
       do n=1,ncnst
          call nested_grid_BC_recv(neststruct%nest_domain, 0, 0, npz_coarse, bd, &
               q_buf(n))
       enddo
#ifndef SW_DYNAMICS
       call nested_grid_BC_recv(neststruct%nest_domain, 0, 0, npz_coarse, bd, &
            pt_buf)

       if (.not. flagstruct%hydrostatic) then
          call nested_grid_BC_recv(neststruct%nest_domain, 0, 0,  npz_coarse, bd, &
               w_buf)
          call nested_grid_BC_recv(neststruct%nest_domain, 0, 0,  npz_coarse, bd, &
               delz_buf)
       endif
#endif
       call nested_grid_BC_recv(neststruct%nest_domain, 0, 1,  npz_coarse+1, bd, &
            pe_u_buf)
       call nested_grid_BC_recv(neststruct%nest_domain, 1, 0,  npz_coarse+1, bd, &
            pe_v_buf)
       call nested_grid_BC_recv(neststruct%nest_domain, 1, 1,  npz_coarse+1, bd, &
            pe_b_buf)

       call nested_grid_BC_recv(neststruct%nest_domain, 0, 1,  npz_coarse, bd, &
            u_buf)
       call nested_grid_BC_recv(neststruct%nest_domain, 0, 1,  npz_coarse, bd, &
            vc_buf)
       call nested_grid_BC_recv(neststruct%nest_domain, 1, 0,  npz_coarse, bd, &
            v_buf)
       call nested_grid_BC_recv(neststruct%nest_domain, 1, 0,  npz_coarse, bd, &
            uc_buf)
       call nested_grid_BC_recv(neststruct%nest_domain, 1, 1,  npz_coarse, bd, &
            divg_buf)
    endif


!! Coarse grid: send to child grids (Lagrangian coordinate, npz_coarse)

    do p=1,size(child_grids)
       if (child_grids(p)) then
          call nested_grid_BC_send(delp, neststruct%nest_domain_all(p), 0, 0)
          do n=1,ncnst
             call nested_grid_BC_send(q(:,:,:,n), neststruct%nest_domain_all(p), 0, 0)
          enddo
#ifndef SW_DYNAMICS
          call nested_grid_BC_send(pt, neststruct%nest_domain_all(p), 0, 0)

          if (.not. flagstruct%hydrostatic) then
             call nested_grid_BC_send(w, neststruct%nest_domain_all(p), 0, 0)
             call nested_grid_BC_send(delz, neststruct%nest_domain_all(p), 0, 0)
          endif          
#endif
          !Compute and send staggered pressure
             !u points
             do j=js,je+1
             do i=is,ie
                pe_ustag(i,j,1) = ak(1)
             enddo
             enddo
             do k=1,npz
             do j=js,je+1
             do i=is,ie
                pe_ustag(i,j,k+1) = pe_ustag(i,j,k) + 0.5*(delp(i,j,k)+delp(i,j-1,k)) 
             enddo
             enddo
             enddo
             call nested_grid_BC_send(pe_ustag, neststruct%nest_domain_all(p), 0, 1)

             !v points
             do j=js,je
             do i=is,ie+1
                pe_vstag(i,j,1) = ak(1)
             enddo
             enddo
             do k=1,npz
             do j=js,je
             do i=is,ie+1
                pe_vstag(i,j,k+1) = pe_vstag(i,j,k) + 0.5*(delp(i,j,k)+delp(i-1,j,k)) 
             enddo
             enddo
             enddo
             call nested_grid_BC_send(pe_vstag, neststruct%nest_domain_all(p), 1, 0)

             !b points
             do j=js,je+1
             do i=is,ie+1
                pe_bstag(i,j,1) = ak(1)
             enddo
             enddo
             !Sets up so 3-point average is automatically done at the corner
             if (is == 1 .and. js == 1) then
                do k=1,npz
                   delp(0,0,k) = a13*(delp(1,1,k) + delp(0,1,k) + delp(1,0,k))
                enddo
             endif
             if (ie == npx-1 .and. js == 1) then
                do k=1,npz
                   delp(npx,0,k) = a13*(delp(npx-1,1,k) + delp(npx,1,k) + delp(npx-1,0,k))
                enddo
             endif
             if (is == 1 .and. je == npy-1) then
                do k=1,npz
                   delp(0,npy,k) = a13*(delp(1,npy-1,k) + delp(0,npy-1,k) + delp(1,npy,k))
                enddo
             endif
             if (ie == npx-1 .and. je == npy-1) then
                do k=1,npz
                   delp(npx,npy,k) = a13*(delp(npx-1,npy-1,k) + delp(npx,npy-1,k) + delp(npx-1,npy,k))
                enddo
             endif

             do k=1,npz
             do j=js,je+1
             do i=is,ie+1
                pe_bstag(i,j,k+1) = pe_bstag(i,j,k) + & 
                     0.25*(delp(i,j,k)+delp(i-1,j,k)+delp(i,j-1,k)+delp(i-1,j-1,k))
             enddo
             enddo
             enddo
             call nested_grid_BC_send(pe_bstag, neststruct%nest_domain_all(p), 1, 1)

          call nested_grid_BC_send(u, neststruct%nest_domain_all(p), 0, 1)
          call nested_grid_BC_send(vc, neststruct%nest_domain_all(p), 0, 1)
          call nested_grid_BC_send(v, neststruct%nest_domain_all(p), 1, 0)
          call nested_grid_BC_send(uc, neststruct%nest_domain_all(p), 1, 0)
          call nested_grid_BC_send(divg, neststruct%nest_domain_all(p), 1, 1)
       endif
    enddo
    
    !Nested grid: do computations
    ! Lag: coarse grid, npz_coarse, lagrangian coordinate---receive and use save_proc to copy into lag_BCs
    ! Eul: nested grid, npz, Eulerian (reference) coordinate
    ! Remapping from Lag to Eul
    if (nested) then
       call allocate_fv_nest_BC_type(delp_lag_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse,ng,0,0,0,.false.)
       call allocate_fv_nest_BC_type(lag_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse,ng,0,0,0,.false.)
       call allocate_fv_nest_BC_type(pe_lag_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse+1,ng,0,0,0,.false.)
       call allocate_fv_nest_BC_type(pe_eul_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz+1,ng,0,0,0,.false.)

       call nested_grid_BC_save_proc(neststruct%nest_domain, &
            neststruct%ind_h, neststruct%wt_h, 0, 0,  npx, npy, npz_coarse, bd, &
            delp_lag_BC, delp_buf, pd_in=do_pd)
       !The incoming delp is on the coarse grid's lagrangian coordinate. Re-create the reference coordinate
       call setup_eul_delp_BC(delp_lag_BC, neststruct%delp_BC, pe_lag_BC, pe_eul_BC, ak, bk, npx, npy, npz, npz_coarse, parent_grid%ptop, bd)

       do n=1,ncnst
          call nested_grid_BC_save_proc(neststruct%nest_domain, &
               neststruct%ind_h, neststruct%wt_h, 0, 0, npx,  npy,  npz_coarse, bd, &
               lag_BC, q_buf(n), pd_in=do_pd)
          call remap_BC(pe_lag_BC, pe_eul_BC, lag_BC, neststruct%q_BC(n), npx, npy, npz, npz_coarse, bd, 0, 0, 0, flagstruct%kord_tr)
       enddo
#ifndef SW_DYNAMICS
       call nested_grid_BC_save_proc(neststruct%nest_domain, &
            neststruct%ind_h, neststruct%wt_h, 0, 0, npx,  npy,  npz_coarse, bd, &
            lag_BC, pt_buf) 
       !NOTE: need to remap using peln, not pe
       call remap_BC(pe_lag_BC, pe_eul_BC, lag_BC, neststruct%pt_BC, npx, npy, npz, npz_coarse, bd, 0, 0, 1, abs(flagstruct%kord_tm), do_log_pe=.true.)

       sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
       if (flagstruct%hydrostatic) then
          call mpp_error(FATAL, " Hydrostatic nesting temporarily disabled during code refactoring")
          call setup_pt_BC(neststruct%pt_BC, neststruct%delp_BC, pe_eul_BC, neststruct%q_BC(sphum), npx, npy, npz, zvir, bd)
       else
          call nested_grid_BC_save_proc(neststruct%nest_domain, &
               neststruct%ind_h, neststruct%wt_h, 0, 0,  npx,  npy,  npz_coarse, bd, &
               lag_BC, w_buf)
          call remap_BC(pe_lag_BC, pe_eul_BC, lag_BC, neststruct%w_BC, npx, npy, npz, npz_coarse, bd, 0, 0, -1, flagstruct%kord_wz)
          call nested_grid_BC_save_proc(neststruct%nest_domain, &
               neststruct%ind_h, neststruct%wt_h, 0, 0,  npx,  npy, npz_coarse, bd, &
               lag_BC, delz_buf) !Need a negative-definite method? 
          call remap_delz_BC(pe_lag_BC, pe_eul_BC, delp_lag_BC, lag_BC, neststruct%delp_BC, neststruct%delz_BC, npx, npy, npz, npz_coarse, bd, 0, 0, 1, flagstruct%kord_wz)
          
          call setup_pt_NH_BC(neststruct%pt_BC, neststruct%delp_BC, neststruct%delz_BC, &
               neststruct%q_BC(sphum), neststruct%q_BC, ncnst, &
#ifdef USE_COND
               neststruct%q_con_BC, &
#ifdef MOIST_CAPPA
               neststruct%cappa_BC, &
#endif
#endif
               npx, npy, npz, zvir, bd)
       endif

#endif

       !!!NOTE: The following require remapping on STAGGERED grids, which requires additional pressure data

       call allocate_fv_nest_BC_type(pe_u_lag_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse+1,ng,0,0,1,.false.)
       call allocate_fv_nest_BC_type(pe_u_eul_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz+1       ,ng,0,0,1,.false.)
       call allocate_fv_nest_BC_type(lag_u_BC,   is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse  ,ng,0,0,1,.false.)
       call allocate_fv_nest_BC_type(pe_v_lag_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse+1,ng,0,1,0,.false.)
       call allocate_fv_nest_BC_type(pe_v_eul_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz+1       ,ng,0,1,0,.false.)
       call allocate_fv_nest_BC_type(lag_v_BC,   is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse  ,ng,0,1,0,.false.)
       call allocate_fv_nest_BC_type(pe_b_lag_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse+1,ng,0,1,1,.false.)
       call allocate_fv_nest_BC_type(pe_b_eul_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz+1       ,ng,0,1,1,.false.)
       call allocate_fv_nest_BC_type(lag_b_BC,   is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse  ,ng,0,1,1,.false.)

       call nested_grid_BC_save_proc(neststruct%nest_domain, &
            neststruct%ind_u, neststruct%wt_u, 0, 1,  npx,  npy,  npz_coarse+1, bd, &
            pe_u_lag_BC, pe_u_buf)
       call setup_eul_pe_BC(pe_u_lag_BC, pe_u_eul_BC, ak, bk, npx, npy, npz, npz_coarse, 0, 1, bd)
       call nested_grid_BC_save_proc(neststruct%nest_domain, &
            neststruct%ind_v, neststruct%wt_v, 1, 0,  npx,  npy,  npz_coarse+1, bd, &
            pe_v_lag_BC, pe_v_buf)
       call setup_eul_pe_BC(pe_v_lag_BC, pe_v_eul_BC, ak, bk, npx, npy, npz, npz_coarse, 1, 0, bd)
       call nested_grid_BC_save_proc(neststruct%nest_domain, &
            neststruct%ind_b, neststruct%wt_b, 1, 1,  npx,  npy,  npz_coarse+1, bd, &
            pe_b_lag_BC, pe_b_buf)
       call setup_eul_pe_BC(pe_b_lag_BC, pe_b_eul_BC, ak, bk, npx, npy, npz, npz_coarse, 1, 1, bd)

       call nested_grid_BC_save_proc(neststruct%nest_domain, &
            neststruct%ind_u, neststruct%wt_u, 0, 1,  npx,  npy,  npz_coarse, bd, &
            lag_u_BC, u_buf)
       call remap_BC(pe_u_lag_BC, pe_u_eul_BC, lag_u_BC, neststruct%u_BC, npx, npy, npz, npz_coarse, bd, 0, 1, -1, flagstruct%kord_mt)
       call nested_grid_BC_save_proc(neststruct%nest_domain, &
            neststruct%ind_u, neststruct%wt_u, 0, 1,  npx,  npy,  npz_coarse, bd, &
            lag_u_BC, vc_buf)
       call remap_BC(pe_u_lag_BC, pe_u_eul_BC, lag_u_BC, neststruct%vc_BC, npx, npy, npz, npz_coarse, bd, 0, 1, -1, flagstruct%kord_mt)
       call nested_grid_BC_save_proc(neststruct%nest_domain, &
            neststruct%ind_v, neststruct%wt_v, 1, 0,  npx,  npy,  npz_coarse, bd, &
            lag_v_BC, v_buf)
       call remap_BC(pe_v_lag_BC, pe_v_eul_BC, lag_v_BC, neststruct%v_BC, npx, npy, npz, npz_coarse, bd, 1, 0, -1, flagstruct%kord_mt)
       call nested_grid_BC_save_proc(neststruct%nest_domain, &
            neststruct%ind_v, neststruct%wt_v, 1, 0,  npx,  npy,  npz_coarse, bd, &
            lag_v_BC, uc_buf)
       call remap_BC(pe_v_lag_BC, pe_v_eul_BC, lag_v_BC, neststruct%uc_BC, npx, npy, npz, npz_coarse, bd, 1, 0, -1, flagstruct%kord_mt)
       call nested_grid_BC_save_proc(neststruct%nest_domain, &
            neststruct%ind_b, neststruct%wt_b, 1, 1,  npx,  npy,  npz_coarse, bd, &
            lag_b_BC, divg_buf)
       call remap_BC(pe_b_lag_BC, pe_b_eul_BC, lag_b_BC, neststruct%divg_BC, npx, npy, npz, npz_coarse, bd, 1, 1, -1, flagstruct%kord_mt)

       call deallocate_fv_nest_BC_type(delp_lag_BC)
       call deallocate_fv_nest_BC_type(lag_BC)
       call deallocate_fv_nest_BC_type(pe_lag_BC)
       call deallocate_fv_nest_BC_type(pe_eul_BC)

       call deallocate_fv_nest_BC_type(pe_u_lag_BC)
       call deallocate_fv_nest_BC_type(pe_u_eul_BC)
       call deallocate_fv_nest_BC_type(lag_u_BC)
       call deallocate_fv_nest_BC_type(pe_v_lag_BC)
       call deallocate_fv_nest_BC_type(pe_v_eul_BC)
       call deallocate_fv_nest_BC_type(lag_v_BC)
       call deallocate_fv_nest_BC_type(pe_b_lag_BC)
       call deallocate_fv_nest_BC_type(pe_b_eul_BC)
       call deallocate_fv_nest_BC_type(lag_b_BC)

       !Correct halo values have now been set up for BCs; we can go ahead and apply them too
       call nested_grid_BC_apply_intT(delp, &
            0, 0, npx, npy, npz, bd, 1., 1., &
            neststruct%delp_BC, bctype=neststruct%nestbctype  )
       do n=1,ncnst
          call nested_grid_BC_apply_intT(q(:,:,:,n), &
               0, 0, npx, npy, npz, bd, 1., 1., &
               neststruct%q_BC(n), bctype=neststruct%nestbctype  )          
       enddo
#ifndef SW_DYNAMICS
       call nested_grid_BC_apply_intT(pt, &
            0, 0, npx, npy, npz, bd, 1., 1., &
            neststruct%pt_BC, bctype=neststruct%nestbctype  )
       if (.not. flagstruct%hydrostatic) then
          call nested_grid_BC_apply_intT(w, &
               0, 0, npx, npy, npz, bd, 1., 1., &
               neststruct%w_BC, bctype=neststruct%nestbctype  )
          call nested_grid_BC_apply_intT(delz, &
               0, 0, npx, npy, npz, bd, 1., 1., &
               neststruct%delz_BC, bctype=neststruct%nestbctype  )
       endif
#ifdef USE_COND
       call nested_grid_BC_apply_intT(q_con, &
            0, 0, npx, npy, npz, bd, 1., 1., &
            neststruct%q_con_BC, bctype=neststruct%nestbctype  )            
#ifdef MOIST_CAPPA
       call nested_grid_BC_apply_intT(cappa, &
            0, 0, npx, npy, npz, bd, 1., 1., &
            neststruct%cappa_BC, bctype=neststruct%nestbctype  )            
#endif
#endif
#endif
       call nested_grid_BC_apply_intT(u, &
            0, 1, npx, npy, npz, bd, 1., 1., &
            neststruct%u_BC, bctype=neststruct%nestbctype  )            
       call nested_grid_BC_apply_intT(vc, &
            0, 1, npx, npy, npz, bd, 1., 1., &
            neststruct%vc_BC, bctype=neststruct%nestbctype  )            
       call nested_grid_BC_apply_intT(v, &
            1, 0, npx, npy, npz, bd, 1., 1., &
            neststruct%v_BC, bctype=neststruct%nestbctype  )            
       call nested_grid_BC_apply_intT(uc, &
            1, 0, npx, npy, npz, bd, 1., 1., &
            neststruct%uc_BC, bctype=neststruct%nestbctype  )            
       !!!NOTE: Divg not available here but not needed
       !!! until dyn_core anyway.
!!$       call nested_grid_BC_apply_intT(divg, &
!!$            1, 1, npx, npy, npz, bd, 1., 1., &
!!$            neststruct%divg_BC, bctype=neststruct%nestbctype  )            

       !Update domains needed for Rayleigh damping
       call mpp_update_domains(u, v, domain, gridtype=DGRID_NE)
       call mpp_update_domains(w, domain, complete=.true.) 

    endif

    if (neststruct%first_step) then
       if (neststruct%nested) call set_BCs_t0(ncnst, flagstruct%hydrostatic, neststruct)
       neststruct%first_step = .false.
       if (.not. flagstruct%hydrostatic) flagstruct%make_nh= .false. 
    else if (flagstruct%make_nh) then
       if (neststruct%nested) call set_NH_BCs_t0(neststruct)
       flagstruct%make_nh= .false. 
    endif

    !Unnecessary?
!!$    if ( neststruct%nested .and. .not. neststruct%divg_BC%initialized) then
!!$       neststruct%divg_BC%east_t0  = neststruct%divg_BC%east_t1
!!$       neststruct%divg_BC%west_t0  = neststruct%divg_BC%west_t1
!!$       neststruct%divg_BC%north_t0 = neststruct%divg_BC%north_t1
!!$       neststruct%divg_BC%south_t0 = neststruct%divg_BC%south_t1 
!!$       neststruct%divg_BC%initialized = .true.
!!$    endif


    call mpp_sync_self

 end subroutine setup_nested_grid_BCs

 subroutine set_physics_BCs(ps, u_dt, v_dt, flagstruct, gridstruct, neststruct, npx, npy, npz, ng, ak, bk, bd)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_flags_type), intent(IN) :: flagstruct
   type(fv_nest_type), intent(INOUT), target :: neststruct
   type(fv_grid_type) :: gridstruct
   integer, intent(IN) :: npx, npy, npz, ng
   real, intent(IN), dimension(npz+1) :: ak, bk
   real, intent(INOUT), dimension(bd%isd:bd%ied,bd%jsd:bd%jed) :: ps
   real, intent(INOUT), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz) :: u_dt, v_dt
   real, dimension(1,1) :: parent_ps ! dummy variable for nesting
   type(fv_nest_BC_type_3d) :: u_dt_buf, v_dt_buf, pe_src_BC, pe_dst_BC!, var_BC

   integer :: n, npz_coarse
   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed

   if (gridstruct%nested) then

       npz_coarse = neststruct%parent_grid%npz

      !Both nested and coarse grids assumed on Eulerian coordinates at this point
      !Only need to fetch ps to form pressure levels
      !Note also u_dt and v_dt are unstaggered
      call nested_grid_BC(ps, parent_ps, neststruct%nest_domain, neststruct%ind_h, neststruct%wt_h, 0, 0, &
           npx, npy, bd, 1, npx-1, 1, npy-1)
      call nested_grid_BC_recv(neststruct%nest_domain, 0, 0, npz_coarse, bd, u_dt_buf) 
      call nested_grid_BC_recv(neststruct%nest_domain, 0, 0, npz_coarse, bd, v_dt_buf) 

      call allocate_fv_nest_BC_type(pe_src_BC, is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse+1,ng,0,0,0,.false.)
      call allocate_fv_nest_BC_type(pe_dst_BC, is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz+1,ng,0,0,0,.false.)

      call copy_ps_BC(ps, pe_src_BC, npx, npy, npz_coarse, 0, 0, bd)
      call setup_eul_pe_BC(pe_src_BC, pe_dst_BC, ak, bk, npx, npy, npz, npz_coarse, 0, 0, bd, &
           make_src_in=.true., ak_src=neststruct%parent_grid%ak, bk_src=neststruct%parent_grid%bk)

      !Note that iv=-1 is used for remapping winds, which sets the lower reconstructed values to 0 if
      ! there is a 2dx signal. Is this the best for **tendencies** though?? Probably not---so iv=1 here
      call set_BC_direct( pe_src_BC, pe_dst_BC, u_dt_buf, u_dt, neststruct, npx, npy, npz, npz_coarse, ng, bd, 0, 0, 1, flagstruct%kord_mt)
      call set_BC_direct( pe_src_BC, pe_dst_BC, v_dt_buf, v_dt, neststruct, npx, npy, npz, npz_coarse, ng, bd, 0, 0, 1, flagstruct%kord_mt)

      call deallocate_fv_nest_BC_type(pe_src_BC)
      call deallocate_fv_nest_BC_type(pe_dst_BC)

   endif
   do n=1,size(neststruct%child_grids)
      if (neststruct%child_grids(n)) then
         call nested_grid_BC(ps, neststruct%nest_domain_all(n), 0, 0)
         call nested_grid_BC_send(u_dt, neststruct%nest_domain_all(n), 0, 0)
         call nested_grid_BC_send(v_dt, neststruct%nest_domain_all(n), 0, 0)
      endif
   enddo


 end subroutine set_physics_BCs

 subroutine set_BC_direct( pe_src_BC, pe_dst_BC, buf, var, neststruct, npx, npy, npz, npz_coarse, ng, bd, istag, jstag, iv, kord)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_nest_type), intent(INOUT) :: neststruct
   integer, intent(IN) :: npx, npy, npz, npz_coarse, ng, istag, jstag, iv, kord
   real, intent(INOUT), dimension(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,npz) :: var
   type(fv_nest_BC_type_3d), intent(INOUT) :: buf, pe_src_BC, pe_dst_BC
   type(fv_nest_BC_type_3d) :: var_BC
   

   call allocate_fv_nest_BC_type(var_BC,is,ie,js,je,isd,ied,jsd,jed,npx,npy,npz_coarse,ng,0,istag,jstag,.false.)

   call nested_grid_BC_save_proc(neststruct%nest_domain, neststruct%ind_h, neststruct%wt_h, istag, jstag, &
        npx, npy, npz_coarse, bd, var_BC, buf)
   call remap_BC_direct(pe_src_BC, pe_dst_BC, var_BC, var, npx, npy, npz, npz_coarse, bd, istag, jstag, iv, kord)
   
   call deallocate_fv_nest_BC_type(var_BC)


 end subroutine set_BC_direct

 subroutine setup_pt_BC(pt_BC, delp_BC, pe_eul_BC, sphum_BC, npx, npy, npz, zvir, bd)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_nest_BC_type_3d), intent(IN), target    :: delp_BC, pe_eul_BC, sphum_BC
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pt_BC
   integer, intent(IN) :: npx, npy, npz
   real, intent(IN) :: zvir

   real, dimension(:,:,:), pointer :: ptBC, sphumBC

   real :: peln, pkz

   integer :: i,j,k, istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   if (is == 1) then
      ptBC    =>    pt_BC%west_t1
      sphumBC => sphum_BC%west_t1
!$OMP parallel do default(none) shared(npz,isd,jsd,jed,ptBC,zvir,sphumBC) private(pkz)
      do k=1,npz
      do j=jsd,jed
      do i=isd,0
         ptBC(i,j,k) = ptBC(i,j,k)/pkz*(1.+zvir*sphumBC(i,j,k))
      end do
      end do
      end do
   end if

   if (js == 1) then
      ptBC    =>    pt_BC%south_t1
      sphumBC => sphum_BC%south_t1
      if (is == 1) then
         istart = is
      else
         istart = isd
      end if
      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

!$OMP parallel do default(none) shared(npz,jsd,istart,iend,ptBC,zvir,sphumBC) private(pkz)
      do k=1,npz
      do j=jsd,0
      do i=istart,iend
         ptBC(i,j,k) = ptBC(i,j,k)/pkz * &
              (1.+zvir*sphumBC(i,j,k))
      end do
      end do
      end do
   end if


   if (ie == npx-1) then
      ptBC    =>    pt_BC%east_t1
      sphumBC => sphum_BC%east_t1
!$OMP parallel do default(none) shared(npz,jsd,jed,npx,ied,ptBC,zvir,sphumBC) private(pkz)
      do k=1,npz
      do j=jsd,jed
      do i=npx,ied
         ptBC(i,j,k) = ptBC(i,j,k)/pkz * &
              (1.+zvir*sphumBC(i,j,k))
      end do
      end do
      end do
   end if

   if (je == npy-1) then
      ptBC    =>    pt_BC%north_t1
      sphumBC => sphum_BC%north_t1
      if (is == 1) then
         istart = is
      else
         istart = isd
      end if
      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

!$OMP parallel do default(none) shared(npz,npy,jed,npx,istart,iend,ptBC,zvir,sphumBC) private(pkz)
      do k=1,npz
      do j=npy,jed
      do i=istart,iend
         ptBC(i,j,k) = ptBC(i,j,k)/pkz * &
              (1.+zvir*sphumBC(i,j,k))
      end do
      end do
      end do
   end if
   
 end subroutine setup_pt_BC


 subroutine setup_pt_BC_k


 end subroutine setup_pt_BC_k

 subroutine setup_eul_delp_BC(delp_lag_BC, delp_eul_BC, pe_lag_BC, pe_eul_BC, ak_dst, bk_dst, npx, npy, npz, npz_coarse, ptop_src, bd)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_nest_BC_type_3d), intent(INOUT), target :: delp_lag_BC
   type(fv_nest_BC_type_3d), intent(INOUT), target :: delp_eul_BC, pe_lag_BC, pe_eul_BC
   integer, intent(IN) :: npx, npy, npz, npz_coarse
   real, intent(IN), dimension(npz+1) :: ak_dst, bk_dst
   real, intent(IN) :: ptop_src

   integer :: i,j,k, istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   if (is == 1) then
      call setup_eul_delp_BC_k(delp_lag_BC%west_t1, delp_eul_BC%west_t1, pe_lag_BC%west_t1, pe_eul_BC%west_t1, &
           ptop_src, ak_dst, bk_dst, isd, 0, isd, 0, jsd, jed, npz, npz_coarse)
   end if

   if (ie == npx-1) then
      call setup_eul_delp_BC_k(delp_lag_BC%east_t1, delp_eul_BC%east_t1, pe_lag_BC%east_t1, pe_eul_BC%east_t1, &
           ptop_src, ak_dst, bk_dst, npx, ied, npx, ied, jsd, jed, npz, npz_coarse)
   end if

   if (is == 1) then
      istart = is
   else
      istart = isd
   end if
   if (ie == npx-1) then
      iend = ie
   else
      iend = ied
   end if

   if (js == 1) then
      call setup_eul_delp_BC_k(delp_lag_BC%south_t1, delp_eul_BC%south_t1, pe_lag_BC%south_t1, pe_eul_BC%south_t1, &
           ptop_src, ak_dst, bk_dst, isd, ied, istart, iend, jsd, 0, npz, npz_coarse)
   end if

   if (je == npy-1) then
      call setup_eul_delp_BC_k(delp_lag_BC%north_t1, delp_eul_BC%north_t1, pe_lag_BC%north_t1, pe_eul_BC%north_t1, &
           ptop_src, ak_dst, bk_dst, isd, ied, istart, iend, npy, jed, npz, npz_coarse)
   end if
   
 end subroutine setup_eul_delp_BC

 subroutine setup_eul_delp_BC_k(delplagBC, delpeulBC, pelagBC, peeulBC, ptop_src, ak_dst, bk_dst, isd, ied, istart, iend, jstart, jend, npz, npz_coarse)

   integer, intent(IN) :: isd, ied, istart, iend, jstart, jend, npz, npz_coarse
   real, intent(INOUT) :: delplagBC(isd:ied,jstart:jend,npz_coarse), pelagBC(isd:ied,jstart:jend,npz_coarse+1)
   real, intent(INOUT) :: delpeulBC(isd:ied,jstart:jend,npz), peeulBC(isd:ied,jstart:jend,npz+1)
   real, intent(IN) :: ptop_src, ak_dst(npz+1), bk_dst(npz+1)

   integer :: i,j,k

   character(len=120) :: errstring

!$no-OMP parallel do default(none) shared(istart,iend,jstart,jend,psBC,ptop_src)
   do j=jstart,jend
   do i=istart,iend
      pelagBC(i,j,1) = ptop_src 
   enddo
   enddo
!$no-OMP parallel do default(none) shared(istart,iend,jstart,jend,npz,psBC,delplagBC)
   do k=1,npz_coarse
   do j=jstart,jend
   do i=istart,iend
      pelagBC(i,j,k+1) = pelagBC(i,j,k) + delplagBC(i,j,k)
   end do
   end do
   end do
!$no-OMP parallel do default(none) shared(istart,iend,jstart,jend,npz,psBC,delpeulBC,ak_dst,bk_dst)
   do k=1,npz+1
   do j=jstart,jend
   do i=istart,iend
      peeulBC(i,j,k) = ak_dst(k) + pelagBC(i,j,npz_coarse+1)*bk_dst(k)
   enddo
   enddo
   enddo
   do k=1,npz
   do j=jstart,jend
   do i=istart,iend
      delpeulBC(i,j,k) = peeulBC(i,j,k+1) - peeulBC(i,j,k)
   enddo
   enddo
   enddo

!!$!!! DEBUG CODE
!!$   !If more than a few percent difference then log the error
!!$   do k=1,npz
!!$   do j=jstart,jend
!!$   do i=istart,iend
!!$      if (delpeulBC(i,j,k) <= 0.) then
!!$         write(errstring,'(3I5, 3(2x, G))'), i, j, k, pelagBC(i,j,k), peeulBC(i,j,k)
!!$         call mpp_error(WARNING, ' Invalid pressure BC at '//errstring)
!!$      else if (abs( peeulBC(i,j,k) - pelagBC(i,j,k)) > 100.0 ) then
!!$         write(errstring,'(3I5, 3(2x, G))'), i, j, k, pelagBC(i,j,k), peeulBC(i,j,k)
!!$         call mpp_error(WARNING, ' Remap BC: pressure deviation at '//errstring)
!!$      endif
!!$   enddo
!!$   enddo
!!$   enddo
!!$!!! END DEBUG CODE

 end subroutine setup_eul_delp_BC_k

 subroutine copy_ps_BC(ps, pe_BC, npx, npy, npz, istag, jstag, bd)

   integer, intent(IN) :: npx, npy, npz, istag, jstag
   type(fv_grid_bounds_type), intent(IN) :: bd
   real, intent(IN) :: ps(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag)
   type(fv_nest_BC_type_3d), intent(INOUT) :: pe_BC

   integer :: i,j,k, istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   if (is == 1) then
      do j=jsd,jed+jstag
      do i=isd,0
         pe_BC%west_t1(i,j,npz+1) = ps(i,j)
      enddo
      enddo
   end if

   if (ie == npx-1) then
      do j=jsd,jed+jstag
      do i=npx+istag,ied+istag
         pe_BC%east_t1(i,j,npz+1) = ps(i,j)
      enddo
      enddo
   end if

   if (is == 1) then
      istart = is
   else
      istart = isd
   end if
   if (ie == npx-1) then
      iend = ie
   else
      iend = ied
   end if

   if (js == 1) then
      do j=jsd,0
      do i=isd,ied+istag
         pe_BC%south_t1(i,j,npz+1) = ps(i,j)
      enddo
      enddo
   end if

   if (je == npy-1) then
      do j=npy+jstag,jed+jstag
      do i=isd,ied+istag
         pe_BC%north_t1(i,j,npz+1) = ps(i,j)
      enddo
      enddo
   end if

 end subroutine copy_ps_BC

 subroutine setup_eul_pe_BC(pe_src_BC, pe_eul_BC, ak_dst, bk_dst, npx, npy, npz, npz_src, istag, jstag, bd, make_src_in, ak_src, bk_src)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pe_src_BC, pe_eul_BC
   integer, intent(IN) :: npx, npy, npz, npz_src, istag, jstag
   real, intent(IN), dimension(npz+1) :: ak_dst, bk_dst
   logical, intent(IN), OPTIONAL :: make_src_in
   real, intent(IN), OPTIONAL :: ak_src(npz_src), bk_src(npz_src)

   logical :: make_src

   integer :: i,j,k, istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   make_src = .false.
   if (present(make_src_in)) make_src = make_src_in

   if (is == 1) then
      call setup_eul_pe_BC_k(pe_src_BC%west_t1, pe_eul_BC%west_t1, ak_dst, bk_dst, isd, 0, isd, 0, jsd, jed+jstag, npz, npz_src, &
           make_src, ak_src, bk_src)
   end if

   if (ie == npx-1) then
      call setup_eul_pe_BC_k(pe_src_BC%east_t1, pe_eul_BC%east_t1, ak_dst, bk_dst, npx+istag, ied+istag, npx+istag, ied+istag, jsd, jed+jstag, npz, npz_src, &
           make_src, ak_src, bk_src)
   end if

   if (is == 1) then
      istart = is
   else
      istart = isd
   end if
   if (ie == npx-1) then
      iend = ie
   else
      iend = ied
   end if

   if (js == 1) then
      call setup_eul_pe_BC_k(pe_src_BC%south_t1, pe_eul_BC%south_t1, ak_dst, bk_dst, isd, ied+istag, istart, iend+istag, jsd, 0, npz, npz_src, &
           make_src, ak_src, bk_src)
!!$         !!! DEBUG CODE
!!$         i = 76
!!$         j = 0
!!$         if ( istart <= i .and. iend >= i ) then
!!$            do k=1,npz+1
!!$               write(mpp_pe()+2000,*) k, pe_src_BC%south_t1(i,j,k), pe_eul_BC%south_t1(i,j,k)
!!$            enddo
!!$            write(mpp_pe()+2000,*) 
!!$         endif
!!$         !!! END DEBUG CODE
   end if

   if (je == npy-1) then
      call setup_eul_pe_BC_k(pe_src_BC%north_t1, pe_eul_BC%north_t1, ak_dst, bk_dst, isd, ied+istag, istart, iend+istag, npy+jstag, jed+jstag, npz, npz_src, &
           make_src, ak_src, bk_src)
   end if
   
 end subroutine setup_eul_pe_BC

 subroutine setup_eul_pe_BC_k(pesrcBC, peeulBC, ak_dst, bk_dst, isd, ied, istart, iend, jstart, jend, npz, npz_src, make_src, ak_src, bk_src)

   integer, intent(IN) :: isd, ied, istart, iend, jstart, jend, npz, npz_src
   real, intent(INOUT) :: pesrcBC(isd:ied,jstart:jend,npz_src+1)
   real, intent(INOUT) :: peeulBC(isd:ied,jstart:jend,npz+1)
   real, intent(IN) :: ak_dst(npz+1), bk_dst(npz+1)
   logical, intent(IN) :: make_src
   real, intent(IN) :: ak_src(npz_src+1), bk_src(npz_src+1)

   integer :: i,j,k

   character(len=120) :: errstring

   do k=1,npz+1
   do j=jstart,jend
   do i=istart,iend
      peeulBC(i,j,k) = ak_dst(k) + pesrcBC(i,j,npz_src+1)*bk_dst(k)
   enddo
   enddo
   enddo
   
   if (make_src) then
   do k=1,npz_src+1
   do j=jstart,jend
   do i=istart,iend
      pesrcBC(i,j,k) = ak_src(k) + pesrcBC(i,j,npz_src+1)*bk_src(k)
   enddo
   enddo
   enddo
   endif


 end subroutine setup_eul_pe_BC_k

 subroutine remap_BC(pe_lag_BC, pe_eul_BC, var_lag_BC, var_eul_BC, npx, npy, npz, npz_coarse, bd, istag, jstag, iv, kord, do_log_pe)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pe_lag_BC, var_lag_BC
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pe_eul_BC, var_eul_BC
   integer, intent(IN) :: npx, npy, npz, npz_coarse, istag, jstag, iv, kord
   logical, intent(IN), OPTIONAL :: do_log_pe

   logical :: log_pe = .false.

   integer :: i,j,k, istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   if (present(do_log_pe)) log_pe = do_log_pe
   
   if (is == 1) then
      call remap_BC_k(pe_lag_BC%west_t1, pe_eul_BC%west_t1, var_lag_BC%west_t1, var_eul_BC%west_t1, isd, 0, isd, 0, jsd, jed+jstag, npz, npz_coarse, iv, kord, log_pe)
   end if

   if (ie == npx-1) then
      call remap_BC_k(pe_lag_BC%east_t1, pe_eul_BC%east_t1, var_lag_BC%east_t1, var_eul_BC%east_t1, npx+istag, ied+istag, npx+istag, ied+istag, jsd, jed+jstag, npz, npz_coarse, iv, kord, log_pe)
   end if

   if (is == 1) then
      istart = is
   else
      istart = isd
   end if
   if (ie == npx-1) then
      iend = ie
   else
      iend = ied
   end if

   if (js == 1) then
      call remap_BC_k(pe_lag_BC%south_t1, pe_eul_BC%south_t1, var_lag_BC%south_t1, var_eul_BC%south_t1, isd, ied+istag, istart, iend+istag, jsd, 0, npz, npz_coarse, iv, kord, log_pe)
!!$         !!! DEBUG CODE
!!$         i = 76
!!$         j = 0
!!$         if ( istart <= i .and. iend >= i ) then
!!$            do k=1,npz
!!$               write(mpp_pe()+1000,*) k, pe_lag_BC%south_t1(i,j,k), pe_eul_BC%south_t1(i,j,k), var_lag_BC%south_t1(i,j,k), var_eul_BC%south_t1(i,j,k)
!!$            enddo
!!$            write(mpp_pe()+1000,*) 
!!$         endif
!!$         !!! END DEBUG CODE
   end if

   if (je == npy-1) then
      call remap_BC_k(pe_lag_BC%north_t1, pe_eul_BC%north_t1, var_lag_BC%north_t1, var_eul_BC%north_t1, isd, ied+istag, istart, iend+istag, npy+jstag, jed+jstag, npz, npz_coarse, iv, kord, log_pe)
   end if
   
 end subroutine remap_BC

 subroutine remap_BC_direct(pe_lag_BC, pe_eul_BC, var_lag_BC, var, npx, npy, npz, npz_coarse, bd, istag, jstag, iv, kord, do_log_pe)

   type(fv_grid_bounds_type), intent(IN) :: bd
   integer, intent(IN) :: npx, npy, npz, npz_coarse, istag, jstag, iv, kord
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pe_lag_BC, var_lag_BC
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pe_eul_BC
   real, intent(INOUT) ::  var(bd%isd:bd%ied+istag,bd%jsd:bd%jed+jstag,npz)
   logical, intent(IN), OPTIONAL :: do_log_pe

   logical :: log_pe = .false.

   integer :: i,j,k, istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   if (present(do_log_pe)) log_pe = do_log_pe
   
   if (is == 1) then
      !I was unable how to do pass-by-memory referencing on parts of the 3D var array,
      ! so instead I am doing an inefficient copy and copy-back. --- lmh 14jun17
      call remap_BC_k(pe_lag_BC%west_t1, pe_eul_BC%west_t1, var_lag_BC%west_t1, var(isd:0,jsd:jed+jstag,:), isd, 0, isd, 0, jsd, jed+jstag, npz, npz_coarse, iv, kord, log_pe)
!!$!!! DEBUG CODE
!!$      write(mpp_pe()+1000,*) j, pe_lag_BC%west_t1(isd,jsd,npz+1), pe_eul_BC%west_t1(isd,jsd,npz+1), var_lag_BC%west_t1(isd,jsd,npz-2), var(isd,jsd,npz-2)
!!$!!! END DEBUG CODE
   end if

   if (ie == npx-1) then
!      var(npx+istag:ied+istag,jsd:jed+jstag,:) = -999.
      call remap_BC_k(pe_lag_BC%east_t1, pe_eul_BC%east_t1, var_lag_BC%east_t1, var(npx+istag:ied+istag,jsd:jed+jstag,:), npx+istag, ied+istag, npx+istag, ied+istag, jsd, jed+jstag, npz, npz_coarse, iv, kord, log_pe)
   end if

   if (is == 1) then
      istart = is
   else
      istart = isd
   end if
   if (ie == npx-1) then
      iend = ie
   else
      iend = ied
   end if

   if (js == 1) then
!      var(isd:ied+istag,jsd:0,:) = -999.
      call remap_BC_k(pe_lag_BC%south_t1, pe_eul_BC%south_t1, var_lag_BC%south_t1, var(isd:ied+istag,jsd:0,:), isd, ied+istag, istart, iend+istag, jsd, 0, npz, npz_coarse, iv, kord, log_pe)
   end if

   if (je == npy-1) then
!      var(isd:ied+istag,npy+jstag:jed+jstag,:) = -999.
      call remap_BC_k(pe_lag_BC%north_t1, pe_eul_BC%north_t1, var_lag_BC%north_t1, var(isd:ied+istag,npy+jstag:jed+jstag,:), isd, ied+istag, istart, iend+istag, npy+jstag, jed+jstag, npz, npz_coarse, iv, kord, log_pe)
   end if
   
 end subroutine remap_BC_direct

 subroutine remap_BC_k(pe_lagBC, pe_eulBC, var_lagBC, var_eulBC, isd, ied, istart, iend, jstart, jend, npz, npz_coarse, iv, kord, log_pe)

   integer, intent(IN) :: isd, ied, istart, iend, jstart, jend, npz, npz_coarse, iv, kord
   logical, intent(IN) :: log_pe
   real, intent(INOUT) :: pe_lagBC(isd:ied,jstart:jend,npz_coarse+1), var_lagBC(isd:ied,jstart:jend,npz_coarse)
   real, intent(INOUT) :: pe_eulBC(isd:ied,jstart:jend,npz+1), var_eulBC(isd:ied,jstart:jend,npz)

   integer :: i, j, k
   real peln_lag(istart:iend,npz_coarse+1)
   real peln_eul(istart:iend,npz+1)
   character(120) :: errstring
   
   do j=jstart,jend

      if (log_pe) then

         do k=1,npz_coarse+1
         do i=istart,iend
!!$!!! DEBUG CODE
!!$            if (pe_lagBC(i,j,k) <= 0.) then
!!$               write(errstring,'(3I5, 2x, G)'), i, j, k, pe_lagBC(i,j,k)
!!$               call mpp_error(WARNING, ' Remap BC: invalid pressure at at '//errstring)               
!!$            endif
!!$!!! END DEBUG CODE
            peln_lag(i,k) = log(pe_lagBC(i,j,k))
         enddo
         enddo

         do k=1,npz+1
         do i=istart,iend
!!$!!! DEBUG CODE
!!$            if (pe_lagBC(i,j,k) <= 0.) then
!!$               write(errstring,'(3I5, 2x, G)'), i, j, k, pe_lagBC(i,j,k)
!!$               call mpp_error(WARNING, ' Remap BC: invalid pressure at at '//errstring)               
!!$            endif
!!$!!! END DEBUG CODE
            peln_eul(i,k) = log(pe_eulBC(i,j,k))
         enddo
         enddo

         call remap_2d(npz_coarse, peln_lag, var_lagBC(istart:iend,j:j,:), &
                       npz, peln_eul, var_eulBC(istart:iend,j:j,:), &
                       istart, iend, iv, kord)

      else

         call remap_2d(npz_coarse, pe_lagBC(istart:iend,j:j,:), var_lagBC(istart:iend,j:j,:), &
                       npz, pe_eulBC(istart:iend,j:j,:), var_eulBC(istart:iend,j:j,:), &
                       istart, iend, iv, kord)

      endif

   enddo

 end subroutine remap_BC_k

 subroutine remap_delz_BC(pe_lag_BC, pe_eul_BC, delp_lag_BC, delz_lag_BC, delp_eul_BC, delz_eul_BC, npx, npy, npz, npz_coarse, bd, istag, jstag, iv, kord)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pe_lag_BC, delp_lag_BC, delz_lag_BC
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pe_eul_BC, delp_eul_BC, delz_eul_BC
   integer, intent(IN) :: npx, npy, npz, npz_coarse, istag, jstag, iv, kord

   integer :: i,j,k, istart, iend

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   if (is == 1) then
      call compute_specific_volume_BC_k(delp_lag_BC%west_t1, delz_lag_BC%west_t1, isd, 0, isd, 0, jsd, jed, npz_coarse)
      call remap_BC_k(pe_lag_BC%west_t1, pe_eul_BC%west_t1, delz_lag_BC%west_t1, delz_eul_BC%west_t1, isd, 0, isd, 0, jsd, jed+jstag, &
           npz, npz_coarse, iv, kord, log_pe=.false.)
      call compute_delz_BC_k(delp_eul_BC%west_t1, delz_eul_BC%west_t1, isd, 0, isd, 0, jsd, jed, npz)
   end if

   if (ie == npx-1) then
      call compute_specific_volume_BC_k(delp_lag_BC%east_t1, delz_lag_BC%east_t1, npx+istag, ied+istag, npx+istag, ied+istag, jsd, jed+jstag, npz_coarse)
      call remap_BC_k(pe_lag_BC%east_t1, pe_eul_BC%east_t1, delz_lag_BC%east_t1, delz_eul_BC%east_t1, npx+istag, ied+istag, npx+istag, ied+istag, jsd, jed+jstag, &
           npz, npz_coarse, iv, kord, log_pe=.false.)
      call compute_delz_BC_k(delp_eul_BC%east_t1, delz_eul_BC%east_t1, npx+istag, ied+istag, npx+istag, ied+istag, jsd, jed+jstag, npz)
   end if

   if (is == 1) then
      istart = is
   else
      istart = isd
   end if
   if (ie == npx-1) then
      iend = ie
   else
      iend = ied
   end if

   if (js == 1) then
      call compute_specific_volume_BC_k(delp_lag_BC%south_t1, delz_lag_BC%south_t1, isd, ied+istag, istart, iend+istag, jsd, 0, npz_coarse)
      call remap_BC_k(pe_lag_BC%south_t1, pe_eul_BC%south_t1, delz_lag_BC%south_t1, delz_eul_BC%south_t1, isd, ied+istag, istart, iend+istag, jsd, 0, npz, npz_coarse, &
           iv, kord, log_pe=.false.)
      call compute_delz_BC_k(delp_eul_BC%south_t1, delz_eul_BC%south_t1, isd, ied+istag, istart, iend+istag, jsd, 0, npz)
   end if

   if (je == npy-1) then
      call compute_specific_volume_BC_k(delp_lag_BC%north_t1, delz_lag_BC%north_t1, isd, ied+istag, istart, iend+istag, npy+jstag, jed+jstag, npz_coarse)
      call remap_BC_k(pe_lag_BC%north_t1, pe_eul_BC%north_t1, delz_lag_BC%north_t1, delz_eul_BC%north_t1, &
           isd, ied+istag, istart, iend+istag, npy+jstag, jed+jstag, npz, npz_coarse, iv, kord, log_pe=.false.)
      call compute_delz_BC_k(delp_eul_BC%north_t1, delz_eul_BC%north_t1, isd, ied+istag, istart, iend+istag, npy+jstag, jed+jstag, npz)
   end if
   
 end subroutine remap_delz_BC

 subroutine compute_specific_volume_BC_k(delpBC, delzBC, isd, ied, istart, iend, jstart, jend, npz)

   integer, intent(IN) :: isd, ied, istart, iend, jstart, jend, npz
   real, intent(IN)    :: delpBC(isd:ied,jstart:jend,npz)
   real, intent(INOUT) :: delzBC(isd:ied,jstart:jend,npz)

   character(len=120) :: errstring
   integer :: i,j,k
   
   do k=1,npz
   do j=jstart,jend
   do i=istart,iend
      delzBC(i,j,k) = -delzBC(i,j,k)/delpBC(i,j,k)
!!$!!! DEBUG CODE 
!!$      if (delzBC(i,j,k) <= 0. ) then
!!$         write(errstring,'(3I5, 2(2x, G))'), i, j, k, delzBC(i,j,k), delpBC(i,j,k)
!!$         call mpp_error(WARNING, ' Remap BC (sfc volume): invalid delz at '//errstring)               
!!$      endif
!!$!!! END DEBUG CODE
   end do
   end do
   end do

 end subroutine compute_specific_volume_BC_k

 subroutine compute_delz_BC_k(delpBC, delzBC, isd, ied, istart, iend, jstart, jend, npz)

   integer, intent(IN) :: isd, ied, istart, iend, jstart, jend, npz
   real, intent(IN)    :: delpBC(isd:ied,jstart:jend,npz)
   real, intent(INOUT) :: delzBC(isd:ied,jstart:jend,npz)
   
   character(len=120) :: errstring
   integer :: i,j,k
   
   do k=1,npz
   do j=jstart,jend
   do i=istart,iend
      delzBC(i,j,k) = -delzBC(i,j,k)*delpBC(i,j,k) 
!!$!!! DEBUG CODE
!!$      if (delzBC(i,j,k) >=0. ) then
!!$         write(errstring,'(3I5, 2(2x, G))'), i, j, k, delzBC(i,j,k), delpBC(i,j,k)
!!$         call mpp_error(WARNING, ' Remap BC (compute delz): invalid delz at '//errstring)               
!!$      endif
!!$!!! END DEBUG CODE
   end do
   end do
   end do

 end subroutine compute_delz_BC_k


 subroutine setup_pt_NH_BC(pt_BC, delp_BC, delz_BC, sphum_BC, q_BC, nq, &
#ifdef USE_COND
      q_con_BC, &
#ifdef MOIST_CAPPA
      cappa_BC, &
#endif
#endif
      npx, npy, npz, zvir, bd)

   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_nest_BC_type_3d), intent(IN), target    :: delp_BC, delz_BC, sphum_BC
   type(fv_nest_BC_type_3d), intent(INOUT), target :: pt_BC
   integer, intent(IN) :: nq
   type(fv_nest_BC_type_3d), intent(IN), target :: q_BC(nq)
#ifdef USE_COND
   type(fv_nest_BC_type_3d), intent(INOUT), target :: q_con_BC
#ifdef MOIST_CAPPA
   type(fv_nest_BC_type_3d), intent(INOUT), target :: cappa_BC
#endif
#endif
   integer, intent(IN) :: npx, npy, npz
   real, intent(IN) :: zvir

    real, parameter:: c_liq = 4185.5      ! heat capacity of water at 0C
    real, parameter:: c_ice = 1972.       ! heat capacity of ice at 0C: c=c_ice+7.3*(T-Tice) 
    real, parameter:: cv_vap = cp_vapor - rvgas  ! 1384.5

   real, dimension(:,:,:), pointer :: liq_watBC_west, ice_watBC_west, rainwatBC_west, snowwatBC_west, graupelBC_west
   real, dimension(:,:,:), pointer :: liq_watBC_east, ice_watBC_east, rainwatBC_east, snowwatBC_east, graupelBC_east
   real, dimension(:,:,:), pointer :: liq_watBC_north, ice_watBC_north, rainwatBC_north, snowwatBC_north, graupelBC_north
   real, dimension(:,:,:), pointer :: liq_watBC_south, ice_watBC_south, rainwatBC_south, snowwatBC_south, graupelBC_south

   real :: dp1, q_liq, q_sol, q_con = 0., cvm, pkz, rdg, cv_air

   integer :: i,j,k, istart, iend
   integer :: liq_wat, ice_wat, rainwat, snowwat, graupel
   real, parameter:: tice = 273.16 ! For GFS Partitioning
   real, parameter:: t_i0 = 15.

   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
   
   rdg = -rdgas / grav
   cv_air =  cp_air - rdgas

   liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
   ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
   rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
   snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
   graupel = get_tracer_index (MODEL_ATMOS, 'graupel')   

   if (is == 1) then
      if (.not. allocated(dum_West)) then
         allocate(dum_West(isd:0,jsd:jed,npz))
!$OMP parallel do default(none) shared(npz,isd,jsd,jed,dum_West)
         do k=1,npz
         do j=jsd,jed
         do i=isd,0
            dum_West(i,j,k) = 0.
         enddo
         enddo
         enddo
      endif
   endif
   if (js == 1) then
      if (.not. allocated(dum_South)) then
         allocate(dum_South(isd:ied,jsd:0,npz))
!$OMP parallel do default(none) shared(npz,isd,ied,jsd,dum_South)
         do k=1,npz
         do j=jsd,0
         do i=isd,ied
            dum_South(i,j,k) = 0.
         enddo
         enddo
         enddo
      endif
   endif
   if (ie == npx-1) then
      if (.not. allocated(dum_East)) then
         allocate(dum_East(npx:ied,jsd:jed,npz))
!$OMP parallel do default(none) shared(npx,npz,ied,jsd,jed,dum_East)
         do k=1,npz
         do j=jsd,jed
         do i=npx,ied
            dum_East(i,j,k) = 0.
         enddo
         enddo
         enddo
      endif
   endif
   if (je == npy-1) then
      if (.not. allocated(dum_North)) then
         allocate(dum_North(isd:ied,npy:jed,npz))
!$OMP parallel do default(none) shared(npy,npz,isd,ied,jed,dum_North)
         do k=1,npz
         do j=npy,jed
         do i=isd,ied
            dum_North(i,j,k) = 0.
         enddo
         enddo
         enddo
      endif
   endif

   if (liq_wat > 0) then
      liq_watBC_west  => q_BC(liq_wat)%west_t1
      liq_watBC_east  => q_BC(liq_wat)%east_t1
      liq_watBC_north => q_BC(liq_wat)%north_t1
      liq_watBC_south => q_BC(liq_wat)%south_t1
   else
      liq_watBC_west  => dum_west
      liq_watBC_east  => dum_east
      liq_watBC_north => dum_north
      liq_watBC_south => dum_south
   endif
   if (ice_wat > 0) then
      ice_watBC_west  => q_BC(ice_wat)%west_t1
      ice_watBC_east  => q_BC(ice_wat)%east_t1
      ice_watBC_north => q_BC(ice_wat)%north_t1
      ice_watBC_south => q_BC(ice_wat)%south_t1
   else
      ice_watBC_west  => dum_west
      ice_watBC_east  => dum_east
      ice_watBC_north => dum_north
      ice_watBC_south => dum_south
   endif
   if (rainwat > 0) then
      rainwatBC_west  => q_BC(rainwat)%west_t1
      rainwatBC_east  => q_BC(rainwat)%east_t1
      rainwatBC_north => q_BC(rainwat)%north_t1
      rainwatBC_south => q_BC(rainwat)%south_t1
   else
      rainwatBC_west  => dum_west
      rainwatBC_east  => dum_east
      rainwatBC_north => dum_north
      rainwatBC_south => dum_south
   endif
   if (snowwat > 0) then
      snowwatBC_west  => q_BC(snowwat)%west_t1
      snowwatBC_east  => q_BC(snowwat)%east_t1
      snowwatBC_north => q_BC(snowwat)%north_t1
      snowwatBC_south => q_BC(snowwat)%south_t1
   else
      snowwatBC_west  => dum_west
      snowwatBC_east  => dum_east
      snowwatBC_north => dum_north
      snowwatBC_south => dum_south
   endif
   if (graupel > 0) then
      graupelBC_west  => q_BC(graupel)%west_t1
      graupelBC_east  => q_BC(graupel)%east_t1
      graupelBC_north => q_BC(graupel)%north_t1
      graupelBC_south => q_BC(graupel)%south_t1
   else
      graupelBC_west  => dum_west
      graupelBC_east  => dum_east
      graupelBC_north => dum_north
      graupelBC_south => dum_south
   endif

   if (is == 1) then

      call setup_pt_NH_BC_k(pt_BC%west_t1, sphum_BC%west_t1, delp_BC%west_t1, delz_BC%west_t1, &
           liq_watBC_west, rainwatBC_west, ice_watBC_west, snowwatBC_west, graupelBC_west, &
#ifdef USE_COND
           q_con_BC%west_t1, &
#ifdef MOIST_CAPPA
           cappa_BC%west_t1, &
#endif
#endif
           zvir, isd, 0, isd, 0, jsd, jed, npz)
   end if


   if (js == 1) then
      if (is == 1) then
         istart = is
      else
         istart = isd
      end if
      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

      call setup_pt_NH_BC_k(pt_BC%south_t1, sphum_BC%south_t1, delp_BC%south_t1, delz_BC%south_t1, &
           liq_watBC_south, rainwatBC_south, ice_watBC_south, snowwatBC_south, graupelBC_south, &
#ifdef USE_COND
           q_con_BC%south_t1, &
#ifdef MOIST_CAPPA
           cappa_BC%south_t1, &
#endif
#endif
           zvir, isd, ied, istart, iend, jsd, 0, npz)
   end if


   if (ie == npx-1) then

      call setup_pt_NH_BC_k(pt_BC%east_t1, sphum_BC%east_t1, delp_BC%east_t1, delz_BC%east_t1, &
           liq_watBC_east, rainwatBC_east, ice_watBC_east, snowwatBC_east, graupelBC_east, &
#ifdef USE_COND
           q_con_BC%east_t1, &
#ifdef MOIST_CAPPA
           cappa_BC%east_t1, &
#endif
#endif
           zvir, npx, ied, npx, ied, jsd, jed, npz)
   end if

   if (je == npy-1) then
      if (is == 1) then
         istart = is
      else
         istart = isd
      end if
      if (ie == npx-1) then
         iend = ie
      else
         iend = ied
      end if

      call setup_pt_NH_BC_k(pt_BC%north_t1, sphum_BC%north_t1, delp_BC%north_t1, delz_BC%north_t1, &
           liq_watBC_north, rainwatBC_north, ice_watBC_north, snowwatBC_north, graupelBC_north, &
#ifdef USE_COND
           q_con_BC%north_t1, &
#ifdef MOIST_CAPPA
           cappa_BC%north_t1, &
#endif
#endif
           zvir, isd, ied, istart, iend, npy, jed, npz)
   end if

 end subroutine setup_pt_NH_BC


 subroutine setup_pt_NH_BC_k(ptBC,sphumBC,delpBC,delzBC, &
                             liq_watBC,rainwatBC,ice_watBC,snowwatBC,graupelBC, &
#ifdef USE_COND
                             q_conBC, &
#ifdef MOIST_CAPPA
                             cappaBC, &
#endif
#endif
                             zvir, isd, ied, istart, iend, jstart, jend, npz)

   integer, intent(IN) :: isd, ied, istart, iend, jstart, jend, npz
   real, intent(OUT), dimension(isd:ied,jstart:jend,npz) :: ptBC
   real, intent(IN),  dimension(isd:ied,jstart:jend,npz) :: sphumBC, delpBC, delzBC
   real, intent(IN),  dimension(isd:ied,jstart:jend,npz) :: liq_watBC,rainwatBC,ice_watBC,snowwatBC,graupelBC
#ifdef USE_COND
   real, intent(OUT), dimension(isd:ied,jstart:jend,npz) ::   q_conBC
#ifdef MOIST_CAPPA
   real, intent(OUT), dimension(isd:ied,jstart:jend,npz) ::   cappaBC
#endif
#endif
   real, intent(IN) :: zvir

   integer :: i,j,k
   real :: dp1, q_con, q_sol, q_liq, cvm, pkz, rdg, cv_air

   real, parameter:: c_liq = 4185.5      ! heat capacity of water at 0C
   real, parameter:: c_ice = 1972.       ! heat capacity of ice at 0C: c=c_ice+7.3*(T-Tice) 
   real, parameter:: cv_vap = cp_vapor - rvgas  ! 1384.5
   real, parameter:: tice = 273.16 ! For GFS Partitioning
   real, parameter:: t_i0 = 15.

   rdg = -rdgas / grav
   cv_air =  cp_air - rdgas

   do k=1,npz
   do j=jstart,jend
   do i=istart,iend
         dp1 = zvir*sphumBC(i,j,k)
#ifdef USE_COND
         q_liq = liq_watBC(i,j,k) + rainwatBC(i,j,k)
         q_sol = ice_watBC(i,j,k) + snowwatBC(i,j,k) + graupelBC(i,j,k)
         q_con = q_liq + q_sol
         q_conBC(i,j,k) = q_con
#ifdef MOIST_CAPPA
         cvm = (1.-(sphumBC(i,j,k)+q_con))*cv_air+sphumBC(i,j,k)*cv_vap+q_liq*c_liq+q_sol*c_ice
         cappaBC(i,j,k) = rdgas/(rdgas + cvm/(1.+dp1))
         pkz = exp( cappaBC(i,j,k)*log(rdg*delpBC(i,j,k)*ptBC(i,j,k) * &
              (1.+dp1)*(1.-q_con)/delzBC(i,j,k)))         
#else
         pkz = exp( kappa*log(rdg*delpBC(i,j,k)*ptBC(i,j,k) * &
              (1.+dp1)*(1.-q_con)/delzBC(i,j,k)))
#endif
         ptBC(i,j,k) = ptBC(i,j,k)*(1.+dp1)*(1.-q_con)/pkz
#else
         pkz = exp( kappa*log(rdg*delpBC(i,j,k)*ptBC(i,j,k) * &
              (1.+dp1)/delzBC(i,j,k)))
         ptBC(i,j,k) = ptBC(i,j,k)*(1.+dp1)/pkz
#endif
   end do
   end do
   end do

 end subroutine setup_pt_NH_BC_k

 subroutine set_NH_BCs_t0(neststruct)

   type(fv_nest_type), intent(INOUT) :: neststruct

#ifndef SW_DYNAMICS
   neststruct%delz_BC%east_t0  = neststruct%delz_BC%east_t1
   neststruct%delz_BC%west_t0  = neststruct%delz_BC%west_t1
   neststruct%delz_BC%north_t0 = neststruct%delz_BC%north_t1
   neststruct%delz_BC%south_t0 = neststruct%delz_BC%south_t1

   neststruct%w_BC%east_t0  = neststruct%w_BC%east_t1
   neststruct%w_BC%west_t0  = neststruct%w_BC%west_t1
   neststruct%w_BC%north_t0 = neststruct%w_BC%north_t1
   neststruct%w_BC%south_t0 = neststruct%w_BC%south_t1
#endif

 end subroutine set_NH_BCs_t0

 subroutine set_BCs_t0(ncnst, hydrostatic, neststruct)

   integer, intent(IN) :: ncnst
   logical, intent(IN) :: hydrostatic
   type(fv_nest_type), intent(INOUT) :: neststruct

   integer :: n

   neststruct%delp_BC%east_t0  = neststruct%delp_BC%east_t1
   neststruct%delp_BC%west_t0  = neststruct%delp_BC%west_t1
   neststruct%delp_BC%north_t0 = neststruct%delp_BC%north_t1
   neststruct%delp_BC%south_t0 = neststruct%delp_BC%south_t1
   do n=1,ncnst
      neststruct%q_BC(n)%east_t0  = neststruct%q_BC(n)%east_t1
      neststruct%q_BC(n)%west_t0  = neststruct%q_BC(n)%west_t1
      neststruct%q_BC(n)%north_t0 = neststruct%q_BC(n)%north_t1
      neststruct%q_BC(n)%south_t0 = neststruct%q_BC(n)%south_t1
   enddo
#ifndef SW_DYNAMICS
   neststruct%pt_BC%east_t0    = neststruct%pt_BC%east_t1
   neststruct%pt_BC%west_t0    = neststruct%pt_BC%west_t1
   neststruct%pt_BC%north_t0   = neststruct%pt_BC%north_t1
   neststruct%pt_BC%south_t0   = neststruct%pt_BC%south_t1
   neststruct%pt_BC%east_t0    = neststruct%pt_BC%east_t1
   neststruct%pt_BC%west_t0    = neststruct%pt_BC%west_t1
   neststruct%pt_BC%north_t0   = neststruct%pt_BC%north_t1
   neststruct%pt_BC%south_t0   = neststruct%pt_BC%south_t1

#ifdef USE_COND
   neststruct%q_con_BC%east_t0    = neststruct%q_con_BC%east_t1
   neststruct%q_con_BC%west_t0    = neststruct%q_con_BC%west_t1
   neststruct%q_con_BC%north_t0   = neststruct%q_con_BC%north_t1
   neststruct%q_con_BC%south_t0   = neststruct%q_con_BC%south_t1
#ifdef MOIST_CAPPA
   neststruct%cappa_BC%east_t0    = neststruct%cappa_BC%east_t1
   neststruct%cappa_BC%west_t0    = neststruct%cappa_BC%west_t1
   neststruct%cappa_BC%north_t0   = neststruct%cappa_BC%north_t1
   neststruct%cappa_BC%south_t0   = neststruct%cappa_BC%south_t1
#endif
#endif

   if (.not. hydrostatic) then
      call set_NH_BCs_t0(neststruct)
   endif
#endif
   neststruct%u_BC%east_t0  = neststruct%u_BC%east_t1
   neststruct%u_BC%west_t0  = neststruct%u_BC%west_t1
   neststruct%u_BC%north_t0 = neststruct%u_BC%north_t1
   neststruct%u_BC%south_t0 = neststruct%u_BC%south_t1
   neststruct%v_BC%east_t0  = neststruct%v_BC%east_t1
   neststruct%v_BC%west_t0  = neststruct%v_BC%west_t1
   neststruct%v_BC%north_t0 = neststruct%v_BC%north_t1
   neststruct%v_BC%south_t0 = neststruct%v_BC%south_t1


   neststruct%vc_BC%east_t0  = neststruct%vc_BC%east_t1
   neststruct%vc_BC%west_t0  = neststruct%vc_BC%west_t1
   neststruct%vc_BC%north_t0 = neststruct%vc_BC%north_t1
   neststruct%vc_BC%south_t0 = neststruct%vc_BC%south_t1
   neststruct%uc_BC%east_t0  = neststruct%uc_BC%east_t1
   neststruct%uc_BC%west_t0  = neststruct%uc_BC%west_t1
   neststruct%uc_BC%north_t0 = neststruct%uc_BC%north_t1
   neststruct%uc_BC%south_t0 = neststruct%uc_BC%south_t1

   neststruct%divg_BC%east_t0  = neststruct%divg_BC%east_t1
   neststruct%divg_BC%west_t0  = neststruct%divg_BC%west_t1
   neststruct%divg_BC%north_t0 = neststruct%divg_BC%north_t1
   neststruct%divg_BC%south_t0 = neststruct%divg_BC%south_t1

 end subroutine set_BCs_t0


!! nestupdate types
!! 1 - Interpolation update on all variables
!! 2 - Conserving update (over areas on cell-
!!     centered variables, over faces on winds) on all variables
!! 3 - Interpolation update on winds only
!! 4 - Interpolation update on all variables except delp (mass conserving)
!! 5 - Remap interpolating update, delp not updated
!! 6 - Remap conserving update, delp not updated
!! 7 - Remap conserving update, delp and q not updated
!! 8 - Remap conserving update, only winds updated

!! Note that nestupdate > 3 will not update delp.

!! "Remap update" remaps updated variables from the nested grid's
!!  vertical coordinate to that of the coarse grid. When delp is not
!!  updated (nestbctype >= 3) the vertical coordinates differ on
!!  the two grids, because the surface pressure will be different
!!  on the two grids.
!! Note: "conserving updates" do not guarantee global conservation
!!  unless flux nested grid BCs are specified, or if a quantity is
!!  not updated at all. This ability has not been implemented.

subroutine twoway_nesting(Atm, ngrids, grids_on_this_pe, zvir)

   type(fv_atmos_type), intent(INOUT) :: Atm(ngrids)
   integer, intent(IN) :: ngrids
   logical, intent(IN) :: grids_on_this_pe(ngrids)
   real, intent(IN) :: zvir

   integer :: n, p, sphum

   
   if (ngrids > 1) then

      do n=ngrids,2,-1 !loop backwards to allow information to propagate from finest to coarsest grids

         !two-way updating    
         if (Atm(n)%neststruct%twowaynest ) then
            if  (grids_on_this_pe(n) .or. grids_on_this_pe(Atm(n)%parent_grid%grid_number)) then
               sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
               call twoway_nest_update(Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, zvir, &
                    Atm(n)%ncnst, sphum, Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%omga, &
                    Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%uc, Atm(n)%vc, &
                    Atm(n)%pkz, Atm(n)%delz, Atm(n)%ps, Atm(n)%ptop, &
                    Atm(n)%gridstruct, Atm(n)%flagstruct, Atm(n)%neststruct, Atm(n)%parent_grid, Atm(N)%bd, .false.)
            endif
         endif

      end do

      !NOTE: these routines need to be used with any grid which has been updated to, not just the coarsest grid.
      do n=1,ngrids
         if (Atm(n)%neststruct%parent_of_twoway .and. grids_on_this_pe(n)) then
            call after_twoway_nest_update( Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%ng, Atm(n)%ncnst,   &
                 Atm(n)%u,  Atm(n)%v,  Atm(n)%w,  Atm(n)%delz, &
                 Atm(n)%pt,  Atm(n)%delp,  Atm(n)%q,   &
                 Atm(n)%ps,  Atm(n)%pe,  Atm(n)%pk,  Atm(n)%peln,  Atm(n)%pkz, &
                 Atm(n)%phis,  Atm(n)%ua,  Atm(n)%va,  &
                 Atm(n)%ptop, Atm(n)%gridstruct, Atm(n)%flagstruct, &
                 Atm(n)%domain, Atm(n)%bd)
         endif
      enddo

   endif ! ngrids > 1




  end subroutine twoway_nesting

!!!CLEANUP: this routine assumes that the PARENT GRID has pt = (regular) temperature,
!!!not potential temperature; which may cause problems when updating if this is not the case.
 subroutine twoway_nest_update(npx, npy, npz, zvir, ncnst, sphum,     &
                        u, v, w, omga, pt, delp, q,   &
                        uc, vc, pkz, delz, ps, ptop, &
                        gridstruct, flagstruct, neststruct, &
                        parent_grid, bd, conv_theta_in)

    real, intent(IN) :: zvir, ptop

    integer, intent(IN) :: npx, npy, npz
    integer, intent(IN) :: ncnst, sphum
    logical, intent(IN), OPTIONAL :: conv_theta_in

    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u ! D grid zonal wind (m/s)
    real, intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   bd%isd:        ,bd%jsd:        ,1: )  !  W (m/s)
    real, intent(inout) :: omga(bd%isd:bd%ied,bd%jsd:bd%jed,npz)      ! Vertical pressure velocity (pa/s)
    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) ! specific humidity and constituents
    real, intent(inout) :: uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) ! (uc,vc) C grid winds
    real, intent(inout) :: vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)

    real, intent(inout) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)             ! finite-volume mean pk
    real, intent(inout) :: delz(bd%isd:      ,bd%jsd:      ,1: )   ! delta-height (m); non-hydrostatic only
    real, intent(inout) :: ps  (bd%isd:bd%ied  ,bd%jsd:bd%jed)           ! Surface pressure (pascal)

    type(fv_grid_type), intent(INOUT) :: gridstruct
    type(fv_flags_type), intent(INOUT) :: flagstruct
    type(fv_nest_type), intent(INOUT) :: neststruct

    type(fv_atmos_type), intent(INOUT) :: parent_grid

    real, allocatable :: t_nest(:,:,:), ps0(:,:)
    integer :: i,j,k,n
    integer :: isd_p, ied_p, jsd_p, jed_p, isc_p, iec_p, jsc_p, jec_p
    integer :: isg, ieg, jsg,jeg, npx_p, npy_p
    integer :: istart, iend
    real :: qmass_b, qmass_a, fix = 1.
    logical :: used, conv_theta=.true.

    real :: qdp(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
    real, allocatable :: qdp_coarse(:,:,:)
    real(kind=f_p), allocatable :: q_diff(:,:,:)
    real :: L_sum_b(npz), L_sum_a(npz), blend_wt(npz)
    
    integer :: upoff
    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed
    integer :: isu, ieu, jsu, jeu

    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    isu = neststruct%isu
    ieu = neststruct%ieu
    jsu = neststruct%jsu
    jeu = neststruct%jeu

    upoff = neststruct%upoff

    !We update actual temperature, not theta.
    !If pt is actual temperature, set conv_theta to .false.
    if (present(conv_theta_in)) conv_theta = conv_theta_in

    if ((.not. neststruct%parent_proc) .and. (.not. neststruct%child_proc)) return

    call mpp_get_data_domain( parent_grid%domain, &
         isd_p,  ied_p,  jsd_p,  jed_p  )
    call mpp_get_compute_domain( parent_grid%domain, &
         isc_p,  iec_p,  jsc_p,  jec_p  )

    !!!FOR NOW:: Set blend_wt to 1.
    blend_wt(:) = 1.


   !delp/ps

   if (neststruct%nestupdate < 3) then

         call update_coarse_grid(parent_grid%delp, delp, neststruct%nest_domain,&
              neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
              isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, &
              neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, &
              npx, npy, npz, 0, 0, &
              neststruct%refinement, neststruct%nestupdate, upoff, 0, &
              neststruct%parent_proc, neststruct%child_proc, parent_grid, blend_wt)

      call mpp_sync!self

#ifdef SW_DYNAMICS
      if (neststruct%parent_proc) then
         do j=jsd_p,jed_p
            do i=isd_p,ied_p

               parent_grid%ps(i,j) = &
                    parent_grid%delp(i,j,1)/grav 

            end do
         end do
      endif
#endif

   end if

   !if (neststruct%nestupdate /= 3 .and. neststruct%nestbctype /= 3) then
   if (neststruct%nestupdate /= 3 .and. neststruct%nestupdate /= 7 .and. neststruct%nestupdate /= 8) then

      allocate(qdp_coarse(isd_p:ied_p,jsd_p:jed_p,npz))
      if (parent_grid%flagstruct%nwat > 0) then
         allocate(q_diff(isd_p:ied_p,jsd_p:jed_p,npz))
         q_diff = 0.
      endif

      do n=1,parent_grid%flagstruct%nwat

         qdp_coarse = 0.
         if (neststruct%child_proc) then
            do k=1,npz
            do j=jsd,jed
            do i=isd,ied
               qdp(i,j,k) = q(i,j,k,n)*delp(i,j,k)
            enddo
            enddo
            enddo
         else
            qdp = 0.
         endif

         if (neststruct%parent_proc) then
            !Add up ONLY region being replaced by nested grid
            do k=1,npz
            do j=jsu,jeu
            do i=isu,ieu
               qdp_coarse(i,j,k) = parent_grid%q(i,j,k,n)*parent_grid%delp(i,j,k)
            enddo
            enddo
            enddo
            call level_sum(qdp_coarse, parent_grid%gridstruct%area, parent_grid%domain, &
                 parent_grid%bd, npz, L_sum_b)
         else
            qdp_coarse = 0.
         endif
         if (neststruct%parent_proc) then
            if (n <= parent_grid%flagstruct%nwat) then
            do k=1,npz
            do j=jsu,jeu
            do i=isu,ieu
               q_diff(i,j,k) = q_diff(i,j,k) - qdp_coarse(i,j,k)
            enddo
            enddo
            enddo
            endif
         endif

            call update_coarse_grid(qdp_coarse, qdp, neststruct%nest_domain, &
                 neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
                 isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, &
                 neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, &
                 npx, npy, npz, 0, 0, &
                 neststruct%refinement, neststruct%nestupdate, upoff, 0, &
                 neststruct%parent_proc, neststruct%child_proc, parent_grid, blend_wt)

               call mpp_sync!self

         if (neststruct%parent_proc) then
            call level_sum(qdp_coarse, parent_grid%gridstruct%area, parent_grid%domain, &
                 parent_grid%bd, npz, L_sum_a)
            do k=1,npz
               if (L_sum_a(k) > 0.) then
                  fix = L_sum_b(k)/L_sum_a(k)
               do j=jsu,jeu
               do i=isu,ieu
                  !Normalization mass fixer
                  parent_grid%q(i,j,k,n) = qdp_coarse(i,j,k)*fix
            enddo
            enddo
               endif
            enddo
               if (n == 1) sphum_ll_fix = 1. - fix
         endif
         if (neststruct%parent_proc) then
            if (n <= parent_grid%flagstruct%nwat) then
            do k=1,npz
            do j=jsu,jeu
            do i=isu,ieu
               q_diff(i,j,k) = q_diff(i,j,k) + parent_grid%q(i,j,k,n)
            enddo
            enddo
            enddo
            endif
         endif

      end do

         if (neststruct%parent_proc) then
            if (parent_grid%flagstruct%nwat > 0) then
               do k=1,npz
            do j=jsu,jeu
            do i=isu,ieu
               parent_grid%delp(i,j,k) = parent_grid%delp(i,j,k) + q_diff(i,j,k)
            enddo
            enddo
            enddo
         endif

         do n=1,parent_grid%flagstruct%nwat
               do k=1,npz
         do j=jsu,jeu
         do i=isu,ieu
            parent_grid%q(i,j,k,n) = parent_grid%q(i,j,k,n)/parent_grid%delp(i,j,k)
         enddo
         enddo
         enddo               
         enddo
         endif

      deallocate(qdp_coarse)
      if  (allocated(q_diff)) deallocate(q_diff)

   endif

#ifndef SW_DYNAMICS
   if (neststruct%nestupdate /= 3 .and. neststruct%nestupdate /= 8) then

      if (conv_theta) then

         if (neststruct%child_proc) then
            !pt is potential temperature on the nested grid, but actual
            !temperature on the coarse grid. Compute actual temperature
            !on the nested grid, then gather.
            allocate(t_nest(isd:ied,jsd:jed,1:npz))
!$OMP parallel do default(none) shared(npz,js,je,is,ie,t_nest,pt,pkz,zvir,q,sphum)
            do k=1,npz
               do j=js,je
                  do i=is,ie
                     t_nest(i,j,k) = pt(i,j,k)*pkz(i,j,k)/(1.+zvir*q(i,j,k,sphum))
                  enddo
               enddo
            enddo
            deallocate(t_nest)
         endif

            call update_coarse_grid(parent_grid%pt, &
                 t_nest, neststruct%nest_domain, &
                 neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
                 isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, &
                 neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, &
                 npx, npy, npz, 0, 0, &
                 neststruct%refinement, neststruct%nestupdate, upoff, 0, &
                 neststruct%parent_proc, neststruct%child_proc, parent_grid, blend_wt)
      else

            call update_coarse_grid(parent_grid%pt, &
              pt, neststruct%nest_domain, &
              neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
              isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, &
              neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, &
              npx, npy, npz, 0, 0, &
              neststruct%refinement, neststruct%nestupdate, upoff, 0, &
              neststruct%parent_proc, neststruct%child_proc, parent_grid, blend_wt)

      endif !conv_theta

      call mpp_sync!self

      if (.not. flagstruct%hydrostatic) then

            call update_coarse_grid(parent_grid%w, w, neststruct%nest_domain, &
                 neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
                 isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, &
                 neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, &
                 npx, npy, npz, 0, 0, &
                 neststruct%refinement, neststruct%nestupdate, upoff, 0, &
                 neststruct%parent_proc, neststruct%child_proc, parent_grid, blend_wt)
            !Updating for delz not yet implemented; 
            ! may need to think very carefully how one would do this!!!
            ! consider updating specific volume instead?
!!$            call update_coarse_grid(parent_grid%delz, delz, neststruct%nest_domain, &
!!$                 neststruct%ind_update_h, &
!!$                 isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, npz, 0, 0, &
!!$                 neststruct%refinement, neststruct%nestupdate, upoff, 0, neststruct%parent_proc, neststruct%child_proc)

         call mpp_sync!self

      end if
      
   end if !Neststruct%nestupdate /= 3

#endif

      call update_coarse_grid(parent_grid%u, u, neststruct%nest_domain, &
           neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
           isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, &
           neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, &
           npx, npy, npz, 0, 1, &
           neststruct%refinement, neststruct%nestupdate, upoff, 0, &
           neststruct%parent_proc, neststruct%child_proc, parent_grid, blend_wt)

      call update_coarse_grid(parent_grid%v, v, neststruct%nest_domain, &
           neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
           isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, &
           neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, &
           npx, npy, npz, 1, 0, &
           neststruct%refinement, neststruct%nestupdate, upoff, 0, &
           neststruct%parent_proc, neststruct%child_proc, parent_grid, blend_wt)

   call mpp_sync!self


#ifndef SW_DYNAMICS
   if (neststruct%nestupdate >= 5 .and. npz > 4) then

      !Use PS0 from nested grid, NOT the full delp. Also we assume the same number of levels on both grids.
      !PS0 should be initially set to be ps so that this routine does NOTHING outside of the update region

      !Re-compute nested (AND COARSE) grid ps

      allocate(ps0(isd_p:ied_p,jsd_p:jed_p))
      if (neststruct%parent_proc) then

         parent_grid%ps = parent_grid%ptop
!This loop appears to cause problems with OMP
!$OMP parallel do default(none) shared(npz,jsd_p,jed_p,isd_p,ied_p,parent_grid)
         do j=jsd_p,jed_p
            do k=1,npz
               do i=isd_p,ied_p
                  parent_grid%ps(i,j) = parent_grid%ps(i,j) + &
                       parent_grid%delp(i,j,k)
               end do
            end do
         end do

         ps0 = parent_grid%ps
      endif

      if (neststruct%child_proc) then

         ps = ptop
!$OMP parallel do default(none) shared(npz,jsd,jed,isd,ied,ps,delp)
         do j=jsd,jed
            do k=1,npz
               do i=isd,ied
                  ps(i,j) = ps(i,j) + delp(i,j,k)
               end do
            end do
         end do
      endif

      call update_coarse_grid(ps0, ps, neststruct%nest_domain, &
              neststruct%ind_update_h, gridstruct%dx, gridstruct%dy, gridstruct%area, &
              isd_p, ied_p, jsd_p, jed_p, isd, ied, jsd, jed, &
              neststruct%isu, neststruct%ieu, neststruct%jsu, neststruct%jeu, &
              npx, npy, 0, 0, &
              neststruct%refinement, neststruct%nestupdate, upoff, 0, &
              neststruct%parent_proc, neststruct%child_proc, parent_grid, 1.)

      !!! The mpp version of update_coarse_grid does not return a consistent value of ps
      !!! across PEs, as it does not go into the haloes of a given coarse-grid PE. This
      !!! update_domains call takes care of the problem.

   if (neststruct%parent_proc) then
     call mpp_update_domains(parent_grid%ps, parent_grid%domain, complete=.false.)
     call mpp_update_domains(ps0, parent_grid%domain, complete=.true.)
   endif


      call mpp_sync!self

      if (parent_grid%tile == neststruct%parent_tile) then 

         if (neststruct%parent_proc) then

         !comment out if statement to always remap theta instead of t in the remap-update.
         !(In LtE typically we use remap_t = .true.: remapping t is better (except in
         !idealized simulations with a background uniform theta) since near the top
         !boundary theta is exponential, which is hard to accurately interpolate with a spline
         if (.not. parent_grid%flagstruct%remap_t) then
!$OMP parallel do default(none) shared(npz,jsc_p,jec_p,isc_p,iec_p,parent_grid,zvir,sphum)
            do k=1,npz
               do j=jsc_p,jec_p
                  do i=isc_p,iec_p
                     parent_grid%pt(i,j,k) = &
                          parent_grid%pt(i,j,k)/parent_grid%pkz(i,j,k)*&
                          (1.+zvir*parent_grid%q(i,j,k,sphum))
                  end do
               end do
            end do
         end if
         call update_remap_tqw(npz, parent_grid%ak, parent_grid%bk, &
              parent_grid%ps, parent_grid%delp, &
              parent_grid%pt, parent_grid%q, parent_grid%w, &
              parent_grid%flagstruct%hydrostatic, &
              npz, ps0, zvir, parent_grid%ptop, ncnst, &
              parent_grid%flagstruct%kord_tm, parent_grid%flagstruct%kord_tr, &
              parent_grid%flagstruct%kord_wz, &
              isc_p, iec_p, jsc_p, jec_p, isd_p, ied_p, jsd_p, jed_p, .false. ) !neststruct%nestupdate < 7)
         if (.not. parent_grid%flagstruct%remap_t) then
!$OMP parallel do default(none) shared(npz,jsc_p,jec_p,isc_p,iec_p,parent_grid,zvir,sphum)
            do k=1,npz
               do j=jsc_p,jec_p
                  do i=isc_p,iec_p
                     parent_grid%pt(i,j,k) = &
                          parent_grid%pt(i,j,k)*parent_grid%pkz(i,j,k) / &
                          (1.+zvir*parent_grid%q(i,j,k,sphum))
                  end do
               end do
            end do
         end if

         call update_remap_uv(npz, parent_grid%ak, parent_grid%bk, &
              parent_grid%ps, &
              parent_grid%u, &
              parent_grid%v, npz, ps0, parent_grid%flagstruct%kord_mt, &
              isc_p, iec_p, jsc_p, jec_p, isd_p, ied_p, jsd_p, jed_p, parent_grid%ptop)

         endif !neststruct%parent_proc

      end if

      if (allocated(ps0)) deallocate(ps0)

   end if

#endif

 end subroutine twoway_nest_update

 subroutine level_sum(q, area, domain, bd, npz, L_sum)

    integer, intent(IN) :: npz
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(in) :: area(   bd%isd:bd%ied  ,bd%jsd:bd%jed)
    real, intent(in) ::    q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)
    real, intent(OUT) :: L_sum( npz ) 
    type(domain2d), intent(IN) :: domain
   
    integer :: i, j, k, n
    real :: qA!(bd%is:bd%ie, bd%js:bd%je)

    do k=1,npz
       qA = 0.
       do j=bd%js,bd%je
       do i=bd%is,bd%ie
          !qA(i,j) = q(i,j,k)*area(i,j)
          qA = qA + q(i,j,k)*area(i,j)
       enddo
       enddo
       call mp_reduce_sum(qA)
       L_sum(k) = qA
!       L_sum(k) = mpp_global_sum(domain, qA, flags=BITWISE_EXACT_SUM)
!       L_sum(k) = mpp_global_sum(domain, qA, flags=BITWISE_EFP_SUM) ! doesn't work??
    enddo

 end subroutine level_sum


 subroutine after_twoway_nest_update(npx, npy, npz, ng, ncnst,   &
                        u, v, w, delz, pt, delp, q,              &
                        ps, pe, pk, peln, pkz, phis, ua, va,     &
                        ptop, gridstruct, flagstruct,            &
                        domain, bd)

   type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(IN) :: ptop

    integer, intent(IN) :: ng, npx, npy, npz
    integer, intent(IN) :: ncnst

    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u ! D grid zonal wind (m/s)
    real, intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v ! D grid meridional wind (m/s)
    real, intent(inout) :: w(   bd%isd:        ,bd%jsd:        ,1: )  !  W (m/s)
    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! temperature (K)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! pressure thickness (pascal)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) ! specific humidity and constituents
    real, intent(inout) :: delz(bd%isd:        ,bd%jsd:        ,1: )   ! delta-height (m); non-hydrostatic only

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout) :: ps  (bd%isd:bd%ied  ,bd%jsd:bd%jed)           ! Surface pressure (pascal)
    real, intent(inout) :: pe  (bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1)  ! edge pressure (pascal)
    real, intent(inout) :: pk  (bd%is:bd%ie,bd%js:bd%je, npz+1)          ! pe**cappa
    real, intent(inout) :: peln(bd%is:bd%ie,npz+1,bd%js:bd%je)           ! ln(pe)
    real, intent(inout) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)             ! finite-volume mean pk
    
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real, intent(inout) :: phis(bd%isd:bd%ied,bd%jsd:bd%jed)       ! Surface geopotential (g*Z_surf)

    real, intent(inout), dimension(bd%isd:bd%ied ,bd%jsd:bd%jed ,npz):: ua, va
    type(fv_grid_type), intent(IN) :: gridstruct
    type(fv_flags_type), intent(IN) :: flagstruct
    type(domain2d), intent(INOUT) :: domain

    logical :: bad_range

    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed
    
    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    call cubed_to_latlon(u, v, ua, va, &
         gridstruct, npx, npy, npz, &
         1, gridstruct%grid_type, domain, &
         gridstruct%nested, flagstruct%c2l_ord, bd)

#ifndef SW_DYNAMICS

   !To get coarse grid pkz, etc right after a two-way update so
   !that it is consistent across a restart:
   !(should only be called after doing such an update)

    !! CLEANUP: move to twoway_nest_update??
   call p_var(npz, is, ie, js, je, ptop, ptop_min,  &
        delp, delz, &
        pt, ps, &
        pe, peln,   &
        pk,   pkz, kappa, &
        q, ng, flagstruct%ncnst,  gridstruct%area_64, 0.,  &
        .false.,  .false., & !mountain argument not used
        flagstruct%moist_phys,  flagstruct%hydrostatic, &
        flagstruct%nwat, domain, .false.)

#endif

      if (flagstruct%range_warn) then
         call range_check('TA update', pt, is, ie, js, je, ng, npz, gridstruct%agrid, 130., 350., bad_range)
         call range_check('UA update', ua, is, ie, js, je, ng, npz, gridstruct%agrid, -220., 250., bad_range)
         call range_check('VA update', va, is, ie, js, je, ng, npz, gridstruct%agrid, -220., 220., bad_range)
         if (.not. flagstruct%hydrostatic) then
            call range_check('W update', w, is, ie, js, je, ng, npz, gridstruct%agrid, -50., 100., bad_range)
         endif
      endif



 end subroutine after_twoway_nest_update

 !Routines for remapping (interpolated) nested-grid data to the coarse-grid's vertical coordinate.

 !This does not yet do anything for the tracers
 subroutine update_remap_tqw( npz, ak,  bk,  ps, delp,  t,  q, w, hydrostatic, &
                      kmd, ps0, zvir, ptop, nq, kord_tm, kord_tr, kord_wz, &
                      is, ie, js, je, isd, ied, jsd, jed, do_q)
  integer, intent(in):: npz, kmd, nq, kord_tm, kord_tr, kord_wz
  real,    intent(in):: zvir, ptop
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(in), dimension(isd:ied,jsd:jed):: ps0
  real,    intent(in), dimension(isd:ied,jsd:jed):: ps
  real, intent(in), dimension(isd:ied,jsd:jed,npz):: delp
  real,    intent(inout), dimension(isd:ied,jsd:jed,npz):: t, w
  real,    intent(inout), dimension(isd:ied,jsd:jed,npz,nq):: q
  integer,  intent(in) ::  is, ie, js, je, isd, ied, jsd, jed
  logical,   intent(in) :: hydrostatic, do_q
! local:
  real, dimension(is:ie,kmd):: tp, qp
  real, dimension(is:ie,kmd+1):: pe0, pn0
  real, dimension(is:ie,npz):: qn1
  real, dimension(is:ie,npz+1):: pe1, pn1
  integer i,j,k,iq

!$OMP parallel do default(none) shared(js,je,kmd,is,ie,ak,bk,ps0,q,npz,ptop,do_q,&
!$OMP          t,w,ps,nq,hydrostatic,kord_tm,kord_tr,kord_wz) &
!$OMP          private(pe0,pn0,pe1,pn1,qp,tp,qn1)
  do 5000 j=js,je

     do k=1,kmd+1
        do i=is,ie
           pe0(i,k) = ak(k) + bk(k)*ps0(i,j)
           pn0(i,k) = log(pe0(i,k))
       enddo
     enddo 
     do k=1,kmd+1
        do i=is,ie
           pe1(i,k) = ak(k) + bk(k)*ps(i,j)
           pn1(i,k) = log(pe1(i,k))
       enddo
     enddo 
     if (do_q) then
        do iq=1,nq
        do k=1,kmd
        do i=is,ie
           qp(i,k) = q(i,j,k,iq)
        enddo
        enddo
        call mappm(kmd, pe0, qp, npz, pe1,  qn1, is,ie, 0, kord_tr, ptop)
        do k=1,npz
           do i=is,ie
              q(i,j,k,iq) = qn1(i,k)
           enddo
        enddo
        enddo
     endif

     do k=1,kmd
        do i=is,ie
           tp(i,k) = t(i,j,k)
        enddo
     enddo
     !Remap T using logp
     call mappm(kmd, pn0, tp, npz, pn1, qn1, is,ie, 1, abs(kord_tm), ptop)
     
     do k=1,npz
        do i=is,ie
           t(i,j,k) = qn1(i,k)
        enddo
     enddo

     if (.not. hydrostatic) then
        do k=1,kmd
           do i=is,ie
              tp(i,k) = w(i,j,k)
           enddo
        enddo
        !Remap w using p
        !Using iv == -1 instead of -2
        call mappm(kmd, pe0, tp, npz, pe1, qn1, is,ie, -1, kord_wz, ptop)

        do k=1,npz
           do i=is,ie
              w(i,j,k) = qn1(i,k)
           enddo
        enddo
     endif

5000 continue

 end subroutine update_remap_tqw

 !remap_uv as-is remaps only a-grid velocities. A new routine has been written to handle staggered grids.
 subroutine update_remap_uv(npz, ak, bk, ps, u, v, kmd, ps0, kord_mt, &
                      is, ie, js, je, isd, ied, jsd, jed, ptop)
  integer, intent(in):: npz
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(in):: ps(isd:ied,jsd:jed)
  real,    intent(inout), dimension(isd:ied,jsd:jed+1,npz):: u
  real,    intent(inout), dimension(isd:ied+1,jsd:jed,npz):: v
!
  integer, intent(in):: kmd, kord_mt
  real,    intent(IN) :: ptop
  real,    intent(in):: ps0(isd:ied,jsd:jed)
  integer,  intent(in) ::  is, ie, js, je, isd, ied, jsd, jed
!
! local:
  real, dimension(is:ie+1,kmd+1):: pe0
  real, dimension(is:ie+1,npz+1):: pe1
  real, dimension(is:ie+1,kmd):: qt
  real, dimension(is:ie+1,npz):: qn1
  integer i,j,k

!------
! map u
!------
!$OMP parallel do default(none) shared(js,je,kmd,is,ie,ak,bk,ps,ps0,npz,u,ptop,kord_mt) &
!$OMP          private(pe0,pe1,qt,qn1)
  do j=js,je+1
!------
! Data
!------
     do k=1,kmd+1
       do i=is,ie
          pe0(i,k) = ak(k) + bk(k)*0.5*(ps0(i,j)+ps0(i,j-1))
       enddo
     enddo
!------
! Model
!------
     do k=1,kmd+1
        do i=is,ie
          pe1(i,k) = ak(k) + bk(k)*0.5*(ps(i,j)+ps(i,j-1))
       enddo
     enddo
!------
!Do map
!------
     qt = 0.
      do k=1,kmd
         do i=is,ie
            qt(i,k) = u(i,j,k)
         enddo
      enddo
      qn1 = 0.
      call mappm(kmd, pe0(is:ie,:), qt(is:ie,:), npz, pe1(is:ie,:), qn1(is:ie,:), is,ie, -1, kord_mt, ptop)
      do k=1,npz
         do i=is,ie
            u(i,j,k) = qn1(i,k)
         enddo
      enddo

   end do

!------
! map v
!------
!$OMP parallel do default(none) shared(js,je,kmd,is,ie,ak,bk,ps,ps0,npz,v,ptop) &
!$OMP          private(pe0,pe1,qt,qn1)
   do j=js,je
!------
! Data
!------
     do k=1,kmd+1
        do i=is,ie+1
          pe0(i,k) = ak(k) + bk(k)*0.5*(ps0(i,j)+ps0(i-1,j))
       enddo
     enddo
!------
! Model
!------
     do k=1,kmd+1
        do i=is,ie+1
          pe1(i,k) = ak(k) + bk(k)*0.5*(ps(i,j)+ps(i-1,j))
       enddo
     enddo
!------
!Do map
!------
     qt = 0.
      do k=1,kmd
         do i=is,ie+1
            qt(i,k) = v(i,j,k)
         enddo
      enddo
      qn1 = 0.
      call mappm(kmd, pe0(is:ie+1,:), qt(is:ie+1,:), npz, pe1(is:ie+1,:), qn1(is:ie+1,:), is,ie+1, -1, 8, ptop)
      do k=1,npz
         do i=is,ie+1
            v(i,j,k) = qn1(i,k)
         enddo
      enddo
   end do

 end subroutine update_remap_uv


end module fv_nesting_mod
