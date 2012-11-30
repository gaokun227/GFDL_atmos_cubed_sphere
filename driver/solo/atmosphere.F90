!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module atmosphere_mod

!-----------------------------------------------------------------------
!
!    interface for FV dynamical core with Held-Suarez forcing
!
!-----------------------------------------------------------------------


use constants_mod, only: grav, kappa, cp_air, pi, rdgas, rvgas, SECONDS_PER_DAY, radius
use fms_mod,       only: file_exist, open_namelist_file,   &
                         error_mesg, FATAL,                &
                         check_nml_error, stdlog, stdout,  &
                         write_version_number,             &
                         close_file, set_domain, nullify_domain, mpp_pe, mpp_root_pe, &
                         mpp_error, NOTE
use time_manager_mod, only: time_type, get_time, set_time, operator(+)
use mpp_domains_mod,  only: domain2d
use mpp_mod,          only: input_nml_file
!------------------
! FV specific codes:
!------------------
use fv_arrays_mod, only: fv_atmos_type, Atm
use fv_current_grid_mod, only: domain
use fv_control_mod,only: fv_init, fv_end, adiabatic, p_ref
use fv_phys_mod,   only: fv_phys, fv_nudge
use fv_diagnostics_mod, only: fv_diag_init, fv_diag, fv_time
use fv_timing_mod,   only: timing_on, timing_off
use fv_restart_mod, only: fv_restart
use fv_dynamics_mod, only: fv_dynamics
use fv_nesting_mod, only: twoway_nest_update, before_twoway_nest_update, after_twoway_nest_update
use fv_grid_tools_mod, only: grid_type
use lin_cld_microphys_mod, only: lin_cld_microphys_init, lin_cld_microphys_end
use fv_nwp_nudge_mod,   only: fv_nwp_nudge_init, fv_nwp_nudge_end
use fv_mp_mod, only: switch_current_Atm, domain
use fv_mp_mod, only: grids_on_this_pe, concurrent, gid
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index
use fv_grid_utils_mod,  only: cubed_to_latlon
!-----------------------------------------------------------------------

implicit none
private

public   atmosphere_init, atmosphere,  atmosphere_end, atmosphere_domain

!-----------------------------------------------------------------------
!---- private data ----

type        (time_type) :: Time_step_atmos
real                    :: dt_atmos
integer :: sec
integer days, seconds

logical :: cold_start      = .false.       ! read in initial condition
integer :: ntiles=1
!-----------------------------------------------------------------------

!---- version number -----
character(len=128) :: version = '$Id: atmosphere.F90,v 17.0.2.1.2.4.2.8.4.1.2.1 2012/09/28 17:38:04 Rusty.Benson Exp $'
character(len=128) :: tagname = '$Name: siena_201211 $'

contains

!#######################################################################

  subroutine atmosphere_init ( Time_init, Time, Time_step )

    type (time_type), intent(in) :: Time_step
    type (time_type), intent(in) :: Time_init
    type (time_type), intent(in) :: Time

    ! local:
    integer :: axes(4)
    integer :: ss, ds
    integer i,j, isc, iec, jsc, jec
    real pp(2)
    logical:: tracers_mca(7)
    real:: zvir
    integer :: f_unit, ios, ierr, n

    namelist /nest_nml/ ntiles !This is an OPTIONAL namelist, that needs to be read before everything else

  !----- write version and namelist to log file -----

    call write_version_number ( version, tagname )

  !---- compute physics/atmos time step in seconds ----

    Time_step_atmos = Time_step
    call get_time (Time_step_atmos, sec)
    dt_atmos = real(sec)

  !----- initialize FV dynamical core -----
!   cold_start = (.not.file_exist('INPUT/fv_rst.res.nc').and. .not.file_exist('INPUT/'//fms_tracers_file))
    cold_start = (.not.file_exist('INPUT/fv_core.res.nc') .and. .not.file_exist('INPUT/fv_core.res.tile1.nc'))

#ifdef INTERNAL_FILE_NML
      read (input_nml_file,nest_nml,iostat=ios)
      ierr = check_nml_error(ios,'nest_nml')
#else
      f_unit=open_namelist_file()
      rewind (f_unit)
      read (f_unit,nest_nml,iostat=ios)
      ierr = check_nml_error(ios,'nest_nml')
      call close_file(f_unit)
#endif

    allocate(Atm(ntiles))

    call fv_init(Atm, dt_atmos)  ! allocates Atm components

    do n=1,ntiles
       Atm(n)%moist_phys = .false.
    end do

    isc = Atm(1)%isc
    iec = Atm(1)%iec
    jsc = Atm(1)%jsc
    jec = Atm(1)%jec

         call timing_on('fv_restart')
    call fv_restart(domain, Atm, dt_atmos, seconds, days, cold_start, grid_type)
         call timing_off('fv_restart')

!   Atm(1)%phis(:,:) = Atm(1)%phis(:,:)*radius / 6371.0e3

     fv_time = time

     do n=1,ntiles
       if (.not. grids_on_this_pe(n)) then
          cycle
       endif
       call switch_current_Atm(Atm(n))
       call fv_diag_init(Atm(n:n), axes, Time, Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, p_ref)

    if ( adiabatic .or. Atm(n)%do_Held_Suarez ) then
         zvir = 0.         ! no virtual effect
    else
         zvir = rvgas/rdgas - 1.
    endif

    !!! The following routines only work with single grids as of yet
    if ( Atm(n)%nwat==6 )    & 
    call lin_cld_microphys_init(iec-isc+1, jec-jsc+1, Atm(n)%npz, axes, Time)

    if ( Atm(n)%nudge )    &
         call fv_nwp_nudge_init( Time, axes, Atm(n)%npz, zvir, Atm(n)%ak, Atm(n)%bk, Atm(n)%ts, Atm(n)%phis)

!   if( nlev > 1 ) call hs_forcing_init ( axes, Time )

    enddo

!-----------------------------------------------------------------------

  end subroutine atmosphere_init


!#######################################################################

  subroutine atmosphere (Time)
    type(time_type), intent(in) :: Time

    real:: zvir
    real:: time_total
    real:: tau_winds, tau_press, tau_temp
    logical:: strat = .false.
    integer :: n, sphum, p, nc

#ifdef NUDGE_IC
    tau_winds =  3600.
    tau_press = -1.
    tau_temp  = -1.
#else
    tau_winds = -1.
    tau_press = -1.
    tau_temp  = -1.
#endif

    call switch_current_Atm(Atm(1)) 

    fv_time = Time + Time_step_atmos
    call get_time (fv_time, seconds,  days)

    time_total = days*SECONDS_PER_DAY + seconds

    do n=1,ntiles

    call switch_current_Atm(Atm(n)) 
       
    if (.not. grids_on_this_pe(n)) then
       cycle
    endif

    call set_domain(Atm(n)%domain)  ! needed for diagnostic output done in fv_dynamics
    if ( tau_winds>0. .or. tau_press>0. .or. tau_temp>0. )     &
    call  fv_nudge(Atm(n)%npz, Atm(n)%isc, Atm(n)%iec, Atm(n)%jsc, Atm(n)%jec, Atm(n)%ng, &
                   Atm(n)%u, Atm(n)%v, Atm(n)%delp, Atm(n)%pt, dt_atmos,    &
                   tau_winds, tau_press, tau_temp)

  !---- call fv dynamics -----
    if ( adiabatic .or. Atm(n)%do_Held_Suarez ) then
         zvir = 0.         ! no virtual effect
    else
         zvir = rvgas/rdgas - 1.
    endif

    call set_domain(Atm(n)%domain)  ! needed for diagnostic output done in fv_dynamics
    call timing_on('fv_dynamics')
    call fv_dynamics(Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%ncnst, Atm(n)%ng,   & 
                     dt_atmos, Atm(n)%consv_te, Atm(n)%fill, Atm(n)%reproduce_sum, kappa,   &
                     cp_air, zvir, Atm(n)%ks, Atm(n)%ncnst, Atm(n)%n_split, Atm(n)%q_split, &
                     Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%delz,       &
                     Atm(n)%hydrostatic, Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%ps,       &
                     Atm(n)%pe, Atm(n)%pk, Atm(n)%peln, Atm(n)%pkz,                         &
                     Atm(n)%phis, Atm(n)%omga, Atm(n)%ua, Atm(n)%va, Atm(n)%uc, Atm(n)%vc,  &
                     Atm(n)%ak, Atm(n)%bk, Atm(n)%mfx, Atm(n)%mfy, Atm(n)%cx, Atm(n)%cy,    &
                     Atm(n)%ze0, Atm(n)%hybrid_z, time_total)
    call timing_off('fv_dynamics')

 end do
    do n=1,ntiles

       if (.not. grids_on_this_pe(n)) then
          cycle
       endif
       call switch_current_Atm(Atm(n))

    if(Atm(n)%npz /=1 .and. .not. adiabatic)then

                                                         call timing_on('FV_PHYS')
       call fv_phys(Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%isc, Atm(n)%iec,  &
                    Atm(n)%jsc, Atm(n)%jec, Atm(n)%ng, Atm(n)%ncnst,             &
                    Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%pt, Atm(n)%q, Atm(n)%pe,  &
                    Atm(n)%delp, Atm(n)%peln, Atm(n)%pkz, dt_atmos,              &
                    Atm(n)%ua, Atm(n)%va, Atm(n)%phis, Atm(n)%agrid,             &
                    Atm(n)%ak, Atm(n)%bk, Atm(n)%ks, Atm(n)%ps, Atm(n)%pk,       &
                    Atm(n)%u_srf, Atm(n)%v_srf,  Atm(n)%ts, Atm(n)%delz,         &
                    Atm(n)%hydrostatic, Atm(n)%oro, strat, .false., p_ref,     &
                    Atm(n)%fv_sg_adj, (mpp_pe()==mpp_root_pe()), Atm(n)%do_Held_Suarez,  &
                    fv_time, time_total)
                                                        call timing_off('FV_PHYS')
    endif

    call nullify_domain()
    end do

    do n=2,ntiles
       if (Atm(n)%parent_grid%grid_number == 1 .and. Atm(n)%twowaynest) then
          !Call this routine only if the coarsest grid has a two-way nested child grid
          call switch_current_Atm(Atm(1)) 
          if (grids_on_this_pe(1)) &
               call before_twoway_nest_update(Atm(1)%npx, Atm(1)%npy, Atm(1)%npz, Atm(1)%ng, Atm(1)%consv_te,               &
               kappa, cp_air, zvir, Atm(1)%ncnst,   &
               Atm(1)%u, Atm(1)%v, Atm(1)%w, Atm(1)%delz, Atm(1)%hydrostatic, Atm(1)%pt, Atm(1)%delp, Atm(1)%q,   &
               Atm(1)%ps, Atm(1)%pe, Atm(1)%pk, Atm(1)%peln, Atm(1)%pkz, Atm(1)%phis, Atm(1)%ua, Atm(1)%va, &
               Atm(1)%dry_mass, Atm(1)%grid_number, Atm(1)%mountain, Atm(1)%make_nh)
          exit
       endif
    enddo

    do n=ntiles,1,-1 !loop backwards to allow information to propagate from finest to coarsest grids

       !two-way updating    
       if (Atm(n)%nested .and. Atm(n)%twowaynest ) then
          if  ((.not. concurrent) .or. grids_on_this_pe(n) .or. ANY(gid == Atm(n)%parent_grid%pelist) ) then
             call switch_current_Atm(Atm(n))
             sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
             call twoway_nest_update(Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, cp_air, zvir, &
                  Atm(n)%ncnst, sphum, Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%omga, &
                  Atm(n)%hydrostatic, Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%uc, Atm(n)%vc, &
                  kappa, Atm(n)%pkz, Atm(n)%delz, Atm(n)%ps, .false.)
          endif
       endif

    end do

 !NOTE: these routines need to be used with any grid which has been updated to, not just the coarsest grid.
    do n=1,ntiles
       do p=1,size(Atm(n)%child_grids)
          if (.not. Atm(n)%child_grids(p)) cycle
          !Call this routine only if the coarsest grid has a two-way nested child grid
          if (Atm(p)%twowaynest) then
             call switch_current_Atm(Atm(n)) 
             if (grids_on_this_pe(n)) &
                  call after_twoway_nest_update( Atm(n)%npz, Atm(n)%ncnst,  Atm(n)%ng, dt_atmos,  &
                  Atm(n)%consv_te,  Atm(n)%fill, &
                  Atm(n)%reproduce_sum, kappa, cp_air, zvir, Atm(n)%ks,  Atm(n)%ncnst,   &
                  Atm(n)%u,  Atm(n)%v,  Atm(n)%w,  Atm(n)%delz,  Atm(n)%hydrostatic, &
                  Atm(n)%pt,  Atm(n)%delp,  Atm(n)%q,   &
                  Atm(n)%ps,  Atm(n)%pe,  Atm(n)%pk,  Atm(n)%peln,  Atm(n)%pkz,  &
                  Atm(n)%phis,  Atm(n)%omga,  Atm(n)%ua,  Atm(n)%va,  Atm(n)%uc,  Atm(n)%vc,          &
                  Atm(n)%ak,  Atm(n)%bk,  Atm(n)%ze0,  Atm(n)%hybrid_z,  Atm(n)%dry_mass, Atm(n)%adjust_dry_mass,&
                  Atm(n)%grid_number,  Atm(n)%mountain,  Atm(n)%make_nh)
             exit
          endif
       enddo
    enddo

  !---- diagnostics for FV dynamics -----

    do n=1,ntiles

       if (.not. grids_on_this_pe(n)) then
          cycle
       endif

       call switch_current_Atm(Atm(n)) 
       call nullify_domain()
       call timing_on('FV_DIAG')
       call cubed_to_latlon(Atm(n)%u, Atm(n)%v, Atm(n)%ua, Atm(n)%va, &
            Atm(n)%dx, Atm(n)%dy, Atm(n)%rdxa, Atm(n)%rdya, Atm(n)%npz, 1)
       call fv_diag(Atm(n:n), zvir, fv_time, Atm(n)%print_freq)
       
       call timing_off('FV_DIAG')
    end do

 end subroutine atmosphere


 subroutine atmosphere_end

    call get_time (fv_time, seconds,  days)

    if ( Atm(1)%nwat==6 ) call lin_cld_microphys_end

    call fv_end(Atm)
    deallocate(Atm)

  end subroutine atmosphere_end

 subroutine atmosphere_domain ( fv_domain )
 type(domain2d), intent(out) :: fv_domain

!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos
        
   fv_domain = domain
        
 end subroutine atmosphere_domain

end module atmosphere_mod
