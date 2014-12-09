module fv_phys_mod

use constants_mod,         only: grav, rdgas, rvgas, pi, cp_air, cp_vapor, hlv, radius, kappa
use time_manager_mod,      only: time_type, get_time
#ifdef MARS_GCM
use  Mars_phys_mod,     only:   Mars_phys
#else
use lin_cld_microphys_mod, only: lin_cld_microphys_driver, sg_conv, qsmith
use hswf_mod,              only: Held_Suarez_Tend
#endif
use fv_sg_mod,             only: fv_dry_conv
use fv_update_phys_mod,    only: fv_update_phys
use fv_timing_mod,         only: timing_on, timing_off
use monin_obukhov_mod,     only: mon_obkv
use tracer_manager_mod,    only: get_tracer_index
use field_manager_mod,     only: MODEL_ATMOS
use fms_mod,               only: error_mesg, FATAL, file_exist, open_namelist_file,  &
                                 check_nml_error, mpp_pe, mpp_root_pe, close_file, &
                                 write_version_number, stdlog, mpp_error
use fv_mp_mod,             only: is_master
use fv_diagnostics_mod,    only: prt_maxmin
use diag_manager_mod,  only: register_diag_field, register_static_field, send_data

use fv_arrays_mod,          only: fv_grid_type, fv_flags_type, fv_nest_type, fv_grid_bounds_type
use mpp_domains_mod,       only: domain2d

implicit none

  logical:: nudge_initialized = .false.
  logical:: sim_phys_initialized = .false.
  logical:: g0_sum_initialized = .false.
  real, allocatable, dimension(:,:) :: l_area
  logical:: master
  real, allocatable:: u0(:,:,:), v0(:,:,:), t0(:,:,:), dp(:,:,:), u_star(:,:)

public :: fv_phys, fv_nudge

!---- version number -----
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

  integer:: sphum, liq_wat, rainwat, snowwat, graupel, ice_wat, cld_amt
! For nudging the IC to a steady state:
  real :: tau_winds = 25.
  real :: tau_temp  = -1.
  real :: tau_press = -1.
!-----
  real, parameter:: c_liq = 4218.
  real, parameter:: c_ice = 2106.
  real, parameter:: cv_vap = 1389.0     ! Satoh value
  real, parameter:: dc_vap = cp_vapor - c_liq     ! = -2368.
  real, parameter:: hlv0 = 2.501e6   ! Emanual Appendix-2
  real, parameter:: Lv0 =  hlv0 - dc_vap*273.15   ! = 3.147782e6
  real::  zvir
  real:: tau_difz = 600.
  real:: area_ref = 6.25E8
  real:: solar_constant = 1367.
! Namelist for Sim_Phys
  logical:: diurnal_cycle    = .false.
  logical:: mixed_layer      = .false.
  logical:: gray_rad         = .false.
  logical:: strat_rad        = .false.
  logical:: do_abl = .false.
  logical:: do_mon_obkv = .true.
  logical:: do_lin_microphys = .false.
  logical:: do_k_warm_rain = .false.
  logical:: do_strat_forcing = .true.
  logical:: do_sg_conv       = .false.
  logical:: prog_cloud       = .true.
  logical:: zero_winds       = .false.  ! use this only for the doubly periodic domain
  logical:: do_reed_phys     = .false.
  logical:: reed_cond_only   = .false.
  real:: tau_zero  = 1.   ! time-scale to "zero" domain-averaged winds (days)
  real:: s_fac  =  0.1
  real:: heating_rate =  0.5    ! deg K per day, stratosphere, for gray-radiation
  real:: cooling_rate = -1.0    ! deg K per day
  real:: c0 = 1.E7              ! Ocean heat capabicity
  real:: low_c = 0.3            ! global mean *low* cloud fraction
  real:: sw_abs = 0.            ! fraction of the solar absorbed/reflected by the atm
  logical:: uniform_sst = .true.
  real:: sst0 = 300.    ! 302.15 for DCIMP case-51
  real:: shift_n = 0.
  logical:: do_t_strat = .false.
  real:: p_strat = 50.E2
  real:: t_strat = 200.
  real:: tau_strat = 10.     ! days
  real:: t_fac = 1.
  integer:: print_freq = 86400   ! sec

namelist /sim_phys_nml/mixed_layer, gray_rad, strat_rad, do_lin_microphys,   &
                       heating_rate, cooling_rate, uniform_sst, sst0, c0,    &
                       sw_abs, prog_cloud, low_c, do_sg_conv, diurnal_cycle, &
                       do_mon_obkv, do_t_strat, p_strat, t_strat, tau_strat, &
                       do_abl, t_fac, tau_difz, do_strat_forcing, s_fac,     &
                       shift_n, print_freq, zero_winds, tau_zero, tau_winds, &
                       tau_temp, tau_press, do_K_warm_rain,                  &
                       do_reed_phys, reed_cond_only


    integer :: id_rain, id_dqdt, id_dTdt, id_dudt, id_dvdt
    real             :: missing_value = -1.e10

contains

 subroutine fv_phys(npx, npy, npz, is, ie, js, je, ng, nq,  &
                    u, v, w, pt, q, pe, delp, peln, pkz, pdt,    &
                    ua, va, phis, grid, ptop, ak, bk, ks, ps, pk,&
                    u_srf, v_srf, ts, delz, hydrostatic,         &
                    oro, rayf, p_ref, fv_sg_adj,                 &
                    do_Held_Suarez, gridstruct, flagstruct,      &
                    neststruct, nwat, bd, domain, axes,          &
                    Time, time_total)

    integer, INTENT(IN   ) :: npx, npy, npz
    integer, INTENT(IN   ) :: is, ie, js, je, ng, nq, nwat
    integer, INTENT(IN   ) :: fv_sg_adj, axes(4)
    real, INTENT(IN) :: p_ref, ptop
    real, INTENT(IN) :: oro(is:ie,js:je)

    real   , INTENT(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
    real   , INTENT(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) ::    w(is-ng:       ,js-ng:       ,1: )
    real   , INTENT(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) ::    q(is-ng:ie+  ng,js-ng:je+  ng,npz, nq)
    real   , INTENT(INOUT) ::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
    real   , INTENT(INOUT) :: peln(is   :ie     ,1:npz+1,js   :je     )
    real   , INTENT(INOUT) ::  pkz(is   :ie     ,js   :je     ,1:npz)
    real   , INTENT(INOUT) ::  pk (is   :ie     ,js   :je     ,npz+1)
    real   , INTENT(INOUT) ::   ua(is-ng:ie+ng,js-ng:je+ng,npz)
    real   , INTENT(INOUT) ::   va(is-ng:ie+ng,js-ng:je+ng,npz)
    real   , INTENT(INOUT) ::   ps(is-ng:ie+ng,js-ng:je+ng)

    real   , INTENT(IN   ) :: phis(is-ng:ie+ng,js-ng:je+ng)
    real   , INTENT(IN   ) :: grid(is-ng:ie+ng,js-ng:je+ng, 2)
    real   , INTENT(IN   ) :: ak(npz+1), bk(npz+1)
    integer, INTENT(IN   ) :: ks

    real   , INTENT(IN   )    :: pdt
    logical, INTENT(IN   )    :: rayf, do_Held_Suarez
    real, INTENT(inout):: u_srf(is:ie,js:je)
    real, INTENT(inout):: v_srf(is:ie,js:je)
    real, INTENT(inout)::    ts(is:ie,js:je)

    type(fv_grid_type) :: gridstruct
    type(fv_flags_type) :: flagstruct
    type(fv_nest_type) :: neststruct
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(domain2d), intent(INOUT) :: domain

    type (time_type), intent(in) :: Time
    real, INTENT(IN), optional:: time_total
    logical, intent(in) ::  hydrostatic
    real, intent(inout) ::  delz(is-ng:     ,js-ng:     ,1: )
    real, allocatable:: u_dt(:,:,:), v_dt(:,:,:), t_dt(:,:,:), q_dt(:,:,:,:)
    real, allocatable :: rain(:,:)
    real, allocatable ::   pz(:,:), rdelp(:,:)
    real :: qdiag(1,1,1)
    logical moist_phys
    integer  isd, ied, jsd, jed
    integer  i, j, k, n
    real rdg
    logical :: first_step = .true.
    logical used

    master = is_master()
    if (.not. sim_phys_initialized) call fv_phys_init(nwat, axes, Time)

    isd = is-ng;   ied = ie + ng
    jsd = js-ng;   jed = je + ng

     allocate ( u_dt(isd:ied,jsd:jed,npz) )
     allocate ( v_dt(isd:ied,jsd:jed,npz) )
     allocate ( t_dt(is:ie,js:je,npz) )
     allocate ( q_dt(is:ie,js:je,npz,nq) )

! Place the memory in the optimal shared mem space!
!$omp parallel do default(shared)
     do k=1, npz
        do j=jsd, jed
           do i=isd, ied
              u_dt(i,j,k) = 0.
              v_dt(i,j,k) = 0.
           enddo
        enddo

        do j=js, je
           do i=is, ie
              t_dt(i,j,k) = 0.
           enddo
        enddo
        do n=1,nq
        do j=js, je
           do i=is, ie
              q_dt(i,j,k,n) = 0.
           enddo
        enddo
        enddo
     enddo

     if ( fv_sg_adj > 0 ) then
          call fv_dry_conv(isd, ied, jsd, jed, is, ie, js, je, npz, min(6,nq), pdt,  &
                           fv_sg_adj, delp, pe, peln, pkz, pt, q, ua, va,  &
                           hydrostatic, w, delz, u_dt, v_dt, t_dt, q_dt )
     endif

    if( do_Held_Suarez ) then
       moist_phys = .false.
       call Held_Suarez_Tend(npx, npy, npz, is, ie, js, je, ng, nq,   &
                              u, v, pt, q, pe, delp, peln, pkz, pdt,  &
                              ua, va, u_dt, v_dt, t_dt, q_dt, grid,   &
                              delz, phis, hydrostatic, ak, bk, ks,    &
                              do_strat_forcing, rayf, master, Time, time_total)

    elseif (do_reed_phys) then

       moist_phys = .true.
       if (.not. hydrostatic) rdg = -rdgas / grav
       if (.not. allocated(pz)) allocate(pz(is:ie,npz))
       if (.not. allocated(rdelp)) allocate(rdelp(is:ie,npz))
       sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
       zvir = rvgas/rdgas - 1.

       do j=js,je

          if (hydrostatic) then

             do k=1,npz
                do i=is,ie
                   !pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
                   pz(i,k) = (pe(i,k+1,j)-pe(i,k,j))/(peln(i,k+1,j)-peln(i,k,j))
                enddo
             enddo

          else
             do k=1,npz
             do i=is,ie
#ifdef USE_COND
                pz(i,k) = rdg*delp(i,j,k)*pt(i,j,k)*    &
                     (1.+zvir*q(i,j,k,sphum)-q_con(i,j,k))/delz(i,j,k) 
!!$                pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
!!$                     (1.+dp1(i,j,k)-q_con(i,j,k))/delz(i,j,k)) )
#else
                pz(i,k) = rdg*delp(i,j,k)*pt(i,j,k)*    &
                     (1.+zvir*q(i,j,k,sphum))/delz(i,j,k)
!!$                pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
!!$                     (1.+dp1(i,j,k))/delz(i,j,k)) )
#endif
             enddo
             enddo
                
          endif

          if (master .and. first_step) then
             do k=1,npz
                print*, pe(i,k,j), pz(is,k)
             enddo
             print*, pe(i,npz+1,j)
             first_step = .false.
          endif

          do k=1,npz
          do i=is,ie
             rdelp(i,k) = 1./delp(i,j,k)
          enddo
          enddo

          if (.not. allocated(rain)) allocate(rain(is:ie,js:je))
!!$          call reed_simple_physics_v5(ie-is+1, npz, pdt, gridstruct%agrid(is:ie,j,2), pt(is:ie,j,:), q(is:ie,j,:,sphum), &
!!$               ua(is:ie,j,:), va(is:ie,j,:), pz, pe(is:ie,:,j), delp(is:ie,j,:), &
!!$               rdelp, ps(is:ie,j), rain(is:ie,j), 0)
          call reed_simple_physics(ie-is+1, npz, pdt, pt(is:ie,j,:), q(is:ie,j,:,sphum), &
               ua(is:ie,j,:), va(is:ie,j,:), pz, pe(is:ie,:,j), delp(is:ie,j,:), &
               rdelp, pe(is:ie,npz+1,j), rain(is:ie,j), reed_cond_only, & 
               q_dt(is:ie,j,:,sphum), t_dt(is:ie,j,:), u_dt(is:ie,j,:), v_dt(is:ie,j,:))

       enddo

       if (id_rain > 0) then
          used=send_data(id_rain, rain, time)
       endif

    else
       moist_phys = .true.
                                             call timing_on('SIM_PHYS')
       call sim_phys(npx, npy, npz, is, ie, js, je, ng, nq, nwat, pk, pkz, &
                     u_dt, v_dt, t_dt, q_dt, u, v, w, ua, va, pt, delz, q, &
                     pe, delp, peln, ts, oro, hydrostatic, pdt, grid, ak, bk, &
                     p_ref, fv_sg_adj, Time, time_total, gridstruct)
                                            call timing_off('SIM_PHYS')
    endif

    if (id_dudt > 0) then
       used=send_data(id_dudt, u_dt(is:ie,js:je,npz), time)
    endif
    if (id_dvdt > 0) then
       used=send_data(id_dvdt, v_dt(is:ie,js:je,npz), time)
    endif
    if (id_dtdt > 0) then
       used=send_data(id_dtdt, t_dt(:,:,:), time)
    endif
    if (id_dqdt > 0) then
       used=send_data(id_dqdt, q_dt(:,:,:,sphum), time)
    endif

                        call timing_on('UPDATE_PHYS')
    call fv_update_phys (pdt, is, ie, js, je, isd, ied, jsd, jed, ng, nq,   &
                         u, v, w, delp, pt, q, qdiag, ua, va, ps, pe, peln, pk, pkz,  &
                         ak, bk, phis, u_srf, v_srf, ts,  &
                         delz, hydrostatic, u_dt, v_dt, t_dt, q_dt, &
                         moist_phys, Time, .false., gridstruct, &
                         gridstruct%agrid(:,:,1), gridstruct%agrid(:,:,2), &
                         npx, npy, npz, flagstruct, neststruct, bd, domain, ptop)
                        call timing_off('UPDATE_PHYS')
    deallocate ( u_dt )
    deallocate ( v_dt )
    deallocate ( t_dt )
    deallocate ( q_dt )

 end subroutine fv_phys


 subroutine sim_phys(npx, npy, npz, is, ie, js, je, ng, nq, nwat, pk, pkz,  &
                     u_dt, v_dt, t_dt, q_dt, u, v, w, ua, va,   &
                     pt, delz, q, pe, delp, peln, sst, oro, hydrostatic, &
                     pdt, agrid, ak, bk, p_ref, fv_sg_adj,    &
                     Time, time_total, gridstruct)
!-----------------------
! A simple moist physics
!-----------------------


 integer, INTENT(IN) :: npx, npy, npz
 integer, INTENT(IN) :: is, ie, js, je, ng, nq, nwat
 integer, INTENT(IN) :: fv_sg_adj
 real   , INTENT(IN) :: pdt
 real   , INTENT(IN) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
 real   , INTENT(IN) :: ak(npz+1), bk(npz+1)
 real, INTENT(IN):: p_ref
 real, INTENT(IN):: oro(is:ie,js:je)       ! land fraction
 logical, INTENT(IN):: hydrostatic

 type(time_type), intent(in) :: Time
 real, INTENT(IN), optional:: time_total

 real   , INTENT(INOUT) :: u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
 real   , INTENT(INOUT) :: v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
 real, INTENT(INOUT)::  pk (is:ie,js:je,npz+1)
 real, INTENT(INOUT)::  pkz(is:ie,js:je,npz)
 real, INTENT(INOUT)::   pt(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: delp(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT)::    q(is-ng:ie+ng,js-ng:je+ng,npz, nq)
 real, INTENT(INOUT)::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
 real, INTENT(INOUT):: peln(is  :ie   ,1:npz+1,js  :je  )
 real, INTENT(INOUT):: delz(is-ng :       ,js-ng  :       ,1: )
 real, intent(inout):: w(is-ng:     ,js-ng:     ,1: )
 real, INTENT(INOUT):: sst(is:ie,js:je)

 type(fv_grid_type) :: gridstruct

! Tendencies:
 real, INTENT(INOUT):: u_dt(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: v_dt(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
 real, INTENT(INOUT):: q_dt(is:ie,js:je,npz,nq)
 real, INTENT(INOUT):: ua(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: va(is-ng:ie+ng,js-ng:je+ng,npz)
! Local
 real, dimension(is:ie,js:je):: flux_t, flux_q, flux_u, flux_v, delm
 logical:: phys_hydrostatic = .true.
 real, parameter:: tice = 273.16
 real, parameter:: f1 = 2./3.
 real, parameter:: f2 = 1./3.
 real, dimension(is:ie,js:je,npz):: u3, v3, t3, p3, dz, zfull
 real, dimension(is:ie,js:je,npz+1):: zhalf, phalf
 real, dimension(is:ie,js:je,npz,nq):: q3
 real, dimension(is:ie,js:je):: rain, snow, ice, graup, land, mu
 real, dimension(is:ie,js:je):: ps, qs, rho, clouds
 real, dimension(is:ie,js:je):: olr, lwu, lwd, sw_surf
 real, dimension(is:ie,npz):: den
 real, dimension(is:ie,npz):: t_dt_rad
 real:: sday, rrg, tvm, olrm, swab, sstm, clds, deg2rad
 real:: rdt, tmp, cooling, heating,utmp, vtmp
 real:: fac_sm, rate_u, rate_v, rate_t, rate_q
 integer  i,j,k, km, iq, k_mp
 integer  isd, ied, jsd, jed
 integer  seconds, days
 logical print_diag

!  if (.not. sim_phys_initialized) call fv_phys_init(nwat)

   call get_time (time, seconds, days)

   if ( mod(seconds, print_freq)==0 ) then
        print_diag = .true.
   else
        print_diag = .false.
   endif

   km = npz
   isd = is-ng;   ied = ie + ng
   jsd = js-ng;   jed = je + ng

   zvir = rvgas/rdgas - 1.
! Factor for Small-Earth Approx.
   fac_sm = radius / 6371.0e3
   rrg  = rdgas / grav
   sday  = 24.*3600.*fac_sm
   rdt = 1. / pdt

   deg2rad = pi/180.

     sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
   liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
   ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
   cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')

   if ( nq.ge.6 ) then
        rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
        snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
        graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
   endif


   do k=1,km+1
      do j=js,je
         do i=is,ie
            phalf(i,j,k) = pe(i,k,j)
         enddo
      enddo
   enddo

! Compute zfull, zhalf
   do j=js,je
      do i=is,ie
         zhalf(i,j,npz+1) = 0.
         ps(i,j) = pe(i,km+1,j)
      enddo
   enddo

   do k=km,1,-1
      do j=js,je
         do i=is,ie
                     tvm = rrg*pt(i,j,k)*(1.+zvir*q(i,j,k,1))
               dz(i,j,k) = -tvm*(peln(i,k+1,j)-peln(i,k,j))
               p3(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
            zfull(i,j,k) = zhalf(i,j,k+1) + tvm*(1.-pe(i,k,j)/p3(i,j,k))
            zhalf(i,j,k) = zhalf(i,j,k+1) - dz(i,j,k)
         enddo
      enddo
   enddo

   do k=1, km
      do j=js,je
         do i=is,ie
            u3(i,j,k) = ua(i,j,k) + pdt*u_dt(i,j,k)
            v3(i,j,k) = va(i,j,k) + pdt*v_dt(i,j,k)
            t3(i,j,k) = pt(i,j,k) + pdt*t_dt(i,j,k)
         enddo 
      enddo 
   enddo

   do j=js,je
      do i=is,ie
         delm(i,j) =  delp(i,j,km)/grav
          rho(i,j) = -delm(i,j)/dz(i,j,km)
      enddo
   enddo
   
!----------
! Setup SST:
!----------

! Need to save sst in a restart file if mixed_layer = .T.
 if ( .not. mixed_layer ) then
   if ( uniform_sst ) then
        sst(:,:) = sst0
   else
      do j=js,je
         do i=is,ie
! Neale & Hoskins 2001
            if ( abs(agrid(i,j,2)-shift_n*deg2rad) <  pi/3. ) then
                 tmp = sin( 1.5*(agrid(i,j,2)-shift_n*deg2rad) ) ** 2
! control:
!                sst(i,j) = 273.16 + 27.*(1.-tmp)
! QOBS:
                 sst(i,j) = 273.16 + 27.*(1.-0.5*(tmp + tmp**2))
            else
                 sst(i,j) = 273.16
            endif
         enddo
      enddo
   endif
 endif


   do iq=1,nq
      do k=1,npz
         do j=js,je
            do i=is,ie
               q3(i,j,k,iq) = q(i,j,k,iq) + pdt*q_dt(i,j,k,iq)
            enddo
         enddo
      enddo
   enddo


!----------------------------
! Apply net radiative cooling
!----------------------------
 if ( gray_rad ) then
    if ( prog_cloud ) then
         call get_low_clouds( is,ie, js,je, km, q3(is,js,1,liq_wat), q3(is,js,1,ice_wat),   &
                                                q3(is,js,1,cld_amt), clouds )
    else
         clouds(:,:) = low_c
    endif

    do j=js,je
       do k=1, km
          do i=is,ie
             den(i,k) = -delp(i,j,k)/(grav*dz(i,j,k))
          enddo 
       enddo 
       call gray_radiation(seconds, is, ie, km, agrid(is,j,1), agrid(is,j,2), clouds(is,j), &
                           sst(is,j), pt(is:ie,j,1:km), ps(is,j), phalf(is:ie,j,1:km+1), &
                           dz(is:ie,j,1:km), den, t_dt_rad,              &
                           olr(is,j), lwu(is,j), lwd(is,j), sw_surf(is,j), ng)
       do k=1, km
          do i=is,ie
             t_dt(i,j,k) = t_dt(i,j,k) + t_dt_rad(i,k)
          enddo 
       enddo 
    enddo
    if ( print_diag ) then
         olrm = g0_sum(olr, is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
         swab = g0_sum(sw_surf, is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)

         if ( mixed_layer ) then
              sstm = g0_sum(sst, is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
              call prt_maxmin('TS', sst, is, ie, js, je, 0,  1, 1.0)
         endif

         if ( prog_cloud ) &
         clds = g0_sum(clouds, is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)

         if( master ) then
            write(*,*) 'Domain mean OLR =', olrm
            write(*,*) 'Domain mean SWA =', swab
            if ( mixed_layer ) &
            write(*,*) 'Domain mean SST (mixed-layer)=', sstm
            if ( prog_cloud ) &
            write(*,*) 'Domain mean low clouds fraction =', clds
         endif
    endif

 else

! Prescribed (non-interating) heating/cooling rate:
   rate_t = 1. / (tau_strat*86400. + pdt )
!  cooling = cooling_rate/sday


   do k=1, km
      do j=js, je
         do i=is, ie
            cooling = cooling_rate*(f1+f2*cos(agrid(i,j,2)))/sday
            if ( p3(i,j,k) >= 100.E2 ) then
               t_dt(i,j,k) = t_dt(i,j,k) + cooling*min(1.0, (p3(i,j,k)-100.E2)/200.E2)
            elseif ( do_t_strat ) then
               if ( p3(i,j,k) < p_strat) then
                    t_dt(i,j,k) = t_dt(i,j,k) + (t_strat - t3(i,j,k))*rate_t
               endif
            endif
         enddo 
      enddo 
   enddo      ! k-loop

 endif


 if ( strat_rad ) then
! Solar heating above 100 mb:
      heating = heating_rate / 86400.
      do k=1, km
         do j=js,je
            do i=is,ie
               if ( p3(i,j,k) < 100.E2 ) then
                  t_dt(i,j,k) = t_dt(i,j,k) + heating*(100.E2-p3(i,j,k))/100.E2
               endif
            enddo 
         enddo 
      enddo 
 endif


 do k=1, km
    do j=js,je
       do i=is,ie
          t3(i,j,k) = pt(i,j,k) + pdt*t_dt(i,j,k)
       enddo 
    enddo 
 enddo



!---------------
! Surface fluxes
!---------------

if( do_mon_obkv ) then

  call qsmith(ie-is+1, je-js+1, 1, sst, ps, q3(is:ie,js:je,km,sphum), qs)
  call qsmith(ie-is+1, je-js+1, 1, sst, ps, qs, qs) ! Iterate once

! Need to save ustar in a restart file (sim_phys)
  if ( .not. allocated(u_star) ) then
       allocate ( u_star(is:ie,js:je) )
!      u_star(:,:) = 1.E25  ! large enough to cause blowup if used
       u_star(:,:) = 1.E-3
  endif

                                                                             call timing_on('mon_obkv')
  call mon_obkv(zvir, ps, t3(is:ie,js:je, km), zfull(is:ie,js:je,km),     &
                rho, p3(is:ie,js:je,km), u3(is:ie,js:je,km), v3(is:ie,js:je,km), sst, &
                qs, q3(is:ie,js:je,km,sphum), flux_t, flux_q, flux_u, flux_v, u_star, &
                delm, pdt, mu, t_fac, master)
                                                                             call timing_off('mon_obkv')
!---------------------------------------------------
! delp/grav = delm = kg/m**2
! watts = J/s = N*m/s = kg * m**2 / s**3
! CP = J/kg/deg
! flux_t = w/m**2 = J / (s*m**2)
!---------------------------------------------------
   do j=js,je
      do i=is,ie
               rate_u = flux_u(i,j)/delm(i,j)
           u3(i,j,km) = u3(i,j,km) + pdt*rate_u
         u_dt(i,j,km) = u_dt(i,j,km) + rate_u
               rate_v = flux_v(i,j)/delm(i,j)
           v3(i,j,km) =   v3(i,j,km) + pdt*rate_v
         v_dt(i,j,km) = v_dt(i,j,km) + rate_v
               rate_t = flux_t(i,j)/(cp_air*delm(i,j))
         t_dt(i,j,km) = t_dt(i,j,km) + rate_t
           t3(i,j,km) =   t3(i,j,km) + rate_t*pdt
               rate_q = flux_q(i,j)/delm(i,j)
         q_dt(i,j,km,sphum) = q_dt(i,j,km,sphum) + rate_q
           q3(i,j,km,sphum) =   q3(i,j,km,sphum) + rate_q*pdt
      enddo
   enddo
endif

if ( zero_winds ) then

! Zero out (doubly periodic) domain averaged winds with time-scale tau_zero:
  rate_u = min(1., 1./(tau_zero*86400.))
!$omp parallel do default(shared) private(utmp, vtmp)
 do k=1, km
    utmp = g0_sum(u3(is,js,k), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
    vtmp = g0_sum(v3(is,js,k), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
    if ( master ) then
        if( sqrt(utmp**2+vtmp**2) > 1. ) then 
            write(*,*) k, 'Domain avg winds=', utmp, vtmp
        endif
    endif
    do j=js,je
       do i=is,ie
          u_dt(i,j,k) = u_dt(i,j,k) - utmp * rate_u
          v_dt(i,j,k) = v_dt(i,j,k) - vtmp * rate_u
       enddo 
    enddo 
 enddo

endif

 if ( do_abl ) then
! ABL: vertical diffusion:
! do j=js,je
!    do i=is,ie
!       mu(i,j) = -sqrt(area(i,j)/area_ref)*dz(i,j,km)/tau_difz
!    enddo
! enddo

  call pbl_diff(hydrostatic, pdt, is, ie, js, je, ng, km, nq, u3, v3, t3,      &
                w, q3, delp, p3, pe, sst, mu, dz, u_dt, v_dt, t_dt, q_dt, gridstruct%area, print_diag )
 endif

  if ( mixed_layer .and. gray_rad ) then
     do j=js, je
        do i=is, ie
           sst(i,j) = sst(i,j)+pdt*(sw_surf(i,j)+lwd(i,j)-flux_t(i,j)-hlv*flux_q(i,j)-lwu(i,j))/c0
        enddo 
     enddo
  endif

! if ( do_sg_conv .and. fv_sg_adj > 0 ) then
!       call sg_conv(is, ie, js, je, isd, ied, jsd, jed, km, nq, pdt, fv_sg_adj,   &
!                    delp, phalf, p3, zfull, zhalf, t3, q3, u3, v3, w, &
!                    u_dt, v_dt, t_dt, q_dt, sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, &
!                    hydrostatic, phys_hydrostatic)
! endif

  if ( do_lin_microphys ) then
      land(:,:) = 0.
      k_mp = 1
      do k=2, km
         tmp = ak(k) + bk(k)*1000.E2
         if ( tmp > 20.E2 ) then
              k_mp = k
              exit
         endif
      enddo
      if(master .and. print_diag) write(*,*) 'k_mp=', k_mp, tmp*0.01
      k_mp = 1

                                                                                          call timing_on('lin_cld_mp')
      call lin_cld_microphys_driver(q3(is:ie,js:je,1:km,sphum), q3(is:ie,js:je,1:km,liq_wat), q3(is:ie,js:je,1:km,rainwat),  &
                                    q3(is:ie,js:je,1:km,ice_wat), q3(is:ie,js:je,1:km,snowwat), q3(is:ie,js:je,1:km,graupel),  &
                                    q3(is:ie,js:je,1:km,cld_amt), q3(is:ie,js:je,1:km,cld_amt),   &
                                  q_dt(is:ie,js:je,1:km,sphum), q_dt(is:ie,js:je,1:km,liq_wat), q_dt(is:ie,js:je,1:km,rainwat), &
                                  q_dt(is:ie,js:je,1:km,ice_wat), q_dt(is:ie,js:je,1:km,snowwat), q_dt(is:ie,js:je,1:km,graupel), &
                                  q_dt(is:ie,js:je,1:km,cld_amt), t_dt, t3, w(is:,js:,1:), u3, v3,     &
                                  u_dt(is:ie,js:je,1:km),v_dt(is:ie,js:je,1:km), dz,      &
                                  delp(is:ie,js:je,1:km), gridstruct%area(is:ie,js:je),  &
                                  pdt, land, rain, snow, ice, graup, hydrostatic, phys_hydrostatic, &
                                  1,ie-is+1, 1,je-js+1, 1,km, k_mp,npz, Time )
                                                                                          call timing_off('lin_cld_mp')
  elseif ( do_K_warm_rain ) then
! This is a stand alone routine (modified from NCAR code). It updates q delp, and pt directly.
! Do not use other components in this module with this code (because NO tendencies)
      call K_warm_rain(pdt, is, ie, js, je, ng, km, zvir, sphum, liq_wat, rainwat, q, pt, q_dt, t_dt,   &
                       delp, delz, rain)
  endif

 end subroutine sim_phys


 subroutine pbl_diff(hydrostatic, dt, is, ie, js, je, ng, km, nq, ua, va, &
                     ta, w, q, delp, pm,  pe, ts, mu, dz, udt, vdt, tdt, qdt, area, print_diag )
 logical, intent(in):: hydrostatic
 integer, intent(in):: is, ie, js, je, ng, km, nq
 real, intent(in):: dt
 real, intent(in), dimension(is:ie,js:je):: ts, mu
 real, intent(in), dimension(is:ie,js:je,km):: dz, pm
 real, intent(in)::   pe(is-1:ie+1 ,1:km+1,js-1:je+1)
 real, intent(inout), dimension(is:ie,js:je,km):: ua, va, ta, tdt
 real, intent(inout), dimension(is-ng:     ,js-ng:     ,1:):: w
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km):: udt, vdt, delp
 real, intent(inout), dimension(is:ie,js:je,km,nq):: q, qdt
 logical, intent(in):: print_diag
 real, intent(in) :: area(is-ng:ie+ng,js-ng:je+ng)
! Local:
 real, dimension(is:ie,js:je):: pblh
 real, dimension(is:ie,km+1):: gh
 real, dimension(is:ie,km):: nu, a, f, r, q2, den
 real, dimension(is:ie,km):: gz, hd, te
 real, dimension(is:ie):: kz
 real, parameter:: mu_min = 1.E-5    ! ~ molecular viscosity
 real, parameter:: ustar2 = 1.E-4
 integer:: n, i, j, k
 real:: cv, rcv, rdt, tmp, tvm, tv_surf, surf_h, rin

 cv = cp_air - rdgas
 rcv = 1./cv
 rdt = 1./dt

! Put opnmp loop here
 do 1000 j=js, je

    do i=is, ie
       gh(i,km+1) = 0.
       pblh(i,j) = 0.
    enddo
    do k=km, 1, -1
       do i=is, ie
          gh(i,k) = gh(i,k+1) - dz(i,j,k)
       enddo
    enddo

! Locate PBL top (m)
    do 125 i=is, ie
       tv_surf = ts(i,j)*(1.+zvir*q(i,j,km,sphum))
       do k=km, 3, -1 
          tvm = ta(i,j,k)*(1.+zvir*q(i,j,k,sphum)-q(i,j,k,liq_wat)-    &
                q(i,j,k,ice_wat)-q(i,j,k,snowwat)-q(i,j,k,rainwat)-q(i,j,k,graupel))
!         tvm = ta(i,j,k)*(1.+zvir*q(i,j,k,sphum)-q(i,j,k,liq_wat)-q(i,j,k,ice_wat))
          tvm = tvm*(pe(i,km+1,j)/pm(i,j,k))**kappa
          rin = grav*(gh(i,k+1)-0.5*dz(i,j,k))*(tvm-tv_surf) / ( 0.5*(tv_surf+tvm)*  &
                (ua(i,j,k)**2+va(i,j,k)**2+ustar2) )
          if ( rin > 1. ) then
               pblh(i,j) = gh(i,k+1)
               goto 125
          endif
       enddo
125 continue

  do k=km, 1, -1
     do i=is, ie
        if ( gh(i,k)>6.E3 .or. (pblh(i,j) < -0.5*dz(i,j,km)) ) then
             nu(i,k) = mu_min
        else
            kz(i) = 0.5*mu(i,j)*gh(i,km)
            surf_h = s_fac*pblh(i,j)
            if ( gh(i,k) <= surf_h) then
                 kz(i) = mu(i,j)*gh(i,k)
                 nu(i,k) = kz(i)
            elseif (gh(i,k) <= pblh(i,j)) then
! Use Dargan's form:
                 nu(i,k) = kz(i)*gh(i,k)/surf_h*(1.-(gh(i,k)-surf_h)/(pblh(i,j)-surf_h))**2
            else
                 nu(i,k) = mu_min
            endif
        endif
     enddo
  enddo

 do k=1,km-1
    do i=is, ie
       a(i,k) = -nu(i,k) / ( 0.5*(dz(i,j,k)+dz(i,j,k+1)) )
    enddo
 enddo
 do i=is, ie
    a(i,km) = 0.
 enddo

 if ( .not. hydrostatic ) then
    do k=1, km
       do i=is,ie
          gz(i,k) = grav*(gh(i,k+1) - 0.5*dz(i,j,k))
             tmp  = gz(i,k) + 0.5*(ua(i,j,k)**2+va(i,j,k)**2+w(i,j,k)**2)
             tvm  = ta(i,j,k)*(1.+zvir*q(i,j,k,sphum)) 
          hd(i,k) = cp_air*tvm + tmp
          te(i,k) =     cv*tvm + tmp
       enddo
    enddo
 endif

! u-diffusion
    do k=1,km
       do i=is, ie
          den(i,k) = delp(i,j,k)/dz(i,j,k)
           q2(i,k) = ua(i,j,k)*den(i,k)
            f(i,k) = -dz(i,j,k)*rdt
            r(i,k) = -f(i,k)*q2(i,k)
       enddo
    enddo
    call trid_dif2(a, f, r, q2, is, ie, km)
    do k=1,km
       do i=is, ie
           q2(i,k) = q2(i,k) / den(i,k)
          udt(i,j,k) = udt(i,j,k) + rdt*(q2(i,k) - ua(i,j,k))
           ua(i,j,k) = q2(i,k)
       enddo
    enddo
!---
! v-diffusion
    do k=1,km
       do i=is, ie
           q2(i,k) = va(i,j,k)*den(i,k)
            r(i,k) = -f(i,k)*q2(i,k)
       enddo
    enddo
    call trid_dif2(a, f, r, q2, is, ie, km)
    do k=1,km
       do i=is, ie
             q2(i,k) = q2(i,k) / den(i,k)
          vdt(i,j,k) = vdt(i,j,k) + rdt*(q2(i,k) - va(i,j,k))
           va(i,j,k) = q2(i,k)
       enddo
    enddo
!---
 if ( .not. hydrostatic ) then
! w-diffusion
    do k=1,km
       do i=is, ie
           q2(i,k) = w(i,j,k) * den(i,k)
            r(i,k) = -f(i,k)*q2(i,k)
       enddo
    enddo
    call trid_dif2(a, f, r, q2, is, ie, km)
    do k=1,km
       do i=is, ie
          w(i,j,k) = q2(i,k) / den(i,k)
       enddo
    enddo
 endif
!--- micro-physics
 do n=1, nq
! if ( (n.ne.rainwat) .and. (n.ne.snowwat) .and. (n.ne.graupel) .and. (n.ne.cld_amt) ) then
  if ( n.ne.cld_amt ) then
    do k=1,km
       do i=is, ie
          q2(i,k) = q(i,j,k,n)*den(i,k)
           r(i,k) = -f(i,k)*q2(i,k)
       enddo
    enddo
    call trid_dif2(a, f, r, q2, is, ie, km)
    do k=1,km
       do i=is, ie
          q2(i,k) = q2(i,k) / den(i,k)
          qdt(i,j,k,n) = qdt(i,j,k,n) + rdt*(q2(i,k) - q(i,j,k,n))
            q(i,j,k,n) = q2(i,k)
       enddo
    enddo
  endif
 enddo

! Diffusion of dry static energy
 if ( .not. hydrostatic ) then
    do k=1,km
       do i=is, ie
           q2(i,k) = hd(i,k)*den(i,k)
            r(i,k) = -f(i,k)*q2(i,k)
       enddo
    enddo
    call trid_dif2(a, f, r, q2, is, ie, km)
    do k=1,km
       do i=is, ie
          te(i,k) = te(i,k) + q2(i,k)/den(i,k) - hd(i,k)
          q2(i,k) = rcv*(te(i,k)-gz(i,k)-0.5*(ua(i,j,k)**2+va(i,j,k)**2+w(i,j,k)**2)) &
                  / (1.+zvir*q(i,j,k,sphum))
          tdt(i,j,k) = tdt(i,j,k) + rdt*(q2(i,k) - ta(i,j,k))
           ta(i,j,k) = q2(i,k)
       enddo
    enddo
 endif

1000 continue

 if ( print_diag ) then
     tmp = g0_sum(pblh, is, ie, js, je, 0, area(is:ie,js:je), 1)
     if (master) write(*,*) 'Mean PBL H (km)=', tmp*0.001
     call prt_maxmin('PBLH(km)', pblh, is, ie, js, je, 0,  1, 0.001)
!    call prt_maxmin('K_ABL (m^2/s)', mu, is, ie, js, je, 0,  1, 1.)
 endif

 end subroutine pbl_diff

 subroutine trid_dif2(a, v, r, q, i1, i2, km)
 integer, intent(in):: i1, i2
 integer, intent(in):: km      ! vertical dimension
 real, intent(in), dimension(i1:i2,km):: a, v, r
 real, intent(out),dimension(i1:i2,km):: q
! Local:
 real:: gam(i1:i2,km)
 real:: bet(i1:i2)
 integer:: i, k

! Zero diffusive fluxes at top and bottom:
! top: k=1
  do i=i1,i2
     bet(i) = -(v(i,1) + a(i,1))
     q(i,1) = r(i,1) / bet(i)   
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) =  a(i,k-1) / bet(i)
          bet(i) = -( a(i,k-1)+v(i,k)+a(i,k) + a(i,k-1)*gam(i,k))
! a(i,km) = 0
          q(i,k) = ( r(i,k) - a(i,k-1)*q(i,k-1) ) / bet(i)
     enddo   
  enddo       
               
  do k=km-1,1,-1
     do i=i1,i2  
        q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
     enddo                                                                                              
  enddo                                                                                                 
                                                                                                        
 end subroutine trid_dif2


 subroutine fv_nudge( npz, is, ie, js, je, ng, u, v, w, delz, delp, pt, dt, hydrostatic)
!
! Nudge the prognostic varaibles toward the IC
! This is only useful for generating balanced steady state IC

    real   , INTENT(IN   ) :: dt
    integer, INTENT(IN   ) :: npz, is, ie, js, je, ng
    logical, INTENT(IN   ) :: hydrostatic
    real   , INTENT(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
    real   , INTENT(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) ::    w(is-ng:       ,js-ng:       ,1: )
    real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) :: delz(is-ng:       ,js-ng:       ,1: )
    real   , INTENT(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
    real c_w, c_p, c_t
    integer  i, j, k

    if ( .not. nudge_initialized ) then

         if( tau_winds > 0. ) then
             allocate ( u0(is:ie,  js:je+1,npz) )
             allocate ( v0(is:ie+1,js:je  ,npz) )
             do k=1,npz
                do j=js,je+1
                   do i=is,ie
                      u0(i,j,k) = u(i,j,k) 
                   enddo
                enddo
                do j=js,je
                   do i=is,ie+1
                      v0(i,j,k) = v(i,j,k)
                   enddo
                enddo
             enddo
         endif

         if( tau_press > 0. ) then
             allocate ( dp(is:ie,js:je,npz) )
             do k=1,npz
                do j=js,je
                   do i=is,ie
                      dp(i,j,k) = delp(i,j,k)
                   enddo
                enddo
             enddo
         endif

         if( tau_temp > 0. ) then
             allocate ( t0(is:ie,js:je,npz) )
             do k=1,npz
                do j=js,je
                   do i=is,ie
                      t0(i,j,k) = pt(i,j,k)
                   enddo
                enddo
             enddo
         endif

         nudge_initialized = .true.
         if ( is_master() ) write(*,*) 'Nudging of IC initialized.'
         return
    endif

! Nudge winds to initial condition:

      do k=1,npz
         if( tau_winds > 0. ) then
            c_w = dt/tau_winds
            do j=js,je+1
               do i=is,ie
                  u(i,j,k) = (u(i,j,k)+c_w*u0(i,j,k)) / (1.+c_w)
               enddo
            enddo
            do j=js,je
               do i=is,ie+1
                  v(i,j,k) = (v(i,j,k)+c_w*v0(i,j,k)) / (1.+c_w)
               enddo
            enddo
            if ( .not.hydrostatic ) then
            do j=js,je
               do i=is,ie+1
                  w(i,j,k) = w(i,j,k) / (1.+c_w)
               enddo
            enddo
            endif
         endif
         if ( tau_temp > 0. ) then
             c_t = dt/tau_temp
             do j=js,je
                do i=is,ie
                   pt(i,j,k) = (pt(i,j,k)+c_t*t0(i,j,k)) / (1.+c_t)
                enddo
             enddo
             if ( .not.hydrostatic ) then
! delz
             endif
         endif

         if ( tau_press > 0. ) then
             c_p = dt/tau_press
             do j=js,je
                do i=is,ie
                   delp(i,j,k) = (delp(i,j,k)+c_p*dp(i,j,k)) / (1.+c_p)
                enddo
             enddo
         endif
      enddo

 end subroutine fv_nudge

 subroutine gray_radiation(sec, is, ie, km, lon, lat, clouds, ts, temp, ps, phalf, &
                           delz, rho, t_dt, olr, lwu, lwd, sw_surf, ng)

! Gray-Radiation algorithms based on Frierson, Held, and Zurita-Gotor, 2006 JAS
! Note: delz is negative
! Coded by S.-J. Lin, June 20, 2012
      integer, intent(in):: sec
      integer, intent(in):: is, ie, km, ng
      real, dimension(is:ie):: ts
      real, intent(in), dimension(is:ie):: lon, lat
      real, intent(in), dimension(is-ng:ie+ng,km):: temp, phalf, rho
      real, intent(in), dimension(is-ng:     ,1:):: delz
      real, intent(in), dimension(is:ie):: ps, clouds
      real, intent(out), dimension(is:ie,km):: t_dt
      real, intent(out), dimension(is:ie):: olr, lwu, lwd, sw_surf
! local:
      real, dimension(is:ie,km+1):: ur, dr
      real, dimension(is:ie,km+1):: tau, lw
      real, dimension(is:ie,km):: delt, b
      real, dimension(is:ie):: tau0
      real, parameter:: t0e =  8.    ! Dargan value= 6.
      real, parameter:: t0p = 1.5
      real, parameter::   fl = 0.1
      real, parameter:: sbc = 5.6734E-8
      real, parameter::  cp = 1004.64
      real:: rs0 = 938.4            ! W m-2
      real:: dels  = 1.4
      real:: sig, sw_rad, solar_ang
      integer:: i, k

      do i=is, ie
         tau0(i) = t0e + (t0p-t0e) * sin(lat(i))**2
      enddo
! Annual & global mean solar_abs ~ 0.25 * 1367 * ( 1-0.3) ~ 240
! Earth cross section/total_area = 0.25; 0.3 = Net cloud reflection and atm abs

      if ( diurnal_cycle ) then
           sw_rad = solar_constant*(1. - sw_abs) 
           do i=is, ie
               solar_ang = 2*pi*real(sec)/86400. + lon(i)
              sw_surf(i) = sw_rad*(1.-clouds(i))*cos(lat(i))*max(0.,cos(solar_ang))
           enddo
      else
           sw_rad = (1./pi) * solar_constant*(1. - sw_abs)
           do i=is, ie
              sw_surf(i) = sw_rad*(1.-clouds(i))*cos(lat(i))
           enddo
      endif

      do k=1, km+1
         do i=is, ie
#ifndef STRAT_OFF
! Dargan version:
                 sig = phalf(i,k)/ps(i) 
            tau(i,k) = tau0(i)*( sig*fl + (1.-fl)*sig**4 )
#else
! SJL: less cooling for the stratosphere
            tau(i,k) = tau0(i) * (phalf(i,k)/ps(i))**4
#endif
         enddo
      enddo
      do k=1, km
         do i=is, ie
            delt(i,k) = tau(i,k+1) - tau(i,k)
               b(i,k) =  sbc*temp(i,k)**4
         enddo
      enddo

! top down integration:
      do i=is, ie
         dr(i,1) = 0.
      enddo
      do k=1, km
         do i=is, ie
            dr(i,k+1) = (dr(i,k)+delt(i,k)*(b(i,k)-0.5*dr(i,k)))/(1.+0.5*delt(i,k)) 
         enddo
      enddo
! Bottom up
      do i=is, ie
         ur(i,km+1) = sbc*ts(i)**4
             lwu(i) = ur(i,km+1)
      enddo
      do k=km, 1, -1
         do i=is, ie
            ur(i,k) = (ur(i,k+1)+delt(i,k)*(b(i,k)-0.5*ur(i,k+1)))/(1.+0.5*delt(i,k)) 
         enddo
      enddo

! Compute net long wave cooling rate:
      do k=1, km+1
         do i=is, ie
            lw(i,k) = ur(i,k) - dr(i,k)
         enddo
      enddo
      do k=1, km
         do i=is, ie
            t_dt(i,k) = (lw(i,k) - lw(i,k+1))/(cp*rho(i,k)*delz(i,k))
         enddo
      enddo

      do i=is, ie
         olr(i) = ur(i,1)
         lwd(i) = dr(i,km+1)
      enddo

 end subroutine gray_radiation



 subroutine get_low_clouds( is,ie, js,je, km, ql, qi, qa, clouds )
 integer, intent(in):: is,ie, js,je, km
 real, intent(in), dimension(is:ie,js:je,km):: ql, qi, qa 
 real, intent(out), dimension(is:ie,js:je):: clouds
 integer:: i, j, k

 do j=js, je
    do i=is, ie
       clouds(i,j) = 0.
       do k=2,km
          if ((ql(i,j,k)>1.E-6 .or. qi(i,j,k)>1.e-5) .and.  qa(i,j,k)>1.E-4 ) then
! Maximum overlap
               clouds(i,j) = max(clouds(i,j), qa(i,j,k))
          endif
       enddo
       clouds(i,j) = min( 1., clouds(i,j) )
    enddo
 enddo

 end subroutine get_low_clouds

 subroutine fv_phys_init(nwat, axes, time)

 integer, intent(IN) :: nwat, axes(4)
 type (time_type), intent(in) :: Time
 integer :: unit, ierr, io

!   ----- read and write namelist -----
    if ( file_exist('input.nml')) then
         unit = open_namelist_file ('input.nml')
         read  (unit, nml=sim_phys_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'sim_phys_nml')
 10     call close_file (unit)
    endif

    if (nwat /= 6 .and. do_lin_microphys) then
       call mpp_error(FATAL, "Need nwat == 6 to run Lin Microphysics.")
    endif

    if (do_reed_phys) then
       id_rain = register_diag_field ( 'dynamics', 'rain', axes(1:2), time,        &
            'rainfall_rate', 'mm H20/d', missing_value=missing_value )
    endif
    id_dqdt    = register_diag_field ( 'dynamics', 'dqdt', axes(1:3), time, &
         'humidity tendency', 'kg/kg/s', missing_value=missing_value )
    id_dTdt    = register_diag_field ( 'dynamics', 'dtdt', axes(1:3), time, &
         'temperature tendency', 'K/s', missing_value=missing_value )
    id_dudt    = register_diag_field ( 'dynamics', 'dudt', axes(1:3), time, &
         'zonal wind tendency', 'm/s/s', missing_value=missing_value )
    id_dvdt    = register_diag_field ( 'dynamics', 'dvdt', axes(1:3), time, &
         'meridional wind tendency', 'm/s/s', missing_value=missing_value )
    !To do: surface fluxes, surface winds, SST

    sim_phys_initialized = .true.

 end subroutine fv_phys_init


 real function g0_sum(p, ifirst, ilast, jfirst, jlast, ngc, area, mode)
 use mpp_mod,           only: mpp_sum
      real, save :: global_area

! Fast version of globalsum
      integer, intent(IN) :: ifirst, ilast
      integer, intent(IN) :: jfirst, jlast, ngc
      integer, intent(IN) :: mode  ! if ==1 divided by area
      real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
      real, intent(IN) :: area(ifirst-ngc:ilast+ngc,jfirst-ngc:jlast+ngc)
      integer :: i,j
      real gsum

!-------------------------
! Quick local sum algorithm
!-------------------------
      if ( .not. g0_sum_initialized ) then
         allocate (l_area(ifirst:ilast,jfirst:jlast))
         global_area = 0.
         do j=jfirst,jlast
           do i=ifirst,ilast
             global_area = global_area + area(i,j)
             l_area(i,j) = area(i,j)
           enddo
         enddo
         call mpp_sum(global_area)
!        if ( mpp_pe().eq.mpp_root_pe() ) write(*,*) 'Global Area=',global_area
         g0_sum_initialized = .true.
      end if

      gsum = 0.
      do j=jfirst,jlast
        do i=ifirst,ilast
          gsum = gsum + p(i,j)*l_area(i,j)
        enddo
      enddo
      call mpp_sum(gsum)

      if ( mode==1 ) then
        g0_sum = gsum / global_area
      else
        g0_sum = gsum
      endif

! Make it reproducible by truncating to 32-bit precision:
      g0_sum = real(g0_sum, 4)

 end function g0_sum

 subroutine K_warm_rain(dt, is, ie, js, je, ng, km, zvir, sphum, liq_wat, rainwat, q, pt, q_dt, t_dt, dp, delz, rain)
 real, intent(in):: dt ! time step
 real, intent(in):: zvir
 integer, intent(in):: is, ie, js, je, km, ng
 integer, intent(in):: sphum, liq_wat, rainwat
 real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng, km):: delz
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km):: dp, pt
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng,km,*):: q
 real, intent(inout),dimension(is:ie,js:je, km,*):: q_dt
 real, intent(inout),dimension(is:ie,js:je, km):: t_dt
 real, intent(out), dimension(is:ie,js:je):: rain
! Local:
 real, dimension(km):: t1, q1, q2, q3, zm, rho, drym, dm
 real, dimension(km):: qa, qb, qc
 real:: ze(km+1)
 real:: rdt
 integer i,j,k

 rdt = 1. / dt
!$omp parallel do default(shared) private(qa,qb,qc,t1, ze, dm, q1, q2, q3, zm, rho, drym)
 do j=js, je
    do i=is, ie
! Compute zm, & rho
       ze(km+1) = 0.
       do k=km, 1, -1
          ze(k) = ze(k+1) - delz(i,j,k)
       enddo
       do k=1,km
          zm(k)  = 0.5*(ze(k)+ze(k+1))
! Moist air mass:
          dm(k) = dp(i,j,k)
! Tracer mass:
          qa(k) = q(i,j,k,sphum)
          qb(k) = q(i,j,k,liq_wat)
          qc(k) = q(i,j,k,rainwat)
          q1(k) = qa(k) * dm(k)
          q2(k) = qb(k) * dm(k)
          q3(k) = qc(k) * dm(k)
! Dry air mass per unit area
          drym(k) = dm(k) - (q1(k)+q2(k)+q3(k))
! Dry air density
           rho(k) = -drym(k)/(grav*delz(i,j,k))
! Convert to dry mixing ratios:
          q1(k) = q1(k) / drym(k)
          q2(k) = q2(k) / drym(k)
          q3(k) = q3(k) / drym(k)
          t1(k) = pt(i,j,k)
       enddo

       call KESSLER1( t1, q1, q2, q3, rho, dt, zm, km, rain(i,j) )

! Convert back to moist mixing ratios
! Update 3D fields:
       do k=1,km
!         t_dt(i,j,k) = t_dt(i,j,k) + (t1(k) - pt(i,j,k))/dt
          pt(i,j,k) = t1(k)
! Tracer mass:
          q1(k) = q1(k) * drym(k)
          q2(k) = q2(k) * drym(k)
          q3(k) = q3(k) * drym(k)
! Moist air mass:
#ifdef LOCAL_UPDATE
          dp(i,j,k) = drym(k) + q1(k)+q2(k)+q3(k)
          q(i,j,k,  sphum) = q1(k) / dp(i,j,k)
          q(i,j,k,liq_wat) = q2(k) / dp(i,j,k)
          q(i,j,k,rainwat) = q3(k) / dp(i,j,k)
#else
! This allows sedimentation transport to be done within fv_update_phys
          t_dt(i,j,k) = 0.
          q_dt(i,j,k,  sphum) = (q1(k)/dm(k) - qa(k))*rdt
          q_dt(i,j,k,liq_wat) = (q2(k)/dm(k) - qb(k))*rdt
          q_dt(i,j,k,rainwat) = (q3(k)/dm(k) - qc(k))*rdt
#endif
       enddo
    enddo
 enddo

 end subroutine K_warm_rain

 SUBROUTINE KESSLER1( T, QV, QC, QR, R, dt, Z, NZ, precip )
! T  - TEMPERATURE (K)
! QV - WATER VAPOR MIXING RATIO (GM/GM)
! QC - CLOUD WATER MIXING RATIO (GM/GM)
! QR - RAIN WATER MIXING RATIO (GM/GM)
! R  - DRY AIR DENSITY (GM/M^3)
! dt - TIME STEP (S)
! Z - HEIGHTS OF THERMODYNAMIC LEVELS IN THE GRID COLUMN (M)
! NZ - NUMBER OF THERMODYNAMIC LEVELS IN THE COLUMN
! VARIABLES IN THE GRID COLUMN ARE ORDERED FROM THE SURFACE TO THE TOP.
! k=1 is the top layer
! Dry mixing ratios?
! OUTPUT VARIABLES:
   integer, intent(in):: nz
   real, intent(in):: dt
   REAL, intent(in)   :: R(NZ), Z(NZ)
   REAL, intent(inout):: T(NZ), QV(NZ), QC(NZ), QR(NZ)
   real, intent(out):: precip
! Local:
   real, parameter:: c_k = 1003.    ! heat capacity of the K scheme
   real, parameter:: F2X = 17.27
   real, parameter:: F5 = 237.3*F2X*2.5E6/1003.
   real, parameter:: XK = .2875
   real, parameter:: cvd = 717.56 
   REAL RHALF(NZ), VELQR(NZ), SED(NZ), pc(NZ)
   REAL ERN, QRPROD, PROD, QVS, cvm
   INTEGER K

   DO K=1,NZ
      RHALF(K) = SQRT(R(1)/R(K))
      pc(K) = 3.8e2 / (R(k)*rdgas*T(k))
! Liquid water terminal velocity (m/s) following K&W78 Eq. 2.15
      VELQR(K) = 36.34*(QR(K)*R(K))**0.1364*RHALF(K)
   END DO

! Sedimentation term using upstream differencing
!
   SED(1) = -dt*QR(1)*VELQR(1)/(.5*(Z(1)-Z(2)))
   DO K=2,NZ
      SED(K) = dt*(R(K-1)*QR(K-1)*VELQR(K-1)         &
              -R(K)*QR(K)*VELQR(K ))/(R(K)*(Z(K-1)-Z(K)))
   END DO
   precip = SED(nz)

   DO K=1,NZ
! Autoconversion and accretion rates following K&W78 Eq. 2.13a,b
! Threshold = 1.e-3
      QRPROD = QC(K) - (QC(K)-dt*max(.001*(QC(K)-.001), 0.))    &
             / (1.+dt*2.2*QR(K)**.875)
      QC(K) = max(QC(K)-QRPROD,0.)
      QR(K) = max(QR(K)+QRPROD+SED(K),0.)

! Saturation vapor mixing ratio (gm/gm) following K&W78 Eq. 2.11
      QVS = pc(K)*EXP(F2X*(T(K)-273.) /(T(K)- 36.))

      PROD = (QV(K)-QVS)/(1.+QVS*F5/(T(K)-36.)**2)
! Evaporation rate following K&W78 Eq3. 3.8-3.10
      ERN = min(dt*(((1.6+124.9*(R(K)*QR(K))**.2046)   &
          *(R(K)*QR(K))**.525)/(2.55E6*pc(K)             &
          /(3.8 *QVS)+5.4E5))*(DIM(QVS,QV(K))            &
          /(R(K)*QVS)),max(-PROD-QC(K),0.),QR(K))
! NICAM way:
!!!      cvm = (1.-qv(k)-qc(k))*cvd + qv(k)*cv_vap + qc(k)*c_liq
!!!     T(k) = T(k) + (Lv0+dc_vap*T(k))*(max(PROD,-QC(k))-ERN) / cvm
! Saturation adjustment following K&W78 Eq.2.14a,b
! NCAR Way:
        T(K) = T(K) + 2.5E6*(max(PROD,-QC(K))-ERN) / c_k
       QV(K) = max(QV(K)-max(PROD,-QC(K))+ERN,0.)
       QC(K) = QC(K)+max(PROD,-QC(K))
       QR(K) = QR(K)-ERN
   END DO

 END SUBROUTINE KESSLER1

 subroutine reed_simple_physics (pcols, pver, dtime, t, q, u, v, pmid, pint, pdel, rpdel, ps, precl, cond_only, dqdt, dTdt, dudt, dvdt)
   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: Simple Physics Package
   !
   ! Author: K. A. Reed (University of Michigan, kareed@umich.edu)
   !         version 4 
   !         May/31/2012
   !         uploaded on June/8/2012
   !
   !         cosmetic changes in comparison to version 3:
   !         parameter list now includes pcols and pver
   !         ncol has been replaced by pcols in the code
   ! 
   ! Description: Includes large-scale precipitation, surface fluxes and
   !              boundary-leyer mixing. The processes are time-split
   !              in that order. A partially implicit formulation is
   !              used to foster numerical stability.
   !              The routine assumes that the model levels are ordered
   !              in a top-down approach, e.g. level 1 denotes the uppermost
   !              full model level.
   !
   !              This routine is based on an implementation which was
   !              developed for the NCAR Community Atmosphere Model (CAM).
   !              Adjustments for other models will be necessary.
   !
   !              The routine provides both updates of the state variables
   !              u, v, T, q (these are local copies of u,v,T,q within this physics
   !              routine) and also collects their time tendencies.
   !              The latter might be used to couple the physics and dynamics
   !              in a process-split way. For a time-split coupling, the final
   !              state should be given to the dynamical core for the next time step.
   !
   !
   ! Reference: Reed, K. A. and C. Jablonowski (2012), Idealized tropical cyclone 
   !            simulations of intermediate complexity: A test case for AGCMs, 
   !            J. Adv. Model. Earth Syst., Vol. 4, M04001, doi:10.1029/2011MS000099
   ! 
   ! NOTE: 17 jul 12 lmh
   !       This has been modified to return only the tendencies and not
   !       to update the model; the model variables are updated in
   !       fv_update_phys.F90. However we do need to use the
   !       tendency computed by the surface fluxes to advance the
   !       boundary layer mixing to be in line with the original
   !       formulation of these physics routines.
   !--------------------------------------------------------------------

   ! use physics_types     , only: physics_dme_adjust   ! This is for CESM/CAM
   ! use cam_diagnostics,    only: diag_phys_writeout   ! This is for CESM/CAM

   implicit none

   integer, parameter :: r8 = selected_real_kind(12)

   !
   ! Input arguments - MODEL DEPENDENT
   !
   integer, intent(in)  :: pcols        ! Set number of atmospheric columns       
   integer, intent(in)  :: pver         ! Set number of model levels
   real(r8), intent(in) :: dtime        ! Set model physics timestep
   logical, intent(IN) :: cond_only     ! Only do large-scale condensation
   !
   ! Input/Output arguments 
   !
   !  pcols is the maximum number of vertical columns per 'chunk' of atmosphere
   !
   real(r8), intent(inout) :: t(pcols,pver)      ! Temperature at full-model level (K)
   real(r8), intent(inout) :: q(pcols,pver)      ! Specific Humidity at full-model level (kg/kg)
   real(r8), intent(inout) :: u(pcols,pver)      ! Zonal wind at full-model level (m/s)
   real(r8), intent(inout) :: v(pcols,pver)      ! Meridional wind at full-model level (m/s)
   real(r8), intent(inout) :: pmid(pcols,pver)   ! Pressure is full-model level (Pa)
   real(r8), intent(inout) :: pint(pcols,pver+1) ! Pressure at model interfaces (Pa)
   real(r8), intent(inout) :: pdel(pcols,pver)   ! Layer thickness (Pa)
   real(r8), intent(inout) :: rpdel(pcols,pver)  ! Reciprocal of layer thickness (1/Pa)
   real(r8), intent(inout) :: ps(pcols)          ! Surface Pressue (Pa)
   !
   ! Output arguments 
   !
   real(r8), intent(out) :: precl(pcols)         ! Precipitation rate (m_water / s)
   !   real(r8), intent(OUT)   :: RH(pcols,pver)     !Output RH

   !
   !---------------------------Local workspace-----------------------------
   !

   ! Integers for loops

   integer  i,k                         ! Longitude, level indices

   ! Physical Constants - Many of these may be model dependent

   real(r8) gravit                      ! Gravity
   real(r8) rair                        ! Gas constant for dry air
   real(r8) cpair                       ! Specific heat of dry air 
   real(r8) latvap                      ! Latent heat of vaporization
   real(r8) rh2o                        ! Gas constant for water vapor
   real(r8) epsilo                      ! Ratio of gas constant for dry air to that for vapor
   real(r8) zvir                        ! Constant for virtual temp. calc. =(rh2o/rair) - 1

   ! Simple Physics Specific Constants 

   real(r8) T0                          ! Control temp for calculation of qsat
   real(r8) e0                          ! Saturation vapor pressure at T0 for calculation of qsat
   real(r8) rhow                        ! Density of Liquid Water
   real(r8) p0                          ! Constant for calculation of potential temperature
   real(r8) Cd0                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cd1                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cm                          ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) v20                         ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
   real(r8) C                           ! Drag coefficient for sensible heat and evaporation

   ! Physics Tendency Arrays
   real(r8), intent(OUT) ::  dtdt(pcols,pver)             ! Temperature tendency 
   real(r8), intent(OUT) ::  dqdt(pcols,pver)             ! Specific humidity tendency
   real(r8), intent(OUT) ::  dudt(pcols,pver)             ! Zonal wind tendency
   real(r8), intent(OUT) ::  dvdt(pcols,pver)             ! Meridional wind tendency

!!$   !Temporary variables for continual updates
!!$   real(r8) :: tt(pcols,pver)      ! Temperature at full-model level (K)
!!$   real(r8) :: qq(pcols,pver)      ! Specific Humidity at full-model level (kg/kg)
!!$   real(r8) :: uu(pcols,pver)      ! Zonal wind at full-model level (m/s)
!!$   real(r8) :: vv(pcols,pver)      ! Meridional wind at full-model level (m/s)

   real(r8) Tsurf(pcols) 

   ! Temporary variables for tendency calculations

   real(r8) tmp                         ! Temporary
   real(r8) qsat                        ! Saturation vapor pressure
   real(r8) qsats                       ! Saturation vapor pressure of SST

   ! Variables for Boundary Layer Calculation

   real(r8) wind(pcols)                 ! Magnitude of Wind
   real(r8) Cd(pcols)                   ! Drag coefficient for momentum
   real(r8) Km(pcols,pver+1)            ! Eddy diffusivity for boundary layer calculations 
   real(r8) Ke(pcols,pver+1)            ! Eddy diffusivity for boundary layer calculations
   real(r8) rho                         ! Density at lower/upper interface
   real(r8) za(pcols)                   ! Heights at midpoints of first model level
   real(r8) dlnpint                     ! Used for calculation of heights
   real(r8) pbltop                      ! Top of boundary layer
   real(r8) pblconst                    ! Constant for the calculation of the decay of diffusivity 
   real(r8) CA(pcols,pver)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CC(pcols,pver)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CE(pcols,pver+1)            ! Matrix Coefficents for PBL Scheme
   real(r8) CAm(pcols,pver)             ! Matrix Coefficents for PBL Scheme
   real(r8) CCm(pcols,pver)             ! Matrix Coefficents for PBL Scheme
   real(r8) CEm(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFu(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFv(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFt(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFq(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme

   !===============================================================================
   !
   ! Physical Constants - MAY BE MODEL DEPENDENT 
   !
   !===============================================================================
   gravit = 9.80616_r8                   ! Gravity
   rair   = 287.0_r8                     ! Gas constant for dry air
   cpair  = 1.004e3_r8                   ! Specific heat of dry air 
   latvap = 2.5e6_r8                     ! Latent heat of vaporization
   rh2o   = 461.5_r8                     ! Gas constant for water vapor
   epsilo = rair/rh2o                    ! Ratio of gas constant for dry air to that for vapor
   zvir   = (rh2o/rair) - 1._r8          ! Constant for virtual temp. calc. =(rh2o/rair) - 1 is approx. 0.608

   !===============================================================================
   !
   ! Local Constants for Simple Physics
   !
   !===============================================================================
   C        = 0.0011_r8      ! From Smith and Vogl 2008
   Tsurf    = 302.15_r8      ! Constant Value for SST
   T0       = 273.16_r8      ! control temp for calculation of qsat
   e0       = 610.78_r8      ! saturation vapor pressure at T0 for calculation of qsat
   rhow     = 1000.0_r8      ! Density of Liquid Water 
   Cd0      = 0.0007_r8      ! Constant for Cd calc. Smith and Vogl 2008
   Cd1      = 0.000065_r8    ! Constant for Cd calc. Smith and Vogl 2008
   Cm       = 0.002_r8       ! Constant for Cd calc. Smith and Vogl 2008
   v20      = 20.0_r8        ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
   p0       = 100000.0_r8    ! Constant for potential temp calculation
   pbltop   = 85000._r8      ! Top of boundary layer
   pblconst = 10000._r8      ! Constant for the calculation of the decay of diffusivity

   !===============================================================================
   !
   ! Definition of local arrays
   !
   !===============================================================================


   !
   ! Calculate hydrostatic height za of the lowest model level
   !
   do i=1,pcols 
      dlnpint = log(ps(i)) - log(pint(i,pver))                 ! ps(i) is identical to pint(i,pver+1), note: this is the correct sign (corrects typo in paper)
      za(i) = rair/gravit*t(i,pver)*(1._r8+zvir*q(i,pver))*0.5_r8*dlnpint
   end do

   !===============================================================================
   !
   ! Set initial physics time tendencies and precipitation field to zero
   !
   !===============================================================================
   do k=1,pver
      do i=1,pcols
         dtdt(i,k)  = 0._r8            ! initialize temperature tendency with zero
         dqdt(i,k)  = 0._r8            ! initialize specific humidity tendency with zero
         dudt(i,k)  = 0._r8            ! initialize zonal wind tendency with zero
         dvdt(i,k)  = 0._r8            ! initialize meridional wind tendency with zero
      enddo
   enddo
   do i=1,pcols
      precl(i) = 0._r8                  ! initialize precipitation rate with zero
   enddo
   !===============================================================================
   !
   ! Large-Scale Condensation and Precipitation Rate
   !
   !===============================================================================
   !
   ! Calculate Tendencies
   !
   do k=1,pver
      do i=1,pcols
         qsat = epsilo*e0/pmid(i,k)*exp(-latvap/rh2o*((1._r8/t(i,k))-1._r8/T0))
         if (q(i,k) > qsat) then                                                 ! saturated?
            tmp  = 1._r8/dtime*(q(i,k)-qsat)/(1._r8+(latvap/cpair)*(epsilo*latvap*qsat/(rair*t(i,k)**2)))
            dtdt(i,k) = dtdt(i,k)+latvap/cpair*tmp
            dqdt(i,k) = dqdt(i,k)-tmp
            precl(i) = precl(i) + tmp*pdel(i,k)/(gravit*rhow)*1000.*86400.                    ! precipitation rate, computed via a vertical integral, converted to m/d
            ! corrected in version 1.3
         end if
         !rh(i,k) = q(i,k)/qsat*100.
      end do
   end do

   !===============================================================================
   !
   ! Turbulent mixing coefficients for the PBL mixing of horizontal momentum,
   ! sensible heat and latent heat
   !
   ! We are using Simplified Ekman theory to compute the diffusion coefficients 
   ! Kx for the boundary-layer mixing. The Kx values are calculated at each time step
   ! and in each column.
   !
   !===============================================================================

   !if (.not. cond_only) then
      !
      ! Compute magnitude of the wind and drag coeffcients for turbulence scheme:
      ! they depend on the conditions at the lowest model level and stay constant
      ! up to the 850 hPa level. Above this level the coefficients are decreased
      ! and tapered to zero. At the 700 hPa level the strength of the K coefficients
      ! is about 10% of the maximum strength. 
      !
      do i=1,pcols
         wind(i) = sqrt(u(i,pver)**2+v(i,pver)**2)    ! wind magnitude at the lowest level
      end do
      do i=1,pcols
         Ke(i,pver+1) = C*wind(i)*za(i)
         if( wind(i) .lt. v20) then
            Cd(i) = Cd0+Cd1*wind(i) 
            Km(i,pver+1) = Cd(i)*wind(i)*za(i)
         else
            Cd(i) = Cm
            Km(i,pver+1) = Cm*wind(i)*za(i)
         endif
      end do

      do k=1,pver
         do i=1,pcols
            if( pint(i,k) .ge. pbltop) then
               Km(i,k) = Km(i,pver+1)                 ! constant Km below 850 hPa level
               Ke(i,k) = Ke(i,pver+1)                 ! constant Ke below 850 hPa level
            else
               Km(i,k) = Km(i,pver+1)*exp(-(pbltop-pint(i,k))**2/(pblconst)**2)  ! Km tapered to 0
               Ke(i,k) = Ke(i,pver+1)*exp(-(pbltop-pint(i,k))**2/(pblconst)**2)  ! Ke tapered to 0
            end if
         end do
      end do


      !===============================================================================
      ! Update the state variables u, v, t, q with the surface fluxes at the
      ! lowest model level, this is done with an implicit approach
      ! see Reed and Jablonowski (JAMES, 2012)
      !
      ! Sea Surface Temperature Tsurf is constant for tropical cyclone test 5-1
      ! Tsurf needs to be dependent on latitude for the
      ! moist baroclinic wave test 4-3 with simple-physics, adjust
      !===============================================================================

      !!! DEBUG CODE
      !write(1000,*) Tsurf

      do i=1,pcols 
         !FIXME: Need to continually update during this routine instead of applying tendencies later?
         qsats = epsilo*e0/ps(i)*exp(-latvap/rh2o*((1._r8/Tsurf(i))-1._r8/T0))
         dudt(i,pver) = dudt(i,pver) + (u(i,pver) &
              /(1._r8+Cd(i)*wind(i)*dtime/za(i))-u(i,pver))/dtime
         dvdt(i,pver) = dvdt(i,pver) + (v(i,pver) &
              /(1._r8+Cd(i)*wind(i)*dtime/za(i))-v(i,pver))/dtime
         !u(i,pver)   = u(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))
         !v(i,pver)   = v(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))
         dtdt(i,pver) = dtdt(i,pver) +((t(i,pver)+C*wind(i)*Tsurf(i)*dtime/za(i)) &
              /(1._r8+C*wind(i)*dtime/za(i))-t(i,pver))/dtime 
         !t(i,pver)   = (t(i,pver)+C*wind(i)*Tsurf(i)*dtime/za(i)) &
         !     /(1._r8+C*wind(i)*dtime/za(i))  
         dqdt(i,pver) = dqdt(i,pver) +((q(i,pver)+C*wind(i)*qsats*dtime/za(i)) &
              /(1._r8+C*wind(i)*dtime/za(i))-q(i,pver))/dtime
         !q(i,pver) = (q(i,pver)+C*wind(i)*qsats*dtime/za(i))/(1._r8+C*wind(i)*dtime/za(i))
      end do
      !===============================================================================


      !===============================================================================
      ! Boundary layer mixing, see Reed and Jablonowski (JAMES, 2012)
      !===============================================================================
      ! Calculate Diagonal Variables for Implicit PBL Scheme
      !
      do k=1,pver-1
         do i=1,pcols
            rho = (pint(i,k+1)/(rair*(t(i,k+1)+t(i,k))/2.0_r8))
            CAm(i,k)   = rpdel(i,k)*dtime*gravit*gravit*Km(i,k+1)*rho*rho   &
                 /(pmid(i,k+1)-pmid(i,k))    
            CCm(i,k+1) = rpdel(i,k+1)*dtime*gravit*gravit*Km(i,k+1)*rho*rho &
                 /(pmid(i,k+1)-pmid(i,k))
            CA(i,k)    = rpdel(i,k)*dtime*gravit*gravit*Ke(i,k+1)*rho*rho   &
                 /(pmid(i,k+1)-pmid(i,k))
            CC(i,k+1)  = rpdel(i,k+1)*dtime*gravit*gravit*Ke(i,k+1)*rho*rho &
                 /(pmid(i,k+1)-pmid(i,k))
         end do
      end do
      do i=1,pcols
         CAm(i,pver) = 0._r8
         CCm(i,1) = 0._r8
         CEm(i,pver+1) = 0._r8
         CA(i,pver) = 0._r8
         CC(i,1) = 0._r8
         CE(i,pver+1) = 0._r8
         CFu(i,pver+1) = 0._r8
         CFv(i,pver+1) = 0._r8
         CFt(i,pver+1) = 0._r8
         CFq(i,pver+1) = 0._r8 
      end do
      do i=1,pcols
         do k=pver,1,-1
            CE(i,k)  = CC(i,k)/(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1)) 
            CEm(i,k) = CCm(i,k)/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFu(i,k) = (u(i,k)+CAm(i,k)*CFu(i,k+1)) &
                 /(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFv(i,k) = (v(i,k)+CAm(i,k)*CFv(i,k+1)) &
                 /(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFt(i,k) = ((p0/pmid(i,k))**(rair/cpair)*t(i,k)+CA(i,k)*CFt(i,k+1)) &
                 /(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1)) 
            CFq(i,k) = (q(i,k)+CA(i,k)*CFq(i,k+1)) &
                 /(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1))
         end do
      end do

      !
      ! Calculate the updated temperature, specific humidity and horizontal wind
      !
      ! First we need to calculate the updates at the top model level
      !
      do i=1,pcols
         dudt(i,1)  = dudt(i,1)+(CFu(i,1)-u(i,1))/dtime
         dvdt(i,1)  = dvdt(i,1)+(CFv(i,1)-v(i,1))/dtime
         !u(i,1)    = CFu(i,1)
         !v(i,1)    = CFv(i,1)
         dtdt(i,1)  = dtdt(i,1)+(CFt(i,1)*(pmid(i,1)/p0)**(rair/cpair)-t(i,1))/dtime  ! corrected in version 1.3
         !t(i,1)    = CFt(i,1)*(pmid(i,1)/p0)**(rair/cpair)
         dqdt(i,1)  = dqdt(i,1)+(CFq(i,1)-q(i,1))/dtime
         !q(i,1)  = CFq(i,1)
      end do
      !
      ! Loop over the remaining level
      !
      do i=1,pcols
         do k=2,pver
            dudt(i,k)  = dudt(i,k)+(CEm(i,k)*u(i,k-1)+CFu(i,k)-u(i,k))/dtime
            dvdt(i,k)  = dvdt(i,k)+(CEm(i,k)*v(i,k-1)+CFv(i,k)-v(i,k))/dtime
            !uu(i,k)    = CEm(i,k)*u(i,k-1)+CFu(i,k) 
            !vv(i,k)    = CEm(i,k)*v(i,k-1)+CFv(i,k)
            dtdt(i,k)  = dtdt(i,k)+((CE(i,k)*t(i,k-1) &
                 *(p0/pmid(i,k-1))**(rair/cpair)+CFt(i,k)) &
                 *(pmid(i,k)/p0)**(rair/cpair)-t(i,k))/dtime 
            !tt(i,k)    = (CE(i,k)*t(i,k-1)*(p0/pmid(i,k-1))**(rair/cpair)+CFt(i,k)) &
            !     *(pmid(i,k)/p0)**(rair/cpair)
            dqdt(i,k)  = dqdt(i,k)+(CE(i,k)*q(i,k-1)+CFq(i,k)-q(i,k))/dtime
            !qq(i,k)  = CE(i,k)*q(i,k-1)+CFq(i,k)
         end do
      end do

   !endif !.not. cond_only

   return
 end subroutine reed_simple_physics


subroutine reed_simple_physics_v5 (pcols, pver, dtime, lat, t, q, u, v, pmid, pint, pdel, rpdel, ps, precl, test)
!----------------------------------------------------------------------- 
! 
! Purpose: Simple Physics Package
!
! Author: K. A. Reed (University of Michigan, kareed@umich.edu)
!         version 5 
!         July/8/2012
!
!  Change log:
!  v2: removal of some NCAR CAM-specific 'use' associations
!  v3: corrected precl(i) computation, the precipitation rate is now computed via a vertical integral, the previous single-level computation in v2 was a bug
!  v3: corrected dtdt(i,1) computation, the term '-(i,1)' was missing the temperature variable: '-t(i,1)'
!  v4: modified and enhanced parameter list to make the routine truly standalone, the number of columns and vertical levels have been added: pcols, pver
!  v4: 'ncol' has been removed, 'pcols' is used instead
!  v5: the sea surface temperature (SST) field Tsurf is now an array, the SST now depends on the latitude
!  v5: addition of the latitude array 'lat' and the flag 'test' in the parameter list
!      if test = 0: constant SST is used, correct setting for the tropical cyclone test case 5-1
!      if test = 1: newly added latitude-dependent SST is used, correct setting for the moist baroclinic wave test with simple-physics (test 4-3)
! 
! Description: Includes large-scale precipitation, surface fluxes and
!              boundary-leyer mixing. The processes are time-split
!              in that order. A partially implicit formulation is
!              used to foster numerical stability.
!              The routine assumes that the model levels are ordered
!              in a top-down approach, e.g. level 1 denotes the uppermost
!              full model level.
!
!              This routine is based on an implementation which was
!              developed for the NCAR Community Atmosphere Model (CAM).
!              Adjustments for other models will be necessary.
!
!              The routine provides both updates of the state variables
!              u, v, T, q (these are local copies of u,v,T,q within this physics
!              routine) and also collects their time tendencies.
!              The latter might be used to couple the physics and dynamics
!              in a process-split way. For a time-split coupling, the final
!              state should be given to the dynamical core for the next time step.
! Test:      0 = Reed and Jablonowski (2011) tropical cyclone test case (test 5-1)
!            1 = Moist baroclinic instability test (test 4-3)
!
!
! Reference: Reed, K. A. and C. Jablonowski (2012), Idealized tropical cyclone 
!            simulations of intermediate complexity: A test case for AGCMs, 
!            J. Adv. Model. Earth Syst., Vol. 4, M04001, doi:10.1029/2011MS000099
!-----------------------------------------------------------------------
  ! use physics_types     , only: physics_dme_adjust   ! This is for CESM/CAM
  ! use cam_diagnostics,    only: diag_phys_writeout   ! This is for CESM/CAM

   implicit none

   integer, parameter :: r8 = selected_real_kind(12)

!
! Input arguments - MODEL DEPENDENT
!
   integer, intent(in)  :: pcols        ! Set number of atmospheric columns       
   integer, intent(in)  :: pver         ! Set number of model levels
   real(r8), intent(in) :: dtime        ! Set model physics timestep
   real(r8), intent(in) :: lat(pcols)   ! Latitude 
   integer, intent(in)  :: test         ! Test number
   
!
! Input/Output arguments 
!
!  pcols is the maximum number of vertical columns per 'chunk' of atmosphere
!
   real(r8), intent(inout) :: t(pcols,pver)      ! Temperature at full-model level (K)
   real(r8), intent(inout) :: q(pcols,pver)      ! Specific Humidity at full-model level (kg/kg)
   real(r8), intent(inout) :: u(pcols,pver)      ! Zonal wind at full-model level (m/s)
   real(r8), intent(inout) :: v(pcols,pver)      ! Meridional wind at full-model level (m/s)
   real(r8), intent(inout) :: pmid(pcols,pver)   ! Pressure is full-model level (Pa)
   real(r8), intent(inout) :: pint(pcols,pver+1) ! Pressure at model interfaces (Pa)
   real(r8), intent(inout) :: pdel(pcols,pver)   ! Layer thickness (Pa)
   real(r8), intent(inout) :: rpdel(pcols,pver)  ! Reciprocal of layer thickness (1/Pa)
   real(r8), intent(inout) :: ps(pcols)          ! Surface Pressue (Pa)

!
! Output arguments 
!
   real(r8), intent(out) :: precl(pcols)         ! Precipitation rate (m_water / s)

!
!---------------------------Local workspace-----------------------------
!

! Integers for loops

   integer  i,k                         ! Longitude, level indices

! Physical Constants - Many of these may be model dependent

   real(r8) gravit                      ! Gravity
   real(r8) rair                        ! Gas constant for dry air
   real(r8) cpair                       ! Specific heat of dry air 
   real(r8) latvap                      ! Latent heat of vaporization
   real(r8) rh2o                        ! Gas constant for water vapor
   real(r8) epsilo                      ! Ratio of gas constant for dry air to that for vapor
   real(r8) zvir                        ! Constant for virtual temp. calc. =(rh2o/rair) - 1
   real(r8) a                           ! Reference Earth's Radius (m)
   real(r8) omega                       ! Reference rotation rate of the Earth (s^-1)
   real(r8) pi                          ! pi

! Simple Physics Specific Constants 

!++++++++                     
   real(r8) Tsurf(pcols)                ! Sea Surface Temperature (constant for tropical cyclone)
!++++++++                                 Tsurf needs to be dependent on latitude for the
                                        ! moist baroclinic wave test 4-3 with simple-physics, adjust

   real(r8) SST_tc                      ! Sea Surface Temperature for tropical cyclone test
   real(r8) T0                          ! Control temp for calculation of qsat
   real(r8) e0                          ! Saturation vapor pressure at T0 for calculation of qsat
   real(r8) rhow                        ! Density of Liquid Water
   real(r8) p0                          ! Constant for calculation of potential temperature
   real(r8) Cd0                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cd1                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cm                          ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) v20                         ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
   real(r8) C                           ! Drag coefficient for sensible heat and evaporation
   real(r8) T00                         ! Horizontal mean T at surface for moist baro test
   real(r8) u0                          ! Zonal wind constant for moist baro test
   real(r8) latw                        ! halfwidth for  for baro test
   real(r8) eta0                        ! Center of jets (hybrid) for baro test
   real(r8) etav                        ! Auxiliary variable for baro test
   real(r8) q0                          ! Maximum specific humidity for baro test

! Physics Tendency Arrays
  real(r8) dtdt(pcols,pver)             ! Temperature tendency 
  real(r8) dqdt(pcols,pver)             ! Specific humidity tendency
  real(r8) dudt(pcols,pver)             ! Zonal wind tendency
  real(r8) dvdt(pcols,pver)             ! Meridional wind tendency

! Temporary variables for tendency calculations

   real(r8) tmp                         ! Temporary
   real(r8) qsat                        ! Saturation vapor pressure
   real(r8) qsats                       ! Saturation vapor pressure of SST

! Variables for Boundary Layer Calculation

   real(r8) wind(pcols)                 ! Magnitude of Wind
   real(r8) Cd(pcols)                   ! Drag coefficient for momentum
   real(r8) Km(pcols,pver+1)            ! Eddy diffusivity for boundary layer calculations 
   real(r8) Ke(pcols,pver+1)            ! Eddy diffusivity for boundary layer calculations
   real(r8) rho                         ! Density at lower/upper interface
   real(r8) za(pcols)                   ! Heights at midpoints of first model level
   real(r8) dlnpint                     ! Used for calculation of heights
   real(r8) pbltop                      ! Top of boundary layer
   real(r8) pblconst                    ! Constant for the calculation of the decay of diffusivity 
   real(r8) CA(pcols,pver)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CC(pcols,pver)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CE(pcols,pver+1)            ! Matrix Coefficents for PBL Scheme
   real(r8) CAm(pcols,pver)             ! Matrix Coefficents for PBL Scheme
   real(r8) CCm(pcols,pver)             ! Matrix Coefficents for PBL Scheme
   real(r8) CEm(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFu(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFv(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFt(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFq(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme


! Variable for Dry Mass Adjustment, this dry air adjustment is necessary to
! conserve the mass of the dry air

   real(r8) qini(pcols,pver)            ! Initial specific humidity

!===============================================================================
!
! Physical Constants - MAY BE MODEL DEPENDENT 
!
!===============================================================================
   gravit = 9.80616_r8                   ! Gravity (9.80616 m/s^2)
   rair   = 287.0_r8                     ! Gas constant for dry air: 287 J/(kg K)
   cpair  = 1.0045e3_r8                  ! Specific heat of dry air: here we use 1004.5 J/(kg K)
   latvap = 2.5e6_r8                     ! Latent heat of vaporization (J/kg)
   rh2o   = 461.5_r8                     ! Gas constant for water vapor: 461.5 J/(kg K)
   epsilo = rair/rh2o                    ! Ratio of gas constant for dry air to that for vapor
   zvir   = (rh2o/rair) - 1._r8          ! Constant for virtual temp. calc. =(rh2o/rair) - 1 is approx. 0.608
   a      = 6371220.0_r8                 ! Reference Earth's Radius (m)
   omega  = 7.29212d-5                   ! Reference rotation rate of the Earth (s^-1)
   pi     = 4._r8*atan(1._r8)            ! pi

!===============================================================================
!
! Local Constants for Simple Physics
!
!===============================================================================
      C        = 0.0011_r8      ! From Smith and Vogl 2008
      SST_tc   = 302.15_r8      ! Constant Value for SST for tropical cyclone test
      T0       = 273.16_r8      ! control temp for calculation of qsat
      e0       = 610.78_r8      ! saturation vapor pressure at T0 for calculation of qsat
      rhow     = 1000.0_r8      ! Density of Liquid Water 
      Cd0      = 0.0007_r8      ! Constant for Cd calc. Smith and Vogl 2008
      Cd1      = 0.000065_r8    ! Constant for Cd calc. Smith and Vogl 2008
      Cm       = 0.002_r8       ! Constant for Cd calc. Smith and Vogl 2008
      v20      = 20.0_r8        ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
      p0       = 100000.0_r8    ! Constant for potential temp calculation
      pbltop   = 85000._r8      ! Top of boundary layer
      pblconst = 10000._r8      ! Constant for the calculation of the decay of diffusivity
      T00      = 288.0_r8         ! Horizontal mean T at surface for moist baro test
      u0       = 35.0_r8          ! Zonal wind constant for moist baro test
      latw     = 2.0_r8*pi/9.0_r8 ! Halfwidth for  for baro test
      eta0     = 0.252_r8         ! Center of jets (hybrid) for baro test
      etav     = (1._r8-eta0)*0.5_r8*pi ! Auxiliary variable for baro test
      q0       = 0.021            ! Maximum specific humidity for baro test

!===============================================================================
!
! Definition of local arrays
!
!===============================================================================
!
! Calculate hydrostatic height za of the lowest model level
!
     do i=1,pcols 
        dlnpint = log(ps(i)) - log(pint(i,pver))                 ! ps(i) is identical to pint(i,pver+1), note: this is the correct sign (corrects typo in JAMES paper)
        za(i) = rair/gravit*t(i,pver)*(1._r8+zvir*q(i,pver))*0.5_r8*dlnpint
     end do
!
! Set Initial Specific Humidity - For dry mass adjustment at the end
!
     qini(:pcols,:pver) = q(:pcols,:pver)
!--------------------------------------------------------------
! Set Sea Surface Temperature (constant for tropical cyclone)
! Tsurf needs to be dependent on latitude for the
! moist baroclinic wave test 4-3 with simple-physics 
!--------------------------------------------------------------
     if (test .eq. 1) then     ! moist baroclinic wave with simple-physics
        do i=1,pcols
           Tsurf(i) = (T00 + pi*u0/rair * 1.5_r8 * sin(etav) * (cos(etav))**0.5_r8 *                 &
                     ((-2._r8*(sin(lat(i)))**6 * ((cos(lat(i)))**2 + 1._r8/3._r8) + 10._r8/63._r8)* &
                     u0 * (cos(etav))**1.5_r8  +                                                    &
                     (8._r8/5._r8*(cos(lat(i)))**3 * ((sin(lat(i)))**2 + 2._r8/3._r8) - pi/4._r8)*a*omega*0.5_r8 ))/ &
                     (1._r8+zvir*q0*exp(-(lat(i)/latw)**4))

        end do
     else
        do i=1,pcols          ! constant SST for the tropical cyclone test case
           Tsurf(i) = SST_tc
        end do
     end if

!===============================================================================
!
! Set initial physics time tendencies and precipitation field to zero
!
!===============================================================================
     dtdt(:pcols,:pver)  = 0._r8            ! initialize temperature tendency with zero
     dqdt(:pcols,:pver)  = 0._r8            ! initialize specific humidity tendency with zero
     dudt(:pcols,:pver)  = 0._r8            ! initialize zonal wind tendency with zero
     dvdt(:pcols,:pver)  = 0._r8            ! initialize meridional wind tendency with zero
     precl(:pcols) = 0._r8                  ! initialize precipitation rate with zero

!===============================================================================
!
! Large-Scale Condensation and Precipitation Rate
!
!===============================================================================
!
! Calculate Tendencies
!
      do k=1,pver
         do i=1,pcols
            qsat = epsilo*e0/pmid(i,k)*exp(-latvap/rh2o*((1._r8/t(i,k))-1._r8/T0))  ! saturation specific humidity
            if (q(i,k) > qsat) then                                                 ! saturated?
               tmp  = 1._r8/dtime*(q(i,k)-qsat)/(1._r8+(latvap/cpair)*(epsilo*latvap*qsat/(rair*t(i,k)**2)))
               dtdt(i,k) = dtdt(i,k)+latvap/cpair*tmp
               dqdt(i,k) = dqdt(i,k)-tmp
               precl(i) = precl(i) + tmp*pdel(i,k)/(gravit*rhow)                    ! precipitation rate, computed via a vertical integral
                                                                                    ! corrected in version 1.3
            end if
         end do
      end do
!
! Update moisture and temperature fields from Large-Scale Precipitation Scheme
!
      do k=1,pver
         do i=1,pcols
            t(i,k) =  t(i,k) + dtdt(i,k)*dtime    ! update the state variables T and q
            q(i,k) =  q(i,k) + dqdt(i,k)*dtime
         end do
      end do

!===============================================================================
! Send variables to history file - THIS PROCESS WILL BE MODEL SPECIFIC
!
! note: The variables, as done in many GCMs, are written to the history file
!       after the moist physics process.  This ensures that the moisture fields
!       are somewhat in equilibrium.
!===============================================================================
  !  call diag_phys_writeout(state)   ! This is for CESM/CAM

!===============================================================================
!
! Turbulent mixing coefficients for the PBL mixing of horizontal momentum,
! sensible heat and latent heat
!
! We are using Simplified Ekman theory to compute the diffusion coefficients 
! Kx for the boundary-layer mixing. The Kx values are calculated at each time step
! and in each column.
!
!===============================================================================
!
! Compute magnitude of the wind and drag coeffcients for turbulence scheme:
! they depend on the conditions at the lowest model level and stay constant
! up to the 850 hPa level. Above this level the coefficients are decreased
! and tapered to zero. At the 700 hPa level the strength of the K coefficients
! is about 10% of the maximum strength. 
!
     do i=1,pcols
        wind(i) = sqrt(u(i,pver)**2+v(i,pver)**2)    ! wind magnitude at the lowest level
     end do
     do i=1,pcols
        Ke(i,pver+1) = C*wind(i)*za(i)
        if( wind(i) .lt. v20) then
           Cd(i) = Cd0+Cd1*wind(i) 
           Km(i,pver+1) = Cd(i)*wind(i)*za(i)
        else
           Cd(i) = Cm
           Km(i,pver+1) = Cm*wind(i)*za(i)
        endif
     end do

      do k=1,pver
         do i=1,pcols
            if( pint(i,k) .ge. pbltop) then
               Km(i,k) = Km(i,pver+1)                 ! constant Km below 850 hPa level
               Ke(i,k) = Ke(i,pver+1)                 ! constant Ke below 850 hPa level
            else
               Km(i,k) = Km(i,pver+1)*exp(-(pbltop-pint(i,k))**2/(pblconst)**2)  ! Km tapered to 0
               Ke(i,k) = Ke(i,pver+1)*exp(-(pbltop-pint(i,k))**2/(pblconst)**2)  ! Ke tapered to 0
            end if 
         end do
      end do     


!===============================================================================
! Update the state variables u, v, t, q with the surface fluxes at the
! lowest model level, this is done with an implicit approach
! see Reed and Jablonowski (JAMES, 2012)
!
! Sea Surface Temperature Tsurf is constant for tropical cyclone test 5-1
! Tsurf needs to be dependent on latitude for the
! moist baroclinic wave test 4-3 with simple-physics, adjust
!===============================================================================

     do i=1,pcols
        qsats = epsilo*e0/ps(i)*exp(-latvap/rh2o*((1._r8/Tsurf(i))-1._r8/T0))  ! saturation specific humidity at the surface
        dudt(i,pver) = dudt(i,pver) + (u(i,pver) &
                            /(1._r8+Cd(i)*wind(i)*dtime/za(i))-u(i,pver))/dtime
        dvdt(i,pver) = dvdt(i,pver) + (v(i,pver) &
                            /(1._r8+Cd(i)*wind(i)*dtime/za(i))-v(i,pver))/dtime
        u(i,pver)   = u(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))
        v(i,pver)   = v(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))
        dtdt(i,pver) = dtdt(i,pver) +((t(i,pver)+C*wind(i)*Tsurf(i)*dtime/za(i)) &
                            /(1._r8+C*wind(i)*dtime/za(i))-t(i,pver))/dtime 
        t(i,pver)   = (t(i,pver)+C*wind(i)*Tsurf(i)*dtime/za(i)) &
                            /(1._r8+C*wind(i)*dtime/za(i))  
        dqdt(i,pver) = dqdt(i,pver) +((q(i,pver)+C*wind(i)*qsats*dtime/za(i)) &
                            /(1._r8+C*wind(i)*dtime/za(i))-q(i,pver))/dtime
        q(i,pver) = (q(i,pver)+C*wind(i)*qsats*dtime/za(i))/(1._r8+C*wind(i)*dtime/za(i))
     end do
!===============================================================================


!===============================================================================
! Boundary layer mixing, see Reed and Jablonowski (JAMES, 2012)
!===============================================================================
! Calculate Diagonal Variables for Implicit PBL Scheme
!
      do k=1,pver-1
         do i=1,pcols
            rho = (pint(i,k+1)/(rair*(t(i,k+1)+t(i,k))/2.0_r8))
            CAm(i,k)   = rpdel(i,k)*dtime*gravit*gravit*Km(i,k+1)*rho*rho   &
                         /(pmid(i,k+1)-pmid(i,k))    
            CCm(i,k+1) = rpdel(i,k+1)*dtime*gravit*gravit*Km(i,k+1)*rho*rho &
                         /(pmid(i,k+1)-pmid(i,k))
            CA(i,k)    = rpdel(i,k)*dtime*gravit*gravit*Ke(i,k+1)*rho*rho   &
                         /(pmid(i,k+1)-pmid(i,k))
            CC(i,k+1)  = rpdel(i,k+1)*dtime*gravit*gravit*Ke(i,k+1)*rho*rho &
                         /(pmid(i,k+1)-pmid(i,k))
         end do
      end do
      do i=1,pcols
         CAm(i,pver) = 0._r8
         CCm(i,1) = 0._r8
         CEm(i,pver+1) = 0._r8
         CA(i,pver) = 0._r8
         CC(i,1) = 0._r8
         CE(i,pver+1) = 0._r8
         CFu(i,pver+1) = 0._r8
         CFv(i,pver+1) = 0._r8
         CFt(i,pver+1) = 0._r8
         CFq(i,pver+1) = 0._r8 
      end do
      do i=1,pcols
         do k=pver,1,-1
            CE(i,k)  = CC(i,k)/(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1)) 
            CEm(i,k) = CCm(i,k)/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFu(i,k) = (u(i,k)+CAm(i,k)*CFu(i,k+1)) &
                       /(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFv(i,k) = (v(i,k)+CAm(i,k)*CFv(i,k+1)) &
                       /(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFt(i,k) = ((p0/pmid(i,k))**(rair/cpair)*t(i,k)+CA(i,k)*CFt(i,k+1)) &
                       /(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1)) 
            CFq(i,k) = (q(i,k)+CA(i,k)*CFq(i,k+1)) &
                       /(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1))
        end do
      end do

!
! Calculate the updated temperature, specific humidity and horizontal wind
!
! First we need to calculate the updates at the top model level
!
      do i=1,pcols
            dudt(i,1)  = dudt(i,1)+(CFu(i,1)-u(i,1))/dtime
            dvdt(i,1)  = dvdt(i,1)+(CFv(i,1)-v(i,1))/dtime
            u(i,1)    = CFu(i,1)
            v(i,1)    = CFv(i,1)
            dtdt(i,1)  = dtdt(i,1)+(CFt(i,1)*(pmid(i,1)/p0)**(rair/cpair)-t(i,1))/dtime  ! corrected in version 1.3
            t(i,1)    = CFt(i,1)*(pmid(i,1)/p0)**(rair/cpair)
            dqdt(i,1)  = dqdt(i,1)+(CFq(i,1)-q(i,1))/dtime
            q(i,1)  = CFq(i,1)
      end do
!
! Loop over the remaining level
!
      do i=1,pcols
         do k=2,pver
            dudt(i,k)  = dudt(i,k)+(CEm(i,k)*u(i,k-1)+CFu(i,k)-u(i,k))/dtime
            dvdt(i,k)  = dvdt(i,k)+(CEm(i,k)*v(i,k-1)+CFv(i,k)-v(i,k))/dtime
            u(i,k)    = CEm(i,k)*u(i,k-1)+CFu(i,k) 
            v(i,k)    = CEm(i,k)*v(i,k-1)+CFv(i,k)
            dtdt(i,k)  = dtdt(i,k)+((CE(i,k)*t(i,k-1) &
                              *(p0/pmid(i,k-1))**(rair/cpair)+CFt(i,k)) &
                              *(pmid(i,k)/p0)**(rair/cpair)-t(i,k))/dtime 
            t(i,k)    = (CE(i,k)*t(i,k-1)*(p0/pmid(i,k-1))**(rair/cpair)+CFt(i,k)) &
                              *(pmid(i,k)/p0)**(rair/cpair)
            dqdt(i,k)  = dqdt(i,k)+(CE(i,k)*q(i,k-1)+CFq(i,k)-q(i,k))/dtime
            q(i,k)  = CE(i,k)*q(i,k-1)+CFq(i,k)
         end do
      end do

!===============================================================================
!
! Dry Mass Adjustment - THIS PROCESS WILL BE MODEL SPECIFIC 
!
! note: Care needs to be taken to ensure that the model conserves the total
!       dry air mass. Add your own routine here.
!===============================================================================
  !  call physics_dme_adjust(state, tend, qini, dtime)   ! This is for CESM/CAM

      

   return
 end subroutine reed_simple_physics_v5





end module fv_phys_mod
