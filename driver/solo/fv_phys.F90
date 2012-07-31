module fv_phys_mod

use constants_mod,         only: grav, rdgas, rvgas, pi
use time_manager_mod,      only: time_type
use lin_cld_microphys_mod, only: lin_cld_microphys_driver, sg_conv, qsmith
use fv_grid_tools_mod,     only: area
use hswf_mod,              only: Held_Suarez_Strat, Held_Suarez_Tend
use fv_sg_mod,             only: fv_dry_conv
use fv_update_phys_mod,    only: fv_update_phys
use fv_timing_mod,         only: timing_on, timing_off
use fv_current_grid_mod, only: nudge_initialized, u0, v0, t0, dp

implicit none

public :: fv_phys, fv_nudge

contains

!-----------------------------------------------------------------------

 subroutine fv_phys(npx, npy, npz, is, ie, js, je, ng, nq,       &
                    u, v, w, pt, q, pe, delp, peln, pkz, pdt,    &
                    ua, va, phis, grid, ak, bk, ks, ps, pk,      &
                    u_srf, v_srf, ts, delz, hydrostatic,         &
                    oro, strat, rayf, p_ref, fv_sg_adj, master,  &
                    do_Held_Suarez, Time, time_total)


    integer, INTENT(IN   ) :: npx, npy, npz
    integer, INTENT(IN   ) :: is, ie, js, je, ng, nq
    integer, INTENT(IN   ) :: fv_sg_adj
    real, INTENT(IN) :: p_ref
    real, INTENT(IN) :: oro(is:ie,js:je)

    real   , INTENT(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
    real   , INTENT(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) ::    w(is-ng:ie+  ng,js-ng:je+  ng,npz)
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
    real   , INTENT(IN   ) :: grid(is-ng:ie+ng,js-ng:je+ng, 1:2)
    real   , INTENT(IN   ) :: ak(npz+1), bk(npz+1)
    integer, INTENT(IN   ) :: ks

    real   , INTENT(IN   )    :: pdt
    logical, INTENT(IN   )    :: strat, rayf, master, do_Held_Suarez
    real, INTENT(inout):: u_srf(is:ie,js:je)
    real, INTENT(inout):: v_srf(is:ie,js:je)
    real, INTENT(inout)::    ts(is:ie,js:je)

    type (time_type), intent(in) :: Time
    real, INTENT(IN), optional:: time_total
    logical, intent(in) ::  hydrostatic
    real, intent(inout) ::  delz(is:ie,js:je,npz)
    real, allocatable:: u_dt(:,:,:), v_dt(:,:,:), t_dt(:,:,:), q_dt(:,:,:,:)
    logical moist_phys

    integer  isd, ied, jsd, jed

#ifdef MARS_GCM
    real, allocatable:: qratio(:,:,:)
    integer :: i, j, k, m
#endif MARS_GCM

    isd = is-ng;   ied = ie + ng
    jsd = js-ng;   jed = je + ng

       allocate ( u_dt(isd:ied,jsd:jed,npz) )
       allocate ( v_dt(isd:ied,jsd:jed,npz) )
       allocate ( t_dt(is:ie,js:je,npz) )
       allocate ( q_dt(is:ie,js:je,npz,nq) )
       u_dt = 0.
       v_dt = 0.
       t_dt = 0.
       q_dt = 0.

#ifdef MARS_GCM
       allocate ( qratio(is:ie,js:je,npz) )

       call timing_on('Mars_PHYS')
       call Mars_phys(npx, npy, npz, is, ie, js, je, ng, nq,   &
                     u_dt, v_dt, t_dt, q_dt, ua, va, pt, q,   &
                     phis, pe, delp, peln, pdt, grid, ak, bk,       &
                     qratio, rayf, master, Time, time_total )

       call timing_off('Mars_PHYS')

!!!      Be sure to set moist_phys to true
       call timing_on('UPDATE_PHYS')

       call fv_update_phys (pdt, is, ie, js, je, isd, ied, jsd, jed, ng, nq,   &
                            u, v, delp, pt, q, ua, va, ps, pe, peln, pk, pkz,  &
                            ak, bk, phis, u_srf, v_srf, ts, delz, hydrostatic, &
                            u_dt, v_dt, t_dt, q_dt, .true., Time )

!-----------------------------------------
! Adjust mass mixing ratio of all tracers:
!!!! May NEED TO FIX THIS?!?!  Perhaps do this for prognostic tracers only, and
!!!!      let the physics code figure out how to handle diagnostic tracers
!-----------------------------------------
     do m=1,nq
       do k=1,npz
         do j=js,je
             do i= is,ie
                q(i,j,k,m) = q(i,j,k,m) / qratio(i,j,k)
             enddo
          enddo
       enddo
     enddo

       deallocate ( u_dt )
       deallocate ( v_dt )
       deallocate ( t_dt )
       deallocate ( q_dt )
       deallocate ( qratio )
       call timing_off('UPDATE_PHYS')

! ------------------------------------------------
#else MARS_GCM

    if( do_Held_Suarez ) then
       moist_phys = .false.
       if ( fv_sg_adj > 0 ) then
          call fv_dry_conv( isd, ied, jsd, jed, is, ie, js, je, npz, nq, pdt,  &
                            fv_sg_adj, delp, pe, peln, pkz, pt, q, ua, va,  &
                            hydrostatic, w, delz, u_dt, v_dt, t_dt, q_dt )

       endif
#ifdef USE_ADJUST_FLD
       call Held_Suarez_Strat(npx, npy, npz, is, ie, js, je, ng, nq,  &
                              u, v, pt, q, pe, delp, peln, pkz, pdt,  &
                              ua, va, grid, ak, bk, ks, strat,  &
                              rayf, master, Time, time_total)
#else
       call Held_Suarez_Tend(npx, npy, npz, is, ie, js, je, ng, nq,   &
                              u, v, pt, q, pe, delp, peln, pkz, pdt,  &
                              ua, va, u_dt, v_dt, t_dt, q_dt, grid,   &
                              delz, hydrostatic, ak, bk, ks,          &
                              strat, rayf, master, Time, time_total)
#endif
    else
       moist_phys = .true.
                                             call timing_on('SIM_PHYS')
       call sim_phys(npx, npy, npz, is, ie, js, je, ng, nq,       &
                     u_dt, v_dt, t_dt, q_dt, u, v, w, ua, va, pt, delz, q, &
                     pe, delp, peln, oro, hydrostatic, pdt, grid, ak, bk, &
                     p_ref, fv_sg_adj, master, Time, time_total)
                                            call timing_off('SIM_PHYS')
    endif

                        call timing_on('UPDATE_PHYS')
    call fv_update_phys (pdt, is, ie, js, je, isd, ied, jsd, jed, ng, nq,   &
                         u, v, delp, pt, q, ua, va, ps, pe, peln, pk, pkz,  &
                         ak, bk, phis, u_srf, v_srf, ts,  &
                         delz, hydrostatic, u_dt, v_dt, t_dt, q_dt, moist_phys, Time, .false.)
                        call timing_off('UPDATE_PHYS')
    deallocate ( u_dt )
    deallocate ( v_dt )
    deallocate ( t_dt )
    deallocate ( q_dt )
#endif MARS_GCM

 end subroutine fv_phys



 subroutine sim_phys(npx, npy, npz, is, ie, js, je, ng, nq,  &
                     u_dt, v_dt, t_dt, q_dt, u, v, w, ua, va,   &
                     pt, delz, q, pe, delp, peln, oro, hydrostatic, &
                     pdt, agrid, ak, bk, p_ref, fv_sg_adj,    &
                     master, Time, time_total)
!-----------------------
! A simple moist physics
!-----------------------


 integer, INTENT(IN) :: npx, npy, npz
 integer, INTENT(IN) :: is, ie, js, je, ng, nq
 integer, INTENT(IN) :: fv_sg_adj
 real   , INTENT(IN) :: pdt
 real   , INTENT(IN) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
 real   , INTENT(IN) :: ak(npz+1), bk(npz+1)
 logical, INTENT(IN) :: master
 real, INTENT(IN):: p_ref
 real, INTENT(IN):: oro(is:ie,js:je)       ! land fraction
 logical, INTENT(IN):: hydrostatic

 type(time_type), intent(in) :: Time
 real, INTENT(IN), optional:: time_total

 real   , INTENT(INOUT) :: u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
 real   , INTENT(INOUT) :: v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
 real, INTENT(INOUT)::   pt(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: delp(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT)::    q(is-ng:ie+ng,js-ng:je+ng,npz, nq)
 real, INTENT(INOUT)::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
 real, INTENT(INOUT):: peln(is  :ie   ,1:npz+1,js  :je  )
 real, INTENT(INOUT):: delz(is   :ie   ,js   :je   ,npz)
 real, intent(inout):: w(is-ng:ie+ng,js-ng:je+ng,npz)

! Tendencies:
 real, INTENT(INOUT):: u_dt(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: v_dt(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
 real, INTENT(INOUT):: q_dt(is:ie,js:je,npz,nq)
 real, INTENT(INOUT):: ua(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: va(is-ng:ie+ng,js-ng:je+ng,npz)
! Local
 logical:: phys_hydrostatic
 real, parameter:: tice = 273.16
 real, dimension(is:ie,js:je,npz):: t3, p3, dz
 real, dimension(is:ie,js:je,npz+1):: zhalf, phalf
 real, dimension(is:ie,js:je,npz  ):: zfull, pfull
 real, dimension(is:ie,js:je,npz,nq):: q3
 real, dimension(is:ie,js:je):: rain, snow, ice, graupel, prec_mp, land
 real, dimension(is:ie,js:je,npz,nq):: qtmp
 real qs(is:ie)
 real pedge(npz+1)
 real pref(npz)
 real sst, sday, shour, rkv, rka, rks, rkt, zvir, rrg, tvm
 real rk1, tmp, rkq, rkd, cooling, frac
 real t00, q00, sst0, rdt, convt, tot_prec
 real fac_sm
 real:: p_conv = 200.E2       ! above which the flow is adiabatic
 logical used
 logical do_lin_cld_microphys
 integer  i,j,k, iq, nq_con, k_mp, k_sg
 integer  seconds, days
 integer  nqv, nql, nqi

   t00  = 200.      ! minimum allowable temperature
   q00  = 1.E-6
   sst0 = 302.      ! K

   phys_hydrostatic = .true.
   do_lin_cld_microphys = .true.

! Factor for Small-Earth Approx.
!  fac_sm = RADIUS / 6371.0e3
   fac_sm = 1.
   rrg  = rdgas / grav
   sday  = 24.*3600.*fac_sm
   shour = 3600.*fac_sm

   rdt = 1. / pdt

   cooling = 2.5 * pdt / sday
!  cooling = 1.5 * pdt / sday

! Drying term:
   rkd = pdt / (120.*shour)      ! sphum
   rkt = pdt / (6.*sday)       ! temp relaxation

#ifdef OLD_PARAM
   rkq = pdt / (3.*shour)      ! sphum
   rk1 = pdt / (3.*shour)      ! sensible heat (from prescribed SST)
   rkv = pdt / (6.*shour)
#else
! TESTING !!!!
   rkq = pdt / 600.      ! sphum
   rk1 = pdt / 600.      ! sensible heat (from prescribed SST)
   rkv = pdt / 10800.
! TESTING !!!!
#endif

   t_dt = 0.;  q_dt = 0.
   u_dt = 0.;  v_dt = 0.

   do k=1,npz+1
      pedge(k) = ak(k) + bk(k)*p_ref
   enddo
   do k=1,npz
      pref(k) = (pedge(k+1)-pedge(k)) / log(pedge(k+1)/pedge(k))
   enddo

   k_sg = 1
   do k=1,npz
      if ( pedge(k) > p_conv ) then
           k_sg = k
           exit
      endif
   enddo

   do k=1,npz+1
      do j=js,je
         do i=is,ie
            phalf(i,j,k) = pe(i,k,j)
         enddo
      enddo
   enddo


!----------------------------------------------
! Apply cooling (i.e., net radiative cooling):
!----------------------------------------------
   do k=1,npz
      do j=js,je
         do i=is,ie
            if ( tmp > t00 ) then
                 t3(i,j,k) = pt(i,j,k) - cooling
            else
! Nudging to t00
                 t3(i,j,k) = (pt(i,j,k)+rkt*t00)/(1.+rkt)
            endif
         enddo 
      enddo 
   enddo      ! k-loop

!----------
! Setup SST:
!----------
   do j=js,je
      do i=is,ie
#ifdef UNIFORM_SST
         sst = sst0
#else
! Neale & Hoskins 2001
        if ( abs(agrid(i,j,2)) <  pi/3. ) then
             tmp = sin( 1.5*agrid(i,j,2) ) ** 2

! QOBS
             sst = 273.16 + 27.*(1.-0.5*(tmp + tmp**2))
! FLAT_SST
!            sst = 273.16 + 27.*(1.-tmp**2)
!----------------------
! Default control case:
!----------------------
!            sst = 273.16 + 27.*(1.-tmp)
         else
             sst = 273.16
         endif
#endif
!------------------------------------------------------
! The following is effectively the sensible heat flux:
!------------------------------------------------------
#ifdef QUICK_SPIN_UP
         t3(i,j,npz) = sst
#else
         t3(i,j,npz) = (t3(i,j,npz)+rk1*sst)/(1.+rk1)
#endif
       enddo
   enddo

   do iq=1,nq
      do k=1,npz
         do j=js,je
            do i=is,ie
               q3(i,j,k,iq) = q(i,j,k,iq)
            enddo
         enddo
      enddo
   enddo

         
! Compute zfull, zhalf
   do j=js,je
      do i=is,ie
         zhalf(i,j,npz+1) = 0.
      enddo
   enddo

   zvir = rvgas/rdgas - 1.

   do k=npz,1,-1
      do j=js,je
         do i=is,ie
#ifdef PHYS_NON_HYDRO
               dz(i,j,k) =  delz(i,j,k)
               p3(i,j,k) = -rrg*delp(i,j,k)/delz(i,j,k)*t3(i,j,k)*(1.+zvir*q3(i,j,k,1))
            zhalf(i,j,k) = zhalf(i,j,k+1) - dz(i,j,k)
            zfull(i,j,k) = 0.5*(zhalf(i,j,k)+zhalf(i,j,k+1))
#else
                     tvm = rrg*t3(i,j,k)*(1.+zvir*q3(i,j,k,1))
               dz(i,j,k) = -tvm*(peln(i,k+1,j)-peln(i,k,j))
               p3(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
            zfull(i,j,k) = zhalf(i,j,k+1) + tvm*(1.-phalf(i,j,k)/p3(i,j,k))
            zhalf(i,j,k) = zhalf(i,j,k+1) - dz(i,j,k)
#endif
         enddo
      enddo
   enddo

!--------------------------
! Moist Physical processes:
!--------------------------

! *** Drying (relaxation) term, at the model top ***
      do j=js,je
         do i=is,ie
            q3(i,j,1,1) =  q3(i,j,1,1)/2.
            q3(i,j,2,1) = (q3(i,j,2,1)+rkd*q00)/(1.+rkd)
         enddo 
      enddo 

!--------------------
! Latent heat fluxes:
!--------------------
      do j=js,je
         call qsmith(ie-is+1, 1, 1, t3(is:ie,j,npz), p3(is:ie,j,npz), q3(is:ie,j,npz,1), qs)
         do i=is,ie
            rkt = rkq * max( 0.25, sqrt(ua(i,j,npz)**2+va(i,j,npz)**2) )
#ifdef QUICK_SPIN_UP
           q3(i,j,npz,1) = qs(i)
           q3(i,j,npz-1,1) = qs(i)
           q3(i,j,npz-2,1) = qs(i)
#else
           q3(i,j,npz,1) = (q3(i,j,npz,1) + rkt*qs(i)) / (1.+rkt)
! Top layer "sink"
!          q3(i,j,1,1) =  q3(i,j,1,1) / (1.+rkd)
#endif
         enddo
      enddo   ! j-loop

   land = 0.

   if ( fv_sg_adj > 0 ) then
         nq_con = nq
         nqv = 1
         nql = 2
         nqi = 3
!        k_sg = 18
         call sg_conv(is, ie, js, je,  is, ie, js, je,    &
                      is, ie, js, je, npz, nq_con, pdt, fv_sg_adj,  &
                      delp(is:ie,js:je,1:npz), phalf,    &
                      p3, zfull, zhalf, t3, q3,    &
                      ua(is:ie,js:je,1:npz),       &
                      va(is:ie,js:je,1:npz),       &
                       w(is:ie,js:je,1:npz),       &
                      u_dt(is:ie,js:je,1:npz),     &
                      v_dt(is:ie,js:je,1:npz),          &
                      t_dt, q_dt, k_sg, nqv, nql, nqi,  &
                      hydrostatic, .true.)
       do k=1,npz
          do j=js,je
             do i=is,ie
                u_dt(i,j,k) = (u_dt(i,j,k) - ua(i,j,k))*rdt                
                v_dt(i,j,k) = (v_dt(i,j,k) - va(i,j,k))*rdt                
                  t3(i,j,k) = t_dt(i,j,k)   ! updated temp
                t_dt(i,j,k) = (t3(i,j,k) - pt(i,j,k)) * rdt
                  dz(i,j,k) = dz(i,j,k) * t3(i,j,k)/pt(i,j,k)  ! isobaric heating
             enddo
          enddo
       enddo

       do iq=1,nq
       do k=1,npz
          do j=js,je
             do i=is,ie
                  q3(i,j,k,iq) = q_dt(i,j,k,iq)
                q_dt(i,j,k,iq) = (q3(i,j,k,iq) - q(i,j,k,iq))*rdt
             enddo
          enddo
       enddo
       enddo
 else
       do k=1,npz
          do j=js,je
             do i=is,ie
                u_dt(i,j,k) = 0.
                v_dt(i,j,k) = 0.
                t_dt(i,j,k) = (t3(i,j,k)-pt(i,j,k))*rdt
             enddo
          enddo
       enddo

       do iq=1,nq
       do k=1,npz
          do j=js,je
             do i=is,ie
                q_dt(i,j,k,iq) = (q3(i,j,k,iq)-q(i,j,k,iq))*rdt
             enddo
          enddo
       enddo
       enddo
 endif    ! fv_sg_adj


   if ( do_lin_cld_microphys ) then
!---------------------------------------
! A 6-class cloud microphysics 
!---------------------------------------)
      k_mp = 1
      call timing_on ('LIN_CLD_MICROPHYS')
! ice 3
! Rain: 4
! snow: 5
      call lin_cld_microphys_driver(q3(is:ie,js:je,1:npz,1), q3(is:ie,js:je,1:npz,2), q3(is:ie,js:je,1:npz,4),  &
                                    q3(is:ie,js:je,1:npz,3), q3(is:ie,js:je,1:npz,5), q3(is:ie,js:je,1:npz,6),  &
                                    q3(is:ie,js:je,1:npz,7),                           &
                                  q_dt(is:ie,js:je,1:npz,1), q_dt(is:ie,js:je,1:npz,2), q_dt(is:ie,js:je,1:npz,4), &
                                  q_dt(is:ie,js:je,1:npz,3), q_dt(is:ie,js:je,1:npz,5), q_dt(is:ie,js:je,1:npz,6), &
                                  q_dt(is:ie,js:je,1:npz,7), t_dt, t3, p3, dz,      &
                                  delp(is:ie,js:je,1:npz), area(is:ie,js:je),  &
                                  pdt, land, rain, snow, ice, graupel, hydrostatic, phys_hydrostatic, &
                                  is,ie, js,je, 1,npz, k_mp,npz, Time)
      call timing_off('LIN_CLD_MICROPHYS')

#ifdef DEBUG_MP
      prec_mp(:,:) = rain(:,:) + snow(:,:) + ice(:,:) + graupel(:,:)
      call prt_maxmin('prec_mp', prec_mp, is, ie, js, je, 0, 1, 1., master)
      call prt_maxmin('Rain', rain, is, ie, js, je, 0, 1, 1., master)
      call prt_maxmin('Snow', snow, is, ie, js, je, 0, 1, 1., master)
      call prt_maxmin('ice',   ice, is, ie, js, je, 0, 1, 1., master)
      call prt_maxmin('Graupel', graupel, is, ie, js, je, 0, 1, 1., master)
#endif
    endif

!-------------------------
! one-layer (surface) drag 
!-------------------------
   k=npz
   do j=js,je
      do i=is,ie
                frac = rkv*max(0.5, sqrt(ua(i,j,npz)**2+va(i,j,npz)**2))
         u_dt(i,j,k) = u_dt(i,j,k) - ua(i,j,k)*frac/(1.+frac) * rdt
         v_dt(i,j,k) = v_dt(i,j,k) - va(i,j,k)*frac/(1.+frac) * rdt
      enddo
   enddo

 end subroutine sim_phys

 subroutine fv_nudge( npz, is, ie, js, je, ng, u, v, delp, pt, dt,    &
                      tau_winds, tau_press, tau_temp )
!
! Nudge the prognostic varaibles toward the IC
! This is only useful for generating balanced steady state IC

    real   , INTENT(IN   ) :: tau_winds, tau_press, tau_temp, dt
    integer, INTENT(IN   ) :: npz, is, ie, js, je, ng
    real   , INTENT(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
    real   , INTENT(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
    real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
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
         endif

         if ( tau_press > 0. ) then
             c_p = dt/tau_press
             do j=js,je
                do i=is,ie
                   delp(i,j,k) = (delp(i,j,k)+c_p*dp(i,j,k)) / (1.+c_p)
                enddo
             enddo
         endif

         if ( tau_temp > 0. ) then
             c_t = dt/tau_temp
             do j=js,je
                do i=is,ie
                   pt(i,j,k) = (pt(i,j,k)+c_t*t0(i,j,k)) / (1.+c_t)
                enddo
             enddo
         endif
      enddo

 end subroutine fv_nudge

end module fv_phys_mod
