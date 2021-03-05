!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

! =======================================================================
! Fast Physics Interface
! Developer: Linjiong Zhou
! Last Update: 3/5/2021
! =======================================================================

module fast_phys_mod

    use constants_mod, only: rdgas, grav
    use fv_grid_utils_mod, only: cubed_to_latlon, update_dwinds_phys
    use fv_arrays_mod, only: fv_grid_type, fv_grid_bounds_type, inline_mp_type
    use mpp_domains_mod, only: domain2d, mpp_update_domains
    use fv_timing_mod, only: timing_on, timing_off
    use tracer_manager_mod, only: get_tracer_index
    use field_manager_mod, only: model_atmos
#ifndef DYCORE_SOLO
    use gfdl_mp_mod, only: gfdl_mp_driver, fast_sat_adj
#endif
    
    implicit none
    
    private

    real, parameter :: consv_min = 0.001

    public :: fast_phys

contains

subroutine fast_phys (is, ie, js, je, isd, ied, jsd, jed, km, npx, npy, &
               c2l_ord, mdt, consv, akap, pfull, hs, te0_2d, ua, va, u, &
               v, w, pt, delp, delz, q_con, cappa, q, pkz, te, peln, pe, pk, ps, &
               inline_mp, gridstruct, domain, bd, hydrostatic, do_adiabatic_init, &
               do_inline_mp, do_sat_adj, last_step)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: is, ie, js, je, isd, ied, jsd, jed, km, npx, npy, c2l_ord

    logical, intent (in) :: hydrostatic, do_adiabatic_init, do_inline_mp, do_sat_adj, last_step

    real, intent (in) :: consv, mdt, akap

    real, intent (in), dimension (km) :: pfull

    real, intent (in), dimension (isd:ied, jsd:jed) :: hs

    real, intent (inout), dimension (isd:ied, jsd:jed) :: ps

    real, intent (inout), dimension (is:ie, js:je) :: te0_2d

    real, intent (inout), dimension (is:ie, js:je, km+1) :: pk

    real, intent (inout), dimension (is:ie, km+1, js:je) :: peln

    real, intent (inout), dimension (is:, js:, 1:) :: delz
    
    real, intent (inout), dimension (isd:, jsd:, 1:) :: q_con, cappa, w
    
    real, intent (inout), dimension (isd:ied, jsd:jed, km) :: pt, ua, va, delp

    real, intent (inout), dimension (isd:ied, jsd:jed, km, *) :: q

    real, intent (inout), dimension (isd:ied, jsd:jed+1, km) :: u

    real, intent (inout), dimension (isd:ied+1, jsd:jed, km) :: v

    real, intent (inout), dimension (is-1:ie+1, km+1, js-1:je+1) :: pe

    real, intent (out), dimension (is:ie, js:je, km) :: pkz

    real, intent (out), dimension (isd:ied, jsd:jed, km) :: te

    type (fv_grid_type), intent (in), target :: gridstruct

    type (fv_grid_bounds_type), intent (in) :: bd

    type (domain2d), intent (inout) :: domain

    type (inline_mp_type), intent (inout) :: inline_mp

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i, j, k, kmp
    integer :: sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, cld_amt, ccn_cm3, cin_cm3

    real :: rrg

    real, dimension (is:ie) :: gsize

    real, dimension (is:ie, km) :: q2, q3

    real, dimension (is:ie, km+1) :: phis

    real, allocatable, dimension (:,:) :: dz, wa

    real, allocatable, dimension (:,:,:) :: u_dt, v_dt, dp0, u0, v0
    
    sphum = get_tracer_index (model_atmos, 'sphum')
    liq_wat = get_tracer_index (model_atmos, 'liq_wat')
    ice_wat = get_tracer_index (model_atmos, 'ice_wat')
    rainwat = get_tracer_index (model_atmos, 'rainwat')
    snowwat = get_tracer_index (model_atmos, 'snowwat')
    graupel = get_tracer_index (model_atmos, 'graupel')
    cld_amt = get_tracer_index (model_atmos, 'cld_amt')
    ccn_cm3 = get_tracer_index (model_atmos, 'ccn_cm3')
    cin_cm3 = get_tracer_index (model_atmos, 'cin_cm3')

    rrg = - rdgas / grav

    ! time saving trick
    if (last_step) then
        kmp = 1
    else
        do k = 1, km
            kmp = k
            if (pfull (k) .gt. 50.E2) exit
        enddo
    endif

    !-----------------------------------------------------------------------
    ! Fast Saturation Adjustment >>>
    !-----------------------------------------------------------------------

    ! Note: pt at this stage is T_v
    if (do_adiabatic_init .or. do_sat_adj) then

        call timing_on ('fast_sat_adj')

        allocate (dz (is:ie, kmp:km))

!$OMP parallel do default (none) shared (is, ie, js, je, isd, jsd, kmp, km, te, &
!$OMP                                    delp, hydrostatic, hs, pt, peln, delz, rainwat, &
!$OMP                                    liq_wat, ice_wat, snowwat, graupel, q_con, &
!$OMP                                    sphum, pkz, last_step, consv, te0_2d, gridstruct, &
!$OMP                                    q, mdt, cld_amt, cappa, rrg, akap, ccn_cm3, &
!$OMP                                    cin_cm3, inline_mp) &
!$OMP                           private (q2, q3, gsize, dz)

        do j = js, je

            gsize (is:ie) = sqrt (gridstruct%area_64 (is:ie, j))

            if (ccn_cm3 .gt. 0) then
                q2 (is:ie, kmp:km) = q (is:ie, j, kmp:km, ccn_cm3)
            else
                q2 (is:ie, kmp:km) = 0.0
            endif
            if (cin_cm3 .gt. 0) then
                q3 (is:ie, kmp:km) = q (is:ie, j, kmp:km, cin_cm3)
            else
                q3 (is:ie, kmp:km) = 0.0
            endif
 
            if (.not. hydrostatic) then
                dz (is:ie, kmp:km) = delz (is:ie, j, kmp:km)
            else
                dz (is:ie, kmp:km) = (peln (is:je, kmp:km, j) - peln (is:ie, kmp+1:km+1, j)) * &
                    rdgas * pt (is:ie, j, kmp:km) / grav
            endif

            call fast_sat_adj (abs (mdt), is, ie, kmp, km, hydrostatic, consv .gt. consv_min, &
                     te (is:ie, j, kmp:km), q (is:ie, j, kmp:km, sphum), q (is:ie, j, kmp:km, liq_wat), &
                     q (is:ie, j, kmp:km, rainwat), q (is:ie, j, kmp:km, ice_wat), &
                     q (is:ie, j, kmp:km, snowwat), q (is:ie, j, kmp:km, graupel), &
                     q (is:ie, j, kmp:km, cld_amt), q2 (is:ie, kmp:km), q3 (is:ie, kmp:km), &
                     hs (is:ie, j), dz (is:ie, kmp:km), pt (is:ie, j, kmp:km), &
                     delp (is:ie, j, kmp:km), &
#ifdef USE_COND
                     q_con (is:ie, j, kmp:km), &
#else
                     q_con (isd:, jsd,1:), &
#endif
#ifdef MOIST_CAPPA
                     cappa (is:ie, j, kmp:km), &
#else
                     cappa (isd:, jsd,1:), &
#endif
                     gsize, last_step, inline_mp%cond (is:ie, j), inline_mp%reevap (is:ie, j), &
                     inline_mp%dep (is:ie, j), inline_mp%sub (is:ie, j))

            ! update pkz
            if (.not. hydrostatic) then
#ifdef MOIST_CAPPA
                pkz (is:ie, j, kmp:km) = exp (cappa (is:ie, j, kmp:km) * &
                    log (rrg * delp (is:ie, j, kmp:km) / &
                    delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
#else
                pkz (is:ie, j, kmp:km) = exp (akap * log (rrg * delp (is:ie, j, kmp:km) / &
                    delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
#endif
            endif
 
            if (consv .gt. consv_min) then
                do i = is, ie
                    do k = kmp, km
                        te0_2d (i, j) = te0_2d (i, j) + te (i, j, k)
                    enddo
                enddo
            endif

        enddo

        deallocate (dz)

        call timing_off ('fast_sat_adj')

    endif

    !-----------------------------------------------------------------------
    ! <<< Fast Saturation Adjustment
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Inline GFDL MP >>>
    !-----------------------------------------------------------------------

    if ((.not. do_adiabatic_init) .and. do_inline_mp) then

        call timing_on ('gfdl_mp')

        allocate (u_dt (isd:ied, jsd:jed, km))
        allocate (v_dt (isd:ied, jsd:jed, km))

        do k = 1, km
            do j = jsd, jed
                do i = isd, ied
                    u_dt (i, j, k) = 0.
                    v_dt (i, j, k) = 0.
                enddo
            enddo
        enddo

        ! save D grid u and v
        if (consv .gt. consv_min) then
            allocate (u0 (isd:ied, jsd:jed+1, km))
            allocate (v0 (isd:ied+1, jsd:jed, km))
            u0 = u
            v0 = v
        endif

        ! D grid wind to A grid wind remap
        call cubed_to_latlon (u, v, ua, va, gridstruct, npx, npy, km, 1, gridstruct%grid_type, &
                 domain, gridstruct%bounded_domain, c2l_ord, bd)

        ! save delp
        if (consv .gt. consv_min) then
            allocate (dp0 (isd:ied, jsd:jed, km))
            dp0 = delp
        endif

        allocate (dz (is:ie, kmp:km))
        allocate (wa (is:ie, kmp:km))

!$OMP parallel do default (none) shared (is, ie, js, je, isd, jsd, kmp, km, pe, ua, va, &
!$OMP                                    te, delp, hydrostatic, hs, pt, peln, delz, &
!$OMP                                    rainwat, liq_wat, ice_wat, snowwat, graupel, q_con, &
!$OMP                                    sphum, w, pk, pkz, last_step, consv, te0_2d, &
!$OMP                                    gridstruct, q, mdt, cld_amt, cappa, rrg, akap, &
!$OMP                                    ccn_cm3, cin_cm3, inline_mp, do_inline_mp, ps) &
!$OMP                           private (u_dt, v_dt, q2, q3, gsize, dz, wa)

        do j = js, je

            gsize (is:ie) = sqrt (gridstruct%area_64 (is:ie, j))

            if (ccn_cm3 .gt. 0) then
                q2 (is:ie, kmp:km) = q (is:ie, j, kmp:km, ccn_cm3)
            else
                q2 (is:ie, kmp:km) = 0.0
            endif
            if (cin_cm3 .gt. 0) then
                q3 (is:ie, kmp:km) = q (is:ie, j, kmp:km, cin_cm3)
            else
                q3 (is:ie, kmp:km) = 0.0
            endif
 
            ! note: ua and va are A-grid variables
            ! note: pt is virtual temperature at this point
            ! note: w is vertical velocity (m/s)
            ! note: delz is negative, delp is positive, delz doesn't change in constant volume situation
            ! note: hs is geopotential height (m^2/s^2)
            ! note: the unit of q2 or q3 is #/cm^3
            ! note: the unit of area is m^2
            ! note: the unit of prer, prei, pres, preg is mm/day
            ! note: the unit of cond, dep, reevap, sub is mm/day

            ! save ua, va for wind tendency calculation
            u_dt (is:ie, j, kmp:km) = ua (is:ie, j, kmp:km)
            v_dt (is:ie, j, kmp:km) = va (is:ie, j, kmp:km)

            if (allocated (inline_mp%liq_wat_dt)) inline_mp%liq_wat_dt (is:ie, j, kmp:km) = &
                inline_mp%liq_wat_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, liq_wat)
            if (allocated (inline_mp%ice_wat_dt)) inline_mp%ice_wat_dt (is:ie, j, kmp:km) = &
                inline_mp%ice_wat_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, ice_wat)
            if (allocated (inline_mp%qv_dt)) inline_mp%qv_dt (is:ie, j, kmp:km) = &
                inline_mp%qv_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, sphum)
            if (allocated (inline_mp%ql_dt)) inline_mp%ql_dt (is:ie, j, kmp:km) = &
                inline_mp%ql_dt (is:ie, j, kmp:km) - (q (is:ie, j, kmp:km, liq_wat) + &
                q (is:ie, j, kmp:km, rainwat))
            if (allocated (inline_mp%qi_dt)) inline_mp%qi_dt (is:ie, j, kmp:km) = &
                inline_mp%qi_dt (is:ie, j, kmp:km) - (q (is:ie, j, kmp:km, ice_wat) + &
                q (is:ie, j, kmp:km, snowwat) + q (is:ie, j, kmp:km, graupel))
            if (allocated (inline_mp%qr_dt)) inline_mp%qr_dt (is:ie, j, kmp:km) = &
                inline_mp%qr_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, rainwat)
            if (allocated (inline_mp%qs_dt)) inline_mp%qs_dt (is:ie, j, kmp:km) = &
                inline_mp%qs_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, snowwat)
            if (allocated (inline_mp%qg_dt)) inline_mp%qg_dt (is:ie, j, kmp:km) = &
                inline_mp%qg_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, graupel)
            if (allocated (inline_mp%t_dt)) inline_mp%t_dt (is:ie, j, kmp:km) = &
                inline_mp%t_dt (is:ie, j, kmp:km) - pt (is:ie, j, kmp:km)
            if (allocated (inline_mp%u_dt)) inline_mp%u_dt (is:ie, j, kmp:km) = &
                inline_mp%u_dt (is:ie, j, kmp:km) - ua (is:ie, j, kmp:km)
            if (allocated (inline_mp%v_dt)) inline_mp%v_dt (is:ie, j, kmp:km) = &
                inline_mp%v_dt (is:ie, j, kmp:km) - va (is:ie, j, kmp:km)

            if (.not. hydrostatic) then
                wa (is:ie, kmp:km) = w (is:ie, j, kmp:km)
                dz (is:ie, kmp:km) = delz (is:ie, j, kmp:km)
            else
                dz (is:ie, kmp:km) = (peln (is:je, kmp:km, j) - peln (is:ie, kmp+1:km+1, j)) * &
                    rdgas * pt (is:ie, j, kmp:km) / grav
            endif

#ifndef DYCORE_SOLO
            call gfdl_mp_driver (q (is:ie, j, kmp:km, sphum), q (is:ie, j, kmp:km, liq_wat), &
                     q (is:ie, j, kmp:km, rainwat), q (is:ie, j, kmp:km, ice_wat), &
                     q (is:ie, j, kmp:km, snowwat), q (is:ie, j, kmp:km, graupel), &
                     q (is:ie, j, kmp:km, cld_amt), q2 (is:ie, kmp:km), &
                     q3 (is:ie, kmp:km), pt (is:ie, j, kmp:km), wa (is:ie, kmp:km), &
                     ua (is:ie, j, kmp:km), va (is:ie, j, kmp:km), dz (is:ie, kmp:km), &
                     delp (is:ie, j, kmp:km), gsize, abs (mdt), hs (is:ie, j), &
                     inline_mp%prer (is:ie, j), inline_mp%pres (is:ie, j), &
                     inline_mp%prei (is:ie, j), inline_mp%preg (is:ie, j), hydrostatic, &
                     is, ie, kmp, km, &
#ifdef USE_COND
                     q_con (is:ie, j, kmp:km), &
#else
                     q_con (isd:, jsd,1:), &
#endif
#ifdef MOIST_CAPPA
                     cappa (is:ie, j, kmp:km), &
#else
                     cappa (isd:, jsd,1:), &
#endif
                     consv .gt. consv_min, te (is:ie, j, kmp:km), inline_mp%cond (is:ie, j), &
                     inline_mp%dep (is:ie, j), inline_mp%reevap (is:ie, j), inline_mp%sub (is:ie, j), &
                     last_step, do_inline_mp)
#endif

            if (.not. hydrostatic) then
                w (is:ie, j, kmp:km) = wa (is:ie, kmp:km)
            endif

            ! compute wind tendency at A grid fori D grid wind update
            u_dt (is:ie, j, kmp:km) = (ua (is:ie, j, kmp:km) - u_dt (is:ie, j, kmp:km)) / abs (mdt)
            v_dt (is:ie, j, kmp:km) = (va (is:ie, j, kmp:km) - v_dt (is:ie, j, kmp:km)) / abs (mdt)

            if (allocated (inline_mp%liq_wat_dt)) inline_mp%liq_wat_dt (is:ie, j, kmp:km) = &
                inline_mp%liq_wat_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, liq_wat)
            if (allocated (inline_mp%ice_wat_dt)) inline_mp%ice_wat_dt (is:ie, j, kmp:km) = &
                inline_mp%ice_wat_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, ice_wat)
            if (allocated (inline_mp%qv_dt)) inline_mp%qv_dt (is:ie, j, kmp:km) = &
                inline_mp%qv_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, sphum)
            if (allocated (inline_mp%ql_dt)) inline_mp%ql_dt (is:ie, j, kmp:km) = &
                inline_mp%ql_dt (is:ie, j, kmp:km) + (q (is:ie, j, kmp:km, liq_wat) + &
                q (is:ie, j, kmp:km, rainwat))
            if (allocated (inline_mp%qi_dt)) inline_mp%qi_dt (is:ie, j, kmp:km) = &
                inline_mp%qi_dt (is:ie, j, kmp:km) + (q (is:ie, j, kmp:km, ice_wat) + &
                q (is:ie, j, kmp:km, snowwat) + q (is:ie, j, kmp:km, graupel))
            if (allocated (inline_mp%qr_dt)) inline_mp%qr_dt (is:ie, j, kmp:km) = &
                inline_mp%qr_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, rainwat)
            if (allocated (inline_mp%qs_dt)) inline_mp%qs_dt (is:ie, j, kmp:km) = &
                inline_mp%qs_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, snowwat)
            if (allocated (inline_mp%qg_dt)) inline_mp%qg_dt (is:ie, j, kmp:km) = &
                inline_mp%qg_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, graupel)
            if (allocated (inline_mp%t_dt)) inline_mp%t_dt (is:ie, j, kmp:km) = &
                inline_mp%t_dt (is:ie, j, kmp:km) + pt (is:ie, j, kmp:km)
            if (allocated (inline_mp%u_dt)) inline_mp%u_dt (is:ie, j, kmp:km) = &
                inline_mp%u_dt (is:ie, j, kmp:km) + ua (is:ie, j, kmp:km)
            if (allocated (inline_mp%v_dt)) inline_mp%v_dt (is:ie, j, kmp:km) = &
                inline_mp%v_dt (is:ie, j, kmp:km) + va (is:ie, j, kmp:km)

            ! update pe, peln, pk, ps
            do k = kmp + 1, km + 1
                pe (is:ie, k, j) = pe (is:ie, k-1, j) + delp (is:ie, j, k-1)
                peln (is:ie, k, j) = log (pe (is:ie, k, j))
                pk (is:ie, j, k) = exp (akap * peln (is:ie, k, j))
            enddo

            ps (is:ie, j) = pe (is:ie, km+1, j)

            ! update pkz
            if (.not. hydrostatic) then
#ifdef MOIST_CAPPA
                pkz (is:ie, j, kmp:km) = exp (cappa (is:ie, j, kmp:km) * &
                    log (rrg * delp (is:ie, j, kmp:km) / &
                    delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
#else
                pkz (is:ie, j, kmp:km) = exp (akap * log (rrg * delp (is:ie, j, kmp:km) / &
                    delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
#endif
            endif
 
            if (consv .gt. consv_min) then
                do i = is, ie
                    do k = kmp, km
                        te0_2d (i, j) = te0_2d (i, j) + te (i, j, k)
                    enddo
                enddo
            endif

        enddo

        deallocate (dz)
        deallocate (wa)

        ! Note: (ua, va) are *lat-lon* wind tendenies on cell centers
        if ( gridstruct%square_domain ) then
            call mpp_update_domains (u_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
            call mpp_update_domains (v_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
        else
            call mpp_update_domains (u_dt, domain, complete=.false.)
            call mpp_update_domains (v_dt, domain, complete=.true.)
        endif
        ! update u_dt and v_dt in halo
        call mpp_update_domains (u_dt, v_dt, domain)

        ! update D grid wind
        call update_dwinds_phys (is, ie, js, je, isd, ied, jsd, jed, abs (mdt), u_dt, v_dt, u, v, &
                 gridstruct, npx, npy, km, domain)

        ! update dry total energy
        if (consv .gt. consv_min) then
!$OMP parallel do default (none) shared (is, ie, js, je, km, te0_2d, hydrostatic, delp, &
!$OMP                                    gridstruct, u, v, dp0, u0, v0, hs, delz, w) &
!$OMP                           private (phis)
            do j = js, je
                if (hydrostatic) then
                    do k = 1, km
                        do i = is, ie
                            te0_2d (i, j) = te0_2d (i, j) + delp (i, j, k) * &
                                (0.25 * gridstruct%rsin2 (i, j) * (u (i, j, k) ** 2 + &
                                u (i, j+1, k) ** 2 + v (i, j, k) ** 2 + v (i+1, j, k) ** 2 - &
                                (u (i, j, k) + u (i, j+1, k)) * (v (i, j, k) + v (i+1, j, k)) * &
                                gridstruct%cosa_s (i, j))) - dp0 (i, j, k) * &
                                (0.25 * gridstruct%rsin2 (i, j) * (u0 (i, j, k) ** 2 + &
                                u0 (i, j+1, k) ** 2 + v0 (i, j, k) ** 2 + v0 (i+1, j, k) ** 2 - &
                                (u0 (i, j, k) + u0 (i, j+1, k)) * (v0 (i, j, k) + v0 (i+1, j, k)) * &
                                gridstruct%cosa_s (i, j)))
                        enddo
                    enddo
                else
                    do i = is, ie
                        phis (i, km+1) = hs (i, j)
                    enddo
                    do k = km, 1, -1
                        do i = is, ie
                            phis (i, k) = phis (i, k+1) - grav * delz (i, j, k)
                        enddo
                    enddo
                    do k = 1, km
                        do i = is, ie
                            te0_2d (i, j) = te0_2d (i, j) + delp (i, j, k) * &
                                (0.5 * (phis (i, k) + phis (i, k+1) + w (i, j, k) ** 2 + 0.5 * &
                                gridstruct%rsin2 (i, j) * (u (i, j, k) ** 2 + u (i, j+1, k) ** 2 + &
                                v (i, j, k) ** 2 + v (i+1, j, k) ** 2 - (u (i, j, k) + &
                                u (i, j+1, k)) * (v (i, j, k) + v (i+1, j, k)) * &
                                gridstruct%cosa_s (i, j)))) - dp0 (i, j, k) * &
                                (0.5 * (phis (i, k) + phis (i, k+1) + w (i, j, k) ** 2 + &
                                0.5 * gridstruct%rsin2 (i, j) * (u0 (i, j, k) ** 2 + &
                                u0 (i, j+1, k) ** 2 + v0 (i, j, k) ** 2 + v0 (i+1, j, k) ** 2 - &
                                (u0 (i, j, k) + u0 (i, j+1, k)) * (v0 (i, j, k) + v0 (i+1, j, k)) * &
                                gridstruct%cosa_s (i, j))))
                        enddo
                    enddo
                endif
            enddo
        end if

        deallocate (u_dt)
        deallocate (v_dt)
        if (consv .gt. consv_min) then
            deallocate (u0)
            deallocate (v0)
            deallocate (dp0)
        endif

        call timing_off ('gfdl_mp')

    endif

    !-----------------------------------------------------------------------
    ! <<< Inline GFDL MP
    !-----------------------------------------------------------------------

end subroutine fast_phys

end module fast_phys_mod
