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
!* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
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

    use constants_mod, only: rdgas, grav, kappa, cp_air
    use fv_grid_utils_mod, only: cubed_to_latlon, update_dwinds_phys
    use fv_arrays_mod, only: fv_grid_type, fv_grid_bounds_type, inline_mp_type, &
                             inline_edmf_type, inline_sas_type, inline_gwd_type
    use mpp_domains_mod, only: domain2d, mpp_update_domains
    use fv_timing_mod, only: timing_on, timing_off
    use tracer_manager_mod, only: get_tracer_index
    use field_manager_mod, only: model_atmos
    use gfdl_mp_mod, only: gfdl_mp_driver, fast_sat_adj, c_liq, c_ice, cv_air, cv_vap
    use sa_sas_mod, only: sa_sas_deep, sa_sas_shal
    use sa_tke_edmf_mod, only: sa_tke_edmf_sfc, sa_tke_edmf_pbl
    use sa_gwd_mod, only: sa_gwd_oro, sa_gwd_cnv
    
    implicit none
    
    private

    real, parameter :: consv_min = 0.001

    public :: fast_phys

contains

subroutine fast_phys (is, ie, js, je, isd, ied, jsd, jed, km, npx, npy, nq, &
               c2l_ord, mdt, consv, akap, ptop, pfull, hs, te0_2d, ua, va, u, &
               v, w, omga, pt, delp, delz, q_con, cappa, q, pkz, te, peln, pe, pk, ps, r_vir, &
               inline_mp, inline_edmf, inline_sas, inline_gwd, gridstruct, domain, bd, hydrostatic, do_adiabatic_init, &
               do_inline_mp, do_inline_edmf, do_inline_sas, do_inline_gwd, do_sat_adj, last_step)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: is, ie, js, je, isd, ied, jsd, jed, km, npx, npy, nq, c2l_ord

    logical, intent (in) :: hydrostatic, do_adiabatic_init, do_inline_mp, do_inline_sas
    logical, intent (in) :: do_inline_edmf, do_inline_gwd, do_sat_adj, last_step

    real, intent (in) :: consv, mdt, akap, r_vir, ptop

    real, intent (in), dimension (km) :: pfull

    real, intent (in), dimension (isd:ied, jsd:jed) :: hs

    real, intent (inout), dimension (isd:ied, jsd:jed) :: ps

    real, intent (inout), dimension (is:ie, js:je) :: te0_2d

    real, intent (inout), dimension (is:ie, js:je, km+1) :: pk

    real, intent (inout), dimension (is:ie, km+1, js:je) :: peln

    real, intent (inout), dimension (is:, js:, 1:) :: delz
    
    real, intent (inout), dimension (isd:, jsd:, 1:) :: q_con, cappa, w
    
    real, intent (inout), dimension (isd:ied, jsd:jed, km) :: pt, ua, va, delp, omga

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
    type (inline_edmf_type), intent (inout) :: inline_edmf
    type (inline_sas_type), intent (inout) :: inline_sas
    type (inline_gwd_type), intent (inout) :: inline_gwd

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i, j, k, kr, kmp, ncld, ntke
    integer :: sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, cld_amt, ccn_cm3, cin_cm3, aerosol

    real :: rrg

    integer, dimension (is:ie, js:je) :: ktop, kbot, kcnv

    real, dimension (is:ie) :: gsize, dqv, dql, dqi, dqr, dqs, dqg, ps_dt, q_liq, q_sol, c_moist

    real, dimension (is:ie, js:je) :: cumabs

    real, dimension (is:ie, km) :: q2, q3, qliq, qsol, cvm

    real, dimension (is:ie, km+1) :: phis

    integer, allocatable, dimension (:) :: lsm, kinver, vegtype

    real, allocatable, dimension (:) :: rn, rb, u10m, v10m, sigmaf
    real, allocatable, dimension (:) :: stress, wind, tmp, qsurf

    real, allocatable, dimension (:,:) :: dz, zm, zi, wa, dp, pm, pi, pmk, pik, qv, ql
    real, allocatable, dimension (:,:) :: ta, uu, vv, ww, radh

    real, allocatable, dimension (:,:,:) :: u_dt, v_dt, dp0, u0, v0, qa
    
    sphum = get_tracer_index (model_atmos, 'sphum')
    liq_wat = get_tracer_index (model_atmos, 'liq_wat')
    ice_wat = get_tracer_index (model_atmos, 'ice_wat')
    rainwat = get_tracer_index (model_atmos, 'rainwat')
    snowwat = get_tracer_index (model_atmos, 'snowwat')
    graupel = get_tracer_index (model_atmos, 'graupel')
    cld_amt = get_tracer_index (model_atmos, 'cld_amt')
    ccn_cm3 = get_tracer_index (model_atmos, 'ccn_cm3')
    cin_cm3 = get_tracer_index (model_atmos, 'cin_cm3')
    aerosol = get_tracer_index (model_atmos, 'aerosol')

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
    if (do_adiabatic_init .or. (.not. do_inline_mp) .or. do_sat_adj) then

        call timing_on ('fast_sat_adj')

        allocate (dz (is:ie, kmp:km))

!$OMP parallel do default (none) shared (is, ie, js, je, isd, jsd, kmp, km, te, &
!$OMP                                    delp, hydrostatic, hs, pt, peln, delz, rainwat, &
!$OMP                                    liq_wat, ice_wat, snowwat, graupel, q_con, &
!$OMP                                    sphum, pkz, last_step, consv, te0_2d, gridstruct, &
!$OMP                                    q, mdt, cld_amt, cappa, rrg, akap, ccn_cm3, &
!$OMP                                    cin_cm3, aerosol, inline_mp, do_sat_adj) &
!$OMP                           private (q2, q3, gsize, dz)

        do j = js, je

            gsize (is:ie) = sqrt (gridstruct%area_64 (is:ie, j))

            if (aerosol .gt. 0) then
                q2 (is:ie, kmp:km) = q (is:ie, j, kmp:km, aerosol)
            elseif (ccn_cm3 .gt. 0) then
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
                dz (is:ie, kmp:km) = (peln (is:ie, kmp:km, j) - peln (is:ie, kmp+1:km+1, j)) * &
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
                     q_con (isd:, jsd, 1:), &
#endif
#ifdef MOIST_CAPPA
                     cappa (is:ie, j, kmp:km), &
#else
                     cappa (isd:, jsd, 1:), &
#endif
                     gsize, last_step, inline_mp%cond (is:ie, j), inline_mp%reevap (is:ie, j), &
                     inline_mp%dep (is:ie, j), inline_mp%sub (is:ie, j), do_sat_adj)

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
    ! Inline SA-TKE-EDMF >>>
    !-----------------------------------------------------------------------

    if ((.not. do_adiabatic_init) .and. do_inline_edmf) then

        call timing_on ('sa_tke_edmf')

        ntke = get_tracer_index (model_atmos, 'sgs_tke')

        allocate (lsm (is:ie))
        allocate (kinver (is:ie))

        allocate (dz (is:ie, 1:km))
        allocate (zm (is:ie, 1:km))
        allocate (zi (is:ie, 1:km+1))
        allocate (dp (is:ie, 1:km))
        allocate (pm (is:ie, 1:km))
        allocate (pi (is:ie, 1:km+1))
        allocate (pmk (is:ie, 1:km))
        allocate (pik (is:ie, 1:km+1))

        allocate (ta (is:ie, 1:km))
        allocate (uu (is:ie, 1:km))
        allocate (vv (is:ie, 1:km))
        allocate (qa (is:ie, 1:km, 1:nq))

        allocate (radh (is:ie, 1:km))
        allocate (rb (is:ie))
        allocate (u10m (is:ie))
        allocate (v10m (is:ie))
        allocate (stress (is:ie))
        allocate (wind (is:ie))
        allocate (sigmaf (is:ie))
        allocate (vegtype (is:ie))
        allocate (qsurf (is:ie))

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

!$OMP parallel do default (none) shared (is, ie, js, je, isd, jsd, km, nq, pe, ua, va, &
!$OMP                                    te, delp, hydrostatic, hs, pt, peln, delz, &
!$OMP                                    rainwat, liq_wat, ice_wat, snowwat, graupel, &
!$OMP                                    sphum, pk, pkz, consv, te0_2d, gridstruct, q, &
!$OMP                                    mdt, cappa, rrg, akap, r_vir, ps, &
!$OMP                                    u_dt, v_dt, ptop, ntke, inline_edmf) &
!$OMP                           private (gsize, dz, lsm, zi, pi, pik, pmk, &
!$OMP                                    zm, dp, pm, ta, uu, vv, qliq, qsol, qa, &
!$OMP                                    radh, rb, u10m, v10m, sigmaf, vegtype, q_liq, &
!$OMP                                    stress, wind, kinver, qsurf, q_sol, c_moist, &
!$OMP                                    cvm, kr, dqv, dql, dqi, dqr, dqs, dqg, ps_dt)

        do j = js, je
 
            gsize (is:ie) = sqrt (gridstruct%area_64 (is:ie, j))

            ! save ua, va for wind tendency calculation
            u_dt (is:ie, j, 1:km) = ua (is:ie, j, 1:km)
            v_dt (is:ie, j, 1:km) = va (is:ie, j, 1:km)

            lsm = 0
            kinver = km

            ! These need to be reviewed later
            qsurf = 0.0

            if (consv .gt. consv_min) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                cvm = (1 - (q (is:ie, j, 1:km, sphum) + qliq + qsol)) * cv_air + &
                    q (is:ie, j, 1:km, sphum) * cv_vap + qliq * c_liq + qsol * c_ice
                te (is:ie, j, 1:km) = - cvm * pt (is:ie, j, 1:km) / ((1. + r_vir * q (is:ie, j, 1:km, sphum)) * &
                    (1. - (qliq + qsol))) * delp (is:ie, j, 1:km)
            endif

            zi (is:ie, 1) = 0.0
            pi (is:ie, km+1) = ptop
            pik (is:ie, km+1) = exp (kappa * log (pi (is:ie, km+1) * 1.e-5))
            do k = 1, km
                kr = km - k + 1
                dp (is:ie, k) = delp (is:ie, j, kr)
                pi (is:ie, kr) = pi (is:ie, kr+1) + delp (is:ie, j, k)
                pik (is:ie, kr) = exp (kappa * log (pi (is:ie, kr) * 1.e-5))
                if (.not. hydrostatic) then
                    pm (is:ie, k) = - dp (is:ie, k) / delz (is:ie, j, kr) * &
                        rdgas * pt (is:ie, j, kr) / grav
                    dz (is:ie, k) = delz (is:ie, j, kr)
                else
                    pm (is:ie, k) = dp (is:ie, k) / (peln (is:ie, kr+1, j) - peln (is:ie, kr, j))
                    dz (is:ie, k) = (peln (is:ie, kr, j) - peln (is:ie, kr+1, j)) * &
                        rdgas * pt (is:ie, j, kr) / grav
                endif
                pmk (is:ie, k) = exp (kappa * log (pm (is:ie, k) * 1.e-5))
                zi (is:ie, k+1) = zi (is:ie, k) - dz (is:ie, k) * grav
                if (k .eq. 1) then
                    zm (is:ie, k) = - 0.5 * dz (is:ie, k) * grav
                else
                    zm (is:ie, k) = zm (is:ie, k-1) - 0.5 * (dz (is:ie, k-1) + dz (is:ie, k)) * grav
                endif
                ta (is:ie, k) = pt (is:ie, j, kr)
                uu (is:ie, k) = ua (is:ie, j, kr)
                vv (is:ie, k) = va (is:ie, j, kr)
                qa (is:ie, k, 1:nq) = q (is:ie, j, kr, 1:nq)
                radh (is:ie, k) = inline_edmf%radh (is:ie, j, kr)
            enddo

            do i = is, ie
                if (hs (i, j) .gt. 0) lsm (i) = 1
                sigmaf (i) = max (inline_edmf%vfrac (i, j), 0.01)
                vegtype (i) = int (inline_edmf%vtype (i, j) + 0.5)
            enddo

            if (inline_edmf%sfc_cpl) then
                u10m = inline_edmf%u10m (is:ie, j)
                v10m = inline_edmf%v10m (is:ie, j)
                rb = inline_edmf%rb (is:ie, j)
                stress = inline_edmf%stress (is:ie, j)
                wind = inline_edmf%wind (is:ie, j)
            else
                call sa_tke_edmf_sfc (ie-is+1, pi (is:ie, 1), uu (is:ie, 1), vv (is:ie, 1), &
                    ta (is:ie, 1), qa (is:ie, 1, 1), inline_edmf%tsfc (is:ie, j), qsurf, &
                    pm (is:ie, 1), pik (is:ie, 1) / pmk (is:ie, 1), &
                    inline_edmf%evap (is:ie, j), inline_edmf%ffmm (is:ie, j), inline_edmf%ffhh (is:ie, j), &
                    zm (is:ie, 1) / grav, inline_edmf%snowd (is:ie, j), &
                    inline_edmf%zorl (is:ie, j), lsm, inline_edmf%uustar (is:ie, j), &
                    sigmaf, vegtype, inline_edmf%shdmax (is:ie, j), &
                    u10m_out = u10m, v10m_out = v10m, rb_out = rb, &
                    stress_out = stress, wind_out = wind)
            endif

            call sa_tke_edmf_pbl (ie-is+1, km, nq, liq_wat, ice_wat, ntke, &
                abs (mdt), uu, vv, ta, qa, gsize, lsm, &
                radh, rb, inline_edmf%zorl (is:ie, j), u10m, v10m, &
                inline_edmf%ffmm (is:ie, j), inline_edmf%ffhh (is:ie, j), &
                inline_edmf%tsfc (is:ie, j), inline_edmf%hflx (is:ie, j), &
                inline_edmf%evap (is:ie, j), stress, wind, kinver, &
                pik (is:ie, 1), dp, pi, pm, pmk, zi, zm, &
                inline_edmf%hpbl (is:ie, j), inline_edmf%kpbl (is:ie, j), &
                inline_edmf%dusfc (is:ie, j), inline_edmf%dvsfc (is:ie, j), &
                inline_edmf%dtsfc (is:ie, j), inline_edmf%dqsfc (is:ie, j))

            do k = 1, km
                kr = km - k + 1
                q (is:ie, j, kr, ntke) = qa (is:ie, k, ntke)
                dqv = qa (is:ie, k, sphum  ) - q (is:ie, j, kr, sphum  )
                dql = qa (is:ie, k, liq_wat) - q (is:ie, j, kr, liq_wat)
                dqi = qa (is:ie, k, ice_wat) - q (is:ie, j, kr, ice_wat)
                dqr = qa (is:ie, k, rainwat) - q (is:ie, j, kr, rainwat)
                dqs = qa (is:ie, k, snowwat) - q (is:ie, j, kr, snowwat)
                dqg = qa (is:ie, k, graupel) - q (is:ie, j, kr, graupel)
                ps_dt = 1 + dqv + dql + dqi + dqr + dqs + dqg
                q (is:ie, j, kr, sphum  ) = qa (is:ie, k, sphum  ) / ps_dt
                q (is:ie, j, kr, liq_wat) = qa (is:ie, k, liq_wat) / ps_dt
                q (is:ie, j, kr, ice_wat) = qa (is:ie, k, ice_wat) / ps_dt
                q (is:ie, j, kr, rainwat) = qa (is:ie, k, rainwat) / ps_dt
                q (is:ie, j, kr, snowwat) = qa (is:ie, k, snowwat) / ps_dt
                q (is:ie, j, kr, graupel) = qa (is:ie, k, graupel) / ps_dt
                delp (is:ie, j, kr) = delp (is:ie, j, kr) * ps_dt
                q_liq = q (is:ie, j, kr, liq_wat) + q (is:ie, j, kr, rainwat)
                q_sol = q (is:ie, j, kr, ice_wat) + q (is:ie, j, kr, snowwat) + q (is:ie, j, kr, graupel)
                c_moist = (1 - (q (is:ie, j, kr, sphum) + q_liq + q_sol)) * cv_air + &
                    q (is:ie, j, kr, sphum) * cv_vap + q_liq * c_liq + q_sol * c_ice
                ps_dt = pt (is:ie, j, kr)
                pt (is:ie, j, kr) = pt (is:ie, j, kr) + (ta (is:ie, k) - pt (is:ie, j, kr)) * cp_air / c_moist
                ua (is:ie, j, kr) = uu (is:ie, k)
                va (is:ie, j, kr) = vv (is:ie, k)
            enddo

            ! compute wind tendency at A grid fori D grid wind update
            u_dt (is:ie, j, 1:km) = (ua (is:ie, j, 1:km) - u_dt (is:ie, j, 1:km)) / abs (mdt)
            v_dt (is:ie, j, 1:km) = (va (is:ie, j, 1:km) - v_dt (is:ie, j, 1:km)) / abs (mdt)

            ! update pe, peln, pk, ps
            do k = 2, km + 1
                pe (is:ie, k, j) = pe (is:ie, k-1, j) + delp (is:ie, j, k-1)
                peln (is:ie, k, j) = log (pe (is:ie, k, j))
                pk (is:ie, j, k) = exp (akap * peln (is:ie, k, j))
            enddo

            ps (is:ie, j) = pe (is:ie, km+1, j)

            ! update pkz
            if (.not. hydrostatic) then
#ifdef MOIST_CAPPA
                pkz (is:ie, j, 1:km) = exp (cappa (is:ie, j, 1:km) * &
                    log (rrg * delp (is:ie, j, 1:km) / &
                    delz (is:ie, j, 1:km) * pt (is:ie, j, 1:km)))
#else
                pkz (is:ie, j, 1:km) = exp (akap * log (rrg * delp (is:ie, j, 1:km) / &
                    delz (is:ie, j, 1:km) * pt (is:ie, j, 1:km)))
#endif
            endif
 
            if (consv .gt. consv_min) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                cvm = (1 - (q (is:ie, j, 1:km, sphum) + qliq + qsol)) * cv_air + &
                    q (is:ie, j, 1:km, sphum) * cv_vap + qliq * c_liq + qsol * c_ice
                te (is:ie, j, 1:km) = te (is:ie, j, 1:km) + &
                    cvm * pt (is:ie, j, 1:km) / ((1. + r_vir * q (is:ie, j, 1:km, sphum)) * &
                    (1. - (qliq + qsol))) * delp (is:ie, j, 1:km)
                do k = 1, km
                    te0_2d (is:ie, j) = te0_2d (is:ie, j) + te (is:ie, j, k)
                enddo
            endif

        enddo

        deallocate (lsm)
        deallocate (kinver)

        deallocate (dz)
        deallocate (zm)
        deallocate (zi)
        deallocate (dp)
        deallocate (pm)
        deallocate (pi)
        deallocate (pmk)
        deallocate (pik)

        deallocate (ta)
        deallocate (uu)
        deallocate (vv)
        deallocate (qa)

        deallocate (radh)
        deallocate (rb)
        deallocate (u10m)
        deallocate (v10m)
        deallocate (stress)
        deallocate (wind)
        deallocate (sigmaf)
        deallocate (vegtype)
        deallocate (qsurf)

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

        deallocate (u_dt)
        deallocate (v_dt)

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

        if (consv .gt. consv_min) then
            deallocate (u0)
            deallocate (v0)
            deallocate (dp0)
        endif

        call timing_off ('sa_tke_edmf')

    endif

    !-----------------------------------------------------------------------
    ! <<< Inline SA-TKE-EDMF
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Inline SA-SAS >>>
    !-----------------------------------------------------------------------

    if ((.not. do_adiabatic_init) .and. do_inline_sas) then

        call timing_on ('sa_sas')

        allocate (lsm (is:ie))

        allocate (rn (is:ie))
        allocate (tmp (is:ie))

        allocate (dz (is:ie, 1:km))
        allocate (zm (is:ie, 1:km))
        allocate (dp (is:ie, 1:km))
        allocate (pm (is:ie, 1:km))

        allocate (ta (is:ie, 1:km))
        allocate (qv (is:ie, 1:km))
        allocate (ql (is:ie, 1:km))
        allocate (uu (is:ie, 1:km))
        allocate (vv (is:ie, 1:km))
        allocate (ww (is:ie, 1:km))

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

!$OMP parallel do default (none) shared (is, ie, js, je, isd, jsd, km, pe, ua, va, &
!$OMP                                    te, delp, hydrostatic, hs, pt, peln, delz, omga, &
!$OMP                                    rainwat, liq_wat, ice_wat, snowwat, graupel, &
!$OMP                                    sphum, pk, pkz, consv, te0_2d, gridstruct, q, &
!$OMP                                    mdt, cappa, rrg, akap, r_vir, inline_sas, ps, &
!$OMP                                    u_dt, v_dt, kbot, ktop, kcnv, cumabs, inline_edmf) &
!$OMP                           private (gsize, dz, lsm, rn, tmp, q_liq, q_sol, &
!$OMP                                    zm, dp, pm, qv, ql, ta, uu, vv, ww, ncld, qliq, qsol, &
!$OMP                                    cvm, kr, dqv, dql, ps_dt, c_moist)

        do j = js, je
 
            gsize (is:ie) = sqrt (gridstruct%area_64 (is:ie, j))

            ! save ua, va for wind tendency calculation
            u_dt (is:ie, j, 1:km) = ua (is:ie, j, 1:km)
            v_dt (is:ie, j, 1:km) = va (is:ie, j, 1:km)

            rn = 0.0
            ktop (is:ie, j) = 1
            kbot (is:ie, j) = km
            kcnv (is:ie, j) = 0
            lsm = 0
            ncld = 1

            if (consv .gt. consv_min) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                cvm = (1 - (q (is:ie, j, 1:km, sphum) + qliq + qsol)) * cv_air + &
                    q (is:ie, j, 1:km, sphum) * cv_vap + qliq * c_liq + qsol * c_ice
                te (is:ie, j, 1:km) = - cvm * pt (is:ie, j, 1:km) / ((1. + r_vir * q (is:ie, j, 1:km, sphum)) * &
                    (1. - (qliq + qsol))) * delp (is:ie, j, 1:km)
            endif

            do k = 1, km
                kr = km - k + 1
                dp (is:ie, k) = delp (is:ie, j, kr)
                if (.not. hydrostatic) then
                    pm (is:ie, k) = - dp (is:ie, k) / delz (is:ie, j, kr) * &
                        rdgas * pt (is:ie, j, kr) / grav
                    dz (is:ie, k) = delz (is:ie, j, kr)
                else
                    pm (is:ie, k) = dp (is:ie, k) / (peln (is:ie, kr+1, j) - peln (is:ie, kr, j))
                    dz (is:ie, k) = (peln (is:ie, kr, j) - peln (is:ie, kr+1, j)) * &
                        rdgas * pt (is:ie, j, kr) / grav
                endif
                if (k .eq. 1) then
                    zm (is:ie, k) = - 0.5 * dz (is:ie, k) * grav
                else
                    zm (is:ie, k) = zm (is:ie, k-1) - 0.5 * (dz (is:ie, k-1) + dz (is:ie, k)) * grav
                endif
                qv (is:ie, k) = q (is:ie, j, kr, sphum)
                ql (is:ie, k) = q (is:ie, j, kr, liq_wat)
                ta (is:ie, k) = pt (is:ie, j, kr)
                uu (is:ie, k) = ua (is:ie, j, kr)
                vv (is:ie, k) = va (is:ie, j, kr)
                ww (is:ie, k) = omga (is:ie, j, kr)
            enddo
  
            do i = is, ie
                if (hs (i, j) .gt. 0) lsm (i) = 1
            enddo

            call sa_sas_deep (ie-is+1, km, abs (mdt), dp, pm, pe (is:ie, km+1, j), zm, ql, &
                qv, ta, uu, vv, rn, kbot (is:ie, j), ktop (is:ie, j), kcnv (is:ie, j), &
                lsm, gsize, ww, ncld)

            inline_sas%prec (is:ie, j) = inline_sas%prec (is:ie, j) + rn
  
            call sa_sas_shal (ie-is+1, km, abs (mdt), dp, pm, pe (is:ie, km+1, j), zm, ql, &
                qv, ta, uu, vv, rn, kbot (is:ie, j), ktop (is:ie, j), kcnv (is:ie, j), &
                lsm, gsize, ww, ncld, inline_edmf%hpbl (is:ie, j))
  
            inline_sas%prec (is:ie, j) = inline_sas%prec (is:ie, j) + rn

            cumabs (is:ie, j) = 0.0
            tmp (is:ie) = 0.0
            do k = 1, km
                kr = km - k + 1
                do i = is, ie
                    if (k .ge. kbot (i, j) .and. k .le. ktop (i, j)) then
                        cumabs (i, j) = cumabs (i, j) + (ta (i, k) - pt (i, j, kr)) * dp (i, k)
                        tmp (i) = tmp (i) + dp (i, k)
                    endif
                enddo
            enddo
            do i = is, ie
                if (tmp (i) .gt. 0.0) cumabs (i, j) = cumabs (i, j) / (abs (mdt) * tmp (i))
            enddo
  
            do k = 1, km
                kr = km - k + 1
                dqv = qv (is:ie, k) - q (is:ie, j, kr, sphum)
                dql = ql (is:ie, k) - q (is:ie, j, kr, liq_wat)
                ps_dt = 1 + dqv + dql
                q (is:ie, j, kr, sphum) = qv (is:ie, k) / ps_dt
                q (is:ie, j, kr, liq_wat) = ql (is:ie, k) / ps_dt
                delp (is:ie, j, kr) = delp (is:ie, j, kr) * ps_dt
                q_liq = q (is:ie, j, kr, liq_wat) + q (is:ie, j, kr, rainwat)
                q_sol = q (is:ie, j, kr, ice_wat) + q (is:ie, j, kr, snowwat) + q (is:ie, j, kr, graupel)
                c_moist = (1 - (q (is:ie, j, kr, sphum) + q_liq + q_sol)) * cv_air + &
                    q (is:ie, j, kr, sphum) * cv_vap + q_liq * c_liq + q_sol * c_ice
                ps_dt = pt (is:ie, j, kr)
                pt (is:ie, j, kr) = pt (is:ie, j, kr) + (ta (is:ie, k) - pt (is:ie, j, kr)) * cp_air / c_moist
                ua (is:ie, j, kr) = uu (is:ie, k)
                va (is:ie, j, kr) = vv (is:ie, k)
            enddo
 
            ! compute wind tendency at A grid fori D grid wind update
            u_dt (is:ie, j, 1:km) = (ua (is:ie, j, 1:km) - u_dt (is:ie, j, 1:km)) / abs (mdt)
            v_dt (is:ie, j, 1:km) = (va (is:ie, j, 1:km) - v_dt (is:ie, j, 1:km)) / abs (mdt)

            ! update pe, peln, pk, ps
            do k = 2, km + 1
                pe (is:ie, k, j) = pe (is:ie, k-1, j) + delp (is:ie, j, k-1)
                peln (is:ie, k, j) = log (pe (is:ie, k, j))
                pk (is:ie, j, k) = exp (akap * peln (is:ie, k, j))
            enddo

            ps (is:ie, j) = pe (is:ie, km+1, j)

            ! update pkz
            if (.not. hydrostatic) then
#ifdef MOIST_CAPPA
                pkz (is:ie, j, 1:km) = exp (cappa (is:ie, j, 1:km) * &
                    log (rrg * delp (is:ie, j, 1:km) / &
                    delz (is:ie, j, 1:km) * pt (is:ie, j, 1:km)))
#else
                pkz (is:ie, j, 1:km) = exp (akap * log (rrg * delp (is:ie, j, 1:km) / &
                    delz (is:ie, j, 1:km) * pt (is:ie, j, 1:km)))
#endif
            endif
 
            if (consv .gt. consv_min) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                cvm = (1 - (q (is:ie, j, 1:km, sphum) + qliq + qsol)) * cv_air + &
                    q (is:ie, j, 1:km, sphum) * cv_vap + qliq * c_liq + qsol * c_ice
                te (is:ie, j, 1:km) = te (is:ie, j, 1:km) + &
                    cvm * pt (is:ie, j, 1:km) / ((1. + r_vir * q (is:ie, j, 1:km, sphum)) * &
                    (1. - (qliq + qsol))) * delp (is:ie, j, 1:km)
                do k = 1, km
                    te0_2d (is:ie, j) = te0_2d (is:ie, j) + te (is:ie, j, k)
                enddo
            endif

        enddo

        deallocate (lsm)

        deallocate (rn)
        deallocate (tmp)

        deallocate (dz)
        deallocate (zm)
        deallocate (dp)
        deallocate (pm)

        deallocate (ta)
        deallocate (qv)
        deallocate (ql)
        deallocate (uu)
        deallocate (vv)
        deallocate (ww)

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

        deallocate (u_dt)
        deallocate (v_dt)

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

        if (consv .gt. consv_min) then
            deallocate (u0)
            deallocate (v0)
            deallocate (dp0)
        endif

        call timing_off ('sa_sas')

    endif

    !-----------------------------------------------------------------------
    ! <<< Inline SA-SAS
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Inline SA-GWD >>>
    !-----------------------------------------------------------------------

    if ((.not. do_adiabatic_init) .and. do_inline_gwd) then

        call timing_on ('sa_gwd')

        allocate (dz (is:ie, 1:km))
        allocate (zm (is:ie, 1:km))
        allocate (zi (is:ie, 1:km+1))
        allocate (dp (is:ie, 1:km))
        allocate (pm (is:ie, 1:km))
        allocate (pi (is:ie, 1:km+1))
        allocate (pmk (is:ie, 1:km))

        allocate (ta (is:ie, 1:km))
        allocate (qv (is:ie, 1:km))
        allocate (uu (is:ie, 1:km))
        allocate (vv (is:ie, 1:km))

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

!$OMP parallel do default (none) shared (is, ie, js, je, isd, jsd, km, pe, ua, va, &
!$OMP                                    te, delp, hydrostatic, pt, peln, delz, &
!$OMP                                    rainwat, liq_wat, ice_wat, snowwat, graupel, &
!$OMP                                    sphum, pk, pkz, consv, te0_2d, gridstruct, q, &
!$OMP                                    mdt, cappa, rrg, akap, r_vir, ps, inline_gwd, &
!$OMP                                    kbot, ktop, kcnv, ptop, cumabs, inline_edmf, &
!$OMP                                    u_dt, v_dt) &
!$OMP                           private (gsize, dz, pi, pmk, zi, q_liq, q_sol, &
!$OMP                                    zm, dp, pm, qv, ta, uu, vv, qliq, qsol, &
!$OMP                                    cvm, kr, dqv, ps_dt, c_moist)

        do j = js, je
 
            gsize (is:ie) = sqrt (gridstruct%area_64 (is:ie, j))

            ! save ua, va for wind tendency calculation
            u_dt (is:ie, j, 1:km) = ua (is:ie, j, 1:km)
            v_dt (is:ie, j, 1:km) = va (is:ie, j, 1:km)

            if (consv .gt. consv_min) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                cvm = (1 - (q (is:ie, j, 1:km, sphum) + qliq + qsol)) * cv_air + &
                    q (is:ie, j, 1:km, sphum) * cv_vap + qliq * c_liq + qsol * c_ice
                te (is:ie, j, 1:km) = - cvm * pt (is:ie, j, 1:km) / ((1. + r_vir * q (is:ie, j, 1:km, sphum)) * &
                    (1. - (qliq + qsol))) * delp (is:ie, j, 1:km)
            endif

            zi (is:ie, 1) = 0.0
            pi (is:ie, km+1) = ptop
            do k = 1, km
                kr = km - k + 1
                dp (is:ie, k) = delp (is:ie, j, kr)
                pi (is:ie, kr) = pi (is:ie, kr+1) + delp (is:ie, j, k)
                if (.not. hydrostatic) then
                    pm (is:ie, k) = - dp (is:ie, k) / delz (is:ie, j, kr) * &
                        rdgas * pt (is:ie, j, kr) / grav
                    dz (is:ie, k) = delz (is:ie, j, kr)
                else
                    pm (is:ie, k) = dp (is:ie, k) / (peln (is:ie, kr+1, j) - peln (is:ie, kr, j))
                    dz (is:ie, k) = (peln (is:ie, kr, j) - peln (is:ie, kr+1, j)) * &
                        rdgas * pt (is:ie, j, kr) / grav
                endif
                pmk (is:ie, k) = exp (kappa * log (pm (is:ie, k) * 1.e-5))
                zi (is:ie, k+1) = zi (is:ie, k) - dz (is:ie, k) * grav
                if (k .eq. 1) then
                    zm (is:ie, k) = - 0.5 * dz (is:ie, k) * grav
                else
                    zm (is:ie, k) = zm (is:ie, k-1) - 0.5 * (dz (is:ie, k-1) + dz (is:ie, k)) * grav
                endif
                qv (is:ie, k) = q (is:ie, j, kr, sphum)
                ta (is:ie, k) = pt (is:ie, j, kr)
                uu (is:ie, k) = ua (is:ie, j, kr)
                vv (is:ie, k) = va (is:ie, j, kr)
            enddo

            call sa_gwd_oro (ie-is+1, km, uu, vv, ta, qv, abs (mdt), gsize, &
                inline_edmf%kpbl (is:ie, j), pi, dp, pm, pmk, zi, zm, &
                inline_gwd%hprime (is:ie, j), inline_gwd%oc (is:ie, j), inline_gwd%oa (is:ie, j, :), &
                inline_gwd%ol (is:ie, j, :), inline_gwd%theta (is:ie, j), inline_gwd%sigma (is:ie, j), &
                inline_gwd%gamma (is:ie, j), inline_gwd%elvmax (is:ie, j))

            call sa_gwd_cnv (ie-is+1, km, uu, vv, ta, qv, abs (mdt), gsize, cumabs (is:ie, j), &
                pm, pi, dp, ktop (is:ie, j), kbot (is:ie, j), kcnv (is:ie, j))
  
            do k = 1, km
                kr = km - k + 1
                dqv = qv (is:ie, k) - q (is:ie, j, kr, sphum)
                ps_dt = 1 + dqv
                q (is:ie, j, kr, sphum) = qv (is:ie, k) / ps_dt
                delp (is:ie, j, kr) = delp (is:ie, j, kr) * ps_dt
                q_liq = q (is:ie, j, kr, liq_wat) + q (is:ie, j, kr, rainwat)
                q_sol = q (is:ie, j, kr, ice_wat) + q (is:ie, j, kr, snowwat) + q (is:ie, j, kr, graupel)
                c_moist = (1 - (q (is:ie, j, kr, sphum) + q_liq + q_sol)) * cv_air + &
                    q (is:ie, j, kr, sphum) * cv_vap + q_liq * c_liq + q_sol * c_ice
                ps_dt = pt (is:ie, j, kr)
                pt (is:ie, j, kr) = pt (is:ie, j, kr) + (ta (is:ie, k) - pt (is:ie, j, kr)) * cp_air / c_moist
                ua (is:ie, j, kr) = uu (is:ie, k)
                va (is:ie, j, kr) = vv (is:ie, k)
            enddo
 
            ! compute wind tendency at A grid fori D grid wind update
            u_dt (is:ie, j, 1:km) = (ua (is:ie, j, 1:km) - u_dt (is:ie, j, 1:km)) / abs (mdt)
            v_dt (is:ie, j, 1:km) = (va (is:ie, j, 1:km) - v_dt (is:ie, j, 1:km)) / abs (mdt)

            ! update pe, peln, pk, ps
            do k = 2, km + 1
                pe (is:ie, k, j) = pe (is:ie, k-1, j) + delp (is:ie, j, k-1)
                peln (is:ie, k, j) = log (pe (is:ie, k, j))
                pk (is:ie, j, k) = exp (akap * peln (is:ie, k, j))
            enddo

            ps (is:ie, j) = pe (is:ie, km+1, j)

            ! update pkz
            if (.not. hydrostatic) then
#ifdef MOIST_CAPPA
                pkz (is:ie, j, 1:km) = exp (cappa (is:ie, j, 1:km) * &
                    log (rrg * delp (is:ie, j, 1:km) / &
                    delz (is:ie, j, 1:km) * pt (is:ie, j, 1:km)))
#else
                pkz (is:ie, j, 1:km) = exp (akap * log (rrg * delp (is:ie, j, 1:km) / &
                    delz (is:ie, j, 1:km) * pt (is:ie, j, 1:km)))
#endif
            endif
 
            if (consv .gt. consv_min) then
                qliq = q (is:ie, j, 1:km, liq_wat) + q (is:ie, j, 1:km, rainwat)
                qsol = q (is:ie, j, 1:km, ice_wat) + q (is:ie, j, 1:km, snowwat) + q (is:ie, j, 1:km, graupel)
                cvm = (1 - (q (is:ie, j, 1:km, sphum) + qliq + qsol)) * cv_air + &
                    q (is:ie, j, 1:km, sphum) * cv_vap + qliq * c_liq + qsol * c_ice
                te (is:ie, j, 1:km) = te (is:ie, j, 1:km) + &
                    cvm * pt (is:ie, j, 1:km) / ((1. + r_vir * q (is:ie, j, 1:km, sphum)) * &
                    (1. - (qliq + qsol))) * delp (is:ie, j, 1:km)
                do k = 1, km
                    te0_2d (is:ie, j) = te0_2d (is:ie, j) + te (is:ie, j, k)
                enddo
            endif

        enddo

        deallocate (dz)
        deallocate (zm)
        deallocate (dp)
        deallocate (pm)
        deallocate (pi)
        deallocate (pmk)

        deallocate (ta)
        deallocate (qv)
        deallocate (uu)
        deallocate (vv)

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

        deallocate (u_dt)
        deallocate (v_dt)

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

        if (consv .gt. consv_min) then
            deallocate (u0)
            deallocate (v0)
            deallocate (dp0)
        endif

        call timing_off ('sa_gwd')

    endif

    !-----------------------------------------------------------------------
    ! <<< Inline SA-GWD
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
!$OMP                                    ccn_cm3, cin_cm3, inline_mp, do_inline_mp, ps, &
!$OMP                                    u_dt, v_dt, aerosol) &
!$OMP                           private (q2, q3, gsize, dz, wa)

        do j = js, je

            gsize (is:ie) = sqrt (gridstruct%area_64 (is:ie, j))

            if (aerosol .gt. 0) then
                q2 (is:ie, kmp:km) = q (is:ie, j, kmp:km, aerosol)
            elseif (ccn_cm3 .gt. 0) then
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
            ! note: the unit of prew, prer, prei, pres, preg is mm/day
            ! note: the unit of prefluxw, prefluxr, prefluxi, prefluxs, prefluxg is mm/day
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
                dz (is:ie, kmp:km) = (peln (is:ie, kmp:km, j) - peln (is:ie, kmp+1:km+1, j)) * &
                    rdgas * pt (is:ie, j, kmp:km) / grav
            endif

            call gfdl_mp_driver (q (is:ie, j, kmp:km, sphum), q (is:ie, j, kmp:km, liq_wat), &
                     q (is:ie, j, kmp:km, rainwat), q (is:ie, j, kmp:km, ice_wat), &
                     q (is:ie, j, kmp:km, snowwat), q (is:ie, j, kmp:km, graupel), &
                     q (is:ie, j, kmp:km, cld_amt), q2 (is:ie, kmp:km), &
                     q3 (is:ie, kmp:km), pt (is:ie, j, kmp:km), wa (is:ie, kmp:km), &
                     ua (is:ie, j, kmp:km), va (is:ie, j, kmp:km), dz (is:ie, kmp:km), &
                     delp (is:ie, j, kmp:km), gsize, abs (mdt), hs (is:ie, j), &
                     inline_mp%prew (is:ie, j), inline_mp%prer (is:ie, j), &
                     inline_mp%prei (is:ie, j), inline_mp%pres (is:ie, j), &
                     inline_mp%preg (is:ie, j), hydrostatic, is, ie, kmp, km, &
#ifdef USE_COND
                     q_con (is:ie, j, kmp:km), &
#else
                     q_con (isd:, jsd, 1:), &
#endif
#ifdef MOIST_CAPPA
                     cappa (is:ie, j, kmp:km), &
#else
                     cappa (isd:, jsd, 1:), &
#endif
                     consv .gt. consv_min, te (is:ie, j, kmp:km), &
                     inline_mp%pcw (is:ie, j, kmp:km), inline_mp%edw (is:ie, j, kmp:km), &
                     inline_mp%oew (is:ie, j, kmp:km), &
                     inline_mp%rrw (is:ie, j, kmp:km), inline_mp%tvw (is:ie, j, kmp:km), &
                     inline_mp%pci (is:ie, j, kmp:km), inline_mp%edi (is:ie, j, kmp:km), &
                     inline_mp%oei (is:ie, j, kmp:km), &
                     inline_mp%rri (is:ie, j, kmp:km), inline_mp%tvi (is:ie, j, kmp:km), &
                     inline_mp%pcr (is:ie, j, kmp:km), inline_mp%edr (is:ie, j, kmp:km), &
                     inline_mp%oer (is:ie, j, kmp:km), &
                     inline_mp%rrr (is:ie, j, kmp:km), inline_mp%tvr (is:ie, j, kmp:km), &
                     inline_mp%pcs (is:ie, j, kmp:km), inline_mp%eds (is:ie, j, kmp:km), &
                     inline_mp%oes (is:ie, j, kmp:km), &
                     inline_mp%rrs (is:ie, j, kmp:km), inline_mp%tvs (is:ie, j, kmp:km), &
                     inline_mp%pcg (is:ie, j, kmp:km), inline_mp%edg (is:ie, j, kmp:km), &
                     inline_mp%oeg (is:ie, j, kmp:km), &
                     inline_mp%rrg (is:ie, j, kmp:km), inline_mp%tvg (is:ie, j, kmp:km), &
                     inline_mp%prefluxw(is:ie, j, kmp:km), &
                     inline_mp%prefluxr(is:ie, j, kmp:km), inline_mp%prefluxi(is:ie, j, kmp:km), &
                     inline_mp%prefluxs(is:ie, j, kmp:km), inline_mp%prefluxg(is:ie, j, kmp:km), &
                     inline_mp%cond (is:ie, j), inline_mp%dep (is:ie, j), inline_mp%reevap (is:ie, j), &
                     inline_mp%sub (is:ie, j), last_step, do_inline_mp)

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

        deallocate (u_dt)
        deallocate (v_dt)

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
