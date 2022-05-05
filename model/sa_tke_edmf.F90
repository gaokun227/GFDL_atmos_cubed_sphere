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
! Scale-Aware Turbulent-Kinetic-Energy based Moist-Eddy-Diffusivity-Mass-Flux
! (SA-TKE-EDMF) Subgrid Vertical Turbulence Mixing Scheme
! For the convective boundary layer, the scheme adopts EDMF parameterization
! (Siebesma et al. 2007) to take into account non-local transport by
! large eddies (mfpblt.f).
! A new mass-flux parameterizaiton for stratocumulus-top-induced turbulence
! mixing has been introduced (previously, it was eddy diffusion form) (mfscu.f).
! For local turbulence mixing, a TKE closure model is used.
! Developers: Jongil Han, Kun Gao, Linjiong Zhou, and the GFDL FV3 Team
! References: Han et al. (2016), Han and Bretherton (2019)
! =======================================================================

! =======================================================================
! Updates at GFDL:
! 1) Jul 2019 by Kun Gao
!    goal: to allow for tke advection
!    change: rearange tracers (q1g)
!    TKE no longer needs to be the last tracer
! 2) Nov 2019 by Kun Gao
!    turn off non-local mixing for hydrometers to avoid unphysical negative values
! 3) Jun 2020 by Kun Gao
!    a) add option for turning off upper-limter on background diff. in inversion layer
!       over land/ice points (cap_k0_land)
!    b) use different xkzm_m, xkzm_h for land, ocean and sea ice points
!    c) add option for turning off hb19 formula for surface backgroud diff. (do_dk_hb19)
! 4) May 2022 by Linjiong Zhou
!    put it into the FV3 dynamical core and revise accordingly
! =======================================================================

module sa_tke_edmf_mod

    use gfdl_mp_mod, only: mqs

    implicit none

    private

    ! -----------------------------------------------------------------------
    ! public subroutines, functions, and variables
    ! -----------------------------------------------------------------------
    
    public :: sa_tke_edmf_init
    public :: sa_tke_edmf

    ! -----------------------------------------------------------------------
    ! physics constants
    ! -----------------------------------------------------------------------
    
    real, parameter :: grav = 9.80665 ! acceleration due to gravity (m/s^2), ref: IFS
    
    real, parameter :: rdgas = 287.05 ! gas constant for dry air (J/kg/K): ref: GFDL, GFS
    real, parameter :: rvgas = 461.50 ! gas constant for water vapor (J/kg/K): ref: GFDL, GFS

    real, parameter :: zvir = rvgas / rdgas - 1. ! 0.6077667316114637
    real, parameter :: eps = rdgas / rvgas ! 0.6219934994582882
    real, parameter :: epsm1 = rdgas / rvgas - 1. ! -0.3780065005417118

    real, parameter :: tice = 273.15 ! freezing temperature (K): ref: GFDL, GFS

    real, parameter :: cp_air = 1004.6 ! heat capacity of dry air at constant pressure (J/kg/K): ref: GFDL, GFS
    real, parameter :: cp_vap = 4.0 * rvgas ! 1846.0885419672554, heat capacity of water vapor at constnat pressure (J/kg/K)

    real, parameter :: c_liq = 4.218e3 ! heat capacity of water at 0 deg C (J/kg/K), ref: IFS

    real, parameter :: hlv = 2.5e6 ! latent heat of evaporation at 0 deg C (J/kg): ref: GFDL, GFS
    real, parameter :: hlf = 3.3358e5 ! latent heat of fusion at 0 deg C (J/kg): ref: GFDL, GFS
    
    ! -----------------------------------------------------------------------
    ! namelist parameters
    ! -----------------------------------------------------------------------

    logical :: cap_k0_land = .true.  ! flag for applying limter on background diff in inversion
    logical :: do_dk_hb19  = .false. ! flag for using hb19 formula for background diff
    logical :: dspheat     = .false. ! flag for tke dissipative heating
    
    real :: xkzm_mo  = 1.0   ! bkgd_vdif_m background vertical diffusion for momentum over ocean
    real :: xkzm_ho  = 1.0   ! bkgd_vdif_h background vertical diffusion for heat q over ocean
    real :: xkzm_ml  = 1.0   ! bkgd_vdif_m background vertical diffusion for momentum over land
    real :: xkzm_hl  = 1.0   ! bkgd_vdif_h background vertical diffusion for heat q over land
    real :: xkzm_mi  = 1.0   ! bkgd_vdif_m background vertical diffusion for momentum over ice
    real :: xkzm_hi  = 1.0   ! bkgd_vdif_h background vertical diffusion for heat q over ice
    real :: xkzm_s   = 1.0   ! bkgd_vdif_s sigma threshold for background mom. diffusion
    real :: xkzm_lim = 0.01  ! background vertical diffusion limit
    real :: xkzm_fac = 1.0   ! background vertical diffusion factor
    real :: xkzinv   = 0.15  ! diffusivity in inversion layers
    real :: xkgdx    = 25.e3 ! background vertical diffusion threshold
    real :: rlmn     = 30.   ! lower-limter on asymtotic mixing length in satmedmfdiff.f
    real :: rlmx     = 300.  ! upper-limter on asymtotic mixing length in satmedmfdiff.f

    ! -----------------------------------------------------------------------
    ! namelist
    ! -----------------------------------------------------------------------
    
    namelist / sa_tke_edmf_nml / &
        xkzm_mo, xkzm_ho, xkzm_ml, xkzm_hl, xkzm_mi, xkzm_hi, xkzm_s, &
        xkzm_lim, xkzm_fac, xkzinv, xkgdx, rlmn, rlmx, &
        cap_k0_land, do_dk_hb19, dspheat

contains

! =======================================================================
! SA-TKE-EDMF initialization
! =======================================================================

subroutine sa_tke_edmf_init (input_nml_file, logunit)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: logunit
    
    character (len = *), intent (in) :: input_nml_file (:)
    
    ! -----------------------------------------------------------------------
    ! read namelist
    ! -----------------------------------------------------------------------
    
    read (input_nml_file, nml = sa_tke_edmf_nml)
    
    ! -----------------------------------------------------------------------
    ! write namelist to log file
    ! -----------------------------------------------------------------------
    
    write (logunit, *) " ================================================================== "
    write (logunit, *) "sa_tke_edmf_mod"
    write (logunit, nml = sa_tke_edmf_nml)
    
end subroutine sa_tke_edmf_init

! =======================================================================
! SA-TKE-EDMF scheme
! =======================================================================

subroutine sa_tke_edmf (im, km, ntrac, ntcw, ntiw, ntke, &
        delt, u1, v1, t1, q1, gsize, islimsk, &
        swh, hlw, xmu, rbsoil, zorl, u10m, v10m, fm, fh, &
        tsea, heat, evap, stress, spd1, kinver, &
        psk, del, prsi, prsl, prslk, phii, phil, &
        hpbl, kpbl, dusfc, dvsfc, dtsfc, dqsfc, dkt_out)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: im, km, ntrac, ntcw, ntiw, ntke
    integer, intent (in) :: kinver (im), islimsk (im)

    real, intent (in) :: delt
    real, intent (in) :: swh (im, km), hlw (im, km), &
        xmu (im), gsize (im), &
        psk (im), rbsoil (im), &
        zorl (im), tsea (im), &
        u10m (im), v10m (im), &
        fm (im), fh (im), &
        evap (im), heat (im), &
        stress (im), spd1 (im), &
        prsi (im, km + 1), del (im, km), &
        prsl (im, km), prslk (im, km), &
        phii (im, km + 1), phil (im, km)

    real, intent (inout) :: u1 (im, km), v1 (im, km), &
        t1 (im, km), q1 (im, km, ntrac)

    integer, intent (out), optional :: kpbl (im)
    
    real, intent (out) :: hpbl (im)
    
    real, intent (out), optional :: dusfc (im), dvsfc (im), dtsfc (im), dqsfc (im), &
        dkt_out (im, km)
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i, is, k, kk, n, km1, kmpbl, kmscu, ntrac1, ntcw_new
    integer :: lcld (im), kcld (im), krad (im), mrad (im)
    integer :: kx1 (im), kpblx (im)
    
    real :: tke (im, km), tkeh (im, km - 1)
    
    real :: theta (im, km), thvx (im, km), thlvx (im, km), &
        qlx (im, km), thetae (im, km), thlx (im, km), &
        slx (im, km), svx (im, km), qtx (im, km), &
        tvx (im, km), pix (im, km), radx (im, km - 1), &
        dku (im, km - 1), dkt (im, km - 1), dkq (im, km - 1), &
        cku (im, km - 1), ckt (im, km - 1), q1g (im, km, ntrac), &
        vdt (im, km), udt (im, km), tdt (im, km), qdt (im, km)
    
    real :: plyr (im, km), rhly (im, km), cfly (im, km), &
        qstl (im, km)
    
    real :: dtdz1 (im), gdx (im), &
        phih (im), phim (im), prn (im, km - 1), &
        rbdn (im), rbup (im), thermal (im), &
        ustar (im), wstar (im), hpblx (im), &
        ust3 (im), wst3 (im), &
        z0 (im), crb (im), &
        hgamt (im), hgamq (im), &
        wscale (im), vpert (im), &
        zol (im), sflux (im), radj (im), &
        tx1 (im), tx2 (im)
    
    real :: radmin (im)
    
    real :: zi (im, km + 1), zl (im, km), zm (im, km), &
        xkzo (im, km - 1), xkzmo (im, km - 1), &
        xkzm_hx (im), xkzm_mx (im), &
        rdzt (im, km - 1), &
        al (im, km - 1), ad (im, km), au (im, km - 1), &
        f1 (im, km), f2 (im, km * (ntrac - 1))
    
    real :: elm (im, km), ele (im, km), rle (im, km - 1), &
        ckz (im, km), chz (im, km), &
        diss (im, km - 1), prod (im, km - 1), &
        bf (im, km - 1), shr2 (im, km - 1), &
        xlamue (im, km - 1), xlamde (im, km - 1), &
        gotvx (im, km), rlam (im, km - 1)
    
    ! variables for updrafts (thermals)
    real :: tcko (im, km), qcko (im, km, ntrac), &
        ucko (im, km), vcko (im, km), &
        buou (im, km), xmf (im, km)
    
    ! variables for stratocumulus - top induced downdrafts
    real :: tcdo (im, km), qcdo (im, km, ntrac), &
        ucdo (im, km), vcdo (im, km), &
        buod (im, km), xmfd (im, km)
    
    logical :: pblflg (im), sfcflg (im), flg (im)
    logical :: scuflg (im), pcnvflg (im)
    logical :: mlenflg
    
    ! pcnvflg: true for unstable pbl
    real :: aphi16, aphi5, &
        wfac, cfac, &
        gamcrt, gamcrq, sfcfrac, &
        conq, cont, conw, &
        dsdz2, dsdzt, dkmax, &
        dsig, dt2, dtodsd, &
        dtodsu, g, factor, dz, &
        gocp, gravi, zol1, zolcru, &
        buop, shrp, dtn, cdtn, &
        prnum, prmax, prmin, prtke, &
        prscu, dw2, dw2min, zk, &
        elmfac, elefac, dspmax, &
        alp, clwt, cql, &
        f0, robn, crbmin, crbmax, &
        es, qs, value, onemrh, &
        cfh, gamma, elocp, el2orc, &
        epsi, beta, chx, cqx, &
        rdt, rdz, qmin, qlmin, &
        ri, rimin, &
        rbcr, rbint, tdzmin, &
        elmx, &
        ttend, utend, qtend, &
        zfac, zfmin, vk, spdk2, &
        tkmin, dspfac, &
        zlup, zldn, bsum, &
        tem, tem1, tem2, &
        ptem, ptem0, ptem1, ptem2
    
    real :: ck0, ck1, ch0, ch1, ce0, rchck
    
    real :: qlcr, zstblmax
    
    real :: h1
    
    parameter (gravi = 1.0 / grav)
    parameter (g = grav)
    parameter (gocp = g / cp_air)
    parameter (cont = cp_air / g, conq = hlv / g, conw = 1.0 / g)
    parameter (elocp = hlv / cp_air, el2orc = hlv * hlv / (rvgas * cp_air))
    parameter (wfac = 7.0, cfac = 4.5)
    parameter (gamcrt = 3., gamcrq = 0., sfcfrac = 0.1)
    parameter (vk = 0.4, rimin = - 100.)
    parameter (rbcr = 0.25, zolcru = - 0.02, tdzmin = 1.e-3)
    parameter (prmin = 0.25, prmax = 4.0, prtke = 1.0, prscu = 0.67)
    parameter (f0 = 1.e-4, crbmin = 0.15, crbmax = 0.35)
    parameter (tkmin = 1.e-9, dspfac = 0.5, dspmax = 10.0)
    parameter (qmin = 1.e-8, qlmin = 1.e-12, zfmin = 1.e-8)
    parameter (aphi5 = 5., aphi16 = 16.)
    parameter (elmfac = 1.0, elefac = 1.0, cql = 100.)
    parameter (dw2min = 1.e-4, dkmax = 1000.)
    parameter (qlcr = 3.5e-5, zstblmax = 2500.)
    parameter (h1 = 0.33333333)
    parameter (ck0 = 0.4, ck1 = 0.15, ch0 = 0.4, ch1 = 0.15, ce0 = 0.4)
    parameter (rchck = 1.5, cdtn = 25.)
    
    elmx = rlmx

    ! -----------------------------------------------------------------------
    ! kgao note (jul 2019)
    ! the code was originally written assuming ntke = ntrac
    ! in this version ntke does not need to be equal to ntrac
    ! in the following we rearrange q1g so that tke is the last tracer
    ! -----------------------------------------------------------------------
    
    !if (ntrac >= 3) then
    if (ntke == ntrac) then ! tke is the last tracer
        q1g (:, :, :) = q1 (:, :, :)
    else ! tke is not
        do kk = 1, ntke - 1
            q1g (:, :, kk) = q1 (:, :, kk)
        enddo
        do kk = ntke + 1, ntrac
            q1g (:, :, kk - 1) = q1 (:, :, kk)
        enddo
        q1g (:, :, ntrac) = q1 (:, :, ntke)
    endif
    !endif
    
    dt2 = delt
    rdt = 1. / dt2
    
    ntrac1 = ntrac - 1
    km1 = km - 1
    kmpbl = km / 2
    kmscu = km / 2
    
    do k = 1, km
        do i = 1, im
            zi (i, k) = phii (i, k) * gravi
            zl (i, k) = phil (i, k) * gravi
            xmf (i, k) = 0.
            xmfd (i, k) = 0.
            buou (i, k) = 0.
            buod (i, k) = 0.
            ckz (i, k) = ck1
            chz (i, k) = ch1
        enddo
    enddo

    do i = 1, im
        zi (i, km + 1) = phii (i, km + 1) * gravi
    enddo
    do k = 1, km
        do i = 1, im
            zm (i, k) = zi (i, k + 1)
        enddo
    enddo

    ! -----------------------------------------------------------------------
    ! horizontal grid size
    ! -----------------------------------------------------------------------

    do i = 1, im
        gdx (i) = gsize (i)
    enddo
    
    do k = 1, km
        do i = 1, im
            tke (i, k) = max (q1g (i, k, ntrac), tkmin) ! tke at layer centers
        enddo
    enddo
    do k = 1, km1
        do i = 1, im
            tkeh (i, k) = 0.5 * (tke (i, k) + tke (i, k + 1)) ! tke at interfaces
        enddo
    enddo
    
    do k = 1, km1
        do i = 1, im
            rdzt (i, k) = 1.0 / (zl (i, k + 1) - zl (i, k))
            prn (i, k) = 1.0
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! han and bretherton, 2019
    ! set background diffusivities as a function of
    ! horizontal grid size with xkzm_h & xkzm_m for gdx >= xkgdx
    ! and 0.01 for gdx = 5m, i.e.,
    ! xkzm_hx = 0.01 + (xkzm_h - 0.01) / (xkgdx - 5.) * (gdx - 5.)
    ! xkzm_mx = 0.01 + (xkzm_h - 0.01) / (xkgdx - 5.) * (gdx - 5.)
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        kx1 (i) = 1
        tx1 (i) = 1.0 / prsi (i, 1)
        tx2 (i) = tx1 (i)
        
        ! -----------------------------------------------------------------------
        ! kgao change - set surface value of background diff (dk) below
        ! -----------------------------------------------------------------------

        if (do_dk_hb19) then ! use eq43 in hb2019
            
            if (gdx (i) >= xkgdx) then ! resolution coarser than xkgdx
                if (islimsk (i) == 1) then ! land points
                    xkzm_hx (i) = xkzm_hl
                    xkzm_mx (i) = xkzm_ml
                elseif (islimsk (i) == 2) then! sea ice points
                    xkzm_hx (i) = xkzm_hi
                    xkzm_mx (i) = xkzm_mi
                else ! ocean points
                    xkzm_hx (i) = xkzm_ho
                    xkzm_mx (i) = xkzm_mo
                endif
            else ! resolution finer than xkgdx
                tem = 1. / (xkgdx - 5.)
                if (islimsk (i) == 1) then ! land points
                    tem1 = (xkzm_hl - xkzm_lim) * tem
                    tem2 = (xkzm_ml - xkzm_lim) * tem
                elseif (islimsk (i) == 2) then! sea ice points
                    tem1 = (xkzm_hi - xkzm_lim) * tem
                    tem2 = (xkzm_mi - xkzm_lim) * tem
                else ! ocean points
                    tem1 = (xkzm_ho - xkzm_lim) * tem
                    tem2 = (xkzm_mo - xkzm_lim) * tem
                endif
                ptem = gdx (i) - 5.
                xkzm_hx (i) = xkzm_lim + tem1 * ptem
                xkzm_mx (i) = xkzm_lim + tem2 * ptem
            endif
            
        else ! use values in the namelist; no res dependency
            
            if (islimsk (i) == 1) then ! land points
                xkzm_hx (i) = xkzm_hl
                xkzm_mx (i) = xkzm_ml
            elseif (islimsk (i) == 2) then ! sea ice points
                xkzm_hx (i) = xkzm_hi
                xkzm_mx (i) = xkzm_mi
            else ! ocean points
                xkzm_hx (i) = xkzm_ho
                xkzm_mx (i) = xkzm_mo
            endif
            
        endif
    enddo
    
    do k = 1, km1
        do i = 1, im
            xkzo (i, k) = 0.0
            xkzmo (i, k) = 0.0
            if (k < kinver (i)) then
                ! -----------------------------------------------------------------------
                ! vertical background diffusivity
                ! -----------------------------------------------------------------------
                ptem = prsi (i, k + 1) * tx1 (i)
                tem1 = 1.0 - ptem
                tem1 = tem1 * tem1 * 10.0
                xkzo (i, k) = xkzm_hx (i) * min (1.0, exp (- tem1))
                ! -----------------------------------------------------------------------
                ! vertical background diffusivity for momentum
                ! -----------------------------------------------------------------------
                if (ptem >= xkzm_s) then
                    xkzmo (i, k) = xkzm_mx (i)
                    kx1 (i) = k + 1
                else
                    if (k == kx1 (i) .and. k > 1) tx2 (i) = 1.0 / prsi (i, k)
                    tem1 = 1.0 - prsi (i, k + 1) * tx2 (i)
                    tem1 = tem1 * tem1 * 5.0
                    xkzmo (i, k) = xkzm_mx (i) * min (1.0, exp (- tem1))
                endif
            endif
        enddo
    enddo
    
    do i = 1, im
        z0 (i) = 0.01 * zorl (i)
        if (present (dusfc)) dusfc (i) = 0.
        if (present (dvsfc)) dvsfc (i) = 0.
        if (present (dtsfc)) dtsfc (i) = 0.
        if (present (dqsfc)) dqsfc (i) = 0.
        kpbl (i) = 1
        hpbl (i) = 0.
        kpblx (i) = 1
        hpblx (i) = 0.
        pblflg (i) = .true.
        sfcflg (i) = .true.
        if (rbsoil (i) > 0.) sfcflg (i) = .false.
        pcnvflg (i) = .false.
        scuflg (i) = .true.
        if (scuflg (i)) then
            radmin (i) = 0.
            mrad (i) = km1
            krad (i) = 1
            lcld (i) = km1
            kcld (i) = km1
        endif
    enddo
    
    do k = 1, km
        do i = 1, im
            pix (i, k) = psk (i) / prslk (i, k)
            theta (i, k) = t1 (i, k) * pix (i, k)
            if (ntiw > 0) then
                tem = max (q1 (i, k, ntcw), qlmin)
                tem1 = max (q1 (i, k, ntiw), qlmin)
                qlx (i, k) = tem + tem1
                ptem = hlv * tem + (hlv + hlf) * tem1
                slx (i, k) = cp_air * t1 (i, k) + phil (i, k) - ptem
            else
                qlx (i, k) = max (q1 (i, k, ntcw), qlmin)
                slx (i, k) = cp_air * t1 (i, k) + phil (i, k) - hlv * qlx (i, k)
            endif
            tem2 = 1. + zvir * max (q1g (i, k, 1), qmin) - qlx (i, k)
            thvx (i, k) = theta (i, k) * tem2
            tvx (i, k) = t1 (i, k) * tem2
            qtx (i, k) = max (q1g (i, k, 1), qmin) + qlx (i, k)
            thlx (i, k) = theta (i, k) - pix (i, k) * elocp * qlx (i, k)
            thlvx (i, k) = thlx (i, k) * (1. + zvir * qtx (i, k))
            svx (i, k) = cp_air * tvx (i, k)
            ptem1 = elocp * pix (i, k) * max (q1g (i, k, 1), qmin)
            thetae (i, k) = theta (i, k) + ptem1
            gotvx (i, k) = g / tvx (i, k)
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! the background vertical diffusivities in the inversion layers are limited
    ! to be less than or equal to xkzminv
    ! -----------------------------------------------------------------------
    
    do k = 1, km1
        do i = 1, im
            tem1 = (tvx (i, k + 1) - tvx (i, k)) * rdzt (i, k)
            
            if (cap_k0_land) then
                if (tem1 > 1.e-5) then
                    xkzo (i, k) = min (xkzo (i, k), xkzinv)
                    xkzmo (i, k) = min (xkzmo (i, k), xkzinv)
                endif
            else
                ! -----------------------------------------------------------------------
                ! kgao note: do not apply upper - limiter over land and sea ice points
                ! (consistent with change in satmedmfdifq.f in jun 2020)
                ! -----------------------------------------------------------------------
                if (tem1 > 0. .and. islimsk (i) == 0) then
                    xkzo (i, k) = min (xkzo (i, k), xkzinv)
                    xkzmo (i, k) = min (xkzmo (i, k), xkzinv)
                endif
            endif
            
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute an empirical cloud fraction based on
    ! xu & randall's (1996, jas) study
    ! -----------------------------------------------------------------------
    
    do k = 1, km
        do i = 1, im
            plyr (i, k) = 0.01 * prsl (i, k) ! pa to mb (hpa)
            ! compute relative humidity
            es = 0.01 * mqs (t1 (i, k)) ! mqs in pa
            qs = max (qmin, eps * es / (plyr (i, k) + epsm1 * es))
            rhly (i, k) = max (0.0, min (1.0, max (qmin, q1g (i, k, 1)) / qs))
            qstl (i, k) = qs
        enddo
    enddo
    
    do k = 1, km
        do i = 1, im
            cfly (i, k) = 0.
            clwt = 1.0e-6 * (plyr (i, k) * 0.001)
            if (qlx (i, k) > clwt) then
                onemrh = max (1.e-10, 1.0 - rhly (i, k))
                tem1 = min (max ((onemrh * qstl (i, k)) ** 0.49, 0.0001), 1.0)
                tem1 = cql / tem1
                value = max (min (tem1 * qlx (i, k), 50.0), 0.0)
                tem2 = sqrt (sqrt (rhly (i, k)))
                cfly (i, k) = min (max (tem2 * (1.0 - exp (- value)), 0.0), 1.0)
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute buoyancy modified by clouds
    ! -----------------------------------------------------------------------
    
    do k = 1, km1
        do i = 1, im
            tem = 0.5 * (svx (i, k) + svx (i, k + 1))
            tem1 = 0.5 * (t1 (i, k) + t1 (i, k + 1))
            tem2 = 0.5 * (qstl (i, k) + qstl (i, k + 1))
            cfh = min (cfly (i, k + 1), 0.5 * (cfly (i, k) + cfly (i, k + 1)))
            alp = g / tem
            gamma = el2orc * tem2 / (tem1 ** 2)
            epsi = tem1 / elocp
            beta = (1. + gamma * epsi * (1. + zvir)) / (1. + gamma)
            chx = cfh * alp * beta + (1. - cfh) * alp
            cqx = cfh * alp * hlv * (beta - epsi)
            cqx = cqx + (1. - cfh) * zvir * g
            ptem1 = (slx (i, k + 1) - slx (i, k)) * rdzt (i, k)
            ptem2 = (qtx (i, k + 1) - qtx (i, k)) * rdzt (i, k)
            bf (i, k) = chx * ptem1 + cqx * ptem2
        enddo
    enddo

    do k = 1, km1
        do i = 1, im
            dku (i, k) = 0.
            dkt (i, k) = 0.
            dkq (i, k) = 0.
            cku (i, k) = 0.
            ckt (i, k) = 0.
            tem = zi (i, k + 1) - zi (i, k)
            radx (i, k) = tem * (swh (i, k) * xmu (i) + hlw (i, k))
        enddo
    enddo
    
    do i = 1, im
        sflux (i) = heat (i) + evap (i) * zvir * theta (i, 1)
        if (.not.sfcflg (i) .or. sflux (i) <= 0.) pblflg (i) = .false.
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute critical bulk richardson number
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        if (pblflg (i)) then
            ! thermal (i) = thvx (i, 1)
            thermal (i) = thlvx (i, 1)
            crb (i) = rbcr
        else
            thermal (i) = tsea (i) * (1. + zvir * max (q1g (i, 1, 1), qmin))
            tem = sqrt (u10m (i) ** 2 + v10m (i) ** 2)
            tem = max (tem, 1.)
            robn = tem / (f0 * z0 (i))
            tem1 = 1.e-7 * robn
            crb (i) = 0.16 * (tem1 ** (- 0.18))
            crb (i) = max (min (crb (i), crbmax), crbmin)
        endif
    enddo
    
    do i = 1, im
        dtdz1 (i) = dt2 / (zi (i, 2) - zi (i, 1))
    enddo
    
    do i = 1, im
        ustar (i) = sqrt (stress (i))
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute buoyancy (bf) and winshear square
    ! -----------------------------------------------------------------------
    
    do k = 1, km1
        do i = 1, im
            rdz = rdzt (i, k)
            ! bf (i, k) = gotvx (i, k) * (thvx (i, k + 1) - thvx (i, k)) * rdz
            dw2 = (u1 (i, k) - u1 (i, k + 1)) ** 2 &
                 + (v1 (i, k) - v1 (i, k + 1)) ** 2
            shr2 (i, k) = max (dw2, dw2min) * rdz * rdz
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! find pbl height based on bulk richardson number (mrf pbl scheme)
    ! and also for diagnostic purpose
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        flg (i) = .false.
        rbup (i) = rbsoil (i)
    enddo
    
    do k = 1, kmpbl
        do i = 1, im
            if (.not.flg (i)) then
                rbdn (i) = rbup (i)
                spdk2 = max ((u1 (i, k) ** 2 + v1 (i, k) ** 2), 1.)
                ! rbup (i) = (thvx (i, k) - thermal (i)) * &
                ! (g * zl (i, k) / thvx (i, 1)) / spdk2
                rbup (i) = (thlvx (i, k) - thermal (i)) * &
                     (g * zl (i, k) / thlvx (i, 1)) / spdk2
                kpblx (i) = k
                flg (i) = rbup (i) > crb (i)
            endif
        enddo
    enddo

    do i = 1, im
        if (kpblx (i) > 1) then
            k = kpblx (i)
            if (rbdn (i) >= crb (i)) then
                rbint = 0.
            elseif (rbup (i) <= crb (i)) then
                rbint = 1.
            else
                rbint = (crb (i) - rbdn (i)) / (rbup (i) - rbdn (i))
            endif
            hpblx (i) = zl (i, k - 1) + rbint * (zl (i, k) - zl (i, k - 1))
            if (hpblx (i) < zi (i, kpblx (i))) kpblx (i) = kpblx (i) - 1
        else
            hpblx (i) = zl (i, 1)
            kpblx (i) = 1
        endif
        hpbl (i) = hpblx (i)
        kpbl (i) = kpblx (i)
        if (kpbl (i) <= 1) pblflg (i) = .false.
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute similarity parameters
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        zol (i) = max (rbsoil (i) * fm (i) * fm (i) / fh (i), rimin)
        if (sfcflg (i)) then
            zol (i) = min (zol (i), - zfmin)
        else
            zol (i) = max (zol (i), zfmin)
        endif
        
        zol1 = zol (i) * sfcfrac * hpbl (i) / zl (i, 1)
        if (sfcflg (i)) then
            tem = 1.0 / (1. - aphi16 * zol1)
            phih (i) = sqrt (tem)
            phim (i) = sqrt (phih (i))
        else
            phim (i) = 1. + aphi5 * zol1
            phih (i) = phim (i)
        endif
    enddo
    
    do i = 1, im
        if (pblflg (i)) then
            if (zol (i) < zolcru) then
                pcnvflg (i) = .true.
            endif
            wst3 (i) = gotvx (i, 1) * sflux (i) * hpbl (i)
            wstar (i) = wst3 (i) ** h1
            ust3 (i) = ustar (i) ** 3.
            wscale (i) = (ust3 (i) + wfac * vk * wst3 (i) * sfcfrac) ** h1
            ptem = ustar (i) / aphi5
            wscale (i) = max (wscale (i), ptem)
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute a thermal excess
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        if (pcnvflg (i)) then
            hgamt (i) = heat (i) / wscale (i)
            hgamq (i) = evap (i) / wscale (i)
            vpert (i) = hgamt (i) + hgamq (i) * zvir * theta (i, 1)
            vpert (i) = max (vpert (i), 0.)
            tem = min (cfac * vpert (i), gamcrt)
            thermal (i) = thermal (i) + tem
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! enhance the pbl height by considering the thermal excess
    ! (overshoot pbl top)
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        flg (i) = .true.
        if (pcnvflg (i)) then
            flg (i) = .false.
            rbup (i) = rbsoil (i)
        endif
    enddo
    do k = 2, kmpbl
        do i = 1, im
            if (.not.flg (i)) then
                rbdn (i) = rbup (i)
                spdk2 = max ((u1 (i, k) ** 2 + v1 (i, k) ** 2), 1.)
                rbup (i) = (thlvx (i, k) - thermal (i)) * &
                     (g * zl (i, k) / thlvx (i, 1)) / spdk2
                kpbl (i) = k
                flg (i) = rbup (i) > crb (i)
            endif
        enddo
    enddo

    do i = 1, im
        if (pcnvflg (i)) then
            k = kpbl (i)
            if (rbdn (i) >= crb (i)) then
                rbint = 0.
            elseif (rbup (i) <= crb (i)) then
                rbint = 1.
            else
                rbint = (crb (i) - rbdn (i)) / (rbup (i) - rbdn (i))
            endif
            hpbl (i) = zl (i, k - 1) + rbint * (zl (i, k) - zl (i, k - 1))
            if (hpbl (i) < zi (i, kpbl (i))) then
                kpbl (i) = kpbl (i) - 1
            endif
            if (kpbl (i) <= 1) then
                pcnvflg (i) = .false.
                pblflg (i) = .false.
            endif
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! look for stratocumulus
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        flg (i) = scuflg (i)
    enddo
    do k = 1, km1
        do i = 1, im
            if (flg (i) .and.zl (i, k) >= zstblmax) then
                lcld (i) = k
                flg (i) = .false.
            endif
        enddo
    enddo
    do i = 1, im
        flg (i) = scuflg (i)
    enddo
    do k = kmscu, 1, - 1
        do i = 1, im
            if (flg (i) .and. k <= lcld (i)) then
                if (qlx (i, k) >= qlcr) then
                    kcld (i) = k
                    flg (i) = .false.
                endif
            endif
        enddo
    enddo
    do i = 1, im
        if (scuflg (i) .and. kcld (i) == km1) scuflg (i) = .false.
    enddo
    
    do i = 1, im
        flg (i) = scuflg (i)
    enddo
    do k = kmscu, 1, - 1
        do i = 1, im
            if (flg (i) .and. k <= kcld (i)) then
                if (qlx (i, k) >= qlcr) then
                    if (radx (i, k) < radmin (i)) then
                        radmin (i) = radx (i, k)
                        krad (i) = k
                    endif
                else
                    flg (i) = .false.
                endif
            endif
        enddo
    enddo
    do i = 1, im
        if (scuflg (i) .and. krad (i) <= 1) scuflg (i) = .false.
        if (scuflg (i) .and. radmin (i) >= 0.) scuflg (i) = .false.
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute components for mass flux mixing by large thermals
    ! -----------------------------------------------------------------------
    
    do k = 1, km
        do i = 1, im
            if (pcnvflg (i)) then
                tcko (i, k) = t1 (i, k)
                ucko (i, k) = u1 (i, k)
                vcko (i, k) = v1 (i, k)
            endif
            if (scuflg (i)) then
                tcdo (i, k) = t1 (i, k)
                ucdo (i, k) = u1 (i, k)
                vcdo (i, k) = v1 (i, k)
            endif
        enddo
    enddo
    do kk = 1, ntrac1
        do k = 1, km
            do i = 1, im
                if (pcnvflg (i)) then
                    qcko (i, k, kk) = q1g (i, k, kk)
                endif
                if (scuflg (i)) then
                    qcdo (i, k, kk) = q1g (i, k, kk)
                endif
            enddo
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! kgao note - change ntcw if q1g is rearranged
    ! -----------------------------------------------------------------------

    if (ntke > ntcw) then
        ntcw_new = ntcw
    else
        ntcw_new = ntcw - 1
    endif

    ! -----------------------------------------------------------------------
    ! edmf parameterization siebesma et al. (2007)
    ! -----------------------------------------------------------------------

    call mfpblt (im, km, kmpbl, ntcw_new, ntrac1, dt2, &
        pcnvflg, zl, zm, q1, t1, u1, v1, plyr, pix, thlx, thvx, &
        gdx, hpbl, kpbl, vpert, buou, xmf, &
        tcko, qcko, ucko, vcko, xlamue)

    ! -----------------------------------------------------------------------
    ! mass - flux parameterization for stratocumulus - top - induced turbulence mixing
    ! -----------------------------------------------------------------------

    call mfscu (im, km, kmscu, ntcw_new, ntrac1, dt2, &
        scuflg, zl, zm, q1, t1, u1, v1, plyr, pix, &
        thlx, thvx, thlvx, gdx, thetae, radj, &
        krad, mrad, radmin, buod, xmfd, &
        tcdo, qcdo, ucdo, vcdo, xlamde)
    
    ! -----------------------------------------------------------------------
    ! compute prandtl number and exchange coefficient varying with height
    ! -----------------------------------------------------------------------
    
    do k = 1, kmpbl
        do i = 1, im
            if (k < kpbl (i)) then
                tem = phih (i) / phim (i)
                ptem = - 3. * (max (zi (i, k + 1) - sfcfrac * hpbl (i), 0.)) ** 2. &
                     / hpbl (i) ** 2.
                if (pcnvflg (i)) then
                    prn (i, k) = 1. + (tem - 1.) * exp (ptem)
                else
                    prn (i, k) = tem
                endif
                prn (i, k) = min (prn (i, k), prmax)
                prn (i, k) = max (prn (i, k), prmin)
                
                ckz (i, k) = ck1 + (ck0 - ck1) * exp (ptem)
                ckz (i, k) = min (ckz (i, k), ck0)
                ckz (i, k) = max (ckz (i, k), ck1)
                chz (i, k) = ch1 + (ch0 - ch1) * exp (ptem)
                chz (i, k) = min (chz (i, k), ch0)
                chz (i, k) = max (chz (i, k), ch1)
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute an asymtotic mixing length
    ! -----------------------------------------------------------------------
    
    do k = 1, km1
        do i = 1, im
            zlup = 0.0
            bsum = 0.0
            mlenflg = .true.
            do n = k, km1
                if (mlenflg) then
                    dz = zl (i, n + 1) - zl (i, n)
                    ptem = gotvx (i, n) * (thvx (i, n + 1) - thvx (i, k)) * dz
                    ! ptem = gotvx (i, n) * (thlvx (i, n + 1) - thlvx (i, k)) * dz
                    bsum = bsum + ptem
                    zlup = zlup + dz
                    if (bsum >= tke (i, k)) then
                        if (ptem >= 0.) then
                            tem2 = max (ptem, zfmin)
                        else
                            tem2 = min (ptem, - zfmin)
                        endif
                        ptem1 = (bsum - tke (i, k)) / tem2
                        zlup = zlup - ptem1 * dz
                        zlup = max (zlup, 0.)
                        mlenflg = .false.
                    endif
                endif
            enddo
            zldn = 0.0
            bsum = 0.0
            mlenflg = .true.
            do n = k, 1, - 1
                if (mlenflg) then
                    if (n == 1) then
                        dz = zl (i, 1)
                        tem1 = tsea (i) * (1. + zvir * max (q1g (i, 1, 1), qmin))
                    else
                        dz = zl (i, n) - zl (i, n - 1)
                        tem1 = thvx (i, n - 1)
                        ! tem1 = thlvx (i, n - 1)
                    endif
                    ptem = gotvx (i, n) * (thvx (i, k) - tem1) * dz
                    ! ptem = gotvx (i, n) * (thlvx (i, k) - tem1) * dz
                    bsum = bsum + ptem
                    zldn = zldn + dz
                    if (bsum >= tke (i, k)) then
                        if (ptem >= 0.) then
                            tem2 = max (ptem, zfmin)
                        else
                            tem2 = min (ptem, - zfmin)
                        endif
                        ptem1 = (bsum - tke (i, k)) / tem2
                        zldn = zldn - ptem1 * dz
                        zldn = max (zldn, 0.)
                        mlenflg = .false.
                    endif
                endif
            enddo
            
            tem = 0.5 * (zi (i, k + 1) - zi (i, k))
            tem1 = min (tem, rlmn)
            
            ptem2 = min (zlup, zldn)
            rlam (i, k) = elmfac * ptem2
            rlam (i, k) = max (rlam (i, k), tem1)
            rlam (i, k) = min (rlam (i, k), rlmx)
            
            ptem2 = sqrt (zlup * zldn)
            ele (i, k) = elefac * ptem2
            ele (i, k) = max (ele (i, k), tem1)
            ele (i, k) = min (ele (i, k), elmx)
            
        enddo
    enddo
    
    do k = 1, km1
        do i = 1, im
            tem = vk * zl (i, k)
            if (zol (i) < 0.) then
                ptem = 1. - 100. * zol (i)
                ptem1 = ptem ** 0.2
                zk = tem * ptem1
            elseif (zol (i) >= 1.) then
                zk = tem / 3.7
            else
                ptem = 1. + 2.7 * zol (i)
                zk = tem / ptem
            endif
            elm (i, k) = zk * rlam (i, k) / (rlam (i, k) + zk)
            
            dz = zi (i, k + 1) - zi (i, k)
            tem = max (gdx (i), dz)
            elm (i, k) = min (elm (i, k), tem)
            ele (i, k) = min (ele (i, k), tem)
            
        enddo
    enddo
    do i = 1, im
        elm (i, km) = elm (i, km1)
        ele (i, km) = ele (i, km1)
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute eddy diffusivities
    ! -----------------------------------------------------------------------
    
    do k = 1, km1
        do i = 1, im
            tem = 0.5 * (elm (i, k) + elm (i, k + 1))
            tem = tem * sqrt (tkeh (i, k))
            if (k < kpbl (i)) then
                if (pblflg (i)) then
                    dku (i, k) = ckz (i, k) * tem
                    dkt (i, k) = dku (i, k) / prn (i, k)
                else
                    dkt (i, k) = chz (i, k) * tem
                    dku (i, k) = dkt (i, k) * prn (i, k)
                endif
            else
                ri = max (bf (i, k) / shr2 (i, k), rimin)
                if (ri < 0.) then ! unstable regime
                    dku (i, k) = ck1 * tem
                    dkt (i, k) = rchck * dku (i, k)
                else ! stable regime
                    dkt (i, k) = ch1 * tem
                    prnum = 1.0 + 2.1 * ri
                    prnum = min (prnum, prmax)
                    dku (i, k) = dkt (i, k) * prnum
                endif
            endif
            
            if (scuflg (i)) then
                if (k >= mrad (i) .and. k < krad (i)) then
                    tem1 = ckz (i, k) * tem
                    ptem1 = tem1 / prscu
                    dku (i, k) = max (dku (i, k), tem1)
                    dkt (i, k) = max (dkt (i, k), ptem1)
                endif
            endif
            
            dkq (i, k) = prtke * dkt (i, k)
            
            dkt (i, k) = min (dkt (i, k), dkmax)
            dkt (i, k) = max (dkt (i, k), xkzo (i, k))
            dkq (i, k) = min (dkq (i, k), dkmax)
            dkq (i, k) = max (dkq (i, k), xkzo (i, k))
            dku (i, k) = min (dku (i, k), dkmax)
            dku (i, k) = max (dku (i, k), xkzmo (i, k))
            
        enddo
    enddo
    
    do i = 1, im
        if (scuflg (i)) then
            k = krad (i)
            tem = bf (i, k) / gotvx (i, k)
            tem1 = max (tem, tdzmin)
            ptem = radj (i) / tem1
            dkt (i, k) = dkt (i, k) + ptem
            dku (i, k) = dku (i, k) + ptem
            dkq (i, k) = dkq (i, k) + ptem
        endif
    enddo
    
    if (present (dkt_out)) then
        do k = 1, km1
            do i = 1, im
                dkt_out (i, k) = dkt (i, k)
            enddo
        enddo
    endif
    
    ! -----------------------------------------------------------------------
    ! compute buoyancy and shear productions of tke
    ! -----------------------------------------------------------------------
    
    do k = 1, km1
        do i = 1, im
            if (k == 1) then
                tem = - dkt (i, 1) * bf (i, 1)
                ! if (pcnvflg (i)) then
                ! ptem1 = xmf (i, 1) * buou (i, 1)
                ! else
                ptem1 = 0.
                ! endif
                if (scuflg (i) .and. mrad (i) == 1) then
                    ptem2 = xmfd (i, 1) * buod (i, 1)
                else
                    ptem2 = 0.
                endif
                tem = tem + ptem1 + ptem2
                buop = 0.5 * (gotvx (i, 1) * sflux (i) + tem)
                
                tem1 = dku (i, 1) * shr2 (i, 1)
                
                tem = (u1 (i, 2) - u1 (i, 1)) * rdzt (i, 1)
                ! if (pcnvflg (i)) then
                ! ptem = xmf (i, 1) * tem
                ! ptem1 = 0.5 * ptem * (u1 (i, 2) - ucko (i, 2))
                ! else
                ptem1 = 0.
                ! endif
                if (scuflg (i) .and. mrad (i) == 1) then
                    ptem = ucdo (i, 1) + ucdo (i, 2) - u1 (i, 1) - u1 (i, 2)
                    ptem = 0.5 * tem * xmfd (i, 1) * ptem
                else
                    ptem = 0.
                endif
                ptem1 = ptem1 + ptem
                
                tem = (v1 (i, 2) - v1 (i, 1)) * rdzt (i, 1)
                ! if (pcnvflg (i)) then
                ! ptem = xmf (i, 1) * tem
                ! ptem2 = 0.5 * ptem * (v1 (i, 2) - vcko (i, 2))
                ! else
                ptem2 = 0.
                ! endif
                if (scuflg (i) .and. mrad (i) == 1) then
                    ptem = vcdo (i, 1) + vcdo (i, 2) - v1 (i, 1) - v1 (i, 2)
                    ptem = 0.5 * tem * xmfd (i, 1) * ptem
                else
                    ptem = 0.
                endif
                ptem2 = ptem2 + ptem
                
                ! tem2 = stress (i) * spd1 (i) / zl (i, 1)
                tem2 = stress (i) * ustar (i) * phim (i) / (vk * zl (i, 1))
                shrp = 0.5 * (tem1 + ptem1 + ptem2 + tem2)
            else
                tem1 = - dkt (i, k - 1) * bf (i, k - 1)
                tem2 = - dkt (i, k) * bf (i, k)
                tem = 0.5 * (tem1 + tem2)
                if (pcnvflg (i) .and. k <= kpbl (i)) then
                    ptem = 0.5 * (xmf (i, k - 1) + xmf (i, k))
                    ptem1 = ptem * buou (i, k)
                else
                    ptem1 = 0.
                endif
                if (scuflg (i)) then
                    if (k >= mrad (i) .and. k < krad (i)) then
                        ptem0 = 0.5 * (xmfd (i, k - 1) + xmfd (i, k))
                        ptem2 = ptem0 * buod (i, k)
                    else
                        ptem2 = 0.
                    endif
                else
                    ptem2 = 0.
                endif
                buop = tem + ptem1 + ptem2
                
                tem1 = dku (i, k - 1) * shr2 (i, k - 1)
                tem2 = dku (i, k) * shr2 (i, k)
                tem = 0.5 * (tem1 + tem2)
                tem1 = (u1 (i, k + 1) - u1 (i, k)) * rdzt (i, k)
                tem2 = (u1 (i, k) - u1 (i, k - 1)) * rdzt (i, k - 1)
                if (pcnvflg (i) .and. k <= kpbl (i)) then
                    ptem = xmf (i, k) * tem1 + xmf (i, k - 1) * tem2
                    ptem1 = 0.5 * ptem * (u1 (i, k) - ucko (i, k))
                else
                    ptem1 = 0.
                endif
                if (scuflg (i)) then
                    if (k >= mrad (i) .and. k < krad (i)) then
                        ptem0 = xmfd (i, k) * tem1 + xmfd (i, k - 1) * tem2
                        ptem2 = 0.5 * ptem0 * (ucdo (i, k) - u1 (i, k))
                    else
                        ptem2 = 0.
                    endif
                else
                    ptem2 = 0.
                endif
                shrp = tem + ptem1 + ptem2
                tem1 = (v1 (i, k + 1) - v1 (i, k)) * rdzt (i, k)
                tem2 = (v1 (i, k) - v1 (i, k - 1)) * rdzt (i, k - 1)
                if (pcnvflg (i) .and. k <= kpbl (i)) then
                    ptem = xmf (i, k) * tem1 + xmf (i, k - 1) * tem2
                    ptem1 = 0.5 * ptem * (v1 (i, k) - vcko (i, k))
                else
                    ptem1 = 0.
                endif
                if (scuflg (i)) then
                    if (k >= mrad (i) .and. k < krad (i)) then
                        ptem0 = xmfd (i, k) * tem1 + xmfd (i, k - 1) * tem2
                        ptem2 = 0.5 * ptem0 * (vcdo (i, k) - v1 (i, k))
                    else
                        ptem2 = 0.
                    endif
                else
                    ptem2 = 0.
                endif
                shrp = shrp + ptem1 + ptem2
            endif
            prod (i, k) = buop + shrp
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! first predict tke due to tke production & dissipation (diss)
    ! -----------------------------------------------------------------------
    
    do k = 1, km1
        do i = 1, im
            rle (i, k) = ce0 / ele (i, k)
        enddo
    enddo
    kk = max (nint (dt2 / cdtn), 1)
    dtn = dt2 / float (kk)
    do n = 1, kk
        do k = 1, km1
            do i = 1, im
                tem = sqrt (tke (i, k))
                diss (i, k) = rle (i, k) * tke (i, k) * tem
                tem1 = prod (i, k) + tke (i, k) / dtn
                diss (i, k) = max (min (diss (i, k), tem1), 0.)
                tke (i, k) = tke (i, k) + dtn * (prod (i, k) - diss (i, k)) ! no diffusion yet
                tke (i, k) = max (tke (i, k), tkmin)
            enddo
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute updraft & downdraft properties for tke
    ! -----------------------------------------------------------------------
    
    do k = 1, km
        do i = 1, im
            if (pcnvflg (i)) then
                ! kgao change
                ! qcko (i, k, ntke) = tke (i, k)
                qcko (i, k, ntrac) = tke (i, k)
            endif
            if (scuflg (i)) then
                ! kgao change
                ! qcdo (i, k, ntke) = tke (i, k)
                qcdo (i, k, ntrac) = tke (i, k)
            endif
        enddo
    enddo
    do k = 2, kmpbl
        do i = 1, im
            if (pcnvflg (i) .and. k <= kpbl (i)) then
                dz = zl (i, k) - zl (i, k - 1)
                tem = 0.5 * xlamue (i, k - 1) * dz
                factor = 1. + tem
                ! kgao change
                ! qcko (i, k, ntke) = ((1. - tem) * qcko (i, k - 1, ntke) + tem * &
                ! (tke (i, k) + tke (i, k - 1))) / factor
                qcko (i, k, ntrac) = ((1. - tem) * qcko (i, k - 1, ntrac) + tem * &
                     (tke (i, k) + tke (i, k - 1))) / factor
                
            endif
        enddo
    enddo
    do k = kmscu, 1, - 1
        do i = 1, im
            if (scuflg (i) .and. k < krad (i)) then
                if (k >= mrad (i)) then
                    dz = zl (i, k + 1) - zl (i, k)
                    tem = 0.5 * xlamde (i, k) * dz
                    factor = 1. + tem
                    ! kgao change
                    ! qcdo (i, k, ntke) = ((1. - tem) * qcdo (i, k + 1, ntke) + tem * &
                    ! (tke (i, k) + tke (i, k + 1))) / factor
                    qcdo (i, k, ntrac) = ((1. - tem) * qcdo (i, k + 1, ntrac) + tem * &
                         (tke (i, k) + tke (i, k + 1))) / factor
                endif
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute tridiagonal matrix elements for turbulent kinetic energy
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        ad (i, 1) = 1.0
        f1 (i, 1) = tke (i, 1)
    enddo
    
    do k = 1, km1
        do i = 1, im
            dtodsd = dt2 / del (i, k)
            dtodsu = dt2 / del (i, k + 1)
            dsig = prsl (i, k) - prsl (i, k + 1)
            rdz = rdzt (i, k)
            tem1 = dsig * dkq (i, k) * rdz
            dsdz2 = tem1 * rdz
            au (i, k) = - dtodsd * dsdz2
            al (i, k) = - dtodsu * dsdz2
            ad (i, k) = ad (i, k) - au (i, k)
            ad (i, k + 1) = 1. - al (i, k)
            tem2 = dsig * rdz
            
            if (pcnvflg (i) .and. k < kpbl (i)) then
                ptem = 0.5 * tem2 * xmf (i, k)
                ptem1 = dtodsd * ptem
                ptem2 = dtodsu * ptem
                tem = tke (i, k) + tke (i, k + 1)
                ! kgao change
                ! ptem = qcko (i, k, ntke) + qcko (i, k + 1, ntke)
                ptem = qcko (i, k, ntrac) + qcko (i, k + 1, ntrac)
                f1 (i, k) = f1 (i, k) - (ptem - tem) * ptem1
                f1 (i, k + 1) = tke (i, k + 1) + (ptem - tem) * ptem2
            else
                f1 (i, k + 1) = tke (i, k + 1)
            endif
            
            if (scuflg (i)) then
                if (k >= mrad (i) .and. k < krad (i)) then
                    ptem = 0.5 * tem2 * xmfd (i, k)
                    ptem1 = dtodsd * ptem
                    ptem2 = dtodsu * ptem
                    tem = tke (i, k) + tke (i, k + 1)
                    ! kgao change
                    ! ptem = qcdo (i, k, ntke) + qcdo (i, k + 1, ntke)
                    ptem = qcdo (i, k, ntrac) + qcdo (i, k + 1, ntrac)
                    f1 (i, k) = f1 (i, k) + (ptem - tem) * ptem1
                    f1 (i, k + 1) = f1 (i, k + 1) - (ptem - tem) * ptem2
                endif
            endif
            
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! solve tridiagonal problem for tke
    ! -----------------------------------------------------------------------
    
    call tridit (im, km, 1, al, ad, au, f1, au, f1)
    
    ! -----------------------------------------------------------------------
    ! recover tendency of tke
    ! -----------------------------------------------------------------------
    
    do k = 1, km
        do i = 1, im
            ! fix negative tke
            f1 (i, k) = max (f1 (i, k), tkmin)
            q1g (i, k, ntrac) = f1 (i, k)
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute tridiagonal matrix elements for heat and moisture (and other tracers, except tke)
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        ad (i, 1) = 1.
        f1 (i, 1) = t1 (i, 1) + dtdz1 (i) * heat (i)
        f2 (i, 1) = q1g (i, 1, 1) + dtdz1 (i) * evap (i)
    enddo
    if (ntrac1 >= 2) then
        do kk = 2, ntrac1
            is = (kk - 1) * km
            do i = 1, im
                f2 (i, 1 + is) = q1g (i, 1, kk)
            enddo
        enddo
    endif
    
    do k = 1, km1
        do i = 1, im
            dtodsd = dt2 / del (i, k)
            dtodsu = dt2 / del (i, k + 1)
            dsig = prsl (i, k) - prsl (i, k + 1)
            rdz = rdzt (i, k)
            tem1 = dsig * dkt (i, k) * rdz
            dsdzt = tem1 * gocp
            dsdz2 = tem1 * rdz
            au (i, k) = - dtodsd * dsdz2
            al (i, k) = - dtodsu * dsdz2
            ad (i, k) = ad (i, k) - au (i, k)
            ad (i, k + 1) = 1. - al (i, k)
            tem2 = dsig * rdz
            
            if (pcnvflg (i) .and. k < kpbl (i)) then
                ptem = 0.5 * tem2 * xmf (i, k)
                ptem1 = dtodsd * ptem
                ptem2 = dtodsu * ptem
                tem = t1 (i, k) + t1 (i, k + 1)
                ptem = tcko (i, k) + tcko (i, k + 1)
                f1 (i, k) = f1 (i, k) + dtodsd * dsdzt - (ptem - tem) * ptem1
                f1 (i, k + 1) = t1 (i, k + 1) - dtodsu * dsdzt + (ptem - tem) * ptem2
                tem = q1g (i, k, 1) + q1g (i, k + 1, 1)
                ptem = qcko (i, k, 1) + qcko (i, k + 1, 1)
                f2 (i, k) = f2 (i, k) - (ptem - tem) * ptem1
                f2 (i, k + 1) = q1g (i, k + 1, 1) + (ptem - tem) * ptem2
            else
                f1 (i, k) = f1 (i, k) + dtodsd * dsdzt
                f1 (i, k + 1) = t1 (i, k + 1) - dtodsu * dsdzt
                f2 (i, k + 1) = q1g (i, k + 1, 1)
            endif
            
            if (scuflg (i)) then
                if (k >= mrad (i) .and. k < krad (i)) then
                    ptem = 0.5 * tem2 * xmfd (i, k)
                    ptem1 = dtodsd * ptem
                    ptem2 = dtodsu * ptem
                    ptem = tcdo (i, k) + tcdo (i, k + 1)
                    tem = t1 (i, k) + t1 (i, k + 1)
                    f1 (i, k) = f1 (i, k) + (ptem - tem) * ptem1
                    f1 (i, k + 1) = f1 (i, k + 1) - (ptem - tem) * ptem2
                    tem = q1g (i, k, 1) + q1g (i, k + 1, 1)
                    ptem = qcdo (i, k, 1) + qcdo (i, k + 1, 1)
                    f2 (i, k) = f2 (i, k) + (ptem - tem) * ptem1
                    f2 (i, k + 1) = f2 (i, k + 1) - (ptem - tem) * ptem2
                endif
            endif
        enddo
    enddo
    
    if (ntrac1 >= 2) then
        do kk = 2, ntrac1
            is = (kk - 1) * km
            do k = 1, km1
                do i = 1, im
                    if (pcnvflg (i) .and. k < kpbl (i)) then
                        dtodsd = dt2 / del (i, k)
                        dtodsu = dt2 / del (i, k + 1)
                        dsig = prsl (i, k) - prsl (i, k + 1)
                        tem = dsig * rdzt (i, k)
                        ptem = 0.5 * tem * xmf (i, k)
                        ptem1 = dtodsd * ptem
                        ptem2 = dtodsu * ptem
                        tem1 = qcko (i, k, kk) + qcko (i, k + 1, kk)
                        tem2 = q1g (i, k, kk) + q1g (i, k + 1, kk)
                        ! kgao note - turn off non - local mixing
                        f2 (i, k + is) = f2 (i, k + is) ! - (tem1 - tem2) * ptem1
                        f2 (i, k + 1 + is) = q1g (i, k + 1, kk) ! + (tem1 - tem2) * ptem2
                    else
                        f2 (i, k + 1 + is) = q1g (i, k + 1, kk)
                    endif
                    
                    if (scuflg (i)) then
                        if (k >= mrad (i) .and. k < krad (i)) then
                            dtodsd = dt2 / del (i, k)
                            dtodsu = dt2 / del (i, k + 1)
                            dsig = prsl (i, k) - prsl (i, k + 1)
                            tem = dsig * rdzt (i, k)
                            ptem = 0.5 * tem * xmfd (i, k)
                            ptem1 = dtodsd * ptem
                            ptem2 = dtodsu * ptem
                            tem1 = qcdo (i, k, kk) + qcdo (i, k + 1, kk)
                            tem2 = q1g (i, k, kk) + q1g (i, k + 1, kk)
                            ! kgao note - turn off non - local mixing
                            f2 (i, k + is) = f2 (i, k + is) ! + (tem1 - tem2) * ptem1
                            f2 (i, k + 1 + is) = f2 (i, k + 1 + is) ! - (tem1 - tem2) * ptem2
                        endif
                    endif
                    
                enddo
            enddo
        enddo
    endif
    
    ! -----------------------------------------------------------------------
    ! solve tridiagonal problem for heat and moisture
    ! -----------------------------------------------------------------------
    
    call tridin (im, km, ntrac1, al, ad, au, f1, f2, au, f1, f2)
    
    ! -----------------------------------------------------------------------
    ! recover tendencies of heat and moisture
    ! -----------------------------------------------------------------------
    
    do k = 1, km
        do i = 1, im
            tdt (i, k) = (f1 (i, k) - t1 (i, k)) * rdt
            qdt (i, k) = (f2 (i, k) - q1g (i, k, 1)) * rdt
            if (present (dtsfc)) dtsfc (i) = dtsfc (i) + cont * del (i, k) * tdt (i, k)
            if (present (dqsfc)) dqsfc (i) = dqsfc (i) + conq * del (i, k) * qdt (i, k)
            t1 (i, k) = f1 (i, k)
            q1g (i, k, 1) = f2 (i, k)
        enddo
    enddo
    
    if (ntrac1 >= 2) then
        do kk = 2, ntrac1
            is = (kk - 1) * km
            do k = 1, km
                do i = 1, im
                    q1g (i, k, kk) = f2 (i, k + is)
                enddo
            enddo
        enddo
    endif
    
    ! -----------------------------------------------------------------------
    ! kgao note - rearrange tracer tendencies
    ! -----------------------------------------------------------------------
    
    !if (ntrac >= 3) then
    if (ntke == ntrac) then ! tke is the last tracer
        q1 (:, :, :) = q1g (:, :, :)
    else ! tke is not
        do kk = 1, ntke - 1
            q1 (:, :, kk) = q1g (:, :, kk)
        enddo
        q1 (:, :, ntke) = q1g (:, :, ntrac)
        do kk = ntke + 1, ntrac
            q1 (:, :, kk) = q1g (:, :, kk - 1)
        enddo
    endif
    !endif
    
    ! -----------------------------------------------------------------------
    ! add tke dissipative heating to temperature tendency
    ! -----------------------------------------------------------------------
    
    if (dspheat) then
        do k = 1, km1
            do i = 1, im
                ! tem = min (diss (i, k), dspmax)
                ! ttend = tem / cp_air
                ttend = diss (i, k) / cp_air
                t1 (i, k) = t1 (i, k) + dspfac * ttend * dt2
            enddo
        enddo
    endif
    
    ! -----------------------------------------------------------------------
    ! compute tridiagonal matrix elements for momentum
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        ad (i, 1) = 1.0 + dtdz1 (i) * stress (i) / spd1 (i)
        f1 (i, 1) = u1 (i, 1)
        f2 (i, 1) = v1 (i, 1)
    enddo
    
    do k = 1, km1
        do i = 1, im
            dtodsd = dt2 / del (i, k)
            dtodsu = dt2 / del (i, k + 1)
            dsig = prsl (i, k) - prsl (i, k + 1)
            rdz = rdzt (i, k)
            tem1 = dsig * dku (i, k) * rdz
            dsdz2 = tem1 * rdz
            au (i, k) = - dtodsd * dsdz2
            al (i, k) = - dtodsu * dsdz2
            ad (i, k) = ad (i, k) - au (i, k)
            ad (i, k + 1) = 1. - al (i, k)
            tem2 = dsig * rdz
            
            if (pcnvflg (i) .and. k < kpbl (i)) then
                ptem = 0.5 * tem2 * xmf (i, k)
                ptem1 = dtodsd * ptem
                ptem2 = dtodsu * ptem
                tem = u1 (i, k) + u1 (i, k + 1)
                ptem = ucko (i, k) + ucko (i, k + 1)
                f1 (i, k) = f1 (i, k) - (ptem - tem) * ptem1
                f1 (i, k + 1) = u1 (i, k + 1) + (ptem - tem) * ptem2
                tem = v1 (i, k) + v1 (i, k + 1)
                ptem = vcko (i, k) + vcko (i, k + 1)
                f2 (i, k) = f2 (i, k) - (ptem - tem) * ptem1
                f2 (i, k + 1) = v1 (i, k + 1) + (ptem - tem) * ptem2
            else
                f1 (i, k + 1) = u1 (i, k + 1)
                f2 (i, k + 1) = v1 (i, k + 1)
            endif
            
            if (scuflg (i)) then
                if (k >= mrad (i) .and. k < krad (i)) then
                    ptem = 0.5 * tem2 * xmfd (i, k)
                    ptem1 = dtodsd * ptem
                    ptem2 = dtodsu * ptem
                    tem = u1 (i, k) + u1 (i, k + 1)
                    ptem = ucdo (i, k) + ucdo (i, k + 1)
                    f1 (i, k) = f1 (i, k) + (ptem - tem) * ptem1
                    f1 (i, k + 1) = f1 (i, k + 1) - (ptem - tem) * ptem2
                    tem = v1 (i, k) + v1 (i, k + 1)
                    ptem = vcdo (i, k) + vcdo (i, k + 1)
                    f2 (i, k) = f2 (i, k) + (ptem - tem) * ptem1
                    f2 (i, k + 1) = f2 (i, k + 1) - (ptem - tem) * ptem2
                endif
            endif
            
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! solve tridiagonal problem for momentum
    ! -----------------------------------------------------------------------
    
    call tridi2 (im, km, al, ad, au, f1, f2, au, f1, f2)
    
    ! -----------------------------------------------------------------------
    ! recover tendencies of momentum
    ! -----------------------------------------------------------------------
    
    do k = 1, km
        do i = 1, im
            udt (i, k) = (f1 (i, k) - u1 (i, k)) * rdt
            vdt (i, k) = (f2 (i, k) - v1 (i, k)) * rdt
            if (present (dusfc)) dusfc (i) = dusfc (i) + conw * del (i, k) * udt (i, k)
            if (present (dvsfc)) dvsfc (i) = dvsfc (i) + conw * del (i, k) * vdt (i, k)
            u1 (i, k) = f1 (i, k)
            v1 (i, k) = f2 (i, k)
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! pbl height for diagnostic purpose
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        hpbl (i) = hpblx (i)
        kpbl (i) = kpblx (i)
    enddo
    
end subroutine sa_tke_edmf
    
! =======================================================================
! solve tridiagonal problem for tke
! =======================================================================

subroutine tridit (l, n, nt, cl, cm, cu, rt, au, at)

    implicit none

    integer :: is, k, kk, n, nt, l, i

    real :: fk (l)
    
    real :: cl (l, 2:n), cm (l, n), cu (l, n - 1), &
        rt (l, n * nt), &
        au (l, n - 1), at (l, n * nt), &
        fkk (l, 2:n - 1)

    do i = 1, l
        fk (i) = 1. / cm (i, 1)
        au (i, 1) = fk (i) * cu (i, 1)
    enddo
    do k = 1, nt
        is = (k - 1) * n
        do i = 1, l
            at (i, 1 + is) = fk (i) * rt (i, 1 + is)
        enddo
    enddo
    do k = 2, n - 1
        do i = 1, l
            fkk (i, k) = 1. / (cm (i, k) - cl (i, k) * au (i, k - 1))
            au (i, k) = fkk (i, k) * cu (i, k)
        enddo
    enddo
    do kk = 1, nt
        is = (kk - 1) * n
        do k = 2, n - 1
            do i = 1, l
                at (i, k + is) = fkk (i, k) * (rt (i, k + is) - cl (i, k) * at (i, k + is - 1))
            enddo
        enddo
    enddo
    do i = 1, l
        fk (i) = 1. / (cm (i, n) - cl (i, n) * au (i, n - 1))
    enddo
    do k = 1, nt
        is = (k - 1) * n
        do i = 1, l
            at (i, n + is) = fk (i) * (rt (i, n + is) - cl (i, n) * at (i, n + is - 1))
        enddo
    enddo
    do kk = 1, nt
        is = (kk - 1) * n
        do k = n - 1, 1, - 1
            do i = 1, l
                at (i, k + is) = at (i, k + is) - au (i, k) * at (i, k + is + 1)
            enddo
        enddo
    enddo

end subroutine tridit
    
! =======================================================================
! edmf parameterization siebesma et al. (2007)
! =======================================================================

subroutine mfpblt (im, km, kmpbl, ntcw, ntrac1, delt, &
        cnvflg, zl, zm, q1, t1, u1, v1, plyr, pix, thlx, thvx, &
        gdx, hpbl, kpbl, vpert, buo, xmf, &
        tcko, qcko, ucko, vcko, xlamue)
    
    implicit none
    
    integer, intent (in) :: im, km, kmpbl, ntcw, ntrac1
    integer :: kpbl (im)

    logical :: cnvflg (im)

    real :: delt
    real :: q1 (im, km, ntrac1), &
        t1 (im, km), u1 (im, km), v1 (im, km), &
        plyr (im, km), pix (im, km), thlx (im, km), &
        thvx (im, km), zl (im, km), zm (im, km), &
        gdx (im), &
        hpbl (im), vpert (im), &
        buo (im, km), xmf (im, km), &
        tcko (im, km), qcko (im, km, ntrac1), &
        ucko (im, km), vcko (im, km), &
        xlamue (im, km - 1)
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i, j, k, n, ndc
    integer :: kpblx (im), kpbly (im)
    
    real :: dt2, dz, ce0, cm, &
        factor, gocp, &
        g, b1, f1, &
        bb1, bb2, &
        alp, a1, pgcon, &
        qmin, qlmin, xmmx, rbint, &
        tem, tem1, tem2, &
        ptem, ptem1, ptem2
    
    real :: elocp, el2orc, qs, es, &
        tlu, gamma, qlu, &
        thup, thvu, dq
    
    real :: rbdn (im), rbup (im), hpblx (im), &
        xlamuem (im, km - 1)
    
    real :: wu2 (im, km), thlu (im, km), &
        qtx (im, km), qtu (im, km)
    
    real :: xlamavg (im), sigma (im), &
        scaldfunc (im), sumx (im)
    
    logical :: totflg, flg (im)
    
    ! physical parameters
    parameter (g = grav)
    parameter (gocp = g / cp_air)
    parameter (elocp = hlv / cp_air, el2orc = hlv * hlv / (rvgas * cp_air))
    parameter (ce0 = 0.4, cm = 1.0)
    parameter (qmin = 1.e-8, qlmin = 1.e-12)
    parameter (alp = 1.0, pgcon = 0.55)
    parameter (a1 = 0.13, b1 = 0.5, f1 = 0.15)
    
    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo
    if (totflg) return
    
    dt2 = delt
    
    do k = 1, km
        do i = 1, im
            if (cnvflg (i)) then
                buo (i, k) = 0.
                wu2 (i, k) = 0.
                qtx (i, k) = q1 (i, k, 1) + q1 (i, k, ntcw)
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute thermal excess
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        if (cnvflg (i)) then
            ptem = alp * vpert (i)
            ptem = min (ptem, 3.0)
            thlu (i, 1) = thlx (i, 1) + ptem
            qtu (i, 1) = qtx (i, 1)
            buo (i, 1) = g * ptem / thvx (i, 1)
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute entrainment rate
    ! -----------------------------------------------------------------------
    
    do k = 1, kmpbl
        do i = 1, im
            if (cnvflg (i)) then
                dz = zl (i, k + 1) - zl (i, k)
                if (k < kpbl (i)) then
                    ptem = 1. / (zm (i, k) + dz)
                    tem = max ((hpbl (i) - zm (i, k) + dz), dz)
                    ptem1 = 1. / tem
                    xlamue (i, k) = ce0 * (ptem + ptem1)
                else
                    xlamue (i, k) = ce0 / dz
                endif
                xlamuem (i, k) = cm * xlamue (i, k)
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute buoyancy for updraft air parcel
    ! -----------------------------------------------------------------------
    
    do k = 2, kmpbl
        do i = 1, im
            if (cnvflg (i)) then
                dz = zl (i, k) - zl (i, k - 1)
                tem = 0.5 * xlamue (i, k - 1) * dz
                factor = 1. + tem
                
                thlu (i, k) = ((1. - tem) * thlu (i, k - 1) + tem * &
                     (thlx (i, k - 1) + thlx (i, k))) / factor
                qtu (i, k) = ((1. - tem) * qtu (i, k - 1) + tem * &
                     (qtx (i, k - 1) + qtx (i, k))) / factor
                
                tlu = thlu (i, k) / pix (i, k)
                es = 0.01 * mqs (tlu) ! mqs in pa
                qs = max (qmin, eps * es / (plyr (i, k) + epsm1 * es))
                dq = qtu (i, k) - qs
                
                if (dq > 0.) then
                    gamma = el2orc * qs / (tlu ** 2)
                    qlu = dq / (1. + gamma)
                    qtu (i, k) = qs + qlu
                    tem1 = 1. + zvir * qs - qlu
                    thup = thlu (i, k) + pix (i, k) * elocp * qlu
                    thvu = thup * tem1
                else
                    tem1 = 1. + zvir * qtu (i, k)
                    thvu = thlu (i, k) * tem1
                endif
                buo (i, k) = g * (thvu / thvx (i, k) - 1.)
                
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute updraft velocity square (wu2)
    ! -----------------------------------------------------------------------
    
    ! tem = 1. - 2. * f1
    ! bb1 = 2. * b1 / tem
    ! bb2 = 2. / tem
    ! from soares et al. (2004, qjrms)
    ! bb1 = 2.
    ! bb2 = 4.
    
    ! from bretherton et al. (2004, mwr)
    ! bb1 = 4.
    ! bb2 = 2.
    
    ! from our tuning
    bb1 = 2.0
    bb2 = 4.0
    
    do i = 1, im
        if (cnvflg (i)) then
            dz = zm (i, 1)
            tem = 0.5 * bb1 * xlamue (i, 1) * dz
            tem1 = bb2 * buo (i, 1) * dz
            ptem1 = 1. + tem
            wu2 (i, 1) = tem1 / ptem1
        endif
    enddo
    do k = 2, kmpbl
        do i = 1, im
            if (cnvflg (i)) then
                dz = zm (i, k) - zm (i, k - 1)
                tem = 0.25 * bb1 * (xlamue (i, k) + xlamue (i, k - 1)) * dz
                tem1 = bb2 * buo (i, k) * dz
                ptem = (1. - tem) * wu2 (i, k - 1)
                ptem1 = 1. + tem
                wu2 (i, k) = (ptem + tem1) / ptem1
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! update pbl height as the height where updraft velocity vanishes
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        flg (i) = .true.
        kpbly (i) = kpbl (i)
        if (cnvflg (i)) then
            flg (i) = .false.
            rbup (i) = wu2 (i, 1)
        endif
    enddo
    do k = 2, kmpbl
        do i = 1, im
            if (.not.flg (i)) then
                rbdn (i) = rbup (i)
                rbup (i) = wu2 (i, k)
                kpblx (i) = k
                flg (i) = rbup (i) .le.0.
            endif
        enddo
    enddo
    do i = 1, im
        if (cnvflg (i)) then
            k = kpblx (i)
            if (rbdn (i) <= 0.) then
                rbint = 0.
            elseif (rbup (i) >= 0.) then
                rbint = 1.
            else
                rbint = rbdn (i) / (rbdn (i) - rbup (i))
            endif
            hpblx (i) = zm (i, k - 1) + rbint * (zm (i, k) - zm (i, k - 1))
        endif
    enddo
    
    do i = 1, im
        if (cnvflg (i)) then
            if (kpbl (i) > kpblx (i)) then
                kpbl (i) = kpblx (i)
                hpbl (i) = hpblx (i)
            endif
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! update entrainment rate
    ! -----------------------------------------------------------------------
    
    do k = 1, kmpbl
        do i = 1, im
            if (cnvflg (i) .and. kpbly (i) > kpblx (i)) then
                dz = zl (i, k + 1) - zl (i, k)
                if (k < kpbl (i)) then
                    ptem = 1. / (zm (i, k) + dz)
                    tem = max ((hpbl (i) - zm (i, k) + dz), dz)
                    ptem1 = 1. / tem
                    xlamue (i, k) = ce0 * (ptem + ptem1)
                else
                    xlamue (i, k) = ce0 / dz
                endif
                xlamuem (i, k) = cm * xlamue (i, k)
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute entrainment rate averaged over the whole pbl
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        xlamavg (i) = 0.
        sumx (i) = 0.
    enddo
    do k = 1, kmpbl
        do i = 1, im
            if (cnvflg (i) .and. k < kpbl (i)) then
                dz = zl (i, k + 1) - zl (i, k)
                xlamavg (i) = xlamavg (i) + xlamue (i, k) * dz
                sumx (i) = sumx (i) + dz
            endif
        enddo
    enddo
    do i = 1, im
        if (cnvflg (i)) then
            xlamavg (i) = xlamavg (i) / sumx (i)
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! updraft mass flux as a function of updraft velocity profile
    ! -----------------------------------------------------------------------
    
    do k = 1, kmpbl
        do i = 1, im
            if (cnvflg (i) .and. k < kpbl (i)) then
                if (wu2 (i, k) > 0.) then
                    tem = sqrt (wu2 (i, k))
                else
                    tem = 0.
                endif
                xmf (i, k) = a1 * tem
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute updraft fraction as a function of mean entrainment rate
    ! (grell & freitas, 2014)
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        if (cnvflg (i)) then
            tem = 0.2 / xlamavg (i)
            tem1 = 3.14 * tem * tem
            sigma (i) = tem1 / (gdx (i) * gdx (i))
            sigma (i) = max (sigma (i), 0.001)
            sigma (i) = min (sigma (i), 0.999)
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute scale - aware function based on arakawa & wu (2013)
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        if (cnvflg (i)) then
            if (sigma (i) > a1) then
                scaldfunc (i) = (1. - sigma (i)) * (1. - sigma (i))
                scaldfunc (i) = max (min (scaldfunc (i), 1.0), 0.)
            else
                scaldfunc (i) = 1.0
            endif
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! final scale - aware updraft mass flux
    ! -----------------------------------------------------------------------
    
    do k = 1, kmpbl
        do i = 1, im
            if (cnvflg (i) .and. k < kpbl (i)) then
                xmf (i, k) = scaldfunc (i) * xmf (i, k)
                dz = zl (i, k + 1) - zl (i, k)
                xmmx = dz / dt2
                xmf (i, k) = min (xmf (i, k), xmmx)
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute updraft property using updated entranment rate
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        if (cnvflg (i)) then
            thlu (i, 1) = thlx (i, 1)
        endif
    enddo
    
    ! do i = 1, im
    ! if (cnvflg (i)) then
    ! ptem1 = max (qcko (i, 1, ntcw), 0.)
    ! tlu = thlu (i, 1) / pix (i, 1)
    ! tcko (i, 1) = tlu + elocp * ptem1
    ! endif
    ! enddo
    
    do k = 2, kmpbl
        do i = 1, im
            if (cnvflg (i) .and. k <= kpbl (i)) then
                dz = zl (i, k) - zl (i, k - 1)
                tem = 0.5 * xlamue (i, k - 1) * dz
                factor = 1. + tem
                
                thlu (i, k) = ((1. - tem) * thlu (i, k - 1) + tem * &
                     (thlx (i, k - 1) + thlx (i, k))) / factor
                qtu (i, k) = ((1. - tem) * qtu (i, k - 1) + tem * &
                     (qtx (i, k - 1) + qtx (i, k))) / factor
                
                tlu = thlu (i, k) / pix (i, k)
                es = 0.01 * mqs (tlu) ! mqs in pa
                qs = max (qmin, eps * es / (plyr (i, k) + epsm1 * es))
                dq = qtu (i, k) - qs
                
                if (dq > 0.) then
                    gamma = el2orc * qs / (tlu ** 2)
                    qlu = dq / (1. + gamma)
                    qtu (i, k) = qs + qlu
                    qcko (i, k, 1) = qs
                    qcko (i, k, ntcw) = qlu
                    tcko (i, k) = tlu + elocp * qlu
                else
                    qcko (i, k, 1) = qtu (i, k)
                    qcko (i, k, ntcw) = 0.
                    tcko (i, k) = tlu
                endif
                
            endif
        enddo
    enddo
    
    do k = 2, kmpbl
        do i = 1, im
            if (cnvflg (i) .and. k <= kpbl (i)) then
                dz = zl (i, k) - zl (i, k - 1)
                tem = 0.5 * xlamuem (i, k - 1) * dz
                factor = 1. + tem
                ptem = tem + pgcon
                ptem1 = tem - pgcon
                ucko (i, k) = ((1. - tem) * ucko (i, k - 1) + ptem * u1 (i, k) &
                     + ptem1 * u1 (i, k - 1)) / factor
                vcko (i, k) = ((1. - tem) * vcko (i, k - 1) + ptem * v1 (i, k) &
                     + ptem1 * v1 (i, k - 1)) / factor
            endif
        enddo
    enddo
    
    if (ntcw > 2) then
        
        do n = 2, ntcw - 1
            do k = 2, kmpbl
                do i = 1, im
                    if (cnvflg (i) .and. k <= kpbl (i)) then
                        dz = zl (i, k) - zl (i, k - 1)
                        tem = 0.5 * xlamue (i, k - 1) * dz
                        factor = 1. + tem
                        
                        qcko (i, k, n) = ((1. - tem) * qcko (i, k - 1, n) + tem * &
                             (q1 (i, k, n) + q1 (i, k - 1, n))) / factor
                    endif
                enddo
            enddo
        enddo
        
    endif
    
    ndc = ntrac1 - ntcw
    
    if (ndc > 0) then
        
        do n = ntcw + 1, ntrac1
            do k = 2, kmpbl
                do i = 1, im
                    if (cnvflg (i) .and. k <= kpbl (i)) then
                        dz = zl (i, k) - zl (i, k - 1)
                        tem = 0.5 * xlamue (i, k - 1) * dz
                        factor = 1. + tem
                        
                        qcko (i, k, n) = ((1. - tem) * qcko (i, k - 1, n) + tem * &
                             (q1 (i, k, n) + q1 (i, k - 1, n))) / factor
                    endif
                enddo
            enddo
        enddo
        
    endif
    
end subroutine mfpblt  
    
! =======================================================================
! mass - flux parameterization for stratocumulus - top - induced turbulence mixing
! =======================================================================
    
subroutine mfscu (im, km, kmscu, ntcw, ntrac1, delt, &
        cnvflg, zl, zm, q1, t1, u1, v1, plyr, pix, &
        thlx, thvx, thlvx, gdx, thetae, radj, &
        krad, mrad, radmin, buo, xmfd, &
        tcdo, qcdo, ucdo, vcdo, xlamde)
    
    implicit none
    
    integer, intent (in) :: im, km, kmscu, ntcw, ntrac1
    integer :: krad (im), mrad (im)
    
    logical :: cnvflg (im)

    real :: delt
    real :: q1 (im, km, ntrac1), t1 (im, km), &
        u1 (im, km), v1 (im, km), &
        plyr (im, km), pix (im, km), &
        thlx (im, km), &
        thvx (im, km), thlvx (im, km), &
        gdx (im), radj (im), &
        zl (im, km), zm (im, km), &
        thetae (im, km), radmin (im), &
        buo (im, km), xmfd (im, km), &
        tcdo (im, km), qcdo (im, km, ntrac1), &
        ucdo (im, km), vcdo (im, km), &
        xlamde (im, km - 1)
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i, j, indx, k, n, kk, ndc

    integer :: krad1 (im), mradx (im), mrady (im)
    
    real :: dt2, dz, ce0, cm, &
        gocp, factor, g, tau, &
        b1, f1, bb1, bb2, &
        a1, a2, a11, a22, &
        cteit, pgcon, &
        qmin, qlmin, &
        xmmx, tem, tem1, tem2, &
        ptem, ptem1, ptem2
    
    real :: elocp, el2orc, qs, es, &
        tld, gamma, qld, thdn, &
        thvd, dq
    
    real :: wd2 (im, km), thld (im, km), &
        qtx (im, km), qtd (im, km), &
        thlvd (im), hrad (im), &
        xlamdem (im, km - 1), ra1 (im), ra2 (im)
    
    real :: xlamavg (im), sigma (im), &
        scaldfunc (im), sumx (im)
    
    logical :: totflg, flg (im)
    
    real :: actei, cldtime
    
    ! physical parameters
    parameter (g = grav)
    parameter (gocp = g / cp_air)
    parameter (elocp = hlv / cp_air, el2orc = hlv * hlv / (rvgas * cp_air))
    parameter (ce0 = 0.4, cm = 1.0, pgcon = 0.55)
    parameter (qmin = 1.e-8, qlmin = 1.e-12)
    parameter (b1 = 0.45, f1 = 0.15)
    parameter (a1 = 0.12, a2 = 0.5)
    parameter (a11 = 0.2, a22 = 1.0)
    parameter (cldtime = 500.)
    parameter (actei = 0.7)
    ! parameter (actei = 0.23)
    
    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo
    if (totflg) return
    
    dt2 = delt
    
    do k = 1, km
        do i = 1, im
            if (cnvflg (i)) then
                buo (i, k) = 0.
                wd2 (i, k) = 0.
                qtx (i, k) = q1 (i, k, 1) + q1 (i, k, ntcw)
            endif
        enddo
    enddo
    
    do i = 1, im
        if (cnvflg (i)) then
            hrad (i) = zm (i, krad (i))
            krad1 (i) = krad (i) - 1
        endif
    enddo
    
    do i = 1, im
        if (cnvflg (i)) then
            k = krad (i)
            tem = zm (i, k + 1) - zm (i, k)
            tem1 = cldtime * radmin (i) / tem
            tem1 = max (tem1, - 3.0)
            thld (i, k) = thlx (i, k) + tem1
            qtd (i, k) = qtx (i, k)
            thlvd (i) = thlvx (i, k) + tem1
            buo (i, k) = - g * tem1 / thvx (i, k)
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! specify downdraft fraction
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        if (cnvflg (i)) then
            ra1 (i) = a1
            ra2 (i) = a11
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! if the condition for cloud - top instability is met,
    ! increase downdraft fraction
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        if (cnvflg (i)) then
            k = krad (i)
            tem = thetae (i, k) - thetae (i, k + 1)
            tem1 = qtx (i, k) - qtx (i, k + 1)
            if (tem > 0. .and. tem1 > 0.) then
                cteit = cp_air * tem / (hlv * tem1)
                if (cteit > actei) then
                    ra1 (i) = a2
                    ra2 (i) = a22
                endif
            endif
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute radiative flux jump at stratocumulus top
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        if (cnvflg (i)) then
            radj (i) = - ra2 (i) * radmin (i)
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! first - quess level of downdraft extension (mrad)
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        flg (i) = cnvflg (i)
        mrad (i) = krad (i)
    enddo
    do k = kmscu, 1, - 1
        do i = 1, im
            if (flg (i) .and. k < krad (i)) then
                if (thlvd (i) <= thlvx (i, k)) then
                    mrad (i) = k
                else
                    flg (i) = .false.
                endif
            endif
        enddo
    enddo
    do i = 1, im
        if (cnvflg (i)) then
            kk = krad (i) - mrad (i)
            if (kk < 1) cnvflg (i) = .false.
        endif
    enddo
    
    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo
    if (totflg) return
    
    ! -----------------------------------------------------------------------
    ! compute entrainment rate
    ! -----------------------------------------------------------------------
    
    do k = 1, kmscu
        do i = 1, im
            if (cnvflg (i)) then
                dz = zl (i, k + 1) - zl (i, k)
                if (k >= mrad (i) .and. k < krad (i)) then
                    if (mrad (i) == 1) then
                        ptem = 1. / (zm (i, k) + dz)
                    else
                        ptem = 1. / (zm (i, k) - zm (i, mrad (i) - 1) + dz)
                    endif
                    tem = max ((hrad (i) - zm (i, k) + dz), dz)
                    ptem1 = 1. / tem
                    xlamde (i, k) = ce0 * (ptem + ptem1)
                else
                    xlamde (i, k) = ce0 / dz
                endif
                xlamdem (i, k) = cm * xlamde (i, k)
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute buoyancy for downdraft air parcel
    ! -----------------------------------------------------------------------
    
    do k = kmscu, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. k < krad (i)) then
                dz = zl (i, k + 1) - zl (i, k)
                tem = 0.5 * xlamde (i, k) * dz
                factor = 1. + tem
                
                thld (i, k) = ((1. - tem) * thld (i, k + 1) + tem * &
                     (thlx (i, k) + thlx (i, k + 1))) / factor
                qtd (i, k) = ((1. - tem) * qtd (i, k + 1) + tem * &
                     (qtx (i, k) + qtx (i, k + 1))) / factor
                
                tld = thld (i, k) / pix (i, k)
                es = 0.01 * mqs (tld) ! mqs in pa
                qs = max (qmin, eps * es / (plyr (i, k) + epsm1 * es))
                dq = qtd (i, k) - qs
                
                if (dq > 0.) then
                    gamma = el2orc * qs / (tld ** 2)
                    qld = dq / (1. + gamma)
                    qtd (i, k) = qs + qld
                    tem1 = 1. + zvir * qs - qld
                    thdn = thld (i, k) + pix (i, k) * elocp * qld
                    thvd = thdn * tem1
                else
                    tem1 = 1. + zvir * qtd (i, k)
                    thvd = thld (i, k) * tem1
                endif
                buo (i, k) = g * (1. - thvd / thvx (i, k))
                
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute downdraft velocity square (wd2)
    ! -----------------------------------------------------------------------
    
    ! tem = 1. - 2. * f1
    ! bb1 = 2. * b1 / tem
    ! bb2 = 2. / tem
    ! from soares et al. (2004, qjrms)
    ! bb1 = 2.
    ! bb2 = 4.
    
    ! from bretherton et al. (2004, mwr)
    ! bb1 = 4.
    ! bb2 = 2.
    
    ! from our tuning
    bb1 = 2.0
    bb2 = 4.0
    
    do i = 1, im
        if (cnvflg (i)) then
            k = krad1 (i)
            dz = zm (i, k + 1) - zm (i, k)
            ! tem = 0.25 * bb1 * (xlamde (i, k) + xlamde (i, k + 1)) * dz
            tem = 0.5 * bb1 * xlamde (i, k) * dz
            tem1 = bb2 * buo (i, k + 1) * dz
            ptem1 = 1. + tem
            wd2 (i, k) = tem1 / ptem1
        endif
    enddo
    do k = kmscu, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. k < krad1 (i)) then
                dz = zm (i, k + 1) - zm (i, k)
                tem = 0.25 * bb1 * (xlamde (i, k) + xlamde (i, k + 1)) * dz
                tem1 = bb2 * buo (i, k + 1) * dz
                ptem = (1. - tem) * wd2 (i, k + 1)
                ptem1 = 1. + tem
                wd2 (i, k) = (ptem + tem1) / ptem1
            endif
        enddo
    enddo
    
    do i = 1, im
        flg (i) = cnvflg (i)
        mrady (i) = mrad (i)
        if (flg (i)) mradx (i) = krad (i)
    enddo
    do k = kmscu, 1, - 1
        do i = 1, im
            if (flg (i) .and. k < krad (i)) then
                if (wd2 (i, k) > 0.) then
                    mradx (i) = k
                else
                    flg (i) = .false.
                endif
            endif
        enddo
    enddo
    
    do i = 1, im
        if (cnvflg (i)) then
            if (mrad (i) < mradx (i)) then
                mrad (i) = mradx (i)
            endif
        endif
    enddo
    
    do i = 1, im
        if (cnvflg (i)) then
            kk = krad (i) - mrad (i)
            if (kk < 1) cnvflg (i) = .false.
        endif
    enddo
    
    totflg = .true.
    do i = 1, im
        totflg = totflg .and. (.not. cnvflg (i))
    enddo
    if (totflg) return
    
    ! -----------------------------------------------------------------------
    ! update entrainment rate
    ! -----------------------------------------------------------------------
    
    do k = 1, kmscu
        do i = 1, im
            if (cnvflg (i) .and. mrady (i) < mradx (i)) then
                dz = zl (i, k + 1) - zl (i, k)
                if (k >= mrad (i) .and. k < krad (i)) then
                    if (mrad (i) == 1) then
                        ptem = 1. / (zm (i, k) + dz)
                    else
                        ptem = 1. / (zm (i, k) - zm (i, mrad (i) - 1) + dz)
                    endif
                    tem = max ((hrad (i) - zm (i, k) + dz), dz)
                    ptem1 = 1. / tem
                    xlamde (i, k) = ce0 * (ptem + ptem1)
                else
                    xlamde (i, k) = ce0 / dz
                endif
                xlamdem (i, k) = cm * xlamde (i, k)
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute entrainment rate averaged over the whole downdraft layers
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        xlamavg (i) = 0.
        sumx (i) = 0.
    enddo
    do k = kmscu, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. (k >= mrad (i) .and. k < krad (i))) then
                dz = zl (i, k + 1) - zl (i, k)
                xlamavg (i) = xlamavg (i) + xlamde (i, k) * dz
                sumx (i) = sumx (i) + dz
            endif
        enddo
    enddo
    do i = 1, im
        if (cnvflg (i)) then
            xlamavg (i) = xlamavg (i) / sumx (i)
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute downdraft mass flux
    ! -----------------------------------------------------------------------
    
    do k = kmscu, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. (k >= mrad (i) .and. k < krad (i))) then
                if (wd2 (i, k) > 0.) then
                    tem = sqrt (wd2 (i, k))
                else
                    tem = 0.
                endif
                xmfd (i, k) = ra1 (i) * tem
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute downdraft fraction as a function of mean entrainment rate
    ! (grell & freitas, 2014)
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        if (cnvflg (i)) then
            tem = 0.2 / xlamavg (i)
            tem1 = 3.14 * tem * tem
            sigma (i) = tem1 / (gdx (i) * gdx (i))
            sigma (i) = max (sigma (i), 0.001)
            sigma (i) = min (sigma (i), 0.999)
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute scale - aware function based on arakawa & wu (2013)
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        if (cnvflg (i)) then
            if (sigma (i) > ra1 (i)) then
                scaldfunc (i) = (1. - sigma (i)) * (1. - sigma (i))
                scaldfunc (i) = max (min (scaldfunc (i), 1.0), 0.)
            else
                scaldfunc (i) = 1.0
            endif
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! final scale - aware downdraft mass flux
    ! -----------------------------------------------------------------------
    
    do k = kmscu, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. (k >= mrad (i) .and. k < krad (i))) then
                xmfd (i, k) = scaldfunc (i) * xmfd (i, k)
                dz = zl (i, k + 1) - zl (i, k)
                xmmx = dz / dt2
                xmfd (i, k) = min (xmfd (i, k), xmmx)
            endif
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute downdraft property using updated entranment rate
    ! -----------------------------------------------------------------------
    
    do i = 1, im
        if (cnvflg (i)) then
            k = krad (i)
            thld (i, k) = thlx (i, k)
        endif
    enddo
    
    ! do i = 1, im
    ! if (cnvflg (i)) then
    ! k = krad (i)
    ! ptem1 = max (qcdo (i, k, ntcw), 0.)
    ! tld = thld (i, k) / pix (i, k)
    ! tcdo (i, k) = tld + elocp * ptem1
    ! qcdo (i, k, 1) = qcdo (i, k, 1) + 0.2 * qcdo (i, k, 1)
    ! qcdo (i, k, ntcw) = qcdo (i, k, ntcw) + 0.2 * qcdo (i, k, ntcw)
    ! endif
    ! enddo
    
    do k = kmscu, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. (k >= mrad (i) .and. k < krad (i))) then
                dz = zl (i, k + 1) - zl (i, k)
                tem = 0.5 * xlamde (i, k) * dz
                factor = 1. + tem
                
                thld (i, k) = ((1. - tem) * thld (i, k + 1) + tem * &
                     (thlx (i, k) + thlx (i, k + 1))) / factor
                qtd (i, k) = ((1. - tem) * qtd (i, k + 1) + tem * &
                     (qtx (i, k) + qtx (i, k + 1))) / factor
                
                tld = thld (i, k) / pix (i, k)
                es = 0.01 * mqs (tld) ! mqs in pa
                qs = max (qmin, eps * es / (plyr (i, k) + epsm1 * es))
                dq = qtd (i, k) - qs
                
                if (dq > 0.) then
                    gamma = el2orc * qs / (tld ** 2)
                    qld = dq / (1. + gamma)
                    qtd (i, k) = qs + qld
                    qcdo (i, k, 1) = qs
                    qcdo (i, k, ntcw) = qld
                    tcdo (i, k) = tld + elocp * qld
                else
                    qcdo (i, k, 1) = qtd (i, k)
                    qcdo (i, k, ntcw) = 0.
                    tcdo (i, k) = tld
                endif
                
            endif
        enddo
    enddo
    
    do k = kmscu, 1, - 1
        do i = 1, im
            if (cnvflg (i) .and. k < krad (i)) then
                if (k >= mrad (i)) then
                    dz = zl (i, k + 1) - zl (i, k)
                    tem = 0.5 * xlamdem (i, k) * dz
                    factor = 1. + tem
                    ptem = tem - pgcon
                    ptem1 = tem + pgcon
                    
                    ucdo (i, k) = ((1. - tem) * ucdo (i, k + 1) + ptem * u1 (i, k + 1) &
                         + ptem1 * u1 (i, k)) / factor
                    vcdo (i, k) = ((1. - tem) * vcdo (i, k + 1) + ptem * v1 (i, k + 1) &
                         + ptem1 * v1 (i, k)) / factor
                endif
            endif
        enddo
    enddo
    
    if (ntcw > 2) then
        
        do n = 2, ntcw - 1
            do k = kmscu, 1, - 1
                do i = 1, im
                    if (cnvflg (i) .and. k < krad (i)) then
                        if (k >= mrad (i)) then
                            dz = zl (i, k + 1) - zl (i, k)
                            tem = 0.5 * xlamde (i, k) * dz
                            factor = 1. + tem
                            
                            qcdo (i, k, n) = ((1. - tem) * qcdo (i, k + 1, n) + tem * &
                                 (q1 (i, k, n) + q1 (i, k + 1, n))) / factor
                        endif
                    endif
                enddo
            enddo
        enddo
        
    endif
    
    ndc = ntrac1 - ntcw
    
    if (ndc > 0) then
        
        do n = ntcw + 1, ntrac1
            do k = kmscu, 1, - 1
                do i = 1, im
                    if (cnvflg (i) .and. k < krad (i)) then
                        if (k >= mrad (i)) then
                            dz = zl (i, k + 1) - zl (i, k)
                            tem = 0.5 * xlamde (i, k) * dz
                            factor = 1. + tem
                            
                            qcdo (i, k, n) = ((1. - tem) * qcdo (i, k + 1, n) + tem * &
                                 (q1 (i, k, n) + q1 (i, k + 1, n))) / factor
                        endif
                    endif
                enddo
            enddo
        enddo
        
    endif
    
end subroutine mfscu

! =======================================================================
! routine to solve the tridiagonal system to calculate temperature and
! moisture at \f$ t + \delta t \f$; part of two - part process to 
! calculate time tendencies due to vertical diffusion.
! =======================================================================

subroutine tridi2 (l, n, cl, cm, cu, r1, r2, au, a1, a2)
    
    implicit none

    integer :: k, n, l, i

    real :: fk
    
    real :: cl (l, 2:n), cm (l, n), cu (l, n - 1), r1 (l, n), r2 (l, n), &
        au (l, n - 1), a1 (l, n), a2 (l, n)

    do i = 1, l
        fk = 1. / cm (i, 1)
        au (i, 1) = fk * cu (i, 1)
        a1 (i, 1) = fk * r1 (i, 1)
        a2 (i, 1) = fk * r2 (i, 1)
    enddo
    do k = 2, n - 1
        do i = 1, l
            fk = 1. / (cm (i, k) - cl (i, k) * au (i, k - 1))
            au (i, k) = fk * cu (i, k)
            a1 (i, k) = fk * (r1 (i, k) - cl (i, k) * a1 (i, k - 1))
            a2 (i, k) = fk * (r2 (i, k) - cl (i, k) * a2 (i, k - 1))
        enddo
    enddo
    do i = 1, l
        fk = 1. / (cm (i, n) - cl (i, n) * au (i, n - 1))
        a1 (i, n) = fk * (r1 (i, n) - cl (i, n) * a1 (i, n - 1))
        a2 (i, n) = fk * (r2 (i, n) - cl (i, n) * a2 (i, n - 1))
    enddo
    do k = n - 1, 1, - 1
        do i = 1, l
            a1 (i, k) = a1 (i, k) - au (i, k) * a1 (i, k + 1)
            a2 (i, k) = a2 (i, k) - au (i, k) * a2 (i, k + 1)
        enddo
    enddo

end subroutine tridi2

! =======================================================================
! routine to solve the tridiagonal system to calculate u - and v - 
! momentum at \f$ t + \delta t \f$; part of two - part process to 
! calculate time tendencies due to vertical diffusion.
! =======================================================================
    
subroutine tridin (l, n, nt, cl, cm, cu, r1, r2, au, a1, a2)
    
    implicit none

    integer :: is, k, kk, n, nt, l, i

    real :: fk (l)
    
    real :: cl (l, 2:n), cm (l, n), cu (l, n - 1), &
        r1 (l, n), r2 (l, n * nt), &
        au (l, n - 1), a1 (l, n), a2 (l, n * nt), &
        fkk (l, 2:n - 1)

    do i = 1, l
        fk (i) = 1. / cm (i, 1)
        au (i, 1) = fk (i) * cu (i, 1)
        a1 (i, 1) = fk (i) * r1 (i, 1)
    enddo
    do k = 1, nt
        is = (k - 1) * n
        do i = 1, l
            a2 (i, 1 + is) = fk (i) * r2 (i, 1 + is)
        enddo
    enddo
    do k = 2, n - 1
        do i = 1, l
            fkk (i, k) = 1. / (cm (i, k) - cl (i, k) * au (i, k - 1))
            au (i, k) = fkk (i, k) * cu (i, k)
            a1 (i, k) = fkk (i, k) * (r1 (i, k) - cl (i, k) * a1 (i, k - 1))
        enddo
    enddo
    do kk = 1, nt
        is = (kk - 1) * n
        do k = 2, n - 1
            do i = 1, l
                a2 (i, k + is) = fkk (i, k) * (r2 (i, k + is) - cl (i, k) * a2 (i, k + is - 1))
            enddo
        enddo
    enddo
    do i = 1, l
        fk (i) = 1. / (cm (i, n) - cl (i, n) * au (i, n - 1))
        a1 (i, n) = fk (i) * (r1 (i, n) - cl (i, n) * a1 (i, n - 1))
    enddo
    do k = 1, nt
        is = (k - 1) * n
        do i = 1, l
            a2 (i, n + is) = fk (i) * (r2 (i, n + is) - cl (i, n) * a2 (i, n + is - 1))
        enddo
    enddo
    do k = n - 1, 1, - 1
        do i = 1, l
            a1 (i, k) = a1 (i, k) - au (i, k) * a1 (i, k + 1)
        enddo
    enddo
    do kk = 1, nt
        is = (kk - 1) * n
        do k = n - 1, 1, - 1
            do i = 1, l
                a2 (i, k + is) = a2 (i, k + is) - au (i, k) * a2 (i, k + is + 1)
            enddo
        enddo
    enddo

end subroutine tridin

end module sa_tke_edmf_mod
