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
! GFDL Cloud Microphysics Package (GFDL MP)
! The algorithms are originally derived from Lin et al. (1983).
! Most of the key elements have been simplified / improved.
! This code at this stage bears little to no similarity to the original Lin MP in ZETAC.
! Developers: Shian-Jiann Lin, Linjiong Zhou, and the GFDL FV3 Team
! References: Chen and Lin (2011, 2013), Zhou et al. (2019), Harris et al. (2020)
! =======================================================================

module gfdl_cld_mp_mod
    
    use machine, only: r_grid => kind_phys
    
    implicit none
    
    private
    
    ! -----------------------------------------------------------------------
    ! public subroutines, functions, and variables
    ! -----------------------------------------------------------------------
    
    public :: gfdl_cld_mp_init
    public :: gfdl_cld_mp_driver
    public :: gfdl_cld_mp_end
    public :: fast_sat_adj, cld_eff_rad, rad_ref
    public :: qsmith_init, qs_moist, qsmith, wqs, wqs_moist
    public :: c_liq, c_ice, rhow, wet_bulb
    
    ! -----------------------------------------------------------------------
    ! initialization conditions
    ! -----------------------------------------------------------------------
    
    logical :: tables_are_initialized = .false. ! initialize satuation tables
    
    ! -----------------------------------------------------------------------
    ! physics constants
    ! -----------------------------------------------------------------------
    
    real, parameter :: grav = 9.80665 ! acceleration due to gravity (m/s^2)

    real, parameter :: rgrav = 1.0 / grav ! inversion of gravity acceleration (s^2/m)

    real, parameter :: pi = 4.0 * atan (1.0) ! ratio of circle circumference to diameter

    real, parameter :: rdgas = 287.05 ! gas constant for dry air (J/kg/K)
    real, parameter :: rvgas = 461.50 ! gas constant for water vapor (J/kg/K)

    real, parameter :: eps = rdgas / rvgas ! 0.6219934995
    real, parameter :: zvir = rvgas / rdgas - 1. ! 0.6077338443
    
    real, parameter :: tice = 273.16 ! freezing temperature (K)
    
    real, parameter :: cp_air = 1.0046e3 ! heat capacity of dry air at constant pressure (J/kg/K)
    real, parameter :: cp_vap = 4.0 * rvgas ! 1846.0, heat capacity of water vapore at constnat pressure (J/kg/K)
    real, parameter :: cv_air = cp_air - rdgas ! 717.55, heat capacity of dry air at constant volume (J/kg/K)
    real, parameter :: cv_vap = 3.0 * rvgas ! 1384.5, heat capacity of water vapor at constant volume (J/kg/K)
    
    real, parameter :: c_ice = 2.106e3 ! heat capacity of ice at 0 deg C (J/kg/K)
    real, parameter :: c_liq = 4.218e3 ! heat capacity of water at 0 deg C (J/kg/K)

    real, parameter :: dc_vap = cp_vap - c_liq ! - 2.372e3, isobaric heating / cooling (J/kg/K)
    real, parameter :: dc_ice = c_liq - c_ice ! 2.112e3, isobaric heating / colling (J/kg/K)
    
    real, parameter :: hlv = 2.5e6 ! latent heat of evaporation (J/kg)
    real, parameter :: hlf = 3.3358e5 ! latent heat of fusion (J/kg)

    real, parameter :: hlv0 = hlv ! evaporation latent heat coefficient at 0 deg C (J/kg)
    real, parameter :: hlf0 = hlf ! fussion latent heat coefficient at 0 deg C (J/kg)
    
    real, parameter :: lv0 = hlv0 - dc_vap * tice ! 3.14893552e6, evaporation latent heat coeff. at 0 deg K (J/kg)
    real, parameter :: li0 = hlf0 - dc_ice * tice ! - 2.2691392e5, fussion latent heat coeff. at 0 deg K (J/kg)
    
    real (kind = r_grid), parameter :: gam263 = 1.4625080991084523 ! Gamma function (2.63)
    real (kind = r_grid), parameter :: gam275 = 1.6083594219855455 ! Gamma function (2.75)
    real (kind = r_grid), parameter :: gam290 = 1.8273550806240362 ! Gamma function (2.9)
    real (kind = r_grid), parameter :: gam325 = 2.549256966718529 ! Gamma function (3.25)
    real (kind = r_grid), parameter :: gam350 = 3.323350970447843 ! Gamma function (3.5)
    real (kind = r_grid), parameter :: gam380 = 4.694174205740422 ! Gamma function (3.8)
    real (kind = r_grid), parameter :: gam425 = 8.28508514183522 ! Gamma function (4.25)
    real (kind = r_grid), parameter :: gam450 = 11.63172839656745 ! Gamma function (4.5)
    real (kind = r_grid), parameter :: gam480 = 17.837861981813603 ! Gamma function (4.8)
    real (kind = r_grid), parameter :: gam625 = 184.86096222719834 ! Gamma functio (6.25)
    real (kind = r_grid), parameter :: gam680 = 496.60607757369075 ! Gamma functio (6.8)
    
    real, parameter :: visk = 1.259e-5 ! kinematic viscosity of air (cm^2/s)
    real, parameter :: vdifu = 2.11e-5 ! diffusivity of water vapor in air (cm^2/s)
    real, parameter :: tcond = 2.36e-2 ! thermal conductivity of air (J/m/s/K)

    real (kind = r_grid), parameter :: d2ice = cp_vap - c_ice ! - 260.0, isobaric heating / cooling (J/kg/K)

    real (kind = r_grid), parameter :: li2 = lv0 + li0 ! 2.9220216e6, sublimation latent heat coeff. at 0 deg K (J/kg)
    
    real (kind = r_grid), parameter :: e00 = 611.21 ! saturation vapor pressure at 0 deg C (Pa)
    
    ! -----------------------------------------------------------------------
    ! predefined parameters
    ! -----------------------------------------------------------------------
    
    real, parameter :: qcmin = 1.0e-12 ! min value for cloud condensates (kg/kg)
    real, parameter :: qfmin = 1.0e-8 ! min value for sedimentation (kg/kg)
    
    real, parameter :: dz_min = 1.0e-2 ! used for correcting flipped height (m)
    
    real, parameter :: sfcrho = 1.2 ! surface air density (kg/m^3)
    
    real, parameter :: rnzr = 8.0e6 ! intercept parameter of rain (Lin et al. 1983) (1/m^4)
    real, parameter :: rnzs = 3.0e6 ! intercept parameter of snow (Lin et al. 1983) (1/m^4)
    real, parameter :: rnzg = 4.0e6 ! intercept parameter of graupel (Rutledge and Hobbs 1984) (1/m^4)
    real, parameter :: rnzh = 4.0e4 ! intercept parameter of hail (Lin et al. 1983) (1/m^4)
    
    real, parameter :: rhow = 1.0e3 ! density of cloud water (kg/m^3)
    real, parameter :: rhor = 1.0e3 ! density of rain (Lin et al. 1983) (kg/m^3)
    real, parameter :: rhos = 0.1e3 ! density of snow (Lin et al. 1983) (kg/m^3)
    real, parameter :: rhog = 0.4e3 ! density of graupel (Rutledge and Hobbs 1984) (kg/m^3)
    real, parameter :: rhoh = 9.17e2 ! density of hail (Lin et al. 1983) (kg/m^3)
    
    real, parameter :: dt_fr = 8.0 ! t_wfr - dt_fr: minimum temperature water can exist (Moore and Molinero 2011)
    
    real, parameter :: p0_min = 100.0 ! minimum pressure for mp to operate (Pa)
    
    real, parameter :: acc (3) = (/ 5.0, 2.0, 0.5 /) ! accretion coefficients
    
    real (kind = r_grid), parameter :: one_r8 = 1.0 ! constant 1
    
    ! -----------------------------------------------------------------------
    ! namelist parameters
    ! -----------------------------------------------------------------------
    
    integer :: ntimes = 1 ! cloud microphysics sub cycles
    
    integer :: cfflag = 1 ! cloud fraction scheme
    ! 1: GFDL cloud scheme
    ! 2: Xu and Randall (1996)
    ! 3: Park et al. (2016)
    ! 4: Pultepe and Isaac (2007)
    
    integer :: icloud_f = 0 ! GFDL cloud scheme
    ! 0: subgrid variability based scheme
    ! 1: same as 0, but for old fvgfs implementation
    ! 2: binary cloud scheme
    ! 3: extension of  0

    integer :: irain_f = 0 ! cloud water to rain auto conversion scheme
    ! 0: subgrid variability based scheme
    ! 1: no subgrid varaibility
    
    integer :: inflag = 1 ! ice nucleation scheme
    ! 1: Hong et al. (2004)
    ! 2: Meyers et al. (1992)
    ! 3: Meyers et al. (1992)
    ! 4: Cooper (1986)
    ! 5: Flecther (1962)

    integer :: igflag = 3 ! ice generation scheme
    ! 1: WSM6
    ! 2: WSM6 with 0 at 0 C
    ! 3: WSM6 with 0 at 0 C and fixed value at - 10 C
    ! 4: combination of 1 and 3

    integer :: ifflag = 1 ! ice fall scheme
    ! 1: Deng and Mace (2008)
    ! 2: Heymsfield and Donner (1990)

    integer :: rewflag = 1 ! cloud water effective radius scheme
    ! 1: Martin et al. (1994)
    ! 2: Martin et al. (1994), GFDL revision
    ! 3: Kiehl et al. (1994)

    integer :: reiflag = 1 ! cloud ice effective radius scheme
    ! 1: Heymsfield and Mcfarquhar (1996)
    ! 2: Donner et al. (1997)
    ! 3: Fu (2007)
    ! 4: Kristjansson et al. (2000)
    ! 5: Wyser (1998)
    
    logical :: do_sedi_uv = .true. ! transport of horizontal momentum in sedimentation
    logical :: do_sedi_w = .true. ! transport of vertical momentum in sedimentation
    logical :: do_sedi_heat = .true. ! transport of heat in sedimentation
    logical :: do_disp_heat = .true. ! dissipative heating due to sedimentation

    logical :: do_qa = .true. ! do inline cloud fraction
    logical :: rad_snow = .true. ! include snow in cloud fraciton calculation
    logical :: rad_graupel = .true. ! include graupel in cloud fraction calculation
    logical :: rad_rain = .true. ! include rain in cloud fraction calculation
    logical :: do_cld_adj = .false. ! do cloud fraction adjustment

    logical :: use_ppm = .false. ! use ppm fall scheme
    logical :: mono_prof = .true. ! perform terminal fall with mono ppm scheme

    logical :: z_slope_liq = .true. ! use linear mono slope for autocconversions
    logical :: z_slope_ice = .true. ! use linear mono slope for autocconversions

    logical :: use_rhc_cevap = .false. ! cap of rh for cloud water evaporation
    logical :: use_rhc_revap = .false. ! cap of rh for rain evaporation

    logical :: const_vi = .false. ! if .ture., the constants are specified by v * _fac
    logical :: const_vs = .false. ! if .ture., the constants are specified by v * _fac
    logical :: const_vg = .false. ! if .ture., the constants are specified by v * _fac
    logical :: const_vr = .false. ! if .ture., the constants are specified by v * _fac
    
    logical :: liq_ice_combine = .false. ! combine all liquid water, combine all solid water
    logical :: snow_grauple_combine = .true. ! combine snow and graupel
    
    logical :: prog_ccn = .false. ! do prognostic ccn (Yi Ming's method)

    logical :: fix_negative = .false. ! fix negative water species

    logical :: do_cond_timescale = .false. ! whether to apply a timescale to condensation

    logical :: do_sat_adj = .false. ! do fast saturation adjustments

    logical :: do_hail = .false. ! use hail parameters instead of graupel

    logical :: consv_checker = .false. ! turn on energy and water conservation checker

    logical :: do_warm_rain_mp = .false. ! do warm rain cloud microphysics only

    real :: mp_time = 150.0 ! maximum microphysics time step (s)
    
    real :: tice_mlt = 273.16 ! can set ice melting temperature to 268 based on observation (Kay et al. 2016) (K)
    
    real :: t_min = 178.0 ! minimum temperature to freeze - dry all water vapor (K)
    real :: t_sub = 184.0 ! minimum temperature for sublimation of cloud ice (K)

    real :: rh_inc = 0.25 ! rh increment for complete evaporation of cloud water and cloud ice
    real :: rh_inr = 0.25 ! rh increment for minimum evaporation of rain
    real :: rh_ins = 0.25 ! rh increment for sublimation of snow
    
    real :: tau_r2g = 900.0 ! rain freezing to graupel time scale (s)
    real :: tau_g2r = 600.0 ! graupel melting to rain time scale (s)
    real :: tau_i2s = 1000.0 ! cloud ice to snow autoconversion time scale (s)
    real :: tau_l2r = 900.0 ! cloud water to rain autoconversion time scale (s)
    real :: tau_v2l = 150.0 ! water vapor to cloud water condensation time scale (s)
    real :: tau_l2v = 300.0 ! cloud water to water vapor evaporation time scale (s)
    real :: tau_revp = 0.0 ! rain evaporation time scale (s)
    real :: tau_imlt = 1200.0 ! cloud ice melting time scale (s)
    real :: tau_smlt = 900.0 ! snow melting time scale (s)
    
    real :: dw_land = 0.20 ! base value for subgrid deviation / variability over land
    real :: dw_ocean = 0.10 ! base value for subgrid deviation / variability over ocean
    
    real :: ccn_o = 90.0 ! ccn over ocean (1/cm^3)
    real :: ccn_l = 270.0 ! ccn over land (1/cm^3)
    
    real :: rthresh = 10.0e-6 ! critical cloud drop radius (micron) for autoconversion
    
    real :: cld_min = 0.05 ! minimum cloud fraction

    real :: sat_adj0 = 0.90 ! adjustment factor (0: no, 1: full) during fast_sat_adj
    
    real :: qi_lim = 1.0 ! cloud ice limiter (0: no, 1: full, >1: extra) to prevent large ice build up
    
    real :: ql_mlt = 2.0e-3 ! maximum cloud water allowed from melted cloud ice (kg/kg)
    real :: qs_mlt = 1.0e-6 ! maximum cloud water allowed from melted snow (kg/kg)
    
    real :: ql_gen = 1.0e-3 ! maximum cloud water generation during remapping step if do_sat_adj = .true. (kg/kg)
    
    real :: ql0_max = 2.0e-3 ! maximum cloud water value (autoconverted to rain) (kg/kg)
    real :: qi0_max = 1.0e-4 ! maximum cloud ice value (by other sources) (kg/m^3)
    
    real :: qi0_crt = 1.0e-4 ! cloud ice to snow autoconversion threshold (kg/m^3)
    real :: qs0_crt = 1.0e-3 ! snow to graupel autoconversion threshold (0.6e-3 in Purdue Lin scheme) (kg/m^3)
    
    real :: c_paut = 0.55 ! cloud water to rain autoconversion efficiency
    real :: c_psaci = 0.02 ! cloud ice to snow accretion efficiency (was 0.1 in ZETAC)
    real :: c_pracw = 0.9 ! cloud water to rain accretion efficiency
    real :: c_pgacs = 2.0e-3 ! snow to graupel accretion efficiency (was 0.1 in ZETAC)
    real :: c_pgaci = 0.05 ! cloud ice to graupel accretion efficiency (was 0.1 in ZETAC)
    
    real :: alin = 842.0 ! "a" in Lin et al. (1983)
    real :: clin = 4.8 ! "c" in Lin et al. (1983), 4.8 -- > 6. (to ehance ql -- > qs)
    
    real :: vi_fac = 1.0 ! IFS: if const_vi: 1 / 3
    real :: vs_fac = 1.0 ! IFS: if const_vs: 1.
    real :: vg_fac = 1.0 ! IFS: if const_vg: 2.
    real :: vr_fac = 1.0 ! IFS: if const_vr: 4.
    
    real :: vi_max = 0.5 ! maximum fall speed for cloud ice (m/s)
    real :: vs_max = 5.0 ! maximum fall speed for snow (m/s)
    real :: vg_max = 8.0 ! maximum fall speed for graupel (m/s)
    real :: vr_max = 12.0 ! maximum fall speed for rain (m/s)
    
    real :: xr_a = 0.25 ! p value in Xu and Randall (1996)
    real :: xr_b = 100.0 ! alpha_0 value in Xu and Randall (1996)
    real :: xr_c = 0.49 ! gamma value in Xu and Randall (1996)
    
    real :: te_err = 1.e-14 ! 64bit: 1.e-14, 32bit: 1.e-7; turn off to save computer time
    
    real :: rh_thres = 0.75 ! minimum relative humidity for cloud fraction
    real :: rhc_cevap = 0.85 ! maximum relative humidity for cloud water evaporation
    real :: rhc_revap = 0.85 ! maximum relative humidity for rain evaporation

    real :: f_dq_p = 1.0 ! cloud fraction adjustment for supersaturation
    real :: f_dq_m = 1.0 ! cloud fraction adjustment for undersaturation

    real :: fi2s_fac = 0.75 ! maximum sink of cloud ice to form snow: 0-1
    real :: fs2g_fac = 1.0 ! maximum sink of snow to form graupel: 0-1
    
    real :: qi0_rei = 0.8e-4 ! maximum cloud ice value (by other sources) (kg/kg)
    
    real :: beta = 1.22 ! defined in Heymsfield and Mcfarquhar (1996)
    
    real :: rewmin = 5.0, rewmax = 15.0 ! minimum and maximum effective radius for cloud water (micron)
    real :: reimin = 10.0, reimax = 150.0 ! minimum and maximum effective radius for cloud ice (micron)
    real :: rermin = 15.0, rermax = 10000.0 ! minimum and maximum effective radius for rain (micron)
    real :: resmin = 150.0, resmax = 10000.0 ! minimum and maximum effective radius for snow (micron)
    real :: regmin = 0.0, regmax = 10000.0 ! minimum and maximum effective radius for graupel
    !real :: rewmax = 15.0, rermin = 15.0 ! Kokhanovsky (2004)
    
    ! -----------------------------------------------------------------------
    ! local shared variables
    ! -----------------------------------------------------------------------
    
    real :: acco (3, 4)
    real :: cracs, csacr, cgacr, cgacs, csacw, craci, csaci, cgacw, cgaci, cracw
    real :: cssub (5), cgsub (5), crevp (5), cgfr (2), csmlt (5), cgmlt (5)
    
    real :: t_wfr, p_min, fac_rc, c_air, c_vap, d0_vap

    real (kind = r_grid) :: lv00, li00, li20
    real (kind = r_grid) :: d1_vap, d1_ice, c1_vap, c1_liq, c1_ice
    real (kind = r_grid) :: vconr, vcons, vcong, vconh, normr, norms, normg, normh
    
    real, allocatable :: table1 (:), table2 (:), table3 (:), table4 (:)
    real, allocatable :: des1 (:), des2 (:), des3 (:), des4 (:)
    
    ! -----------------------------------------------------------------------
    ! namelist
    ! -----------------------------------------------------------------------
    
    namelist / gfdl_mp_nml / &
        t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, dw_land, dw_ocean, vi_fac, &
        vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vi_max, vs_max, &
        vg_max, vr_max, qs_mlt, qs0_crt, ql0_max, qi0_max, qi0_crt, ifflag, &
        do_sat_adj, rh_inc, rh_ins, rh_inr, const_vi, const_vs, const_vg, &
        const_vr, rthresh, ccn_l, ccn_o, sat_adj0, igflag, c_paut, tau_imlt, &
        tau_v2l, tau_l2v, tau_i2s, tau_l2r, qi_lim, ql_gen, do_hail, inflag, &
        c_psaci, c_pgacs, c_pgaci, z_slope_liq, z_slope_ice, prog_ccn, &
        c_pracw, alin, clin, rad_snow, rad_graupel, rad_rain, cld_min, &
        use_ppm, mono_prof, do_sedi_heat, do_sedi_uv, do_sedi_w, icloud_f, &
        irain_f, xr_a, xr_b, xr_c, ntimes, do_disp_heat, tau_revp, tice_mlt, &
        do_cond_timescale, mp_time, consv_checker, te_err, use_rhc_cevap, &
        use_rhc_revap, do_warm_rain_mp, rh_thres, f_dq_p, f_dq_m, do_cld_adj, &
        rhc_cevap, rhc_revap, qi0_rei, beta, liq_ice_combine, rewflag, &
        reiflag, rewmin, rewmax, reimin, reimax, rermin, rermax, resmin, resmax, &
        regmin, regmax, fs2g_fac, fi2s_fac
    
contains

! =======================================================================
! GFDL cloud microphysics initialization
! =======================================================================

subroutine gfdl_cld_mp_init (me, master, nlunit, input_nml_file, logunit, fn_nml)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: me, master, nlunit, logunit
    
    character (len = 64), intent (in) :: fn_nml

    character (len = *), intent (in) :: input_nml_file (:)
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    logical :: exists
    
    ! -----------------------------------------------------------------------
    ! read namelist
    ! -----------------------------------------------------------------------
    
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml = gfdl_mp_nml)
#else
    inquire (file = trim (fn_nml), exist = exists)
    if (.not. exists) then
        write (6, *) 'gfdl - mp :: namelist file: ', trim (fn_nml), ' does not exist'
        stop
    else
        open (unit = nlunit, file = fn_nml, readonly, status = 'old')
    endif
    rewind (nlunit)
    read (nlunit, nml = gfdl_mp_nml)
    close (nlunit)
#endif
    
    ! -----------------------------------------------------------------------
    ! write namelist to log file
    ! -----------------------------------------------------------------------
    
    if (me .eq. master) then
        write (logunit, *) " ================================================================== "
        write (logunit, *) "gfdl_mp_mod"
        write (logunit, nml = gfdl_mp_nml)
    endif
    
    ! -----------------------------------------------------------------------
    ! initialize microphysics variables
    ! -----------------------------------------------------------------------
    
    if (.not. tables_are_initialized) call qsmith_init
    
    call setup_mp
    
end subroutine gfdl_cld_mp_init

! =======================================================================
! GFDL cloud microphysics driver
! =======================================================================

subroutine gfdl_cld_mp_driver (qv, ql, qr, qi, qs, qg, qa, qnl, qni, pt, w, &
        ua, va, dz, delp, gsize, dts, hs, rain, snow, ice, graupel, &
        hydrostatic, is, ie, ks, ke, q_con, cappa, consv_te, te, &
        condensation, deposition, evaporation, sublimation, last_step, &
        do_inline_mp)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: is, ie, ks, ke
    
    logical, intent (in) :: hydrostatic, last_step, consv_te, do_inline_mp
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (is:ie) :: hs, gsize
    
    real, intent (in), dimension (is:ie, ks:ke) :: dz, qnl, qni
    
    real, intent (inout), dimension (is:ie, ks:ke) :: delp, pt, ua, va, w, te
    real, intent (inout), dimension (is:ie, ks:ke) :: qv, ql, qr, qi, qs, qg, qa

    real, intent (inout), dimension (is:, ks:) :: q_con, cappa

    real, intent (inout), dimension (is:ie) :: rain, snow, ice, graupel
    real, intent (inout), dimension (is:ie) :: condensation, deposition
    real, intent (inout), dimension (is:ie) :: evaporation, sublimation

    ! -----------------------------------------------------------------------
    ! top level for microphysical processes
    ! -----------------------------------------------------------------------
    
    if (last_step) then
        p_min = p0_min ! final cleanup
    else
        p_min = 30.e2 ! time saving trick
    endif
    
    ! -----------------------------------------------------------------------
    ! define various heat capacities and latent heat coefficients at 0 deg K
    ! -----------------------------------------------------------------------
    
    if (hydrostatic) then
        c_air = cp_air
        c_vap = cp_vap
        do_sedi_w = .false.
    else
        c_air = cv_air
        c_vap = cv_vap
    endif
    d0_vap = c_vap - c_liq
    
    ! scaled constants (to reduce float point errors for 32-bit)

    d1_vap = d0_vap / c_air
    d1_ice = dc_ice / c_air
    
    lv00 = (hlv0 - d0_vap * tice) / c_air
    li00 = (hlf0 - dc_ice * tice) / c_air
    li20 = lv00 + li00
    
    c1_vap = c_vap / c_air
    c1_liq = c_liq / c_air
    c1_ice = c_ice / c_air
    
    ! -----------------------------------------------------------------------
    ! major cloud microphysics driver
    ! -----------------------------------------------------------------------
    
    call mpdrv (hydrostatic, ua, va, w, delp, pt, qv, ql, qr, qi, qs, qg, qa, &
        qnl, qni, dz, is, ie, ks, ke, dts, rain, snow, graupel, ice, gsize, &
        hs, q_con, cappa, consv_te, te, condensation, deposition, evaporation, &
        sublimation, last_step, do_inline_mp)
    
end subroutine gfdl_cld_mp_driver

! =======================================================================
! GFDL cloud microphysics end
! =======================================================================

subroutine gfdl_cld_mp_end
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! free up memory
    ! -----------------------------------------------------------------------
    
    deallocate (table1)
    deallocate (table2)
    deallocate (table3)
    deallocate (table4)
    deallocate (des1)
    deallocate (des2)
    deallocate (des3)
    deallocate (des4)
    
    tables_are_initialized = .false.
    
end subroutine gfdl_cld_mp_end

! =======================================================================
! setup GFDL cloud microphysics parameters
! =======================================================================

subroutine setup_mp
    
    implicit none
    
    integer :: i, k
    
    real :: gcon, scm3, pisq, act (8)
    
    ! -----------------------------------------------------------------------
    ! complete freezing temperature
    ! -----------------------------------------------------------------------
    
    if (do_warm_rain_mp) then
        t_wfr = t_min
    else
        t_wfr = tice - 40.0
    endif ! do_warm_rain_mp
    
    ! -----------------------------------------------------------------------
    ! cloud water autoconversion threshold in mass
    ! -----------------------------------------------------------------------
    
    fac_rc = (4. / 3.) * pi * rhor * rthresh ** 3
    
    ! -----------------------------------------------------------------------
    ! terminal fall speed for rain, snow, and graupel or hail
    ! -----------------------------------------------------------------------
    
    gcon = 40.74 * sqrt (sfcrho)
    
    vconr = alin * gam480 / 6.0
    vcons = clin * gam425 / 6.0
    vcong = gam450 / 6.0 * gcon
    vconh = vcong * sqrt (rhoh / rhog)

    normr = pi * rhor * rnzr
    norms = pi * rhos * rnzs
    normg = pi * rhog * rnzg
    normh = pi * rhoh * rnzh

    ! -----------------------------------------------------------------------
    ! physics constants
    ! -----------------------------------------------------------------------
    
    scm3 = (visk / vdifu) ** (1. / 3.)
    
    ! -----------------------------------------------------------------------
    ! accretion between rain, snow, and graupel or hail
    ! -----------------------------------------------------------------------
    
    pisq = pi * pi
    
    cracs = pisq * rnzr * rnzs * rhos
    csacr = pisq * rnzr * rnzs * rhor
    if (do_hail) then
        cgacr = pisq * rnzr * rnzh * rhor
        cgacs = pisq * rnzh * rnzs * rhos
    else
        cgacr = pisq * rnzr * rnzg * rhor
        cgacs = pisq * rnzg * rnzs * rhos
    endif
    cgacs = cgacs * c_pgacs
    
    ! act:
    ! 1-2: racs (s-r)
    ! 3-4: sacr (r-s)
    ! 5-6: gacr (r-g)
    ! 7-8: gacs (s-g)
    
    act (1) = norms
    act (2) = normr
    if (do_hail) then
        act (6) = normh
    else
        act (6) = normg
    endif
    act (3) = act (2)
    act (4) = act (1)
    act (5) = act (2)
    act (7) = act (1)
    act (8) = act (6)
    
    do i = 1, 3
        do k = 1, 4
            acco (i, k) = acc (i) / (act (2 * k - 1) ** ((7 - i) * 0.25) * act (2 * k) ** (i * 0.25))
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! accretion between cloud water, cloud ice, rain, and snow
    ! -----------------------------------------------------------------------
    
    csacw = pi * rnzs * clin * gam325 / (4. * act (1) ** 0.8125)
    craci = pi * rnzr * alin * gam380 / (4. * act (2) ** 0.95)
    csaci = csacw * c_psaci
    cracw = craci * c_pracw
    
    ! -----------------------------------------------------------------------
    ! accretion between cloud water, cloud ice, and graupel or hail
    ! -----------------------------------------------------------------------
    
    if (do_hail) then
        cgacw = pi * rnzh * gam350 * gcon / (4. * act (6) ** 0.875)
    else
        cgacw = pi * rnzg * gam350 * gcon / (4. * act (6) ** 0.875)
    endif
    cgaci = cgacw * c_pgaci
    
    ! -----------------------------------------------------------------------
    ! snow sublimation
    ! -----------------------------------------------------------------------

    cssub (1) = 2. * pi * vdifu * tcond * rvgas * rnzs
    cssub (2) = 0.78 / sqrt (act (1))
    cssub (3) = 0.31 * scm3 * gam263 * sqrt (clin / visk) / act (1) ** 0.65625
    cssub (4) = tcond * rvgas
    cssub (5) = (hlv + hlf) ** 2 * vdifu

    ! -----------------------------------------------------------------------
    ! graupel or hail sublimation
    ! -----------------------------------------------------------------------
    
    if (do_hail) then
        cgsub (1) = 2. * pi * vdifu * tcond * rvgas * rnzh
    else
        cgsub (1) = 2. * pi * vdifu * tcond * rvgas * rnzg
    endif
    cgsub (2) = 0.78 / sqrt (act (6))
    cgsub (3) = 0.31 * scm3 * gam275 * sqrt (gcon / visk) / act (6) ** 0.6875
    cgsub (4) = cssub (4)
    cgsub (5) = cssub (5)

    ! -----------------------------------------------------------------------
    ! rain evaporation
    ! -----------------------------------------------------------------------
    
    crevp (1) = 2. * pi * vdifu * tcond * rvgas * rnzr
    crevp (2) = 0.78 / sqrt (act (2))
    crevp (3) = 0.31 * scm3 * gam290 * sqrt (alin / visk) / act (2) ** 0.725
    crevp (4) = cssub (4)
    crevp (5) = hlv ** 2 * vdifu
    
    ! -----------------------------------------------------------------------
    ! rain freezing
    ! -----------------------------------------------------------------------
    
    cgfr (1) = 20.e2 * pisq * rnzr * rhor / act (2) ** 1.75
    cgfr (2) = 0.66
    
    ! -----------------------------------------------------------------------
    ! snow melting
    ! -----------------------------------------------------------------------
    
    csmlt (1) = 2. * pi * tcond * rnzs / hlf
    csmlt (2) = 2. * pi * vdifu * rnzs * hlv / hlf
    csmlt (3) = cssub (2)
    csmlt (4) = cssub (3)
    csmlt (5) = c_liq / hlf
    
    ! -----------------------------------------------------------------------
    ! graupel or hail melting
    ! -----------------------------------------------------------------------
    
    if (do_hail) then
        cgmlt (1) = 2. * pi * tcond * rnzh / hlf
        cgmlt (2) = 2. * pi * vdifu * rnzh * hlv / hlf
    else
        cgmlt (1) = 2. * pi * tcond * rnzg / hlf
        cgmlt (2) = 2. * pi * vdifu * rnzg * hlv / hlf
    endif
    cgmlt (3) = cgsub (2)
    cgmlt (4) = cgsub (3)
    cgmlt (5) = c_liq / hlf
    
end subroutine setup_mp

! =======================================================================
! major cloud microphysics driver
! =======================================================================

subroutine mpdrv (hydrostatic, ua, va, w, delp, pt, qv, ql, qr, qi, qs, &
        qg, qa, qnl, qni, dz, is, ie, ks, ke, dt_in, rain, snow, graupel, &
        ice, gsize, hs, q_con, cappa, consv_te, te, condensation, &
        deposition, evaporation, sublimation, last_step, do_inline_mp)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: is, ie, ks, ke

    logical, intent (in) :: hydrostatic, last_step, consv_te, do_inline_mp

    real, intent (in) :: dt_in

    real, intent (in), dimension (is:ie) :: gsize, hs

    real, intent (in), dimension (is:ie, ks:ke) :: dz, qnl, qni
    
    real, intent (inout), dimension (is:ie, ks:ke) :: delp, pt, ua, va, w
    real, intent (inout), dimension (is:ie, ks:ke) :: qv, ql, qr, qi, qs, qg, qa

    real, intent (inout), dimension (is:, ks:) :: q_con, cappa

    real, intent (inout), dimension (is:ie) :: rain, snow, ice, graupel
    real, intent (inout), dimension (is:ie) :: condensation, deposition
    real, intent (inout), dimension (is:ie) :: evaporation, sublimation
    
    real, intent (out), dimension (is:ie, ks:ke) :: te

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i, k, n
    
    real :: cpaut, rh_adj, rh_rain, r1, s1, i1, g1, rdt, ccn0, cond, dep, reevap, sub
    real :: dt_rain, convt, dts, q_cond, s_leng, t_land, t_ocean, h_var, tmp
    
    real, dimension (ks:ke) :: q_liq, q_sol, vtiz, vtsz, vtgz, vtrz, dp1, dz1
    real, dimension (ks:ke) :: qvz, qlz, qrz, qiz, qsz, qgz, qaz, c_praut, m1
    real, dimension (ks:ke) :: den, p1, denfac, ccn, cin, m1_rain, m1_sol
    real, dimension (ks:ke) :: u0, v0, u1, v1, w1

    real, dimension (is:ie, ks:ke) :: m2_rain, m2_sol
    
    real (kind = r_grid) :: con_r8, c8

    real (kind = r_grid), dimension (is:ie, ks:ke) :: te_beg, te_end, tw_beg, tw_end
    real (kind = r_grid), dimension (is:ie, ks:ke) :: te_beg_0, te_end_0, tw_beg_0, tw_end_0

    real (kind = r_grid), dimension (is:ie) :: te_b_beg, te_b_end, tw_b_beg, tw_b_end, dte, te_loss
    real (kind = r_grid), dimension (is:ie) :: te_b_beg_0, te_b_end_0, tw_b_beg_0, tw_b_end_0

    real (kind = r_grid), dimension (ks:ke) :: dp0, tz, cvm, te1, te2
    
    ! -----------------------------------------------------------------------
    ! time steps
    ! -----------------------------------------------------------------------

    ntimes = max (ntimes, int (dt_in / min (dt_in, mp_time)))
    dts = dt_in / real (ntimes)
    
    dt_rain = dts * 0.5
    rdt = one_r8 / dts
    
    ! -----------------------------------------------------------------------
    ! initialization of total energy difference
    ! -----------------------------------------------------------------------

    dte = 0.0
    
    ! -----------------------------------------------------------------------
    ! unit convert to mm/day
    ! -----------------------------------------------------------------------

    convt = 86400. * rdt * rgrav
    
    do i = is, ie
        
        ! -----------------------------------------------------------------------
        ! conversion of temperature
        ! -----------------------------------------------------------------------
    
        if (do_inline_mp) then
            do k = ks, ke
#ifdef MOIST_CAPPA
                tz (k) = pt (i, k) / ((1. + zvir * qv (i, k)) * &
                   (1. - (ql (i, k) + qr (i, k) + qi (i, k) + qs (i, k) + qg (i, k))))
#else
                tz (k) = pt (i, k) / (1. + zvir * qv (i, k))
#endif
            enddo
        else
            do k = ks, ke
                tz (k) = pt (i, k)
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do k = ks, ke
                q_liq (k) = ql (i, k) + qr (i, k)
                q_sol (k) = qi (i, k) + qs (i, k) + qg (i, k)
                cvm (k) = c_air * (1.0 - qv (i, k) - q_liq (k) - q_sol (k)) + &
                    qv (i, k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                te_beg_0 (i, k) = cvm (k) * tz (k) + lv00 * c_air * qv (i, k) - li00 * c_air * q_sol (k)
                if (hydrostatic) then
                    te_beg_0 (i, k) = te_beg_0 (i, k) + 0.5 * (ua (i, k) ** 2 + va (i, k) ** 2)
                else
                    te_beg_0 (i, k) = te_beg_0 (i, k) + 0.5 * (ua (i, k) ** 2 + va (i, k) ** 2 + w (i, k) ** 2)
                endif
                te_beg_0 (i, k) = rgrav * te_beg_0 (i, k) * delp (i, k) * gsize (i) ** 2.0
                tw_beg_0 (i, k) = rgrav * (qv (i, k) + q_liq (k) + q_sol (k)) * delp (i, k) * gsize (i) ** 2.0
            enddo
            te_b_beg_0 (i) = (dte (i) - li00 * c_air * (ice (i) + snow (i) + graupel (i)) * dt_in / 86400) * &
                gsize (i) ** 2.0
            tw_b_beg_0 (i) = (rain (i) + ice (i) + snow (i) + graupel (i)) * dt_in / 86400 * gsize (i) ** 2.0
        endif
        
        do k = ks, ke

            ! -----------------------------------------------------------------------
            ! initialization
            ! -----------------------------------------------------------------------

            m1 (k) = 0.
            qaz (k) = 0.
            dp0 (k) = delp (i, k)
            dz1 (k) = dz (i, k)

            ! -----------------------------------------------------------------------
            ! convert moist mixing ratios to dry mixing ratios
            ! -----------------------------------------------------------------------

            qvz (k) = qv (i, k)
            qlz (k) = ql (i, k)
            qrz (k) = qr (i, k)
            qiz (k) = qi (i, k)
            qsz (k) = qs (i, k)
            qgz (k) = qg (i, k)

            q_liq (k) = qlz (k) + qrz (k)
            q_sol (k) = qiz (k) + qsz (k) + qgz (k)
            q_cond = q_liq (k) + q_sol (k)
            con_r8 = one_r8 - (qvz (k) + q_cond)

            dp1 (k) = dp0 (k) * con_r8
            con_r8 = one_r8 / con_r8
            qvz (k) = qvz (k) * con_r8
            qlz (k) = qlz (k) * con_r8
            qrz (k) = qrz (k) * con_r8
            qiz (k) = qiz (k) * con_r8
            qsz (k) = qsz (k) * con_r8
            qgz (k) = qgz (k) * con_r8
            
            ! -----------------------------------------------------------------------
            ! dry air density and layer-mean pressure thickness
            ! -----------------------------------------------------------------------
            
            den (k) = - dp1 (k) / (grav * dz1 (k))
            p1 (k) = den (k) * rdgas * tz (k)
            denfac (k) = sqrt (sfcrho / den (k))
            
            ! -----------------------------------------------------------------------
            ! for sedi_momentum transport
            ! -----------------------------------------------------------------------
            
            u0 (k) = ua (i, k)
            v0 (k) = va (i, k)
            if (.not. hydrostatic) then
                w1 (k) = w (i, k)
            endif
            u1 (k) = u0 (k)
            v1 (k) = v0 (k)

        enddo ! k loop
        
        ! -----------------------------------------------------------------------
        ! calculate total energy
        ! -----------------------------------------------------------------------
        
        if (consv_te) then
            if (hydrostatic) then
                do k = ks, ke
                    te (i, k) = - c_air * tz (k) * delp (i, k)
                enddo
            else
                do k = ks, ke
#ifdef MOIST_CAPPA
                    q_liq (k) = ql (i, k) + qr (i, k)
                    q_sol (k) = qi (i, k) + qs (i, k) + qg (i, k)
                    q_cond = q_liq (k) + q_sol (k)
                    cvm (k) = (one_r8 - (qv (i, k) + q_cond)) * c_air + &
                        qv (i, k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                    te (i, k) = - cvm (k) * tz (k) * delp (i, k)
#else
                    te (i, k) = - c_air * tz (k) * delp (i, k)
#endif
                enddo
            endif
        endif
        
        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do k = ks, ke
                q_liq (k) = qlz (k) + qrz (k)
                q_sol (k) = qiz (k) + qsz (k) + qgz (k)
                cvm (k) = c_air + qvz (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                te_beg (i, k) = cvm (k) * tz (k) + lv00 * c_air * qvz (k) - li00 * c_air * q_sol (k)
                if (hydrostatic) then
                    te_beg (i, k) = te_beg (i, k) + 0.5 * (u1 (k) ** 2 + v1 (k) ** 2)
                else
                    te_beg (i, k) = te_beg (i, k) + 0.5 * (u1 (k) ** 2 + v1 (k) ** 2 + w1 (k) ** 2)
                endif
                te_beg (i, k) = rgrav * te_beg (i, k) * dp1 (k) * gsize (i) ** 2.0
                tw_beg (i, k) = rgrav * (qvz (k) + q_liq (k) + q_sol (k)) * dp1 (k) * gsize (i) ** 2.0
            enddo
            te_b_beg (i) = (dte (i) - li00 * c_air * (ice (i) + snow (i) + graupel (i)) * dt_in / 86400) * &
                gsize (i) ** 2.0
            tw_b_beg (i) = (rain (i) + ice (i) + snow (i) + graupel (i)) * dt_in / 86400 * gsize (i) ** 2.0
        endif
        
        ! -----------------------------------------------------------------------
        ! calculate cloud droplet concentration based on cloud condensation nuclei (CCN)
        ! -----------------------------------------------------------------------
        
        cpaut = c_paut * 0.104 * grav / 1.717e-5
        if (prog_ccn) then
            do k = ks, ke
                ccn (k) = max (10.0, qnl (i, k)) * 1.e6
                cin (k) = max (10.0, qni (i, k)) * 1.e6
                ccn (k) = ccn (k) / den (k)
                c_praut (k) = cpaut * (ccn (k) * rhor) ** (- 1. / 3.)
            enddo
        else
            ccn0 = (ccn_l * min (1., abs (hs (i)) / (10. * grav)) + &
                ccn_o * (1. - min (1., abs (hs (i)) / (10. * grav)))) * 1.e6
            do k = ks, ke
                ccn (k) = ccn0 / den (k)
                c_praut (k) = cpaut * (ccn (k) * rhor) ** (- 1. / 3.)
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! calculate horizontal subgrid variability
        ! subgrid deviation in horizontal direction
        ! default area dependent form: use dx ~ 100 km as the base
        ! -----------------------------------------------------------------------
        
        s_leng = sqrt (gsize (i) / 1.e5)
        t_land = dw_land * s_leng
        t_ocean = dw_ocean * s_leng
        tmp = min (1., abs (hs (i)) / (10. * grav))
        h_var = t_land * tmp + t_ocean * (1. - tmp)
        h_var = min (0.20, max (0.01, h_var))
        
        ! -----------------------------------------------------------------------
        ! relative humidity thresholds
        ! -----------------------------------------------------------------------
        
        rh_adj = 1. - h_var - rh_inc
        rh_rain = max (0.35, rh_adj - rh_inr)
        
        ! -----------------------------------------------------------------------
        ! fix negative water species
        ! -----------------------------------------------------------------------
        
        if (fix_negative) &
            call neg_adj (ks, ke, tz, dp1, qvz, qlz, qrz, qiz, qsz, qgz, cond)
        
        condensation (i) = condensation (i) + cond * convt * ntimes
            
        ! -----------------------------------------------------------------------
        ! initialization
        ! -----------------------------------------------------------------------

        m2_rain (i, :) = 0.
        m2_sol (i, :) = 0.
        
        do n = 1, ntimes
            
            ! -----------------------------------------------------------------------
            ! time-split warm rain processes: 1st pass
            ! -----------------------------------------------------------------------
            
            call warm_rain (dt_rain, ks, ke, dp1, dz1, tz, qvz, qlz, qrz, qiz, qsz, &
                qgz, den, denfac, ccn, c_praut, rh_rain, vtrz, r1, m1_rain, w1, h_var, reevap, dte (i))
            
            evaporation (i) = evaporation (i) + reevap * convt
            rain (i) = rain (i) + r1 * convt
            
            do k = ks, ke
                m2_rain (i, k) = m2_rain (i, k) + m1_rain (k)
                m1 (k) = m1 (k) + m1_rain (k)
            enddo
            
            ! -----------------------------------------------------------------------
            ! sedimentation of cloud ice, snow, and graupel
            ! -----------------------------------------------------------------------
            
            call fall_speed (ks, ke, den, qsz, qiz, qgz, qlz, tz, vtsz, vtiz, vtgz)
            
            call terminal_fall (dts, ks, ke, tz, qvz, qlz, qrz, qgz, qsz, qiz, &
                dz1, dp1, den, vtgz, vtsz, vtiz, r1, g1, s1, i1, m1_sol, w1, dte (i))
            
            rain (i) = rain (i) + r1 * convt
            snow (i) = snow (i) + s1 * convt
            graupel (i) = graupel (i) + g1 * convt
            ice (i) = ice (i) + i1 * convt
            
            ! -----------------------------------------------------------------------
            ! energy loss during sedimentation heating
            ! -----------------------------------------------------------------------
            
            if (consv_checker) then
                do k = ks, ke
                    te1 (k) = one_r8 + qvz (k) * c1_vap + (qlz (k) + qrz (k)) * c1_liq + &
                        (qiz (k) + qsz (k) + qgz (k)) * c1_ice
                    te1 (k) = rgrav * te1 (k) * c_air * tz (k) * dp1 (k)
                enddo
            endif
            
            ! -----------------------------------------------------------------------
            ! heat transportation during sedimentation
            ! -----------------------------------------------------------------------
            
            if (do_sedi_heat) then
                call sedi_heat (ks, ke, dp1, m1_sol, dz1, tz, qvz, qlz, qrz, qiz, &
                    qsz, qgz, c_ice)
            endif
            
            ! -----------------------------------------------------------------------
            ! energy loss during sedimentation heating
            ! -----------------------------------------------------------------------
            
            if (consv_checker) then
                do k = ks, ke
                    te2 (k) = one_r8 + qvz (k) * c1_vap + (qlz (k) + qrz (k)) * c1_liq + &
                        (qiz (k) + qsz (k) + qgz (k)) * c1_ice
                    te2 (k) = rgrav * te2 (k) * c_air * tz (k) * dp1 (k)
                enddo
                dte (i) = dte (i) + sum (te1) - sum (te2)
            endif
            
            ! -----------------------------------------------------------------------
            ! time-split warm rain processes: 2nd pass
            ! -----------------------------------------------------------------------
            
            call warm_rain (dt_rain, ks, ke, dp1, dz1, tz, qvz, qlz, qrz, qiz, qsz, &
                qgz, den, denfac, ccn, c_praut, rh_rain, vtrz, r1, m1_rain, w1, h_var, reevap, dte (i))
            
            evaporation (i) = evaporation (i) + reevap * convt
            rain (i) = rain (i) + r1 * convt
            
            do k = ks, ke
                m2_rain (i, k) = m2_rain (i, k) + m1_rain (k)
                m2_sol (i, k) = m2_sol (i, k) + m1_sol (k)
                m1 (k) = m1 (k) + m1_rain (k) + m1_sol (k)
            enddo
            
            ! -----------------------------------------------------------------------
            ! ice-phase microphysics
            ! -----------------------------------------------------------------------
            
            call icloud (ks, ke, tz, p1, qvz, qlz, qrz, qiz, qsz, qgz, dp1, den, ccn, &
                cin, denfac, vtsz, vtgz, vtrz, qaz, rh_adj, dts, h_var, gsize (i), &
                cond, dep, reevap, sub, last_step)

            condensation (i) = condensation (i) + cond * convt
            deposition (i) = deposition (i) + dep * convt
            evaporation (i) = evaporation (i) + reevap * convt
            sublimation (i) = sublimation (i) + sub * convt
            
        enddo ! n loop
        
        ! -----------------------------------------------------------------------
        ! momentum transportation during sedimentation
        ! note: dp1 is dry mass; dp0 is the old moist (total) mass
        ! -----------------------------------------------------------------------
        
        if (do_sedi_uv) then
            do k = ks + 1, ke
                u1 (k) = (dp0 (k) * u1 (k) + m1 (k - 1) * u1 (k - 1)) / (dp0 (k) + m1 (k - 1))
                v1 (k) = (dp0 (k) * v1 (k) + m1 (k - 1) * v1 (k - 1)) / (dp0 (k) + m1 (k - 1))
                ua (i, k) = u1 (k)
                va (i, k) = v1 (k)
            enddo
            if (do_disp_heat) then
                do k = ks + 1, ke
#ifdef MOIST_CAPPA
                    c8 = c_air + qvz (k) * c_vap + (qrz (k) + qlz (k)) * c_liq + &
                        (qiz (k) + qsz (k) + qgz (k)) * c_ice
                    tz (k) = tz (k) + 0.5 * (u0 (k) ** 2 + v0 (k) ** 2 - (u1 (k) ** 2 + v1 (k) ** 2)) / c8
#else
                    tz (k) = tz (k) + 0.5 * (u0 (k) ** 2 + v0 (k) ** 2 - (u1 (k) ** 2 + v1 (k) ** 2)) / c_air
#endif
                enddo
            endif
        endif
        
        if (do_sedi_w) then
            if (do_disp_heat) then
                do k = ks, ke
#ifdef MOIST_CAPPA
                    c8 = c_air + qvz (k) * c_vap + (qrz (k) + qlz (k)) * c_liq + &
                        (qiz (k) + qsz (k) + qgz (k)) * c_ice
                    tz (k) = tz (k) + 0.5 * (w (i, k) ** 2 - w1 (k) ** 2) / c8
#else
                    tz (k) = tz (k) + 0.5 * (w (i, k) ** 2 - w1 (k) ** 2) / c_air
#endif
                enddo
            endif
            do k = ks, ke
                w (i, k) = w1 (k)
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do k = ks, ke
                q_liq (k) = qlz (k) + qrz (k)
                q_sol (k) = qiz (k) + qsz (k) + qgz (k)
                cvm (k) = c_air + qvz (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                te_end (i, k) = cvm (k) * tz (k) + lv00 * c_air * qvz (k) - li00 * c_air * q_sol (k)
                if (hydrostatic) then
                    te_end (i, k) = te_end (i, k) + 0.5 * (u1 (k) ** 2 + v1 (k) ** 2)
                else
                    te_end (i, k) = te_end (i, k) + 0.5 * (u1 (k) ** 2 + v1 (k) ** 2 + w1 (k) ** 2)
                endif
                te_end (i, k) = rgrav * te_end (i, k) * dp1 (k) * gsize (i) ** 2.0
                tw_end (i, k) = rgrav * (qvz (k) + q_liq (k) + q_sol (k)) * dp1 (k) * gsize (i) ** 2.0
            enddo
            te_b_end (i) = (dte (i) - li00 * c_air * (ice (i) + snow (i) + graupel (i)) * dt_in / 86400) * &
                gsize (i) ** 2.0
            tw_b_end (i) = (rain (i) + ice (i) + snow (i) + graupel (i)) * dt_in / 86400 * gsize (i) ** 2.0
            ! total energy loss due to sedimentation and its heating
            te_loss (i) = dte (i) * gsize (i) ** 2.0
        endif
        
        ! -----------------------------------------------------------------------
        ! update moist air mass (actually hydrostatic pressure)
        ! convert to dry mixing ratios
        ! -----------------------------------------------------------------------
        
        do k = ks, ke

            ! -----------------------------------------------------------------------
            ! convert dry mixing ratios back to moist mixing ratios
            ! -----------------------------------------------------------------------

            con_r8 = one_r8 + qvz (k) + qlz (k) + qrz (k) + qiz (k) + qsz (k) + qgz (k)
            delp (i, k) = dp1 (k) * con_r8
            con_r8 = one_r8 / con_r8
            qvz (k) = qvz (k) * con_r8
            qlz (k) = qlz (k) * con_r8
            qrz (k) = qrz (k) * con_r8
            qiz (k) = qiz (k) * con_r8
            qsz (k) = qsz (k) * con_r8
            qgz (k) = qgz (k) * con_r8

            qv (i, k) = qvz (k)
            ql (i, k) = qlz (k)
            qr (i, k) = qrz (k)
            qi (i, k) = qiz (k)
            qs (i, k) = qsz (k)
            qg (i, k) = qgz (k)

            ! -----------------------------------------------------------------------
            ! conversion of temperature
            ! -----------------------------------------------------------------------
    
            q_liq (k) = qlz (k) + qrz (k)
            q_sol (k) = qiz (k) + qsz (k) + qgz (k)
            q_cond = q_liq (k) + q_sol (k)
            cvm (k) = (one_r8 - (qvz (k) + q_cond)) * c_air + qvz (k) * c_vap + q_liq (k) * c_liq + &
                q_sol (k) * c_ice

#ifdef MOIST_CAPPA
            q_con (i, k) = q_cond
            tmp = rdgas * (1. + zvir * qvz (k))
            cappa (i, k) = tmp / (tmp + cvm (k))
#endif

            if (do_inline_mp) then
#ifdef MOIST_CAPPA
                pt (i, k) = tz (k) * (1. + zvir * qvz (k)) * (1. - q_cond)
#else
                pt (i, k) = tz (k) * (1. + zvir * qvz (k))
#endif
            else
                pt (i, k) = pt (i, k) + (tz (k) - pt (i, k)) * cvm (k) / cp_air
            endif

        enddo ! k loop
        
        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do k = ks, ke
                q_liq (k) = ql (i, k) + qr (i, k)
                q_sol (k) = qi (i, k) + qs (i, k) + qg (i, k)
                cvm (k) = c_air * (1.0 - qv (i, k) - q_liq (k) - q_sol (k)) + &
                    qv (i, k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                te_end_0 (i, k) = cvm (k) * tz (k) + lv00 * c_air * qv (i, k) - li00 * c_air * q_sol (k)
                te_end_0 (i, k) = te_end_0 (i, k) + 0.5 * (ua (i, k) ** 2 + va (i, k) ** 2 + w (i, k) ** 2)
                te_end_0 (i, k) = rgrav * te_end_0 (i, k) * delp (i, k) * gsize (i) ** 2.0
                tw_end_0 (i, k) = rgrav * (qv (i, k) + q_liq (k) + q_sol (k)) * delp (i, k) * gsize (i) ** 2.0
            enddo
            te_b_end_0 (i) = (dte (i) - li00 * c_air * (ice (i) + snow (i) + graupel (i)) * dt_in / 86400) * &
                gsize (i) ** 2.0
            tw_b_end_0 (i) = (rain (i) + ice (i) + snow (i) + graupel (i)) * dt_in / 86400 * gsize (i) ** 2.0
        endif
        
        ! -----------------------------------------------------------------------
        ! calculate total energy loss or gain
        ! -----------------------------------------------------------------------
        
        if (consv_te) then
            if (hydrostatic) then
                do k = ks, ke
                    te (i, k) = te (i, k) + c_air * tz (k) * delp (i, k)
                enddo
            else
                do k = ks, ke
#ifdef MOIST_CAPPA
                    te (i, k) = te (i, k) + cvm (k) * tz (k) * delp (i, k)
#else
                    te (i, k) = te (i, k) + c_air * tz (k) * delp (i, k)
#endif
                enddo
            endif
        endif
        
        ! -----------------------------------------------------------------------
        ! update cloud fraction tendency
        ! -----------------------------------------------------------------------
        
        do k = ks, ke
            qa (i, k) = qaz (k)
        enddo
        
    enddo ! i loop
    
    ! -----------------------------------------------------------------------
    ! total energy checker
    ! -----------------------------------------------------------------------
    
    if (consv_checker) then
        if (abs (sum (te_end) + sum (te_b_end) - sum (te_beg) - sum (te_b_beg)) / &
            (sum (te_beg) + sum (te_b_beg)) .gt. te_err) then
            print *, "gfdl_mp te: ", sum (te_beg) / sum (gsize ** 2) + sum (te_b_beg) / sum (gsize ** 2), &
                sum (te_end) / sum (gsize ** 2) + sum (te_b_end) / sum (gsize ** 2), &
                 (sum (te_end) + sum (te_b_end) - sum (te_beg) - sum (te_b_beg)) / (sum (te_beg) + sum (te_b_beg))
        endif
        if (abs (sum (tw_end) + sum (tw_b_end) - sum (tw_beg) - sum (tw_b_beg)) / &
            (sum (tw_beg) + sum (tw_b_beg)) .gt. te_err) then
            print *, "gfdl_mp tw: ", sum (tw_beg) / sum (gsize ** 2) + sum (tw_b_beg) / sum (gsize ** 2), &
                sum (tw_end) / sum (gsize ** 2) + sum (tw_b_end) / sum (gsize ** 2), &
                 (sum (tw_end) + sum (tw_b_end) - sum (tw_beg) - sum (tw_b_beg)) / (sum (tw_beg) + sum (tw_b_beg))
        endif
        ! print *, "gfdl_mp te loss (%) : ", sum (te_loss) / (sum (te_beg) + sum (te_b_beg)) * 100.0
    endif
    
end subroutine mpdrv

! =======================================================================
! warm rain cloud microphysics
! includes: rain evaporation, accretion, autoconversion, and sedimentation
! =======================================================================

subroutine warm_rain (dt, ks, ke, dp, dz, tz, qv, ql, qr, qi, qs, qg, &
        den, denfac, ccn, c_praut, rh_rain, vtr, r1, m1_rain, w1, h_var, &
        reevap, dte)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke

    real, intent (in) :: dt, rh_rain, h_var

    real, intent (in), dimension (ks:ke) :: dp, dz, den, denfac, ccn, c_praut
    
    real, intent (inout), dimension (ks:ke) :: vtr, qv, ql, qr, qi, qs, qg, m1_rain, w1

    real (kind = r_grid), intent (inout) :: dte

    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz

    real, intent (out) :: r1, reevap

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    real, parameter :: so3 = 7.0 / 3.0
    
    integer :: k
    
    logical :: no_fall

    real :: sink, dq, qc, qden, zs, dt5
    
    real, dimension (ks:ke) :: dl, dm

    real, dimension (ks:ke + 1) :: ze, zt

    real (kind = r_grid), dimension (ks:ke) :: te1, te2

    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
    
    zs = 0.0
    dt5 = 0.5 * dt
    
    reevap = 0
    
    ! -----------------------------------------------------------------------
    ! rain evaporation and accretion for the first 1 / 2 time step
    ! -----------------------------------------------------------------------
    
    call revap_racc (ks, ke, dt5, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var, dp, reevap)

    ! -----------------------------------------------------------------------
    ! sedimentation of rain
    ! -----------------------------------------------------------------------
    
    m1_rain (:) = 0.0
    
    call check_column (ks, ke, qr, no_fall)

    if (no_fall) then

        vtr (:) = 0.0
        r1 = 0.0

    else
        
        ! -----------------------------------------------------------------------
        ! fall speed of rain
        ! -----------------------------------------------------------------------
        
        if (const_vr) then
            vtr (:) = vr_fac
        else
            do k = ks, ke
                qden = qr (k) * den (k)
                if (qr (k) .lt. qcmin) then
                    vtr (k) = 0.0
                else
                    vtr (k) = vr_fac * vconr * sqrt (min (10., sfcrho / den (k))) * &
                        exp (0.2 * log (qden / normr))
                    vtr (k) = min (vr_max, max (0.0, vtr (k)))
                endif
            enddo
        endif
        
        if (do_sedi_w) then
            do k = ks, ke
                dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! energy change during sedimentation
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do k = ks, ke
                te1 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te1 (k) = rgrav * te1 (k) * c_air * tz (k) * dp (k)
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! sedimentation mass flux
        ! -----------------------------------------------------------------------
        
        ze (ke + 1) = zs
        do k = ke, ks, - 1
            ze (k) = ze (k + 1) - dz (k)
        enddo
        
        if (use_ppm) then
            zt (ks) = ze (ks)
            do k = ks + 1, ke
                zt (k) = ze (k) - dt5 * (vtr (k - 1) + vtr (k))
            enddo
            zt (ke + 1) = zs - dt * vtr (ke)
            
            do k = ks, ke
                if (zt (k + 1) .ge. zt (k)) zt (k + 1) = zt (k) - dz_min
            enddo
            call lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, qr, r1, m1_rain, mono_prof)
        else
            call implicit_fall (dt, ks, ke, ze, vtr, dp, qr, r1, m1_rain)
        endif
        
        ! -----------------------------------------------------------------------
        ! energy change during sedimentation
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do k = ks, ke
                te2 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te2 (k) = rgrav * te2 (k) * c_air * tz (k) * dp (k)
            enddo
            dte = dte + sum (te1) - sum (te2)
        endif
        
        ! -----------------------------------------------------------------------
        ! vertical velocity transportation during sedimentation
        ! -----------------------------------------------------------------------
        
        if (do_sedi_w) then
            w1 (ks) = w1 (ks) + m1_rain (ks) * vtr (ks) / dm (ks)
            do k = ks + 1, ke
                w1 (k) = (dm (k) * w1 (k) + m1_rain (k - 1) * (w1 (k - 1) - vtr (k - 1)) + m1_rain (k) * vtr (k)) &
                     / (dm (k) + m1_rain (k - 1))
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! energy change during sedimentation heating
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do k = ks, ke
                te1 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te1 (k) = rgrav * te1 (k) * c_air * tz (k) * dp (k)
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! heat exchanges during sedimentation
        ! -----------------------------------------------------------------------
        
        if (do_sedi_heat) then
            call sedi_heat (ks, ke, dp, m1_rain, dz, tz, qv, ql, qr, qi, qs, qg, c_liq)
        endif
        
        ! -----------------------------------------------------------------------
        ! energy change during sedimentation heating
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do k = ks, ke
                te2 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te2 (k) = rgrav * te2 (k) * c_air * tz (k) * dp (k)
            enddo
            dte = dte + sum (te1) - sum (te2)
        endif
        
    endif ! no_fall
    
    ! -----------------------------------------------------------------------
    ! rain evaporation and accretion for the remaing 1 / 2 time step
    ! -----------------------------------------------------------------------
    
    call revap_racc (ks, ke, dt5, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var, dp, reevap)
        
    ! -----------------------------------------------------------------------
    ! autoconversion
    ! -----------------------------------------------------------------------
    
    if (irain_f .eq. 0) then
        call linear_prof (ke - ks + 1, ql (ks), dl (ks), z_slope_liq, h_var)
        do k = ks, ke
            if (tz (k) .gt. t_wfr .and. ql (k) .gt. qcmin) then
                qc = fac_rc * ccn (k)
                dl (k) = min (max (qcmin, dl (k)), 0.5 * ql (k))
                dq = 0.5 * (ql (k) + dl (k) - qc)
                if (dq .gt. 0.) then
                    sink = min (1., dq / dl (k)) * dt * c_praut (k) * den (k) * exp (so3 * log (ql (k)))
                    sink = min (ql (k), sink)
                    ql (k) = ql (k) - sink
                    qr (k) = qr (k) + sink
                endif
            endif
        enddo
    endif
    
    if (irain_f .eq. 1) then
        do k = ks, ke
            if (tz (k) .gt. t_wfr .and. ql (k) .gt. qcmin) then
                qc = fac_rc * ccn (k)
                dq = ql (k) - qc
                if (dq .gt. 0.) then
                    sink = min (dq, dt * c_praut (k) * den (k) * exp (so3 * log (ql (k))))
                    sink = min (ql (k), sink)
                    ql (k) = ql (k) - sink
                    qr (k) = qr (k) + sink
                endif
            endif
        enddo
    endif

end subroutine warm_rain

! =======================================================================
! rain evaporation and accretion
! =======================================================================

subroutine revap_racc (ks, ke, dt, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var, dp, reevap)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke

    real, intent (in) :: dt, rh_rain, h_var

    real, intent (in), dimension (ks:ke) :: den, denfac, dp

    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz

    real, intent (inout), dimension (ks:ke) :: qv, qr, ql, qi, qs, qg

    real, intent (out) :: reevap

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: k
    
    real :: dqv, qsat, dqsdt, tmp, t2, qden, q_plus, q_minus, sink
    real :: qpz, dq, dqh, tin, fac_revp, rh_tem
    
    real, dimension (ks:ke) :: q_liq, q_sol, lcpk

    real (kind = r_grid), dimension (ks:ke) :: cvm, te8

    ! -----------------------------------------------------------------------
    ! time-scale factor
    ! -----------------------------------------------------------------------
            
    if (tau_revp .gt. 1.e-6) then
        fac_revp = 1. - exp (- dt / tau_revp)
    else
        fac_revp = 1.
    endif
    
    do k = ks, ke
        
        ! -----------------------------------------------------------------------
        ! calculate heat capacities and latent heat coefficients
        ! -----------------------------------------------------------------------
        
        q_liq (k) = ql (k) + qr (k)
        q_sol (k) = qi (k) + qs (k) + qg (k)
        
        cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
        te8 (k) = cvm (k) * tz (k) + lv00 * qv (k) - li00 * q_sol (k)
        lcpk (k) = (lv00 + d1_vap * tz (k)) / cvm (k)
        tin = (tz (k) * cvm (k) - lv00 * ql (k)) / (1. + (qv (k) + ql (k)) * c1_vap + qr (k) * c1_liq + q_sol (k) * c1_ice)
        
        ! -----------------------------------------------------------------------
        ! calculate supersaturation and subgrid variability of water
        ! -----------------------------------------------------------------------
        
        qpz = qv (k) + ql (k)
        qsat = wqs (tin, den (k), dqsdt)
        dqv = qsat - qv (k)

        dqh = max (ql (k), h_var * max (qpz, qcmin))
        dqh = min (dqh, 0.2 * qpz)
        q_minus = qpz - dqh
        q_plus = qpz + dqh
        
        ! -----------------------------------------------------------------------
        ! rain evaporation
        ! -----------------------------------------------------------------------
        
        rh_tem = qpz / iqs (tin, den (k))
            
        if (tz (k) .gt. t_wfr .and. qr (k) .gt. qcmin .and. dqv .gt. 0.0 .and. qsat .gt. q_minus) then

            if (qsat .gt. q_plus) then
                dq = qsat - qpz
            else
                dq = 0.25 * (qsat - q_minus) ** 2 / dqh
            endif
            qden = qr (k) * den (k)
            t2 = tin * tin
            if (use_rhc_revap) then
                sink = 0.0
                if (rh_tem .lt. rhc_revap) then
                    sink = crevp (1) * t2 * dq * (crevp (2) * sqrt (qden) + crevp (3) * &
                        exp (0.725 * log (qden)) * sqrt (denfac (k))) / (crevp (4) * t2 + crevp (5) * qsat * den (k))
                    sink = min (qr (k), dt * fac_revp * sink, dqv / (1. + lcpk (k) * dqsdt))
                endif
            else
                sink = crevp (1) * t2 * dq * (crevp (2) * sqrt (qden) + crevp (3) * &
                    exp (0.725 * log (qden))) / (crevp (4) * t2 + crevp (5) * qsat * den (k))
                sink = min (qr (k), dt * fac_revp * sink, dqv / (1. + lcpk (k) * dqsdt))
            endif

            ! -----------------------------------------------------------------------
            ! alternative minimum evaporation in dry environmental air
            ! -----------------------------------------------------------------------
            ! tmp = min (qr (k), dim (rh_rain * qsat, qv (k)) / (1. + lcpk (k) * dqsdt))
            ! sink = max (sink, tmp)

            reevap = reevap + sink * dp (k)
            qr (k) = qr (k) - sink
            qv (k) = qv (k) + sink
            q_liq (k) = q_liq (k) - sink

            cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
            tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)

        endif
            
        ! -----------------------------------------------------------------------
        ! rain accretion
        ! -----------------------------------------------------------------------
            
        if (tz (k) .gt. t_wfr .and. qr (k) .gt. qcmin .and. ql (k) .gt. qcmin) then

            sink = dt * denfac (k) * cracw * exp (0.95 * log (qr (k) * den (k)))
            sink = sink / (1. + sink) * ql (k)
            ql (k) = ql (k) - sink
            qr (k) = qr (k) + sink

        endif

    enddo ! k loop
    
end subroutine revap_racc

! =======================================================================
! vertical subgrid variability used for cloud ice and cloud water autoconversion
! edges: qe == qbar + / - dm
! =======================================================================

subroutine linear_prof (km, q, dm, z_var, h_var)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: km

    logical, intent (in) :: z_var

    real, intent (in) :: q (km), h_var

    real, intent (out) :: dm (km)

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: k
    
    real :: dq (km)

    if (z_var) then
        do k = 2, km
            dq (k) = 0.5 * (q (k) - q (k - 1))
        enddo
        dm (1) = 0.
        ! -----------------------------------------------------------------------
        ! use twice the strength of the positive definiteness limiter (Lin et al. 1994)
        ! -----------------------------------------------------------------------
        do k = 2, km - 1
            dm (k) = 0.5 * min (abs (dq (k) + dq (k + 1)), 0.5 * q (k))
            if (dq (k) * dq (k + 1) .le. 0.) then
                if (dq (k) .gt. 0.) then
                    dm (k) = min (dm (k), dq (k), - dq (k + 1))
                else
                    dm (k) = 0.
                endif
            endif
        enddo
        dm (km) = 0.
        ! -----------------------------------------------------------------------
        ! impose a presumed background horizontal variability that is proportional to the value itself
        ! -----------------------------------------------------------------------
        do k = 1, km
            dm (k) = max (dm (k), 0.0, h_var * q (k))
        enddo
    else
        do k = 1, km
            dm (k) = max (0.0, h_var * q (k))
        enddo
    endif
    
end subroutine linear_prof

! =======================================================================
! ice cloud microphysics processes
! =======================================================================

subroutine icloud (ks, ke, tz, p1, qv, ql, qr, qi, qs, qg, dp1, den, &
        ccn, cin, denfac, vts, vtg, vtr, qa, rh_adj, dts, h_var, &
        gsize, cond, dep, reevap, sub, last_step)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke

    logical, intent (in) :: last_step

    real, intent (in) :: rh_adj, dts, h_var, gsize

    real, intent (in), dimension (ks:ke) :: p1, dp1, den, denfac, vts, vtg, vtr, ccn

    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, qa, cin

    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz

    real, intent (out) :: cond, dep, reevap, sub

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: k

    real, dimension (ks:ke) :: icpk, di, q_liq, q_sol

    real :: pracs, psacw, pgacw, psacr, pgacr, pgaci, psaci, pgacs
    real :: pgmlt, psmlt, pgfr, psaut, pgaut
    real :: tmp, dq, dtmp, q_plus, sink, factor, qim, tin
    real :: rdts, fac_i2s, fac_imlt, tc, qden, qsm
    
    real (kind = r_grid), dimension (ks:ke) :: cvm, te8

    rdts = 1. / dts
    
    ! -----------------------------------------------------------------------
    ! time-scale factor
    ! -----------------------------------------------------------------------
    
    fac_i2s = 1. - exp (- dts / tau_i2s)
    fac_imlt = 1. - exp (- dts / tau_imlt)
    
    ! -----------------------------------------------------------------------
    ! calculate heat capacities and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        q_liq (k) = ql (k) + qr (k)
        q_sol (k) = qi (k) + qs (k) + qg (k)
        cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
        te8 (k) = cvm (k) * tz (k) + lv00 * qv (k) - li00 * q_sol (k)
        icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)
    enddo
    
    if (.not. do_warm_rain_mp) then

        do k = ks, ke

            if (tz (k) .gt. tice_mlt .and. qi (k) .gt. qcmin) then
                
                ! -----------------------------------------------------------------------
                ! melting of cloud ice
                ! -----------------------------------------------------------------------
                
                sink = min (qi (k), fac_imlt * (tz (k) - tice_mlt) / icpk (k))
                tmp = min (sink, dim (ql_mlt, ql (k)))
                qi (k) = qi (k) - sink
                ql (k) = ql (k) + tmp
                qr (k) = qr (k) + sink - tmp
                q_sol (k) = q_sol (k) - sink
                q_liq (k) = q_liq (k) + sink

            elseif (tz (k) .lt. t_wfr .and. ql (k) .gt. qcmin) then
                
                ! -----------------------------------------------------------------------
                ! homogeneous freezing of cloud water
                ! -----------------------------------------------------------------------
                
                dtmp = t_wfr - tz (k)
                factor = min (1., dtmp / dt_fr)
                sink = min (ql (k) * factor, dtmp / icpk (k))
                qim = qi0_crt / den (k)
                tmp = min (sink, dim (qim, qi (k)))
                ql (k) = ql (k) - sink
                qs (k) = qs (k) + sink - tmp
                qi (k) = qi (k) + tmp
                q_liq (k) = q_liq (k) - sink
                q_sol (k) = q_sol (k) + sink

            endif
        
            cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
            tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
            icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)

        enddo
        
        ! -----------------------------------------------------------------------
        ! vertical subgrid variability
        ! -----------------------------------------------------------------------
        
        call linear_prof (ke - ks + 1, qi (ks), di (ks), z_slope_ice, h_var)
        
        do k = ks, ke
            
            ! -----------------------------------------------------------------------
            ! do nothing above p_min
            ! -----------------------------------------------------------------------
            
            if (p1 (k) .lt. p_min) cycle
            
            tc = tz (k) - tice
            
            if (tc .ge. 0.) then
                
                ! -----------------------------------------------------------------------
                ! melting of snow (includes snow accretion with cloud water and rain)
                ! -----------------------------------------------------------------------
                
                if (qs (k) .gt. qcmin) then
                    
                    psacw = 0.
                    if (ql (k) .gt. qcmin) then
                        factor = denfac (k) * csacw * exp (0.8125 * log (qs (k) * den (k)))
                        psacw = factor / (1. + dts * factor) * ql (k)
                    endif
                    
                    psacr = 0.
                    pracs = 0.
                    if (qr (k) .gt. qcmin) then
                        psacr = min (acr3d (vts (k), vtr (k), qr (k), qs (k), csacr, acco (1, 2), &
                            den (k)), qr (k) * rdts)
                        pracs = acr3d (vtr (k), vts (k), qs (k), qr (k), cracs, acco (1, 1), den (k))
                    endif
                    
                    tin = tz (k)
                    dq = wqs (tin, den (k)) - qv (k)
                    psmlt = max (0., smlt (tc, dq, qs (k) * den (k), psacw, psacr, csmlt, &
                        den (k), denfac (k)))

                    sink = min (qs (k), dts * (psmlt + pracs), tc / icpk (k))
                    tmp = min (sink, dim (qs_mlt, ql (k)))
                    qs (k) = qs (k) - sink
                    ql (k) = ql (k) + tmp
                    qr (k) = qr (k) + sink - tmp
                    q_sol (k) = q_sol (k) - sink
                    q_liq (k) = q_liq (k) + sink
                    
                    cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                    tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
                    icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)

                endif
                
                ! -----------------------------------------------------------------------
                ! melting of graupel (includes graupel accretion with cloud water and rain)
                ! -----------------------------------------------------------------------
                
                if (qg (k) .gt. qcmin) then
                    
                    pgacr = 0.
                    if (qr (k) .gt. qcmin) then
                        pgacr = min (acr3d (vtg (k), vtr (k), qr (k), qg (k), cgacr, acco (1, 3), &
                        den (k)), rdts * qr (k))
                    endif
                    
                    pgacw = 0.
                    qden = qg (k) * den (k)
                    if (ql (k) .gt. qcmin) then
                        factor = cgacw * qden / sqrt (den (k) * sqrt (sqrt (qden)))
                        pgacw = factor / (1. + dts * factor) * ql (k)
                    endif
                    
                    tin = tz (k)
                    dq = wqs (tin, den (k)) - qv (k)
                    pgmlt = dts * gmlt (tc, dq, qden, pgacw, pgacr, cgmlt, den (k))

                    sink = min (max (0., pgmlt), qg (k), tc / icpk (k))
                    qg (k) = qg (k) - sink
                    qr (k) = qr (k) + sink
                    q_sol (k) = q_sol (k) - sink
                    q_liq (k) = q_liq (k) + sink

                    cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                    tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)

                endif
                

            else
                
                ! -----------------------------------------------------------------------
                ! cloud ice sink terms (includes cloud ice accretion with snow and graupel, autoconversion)
                ! -----------------------------------------------------------------------
        
                if (qi (k) .gt. qcmin) then
                    
                    psaci = 0.
                    if (qs (k) .gt. qcmin) then
                        factor = dts * denfac (k) * csaci * exp (0.05 * tc + 0.8125 * log (qs (k) * den (k)))
                        psaci = factor / (1. + factor) * qi (k)
                    endif
                    
                    psaut = 0.
                    tmp = fac_i2s * exp (0.025 * tc)
                    di (k) = max (di (k), qcmin)
                    q_plus = qi (k) + di (k)
                    qim = qi0_crt / den (k)
                    if (q_plus .gt. (qim + qcmin)) then
                        if (qim .gt. (qi (k) - di (k))) then
                            dq = (0.25 * (q_plus - qim) ** 2) / di (k)
                        else
                            dq = qi (k) - qim
                        endif
                        psaut = tmp * dq
                    endif

                    sink = min (fi2s_fac * qi (k), psaci + psaut)
                    qi (k) = qi (k) - sink
                    qs (k) = qs (k) + sink
                    
                    pgaci = 0.
                    if (qg (k) .gt. qcmin) then
                        factor = dts * cgaci / sqrt (den (k)) * exp (0.875 * log (qg (k) * den (k)))
                        pgaci = factor / (1. + factor) * qi (k)
                    endif

                    sink = min (qi (k), pgaci)
                    qi (k) = qi (k) - sink
                    qg (k) = qg (k) + sink
                    
                endif
                
                ! -----------------------------------------------------------------------
                ! rain sink terms (includes rain accretion with snow, freezing)
                ! -----------------------------------------------------------------------
                
                if (qr (k) .gt. qcmin) then
                    
                    psacr = 0.
                    if (qs (k) .gt. qcmin) then
                        psacr = dts * acr3d (vts (k), vtr (k), qr (k), qs (k), csacr, acco (1, 2), den (k))
                    endif
                    
                    pgfr = dts * cgfr (1) / den (k) * (exp (- cgfr (2) * tc) - 1.) * &
                        exp (1.75 * log (qr (k) * den (k)))
                    
                    sink = psacr + pgfr
                    factor = min (sink, qr (k), - tc / icpk (k)) / max (sink, qcmin)
                    psacr = factor * psacr
                    pgfr = factor * pgfr
                    
                    sink = psacr + pgfr
                    qr (k) = qr (k) - sink
                    qs (k) = qs (k) + psacr
                    qg (k) = qg (k) + pgfr
                    q_liq (k) = q_liq (k) - sink
                    q_sol (k) = q_sol (k) + sink
                    
                    cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                    tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
                    icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)

                endif
                
                ! -----------------------------------------------------------------------
                ! graupel production terms (includes graupel accretion with cloud water and autoconversion)
                ! -----------------------------------------------------------------------
                
                if (qg (k) .gt. qcmin) then
                    
                    pgacs = 0
                    if (qs (k) .gt. qcmin) then
                        pgacs = dts * acr3d (vtg (k), vts (k), qs (k), qg (k), cgacs, acco (1, 4), den (k))
                    endif
                    
                    pgaut = 0
                    qsm = qs0_crt / den (k)
                    if (qs (k) .gt. qsm) then
                        factor = dts * 1.e-3 * exp (0.09 * (tz (k) - tice))
                        pgaut = factor / (1. + factor) * (qs (k) - qsm)
                    endif

                    sink = min (fs2g_fac * qs (k), pgacs + pgaut)
                    qs (k) = qs (k) - sink
                    qg (k) = qg (k) + sink
                    
                endif
                
                ! -----------------------------------------------------------------------
                ! graupel production terms (includes graupel accretion with cloud water and rain)
                ! -----------------------------------------------------------------------
                
                if (qg (k) .gt. qcmin) then
                    
                    pgacw = 0.
                    if (ql (k) .gt. qcmin) then
                        qden = qg (k) * den (k)
                        factor = dts * cgacw * qden / sqrt (den (k) * sqrt (sqrt (qden)))
                        pgacw = factor / (1. + factor) * ql (k)
                    endif
                    
                    pgacr = 0.
                    if (qr (k) .gt. qcmin) then
                        pgacr = min (dts * acr3d (vtg (k), vtr (k), qr (k), qg (k), cgacr, acco (1, 3), &
                            den (k)), qr (k))
                    endif
                    
                    sink = pgacr + pgacw
                    factor = min (sink, dim (tice, tz (k)) / icpk (k)) / max (sink, qcmin)
                    pgacr = factor * pgacr
                    pgacw = factor * pgacw
                    
                    sink = pgacr + pgacw
                    qg (k) = qg (k) + sink
                    qr (k) = qr (k) - pgacr
                    ql (k) = ql (k) - pgacw
                    q_liq (k) = q_liq (k) - sink
                    q_sol (k) = q_sol (k) + sink

                    cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                    tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)

                endif
                
            endif ! tc .ge. 0.
            
        enddo ! k loop

    endif ! do_warm_rain_mp
    
    call subgrid_z_proc (ks, ke, p1, den, denfac, dts, rh_adj, tz, qv, ql, &
        qr, qi, qs, qg, qa, dp1, h_var, te8, ccn, cin, gsize, cond, dep, &
        reevap, sub, last_step)
    
end subroutine icloud

! =======================================================================
! temperature sentive high vertical resolution processes
! =======================================================================

subroutine subgrid_z_proc (ks, ke, p1, den, denfac, dts, rh_adj, tz, qv, ql, qr, &
        qi, qs, qg, qa, dp1, h_var, te8, ccn, cin, gsize, cond, dep, reevap, sub, last_step)
    
    implicit none
    
    integer, intent (in) :: ks, ke
    real, intent (in) :: dts, rh_adj, h_var, gsize
    real, intent (in), dimension (ks:ke) :: p1, den, denfac, ccn, dp1
    real (kind = r_grid), intent (in), dimension (ks:ke) :: te8
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    real, intent (inout), dimension (ks:ke) :: cin
    logical, intent (in) :: last_step
    real, intent (out) :: cond, dep, reevap, sub
    ! local:
    real, dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    real, dimension (ks:ke) :: q_liq, q_sol, q_cond
    real (kind = r_grid), dimension (ks:ke) :: cvm
    real :: pidep, qi_gen, qi_crt
    real :: sigma, gam
    ! -----------------------------------------------------------------------
    ! qstar over water may be accurate only down to - 80 deg c with ~10% uncertainty
    ! must not be too large to allow psc
    ! -----------------------------------------------------------------------
    real :: rh, rqi, tin, qsw, qsi, qpz, qstar, rh_tem
    real :: dqsdt, dwsdt, dq, dq0, factor, tmp, liq, ice
    real :: q_plus, q_minus, dt_evap, dt_pisub
    real :: evap, sink, tc, dtmp, qa10, qa100
    real :: pssub, pgsub, tsq, qden
    real :: fac_l2v, fac_v2l
    integer :: k
    
    if (do_sat_adj) then
        dt_evap = 0.5 * dts
    else
        dt_evap = dts
    endif
    
    ! -----------------------------------------------------------------------
    ! define conversion scalar / factor
    ! -----------------------------------------------------------------------
    
    fac_l2v = 1. - exp (- dt_evap / tau_l2v)
    fac_v2l = 1. - exp (- dt_evap / tau_v2l)
    
    ! -----------------------------------------------------------------------
    ! define heat capacity and latent heat coefficient
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        q_liq (k) = ql (k) + qr (k)
        q_sol (k) = qi (k) + qs (k) + qg (k)
        cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
        lcpk (k) = (lv00 + d1_vap * tz (k)) / cvm (k)
        icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)
        tcp3 (k) = lcpk (k) + icpk (k) * min (1., dim (tice, tz (k)) / (tice - t_wfr))
    enddo

    cond = 0
    dep = 0
    reevap = 0
    sub = 0
    
    do k = ks, ke
        
        if (p1 (k) .lt. p_min) cycle
        
        if (.not. do_warm_rain_mp) then

            ! -----------------------------------------------------------------------
            ! instant deposit all water vapor to cloud ice when temperature is super low
            ! -----------------------------------------------------------------------
            
            if (tz (k) .lt. t_min) then
                sink = dim (qv (k), qcmin)
                dep = dep + sink * dp1 (k)
                qv (k) = qv (k) - sink
                qi (k) = qi (k) + sink
                q_sol (k) = q_sol (k) + sink
                tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / &
                     (one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice)
                if (do_qa) qa (k) = 1. ! air fully saturated; 100 % cloud cover
                cycle
            endif
            
            ! -----------------------------------------------------------------------
            ! instant evaporation / sublimation of all clouds if rh < rh_adj -- > cloud free
            ! -----------------------------------------------------------------------
            ! rain water is handled in warm - rain process.
            qpz = qv (k) + ql (k) + qi (k)
            tin = (te8 (k) - lv00 * qpz + li00 * (qs (k) + qg (k))) / &
                 (one_r8 + qpz * c1_vap + qr (k) * c1_liq + (qs (k) + qg (k)) * c1_ice)
            if (tin .gt. t_sub + 6.) then
                rh = qpz / iqs (tin, den (k))
                if (rh .lt. rh_adj) then ! qpz / rh_adj < qs
                    reevap = reevap + ql (k) * dp1 (k)
                    sub = sub + qi (k) * dp1 (k)
                    tz (k) = tin
                    qv (k) = qpz
                    ql (k) = 0.
                    qi (k) = 0.
                    cycle ! cloud free
                endif
            endif

        endif
        
        ! -----------------------------------------------------------------------
        ! cloud water < -- > vapor adjustment:
        ! -----------------------------------------------------------------------
        
        tin = tz (k)
        rh_tem = qpz / iqs (tin, den (k))
        qsw = wqs (tin, den (k), dwsdt)
        dq0 = qsw - qv (k)
        if (use_rhc_cevap) then
            evap = 0.
            if (rh_tem .lt. rhc_cevap) then
                if (dq0 .gt. 0.) then ! evaporation
                    factor = min (1., fac_l2v * (10. * dq0 / qsw)) ! the rh dependent factor = 1 at 90%
                    evap = min (ql (k), factor * dq0 / (1. + tcp3 (k) * dwsdt))
                    reevap = reevap + evap * dp1 (k)
                elseif (do_cond_timescale) then
                    factor = min (1., fac_v2l * (10. * (- dq0) / qsw))
                    evap = - min (qv (k), factor * (- dq0) / (1. + tcp3 (k) * dwsdt))
                    cond = cond - evap * dp1 (k)
                else ! condensate all excess vapor into cloud water
                    evap = dq0 / (1. + tcp3 (k) * dwsdt)
                    cond = cond - evap * dp1 (k)
                endif
            endif
        else
            if (dq0 .gt. 0.) then ! evaporation
                factor = min (1., fac_l2v * (10. * dq0 / qsw)) ! the rh dependent factor = 1 at 90%
                evap = min (ql (k), factor * dq0 / (1. + tcp3 (k) * dwsdt))
                reevap = reevap + evap * dp1 (k)
            elseif (do_cond_timescale) then
                factor = min (1., fac_v2l * (10. * (- dq0) / qsw))
                evap = - min (qv (k), factor * (- dq0) / (1. + tcp3 (k) * dwsdt))
                cond = cond - evap * dp1 (k)
            else ! condensate all excess vapor into cloud water
                evap = dq0 / (1. + tcp3 (k) * dwsdt)
                cond = cond - evap * dp1 (k)
            endif
        endif
        ! sjl on jan 23 2018: reversible evap / condensation:
        qv (k) = qv (k) + evap
        ql (k) = ql (k) - evap
        q_liq (k) = q_liq (k) - evap
        
        cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
        tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
        
        ! -----------------------------------------------------------------------
        ! update heat capacity and latent heat coefficient
        ! -----------------------------------------------------------------------
        
        icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)
        
        if (.not. do_warm_rain_mp) then

            ! -----------------------------------------------------------------------
            ! enforce complete freezing below - 48 c
            ! -----------------------------------------------------------------------
            
            dtmp = t_wfr - tz (k) ! [ - 40, - 48]
            if (dtmp .gt. 0. .and. ql (k) .gt. qcmin) then
                sink = min (ql (k), ql (k) * dtmp * 0.125, dtmp / icpk (k))
                ql (k) = ql (k) - sink
                qi (k) = qi (k) + sink
                q_liq (k) = q_liq (k) - sink
                q_sol (k) = q_sol (k) + sink
                cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
                icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)
            endif
            
            ! -----------------------------------------------------------------------
            ! bigg mechanism
            ! -----------------------------------------------------------------------
            
            if (do_sat_adj) then
                dt_pisub = 0.5 * dts
            else
                dt_pisub = dts
                tc = tice - tz (k)
                if (ql (k) .gt. qcmin .and. tc .gt. 0.1) then
                    sink = 100. / (rhow * ccn (k)) * dts * (exp (0.66 * tc) - 1.) * ql (k) ** 2
                    sink = min (ql (k), tc / icpk (k), sink)
                    ql (k) = ql (k) - sink
                    qi (k) = qi (k) + sink
                    q_liq (k) = q_liq (k) - sink
                    q_sol (k) = q_sol (k) + sink
                    cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                    tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
                endif ! significant ql existed
            endif
            
            ! -----------------------------------------------------------------------
            ! update capacity heat and latent heat coefficient
            ! -----------------------------------------------------------------------
            
            tcpk (k) = (li20 + (d1_vap + d1_ice) * tz (k)) / cvm (k)
            
            ! -----------------------------------------------------------------------
            ! sublimation / deposition of ice
            ! -----------------------------------------------------------------------
            
            if (tz (k) .lt. tice) then
                tin = tz (k)
                qsi = iqs (tin, den (k), dqsdt)
                dq = qv (k) - qsi
                sink = dq / (1. + tcpk (k) * dqsdt)
                if (qi (k) .gt. qcmin) then
                    if (.not. prog_ccn) then
                        if (inflag .eq. 1) &
                            ! hong et al., 2004
                            cin (k) = 5.38e7 * exp (0.75 * log (qi (k) * den (k)))
                        if (inflag .eq. 2) &
                            ! meyers et al., 1992
                            cin (k) = exp (-2.80 + 0.262 * (tice - tz (k))) * 1000.0 ! convert from L^-1 to m^-3
                        if (inflag .eq. 3) &
                            ! meyers et al., 1992
                            cin (k) = exp (-0.639 + 12.96 * (qv (k) / qsi - 1.0)) * 1000.0 ! convert from L^-1 to m^-3
                        if (inflag .eq. 4) &
                            ! cooper, 1986
                            cin (k) = 5.e-3 * exp (0.304 * (tice - tz (k))) * 1000.0 ! convert from L^-1 to m^-3
                        if (inflag .eq. 5) &
                            ! flecther, 1962
                            cin (k) = 1.e-5 * exp (0.5 * (tice - tz (k))) * 1000.0 ! convert from L^-1 to m^-3
                    endif
                    pidep = dt_pisub * dq * 4.0 * 11.9 * exp (0.5 * log (qi (k) * den (k) * cin (k))) &
                         / (qsi * den (k) * (hlv + hlf) ** 2 / (0.0243 * rvgas * tz (k) ** 2) + 4.42478e4)
                else
                    pidep = 0.
                endif
                if (dq .gt. 0.) then ! vapor - > ice
                    tmp = tice - tz (k)
                    ! -----------------------------------------------------------------------
                    ! WRF / WSM6 scheme: qi_gen = 4.92e-11 * (1.e3 * exp (0.1 * tmp)) ** 1.33
                    ! optimized: qi_gen = 4.92e-11 * exp (1.33 * log (1.e3 * exp (0.1 * tmp)))
                    ! qi_gen ~ 4.808e-7 at 0 c; 1.818e-6 at - 10 c, 9.827e-5 at - 40 c
                    ! the following value is constructed such that qi_crt = 0 at 0 c and at - 10 c matches
                    ! WRF / WSM6 ice initiation scheme; qi_crt = qi_gen * min (qi_lim, 0.1 * tmp) / den
                    ! -----------------------------------------------------------------------
                    qi_gen = 4.92e-11 * exp (1.33 * log (1.e3 * exp (0.1 * tmp)))
                    if (igflag .eq. 1) &
                        qi_crt = qi_gen / den (k)
                    if (igflag .eq. 2) &
                        qi_crt = qi_gen * min (qi_lim, 0.1 * tmp) / den (k)
                    if (igflag .eq. 3) &
                        qi_crt = 1.82e-6 * min (qi_lim, 0.1 * tmp) / den (k)
                    if (igflag .eq. 4) &
                        qi_crt = max (qi_gen, 1.82e-6) * min (qi_lim, 0.1 * tmp) / den (k)
                    sink = min (sink, max (qi_crt - qi (k), pidep), tmp / tcpk (k))
                    dep = dep + sink * dp1 (k)
                else ! ice -- > vapor
                    pidep = pidep * min (1., dim (tz (k), t_sub) * 0.2)
                    sink = max (pidep, sink, - qi (k))
                    sub = sub - sink * dp1 (k)
                endif
                qv (k) = qv (k) - sink
                qi (k) = qi (k) + sink
                q_sol (k) = q_sol (k) + sink
                cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
            endif
            
            ! -----------------------------------------------------------------------
            ! update capacity heat and latent heat coefficient
            ! -----------------------------------------------------------------------
            
            tcpk (k) = (li20 + (d1_vap + d1_ice) * tz (k)) / cvm (k)
            
            ! -----------------------------------------------------------------------
            ! sublimation / deposition of snow
            ! this process happens for all temp rage
            ! -----------------------------------------------------------------------
            
            if (qs (k) .gt. qcmin) then
                tin = tz (k)
                qsi = iqs (tin, den (k), dqsdt)
                qden = qs (k) * den (k)
                tmp = exp (0.65625 * log (qden))
                tsq = tz (k) * tz (k)
                dq = (qsi - qv (k)) / (1. + tcpk (k) * dqsdt)
                pssub = cssub (1) * tsq * (cssub (2) * sqrt (qden) + cssub (3) * tmp * &
                    sqrt (denfac (k))) / (cssub (4) * tsq + cssub (5) * qsi * den (k))
                pssub = (qsi - qv (k)) * dts * pssub
                if (pssub .gt. 0.) then ! qs -- > qv, sublimation
                    pssub = min (pssub * min (1., dim (tz (k), t_sub) * 0.2), qs (k))
                    sub = sub + pssub * dp1 (k)
                else
                    if (tz (k) .gt. tice) then
                        pssub = 0. ! no deposition
                    else
                        pssub = max (pssub, dq, (tz (k) - tice) / tcpk (k))
                    endif
                    dep = dep - pssub * dp1 (k)
                endif
                qs (k) = qs (k) - pssub
                qv (k) = qv (k) + pssub
                q_sol (k) = q_sol (k) - pssub
                cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
                tcpk (k) = (li20 + (d1_vap + d1_ice) * tz (k)) / cvm (k)
            endif
            
            ! -----------------------------------------------------------------------
            ! sublimation / deposition of graupel
            ! this process happens for all temp rage
            ! -----------------------------------------------------------------------
            
            if (qg (k) .gt. qcmin) then
                tin = tz (k)
                qsi = iqs (tin, den (k), dqsdt)
                qden = qg (k) * den (k)
                tmp = exp (0.6875 * log (qden))
                tsq = tz (k) * tz (k)
                dq = (qsi - qv (k)) / (1. + tcpk (k) * dqsdt)
                pgsub = cgsub (1) * tsq * (cgsub (2) * sqrt (qden) + cgsub (3) * tmp / &
                    sqrt (sqrt (den (k)))) / (cgsub (4) * tsq + cgsub (5) * qsi * den (k))
                pgsub = (qsi - qv (k)) * dts * pgsub
                if (pgsub .gt. 0.) then ! qs -- > qv, sublimation
                    pgsub = min (pgsub * min (1., dim (tz (k), t_sub) * 0.2), qg (k))
                    sub = sub + pgsub * dp1 (k)
                else
                    if (tz (k) .gt. tice) then
                        pgsub = 0. ! no deposition
                    else
                        pgsub = max (pgsub, dq, (tz (k) - tice) / tcpk (k))
                    endif
                    dep = dep - pgsub * dp1 (k)
                endif
                qg (k) = qg (k) - pgsub
                qv (k) = qv (k) + pgsub
                q_sol (k) = q_sol (k) - pgsub
                cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)
                tcpk (k) = (li20 + (d1_vap + d1_ice) * tz (k)) / cvm (k)
            endif
            
        endif
        
        ! -----------------------------------------------------------------------
        ! compute cloud fraction
        ! -----------------------------------------------------------------------
        
        ! -----------------------------------------------------------------------
        ! combine water species
        ! -----------------------------------------------------------------------
        
        if (.not. (do_qa .and. last_step)) cycle
        
        ice = q_sol (k)
        if (rad_snow) then
            if (rad_graupel) then
                q_sol (k) = qi (k) + qs (k) + qg (k)
            else
                q_sol (k) = qi (k) + qs (k)
            endif
        else
            q_sol (k) = qi (k)
        endif
        liq = q_liq (k)
        if (rad_rain) then
            q_liq (k) = ql (k) + qr (k)
        else
            q_liq (k) = ql (k)
        endif
        
        q_cond (k) = q_liq (k) + q_sol (k)
        qpz = qv (k) + q_cond (k)
        
        ! -----------------------------------------------------------------------
        ! use the "liquid - frozen water temperature" (tin) to compute saturated specific humidity
        ! -----------------------------------------------------------------------
        ! tin = tz (k) - (lcpk (k) * q_cond (k) + icpk (k) * q_sol (k)) ! minimum temperature
        !! tin = (tz (k) * cvm (i) + li00 * q_sol (k) - lv00 * q_cond (k)) / &
        !! (one_r8 + (qv (k) + q_cond (k)) * c1_vap)
        ice = ice - q_sol (k)
        liq = liq - q_liq (k)
        tin = (te8 (k) - lv00 * qpz + li00 * ice) / (one_r8 + qpz * c1_vap + liq * c1_liq + ice * c1_ice)
        ! -----------------------------------------------------------------------
        ! determine saturated specific humidity
        ! -----------------------------------------------------------------------
        
        if (tin .le. t_wfr) then
            ! ice phase:
            qstar = iqs (tin, den (k))
        elseif (tin .ge. tice) then
            ! liquid phase:
            qstar = wqs (tin, den (k))
        else
            ! mixed phase:
            qsi = iqs (tin, den (k))
            qsw = wqs (tin, den (k))
            if (q_cond (k) .gt. qcmin) then
                rqi = q_sol (k) / q_cond (k)
            else
                ! -----------------------------------------------------------------------
                ! mostly liquid water q_cond (k) at initial cloud development stage
                ! -----------------------------------------------------------------------
                rqi = (tice - tin) / (tice - t_wfr)
            endif
            qstar = rqi * qsi + (1. - rqi) * qsw
        endif
        
        ! -----------------------------------------------------------------------
        ! assuming subgrid linear distribution in horizontal; this is effectively a smoother for the
        ! binary cloud scheme
        ! -----------------------------------------------------------------------
        
        ! -----------------------------------------------------------------------
        ! partial cloudiness by pdf:
        ! assuming subgrid linear distribution in horizontal; this is effectively a smoother for the
        ! binary cloud scheme; qa = 0.5 if qstar == qpz
        ! -----------------------------------------------------------------------
        
        rh = qpz / qstar
        
        if (cfflag .eq. 1) then
            if (rh .gt. rh_thres .and. qpz .gt. qcmin) then
                
                dq = h_var * qpz
                if (do_cld_adj) then
                    q_plus = qpz + dq * f_dq_p * min(1.0, max(0.0, (p1 (k) - 200.e2) / (1000.e2 - 200.e2)))
                else
                    q_plus = qpz + dq * f_dq_p
                endif
                q_minus = qpz - dq * f_dq_m
                
                if (icloud_f .eq. 2) then
                    if (qstar .lt. qpz) then
                        qa (k) = 1.
                    else
                        qa (k) = 0.
                    endif
                elseif (icloud_f .eq. 3) then
                    if (qstar .lt. qpz) then
                        qa (k) = 1.
                    else
                        if (qstar .lt. q_plus) then
                            qa (k) = (q_plus - qstar) / (dq * f_dq_p)
                        else
                            qa (k) = 0.
                        endif
                        ! impose minimum cloudiness if substantial q_cond (k) exist
                        if (q_cond (k) .gt. qcmin) then
                            qa (k) = max (cld_min, qa (k))
                        endif
                        qa (k) = min (1., qa (k))
                    endif
                else
                    if (qstar .lt. q_minus) then
                        qa (k) = 1.
                    else
                        if (qstar .lt. q_plus) then
                            if (icloud_f .eq. 0) then
                                qa (k) = (q_plus - qstar) / (dq * f_dq_p + dq * f_dq_m)
                            else
                                qa (k) = (q_plus - qstar) / ((dq * f_dq_p + dq * f_dq_m) * (1. - q_cond (k)))
                            endif
                        else
                            qa (k) = 0.
                        endif
                        ! impose minimum cloudiness if substantial q_cond (k) exist
                        if (q_cond (k) .gt. qcmin) then
                            qa (k) = max (cld_min, qa (k))
                        endif
                        qa (k) = min (1., qa (k))
                    endif
                endif
            else
                qa (k) = 0.
            endif
        endif
        
        if (cfflag .eq. 2) then
            if (rh .ge. 1.0) then
                qa (k) = 1.0
            elseif (rh .gt. rh_thres .and. q_cond (k) .gt. qcmin) then
                qa (k) = rh ** xr_a * (1.0 - exp (- xr_b * max (0.0, q_cond (k)) / &
                    max (1.e-5, (max (1.e-10, 1.0 - rh) * qstar) ** xr_c)))
                qa (k) = max (0.0, min (1., qa (k)))
            else
                qa (k) = 0.0
            endif
        endif

        if (cfflag .eq. 3) then
            if (q_cond (k) .gt. qcmin) then
                qa (k) = 1. / 50. * (5.77 * (100. - gsize / 1000.) * max (0.0, q_cond (k) * 1000.) ** 1.07 + &
                    4.82 * (gsize / 1000. - 50.) * max (0.0, q_cond (k) * 1000.) ** 0.94)
                qa (k) = qa (k) * (0.92 / 0.96 * q_liq (k) / q_cond (k) + 1.0 / 0.96 * q_sol (k) / q_cond (k))
                qa (k) = max (0.0, min (1., qa (k)))
            else
                qa (k) = 0.0
            endif
        endif

        if (cfflag .eq. 4) then
            sigma = 0.28 + max (0.0, q_cond (k) * 1000.) ** 0.49
            gam = max (0.0, q_cond (k) * 1000.) / sigma
            if (gam .lt. 0.18) then
                qa10 = 0.
            elseif (gam .gt. 2.0) then
                qa10 = 1.0
            else
                qa10 = - 0.1754 + 0.9811 * gam - 0.2223 * gam ** 2 + 0.0104 * gam ** 3
                qa10 = max (0.0, min (1., qa10))
            endif
            if (gam .lt. 0.12) then
                qa100 = 0.
            elseif (gam .gt. 1.85) then
                qa100 = 1.0
            else
                qa100 = - 0.0913 + 0.7213 * gam + 0.1060 * gam ** 2 - 0.0946 * gam ** 3
                qa100 = max (0.0, min (1., qa100))
            endif
            qa (k) = qa10 + (log10 (gsize / 1000.) - 1) * (qa100 - qa10)
            qa (k) = max (0.0, min (1., qa (k)))
        endif

    enddo
    
end subroutine subgrid_z_proc

! =======================================================================
! compute terminal fall speed
! consider cloud ice, snow, and graupel's melting during fall
! =======================================================================

subroutine terminal_fall (dtm, ks, ke, tz, qv, ql, qr, qg, qs, qi, dz, dp, &
        den, vtg, vts, vti, r1, g1, s1, i1, m1_sol, w1, dte)
    
    implicit none
    
    integer, intent (in) :: ks, ke
    real, intent (in) :: dtm ! time step (s)
    real, intent (in), dimension (ks:ke) :: vtg, vts, vti, den, dp, dz
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qg, qs, qi, m1_sol, w1
    real (kind = r_grid), intent (inout) :: dte
    real, intent (out) :: r1, g1, s1, i1
    ! local:
    real, dimension (ks:ke + 1) :: ze, zt
    real :: qsat, dqsdt, dt5, evap, dtime
    real :: factor, frac
    real :: tmp, precip, tc, sink
    real, dimension (ks:ke) :: lcpk, icpk, cvm, q_liq, q_sol
    real, dimension (ks:ke) :: m1, dm
    real (kind = r_grid), dimension (ks:ke) :: te1, te2
    real :: zs = 0.
    real :: fac_imlt
    
    integer :: k, k0, m
    logical :: no_fall
    
    dt5 = 0.5 * dtm
    fac_imlt = 1. - exp (- dtm / tau_imlt)
    
    ! -----------------------------------------------------------------------
    ! define heat capacity and latent heat coefficient
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        m1_sol (k) = 0.
        q_liq (k) = ql (k) + qr (k)
        q_sol (k) = qi (k) + qs (k) + qg (k)
        cvm (k) = 1. + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
        lcpk (k) = (lv00 + d1_vap * tz (k)) / cvm (k)
        icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)
    enddo
    
    ! -----------------------------------------------------------------------
    ! find significant melting level
    ! -----------------------------------------------------------------------
    
    k0 = ke
    do k = ks, ke - 1
        if (tz (k) .gt. tice) then
            k0 = k
            exit
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! melting of cloud_ice (before fall) :
    ! -----------------------------------------------------------------------
    
    do k = k0, ke
        tc = tz (k) - tice
        if (qi (k) .gt. qcmin .and. tc .gt. 0.) then
            sink = min (qi (k), fac_imlt * tc / icpk (k))
            tmp = min (sink, dim (ql_mlt, ql (k)))
            ql (k) = ql (k) + tmp
            qr (k) = qr (k) + sink - tmp
            qi (k) = qi (k) - sink
            q_liq (k) = q_liq (k) + sink
            q_sol (k) = q_sol (k) - sink
            tz (k) = tz (k) * cvm (k) - li00 * sink
            cvm (k) = 1. + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
            tz (k) = tz (k) / cvm (k)
            tc = tz (k) - tice
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! turn off melting when cloud microphysics time step is small
    ! -----------------------------------------------------------------------
    
    ! sjl, turn off melting of falling cloud ice, snow and graupel
    ! if (dtm < 60.) k0 = ke
    k0 = ke
    ! sjl, turn off melting of falling cloud ice, snow and graupel
    
    ze (ke + 1) = zs
    do k = ke, ks, - 1
        ze (k) = ze (k + 1) - dz (k) ! dz < 0
    enddo
    
    zt (ks) = ze (ks)
    
    ! -----------------------------------------------------------------------
    ! update capacity heat and latent heat coefficient
    ! -----------------------------------------------------------------------
    
    do k = k0, ke
        icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)
    enddo
    
    ! -----------------------------------------------------------------------
    ! melting of falling cloud ice into cloud water and rain
    ! -----------------------------------------------------------------------
    
    call check_column (ks, ke, qi, no_fall)
    
    if (no_fall) then
        i1 = 0.
    else
        
        do k = ks + 1, ke
            zt (k) = ze (k) - dt5 * (vti (k - 1) + vti (k))
        enddo
        zt (ke + 1) = zs - dtm * vti (ke)
        
        do k = ks, ke
            if (zt (k + 1) .ge. zt (k)) zt (k + 1) = zt (k) - dz_min
        enddo
        
        if (k0 .lt. ke) then
            do k = ke - 1, k0, - 1
                if (qi (k) .gt. qcmin) then
                    do m = k + 1, ke
                        if (zt (k + 1) .ge. ze (m)) exit
                        if (zt (k) .lt. ze (m + 1) .and. tz (m) .gt. tice) then
                            dtime = min (1.0, (ze (m) - ze (m + 1)) / (max (0.0, vti (k)) * tau_imlt))
                            sink = min (qi (k) * dp (k) / dp (m), dtime * (tz (m) - tice) / icpk (m))
                            tmp = min (sink, dim (ql_mlt, ql (m)))
                            ql (m) = ql (m) + tmp
                            qr (m) = qr (m) - tmp + sink
                            qi (k) = qi (k) - sink * dp (m) / dp (k)
                            tz (m) = (tz (m) * cvm (m) - li00 * sink) / &
                                 (1. + qv (m) * c1_vap + (ql (m) + qr (m)) * c1_liq + (qi (m) + qs (m) + qg (m)) * c1_ice)
                        endif
                    enddo
                endif
            enddo
        endif
        
        if (do_sedi_w) then
            do k = ks, ke
                dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do k = ks, ke
                te1 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te1 (k) = rgrav * te1 (k) * c_air * tz (k) * dp (k)
            enddo
        endif
        
        if (use_ppm) then
            call lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, qi, i1, m1_sol, mono_prof)
        else
            call implicit_fall (dtm, ks, ke, ze, vti, dp, qi, i1, m1_sol)
        endif
        
        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do k = ks, ke
                te2 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te2 (k) = rgrav * te2 (k) * c_air * tz (k) * dp (k)
            enddo
            dte = dte + sum (te1) - sum (te2)
        endif
        
        if (do_sedi_w) then
            w1 (ks) = w1 (ks) + m1_sol (ks) * vti (ks) / dm (ks)
            do k = ks + 1, ke
                w1 (k) = (dm (k) * w1 (k) + m1_sol (k - 1) * (w1 (k - 1) - vti (k - 1)) + m1_sol (k) * vti (k)) &
                     / (dm (k) + m1_sol (k - 1))
            enddo
        endif
        
    endif
    
    ! -----------------------------------------------------------------------
    ! melting of falling snow into rain
    ! -----------------------------------------------------------------------
    
    r1 = 0.
    
    call check_column (ks, ke, qs, no_fall)
    
    if (no_fall) then
        s1 = 0.
    else
        
        do k = ks + 1, ke
            zt (k) = ze (k) - dt5 * (vts (k - 1) + vts (k))
        enddo
        zt (ke + 1) = zs - dtm * vts (ke)
        
        do k = ks, ke
            if (zt (k + 1) .ge. zt (k)) zt (k + 1) = zt (k) - dz_min
        enddo
        
        if (k0 .lt. ke) then
            do k = ke - 1, k0, - 1
                if (qs (k) .gt. qcmin) then
                    do m = k + 1, ke
                        if (zt (k + 1) .ge. ze (m)) exit
                        dtime = min (dtm, (ze (m) - ze (m + 1)) / (0.0 + vts (k)))
                        if (zt (k) .lt. ze (m + 1) .and. tz (m) .gt. tice) then
                            dtime = min (1.0, dtime / tau_smlt)
                            sink = min (qs (k) * dp (k) / dp (m), dtime * (tz (m) - tice) / icpk (m))
                            tz (m) = tz (m) - sink * icpk (m)
                            qs (k) = qs (k) - sink * dp (m) / dp (k)
                            if (zt (k) .lt. zs) then
                                r1 = r1 + sink * dp (m) ! precip as rain
                            else
                                ! qr source here will fall next time step (therefore, can evap)
                                qr (m) = qr (m) + sink
                            endif
                        endif
                        if (qs (k) .lt. qcmin) exit
                    enddo
                endif
            enddo
        endif
        
        if (do_sedi_w) then
            do k = ks, ke
                dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do k = ks, ke
                te1 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te1 (k) = rgrav * te1 (k) * c_air * tz (k) * dp (k)
            enddo
        endif
        
        if (use_ppm) then
            call lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, qs, s1, m1, mono_prof)
        else
            call implicit_fall (dtm, ks, ke, ze, vts, dp, qs, s1, m1)
        endif
        
        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do k = ks, ke
                te2 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te2 (k) = rgrav * te2 (k) * c_air * tz (k) * dp (k)
            enddo
            dte = dte + sum (te1) - sum (te2)
        endif
        
        do k = ks, ke
            m1_sol (k) = m1_sol (k) + m1 (k)
        enddo
        
        if (do_sedi_w) then
            w1 (ks) = w1 (ks) + m1 (ks) * vts (ks) / dm (ks)
            do k = ks + 1, ke
                w1 (k) = (dm (k) * w1 (k) + m1 (k - 1) * (w1 (k - 1) - vts (k - 1)) + m1 (k) * vts (k)) &
                     / (dm (k) + m1 (k - 1))
            enddo
        endif
        
    endif
    
    ! ----------------------------------------------
    ! melting of falling graupel into rain
    ! ----------------------------------------------
    
    call check_column (ks, ke, qg, no_fall)
    
    if (no_fall) then
        g1 = 0.
    else
        
        do k = ks + 1, ke
            zt (k) = ze (k) - dt5 * (vtg (k - 1) + vtg (k))
        enddo
        zt (ke + 1) = zs - dtm * vtg (ke)
        
        do k = ks, ke
            if (zt (k + 1) .ge. zt (k)) zt (k + 1) = zt (k) - dz_min
        enddo
        
        if (k0 .lt. ke) then
            do k = ke - 1, k0, - 1
                if (qg (k) .gt. qcmin) then
                    do m = k + 1, ke
                        if (zt (k + 1) .ge. ze (m)) exit
                        dtime = min (dtm, (ze (m) - ze (m + 1)) / vtg (k))
                        if (zt (k) .lt. ze (m + 1) .and. tz (m) .gt. tice) then
                            dtime = min (1., dtime / tau_g2r)
                            sink = min (qg (k) * dp (k) / dp (m), dtime * (tz (m) - tice) / icpk (m))
                            tz (m) = tz (m) - sink * icpk (m)
                            qg (k) = qg (k) - sink * dp (m) / dp (k)
                            if (zt (k) .lt. zs) then
                                r1 = r1 + sink * dp (m)
                            else
                                qr (m) = qr (m) + sink
                            endif
                        endif
                        if (qg (k) .lt. qcmin) exit
                    enddo
                endif
            enddo
        endif
        
        if (do_sedi_w) then
            do k = ks, ke
                dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do k = ks, ke
                te1 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te1 (k) = rgrav * te1 (k) * c_air * tz (k) * dp (k)
            enddo
        endif
        
        if (use_ppm) then
            call lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, qg, g1, m1, mono_prof)
        else
            call implicit_fall (dtm, ks, ke, ze, vtg, dp, qg, g1, m1)
        endif
        
        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do k = ks, ke
                te2 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te2 (k) = rgrav * te2 (k) * c_air * tz (k) * dp (k)
            enddo
            dte = dte + sum (te1) - sum (te2)
        endif
        
        do k = ks, ke
            m1_sol (k) = m1_sol (k) + m1 (k)
        enddo
        
        if (do_sedi_w) then
            w1 (ks) = w1 (ks) + m1 (ks) * vtg (ks) / dm (ks)
            do k = ks + 1, ke
                w1 (k) = (dm (k) * w1 (k) + m1 (k - 1) * (w1 (k - 1) - vtg (k - 1)) + m1 (k) * vtg (k)) &
                     / (dm (k) + m1 (k - 1))
            enddo
        endif
        
    endif
    
end subroutine terminal_fall

! =======================================================================
! time - implicit monotonic scheme
! developed by sj lin, 2016
! =======================================================================

subroutine implicit_fall (dt, ks, ke, ze, vt, dp, q, precip, m1)
    
    implicit none
    
    integer, intent (in) :: ks, ke
    real, intent (in) :: dt
    real, intent (in), dimension (ks:ke + 1) :: ze
    real, intent (in), dimension (ks:ke) :: vt, dp
    real, intent (inout), dimension (ks:ke) :: q
    real, intent (out), dimension (ks:ke) :: m1
    real, intent (out) :: precip
    real, dimension (ks:ke) :: dz, qm, dd
    integer :: k
    
    do k = ks, ke
        dz (k) = ze (k) - ze (k + 1)
        dd (k) = dt * vt (k)
        q (k) = q (k) * dp (k)
    enddo
    
    ! -----------------------------------------------------------------------
    ! sedimentation: non - vectorizable loop
    ! -----------------------------------------------------------------------
    
    qm (ks) = q (ks) / (dz (ks) + dd (ks))
    do k = ks + 1, ke
        qm (k) = (q (k) + dd (k - 1) * qm (k - 1)) / (dz (k) + dd (k))
    enddo
    
    ! -----------------------------------------------------------------------
    ! qm is density at this stage
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        qm (k) = qm (k) * dz (k)
    enddo
    
    ! -----------------------------------------------------------------------
    ! output mass fluxes: non - vectorizable loop
    ! -----------------------------------------------------------------------
    
    m1 (ks) = q (ks) - qm (ks)
    do k = ks + 1, ke
        m1 (k) = m1 (k - 1) + q (k) - qm (k)
    enddo
    precip = m1 (ke)
    
    ! -----------------------------------------------------------------------
    ! update:
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        q (k) = qm (k) / dp (k) !dry dp used inside MP
    enddo
    
end subroutine implicit_fall

! =======================================================================
! lagrangian scheme
! developed by sj lin, around 2006
! =======================================================================

subroutine lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, q, precip, m1, mono)
    
    implicit none
    
    integer, intent (in) :: ks, ke
    real, intent (in) :: zs
    logical, intent (in) :: mono
    real, intent (in), dimension (ks:ke + 1) :: ze, zt
    real, intent (in), dimension (ks:ke) :: dp
    
    ! m1: flux
    real, intent (inout), dimension (ks:ke) :: q, m1
    real, intent (out) :: precip
    real, dimension (ks:ke) :: qm, dz
    
    real :: a4 (4, ks:ke)
    real :: pl, pr, delz, esl
    integer :: k, k0, n, m
    real, parameter :: r3 = 1. / 3., r23 = 2. / 3.
    
    ! -----------------------------------------------------------------------
    ! density:
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        dz (k) = zt (k) - zt (k + 1) ! note: dz is positive
        q (k) = q (k) * dp (k)
        a4 (1, k) = q (k) / dz (k)
        qm (k) = 0.
    enddo
    
    ! -----------------------------------------------------------------------
    ! construct vertical profile with zt as coordinate
    ! -----------------------------------------------------------------------
    
    call cs_profile (a4 (1, ks), dz (ks), ke - ks + 1, mono)
    
    k0 = ks
    do k = ks, ke
        do n = k0, ke
            if (ze (k) .le. zt (n) .and. ze (k) .ge. zt (n + 1)) then
                pl = (zt (n) - ze (k)) / dz (n)
                if (zt (n + 1) .le. ze (k + 1)) then
                    ! entire new grid is within the original grid
                    pr = (zt (n) - ze (k + 1)) / dz (n)
                    qm (k) = a4 (2, n) + 0.5 * (a4 (4, n) + a4 (3, n) - a4 (2, n)) * (pr + pl) - &
                        a4 (4, n) * r3 * (pr * (pr + pl) + pl ** 2)
                    qm (k) = qm (k) * (ze (k) - ze (k + 1))
                    k0 = n
                    goto 555
                else
                    qm (k) = (ze (k) - zt (n + 1)) * (a4 (2, n) + 0.5 * (a4 (4, n) + &
                        a4 (3, n) - a4 (2, n)) * (1. + pl) - a4 (4, n) * (r3 * (1. + pl * (1. + pl))))
                    if (n .lt. ke) then
                        do m = n + 1, ke
                            ! locate the bottom edge: ze (k + 1)
                            if (ze (k + 1) .lt. zt (m + 1)) then
                                qm (k) = qm (k) + q (m)
                            else
                                delz = zt (m) - ze (k + 1)
                                esl = delz / dz (m)
                                qm (k) = qm (k) + delz * (a4 (2, m) + 0.5 * esl * &
                                     (a4 (3, m) - a4 (2, m) + a4 (4, m) * (1. - r23 * esl)))
                                k0 = m
                                goto 555
                            endif
                        enddo
                    endif
                    goto 555
                endif
            endif
        enddo
        555 continue
    enddo
    
    m1 (ks) = q (ks) - qm (ks)
    do k = ks + 1, ke
        m1 (k) = m1 (k - 1) + q (k) - qm (k)
    enddo
    precip = m1 (ke)
    
    ! convert back to * dry * mixing ratio:
    ! dp must be dry air_mass (because moist air mass will be changed due to terminal fall) .
    
    do k = ks, ke
        q (k) = qm (k) / dp (k)
    enddo
    
end subroutine lagrangian_fall_ppm

! =======================================================================
! =======================================================================

subroutine cs_profile (a4, del, km, do_mono)
    
    implicit none
    
    integer, intent (in) :: km ! vertical dimension
    real, intent (in) :: del (km)
    logical, intent (in) :: do_mono
    real, intent (inout) :: a4 (4, km)
    real, parameter :: qp_min = 1.e-6
    real :: gam (km)
    real :: q (km + 1)
    real :: d4, bet, a_bot, grat, pmp, lac
    real :: pmp_1, lac_1, pmp_2, lac_2
    real :: da1, da2, a6da
    
    integer :: k
    
    logical extm (km)
    
    grat = del (2) / del (1) ! grid ratio
    bet = grat * (grat + 0.5)
    q (1) = (2. * grat * (grat + 1.) * a4 (1, 1) + a4 (1, 2)) / bet
    gam (1) = (1. + grat * (grat + 1.5)) / bet
    
    do k = 2, km
        d4 = del (k - 1) / del (k)
        bet = 2. + 2. * d4 - gam (k - 1)
        q (k) = (3. * (a4 (1, k - 1) + d4 * a4 (1, k)) - q (k - 1)) / bet
        gam (k) = d4 / bet
    enddo
    
    a_bot = 1. + d4 * (d4 + 1.5)
    q (km + 1) = (2. * d4 * (d4 + 1.) * a4 (1, km) + a4 (1, km - 1) - a_bot * q (km)) &
         / (d4 * (d4 + 0.5) - a_bot * gam (km))
    
    do k = km, 1, - 1
        q (k) = q (k) - gam (k) * q (k + 1)
    enddo
    
    ! -----------------------------------------------------------------------
    ! apply constraints
    ! -----------------------------------------------------------------------
    
    do k = 2, km
        gam (k) = a4 (1, k) - a4 (1, k - 1)
    enddo
    
    ! -----------------------------------------------------------------------
    ! apply large - scale constraints to all fields if not local max / min
    ! -----------------------------------------------------------------------
    
    ! -----------------------------------------------------------------------
    ! top:
    ! -----------------------------------------------------------------------
    
    q (1) = max (q (1), 0.)
    q (2) = min (q (2), max (a4 (1, 1), a4 (1, 2)))
    q (2) = max (q (2), min (a4 (1, 1), a4 (1, 2)), 0.)
    
    ! -----------------------------------------------------------------------
    ! interior:
    ! -----------------------------------------------------------------------
    
    do k = 3, km - 1
        if (gam (k - 1) * gam (k + 1) .gt. 0.) then
            q (k) = min (q (k), max (a4 (1, k - 1), a4 (1, k)))
            q (k) = max (q (k), min (a4 (1, k - 1), a4 (1, k)))
        else
            if (gam (k - 1) .gt. 0.) then
                ! there exists a local max
                q (k) = max (q (k), min (a4 (1, k - 1), a4 (1, k)))
            else
                ! there exists a local min
                q (k) = min (q (k), max (a4 (1, k - 1), a4 (1, k)))
                q (k) = max (q (k), 0.0)
            endif
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! bottom :
    ! -----------------------------------------------------------------------
    
    q (km) = min (q (km), max (a4 (1, km - 1), a4 (1, km)))
    q (km) = max (q (km), min (a4 (1, km - 1), a4 (1, km)), 0.)
    ! q (km + 1) = max (q (km + 1), 0.)
    
    ! -----------------------------------------------------------------------
    ! f (s) = al + s * [ (ar - al) + a6 * (1 - s) ] (0 <= s <= 1)
    ! -----------------------------------------------------------------------
    
    do k = 1, km - 1
        a4 (2, k) = q (k)
        a4 (3, k) = q (k + 1)
    enddo
    
    do k = 2, km - 1
        if (gam (k) * gam (k + 1) .gt. 0.0) then
            extm (k) = .false.
        else
            extm (k) = .true.
        endif
    enddo
    
    if (do_mono) then
        do k = 3, km - 2
            if (extm (k)) then
                ! positive definite constraint only if true local extrema
                if (a4 (1, k) .lt. qp_min .or. extm (k - 1) .or. extm (k + 1)) then
                    a4 (2, k) = a4 (1, k)
                    a4 (3, k) = a4 (1, k)
                endif
            else
                a4 (4, k) = 6. * a4 (1, k) - 3. * (a4 (2, k) + a4 (3, k))
                if (abs (a4 (4, k)) .gt. abs (a4 (2, k) - a4 (3, k))) then
                    ! check within the smooth region if subgrid profile is non - monotonic
                    pmp_1 = a4 (1, k) - 2.0 * gam (k + 1)
                    lac_1 = pmp_1 + 1.5 * gam (k + 2)
                    a4 (2, k) = min (max (a4 (2, k), min (a4 (1, k), pmp_1, lac_1)), &
                        max (a4 (1, k), pmp_1, lac_1))
                    pmp_2 = a4 (1, k) + 2.0 * gam (k)
                    lac_2 = pmp_2 - 1.5 * gam (k - 1)
                    a4 (3, k) = min (max (a4 (3, k), min (a4 (1, k), pmp_2, lac_2)), &
                        max (a4 (1, k), pmp_2, lac_2))
                endif
            endif
        enddo
    else
        do k = 3, km - 2
            if (extm (k)) then
                if (a4 (1, k) .lt. qp_min .or. extm (k - 1) .or. extm (k + 1)) then
                    a4 (2, k) = a4 (1, k)
                    a4 (3, k) = a4 (1, k)
                endif
            endif
        enddo
    endif
    
    do k = 1, km - 1
        a4 (4, k) = 6. * a4 (1, k) - 3. * (a4 (2, k) + a4 (3, k))
    enddo
    
    k = km - 1
    if (extm (k)) then
        a4 (2, k) = a4 (1, k)
        a4 (3, k) = a4 (1, k)
        a4 (4, k) = 0.
    else
        da1 = a4 (3, k) - a4 (2, k)
        da2 = da1 ** 2
        a6da = a4 (4, k) * da1
        if (a6da .lt. - da2) then
            a4 (4, k) = 3. * (a4 (2, k) - a4 (1, k))
            a4 (3, k) = a4 (2, k) - a4 (4, k)
        elseif (a6da .gt. da2) then
            a4 (4, k) = 3. * (a4 (3, k) - a4 (1, k))
            a4 (2, k) = a4 (3, k) - a4 (4, k)
        endif
    endif
    
    call cs_limiters (km - 1, a4)
    
    ! -----------------------------------------------------------------------
    ! bottom layer:
    ! -----------------------------------------------------------------------
    
    a4 (2, km) = a4 (1, km)
    a4 (3, km) = a4 (1, km)
    a4 (4, km) = 0.
    
end subroutine cs_profile

! =======================================================================
! =======================================================================

subroutine cs_limiters (km, a4)
    
    implicit none
    
    integer, intent (in) :: km
    
    real, intent (inout) :: a4 (4, km) ! ppm array
    
    real, parameter :: r12 = 1. / 12.
    
    integer :: k
    
    ! -----------------------------------------------------------------------
    ! positive definite constraint
    ! -----------------------------------------------------------------------
    
    do k = 1, km
        if (abs (a4 (3, k) - a4 (2, k)) .lt. - a4 (4, k)) then
            if ((a4 (1, k) + 0.25 * (a4 (3, k) - a4 (2, k)) ** 2 / a4 (4, k) + a4 (4, k) * r12) .lt. 0.) then
                if (a4 (1, k) .lt. a4 (3, k) .and. a4 (1, k) .lt. a4 (2, k)) then
                    a4 (3, k) = a4 (1, k)
                    a4 (2, k) = a4 (1, k)
                    a4 (4, k) = 0.
                elseif (a4 (3, k) .gt. a4 (2, k)) then
                    a4 (4, k) = 3. * (a4 (2, k) - a4 (1, k))
                    a4 (3, k) = a4 (2, k) - a4 (4, k)
                else
                    a4 (4, k) = 3. * (a4 (3, k) - a4 (1, k))
                    a4 (2, k) = a4 (3, k) - a4 (4, k)
                endif
            endif
        endif
    enddo
    
end subroutine cs_limiters

! =======================================================================
! calculation of vertical fall speed
! =======================================================================

subroutine fall_speed (ks, ke, den, qs, qi, qg, ql, tk, vts, vti, vtg)
    
    implicit none
    
    integer, intent (in) :: ks, ke
    
    real (kind = r_grid), intent (in), dimension (ks:ke) :: tk
    real, intent (in), dimension (ks:ke) :: den, qs, qi, qg, ql
    real, intent (out), dimension (ks:ke) :: vts, vti, vtg
    
    ! fall velocity constants:
    
    real, parameter :: aa = - 4.14122e-5
    real, parameter :: bb = - 0.00538922
    real, parameter :: cc = - 0.0516344
    real, parameter :: dd = 0.00216078
    real, parameter :: ee = 1.9714
    
    real, dimension (ks:ke) :: qden, tc, rhof
    
    real :: vi0
    
    integer :: k
    
    ! -----------------------------------------------------------------------
    ! marshall - palmer formula
    ! -----------------------------------------------------------------------
    
    ! -----------------------------------------------------------------------
    ! try the local air density -- for global model; the true value could be
    ! much smaller than sfcrho over high mountains
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        rhof (k) = sqrt (min (10., sfcrho / den (k)))
    enddo
    
    ! -----------------------------------------------------------------------
    ! ice:
    ! -----------------------------------------------------------------------
    
    if (const_vi) then
        vti (:) = vi_fac
    else
        ! -----------------------------------------------------------------------
        ! use deng and mace (2008, grl), which gives smaller fall speed than hd90 formula
        ! -----------------------------------------------------------------------
        vi0 = 0.01 * vi_fac
        do k = ks, ke
            if (qi (k) .lt. qfmin) then ! this is needed as the fall - speed maybe problematic for small qi
                vti (k) = 0.0
            else
                tc (k) = tk (k) - tice
                if (ifflag .eq. 1) then
                    vti (k) = (3. + log10 (qi (k) * den (k))) * (tc (k) * (aa * tc (k) + bb) + cc) + dd * tc (k) + ee
                    vti (k) = vi0 * exp (log (10.) * vti (k))
                endif
                if (ifflag .eq. 2) &
                    vti (k) = vi_fac * 3.29 * (qi (k) * den (k)) ** 0.16
                vti (k) = min (vi_max, max (0.0, vti (k)))
            endif
        enddo
    endif
    
    ! -----------------------------------------------------------------------
    ! snow:
    ! -----------------------------------------------------------------------
    
    if (const_vs) then
        vts (:) = vs_fac ! 1. ifs_2016
    else
        do k = ks, ke
            if (qs (k) .lt. qfmin) then
                vts (k) = 0.0
            else
                vts (k) = vs_fac * vcons * rhof (k) * exp (0.0625 * log (qs (k) * den (k) / norms))
                vts (k) = min (vs_max, max (0.0, vts (k)))
            endif
        enddo
    endif
    
    ! -----------------------------------------------------------------------
    ! graupel:
    ! -----------------------------------------------------------------------
    
    if (const_vg) then
        vtg (:) = vg_fac ! 2.
    else
        if (do_hail) then
            do k = ks, ke
                if (qg (k) .lt. qfmin) then
                    vtg (k) = 0.0
                else
                    vtg (k) = vg_fac * vconh * rhof (k) * sqrt (sqrt (sqrt (qg (k) * den (k) / normh))) / sqrt (den (k))
                    vtg (k) = min (vg_max, max (0.0, vtg (k)))
                endif
            enddo
        else
            do k = ks, ke
                if (qg (k) .lt. qfmin) then
                    vtg (k) = 0.0
                else
                    vtg (k) = vg_fac * vcong * rhof (k) * sqrt (sqrt (sqrt (qg (k) * den (k) / normg))) / sqrt (den (k))
                    vtg (k) = min (vg_max, max (0.0, vtg (k)))
                endif
            enddo
        endif
    endif
    
end subroutine fall_speed

! =======================================================================
! fast saturation adjustments
! this is designed for single - moment 6 - class cloud microphysics schemes
! handles the heat release due to in situ phase changes.
! it mainly consists of melting / freezing, condensation / evaporation,
! sublimation / deposition, and autoconversion processes.
! =======================================================================

subroutine fast_sat_adj (mdt, is, ie, js, je, ng, hydrostatic, consv_te, &
        te, qv, ql, qi, qr, qs, qg, qa, qnl, qni, hs, dpln, delz, pt, delp, &
        q_con, cappa, gsize, dtdt, out_dt, last_step)
    
    implicit none
    
    logical, intent (in) :: hydrostatic, consv_te, out_dt, last_step
    
    integer, intent (in) :: is, ie, js, je, ng
    
    real, intent (in) :: mdt
    
    real, intent (in), dimension (is - ng:ie + ng, js - ng:je + ng) :: delp, hs
    real, intent (in), dimension (is:ie, js:je) :: dpln
    real, intent (in), dimension (is:ie, js:je) :: delz
    
    real (kind = r_grid), intent (in), dimension (is:ie, js:je) :: gsize
    
    real, intent (inout), dimension (is - ng:ie + ng, js - ng:je + ng) :: pt, qv, ql, qr
    real, intent (inout), dimension (is - ng:ie + ng, js - ng:je + ng) :: qi, qs, qg
    real, intent (inout), dimension (is - ng:, js - ng:) :: q_con, cappa
    real, intent (inout), dimension (is:ie, js:je) :: dtdt
    
    real, intent (inout), dimension (is - ng:ie + ng, js - ng:je + ng) :: qa, te, qnl, qni
    
    real (kind = r_grid), dimension (is:ie, js:je) :: te_beg, te_end, tw_beg, tw_end
    
    real (kind = r_grid), dimension (is:ie) :: pt1
    real, dimension (is:ie) :: wqsat, dq2dt, qpz, cvm, t0, qstar
    real, dimension (is:ie) :: icp2, lcp2, tcp2, tcp3
    real, dimension (is:ie) :: den, q_liq, q_sol, q_cond, src, sink, hvar
    real, dimension (is:ie) :: mc_air, lhl, lhi, ccn, cin
    
    real :: d0_vap ! the same as dc_vap, except that cp_vap can be cp_vap or cv_vap
    real :: lv00 ! the same as lv0, except that cp_vap can be cp_vap or cv_vap
    real :: li00
    real :: qsw, rh, ccn0
    real :: tc, qsi, dqsdt, dq, dq0, pidep, qi_gen, qi_crt, tmp, dtmp
    real :: rqi, q_plus, q_minus
    real :: tin, sdt, dt_bigg, adj_fac
    real :: fac_smlt, fac_r2g, fac_i2s, fac_imlt, fac_l2r, fac_v2l, fac_l2v
    real :: factor, qim, c_air, c_vap, dw
    
    integer :: i, j
    
    sdt = 0.5 * mdt
    dt_bigg = mdt
    
    ! -----------------------------------------------------------------------
    ! conversion scalar / factor
    ! -----------------------------------------------------------------------
    
    fac_i2s = 1. - exp (- mdt / tau_i2s)
    fac_r2g = 1. - exp (- mdt / tau_r2g)
    fac_l2r = 1. - exp (- mdt / tau_l2r)
    fac_v2l = 1. - exp (- sdt / tau_v2l)
    
    fac_l2v = 1. - exp (- sdt / tau_l2v)
    fac_l2v = min (sat_adj0, fac_l2v)
    
    fac_imlt = 1. - exp (- mdt / tau_imlt)
    fac_smlt = 1. - exp (- mdt / tau_smlt)
    
    ! -----------------------------------------------------------------------
    ! heat capacity of dry air and water vapor based on hydrostatical property
    ! -----------------------------------------------------------------------
    
    if (hydrostatic) then
        c_air = cp_air
        c_vap = cp_vap
    else
        c_air = cv_air
        c_vap = cv_vap
    endif
    d0_vap = c_vap - c_liq
    lv00 = hlv - d0_vap * tice
    li00 = hlf - dc_ice * tice
    
    do j = js, je
        
        ! -----------------------------------------------------------------------
        ! compute true temperature
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            q_liq (i) = ql (i, j) + qr (i, j)
            q_sol (i) = qi (i, j) + qs (i, j) + qg (i, j)
            qpz (i) = q_liq (i) + q_sol (i)
#ifdef MOIST_CAPPA
            pt1 (i) = pt (i, j) / ((1 + zvir * qv (i, j)) * (1 - qpz (i)))
#else
            pt1 (i) = pt (i, j) / (1 + zvir * qv (i, j))
#endif
            t0 (i) = pt1 (i)
            qpz (i) = qpz (i) + qv (i, j)
        enddo
        
        ! -----------------------------------------------------------------------
        ! moist air density based on hydrostatical property
        ! -----------------------------------------------------------------------
        
        if (hydrostatic) then
            do i = is, ie
                den (i) = delp (i, j) / (dpln (i, j) * rdgas * pt (i, j))
            enddo
        else
            do i = is, ie
                den (i) = - delp (i, j) / (grav * delz (i, j))
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! calculate cloud condensation nuclei (ccn)
        ! the following is based on klein eq. 15
        ! -----------------------------------------------------------------------
        
        if (prog_ccn) then
            do i = is, ie
                ccn (i) = max (10.0, qnl (i, j)) * 1.e6
                cin (i) = max (10.0, qni (i, j)) * 1.e6
                ccn (i) = ccn (i) / den (i)
            enddo
        else
            do i = is, ie
                ccn0 = (ccn_l * min (1., abs (hs (i, j)) / (10. * grav)) + &
                    ccn_o * (1. - min (1., abs (hs (i, j)) / (10. * grav)))) * 1.e6
                ccn (i) = ccn0 / den (i)
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! moist heat capacity and latent heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            mc_air (i) = (1. - qpz (i)) * c_air
            cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! for energy fixer
        ! -----------------------------------------------------------------------
        
        if (consv_te) then
            if (hydrostatic) then
                do i = is, ie
                    te (i, j) = - c_air * t0 (i)
                enddo
            else
                do i = is, ie
#ifdef MOIST_CAPPA
                    te (i, j) = - cvm (i) * t0 (i)
#else
                    te (i, j) = - c_air * t0 (i)
#endif
                enddo
            endif
        endif
        
        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do i = is, ie
                te_beg (i, j) = cvm (i) * pt1 (i) + lv00 * qv (i, j) - li00 * q_sol (i)
                te_beg (i, j) = rgrav * te_beg (i, j) * delp (i, j) * gsize (i, j) ** 2.0
                tw_beg (i, j) = rgrav * (qv (i, j) + q_liq (i) + q_sol (i)) * delp (i, j) * gsize (i, j) ** 2.0
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! fix negative cloud ice with snow
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (qi (i, j) .lt. 0.) then
                qs (i, j) = qs (i, j) + qi (i, j)
                qi (i, j) = 0.
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! melting of cloud ice to cloud water and rain
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (qi (i, j) .gt. qcmin .and. pt1 (i) .gt. tice) then
                sink (i) = min (qi (i, j), fac_imlt * (pt1 (i) - tice) / icp2 (i))
                qi (i, j) = qi (i, j) - sink (i)
                tmp = min (sink (i), dim (ql_mlt, ql (i, j)))
                ql (i, j) = ql (i, j) + tmp
                qr (i, j) = qr (i, j) + sink (i) - tmp
                q_liq (i) = q_liq (i) + sink (i)
                q_sol (i) = q_sol (i) - sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) - sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latent heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! fix negative snow with graupel or graupel with available snow
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (qs (i, j) .lt. 0.) then
                qg (i, j) = qg (i, j) + qs (i, j)
                qs (i, j) = 0.
            elseif (qg (i, j) .lt. 0.) then
                tmp = min (- qg (i, j), max (0., qs (i, j)))
                qg (i, j) = qg (i, j) + tmp
                qs (i, j) = qs (i, j) - tmp
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! fix negative cloud water with rain or rain with available cloud water
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (ql (i, j) .lt. 0.) then
                tmp = min (- ql (i, j), max (0., qr (i, j)))
                ql (i, j) = ql (i, j) + tmp
                qr (i, j) = qr (i, j) - tmp
            elseif (qr (i, j) .lt. 0.) then
                tmp = min (- qr (i, j), max (0., ql (i, j)))
                ql (i, j) = ql (i, j) - tmp
                qr (i, j) = qr (i, j) + tmp
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! enforce complete freezing of cloud water to cloud ice below - 48 c
        ! it can be - 50 c, straka, 2009
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            dtmp = tice - 48. - pt1 (i)
            if (ql (i, j) .gt. qcmin .and. dtmp .gt. 0.) then
                sink (i) = min (ql (i, j), dtmp / icp2 (i))
                ql (i, j) = ql (i, j) - sink (i)
                qi (i, j) = qi (i, j) + sink (i)
                q_liq (i) = q_liq (i) - sink (i)
                q_sol (i) = q_sol (i) + sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) + sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latent heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhl (i) = lv00 + d0_vap * pt1 (i)
            lhi (i) = li00 + dc_ice * pt1 (i)
            lcp2 (i) = lhl (i) / cvm (i)
            icp2 (i) = lhi (i) / cvm (i)
            tcp3 (i) = lcp2 (i) + icp2 (i) * min (1., dim (tice, pt1 (i)) / 48.)
        enddo
        
        ! -----------------------------------------------------------------------
        ! condensation / evaporation between water vapor and cloud water
        ! -----------------------------------------------------------------------
        
        adj_fac = sat_adj0
        do i = is, ie
            tin = pt1 (i)
            wqsat (i) = wqs (tin, den (i), dq2dt (i))
            dq0 = (qv (i, j) - wqsat (i)) / (1. + tcp3 (i) * dq2dt (i))
            if (dq0 .gt. 0.) then
                src (i) = min (adj_fac * dq0, max (ql_gen - ql (i, j), fac_v2l * dq0))
            else
                ! sjl, 20170703
                ! factor = - min (1., fac_l2v * sqrt (max (0., ql (i, j)) / 1.e-5) * &
                ! 10. * (1. - qv (i, j) / wqsat (i)))
                ! factor = - fac_l2v
                ! factor = - 1
                factor = - min (1., fac_l2v * 10. * (1. - qv (i, j) / wqsat (i)))
                src (i) = - min (ql (i, j), factor * dq0)
            endif
            qv (i, j) = qv (i, j) - src (i)
            ql (i, j) = ql (i, j) + src (i)
            q_liq (i) = q_liq (i) + src (i)
            cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
            pt1 (i) = pt1 (i) + src (i) * lhl (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latent heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhl (i) = lv00 + d0_vap * pt1 (i)
            lhi (i) = li00 + dc_ice * pt1 (i)
            lcp2 (i) = lhl (i) / cvm (i)
            icp2 (i) = lhi (i) / cvm (i)
            tcp3 (i) = lcp2 (i) + icp2 (i) * min (1., dim (tice, pt1 (i)) / 48.)
        enddo
        
        if (last_step) then
            
            ! -----------------------------------------------------------------------
            ! condensation / evaporation between water vapor and cloud water at last time step
            ! enforce upper (no super_sat) & lower (critical rh) bounds
            ! final iteration:
            ! -----------------------------------------------------------------------
            
            do i = is, ie
                tin = pt1 (i)
                wqsat (i) = wqs (tin, den (i), dq2dt (i))
                dq0 = (qv (i, j) - wqsat (i)) / (1. + tcp3 (i) * dq2dt (i))
                if (dq0 .gt. 0.) then
                    src (i) = dq0
                else
                    ! sjl, 20170703
                    ! factor = - min (1., fac_l2v * sqrt (max (0., ql (i, j)) / 1.e-5) * &
                    ! 10. * (1. - qv (i, j) / wqsat (i)))
                    ! factor = - fac_l2v
                    ! factor = - 1
                    factor = - min (1., fac_l2v * 10. * (1. - qv (i, j) / wqsat (i)))
                    src (i) = - min (ql (i, j), factor * dq0)
                endif
                adj_fac = 1.
                qv (i, j) = qv (i, j) - src (i)
                ql (i, j) = ql (i, j) + src (i)
                q_liq (i) = q_liq (i) + src (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) + src (i) * lhl (i) / cvm (i)
            enddo
            
            ! -----------------------------------------------------------------------
            ! update latent heat coefficient
            ! -----------------------------------------------------------------------
            
            do i = is, ie
                lhl (i) = lv00 + d0_vap * pt1 (i)
                lhi (i) = li00 + dc_ice * pt1 (i)
                lcp2 (i) = lhl (i) / cvm (i)
                icp2 (i) = lhi (i) / cvm (i)
            enddo
            
        endif
        
        ! -----------------------------------------------------------------------
        ! homogeneous freezing of cloud water to cloud ice, - 40 c to - 48 c
        ! it can be - 50 c, straka, 2009
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            dtmp = t_wfr - pt1 (i)
            if (ql (i, j) .gt. qcmin .and. dtmp .gt. 0.) then
                sink (i) = min (ql (i, j), ql (i, j) * dtmp * 0.125, dtmp / icp2 (i))
                ql (i, j) = ql (i, j) - sink (i)
                qi (i, j) = qi (i, j) + sink (i)
                q_liq (i) = q_liq (i) - sink (i)
                q_sol (i) = q_sol (i) + sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) + sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latent heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! bigg mechanism (heterogeneous freezing of cloud water to cloud ice)
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            tc = tice - pt1 (i)
            if (ql (i, j) .gt. qcmin .and. tc .gt. 0.) then
                sink (i) = 100. / (rhow * ccn (i)) * dt_bigg * (exp (0.66 * tc) - 1.) * ql (i, j) ** 2
                sink (i) = min (ql (i, j), tc / icp2 (i), sink (i))
                ql (i, j) = ql (i, j) - sink (i)
                qi (i, j) = qi (i, j) + sink (i)
                q_liq (i) = q_liq (i) - sink (i)
                q_sol (i) = q_sol (i) + sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) + sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latent heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! freezing of rain to graupel, complete freezing below - 40 c
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            dtmp = (tice - 0.1) - pt1 (i)
            if (qr (i, j) .gt. qcmin .and. dtmp .gt. 0.) then
                tmp = min (1., (dtmp * 0.025) ** 2) * qr (i, j)
                sink (i) = min (tmp, fac_r2g * dtmp / icp2 (i))
                qr (i, j) = qr (i, j) - sink (i)
                qg (i, j) = qg (i, j) + sink (i)
                q_liq (i) = q_liq (i) - sink (i)
                q_sol (i) = q_sol (i) + sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) + sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latent heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! melting of snow to rain or cloud water, complete melting above 10 c
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            dtmp = pt1 (i) - (tice + 0.1)
            if (qs (i, j) .gt. qcmin .and. dtmp .gt. 0.) then
                tmp = min (1., (dtmp * 0.1) ** 2) * qs (i, j)
                sink (i) = min (tmp, fac_smlt * dtmp / icp2 (i))
                tmp = min (sink (i), dim (qs_mlt, ql (i, j)))
                qs (i, j) = qs (i, j) - sink (i)
                ql (i, j) = ql (i, j) + tmp
                qr (i, j) = qr (i, j) + sink (i) - tmp
                ! ljz, 20190716
                ! qr (i, j) = qr (i, j) + sink (i)
                q_liq (i) = q_liq (i) + sink (i)
                q_sol (i) = q_sol (i) - sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) - sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! autoconversion from cloud water to rain
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (ql (i, j) .gt. ql0_max) then
                sink (i) = fac_l2r * (ql (i, j) - ql0_max)
                qr (i, j) = qr (i, j) + sink (i)
                ql (i, j) = ql (i, j) - sink (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latent heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhl (i) = lv00 + d0_vap * pt1 (i)
            lhi (i) = li00 + dc_ice * pt1 (i)
            lcp2 (i) = lhl (i) / cvm (i)
            icp2 (i) = lhi (i) / cvm (i)
            tcp2 (i) = lcp2 (i) + icp2 (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! sublimation / deposition between water vapor and cloud ice
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            src (i) = 0.
            if (pt1 (i) .lt. t_sub) then
                src (i) = dim (qv (i, j), qcmin)
            elseif (pt1 (i) .lt. tice) then
                tin = pt1 (i)
                qsi = iqs (tin, den (i), dqsdt)
                dq = qv (i, j) - qsi
                sink (i) = adj_fac * dq / (1. + tcp2 (i) * dqsdt)
                if (qi (i, j) .gt. qcmin) then
                    if (.not. prog_ccn) then
                        if (inflag .eq. 1) &
                            ! hong et al., 2004
                            cin (i) = 5.38e7 * exp (0.75 * log (qi (i, j) * den (i)))
                        if (inflag .eq. 2) &
                            ! meyers et al., 1992
                            cin (i) = exp (-2.80 + 0.262 * (tice - pt1 (i))) * 1000.0 ! convert from L^-1 to m^-3
                        if (inflag .eq. 3) &
                            ! meyers et al., 1992
                            cin (i) = exp (-0.639 + 12.96 * (qv (i, j) / qsi - 1.0)) * 1000.0 ! convert from L^-1 to m^-3
                        if (inflag .eq. 4) &
                            ! cooper, 1986
                            cin (i) = 5.e-3 * exp (0.304 * (tice - pt1 (i))) * 1000.0 ! convert from L^-1 to m^-3
                        if (inflag .eq. 5) &
                            ! flecther, 1962
                            cin (i) = 1.e-5 * exp (0.5 * (tice - pt1 (i))) * 1000.0 ! convert from L^-1 to m^-3
                    endif
                    pidep = sdt * dq * 4.0 * 11.9 * exp (0.5 * log (qi (i, j) * den (i) * cin (i))) &
                         / (qsi * den (i) * (hlv + hlf) ** 2 / (0.0243 * rvgas * pt1 (i) ** 2) + 4.42478e4)
                else
                    pidep = 0.
                endif
                if (dq .gt. 0.) then
                    tmp = tice - pt1 (i)
                    ! -----------------------------------------------------------------------
                    ! WRF / WSM6 scheme: qi_gen = 4.92e-11 * (1.e3 * exp (0.1 * tmp)) ** 1.33
                    ! optimized: qi_gen = 4.92e-11 * exp (1.33 * log (1.e3 * exp (0.1 * tmp)))
                    ! qi_gen ~ 4.808e-7 at 0 c; 1.818e-6 at - 10 c, 9.827e-5 at - 40 c
                    ! the following value is constructed such that qi_crt = 0 at 0 c and at - 10 c matches
                    ! WRF / WSM6 ice initiation scheme; qi_crt = qi_gen * min (qi_lim, 0.1 * tmp) / den
                    ! -----------------------------------------------------------------------
                    qi_gen = 4.92e-11 * exp (1.33 * log (1.e3 * exp (0.1 * tmp)))
                    if (igflag .eq. 1) &
                        qi_crt = qi_gen / den (i)
                    if (igflag .eq. 2) &
                        qi_crt = qi_gen * min (qi_lim, 0.1 * tmp) / den (i)
                    if (igflag .eq. 3) &
                        qi_crt = 1.82e-6 * min (qi_lim, 0.1 * tmp) / den (i)
                    if (igflag .eq. 4) &
                        qi_crt = max (qi_gen, 1.82e-6) * min (qi_lim, 0.1 * tmp) / den (i)
                    src (i) = min (sink (i), max (qi_crt - qi (i, j), pidep), tmp / tcp2 (i))
                else
                    pidep = pidep * min (1., dim (pt1 (i), t_sub) * 0.2)
                    src (i) = max (pidep, sink (i), - qi (i, j))
                endif
            endif
            qv (i, j) = qv (i, j) - src (i)
            qi (i, j) = qi (i, j) + src (i)
            q_sol (i) = q_sol (i) + src (i)
            cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
            pt1 (i) = pt1 (i) + src (i) * (lhl (i) + lhi (i)) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! fix negative graupel with available cloud ice
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (qg (i, j) .lt. 0.) then
                tmp = min (- qg (i, j), max (0., qi (i, j)))
                qg (i, j) = qg (i, j) + tmp
                qi (i, j) = qi (i, j) - tmp
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! autoconversion from cloud ice to snow
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            qim = qi0_max / den (i)
            if (qi (i, j) .gt. qim) then
                sink (i) = fac_i2s * (qi (i, j) - qim)
                qi (i, j) = qi (i, j) - sink (i)
                qs (i, j) = qs (i, j) + sink (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            do i = is, ie
                te_end (i, j) = cvm (i) * pt1 (i) + lv00 * qv (i, j) - li00 * q_sol (i)
                te_end (i, j) = rgrav * te_end (i, j) * delp (i, j) * gsize (i, j) ** 2.0
                tw_end (i, j) = rgrav * (qv (i, j) + q_liq (i) + q_sol (i)) * delp (i, j) * gsize (i, j) ** 2.0
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! update virtual temperature
        ! -----------------------------------------------------------------------
        
        do i = is, ie
#ifdef MOIST_CAPPA
            q_con (i, j) = q_liq (i) + q_sol (i)
            tmp = 1. + zvir * qv (i, j)
            pt (i, j) = pt1 (i) * tmp * (1. - q_con (i, j))
            tmp = rdgas * tmp
            cappa (i, j) = tmp / (tmp + cvm (i))
#else
            pt (i, j) = pt1 (i) * (1. + zvir * qv (i, j))
#endif
        enddo
        
        if (out_dt) then
            do i = is, ie
                dtdt (i, j) = dtdt (i, j) + pt1 (i) - t0 (i)
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! for energy fixer
        ! -----------------------------------------------------------------------
        
        if (consv_te) then
            do i = is, ie
                if (hydrostatic) then
                    te (i, j) = delp (i, j) * (te (i, j) + c_air * pt1 (i))
                else
#ifdef MOIST_CAPPA
                    te (i, j) = delp (i, j) * (te (i, j) + cvm (i) * pt1 (i))
#else
                    te (i, j) = delp (i, j) * (te (i, j) + c_air * pt1 (i))
#endif
                endif
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! update latent heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhl (i) = lv00 + d0_vap * pt1 (i)
            lhi (i) = li00 + dc_ice * pt1 (i)
            cvm (i) = mc_air (i) + (qv (i, j) + q_liq (i) + q_sol (i)) * c_vap
            lcp2 (i) = lhl (i) / cvm (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! compute cloud fraction
        ! -----------------------------------------------------------------------
        
        if (do_qa .and. last_step) then
            
            ! -----------------------------------------------------------------------
            ! combine water species
            ! -----------------------------------------------------------------------
            
            if (rad_snow) then
                if (rad_graupel) then
                    do i = is, ie
                        q_sol (i) = qi (i, j) + qs (i, j) + qg (i, j)
                    enddo
                else
                    do i = is, ie
                        q_sol (i) = qi (i, j) + qs (i, j)
                    enddo
                endif
            else
                do i = is, ie
                    q_sol (i) = qi (i, j)
                enddo
            endif
            if (rad_rain) then
                do i = is, ie
                    q_liq (i) = ql (i, j) + qr (i, j)
                enddo
            else
                do i = is, ie
                    q_liq (i) = ql (i, j)
                enddo
            endif
            do i = is, ie
                q_cond (i) = q_sol (i) + q_liq (i)
            enddo
            
            ! -----------------------------------------------------------------------
            ! use the "liquid - frozen water temperature" (tin) to compute saturated
            ! specific humidity
            ! -----------------------------------------------------------------------
            
            do i = is, ie
                
                tin = pt1 (i) - (lcp2 (i) * q_cond (i) + icp2 (i) * q_sol (i))
                
                ! -----------------------------------------------------------------------
                ! compute saturated specific humidity
                ! -----------------------------------------------------------------------
                
                if (tin .le. t_wfr) then
                    qstar (i) = iqs (tin, den (i))
                elseif (tin .ge. tice) then
                    qstar (i) = wqs (tin, den (i))
                else
                    qsi = iqs (tin, den (i))
                    qsw = wqs (tin, den (i))
                    if (q_cond (i) .gt. qcmin) then
                        rqi = q_sol (i) / q_cond (i)
                    else
                        rqi = ((tice - tin) / (tice - t_wfr))
                    endif
                    qstar (i) = rqi * qsi + (1. - rqi) * qsw
                endif
                
                ! -----------------------------------------------------------------------
                ! compute sub - grid variability
                ! -----------------------------------------------------------------------
                
                dw = dw_ocean + (dw_land - dw_ocean) * min (1., abs (hs (i, j)) / (10. * grav))
                hvar (i) = min (0.2, max (0.01, dw * sqrt (gsize (i, j) / 100.e3)))
                
                ! -----------------------------------------------------------------------
                ! partial cloudiness by pdf:
                ! assuming subgrid linear distribution in horizontal;
                ! this is effectively a smoother for the binary cloud scheme;
                ! qa = 0.5 if qstar == qpz;
                ! -----------------------------------------------------------------------
                
                rh = qpz (i) / qstar (i)
                
                if (rh .gt. 0.75 .and. qpz (i) .gt. qcmin) then
                    dq = hvar (i) * qpz (i)
                    q_plus = qpz (i) + dq
                    q_minus = qpz (i) - dq
                    if (icloud_f .eq. 2) then
                        if (qpz (i) .gt. qstar (i)) then
                            qa (i, j) = 1.
                        elseif (qstar (i) .lt. q_plus .and. q_cond (i) .gt. qcmin) then
                            qa (i, j) = ((q_plus - qstar (i)) / dq) ** 2
                            qa (i, j) = min (1., qa (i, j))
                        else
                            qa (i, j) = 0.
                        endif
                    else
                        if (qstar (i) .lt. q_minus) then
                            qa (i, j) = 1.
                        else
                            if (qstar (i) .lt. q_plus) then
                                if (icloud_f .eq. 0) then
                                    qa (i, j) = (q_plus - qstar (i)) / (dq + dq)
                                else
                                    qa (i, j) = (q_plus - qstar (i)) / &
                                         (2. * dq * (1. - q_cond (i)))
                                endif
                            else
                                qa (i, j) = 0.
                            endif
                            if (q_cond (i) .gt. qcmin) then
                                qa (i, j) = max (cld_min, qa (i, j))
                            endif
                            qa (i, j) = min (1., qa (i, j))
                        endif
                    endif
                else
                    qa (i, j) = 0.
                endif
                
            enddo
            
        endif
        
    enddo
    
    ! -----------------------------------------------------------------------
    ! total energy checker
    ! -----------------------------------------------------------------------
    
    if (consv_checker) then
        if (abs (sum (te_end) - sum (te_beg)) / sum (te_beg) .gt. te_err) then
            print *, "fast_sat_adj te: ", sum (te_beg) / sum (gsize ** 2.0), &
                sum (te_end) / sum (gsize ** 2.0), &
                 (sum (te_end) - sum (te_beg)) / sum (te_beg)
        endif
        if (abs (sum (tw_end) - sum (tw_beg)) / sum (tw_beg) .gt. te_err) then
            print *, "fast_sat_adj tw: ", sum (tw_beg) / sum (gsize ** 2.0), &
                sum (tw_end) / sum (gsize ** 2.0), &
                 (sum (tw_end) - sum (tw_beg)) / sum (tw_beg)
        endif
    endif
    
end subroutine fast_sat_adj

! =======================================================================
! cloud radii diagnosis built for gfdl cloud microphysics
! radius of cloud species diagnosis
! =======================================================================

subroutine cld_eff_rad (is, ie, ks, ke, lsm, p, delp, t, qw, qi, qr, qs, qg, &
        qcw, qci, qcr, qcs, qcg, rew, rei, rer, res, reg, &
        cld, cloud, snowd, cnvw, cnvi, cnvc)
    
    implicit none
    
    integer, intent (in) :: is, ie
    integer, intent (in) :: ks, ke
    
    real, intent (in), dimension (is:ie) :: lsm ! land sea mask, 0: ocean, 1: land, 2: sea ice
    real, intent (in), dimension (is:ie) :: snowd ! snow depth (mm)
    
    real, intent (in), dimension (is:ie, ks:ke) :: delp, t, p
    real, intent (in), dimension (is:ie, ks:ke) :: cloud ! cloud fraction
    real, intent (in), dimension (is:ie, ks:ke) :: qw, qi, qr, qs, qg ! mass mixing ratio (kg / kg)
    
    real, intent (in), dimension (is:ie, ks:ke), optional :: cnvw, cnvi ! convective cloud water / ice mass mixing ratio (kg / kg)
    real, intent (in), dimension (is:ie, ks:ke), optional :: cnvc ! convective cloud fraction
    
    real, intent (inout), dimension (is:ie, ks:ke) :: qcw, qci, qcr, qcs, qcg ! units: g / m^2
    real, intent (inout), dimension (is:ie, ks:ke) :: rew, rei, rer, res, reg ! radii (micron)
    real, intent (inout), dimension (is:ie, ks:ke) :: cld ! total cloud fraction
    
    ! local variables
    
    integer :: i, k, ind
    
    real, dimension (is:ie, ks:ke) :: qmw, qmr, qmi, qms, qmg ! mass mixing ratio (kg / kg)
    
    real :: dpg ! dp / g
    real :: rho ! density (kg / m^3)
    real :: ccnw ! cloud condensate nuclei for cloud water (cm^ - 3)
    real :: mask
    real :: cor
    real :: tc0
    real :: bw
    
    real :: lambdar, lambdas, lambdag
    real :: rei_fac
    
    real :: ccno = 90. ! ccn over ocean (cm^ - 3)
    real :: ccnl = 270. ! ccn over land (cm^ - 3)
    
    real :: retab (138) = (/ &
        0.05000, 0.05000, 0.05000, 0.05000, 0.05000, 0.05000, &
        0.05500, 0.06000, 0.07000, 0.08000, 0.09000, 0.10000, &
        0.20000, 0.30000, 0.40000, 0.50000, 0.60000, 0.70000, &
        0.80000, 0.90000, 1.00000, 1.10000, 1.20000, 1.30000, &
        1.40000, 1.50000, 1.60000, 1.80000, 2.00000, 2.20000, &
        2.40000, 2.60000, 2.80000, 3.00000, 3.20000, 3.50000, &
        3.80000, 4.10000, 4.40000, 4.70000, 5.00000, 5.30000, &
        5.60000, 5.92779, 6.26422, 6.61973, 6.99539, 7.39234, &
        7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930, &
        10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319, &
        15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955, &
        20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125, &
        27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943, &
        31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601, &
        34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078, &
        38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635, &
        42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221, &
        50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898, &
        65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833, &
        93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424, &
        124.954, 130.630, 136.457, 142.446, 148.608, 154.956, &
        161.503, 168.262, 175.248, 182.473, 189.952, 197.699, &
        205.728, 214.055, 222.694, 231.661, 240.971, 250.639 /)
    
    qmw = qw
    qmi = qi
    qmr = qr
    qms = qs
    qmg = qg
    cld = cloud
    
    if (present (cnvw)) then
        qmw = qmw + cnvw
    endif
    if (present (cnvi)) then
        qmi = qmi + cnvi
    endif
    if (present (cnvc)) then
        cld = cnvc + (1 - cnvc) * cld
    endif
    
    if (liq_ice_combine) then
        do k = ks, ke
            do i = is, ie
                qmw (i, k) = qmw (i, k) + qmr (i, k)
                qmr (i, k) = 0.0
                qmi (i, k) = qmi (i, k) + qms (i, k) + qmg (i, k)
                qms (i, k) = 0.0
                qmg (i, k) = 0.0
            enddo
        enddo
    endif
    
    if (snow_grauple_combine) then
        do k = ks, ke
            do i = is, ie
                qms (i, k) = qms (i, k) + qmg (i, k)
                qmg (i, k) = 0.0
            enddo
        enddo
    endif
    
    do k = ks, ke
        
        do i = is, ie
            
            qmw (i, k) = max (qmw (i, k), 0.0)
            qmi (i, k) = max (qmi (i, k), 0.0)
            qmr (i, k) = max (qmr (i, k), 0.0)
            qms (i, k) = max (qms (i, k), 0.0)
            qmg (i, k) = max (qmg (i, k), 0.0)
            
            cld (i, k) = min (max (cld (i, k), 0.0), 1.0)
            
            mask = min (max (lsm (i), 0.0), 2.0)
            
            dpg = abs (delp (i, k)) / grav
            ! rho = p (i, k) / (rdgas * t (i, k) * (1. + zvir * qv)) ! needs qv
            rho = p (i, k) / (rdgas * t (i, k))
            ! use rho = dpg / delz ! needs delz
            
            tc0 = t (i, k) - tice
            
            if (rewflag .eq. 1) then
                
                ! -----------------------------------------------------------------------
                ! cloud water (martin et al., 1994)
                ! -----------------------------------------------------------------------
                
#ifndef MARTIN_CCN
                ccnw = ccno * abs (mask - 1.0) + ccnl * (1.0 - abs (mask - 1.0))
#else
                ccnw = 0.80 * (- 1.15e-3 * (ccno ** 2) + 0.963 * ccno + 5.30) * abs (mask - 1.0) + &
                    0.67 * (- 2.10e-4 * (ccnl ** 2) + 0.568 * ccnl - 27.9) * (1.0 - abs (mask - 1.0))
#endif
                
                if (qmw (i, k) .gt. qcmin) then
                    qcw (i, k) = dpg * qmw (i, k) * 1.0e3
                    rew (i, k) = exp (1.0 / 3.0 * log ((3.0 * qmw (i, k) * rho) / (4.0 * pi * rhow * ccnw))) * 1.0e4
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif
                
            endif
            
            if (rewflag .eq. 2) then
                
                ! -----------------------------------------------------------------------
                ! cloud water (martin et al., 1994, gfdl revision)
                ! -----------------------------------------------------------------------
                
                ccnw = 1.077 * ccno * abs (mask - 1.0) + 1.143 * ccnl * (1.0 - abs (mask - 1.0))
                
                if (qmw (i, k) .gt. qcmin) then
                    qcw (i, k) = dpg * qmw (i, k) * 1.0e3
                    rew (i, k) = exp (1.0 / 3.0 * log ((3.0 * qmw (i, k) * rho) / (4.0 * pi * rhow * ccnw))) * 1.0e4
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif
                
            endif
            
            if (rewflag .eq. 3) then
                
                ! -----------------------------------------------------------------------
                ! cloud water (kiehl et al., 1994)
                ! -----------------------------------------------------------------------
                
                if (qmw (i, k) .gt. qcmin) then
                    qcw (i, k) = dpg * qmw (i, k) * 1.0e3
                    rew (i, k) = 14.0 * abs (mask - 1.0) + &
                         (8.0 + (14.0 - 8.0) * min (1.0, max (0.0, - tc0 / 30.0))) * (1.0 - abs (mask - 1.0))
                    rew (i, k) = rew (i, k) + (14.0 - rew (i, k)) * min (1.0, max (0.0, snowd (i) / 1000.0))
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif
                
            endif
            
            if (reiflag .eq. 1) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (heymsfield and mcfarquhar, 1996)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    rei_fac = log (1.0e3 * qmi (i, k) * rho)
                    if (tc0 .lt. - 50) then
                        rei (i, k) = beta / 9.917 * exp (0.109 * rei_fac) * 1.0e3
                    elseif (tc0 .lt. - 40) then
                        rei (i, k) = beta / 9.337 * exp (0.080 * rei_fac) * 1.0e3
                    elseif (tc0 .lt. - 30) then
                        rei (i, k) = beta / 9.208 * exp (0.055 * rei_fac) * 1.0e3
                    else
                        rei (i, k) = beta / 9.387 * exp (0.031 * rei_fac) * 1.0e3
                    endif
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (reiflag .eq. 2) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (donner et al., 1997)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    if (tc0 .le. - 55) then
                        rei (i, k) = 15.41627
                    elseif (tc0 .le. - 50) then
                        rei (i, k) = 16.60895
                    elseif (tc0 .le. - 45) then
                        rei (i, k) = 32.89967
                    elseif (tc0 .le. - 40) then
                        rei (i, k) = 35.29989
                    elseif (tc0 .le. - 35) then
                        rei (i, k) = 55.65818
                    elseif (tc0 .le. - 30) then
                        rei (i, k) = 85.19071
                    elseif (tc0 .le. - 25) then
                        rei (i, k) = 72.35392
                    else
                        rei (i, k) = 92.46298
                    endif
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (reiflag .eq. 3) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (fu, 2007)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    rei (i, k) = 47.05 + tc0 * (0.6624 + 0.001741 * tc0)
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (reiflag .eq. 4) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (kristjansson et al., 2000)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    ind = min (max (int (t (i, k) - 136.0), 44), 138 - 1)
                    cor = t (i, k) - int (t (i, k))
                    rei (i, k) = retab (ind) * (1. - cor) + retab (ind + 1) * cor
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (reiflag .eq. 5) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (wyser, 1998)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    bw = - 2. + 1.e-3 * log10 (rho * qmi (i, k) / 50.e-3) * max (0.0, - tc0) ** 1.5
                    rei (i, k) = 377.4 + bw * (203.3 + bw * (37.91 + 2.3696 * bw))
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            ! -----------------------------------------------------------------------
            ! rain (lin et al., 1983)
            ! -----------------------------------------------------------------------
            
            if (qmr (i, k) .gt. qcmin) then
                qcr (i, k) = dpg * qmr (i, k) * 1.0e3
                lambdar = exp (0.25 * log (normr / qmr (i, k) / rho))
                rer (i, k) = 0.5 * exp (log (gam480 / 6) / 0.80) / lambdar * 1.0e6
                rer (i, k) = max (rermin, min (rermax, rer (i, k)))
            else
                qcr (i, k) = 0.0
                rer (i, k) = rermin
            endif
            
            ! -----------------------------------------------------------------------
            ! snow (lin et al., 1983)
            ! -----------------------------------------------------------------------
            
            if (qms (i, k) .gt. qcmin) then
                qcs (i, k) = dpg * qms (i, k) * 1.0e3
                lambdas = exp (0.25 * log (norms / qms (i, k) / rho))
                res (i, k) = 0.5 * exp (log (gam425 / 6) / 0.25) / lambdas * 1.0e6
                res (i, k) = max (resmin, min (resmax, res (i, k)))
            else
                qcs (i, k) = 0.0
                res (i, k) = resmin
            endif
            
            ! -----------------------------------------------------------------------
            ! graupel (lin et al., 1983)
            ! -----------------------------------------------------------------------
            
            if (qmg (i, k) .gt. qcmin) then
                qcg (i, k) = dpg * qmg (i, k) * 1.0e3
                lambdag = exp (0.25 * log (normg / qmg (i, k) / rho))
                reg (i, k) = 0.5 * exp (log (gam450 / 6) / 0.50) / lambdag * 1.0e6
                reg (i, k) = max (regmin, min (regmax, reg (i, k)))
            else
                qcg (i, k) = 0.0
                reg (i, k) = regmin
            endif
            
        enddo
        
    enddo
    
end subroutine cld_eff_rad

! =======================================================================
! radar reflectivity
! =======================================================================

subroutine rad_ref (is, ie, js, je, isd, ied, jsd, jed, &
        q, pt, delp, peln, delz, dbz, maxdbz, allmax, &
        npz, ncnst, hydrostatic, zvir, in0r, in0s, in0g, iliqskin, do_inline_mp, &
        sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, mp_top)
    
    ! code from mark stoelinga's dbzcalc.f from the rip package.
    ! currently just using values taken directly from that code, which is
    ! consistent for the mm5 reisner - 2 microphysics. from that file:
    
    ! this routine computes equivalent reflectivity factor (in dbz) at
    ! each model grid point. in calculating ze, the rip algorithm makes
    ! assumptions consistent with those made in an early version
    ! (ca. 1996) of the bulk mixed - phase microphysical scheme in the mm5
    ! model (i.e., the scheme known as "resiner - 2") . for each species:
    !
    ! 1. particles are assumed to be spheres of constant density. the
    ! densities of rain drops, snow particles, and graupel particles are
    ! taken to be rho_r = rho_l = 1000 kg m^ - 3, rho_s = 100 kg m^ - 3, and
    ! rho_g = 400 kg m^ - 3, respectively. (l refers to the density of
    ! liquid water.)
    !
    ! 2. the size distribution (in terms of the actual diameter of the
    ! particles, rather than the melted diameter or the equivalent solid
    ! ice sphere diameter) is assumed to follow an exponential
    ! distribution of the form n (d) = n_0 * exp (lambda * d) .
    !
    ! 3. if in0x = 0, the intercept parameter is assumed constant (as in
    ! early reisner - 2), with values of 8x10^6, 2x10^7, and 4x10^6 m^ - 4,
    ! for rain, snow, and graupel, respectively. various choices of
    ! in0x are available (or can be added) . currently, in0x = 1 gives the
    ! variable intercept for each species that is consistent with
    ! thompson, rasmussen, and manning (2004, monthly weather review,
    ! vol. 132, no. 2, pp. 519 - 542.)
    !
    ! 4. if iliqskin = 1, frozen particles that are at a temperature above
    ! freezing are assumed to scatter as a liquid particle.
    !
    ! more information on the derivation of simulated reflectivity in rip
    ! can be found in stoelinga (2005, unpublished write - up) . contact
    ! mark stoelinga (stoeling@atmos.washington.edu) for a copy.
    
    ! 22sep16: modifying to use the gfdl mp parameters. if doing so remember
    ! that the gfdl mp assumes a constant intercept (in0x = .false.)
    ! ferrier - aligo has an option for fixed slope (rather than fixed intercept) .
    ! thompson presumably is an extension of reisner mp.
    
    implicit none
    
    logical, intent (in) :: hydrostatic, in0r, in0s, in0g, iliqskin, do_inline_mp
    
    integer, intent (in) :: is, ie, js, je, isd, ied, jsd, jed
    integer, intent (in) :: npz, ncnst, mp_top
    integer, intent (in) :: sphum, liq_wat, ice_wat, rainwat, snowwat, graupel
    
    real, intent (in), dimension (isd:ied, jsd:jed, npz) :: pt, delp
    real, intent (in), dimension (is:, js:, 1:) :: delz
    real, intent (in), dimension (isd:ied, jsd:jed, npz, ncnst) :: q
    real, intent (in), dimension (is :ie, npz + 1, js:je) :: peln
    real, intent (out), dimension (is :ie, js :je, npz) :: dbz
    real, intent (out), dimension (is :ie, js :je) :: maxdbz
    
    real, intent (in) :: zvir
    real, intent (out) :: allmax
    
    ! constants for variable intercepts
    ! will need to be changed based on mp scheme
    
    real, parameter :: ron = 8.e6
    real, parameter :: ron2 = 1.e10
    real, parameter :: son = 2.e7
    real, parameter :: gon = 5.e7
    real, parameter :: ron_min = 8.e6
    real, parameter :: ron_qr0 = 0.00010
    real, parameter :: ron_delqr0 = 0.25 * ron_qr0
    real, parameter :: ron_const1r = (ron2 - ron_min) * 0.5
    real, parameter :: ron_const2r = (ron2 + ron_min) * 0.5
    
    ! other constants
    
    real, parameter :: gamma_seven = 720.
    real, parameter :: alpha = 0.224
    real (kind = r_grid), parameter :: factor_s = gamma_seven * 1.e18 * (1. / (pi * rhos)) ** 1.75 &
         * (rhos / rhor) ** 2 * alpha
    
    ! double precision
    
    real (kind = r_grid), dimension (is:ie) :: rhoair, denfac, z_e
    real (kind = r_grid) :: qr1, qs1, qg1, t1, t2, t3, rwat, vtr, vtg, vts
    real (kind = r_grid) :: factorb_s, factorb_g
    real (kind = r_grid) :: temp_c, pres, sonv, gonv, ronv
    
    real :: rhogh, vcongh, normgh
    
    integer :: i, j, k
    
    if (rainwat .lt. 1) return
    
    dbz (:, :, 1:mp_top) = - 20.
    maxdbz (:, :) = - 20. ! minimum value
    allmax = - 20.
    
    if (do_hail .and. .not. do_inline_mp) then
        rhogh = rhoh
        vcongh = vconh
        normgh = normh
    else
        rhogh = rhog
        vcongh = vcong
        normgh = normg
    endif
    
    !$omp parallel do default (shared) private (rhoair, t1, t2, t3, denfac, vtr, vtg, vts, z_e)
    do k = mp_top + 1, npz
        do j = js, je
            if (hydrostatic) then
                do i = is, ie
                    rhoair (i) = delp (i, j, k) / ((peln (i, k + 1, j) - peln (i, k, j)) * &
                        rdgas * pt (i, j, k) * (1. + zvir * q (i, j, k, sphum)))
                    denfac (i) = sqrt (min (10., 1.2 / rhoair (i)))
                    z_e (i) = 0.
                enddo
            else
                do i = is, ie
                    rhoair (i) = - delp (i, j, k) / (grav * delz (i, j, k)) ! moist air density
                    denfac (i) = sqrt (min (10., 1.2 / rhoair (i)))
                    z_e (i) = 0.
                enddo
            endif
            if (rainwat .gt. 0) then
                do i = is, ie
                    ! the following form vectorizes better & more consistent with gfdl_mp
                    ! sjl notes: marshall - palmer, dbz = 200 * precip ** 1.6, precip = 3.6e6 * t1 / rhor * vtr ! [mm / hr]
                    ! gfdl_mp terminal fall speeds are used
                    ! date modified 20170701
                    ! account for excessively high cloud water - > autoconvert (diag only) excess cloud water
                    t1 = rhoair (i) * max (qcmin, q (i, j, k, rainwat) + dim (q (i, j, k, liq_wat), 1.0e-3))
                    vtr = max (1.e-3, vconr * denfac (i) * exp (0.2 * log (t1 / normr)))
                    z_e (i) = 200. * exp (1.6 * log (3.6e6 * t1 / rhor * vtr))
                    ! z_e = 200. * (exp (1.6 * log (3.6e6 * t1 / rhor * vtr)) + &
                    ! exp (1.6 * log (3.6e6 * t3 / rhogh * vtg)) + &
                    ! exp (1.6 * log (3.6e6 * t2 / rhos * vts)))
                enddo
            endif
            if (graupel .gt. 0) then
                do i = is, ie
                    t3 = rhoair (i) * max (qcmin, q (i, j, k, graupel))
                    vtg = max (1.e-3, vcongh * denfac (i) * exp (0.125 * log (t3 / normgh))) / sqrt (rhoair (i))
                    z_e (i) = z_e (i) + 200. * exp (1.6 * log (3.6e6 * t3 / rhogh * vtg))
                enddo
            endif
            if (snowwat .gt. 0) then
                do i = is, ie
                    t2 = rhoair (i) * max (qcmin, q (i, j, k, snowwat))
                    ! vts = max (1.e-3, vcons * denfac * exp (0.0625 * log (t2 / norms)))
                    z_e (i) = z_e (i) + (factor_s / alpha) * t2 * exp (0.75 * log (t2 / rnzs))
                    ! z_e = 200. * (exp (1.6 * log (3.6e6 * t1 / rhor * vtr)) + &
                    ! exp (1.6 * log (3.6e6 * t3 / rhogh * vtg)) + &
                    ! exp (1.6 * log (3.6e6 * t2 / rhos * vts)))
                enddo
            endif
            do i = is, ie
                dbz (i, j, k) = 10. * log10 (max (0.01, z_e (i)))
            enddo
        enddo
    enddo
    
    !$omp parallel do default (shared)
    do j = js, je
        do k = mp_top + 1, npz
            do i = is, ie
                maxdbz (i, j) = max (dbz (i, j, k), maxdbz (i, j))
            enddo
        enddo
    enddo
    
    do j = js, je
        do i = is, ie
            allmax = max (maxdbz (i, j), allmax)
        enddo
    enddo
    
end subroutine rad_ref

! =======================================================================
! accretion function (lin et al. 1983)
! =======================================================================

real function acr3d (v1, v2, q1, q2, c, cac, rho)
    
    implicit none
    
    real, intent (in) :: v1, v2, c, rho
    real, intent (in) :: q1, q2 ! mixing ratio!!!
    real, intent (in) :: cac (3)
    
    real :: t1, s1, s2
    
    ! integer :: k
    !
    ! real :: a
    !
    ! a = 0.0
    ! do k = 1, 3
    ! a = a + cac (k) * ((q1 * rho) ** ((7 - k) * 0.25) * (q2 * rho) ** (k * 0.25))
    ! enddo
    ! acr3d = c * abs (v1 - v2) * a / rho
    
    ! optimized
    
    t1 = sqrt (q1 * rho)
    s1 = sqrt (q2 * rho)
    s2 = sqrt (s1) ! s1 = s2 ** 2
    acr3d = c * abs (v1 - v2) * q1 * s2 * (cac (1) * t1 + cac (2) * sqrt (t1) * s2 + cac (3) * s1)
    
end function acr3d

! =======================================================================
! melting of snow function (lin et al. 1983)
! note: psacw and psacr must be calc before smlt is called
! =======================================================================

real function smlt (tc, dqs, qsrho, psacw, psacr, c, rho, rhofac)
    
    implicit none
    
    real, intent (in) :: tc, dqs, qsrho, psacw, psacr, c (5), rho, rhofac
    
    smlt = (c (1) * tc / rho - c (2) * dqs) * (c (3) * sqrt (qsrho) + &
        c (4) * qsrho ** 0.65625 * sqrt (rhofac)) + c (5) * tc * (psacw + psacr)
    
end function smlt

! =======================================================================
! melting of graupel function (lin et al. 1983)
! note: pgacw and pgacr must be calc before gmlt is called
! =======================================================================

real function gmlt (tc, dqs, qgrho, pgacw, pgacr, c, rho)
    
    implicit none
    
    real, intent (in) :: tc, dqs, qgrho, pgacw, pgacr, c (5), rho
    
    gmlt = (c (1) * tc / rho - c (2) * dqs) * (c (3) * sqrt (qgrho) + &
        c (4) * qgrho ** 0.6875 / rho ** 0.25) + c (5) * tc * (pgacw + pgacr)
    
end function gmlt

! =======================================================================
! sedimentation of heat
! =======================================================================

subroutine sedi_heat (ks, ke, dm, m1, dz, tz, qv, ql, qr, qi, qs, qg, cw)
    ! revised with a precise energy conserving form: s. - j. lin, jan 22, 2018
    ! input q fields are dry mixing ratios, and dm is dry air mass
    implicit none
    integer, intent (in) :: ks, ke
    real, intent (in), dimension (ks:ke) :: dm, m1, dz, qv, ql, qr, qi, qs, qg
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (in) :: cw ! heat capacity
    ! local:
    real, dimension (ks:ke) :: dgz, cv0
    integer :: k
    
    ! this is the vectorized loop
    do k = ks + 1, ke
        dgz (k) = - 0.5 * grav * (dz (k - 1) + dz (k))
        cv0 (k) = dm (k) * (cv_air + qv (k) * cv_vap + (qr (k) + ql (k)) * c_liq + &
             (qi (k) + qs (k) + qg (k)) * c_ice) + cw * (m1 (k) - m1 (k - 1))
        ! cvm_new + cw * m1 (k) = cvm_old + cw * m1 (k - 1)
    enddo
    ! -----------------------------------------------------------------------
    ! implicit algorithm: can't be vectorized
    ! needs an inner i - loop for vectorization
    ! -----------------------------------------------------------------------
    ! top layer: cv0 = cvn + cw * m1 (k)
    ! tz (k) = cv0 (k) * tz (k) / (cvn (k) + cw * m1 (k)) = tz (k) -- > no change
    do k = ks + 1, ke
        tz (k) = (cv0 (k) * tz (k) + m1 (k - 1) * (cw * tz (k - 1) + dgz (k))) / (cv0 (k) + cw * m1 (k - 1))
    enddo
    
end subroutine sedi_heat

! =======================================================================
! check if water species is large enough to fall
! =======================================================================

subroutine check_column (ks, ke, q, no_fall)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke

    real, intent (in) :: q (ks:ke)

    logical, intent (out) :: no_fall

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: k
    
    no_fall = .true.
    
    do k = ks, ke
        if (q (k) .gt. qfmin) then
            no_fall = .false.
            exit
        endif
    enddo
    
end subroutine check_column

! =======================================================================
! fix negative water species
! =======================================================================

subroutine neg_adj (ks, ke, pt, dp, qv, ql, qr, qi, qs, qg, cond)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke

    real, intent (in), dimension (ks:ke) :: dp

    real (kind = r_grid), intent (inout), dimension (ks:ke) :: pt

    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg

    real, intent (out) :: cond
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: k
    
    real :: dq, cvm
    
    real, dimension (ks:ke) :: lcpk, icpk
    
    ! -----------------------------------------------------------------------
    ! calculate moist heat capacity and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        cvm = 1. + qv (k) * c1_vap + (qr (k) + ql (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
        lcpk (k) = (lv00 + d1_vap * pt (k)) / cvm
        icpk (k) = (li00 + d1_ice * pt (k)) / cvm
    enddo

    cond = 0
    
    do k = ks, ke
        
        ! -----------------------------------------------------------------------
        ! fix negative solid-phase hydrometeors
        ! -----------------------------------------------------------------------
        
        ! if cloud ice < 0, borrow from snow
        if (qi (k) .lt. 0.) then
            qs (k) = qs (k) + qi (k)
            qi (k) = 0.
        endif

        ! if snow < 0, borrow from graupel
        if (qs (k) .lt. 0.) then
            qg (k) = qg (k) + qs (k)
            qs (k) = 0.
        endif

        ! if graupel < 0, borrow from rain
        if (qg (k) .lt. 0.) then
            qr (k) = qr (k) + qg (k)
            pt (k) = pt (k) - qg (k) * icpk (k)
            qg (k) = 0.
        endif
        
        ! -----------------------------------------------------------------------
        ! fix negative liquid-phase hydrometeors
        ! -----------------------------------------------------------------------
        
        ! if rain < 0, borrow from cloud water
        if (qr (k) .lt. 0.) then
            ql (k) = ql (k) + qr (k)
            qr (k) = 0.
        endif

        ! if cloud water < 0, borrow from water vapor
        if (ql (k) .lt. 0.) then
            cond = cond - ql (k) * dp (k)
            qv (k) = qv (k) + ql (k)
            pt (k) = pt (k) - ql (k) * lcpk (k)
            ql (k) = 0.
        endif
        
    enddo
    
    ! -----------------------------------------------------------------------
    ! fix water vapor
    ! -----------------------------------------------------------------------
    
    ! borrow water vapor from below
    do k = ks, ke - 1
        if (qv (k) .lt. 0.) then
            qv (k + 1) = qv (k + 1) + qv (k) * dp (k) / dp (k + 1)
            qv (k) = 0.
        endif
    enddo
    
    ! borrow water vapor from above
    if (qv (ke) .lt. 0. .and. qv (ke - 1) .gt. 0.) then
        dq = min (- qv (ke) * dp (ke), qv (ke - 1) * dp (ke - 1))
        qv (ke - 1) = qv (ke - 1) - dq / dp (ke - 1)
        qv (ke) = qv (ke) + dq / dp (ke)
    endif
    
end subroutine neg_adj

! =======================================================================
! qsmith table initialization
! prepare saturation water vapor pressure tables
! =======================================================================

subroutine qsmith_init
    
    implicit none
    
    integer, parameter :: length = 2621
    
    integer :: i
    
    if (.not. tables_are_initialized) then
        
        allocate (table1 (length))
        allocate (table2 (length))
        allocate (table3 (length))
        allocate (table4 (length))

        allocate (des1 (length))
        allocate (des2 (length))
        allocate (des3 (length))
        allocate (des4 (length))
        
        call qs_table1 (length)
        call qs_table2 (length)
        call qs_table3 (length)
        call qs_table4 (length)
        
        do i = 1, length - 1
            des1 (i) = max (0., table1 (i + 1) - table1 (i))
            des2 (i) = max (0., table2 (i + 1) - table2 (i))
            des3 (i) = max (0., table3 (i + 1) - table3 (i))
            des4 (i) = max (0., table4 (i + 1) - table4 (i))
        enddo
        des1 (length) = des1 (length - 1)
        des2 (length) = des2 (length - 1)
        des3 (length) = des3 (length - 1)
        des4 (length) = des4 (length - 1)
        
        tables_are_initialized = .true.
        
    endif
    
end subroutine qsmith_init

! =======================================================================
! saturation water vapor pressure table 1
! 3 - phase table, blended between -20 deg C and 0 deg C
! =======================================================================

subroutine qs_table1 (n)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i
    
    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: tmin, tem, esh20
    real (kind = r_grid) :: wice, wh2o, fac0, fac1, fac2
    real (kind = r_grid) :: esupc (200)
    
    tmin = tice - 160.
    
    ! -----------------------------------------------------------------------
    ! compute es over ice between - 160 deg C and 0 deg C
    ! -----------------------------------------------------------------------
    
    do i = 1, 1600
        tem = tmin + delt * real (i - 1)
        fac0 = (tem - tice) / (tem * tice)
        fac1 = fac0 * li2
        fac2 = (d2ice * log (tem / tice) + fac1) / rvgas
        table1 (i) = e00 * exp (fac2)
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute es over water between - 20 deg C and 102 deg C
    ! -----------------------------------------------------------------------
    
    do i = 1, 1221
        tem = 253.16 + delt * real (i - 1)
        fac0 = (tem - tice) / (tem * tice)
        fac1 = fac0 * lv0
        fac2 = (dc_vap * log (tem / tice) + fac1) / rvgas
        esh20 = e00 * exp (fac2)
        if (i .le. 200) then
            esupc (i) = esh20
        else
            table1 (i + 1400) = esh20
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! derive blended es over ice and supercooled water between - 20 deg C and 0 deg C
    ! -----------------------------------------------------------------------
    
    do i = 1, 200
        tem = 253.16 + delt * real (i - 1)
        wice = 0.05 * (tice - tem)
        wh2o = 0.05 * (tem - 253.16)
        table1 (i + 1400) = wice * table1 (i + 1400) + wh2o * esupc (i)
    enddo
    
end subroutine qs_table1

! =======================================================================
! saturation water vapor pressure table 2
! 2 - phase table, smoothed around 0 deg C
! =======================================================================

subroutine qs_table2 (n)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i, i0, i1
    
    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: tmin, tem0, tem1, fac0, fac1, fac2
    
    tmin = tice - 160.
    
    do i = 1, n
        tem0 = tmin + delt * real (i - 1)
        fac0 = (tem0 - tice) / (tem0 * tice)
        if (i .le. 1600) then
            ! -----------------------------------------------------------------------
            ! compute es over ice between - 160 deg C and 0 deg C
            ! -----------------------------------------------------------------------
            fac1 = fac0 * li2
            fac2 = (d2ice * log (tem0 / tice) + fac1) / rvgas
        else
            ! -----------------------------------------------------------------------
            ! compute es over water between 0 deg C and 102 deg C
            ! -----------------------------------------------------------------------
            fac1 = fac0 * lv0
            fac2 = (dc_vap * log (tem0 / tice) + fac1) / rvgas
        endif
        table2 (i) = e00 * exp (fac2)
    enddo
    
    ! -----------------------------------------------------------------------
    ! smoother around 0 deg C
    ! -----------------------------------------------------------------------
    
    i0 = 1600
    i1 = 1601
    tem0 = 0.25 * (table2 (i0 - 1) + 2. * table1 (i0) + table2 (i0 + 1))
    tem1 = 0.25 * (table2 (i1 - 1) + 2. * table1 (i1) + table2 (i1 + 1))
    table2 (i0) = tem0
    table2 (i1) = tem1
    
end subroutine qs_table2

! =======================================================================
! saturation water vapor pressure table 3
! 2 - phase table, smoothed around - 2 deg C
! see smithsonian meteorological tables page 350
! =======================================================================

subroutine qs_table3 (n)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i, i0, i1
    
    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: esbasw, tbasw, esbasi, tmin, tem, a, b, c, d, e
    real (kind = r_grid) :: tem0, tem1
    
    esbasw = 1013246.0
    tbasw = tice + 100.
    esbasi = 6107.1
    tmin = tice - 160.
    
    do i = 1, n
        tem = tmin + delt * real (i - 1)
        if (i .le. 1580) then ! change from 0 deg C to - 2 deg C
            ! -----------------------------------------------------------------------
            ! compute es over ice between - 160 deg C and - 2 deg C
            ! -----------------------------------------------------------------------
            a = - 9.09718 * (tice / tem - 1.)
            b = - 3.56654 * log10 (tice / tem)
            c = 0.876793 * (1. - tem / tice)
            e = log10 (esbasi)
            table3 (i) = 0.1 * 10 ** (a + b + c + e)
        else
            ! -----------------------------------------------------------------------
            ! compute es over water between - 2 deg C and 102 deg C
            ! -----------------------------------------------------------------------
            a = - 7.90298 * (tbasw / tem - 1.)
            b = 5.02808 * log10 (tbasw / tem)
            c = - 1.3816e-7 * (10 ** ((1. - tem / tbasw) * 11.344) - 1.)
            d = 8.1328e-3 * (10 ** ((tbasw / tem - 1.) * (- 3.49149)) - 1.)
            e = log10 (esbasw)
            table3 (i) = 0.1 * 10 ** (a + b + c + d + e)
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! smoother around - 2 deg C
    ! -----------------------------------------------------------------------
    
    i0 = 1580
    i1 = 1581
    tem0 = 0.25 * (table3 (i0 - 1) + 2. * table1 (i0) + table3 (i0 + 1))
    tem1 = 0.25 * (table3 (i1 - 1) + 2. * table1 (i1) + table3 (i1 + 1))
    table3 (i0) = tem0
    table3 (i1) = tem1
    
end subroutine qs_table3

! =======================================================================
! saturation water vapor pressure table 4
! 1 - phase table
! =======================================================================

subroutine qs_table4 (n)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i
    
    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: tmin, tem, fac0, fac1, fac2
    
    tmin = tice - 160.
    
    ! -----------------------------------------------------------------------
    ! compute es over water only
    ! -----------------------------------------------------------------------
    
    do i = 1, n
        tem = tmin + delt * real (i - 1)
        fac0 = (tem - tice) / (tem * tice)
        fac1 = fac0 * lv0
        fac2 = (dc_vap * log (tem / tice) + fac1) / rvgas
        table4 (i) = e00 * exp (fac2)
    enddo
    
end subroutine qs_table4

! =======================================================================
! compute the saturated specific humidity and its gradient for table 4
! =======================================================================

real function wqs (ta, den, dqdt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: ta, den

    real, intent (out), optional :: dqdt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: it
    
    real :: es, ap1, tmin
    
    if (.not. tables_are_initialized) call qsmith_init
    
    tmin = tice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table4 (it) + (ap1 - it) * des4 (it)
    wqs = es / (rvgas * ta * den)
    if (present (dqdt)) then
        it = ap1 - 0.5
        ! finite diff, del_t = 0.1:
        dqdt = 10. * (des4 (it) + (ap1 - it) * (des4 (it + 1) - des4 (it))) / (rvgas * ta * den)
    endif
    
end function wqs

! =======================================================================
! compute the saturated specific humidity and its gradient for table 2
! =======================================================================

real function iqs (ta, den, dqdt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: ta, den

    real, intent (out), optional :: dqdt

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: it
    
    real :: es, ap1, tmin

    if (.not. tables_are_initialized) call qsmith_init
    
    tmin = tice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    iqs = es / (rvgas * ta * den)
    if (present (dqdt)) then
        it = ap1 - 0.5
        dqdt = 10. * (des2 (it) + (ap1 - it) * (des2 (it + 1) - des2 (it))) / (rvgas * ta * den)
    endif
    
end function iqs

! =======================================================================
! compute the saturated specific humidity and its gradient for table 4
! it is the same as "wqs", but includes the moist effect
! =======================================================================

real function wqs_moist (ta, pa, qv, dqdt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: ta, pa, qv
    
    real, intent (out), optional :: dqdt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: it
    
    real :: es, ap1, tmin, eps10
    
    if (.not. tables_are_initialized) call qsmith_init
    
    tmin = tice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table4 (it) + (ap1 - it) * des4 (it)
    wqs_moist = eps * es * (1. + zvir * qv) / pa
    if (present (dqdt)) then
        eps10 = 10. * eps
        it = ap1 - 0.5
        dqdt = eps10 * (des4 (it) + (ap1 - it) * (des4 (it + 1) - des4 (it))) * (1. + zvir * qv) / pa
    endif
    
end function wqs_moist

! =======================================================================
! compute the saturated specific humidity and its gradient for table 1
! it is the same as "wqs_moist", but for table 1
! =======================================================================

real function qs_moist (ta, pa, qv, dqdt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: ta, pa, qv
    
    real, intent (out), optional :: dqdt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: it
    
    real :: es, ap1, tmin, eps10
    
    if (.not. tables_are_initialized) call qsmith_init
    
    tmin = tice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table1 (it) + (ap1 - it) * des1 (it)
    qs_moist = eps * es * (1. + zvir * qv) / pa
    if (present (dqdt)) then
        eps10 = 10. * eps
        it = ap1 - 0.5
        dqdt = eps10 * (des1 (it) + (ap1 - it) * (des1 (it + 1) - des1 (it))) * (1. + zvir * qv) / pa
    endif
    
end function qs_moist

! =======================================================================
! compute the saturated specific humidity and its gradient for table 1
! it is the same as "qs_moist", but written as vector function
! =======================================================================

subroutine qsmith (im, km, ks, ta, pa, qv, qs, dqdt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: im, km, ks
    
    real, intent (in), dimension (im, km) :: ta, pa, qv
    
    real, intent (out), dimension (im, km) :: qs
    
    real, intent (out), dimension (im, km), optional :: dqdt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i, k, it
    
    real :: eps10, ap1, tmin
    
    real, dimension (im, km) :: es
    
    if (.not. tables_are_initialized) call qsmith_init
    
    tmin = tice - 160.
    do k = ks, km
        do i = 1, im
            ap1 = 10. * dim (ta (i, k), tmin) + 1.
            ap1 = min (2621., ap1)
            it = ap1
            es (i, k) = table1 (it) + (ap1 - it) * des1 (it)
            qs (i, k) = eps * es (i, k) * (1. + zvir * qv (i, k)) / pa (i, k)
        enddo
    enddo
    
    if (present (dqdt)) then
        eps10 = 10. * eps
        do k = ks, km
            do i = 1, im
                ap1 = 10. * dim (ta (i, k), tmin) + 1.
                ap1 = min (2621., ap1) - 0.5
                it = ap1
                dqdt (i, k) = eps10 * (des1 (it) + (ap1 - it) * (des1 (it + 1) - des1 (it))) * &
                    (1. + zvir * qv (i, k)) / pa (i, k)
            enddo
        enddo
    endif
    
end subroutine qsmith

! =======================================================================
! compute wet buld temperature
! =======================================================================

real function wet_bulb (qv, ta, den)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: qv, ta, den
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    real :: qs, tp, dqdt, lcp
    
    lcp = hlv / cp_air
    
    wet_bulb = ta
    qs = wqs (wet_bulb, den, dqdt)
    tp = 0.5 * (qs - qv) / (1. + lcp * dqdt) * lcp
    wet_bulb = wet_bulb - tp
    
    ! tp is negative if supersaturated
    if (tp .gt. 0.01) then
        qs = wqs (wet_bulb, den, dqdt)
        tp = (qs - qv) / (1. + lcp * dqdt) * lcp
        wet_bulb = wet_bulb - tp
    endif
    
end function wet_bulb

end module gfdl_cld_mp_mod
