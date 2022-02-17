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

module gfdl_mp_mod

    use fv_arrays_mod, only: r_grid

    implicit none

    private
<<<<<<< HEAD
    
    ! -----------------------------------------------------------------------
    ! interface functions
    ! -----------------------------------------------------------------------
    
    interface wqs
        procedure wqs_trho
        procedure wqs_ptqv
    end interface wqs
    
    interface mqs
        procedure mqs_trho
        procedure mqs_ptqv
    end interface mqs
    
    interface iqs
        procedure iqs_trho
        procedure iqs_ptqv
    end interface iqs
    
    interface mhc
        procedure mhc3
        procedure mhc4
        procedure mhc6
    end interface mhc
    
    interface wet_bulb
        procedure wet_bulb_dry
        procedure wet_bulb_moist
    end interface wet_bulb
    
    ! -----------------------------------------------------------------------
    ! public subroutines, functions, and variables
    ! -----------------------------------------------------------------------
    
    public :: gfdl_mp_init
    public :: gfdl_mp_driver
    public :: gfdl_mp_end
    public :: fast_sat_adj, cld_eff_rad, rad_ref
    public :: qs_init, wqs, mqs, mqs3d
    public :: c_liq, c_ice, rhow, wet_bulb
    
    ! -----------------------------------------------------------------------
    ! precision definition
    ! -----------------------------------------------------------------------
=======

    public gfdl_mp_driver, gfdl_mp_init, gfdl_mp_end
    public wqs1, wqs2, iqs1, iqs2, mpdrv, sedi_heat, warm_rain, revap_racc, &
        linear_prof, icloud, subgrid_z_proc, terminal_fall, check_column, implicit_fall, &
        lagrangian_fall_ppm, cs_profile, cs_limiters, fall_speed, setupm, setup_con, &
        qsmith_init, qs_tablew, qs_table2, qs_table3, qs_table, neg_adj, acr3d, smlt, gmlt, &
        wet_bulb, qsmith, qs_blend, es3_table1d, es2_table1d, esw_table1d, es2_table, &
        esw_table, d_sat, qs1d_m, wqsat_moist, wqsat2_moist, qs1d_moist, revap_rac1, &
        wqs2_vect, rhow, rhor, rhos, rhog, rhoh, rnzr, rnzs, rnzg, rnzh, rvgas, rdgas, &
        grav, hlv, hlf, cp_air, cp_vap, cv_air, cv_vap, c_ice, c_liq, dc_vap, dc_ice, &
        t_ice, t_wfr, e00, pi, zvir, rgrav
>>>>>>> user/lnz/shield2022
    
    integer, parameter :: r8 = 8 ! double precision
    
    ! -----------------------------------------------------------------------
    ! initialization conditions
    ! -----------------------------------------------------------------------
    
    logical :: tables_are_initialized = .false. ! initialize satuation tables
    
    ! -----------------------------------------------------------------------
    ! physics constants
    ! -----------------------------------------------------------------------
    
    real, parameter :: grav = 9.80665 ! acceleration due to gravity (m/s^2), ref: IFS
    
    real, parameter :: rgrav = 1.0 / grav ! inversion of gravity acceleration (s^2/m)
    
    real, parameter :: pi = 4.0 * atan (1.0) ! ratio of circle circumference to diameter

<<<<<<< HEAD
    real, parameter :: boltzmann = 1.38064852e-23 ! boltzmann constant (J/K)
    real, parameter :: avogadro = 6.02214076e23 ! avogadro number (1/mol)
    real, parameter :: runiver = avogadro * boltzmann ! 8.314459727525675, universal gas constant (J/K/mol)
    real, parameter :: mmd = 2.89644e-2 ! dry air molar mass (kg/mol), ref: IFS
    real, parameter :: mmv = 1.80153e-2 ! water vapor molar mass (kg/mol), ref: IFS
    
    real, parameter :: rdgas = 287.05 ! gas constant for dry air (J/kg/K): ref: GFDL, GFS
    real, parameter :: rvgas = 461.50 ! gas constant for water vapor (J/kg/K): ref: GFDL, GFS
    !real, parameter :: rdgas = runiver / mmd ! 287.0578961596192, gas constant for dry air (J/kg/K)
    !real, parameter :: rvgas = runiver / mmv ! 461.52213549181386, gas constant for water vapor (J/kg/K)

    real, parameter :: zvir = rvgas / rdgas - 1. ! 0.6077667316114637
    
    real, parameter :: tice = 273.15 ! freezing temperature (K): ref: GFDL, GFS
    !real, parameter :: tice = 273.16 ! freezing temperature (K), ref: IFS
    
    real, parameter :: cp_air = 1004.6 ! heat capacity of dry air at constant pressure (J/kg/K): ref: GFDL, GFS
    real, parameter :: cv_air = cp_air - rdgas ! 717.55, heat capacity of dry air at constant volume (J/kg/K): ref: GFDL, GFS
    !real, parameter :: cp_air = 7. / 2. * rdgas ! 1004.7026365586671, heat capacity of dry air at constant pressure (J/kg/K)
    !real, parameter :: cv_air = 5. / 2. * rdgas ! 717.644740399048, heat capacity of dry air at constant volume (J/kg/K)
    real, parameter :: cp_vap = 4.0 * rvgas ! 1846.0885419672554, heat capacity of water vapor at constnat pressure (J/kg/K)
    real, parameter :: cv_vap = 3.0 * rvgas ! 1384.5664064754415, heat capacity of water vapor at constant volume (J/kg/K)
    
    real, parameter :: c_ice = 2.106e3 ! heat capacity of ice at 0 deg C (J/kg/K), ref: IFS
    real, parameter :: c_liq = 4.218e3 ! heat capacity of water at 0 deg C (J/kg/K), ref: IFS
    
    real, parameter :: dc_vap = cp_vap - c_liq ! - 2371.9114580327446, isobaric heating / cooling (J/kg/K)
    real, parameter :: dc_ice = c_liq - c_ice ! 2112.0, isobaric heating / colling (J/kg/K)
    real, parameter :: d2_ice = cp_vap - c_ice ! - 259.9114580327446, isobaric heating / cooling (J/kg/K)
    
    real, parameter :: hlv = 2.5e6 ! latent heat of evaporation at 0 deg C (J/kg): ref: GFDL, GFS
    real, parameter :: hlf = 3.3358e5 ! latent heat of fusion at 0 deg C (J/kg): ref: GFDL, GFS
    !real, parameter :: hlv = 2.5008e6 ! latent heat of evaporation at 0 deg C (J/kg), ref: IFS
    !real, parameter :: hlf = 3.345e5 ! latent heat of fusion at 0 deg C (J/kg), ref: IFS
    
    real, parameter :: visd = 1.717e-5 ! dynamics viscosity of air at 0 deg C and 1000 hPa (Mason, 1971) (kg/m/s)
    real, parameter :: visk = 1.35e-5 ! kinematic viscosity of air at 0 deg C  and 1000 hPa (Mason, 1971) (m^2/s)
    real, parameter :: vdifu = 2.25e-5 ! diffusivity of water vapor in air at 0 deg C  and 1000 hPa (Mason, 1971) (m^2/s)
    real, parameter :: tcond = 2.40e-2 ! thermal conductivity of air at 0 deg C  and 1000 hPa (Mason, 1971) (J/m/s/K)

    real, parameter :: rho0 = 1.0 ! reference air density (kg/m^3), ref: IFS
    real, parameter :: cdg = 3.15121 ! drag coefficient of graupel (Locatelli and Hobbs, 1974)
    real, parameter :: cdh = 0.5 ! drag coefficient of hail (Heymsfield and Wright, 2014)
    
    real (kind = r8), parameter :: lv0 = hlv - dc_vap * tice ! 3148711.3338762247, evaporation latent heat coeff. at 0 deg K (J/kg)
    real (kind = r8), parameter :: li0 = hlf - dc_ice * tice ! - 242413.92000000004, fussion latent heat coeff. at 0 deg K (J/kg)
    real (kind = r8), parameter :: li2 = lv0 + li0 ! 2906297.413876225, sublimation latent heat coeff. at 0 deg K (J/kg)
    
    real (kind = r8), parameter :: e00 = 611.21 ! saturation vapor pressure at 0 deg C (Pa), ref: IFS
    
    ! -----------------------------------------------------------------------
    ! predefined parameters
    ! -----------------------------------------------------------------------
    
    integer, parameter :: length = 2621 ! length of the saturation table
    
    real, parameter :: qcmin = 1.0e-15 ! min value for cloud condensates (kg/kg)
    real, parameter :: qfmin = 1.0e-8 ! min value for sedimentation (kg/kg)
    
    real, parameter :: dz_min = 1.0e-2 ! used for correcting flipped height (m)
    
    real, parameter :: rhow = 1.0e3 ! density of cloud water (kg/m^3)
    real, parameter :: rhoi = 9.17e2 ! density of cloud ice (kg/m^3)
    real, parameter :: rhor = 1.0e3 ! density of rain (Lin et al. 1983) (kg/m^3)
    real, parameter :: rhos = 1.0e2 ! density of snow (Lin et al. 1983) (kg/m^3)
    real, parameter :: rhog = 4.0e2 ! density of graupel (Rutledge and Hobbs 1984) (kg/m^3)
    real, parameter :: rhoh = 9.17e2 ! density of hail (Lin et al. 1983) (kg/m^3)
    
    real, parameter :: dt_fr = 8.0 ! t_wfr - dt_fr: minimum temperature water can exist (Moore and Molinero 2011)
    
    real (kind = r8), parameter :: one_r8 = 1.0 ! constant 1
    
=======
    logical :: module_is_initialized = .false.
    logical :: qsmith_tables_initialized = .false.

    real, parameter :: grav = 9.80665 ! gfs: acceleration due to gravity
    real, parameter :: rdgas = 287.05 ! gfs: gas constant for dry air
    real, parameter :: rvgas = 461.50 ! gfs: gas constant for water vapor
    real, parameter :: cp_air = 1.0046e3 ! gfs: heat capacity of dry air at constant pressure
    real, parameter :: hlv = 2.5e6 ! gfs: latent heat of evaporation
    real, parameter :: hlf = 3.3358e5 ! gfs: latent heat of fusion
    real, parameter :: pi = 3.1415926535897931 ! gfs: ratio of circle circumference to diameter

    ! real, parameter :: cp_air = rdgas * 7. / 2. ! 1004.675, heat capacity of dry air at constant pressure
    real, parameter :: cp_vap = 4.0 * rvgas ! 1846.0, heat capacity of water vapore at constnat pressure
    ! real, parameter :: cv_air = 717.56 ! satoh value, heat capacity of dry air at constant volume
    real, parameter :: cv_air = cp_air - rdgas ! 717.55, heat capacity of dry air at constant volume
    ! real, parameter :: cv_vap = 1410.0 ! emanuel value, heat capacity of water vapor at constant volume
    real, parameter :: cv_vap = 3.0 * rvgas ! 1384.5, heat capacity of water vapor at constant volume

    ! http: // www.engineeringtoolbox.com / ice - thermal - properties - d_576.html
    ! c_ice = 2050.0 at 0 deg c
    ! c_ice = 2000.0 at - 10 deg c
    ! c_ice = 1943.0 at - 20 deg c
    ! c_ice = 1882.0 at - 30 deg c
    ! c_ice = 1818.0 at - 40 deg c

    ! https: // www.engineeringtoolbox.com / specific - heat - capacity - water - d_660.html
    ! c_liq = 4219.9 at 0.01 deg c
    ! c_liq = 4195.5 at 10 deg c
    ! c_liq = 4184.4 at 20 deg c
    ! c_liq = 4180.1 at 30 deg c
    ! c_liq = 4179.6 at 40 deg c

    ! the following two are from emanuel's book "atmospheric convection"
    ! real, parameter :: c_ice = 2.106e3 ! heat capacity of ice at 0 deg c: c = c_ice + 7.3 * (t - tice)
    ! real, parameter :: c_liq = 4.190e3 ! heat capacity of water at 0 deg c
    ! real, parameter :: c_ice = 1.972e3 ! gfdl: heat capacity of ice at - 15 deg c
    ! real, parameter :: c_liq = 4.1855e3 ! gfdl: heat capacity of water at 15 deg c
    ! real, parameter :: c_ice = 2.106e3 ! gfs: heat capacity of ice at 0 deg c
    ! real, parameter :: c_liq = 4.1855e3 ! gfs: heat capacity of liquid at 15 deg c
    real, parameter :: c_ice = 2.106e3 ! ifs: heat capacity of ice at 0 deg c
    real, parameter :: c_liq = 4.218e3 ! ifs: heat capacity of water at 0 deg c

    real, parameter :: eps = rdgas / rvgas ! 0.6219934995
    real, parameter :: zvir = rvgas / rdgas - 1. ! 0.6077338443

    real, parameter :: dc_vap = cp_vap - c_liq ! - 2.372e3, isobaric heating / cooling
    real, parameter :: dc_ice = c_liq - c_ice ! 2.112e3, isobaric heating / colling

    real, parameter :: t_ice = 273.16 ! freezing temperature
    real, parameter :: table_ice = 273.16 ! freezing point for qs table
    real :: t_wfr ! complete freezing temperature

    real (kind = r_grid), parameter :: e00 = 611.21 ! ifs: saturation vapor pressure at 0 deg c
    ! real (kind = r_grid), parameter :: e00 = 610.71 ! gfdl: saturation vapor pressure at 0 deg c

    real, parameter :: hlv0 = hlv ! gfs: evaporation latent heat coefficient at 0 deg c
    ! real, parameter :: hlv0 = 2.501e6 ! emanuel value
    real, parameter :: hlf0 = hlf ! gfs: fussion latent heat coefficient at 0 deg c
    ! real, parameter :: hlf0 = 3.337e5 ! emanuel value

    real, parameter :: lv0 = hlv0 - dc_vap * t_ice ! 3.14893552e6, evaporation latent heat coefficient at 0 deg k
    real, parameter :: li0 = hlf0 - dc_ice * t_ice ! - 2.2691392e5, fussion latend heat coefficient at 0 deg k

    real (kind = r_grid), parameter :: d2ice = cp_vap - c_ice ! - 260.0, isobaric heating / cooling
    real (kind = r_grid), parameter :: li2 = lv0 + li0 ! 2.9220216e6, sublimation latent heat coefficient at 0 deg k

    real, parameter :: qrmin = 1.e-8 ! min value for cloud condensates
    real, parameter :: qvmin = 1.e-20 ! min value for water vapor (treated as zero)
    real, parameter :: qcmin = 1.e-12 ! min value for cloud condensates

    real, parameter :: vr_min = 1.e-3 ! min fall speed for rain
    real, parameter :: vf_min = 1.e-5 ! min fall speed for cloud ice, snow, graupel

    real, parameter :: dz_min = 1.e-2 ! used for correcting flipped height

    real, parameter :: sfcrho = 1.2 ! surface air density

    real, parameter :: rnzr = 8.0e6 ! lin et al. 1983
    real, parameter :: rnzs = 3.0e6 ! lin et al. 1983
    real, parameter :: rnzg = 4.0e6 ! rutledge and hobbs 1984
    ! lmh, 20170929
    real, parameter :: rnzh = 4.0e4 ! lin et al. 1983

    real, parameter :: rhow = 1.0e3 ! density of cloud water
    real, parameter :: rhor = 1.0e3 ! lin et al. 1983
    real, parameter :: rhos = 0.1e3 ! lin et al. 1983
    real, parameter :: rhog = 0.4e3 ! rutledge and hobbs 1984
    ! lmh, 20170929
    real, parameter :: rhoh = 0.917e3 ! lin et al. 1983

    real, parameter :: rgrav = 1. / grav

    real :: cracs, csacr, cgacr, cgacs, csacw, craci, csaci, cgacw, cgaci, cracw ! constants for accretions
    real :: acco (3, 4) ! constants for accretions
    ! constants for sublimation / deposition, freezing / melting, condensation / evaporation
    real :: cssub (5), cgsub (5), crevp (5), cgfr (2), csmlt (5), cgmlt (5)

    real :: es0, ces0
    real :: pie, fac_rc
    real :: c_air, c_vap

    real :: lat2, lcp, icp, tcp ! used in bigg mechanism and wet bulk

    real :: d0_vap ! the same as dc_vap, except that cp_vap can be cp_vap or cv_vap
    real (kind = r_grid) :: lv00, li00, li20
    real (kind = r_grid) :: d1_vap, d1_ice, c1_vap, c1_liq, c1_ice
    real (kind = r_grid), parameter :: one_r8 = 1.

    real, allocatable :: table (:), table2 (:), table3 (:), tablew (:)
    real, allocatable :: des (:), des2 (:), des3 (:), desw (:)

    logical :: tables_are_initialized = .false.

    real, parameter :: dt_fr = 8. ! homogeneous freezing of all cloud water at t_wfr - dt_fr
    ! minimum temperature water can exist (moore & molinero nov. 2011, nature)
    ! dt_fr can be considered as the error bar

    real, parameter :: p0_min = 100. ! minimum pressure (pascal) for mp to operate
    real :: p_min

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! namelist parameters
    ! -----------------------------------------------------------------------

    integer :: ntimes = 1 ! cloud microphysics sub cycles
<<<<<<< HEAD
    
    integer :: cfflag = 1 ! cloud fraction scheme
    ! 1: GFDL cloud scheme
    ! 2: Xu and Randall (1996)
    ! 3: Park et al. (2016)
    ! 4: Gultepe and Isaac (2007)
    
    integer :: icloud_f = 0 ! GFDL cloud scheme
    ! 0: subgrid variability based scheme
    ! 1: same as 0, but for old fvgfs implementation
    ! 2: binary cloud scheme
    ! 3: extension of 0
    
    integer :: irain_f = 0 ! cloud water to rain auto conversion scheme
    ! 0: subgrid variability based scheme
    ! 1: no subgrid varaibility
    
    integer :: inflag = 1 ! ice nucleation scheme
    ! 1: Hong et al. (2004)
    ! 2: Meyers et al. (1992)
    ! 3: Meyers et al. (1992)
    ! 4: Cooper (1986)
    ! 5: Fletcher (1962)
    
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
    ! 4: effective radius
    ! 5: mass-weighted effective radius
    
    integer :: reiflag = 5 ! cloud ice effective radius scheme
    ! 1: Heymsfield and Mcfarquhar (1996)
    ! 2: Donner et al. (1997)
    ! 3: Fu (2007)
    ! 4: Kristjansson et al. (2000)
    ! 5: Wyser (1998)
    ! 6: Sun and Rikus (1999), Sun (2001)
    ! 7: effective radius
    ! 8: mass-weighted effective radius
    
    integer :: rerflag = 1 ! rain effective radius scheme
    ! 1: effective radius
    ! 2: mass-weighted effective radius
    
    integer :: resflag = 1 ! snow effective radius scheme
    ! 1: effective radius
    ! 2: mass-weighted effective radius
    
    integer :: regflag = 1 ! graupel effective radius scheme
    ! 1: effective radius
    ! 2: mass-weighted effective radius
    
    integer :: radr_flag = 1 ! radar reflectivity for rain
    ! 1: Mark Stoelinga (2005)
    ! 2: Smith et al. (1975), Tong and Xue (2005)
    ! 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
    
    integer :: rads_flag = 1 ! radar reflectivity for snow
    ! 1: Mark Stoelinga (2005)
    ! 2: Smith et al. (1975), Tong and Xue (2005)
    ! 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
    
    integer :: radg_flag = 1 ! radar reflectivity for graupel
    ! 1: Mark Stoelinga (2005)
    ! 2: Smith et al. (1975), Tong and Xue (2005)
    ! 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
    
    logical :: do_sedi_uv = .true. ! transport of horizontal momentum in sedimentation
    logical :: do_sedi_w = .true. ! transport of vertical momentum in sedimentation
=======

    integer :: icloud_f = 0 ! cloud scheme
    integer :: irain_f = 0 ! cloud water to rain auto conversion scheme

    logical :: sedi_transport = .true. ! transport of momentum in sedimentation
    logical :: do_sedi_w = .true. ! transport of vertical momentum during sedimentation
>>>>>>> user/lnz/shield2022
    logical :: do_sedi_heat = .true. ! transport of heat in sedimentation
    logical :: do_sedi_melt = .true. ! melt cloud ice, snow, and graupel during sedimentation
    
    logical :: do_qa = .true. ! do inline cloud fraction
    logical :: rad_snow = .true. ! include snow in cloud fraciton calculation
    logical :: rad_graupel = .true. ! include graupel in cloud fraction calculation
    logical :: rad_rain = .true. ! include rain in cloud fraction calculation
    logical :: do_cld_adj = .false. ! do cloud fraction adjustment
    
    logical :: use_ppm = .false. ! use ppm fall scheme
    
    logical :: use_implicit_fall = .true. ! use implicit fall scheme
    logical :: use_explicit_fall = .false. ! use explicit fall scheme
    
    logical :: z_slope_liq = .true. ! use linear mono slope for autocconversions
    logical :: z_slope_ice = .true. ! use linear mono slope for autocconversions
    
    logical :: use_rhc_cevap = .false. ! cap of rh for cloud water evaporation
    logical :: use_rhc_revap = .false. ! cap of rh for rain evaporation
    
    logical :: const_vw = .false. ! if .ture., the constants are specified by v * _fac
    logical :: const_vi = .false. ! if .ture., the constants are specified by v * _fac
    logical :: const_vs = .false. ! if .ture., the constants are specified by v * _fac
    logical :: const_vg = .false. ! if .ture., the constants are specified by v * _fac
    logical :: const_vr = .false. ! if .ture., the constants are specified by v * _fac
    
    logical :: liq_ice_combine = .false. ! combine all liquid water, combine all solid water
    logical :: snow_grauple_combine = .true. ! combine snow and graupel
    
    logical :: prog_ccn = .false. ! do prognostic ccn (Yi Ming's method)
    
    logical :: fix_negative = .true. ! fix negative water species
    
    logical :: do_cond_timescale = .false. ! whether to apply a timescale to condensation
<<<<<<< HEAD
    
    logical :: do_hail = .false. ! use hail parameters instead of graupel
    
    logical :: consv_checker = .false. ! turn on energy and water conservation checker
    
    logical :: do_warm_rain_mp = .false. ! do warm rain cloud microphysics only
    
    logical :: do_wbf = .false. ! do Wegener Bergeron Findeisen process

    logical :: do_psd_water_fall = .false. ! calculate cloud water terminal velocity based on PSD
    logical :: do_psd_ice_fall = .false. ! calculate cloud ice terminal velocity based on PSD

    logical :: do_psd_water_num = .false. ! calculate cloud water number concentration based on PSD
    logical :: do_psd_ice_num = .false. ! calculate cloud ice number concentration based on PSD

    logical :: do_new_acc_water = .false. ! perform the new accretion for cloud water
    logical :: do_new_acc_ice = .false. ! perform the new accretion for cloud ice
    
    real :: mp_time = 150.0 ! maximum microphysics time step (s)
    
    real :: n0w_sig = 1.1 ! intercept parameter (significand) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994)
    !real :: n0w_sig = 1.4 ! intercept parameter (significand) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994)
    real :: n0i_sig = 1.3 ! intercept parameter (significand) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015)
    !real :: n0i_sig = 9.4 ! intercept parameter (significand) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015)
    real :: n0r_sig = 8.0 ! intercept parameter (significand) of rain (Lin et al. 1983) (1/m^4) (Marshall and Palmer 1948)
    real :: n0s_sig = 3.0 ! intercept parameter (significand) of snow (Lin et al. 1983) (1/m^4) (Gunn and Marshall 1958)
    real :: n0g_sig = 4.0 ! intercept parameter (significand) of graupel (Rutledge and Hobbs 1984) (1/m^4) (Houze et al. 1979)
    real :: n0h_sig = 4.0 ! intercept parameter (significand) of hail (Lin et al. 1983) (1/m^4) (Federer and Waldvogel 1975)
    
    real :: n0w_exp = 41 ! intercept parameter (exponent) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994)
    !real :: n0w_exp = 91 ! intercept parameter (exponent) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994)
    real :: n0i_exp = 18 ! intercept parameter (exponent) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015)
    !real :: n0i_exp = 17 ! intercept parameter (exponent) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015)
    real :: n0r_exp = 6 ! intercept parameter (exponent) of rain (Lin et al. 1983) (1/m^4) (Marshall and Palmer 1948)
    real :: n0s_exp = 6 ! intercept parameter (exponent) of snow (Lin et al. 1983) (1/m^4) (Gunn and Marshall 1958)
    real :: n0g_exp = 6 ! intercept parameter (exponent) of graupel (Rutledge and Hobbs 1984) (1/m^4) (Houze et al. 1979)
    real :: n0h_exp = 4 ! intercept parameter (exponent) of hail (Lin et al. 1983) (1/m^4) (Federer and Waldvogel 1975)
    
    real :: muw = 6.0 ! shape parameter of cloud water in Gamma distribution (Martin et al. 1994)
    !real :: muw = 16.0 ! shape parameter of cloud water in Gamma distribution (Martin et al. 1994)
    real :: mui = 3.35 ! shape parameter of cloud ice in Gamma distribution (McFarquhar et al. 2015)
    !real :: mui = 3.54 ! shape parameter of cloud ice in Gamma distribution (McFarquhar et al. 2015)
    real :: mur = 1.0 ! shape parameter of rain in Gamma distribution (Marshall and Palmer 1948)
    real :: mus = 1.0 ! shape parameter of snow in Gamma distribution (Gunn and Marshall 1958)
    real :: mug = 1.0 ! shape parameter of graupel in Gamma distribution (Houze et al. 1979)
    real :: muh = 1.0 ! shape parameter of hail in Gamma distribution (Federer and Waldvogel 1975)
    
    real :: alinw = 3.e7 ! "a" in Lin et al. (1983) for cloud water (Ikawa and Saito 1990)
    real :: alini = 7.e2 ! "a" in Lin et al. (1983) for cloud ice (Ikawa and Saita 1990)
    real :: alinr = 842.0 ! "a" in Lin et al. (1983) for rain (Liu and Orville 1969)
    real :: alins = 4.8 ! "a" in Lin et al. (1983) for snow (straka 2009)
    real :: aling = 1.0 ! "a" in Lin et al. (1983), similar to a, but for graupel (Pruppacher and Klett 2010)
    real :: alinh = 1.0 ! "a" in Lin et al. (1983), similar to a, but for hail (Pruppacher and Klett 2010)

    real :: blinw = 2.0 ! "b" in Lin et al. (1983) for cloud water (Ikawa and Saito 1990)
    real :: blini = 1.0 ! "b" in Lin et al. (1983) for cloud ice (Ikawa and Saita 1990)
    real :: blinr = 0.8 ! "b" in Lin et al. (1983) for rain (Liu and Orville 1969)
    real :: blins = 0.25 ! "b" in Lin et al. (1983) for snow (straka 2009)
    real :: bling = 0.5 ! "b" in Lin et al. (1983), similar to b, but for graupel (Pruppacher and Klett 2010)
    real :: blinh = 0.5 ! "b" in Lin et al. (1983), similar to b, but for hail (Pruppacher and Klett 2010)
    
    real :: tice_mlt = 273.16 ! can set ice melting temperature to 268 based on observation (Kay et al. 2016) (K)
    
    real :: t_min = 178.0 ! minimum temperature to freeze - dry all water vapor (K)
    real :: t_sub = 184.0 ! minimum temperature for sublimation of cloud ice (K)
    
    real :: rh_inc = 0.25 ! rh increment for complete evaporation of cloud water and cloud ice
    real :: rh_inr = 0.25 ! rh increment for minimum evaporation of rain
    real :: rh_ins = 0.25 ! rh increment for sublimation of snow
    
    real :: tau_r2g = 900.0 ! rain freezing to graupel time scale (s)
    real :: tau_i2s = 1000.0 ! cloud ice to snow autoconversion time scale (s)
    real :: tau_l2r = 900.0 ! cloud water to rain autoconversion time scale (s)
    real :: tau_v2l = 150.0 ! water vapor to cloud water condensation time scale (s)
    real :: tau_l2v = 300.0 ! cloud water to water vapor evaporation time scale (s)
    real :: tau_revp = 0.0 ! rain evaporation time scale (s)
    real :: tau_imlt = 1200.0 ! cloud ice melting time scale (s)
    real :: tau_smlt = 900.0 ! snow melting time scale (s)
    real :: tau_gmlt = 600.0 ! graupel melting time scale (s)
    real :: tau_wbf = 300.0 ! graupel melting time scale (s)
    
    real :: dw_land = 0.20 ! base value for subgrid deviation / variability over land
    real :: dw_ocean = 0.10 ! base value for subgrid deviation / variability over ocean
    
    real :: ccn_o = 90.0 ! ccn over ocean (1/cm^3)
    real :: ccn_l = 270.0 ! ccn over land (1/cm^3)
    
    real :: rthresh = 10.0e-6 ! critical cloud drop radius (micron) for autoconversion
    
    real :: cld_min = 0.05 ! minimum cloud fraction
    
    real :: qi_lim = 1.0 ! cloud ice limiter (0: no, 1: full, >1: extra) to prevent large ice build up
    
    real :: ql_mlt = 2.0e-3 ! maximum cloud water allowed from melted cloud ice (kg/kg)
    real :: qs_mlt = 1.0e-6 ! maximum cloud water allowed from melted snow (kg/kg)
    
    real :: ql_gen = 1.0e-3 ! maximum cloud water generation during remapping step (kg/kg)
    
    real :: ql0_max = 2.0e-3 ! maximum cloud water value (autoconverted to rain) (kg/kg)
    real :: qi0_max = 1.0e-4 ! maximum cloud ice value (autoconverted to snow) (kg/m^3)
    
    real :: qi0_crt = 1.0e-4 ! cloud ice to snow autoconversion threshold (kg/m^3)
    real :: qs0_crt = 1.0e-3 ! snow to graupel autoconversion threshold (0.6e-3 in Purdue Lin scheme) (kg/m^3)
    
    real :: c_paut = 0.55 ! cloud water to rain autoconversion efficiency
    real :: c_psacw = 1.0 ! cloud water to snow accretion efficiency
    real :: c_psaci = 0.02 ! cloud ice to snow accretion efficiency (was 0.1 in ZETAC)
    real :: c_pracw = 0.9 ! cloud water to rain accretion efficiency
    real :: c_praci = 1.0 ! cloud ice to rain accretion efficiency
    real :: c_pgacw = 1.0 ! cloud water to graupel accretion efficiency
    real :: c_pgaci = 0.05 ! cloud ice to graupel accretion efficiency (was 0.1 in ZETAC)
    real :: c_pracs = 1.0 ! snow to rain accretion efficiency
    real :: c_psacr = 1.0 ! rain to snow accretion efficiency
    real :: c_pgacr = 1.0 ! rain to graupel accretion efficiency
    real :: c_pgacs = 0.01 ! snow to graupel accretion efficiency (was 0.1 in ZETAC)

    real :: is_fac = 0.2 ! cloud ice sublimation temperature factor
    real :: ss_fac = 0.2 ! snow sublimation temperature factor
    real :: gs_fac = 0.2 ! graupel sublimation temperature factor

    real :: rh_fac = 10.0 ! cloud water condensation / evaporation relative humidity factor
    
    real :: vw_fac = 1.0
    real :: vi_fac = 1.0 ! IFS: if const_vi: 1 / 3
    real :: vs_fac = 1.0 ! IFS: if const_vs: 1.
    real :: vg_fac = 1.0 ! IFS: if const_vg: 2.
    real :: vr_fac = 1.0 ! IFS: if const_vr: 4.
    
    real :: vw_max = 0.01 ! maximum fall speed for cloud water (m/s)
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
    
    real :: fi2s_fac = 1.0 ! maximum sink of cloud ice to form snow: 0-1
    real :: fi2g_fac = 1.0 ! maximum sink of cloud ice to form graupel: 0-1
    real :: fs2g_fac = 1.0 ! maximum sink of snow to form graupel: 0-1
    
    real :: beta = 1.22 ! defined in Heymsfield and Mcfarquhar (1996)
    
    real :: rewmin = 5.0, rewmax = 15.0 ! minimum and maximum effective radius for cloud water (micron)
    real :: reimin = 10.0, reimax = 150.0 ! minimum and maximum effective radius for cloud ice (micron)
    real :: rermin = 15.0, rermax = 10000.0 ! minimum and maximum effective radius for rain (micron)
    real :: resmin = 150.0, resmax = 10000.0 ! minimum and maximum effective radius for snow (micron)
    real :: regmin = 150.0, regmax = 10000.0 ! minimum and maximum effective radius for graupel
    !real :: rewmax = 15.0, rermin = 15.0 ! Kokhanovsky (2004)
    
    ! -----------------------------------------------------------------------
    ! local shared variables
    ! -----------------------------------------------------------------------
    
    real :: acco (3, 10), acc (20)
    real :: cracs, csacr, cgacr, cgacs, csacw, craci, csaci, cgacw, cgaci, cracw
    real :: cssub (5), cgsub (5), crevp (5), cgfr (2), csmlt (4), cgmlt (4)
    
    real :: t_wfr, fac_rc, c_air, c_vap, d0_vap
    
    real (kind = r8) :: lv00, li00, li20, cpaut
    real (kind = r8) :: d1_vap, d1_ice, c1_vap, c1_liq, c1_ice
    real (kind = r8) :: vconw, vconr, vconi, vcons, vcong, vconh
    real (kind = r8) :: normw, normr, normi, norms, normg, normh
    real (kind = r8) :: expow, expor, expoi, expos, expog, expoh
    real (kind = r8) :: coeaw, coear, coeai, coeas, coeag, coeah
    real (kind = r8) :: coebw, coebr, coebi, coebs, coebg, coebh
    
    real, allocatable :: table0 (:), table1 (:), table2 (:), table3 (:), table4 (:)
    real, allocatable :: des0 (:), des1 (:), des2 (:), des3 (:), des4 (:)
    
=======

    real :: cld_fac = 1.0 ! multiplication factor for cloud fraction
    real :: cld_min = 0.05 ! minimum cloud fraction
    real :: tice = 273.16 ! set tice = 165. to trun off ice - phase phys (kessler emulator)
    real :: tice_mlt = 273.16 ! set ice melting temperature to 268.0 based on observation (kay et al., 2016, jc)

    real :: t_min = 178. ! min temp to freeze - dry all water vapor
    real :: t_sub = 184. ! min temp for sublimation of cloud ice
    real :: mp_time = 150. ! maximum micro - physics time step (sec)

    real :: rh_inc = 0.25 ! rh increment for complete evaporation of cloud water and cloud ice
    real :: rh_inr = 0.25 ! rh increment for minimum evaporation of rain
    real :: rh_ins = 0.25 ! rh increment for sublimation of snow

    real :: tau_r2g = 900. ! rain freezing during fast_sat
    real :: tau_smlt = 900. ! snow melting
    real :: tau_g2r = 600. ! graupel melting to rain
    real :: tau_imlt = 600. ! cloud ice melting
    real :: tau_i2s = 1000. ! cloud ice to snow auto - conversion
    real :: tau_l2r = 900. ! cloud water to rain auto - conversion
    real :: tau_v2l = 150. ! water vapor to cloud water (condensation)
    real :: tau_l2v = 300. ! cloud water to water vapor (evaporation)
    real :: tau_g2v = 900. ! grapuel sublimation
    real :: tau_v2g = 21600. ! grapuel deposition -- make it a slow process
    real :: tau_revp = 0. ! rain evaporation

    real :: dw_land = 0.20 ! base value for subgrid deviation / variability over land
    real :: dw_ocean = 0.10 ! base value for ocean

    real :: ccn_o = 90. ! ccn over ocean (cm^ - 3)
    real :: ccn_l = 270. ! ccn over land (cm^ - 3)

    real :: rthresh = 10.0e-6 ! critical cloud drop radius (micron)

    ! -----------------------------------------------------------------------
    ! wrf / wsm6 scheme: qi_gen = 4.92e-11 * (1.e3 * exp (0.1 * tmp)) ** 1.33
    ! optimized: qi_gen = 4.92e-11 * exp (1.33 * log (1.e3 * exp (0.1 * tmp)))
    ! qi_gen ~ 4.808e-7 at 0 c; 1.818e-6 at - 10 c, 9.82679e-5 at - 40c
    ! the following value is constructed such that qc_crt = 0 at zero c and @ - 10c matches
    ! wrf / wsm6 ice initiation scheme; qi_crt = qi_gen * min (qi_lim, 0.1 * tmp) / den
    ! -----------------------------------------------------------------------

    real :: sat_adj0 = 0.90 ! adjustment factor (0: no, 1: full) during fast_sat_adj

    real :: qc_crt = 5.0e-8 ! mini condensate mixing ratio to allow partial cloudiness

    real :: qi_lim = 1. ! cloud ice limiter to prevent large ice build up

    real :: ql_mlt = 2.0e-3 ! max value of cloud water allowed from melted cloud ice
    real :: qs_mlt = 1.0e-6 ! max cloud water due to snow melt

    real :: ql_gen = 1.0e-3 ! max cloud water generation during remapping step if do_sat_adj = .t.
    real :: qi_gen = 1.82e-6 ! max cloud ice generation during remapping step

    ! cloud condensate upper bounds: "safety valves" for ql & qi

    real :: ql0_max = 2.0e-3 ! max cloud water value (auto converted to rain)
    real :: qi0_max = 1.0e-4 ! max cloud ice value (by other sources)

    real :: qi0_crt = 1.0e-4 ! cloud ice to snow autoconversion threshold (was 1.e-4)
    ! qi0_crt if negative, its magnitude is used as the mixing ration threshold; otherwise, used as density
    real :: qr0_crt = 1.0e-4 ! rain to snow or graupel / hail threshold
    ! lin et al. (1983) used * mixing ratio * = 1.e-4 (hail)
    real :: qs0_crt = 1.0e-3 ! snow to graupel density threshold (0.6e-3 in purdue lin scheme)

    real :: c_paut = 0.55 ! autoconversion cloud water to rain (use 0.5 to reduce autoconversion)
    real :: c_psaci = 0.02 ! accretion: cloud ice to snow (was 0.1 in zetac)
    real :: c_piacr = 5.0 ! accretion: rain to ice: (not used)
    real :: c_cracw = 0.9 ! rain accretion efficiency
    real :: c_pgacs = 2.0e-3 ! snow to graupel "accretion" eff. (was 0.1 in zetac)

    ! decreasing clin to reduce csacw (so as to reduce cloud water --- > snow)

    real :: alin = 842.0 ! "a" in lin et al. (1983)
    real :: clin = 4.8 ! "c" in lin et al. (1983), 4.8 -- > 6. (to ehance ql -- > qs)

    logical :: const_vi = .false. ! if .t. the constants are specified by v * _fac
    logical :: const_vs = .false. ! if .t. the constants are specified by v * _fac
    logical :: const_vg = .false. ! if .t. the constants are specified by v * _fac
    logical :: const_vr = .false. ! if .t. the constants are specified by v * _fac

    real :: vi_fac = 1. ! ifs: if const_vi: 1 / 3
    real :: vs_fac = 1. ! ifs: if const_vs: 1.
    real :: vg_fac = 1. ! ifs: if const_vg: 2.
    real :: vr_fac = 1. ! ifs: if const_vr: 4.

    real :: vi_max = 0.5 ! max fall speed for ice
    real :: vs_max = 5.0 ! max fall speed for snow
    real :: vg_max = 8.0 ! max fall speed for graupel
    real :: vr_max = 12. ! max fall speed for rain

    real :: xr_a = 0.25 ! p value in xu and randall, 1996
    real :: xr_b = 100. ! alpha_0 value in xu and randall, 1996
    real :: xr_c = 0.49 ! gamma value in xu and randall, 1996

    real :: te_err = 1.e-14 ! 64bit: 1.e-14, 32bit: 1.e-7

    logical :: do_sat_adj = .false. ! has fast saturation adjustments
    logical :: z_slope_liq = .true. ! use linear mono slope for autocconversions
    logical :: z_slope_ice = .false. ! use linear mono slope for autocconversions
    logical :: use_ccn = .false. ! must be true when prog_ccn is false
    logical :: use_ppm = .false. ! use ppm fall scheme
    logical :: use_ppm_ice = .false. ! use ppm fall scheme for cloud ice
    logical :: mono_prof = .true. ! perform terminal fall with mono ppm scheme
    logical :: do_hail = .false. ! use hail parameters instead of graupel
    logical :: hd_icefall = .false. ! use heymsfield and donner, 1990's fall speed of cloud ice
    logical :: use_xr_cloud = .false. ! use xu and randall, 1996's cloud diagnosis
    logical :: use_park_cloud = .false. ! park et al. 2016
    logical :: use_gi_cloud = .false. ! gultepe and isaac (2007, grl)
    logical :: use_rhc_cevap = .false. ! cap of rh for cloud water evaporation
    logical :: use_rhc_revap = .false. ! cap of rh for rain evaporation
    logical :: consv_checker = .false. ! turn on energy and water conservation checker
    logical :: do_warm_rain_mp = .false. ! do warm rain cloud microphysics only
    ! turn off to save time, turn on only in c48 64bit

    real :: g2, log_10

    real :: rh_thres = 0.75
    real :: rhc_cevap = 0.85 ! cloud water
    real :: rhc_revap = 0.85 ! cloud water

    real :: f_dq_p = 1.0
    real :: f_dq_m = 1.0
    logical :: do_cld_adj = .false.

    integer :: inflag = 1 ! ice nucleation scheme
    ! 1: hong et al., 2004
    ! 2: meyers et al., 1992
    ! 3: meyers et al., 1992
    ! 4: cooper, 1986
    ! 5: flecther, 1962

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! namelist
    ! -----------------------------------------------------------------------

    namelist / gfdl_mp_nml / &
<<<<<<< HEAD
        t_min, t_sub, tau_r2g, tau_smlt, tau_gmlt, dw_land, dw_ocean, vw_fac, vi_fac, &
        vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vw_max, vi_max, vs_max, &
        vg_max, vr_max, qs_mlt, qs0_crt, ql0_max, qi0_max, qi0_crt, ifflag, &
        rh_inc, rh_ins, rh_inr, const_vw, const_vi, const_vs, const_vg, const_vr, rthresh, &
        ccn_l, ccn_o, igflag, c_paut, tau_imlt, tau_v2l, tau_l2v, tau_i2s, &
        tau_l2r, qi_lim, ql_gen, do_hail, inflag, c_psacw, c_psaci, c_pracs, &
        c_psacr, c_pgacr, c_pgacs, c_pgacw, c_pgaci, z_slope_liq, z_slope_ice, &
        prog_ccn, c_pracw, c_praci, rad_snow, rad_graupel, rad_rain, cld_min, &
        use_ppm, do_sedi_uv, do_sedi_w, do_sedi_heat, icloud_f, &
        irain_f, xr_a, xr_b, xr_c, ntimes, tau_revp, tice_mlt, do_cond_timescale, &
        mp_time, consv_checker, te_err, use_rhc_cevap, use_rhc_revap, tau_wbf, &
        do_warm_rain_mp, rh_thres, f_dq_p, f_dq_m, do_cld_adj, rhc_cevap, &
        rhc_revap, beta, liq_ice_combine, rewflag, reiflag, rerflag, resflag, &
        regflag, rewmin, rewmax, reimin, reimax, rermin, rermax, resmin, &
        resmax, regmin, regmax, fs2g_fac, fi2s_fac, fi2g_fac, do_sedi_melt, &
        radr_flag, rads_flag, radg_flag, do_wbf, do_psd_water_fall, do_psd_ice_fall, &
        n0w_sig, n0i_sig, n0r_sig, n0s_sig, n0g_sig, n0h_sig, n0w_exp, n0i_exp, &
        n0r_exp, n0s_exp, n0g_exp, n0h_exp, muw, mui, mur, mus, mug, muh, &
        alinw, alini, alinr, alins, aling, alinh, blinw, blini, blinr, blins, bling, blinh, &
        do_new_acc_water, do_new_acc_ice, is_fac, ss_fac, gs_fac, rh_fac, &
        snow_grauple_combine, do_psd_water_num, do_psd_ice_num, use_implicit_fall, &
        use_explicit_fall
    
=======
        t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, dw_land, dw_ocean, &
        vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vi_max, &
        vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max, &
        qi0_crt, do_sat_adj, rh_inc, rh_ins, rh_inr, const_vi, &
        const_vs, const_vg, const_vr, use_ccn, rthresh, ccn_l, ccn_o, qc_crt, &
        tau_g2v, tau_v2g, sat_adj0, tau_imlt, tau_v2l, tau_l2v, &
        tau_i2s, tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, c_pgacs, &
        z_slope_liq, z_slope_ice, prog_ccn, c_cracw, alin, clin, tice, &
        rad_snow, rad_graupel, rad_rain, cld_fac, cld_min, use_ppm, use_ppm_ice, mono_prof, &
        do_sedi_heat, sedi_transport, do_sedi_w, icloud_f, irain_f, &
        ntimes, disp_heat, do_hail, use_xr_cloud, xr_a, xr_b, xr_c, tau_revp, tice_mlt, hd_icefall, &
        do_cond_timescale, mp_time, consv_checker, te_err, use_park_cloud, &
        use_gi_cloud, use_rhc_cevap, use_rhc_revap, inflag, do_warm_rain_mp, &
        rh_thres, f_dq_p, f_dq_m, do_cld_adj

    public &
        t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, dw_land, dw_ocean, &
        vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vi_max, &
        vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max, &
        qi0_crt, do_sat_adj, rh_inc, rh_ins, rh_inr, const_vi, &
        const_vs, const_vg, const_vr, use_ccn, rthresh, ccn_l, ccn_o, qc_crt, &
        tau_g2v, tau_v2g, sat_adj0, tau_imlt, tau_v2l, tau_l2v, &
        tau_i2s, tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, c_pgacs, &
        z_slope_liq, z_slope_ice, prog_ccn, c_cracw, alin, clin, tice, &
        rad_snow, rad_graupel, rad_rain, cld_fac, cld_min, use_ppm, use_ppm_ice, mono_prof, &
        do_sedi_heat, sedi_transport, do_sedi_w, icloud_f, irain_f, &
        ntimes, disp_heat, do_hail, use_xr_cloud, xr_a, xr_b, xr_c, tau_revp, tice_mlt, hd_icefall, &
        do_cond_timescale, mp_time, consv_checker, te_err, use_park_cloud, &
        use_gi_cloud, use_rhc_cevap, use_rhc_revap, inflag, do_warm_rain_mp, &
        rh_thres, f_dq_p, f_dq_m, do_cld_adj

>>>>>>> user/lnz/shield2022
contains

! =======================================================================
! GFDL cloud microphysics initialization
! =======================================================================

subroutine gfdl_mp_init (me, master, nlunit, input_nml_file, logunit, fn_nml)
    
    implicit none
<<<<<<< HEAD
    
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
=======

    logical, intent (in) :: hydrostatic
    logical, intent (in) :: last_step
    logical, intent (in) :: consv_te
    logical, intent (in) :: do_inline_mp

    integer, intent (in) :: is, ie ! physics window
    integer, intent (in) :: ks, ke ! vertical dimension

    real, intent (in) :: dts ! physics time step

    real, intent (in), dimension (is:ie) :: hs, gsize

    real, intent (in), dimension (is:ie, ks:ke) :: dz
    real, intent (in), dimension (is:ie, ks:ke) :: qnl, qni

    real, intent (inout), dimension (is:ie, ks:ke) :: delp
    real, intent (inout), dimension (is:ie, ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    real, intent (inout), dimension (is:ie, ks:ke) :: pt, ua, va, w
    real, intent (inout), dimension (is:ie, ks:ke) :: prefluxr, prefluxi, prefluxs, prefluxg
    real, intent (inout), dimension (is:, ks:) :: q_con, cappa
    real, intent (inout), dimension (is:ie) :: rain, snow, ice, graupel
    real, intent (inout), dimension (is:ie) :: condensation, deposition
    real, intent (inout), dimension (is:ie) :: evaporation, sublimation

    real, intent (inout), dimension (is:ie, ks:ke) :: te
    ! logical :: used
    real, dimension (is:ie) :: w_var
    real, dimension (is:ie, ks:ke) :: vt_r, vt_s, vt_g, vt_i
    real, dimension (is:ie, ks:ke) :: m2_rain, m2_sol

    if (last_step) then
        p_min = p0_min ! final clean - up
>>>>>>> user/lnz/shield2022
    else
        open (unit = nlunit, file = fn_nml, readonly, status = 'old')
    endif
<<<<<<< HEAD
    rewind (nlunit)
    read (nlunit, nml = gfdl_mp_nml)
    close (nlunit)
#endif
    
=======

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! write namelist to log file
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    if (me .eq. master) then
        write (logunit, *) " ================================================================== "
        write (logunit, *) "gfdl_mp_mod"
        write (logunit, nml = gfdl_mp_nml)
    endif
    
    ! -----------------------------------------------------------------------
    ! initialize microphysics variables
    ! -----------------------------------------------------------------------
    
    if (.not. tables_are_initialized) call qs_init
    
    call setup_mp
    
end subroutine gfdl_mp_init

! =======================================================================
! GFDL cloud microphysics driver
! =======================================================================

subroutine gfdl_mp_driver (qv, ql, qr, qi, qs, qg, qa, qnl, qni, pt, wa, &
        ua, va, delz, delp, gsize, dtm, hs, water, rain, ice, snow, graupel, &
        hydrostatic, is, ie, ks, ke, q_con, cappa, consv_te, te, &
        prefluxw, prefluxr, prefluxi, prefluxs, prefluxg, condensation, &
        deposition, evaporation, sublimation, last_step, do_inline_mp)
    
    implicit none
    
=======

    if (hydrostatic) then
        c_air = cp_air
        c_vap = cp_vap
        do_sedi_w = .false.
    else
        c_air = cv_air
        c_vap = cv_vap
    endif
    d0_vap = c_vap - c_liq

    ! scaled constants (to reduce fp errors for 32 - bit) :
    d1_vap = d0_vap / c_air
    d1_ice = dc_ice / c_air

    ! lv0 = hlv0 - (c_vap - c_liq) * t_ice! 3.13905782e6, evaporation latent heat coefficient at 0 deg k
    lv00 = (hlv0 - d0_vap * t_ice) / c_air
    li00 = (hlf0 - dc_ice * t_ice) / c_air
    li20 = lv00 + li00

    c1_vap = c_vap / c_air
    c1_liq = c_liq / c_air
    c1_ice = c_ice / c_air

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    integer, intent (in) :: is, ie, ks, ke
    
    logical, intent (in) :: hydrostatic, last_step, consv_te, do_inline_mp
    
    real, intent (in) :: dtm
    
    real, intent (in), dimension (is:ie) :: hs, gsize
    
    real, intent (in), dimension (is:ie, ks:ke) :: delz, qnl, qni
    
    real, intent (inout), dimension (is:ie, ks:ke) :: delp, pt, ua, va, wa, te
    real, intent (inout), dimension (is:ie, ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    real, intent (inout), dimension (is:ie, ks:ke) :: prefluxw, prefluxr, prefluxi, prefluxs, prefluxg
    
    real, intent (inout), dimension (is:, ks:) :: q_con, cappa
    
    real, intent (inout), dimension (is:ie) :: water, rain, ice, snow, graupel
    real, intent (inout), dimension (is:ie) :: condensation, deposition
    real, intent (inout), dimension (is:ie) :: evaporation, sublimation
    
    ! -----------------------------------------------------------------------
    ! define various heat capacities and latent heat coefficients at 0 deg K
    ! -----------------------------------------------------------------------
    
    call setup_mhc_lhc (hydrostatic)
    
=======

    lat2 = (hlv + hlf) ** 2

    lcp = hlv / cp_air
    icp = hlf / cp_air
    tcp = (hlv + hlf) / cp_air

    ! tendency zero out for am moist processes should be done outside the driver

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! major cloud microphysics driver
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    call mpdrv (hydrostatic, ua, va, wa, delp, pt, qv, ql, qr, qi, qs, qg, qa, &
        qnl, qni, delz, is, ie, ks, ke, dtm, water, rain, ice, snow, graupel, &
        gsize, hs, q_con, cappa, consv_te, te, prefluxw, prefluxr, prefluxi, &
        prefluxs, prefluxg, condensation, deposition, evaporation, sublimation, &
        last_step, do_inline_mp, .false., .true.)
=======

    call mpdrv (hydrostatic, ua, va, w, delp, pt, qv, ql, qr, qi, qs, qg, &
        qa, qnl, qni, dz, is, ie, ks, ke, dts, &
        rain, snow, graupel, ice, m2_rain, m2_sol, gsize, hs, &
        w_var, vt_r, vt_s, vt_g, vt_i, q_con, cappa, consv_te, te, &
        prefluxr, prefluxi, prefluxs, prefluxg, condensation, deposition, &
        evaporation, sublimation, last_step, do_inline_mp)
>>>>>>> user/lnz/shield2022
    
end subroutine gfdl_mp_driver

! =======================================================================
! GFDL cloud microphysics end
! =======================================================================

subroutine gfdl_mp_end
    
    implicit none
<<<<<<< HEAD
    
    ! -----------------------------------------------------------------------
    ! free up memory
    ! -----------------------------------------------------------------------
    
    deallocate (table0)
    deallocate (table1)
    deallocate (table2)
    deallocate (table3)
    deallocate (table4)
    deallocate (des0)
    deallocate (des1)
    deallocate (des2)
    deallocate (des3)
    deallocate (des4)
    
    tables_are_initialized = .false.
    
end subroutine gfdl_mp_end

! =======================================================================
! setup cloud microphysics parameters
! =======================================================================

subroutine setup_mp
    
    implicit none
    
    integer :: i, k
    
    real :: gcon, hcon, scm3, pisq, act (20), ace (20), occ (3), aone
    
    ! -----------------------------------------------------------------------
    ! complete freezing temperature
    ! -----------------------------------------------------------------------
    
    if (do_warm_rain_mp) then
        t_wfr = t_min
    else
        t_wfr = tice - 40.0
    endif
    
    ! -----------------------------------------------------------------------
    ! cloud water autoconversion, Hong et al. (2004)
    ! -----------------------------------------------------------------------
    
    fac_rc = (4. / 3.) * pi * rhow * rthresh ** 3
    
    aone = 2. / 9. * (3. / 4.) ** (4. / 3.) / pi ** (1. / 3.)
    cpaut = c_paut * aone * grav / visd
    
=======

    logical, intent (in) :: hydrostatic
    logical, intent (in) :: last_step
    logical, intent (in) :: consv_te
    logical, intent (in) :: do_inline_mp
    integer, intent (in) :: is, ie, ks, ke
    real, intent (in) :: dt_in
    real, intent (in), dimension (is:ie) :: gsize
    real, intent (in), dimension (is:ie) :: hs
    real, intent (in), dimension (is:ie, ks:ke) :: dz
    real, intent (in), dimension (is:ie, ks:ke) :: qnl, qni

    real, intent (inout), dimension (is:ie, ks:ke) :: delp
    real, intent (inout), dimension (is:ie, ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    real, intent (inout), dimension (is:ie, ks:ke) :: pt, ua, va, w
    real, intent (inout), dimension (is:ie, ks:ke) :: prefluxr, prefluxi, prefluxs, prefluxg
    real, intent (inout), dimension (is:, ks:) :: q_con, cappa
    real, intent (inout), dimension (is:ie) :: rain, snow, ice, graupel
    real, intent (inout), dimension (is:ie) :: condensation, deposition
    real, intent (inout), dimension (is:ie) :: evaporation, sublimation

    real, intent (out), dimension (is:ie) :: w_var
    real, intent (out), dimension (is:ie, ks:ke) :: vt_r, vt_s, vt_g, vt_i
    real, intent (out), dimension (is:ie, ks:ke) :: m2_rain, m2_sol
    real, intent (out), dimension (is:ie, ks:ke) :: te
    ! local:
    real, dimension (ks:ke) :: q_liq, q_sol
    real, dimension (ks:ke) :: qvz, qlz, qrz, qiz, qsz, qgz, qaz
    real, dimension (ks:ke) :: vtiz, vtsz, vtgz, vtrz, pfr, pfi, pfs, pfg
    real, dimension (ks:ke) :: dp1, dz1
    real, dimension (ks:ke) :: den, p1, denfac
    real, dimension (ks:ke) :: ccn, cin, c_praut, m1_rain, m1_sol, m1
    real, dimension (ks:ke) :: u0, v0, u1, v1, w1

    real (kind = r_grid), dimension (is:ie, ks:ke) :: te_beg, te_end, tw_beg, tw_end
    real (kind = r_grid), dimension (is:ie, ks:ke) :: te_beg_0, te_end_0, tw_beg_0, tw_end_0
    real (kind = r_grid), dimension (is:ie) :: te_b_beg, te_b_end, tw_b_beg, tw_b_end, dte, te_loss
    real (kind = r_grid), dimension (is:ie) :: te_b_beg_0, te_b_end_0, tw_b_beg_0, tw_b_end_0
    real (kind = r_grid), dimension (ks:ke) :: te1, te2

    real :: cpaut, rh_adj, rh_rain
    real :: r1, s1, i1, g1, rdt, ccn0
    real :: dt_rain
    real :: s_leng, t_land, t_ocean, h_var, tmp
    real (kind = r_grid), dimension (ks:ke) :: dp0, tz, cvm
    real (kind = r_grid) :: con_r8, c8
    real :: convt
    real :: dts, q_cond
    real :: nl, ni
    real :: cond, dep, reevap, sub

    integer :: i, k, n

    ntimes = max (ntimes, int (dt_in / min (dt_in, mp_time)))
    dts = dt_in / real (ntimes)

    dt_rain = dts * 0.5
    rdt = one_r8 / dts

    dte = 0.0

    ! convert to mm / day
    convt = 86400. * rdt * rgrav

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! terminal velocities parameters, Lin et al. (1983)
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    gcon = (4. * grav * rhog / (3. * cdg * rho0)) ** 0.5
    hcon = (4. * grav * rhoh / (3. * cdh * rho0)) ** 0.5
    
    vconw = alinw * gamma (3 + muw + blinw) / gamma (3 + muw)
    vconi = alini * gamma (3 + mui + blini) / gamma (3 + mui)
    vconr = alinr * gamma (3 + mur + blinr) / gamma (3 + mur)
    vcons = alins * gamma (3 + mus + blins) / gamma (3 + mus)
    vcong = aling * gamma (3 + mug + bling) / gamma (3 + mug) * gcon
    vconh = alinh * gamma (3 + muh + blinh) / gamma (3 + muh) * hcon
    
    ! -----------------------------------------------------------------------
    ! slope parameters, Lin et al. (1983)
    ! -----------------------------------------------------------------------
    
    normw = pi * rhow * n0w_sig * gamma (muw + 3)
    normi = pi * rhoi * n0i_sig * gamma (mui + 3)
    normr = pi * rhor * n0r_sig * gamma (mur + 3)
    norms = pi * rhos * n0s_sig * gamma (mus + 3)
    normg = pi * rhog * n0g_sig * gamma (mug + 3)
    normh = pi * rhoh * n0h_sig * gamma (muh + 3)

    expow = exp (n0w_exp / (muw + 3) * log (10.))
    expoi = exp (n0i_exp / (mui + 3) * log (10.))
    expor = exp (n0r_exp / (mur + 3) * log (10.))
    expos = exp (n0s_exp / (mus + 3) * log (10.))
    expog = exp (n0g_exp / (mug + 3) * log (10.))
    expoh = exp (n0h_exp / (muh + 3) * log (10.))

    coeaw = exp (3 / (muw + 3) * log (n0w_sig)) * gamma (muw) * exp (3 * n0w_exp / (muw + 3) * log (10.))
    coeai = exp (3 / (mui + 3) * log (n0i_sig)) * gamma (muw) * exp (3 * n0i_exp / (mui + 3) * log (10.))
    coear = exp (3 / (mur + 3) * log (n0r_sig)) * gamma (muw) * exp (3 * n0r_exp / (mur + 3) * log (10.))
    coeas = exp (3 / (mus + 3) * log (n0s_sig)) * gamma (muw) * exp (3 * n0s_exp / (mus + 3) * log (10.))
    coeag = exp (3 / (mug + 3) * log (n0g_sig)) * gamma (muw) * exp (3 * n0g_exp / (mug + 3) * log (10.))
    coeah = exp (3 / (muh + 3) * log (n0h_sig)) * gamma (muw) * exp (3 * n0h_exp / (muh + 3) * log (10.))

    coebw = exp (muw / (muw + 3) * log (pi / 6 * rhow * gamma (muw + 3)))
    coebi = exp (mui / (mui + 3) * log (pi / 6 * rhoi * gamma (mui + 3)))
    coebr = exp (mur / (mur + 3) * log (pi / 6 * rhor * gamma (mur + 3)))
    coebs = exp (mus / (mus + 3) * log (pi / 6 * rhos * gamma (mus + 3)))
    coebg = exp (mug / (mug + 3) * log (pi / 6 * rhog * gamma (mug + 3)))
    coebh = exp (muh / (muh + 3) * log (pi / 6 * rhoh * gamma (muh + 3)))
    
    ! -----------------------------------------------------------------------
    ! Schmidt number, Sc ** (1 / 3) in Lin et al. (1983)
    ! -----------------------------------------------------------------------
    
    scm3 = exp (1. / 3. * log (visk / vdifu))
    
    pisq = pi * pi
    
    ! -----------------------------------------------------------------------
    ! accretion between cloud water, cloud ice, rain, snow, and graupel or hail, Lin et al. (1983)
    ! -----------------------------------------------------------------------

    cracw = pi * n0r_sig * alinr * gamma (2 + mur + blinr) / &
         (4. * exp (1. / (mur + 3) * (2 + mur + blinr) * log (normr))) * &
         exp ((1 - blinr) * log (expor))
    craci = pi * n0r_sig * alinr * gamma (2 + mur + blinr) / &
         (4. * exp (1. / (mur + 3) * (2 + mur + blinr) * log (normr))) * &
         exp ((1 - blinr) * log (expor))
    csacw = pi * n0s_sig * alins * gamma (2 + mus + blins) / &
         (4. * exp (1. / (mus + 3) * (2 + mus + blins) * log (norms))) * &
         exp ((1 - blins) * log (expos))
    csaci = pi * n0s_sig * alins * gamma (2 + mus + blins) / &
         (4. * exp (1. / (mus + 3) * (2 + mus + blins) * log (norms))) * &
         exp ((1 - blins) * log (expos))
    if (do_hail) then
        cgacw = pi * n0h_sig * alinh * gamma (2 + muh + blinh) * hcon / &
             (4. * exp (1. / (muh + 3) * (2 + muh + blinh) * log (normh))) * &
             exp ((1 - blinh) * log (expoh))
        cgaci = pi * n0h_sig * alinh * gamma (2 + muh + blinh) * hcon / &
             (4. * exp (1. / (muh + 3) * (2 + muh + blinh) * log (normh))) * &
             exp ((1 - blinh) * log (expoh))
    else
        cgacw = pi * n0g_sig * aling * gamma (2 + mug + bling) * gcon / &
             (4. * exp (1. / (mug + 3) * (2 + mug + bling) * log (normg))) * &
             exp ((1 - bling) * log (expog))
        cgaci = pi * n0g_sig * aling * gamma (2 + mug + bling) * gcon / &
             (4. * exp (1. / (mug + 3) * (2 + mug + bling) * log (normg))) * &
             exp ((1 - bling) * log (expog))
    endif

    if (do_new_acc_water) then
    
        cracw = pisq * n0r_sig * n0w_sig * rhow / 24.
        csacw = pisq * n0s_sig * n0w_sig * rhow / 24.
        if (do_hail) then
            cgacw = pisq * n0h_sig * n0w_sig * rhow / 24.
        else
            cgacw = pisq * n0g_sig * n0w_sig * rhow / 24.
        endif

    endif

    if (do_new_acc_ice) then
    
        craci = pisq * n0r_sig * n0i_sig * rhoi / 24.
        csaci = pisq * n0s_sig * n0i_sig * rhoi / 24.
        if (do_hail) then
            cgaci = pisq * n0h_sig * n0i_sig * rhoi / 24.
        else
            cgaci = pisq * n0g_sig * n0i_sig * rhoi / 24.
        endif

    endif

    cracw = cracw * c_pracw
    craci = craci * c_praci
    csacw = csacw * c_psacw
    csaci = csaci * c_psaci
    cgacw = cgacw * c_pgacw
    cgaci = cgaci * c_pgaci

    ! -----------------------------------------------------------------------
    ! accretion between rain, snow, and graupel or hail, Lin et al. (1983)
    ! -----------------------------------------------------------------------
    
    cracs = pisq * n0r_sig * n0s_sig * rhos / 24.
    csacr = pisq * n0s_sig * n0r_sig * rhor / 24.
    if (do_hail) then
        cgacr = pisq * n0h_sig * n0r_sig * rhor / 24.
        cgacs = pisq * n0h_sig * n0s_sig * rhos / 24.
    else
        cgacr = pisq * n0g_sig * n0r_sig * rhor / 24.
        cgacs = pisq * n0g_sig * n0s_sig * rhos / 24.
    endif
    
    cracs = cracs * c_pracs
    csacr = csacr * c_psacr
    cgacr = cgacr * c_pgacr
    cgacs = cgacs * c_pgacs
    
    ! act / ace / acc:
    !  1 -  2: racs (s - r)
    !  3 -  4: sacr (r - s)
    !  5 -  6: gacr (r - g)
    !  7 -  8: gacs (s - g)
    !  9 - 10: racw (w - r)
    ! 11 - 12: raci (i - r)
    ! 13 - 14: sacw (w - s)
    ! 15 - 16: saci (i - s)
    ! 17 - 18: sacw (w - g)
    ! 19 - 20: saci (i - g)
    
    act (1) = norms
    act (2) = normr
    act (3) = act (2)
    act (4) = act (1)
    act (5) = act (2)
    if (do_hail) then
        act (6) = normh
    else
        act (6) = normg
    endif
    act (7) = act (1)
    act (8) = act (6)
    act (9) = normw
    act (10) = act (2)
    act (11) = normi
    act (12) = act (2)
    act (13) = act (9)
    act (14) = act (1)
    act (15) = act (11)
    act (16) = act (1)
    act (17) = act (9)
    act (18) = act (6)
    act (19) = act (11)
    act (20) = act (6)
    
    ace (1) = expos
    ace (2) = expor
    ace (3) = ace (2)
    ace (4) = ace (1)
    ace (5) = ace (2)
    if (do_hail) then
        ace (6) = expoh
    else
        ace (6) = expog
    endif
    ace (7) = ace (1)
    ace (8) = ace (6)
    ace (9) = expow
    ace (10) = ace (2)
    ace (11) = expoi
    ace (12) = ace (2)
    ace (13) = ace (9)
    ace (14) = ace (1)
    ace (15) = ace (11)
    ace (16) = ace (1)
    ace (17) = ace (9)
    ace (18) = ace (6)
    ace (19) = ace (11)
    ace (20) = ace (6)
    
    acc (1) = mus
    acc (2) = mur
    acc (3) = acc (2)
    acc (4) = acc (1)
    acc (5) = acc (2)
    if (do_hail) then
        acc (6) = muh
    else
        acc (6) = mug
    endif
    acc (7) = acc (1)
    acc (8) = acc (6)
    acc (9) = muw
    acc (10) = acc (2)
    acc (11) = mui
    acc (12) = acc (2)
    acc (13) = acc (9)
    acc (14) = acc (1)
    acc (15) = acc (11)
    acc (16) = acc (1)
    acc (17) = acc (9)
    acc (18) = acc (6)
    acc (19) = acc (11)
    acc (20) = acc (6)
    
    occ (1) = 1.
    occ (2) = 2.
    occ (3) = 1.
    
    do i = 1, 3
        do k = 1, 10
            acco (i, k) = occ (i) * gamma (6 + acc (2 * k - 1) - i) * gamma (acc (2 * k) + i - 1) / &
                 (exp (1. / (acc (2 * k - 1) + 3) * (6 + acc (2 * k - 1) - i) * log (act (2 * k - 1))) * &
                exp (1. / (acc (2 * k) + 3) * (acc (2 * k) + i - 1) * log (act (2 * k)))) * &
                exp ((i - 3) * log (ace (2 * k - 1))) * exp ((4 - i) * log (ace (2 * k)))
        enddo
    enddo
    
    ! -----------------------------------------------------------------------
    ! rain evaporation, snow sublimation, and graupel or hail sublimation, Lin et al. (1983)
    ! -----------------------------------------------------------------------
    
    crevp (1) = 2. * pi * vdifu * tcond * rvgas * n0r_sig * gamma (1 + mur) / &
        exp (1. / (mur + 3) * (1 + mur) * log (normr)) * exp (2.0 * log (expor))
    crevp (2) = 0.78
    crevp (3) = 0.31 * scm3 * sqrt (alinr / visk) * gamma ((3 + 2 * mur + blinr) / 2) / &
        exp (1. / (mur + 3) * (3 + 2 * mur + blinr) / 2 * log (normr)) * &
        exp (1. / (mur + 3) * (1 + mur) * log (normr)) / gamma (1 + mur) * &
        exp ((- 1 - blinr) / 2. * log (expor))
    crevp (4) = tcond * rvgas
    crevp (5) = vdifu
    
    cssub (1) = 2. * pi * vdifu * tcond * rvgas * n0s_sig * gamma (1 + mus) / &
        exp (1. / (mus + 3) * (1 + mus) * log (norms)) * exp (2.0 * log (expos))
    cssub (2) = 0.78
    cssub (3) = 0.31 * scm3 * sqrt (alins / visk) * gamma ((3 + 2 * mus + blins) / 2) / &
        exp (1. / (mus + 3) * (3 + 2 * mus + blins) / 2 * log (norms)) * &
        exp (1. / (mus + 3) * (1 + mus) * log (norms)) / gamma (1 + mus) * &
        exp ((- 1 - blins) / 2. * log (expos))
    cssub (4) = tcond * rvgas
    cssub (5) = vdifu
    
    if (do_hail) then
        cgsub (1) = 2. * pi * vdifu * tcond * rvgas * n0h_sig * gamma (1 + muh) / &
            exp (1. / (muh + 3) * (1 + muh) * log (normh)) * exp (2.0 * log (expoh))
        cgsub (2) = 0.78
        cgsub (3) = 0.31 * scm3 * sqrt (alinh * hcon / visk) * gamma ((3 + 2 * muh + blinh) / 2) / &
            exp (1. / (muh + 3) * (3 + 2 * muh + blinh) / 2 * log (normh)) * &
            exp (1. / (muh + 3) * (1 + muh) * log (normh)) / gamma (1 + muh) * &
            exp ((- 1 - blinh) / 2. * log (expoh))
    else
        cgsub (1) = 2. * pi * vdifu * tcond * rvgas * n0g_sig * gamma (1 + mug) / &
            exp (1. / (mug + 3) * (1 + mug) * log (normg)) * exp (2.0 * log (expog))
        cgsub (2) = 0.78
        cgsub (3) = 0.31 * scm3 * sqrt (aling * gcon / visk) * gamma ((3 + 2 * mug + bling) / 2) / &
            exp (1. / (mug + 3) * (3 + 2 * mug + bling) / 2 * log (normg)) * &
            exp (1. / (mug + 3) * (1 + mug) * log (normg)) / gamma (1 + mug) * &
            exp ((- 1 - bling) / 2. * log (expog))
    endif
    cgsub (4) = tcond * rvgas
    cgsub (5) = vdifu
    
    ! -----------------------------------------------------------------------
    ! snow melting, Lin et al. (1983)
    ! -----------------------------------------------------------------------
    
    csmlt (1) = 2. * pi * tcond * n0s_sig * gamma (1 + mus) / &
        exp (1. / (mus + 3) * (1 + mus) * log (norms)) * exp (2.0 * log (expos))
    csmlt (2) = 2. * pi * vdifu * n0s_sig * gamma (1 + mus) / &
        exp (1. / (mus + 3) * (1 + mus) * log (norms)) * exp (2.0 * log (expos))
    csmlt (3) = cssub (2)
    csmlt (4) = cssub (3)
    
    ! -----------------------------------------------------------------------
    ! graupel or hail melting, Lin et al. (1983)
    ! -----------------------------------------------------------------------
    
    if (do_hail) then
        cgmlt (1) = 2. * pi * tcond * n0h_sig * gamma (1 + muh) / &
            exp (1. / (muh + 3) * (1 + muh) * log (normh)) * exp (2.0 * log (expoh))
        cgmlt (2) = 2. * pi * vdifu * n0h_sig * gamma (1 + muh) / &
            exp (1. / (muh + 3) * (1 + muh) * log (normh)) * exp (2.0 * log (expoh))
    else
        cgmlt (1) = 2. * pi * tcond * n0g_sig * gamma (1 + mug) / &
            exp (1. / (mug + 3) * (1 + mug) * log (normg)) * exp (2.0 * log (expog))
        cgmlt (2) = 2. * pi * vdifu * n0g_sig * gamma (1 + mug) / &
            exp (1. / (mug + 3) * (1 + mug) * log (normg)) * exp (2.0 * log (expog))
    endif
    cgmlt (3) = cgsub (2)
    cgmlt (4) = cgsub (3)
    
    ! -----------------------------------------------------------------------
    ! rain freezing, Lin et al. (1983)
    ! -----------------------------------------------------------------------
    
    cgfr (1) = 1.e2 / 36 * pisq * n0r_sig * rhor * gamma (6 + mur) / &
        exp (1. / (mur + 3) * (6 + mur) * log (normr)) * exp (- 3.0 * log (expor))
    cgfr (2) = 0.66
    
end subroutine setup_mp

! =======================================================================
! define various heat capacities and latent heat coefficients at 0 deg K
! =======================================================================

subroutine setup_mhc_lhc (hydrostatic)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    logical, intent (in) :: hydrostatic
    
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
    
    lv00 = (hlv - d0_vap * tice) / c_air
    li00 = (hlf - dc_ice * tice) / c_air
    li20 = lv00 + li00
    
    c1_vap = c_vap / c_air
    c1_liq = c_liq / c_air
    c1_ice = c_ice / c_air
    
end subroutine setup_mhc_lhc

! =======================================================================
! major cloud microphysics driver
! =======================================================================

subroutine mpdrv (hydrostatic, ua, va, wa, delp, pt, qv, ql, qr, qi, qs, &
        qg, qa, qnl, qni, delz, is, ie, ks, ke, dtm, water, rain, ice, snow, &
        graupel, gsize, hs, q_con, cappa, consv_te, te, prefluxw, prefluxr, &
        prefluxi, prefluxs, prefluxg, condensation, deposition, evaporation, &
        sublimation, last_step, do_inline_mp, do_mp_fast, do_mp_full)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: is, ie, ks, ke
    
    logical, intent (in) :: hydrostatic, last_step, consv_te, do_inline_mp
    logical, intent (in) :: do_mp_fast, do_mp_full
    
    real, intent (in) :: dtm
    
    real, intent (in), dimension (is:ie) :: gsize, hs
    
    real, intent (in), dimension (is:ie, ks:ke) :: delz, qnl, qni
    
    real, intent (inout), dimension (is:ie, ks:ke) :: delp, pt, ua, va, wa
    real, intent (inout), dimension (is:ie, ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    real, intent (inout), dimension (is:ie, ks:ke) :: prefluxw, prefluxr, prefluxi, prefluxs, prefluxg
    
    real, intent (inout), dimension (is:, ks:) :: q_con, cappa
    
    real, intent (inout), dimension (is:ie) :: water, rain, ice, snow, graupel
    real, intent (inout), dimension (is:ie) :: condensation, deposition
    real, intent (inout), dimension (is:ie) :: evaporation, sublimation
    
    real, intent (out), dimension (is:ie, ks:ke) :: te
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i, k, n
    
    real :: rh_adj, rh_rain, ccn0, cin0, cond
    real :: convt, dts, q_cond, t_lnd, t_ocn, h_var, tmp, nl, ni
    
    real, dimension (ks:ke) :: q_liq, q_sol, dp, dz
    real, dimension (ks:ke) :: qvz, qlz, qrz, qiz, qsz, qgz, qaz
    real, dimension (ks:ke) :: den, pz, denfac, ccn, cin
    real, dimension (ks:ke) :: u, v, w
    
    real (kind = r8) :: con_r8, c8
    
    real (kind = r8), dimension (is:ie, ks:ke) :: te_beg, te_end, tw_beg, tw_end
    real (kind = r8), dimension (is:ie, ks:ke) :: te_beg_0, te_end_0, tw_beg_0, tw_end_0
    
    real (kind = r8), dimension (is:ie) :: te_b_beg, te_b_end, tw_b_beg, tw_b_end, dte, te_loss
    real (kind = r8), dimension (is:ie) :: te_b_beg_0, te_b_end_0, tw_b_beg_0, tw_b_end_0
    
    real (kind = r8), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! time steps
    ! -----------------------------------------------------------------------
    
    ntimes = max (ntimes, int (dtm / min (dtm, mp_time)))
    dts = dtm / real (ntimes)
    
    ! -----------------------------------------------------------------------
    ! initialization of total energy difference and condensation diag
    ! -----------------------------------------------------------------------
    
    dte = 0.0
    cond = 0.0
    
    ! -----------------------------------------------------------------------
    ! unit convert to mm/day
    ! -----------------------------------------------------------------------
    
    convt = 86400. * rgrav / dts
    
    do i = is, ie
        
        ! -----------------------------------------------------------------------
        ! conversion of temperature
        ! -----------------------------------------------------------------------
        
        if (do_inline_mp) then
            do k = ks, ke
                q_cond = ql (i, k) + qr (i, k) + qi (i, k) + qs (i, k) + qg (i, k)
                tz (k) = pt (i, k) / ((1. + zvir * qv (i, k)) * (1. - q_cond))
            enddo
        else
            do k = ks, ke
                tz (k) = pt (i, k)
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! calculate base total energy
        ! -----------------------------------------------------------------------
        
        if (consv_te) then
            if (hydrostatic) then
                do k = ks, ke
                    te (i, k) = - c_air * tz (k) * delp (i, k)
                enddo
            else
                do k = ks, ke
                    te (i, k) = - mte (qv (i, k), ql (i, k), qr (i, k), qi (i, k), &
                        qs (i, k), qg (i, k), tz (k), delp (i, k), .true.) * grav
                enddo
            endif
        endif
        
        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------
        
        if (consv_checker) then
            call mtetw (ks, ke, qv (i, :), ql (i, :), qr (i, :), qi (i, :), &
                qs (i, :), qg (i, :), tz, ua (i, :), va (i, :), wa (i, :), &
                delp (i, :), gsize (i), dte (i), water (i), rain (i), ice (i), snow (i), &
                graupel (i), dtm, te_beg_0 (i, :), tw_beg_0 (i, :), &
                te_b_beg_0 (i), tw_b_beg_0 (i), .true., hydrostatic)
        endif
        
        do k = ks, ke
            
            ! -----------------------------------------------------------------------
            ! convert moist mixing ratios to dry mixing ratios
            ! -----------------------------------------------------------------------
            
=======

    do i = is, ie

        do k = ks, ke
            if (do_inline_mp) then
#ifdef MOIST_CAPPA
                tz (k) = pt (i, k) / ((1. + zvir * qv (i, k)) * (1. - (ql (i, k) + qr (i, k) + qi (i, k) + qs (i, k) + qg (i, k))))
#else
                tz (k) = pt (i, k) / (1. + zvir * qv (i, k))
#endif
            else
                tz (k) = pt (i, k)
            endif
        enddo

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
            te_b_beg_0 (i) = (dte (i) - li00 * c_air * (ice (i) + snow (i) + graupel (i)) * dt_in / 86400) * gsize (i) ** 2.0
            tw_b_beg_0 (i) = (rain (i) + ice (i) + snow (i) + graupel (i)) * dt_in / 86400 * gsize (i) ** 2.0
        endif

        do k = ks, ke
            dp0 (k) = delp (i, k)
            ! -----------------------------------------------------------------------
            ! convert moist mixing ratios to dry mixing ratios
            ! -----------------------------------------------------------------------
>>>>>>> user/lnz/shield2022
            qvz (k) = qv (i, k)
            qlz (k) = ql (i, k)
            qrz (k) = qr (i, k)
            qiz (k) = qi (i, k)
            qsz (k) = qs (i, k)
            qgz (k) = qg (i, k)
            qaz (k) = qa (i, k)
            
            q_cond = qlz (k) + qrz (k) + qiz (k) + qsz (k) + qgz (k)
            con_r8 = one_r8 - (qvz (k) + q_cond)
            
            dp (k) = delp (i, k) * con_r8
            con_r8 = one_r8 / con_r8
            qvz (k) = qvz (k) * con_r8
            qlz (k) = qlz (k) * con_r8
            qrz (k) = qrz (k) * con_r8
            qiz (k) = qiz (k) * con_r8
            qsz (k) = qsz (k) * con_r8
            qgz (k) = qgz (k) * con_r8
<<<<<<< HEAD
            
            ! -----------------------------------------------------------------------
            ! dry air density and layer-mean pressure thickness
            ! -----------------------------------------------------------------------
            
            dz (k) = delz (i, k)
            den (k) = - dp (k) / (grav * dz (k))
            pz (k) = den (k) * rdgas * tz (k)
            
=======

            den (k) = - dp1 (k) / (grav * dz1 (k)) ! density of dry air
            p1 (k) = den (k) * rdgas * tz (k) ! dry air pressure

>>>>>>> user/lnz/shield2022
            ! -----------------------------------------------------------------------
            ! for sedi_momentum transport
            ! -----------------------------------------------------------------------
<<<<<<< HEAD
            
            u (k) = ua (i, k)
            v (k) = va (i, k)
=======

            m1 (k) = 0.
            u0 (k) = ua (i, k)
            v0 (k) = va (i, k)
>>>>>>> user/lnz/shield2022
            if (.not. hydrostatic) then
                w (k) = wa (i, k)
            endif
            
        enddo
<<<<<<< HEAD
        
        do k = ks, ke
            denfac (k) = sqrt (den (ke) / den (k))
        enddo
        
=======

        ! -----------------------------------------------------------------------
        ! fix energy conservation
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

>>>>>>> user/lnz/shield2022
        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            call mtetw (ks, ke, qvz, qlz, qrz, qiz, qsz, qgz, tz, u, v, w, &
                dp, gsize (i), dte (i), water (i), rain (i), ice (i), snow (i), &
                graupel (i), dtm, te_beg (i, :), tw_beg (i, :), &
                te_b_beg (i), tw_b_beg (i), .false., hydrostatic)
        endif

        ! -----------------------------------------------------------------------
        ! cloud condensation nuclei (CCN), cloud ice nuclei (CIN)
        ! -----------------------------------------------------------------------
<<<<<<< HEAD
        
=======

        cpaut = c_paut * 0.104 * grav / 1.717e-5

>>>>>>> user/lnz/shield2022
        if (prog_ccn) then
            do k = ks, ke
                ! boucher and lohmann (1995)
                nl = min (1., abs (hs (i)) / (10. * grav)) * &
                     (10. ** 2.24 * (qnl (i, k) * den (k) * 1.e9) ** 0.257) + &
                     (1. - min (1., abs (hs (i)) / (10. * grav))) * &
                     (10. ** 2.06 * (qnl (i, k) * den (k) * 1.e9) ** 0.48)
                ni = qni (i, k)
                ccn (k) = max (10.0, nl) * 1.e6
                cin (k) = max (10.0, ni) * 1.e6
                ccn (k) = ccn (k) / den (k)
                cin (k) = cin (k) / den (k)
            enddo
        else
            ccn0 = (ccn_l * min (1., abs (hs (i)) / (10. * grav)) + &
                ccn_o * (1. - min (1., abs (hs (i)) / (10. * grav)))) * 1.e6
            cin0 = 0.0
            do k = ks, ke
                ccn (k) = ccn0 / den (k)
                cin (k) = cin0 / den (k)
            enddo
        endif

        ! -----------------------------------------------------------------------
        ! subgrid deviation in horizontal direction
        ! default area dependent form: use dx ~ 100 km as the base
        ! -----------------------------------------------------------------------
<<<<<<< HEAD
        
        t_lnd = dw_land * sqrt (gsize (i) / 1.e5)
        t_ocn = dw_ocean * sqrt (gsize (i) / 1.e5)
=======

        s_leng = sqrt (gsize (i) / 1.e5)
        t_land = dw_land * s_leng
        t_ocean = dw_ocean * s_leng
>>>>>>> user/lnz/shield2022
        tmp = min (1., abs (hs (i)) / (10. * grav))
        h_var = t_lnd * tmp + t_ocn * (1. - tmp)
        h_var = min (0.20, max (0.01, h_var))

        ! -----------------------------------------------------------------------
        ! relative humidity thresholds
        ! -----------------------------------------------------------------------

        rh_adj = 1. - h_var - rh_inc
<<<<<<< HEAD
        rh_rain = max (0.35, rh_adj - rh_inr)
        
=======
        rh_rain = max (0.35, rh_adj - rh_inr) ! rh_inr = 0.25

>>>>>>> user/lnz/shield2022
        ! -----------------------------------------------------------------------
        ! fix negative water species from outside
        ! -----------------------------------------------------------------------

        if (fix_negative) &
<<<<<<< HEAD
            call neg_adj (ks, ke, tz, dp, qvz, qlz, qrz, qiz, qsz, qgz, cond)
        
        condensation (i) = condensation (i) + cond * convt * ntimes
        
        ! -----------------------------------------------------------------------
        ! fast microphysics loop
        ! -----------------------------------------------------------------------
        
        if (do_mp_fast) then
            
            call mp_fast (ks, ke, tz, qvz, qlz, qrz, qiz, qsz, qgz, dtm, dp, den, &
                ccn, cin, condensation (i), deposition (i), evaporation (i), &
                sublimation (i), convt)
            
        endif
        
        ! -----------------------------------------------------------------------
        ! full microphysics loop
        ! -----------------------------------------------------------------------
        
        if (do_mp_full) then
=======
            call neg_adj (ks, ke, tz, dp1, qvz, qlz, qrz, qiz, qsz, qgz, cond)

        condensation (i) = condensation (i) + cond * convt * ntimes

        m2_rain (i, :) = 0.
        m2_sol (i, :) = 0.

        do n = 1, ntimes

            ! -----------------------------------------------------------------------
            ! time - split warm rain processes: 1st pass
            ! -----------------------------------------------------------------------

            call warm_rain (dt_rain, ks, ke, dp1, dz1, tz, qvz, qlz, qrz, qiz, qsz, &
                qgz, den, denfac, ccn, c_praut, rh_rain, vtrz, r1, pfr, m1_rain, w1, h_var, reevap, dte (i))
            
            evaporation (i) = evaporation (i) + reevap * convt
            rain (i) = rain (i) + r1 * convt
            prefluxr (i, :) = prefluxr (i, :) + pfr * convt
            
            do k = ks, ke
                m2_rain (i, k) = m2_rain (i, k) + m1_rain (k)
                m1 (k) = m1 (k) + m1_rain (k)
            enddo

            ! -----------------------------------------------------------------------
            ! sedimentation of cloud ice, snow, and graupel
            ! -----------------------------------------------------------------------

            call fall_speed (ks, ke, den, qsz, qiz, qgz, qlz, tz, vtsz, vtiz, vtgz)

            call terminal_fall (dts, ks, ke, tz, qvz, qlz, qrz, qgz, qsz, qiz, &
                dz1, dp1, den, vtgz, vtsz, vtiz, r1, g1, s1, i1, pfr, pfi, pfs, pfg, m1_sol, w1, dte (i))
>>>>>>> user/lnz/shield2022
            
            call mp_full (ks, ke, ntimes, tz, qvz, qlz, qrz, qiz, qsz, qgz, dp, dz, &
                u, v, w, den, denfac, ccn, cin, dts, rh_adj, rh_rain, h_var, dte (i), &
                water (i), rain (i), ice (i), snow (i), graupel (i), prefluxw (i, :), &
                prefluxr (i, :), prefluxi (i, :), prefluxs (i, :), prefluxg (i, :), &
                condensation (i), deposition (i), evaporation (i), sublimation (i), convt)
            
<<<<<<< HEAD
        endif
        
        ! -----------------------------------------------------------------------
        ! cloud fraction diagnostic
        ! -----------------------------------------------------------------------
        
        if (do_qa .and. last_step) then
            call cloud_fraction (ks, ke, pz, den, qvz, qlz, qrz, qiz, qsz, qgz, qaz, &
                tz, h_var, gsize (i))
        endif
        
=======
            ! -----------------------------------------------------------------------
            ! energy loss during sedimentation heating
            ! -----------------------------------------------------------------------

            if (consv_checker) then
                do k = ks, ke
                    te1 (k) = one_r8 + qvz (k) * c1_vap + (qlz (k) + qrz (k)) * c1_liq + (qiz (k) + qsz (k) + qgz (k)) * c1_ice
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
                    te2 (k) = one_r8 + qvz (k) * c1_vap + (qlz (k) + qrz (k)) * c1_liq + (qiz (k) + qsz (k) + qgz (k)) * c1_ice
                    te2 (k) = rgrav * te2 (k) * c_air * tz (k) * dp1 (k)
                enddo
                dte (i) = dte (i) + sum (te1) - sum (te2)
            endif

            ! -----------------------------------------------------------------------
            ! time - split warm rain processes: 2nd pass
            ! -----------------------------------------------------------------------

            call warm_rain (dt_rain, ks, ke, dp1, dz1, tz, qvz, qlz, qrz, qiz, qsz, &
                qgz, den, denfac, ccn, c_praut, rh_rain, vtrz, r1, pfr, m1_rain, w1, h_var, reevap, dte (i))
            
            evaporation (i) = evaporation (i) + reevap * convt
            rain (i) = rain (i) + r1 * convt
            prefluxr (i, :) = prefluxr (i, :) + pfr * convt
            
            do k = ks, ke
                m2_rain (i, k) = m2_rain (i, k) + m1_rain (k)
                m2_sol (i, k) = m2_sol (i, k) + m1_sol (k)
                m1 (k) = m1 (k) + m1_rain (k) + m1_sol (k)
            enddo

            ! -----------------------------------------------------------------------
            ! ice - phase microphysics
            ! -----------------------------------------------------------------------

            call icloud (ks, ke, tz, p1, qvz, qlz, qrz, qiz, qsz, qgz, dp1, den, ccn, &
                cin, denfac, vtsz, vtgz, vtrz, qaz, rh_adj, rh_rain, dts, h_var, gsize (i), &
                cond, dep, reevap, sub, last_step)

            condensation (i) = condensation (i) + cond * convt
            deposition (i) = deposition (i) + dep * convt
            evaporation (i) = evaporation (i) + reevap * convt
            sublimation (i) = sublimation (i) + sub * convt

        enddo

>>>>>>> user/lnz/shield2022
        ! -----------------------------------------------------------------------
        ! momentum transportation during sedimentation
        ! -----------------------------------------------------------------------
<<<<<<< HEAD
        
        if (do_sedi_uv) then
=======

        if (sedi_transport) then
>>>>>>> user/lnz/shield2022
            do k = ks + 1, ke
                c8 = mhc (qvz (k), qlz (k), qrz (k), qiz (k), qsz (k), qgz (k)) * c_air
                tz (k) = tz (k) + 0.5 * (ua (i, k) ** 2 + va (i, k) ** 2 - &
                     (u (k) ** 2 + v (k) ** 2)) / c8
            enddo
            do k = ks + 1, ke
                ua (i, k) = u (k)
                va (i, k) = v (k)
            enddo
        endif

        if (do_sedi_w) then
            do k = ks, ke
                c8 = mhc (qvz (k), qlz (k), qrz (k), qiz (k), qsz (k), qgz (k)) * c_air
                tz (k) = tz (k) + 0.5 * (wa (i, k) ** 2 - w (k) ** 2) / c8
            enddo
            do k = ks, ke
                wa (i, k) = w (k)
            enddo
        endif

        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            call mtetw (ks, ke, qvz, qlz, qrz, qiz, qsz, qgz, tz, u, v, w, &
                dp, gsize (i), dte (i), water (i), rain (i), ice (i), snow (i), &
                graupel (i), dtm, te_end (i, :), tw_end (i, :), &
                te_b_end (i), tw_b_end (i), .false., hydrostatic, te_loss (i))
        endif
<<<<<<< HEAD
        
=======

        ! -----------------------------------------------------------------------
        ! update moist air mass (actually hydrostatic pressure)
        ! convert to dry mixing ratios
        ! -----------------------------------------------------------------------

>>>>>>> user/lnz/shield2022
        do k = ks, ke
            
            ! -----------------------------------------------------------------------
            ! convert dry mixing ratios back to moist mixing ratios
            ! -----------------------------------------------------------------------
            
            q_cond = qlz (k) + qrz (k) + qiz (k) + qsz (k) + qgz (k)
            con_r8 = one_r8 + qvz (k) + q_cond
            
            delp (i, k) = dp (k) * con_r8
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
            qa (i, k) = qaz (k)
            
            ! -----------------------------------------------------------------------
            ! calculate some more variables needed outside
            ! -----------------------------------------------------------------------
            
            q_liq (k) = qlz (k) + qrz (k)
            q_sol (k) = qiz (k) + qsz (k) + qgz (k)
            q_cond = q_liq (k) + q_sol (k)
            con_r8 = one_r8 - (qvz (k) + q_cond)
            c8 = mhc (con_r8, qvz (k), q_liq (k), q_sol (k)) * c_air
            
#ifdef USE_COND
            q_con (i, k) = q_cond
#endif
#ifdef MOIST_CAPPA
            tmp = rdgas * (1. + zvir * qvz (k))
            cappa (i, k) = tmp / (tmp + c8)
#endif
            
        enddo

        ! -----------------------------------------------------------------------
        ! total energy checker
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            call mtetw (ks, ke, qv (i, :), ql (i, :), qr (i, :), qi (i, :), &
                qs (i, :), qg (i, :), tz, ua (i, :), va (i, :), wa (i, :), &
                delp (i, :), gsize (i), dte (i), water (i), rain (i), ice (i), snow (i), &
                graupel (i), dtm, te_end_0 (i, :), tw_end_0 (i, :), &
                te_b_end_0 (i), tw_b_end_0 (i), .true., hydrostatic)
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
                    te (i, k) = te (i, k) + mte (qv (i, k), ql (i, k), qr (i, k), qi (i, k), &
                        qs (i, k), qg (i, k), tz (k), delp (i, k), .true.) * grav
                enddo
            endif
        endif

        ! -----------------------------------------------------------------------
        ! conversion of temperature
        ! -----------------------------------------------------------------------
<<<<<<< HEAD
        
        if (do_inline_mp) then
            do k = ks, ke
                q_cond = qlz (k) + qrz (k) + qiz (k) + qsz (k) + qgz (k)
                pt (i, k) = tz (k) * (1. + zvir * qvz (k)) * (1. - q_cond)
            enddo
        else
            do k = ks, ke
                q_liq (k) = qlz (k) + qrz (k)
                q_sol (k) = qiz (k) + qsz (k) + qgz (k)
                q_cond = q_liq (k) + q_sol (k)
                con_r8 = one_r8 - (qvz (k) + q_cond)
                c8 = mhc (con_r8, qvz (k), q_liq (k), q_sol (k)) * c_air
                pt (i, k) = pt (i, k) + (tz (k) - pt (i, k)) * c8 / cp_air
            enddo
        endif
        
    enddo ! i loop
    
=======

        do k = ks, ke
            qa (i, k) = qaz (k)
        enddo

    enddo

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! total energy checker
    ! -----------------------------------------------------------------------

    if (consv_checker) then
        if (abs (sum (te_end) + sum (te_b_end) - sum (te_beg) - sum (te_b_beg)) / &
             (sum (te_beg) + sum (te_b_beg)) .gt. te_err) then
            print *, "GFDL MP TE: ", (sum (te_beg) + sum (te_b_beg)) / sum (gsize ** 2), &
                (sum (te_end) + sum (te_b_end)) / sum (gsize ** 2), &
                (sum (te_end) + sum (te_b_end) - sum (te_beg) - sum (te_b_beg)) / &
                (sum (te_beg) + sum (te_b_beg))
        endif
        if (abs (sum (tw_end) + sum (tw_b_end) - sum (tw_beg) - sum (tw_b_beg)) / &
             (sum (tw_beg) + sum (tw_b_beg)) .gt. te_err) then
            print *, "GFDL MP TW: ", (sum (tw_beg) + sum (tw_b_beg)) / sum (gsize ** 2), &
                (sum (tw_end) + sum (tw_b_end)) / sum (gsize ** 2), &
                (sum (tw_end) + sum (tw_b_end) - sum (tw_beg) - sum (tw_b_beg)) / &
                (sum (tw_beg) + sum (tw_b_beg))
        endif
        ! print *, "GFDL MP TE LOSS (%) : ", sum (te_loss) / (sum (te_beg) + sum (te_b_beg)) * 100.0
    endif

end subroutine mpdrv

! =======================================================================
! fix negative water species
! =======================================================================

subroutine neg_adj (ks, ke, tz, dp, qv, ql, qr, qi, qs, qg, cond)
    
    implicit none
<<<<<<< HEAD
    
=======
    integer, intent (in) :: ks, ke
    real, intent (in), dimension (ks:ke) :: dm, m1, dz, qv, ql, qr, qi, qs, qg
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (in) :: cw ! heat capacity
    ! local:
    real, dimension (ks:ke) :: dgz, cv0
    integer :: k

    ! this is the vectorized loop
    do k = ks + 1, ke
        dgz (k) = - g2 * (dz (k - 1) + dz (k))
        cv0 (k) = dm (k) * (cv_air + qv (k) * cv_vap + (qr (k) + ql (k)) * c_liq + &
             (qi (k) + qs (k) + qg (k)) * c_ice) + cw * (m1 (k) - m1 (k - 1))
        ! cvm_new + cw * m1 (k) = cvm_old + cw * m1 (k - 1)
    enddo
>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    integer, intent (in) :: ks, ke
    
    real, intent (in), dimension (ks:ke) :: dp
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real, intent (out) :: cond
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: dq, sink
    
    real, dimension (ks:ke) :: q_liq, q_sol, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), dimension (ks:ke) :: cvm, te8
=======
    ! top layer: cv0 = cvn + cw * m1 (k)
    ! tz (k) = cv0 (k) * tz (k) / (cvn (k) + cw * m1 (k)) = tz (k) -- > no change
    do k = ks + 1, ke
        tz (k) = (cv0 (k) * tz (k) + m1 (k - 1) * (cw * tz (k - 1) + dgz (k))) / (cv0 (k) + cw * m1 (k - 1))
    enddo

end subroutine sedi_heat

! -----------------------------------------------------------------------
! warm rain cloud microphysics
! -----------------------------------------------------------------------

subroutine warm_rain (dt, ks, ke, dp, dz, tz, qv, ql, qr, qi, qs, qg, &
        den, denfac, ccn, c_praut, rh_rain, vtr, r1, pfr, m1_rain, w1, h_var, reevap, dte)
    
    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: dt ! time step (s)
    real, intent (in) :: rh_rain, h_var
    real, intent (in), dimension (ks:ke) :: dp, dz, den
    real, intent (in), dimension (ks:ke) :: denfac, ccn, c_praut

    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (inout), dimension (ks:ke) :: vtr, qv, ql, qr, qi, qs, qg, m1_rain, w1, pfr
    real (kind = r_grid), intent (inout) :: dte
    real, intent (out) :: r1
    real, intent (out) :: reevap
    real, parameter :: so3 = 7. / 3.
    ! fall velocity constants:
    real, parameter :: vconr = 2503.23638966667
    real, parameter :: normr = 25132741228.7183
    real, parameter :: thr = 1.e-8

    real, dimension (ks:ke) :: dl, dm
    real (kind = r_grid), dimension (ks:ke) :: te1, te2
    real, dimension (ks:ke + 1) :: ze, zt
    real :: sink, dq, qc
    real :: qden
    real :: zs = 0.
    real :: dt5
    integer :: k

    logical :: no_fall

    dt5 = 0.5 * dt
    pfr = 0.
>>>>>>> user/lnz/shield2022
    
    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    cond = 0
    
    ! -----------------------------------------------------------------------
    ! calculate moist heat capacity and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    call cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, cvm, te8, tz, &
        lcpk, icpk, tcpk, tcp3)
    
    do k = ks, ke
        
=======

    m1_rain (:) = 0.

    call check_column (ks, ke, qr, no_fall)

    reevap = 0

    if (no_fall) then
        vtr (:) = vf_min
        r1 = 0.
    else

>>>>>>> user/lnz/shield2022
        ! -----------------------------------------------------------------------
        ! fix negative solid-phase hydrometeors
        ! -----------------------------------------------------------------------
<<<<<<< HEAD
        
        ! if cloud ice < 0, borrow from snow
        if (qi (k) .lt. 0.) then
            sink = min (- qi (k), max (0., qs (k)))
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., sink, - sink, 0.)
        endif
        
        ! if snow < 0, borrow from graupel
        if (qs (k) .lt. 0.) then
            sink = min (- qs (k), max (0., qg (k)))
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., 0., sink, - sink)
        endif
        
        ! if graupel < 0, borrow from rain
        if (qg (k) .lt. 0.) then
            sink = min (- qg (k), max (0., qr (k)))
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., - sink, 0., 0., sink, te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
=======

        if (const_vr) then
            vtr (:) = vr_fac ! ifs_2016: 4.0
        else
            do k = ks, ke
                qden = qr (k) * den (k)
                if (qr (k) < thr) then
                    vtr (k) = vr_min
                else
                    vtr (k) = vr_fac * vconr * sqrt (min (10., sfcrho / den (k))) * &
                        exp (0.2 * log (qden / normr))
                    vtr (k) = min (vr_max, max (vr_min, vtr (k)))
                endif
            enddo
        endif

        ze (ke + 1) = zs
        do k = ke, ks, - 1
            ze (k) = ze (k + 1) - dz (k) ! dz < 0
        enddo

        ! -----------------------------------------------------------------------
        ! evaporation and accretion of rain for the first 1 / 2 time step
        ! -----------------------------------------------------------------------

        call revap_racc (ks, ke, dt5, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var, dp, reevap)

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
>>>>>>> user/lnz/shield2022
        endif

        ! -----------------------------------------------------------------------
        ! fix negative liquid-phase hydrometeors
        ! -----------------------------------------------------------------------
<<<<<<< HEAD
        
        ! if rain < 0, borrow from cloud water
        if (qr (k) .lt. 0.) then
            sink = min (- qr (k), max (0., ql (k)))
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - sink, sink, 0., 0., 0.)
        endif
        
        ! if cloud water < 0, borrow from water vapor
        if (ql (k) .lt. 0.) then
            sink = min (- ql (k), max (0., qv (k)))
            cond = cond + sink * dp (k)
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                 - sink, sink, 0., 0., 0., 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
        endif
        
    enddo
    
    ! -----------------------------------------------------------------------
    ! fix negative water vapor
    ! -----------------------------------------------------------------------
    
    ! if water vapor < 0, borrow water vapor from below
    do k = ks, ke - 1
        if (qv (k) .lt. 0.) then
            qv (k + 1) = qv (k + 1) + qv (k) * dp (k) / dp (k + 1)
            qv (k) = 0.
        endif
    enddo
    
    ! if water vapor < 0, borrow water vapor from above
    if (qv (ke) .lt. 0. .and. qv (ke - 1) .gt. 0.) then
        dq = min (- qv (ke) * dp (ke), qv (ke - 1) * dp (ke - 1))
        qv (ke - 1) = qv (ke - 1) - dq / dp (ke - 1)
        qv (ke) = qv (ke) + dq / dp (ke)
    endif
    
end subroutine neg_adj

! =======================================================================
! full microphysics loop
! =======================================================================

subroutine mp_full (ks, ke, ntimes, tz, qv, ql, qr, qi, qs, qg, dp, dz, u, v, w, &
        den, denfac, ccn, cin, dts, rh_adj, rh_rain, h_var, dte, water, rain, ice, &
        snow, graupel, prefluxw, prefluxr, prefluxi, prefluxs, prefluxg, &
        condensation, deposition, evaporation, sublimation, convt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke, ntimes
    
    real, intent (in) :: dts, rh_adj, rh_rain, h_var, convt
    
    real, intent (in), dimension (ks:ke) :: dp, dz, den, denfac
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, u, v, w, ccn, cin
    real, intent (inout), dimension (ks:ke) :: prefluxw, prefluxr, prefluxi, prefluxs, prefluxg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    real, intent (inout) :: water, rain, ice, snow, graupel
    real, intent (inout) :: condensation, deposition
    real, intent (inout) :: evaporation, sublimation
    
    real (kind = r8), intent (inout) :: dte
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: n
    
    real :: w1, r1, i1, s1, g1, cond, dep, reevap, sub
    
    real, dimension (ks:ke) :: vtw, vtr, vti, vts, vtg, pfw, pfr, pfi, pfs, pfg
    
    do n = 1, ntimes
        
=======

        if (use_ppm) then
            zt (ks) = ze (ks)
            do k = ks + 1, ke
                zt (k) = ze (k) - dt5 * (vtr (k - 1) + vtr (k))
            enddo
            zt (ke + 1) = zs - dt * vtr (ke)

            do k = ks, ke
                if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
            enddo
            call lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, qr, r1, m1_rain, mono_prof)
        else
            call implicit_fall (dt, ks, ke, ze, vtr, dp, qr, r1, m1_rain)
        endif
        
        pfr (ks) = max (0.0, m1_rain (ks))
        do k = ke, ks + 1, -1
            pfr (k) = max (0.0, m1_rain (k) - m1_rain (k - 1))
        enddo

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

>>>>>>> user/lnz/shield2022
        ! -----------------------------------------------------------------------
        ! sedimentation of cloud ice, snow, graupel or hail, and rain
        ! -----------------------------------------------------------------------
<<<<<<< HEAD
        
        call sedimentation (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, &
            dz, dp, vtw, vtr, vti, vts, vtg, w1, r1, i1, s1, g1, pfw, pfr, pfi, pfs, pfg, &
            u, v, w, den, denfac, dte)
        
        water = water + w1 * convt
        rain = rain + r1 * convt
        ice = ice + i1 * convt
        snow = snow + s1 * convt
        graupel = graupel + g1 * convt
        
        prefluxw = prefluxw + pfw * convt
        prefluxr = prefluxr + pfr * convt
        prefluxi = prefluxi + pfi * convt
        prefluxs = prefluxs + pfs * convt
        prefluxg = prefluxg + pfg * convt
        
=======

        if (do_sedi_w) then
            ! conservation of vertical momentum:
            w1 (ks) = w1 (ks) + m1_rain (ks) * vtr (ks) / dm (ks)
            do k = ks + 1, ke
                w1 (k) = (dm (k) * w1 (k) + m1_rain (k - 1) * (w1 (k - 1) - vtr (k - 1)) + m1_rain (k) * vtr (k)) &
                     / (dm (k) + m1_rain (k - 1))
            enddo
        endif

        ! -----------------------------------------------------------------------
        ! energy loss during sedimentation heating
        ! -----------------------------------------------------------------------

        if (consv_checker) then
            do k = ks, ke
                te1 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te1 (k) = rgrav * te1 (k) * c_air * tz (k) * dp (k)
            enddo
        endif

>>>>>>> user/lnz/shield2022
        ! -----------------------------------------------------------------------
        ! warm rain cloud microphysics
        ! -----------------------------------------------------------------------
<<<<<<< HEAD
        
        call warm_rain (dts, ks, ke, dp, dz, tz, qv, ql, qr, qi, qs, qg, &
            den, denfac, vtw, vtr, ccn, rh_rain, h_var, reevap)
        
        evaporation = evaporation + reevap * convt
        
=======

        if (do_sedi_heat) then
            call sedi_heat (ks, ke, dp, m1_rain, dz, tz, qv, ql, qr, qi, qs, qg, c_liq)
        endif

>>>>>>> user/lnz/shield2022
        ! -----------------------------------------------------------------------
        ! ice cloud microphysics
        ! -----------------------------------------------------------------------
<<<<<<< HEAD
        
        call ice_cloud (ks, ke, tz, qv, ql, qr, qi, qs, qg, den, &
            denfac, vtw, vtr, vti, vts, vtg, dts, h_var)
        
=======

        if (consv_checker) then
            do k = ks, ke
                te2 (k) = one_r8 + qv (k) * c1_vap + (ql (k) + qr (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
                te2 (k) = rgrav * te2 (k) * c_air * tz (k) * dp (k)
            enddo
            dte = dte + sum (te1) - sum (te2)
        endif


>>>>>>> user/lnz/shield2022
        ! -----------------------------------------------------------------------
        ! temperature sentive high vertical resolution processes
        ! -----------------------------------------------------------------------
<<<<<<< HEAD
        
        call subgrid_z_proc (ks, ke, den, denfac, dts, rh_adj, tz, qv, ql, &
            qr, qi, qs, qg, dp, ccn, cin, cond, dep, reevap, sub)
        
        condensation = condensation + cond * convt
        deposition = deposition + dep * convt
        evaporation = evaporation + reevap * convt
        sublimation = sublimation + sub * convt
        
    enddo
    
end subroutine mp_full

! =======================================================================
! fast microphysics loop
! =======================================================================

subroutine mp_fast (ks, ke, tz, qv, ql, qr, qi, qs, qg, dtm, dp, den, &
        ccn, cin, condensation, deposition, evaporation, sublimation, convt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dtm, convt
    
    real, intent (in), dimension (ks:ke) :: dp, den
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, ccn, cin
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    real, intent (inout) :: condensation, deposition
    real, intent (inout) :: evaporation, sublimation
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real :: cond, dep, reevap, sub
    
    real, dimension (ks:ke) :: q_liq, q_sol, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), dimension (ks:ke) :: cvm, te8
    
    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
    
    cond = 0
    dep = 0
    reevap = 0
    sub = 0
    
    ! -----------------------------------------------------------------------
    ! calculate heat capacities and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    call cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, cvm, te8, tz, &
        lcpk, icpk, tcpk, tcp3)
    
    if (.not. do_warm_rain_mp) then
        
=======

        call revap_racc (ks, ke, dt5, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var, dp, reevap)

    endif

    ! -----------------------------------------------------------------------
    ! auto - conversion
    ! assuming linear subgrid vertical distribution of cloud water
    ! following lin et al. 1994, mwr
    ! -----------------------------------------------------------------------

    if (irain_f /= 0) then

>>>>>>> user/lnz/shield2022
        ! -----------------------------------------------------------------------
        ! cloud ice melting to form cloud water and rain
        ! -----------------------------------------------------------------------
<<<<<<< HEAD
        
        call pimlt (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, cvm, te8, &
            lcpk, icpk, tcpk, tcp3)
        
=======

        do k = ks, ke
            qc = fac_rc * ccn (k)
            if (tz (k) > t_wfr) then
                dq = ql (k) - qc
                if (dq > 0.) then
                    sink = min (dq, dt * c_praut (k) * den (k) * exp (so3 * log (ql (k))))
                    ql (k) = ql (k) - sink
                    qr (k) = qr (k) + sink
                endif
            endif
        enddo

    else

>>>>>>> user/lnz/shield2022
        ! -----------------------------------------------------------------------
        ! enforce complete freezing below t_wfr
        ! -----------------------------------------------------------------------
<<<<<<< HEAD
        
        call pcomp (ks, ke, qv, ql, qr, qi, qs, qg, tz, cvm, te8, &
            lcpk, icpk, tcpk, tcp3)
        
    endif
    
    ! -----------------------------------------------------------------------
    ! cloud water condensation and evaporation
    ! -----------------------------------------------------------------------
    
    call pcond_pevap (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
        lcpk, icpk, tcpk, tcp3, cond, reevap)
    
    condensation = condensation + cond * convt
    evaporation = evaporation + reevap * convt
        
    if (.not. do_warm_rain_mp) then
 
        ! -----------------------------------------------------------------------
        ! cloud water freezing to form cloud ice and snow
        ! -----------------------------------------------------------------------
        
        call pifr (ks, ke, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, &
            lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! Wegener Bergeron Findeisen process
        ! -----------------------------------------------------------------------
        
        call pwbf (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, &
            lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! Bigg freezing mechanism
        ! -----------------------------------------------------------------------
        
        call pbigg (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, ccn, &
            lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! rain freezing to form graupel
        ! -----------------------------------------------------------------------
        
        call pgfr_simp (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, cvm, te8, &
            lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! snow melting to form cloud water and rain
        ! -----------------------------------------------------------------------
        
        call psmlt_simp (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, cvm, te8, &
            lcpk, icpk, tcpk, tcp3)
        
    endif
    
    ! -----------------------------------------------------------------------
    ! cloud water to rain autoconversion
    ! -----------------------------------------------------------------------
    
    call praut_simp (ks, ke, dtm, tz, qv, ql, qr, qi, qs, qg)
    
    if (.not. do_warm_rain_mp) then
        
=======

        call linear_prof (ke - ks + 1, ql (ks), dl (ks), z_slope_liq, h_var)

        do k = ks, ke
            qc = fac_rc * ccn (k)
            if (tz (k) > t_wfr + dt_fr) then
                dl (k) = min (max (1.e-6, dl (k)), 0.5 * ql (k))
                ! --------------------------------------------------------------------
                ! as in klein's gfdl am2 stratiform scheme (with subgrid variations)
                ! --------------------------------------------------------------------
                dq = 0.5 * (ql (k) + dl (k) - qc)
                ! --------------------------------------------------------------------
                ! dq = dl if qc == q_minus = ql - dl
                ! dq = 0 if qc == q_plus = ql + dl
                ! --------------------------------------------------------------------
                if (dq > 0.) then ! q_plus > qc
                    ! --------------------------------------------------------------------
                    ! revised continuous form: linearly decays (with subgrid dl) to zero at qc == ql + dl
                    ! --------------------------------------------------------------------
                    sink = min (1., dq / dl (k)) * dt * c_praut (k) * den (k) * exp (so3 * log (ql (k)))
                    ql (k) = ql (k) - sink
                    qr (k) = qr (k) + sink
                endif
            endif
        enddo
    endif

end subroutine warm_rain

! -----------------------------------------------------------------------
! evaporation of rain
! -----------------------------------------------------------------------

subroutine revap_racc (ks, ke, dt, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var, dp, reevap)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: dt ! time step (s)
    real, intent (in) :: rh_rain, h_var
    real, intent (in), dimension (ks:ke) :: den, denfac, dp
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (inout), dimension (ks:ke) :: qv, qr, ql, qi, qs, qg
    real, intent (out) :: reevap
    ! local:
    real (kind = r_grid), dimension (ks:ke) :: cvm
    real, dimension (ks:ke) :: q_liq, q_sol, lcpk
    real :: dqv, qsat, dqsdt, evap, t2, qden, q_plus, q_minus, sink
    real :: qpz, dq, dqh, tin
    real :: fac_revp, rh_tem

    integer :: k

    if (tau_revp .gt. 1.e-6) then
        fac_revp = 1. - exp (- dt / tau_revp)
    else
        fac_revp = 1.
    endif

    do k = ks, ke

        if (tz (k) > t_wfr .and. qr (k) > qrmin) then

            ! -----------------------------------------------------------------------
            ! define heat capacity and latent heat coefficient
            ! -----------------------------------------------------------------------

            q_liq (k) = ql (k) + qr (k)
            q_sol (k) = qi (k) + qs (k) + qg (k)

            cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
            lcpk (k) = (lv00 + d1_vap * tz (k)) / cvm (k)
            tin = (tz (k) * cvm (k) - lv00 * ql (k)) / (1. + (qv (k) + ql (k)) * c1_vap + qr (k) * c1_liq + q_sol (k) * c1_ice)

            qpz = qv (k) + ql (k)
            qsat = wqs2 (tin, den (k), dqsdt)
            dqh = max (ql (k), h_var * max (qpz, qcmin))
            dqh = min (dqh, 0.2 * qpz) ! new limiter
            dqv = qsat - qv (k) ! use this to prevent super - sat the gird box
            q_minus = qpz - dqh
            q_plus = qpz + dqh

            ! -----------------------------------------------------------------------
            ! qsat must be > q_minus to activate evaporation
            ! qsat must be < q_plus to activate accretion
            ! -----------------------------------------------------------------------

            ! -----------------------------------------------------------------------
            ! rain evaporation
            ! -----------------------------------------------------------------------

            rh_tem = qpz / iqs1 (tin, den (k))

            if (dqv > qvmin .and. qsat > q_minus) then
                if (qsat > q_plus) then
                    dq = qsat - qpz
                else
                    ! -----------------------------------------------------------------------
                    ! q_minus < qsat < q_plus
                    ! dq == dqh if qsat == q_minus
                    ! -----------------------------------------------------------------------
                    dq = 0.25 * (qsat - q_minus) ** 2 / dqh
                endif
                qden = qr (k) * den (k)
                t2 = tin * tin
                if (use_rhc_revap) then
                    evap = 0.0
                    if (rh_tem < rhc_revap) then
                        evap = crevp (1) * t2 * dq * (crevp (2) * sqrt (qden) + crevp (3) * &
                            exp (0.725 * log (qden)) * sqrt (denfac (k))) / (crevp (4) * t2 + crevp (5) * qsat * den (k))
                        evap = min (qr (k), dt * fac_revp * evap, dqv / (1. + lcpk (k) * dqsdt))
                    endif
                else
                    evap = crevp (1) * t2 * dq * (crevp (2) * sqrt (qden) + crevp (3) * &
                        exp (0.725 * log (qden))) / (crevp (4) * t2 + crevp (5) * qsat * den (k))
                    evap = min (qr (k), dt * fac_revp * evap, dqv / (1. + lcpk (k) * dqsdt))
                endif
                reevap = reevap + evap * dp (k)

                ! -----------------------------------------------------------------------
                ! alternative minimum evap in dry environmental air
                ! sink = min (qr (k), dim (rh_rain * qsat, qv (k)) / (1. + lcpk (k) * dqsdt))
                ! evap = max (evap, sink)
                ! -----------------------------------------------------------------------
                qr (k) = qr (k) - evap
                qv (k) = qv (k) + evap
                q_liq (k) = q_liq (k) - evap
                tz (k) = (cvm (k) * tz (k) - lv00 * evap) / (one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice)
            endif

            ! -----------------------------------------------------------------------
            ! accretion: pracc
            ! -----------------------------------------------------------------------

            ! if (qr (k) > qrmin .and. ql (k) > 1.e-7 .and. qsat < q_plus) then
            if (qr (k) > qrmin .and. ql (k) > 1.e-6 .and. qsat < q_minus) then
                sink = dt * denfac (k) * cracw * exp (0.95 * log (qr (k) * den (k)))
                sink = sink / (1. + sink) * ql (k)
                ql (k) = ql (k) - sink
                qr (k) = qr (k) + sink
            endif

        endif ! warm - rain
    enddo

end subroutine revap_racc

! -----------------------------------------------------------------------
! definition of vertical subgrid variability
! used for cloud ice and cloud water autoconversion
! qi -- > ql & ql -- > qr
! edges: qe == qbar + / - dm
! -----------------------------------------------------------------------

subroutine linear_prof (km, q, dm, z_var, h_var)

    implicit none

    integer, intent (in) :: km
    real, intent (in) :: q (km), h_var
    real, intent (out) :: dm (km)
    logical, intent (in) :: z_var
    real :: dq (km)
    integer :: k

    if (z_var) then
        do k = 2, km
            dq (k) = 0.5 * (q (k) - q (k - 1))
        enddo
        dm (1) = 0.

>>>>>>> user/lnz/shield2022
        ! -----------------------------------------------------------------------
        ! cloud ice deposition and sublimation
        ! -----------------------------------------------------------------------
<<<<<<< HEAD
        
        call pidep_pisub (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
            lcpk, icpk, tcpk, tcp3, cin, dep, sub)
        
        deposition = deposition + dep * convt
        sublimation = sublimation + sub * convt
        
=======

        do k = 2, km - 1
            dm (k) = 0.5 * min (abs (dq (k) + dq (k + 1)), 0.5 * q (k))
            if (dq (k) * dq (k + 1) <= 0.) then
                if (dq (k) > 0.) then ! local max
                    dm (k) = min (dm (k), dq (k), - dq (k + 1))
                else
                    dm (k) = 0.
                endif
            endif
        enddo
        dm (km) = 0.

>>>>>>> user/lnz/shield2022
        ! -----------------------------------------------------------------------
        ! cloud ice to snow autoconversion
        ! -----------------------------------------------------------------------
<<<<<<< HEAD
        
        call psaut_simp (ks, ke, dtm, qv, ql, qr, qi, qs, qg, tz, den)
        
    endif
    
end subroutine mp_fast
=======

        do k = 1, km
            dm (k) = max (dm (k), qvmin, h_var * q (k))
        enddo
    else
        do k = 1, km
            dm (k) = max (qvmin, h_var * q (k))
        enddo
    endif

end subroutine linear_prof
>>>>>>> user/lnz/shield2022

! =======================================================================
! sedimentation of cloud ice, snow, graupel or hail, and rain
! =======================================================================

<<<<<<< HEAD
subroutine sedimentation (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
        vtw, vtr, vti, vts, vtg, w1, r1, i1, s1, g1, pfw, pfr, pfi, pfs, pfg, &
        u, v, w, den, denfac, dte)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: dp, dz, den, denfac
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, u, v, w
    
    real, intent (out) :: w1, r1, i1, s1, g1
    
    real, intent (out), dimension (ks:ke) :: vtw, vtr, vti, vts, vtg, pfw, pfr, pfi, pfs, pfg
    
    real (kind = r8), intent (inout) :: dte
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real, dimension (ks:ke) :: q_liq, q_sol, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), dimension (ks:ke) :: te8, cvm
    
    w1 = 0.
    r1 = 0.
    i1 = 0.
    s1 = 0.
    g1 = 0.
    
    vtw = 0.
    vtr = 0.
    vti = 0.
    vts = 0.
    vtg = 0.
    
    pfw = 0.
    pfr = 0.
    pfi = 0.
    pfs = 0.
    pfg = 0.
=======
subroutine icloud (ks, ke, tzk, p1, qvk, qlk, qrk, qik, qsk, qgk, dp1, den, &
        ccn, cin, denfac, vts, vtg, vtr, qak, rh_adj, rh_rain, dts, h_var, &
        gsize, cond, dep, reevap, sub, last_step)

    implicit none

    logical, intent (in) :: last_step
    integer, intent (in) :: ks, ke
    real, intent (in), dimension (ks:ke) :: p1, dp1, den, denfac, vts, vtg, vtr, ccn
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tzk
    real, intent (inout), dimension (ks:ke) :: qvk, qlk, qrk, qik, qsk, qgk, qak
    real, intent (inout), dimension (ks:ke) :: cin
    real, intent (in) :: rh_adj, rh_rain, dts, h_var, gsize
    real, intent (out) :: cond, dep, reevap, sub
    ! local:
    real, dimension (ks:ke) :: icpk, di, qim
    real, dimension (ks:ke) :: q_liq, q_sol
    real (kind = r_grid), dimension (ks:ke) :: cvm, te8
    real (kind = r_grid) :: tz
    real :: rdts, fac_g2v, fac_v2g, fac_i2s, fac_imlt
    real :: qv, ql, qr, qi, qs, qg, melt
    real :: pracs, psacw, pgacw, psacr, pgacr, pgaci, praci, psaci
    real :: pgmlt, psmlt, pgfr, psaut
    real :: tc, dqs0, qden, qsm
    real :: dt5, factor, sink, qi_crt
    real :: tmp, qsw, qsi, dqsdt, dq
    real :: dtmp, qc, q_plus, q_minus
    integer :: k

    dt5 = 0.5 * dts
    rdts = 1. / dts
>>>>>>> user/lnz/shield2022

    ! -----------------------------------------------------------------------
    ! calculate heat capacities and latent heat coefficients
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    call cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, cvm, te8, tz, &
        lcpk, icpk, tcpk, tcp3)
    
=======

    fac_i2s = 1. - exp (- dts / tau_i2s)
    fac_g2v = 1. - exp (- dts / tau_g2v)
    fac_v2g = 1. - exp (- dts / tau_v2g)
    fac_imlt = 1. - exp (- dt5 / tau_imlt)

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! terminal fall and melting of falling cloud ice into rain
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    if (do_psd_ice_fall) then
        call term_rsg (ks, ke, qi, den, denfac, vi_fac, vconi, blini, normi, expoi, mui, vi_max, const_vi, vti)
    else
        call term_ice (ks, ke, tz, qi, den, vi_fac, vi_max, const_vi, vti)
    endif
    
    if (do_sedi_melt) then
        call sedi_melt (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
            vti, r1, tau_imlt, cvm, icpk, "qi")
    endif
    
    call terminal_fall (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
        vti, i1, pfi, u, v, w, dte, "qi")
    
    pfi (ks) = max (0.0, pfi (ks))
    do k = ke, ks + 1, -1
        pfi (k) = max (0.0, pfi (k) - pfi (k - 1))
    enddo

    ! -----------------------------------------------------------------------
    ! terminal fall and melting of falling snow into rain
    ! -----------------------------------------------------------------------
    
    call term_rsg (ks, ke, qs, den, denfac, vs_fac, vcons, blins, norms, expos, mus, vs_max, const_vs, vts)
    
    if (do_sedi_melt) then
        call sedi_melt (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
            vts, r1, tau_smlt, cvm, icpk, "qs")
    endif
    
    call terminal_fall (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
        vts, s1, pfs, u, v, w, dte, "qs")
    
    pfs (ks) = max (0.0, pfs (ks))
    do k = ke, ks + 1, -1
        pfs (k) = max (0.0, pfs (k) - pfs (k - 1))
    enddo

=======

    do k = ks, ke
        q_liq (k) = qlk (k) + qrk (k)
        q_sol (k) = qik (k) + qsk (k) + qgk (k)
        cvm (k) = one_r8 + qvk (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
        te8 (k) = cvm (k) * tzk (k) + lv00 * qvk (k) - li00 * q_sol (k)
        icpk (k) = (li00 + d1_ice * tzk (k)) / cvm (k)
    enddo

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! terminal fall and melting of falling graupel into rain
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    if (do_hail) then
        call term_rsg (ks, ke, qg, den, denfac, vg_fac, vconh, blinh, normh, expoh, muh, vg_max, const_vg, vtg)
    else
        call term_rsg (ks, ke, qg, den, denfac, vg_fac, vcong, bling, normg, expog, mug, vg_max, const_vg, vtg)
    endif
    
    if (do_sedi_melt) then
        call sedi_melt (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
            vtg, r1, tau_gmlt, cvm, icpk, "qg")
    endif
    
    call terminal_fall (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
        vtg, g1, pfg, u, v, w, dte, "qg")
    
    pfg (ks) = max (0.0, pfg (ks))
    do k = ke, ks + 1, -1
        pfg (k) = max (0.0, pfg (k) - pfg (k - 1))
    enddo
    
    ! -----------------------------------------------------------------------
    ! terminal fall of cloud water
    ! -----------------------------------------------------------------------
    
    if (do_psd_water_fall) then

        call term_rsg (ks, ke, ql, den, denfac, vw_fac, vconw, blinw, normw, expow, muw, vw_max, const_vw, vtw)
    
        call terminal_fall (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
            vtw, w1, pfw, u, v, w, dte, "ql")

        pfw (ks) = max (0.0, pfw (ks))
        do k = ke, ks + 1, -1
            pfw (k) = max (0.0, pfw (k) - pfw (k - 1))
        enddo

    endif
    
    ! -----------------------------------------------------------------------
    ! terminal fall of rain
    ! -----------------------------------------------------------------------
    
    call term_rsg (ks, ke, qr, den, denfac, vr_fac, vconr, blinr, normr, expor, mur, vr_max, const_vr, vtr)
    
    call terminal_fall (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
        vtr, r1, pfr, u, v, w, dte, "qr")
    
    pfr (ks) = max (0.0, pfr (ks))
    do k = ke, ks + 1, -1
        pfr (k) = max (0.0, pfr (k) - pfr (k - 1))
    enddo

end subroutine sedimentation
=======

    do k = ks, ke
        if (qi0_crt < 0.) then
            qim (k) = - qi0_crt
        else
            qim (k) = qi0_crt / den (k)
        endif
    enddo

    if (.not. do_warm_rain_mp) then

        ! -----------------------------------------------------------------------
        ! sources of cloud ice: pihom, cold rain, and the sat_adj
        ! (initiation plus deposition)
        ! sources of snow: cold rain, auto conversion + accretion (from cloud ice)
        ! sat_adj (deposition; requires pre - existing snow) ; initial snow comes from auto conversion
        ! -----------------------------------------------------------------------

        do k = ks, ke
            if (tzk (k) > tice_mlt .and. qik (k) > qcmin) then

                ! -----------------------------------------------------------------------
                ! pimlt: instant melting of cloud ice
                ! -----------------------------------------------------------------------

                melt = min (qik (k), fac_imlt * (tzk (k) - tice_mlt) / icpk (k))
                tmp = min (melt, dim (ql_mlt, qlk (k))) ! max ql amount
                qlk (k) = qlk (k) + tmp
                qrk (k) = qrk (k) + melt - tmp
                qik (k) = qik (k) - melt
                q_liq (k) = q_liq (k) + melt
                q_sol (k) = q_sol (k) - melt
            elseif (tzk (k) < t_wfr .and. qlk (k) > qcmin) then

                ! -----------------------------------------------------------------------
                ! pihom: homogeneous freezing of cloud water into cloud ice
                ! -----------------------------------------------------------------------

                dtmp = t_wfr - tzk (k)
                factor = min (1., dtmp / dt_fr)
                sink = min (qlk (k) * factor, dtmp / icpk (k))
                tmp = min (sink, dim (qim (k), qik (k)))
                qlk (k) = qlk (k) - sink
                qsk (k) = qsk (k) + sink - tmp
                qik (k) = qik (k) + tmp
                q_liq (k) = q_liq (k) - sink
                q_sol (k) = q_sol (k) + sink
            endif
        enddo

        ! -----------------------------------------------------------------------
        ! vertical subgrid variability
        ! -----------------------------------------------------------------------

        call linear_prof (ke - ks + 1, qik (ks), di (ks), z_slope_ice, h_var)

        ! -----------------------------------------------------------------------
        ! update capacity heat and latend heat coefficient
        ! -----------------------------------------------------------------------

        do k = ks, ke
            cvm (k) = one_r8 + qvk (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
            tzk (k) = (te8 (k) - lv00 * qvk (k) + li00 * q_sol (k)) / cvm (k)
            icpk (k) = (li00 + d1_ice * tzk (k)) / cvm (k)
        enddo

        do k = ks, ke

            ! -----------------------------------------------------------------------
            ! do nothing above p_min
            ! -----------------------------------------------------------------------

            if (p1 (k) < p_min) cycle

            tz = tzk (k)
            qv = qvk (k)
            ql = qlk (k)
            qi = qik (k)
            qr = qrk (k)
            qs = qsk (k)
            qg = qgk (k)

            pgacr = 0.
            pgacw = 0.
            tc = tz - tice

            if (tc .ge. 0.) then

                ! -----------------------------------------------------------------------
                ! melting of snow
                ! -----------------------------------------------------------------------

                dqs0 = ces0 / p1 (k) - qv ! not sure if this is correct; check again

                if (qs > qcmin) then

                    ! -----------------------------------------------------------------------
                    ! psacw: accretion of cloud water by snow
                    ! only rate is used (for snow melt) since tc > 0.
                    ! -----------------------------------------------------------------------

                    if (ql > qrmin) then
                        factor = denfac (k) * csacw * exp (0.8125 * log (qs * den (k)))
                        psacw = factor / (1. + dts * factor) * ql ! rate
                    else
                        psacw = 0.
                    endif

                    ! -----------------------------------------------------------------------
                    ! psacr: accretion of rain by melted snow
                    ! pracs: accretion of snow by rain
                    ! -----------------------------------------------------------------------

                    if (qr > qrmin) then
                        psacr = min (acr3d (vts (k), vtr (k), qr, qs, csacr, acco (1, 2), &
                            den (k)), qr * rdts)
                        pracs = acr3d (vtr (k), vts (k), qs, qr, cracs, acco (1, 1), den (k))
                    else
                        psacr = 0.
                        pracs = 0.
                    endif

                    ! -----------------------------------------------------------------------
                    ! total snow sink:
                    ! psmlt: snow melt (due to rain accretion)
                    ! -----------------------------------------------------------------------

                    psmlt = max (0., smlt (tc, dqs0, qs * den (k), psacw, psacr, csmlt, &
                        den (k), denfac (k)))
                    sink = min (qs, dts * (psmlt + pracs), tc / icpk (k))
                    qs = qs - sink
                    tmp = min (sink, dim (qs_mlt, ql)) ! max ql due to snow melt
                    ql = ql + tmp
                    qr = qr + sink - tmp
                    q_liq (k) = q_liq (k) + sink
                    q_sol (k) = q_sol (k) - sink

                    cvm (k) = one_r8 + qv * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                    tz = (te8 (k) - lv00 * qv + li00 * q_sol (k)) / cvm (k)
                    tc = tz - tice
                    icpk (k) = (li00 + d1_ice * tz) / cvm (k)

                endif

                ! -----------------------------------------------------------------------
                ! melting of graupel
                ! -----------------------------------------------------------------------

                if (qg > qcmin .and. tc > 0.) then

                    ! -----------------------------------------------------------------------
                    ! pgacr: accretion of rain by graupel
                    ! -----------------------------------------------------------------------

                    if (qr > qrmin) &
                        pgacr = min (acr3d (vtg (k), vtr (k), qr, qg, cgacr, acco (1, 3), &
                        den (k)), rdts * qr)

                    ! -----------------------------------------------------------------------
                    ! pgacw: accretion of cloud water by graupel
                    ! -----------------------------------------------------------------------

                    qden = qg * den (k)
                    if (ql > qrmin) then
                        factor = cgacw * qden / sqrt (den (k) * sqrt (sqrt (qden)))
                        pgacw = factor / (1. + dts * factor) * ql ! rate
                    endif

                    ! -----------------------------------------------------------------------
                    ! pgmlt: graupel melt
                    ! -----------------------------------------------------------------------

                    pgmlt = dts * gmlt (tc, dqs0, qden, pgacw, pgacr, cgmlt, den (k))
                    pgmlt = min (max (0., pgmlt), qg, tc / icpk (k))
                    qg = qg - pgmlt
                    qr = qr + pgmlt
                    q_liq (k) = q_liq (k) + pgmlt
                    q_sol (k) = q_sol (k) - pgmlt
                    cvm (k) = one_r8 + qv * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                    tz = (te8 (k) - lv00 * qv + li00 * q_sol (k)) / cvm (k)
                endif

            else

                ! -----------------------------------------------------------------------
                ! cloud ice proc:
                ! -----------------------------------------------------------------------

                ! -----------------------------------------------------------------------
                ! psaci: accretion of cloud ice by snow
                ! -----------------------------------------------------------------------

                if (qi > 3.e-7) then ! cloud ice sink terms

                    if (qs > 1.e-7) then
                        ! -----------------------------------------------------------------------
                        ! sjl added (following lin eq. 23) the temperature dependency
                        ! to reduce accretion, use esi = exp (0.05 * tc) as in hong et al 2004
                        ! -----------------------------------------------------------------------
                        factor = dts * denfac (k) * csaci * exp (0.05 * tc + 0.8125 * log (qs * den (k)))
                        psaci = factor / (1. + factor) * qi
                    else
                        psaci = 0.
                    endif

                    ! -----------------------------------------------------------------------
                    ! assuming linear subgrid vertical distribution of cloud ice
                    ! the mismatch computation following lin et al. 1994, mwr
                    ! -----------------------------------------------------------------------

                    if (const_vi) then
                        tmp = fac_i2s
                    else
                        tmp = fac_i2s * exp (0.025 * tc)
                    endif

                    di (k) = max (di (k), qrmin)
                    q_plus = qi + di (k)
                    if (q_plus > (qim (k) + qrmin)) then
                        if (qim (k) > (qi - di (k))) then
                            dq = (0.25 * (q_plus - qim (k)) ** 2) / di (k)
                        else
                            dq = qi - qim (k)
                        endif
                        psaut = tmp * dq
                    else
                        psaut = 0.
                    endif
                    ! -----------------------------------------------------------------------
                    ! sink is no greater than 75% of qi
                    ! -----------------------------------------------------------------------
                    sink = min (0.75 * qi, psaci + psaut)
                    qi = qi - sink
                    qs = qs + sink

                    ! -----------------------------------------------------------------------
                    ! pgaci: accretion of cloud ice by graupel
                    ! -----------------------------------------------------------------------

                    if (qg > 1.e-6) then
                        ! -----------------------------------------------------------------------
                        ! factor = dts * cgaci / sqrt (den (k)) * exp (0.05 * tc + 0.875 * log (qg * den (k)))
                        ! simplified form: remove temp dependency & set the exponent "0.875" -- > 1
                        ! -----------------------------------------------------------------------
                        factor = dts * cgaci / sqrt (den (k)) * exp (0.875 * log (qg * den (k)))
                        pgaci = factor / (1. + factor) * qi
                        qi = qi - pgaci
                        qg = qg + pgaci
                    endif

                endif

                ! -----------------------------------------------------------------------
                ! cold - rain proc:
                ! -----------------------------------------------------------------------

                ! -----------------------------------------------------------------------
                ! rain to ice, snow, graupel processes:
                ! -----------------------------------------------------------------------

                tc = tz - tice

                if (qr > 1.e-7 .and. tc < 0.) then

                    ! -----------------------------------------------------------------------
                    ! * sink * terms to qr: psacr + pgfr
                    ! source terms to qs: psacr
                    ! source terms to qg: pgfr
                    ! -----------------------------------------------------------------------

                    ! -----------------------------------------------------------------------
                    ! psacr accretion of rain by snow
                    ! -----------------------------------------------------------------------

                    if (qs > 1.e-7) then ! if snow exists
                        psacr = dts * acr3d (vts (k), vtr (k), qr, qs, csacr, acco (1, 2), den (k))
                    else
                        psacr = 0.
                    endif

                    ! -----------------------------------------------------------------------
                    ! pgfr: rain freezing -- > graupel
                    ! -----------------------------------------------------------------------

                    pgfr = dts * cgfr (1) / den (k) * (exp (- cgfr (2) * tc) - 1.) * &
                        exp (1.75 * log (qr * den (k)))

                    ! -----------------------------------------------------------------------
                    ! total sink to qr
                    ! -----------------------------------------------------------------------

                    sink = psacr + pgfr
                    factor = min (sink, qr, - tc / icpk (k)) / max (sink, qrmin)

                    psacr = factor * psacr
                    pgfr = factor * pgfr

                    sink = psacr + pgfr
                    qr = qr - sink
                    qs = qs + psacr
                    qg = qg + pgfr
                    q_liq (k) = q_liq (k) - sink
                    q_sol (k) = q_sol (k) + sink

                    cvm (k) = one_r8 + qv * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
                    tz = (te8 (k) - lv00 * qv + li00 * q_sol (k)) / cvm (k)
                    icpk (k) = (li00 + d1_ice * tz) / cvm (k)
                endif

                ! -----------------------------------------------------------------------
                ! graupel production terms:
                ! -----------------------------------------------------------------------

                if (qs > 1.e-7) then

                    ! -----------------------------------------------------------------------
                    ! accretion: snow -- > graupel
                    ! -----------------------------------------------------------------------

                    if (qg > qrmin) then
                        sink = dts * acr3d (vtg (k), vts (k), qs, qg, cgacs, acco (1, 4), den (k))
                    else
                        sink = 0.
                    endif

                    ! -----------------------------------------------------------------------
                    ! autoconversion snow -- > graupel
                    ! -----------------------------------------------------------------------

                    qsm = qs0_crt / den (k)
                    if (qs > qsm) then
                        factor = dts * 1.e-3 * exp (0.09 * (tz - tice))
                        sink = sink + factor / (1. + factor) * (qs - qsm)
                    endif
                    sink = min (qs, sink)
                    qs = qs - sink
                    qg = qg + sink

                endif ! snow existed

                if (qg > 1.e-7 .and. tz < tice) then

                    ! -----------------------------------------------------------------------
                    ! pgacw: accretion of cloud water by graupel
                    ! -----------------------------------------------------------------------

                    if (ql > 1.e-6) then
                        qden = qg * den (k)
                        factor = dts * cgacw * qden / sqrt (den (k) * sqrt (sqrt (qden)))
                        pgacw = factor / (1. + factor) * ql
                    else
                        pgacw = 0.
                    endif

                    ! -----------------------------------------------------------------------
                    ! pgacr: accretion of rain by graupel
                    ! -----------------------------------------------------------------------

                    if (qr > 1.e-6) then
                        pgacr = min (dts * acr3d (vtg (k), vtr (k), qr, qg, cgacr, acco (1, 3), &
                            den (k)), qr)
                    else
                        pgacr = 0.
                    endif

                    sink = pgacr + pgacw
                    factor = min (sink, dim (tice, tz) / icpk (k)) / max (sink, qrmin)
                    pgacr = factor * pgacr
                    pgacw = factor * pgacw

                    sink = pgacr + pgacw
                    qg = qg + sink
                    qr = qr - pgacr
                    ql = ql - pgacw

                    q_liq (k) = q_liq (k) - sink
                    q_sol (k) = q_sol (k) + sink
                    tz = (te8 (k) - lv00 * qv + li00 * q_sol (k)) / (one_r8 + qv * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice)
                endif

            endif

            tzk (k) = tz
            qvk (k) = qv
            qlk (k) = ql
            qik (k) = qi
            qrk (k) = qr
            qsk (k) = qs
            qgk (k) = qg

        enddo

    endif

    call subgrid_z_proc (ks, ke, p1, den, denfac, dts, rh_adj, tzk, qvk, qlk, &
        qrk, qik, qsk, qgk, qak, dp1, h_var, rh_rain, te8, ccn, cin, gsize, &
        cond, dep, reevap, sub, last_step)

end subroutine icloud
>>>>>>> user/lnz/shield2022

! =======================================================================
! terminal velocity for cloud ice
! =======================================================================

<<<<<<< HEAD
subroutine term_ice (ks, ke, tz, q, den, v_fac, v_max, const_v, vt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    logical, intent (in) :: const_v
    
    real, intent (in) :: v_fac, v_max
    
    real, intent (in), dimension (ks:ke) :: q, den
    
    real (kind = r8), intent (in), dimension (ks:ke) :: tz
    
    real, intent (out), dimension (ks:ke) :: vt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: qden
    
    real, parameter :: aa = - 4.14122e-5
    real, parameter :: bb = - 0.00538922
    real, parameter :: cc = - 0.0516344
    real, parameter :: dd = 0.00216078
    real, parameter :: ee = 1.9714
    
    real, dimension (ks:ke) :: tc
    
    if (const_v) then
        vt (:) = v_fac
=======
subroutine subgrid_z_proc (ks, ke, p1, den, denfac, dts, rh_adj, tz, qv, ql, qr, &
        qi, qs, qg, qa, dp1, h_var, rh_rain, te8, ccn, cin, gsize, cond, dep, reevap, sub, last_step)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: dts, rh_adj, h_var, rh_rain, gsize
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
    real :: pidep, qi_crt
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
    real :: fac_l2v, fac_v2l, fac_g2v, fac_v2g
    integer :: k

    if (do_sat_adj) then
        dt_evap = 0.5 * dts
>>>>>>> user/lnz/shield2022
    else
        do k = ks, ke
            qden = q (k) * den (k)
            if (q (k) .lt. qfmin) then
                vt (k) = 0.0
            else
                tc (k) = tz (k) - tice
                if (ifflag .eq. 1) then
                    vt (k) = (3. + log10 (qden)) * (tc (k) * (aa * tc (k) + bb) + cc) + &
                        dd * tc (k) + ee
                    vt (k) = 0.01 * v_fac * exp (vt (k) * log (10.))
                endif
                if (ifflag .eq. 2) &
                    vt (k) = v_fac * 3.29 * exp (0.16 * log (qden))
                vt (k) = min (v_max, max (0.0, vt (k)))
            endif
        enddo
    endif
<<<<<<< HEAD
    
end subroutine term_ice

! =======================================================================
! terminal velocity for rain, snow, and graupel, Lin et al. (1983)
! =======================================================================

subroutine term_rsg (ks, ke, q, den, denfac, v_fac, vcon, blin, norm, expo, mu, v_max, const_v, vt)
    
    implicit none
    
=======

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    integer, intent (in) :: ks, ke
    
    logical, intent (in) :: const_v
    
    real, intent (in) :: v_fac, blin, v_max, mu
    
    real (kind = r8), intent (in) :: vcon, norm, expo
    
    real, intent (in), dimension (ks:ke) :: q, den, denfac
    
    real, intent (out), dimension (ks:ke) :: vt
    
=======

    fac_l2v = 1. - exp (- dt_evap / tau_l2v)
    fac_v2l = 1. - exp (- dt_evap / tau_v2l)
    fac_g2v = 1. - exp (- dts / tau_g2v)
    fac_v2g = 1. - exp (- dts / tau_v2g)

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    integer :: k
    
    real :: qden
    
    if (const_v) then
        vt (:) = v_fac
    else
        do k = ks, ke
            qden = q (k) * den (k)
            if (q (k) .lt. qfmin) then
                vt (k) = 0.0
            else
                vt (k) = v_fac * vcon * denfac (k) * exp (1. / (mu + 3) * blin * log (6 * qden / norm)) * &
                    exp (- blin * log (expo))
                vt (k) = min (v_max, max (0.0, vt (k)))
=======

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

        if (p1 (k) < p_min) cycle

        if (.not. do_warm_rain_mp) then

            ! -----------------------------------------------------------------------
            ! instant deposit all water vapor to cloud ice when temperature is super low
            ! -----------------------------------------------------------------------

            if (tz (k) < t_min) then
                sink = dim (qv (k), 1.e-7)
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
            if (tin > t_sub + 6.) then
                rh = qpz / iqs1 (tin, den (k))
                if (rh < rh_adj) then ! qpz / rh_adj < qs
                    reevap = reevap + ql (k) * dp1 (k)
                    sub = sub + qi (k) * dp1 (k)
                    tz (k) = tin
                    qv (k) = qpz
                    ql (k) = 0.
                    qi (k) = 0.
                    cycle ! cloud free
                endif
>>>>>>> user/lnz/shield2022
            endif
        enddo
    endif
    
end subroutine term_rsg

<<<<<<< HEAD
! =======================================================================
! melting during sedimentation
! =======================================================================

subroutine sedi_melt (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
        vt, r1, tau_mlt, cvm, icpk, qflag)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts, tau_mlt
    
    real, intent (in), dimension (ks:ke) :: vt, dp, dz, icpk
    
    real (kind = r8), intent (in), dimension (ks:ke) :: cvm
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real, intent (inout) :: r1
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    character (len = 2), intent (in) :: qflag
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k, m
    
    real :: dtime, sink, zs
    
    real, dimension (ks:ke) :: q
    
    real, dimension (ks:ke + 1) :: ze, zt
    
    call zezt (ks, ke, dts, zs, dz, vt, ze, zt)
    
    select case (qflag)
        case ("qi")
            q = qi
        case ("qs")
            q = qs
        case ("qg")
            q = qg
        case default
            print *, "gfdl_mp: qflag error!"
    end select
    
    ! -----------------------------------------------------------------------
    ! melting to rain
    ! -----------------------------------------------------------------------
    
    do k = ke - 1, ks, - 1
        if (vt (k) .lt. 1.e-10) cycle
        if (q (k) .gt. qcmin) then
            do m = k + 1, ke
                if (zt (k + 1) .ge. ze (m)) exit
                if (zt (k) .lt. ze (m + 1) .and. tz (m) .gt. tice) then
                    dtime = min (dts, (ze (m) - ze (m + 1)) / vt (k))
                    dtime = min (1.0, dtime / tau_mlt)
                    sink = min (q (k) * dp (k) / dp (m), dtime * (tz (m) - tice) / icpk (m))
                    q (k) = q (k) - sink * dp (m) / dp (k)
                    if (zt (k) .lt. zs) then
                        r1 = r1 + sink * dp (m)
                    else
                        qr (m) = qr (m) + sink
                    endif
                    select case (qflag)
                        case ("qi")
                            qi (k) = q (k)
                        case ("qs")
                            qs (k) = q (k)
                        case ("qg")
                            qg (k) = q (k)
                        case default
                            print *, "gfdl_mp: qflag error!"
                    end select
                    tz (m) = (tz (m) * cvm (m) - li00 * sink) / &
                        mhc (qv (k), ql (m), qr (m), qi (m), qs (m), qg (m))
=======
        endif

        ! -----------------------------------------------------------------------
        ! cloud water < -- > vapor adjustment:
        ! -----------------------------------------------------------------------

        tin = tz (k)
        rh_tem = qpz / iqs1 (tin, den (k))
        qsw = wqs2 (tin, den (k), dwsdt)
        dq0 = qsw - qv (k)
        if (use_rhc_cevap) then
            evap = 0.
            if (rh_tem .lt. rhc_cevap) then
                if (dq0 > 0.) then ! evaporation
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
>>>>>>> user/lnz/shield2022
                endif
                if (q (k) .lt. qcmin) exit
            enddo
        endif
<<<<<<< HEAD
    enddo
    
end subroutine sedi_melt
=======
        ! sjl on jan 23 2018: reversible evap / condensation:
        qv (k) = qv (k) + evap
        ql (k) = ql (k) - evap
        q_liq (k) = q_liq (k) - evap

        cvm (k) = one_r8 + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
        tz (k) = (te8 (k) - lv00 * qv (k) + li00 * q_sol (k)) / cvm (k)

        ! -----------------------------------------------------------------------
        ! update heat capacity and latend heat coefficient
        ! -----------------------------------------------------------------------

        icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)

        if (.not. do_warm_rain_mp) then

            ! -----------------------------------------------------------------------
            ! enforce complete freezing below - 48 c
            ! -----------------------------------------------------------------------

            dtmp = t_wfr - tz (k) ! [ - 40, - 48]
            if (dtmp > 0. .and. ql (k) > qcmin) then
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
                if (ql (k) > qrmin .and. tc > 0.1) then
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
            ! update capacity heat and latend heat coefficient
            ! -----------------------------------------------------------------------

            tcpk (k) = (li20 + (d1_vap + d1_ice) * tz (k)) / cvm (k)

            ! -----------------------------------------------------------------------
            ! sublimation / deposition of ice
            ! -----------------------------------------------------------------------

            if (tz (k) < tice) then
                qsi = iqs2 (tz (k), den (k), dqsdt)
                dq = qv (k) - qsi
                sink = dq / (1. + tcpk (k) * dqsdt)
                if (qi (k) > qrmin) then
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
                         / (qsi * den (k) * lat2 / (0.0243 * rvgas * tz (k) ** 2) + 4.42478e4)
                else
                    pidep = 0.
                endif
                if (dq > 0.) then ! vapor - > ice
                    tmp = tice - tz (k)
                    ! 20160912: the following should produce more ice at higher altitude
                    ! qi_crt = 4.92e-11 * exp (1.33 * log (1.e3 * exp (0.1 * tmp))) / den (k)
                    qi_crt = qi_gen * min (qi_lim, 0.1 * tmp) / den (k)
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
            ! update capacity heat and latend heat coefficient
            ! -----------------------------------------------------------------------

            tcpk (k) = (li20 + (d1_vap + d1_ice) * tz (k)) / cvm (k)

            ! -----------------------------------------------------------------------
            ! sublimation / deposition of snow
            ! this process happens for all temp rage
            ! -----------------------------------------------------------------------

            if (qs (k) > qrmin) then
                qsi = iqs2 (tz (k), den (k), dqsdt)
                qden = qs (k) * den (k)
                tmp = exp (0.65625 * log (qden))
                tsq = tz (k) * tz (k)
                dq = (qsi - qv (k)) / (1. + tcpk (k) * dqsdt)
                pssub = cssub (1) * tsq * (cssub (2) * sqrt (qden) + cssub (3) * tmp * &
                    sqrt (denfac (k))) / (cssub (4) * tsq + cssub (5) * qsi * den (k))
                pssub = (qsi - qv (k)) * dts * pssub
                if (pssub > 0.) then ! qs -- > qv, sublimation
                    pssub = min (pssub * min (1., dim (tz (k), t_sub) * 0.2), qs (k))
                    sub = sub + pssub * dp1 (k)
                else
                    if (tz (k) > tice) then
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

            if (qg (k) > qrmin) then
                qsi = iqs2 (tz (k), den (k), dqsdt)
                qden = qg (k) * den (k)
                tmp = exp (0.6875 * log (qden))
                tsq = tz (k) * tz (k)
                dq = (qsi - qv (k)) / (1. + tcpk (k) * dqsdt)
                pgsub = cgsub (1) * tsq * (cgsub (2) * sqrt (qden) + cgsub (3) * tmp / &
                    sqrt (sqrt (den (k)))) / (cgsub (4) * tsq + cgsub (5) * qsi * den (k))
                pgsub = (qsi - qv (k)) * dts * pgsub
                if (pgsub > 0.) then ! qs -- > qv, sublimation
                    pgsub = min (pgsub * min (1., dim (tz (k), t_sub) * 0.2), qg (k))
                    sub = sub + pgsub * dp1 (k)
                else
                    if (tz (k) > tice) then
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

        if (tin <= t_wfr) then
            ! ice phase:
            qstar = iqs1 (tin, den (k))
        elseif (tin >= tice) then
            ! liquid phase:
            qstar = wqs1 (tin, den (k))
        else
            ! mixed phase:
            qsi = iqs1 (tin, den (k))
            qsw = wqs1 (tin, den (k))
            if (q_cond (k) > 3.e-6) then
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

        qpz = cld_fac * qpz
        rh = qpz / qstar

        ! -----------------------------------------------------------------------
        ! icloud_f = 0: bug - fixed
        ! icloud_f = 1: old fvgfs gfdl) mp implementation
        ! icloud_f = 2: binary cloud scheme (0 / 1)
        ! icloud_f = 3: revision of icloud = 0
        ! -----------------------------------------------------------------------

        if (use_xr_cloud) then ! xu and randall cloud scheme (1996)
            if (rh >= 1.0) then
                qa (k) = 1.0
            elseif (rh > rh_thres .and. q_cond (k) > 1.e-6) then
                qa (k) = rh ** xr_a * (1.0 - exp (- xr_b * max (0.0, q_cond (k)) / &
                    max (1.e-5, (max (1.e-10, 1.0 - rh) * qstar) ** xr_c)))
                qa (k) = max (0.0, min (1., qa (k)))
            else
                qa (k) = 0.0
            endif
        elseif (use_park_cloud) then ! park et al. 2016 (mon. wea. review)
            if (q_cond (k) > 1.e-6) then
                qa (k) = 1. / 50. * (5.77 * (100. - gsize / 1000.) * max (0.0, q_cond (k) * 1000.) ** 1.07 + &
                    4.82 * (gsize / 1000. - 50.) * max (0.0, q_cond (k) * 1000.) ** 0.94)
                qa (k) = qa (k) * (0.92 / 0.96 * q_liq (k) / q_cond (k) + 1.0 / 0.96 * q_sol (k) / q_cond (k))
                qa (k) = max (0.0, min (1., qa (k)))
            else
                qa (k) = 0.0
            endif
        elseif (use_gi_cloud) then ! gultepe and isaac (2007)
            sigma = 0.28 + max (0.0, q_cond (k) * 1000.) ** 0.49
            gam = max (0.0, q_cond (k) * 1000.) / sigma
            if (gam < 0.18) then
                qa10 = 0.
            elseif (gam > 2.0) then
                qa10 = 1.0
            else
                qa10 = - 0.1754 + 0.9811 * gam - 0.2223 * gam ** 2 + 0.0104 * gam ** 3
                qa10 = max (0.0, min (1., qa10))
            endif
            if (gam < 0.12) then
                qa100 = 0.
            elseif (gam > 1.85) then
                qa100 = 1.0
            else
                qa100 = - 0.0913 + 0.7213 * gam + 0.1060 * gam ** 2 - 0.0946 * gam ** 3
                qa100 = max (0.0, min (1., qa100))
            endif
            qa (k) = qa10 + (log10 (gsize / 1000.) - 1) * (qa100 - qa10)
            qa (k) = max (0.0, min (1., qa (k)))
        else
            if (rh > rh_thres .and. qpz > 1.e-6) then

                dq = h_var * qpz
                if (do_cld_adj) then
                    q_plus = qpz + dq * f_dq_p * min(1.0, max(0.0, (p1 (k) - 200.e2) / (1000.e2 - 200.e2)))
                else
                    q_plus = qpz + dq * f_dq_p
                endif
                q_minus = qpz - dq * f_dq_m

                if (icloud_f .eq. 2) then
                    if (qstar < qpz) then
                        qa (k) = 1.
                    else
                        qa (k) = 0.
                    endif
                elseif (icloud_f .eq. 3) then
                    if (qstar < qpz) then
                        qa (k) = 1.
                    else
                        if (qstar < q_plus) then
                            qa (k) = (q_plus - qstar) / (dq * f_dq_p)
                        else
                            qa (k) = 0.
                        endif
                        ! impose minimum cloudiness if substantial q_cond (k) exist
                        if (q_cond (k) > 1.e-6) then
                            qa (k) = max (cld_min, qa (k))
                        endif
                        qa (k) = min (1., qa (k))
                    endif
                else
                    if (qstar < q_minus) then
                        qa (k) = 1.
                    else
                        if (qstar < q_plus) then
                            if (icloud_f .eq. 0) then
                                qa (k) = (q_plus - qstar) / (dq * f_dq_p + dq * f_dq_m)
                            else
                                qa (k) = (q_plus - qstar) / ((dq * f_dq_p + dq * f_dq_m) * (1. - q_cond (k)))
                            endif
                        else
                            qa (k) = 0.
                        endif
                        ! impose minimum cloudiness if substantial q_cond (k) exist
                        if (q_cond (k) > 1.e-6) then
                            qa (k) = max (cld_min, qa (k))
                        endif
                        qa (k) = min (1., qa (k))
                    endif
                endif
            else
                qa (k) = 0.
            endif
        endif

    enddo

end subroutine subgrid_z_proc
>>>>>>> user/lnz/shield2022

! =======================================================================
! melting during sedimentation
! =======================================================================

<<<<<<< HEAD
subroutine terminal_fall (dts, ks, ke, tz, qv, ql, qr, qi, qs, qg, dz, dp, &
        vt, x1, m1, u, v, w, dte, qflag)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: vt, dp, dz
    
    character (len = 2), intent (in) :: qflag
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, u, v, w
    
    real, intent (inout) :: x1
    
    real (kind = r8), intent (inout) :: dte
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    real, intent (out), dimension (ks:ke) :: m1
    
=======
subroutine revap_rac1 (hydrostatic, is, ie, dt, tz, qv, ql, qr, qi, qs, qg, den, hvar)

    implicit none

    logical, intent (in) :: hydrostatic

    integer, intent (in) :: is, ie

    real, intent (in) :: dt ! time step (s)

    real, intent (in), dimension (is:ie) :: den, hvar, qi, qs, qg

    real, intent (inout), dimension (is:ie) :: tz, qv, qr, ql

    real, dimension (is:ie) :: lcp2, denfac, q_liq, q_sol, cvm, lhl

    real :: dqv, qsat, dqsdt, evap, qden, q_plus, q_minus, sink
    real :: tin, t2, qpz, dq, dqh

    integer :: i

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    integer :: k
    
    logical :: no_fall
    
    real :: zs
    
    real, dimension (ks:ke) :: dm, q
    
    real, dimension (ks:ke + 1) :: ze, zt
    
    real (kind = r8), dimension (ks:ke) :: te1, te2

    m1 = 0.0
    
    call zezt (ks, ke, dts, zs, dz, vt, ze, zt)
    
    select case (qflag)
        case ("ql")
            q = ql
        case ("qr")
            q = qr
        case ("qi")
            q = qi
        case ("qs")
            q = qs
        case ("qg")
            q = qg
        case default
            print *, "gfdl_mp: qflag error!"
    end select
    
    call check_column (ks, ke, q, no_fall)
    
    if (no_fall) return
    
=======

    do i = is, ie
        lhl (i) = lv00 + d0_vap * tz (i)
        q_liq (i) = ql (i) + qr (i)
        q_sol (i) = qi (i) + qs (i) + qg (i)
        cvm (i) = c_air + qv (i) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
        lcp2 (i) = lhl (i) / cvm (i)
        ! denfac (i) = sqrt (sfcrho / den (i))
    enddo

    do i = is, ie
        if (qr (i) > qrmin .and. tz (i) > t_wfr) then
            qpz = qv (i) + ql (i)
            tin = tz (i) - lcp2 (i) * ql (i) ! presence of clouds suppresses the rain evap
            qsat = wqs2 (tin, den (i), dqsdt)
            dqh = max (ql (i), hvar (i) * max (qpz, qcmin))
            dqv = qsat - qv (i)
            q_minus = qpz - dqh
            q_plus = qpz + dqh

            ! -----------------------------------------------------------------------
            ! qsat must be > q_minus to activate evaporation
            ! qsat must be < q_plus to activate accretion
            ! -----------------------------------------------------------------------

            ! -----------------------------------------------------------------------
            ! rain evaporation
            ! -----------------------------------------------------------------------

            if (dqv > qvmin .and. qsat > q_minus) then
                if (qsat > q_plus) then
                    dq = qsat - qpz
                else
                    ! q_minus < qsat < q_plus
                    ! dq == dqh if qsat == q_minus
                    dq = 0.25 * (q_minus - qsat) ** 2 / dqh
                endif
                qden = qr (i) * den (i)
                t2 = tin * tin
                evap = crevp (1) * t2 * dq * (crevp (2) * sqrt (qden) + crevp (3) * exp (0.725 * log (qden))) &
                     / (crevp (4) * t2 + crevp (5) * qsat * den (i))
                evap = min (qr (i), dt * evap, dqv / (1. + lcp2 (i) * dqsdt))
                qr (i) = qr (i) - evap
                qv (i) = qv (i) + evap
                q_liq (i) = q_liq (i) - evap
                cvm (i) = c_air + qv (i) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                tz (i) = tz (i) - evap * lhl (i) / cvm (i)
            endif

            ! -----------------------------------------------------------------------
            ! accretion: pracc
            ! -----------------------------------------------------------------------

            if (qr (i) > qrmin .and. ql (i) > 1.e-8 .and. qsat < q_plus) then
                denfac (i) = sqrt (sfcrho / den (i))
                sink = dt * denfac (i) * cracw * exp (0.95 * log (qr (i) * den (i)))
                sink = sink / (1. + sink) * ql (i)
                ql (i) = ql (i) - sink
                qr (i) = qr (i) + sink
            endif
        endif
    enddo

end subroutine revap_rac1

! =======================================================================
! compute terminal fall speed
! consider cloud ice, snow, and graupel's melting during fall
! =======================================================================

subroutine terminal_fall (dtm, ks, ke, tz, qv, ql, qr, qg, qs, qi, dz, dp, &
        den, vtg, vts, vti, r1, g1, s1, i1, pfr, pfi, pfs, pfg, m1_sol, w1, dte)
    
    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in) :: dtm ! time step (s)
    real, intent (in), dimension (ks:ke) :: vtg, vts, vti, den, dp, dz
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: tz
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qg, qs, qi, m1_sol, w1, pfr, pfi, pfs, pfg
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
    fac_imlt = 1. - exp (- dt5 / tau_imlt)
    pfr = 0.0
    pfi = 0.0
    pfs = 0.0
    pfg = 0.0
    
    ! -----------------------------------------------------------------------
    ! define heat capacity and latend heat coefficient
    ! -----------------------------------------------------------------------

    do k = ks, ke
        m1_sol (k) = 0.
        q_liq (k) = ql (k) + qr (k)
        q_sol (k) = qi (k) + qs (k) + qg (k)
        cvm (k) = 1. + qv (k) * c1_vap + q_liq (k) * c1_liq + q_sol (k) * c1_ice
        lcpk (k) = (lv00 + d1_vap * tz (k)) / cvm (k)
        icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)
    enddo

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! momentum transportation during sedimentation
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    if (do_sedi_w) then
        do k = ks, ke
            dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
        enddo
    endif
    
=======

    k0 = ke
    do k = ks, ke - 1
        if (tz (k) > tice) then
            k0 = k
            exit
        endif
    enddo

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! energy change during sedimentation
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    if (consv_checker) then
        do k = ks, ke
            te1 (k) = mte (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), tz (k), dp (k), .false.)
        enddo
    endif
    
=======

    do k = k0, ke
        tc = tz (k) - tice
        if (qi (k) > qcmin .and. tc > 0.) then
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

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! sedimentation
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    select case (qflag)
        case ("ql")
            q = ql
        case ("qr")
            q = qr
        case ("qi")
            q = qi
        case ("qs")
            q = qs
        case ("qg")
            q = qg
        case default
            print *, "gfdl_mp: qflag error!"
    end select

    if (use_implicit_fall .eq. use_explicit_fall) then
        write (6, *) 'gfdl_mp: use_implicit_fall and use_explicit_fall cannot be the same.'
        stop
    endif
    
    if (use_ppm) then
        call lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, q, x1, m1)
    else
        if (use_implicit_fall) then
            call implicit_fall (dts, ks, ke, ze, vt, dp, q, x1, m1)
        elseif (use_explicit_fall) then
            call explicit_fall (dts, ks, ke, ze, vt, dp, q, x1, m1)
        endif
    endif
    
    select case (qflag)
        case ("ql")
            ql = q
        case ("qr")
            qr = q
        case ("qi")
            qi = q
        case ("qs")
            qs = q
        case ("qg")
            qg = q
        case default
            print *, "gfdl_mp: qflag error!"
    end select
    
=======

    ! sjl, turn off melting of falling cloud ice, snow and graupel
    ! if (dtm < 60.) k0 = ke
    k0 = ke
    ! sjl, turn off melting of falling cloud ice, snow and graupel

    ze (ke + 1) = zs
    do k = ke, ks, - 1
        ze (k) = ze (k + 1) - dz (k) ! dz < 0
    enddo

    zt (ks) = ze (ks)

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! energy change during sedimentation
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    if (consv_checker) then
        do k = ks, ke
            te2 (k) = mte (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), tz (k), dp (k), .false.)
        enddo
        dte = dte + sum (te1) - sum (te2)
    endif
    
=======

    do k = k0, ke
        icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)
    enddo

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! momentum transportation during sedimentation
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    if (do_sedi_uv) then
        call sedi_uv (ks, ke, m1, dp, u, v)
    endif
    
    if (do_sedi_w) then
        call sedi_w (ks, ke, m1, w, vt, dm)
=======

    call check_column (ks, ke, qi, no_fall)

    if (vi_fac < 1.e-5 .or. no_fall) then
        i1 = 0.
    else

        do k = ks + 1, ke
            zt (k) = ze (k) - dt5 * (vti (k - 1) + vti (k))
        enddo
        zt (ke + 1) = zs - dtm * vti (ke)

        do k = ks, ke
            if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
        enddo

        if (k0 < ke) then
            do k = ke - 1, k0, - 1
                if (qi (k) > qrmin) then
                    do m = k + 1, ke
                        if (zt (k + 1) >= ze (m)) exit
                        if (zt (k) < ze (m + 1) .and. tz (m) > tice) then
                            dtime = min (1.0, (ze (m) - ze (m + 1)) / (max (vr_min, vti (k)) * tau_imlt))
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

        if (use_ppm_ice) then
            call lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, qi, i1, m1_sol, mono_prof)
        else
            call implicit_fall (dtm, ks, ke, ze, vti, dp, qi, i1, m1_sol)
        endif
        
        pfi (ks) = max (0.0, m1_sol (ks))
        do k = ke, ks + 1, -1
            pfi (k) = max (0.0, m1_sol (k) - m1_sol (k - 1))
        enddo

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

>>>>>>> user/lnz/shield2022
    endif

    ! -----------------------------------------------------------------------
    ! energy change during sedimentation heating
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    if (consv_checker) then
=======

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
            if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
        enddo

        if (k0 < ke) then
            do k = ke - 1, k0, - 1
                if (qs (k) > qrmin) then
                    do m = k + 1, ke
                        if (zt (k + 1) >= ze (m)) exit
                        dtime = min (dtm, (ze (m) - ze (m + 1)) / (vr_min + vts (k)))
                        if (zt (k) < ze (m + 1) .and. tz (m) > tice) then
                            dtime = min (1.0, dtime / tau_smlt)
                            sink = min (qs (k) * dp (k) / dp (m), dtime * (tz (m) - tice) / icpk (m))
                            tz (m) = tz (m) - sink * icpk (m)
                            qs (k) = qs (k) - sink * dp (m) / dp (k)
                            if (zt (k) < zs) then
                                r1 = r1 + sink * dp (m) ! precip as rain
                                pfr (k) = sink * dp (m)
                            else
                                ! qr source here will fall next time step (therefore, can evap)
                                qr (m) = qr (m) + sink
                            endif
                        endif
                        if (qs (k) < qrmin) exit
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
        
        pfs (ks) = max (0.0, m1 (ks))
        do k = ke, ks + 1, -1
            pfs (k) = max (0.0, m1 (k) - m1 (k - 1))
        enddo

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

>>>>>>> user/lnz/shield2022
        do k = ks, ke
            te1 (k) = mte (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), tz (k), dp (k), .false.)
        enddo
<<<<<<< HEAD
    endif
    
    ! -----------------------------------------------------------------------
    ! heat exchanges during sedimentation
    ! -----------------------------------------------------------------------
    
    if (do_sedi_heat) then
        call sedi_heat (ks, ke, dp, m1, dz, tz, qv, ql, qr, qi, qs, qg, c_ice)
    endif
    
    ! -----------------------------------------------------------------------
    ! energy change during sedimentation heating
    ! -----------------------------------------------------------------------
    
    if (consv_checker) then
=======

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
            if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
        enddo

        if (k0 < ke) then
            do k = ke - 1, k0, - 1
                if (qg (k) > qrmin) then
                    do m = k + 1, ke
                        if (zt (k + 1) >= ze (m)) exit
                        dtime = min (dtm, (ze (m) - ze (m + 1)) / vtg (k))
                        if (zt (k) < ze (m + 1) .and. tz (m) > tice) then
                            dtime = min (1., dtime / tau_g2r)
                            sink = min (qg (k) * dp (k) / dp (m), dtime * (tz (m) - tice) / icpk (m))
                            tz (m) = tz (m) - sink * icpk (m)
                            qg (k) = qg (k) - sink * dp (m) / dp (k)
                            if (zt (k) < zs) then
                                r1 = r1 + sink * dp (m)
                                pfr (k) = sink * dp (m)
                            else
                                qr (m) = qr (m) + sink
                            endif
                        endif
                        if (qg (k) < qrmin) exit
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
        
        pfg (ks) = max (0.0, m1 (ks))
        do k = ke, ks + 1, -1
            pfg (k) = max (0.0, m1 (k) - m1 (k - 1))
        enddo
    
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

>>>>>>> user/lnz/shield2022
        do k = ks, ke
            te2 (k) = mte (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), tz (k), dp (k), .false.)
        enddo
<<<<<<< HEAD
        dte = dte + sum (te1) - sum (te2)
=======

        if (do_sedi_w) then
            w1 (ks) = w1 (ks) + m1 (ks) * vtg (ks) / dm (ks)
            do k = ks + 1, ke
                w1 (k) = (dm (k) * w1 (k) + m1 (k - 1) * (w1 (k - 1) - vtg (k - 1)) + m1 (k) * vtg (k)) &
                     / (dm (k) + m1 (k - 1))
            enddo
        endif

>>>>>>> user/lnz/shield2022
    endif

end subroutine terminal_fall

! =======================================================================
! calculate ze zt for sedimentation
! =======================================================================

subroutine zezt (ks, ke, dts, zs, dz, vt, ze, zt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: dz, vt
    
    real, intent (out) :: zs
    
    real, intent (out), dimension (ks:ke + 1) :: ze, zt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: dt5
    
    dt5 = 0.5 * dts
    zs = 0.0
    ze (ke + 1) = zs
    do k = ke, ks, - 1
        ze (k) = ze (k + 1) - dz (k)
    enddo
    zt (ks) = ze (ks)
    do k = ks + 1, ke
        zt (k) = ze (k) - dt5 * (vt (k - 1) + vt (k))
    enddo
    zt (ke + 1) = zs - dts * vt (ke)
    do k = ks, ke
        if (zt (k + 1) .ge. zt (k)) zt (k + 1) = zt (k) - dz_min
    enddo
    
end subroutine zezt

! =======================================================================
! check if water species is large enough to fall
! =======================================================================

subroutine check_column (ks, ke, q, no_fall)

    implicit none
<<<<<<< HEAD
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
=======

>>>>>>> user/lnz/shield2022
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
! warm rain cloud microphysics
! =======================================================================

<<<<<<< HEAD
subroutine warm_rain (dts, ks, ke, dp, dz, tz, qv, ql, qr, qi, qs, qg, &
        den, denfac, vtw, vtr, ccn, rh_rain, h_var, reevap)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts, rh_rain, h_var
    
    real, intent (in), dimension (ks:ke) :: dp, dz, den, denfac, vtw, vtr
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, ccn
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    real, intent (out) :: reevap
    
=======
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

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    reevap = 0
    
=======

    qm (ks) = q (ks) / (dz (ks) + dd (ks))
    do k = ks + 1, ke
        qm (k) = (q (k) + dd (k - 1) * qm (k - 1)) / (dz (k) + dd (k))
    enddo

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! rain evaporation to form water vapor
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    call prevp (ks, ke, dts, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var, dp, reevap)
    
=======

    do k = ks, ke
        qm (k) = qm (k) * dz (k)
    enddo

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! rain accretion with cloud water
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    call pracw (ks, ke, dts, tz, qv, ql, qr, qi, qs, qg, den, denfac, vtw, vtr)
    
=======

    m1 (ks) = q (ks) - qm (ks)
    do k = ks + 1, ke
        m1 (k) = m1 (k - 1) + q (k) - qm (k)
    enddo
    precip = m1 (ke)

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! cloud water to rain autoconversion
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    call praut (ks, ke, dts, tz, qv, ql, qr, qi, qs, qg, den, ccn, h_var)
    
end subroutine warm_rain
=======

    do k = ks, ke
        q (k) = qm (k) / dp (k) !dry dp used inside MP
    enddo

end subroutine implicit_fall
>>>>>>> user/lnz/shield2022

! =======================================================================
! rain evaporation to form water vapor, Lin et al. (1983)
! =======================================================================

<<<<<<< HEAD
subroutine prevp (ks, ke, dts, tz, qv, ql, qr, qi, qs, qg, den, denfac, rh_rain, h_var, dp, reevap)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts, rh_rain, h_var
    
    real, intent (in), dimension (ks:ke) :: den, denfac, dp
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    real, intent (inout), dimension (ks:ke) :: qv, qr, ql, qi, qs, qg
    
    real, intent (out) :: reevap
    
=======
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

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    integer :: k
    
    real :: dqv, qsat, dqdt, tmp, t2, qden, q_plus, q_minus, sink
    real :: qpz, dq, dqh, tin, fac_revp, rh_tem
    
    real, dimension (ks:ke) :: q_liq, q_sol, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), dimension (ks:ke) :: cvm, te8
    
=======

    do k = ks, ke
        dz (k) = zt (k) - zt (k + 1) ! note: dz is positive
        q (k) = q (k) * dp (k)
        a4 (1, k) = q (k) / dz (k)
        qm (k) = 0.
    enddo

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    reevap = 0
    
    ! -----------------------------------------------------------------------
    ! time-scale factor
    ! -----------------------------------------------------------------------
    
    fac_revp = 1.
    if (tau_revp .gt. 1.e-6) then
        fac_revp = 1. - exp (- dts / tau_revp)
    endif
    
    ! -----------------------------------------------------------------------
    ! calculate heat capacities and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    call cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, cvm, te8, tz, &
        lcpk, icpk, tcpk, tcp3)
    
    do k = ks, ke
        
        tin = (tz (k) * cvm (k) - lv00 * ql (k)) / mhc (qv (k) + ql (k), qr (k), q_sol (k))
        
        ! -----------------------------------------------------------------------
        ! calculate supersaturation and subgrid variability of water
        ! -----------------------------------------------------------------------
        
        qpz = qv (k) + ql (k)
        qsat = wqs (tin, den (k), dqdt)
        dqv = qsat - qv (k)
        
        dqh = max (ql (k), h_var * max (qpz, qcmin))
        dqh = min (dqh, 0.2 * qpz)
        q_minus = qpz - dqh
        q_plus = qpz + dqh
        
        ! -----------------------------------------------------------------------
        ! rain evaporation
        ! -----------------------------------------------------------------------
        
        rh_tem = qpz / qsat
            
        if (tz (k) .gt. t_wfr .and. qr (k) .gt. qcmin .and. dqv .gt. 0.0 .and. qsat .gt. q_minus) then

            if (qsat .gt. q_plus) then
                dq = qsat - qpz
            else
                dq = 0.25 * (qsat - q_minus) ** 2 / dqh
            endif
            qden = qr (k) * den (k)
            t2 = tin * tin
            sink = psub (t2, dq, qden, qsat, crevp, den (k), denfac (k), blinr, mur, lcpk (k), cvm (k))
            sink = min (qr (k), dts * fac_revp * sink, dqv / (1. + lcpk (k) * dqdt))
            if (use_rhc_revap .and. rh_tem .ge. rhc_revap) then
                sink = 0.0
            endif
            
            ! -----------------------------------------------------------------------
            ! alternative minimum evaporation in dry environmental air
            ! -----------------------------------------------------------------------
            ! tmp = min (qr (k), dim (rh_rain * qsat, qv (k)) / (1. + lcpk (k) * dqdt))
            ! sink = max (sink, tmp)
            
            reevap = reevap + sink * dp (k)
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                sink, 0., - sink, 0., 0., 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo ! k loop
    
end subroutine prevp

! =======================================================================
! rain accretion with cloud water, Lin et al. (1983)
! =======================================================================

subroutine pracw (ks, ke, dts, tz, qv, ql, qr, qi, qs, qg, den, denfac, vtw, vtr)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vtw, vtr
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    real, intent (inout), dimension (ks:ke) :: qv, qr, ql, qi, qs, qg
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: qden, sink
    
    do k = ks, ke
        
        if (tz (k) .gt. t_wfr .and. qr (k) .gt. qcmin .and. ql (k) .gt. qcmin) then
            
            qden = qr (k) * den (k)
            if (do_new_acc_water) then
                sink = dts * acr3d (vtr (k), vtw (k), ql (k), qr (k), cracw, acco (:, 5), &
                    acc (9), acc (10), den (k))
            else
                sink = dts * acr2d (qden, cracw, denfac (k), blinr, mur)
                sink = sink / (1. + sink) * ql (k)
            endif
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - sink, sink, 0., 0., 0.)
            
        endif
        
    enddo
    
end subroutine pracw

! =======================================================================
! cloud water to rain autoconversion, Hong et al. (2004)
! =======================================================================

subroutine praut (ks, ke, dts, tz, qv, ql, qr, qi, qs, qg, den, ccn, h_var)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts, h_var
    
    real, intent (in), dimension (ks:ke) :: den
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, ccn
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real, parameter :: so3 = 7.0 / 3.0
    real, parameter :: so1 = - 1.0 / 3.0
    
    integer :: k
    
    real :: sink, dq, qc
    
    real, dimension (ks:ke) :: dl, c_praut

    if (irain_f .eq. 0) then
        
        call linear_prof (ke - ks + 1, ql (ks), dl (ks), z_slope_liq, h_var)
        
        do k = ks, ke
            
            if (tz (k) .gt. t_wfr .and. ql (k) .gt. qcmin) then
                
                if (do_psd_water_num) then
                    ccn (k) = coeaw / coebw * exp (muw / (muw + 3) * log (den (k) * ql (k))) / den (k)
                endif

                qc = fac_rc * ccn (k)
                dl (k) = min (max (qcmin, dl (k)), 0.5 * ql (k))
                dq = 0.5 * (ql (k) + dl (k) - qc)
                
                if (dq .gt. 0.) then
                    
                    c_praut (k) = cpaut * exp (so1 * log (ccn (k) * rhow))
                    sink = min (1., dq / dl (k)) * dts * c_praut (k) * den (k) * &
                        exp (so3 * log (ql (k)))
                    sink = min (ql (k), sink)
                    
                    call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                        0., - sink, sink, 0., 0., 0.)
                    
                endif
                
            endif
            
        enddo
        
    endif
    
    if (irain_f .eq. 1) then
        
        do k = ks, ke
            
            if (tz (k) .gt. t_wfr .and. ql (k) .gt. qcmin) then
                
                if (do_psd_water_num) then
                    ccn (k) = coeaw / coebw * exp (muw / (muw + 3) * log (den (k) * ql (k))) / den (k)
                endif

                qc = fac_rc * ccn (k)
                dq = ql (k) - qc
                
                if (dq .gt. 0.) then
                    
                    c_praut (k) = cpaut * exp (so1 * log (ccn (k) * rhow))
                    sink = min (dq, dts * c_praut (k) * den (k) * exp (so3 * log (ql (k))))
                    sink = min (ql (k), sink)
                    
                    call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                        0., - sink, sink, 0., 0., 0.)
                    
                endif
                
            endif
            
        enddo
        
    endif

end subroutine praut
    
! =======================================================================
! ice cloud microphysics
! =======================================================================

subroutine ice_cloud (ks, ke, tz, qv, ql, qr, qi, qs, qg, den, &
        denfac, vtw, vtr, vti, vts, vtg, dts, h_var)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts, h_var
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vtw, vtr, vti, vts, vtg
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real, dimension (ks:ke) :: di, q_liq, q_sol, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), dimension (ks:ke) :: cvm, te8
    
    ! -----------------------------------------------------------------------
    ! calculate heat capacities and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    call cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, cvm, te8, tz, &
        lcpk, icpk, tcpk, tcp3)
    
    if (.not. do_warm_rain_mp) then
        
        ! -----------------------------------------------------------------------
        ! cloud ice melting to form cloud water and rain
        ! -----------------------------------------------------------------------
        
        call pimlt (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! cloud water freezing to form cloud ice and snow
        ! -----------------------------------------------------------------------
        
        call pifr (ks, ke, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! vertical subgrid variability
        ! -----------------------------------------------------------------------
        
        call linear_prof (ke - ks + 1, qi, di, z_slope_ice, h_var)
        
        ! -----------------------------------------------------------------------
        ! snow melting (includes snow accretion with cloud water and rain) to form cloud water and rain
        ! -----------------------------------------------------------------------
        
        call psmlt (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
            vtw, vtr, vts, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! graupel melting (includes graupel accretion with cloud water and rain) to form rain
        ! -----------------------------------------------------------------------
        
        call pgmlt (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
            vtw, vtr, vtg, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! snow accretion with cloud ice
        ! -----------------------------------------------------------------------
        
        call psaci (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, denfac, vti, vts)
        
        ! -----------------------------------------------------------------------
        ! cloud ice to snow autoconversion
        ! -----------------------------------------------------------------------
        
        call psaut (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, di)
        
        ! -----------------------------------------------------------------------
        ! graupel accretion with cloud ice
        ! -----------------------------------------------------------------------
        
        call pgaci (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, denfac, vti, vtg)
        
        ! -----------------------------------------------------------------------
        ! snow accretion with rain and rain freezing to form graupel
        ! -----------------------------------------------------------------------
        
        call psacr_pgfr (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
            vtr, vts, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! graupel accretion with snow
        ! -----------------------------------------------------------------------
        
        call pgacs (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, vts, vtg)
        
        ! -----------------------------------------------------------------------
        ! snow to graupel autoconversion
        ! -----------------------------------------------------------------------
        
        call pgaut (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den)
        
        ! -----------------------------------------------------------------------
        ! graupel accretion with cloud water and rain
        ! -----------------------------------------------------------------------
        
        call pgacw_pgacr (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
            vtr, vtg, lcpk, icpk, tcpk, tcp3)
        
    endif ! do_warm_rain_mp
    
end subroutine ice_cloud

! =======================================================================
! cloud ice melting to form cloud water and rain, Lin et al. (1983)
! =======================================================================

subroutine pimlt (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, tmp, sink, fac_imlt
    
    fac_imlt = 1. - exp (- dts / tau_imlt)
    
    do k = ks, ke
        
        tc = tz (k) - tice_mlt
        
        if (tc .gt. 0 .and. qi (k) .gt. qcmin) then
            
            sink = fac_imlt * tc / icpk (k)
            sink = min (qi (k), sink)
            tmp = min (sink, dim (ql_mlt, ql (k)))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., tmp, sink - tmp, - sink, 0., 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pimlt

! =======================================================================
! cloud water freezing to form cloud ice and snow, Lin et al. (1983)
! =======================================================================

subroutine pifr (ks, ke, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in), dimension (ks:ke) :: den
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, tmp, sink, qim
    
    do k = ks, ke
        
        tc = t_wfr - tz (k)
        
        if (tc .gt. 0. .and. ql (k) .gt. qcmin) then
            
            sink = ql (k) * tc / dt_fr
            sink = min (ql (k), sink, tc / icpk (k))
            qim = qi0_crt / den (k)
            tmp = min (sink, dim (qim, qi (k)))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - sink, 0., tmp, sink - tmp, 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pifr

! =======================================================================
! snow melting (includes snow accretion with cloud water and rain) to form cloud water and rain
! Lin et al. (1983)
! =======================================================================

subroutine psmlt (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
        vtw, vtr, vts, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vtw, vtr, vts
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, factor, tmp, sink, qden, dqdt, tin, dq
    real :: psacw, psacr, pracs
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .ge. 0. .and. qs (k) .gt. qcmin) then
            
            psacw = 0.
            qden = qs (k) * den (k)
            if (ql (k) .gt. qcmin) then
                if (do_new_acc_water) then
                    psacw = acr3d (vts (k), vtw (k), ql (k), qs (k), csacw, acco (:, 7), &
                        acc (13), acc (14), den (k))
                else
                    factor = acr2d (qden, csacw, denfac (k), blins, mus)
                    psacw = factor / (1. + dts * factor) * ql (k)
                endif
            endif
            
            psacr = 0.
            pracs = 0.
            if (qr (k) .gt. qcmin) then
                psacr = min (acr3d (vts (k), vtr (k), qr (k), qs (k), csacr, acco (:, 2), &
                    acc (3), acc (4), den (k)), qr (k) / dts)
                pracs = acr3d (vtr (k), vts (k), qs (k), qr (k), cracs, acco (:, 1), &
                    acc (1), acc (2), den (k))
            endif
            
            tin = tz (k)
            dq = iqs (tin, den (k), dqdt) - qv (k)
            sink = max (0., pmlt (tc, dq, qden, psacw, psacr, csmlt, den (k), denfac (k), blins, mus, &
                lcpk (k), icpk (k), cvm (k)))
            
            sink = min (qs (k), (sink + pracs) * dts, tc / icpk (k))
            tmp = min (sink, dim (qs_mlt, ql (k)))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., tmp, sink - tmp, 0., - sink, 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine psmlt

! =======================================================================
! graupel melting (includes graupel accretion with cloud water and rain) to form rain
! Lin et al. (1983)
! =======================================================================

subroutine pgmlt (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
        vtw, vtr, vtg, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vtw, vtr, vtg
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, factor, sink, qden, dqdt, tin, dq
    real :: pgacw, pgacr
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .ge. 0. .and. qg (k) .gt. qcmin) then
            
            pgacw = 0.
            qden = qg (k) * den (k)
            if (ql (k) .gt. qcmin) then
                if (do_new_acc_water) then
                    pgacw = acr3d (vtg (k), vtw (k), ql (k), qg (k), cgacw, acco (:, 9), &
                        acc (17), acc (18), den (k))
                else
                    if (do_hail) then
                        factor = acr2d (qden, cgacw, denfac (k), blinh, muh)
                    else
                        factor = acr2d (qden, cgacw, denfac (k), bling, mug)
                    endif
                    pgacw = factor / (1. + dts * factor) * ql (k)
                endif
            endif
            
            pgacr = 0.
            if (qr (k) .gt. qcmin) then
                pgacr = min (acr3d (vtg (k), vtr (k), qr (k), qg (k), cgacr, acco (:, 3), &
                    acc (5), acc (6), den (k)), qr (k) / dts)
            endif
            
            tin = tz (k)
            dq = iqs (tin, den (k), dqdt) - qv (k)
            if (do_hail) then
                sink = max (0., pmlt (tc, dq, qden, pgacw, pgacr, cgmlt, den (k), denfac (k), blinh, muh, &
                    lcpk (k), icpk (k), cvm (k)))
            else
                sink = max (0., pmlt (tc, dq, qden, pgacw, pgacr, cgmlt, den (k), denfac (k), bling, mug, &
                    lcpk (k), icpk (k), cvm (k)))
            endif
            
            sink = min (qg (k), sink * dts, tc / icpk (k))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., sink, 0., 0., - sink, te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pgmlt

! =======================================================================
! snow accretion with cloud ice, Lin et al. (1983)
! =======================================================================

subroutine psaci (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, denfac, vti, vts)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vti, vts
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, factor, sink, qden
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .lt. 0. .and. qi (k) .gt. qcmin) then
            
            sink = 0.
            qden = qs (k) * den (k)
            if (qs (k) .gt. qcmin) then
                if (do_new_acc_ice) then
                    sink = dts * acr3d (vts (k), vti (k), qi (k), qs (k), csaci, acco (:, 8), &
                        acc (15), acc (16), den (k))
                else
                    factor = dts * acr2d (qden, csaci, denfac (k), blins, mus)
                    sink = factor / (1. + factor) * qi (k)
                endif
            endif
            
            sink = min (fi2s_fac * qi (k), sink)
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., - sink, sink, 0.)
            
        endif
        
    enddo
    
end subroutine psaci

! =======================================================================
! cloud ice to snow autoconversion, Lin et al. (1983)
! =======================================================================

subroutine psaut (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, di)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, di
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, sink, fac_i2s, q_plus, qim, dq, tmp
    
    fac_i2s = 1. - exp (- dts / tau_i2s)
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .lt. 0. .and. qi (k) .gt. qcmin) then
            
            sink = 0.
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
                sink = tmp * dq
            endif
            
            sink = min (fi2s_fac * qi (k), sink)
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., - sink, sink, 0.)
            
        endif
        
    enddo
    
end subroutine psaut

! =======================================================================
! graupel accretion with cloud ice, Lin et al. (1983)
! =======================================================================

subroutine pgaci (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, denfac, vti, vtg)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vti, vtg
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, factor, sink, qden
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .lt. 0. .and. qi (k) .gt. qcmin) then
            
            sink = 0.
            qden = qg (k) * den (k)
            if (qg (k) .gt. qcmin) then
                if (do_new_acc_ice) then
                    sink = dts * acr3d (vtg (k), vti (k), qi (k), qg (k), cgaci, acco (:, 10), &
                        acc (19), acc (20), den (k))
                else
                    if (do_hail) then
                        factor = dts * acr2d (qden, cgaci, denfac (k), blinh, muh)
                    else
                        factor = dts * acr2d (qden, cgaci, denfac (k), bling, mug)
                    endif
                    sink = factor / (1. + factor) * qi (k)
                endif
            endif
            
            sink = min (fi2g_fac * qi (k), sink)
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., - sink, 0., sink)
            
        endif
        
    enddo
    
end subroutine pgaci

! =======================================================================
! snow accretion with rain and rain freezing to form graupel, Lin et al. (1983)
! =======================================================================

subroutine psacr_pgfr (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
        vtr, vts, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vtr, vts
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, factor, sink
    real :: psacr, pgfr
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .lt. 0. .and. qr (k) .gt. qcmin) then
            
            psacr = 0.
            if (qs (k) .gt. qcmin) then
                psacr = dts * acr3d (vts (k), vtr (k), qr (k), qs (k), csacr, acco (:, 2), &
                    acc (3), acc (4), den (k))
            endif
            
            pgfr = dts * cgfr (1) / den (k) * (exp (- cgfr (2) * tc) - 1.) * &
                exp (1. / (mur + 3) * (6 + mur) * log (6 * qr (k) * den (k)))
            
            sink = psacr + pgfr
            factor = min (sink, qr (k), - tc / icpk (k)) / max (sink, qcmin)
            psacr = factor * psacr
            pgfr = factor * pgfr
            
            sink = min (qr (k), psacr + pgfr)
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., - sink, 0., psacr, pgfr, te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine psacr_pgfr

! =======================================================================
! graupel accretion with snow, Lin et al. (1983)
! =======================================================================

subroutine pgacs (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den, vts, vtg)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, vts, vtg
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: sink
    
    do k = ks, ke
        
        if (tz (k) .lt. tice .and. qs (k) .gt. qcmin .and. qg (k) .gt. qcmin) then
            
            sink = dts * acr3d (vtg (k), vts (k), qs (k), qg (k), cgacs, acco (:, 4), &
                acc (7), acc (8), den (k))
            sink = min (fs2g_fac * qs (k), sink)
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., 0., - sink, sink)
            
        endif
        
    enddo
    
end subroutine pgacs

! =======================================================================
! snow to graupel autoconversion, Lin et al. (1983)
! =======================================================================

subroutine pgaut (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, factor, sink, qsm
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .lt. 0. .and. qs (k) .gt. qcmin) then
            
            sink = 0
            qsm = qs0_crt / den (k)
            if (qs (k) .gt. qsm) then
                factor = dts * 1.e-3 * exp (0.09 * (tz (k) - tice))
                sink = factor / (1. + factor) * (qs (k) - qsm)
            endif
            
            sink = min (fs2g_fac * qs (k), sink)
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., 0., - sink, sink)
            
        endif
        
    enddo
    
end subroutine pgaut

! =======================================================================
! graupel accretion with cloud water and rain, Lin et al. (1983)
! =======================================================================

subroutine pgacw_pgacr (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, denfac, &
        vtr, vtg, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, denfac, vtr, vtg
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, factor, sink, qden
    real :: pgacw, pgacr
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .lt. 0. .and. qg (k) .gt. qcmin) then
            
            pgacw = 0.
            if (ql (k) .gt. qcmin) then
                qden = qg (k) * den (k)
                if (do_hail) then
                    factor = dts * acr2d (qden, cgacw, denfac (k), blinh, muh)
                else
                    factor = dts * acr2d (qden, cgacw, denfac (k), bling, mug)
                endif
                pgacw = factor / (1. + factor) * ql (k)
            endif
            
            pgacr = 0.
            if (qr (k) .gt. qcmin) then
                pgacr = min (dts * acr3d (vtg (k), vtr (k), qr (k), qg (k), cgacr, acco (:, 3), &
                    acc (5), acc (6), den (k)), qr (k))
            endif
            
            sink = pgacr + pgacw
            factor = min (sink, dim (tice, tz (k)) / icpk (k)) / max (sink, qcmin)
            pgacr = factor * pgacr
            pgacw = factor * pgacw
            
            sink = pgacr + pgacw
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - pgacw, - pgacr, 0., 0., sink, te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pgacw_pgacr

! =======================================================================
! temperature sentive high vertical resolution processes
! =======================================================================

subroutine subgrid_z_proc (ks, ke, den, denfac, dts, rh_adj, tz, qv, ql, qr, &
        qi, qs, qg, dp, ccn, cin, cond, dep, reevap, sub)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts, rh_adj
    
    real, intent (in), dimension (ks:ke) :: den, denfac, dp
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, ccn, cin
    
    real, intent (out) :: cond, dep, reevap, sub
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real, dimension (ks:ke) :: q_liq, q_sol, q_cond, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), dimension (ks:ke) :: cvm, te8
    
    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
    
    cond = 0
    dep = 0
    reevap = 0
    sub = 0
    
    ! -----------------------------------------------------------------------
    ! calculate heat capacities and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    call cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, cvm, te8, tz, &
        lcpk, icpk, tcpk, tcp3)
    
    ! -----------------------------------------------------------------------
    ! instant processes (include deposition, evaporation, and sublimation)
    ! -----------------------------------------------------------------------
    
    if (.not. do_warm_rain_mp) then
        
        call pinst (ks, ke, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
            lcpk, icpk, tcpk, tcp3, rh_adj, dep, sub, reevap)
        
    endif
    
    ! -----------------------------------------------------------------------
    ! cloud water condensation and evaporation
    ! -----------------------------------------------------------------------
    
    call pcond_pevap (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
        lcpk, icpk, tcpk, tcp3, cond, reevap)
    
    if (.not. do_warm_rain_mp) then
        
        ! -----------------------------------------------------------------------
        ! enforce complete freezing below t_wfr
        ! -----------------------------------------------------------------------
        
        call pcomp (ks, ke, qv, ql, qr, qi, qs, qg, tz, cvm, te8, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! Wegener Bergeron Findeisen process
        ! -----------------------------------------------------------------------
        
        call pwbf (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! Bigg freezing mechanism
        ! -----------------------------------------------------------------------
        
        call pbigg (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, ccn, lcpk, icpk, tcpk, tcp3)
        
        ! -----------------------------------------------------------------------
        ! cloud ice deposition and sublimation
        ! -----------------------------------------------------------------------
        
        call pidep_pisub (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
            lcpk, icpk, tcpk, tcp3, cin, dep, sub)
        
        ! -----------------------------------------------------------------------
        ! snow deposition and sublimation
        ! -----------------------------------------------------------------------
        
        call psdep_pssub (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
            denfac, lcpk, icpk, tcpk, tcp3, dep, sub)
        
        ! -----------------------------------------------------------------------
        ! graupel deposition and sublimation
        ! -----------------------------------------------------------------------
        
        call pgdep_pgsub (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
            denfac, lcpk, icpk, tcpk, tcp3, dep, sub)
        
    endif
    
end subroutine subgrid_z_proc

! =======================================================================
! instant processes (include deposition, evaporation, and sublimation)
! =======================================================================

subroutine pinst (ks, ke, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
        lcpk, icpk, tcpk, tcp3, rh_adj, dep, sub, reevap)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: rh_adj
    
    real, intent (in), dimension (ks:ke) :: den, dp
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    real, intent (out) :: dep, reevap, sub
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: sink, tin, qpz, rh, dqdt, tmp
    
    do k = ks, ke
        
        ! -----------------------------------------------------------------------
        ! instant deposit all water vapor to cloud ice when temperature is super low
        ! -----------------------------------------------------------------------
        
        if (tz (k) .lt. t_min) then
            
            sink = dim (qv (k), qcmin)
            dep = dep + sink * dp (k)
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                 - sink, 0., 0., sink, 0., 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
        ! -----------------------------------------------------------------------
        ! instant evaporation / sublimation of all clouds when rh < rh_adj
        ! -----------------------------------------------------------------------
        
        qpz = qv (k) + ql (k) + qi (k)
        tin = (te8 (k) - lv00 * qpz + li00 * (qs (k) + qg (k))) / &
            mhc (qpz, qr (k), qs (k) + qg (k))
        
        if (tin .gt. t_sub + 6.) then
            
            rh = qpz / iqs (tin, den (k), dqdt)
            if (rh .lt. rh_adj) then
                
                sink = ql (k)
                tmp = qi (k)
                
                reevap = reevap + sink * dp (k)
                sub = sub + tmp * dp (k)
                
                call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                    sink + tmp, - sink, 0., - tmp, 0., 0., te8 (k), cvm (k), tz (k), &
                    lcpk (k), icpk (k), tcpk (k), tcp3 (k))
                
            endif
            
        endif
        
    enddo
    
end subroutine pinst

! =======================================================================
! cloud water condensation and evaporation, Hong and Lim (2006)
! =======================================================================

subroutine pcond_pevap (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
        lcpk, icpk, tcpk, tcp3, cond, reevap)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, dp
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    real, intent (out) :: cond, reevap
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: sink, tin, qpz, dqdt, qsw, rh_tem, dq, factor, fac_l2v, fac_v2l
    
    fac_l2v = 1. - exp (- dts / tau_l2v)
    fac_v2l = 1. - exp (- dts / tau_v2l)
    
    do k = ks, ke
        
        tin = tz (k)
        qsw = wqs (tin, den (k), dqdt)
        qpz = qv (k) + ql (k) + qi (k)
        rh_tem = qpz / qsw
        dq = qsw - qv (k)
        if (dq .gt. 0.) then
            factor = min (1., fac_l2v * (rh_fac * dq / qsw))
            sink = min (ql (k), factor * dq / (1. + tcp3 (k) * dqdt))
            if (use_rhc_cevap .and. rh_tem .ge. rhc_cevap) then
                sink = 0.
            endif
            reevap = reevap + sink * dp (k)
        elseif (do_cond_timescale) then
            factor = min (1., fac_v2l * (rh_fac * (- dq) / qsw))
            sink = - min (qv (k), factor * (- dq) / (1. + tcp3 (k) * dqdt))
            cond = cond - sink * dp (k)
        else
            sink = - min (qv (k), - dq / (1. + tcp3 (k) * dqdt))
            cond = cond - sink * dp (k)
        endif
        
        call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
            sink, - sink, 0., 0., 0., 0., te8 (k), cvm (k), tz (k), &
            lcpk (k), icpk (k), tcpk (k), tcp3 (k))
        
    enddo
    
end subroutine pcond_pevap

! =======================================================================
! enforce complete freezing below t_wfr, Lin et al. (1983)
! =======================================================================

subroutine pcomp (ks, ke, qv, ql, qr, qi, qs, qg, tz, cvm, te8, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, sink
    
    do k = ks, ke
        
        tc = t_wfr - tz (k)
        
        if (tc .gt. 0. .and. ql (k) .gt. qcmin) then
            
            sink = ql (k) * tc / dt_fr
            sink = min (ql (k), sink, tc / icpk (k))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - sink, 0., sink, 0., 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pcomp

! =======================================================================
! Wegener Bergeron Findeisen process, Storelvmo and Tan (2015)
! =======================================================================

subroutine pwbf (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, tin, sink, dqdt, qsw, qsi, qim, tmp, fac_wbf

    if (.not. do_wbf) return
    
    fac_wbf = 1. - exp (- dts / tau_wbf)
    
    do k = ks, ke
        
        tc = tice - tz (k)
        
        tin = tz (k)
        qsw = wqs (tin, den (k), dqdt)
        qsi = iqs (tin, den (k), dqdt)

        if (tc .gt. 0. .and. ql (k) .gt. qcmin .and. qi (k) .gt. qcmin .and. &
            qv (k) .gt. qsi .and. qv (k) .lt. qsw) then

            sink = min (fac_wbf * ql (k), tc / icpk (k))
            qim = qi0_crt / den (k)
            tmp = min (sink, dim (qim, qi (k)))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - sink, 0., tmp, sink - tmp, 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pwbf

! =======================================================================
! Bigg freezing mechanism, Bigg (1953)
! =======================================================================

subroutine pbigg (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, den, ccn, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, ccn
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: sink, tc
    
    do k = ks, ke
        
        tc = tice - tz (k)
        
        if (tc .gt. 0 .and. ql (k) .gt. qcmin) then
            
            if (do_psd_water_num) then
                ccn (k) = coeaw / coebw * exp (muw / (muw + 3) * log (den (k) * ql (k))) / den (k)
            endif

            sink = 100. / (rhow * ccn (k)) * dts * (exp (0.66 * tc) - 1.) * ql (k) ** 2
            sink = min (ql (k), sink, tc / icpk (k))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - sink, 0., sink, 0., 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo

end subroutine pbigg
        
! =======================================================================
! cloud ice deposition and sublimation, Hong et al. (2004)
! =======================================================================

subroutine pidep_pisub (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
        lcpk, icpk, tcpk, tcp3, cin, dep, sub)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, dp
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, cin
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    real, intent (out) :: dep, sub
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: sink, tin, dqdt, qsi, dq, pidep, tmp, tc, qi_gen, qi_crt
    
    do k = ks, ke
        
        if (tz (k) .lt. tice) then
            
            pidep = 0.
            tin = tz (k)
            qsi = iqs (tin, den (k), dqdt)
            dq = qv (k) - qsi
            tmp = dq / (1. + tcpk (k) * dqdt)
            
            if (qi (k) .gt. qcmin) then
                 if (do_psd_ice_num) then
                     cin (k) = coeai / coebi * exp (mui / (mui + 3) * log (den (k) * qi (k))) / den (k)
                 endif
                if (.not. prog_ccn) then
                    if (inflag .eq. 1) &
                        cin (k) = 5.38e7 * exp (0.75 * log (qi (k) * den (k)))
                    if (inflag .eq. 2) &
                        cin (k) = exp (- 2.80 + 0.262 * (tice - tz (k))) * 1000.0
                    if (inflag .eq. 3) &
                        cin (k) = exp (- 0.639 + 12.96 * (qv (k) / qsi - 1.0)) * 1000.0
                    if (inflag .eq. 4) &
                        cin (k) = 5.e-3 * exp (0.304 * (tice - tz (k))) * 1000.0
                    if (inflag .eq. 5) &
                        cin (k) = 1.e-5 * exp (0.5 * (tice - tz (k))) * 1000.0
                endif
                pidep = dts * dq * 4.0 * 11.9 * exp (0.5 * log (qi (k) * den (k) * cin (k))) / &
                     (qsi * den (k) * (tcpk (k) * cvm (k)) ** 2 / (tcond * rvgas * tz (k) ** 2) + &
                    1 / vdifu)
            endif
            
            if (dq .gt. 0.) then
                tc = tice - tz (k)
                qi_gen = 4.92e-11 * exp (1.33 * log (1.e3 * exp (0.1 * tc)))
                if (igflag .eq. 1) &
                    qi_crt = qi_gen / den (k)
                if (igflag .eq. 2) &
                    qi_crt = qi_gen * min (qi_lim, 0.1 * tc) / den (k)
                if (igflag .eq. 3) &
                    qi_crt = 1.82e-6 * min (qi_lim, 0.1 * tc) / den (k)
                if (igflag .eq. 4) &
                    qi_crt = max (qi_gen, 1.82e-6) * min (qi_lim, 0.1 * tc) / den (k)
                sink = min (tmp, max (qi_crt - qi (k), pidep), tc / tcpk (k))
                dep = dep + sink * dp (k)
            else
                pidep = pidep * min (1., dim (tz (k), t_sub) * is_fac)
                sink = max (pidep, tmp, - qi (k))
                sub = sub - sink * dp (k)
            endif
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                 - sink, 0., 0., sink, 0., 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pidep_pisub

! =======================================================================
! snow deposition and sublimation, Lin et al. (1983)
! =======================================================================

subroutine psdep_pssub (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
        denfac, lcpk, icpk, tcpk, tcp3, dep, sub)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, dp, denfac
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    real, intent (out) :: dep, sub
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: sink, tin, dqdt, qsi, qden, t2, dq, pssub
    
    do k = ks, ke
        
        if (qs (k) .gt. qcmin) then
            
            tin = tz (k)
            qsi = iqs (tin, den (k), dqdt)
            qden = qs (k) * den (k)
            t2 = tz (k) * tz (k)
            dq = qsi - qv (k)
            pssub = psub (t2, dq, qden, qsi, cssub, den (k), denfac (k), blins, mus, tcpk (k), cvm (k))
            pssub = dts * pssub
            dq = dq / (1. + tcpk (k) * dqdt)
            if (pssub .gt. 0.) then
                sink = min (pssub * min (1., dim (tz (k), t_sub) * ss_fac), qs (k))
                sub = sub + sink * dp (k)
            else
                sink = 0.
                if (tz (k) .le. tice) then
                    sink = max (pssub, dq, (tz (k) - tice) / tcpk (k))
                endif
                dep = dep - sink * dp (k)
            endif
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                sink, 0., 0., 0., - sink, 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine psdep_pssub

! =======================================================================
! graupel deposition and sublimation, Lin et al. (1983)
! =======================================================================

subroutine pgdep_pgsub (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, dp, cvm, te8, den, &
        denfac, lcpk, icpk, tcpk, tcp3, dep, sub)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den, dp, denfac
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    real, intent (out) :: dep, sub
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: sink, tin, dqdt, qsi, qden, t2, dq, pgsub
    
    do k = ks, ke
        
        if (qg (k) .gt. qcmin) then
            
            tin = tz (k)
            qsi = iqs (tin, den (k), dqdt)
            qden = qg (k) * den (k)
            t2 = tz (k) * tz (k)
            dq = qsi - qv (k)
            if (do_hail) then
                pgsub = psub (t2, dq, qden, qsi, cgsub, den (k), denfac (k), blinh, muh, tcpk (k), cvm (k))
            else
                pgsub = psub (t2, dq, qden, qsi, cgsub, den (k), denfac (k), bling, mug, tcpk (k), cvm (k))
            endif
            pgsub = dts * pgsub
            dq = dq / (1. + tcpk (k) * dqdt)
            if (pgsub .gt. 0.) then
                sink = min (pgsub * min (1., dim (tz (k), t_sub) * gs_fac), qg (k))
                sub = sub + sink * dp (k)
            else
                sink = 0.
                if (tz (k) .le. tice) then
                    sink = max (pgsub, dq, (tz (k) - tice) / tcpk (k))
                endif
                dep = dep - sink * dp (k)
            endif
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                sink, 0., 0., 0., 0., - sink, te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pgdep_pgsub

! =======================================================================
! cloud fraction diagnostic
! =======================================================================

subroutine cloud_fraction (ks, ke, pz, den, qv, ql, qr, qi, qs, qg, qa, tz, h_var, gsize)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: h_var, gsize
    
    real, intent (in), dimension (ks:ke) :: pz, den
    
    real (kind = r8), intent (in), dimension (ks:ke) :: tz
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: q_plus, q_minus
    real :: rh, rqi, tin, qsw, qsi, qpz, qstar, sigma, gam
    real :: dqdt, dq, liq, ice
    real :: qa10, qa100
    
    real, dimension (ks:ke) :: q_liq, q_sol, q_cond, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), dimension (ks:ke) :: cvm, te8
    
    ! -----------------------------------------------------------------------
    ! calculate heat capacities and latent heat coefficients
    ! -----------------------------------------------------------------------
    
    call cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, cvm, te8, tz, &
        lcpk, icpk, tcpk, tcp3)
    
    do k = ks, ke
        
        ! combine water species
        
        ice = q_sol (k)
        q_sol (k) = qi (k)
        if (rad_snow) then
            q_sol (k) = qi (k) + qs (k)
            if (rad_graupel) then
                q_sol (k) = qi (k) + qs (k) + qg (k)
            endif
        endif
        
        liq = q_liq (k)
        q_liq (k) = ql (k)
        if (rad_rain) then
            q_liq (k) = ql (k) + qr (k)
        endif
        
        q_cond (k) = q_liq (k) + q_sol (k)
        qpz = qv (k) + q_cond (k)
        
        ! use the "liquid - frozen water temperature" (tin) to compute saturated specific humidity
        
        ice = ice - q_sol (k)
        liq = liq - q_liq (k)
        tin = (te8 (k) - lv00 * qpz + li00 * ice) / mhc (qpz, liq, ice)
        
        ! calculate saturated specific humidity
        
        if (tin .le. t_wfr) then
            qstar = iqs (tin, den (k), dqdt)
        elseif (tin .ge. tice) then
            qstar = wqs (tin, den (k), dqdt)
        else
            qsi = iqs (tin, den (k), dqdt)
            qsw = wqs (tin, den (k), dqdt)
            if (q_cond (k) .gt. qcmin) then
                rqi = q_sol (k) / q_cond (k)
            else
                rqi = (tice - tin) / (tice - t_wfr)
            endif
            qstar = rqi * qsi + (1. - rqi) * qsw
        endif
        
        ! cloud schemes
        
        rh = qpz / qstar
        
        if (cfflag .eq. 1) then
            if (rh .gt. rh_thres .and. qpz .gt. qcmin) then
                
                dq = h_var * qpz
                if (do_cld_adj) then
                    q_plus = qpz + dq * f_dq_p * min (1.0, max (0.0, (pz (k) - 200.e2) / &
                         (1000.e2 - 200.e2)))
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
                                qa (k) = (q_plus - qstar) / ((dq * f_dq_p + dq * f_dq_m) * &
                                     (1. - q_cond (k)))
                            endif
                        else
                            qa (k) = 0.
                        endif
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
                qa (k) = exp (xr_a * log (rh)) * (1.0 - exp (- xr_b * max (0.0, q_cond (k)) / &
                    max (1.e-5, exp (xr_c * log (max (1.e-10, 1.0 - rh) * qstar)))))
                qa (k) = max (0.0, min (1., qa (k)))
            else
                qa (k) = 0.0
            endif
        endif
        
        if (cfflag .eq. 3) then
            if (q_cond (k) .gt. qcmin) then
                qa (k) = 1. / 50. * (5.77 * (100. - gsize / 1000.) * &
                    exp (1.07 * log (max (qcmin * 1000., q_cond (k) * 1000.))) + &
                    4.82 * (gsize / 1000. - 50.) * &
                    exp (0.94 * log (max (qcmin * 1000., q_cond (k) * 1000.))))
                qa (k) = qa (k) * (0.92 / 0.96 * q_liq (k) / q_cond (k) + &
                    1.0 / 0.96 * q_sol (k) / q_cond (k))
                qa (k) = max (0.0, min (1., qa (k)))
            else
                qa (k) = 0.0
            endif
        endif
        
        if (cfflag .eq. 4) then
            sigma = 0.28 + exp (0.49 * log (max (qcmin * 1000., q_cond (k) * 1000.)))
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
    
end subroutine cloud_fraction

! =======================================================================
! piecewise parabolic lagrangian scheme
! this subroutine is the same as map1_q2 in fv_mapz_mod.
! =======================================================================

subroutine lagrangian_fall_ppm (ks, ke, zs, ze, zt, dp, q, precip, m1)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: zs
    
    real, intent (in), dimension (ks:ke + 1) :: ze, zt
    
    real, intent (in), dimension (ks:ke) :: dp
    
    real, intent (inout), dimension (ks:ke) :: q
    
    real, intent (inout) :: precip
    
    real, intent (out), dimension (ks:ke) :: m1
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k, k0, n, m
    
    real :: a4 (4, ks:ke), pl, pr, delz, esl
    
    real, parameter :: r3 = 1. / 3., r23 = 2. / 3.
    
    real, dimension (ks:ke) :: qm, dz
    
    ! -----------------------------------------------------------------------
    ! density:
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        dz (k) = zt (k) - zt (k + 1)
        q (k) = q (k) * dp (k)
        a4 (1, k) = q (k) / dz (k)
        qm (k) = 0.
    enddo
    
    ! -----------------------------------------------------------------------
    ! construct vertical profile with zt as coordinate
    ! -----------------------------------------------------------------------
    
    call cs_profile (a4 (1, ks), dz (ks), ke - ks + 1)
    
    k0 = ks
    do k = ks, ke
        do n = k0, ke
            if (ze (k) .le. zt (n) .and. ze (k) .ge. zt (n + 1)) then
                pl = (zt (n) - ze (k)) / dz (n)
                if (zt (n + 1) .le. ze (k + 1)) then
=======

    call cs_profile (a4 (1, ks), dz (ks), ke - ks + 1, mono)

    k0 = ks
    do k = ks, ke
        do n = k0, ke
            if (ze (k) <= zt (n) .and. ze (k) >= zt (n + 1)) then
                pl = (zt (n) - ze (k)) / dz (n)
                if (zt (n + 1) <= ze (k + 1)) then
>>>>>>> user/lnz/shield2022
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
    precip = precip + m1 (ke)
    
    ! -----------------------------------------------------------------------
    ! convert back to * dry * mixing ratio:
    ! dp must be dry air_mass (because moist air mass will be changed due to terminal fall) .
    ! -----------------------------------------------------------------------
    
    do k = ks, ke
        q (k) = qm (k) / dp (k)
    enddo
    
end subroutine lagrangian_fall_ppm

! =======================================================================
! vertical profile reconstruction
! this subroutine is the same as cs_profile in fv_mapz_mod where iv = 0 and kord = 9
! =======================================================================

subroutine cs_profile (a4, del, km)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: km
    
    real, intent (in) :: del (km)
    
    real, intent (inout) :: a4 (4, km)
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    logical :: extm (km)
    
    real :: gam (km), q (km + 1), d4, bet, a_bot, grat, pmp, lac
    real :: pmp_1, lac_1, pmp_2, lac_2, da1, da2, a6da
    
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
            ! apply large - scale constraints to all fields if not local max / min
            q (k) = min (q (k), max (a4 (1, k - 1), a4 (1, k)))
            q (k) = max (q (k), min (a4 (1, k - 1), a4 (1, k)))
        else
            if (gam (k - 1) .gt. 0.) then
                ! there exists a local max
                q (k) = max (q (k), min (a4 (1, k - 1), a4 (1, k)))
            else
                ! there exists a local min
                q (k) = min (q (k), max (a4 (1, k - 1), a4 (1, k)))
                ! positive-definite
                q (k) = max (q (k), 0.0)
            endif
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! bottom:
    ! -----------------------------------------------------------------------
    
    q (km) = min (q (km), max (a4 (1, km - 1), a4 (1, km)))
    q (km) = max (q (km), min (a4 (1, km - 1), a4 (1, km)), 0.)
    q (km + 1) = max (q (km + 1), 0.)
    
    do k = 1, km
        a4 (2, k) = q (k)
        a4 (3, k) = q (k + 1)
    enddo
    
    do k = 1, km
        if (k .eq. 1 .or. k .eq. km) then
            extm (k) = (a4 (2, k) - a4 (1, k)) * (a4 (3, k) - a4 (1, k)) .gt. 0.
        else
            extm (k) = gam (k) * gam (k + 1) .lt. 0.
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! apply constraints
    ! f (s) = al + s * [ (ar - al) + a6 * (1 - s) ] (0 <= s <= 1)
    ! always use monotonic mapping
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! top:
    ! -----------------------------------------------------------------------

    a4 (2, 1) = max (0., a4 (2, 1))
    
    ! -----------------------------------------------------------------------
    ! Huynh's 2nd constraint for interior:
    ! -----------------------------------------------------------------------

    do k = 3, km - 2
        if (extm (k)) then
            ! positive definite constraint only if true local extrema
            if (a4 (1, k) .lt. qcmin .or. extm (k - 1) .or. extm (k + 1)) then
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
    ! bottom:
    ! -----------------------------------------------------------------------
    
    a4 (2, km) = a4 (1, km)
    a4 (3, km) = a4 (1, km)
    a4 (4, km) = 0.
    
end subroutine cs_profile

! =======================================================================
! cubic spline (cs) limiters or boundary conditions
! a positive-definite constraint (iv = 0) is applied to tracers in every layer, 
! adjusting the top-most and bottom-most interface values to enforce positive.
! this subroutine is the same as cs_limiters in fv_mapz_mod where iv = 0.
! =======================================================================

subroutine cs_limiters (km, a4)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: km
    
    real, intent (inout) :: a4 (4, km) ! ppm array
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real, parameter :: r12 = 1. / 12.
    
    do k = 1, km
        if (a4 (1, k) .le. 0.) then
            a4 (2, k) = a4 (1, k)
            a4 (3, k) = a4 (1, k)
            a4 (4, k) = 0.
        else
            if (abs (a4 (3, k) - a4 (2, k)) .lt. - a4 (4, k)) then
                if ((a4 (1, k) + 0.25 * (a4 (3, k) - a4 (2, k)) ** 2 / a4 (4, k) + &
                    a4 (4, k) * r12) .lt. 0.) then
                    ! local minimum is negative
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
        endif
    enddo
    
end subroutine cs_limiters

! =======================================================================
! time-implicit monotonic scheme
! =======================================================================

subroutine implicit_fall (dts, ks, ke, ze, vt, dp, q, precip, m1)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke + 1) :: ze
    
    real, intent (in), dimension (ks:ke) :: vt, dp
    
    real, intent (inout), dimension (ks:ke) :: q
    
    real, intent (inout) :: precip
    
    real, intent (out), dimension (ks:ke) :: m1
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real, dimension (ks:ke) :: dz, qm, dd
    
    do k = ks, ke
        dz (k) = ze (k) - ze (k + 1)
        dd (k) = dts * vt (k)
        q (k) = q (k) * dp (k)
    enddo
    
    qm (ks) = q (ks) / (dz (ks) + dd (ks))
    do k = ks + 1, ke
        qm (k) = (q (k) + qm (k - 1) * dd (k - 1)) / (dz (k) + dd (k))
    enddo
    
    do k = ks, ke
        qm (k) = qm (k) * dz (k)
    enddo
    
    m1 (ks) = q (ks) - qm (ks)
    do k = ks + 1, ke
        m1 (k) = m1 (k - 1) + q (k) - qm (k)
    enddo
    precip = precip + m1 (ke)
    
    do k = ks, ke
        q (k) = qm (k) / dp (k)
    enddo
    
end subroutine implicit_fall

! =======================================================================
! time-explicit monotonic scheme
! =======================================================================

subroutine explicit_fall (dts, ks, ke, ze, vt, dp, q, precip, m1)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke + 1) :: ze
    
    real, intent (in), dimension (ks:ke) :: vt, dp
    
    real, intent (inout), dimension (ks:ke) :: q
    
    real, intent (inout) :: precip
    
    real, intent (out), dimension (ks:ke) :: m1
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: n, k, nstep
    
    real, dimension (ks:ke) :: dz, qm, q0, dd
    
    do k = ks, ke
        dz (k) = ze (k) - ze (k + 1)
        dd (k) = dts * vt (k)
        q0 (k) = q (k) * dp (k)
    enddo
    
    nstep = 1 + int (maxval (dd / dz))
    do k = ks, ke
        dd (k) = dd (k) / nstep
        q (k) = q0 (k)
    enddo
    
    do n = 1, nstep
        qm (ks) = q (ks) - q (ks) * dd (ks) / dz (ks)
        do k = ks + 1, ke
            qm (k) = q (k) - q (k) * dd (k) / dz (k) + q (k - 1) * dd (k - 1) / dz (k)
        enddo
        q = qm
    enddo
    
    m1 (ks) = q0 (ks) - qm (ks)
    do k = ks + 1, ke
        m1 (k) = m1 (k - 1) + q0 (k) - qm (k)
    enddo
    precip = precip + m1 (ke)
    
    do k = ks, ke
        q (k) = qm (k) / dp (k)
    enddo
    
end subroutine explicit_fall

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
! accretion function, Lin et al. (1983)
! =======================================================================

function acr2d (qden, c, denfac, blin, mu)
    
    implicit none
    
    real :: acr2d
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: qden, c, denfac, blin, mu
    
    acr2d = denfac * c * exp (1. / (mu + 3) * (2 + mu + blin) * log (6 * qden))
    
end function acr2d

! =======================================================================
! accretion function, Lin et al. (1983)
! =======================================================================

function acr3d (v1, v2, q1, q2, c, acco, acc1, acc2, den)
    
    implicit none
    
    real :: acr3d
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: v1, v2, c, den, q1, q2, acco (3), acc1, acc2
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i
    
    real :: t1, t2, tmp
    
    t1 = exp (1. / (acc1 + 3) * log (6 * q1 * den))
    t2 = exp (1. / (acc2 + 3) * log (6 * q2 * den))
    
    acr3d = c * abs (v1 - v2) / den
    
    tmp = 0
    do i = 1, 3
        tmp = tmp + acco (i) * exp ((6 + acc1 - i) * log (t1)) * exp ((acc2 + i - 1) * log (t2))
    enddo
    
    acr3d = acr3d * tmp
    
end function acr3d

! =======================================================================
! ventilation coefficient, Lin et al. (1983)
! =======================================================================

function vent_coeff (qden, c1, c2, denfac, blin, mu)
    
    implicit none
    
    real :: vent_coeff
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: qden, c1, c2, denfac, blin, mu
    
    vent_coeff = c1 + c2 * exp (1. / (mu + 3) * (3 + 2 * mu + blin) / 2 * log (6 * qden)) * &
        sqrt (denfac) / exp (1. / (mu + 3) * (1 + mu) * log (6 * qden))
    
end function vent_coeff

! =======================================================================
! sublimation or evaporation function, Lin et al. (1983)
! =======================================================================

function psub (t2, dq, qden, qsat, c, den, denfac, blin, mu, cpk, cvm)
    
    implicit none
    
    real :: psub
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: t2, dq, qden, qsat, c (5), den, denfac, blin, cpk, mu
    
    real (kind = r8), intent (in) :: cvm
    
    psub = c (1) * t2 * dq * exp (1. / (mu + 3) * (1 + mu) * log (6 * qden)) * &
        vent_coeff (qden, c (2), c (3), denfac, blin, mu) / &
         (c (4) * t2 + c (5) * (cpk * cvm) ** 2 * qsat * den)
    
end function psub

! =======================================================================
! melting function, Lin et al. (1983)
! =======================================================================

function pmlt (tc, dq, qden, pxacw, pxacr, c, den, denfac, blin, mu, lcpk, icpk, cvm)
    
    implicit none
    
    real :: pmlt
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tc, dq, qden, pxacw, pxacr, c (4), den, denfac, blin, lcpk, icpk, mu
    
    real (kind = r8), intent (in) :: cvm
    
    pmlt = (c (1) / (icpk * cvm) * tc / den - c (2) * lcpk / icpk * dq) * &
        exp (1. / (mu + 3) * (1 + mu) * log (6 * qden)) * &
        vent_coeff (qden, c (3), c (4), denfac, blin, mu) + &
        c_liq / (icpk * cvm) * tc * (pxacw + pxacr)
    
end function pmlt

! =======================================================================
! sedimentation of horizontal momentum
! =======================================================================

subroutine sedi_uv (ks, ke, m1, dp, u, v)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in), dimension (ks:ke) :: m1, dp
    
    real, intent (inout), dimension (ks:ke) :: u, v
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    do k = ks + 1, ke
        u (k) = (dp (k) * u (k) + m1 (k - 1) * u (k - 1)) / (dp (k) + m1 (k - 1))
        v (k) = (dp (k) * v (k) + m1 (k - 1) * v (k - 1)) / (dp (k) + m1 (k - 1))
    enddo
    
end subroutine sedi_uv

! =======================================================================
! sedimentation of vertical momentum
! =======================================================================

subroutine sedi_w (ks, ke, m1, w, vt, dm)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in), dimension (ks:ke) :: m1, vt, dm
    
    real, intent (inout), dimension (ks:ke) :: w
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    w (ks) = w (ks) + m1 (ks) * vt (ks) / dm (ks)
    do k = ks + 1, ke
        w (k) = (dm (k) * w (k) + m1 (k - 1) * (w (k - 1) - vt (k - 1)) + m1 (k) * vt (k)) / &
             (dm (k) + m1 (k - 1))
    enddo
    
end subroutine sedi_w

! =======================================================================
! sedimentation of heat
! =======================================================================

subroutine sedi_heat (ks, ke, dm, m1, dz, tz, qv, ql, qr, qi, qs, qg, cw)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: cw
    
    real, intent (in), dimension (ks:ke) :: dm, m1, dz, qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real, dimension (ks:ke) :: dgz, cv0
    
    do k = ks + 1, ke
        dgz (k) = - 0.5 * grav * (dz (k - 1) + dz (k))
        cv0 (k) = dm (k) * (cv_air + qv (k) * cv_vap + (qr (k) + ql (k)) * c_liq + &
             (qi (k) + qs (k) + qg (k)) * c_ice) + cw * (m1 (k) - m1 (k - 1))
    enddo
    
    do k = ks + 1, ke
        tz (k) = (cv0 (k) * tz (k) + m1 (k - 1) * (cw * tz (k - 1) + dgz (k))) / &
             (cv0 (k) + cw * m1 (k - 1))
    enddo
    
end subroutine sedi_heat

! =======================================================================
! fast saturation adjustments
! =======================================================================

subroutine fast_sat_adj (dtm, is, ie, ks, ke, hydrostatic, consv_te, &
        te, qv, ql, qr, qi, qs, qg, qa, qnl, qni, hs, delz, pt, delp, &
        q_con, cappa, gsize, last_step, condensation, evaporation, &
        deposition, sublimation, do_sat_adj)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: is, ie, ks, ke
    
    logical, intent (in) :: hydrostatic, last_step, consv_te, do_sat_adj
    
    real, intent (in) :: dtm
    
    real, intent (in), dimension (is:ie) :: hs, gsize
    
    real, intent (in), dimension (is:ie, ks:ke) :: delz, qnl, qni
    
    real, intent (inout), dimension (is:ie, ks:ke) :: delp, pt, te
    real, intent (inout), dimension (is:ie, ks:ke) :: qv, ql, qr, qi, qs, qg, qa
    
    real, intent (inout), dimension (is:, ks:) :: q_con, cappa
    
    real, intent (inout), dimension (is:ie) :: condensation, deposition
    real, intent (inout), dimension (is:ie) :: evaporation, sublimation
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real, dimension (is:ie, ks:ke) :: ua, va, wa, prefluxw, prefluxr, prefluxi, prefluxs, prefluxg
    
    real, dimension (is:ie) :: water, rain, ice, snow, graupel
    
    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
    
    ua = 0.0
    va = 0.0
    wa = 0.0
    
    water = 0.0
    rain = 0.0
    ice = 0.0
    snow = 0.0
    graupel = 0.0
    
    prefluxw = 0.0
    prefluxr = 0.0
    prefluxi = 0.0
    prefluxs = 0.0
    prefluxg = 0.0
    
    ! -----------------------------------------------------------------------
    ! define various heat capacities and latent heat coefficients at 0 deg K
    ! -----------------------------------------------------------------------
    
    call setup_mhc_lhc (hydrostatic)
    
    ! -----------------------------------------------------------------------
    ! major cloud microphysics driver
    ! -----------------------------------------------------------------------
    
    call mpdrv (hydrostatic, ua, va, wa, delp, pt, qv, ql, qr, qi, qs, qg, qa, &
        qnl, qni, delz, is, ie, ks, ke, dtm, water, rain, ice, snow, graupel, &
        gsize, hs, q_con, cappa, consv_te, te, prefluxw, prefluxr, prefluxi, &
        prefluxs, prefluxg, condensation, deposition, evaporation, sublimation, &
        last_step, .true., do_sat_adj, .false.)
    
end subroutine fast_sat_adj

! =======================================================================
! rain freezing to form graupel, simple version
! =======================================================================

subroutine pgfr_simp (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, &
        lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, sink, fac_r2g
    
    fac_r2g = 1. - exp (- dts / tau_r2g)
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .lt. 0. .and. qr (k) .gt. qcmin) then
            
            sink = (- tc * 0.025) ** 2 * qr (k)
            sink = min (qr (k), sink, - fac_r2g * tc / icpk (k))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., - sink, 0., 0., sink, te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine pgfr_simp

! =======================================================================
! snow melting to form cloud water and rain, simple version
! =======================================================================

subroutine psmlt_simp (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, cvm, te8, &
        lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real (kind = r8), intent (in), dimension (ks:ke) :: te8
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (inout), dimension (ks:ke) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: cvm, tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, tmp, sink, fac_smlt
    
    fac_smlt = 1. - exp (- dts / tau_smlt)
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        if (tc .ge. 0. .and. qs (k) .gt. qcmin) then
            
            sink = (tc * 0.1) ** 2 * qs (k)
            sink = min (qs (k), sink, fac_smlt * tc / icpk (k))
            tmp = min (sink, dim (qs_mlt, ql (k)))
            
            call update_qt (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., tmp, sink - tmp, 0., - sink, 0., te8 (k), cvm (k), tz (k), &
                lcpk (k), icpk (k), tcpk (k), tcp3 (k))
            
        endif
        
    enddo
    
end subroutine psmlt_simp

! =======================================================================
! cloud water to rain autoconversion, simple version
! =======================================================================

subroutine praut_simp (ks, ke, dts, tz, qv, ql, qr, qi, qs, qg)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, sink, fac_l2r
    
    fac_l2r = 1. - exp (- dts / tau_l2r)
    
    do k = ks, ke
        
        tc = tz (k) - t_wfr
        
        if (tc .gt. 0 .and. ql (k) .gt. ql0_max) then
            
            sink = fac_l2r * (ql (k) - ql0_max)
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., - sink, sink, 0., 0., 0.)
            
        endif
        
    enddo

end subroutine praut_simp
    
! =======================================================================
! cloud ice to snow autoconversion, simple version
! =======================================================================

subroutine psaut_simp (ks, ke, dts, qv, ql, qr, qi, qs, qg, tz, den)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in) :: dts
    
    real, intent (in), dimension (ks:ke) :: den
    
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (inout), dimension (ks:ke) :: tz
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: tc, sink, fac_i2s, qim
    
    fac_i2s = 1. - exp (- dts / tau_i2s)
    
    do k = ks, ke
        
        tc = tz (k) - tice
        
        qim = qi0_max / den (k)
        
        if (tc .lt. 0. .and. qi (k) .gt. qim) then
            
            sink = fac_i2s * (qi (k) - qim)
            
            call update_qq (qv (k), ql (k), qr (k), qi (k), qs (k), qg (k), &
                0., 0., 0., - sink, sink, 0.)
            
        endif
        
    enddo
    
end subroutine psaut_simp

! =======================================================================
! cloud radii diagnosis built for gfdl cloud microphysics
! =======================================================================

subroutine cld_eff_rad (is, ie, ks, ke, lsm, p, delp, t, qw, qi, qr, qs, qg, qa, &
        qcw, qci, qcr, qcs, qcg, rew, rei, rer, res, reg, cld, cloud, snowd, &
        cnvw, cnvi, cnvc)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: is, ie, ks, ke
    
    real, intent (in), dimension (is:ie) :: lsm, snowd
    
    real, intent (in), dimension (is:ie, ks:ke) :: delp, t, p, cloud
    real, intent (in), dimension (is:ie, ks:ke) :: qw, qi, qr, qs, qg, qa
    
    real, intent (in), dimension (is:ie, ks:ke), optional :: cnvw, cnvi, cnvc
    
    real, intent (inout), dimension (is:ie, ks:ke) :: qcw, qci, qcr, qcs, qcg
    real, intent (inout), dimension (is:ie, ks:ke) :: rew, rei, rer, res, reg
    real, intent (inout), dimension (is:ie, ks:ke) :: cld
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i, k, ind
    
    real, dimension (is:ie, ks:ke) :: qmw, qmr, qmi, qms, qmg
    
    real :: dpg, rho, ccnw, mask, cor, tc, bw
    real :: lambdaw, lambdar, lambdai, lambdas, lambdag, rei_fac
    
    real :: ccno = 90.
    real :: ccnl = 270.
    
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
    
    ! -----------------------------------------------------------------------
    ! merge convective cloud to total cloud
    ! -----------------------------------------------------------------------
    
    if (present (cnvw)) then
        qmw = qmw + cnvw
    endif
    if (present (cnvi)) then
        qmi = qmi + cnvi
    endif
    if (present (cnvc)) then
        cld = cnvc + (1 - cnvc) * cld
    endif
    
    ! -----------------------------------------------------------------------
    ! combine liquid and solid phases
    ! -----------------------------------------------------------------------
    
    if (liq_ice_combine) then
        do i = is, ie
            do k = ks, ke
                qmw (i, k) = qmw (i, k) + qmr (i, k)
                qmr (i, k) = 0.0
                qmi (i, k) = qmi (i, k) + qms (i, k) + qmg (i, k)
                qms (i, k) = 0.0
                qmg (i, k) = 0.0
            enddo
        enddo
    endif
    
    ! -----------------------------------------------------------------------
    ! combine snow and graupel
    ! -----------------------------------------------------------------------
    
    if (snow_grauple_combine) then
        do i = is, ie
            do k = ks, ke
                qms (i, k) = qms (i, k) + qmg (i, k)
                qmg (i, k) = 0.0
            enddo
        enddo
    endif
    
        
    do i = is, ie

        do k = ks, ke
            
            qmw (i, k) = max (qmw (i, k), 0.0)
            qmi (i, k) = max (qmi (i, k), 0.0)
            qmr (i, k) = max (qmr (i, k), 0.0)
            qms (i, k) = max (qms (i, k), 0.0)
            qmg (i, k) = max (qmg (i, k), 0.0)
            
            cld (i, k) = min (max (cld (i, k), 0.0), 1.0)
            
            mask = min (max (lsm (i), 0.0), 2.0)
            
            dpg = abs (delp (i, k)) / grav
            rho = p (i, k) / (rdgas * t (i, k))
            
            tc = t (i, k) - tice
            
            if (rewflag .eq. 1) then
                
                ! -----------------------------------------------------------------------
                ! cloud water (Martin et al. 1994)
                ! -----------------------------------------------------------------------
                
                if (prog_ccn) then
                    ! boucher and lohmann (1995)
                    ccnw = (1.0 - abs (mask - 1.0)) * &
                         (10. ** 2.24 * (qa (i, k) * rho * 1.e9) ** 0.257) + &
                        abs (mask - 1.0) * &
                         (10. ** 2.06 * (qa (i, k) * rho * 1.e9) ** 0.48)
                else
                    ccnw = ccno * abs (mask - 1.0) + ccnl * (1.0 - abs (mask - 1.0))
                endif
                
                if (qmw (i, k) .gt. qcmin) then
                    qcw (i, k) = dpg * qmw (i, k) * 1.0e3
                    rew (i, k) = exp (1.0 / 3.0 * log ((3.0 * qmw (i, k) * rho) / &
                         (4.0 * pi * rhow * ccnw))) * 1.0e4
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif
                
            endif
            
            if (rewflag .eq. 2) then
                
                ! -----------------------------------------------------------------------
                ! cloud water (Martin et al. 1994, gfdl revision)
                ! -----------------------------------------------------------------------
                
                if (prog_ccn) then
                    ! boucher and lohmann (1995)
                    ccnw = (1.0 - abs (mask - 1.0)) * &
                         (10. ** 2.24 * (qa (i, k) * rho * 1.e9) ** 0.257) + &
                        abs (mask - 1.0) * &
                         (10. ** 2.06 * (qa (i, k) * rho * 1.e9) ** 0.48)
                else
                    ccnw = 1.077 * ccno * abs (mask - 1.0) + 1.143 * ccnl * (1.0 - abs (mask - 1.0))
                endif
                
                if (qmw (i, k) .gt. qcmin) then
                    qcw (i, k) = dpg * qmw (i, k) * 1.0e3
                    rew (i, k) = exp (1.0 / 3.0 * log ((3.0 * qmw (i, k) * rho) / &
                         (4.0 * pi * rhow * ccnw))) * 1.0e4
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif
                
            endif
            
            if (rewflag .eq. 3) then
                
                ! -----------------------------------------------------------------------
                ! cloud water (Kiehl et al. 1994)
                ! -----------------------------------------------------------------------
                
                if (qmw (i, k) .gt. qcmin) then
                    qcw (i, k) = dpg * qmw (i, k) * 1.0e3
                    rew (i, k) = 14.0 * abs (mask - 1.0) + &
                         (8.0 + (14.0 - 8.0) * min (1.0, max (0.0, - tc / 30.0))) * &
                         (1.0 - abs (mask - 1.0))
                    rew (i, k) = rew (i, k) + (14.0 - rew (i, k)) * &
                        min (1.0, max (0.0, snowd (i) / 1000.0))
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif
                
            endif
            
            if (rewflag .eq. 4) then
                
                ! -----------------------------------------------------------------------
                ! cloud water (Lin et al. 1983)
                ! -----------------------------------------------------------------------
                
                if (qmw (i, k) .gt. qcmin) then
                    qcw (i, k) = dpg * qmw (i, k) * 1.0e3
                    lambdaw = exp (1. / (muw + 3) * log (normw / (6 * qmw (i, k) * rho))) * expow
                    rew (i, k) = 0.5 * (muw + 2) / lambdaw * 1.0e6
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif
                
            endif
            
            if (rewflag .eq. 5) then
                
                ! -----------------------------------------------------------------------
                ! cloud water (Lin et al. 1983)
                ! -----------------------------------------------------------------------
                
                if (qmw (i, k) .gt. qcmin) then
                    qcw (i, k) = dpg * qmw (i, k) * 1.0e3
                    lambdaw = exp (1. / (muw + 3) * log (normw / (6 * qmw (i, k) * rho))) * expow
                    rew (i, k) = 0.5 * exp (log (gamma (4 + blinw) / 6) / blinw) / lambdaw * 1.0e6
                    rew (i, k) = max (rewmin, min (rewmax, rew (i, k)))
                else
                    qcw (i, k) = 0.0
                    rew (i, k) = rewmin
                endif
                
            endif
            
            if (reiflag .eq. 1) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (Heymsfield and Mcfarquhar 1996)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    rei_fac = log (1.0e3 * qmi (i, k) * rho)
                    if (tc .lt. - 50) then
                        rei (i, k) = beta / 9.917 * exp (0.109 * rei_fac) * 1.0e3
                    elseif (tc .lt. - 40) then
                        rei (i, k) = beta / 9.337 * exp (0.080 * rei_fac) * 1.0e3
                    elseif (tc .lt. - 30) then
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
                ! cloud ice (Donner et al. 1997)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    if (tc .le. - 55) then
                        rei (i, k) = 15.41627
                    elseif (tc .le. - 50) then
                        rei (i, k) = 16.60895
                    elseif (tc .le. - 45) then
                        rei (i, k) = 32.89967
                    elseif (tc .le. - 40) then
                        rei (i, k) = 35.29989
                    elseif (tc .le. - 35) then
                        rei (i, k) = 55.65818
                    elseif (tc .le. - 30) then
                        rei (i, k) = 85.19071
                    elseif (tc .le. - 25) then
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
                ! cloud ice (Fu 2007)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    rei (i, k) = 47.05 + tc * (0.6624 + 0.001741 * tc)
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (reiflag .eq. 4) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (Kristjansson et al. 2000)
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
                ! cloud ice (Wyser 1998)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    bw = - 2. + 1.e-3 * log10 (rho * qmi (i, k) / 50.e-3) * &
                        exp (1.5 * log (max (1.e-10, - tc)))
                    rei (i, k) = 377.4 + bw * (203.3 + bw * (37.91 + 2.3696 * bw))
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (reiflag .eq. 6) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (Sun and Rikus 1999, Sun 2001)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    rei_fac = log (1.0e3 * qmi (i, k) * rho)
                    rei (i, k) = 45.8966 * exp (0.2214 * rei_fac) + &
                        0.7957 * exp (0.2535 * rei_fac) * (tc + 190.0)
                    rei (i, k) = (1.2351 + 0.0105 * tc) * rei (i, k)
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (reiflag .eq. 7) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (Lin et al. 1983)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    lambdai = exp (1. / (mui + 3) * log (normi / (6 * qmi (i, k) * rho))) * expoi
                    rei (i, k) = 0.5 * (mui + 2) / lambdai * 1.0e6
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (reiflag .eq. 8) then
                
                ! -----------------------------------------------------------------------
                ! cloud ice (Lin et al. 1983)
                ! -----------------------------------------------------------------------
                
                if (qmi (i, k) .gt. qcmin) then
                    qci (i, k) = dpg * qmi (i, k) * 1.0e3
                    lambdai = exp (1. / (mui + 3) * log (normi / (6 * qmi (i, k) * rho))) * expoi
                    rei (i, k) = 0.5 * exp (log (gamma (4 + blini) / 6) / blini) / lambdai * 1.0e6
                    rei (i, k) = max (reimin, min (reimax, rei (i, k)))
                else
                    qci (i, k) = 0.0
                    rei (i, k) = reimin
                endif
                
            endif
            
            if (rerflag .eq. 1) then
                
                ! -----------------------------------------------------------------------
                ! rain (Lin et al. 1983)
                ! -----------------------------------------------------------------------
                
                if (qmr (i, k) .gt. qcmin) then
                    qcr (i, k) = dpg * qmr (i, k) * 1.0e3
                    lambdar = exp (1. / (mur + 3) * log (normr / (6 * qmr (i, k) * rho))) * expor
                    rer (i, k) = 0.5 * (mur + 2) / lambdar * 1.0e6
                    rer (i, k) = max (rermin, min (rermax, rer (i, k)))
                else
                    qcr (i, k) = 0.0
                    rer (i, k) = rermin
                endif
                
            endif
            
            if (rerflag .eq. 2) then
                
                ! -----------------------------------------------------------------------
                ! rain (Lin et al. 1983)
                ! -----------------------------------------------------------------------
                
                if (qmr (i, k) .gt. qcmin) then
                    qcr (i, k) = dpg * qmr (i, k) * 1.0e3
                    lambdar = exp (1. / (mur + 3) * log (normr / (6 * qmr (i, k) * rho))) * expor
                    rer (i, k) = 0.5 * exp (log (gamma (4 + blinr) / 6) / blinr) / lambdar * 1.0e6
                    rer (i, k) = max (rermin, min (rermax, rer (i, k)))
                else
                    qcr (i, k) = 0.0
                    rer (i, k) = rermin
                endif
                
            endif
            
            if (resflag .eq. 1) then
                
                ! -----------------------------------------------------------------------
                ! snow (Lin et al. 1983)
                ! -----------------------------------------------------------------------
                
                if (qms (i, k) .gt. qcmin) then
                    qcs (i, k) = dpg * qms (i, k) * 1.0e3
                    lambdas = exp (1. / (mus + 3) * log (norms / (6 * qms (i, k) * rho))) * expos
                    res (i, k) = 0.5 * (mus + 2) / lambdas * 1.0e6
                    res (i, k) = max (resmin, min (resmax, res (i, k)))
                else
                    qcs (i, k) = 0.0
                    res (i, k) = resmin
                endif
                
            endif
            
            if (resflag .eq. 2) then
                
                ! -----------------------------------------------------------------------
                ! snow (Lin et al. 1983)
                ! -----------------------------------------------------------------------
                
                if (qms (i, k) .gt. qcmin) then
                    qcs (i, k) = dpg * qms (i, k) * 1.0e3
                    lambdas = exp (1. / (mus + 3) * log (norms / (6 * qms (i, k) * rho))) * expos
                    res (i, k) = 0.5 * exp (log (gamma (4 + blins) / 6) / blins) / lambdas * 1.0e6
                    res (i, k) = max (resmin, min (resmax, res (i, k)))
                else
                    qcs (i, k) = 0.0
                    res (i, k) = resmin
                endif
                
            endif
            
            if (regflag .eq. 1) then
                
                ! -----------------------------------------------------------------------
                ! graupel (Lin et al. 1983)
                ! -----------------------------------------------------------------------
                
                if (qmg (i, k) .gt. qcmin) then
                    qcg (i, k) = dpg * qmg (i, k) * 1.0e3
                    if (do_hail) then
                        lambdag = exp (1. / (muh + 3) * log (normh / (6 * qmg (i, k) * rho))) * expoh
                        reg (i, k) = 0.5 * (muh + 2) / lambdag * 1.0e6
                    else
                        lambdag = exp (1. / (mug + 3) * log (normg / (6 * qmg (i, k) * rho))) * expog
                        reg (i, k) = 0.5 * (mug + 2) / lambdag * 1.0e6
                    endif
                    reg (i, k) = max (regmin, min (regmax, reg (i, k)))
                else
                    qcg (i, k) = 0.0
                    reg (i, k) = regmin
                endif
                
            endif
            
            if (regflag .eq. 2) then
                
                ! -----------------------------------------------------------------------
                ! graupel (Lin et al. 1983)
                ! -----------------------------------------------------------------------
                
                if (qmg (i, k) .gt. qcmin) then
                    qcg (i, k) = dpg * qmg (i, k) * 1.0e3
                    if (do_hail) then
                        lambdag = exp (1. / (muh + 3) * log (normh / (6 * qmg (i, k) * rho))) * expoh
                        reg (i, k) = 0.5 * exp (log (gamma (4 + blinh) / 6) / blinh) / lambdag * 1.0e6
                    else
                        lambdag = exp (1. / (mug + 3) * log (normg / (6 * qmg (i, k) * rho))) * expog
                        reg (i, k) = 0.5 * exp (log (gamma (4 + bling) / 6) / bling) / lambdag * 1.0e6
                    endif
                    reg (i, k) = max (regmin, min (regmax, reg (i, k)))
                else
                    qcg (i, k) = 0.0
                    reg (i, k) = regmin
                endif
                
            endif
            
        enddo
        
    enddo
<<<<<<< HEAD
    
end subroutine cld_eff_rad

! =======================================================================
! radar reflectivity
! =======================================================================

subroutine rad_ref (is, ie, js, je, isd, ied, jsd, jed, q, pt, delp, peln, &
        delz, dbz, maxdbz, allmax, npz, ncnst, hydrostatic, zvir, &
        do_inline_mp, sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, mp_top)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    logical, intent (in) :: hydrostatic, do_inline_mp
    
    integer, intent (in) :: is, ie, js, je, isd, ied, jsd, jed
    integer, intent (in) :: npz, ncnst, mp_top
    integer, intent (in) :: sphum, liq_wat, ice_wat, rainwat, snowwat, graupel
    
    real, intent (in) :: zvir
    
    real, intent (in), dimension (is:, js:, 1:) :: delz
    
    real, intent (in), dimension (isd:ied, jsd:jed, npz) :: pt, delp
    
    real, intent (in), dimension (isd:ied, jsd:jed, npz, ncnst) :: q
    
    real, intent (in), dimension (is:ie, npz + 1, js:je) :: peln
    
    real, intent (out) :: allmax
    
    real, intent (out), dimension (is:ie, js:je) :: maxdbz
    
    real, intent (out), dimension (is:ie, js:je, npz) :: dbz
    
=======

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

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    integer :: i, j, k
    
    real, parameter :: alpha = 0.224, mp_const = 200 * exp (1.6 * log (3.6e6))
    
    real (kind = r8) :: qden, z_e, fac_r, fac_s, fac_g, fac_sw, fac_sd
    
    real, dimension (npz) :: den, denfac, qmr, qms, qmg, vtr, vts, vtg
    
=======

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

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! return if the microphysics scheme doesn't include rain
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    if (rainwat .lt. 1) return
    
=======

    do k = 3, km - 1
        if (gam (k - 1) * gam (k + 1) > 0.) then
            q (k) = min (q (k), max (a4 (1, k - 1), a4 (1, k)))
            q (k) = max (q (k), min (a4 (1, k - 1), a4 (1, k)))
        else
            if (gam (k - 1) > 0.) then
                ! there exists a local max
                q (k) = max (q (k), min (a4 (1, k - 1), a4 (1, k)))
            else
                ! there exists a local min
                q (k) = min (q (k), max (a4 (1, k - 1), a4 (1, k)))
                q (k) = max (q (k), 0.0)
            endif
        endif
    enddo

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! initialization
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    dbz = - 20.
    maxdbz = - 20.
    allmax = - 20.
    
=======

    q (km) = min (q (km), max (a4 (1, km - 1), a4 (1, km)))
    q (km) = max (q (km), min (a4 (1, km - 1), a4 (1, km)), 0.)
    ! q (km + 1) = max (q (km + 1), 0.)

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! constants for radar reflectivity
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    if (radr_flag .eq. 1 .or. radr_flag .eq. 2) &
        fac_r = gamma (mur + 6) * exp (- 1. / (mur + 3) * (mur + 6) * log (normr)) * n0r_sig * &
            exp (- 3.0 * log (expor)) * 1.e18
    
    if (rads_flag .eq. 1) &
        fac_s = gamma (mus + 6) * exp (- 1. / (mus + 3) * (mus + 6) * log (norms)) * n0s_sig * &
            exp (- 3.0 * log (expos)) * 1.e18 * alpha * (rhos / rhor) ** 2
    if (rads_flag .eq. 2) then
        fac_sw = gamma (mus + 6) * exp (- 1. / (mus + 3) * (mus + 6) * log (norms)) * n0s_sig * &
            exp (- 3.0 * log (expos)) * 1.e18
        fac_sd = gamma (mus + 6) * exp (- 1. / (mus + 3) * (mus + 6) * log (norms)) * n0s_sig * &
            exp (- 3.0 * log (expos)) * 1.e18 * alpha * (rhos / rhoi) ** 2
    endif
    
    if (radg_flag .eq. 1) then
        if (do_hail .and. .not. do_inline_mp) then
            fac_g = gamma (muh + 6) * exp (- 1. / (muh + 3) * (muh + 6) * log (normh)) * n0h_sig * &
                exp (- 3.0 * log (expoh)) * 1.e18 * alpha * (rhoh / rhor) ** 2
=======

    do k = 1, km - 1
        a4 (2, k) = q (k)
        a4 (3, k) = q (k + 1)
    enddo

    do k = 2, km - 1
        if (gam (k) * gam (k + 1) > 0.0) then
            extm (k) = .false.
>>>>>>> user/lnz/shield2022
        else
            fac_g = gamma (mug + 6) * exp (- 1. / (mug + 3) * (mug + 6) * log (normg)) * n0g_sig * &
                exp (- 3.0 * log (expog)) * 1.e18 * alpha * (rhog / rhor) ** 2
        endif
<<<<<<< HEAD
    endif
    if (radg_flag .eq. 2) then
        if (do_hail .and. .not. do_inline_mp) then
            fac_g = gamma (muh + 6) * exp (- 1. / (muh + 3) * (muh + 6) * log (normh)) * n0h_sig * &
                exp (- 3.0 * log (expoh)) * 1.e18
        else
            fac_g = gamma (mug + 6) * exp (- 1. / (mug + 3) * (mug + 6) * log (normg)) * n0g_sig * &
                exp (- 3.0 * log (expog)) * 1.e18
        endif
    endif
    
=======
    enddo

    if (do_mono) then
        do k = 3, km - 2
            if (extm (k)) then
                ! positive definite constraint only if true local extrema
                if (a4 (1, k) < qp_min .or. extm (k - 1) .or. extm (k + 1)) then
                    a4 (2, k) = a4 (1, k)
                    a4 (3, k) = a4 (1, k)
                endif
            else
                a4 (4, k) = 6. * a4 (1, k) - 3. * (a4 (2, k) + a4 (3, k))
                if (abs (a4 (4, k)) > abs (a4 (2, k) - a4 (3, k))) then
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
                if (a4 (1, k) < qp_min .or. extm (k - 1) .or. extm (k + 1)) then
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
        if (a6da < - da2) then
            a4 (4, k) = 3. * (a4 (2, k) - a4 (1, k))
            a4 (3, k) = a4 (2, k) - a4 (4, k)
        elseif (a6da > da2) then
            a4 (4, k) = 3. * (a4 (3, k) - a4 (1, k))
            a4 (2, k) = a4 (3, k) - a4 (4, k)
        endif
    endif

    call cs_limiters (km - 1, a4)

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! calculate radar reflectivity
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    do j = js, je
        do i = is, ie
            
            ! -----------------------------------------------------------------------
            ! air density
            ! -----------------------------------------------------------------------
            
            do k = 1, npz
                if (hydrostatic) then
                    den (k) = delp (i, j, k) / ((peln (i, k + 1, j) - peln (i, k, j)) * &
                        rdgas * pt (i, j, k) * (1. + zvir * q (i, j, k, sphum)))
                else
                    den (k) = - delp (i, j, k) / (grav * delz (i, j, k))
                endif
                qmr (k) = max (qcmin, q (i, j, k, rainwat))
                qms (k) = max (qcmin, q (i, j, k, snowwat))
                qmg (k) = max (qcmin, q (i, j, k, graupel))
            enddo
            
            do k = 1, npz
                denfac (k) = sqrt (den (npz) / den (k))
            enddo
            
            ! -----------------------------------------------------------------------
            ! fall speed
            ! -----------------------------------------------------------------------
            
            if (radr_flag .eq. 3) then
                call term_rsg (1, npz, qmr, den, denfac, vr_fac, vconr, blinr, normr, expor, &
                    mur, vr_max, const_vr, vtr)
                vtr = vtr / rhor
            endif
            
            if (rads_flag .eq. 3) then
                call term_rsg (1, npz, qms, den, denfac, vs_fac, vcons, blins, norms, expos, &
                    mus, vs_max, const_vs, vts)
                vts = vts / rhos
            endif
            
            if (radg_flag .eq. 3) then
                if (do_hail .and. .not. do_inline_mp) then
                    call term_rsg (1, npz, qmg, den, denfac, vg_fac, vconh, blinh, normh, expoh, &
                        muh, vg_max, const_vg, vtg)
                    vtg = vtg / rhoh
                else
                    call term_rsg (1, npz, qmg, den, denfac, vg_fac, vcong, bling, normg, expog, &
                        mug, vg_max, const_vg, vtg)
                    vtg = vtg / rhog
                endif
            endif
            
            ! -----------------------------------------------------------------------
            ! radar reflectivity
            ! -----------------------------------------------------------------------
            
            do k = mp_top + 1, npz
                z_e = 0.
                
                if (rainwat .gt. 0) then
                    qden = den (k) * qmr (k)
                    if (radr_flag .eq. 1 .or. radr_flag .eq. 2) then
                        z_e = z_e + fac_r * exp (1. / (mur + 3) * (mur + 6) * log (6 * qden))
                    endif
                    if (radr_flag .eq. 3) then
                        z_e = z_e + mp_const * exp (1.6 * log (qden * vtr (k)))
                    endif
                endif
                
                if (snowwat .gt. 0) then
                    qden = den (k) * qms (k)
                    if (rads_flag .eq. 1) then
                        if (pt (i, j, k) .lt. tice) then
                            z_e = z_e + fac_s * exp (1. / (mus + 3) * (mus + 6) * log (6 * qden))
                        else
                            z_e = z_e + fac_s * exp (1. / (mus + 3) * (mus + 6) * log (6 * qden)) / alpha
                        endif
                    endif
                    if (rads_flag .eq. 2) then
                        if (pt (i, j, k) .lt. tice) then
                            z_e = z_e + fac_sd * exp (1. / (mus + 3) * (mus + 6) * log (6 * qden))
                        else
                            z_e = z_e + fac_sw * exp (1. / (mus + 3) * (mus + 6) * log (6 * qden))
                        endif
                    endif
                    if (rads_flag .eq. 3) then
                        z_e = z_e + mp_const * exp (1.6 * log (qden * vts (k)))
                    endif
                endif
                
                if (graupel .gt. 0) then
                    qden = den (k) * qmg (k)
                    if (do_hail .and. .not. do_inline_mp) then
                        if (radg_flag .eq. 1) then
                            if (pt (i, j, k) .lt. tice) then
                                z_e = z_e + fac_g * exp (1. / (muh + 3) * (muh + 6) * log (6 * qden))
                            else
                                z_e = z_e + fac_g * exp (1. / (muh + 3) * (muh + 6) * log (6 * qden)) / alpha
                            endif
                        endif
                        if (radg_flag .eq. 2) then
                            z_e = z_e + fac_g * exp (1. / (muh + 3) * (muh + 6) * log (6 * qden))
                        endif
                    else
                        if (radg_flag .eq. 1) then
                            if (pt (i, j, k) .lt. tice) then
                                z_e = z_e + fac_g * exp (1. / (mug + 3) * (mug + 6) * log (6 * qden))
                            else
                                z_e = z_e + fac_g * exp (1. / (mug + 3) * (mug + 6) * log (6 * qden)) / alpha
                            endif
                        endif
                        if (radg_flag .eq. 2) then
                            z_e = z_e + fac_g * exp (1. / (mug + 3) * (mug + 6) * log (6 * qden))
                        endif
                    endif
                    if (radg_flag .eq. 3) then
                        z_e = z_e + mp_const * exp (1.6 * log (qden * vtg (k)))
                    endif
                endif
                
                dbz (i, j, k) = 10. * log10 (max (0.01, z_e))
            enddo
            
            do k = mp_top + 1, npz
                maxdbz (i, j) = max (dbz (i, j, k), maxdbz (i, j))
            enddo
            
            allmax = max (maxdbz (i, j), allmax)
            
        enddo
    enddo
    
end subroutine rad_ref

! =======================================================================
! moist heat capacity, 3 input variables
! =======================================================================

function mhc3 (qv, q_liq, q_sol)
    
    implicit none
    
    real (kind = r8) :: mhc3
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: qv, q_liq, q_sol
    
=======

    a4 (2, km) = a4 (1, km)
    a4 (3, km) = a4 (1, km)
    a4 (4, km) = 0.

end subroutine cs_profile

subroutine cs_limiters (km, a4)

    implicit none

    integer, intent (in) :: km

    real, intent (inout) :: a4 (4, km) ! ppm array

    real, parameter :: r12 = 1. / 12.

    integer :: k

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    mhc3 = one_r8 + qv * c1_vap + q_liq * c1_liq + q_sol * c1_ice
    
end function mhc3
=======

    do k = 1, km
        if (abs (a4 (3, k) - a4 (2, k)) < - a4 (4, k)) then
            if ((a4 (1, k) + 0.25 * (a4 (3, k) - a4 (2, k)) ** 2 / a4 (4, k) + a4 (4, k) * r12) < 0.) then
                if (a4 (1, k) < a4 (3, k) .and. a4 (1, k) < a4 (2, k)) then
                    a4 (3, k) = a4 (1, k)
                    a4 (2, k) = a4 (1, k)
                    a4 (4, k) = 0.
                elseif (a4 (3, k) > a4 (2, k)) then
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
>>>>>>> user/lnz/shield2022

! =======================================================================
! moist heat capacity, 4 input variables
! =======================================================================

<<<<<<< HEAD
function mhc4 (qd, qv, q_liq, q_sol)
    
    implicit none
    
    real (kind = r8) :: mhc4
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: qv, q_liq, q_sol
    
    real (kind = r8), intent (in) :: qd
    
=======
subroutine fall_speed (ks, ke, den, qs, qi, qg, ql, tk, vts, vti, vtg)

    implicit none

    integer, intent (in) :: ks, ke

    real (kind = r_grid), intent (in), dimension (ks:ke) :: tk
    real, intent (in), dimension (ks:ke) :: den, qs, qi, qg, ql
    real, intent (out), dimension (ks:ke) :: vts, vti, vtg

    ! fall velocity constants:

    real, parameter :: thi = 1.0e-8 ! cloud ice threshold for terminal fall
    real, parameter :: thg = 1.0e-8
    real, parameter :: ths = 1.0e-8

    real, parameter :: aa = - 4.14122e-5
    real, parameter :: bb = - 0.00538922
    real, parameter :: cc = - 0.0516344
    real, parameter :: dd = 0.00216078
    real, parameter :: ee = 1.9714

    ! marshall - palmer constants

    real, parameter :: vcons = 6.6280504
    real, parameter :: vcong = 87.2382675
    real, parameter :: vconh = vcong * sqrt (rhoh / rhog) ! 132.087495104005
    real, parameter :: norms = 942477796.076938
    real, parameter :: normg = 5026548245.74367
    real, parameter :: normh = pi * rhoh * rnzh ! 115233618.533674

    real, dimension (ks:ke) :: qden, tc, rhof

    real :: vi0

    integer :: k

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    mhc4 = qd + qv * c1_vap + q_liq * c1_liq + q_sol * c1_ice
    
end function mhc4

! =======================================================================
! moist heat capacity, 6 input variables
! =======================================================================

function mhc6 (qv, ql, qr, qi, qs, qg)
    
    implicit none
    
    real (kind = r8) :: mhc6
    
=======

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
            if (qi (k) < thi) then ! this is needed as the fall - speed maybe problematic for small qi
                vti (k) = vf_min
            else
                tc (k) = tk (k) - tice
                if (hd_icefall) then
                    ! heymsfield and donner, 1990, jas
                    vti (k) = vi_fac * 3.29 * (qi (k) * den (k)) ** 0.16
                else
                    ! deng and mace, 2008, grl
                    vti (k) = (3. + log10 (qi (k) * den (k))) * (tc (k) * (aa * tc (k) + bb) + cc) + dd * tc (k) + ee
                    vti (k) = vi0 * exp (log_10 * vti (k))
                endif
                vti (k) = min (vi_max, max (vf_min, vti (k)))
            endif
        enddo
    endif

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    real, intent (in) :: qv, ql, qr, qi, qs, qg
    
=======

    if (const_vs) then
        vts (:) = vs_fac ! 1. ifs_2016
    else
        do k = ks, ke
            if (qs (k) < ths) then
                vts (k) = vf_min
            else
                vts (k) = vs_fac * vcons * rhof (k) * exp (0.0625 * log (qs (k) * den (k) / norms))
                vts (k) = min (vs_max, max (vf_min, vts (k)))
            endif
        enddo
    endif

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    real :: q_liq, q_sol
    
    q_liq = ql + qr
    q_sol = qi + qs + qg
    mhc6 = mhc (qv, q_liq, q_sol)
    
end function mhc6
=======

    if (const_vg) then
        vtg (:) = vg_fac ! 2.
    else
        if (do_hail) then
            do k = ks, ke
                if (qg (k) < thg) then
                    vtg (k) = vf_min
                else
                    vtg (k) = vg_fac * vconh * rhof (k) * sqrt (sqrt (sqrt (qg (k) * den (k) / normh)))
                    vtg (k) = min (vg_max, max (vf_min, vtg (k)))
                endif
            enddo
        else
            do k = ks, ke
                if (qg (k) < thg) then
                    vtg (k) = vf_min
                else
                    vtg (k) = vg_fac * vcong * rhof (k) * sqrt (sqrt (sqrt (qg (k) * den (k) / normg)))
                    vtg (k) = min (vg_max, max (vf_min, vtg (k)))
                endif
            enddo
        endif
    endif

end subroutine fall_speed
>>>>>>> user/lnz/shield2022

! =======================================================================
! moist total energy
! =======================================================================

<<<<<<< HEAD
function mte (qv, ql, qr, qi, qs, qg, tk, dp, moist_q)
    
    implicit none
    
    real (kind = r8) :: mte
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    logical, intent (in) :: moist_q
    
    real, intent (in) :: qv, ql, qr, qi, qs, qg, dp
    
    real (kind = r8), intent (in) :: tk
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real :: q_liq, q_sol, q_cond
    
    real (kind = r8) :: cvm, con_r8
    
    q_liq = ql + qr
    q_sol = qi + qs + qg
    q_cond = q_liq + q_sol
    con_r8 = one_r8 - (qv + q_cond)
    if (moist_q) then
        cvm = mhc (con_r8, qv, q_liq, q_sol)
=======
subroutine setupm

    implicit none

    real :: gcon, cd, scm3, pisq, act (8)
    real :: vdifu, tcond
    real :: visk
    real :: ch2o, hltf
    real :: hlts, hltc, ri50

    real, parameter :: gam263 = 1.456943, gam275 = 1.608355, gam290 = 1.827363, &
        gam325 = 2.54925, gam350 = 3.323363, gam380 = 4.694155, &
        gam425 = 8.285063, gam450 = 11.631769, gam480 = 17.837789, &
        gam625 = 184.860962, gam680 = 496.604067

    real, parameter :: acc (3) = (/ 5.0, 2.0, 0.5 /)

    real den_rc

    integer :: i, k

    pie = 4. * atan (1.0)

    ! s. klein's formular (eq 16) from am2

    fac_rc = (4. / 3.) * pie * rhor * rthresh ** 3

    if (prog_ccn) then
        ! if (master) write (*, *) 'prog_ccn option is .t.'
    else
        den_rc = fac_rc * ccn_o * 1.e6
        ! if (master) write (*, *) 'mp: for ccn_o = ', ccn_o, 'ql_rc = ', den_rc
        den_rc = fac_rc * ccn_l * 1.e6
        ! if (master) write (*, *) 'mp: for ccn_l = ', ccn_l, 'ql_rc = ', den_rc
    endif

    vdifu = 2.11e-5
    tcond = 2.36e-2

    visk = 1.259e-5
    hlts = hlv + hlf
    hltc = hlv
    hltf = hlf

    ch2o = c_liq
    ri50 = 1.e-4

    pisq = pie * pie
    scm3 = (visk / vdifu) ** (1. / 3.)

    cracs = pisq * rnzr * rnzs * rhos
    csacr = pisq * rnzr * rnzs * rhor
    if (do_hail) then
        cgacr = pisq * rnzr * rnzh * rhor
        cgacs = pisq * rnzh * rnzs * rhos
>>>>>>> user/lnz/shield2022
    else
        cvm = mhc (qv, q_liq, q_sol)
    endif
<<<<<<< HEAD
    mte = rgrav * cvm * c_air * tk * dp
    
end function mte

! =======================================================================
! moist total energy and total water
! =======================================================================

subroutine mtetw (ks, ke, qv, ql, qr, qi, qs, qg, tz, ua, va, wa, delp, &
        gsize, dte, water, rain, ice, snow, graupel, dts, te, tw, te_b, tw_b, &
        moist_q, hydrostatic, te_loss)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    logical, intent (in) :: moist_q, hydrostatic
    
    real, intent (in) :: gsize, water, rain, ice, snow, graupel, dts
    
    real, intent (in), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg, ua, va, wa, delp
    
    real (kind = r8), intent (in) :: dte
    
    real (kind = r8), intent (in), dimension (ks:ke) :: tz
    
    real (kind = r8), intent (out) :: te_b, tw_b
    
    real (kind = r8), intent (out), optional :: te_loss
    
    real (kind = r8), intent (out), dimension (ks:ke) :: te, tw
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    real :: q_cond
    
    real (kind = r8) :: con_r8
    
    real, dimension (ks:ke) :: q_liq, q_sol
    
    real (kind = r8), dimension (ks:ke) :: cvm

     do k = ks, ke
         q_liq (k) = ql (k) + qr (k)
         q_sol (k) = qi (k) + qs (k) + qg (k)
         q_cond = q_liq (k) + q_sol (k)
         con_r8 = one_r8 - (qv (k) + q_cond)
         if (moist_q) then
             cvm (k) = mhc (con_r8, qv (k), q_liq (k), q_sol (k))
         else
             cvm (k) = mhc (qv (k), q_liq (k), q_sol (k))
         endif
         te (k) = (cvm (k) * tz (k) + lv00 * qv (k) - li00 * q_sol (k)) * c_air
         if (hydrostatic) then
             te (k) = te (k) + 0.5 * (ua (k) ** 2 + va (k) ** 2)
         else
             te (k) = te (k) + 0.5 * (ua (k) ** 2 + va (k) ** 2 + wa (k) ** 2)
         endif
         te (k) = rgrav * te (k) * delp (k) * gsize ** 2.0
         tw (k) = rgrav * (qv (k) + q_cond) * delp (k) * gsize ** 2.0
     enddo
     te_b = (dte - li00 * c_air * (ice + snow + graupel) * dts / 86400) * gsize ** 2.0
     tw_b = (water + rain + ice + snow + graupel) * dts / 86400 * gsize ** 2.0

     if (present (te_loss)) then
          ! total energy change due to sedimentation and its heating
          te_loss = dte * gsize ** 2.0
     endif
=======
    cgacs = cgacs * c_pgacs

    ! act: 1 - 2:racs (s - r) ; 3 - 4:sacr (r - s) ;
    ! 5 - 6:gacr (r - g) ; 7 - 8:gacs (s - g)

    act (1) = pie * rnzs * rhos
    act (2) = pie * rnzr * rhor
    if (do_hail) then
        act (6) = pie * rnzh * rhoh
    else
        act (6) = pie * rnzg * rhog
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

    gcon = 40.74 * sqrt (sfcrho) ! 44.628

    csacw = pie * rnzs * clin * gam325 / (4. * act (1) ** 0.8125)
    ! decreasing csacw to reduce cloud water --- > snow

    craci = pie * rnzr * alin * gam380 / (4. * act (2) ** 0.95)
    csaci = csacw * c_psaci

    if (do_hail) then
        cgacw = pie * rnzh * gam350 * gcon / (4. * act (6) ** 0.875)
    else
        cgacw = pie * rnzg * gam350 * gcon / (4. * act (6) ** 0.875)
    endif
    ! cgaci = cgacw * 0.1

    ! sjl, may 28, 2012
    cgaci = cgacw * 0.05
    ! sjl, may 28, 2012

    cracw = craci ! cracw = 3.27206196043822
    cracw = c_cracw * cracw

    ! subl and revp: five constants for three separate processes

    cssub (1) = 2. * pie * vdifu * tcond * rvgas * rnzs
    if (do_hail) then
        cgsub (1) = 2. * pie * vdifu * tcond * rvgas * rnzh
    else
        cgsub (1) = 2. * pie * vdifu * tcond * rvgas * rnzg
    endif
    crevp (1) = 2. * pie * vdifu * tcond * rvgas * rnzr
    cssub (2) = 0.78 / sqrt (act (1))
    cgsub (2) = 0.78 / sqrt (act (6))
    crevp (2) = 0.78 / sqrt (act (2))
    cssub (3) = 0.31 * scm3 * gam263 * sqrt (clin / visk) / act (1) ** 0.65625
    cgsub (3) = 0.31 * scm3 * gam275 * sqrt (gcon / visk) / act (6) ** 0.6875
    crevp (3) = 0.31 * scm3 * gam290 * sqrt (alin / visk) / act (2) ** 0.725
    cssub (4) = tcond * rvgas
    cssub (5) = hlts ** 2 * vdifu
    cgsub (4) = cssub (4)
    crevp (4) = cssub (4)
    cgsub (5) = cssub (5)
    crevp (5) = hltc ** 2 * vdifu

    cgfr (1) = 20.e2 * pisq * rnzr * rhor / act (2) ** 1.75
    cgfr (2) = 0.66

    ! smlt: five constants (lin et al. 1983)

    csmlt (1) = 2. * pie * tcond * rnzs / hltf
    csmlt (2) = 2. * pie * vdifu * rnzs * hltc / hltf
    csmlt (3) = cssub (2)
    csmlt (4) = cssub (3)
    csmlt (5) = ch2o / hltf

    ! gmlt: five constants

    if (do_hail) then
        cgmlt (1) = 2. * pie * tcond * rnzh / hltf
        cgmlt (2) = 2. * pie * vdifu * rnzh * hltc / hltf
    else
        cgmlt (1) = 2. * pie * tcond * rnzg / hltf
        cgmlt (2) = 2. * pie * vdifu * rnzg * hltc / hltf
    endif
    cgmlt (3) = cgsub (2)
    cgmlt (4) = cgsub (3)
    cgmlt (5) = ch2o / hltf

    es0 = e00
    ces0 = eps * es0

end subroutine setupm
>>>>>>> user/lnz/shield2022

end subroutine mtetw
    
! =======================================================================
! calculate heat capacities and latent heat coefficients
! =======================================================================

<<<<<<< HEAD
subroutine cal_mhc_lhc (ks, ke, qv, ql, qr, qi, qs, qg, q_liq, q_sol, &
        cvm, te8, tz, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: ks, ke
    
    real, intent (in), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    
    real (kind = r8), intent (in), dimension (ks:ke) :: tz
    
    real, intent (out), dimension (ks:ke) :: q_liq, q_sol, lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (out), dimension (ks:ke) :: cvm, te8
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: k
    
    do k = ks, ke
        q_liq (k) = ql (k) + qr (k)
        q_sol (k) = qi (k) + qs (k) + qg (k)
        cvm (k) = mhc (qv (k), q_liq (k), q_sol (k))
        te8 (k) = cvm (k) * tz (k) + lv00 * qv (k) - li00 * q_sol (k)
        lcpk (k) = (lv00 + d1_vap * tz (k)) / cvm (k)
        icpk (k) = (li00 + d1_ice * tz (k)) / cvm (k)
        tcpk (k) = (li20 + (d1_vap + d1_ice) * tz (k)) / cvm (k)
        tcp3 (k) = lcpk (k) + icpk (k) * min (1., dim (tice, tz (k)) / (tice - t_wfr))
    enddo

end subroutine cal_mhc_lhc
    
=======
subroutine gfdl_mp_init (input_nml_file, logunit)

    implicit none

    character (len = *), intent (in) :: input_nml_file (:)
    integer, intent (in) :: logunit

    logical :: exists

    read (input_nml_file, nml = gfdl_mp_nml)

    ! write version number and namelist to log file
    write (logunit, *) " ================================================================== "
    write (logunit, *) "gfdl_mp_mod"
    write (logunit, nml = gfdl_mp_nml)

    if (do_setup) then
        call setup_con
        call setupm
        do_setup = .false.
    endif

    g2 = 0.5 * grav
    log_10 = log (10.)

    if (do_warm_rain_mp) then
        t_wfr = t_min
    else
        t_wfr = t_ice - 40.0
    endif

    module_is_initialized = .true.

end subroutine gfdl_mp_init

! =======================================================================
! end of gfdl cloud microphysics
! =======================================================================

subroutine gfdl_mp_end

    implicit none

    deallocate (table)
    deallocate (table2)
    deallocate (table3)
    deallocate (tablew)
    deallocate (des)
    deallocate (des2)
    deallocate (des3)
    deallocate (desw)

    tables_are_initialized = .false.

end subroutine gfdl_mp_end

>>>>>>> user/lnz/shield2022
! =======================================================================
! update hydrometeors
! =======================================================================

<<<<<<< HEAD
subroutine update_qq (qv, ql, qr, qi, qs, qg, dqv, dql, dqr, dqi, dqs, dqg)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: dqv, dql, dqr, dqi, dqs, dqg
    
    real, intent (inout) :: qv, ql, qr, qi, qs, qg
    
    qv = qv + dqv
    ql = ql + dql
    qr = qr + dqr
    qi = qi + dqi
    qs = qs + dqs
    qg = qg + dqg
    
end subroutine update_qq
=======
subroutine setup_con

    implicit none

    ! master = (mpp_pe () .eq.mpp_root_pe ())

    if (.not. qsmith_tables_initialized) call qsmith_init

    qsmith_tables_initialized = .true.

end subroutine setup_con
>>>>>>> user/lnz/shield2022

! =======================================================================
! update hydrometeors and temperature
! =======================================================================

<<<<<<< HEAD
subroutine update_qt (qv, ql, qr, qi, qs, qg, dqv, dql, dqr, dqi, dqs, dqg, te8, &
        cvm, tk, lcpk, icpk, tcpk, tcp3)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: dqv, dql, dqr, dqi, dqs, dqg
    
    real (kind = r8), intent (in) :: te8
    
    real, intent (inout) :: qv, ql, qr, qi, qs, qg
    
    real, intent (out) :: lcpk, icpk, tcpk, tcp3
    
    real (kind = r8), intent (out) :: cvm, tk
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    qv = qv + dqv
    ql = ql + dql
    qr = qr + dqr
    qi = qi + dqi
    qs = qs + dqs
    qg = qg + dqg
    
    cvm = mhc (qv, ql, qr, qi, qs, qg)
    tk = (te8 - lv00 * qv + li00 * (qi + qs + qg)) / cvm
    
    lcpk = (lv00 + d1_vap * tk) / cvm
    icpk = (li00 + d1_ice * tk) / cvm
    tcpk = (li20 + (d1_vap + d1_ice) * tk) / cvm
    tcp3 = lcpk + icpk * min (1., dim (tice, tk) / (tice - t_wfr))
    
end subroutine update_qt
=======
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
>>>>>>> user/lnz/shield2022

! =======================================================================
! prepare saturation water vapor pressure tables
! =======================================================================

<<<<<<< HEAD
subroutine qs_init
    
    implicit none
    
=======
subroutine qsmith_init

    implicit none

    integer, parameter :: length = 2621

>>>>>>> user/lnz/shield2022
    integer :: i

    if (.not. tables_are_initialized) then
<<<<<<< HEAD
        
        allocate (table0 (length))
        allocate (table1 (length))
=======

        allocate (table (length))
>>>>>>> user/lnz/shield2022
        allocate (table2 (length))
        allocate (table3 (length))
        allocate (table4 (length))
        
        allocate (des0 (length))
        allocate (des1 (length))
        allocate (des2 (length))
        allocate (des3 (length))
<<<<<<< HEAD
        allocate (des4 (length))
        
        call qs_table0 (length)
        call qs_table1 (length)
        call qs_table2 (length)
        call qs_table3 (length)
        call qs_table4 (length)
        
=======
        allocate (desw (length))

        call qs_table (length)
        call qs_table2 (length)
        call qs_table3 (length)
        call qs_tablew (length)

>>>>>>> user/lnz/shield2022
        do i = 1, length - 1
            des0 (i) = max (0., table0 (i + 1) - table0 (i))
            des1 (i) = max (0., table1 (i + 1) - table1 (i))
            des2 (i) = max (0., table2 (i + 1) - table2 (i))
            des3 (i) = max (0., table3 (i + 1) - table3 (i))
            des4 (i) = max (0., table4 (i + 1) - table4 (i))
        enddo
        des0 (length) = des0 (length - 1)
        des1 (length) = des1 (length - 1)
        des2 (length) = des2 (length - 1)
        des3 (length) = des3 (length - 1)
<<<<<<< HEAD
        des4 (length) = des4 (length - 1)
        
        tables_are_initialized = .true.

    endif
    
end subroutine qs_init
=======
        desw (length) = desw (length - 1)

        tables_are_initialized = .true.

    endif

end subroutine qsmith_init

! =======================================================================
! compute the saturated specific humidity for table ii
! =======================================================================

real function wqs1 (ta, den)

    implicit none

    ! pure water phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqs1 = es / (rvgas * ta * den)

end function wqs1

! =======================================================================
! compute the gradient of saturated specific humidity for table ii
! =======================================================================

real function wqs2 (ta, den, dqdt)

    implicit none

    ! pure water phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    real, intent (in) :: ta, den
    real, intent (out) :: dqdt
    real :: es, ap1, tmin
    integer :: it

    tmin = table_ice - 160.

    if (.not. tables_are_initialized) call qsmith_init

    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqs2 = es / (rvgas * ta * den)
    it = ap1 - 0.5
    ! finite diff, del_t = 0.1:
    dqdt = 10. * (desw (it) + (ap1 - it) * (desw (it + 1) - desw (it))) / (rvgas * ta * den)

end function wqs2
>>>>>>> user/lnz/shield2022

! =======================================================================
! saturation water vapor pressure table, core function
! =======================================================================

<<<<<<< HEAD
subroutine qs_table_core (n, n_blend, do_smith_table, table)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n, n_blend
    
    logical, intent (in) :: do_smith_table
    
    real, intent (out), dimension (n) :: table
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i
    integer, parameter :: n_min = 1600
    
    real (kind = r8) :: delt = 0.1
    real (kind = r8) :: tmin, tem, esh
    real (kind = r8) :: wice, wh2o, fac0, fac1, fac2
    real (kind = r8) :: esbasw, tbasw, esbasi, a, b, c, d, e
    real (kind = r8) :: esupc (n_blend)
    
    esbasw = 1013246.0
    tbasw = tice + 100.
    esbasi = 6107.1
    tmin = tice - n_min * delt
    
    ! -----------------------------------------------------------------------
    ! compute es over ice between - (n_min * delt) deg C and 0 deg C
    ! -----------------------------------------------------------------------
    
    if (do_smith_table) then
        do i = 1, n_min
            tem = tmin + delt * real (i - 1)
            a = - 9.09718 * (tice / tem - 1.)
            b = - 3.56654 * log10 (tice / tem)
            c = 0.876793 * (1. - tem / tice)
            e = log10 (esbasi)
            table (i) = 0.1 * exp ((a + b + c + e) * log (10.))
        enddo
    else
        do i = 1, n_min
            tem = tmin + delt * real (i - 1)
            fac0 = (tem - tice) / (tem * tice)
            fac1 = fac0 * li2
            fac2 = (d2_ice * log (tem / tice) + fac1) / rvgas
            table (i) = e00 * exp (fac2)
        enddo
    endif
    
    ! -----------------------------------------------------------------------
    ! compute es over water between - (n_blend * delt) deg C and [ (n - n_min - 1) * delt] deg C
    ! -----------------------------------------------------------------------
    
    if (do_smith_table) then
        do i = 1, n - n_min + n_blend
            tem = tice + delt * (real (i - 1) - n_blend)
            a = - 7.90298 * (tbasw / tem - 1.)
            b = 5.02808 * log10 (tbasw / tem)
            c = - 1.3816e-7 * (exp ((1. - tem / tbasw) * 11.344 * log (10.)) - 1.)
            d = 8.1328e-3 * (exp ((tbasw / tem - 1.) * (- 3.49149) * log (10.)) - 1.)
            e = log10 (esbasw)
            esh = 0.1 * exp ((a + b + c + d + e) * log (10.))
            if (i .le. n_blend) then
                esupc (i) = esh
            else
                table (i + n_min - n_blend) = esh
            endif
        enddo
    else
        do i = 1, n - n_min + n_blend
            tem = tice + delt * (real (i - 1) - n_blend)
            fac0 = (tem - tice) / (tem * tice)
            fac1 = fac0 * lv0
            fac2 = (dc_vap * log (tem / tice) + fac1) / rvgas
            esh = e00 * exp (fac2)
            if (i .le. n_blend) then
                esupc (i) = esh
            else
                table (i + n_min - n_blend) = esh
            endif
        enddo
    endif
    
    ! -----------------------------------------------------------------------
    ! derive blended es over ice and supercooled water between - (n_blend * delt) deg C and 0 deg C
    ! -----------------------------------------------------------------------
    
    do i = 1, n_blend
        tem = tice + delt * (real (i - 1) - n_blend)
        wice = 1.0 / (delt * n_blend) * (tice - tem)
        wh2o = 1.0 / (delt * n_blend) * (tem - tice + delt * n_blend)
        table (i + n_min - n_blend) = wice * table (i + n_min - n_blend) + wh2o * esupc (i)
    enddo
    
end subroutine qs_table_core
=======
subroutine wqs2_vect (is, ie, ta, den, wqsat, dqdt)

    implicit none

    ! pure water phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    integer, intent (in) :: is, ie

    real, intent (in), dimension (is:ie) :: ta, den

    real, intent (out), dimension (is:ie) :: wqsat, dqdt

    real :: es, ap1, tmin

    integer :: i, it

    tmin = t_ice - 160.

    do i = is, ie
        ap1 = 10. * dim (ta (i), tmin) + 1.
        ap1 = min (2621., ap1)
        it = ap1
        es = tablew (it) + (ap1 - it) * desw (it)
        wqsat (i) = es / (rvgas * ta (i) * den (i))
        it = ap1 - 0.5
        ! finite diff, del_t = 0.1:
        dqdt (i) = 10. * (desw (it) + (ap1 - it) * (desw (it + 1) - desw (it))) / (rvgas * ta (i) * den (i))
    enddo

end subroutine wqs2_vect

! =======================================================================
! compute wet buld temperature
! =======================================================================

real function wet_bulb (q, t, den)

    implicit none

    real, intent (in) :: t, q, den

    real :: qs, tp, dqdt

    wet_bulb = t
    qs = wqs2 (wet_bulb, den, dqdt)
    tp = 0.5 * (qs - q) / (1. + lcp * dqdt) * lcp
    wet_bulb = wet_bulb - tp

    ! tp is negative if super - saturated
    if (tp > 0.01) then
        qs = wqs2 (wet_bulb, den, dqdt)
        tp = (qs - q) / (1. + lcp * dqdt) * lcp
        wet_bulb = wet_bulb - tp
    endif

end function wet_bulb

! =======================================================================
! compute the saturated specific humidity for table iii
! =======================================================================

real function iqs1 (ta, den)

    implicit none

    ! water - ice phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    iqs1 = es / (rvgas * ta * den)

end function iqs1
>>>>>>> user/lnz/shield2022

! =======================================================================
! saturation water vapor pressure table 0, water only
! useful for idealized experiments
! it can also be used in warm rain microphyscis only
! =======================================================================

<<<<<<< HEAD
subroutine qs_table0 (n)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i
    
    real (kind = r8) :: delt = 0.1
    real (kind = r8) :: tmin, tem, fac0, fac1, fac2
    
    tmin = tice - 160.
    
    ! -----------------------------------------------------------------------
    ! compute es over water only
    ! -----------------------------------------------------------------------
    
    do i = 1, n
        tem = tmin + delt * real (i - 1)
        fac0 = (tem - tice) / (tem * tice)
        fac1 = fac0 * lv0
        fac2 = (dc_vap * log (tem / tice) + fac1) / rvgas
        table0 (i) = e00 * exp (fac2)
    enddo
    
end subroutine qs_table0
=======
real function iqs2 (ta, den, dqdt)

    implicit none

    ! water - ice phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    real (kind = r_grid), intent (in) :: ta
    real, intent (in) :: den
    real, intent (out) :: dqdt
    real (kind = r_grid) :: tmin, es, ap1
    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    iqs2 = es / (rvgas * ta * den)
    it = ap1 - 0.5
    dqdt = 10. * (des2 (it) + (ap1 - it) * (des2 (it + 1) - des2 (it))) / (rvgas * ta * den)

end function iqs2

! =======================================================================
! compute the gradient of saturated specific humidity for table iii
! =======================================================================

real function qs1d_moist (ta, qv, pa, dqdt)

    implicit none

    real, intent (in) :: ta, pa, qv

    real, intent (out) :: dqdt

    real :: es, ap1, tmin, eps10

    integer :: it

    tmin = table_ice - 160.
    eps10 = 10. * eps
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    qs1d_moist = eps * es * (1. + zvir * qv) / pa
    it = ap1 - 0.5
    dqdt = eps10 * (des2 (it) + (ap1 - it) * (des2 (it + 1) - des2 (it))) * (1. + zvir * qv) / pa

end function qs1d_moist
>>>>>>> user/lnz/shield2022

! =======================================================================
! saturation water vapor pressure table 1, water and ice
! blended between -20 deg C and 0 deg C
! the most realistic saturation water vapor pressure for the full temperature range
! =======================================================================

<<<<<<< HEAD
subroutine qs_table1 (n)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n
    
    call qs_table_core (n, 200, .false., table1)
    
end subroutine qs_table1
=======
real function wqsat2_moist (ta, qv, pa, dqdt)

    implicit none

    real, intent (in) :: ta, pa, qv

    real, intent (out) :: dqdt

    real :: es, ap1, tmin, eps10

    integer :: it

    tmin = table_ice - 160.
    eps10 = 10. * eps
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqsat2_moist = eps * es * (1. + zvir * qv) / pa
    it = ap1 - 0.5
    dqdt = eps10 * (desw (it) + (ap1 - it) * (desw (it + 1) - desw (it))) * (1. + zvir * qv) / pa

end function wqsat2_moist
>>>>>>> user/lnz/shield2022

! =======================================================================
! saturation water vapor pressure table 2, water and ice
! same as table 1, but the blending is replaced with smoothing around 0 deg C
! it is not designed for mixed-phase cloud microphysics
! used for ice microphysics (< 0 deg C) or warm rain microphysics (> 0 deg C)
! =======================================================================

<<<<<<< HEAD
subroutine qs_table2 (n)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n
    
    call qs_table_core (n, 0, .false., table2)
    
end subroutine qs_table2
=======
real function wqsat_moist (ta, qv, pa)

    implicit none

    real, intent (in) :: ta, pa, qv

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqsat_moist = eps * es * (1. + zvir * qv) / pa

end function wqsat_moist
>>>>>>> user/lnz/shield2022

! =======================================================================
! saturation water vapor pressure table 3, water and ice
! blended between -20 deg C and 0 deg C
! the same as table 1, but from smithsonian meteorological tables page 350
! =======================================================================

<<<<<<< HEAD
subroutine qs_table3 (n)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n
    
    call qs_table_core (n, 200, .true., table3)
    
end subroutine qs_table3
=======
real function qs1d_m (ta, qv, pa)

    implicit none

    real, intent (in) :: ta, pa, qv

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    qs1d_m = eps * es * (1. + zvir * qv) / pa

end function qs1d_m
>>>>>>> user/lnz/shield2022

! =======================================================================
! saturation water vapor pressure table 4, water and ice
! same as table 3, but the blending is replaced with smoothing around 0 deg C
! the same as table 2, but from smithsonian meteorological tables page 350
! =======================================================================

<<<<<<< HEAD
subroutine qs_table4 (n)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: n
    
    call qs_table_core (n, 0, .true., table4)
    
end subroutine qs_table4
=======
real function d_sat (ta, den)

    implicit none

    real, intent (in) :: ta, den

    real :: es_w, es_i, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es_w = tablew (it) + (ap1 - it) * desw (it)
    es_i = table2 (it) + (ap1 - it) * des2 (it)
    d_sat = dim (es_w, es_i) / (rvgas * ta * den) ! take positive difference

end function d_sat
>>>>>>> user/lnz/shield2022

! =======================================================================
! compute the saturated specific humidity, core function
! =======================================================================

<<<<<<< HEAD
function qs_core (length, tk, den, dqdt, table, des)
    
    implicit none
    
    real :: qs_core
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: length
    
    real, intent (in) :: tk, den
    
    real, intent (in), dimension (length) :: table, des
    
    real, intent (out) :: dqdt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: it
    
    real :: es, ap1, tmin
    
    if (.not. tables_are_initialized) call qs_init
    
    tmin = tice - 160.
    ap1 = 10. * dim (tk, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table (it) + (ap1 - it) * des (it)
    qs_core = es / (rvgas * tk * den)
    it = ap1 - 0.5
    dqdt = 10. * (des (it) + (ap1 - it) * (des (it + 1) - des (it))) / (rvgas * tk * den)
    
end function qs_core
=======
real function esw_table (ta)

    implicit none

    real, intent (in) :: ta

    real :: ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    esw_table = tablew (it) + (ap1 - it) * desw (it)

end function esw_table

! =======================================================================
! compute the saturated water vapor pressure for table iii
! =======================================================================

real function es2_table (ta)

    implicit none

    real, intent (in) :: ta

    real :: ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es2_table = table2 (it) + (ap1 - it) * des2 (it)

end function es2_table
>>>>>>> user/lnz/shield2022

! =======================================================================
! compute the saturated specific humidity based on table 0, water only
! useful for idealized experiments
! it can also be used in warm rain microphyscis only
! =======================================================================

<<<<<<< HEAD
function wqs_trho (tk, den, dqdt)
    
    implicit none
    
    real :: wqs_trho
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk, den
    
    real, intent (out) :: dqdt
    
    wqs_trho = qs_core (length, tk, den, dqdt, table0, des0)
    
end function wqs_trho
=======
subroutine esw_table1d (ta, es, n)

    implicit none

    integer, intent (in) :: n

    real, intent (in) :: ta (n)

    real, intent (out) :: es (n)

    real :: ap1, tmin

    integer :: i, it

    tmin = table_ice - 160.

    do i = 1, n
        ap1 = 10. * dim (ta (i), tmin) + 1.
        ap1 = min (2621., ap1)
        it = ap1
        es (i) = tablew (it) + (ap1 - it) * desw (it)
    enddo

end subroutine esw_table1d
>>>>>>> user/lnz/shield2022

! =======================================================================
! compute the saturated specific humidity based on table 1, water and ice
! the most realistic saturation water vapor pressure for the full temperature range
! =======================================================================

<<<<<<< HEAD
function mqs_trho (tk, den, dqdt)
    
    implicit none
    
    real :: mqs_trho
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk, den
    
    real, intent (out) :: dqdt
    
    mqs_trho = qs_core (length, tk, den, dqdt, table1, des1)
    
end function mqs_trho
=======
subroutine es2_table1d (ta, es, n)

    implicit none

    integer, intent (in) :: n

    real, intent (in) :: ta (n)

    real, intent (out) :: es (n)

    real :: ap1, tmin

    integer :: i, it

    tmin = table_ice - 160.

    do i = 1, n
        ap1 = 10. * dim (ta (i), tmin) + 1.
        ap1 = min (2621., ap1)
        it = ap1
        es (i) = table2 (it) + (ap1 - it) * des2 (it)
    enddo

end subroutine es2_table1d
>>>>>>> user/lnz/shield2022

! =======================================================================
! compute the saturated specific humidity based on table 2, water and ice
! it is not designed for mixed-phase cloud microphysics
! used for ice microphysics (< 0 deg C) or warm rain microphysics (> 0 deg C)
! =======================================================================

<<<<<<< HEAD
function iqs_trho (tk, den, dqdt)
    
    implicit none
    
    real :: iqs_trho
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk, den
    
    real, intent (out) :: dqdt
    
    iqs_trho = qs_core (length, tk, den, dqdt, table2, des2)
    
end function iqs_trho
=======
subroutine es3_table1d (ta, es, n)

    implicit none

    integer, intent (in) :: n

    real, intent (in) :: ta (n)

    real, intent (out) :: es (n)

    real :: ap1, tmin

    integer :: i, it

    tmin = table_ice - 160.

    do i = 1, n
        ap1 = 10. * dim (ta (i), tmin) + 1.
        ap1 = min (2621., ap1)
        it = ap1
        es (i) = table3 (it) + (ap1 - it) * des3 (it)
    enddo

end subroutine es3_table1d
>>>>>>> user/lnz/shield2022

! =======================================================================
! compute the saturated specific humidity based on table 0, water only
! useful for idealized experiments
! it can also be used in warm rain microphyscis only
! =======================================================================

<<<<<<< HEAD
function wqs_ptqv (tk, pa, qv, dqdt)
    
    implicit none
    
    real :: wqs_ptqv
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk, pa, qv
    
    real, intent (out) :: dqdt
    
=======
subroutine qs_tablew (n)

    implicit none

    integer, intent (in) :: n

    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: tmin, tem, fac0, fac1, fac2

    integer :: i

    tmin = table_ice - 160.

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    real :: den
    
    den = pa / (rdgas * tk * (1. + zvir * qv))
    
    wqs_ptqv = wqs (tk, den, dqdt)
    
end function wqs_ptqv
=======

    do i = 1, n
        tem = tmin + delt * real (i - 1)
        fac0 = (tem - t_ice) / (tem * t_ice)
        fac1 = fac0 * lv0
        fac2 = (dc_vap * log (tem / t_ice) + fac1) / rvgas
        tablew (i) = e00 * exp (fac2)
    enddo

end subroutine qs_tablew
>>>>>>> user/lnz/shield2022

! =======================================================================
! compute the saturated specific humidity based on table 1, water and ice
! the most realistic saturation water vapor pressure for the full temperature range
! =======================================================================

<<<<<<< HEAD
function mqs_ptqv (tk, pa, qv, dqdt)
    
    implicit none
    
    real :: mqs_ptqv
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk, pa, qv
    
    real, intent (out) :: dqdt
    
=======
subroutine qs_table2 (n)

    implicit none

    integer, intent (in) :: n

    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: tmin, tem0, tem1, fac0, fac1, fac2

    integer :: i, i0, i1

    tmin = table_ice - 160.

    do i = 1, n
        tem0 = tmin + delt * real (i - 1)
        fac0 = (tem0 - t_ice) / (tem0 * t_ice)
        if (i <= 1600) then
            ! -----------------------------------------------------------------------
            ! compute es over ice between - 160 deg c and 0 deg c.
            ! -----------------------------------------------------------------------
            fac1 = fac0 * li2
            fac2 = (d2ice * log (tem0 / t_ice) + fac1) / rvgas
        else
            ! -----------------------------------------------------------------------
            ! compute es over water between 0 deg c and 102 deg c.
            ! -----------------------------------------------------------------------
            fac1 = fac0 * lv0
            fac2 = (dc_vap * log (tem0 / t_ice) + fac1) / rvgas
        endif
        table2 (i) = e00 * exp (fac2)
    enddo

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    real :: den
    
    den = pa / (rdgas * tk * (1. + zvir * qv))
    
    mqs_ptqv = mqs (tk, den, dqdt)
    
end function mqs_ptqv
=======

    i0 = 1600
    i1 = 1601
    tem0 = 0.25 * (table2 (i0 - 1) + 2. * table (i0) + table2 (i0 + 1))
    tem1 = 0.25 * (table2 (i1 - 1) + 2. * table (i1) + table2 (i1 + 1))
    table2 (i0) = tem0
    table2 (i1) = tem1

end subroutine qs_table2
>>>>>>> user/lnz/shield2022

! =======================================================================
! compute the saturated specific humidity based on table 2, water and ice
! it is not designed for mixed-phase cloud microphysics
! used for ice microphysics (< 0 deg C) or warm rain microphysics (> 0 deg C)
! =======================================================================

<<<<<<< HEAD
function iqs_ptqv (tk, pa, qv, dqdt)
    
    implicit none
    
    real :: iqs_ptqv
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: tk, pa, qv
    
    real, intent (out) :: dqdt
    
=======
subroutine qs_table3 (n)

    implicit none

    integer, intent (in) :: n

    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: esbasw, tbasw, esbasi, tmin, tem, aa, b, c, d, e
    real (kind = r_grid) :: tem0, tem1

    integer :: i, i0, i1

    esbasw = 1013246.0
    tbasw = table_ice + 100.
    esbasi = 6107.1
    tmin = table_ice - 160.

    do i = 1, n
        tem = tmin + delt * real (i - 1)
        ! if (i <= 1600) then
        if (i <= 1580) then ! change to - 2 c
            ! -----------------------------------------------------------------------
            ! compute es over ice between - 160 deg c and 0 deg c.
            ! see smithsonian meteorological tables page 350.
            ! -----------------------------------------------------------------------
            aa = - 9.09718 * (table_ice / tem - 1.)
            b = - 3.56654 * log10 (table_ice / tem)
            c = 0.876793 * (1. - tem / table_ice)
            e = log10 (esbasi)
            table3 (i) = 0.1 * 10 ** (aa + b + c + e)
        else
            ! -----------------------------------------------------------------------
            ! compute es over water between - 2 deg c and 102 deg c.
            ! see smithsonian meteorological tables page 350.
            ! -----------------------------------------------------------------------
            aa = - 7.90298 * (tbasw / tem - 1.)
            b = 5.02808 * log10 (tbasw / tem)
            c = - 1.3816e-7 * (10 ** ((1. - tem / tbasw) * 11.344) - 1.)
            d = 8.1328e-3 * (10 ** ((tbasw / tem - 1.) * (- 3.49149)) - 1.)
            e = log10 (esbasw)
            table3 (i) = 0.1 * 10 ** (aa + b + c + d + e)
        endif
    enddo

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    real :: den
    
    den = pa / (rdgas * tk * (1. + zvir * qv))
    
    iqs_ptqv = iqs (tk, den, dqdt)
    
end function iqs_ptqv
=======

    i0 = 1580
    i1 = 1581
    tem0 = 0.25 * (table3 (i0 - 1) + 2. * table (i0) + table3 (i0 + 1))
    tem1 = 0.25 * (table3 (i1 - 1) + 2. * table (i1) + table3 (i1 + 1))
    table3 (i0) = tem0
    table3 (i1) = tem1

end subroutine qs_table3
>>>>>>> user/lnz/shield2022

! =======================================================================
! compute the saturated specific humidity based on table 1, water and ice
! the most realistic saturation water vapor pressure for the full temperature range
! it is the 3d version of "mqs"
! =======================================================================

<<<<<<< HEAD
subroutine mqs3d (im, km, ks, tk, pa, qv, qs, dqdt)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    integer, intent (in) :: im, km, ks
    
    real, intent (in), dimension (im, ks:km) :: tk, pa, qv
    
    real, intent (out), dimension (im, ks:km) :: qs
    
    real, intent (out), dimension (im, ks:km), optional :: dqdt
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    integer :: i, k
    
    real :: dqdt0
    
    if (present (dqdt)) then
        do k = ks, km
            do i = 1, im
                qs (i, k) = mqs (tk (i, k), pa (i, k), qv (i, k), dqdt (i, k))
            enddo
        enddo
    else
        do k = ks, km
            do i = 1, im
                qs (i, k) = mqs (tk (i, k), pa (i, k), qv (i, k), dqdt0)
            enddo
        enddo
    endif
    
end subroutine mqs3d
=======
real function qs_blend (t, p, q)

    implicit none

    real, intent (in) :: t, p, q

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (t, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table (it) + (ap1 - it) * des (it)
    qs_blend = eps * es * (1. + zvir * q) / p

end function qs_blend
>>>>>>> user/lnz/shield2022

! =======================================================================
! compute wet buld temperature, core function
! Knox et al. (2017)
! =======================================================================

<<<<<<< HEAD
function wet_bulb_core (qv, tk, den, lcp)
    
    implicit none
    
    real :: wet_bulb_core
    
=======
subroutine qs_table (n)

    implicit none

    integer, intent (in) :: n

    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: tmin, tem, esh20
    real (kind = r_grid) :: wice, wh2o, fac0, fac1, fac2
    real (kind = r_grid) :: esupc (200)

    integer :: i

    tmin = table_ice - 160.

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    real, intent (in) :: qv, tk, den, lcp
    
=======

    do i = 1, 1600
        tem = tmin + delt * real (i - 1)
        fac0 = (tem - t_ice) / (tem * t_ice)
        fac1 = fac0 * li2
        fac2 = (d2ice * log (tem / t_ice) + fac1) / rvgas
        table (i) = e00 * exp (fac2)
    enddo

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    logical :: do_adjust = .false.
    
    real :: factor = 1 / 3.
    real :: qsat, tp, dqdt
    
    wet_bulb_core = tk
    qsat = wqs (wet_bulb_core, den, dqdt)
    tp = factor * (qsat - qv) / (1. + lcp * dqdt) * lcp
    wet_bulb_core = wet_bulb_core - tp
    
    if (do_adjust .and. tp .gt. 0.0) then
        qsat = wqs (wet_bulb_core, den, dqdt)
        tp = (qsat - qv) / (1. + lcp * dqdt) * lcp
        wet_bulb_core = wet_bulb_core - tp
    endif
    
end function wet_bulb_core
=======

    do i = 1, 1221
        tem = 253.16 + delt * real (i - 1)
        fac0 = (tem - t_ice) / (tem * t_ice)
        fac1 = fac0 * lv0
        fac2 = (dc_vap * log (tem / t_ice) + fac1) / rvgas
        esh20 = e00 * exp (fac2)
        if (i <= 200) then
            esupc (i) = esh20
        else
            table (i + 1400) = esh20
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! derive blended es over ice and supercooled water between - 20 deg c and 0 deg c
    ! -----------------------------------------------------------------------

    do i = 1, 200
        tem = 253.16 + delt * real (i - 1)
        wice = 0.05 * (table_ice - tem)
        wh2o = 0.05 * (tem - 253.16)
        table (i + 1400) = wice * table (i + 1400) + wh2o * esupc (i)
    enddo

end subroutine qs_table
>>>>>>> user/lnz/shield2022

! =======================================================================
! compute wet buld temperature, dry air case
! =======================================================================

<<<<<<< HEAD
function wet_bulb_dry (qv, tk, den)
    
    implicit none
    
    real :: wet_bulb_dry
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
    
    real, intent (in) :: qv, tk, den
    
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
    
    real :: lcp
    
    lcp = hlv / cp_air
    
    wet_bulb_dry = wet_bulb_core (qv, tk, den, lcp)
    
end function wet_bulb_dry
=======
subroutine qsmith (im, km, ks, t, p, q, qs, dqdt)

    implicit none

    integer, intent (in) :: im, km, ks

    real, intent (in), dimension (im, km) :: t, p, q

    real, intent (out), dimension (im, km) :: qs

    real, intent (out), dimension (im, km), optional :: dqdt

    real :: eps10, ap1, tmin

    real, dimension (im, km) :: es

    integer :: i, k, it

    tmin = table_ice - 160.
    eps10 = 10. * eps

    if (.not. tables_are_initialized) then
        call qsmith_init
    endif

    do k = ks, km
        do i = 1, im
            ap1 = 10. * dim (t (i, k), tmin) + 1.
            ap1 = min (2621., ap1)
            it = ap1
            es (i, k) = table (it) + (ap1 - it) * des (it)
            qs (i, k) = eps * es (i, k) * (1. + zvir * q (i, k)) / p (i, k)
        enddo
    enddo

    if (present (dqdt)) then
        do k = ks, km
            do i = 1, im
                ap1 = 10. * dim (t (i, k), tmin) + 1.
                ap1 = min (2621., ap1) - 0.5
                it = ap1
                dqdt (i, k) = eps10 * (des (it) + (ap1 - it) * (des (it + 1) - des (it))) * (1. + zvir * q (i, k)) / p (i, k)
            enddo
        enddo
    endif

end subroutine qsmith
>>>>>>> user/lnz/shield2022

! =======================================================================
! compute wet buld temperature, moist air case
! =======================================================================

<<<<<<< HEAD
function wet_bulb_moist (qv, ql, qi, qr, qs, qg, tk, den)
    
    implicit none
    
    real :: wet_bulb_moist
    
=======
subroutine neg_adj (ks, ke, pt, dp, qv, ql, qr, qi, qs, qg, cond)

    implicit none

    integer, intent (in) :: ks, ke
    real, intent (in), dimension (ks:ke) :: dp
    real (kind = r_grid), intent (inout), dimension (ks:ke) :: pt
    real, intent (inout), dimension (ks:ke) :: qv, ql, qr, qi, qs, qg
    real, intent (out) :: cond

    real, dimension (ks:ke) :: lcpk, icpk

    real :: dq, cvm

    integer :: k

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    real, intent (in) :: qv, ql, qi, qr, qs, qg, tk, den
    
=======

    do k = ks, ke
        cvm = 1. + qv (k) * c1_vap + (qr (k) + ql (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
        lcpk (k) = (lv00 + d1_vap * pt (k)) / cvm
        icpk (k) = (li00 + d1_ice * pt (k)) / cvm
    enddo

    cond = 0

    do k = ks, ke

        ! -----------------------------------------------------------------------
        ! ice phase:
        ! -----------------------------------------------------------------------

        ! if cloud ice < 0, borrow from snow
        if (qi (k) < 0.) then
            qs (k) = qs (k) + qi (k)
            qi (k) = 0.
        endif
        ! if snow < 0, borrow from graupel
        if (qs (k) < 0.) then
            qg (k) = qg (k) + qs (k)
            qs (k) = 0.
        endif
        ! if graupel < 0, borrow from rain
        if (qg (k) < 0.) then
            qr (k) = qr (k) + qg (k)
            pt (k) = pt (k) - qg (k) * icpk (k) ! heating
            qg (k) = 0.
        endif

        ! -----------------------------------------------------------------------
        ! liquid phase:
        ! -----------------------------------------------------------------------

        ! if rain < 0, borrow from cloud water
        if (qr (k) < 0.) then
            ql (k) = ql (k) + qr (k)
            qr (k) = 0.
        endif
        ! if cloud water < 0, borrow from water vapor
        if (ql (k) < 0.) then
            cond = cond - ql (k) * dp (k)
            qv (k) = qv (k) + ql (k)
            pt (k) = pt (k) - ql (k) * lcpk (k) ! heating
            ql (k) = 0.
        endif

    enddo

>>>>>>> user/lnz/shield2022
    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------
<<<<<<< HEAD
    
    real :: lcp, q_liq, q_sol
    
    real (kind = r8) :: cvm
    
    q_liq = ql + qr
    q_sol = qi + qs + qg
    cvm = mhc (qv, q_liq, q_sol)
    lcp = (lv00 + d1_vap * tk) / cvm
    
    wet_bulb_moist = wet_bulb_core (qv, tk, den, lcp)
    
end function wet_bulb_moist
=======

    do k = ks, ke - 1
        if (qv (k) < 0.) then
            qv (k + 1) = qv (k + 1) + qv (k) * dp (k) / dp (k + 1)
            qv (k) = 0.
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! bottom layer; borrow from above
    ! -----------------------------------------------------------------------

    if (qv (ke) < 0. .and. qv (ke - 1) > 0.) then
        dq = min (- qv (ke) * dp (ke), qv (ke - 1) * dp (ke - 1))
        qv (ke - 1) = qv (ke - 1) - dq / dp (ke - 1)
        qv (ke) = qv (ke) + dq / dp (ke)
    endif

end subroutine neg_adj
>>>>>>> user/lnz/shield2022

end module gfdl_mp_mod
