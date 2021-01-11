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

module external_sst_mod

  use fms_io_mod,          only: restart_file_type, register_restart_field
  use fms_io_mod,          only: restore_state
  use fms_io_mod,          only: read_data, get_var_att_value, field_size
  use fms_mod,             only: file_exist, mpp_error, FATAL
  use mpp_mod,             only: mpp_pe, mpp_root_pe
  use fv_arrays_mod,       only: fv_atmos_type
  use time_manager_mod,    only: time_type, get_date, operator(==), &
                                 operator(>), operator(<)
  use get_cal_time_mod,    only: get_cal_time

#ifdef NO_GFDL_SHARED
!----------------- Public Data -----------------------------------
integer :: i_sst = -1
integer :: j_sst = -1
logical :: forecast_mode = .false.
logical :: use_ncep_sst  = .false.
real, allocatable, dimension(:,:) ::  sst_ncep, sst_anom
#else
use amip_interp_mod, only: i_sst, j_sst, sst_ncep, sst_anom, &
                           forecast_mode, use_ncep_sst
#endif

public i_sst, j_sst, sst_ncep, sst_anom, forecast_mode, use_ncep_sst

! Note: Added by Linjiong Zhou to use EC SST dataset for forcast
integer :: time_len
character(len=64) :: fn_ec_sst = "ec_sst_data.nc"
character(len=64) :: fn_ec_time = "INPUT/ec_sst.time.nc"
real, allocatable, dimension(:,:,:) :: ec_sst
type (time_type), allocatable, dimension(:) :: sst_time

public :: load_ec_sst, get_ec_sst

contains

! Note: Added by Linjiong Zhou to use EC SST dataset for forcast
subroutine load_ec_sst(Atm)

    implicit none

    type(fv_atmos_type), intent(in), target :: Atm
    type (restart_file_type) :: SST_restart
    type (restart_file_type) :: time_restart
    integer :: id_res
    integer :: is, ie, js, je, n
    integer :: siz(4)
    character(len = 256) :: units
    real, allocatable, dimension(:) :: time

    is = Atm%bd%is
    ie = Atm%bd%ie
    js = Atm%bd%js
    je = Atm%bd%je

    if (mpp_pe() .eq. mpp_root_pe()) then
        write(*,*) "EC SST dataset is used for forecast."
        write(*,*) "EC SST update frequency depends on the time interval of EC SST dataset."
    endif

    if (file_exist(fn_ec_time)) then

        call field_size(fn_ec_time, "tnew", siz)
 
        time_len = siz(1)
 
        if (time_len .lt. 1) then
            call mpp_error("external_sst_mod", "Invalid number of time records in file: "//trim(fn_ec_time), FATAL)
        endif

        allocate(time(time_len))
        allocate(sst_time(time_len))
 
        call read_data(fn_ec_time, "tnew", time, no_domain=.true.)
        call get_var_att_value(fn_ec_time, "tnew", "units", units)
 
        do n = 1, time_len
            sst_time(n) = get_cal_time(time(n), units, "gregorian")
        enddo
 
        allocate(ec_sst(is:ie,js:je,time_len))
 
        id_res = register_restart_field (SST_restart, fn_ec_sst, "sst", ec_sst, domain=Atm%domain)
        call restore_state (SST_restart)

    else

        call mpp_error("external_sst_mod","file: "//trim(fn_ec_time)//" does not exist.", FATAL)

    endif

end subroutine load_ec_sst

! Note: Added by Linjiong Zhou to use EC SST dataset for forcast
subroutine get_ec_sst(time, is, ie, js, je, sst)

    implicit none

    type (time_type), intent(in) :: Time

    integer, intent(in) :: is, ie, js, je

    real, intent(inout), dimension(is:ie,js:je) :: sst

    integer :: n
    integer :: year, month, day, hour, minute, second

    call get_date(Time, year, month, day, hour, minute, second)
    if (.not. (hour .eq. 0 .and. minute .eq. 0 .and. second .eq. 0)) return ! this requires daily data

    if (time < sst_time(1) .or. time > sst_time(time_len)) then
        call mpp_error("external_sst_mod","model time falls outside the range of SST data time.", FATAL)
    endif

    if (.not. (allocated(sst_time) .and. allocated(ec_sst))) then
        call mpp_error("external_sst_mod","neither sst_time nor ec_sst is allocated.", FATAL)
    endif

    do n = 1, time_len
        if (time == sst_time(n)) then
            sst = ec_sst(:,:,n)
            if (mpp_pe() .eq. mpp_root_pe()) then
                write(*,'(A,i4,A,i2,A,i2)') "Using EC SST on: ",year,"-",month,"-",day
            endif
            exit
        endif
    enddo

end subroutine get_ec_sst

end module external_sst_mod
