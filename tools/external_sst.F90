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
  use fms_mod,             only: file_exist, mpp_error, FATAL, NOTE
  use mpp_mod,             only: mpp_pe, mpp_root_pe
  use fv_arrays_mod,       only: fv_atmos_type
  use time_manager_mod,    only: time_type, get_date, operator(==), &
                                 operator(>), operator(<)
  use get_cal_time_mod,    only: get_cal_time
  use data_override_mod,   only: data_override

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

public :: i_sst, j_sst, sst_ncep, sst_anom, forecast_mode, use_ncep_sst
public :: get_ec_sst

contains

! Here is a sample data_table that will enable reading in
! EC SSTs and sea-ice from an external file. If read_ec_sst = .true.
! then the 'sst' entry MUST be present. The 'ci' entry is optional.
!
!"ATM", "sst", "sst", "INPUT/ec_sst.nc", "bilinear", 1.0
!"ATM", "ci", "ci", "INPUT/ec_sst.nc", "bilinear", 1.0

  
subroutine get_ec_sst(time, is, ie, js, je, sst, ci)

    implicit none

    type (time_type), intent(in) :: Time

    integer, intent(in) :: is, ie, js, je
    real, intent(inout), dimension(is:ie,js:je) :: sst, ci
    logical :: used

    call data_override('ATM', 'sst', sst, time, override=used)
    if (.not. used) then
       call mpp_error(FATAL, " SST dataset not specified in data_table.")
    endif
    call data_override('ATM', 'ci', ci, time, override=used)
    if (.not. used) then
       call mpp_error(NOTE, " Sea ice fraction dataset not specified in data_table. No override will occur.")
       ci(:,:) = -999.
    endif
    
end subroutine get_ec_sst

end module external_sst_mod
