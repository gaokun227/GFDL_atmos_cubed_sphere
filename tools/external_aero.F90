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
! this module is designed to read 12 months climatology aerosol and
! interpolate to daily aerosol
! developer: linjiong zhou
! =======================================================================

module external_aero_mod

	use fv_arrays_mod, only: fv_atmos_type
	use mpp_mod, only: mpp_pe, mpp_root_pe

	public :: load_aero, read_aero, clean_aero

	integer :: nmon = 12, nlev = 72
	integer :: id_aero

	real, allocatable, dimension(:,:,:) :: aero_ps
	real, allocatable, dimension(:,:,:,:) :: aero_dp
	real, allocatable, dimension(:,:,:,:) :: aerosol

contains

! =======================================================================
! load aerosol 12 months climatological dataset

subroutine load_aero(Atm, fv_domain)

	use fms_io_mod, only: restart_file_type, register_restart_field
	use fms_io_mod, only: restore_state
	use fms_mod, only: file_exist, mpp_error, FATAL
	use mpp_domains_mod, only: domain2d
	use diag_manager_mod, only: register_static_field

	implicit none

	type(fv_atmos_type), intent(in), target :: Atm
	type(domain2d), intent(in) :: fv_domain
	type(restart_file_type) :: aero_restart

	integer :: is, ie, js, je
	integer :: id_res

	character(len=64) :: file_name = "MERRA2_400.inst3_3d_aer_Nv.climatology.nc"

	is = Atm%bd%is
	ie = Atm%bd%ie
	js = Atm%bd%js
	je = Atm%bd%je

	if (mpp_pe() .eq. mpp_root_pe()) then
		write(*,*) "aerosol 12 months climatological dataset is used for forecast."
	endif

	if (file_exist(file_name)) then
		if (.not. allocated(aero_ps)) allocate(aero_ps(is:ie,js:je,nmon))
		if (.not. allocated(aero_dp)) allocate(aero_dp(is:ie,js:je,nlev,nmon))
		if (.not. allocated(aerosol)) allocate(aerosol(is:ie,js:je,nlev,nmon))
		id_res = register_restart_field(aero_restart,file_name,"PS",aero_ps,domain=fv_domain)
		id_res = register_restart_field(aero_restart,file_name,"DELP",aero_dp,domain=fv_domain)
		id_res = register_restart_field(aero_restart,file_name,"SO4",aerosol,domain=fv_domain)
		call restore_state(aero_restart)
	else
		call mpp_error("external_aero_mod","file: "//trim(file_name)//" does not exist.",FATAL)
	endif

	id_aero = register_static_field('dynamics','aerosol',Atm%atmos_axes(1:2),'none','none')

end subroutine load_aero

! =======================================================================
! read aerosol climatological dataset

subroutine read_aero(Atm, Time)

	use constants_mod, only: grav
	use diag_manager_mod, only: send_data
	use time_manager_mod, only: time_type

	implicit none

	type(fv_atmos_type), intent(in), target :: Atm
	type(time_type), intent(in) :: Time

	integer :: k, n
	integer :: is, ie, js, je

	real, allocatable, dimension(:,:) :: vi_aero

	logical :: used

	is = Atm%bd%is
	ie = Atm%bd%ie
	js = Atm%bd%js
	je = Atm%bd%je

	if (.not. allocated(vi_aero)) allocate(vi_aero(is:ie,js:je))

	vi_aero = 0.0
	do n = 1, nmon
		do k = 1, nlev
			vi_aero = vi_aero+aerosol(:,:,k,n)*aero_dp(:,:,k,n)/grav*1.e6
		enddo
	enddo
	vi_aero = vi_aero/nmon

	if (id_aero > 0) used = send_data(id_aero,vi_aero,Time)

	if (allocated(vi_aero)) deallocate(vi_aero)

end subroutine read_aero

! =======================================================================
! clean aerosol climatological dataset

subroutine clean_aero()

	implicit none

	if (allocated(aero_ps)) deallocate(aero_ps)
	if (allocated(aero_dp)) deallocate(aero_dp)
	if (allocated(aerosol)) deallocate(aerosol)

end subroutine clean_aero

end module external_aero_mod
