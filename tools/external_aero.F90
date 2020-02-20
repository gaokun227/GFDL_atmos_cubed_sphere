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

	use mpp_mod, only: mpp_pe, mpp_root_pe

	public :: load_aero, read_aero, clean_aero

	real, allocatable, dimension(:,:,:,:) :: aerosol

contains

! =======================================================================
! load aerosol 12 months climatological dataset

subroutine load_aero(Atm)

	use fms_io_mod, only: restart_file_type, register_restart_field
	use fms_io_mod, only: restore_state
	use fv_arrays_mod, only: fv_atmos_type
  use fms_mod, only: file_exist, mpp_error, FATAL

	implicit none

	type(fv_atmos_type), intent(in), target :: Atm(:)
	type(restart_file_type) :: aero_restart

	integer :: is, ie, js, je, n
	integer :: nmon = 12, nlev = 72
	integer :: id_res

	character(len=64) :: file_name = "MERRA2_400.inst3_3d_aer_Nv.climatology.nc"

	is = Atm(1)%bd%is
	ie = Atm(1)%bd%ie
	js = Atm(1)%bd%js
	je = Atm(1)%bd%je

	if (mpp_pe() .eq. mpp_root_pe()) then
		write(*,*) "aerosol 12 months climatological dataset is used for forecast."
	endif

	if (file_exist(file_name)) then
		if (.not. allocated(aerosol)) allocate(aerosol(is:ie,js:je,nlev,nmon))
		do n = 1, size(Atm(:))
			id_res = register_restart_field(aero_restart,file_name,"SO4",aerosol,domain=Atm(n)%domain)
			call restore_state(aero_restart)
		enddo
	else
		call mpp_error("external_aero_mod","file: "//trim(file_name)//" does not exist.",FATAL)
	endif

end subroutine load_aero

! =======================================================================
! read aerosol climatological dataset

subroutine read_aero()

	implicit none

end subroutine read_aero

! =======================================================================
! clean aerosol climatological dataset

subroutine clean_aero()

	implicit none

	if (allocated(aerosol)) deallocate(aerosol)

end subroutine clean_aero

end module external_aero_mod
