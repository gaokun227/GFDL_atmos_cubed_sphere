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

module sa_3d_tke_mod

  use fv_arrays_mod,      only: fv_grid_bounds_type, fv_grid_type
  !use nh_utils_mod,       only: edge_profile1  
  !use fv_mp_mod,          only: is_master ! KGao - debug

  implicit none
  private
  public :: cal_3d_tke_budget, cal_3d_tke_budget_agrid 

contains

! =======================================================================
! cal_3d_tke_budget :: calculate 3D TKE shear production
! - this version is created based on Ping Zhu's anisotropic 
!   turbulence shear production formulation (REF to be added) and also
!   his code for 3D TKE shear production  
! - in this verion we use D-grid wind, instead of A-grid wind 
! - the algorithms are updated to be consistent with FV3
! - created by Kun Gao (kun.gao@noaa.gov)
! =======================================================================

  subroutine cal_3d_tke_budget(u, v, ua, va, w, tke,  &
                   delz, npz, ak, bk, gridstruct, bd, &
                   deform_1h, deform_1v)

    integer, intent(in) :: npz
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(fv_grid_type), intent(IN), target :: gridstruct
    real, intent(in) ::     ak(npz+1), bk(npz+1)
    real, intent(in) ::      u(bd%isd:bd%ied,  bd%jsd:bd%jed+1,npz)
    real, intent(in) ::      v(bd%isd:bd%ied+1,bd%jsd:bd%jed,  npz)
    real, intent(in) ::     ua(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    real, intent(in) ::     va(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    real, intent(in) ::      w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    real, intent(in) ::    tke(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    real, intent(in) ::   delz(bd%is:bd%ie,   bd%js:bd%je,   npz)

    real, intent(out) ::    deform_1h(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    real, intent(out) ::    deform_1v(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    !real, intent(out) ::    deform_2(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    !real, intent(out) ::    scl(bd%isd:bd%ied, bd%jsd:bd%jed)

! local

    real, pointer, dimension(:,:) :: dx, dy, rarea, dxc, dyc, rarea_c, dxa, dya
    integer :: is, ie, js, je, isd, ied, jsd, jed, i, j, k
    real:: dudx(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dudy(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dvdx(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dvdy(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dwdx(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dwdy(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dudz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dvdz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dwdz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: ut(bd%isd:bd%ied+1,bd%jsd:bd%jed)
    real:: vt(bd%isd:bd%ied,  bd%jsd:bd%jed+1)
    real:: wk1(bd%isd:bd%ied,bd%jsd:bd%jed) 
    real:: wk2(bd%isd:bd%ied,bd%jsd:bd%jed)
    real:: u_e(bd%is:bd%ie,npz+1)
    real:: v_e(bd%is:bd%ie,npz+1)
    real:: w_e(bd%is:bd%ie,npz+1)
    
    !real :: tke_1(bd%isd:bd%ied+2,  bd%jsd:bd%jed+2)
    !real :: dedy_1(bd%isd:bd%ied+1,  bd%jsd:bd%jed+1)
    !real :: dedy_2(bd%isd:bd%ied,  bd%jsd:bd%jed,npz)
    !real :: dedx_1(bd%isd:bd%ied+1,  bd%jsd:bd%jed+1)
    !real :: dedx_2(bd%isd:bd%ied,  bd%jsd:bd%jed,npz)
    !real :: tke_2(npz)
    !real :: dedz_1(npz)
    !real :: dedz_2(bd%isd:bd%ied,  bd%jsd:bd%jed,npz)

    real :: dp_ref(npz)
    real :: tmp, tkemax, dz
    integer:: l_tkemax,kscl

    real, parameter:: p1 =  7./12.
    real, parameter:: p2 = -1./12.

    dxc => gridstruct%dxc
    dyc => gridstruct%dyc
    dx  => gridstruct%dx
    dy  => gridstruct%dy
    dxa  => gridstruct%dxa
    dya  => gridstruct%dya
    rarea   => gridstruct%rarea
    rarea_c => gridstruct%rarea_c

    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd  = bd%isd
    ied  = bd%ied
    jsd  = bd%jsd
    jed  = bd%jed

!===========================================================
! Calculate deform_1h and deform_1v
!
!      shr_prod = kh * deform_1h + kv * deform_1v
!
!     deform_1h = (2*(du/dx**2 + dv/dy**2)
!                  + (du/dy + dv/dx)**2
!                  + (dw/dx**2 + dw/dy**2)
!                  + (dw/dx*du/dz + dw/dy*dv/dz))
!
!     deform_1v = (2*dw/dz**2 
!                  + (du/dz**2 + dv/dz**2)
!                  + (dw/dx*du/dz + dw/dy*dv/dz))
!===========================================================

!$OMP parallel do default(none) shared(npz,is,ie,js,je,u,v,w, &
!$OMP                                  dx,dy,rarea,dxc,dyc,rarea_c,dxa,dya, &
!$OMP                                  dudx,dudy,dvdx,dvdy,dwdx,dwdy) &
!$OMP                          private(ut,vt,wk1,wk2)
   do k=1,npz

!---------------------------------------------
! Calculate du/dx and dv/dy:
!  get corner values first, and 
!  then average the 4 corners to get cell-mean
!---------------------------------------------

      do j=js-1,je+2
         do i=is-2,ie+2
            ut(i,j) = u(i,j,k)*dyc(i,j)
         enddo
      enddo
      do j=js-2,je+2
         do i=is-1,ie+2
            vt(i,j) = v(i,j,k)*dxc(i,j)
         enddo
      enddo
      do j=js,je+1
         do i=is,ie+1
            wk1(i,j) = rarea_c(i,j)*( -ut(i-1,j) + ut(i,j) ) ! du/dx at corners
            wk2(i,j) = rarea_c(i,j)*( -vt(i,j-1) + vt(i,j) ) ! dv/dy at corners
         enddo
      enddo
! Fix the corners if grid_type /= 4
      do j=js,je
         do i=is,ie
            ! mean of 4 corners
            dudx(i,j,k) = 0.25*(wk1(i,j) + wk1(i,j+1) + wk1(i+1,j) + wk1(i+1,j+1))
            dvdy(i,j,k) = 0.25*(wk2(i,j) + wk2(i,j+1) + wk2(i+1,j) + wk2(i+1,j+1))
         enddo
      enddo

!-------------------------------------
! Calculate du/dy and dv/dx
!-------------------------------------

      do j=js,je+1
         do i=is,ie
            vt(i,j) = u(i,j,k)*dx(i,j)
         enddo
      enddo
      do j=js,je
         do i=is,ie+1
            ut(i,j) = v(i,j,k)*dy(i,j)
         enddo
      enddo

      do j=js,je
         do i=is,ie
            dudy(i,j,k) = rarea(i,j)*( -vt(i,j) + vt(i,j+1) )
            dvdx(i,j,k) = rarea(i,j)*( -ut(i,j) + ut(i+1,j) )
         enddo
      enddo

!-------------------------------------------
! Calculate dw/dx and dw/dy:
!  get edge values of w first, and then 
!  use the two edge values to get cell-mean 
!  dw/dx and dw/dy
!-------------------------------------------

      do j=js,je
         do i=is,ie
           ! using the 1d ppm method to get w values at cell edges
           ! q(i+0.5) = p1*(q(i-1)+q(i)) + p2*(q(i-2)+q(i+1))
           wk1(i+1,j) = p1*(w(i-1,j,k)+w(i,j,k)) + p2*(w(i-2,j,k)+w(i+1,j,k))
           wk2(i,j+1) = p1*(w(i,j-1,k)+w(i,j,k)) + p2*(w(i,j-2,k)+w(i,j+1,k))
         enddo
      enddo

      do j=js,je
         do i=is,ie
            dwdx(i,j,k) = 1./dxa(i,j) * (-wk1(i,j) + wk1(i+1,j))
            dwdy(i,j,k) = 1./dya(i,j) * (-wk2(i,j) + wk2(i,j+1))
         enddo
      enddo

   enddo   !z loop

!-------------------------------------
! get du/dz, dv/dz and dw/dz
!-------------------------------------

! TODO: use $OMP 
   do k=1,npz
      dp_ref(k) = ak(k+1)-ak(k) + (bk(k+1)-bk(k))*1.E5
   enddo

   do j=js,je
      call edge_profile1(ua(is:ie,j,:), u_e, is, ie, npz, dp_ref, 1)
      call edge_profile1(va(is:ie,j,:), v_e, is, ie, npz, dp_ref, 1)
      call edge_profile1(w(is:ie,j,:),  w_e, is, ie, npz, dp_ref, 1)
      do k=1,npz
         do i=is,ie
            dz = delz(i,j,k) !zh(i,j,k) - zh(i,j,k+1) 
            dudz(i,j,k) = (u_e(i,k)-u_e(i,k+1))/dz
            dvdz(i,j,k) = (v_e(i,k)-v_e(i,k+1))/dz
            dwdz(i,j,k) = (w_e(i,k)-w_e(i,k+1))/dz
         enddo
      enddo
   enddo

!-----------------------------------------------
! get deform_1h and deform_1v based on all terms
!-----------------------------------------------

!$OMP parallel do default(none) shared(npz,is,ie,js,je, &
!$OMP                                  deform_1h, deform_1v, &
!$OMP                                  dudx,dudy,dudz,dvdx, &
!$OMP                                  dvdy,dvdz,dwdx,dwdy,dwdz) &
!$OMP                          private(tmp)
   do k=1,npz
       do j=js,je
          do i=is,ie
             !deform_1(i,j,k)=2*(dudx(i,j,k)**2+dvdy(i,j,k)**2+  &
             !   dwdz(i,j,k)**2)+(dudy(i,j,k)+dvdx(i,j,k))**2+   &
             !   (dudz(i,j,k)+dwdx(i,j,k))**2+                   &
             !   (dvdz(i,j,k)+dwdy(i,j,k))**2

             tmp = dwdx(i,j,k)*dudz(i,j,k)+dwdy(i,j,k)*dvdz(i,j,k)

             deform_1h(i,j,k)=                            &
               (2*(dudx(i,j,k)**2 + dvdy(i,j,k)**2)       &
                + (dudy(i,j,k) + dvdx(i,j,k))**2+         &
                + dwdx(i,j,k)**2 + dwdy(i,j,k)**2 + tmp)


             deform_1v(i,j,k)=                            &
               (2*dwdz(i,j,k)**2                          &
                + dudz(i,j,k)**2 + dvdz(i,j,k)**2 + tmp)

             ! apply limiters
             deform_1h(i,j,k) = min( 1.e3, deform_1h(i,j,k) )
             !deform_1v(i,j,k) = min( 1.e3, deform_1v(i,j,k) )
          enddo
       enddo
   enddo

   ! KGao: debug code
   !if (is_master())  then
   !   write(*,*) 'KGao debug - max deform_1 ', maxval(abs(deform_1(is:ie,js:je,:)))
   !   write(*,*) 'KGao debug - max dudx ', maxval(abs(dudx(is:ie,js:je,:)))
   !   write(*,*) 'KGao debug - max dudz ', maxval(abs(dudz(is:ie,js:je,:)))
   !   write(*,*) 'KGao debug - max dz ', maxval(abs(delz))
   !   write(*,*) 'KGao debug - max dx ', maxval(abs(dx))
   !endif

  end subroutine cal_3d_tke_budget

! =======================================================================
! cal_3d_tke_budget_agrid :: Ping Zhu's scheme for 3D TKE budget terms
! This version calculate the shear production terms based on A-grid wind 
! =======================================================================

  subroutine cal_3d_tke_budget_agrid(ua, va, w, tke, delz, npz, ak, bk, gridstruct, bd, &
                   deform_1) !, deform_2)


    integer, intent(in) :: npz
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(fv_grid_type), intent(IN), target :: gridstruct
    real, intent(in) ::     ak(npz+1), bk(npz+1)
    real, intent(in) ::     ua(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    real, intent(in) ::     va(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    real, intent(in) ::      w(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    real, intent(in) ::    tke(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    real, intent(in) ::   delz(bd%is:bd%ie,   bd%js:bd%je,   npz)

    real, intent(out) ::    deform_1(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    !real, intent(out) ::    deform_2(bd%isd:bd%ied, bd%jsd:bd%jed, npz)
    !real, intent(out) ::    scl(bd%isd:bd%ied, bd%jsd:bd%jed)

! local

    real, pointer, dimension(:,:) :: dx, dy, rarea
    integer :: is, ie, js, je, isd, ied, jsd, jed, i, j, k
    real:: dudx(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dudy(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dvdx(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dvdy(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dwdx(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dwdy(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dudz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dvdz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: dwdz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
    real:: ut(bd%isd:bd%ied+1,bd%jsd:bd%jed)
    real:: vt(bd%isd:bd%ied,  bd%jsd:bd%jed+1)
    real:: u_e(bd%is:bd%ie,npz+1)
    real:: v_e(bd%is:bd%ie,npz+1)
    real:: w_e(bd%is:bd%ie,npz+1)
    
    real :: tke_1(bd%isd:bd%ied+2,  bd%jsd:bd%jed+2)
    real :: dedy_1(bd%isd:bd%ied+1,  bd%jsd:bd%jed+1)
    real :: dedy_2(bd%isd:bd%ied,  bd%jsd:bd%jed,npz)
    real :: dedx_1(bd%isd:bd%ied+1,  bd%jsd:bd%jed+1)
    real :: dedx_2(bd%isd:bd%ied,  bd%jsd:bd%jed,npz)
    real :: tke_2(npz)
    real :: dedz_1(npz)
    real :: dedz_2(bd%isd:bd%ied,  bd%jsd:bd%jed,npz)

    real :: dp_ref(npz)
    real :: tkemax, dz
    integer:: l_tkemax,kscl

    dx  => gridstruct%dx
    dy  => gridstruct%dy
    rarea   => gridstruct%rarea

    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd  = bd%isd
    ied  = bd%ied
    jsd  = bd%jsd
    jed  = bd%jed

!===========================================================
! Calculate deform_1
! deform_1 = 2*(du/dx**2 + dv/dy**2 + dw/dz**2)
!          + (du/dy + dv/dx) ** 2
!          + (du/dz + dw/dx) ** 2
!          + (dv/dz + dw/dy) ** 2
!===========================================================

! KGao: make ut and vt private as suggested by Lucas

!$OMP parallel do default(none) shared(npz,is,ie,js,je,ua,va,w,dx,dy,rarea, &
!$OMP                           dudx,dudy,dvdx,dvdy,dwdx,dwdy)              &
!$OMP                           private(ut,vt)
   do k=1,npz

!-------------------------------------
! get du/dy and dv/dx

       do j=js,je+1
          do i=is,ie
             vt(i,j) = ua(i,j,k)*dx(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             ut(i,j) = va(i,j,k)*dy(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie
             dudy(i,j,k) = rarea(i,j)*(vt(i,j+1)-vt(i,j))
             dvdx(i,j,k) = rarea(i,j)*(ut(i+1,j)-ut(i,j))
          enddo
       enddo

!-------------------------------------
! get du/dx and dv/dy

       do j=js,je
          do i=is,ie+1
             ut(i,j) = ua(i,j,k)*dy(i,j)
          enddo
       enddo
       do j=js,je+1
          do i=is,ie
             vt(i,j) = va(i,j,k)*dx(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie
             dudx(i,j,k) = rarea(i,j)*(ut(i+1,j)-ut(i,j))
             dvdy(i,j,k) = rarea(i,j)*(vt(i,j+1)-vt(i,j))
          enddo
       enddo

!-------------------------------------
! get dw/dx and dw/dy

       do j=js,je
          do i=is,ie+1
             ut(i,j) = w(i,j,k)*dy(i,j)
          enddo
       enddo
       do j=js,je+1
          do i=is,ie
             vt(i,j) = w(i,j,k)*dx(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie
             dwdx(i,j,k) = rarea(i,j)*(ut(i+1,j)-ut(i,j))
             dwdy(i,j,k) = rarea(i,j)*(vt(i,j+1)-vt(i,j))
          enddo
       enddo

   enddo   !z loop

!-------------------------------------
! get du/dz, dv/dz and dw/dz

! KGao: use Lucas's method as in compute_dudz()
   do k=1,npz
      dp_ref(k) = ak(k+1)-ak(k) + (bk(k+1)-bk(k))*1.E5
   enddo

   do j=js,je
      call edge_profile1(ua(is:ie,j,:), u_e, is, ie, npz, dp_ref, 1)
      call edge_profile1(va(is:ie,j,:), v_e, is, ie, npz, dp_ref, 1)
      call edge_profile1(w(is:ie,j,:),  w_e, is, ie, npz, dp_ref, 1)
      do k=1,npz
         do i=is,ie
            dz = delz(i,j,k) !zh(i,j,k) - zh(i,j,k+1) ! can use delz too
            dudz(i,j,k) = (u_e(i,k)-u_e(i,k+1))/dz
            dvdz(i,j,k) = (v_e(i,k)-v_e(i,k+1))/dz
            dwdz(i,j,k) = (w_e(i,k)-w_e(i,k+1))/dz
         enddo
      enddo
   enddo

!!$OMP parallel do default(none) shared(npz,is,ie,js,je,ua,va,w,delz, &
!!$OMP                          dudz,dvdz,dwdz)

!   do j=js,je
!       do i=is,ie
!          do k=1,npz-1
!             dudz(i,j,k)=(ua(i,j,k+1)-ua(i,j,k))/delz(i,j,k) ! KGao: this is problematic
!             dvdz(i,j,k)=(va(i,j,k+1)-va(i,j,k))/delz(i,j,k)
!             dwdz(i,j,k)=(w(i,j,k+1)-w(i,j,k))/delz(i,j,k)
!          enddo
!          dudz(i,j,npz)=0.0
!          dvdz(i,j,npz)=0.0
!          dwdz(i,j,npz)=0.0
!       enddo
!   enddo

!-------------------------------------
! get deform_1 based on all terms

!$OMP parallel do default(none) shared(npz,is,ie,js,je,deform_1, &
!$OMP        dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
   do k=1,npz
       do j=js,je
          do i=is,ie
             deform_1(i,j,k)=2*(dudx(i,j,k)**2+dvdy(i,j,k)**2+  &
                dwdz(i,j,k)**2)+(dudy(i,j,k)+dvdx(i,j,k))**2+   &
                (dudz(i,j,k)+dwdx(i,j,k))**2+                     &
                (dvdz(i,j,k)+dwdy(i,j,k))**2
             ! KGao - enforce a limiter
             deform_1(i,j,k) = min( 1.e3, deform_1(i,j,k) )
          enddo
       enddo
   enddo

   !if (is_master())  then
   !   write(*,*) 'KGao debug - max deform_1 ', maxval(abs(deform_1(is:ie,js:je,:)))
   !   write(*,*) 'KGao debug - max dudx ', maxval(abs(dudx(is:ie,js:je,:)))
   !   write(*,*) 'KGao debug - max dudz ', maxval(abs(dudz(is:ie,js:je,:)))
   !   write(*,*) 'KGao debug - max dz ', maxval(abs(delz))
   !   write(*,*) 'KGao debug - max dx ', maxval(abs(dx))
   !endif

!===========================================================
! Calculate deform_2
! deform_2 = (d^2e/dx^2 + d^2e/dy^2 + d^2e/dz^2)
! where e is tke
!===========================================================

!find scale of energy containin eddies using TKE

!!$OMP parallel do default(none) shared(npz,is,ie,js,je,tke, &
!!$OMP     tkemax,l_tkemax,kscl,scl,zh)
!   do j=js,je
!      do i=is,ie
!         !scl(i,j)=1000.0
!         l_tkemax=10
!         kscl=10
!         tkemax=0.0
!         do k=1,npz
!            tkemax=max(tkemax,tke(i,j,k))
!         enddo
!         do k=1,npz
!            if (abs(tke(i,j,k)-tkemax)/tkemax .lt. 1.0e-9) then
!               l_tkemax=k
!            endif
!         enddo
!         do k=l_tkemax,npz
!            if (tke(i,j,k)-0.5*tkemax .gt. 0.0) then
!               kscl=k
!            endif
!         enddo
!         kscl=min(kscl,npz-10)
!         !scl(i,j)=zh(i,j,kscl+1) ! KGao: use zh, not gz; zh is needed as input !!!
!      enddo
!   enddo

! KGao: make the 2d temporay arrays private - as suggested by Lucas

!!$OMP parallel do default(none) shared(npz,is,ie,js,je,ua,va,w,dx,dy,rarea, &
!!$OMP                                  tke,dedy_2,dedx_2)                &
!!$OMP                           private(ut,vt,tke_1,dedx_1,dedy_1)
!   do k=1,npz

!-------------------------------------
! get d^2e/dy^2

!       do j=js,je+2
!          do i=is,ie+2
!             tke_1(i,j)=tke(i,j,k)
!          enddo
!       enddo
!       do j=js,je+2
!          do i=is,ie
!             vt(i,j)=tke_1(i,j)*dx(i,j)
!          enddo
!       enddo
!       do j=js,je+1
!          do i=is,ie
!             dedy_1(i,j)=rarea(i,j)*(vt(i,j+1)-vt(i,j))
!          enddo
!       enddo
!       do j=js,je+1
!          do i=is,ie
!             vt(i,j)=dedy_1(i,j)*dx(i,j)
!          enddo
!       enddo
!       do j=js,je
!          do i=is,ie
!             dedy_2(i,j,k)=rarea(i,j)*(vt(i,j+1)-vt(i,j))
!          enddo
!       enddo

!-------------------------------------
! get d^2e/dx^2

!       do j=js,je
!          do i=is,ie+2
!             ut(i,j)=tke_1(i,j)*dy(i,j)
!          enddo
!       enddo
!       do j=js,je
!          do i=is,ie+1
!             dedx_1(i,j)=rarea(i,j)*(ut(i+1,j)-ut(i,j))
!          enddo
!       enddo
!       do j=js,je
!          do i=is,ie+1
!             ut(i,j)=dedx_1(i,j)*dy(i,j)
!          enddo
!       enddo
!       do j=js,je
!          do i=is,ie
!             dedx_2(i,j,k)=rarea(i,j)*(ut(i+1,j)-ut(i,j))
!          enddo
!       enddo
!   enddo ! z loop

!-------------------------------------
! get d^2e/dz^2

!!$OMP parallel do default(none) shared(npz,is,ie,js,je,q,delz, &
!!$OMP     tke_2,dedz_1,dedz_2)
!   do j=js,je
!       do i=is,ie
!          do k=1,npz
!             tke_2(k)=tke(i,j,k)
!          enddo
!          do k=1,npz-1
!             dedz_1(k)=(tke_2(k+1)-tke_2(k))/delz(i,j,k)
!          enddo
!          do k=1,npz-2
!             dedz_2(i,j,k)=(dedz_1(k+1)-dedz_1(k))/delz(i,j,k)
!          enddo
!          dedz_2(i,j,npz-1)=0.0
!          dedz_2(i,j,npz)=0.0
!       enddo
!   enddo

!-------------------------------------
! get deform_2 based on all terms

!!$OMP parallel do default(none) shared(npz,is,ie,js,je,deform_2, &
!!$OMP                                  dedx_2,dedy_2)
!   do k=1,npz
!       do j=js,je
!          do i=is,ie
!   !          deform_2(i,j,k)=dedx_2(i,j,k)+dedy_2(i,j,k)+dedz_2(i,j,k)
!             deform_2(i,j,k)=dedx_2(i,j,k)+dedy_2(i,j,k)
!          enddo
!       enddo
!   enddo

  end subroutine cal_3d_tke_budget_agrid

  subroutine edge_profile1(q1, q1e, i1, i2, km, dp0, limiter)
! Edge profiles for a single scalar quantity
 integer, intent(in):: i1, i2
 integer, intent(in):: km
 integer, intent(in):: limiter
 real, intent(in):: dp0(km)
 real, intent(in),  dimension(i1:i2,km):: q1
 real, intent(out), dimension(i1:i2,km+1):: q1e
!-----------------------------------------------------------------------
 real, dimension(i1:i2,km+1):: qe1, gam  ! edge values
 real  gak(km)
 real  bet, r2o3, r4o3
 real  g0, gk, xt1, xt2, a_bot
 integer i, k

! Assuming grid varying in vertical only
   g0 = dp0(2) / dp0(1)
  xt1 = 2.*g0*(g0+1. )
  bet =    g0*(g0+0.5)
  do i=i1,i2
      qe1(i,1) = ( xt1*q1(i,1) + q1(i,2) ) / bet
      gam(i,1) = ( 1. + g0*(g0+1.5) ) / bet
  enddo

  do k=2,km
     gk = dp0(k-1) / dp0(k)
     do i=i1,i2
             bet =  2. + 2.*gk - gam(i,k-1)
        qe1(i,k) = ( 3.*(q1(i,k-1)+gk*q1(i,k)) - qe1(i,k-1) ) / bet
        gam(i,k) = gk / bet
     enddo
  enddo

  a_bot = 1. + gk*(gk+1.5)
    xt1 =   2.*gk*(gk+1.)
  do i=i1,i2
             xt2 = gk*(gk+0.5) - a_bot*gam(i,km)
     qe1(i,km+1) = ( xt1*q1(i,km) + q1(i,km-1) - a_bot*qe1(i,km) ) / xt2
  enddo

  do k=km,1,-1
     do i=i1,i2
        qe1(i,k) = qe1(i,k) - gam(i,k)*qe1(i,k+1)
     enddo
  enddo

!------------------
! Apply constraints
!------------------
    if ( limiter/=0 ) then   ! limit the top & bottom winds
         do i=i1,i2
! Top
            if ( q1(i,1)*qe1(i,1) < 0. ) qe1(i,1) = 0.
! Surface:
            if ( q1(i,km)*qe1(i,km+1) < 0. ) qe1(i,km+1) = 0.
         enddo
    endif

    do k=1,km+1
       do i=i1,i2
          q1e(i,k) = qe1(i,k)
       enddo
    enddo

  end subroutine edge_profile1

end module sa_3d_tke_mod
