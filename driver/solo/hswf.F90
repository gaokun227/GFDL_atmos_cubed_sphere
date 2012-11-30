module hswf_mod

 use constants_mod,      only: grav, rdgas, rvgas, RADIAN, kappa, radius
 use fv_grid_utils_mod,  only: g_sum
 use mpp_domains_mod,    only: mpp_update_domains
 use fv_mp_mod,          only: domain
 use time_manager_mod,   only: time_type, get_date, get_time
 use diag_manager_mod,   only: send_data
 use fv_timing_mod,      only: timing_on, timing_off

#ifdef MARS_GCM
      use fv_current_grid_mod, only: tmars_initialized, tmars
      use fms_mod, only: file_exist, read_data, field_size
      use mpp_mod, only: mpp_error, FATAL
      use horiz_interp_mod, only: horiz_interp
#endif

      implicit none
!-----------------------------------------------------------------------
      logical :: rf_initialized = .false.

#ifdef MARS_GCM
      logical :: tmars_initialized = .false.
      real,  allocatable, dimension(:,:,:):: tmars
#endif


      private
      public :: Held_Suarez_Strat, Held_Suarez_Tend, age_of_air

!---- version number -----
      character(len=128) :: version = '$Id: hswf.F90,v 17.0.2.2.2.3.2.2 2012/04/30 17:08:46 Lucas.Harris Exp $'
      character(len=128) :: tagname = '$Name: siena_201211 $'

contains

!-----------------------------------------------------------------------
 subroutine Held_Suarez_Tend(npx, npy, npz, is, ie, js, je, ng, nq,   &
                              u, v, pt, q, pe, delp, peln, pkz, pdt,  &
                              ua, va, u_dt, v_dt, t_dt, q_dt, agrid,  &
                              delz, hydrostatic, ak, bk, ks,   &
                              strat, rayf, master, Time, time_total)

      integer, INTENT(IN   ) :: npx, npy, npz
      integer, INTENT(IN   ) :: is, ie, js, je, ng, nq
      logical, intent(IN)    :: hydrostatic
      real   , INTENT(IN   ) ::  delz(is:ie,js:je,npz)

      real   , INTENT(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
      real   , INTENT(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) ::    q(is-ng:ie+  ng,js-ng:je+  ng,npz, nq)
      real   , INTENT(INOUT) ::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
      real   , INTENT(INOUT) :: peln(is  :ie     ,1:npz+1,js   :je     )
      real   , INTENT(INOUT) ::  pkz(is  :ie     ,js   :je     ,1:npz)

      real   , INTENT(INOUT) ::   ua(is-ng:ie+ng,js-ng:je+ng,npz)
      real   , INTENT(INOUT) ::   va(is-ng:ie+ng,js-ng:je+ng,npz)

! Tendencies:
      real, INTENT(INOUT):: u_dt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: v_dt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
      real, INTENT(INOUT):: q_dt(is:ie,js:je,npz,nq)


      real   , INTENT(IN   ) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
      real   , INTENT(IN   ) :: ak(npz+1), bk(npz+1)
      integer, INTENT(IN   ) :: ks

      real   , INTENT(IN   ) :: pdt
      logical, INTENT(IN   ) :: strat, rayf, master

      type(time_type), intent(in) :: Time
      real, INTENT(IN), optional:: time_total

! Local
      real pref(npz)
      integer  i,j,k
      real  ty, tz, akap 
      real  p0, t0, sday, rkv, rka, rks, rkt, sigb, rsgb
      real  tmp
      real  ap0k, algpk
      real  tey, tez, fac, pw, sigl
      real  h0, dz
      real  dt_tropic
      real  rmr, rms
      real  relx, tau
      real  t_st, t_ms
      real  rdt, f1
      real  pc, c1

      real, allocatable, SAVE ::  rf(:)
      real :: frac(is-ng:ie+ng,js-ng:je+ng)
      real ::  pl(is:ie, js:je, 1:npz)
      real :: teq(is:ie, js:je, 1:npz)
      real rdg, rad_ratio, kf_day

      rdg = -rdgas / grav

      ty = 60.0
      tz = 10.0             ! Original value from H-S was 10.
      akap = 2./7.

      p0 = 100000.
      t0 = 200.
      h0 = 7.
      sday = 24.*3600.
      rdt = 1. / pdt

!--------------------------
      rad_ratio = radius / 6371.0e3

      kf_day = sday * rad_ratio
      rkv = pdt / kf_day
      rka = pdt      / (40.*kf_day)
      rks = pdt      / (4.0*kf_day)

!     rkv = pdt/sday
!     rka = pdt/ (40.*sday)      ! was 40 days
!     rks = pdt/ (4.*sday)       ! was 4 days
!--------------------------

! For strat-mesosphere
      t_ms = 10.*rad_ratio
      t_st = 40.*rad_ratio

      tau = (t_st - t_ms) / log(100.)
      rms = pdt/(t_ms*sday)
      rmr =  1./(1.+rms)

      sigb = 0.7
      rsgb = 1./(1.-sigb)
      ap0k = 1./p0**akap
      algpk = log(ap0k)

      do k=1,npz
         pref(k) = ak(k) + bk(k)*1.E5
      enddo

! Setup...
      if ( rayf .and. (.not. rf_initialized) ) then
          allocate( rf(npz) )
          c1 = 1. / (36.*3600)
          pc = 1.
          if(master) write(6,*) 'HSWF Forcing ...' 
          do k=1,ks
             tmp = (ak(k+1)-ak(k))/log(ak(k+1)/ak(k))
             rf(k) = c1*(1.+tanh(log10(pc/tmp)))
             if(master) write(6,*) k, 0.01*tmp, 1./(rf(k)*sday)
             rf(k) = 1./(1.+pdt*rf(k))
          enddo
          if(master) write(6,*) ' '
          rf_initialized = .true.
      endif

        do k=1,npz
           do j=js,je
              do i=is,ie
                 pl(i,j,k) = delp(i,j,k) / ( peln(i,k+1,j)-peln(i,k,j))
              enddo
           enddo 
        enddo

#ifdef TEST_HYDRO
        if ( .not. hydrostatic ) then
           do k=1,npz
              do j=js,je
                 do i=is,ie
                    pkz(i,j,k) = pl(i,j,k)**akap
                 enddo
              enddo 
           enddo
         endif
#endif
        
! Temperature forcing...
        do k=npz,1,-1
           do j=js,je
              do i=is,ie
                 tey = ap0k*( 315.0 - ty*SIN(agrid(i,j,2))*SIN(agrid(i,j,2)) )
                 tez = tz*( ap0k/akap )*COS(agrid(i,j,2))*COS(agrid(i,j,2)) 

                 if (strat .and. pl(i,j,k) < 10000.    &
                           .and. pl(i,j,k) > 100.  )  then
                   dz = h0 * log(pl(i,j,k+1)/pl(i,j,k))
!
! Lapse rate above tropic stratopause is 2.25 deg/km
! Relaxation time is t_st days at 100 mb (as H-S) and gradually
! decreases to t_ms Days at and above the stratopause
!
                   relx =  t_ms + tau*log(0.01*pl(i,j,k))
                   relx = pdt/(relx*sday)
                   dt_tropic = 2.25*COS(agrid(i,j,2)) * dz
                   teq(i,j,k)  =  teq(i,j,k+1) + dt_tropic
                   t_dt(i,j,k) = t_dt(i,j,k) + relx*(teq(i,j,k)-pt(i,j,k))/(1.+relx) * rdt
                 elseif (strat .and. pl(i,j,k) <= 100.)  then
!
! Mesosphere: defined as the region above 1 mb
!
                   dz = h0 * log(pl(i,j,k+1)/pl(i,j,k))
                   dt_tropic = -2.25*COS(agrid(i,j,2)) * dz
                   teq(i,j,k) = teq(i,j,k+1) + dt_tropic
                   t_dt(i,j,k) = t_dt(i,j,k) + ((pt(i,j,k)+rms*teq(i,j,k))*rmr - pt(i,j,k))*rdt
                 else

! Trop:  strictly Held-Suarez

                   sigl = pl(i,j,k)/pe(i,npz+1,j)
                   f1 = max(0., (sigl-sigb) * rsgb )
                   teq(i,j,k) = tey - tez*(log(pkz(i,j,k))+algpk)
                   teq(i,j,k) = max(t0, teq(i,j,k)*pkz(i,j,k))
                   rkt = rka + (rks-rka)*f1*(COS(agrid(i,j,2))**4.0)
                   t_dt(i,j,k) = t_dt(i,j,k) + rkt*(teq(i,j,k)-pt(i,j,k))/(1.+rkt) * rdt
                 endif
              enddo     !i-loop
           enddo     !j-loop
        enddo     !k-loop

#ifdef TEST_HYDRO
        if ( .not. hydrostatic ) then
!------------------------------------------------------------------------
! The nonhydrostatic pressure changes if there is heating (under constant
! volume and mass is locally conserved).
!------------------------------------------------------------------------
           do k=1,npz
              do j=js,je
                 do i=is,ie
                      pkz(i,j,k) = (rdg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k))**akap
                 enddo
              enddo 
           enddo
        endif
#endif

!!! NOTE: Not sure what to do in the boundary for a nested grid
        u_dt = 0.
        v_dt = 0.
! Velocity dissipation damping
      do 2000 k=1,npz

        if (rayf .and. k <= ks) then
! Apply Rayleigh friction
          do j=js,je
             do i=is,ie
                u_dt(i,j,k) = u_dt(i,j,k) + ua(i,j,k)*(rf(k) - 1.) * rdt
                v_dt(i,j,k) = v_dt(i,j,k) + va(i,j,k)*(rf(k) - 1.) * rdt
             enddo
          enddo
        else
! Surface Rayleigh friction according to Held-Suarez
           do j=js,je
              do i=is,ie
                 sigl = 0.5*(pe(i,k,j)+pe(i,k+1,j)) / pe(i,npz+1,j)
                 frac(i,j) = max(0., (sigl-sigb)*rsgb ) * rkv
                 if (frac(i,j) > 0.) then
                     u_dt(i,j,k) = u_dt(i,j,k) - ua(i,j,k)*frac(i,j)/(1.+frac(i,j)) * rdt
                     v_dt(i,j,k) = v_dt(i,j,k) - va(i,j,k)*frac(i,j)/(1.+frac(i,j)) * rdt
                 endif
               enddo
           enddo
        endif

2000  continue

#ifdef DO_AGE
      if( nq/=0 )     &
      call age_of_air(is, ie, js, je, npz, ng, time_total, pe, q(is-ng,js-ng,1,nq))
#endif


 end subroutine Held_Suarez_Tend

 subroutine Held_Suarez_Strat(npx, npy, npz, is, ie, js, je, ng, nq,  &
                              u, v, pt, q, pe, delp, peln, pkz, pdt,  &
                              ua, va, agrid, ak, bk, ks, strat,  &
                              rayf, master, Time, time_total)

      integer, INTENT(IN   ) :: npx, npy, npz
      integer, INTENT(IN   ) :: is, ie, js, je, ng, nq

	      real   , INTENT(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
	      real   , INTENT(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
	      real   , INTENT(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
	      real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
	      real   , INTENT(INOUT) ::    q(is-ng:ie+  ng,js-ng:je+  ng,npz, nq)
	      real   , INTENT(INOUT) ::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
	      real   , INTENT(INOUT) :: peln(is   :ie     ,1:npz+1,js   :je     )
	      real   , INTENT(INOUT) ::  pkz(is   :ie     ,js   :je     ,1:npz)

	      real   , INTENT(INOUT) ::   ua(is-ng:ie+ng,js-ng:je+ng,npz)
	      real   , INTENT(INOUT) ::   va(is-ng:ie+ng,js-ng:je+ng,npz)

	      real   , INTENT(IN   ) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
	      real   , INTENT(IN   ) :: ak(npz+1), bk(npz+1)
              integer, INTENT(IN   ) :: ks

	      real   , INTENT(IN   ) :: pdt
	      logical, INTENT(IN   ) :: strat, rayf, master

              type(time_type), intent(in) :: Time
	      real, INTENT(IN), optional:: time_total

! Local
	      real pref(npz)
              integer  i,j,k
	      real  ty, tz, akap 
	      real  p0, t0, sday, rkv, rka, rks, rkt, sigb, rsgb
	      real  tmp
	      real  ap0k, algpk
	      real  tey, tez, fac, pw, sigl
	      real  h0, dz
	      real  dt_tropic
	      real  rmr, rms
	      real  relx, tau
	      real  t_st, t_ms
	      real  f1
	      real  pc, c1, kf_day
              real  rad_ratio        ! small earth ratio

	      real, allocatable, SAVE ::  rf(:)
	      real :: frac(is-ng:ie+ng,js-ng:je+ng)
	      real ::  pl(is:ie, js:je, 1:npz)
	      real :: teq(is:ie, js:je, 1:npz)

      ty = 60.0   ! origina; setting
      tz = 10.0   ! Original value from H-S was 10.
      akap = 2./7.

      p0 = 100000.
      t0 = 200.
      h0 = 7.

      rad_ratio = radius / 6371.0e3

        sday = 24.*3600.
      kf_day = sday * rad_ratio

      rkv = (0.5*pdt) / kf_day
      rka = pdt      / (40.*kf_day)
      rks = pdt      / (4.0*kf_day)

! For strat-mesosphere
      t_ms = 10.*rad_ratio
      t_st = 40.*rad_ratio

      tau = (t_st - t_ms) / log(100.)
      rms = pdt/(t_ms*sday)
      rmr =  1./(1.+rms)

      sigb = 0.7
      rsgb = 1./(1.-sigb)
      ap0k = 1./p0**akap
      algpk = log(ap0k)

      do k=1,npz
         pref(k) = ak(k) + bk(k)*1.E5
      enddo

! Setup...
      if ( rayf .and. (.not. rf_initialized) ) then
           allocate( rf(npz) )
           c1 = 1. / (36.*3600)
           pc = 1.
           if(master) write(6,*) 'HSWF Forcing ...' 
           do k=1,ks
              tmp = (ak(k+1)-ak(k))/log(ak(k+1)/ak(k))
              rf(k) = c1*(1.+tanh(log10(pc/tmp)))
             if(master) write(6,*) k, 0.01*tmp, 1./(rf(k)*sday)
             rf(k) = 1./(1.+pdt*rf(k))
           enddo
           if(master) write(6,*) ' '
           rf_initialized = .true.
      endif

   do k=1,npz
      do j=js,je
         do i=is,ie
            pl(i,j,k) = delp(i,j,k) / ( peln(i,k+1,j)-peln(i,k,j))
         enddo
      enddo 
   enddo

! Temperature forcing...
   do k=npz,1,-1
      do j=js,je
      do i=is,ie
         tey = ap0k*( 315.0 - ty*SIN(agrid(i,j,2))*SIN(agrid(i,j,2)) )
         tez = tz*( ap0k/akap )*COS(agrid(i,j,2))*COS(agrid(i,j,2)) 

      if (strat .and. pl(i,j,k) < 10000.    &
          .and. pl(i,j,k) > 100.  )  then
         dz = h0 * log(pl(i,j,k+1)/pl(i,j,k))
! Lapse rate above tropic stratopause is 2.25 deg/km
! Relaxation time is t_st days at 100 mb (as H-S) and gradually
! decreases to t_ms Days at and above the stratopause
         relx =  t_ms + tau*log(0.01*pl(i,j,k))
         relx = pdt/(relx*sday)
         dt_tropic = 2.25*COS(agrid(i,j,2)) * dz
         teq(i,j,k)  =  teq(i,j,k+1) + dt_tropic
         pt(i,j,k) = (pt(i,j,k)+relx*teq(i,j,k))/(1.+relx)
      elseif (strat .and. pl(i,j,k) <= 100.)  then
!
! Mesosphere: defined as the region above 1 mb
!
         dz = h0 * log(pl(i,j,k+1)/pl(i,j,k))
         dt_tropic = -2.25*COS(agrid(i,j,2)) * dz
         teq(i,j,k) = teq(i,j,k+1) + dt_tropic
         pt(i,j,k) = (pt(i,j,k)+rms*teq(i,j,k))*rmr
      else

! Trop:  strictly Held-Suarez

         sigl = pl(i,j,k)/pe(i,npz+1,j)
         f1 = max(0., (sigl-sigb) * rsgb )
         teq(i,j,k) = tey - tez*(log(pkz(i,j,k))+algpk)
         teq(i,j,k) = max(t0, teq(i,j,k)*pkz(i,j,k))
         rkt = rka + (rks-rka)*f1*(COS(agrid(i,j,2))**4.0)
         pt(i,j,k) = (pt(i,j,k)+rkt*teq(i,j,k))/(1.+rkt)
      endif
      enddo     !i-loop
   enddo     !j-loop
   enddo     !k-loop

! Velocity dissipation damping
  do 2000 k=1,npz

     frac = 0. !Avoids severe errors at nested grid boundary
      if (rayf .and. k <= ks) then
! Apply Rayleigh friction
          do j=js,je+1
             do i=is,ie
                u(i,j,k) = u(i,j,k)*rf(k)
              enddo
          enddo
          do j=js,je
             do i=is,ie+1
                v(i,j,k) = v(i,j,k)*rf(k)
             enddo
          enddo
       else
! Surface Rayleigh friction according to Held-Suarez
          do j=js,je
             do i=is,ie
                sigl = 0.5*(pe(i,k,j)+pe(i,k+1,j)) / pe(i,npz+1,j)
                frac(i,j) = max(0., (sigl-sigb)*rsgb )
             enddo
          enddo

     call mpp_update_domains( frac, domain )

! Backward adjustment
  do j=js,je+1
    do i=is,ie+1
      fac = frac(i,j)+frac(i,j-1)
      if (fac > 0.) then
         u(i,j,k) = u(i,j,k)/(1.+rkv*fac)
      endif
    enddo
  enddo

  do j=js,je
    do i=is,ie+1
       fac = frac(i,j)+frac(i-1,j)
       if (fac > 0.) then
           v(i,j,k) = v(i,j,k)/(1.+rkv*fac)
       endif
    enddo
  enddo

  endif

2000  continue

#ifdef DO_AGE
       if( nq/=0 )     &
        call age_of_air(is, ie, js, je, npz, ng, time_total, pe, q(is-ng,js-ng,1,nq))
#endif


  end subroutine Held_Suarez_Strat


#ifdef MARS_GCM
! Forcing for MARS_GCM
 subroutine Mars_phys(npx, npy, npz, is, ie, js, je, ng, nq,  &
		     u_dt, v_dt, t_dt, q_dt, u, v, ua, va,   &
		     pt, q, pe, delp, peln, oro, hydrostatic, &
		     pdt, agrid, ak, bk, rayf, p_ref, master, Time, time_total)

	      integer, INTENT(IN) :: npx, npy, npz
	      integer, INTENT(IN) :: is, ie, js, je, ng, nq
	      real   , INTENT(IN) :: pdt
	      real   , INTENT(IN) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
	      real   , INTENT(IN) :: ak(npz+1), bk(npz+1)
	      logical, INTENT(IN) :: rayf, master
	      real, INTENT(IN):: oro(is:ie,js:je)       ! land fraction
	      logical, INTENT(IN):: hydrostatic
	      real, INTENT(IN):: p_ref

	      type(time_type), intent(in) :: Time
	      real, INTENT(IN), optional:: time_total

	      real   , INTENT(INOUT) :: u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
	      real   , INTENT(INOUT) :: v(is-ng:ie+1+ng,js-ng:je+  ng,npz)

	      real, INTENT(INOUT)::   pt(is-ng:ie+ng,js-ng:je+ng,npz)
	      real, INTENT(INOUT):: delp(is-ng:ie+ng,js-ng:je+ng,npz)
	      real, INTENT(INOUT)::    q(is-ng:ie+ng,js-ng:je+ng,npz, nq)
	      real, INTENT(INOUT)::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
	      real, INTENT(INOUT):: peln(is  :ie   ,1:npz+1,js  :je  )

	! Tendencies:
	      real, INTENT(INOUT):: u_dt(is-ng:ie+ng,js-ng:je+ng,npz)
	      real, INTENT(INOUT):: v_dt(is-ng:ie+ng,js-ng:je+ng,npz)
	      real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
	      real, INTENT(INOUT):: q_dt(is:ie,js:je,npz,nq)

	      real, INTENT(INOUT):: ua(is-ng:ie+ng,js-ng:je+ng,npz)
	      real, INTENT(INOUT):: va(is-ng:ie+ng,js-ng:je+ng,npz)

	! Local
	      real pedge(npz+1)
	      real pref(npz)

	      real  sday, rkv, rkt, sigb, rsgb, sigl, cs
	      real frac
	      real rdt                    ! Reciprocal of dt
	      character (len=128) :: filename
	      integer  i,j,k

	      sday = 24.*3600.

	      do k=1,npz+1
		 pedge(k) = ak(k) + bk(k)*p_ref
	      enddo

	      do k=1,npz
		 pref(k) = (pedge(k+1)-pedge(k)) / log(pedge(k+1)/pedge(k))
	      enddo

	      if ( .not. tmars_initialized ) then
		  allocate( tmars(is:ie,js:je,npz) )

		  filename= 'INPUT/teq.nc' 
		  if( file_exist( trim( filename ) ) ) then 
		      call read_teq( filename, npz, agrid(is:ie,js:je,1), agrid(is:ie,js:je,2),  &
				     pref, tmars(is:ie,js:je,:)  )
		      if(master) write(6,*) 'TEQ for Mars initialized.'
		  else
#ifdef FAIL_SAFE
	      call mpp_error(FATAL,'Mars_GCM: TEQ data not found')
#else
		      do k=1,npz
			 do j=js,je
			    do i=is,ie
			       tmars(i,j,k) = 100.+ 100.*max(0.25, (1.-sin(agrid(i,j,2)))*pref(k)/pedge(npz+1))
	!                      tmars(i,j,k) = 120. + 100.*    &
	!                      max(0.2, (1.-sin(agrid(i,j,2)))*log(pref(k))/log(pedge(npz+1)))
			    enddo 
			 enddo 
		      enddo
		      if(master) write(6,*) 'Data for Mars not found; using analytic profile'
#endif
		  endif 
		  cs = sqrt( rdgas * 273./(1.-kappa) ) 
		  if(master) write(6,*) 'Sound speed (T=273)=', cs
		  tmars_initialized = .true.
	      endif

	! ***  Newtonian cooling/relaxation
	      rdt = 1. / pdt
	!     rkt = pdt / (4.*sday)
	      rkt = pdt / (8.*sday)

	      do k=1,npz
		 do j=js,je
		    do i=is,ie
		       t_dt(i,j,k) = rkt*(tmars(i,j,k)-pt(i,j,k))/(1.+rkt) * rdt
		    enddo 
		 enddo 
	      enddo

	! *** Surface Rayleigh friction according to Held-Suarez
	      rkv = pdt / (1.*sday)          ! 1 day
	!     rkv = pdt / (2.*sday)          ! 2 day
	      sigb = 0.7
	      rsgb = 1./(1.-sigb)

	      do k=1,npz
		  do j=js,je
		     do i=is,ie
			sigl = 0.5*(pe(i,k,j)+pe(i,k+1,j)) / pe(i,npz+1,j)
			frac = rkv * (sigl-sigb)*rsgb
			if (frac > 0.) then
			    u_dt(i,j,k) = -ua(i,j,k)*frac/(1.+frac) * rdt
			    v_dt(i,j,k) = -va(i,j,k)*frac/(1.+frac) * rdt
			endif
		     enddo
		  enddo
	      enddo



	 end subroutine Mars_phys


	 subroutine read_teq ( filename, nlevels, lon, lat, pstd, tout )

	!-----------------------------------------------------------------------
	!
	! routine for initializing the radiative-convective temperature cross-section
	!
	!-----------------------------------------------------------------------



	   character(len=128), intent(in)            :: filename
	   integer, intent(in):: nlevels
	   real,    intent(in),  dimension(:,:)      :: lon, lat 
	   real,    intent(in)       ::  pstd(nlevels)

	   real,    intent(out),  dimension(:,:,:)   ::  tout

	!-----------------------------------------------------------------------
	   integer  unit, io, ierr, i, j, k, klev, kk

	   integer  im, jm, km, fld_dims(4)
	   real    ::  frac

	   real, dimension(:,:,:),  allocatable  ::   teq_inpt
	   real, dimension(:,:),    allocatable  ::   txy

	   real, dimension(:),  allocatable  ::   lat_inpt, pres_inpt,  presh_inpt
	   real, dimension(:),  allocatable  ::   lonb_inpt, latb_inpt 



	!       Get input field teq
	   call field_size( trim(filename), 'lat', fld_dims )

	   allocate( lat_inpt (fld_dims(1)  ) )
	   allocate( latb_inpt(fld_dims(1)+1) )

	   call read_data( trim(filename), 'lat',  lat_inpt,  no_domain=.true. )
	   call read_data( trim(filename), 'latb', latb_inpt, no_domain=.true. )

	   call field_size( trim(filename), 'lonb', fld_dims )

	   allocate( lonb_inpt (fld_dims(1)  ) )

	   call read_data( trim(filename), 'lonb', lonb_inpt, no_domain=.true. )


	   call field_size( trim(filename), 'pfull', fld_dims )
	   allocate( pres_inpt( fld_dims(1) ) )
	   call read_data( trim(filename), 'pfull', pres_inpt, no_domain=.true. )
	   print *, 'pfull dims:  ', fld_dims 
	   pres_inpt= 100.0 * pres_inpt



	   call field_size( trim(filename), 'teq', fld_dims )

	   im= fld_dims(1);  jm= fld_dims(2);  km= fld_dims(3) 
	       print *, 'Input Teq dims:  ', fld_dims 

	   allocate( teq_inpt( im,jm,km ) )
	   allocate( txy     ( im,jm )    )

	   call read_data( trim(filename), 'teq', teq_inpt, no_domain=.true. )

	   latb_inpt(:)= latb_inpt(:)/RADIAN
	   lonb_inpt(:)= lonb_inpt(:)/RADIAN


	!           If km != nlevels  then require vertical interpolation 
	  if( nlevels > km .or. nlevels < km )  then

	     DO k= 1, nlevels
		if( pres_inpt(1) > pstd(k) ) then
		   klev= 2;  frac= 1.0
		else
		  DO kk= 2, km
		    if( pres_inpt(kk) > pstd(k) ) then
		       frac=  ( pres_inpt(kk) - pstd(k) )/(pres_inpt(kk)-pres_inpt(kk-1))
		       klev= kk
		       exit
		    endif
		  ENDDO

		  if( kk > km )  then
		       klev= km;  frac= 0.0
		  endif 
		endif

	!             Complete pressure interpolation 
		DO i= 1, im
		  DO j= 1, jm
		    txy(i,j)= (1.0-frac)*teq_inpt(i,j,klev) + frac *teq_inpt(i,j,klev-1)
		  ENDDO
		ENDDO

	!             Carry out horizontal interpolation 
		 call horiz_interp( txy(:,:), lonb_inpt, latb_inpt, lon, lat,    &
					      tout(:,:,k), interp_method= 'bilinear' )
	     ENDDO      ! -----------  end loop over k 

	  else

	       DO k= 1, nlevels
		 call horiz_interp( teq_inpt(:,:,k), lonb_inpt, latb_inpt, lon, lat, &
				     tout(:,:,k), interp_method= 'bilinear' )
	       ENDDO
	  endif

	   deallocate ( teq_inpt )
	   deallocate ( txy )
	   deallocate ( lat_inpt )
	   deallocate ( latb_inpt )
	   deallocate ( lonb_inpt )
	   deallocate ( pres_inpt )


	 end subroutine read_teq

#endif MARS_GCM

      subroutine age_of_air(is, ie, js, je, km, ng, time, pe, q)

      integer is, ie, js, je
      integer km
      integer ng

! q is the age tracer
! Need to be converted to mixing ratio (mass of tracer / dry_air-mass)
! Ignore this inconsistency for now.

      real, intent(inout):: pe(is-1:ie+1, km+1, js-1:je+1)
      real, intent(in):: time        ! accumulated time since init
      real, intent(inout):: q(is-ng:ie+ng,js-ng:je+ng,km)

! Local
      integer i, j, k
      real p_source      ! source level (pa)
      real ascale
      real tiny
      parameter ( tiny = 1.e-6 )
      parameter ( p_source = 75000. )
      parameter ( ascale = 5.e-6 / 60. )

!$omp parallel do private(i, j, k)
      do k=1,km
        do j=js,je
            do i=is,ie
               if( time < tiny ) then
                   q(i,j,k) = 0.
               elseif( pe(i,k,j) >= p_source ) then
                   q(i,j,k) = ascale * time
               endif
            enddo
        enddo          ! j-loop
      enddo             ! k-loop

      end subroutine age_of_air

end module hswf_mod
