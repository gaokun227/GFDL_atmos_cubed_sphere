subroutine reed_simple_physics (pcols, pver, dtime, t, q, u, v, pmid, pint, pdel, rpdel, ps, precl, TSurf, rh, cond_only, dqdt, dTdt, dudt, dvdt)
!----------------------------------------------------------------------- 
! 
! Purpose: Simple Physics Package
!
! Author: K. A. Reed (University of Michigan, kareed@umich.edu)
!         version 4 
!         May/31/2012
!         uploaded on June/8/2012
!
!         cosmetic changes in comparison to version 3:
!         parameter list now includes pcols and pver
!         ncol has been replaced by pcols in the code
! 
! Description: Includes large-scale precipitation, surface fluxes and
!              boundary-leyer mixing. The processes are time-split
!              in that order. A partially implicit formulation is
!              used to foster numerical stability.
!              The routine assumes that the model levels are ordered
!              in a top-down approach, e.g. level 1 denotes the uppermost
!              full model level.
!
!              This routine is based on an implementation which was
!              developed for the NCAR Community Atmosphere Model (CAM).
!              Adjustments for other models will be necessary.
!
!              The routine provides both updates of the state variables
!              u, v, T, q (these are local copies of u,v,T,q within this physics
!              routine) and also collects their time tendencies.
!              The latter might be used to couple the physics and dynamics
!              in a process-split way. For a time-split coupling, the final
!              state should be given to the dynamical core for the next time step.
!
!
! Reference: Reed, K. A. and C. Jablonowski (2012), Idealized tropical cyclone 
!            simulations of intermediate complexity: A test case for AGCMs, 
!            J. Adv. Model. Earth Syst., Vol. 4, M04001, doi:10.1029/2011MS000099
! 
! NOTE: 17 jul 12 lmh
!       This has been modified to return only the tendencies and not
!       to update the model; the model variables are updated in
!       fv_update_phys.F90. However we do need to use the
!       tendency computed by the surface fluxes to advance the
!       boundary layer mixing to be in line with the original
!       formulation of these physics routines.
!--------------------------------------------------------------------

  ! use physics_types     , only: physics_dme_adjust   ! This is for CESM/CAM
  ! use cam_diagnostics,    only: diag_phys_writeout   ! This is for CESM/CAM

  use fv_sg_mod, only: qsmith
  implicit none

!---- version number -----
   character(len=128) :: version = '$Id: simple_physics.F90,v 20.0 2013/12/13 23:04:13 fms Exp $'
   character(len=128) :: tagname = '$Name: tikal_201409 $'

   integer, parameter :: r8 = selected_real_kind(12)

!
! Input arguments - MODEL DEPENDENT
!
   integer, intent(in)  :: pcols        ! Set number of atmospheric columns       
   integer, intent(in)  :: pver         ! Set number of model levels
   real(r8), intent(in) :: dtime        ! Set model physics timestep
   logical, intent(IN) :: cond_only     ! Only do large-scale condensation
!
! Input/Output arguments 
!
!  pcols is the maximum number of vertical columns per 'chunk' of atmosphere
!
   real(r8), intent(inout) :: t(pcols,pver)      ! Temperature at full-model level (K)
   real(r8), intent(inout) :: q(pcols,pver)      ! Specific Humidity at full-model level (kg/kg)
   real(r8), intent(inout) :: u(pcols,pver)      ! Zonal wind at full-model level (m/s)
   real(r8), intent(inout) :: v(pcols,pver)      ! Meridional wind at full-model level (m/s)
   real(r8), intent(inout) :: pmid(pcols,pver)   ! Pressure is full-model level (Pa)
   real(r8), intent(inout) :: pint(pcols,pver+1) ! Pressure at model interfaces (Pa)
   real(r8), intent(inout) :: pdel(pcols,pver)   ! Layer thickness (Pa)
   real(r8), intent(inout) :: rpdel(pcols,pver)  ! Reciprocal of layer thickness (1/Pa)
   real(r8), intent(inout) :: ps(pcols)          ! Surface Pressue (Pa)
   real(r8), intent(IN)    :: Tsurf(pcols)              ! Specified SST
!
! Output arguments 
!
   real(r8), intent(out) :: precl(pcols)         ! Precipitation rate (m_water / s)
   real(r8), intent(OUT)   :: RH(pcols,pver)     !Output RH

!
!---------------------------Local workspace-----------------------------
!

! Integers for loops

   integer  i,k                         ! Longitude, level indices

! Physical Constants - Many of these may be model dependent

   real(r8) gravit                      ! Gravity
   real(r8) rair                        ! Gas constant for dry air
   real(r8) cpair                       ! Specific heat of dry air 
   real(r8) latvap                      ! Latent heat of vaporization
   real(r8) rh2o                        ! Gas constant for water vapor
   real(r8) epsilo                      ! Ratio of gas constant for dry air to that for vapor
   real(r8) zvir                        ! Constant for virtual temp. calc. =(rh2o/rair) - 1

! Simple Physics Specific Constants 

   real(r8) T0                          ! Control temp for calculation of qsat
   real(r8) e0                          ! Saturation vapor pressure at T0 for calculation of qsat
   real(r8) rhow                        ! Density of Liquid Water
   real(r8) p0                          ! Constant for calculation of potential temperature
   real(r8) Cd0                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cd1                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cm                          ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) v20                         ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
   real(r8) C                           ! Drag coefficient for sensible heat and evaporation

! Physics Tendency Arrays
  real(r8), intent(OUT) ::  dtdt(pcols,pver)             ! Temperature tendency 
  real(r8), intent(OUT) ::  dqdt(pcols,pver)             ! Specific humidity tendency
  real(r8), intent(OUT) ::  dudt(pcols,pver)             ! Zonal wind tendency
  real(r8), intent(OUT) ::  dvdt(pcols,pver)             ! Meridional wind tendency

  real(r8) q_orig(pcols,pver)

! Temporary variables for tendency calculations

   real(r8) tmp                         ! Temporary
   real(r8) qsat(pcols,pver)                        ! Saturation vapor pressure
   !real(r8) qsat                        ! Saturation vapor pressure
   real(r8) qsats(pcols,1)                       ! Saturation vapor pressure of SST
   !real(r8) qsats                       ! Saturation vapor pressure of SST

   !Updated variables, using the tendencies from LS precipitation and
   !surface fluxes after each is computed to yield the same answer as
   !if the tendencies were added into the original variables
   !immediately after each step of the physics was called.
!!$   real(r8) :: tup(pcols,pver)      ! Temperature at full-model level (K)
!!$   real(r8) :: qup(pcols,pver)      ! Specific Humidity at full-model level (kg/kg)
!!$   real(r8) :: uup(pcols,pver)      ! Zonal wind at full-model level (m/s)
!!$   real(r8) :: vup(pcols,pver)      ! Meridional wind at full-model level (m/s)

! Variables for Boundary Layer Calculation

   real(r8) wind(pcols)                 ! Magnitude of Wind
   real(r8) Cd(pcols)                   ! Drag coefficient for momentum
   real(r8) Km(pcols,pver+1)            ! Eddy diffusivity for boundary layer calculations 
   real(r8) Ke(pcols,pver+1)            ! Eddy diffusivity for boundary layer calculations
   real(r8) rho                         ! Density at lower/upper interface
   real(r8) za(pcols)                   ! Heights at midpoints of first model level
   real(r8) dlnpint                     ! Used for calculation of heights
   real(r8) pbltop                      ! Top of boundary layer
   real(r8) pblconst                    ! Constant for the calculation of the decay of diffusivity 
   real(r8) CA(pcols,pver)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CC(pcols,pver)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CE(pcols,pver+1)            ! Matrix Coefficents for PBL Scheme
   real(r8) CAm(pcols,pver)             ! Matrix Coefficents for PBL Scheme
   real(r8) CCm(pcols,pver)             ! Matrix Coefficents for PBL Scheme
   real(r8) CEm(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFu(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFv(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFt(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFq(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme


! Variable for Dry Mass Adjustment, this dry air adjustment is necessary to
! conserve the mass of the dry air

   real(r8) qini(pcols,pver)            ! Initial specific humidity
   real(r8) ddelp(pcols)                ! Column delp change

!===============================================================================
!
! Physical Constants - MAY BE MODEL DEPENDENT 
!
!===============================================================================
   gravit = 9.80616_r8                   ! Gravity
   rair   = 287.0_r8                     ! Gas constant for dry air
   cpair  = 1.004e3_r8                   ! Specific heat of dry air 
   latvap = 2.5e6_r8                     ! Latent heat of vaporization
   rh2o   = 461.5_r8                     ! Gas constant for water vapor
   epsilo = rair/rh2o                    ! Ratio of gas constant for dry air to that for vapor
   zvir   = (rh2o/rair) - 1._r8          ! Constant for virtual temp. calc. =(rh2o/rair) - 1 is approx. 0.608

!===============================================================================
!
! Local Constants for Simple Physics
!
!===============================================================================
!!$      C        = 0.0011_r8*2.      ! From Smith and Vogl 2008
      C        = 0.0011_r8      ! From Smith and Vogl 2008
      !Tsurf    = 302.15_r8      ! Constant Value for SST
      T0       = 273.16_r8      ! control temp for calculation of qsat
      e0       = 610.78_r8      ! saturation vapor pressure at T0 for calculation of qsat
      rhow     = 1000.0_r8      ! Density of Liquid Water 
!!$      Cd0      = 0.007_r8      ! Constant for Cd calc. Smith and Vogl 2008
!!$      Cd1      = 0.00065_r8    ! Constant for Cd calc. Smith and Vogl 2008
!!$      Cm       = 0.02_r8       ! Constant for Cd calc. Smith and Vogl 2008
      Cd0      = 0.0007_r8      ! Constant for Cd calc. Smith and Vogl 2008
      Cd1      = 0.000065_r8    ! Constant for Cd calc. Smith and Vogl 2008
      Cm       = 0.002_r8       ! Constant for Cd calc. Smith and Vogl 2008
      v20      = 20.0_r8        ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
      p0       = 100000.0_r8    ! Constant for potential temp calculation
      pbltop   = 85000._r8      ! Top of boundary layer
      pblconst = 10000._r8      ! Constant for the calculation of the decay of diffusivity

!===============================================================================
!
! Definition of local arrays
!
!===============================================================================

      q_orig = q

!
! Calculate hydrostatic height za of the lowest model level
!
     do i=1,pcols 
        dlnpint = log(ps(i)) - log(pint(i,pver))                 ! ps(i) is identical to pint(i,pver+1), note: this is the correct sign (corrects typo in paper)
        za(i) = rair/gravit*t(i,pver)*(1._r8+zvir*q(i,pver))*0.5_r8*dlnpint
     end do
!
! Set Initial Specific Humidity - For dry mass adjustment at the end
!
     qini(:pcols,:pver) = q(:pcols,:pver)

!===============================================================================
!
! Set initial physics time tendencies and precipitation field to zero
!
!===============================================================================
     dtdt(:pcols,:pver)  = 0._r8            ! initialize temperature tendency with zero
     dqdt(:pcols,:pver)  = 0._r8            ! initialize specific humidity tendency with zero
     dudt(:pcols,:pver)  = 0._r8            ! initialize zonal wind tendency with zero
     dvdt(:pcols,:pver)  = 0._r8            ! initialize meridional wind tendency with zero
     precl(:pcols) = 0._r8                  ! initialize precipitation rate with zero

!===============================================================================
!
! Large-Scale Condensation and Precipitation Rate
!
!===============================================================================
!
! Calculate Tendencies
!
     call qsmith(pcols,pver,1, t, pmid, q, qsat)
      do k=1,pver
         do i=1,pcols
            !qsat = epsilo*e0/pmid(i,k)*exp(-latvap/rh2o*((1._r8/t(i,k))-1._r8/T0))  ! saturation specific humidity
            !if (q(i,k) > qsat) then                                                 ! saturated?
            if (q(i,k) > qsat(i,k)) then                                                 ! saturated?
               !tmp  = 1._r8/dtime*(q(i,k)-qsat)/(1._r8+(latvap/cpair)*(epsilo*latvap*qsat/(rair*t(i,k)**2)))
               tmp  = 1._r8/dtime*(q(i,k)-qsat(i,k))/(1._r8+(latvap/cpair)*(epsilo*latvap*qsat(i,k)/(rair*t(i,k)**2)))
               dtdt(i,k) = dtdt(i,k)+latvap/cpair*tmp
               dqdt(i,k) = dqdt(i,k)-tmp
               precl(i) = precl(i) + tmp*pdel(i,k)/(gravit*rhow)                    ! precipitation rate, computed via a vertical integral
                                                                                    ! corrected in version 1.3
            end if
            !rh(i,k) = q(i,k)/qsat*100.
            rh(i,k) = q(i,k)/qsat(i,k)*100.
         end do
      end do
!
! Update moisture and temperature fields from Large-Scale Precipitation Scheme
!
!!$      do k=1,pver
!!$         do i=1,pcols
!!$            tup(i,k) =  t(i,k) + dtdt(i,k)*dtime    ! update the state variables T and q
!!$            qup(i,k) =  q(i,k) + dqdt(i,k)*dtime
!!$         end do
!!$      end do

!===============================================================================
! Send variables to history file - THIS PROCESS WILL BE MODEL SPECIFIC
!
! note: The variables, as done in many GCMs, are written to the history file
!       after the moist physics process.  This ensures that the moisture fields
!       are somewhat in equilibrium.
!===============================================================================
  !  call diag_phys_writeout(state)   ! This is for CESM/CAM

!===============================================================================
!
! Turbulent mixing coefficients for the PBL mixing of horizontal momentum,
! sensible heat and latent heat
!
! We are using Simplified Ekman theory to compute the diffusion coefficients 
! Kx for the boundary-layer mixing. The Kx values are calculated at each time step
! and in each column.
!
!===============================================================================

if (.not. cond_only) then
!
! Compute magnitude of the wind and drag coeffcients for turbulence scheme:
! they depend on the conditions at the lowest model level and stay constant
! up to the 850 hPa level. Above this level the coefficients are decreased
! and tapered to zero. At the 700 hPa level the strength of the K coefficients
! is about 10% of the maximum strength. 
!
     do i=1,pcols
        wind(i) = sqrt(u(i,pver)**2+v(i,pver)**2)    ! wind magnitude at the lowest level
     end do
     do i=1,pcols
        Ke(i,pver+1) = C*wind(i)*za(i)
        if( wind(i) .lt. v20) then
           Cd(i) = Cd0+Cd1*wind(i) 
           Km(i,pver+1) = Cd(i)*wind(i)*za(i)
        else
           Cd(i) = Cm
           Km(i,pver+1) = Cm*wind(i)*za(i)
        endif
     end do

      do k=1,pver
         do i=1,pcols
            if( pint(i,k) .ge. pbltop) then
               Km(i,k) = Km(i,pver+1)                 ! constant Km below 850 hPa level
               Ke(i,k) = Ke(i,pver+1)                 ! constant Ke below 850 hPa level
            else
               Km(i,k) = Km(i,pver+1)*exp(-(pbltop-pint(i,k))**2/(pblconst)**2)  ! Km tapered to 0
               Ke(i,k) = Ke(i,pver+1)*exp(-(pbltop-pint(i,k))**2/(pblconst)**2)  ! Ke tapered to 0
            end if 
         end do
      end do     


!===============================================================================
! Update the state variables u, v, t, q with the surface fluxes at the
! lowest model level, this is done with an implicit approach
! see Reed and Jablonowski (JAMES, 2012)
!
! Sea Surface Temperature Tsurf is constant for tropical cyclone test 5-1
! Tsurf needs to be dependent on latitude for the
! moist baroclinic wave test 4-3 with simple-physics, adjust
!===============================================================================

      !Here qsats is only one in the y-direction
      call qsmith(pcols,1,1, reshape(tsurf, (/pcols,1/)), reshape(ps, (/pcols,1/)), q(:,pver:pver), qsats)


      !uup = u
      !vup = v

     do i=1,pcols 
        !qsats = epsilo*e0/ps(i)*exp(-latvap/rh2o*((1._r8/Tsurf(i))-1._r8/T0))  ! saturation specific humidity at the surface
        dudt(i,pver) = dudt(i,pver) + (u(i,pver) &
                            /(1._r8+Cd(i)*wind(i)*dtime/za(i))-u(i,pver))/dtime
        dvdt(i,pver) = dvdt(i,pver) + (v(i,pver) &
                            /(1._r8+Cd(i)*wind(i)*dtime/za(i))-v(i,pver))/dtime
        !uup(i,pver)   = u(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))
        !vup(i,pver)   = v(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))
        dtdt(i,pver) = dtdt(i,pver) +((t(i,pver)+C*wind(i)*Tsurf(i)*dtime/za(i)) &
                            /(1._r8+C*wind(i)*dtime/za(i))-t(i,pver))/dtime 
!!$        dtdt(i,pver) = dtdt(i,pver) +((tup(i,pver)+C*wind(i)*Tsurf(i)*dtime/za(i)) &
!!$                            /(1._r8+C*wind(i)*dtime/za(i))-tup(i,pver))/dtime 
!!$        tup(i,pver)   = (tup(i,pver)+C*wind(i)*Tsurf(i)*dtime/za(i)) &
!!$                            /(1._r8+C*wind(i)*dtime/za(i))  
        dqdt(i,pver) = dqdt(i,pver) +((q(i,pver)+C*wind(i)*qsats(i,1)*dtime/za(i)) &
                            /(1._r8+C*wind(i)*dtime/za(i))-q(i,pver))/dtime
!!$        dqdt(i,pver) = dqdt(i,pver) +((qup(i,pver)+C*wind(i)*qsats(i,1)*dtime/za(i)) &
!!$                            /(1._r8+C*wind(i)*dtime/za(i))-qup(i,pver))/dtime
!!$        qup(i,pver) = (qup(i,pver)+C*wind(i)*qsats(i,1)*dtime/za(i))/(1._r8+C*wind(i)*dtime/za(i))
!!! DEBUG CODE
        !if (wind(i) > 10.) print*, 'dudt 1:', dudt(i,pver), dvdt(i,pver)
        !rh(i,pver) = q(i,pver)/qsats*100.
        rh(i,pver) = q(i,pver)/qsats(i,1)*100.
!!! END DEBUG CODE
     end do
!===============================================================================


!===============================================================================
! Boundary layer mixing, see Reed and Jablonowski (JAMES, 2012)
!===============================================================================
! Calculate Diagonal Variables for Implicit PBL Scheme
!
      do k=1,pver-1
         do i=1,pcols
            rho = (pint(i,k+1)/(rair*(t(i,k+1)+t(i,k))/2.0_r8))
!!$            rho = (pint(i,k+1)/(rair*(tup(i,k+1)+tup(i,k))/2.0_r8))
            CAm(i,k)   = rpdel(i,k)*dtime*gravit*gravit*Km(i,k+1)*rho*rho   &
                         /(pmid(i,k+1)-pmid(i,k))    
            CCm(i,k+1) = rpdel(i,k+1)*dtime*gravit*gravit*Km(i,k+1)*rho*rho &
                         /(pmid(i,k+1)-pmid(i,k))
            CA(i,k)    = rpdel(i,k)*dtime*gravit*gravit*Ke(i,k+1)*rho*rho   &
                         /(pmid(i,k+1)-pmid(i,k))
            CC(i,k+1)  = rpdel(i,k+1)*dtime*gravit*gravit*Ke(i,k+1)*rho*rho &
                         /(pmid(i,k+1)-pmid(i,k))
         end do
      end do
      do i=1,pcols
         CAm(i,pver) = 0._r8
         CCm(i,1) = 0._r8
         CEm(i,pver+1) = 0._r8
         CA(i,pver) = 0._r8
         CC(i,1) = 0._r8
         CE(i,pver+1) = 0._r8
         CFu(i,pver+1) = 0._r8
         CFv(i,pver+1) = 0._r8
         CFt(i,pver+1) = 0._r8
         CFq(i,pver+1) = 0._r8 
      end do
      do i=1,pcols
         do k=pver,1,-1
            CE(i,k)  = CC(i,k)/(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1)) 
            CEm(i,k) = CCm(i,k)/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFu(i,k) = (u(i,k)+CAm(i,k)*CFu(i,k+1)) &
                       /(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFv(i,k) = (v(i,k)+CAm(i,k)*CFv(i,k+1)) &
                       /(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFt(i,k) = ((p0/pmid(i,k))**(rair/cpair)*t(i,k)+CA(i,k)*CFt(i,k+1)) &
                       /(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1)) 
            CFq(i,k) = (q(i,k)+CA(i,k)*CFq(i,k+1)) &
                       /(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1))
!!$            CFu(i,k) = (uup(i,k)+CAm(i,k)*CFu(i,k+1)) &
!!$                       /(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
!!$            CFv(i,k) = (vup(i,k)+CAm(i,k)*CFv(i,k+1)) &
!!$                       /(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
!!$            CFt(i,k) = ((p0/pmid(i,k))**(rair/cpair)*tup(i,k)+CA(i,k)*CFt(i,k+1)) &
!!$                       /(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1)) 
!!$            CFq(i,k) = (qup(i,k)+CA(i,k)*CFq(i,k+1)) &
!!$                       /(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1))
        end do
      end do

!
! Calculate the updated temperature, specific humidity and horizontal wind
!
! First we need to calculate the updates at the top model level
!
      do i=1,pcols
            dudt(i,1)  = dudt(i,1)+(CFu(i,1)-u(i,1))/dtime
            dvdt(i,1)  = dvdt(i,1)+(CFv(i,1)-v(i,1))/dtime
!!$            dudt(i,1)  = dudt(i,1)+(CFu(i,1)-uup(i,1))/dtime
!!$            dvdt(i,1)  = dvdt(i,1)+(CFv(i,1)-vup(i,1))/dtime
            !u(i,1)    = CFu(i,1)
            !v(i,1)    = CFv(i,1)
            dtdt(i,1)  = dtdt(i,1)+(CFt(i,1)*(pmid(i,1)/p0)**(rair/cpair)-t(i,1))/dtime  ! corrected in version 1.3
!!$            dtdt(i,1)  = dtdt(i,1)+(CFt(i,1)*(pmid(i,1)/p0)**(rair/cpair)-tup(i,1))/dtime  ! corrected in version 1.3
            !t(i,1)    = CFt(i,1)*(pmid(i,1)/p0)**(rair/cpair)
            dqdt(i,1)  = dqdt(i,1)+(CFq(i,1)-q(i,1))/dtime
!!$            dqdt(i,1)  = dqdt(i,1)+(CFq(i,1)-qup(i,1))/dtime
            !q(i,1)  = CFq(i,1)
      end do
!
! Loop over the remaining level
!
      do i=1,pcols
         do k=2,pver
            dudt(i,k)  = dudt(i,k)+(CEm(i,k)*u(i,k-1)+CFu(i,k)-u(i,k))/dtime
            dvdt(i,k)  = dvdt(i,k)+(CEm(i,k)*v(i,k-1)+CFv(i,k)-v(i,k))/dtime
!!$            dudt(i,k)  = dudt(i,k)+(CEm(i,k)*uup(i,k-1)+CFu(i,k)-uup(i,k))/dtime
!!$            dvdt(i,k)  = dvdt(i,k)+(CEm(i,k)*vup(i,k-1)+CFv(i,k)-vup(i,k))/dtime
            !u(i,k)    = CEm(i,k)*u(i,k-1)+CFu(i,k) 
            !v(i,k)    = CEm(i,k)*v(i,k-1)+CFv(i,k)
            dtdt(i,k)  = dtdt(i,k)+((CE(i,k)*t(i,k-1) &
                              *(p0/pmid(i,k-1))**(rair/cpair)+CFt(i,k)) &
                              *(pmid(i,k)/p0)**(rair/cpair)-t(i,k))/dtime 
!!$            dtdt(i,k)  = dtdt(i,k)+((CE(i,k)*tup(i,k-1) &
!!$                              *(p0/pmid(i,k-1))**(rair/cpair)+CFt(i,k)) &
!!$                              *(pmid(i,k)/p0)**(rair/cpair)-tup(i,k))/dtime 
            !t(i,k)    = (CE(i,k)*t(i,k-1)*(p0/pmid(i,k-1))**(rair/cpair)+CFt(i,k)) &
            !                  *(pmid(i,k)/p0)**(rair/cpair)
            dqdt(i,k)  = dqdt(i,k)+(CE(i,k)*q(i,k-1)+CFq(i,k)-q(i,k))/dtime
!!$            dqdt(i,k)  = dqdt(i,k)+(CE(i,k)*qup(i,k-1)+CFq(i,k)-qup(i,k))/dtime
            !q(i,k)  = CE(i,k)*q(i,k-1)+CFq(i,k)
         end do
      end do

   endif !.not. cond_only
!===============================================================================
!
! Dry Mass Adjustment - THIS PROCESS WILL BE MODEL SPECIFIC 
!
! note: Care needs to be taken to ensure that the model conserves the total
!       dry air mass. 
!===============================================================================
  !  call physics_dme_adjust(state, tend, qini, dtime)   ! This is for CESM/CAM

   !Adjustment done in fv_update_physics

!!$   ddelp = 0.
!!$      do k=1,pver
!!$      do i=1,pcols
!!$         ddelp(i) = ddelp(i) + (q(i,k) - q_orig(i,k))*pdel(i,k)
!!$         pdel(i,k) = pdel(i,k)*(1. + q(i,k) - q_orig(i,k))
!!$         q(i,k) = q(i,k)/(1. + q(i,k) - q_orig(i,k))
!!$      enddo
!!$      enddo
!!$
!!$      do i=1,pcols
!!$         ps(i) = ps(i) + ddelp(i)
!!$      enddo

   return
 end subroutine reed_simple_physics

