!!
!! Dean M. Townsley 2009
!!
!! This subroutine interpolates within a 1-dimensional WD profile
!! by averaging over a cell width.  The WD profile, stored in the
!! Simulation Unit's module data area, is averaged over an interval of
!! size "dr" at radius "radius".  Returned are the density, temperature,
!! and c12 and ne22 mass fractions.  Temperature averaging is just
!! multiplicatively mass-weighted for lack of a better thing to do.
!!
!! Averaging over an interval is used to allow a graceful mapping onto
!! the unrefined grid during initial refinement.
!!
!! JAH: Added extra dimension to tabulated profile to simplify logical for including fluff
!! JAH: No longer assumes evenly spaced radial grid
!! JAH: Interpolate instead of average if grid resolution finer than profile
!!

subroutine sim_interpolate1dWd( volume, r_inner, r_outer, dens, temp, xhe4, xc12, xo16 )

  use Simulation_data, ONLY : sim_wd_dens_tab, sim_wd_temp_tab, sim_wd_he4_tab, &
                              sim_wd_c12_tab, sim_wd_o16_tab, sim_wd_npnts, &
                              sim_wd_radius, sim_wd_dr, sim_wd_rad_tab, &
                              sim_wd_vol_tab, sim_wd_radius, sim_wd_volume, &
                              sim_smallrho, sim_smallt, sim_globalMe
  use sim_local_interface, ONLY : interp1d_linear, locate
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"

  real, intent(in) :: volume, r_inner, r_outer
  real, intent(out):: dens, temp, xhe4, xc12, xo16

  real :: r1, r2, dr, dx, vol, mass, dvol, dmass, masstemp, masshe4, massc12, massne22
  integer :: imin, imax, i

  ! check if within the profile
  if ( r_inner > sim_wd_radius ) then
        ! if out of the profile, assign the last value in it, fluff will be handled by initBlock
        dens = sim_wd_dens_tab(sim_wd_npnts)
        temp = sim_wd_temp_tab(sim_wd_npnts)
        xhe4 = sim_wd_he4_tab(sim_wd_npnts)
        xc12 = sim_wd_c12_tab(sim_wd_npnts)
        xo16 = sim_wd_o16_tab(sim_wd_npnts)
  else
        ! Determine the indices in the WD profile corresponding to r_inner and r_outer
        imin = locate( r_inner, sim_wd_npnts, sim_wd_rad_tab(1:sim_wd_npnts) )
        imax = locate( r_outer, sim_wd_npnts, sim_wd_rad_tab(1:sim_wd_npnts) )
        !dx = (volume - sim_wd_vol_tab(imin))/(sim_wd_vol_tab(imax)-sim_wd_vol_tab(imin))

        ! If grid resolution is coarser than WD profile, construct an average
        ! If grid resolution is finer than WD profile, do simple interpolation
        if ( imax > imin ) then
            ! Average cells by mass
            vol = 0.0
            mass = 0.0
            masstemp = 0.0
            masshe4 = 0.0
            massc12 = 0.0
            massne22 = 0.0
            do i = imin, imax
                r1 = max( sim_wd_rad_tab(i-1),   r_inner )
                r2 = min( sim_wd_rad_tab(i), r_outer )
                dr = r2 - r1
                dvol = dr * ( 3.0*r1*r2 + dr*dr ) * 4.0*PI/3.0
                dmass = dvol * sim_wd_dens_tab(i)
                vol = vol + dvol
                mass = mass + dmass
                masstemp = masstemp + dmass*sim_wd_temp_tab(i)
                masshe4  = masshe4 + dmass*sim_wd_he4_tab(i)
                massc12  = massc12  + dmass*sim_wd_c12_tab(i)
                massne22 = massne22 + dmass*sim_wd_o16_tab(i)
            end do
            dens = mass / vol 
            temp = masstemp / mass
            xhe4 = masshe4 / mass
            xc12 = massc12 / mass
            xo16 = massne22 / mass
         else if ( volume > sim_wd_volume ) then 
            dens = sim_wd_dens_tab(sim_wd_npnts)
            temp = sim_wd_temp_tab(sim_wd_npnts)
            xhe4 = sim_wd_he4_tab(sim_wd_npnts)
            xc12 = sim_wd_c12_tab(sim_wd_npnts)
            xo16 = sim_wd_o16_tab(sim_wd_npnts)
         else
            call interp1d_linear(sim_wd_vol_tab(1:sim_wd_npnts), sim_wd_dens_tab(1:sim_wd_npnts), volume, dens)
            call interp1d_linear(sim_wd_vol_tab(1:sim_wd_npnts), sim_wd_dens_tab(1:sim_wd_npnts)*&
                 sim_wd_temp_tab(1:sim_wd_npnts), volume, temp)
            call interp1d_linear(sim_wd_vol_tab(1:sim_wd_npnts), sim_wd_dens_tab(1:sim_wd_npnts)*&
                 sim_wd_he4_tab(1:sim_wd_npnts),  volume, xhe4)
            call interp1d_linear(sim_wd_vol_tab(1:sim_wd_npnts), sim_wd_dens_tab(1:sim_wd_npnts)*&
                 sim_wd_c12_tab(1:sim_wd_npnts),  volume, xc12)
            call interp1d_linear(sim_wd_vol_tab(1:sim_wd_npnts), sim_wd_dens_tab(1:sim_wd_npnts)*&
                 sim_wd_o16_tab(1:sim_wd_npnts),  volume, xo16)
            temp = temp / dens
            xhe4 = xhe4 / dens
            xc12 = xc12 / dens
            xo16 = xo16 / dens
         end if
  end if

  if ( dens < sim_smallrho .or. temp < sim_smallt ) then
     write(*,'(a,i5)')       '[sim_interpolate1dWd] Bad value(s) on PE=', sim_globalMe
     write(*,'(a,2i5)')      '                  imin, imax = ', imin, imax
     write(*,'(a,2es15.7)')  '            r_inner, r_outer = ', r_inner, r_outer
     write(*,'(a,30es15.7)') ' sim_wd_rad_tab(imin-1:imax) = ', sim_wd_rad_tab(imin-1:imax)
     write(*,'(a,1es15.7)')  '                      volume = ', volume
     write(*,'(a,4es15.7)')  '          dens,temp,mass,vol = ', dens, temp, mass, vol
     call Driver_abortFlash("[sim_interpolate1dWd] Bad values in initialization")
  end if

  return
end subroutine
