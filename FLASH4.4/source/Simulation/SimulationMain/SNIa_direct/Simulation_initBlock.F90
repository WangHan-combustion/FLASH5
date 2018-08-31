!!
!! Adapted from Dean M. Townsley (2009) SNIa_ddt setup
!!
!! This initializes the white dwarf and flame for Type Ia simulation
!!

subroutine Simulation_initBlock(solnData, block)
  
  use Simulation_data
  use sim_local_interface, ONLY : sim_interpolate1dWd
  use Driver_interface, ONLY : Driver_abortFlash
  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
    Multispecies_getSumFrac
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkBoundBox, Grid_getDeltas, Grid_putPointData, &
    Grid_getCellCoords, Grid_getBlkBoundBox
  use Eos_interface, ONLY : Eos
  use block_metadata, ONLY : block_metadata_t
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Multispecies.h"
  
  real,dimension(:,:,:,:),pointer :: solnData
  type(block_metadata_t), intent(in) :: block

  logical, parameter :: useGuardCell = .true.
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: iSizeGC, jSizeGC, kSizeGC
  real, allocatable, dimension(:) :: xCenter, yCenter, zCenter
  real, allocatable, dimension(:) :: xLeft, yLeft, zLeft
  real, allocatable, dimension(:) :: xRight, yRight, zRight
  real, dimension(MDIM) :: delta
  integer, dimension(MDIM) :: cell
  integer :: meshGeom
  !real, dimension(2,MDIM) :: Bspan

  real, dimension(EOS_NUM) :: eosData
  real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction

  real :: radCenter, thtCenter, radCenterVol, dVol, ign_dist, dx, dy, dz, dr
  real :: radMin, radMax, x2Min, x2Max, y2Min, y2Max, z2Min, z2Max
  real :: xhe4initial, xc12initial, xo16initial
  real :: ytot, atot, ztot, abar, zbar
  real :: velx, vely, velz, dens, temp, pres, eint, etot, gamc, game
  integer :: i, j, k, n, level

!==============================================================================

  call Grid_getGeometry(meshGeom)

  ! Get the indices of the blocks
  blkLimits=block%Limits
  blkLimitsGC=block%LimitsGC
  iSizeGC = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  jSizeGC = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  kSizeGC = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(xCenter(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
  allocate(xLeft(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
  allocate(xRight(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
  allocate(yCenter(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
  allocate(yLeft(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
  allocate(yRight(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
  allocate(zCenter(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
  allocate(zLeft(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
  allocate(zRight(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))

  call Grid_getDeltas(block%level, delta)
  dx = delta(IAXIS)
  dy = delta(JAXIS)
  dz = delta(KAXIS)

  call Grid_getCellCoords(IAXIS,block,CENTER,    useGuardCell,xCenter,iSizeGC)
  call Grid_getCellCoords(IAXIS,block,LEFT_EDGE, useGuardCell,xLeft,  iSizeGC)
  call Grid_getCellCoords(IAXIS,block,RIGHT_EDGE,useGuardCell,xRight, iSizeGC)

  call Grid_getCellCoords(JAXIS,block,CENTER,    useGuardCell,yCenter,jSizeGC)
  call Grid_getCellCoords(JAXIS,block,LEFT_EDGE, useGuardCell,yLeft,  jSizeGC)
  call Grid_getCellCoords(JAXIS,block,RIGHT_EDGE,useGuardCell,yRight, jSizeGC)

  call Grid_getCellCoords(KAXIS,block,CENTER,    useGuardCell,zCenter,kSizeGC)
  call Grid_getCellCoords(KAXIS,block,LEFT_EDGE, useGuardCell,zLeft,  kSizeGC)
  call Grid_getCellCoords(KAXIS,block,RIGHT_EDGE,useGuardCell,zRight, kSizeGC)
  
  !call Grid_getBlkBoundBox(blockID, Bspan)
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

           !-----------------------------------------------
           !  determine state of material at this radius if unburned from external 1-d hyrdostatic model
           !-----------------------------------------------

           if ( meshGeom == SPHERICAL ) then
              radCenter = xCenter(i)
              thtCenter = yCenter(j)
              dVol = (4.0*PI/3.0)*(xRight(i)-xLeft(i))*(3.0*xLeft(i)*xRight(i) + (xRight(i)-xLeft(i))**2)
              radCenterVol = (4.0*PI/3.0)*xLeft(i)**3 + 0.5*dVol

              radMin = xLeft(i)
              radMax = xRight(i)
           else if ( meshGeom == CYLINDRICAL ) then
              radCenter = sqrt(xCenter(i)**2 + yCenter(j)**2)
              if ( yCenter(j) /= 0.0 ) then
                 thtCenter = atan( xCenter(i) / yCenter(j) )
              else
                 thtCenter = 0.5 * PI
              end if
              if ( thtCenter < 0.0 ) then
                 thtCenter = thtCenter + PI
              end if
              radCenterVol = (4.0*PI/3.0) * radCenter**3

              x2Min = xLeft(i)**2
              x2Max = xRight(i)**2

              if ( yLeft(j)*yRight(j) > 0.0 ) then
                 y2Min = min(yLeft(j)**2,yRight(j)**2)
              else
                 y2Min = 0.0
              end if
              y2Max = max(yLeft(j)**2,yRight(j)**2)

              radMin = sqrt( x2Min + y2Min )
              radMax = sqrt( x2Max + y2Max )
           else if ( meshGeom == CARTESIAN ) then
              radCenter = sqrt(xCenter(i)**2 + yCenter(j)**2 + zCenter(k)**2)
              thtCenter = acos( zCenter(k) / radCenter )
              radCenterVol = (4.0*PI/3.0) * radCenter**3

              if ( xLeft(i)*xRight(i) > 0.0 ) then
                 x2Min = min(xLeft(i)**2,xRight(i)**2)
              else
                 x2Min = 0.0
              end if
              x2Max = max(xLeft(i)**2,xRight(i)**2)

              if ( yLeft(j)*yRight(j) > 0.0 ) then
                 y2Min = min(yLeft(j)**2,yRight(j)**2)
              else
                 y2Min = 0.0
              end if
              y2Max = max(yLeft(j)**2,yRight(j)**2)

              if ( zLeft(k)*zRight(k) > 0.0 ) then
                 z2Min = min(zLeft(k)**2,zRight(k)**2)
              else
                 z2Min = 0.0
              end if
              z2Max = max(zLeft(k)**2,zRight(k)**2)

              radMin = sqrt( x2Min + y2Min + z2Min )
              radMax = sqrt( x2Max + y2Max + z2Max )
           else
              call Driver_abortFlash("Geometry not supported")
           end if
           
           call sim_interpolate1dWd(radCenterVol, radMin, radMax, dens, temp, xhe4initial, xc12initial, xo16initial)
           
           ! interp1D sorts the edge of the wd and shell, so just assign the values and the desired filler (Ne20 here)
           if ( radCenterVol < sim_wd_vol_tab(sim_wd_npnts) ) then
               massFraction(:) = sim_smallx
               if ( HE4_SPEC > 0 ) massFraction(HE4_SPEC) = max(xhe4initial,sim_smallx)
               if ( C12_SPEC > 0 ) massFraction(C12_SPEC) = max(xc12initial,sim_smallx)
               if ( O16_SPEC > 0 ) massFraction(O16_SPEC) = max(xo16initial,sim_smallx)
               if ( NE20_SPEC > 0 ) massFraction(NE20_SPEC) = max(1.0-xhe4initial-xc12initial-xo16initial,sim_smallx)
           else  ! on the fluff, fill with helium
               !dens = sim_densFluff
               !temp = sim_tempFluff
               massFraction(:) = sim_smallx
               if ( C12_SPEC > 0 ) massFraction(C12_SPEC) = max(sim_xc12Fluff,sim_smallx)
               if ( O16_SPEC > 0 ) massFraction(O16_SPEC) = max(sim_xo16Fluff,sim_smallx)
               if ( NI56_SPEC > 0 ) massFraction(NI56_SPEC) = max(sim_xni56Fluff,sim_smallx)
               if ( HE4_SPEC > 0 ) massFraction(HE4_SPEC) = max(1.0-sim_xo16Fluff-sim_xc12Fluff-sim_xni56Fluff, sim_smallx)
           endif
           
           ! get z and a bar for gamma calc below
           call Multispecies_getSumInv(A,ytot,massFraction)
           ! atot must be 1.0 so this call is redundant
		   !call Multispecies_getSumFrac(A,atot,massFraction)
		   abar = 1.0 / ytot
           call Multispecies_getSumFrac(Z,ztot,massFraction)
           zbar = ztot * abar

           if (sim_ignite) then
              !-----------------------------------------------
              ! initialize burned region
              !-----------------------------------------------
              ! default to a spherical region centered at specified coordinates
              ! distance from center of ignition region
              ign_dist = (xCenter(i) - sim_ignX)**2
              if ( NDIM >= 2 ) ign_dist = ign_dist + (yCenter(j) - sim_ignY)**2
              if ( NDIM == 3 ) ign_dist = ign_dist + (zCenter(k) - sim_ignZ)**2
              ign_dist = sqrt(ign_dist)

              ! heat to ignition temp if zone center is within half zone-width of match
              if ( ign_dist <= sim_ignROuter + 0.5*dx .and. ign_dist >= sim_ignRInner - 0.5*dx ) then
                 temp = sim_ignTOuter
              end if

           end if ! sim_ignite
           
           
           eosData(:) = 0.0
           eosData(EOS_ABAR) = abar
           eosData(EOS_ZBAR) = zbar
           eosData(EOS_DENS) = dens
           eosData(EOS_TEMP) = temp


           ! Giant traffic cone
           if ( dens < sim_smallrho .or. temp < sim_smallt) then
!               write(*,'(2(a,i5),a,3i3)') '[Before EOS] Bad value(s) on PE=',myPE,', blockID=',blockID,', (i,j,k)=',i,j,k
              write(*,'(a,2es15.7)') '  radCenter, thtCenter = ', radCenter, thtCenter
              write(*,'(a,2es15.7)') '     radMin,    radMax = ', radMin, radMax
              write(*,'(a,1es15.7)') '          radCenterVol = ', radCenterVol
              write(*,'(a,3es15.7)') '        x,y,z (center) = ', xCenter(i), yCenter(j), zCenter(k)
              write(*,'(a,3es15.7)') '        x,y,z   (left) = ', xLeft(i),   yLeft(j),   zLeft(k)
              write(*,'(a,3es15.7)') '        x,y,z  (right) = ', xRight(i),  yRight(j),  zRight(k)
              write(*,'(a,3es15.7)') '              dx,dy,dz = ', dx, dy, dz
              write(*,'(a,4es15.7)') '             dens,temp = ', dens, temp, pres, eint
              call Driver_abortFlash("[Simulation_initBlock] Bad values BEFORE EOS call.")
           end if
           
           ! Update the rest of the thermodynamic state
           call Eos(MODE_DENS_TEMP, 1, eosData, massFraction)
           
           !-----------------------------------------------
           !  Now store all this info on the grid
           !-----------------------------------------------
           velx = 0.0
           vely = 0.0
           velz = 0.0
           pres = eosData(EOS_PRES)
           eint = eosData(EOS_EINT)
           etot = eint + 0.5*(velx**2+vely**2+velz**2)
           gamc = eosData(EOS_GAMC)
           game = 1.0 + pres/eint/dens
           
           
           ! Giant traffic cone
           if ( dens < sim_smallrho .or. temp < sim_smallt .or. pres < sim_smallp .or. eint < sim_smalle ) then
!               write(*,'(2(a,i5),a,3i3)') '[After EOS] Bad value(s) on PE=',myPE,', blockID=',blockID,', (i,j,k)=',i,j,k
              write(*,'(a,2es15.7)') '  radCenter, thtCenter = ', radCenter, thtCenter
              write(*,'(a,2es15.7)') '     radMin,    radMax = ', radMin, radMax
              write(*,'(a,1es15.7)') '          radCenterVol = ', radCenterVol
              write(*,'(a,3es15.7)') '        x,y,z (center) = ', xCenter(i), yCenter(j), zCenter(k)
              write(*,'(a,3es15.7)') '        x,y,z   (left) = ', xLeft(i),   yLeft(j),   zLeft(k)
              write(*,'(a,3es15.7)') '        x,y,z  (right) = ', xRight(i),  yRight(j),  zRight(k)
              write(*,'(a,3es15.7)') '              dx,dy,dz = ', dx, dy, dz
              write(*,'(a,4es15.7)') '   dens,temp,pres,eint = ', dens, temp, pres, eint
              call Driver_abortFlash("[Simulation_initBlock] Bad values AFTER EOS call.")
           end if
           
           
           
           cell(IAXIS) = i-blkLimits(LOW,IAXIS)+1
           cell(JAXIS) = j-blkLimits(LOW,JAXIS)+1
           cell(KAXIS) = k-blkLimits(LOW,KAXIS)+1

           ! change state to eosData below
           call Grid_putPointData(block, CENTER, VELX_VAR, INTERIOR, cell, velx)
           call Grid_putPointData(block, CENTER, VELY_VAR, INTERIOR, cell, vely)
           call Grid_putPointData(block, CENTER, VELZ_VAR, INTERIOR, cell, velz)
           call Grid_putPointData(block, CENTER, DENS_VAR, INTERIOR, cell, dens)
           call Grid_putPointData(block, CENTER, TEMP_VAR, INTERIOR, cell, temp)
           call Grid_putPointData(block, CENTER, PRES_VAR, INTERIOR, cell, pres)
           call Grid_putPointData(block, CENTER, EINT_VAR, INTERIOR, cell, eint)
           call Grid_putPointData(block, CENTER, ENER_VAR, INTERIOR, cell, etot)
           call Grid_putPointData(block, CENTER, GAMC_VAR, INTERIOR, cell, gamc)
           call Grid_putPointData(block, CENTER, GAME_VAR, INTERIOR, cell, game)
           do n = SPECIES_BEGIN, SPECIES_END
              call Grid_putPointData(block, CENTER, n, INTERIOR, cell, massFraction(n))
           end do
        end do
     end do
  end do

  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yLeft)
  deallocate(yRight)
  deallocate(yCenter)
  deallocate(zLeft)
  deallocate(zRight)
  deallocate(zCenter)
  return
end subroutine Simulation_initBlock
