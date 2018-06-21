!!****if* source/Simulation/SimulationMain/StreamingSineWave/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(real,pointer :: solnData(:,:,:,:),
!!                            integer(IN)  :: blockDesc  )
!!
!!
!! DESCRIPTION
!!
!!   Initialize solution data in one block for a streaming sine wave
!!
!! ARGUMENTS
!!
!!  solnData  -        pointer to solution data
!!  blockDesc -        describes the block to initialize
!!
!! PARAMETERS
!!
!!  
!!***

!!REORDER(4): solnData


subroutine Simulation_initBlock(solnData,block)

  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getCellCoords
  use ProgramHeaderModule, ONLY : nE, nDOF, nNodesX, nNodesE
  use RadiationFieldsModule, ONLY : nSpecies, nCR

  use block_metadata, ONLY : block_metadata_t

  implicit none

#include "constants.h"
#include "Flash.h"
  
  real, dimension(:,:,:,:), pointer :: solnData
  type(block_metadata_t), intent(in) :: block

  logical, parameter :: useGuardCell = .TRUE.

  real, allocatable, dimension(:) :: xLeft, xCenter, xRight
  real, allocatable, dimension(:) :: yLeft, yCenter, yRight
  real, allocatable, dimension(:) :: zLeft, zCenter, zRight
  real, dimension(MDIM) :: delta
  real :: dx, dy, dz

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: iSize, iSizeGC
  integer :: jSize, jSizeGC
  integer :: kSize, kSizeGC

  integer :: meshGeom

  integer :: i, j, k, n
  integer :: is, im, ie, ixnode, iynode, iznode, ienode, ii_0, id, ii
  integer :: nx, ny, nz
  real :: xnode, ynode, znode, ss

  blkLimits = block%limits
  blkLimitsGC = block%limitsGC

  iSizeGC = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  jSizeGC = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  kSizeGC = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
  jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
  kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

  !! allocate all needed space
  allocate(xLeft(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
  allocate(xCenter(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
  allocate(xRight(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
  allocate(yLeft(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
  allocate(yCenter(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
  allocate(yRight(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
  allocate(zLeft(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
  allocate(zCenter(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
  allocate(zRight(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))

  xLeft = 0.0 ; xCenter(:) = 0.0 ; xRight(:) = 0.0
  yLeft = 0.0 ; yCenter(:) = 0.0 ; yRight(:) = 0.0
  zLeft = 0.0 ; zCenter(:) = 0.0 ; zRight(:) = 0.0

  call Grid_getDeltas(block%level, delta)
  dx = delta(IAXIS)
  dy = delta(JAXIS)
  dz = delta(KAXIS)

  call Grid_getGeometry(meshGeom)

  call Grid_getCellCoords(IAXIS, block, LEFT_EDGE,  useGuardCell, xLeft,   iSizeGC)
  call Grid_getCellCoords(IAXIS, block, CENTER,     useGuardCell, xCenter, iSizeGC)
  call Grid_getCellCoords(IAXIS, block, RIGHT_EDGE, useGuardCell, xRight,  iSizeGC)

  call Grid_getCellCoords(JAXIS, block, LEFT_EDGE,  useGuardCell, yLeft,   jSizeGC)
  call Grid_getCellCoords(JAXIS, block, CENTER,     useGuardCell, yCenter, jSizeGC)
  call Grid_getCellCoords(JAXIS, block, RIGHT_EDGE, useGuardCell, yRight,  jSizeGC)

  call Grid_getCellCoords(KAXIS, block, LEFT_EDGE,  useGuardCell, zLeft,   kSizeGC)
  call Grid_getCellCoords(KAXIS, block, CENTER,     useGuardCell, zCenter, kSizeGC)
  call Grid_getCellCoords(KAXIS, block, RIGHT_EDGE, useGuardCell, zRight,  kSizeGC)

  nz = nNodesE*nNodesX(1)*nNodesX(2)
  ny = nNodesE*nNodesX(1)
  nx = nNodesE

  if ( nDOF /= nz*nNodesX(3) ) then
     call Driver_abortFlash("nDOE /= nNodesX(1)*nNodesX(2)*nNodesX(3)*nNodesE")
  endif

  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           ! Initialize hydro data
           solnData(VELX_VAR,i,j,k) = sim_velx_i
           solnData(VELY_VAR,i,j,k) = sim_vely_i
           solnData(VELZ_VAR,i,j,k) = sim_velz_i
           solnData(DENS_VAR,i,j,k) = sim_dens_i
           solnData(TEMP_VAR,i,j,k) = sim_temp_i
           solnData(PRES_VAR,i,j,k) = sim_pres_i
           solnData(EINT_VAR,i,j,k) = sim_eint_i
           solnData(ENER_VAR,i,j,k) = sim_etot_i
           solnData(GAMC_VAR,i,j,k) = sim_gamc_i
           solnData(GAME_VAR,i,j,k) = sim_game_i
           do n = SPECIES_BEGIN,SPECIES_END
              solnData(n,i,j,k) = sim_xn_i(n)
           enddo

           ! Initialize neutrino data
           do is = 1, nSpecies ; do im = 1, nCR ; do ie = 1, nE

              ii_0 = (is-1)*(nCR*nE*nDOF) + (im-1)*(nE*nDOF) + (ie-1)*nDOF

              do iznode = 1, nNodesX(3) ; do iynode = 1, nNodesX(2) ; do ixnode = 1, nNodesX(1) ; do ienode = 1, nNodesE

                 ! calculate the indices
                 id = (ienode-1) + nx*(ixnode-1) + ny*(iynode-1) + nz*(iznode-1)
                 ii = THORNADO_BEGIN + ii_0 + id

                 ! calculate actual positions of the nodes used for the gaussian quadrature
                 xnode = xCenter(i) + (real(ixnode)-1.5)*dx/sqrt(3.0)
                 ynode = yCenter(j) + (real(iynode)-1.5)*dy/sqrt(3.0)
                 znode = zCenter(k) + (real(iznode)-1.5)*dz/sqrt(3.0)

                 ss = 1.0e0 + 0.9999e0 * sin(2.0e0*PI*xnode)

                 ! J moment, im = 1
                 if (im == 1) solnData(ii,i,j,k) = ss

                 ! H_x moment, im = 2
                 if (im == 2) solnData(ii,i,j,k) = 3.0e10 * 0.9999e0 * ss

                 ! H_y moment, im = 3
                 if (im == 3) solnData(ii,i,j,k) = 0.0e0

                 ! H_z moment, im = 4
                 if (im == 4) solnData(ii,i,j,k) = 0.0e0

              enddo ; enddo ; enddo ; enddo
           enddo ; enddo ; enddo

        enddo
     enddo
  enddo

  ! cleanup
  deallocate(xLeft,xCenter,xRight)
  deallocate(yLeft,yCenter,yRight)
  deallocate(zLeft,zCenter,zRight)

  return
end subroutine Simulation_initBlock
