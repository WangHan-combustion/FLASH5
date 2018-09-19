!!****if* source/Simulation/SimulationMain/unitTest/PFFT_PoissonFD/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!!  Reference:
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!  
!!
!! 
!!
!!***

!!-----!! Do not REORDER(4): solnData

subroutine Simulation_initBlock(solnData,block)

!  use Simulation_data
  use Simulation_data, ONLY :sim_xMin,sim_xMax,sim_yMin,sim_yMax,sim_zMin,sim_zMax
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkBoundBox, Grid_getBlkCenterCoords, Grid_getDeltas
  use block_metadata, ONLY : block_metadata_t
  use amrex_fort_module,         ONLY : wp => amrex_real

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  real,dimension(:,:,:,:),pointer :: solnData
  type(block_metadata_t), intent(in) :: block
  integer :: blockID
  !!$ ---------------------------------
 
  integer :: i, j, k
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  real,allocatable, dimension(:) ::xCenter,yCenter,zCenter
  integer :: sizeX,sizeY,sizeZ

  real :: Lx, Ly, Lz, xi, yi, zi, Phi_ijk, F_ijk, xf, yf, zf


  real, parameter :: pfb_waven_x = 2.
  real, parameter :: pfb_waven_y = 1.
  real, parameter :: pfb_waven_z = 2.
  real, parameter :: pfb_alpha_x = 0.

  logical :: gcell = .true.

  real(wp), contiguous, pointer :: facexData(:,:,:,:), faceyData(:,:,:,:), facezData(:,:,:,:)

  real :: del(MDIM)

  !----------------------------------------------------------------------
  blkLimits = block%limits
  blkLimitsGC = block%limitsGC
  allocate(xCenter(blkLimitsGC(LOW, IAXIS):blkLimitsGC(HIGH, IAXIS)))
  allocate(yCenter(blkLimitsGC(LOW, JAXIS):blkLimitsGC(HIGH, JAXIS)))
  allocate(zCenter(blkLimitsGC(LOW, KAXIS):blkLimitsGC(HIGH, KAXIS)))
  xCenter = 0.0
  yCenter = 0.0
  zCenter = 0.0

  sizeX = SIZE(xCenter)
  sizeY = SIZE(yCenter)
  sizeZ = SIZE(zCenter)
 
  call Grid_getDeltas(block%level,del)
 
  call Grid_getCellCoords(IAXIS, block, CENTER, gcell, xCenter, sizeX)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS, block, CENTER, gcell, yCenter, sizeY)
  if (NDIM == 3) call Grid_getCellCoords(KAXIS, block, CENTER, gcell, zCenter, sizeZ)

#ifdef DEBUG_SIMULATION
98 format('initBlock:',A4,'(',I3,':   ,',   I3,':   ,',   I3,':   ,',   I3,':   )')
99 format('initBlock:',A4,'(',I3,':',I3,',',I3,':',I3,',',I3,':',I3,',',I3,':',I3,')')
  print 99,"solnData" ,(lbound(solnData ,i),ubound(solnData ,i),i=1,4)
  print*,'blkLim  :',blkLimits
  print*,'blkLimGC:',blkLimitsGC
#endif
!------------------------------------------------------------------------------

  Lx = sim_xMax - sim_xMin
  Ly = sim_yMax - sim_yMin  
  Lz = sim_zMax - sim_zMin

  call Grid_getBlkPtr(block,facexData,FACEX)
  call Grid_getBlkPtr(block,faceyData,FACEY)
  call Grid_getBlkPtr(block,facezData,FACEZ)

  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

           xi=xCenter(i)
           yi=yCenter(j)
           zi=zCenter(k)

           xf = xi - 0.5*del(IAXIS)
           yf = yi - 0.5*del(JAXIS)
           zf = zi - 0.5*del(KAXIS)

           facexData(i,j,k,VELC_FACE_VAR) =  1.0*cos(2*PI*xf)*sin(2*PI*yi)*sin(2*PI*zi)
           faceyData(i,j,k,VELC_FACE_VAR) = -0.5*sin(2*PI*xi)*cos(2*PI*yf)*sin(2*PI*zi)
           facezData(i,j,k,VELC_FACE_VAR) = -0.5*sin(2*PI*xi)*sin(2*PI*yi)*cos(2*PI*zf)

           facexData(i,j,k,VELA_FACE_VAR) = facexData(i,j,k,VELC_FACE_VAR)
           faceyData(i,j,k,VELA_FACE_VAR) = faceyData(i,j,k,VELC_FACE_VAR)
           facezData(i,j,k,VELA_FACE_VAR) = facezData(i,j,k,VELC_FACE_VAR) 

        enddo
     enddo
  enddo


  ! set values for u,v velocities and pressure
  solnData(:,:,:,DENS_VAR) = 1.0

  call Grid_releaseBlkPtr(block,facexData,FACEX)
  call Grid_releaseBlkPtr(block,faceyData,FACEY)
  call Grid_releaseBlkPtr(block,facezData,FACEZ)


!!$  write(*,*) 'BlockID=',blockID
!!$  write(*,*) 'Center coordinates=',coord
!!$  write(*,*) 'Size =',bsize
 

  ! Release pointer
!  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  return

end subroutine Simulation_initBlock
