!!****if* source/physics/Hydro/HydroMain/Spark/Hydro_correctGravSrc
!!
!!
!! NAME
!!
!!  Hydro_correctGravSrc
!!
!!
!! SYNOPSIS
!!
!!  Hydro_correctGravSrc(integer(IN) :: blockCount,
!!                       integer(IN) :: blockList(blockCount)
!!                       real(IN)    :: dt )
!!
!!
!! DESCRIPTION
!!
!!  Corrects gravitational source terms with n+1 potential to make source
!!  source terms second-order accurate
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  dt         - timestep
!!
!!***
subroutine Hydro_correctGravSrc(nblk,blklst,dt)
  use Hydro_data, ONLY : hy_grav
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr
  implicit none

#include "Flash.h"
#include "constants.h"
#include "Spark.h"

  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)
  real,    intent(in) :: dt

  integer :: blockID
  integer :: n, i, j, k
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real, pointer :: solnData(:,:,:,:)
  real, dimension(MDIM) :: momOld, momNew
  real :: ekin

  do n=1,nblk
     blockID = blklst(n)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData)
     call hy_rk_getGravAccel(blockID,blkLimitsGC)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              momOld = solnData(DENS_VAR,i,j,k)*solnData(VELX_VAR:VELZ_VAR,i,j,k)
              momNew = momOld + 0.5*dt*solnData(DENS_VAR,i,j,k)*hy_grav(:,i,j,k)
              solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) &
                   + 0.5*dt*dot_product(momNew, hy_grav(:,i,j,k))/solnData(DENS_VAR,i,j,k)
              solnData(VELX_VAR:VELZ_VAR,i,j,k) = momNew/solnData(DENS_VAR,i,j,k)
              ! ekin = 0.5*dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k), &
              !      solnData(VELX_VAR:VELZ_VAR,i,j,k))
              ! print *, (solnData(EINT_VAR,i,j,k) - (solnData(ENER_VAR,i,j,k)-ekin))/solnData(EINT_VAR,i,j,k)
           end do
        end do
     end do

  end do

  return
end subroutine Hydro_correctGravSrc
