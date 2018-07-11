subroutine hy_rk_getGravAccel(blockID,limits)
  use Gravity_interface, ONLY : Gravity_accelOneRow
  use Hydro_data, ONLY : hy_grav

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: blockID
  integer, intent(IN) :: limits(LOW:HIGH,MDIM)
  integer :: dir, i,j,k
  hy_grav = 0.
#ifdef GRAVITY
  do k=limits(LOW,KAXIS),limits(HIGH,KAXIS)
     do j=limits(LOW,JAXIS),limits(HIGH,JAXIS)
#ifdef FLASH_GRAVITY_TIMEDEP
        ! For time-dependent gravity, we call the local acceleration routine
        ! since the GPOT will be extrapolated forward in time for the RK
        ! sub-stages. A cleaner way might be to pass a pointer to Gravity_accelOneRow
        ! telling it what data structure to use for computing the acceleration.
        call accelOneRow((/j,k/),IAXIS,blockID,&
             GRID_IHI_GC,hy_grav(IAXIS,:,j,k))
#else
        ! For time-independent gravity, just call the regular Gravity routine.
        call Gravity_accelOneRow((/j,k/),IAXIS,blockID,&
             GRID_IHI_GC,hy_grav(IAXIS,:,j,k))
#endif
     enddo
  enddo
#if NDIM>1
  do k=limits(LOW,KAXIS),limits(HIGH,KAXIS)
     do i=limits(LOW,IAXIS),limits(HIGH,IAXIS)
#ifdef FLASH_GRAVITY_TIMEDEP
        call accelOneRow((/i,k/),JAXIS,blockID,&
             GRID_JHI_GC,hy_grav(JAXIS,i,:,k))
#else
        call Gravity_accelOneRow((/i,k/),JAXIS,blockID,&
             GRID_JHI_GC,hy_grav(JAXIS,i,:,k))
#endif
     enddo
  enddo
#if NDIM==3
  do j=limits(LOW,JAXIS),limits(HIGH,JAXIS)
     do i=limits(LOW,IAXIS),limits(HIGH,IAXIS)
#ifdef FLASH_GRAVITY_TIMEDEP
        call accelOneRow((/i,j/),KAXIS,blockID,&
             GRID_KHI_GC,hy_grav(KAXIS,i,j,:))
#else
        call accelOneRow((/i,j/),KAXIS,blockID,&
             GRID_KHI_GC,hy_grav(KAXIS,i,j,:))
#endif
     enddo
  enddo
#endif
#endif
#endif /* GRAVITY */

contains

  subroutine accelOneRow(pos, sweepDir, blockID, numCells, grav)
    use Grid_interface, ONLY : Grid_getBlkPhysicalSize, Grid_getBlkIndexLimits
    use Hydro_data, ONLY : hy_starState

    implicit none

    integer, dimension(2), intent(in) :: pos
    integer, intent(in)               :: sweepDir, blockID,  numCells
    real, intent(inout)               :: grav(numCells)

    real            :: blockSize(MDIM)
    integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
    integer         :: ii, iimin, iimax
    real            :: gpot(numCells), delxinv
    integer         :: nxbBlock, nybBlock, nzbBlock

    !==================================================

#ifdef GPOT_VAR
    call Grid_getBlkPhysicalSize(blockID, blockSize)
    call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC)
    nxbBlock = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
    nybBlock = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
    nzbBlock = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1


    iimin   = 1
    iimax   = numCells
    grav(iimin:iimax) = 0.0

    !Get row of potential values and compute inverse of zone spacing
    if (sweepDir == SWEEP_X) then                     ! x-direction
       delxinv = real(nxbBlock) / blockSize(IAXIS)
       gpot(:) = hy_starState(GPOT_VAR,:,pos(1),pos(2))
    elseif (sweepDir == SWEEP_Y) then                 ! y-direction
       delxinv = real(nybBlock) / blockSize(JAXIS)
       gpot(:) = hy_starState(GPOT_VAR,pos(1),:,pos(2))
    else                                            ! z-direction
       delxinv = real(nzbBlock) / blockSize(KAXIS)
       gpot(:) = hy_starState(GPOT_VAR,pos(1),pos(2),:)
    endif

    !----------------------------------------------------------------------
    !               Compute gravitational acceleration
    !**************** first-order differences
    !                 preserves conservation
    delxinv = 0.5e0 * delxinv
    do ii = iimin+1, iimax-1
       grav(ii) = grav(ii) + delxinv * (gpot(ii-1) - gpot(ii+1))
    enddo

    grav(iimin) = grav(iimin+1)     ! this is invalid data - must not be used
    grav(iimax) = grav(iimax-1)
#else
    ! Assume constant gravitational acceleration

#endif /* GPOT_VAR */
  end subroutine accelOneRow

end subroutine hy_rk_getGravAccel
