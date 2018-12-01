!!****if* source/Grid/GridMain/AMR/Amrex/Grid_zeroFluxData
!!
!! NAME
!!  Grid_zeroFluxData
!!
!! SYNOPSIS
!!  call Grid_zeroFluxData
!!
!! DESCRIPTION 
!!  Request that the Grid unit zero all flux data managed by the unit.
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Grid_zeroFluxData()
  use Grid_interface,       ONLY : Grid_getTileIterator, &
                                   Grid_releaseTileIterator
  use flash_iterator,       ONLY : flash_iterator_t
  use flash_tile,           ONLY : flash_tile_t

  implicit none

  type(flash_iterator_t) :: itor
  type(flash_tile_t)     :: tileDesc

  real, pointer :: fluxData(:, :, :, :)

  if(NFLUXES < 1)   RETURN

  nullify(fluxData)

  call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.)
  do while (itor%isValid())
    call itor%currentTile(tileDesc)

    call tileDesc%getDataPtr(fluxData, FLUXX)
    if (associated(fluxData)) then
        fluxData(:,:,:,:) = 0.0
    end if
    call tileDesc%releaseDataPtr(fluxData, FLUXX)

    call tileDesc%getDataPtr(fluxData, FLUXY)
    if (associated(fluxData)) then
        fluxData(:,:,:,:) = 0.0
    end if
    call tileDesc%releaseDataPtr(fluxData, FLUXY)

    call tileDesc%getDataPtr(fluxData, FLUXZ)
    if (associated(fluxData)) then
        fluxData(:,:,:,:) = 0.0
    end if
    call tileDesc%releaseDataPtr(fluxData, FLUXZ)

    call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
end subroutine Grid_zeroFluxData

