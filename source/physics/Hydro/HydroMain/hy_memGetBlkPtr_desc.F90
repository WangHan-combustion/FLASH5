!!****if* source/physics/Hydro/HydroMain/hy_memGetBlkPtr
!!
!! NAME
!!  hy_memGetBlkPtr
!!
!! SYNOPSIS
!!
!!  call hy_memGetBlkPtr(Grid_tile_t(IN)  :: tileDesc,
!!                 real(pointer)(:,:,:,:) :: dataPtr,
!!                 integer(IN),optional   :: gridDataStruct)
!!
!! DESCRIPTION
!!
!!  Gets a pointer to a single block of allocated scratch data that conforms with
!!  the specified Grid data structure.
!!
!!  The block data may include zero, one, or several layers of guard cells,
!!  dependent on how hy_memAllocScratch was called for the gridDataStruct.
!!  If the optional argument "gridDataStruct" is not specified,
!!  it returns a block from cell centered data structure.
!!
!! ARGUMENTS
!!
!!  blockDesc : describes the local block
!!
!!  dataPtr : Pointer to the data block
!!
!!  gridDataStruct : optional integer value specifying data structure.
!!                   The options are defined in constants.h and they are :
!!                   SCRATCH scratch space that can fit cell and face centered variables
!!                   SCRATCH_CTR scratch space for cell centered variables
!!                   SCRATCH_FACEX scratch space for facex variables
!!                   SCRATCH_FACEY scratch space for facey variables
!!                   SCRATCH_FACEZ scratch space for facez variables
!!
!!
!!
!! NOTES
!!
!!  hy_memGetBlkPtr is an accessor function that passes a pointer
!!  as an argument and requires an explicit interface.
!!
!!  If you call hy_memGetBlkPtr, you should also  call hy_memReleaseBlkPtr when you
!!  are done using the pointer.
!!
!!***

!!REORDER(5):hy_memArrayScratch, hy_memArrayScratch_ctr, hy_memArrayScratch_facevar[xyz]

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine hy_memGetBlkPtr_desc(tileDesc,dataPtr, gridDataStruct)

#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"

  use hy_memInterface, ONLY : hy_memGetBlkPtr
  use Grid_tile, ONLY : Grid_tile_t
  implicit none

  type(Grid_tile_t), intent(IN) :: tileDesc
  real, POINTER_INTENT_OUT :: dataPtr(:,:,:,:)
  integer, optional,intent(in) :: gridDataStruct

  real, dimension(:,:,:,:), pointer :: medPtr
  integer,dimension(MDIM+1) :: lo
#ifdef INDEXREORDER
  integer, parameter :: iX = 1
#else
  integer, parameter :: iX = 2
#endif

  call hy_memGetBlkPtr(tileDesc%id, medPtr, gridDataStruct)

  ! DEV: How to set this if we eventually have tiling with Paramesh?
  lo = lbound(medPtr)
  lo(iX:ix+MDIM-1) = lo(iX:ix+MDIM-1) + tileDesc%grownLimits(LOW,:) - 1

  dataPtr(lo(1):,lo(2):,lo(3):,lo(4):) => medPtr

  return
end subroutine hy_memGetBlkPtr_desc

