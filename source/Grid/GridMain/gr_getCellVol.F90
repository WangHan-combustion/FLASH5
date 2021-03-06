!!****if* source/Grid/GridMain/gr_getCellVol
!!
!! NAME
!!  gr_getCellVol
!!
!! SYNOPSIS
!!
!!  gr_getCellVol(integer(IN) :: xb,
!!                integer(IN) :: xe,
!!                integer(IN) :: yb,
!!                integer(IN) :: ye,
!!                integer(IN) :: zb,
!!                integer(IN) :: ze,
!!                integer(IN) :: blockID,
!!                real(xb:xe,yb:ye,zb:ze) (OUT) :: dataBlock,
!!                integer (IN) :: indexing)
!!  
!! DESCRIPTION 
!!  
!!  This routine calculate the volume of 
!!  the specified cells (bounded by xb:xe,yb:ye,zb:ze) of a given block.
!!
!! ARGUMENTS 
!!
!!  xb        : starting cell along IAXIS
!!  xe        : last cell along IAXIS
!!  yb        : starting cell along JAXIS
!!  ye        : last cell along JAXIS
!!  zb        : starting cell along KAXIS
!!  ze        : last cell along KAXIS
!!  blockID   : my block number
!!  dataBlock : storage for returning calculated values
!!  indexing  : indexing convention used for the first six arguments;
!!              should be GLOBALIDX1 or EXTERIOR
!!
!!***

subroutine gr_getCellVol(xb,xe,yb,ye,zb,ze,block,dataBlock,indexing)
  use Driver_interface, ONLY : Driver_abortFlash
!  use Grid_interface, ONLY : Grid_getDeltas, Grid_getSingleCellVol, &
!                             Grid_getGeometry
  use block_metadata, ONLY : block_metadata_t

  implicit none
#include "Flash.h"
#include "constants.h"

  type(block_metadata_t),intent(IN) :: block
  integer,intent(IN) :: xb,xe,yb,ye,zb,ze
  real,dimension(xb:xe,yb:ye,zb:ze),intent(OUT)::dataBlock
  integer,intent(IN) :: indexing

  dataBlock(:,:,:) = 0.0
  call Driver_abortFlash("[gr_getCellVol] DEPRECATED. Use Grid_getCellVolumes")
 
!  integer :: geometry
!  integer,dimension(MDIM) :: point
!  real,dimension(MDIM) :: del
!  integer :: i, j, k
!  integer :: beginCount
!
!#ifdef DEBUG_GRID
!  print*,xb,xe,yb,ye,zb,ze
!#endif
!  call Grid_getGeometry(geometry)
!  if (geometry==CARTESIAN) then
!     call Grid_getDeltas(block%level, del)
!     dataBlock=del(IAXIS)
!     do i = 2,NDIM
!        dataBlock = dataBlock*del(i)
!     end do
!#ifdef DEBUG_GRID
!     print*,'in here it is ',maxval(dataBlock),del
!#endif
!  else
!     do k=zb,ze
!        point(KAXIS) = k
!        do j=yb,ye
!           point(JAXIS) = j
!           do i=xb,xe
!              point(IAXIS) = i
!              call Grid_getSingleCellVol(point, block%level, dataBlock(i,j,k))
!           end do
!        end do
!     end do
!  end if
end subroutine gr_getCellVol
