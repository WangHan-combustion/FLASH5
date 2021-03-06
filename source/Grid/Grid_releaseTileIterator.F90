!!****if* source/Grid/Grid_releaseTileIterator
!!
!! NAME
!!  Grid_releaseTileIterator
!!
!! SYNOPSIS
!!  Grid_releaseTileIterator(Grid_iterator_t(INOUT) :: itor)
!!  
!! DESCRIPTION 
!!  Destroy given block/tile iterator.
!!
!! ARGUMENTS 
!!  itor - the block/tile iterator to destroy.
!!
!! SEE ALSO
!!  Grid_getTileIterator
!!
!!***

subroutine Grid_releaseTileIterator(itor)
  use Grid_iterator, ONLY : Grid_iterator_t

  implicit none

  type(Grid_iterator_t), intent(INOUT) :: itor

  return
end subroutine Grid_releaseTileIterator

