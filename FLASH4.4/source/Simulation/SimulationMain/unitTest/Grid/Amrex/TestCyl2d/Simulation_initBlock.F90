!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestInit/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer (IN) ::blockId, 
!!
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes the Grid with a composit number which is a combination
!!  of the block number and the indices of the cell
!! 
!! ARGUMENTS
!!
!!  blockId -          the blockId to update
!!  
!!
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Simulation_initBlock(initData, tileDesc)
    use flash_tile, ONLY : flash_tile_t

    implicit none

    real,                           pointer :: initData(:, :, :, :)
    type(flash_tile_t), intent(IN)          :: tileDesc

    integer :: i
    integer :: j
    integer :: k
    integer :: var

    associate(lo => tileDesc%limitsGC(LOW,  :), &
              hi => tileDesc%limitsGC(HIGH, :))
        do           var = UNK_VARS_BEGIN, UNK_VARS_END
            do         k = lo(KAXIS), hi(KAXIS)
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        initData(i, j, k, var) = 1.1 * var
                    end do
                end do
            end do
        end do
    end associate
end subroutine Simulation_initBlock

