!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestInit/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!  A subset of simulation configuration data is loaded into AMReX at
!!  initialization and is therefore owned by AMReX.  As a result, AMReX is used
!!  to provide these data values to client code through the associated public
!!  Grid_* and local gr_* interface accessor routines.
!!
!!  This code tests that AMReX is properly initialized for a Cartesian domain by
!!  verifying correct results as obtained through the accessor routines.
!!
!! NOTES
!!  This simulation *must* be configured with at least the following
!!  2D run:
!!     ./setup -auto -2d -nxb=8 -nyb=4 
!!              unitTest/Grid/Amrex/TestInit 
!!             +noio -index-reorder
!!  3D run:
!!     ./setup -auto -3d -nxb=8 -nyb=4 -nzb=2
!!              unitTest/Grid/Amrex/TestInit 
!!             +noio -index-reorder
!!
!!  For the future:
!!             -unit=IO/IOMain/hdf5/serial/AM
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Driver_evolveFlash()
    use flash_iterator,        ONLY : flash_iterator_t
    use flash_tile,            ONLY : flash_tile_t
    use Grid_interface,        ONLY : Grid_getTileIterator, &
                                      Grid_releaseTileIterator
    use Grid_data,             ONLY : gr_enableTiling, &
                                      gr_tileSize
    use ut_testDriverMod

    implicit none

    integer :: cnt

    type(flash_iterator_t) :: itor
    type(flash_tile_t)     :: tileDesc

    call start_test_run

    ! Confirm contents of parfile
    call assertTrue(gr_enableTiling, "Tiling not enabled")
    call assertEqual(gr_tileSize(IAXIS), 4, "Incorrect X tile size")
    call assertEqual(gr_tileSize(JAXIS), 2, "Incorrect Y tile size")
    call assertEqual(gr_tileSize(KAXIS), 1, "Incorrect Z tile size")

    call assertEqual(NXB, 8, "Invalid number of cells/block in X")
    call assertEqual(NYB, 8, "Invalid number of cells/block in Y")
    call assertEqual(NZB, 1, "Invalid number of cells/block in Z")

    call assertEqual(NGUARD, 2, "Invalid number of guardcells")

    ! Get total blocks by explicitly turning off tiling at iterator creation
    cnt = 0
    call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.)
    do while(itor%isValid())
        call itor%currentTile(tileDesc)

        associate(lo      => tileDesc%limits(LOW,  :), &
                  hi      => tileDesc%limits(HIGH, :), &
                  loGC    => tileDesc%limitsGC(LOW,  :), &
                  hiGC    => tileDesc%limitsGC(HIGH, :), &
                  blkloGC => tileDesc%blkLimitsGC(LOW,  :), &
                  blkhiGC => tileDesc%blkLimitsGC(HIGH, :))
            call assertEqual(NXB, hi(IAXIS) - lo(IAXIS) + 1, "Invalid tile x-length")
            call assertEqual(NYB, hi(JAXIS) - lo(JAXIS) + 1, "Invalid tile y-length")
            call assertEqual(NZB, hi(KAXIS) - lo(KAXIS) + 1, "Invalid tile z-length")

            call assertEqual(NXB+2*NGUARD, hiGC(IAXIS) - loGC(IAXIS) + 1, "Invalid tile x-length")
            call assertEqual(NYB+2*NGUARD, hiGC(JAXIS) - loGC(JAXIS) + 1, "Invalid tile y-length")
            call assertEqual(NZB,          hiGC(KAXIS) - loGC(KAXIS) + 1, "Invalid tile z-length")

            call assertTrue(ALL(loGC == blkLoGC), "blkLimitsGC low != limitsGC low for block")
            call assertTrue(ALL(hiGC == blkHiGC), "blkLimitsGC high != limitsGC high for block")
        end associate

        cnt = cnt + 1

        call itor%next()
    end do
    call Grid_releaseTileIterator(itor)
    call assertEqual(cnt, 2, "Incorrect number of blocks")

    cnt = 0
    call Grid_getTileIterator(itor, ALL_BLKS, tiling=.TRUE.)
    do while(itor%isValid())
        call itor%currentTile(tileDesc)

        associate(lo      => tileDesc%limits(LOW,  :), &
                  hi      => tileDesc%limits(HIGH, :), &
                  loGC    => tileDesc%limitsGC(LOW,  :), &
                  hiGC    => tileDesc%limitsGC(HIGH, :), &
                  blkloGC => tileDesc%blkLimitsGC(LOW,  :), &
                  blkhiGC => tileDesc%blkLimitsGC(HIGH, :))

            ! Confirm appropriate size of tile (the interior)
            call assertEqual(4, hi(IAXIS) - lo(IAXIS) + 1, "Invalid tile x-length")
            call assertEqual(2, hi(JAXIS) - lo(JAXIS) + 1, "Invalid tile y-length")
            call assertEqual(1, hi(KAXIS) - lo(KAXIS) + 1, "Invalid tile z-length")

            ! Confirm that the grown tile has been grown appropriately
            if      (lo(IAXIS) == 1) then
                call assertEqual(-1, loGC(IAXIS), "Invalid grown tile lo x")
                call assertEqual( 4, hiGC(IAXIS), "Invalid grown tile hi x")
            else if (lo(IAXIS) == 5) then
                call assertEqual( 5, loGC(IAXIS), "Invalid grown tile lo x")
                call assertEqual(10, hiGC(IAXIS), "Invalid grown tile hi x")
            end if

            if      (lo(JAXIS) == 1) then
                call assertEqual(-1, loGC(JAXIS), "Invalid grown tile lo y")
                call assertEqual( 2, hiGC(JAXIS), "Invalid grown tile hi y")
            else if (lo(JAXIS) == 3) then
                call assertEqual( 3, loGC(JAXIS), "Invalid grown tile lo y")
                call assertEqual( 4, hiGC(JAXIS), "Invalid grown tile hi y")
            else if (lo(JAXIS) == 5) then
                call assertEqual( 5, loGC(JAXIS), "Invalid grown tile lo y")
                call assertEqual( 6, hiGC(JAXIS), "Invalid grown tile hi y")
            else if (lo(JAXIS) == 7) then
                call assertEqual( 7, loGC(JAXIS), "Invalid grown tile lo y")
                call assertEqual(10, hiGC(JAXIS), "Invalid grown tile hi y")
            end if

            call assertEqual(1, loGC(KAXIS), "Invalid grown tile lo z")
            call assertEqual(1, hiGC(KAXIS), "Invalid grown tile lo z")

            ! The limits of the enclosing block should be the same
            call assertEqual(NXB+2*NGUARD, blkhiGC(IAXIS) - blkloGC(IAXIS) + 1, "Invalid tile x-length")
            call assertEqual(NYB+2*NGUARD, blkhiGC(JAXIS) - blkloGC(JAXIS) + 1, "Invalid tile y-length")
            call assertEqual(NZB,          blkhiGC(KAXIS) - blkloGC(KAXIS) + 1, "Invalid tile z-length")
        end associate

        cnt = cnt + 1

        call itor%next()
    end do
    call Grid_releaseTileIterator(itor)
    call assertEqual(cnt, 16, "Incorrect number of blocks")

    ! OVERLOAD RUNTIME PARAMETERS
    ! This is not an expected use, but it is available
    gr_enableTiling = .FALSE. 

    ! Get total blocks by disabling all tiling.  This should ignore
    ! our tiling parameter value
    cnt = 0
    call Grid_getTileIterator(itor, ALL_BLKS, tiling=.TRUE.)
    do while(itor%isValid())
        cnt = cnt + 1
        call itor%next()
    end do
    call Grid_releaseTileIterator(itor)
    call assertEqual(cnt, 2, "Incorrect number of blocks")
    
    call finish_test_run

end subroutine Driver_evolveFlash

