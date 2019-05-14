!!****if* source/physics/Hydro/HydroMain/Spark/rk3/Hydro
!!
!!
!! NAME
!!
!!  Hydro
!!
!!
!! SYNOPSIS
!!
!!  Hydro(integer(IN) :: blockCount,
!!        integer(IN) :: blockList(blockCount)
!!        real(IN)    :: timeEndAdv,
!!        real(IN)    :: dt,
!!        real(IN)    :: dtOld,
!!        integer(IN) :: sweepOrder)
!!
!!
!! DESCRIPTION
!!
!!  Performs physics update in a directionally unsplit fashion.
!!
!!  The blockList and blockCount arguments tell this routine on
!!  which blocks and on how many to operate.  blockList is an
!!  integer array of size blockCount that contains the local
!!  block numbers of blocks on which to advance.
!!
!!  dt gives the timestep through which this update should advance,
!!  and timeEndAdv tells the time that this update will reach when
!!  it finishes.  dtOld gives the previously taken timestep.
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  timeEndAdv - end time
!!  dt         - timestep
!!  dtOld      - old timestep
!!  sweepOrder - dummy argument for the unsplit scheme, just a dummy
!!               variable to be consistent with a toplayer stub function
!!
!!***

subroutine Hydro( nblk, blklst, &
     timeEndAdv, dt,  dtOld,&
     sweepOrder)

  use Hydro_data, ONLY : hy_useHydro, hy_fluxCorrect, hy_gcMask, hy_lChyp, &
       hy_C_hyp, hy_globalComm
  use hy_rk_interface, ONLY : hy_rk_eos, hy_rk_getFaceFlux, hy_rk_getGravAccel,&
       hy_rk_updateSoln
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_conserveFluxes, &
       Grid_fillGuardCells, Grid_getBlkIndexLimits
  use Eos_interface, ONLY : Eos_wrapped
  use IO_interface, ONLY : IO_setScalar

  implicit none

#include "Flash.h"
#include "constants.h"
  include "Flash_mpi.h"

  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)
  real,    intent(in) :: timeEndAdv, dt, dtOld
  integer, intent(IN) :: sweepOrder

  integer :: blockID
  integer :: n, error
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(LOW:HIGH,MDIM) :: limits
  real :: hdt
  real, dimension(3) :: coeffs
  real, parameter :: onethird = 1./3.
  real, parameter :: twothird = 2./3.
  real, parameter :: oneSixth = 1./6.

  if (.NOT. hy_useHydro) return

  call Timers_start("Hydro")

  hdt = 0.5*dt

  ! Find the global maximum hyperbolic speed. hy_lChyp from Hydro_computeDt
#ifdef SPARK_GLM
  call MPI_AllReduce (hy_lChyp, hy_C_hyp, 1, &
       FLASH_REAL, MPI_MAX, hy_globalComm, error)
  call IO_setScalar("C_hyp", hy_lChyp)
#endif

  call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.false.,maskSize=NUNK_VARS,mask=hy_gcMask)

  ! Loop over blocks and compute M1 update block-by-block
  do n=1,nblk
     blockID = blklst(n)

     ! DivB will technically be lagged by 1 step, but we need ghost zones to
     ! compute the gradients. I ain't doing more communication for a diagnostic...
     call calcDivB(blockID)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     call setLims(limits, NGUARD-1)
     call shockDetect(blockID,limits)

     ! Setup scratch storage of block data
     call Timers_start("scratch")
     call saveState(blockID)
     call Timers_stop("scratch")

     !*********************
     ! Stage 1:
     !*********************

     ! calculate gravitational acceleration based on current value of GPOT_VAR
     ! This is stored in module-scope variable hy_grav
     call hy_rk_getGravAccel(blockID,limits)

     ! Perform reconstruction and flux calculation
     ! In Stage 1, compute low-side fluxes and update for 2 guardcells
     call setLims(limits, 2*NSTENCIL)
     call Timers_start("getFaceFlux")
     call hy_rk_getFaceFlux(blockID, limits)
     call Timers_stop("getFaceFlux")

     if (hy_fluxCorrect) call addFluxes(oneSixth, .false.)

     ! Now update solution based on conservative fluxes
     ! U2 = C1*U0 + C2*U1 + C3*dt*L(U1)
     coeffs = (/1.0, 0.0, 1.0/)
     call Timers_start("updateSoln")
     call hy_rk_updateSoln(blockID,dt,dtOld,limits,coeffs)
     call Timers_stop("updateSoln")

     ! Update EOS based on intermediate solution
     call Timers_start("eos")
     call hy_rk_eos(limits)
     call Timers_stop("eos")

     !*********************
     ! Stage 2:
     !*********************

     ! calculate gravitational acceleration based on current value of GPOT_VAR
     ! This is stored in module-scope variable hy_grav
     call hy_rk_getGravAccel(blockID,limits)

     ! Perform reconstruction and flux calculation
     ! In Stage 1, compute low-side fluxes and update for 0 guardcells
     call setLims(limits, NSTENCIL)
     call Timers_start("getFaceFlux")
     call hy_rk_getFaceFlux(blockID, limits)
     call Timers_stop("getFaceFlux")

     if (hy_fluxCorrect) call addFluxes(oneSixth, .true.)

     ! Now update solution based on conservative fluxes
     ! U2 = C1*U0 + C2*U1 + C3*dt*L(U1)
     coeffs = (/0.75, 0.25, 0.25/)
     call Timers_start("updateSoln")
     call hy_rk_updateSoln(blockID,dt,dtOld,limits,coeffs)
     call Timers_stop("updateSoln")

     ! Update EOS based on intermediate solution
     call Timers_start("eos")
     call hy_rk_eos(limits)
     call Timers_stop("eos")

     !*********************
     ! Stage 3:
     !*********************

     ! calculate gravitational acceleration based on current value of GPOT_VAR
     ! This is stored in module-scope variable hy_grav
     call hy_rk_getGravAccel(blockID,limits)

     ! Perform reconstruction and flux calculation
     ! In Stage 1, compute low-side fluxes and update for 0 guardcells
     call setLims(limits, 0)
     call Timers_start("getFaceFlux")
     call hy_rk_getFaceFlux(blockID, limits)
     call Timers_stop("getFaceFlux")

     if (hy_fluxCorrect) call addFluxes(twothird, .true.)

     ! Now update solution based on conservative fluxes
     ! U2 = C1*U0 + C2*U1 + C3*dt*L(U1)
     coeffs = (/onethird, twothird, twothird/)
     call Timers_start("updateSoln")
     call hy_rk_updateSoln(blockID,dt,dtOld,limits,coeffs)
     call Timers_stop("updateSoln")

     ! Update EOS based on intermediate solution
     call Timers_start("eos")
     call hy_rk_eos(limits)
     call Timers_stop("eos")

     ! Finally, store the output and free up the scratch array
     call updateState(blockID)

  end do

  if (hy_fluxCorrect) then
     ! Call the PARAMESH routine to correct fluxes at fine-coarse boundaries
     call Timers_start("flux correct")
     call Grid_conserveFluxes(ALLDIR,0)
     ! Loop over blocks and correct block-edge solution
     do n=1,nblk
        blockID = blklst(n)
        call hy_rk_correctFluxes(blockID,dt)
      end do
     call Timers_stop("flux correct")
  end if

  ! Reset local maximum hyperbolic speed. This will be updated in Hydro_computeDt.
  hy_lChyp = TINY(1.0)

  call Timers_stop("Hydro")

contains

#include "Hydro_funcs.F90"

end subroutine Hydro
