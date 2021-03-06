!!****if* source/Simulation/SimulationMain/Sedov/gr_expandDomain
!!
!!  NAME
!!     gr_expandDomain
!!
!!  SYNOPSIS
!!     call gr_expandDomain(logical(OUT) :: particlesInitialized)
!!
!!  DESCRIPTION
!!
!!    The grid is initialized in gr_createDomain, with a specified
!!    number of root blocks, typically one single block. This routine
!!    refines appropriate portions of the initialized physical
!!    domain according to the given refinement criteria, and applies
!!    initial conditions to the AMR domain.
!!
!!    In simulations with particles, under certain conditions particle
!!    positions will also be initialized.  Currently this is the case
!!    if and only if the runtime parameter refine_on_particle_count is
!!    true.
!!
!!  ARGUMENTS
!!    particlesInitialized : is true if this routine initialized particles positions
!!
!!  SIDE EFFECTS
!!
!!    Particle positions may be initialized, see DESCRIPTION above.
!!***

#define DEBUG_PARTICLES

subroutine gr_expandDomain (particlesInitialized)

  use Grid_interface, ONLY : Grid_copyF4DataToMultiFabs
  use Grid_data, ONLY : gr_domainBC,gr_eosModeInit,gr_refineOnParticleCount,&
       gr_meshMe,&
       gr_meshNumProcs, gr_lrefineMinInit, gr_gcellsUpToDate
  use Driver_data,         ONLY : dr_simGeneration
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_stampVarMask
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getLocalNumBlks, Grid_getListOfBlocks, &
       Grid_markRefineDerefine, Grid_getBlkPtr, Grid_releaseBlkPtr
  use gr_amrextInterface, ONLY : gr_fillMetaData, gr_amrextBuildMultiFabsFromF4Grid
  use tree, ONLY : lrefine, lrefine_min, lrefine_max, grid_changed
  use paramesh_interfaces, ONLY : amr_refine_derefine, amr_restrict
  use Eos_interface, ONLY : Eos_wrapped_blkid

#include "Flash.h"

#ifndef FLASH_GRID_PARAMESH2
  use physicaldata, ONLY: no_permanent_guardcells, mpi_pattern_id
#endif
  use Simulation_interface, ONLY : Simulation_initBlock
  use Particles_interface, ONLY : Particles_accumCount, &
    Particles_initPositions, &
    Particles_updateRefinement
  use Driver_interface, ONLY : Driver_abortFlash
  use RadTrans_interface, ONLY: RadTrans_sumEnergy
  use block_metadata, ONLY : block_metadata_t

  implicit none

#include "constants.h"

  include 'Flash_mpi.h'

  real, pointer:: solnData(:,:,:,:)
  logical, intent(out) :: particlesInitialized
  integer :: lnblocks, lrefineMinSave


  !!          Local variables and functions

  integer :: ntimes, i

  integer, dimension(2,MDIM) :: blkLimits
  integer, dimension(2,MDIM) :: blkLimitsGC
  integer :: count, cur_treedepth, grid_changed_anytime
  logical :: restart = .false.
  logical :: particlesPosnsDone, retainParticles
  integer ,dimension(MAXBLOCKS) :: blkList
  character(len=32), dimension(2,2) :: block_buff
  character(len=32)                 :: int_to_str
  integer :: gridDataStruct, whichBlocks
  type(block_metadata_t) :: blockDesc

  !!============================================================================



  !!============================================================================

  !! The beginning timestep number, time, and timestep.
  !! If the initial redshift (zinitial) is physical (>= 0),
  !! we use it to initialize the time; otherwise we set the
  !! redshift to zero and get the initial time from tinitial.
  !! The latter case (no cosmology) is the default.

  ! initialize the step counter and the simulation time
  ! the timestep initialization is moved to after the initialization,
  ! so we can check whether it is > t_cfl

  particlesInitialized=.false.
  call Grid_getBlkIndexLimits(1, blkLimits, blkLimitsGC)

  call gr_initParameshArrays(restart,        &
       gr_domainBC(LOW,IAXIS),gr_domainBC(HIGH,IAXIS), &
       gr_domainBC(LOW,JAXIS),gr_domainBC(HIGH,JAXIS), &
       gr_domainBC(LOW,KAXIS),gr_domainBC(HIGH,KAXIS))

  ! The Paramesh call above may have already resulted in some block refining,
  ! so try get the current max level from the lrefine array. This is only used for
  ! diagnostic output. Note that this assumes that lrefine on the current
  ! processor is representative of the grid as a whole.
  cur_treedepth = maxval(lrefine)

  gridDataStruct = CENTER
#if NFACE_VARS > 0
  gridDataStruct = CENTER_FACES
#endif

  lrefineMinSave = lrefine_min
  lrefine_min = min(gr_lrefineMinInit,lrefine_max)

  grid_changed_anytime = grid_changed ! save value established by previous Paramesh initialization

  retainParticles=.false.

  do ntimes = 1, lrefine_max+2
     if (ntimes .EQ. gr_lrefineMinInit) then
        lrefine_min = lrefineMinSave
     end if
     write (block_buff(1,1), '(a)') 'iteration'
     write (int_to_str, '(i7,a1)') ntimes, ','
     write (block_buff(1,2), '(a,1x,a)') trim(adjustl(int_to_str))

     write (block_buff(2,1), '(a)') 'create level'
     write (int_to_str, '(i7)') min(cur_treedepth+1,lrefine_max)
     write (block_buff(2,2), '(a)') trim(adjustl(int_to_str))

     call Logfile_stamp( block_buff, 2, 2, '[GRID gr_expandDomain]')

     call gr_updateData()
     call Grid_getLocalNumBlks(lnblocks)
     whichBlocks = LEAF
     ! Paramesh may have already refined the original root block(s) once by this point.
     ! (That can happen in particular when lrefine_min > 1.)
     ! So initialize all existing blocks in the first iteration here.  This makes sure
     ! that root blocks are not left with uninitialized contents. (Normally that
     ! situation would only last until Grid_fillGuardCells is called with LEAF blocks
     ! at refinement level 2 anyway, since PARENT blocks are then updated as a side
     ! effect by restriction.) - KW
     if (ntimes == 1) whichBlocks = ALL_BLKS
     call Grid_getListOfBlocks(whichBlocks, blkList,count)


#ifndef FLASH_GRID_PARAMESH2
     if (no_permanent_guardcells) then
        call gr_commSetup(gridDataStruct)
     end if
#endif

     do i = 1, count
        !  We need to zero data in case we reuse blocks from previous levels
        !  but don't initialize all data in Simulation_initBlock... in particular
        !  the total vs. internal energies can cause problems in the eos call that
        !  follows.
        call Grid_getBlkPtr(blkList(i), solnData)
        solnData = 0.0
        !      Now reinitialize the solution on the new grid so that it has
        !      the exact solution.
        call gr_fillMetaData (blkList(i), blockDesc)
        call Simulation_initBlock (solnData, blockDesc)
        call Grid_releaseBlkPtr(blkList(i), solnData)
     end do

#ifdef ERAD_VAR
     ! Sum radiation energy density over all meshes. This call is
     ! needed for mesh replication.

     ! DEVNOTE: The interface of this function has not yet been changed
     ! It should create its own iterator.  How to get count?
     call RadTrans_sumEnergy(ERAD_VAR, count, whichBlocks)
#endif

     ! This is here for safety, in case the user did not take care to make things
     ! thermodynamically consistent in the initial state.- KW
     call Timers_start("eos")
     do i = 1, count
        call Eos_wrapped_blkid(gr_eosModeInit,blkLimits,blkList(i))
     end do

     call Timers_stop("eos")

     if(gr_refineOnParticleCount ) then

        !!   This loop initializes the particle positions if
        !!   their count is one of the refinement criteria.
        !!   If the initialization routine intends to keep
        !!   the already initialized particles around, instead
        !!   of reinitializing them as the grid is refined,
        !!   it should return retainParticles true.
        !!   In case of initializing particles from a file,
        !!   if the whole file has been read in
        !!   then particlesPosnsDone should be true, otherwise false.
        !!
        if(.not.retainParticles) then
           particlesPosnsDone=.false.
        end if
        call Particles_initPositions(particlesPosnsDone,retainParticles)
#ifdef DEBUG_PARTICLES
        if (gr_meshMe == MASTER_PE .OR. gr_meshNumProcs .LE. 4) then
           print*,'gr_expandDomain after Particles_initPositions on',gr_meshMe,':',particlesPosnsDone,retainParticles
        end if
#endif
     end if
     dr_simGeneration = dr_simGeneration + 1
     if (ntimes .le. lrefine_max+1) then
        ! Guard cell filling and Eos_wrapped are done in Grid_markRefineDerefine as needed.
        call gr_amrextBuildMultiFabsFromF4Grid(CENTER,lrefine_max,ACTIVE_BLKS)
        call Grid_copyF4DataToMultiFabs(CENTER, nodetype=ACTIVE_BLKS)
        call Grid_markRefineDerefine()
        grid_changed_anytime = max(grid_changed, grid_changed_anytime)
        grid_changed = 0              ! will be 1 after amr_refine_derefine if the grid actually changed
        call amr_refine_derefine()
#ifndef FLASH_GRID_PARAMESH2
        if (grid_changed .NE. 0) mpi_pattern_id = -abs(mpi_pattern_id) !make it different from recognized values
#endif
        if(gr_refineOnParticleCount.and.retainParticles) call Particles_updateRefinement(lnblocks)
        cur_treedepth = max(maxval(lrefine),min(cur_treedepth+1,lrefine_max))

        gr_gcellsUpToDate = .false.
     end if

  end do !ntimes

  grid_changed = max(grid_changed, grid_changed_anytime) !leave global flag true if grid changed in ANY iteration

  if(gr_refineOnParticleCount) then
     if(.not.particlesPosnsDone) call Driver_abortFlash(&
       "This distribution of particles will not fit on the grid. Increase pt_maxPerProc, or decrease the particle count.")
     particlesInitialized=.true.
  end if

  lrefine_min = lrefineMinSave

  call gr_ensureValidNeighborInfo(10)

  call gr_amrextBuildMultiFabsFromF4Grid(CENTER,lrefine_max,ACTIVE_BLKS)
  call Grid_copyF4DataToMultiFabs(CENTER, nodetype=ACTIVE_BLKS)

end subroutine gr_expandDomain
