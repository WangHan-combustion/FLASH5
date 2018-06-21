!!****if* source/Simulation/SimulationMain/StreamingSineWave/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!   Initialize general solution data for streaming sine wave radiation
!!   problem (steady-state)
!!
!! ARGUMENTS
!!
!!
!!
!! PARAMETERS
!!
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use Eos_interface, ONLY : Eos
  use Logfile_interface, ONLY : Logfile_stamp
  use ProgramHeaderModule, ONLY : nE, nDOF
  use RadiationFieldsModule, ONLY : nCR, nSpecies
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use ThornadoInitializationModule, ONLY : InitThornado

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Multispecies.h"

  real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction
  real, dimension(EOS_NUM) :: eosData

  call Driver_getMype(MESH_COMM, sim_meshMe)
  call Driver_getMype(GLOBAL_COMM, sim_globalMe)

  call RuntimeParameters_get('dens_i', sim_dens_i)
  call RuntimeParameters_get('temp_i', sim_temp_i)
  sim_xn_i(:) = 0.0e0
  sim_xn_i(SPECIES_BEGIN) = 1.0e0

  massFraction(:) = sim_xn_i(:)

  eosData(EOS_DENS) = sim_dens_i
  eosData(EOS_TEMP) = sim_temp_i

  call Eos(MODE_DENS_TEMP,1,eosData,massFraction)

  sim_velx_i = 0.0e0
  sim_vely_i = 0.0e0
  sim_velz_i = 0.0e0
  sim_pres_i = eosData(EOS_PRES)
  sim_eint_i = eosData(EOS_EINT)
  sim_etot_i = sim_eint_i + 0.5*(sim_velx_i**2 + sim_vely_i**2 + sim_velz_i**2)
  sim_gamc_i = eosData(EOS_GAMC)
  sim_game_i = sim_pres_i/(sim_eint_i*sim_dens_i) + 1.0e0

  call InitThornado( NDIM, THORNADO_NE, THORNADO_NSPECIES )

  sim_nComp = nSpecies * nCR * nE * nDOF

  return

end subroutine Simulation_init
