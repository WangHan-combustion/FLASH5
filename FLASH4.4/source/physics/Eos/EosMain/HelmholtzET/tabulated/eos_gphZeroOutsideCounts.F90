!!****if* source/physics/Eos/EosMain/Tabulated/eos_gphZeroOutsideCounts
!!
!! NAME
!!
!!  eos_gphZeroOutsideCounts
!!
!! SYNOPSIS
!!
!!  call eos_gphZeroOutsideCounts()
!!
!! DESCRIPTION
!! 
!! 
!! ARGUMENTS
!!
!!   No arguments
!!
!!
!!
!!***

subroutine eos_gphZeroOutsideCounts()

  use eos_gphData,                ONLY : eos_gphAllDiag

  implicit none

#include "Flash.h"

  integer :: species

  do species = 1,NSPECIES
     eos_gphAllDiag(species)%highTempCount = 0
     eos_gphAllDiag(species)%highDensCount = 0
     eos_gphAllDiag(species)%highestTemp = -999.0
     eos_gphAllDiag(species)%highestDens = -999.0
     eos_gphAllDiag(species)%highTempVarsLookedUp(:) = .FALSE.
     eos_gphAllDiag(species)%highDensVarsLookedUp(:) = .FALSE.
  end do
  

end subroutine eos_gphZeroOutsideCounts
