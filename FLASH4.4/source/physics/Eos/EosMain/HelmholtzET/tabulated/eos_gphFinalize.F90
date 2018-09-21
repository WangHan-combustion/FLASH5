!!****if* source/physics/Eos/EosMain/Tabulated/eos_gphFinalize
!!
!! NAME
!!
!!  eos_gphFinalize
!!
!! SYNOPSIS
!!
!!  call eos_gphFinalize ()
!!
!! DESCRIPTION
!!
!!  Clean up the Opacity unit.
!!
!! ARGUMENTS
!!
!!***
#include "Flash.h"

subroutine eos_gphFinalize ()

  use eos_gphData,   ONLY : eos_gphleKind,              &
                             eos_gphleName,              &
                             eos_groupName,              &
                             eos_gphIonizationKind,         &
                             eos_gphIntEnergyKind,           &
                             eos_gphHeatCpKind,          &
           eos_allTab, &
           eos_gphAllDiag, &
                             EOS_TAB_NALLTAB
  implicit none

  integer :: species, i, j

  do species = 1,NSPECIES
     do i = 1,EOS_TAB_NALLTAB
        ! Currently, the only way these pointers get different from null is by being allocated. - KW
        if (associated(eos_allTab(species)%tg(i)%table)) then
           do j = LBOUND(eos_allTab(species)%tg(i)%table,1), &
                UBOUND(eos_allTab(species)%tg(i)%table,1)
              if (associated(eos_allTab(species)%tg(i)%table(j)%table)) deallocate(eos_allTab(species)%tg(i)%table(j)%table)
           end do
           deallocate(eos_allTab(species)%tg(i)%table)
        end if
        if (associated(eos_allTab(species)%tg(i)%mgTable)) deallocate(eos_allTab(species)%tg(i)%mgTable)
        if (associated(eos_allTab(species)%tg(i)%td%Temperatures)) deallocate(eos_allTab(species)%tg(i)%td%Temperatures)
        if (associated(eos_allTab(species)%tg(i)%td%Densities)) deallocate(eos_allTab(species)%tg(i)%td%Densities)
     end do
  end do
  deallocate(eos_allTab)
  deallocate(eos_gphAllDiag)

  deallocate (eos_gphleKind)
  deallocate (eos_gphleName)
  deallocate (eos_gphIonizationKind)
  deallocate (eos_gphIntEnergyKind)
  deallocate (eos_gphHeatCpKind)

end subroutine eos_gphFinalize
