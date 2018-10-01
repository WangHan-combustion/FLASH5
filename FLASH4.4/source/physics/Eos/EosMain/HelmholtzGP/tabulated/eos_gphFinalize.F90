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

  use eos_gphData,   ONLY : eos_gphKind,              &
                             eos_gphFileName,              &
                             eos_gphIonizationKind,         &
                             eos_gphIntEnergyKind,           &
                             eos_gphHeatCpKind,          &
                             TheGphTable => eos_gphTheTable, &
           eos_gphAllDiag, &
                             EOS_GPH_NALLTAB
  implicit none

  integer :: i, j

     do i = 1,EOS_GPH_NALLTAB
        ! Currently, the only way these pointers get different from null is by being allocated. - KW
        nullify(theGphTable % tg(1) % table(1) % table)
        nullify(theGphTable % tg(1) % table(1) % ells)
           do j = LBOUND(theGphTable%tg(i)%table,1), &
                UBOUND(theGphTable%tg(i)%table,1)
              if (associated(theGphTable%tg(i)%table(j)%table)) deallocate(theGphTable%tg(i)%table(j)%table)
              if (associated(theGphTable%tg(i)%table(j)%ells )) deallocate(theGphTable%tg(i)%table(j)%ells)
              if (associated(theGphTable%tg(i)%table(j)%derDefs)) deallocate(theGphTable%tg(i)%table(j)%derDefs)
           end do
!!$        if (associated(theGphTable%tg(i)%td%Temperatures)) deallocate(theGphTable%tg(i)%td%Temperatures)
!!$        if (associated(theGphTable%tg(i)%td%Densities)) deallocate(theGphTable%tg(i)%td%Densities)
     end do

  deallocate(theGphTable)
!!$  deallocate(eos_gphAllDiag)

  deallocate (eos_gphKind)
  deallocate (eos_gphFileName)
  deallocate (eos_gphIonizationKind)
  deallocate (eos_gphIntEnergyKind)
  deallocate (eos_gphHeatCpKind)

end subroutine eos_gphFinalize
