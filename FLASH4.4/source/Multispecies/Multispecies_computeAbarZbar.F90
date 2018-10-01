!  This function is a callback for the Eos/Helmholtz/ExternelAbarZbar
!  EOS implementation and similar Eos implementations.
!  The local Abar and Zbar are calculated from the species mass
!  fractions from the provided solution state.

!!REORDER(2): solnScalars

subroutine Multispecies_computeAbarZbar(solnScalars, abarData, zbarData)

  implicit none

#include "Flash.h"

  real, intent(in), dimension(SPECIES_BEGIN:,:)  :: solnScalars
  real, intent(inout), dimension(:) :: abarData, zbarData

end subroutine Multispecies_computeAbarZbar
