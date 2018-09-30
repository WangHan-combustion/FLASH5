!  This function is a callback for the Eos/Helmholtz/ExternelAbarZbar
!  EOS implementation and similar Eos implementations.
!  The local Abar and Zbar are calculated from the species mass
!  fractions from the provided solution state.

!!REORDER(2): solnScalars

subroutine Multispecies_computeAbarZbar(solnScalars, abarData, zbarData)

  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
       Multispecies_getSumFrac

  implicit none

#include "Flash.h"
#include "Multispecies.h"

  real, intent(in), dimension(SPECIES_BEGIN:,:)  :: solnScalars
  real, intent(inout), dimension(:)                      :: abarData, zbarData
  real :: abarInv,zbarFrac

!  integer,parameter :: sumyi = SUMY_MSCALAR-SPECIES_BEGIN+1
!  integer,parameter :: yei = YE_MSCALAR-SPECIES_BEGIN+1

  integer :: i, ubo

#ifdef INDEXREORDER
  ubo = ubound(solnScalars,1)
#else
  ubo = ubound(solnScalars,2)
#endif

  do i = 1, ubo

     call Multispecies_getSumInv(A, abarInv, solnScalars(SPECIES_BEGIN:SPECIES_END,i))
     abarData(i) = 1.e0 / abarInv

     call Multispecies_getSumFrac(Z, zbarFrac, solnScalars(SPECIES_BEGIN:SPECIES_END,i))
     zbarData(i) = abarData(i) * zbarFrac

     ! check
!     if ( .not. (abarData(i) > 0.0 ) ) then
!        !print but don't die, EOS will die and report density and temperature
!        write (6,*) "[Multispecies_computeAbarZbar] negative abar =", abarData(i)
!        write (6,*) "grid data sumy, ye :"
!        write (6,*) solnScalars(sumyi,i), solnScalars(yei,i)
!     endif
  enddo

end subroutine Multispecies_computeAbarZbar
