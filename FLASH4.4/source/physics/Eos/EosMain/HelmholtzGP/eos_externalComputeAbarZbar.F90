!!****if* source/physics/Eos/EosMain/Helmholtz/ExternalAbarZbar/eos_externalComputeAbarZbar
!!
!! NAME
!!
!!  eos_externalComputeAbarZbar
!!
!! SYNOPSIS
!!
!!  call eos_externalComputeAbarZbar(real(in), dimension(:,:)  :: solnscalars,
!!                                   real(out), dimension(:)  :: abardata,
!!                                   real(out), dimension(:)  :: zbardata)
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!! Klaus Weide   2013
!!
!!  This routine private to the Eos unit serves to implement callbacks
!!  for the Eos/Helmholtz/ExternelAbarZbar EOS implementation.
!!  Code units that implement ways for computing Abar and Zbar from
!!  the solutions state (or otherwise) are polled here by calling their
!!  public UNITNAME_computeAbarZbar interfaces.
!!
!!  It is assumed that not more than one polled units have implementations
!!  included in the simulation that actually provide Abar and Zbar.
!!  If several units do provide the data, the last one polled here will win.
!!
!! ARGUMENTS
!!
!!   solnscalars : scalars of the solution 
!!
!!   abardata : abar info
!!
!!   zbardata : zbar info
!!
!!
!!
!!***

#include "Flash.h"
subroutine eos_externalComputeAbarZbar(solnScalars, abarData, zbarData)

  use Burn_interface,  ONLY: Burn_computeAbarZbar
  use Multispecies_interface,  ONLY: Multispecies_computeAbarZbar
  use Eos_data, ONLY: eos_singleSpeciesA, eos_singleSpeciesZ

  implicit none

  real, intent(in),  dimension(:,:)  :: solnScalars
  real, intent(out), dimension(:)    :: abarData, zbarData

#if defined(FLASH_SOURCEBURN)
  call Burn_computeAbarZbar(solnScalars, abarData, zbarData)
#elif defined(FLASH_MULTISPECIES) && (NSPECIES > 0)
  call Multispecies_computeAbarZbar(solnScalars, abarData, zbarData)
#else
  abarData(:) = eos_singleSpeciesA
  zbarData(:) = eos_singleSpeciesZ
  print*,'abarData:',abarData
  print*,'zbarData:',zbarData
#endif
     if ( .not. (zbarData(1) > 0.0 ) ) then
        !print but don't die, EOS will die and report density and temperature
        write (6,*) "[eos_externalComputeAbarZbar] nonpositive zbar =", zbarData(1)
!!$        write (6,*) "grid data sumy, ye :"
!!$        write (6,*) solnScalars(sumyi,i), solnScalars(yei,i)
     endif
     if ( .not. (abarData(1) > 0.0 ) ) then
        !print but don't die, EOS will die and report density and temperature
        write (6,*) "[eos_externalComputeAbarZbar] nonpositive abar =", abarData(1)
!!$        write (6,*) "grid data sumy, ye :"
!!$        write (6,*) solnScalars(sumyi,i), solnScalars(yei,i)
     endif

end subroutine eos_externalComputeAbarZbar
