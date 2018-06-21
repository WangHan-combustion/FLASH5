!!****f* source/physics/RadTrans/RadTransTwoMoment/RadTrans
!!
!!  NAME 
!!
!!  RadTrans
!!
!!  SYNOPSIS
!!
!!  call RadTrans( real(IN)    :: dt, 
!!       optional, integer(IN) :: pass)
!!
!!  DESCRIPTION 
!!      This subroutine performs the radiatiative transfer calculation
!!      for this step. 
!!
!! ARGUMENTS
!!
!!   dt     : The time step
!!   pass   : Reverses direction of solve
!!
!!***

!!REORDER(4): solnData

subroutine RadTrans_desc( dt, pass )

  use Driver_interface, ONLY : Driver_abortFlash
  use FluidFieldsModule, ONLY : uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
  use Grid_interface, ONLY : Grid_fillGuardCells, &
     Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getBlkPtr, &
     Grid_releaseBlkPtr, Grid_getMaxRefinement, Grid_getLeafIterator, &
     Grid_releaseLeafIterator, Grid_getBlkBoundBox
  use ProgramHeaderModule, ONLY : nE, nDOF, nDOFX
  use RadiationFieldsModule, ONLY : nSpecies, uCR, nCR
  use RadTrans_data, ONLY : rt_useRadTrans
  use ThornadoInitializationModule, ONLY : InitThornado_Patch, FreeThornado_Patch
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use TimeSteppingModule_Castro, ONLY : Update_IMEX_PC2
  use UnitsModule, ONLY : Gram, Centimeter, Second

  use leaf_iterator, ONLY : leaf_iterator_t
  use block_metadata, ONLY : block_metadata_t

  implicit none

#include "Flash.h"
#include "constants.h"

  real,    intent(in) :: dt
  integer, intent(in), optional :: pass

  integer :: blockID, thisBlock, blockCount
  real, pointer, dimension(:,:,:,:) :: solnData
  real, dimension(LOW:HIGH,MDIM) :: boundBox
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  logical, parameter :: getGuardCells = .true.
  integer :: iSize, jSize, kSize, iSizeGC, jSizeGC, kSizeGC

  integer :: i, j, k, ic, jc, kc, is, im, ie, id, ii

  integer :: level, maxLev
  type(leaf_iterator_t) :: itor
  type(block_metadata_t) :: blockDesc

  integer, parameter :: my_ngrow = 2
  integer :: swE
  real :: eL, eR
  integer :: nX(MDIM), swX(MDIM)
  real :: xL(MDIM), xR(MDIM)

  ! Old and new radiation state
  real, allocatable, dimension(:,:,:,:) :: U_R_o, U_R_n

  real, parameter :: conv_dens = Gram / Centimeter**3
  real, parameter :: conv_mom  = Gram / Centimeter**2 / Second
  real, parameter :: conv_enr  = Gram / Centimeter / Second**2
  real, parameter :: conv_ne   = 1.0 / Centimeter**3
  real, parameter :: conv_J    = Gram/Second**2/Centimeter
  real, parameter :: conv_H    = Gram/Second**3

  if (.NOT. rt_useRadTrans) return

  if ( my_ngrow /= 2 ) then
     call Driver_abortFlash("Need two ghost cells in call_to_thornado!")
  end if

  call Timers_start("RadTrans")

  swX(1:3) = my_ngrow
  swE = 0
  eL = 0.0
  eR = 0.0

#ifdef FLASH_GRID_UG
  maxLev = 1
#else
  call Grid_getMaxRefinement(maxLev,mode=1) !mode=1 means lrefine_max, which does not change during sim.
#endif

  do level = 1, maxLev
     thisBlock = 0
     call Grid_getLeafIterator(itor, level=level)
     do while(itor%is_valid())
        call itor%blkMetaData(blockDesc)
        thisBlock = thisBlock + 1

        ! get dimensions/limits and coordinates
        blkLimitsGC = blockDesc%limitsGC
        iSizeGC = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
        jSizeGC = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
        kSizeGC = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

        blkLimits = blockDesc%limits
        iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
        jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
        kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1
        nX(1:3) = [ iSize, jSize, kSize ]

        call Grid_getBlkBoundBox(blockDesc, boundBox)

        ! Convert cm to m for Thornado
        xL(1:NDIM) = boundBox(LOW,1:NDIM) * 1.0e-2
        xR(1:NDIM) = boundBox(HIGH,1:NDIM) * 1.0e-2

        ! Get a pointer to solution data
        call Grid_getBlkPtr(blockDesc, solnData)

        call InitThornado_Patch(nX, swX, xL, xR, swE, eL, eR)

        ! Copy from the Flash arrays into Thornado arrays from InitThornado_Patch
        do kc = blkLimits(LOW,KAXIS)-swX(3), blkLimits(HIGH,KAXIS)+swX(3)
           do jc = blkLimits(LOW,JAXIS)-swX(2), blkLimits(HIGH,JAXIS)+swX(2)
              do ic = blkLimits(LOW,IAXIS)-swX(1), blkLimits(HIGH,IAXIS)+swX(1)

                 ! U_R_o spatial indices start at lo - (number of ghost zones)
                 ! uCR spatial indices start at 1 - (number of ghost zones)
                 i = ic - blkLimits(LOW,IAXIS) + 1
                 j = jc - blkLimits(LOW,JAXIS) + 1
                 k = kc - blkLimits(LOW,KAXIS) + 1

                 ! Thornado uses units where c = G = k = 1, Meter = 1
                 uCF(1:nDOFX,i,j,k,iCF_D)  = solnData(DENS_VAR,ic,jc,kc) * conv_dens
                 uCF(1:nDOFX,i,j,k,iCF_S1) = solnData(DENS_VAR,ic,jc,kc) * solnData(VELX_VAR,  ic,jc,kc) * conv_mom
                 uCF(1:nDOFX,i,j,k,iCF_S2) = solnData(DENS_VAR,ic,jc,kc) * solnData(VELY_VAR,  ic,jc,kc) * conv_mom
                 uCF(1:nDOFX,i,j,k,iCF_S3) = solnData(DENS_VAR,ic,jc,kc) * solnData(VELZ_VAR,  ic,jc,kc) * conv_mom
                 uCF(1:nDOFX,i,j,k,iCF_E)  = solnData(DENS_VAR,ic,jc,kc) * solnData(ENER_VAR,  ic,jc,kc) * conv_enr
                 uCF(1:nDOFX,i,j,k,iCF_Ne) = solnData(DENS_VAR,ic,jc,kc) * solnData(YE_MSCALAR,ic,jc,kc) * conv_ne

                 do is = 1, nSpecies ; do im = 1, nCR ; do ie = 1, nE ; do id = 1, nDOF
                    ii = THORNADO_BEGIN + (is-1)*(nCR*nE*nDOF) + (im-1)*(nE*nDOF) + (ie-1)*nDOF + (id-1)

                    if ( im == 1 ) then
                       uCR(id,ie,i,j,k,im,is) = solnData(ii,ic,jc,kc) * conv_J
                    else if ( im > 1 ) then
                       uCR(id,ie,i,j,k,im,is) = solnData(ii,ic,jc,kc) * conv_H
                    end if

                 end do ; end do ; end do ; end do

              end do
           end do
        end do

        ! Call the Fortran interface that lives in the Thornado repo
        call Update_IMEX_PC2(dt*Second, uCF, uCR)

        ! Copy back from the thornado arrays into Castro arrays
        do kc = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do jc = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do ic = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 ! uCR spatial indices start at 1 - ng
                 ! U_R_n spatial indices start at lo
                 i = ic - blkLimits(LOW,IAXIS) + 1
                 j = jc - blkLimits(LOW,JAXIS) + 1
                 k = kc - blkLimits(LOW,KAXIS) + 1

                 !solnData(DENS_VAR,  ic,jc,kc) = uCF(1,i,j,k,iCF_D)  / conv_dens
                 !solnData(VELX_VAR,  ic,jc,kc) = uCF(1,i,j,k,iCF_S1) / uCF(1,i,j,k,iCF_D) / conv_mom
                 !solnData(VELY_VAR,  ic,jc,kc) = uCF(1,i,j,k,iCF_S2) / uCF(1,i,j,k,iCF_D) / conv_mom
                 !solnData(VELZ_VAR,  ic,jc,kc) = uCF(1,i,j,k,iCF_S3) / uCF(1,i,j,k,iCF_D) / conv_mom
                 !solnData(ENER_VAR,  ic,jc,kc) = uCF(1,i,j,k,iCF_E)  / uCF(1,i,j,k,iCF_D) / conv_enr
                 !solnData(YE_MSCALAR,ic,jc,kc) = uCF(1,i,j,k,iCF_Ne) / uCF(1,i,j,k,iCF_D) / conv_ne

                 do is = 1, nSpecies ; do im = 1, nCR ; do ie = 1, nE ; do id = 1, nDOF
                    ii = THORNADO_BEGIN + (is-1)*(nCR*nE*nDOF) + (im-1)*(nE*nDOF) + (ie-1)*nDOF + (id-1)

                    if ( im == 1 ) then
                       solnData(ii,ic,jc,kc) = uCR(id,ie,i,j,k,im,is) / conv_J
                    else if ( im > 1 ) then
                       solnData(ii,ic,jc,kc) = uCR(id,ie,i,j,k,im,is) / conv_H
                    end if

                 end do ; end do ; end do ; end do

              end do
           end do
        end do

        call FreeThornado_Patch()

        call Grid_releaseBlkPtr(blockDesc,solnData)
        nullify(solnData)

        call itor%next()

     end do
     call Grid_releaseLeafIterator(itor)
  end do

  call Timers_stop("RadTrans")

  return

end subroutine RadTrans_desc
