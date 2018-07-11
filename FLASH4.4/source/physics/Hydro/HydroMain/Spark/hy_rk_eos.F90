!!****if* source/physics/Hydro/HydroMain/Spark/hy_rk_eos
!!
!!  NAME
!!
!!  hy_rk_eos
!!
!!  SYNOPSIS
!!
!!  call hy_rk_eos ( integer(IN) :: blockID )
!!
!!  DESCRIPTION
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine hy_rk_eos(limits)

  use Hydro_data, ONLY : hy_starState, hy_threadWithinBlock
  use Eos_interface, ONLY : Eos_putData, Eos_getData, Eos

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Spark.h"
#include "Eos.h"

  integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
  integer :: i,j,k,vecLen
  integer,dimension(MDIM) :: pos

  real, dimension(NSPECIES*MAXCELLS) :: massFraction
  real, dimension(EOS_NUM*MAXCELLS) :: eosData
  real, pointer :: tempData(:,:,:,:)

  tempData   => hy_starState

  !$omp parallel if (hy_threadWithinBlock .AND. NDIM > 1) &
  !$omp default(none) &
  !$omp firstprivate(vecLen) &
  !$omp private(i,j,k,pos,eosData,massFraction)&
  !$omp shared(limits,hy_starState,tempData)

  pos(IAXIS) = limits(LOW,IAXIS)
  vecLen     = limits(HIGH,IAXIS)-limits(LOW,IAXIS)+1

  !  Begin loop over zones
  !$omp do schedule(static) collapse(2)
  do k = limits(LOW,KAXIS), limits(HIGH,KAXIS)
     do j = limits(LOW,JAXIS), limits(HIGH,JAXIS)

       pos(JAXIS) = j
       pos(KAXIS) = k
       call Eos_getData(IAXIS,pos,vecLen,tempData,CENTER,eosData,massFraction)
       call Eos(MODE_DENS_EI,vecLen,eosData,massFraction)
       call Eos_putData(IAXIS,pos,vecLen,tempData,CENTER,eosData)

     end do
  end do
  !$omp end do nowait
  !$omp end parallel

  nullify(tempData)

end subroutine hy_rk_eos
