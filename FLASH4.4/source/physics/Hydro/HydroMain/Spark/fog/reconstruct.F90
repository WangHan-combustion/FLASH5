!!****if* source/physics/Hydro/HydroMain/Spark/reconstruct
!!
!!  NAME
!!
!!  hy_rk_reconstruct
!!
!!  SYNOPSIS
!!
!!  call hy_rk_reconstruct (  )
!!
!!  DESCRIPTION
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine reconstruct(uPlus, uMinus, data1d, flat, size1d, ind, dx)
  use Hydro_data, ONLY : hy_limRad
  use Timers_interface, ONLY : Timers_stop, Timers_start

  implicit none

#include "Flash.h"
#include "Spark.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS

  integer, intent(IN) :: size1d, ind
  real, intent(IN) :: data1d(NRECON,size1d), flat, dx
  real, dimension(NRECON), intent(OUT) :: uPlus, uMinus

  uPlus  = data1d(:,ind)
  uMinus = data1d(:,ind)

end subroutine reconstruct
