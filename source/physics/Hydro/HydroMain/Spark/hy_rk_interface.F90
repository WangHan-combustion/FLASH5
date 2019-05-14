!!****if* source/physics/Hydro/HydroMain/Spark/hy_rk_interface
!!
!! NAME
!!   hy_rk_interface
!!
!! SYNOPSIS
!!   use hy_rk_interface : ONLY
!!
!!  DESCRIPTION
!!    Interface for internal Spark Hydro subroutines
!!
!!***
module hy_rk_interface

#include "Flash.h"
#include "constants.h"

  interface
     subroutine hy_rk_getFaceFlux (blockID,limits)
       implicit none
       integer, intent(in) :: blockID
       integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
     end subroutine hy_rk_getFaceFlux
  end interface

  interface
     subroutine hy_rk_updateSoln (blockID, dt, dtOld, limits, coeffs)
       implicit none
       integer, intent(IN) :: blockID
       integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
       real, intent(IN) :: dt, dtOld
       real, dimension(3), intent(IN) :: coeffs
     end subroutine hy_rk_updateSoln
  end interface

  interface
     subroutine hy_rk_eos(limits)
       implicit none
       integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
     end subroutine hy_rk_eos
  end interface

  interface
     subroutine hy_rk_getGravAccel(blockID,limits)
       implicit none
       integer, intent(IN) :: blockID
       integer, intent(IN) :: limits(LOW:HIGH,MDIM)
     end subroutine hy_rk_getGravAccel
  end interface

end module hy_rk_interface
