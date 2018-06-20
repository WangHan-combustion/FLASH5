!!****if* source/physics/RadTrans/RadTransTwoMoment/RadTrans_init
!!
!!  NAME 
!!
!!  RadTrans_init
!!
!!  SYNOPSIS
!!
!!  call RadTrans_init()
!!
!!  DESCRIPTION 
!!    Initialize radiative transfer unit
!!
!! ARGUMENTS
!!
!!
!!***

#include "constants.h"

subroutine RadTrans_init()
  use RadTrans_data
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype, Driver_getComm
  implicit none

  call Driver_getMype(MESH_COMM,rt_meshMe)
  call RuntimeParameters_get ("useRadTrans", rt_useRadTrans)
  
  return

end subroutine RadTrans_init
