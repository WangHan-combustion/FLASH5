!!****if* source/physics/RadTrans/RadTransTwoMoment/RadTrans_data
!!
!!  NAME 
!!
!!  RadTrans_data
!!
!!  SYNOPSIS
!!   use RadTrans_data
!!
!!  DESCRIPTION 
!!    Stores local data for the RadTrans unit
!!
!!***
module RadTrans_data
  implicit none
  
  integer, save :: rt_meshMe ! Process rank
  logical, save :: rt_useRadTrans
  
end module RadTrans_data
