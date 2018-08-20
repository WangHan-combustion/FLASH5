!!****if* source/physics/Hydro/localAPI/hy_getRiemannState_loop1
!!
!! NAME
!!
!!  hy_getRiemannState_loop1
!!
!! SYNOPSIS
!!
!!  hy_getRiemannState_loop1( integer(IN) :: dir,
!!                          real, pointer, dimension (:,:,:,:) :: U
!!                          integer(IN) :: i,j,k
!!                          real(INOUT), dimension(5) :: tempState)
!!
!! DESCRIPTION
!!
!!  This routine is called from a nested loop from hy_getRiemannState
!!  and does computations for flattening tilde.
!!
!! ARGUMENTS
!!
!!  dir     - firection x, y or z
!!  U   - Pointer to blockdata for the block being computed
!!  i,j,k - indeces of the cell being computed from the loop for this block
!!  tempState - array containing 5 outputs to the calling function
!!                    - dp1, dp2, dv1, presL, presR
!!   
!! NOTE : :
!! Will fail for NDIM <3 in current format because of missing #if NDIM directive
!!***

!!REORDER(4):U

subroutine  hy_getReimannState_loop1(dir, U, i,j,k, temp)

#include "UHD.h"
implicit none
    integer, intent(in) :: dir
    real, pointer, dimension(:,:,:,:) :: U
    integer, intent(in) :: i,j,k
    real, dimension(5), intent(inout) :: temp
    
return
end subroutine hy_getReimannState_loop1
