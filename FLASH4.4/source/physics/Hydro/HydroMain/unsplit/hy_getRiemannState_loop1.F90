!!****if* source/physics/Hydro/HydroMain/unsplit/hy_getRiemannState_loop1
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
 	   
 	   integer,dimension(5) :: ii,jj,kk
	   integer :: var
	   ii=i
	   jj=j
	   kk=k
	   select case(dir)
	   case(DIR_X)
             ii(1)=i+1
             ii(2)=i-1
             ii(3)=i+2
             ii(4)=i-2
	     var=VELX_VAR
	   case(DIR_Y)
             jj(1)=j+1
             jj(2)=j-1
             jj(3)=j+2
             jj(4)=j-2
	     var=VELY_VAR
	   case(DIR_Z)
             kk(1)=k+1
             kk(2)=k-1
             kk(3)=k+2
             kk(4)=k-2
	     var=VELY_VAR
          end select
        temp(1)=U(ii(1),jj(1),kk(1), PRES_VAR)-U(ii(2),jj(2),kk(2),PRES_VAR)
        temp(2)=U(ii(3),jj(3),kk(3), PRES_VAR)-U(ii(4),jj(4),kk(4),PRES_VAR)
	temp(3)=U(ii(1),jj(1),kk(1), PRES_VAR)-U(ii(2),jj(2),kk(2), PRES_VAR)
	temp(4)=U(ii(2),jj(2),kk(2), PRES_VAR)
	temp(5)=U(ii(2),jj(2),kk(1), PRES_VAR)

end subroutine hy_getReimannState_loop1