!!****if* source/physics/Hydro/HydroMain/Spark/reconstruct
!!
!!  NAME
!!
!!  reconstruct
!!
!!  SYNOPSIS
!!
!!  call reconstruct (  )
!!
!!  DESCRIPTION
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine reconstruct(uPlus, uMinus, data1d, flat, size1d, ind, dx)

  implicit none

#include "Flash.h"
#include "Spark.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS

  integer, intent(IN) :: size1d, ind
  real, intent(IN) :: data1d(NRECON,size1d), flat, dx
  real, dimension(NRECON), intent(OUT) :: uPlus, uMinus
  real, dimension(NRECON) :: delbar, del_p, del_m, invSumAlpha
  real, dimension(NRECON,3) :: W5p, W5m, betaWeno, Alpha5, omega
  integer :: g
  real, parameter :: epsilon = 1.e-36
  real, parameter :: n13o12 = 13./12.
  real, dimension(NRECON) :: mags, betaSum
  !! Set WENO5 coefficients once and for all
  !u_{1,i+1/2}= 2/6*u_{i-2} -7/6*u_{i-1} +11/6*u_{i}
  real, dimension(3), parameter :: coeff1p1(1:3) = (/ 2./6., -7./6., 11./6./)
  !u_{2,i+1/2}=-1/6*u_{i-2} +5/6*u_{i-1} + 2/6*u_{i}
  real, dimension(3), parameter :: coeff1p2(1:3) = (/-1./6.,  5./6.,  2./6./)
  !u_{3,i+1/2}= 2/6*u_{i-2} +5/6*u_{i-1} - 1/6*u_{i}
  real, dimension(3), parameter :: coeff1p3(1:3) = (/ 2./6.,  5./6., -1./6./)
  !(gamma1,gamma2,gamma3)
  real, dimension(3), parameter :: coeff2p(1:3)   = (/0.1, 0.6, 0.3/)
  !u_{1,i-1/2}=-1/6*u_{i-2} +5/6*u_{i-1} + 2/6*u_{i}
  real, dimension(3), parameter :: coeff1m1(1:3) = (/-1./6.,  5./6.,  2./6./)
  !u_{2,i-1/2}= 2/6*u_{i-2} +5/6*u_{i-1} - 1/6*u_{i}
  real, dimension(3), parameter :: coeff1m2(1:3) = (/ 2./6.,  5./6., -1./6./)
  !u_{3,i-1/2}=11/6*u_{i-2} -7/6*u_{i-1} + 2/6*u_{i}
  real, dimension(3), parameter :: coeff1m3(1:3) = (/ 11./6.,-7./6.,  2./6./)
  !(gamma1,gamma2,gamma3)
  real, dimension(3), parameter :: coeff2m(1:3)   = (/0.3, 0.6, 0.1/)

  ! Interpolation stencil for weno

  !! Calculate interface values at i+1/2
  W5p(:,1) = coeff1p1(1)*data1d(:,ind-2) + coeff1p1(2)*data1d(:,ind-1) + coeff1p1(3)*data1d(:,ind)
  W5p(:,2) = coeff1p2(1)*data1d(:,ind-1) + coeff1p2(2)*data1d(:,ind)   + coeff1p2(3)*data1d(:,ind+1)
  W5p(:,3) = coeff1p3(1)*data1d(:,ind)   + coeff1p3(2)*data1d(:,ind+1) + coeff1p3(3)*data1d(:,ind+2)

  !! Calculate interface values at i-1/2
  W5m(:,1) = coeff1m1(1)*data1d(:,ind-2) + coeff1m1(2)*data1d(:,ind-1) + coeff1m1(3)*data1d(:,ind)
  W5m(:,2) = coeff1m2(1)*data1d(:,ind-1) + coeff1m2(2)*data1d(:,ind)   + coeff1m2(3)*data1d(:,ind+1)
  W5m(:,3) = coeff1m3(1)*data1d(:,ind)   + coeff1m3(2)*data1d(:,ind+1) + coeff1m3(3)*data1d(:,ind+2)

  !! Calculate smoothness indicators at i+1/2
  betaWeno(:,1) = n13o12*(data1d(:,ind-2) - 2.*data1d(:,ind-1) +    data1d(:,ind)  )**2 &
       +            0.25*(data1d(:,ind-2) - 4.*data1d(:,ind-1) + 3.*data1d(:,ind)  )**2
  betaWeno(:,2) = n13o12*(data1d(:,ind-1) - 2.*data1d(:,ind)   +    data1d(:,ind+1))**2 &
       +            0.25*(data1d(:,ind-1)                      -    data1d(:,ind+1))**2
  betaWeno(:,3) = n13o12*(data1d(:,ind)   - 2.*data1d(:,ind+1) +    data1d(:,ind+2))**2 &
       +            0.25*(3.*data1d(:,ind)- 4.*data1d(:,ind+1) +    data1d(:,ind+2))**2

  !! Use problem-adaptive epsilong as in Tchekovskoy+07, A3
  ! This does not seem to work with the WENO-Z indicators of Borges+08
  ! mags(:) = data1d(:,ind-2)**2 + data1d(:,ind-1)**2 + data1d(:,ind)**2 &
  !      + data1d(:,ind+1)**2 + data1d(:,ind+2)**2
  ! betaWeno(:,1) = betaWeno(:,1) + epsilon*mags(:) + TINY(1.0)
  ! betaWeno(:,2) = betaWeno(:,2) + epsilon*mags(:) + TINY(1.0)
  ! betaWeno(:,3) = betaWeno(:,3) + epsilon*mags(:) + TINY(1.0)

  !! This is WENO-Z: this is very similar to weno5 with wenoExp=1
  Alpha5(:,1) = coeff2p(1)*(1.+(abs(betaWeno(:,1)-betaWeno(:,3))/(betaWeno(:,1)+epsilon)))
  Alpha5(:,2) = coeff2p(2)*(1.+(abs(betaWeno(:,1)-betaWeno(:,3))/(betaWeno(:,2)+epsilon)))
  Alpha5(:,3) = coeff2p(3)*(1.+(abs(betaWeno(:,1)-betaWeno(:,3))/(betaWeno(:,3)+epsilon)))

  !! Normalize weights at i+1/2
  invSumAlpha = 1./(Alpha5(:,1)+Alpha5(:,2)+Alpha5(:,3))
  omega(:,1)  = Alpha5(:,1)*invSumAlpha
  omega(:,2)  = Alpha5(:,2)*invSumAlpha
  omega(:,3)  = Alpha5(:,3)*invSumAlpha

  !! Compute interface value at i+1/2
  uPlus = omega(:,1)*W5p(:,1) + omega(:,2)*W5p(:,2) + omega(:,3)*W5p(:,3)
  !! Apply flattening
  uPlus = flat*uPlus + (1.-flat)*data1d(:,ind)

  !! Now move on to i-1/2
  !! This is WENO-Z
  Alpha5(:,1) = coeff2m(1)*(1.+(abs(betaWeno(:,1)-betaWeno(:,3))/(betaWeno(:,1)+epsilon)))
  Alpha5(:,2) = coeff2m(2)*(1.+(abs(betaWeno(:,1)-betaWeno(:,3))/(betaWeno(:,2)+epsilon)))
  Alpha5(:,3) = coeff2m(3)*(1.+(abs(betaWeno(:,1)-betaWeno(:,3))/(betaWeno(:,3)+epsilon)))

  !! Normalize weights at i-1/2
  invSumAlpha = 1./(Alpha5(:,1)+Alpha5(:,2)+Alpha5(:,3))
  omega(:,1)  = Alpha5(:,1)*invSumAlpha
  omega(:,2)  = Alpha5(:,2)*invSumAlpha
  omega(:,3)  = Alpha5(:,3)*invSumAlpha

  !! Compute interface value at i-1/2
  uMinus = omega(:,1)*W5m(:,1) + omega(:,2)*W5m(:,2) + omega(:,3)*W5m(:,3)
  !! Apply flattening
  uMinus = flat*uMinus + (1.-flat)*data1d(:,ind)

  !! Check for monotonicity
  where ( (uPlus-data1d(:,ind))*(data1d(:,ind)-uMinus) <= 0. )
    uPlus  = data1d(:,ind)
    uMinus = data1d(:,ind)
  end where

end subroutine reconstruct

real function minmod(a,b)
  implicit none
  real :: a,b
  minmod=.5 * (sign(1.,a) + sign(1.,b))*min(abs(a),abs(b))
end function minmod
