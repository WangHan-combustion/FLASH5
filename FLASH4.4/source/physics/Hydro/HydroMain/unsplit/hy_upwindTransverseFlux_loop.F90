!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_upwindTransverseFlux
!!
!! NAME
!!
!!  hy_uhd_upwindTransverseFlux
!!
!! SYNOPSIS
!!
!!  call hy_uhd_upwindTransverseFlux( integer(IN) :: dir,
!!                                    integer(IN) :: order,
!!                                    real(IN)    :: vm1(:),
!!                                    real(IN)    :: vc0(:),
!!                                    real(IN)    :: vp1(:),
!!                                    real(IN)    :: lambda(HY_WAVENUM),
!!                                    real(IN)    :: leig(HY_VARINUM,HY_WAVENUM),
!!                                    real(IN)    :: reig(HY_VARINUM,HY_WAVENUM),
!!                                    integer(IN) :: sigSize,
!!                                    real(OUT)   :: sig(sigSize),
!!                                    logical(IN),optional :: speciesScalar)
!!
!!
!! ARGUMENTS
!!
!!  vm1      - input data at i-1 cell
!!  vc0      - input data at i cell
!!  vp1      - input data at i+1 cell
!!  lambda   - eigen values
!!  leig     - left eigen vectors
!!  reig     - right eigen vectors
!!  sigSize  - size of the transverse flux vector
!!  sig      - transverse flux vector
!!  TransFlux - flux vector
!!  delbar    -  gradient of flux
!!  hyBeg    -  Begining indices for range Flux vec
!!  hyConEnd -  Ending indices for range in Flux vec
!!  
!!
!! DESCRIPTION
!!
!!  This routine calculate upwind transverse fluxes for conservative variables.
!!  This is a loop subroutine i.e. a subroutine to encapsulate a compute-intensive loop
!!*** 

#include "Flash.h"

Subroutine hy_uhd_upwindTransverseFlux_loop&
     (vm1,vc0,vp1,lambda,leig,reig,sigSize,sig,TransFlux,delbar,hyBeg,hyConEnd)

  implicit none

#include "UHD.h"

  !!-----Arguments---------------------------------------------------------
  integer, intent(IN) :: hyBeg,hyConEnd
  real,intent(INOUT),dimension(hyBeg:hyConEnd)  :: vm1,vc0 ,vp1
  real,intent(IN),dimension(HY_WAVENUM) :: lambda
  real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: leig
  real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: reig
  integer, intent(IN) :: sigSize
  real,intent(INOUT), dimension(sigSize) :: sig
  real,intent(INOUT),dimension(HY_END_VARS) :: TransFlux
  real,intent(INOUT),dimension(HY_WAVENUM)     :: delbar
  !!------------------------------------------------------------------------

  integer :: n
    !! (1) For hydro/MHD variables:
     !!    (1a) Calculate upwind transverse fluxes for conservative variables
     do n=1,HY_WAVENUM
        ! Upwinding
        ! (NOTE: Using this If-else-endif is much faster than using the signum approach as before)
        if (lambda(n) > 0.) then
           TransFlux(hyBeg:hyConEnd) = vc0(hyBeg:hyConEnd)-vm1(hyBeg:hyConEnd)
        else
           TransFlux(hyBeg:hyConEnd) = vp1(hyBeg:hyConEnd)-vc0(hyBeg:hyConEnd)
        endif

        ! Make sigma sums (or the transverse fluxes) for primitive variables 
        ! except for gamc, game, eint, gravity, and 3T variables.
        ! gamc, game, eint, grav, and 3T are treated separately in the below, (1b).
        delbar(n) = dot_product(leig(hyBeg:hyConEnd,n),TransFlux(hyBeg:hyConEnd))
        TransFlux(hyBeg:hyConEnd) = lambda(n)*reig(hyBeg:hyConEnd,n)*delbar(n)
        sig(hyBeg:hyConEnd) = sig(hyBeg:hyConEnd) + TransFlux(hyBeg:hyConEnd)
     enddo ! End of do n=1,HY_WAVENUM
 
End Subroutine hy_uhd_upwindTransverseFlux_loop