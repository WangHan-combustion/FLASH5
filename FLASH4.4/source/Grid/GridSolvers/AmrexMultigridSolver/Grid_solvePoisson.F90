!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/Grid_solvePoisson
!!
!!  NAME 
!!
!! Grid_solvePoisson
!!
!!  SYNOPSIS
!!
!!  call Grid_solvePoisson()
!!
!!
!!  DESCRIPTION 
!! This routine solves the Poisson equation from the 
!! Amrex Linear Solvers using the variables from Unk multifab
!! for rhs and unknown phi 
!!
!!
!! ARGUMENTS
!!
!!  iSoln    - the index for the solution variable (potential when used for self-gravity)
!!  iSrc     - the index of the source variable (density when used for self-gravity)
!!  bcTypes  - the boundary condition type; only the first entry is used.
!!             Only the first 2*NDIM elements are significant. They are interpreted
!!             in the order (X left, X right, Y left, Y right, Z left, Z right).
!!             Valid values are:
!!               GRID_PDE_BND_PERIODIC (1)
!!               GRID_PDE_BND_DIRICHLET (2) (homogeneous or constant Dirichlet)
!!               GRID_PDE_BND_NEUMANN (3) (homogeneous or constant Neumann)
!!               GRID_PDE_BND_ISOLATED (0)
!!
!!  bcValues - the values to boundary conditions, currently not used (treated as 0)
!!  poisfact      - scaling factor to be used in calculation
!!  iAlpha   - index for alpha matrix in ABECLaplacian equation (cell centered)
!!  iBeta    - index for beta matrix (the coefficients) in ABECLaplacian equation (Face-centered)
!!  ascalar  - ascalar in ABECLaplacian equation 
!!  bscalar  - bscalar in ABECLaplacian equation 
!!   
!!  THE EQUATION TO BE SOLVED   ::  
!!  ascalar*alpha*Soln - bscalar*del dot (beta grad(Soln)) = Src
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  Currently, solver only works for GRID_PDE_BND_PERIODIC i.e. periodic boundary conditions
!!  Other BCs to be implemented later
!!  Relative and absolute tolerences for multigrid sovle - 1.e-10, 0.0
!!
!!***

subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, bcValues, poisfact, iAlpha, iBeta, ascalar, bscalar)
  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none
  
    integer, intent(in)    :: iSoln, iSrc
    integer, intent(in)    :: bcTypes(6)
    real, intent(in)       :: bcValues(2,6)
    real, intent(inout)    :: poisfact
    integer, intent(in), optional :: iAlpha, iBeta
    real, intent(inout), optional :: ascalar, bscalar

#include "Flash.h"
#include "constants.h"   
  
  call Timers_start("Grid_solvePoisson")
  if(present(iBeta) .AND. present(iBeta) ) then
    !! call variable coefficient poisson solver
    if(.NOT.present(ascalar)) ascalar=1.
    if(.NOT.present(bscalar)) bscalar=1.
    call gr_amrexLsVarPoisson(iSoln, iSrc, bcTypes, bcValues, poisfact, iAlpha, iBeta, ascalar, bscalar)
  else
    !! call fixed coefficient poisson solver
    call gr_amrexLsFixedPoisson(iSoln, iSrc, bcTypes, bcValues, poisfact)
  endif
       
  call Timers_stop("Grid_solvePoisson")
end subroutine Grid_solvePoisson
