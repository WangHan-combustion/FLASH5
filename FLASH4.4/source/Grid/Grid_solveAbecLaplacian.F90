!!****if* source/Grid/Grid_solveAbecLaplacian
!!
!!  NAME 
!!
!! Grid_solveAbecLaplacian
!!
!!  SYNOPSIS
!!
!!  call Grid_solveAbecLaplacian()
!!
!!
!!  DESCRIPTION 
!! This routine solves the Abec Laplacian equation from the 
!! Amrex Linear Solvers using the variables from Unk multifab
!! for rhs, unknown coefficients beta and alpha, and scalars a,b
!!Equation to be solved (known as Abec Laplacian) is :: 
!!    !! ascalar * alpha * phi - bscalar * del dot (beta grad phi) = rhs !!
!!
!!
!! ARGUMENTS
!!
!!  iSoln    - the index for the solution (phi) variable (potential when used for self-gravity)
!!  iSrc      - the index of the source (rhs) variable (density when used for self-gravity)
!!  iBeta    - the index of the variable coefficients beta of the Abec Laplacian Eq.
!!  iAlpha  - the index of the coefficients alpha of the Abec Laplacian Eq.
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
!!  ascalar    - ascalar value in ABecLaplacian
!!  bscalar    - bscalar value in ABecLaplacian
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  Currently, solver only works for GRID_PDE_BND_NEUMANN i.e. neumann boundary conditions
!!  Other BCs to be implemented later
!!  ascalar=0, alpha=0, bscalar=-1, beta=1 makes this equivalent to Poisson solver
!!
!!***

subroutine Grid_solveAbecLaplacian (iSoln, iSrc, iAlpha, iBeta, bcTypes, bcValues, poisfact,ascalar,bscalar)
  implicit none
  
    integer, intent(in)    :: iSoln, iSrc, iAlpha, iBeta
    integer, intent(in)    :: bcTypes(6)
    real, intent(in)       :: bcValues(2,6)
    real, intent(in)       :: ascalar, bscalar
    real, intent(inout)    :: poisfact

  return
end subroutine Grid_solveAbecLaplacian
