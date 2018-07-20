!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexLsVarPoisson
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

subroutine gr_amrexLsVarPoisson (iSoln, iSrc, bcTypes, bcValues, iAlpha, iBeta, ascalar, bscalar)
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
       GRID_PDE_BND_NEUMANN,   &
       GRID_PDE_BND_DIRICHLET
  use amrex_multigrid_module, ONLY : amrex_multigrid, amrex_multigrid_build, amrex_multigrid_destroy, &
                                     amrex_bottom_default, amrex_bottom_hypre, amrex_bottom_bicgstab, &
                                     amrex_bottom_smoother, amrex_bottom_cg
  use amrex_abeclaplacian_module, ONLY : amrex_abeclaplacian, amrex_abeclaplacian_build, amrex_abeclaplacian_destroy
  use amrex_lo_bctypes_module, ONLY : amrex_lo_periodic, amrex_lo_dirichlet, amrex_lo_neumann
  use amrex_amr_module, ONLY : amrex_geom, amrex_get_finest_level, amrex_max_level
  use amrex_fort_module,     ONLY : amrex_real
  use gr_amrexLsData, ONLY : gr_amrexLs_agglomeration, gr_amrexLs_consolidation, &
                                    gr_amrexLs_linop_maxorder, gr_amrexLs_verbose, gr_amrexLs_cg_verbose, &
                                    gr_amrexLs_max_iter, gr_amrexLs_max_fmg_iter
  use gr_physicalMultifabs,  ONLY : unk, facevarx, facevary, facevarz
  !!
  use amrex_multifab_module, ONLY : amrex_multifab, amrex_multifab_destroy, amrex_multifab_build_alias, &
                                  amrex_multifab_build
  use amrex_multifabutil_module, ONLY :  amrex_average_cellcenter_to_face

  use Grid_data, ONLY : gr_meshMe

  implicit none
  
    integer, intent(in)    :: iSoln, iSrc
    integer, intent(in)    :: bcTypes(6)
    real, intent(in)       :: bcValues(2,6)
    integer, intent(in)    :: iAlpha, iBeta
    real, intent(in)       :: ascalar, bscalar
    integer                :: amrexPoissonBcTypes(6)
    integer                :: i
    
    type(amrex_abeclaplacian) :: abeclap
    type(amrex_multigrid) :: multigrid
    integer :: ilev, maxLevel
    real(amrex_real) :: err
    type(amrex_multifab), allocatable, save :: solution(:)
    type(amrex_multifab), allocatable, save :: rhs(:)
    type(amrex_multifab), allocatable, save :: alpha(:)
    type(amrex_multifab), allocatable, save :: beta(:,:)
    type(amrex_multifab) :: null 
type(amrex_multifab), allocatable, save :: beta2(:,:)
logical :: nodal(3)
!    type(amrex_boxarray), allocatable, save :: ba(:)
!    type(amrex_distromap), allocatable, save :: dm(:)
    integer :: idim

!This temp_cc_beta is a cell-centered approx of beta. It is there only because we dont have facevars in flash5 yet.
!this will eventually be delated and instead of interpolating cc_beta to face-cc beta(:,:) we just assign it as
!beta(:,:)=[facevarx(iBeta) facevary(iBeta) facevarz(iBeta)] or beta = facevar(NDIM, iBeta)
!Above assignment should be done by looping over dimension and levels calling amrex_multifab_build_alias(...)
    type(amrex_multifab), allocatable, save :: temp_cc_beta(:)

#include "Flash.h"
#include "constants.h"   
  
  call Timers_start("gr_amrexLsVarPoisson")
     maxLevel = amrex_get_finest_level() !! TODO :: Check with Jared if this is the best wat to get max level of current 
                                                                   !! grid. There is likely a flash ssubroutine for same Grid_getMaxRefinement?
!   Allocate space for multifab array storing phi (solution) and rhs
    allocate(solution(0:maxLevel))
    allocate(rhs(0:maxLevel))
    allocate(alpha(0:maxLevel))
    allocate(beta(1:3,0:maxLevel))
!!!!!!!!!!TEST BEGIN!!!!!!
    allocate(beta2(NDIM,0:maxLevel))
    allocate(temp_cc_beta(0:maxLevel))
!    allocate(ba(0:maxLevel))
!    allocate(dm(0:maxLevel))
    do ilev = 0, maxLevel
!       ba(ilev)=unk(ilev)%ba
!       dm(ilev)=unk(ilev)%dm
       call amrex_multifab_build_alias(temp_cc_beta(ilev), unk(ilev), DENS_VAR, 1)
       do idim = 1, NDIM
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(beta2(idim,ilev), unk(ilev)%ba, unk(ilev)%dm, 1, 0, nodal)
       end do
       call amrex_average_cellcenter_to_face(beta2(:,ilev), temp_cc_beta(ilev), amrex_geom(ilev))
    end do
!!!!!!!!!!TEST END!!!!!!!!!!!

!     Create alias from multifab unk to rhs and solution with respective components in unk
    do ilev = 0, maxLevel
        call amrex_multifab_build_alias(solution(ilev), unk(ilev), iSoln, 1)
        call amrex_multifab_build_alias(rhs(ilev), unk(ilev), iSrc, 1)
        call amrex_multifab_build_alias(alpha(ilev), unk(ilev), iAlpha, 1)
        call amrex_multifab_build_alias(beta(1,ilev), facevarx(ilev), iBeta, 1)
        call amrex_multifab_build_alias(beta(2,ilev), facevary(ilev), iBeta, 1)
#ifdef NDIM > 2
        call amrex_multifab_build_alias(beta(3,ilev), facevarz(ilev), iBeta, 1)
#endif
        call solution(ilev)%setVal(0.0_amrex_real)
    end do

!   Build abeclaplacian object with the geometry amrex_geom, boxarray unk%ba  and distromap unk%dm
       call amrex_abeclaplacian_build(abeclap,amrex_geom(0:maxLevel), rhs%ba, rhs%dm, &
            metric_term=.false., agglomeration=gr_amrexLs_agglomeration, consolidation=gr_amrexLs_consolidation)
!            max_coarsening_level=max_coarsening_level)
       call abeclap % set_maxorder(gr_amrexLs_linop_maxorder)

!  Select BCs to send to AMReX solver
     do i=1,6
       select case (bcTypes(i))
       case (GRID_PDE_BND_PERIODIC)
          amrexPoissonBcTypes(i)=amrex_lo_periodic
       case (GRID_PDE_BND_NEUMANN)
          amrexPoissonBcTypes(i)=amrex_lo_neumann
       case (GRID_PDE_BND_DIRICHLET)
          amrexPoissonBcTypes(i)=amrex_lo_dirichlet
       case default
          call Driver_abortFlash('Invalid or not implemented BC for AMReX ABEC Laplacian solver!')
       end select
     end do
     call abeclap % set_domain_bc([amrexPoissonBcTypes(1),amrexPoissonBcTypes(3),amrexPoissonBcTypes(5)], &
          &                       [amrexPoissonBcTypes(2),amrexPoissonBcTypes(4),amrexPoissonBcTypes(6)])

       do ilev = 0, maxLevel
          ! for problem with pure homogeneous Neumann BC, we could pass an empty multifab otherwise pass solution (ilev)
          if(ALL(amrexPoissonBcTypes==amrex_lo_neumann)) then
             call abeclap % set_level_bc(ilev, null)
          else 
             call abeclap  % set_level_bc(ilev, solution(ilev))
          end if
       end do

       call abeclap % set_scalars(ascalar, bscalar)
       do ilev = 0, maxLevel
          call abeclap % set_acoeffs(ilev, alpha(ilev))
          call abeclap % set_bcoeffs(ilev, beta(:,ilev))
       end do



       call amrex_multigrid_build(multigrid, abeclap)
       call multigrid % set_verbose(gr_amrexLs_verbose)
       call multigrid % set_cg_verbose(gr_amrexLs_cg_verbose)
       call multigrid % set_max_iter(gr_amrexLs_max_iter)
       call multigrid % set_max_fmg_iter(gr_amrexLs_max_fmg_iter)
       !uncomment to invoke following when amrex lib has been updated
!       call multigrid % set_bottom_solver(amrex_bottom_cg)

       err = multigrid % solve(solution, rhs, 1.e-10_amrex_real, 1.e-15_amrex_real)
       if(gr_meshMe==MASTER_PE) then
          print*, err
       endif
!      !!Finalize objects
        do ilev = 0, maxLevel
            call amrex_multifab_destroy(solution(ilev))
            call amrex_multifab_destroy(rhs(ilev))
            call amrex_multifab_destroy(alpha(ilev))
            call amrex_multifab_destroy(beta(1,ilev))
            call amrex_multifab_destroy(beta(2,ilev))
            call amrex_multifab_destroy(beta(3,ilev))

!            call amrex_boxarray_destroy(ba(ilev))
!            call amrex_distromap_destroy(dm(ilev))
            call amrex_multifab_destroy(beta2(1,ilev))
            call amrex_multifab_destroy(beta2(2,ilev))
            call amrex_multifab_destroy(beta2(3,ilev))
            call amrex_multifab_destroy(temp_cc_beta(ilev))
        end do

       call amrex_multigrid_destroy(multigrid)
       call amrex_abeclaplacian_destroy(abeclap)
       
  call Timers_stop("gr_amrexLsVarPoisson")
end subroutine gr_amrexLsVarPoisson
