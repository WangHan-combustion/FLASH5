!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/Grid_solveAbecLaplacian
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
!!  Hardcoded relative and absolute tolerences for multigrid sovle - 1.e-10, 0.0 
!!  (TODO :: should be changed)
!!  ascalar=0, alpha=0, bscalar=-1, beta=1 makes this equivalent to Poisson solver
!!***

subroutine Grid_solveAbecLaplacian (iSoln, iSrc, iAlpha, iBeta, bcTypes, bcValues, poisfact,ascalar,bscalar)
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
       GRID_PDE_BND_NEUMANN,   &
       GRID_PDE_BND_DIRICHLET
  use amrex_multigrid_module, ONLY : amrex_multigrid, amrex_multigrid_build, amrex_multigrid_destroy, amrex_bottom_default
  use amrex_abeclaplacian_module, ONLY : amrex_abeclaplacian, amrex_abeclaplacian_build, amrex_abeclaplacian_destroy
  use amrex_lo_bctypes_module, ONLY : amrex_lo_periodic, amrex_lo_dirichlet, amrex_lo_neumann
  use amrex_amr_module, ONLY : amrex_geom, amrex_get_finest_level, amrex_max_level
  use amrex_fort_module,     ONLY : amrex_real
  use gr_amrexLsData, ONLY : gr_amrexLs_agglomeration, gr_amrexLs_consolidation, &
                                    gr_amrexLs_linop_maxorder, gr_amrexLs_verbose, gr_amrexLs_cg_verbose, &
                                    gr_amrexLs_max_iter,gr_amrexLs_max_fmg_iter,&
                                    gr_amrexLs_composite_solve, gr_amrexLs_ref_ratio,&
                                    gr_amrexLs_max_coarsening_level, gr_amrexLs_bottom_solver
  use gr_physicalMultifabs,  ONLY : unk
  !!
  use amrex_multifab_module, ONLY : amrex_multifab, amrex_multifab_destroy, & 
                                                              amrex_multifab_build_alias, amrex_multifab_build
  use amrex_multifabutil_module, ONLY: amrex_average_cellcenter_to_face
  use amrex_boxarray_module, ONLY : amrex_boxarray, amrex_boxarray_destroy
  use amrex_distromap_module, ONLY : amrex_distromap, amrex_distromap_destroy
  use amrex_geometry_module, ONLY : amrex_geometry, amrex_geometry_destroy

  use Grid_data, only: gr_globalMe

  implicit none
  
    integer, intent(in)    :: iSoln, iSrc, iAlpha, iBeta
    integer, intent(in)    :: bcTypes(6)
    real, intent(in)       :: bcValues(2,6)
    real, intent(in)       :: ascalar, bscalar
    real, intent(inout)    :: poisfact

    integer                :: amrexPoissonBcTypes(6)
    integer                :: i, ilev, maxLevel, idim
    logical :: nodal(3)

    real(amrex_real) :: err
    type(amrex_geometry), allocatable, save :: geom(:)
    type(amrex_boxarray), allocatable, save :: ba(:)
    type(amrex_distromap), allocatable, save :: dm(:)
    type(amrex_abeclaplacian) :: abeclap
    type(amrex_multigrid) :: multigrid
    type(amrex_multifab) :: null
    type(amrex_multifab), allocatable, save :: solution(:)
    type(amrex_multifab), allocatable, save :: rhs(:)
    type(amrex_multifab), allocatable, save :: acoef(:)
    type(amrex_multifab), allocatable, save :: bcoef(:)
    type(amrex_multifab), allocatable :: beta(:,:)

#include "Flash.h"
#include "constants.h"   
  
  call Timers_start("Grid_solveAbecLaplacian")
     maxLevel = amrex_get_finest_level() !! TODO :: Check with Jared if this is the best wat to get max level of current 
                                                                   !! grid. There is likely a flash ssubroutine for same Grid_getMaxRefinement?
!   Allocate space for multifab array storing phi (solution) and rhs
    allocate(solution(0:maxLevel))
    allocate(rhs(0:maxLevel))
    allocate(geom(0:maxLevel))
    allocate(ba(0:maxLevel))
    allocate(dm(0:maxLevel))
    allocate(acoef(0:maxLevel))
    allocate(bcoef(0:maxLevel))

!!initialization
    do ilev = 0, maxLevel
       call amrex_multifab_build_alias(solution(ilev), unk(ilev), NSOL_VAR, 1)
       call amrex_multifab_build_alias(rhs(ilev), unk(ilev), RHS_VAR, 1)
       call amrex_multifab_build_alias(acoef(ilev), unk(ilev), ALPHACC_VAR, 1)
       call amrex_multifab_build_alias(bcoef(ilev), unk(ilev), BETACC_VAR, 1)
       ba(ilev) = rhs(ilev)%ba
       dm(ilev) = rhs(ilev)%dm
       geom(ilev) = amrex_geom(ilev)
      ! This will be used to provide bc and initial guess for the solver.
       call solution(ilev)%setVal(0.0_amrex_real)
    end do

!! For ABecLaplacian, the beta coefficents are on faces. 
!! For the purpose of current version of solver, the beta coefficients are assumed
!! to bedefined on cell center in unk. Hence they are temporarily interpolated to facecentered multifab, 'beta''. 
!!Later, these should be defined on face centers and obtained from facevarx, facevary, facevarz.
    allocate(beta(NDIM,0:maxLevel))
    do ilev = 0, maxLevel
       do idim = 1, NDIM
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(beta(idim,ilev), ba(ilev), dm(ilev), 1, 0, nodal)
       end do
       call amrex_average_cellcenter_to_face(beta(:,ilev), bcoef(ilev), geom(ilev))
!       call amrex_multifab_build_alias(beta(1,ilev), facevarx(ilev), BETA_FACE_VAR, 1)
!       call amrex_multifab_build_alias(beta(2,ilev), facevary(ilev), BETA_FACE_VAR, 1)
!       call amrex_multifab_build_alias(beta(3,ilev), facevarz(ilev), BETA_FACE_VAR, 1)      
    end do
     
!  Select BCs to send to AMReX poisson solver
     do i=1,6
       select case (bcTypes(i))
       case (GRID_PDE_BND_PERIODIC)
          amrexPoissonBcTypes(i)=amrex_lo_periodic
       case (GRID_PDE_BND_NEUMANN)
          amrexPoissonBcTypes(i)=amrex_lo_neumann
       case (GRID_PDE_BND_DIRICHLET)
          amrexPoissonBcTypes(i)=amrex_lo_dirichlet
       case default
          call Driver_abortFlash('BC unsupported or not implemented for AMReX AbecLaplacian solver!')
       end select
     end do
     gr_amrexLs_composite_solve=.TRUE.

  if(gr_amrexLs_composite_solve) then
       call amrex_abeclaplacian_build(abeclap, amrex_geom(0:maxLevel), ba, dm, &
            metric_term=.false., agglomeration=gr_amrexLs_agglomeration, &
            consolidation=gr_amrexLs_consolidation, max_coarsening_level=gr_amrexLs_max_coarsening_level)
       call abeclap % set_maxorder(gr_amrexLs_linop_maxorder)

       ! This is set up to have homogeneous Neumann BC
       call abeclap % set_domain_bc([amrexPoissonBcTypes(1),amrexPoissonBcTypes(3),amrexPoissonBcTypes(5)], &
          &                       [amrexPoissonBcTypes(2),amrexPoissonBcTypes(4),amrexPoissonBcTypes(6)])

       if (ilev > 0) then
          ! use coarse level data to set up bc at corase/fine boundary
          call abeclap % set_coarse_fine_bc(solution(ilev-1), gr_amrexLs_ref_ratio)
       end if

       do ilev = 0, maxLevel
          ! for problem with pure homogeneous Neumann BC, we could pass an empty multifab
          call abeclap % set_level_bc(ilev, solution (ilev))
       end do

       call abeclap % set_scalars(ascalar, bscalar)
       do ilev = 0, maxLevel
          call abeclap % set_acoeffs(ilev, acoef(ilev))
          call abeclap % set_bcoeffs(ilev, beta(:,ilev))
       end do

       call amrex_multigrid_build(multigrid, abeclap)
       call multigrid % set_verbose(gr_amrexLs_verbose)
       call multigrid % set_cg_verbose(gr_amrexLs_cg_verbose)
       call multigrid % set_max_iter(gr_amrexLs_max_iter)
       call multigrid % set_max_fmg_iter(gr_amrexLs_max_fmg_iter)
       gr_amrexLs_bottom_solver = amrex_bottom_default
       call multigrid % set_bottom_solver(gr_amrexLs_bottom_solver)

       err = multigrid % solve(solution, rhs, 1.e-10_amrex_real, 0.0_amrex_real)
       if (gr_globalMe .eq. 0) print*, err
       call amrex_multigrid_destroy(multigrid)
       call amrex_abeclaplacian_destroy(abeclap)
  else
!else do level by level solve instead of composite. Lbl seems to be always ~2x faster than composite
    call Driver_abortFlash('Level by level solve for AbecLaplacian not implemented. Try gr_amrexLs_composite_solve=.TRUE.!')
    do ilev = 0, maxLevel
    end do
  endif
 !!Finalize temporary alias objects
  do ilev = 0, maxLevel
       call amrex_geometry_destroy(geom(ilev))
       call amrex_boxarray_destroy(ba(ilev))
       call amrex_distromap_destroy(dm(ilev))
       call amrex_multifab_destroy(solution(ilev))
       call amrex_multifab_destroy(rhs(ilev))
       call amrex_multifab_destroy(acoef(ilev))
       call amrex_multifab_destroy(bcoef(ilev))
  end do
       
  call Timers_stop("Grid_solveAbecLaplacian")
end subroutine Grid_solveAbecLaplacian
