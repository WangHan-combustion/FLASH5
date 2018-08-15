!!****if* source/Simulation/SimulationMain/unitTest/Multigrid_Amrex
!!
!! NAME
!!
!!  Grid_unitTest
!!
!! SYNOPSIS
!!
!!  call Grid_unitTest(integer(in):: fileUnit,
!!                     logical(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  This unit test is to test the multigrid solver of AMReX library.
!!  It uses known analytical function as well as harmonic
!!
!! ARGUMENTS
!!
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!!
!! NOTES
!!
!!***

!!------!! Do not REORDER(4): solnData

subroutine Grid_unitTest(fileUnit,perfect)

  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET, GRID_PDE_BND_NEUMANN,&
!       Grid_getBlkIndexLimits, &
       Grid_solvePoisson, &
       Grid_getBlkPtr,Grid_releaseBlkPtr, &
       Grid_getDeltas, Grid_fillGuardCells, &
       Grid_getLeafIterator, Grid_releaseLeafIterator
  use gr_interface ,ONLY : gr_findMean
  use Grid_data, ONLY : gr_meshMe, gr_meshComm
  use leaf_iterator, ONLY : leaf_iterator_t
  use block_metadata, ONLY : block_metadata_t
use amrex_amr_module
use amrex_multifab_module
use amrex_boxarray_module
use amrex_distromap_module
use amrex_geometry_module
use amrex_fort_module
use amrex_multigrid_module
use amrex_abeclaplacian_module
use amrex_lo_bctypes_module
  use gr_amrexLsData, ONLY : gr_amrexLs_agglomeration, gr_amrexLs_consolidation, &
                                    gr_amrexLs_linop_maxorder, gr_amrexLs_verbose, gr_amrexLs_cg_verbose, &
                                    gr_amrexLs_max_iter,gr_amrexLs_max_fmg_iter,&
                             gr_amrexLs_composite_solve, gr_amrexLs_ref_ratio
  use gr_physicalMultifabs,      ONLY : unk!, facevarx, facevary, facevarz
 
!,     ONLY : amrex_init_from_scratch, &
!                                   amrex_max_level
!  use gr_amrexLsInterface, ONLY : Grid_amrexLsSolvePoissonUnk

#include "Flash.h"
#include "constants.h"


  implicit none
#include "Flash_mpi.h"

  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors

  real,dimension(:,:,:,:),pointer :: solnData
  type(leaf_iterator_t) :: itor
  type(block_metadata_t) :: block

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(2*MDIM) :: bcTypes
  real,dimension(2,2*MDIM) :: bcValues

  logical :: gcMask(NUNK_VARS), evaluateData

  real :: poisfact
  integer,dimension(MAXBLOCKS) :: blkList
  integer :: blkCount=0,lb,i,j,k
  real:: del(MDIM)
  real meanASOL,meanPFFT
  integer nx,ny,nz
  real, allocatable, dimension(:,:,:) :: fsrc, usol
  integer blkpoints, blkpointsaux,blkCountaux
  real L2_err, L2_erraux, Linf_err, Linf_erraux, Tvol, Tvolaux, vcell
  integer :: refinelevel
  real, parameter :: tol = 1.e-6
  real, parameter :: tolInf = 4.1e-2
  real, parameter :: tol2  = 2.24e-2

  integer TA(2),count_rate,ierr
  real :: ET
!!example test vars..!!
 integer, save :: max_level = 1
  integer, save :: ref_ratio = 2
  integer, save :: n_cell = 32 !128
  integer, save :: max_grid_size = 64
  integer, save :: prob_type = 2
  integer, save :: verbose = 2
  integer, save :: cg_verbose = 0
  integer, save :: max_iter = 100
  integer, save :: max_fmg_iter = 0
  integer, save :: bottom_solver = amrex_bottom_default
  integer, save :: linop_maxorder = 2
  logical, save :: agglomeration = .true.
  logical, save :: consolidation = .true.
  integer, save :: max_coarsening_level = 30

  ! data
  type(amrex_geometry), allocatable, save :: geom(:)
  type(amrex_boxarray), allocatable, save :: ba(:)
  type(amrex_distromap), allocatable, save :: dm(:)
  type(amrex_multifab), allocatable, save :: solution(:)
  type(amrex_multifab), allocatable, save :: rhs(:)
  type(amrex_multifab), allocatable, save :: exact_solution(:)
  type(amrex_multifab), allocatable, save :: acoef(:)
  type(amrex_multifab), allocatable, save :: bcoef(:)
  real(amrex_real), save :: ascalar, bscalar

    integer ilev,idim
    type(amrex_box) :: domain
    type(amrex_box) :: dom

    integer :: rlo(4), rhi(4), elo(4), ehi(4), alo(4), ahi(4), blo(4), bhi(4)
    type(amrex_box) :: bx, gbx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: prhs, pexact, pa, pb

    type(amrex_abeclaplacian) :: abeclap
    type(amrex_multigrid) :: multigrid
    integer(c_long) :: npts
    real(amrex_real) :: err, avg1, avg2, offset
    type(amrex_multifab), allocatable :: beta(:,:)
    logical :: nodal(3)
    type(amrex_multifab) :: null

  ! -------------------------------------------------------------------
  bcTypes(:)=GRID_PDE_BND_NEUMANN  
  bcTypes(3:4)=GRID_PDE_BND_DIRICHLET
  bcValues(:,:)=0.
!!-----------------------------------------------------------!!
!!allocate
    max_level = amrex_get_finest_level()
    allocate(geom(0:max_level))
    allocate(ba(0:max_level))
    allocate(dm(0:max_level))
    allocate(solution(0:max_level))
    allocate(rhs(0:max_level))
    allocate(exact_solution(0:max_level))
       allocate(acoef(0:max_level))
       allocate(bcoef(0:max_level))

!!init mf
    do ilev = 0, max_level
       call amrex_multifab_build_alias(solution(ilev), unk(ilev), NSOL_VAR, 1)
       call amrex_multifab_build_alias(rhs(ilev), unk(ilev), RHS_VAR, 1)
       call amrex_multifab_build_alias(exact_solution(ilev), unk(ilev), ASOL_VAR, 1)
       call amrex_multifab_build_alias(acoef(ilev), unk(ilev), ALPHACC_VAR, 1)
       call amrex_multifab_build_alias(bcoef(ilev), unk(ilev), BETACC_VAR, 1)
       ba(ilev) = rhs(ilev)%ba
       dm(ilev) = rhs(ilev)%dm
       geom(ilev) = amrex_geom(ilev)
    end do

!!init_prob_abeclap
    ascalar=0.
    bscalar=-1.
    do ilev = 0, size(rhs)-1
       !$omp parallel private(rlo,rhi,elo,ehi,alo,ahi,blo,bhi,bx,mfi,prhs,pexact,pa,pb)
       call amrex_mfiter_build(mfi, rhs(ilev), tiling=.true.)

       do while (mfi%next())
          bx = mfi%tilebox()
          gbx = mfi%growntilebox(1)
          prhs   =>            rhs(ilev) % dataptr(mfi)
          pexact => exact_solution(ilev) % dataptr(mfi)
          pa     =>          acoef(ilev) % dataptr(mfi)
          pb     =>          bcoef(ilev) % dataptr(mfi)
          rlo = lbound(prhs)
          rhi = ubound(prhs)
          elo = lbound(pexact)
          ehi = ubound(pexact)
          alo = lbound(pa)
          ahi = ubound(pa)
          blo = lbound(pb)
          bhi = ubound(pb)
          call actual_init_abeclapcian(bx%lo, bx%hi, gbx%lo, gbx%hi, &
               prhs, rlo(1:3), rhi(1:3), pexact, elo(1:3), ehi(1:3), &
               pa, alo(1:3), ahi(1:3), pb, blo(1:3), bhi(1:3), &
               amrex_problo, amrex_probhi, geom(ilev)%dx)
       end do

       call amrex_mfiter_destroy(mfi)
       !$omp end parallel

       ! This will be used to provide bc and initial guess for the solver.
       call solution(ilev)%setVal(0.0_amrex_real)
    end do

!! solve abeclaplacian
    ! For ABecLaplacian, the b coefficents are on faces
    allocate(beta(amrex_spacedim,0:max_level))
    do ilev = 0, max_level
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(beta(idim,ilev), ba(ilev), dm(ilev), 1, 0, nodal)
       end do
       call amrex_average_cellcenter_to_face(beta(:,ilev), bcoef(ilev), geom(ilev))
!       call amrex_multifab_build_alias(beta(1,ilev), facevarx(ilev), BETA_FACE_VAR, 1)
!       call amrex_multifab_build_alias(beta(2,ilev), facevary(ilev), BETA_FACE_VAR, 1)
!       call amrex_multifab_build_alias(beta(3,ilev), facevarz(ilev), BETA_FACE_VAR, 1)      
    end do
     
       call amrex_abeclaplacian_build(abeclap, geom, ba, dm, &
            metric_term=.false., agglomeration=agglomeration, consolidation=consolidation, &
            max_coarsening_level=max_coarsening_level)

       call abeclap % set_maxorder(linop_maxorder)

       ! This is set up to have homogeneous Neumann BC
       call abeclap % set_domain_bc([amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann], &
            &                       [amrex_lo_neumann, amrex_lo_neumann, amrex_lo_neumann])

       if (ilev > 0) then
          ! use coarse level data to set up bc at corase/fine boundary
          call abeclap % set_coarse_fine_bc(solution(ilev-1), ref_ratio)
       end if

       do ilev = 0, max_level
          ! for problem with pure homogeneous Neumann BC, we could pass an empty multifab
          call abeclap % set_level_bc(ilev, solution (ilev))
       end do

       call abeclap % set_scalars(ascalar, bscalar)
       do ilev = 0, max_level
          call abeclap % set_acoeffs(ilev, acoef(ilev))
          call abeclap % set_bcoeffs(ilev, beta(:,ilev))
       end do

       call amrex_multigrid_build(multigrid, abeclap)
       call multigrid % set_verbose(verbose)
       call multigrid % set_cg_verbose(cg_verbose)
       call multigrid % set_max_iter(max_iter)
       call multigrid % set_max_fmg_iter(max_fmg_iter)
       call multigrid % set_bottom_solver(bottom_solver)

  call mpi_barrier(gr_meshComm,ierr)
  if (gr_meshMe .eq. 0) CALL SYSTEM_CLOCK(TA(1),count_rate)  
       err = multigrid % solve(solution, rhs, 1.e-10_amrex_real, 0.0_amrex_real)
  call mpi_barrier(gr_meshComm,ierr)
  if (gr_meshMe .eq. 0) then
     CALL SYSTEM_CLOCK(TA(2),count_rate)
     ET=REAL(TA(2)-TA(1))/count_rate
     write(*,*) ' ' 
     write(*,*) '3 PERIODIC Poisson Solver time = ',ET,' sec.'
  endif
       call amrex_abeclaplacian_destroy(abeclap)
       call amrex_multigrid_destroy(multigrid)

!!Destroy-----------------------------------------------------------!!
    do ilev = 0, max_level
       call amrex_geometry_destroy(geom(ilev))
       call amrex_boxarray_destroy(ba(ilev))
       call amrex_distromap_destroy(dm(ilev))
       call amrex_multifab_destroy(solution(ilev))
       call amrex_multifab_destroy(rhs(ilev))
       call amrex_multifab_destroy(exact_solution(ilev))
       if (allocated(acoef)) then
          call amrex_multifab_destroy(acoef(ilev))
          call amrex_multifab_destroy(bcoef(ilev))
       end if
    end do

!!-----------------------------------------------------------!!
  ! Check error in the solution:
  L2_err = 0.
  blkpoints = 0.
  Linf_err = 0.
  Tvol = 0.
  ! Get Block iterator
!   itor = block_iterator_t(LEAF)
  call Grid_getLeafIterator(itor)
  do while (itor%is_valid())
     call itor%blkMetaData(block)
     !get the index limits of the block
     blkLimits   = block%limits
     blkLimitsGC = block%limitsGC

     ! get a pointer to the current block of data
!     call Grid_getBlkPtr(block, solnData)
     call Grid_getBlkPtr(block,solnData,CENTER)

     call Grid_getDeltas(block%level,del)

     select case (NDIM)
     case(1)
        vcell = del(IAXIS)
     case(2)
        vcell = del(IAXIS)*del(JAXIS)
     case(3)
        vcell = del(IAXIS)*del(JAXIS)*del(KAXIS)
     end select

     blkCount = blkCount + 1
     blkpoints = blkpoints + &
          (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1) * &
          (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1) * &
          (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)

     Tvol = Tvol + vcell*real( (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1) * &
          (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1) * &
          (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1))

     solnData(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),           &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS),DIFF_VAR)   =       &
          abs(  solnData(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS),NSOL_VAR)   - &
          solnData(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS),ASOL_VAR)   )

     ! L2 norm of error:
     L2_err = L2_err + sum( vcell*solnData(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),           &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS),DIFF_VAR)**2.)

     ! Linf norm of error:
     Linf_err = max(Linf_err,maxval(solnData(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS),DIFF_VAR) ))

     call Grid_releaseBlkPtr(block,solnData,CENTER)
     call itor%next()
  enddo
 call Grid_releaseLeafIterator(itor)
! #if defined(__GFORTRAN__) && (__GNUC__ <= 4)
!   call destroy_iterator(itor)
! #endif


  ! Sum processors Volumes
  Tvolaux = Tvol
  call MPI_Allreduce(Tvolaux, Tvol, 1, FLASH_REAL,&
       MPI_SUM, gr_meshComm, ierr)

  ! Sum processors points
  blkpointsaux = blkpoints
  call MPI_Allreduce(blkpointsaux, blkpoints, 1, FLASH_INTEGER,&
       MPI_SUM, gr_meshComm, ierr)

  ! Sum processors L2 norm of error squared
  L2_erraux = L2_err
  call MPI_Allreduce(L2_erraux, L2_err, 1, FLASH_REAL,&
       MPI_SUM, gr_meshComm, ierr)


  ! Compute L2 norm for whole domain
  L2_err = sqrt(L2_err/Tvol)

  ! Sum processors Linf norm of error
  Linf_erraux = Linf_err
  call MPI_Allreduce(Linf_erraux, Linf_err, 1, FLASH_REAL,&
       MPI_MAX, gr_meshComm, ierr)


  ! Fill GuardCells:
  gcMask = .TRUE.                         
  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS,mask=gcMask)          

  ! Export to Tecplot:
  !call outtotecplot(gr_meshMe,0.0,1.,1,0,0.0,blkList,blkCount,0)


   call gr_findMean(ASOL_VAR,2,.false.,meanASOL)
   call gr_findMean(NSOL_VAR,2,.false.,meanPFFT)

  !Unit test gives a mean analytical solution of zero.  Ensure the absolute
  !value of the mean numerical solution is between 0 and the tolerance value.
  perfect = abs(meanPFFT) < MAX(TINY(1.), tol)

  if (perfect) perfect = (abs(Linf_err) .LE. tolInf)
  if (perfect) perfect = (L2_err .LE. tol2)

  ! Sum processors points
  call MPI_Allreduce(blkCount, blkCountaux, 1, FLASH_INTEGER,&
       MPI_SUM, gr_meshComm, ierr)


  if (gr_meshMe .eq. 0) then
     write(*,*) ' ' 
     write(*,'(A,2g16.8)') ' Mean Analytical, Numerical Sol=',meanASOL,meanPFFT
     write(*,'(A,1g16.8)') " ||Phi - PhiAnalytical||inf =" ,Linf_err
     write(*,'(A,1g16.8)') " ||Phi - PhiAnalytical||2   =" ,L2_err
     write(*,*) " Total Volume =",Tvol
     write(*,*) " Total Number of Leaf Blocks=", blkCountaux
     write(*,*) ' ' 
  endif
  return


contains
  subroutine actual_init_abeclapcian (lo, hi, glo, ghi, rhs, rlo, rhi, exact, elo, ehi, &
       alpha, alo, ahi, beta, blo, bhi, prob_lo, prob_hi, dx)
    integer, dimension(3), intent(in) :: lo, hi, glo, ghi, rlo, rhi, elo, ehi, alo, ahi, blo, bhi
    real(amrex_real), intent(inout) :: rhs  (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(inout) :: exact(elo(1):ehi(1),elo(2):ehi(2),elo(3):ehi(3))
    real(amrex_real), intent(inout) :: alpha(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
    real(amrex_real), intent(inout) :: beta (blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3))
    real(amrex_real), dimension(3), intent(in) :: prob_lo, prob_hi, dx

    integer :: i,j,k
    real(amrex_real) x, y, z, xc, yc, zc
    real(amrex_real) r, theta, dbdrfac
    real(amrex_real) pi, fpi, tpi, fac
    real(amrex_real), parameter :: w = 0.05d0
    real(amrex_real), parameter :: sigma = 10.d0
    integer :: bc_type=1
    type(amrex_parmparse) pp

    call amrex_parmparse_build(pp)
    call pp % query("bc_type", bc_type)
    call amrex_parmparse_destroy(pp)

    pi = 4.d0 * atan(1.d0)
    tpi = 2.0d0 * pi
    fpi = 4.0d0 * pi
    fac = 12.d0 * pi**2

    xc = (prob_hi(1) + prob_lo(1))/2.d0
    yc = (prob_hi(2) + prob_lo(2))/2.d0
    zc = (prob_hi(3) + prob_lo(3))/2.d0

    theta = 0.5d0*log(3.d0) / (w + 1.d-50)

    do k = glo(3), ghi(3)
       z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
       do j = glo(2), ghi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i = glo(1), ghi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)

             r = sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)

             beta(i,j,k) = 1.0 !(sigma-1.d0)/2.d0*tanh(theta*(r-0.25d0)) + (sigma+1.d0)/2.d0
          end do
       end do
    end do

    do k = lo(3), hi(3)
       z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
       do j = lo(2), hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i = lo(1), hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
             r = sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)

             dbdrfac = (sigma-1.d0)/2.d0/(cosh(theta*(r-0.25d0)))**2 * theta/r
             dbdrfac = dbdrfac * bscalar

             alpha(i,j,k) = 0. !1.d0

!             exact(i,j,k) = 1.d0 * cos(tpi*x) * cos(tpi*y) * cos(tpi*z)   &
!                  &      + .25d0 * cos(fpi*x) * cos(fpi*y) * cos(fpi*z)
!
!             rhs(i,j,k) = beta(i,j,k)*bscalar*fac*(cos(tpi*x) * cos(tpi*y) * cos(tpi*z)   &
!                  &                        + cos(fpi*x) * cos(fpi*y) * cos(fpi*z))  &
!                  &   + dbdrfac*((x-xc)*(tpi*sin(tpi*x) * cos(tpi*y) * cos(tpi*z)   &
!                  &                     + pi*sin(fpi*x) * cos(fpi*y) * cos(fpi*z))  &
!                  &            + (y-yc)*(tpi*cos(tpi*x) * sin(tpi*y) * cos(tpi*z)   &
!                  &                     + pi*cos(fpi*x) * sin(fpi*y) * cos(fpi*z))  &
!                  &            + (z-zc)*(tpi*cos(tpi*x) * cos(tpi*y) * sin(tpi*z)   &
!                  &                     + pi*cos(fpi*x) * cos(fpi*y) * sin(fpi*z))) &
!                  &                   + ascalar * (cos(tpi*x) * cos(tpi*y) * cos(tpi*z)   &
!                  &               + 0.25d0 * cos(fpi*x) * cos(fpi*y) * cos(fpi*z))
!!!!!!!!!!!!!!
              exact(i,j,k) = cos(tpi*x)*cos(tpi*y)*cos(tpi*z)
              rhs(i,j,k)   = -3.0*tpi*tpi*exact(i,j,k)
!!!!!!!!
          end do
       end do
    end do

  end subroutine actual_init_abeclapcian



end subroutine Grid_unitTest

