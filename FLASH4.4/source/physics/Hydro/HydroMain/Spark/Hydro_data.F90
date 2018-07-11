!!****if* source/physics/Hydro/HydroMain/Spark/Hydro_data
!!
!!  NAME
!!    Hydro_data
!!
!!  SYNOPSIS
!!    use Hydro_data
!!
!!  DESCRIPTION
!!    Stores data for Spark Hydro
!!
!!  NOTES
!!
!!***
module Hydro_data
#include "constants.h"
#include "Flash.h"
#include "Spark.h"

  implicit none
  save

  integer :: hy_meshMe
  real :: hy_cfl
  logical :: hy_hydroComputeDtFirstCall
  logical :: hy_updateHydroFluxes
  real :: hy_dt, hy_dtmin
  integer :: hy_gcMaskSize
  integer :: hy_globalComm
  integer :: hy_meshNumProcs
  logical :: hy_restart
  logical :: hy_shockDetectOn
  real :: hy_smalldens, hy_smallE, hy_smallpres, hy_smallX

  ! Note that we are making the assumption of fixed-block-size mode!!
  ! Storage for one block's worth of fluxes
  real, target, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: hy_flx, hy_fly, hy_flz

  ! Storage for cell-centered gravitational acceleration
  real, dimension(MDIM,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: hy_grav

  logical :: hy_useHydro
  logical :: hy_fluxCorrect
  integer, dimension(NFLUXES) :: hy_fluxCorVars
  integer :: hy_geometry

  logical :: hy_threadWithinBlock
  logical, dimension(NUNK_VARS) :: hy_gcMask
  ! Additional scratch storage for RK time stepping
  real, allocatable, target :: hy_starState(:,:,:,:)

  ! Limiter info
  real :: hy_limRad
  real :: hy_cvisc

  real :: hy_tiny=1.e-32
  real :: hy_gravConst, hy_4piGinv

  logical :: hy_hybridRiemann, hy_flattening

  real :: hy_C_hyp, hy_alphaGLM, hy_lChyp
  real :: hy_bref
  ! System of units used
  character(4) :: hy_units

end module Hydro_data
