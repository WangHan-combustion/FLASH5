!!****if* source/physics/Hydro/HydroMain/M1/hy_rk_getFaceFlux
!!
!!  NAME
!!
!!  hy_rk_getFaceFlux
!!
!!  SYNOPSIS
!!
!!  call hy_rk_getFaceFlux ( integer(IN) :: blockID )
!!
!!  DESCRIPTION
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine hy_rk_getFaceFlux (blockID, limits)

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getDeltas, &
       Grid_getBlkData, Grid_putFluxData
  use Hydro_data, ONLY : hy_fluxCorrect, hy_geometry, hy_fluxCorVars, &
       hy_flx, hy_fly, hy_flz, hy_threadWithinBlock, hy_limRad, &
       hy_starState, hy_grav, hy_flattening
  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Spark.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS

  integer, intent(IN) :: blockID
  integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  integer :: i1,i2,i3, n, g, dir, ierr

  real, dimension(NFLUXES) :: flux

  real :: faceAreas(GRID_ILO_GC:GRID_IHI_GC,     &
       GRID_JLO_GC:GRID_JHI_GC,     &
       GRID_KLO_GC:GRID_KHI_GC)

  integer, dimension(MDIM) :: datasize
  real, dimension(MDIM) :: del

  ! 1D pencil of data
  integer :: size1d
  real, allocatable :: pencil(:,:)
  integer, dimension(LOW:HIGH,MDIM) :: dirLims
  real, dimension(NRECON) :: leftState, rightState
  real, dimension(NRECON) :: uPlus, uMinus
  real :: speed
  real, allocatable :: grv(:), shck(:), flat(:)
  logical :: inShock
  real :: flat3d(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getDeltas(blockID,del)

  datasize(1:MDIM) = blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
  size1d = maxval(datasize)
  allocate(pencil(NRECON,size1d))
  allocate(grv(size1d))
  allocate(shck(size1d))
  allocate(flat(size1d))

  if (hy_flattening) then
     call flattening(flat3d, limits)
  else
     flat3d = 1.0
  end if

  !$omp parallel if(hy_threadWithinBlock .AND. NDIM > 1) &
  !$omp default(none) &
  !$omp private(i1,i2,i3,n,grv,shck,flat,inShock,flux,ierr,speed,&
  !$omp         leftState,rightState,uPlus,uMinus,pencil,dir,dirLims)&
  !$omp shared(hy_starState,blockID,size1d,del,flat3d)

  !  Begin loop over zones
  do dir = 1, NDIM
     call setLoop(dir, dirLims)
     !$omp do schedule(guided) collapse(2)
     do i3 = dirLims(LOW,3), dirLims(HIGH,3)
        do i2 = dirLims(LOW,2), dirLims(HIGH,2)
           call setPencil(pencil,grv,shck,flat,i2,i3,dir)
           !call setHSE(pencil,grv,dirLims(LOW,1)-1,del(dir))
           i1 = dirLims(LOW,1)-1
           call reconstruct(uPlus, uMinus, &
                pencil, flat(i1), size1d, i1, del(dir))
           !call resetHSE(pencil,grv,dirLims(LOW,1)-1,del(dir))
           do i1 = dirLims(LOW,1), dirLims(HIGH,1)+1
              ! cycle left-right state vectors
              leftState = uPlus
              ! Now do the reconstruction for ALL variables in one call
              !call setHSE(pencil,grv,i1,del(dir))
              call reconstruct(uPlus, uMinus, &
                   pencil, flat(i1), size1d, i1, del(dir))
              !call resetHSE(pencil,grv,i1,del(dir))
              rightState = uMinus

              ! Make sure that reconstruction has not introduced unphysical states
              call ensurePhysicalState(leftState)
              call ensurePhysicalState(rightState)

              ! Check for shocks in the zones involved in flux calculation
              inShock = any(shck(i1-1:i1) /= 0.0)
              ! Now call the Riemann solver to compute fluxes
              call riemann(dir,leftState(1:HY_NUM_VARS),&
                   rightState(1:HY_NUM_VARS),inShock,&
                   flux(1:HY_NUM_FLUX),speed,ierr)
              ! Add artificial viscosity for strong-shock capturing
              call avisc(flux,pencil,i1,NRECON,dir)

              ! Here, we compute the species and mass scalar
              ! fluxes based on the density flux and the reconstructed
              ! mass scalar interface values
              call mscalarFluxes(flux,leftState(HY_NUM_VARS+1:NRECON),&
                   rightState(HY_NUM_VARS+1:NRECON))

              ! ***************
              ! Fluxes computed for one face of this zone
              ! Save the fluxes
              ! ***************
              call saveFluxes(dir,i1,i2,i3,flux)
           end do  ! i
        end do ! j
     end do ! k
     !$omp end do nowait
  end do ! dir
  !$omp end parallel

  deallocate(pencil)
  deallocate(grv)
  deallocate(shck)
  deallocate(flat)


contains

  subroutine setLoop(dir, dirLims)
    implicit none
    integer, intent(IN) :: dir
    integer, intent(OUT), dimension(LOW:HIGH,MDIM) :: dirLims
    select case(dir)
    case (IAXIS)
       dirLims(:,:) = limits(:,:)
    case (JAXIS)
       dirLims(:,1) = limits(:,JAXIS)
       dirLims(:,2) = limits(:,IAXIS)
       dirLims(:,3) = limits(:,KAXIS)
    case (KAXIS)
       dirLims(:,1) = limits(:,KAXIS)
       dirLims(:,2) = limits(:,IAXIS)
       dirLims(:,3) = limits(:,JAXIS)
    end select
  end subroutine setLoop

  !! the 'pencil' holds a 1D array of the solution data to be operated
  !! on.  It is unrolled this way so that all data that are needed for
  !! interpolation and flux calculation are truly contiguous in memory
  !! space.  And so that it all fits in cache at once...
  !! The maximum amount of calculation is done on these data prior to
  !! reseting the pencil data to a new ray through the block.
  !! If I were talented at ASCII art, I would make a diagram...
  subroutine setPencil(pencil,grv,shck,flat,i2,i3,dir)
    implicit none
    real, intent(INOUT) :: pencil(HY_NUM_VARS+NSPECIES+NMASS_SCALARS,size1d)
    real, intent(INOUT) :: grv(size1d)
    real, intent(INOUT) :: shck(size1d)
    real, intent(INOUT) :: flat(size1d)
    integer, intent(IN) :: i2,i3,dir
    integer :: n
    select case(dir)
    case (IAXIS)
       pencil(HY_DENS,:) = hy_starState(DENS_VAR,:,i2,i3)
       pencil(HY_VELX,:) = hy_starState(VELX_VAR,:,i2,i3)
       pencil(HY_VELY,:) = hy_starState(VELY_VAR,:,i2,i3)
       pencil(HY_VELZ,:) = hy_starState(VELZ_VAR,:,i2,i3)
       pencil(HY_PRES,:) = hy_starState(PRES_VAR,:,i2,i3)
       pencil(HY_GAMC,:) = hy_starState(GAMC_VAR,:,i2,i3)
       pencil(HY_RHOE,:) = hy_starState(DENS_VAR,:,i2,i3)*hy_starState(EINT_VAR,:,i2,i3)
#ifdef SPARK_GLM
       pencil(HY_MAGX,:) = hy_starState(MAGX_VAR,:,i2,i3)
       pencil(HY_MAGY,:) = hy_starState(MAGY_VAR,:,i2,i3)
       pencil(HY_MAGZ,:) = hy_starState(MAGZ_VAR,:,i2,i3)
       pencil(HY_PSIB,:) = hy_starState(PSIB_VAR,:,i2,i3)
#endif
#if NSPECIES+NMASS_SCALARS>0
       do n=SPECIES_BEGIN, MASS_SCALARS_END
          pencil(HY_NUM_VARS+1+n-SPECIES_BEGIN,:)    = hy_starState(n,:,i2,i3)
       enddo
#endif
#ifdef GRAVITY
#ifdef GPOT_VAR
       grv(:) = hy_starState(GPOT_VAR,:,i2,i3)
#else
       grv(:) = hy_grav(IAXIS,:,i2,i3)
#endif
#endif
#ifdef SHOK_VAR
       shck(:) = hy_starState(SHOK_VAR,:,i2,i3)
#else
       shck(:) = 0.0
#endif
       flat(:) = flat3d(:,i2,i3)
    case (JAXIS)
       pencil(HY_DENS,:) = hy_starState(DENS_VAR,i2,:,i3)
       pencil(HY_VELX,:) = hy_starState(VELX_VAR,i2,:,i3)
       pencil(HY_VELY,:) = hy_starState(VELY_VAR,i2,:,i3)
       pencil(HY_VELZ,:) = hy_starState(VELZ_VAR,i2,:,i3)
       pencil(HY_PRES,:) = hy_starState(PRES_VAR,i2,:,i3)
       pencil(HY_GAMC,:) = hy_starState(GAMC_VAR,i2,:,i3)
       pencil(HY_RHOE,:) = hy_starState(DENS_VAR,i2,:,i3)*hy_starState(EINT_VAR,i2,:,i3)
#ifdef SPARK_GLM
       pencil(HY_MAGX,:) = hy_starState(MAGX_VAR,i2,:,i3)
       pencil(HY_MAGY,:) = hy_starState(MAGY_VAR,i2,:,i3)
       pencil(HY_MAGZ,:) = hy_starState(MAGZ_VAR,i2,:,i3)
       pencil(HY_PSIB,:) = hy_starState(PSIB_VAR,i2,:,i3)
#endif
#if NSPECIES+NMASS_SCALARS>0
       do n=SPECIES_BEGIN, MASS_SCALARS_END
          pencil(HY_NUM_VARS+1+n-SPECIES_BEGIN,:)    = hy_starState(n,i2,:,i3)
       enddo
#endif
#ifdef GRAVITY
#ifdef GPOT_VAR
       grv(:) = hy_starState(GPOT_VAR,i2,:,i3)
#else
       grv(:) = hy_grav(JAXIS,i2,:,i3)
#endif
#endif
#ifdef SHOK_VAR
       shck(:) = hy_starState(SHOK_VAR,i2,:,i3)
#else
       shck(:) = 0.0
#endif
       flat(:) = flat3d(i2,:,i3)
    case (KAXIS)
       pencil(HY_DENS,:) = hy_starState(DENS_VAR,i2,i3,:)
       pencil(HY_VELX,:) = hy_starState(VELX_VAR,i2,i3,:)
       pencil(HY_VELY,:) = hy_starState(VELY_VAR,i2,i3,:)
       pencil(HY_VELZ,:) = hy_starState(VELZ_VAR,i2,i3,:)
       pencil(HY_PRES,:) = hy_starState(PRES_VAR,i2,i3,:)
       pencil(HY_GAMC,:) = hy_starState(GAMC_VAR,i2,i3,:)
       pencil(HY_RHOE,:) = hy_starState(DENS_VAR,i2,i3,:)*hy_starState(EINT_VAR,i2,i3,:)
#ifdef SPARK_GLM
       pencil(HY_MAGX,:) = hy_starState(MAGX_VAR,i2,i3,:)
       pencil(HY_MAGY,:) = hy_starState(MAGY_VAR,i2,i3,:)
       pencil(HY_MAGZ,:) = hy_starState(MAGZ_VAR,i2,i3,:)
       pencil(HY_PSIB,:) = hy_starState(PSIB_VAR,i2,i3,:)
#endif
#if NSPECIES+NMASS_SCALARS>0
       do n=SPECIES_BEGIN, MASS_SCALARS_END
          pencil(HY_NUM_VARS+1+n-SPECIES_BEGIN,:)    = hy_starState(n,i2,i3,:)
       enddo
#endif
#ifdef GRAVITY
#ifdef GPOT_VAR
       grv(:) = hy_starState(GPOT_VAR,i2,i3,:)
#else
       grv(:) = hy_grav(KAXIS,i2,i3,:)
#endif
#endif
#ifdef SHOK_VAR
       shck(:) = hy_starState(SHOK_VAR,i2,i3,:)
#else
       shck(:) = 0.0
#endif
       flat(:) = flat3d(i2,i3,:)
    end select
  end subroutine setPencil


  subroutine saveFluxes(dir,i1,i2,i3,flux)
    use Hydro_data, ONLY : hy_flx, hy_fly, hy_flz
    integer, intent(IN) :: dir,i1,i2,i3
    real, dimension(NFLUXES), intent(IN) :: flux
    select case(dir)
    case(IAXIS)
       hy_flx(:,i1,i2,i3) = flux(:)
    case (JAXIS)
       hy_fly(:,i2,i1,i3) = flux(:)
    case (KAXIS)
       hy_flz(:,i2,i3,i1) = flux(:)
    end select
  end subroutine saveFluxes

  subroutine mscalarFluxes(flux,XL,XR)
    implicit none
    real, intent(INOUT) :: flux(NFLUXES)
    real, target, dimension(NSPECIES+NMASS_SCALARS), intent(IN) :: XL, XR
    real, pointer :: Xstar(:)
#if NSPECIES+NMASS_SCALARS==0
    return
#else
    if (flux(HY_MASS) > 0.) then
       Xstar => XL
    else
       Xstar => XR
    endif
    flux(HY_NUM_FLUX+1:NFLUXES) = Xstar*flux(HY_MASS)
#endif
  end subroutine mscalarFluxes

  subroutine setHSE(V,grv,i1,dx)
    implicit none
    real, intent(INOUT) :: V(HY_NUM_VARS+NSPECIES+NMASS_SCALARS,size1d)
    real, intent(IN) :: grv(size1d)
    integer, intent(IN) :: i1
    real, intent(IN) :: dx
    real :: accelM, accelP
#ifndef GRAVITY
    return
#else
#ifdef GPOT_VAR
    accelM = grv(i1)   - grv(i1-1)
    accelP = grv(i1+1) - grv(i1)
#else
    accelM = 0.5*(grv(i1-1) + grv(i1)  )*dx
    accelP = 0.5*(grv(i1)   + grv(i1+1))*dx
#endif
    V(HY_PRES,i1-1) = V(HY_PRES,i1-1) - &
         0.5*(V(HY_DENS,i1-1)+V(HY_DENS,i1)) * accelM
    V(HY_PRES,i1+1) = V(HY_PRES,i1+1) + &
         0.5*(V(HY_DENS,i1)+V(HY_DENS,i1+1)) * accelP
#endif /* GRAVITY */
  end subroutine setHSE

  subroutine resetHSE(V,grv,i1,dx)
    implicit none
    real, intent(INOUT) :: V(HY_NUM_VARS+NSPECIES+NMASS_SCALARS,size1d)
    real, intent(IN) :: grv(size1d)
    integer, intent(IN) :: i1
    real, intent(IN) :: dx
    real :: accelM, accelP
#ifndef GRAVITY
    return
#else
#ifdef GPOT_VAR
    accelM = grv(i1)   - grv(i1-1)
    accelP = grv(i1+1) - grv(i1)
#else
    accelM = 0.5*(grv(i1-1) + grv(i1)  )*dx
    accelP = 0.5*(grv(i1)   + grv(i1+1))*dx
#endif
    V(HY_PRES,i1-1) = V(HY_PRES,i1-1) + &
         0.5*(V(HY_DENS,i1-1)+V(HY_DENS,i1)) * accelM
    V(HY_PRES,i1+1) = V(HY_PRES,i1+1) - &
         0.5*(V(HY_DENS,i1)+V(HY_DENS,i1+1)) * accelP
#endif /* GRAVITY */
  end subroutine resetHSE

  !! Artificial viscosity as in Colella and Woodward
  !! Calculation of cvisc here is different from the UHD solver
  !! I am not including the transverse velocity terms to keep everything
  !! neat and 1D for the pencil data
  subroutine avisc(flux,V,i1,nvars,dir)
    use Hydro_data, ONLY : hy_cvisc
    implicit none
    real, intent(INOUT) :: flux(NFLUXES)
    integer, intent(IN) :: nvars
    real, intent(IN) :: V(nvars,size1d)
    integer, intent(IN) :: i1,dir
    real :: cvisc, VenerLo, VenerHi

    cvisc = hy_cvisc*max(-(V(HY_VELX+dir-1,i1) - V(HY_VELX+dir-1,i1-1)),0.)

    ! Construct minus and plus TOTAL energy densities
    VenerLo = V(HY_DENS,i1-1)*0.5*(dot_product(V(HY_VELX:HY_VELZ,i1-1),V(HY_VELX:HY_VELZ,i1-1)))&
            + V(HY_RHOE,i1-1)
    VenerHi = V(HY_DENS,i1)*0.5*(dot_product(V(HY_VELX:HY_VELZ,i1),V(HY_VELX:HY_VELZ,i1)))&
            + V(HY_RHOE,i1)

    flux(HY_MASS:HY_ENER) = &
         flux(HY_MASS:HY_ENER) &
         +cvisc*(/V(HY_DENS,i1-1)                 - V(HY_DENS,i1)&
         ,        V(HY_DENS,i1-1)*V(HY_VELX,i1-1) - V(HY_DENS,i1)*V(HY_VELX,i1)&
         ,        V(HY_DENS,i1-1)*V(HY_VELY,i1-1) - V(HY_DENS,i1)*V(HY_VELY,i1)&
         ,        V(HY_DENS,i1-1)*V(HY_VELZ,i1-1) - V(HY_DENS,i1)*V(HY_VELZ,i1)&
         ,        VenerLo                         - VenerHi/)
#ifdef SPARK_GLM
    flux(HY_FMGX:HY_FPSI) = &
         flux(HY_FMGX:HY_FPSI) &
         +cvisc*(/V(HY_MAGX,i1-1)                 - V(HY_MAGX,i1)&
         ,        V(HY_MAGY,i1-1)                 - V(HY_MAGY,i1)&
         ,        V(HY_MAGZ,i1-1)                 - V(HY_MAGZ,i1)&
         ,        V(HY_PSIB,i1-1)                 - V(HY_PSIB,i1)/)
#endif
  end subroutine avisc

  subroutine flattening(flat3d,limits)
    !! This follows Miller & Colella 2002
    use Hydro_data, ONLY : hy_starState
    implicit none
    integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
    real, intent(OUT) :: flat3d(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
    real :: flatTilde(NDIM,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
    real :: beta, Z
    real, parameter :: betaMin = 0.75, betaMax = 0.85
    real, parameter :: Zmin = 0.25, Zmax = 0.75
    integer :: i,j,k, kx, ky, kz
    kx = 1
    ky = 0
    kz = 0
#if NDIM>1
    ky = 1
#if NDIM==3
    kz = 1
#endif
#endif
    do k=limits(LOW,KAXIS)-kz, limits(HIGH,KAXIS)+kz
       do j=limits(LOW,JAXIS)-ky, limits(HIGH,JAXIS)+ky
          do i=limits(LOW,IAXIS)-kx, limits(HIGH,IAXIS)+kx
             beta = abs(hy_starState(PRES_VAR,i+1,j,k)-hy_starState(PRES_VAR,i-1,j,k)) &
                  / max(TINY(1.0), abs(hy_starState(PRES_VAR,i+2,j,k)-hy_starState(PRES_VAR,i-2,j,k)))
             Z    = abs(hy_starState(PRES_VAR,i+1,j,k)-hy_starState(PRES_VAR,i-1,j,k)) &
                  / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
             flatTilde(IAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
             if (hy_starState(VELX_VAR,i+1,j,k)<hy_starState(VELX_VAR,i,j,k)) then
                flatTilde(IAXIS,i,j,k) = max(flatTilde(IAXIS,i,j,k), &
                     &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
             else
                flatTilde(IAXIS,i,j,k) = 1.0
             end if
#if NDIM>1
             beta = abs(hy_starState(PRES_VAR,i,j+1,k)-hy_starState(PRES_VAR,i,j-1,k)) &
                  / max(TINY(1.0),abs(hy_starState(PRES_VAR,i,j+2,k)-hy_starState(PRES_VAR,i,j-2,k)))
             Z    = abs(hy_starState(PRES_VAR,i,j+1,k)-hy_starState(PRES_VAR,i,j-1,k)) &
                  / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
             flatTilde(JAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
             if (hy_starState(VELY_VAR,i,j+1,k)<hy_starState(VELY_VAR,i,j,k)) then
                flatTilde(JAXIS,i,j,k) = max(flatTilde(JAXIS,i,j,k), &
                     &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
             else
                flatTilde(JAXIS,i,j,k) = 1.0
             end if
#if NDIM==3
             beta = abs(hy_starState(PRES_VAR,i,j,k+1)-hy_starState(PRES_VAR,i,j,k-1)) &
                  / max(TINY(1.0),abs(hy_starState(PRES_VAR,i,j,k+2)-hy_starState(PRES_VAR,i,j,k-2)))
             Z    = abs(hy_starState(PRES_VAR,i,j,k+1)-hy_starState(PRES_VAR,i,j,k-1)) &
                  / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
             flatTilde(KAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
             if (hy_starState(VELZ_VAR,i,j,k+1)<hy_starState(VELZ_VAR,i,j,k)) then
                flatTilde(KAXIS,i,j,k) = max(flatTilde(KAXIS,i,j,k), &
                     &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
             else
                flatTilde(KAXIS,i,j,k) = 1.0
             end if
#endif
#endif
          end do
       end do
    end do
    do k=limits(LOW,KAXIS)-kz, limits(HIGH,KAXIS)+kz
       do j=limits(LOW,JAXIS)-ky, limits(HIGH,JAXIS)+ky
          do i=limits(LOW,IAXIS)-kx, limits(HIGH,IAXIS)+kx
             flat3d(i,j,k) = minval(flatTilde(1:NDIM,i,j,k))
#ifdef FLAT_VAR
             hy_starState(FLAT_VAR,i,j,k) = flat3d(i,j,k)
#endif
          end do
       end do
    end do
  end subroutine flattening

  subroutine ensurePhysicalState(state)
    use Hydro_data, ONLY : hy_smalldens, hy_smallE, hy_smallpres, hy_smallX
    implicit none
    real, dimension(NRECON), intent(INOUT), target  :: state
    integer :: s
    real :: spcSumInv
    real, pointer :: spc(:)
    state(HY_DENS) = max(hy_smalldens, state(HY_DENS))
    state(HY_PRES) = max(hy_smallpres, state(HY_PRES))
#if NSPECIES>0
    ! Limit and renormalize the species.
    spc => state(HY_NUM_VARS+1:HY_NUM_VARS+NSPECIES)
    do s = 1, NSPECIES
       spc(s) = max(hy_smallX,min(1.0,spc(s)))
    end do
    spcSumInv = 1./sum(spc(1:NSPECIES))
    spc = spc*spcSumInv
#endif
  end subroutine ensurePhysicalState

#include "reconstruct.F90"
#include "riemann.F90"

end subroutine hy_rk_getFaceFlux
