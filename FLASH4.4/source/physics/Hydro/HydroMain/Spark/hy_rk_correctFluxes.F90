!!****if* source/physics/Hydro/HydroMain/Spark/hy_rk_correctFluxes
!!
!!  NAME
!!
!!  hy_rk_correctFluxes
!!
!!  SYNOPSIS
!!
!!  call hy_rk_correctFluxes ( integer(IN) :: blockID )
!!
!!  DESCRIPTION
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine hy_rk_correctFluxes(blockID, dt)

  use Hydro_data, ONLY : hy_threadWithinBlock, hy_starState, &
       hy_flx, hy_fly, hy_flz, hy_smallE, hy_smalldens, hy_geometry,&
       hy_grav, hy_4piGinv, hy_alphaGLM, hy_C_hyp, hy_fluxCorVars
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_getDeltas, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getBlkData, Grid_getFluxData
  use Eos_interface, ONLY : Eos_putData, Eos_getData, Eos


  implicit none

#include "Flash.h"
#include "constants.h"
#include "Spark.h"
#include "Eos.h"

  integer, intent(IN) :: blockID
  real, intent(IN) :: dt

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(GRID_IHI_GC) :: xCenter, xLeft, xRight
  integer :: i,j,k,n,g

  real, pointer :: solnData(:,:,:,:)
  real, pointer :: Vstar(:)

  real :: dx, dy, dz, del(MDIM)
  real :: dFlux(NFLUXES)

  real :: eint, ekin, emag
  ! Geometry factors
  real :: facM, facP
  integer :: isize, jsize, ksize
  integer, dimension(MDIM) :: datasize
  real :: dhdt, fac
  real, dimension(GRID_ILO_GC:GRID_IHI_GC,     &
       GRID_JLO_GC:GRID_JHI_GC,     &
       GRID_KLO_GC:GRID_KHI_GC) :: faceAreas, cellVolumes

  ! For EOS call
  integer,dimension(MDIM) :: pos
  real, dimension(NSPECIES*MAXCELLS) :: massFraction
  real, dimension(EOS_NUM*MAXCELLS) :: eosData

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  iSize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  jSize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  kSize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
  datasize(1:MDIM) = blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

  call Grid_getBlkPtr(blockID,solnData)

  call Grid_getDeltas(blockID,del)
  dhdt = minval(del(1:NDIM))/dt

  !! Grab the conservative flux differences
  if (hy_geometry /= CARTESIAN) then
     call Grid_getBlkData(blockID,CELL_FACEAREA,ILO_FACE, EXTERIOR, &
          (/1,1,1/), faceAreas, datasize)
     call Grid_getFluxData(blockID,IAXIS,hy_flx,datasize,hy_fluxCorVars,faceAreas)
  else
     call Grid_getFluxData(blockID,IAXIS,hy_flx,datasize)
  endif
#if NDIM > 1
  if (hy_geometry /= CARTESIAN) then
     call Grid_getBlkData(blockID,CELL_FACEAREA,JLO_FACE, EXTERIOR, &
          (/1,1,1/), faceAreas, datasize)
     call Grid_getFluxData(blockID,JAXIS,hy_fly,datasize,hy_fluxCorVars,faceAreas)
  else
     call Grid_getFluxData(blockID,JAXIS,hy_fly,datasize)
  endif
#if NDIM > 2
  if (hy_geometry /= CARTESIAN) then
     call Grid_getBlkData(blockID,CELL_FACEAREA,KLO_FACE, EXTERIOR, &
          (/1,1,1/), faceAreas, datasize)
     call Grid_getFluxData(blockID,KAXIS,hy_flz,datasize,hy_fluxCorVars,faceAreas)
  else
     call Grid_getFluxData(blockID,KAXIS,hy_flz,datasize)
  endif
#endif
#endif

  if (hy_geometry /= CARTESIAN) then
     ! most of the following could and should be moved to geoFacs
     call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, &
          (/1,1,1/), faceAreas, datasize )
     call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR, &
          (/1,1,1/), cellVolumes, (/isize, jsize, ksize/) )
     call Grid_getCellCoords(IAXIS,blockID, CENTER,.true.,xCenter, blkLimitsGC(HIGH,IAXIS))
     call Grid_getCellCoords(IAXIS,blockID, LEFT_EDGE,.true.,xLeft, blkLimitsGC(HIGH,IAXIS))
     call Grid_getCellCoords(IAXIS,blockID, RIGHT_EDGE,.true.,xRight, blkLimitsGC(HIGH,IAXIS))
  endif

  !$omp parallel if (.FALSE.) &
  !$omp default(none) &
  !$omp private(n,i,j,k,Vstar,facM,facP,pos,eosData,massFraction,&
  !$omp         emag,ekin,dx,dy,dz,fac,dFlux)&
  !$omp shared(dt,solnData,hy_starState,hy_flx,hy_fly,hy_flz,hy_grav,&
  !$omp        xCenter,xLeft,xRight,blockID,&
  !$omp        blkLimits,blkLimitsGC,hy_alphaGLM, hy_C_hyp,&
  !$omp        dhdt, hy_smalldens, hy_smallE,del)


  ! Correct IAXIS sides
  !$omp do schedule(static) collapse(2)
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        !! LOW side
        i = blkLimits(LOW,IAXIS)
        ! Point to old/intermediate states
        Vstar => solnData(:,i,j,k)
        ! Point to the correct fluxes
        dFlux = -hy_flx(:,i  ,j  ,k  )
        dx = del(IAXIS)
        ! Get gemoetric factors and sources
        call geoFacs(i,j,k,facM,facP)
        fac = facM
        ! if (dFlux(HY_ENER) /= 0.) print *, 'b', Vstar(TEMP_VAR), dt/dx*fac*dFlux(HY_ENER)/(Vstar(ENER_VAR)*Vstar(DENS_VAR)), Vstar(VELX_VAR)
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        pos = (/i,j,k/)
        call Eos_getData(IAXIS,pos,1,solnData,CENTER,eosData,massFraction)
        call Eos(MODE_DENS_EI,1,eosData,massFraction)
        call Eos_putData(IAXIS,pos,1,solnData,CENTER,eosData)
        ! if (dFlux(HY_ENER) /= 0.) print *, 'a', Vstar(TEMP_VAR), dt/dx*fac*dFlux(HY_ENER)/(Vstar(ENER_VAR)*Vstar(DENS_VAR)), Vstar(VELX_VAR)
        ! Release pointers
        nullify(Vstar)

        !! HIGH side
        i = blkLimits(HIGH,IAXIS)
        ! Point to old/intermediate states
        Vstar => solnData(:,i,j,k)
        ! Point to the correct fluxes
        dFlux = hy_flx(:,i+1  ,j  ,k  )
        dx = del(IAXIS)
        ! Get gemoetric factors and sources
        call geoFacs(i,j,k,facM,facP)
        fac = facP
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        pos = (/i,j,k/)
        call Eos_getData(IAXIS,pos,1,solnData,CENTER,eosData,massFraction)
        call Eos(MODE_DENS_EI,1,eosData,massFraction)
        call Eos_putData(IAXIS,pos,1,solnData,CENTER,eosData)
        ! Release pointers
        nullify(Vstar)
     enddo !j
  enddo !k
  !$omp end do
#if NDIM>1
  ! Correct JAXIS sides
  !$omp do schedule(static) collapse(2)
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
        !! LOW side
        j = blkLimits(LOW,JAXIS)
        ! Point to old/intermediate states
        Vstar => solnData(:,i,j,k)
        ! Point to the correct fluxes
        dFlux = -hy_fly(:,i  ,j  ,k  )
        fac = 1.0
        dx = del(JAXIS)
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        pos = (/i,j,k/)
        call Eos_getData(IAXIS,pos,1,solnData,CENTER,eosData,massFraction)
        call Eos(MODE_DENS_EI,1,eosData,massFraction)
        call Eos_putData(IAXIS,pos,1,solnData,CENTER,eosData)
        ! Release pointers
        nullify(Vstar)

        !! HIGH side
        j = blkLimits(HIGH,JAXIS)
        ! Point to old/intermediate states
        Vstar => solnData(:,i,j,k)
        ! Point to the correct fluxes
        dFlux = hy_fly(:,i  ,j+1  ,k  )
        fac = 1.0
        dx = del(JAXIS)
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        pos = (/i,j,k/)
        call Eos_getData(IAXIS,pos,1,solnData,CENTER,eosData,massFraction)
        call Eos(MODE_DENS_EI,1,eosData,massFraction)
        call Eos_putData(IAXIS,pos,1,solnData,CENTER,eosData)
        ! Release pointers
        nullify(Vstar)
     enddo !j
  enddo !k
  !$omp end do
#if NDIM==3
  ! Correct KAXIS sides
  !$omp do schedule(static) collapse(2)
  do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
     do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
        !! LOW side
        k = blkLimits(LOW,KAXIS)
        ! Point to old/intermediate states
        Vstar => solnData(:,i,j,k)
        ! Point to the correct fluxes
        dFlux = -hy_flz(:,i  ,j  ,k  )
        fac = 1.0
        dx = del(KAXIS)
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        pos = (/i,j,k/)
        call Eos_getData(IAXIS,pos,1,solnData,CENTER,eosData,massFraction)
        call Eos(MODE_DENS_EI,1,eosData,massFraction)
        call Eos_putData(IAXIS,pos,1,solnData,CENTER,eosData)
        ! Release pointers
        nullify(Vstar)

        !! HIGH side
        k = blkLimits(HIGH,KAXIS)
        ! Point to old/intermediate states
        Vstar => solnData(:,i,j,k)
        ! Point to the correct fluxes
        dFlux = hy_flz(:,i  ,j  ,k+1  )
        fac = 1.0
        dx = del(KAXIS)
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        pos = (/i,j,k/)
        call Eos_getData(IAXIS,pos,1,solnData,CENTER,eosData,massFraction)
        call Eos(MODE_DENS_EI,1,eosData,massFraction)
        call Eos_putData(IAXIS,pos,1,solnData,CENTER,eosData)
        ! Release pointers
        nullify(Vstar)
     enddo !j
  enddo !k
  !$omp end do
#endif
#endif

  !$omp end parallel

  call Grid_releaseBlkPtr(blockID,solnData)

contains

  subroutine correctZone(Vstar,dFlux,dt,dx,fac)
    implicit none
    real, pointer, intent(INOUT) :: Vstar(:)
    real, intent(IN) :: dFlux(NFLUXES), dt, dx, fac
    real :: Ustar(NFLUXES)

    ! Construct vectors of conserved variables
    Ustar(HY_MASS)         = Vstar(DENS_VAR)
    Ustar(HY_XMOM:HY_ZMOM) = Vstar(DENS_VAR)*Vstar(VELX_VAR:VELZ_VAR)
    Ustar(HY_ENER)         = Vstar(DENS_VAR)*Vstar(ENER_VAR)
    Ustar(HY_NUM_FLUX+1:NFLUXES) = Vstar(SPECIES_BEGIN:MASS_SCALARS_END)*Vstar(DENS_VAR)
#ifdef SPARK_GLM
    Ustar(HY_FMGX:HY_FMGZ) = Vstar(MAGX_VAR:MAGZ_VAR)
    Ustar(HY_ENER) = Ustar(HY_ENER)+0.5*dot_product(Vstar(MAGX_VAR:MAGZ_VAR),&
         Vstar(MAGX_VAR:MAGZ_VAR))  ! * density ??? [KC]
    Ustar(HY_FPSI) = Vstar(PSIB_VAR)
#endif

    ! Now correct conserved vector with flux deltas
    ! The facP, facM definition is cracked. Need to know what side we are on
    Ustar = Ustar -dt/dx*(fac*dFlux)

    ! Update primitive variables
    emag = 0.0
#ifdef SPARK_GLM
    Vstar(MAGX_VAR:MAGZ_VAR) = Ustar(HY_FMGX:HY_FMGZ)
    ! Parabolic damping of PSI is applied to flux correction difference above
    Vstar(PSIB_VAR) = Ustar(HY_FPSI)
    emag = 0.5*dot_product(Vstar(MAGX_VAR:MAGZ_VAR),Vstar(MAGX_VAR:MAGZ_VAR))
    Vstar(MAGP_VAR) = emag
    Ustar(HY_ENER) = Ustar(HY_ENER) - emag
#endif
    Vstar(DENS_VAR)          = max(Ustar(HY_MASS),hy_smalldens)
    Vstar(VELX_VAR:VELZ_VAR) = Ustar(HY_XMOM:HY_ZMOM)/Vstar(DENS_VAR)
    Vstar(ENER_VAR)          = max(hy_smallE,Ustar(HY_ENER)/Vstar(DENS_VAR))

    ekin = .5*dot_product(Vstar(VELX_VAR:VELZ_VAR),Vstar(VELX_VAR:VELZ_VAR))
    Vstar(EINT_VAR) = max(hy_smallE,Vstar(ENER_VAR)-ekin)

    ! Divide partial densities by new mass densities to finalize
    ! update of new mass fractions.
    Vstar(SPECIES_BEGIN:MASS_SCALARS_END) = Ustar(HY_NUM_FLUX+1:NFLUXES)/Vstar(DENS_VAR)
  end subroutine correctZone

  subroutine  geoFacs(i,j,k,facM,facP)
    implicit none
    integer, intent(IN) :: i,j,k
    real, intent(OUT) :: facM, facP
    real    :: alpha, dx_sph

    if (hy_geometry == CARTESIAN) then
       facM = 1.0; facP = 1.0
       return
    endif
    facM = faceAreas(i  ,j,k)*dx/cellVolumes(i,j,k)
    facP = faceAreas(i+1,j,k)*dx/cellVolumes(i,j,k)

    if (xCenter(i) < 0.0) then
       facM = 0.
       facP = 0.
    end if
  end subroutine geoFacs

end subroutine hy_rk_correctFluxes
