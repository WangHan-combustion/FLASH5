!!****if* source/physics/Hydro/HydroMain/Spark/hy_rk_updateSoln
!!
!!  NAME
!!
!!  hy_rk_updateSoln
!!
!!  SYNOPSIS
!!
!!  call hy_rk_updateSoln ( integer(IN) :: blockID )
!!
!!  DESCRIPTION
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine hy_rk_updateSoln (blockID, dt, dtOld, limits, coeffs)

  use Hydro_data, ONLY : hy_threadWithinBlock, hy_starState, &
       hy_flx, hy_fly, hy_flz, hy_smallE, hy_smalldens, hy_geometry,&
       hy_grav, hy_4piGinv, hy_alphaGLM, hy_C_hyp
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_getDeltas, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getBlkData

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Spark.h"

  integer, intent(IN) :: blockID
  integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
  real, intent(IN) :: dt, dtOld
  real, dimension(3), intent(IN) :: coeffs

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(GRID_IHI_GC) :: xCenter, xLeft, xRight
  real, dimension(GRID_JHI_GC) :: yCenter
  real, dimension(GRID_KHI_GC) :: zCenter
  integer :: i,j,k,n,g

  real, pointer :: solnData(:,:,:,:)
  real, pointer :: V0(:), Vstar(:)

  real, dimension(NFLUXES) :: U0, Ustar

  real :: dx, dy, dz, del(MDIM)
  real, pointer :: Fm(:), Fp(:), Gm(:), Gp(:), Hm(:), Hp(:)

  real :: eint, ekin, emag
  ! Geometry factors
  real :: facM, facP
  real, dimension(NFLUXES) :: Sgeo, Sgrv, Stot
  real, allocatable, dimension(:,:,:) :: faceAreas
  real, allocatable, dimension(:,:,:) :: cellVolumes
  integer :: isize, jsize, ksize
  real :: dhdt

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  iSize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  jSize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  kSize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  call Grid_getBlkPtr(blockID,solnData)

  call Grid_getDeltas(blockID,del)
  dx = del(IAXIS); dy = del(JAXIS); dz = del(KAXIS)
  dhdt = minval(del(1:NDIM))/(coeffs(3)*dt)

  if (hy_geometry /= CARTESIAN) then
     allocate(faceAreas(iSize,jSize,kSize))
     call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, &
          (/1,1,1/), faceAreas, (/isize, jsize, ksize/) )

     allocate(cellVolumes(iSize, jSize, kSize))
     call Grid_getBlkData(blockID, CELL_VOLUME, 0, EXTERIOR, &
          (/1,1,1/), cellVolumes, (/isize, jsize, ksize/) )

     call Grid_getCellCoords(IAXIS,blockID, CENTER,.true.,xCenter, blkLimitsGC(HIGH,IAXIS))
     call Grid_getCellCoords(IAXIS,blockID, LEFT_EDGE,.true.,xLeft, blkLimitsGC(HIGH,IAXIS))
     call Grid_getCellCoords(IAXIS,blockID, RIGHT_EDGE,.true.,xRight, blkLimitsGC(HIGH,IAXIS))
     call Grid_getCellCoords(JAXIS,blockID, CENTER,.true.,yCenter, blkLimitsGC(HIGH,JAXIS))
     call Grid_getCellCoords(KAXIS,blockID, CENTER,.true.,zCenter, blkLimitsGC(HIGH,KAXIS))
  endif

  !$omp parallel if (hy_threadWithinBlock) &
  !$omp default(none) &
  !$omp private(n,i,j,k,Fm,Fp,Gm,Gp,Hm,Hp,V0,Vstar,U0,Ustar,Sgeo,facM,facP,&
  !$omp         Sgrv,Stot,emag,ekin)&
  !$omp shared(dt,solnData,limits,hy_starState,hy_flx,hy_fly,hy_flz,hy_grav,&
  !$omp        dtOld,xCenter,xLeft,xRight,yCenter,zCenter,blockID,&
  !$omp        blkLimits,blkLimitsGC,coeffs,dx,dy,dz,hy_alphaGLM, hy_C_hyp,&
  !$omp        dhdt, hy_smalldens, hy_smallE)

  !  Begin loop over zones
  !$omp do schedule(guided) collapse(3)
  do k = limits(LOW,KAXIS), limits(HIGH,KAXIS)
     do j = limits(LOW,JAXIS), limits(HIGH,JAXIS)
        do i = limits(LOW,IAXIS), limits(HIGH,IAXIS)

           ! Point to old/intermediate states
           V0    => solnData(:,i,j,k)
           Vstar => hy_starState(:,i,j,k)

           ! Point to the correct fluxes
           Fm => hy_flx(:,i  ,j  ,k  )
           Fp => hy_flx(:,i+1,j  ,k  )
           Gm => hy_fly(:,i  ,j  ,k  )
           Gp => hy_fly(:,i  ,j+1,k  )
           Hm => hy_flz(:,i  ,j  ,k  )
           Hp => hy_flz(:,i  ,j  ,k+1)

           ! Construct vectors of conserved variables
           U0(HY_MASS)            = V0(DENS_VAR)
           U0(HY_XMOM:HY_ZMOM)    = V0(DENS_VAR)*V0(VELX_VAR:VELZ_VAR)
           U0(HY_ENER)            = V0(DENS_VAR)*V0(ENER_VAR)
           U0(HY_NUM_FLUX+1:NFLUXES) = V0(SPECIES_BEGIN:MASS_SCALARS_END)*V0(DENS_VAR)
           Ustar(HY_MASS)         = Vstar(DENS_VAR)
           Ustar(HY_XMOM:HY_ZMOM) = Vstar(DENS_VAR)*Vstar(VELX_VAR:VELZ_VAR)
           Ustar(HY_ENER)         = Vstar(DENS_VAR)*Vstar(ENER_VAR)
           Ustar(HY_NUM_FLUX+1:NFLUXES) = Vstar(SPECIES_BEGIN:MASS_SCALARS_END)*Vstar(DENS_VAR)
#ifdef SPARK_GLM
           U0(HY_FMGX:HY_FMGZ)    = V0(MAGX_VAR:MAGZ_VAR)
           U0(HY_ENER) = U0(HY_ENER)+0.5*dot_product(V0(MAGX_VAR:MAGZ_VAR),&
                V0(MAGX_VAR:MAGZ_VAR))
           U0(HY_FPSI) = V0(PSIB_VAR)
           Ustar(HY_FMGX:HY_FMGZ) = Vstar(MAGX_VAR:MAGZ_VAR)
           Ustar(HY_ENER) = Ustar(HY_ENER)+0.5*dot_product(Vstar(MAGX_VAR:MAGZ_VAR),&
                Vstar(MAGX_VAR:MAGZ_VAR))
           Ustar(HY_FPSI) = Vstar(PSIB_VAR)
#endif

           ! Get gemoetric factors and sources
           call geoFacs(i,j,k,facM,facP,Sgeo,Ustar,Vstar)

           ! Get gravitational source terms
           call gravSources(Ustar,hy_grav(:,i,j,k),Sgrv)

           ! Sum total source terms
           Stot = Sgeo + Sgrv

           ! Now update conserved vector with flux gradients
           Ustar = coeffs(1)*U0 + coeffs(2)*Ustar +coeffs(3)*( &
                -dt/dx*(facP*Fp-facM*Fm) &
#if NDIM > 1
                -dt/dy*(Gp-Gm) &
#if NDIM ==3
                -dt/dz*(Hp-Hm) &
#endif
#endif
                +dt*Stot)

           ! Update primitive variables
           emag = 0.0
#ifdef SPARK_GLM
           !print *, coeffs(1)*U0(HY_FPSI), coeffs(2)*Ustar(HY_FPSI), Fp(HY_FPSI)-Fm(HY_FPSI)
           Vstar(MAGX_VAR:MAGZ_VAR) = Ustar(HY_FMGX:HY_FMGZ)
           Vstar(PSIB_VAR) = Ustar(HY_FPSI)*exp(-hy_alphaGLM*hy_C_hyp/dhdt)
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

#ifdef GPOT_VAR
           ! Now extrapolate gravity to intermediate time state
           ! the star state GPOT_VAR will be reset so that GPOL_VAR isn't screwed up
           Vstar(GPOT_VAR) = coeffs(1)*V0(GPOT_VAR) + coeffs(2)*Vstar(GPOT_VAR) &
                + coeffs(3)*dt*(V0(GPOT_VAR) - V0(GPOL_VAR))/dtOld
#endif
           ! Release pointers
           nullify(V0)
           nullify(Vstar)
           nullify(Fm)
           nullify(Fp)
           nullify(Gm)
           nullify(Gp)
           nullify(Hm)
           nullify(Hp)

        enddo !i
     enddo !j
  enddo !k
  !$omp end do nowait
  !$omp end parallel

  call Grid_releaseBlkPtr(blockID,solnData)
  if (hy_geometry /= CARTESIAN) then
     deallocate(faceAreas)
     deallocate(cellVolumes)
  end if

contains

  subroutine  geoFacs(i,j,k,facM,facP,Sgeo,U,V)
    implicit none
    integer, intent(IN) :: i,j,k
    real, intent(OUT) :: facM, facP
    real, dimension(NFLUXES) :: Sgeo
    real, intent(IN) :: U(:)
    real, pointer, intent(IN) :: V(:)

    real    :: presStar, densStar, pmomStar, tmomStar, xmomStar
    real    :: pmagStar, xmagStar, zmagStar
    integer :: VEL_PHI, MOM_PHI, MOM_PHI_FLUX, MAG_PHI,  MAG_PHI_FLUX
    integer :: VEL_ZI, MOM_ZI, MOM_ZI_FLUX, MAG_ZI,  MAG_ZI_FLUX
    integer :: VEL_THT, MOM_THT, MOM_THT_FLUX
    real    :: alpha, dx_sph

    if (hy_geometry == CARTESIAN) then
       facM = 1.0; facP = 1.0; Sgeo = 0.0
       return
    endif

    select case(hy_geometry) ! First, select whether y or z is phi-direction
    case(CYLINDRICAL)
       MOM_PHI = HY_ZMOM
       MOM_PHI_FLUX = HY_ZMOM
       MOM_ZI       = HY_YMOM
       MOM_ZI_FLUX  = HY_YMOM
#ifdef SPARK_GLM
       MAG_PHI      = HY_MAGZ
#endif

       alpha = 1.
    case(POLAR)
       MOM_PHI      = HY_YMOM
       MOM_PHI_FLUX = HY_YMOM
       MOM_ZI       = HY_ZMOM
       MOM_ZI_FLUX  = HY_ZMOM
#ifdef SPARK_GLM
       MAG_PHI      = HY_MAGY
#endif
       alpha = 1.
    case(SPHERICAL)
       MOM_PHI      = HY_ZMOM
       MOM_PHI_FLUX = HY_ZMOM
       MOM_THT      = HY_YMOM
       MOM_THT_FLUX = HY_YMOM
       dx_sph = (xRight(i)**3 - xLeft(i)**3) / (3.*xCenter(i)**2)
       alpha  = 2.
    end select

    facM = faceAreas(i  ,j,k)*dx/cellVolumes(i,j,k)
    facP = faceAreas(i+1,j,k)*dx/cellVolumes(i,j,k)

    Sgeo = 0.
    !! Calculate geometrical source terms.  See S&O 75.
    Sgeo(HY_XMOM) = (V(DENS_VAR)*V(VELZ_VAR)*V(VELZ_VAR) + alpha*V(PRES_VAR)) / xCenter(i)!T phi,phi
    Sgeo(MOM_PHI) =  V(DENS_VAR)*V(VELZ_VAR)*V(VELX_VAR) / xCenter(i)!T phi,r

#ifdef SPARK_GLM
    ! P* is the total Pressure
    ! This presently does not work for POLAR coordinates
    Sgeo(HY_XMOM) = Sgeo(HY_XMOM) - (V(MAGZ_VAR)**2 - alpha*V(MAGP_VAR))/ xCenter(i)
    Sgeo(MOM_PHI) = Sgeo(MOM_PHI) - V(MAGZ_VAR)*V(MAGX_VAR) / xCenter(i)

    Sgeo(MAG_PHI) = - (V(VELZ_VAR)*V(MAGX_VAR) - V(MAGZ_VAR)*V(VELX_VAR)) / xCenter(i) !O phi,r
#endif
    Sgeo(MOM_PHI) = - Sgeo(MOM_PHI)

    if (hy_geometry == SPHERICAL) then
       Sgeo(HY_XMOM) = Sgeo(HY_XMOM) + U(MOM_THT)**2/V(DENS_VAR) / xCenter(i)
       Sgeo = Sgeo*dx/dx_sph
    endif
    if (xCenter(i) < 0.0) then
       facM = 0.
       facP = 0.
       Sgeo = 0.
    end if
  end subroutine geoFacs

  subroutine gravSources(U,g,S)
    implicit none
    real, intent(IN) :: U(NFLUXES)
    real, intent(IN) :: g(MDIM)
    real, intent(OUT) :: S(NFLUXES)
    S = 0.
#ifdef GRAVITY
    S(HY_XMOM:HY_ZMOM) = U(HY_MASS)*g(:)
    S(HY_ENER) = dot_product(U(HY_XMOM:HY_ZMOM),g(:))
#endif
  end subroutine gravSources

end subroutine hy_rk_updateSoln
