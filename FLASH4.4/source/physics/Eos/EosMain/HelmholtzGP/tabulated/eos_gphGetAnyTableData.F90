!!****if* source/physics/Eos/EosMain/Tabulated/eos_gphGetSpeciesAnyTableData
!!
!! NAME
!!
!!  eos_gphGetAnyTableData
!!
!! SYNOPSIS
!!
!!  call eos_gphGetAnyTableData(integer (in)  :: species,
!!                                    real    (in)  :: speciesTemperature,
!!                                    real    (in)  :: speciesDensity,
!!                                    integer (in)  :: maxComp,
!!                                    integer (in)  :: selTT,
!!                                    integer (in)  :: needDerivs(:,:),
!!                                    real    (out) :: resultTT(:,:))
!!
!! DESCRIPTION
!!
!!  Extracts via interpolation an average ionization, internal energy, or heat
!!  capacity value for all components for
!!  the specified temperature, density, and species. The interpolation method is
!!  bilinear.
!!
!! ARGUMENTS
!!
!!   species            : The species index
!!   speciesTemperature : The species temperature
!!   speciesDensity     : The species density
!!   maxComp :          : The highest component (ions, electrons, or matter) for which
!!                        data are to be returned
!!   selTT              : Select ZF, EN, or HC table type
!!   resultTT          : The value of the determined average ionization etc.
!!
!!***

#include "constants.h"
#include "Eos.h"

subroutine eos_gphGetAnyTableData (xpos,            &
                                        wanted, &
                                        selTT,              &
                                        derDefsIn,         &
                                        needDerivs,         &
                                        ipos,         &
                                        tausIn, &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        resultTT, &
                                        cloIn, chiIn)

  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_data,         ONLY : eos_meshMe
  use eos_gphData,      ONLY : one,ten,                    &
                               EOS_TAB_NCOMP,              &
                               EOS_TABINT_DERIV_0,         &
                               EOS_TABINT_DERIV_DT,        &
                               EOS_TABINT_DERIV_DD,        &
!                               EOS_GPH_NALLTAB,            &
                               EOS_TABVT_EN,               &
                               EOS_TABVT_HC,               &
                               EOS_TABVT_PR,               &
                               EOS_TABVT_ENTR,             &
                               DENS_IVAR, TEMP_IVAR,&
                                eos_useLogTables,          &
                                eosT_varTableGroupT,      &
                                eos_gphTheTable

  
  implicit none

  real,    intent (in)  :: xpos(EOST_MAX_IVARS)
  logical,            intent (in) :: wanted(EOS_TAB_NCOMP)
  integer, intent (in)  :: selTT
  integer, intent (in)  :: derDefsIn(:,0:)
  integer, intent (in)  :: needDerivs(:,:)
  real,    intent (out) :: resultTT(0:,:)

  type(eosT_varTableGroupT),pointer :: thisTypeTable
  integer, pointer      :: derDefs(:,:)

  logical,intent(IN) :: lowerBoundary(1:EOST_MAX_IVARS)
  logical :: upperBoundary(1:EOST_MAX_IVARS)
  logical,intent(IN) :: withinBoundary(1:EOST_MAX_IVARS)
  logical :: found

  integer :: d,t
  integer :: ivi,      iw,ix,iy,iz
  integer :: tableDim, nw,nx,ny,nz
  integer :: ic, ib, der, derW, wrtVar, degW
  integer,intent(IN) :: ipos(EOST_MAX_IVARS)
  integer :: varType
  integer :: ncorners
  integer :: nderiW
  integer :: nderivs
  integer :: nstepsDensity
  integer :: nstepsTemperature
  integer :: g

  real,           dimension(EOST_MAX_IVARS,0:1) :: clohi
!!$  real,           dimension(EOST_MAX_IVARS) :: clo, chi
  real,intent(IN),dimension(EOST_MAX_IVARS) :: cloIn, chiIn
  real :: f1,f2,f3,f4
  real :: o1,o2,o3,o4
  real,allocatable :: o(:,:)
  real,allocatable :: f(:,:,:)
  real :: r
  real :: expo1(0:1, EOST_MAX_IVARS)
  real :: speciesDensityTT
  real :: speciesTemperatureTT
  real,intent(IN) :: tausIn(EOST_MAX_IVARS)
  real            :: taus(EOST_MAX_IVARS)
  real            :: ellFact, ellReno, fact

  real,parameter :: largeExponent = 77.0      ! was 0.25 * log10(HUGE(largeExponent))
  real,parameter :: largeMultiplier = 10.0**77  ! was 10.0**largeExponent

  intrinsic ibits
!
!
!   ...Set the handle to the ZF, EN, PR, or HC tables and associated data arrays.
!
!
  select case (selTT)
#if(0)
  case(EOS_TABULAR_Z)
     varType = EOS_TABVT_ZF
     thisTypeTable => eos_allTab(species)%tg(varType)
  case(EOS_TABULAR_E)
     varType = EOS_TABVT_EN
     thisTypeTable => eos_allTab(species)%tg(varType)
  case(EOS_TABULAR_C)
     varType = EOS_TABVT_HC
     thisTypeTable => eos_allTab(species)%tg(varType)
  case(EOS_TABULAR_P)
     varType = EOS_TABVT_PR
     thisTypeTable => eos_allTab(species)%tg(varType)
#endif
  case(EOS_TABULAR_S)
     varType = EOS_TABVT_ENTR
     thisTypeTable => eos_gphTheTable%tg(varType)
  case default
     varType = 0
  end select

  if (varType == 0) then
      call Driver_abortFlash ('[eos_gphGetAnyTableData] ERROR: no handle to tables')
  end if
!
!
!   ...Get the current temperature and ion number density of the species and find:
!
!        1) the temperature/density quadrant (T1,T2,D1,D2) boundary containing
!           containing the species's temperature and density (x)
!
!        2) the table index quadrant (i,j,k,l) of the boundary, as indicated
!           in the figure below.
!
!        3) the four tabulated values (o1,o2,o3,o4)
!
!
!                    o3------------o4   T2 (j)
!                     |            |
!                     |            |
!                     |        x   |
!                     |            |
!                     |            |
!                    o1 -----------o2   T1 (i)
!
!                  D1 (k)         D2 (l)
!
!
!      In case the temperature and/or density of the species lay outside
!      the tabulated boundaries, take the corresponding boundary values
!      of the tables. The criteria as to when the species's temperature and
!      density belong within the boundary are:
!
!                         T1 =< speciesTemperature =< T2
!                         D1 =<   speciesDensity   =< D2
!
!                         clo(TEMP_IVAR) =< speciesTemperature =< chi(TEMP_IVAR)
!                         clo(DENS_IVAR) =<   speciesDensity   =< chi(DENS_IVAR)
!
!


  taus = tausIn
!!$  clo = cloIn; chi = chiIn
  clohi(:,0) = cloIn
  clohi(:,1) = chiIn

  upperBoundary = .NOT. (lowerBoundary .OR. withinBoundary)


#ifdef DEBUG_EOS
  if (selTT==EOS_TABULAR_S) then
     print*,'selTT,wanted:',selTT,wanted
  end if
#endif
  do g = 1,EOS_TAB_NCOMP
     if (.NOT.wanted(g)) cycle
     if (.NOT.associated(thisTypeTable%table(g)%table)) then
        if (eos_meshMe==MASTER_PE) then
           print*,'Table not associated:g,selTT',g,selTT
        end if
        cycle
     end if
!!$     o1 = thisTypeTable%table(g)%table(i,k)
!!$     o2 = thisTypeTable%table(g)%table(i,l)
!!$     o3 = thisTypeTable%table(g)%table(j,k)
!!$     o4 = thisTypeTable%table(g)%table(j,l)

     print*,'gphGetAnyTableData: sought ipos is',ipos
     ncorners = thisTypeTable % td % ncorners
     nderivs  = thisTypeTable % td % nderivs
     allocate(o(0:ncorners-1, 0:nderivs))
     associate (d => thisTypeTable % table(g) % table)
       do der=0,nderivs
          do ic=0,ncorners-1
             o(ic, der) = d (ic, der, ipos(1), ipos(2), ipos(3), ipos(4) )
             print*,'ic,der, o:',ic,der,o(ic, der)
          end do
       end do
     End associate

#if(0)
     if (.NOT. withinBoundary( ivi ) .AND. &
          (selTT==EOS_TABULAR_E .OR. selTT==EOS_TABULAR_P) ) then
        if (eos_useLogTables) then
           speciesTemperatureTT = log10 (speciesTemperature)
        else
           speciesTemperatureTT = speciesTemperature
        end if
        if ( (lowerBoundary( ivi ) .AND. (eos_useLogTables .OR. chi(TEMP_IVAR) > 0.0)) ) then
           if (eos_useLogTables) then
              clo(TEMP_IVAR) = chi(TEMP_IVAR) - largeExponent
              tau = (speciesTemperatureTT - clo(TEMP_IVAR)) / largeExponent
              o1 = o3 - largeExponent
              o2 = o4 - largeExponent
           else
              clo(TEMP_IVAR) = 0.0
              tau = speciesTemperatureTT / chi(TEMP_IVAR)
              o1 = 0.0
              o2 = 0.0
           end if
        else if ( (upperBoundary( ivi ) .AND. (eos_useLogTables .OR. chi(TEMP_IVAR) > 0.0)) ) then
           if (eos_useLogTables) then
              chi(TEMP_IVAR) = clo(TEMP_IVAR) + largeExponent
              tau = (speciesTemperatureTT - clo(TEMP_IVAR)) / largeExponent
              o3 = o1 + largeExponent
              o4 = o2 + largeExponent
           else
              chi(TEMP_IVAR) = clo(TEMP_IVAR) * largeMultiplier
              tau = (speciesTemperatureTT - clo(TEMP_IVAR)) / (chi(TEMP_IVAR) - clo(TEMP_IVAR))
              o3 = o1 * largeMultiplier
              o4 = o2 * largeMultiplier
           end if
        end if
     end if
     if (.NOT. withinBoundary( ivi ) .AND. &
          (selTT==EOS_TABULAR_P) ) then
        if (eos_useLogTables) then
           speciesDensityTT     = log10 (speciesDensity)
        else
           speciesDensityTT     = speciesDensity
        end if
        if ( (lowerBoundary( ivi ) .AND. (eos_useLogTables .OR. chi(DENS_IVAR) > 0.0)) ) then
           if (eos_useLogTables) then
              clo(DENS_IVAR) = chi(DENS_IVAR) - largeExponent
              delta = (speciesDensityTT - clo(DENS_IVAR)) / largeExponent
              o1 = o2 - largeExponent
              o3 = o4 - largeExponent
           else
              clo(DENS_IVAR) = 0.0
              delta = speciesDensityTT / chi(DENS_IVAR)
              o1 = 0.0
              o3 = 0.0
           end if
        else if ( (upperBoundary( ivi ) .AND. (eos_useLogTables .OR. chi(DENS_IVAR) > 0.0)) ) then
           if (eos_useLogTables) then
              chi(DENS_IVAR) = clo(DENS_IVAR) + largeExponent
              delta = (speciesDensityTT - clo(DENS_IVAR)) / largeExponent
              o2 = o1 + largeExponent
              o4 = o3 + largeExponent
           else
              chi(DENS_IVAR) = clo(DENS_IVAR) * largeMultiplier
              delta = (speciesDensityTT - clo(DENS_IVAR)) / (chi(DENS_IVAR) - clo(DENS_IVAR))
              o2 = o1 * largeMultiplier
              o4 = o3 * largeMultiplier
           end if
        end if
     end if
#endif
     !
     !
     !   ...Do the bilinear interpolation:
     !
     !                result =   o1 * [(1-tau)*(1-delta)]
     !                         + o2 * [delta*(1-tau)]
     !                         + o3 * [tau*(1-delta)]
     !                         + o4 * [delta*tau]
     !
     !
!!$     tau   = (speciesTemperatureTT - T1) / (T2 - T1)
!!$     delta = (speciesDensityTT     - D1) / (D2 - D1)

!!$     tau   = (speciesTemperatureTT - clo(TEMP_IVAR)) / (chi(TEMP_IVAR) - clo(TEMP_IVAR))
!!$     delta = (speciesDensityTT     - clo(DENS_IVAR)) / (chi(DENS_IVAR) - clo(DENS_IVAR))

!!$     f1 = (one - tau) * (one - delta)
!!$     f2 = delta * (one - tau)
!!$     f3 = tau * (one - delta)
!!$     f4 = delta * tau


     associate (tableDim => thisTypeTable % td % n, &
          ell => thisTypeTable % table(g) &
          % ells(:,ipos(1),ipos(2),ipos(3),ipos(4)) )
     nw = 1; nx = 1; ny = 1; nz = 1
     if (tableDim < 1) nw = 0
     if (tableDim < 2) nx = 0
     if (tableDim < 3) ny = 0
     if (tableDim < 4) nz = 0

     print*,'ells:',thisTypeTable % table(g) % ells
     print*,'ell :',ell
     expo1(:,:) = 1.0
     do ivi = 1, tableDim
        ellReno = ell(ivi) / (clohi(ivi,1) - clohi(ivi,0))
        ellFact = 1.0 / ellReno
        print*,'ellerie:',ell(ivi) , (clohi(ivi,1) - clohi(ivi,0)),ellReno,ellFact
        expo1(0, ivi) = exp(-0.5 * (ellFact *      taus(ivi) )**2)
        expo1(1, ivi) = exp(-0.5 * (ellFact * (1.0-taus(ivi)))**2)
        print*,'ivi, expo1(0:1,ivi):',ivi, expo1(0:1,ivi)
     end do

   derDefs => thisTypeTable % table(g) % derDefs
   print*,'lbound(needDerivs):',lbound(needDerivs)
   print*,'ubound(needDerivs):',ubound(needDerivs)
   print*,'lbound( derDefs  ):',lbound( derDefs  )
   print*,'ubound( derDefs  ):',ubound( derDefs  )
     nderiW = size(needDerivs,2)
     allocate(f(0:ncorners-1,0:nderivs,0:nderiW))
!!$     allocate(result(0:nderiW))

     do iz=0,nz
        do iy=0,ny
           do ix=0,nx
              do iw=0,nw
                 ic = iw + 2 * (ix + 2 * (iy + 2*iz))
                 f(ic,0,0) = expo1(iw,1) * expo1(ix,2) * expo1(iy,3) * expo1(iz,4)
 
                 do derW=0,nderiW
                    do der=0,nderivs
                       f(ic,der,derW) = f(ic,0,0)
                       if (derW+der == 0) CYCLE
                       do wrtVar = 1,tableDim
                          associate ( degree => derDefs   (wrtVar,der ) )
                            degW = 0
                            if (derW > 0) degW = needDerivs(wrtVar,derW)

                            if (degW+degree > 0) then
                               ib = ibits(ic, wrtVar-1, 1)
                               select case (degW+degree)
                               case(1) ! Cf. "probsbilists'" Hermite polynomials...
                                  fact = (xpos(wrtVar) - clohi(wrtVar,ib)) / ell(wrtVar)**2
                                  fact = fact * (-1)**degW
!!$                                  print*,'...derW,der,wrtVar,xpos(wrtVar),ib,clohi(wrtVar,ib):', &
!!$                                       derW,der,wrtVar,xpos(wrtVar),ib,clohi(wrtVar,ib),'->',fact
                               case(2)
!!$                               fact = (chi(wrtVar) - clo(wrtVar)) / ell(wrtVar)**2
!!$                               fact = fact**2 - 1. / ell(wrtVar)**2
                                  fact = ((xpos(wrtVar) - clohi(wrtVar,ib))**2 - ell(wrtVar)**2) / &
                                                                                 ell(wrtVar)**4
                                  fact = fact * (-1)**degW
                               case(3)
                                  fact = ((xpos(wrtVar) - clohi(wrtVar,ib))**3                     &
                                           - 3*(xpos(wrtVar) - clohi(wrtVar,ib))*ell(wrtVar)**2) / &
                                                                                 ell(wrtVar)**6
                                  fact = -fact * (-1)**degW
                               case(4)
                                  fact = ((xpos(wrtVar) - clohi(wrtVar,ib))**4                     &
                                           - 6*(xpos(wrtVar) - clohi(wrtVar,ib))**2*ell(wrtVar)**2 &
                                           + 3                                     *ell(wrtVar)**4) / &
                                                                                 ell(wrtVar)**8
                                  fact = fact * (-1)**degW
                               case default
                                  print*,'derW,der,wrtVar,degW,degree,fact,derDefs:',needDerivs(:,derW),derDefs(:,der)
                                  print*, derW,der,wrtVar,degW,degree,fact,derDefs
                                  call Driver_abortFlash("derivative degree is too high!")
                               end select
                               f(ic,der,derW) = f(ic,der,derW) * fact
                            end if
                          end associate
                       end do
                    end do
                 end do

             end do
           end do
        end do
     end do

   end associate


   print*,'lbound(o):',lbound(o)
   print*,'ubound(o):',ubound(o)
   print*,'lbound(f):',lbound(f)
   print*,'ubound(f):',ubound(f)
   
   do derW=0,nderiW
      r = 0.0
      do der=0,nderivs
         do ic=0,ncorners-1
            r = r + o(ic,der) * f(ic,der,derW) 
!!$            print*,'...just added o*f; ic,der,derW,g=',ic,der,derW,g,', o,f,r', o(ic,der),f(ic,der,derW),r
         end do
      end do
      resultTT(derW, g) = r !DOT_PRODUCT(o,f)
   end do

!!$     resultTT(EOS_TABINT_DERIV_0,g) = 0 !DOT_PRODUCT(o,f)
     !   ...Convert logarithmic form to real form (if needed).
     if (eos_useLogTables) then
        resultTT(EOS_TABINT_DERIV_0,g) = ten ** resultTT(EOS_TABINT_DERIV_0,g)
     end if

#if 0
     if (needDerivs .GE. EOS_TABINT_DERIV_DT) then
        if (.NOT.withinBoundary( ivi ) .AND. &
            .NOT.( (selTT==EOS_TABULAR_E .OR. selTT==EOS_TABULAR_P) .AND. &
                   (  (lowerBoundary( ivi ) .AND. (eos_useLogTables .OR. chi(TEMP_IVAR) > 0.0) ) .OR. &
                      (upperBoundary( ivi ) .AND. (eos_useLogTables .OR. chi(TEMP_IVAR) > 0.0) ) ) )) then
           resultTT(EOS_TABINT_DERIV_DT,g) = 0.0
        else
           f1 =   delta - one
           f2 = - delta
           f3 =   one - delta
           f4 =   delta
           resultTT(EOS_TABINT_DERIV_DT,g) = &
                (o1 * f1 + o2 * f2 + o3 * f3 + o4 * f4) / (chi(TEMP_IVAR) - clo(TEMP_IVAR))
           !   ...Convert logarithmic form to real form (if needed).
           if (eos_useLogTables) then
              resultTT(EOS_TABINT_DERIV_DT,g) = &
                   resultTT(EOS_TABINT_DERIV_0,g) * resultTT(EOS_TABINT_DERIV_DT,g) &
                   / speciesTemperature
           end if
        end if
     end if

     if (needDerivs .GE. EOS_TABINT_DERIV_DD) then
        if (.NOT.withinBoundary( ivi ) .AND. &
            .NOT.( (selTT==EOS_TABULAR_P) .AND. &
                   (  (lowerBoundary( ivi ) .AND. (eos_useLogTables .OR. chi(DENS_IVAR) > 0.0) ) .OR. &
                      (upperBoundary( ivi ) .AND. (eos_useLogTables .OR. chi(DENS_IVAR) > 0.0) ) ) )) then
           resultTT(EOS_TABINT_DERIV_DD,g) = 0.0
        else
           f1 =   tau - one
           f2 =   one - tau
           f3 = - tau
           f4 =   tau
           resultTT(EOS_TABINT_DERIV_DD,g) = &
                (o1 * f1 + o2 * f2 + o3 * f3 + o4 * f4) / (chi(DENS_IVAR) - clo(DENS_IVAR))
           !   ...Convert logarithmic form to real form (if needed).
           if (eos_useLogTables) then
              resultTT(EOS_TABINT_DERIV_DD,g) = &
                   resultTT(EOS_TABINT_DERIV_0,g) * resultTT(EOS_TABINT_DERIV_DD,g) &
                   / speciesDensity
           end if
        end if
     end if
#endif
  end do                        ! do g
  print*,'At ',xpos,', returning',resultTT
!
!
!
!


!
!
!   ...Ready! 
!
!
  return
end subroutine eos_gphGetAnyTableData
