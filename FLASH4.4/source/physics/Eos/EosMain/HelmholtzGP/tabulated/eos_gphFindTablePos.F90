!!****if* source/physics/Eos/EosMain/Tabulated/eos_gphFindTablePos
!!
!! NAME
!!
!!  eos_gphFindTablePos
!!
!! SYNOPSIS
!!
!!  call eos_gphFindTablePos(integer (real    (in)  :: xpos(EOST_MAX_IVARS),
!!                                    integer (in)  :: maxComp,
!!                                    integer (in)  :: selTT)
!!
!! DESCRIPTION
!!
!!  Given desired Temperature and Density values and a ((temp(i), dens(k),
!!  i=1,ntemp, k=1,ndens) description of a (temperature, density) table grid,
!!  find and return the index pair (i,k) such that
!!        temp(i) < Temperature <=  temp(j)
!!  and
!!        dens(k) < Density     <=  dens(l)  .
!!
!! ARGUMENTS
!!
!!   species            : The species index
!!   xpos(ivi) : The species temperature
!!   speciesDensity     : The species density
!!   maxComp :          : The highest component (ions, electrons, or matter) for which
!!                        data are to be returned
!!   selTT              : Select ZF, EN, or HC table type
!!
!!***

#include "Flash.h"
#include "Eos.h"

subroutine eos_gphFindTablePos (        xpos, &
                                        selTT,              &
                                        i,j,k,l,&
                                        taus, &
                                        iPrev, kPrev, &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        T1,T2,D1,D2)

  use Driver_interface,  ONLY : Driver_abortFlash
  use eos_gphData,      ONLY : EOS_TABVT_ENTR,             &
                                eos_useLogTables,          &
                                eosT_varTableGroupT,      &
                                eos_gphTheTable


  implicit none

  real,    intent (in)  :: xpos(EOST_MAX_IVARS)
  integer, intent (in)  :: selTT

  type(eosT_varTableGroupT),pointer :: thisTypeTable
  real,    pointer :: thisTypeIvar (:)

#define DENS_IVAR 1
#define TEMP_IVAR 2

  logical,intent(OUT) :: lowerBoundary(1:EOST_MAX_IVARS)
  logical :: upperBoundary(1:EOST_MAX_IVARS)
  logical,intent(OUT) :: withinBoundary(1:EOST_MAX_IVARS)
  logical :: found
  integer :: ivi                ! input var index

  integer :: d,t
  integer,intent(OUT) :: i,j,k,l
  integer,intent(INOUT) :: iPrev,kPrev
  integer :: varType
  integer :: nsteps
  integer :: g

  real,intent(OUT) :: D1,D2,T1,T2
  real :: f1,f2,f3,f4
  real :: o1,o2,o3,o4
  real :: xi
  real,intent(out) :: taus(EOST_MAX_IVARS)

!
!
!   ...Set the handle to the ZF, EN, or HC tables and associated data arrays.
!
!
  select case (selTT)
#if(0)
  case(EOS_TABULAR_Z)
     varType = EOS_TABVT_ZF
  case(EOS_TABULAR_E)
     varType = EOS_TABVT_EN
  case(EOS_TABULAR_C)
     varType = EOS_TABVT_HC
  case(EOS_TABULAR_P)
     varType = EOS_TABVT_PR
#endif
  case(EOS_TABULAR_S)
     varType = EOS_TABVT_ENTR
  case default
     varType = 0
  end select

  if (varType == 0) then
      call Driver_abortFlash ('[eos_gphFindTablePos] ERROR: no handle to tables')
  end if

  thisTypeTable => eos_gphTheTable%tg(varType)

 Do ivi = 1,tableDim
  thisTypeIvar => thisTypeTable%td % c(i) % val
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
!      Generic:
!                         T1 =< xpos(ivi) =< T2
!                         D1 =<   pos(ivi)   =< D2
!
!
!!$  if (eos_useLogTables) then
!!$      xi = log10 (xpos(ivi))
!!$      xi     = log10 (pos(ivi))
!!$  else
      xi = xpos(ivi)
      xi     = pos(ivi)
!!$  end if

  nsteps = thisTypeTable%td % c(i) % nIval + 1

  lowerBoundary( ivi ) = xi < thisTypeIvar (1                )
  if (lowerBoundary( ivi )) then
     upperBoundary( ivi ) = .FALSE.
  else
     upperBoundary( ivi ) = xi > thisTypeIvar (nsteps) .OR. nsteps == 1
  end if
  withinBoundary( ivi ) = (.not.lowerBoundary( ivi )) .and. (.not.upperBoundary( ivi ))

  if (withinBoundary( ivi )) then
     if (iPrev > 0 .AND. iPrev < nsteps) then
        i = iPrev
        j = i + 1

        T1 = thisTypeIvar (i)
        T2 = thisTypeIvar (j)

        found = ((T1-xi)*(T2-xi) .LE. 0.0)
     else
        found = .FALSE.
     end if
     if (.NOT.found) then
        do t = 2,nsteps
           if (thisTypeIvar (t) >= xi) then
              i  = t - 1
              j  = t
              T1 = thisTypeIvar (i)
              T2 = thisTypeIvar (j)
#ifndef FLASH_OPENMP              
              iPrev = i
#endif
              exit
           end if
        end do
     end if
  else
      if (lowerBoundary( ivi )) then
          i  = 1
          j  = 1
          T1 = xi
          T2 = thisTypeIvar (1)
      end if

      if (upperBoundary( ivi )) then
          i  = nsteps
          j  = nsteps
          T1 = thisTypeIvar (nsteps)
          T2 = xi
      end if
  end if



  nsteps = thisTypeTable%td % c(i) % nIval + 1

  ! Check to see if density is off of the table:
  lowerBoundary( ivi ) = xi < thisTypeIvar (1            )
  if (lowerBoundary( ivi )) then
     upperBoundary( ivi ) = .FALSE.
  else
     upperBoundary( ivi ) = xi > thisTypeIvar (nsteps) .OR. nsteps == 1
  end if
  withinBoundary( ivi ) = (.not.lowerBoundary( ivi )) .and. (.not.upperBoundary( ivi ))

  ! if (withinBoundaryDens == .true.) -> Density is within table boundaries
  if (withinBoundary( ivi )) then
     if (kPrev > 0 .AND. kPrev < nsteps) then
        k = kPrev
        l = k + 1

        D1 = thisTypeIvar (k)
        D2 = thisTypeIvar (l)

        found = ((D1-xi)*(D2-xi) .LE. 0.0)
     else
        found = .FALSE.
     end if
     if (.NOT.found) then
        do d = 2,nsteps
           if (thisTypeIvar (d) >= xi) then
              k  = d - 1
              l  = d
              D1 = thisTypeIvar (k)
              D2 = thisTypeIvar (l)
#ifndef FLASH_OPENMP
              kPrev = k
#endif              
              exit
           end if
        end do
     end if
  else
      if (lowerBoundary( ivi )) then
          k  = 1
          l  = 1
          D1 = xi
          D2 = thisTypeIvar (1)
      end if

      if (upperBoundary( ivi )) then
          k  = nsteps
          l  = nsteps
          D1 = thisTypeIvar (nsteps)
          D2 = xi
      end if
  end if

  taus( ivi ) =
  tau   = (xi - T1) / (T2 - T1)
  delta = (xi     - D1) / (D2 - D1)


  return
end subroutine eos_gphFindTablePos
