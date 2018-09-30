!!****if* source/physics/Eos/EosMain/Tabulated/eos_gphUpdateOutsideCount
!!
!! NAME
!!
!!  eos_gphUpdateOutsideCount
!!
!! SYNOPSIS
!!
!!  call eos_gphUpdateOutsideCount(integer(in) :: species,
!!                                 logical(in) :: upperboundarytemp,
!!                                 logical(in) :: tempislog,
!!                                 real(in) :: temp,
!!                                 logical(in) :: upperboundarydens,
!!                                 logical(in) :: densislog,
!!                                 real(in) :: dens)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   species : species
!!
!!   upperboundarytemp : check if upper boundary for temperature
!!
!!   tempislog : check if temperature is log
!!
!!   temp : temperature 
!!
!!   upperboundarydens : check if upper boundary for density
!!
!!   densislog : check if density is log
!!
!!   dens : density
!!
!!
!!
!!***

#include "Eos.h"

subroutine eos_gphUpdateOutsideCount(species, &
                                     upperBoundaryTemp, tempIsLog, &
                                     temp, &
                                     upperBoundaryDens, densIsLog, &
                                     dens)

  use Eos_data,                   ONLY : eos_logLevel
  use eos_gphData,                ONLY : EOS_GPH_NALLTAB,           &
                                         eos_gphAllDiag
  implicit none

  integer,intent(in) :: species
  logical,intent(in) :: upperBoundaryTemp, upperBoundaryDens
  logical,intent(in) :: tempIsLog, densIsLog
  real,   intent(in) :: temp, dens

  real :: physTemp, physDens


  if (eos_logLevel .LT. EOS_LOGLEVEL_WARN_ANY) then
     return                     !RETURN IMMEDIATELY
  end if
  if (.NOT. (upperBoundaryTemp .OR. upperBoundaryDens)) then
     return                     !RETURN IMMEDIATELY
  end if
  
  if (tempIsLog) then
     physTemp = 10.0**temp
  else
     physTemp = temp
  end if
  if (densIsLog) then
     physDens = 10.0**dens
  else
     physDens = dens
  end if

  if (upperBoundaryTemp) then
     if (eos_gphAllDiag(species)%highTempCount == 0) then
        eos_gphAllDiag(species)%firstHighTempEvent%temp = physTemp
        eos_gphAllDiag(species)%firstHighTempEvent%dens = physDens
     end if
     eos_gphAllDiag(species)%highTempCount = &
          eos_gphAllDiag(species)%highTempCount + 1
     if (physTemp > eos_gphAllDiag(species)%highestTemp) &
          eos_gphAllDiag(species)%highestTemp = physTemp
  end if

  if (upperBoundaryDens) then
     if (eos_gphAllDiag(species)%highDensCount == 0) then
        eos_gphAllDiag(species)%firstHighDensEvent%temp = physTemp
        eos_gphAllDiag(species)%firstHighDensEvent%dens = physDens
     end if
     eos_gphAllDiag(species)%highDensCount = &
          eos_gphAllDiag(species)%highDensCount + 1
     if (physDens > eos_gphAllDiag(species)%highestDens) &
          eos_gphAllDiag(species)%highestDens = physDens
  end if
  

end subroutine eos_gphUpdateOutsideCount
