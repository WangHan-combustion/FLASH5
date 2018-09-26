!!****if* source/physics/Eos/EosMain/Tabulated/eos_gphGetSpeciesTabulatedData
!!
!! NAME
!!
!!  eos_gphGetTabulatedData
!!
!! SYNOPSIS
!!
!!  call eos_gphGetTabulatedData (integer (in) :: species,
!!                                        real    (in) :: speciesTemperature,
!!                                        real    (in) :: speciesDensity,
!!                                        integer (in) :: maxComp,
!!                                        integer (in) :: needZFTable,
!!                                        integer (in) :: needENDerivs,
!!                                        integer (in) :: needHCDerivs)
!!
!! DESCRIPTION
!!
!!  Computes absorption, emission and transport opacities from tabulated
!!  values for a particular (temperature, density, energy group, species)
!!  quadruple.
!!
!! ARGUMENTS
!!
!!   species            : The species index
!!   speciesTemperature : The species temperature
!!   speciesDensity     : The species density
!!   needZFDerivs        : Knob to activate data extraction from the average ionization tables.
!!                         Number of highest derivative needed, may be 0, negative for none.
!!   needENDerivs        : Knob to activate data extraction from the internal energy tables
!!                         Number of highest derivative needed, may be 0, negative for none.
!!   needHCDerivs        : Knob to activate data extraction from the specific heat tables
!!                         Number of highest derivative needed, may be 0, negative for none.
!!
!! NOTES
!!
!!  See definitions of EOS_TABINT_DERIV_* in eos_gphData and EOS_TAB_NDERIVS for
!!  the numbering of derivatives.
!!***
subroutine eos_gphGetTabulatedData (species,            &
                                            speciesTemperature, &
                                            speciesDensity,     &
                                            wantedDerivs, &
                                            outData              )

  use eos_gphInterface, ONLY: eos_gphGetAnyTableData, eos_gphFindTablePos, &
                              eos_gphUpdateOutsideCount
  use eos_gphData, ONLY: eos_gphIonizationKind,     &
                          eos_gphIntEnergyKind,       &
                          eos_gphHeatCpKind,      &
                          EOS_TAB_NCOMP,EOS_TAB_NALLTAB, &
                          EOS_TAB_NDERIVS, &
                          EOS_TAB_FOR_ION, EOS_TAB_FOR_ELE, EOS_TAB_FOR_MAT, &
                          EOS_TABVT_ZF, EOS_TABVT_EN, EOS_TABVT_PR, EOS_TABVT_HC, &
                          EOS_TABVT_ENTR, &
                          EOS_TABINT_DERIV_0, EOS_TAB_NDERIVS,      &
                          eos_allTab
  implicit none
  
#include "Eos.h"
  
  integer, intent (in) :: wantedDerivs(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
  real,    intent (in) :: speciesTemperature
  real,    intent (in) :: speciesDensity
  real,intent(OUT) :: outData(0:EOS_TAB_NDERIVS,1:EOS_TABVT_ENTR,EOS_TAB_NCOMP)

  real,    intent (in)  :: xpos(EOST_MAX_IVARS)

  logical :: lowerBoundary(1:EOST_MAX_IVARS)
  logical :: withinBoundary(1:EOST_MAX_IVARS)
  logical :: isLog
  integer :: i,j,k,l
  integer :: varType
  integer,save :: iSave=1,kSave=1
  real :: D1,D2,T1,T2
  real :: taus(EOST_MAX_IVARS)


  iSave = -1; kSave = -1
  call eos_gphFindTablePos (xpos,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        EOS_TABULAR_S,       &
                                        i,j,k,l,&
                                        taus, &
                                        iSave, kSave, &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        T1,T2,D1,D2)

  varType = EOS_TABVT_ZF
!!$  thisTypeTable => eos_allTab(species)%tg(varType)
  isLog = eos_allTab(species)%tg(varType)%td%isLog
  call eos_gphUpdateOutsideCount( &
                                 .NOT.(lowerBoundaryTemp .OR. withinBoundaryTemp), &
                                 isLog, T2, &
                                 .NOT.(lowerBoundary .OR. withinBoundaryDens), &
                                 isLog, D2)

!
!
!   ...Extract only the necessary opacities from the tables.
!
!
  if (ANY(wantedDerivs(:,EOS_TABVT_ZF) .GE. 0)) then

      call eos_gphGetAnyTableData (xpos,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        (wantedDerivs .GE. 0),&
                                        EOS_TABULAR_Z,      &
                                        maxval(wantedDerivs(:,EOS_TABVT_ZF)),       &
                                        i,j,k,l,&
                                        tau, delta, &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        outData(:,EOS_TABVT_ZF,:), &
                                        T1,T2,D1,D2)
  end if

  if (ANY(wantedDerivs(:,EOS_TABVT_EN) .GE. 0)) then

      call eos_gphGetAnyTableData (xpos,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        (wantedDerivs .GE. 0),&
                                        EOS_TABULAR_E,      &
                                        maxval(wantedDerivs(:,EOS_TABVT_EN)),       &
                                        i,j,k,l,&
                                        tau, delta, &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        outData(:,EOS_TABVT_EN,:), &
                                        T1,T2,D1,D2)
  end if

  if (ANY(wantedDerivs(:,EOS_TABVT_HC) .GE. 0)) then

      call eos_gphGetAnyTableData (xpos,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        (wantedDerivs .GE. 0),&
                                        EOS_TABULAR_C,      &
                                        maxval(wantedDerivs(:,EOS_TABVT_HC)),       &
                                        i,j,k,l,&
                                        tau, delta, &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        outData(:,EOS_TABULAR_C,:), &
                                        T1,T2,D1,D2)
  end if

  if (ANY(wantedDerivs(:,EOS_TABVT_PR) .GE. 0)) then

      call eos_gphGetAnyTableData (xpos,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        (wantedDerivs .GE. 0),&
                                        EOS_TABULAR_P,      &
                                        maxval(wantedDerivs(:,EOS_TABVT_PR)),       &
                                        i,j,k,l,&
                                        tau, delta, &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        outData(:,EOS_TABVT_PR,:), &
                                        T1,T2,D1,D2)
  end if

  if (ANY(wantedDerivs(:,EOS_TABVT_ENTR) .GE. 0)) then

      call eos_gphGetAnyTableData (xpos,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        (wantedDerivs .GE. 0),&
                                        EOS_TABULAR_S,      &
                                        maxval(wantedDerivs(:,EOS_TABVT_ENTR)),       &
                                        i,j,k,l,&
                                        tau, delta, &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        outData(:,EOS_TABVT_ENTR,:), &
                                        T1,T2,D1,D2)
  end if

!
!
!   ...Ready! 
!
!
  return
end subroutine eos_gphGetTabulatedData
