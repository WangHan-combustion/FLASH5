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
subroutine eos_gphGetTabulatedData (xpos,   &
                                    derDefs,        &
                                    wantedTabDeriv, &
                                            outData              )

  use eos_gphInterface, ONLY: eos_gphGetAnyTableData, eos_gphFindTablePos, &
                              eos_gphUpdateOutsideCount
  use eos_gphData, ONLY: eos_gphIonizationKind,     &
                          eos_gphIntEnergyKind,       &
                          eos_gphHeatCpKind,      &
                          EOS_TAB_NCOMP,EOS_GPH_NALLTAB, &
                          EOST_MAX_DERIVS, &
                          EOS_TAB_FOR_MAT, &
                          EOS_TABVT_ENTR, &
                          ENER_IVAR, DENS_IVAR, &
                          EOS_TABINT_DERIV_0, EOST_MAX_DERIVS,      &
                                eosT_varTableGroupT,      &
                          eos_gphTheTable
  use eos_gphData, ONLY: EOS_GPHDERIV_D , &
                         EOS_GPHDERIV_E , &
                         EOS_GPHDERIV_E2
  implicit none
  
#include "Eos.h"
  
!!$  integer, intent(in) :: derDefs(EOST_MAX_IVARS,EOST_MAX_DERIVS)
  integer, intent(in) :: derDefs(:,0:)
  logical, intent(in) :: wantedTabDeriv(EOST_MAX_DERIVS,EOS_TAB_NCOMP)
  integer,allocatable:: wantedDerivs(:,:)
  real,intent(OUT) :: outData(0:EOST_MAX_DERIVS,1:EOS_TABVT_ENTR,EOS_TAB_NCOMP)

  real,intent (in)  :: xpos(EOST_MAX_IVARS)

  type(eosT_varTableGroupT),pointer :: thisTypeTable

  logical :: lowerBoundary(1:EOST_MAX_IVARS)
  logical :: withinBoundary(1:EOST_MAX_IVARS)
  logical :: isLog1, isLog2
  logical :: wantedComp(EOS_TAB_NCOMP)
  integer :: ipos(EOST_MAX_IVARS),jpos(EOST_MAX_IVARS)
  integer :: iposPrev(EOST_MAX_IVARS)
  integer :: varType
  integer :: numDerivs
  integer :: i
  integer,save :: iSave(EOST_MAX_IVARS)=1

  real :: taus(EOST_MAX_IVARS)
  real :: clo(EOST_MAX_IVARS),chi(EOST_MAX_IVARS)

!!$  print*,'eos_gphGetTabulatedData: derDefs=',derDefs
  iSave = -1
  call eos_gphFindTablePos (xpos,            &
                                        EOS_TABULAR_S,       &
                                        ipos, jpos,&
                                        taus, &
                                        iSave, &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        clo, chi)

  varType = EOS_TABVT_ENTR

  print*,'taus',taus
  print*,'clo',clo,', chi',chi

  thisTypeTable => eos_gphTheTable%tg(varType)
  isLog1 = thisTypeTable%td%c(ENER_IVAR)%isLog
  isLog2 = thisTypeTable%td%c(DENS_IVAR)%isLog
  call eos_gphUpdateOutsideCount(EOS_TABVT_ENTR, &
                                 .NOT.(lowerBoundary(ENER_IVAR) .OR. withinBoundary(ENER_IVAR)), &
                                 isLog1, chi(ENER_IVAR), &
                                 .NOT.(lowerBoundary(DENS_IVAR) .OR. withinBoundary(DENS_IVAR)), &
                                 isLog2, chi(DENS_IVAR))

  wantedComp = .FALSE.
  wantedComp(EOS_TAB_FOR_MAT) = .TRUE.

  numDerivs = 0
  if (wantedTabDeriv(EOS_GPHDERIV_E  ,EOS_TAB_FOR_MAT)) numDerivs = 1
  if (wantedTabDeriv(EOS_GPHDERIV_D  ,EOS_TAB_FOR_MAT)) numDerivs = numDerivs + 1
  if (wantedTabDeriv(EOS_GPHDERIV_E2 ,EOS_TAB_FOR_MAT)) numDerivs = numDerivs + 1

  allocate(wantedDerivs(EOST_MAX_IVARS,numDerivs))
  wantedDerivs(:,:) = 0

  i = 0
  if (wantedTabDeriv(EOS_GPHDERIV_E  ,EOS_TAB_FOR_MAT)) then
     i = 1
     wantedDerivs(:,i) = derDefs(:,EOS_GPHDERIV_E)
  end if
  if (wantedTabDeriv(EOS_GPHDERIV_D  ,EOS_TAB_FOR_MAT)) then
     i = i + 1
     wantedDerivs(:,i) = derDefs(:,EOS_GPHDERIV_D)
  end if
  if (wantedTabDeriv(EOS_GPHDERIV_E2 ,EOS_TAB_FOR_MAT)) then
     i = i + 1
     wantedDerivs(:,i) = derDefs(:,EOS_GPHDERIV_E2)
  end if

!
!   ...Extract only the necessary variables from the tables.
!
!
#if 0
  if (ANY(wantedDerivs(:,EOS_TABVT_ZF) .GE. 0)) then

      call eos_gphGetAnyTableData (xpos,            &
                                        (wantedDerivs .GE. 0),&
                                        EOS_TABULAR_Z,      &
                                        maxval(wantedDerivs(:,EOS_TABVT_ZF)),       &
                                        ipos, &
                                        taus, &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        outData(:,EOS_TABVT_ZF,:), &
                                        clo, chi)
  end if

  if (ANY(wantedDerivs(:,EOS_TABVT_EN) .GE. 0)) then

      call eos_gphGetAnyTableData (xpos,            &
                                        (wantedDerivs .GE. 0),&
                                        EOS_TABULAR_E,      &
                                        maxval(wantedDerivs(:,EOS_TABVT_EN)),       &
                                        ipos, &
                                        taus, &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        outData(:,EOS_TABVT_EN,:), &
                                        clo, chi)
  end if

  if (ANY(wantedDerivs(:,EOS_TABVT_HC) .GE. 0)) then

      call eos_gphGetAnyTableData (xpos,            &
                                        (wantedDerivs .GE. 0),&
                                        EOS_TABULAR_C,      &
                                        maxval(wantedDerivs(:,EOS_TABVT_HC)),       &
                                        ipos, &
                                        taus, &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        outData(:,EOS_TABULAR_C,:), &
                                        clo, chi)
  end if

  if (ANY(wantedDerivs(:,EOS_TABVT_PR) .GE. 0)) then

      call eos_gphGetAnyTableData (xpos,            &
                                        (wantedDerivs .GE. 0),&
                                        EOS_TABULAR_P,      &
                                        maxval(wantedDerivs(:,EOS_TABVT_PR)),       &
                                        ipos, &
                                        taus, &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        outData(:,EOS_TABVT_PR,:), &
                                        clo, chi)
  end if
#endif

  if (ANY(wantedDerivs(:,EOS_TABVT_ENTR) .GE. 0)) then


     do i=lbound(wantedDerivs,2),ubound(wantedDerivs,2)
        outData(i,1:EOS_TABVT_ENTR,EOS_TAB_NCOMP) = 100 + i
     end do
     do i=0,EOST_MAX_DERIVS
        outData(i,1:EOS_TABVT_ENTR,EOS_TAB_NCOMP) = 200 + i
     end do
     print*,'gphGeTabulatedData: sought xpos is',xpos
      call eos_gphGetAnyTableData (xpos,            &
                                        wantedComp,&
                                        EOS_TABULAR_S,      &
                                        derDefs,            &
                                        wantedDerivs,       &
                                        ipos, &
                                        taus, &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        outData(:,EOS_TABVT_ENTR,:), &
                                        clo, chi)
  end if

!
!
!   ...Ready! 
!
!
  deallocate(wantedDerivs)
  return
end subroutine eos_gphGetTabulatedData
