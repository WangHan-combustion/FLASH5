!!****ih* source/physics/Eos/localAPI/eos_gphInterface
!!
!! NAME
!!     eos_gphInterface
!!
!! SYNOPSIS
!!     use eos_gphInterface
!!
!! DESCRIPTION
!!
!! This is an interface module for internal use of
!! the GP-Helmholtz Tabulated implementation of the Eos unit.
!!
!!***

module eos_gphInterface
  implicit none
#include "Eos.h"
!#include "Flash.h"

  interface
     subroutine eos_gphBrowseIonmixTables (tableName,                   &
                                       needZFTable,                 &
                                       needENTable,                 &
                                       needHCTable,                 &
                                               nstepsDensityZF,     &
                                               nstepsDensityEN,     &
                                               nstepsDensityHC,     &
                                               nstepsTemperatureZF, &
                                               nstepsTemperatureEN, &
                                               nstepsTemperatureHC  )

       logical,            intent (in)  :: needZFTable
       logical,            intent (in)  :: needENTable
       logical,            intent (in)  :: needHCTable
       integer,            intent (out) :: nstepsDensityZF
       integer,            intent (out) :: nstepsDensityEN
       integer,            intent (out) :: nstepsDensityHC
       integer,            intent (out) :: nstepsTemperatureZF
       integer,            intent (out) :: nstepsTemperatureEN
       integer,            intent (out) :: nstepsTemperatureHC
       character (len=80), intent (in)  :: tableName
     end subroutine eos_gphBrowseIonmixTables
  end interface

  interface
     subroutine eos_gphBrowseTables (tableKind,                   &
                                 tableName,                   &
                                 groupName,                   &
                                 needZFTable,                 &
                                 needENTable,                 &
                                 needHCTable,                 &
                                 needEntrTable,                 &
                                         nstepsDensityZF,     &
                                         nstepsDensityEN,     &
                                         nstepsDensityHC,     &
                                         nstepsDensityEntr,     &
                                         nstepsTemperatureZF, &
                                         nstepsTemperatureEN, &
                                         nstepsTemperatureHC, &
                                         nstepsTemperatureEntr  )

       logical,            intent (in)  :: needZFTable
       logical,            intent (in)  :: needENTable
       logical,            intent (in)  :: needHCTable
       logical,            intent (in)  :: needEntrTable
       integer,            intent (out) :: nstepsDensityZF
       integer,            intent (out) :: nstepsDensityEN
       integer,            intent (out) :: nstepsDensityHC
       integer,            intent (out) :: nstepsDensityEntr
       integer,            intent (out) :: nstepsTemperatureZF
       integer,            intent (out) :: nstepsTemperatureEN
       integer,            intent (out) :: nstepsTemperatureHC
       integer,            intent (out) :: nstepsTemperatureEntr
       character (len=80), intent (in)  :: tableKind
       character (len=80), intent (in)  :: tableName
       character (len=80), intent (in)  :: groupName
     end subroutine eos_gphBrowseTables
  end interface

  interface
subroutine eos_gphFindTablePos (xpos   ,                    &
                                        selTT,              &
                                        ipos, jpos, &
                                        taus, &
                                        iposPrev,          &
                                        lowerBoundary, &
                                        withinBoundary, &
                                        clo, chi)

  implicit none
  real,    intent (in)  :: xpos(EOST_MAX_IVARS)
  integer, intent (in)  :: selTT
  integer,intent(OUT) :: ipos(EOST_MAX_IVARS),jpos(EOST_MAX_IVARS)
  integer,intent(INOUT) :: iposPrev(EOST_MAX_IVARS)
  real,intent(out) :: taus(EOST_MAX_IVARS)
  logical,intent(OUT) :: lowerBoundary(1:EOST_MAX_IVARS)
  logical,intent(OUT) :: withinBoundary(1:EOST_MAX_IVARS)
  real,intent(OUT) :: clo(EOST_MAX_IVARS),chi(EOST_MAX_IVARS)
end subroutine eos_gphFindTablePos
  end interface

  interface
     subroutine eos_gphGetAnyTableData(pos,     &
                                             wanted, &
                                             selTT             , &
                                             derDefs           , &
                                             needDerivs        , &
                                        ipos,         &
                                        tauIn, &
                                        lowerBoundary,&
                                        withinBoundary,&
                                             resultTT, &
                                        cloIn, chiIn)
       use eos_gphData,ONLY: EOS_TAB_NCOMP
       implicit none
       real,    intent (in)  :: pos(EOST_MAX_IVARS)
       logical,            intent (in) :: wanted(EOS_TAB_NCOMP)
       integer, intent (in)  :: selTT
  integer, intent (in)  :: derDefs(:,:)
       integer, intent (in)  :: needDerivs(:,:)
  integer,intent(IN) :: ipos(EOST_MAX_IVARS)
  real,intent(IN) :: tauIn(EOST_MAX_IVARS)
  logical,intent(IN) :: lowerBoundary(1:EOST_MAX_IVARS)
  logical,intent(IN) :: withinBoundary(1:EOST_MAX_IVARS)
       real,    intent (out) :: resultTT(0:,:)
  real,intent(IN),dimension(EOST_MAX_IVARS) :: cloIn, chiIn
     end subroutine eos_gphGetAnyTableData
  end interface

  interface
     subroutine eos_gphGetTabulatedData (xpos,            &
                                             derDefs           , &
                                             wantedTabDeriv    , &
                                                 outData              )
       use eos_gphData,ONLY: EOS_TAB_NCOMP,EOS_GPH_NALLTAB, &
                             EOS_TABVT_ENTR, &
                             EOST_MAX_DERIVS
       implicit none
       real,    intent (in)  :: xpos(EOST_MAX_IVARS)
       integer, intent (in)  :: derDefs(:,0:)
       logical, intent(in) :: wantedTabDeriv(EOST_MAX_DERIVS,EOS_TAB_NCOMP)
!!$       integer, intent (in) :: wantedDerivs(EOS_TAB_NCOMP,EOS_GPH_NALLTAB)
       real,intent(OUT) :: outData(0:EOST_MAX_DERIVS,1:EOS_TABVT_ENTR,EOS_TAB_NCOMP)
     end subroutine eos_gphGetTabulatedData
  end interface

  interface
     subroutine eos_gphIonmix4(mode, vecLen, eosData, vecBegin, vecEnd, eosType, subtype, material, mask)
       implicit none
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       integer,optional,INTENT(in) :: vecBegin,vecEnd
       integer,optional,INTENT(in) :: eosType,subtype
       integer,optional,INTENT(in) :: material
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
     end subroutine eos_gphIonmix4
     subroutine eos_gphGpN(mode, vecLen, eosData, vecBegin, vecEnd, eosType, subtype, mask)
       implicit none
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       integer,optional,INTENT(in) :: vecBegin,vecEnd
       integer,optional,INTENT(in) :: eosType,subtype
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
     end subroutine eos_gphGpN
  end interface

  interface
     subroutine eos_tabOpacplot(mode, vecLen, eosData, vecBegin, vecEnd, eosType, subtype, material, mask)
       implicit none
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       integer,optional,INTENT(in) :: vecBegin,vecEnd
       integer,optional,INTENT(in) :: eosType,subtype
       integer,optional,INTENT(in) :: material
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
     end subroutine Eos_tabOpacplot
  end interface


  interface
     subroutine eos_gphLogOutsideCounts(force)
       implicit none
       logical, INTENT(in) :: force
     end subroutine eos_gphLogOutsideCounts
  end interface

  interface
     subroutine eos_gphReadIonmixTables (tableName,   &
                                     wanted,      &
                                     td,& 
                                     tbZF,tbEN,tbHC)
       use eos_gphData,ONLY: EOS_TAB_NCOMP,EOS_GPH_NALLTAB, &
            eosT_tableIvarsetDescT, &
            eosT_oneVarTablePT

       implicit none

       logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_GPH_NALLTAB)
       type(eosT_tableIvarsetDescT),intent(inout) :: td(:)
       type(eosT_oneVarTablePT),pointer,dimension(:) :: tbZF,tbEN,tbHC

       character (len=80), intent (in) :: tableName
     end subroutine eos_gphReadIonmixTables
  end interface

  interface
     subroutine eos_gphReadIonmix4Tables (tableName,   &
                                     wanted,      &
                                     td,& 
                                     tbZF,tbEN,tbPR,tbHC,tbEntr)
       use eos_gphData,ONLY: EOS_TAB_NCOMP,EOS_GPH_NALLTAB, &
            eosT_tableIvarsetDescT, &
            eosT_oneVarTablePT
       implicit none
       logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_GPH_NALLTAB)
       type(eosT_tableIvarsetDescT),intent(inout) :: td(:)
       type(eosT_oneVarTablePT),pointer,dimension(:) :: tbZF,tbEN,tbPR,tbHC,tbEntr

       character (len=80), intent (in) :: tableName
     end subroutine eos_gphReadIonmix4Tables
     subroutine eos_gphReadGpNTables (tableName,   &
                                     wanted,      &
                                     td,& 
                                     tbEntr)
       use eos_gphData,ONLY: EOS_TAB_NCOMP,EOS_GPH_NALLTAB, &
            eosT_tableIvarsetDescT, &
            eosT_oneVarTablePT
       implicit none
       logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_GPH_NALLTAB)
       type(eosT_tableIvarsetDescT),intent(inout) :: td(:)
       type(eosT_oneVarTablePT),intent(INOUT),dimension(:) :: tbEntr

       character (len=80), intent (in) :: tableName
     end subroutine eos_gphReadGpNTables
  end interface


  interface
     subroutine eos_gphReadTables (tableKind,   &
                               tableName,   &
                               wanted,      &
                               td,&
                                     tbEntr)
       use eos_gphData,ONLY: EOS_TAB_NCOMP,EOS_GPH_NALLTAB, &
            eosT_tableIvarsetDescT, &
            eosT_oneVarTablePT

       implicit none

       logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_GPH_NALLTAB)
       type(eosT_tableIvarsetDescT),intent(inout) :: td(:)
       type(eosT_oneVarTablePT),intent(inout) :: tbEntr

       character (len=80), intent (in) :: tableKind
       character (len=80), intent (in) :: tableName
     end subroutine eos_gphReadTables
  end interface

  interface
     subroutine eos_gphUpdateOutsideCount(species, &
                                     upperBoundaryTemp, tempIsLog, &
                                     temp, &
                                     upperBoundaryDens, densIsLog, &
                                     dens)
       implicit none
       integer,intent(in) :: species
       logical,intent(in) :: upperBoundaryTemp, upperBoundaryDens
       logical,intent(in) :: tempIsLog, densIsLog
       real,   intent(in) :: temp, dens
     end subroutine eos_gphUpdateOutsideCount
  end interface

  interface
     subroutine eos_gphZeroOutsideCounts()
     end subroutine eos_gphZeroOutsideCounts
  end interface

end module eos_gphInterface
