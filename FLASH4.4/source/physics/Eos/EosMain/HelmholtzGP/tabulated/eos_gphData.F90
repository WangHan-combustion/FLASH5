!!****if* source/physics/Eos/EosMain/Tabulated/eos_gphData
!!
!! NAME
!!
!!  eos_gphData
!!
!! SYNOPSIS
!!  use eos_gphData
!!
!! DESCRIPTION
!!
!!  Defines and stores the local data for Unit Opacity.
!!  
!!***

#include "Eos.h"
#include "Eos_components.h"

module eos_gphData
  
  implicit none

#if (EOSCOMP_NUM_COMPONENTS==1)
!!notdef  integer, parameter :: EOS_TAB_FOR_ION = 1
!!notdef  integer, parameter :: EOS_TAB_FOR_ELE = 1
  integer, parameter :: EOS_TAB_FOR_MAT = 1
#else
!!notdef  integer, parameter :: EOS_TAB_FOR_ION = 1
!!notdef  integer, parameter :: EOS_TAB_FOR_ELE = 2
!!notdef  integer, parameter :: EOS_TAB_FOR_MAT = 3
#endif

  ! maximum in this code version is 4 input variables (canonical: spec. energy, density, Ye, A or similar)
#define EOST_MAX_IVARS 4
  integer, parameter :: EOS_TAB_NCOMP = EOS_TAB_FOR_MAT

  ! variable types - families of tabulated variables
!!$  integer, parameter :: EOS_TABVT_ZF = 1 !Z_free
  integer, parameter :: EOS_TABVT_EN = 2 !(internal) energies
  integer, parameter :: EOS_TABVT_PR = 3 !pressures
  integer, parameter :: EOS_TABVT_ENTR = 1 !entropies
  integer, parameter :: EOS_TABVT_HC = 5 !heat capacities
  integer, parameter :: EOS_TABVT_OP = 6 !opacities
  integer, parameter :: EOS_GPH_NALLTAB = 1 ! number of output vars (aka saved state functions, aka ovars) supported

  integer, parameter :: EOS_TABINT_DERIV_0 = 0 ! leave as 0
  integer, parameter :: EOS_TABINT_DERIV_DT = 1 ! leave consecutive
  integer, parameter :: EOS_TABINT_DERIV_DD = 2 ! leave consecutive
  integer, parameter :: EOS_TAB_NDERIVS = 2


  logical, save :: eos_useLogTables

  integer, save :: op_maxNstepsDensityZF
  integer, save :: op_maxNstepsDensityEN
  integer, save :: op_maxNstepsDensityHC
  integer, save :: op_maxNstepsTemperatureZF
  integer, save :: op_maxNstepsTemperatureEN
  integer, save :: op_maxNstepsTemperatureHC
  integer, save :: op_maxTablesZF
  integer, save :: op_maxTablesEN
  integer, save :: op_maxTablesHC
  integer, save :: eos_gphTotalNumSpecies


  real, parameter :: zero =  0.0
  real, parameter :: one  =  1.0
  real, parameter :: ten  = 10.0

  character (len=80), allocatable, save ::  eos_gphleKind (:)
  character (len=80), allocatable, save ::  eos_gphleName (:)
  character (len=80), allocatable, save ::  eos_groupName (:)

  integer, allocatable, save :: eos_gphIonizationKind    (:,:)
  integer, allocatable, save :: eos_gphIntEnergyKind     (:,:)
  integer, allocatable, save :: eos_gphHeatCpKind        (:,:)


  type eosT_oneVarTablePT
     type(eosT_varTableGroupT),pointer :: pg ! group to which this table belongs
     real, pointer                      :: table(:,:) ! the data, for one cell
!!$     logical                            :: isLogData
!!$     character(len=80)                  :: fromFile !for debugging
!!$     integer                            :: tableNo  !for debugging
!!$     integer                            :: lineNo   !for debugging
!!$     integer                            :: derivedFrom1 !for debugging ?
!!$     integer                            :: derivedFrom2 !for debugging ?
!!$     integer                            :: varType ! whether (1) Z, (2) eint, (3) pres, (4) hc, ..
!!$     integer                            :: component ! whether (1) EOS_TAB_FOR_ION,
!!$                                                     ! (2) EOS_TAB_FOR_ELE, (3) EOS_TAB_FOR_MAT, ?..
!     integer                            :: iSave     ! cached location of previous lookup
!     integer                            :: kSave     ! cached location of previous lookup
  end type eosT_oneVarTablePT

  integer,dimension(EOST_MAX_IVARS) :: tdims
  tdims(:) = 1
  ASSERT(  TheGphTable % tg(EOS_TABVT_ENTR) % td % N  == tableDim  )
  do i=1,tableDim
     tdims(i) = TheGphTable % tg(EOS_TABVT_ENTR) % td % c(i) % nIval
  end do
  allocate( d (0:ncorners-1, 0:nderivs, tdims(1), tdims(2), tdims(2), tdims(4)) )
  TheGphTable % tg(EOS_TABVT_ENTR) % table => d
  
  
  type eosT_oneVarTablePT
     type(eosT_varTableGroupT),pointer :: pg ! group to which this table belongs
     real, pointer                      :: table(:,:,:,:,:,:) ! the data, (2+EOST_MAX_IVARS)-dimensional
!!$     logical                            :: isLogData
!!$     character(len=80)                  :: fromFile !for debugging
!!$     integer                            :: tableNo  !for debugging
!!$     integer                            :: lineNo   !for debugging
!!$     integer                            :: derivedFrom1 !for debugging ?
!!$     integer                            :: derivedFrom2 !for debugging ?
!!$     integer                            :: varType ! whether (1) Z, (2) eint, (3) pres, (4) hc, ..
!!$     integer                            :: component ! whether (1) EOS_TAB_FOR_ION,
!!$                                                     ! (2) EOS_TAB_FOR_ELE, (3) EOS_TAB_FOR_MAT, ?..
!     integer                            :: iSave     ! cached location of previous lookup
!     integer                            :: kSave     ! cached location of previous lookup
  end type eosT_oneVarTablePT


  type eosT_tableIvarDescT
     integer                           :: nIval  = 0
     real, pointer                     :: val(:) => null() ! to be allocated to length nIval+1
     logical                           :: isLog = .FALSE. !whether *temperatures* and *densities* are stored as logarithms.
  end type eosT_tableIvarDescT

  type eosT_tableIvarsetDescT
     integer                           :: n
     type(eosT_tableIvarDescT)         :: c(EOST_MAX_IVARS) ! up to EOST_MAX_IVARS table "coordinates"
!!$     logical                           :: isLog = .FALSE. !whether any *temperatures* and *densities* etc. are stored as logarithms?
  end type eosT_tableIvarsetDescT

  ! Now a type that can contain several (up-to-)EOST_MAX_IVARS-dimensional data tables, each of them
  ! of the same size an row and columns ranges.
  ! We shall use one object of this type for each species for each of (Z, eint, pres, hc)
  ! where necessary.
  type eosT_varTableGroupT
     type(eosT_oneVarTablePT)           :: table(EOS_TAB_FOR_MAT:EOS_TAB_FOR_MAT) ! tables for one (or, after adjusting code, several) components
     type(eosT_tableIvarsetDescT)        :: td
!     integer                           :: ntemp
!     integer                           :: ndens
!     integer                           :: nmg ! length of vector in multigroup variables
!     real, pointer                     :: Temperatures(:)
!     real, pointer                     :: Densities(:)
  end type eosT_varTableGroupT

  ! Now a type that contains all the tables (or pointers to them) that
  ! pertain to one output var (aka saved state function, aka ovar).
  type eosT_varAllTablesT
     type(eosT_varTableGroupT)  :: tg(1:EOS_GPH_NALLTAB)
  end type eosT_varAllTablesT

  type(eosT_varAllTablesT), allocatable,target,save :: eos_gphTheTable

  type eosT_diagEvent
     real :: temp
     real :: dens
  end type eosT_diagEvent

  type eosT_diagPT
     integer :: highTempCount
     integer :: highDensCount
     real :: highestTemp
     real :: highestDens
     type(eosT_diagEvent) :: firstHighTempEvent
     type(eosT_diagEvent) :: firstHighDensEvent
     logical :: highTempVarsLookedUp(1:EOS_GPH_NALLTAB)
     logical :: highDensVarsLookedUp(1:EOS_GPH_NALLTAB)
  end type eosT_diagPT


  type(eosT_diagPT), allocatable,target,save :: eos_gphAllDiag(:)



end module eos_gphData

