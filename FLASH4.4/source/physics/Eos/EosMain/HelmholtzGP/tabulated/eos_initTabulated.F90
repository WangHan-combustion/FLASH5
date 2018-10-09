!!****if* source/physics/Eos/EosMain/Tabulated/eos_initTabulated
!!
!! NAME
!!
!!  eos_initTabulated
!!
!! SYNOPSIS
!!
!!  call eos_initTabulated()
!!                      
!!                    
!!
!! DESCRIPTION
!!
!!  Initializes data for Tabulated implementations of the Eos unit.
!!
!! ARGUMENTS
!!
!!  none  
!!  
!!  
!!
!!***

#include "Eos.h"
#include "Flash.h"
#include "constants.h"

subroutine eos_initTabulated()

  use eos_gphData,                ONLY : EOS_TAB_FOR_ION => EOS_TAB_FOR_MAT,&
                                         EOS_TAB_FOR_ELE => EOS_TAB_FOR_MAT,&
                                         EOS_TAB_FOR_MAT,           &
                                         EOS_TAB_NALLTAB => EOS_GPH_NALLTAB,           &
                                         EOS_GPH_NALLTAB,           &
                                         EOS_TABVT_ZF =>EOS_TABVT_ENTR  ,           &
                                         EOS_TABVT_EN   ,           &
                                         EOS_TABVT_HC   ,           &
                                         EOS_TABVT_PR   ,           &
                                         EOS_TABVT_ENTR ,           &
                                          op_maxNstepsDensityZF,     &
                                          op_maxNstepsDensityEN,     &
                                          op_maxNstepsDensityHC,     &
                                          op_maxNstepsTemperatureZF, &
                                          op_maxNstepsTemperatureEN, &
                                          op_maxNstepsTemperatureHC, &
                                          op_maxTablesZF,            &
                                          op_maxTablesEN,            &
                                          op_maxTablesHC,            &
                                          EOS_TAB_NCOMP=> EOS_TAB_FOR_MAT,          &
                                          eos_useLogTables,           &
                                          eos_gphKind,              &
                                          eos_gphFileName,              &
                                          eos_gphIonizationKind,         &
                                          eos_gphIntEnergyKind,           &
                                          eos_gphHeatCpKind,          &
                                          TheGphTable => eos_gphTheTable, &
                                          eos_gphAllDiag
  use eos_gphData,                 ONLY : eos_maxFactorUp, eos_maxFactorDown


  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Driver_interface,            ONLY : Driver_abortFlash
!!$  use eos_gphInterface,            ONLY : eos_tabBrowseTables, eos_gphReadTables, &
!!$                                          eos_tabWriteTables
  use eos_gphInterface,            ONLY : eos_gphReadTables, eos_gphReadGpNTables
  
  use Eos_data,                    ONLY : eos_type
  use Eos_data,                    ONLY : eos_globalMe, eos_meshMe, eos_entrEleScaleChoice

  implicit none

  character*80 :: dummyLine
  character*80 :: printoutHeader

  logical :: needZFTable
  logical :: needENTable
  logical :: needHCTable
  logical :: needPRTable
  logical :: needEntrTable
  logical :: needTable
  logical :: isGpNLike
  logical :: isGpN, isIonmix6
  logical :: isOpacplot
  logical :: wanted(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)

  integer :: i, j, ivi

  integer :: fileUnit
  integer :: nstepsDensityZF
  integer :: nstepsDensityEN
  integer :: nstepsDensityHC
  integer :: nstepsDensityEntr
  integer :: nstepsTemperatureZF
  integer :: nstepsTemperatureEN
  integer :: nstepsTemperatureHC
  integer :: nstepsTemperatureEntr
  integer :: ut_getFreeFileUnit
  integer :: tabno


  real,pointer :: temperatures(:)
  real,pointer :: densities(:)


  integer,dimension(EOST_MAX_IVARS) :: tdims

#if 0
  tdims(:) = 1
  ASSERT(  TheGphTable % tg(EOS_TABVT_ENTR) % td % N  == tableDim  )
  ASSERT(  TheGphTable % tg(EOS_TABVT_ENTR) % td % nderivs  == nderivs  )
  TheGphTable % tg(EOS_TABVT_ENTR) % td % ncorners = 2**tableDim
  do i=1,tableDim
     tdims(i) = TheGphTable % tg(EOS_TABVT_ENTR) % td % c(i) % nIval
  end do
  allocate( d (0:ncorners-1, 0:nderivs, tdims(1), tdims(2), tdims(2), tdims(4)) )
  TheGphTable % tg(EOS_TABVT_ENTR) % table => d
#endif
!
!
!    ...Set internal parameters.
!
!
  ! DONT HAVE ANY ?
!
!
!    ...Allocate those arrays that depend on the # of species and # of energy groups alone.
!
!
  allocate (eos_gphKind)
  allocate (eos_gphFileName)
  allocate (eos_gphIonizationKind)
  allocate (eos_gphIntEnergyKind)
  allocate (eos_gphHeatCpKind)

!
!
!    ...Get the external data.
!
!
!!$  call RuntimeParameters_get ("eos_useLogTables",   eos_useLogTables )

  call RuntimeParameters_get('eos_maxFactorUp', eos_maxFactorUp)
  call RuntimeParameters_get('eos_maxFactorDown', eos_maxFactorDown)

  call RuntimeParameters_get("eos_entrEleScaleChoice", eos_entrEleScaleChoice)
  call RuntimeParameters_get("eos_gphFileName", eos_gphFileName)


  ! noLOOP # 0 - Allocation & zeroing of pointers

  allocate(theGphTable)
     do i = 1,EOS_TAB_NALLTAB
        theGphTable%tg(i)%td%N        = 0 ! to change, seen below
        theGphTable%tg(i)%td%ncorners = 0 ! to change, seen below
        theGphTable%tg(i)%td%nderivs  = 0 ! to change, seen below
        !nullify pointers for coordinate arrays
        do ivi = 1,EOST_MAX_IVARS
           theGphTable%tg(i)%td%c(ivi)%nIval = 0 ! to change seen below
           nullify(theGphTable%tg(i)%td%c(ivi)%val)
           theGphTable%tg(i)%td%c(ivi)%isLog = .FALSE.
        end do
        !nullify pointers for tables of all types
        do j = LBOUND(theGphTable%tg(i)%table,1), &
             UBOUND(theGphTable%tg(i)%table,1)
           nullify(theGphTable % tg(1) % table(j)%table)
           nullify(theGphTable % tg(1) % table(j)%ells)
           nullify(theGphTable % tg(1) % table(j)%derDefs)
        end do
!!$        if (i==EOS_TABVT_ENTR) then
!!$           theGphTable%tg(i)%td%isLog = eos_useLogTables
!!$        else
!!$           theGphTable%tg(i)%td%isLog = .FALSE.
!!$        end if
     end do
!
!
!    ...Initialize all data to zero.
!
!
  op_maxTablesZF            = 0
  op_maxTablesEN            = 0
  op_maxTablesHC            = 0
  op_maxNstepsDensityZF     = 0
  op_maxNstepsDensityEN     = 0
  op_maxNstepsDensityHC     = 0
  op_maxNstepsTemperatureZF = 0
  op_maxNstepsTemperatureEN = 0
  op_maxNstepsTemperatureHC = 0

  eos_gphIonizationKind         = 0
  eos_gphIntEnergyKind           = 0
  eos_gphHeatCpKind          = 0

  wanted = .FALSE.
!
!
!    ...Assume that Multispecies_init as initialized the relevant fields in
!    the multispecies database.
!    !!DEV: Some of that initialization should be moved to the Eos unit -
!    filenames, at least. - KW
!
!

  eos_type = EOS_TAB

!
!
!    ...Read in the needed info from the external file and set the handle
!       arrays.
!
!
  ! noLOOP # 1 - Get info from the Multispecies database that
  ! used to be read from the EOS_sources file.

     eos_gphKind = "gphN"
     
     isOpacplot  = .FALSE.
     isGpN   = .FALSE.
     isIonmix6   = .FALSE.

     eos_gphIonizationKind = EOS_APPROX_KIN
     eos_gphIntEnergyKind = EOS_APPROX_KIN
     eos_gphHeatCpKind = EOS_APPROX_KIN

  

!
!
  ! noLOOP # 2 - Browse & Allocate
  !
!    ...Determine first the overall maximal dimensions needed for allocating
!       the tables.
!
!    ...Allocate the tables and associated data for interpolation.
!

     isOpacplot  = (eos_gphKind=='OPACPLOT')
     isGpN   = (eos_gphKind=='gphN')
     isIonmix6   = (eos_gphKind=='IONMIX6')
     isGpNLike = (isGpN .OR. isIonmix6)

     needZFTable  =     ((eos_gphIonizationKind == EOS_TABULAR_Z) &
                    .or. (eos_gphIntEnergyKind == EOS_TABULAR_Z) &
                    .or. (eos_gphHeatCpKind == EOS_TABULAR_Z) )

     needENTable  =     ((eos_gphIonizationKind == EOS_TABULAR_E) &
                    .or. (eos_gphIntEnergyKind == EOS_TABULAR_E) &
                    .or. (eos_gphHeatCpKind == EOS_TABULAR_E))

     needHCTable  =     ((eos_gphIonizationKind == EOS_TABULAR_C) &
                    .or. (eos_gphIntEnergyKind == EOS_TABULAR_C) &
                    .or. (eos_gphHeatCpKind == EOS_TABULAR_C))!!!   .AND. .FALSE. !!!!!!

     needPRTable  =     ((eos_gphIonizationKind == EOS_TABULAR_P) &
                    .or. (eos_gphIntEnergyKind == EOS_TABULAR_P) &
                    .or. (eos_gphHeatCpKind == EOS_TABULAR_P))

     needEntrTable = isIonmix6

     needTable    =       needZFTable &
                    .or.  needENTable &
                    .or.  needHCTable &
                    .or.  needPRTable &
                    .or.  needEntrTable

     if (needZFTable) op_maxTablesZF = op_maxTablesZF + 1
     if (needENTable) op_maxTablesEN = op_maxTablesEN + 1
     if (needHCTable) op_maxTablesHC = op_maxTablesHC + 1
     if (needPRTable) op_maxTablesHC = op_maxTablesHC + 1

#if 0
     if (needTable) then

         call eos_tabBrowseTables (eos_gphKind,      &
                               eos_gphFileName,      &
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
                                       nstepsTemperatureEntr)

         op_maxNstepsDensityZF     = max (nstepsDensityZF     , op_maxNstepsDensityZF    )
         op_maxNstepsDensityEN     = max (nstepsDensityEN     , op_maxNstepsDensityEN    )
         op_maxNstepsDensityHC     = max (nstepsDensityHC     , op_maxNstepsDensityHC    )
         op_maxNstepsTemperatureZF = max (nstepsTemperatureZF , op_maxNstepsTemperatureZF)
         op_maxNstepsTemperatureEN = max (nstepsTemperatureEN , op_maxNstepsTemperatureEN)
         op_maxNstepsTemperatureHC = max (nstepsTemperatureHC , op_maxNstepsTemperatureHC)

         if (needZFTable) then
            theGphTable%tg(EOS_TABVT_ZF)%td%ntemp = nstepsTemperatureZF
            theGphTable%tg(EOS_TABVT_ZF)%td%ndens = nstepsDensityZF
            nullify(theGphTable%tg(EOS_TABVT_ZF)%td%Temperatures)
            nullify(theGphTable%tg(EOS_TABVT_ZF)%td%Densities)
!!$            allocate(theGphTable%tg(EOS_TABVT_ZF)%td%Temperatures(nstepsTemperatureZF))
!!$            allocate(theGphTable%tg(EOS_TABVT_ZF)%td%Densities(nstepsDensityZF))
            allocate(theGphTable%tg(EOS_TABVT_ZF)%table(EOS_TAB_FOR_ELE:EOS_TAB_FOR_ELE))
            theGphTable%tg(EOS_TABVT_ZF)%numTables = 1
            do i = LBOUND(theGphTable%tg(EOS_TABVT_ZF)%table,1), &
                 UBOUND(theGphTable%tg(EOS_TABVT_ZF)%table,1)
               nullify(theGphTable%tg(EOS_TABVT_ZF)%table(i)%table)
               theGphTable%tg(EOS_TABVT_ZF)%table(i)%isLogData = .FALSE.
            end do
         End if
         if (needENTable) then
            theGphTable%tg(EOS_TABVT_EN)%td%ntemp = nstepsTemperatureEN
            theGphTable%tg(EOS_TABVT_EN)%td%ndens = nstepsDensityEN
            nullify(theGphTable%tg(EOS_TABVT_EN)%td%Temperatures)
            nullify(theGphTable%tg(EOS_TABVT_EN)%td%Densities)
!!$            allocate(theGphTable%tg(EOS_TABVT_EN)%td%Temperatures(nstepsTemperatureEN))
!!$            allocate(theGphTable%tg(EOS_TABVT_EN)%td%Densities(nstepsDensityEN))
            if (isGpNLike) then
               allocate(theGphTable%tg(EOS_TABVT_EN)%table(EOS_TAB_NCOMP))
               theGphTable%tg(EOS_TABVT_EN)%numTables = EOS_TAB_NCOMP
            else
               allocate(theGphTable%tg(EOS_TABVT_EN)%table(EOS_TAB_FOR_MAT:EOS_TAB_FOR_MAT))
               theGphTable%tg(EOS_TABVT_EN)%numTables = 1 
            end if
            do i = LBOUND(theGphTable%tg(EOS_TABVT_EN)%table,1), &
                 UBOUND(theGphTable%tg(EOS_TABVT_EN)%table,1)
               nullify(theGphTable%tg(EOS_TABVT_EN)%table(i)%table)
               theGphTable%tg(EOS_TABVT_EN)%table(i)%isLogData = .FALSE.
            end do
         end if
         if (needHCTable) then
            theGphTable%tg(EOS_TABVT_HC)%td%ntemp = nstepsTemperatureHC
            theGphTable%tg(EOS_TABVT_HC)%td%ndens = nstepsDensityHC
            nullify(theGphTable%tg(EOS_TABVT_HC)%td%Temperatures)
            nullify(theGphTable%tg(EOS_TABVT_HC)%td%Densities)
!!$            allocate(theGphTable%tg(EOS_TABVT_HC)%td%Temperatures(nstepsTemperatureHC))
!!$            allocate(theGphTable%tg(EOS_TABVT_HC)%td%Densities(nstepsDensityHC))
            allocate(theGphTable%tg(EOS_TABVT_HC)%table(EOS_TAB_NCOMP))
            theGphTable%tg(EOS_TABVT_HC)%numTables = EOS_TAB_NCOMP
            do i = LBOUND(theGphTable%tg(EOS_TABVT_HC)%table,1), &
                 UBOUND(theGphTable%tg(EOS_TABVT_HC)%table,1)
               nullify(theGphTable%tg(EOS_TABVT_HC)%table(i)%table)
               theGphTable%tg(EOS_TABVT_HC)%table(i)%isLogData = .FALSE.
            end do
         end if
         if (needPRTable) then
            theGphTable%tg(EOS_TABVT_PR)%td%ntemp = nstepsTemperatureHC !DEV: ??
            theGphTable%tg(EOS_TABVT_PR)%td%ndens = nstepsDensityHC !DEV: ??
            nullify(theGphTable%tg(EOS_TABVT_PR)%td%Temperatures)
            nullify(theGphTable%tg(EOS_TABVT_PR)%td%Densities)
!!$            allocate(theGphTable%tg(EOS_TABVT_PR)%td%Temperatures(nstepsTemperatureHC))
!!$            allocate(theGphTable%tg(EOS_TABVT_PR)%td%Densities(nstepsDensityHC))

!!$            allocate(theGphTable%tg(EOS_TABVT_PR)%table(EOS_TAB_FOR_ION:EOS_TAB_FOR_ELE))
            if (isGpNLike) then
               allocate(theGphTable%tg(EOS_TABVT_PR)%table(EOS_TAB_NCOMP))
               theGphTable%tg(EOS_TABVT_PR)%numTables = EOS_TAB_NCOMP
            else
               allocate(theGphTable%tg(EOS_TABVT_PR)%table(EOS_TAB_FOR_MAT:EOS_TAB_FOR_MAT))
               theGphTable%tg(EOS_TABVT_PR)%numTables = 1 
            end if
!!$            theGphTable%tg(EOS_TABVT_PR)%numTables = EOS_TAB_FOR_ELE - EOS_TAB_FOR_ION + 1
            do i = LBOUND(theGphTable%tg(EOS_TABVT_PR)%table,1), &
                 UBOUND(theGphTable%tg(EOS_TABVT_PR)%table,1)
               nullify(theGphTable%tg(EOS_TABVT_PR)%table(i)%table)
               theGphTable%tg(EOS_TABVT_PR)%table(i)%isLogData = .FALSE.
            end do
         end if
         if (needEntrTable) then
            theGphTable%tg(EOS_TABVT_ENTR)%td%ntemp = nstepsTemperatureEntr
            theGphTable%tg(EOS_TABVT_ENTR)%td%ndens = nstepsDensityEntr
            nullify(theGphTable%tg(EOS_TABVT_ENTR)%td%Temperatures)
            nullify(theGphTable%tg(EOS_TABVT_ENTR)%td%Densities)
!!$            allocate(theGphTable%tg(EOS_TABVT_ENTR)%td%Temperatures(nstepsTemperatureEntr))
!!$            allocate(theGphTable%tg(EOS_TABVT_ENTR)%td%Densities(nstepsDensityEntr))

!!$            allocate(theGphTable%tg(EOS_TABVT_ENTR)%table(EOS_TAB_FOR_ELE))
            allocate(theGphTable%tg(EOS_TABVT_ENTR)%table(EOS_TAB_FOR_ELE:EOS_TAB_FOR_ELE))
            theGphTable%tg(EOS_TABVT_ENTR)%numTables = 1
            do i = LBOUND(theGphTable%tg(EOS_TABVT_ENTR)%table,1), &
                 UBOUND(theGphTable%tg(EOS_TABVT_ENTR)%table,1)
               nullify(theGphTable%tg(EOS_TABVT_ENTR)%table(i)%table)
               theGphTable%tg(EOS_TABVT_ENTR)%table(i)%isLogData = .FALSE.
            end do
         end if

     end if
#endif

!
!
!    ...Close the 'EOS_sources.txt' file.
!
!
!
!
  ! PRE-LOOP # 3  (some checks)
  !
!
!
  if (op_maxTablesZF > 0) then

      if (op_maxNstepsDensityZF == 0) then
          call Driver_abortFlash ('[eos_initTabulated] ERROR: no ZF table density grid')
      end if
      if (op_maxNstepsTemperatureZF == 0) then
          call Driver_abortFlash ('[eos_initTabulated] ERROR: no ZF table temperature grid')
      end if

  end if

  if (op_maxTablesEN > 0) then

      if (op_maxNstepsDensityEN == 0) then
          call Driver_abortFlash ('[eos_initTabulated] ERROR: no EN table density grid')
      end if
      if (op_maxNstepsTemperatureEN == 0) then
          call Driver_abortFlash ('[eos_initTabulated] ERROR: no EN table temperature grid')
      end if

  end if

  if (op_maxTablesHC > 0) then

      if (op_maxNstepsDensityHC == 0) then
          call Driver_abortFlash ('[eos_initTabulated] ERROR: no HC table density grid')
      end if
      if (op_maxNstepsTemperatureHC == 0) then
          call Driver_abortFlash ('[eos_initTabulated] ERROR: no HC table temperature grid')
      end if

  end if

!
!
!    ...Read the tables and store the needed opacity values and all associated data
!       into the arrays.
!
!

  ! noLOOP # 3 - Actually call eos_gphReadTables with appropriate arguments
  ! to get the table data that will be needed into memory.
  !

     wanted = .FALSE.
     isOpacplot = (eos_gphKind=='OPACPLOT')
     isGpN = (eos_gphKind=='GPN')
     isIonmix6 = (eos_gphKind=='IONMIX6')
     isGpNLike = (isGpN .OR. isIonmix6)

     needEntrTable = .TRUE.

     needTable    =       needZFTable &
                    .or.  needENTable &
                    .or.  needHCTable &
                    .or.  needPRTable &
                    .or.  needEntrTable

#ifdef SUPPRESS_OLD4
     if (needZFTable) then
         wanted(EOS_TAB_FOR_ELE,EOS_TABVT_ZF) = .TRUE.
     end if

     if (needENTable) then
        wanted(EOS_TAB_FOR_MAT,EOS_TABVT_EN) = .TRUE.
     end if


     nullify(temperatures)
     nullify(densities)
#endif
     if (needEntrTable) then
        wanted(EOS_TAB_FOR_MAT,EOS_TABVT_ENTR) = .TRUE.
     end if

     if (needTable) then

        if(eos_meshMe < 4) then
           print *, "in eos_inittabulated, tableName = ", trim(eos_gphFileName)
        end if

!!$        call eos_gphReadTables (eos_gphKind, &
!!$                                eos_gphFileName, &
!!$                                wanted,                 &
!!$                                theGphTable%tg(:)%td, &
!!$                                theGphTable%tg(EOS_TABVT_ENTR)%table(EOS_TAB_FOR_MAT))
        call eos_gphReadGpNTables (eos_gphFileName, &
                                wanted,                 &
                                theGphTable%tg(:)%td, &
                                theGphTable%tg(EOS_TABVT_ENTR)%table)
     end if

!
!
!    DEV: ?? ...Print out the EOS data constants to see what has been stored.
!
!
#if(0)
  fileUnit = ut_getFreeFileUnit ()

  open (unit = fileUnit, &
        file = "opacity_printout_constants.txt", &
        form = 'formatted')

  printoutHeader = "   EOS_APPROX_KIN CONST PRINTOUT (in ?some? units)"

  call eos_tabWriteConstants (fileUnit,printoutHeader)

  close (fileUnit)
#endif
!
!
!    ...Print out the data tables to see what has been stored.
!
!

#if 0
  if ( eos_globalMe == MASTER_PE ) then

     fileUnit = ut_getFreeFileUnit ()

     open (unit = fileUnit, &
          file = "EOS_printout_tables.txt", &
          form = 'formatted')

     if (eos_useLogTables) then
        printoutHeader = "   EOS TABLES PRINTOUT (logarithmic base 10, energies in units of erg/g)"
     else
        printoutHeader = "   EOS TABLES PRINTOUT (energies in units of erg/g)"
     end if

     call eos_tabWriteTables (fileUnit,printoutHeader)


     close (fileUnit)

  end if
#endif

  ! LOOP # 4 - Initialize for EOS table diagnostics.
  !
  allocate(eos_gphAllDiag(EOS_GPH_NALLTAB))
  do tabno = 1,EOS_GPH_NALLTAB
     eos_gphAllDiag(tabno)%highTempCount = 0
     eos_gphAllDiag(tabno)%highDensCount = 0
     eos_gphAllDiag(tabno)%highestTemp = -999.0
     eos_gphAllDiag(tabno)%highestDens = -999.0
     eos_gphAllDiag(tabno)%highTempVarsLookedUp(:) = .FALSE.
     eos_gphAllDiag(tabno)%highDensVarsLookedUp(:) = .FALSE.
     eos_gphAllDiag(tabno)%firstHighTempEvent%temp = -999.0
     eos_gphAllDiag(tabno)%firstHighTempEvent%dens = -999.0
     eos_gphAllDiag(tabno)%firstHighDensEvent%temp = -999.0
     eos_gphAllDiag(tabno)%firstHighDensEvent%dens = -999.0
  end do
  
!
!
!    ...Ready!
!
!
  return
end subroutine eos_initTabulated
