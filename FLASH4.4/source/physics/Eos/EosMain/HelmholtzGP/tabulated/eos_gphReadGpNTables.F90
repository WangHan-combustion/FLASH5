!!****if* source/physics/Eos/EosMain/Tabulated/eos_gphReadIonmix4Tables
!!
!! NAME
!!
!!  eos_gphReadGpNTables
!!
!! SYNOPSIS
!!
!!  call eos_gphReadGpNTables (character (in) :: tableName (len=80),
!!                            logical   (in) :: needZFTable,
!!                            logical   (in) :: needENTable,
!!                            logical   (in) :: needHCTable)
!!
!! DESCRIPTION
!!
!!  Reads tabulated opacities from an GPN datafile output. The tabulated opacities
!!  will be stored into the 4-dimensional arrays: ! DEV: removed !
!!
!!             eos_gphIonizationTables  (t,d,g,indexZF)  (in cm^2/g)
!!             eos_gphIntEnergyTables    (t,d,g,indexEN)  (in cm^2/g)
!!             eos_gphHeatCpTables        (t,d,g,indexHC)  (in cm^2/g)
!!
!!    where:   t = temperature (eV) index
!!             d = ion number density (# ions/cm^3) index
!!             g = energy group (eV) index
!!       indexXX = table counting index for Opacity kind XX (XX = ZF,EN,HC)
!!
!!  The number of temperature and density indices will be stored in the 2-dimensional arrays: 
!!
!!         op_nstepsDensity (EOS_TABULAR_Z,indexZF)  !!DEV: outdated! - KW
!!         op_nstepsDensity (EOS_TABULAR_E,indexEN)
!!         op_nstepsDensity (EOS_TABULAR_C,indexHC)
!!
!!         op_nstepsTemperature (EOS_TABULAR_Z,indexZF)
!!         op_nstepsTemperature (EOS_TABULAR_E,indexEN)
!!         op_nstepsTemperature (EOS_TABULAR_C,indexHC)
!!
!!  The actual temperatures and densities will be stored in the 2-dimensional arrays: !!DEV: removed !
!!
!!         op_plasmaDensityZF     (d,indexZF) ; d = 1,op_nstepsDensity (EOS_TABULAR_Z,indexZF)       (in # ions/cm^3)
!!         op_plasmaDensityEN     (d,indexEN) ; d = 1,op_nstepsDensity (EOS_TABULAR_E,indexEN)       (in # ions/cm^3)
!!         op_plasmaDensityHC     (d,indexHC) ; d = 1,op_nstepsDensity (EOS_TABULAR_C,indexHC)       (in # ions/cm^3)
!!
!!         op_plasmaTemperatureZF (t,indexZF) ; t = 1,op_nstepsTemperature (EOS_TABULAR_Z,indexZF)   (in eV, 1eV = 11405K)
!!         op_plasmaTemperatureEN (t,indexEN) ; t = 1,op_nstepsTemperature (EOS_TABULAR_E,indexEN)   (in eV, 1eV = 11405K)
!!         op_plasmaTemperatureHC (t,indexHC) ; t = 1,op_nstepsTemperature (EOS_TABULAR_C,indexHC)   (in eV, 1eV = 11405K)
!!
!!
!! ARGUMENTS
!!
!!  tableName   : the name of the GPN file
!!  indexZF     : table counting index where average ionization data will be placed
!!  indexEN     : table counting index where internal energy data will be placed
!!  indexHC     : table counting index where specific heat data will be placed
!!
!!***
subroutine eos_gphReadGpNTables (tableName,   &
                                wanted, &
                                td,&
!!$                                tbZF,tbEN,tbPR,tbHC,&
                                tbEntr)

  use Eos_data,         ONLY : eos_smallT
  use eos_gphData,      ONLY :  eos_useLogTables,          &
                                EOS_TAB_NCOMP,         &
                                EOS_GPH_NALLTAB,         &
!!$                                EOS_TAB_FOR_ION,         &
!!$                                EOS_TAB_FOR_ELE,         &
                                EOS_TAB_FOR_MAT,         &
!!$                                EOS_TABVT_ZF,         &
!!$                                EOS_TABVT_EN,         &
!!$                                EOS_TABVT_PR,         &
!!$                                EOS_TABVT_HC,         &
                                EOS_TABVT_ENTR,       &
                                eosT_tableIvarsetDescT, &
                                eosT_oneVarTablePT

  use Driver_interface,  ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_GPH_NALLTAB)
  type(eosT_tableIvarsetDescT),intent(inout) :: td(:)
!!$  type(eosT_oneVarTablePT),pointer,dimension(:) :: tbZF,tbEN,tbPR,tbHC,tbEntr
  type(eosT_oneVarTablePT),intent(INOUT),dimension(:) :: tbEntr
  character (len=80), intent (in) :: tableName

  character (len=80) :: dummyLine

!!  needZFTable : if yes, average ionization data are needed from the GPN table
!!  needENTable : if yes, internal energy data are needed from the GPN table
!!  needHCTable : if yes,        specific heat data are needed from the GPN table
  logical :: needZFTable
  logical :: needENTables
  logical :: needPRTables
  logical :: needHCTables
  logical :: needEntrTables
  logical :: fileExists, doread

  integer :: tableDim, ncorners, nderivs
  integer,dimension(EOST_MAX_IVARS) :: tdims

  integer :: ntemp, ndens
  real, allocatable :: temperatures(:)
  real, allocatable :: densities(:)
  integer :: fileUnit
  integer :: notneededData
  integer :: step
  integer :: i,n,t,d,g
  integer :: iw,ix,iy,iz
  integer :: ngroupsEnergy
  integer :: nstepsDensity
  integer :: nstepsTemperature
  integer :: ut_getFreeFileUnit

  real    :: dummyData
!!$  real    :: log10Density
  real    :: log10DensityStep
  real    :: log10DensityFirst
!!$  real    :: log10Temperature
  real    :: log10TemperatureStep
  real    :: log10TemperatureFirst
  real    :: maxTempStepFactor,maxDensStepFactor
  real    :: minTempStepFactor,minDensStepFactor

  ! Conversion from eV to K:
  real, parameter :: K = 11604.5221
  real, parameter :: joule2erg = 1.0e7
!
!
!   ...Check and open the opacity file.
!
!
  inquire (file = tableName , exist = fileExists)

  if (.not.fileExists) then
       call Driver_abortFlash ('[eos_gphReadGpNTables] ERROR: no GPN file found')
  end if

  fileUnit = ut_getFreeFileUnit ()
  open (unit = fileUnit , file = tableName)
!
!
!   ...Read the temperature, density and energy group grids. Abort the calculation,
!      if any of the grids is not found.
!
!
  read (fileUnit,'(A80)')  dummyLine
  print*,'Read "',trim(dummyLine),'" from ',trim(tableName)
  read(fileUnit,*) tableDim, nderivs
  ncorners      = 2**tableDim

#ifdef SUPPRESS_OLD4
  read (fileUnit,'(2I10)') nstepsTemperature , nstepsDensity
  ntemp = nstepsTemperature
  ndens = nstepsDensity
  read (fileUnit,'(A80)')  dummyLine
  read (fileUnit,'(A80)')  dummyLine

!!$  read (fileUnit,'(4E12.6,I12)') log10DensityStep,       &
!!$                                 log10DensityFirst,      &
!!$                                 log10TemperatureStep,   &
!!$                                 log10TemperatureFirst,  &
!!$                                 ngroupsEnergy
  read (fileUnit,'(I12)') ngroupsEnergy

  if (nstepsTemperature <= 0) then
      call Driver_abortFlash ('[eos_gphReadGpNTables] ERROR: no GPN temperature grid found')
  end if

  if (nstepsDensity <= 0) then
      call Driver_abortFlash ('[eos_gphReadGpNTables] ERROR: no GPN density grid found')
  end if

  if (ngroupsEnergy <= 0) then
      call Logfile_stampMessage('[eos_gphReadGpNTables] WARNING: no GPN energy group grid found')
  end if


  allocate(temperatures(ntemp))
  allocate(densities(ndens))

  read (fileUnit,'(4E12.6)') (temperatures(t),t=1,ntemp)
  ! Convert the temperatures into K from eV:
  temperatures(:) = temperatures(:) * K

  log10TemperatureFirst = log10(Temperatures(1))
  minTempStepFactor = HUGE(minTempStepFactor)
  maxTempStepFactor = 1
  do t = 1,ntemp-1
     if(temperatures(t) < eos_smallT) then 
        if(temperatures(t+1) > eos_smallT) &
             minTempStepFactor = min(minTempStepFactor,temperatures(t+1)/eos_smallT)
        maxTempStepFactor = max(maxTempStepFactor,temperatures(t+1)/eos_smallT)
     else
        minTempStepFactor = min(minTempStepFactor,temperatures(t+1)/temperatures(t))
        maxTempStepFactor = max(maxTempStepFactor,temperatures(t+1)/temperatures(t))
     end if
  end do


  read (fileUnit,'(4E12.6)') (densities(d),d=1,ndens)

  log10DensityFirst = log10(densities(1))
  minDensStepFactor = HUGE(minDensStepFactor)
  maxDensStepFactor = 1
  do d = 1,ndens-1
     minDensStepFactor = min(maxDensStepFactor,densities(d+1)/densities(d))
     maxDensStepFactor = max(maxDensStepFactor,densities(d+1)/densities(d))
  end do



  needZFTable = wanted(EOS_TAB_FOR_ELE,EOS_TABVT_ZF)
  needENTables = ANY(wanted(EOS_TAB_FOR_ION:EOS_TAB_FOR_ELE,EOS_TABVT_EN))
  needPRTables = ANY(wanted(EOS_TAB_FOR_ION:EOS_TAB_FOR_ELE,EOS_TABVT_PR))
  needHCTables = ANY(wanted(:,EOS_TABVT_HC))
#endif
  needEntrTables = ANY(wanted(EOS_TAB_FOR_MAT:EOS_TAB_FOR_MAT,EOS_TABVT_ENTR))
  print*,'wanted:',wanted
  print*,'wanted(EOS_TAB_FOR_MAT:EOS_TAB_FOR_MAT,EOS_TABVT_ENTR):',wanted(EOS_TAB_FOR_MAT:EOS_TAB_FOR_MAT,EOS_TABVT_ENTR)
#ifdef SUPPRESS_OLD4
!
!
!   ...Establish the temperature and density grid for the average ionization (if needed)
!      and read in the data.
!
!

  if (needZFTable) then

      td(EOS_TABVT_ZF)%ntemp = nstepsTemperature
      td(EOS_TABVT_ZF)%ndens = nstepsDensity
      if (.NOT.associated(td(EOS_TABVT_ZF)%Temperatures)) then
         allocate(td(EOS_TABVT_ZF)%Temperatures(ntemp)) 
      end if
      if (.NOT.associated(td(EOS_TABVT_ZF)%Densities)) then
         allocate(td(EOS_TABVT_ZF)%Densities(ndens))
      end if

      do step = 1,nstepsDensity
         td(EOS_TABVT_ZF)%Densities(step) = densities(step)
      end do

      do step = 1,nstepsTemperature
         td(EOS_TABVT_ZF)%Temperatures(step) = temperatures(step)
      end do

!!$     print*,'LBOUND(tbZF,1):',LBOUND(tbZF,1)
!!$     print*,'UBOUND(tbZF,1):',UBOUND(tbZF,1)
      allocate(tbZF(EOS_TAB_FOR_ELE)%table(nstepsTemperature,nstepsDensity))
      read (fileUnit,'(4E12.6)') ((tbZF(EOS_TAB_FOR_ELE)%table(t,d), t = 1,nstepsTemperature), d = 1,nstepsDensity    )

      if (td(EOS_TABVT_ZF)%isLog) then
         td(EOS_TABVT_ZF)%Densities    (1:nstepsDensity)     = log10(td(EOS_TABVT_ZF)%Densities    (1:nstepsDensity) )
         td(EOS_TABVT_ZF)%Temperatures (1:nstepsTemperature) = log10(td(EOS_TABVT_ZF)%Temperatures (1:nstepsTemperature) )
         tbZF(EOS_TAB_FOR_ELE)%isLogData = .TRUE.
      end if
      if (tbZF(EOS_TAB_FOR_ELE)%isLogData) then
         tbZF(EOS_TAB_FOR_ELE)%table(1:nstepsTemperature,1:nstepsDensity) = &
              log10(tbZF(EOS_TAB_FOR_ELE)%table(1:nstepsTemperature,1:nstepsDensity))
      end if

  else
      !   ...Skip not needed data from the GPN file.
      notneededData =   nstepsDensity * nstepsTemperature
      read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)
  end if

  ! ...Skip dz/dT table
  notneededData =   nstepsDensity * nstepsTemperature
  read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)
!
!
!   ...Establish the temperature and density grid for pressures (if needed)
!      and read in the data.
!
!
  if (needPRTables) then

      td(EOS_TABVT_PR)%ntemp = nstepsTemperature
      td(EOS_TABVT_PR)%ndens = nstepsDensity
      if (.NOT.associated(td(EOS_TABVT_PR)%Temperatures)) then
         allocate(td(EOS_TABVT_PR)%Temperatures(ntemp)) 
      end if
      if (.NOT.associated(td(EOS_TABVT_PR)%Densities)) then
         allocate(td(EOS_TABVT_PR)%Densities(ndens))
      end if

      do step = 1,nstepsDensity
         td(EOS_TABVT_PR)%Densities(step) = densities(step)
      end do

      do step = 1,nstepsTemperature
         td(EOS_TABVT_PR)%Temperatures(step) = temperatures(step)
      end do

      if (td(EOS_TABVT_PR)%isLog) then
         td(EOS_TABVT_PR)%Densities    (1:nstepsDensity)     = log10(td(EOS_TABVT_PR)%Densities    (1:nstepsDensity) )
         td(EOS_TABVT_PR)%Temperatures (1:nstepsTemperature) = log10(td(EOS_TABVT_PR)%Temperatures (1:nstepsTemperature) )
         do n = 1,EOS_TAB_NCOMP
            if (wanted(n,EOS_TABVT_PR)) tbPR(n)%isLogData = .TRUE.
         end do
      end if

      do n=EOS_TAB_FOR_ION,EOS_TAB_FOR_ION+1
         doread = (n .LE. EOS_TAB_NCOMP)
         if (doread) doread = wanted(n,EOS_TABVT_PR)
         if (doread) then
            allocate(tbPR(n)%table(nstepsTemperature,nstepsDensity))
            read (fileUnit,'(4E12.6)') ((tbPR(n)%table(t,d), t = 1,nstepsTemperature), d = 1,nstepsDensity    )

            tbPR(n)%table(1:nstepsTemperature,1:nstepsDensity) = &
                 tbPR(n)%table(1:nstepsTemperature,1:nstepsDensity) * joule2erg
            if (tbPR(n)%isLogData) then
               tbPR(n)%table(1:nstepsTemperature,1:nstepsDensity) = &
                    log10(tbPR(n)%table(1:nstepsTemperature,1:nstepsDensity))
            end if
         else
            !   ...Skip not needed data from the GPN file.
            notneededData =   nstepsDensity * nstepsTemperature
            read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)
         end if
      end do
  else
      !   ...Skip not needed data from the GPN file.
      notneededData =   nstepsDensity * nstepsTemperature
      read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)
      read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)
  end if

  ! ...Skip dPion/dTion table
  notneededData =   nstepsDensity * nstepsTemperature
  read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)

  ! ...Skip dPele/dTele table
  notneededData =   nstepsDensity * nstepsTemperature
  read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)

!
!
!   ...Establish the temperature and density grid for the internal energy (if needed)
!      and read in the data.
!
!
  if (needENTables) then

      td(EOS_TABVT_EN)%ntemp = nstepsTemperature
      td(EOS_TABVT_EN)%ndens = nstepsDensity
      if (.NOT.associated(td(EOS_TABVT_EN)%Temperatures)) then
         allocate(td(EOS_TABVT_EN)%Temperatures(ntemp)) 
      end if
      if (.NOT.associated(td(EOS_TABVT_EN)%Densities)) then
         allocate(td(EOS_TABVT_EN)%Densities(ndens))
      end if

      do step = 1,nstepsDensity
         td(EOS_TABVT_EN)%Densities(step) = densities(step)
      end do

      do step = 1,nstepsTemperature
         td(EOS_TABVT_EN)%Temperatures(step) = temperatures(step)
      end do

      if (td(EOS_TABVT_EN)%isLog) then
         td(EOS_TABVT_EN)%Densities    (1:nstepsDensity)     = log10(td(EOS_TABVT_EN)%Densities    (1:nstepsDensity) )
         td(EOS_TABVT_EN)%Temperatures (1:nstepsTemperature) = log10(td(EOS_TABVT_EN)%Temperatures (1:nstepsTemperature) )
         do n = 1,EOS_TAB_NCOMP
            if (wanted(n,EOS_TABVT_EN)) tbEN(n)%isLogData = .TRUE.
         end do
      end if

      do n=EOS_TAB_FOR_ION,EOS_TAB_FOR_ION+1
         doread = (n .LE. EOS_TAB_NCOMP)
         if (doread) doread = wanted(n,EOS_TABVT_EN)
         if (doread) then
            allocate(tbEN(n)%table(nstepsTemperature,nstepsDensity))
            read (fileUnit,'(4E12.6)') ((tbEN(n)%table(t,d), t = 1,nstepsTemperature), d = 1,nstepsDensity    )

            tbEN(n)%table(1:nstepsTemperature,1:nstepsDensity) = &
                 tbEN(n)%table(1:nstepsTemperature,1:nstepsDensity) * joule2erg
            if (tbEN(n)%isLogData) then
               tbEN(n)%table(1:nstepsTemperature,1:nstepsDensity) = &
                    log10(tbEN(n)%table(1:nstepsTemperature,1:nstepsDensity))
            end if
         else
            !   ...Skip not needed data from the GPN file.
            notneededData =   nstepsDensity * nstepsTemperature
            read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)
         end if
      end do
  else
      !   ...Skip not needed data from the GPN file.
      notneededData =   nstepsDensity * nstepsTemperature
      read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)
      read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)
  end if
!
!
!   ...Establish the temperature and density grid for the specific heat capacity (if needed)
!      and read in the data.
!
!
  if (needHCTables) then

      td(EOS_TABVT_HC)%ntemp = nstepsTemperature
      td(EOS_TABVT_HC)%ndens = nstepsDensity
      if (.NOT.associated(td(EOS_TABVT_HC)%Temperatures)) then
         allocate(td(EOS_TABVT_HC)%Temperatures(ntemp)) 
      end if
      if (.NOT.associated(td(EOS_TABVT_HC)%Densities)) then
         allocate(td(EOS_TABVT_HC)%Densities(ndens))
      end if

      do step = 1,nstepsDensity
         td(EOS_TABVT_HC)%Densities(step) = densities(step)
      end do

      do step = 1,nstepsTemperature
         td(EOS_TABVT_HC)%Temperatures(step) = temperatures(step)
      end do

      if (td(EOS_TABVT_HC)%isLog) then
         td(EOS_TABVT_HC)%Densities    (1:nstepsDensity)     = log10(td(EOS_TABVT_HC)%Densities    (1:nstepsDensity) )
         td(EOS_TABVT_HC)%Temperatures (1:nstepsTemperature) = log10(td(EOS_TABVT_HC)%Temperatures (1:nstepsTemperature) )
         do n = 1,EOS_TAB_NCOMP
            if (wanted(n,EOS_TABVT_HC)) tbHC(n)%isLogData = .TRUE.
         end do
      end if

      do n=EOS_TAB_FOR_ION,EOS_TAB_FOR_ION+1
         doread = (n .LE. EOS_TAB_NCOMP)
         if (doread) doread = wanted(n,EOS_TABVT_HC)
         if (doread) then
            allocate(tbHC(n)%table(nstepsTemperature,nstepsDensity))
            read (fileUnit,'(4E12.6)') ((tbHC(n)%table(t,d), t = 1,nstepsTemperature), d = 1,nstepsDensity    )

            tbHC(n)%table(1:nstepsTemperature,1:nstepsDensity) = &
                 tbHC(n)%table(1:nstepsTemperature,1:nstepsDensity) * joule2erg / K
            if (tbHC(n)%isLogData) then
               tbHC(n)%table(1:nstepsTemperature,1:nstepsDensity) = &
                    log10(tbHC(n)%table(1:nstepsTemperature,1:nstepsDensity))
            end if
         else
            !   ...Skip not needed data from the GPN file.
            notneededData =   nstepsDensity * nstepsTemperature
            read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)
         end if
      end do

  else
      !   ...Skip not needed data from the GPN / IONMIX6 file.
      notneededData =   nstepsDensity * nstepsTemperature
      read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)
      read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)

  end if

  ! ...Skip dEion/dNion table
  notneededData =   nstepsDensity * nstepsTemperature
  read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)

  ! ...Skip dEele/dNion table
  notneededData =   nstepsDensity * nstepsTemperature
  read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)

#endif
!
!
!   ...Establish the temperature and density grid for the specific (electron) entropy (if needed)
!      and read in the data.  This is onlyfor IONMIX6, not GPN.
!
!
  if (needEntrTables) then

    n = EOS_TAB_FOR_MAT
    associate (tdS => td(EOS_TABVT_ENTR), &
               tbS => tbEntr(n))

      tdS % N       = tableDim
      tdS % nderivs = nderivs
      tdS % ncorners = ncorners

      tdims(:) = 1
      read(fileUnit,*) (tdS % c(i) % nIval, i=1,tableDim)
      do i=1,tableDim
         tdims(i) = tdS % c(i) % nIval
      end do
      read(fileUnit,*)

      do i=1,tableDim
         allocate(tdS % c(i) % val(tdims(i)+1))
         read(fileUnit,*)  tdS % c(i) % nIval, notneededData, &
                           tdS % c(i) % val 
      end do

      allocate(tbS % derDefs(tableDim,0:nderivs))
      tbS % derDefs(:,0) = 0    ! for the 'zeroth' derivative
      do i=1,nderivs
         read(fileUnit,*)  tbS % derDefs(:,i)
      end do

      allocate(tbS % table(0:ncorners-1, 0:nderivs, tdims(1), tdims(2), tdims(3), tdims(4)))
      allocate(tbS % ells (tableDim,                tdims(1), tdims(2), tdims(3), tdims(4)))
      do iz=1,tdims(4)
         do iy=1,tdims(3)
            do ix=1,tdims(2)
               do iw=1,tdims(1)
                  read(fileUnit,*)  tbS % table(0:,0:,iw,ix,iy,iz), tbS % ells(:,iw,ix,iy,iz)
               end do
            end do
         end do
      end do
      print*,'ells:',tbS % ells
         

!!$      if (tdS % isLog) then
!!$         tdS % c(2)%val    (1:nstepsDensity)     = log10(tdS % c(2)%val    (1:nstepsDensity) )
!!$         tdS % c(1)%val (1:nstepsTemperature) = log10(tdS % c(1)%val (1:nstepsTemperature) )
!!$         do n = 1,EOS_TAB_NCOMP
!!$            if (wanted(n,EOS_TABVT_ENTR)) tbEntr(n)%isLogData = .TRUE.
!!$         end do
!!$      end if

#ifdef SUPPRESS_OLD4
         doread = (n .LE. EOS_TAB_NCOMP)
         if (doread) doread = wanted(n,EOS_TABVT_ENTR)
         if (doread) then
            allocate(tbS (nstepsTemperature,nstepsDensity))
            read (fileUnit,'(4E12.6)') ((tbS (t,d), t = 1,nstepsTemperature), d = 1,nstepsDensity    )

            tbS (1:nstepsTemperature,1:nstepsDensity) = &
                 tbS (1:nstepsTemperature,1:nstepsDensity) * joule2erg / K
!!$            if (tbEntr(n)%isLogData) then
!!$               tbS (1:nstepsTemperature,1:nstepsDensity) = &
!!$                    log10(tbS (1:nstepsTemperature,1:nstepsDensity))
!!$            end if
         else
            !   ...Skip not needed data from the GPN file.
            notneededData =   nstepsDensity * nstepsTemperature
            read (fileUnit,'(4E12.6)') (dummyData, i = 1,notneededData)
         end if
#endif
        end associate
      !
      !   ...Stop reading here, we do not care about the remainder of the file.
      !
      !


  end if

#ifdef SUPPRESS_OLD4
  deallocate(temperatures)
  deallocate(densities)
#endif

!
!
!   ...Close the GPN file.
!
!
  close (fileUnit)
!
!
!   ...Ready! 
!
!
  return
end subroutine eos_gphReadGpNTables
