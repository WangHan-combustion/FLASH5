!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!!  This is a modification of the standard Driver_evolveFlash for the unitTest for 
!!     the Particles unit.  Instead of using Hydro sweeps in xyz/zyx, the
!!     velocities are generated with a random perturbation overlaid on a 
!!     constant flow in the x-direction.
!!
!! NOTES
!!
!!
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe, dr_nstep
       
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
    Timers_getSummary
  use Gravity_interface, ONLY :Gravity_potential
  use IO_interface, ONLY :IO_writeCheckpoint, IO_writePlotfile
  implicit none

#include "constants.h"
#include "Flash.h"

  integer   :: localNumBlocks, i

  integer, parameter :: stepsPerAdvance = 2


  
  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr


  logical ::  perfect = .true.

  character(len=20) :: fileName
  integer, parameter        :: fileUnit = 2
  integer,dimension(4) :: prNum
  integer :: temp


  ! stays true if no errors are found
  perfect = .true.
  
  temp = dr_globalMe
  
  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do
  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
                                 char(48+prNum(2))//char(48+prNum(1))
  
  open(fileUnit,file=fileName)
  write(fileUnit,'("P",I0)') dr_globalMe



  call Logfile_stamp( 'Starting Calculation' , '[Driver_evolveFlash]')
  print*,'started calculation'
  call Timers_start("calculation")


  call Gravity_potential()
  print*,'got Potential'

  call Timers_stop("calculation")

  call Logfile_stamp( 'Ending Calculation' , '[Driver_evolveFlash]')

  call IO_writeCheckpoint()

  call IO_writePlotfile( )

  call Timers_getSummary( dr_nstep)

  
  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")

  call Logfile_close()


  !finish unit test write out file

  if (perfect) then
    write(fileUnit,'("all results conformed with expected values.")')
  endif

  close(fileUnit)



  return
  
end subroutine Driver_evolveFlash



