!!****if* source/IO/IOMain/hdf5/serial/AM/io_initFile
!!
!! NAME
!!  io_initFile
!!
!! SYNOPSIS
!!
!!  io_initFile(integer(in)      :: filenum,
!!              integer(io_fileID_t)(INOUT)   :: fileID,
!!              character(out)   :: filename(MAX_STRING_LENGTH),
!!              character(in)    :: outputType,
!!              logical(in)      :: forced)
!!
!! DESCRIPTION
!!
!!
!!  Initialized the hdf5 file for serial io
!!  only MASTER_PE opens the file
!!
!! ARGUMENTS
!!
!!
!!  filenum - number of the output file
!!
!!  fileID - file handle returned from hdf5 init file calls
!!
!!  filename - name of the returned file made up of the outputType, filenum and basename
!!
!!  outputType - string indicating output type, usually "chk" or "plt_cnt" or "part"
!!
!!  forced - .true. if this file is considered "forced."
!!
!!***


subroutine io_initFile( filenum, fileID, filename, outputType, forced)

  use IO_data, ONLY : io_comm, io_outputSplitNum, io_globalMe
  use Driver_interface, ONLY : Driver_abortFlash
  use io_intfTypesModule, ONLY : io_fileID_t

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: filenum
  integer(io_fileID_t), intent(INOUT) :: fileID
  character(len=*), intent(in) :: outputType
  character (len=MAX_STRING_LENGTH), intent(inout) :: filename
  logical, intent(in) :: forced


  if (io_globalMe /= MASTER_PE) then
     return
  end if

  call io_getOutputName(filenum, "hdf5", outputType, filename, forced)

  fileID = -1

  if (outputType(1:4) .NE. '_plt') then
     call io_h5init_file(fileID, filename)
     if(fileID == -1) then
        print *, "Error: unable to initialize file"
        call Driver_abortFlash("unable to initialize hdf5 file")
     end if
  end if

end subroutine io_initFile
