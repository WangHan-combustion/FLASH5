#-------------------------------------------------------------------
# FLASH makefile definitions for the bnl bgp
# 
#
# These modules set environment variables pointing to the
# proper library locations. Only python is essential.
#-------------------------------------------------------------------

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

# This is HDF5 1.8.x, so the necessary define for the 1.6 API is included in
# the HDF5 flags below.

# HDF5_PATH = /gpfs/home2/ajackson/local/hdf5/current/bgl/serial
HDF5_PATH = /gpfs/home2/ajackson/local/hdf5/current/bgp/parallel
# HDF5_PATH = /bgl/apps/hdf5

# NCMPI_PATH = /usr/local/tools/parallel-netcdf/parallel-netcdf-1.0.0

ZLIB_PATH =

HPM_PATH =
PNG_PATH = /usr

MPI_PATH = /bgsys/drivers/ppcfloor/comm

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

ifdef PDTDIR
OPT=-optPDBFile=merged.pdb -optTauSelectFile="select.tau" -optReset="" -optVerbose
else
LIB_MPI = 
# LIB_MPI = -L$(MPI_PATH)/lib -lmpich.rts -lmsglayer.rts -ldevices.rts -lrts.rts -ldevices.rts
endif
FCOMP   = $(MPI_PATH)/bin/mpixlf90 -I$(MPI_PATH)/include
CCOMP   = $(MPI_PATH)/bin/mpixlc -I$(MPI_PATH)/include
CPPCOMP = $(MPI_PATH)/bin/mpixlcxx
LINK    = $(MPI_PATH)/bin/mpixlf90

#CONFIG_LIB = -L$(MPI_PATH)/lib -lmpich.rts -lmsglayer.rts -ldevices.rts -lrts.rts -ldevices.rts
#-----------------------------------------------------------------------------
# Compilation flags
#
#  Three sets of compilation/linking flags are defined: one for optimized code
#  code ("-opt"), one for debugging ("-debug"), and one for testing ("-test").
#  Passing these flags to the setup script will cause the value associated with
#  the corresponding keys (i.e. those ending in "_OPT", "_DEBUG", or "_TEST") to
#  be incorporated into the final Makefile. For example, passing "-opt" to the
#  setup script will cause the flags following "FFLAGS_OPT" to be assigned to
#  "FFLAGS" in the final Makefile. If none of these flags are passed, the default
#  behavior will match that of the "-opt" flag.
#  In general, "-opt" is meant to optimize compilation and linking. "-debug"
#  should enable runtime bounds checking, debugger symbols, and other compiler-
#  specific debugging options. "-test" is useful for testing different
#  combinations of compiler flags particular to your individual system.
#----------------------------------------------------------------------------

FFLAGS_OPT   = -O2 -qintsize=4 -qrealsize=8 -qdpc=e \
               -qfixed -qnosave -c \
               -qinline -qmaxmem=131072 \
               -qsuffix=cpp=F -qarch=450 -qtune=450 
               
FFLAGS_TEST  = -O2 -qintsize=4 -qrealsize=8 -qdpc=e \
               -qfixed -qnosave -c \
               -qinline -qmaxmem=131072 \
               -qsuffix=cpp=F -qarch=450 -qtune=450 
              
FFLAGS_DEBUG = -g -C -qintsize=4 -qrealsize=8 -qfixed -qnosave -c \
               -qinline -qmaxmem=131072 \
               -qsuffix=cpp=F -qarch=450 -qtune=450

F90FLAGS     = -qsuffix=f=F90:cpp=F90 -qfree=f90
f90FLAGS     = -qsuffix=f=f90:cpp=F90 -qfree=f90

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5  = -I${HDF5_PATH}/include -DH5_USE_16_API
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include

CFLAGS_OPT   = -O2 -DIBM -DNOUNDERSCORE -c -qlanglvl=extended \
               -qarch=450 -qtune=450 -qmaxmem=131072
CFLAGS_TEST  = -O3 -DIBM -DNOUNDERSCORE -c \
               -qarch=450 -qtune=450 
CFLAGS_DEBUG = -g -C -DIBM -DNOUNDERSCORE -c -qlanglvl=extended\
               -qarch=450 -qtune=450 -qmaxmem=131072

MDEFS = -WF,

.SUFFIXES: .o .c .f .F .h .fh .F90 .f90

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o
LFLAGS_TEST  = -o
LFLAGS_DEBUG = -g -o

#----------------------------------------------------------------------------
# Library specific linking
#
#  If a FLASH module has a 'LIBRARY xxx' line in its Config file, we need to
#  create a macro in this Makefile.h for LIB_xxx, which will be added to the
#  link line when FLASH is built.  This allows us to switch between different
#  (incompatible) libraries.  We also create a _OPT, _DEBUG, and _TEST
#  library macro to add any performance-minded libraries (like fast math),
#  depending on how FLASH was setup.
#----------------------------------------------------------------------------

LIB_HDF5  = -L$(HDF5_PATH)/lib -lhdf5
LIB_NCMPI = -L$(NCMPI_PATH)/lib -lpnetcdf




LIB_OPT   = 
LIB_DEBUG =
LIB_TEST  =

#----------------------------------------------------------------------------
# Additional machine-dependent object files
#
#  Add any machine specific files here -- they will be compiled and linked
#  when FLASH is built.
#----------------------------------------------------------------------------

MACHOBJ = 

#----------------------------------------------------------------------------
# Additional commands
#---------------------------------------------------------------------------- 

MV = mv -f
AR = ar -r
RM = rm -f
CD = cd
RL = ranlib
ECHO = echo

# full optimization does not work on runtime_parameters.F90
#runtime_parameters.o : runtime_parameters.F90
#	${FCOMP} ${FFLAGS_TEST}  $(F90FLAGS) $<

#amr_%.o : amr_%.F90
#	${FCOMP} ${FFLAGS_OPT_NOIPA} ${F90FLAGS} ${FDEFINES} $<
