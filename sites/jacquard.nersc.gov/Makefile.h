# FLASH makefile definitions for NERSC Jacquard Opteron Cluster (PathScale compiler)
#
# module load hdf5_par
#----------------------------------------------------------------------------
# Set the HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HDF5_PATH   = /usr/common/usg/hdf5/1.6.4/parallel
MPI_PATH    = 


PAPI_PATH  = 
PAPI_FLAGS = -c -I$(PAPI_PATH)/include


ZLIB_PATH  =
SZIP_PATH = /usr/common/usg/szip/2.0/lib

NCMPI_PATH = /usr/local/pnetcdf-1.0.0
MPE_PATH   =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP      = mpif90
CCOMP      = mpicc
CPPCOMP    = mpicxx
LINK       = mpif90




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

FFLAGS_OPT   = -c -r8 -i4 
FFLAGS_DEBUG = -c -r8 -i4 -g 
FFLAGS_TEST  = -c -r8 -i4 

CFLAGS_OPT   = -c -O2
CFLAGS_DEBUG = -c -g
CFLAGS_TEST  = -c 

FFLAGS_PAPI  = -I$(PAPI_PATH)/include
FFLAGS_MPI   = -I$(MPI_PATH)/include

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5  = -I$(HDF5_PATH)/include
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include
CFLAGS_MPI   = -I$(MPI_PATH)/include

CFLAGS_SAMRAI = -I$(SAMRAI_PATH)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o
LFLAGS_DEBUG = -g -o
LFLAGS_TEST  = -Mprof=lines -o

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


LIB_HDF5    = -L$(HDF5_PATH)/lib -lhdf5 -lz  -lpthread \
              -L$(SZIP_PATH) -lsz

LIB_MPI     = 
LIB_PAPI    = $(PAPI_PATH)/lib/libpapi.a $(PAPI_PATH)/lib/_fixunssfdi.o
LIB_PNG     = -lpng -lz

LIB_OPT     = 
LIB_DEBUG   =
LIB_TEST    =

LIB_SAMRAI  = -L$(SAMRAI_PATH)/lib -lSAMRAI



LIB_FISHPAK = -L$(FISHPAK_PATH)/lib -lfishpak
LIB_NCMPI   = -L$(NCMPI_PATH)/lib -lpnetcdf
LIB_MPE     =
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

MV    = mv -f
AR    = ar -r
RM    = rm -f
CD    = cd
RL    = ranlib
ECHO  = echo
AWK   = awk
CAT   = cat


