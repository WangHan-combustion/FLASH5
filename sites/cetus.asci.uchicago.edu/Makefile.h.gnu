# FLASH makefile definitions for x86-64 Linux (GNU compilers)
#
# Openmpi-1.3 has tight integration with Valgrind.
#   Run using:
#   export valgrind_log='valgrind --tool=memcheck --log-file=valgrind.log.%p \
#           --time-stamp=yes --track-origins=yes --run-libc-freeres=yes'   
#   mpiexec -n 4 ${valgrind_log} ./flash3
#   (Valgrind installed at /opt/valgrind-3.4.0)
#
# TAU makefiles configured for this software stack.
#   Setup FLASH using:
#   ./setup Sedov -auto \
#   -tau=/opt/tau-2.18.1/x86_64/lib/Makefile.tau-callpath-mpi-pdt -makefile=gnu
#   Analyse using:
#   export PATH=/usr/java/jdk1.6.0_11/bin/java:$PATH
#   /opt/tau-2.18.1/x86_64/bin/paraprof
#
#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------
HDF4_PATH =
HDF5_PATH = /opt/hdf5-1.8.5

ZLIB_PATH  =

PAPI_PATH  =
PAPI_FLAGS =

NCMPI_PATH = /opt/parallel-netcdf-1.2.0pre
MPE_PATH   =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#   Use the MPICH wrappers around the compilers -- these will automatically
#   load the proper libraries and include files.  Version of MPICH prior
#   to 1.2.2 (?) do not recognize .F90 as a valid Fortran file extension.
#   You need to edit mpif90 and add .F90 to the test of filename extensions,
#   or upgrade your MPICH.
#----------------------------------------------------------------------------

FCOMP   = /opt/openmpi-1.4.2/bin/mpif90
CCOMP   = /opt/openmpi-1.4.2/bin/mpicc
CPPCOMP = /opt/openmpi-1.4.2/bin/mpicxx
LINK    = /opt/openmpi-1.4.2/bin/mpif90

# pre-processor flag
PP      = -D

#----------------------------------------------------------------------------
# Compilation flags
#
#  Three sets of compilation/linking flags are defined: one for optimized
#  code, one for testing, and one for debugging.  The default is to use the
#  _OPT version.  Specifying -debug to setup will pick the _DEBUG version,
#  these should enable bounds checking.  Specifying _TEST is used for
#  flash_test, and is set for quick code generation, and (sometimes)
#  profiling.  The Makefile generated by setup will assign the generic token
#  (ex. FFLAGS) to the proper set of flags (ex. FFLAGS_OPT).
#----------------------------------------------------------------------------

FFLAGS_OPT = -g -c -O2 -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none -Wuninitialized

#Unfortunately gfortran is very out of date in the fedora repositories so 
#the following VERY USEFUL debugging options are not available:
#-fbacktrace -fdump-core -finit-real=nan
#Add to FFLAGS_DEBUG when fedora update gfortran to version>=4.3.3.

FFLAGS_DEBUG = -g -c -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none -pedantic -Wall -Wextra -Wconversion -Wunderflow \
-ffpe-trap=invalid,zero,overflow -fbounds-check

FFLAGS_TEST = -g -c -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none

F90FLAGS =


CFLAGS_OPT = -g -c -O2 -Wuninitialized

CFLAGS_DEBUG = -g -c -pedantic -Wall -Wextra -Winit-self -ftree-vrp \
-Wfloat-equal -Wunsafe-loop-optimizations -Wpadded

CFLAGS_TEST = -g -c


# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5 = -I ${HDF5_PATH}/include -DH5_USE_16_API
CFLAGS_NCMPI = -I ${NCMPI_PATH}/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -g -o

#Electric Fence (or eFence) is a memory debugger written by Bruce Perens.
#LFLAGS_DEBUG = -g -L/usr/lib64 -lefence -o

#mtrace is the memory debugger included in the GNU C Library.
#LFLAGS_DEBUG = -g -lmcheck -o

LFLAGS_DEBUG = -g -o
LFLAGS_TEST  = -g -o


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

LIB_OPT   = 
LIB_DEBUG = 
LIB_TEST  =

LIB_HDF4  = 
LIB_HDF5  = -L ${HDF5_PATH}/lib -lhdf5 

LIB_PAPI  =
LIB_MATH  = 

LIB_MPI   = 
LIB_NCMPI = -L ${NCMPI_PATH}/lib -lpnetcdf
LIB_MPE   =

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

#Configure lines for complete software stack (install in this order):
#
#valgrind-3.4.0:
#./configure --prefix=/opt/valgrind-3.4.0 --without-mpicc
#
#openmpi-1.3:
#./configure --prefix=/opt/openmpi-1.3 FC=gfortran F90=gfortran \
#--enable-debug --enable-memchecker --with-valgrind=/opt/valgrind-3.4.0 \
#--with-mpi-f90-size=medium --with-f90-max-array-dim=5 --enable-static
#
#hdf5-1.6.8:
#./configure --prefix=/opt/hdf5-1.6.8 CC=mpicc --enable-parallel \
#--enable-production --disable-stream-vfd --enable-threadsafe \
#--enable-debug=all --enable-parallel --enable-shared
#
#hdf5-1.8.3_srcbuild:
#./configure --prefix=/opt/hdf5-1.8.3_srcbuild CC=mpicc --disable-production \
#--enable-debug=all --enable-trace --enable-parallel --enable-shared \
#--enable-static
#
#pdtoolkit-3.14.1:
#./configure -prefix=/opt/pdtoolkit-3.14.1 -GNU
#
#tau-2.18.1:
#./configure -prefix=/opt/tau-2.18.1 -c++=g++ -cc=gcc -fortran=gnu \
#-pdt=/opt/pdtoolkit-3.14.1 -mpi -mpiinc=/opt/openmpi-1.3/include \
#-mpilib=/opt/openmpi-1.3/lib -PROFILECALLPATH
#
#parallel-netcdf-1.0.3:
#./configure --prefix=/opt/parallel-netcdf-1.0.3 --enable-fortran \
#MPICC=mpicc CFLAGS=-g MPIF77=mpif77
#
