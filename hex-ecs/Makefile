.PHONY: all doc docclean clean allclean pollversion dist distclean

################################################################################
# Editable section                                                             #

# ---------------------------------------------------------------------------- #
# Compiler and linker + flags                                                  #
#  - Hex-ecs requires compiler with C++11 support (GCC >= 4.8, ICC >= 2015)    #
# ---------------------------------------------------------------------------- #

# C++ compiler (g++, clang, ...)
CXX ?= g++

# build mode ("debug" or "release")
TARGET ?= release

# global compiler flags
#
#    Single precision: -DSINGLE
#    Double precision: -fcx-limited-range (optional)
#
CFLAGS ?= -fopenmp -Wall -mtune=native -march=native -fcx-limited-range -Wno-unknown-pragmas
LFLAGS ?= -fopenmp -lstdc++

# ---------------------------------------------------------------------------- #
# OpenMPI environment (optional)                                               #
# ---------------------------------------------------------------------------- #

OPENMPIFLAGS ?= -DWITH_MPI $(shell mpicxx --showme:compile)
OPENMPILIBS  ?= $(shell mpicxx --showme:link)

# ---------------------------------------------------------------------------- #
# BLAS & LAPACK                                                                #
#                                                                              #
# The macro DISABLE_PARALLEL_PRECONDITION will inhibit nested parallelism      #
# e.g. when using threaded OpenBLAS. This may (and may not) deliver            #
# performance benefits in exchange of giving up possibility to precondition    #
# several blocks concurrently.                                                 #
# ---------------------------------------------------------------------------- #

# BLASLAPACKFLAGS = 
# BLASLAPACKLIBS = -lblas -llapack

BLASLAPACKFLAGS ?= -DDISABLE_PARALLEL_PRECONDITION
BLASLAPACKLIBS ?= -lopenblas

# ---------------------------------------------------------------------------- #
# LU decomposition libraries                                                   #
# - strongly recommended, otherwise everything will be done using LAPACK       #
# ---------------------------------------------------------------------------- #

# SuiteSparse/UMFPACK
# - ILP64 mode needed already for medium-size matrices (-D_LONGINT)
SUITESPARSEFLAGS ?= -DWITH_UMFPACK -D_LONGINT
SUITESPARSELIBS ?= -lumfpack -lamd -lsuitesparseconfig

# Sequential SuperLU
# SUPERLUFLAGS = -DWITH_SUPERLU
# SUPERLULIBS = -lsuperlu

# Distributed SuperLU + ParMetis (LP64 mode)
# SUPERLUDISTFLAGS = -DWITH_SUPERLU_DIST -U_LONGINT
# SUPERLUDISTLIBS = -lsuperlu_dist_lp64 -lparmetis_lp64

# Distributed SuperLU + ParMetis (ILP64 mode)
# SUPERLUDISTFLAGS = -DWITH_SUPERLU_DIST -D_LONGINT
# SUPERLUDISTLIBS = -lsuperlu_dist_ilp64 -lparmetis_ilp64

# Distributed MUMPS
# MUMPSFLAGS = -DWITH_MUMPS
# MUMPSLIBS = -lmumps

# ScaLAPACK
# SCALAPACKFLAGS = -DWITH_SCALAPACK
# SCALAPACKLIBS = -lscalapack

# Intel MKL Pardiso
# - incompatible with Pardiso 5.0.0
# PARDISOFLAGS = -DWITH_MKL
# PARDISOLIBS = -mkl

# Pardiso 5.0.0
# - incompatible with Intel MKL Pardiso
# PARDISOFLAGS = -DWITH_PARDISO
# PARDISOLIBS = -lpardiso500-MPI-GNU472-X86-64

# ---------------------------------------------------------------------------- #
# GSL                                                                          #
# ---------------------------------------------------------------------------- #

GSLFLAGS ?= $(shell gsl-config --cflags)
GSLLIBS ?= $(shell gsl-config --libs)

# ---------------------------------------------------------------------------- #
# OpenCL (optional)                                                            #
# - Use the flag -DOPENCL2 if you are sure that your implementation supports   #
#   OpenCL 2.0+. This will avoid some deprecated functions.                    #
# ---------------------------------------------------------------------------- #

# OPENCLFLAGS = -DWITH_OPENCL
# OPENCLLIBS = -lOpenCL

# ---------------------------------------------------------------------------- #
# PNG++ (optional, for debugging purposes)                                     #
# ---------------------------------------------------------------------------- #

# PNGPPFLAGS = -DWITH_PNG -I../libs
# PNGPPLIBS = -lpng
                                                                               #
################################################################################


# compiler flags
CFLAGS += -fPIC -std=c++11 -pedantic-errors -DGIT_COMMIT=\"$(GIT_COMMIT)\" \
          -Isrc -I../common \
          $(SUITESPARSEFLAGS) $(SUPERLUFLAGS) $(SUPERLUDISTFLAGS) $(MUMPSFLAGS) $(PARDISOFLAGS)
          $(OPENMPIFLAGS) $(OPENCLFLAGS) $(GSLFLAGS) $(BLASLAPACKFLAGS) $(PNGPPFLAGS)

# linker flags
LFLAGS += -fPIC -rdynamic \
          $(SUITESPARSELIBS) $(SUPERLULIBS) $(SUPERLUDISTLIBS) $(MUMPSLIBS) $(PARDISOLIBS) \
          $(OPENMPILIBS) $(OPENCLLIBS) $(GSLLIBS) $(BLASLAPACKLIBS) $(PNGPPLIBS)

# git revision number
GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

# DEBUG and RELEASE variants
ifeq ($(TARGET),debug)
CXXFLAGS = -g $(CFLAGS) -O0
else ifeq ($(TARGET),release)
CXXFLAGS = -g $(CFLAGS) -O3 -DNDEBUG
else
$(error "Unknown target $(TARGET)")
endif

# compiler name
BASECXX = $(shell basename $(CXX))

# common source files
COMMON_FILES   = hex-arrays hex-blas hex-hdffile hex-hydrogen \
                 hex-matrix hex-misc hex-special hex-version hex-vtkfile
COMMON_SOURCES = $(patsubst %,../common/%.cpp,$(COMMON_FILES))
COMMON_OBJECTS = $(patsubst %,obj/$(BASECXX)/common/%.o,$(COMMON_FILES))
COMMON_DEPS    = $(patsubst %,obj/$(BASECXX)/common/%.dep,$(COMMON_FILES))

# source files specific to hex-ecs
HEX_ECS_FILES   = amplitudes ang gauss bspline hldata io luft main preconditioners radial slater solver \
                  factorizers/lu-umfpack factorizers/lu-mumps factorizers/lu-superlu factorizers/lu-lapack \
                  factorizers/lu-superlu-dist factorizers/lu-pardiso factorizers/lu-scalapack \
                  preconditioners/NoPreconditioner preconditioners/CGPreconditioner \
                  preconditioners/KPAPreconditioner preconditioners/ILUPreconditioner \
                  preconditioners/GPUPreconditioner preconditioners/HybPreconditioner \
                  preconditioners/DOMPreconditioner preconditioners/CoupledPreconditioner
HEX_ECS_SOURCES = $(patsubst %,src/%.cpp,$(HEX_ECS_FILES))
HEX_ECS_OBJECTS = $(patsubst %,obj/$(BASECXX)/hex-ecs/%.o,$(HEX_ECS_FILES))
HEX_ECS_DEPS    = $(patsubst %,obj/$(BASECXX)/hex-ecs/%.dep,$(HEX_ECS_FILES))

# OpenCL kernels
OCLSOURCES = src/preconditioners/GPUPreconditioner.cl
OCLKERNELS = src/preconditioners/GPUPreconditioner.inc

# build all targets
all: pollversion bin/hex-ecs.$(BASECXX)

# dummy target to force update of git revision
pollversion:
	@touch ../common/hex-version.cpp

# build documentation from source comments using Doxygen
doc:
	doxygen

# remove documentation
docclean:
	rm -rf html

# remove build files
clean:
	rm -f bin/hex-ecs bin/hex-ecs.$(BASECXX) $(COMMON_OBJECTS) $(COMMON_DEPS) $(HEX_ECS_OBJECTS) $(HEX_ECS_DEPS)

# remove everything (build files, documentation and distribution files)
allclean: clean docclean

# main program target
bin/hex-ecs.$(BASECXX): $(OCLKERNELS) $(COMMON_OBJECTS) $(HEX_ECS_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(HEX_ECS_OBJECTS) $(COMMON_OBJECTS) $(LFLAGS) -o $@
	@ln -fs hex-ecs.$(BASECXX) bin/hex-ecs

# CL -> INC rule (convert OpenCL sources to C byte array)
%.inc: %.cl
	xxd -i $< > $@

# CPP -> O rule (compile source files)
obj/$(BASECXX)/common/%.o: ../common/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)
obj/$(BASECXX)/hex-ecs/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

# include generated source dependencied
-include $(COMMON_DEPS) $(HEX_ECS_DEPS)
