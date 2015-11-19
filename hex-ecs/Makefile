.PHONY: all doc docclean clean allclean pollversion dist distclean

################################################################################
# Editable section                                                             #

# ---------------------------------------------------------------------------- #
# Compiler and linker + flags                                                  #
#  - Hex-ecs requires compiler with C++11 support (GCC >= 4.8, ICC >= 2015)    #
# ---------------------------------------------------------------------------- #

# C++ compiler (g++, clang, ...)
CXX = g++

# build mode ("debug" or "release")
TARGET = release

# global compiler flags
CFLAGS = -fopenmp -Wall -mtune=native -march=native -fcx-limited-range
LFLAGS = -fopenmp -lstdc++

# ---------------------------------------------------------------------------- #
# BLAS & LAPACK                                                                #
#                                                                              #
# The macro DISABLE_PARALLEL_PRECONDITION will inhibit nested parallelism      #
# e.g. when using threaded OpenBLAS. This may (and may not) deliver            #
# performance benefits in exchange of giving up possibility to precondition    #
# several blocks concurrently.                                                 #
# ---------------------------------------------------------------------------- #

BLASLAPACKFLAGS = 
BLASLAPACKLIBS = -lblas -llapack

# BLASLAPACKFLAGS = -DDISABLE_PARALLEL_PRECONDITION
# BLASLAPACKLIBS = -lopenblas

# ---------------------------------------------------------------------------- #
# SuiteSparse/UMFPACK                                                          #
#  - needed for sparse matrix manipulation                                     #
#  - default choice for LU factorization (e.g. in ILU preconditioner)          #
#  - supports both LP64 (32-bit integer) and ILP64 (64-bit integer)            #
#  - already the factorization of moderate matrices requires ILP64, so         #
#    the flag -D_LONGINT is included by default                                #
# ---------------------------------------------------------------------------- #

SUITESPARSEFLAGS = -DWITH_UMFPACK -D_LONGINT
SUITESPARSELIBS = -lumfpack -lamd -lsuitesparseconfig

# ---------------------------------------------------------------------------- #
# Other LU factorization libraries                                             #
#  - these are optional                                                        #
#  - use only one at a time                                                    #
# ---------------------------------------------------------------------------- #

# - sequential SuperLU with only 32-bit integer support
# SUPERLUFLAGS = -DWITH_SUPERLU -U_LONGINT
# SUPERLULIBS = -lsuperlu

# - distributed SuperLU + Metis compiled with 32-bit integer (LP64 mode)
# SUPERLUDISTFLAGS = -DWITH_SUPERLU_DIST -U_LONGINT
# SUPERLUDISTLIBS = -lsuperlu_dist_lp64 -lmetis_lp64

# - distributed SuperLU + Metis compiled with 64-bit integer (ILP64 mode)
# SUPERLUDISTFLAGS = -DWITH_SUPERLU_DIST -D_LONGINT
# SUPERLUDISTLIBS = -lsuperlu_dist_ilp64 -lmetis_ilp64

# ---------------------------------------------------------------------------- #
# MUMPS library for the out-of-core coupled solver.                            #
# ---------------------------------------------------------------------------- #

# MUMPSFLAGS = -DWITH_MUMPS
# MUMPSLIBS = -lzmumps_seq

# ---------------------------------------------------------------------------- #
# GSL                                                                          #
# ---------------------------------------------------------------------------- #

GSLFLAGS = $(shell gsl-config --cflags)
GSLLIBS = $(shell gsl-config --libs)

# ---------------------------------------------------------------------------- #
# OpenMPI environment (optional)                                               #
# ---------------------------------------------------------------------------- #

OPENMPIFLAGS = -DWITH_MPI $(shell mpic++ --showme:compile)
OPENMPILIBS  = $(shell mpic++ --showme:link)

# ---------------------------------------------------------------------------- #
# OpenCL (optional)                                                            #
# ---------------------------------------------------------------------------- #

# OPENCLFLAGS = -DWITH_OPENCL
# OPENCLLIBS = -lOpenCL

# ---------------------------------------------------------------------------- #
# PNG++ (optional)                                                             #
# ---------------------------------------------------------------------------- #

# PNGPPFLAGS = -DWITH_PNG -I../libs
# PNGPPLIBS = -lpng
                                                                               #
################################################################################


# compiler flags
CFLAGS += -fPIC -std=c++11 -pedantic-errors -DGIT_COMMIT=\"$(GIT_COMMIT)\" \
          -Isrc -I../common \
          $(SUITESPARSEFLAGS) $(SUPERLUFLAGS) $(SUPERLUDISTFLAGS) $(OPENMPIFLAGS) \
          $(OPENCLFLAGS) $(GSLFLAGS) $(BLASLAPACKFLAGS) $(PNGPPFLAGS) 

# linker flags
LFLAGS += -fPIC -rdynamic \
          $(SUITESPARSELIBS) $(SUPERLULIBS) $(SUPERLUDISTLIBS) $(OPENMPILIBS) \
          $(OPENCLLIBS) $(GSLLIBS) $(BLASLAPACKLIBS) $(PNGPPLIBS)

# git revision number
GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

# DEBUG and RELEASE variants
ifeq ($(TARGET),debug)
CXXFLAGS = -g $(CFLAGS)
else ifeq ($(TARGET),release)
CXXFLAGS = -g $(CFLAGS) -O3 -DNDEBUG
else
$(error "Unknown target $(TARGET)")
endif

# compiler name
BASECXX = $(shell basename $(CXX))

# common source files
COMMON_FILES   = hex-arrays hex-hdffile hex-hydrogen hex-luft hex-matrix \
                 hex-misc hex-special hex-version
COMMON_SOURCES = $(patsubst %,../common/%.cpp,$(COMMON_FILES))
COMMON_OBJECTS = $(patsubst %,obj/$(BASECXX)/common/%.o,$(COMMON_FILES))
COMMON_DEPS    = $(patsubst %,obj/$(BASECXX)/common/%.dep,$(COMMON_FILES))

# source files specific to hex-ecs
HEX_ECS_FILES   = amplitudes ang gauss bspline io main radial slater solver \
                  preconditioners/NoPreconditioner preconditioners/CGPreconditioner \
                  preconditioners/KPAPreconditioner preconditioners/ILUPreconditioner \
                  preconditioners/GPUPreconditioner preconditioners/HybPreconditioner
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
