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
# ---------------------------------------------------------------------------- #

# BLASLAPACKFLAGS = -DNO_LAPACK
# BLASLAPACKLIBS = -lblas

BLASLAPACKFLAGS =
BLASLAPACKLIBS = -lopenblas

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

# PNGPPFLAGS = -DWITH_PNG -I./libs
# PNGPPLIBS = -lpng
                                                                               #
################################################################################


# compiler flags
CFLAGS += -fPIC -std=c++11 -pedantic-errors -DGIT_COMMIT=\"$(GIT_COMMIT)\" \
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

# source files
CXXFILES = amplitudes arrays gauss bspline io hdffile hydrogen luft main \
           matrix misc radial slater solver special version \
           preconditioners/NoPreconditioner preconditioners/CGPreconditioner \
           preconditioners/KPAPreconditioner preconditioners/ILUPreconditioner \
           preconditioners/GPUPreconditioner preconditioners/HybPreconditioner

# sources and derived files (object files, dep-files, etc)
BASECXX      = $(shell basename $(CXX))
CXXSOURCES   = $(patsubst %,src/%.cpp,$(CXXFILES))
OBJECTS      = $(patsubst %,obj/$(BASECXX)/%.o,  $(CXXFILES))
DEPENDENCIES = $(patsubst %,obj/$(BASECXX)/%.dep,$(CXXFILES))
OCLSOURCES   = src/preconditioners/GPUPreconditioner.cl
OCLKERNELS   = src/preconditioners/GPUPreconditioner.inc

# build all targets
all: pollversion bin/hex-ecs.$(BASECXX)

# dummy target to force update of git revision
pollversion:
	@touch src/version.cpp

# build documentation from source comments using Doxygen
doc:
	doxygen

# remove documentation
docclean:
	rm -rf html

# make distribution archive
dist:
	mkdir -p hex-ecs/libs
	cp --recursive --dereference README Makefile Doxyfile *.cbp ../media/hexe-small.gif src hex-ecs/
	cp --recursive libs/png++ hex-ecs/libs
	sed -i "s/\(^GIT_COMMIT = \)\(.*\)/\1$(GIT_COMMIT)/g" hex-ecs/Makefile
	sed -i "s/\(^all:.*\)\(pollversion\)\(.*\)/\1\3/g" hex-ecs/Makefile
	sed -i "s/\.\.\/media\///g" hex-ecs/Doxyfile
	tar --create --gzip --file=hex-ecs-$(GIT_COMMIT).tar.gz hex-ecs
	rm -rf hex-ecs/

# remove build files
clean:
	rm -f bin/hex-ecs.$(BASECXX) $(OBJECTS) $(DEPENDENCIES)

# remove everything (build files, documentation and distribution files)
allclean: clean docclean distclean

# remove distribution files
distclean:
	rm -rf hex-ecs hex-ecs-*.tar.gz

# main program target
bin/hex-ecs.$(BASECXX): $(OCLKERNELS) $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(OBJECTS) $(LFLAGS) -o $@
	@ln -fs hex-ecs.$(BASECXX) bin/hex-ecs

# CL -> INC rule (convert OpenCL sources to C byte array)
src/%.inc: src/%.cl
	xxd -i $< | tail -n+2 | head -n-2 > $@

# CPP -> O rule (compile source files)
obj/$(BASECXX)/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

# include generated source dependencied
-include $(DEPENDENCIES)

