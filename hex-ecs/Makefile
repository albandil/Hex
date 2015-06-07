.PHONY: all doc docclean clean allclean pollversion dist distclean

################################################################################
# Editable section                                                             #

# ---------------------------------------------------------------------------- #
# Compiler and linker + flags                                                  #
# ---------------------------------------------------------------------------- #

CXX = g++

CFLAGS = -fopenmp -Wall -mtune=native -march=native
LFLAGS = -fopenmp

# ---------------------------------------------------------------------------- #
# LU factorization libraries (use at least one)                                #
# ---------------------------------------------------------------------------- #

# - sequential SuperLU
SUPERLUFLAGS = -DWITH_SUPERLU
SUPERLULIBS = -lsuperlu

# - distributed SuperLU (incompatible with SuperLU)
# SUPERLUDISTFLAGS = -DWITH_SUPERLU_DIST
# SUPERLUDISTLIBS = -lsuperlu_dist -lptscotchmetis -lmetis

# - SuiteSparse/UMFPACK
SUITESPARSEFLAGS = -DWITH_UMFPACK
SUITESPARSELIBS = -lumfpack -lamd -lsuitesparseconfig -lopenblas

# ---------------------------------------------------------------------------- #
# OpenMPI environment (optional)                                               #
# ---------------------------------------------------------------------------- #

OPENMPIFLAGS = -DWITH_MPI $(shell mpic++ --showme:compile)
OPENMPILIBS  = $(shell mpic++ --showme:link)

# ---------------------------------------------------------------------------- #
# OpenCL (optional)                                                            #
# ---------------------------------------------------------------------------- #

# - Intel OpenCL
OPENCLFLAGS = -DWITH_OPENCL -I/opt/intel/opencl-1.2-4.5.0.8/include
OPENCLLIBS = -L/opt/intel/opencl-1.2-4.5.0.8/lib64 -lOpenCL

# - AMD OpenCL
# OPENCLFLAGS = -DWITH_OPENCL -I/opt/AMDAPPSDK-3.0-0-Beta/include
# OPENCLLIBS = -L/opt/AMDAPPSDK-3.0-0-Beta/lib/x86_64/sdk -lOpenCL

# ---------------------------------------------------------------------------- #
# GSL                                                                          #
# ---------------------------------------------------------------------------- #

GSLFLAGS = $(shell gsl-config --cflags)
GSLLIBS = $(shell gsl-config --libs)

# ---------------------------------------------------------------------------- #
# FFTW                                                                         #
# ---------------------------------------------------------------------------- #

FFTWFLAGS = 
FFTWLIBS = -lfftw3

# ---------------------------------------------------------------------------- #
# Lapack (optional if not already in OpenBLAS; use -DNO_LAPACK to omit)        #
# ---------------------------------------------------------------------------- #

# LAPACKFLAGS = -DNO_LAPACK
# LAPACKLIBS = -lgfortran -llapack

# ---------------------------------------------------------------------------- #
# PNG++ (optional)                                                             #
# ---------------------------------------------------------------------------- #

PNGPPFLAGS = -DWITH_PNG -I./libs
PNGPPLIBS = -lpng
                                                                               #
################################################################################


# compiler flags
CFLAGS +=  -fPIC -std=c++11 -DGIT_COMMIT=\"$(GIT_COMMIT)\" \
          $(SUPERLUFLAGS) $(SUPERLUDISTFLAGS)  $(SUITESPARSEFLAGS) $(OPENMPIFLAGS) \
          $(OPENCLFLAGS) $(GSLFLAGS) $(FFTWFLAGS) $(LAPACKFLAGS) $(PNGPPFLAGS) 

# linker flags
LFLAGS += -fPIC -rdynamic \
          $(SUPERLULIBS) $(SUPERLUDISTLIBS) $(SUITESPARSELIBS) $(OPENMPILIBS) \
          $(OPENCLLIBS) $(GSLLIBS) $(FFTWLIBS) $(LAPACKLIBS) $(PNGPPLIBS)

# git revision number
GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

# DEBUG and RELEASE variants
CXXDEBUG = -g $(CFLAGS) 
CXXRILIS = -g $(CFLAGS) -O3 -fcx-limited-range -DNDEBUG

# choose between DEBUG and RELEASE compiler flags
# CXXFLAGS = $(CXXDEBUG)
CXXFLAGS = $(CXXRILIS)

# source files
CXXFILES = amplitudes arrays gauss bspline io hdffile hydrogen luft main \
           matrix misc radial slater special version \
           preconditioners/NoPreconditioner preconditioners/CGPreconditioner \
           preconditioners/KPAPreconditioner preconditioners/ILUPreconditioner \
           preconditioners/GPUPreconditioner

# sources and derived files (object files, dep-files, etc)
CXXSOURCES   = $(patsubst %,src/%.cpp,$(CXXFILES))
OBJECTS      = $(patsubst %,obj/%.o,  $(CXXFILES))
DEPENDENCIES = $(patsubst %,obj/%.dep,$(CXXFILES))
OCLSOURCES   = src/preconditioners/GPUPreconditioner.cl
OCLKERNELS   = src/preconditioners/GPUPreconditioner.inc

# build all targets
all: pollversion bin/hex-ecs

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
	rm -f bin/hex-ecs $(OBJECTS) $(DEPENDENCIES)

# remove everything (build files, documentation and distribution files)
allclean: clean docclean distclean

# remove distribution files
distclean:
	rm -rf hex-ecs hex-ecs-*.tar.gz

# main program target
bin/hex-ecs: $(OCLKERNELS) $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(OBJECTS) $(LFLAGS) -o $@

# CL -> INC rule (convert OpenCL sources to C byte array)
src/%.inc: src/%.cl
	xxd -i $< | tail -n+2 | head -n-2 > $@

# CPP -> O rule (compile source files)
obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

# include generated source dependencied
-include $(DEPENDENCIES)
