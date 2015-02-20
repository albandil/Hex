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
# Suite Sparse environment                                                     #
# ---------------------------------------------------------------------------- #

SUITESPARSENAMES = umfpack amd cholmod colamd camd ccolamd

# - for repository release
# SUITESPARSEFLAGS = -I/usr/include/suitesparse
# SUITESPARSELIBS = $(patsubst %,-l%,$(SUITESPARSENAMES))

# - for custom build
SUITESPARSEFLAGS =
SUITESPARSELIBS = $(patsubst %,-l%,$(SUITESPARSENAMES)) \
                  -lmetis -lsuitesparseconfig -lopenblas

# ---------------------------------------------------------------------------- #
# OpenMPI environment (optional, use -DNO_MPI to omit)                         #
# ---------------------------------------------------------------------------- #

# OPENMPIFLAGS = -DNO_MPI
# OPENMPILIBS =

# - for repository release
OPENMPIFLAGS = $(shell mpic++ --showme:compile)
OPENMPILIBS  = $(shell mpic++ --showme:link)

# ---------------------------------------------------------------------------- #
# HDF                                                                          #
# ---------------------------------------------------------------------------- #

HDFFLAGS = 
HDFLIBS = -lhdf5

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
# PNG++ (optional, use -DNO_PNG to omit)                                       #
# ---------------------------------------------------------------------------- #

PNGPPFLAGS = -DNO_PNG
PNGPPFLIBS =

# PNGPPFLAGS = -I./libs
# PNGPPLIBS = -lpng
                                                                               #
################################################################################


# compiler flags
CFLAGS +=  -fPIC -std=c++11 $(SUITESPARSEFLAGS) $(OPENMPIFLAGS) $(OPENCLFLAGS) $(HDFFLAGS) \
          $(GSLFLAGS) $(FFTWFLAGS) $(LAPACKFLAGS) $(PNGPPFLAGS)

# linker flags
LFLAGS += -fPIC $(SUITESPARSELIBS) $(OPENMPILIBS) $(OPENCLLIBS) $(HDFLIBS) \
          $(GSLLIBS) $(FFTWLIBS) $(LAPACKLIBS) $(PNGPPLIBS) -rdynamic

# git revision number
GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

# DEBUG and RELEASE variants
CXXDEBUG = -g $(CFLAGS) -DGIT_COMMIT=\"$(GIT_COMMIT)\"
CXXRILIS = -O3 -fcx-limited-range $(CFLAGS) -DNDEBUG -DGIT_COMMIT=\"$(GIT_COMMIT)\"

# choose between DEBUG and RELEASE compiler flags
# CXXFLAGS = $(CXXDEBUG)
CXXFLAGS = $(CXXRILIS)

# source files
CXXFILES = amplitudes arrays gauss bspline hdffile hydrogen main io \
           matrix radial slater special version \
           preconditioners/NoPreconditioner preconditioners/CGPreconditioner \
           preconditioners/KPAPreconditioner preconditioners/ILUPreconditioner

# sources and derived files (object files, dep-files, etc)
CXXSOURCES   = $(patsubst %,src/%.cpp,$(CXXFILES))
OBJECTS      = $(patsubst %,obj/%.o,  $(CXXFILES))
DEPENDENCIES = $(patsubst %,obj/%.dep,$(CXXFILES))

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
bin/hex-ecs: $(OBJECTS)
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
