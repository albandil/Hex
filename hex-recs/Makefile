.PHONY: all doc docclean clean allclean pollversion dist distclean


################################################################################
# Editable section                                                             #

# ---------------------------------------------------------------------------- #
# Compiler and linker + flags                                                  #
# ---------------------------------------------------------------------------- #

CXX = g++

CFLAGS = -Wall -std=c++11 -fopenmp -fPIC -mtune=native -march=native
LFLAGS = -fopenmp -fPIC

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

OPENMPIFLAGS = -DNO_MPI
OPENMPILIBS =

# - for repository release
# OPENMPIFLAGS = -I/usr/lib64/mpi/gcc/openmpi/include
# OPENMPILIBS = -L/usr/lib64/mpi/gcc/openmpi/lib64 -lmpi -lmpi_cxx

# ---------------------------------------------------------------------------- #
# OpenCL environment (optional, use -DNO_OPENCL to omit)                       #
# ---------------------------------------------------------------------------- #

OPENCLFLAGS = -DNO_OPENCL

# - for AMD OpenCL
# OPENCLFLAGS = -I/opt/AMDAPP/include
# OPENCLLIBS = -lOpenCL

# - for Intel OpenCL
# OPENCLFLAGS = -I/opt/intel/opencl-1.2-3.1.1.11385/include
# OPENCLLIBS = -L/opt/intel/opencl-1.2-3.1.1.11385/lib64 -lOpenCL

# - for NVidia OpenCL
# OPENCLFLAGS = -I/usr/include/nvidia-current/
# OPENCLLIBS = -lOpenCL

# ---------------------------------------------------------------------------- #
# HDF                                                                          #
# ---------------------------------------------------------------------------- #

HDFFLAGS = 
HDFLIBS = -lhdf5_cpp -lhdf5

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
# Lapack (optional, use -DNO_LAPACK to omit)                                   #
# ---------------------------------------------------------------------------- #

LAPACKFLAGS = -DNO_LAPACK
# LAPACKLIBS = -lgfortran -llapack

# ---------------------------------------------------------------------------- #
# PNG++ (optional, use -DNO_PNG to omit)                                       #
# ---------------------------------------------------------------------------- #

PNGPPFLAGS = -DNO_PNG
# PNGPPFLAGS = -I./libs
# PNGPPLIBS = -lpng
                                                                               #
################################################################################


CFLAGS += $(SUITESPARSEFLAGS) $(OPENMPIFLAGS) $(OPENCLFLAGS) $(HDFFLAGS) \
          $(GSLFLAGS) $(FFTWFLAGS) $(LAPACKFLAGS) $(PNGPPFLAGS)

LFLAGS += $(SUITESPARSELIBS) $(OPENMPILIBS) $(OPENCLLIBS) $(HDFLIBS) \
          $(GSLLIBS) $(FFTWLIBS) $(LAPACKLIBS) $(PNGPPLIBS)

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CXXDEBUG = -g $(CFLAGS) -DGIT_COMMIT=\"$(GIT_COMMIT)\"
CXXRILIS = -O3 -fcx-limited-range $(CFLAGS) -DNDEBUG -DGIT_COMMIT=\"$(GIT_COMMIT)\"

# CXXFLAGS = $(CXXDEBUG)
CXXFLAGS = -g $(CXXRILIS)

FILES = amplitudes arrays gauss bspline hdffile hydrogen main input \
	matrix preconditioners radial slater special version

SOURCES = $(patsubst %,src/%.cpp,$(FILES))
OBJECTS = $(patsubst %,obj/%.o,$(FILES))
DEPENDENCIES = $(patsubst %,obj/%.dep,$(FILES))

all: pollversion bin/hex-recs

pollversion:
	@touch src/version.cpp

doc:
	doxygen

docclean:
	rm -rf html

dist:
	mkdir -p hex-recs/libs
	cp --recursive --dereference README Makefile Doxyfile *.cbp ../media/hexe-small.gif src hex-recs/
	cp --recursive libs/png++ hex-recs/libs
	sed -i "s/\(^GIT_COMMIT = \)\(.*\)/\1$(GIT_COMMIT)/g" hex-recs/Makefile
	sed -i "s/\(^all:.*\)\(pollversion\)\(.*\)/\1\3/g" hex-recs/Makefile
	sed -i "s/\.\.\/media\///g" hex-recs/Doxyfile
	tar --create --gzip --file=hex-recs-$(GIT_COMMIT).tar.gz hex-recs
	rm -rf hex-recs/

clean:
	rm -f bin/hex-recs $(OBJECTS) $(DEPENDENCIES)

allclean: clean docclean distclean

distclean:
	rm -rf hex-recs hex-recs-*.tar.gz

bin/hex-recs: $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $^ $(LFLAGS) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

src/kernels_cl.c: src/kernels.cl
	xxd -i src/kernels.cl | tail -n+2 | head -n-2 > src/kernels_cl.c

-include $(wildcard obj/*.dep)
