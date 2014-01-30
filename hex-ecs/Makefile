.PHONY: all doc docclean clean allclean pollversion dist distclean

###################################################################################
# Editable section

# Compiler executable
CXX = g++
# CXX = icpc

# Suite Sparse environment
SUITESPARSENAMES = umfpack amd cholmod colamd camd ccolamd metis suitesparseconfig
SUITESPARSELIBS = $(patsubst %,-l%,$(SUITESPARSENAMES))

# OpenMPI environment
OPENMPIFLAGS = -DNO_MPI
#OPENMPIINCLUDE = -I/usr/lib64/mpi/gcc/openmpi/include
#OPENMPILIBS = -L/usr/lib64/mpi/gcc/openmpi/lib64 -lmpi -lmpi_cxx

# OpenCL environment
OPENCLFLAGS = -DNO_OPENCL
# AMD
#OPENCLINCLUDE = -I/opt/AMDAPP/include
#OPENCLLIBS = -lOpenCL
# Intel
#OPENCLINCLUDE = -I/opt/intel/opencl-1.2-3.1.1.11385/include
#OPENCLLIBS = -L/opt/intel/opencl-1.2-3.1.1.11385/lib64 -lOpenCL
# NVidia
#OPENCLINCLUDE = -I/usr/include/nvidia-current/
#OPENCLLIBS = -lOpenCL

# Lapack
LAPACKFLAGS = -DNO_LAPACK
#LAPACKINCLUDE =
#LAPACKLIBS = -lgfortran -llapack

# PNG++
PNGFLAGS = -DNO_PNG
#PNGPPINCLUDE = -I./libs
#PNGPPLIBS = -lpng
                                                                                    #
#####################################################################################

INCLUDE = $(OPENMPIINCLUDE) $(OPENCLINCLUDE) $(LAPACKINCLUDE) $(PNGPPINCLUDE)
LIBRARY = $(SUITESPARSELIBS) $(OPENMPILIBS) \
	-lrt -lhdf5_cpp -lhdf5 -lgsl -lopenblas -lgslcblas -lpthread \
	-lm -lfftw3 $(OPENMPILIBS) $(OPENCLLIBS) $(LAPACKLIBS) $(PNGPPLIBS)

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CFLAGS_DEBUG = -g -fopenmp -Wall -std=c++11 -fPIC $(INCLUDE) -DGIT_COMMIT=\"$(GIT_COMMIT)\"
CFLAGS_RELEASE = -O3 -fopenmp -Wall -std=c++11 -fPIC $(INCLUDE) -DNDEBUG -DGIT_COMMIT=\"$(GIT_COMMIT)\"

LFLAGS = -fopenmp -fPIC $(LIBRARY)

# CFLAGS = $(CFLAGS_DEBUG) $(LAPACKFLAGS) $(PNGFLAGS) $(OPENCLFLAGS) $(OPENMPIFLAGS)
CFLAGS = $(CFLAGS_RELEASE) $(LAPACKFLAGS) $(PNGFLAGS) $(OPENCLFLAGS) $(OPENMPIFLAGS) \
         -mtune=native -march=native

FILES = amplitudes arrays gauss bspline hdffile main input \
	matrix preconditioners radial slater specf version

SOURCES = $(patsubst %,src/%.cpp,$(FILES))
OBJECTS = $(patsubst %,obj/%.o,$(FILES))
DEPENDENCIES = $(patsubst %,obj/%.dep,$(FILES))

all: pollversion bin/hex-ecs

pollversion:
	@touch src/version.cpp

doc:
	doxygen

docclean:
	rm -rf html

dist:
	mkdir -p hex-ecs/libs
	cp --recursive --dereference README Makefile Doxyfile *.cbp ../media/hexe-small.gif src hex-ecs/
	cp --recursive libs/png++ hex-ecs/libs
	sed -i "s/\(^GIT_COMMIT = \)\(.*\)/\1$(GIT_COMMIT)/g" hex-ecs/Makefile
	sed -i "s/\(^all:.*\)\(pollversion\)\(.*\)/\1\3/g" hex-ecs/Makefile
	sed -i "s/\.\.\/media\///g" hex-ecs/Doxyfile
	tar --create --gzip --file=hex-ecs-$(GIT_COMMIT).tar.gz hex-ecs
	rm -rf hex-ecs/

clean:
	rm -f bin/hex $(OBJECTS) $(DEPENDENCIES)

allclean: clean docclean distclean

distclean:
	rm -rf hex-ecs hex-ecs-*.tar.gz

bin/hex-ecs: $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $^ $(LFLAGS) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

src/kernels_cl.c: src/kernels.cl
	xxd -i src/kernels.cl | tail -n+2 | head -n-2 > src/kernels_cl.c

-include $(wildcard obj/*.dep)
