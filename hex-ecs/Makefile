.PHONY: all doc clean pollversion

CPP = g++

SUITESPARSENAMES = umfpack amd cholmod colamd camd ccolamd metis suitesparseconfig
SUITESPARSELIBS = $(patsubst %,-l%,$(SUITESPARSENAMES))

OPENMPIINCLUDE = -I/usr/lib64/mpi/gcc/openmpi/include
OPENMPILIBS = -L/usr/lib64/mpi/gcc/openmpi/lib64

OPENCLINCLUDE = -I/opt/AMDAPP/include
OPENCLLIBS = -lOpenCL

INCLUDE = $(OPENMPIINCLUDE) -I$(shell pwd)/../libs $(OPENCLINCLUDE)
LIBRARY = $(SUITESPARSELIBS) $(OPENMPILIBS) $(OPENCLLIBS) \
	-lmpi -lmpi_cxx -lrt -lhdf5_cpp -lhdf5 -lgsl -lopenblas -lgslcblas -lpthread \
	-lm -lfftw3 -lgfortran -lpng -llapack

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CFLAGS_DEBUG = -g -fopenmp -Wall -std=c++11 -fPIC $(INCLUDE) -DGIT_COMMIT=\"$(GIT_COMMIT)\"
LFLAGS_DEBUG = -g -fopenmp -fPIC $(LIBRARY)

CFLAGS_RELEASE = -Ofast -fopenmp -Wall -std=c++11 -fPIC $(INCLUDE) -DNDEBUG -DGIT_COMMIT=\"$(GIT_COMMIT)\"
LFLAGS_RELEASE = -Ofast -fopenmp -fPIC $(LIBRARY)

# CFLAGS = $(CFLAGS_DEBUG)
# LFLAGS = $(LFLAGS_DEBUG)
CFLAGS = $(CFLAGS_RELEASE)
LFLAGS = $(LFLAGS_RELEASE)

FILES = amplitudes arrays gauss bspline hdffile hex input \
	matrix opencl preconditioners radial slater specf version

SOURCES = $(patsubst %,src/%.cpp,$(FILES))
OBJECTS = $(patsubst %,obj/%.o,$(FILES))
DEPENDENCIES = $(patsubst %,obj/%.dep,$(FILES))

all: pollversion bin/hex-ecs

pollversion:
	@touch src/version.cpp

doc:
	doxygen

clean:
	rm -f bin/hex $(OBJECTS) $(DEPENDENCIES)

bin/hex-ecs: $(OBJECTS)
	@mkdir -p $(@D)
	$(CPP) $^ $(LFLAGS) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

-include $(wildcard obj/*.dep)
