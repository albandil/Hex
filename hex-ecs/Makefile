.PHONY: all doc clean pollversion

###################################################################################
# Editable section

# Compiler executable
CPP = g++

# Suite Sparse environment
SUITESPARSENAMES = umfpack amd cholmod colamd camd ccolamd metis suitesparseconfig
SUITESPARSELIBS = $(patsubst %,-l%,$(SUITESPARSENAMES))

# OpenMPI environment
OPENMPIFLAGS = -DNO_MPI
# OPENMPIINCLUDE = -I/usr/lib64/mpi/gcc/openmpi/include
# OPENMPILIBS = -L/usr/lib64/mpi/gcc/openmpi/lib64

# AMD OpenCL environment
OPENCLFLAGS = -DNO_OPENCL
# OPENCLINCLUDE = -I/opt/AMDAPP/include
# OPENCLLIBS = -lOpenCL

# Intel OpenCL environment
# OPENCLINCLUDE = -I/opt/intel/opencl-1.2-3.1.1.11385/include
# OPENCLLIBS = -L/opt/intel/opencl-1.2-3.1.1.11385/lib64 -lOpenCL

# Lapack environment
LAPACKFLAGS = -DNO_LAPACK
# LAPACKINCLUDE =
# LAPACKLIBS = -lgfortran -llapack

# PNG++ environment
PNGFLAGS = -DNO_PNG
# PNGPPINCLUDE = -I./libs
# PNGPPLIBS = -lpng
                                                                                    #
#####################################################################################

INCLUDE = $(OPENMPIINCLUDE) $(OPENCLINCLUDE) $(LAPACKINCLUDE) $(PNGPPINCLUDE)
LIBRARY = $(SUITESPARSELIBS) $(OPENMPILIBS) \
	-lmpi -lmpi_cxx -lrt -lhdf5_cpp -lhdf5 -lgsl -lopenblas -lgslcblas -lpthread \
	-lm -lfftw3 $(OPENCLLIBS) $(LAPACKLIBS) $(PNGPPLIBS)

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CFLAGS_DEBUG = -g -fopenmp -Wall -std=c++11 -fPIC $(INCLUDE) -DGIT_COMMIT=\"$(GIT_COMMIT)\"
CFLAGS_RELEASE = -Ofast -fopenmp -Wall -std=c++11 -fPIC $(INCLUDE) -DNDEBUG -DGIT_COMMIT=\"$(GIT_COMMIT)\"

LFLAGS = -fopenmp -fPIC $(LIBRARY)

# CFLAGS = $(CFLAGS_DEBUG) $(LAPACKFLAGS) $(OPENCLFLAGS) $(PNGFLAGS)
CFLAGS = $(CFLAGS_RELEASE) $(LAPACKFLAGS) $(OPENMPIFLAGS) $(OPENCLFLAGS) $(PNGFLAGS)

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

clean:
	rm -f bin/hex $(OBJECTS) $(DEPENDENCIES)

bin/hex-ecs: $(OBJECTS)
	@mkdir -p $(@D)
	$(CPP) $^ $(LFLAGS) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

src/kernels_cl.c: src/kernels.cl
	xxd -i src/kernels.cl | tail -n+2 | head -n-2 > src/kernels_cl.c

-include $(wildcard obj/*.dep)
