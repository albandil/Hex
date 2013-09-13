.PHONY: all docu clean obj/version.o

CPP = g++

SUITESPARSELIBS = /usr/local/lib/libumfpack.a /usr/local/lib/libamd.a \
	/usr/local/lib/libcholmod.a /usr/local/lib/libcolamd.a \
	/usr/local/lib/libcamd.a /usr/local/lib/libccolamd.a \
	/usr/local/lib/libmetis.a /usr/local/lib/libsuitesparseconfig.a

OPENMPILIBS = -L/usr/lib64/mpi/gcc/openmpi/lib64
OPENMPIINCLUDE = -I/usr/lib64/mpi/gcc/openmpi/include

INCLUDE = $(OPENMPIINCLUDE)
LIBRARY = $(SUITESPARSELIBS) $(OPENMPILIBS) -lmpi -lmpi_cxx -lrt \
	-lhdf5_cpp -lhdf5 -lgsl -lopenblas -lgslcblas -lpthread \
	-lm -lfftw3 -lgfortran -lpng

CFLAGS_DEBUG = -g -fopenmp -Wall -std=c++11 -fPIC $(INCLUDE) #-DWITH_PNGPP
LFLAGS_DEBUG = -g -fopenmp -fPIC $(LIBRARY)

CFLAGS_RELEASE = -Ofast -fopenmp -Wall -std=c++11 -fPIC $(INCLUDE) -DNDEBUG
LFLAGS_RELEASE = -Ofast -fopenmp -fPIC $(LIBRARY)

CFLAGS = $(CFLAGS_DEBUG)
LFLAGS = $(LFLAGS_DEBUG)

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

all: bin/hex 

docu:
	doxygen

clean:
	rm -f bin/hex obj/*.o

bin/hex: obj/hex.o obj/amplitudes.o obj/gauss.o obj/bspline.o obj/input.o obj/moments.o obj/slater.o obj/spmatrix.o obj/specf.o obj/arrays.o obj/hdffile.o obj/version.o
	@mkdir -p $(@D)
	$(CPP) $^ $(LFLAGS) -o $@

obj/hex.o: src/hex.cpp src/amplitudes.h src/arrays.h src/bspline.h src/moments.h src/slater.h src/spmatrix.h src/specf.h src/input.h src/hdffile.h src/version.h src/itersolve.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/amplitudes.o: src/amplitudes.cpp src/amplitudes.h src/arrays.h src/bspline.h src/specf.h src/spmatrix.h src/hdffile.h src/moments.h src/clenshawcurtis.h src/chebyshev.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@
	
obj/arrays.o: src/arrays.cpp src/complex.h src/arrays.h src/hdffile.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/bspline.o: src/bspline.cpp src/arrays.h src/bspline.h src/hdffile.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/gauss.o: src/gauss.cpp src/arrays.h src/bspline.h src/hdffile.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@
	
obj/input.o: src/input.cpp src/arrays.h src/hdffile.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/moments.o: src/moments.cpp src/arrays.h src/bspline.h src/gauss.h src/hdffile.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/slater.o: src/slater.cpp src/arrays.h src/bspline.h src/moments.h src/gauss.h src/hdffile.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/specf.o: src/specf.cpp src/arrays.h src/hdffile.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/spmatrix.o: src/spmatrix.cpp src/arrays.h src/spmatrix.h src/hdffile.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@ #-ftree-vectorizer-verbose=5

obj/hdffile.o: src/hdffile.cpp src/hdffile.h src/complex.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/version.o:
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c src/version.cpp -o $@ -DGIT_COMMIT=\"$(GIT_COMMIT)\"
	
