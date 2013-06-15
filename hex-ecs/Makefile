CPP = g++

SUITESPARSELIBS = /usr/local/lib/libumfpack.a /usr/local/lib/libamd.a /usr/local/lib/libcholmod.a /usr/local/lib/libcolamd.a /usr/local/lib/libcholmod.a /usr/local/lib/libcamd.a /usr/local/lib/libccolamd.a /usr/local/lib/libmetis.a /usr/local/lib/libsuitesparseconfig.a
OPENMPILIBS = -L/usr/lib64/mpi/gcc/openmpi/lib64
OPENMPIINCLUDE = -I/usr/lib64/mpi/gcc/openmpi/include

INCLUDE = $(OPENMPIINCLUDE)
LIBRARY = $(SUITESPARSELIBS) $(OPENMPILIBS) -lmpi -lmpi_cxx -lrt -lhdf5_cpp -lhdf5 -lgsl -lopenblas -lgslcblas -lz -lpthread -lm -lpng -lfftw3

CFLAGS = -O3 -fopenmp -Wall -std=c++11 -fPIC $(INCLUDE) #-DHEX_USE_LU_SCRATCH
LFLAGS = -O3 -fopenmp -fPIC $(LIBRARY)

all: bin/hex 

docu:
	doxygen

clean:
	rm -f bin/hex obj/*.o

bin/hex: obj/hex.o obj/amplitudes.o obj/gauss.o obj/bspline.o obj/input.o obj/moments.o obj/slater.o obj/angs.o obj/spmatrix.o obj/special.o obj/gauss_weights.o obj/gauss_roots.o obj/arrays.o
	@mkdir -p $(@D)
	$(CPP) $^ $(LFLAGS) -o $@

obj/hex.o: src/hex.cpp src/amplitudes.h src/angs.h src/arrays.h src/bspline.h src/moments.h src/slater.h src/spmatrix.h src/special.h src/input.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/amplitudes.o: src/amplitudes.cpp src/amplitudes.h src/angs.h src/arrays.h src/bspline.h src/special.h src/spmatrix.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@
	
obj/angs.o: src/angs.cpp
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/arrays.o: src/arrays.cpp src/complex.h src/arrays.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/bspline.o: src/bspline.cpp src/arrays.h src/bspline.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/gauss.o: src/gauss.cpp src/arrays.h src/bspline.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@
	
obj/gauss_weights.o: src/gauss_weights.cpp src/arrays.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@
	
obj/gauss_roots.o: src/gauss_roots.cpp src/arrays.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/input.o: src/input.cpp src/arrays.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/moments.o: src/moments.cpp src/arrays.h src/bspline.h src/gauss.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/slater.o: src/slater.cpp src/arrays.h src/bspline.h src/moments.h src/gauss.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/special.o: src/special.cpp src/arrays.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/spmatrix.o: src/spmatrix.cpp src/arrays.h src/spmatrix.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@ #-ftree-vectorizer-verbose=5

