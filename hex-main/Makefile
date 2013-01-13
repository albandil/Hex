CPP = g++
MPICPP = mpic++ -DWITH_MPI

SUITESPARSELIBS = /usr/local/lib/libumfpack.a /usr/local/lib/libamd.a /usr/local/lib/libcholmod.a /usr/local/lib/libcolamd.a /usr/local/lib/libcholmod.a /usr/local/lib/libcamd.a /usr/local/lib/libccolamd.a /usr/local/lib/libmetis.a /usr/local/lib/libsuitesparseconfig.a

INCLUDE = 
LIBRARY = $(SUITESPARSELIBS) -lrt -lhdf5_cpp -lhdf5 -lgsl -lopenblas -lgslcblas -lz -llapack -lpthread -lm -lpng

CFLAGS = -O3 -fopenmp -Wall -std=c++0x -fPIC $(INCLUDE) #-DHEX_USE_LU_SCRATCH
LFLAGS = -O3 -fopenmp -fPIC $(LIBRARY)

all: bin/hex bin/hex-mpi

docu:
	doxygen

clean:
	rm -f bin/hex obj/*.o

bin/hex: obj/hex.o obj/amplitudes.o obj/gauss.o obj/bspline.o obj/input.o obj/moments.o obj/slater.o obj/angs.o obj/spmatrix.o obj/special.o
	$(CPP) $^ $(LFLAGS) -o $@

bin/hex-mpi: obj/hex-mpi.o obj/amplitudes-mpi.o obj/gauss-mpi.o obj/bspline-mpi.o obj/input-mpi.o obj/moments-mpi.o obj/slater-mpi.o obj/angs-mpi.o obj/spmatrix-mpi.o obj/special-mpi.o
	$(MPICPP) $^ $(LFLAGS) -o $@

obj/hex.o: src/hex.cpp src/amplitudes.h src/angs.h src/arrays.h src/bspline.h src/moments.h src/slater.h src/spmatrix.h src/special.h src/input.h
	$(CPP) $(CFLAGS) -c $< -o $@

obj/hex-mpi.o: src/hex-mpi.cpp src/amplitudes.h src/angs.h src/arrays.h src/bspline.h src/moments.h src/slater.h src/spmatrix.h src/special.h src/input.h
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/amplitudes.o: src/amplitudes.cpp src/amplitudes.h src/angs.h src/arrays.h src/bspline.h src/special.h src/spmatrix.h
	$(CPP) $(CFLAGS) -c $< -o $@
	
obj/amplitudes-mpi.o: src/amplitudes.cpp src/amplitudes.h src/angs.h src/arrays.h src/bspline.h src/special.h src/spmatrix.h
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/angs.o: src/angs.cpp
	$(CPP) $(CFLAGS) -c $< -o $@

obj/angs-mpi.o: src/angs.cpp
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/arrays.o: src/arrays.cpp src/complex.h
	$(CPP) $(CFLAGS) -c $< -o $@

obj/arrays-mpi.o: src/arrays.cpp src/complex.h
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/bspline.o: src/bspline.cpp src/arrays.h src/bspline.h
	$(CPP) $(CFLAGS) -c $< -o $@

obj/bspline-mpi.o: src/bspline.cpp src/arrays.h src/bspline.h
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/gauss.o: src/gauss.cpp src/arrays.h src/bspline.h
	$(CPP) $(CFLAGS) -c $< -o $@

obj/gauss-mpi.o: src/gauss.cpp src/arrays.h src/bspline.h
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/input.o: src/input.cpp src/arrays.h
	$(CPP) $(CFLAGS) -c $< -o $@

obj/input-mpi.o: src/input.cpp src/arrays.h
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/moments.o: src/moments.cpp src/arrays.h src/bspline.h src/gauss.h
	$(CPP) $(CFLAGS) -c $< -o $@

obj/moments-mpi.o: src/moments.cpp src/arrays.h src/bspline.h src/gauss.h
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/slater.o: src/slater.cpp src/arrays.h src/bspline.h src/moments.h src/gauss.h
	$(CPP) $(CFLAGS) -c $< -o $@

obj/slater-mpi.o: src/slater.cpp src/arrays.h src/bspline.h src/moments.h src/gauss.h
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/special.o: src/special.cpp src/arrays.h
	$(CPP) $(CFLAGS) -c $< -o $@

obj/special-mpi.o: src/special.cpp src/arrays.h
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/spmatrix.o: src/spmatrix.cpp src/arrays.h src/spmatrix.h
	$(CPP) $(CFLAGS) -c $< -o $@

obj/spmatrix-mpi.o: src/spmatrix.cpp src/arrays.h src/spmatrix.h
	$(MPICPP) $(CFLAGS) -c $< -o $@
