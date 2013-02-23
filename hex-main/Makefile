CPP = g++
MPICPP = mpic++ -DWITH_MPI

SUITESPARSELIBS = /usr/local/lib/libumfpack.a /usr/local/lib/libamd.a /usr/local/lib/libcholmod.a /usr/local/lib/libcolamd.a /usr/local/lib/libcholmod.a /usr/local/lib/libcamd.a /usr/local/lib/libccolamd.a /usr/local/lib/libmetis.a /usr/local/lib/libsuitesparseconfig.a

INCLUDE = 
LIBRARY = $(SUITESPARSELIBS) -lrt -lhdf5_cpp -lhdf5 -lgsl -lopenblas -lgslcblas -lo2scl -lz -llapack -lpthread -lm -lpng

CFLAGS = -O3 -fopenmp -Wall -std=c++0x -fPIC $(INCLUDE) #-DHEX_USE_LU_SCRATCH
LFLAGS = -O3 -fopenmp -fPIC $(LIBRARY)

all: bin/hex 

docu:
	doxygen

clean:
	rm -f bin/hex obj/*.o

bin/hex: obj/hex.o obj/amplitudes.o obj/gauss.o obj/bspline.o obj/input.o obj/moments.o obj/slater.o obj/angs.o obj/spmatrix.o obj/special.o obj/gauss_weights.o obj/gauss_roots.o
	@mkdir -p $(@D)
	$(CPP) $^ $(LFLAGS) -o $@

bin/hex-mpi: obj/hex-mpi.o obj/amplitudes-mpi.o obj/gauss-mpi.o obj/bspline-mpi.o obj/input-mpi.o obj/moments-mpi.o obj/slater-mpi.o obj/angs-mpi.o obj/spmatrix-mpi.o obj/special-mpi.o
	@mkdir -p $(@D)
	$(MPICPP) $^ $(LFLAGS) -o $@

obj/hex.o: src/hex.cpp src/amplitudes.h src/angs.h src/arrays.h src/bspline.h src/moments.h src/slater.h src/spmatrix.h src/special.h src/input.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/hex-mpi.o: src/hex-mpi.cpp src/amplitudes.h src/angs.h src/arrays.h src/bspline.h src/moments.h src/slater.h src/spmatrix.h src/special.h src/input.h
	@mkdir -p $(@D)
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/amplitudes.o: src/amplitudes.cpp src/amplitudes.h src/angs.h src/arrays.h src/bspline.h src/special.h src/spmatrix.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@
	
obj/amplitudes-mpi.o: src/amplitudes.cpp src/amplitudes.h src/angs.h src/arrays.h src/bspline.h src/special.h src/spmatrix.h
	@mkdir -p $(@D)
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/angs.o: src/angs.cpp
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/angs-mpi.o: src/angs.cpp
	@mkdir -p $(@D)
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/arrays.o: src/arrays.cpp src/complex.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/arrays-mpi.o: src/arrays.cpp src/complex.h
	@mkdir -p $(@D)
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/bspline.o: src/bspline.cpp src/arrays.h src/bspline.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/bspline-mpi.o: src/bspline.cpp src/arrays.h src/bspline.h
	@mkdir -p $(@D)
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/gauss.o: src/gauss.cpp src/arrays.h src/bspline.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@
	
obj/gauss_weights.o: src/gauss_weights.cpp src/arrays.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@
	
obj/gauss_roots.o: src/gauss_roots.cpp src/arrays.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/gauss-mpi.o: src/gauss.cpp src/arrays.h src/bspline.h
	@mkdir -p $(@D)
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/input.o: src/input.cpp src/arrays.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/input-mpi.o: src/input.cpp src/arrays.h
	@mkdir -p $(@D)
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/moments.o: src/moments.cpp src/arrays.h src/bspline.h src/gauss.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/moments-mpi.o: src/moments.cpp src/arrays.h src/bspline.h src/gauss.h
	@mkdir -p $(@D)
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/slater.o: src/slater.cpp src/arrays.h src/bspline.h src/moments.h src/gauss.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/slater-mpi.o: src/slater.cpp src/arrays.h src/bspline.h src/moments.h src/gauss.h
	@mkdir -p $(@D)
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/special.o: src/special.cpp src/arrays.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/special-mpi.o: src/special.cpp src/arrays.h
	@mkdir -p $(@D)
	$(MPICPP) $(CFLAGS) -c $< -o $@

obj/spmatrix.o: src/spmatrix.cpp src/arrays.h src/spmatrix.h
	@mkdir -p $(@D)
	$(CPP) $(CFLAGS) -c $< -o $@

obj/spmatrix-mpi.o: src/spmatrix.cpp src/arrays.h src/spmatrix.h
	@mkdir -p $(@D)
	$(MPICPP) $(CFLAGS) -c $< -o $@
