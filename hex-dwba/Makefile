SUITESPARSELIBS = /usr/local/lib/libumfpack.a /usr/local/lib/libamd.a /usr/local/lib/libcholmod.a /usr/local/lib/libcolamd.a /usr/local/lib/libcholmod.a /usr/local/lib/libcamd.a /usr/local/lib/libccolamd.a /usr/local/lib/libmetis.a /usr/local/lib/libsuitesparseconfig.a

CXXFLAGS = -g -fPIC -O3 -std=c++0x -Wall -fopenmp
LFLAGS = -g -fPIC -O3 $(SUITESPARSELIBS) -lrt -lopenblas -lo2scl -lgsl -lgslcblas -lhdf5 -lhdf5_cpp -fopenmp -lpng -lcln

all: bin/dwba

bin/dwba: obj/main.o obj/potential.o obj/wave_distort.o \
          obj/dwba1.o obj/dwba2.o obj/hydrogen.o obj/arrays.o \
          obj/angs.o obj/specf.o obj/spmatrix.o
	@mkdir -p $(@D)
	g++ $^ $(LFLAGS) -o $@

obj/angs.o: src/angs.cpp src/angs.h
	@mkdir -p $(@D)
	g++ -c $(CXXFLAGS) $< -o $@

obj/main.o: src/main.cpp src/dwba1.h src/complex.h src/arrays.h src/dwba1.h src/dwba2.h src/potential.h src/hydrogen.h
	@mkdir -p $(@D)
	g++ -c $(CXXFLAGS) $< -o $@

obj/potential.o: src/potential.cpp src/potential.h src/arrays.h
	@mkdir -p $(@D)
	g++ -c $(CXXFLAGS) $< -o $@

obj/arrays.o: src/arrays.cpp src/arrays.h src/complex.h
	@mkdir -p $(@D)
	g++ -c $(CXXFLAGS) $< -o $@

obj/spmatrix.o: src/spmatrix.cpp src/spmatrix.h src/arrays.h src/complex.h
	@mkdir -p $(@D)
	g++ -c $(CXXFLAGS) $< -o $@

obj/wave_distort.o: src/wave_distort.cpp src/wave_distort.h src/arrays.h src/complex.h src/potential.h
	@mkdir -p $(@D)
	g++ -c $(CXXFLAGS) $< -o $@
	
# obj/wave_irreg.o: src/wave_irreg.cpp src/wave_irreg.h src/arrays.h src/complex.h
# 	g++ -c $(CXXFLAGS) $< -o $@

# obj/wave_forbid.o: src/wave_forbid.cpp src/wave_forbid.h src/arrays.h src/complex.h
# 	g++ -c $(CXXFLAGS) $< -o $@

# obj/wave_hyperb.o: src/wave_hyperb.cpp src/wave_hyperb.h src/arrays.h src/complex.h
# 	g++ -c $(CXXFLAGS) $< -o $@

obj/dwba1.o: src/dwba1.cpp src/dwba1.h src/complex.h src/integrate.h src/specf.h src/potential.h src/hydrogen.h
	@mkdir -p $(@D)
	g++ -c $(CXXFLAGS) $< -o $@

obj/dwba2.o: src/dwba2.cpp src/dwba2.h src/complex.h src/integrate.h src/specf.h src/potential.h src/angs.h src/spmatrix.h src/hydrogen.h src/chebyshev.h
	@mkdir -p $(@D)
	g++ -c $(CXXFLAGS) $< -o $@

# obj/dwba2_cache.o: src/dwba2_cache.cpp src/dwba2.h src/complex.h src/integrate.h src/specf.h src/potential.h src/angs.h src/spmatrix.h src/hydrogen.h
# 	g++ -c $(CXXFLAGS) $< -o $@

# obj/twoel.o: src/twoel.cpp src/twoel.h src/arrays.h src/complex.h src/integrate.h
# 	g++ -c $(CXXFLAGS) $< -o $@

obj/hydrogen.o: src/hydrogen.cpp src/hydrogen.h
	@mkdir -p $(@D)
	g++ -c $(CXXFLAGS) $< -o $@

obj/specf.o: src/specf.cpp src/specf.h src/complex.h
	@mkdir -p $(@D)
	g++ -c $(CXXFLAGS) $< -o $@

clean:
	rm -f obj/*.o bin/dwba
