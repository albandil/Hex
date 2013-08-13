CXX = g++
CFLAGS = -Ofast -Wall -std=c++11
LFLAGS = -lgsl -lgslcblas -lo2scl -lcln -lhdf5 -lhdf5_cpp

all: bin/pwba

bin/pwba: obj/pwba.o obj/symbolic.o obj/angs.o obj/specf.o obj/arrays.o obj/hdffile.o
	@mkdir -p $(@D)
	$(CXX) $(LFLAGS) $^ -o $@

obj/pwba.o: src/pwba.cpp src/gausskronrod.h src/symbolic.h src/arrays.h src/complex.h src/misc.h src/hdffile.h
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -c $< -o $@

obj/arrays.o: src/arrays.cpp src/arrays.h src/complex.h src/misc.h src/hdffile.h
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -c $< -o $@

obj/symbolic.o: src/symbolic.cpp src/symbolic.h
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -c $< -o $@

obj/specf.o: src/specf.cpp src/specf.h src/complex.h src/misc.h src/arrays.h src/hdffile.h
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -c $< -o $@

obj/hdffile.o: src/hdffile.cpp src/complex.h src/misc.h
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -c $< -o $@

clean:
	rm -f obj/*.o bin/pwba

