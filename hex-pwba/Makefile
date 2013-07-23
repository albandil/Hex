CXX = g++
CFLAGS = -Ofast -Wall -std=c++11
LFLAGS = -lgsl -lgslcblas -lo2scl -lcln

all: bin/pwba

bin/pwba: obj/pwba.o obj/symbolic.o obj/angs.o obj/specf.o
	@mkdir -p $(@D)
	$(CXX) $(LFLAGS) $^ -o $@

obj/pwba.o: src/pwba.cpp src/gausskronrod.h src/symbolic.h
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -c $< -o $@

obj/symbolic.o: src/symbolic.cpp src/symbolic.h
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -c $< -o $@

obj/angs.o: src/angs.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -c $< -o $@
	
obj/specf.o: src/specf.cpp src/specf.h src/complex.h
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -c $< -o $@

clean:
	rm -f obj/*.o bin/pwba

