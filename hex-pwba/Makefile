CXX = g++
CFLAGS = -Ofast -Wall -std=c++0x
LFLAGS = -lgsl -lgslcblas -lo2scl -lcln

all: bin/pwba

bin/pwba: obj/pwba.o obj/symbolic.o
	$(CXX) $(LFLAGS) $^ -o $@

obj/pwba.o: src/pwba.cpp src/integrate.h src/symbolic.h
	$(CXX) $(CFLAGS) -c $< -o $@

obj/symbolic.o: src/symbolic.cpp src/symbolic.h
	$(CXX) $(CFLAGS) -c $< -o $@

clean:
	rm -f obj/*.o bin/pwba

