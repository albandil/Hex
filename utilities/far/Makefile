.PHONY: all clean

all: far

clean:
	rm -f far

far: far.cpp ../../src/specf.h ../../src/specf.cpp ../../src/hydrogen.h ../../src/hydrogen.cpp
	g++ -std=c++11 far.cpp -o far -lgsl -lgslcblas -lo2scl
