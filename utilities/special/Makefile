.PHONY: all

all: hex-special

clean:
	rm -rf hex-special

hex-special: hex-special.cpp ../../src/special.cpp ../../src/special.h
	g++ -O2 -Wall -DNO_LAPACK -std=c++11 $< -o $@ -lgsl -lgslcblas
