.PHONY: all

all: ../bin/hex-special

clean:
	rm -rf ../bin/hex-special

../bin/hex-special: hex-special.cpp ../../src/special.cpp ../../src/special.h
	g++ -O2 -Wall -DNO_LAPACK -std=c++11 $< -o $@ -lgsl -lgslcblas -DNO_MPI
