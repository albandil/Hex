.PHONY: all clean

all: chebeval

chebeval: chebeval.cpp ../../src/arrays.cpp ../../src/arrays.h ../../src/complex.h ../../src/compact.h ../../src/misc.h ../../src/chebyshev.h
	g++ -std=c++11 -O2 chebeval.cpp -o chebeval -lhdf5 -lhdf5_cpp

clean:
	rm -rf chebeval
