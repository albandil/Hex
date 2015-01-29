.PHONY: all clean

all: ../bin/hex-chebeval

../bin/hex-chebeval: chebeval.cpp ../../src/arrithm.h ../../src/arrays.cpp ../../src/arrays.h ../../src/numbers.h ../../src/compact.h ../../src/misc.h ../../src/chebyshev.h ../../src/hdffile.cpp
	g++ -std=c++11 -O2 chebeval.cpp -o ../bin/hex-chebeval -lhdf5 -lhdf5_cpp -DNO_MPI

clean:
	rm -rf ../bin/hex-chebeval
