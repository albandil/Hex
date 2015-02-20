.PHONY: all clean

all: ../bin/hex-hdf2hdf

../bin/hex-hdf2hdf: hdf2hdf.cpp ../../src/arrithm.h ../../src/arrays.h ../../src/arrays.cpp ../../src/numbers.h ../../src/hdffile.cpp
	g++ -std=c++11 -O2 hdf2hdf.cpp -o ../bin/hex-hdf2hdf -lhdf5 -DNO_MPI

clean:
	rm -rf ../bin/hex-hdf2hdf
