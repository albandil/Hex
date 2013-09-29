.PHONY: all clean

all: hdf2hdf

hdf2hdf: hdf2hdf.cpp ../../src/arrithm.h ../../src/arrays.h ../../src/arrays.cpp ../../src/complex.h ../../src/hdffile.cpp
	g++ -std=c++11 -O2 hdf2hdf.cpp -o hdf2hdf -lhdf5 -lhdf5_cpp

clean:
	rm -rf hdf2hdf
