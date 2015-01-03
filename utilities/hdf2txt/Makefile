.PHONY: all clean

all: ../bin/hex-hdf2txt

../bin/hex-hdf2txt: hdf2txt.cpp ../../src/arrithm.h  ../../src/arrays.cpp ../../src/arrays.h ../../src/complex.h ../../src/hdffile.cpp
	g++ -std=c++11 -O2 hdf2txt.cpp -o ../bin/hex-hdf2txt -lhdf5 -lhdf5_cpp -DNO_MPI

clean:
	rm -rf ../bin/hdf2txt
