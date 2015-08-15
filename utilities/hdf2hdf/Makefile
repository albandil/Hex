.PHONY: all clean

all: ../bin/hex-hdf2hdf

../bin/hex-hdf2hdf: hdf2hdf.cpp ../../src/arrithm.h ../../src/arrays.h ../../src/arrays.cpp ../../src/numbers.h ../../src/hdffile.cpp ../../src/hdffile.h ../../src/h5file.cpp ../../src/h5file.h
	$(CXX) -std=c++11 -O2 hdf2hdf.cpp -o ../bin/hex-hdf2hdf -lhdf5 -DNO_MPI -lstdc++ -lm

clean:
	rm -rf ../bin/hex-hdf2hdf
