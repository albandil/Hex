.PHONY: all clean

all: hdf2txt

hdf2txt: hdf2txt.cpp ../../src/arrithm.h  ../../src/arrays.cpp ../../src/arrays.h ../../src/complex.h ../../src/hdffile.cpp
	g++ -std=c++11 -O2 hdf2txt.cpp -o hdf2txt -lhdf5 -lhdf5_cpp

clean:
	rm -rf hdf2txt
