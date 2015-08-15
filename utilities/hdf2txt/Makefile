.PHONY: all clean

all: ../bin/hex-hdf2txt

../bin/hex-hdf2txt: hdf2txt.cpp ../../src/arrithm.h  ../../src/arrays.cpp ../../src/arrays.h ../../src/numbers.h ../../src/hdffile.cpp
	$(CXX) -std=c++11 -O2 hdf2txt.cpp -o ../bin/hex-hdf2txt -lstdc++ -lm

clean:
	rm -rf ../bin/hex-hdf2txt
