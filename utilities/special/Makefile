.PHONY: all

all: ../bin/hex-special

clean:
	rm -rf ../bin/hex-special

../bin/hex-special: hex-special.cpp ../../src/special.cpp ../../src/special.h ../../src/misc.cpp
	$(CXX) -O2 -Wall -DNO_LAPACK -std=c++11 $< -o $@ -lgsl -lgslcblas -lstdc++ -lm
