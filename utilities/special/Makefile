.PHONY: all

all: special

clean:
	rm -rf special

special: hex-special.cpp ../../src/special.cpp ../../src/special.h
	g++ $< -o special -DNO_LAPACK -lgsl -lgslcblas -std=c++11
