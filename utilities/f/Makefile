.PHONY: all

all: f

clean:
	rm -rf f

f: f.cpp ../../src/special.cpp ../../src/special.h
	g++ $< -o f -DNO_LAPACK -lgsl -lgslcblas -std=c++11
