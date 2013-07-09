.PHONY: all

all: f

clean:
	rm -rf f

f: f.cpp ../../src/angs.cpp ../../src/angs.h
	g++ $< -o f -lgsl -lgslcblas -std=c++11
