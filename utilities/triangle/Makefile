.PHONY: all clean

CFLAGS = -std=c++11
LFLAGS = -lgsl -lgslcblas

all: triangle

triangle: triangle.o angs.o
	g++ $^ -o $@ $(LFLAGS)

triangle.o: triangle.cpp angs.h
	g++ -c $< -o $@ $(CFLAGS)

angs.o: angs.cpp angs.h
	g++ -c $< -o $@ $(CFLAGS)

clean:
	rm -rf triangle *.o
