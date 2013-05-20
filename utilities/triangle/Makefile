.PHONY: all clean

LFLAGS = -lgsl -lgslcblas

all: triangle

triangle: triangle.o angs.o
	g++ $^ -o $@ $(LFLAGS)

triangle.o: triangle.cpp angs.h
	g++ -c $< -o $@

angs.o: angs.cpp angs.h
	g++ -c $< -o $@

clean:
	rm -rf triangle *.o
