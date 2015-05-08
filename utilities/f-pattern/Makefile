all: ../bin/hex-f-pattern

../bin/hex-f-pattern: hex-f-pattern.o special.o misc.o
	g++ $^ -o $@ -lgsl -lgslcblas

hex-f-pattern.o: hex-f-pattern.cpp arrays.h arrithm.h misc.h numbers.h special.h
	g++ -Wall -std=c++11 -DNO_HDF -DNO_MPI -c $< -o $@

special.o: special.cpp arrays.h arrithm.h misc.h numbers.h special.h
	g++ -Wall -std=c++11 -DNO_HDF -DNO_MPI -DNO_LAPACK -c $< -o $@

misc.o: misc.cpp misc.h
	g++ -Wall -std=c++11 -DNO_HDF -DNO_MPI -DNO_LAPACK -c $< -o $@
