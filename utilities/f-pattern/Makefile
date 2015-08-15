all: ../bin/hex-f-pattern

clean:
	rm -f *.o ../bin/hex-f-pattern

../bin/hex-f-pattern: hex-f-pattern.o special.o misc.o
	$(CXX) $^ -o $@ -lgsl -lgslcblas -lstdc++ -lm

hex-f-pattern.o: hex-f-pattern.cpp arrays.h arrithm.h misc.h numbers.h special.h
	$(CXX) -Wall -std=c++11 -DNO_HDF -DNO_MPI -c $< -o $@

special.o: special.cpp arrays.h arrithm.h misc.h numbers.h special.h
	$(CXX) -Wall -std=c++11 -DNO_HDF -DNO_MPI -DNO_LAPACK -c $< -o $@

misc.o: misc.cpp misc.h
	$(CXX) -Wall -std=c++11 -DNO_HDF -DNO_MPI -DNO_LAPACK -c $< -o $@
