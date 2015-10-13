SOURCES = hdf2hdf.cpp \
          ../../common/hex-arrays.cpp \
          ../../common/hex-hdffile.cpp \
          ../../common/hex-h5file.cpp \
          ../../common/hex-misc.cpp

HEADERS = ../../common/hex-arrithm.h \
          ../../common/hex-arrays.h \
          ../../common/hex-numbers.h \
          ../../common/hex-hdffile.h \
          ../../common/hex-h5file.h \
          ../../common/hex-misc.h

.PHONY: all clean

all: ../bin/hex-hdf2hdf

../bin/hex-hdf2hdf: $(SOURCES) $(HEADERS)
	$(CXX) -std=c++11 -I../../common -O2 $(SOURCES) -o ../bin/hex-hdf2hdf -lhdf5 -DNO_MPI -lstdc++ -lm

clean:
	rm -rf ../bin/hex-hdf2hdf
