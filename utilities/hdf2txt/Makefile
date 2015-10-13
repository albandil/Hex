SOURCES = hdf2txt.cpp \
         ../../common/hex-arrays.cpp \
         ../../common/hex-hdffile.cpp \
         ../../common/hex-misc.cpp

HEADERS = ../../common/hex-arrithm.h  \
         ../../common/hex-arrays.h \
         ../../common/hex-numbers.h \
         ../../common/hex-misc.h

.PHONY: all clean

all: ../bin/hex-hdf2txt

../bin/hex-hdf2txt: $(SOURCES) $(HEADERS)
	$(CXX) -std=c++11 -O2 -I../../common $(SOURCES) -o ../bin/hex-hdf2txt -lstdc++ -lm

clean:
	rm -rf ../bin/hex-hdf2txt
