SOURCES = chebeval.cpp \
          ../../common/hex-arrays.cpp \
          ../../common/hex-hdffile.cpp \
          ../../common/hex-misc.cpp
HEADERS = ../../common/hex-misc.h \
          ../../common/hex-arrithm.h \
          ../../common/hex-arrays.h \
          ../../common/hex-numbers.h \
          ../../common/hex-compact.h \
          ../../common/hex-chebyshev.h \
          ../../common/hex-hdffile.h

HEADERS = 

.PHONY: all clean

all: ../bin/hex-chebeval

../bin/hex-chebeval: $(SOURCES) $(HEADERS)
	$(CXX) -std=c++11 -O2 -I../../common $(SOURCES) -o ../bin/hex-chebeval -lstdc++ -lm

clean:
	rm -rf ../bin/hex-chebeval
