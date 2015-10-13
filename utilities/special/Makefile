SOURCES = hex-special.cpp \
          ../../common/hex-special.cpp \
          ../../common/hex-misc.cpp

HEADERS = ../../common/hex-special.h \
          ../../common/hex-misc.h

.PHONY: all

all: ../bin/hex-special

clean:
	rm -rf ../bin/hex-special

../bin/hex-special: $(SOURCES) $(HEADERS)
	$(CXX) -O2 -Wall -DNO_LAPACK -I../../common -std=c++11 $(SOURCES) -o $@ -lgsl -lgslcblas -lstdc++ -lm
