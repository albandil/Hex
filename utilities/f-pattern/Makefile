SOURCES = hex-f-pattern.cpp \
          ../../common/hex-special.cpp \
          ../../common/hex-misc.cpp
          
HEADERS = ../../common/hex-special.h \
          ../../common/hex-misc.cpp

all: ../bin/hex-f-pattern

clean:
	rm -f *.o ../bin/hex-f-pattern

../bin/hex-f-pattern: $(SOURCES) $(HEADERS)
	$(CXX) -std=c++11 -I../../common $(SOURCES) -o $@ -DNO_LAPACK -lgsl -lgslcblas -lstdc++ -lm
