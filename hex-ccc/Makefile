.PHONY: all clean allclean dist distclean doc docclean

##############################################################
# Editable part

CPP = g++
# CPP = icpc

##############################################################

CFLAGS = -Wall -std=c++11 -fPIC -fopenmp
LFLAGS =  -fopenmp -lgsl -lopenblas -lgslcblas -lcln

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CXXRILIS = $(CFLAGS) -DNDEBUG -DNO_HDF -DNO_UMFPACK -DNO_PNG -O2 -DGIT_COMMIT=\"$(GIT_COMMIT)\"
CXXDEBUG = $(CFLAGS) -g -DNO_HDF -DNO_UMFPACK -DNO_PNG -O2 -DGIT_COMMIT=\"$(GIT_COMMIT)\"

CXXFLAGS = $(CXXDEBUG)
# CXXFLAGS = $(CXXRILIS)

FILES = arrays basis gauss hdffile main matrix moments potential quadrature specf symbolic version
SOURCES = $(patsubst %,src/%.cpp,$(FILES))
OBJECTS = $(patsubst %,obj/%.o,$(FILES))
DEPENDENCIES = $(patsubst %,obj/%.dep,$(FILES))

all: pollversion bin/hex-ccc

pollversion:
	@touch src/version.cpp

doc:
	doxygen

docclean:
	rm -rf html

dist:
	mkdir -p hex-ccc
	cp --recursive --dereference README Makefile Doxyfile *.cbp ../media/hexe-small.gif src hex-ccc/
	sed -i "s/\(^GIT_COMMIT = \)\(.*\)/\1$(GIT_COMMIT)/g" hex-ccc/Makefile
	sed -i "s/\(^all:.*\)\(pollversion\)\(.*\)/\1\3/g" hex-ccc/Makefile
	sed -i "s/\.\.\/media\///g" hex-ccc/Doxyfile
	tar --create --gzip --file=hex-ccc-$(GIT_COMMIT).tar.gz hex-ccc
	rm -rf hex-ccc/

clean:
	rm -f obj/*.o obj/*.dep bin/hex-ccc

distclean:
	rm -rf hex-ccc hex-ccc-*.tar.gz

allclean: clean distclean docclean

bin/hex-ccc: $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $^ $(LFLAGS) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CPP) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

-include $(wildcard obj/*.dep)
