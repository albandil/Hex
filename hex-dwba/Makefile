.PHONY: all clean allclean doc docclean pollversion dist distclean

##############################################################
# Editable part

CXX = g++
# CXX = icpc

# HDF library
HDFFLAGS = -DNO_HDF
# HDFINCLUDE = 
# HDFLIBS = -lhdf5 -lhdf5_cpp
                                                              #
###############################################################

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CXXDEBUG = -g -fPIC -std=c++11 -Wall $(HDFINCLUDE) $(HDFFLAGS) -DGIT_COMMIT=\"$(GIT_COMMIT)\"
CXXRILIS = -O3 -fPIC -std=c++11 -Wall $(HDFINCLUDE) $(HDFFLAGS) -DGIT_COMMIT=\"$(GIT_COMMIT)\" -DNDEBUG
LIBRARY = -lrt -lgsl -lgslcblas -lfftw3 $(HDFLIBS)

# CXXFLAGS = $(CXXDEBUG)
CXXFLAGS = $(CXXRILIS)
LFLAGS = -fPIC $(LIBRARY)

FILES = arrays dwba1 hdffile hydrogen main potential specf \
        version wave_distort
SOURCES = $(patsubst %,src/%.cpp,$(FILES))
OBJECTS = $(patsubst %,obj/%.o,$(FILES))
DEPENDENCIES = $(patsubst %,obj/%.dep,$(FILES))

all: pollversion bin/hex-dwba

allclean: clean docclean distclean

pollversion:
	touch src/version.cpp

dist:
	mkdir -p hex-dwba
	cp --recursive --dereference README Makefile Doxyfile *.cbp ../media/hexe-small.gif src hex-dwba/
	sed -i "s/\(^GIT_COMMIT = \)\(.*\)/\1$(GIT_COMMIT)/g" hex-dwba/Makefile
	sed -i "s/\(^all:.*\)\(pollversion\)\(.*\)/\1\3/g" hex-dwba/Makefile
	sed -i "s/\.\.\/media\///g" hex-dwba/Doxyfile
	tar --create --gzip --file=hex-dwba-$(GIT_COMMIT).tar.gz hex-dwba
	rm -rf hex-dwba/

distclean:
	rm -rf hex-dwba hex-dwba-*.tar.gz

bin/hex-dwba: $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $^ $(LFLAGS) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(CXXFLAGS) $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

clean:
	rm -f $(OBJECTS) $(DEPENDENCIES) bin/dwba

doc:
	doxygen

docclean:
	rm -rf html

-include $(wildcard obj/*.dep)
