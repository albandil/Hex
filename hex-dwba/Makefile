.PHONY: all clean allclean doc docclean pollversion dist distclean


################################################################################
# Editable section                                                             #

# ---------------------------------------------------------------------------- #
# Compiler and linker + flags                                                  #
# ---------------------------------------------------------------------------- #

CXX = g++

CFLAGS = -Wall -std=c++11 -fPIC
LFLAGS = 

# ---------------------------------------------------------------------------- #
# HDF (optional, use -DNO_HDF to omit)                                         #
# ---------------------------------------------------------------------------- #

HDFFLAGS = -DNO_HDF
# HDFLIBS = -lhdf5 -lhdf5_cpp

# ---------------------------------------------------------------------------- #
# GSL                                                                          #
# ---------------------------------------------------------------------------- #

GSLFLAGS = $(shell gsl-config --cflags)
GSLLIBS = $(shell gsl-config --libs)

# ---------------------------------------------------------------------------- #
# CLN                                                                          #
# ---------------------------------------------------------------------------- #

CLNFLAGS =
CLNLIBS = -lcln
                                                                               #
################################################################################


CFLAGS += $(HDFFLAGS) $(GSLFLAGS) $(CLNFLAGS)
LFLAGS += $(HDFLIBS) $(GSLLIBS) $(CLNLIBS)

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CXXRILIS = $(CFLAGS) -O3 -DNDEBUG -DNO_LAPACK -DGIT_COMMIT=\"$(GIT_COMMIT)\"
CXXDEBUG = $(CFLAGS) -g -DNO_LAPACK -DGIT_COMMIT=\"$(GIT_COMMIT)\"

# CXXFLAGS = $(CXXDEBUG)
CXXFLAGS = $(CXXRILIS)

FILES = arrays dwba hdffile hydrogen potential pwba special symbolic ui \
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
	rm -f $(OBJECTS) $(DEPENDENCIES) bin/hex-dwba

doc:
	doxygen

docclean:
	rm -rf html

-include $(wildcard obj/*.dep)
