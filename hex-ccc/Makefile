.PHONY: all clean allclean dist distclean doc docclean

##############################################################
# Editable part

CPP = g++
# CPP = icpc

##############################################################

CFLAGS = -Wall -std=c++11 -fPIC
LFLAGS = -lgsl -lgslcblas

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CXXRILIS = $(CFLAGS) -DNDEBUG -DNO_HDF -O2 -DGIT_COMMIT=\"$(GIT_COMMIT)\"
CXXDEBUG = $(CFLAGS) -DNO_HDF -g -O2 -DGIT_COMMIT=\"$(GIT_COMMIT)\"

CXXFLAGS = $(CXXDEBUG)
# CXXFLAGS = $(CXXRILIS)

FILES = arrays hdffile main specf version
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
