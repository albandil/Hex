.PHONY: all clean allclean doc docclean pollversion dist distclean


################################################################################
# Editable section                                                             #

# ---------------------------------------------------------------------------- #
# Compiler and linker + flags                                                  #
# ---------------------------------------------------------------------------- #

CXX = g++

CFLAGS = -Wall -std=c++11 -pedantic-errors -fPIC -fopenmp
LFLAGS = -fopenmp -lstdc++

# ---------------------------------------------------------------------------- #
# GSL                                                                          #
# ---------------------------------------------------------------------------- #

GSLFLAGS = $(shell gsl-config --cflags)
GSLLIBS = $(shell gsl-config --libs)

# ---------------------------------------------------------------------------- #
# GiNaC                                                                        #
# ---------------------------------------------------------------------------- #

GINACFLAGS =
GINACLIBS = -lcln -lginac

# ---------------------------------------------------------------------------- #
# Lapack                                                                       #
# ---------------------------------------------------------------------------- #

LAPACKFLAGS = 
LAPACKLIBS = -lopenblas
                                                                               #
################################################################################


CFLAGS += $(GSLFLAGS) $(GINACFLAGS) $(LAPACKFLAGS)
LFLAGS += $(GSLLIBS) $(GINACLIBS) $(LAPACKLIBS) -rdynamic

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CXXRILIS = $(CFLAGS) -O3 -fcx-limited-range -DNO_HDF -DNO_MPI -DGIT_COMMIT=\"$(GIT_COMMIT)\" -DNDEBUG
CXXDEBUG = $(CFLAGS) -g -fcx-limited-range -DNO_HDF -DNO_MPI -DGIT_COMMIT=\"$(GIT_COMMIT)\"

# CXXFLAGS = $(CXXDEBUG)
CXXFLAGS = -g $(CXXRILIS)

FILES = arrays born hydrogen pwba2-pwx pwba2-spg misc radial special spgrid \
        ui vec3d version hyp_2F1
SOURCES = $(patsubst %,src/%.cpp,$(FILES))
OBJECTS = $(patsubst %,obj/%.o,$(FILES))
DEPENDENCIES = $(patsubst %,obj/%.dep,$(FILES))

all: pollversion bin/hex-pwba2

allclean: clean docclean distclean

pollversion:
	@touch src/version.cpp

dist:
	mkdir -p hex-pwba2
	cp --recursive --dereference README Makefile Doxyfile *.cbp ../media/hexe-small.gif src hex-pwba2/
	sed -i "s/\(^GIT_COMMIT = \)\(.*\)/\1$(GIT_COMMIT)/g" hex-pwba2/Makefile
	sed -i "s/\(^all:.*\)\(pollversion\)\(.*\)/\1\3/g" hex-pwba2/Makefile
	sed -i "s/\.\.\/media\///g" hex-pwba2/Doxyfile
	tar --create --gzip --file=hex-pwba2-$(GIT_COMMIT).tar.gz hex-pwba2
	rm -rf hex-pwba2/

distclean:
	rm -rf hex-pwba2 hex-pwba2-*.tar.gz

bin/hex-pwba2: $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $^ $(LFLAGS) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

clean:
	rm -f $(OBJECTS) $(DEPENDENCIES) bin/pwba2

doc:
	doxygen

docclean:
	rm -rf html

-include $(wildcard obj/*.dep)
