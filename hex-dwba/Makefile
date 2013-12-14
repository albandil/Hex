.PHONY: all clean doc pollversion

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
LIBRARY = -lrt -lo2scl -lgsl -lgslcblas -lfftw3 $(HDFLIBS)

# CXXFLAGS = $(CXXDEBUG)
CXXFLAGS = $(CXXRILIS)
LFLAGS = -fPIC $(LIBRARY)

FILES = arrays dwba1 hdffile hydrogen main potential specf \
        version wave_distort
SOURCES = $(patsubst %,src/%.cpp,$(FILES))
OBJECTS = $(patsubst %,obj/%.o,$(FILES))
DEPENDENCIES = $(patsubst %,obj/%.dep,$(FILES))

all: pollversion bin/hex-dwba

pollversion:
	touch src/version.cpp

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
