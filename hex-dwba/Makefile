.PHONY: all clean doc

##############################################################
# Editable part

# HDF library
HDFFLAGS = -DNO_HDF
# HDFINCLUDE = 
# HDFLIBS = -lhdf5 -lhdf5_cpp
                                                              #
###############################################################

CXXDEBUG = -g -fPIC -std=c++11 -Wall $(HDFINCLUDE) $(HDFFLAGS)
CXXRILIS = -O3 -fPIC -std=c++11 -Wall $(HDFINCLUDE) $(HDFFLAGS) -DNDEBUG
LIBRARY = -lrt -lo2scl -lgsl -lgslcblas -lfftw3 $(HDFLIBS)

# CXXFLAGS = $(CXXDEBUG)
CXXFLAGS = $(CXXRILIS)
LFLAGS = -fPIC $(LIBRARY)

FILES = arrays dwba1 hdffile hydrogen main potential specf wave_distort
SOURCES = $(patsubst %,src/%.cpp,$(FILES))
OBJECTS = $(patsubst %,obj/%.o,$(FILES))
DEPENDENCIES = $(patsubst %,obj/%.dep,$(FILES))

all: bin/hex-dwba

bin/hex-dwba: $(OBJECTS)
	@mkdir -p $(@D)
	g++ $^ $(LFLAGS) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	g++ -c $(CXXFLAGS) $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

clean:
	rm -f $(OBJECTS) $(DEPENDENCIES) bin/dwba

doc:
	doxygen
