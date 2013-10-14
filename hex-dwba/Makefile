CXXDEBUG = -g -fPIC -std=c++11 -Wall
CXXRILIS = -O3 -fPIC -std=c++11 -Wall -DNDEBUG

LIB = -lrt -lopenblas -lo2scl -lgsl -lgslcblas -lhdf5 -lhdf5_cpp -lcln -lfftw3

LDEBUG = -g -fPIC $(LIB)
LRILIS = -O3 -fPIC $(LIB)

CXXFLAGS = $(CXXDEBUG)
LFLAGS = $(LDEBUG)
# CXXFLAGS = $(CXXRILIS)
# LFLAGS = $(LRILIS)

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
