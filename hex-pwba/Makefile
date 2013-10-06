CXX = g++
CFLAGS = -Ofast -Wall -std=c++11
LFLAGS = -lgsl -lgslcblas -lo2scl -lcln -lhdf5 -lhdf5_cpp

FILES = pwba symbolic specf arrays hdffile
SOURCES = $(patsubst %,src/%.cpp,$(FILES))
OBJECTS = $(patsubst %,obj/%.o,$(FILES))
DEPENDENCIES = $(patsubst %,obj/%.dep,$(FILES))

all: bin/pwba

bin/pwba: $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(LFLAGS) $^ -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

clean:
	rm -f $(OBJECTS) $(DEPENDENCIES) bin/pwba
