.PHONY: all clean doc

##############################################################
# Editable part

# compiler executable
CXX = g++
# CXX = icpc
                                                             #
##############################################################

CFLAGS = -Ofast -Wall -std=c++11 -DNO_HDF
LFLAGS = -lgsl -lgslcblas -lo2scl -lcln

FILES = pwba symbolic specf arrays
SOURCES = $(patsubst %,src/%.cpp,$(FILES))
OBJECTS = $(patsubst %,obj/%.o,$(FILES))
DEPENDENCIES = $(patsubst %,obj/%.dep,$(FILES))

all: bin/hex-pwba

bin/hex-pwba: $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(LFLAGS) $^ -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

clean:
	rm -f $(OBJECTS) $(DEPENDENCIES) bin/pwba

doc:
	doxygen
