.PHONY: all clean doc pollversion dist distclean

##############################################################
# Editable part

# compiler executable
CXX = g++
# CXX = icpc
                                                             #
##############################################################

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CFLAGS = -Ofast -Wall -std=c++11 -DNO_HDF -DGIT_COMMIT=\"$(GIT_COMMIT)\" -DNDEBUG
LFLAGS = -lgsl -lgslcblas -lcln

FILES = pwba symbolic specf arrays version
SOURCES = $(patsubst %,src/%.cpp,$(FILES))
OBJECTS = $(patsubst %,obj/%.o,$(FILES))
DEPENDENCIES = $(patsubst %,obj/%.dep,$(FILES))

all: bin/hex-pwba

dist:
	mkdir -p hex-pwba
	cp --recursive --dereference README Makefile Doxyfile src hex-pwba/
	sed -i "s/\(^GIT_COMMIT = \)\(.*\)/\1$(GIT_COMMIT)/g" hex-pwba/Makefile
	sed -i "s/\(^all:.*\)\(pollversion\)\(.*\)/\1\3/g" hex-pwba/Makefile
	tar --create --gzip --file=hex-pwba-$(GIT_COMMIT).tar.gz hex-pwba
	rm -rf hex-pwba/

distclean:
	rm -rf hex-dwba hex-dwba-*.tar.gz

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
