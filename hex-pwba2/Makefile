.PHONY: all clean allclean doc docclean pollversion dist distclean

##############################################################
# Editable part

# compiler executable
CXX = g++
# CXX = icpc
                                                             #
##############################################################

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CFLAGS = -O3 -fopenmp -Wall -std=c++11 -DNO_HDF -DGIT_COMMIT=\"$(GIT_COMMIT)\" -DNDEBUG -DNO_ALIGN
LFLAGS = -O3 -fopenmp -lgsl -lcln -lopenblas

FILES = arrays hydrogen main multipot radial specf symbolic version
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
	tar --create --gzip --file=hex-pwba-$(GIT_COMMIT).tar.gz hex-pwba2
	rm -rf hex-pwba2/

distclean:
	rm -rf hex-pwba2 hex-pwba2-*.tar.gz

bin/hex-pwba2: $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $^ $(LFLAGS) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

clean:
	rm -f $(OBJECTS) $(DEPENDENCIES) bin/pwba2

doc:
	doxygen

docclean:
	rm -rf html

-include $(wildcard obj/*.dep)
