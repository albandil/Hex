.PHONY: all clean allclean doc docclean pollversion dist distclean


################################################################################
# Editable section                                                             #

# ---------------------------------------------------------------------------- #
# Compiler and linker + flags                                                  #
# ---------------------------------------------------------------------------- #

CXX ?= g++

CFLAGS ?= -Wall -std=c++11 -pedantic-errors -fPIC
LFLAGS ?= -lstdc++

# ---------------------------------------------------------------------------- #
# GSL                                                                          #
# ---------------------------------------------------------------------------- #

GSLFLAGS ?= $(shell gsl-config --cflags)
GSLLIBS ?= $(shell gsl-config --libs)

# ---------------------------------------------------------------------------- #
# CLN                                                                          #
# ---------------------------------------------------------------------------- #

CLNFLAGS ?=
CLNLIBS ?= -lcln
                                                                               #
################################################################################


CFLAGS += $(GSLFLAGS) $(CLNFLAGS)
LFLAGS += $(GSLLIBS) $(CLNLIBS)

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CXXRILIS = $(CFLAGS) -O3 -DNDEBUG -DNO_LAPACK -DGIT_COMMIT=\"$(GIT_COMMIT)\"
CXXDEBUG = $(CFLAGS) -g -DNO_LAPACK -DGIT_COMMIT=\"$(GIT_COMMIT)\"

# CXXFLAGS = $(CXXDEBUG)
CXXFLAGS = $(CXXRILIS)

BASECXX = $(shell basename $(CXX))

COMMON_FILES   = hex-arrays hex-hdffile hex-hydrogen hex-misc hex-special hex-symbolic hex-version
COMMON_SOURCES = $(patsubst %,../common/%.cpp,$(COMMON_FILES))
COMMON_OBJECTS = $(patsubst %,obj/$(BASECXX)/common/%.o,$(COMMON_FILES))
COMMON_DEPS    = $(patsubst %,obj/$(BASECXX)/common/%.dep,$(COMMON_FILES))

HEX_DWBA_FILES   = dwba potential pwba ui wave_distort
HEX_DWBA_SOURCES = $(patsubst %,src/%.cpp,$(HEX_DWBA_FILES))
HEX_DWBA_OBJECTS = $(patsubst %,obj/$(BASECXX)/hex-dwba/%.o,$(HEX_DWBA_FILES))
HEX_DWBA_DEPS    = $(patsubst %,obj/$(BASECXX)/hex-dwba/%.dep,$(HEX_DWBA_FILES))

all: pollversion bin/hex-dwba

allclean: clean docclean

pollversion:
	touch ../common/hex-version.cpp

bin/hex-dwba: bin/hex-dwba.$(BASECXX)
	@ln -fs hex-dwba.$(BASECXX) bin/hex-dwba

bin/hex-dwba.$(BASECXX): $(COMMON_OBJECTS) $(HEX_DWBA_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $^ $(LFLAGS) -o $@

obj/$(BASECXX)/common/%.o: ../common/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(CXXFLAGS) -I../common $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

obj/$(BASECXX)/hex-dwba/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(CXXFLAGS) -I../common -Isrc $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

clean:
	rm -rf obj/$(BASECXX) bin/hex-dwba.$(BASECXX) bin/hex-dwba

doc:
	doxygen

docclean:
	rm -rf html

-include $(COMMON_DEPS) $(HEX_DWBA_DEPS)
