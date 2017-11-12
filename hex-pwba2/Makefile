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

CXXRILIS = $(CFLAGS) -O3 -fcx-limited-range -DHEX_GIT_COMMIT=\"$(GIT_COMMIT)\" -DNDEBUG
CXXDEBUG = $(CFLAGS) -g -fcx-limited-range -DHEX_GIT_COMMIT=\"$(GIT_COMMIT)\"

# CXXFLAGS = $(CXXDEBUG)
CXXFLAGS = -g $(CXXRILIS)

BASECXX = $(shell basename $(CXX))

COMMON_FILES   = hex-arrays hex-born hex-hydrogen hex-misc hex-special hex-spgrid \
                 hex-vec3d hex-version
COMMON_SOURCES = $(patsubst %,../common/%.cpp,$(COMMON_FILES))
COMMON_OBJECTS = $(patsubst %,obj/$(BASECXX)/common/%.o,$(COMMON_FILES))
COMMON_DEPS    = $(patsubst %,obj/$(BASECXX)/common/%.dep,$(COMMON_FILES))

HEX_PWBA2_FILES   = pwba2-pwx pwba2-spg radial ui hyp_2F1
HEX_PWBA2_SOURCES = $(patsubst %,src/%.cpp,$(HEX_PWBA2_FILES))
HEX_PWBA2_OBJECTS = $(patsubst %,obj/$(BASECXX)/hex-pwba2/%.o,$(HEX_PWBA2_FILES))
HEX_PWBA2_DEPS    = $(patsubst %,obj/$(BASECXX)/hex-pwba2/%.dep,$(HEX_PWBA2_FILES))

all: pollversion bin/hex-pwba2

allclean: clean docclean

pollversion:
	@touch ../common/hex-version.cpp

bin/hex-pwba2: bin/hex-pwba2.$(BASECXX)
	@ln -fs hex-pwba2.$(BASECXX) bin/hex-pwba2

bin/hex-pwba2.$(BASECXX): $(COMMON_OBJECTS) $(HEX_PWBA2_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $^ $(LFLAGS) -o $@

obj/$(BASECXX)/common/%.o: ../common/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -I../common -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

obj/$(BASECXX)/hex-pwba2/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -I../common -Isrc -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

clean:
	rm -f $(OBJECTS) $(DEPENDENCIES) bin/pwba2

doc:
	doxygen

docclean:
	rm -rf html

-include $(COMMON_DEPS) $(HEX_PWBA2_DEPS)
