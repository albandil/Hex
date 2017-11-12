.PHONY: all doc docclean clean allclean pollversion dist distclean

################################################################################
# Editable section                                                             #

# ---------------------------------------------------------------------------- #
# Compiler and linker + flags                                                  #
#  - Hex-pecs requires compiler with C++11 support (GCC >= 5, ICC >= 2016)     #
# ---------------------------------------------------------------------------- #

# C++ compiler (g++, clang, ...)
CXX = g++

# build mode ("debug" or "release")
TARGET = release

# global compiler flags
CFLAGS = -fopenmp -Wall -mtune=native -march=native -fcx-limited-range
LFLAGS = -fopenmp -lstdc++

# ---------------------------------------------------------------------------- #
# BLAS & LAPACK                                                                #
#                                                                              #
# The macro DISABLE_PARALLEL_PRECONDITION will inhibit nested parallelism      #
# e.g. when using threaded OpenBLAS. This may (and may not) deliver            #
# performance benefits in exchange of giving up possibility to precondition    #
# several blocks concurrently.                                                 #
# ---------------------------------------------------------------------------- #

BLASLAPACKFLAGS = 
BLASLAPACKLIBS = -llapack -lblas

# BLASLAPACKFLAGS = -DDISABLE_PARALLEL_PRECONDITION
# BLASLAPACKLIBS = -lopenblas

# ---------------------------------------------------------------------------- #
# GSL                                                                          #
# ---------------------------------------------------------------------------- #

GSLFLAGS = $(shell gsl-config --cflags)
GSLLIBS = $(shell gsl-config --libs)

# ---------------------------------------------------------------------------- #
# OpenMPI environment (optional)                                               #
# ---------------------------------------------------------------------------- #

OPENMPIFLAGS = -DWITH_MPI $(shell mpic++ --showme:compile)
OPENMPILIBS  = $(shell mpic++ --showme:link)

# ---------------------------------------------------------------------------- #
# OpenCL (optional)                                                            #
# ---------------------------------------------------------------------------- #

# OPENCLFLAGS = -DWITH_OPENCL
# OPENCLLIBS = -lOpenCL

# ---------------------------------------------------------------------------- #
# PNG++ (optional)                                                             #
# ---------------------------------------------------------------------------- #

# PNGPPFLAGS = -DWITH_PNG -I../libs
# PNGPPLIBS = -lpng
                                                                               #
################################################################################


# compiler flags
CFLAGS += -fPIC -std=c++11 -pedantic-errors -DHEX_GIT_COMMIT=\"$(GIT_COMMIT)\" \
          -Isrc -I../common \
          $(SUITESPARSEFLAGS) $(SUPERLUFLAGS) $(SUPERLUDISTFLAGS) $(OPENMPIFLAGS) \
          $(OPENCLFLAGS) $(GSLFLAGS) $(BLASLAPACKFLAGS) $(PNGPPFLAGS) 

# linker flags
LFLAGS += -fPIC -rdynamic \
          $(SUITESPARSELIBS) $(SUPERLULIBS) $(SUPERLUDISTLIBS) $(OPENMPILIBS) \
          $(OPENCLLIBS) $(GSLLIBS) $(BLASLAPACKLIBS) $(PNGPPLIBS)

# git revision number
GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

# DEBUG and RELEASE variants
ifeq ($(TARGET),debug)
CXXFLAGS = -g $(CFLAGS)
else ifeq ($(TARGET),release)
CXXFLAGS = -g $(CFLAGS) -O3 -DNDEBUG
else
$(error "Unknown target $(TARGET)")
endif

# compiler name
BASECXX = $(shell basename $(CXX))

# common source files
COMMON_FILES   = hex-arrays hex-hdffile hex-hydrogen \
                 hex-misc hex-special hex-version
COMMON_SOURCES = $(patsubst %,../common/%.cpp,$(COMMON_FILES))
COMMON_OBJECTS = $(patsubst %,obj/$(BASECXX)/common/%.o,$(COMMON_FILES))
COMMON_DEPS    = $(patsubst %,obj/$(BASECXX)/common/%.dep,$(COMMON_FILES))

# source files specific to hex-pecs
HEX_PECS_FILES   = amplitudes ang io main matops numerov os radial
HEX_PECS_SOURCES = $(patsubst %,src/%.cpp,$(HEX_ECS_FILES))
HEX_PECS_OBJECTS = $(patsubst %,obj/$(BASECXX)/hex-pecs/%.o,$(HEX_PECS_FILES))
HEX_PECS_DEPS    = $(patsubst %,obj/$(BASECXX)/hex-pecs/%.dep,$(HEX_PECS_FILES))

# OpenCL kernels
OCLSOURCES = 
OCLKERNELS = 

# build all targets
all: pollversion bin/hex-pecs.$(BASECXX)

# dummy target to force update of git revision
pollversion:
	@touch ../common/hex-version.cpp

# build documentation from source comments using Doxygen
doc:
	doxygen

# remove documentation
docclean:
	rm -rf html

# remove build files
clean:
	rm -f bin/hex-pecs bin/hex-pecs.$(BASECXX) $(COMMON_OBJECTS) $(COMMON_DEPS) $(HEX_PECS_OBJECTS) $(HEX_PECS_DEPS)

# remove everything (build files, documentation and distribution files)
allclean: clean docclean

# main program target
bin/hex-pecs.$(BASECXX): $(OCLKERNELS) $(COMMON_OBJECTS) $(HEX_PECS_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(HEX_PECS_OBJECTS) $(COMMON_OBJECTS) $(LFLAGS) -o $@
	@ln -fs hex-pecs.$(BASECXX) bin/hex-pecs

# CL -> INC rule (convert OpenCL sources to C byte array)
%.inc: %.cl
	xxd -i $< > $@

# CPP -> O rule (compile source files)
obj/$(BASECXX)/common/%.o: ../common/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)
obj/$(BASECXX)/hex-pecs/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

# include generated source dependencied
-include $(COMMON_DEPS) $(HEX_PECS_DEPS)
