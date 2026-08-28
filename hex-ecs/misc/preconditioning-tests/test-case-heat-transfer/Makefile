# The test cases are built against the sources of the surrounding Hex tree and
# against UMFPACK. Override these on the command line if they live elsewhere:
#     make HEX_PATH=~/Hex SUITESPARSE_INCLUDE=/opt/suitesparse-4.5.6/include
HEX_PATH ?= ../../../..
SUITESPARSE_INCLUDE ?= /usr/include/suitesparse

INCLUDES ?= -I$(HEX_PATH)/common \
            -I$(HEX_PATH)/libs \
            -I$(HEX_PATH)/hex-ecs/src \
            -I$(SUITESPARSE_INCLUDE)

LIBS ?= -lumfpack -lopenblas -lgsl -lgslcblas -lpng

CXX ?= g++

CXXFLAGS ?= -std=c++17 -g -O2 $(INCLUDES) -D_LONGINT -DWITH_UMFPACK -DWITH_PNG -DNDEBUG

LFLAGS ?= $(LIBS)

COMMON_DEPS_FILES = hex-arrays hex-blas hex-hdffile hex-misc hex-matrix hex-special hex-vtkfile
COMMON_DEPS_SRCS  = $(patsubst %,$(HEX_PATH)/common/%.cpp,$(COMMON_DEPS_FILES))
COMMON_DEPS_DEPS  = $(patsubst %,obj/common/%.dep,$(COMMON_DEPS_FILES))
COMMON_DEPS_OBJS  = $(patsubst %,obj/common/%.o,$(COMMON_DEPS_FILES))

ECS_DEPS_FILES = bspline gauss luft radial slater factorizers/lu-umfpack
ECS_DEPS_SRCS  = $(patsubst %,$(HEX_PATH)/hex-ecs/src/%.cpp,$(ECS_DEPS_FILES))
ECS_DEPS_DEPS  = $(patsubst %,obj/hex-ecs/%.dep,$(ECS_DEPS_FILES))
ECS_DEPS_OBJS  = $(patsubst %,obj/hex-ecs/%.o,$(ECS_DEPS_FILES))

HT_EXE = bin/ht
HT_SRC = src/ht.cpp
HT_OBJ = obj/ht.o
HT_DEP = obj/ht.dep

all: $(HT_EXE)

clean:
	rm -rf $(COMMON_DEPS_OBJS) $(COMMON_DEPS_DEPS) $(ECS_DEPS_OBJS) $(ECS_DEPS_DEPS) $(HT_OBJ) $(HT_DEP)

$(HT_EXE): $(COMMON_DEPS_OBJS) $(ECS_DEPS_OBJS) $(HT_OBJ)
	@mkdir -p $(@D)
	$(CXX) $(LFLAGS) $^ -o $@

obj/common/%.o: $(HEX_PATH)/common/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

obj/hex-ecs/%.o: $(HEX_PATH)/hex-ecs/src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

$(HT_OBJ): $(HT_SRC)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

-include $(COMMON_DEPS_DEPS) $(ECS_DEPS_DEPS) $(HT_DEP)
