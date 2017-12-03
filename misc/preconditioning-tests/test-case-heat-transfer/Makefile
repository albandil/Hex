HEX_PATH ?= /home/jacob/Dokumenty/prog/Hex

INCLUDES ?= -I$(HEX_PATH)/common \
            -I$(HEX_PATH)/libs \
            -I$(HEX_PATH)/hex-ecs/src \
            -I/opt/suitesparse-4.5.6/include

LIBS ?= -L/opt/suitesparse-4.5.6/lib -lumfpack -Wl,-rpath,/opt/suitesparse-4.5.6/lib \
        -L/opt/openblas-0.2.20/lib -lopenblas_openmp -Wl,-rpath,/opt/openblas-0.2.20/lib \
        -lgsl -lgslcblas -lpng

CXX ?= g++

CXXFLAGS ?= -g -O2 $(INCLUDES) -DWITH_UMFPACK -DWITH_PNG -DNDEBUG

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
