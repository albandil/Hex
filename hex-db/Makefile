.PHONY: all exe lib clean allclean doc docclean dist distclean


################################################################################
# Editable section                                                             #

# ---------------------------------------------------------------------------- #
# Compiler and linker + flags                                                  #
# ---------------------------------------------------------------------------- #

CXX = g++

CFLAGS = -Wall -std=c++11 -pedantic-errors -fPIC
LFLAGS = -lstdc++

# ---------------------------------------------------------------------------- #
# GiNaC                                                                        #
# ---------------------------------------------------------------------------- #

GINACFLAGS =
GINACLIBS = -lginac -lcln

# ---------------------------------------------------------------------------- #
# SQLite3                                                                      #
# ---------------------------------------------------------------------------- #

SQLITEFLAGS =
SQLITELIBS = -lsqlite3

# ---------------------------------------------------------------------------- #
# GSL                                                                          #
# ---------------------------------------------------------------------------- #

GSLFLAGS = 
GSLLIBS = -lgsl -lgslcblas
                                                                               #
################################################################################


CFLAGS += $(SQLITEFLAGS) $(GINACFLAGS) $(GSLFLAGS)
LFLAGS += $(SQLITELIBS) $(GINACLIBS) $(GSLLIBS) -Wl,-rpath,'$$ORIGIN/../lib'

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CXXRILIS = $(CFLAGS) -DNDEBUG -DNO_LAPACK -O2 -DGIT_COMMIT=\"$(GIT_COMMIT)\"
CXXDEBUG = $(CFLAGS) -DNO_LAPACK -g -O2 -DGIT_COMMIT=\"$(GIT_COMMIT)\"

# CXXFLAGS = $(CXXDEBUG)
CXXFLAGS = $(CXXRILIS)

# compiler name
BASECXX = $(shell basename $(CXX))

# SQLite++ library
SQLITEPP_LIBRARY = obj/$(BASECXX)/libsqlitepp.a
SQLITEPP_FILES   = binders exception query session statement string transaction
SQLITEPP_SOURCES = $(patsubst %,../libs/sqlitepp/%.cpp,$(SQLITEPP_FILES))
SQLITEPP_OBJECTS = $(patsubst %,obj/$(BASECXX)/sqlitepp/%.o,$(SQLITEPP_FILES))
CFLAGS += -I../libs
LFLAGS += -Lobj/$(BASECXX) -lsqlitepp

# Common source files
COMMON_FILES   = hex-arrays hex-misc hex-special hex-vec3d hex-version
COMMON_SOURCES = $(patsubst %,../common/%.cpp,$(COMMON_FILES))
COMMON_OBJECTS = $(patsubst %,obj/$(BASECXX)/common/%.o,$(COMMON_FILES))
COMMON_DEPS    = $(patsubst %,obj/$(BASECXX)/common/%.dep,$(COMMON_FILES))

# Source files specific to hex-db.
HEX_DB_FILES   = db ui variables var/asy var/ccs var/colls var/dcs var/xcs var/ics var/momtf \
                 var/scatamp var/scatamp-dir var/tcs var/tmat var/ionf var/ionamp var/tdcs \
                 var/stokes var/spflip
HEX_DB_SOURCES = $(patsubst %,src/%.cpp,$(HEX_DB_FILES))
HEX_DB_OBJECTS = $(patsubst %,obj/$(BASECXX)/hex-db/%.o,$(HEX_DB_FILES))
HEX_DB_DEPS    = $(patsubst %,obj/$(BASECXX)/hex-db/%.dep,$(HEX_DB_FILES))

HEX_DB_NAME = hex-db-$(BASECXX)
HEX_DB_LIB  = lib$(HEX_DB_NAME).so
HEX_DB_EXE  = hex-db

all: pollversion lib exe

pollversion:
	@touch ../common/hex-version.cpp

doc:
	doxygen

docclean:
	rm -rf html

clean:
	rm -rf obj/$(BASECXX)
	rm -f bin/$(HEX_DB_EXE) bin/$(HEX_DB_EXE).$(BASECXX) lib/$(HEX_DB_LIB)*

allclean: clean docclean

exe: bin/$(HEX_DB_EXE)

lib: lib/$(HEX_DB_LIB)

bin/$(HEX_DB_EXE): bin/$(HEX_DB_EXE).$(BASECXX)
	@ln -fs $(HEX_DB_EXE).$(BASECXX) bin/$(HEX_DB_EXE)

bin/$(HEX_DB_EXE).$(BASECXX): obj/$(BASECXX)/hex-db/ui.o lib/$(HEX_DB_LIB).0
	@mkdir -p $(@D)
	$(CXX) $^ -Llib -l$(HEX_DB_NAME) $(LFLAGS) -o $@

lib/$(HEX_DB_LIB).0: lib/$(HEX_DB_LIB)

lib/$(HEX_DB_LIB): $(COMMON_OBJECTS) $(HEX_DB_OBJECTS) $(SQLITEPP_LIBRARY)
	@mkdir -p $(@D)
	$(CXX) $^ $(LFLAGS) -fPIC -shared -Wl,-soname,$(HEX_DB_LIB).0 -o $@.0
	cd lib && ln -fs $(HEX_DB_LIB).0 $(HEX_DB_LIB) && cd ..

obj/$(BASECXX)/common/%.o: ../common/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -I../common -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

obj/$(BASECXX)/hex-db/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -Isrc -I../common -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

# -----------------------------------------------------------------------------
# SQLite++ library

$(SQLITEPP_OBJECTS) : obj/$(BASECXX)/sqlitepp/%.o : ../libs/sqlitepp/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(SQLITEFLAGS) -I../libs -O2 -fPIC -c $^ -o $@

$(SQLITEPP_LIBRARY) : $(SQLITEPP_OBJECTS)
	@mkdir -p $(@D)
	ar cr $@ $^

-include $(COMMON_DEPS) $(HEX_DB_DEPS)
