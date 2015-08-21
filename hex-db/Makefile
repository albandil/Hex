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
GINACLIBS = -lcln -lginac

# ---------------------------------------------------------------------------- #
# SQLite3                                                                      #
# ---------------------------------------------------------------------------- #

SQLITEFLAGS =
SQLITELIBS = -lsqlite3

# ---------------------------------------------------------------------------- #
# SQLite++ (bundled)                                                           #
# ---------------------------------------------------------------------------- #

SQLITEPPFLAGS = -I$(shell pwd)/libs
SQLITEPPLIBS = -L$(shell pwd)/libs -lsqlitepp

# ---------------------------------------------------------------------------- #
# GSL                                                                          #
# ---------------------------------------------------------------------------- #

GSLFLAGS = $(shell gsl-config --cflags)
GSLLIBS = $(shell gsl-config --libs) 
                                                                               #
################################################################################


CFLAGS += $(SQLITEFLAGS) $(GINACFLAGS) $(SQLITEPPFLAGS) $(GSLFLAGS)
LFLAGS += $(SQLITELIBS) $(GINACLIBS) $(SQLITEPPLIBS) $(GSLLIBS) -Wl,-rpath,'$$ORIGIN/../lib'

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CXXRILIS = $(CFLAGS) -DNDEBUG -DNO_LAPACK -DNO_HDF -DNO_MPI -O2 -DGIT_COMMIT=\"$(GIT_COMMIT)\"
CXXDEBUG = $(CFLAGS) -DNO_HDF -DNO_LAPACK -DNO_MPI -g -O2 -DGIT_COMMIT=\"$(GIT_COMMIT)\"

CXXFLAGS = $(CXXDEBUG)
# CXXFLAGS = $(CXXRILIS)

SQLITEPP = libs/libsqlitepp.a
SQLITEPPNAMES = binders exception query session statement string transaction
SQLITEPPSRC = $(patsubst %,libs/sqlitepp/%.cpp,$(SQLITEPPNAMES))
SQLITEPPOBJ = $(patsubst %,obj/%.o,$(SQLITEPPNAMES))

FILES = arrays hex-db born misc special variables vec3d version
FILES_SOURCES = $(patsubst %,src/%.cpp,$(FILES))
FILES_OBJECTS = $(patsubst %,obj/%.o,$(FILES))
FILES_DEPENDENCIES = $(patsubst %,obj/%.dep,$(FILES))

VARS = asy ccs colls dcs xcs ics momtf scatamp scatamp-dir tcs tmat ionf ionamp tdcs stokes bornf tmatb
VARS_SOURCES = $(patsubst %,src/var/%.cpp,$(VARS))
VARS_OBJECTS = $(patsubst %,obj/%.o,$(VARS))
VARS_DEPENDENCIES = $(patsubst %,obj/%.dep,$(VARS))

HEXDBLIB = libhex-db.so
HEXDBEXE = hex-db

all: pollversion lib exe test

pollversion:
	@touch src/version.cpp

doc:
	doxygen

docclean:
	rm -rf html

dist:
	mkdir -p hex-db/libs hex-db/test
	cp --recursive --dereference README Makefile Doxyfile *.cbp ../media/hexe-small.gif src hex-db/
	cp --recursive libs/sqlitepp hex-db/libs
	cp test/test-hexdb-lib-c.c test/test-hexdb-lib-f.f95 test/Makefile hex-db/test
	sed -i "s/\(^GIT_COMMIT = \)\(.*\)/\1$(GIT_COMMIT)/g" hex-db/Makefile
	sed -i "s/\(^all:.*\)\(pollversion\)\(.*\)/\1\3/g" hex-db/Makefile
	sed -i "s/\.\.\/media\///g" hex-db/Doxyfile
	tar --create --gzip --file=hex-db-$(GIT_COMMIT).tar.gz hex-db
	rm -rf hex-db/

clean:
	+make -C test $(MAKECMDGOALS)
	rm -f obj/*.o obj/*.dep bin/hex-db* lib/*

distclean:
	rm -rf hex-db hex-db-*.tar.gz

allclean: clean distclean docclean

exe: bin/$(HEXDBEXE)

lib: lib/$(HEXDBLIB)

test: lib/$(HEXDBLIB)
	+make -C $@ $(MAKECMDGOALS)

bin/$(HEXDBEXE): obj/ui.o lib/$(HEXDBLIB).0
	@mkdir -p $(@D)
	$(CXX) $^ -L$(shell pwd)/lib/ -lhex-db $(LFLAGS) -o $@

lib/$(HEXDBLIB).0: lib/$(HEXDBLIB)

lib/$(HEXDBLIB): $(FILES_OBJECTS) $(VARS_OBJECTS) $(SQLITEPP)
	@mkdir -p $(@D)
	$(CXX) $^ $(LFLAGS) -fPIC -shared -Wl,-soname,$(HEXDBLIB).0 -o $@.0
	cd lib && ln -fs $(HEXDBLIB).0 $(HEXDBLIB) && cd ..

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

obj/%.o: src/var/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

# -----------------------------------------------------------------------------
# SQLite++ library

$(SQLITEPPOBJ) : obj/%.o : libs/sqlitepp/%.cpp
	$(CXX) $(SQLITEFLAGS) -Ilibs -O2 -fPIC -c $^ -o $@

$(SQLITEPP) : $(SQLITEPPOBJ)
	ar cr $@ $^

-include $(wildcard obj/*.dep)
