.PHONY: all exe lib clean allclean doc docclean dist distclean

##############################################################
# Editable part

CPP = g++
# CPP = icpc

##############################################################

CFLAGS = -Wall -std=c++11 -I $(shell pwd)/libs -fPIC
LFLAGS = -L $(shell pwd)/libs -lsqlitepp -lsqlite3 -lgsl -lgslcblas -lfftw3

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

CXXRILIS = $(CFLAGS) -DNDEBUG -DNO_HDF -O2 -DGIT_COMMIT=\"$(GIT_COMMIT)\"
CXXDEBUG = $(CFLAGS) -DNO_HDF -g -O2 -DGIT_COMMIT=\"$(GIT_COMMIT)\"

# CXXFLAGS = $(CXXDEBUG)
CXXFLAGS = $(CXXRILIS)

SQLITEPP = libs/libsqlitepp.a
SQLITEPPNAMES = binders exception query session statement string transaction
SQLITEPPSRC = $(patsubst %,libs/sqlitepp/%.cpp,$(SQLITEPPNAMES))
SQLITEPPOBJ = $(patsubst %,obj/%.o,$(SQLITEPPNAMES))

FILES = arrays hex-db specf variables vec3d version
FILES_SOURCES = $(patsubst %,src/%.cpp,$(FILES))
FILES_OBJECTS = $(patsubst %,obj/%.o,$(FILES))
FILES_DEPENDENCIES = $(patsubst %,obj/%.dep,$(FILES))

VARS = ccs colls dcs xcs ics momtf scatamp tcs tmat ionf ionamp tdcs stokes
VARS_SOURCES = $(patsubst %,src/var/%.cpp,$(VARS))
VARS_OBJECTS = $(patsubst %,obj/%.o,$(VARS))
VARS_DEPENDENCIES = $(patsubst %,obj/%.dep,$(VARS))

all: pollversion exe

pollversion:
	@touch src/version.cpp

doc:
	doxygen

docclean:
	rm -rf html

dist:
	mkdir -p hex-db/libs
	cp --recursive --dereference README Makefile Doxyfile *.cbp ../media/hexe-small.gif src hex-db/
	cp --recursive libs/sqlitepp hex-db/libs
	sed -i "s/\(^GIT_COMMIT = \)\(.*\)/\1$(GIT_COMMIT)/g" hex-db/Makefile
	sed -i "s/\(^all:.*\)\(pollversion\)\(.*\)/\1\3/g" hex-db/Makefile
	sed -i "s/\.\.\/media\///g" hex-db/Doxyfile
	tar --create --gzip --file=hex-db-$(GIT_COMMIT).tar.gz hex-db
	rm -rf hex-db/

clean:
	rm -f obj/*.o obj/*.dep bin/hex-db*

distclean:
	rm -rf hex-db hex-db-*.tar.gz

allclean: clean distclean docclean

exe: bin/hex-db

lib: bin/hex-db.so

bin/hex-db: obj/ui.o $(FILES_OBJECTS) $(VARS_OBJECTS) $(SQLITEPP)
	@mkdir -p $(@D)
	$(CPP) $^ $(LFLAGS) -o $@

bin/hex-db.so: $(FILES_OBJECTS) $(VARS_OBJECTS) $(SQLITEPP)
	@mkdir -p $(@D)
	$(CPP) $^ $(LFLAGS) -fPIC -shared -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CPP) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

obj/%.o: src/var/%.cpp
	@mkdir -p $(@D)
	$(CPP) $(CXXFLAGS) -c $< -o $@ -MMD -MF $(patsubst %.o,%.dep,$@)

# -----------------------------------------------------------------------------
# SQLite++ library

$(SQLITEPPOBJ) : obj/%.o : libs/sqlitepp/%.cpp
	$(CPP) -O2 -fPIC -c $^ -o $@

$(SQLITEPP) : $(SQLITEPPOBJ)
	ar cr $@ $^

-include $(wildcard obj/*.dep)
