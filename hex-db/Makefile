.PHONY: all exe lib clean doc

CPP = g++

CFLAGS = -Wall -std=c++11 -I $(shell pwd)/libs -fPIC
LFLAGS = -L $(shell pwd)/libs -lsqlitepp -lsqlite3 -lgsl -lgslcblas -lo2scl -lfftw3

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

all: pollversion exe lib

pollversion:
	@touch src/version.cpp

doc:
	doxygen

clean:
	rm -f obj/*.o obj/*.dep bin/hex-db*

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
	g++ -O2 -fPIC -c $^ -o $@

$(SQLITEPP) : $(SQLITEPPOBJ)
	ar cr $@ $^
