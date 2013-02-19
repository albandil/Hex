CPP = g++

CFLAGS = -Wall -std=c++11 -I $(shell pwd)/include
LFLAGS = -L $(shell pwd)/lib -lsqlitepp -lsqlite3 -lgsl -lgslcblas

CXXRILIS = $(CFLAGS) -DNDEBUG -DNO_HDF -O2
CXXDEBUG = $(CFLAGS) -DNO_HDF -g -O2

SQLITEPP = lib/libsqlitepp.a
SQLITEPPNAMES = binders exception query session statement string transaction
SQLITEPPSRC = $(patsubst %,lib/sqlitepp/%.cpp,$(SQLITEPPNAMES))
SQLITEPPOBJ = $(patsubst %,obj/%.o,$(SQLITEPPNAMES))

all: bin/hex-db

clean:
	rm -f obj/*.o bin/hex-db*

VARS = ccs colls dcs xcs ics momtf scatamp tcs tmat
VARS_SOURCES = $(patsubst %,src/var/%.cpp,$(VARS))
VARS_OBJECTS = $(patsubst %,obj/%.o,$(VARS))

# -----------------------------------------------------------------------------
# Main targets

bin/hex-db: obj/ui.o obj/hex-db.o obj/specf.o obj/angs.o obj/arrays.o obj/variables.o $(VARS_OBJECTS) $(SQLITEPP)
	@mkdir -p $(@D)
	$(CPP) $^ $(LFLAGS) -o $@

# -----------------------------------------------------------------------------
# Release

$(VARS_OBJECTS):obj/%.o: src/var/%.cpp src/variables.h src/interpolate.h src/arrays.h src/complex.h src/specf.h
	@mkdir -p $(@D)
	$(CPP) $(CXXRILIS) -c $< -o $@

obj/ui.o: src/ui.cpp src/hex-db.h src/variables.h
	@mkdir -p $(@D)
	$(CPP) $(CXXRILIS) -c $< -o $@
	
obj/hex-db.o: src/hex-db.cpp src/hex-db.h src/angs.h src/specf.h src/complex.h src/arrays.h src/interpolate.h src/variables.h
	@mkdir -p $(@D)
	$(CPP) $(CXXRILIS) -c $< -o $@
	
obj/specf.o: src/specf.cpp src/complex.h
	@mkdir -p $(@D)
	$(CPP) $(CXXRILIS) -c $< -o $@
	
obj/angs.o: src/angs.cpp
	@mkdir -p $(@D)
	$(CPP) $(CXXRILIS) -c $< -o $@
	
obj/arrays.o: src/arrays.cpp src/arrays.h src/complex.h
	@mkdir -p $(@D)
	$(CPP) $(CXXRILIS) -c $< -o $@
	
obj/variables.o: src/variables.cpp src/variables.h src/arrays.h src/complex.h
	@mkdir -p $(@D)
	$(CPP) $(CXXRILIS) -c $< -o $@

# -----------------------------------------------------------------------------
# SQLite++ library

$(SQLITEPPOBJ) : obj/%.o : lib/sqlitepp/%.cpp
	g++ -O2 -c $^ -o $@

$(SQLITEPP) : $(SQLITEPPOBJ)
	ar cr $@ $^
