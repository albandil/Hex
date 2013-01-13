CPP = g++

CFLAGS = -Wall -std=c++0x -I ./include
LFLAGS = -L ./lib -lsqlitepp -lsqlite3 -lgsl -lgslcblas -lhdf5 -lhdf5_cpp

CXXRILIS = $(CFLAGS) -DNDEBUG -O2
CXXDEBUG = $(CFLAGS) -g -O2

all: bin/hex-db bin/hex-db_debug

clean:
	rm -f obj/*.o bin/hex-db*

# -----------------------------------------------------------------------------
# Main targets

bin/hex-db: obj/ui.o obj/hex-db.o obj/specf.o obj/angs.o obj/arrays.o
	$(CPP) $^ $(LFLAGS) -o $@
bin/hex-db_debug: obj/ui_debug.o obj/hex-db_debug.o obj/specf_debug.o obj/angs_debug.o obj/arrays_debug.o
	$(CPP) $^ $(LFLAGS) -o $@

# -----------------------------------------------------------------------------
# Release

obj/ui.o: src/ui.cpp src/hex-db.h
	$(CPP) $(CXXRILIS) -c $< -o $@
obj/hex-db.o: src/hex-db.cpp src/hex-db.h src/angs.h src/specf.h src/complex.h src/arrays.h src/interpolate.h
	$(CPP) $(CXXRILIS) -c $< -o $@
obj/specf.o: src/specf.cpp src/complex.h
	$(CPP) $(CXXRILIS) -c $< -o $@
obj/angs.o: src/angs.cpp
	$(CPP) $(CXXRILIS) -c $< -o $@
obj/arrays.o: src/arrays.cpp src/arrays.h src/complex.h
	$(CPP) $(CXXRILIS) -c $< -o $@

# -----------------------------------------------------------------------------
# Debug

obj/ui_debug.o: src/ui.cpp src/hex-db.h
	$(CPP) $(CXXDEBUG) -c $< -o $@
obj/hex-db_debug.o: src/hex-db.cpp src/hex-db.h src/angs.h src/specf.h src/complex.h src/arrays.h src/interpolate.h
	$(CPP) $(CXXDEBUG) -c $< -o $@
obj/specf_debug.o: src/specf.cpp src/complex.h
	$(CPP) $(CXXDEBUG) -c $< -o $@
obj/angs_debug.o: src/angs.cpp
	$(CPP) $(CXXDEBUG) -c $< -o $@
obj/arrays_debug.o: src/arrays.cpp src/arrays.h src/complex.h
	$(CPP) $(CXXDEBUG) -c $< -o $@