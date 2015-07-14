all: test-hexdb-lib-c test-hexdb-lib-f
allclean: clean

clean:
	rm -f *.o test-hexdb-lib-c test-hexdb-lib-f

doc: ;
docclean: ;
dist: ;
distclean: ;

#
# C test
#

test-hexdb-lib-c: test-hexdb-lib-c.o
	gcc -std=c99 test-hexdb-lib-c.o -L../lib -lhex-db -o test-hexdb-lib-c

test-hexdb-lib-c.o: test-hexdb-lib-c.c
	gcc -std=c99 -c test-hexdb-lib-c.c -o test-hexdb-lib-c.o

#
# Fortran test
#

test-hexdb-lib-f: test-hexdb-lib-f.o
	gfortran test-hexdb-lib-f.o -L../lib -lhex-db -o test-hexdb-lib-f

test-hexdb-lib-f.o: test-hexdb-lib-f.f95
	gfortran -c test-hexdb-lib-f.f95 -o test-hexdb-lib-f.o
