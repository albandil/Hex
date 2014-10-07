all: test-hexdb-lib-c test-hexdb-lib-f
clean:
	rm -f *.o test-hexdb-lib-c test-hexdb-lib-f

#
# C test
#

test-hexdb-lib-c: test-hexdb-lib-c.o
	gcc test-hexdb-lib-c.o -L../lib -lhex-db -o test-hexdb-lib-c

test-hexdb-lib-c.o: test-hexdb-lib-c.c
	gcc -c test-hexdb-lib-c.c -o test-hexdb-lib-c.o

#
# Fortran test
#

test-hexdb-lib-f: test-hexdb-lib-f.o
	gfortran test-hexdb-lib-f.o -L../lib -lhex-db -o test-hexdb-lib-f

test-hexdb-lib-f.o: test-hexdb-lib-f.f95
	gfortran -c test-hexdb-lib-f.f95 -o test-hexdb-lib-f.o
