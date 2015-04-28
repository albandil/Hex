.PHONY: all clean

CXXFLAGS = -Wall -std=c++11 -DNO_MPI -fopenmp
LDFLAGS = -fopenmp

TESTS = test-spgrid
OBJS = $(patsubst %, %.o, $(TESTS))

all: $(TESTS)

clean:
	rm -f $(TESTS) $(OBJS)
