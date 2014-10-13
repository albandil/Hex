.PHONY: all clean

all: ../bin/hex-stack

../bin/hex-stack: stack.cpp
	g++ -std=c++11 -O2 $< -o $@

clean:
	rm -rf ../bin/hex-stack
