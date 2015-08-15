.PHONY: all clean

all: ../bin/hex-stack

../bin/hex-stack: stack.cpp
	$(CXX) -std=c++11 -O2 $< -o $@ -lstdc++

clean:
	rm -rf ../bin/hex-stack
