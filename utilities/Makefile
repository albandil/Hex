all:
	make -C stack
	make -C triangle
	mkdir -p bin; cd bin; pwd; ln -fs ../stack/stack ../triangle/triangle .
	
clean:
	make -C stack clean
	make -C triangle clean
	rm -rf bin
	