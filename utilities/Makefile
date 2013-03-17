all:
	make -C far
	make -C stack
	make -C triangle
	mkdir -p bin; cd bin; pwd; ln -fs ../stack/stack ../triangle/triangle ../far/far .
	
clean:
	make -C far clean
	make -C stack clean
	make -C triangle clean
	rm -rf bin
	