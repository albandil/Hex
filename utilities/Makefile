all:
	make -C far
	make -C stack
	make -C triangle
	make -C hdf2txt
	mkdir -p bin; cd bin; pwd; ln -fs ../stack/stack ../triangle/triangle ../far/far ../hdf2txt/hdf2txt .
	
clean:
	make -C far clean
	make -C stack clean
	make -C triangle clean
	make -C hdf2txt clean
	rm -rf bin
	