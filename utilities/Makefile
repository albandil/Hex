all:
	make -C far
	make -C stack
	make -C triangle
	make -C hdf2txt
	make -C hdf2hdf
	mkdir -p bin; cd bin; pwd; ln -fs ../stack/stack ../triangle/triangle ../far/far ../hdf2txt/hdf2txt ../hdf2hdf/hdf2hdf .
	
clean:
	make -C far clean
	make -C stack clean
	make -C triangle clean
	make -C hdf2txt clean
	make -C hdf2hdf clean
	rm -rf bin
	