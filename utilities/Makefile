all:
	make -C stack
	make -C triangle
	make -C hdf2txt
	make -C hdf2hdf
	make -C chebeval
	mkdir -p bin; cd bin; pwd; ln -fs ../stack/stack ../triangle/triangle ../hdf2txt/hdf2txt ../hdf2hdf/hdf2hdf ../chebeval/chebeval .
	
clean:
	make -C stack clean
	make -C triangle clean
	make -C hdf2txt clean
	make -C hdf2hdf clean
	make -C chebeval clean
	rm -rf bin
	