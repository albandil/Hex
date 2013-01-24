all:
	make -C hex-db
	make -C hex-main
	make -C hex-dwba
	make -C hex-pwba
	cd bin; ln -fs ../hex-db/bin/* ../hex-main/bin/* ../hex-dwba/bin/* ../hex-pwba/bin/* .
	
	
clean:
	make -C hex-db clean
	make -C hex-main clean
	make -C hex-dwba clean
	make -C hex-pwba clean
	rm -f bin/*
