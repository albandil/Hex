all:
	make -C hex-db
	make -C hex-ecs
	make -C hex-dwba
	make -C hex-pwba
	make -C utilities
	mkdir -p bin; cd bin; ln -fs ../hex-db/bin/* ../hex-ecs/bin/* ../hex-dwba/bin/* ../hex-pwba/bin/* ../utilities/bin/* .
	
clean:
	make -C hex-db clean
	make -C hex-ecs clean
	make -C hex-dwba clean
	make -C hex-pwba clean
	make -C utilities clean
	rm -rf bin
