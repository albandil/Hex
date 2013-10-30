SUBDIRS = hex-db hex-ecs hex-pwba hex-dwba utilities

.PHONY: $(SUBDIRS)

default : all

$(SUBDIRS)::
	+make -C $@ $(MAKECMDGOALS)
	mkdir -p bin; cd bin; ln -fs ../$@/bin/* .; cd ..
	
all clean doc : $(SUBDIRS)

