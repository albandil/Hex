.PHONY: all clean doc

SUBDIRS = stack f hdf2txt hdf2hdf chebeval

default : all

$(SUBDIRS)::
	+make -C $@ $(MAKECMDGOALS)
	mkdir -p bin; cd bin; ln -fs ../$@/$@ .
	
all clean: $(SUBDIRS)

doc:
	