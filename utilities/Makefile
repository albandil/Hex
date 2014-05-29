.PHONY: all clean allclean doc docclean dist distclean

SUBDIRS = stack special hdf2txt hdf2hdf chebeval

default:
	make all

allclean:
	make clean

all clean : $(SUBDIRS)

doc dist distclean :

$(SUBDIRS)::
	@echo $@ $(MAKECMDGOALS)
	+make -C $@ $(MAKECMDGOALS)
	mkdir -p bin; cd bin; ln -fs ../$@/$@ .
