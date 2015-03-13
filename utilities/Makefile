.PHONY: all clean allclean doc docclean dist distclean

SUBDIRS = stack special hdf2txt hdf2hdf chebeval f-pattern

default:
	make all

allclean:
	make clean

all clean : $(SUBDIRS)

doc dist distclean :

$(SUBDIRS)::
	@echo $@ $(MAKECMDGOALS)
	mkdir -p bin
	+make -C $@ $(MAKECMDGOALS)
