SUBDIRS = hex-db hex-ecs hex-pwba2 hex-dwba hex-recs utilities

GIT_COMMIT = $(shell git rev-parse HEAD | cut -c -8)

.PHONY: $(SUBDIRS)

default : all

$(SUBDIRS)::
	+make -C $@ $(MAKECMDGOALS)
	mkdir -p bin; cd bin; ln -fs ../$@/bin/* .; cd ..
	
all clean allclean doc docclean distclean : $(SUBDIRS)

dist: $(SUBDIRS)
	@mkdir -p release
	cp hex-db/hex-db-$(GIT_COMMIT).tar.gz \
	   hex-dwba/hex-dwba-$(GIT_COMMIT).tar.gz \
	   hex-ecs/hex-ecs-$(GIT_COMMIT).tar.gz \
	   hex-pwba2/hex-pwba2-$(GIT_COMMIT).tar.gz     release/

