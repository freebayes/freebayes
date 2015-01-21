all: vcflib/Makefile log
	cd src && $(MAKE)

log: src/version_git.h
	wget -q http://hypervolu.me/freebayes/build/$(shell cat src/version_git.h | grep v | cut -f 3 -d\  | sed s/\"//g) &

src/version_git.h:
	cd src && $(MAKE) autoversion
	touch src/version_git.h

# If 'vcflib/Makefile' doesn't exist,
# it means git submodules were not checked-out with 'git clone --recursive'.
# manually clone and update them now.
vcflib/Makefile:
	git submodule init
	git submodule update
	( cd vcflib ; git submodule init && git submodule update )


debug:
	cd src && $(MAKE) debug

install:
	cp bin/freebayes bin/bamleftalign /usr/local/bin/

uninstall:
	rm /usr/local/bin/freebayes /usr/local/bin/bamleftalign

test: all
	prove test/*.t

clean:
	cd src && $(MAKE) clean
	rm -f bin/*

.PHONY: all install uninstall clean
