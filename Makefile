all: vcflib/Makefile log
	cd src && $(MAKE)

log: src/version_git.h
	wget -q http://hypervolu.me/freebayes/build/$(shell cat src/version_git.h | grep v | cut -f 3 -d\  | sed s/\"//g) &

src/version_git.h:
	cd src && $(MAKE) autoversion
	touch src/version_git.h

vcflib/Makefile:
	@echo "To build freebayes you must use git to also download its submodules."
	@echo "Do so by downloading freebayes again using this command (note --recursive flag):"
	@echo "    git clone --recursive git://github.com/ekg/freebayes.git"
	@error

debug:
	cd src && $(MAKE) debug

install:
	cp bin/freebayes bin/bamleftalign /usr/local/bin/

uninstall:
	rm /usr/local/bin/freebayes /usr/local/bin/bamleftalign

test:
	cd test && make test

clean:
	cd src && $(MAKE) clean
	rm -fr bin/*

.PHONY: all install uninstall clean test
