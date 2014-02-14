all: vcflib/Makefile
	cd src && $(MAKE)

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

clean:
	cd src && $(MAKE) clean
	rm -f bin/*

.PHONY: all install uninstall clean
