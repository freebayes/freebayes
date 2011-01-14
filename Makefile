all:
	cd src && $(MAKE)

debug:
	cd src && $(MAKE) debug

install:
	cp bin/freebayes bin/bamleftalign bin/bamfiltertech /usr/local/bin/

uninstall:
	rm /usr/local/bin/freebayes /usr/local/bin/bamleftalign /usr/local/bin/bamfiltertech

clean:
	cd src && $(MAKE) clean
	rm -f bin/*

.PHONY: all install uninstall clean
