all:
	cd src && $(MAKE)

install:
	cp bin/freebayes /usr/bin/

uninstall:
	rm /usr/bin/freebayes

clean:
	cd src && $(MAKE) clean
	rm -f bin/*

.PHONY: all install uninstall clean
