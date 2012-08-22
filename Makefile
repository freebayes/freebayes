DESTDIR=/usr/local/bin

all:
	cd src && $(MAKE)

debug:
	cd src && $(MAKE) debug

install:
	install -m 0755 bin/freebayes $(DESTDIR)
	install -m 0755 bin/bamleftalign $(DESTDIR)

uninstall:
	rm /usr/local/bin/freebayes /usr/local/bin/bamleftalign

clean:
	cd src && $(MAKE) clean
	rm -f bin/*

.PHONY: all install uninstall clean
