all:
	cd src && make

install:
	cp bin/freebayes /usr/bin/

uninstall:
	rm /usr/bin/freebayes

clean:
	cd src && make clean
	rm -f bin/*
	rm -rf doc/*

.PHONY: all install uninstall clean
