all:
	cd src && make
	cd src && make install

.PHONY: all

doc:
	doxygen doxygen.cfg

.PHONY: doc

clean:
	cd src && make clean
	rm -f bin/*
	rm -rf doc/*

.PHONY: clean
