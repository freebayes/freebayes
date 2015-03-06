.PHONY: all clean

freebayes=../bin/freebayes

all: test

test: $(freebayes)
	prove -v t

$(freebayes):
	cd .. && $(MAKE)
