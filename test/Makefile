.PHONY: all clean

freebayes=../bin/freebayes
vcfuniq=../vcflib/bin/vcfuniq

all: test

test: $(freebayes) $(vcfuniq)
	prove -v t

$(freebayes):
	cd .. && $(MAKE)

$(vcfuniq):
	cd ../vcflib && make clean && make
