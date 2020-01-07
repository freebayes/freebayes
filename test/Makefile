.PHONY: all clean

freebayes=../bin/freebayes
vcfuniq=build/vcfuniq

all: test

test: $(freebayes) $(vcfuniq)
	prove -v t

$(freebayes):
	cd .. && $(MAKE)

$(vcfuniq):
	cmake -H. -Bbuild -S ../vcflib && cmake --build build -- -j4
