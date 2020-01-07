.PHONY: all clean

freebayes=../bin/freebayes
vcfuniq=build/vcfuniq

all: test

test: $(freebayes) $(vcfuniq)
	prove -v t

$(freebayes):
	cd .. && $(MAKE)

build:
	mkdir -p build

$(vcfuniq): build
	cmake -H. -Bbuild -S ../vcflib && cmake --build build -- -j4

clean:
	rm -rf build
