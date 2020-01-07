# convenience wrapper for make

all: log freebayes vcflib/CMakeLists.txt

freebayes:
	cmake -H. -Bbuild && cmake --build build -- $(MAKEFLAGS)

log: src/version_git.h
	wget -q http://hypervolu.me/freebayes/build/$(shell cat src/version_git.h | grep v | cut -f 3 -d\  | sed s/\"//g) &

vcflib/CMakeLists.txt:
SeqLib/configure:
	@echo "To build freebayes you must use git to also download its submodules."
	@echo "Do so by downloading freebayes again using this command (note --recursive flag):"
	@echo "    git clone --recursive git://github.com/ekg/freebayes.git"
	@error

test:
	cd test && $(MAKE) test

clean:
	rm -rf test/build
	rm -rf build
	rm -fr bin/*

.PHONY: all clean test freebayes
