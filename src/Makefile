################################################################################
# Makefile for FreeBayes
# Erik Garrison, 2010
# Boston College
################################################################################

# Compiler
CXX=g++
C=gcc

# Compiler flags
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -g
#CFLAGS=-O3 -static -D VERBOSE_DEBUG  # enables verbose debugging via --debug2

BAMTOOLS_ROOT=../bamtools
SEQLIB_ROOT=../SeqLib
VCFLIB_ROOT=../vcflib
TABIX_ROOT=$(VCFLIB_ROOT)/tabixpp
HTSLIB_ROOT=$(TABIX_ROOT)/htslib

LIBS = -lz -lm -lpthread
INCLUDE = -I../ttmath -I$(BAMTOOLS_ROOT)/src/ -I$(VCFLIB_ROOT)/src/ -I$(TABIX_ROOT)/ -I$(VCFLIB_ROOT)/smithwaterman/ -I$(VCFLIB_ROOT)/multichoose/ -I$(VCFLIB_ROOT)/filevercmp/ -I$(HTSLIB_ROOT) -I$(SEQLIB_ROOT) -I$(SEQLIB_ROOT)/htslib

all: autoversion ../bin/freebayes ../bin/bamleftalign

static:
	$(MAKE) CFLAGS="$(CFLAGS) -static" all

debug:
	$(MAKE) CFLAGS="$(CFLAGS) -D VERBOSE_DEBUG -g -rdynamic" all

profiling:
	$(MAKE) CFLAGS="$(CFLAGS) -g" all

gprof:
	$(MAKE) CFLAGS="$(CFLAGS) -pg" all

.PHONY: all static debug profiling gprof

# builds bamtools static lib, and copies into root
$(BAMTOOLS_ROOT)/lib/libbamtools.a:
	cd $(BAMTOOLS_ROOT) && mkdir -p build && cd build && cmake .. && $(MAKE)
$(HTSLIB_ROOT)/libhts.a:
	cd $(HTSLIB_ROOT) && make

$(SEQLIB_ROOT)/src/libseqlib.a:
	cd $(SEQLIB_ROOT) && ./configure && make 

OBJECTS=BedReader.o \
		CNV.o \
		fastlz.o \
		Fasta.o \
		Parameters.o \
		Allele.o \
		Sample.o \
		Result.o \
		AlleleParser.o \
		Utility.o \
		Genotype.o \
		DataLikelihood.o \
		Multinomial.o \
		Ewens.o \
		ResultData.o \
		Dirichlet.o \
		Marginals.o \
		split.o \
		LeftAlign.o \
		IndelAllele.o \
		Bias.o \
		Contamination.o \
		NonCall.o \
		SegfaultHandler.o \
		../vcflib/tabixpp/tabix.o \
		../vcflib/tabixpp/htslib/bgzf.o \
		../vcflib/smithwaterman/SmithWatermanGotoh.o \
		../vcflib/smithwaterman/disorder.cpp \
		../vcflib/smithwaterman/LeftAlign.o \
		../vcflib/smithwaterman/Repeats.o \
		../vcflib/smithwaterman/IndelAllele.o \
		Variant.o \
		$(BAMTOOLS_ROOT)/lib/libbamtools.a \
			$(SEQLIB_ROOT)/src/libseqlib.a	\
			$(SEQLIB_ROOT)/bwa/libbwa.a	\
			$(SEQLIB_ROOT)/fermi-lite/libfml.a	\
			$(SEQLIB_ROOT)/htslib/libhts.a

HEADERS=multichoose.h version_git.h

# executables

freebayes ../bin/freebayes: freebayes.o $(OBJECTS) $(HEADERS)
	$(CXX) $(CFLAGS) $(INCLUDE) freebayes.o $(OBJECTS) -o ../bin/freebayes $(LIBS)

alleles ../bin/alleles: alleles.o $(OBJECTS) $(HEADERS)
	$(CXX) $(CFLAGS) $(INCLUDE) alleles.o $(OBJECTS) -o ../bin/alleles $(LIBS)

dummy ../bin/dummy: dummy.o $(OBJECTS) $(HEADERS)
	$(CXX) $(CFLAGS) $(INCLUDE) dummy.o $(OBJECTS) -o ../bin/dummy $(LIBS)

bamleftalign ../bin/bamleftalign: $(BAMTOOLS_ROOT)/lib/libbamtools.a $(SEQLIB_ROOT)/src/libseqlib.a $(SEQLIB_ROOT)/htslib/libhts.a bamleftalign.o Fasta.o LeftAlign.o IndelAllele.o split.o
	$(CXX) $(CFLAGS) $(INCLUDE) bamleftalign.o Fasta.o Utility.o LeftAlign.o IndelAllele.o split.o $(BAMTOOLS_ROOT)/lib/libbamtools.a $(SEQLIB_ROOT)/src/libseqlib.a $(SEQLIB_ROOT)/htslib/libhts.a -o ../bin/bamleftalign $(LIBS)

bamfiltertech ../bin/bamfiltertech: $(BAMTOOLS_ROOT)/lib/libbamtools.a $(SEQLIB_ROOT)/src/libseqlib.a $(SEQLIB_ROOT)/htslib/libhts.a bamfiltertech.o $(OBJECTS) $(HEADERS)
	$(CXX) $(CFLAGS) $(INCLUDE) bamfiltertech.o $(OBJECTS) -o ../bin/bamfiltertech $(LIBS)


# objects

Fasta.o: Fasta.cpp Utility.o
	$(CXX) $(CFLAGS) $(INCLUDE) -c Fasta.cpp

alleles.o: alleles.cpp AlleleParser.o Allele.o
	$(CXX) $(CFLAGS) $(INCLUDE) -c alleles.cpp

dummy.o: dummy.cpp AlleleParser.o Allele.o
	$(CXX) $(CFLAGS) $(INCLUDE) -c dummy.cpp

freebayes.o: freebayes.cpp TryCatch.h $(HTSLIB_ROOT)/libhts.a $(BAMTOOLS_ROOT)/lib/libbamtools.a
	$(CXX) $(CFLAGS) $(INCLUDE) -c freebayes.cpp

fastlz.o: fastlz.c fastlz.h
	$(C) $(CFLAGS) $(INCLUDE) -c fastlz.c	

Parameters.o: Parameters.cpp Parameters.h Version.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c Parameters.cpp

Allele.o: Allele.cpp Allele.h multichoose.h Genotype.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c Allele.cpp

Sample.o: Sample.cpp Sample.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c Sample.cpp

Genotype.o: Genotype.cpp Genotype.h Allele.h multipermute.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c Genotype.cpp

Ewens.o: Ewens.cpp Ewens.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c Ewens.cpp

AlleleParser.o: AlleleParser.cpp AlleleParser.h multichoose.h Parameters.h $(BAMTOOLS_ROOT)/lib/libbamtools.a $(HTSLIB_ROOT)/libhts.a
	$(CXX) $(CFLAGS) $(INCLUDE) -c AlleleParser.cpp

Utility.o: Utility.cpp Utility.h Sum.h Product.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c Utility.cpp

SegfaultHandler.o: SegfaultHandler.cpp SegfaultHandler.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c SegfaultHandler.cpp

Dirichlet.o: Dirichlet.h Dirichlet.cpp Sum.h Product.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c Dirichlet.cpp

Multinomial.o: Multinomial.h Multinomial.cpp Sum.h Product.h Utility.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c Multinomial.cpp

DataLikelihood.o: DataLikelihood.cpp DataLikelihood.h Sum.h Product.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c DataLikelihood.cpp

Marginals.o: Marginals.cpp Marginals.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c Marginals.cpp

ResultData.o: ResultData.cpp ResultData.h Result.h Result.cpp Allele.h Utility.h Genotype.h AlleleParser.h Version.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c ResultData.cpp

Result.o: Result.cpp Result.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c Result.cpp

NonCall.o: NonCall.cpp NonCall.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c NonCall.cpp

BedReader.o: BedReader.cpp BedReader.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c BedReader.cpp

CNV.o: CNV.cpp CNV.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c CNV.cpp

Bias.o: Bias.cpp Bias.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c Bias.cpp

split.o: split.h split.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -c split.cpp

bamleftalign.o: bamleftalign.cpp LeftAlign.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -c bamleftalign.cpp

bamfiltertech.o: bamfiltertech.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -c bamfiltertech.cpp

LeftAlign.o: LeftAlign.h LeftAlign.cpp $(BAMTOOLS_ROOT)/lib/libbamtools.a $(HTSLIB_ROOT)/libhts.a
	$(CXX) $(CFLAGS) $(INCLUDE) -c LeftAlign.cpp

IndelAllele.o: IndelAllele.cpp IndelAllele.h
	$(CXX) $(CFLAGS) $(INCLUDE) -c IndelAllele.cpp

Variant.o: $(VCFLIB_ROOT)/src/Variant.h $(VCFLIB_ROOT)/src/Variant.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -c $(VCFLIB_ROOT)/src/Variant.cpp

../vcflib/tabixpp/tabix.o:
	cd $(TABIX_ROOT)/ && make
../vcflib/tabixpp/htslib/bgzf.o: ../vcflib/tabixpp/htslib/bgzf.c ../vcflib/tabixpp/htslib/htslib/bgzf.h
	cd ../vcflib/tabixpp && $(MAKE)

../vcflib/smithwaterman/SmithWatermanGotoh.o: ../vcflib/smithwaterman/SmithWatermanGotoh.h ../vcflib/smithwaterman/SmithWatermanGotoh.cpp
	cd ../vcflib/smithwaterman && $(MAKE)


VERSION_FILE=./version_git.h
RELEASED_VERSION_FILE=./version_release.txt


## versioning system from BedTools
##
## For freebayes developers (not users):
## When you want to release (and tag) a new version, run:
##   $ make setversion VERSION=v2.17.2
## This will:
##   1. Update the "/src/utils/version/version_release.txt" file
##   2. Commit the file
##   3. Git-Tag the commit with the latest version
##
.PHONY: setversion
setversion:
    ifeq "$(VERSION)" ""
		$(error please set VERSION variable to the new version (e.g "make setversion VERSION=v2.17.2"))
    endif
		@echo "# This file was auto-generated by running \"make setversion VERSION=$(VERSION)\"" > "$(RELEASED_VERSION_FILE)"
		@echo "# on $$(date) ." >> "$(RELEASED_VERSION_FILE)"
		@echo "# Please do not edit or commit this file manually." >> "$(RELEASED_VERSION_FILE)"
		@echo "#" >> "$(RELEASED_VERSION_FILE)"
		@echo "$(VERSION)" >> $(RELEASED_VERSION_FILE)
		@git add $(RELEASED_VERSION_FILE)
		@git commit -q -m "Setting Release-Version $(VERSION)"
		@git tag "$(VERSION)"
		@echo "Version updated to $(VERSION)."
		@echo ""
		@echo "Don't forget to push the commits AND the tags:"
		@echo "  git push --all --tags"
		@echo ""


## Automatic version detection
##
## What's going on here?
## 1. If there's a ".git" repository - use the version from the repository.
##    ignore any released-version file. git repository is authorative.
##
## 2, If there's no ".git" repository,
##    get the "released" version number from the release-version file.
##
## 3. Compare the detected version (from steps 1,2) to the current string
##    in ./src/utils/version/version_git.h .
##    If they differ, update the header file - will cause a recompilation
##    of version.o .
##
.PHONY: autoversion
autoversion:
	@( \
	if [ -d "../.git" ] && which git > /dev/null ; then \
		DETECTED_VERSION=$$(git describe --always --tags --dirty) ; \
	else \
		DETECTED_VERSION=$$(grep -v "^#" "$(RELEASED_VERSION_FILE)") ; \
	fi ; \
	\
	CURRENT_VERSION="" ; \
	[ -e "$(VERSION_FILE)" ] && CURRENT_VERSION=$$(grep "define VERSION_GIT " "$(VERSION_FILE)" | cut -f3 -d" " | sed 's/"//g') ; \
	\
	echo "DETECTED_VERSION = $$DETECTED_VERSION" ; \
	echo "CURRENT_VERSION  = $$CURRENT_VERSION" ; \
	if [ "$${DETECTED_VERSION}" != "$${CURRENT_VERSION}" ] ; then \
		echo "Updating version file." ; \
		echo "#ifndef VERSION_GIT_H" > $(VERSION_FILE) ; \
		echo "#define VERSION_GIT_H" >> $(VERSION_FILE) ; \
		echo "#define VERSION_GIT \"$${DETECTED_VERSION}\"" >> $(VERSION_FILE) ; \
		echo "#endif /* VERSION_GIT_H */" >> $(VERSION_FILE) ; \
	fi )


clean:
	rm -rf *.o *.cgh *~ freebayes alleles ../bin/freebayes ../bin/alleles ../vcflib/*.o ../vcflib/tabixpp/*.{o,a}
	cd $(BAMTOOLS_ROOT)/build && make clean
	cd ../vcflib/smithwaterman && make clean

