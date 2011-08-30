################################################################################
# Makefile for FreeBayes
# Erik Garrison, 2010
# Boston College
################################################################################

# Compiler
CC=g++
C=gcc

# Compiler flags
CFLAGS=-O3
#CFLAGS=-O3 -static -D VERBOSE_DEBUG  # enables verbose debugging via --debug2

BAMTOOLS_ROOT=../bamtools

LIBS = -lz -lm -L./ -L../vcflib/tabixpp/ -L$(BAMTOOLS_ROOT)/lib -ltabix -lbamtools
INCLUDE = -I$(BAMTOOLS_ROOT)/src


all: ../bin/freebayes ../bin/bamleftalign

static:
	$(MAKE) CFLAGS="$(CFLAGS) -static" all

debug:
	$(MAKE) CFLAGS="$(CFLAGS) -D VERBOSE_DEBUG" all

profiling:
	$(MAKE) CFLAGS="$(CFLAGS) -g" all

gprof:
	$(MAKE) CFLAGS="$(CFLAGS) -pg" all

.PHONY: all static debug profiling gprof

# builds bamtools static lib, and copies into root
$(BAMTOOLS_ROOT)/lib/libbamtools.a:
	cd $(BAMTOOLS_ROOT) && mkdir -p build && cd build && cmake .. && $(MAKE)


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
		../vcflib/tabixpp/tabix.o \
		../vcflib/tabixpp/bgzf.o \
		../vcflib/smithwaterman/SmithWatermanGotoh.o \
		Variant.o

HEADERS=multichoose.h

# executables

freebayes ../bin/freebayes: freebayes.o $(BAMTOOLS_ROOT)/lib/libbamtools.a $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDE) freebayes.o $(OBJECTS) $(BAMTOOLS_ROOT)/lib/libbamtools.a -o ../bin/freebayes $(LIBS)

alleles ../bin/alleles: alleles.o $(BAMTOOLS_ROOT)/lib/libbamtools.a $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDE) alleles.o $(OBJECTS) $(BAMTOOLS_ROOT)/lib/libbamtools.a -o ../bin/alleles $(LIBS)

dummy ../bin/dummy: dummy.o $(BAMTOOLS_ROOT)/lib/libbamtools.a $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDE) dummy.o $(OBJECTS) $(BAMTOOLS_ROOT)/lib/libbamtools.a -o ../bin/dummy $(LIBS)

bamleftalign ../bin/bamleftalign: $(BAMTOOLS_ROOT)/lib/libbamtools.a bamleftalign.o Fasta.o LeftAlign.o IndelAllele.o split.o
	$(CC) $(CFLAGS) $(INCLUDE) $(BAMTOOLS_ROOT)/lib/libbamtools.a bamleftalign.o Fasta.o LeftAlign.o IndelAllele.o split.o -o ../bin/bamleftalign $(LIBS)

bamfiltertech ../bin/bamfiltertech: $(BAMTOOLS_ROOT)/lib/libbamtools.a bamfiltertech.o $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDE) bamfiltertech.o $(OBJECTS) $(BAMTOOLS_ROOT)/lib/libbamtools.a -o ../bin/bamfiltertech $(LIBS)


# objects

Fasta.o: Fasta.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c Fasta.cpp

alleles.o: alleles.cpp AlleleParser.o Allele.o
	$(CC) $(CFLAGS) $(INCLUDE) -c alleles.cpp

dummy.o: dummy.cpp AlleleParser.o Allele.o
	$(CC) $(CFLAGS) $(INCLUDE) -c dummy.cpp

freebayes.o: freebayes.cpp TryCatch.h $(BAMTOOLS_ROOT)/lib/libbamtools.a
	$(CC) $(CFLAGS) $(INCLUDE) -c freebayes.cpp

fastlz.o: fastlz.c fastlz.h
	$(C) $(CFLAGS) $(INCLUDE) -c fastlz.c	

Parameters.o: Parameters.cpp Parameters.h Version.h
	$(CC) $(CFLAGS) $(INCLUDE) -c Parameters.cpp

Allele.o: Allele.cpp Allele.h multichoose.h Genotype.h
	$(CC) $(CFLAGS) $(INCLUDE) -c Allele.cpp

Sample.o: Sample.cpp Sample.h
	$(CC) $(CFLAGS) $(INCLUDE) -c Sample.cpp

Genotype.o: Genotype.cpp Genotype.h Allele.h multipermute.h
	$(CC) $(CFLAGS) $(INCLUDE) -c Genotype.cpp

Ewens.o: Ewens.cpp Ewens.h
	$(CC) $(CFLAGS) $(INCLUDE) -c Ewens.cpp

AlleleParser.o: AlleleParser.cpp AlleleParser.h multichoose.h Parameters.h $(BAMTOOLS_ROOT)/lib/libbamtools.a
	$(CC) $(CFLAGS) $(INCLUDE) -c AlleleParser.cpp

Utility.o: Utility.cpp Utility.h Sum.h Product.h
	$(CC) $(CFLAGS) $(INCLUDE) -c Utility.cpp

Dirichlet.o: Dirichlet.h Dirichlet.cpp Sum.h Product.h
	$(CC) $(CFLAGS) $(INCLUDE) -c Dirichlet.cpp

Multinomial.o: Multinomial.h Multinomial.cpp Sum.h Product.h Utility.h
	$(CC) $(CFLAGS) $(INCLUDE) -c Multinomial.cpp

DataLikelihood.o: DataLikelihood.cpp DataLikelihood.h Sum.h Product.h
	$(CC) $(CFLAGS) $(INCLUDE) -c DataLikelihood.cpp

Marginals.o: Marginals.cpp Marginals.h
	$(CC) $(CFLAGS) $(INCLUDE) -c Marginals.cpp

ResultData.o: ResultData.cpp ResultData.h Result.h Result.cpp Allele.h Utility.h Genotype.h AlleleParser.h Version.h
	$(CC) $(CFLAGS) $(INCLUDE) -c ResultData.cpp

Result.o: Result.cpp Result.h
	$(CC) $(CFLAGS) $(INCLUDE) -c Result.cpp

BedReader.o: BedReader.cpp BedReader.h
	$(CC) $(CFLAGS) $(INCLUDE) -c BedReader.cpp

CNV.o: CNV.cpp CNV.h
	$(CC) $(CFLAGS) $(INCLUDE) -c CNV.cpp

split.o: split.h split.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c split.cpp

bamleftalign.o: bamleftalign.cpp LeftAlign.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c bamleftalign.cpp

bamfiltertech.o: bamfiltertech.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c bamfiltertech.cpp

LeftAlign.o: LeftAlign.h LeftAlign.cpp $(BAMTOOLS_ROOT)/lib/libbamtools.a
	$(CC) $(CFLAGS) $(INCLUDE) -c LeftAlign.cpp

IndelAllele.o: IndelAllele.cpp IndelAllele.h
	$(CC) $(CFLAGS) $(INCLUDE) -c IndelAllele.cpp

Variant.o: ../vcflib/Variant.h ../vcflib/Variant.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c ../vcflib/Variant.cpp

../vcflib/tabixpp/tabix.o: ../vcflib/tabixpp/tabix.hpp ../vcflib/tabixpp/tabix.cpp
../vcflib/tabixpp/bgzf.o: ../vcflib/tabixpp/bgzf.c ../vcflib/tabixpp/bgzf.h
	cd ../vcflib/tabixpp && $(MAKE)

../vcflib/smithwaterman/SmithWatermanGotoh.o: ../vcflib/smithwaterman/SmithWatermanGotoh.h ../vcflib/smithwaterman/SmithWatermanGotoh.cpp
	cd ../vcflib/smithwaterman && $(MAKE)

# versioning for git point releases
version:
	echo \#define FREEBAYES_VERSION \"`git log | head -1 | cut -f 2 -d\ `\" >Version.h
	echo \#define FREEBAYES_COMPILE_DATE \"`date +%F`\" >>Version.h

clean:
	rm -rf *.o *.cgh *~ freebayes alleles ../bin/freebayes ../bin/alleles ../vcflib/*.o ../vcflib/tabixpp/*.{o,a}
	cd $(BAMTOOLS_ROOT)/build && make clean
	cd ../vcflib/smithwaterman && make clean

