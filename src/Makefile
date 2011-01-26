################################################################################
# Makefile for FreeBayes
# Erik Garrison, 2010
# Boston College
################################################################################

# optimized + static
# make CFLAGS="-O9 -static"

# to profile
# make CFLAGS="-O9 -static -pg"

# to enable verbose debugging
# make CFLAGS="-O9 -static -pg -D VERBOSE_DEBUG"

# Compiler
CC=g++
C=gcc

# Compiler flags
CFLAGS=-O3 -Wno-deprecated
#CFLAGS=-O3 -static -D VERBOSE_DEBUG  # enables verbose debugging via --debug2
LIBS = -lz -lm


all: ../bin/freebayes ../bin/bamleftalign ../bin/bamfiltertech

static:
	$(MAKE) CFLAGS="$(CFLAGS) -static" all

debug:
	$(MAKE) CFLAGS="$(CFLAGS) -D VERBOSE_DEBUG" all

profiling:
	$(MAKE) CFLAGS="$(CFLAGS) -g" all

gprof:
	$(MAKE) CFLAGS="$(CFLAGS) -pg" all

.PHONY: all static debug profiling gprof

OBJECTS=BamAlignment.o \
		BamMultiReader.o \
		BamReader.o \
		BamWriter.o \
		BamIndex.o \
		BGZF.o \
		BedReader.o \
		CNV.o \
		fastlz.o \
		Fasta.o \
		Parameters.o \
		Allele.o \
		AlleleParser.o \
		Utility.o \
		Genotype.o \
		DataLikelihood.o \
		Multinomial.o \
		GenotypePriors.o \
		ResultData.o \
		Dirichlet.o \
		Marginals.o \
		split.o \
		LeftAlign.o \
		IndelAllele.o

HEADERS=multichoose.h

# executables

freebayes ../bin/freebayes: freebayes.o $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) freebayes.o $(OBJECTS) -o ../bin/freebayes $(LIBS)

alleles ../bin/alleles: alleles.o $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) alleles.o $(OBJECTS) -o ../bin/alleles $(LIBS)

dummy ../bin/dummy: dummy.o $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) dummy.o $(OBJECTS) -o ../bin/dummy $(LIBS)

bamleftalign ../bin/bamleftalign: bamleftalign.o $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) bamleftalign.o $(OBJECTS) -o ../bin/bamleftalign $(LIBS)

bamfiltertech ../bin/bamfiltertech: bamfiltertech.o $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) bamfiltertech.o $(OBJECTS) -o ../bin/bamfiltertech $(LIBS)


# objects

Fasta.o: Fasta.cpp
	$(CC) $(CFLAGS) -c Fasta.cpp

alleles.o: alleles.cpp AlleleParser.o Allele.o
	$(CC) $(CFLAGS) -c alleles.cpp

dummy.o: dummy.cpp AlleleParser.o Allele.o
	$(CC) $(CFLAGS) -c dummy.cpp

freebayes.o: freebayes.cpp TryCatch.h
	$(CC) $(CFLAGS) -c freebayes.cpp

fastlz.o: fastlz.c fastlz.h
	$(C) $(CFLAGS) -c fastlz.c	

BamAlignment.o: BamAlignment.cpp BamAlignment.h BamAux.h
	$(CC) $(CFLAGS) -c BamAlignment.cpp 

BamMultiReader.o: BamMultiReader.cpp BamMultiReader.h BamReader.cpp BamReader.h BamAux.h
	$(CC) $(CFLAGS) -c BamMultiReader.cpp 	

BamReader.o: BamReader.cpp BamReader.h BamAux.h
	$(CC) $(CFLAGS) -c BamReader.cpp 	

BamWriter.o: BamWriter.cpp BamReader.h BamAux.h
	$(CC) $(CFLAGS) -c BamWriter.cpp 	

BamIndex.o: BamIndex.cpp BamIndex.h BamAux.h
	$(CC) $(CFLAGS) -c BamIndex.cpp 	

BGZF.o: BGZF.cpp BGZF.h
	$(CC) $(CFLAGS) -c BGZF.cpp

Parameters.o: Parameters.cpp Parameters.h Version.h
	$(CC) $(CFLAGS) -c Parameters.cpp

Allele.o: Allele.cpp Allele.h multichoose.h Genotype.h
	$(CC) $(CFLAGS) -c Allele.cpp

Genotype.o: Genotype.cpp Genotype.h Allele.h multipermute.h
	$(CC) $(CFLAGS) -c Genotype.cpp

GenotypePriors.o: GenotypePriors.cpp GenotypePriors.h Allele.h
	$(CC) $(CFLAGS) -c GenotypePriors.cpp

AlleleParser.o: AlleleParser.cpp AlleleParser.h multichoose.h Parameters.h
	$(CC) $(CFLAGS) -c AlleleParser.cpp

Utility.o: Utility.cpp Utility.h Sum.h Product.h
	$(CC) $(CFLAGS) -c Utility.cpp

Dirichlet.o: Dirichlet.h Dirichlet.cpp Sum.h Product.h
	$(CC) $(CFLAGS) -c Dirichlet.cpp

Multinomial.o: Multinomial.h Multinomial.cpp Sum.h Product.h Utility.h
	$(CC) $(CFLAGS) -c Multinomial.cpp

DataLikelihood.o: DataLikelihood.cpp DataLikelihood.h Sum.h Product.h
	$(CC) $(CFLAGS) -c DataLikelihood.cpp

Marginals.o: Marginals.cpp Marginals.h
	$(CC) $(CFLAGS) -c Marginals.cpp

ResultData.o: ResultData.cpp ResultData.h Allele.h Utility.h Genotype.h AlleleParser.h Version.h
	$(CC) $(CFLAGS) -c ResultData.cpp

BedReader.o: BedReader.cpp BedReader.h
	$(CC) $(CFLAGS) -c BedReader.cpp

CNV.o: CNV.cpp CNV.h
	$(CC) $(CFLAGS) -c CNV.cpp

split.o: split.h split.cpp
	$(CC) $(CFLAGS) -c split.cpp

bamleftalign.o: bamleftalign.cpp LeftAlign.cpp
	$(CC) $(CFLAGS) -c bamleftalign.cpp

bamfiltertech.o: bamfiltertech.cpp
	$(CC) $(CFLAGS) -c bamfiltertech.cpp

LeftAlign.o: LeftAlign.h LeftAlign.cpp
	$(CC) $(CFLAGS) -c LeftAlign.cpp

IndelAllele.o: IndelAllele.cpp IndelAllele.h
	$(CC) $(CFLAGS) -c IndelAllele.cpp

# versioning for git point releases
version:
	echo \#define FREEBAYES_VERSION \"`git log | head -1 | cut -f 2 -d\ `\" >Version.h
	echo \#define FREEBAYES_COMPILE_DATE \"`date +%F`\" >>Version.h

clean:
	rm -rf *.o *.cgh *~ freebayes alleles ../bin/freebayes ../bin/alleles

