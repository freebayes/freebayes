################################################################################
# Makefile for the GigaBayes package
# Gabor T. Marth 2008, 2009
# Boston College
################################################################################

################################################################################
# To run (examples)
################################################################################

# optimized + static
# make CFLAGS="-O9 -static"

# to profile
# make CFLAGS="-O9 -static -pg"

################################################################################
# Variables
################################################################################

# Compiler
CC=g++
C=gcc

# Compiler flags
#CFLAGS=-g # debugging
#CFLAGS=-O3 #-pg #profiling
#CFLAGS=-pg -O3
CFLAGS=-O3

################################################################################
# Program builds
################################################################################
all: esfbayes alleles

esfbayes: esfbayes.o BamMultiReader.o BamReader.o BamIndex.o BGZF.o Class-BedReader.o Function-Sequence.o fastlz.o Fasta.o Parameters.o Allele.o AlleleParser.o Utility.o multichoose.h Genotype.o DataLikelihood.o Multinomial.o GenotypePriors.o ResultData.o
	$(CC) $(CFLAGS) esfbayes.o BamMultiReader.o BamReader.o BamIndex.o BGZF.o Class-BedReader.o Function-Sequence.o Parameters.o fastlz.o Fasta.o Allele.o AlleleParser.o Utility.o Genotype.o DataLikelihood.o Multinomial.o GenotypePriors.o ResultData.o -o esfbayes -lz -lboost_regex -lpthread -Wno-deprecated

alleles: alleles.o BamMultiReader.o BamReader.o BamIndex.o BGZF.o Class-BedReader.o Function-Sequence.o fastlz.o Fasta.o Parameters.o Allele.o AlleleParser.o Utility.o multichoose.h
	$(CC) $(CFLAGS) alleles.o BamMultiReader.o BamReader.o BamIndex.o BGZF.o Class-BedReader.o Function-Sequence.o Parameters.o fastlz.o Fasta.o Allele.o AlleleParser.o Utility.o  -o alleles -lz -lboost_regex -lpthread -Wno-deprecated

################################################################################
# Compilation
################################################################################

# Main program compilations
Fasta.o: Fasta.cpp
	$(CC) $(CFLAGS) -c Fasta.cpp -Wno-deprecated

alleles.o: alleles.cpp AlleleParser.o Allele.o
	$(CC) $(CFLAGS) -c alleles.cpp -Wno-deprecated

esfbayes.o: esfbayes.cpp TryCatch.h
	$(CC) $(CFLAGS) -c esfbayes.cpp -Wno-deprecated


# Class compilations
fastlz.o: fastlz.c fastlz.h
	$(C) $(CFLAGS) -c fastlz.c	

BamMultiReader.o: BamMultiReader.cpp BamMultiReader.h BamReader.cpp BamReader.h BamAux.h
	$(CC) $(CFLAGS) -c BamMultiReader.cpp  -Wno-deprecated	

BamReader.o: BamReader.cpp BamReader.h BamAux.h
	$(CC) $(CFLAGS) -c BamReader.cpp  -Wno-deprecated	

BamWriter.o: BamWriter.cpp BamReader.h BamAux.h
	$(CC) $(CFLAGS) -c BamWriter.cpp  -Wno-deprecated	

BamIndex.o: BamIndex.cpp BamIndex.h BamAux.h
	$(CC) $(CFLAGS) -c BamIndex.cpp  -Wno-deprecated	

BGZF.o: BGZF.cpp BGZF.h
	$(CC) $(CFLAGS) -c BGZF.cpp -Wno-deprecated

Parameters.o: Parameters.cpp Parameters.h
	$(CC) $(CFLAGS) -c Parameters.cpp -Wno-deprecated

Allele.o: Allele.cpp Allele.h multichoose.h Genotype.h
	$(CC) $(CFLAGS) -c Allele.cpp -Wno-deprecated

Genotype.o: Genotype.cpp Genotype.h Allele.h multipermute.h
	$(CC) $(CFLAGS) -c Genotype.cpp -Wno-deprecated

GenotypePriors.o: GenotypePriors.cpp GenotypePriors.h Allele.h
	$(CC) $(CFLAGS) -c GenotypePriors.cpp -Wno-deprecated

AlleleParser.o: AlleleParser.cpp AlleleParser.h multichoose.h Parameters.cpp Parameters.h
	$(CC) $(CFLAGS) -c AlleleParser.cpp -Wno-deprecated

Utility.o: Utility.cpp Utility.h Sum.h Product.h
	$(CC) $(CFLAGS) -c Utility.cpp -Wno-deprecated

Dirichlet.o: Dirichlet.h Dirichlet.cpp Sum.h Product.h
	$(CC) $(CFLAGS) -c Dirichlet.cpp -Wno-deprecated

Multinomial.o: Multinomial.h Multinomial.cpp Sum.h Product.h Utility.h
	$(CC) $(CFLAGS) -c Multinomial.cpp -Wno-deprecated

DataLikelihood.o: DataLikelihood.cpp DataLikelihood.h Sum.h Product.h
	$(CC) $(CFLAGS) -c DataLikelihood.cpp -Wno-deprecated

ResultData.o: ResultData.cpp ResultData.h
	$(CC) $(CFLAGS) -c ResultData.cpp -Wno-deprecated

# ...
Class-BedReader.o: Class-BedReader.cpp Class-BedReader.h Function-Generic.h
	$(CC) $(CFLAGS) -c Class-BedReader.cpp -Wno-deprecated

# Function compilations
Function-Sequence.o: Function-Sequence.cpp Function-Sequence.h
	$(CC) $(CFLAGS) -c Function-Sequence.cpp -Wno-deprecated

Function-Generic.o: Function-Generic.cpp Function-Generic.h Type-Hash.h
	$(CC) $(CFLAGS) -c Function-Generic.cpp -Wno-deprecated

Function-Math.o: Function-Math.cpp Function-Math.h Type-Hash.h
	$(CC) $(CFLAGS) -c Function-Math.cpp -Wno-deprecated

################################################################################
# Archival
################################################################################
# TODO
#
################################################################################
# Installation
################################################################################

################################################################################
# House-keeping
################################################################################
clean:
	rm -rf *.o *.cgh *~ esfbayes alleles 

cleanEmacs:
	rm -rf *~

