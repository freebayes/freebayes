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
CFLAGS=-g # debugging

################################################################################
# Program builds
################################################################################
all: testMath bayes bamBayes bamIndel bamDup bamDump bamSlice bamStat rptSum aggregateHapmap

testMath: testMath.o Function-Sequence.o Function-Math.o
	$(CC) $(CFLAGS) testMath.o Function-Math.o Function-Sequence.o -o testMath -lboost_regex -lpthread -Wno-deprecated

bamBayes: bamBayes.o ReferenceSequenceReader.o BamReader.o BGZF.o Class-BedReader.o Class-FastaReader.o Function-Sequence.o Function-Generic.o Function-Math.o fastlz.o Fasta.o
	$(CC) $(CFLAGS) bamBayes.o ReferenceSequenceReader.o BamReader.o BGZF.o Class-BedReader.o Class-FastaReader.o Function-Sequence.o Function-Generic.o Function-Math.o  fastlz.o Fasta.o -o bamBayes -lz -lboost_regex -lpthread -Wno-deprecated

bayes: bayes.o BamReader.o BGZF.o Class-BedReader.o Class-FastaReader.o Function-Sequence.o Function-Generic.o Function-Math.o fastlz.o Fasta.o Parameters.o Allele.o Caller.o Utility.o multichoose.h
	$(CC) $(CFLAGS) bayes.o BamReader.o BGZF.o Class-BedReader.o Class-FastaReader.o Function-Sequence.o Function-Generic.o Function-Math.o Parameters.o fastlz.o Fasta.o Allele.o Caller.o Utility.o  -o bayes -lz -lboost_regex -lpthread -Wno-deprecated

tests: tests.o Allele.o
	$(CC) $(CFLAGS) tests.o Allele.o Utility.o  -o tests

bamIndel: bamIndel.o ReferenceSequenceReader.o BamReader.o BGZF.o Class-BedReader.o Class-FastaReader.o Function-Sequence.o Function-Generic.o Function-Math.o fastlz.o
	$(CC) $(CFLAGS) bamIndel.o ReferenceSequenceReader.o BamReader.o BGZF.o Class-BedReader.o Class-FastaReader.o Function-Sequence.o Function-Generic.o Function-Math.o  fastlz.o -o bamIndel -lz -lboost_regex -lpthread -Wno-deprecated

bamDup: bamDup.o BamReader.o BamWriter.o BGZF.o Function-Generic.o
	$(CC) $(CFLAGS) bamDup.o BamReader.o BamWriter.o BGZF.o Function-Generic.o -o bamDup -lz -lboost_regex -lpthread -Wno-deprecated

bamDump: bamDump.o BamReader.o BGZF.o 
	$(CC) $(CFLAGS) bamDump.o BamReader.o BGZF.o -o bamDump -lz -lboost_regex -lpthread -Wno-deprecated

bamSlice: bamSlice.o ReferenceSequenceReader.o BamReader.o BGZF.o Function-Generic.o Function-Sequence.o Function-Math.o fastlz.o
	$(CC) $(CFLAGS) bamSlice.o ReferenceSequenceReader.o BamReader.o BGZF.o Function-Generic.o Function-Sequence.o Function-Math.o fastlz.o -o bamSlice -lz -lboost_regex -lpthread -Wno-deprecated

bamStat: bamStat.o BamReader.o Function-Generic.o Class-BedReader.o
	$(CC) $(CFLAGS) bamStat.o BamReader.o BGZF.o Function-Generic.o Class-BedReader.o -o bamStat -lz -lboost_regex -lpthread -Wno-deprecated

rptSum: rptSum.o Function-Generic.o Class-BedReader.o
	$(CC) $(CFLAGS) rptSum.o Function-Generic.o  Class-BedReader.o -o rptSum -lboost_regex -lpthread -Wno-deprecated

aggregateHapmap: aggregateHapmap.o ReferenceSequenceReader.o Function-Generic.o fastlz.o
	$(CC) $(CFLAGS) aggregateHapmap.o ReferenceSequenceReader.o Function-Generic.o -o aggregateHapmap -lboost_regex  fastlz.o -lpthread -Wno-deprecated

################################################################################
# Compilation
################################################################################

# Main program compilations
testMath.o: testMath.cpp
	$(CC) $(CFLAGS) -c testMath.cpp -Wno-deprecated

Fasta.o: Fasta.cpp
	$(CC) $(CFLAGS) -c Fasta.cpp -Wno-deprecated

bamBayes.o: bamBayes.cpp TryCatch.h
	$(CC) $(CFLAGS) -c bamBayes.cpp -Wno-deprecated

bayes.o: bayes.cpp TryCatch.h Caller.o
	$(CC) $(CFLAGS) -c bayes.cpp -Wno-deprecated

bamIndel.o: bamIndel.cpp
	$(CC) $(CFLAGS) -c bamIndel.cpp -Wno-deprecated

bamDup.o: bamDup.cpp
	$(CC) $(CFLAGS) -c bamDup.cpp -Wno-deprecated

bamDump.o: bamDump.cpp
	$(CC) $(CFLAGS) -c bamDump.cpp -Wno-deprecated

bamSlice.o: bamSlice.cpp
	$(CC) $(CFLAGS) -c bamSlice.cpp -Wno-deprecated

bamStat.o: bamStat.cpp
	$(CC) $(CFLAGS) -c bamStat.cpp -Wno-deprecated

rptSum.o: rptSum.cpp
	$(CC) $(CFLAGS) -c rptSum.cpp -Wno-deprecated

aggregateHapmap.o: aggregateHapmap.cpp
	$(CC) $(CFLAGS) -c aggregateHapmap.cpp -Wno-deprecated

# Class compilations
Class-GigReader.o: Class-GigReader.cpp Class-GigReader.h Function-Sequence.h Function-Generic.h LargeFileSupport.h 
	$(CC) $(CFLAGS) -c Class-GigReader.cpp  -Wno-deprecated

ReferenceSequenceReader.o: ReferenceSequenceReader.cpp ReferenceSequenceReader.h
	$(CC) $(CFLAGS) -c ReferenceSequenceReader.cpp  -Wno-deprecated	

fastlz.o: fastlz.c fastlz.h
	$(C) $(CFLAGS) -c fastlz.c	

BamReader.o: BamReader.cpp BamReader.h BamAux.h
	$(CC) $(CFLAGS) -c BamReader.cpp  -Wno-deprecated	

BamWriter.o: BamWriter.cpp BamReader.h BamAux.h
	$(CC) $(CFLAGS) -c BamWriter.cpp  -Wno-deprecated	

BGZF.o: BGZF.cpp BGZF.h
	$(CC) $(CFLAGS) -c BGZF.cpp -Wno-deprecated

Class-BedReader.o: Class-BedReader.cpp Class-BedReader.h Function-Generic.h
	$(CC) $(CFLAGS) -c Class-BedReader.cpp -Wno-deprecated

Class-FastaReader.o: Class-FastaReader.cpp Class-FastaReader.h Function-Generic.h
	$(CC) $(CFLAGS) -c Class-FastaReader.cpp -Wno-deprecated

Parameters.o: Parameters.cpp Parameters.h
	$(CC) $(CFLAGS) -c Parameters.cpp -Wno-deprecated

tests.o: tests.cpp
	$(CC) $(CFLAGS) -c tests.cpp

Allele.o: Allele.cpp Allele.h multichoose.h
	$(CC) $(CFLAGS) -c Allele.cpp -Wno-deprecated

Caller.o: Caller.cpp Caller.h multichoose.h Allele.cpp Allele.h
	$(CC) $(CFLAGS) -c Caller.cpp -Wno-deprecated

Utility.o: Utility.cpp Utility.h
	$(CC) $(CFLAGS) -c Utility.cpp -Wno-deprecated

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
tarSource: 	Makefile gigaBuild.cpp gigaBayes.cpp gigaDump.cpp testMath.cpp \
		bamBayes.cpp bamIndel.cpp bamSlice.cpp bamStat.cpp bamDup.cpp bamDump.cpp \
		rptSum.cpp aggregateHapmap.cpp \
		Class-GigReader.h Class-GigReader.cpp \
		Function-Sequence.h Function-Sequence.cpp \
		Function-Generic.h Function-Generic.cpp \
		Function-Math.h Function-Math.cpp \
		Class-BedReader.h Class-BedReader.cpp \
		Class-FastaReader.h Class-FastaReader.cpp \
		BamAux.h BamReader.h BamReader.cpp \
		fastlz.h fastlz.c \
		ReferenceSequenceReader.h ReferenceSequenceReader.cpp \
		BamWriter.h BamWriter.cpp \
		LargeFileSupport.h Type-Hash.h README Makefile snpFilter.pl snpRate.pl
	tar cvf BamBayes.tar Makefile gigaBuild.cpp gigaBayes.cpp gigaDump.cpp testMath.cpp \
		bamBayes.cpp bamIndel.cpp bamSlice.cpp bamStat.cpp bamDup.cpp bamDump.cpp \
		rptSum.cpp aggregateHapmap.cpp \
		Class-GigReader.h Class-GigReader.cpp \
		Function-Sequence.h Function-Sequence.cpp \
		Function-Generic.h Function-Generic.cpp \
		Function-Math.h Function-Math.cpp \
		Class-BedReader.h Class-BedReader.cpp \
		Class-FastaReader.h Class-FastaReader.cpp \
		BamAux.h BamReader.h BamReader.cpp \
		fastlz.h  fastlz.c \
		ReferenceSequenceReader.h ReferenceSequenceReader.cpp \
		BamWriter.h BamWriter.cpp \
		LargeFileSupport.h Type-Hash.h README Makefile snpFilter.pl snpRate.pl; \
	gzip BamBayes.tar

tarExecutable:	testMath \
		bamBayes bamIndel bamDup bamDump bamStat bamSlice \
		rptSum aggregateHapmap \
		snpFilter.pl snpRate.pl
	tar cvf BamBayesBin.tar testMath \
		bamBayes bamIndel bamDup bamDump bamStat bamSlice \
		rptSum aggregateHapmap \
		snpFilter.pl snpRate.pl; \
	gzip BamBayesBin.tar

################################################################################
# Installation
################################################################################

install:	testMath \
		bamBayes bamIndel bamDup bamDump bamStat bamSlice \
		rptSum aggregateHapmap \
		snpFilter.pl snpRate.pl
	cp testMath ../bin/; \
	cp bamBayes ../bin/; \
	cp bamIndel ../bin/; \
	cp bamDup ../bin/; \
	cp bamDump ../bin/; \
	cp bamStat ../bin/; \
	cp bamSlice ../bin/; \
	cp rptSum ../bin/; \
	cp aggregateHapmap ../bin/; \
	cp snpFilter.pl ../bin/; \
	cp snpRate.pl ../bin/;

################################################################################
# House-keeping
################################################################################
clean:
	rm -rf *.o *.cgh *~ gigaBuild gigaBayes gigaDump testMath \
		bamBayes bamIndel bamDup bamDump bamStat bamSlice \
		rptStat rptSum aggregateHapmap

cleanEmacs:
	rm -rf *~

