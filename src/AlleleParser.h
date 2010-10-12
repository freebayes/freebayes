#ifndef _ALLELE_PARSER_H
#define _ALLELE_PARSER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <deque>
#include <utility>
#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <boost/bind.hpp>
#include <time.h>
#include "split.h"
#include "BamReader.h"
#include "BedReader.h"
#include "Parameters.h"
#include "Utility.h"
#include "Allele.h"
#include "Fasta.h"
#include "TryCatch.h"
//#include "Function-Math.h"
#include "Function-Sequence.h"
#include "BamMultiReader.h"
#include "Genotype.h"

using namespace std;
using namespace BamTools;
using boost::tuple;
using boost::make_tuple;

// a structure holding information about our parameters

// structure to encapsulate registered reads and alleles
class RegisteredAlignment {
    friend ostream &operator<<(ostream &out, RegisteredAlignment &a);
public:
    BamAlignment alignment;
    vector<Allele*> alleles;
    int mismatches;

    RegisteredAlignment(BamAlignment alignment)
        : alignment(alignment)
        , mismatches(0)
    { }

};

// functor to filter alleles outside of our analysis window
class AlleleFilter {

public:
    AlleleFilter(long unsigned int s, long unsigned int e) : start(s), end(e) {}

    // true of the allele is outside of our window
    bool operator()(Allele& a) { 
        return !(start >= a.position && end < a.position + a.length);
    }

    bool operator()(Allele*& a) { 
        return !(start >= a->position && end < a->position + a->length);
    }

private:
    long unsigned int start, end;

};

class AllelePtrCmp {

public:
    bool operator()(Allele* &a, Allele* &b) {
        return a->type < b->type;
    }

};


class AlleleParser {

public:

    Parameters parameters; // holds operational parameters passed at program invocation
    
    AlleleParser(int argc, char** argv);
    ~AlleleParser(void); 

    void writeVcfHeader(ostream& out);

    vector<string> sampleList; // list of sample names, indexed by sample id
    vector<string> sampleListFromBam; // sample names drawn from BAM file
    map<string, string> readGroupToSampleNames; // maps read groups to samples

    // reference
    FastaReference* reference;
    vector<string> referenceSequenceNames;
    map<string, int> referenceSequenceNameToID;
    
    // target regions
    //vector<vector<BedTarget>> targetRegions;  // beddatas indexed by sequence id
    //map<string, vector<BedTarget> > targetsByRefseq; // same, indexed by sequence name
    vector<BedTarget> targets;

    // bamreader
    BamMultiReader bamMultiReader;

    deque<RegisteredAlignment> registeredAlignmentQueue;
    vector<Allele*> registeredAlleles;
    //list<Allele*> registeredAlleles;
    //map<string, list<Allele*> > allelesBySample;

    // reference names indexed by id
    vector<RefData> referenceSequences;
    // ^^ vector of objects containing:
    //RefName;          //!< Name of reference sequence
    //RefLength;        //!< Length of reference sequence
    //RefHasAlignments; //!< True if BAM file contains alignments mapped to reference sequence
 
    void openBams(void);
    void openLogFile(void);
    void openTraceFile(void);
    void openOutputFile(void);
    void getSampleNames(void);
    void loadBamReferenceSequenceNames(void);
    void loadFastaReference(void);
    void loadReferenceSequence(int seqID);
    void loadReferenceSequence(string seqName, int start, int length);
    void loadReferenceSequence(BedTarget*);
    void loadReferenceSequence(BedTarget*, int, int);
    void extendReferenceSequence(int);
    void loadTargets(void);
    void initializeOutputFiles(void);
    RegisteredAlignment registerAlignment(BamAlignment& alignment, string sampleName);
    void updateAlignmentQueue(void);
    void removeNonOverlappingAlleles(vector<Allele*>& alleles);
    void removeFilteredAlleles(vector<Allele*>& alleles);
    void updateRegisteredAlleles(void);
    vector<BedTarget>* targetsInCurrentRefSeq(void);
    bool toNextRefID(void);
    bool loadTarget(BedTarget*);
    bool toFirstTargetPosition(void);
    bool toNextTargetPosition(void);
    bool dummyProcessNextTarget(void);
    bool toNextTarget(void);
    void setPosition(long unsigned int);
    int currentSequencePosition(const BamAlignment& alignment);
    bool getNextAlleles(map<string, vector<Allele*> >& allelesBySample, int allowedAlleleTypes);
    void getAlleles(map<string, vector<Allele*> >& allelesBySample, int allowedAlleleTypes);
    Allele* referenceAllele(int mapQ, int baseQ);
    Allele* alternateAllele(int mapQ, int baseQ);
    int homopolymerRunLeft(string altbase);
    int homopolymerRunRight(string altbase);

    // gets the genotype alleles we should evaluate among the allele groups and
    // sample groups at the current position, according to our filters
    vector<Allele> genotypeAlleles(vector<vector<Allele*> >& alleleGroups,
            map<string, vector<Allele*> >& sampleGroups,
            vector<Allele>& allGenotypeAlleles);

    // TODO clean these up...
    // p( observedAlleles | genotype ) for all genotypes
    vector<pair<Genotype, long double> > probObservedAllelesGivenGenotypes(vector<Allele*> &observedAlleles, vector< vector<Allele> > &genotypes);
    vector<pair<Genotype, long double> > probObservedAllelesGivenGenotypes_approx(vector<Allele*> &observedAlleles, vector< vector<Allele> > &genotypes);
    vector<pair<Genotype, long double> > probObservedAllelesGivenPossibleGenotypes(vector<Allele*> &observedAlleles, int ploidy);
    long double probObservedAllelesGivenGenotype(vector<vector<Allele*> > &alleleCombo, vector<Allele> &genotype);
    long double probObservedAllelesGivenGenotype(vector<tuple<long double, long double, vector<Allele*> > > &alleleComboProbs, 
            vector<Allele> &genotype);

    // NB: this is the one we're currently using
    long double probObservedAllelesGivenGenotypeAndTrueAlleles(
            vector<tuple<long double, long double, Allele* > >& allelesAndProbs,
            vector<Allele*>& trueAlleles,
            vector<Allele>& genotype);

    // caches probability products of the allele observations
    void calculateAlleleGroupProbabilities(vector<vector<Allele*> >& alleleGroups, 
            vector<tuple<long double, long double, vector<Allele*> > >& alleleGroupsAndQualities);


    // calculates log observation true/false probabilities and caches them in allelesAndProbs
    // true prob = log(1 - phred2float(Q score)), false prob = phred2log(Q score)
    void calculateAlleleBinaryProbabilities(
            vector<Allele*>& alleles, 
            vector<tuple<long double, long double, Allele* > >& allelesAndProbs);

    // pointer to current position in targets
    int fastaReferenceSequenceCount; // number of reference sequences
    BedTarget* currentTarget;
    long unsigned int currentPosition;  // 0-based current position
    char currentReferenceBase;
    string currentSequence;
    char currentReferenceBaseChar();
    string currentReferenceBaseString();
    string::iterator currentReferenceBaseIterator();

    // output files
    ofstream logFile, outputFile, traceFile;
    ostream* output;

private:

    bool justSwitchedTargets;  // to trigger clearing of queues, maps and such holding Allele*'s on jump

    Allele* currentReferenceAllele;
    Allele* currentAlternateAllele;

    int basesBeforeCurrentTarget; // number of bases in sequence we're storing before the current target
    int basesAfterCurrentTarget;  // ........................................  after ...................

    int currentRefID;
    BamAlignment currentAlignment;

};


#endif
