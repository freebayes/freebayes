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
#include <time.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include "split.h"
#include "join.h"
#include "BamReader.h"
#include "BedReader.h"
#include "Parameters.h"
#include "Utility.h"
#include "Allele.h"
#include "Fasta.h"
#include "TryCatch.h"
#include "BamMultiReader.h"
#include "Genotype.h"

#define CACHED_REFERENCE_WINDOW 10

using namespace std;
using namespace BamTools;

// a structure holding information about our parameters

// structure to encapsulate registered reads and alleles
class RegisteredAlignment {
    friend ostream &operator<<(ostream &out, RegisteredAlignment &a);
public:
    //BamAlignment alignment;
    long unsigned int start;
    long unsigned int end;
    int refid;
    string name;
    vector<Allele> alleles;
    int mismatches;

    RegisteredAlignment(BamAlignment& alignment)
        //: alignment(alignment)
        : start(alignment.Position)
        , end(alignment.GetEndPosition())
        , refid(alignment.RefID)
        , name(alignment.Name)
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
    map<int, string> referenceIDToName;
    
    // target regions
    //vector<vector<BedTarget>> targetRegions;  // beddatas indexed by sequence id
    //map<string, vector<BedTarget> > targetsByRefseq; // same, indexed by sequence name
    vector<BedTarget> targets;

    // bamreader
    BamMultiReader bamMultiReader;

    map<long unsigned int, deque<RegisteredAlignment> > registeredAlignments;
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
    void openTraceFile(void);
    void openFailedFile(void);
    void openOutputFile(void);
    void getSampleNames(void);
    void loadBamReferenceSequenceNames(void);
    void loadFastaReference(void);
    void loadReferenceSequence(BedTarget*, int, int);
    void loadReferenceSequence(BamAlignment& alignment);
    void preserveReferenceSequenceWindow(int bp);
    void extendReferenceSequence(int);
    void extendReferenceSequence(BamAlignment& alignment);
    void eraseReferenceSequence(int leftErasure);
    string referenceSubstr(long double position, unsigned int length);
    void loadTargets(void);
    bool getFirstAlignment(void);
    void loadTargetsFromBams(void);
    void initializeOutputFiles(void);
    RegisteredAlignment& registerAlignment(BamAlignment& alignment, RegisteredAlignment& ra, string sampleName);
    void clearRegisteredAlignments(void);
    void updateAlignmentQueue(void);
    void removeNonOverlappingAlleles(vector<Allele*>& alleles);
    void removeFilteredAlleles(vector<Allele*>& alleles);
    void updateRegisteredAlleles(void);
    vector<BedTarget>* targetsInCurrentRefSeq(void);
    bool toNextRefID(void);
    bool loadTarget(BedTarget*);
    bool toFirstTargetPosition(void);
    bool toNextPosition(void);
    bool dummyProcessNextTarget(void);
    bool toNextTarget(void);
    void setPosition(long unsigned int);
    int currentSequencePosition(const BamAlignment& alignment);
    bool getNextAlleles(Samples& allelesBySample, int allowedAlleleTypes);
    void getAlleles(Samples& allelesBySample, int allowedAlleleTypes);
    Allele* referenceAllele(int mapQ, int baseQ);
    Allele* alternateAllele(int mapQ, int baseQ);
    int homopolymerRunLeft(string altbase);
    int homopolymerRunRight(string altbase);

    // gets the genotype alleles we should evaluate among the allele groups and
    // sample groups at the current position, according to our filters
    vector<Allele> genotypeAlleles(map<string, vector<Allele*> >& alleleGroups,
            Samples& samples,
            vector<Allele>& allGenotypeAlleles);

    // pointer to current position in targets
    int fastaReferenceSequenceCount; // number of reference sequences
    BedTarget* currentTarget;
    long double currentPosition;  // 0-based current position
    char currentReferenceBase;
    string currentSequence;
    char currentReferenceBaseChar();
    string currentReferenceBaseString();
    string::iterator currentReferenceBaseIterator();

    // output files
    ofstream logFile, outputFile, traceFile, failedFile;
    ostream* output;

    // utility
    bool isCpG(string& altbase);

private:

    bool justSwitchedTargets;  // to trigger clearing of queues, maps and such holding Allele*'s on jump

    Allele* currentReferenceAllele;
    Allele* currentAlternateAllele;

    //BedTarget currentSequenceBounds;
    long int currentSequenceStart;
    string currentSequenceName;

    bool hasMoreAlignments;

    int basesBeforeCurrentTarget; // number of bases in sequence we're storing before the current target
    int basesAfterCurrentTarget;  // ........................................  after ...................

    int currentRefID;
    BamAlignment currentAlignment;

};

#endif
