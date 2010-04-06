#ifndef _CALLER_H
#define _CALLER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <deque>
#include <boost/regex.hpp>
#include "BamReader.h"
#include "Class-BedReader.h"
#include "Parameters.h"
#include "Utility.h"
#include "Allele.h"
#include "Fasta.h"
#include "TryCatch.h"
#include "Function-Sequence.h"

using namespace std;
using namespace BamTools;

// a structure holding information about our parameters

// structure to encapsulate registered reads and alleles
class RegisteredAlignment {
    friend ostream &operator<<(ostream &out, RegisteredAlignment &a);
public:
    BamAlignment alignment;
    vector<Allele> alleles;
    int mismatches;

    RegisteredAlignment(BamAlignment alignment)
        : alignment(alignment)
        , mismatches(0)
    { }
};


class Caller {

public:

    Parameters* parameters; // holds operational parameters passed at program invocation
    
    Caller(int argc, char** argv);
    ~Caller(void); 

    vector<string> sampleList; // list of sample names, indexed by sample id

    // reference
    FastaReference* reference;
    vector<string> referenceSequenceNames;
    map<string, int> referenceSequenceNameToID;
    
    // target regions
    //vector<vector<BedData>> targetRegions;  // beddatas indexed by sequence id
    //map<string, vector<BedData> > targetsByRefseq; // same, indexed by sequence name
    vector<BedData> targets;

    // bamreader
    BamReader bamReader;

    deque<RegisteredAlignment> registeredAlignmentQueue;

    // reference names indexed by id
    vector<RefData> referenceSequences;
    // ^^ vector of objects containing:
    //RefName;          //!< Name of reference sequence
    //RefLength;        //!< Length of reference sequence
    //RefHasAlignments; //!< True if BAM file contains alignments mapped to reference sequence
 
    void openBam(void);
    void openLogFile(void);
    void getSampleNames(void);
    void loadBamReferenceSequenceNames(void);
    void loadFastaReference(void);
    void loadReferenceSequence(int seqID);
    void loadReferenceSequence(string seqName, int start, int length);
    void loadReferenceSequence(BedData*);
    void loadReferenceSequence(BedData*, int, int);
    void loadTargets(void);
    void initializeOutputFiles(void);
    RegisteredAlignment registerAlignment(BamAlignment& alignment);
    void updateAlignmentQueue(void);
    vector<BedData>* targetsInCurrentRefSeq(void);
    bool toNextRefID(void);
    bool loadTarget(BedData*);
    bool toFirstTargetPosition(void);
    bool toNextTargetPosition(void);
    bool toNextTarget(void);
    void setPosition(long unsigned int);
    int currentSequencePosition(const BamAlignment& alignment);
    bool getNextAlleles(vector<Allele>& alleles);
    void getAlleles(vector<Allele>& alleles);

    // math
    //long double logGenotypeLikelihood(vector<Allele>, string);

private:
    // output files
    ofstream rptFile, vcfFile, logFile;

    string currentSequence;
    int basesBeforeCurrentTarget; // number of bases in sequence we're storing before the current target
    int basesAfterCurrentTarget;  // ........................................  after ...................

    // pointer to current position in targets
    int fastaReferenceSequenceCount; // number of reference sequences
    int currentRefID;
    BedData* currentTarget;
    long unsigned int currentPosition;
    BamAlignment currentAlignment;

    // constants
    long double LOGFACTOR;  // log(10) / -10
    long double LN3;  // log 3

};

#endif
