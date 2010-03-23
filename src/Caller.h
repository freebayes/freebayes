#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include "BamReader.h"
#include "Class-BedReader.h"
#include "Parameters.h"
#include "Utility.h"

using namespace std;
using namespace BamTools;

// a structure holding information about our parameters

// structure to encapsulate registered reads and alleles
class RegisteredAlignment {
public:
    BamAlignment alignment;
    vector<Allele> alleles;
    int mismatches;

    RegisteredAlignment(BamAlignment alignment)
        : alignment = alignment
        , mismatches = 0
    { }
}


class BayesCaller {

public:

    Parameters params; // holds operational parameters passed at program invocation
    BayesCaller(Parameters& params);

    
    vector<string> sampleList;

    // reference
    FastaReference reference;
    vector<string> referenceSequenceNames;
    map<string, int> referenceSequenceNameToID;
    
    // target regions
    //vector<vector<BedData>> targetRegions;  // beddatas indexed by sequence id
    map<string, vector<BedData>> targetsByRefseq; // same, indexed by sequence name

    // bamreader
    BamReader bamReader;
    // bed regions
    BedReader bedReader;

    ofstream logFile;

    //queue<BamAlignment> alignmentQueue;
    //queue<Allele> alleleQueue;
    queue<RegisteredAlignment> registeredAlignmentQueue;

    // reference names indexed by id
    vector<RefData> referenceSequences;
    //map<string, RefData> referenceSequencesByName;
    // ^^ vector of objects containing:
    //RefName;          //!< Name of reference sequence
    //RefLength;        //!< Length of reference sequence
    //RefHasAlignments; //!< True if BAM file contains alignments mapped to reference sequence
    //

private:
    // output files
    ofstream rptFile, vcfFile;

    string currentSequence;

    // pointer to current position in targets
    int referenceSequenceCount; // number of reference sequences
    int currentRefID;
    BedData* currentTarget;
    long unsigned int currentPosition;
    BamAlignment currentAlignment;
    string currentSequence;

}
