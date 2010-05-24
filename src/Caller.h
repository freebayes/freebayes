#ifndef _CALLER_H
#define _CALLER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <deque>
#include <utility>
#include <algorithm>
#include <boost/regex.hpp>
#include <boost/tuple/tuple.hpp>
#include "BamReader.h"
#include "Class-BedReader.h"
#include "Parameters.h"
#include "Utility.h"
#include "Allele.h"
#include "Fasta.h"
#include "TryCatch.h"
//#include "Function-Math.h"
#include "Function-Sequence.h"
#include "BamMultiReader.h"

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


class Caller {

public:

    Parameters parameters; // holds operational parameters passed at program invocation
    
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
    void updateRegisteredAlleles(void);
    vector<BedData>* targetsInCurrentRefSeq(void);
    bool toNextRefID(void);
    bool loadTarget(BedData*);
    bool toFirstTargetPosition(void);
    bool toNextTargetPosition(void);
    bool toNextTarget(void);
    void setPosition(long unsigned int);
    int currentSequencePosition(const BamAlignment& alignment);
    bool getNextAlleles(list<Allele*>& alleles);
    void getAlleles(list<Allele*>& alleles);

    // p( observedAlleles | genotype ) for all genotypes
    vector<pair<Genotype, long double> > probObservedAllelesGivenGenotypes(vector<Allele*> &observedAlleles, vector< vector<Allele> > &genotypes);
    vector<pair<Genotype, long double> > probObservedAllelesGivenPossibleGenotypes(vector<Allele*> &observedAlleles, int ploidy);
    long double probObservedAllelesGivenGenotype(vector<vector<Allele*> > &alleleCombo, vector<Allele> &genotype);
    long double probObservedAllelesGivenGenotype(vector<tuple<long double, long double, vector<Allele*> > > &alleleComboProbs, 
            vector<Allele> &genotype);
    // caches probability products of the allele observations
    void calculateAlleleGroupProbabilities(vector<vector<Allele*> >& alleleGroups, 
            vector<tuple<long double, long double, vector<Allele*> > >& alleleGroupsAndQualities);

    // pointer to current position in targets
    int fastaReferenceSequenceCount; // number of reference sequences
    BedData* currentTarget;
    long unsigned int currentPosition;  // 0-based current position

private:
    // output files
    ofstream rptFile, vcfFile, logFile;

    string currentSequence;
    int basesBeforeCurrentTarget; // number of bases in sequence we're storing before the current target
    int basesAfterCurrentTarget;  // ........................................  after ...................

    int currentRefID;
    BamAlignment currentAlignment;

};


// helper functions:

// probability of drawing an allele from the set of alleles
long double probChooseAlleleFromAlleles(Allele &allele, vector<Allele> &alleles);
int observationsInAlleleCombo(vector<vector<Allele> > &combo);
long double approximateBayesianNormalizationFactor(vector<vector<Allele> > &genotypes, 
        vector<vector<long double> > &probGenotypesGivenSampleObs,
        vector<vector<Allele> > &sampleGroups);
long double bayesianNormalizationFactor(vector<vector<Allele> > &genotypes, 
        vector<vector<long double> > &probGenotypesGivenSampleObs,
        vector<vector<Allele> > &sampleGroups);
vector<pair<long double, vector<Allele> > > mostLikelyGenotypesGivenObservations(vector<vector<Allele> > &genotypeCombos, 
                vector<vector<long double> > &probsBySample, bool normalize = false);
void normalizeGenotypeProbabilities(vector<pair<Genotype, long double> >& genotypeProbabilities);


#endif
