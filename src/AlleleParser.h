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
#include <cmath>
#include "split.h"
#include "join.h"

#include "BedReader.h"
#include "Parameters.h"
#include "Utility.h"
#include "Allele.h"
#include "Sample.h"
#include "Fasta.h"
#include "TryCatch.h"

#include "Genotype.h"
#include "CNV.h"
#include "Result.h"
#include "LeftAlign.h"
#include "Variant.h"
#include "version_git.h"

// the size of the window of the reference which is always cached in memory
#define CACHED_REFERENCE_WINDOW 300

// the window of haplotype basis alleles which we ensure we keep
// increasing this reduces disk access when using haplotype basis alleles, but increases memory usage
#define CACHED_BASIS_HAPLOTYPE_WINDOW 1000

using namespace std;

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
    string readgroup;
    vector<Allele> alleles;
    int mismatches;
    int snpCount;
    int indelCount;
    int alleleTypes;
    Parameters parameters;

    RegisteredAlignment(BAMALIGN& alignment, Parameters parameters)
        //: alignment(alignment)
        : start(alignment.POSITION)
        , end(alignment.ENDPOSITION)
        , refid(alignment.REFID)
        , name(alignment.QNAME)
        , mismatches(0)
        , snpCount(0)
        , indelCount(0)
        , alleleTypes(0)
        , parameters(parameters)
    {
      FILLREADGROUP(readgroup, alignment);
    }

    void addAllele(Allele allele, bool mergeComplex = true,
                   int maxComplexGap = 0, bool boundIndels = false);
    bool fitHaplotype(int pos, int haplotypeLength, Allele*& aptr, bool allowPartials = false);

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

class AllelicPrimitive {
public:
    string alt;
    string ref;
    AllelicPrimitive(string& r, string& a)
	: ref(r)
	, alt(a) { }
};

bool operator<(const AllelicPrimitive& a, const AllelicPrimitive& b);

void capBaseQuality(BAMALIGN& alignment, int baseQualityCap);

class AlleleParser {

public:

    Parameters parameters; // holds operational parameters passed at program invocation

    AlleleParser(int argc, char** argv);
    ~AlleleParser(void);

    vector<string> sampleList; // list of sample names, indexed by sample id
    vector<string> sampleListFromBam; // sample names drawn from BAM file
    vector<string> sampleListFromVCF; // sample names drawn from input VCF
    map<string, string> samplePopulation; // population subdivisions of samples
    map<string, vector<string> > populationSamples; // inversion of samplePopulation
    map<string, string> readGroupToSampleNames; // maps read groups to samples
    map<string, string> readGroupToTechnology; // maps read groups to technologies
    vector<string> sequencingTechnologies;  // a list of the present technologies

    CNVMap sampleCNV;

    // reference
    FastaReference reference;
    vector<string> referenceSequenceNames;
    map<int, string> referenceIDToName;
    string referenceSampleName;

    // target regions
    vector<BedTarget> targets;
    // returns true if we are within a target
    // useful for controlling output when we are reading from stdin
    bool inTarget(void);

    // bamreader
    BAMREADER bamMultiReader;

    // bed reader
    BedReader bedReader;

    // VCF
    vcflib::VariantCallFile variantCallFile;
    vcflib::VariantCallFile variantCallInputFile;   // input variant alleles, to target analysis
    vcflib::VariantCallFile haplotypeVariantInputFile;  // input alleles which will be used to construct haplotype alleles

    // input haplotype alleles
    //
    // as calling progresses, a window of haplotype basis alleles from the flanking sequence
    // map from starting position to length->alle
    map<long int, vector<AllelicPrimitive> > haplotypeBasisAlleles;  // this is in the current reference sequence
    bool usingHaplotypeBasisAlleles;
    bool usingVariantInputAlleles;
    long int rightmostHaplotypeBasisAllelePosition;
    long int rightmostInputAllelePosition;
    void updateHaplotypeBasisAlleles(long int pos, int referenceLength);
    bool allowedHaplotypeBasisAllele(long int pos, string& ref, string& alt);

    Allele makeAllele(RegisteredAlignment& ra,
		      AlleleType type,
		      long int pos,
		      int length,
		      int basesLeft,
		      int basesRight,
		      string& readSequence,
		      string& sampleName,
		      BAMALIGN& alignment,
		      string& sequencingTech,
		      long double qual,
		      string& qualstr);



    vector<Allele*> registeredAlleles;
    map<long unsigned int, deque<RegisteredAlignment> > registeredAlignments;
    map<int, map<long int, vector<Allele> > > inputVariantAlleles; // all variants present in the input VCF, as 'genotype' alleles
    pair<int, long int> nextInputVariantPosition(void);
    void getInputVariantsInRegion(string& seq, long start = 0, long end = 0);
    void getAllInputVariants(void);
    //  position         sample     genotype  likelihood
    map<string, map<long int, map<string, map<string, long double> > > > inputGenotypeLikelihoods; // drawn from input VCF
    map<string, map<long int, map<Allele, int> > > inputAlleleCounts; // drawn from input VCF
    Sample* nullSample;

    bool loadNextPositionWithAlignmentOrInputVariant(BAMALIGN& currentAlignment);
    bool loadNextPositionWithInputVariant(void);
    bool hasMoreInputVariants(void);

    void addCurrentGenotypeLikelihoods(map<int, vector<Genotype> >& genotypesByPloidy,
            vector<vector<SampleDataLikelihood> >& sampleDataLikelihoods);

    void getInputAlleleCounts(vector<Allele>& genotypeAlleles, map<string, int>& inputAFs);

    // reference names indexed by id
    REFVEC referenceSequences;
    // ^^ vector of objects containing:
    //RefName;          //!< Name of reference sequence
    //RefLength;        //!< Length of reference sequence
    //RefHasAlignments; //!< True if BAM file contains alignments mapped to reference sequence

    vector<string> bamHeaderLines;

    void openBams(void);
    void openOutputFile(void);
    void getSampleNames(void);
    void getPopulations(void);
    void getSequencingTechnologies(void);
    void loadSampleCNVMap(void);
    int currentSamplePloidy(string const& sample);
    int copiesOfLocus(Samples& samples);
    vector<int> currentPloidies(Samples& samples);
    void loadBamReferenceSequenceNames(void);
    void loadFastaReference(void);
    void loadReferenceSequence(BAMALIGN& alignment);
    void loadReferenceSequence(string& seqname);
    string referenceSubstr(long int position, unsigned int length);
    void loadTargets(void);
    bool getFirstAlignment(void);
    bool getFirstVariant(void);
    void loadTargetsFromBams(void);
    void initializeOutputFiles(void);
    RegisteredAlignment& registerAlignment(BAMALIGN& alignment, RegisteredAlignment& ra, string& sampleName, string& sequencingTech);
    void clearRegisteredAlignments(void);
    void updateAlignmentQueue(long int position, vector<Allele*>& newAlleles, bool gettingPartials = false);
    void updateInputVariants(long int pos, int referenceLength);
    void updateHaplotypeBasisAlleles(void);
    void removeAllelesWithoutReadSpan(vector<Allele*>& alleles, int probeLength, int haplotypeLength);
    void removeNonOverlappingAlleles(vector<Allele*>& alleles,
                                     int haplotypeLength = 1,
                                     bool getAllAllelesInHaplotype = false);
    void removePreviousAlleles(vector<Allele*>& alleles);
    void removeFilteredAlleles(vector<Allele*>& alleles);
    void removeDuplicateAlleles(Samples& samples, map<string, vector<Allele*> >& alleleGroups,
                                int allowedAlleleTypes, int haplotypeLength, Allele& refallele);
    void updateRegisteredAlleles(void);
    void addToRegisteredAlleles(vector<Allele*>& alleles);
    void updatePriorAlleles(void);
    vector<BedTarget>* targetsInCurrentRefSeq(void);
    bool toNextRefID(void);
    bool loadTarget(BedTarget*);
    bool toFirstTargetPosition(void);
    bool toNextPosition(void);
    void getCompleteObservationsOfHaplotype(Samples& samples, int haplotypeLength, vector<Allele*>& haplotypeObservations);
    void getPartialObservationsOfHaplotype(Samples& samples, int haplotypeLength, vector<Allele*>& partials);
    bool dummyProcessNextTarget(void);
    bool toNextTarget(void);
    void setPosition(long unsigned int);
    int currentSequencePosition(const BAMALIGN& alignment);
    int currentSequencePosition();
    void unsetAllProcessedFlags(void);
    bool getNextAlleles(Samples& allelesBySample, int allowedAlleleTypes);

    // builds up haplotype (longer, e.g. ref+snp+ref) alleles to match the longest allele in genotypeAlleles
    // updates vector<Allele>& alleles with the new alleles
    void buildHaplotypeAlleles(vector<Allele>& alleles,
                               Samples& allelesBySample,
                               map<string, vector<Allele*> >& alleleGroups,
                               // provides observation group counts, counts of partial observations
                               map<string, vector<Allele*> >& partialObservationGroups,
                               map<Allele*, set<Allele*> >& partialObservationSupport,
                               int allowedAlleleTypes);
    void getAlleles(Samples& allelesBySample,
                    int allowedAlleleTypes,
                    int haplotypeLength = 1,
                    bool getAllAllelesInHaplotype = false,
                    bool ignoreProcessedAlleles = true);
    Allele* referenceAllele(int mapQ, int baseQ);
    Allele* alternateAllele(int mapQ, int baseQ);
    int homopolymerRunLeft(string altbase);
    int homopolymerRunRight(string altbase);
    map<string, int> repeatCounts(long int position, const string& sequence, int maxsize);
    map<long int, map<string, int> > cachedRepeatCounts; // cached version of previous
    bool isRepeatUnit(const string& seq, const string& unit);
    void setupVCFOutput(void);
    void setupVCFInput(void);
    string vcfHeader(void);
    bool hasInputVariantAllelesAtCurrentPosition(void);

    // gets the genotype alleles we should evaluate among the allele groups and
    // sample groups at the current position, according to our filters
    vector<Allele> genotypeAlleles(map<string, vector<Allele*> >& alleleGroups,
                                   Samples& samples,
                                   bool useOnlyInputAlleles,
                                   int haplotypeLength = 1);

    // pointer to current position in targets
    int fastaReferenceSequenceCount; // number of reference sequences
    bool hasTarget;
    BedTarget* currentTarget;
    long int currentPosition;  // 0-based current position
    int lastHaplotypeLength;
    char currentReferenceBase;
    string currentSequence;
    char currentReferenceBaseChar();
    string currentReferenceBaseString();
    string::iterator currentReferenceBaseIterator();
    string currentReferenceHaplotype();

    // output files
    ofstream logFile, outputFile;
    ostream* output;

    // utility
    bool isCpG(string& altbase);

    string currentSequenceName;

private:

    bool justSwitchedTargets;  // to trigger clearing of queues, maps and such holding Allele*'s on jump

    Allele* currentReferenceAllele;
    Allele* currentAlternateAllele;

    //BedTarget currentSequenceBounds;
    long int currentSequenceStart;

    bool hasMoreAlignments;
    bool hasMoreVariants;;

    bool oneSampleAnalysis; // if we are analyzing just one sample, and there are no specified read groups

    int basesBeforeCurrentTarget; // number of bases in sequence we're storing before the current target
    int basesAfterCurrentTarget;  // ........................................  after ...................

    int currentRefID;
    BAMALIGN currentAlignment;
    vcflib::Variant* currentVariant;

};

#endif
