#ifndef _ALLELE_H
#define _ALLELE_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <limits>
#include <sstream>
#include <assert.h>
#include "Utility.h"
#include "convert.h"

//#ifdef HAVE_BAMTOOLS
//#include "api/BamAlignment.h"
//using namespace BamTools;
//#endif

using namespace std;

class Allele;

// Allele recycling allocator
// without we spend 30% of our runtime deleting Allele instances

class AlleleFreeList {

public:
    AlleleFreeList()
        : _p(NULL)
        , _size(0)
        , _allocs(0)
        , _min_size(0)
        , _max_size(0)
        , _tick_allocs(1000000)  // attempt realloc every million allocs
    { }

    ~AlleleFreeList();
    void Purge();
    void Resize(int new_size);
    void* NewAllele();
    void Recycle(void* mem);

private:
    Allele* _p;
    int _size;    // number of alleles on list
    int _allocs;  // allocation counter
    int _min_size; // min size within some previous number of calls to new
    int _max_size;
    int _tick_allocs; // GC cycle length

};


// a structure describing an allele

enum AlleleType {
    ALLELE_GENOTYPE = 1,
    ALLELE_REFERENCE = 2,
    ALLELE_MNP = 4,
    ALLELE_SNP = 8,
    ALLELE_INSERTION = 16,
    ALLELE_DELETION = 32,
    ALLELE_COMPLEX = 64,
    ALLELE_NULL = 128
};

// used in making allele type filter vectors
//const int numberOfPossibleAlleleTypes = 7;

enum AlleleStrand {
    STRAND_FORWARD,
    STRAND_REVERSE
};

typedef long int Position;

class Allele {

    friend class AlleleFreeList;

    friend string stringForAllele(const Allele &a);
    friend string stringForAlleles(vector<Allele> &av);

    friend bool operator<(const Allele &a, const Allele &b);
    friend bool operator==(const Allele &a, const Allele &b);
    friend bool operator!=(const Allele &a, const Allele &b);

    friend ostream &operator<<(ostream &out, vector<Allele> &a);
    friend ostream &operator<<(ostream &out, vector<Allele*> &a);
    friend ostream &operator<<(ostream &out, list<Allele*> &a);

    friend ostream &operator<<(ostream &out, Allele &a);
    friend ostream &operator<<(ostream &out, Allele* &a);

    friend string json(vector<Allele*> &alleles, long int &position);
    friend string json(vector<Allele*> &alleles);
    friend string json(Allele &allele, long int &position);
    friend string json(Allele* &allele);

public:

    AlleleType type;        // type of the allele, enumerated above
    string referenceName;   // reference name, for sanity checking
    string referenceSequence; // reference sequence or "" (in case of insertions)
    string alternateSequence; // alternate sequence or "" (in case of deletions and reference alleles)
    string sequencingTechnology; // the technology used to generate this allele
    long int position;      // position 0-based against reference
    long int* currentReferencePosition; // pointer to the current reference position (which may be updated during the life of this allele)
    char* currentReferenceBase;  // pointer to current reference base
    unsigned int length;    // and event length (deletion implies 0, snp implies 1, insertion >1)
    unsigned int referenceLength; // length of the event relative to the reference
    long int repeatRightBoundary;  // if this allele is an indel, and if it is embedded in a tandem repeat 
    // TODO cleanup
    int basesLeft;  // these are the "updated" versions of the above
    int basesRight;
    AlleleStrand strand;          // strand, true = +, false = -
    string sampleID;        // representative sample ID
    string readGroupID;     // read group membership
    string readID;          // id of the read which the allele is drawn from
    vector<short> baseQualities;
    long double quality;          // base quality score associated with this allele, updated every position in the case of reference alleles
    long double lnquality;  // log version of above
    string currentBase;       // current base, meant to be updated every position
    short mapQuality;       // map quality for the originating read
    long double lnmapQuality;       // map quality for the originating read
    double readMismatchRate; // per-base mismatch rate for the read
    double readIndelRate;  // only considering gaps
    double readSNPRate;    // only considering snps/mnps
    bool isProperPair;    // if the allele is supported by a properly paired read
    bool isPaired;  // if the allele is supported by a read that is part of a pair
    bool isMateMapped;  // if the mate in the pair is mapped
    bool genotypeAllele;    // if this is an abstract 'genotype' allele
    bool processed; // flag to mark if we've presented this allele for analysis
    string cigar; // a cigar representation of the allele
    vector<Allele>* alignmentAlleles;
    long int alignmentStart;
    long int alignmentEnd;

    // default constructor, for converting alignments into allele observations
    Allele(AlleleType t, 
           string& refname,
           long int pos, 
           long int* crefpos,
           char* crefbase,
           unsigned int len,
           long int rrbound,
           int bleft,
           int bright,
           string alt,
           string& sampleid,
           string& readid,
           string& readgroupid,
           string& sqtech,
           bool strnd, 
           long double qual,
           string qstr, 
           short mapqual,
           bool ispair,
           bool ismm,
           bool isproppair,
           string cigarstr,
           vector<Allele>* ra,
           long int bas,
           long int bae)
        : type(t)
        , referenceName(refname)
        , position(pos)
        , currentReferencePosition(crefpos)
        , currentReferenceBase(crefbase)
        , length(len)
        , repeatRightBoundary(rrbound)
        , basesLeft(bleft)
        , basesRight(bright)
        , currentBase(alt)
        , alternateSequence(alt)
        , sampleID(sampleid)
        , readID(readid)
        , readGroupID(readgroupid)
        , sequencingTechnology(sqtech)
        , strand(strnd ? STRAND_FORWARD : STRAND_REVERSE)
        , quality((qual == -1) ? averageQuality(qstr) : qual) // passing -1 as quality triggers this calculation
        , lnquality(phred2ln((qual == -1) ? averageQuality(qstr) : qual))
        , mapQuality(mapqual) 
        , lnmapQuality(phred2ln(mapqual))
        , isProperPair(isproppair)
        , isPaired(ispair)
        , isMateMapped(ismm)
        , genotypeAllele(false)
        , processed(false)
        , readMismatchRate(0)
        , readIndelRate(0)
        , readSNPRate(0)
        , cigar(cigarstr)
        , alignmentAlleles(ra)
        , alignmentStart(bas)
        , alignmentEnd(bae)
    {

        baseQualities.resize(qstr.size()); // cache qualities
        transform(qstr.begin(), qstr.end(), baseQualities.begin(), qualityChar2ShortInt);
        referenceLength = referenceLengthFromCigar();

    }

    // for constructing genotype alleles
    Allele(AlleleType t,
	   string alt,
	   unsigned int len,
	   unsigned int reflen,
	   string cigarStr,
	   long int pos=0,
	   long int rrbound=0,
	   bool gallele=true) 
        : type(t)
        , alternateSequence(alt)
        , length(len)
        , referenceLength(reflen)
        , repeatRightBoundary(rrbound)
        , quality(0)
        , lnquality(1)
        , position(pos)
        , genotypeAllele(true)
        , readMismatchRate(0)
        , readIndelRate(0)
        , readSNPRate(0)
        , cigar(cigarStr)
        , alignmentAlleles(NULL)
        , processed(false)
    {
        currentBase = base();
        baseQualities.assign(alternateSequence.size(), 0);
        referenceLength = referenceLengthFromCigar();
    }

    bool equivalent(Allele &a);  // heuristic 'equivalency' between two alleles, which depends on their type
    string typeStr(void) const; // return a string representation of the allele type, for output
    bool isReference(void) const; // true if type == ALLELE_REFERENCE
    bool isSNP(void) const; // true if type == ALLELE_SNP
    bool isInsertion(void) const; // true if type == ALLELE_INSERTION
    bool isDeletion(void) const; // true if type == ALLELE_DELETION
    bool isMNP(void) const; // true if type == ALLELE_MNP
    bool isComplex(void) const; // true if type == ALLELE_COMPLEX
    bool isNull(void) const; // true if type == ALLELE_NULL
    int referenceOffset(void) const;
    const short currentQuality(void) const;  // for getting the quality of a given position in multi-bp alleles
    const long double lncurrentQuality(void) const;
    const int subquality(int startpos, int len) const;
    const long double lnsubquality(int startpos, int len) const;
    const int subquality(const Allele &a) const;
    const long double lnsubquality(const Allele &a) const;
    //const int basesLeft(void) const; // returns the bases left within the read of the current position within the allele
    //const int basesRight(void) const; // returns the bases right within the read of the current position within the allele
    bool sameSample(Allele &other);  // if the other allele has the same sample as this one
    void update(int haplotypeLength = 1); // for reference alleles, updates currentBase and quality
    void setQuality(void); // sets 'current quality' for alleles
    // TODO update this to reflect different insertions (e.g. IATGC instead of I4)
    const string base(void) const;  // the 'current' base of the allele or a string describing the allele, e.g. I10 or D2
                                    //  this is used to update cached data in the allele prior to presenting the allele for analysis
                                    //  for the current base, just use allele.currentBase

    string json(void);
    unsigned int getLengthOnReference(void);
    int referenceLengthFromCigar(void);

    string readSeq(void);
    string read5p(void);
    string read3p(void);
    string read5pNonNull(void);
    string read3pNonNull(void);

    // the number of bases from the 5p edge of the allele until the end or the next null allele
    int read5pNonNullBases(void);
    // the number of bases from the 3p edge of the allele until the end or the next null allele
    int read3pNonNullBases(void);

    // wish list...
    //string readRefRelativeSubstr(long int start, long int end);
    //string readRefStartLenSubstr(long int start, int bp);

    vector<Allele*> extend(int pos, int haplotypeLength);
    void squash(void);
    void subtract(int subtractFromRefStart,
            int subtractFromRefEnd,
            string& substart,
            string& subend,
            vector<pair<int, string> >& cigarstart,
            vector<pair<int, string> >& cigarend,
            vector<short>& qsubstart,
            vector<short>& qsubend);

    void add(string& addToStart,
            string& addToEnd,
            vector<pair<int, string> >& cigarStart,
            vector<pair<int, string> >& cigarEnd,
            vector<short>& qaddToStart,
            vector<short>& qaddToEnd);


    void subtractFromStart(int bp, string& seq, vector<pair<int, string> >& cig, vector<short>& quals);
    void subtractFromEnd(int bp, string& seq, vector<pair<int, string> >& cig, vector<short>& quals);
    void addToStart(string& seq, vector<pair<int, string> >& cig, vector<short>& quals);
    void addToEnd(string& seq, vector<pair<int, string> >& cig, vector<short>& quals);

    void mergeAllele(const Allele& allele, AlleleType newType);

    void updateTypeAndLengthFromCigar(void);
    int bpLeft(void); // how many bases are in the read to the left of the allele
    int bpRight(void); // how many bases are in the read to the left of the allele


};

// for sorting pairs of alleles and ints
class AllelePairIntCompare {
public:
    bool operator()(const pair<Allele, int>& a, const pair<Allele, int>& b) {
        return a.second > b.second;
    }
};

class AllelePositionCompare {
public:
    bool operator()(const Allele& a, const Allele& b) {
        return a.position < b.position;
    }
};


void updateAllelesCachedData(vector<Allele*>& alleles);

map<string, vector<Allele*> > groupAllelesBySample(list<Allele*>& alleles);
void groupAllelesBySample(list<Allele*>& alleles, map<string, vector<Allele*> >& groups);

int allowedAlleleTypes(vector<AlleleType>& allowedEnumeratedTypes);
void filterAlleles(list<Allele*>& alleles, int allowedTypes);
int countAlleles(map<string, vector<Allele*> >& sampleGroups);
int baseCount(vector<Allele*>& alleles, string base, AlleleStrand strand);
pair<pair<int, int>, pair<int, int> >
baseCount(vector<Allele*>& alleles, string refbase, string altbase);
int countAllelesWithBase(vector<Allele*>& alleles, string base);

bool areHomozygous(vector<Allele*>& alleles);

map<Allele, int> countAlleles(vector<Allele*>& alleles);
map<string, int> countAllelesString(vector<Allele*>& alleles);
map<string, int> countAllelesString(vector<Allele>& alleles);
map<Allele, int> countAlleles(vector<Allele>& alleles);
map<Allele, int> countAlleles(list<Allele*>& alleles);

vector<Allele> uniqueAlleles(list<Allele*>& alleles);

bool allelesSameType(Allele* &a, Allele* &b);
bool allelesEquivalent(Allele* &a, Allele* &b);
bool allelesSameSample(Allele* &a, Allele* &b);
bool allelesSameType(Allele &a, Allele &b);
bool allelesEquivalent(Allele &a, Allele &b);
bool allelesSameSample(Allele &a, Allele &b);
bool allelesEqual(Allele &a, Allele &b);

void groupAlleles(map<string, vector<Allele*> >& sampleGroups, map<string, vector<Allele*> >& alleleGroups);
void homogenizeAlleles(map<string, vector<Allele*> >& alleleGroups, string& refseq, Allele& refallele);
void resetProcessedFlag(map<string, vector<Allele*> >& alleleGroups);

vector<Allele> alleleUnion(vector<Allele>& a1, vector<Allele>& a2);

// XXX cleanup
// is there a way to template these?  difficult as the syntax for pointer-based comparisons is different than non-pointer
vector<vector<Allele*> >  groupAlleles(list<Allele*> &alleles, bool (*fncompare)(Allele* &a, Allele* &b));
vector<vector<Allele*> >  groupAlleles(list<Allele> &alleles, bool (*fncompare)(Allele &a, Allele &b));
vector<vector<Allele*> > groupAlleles(vector<Allele*> &alleles, bool (*fncompare)(Allele &a, Allele &b));
vector<vector<Allele*> > groupAlleles(vector<Allele> &alleles, bool (*fncompare)(Allele &a, Allele &b));
vector<vector<Allele*> > groupAlleles(map<string, vector<Allele*> > &alleles, bool (*fncompare)(Allele &a, Allele &b));
vector<vector<Allele> > groupAlleles_copy(vector<Allele> &alleles, bool (*fncompare)(Allele &a, Allele &b));
vector<vector<Allele> > groupAlleles_copy(list<Allele> &alleles, bool (*fncompare)(Allele &a, Allele &b));
vector<vector<Allele> > groupAlleles_copy(vector<Allele> &alleles);
vector<Allele> genotypeAllelesFromAlleleGroups(vector<vector<Allele> > &groups);
vector<Allele> genotypeAllelesFromAlleleGroups(vector<vector<Allele*> > &groups);
vector<Allele> genotypeAllelesFromAlleles(vector<Allele> &alleles);
vector<Allele> genotypeAllelesFromAlleles(vector<Allele*> &alleles);
Allele genotypeAllele(Allele& a);
Allele genotypeAllele(AlleleType type, string alt = "", unsigned int length = 0, string cigar = "", unsigned int reflen = 0, long int position = 0, long int rrbound = 0);

bool isEmptyAllele(const Allele& allele);
bool isDividedIndel(const Allele& allele);
bool isEmptyAlleleOrIsDividedIndel(const Allele& allele);
bool isUnflankedIndel(const Allele& allele);

int referenceLengthFromCigar(string& cigar);

//AlleleFreeList Allele::_freeList;

#endif
