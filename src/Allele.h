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
    ALLELE_COMPLEX = 64
    //ALLELE_CNV = 128,
};

// used in making allele type filter vectors
//const int numberOfPossibleAlleleTypes = 7;

enum AlleleStrand {
    STRAND_FORWARD,
    STRAND_REVERSE
};

typedef long double Position;

class Allele {

    friend class AlleleFreeList;

    friend string stringForAllele(Allele &a);
    friend string stringForAlleles(vector<Allele> &av);

    friend bool operator<(const Allele &a, const Allele &b);
    friend bool operator==(const Allele &a, const Allele &b);
    friend bool operator!=(const Allele &a, const Allele &b);

    friend ostream &operator<<(ostream &out, vector<Allele> &a);
    friend ostream &operator<<(ostream &out, vector<Allele*> &a);
    friend ostream &operator<<(ostream &out, list<Allele*> &a);

    friend ostream &operator<<(ostream &out, Allele &a);
    friend ostream &operator<<(ostream &out, Allele* &a);

    friend string json(vector<Allele*> &alleles, long double &position);
    friend string json(vector<Allele*> &alleles);
    friend string json(Allele &allele, long double &position);
    friend string json(Allele* &allele);

public:

    AlleleType type;        // type of the allele, enumerated above
    string referenceName;   // reference name, for sanity checking
    string referenceSequence; // reference sequence or "" (in case of insertions)
    string alternateSequence; // alternate sequence or "" (in case of deletions and reference alleles)
    string sequencingTechnology; // the technology used to generate this allele
    long double position;      // position 0-based against reference
    long double* currentReferencePosition; // pointer to the current reference position (which may be updated during the life of this allele)
    char* currentReferenceBase;  // pointer to current reference base
    unsigned int length;    // and event length (deletion implies 0, snp implies 1, insertion >1)
    unsigned int referenceLength; // length of the event relative to the reference
    int bpLeft; // how many bases are in the read to the left of the allele
    int bpRight; // how many bases are in the read to the left of the allele
    // TODO cleanup
    int basesLeft;  // these are the "updated" versions of the above
    int basesRight;
    AlleleStrand strand;          // strand, true = +, false = -
    string sampleID;        // representative sample ID
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
    vector<bool> indelMask; // indel mask structure, masks sites within the IDW from indels
    const bool masked(void) const;      // if the allele is masked at the *currentReferencePosition
    bool processed; // flag to mark if we've presented this allele for analysis

    // default constructor, for converting alignments into allele observations
    Allele(AlleleType t, 
                string& refname,
                long double pos, 
                long double* crefpos,
                char* crefbase,
                unsigned int len, 
                int bleft,
                int bright,
                string alt, 
                string& sampleid,
                string& readid,
                string& sqtech,
                bool strnd, 
                long double qual,
                string qstr, 
                short mapqual,
                bool ispair,
                bool ismm,
                bool isproppair)
        : type(t)
        , referenceName(refname)
        , position(pos)
        , currentReferencePosition(crefpos)
        , currentReferenceBase(crefbase)
        , length(len)
        , bpLeft(bleft)
        , basesLeft(bleft)
        , bpRight(bright)
        , basesRight(bright)
        , currentBase(alt)
        , alternateSequence(alt)
        , sampleID(sampleid)
        , readID(readid)
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
    {

        baseQualities.resize(qstr.size()); // cache qualities
        transform(qstr.begin(), qstr.end(), baseQualities.begin(), qualityChar2ShortInt);
        referenceLength = getLengthOnReference();

    }

    // for constructing genotype alleles
    Allele(AlleleType t,
            string alt,
            unsigned int len,
            unsigned int reflen,
            long double pos=0,
            bool gallele=true) 
        : type(t)
        , alternateSequence(alt)
        , length(len)
        , referenceLength(reflen)
        , quality(0)
        , lnquality(1)
        , position(pos)
        , genotypeAllele(true)
        , readMismatchRate(0)
        , readIndelRate(0)
        , readSNPRate(0)
    {
        currentBase = base();
    }

    /*
    Allele(const Allele& other) 
        : type(other.type)
        , referenceName(other.referenceName)
        , position(other.position)
        , currentReferencePosition(other.currentReferencePosition)
        , currentReferenceBase(other.currentReferenceBase)
        , length(other.length)
        , referenceSequence(other.referenceSequence)
        , alternateSequence(other.alternateSequence)
        , sampleID(other.sampleID)
        , readID(other.readID)
        , strand(other.strand)
        , quality(other.quality)
        , lnquality(other.lnquality)
        , currentBase(other.currentBase)
        , baseQualities(other.baseQualities)
        , mapQuality(other.mapQuality) 
        , genotypeAllele(other.genotypeAllele)
        , processed(other.processed)
    { }
    */

    bool equivalent(Allele &a);  // heuristic 'equivalency' between two alleles, which depends on their type
    string typeStr(void); // return a string representation of the allele type, for output
    bool isReference(void); // true if type == ALLELE_REFERENCE
    bool isSNP(void); // true if type == ALLELE_SNP
    bool isInsertion(void); // true if type == ALLELE_INSERTION
    bool isDeletion(void); // true if type == ALLELE_DELETION
    bool isMNP(void); // true if type == ALLELE_MNP
    bool isComplex(void); // true if type == ALLELE_COMPLEX
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
    void update(void); // for reference alleles, updates currentBase and quality
    // TODO update this to reflect different insertions (e.g. IATGC instead of I4)
    const string base(void) const;  // the 'current' base of the allele or a string describing the allele, e.g. I10 or D2
                                    //  this is used to update cached data in the allele prior to presenting the allele for analysis
                                    //  for the current base, just use allele.currentBase

    string json(void);
    unsigned int getLengthOnReference(void);



    void mergeAllele(const Allele& allele);

    // overload new and delete for object recycling pool

    /*
    void* operator new (size_t size) {
        assert (size == sizeof(Allele));
        return _freeList.NewAllele();
    }

    void operator delete (void* mem) {
        if (mem) _freeList.Recycle(mem);  // check allows us to safely subclass Allele
    }

    static void Purge() { _freeList.Purge(); }

    Allele* pNext() const { return _pNext; } // next allocated allele

private:
    static AlleleFreeList _freeList; // for object recycling
    Allele* _pNext; // for use in object recycling on the AlleleFreeList, used to link freed alleles into a list
    */

};

// for sorting pairs of alleles and ints
class AllelePairIntCompare {
public:
    bool operator()(const pair<Allele, int>& a, const pair<Allele, int>& b) {
        return a.second > b.second;
    }
};



void updateAllelesCachedData(vector<Allele*>& alleles);

map<string, vector<Allele*> > groupAllelesBySample(list<Allele*>& alleles);
void groupAllelesBySample(list<Allele*>& alleles, map<string, vector<Allele*> >& groups);

int allowedAlleleTypes(vector<AlleleType>& allowedEnumeratedTypes);
void filterAlleles(list<Allele*>& alleles, int allowedTypes);
void removeIndelMaskedAlleles(list<Allele*>& alleles, long double position);
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
Allele genotypeAllele(AlleleType type, string alt = "", unsigned int length = 0, unsigned int reflen = 0, long double position = 0);


//AlleleFreeList Allele::_freeList;

#endif
