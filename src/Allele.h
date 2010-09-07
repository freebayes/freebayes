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
#include <boost/pool/object_pool.hpp>

using namespace std;

class Allele;

// Allele recycling allocator
// without we spend 30% of our runtime deleting Allele instances

class AlleleFreeList {

public:
    AlleleFreeList() : _p(NULL), _allocs(0) { }
    ~AlleleFreeList();
    void Purge();
    void* NewAllele();
    void Recycle(void* link);

private:
    Allele* _p;
    int _allocs;

};


// a structure describing an allele

enum AlleleType {
    ALLELE_GENOTYPE = 1,
    ALLELE_REFERENCE = 2,
    ALLELE_MISMATCH = 4,
    ALLELE_SNP = 8,
    ALLELE_INSERTION = 16,
    ALLELE_DELETION = 32,
    ALLELE_CNV = 64
};

// used in making allele type filter vectors
//const int numberOfPossibleAlleleTypes = 7;

enum AlleleStrand {
    STRAND_FORWARD,
    STRAND_REVERSE
};

typedef long unsigned int Position;

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

    friend string json(vector<Allele*> &alleles, long unsigned int &position);
    friend string json(vector<Allele*> &alleles);
    friend string json(Allele &allele, long unsigned int &position);
    friend string json(Allele* &allele);
    friend string json(Allele &allele);

public:

    AlleleType type;        // type of the allele, enumerated above
    string referenceName;   // reference name, for sanity checking
    string referenceSequence; // reference sequence or "" (in case of insertions)
    string alternateSequence; // alternate sequence or "" (in case of deletions and reference alleles)
    long unsigned int position;      // position 0-based against reference
    long unsigned int *currentReferencePosition; // pointer to the current reference position (which may be updated during the life of this allele)
    unsigned int length;    // and event length (deletion implies 0, snp implies 1, insertion >1)
    AlleleStrand strand;          // strand, true = +, false = -
    string sampleID;        // representative sample ID
    string readID;          // id of the read which the allele is drawn from
    string qualityString;   // quality string drawn from sequencer
    vector<short> baseQualities;
    short quality;          // base quality score associated with this allele
    short mapQuality;       // map quality for the originating read
    bool genotypeAllele;    // if this is an abstract 'genotype' allele
    vector<bool> indelMask; // indel mask structure, masks sites within the IDW from indels
    const bool masked(void) const;      // if the allele is masked at the *currentReferencePosition
    bool processed; // flag to mark if we've presented this allele for analysis

    // default constructor, for converting alignments into allele observations
    Allele(AlleleType t, 
                string refname, 
                long unsigned int pos, 
                long unsigned int *crefpos,
                unsigned int len, 
                string refallele, 
                string alt, 
                string sampleid,
                string readid,
                bool strnd, 
                short qual,
                string qstr, 
                short mapqual)
        : type(t)
        , referenceName(refname)
        , position(pos)
        , currentReferencePosition(crefpos)
        , length(len)
        //, referenceSequence(refallele)
        , alternateSequence(alt)
        , sampleID(sampleid)
        , readID(readid)
        , strand(strnd ? STRAND_FORWARD : STRAND_REVERSE)
        , quality((qual == -1) ? averageQuality(qstr) : qual) // passing -1 as quality triggers this calculation
        , qualityString(qstr)
        , mapQuality(mapqual) 
        , genotypeAllele(false)
        , processed(false)
    { 
        baseQualities.resize(qstr.size()); // cache qualities
        transform(qstr.begin(), qstr.end(), baseQualities.begin(), qualityChar2ShortInt);
    }

    // for constructing genotype alleles
    Allele(AlleleType t,
            string alt,
            unsigned int len,
            long unsigned int pos=0,
            bool gallele=true) 
        : type(t)
        , alternateSequence(alt)
        , length(len)
        , quality(0)
        , position(pos)
        , genotypeAllele(true)
    { }

    // I'm not sure if this explicit copy constructor is necessary
    // but it improves performance slightly
    Allele(const Allele& other) 
        : type(other.type)
        , referenceName(other.referenceName)
        , position(other.position)
        , currentReferencePosition(other.currentReferencePosition)
        , length(other.length)
        , referenceSequence(other.referenceSequence)
        , alternateSequence(other.alternateSequence)
        , sampleID(other.sampleID)
        , readID(other.readID)
        , strand(other.strand)
        , quality(other.quality)
        , qualityString(other.qualityString)
        , baseQualities(other.baseQualities)
        , mapQuality(other.mapQuality) 
        , genotypeAllele(other.genotypeAllele)
        , processed(other.processed)
    { }

    bool equivalent(Allele &a);  // heuristic 'equivalency' between two alleles, which depends on their type
    string typeStr(void); // return a string representation of the allele type, for output
    int referenceOffset(void) const;
    const short currentQuality(void) const;  // for getting the quality of a given position in multi-bp alleles
    const long double lncurrentQuality(void) const;
    bool sameSample(Allele &other);  // if the other allele has the same sample as this one
    const string base(void) const;  // the 'current' base of the allele or a string describing the allele, e.g. I10 or D2


    // overload new and delete for object recycling pool

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
    //static boost::object_pool<Allele> Pool;

};

map<string, vector<Allele*> > groupAllelesBySample(list<Allele*>& alleles);
void groupAllelesBySample(list<Allele*>& alleles, map<string, vector<Allele*> >& groups);

int allowedAlleleTypes(vector<AlleleType>& allowedEnumeratedTypes);
void filterAlleles(list<Allele*>& alleles, int allowedTypes);
void removeIndelMaskedAlleles(list<Allele*>& alleles, long unsigned int position);
int countAlleles(map<string, vector<Allele*> > sampleGroups);

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
Allele genotypeAllele(AlleleType type, string alt = "", unsigned int length = 0);

// filters... maybe move to its own file?
bool sufficientAlternateObservations(map<string, vector<Allele*> >& observations, int mincount, float minfraction);

//AlleleFreeList Allele::_freeList;

#endif
