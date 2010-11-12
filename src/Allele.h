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

using namespace std;

class Allele;
class Sample;

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
    long double position;      // position 0-based against reference
    long double* currentReferencePosition; // pointer to the current reference position (which may be updated during the life of this allele)
    char* currentReferenceBase;  // pointer to current reference base
    unsigned int length;    // and event length (deletion implies 0, snp implies 1, insertion >1)
    AlleleStrand strand;          // strand, true = +, false = -
    string sampleID;        // representative sample ID
    string readID;          // id of the read which the allele is drawn from
    vector<short> baseQualities;
    long double quality;          // base quality score associated with this allele, updated every position in the case of reference alleles
    long double lnquality;  // log version of above
    string currentBase;       // current base, meant to be updated every position
    short mapQuality;       // map quality for the originating read
    bool genotypeAllele;    // if this is an abstract 'genotype' allele
    vector<bool> indelMask; // indel mask structure, masks sites within the IDW from indels
    const bool masked(void) const;      // if the allele is masked at the *currentReferencePosition
    bool processed; // flag to mark if we've presented this allele for analysis

    // default constructor, for converting alignments into allele observations
    Allele(AlleleType t, 
                string refname, 
                long double pos, 
                long double* crefpos,
                char* crefbase,
                unsigned int len, 
                string refallele, 
                string alt, 
                string sampleid,
                string readid,
                bool strnd, 
                long double qual,
                string qstr, 
                short mapqual)
        : type(t)
        , referenceName(refname)
        , position(pos)
        , currentReferencePosition(crefpos)
        , currentReferenceBase(crefbase)
        , length(len)
        , currentBase(alt)
        , alternateSequence(alt)
        , sampleID(sampleid)
        , readID(readid)
        , strand(strnd ? STRAND_FORWARD : STRAND_REVERSE)
        , quality((qual == -1) ? averageQuality(qstr) : qual) // passing -1 as quality triggers this calculation
        , lnquality(log((qual == -1) ? averageQuality(qstr) : qual))
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
            long double pos=0,
            bool gallele=true) 
        : type(t)
        , alternateSequence(alt)
        , currentBase(alt)
        , length(len)
        , quality(0)
        , lnquality(1)
        , position(pos)
        , genotypeAllele(true)
    { }

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
    int referenceOffset(void) const;
    const short currentQuality(void) const;  // for getting the quality of a given position in multi-bp alleles
    const long double lncurrentQuality(void) const;
    bool sameSample(Allele &other);  // if the other allele has the same sample as this one
    void update(void); // for reference alleles, updates currentBase and quality
    // TODO update this to reflect different insertions (e.g. IATGC instead of I4)
    const string base(void) const;  // the 'current' base of the allele or a string describing the allele, e.g. I10 or D2
                                    //  this is used to update cached data in the allele prior to presenting the allele for analysis
                                    //  for the current base, just use allele.currentBase

    string json(void);

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


// sample tracking and allele sorting
class Sample : public map<string, vector<Allele*> > {

public:

    // the number of observations for this allele
    int observationCount(Allele& allele) {
        return observationCount(allele.currentBase);
    }

    // the number of observations for this base
    int observationCount(const string& base) {
        Sample::iterator g = find(base);
        if (g != end())
            return g->second.size();
        else
            return 0;
    }

    // the total number of observations
    int observationCount(void) {
        int count = 0;
        for (Sample::iterator g = begin(); g != end(); ++g) {
            count += g->second.size();
        }
        return count;
    }

    // puts alleles into the right bins if they have changed their base (as
    // occurs in the case of reference alleles)
    void sortAlleles(void) {
        for (Sample::iterator g = begin(); g != end(); ++g) {
            const string& groupBase = g->first;
            vector<Allele*>& alleles = g->second;
            for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                const string& base = (*a)->currentBase;
                if (base != groupBase) {
                    Sample::iterator g = find(base);
                    if (g != end()) {
                        g->second.push_back(*a);
                    } else {
                        vector<Allele*> alleles;
                        alleles.push_back(*a);
                        insert(begin(), make_pair(base, alleles));
                    }
                    *a = NULL;
                }
            }
            alleles.erase(remove(alleles.begin(), alleles.end(), (Allele*)NULL), alleles.end());
        }
    }

    pair<pair<int, int>, pair<int, int> >
    baseCount(string refbase, string altbase) {

        int forwardRef = 0;
        int reverseRef = 0;
        int forwardAlt = 0;
        int reverseAlt = 0;

        for (Sample::iterator s = begin(); s != end(); ++s) {
            vector<Allele*>& alleles = s->second;
            for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                string base = (*a)->currentBase;
                AlleleStrand strand = (*a)->strand;
                if (base == refbase) {
                    if (strand == STRAND_FORWARD)
                        ++forwardRef;
                    else if (strand == STRAND_REVERSE)
                        ++reverseRef;
                } else if (base == altbase) {
                    if (strand == STRAND_FORWARD)
                        ++forwardAlt;
                    else if (strand == STRAND_REVERSE)
                        ++reverseAlt;
                }
            }
        }

        return make_pair(make_pair(forwardRef, forwardAlt), make_pair(reverseRef, reverseAlt));

    }


    int baseCount(string base, AlleleStrand strand) {

        int count = 0;
        for (Sample::iterator g = begin(); g != end(); ++g) {
            vector<Allele*>& alleles = g->second;
            for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                if ((*a)->currentBase == base && (*a)->strand == strand)
                    ++count;
            }
        }
        return count;

    }


    string json(void) {
        stringstream out;
        out << "[";
        bool first = true;
        for (map<string, vector<Allele*> >::iterator g = this->begin(); g != this->end(); ++g) {
            vector<Allele*>& alleles = g->second;
            for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                if (!first) { out << ","; } else { first = false; }
                out << (*a)->json();
            }
        }
        out << "]";
        return out.str();
    }

};

class Samples : public map<string, Sample> { };


void updateAllelesCachedData(vector<Allele*>& alleles);

map<string, vector<Allele*> > groupAllelesBySample(list<Allele*>& alleles);
void groupAllelesBySample(list<Allele*>& alleles, map<string, vector<Allele*> >& groups);

int allowedAlleleTypes(vector<AlleleType>& allowedEnumeratedTypes);
void filterAlleles(list<Allele*>& alleles, int allowedTypes);
void removeIndelMaskedAlleles(list<Allele*>& alleles, long double position);
int countAlleles(Samples& samples);
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

// using this one...
void groupAlleles(Samples& samples, map<string, vector<Allele*> >& alleleGroups);

void groupAlleles(map<string, vector<Allele*> >& sampleGroups, map<string, vector<Allele*> >& alleleGroups);

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
Allele genotypeAllele(AlleleType type, string alt = "", unsigned int length = 0);

// filters... maybe move to its own file?
bool sufficientAlternateObservations(Samples& observations, int mincount, float minfraction);

//AlleleFreeList Allele::_freeList;

#endif
