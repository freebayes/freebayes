#ifndef __GENOTYPE_H
#define __GENOTYPE_H

#include <iostream>
#include <vector>
#include <utility> // pair
#include <algorithm>
#include <numeric>
#include <vector>
#include <iterator>
#include <math.h>
#include <numeric>
#include "Allele.h"

using namespace std;


// each genotype is a vetor of GenotypeElements, each is a count of alleles
class GenotypeElement {

    friend ostream& operator<<(ostream& out, GenotypeElement& rhs);

public:
    Allele allele;
    int count;
    GenotypeElement(const Allele& a, int c) : allele(a), count(c) { }

};


class Genotype : public vector<GenotypeElement> {

    friend ostream& operator<<(ostream& out, const pair<Allele, int>& rhs);
    friend ostream& operator<<(ostream& out, const Genotype& g);
    friend bool operator<(Genotype& a, Genotype& b);

public:
    
    int ploidy;
    vector<Allele> alleles;
    map<string, int> alleleCounts;
    bool homozygous;

    Genotype(vector<Allele>& ungroupedAlleles) {
        alleles = ungroupedAlleles;
        sort(alleles.begin(), alleles.end());
        vector<vector<Allele> > groups = groupAlleles_copy(alleles);
        for (vector<vector<Allele> >::const_iterator group = groups.begin(); group != groups.end(); ++group) {
            this->push_back(GenotypeElement(group->front(), group->size()));
            alleleCounts[group->front().currentBase] = group->size();
        }
        ploidy = getPloidy();
        homozygous = isHomozygous();
    }

    vector<Allele> uniqueAlleles(void);
    int getPloidy(void);
    int alleleFrequency(const string& base);
    int alleleFrequency(Allele& allele);
    bool containsAllele(Allele& allele);
    bool containsAllele(const string& base);
    vector<Allele> alternateAlleles(string& refbase);
    vector<int> counts(void);
    // the probability of drawing each allele out of the genotype, ordered by allele
    vector<long double> alleleProbabilities(void);
    string str(void);
    string relativeGenotype(string& refbase, string& altbase);
    bool isHomozygous(void);
    vector<int> alleleCountsInObservations(Sample& observations);
    int containedAlleleTypes(void);

};


string IUPAC(Genotype& g);
string IUPAC2GenotypeStr(string iupac);

vector<Genotype> allPossibleGenotypes(int ploidy, vector<Allele> potentialAlleles);

class SampleGenotypeProb {
public:
    string name;
    Genotype* genotype;
    long double prob;
    SampleGenotypeProb(string n, Genotype* g, long double p)
        : name(n)
        , genotype(g)
        , prob(p)
    { }
};

// a combination of genotypes for the population of samples in the analysis
class GenotypeCombo : public vector<SampleGenotypeProb> {
public:
    // GenotypeCombo::prob is equal to the sum of probs in the combo.  We
    // factor it out so that we can construct the probabilities efficiently as
    // we generate the genotype combinations
    long double prob;
    map<string, int> alleleFrequencies;  // *must* be generated at construction time

    GenotypeCombo(void) : prob(0) { }

    int numberOfAlleles(void) {
        int count = 0;
        for (GenotypeCombo::iterator g = this->begin(); g != this->end(); ++g) {
            count += g->genotype->ploidy;
        }
        return count;
    }

    void initAlleleFrequencies(void) {
        alleleFrequencies = countAlleles();
    }

    int alleleFrequency(Allele& allele) {
        map<string, int>::iterator f = alleleFrequencies.find(allele.currentBase);
        if (f == alleleFrequencies.end()) {
            return 0;
        } else {
            return f->second;
        }
    }


    void updateAlleleFrequencies(Genotype* oldGenotype, Genotype* newGenotype) {
        // remove allele frequency information for old genotype
        for (Genotype::iterator g = oldGenotype->begin(); g != oldGenotype->end(); ++g) {
            GenotypeElement& ge = *g;
            alleleFrequencies[ge.allele.currentBase] -= ge.count;
        }
        // add allele frequency information for new genotype
        for (Genotype::iterator g = newGenotype->begin(); g != newGenotype->end(); ++g) {
            GenotypeElement& ge = *g;
            alleleFrequencies[ge.allele.currentBase] += ge.count;
        }
        // remove allele frequencies which are now 0 or below
        for (map<string, int>::iterator af = alleleFrequencies.begin(); af != alleleFrequencies.end(); ++af) {
            if (af->second <= 0) {
                alleleFrequencies.erase(af);
            }
        }
    }

    map<string, int> countAlleles(void) {
        map<string, int> alleleCounts;
        for (GenotypeCombo::iterator g = this->begin(); g != this->end(); ++g) {
            for (Genotype::iterator a = g->genotype->begin(); a != g->genotype->end(); ++a) {
                map<string, int>::iterator c = alleleCounts.find(a->allele.currentBase);
                if (c != alleleCounts.end()) {
                    c->second += a->count;
                } else {
                    alleleCounts.insert(make_pair(a->allele.currentBase, a->count));
                }
            }
        }
        return alleleCounts;
    }

    map<int, int> countFrequencies(void) {
        map<int, int> frequencyCounts;
        for (map<string, int>::iterator a = alleleFrequencies.begin(); a != alleleFrequencies.end(); ++a) {
            map<int, int>::iterator c = frequencyCounts.find(a->second);
            if (c != frequencyCounts.end()) {
                c->second += 1;
            } else {
                frequencyCounts[a->second] = 1;
            }
        }
        return frequencyCounts;
    }

    vector<int> counts(void) {
        map<string, int> alleleCounts = countAlleles();
        vector<int> counts;
        for (map<string, int>::iterator a = alleleCounts.begin(); a != alleleCounts.end(); ++a) {
            counts.push_back(a->second);
        }
        return counts;
    }

    // returns true if the combination is 100% homozygous and equale
    bool isHomozygous(void) {
        GenotypeCombo::iterator g = begin();
        Genotype* genotype = g->genotype;
        if (!genotype->homozygous)
            return false;
        for (; g != end(); ++g) {
            if (g->genotype != genotype)
                return false;
        }
        return true;
    }

};


// combines a genotype combination with probabilities
class GenotypeComboResult {
public:

    GenotypeCombo* combo;
    long double priorComboProb; // derived from the below values
    long double probObsGivenGenotypes;
    long double priorProbGenotypeCombo;
    long double priorProbGenotypeComboG_Af;
    long double priorProbGenotypeComboAf;

    GenotypeComboResult(GenotypeCombo* gc,
            long double cp,
            long double pogg,
            long double ppgc,
            long double ppgcgaf,
            long double ppgcaf)
        : combo(gc)
        , priorComboProb(cp)
        , probObsGivenGenotypes(pogg)
        , priorProbGenotypeCombo(ppgc)
        , priorProbGenotypeComboG_Af(ppgcgaf)
        , priorProbGenotypeComboAf(ppgcaf)
    { }

};

class GenotypeComboResultSorter {
public:
    bool operator()(const GenotypeComboResult& gc1, const GenotypeComboResult& gc2) {
        return gc1.priorComboProb > gc2.priorComboProb;
    }
};

// a set of probabilities for a set of genotypes for a set of samples
typedef vector<pair<string, vector<pair<Genotype*, long double> > > > SampleGenotypesAndProbs;

typedef map<string, pair<Genotype*, long double> > GenotypeComboMap;

void genotypeCombo2Map(GenotypeCombo& gc, GenotypeComboMap& gcm);

void
bandedGenotypeCombinations(
    vector<GenotypeCombo>& combos,
    SampleGenotypesAndProbs& sampleGenotypes,
    int bandwidth, int banddepth);

void
bandedGenotypeCombinationsIncludingBestHomozygousCombo(
    vector<GenotypeCombo>& combos,
    SampleGenotypesAndProbs& sampleGenotypes,
    int bandwidth, int banddepth);

void
bandedGenotypeCombinationsIncludingAllHomozygousCombos(
    vector<GenotypeCombo>& combos,
    SampleGenotypesAndProbs& sampleGenotypes,
    vector<Genotype>& genotypes,
    int bandwidth, int banddepth);

vector<pair<Allele, int> > alternateAlleles(GenotypeCombo& combo, string referenceBase);

pair<int, int> alternateAndReferenceCount(vector<Allele*>& observations, string& refbase, string altbase);

ostream& operator<<(ostream& out, vector<GenotypeCombo>& combo);
ostream& operator<<(ostream& out, GenotypeCombo& g);


#endif
