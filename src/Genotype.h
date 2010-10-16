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

    Genotype(vector<Allele>& ungroupedAlleles) {
        alleles = ungroupedAlleles;
        sort(alleles.begin(), alleles.end());
        vector<vector<Allele> > groups = groupAlleles_copy(alleles);
        for (vector<vector<Allele> >::const_iterator group = groups.begin(); group != groups.end(); ++group) {
            this->push_back(GenotypeElement(group->front(), group->size()));
        }
        ploidy = getPloidy();
    }

    vector<Allele> uniqueAlleles(void);
    int getPloidy(void);
    int alleleFrequency(Allele& allele);
    bool containsAllele(Allele& allele);
    bool containsAlleleOtherThan(string& base);
    vector<Allele> alternateAlleles(string& refbase);
    int alleleCount(string& base);
    // the probability of drawing each allele out of the genotype, ordered by allele
    vector<long double> alleleProbabilities(void);
    string str(void);
    string relativeGenotype(string& refbase, string& altbase);
    bool homozygous(void);
    vector<int> alleleCountsInObservations(vector<Allele*> observations);

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

    GenotypeCombo(void) : prob(0) { }

    int numberOfAlleles(void) {
        int count = 0;
        for (GenotypeCombo::iterator g = this->begin(); g != this->end(); ++g) {
            count += g->genotype->ploidy;
        }
        return count;
    }

    map<Allele, int> countAlleles(void) {
        map<Allele, int> alleleCounts;
        for (GenotypeCombo::iterator g = this->begin(); g != this->end(); ++g) {
            for (Genotype::iterator a = g->genotype->begin(); a != g->genotype->end(); ++a) {
                map<Allele, int>::iterator c = alleleCounts.find(a->allele);
                if (c != alleleCounts.end()) {
                    c->second += a->count;
                } else {
                    alleleCounts.insert(make_pair(a->allele, a->count));
                }
            }
        }
        return alleleCounts;
    }

    map<int, int> countFrequencies(void) {
        map<int, int> frequencyCounts;
        map<Allele, int> alleles = countAlleles();
        for (map<Allele, int>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
            map<int, int>::iterator c = frequencyCounts.find(a->second);
            if (c != frequencyCounts.end()) {
                c->second += 1;
            } else {
                frequencyCounts[a->second] = 1;
            }
        }
        return frequencyCounts;
    }

};


// combines a genotype combination with probabilities
class GenotypeComboResult {
public:

    GenotypeCombo combo;
    long double priorComboProb; // derived from the below values
    long double probObsGivenGenotypes;
    long double priorProbGenotypeCombo;
    long double priorProbGenotypeComboG_Af;
    long double priorProbGenotypeComboAf;

    GenotypeComboResult(GenotypeCombo& gc,
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

bool isHomozygousCombo(GenotypeCombo& combo);

vector<pair<Allele, int> > alternateAlleles(GenotypeCombo& combo, string referenceBase);

pair<int, int> alternateAndReferenceCount(vector<Allele*>& observations, string& refbase, string altbase);

ostream& operator<<(ostream& out, vector<GenotypeCombo>& combo);
ostream& operator<<(ostream& out, GenotypeCombo& g);


#endif
