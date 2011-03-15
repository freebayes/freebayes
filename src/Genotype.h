#ifndef __GENOTYPE_H
#define __GENOTYPE_H

#include <iostream>
#include <vector>
#include <utility> // pair
#include <algorithm>
#include <numeric>
#include <vector>
#include <iterator>
#include <cmath>
#include <numeric>
#include "Allele.h"
#include "Sample.h"

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
    int containedAlleleTypes(void);
    vector<int> alleleObservationCounts(Sample& sample);

};


string IUPAC(Genotype& g);
string IUPAC2GenotypeStr(string iupac);

vector<Genotype> allPossibleGenotypes(int ploidy, vector<Allele> potentialAlleles);

class SampleGenotypeProb {
public:
    string sampleName;
    Genotype* genotype;
    long double prob;
    Sample* sample;
    SampleGenotypeProb(string n, Sample* s, Genotype* g, long double p)
        : sampleName(n)
        , sample(s)
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

    // these *must* be generated at construction time
    // for efficiency they can be updated as each genotype combo is generated
    map<string, int> alleleFrequencies; // frequencies of each allele in the combo
    map<string, pair<int, int> > alleleStrandCounts; // map from allele spec to (forword, reverse) counts
    map<string, pair<int, int> > alleleReadPlacementCounts; // map from allele spec to (left, right) counts
    map<string, pair<int, int> > alleleHetRefAltCounts; // map from allele spec to (ref, alt) counts

    GenotypeCombo(void) : prob(0) { }

    void init(void);

    int numberOfAlleles(void);
    void initAlleleFrequencies(void);
    int alleleFrequency(Allele& allele);
    void updateCachedCounts(Sample* sample, Genotype* oldGenotype, Genotype* newGenotype);
    map<string, int> countAlleles(void);
    map<int, int> countFrequencies(void);
    vector<int> counts(void); // the counts of frequencies of the alleles in the genotype combo
    vector<string> alleles(void);  // the string representations of alleles in the genotype combo
    bool isHomozygous(void); // returns true if the combination is 100% homozygous across all individuals
                             // e.g. if there is no variation

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
    long double priorProbBinomialObservations;

    GenotypeComboResult(GenotypeCombo* gc,
            long double cp,
            long double pogg,
            long double ppgc,
            long double ppgcgaf,
            long double ppgcaf,
            long double ppbo)
        : combo(gc)
        , priorComboProb(cp)
        , probObsGivenGenotypes(pogg)
        , priorProbGenotypeCombo(ppgc)
        , priorProbGenotypeComboG_Af(ppgcgaf)
        , priorProbGenotypeComboAf(ppgcaf)
        , priorProbBinomialObservations(ppbo)
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
    Samples& samples,
    int bandwidth, int banddepth,
    float logStepMax);

void
bandedGenotypeCombinationsIncludingBestHomozygousCombo(
    vector<GenotypeCombo>& combos,
    SampleGenotypesAndProbs& sampleGenotypes,
    Samples& samples,
    int bandwidth, int banddepth,
    float logStepMax);

void
bandedGenotypeCombinationsIncludingAllHomozygousCombos(
    vector<GenotypeCombo>& combos,
    SampleGenotypesAndProbs& sampleGenotypes,
    Samples& samples,
    map<int, vector<Genotype> >& genotypesByPloidy,
    vector<Allele>& genotypeAlleles,
    int bandwidth, int banddepth,
    float logStepMax);

vector<pair<Allele, int> > alternateAlleles(GenotypeCombo& combo, string referenceBase);

pair<int, int> alternateAndReferenceCount(vector<Allele*>& observations, string& refbase, string altbase);

ostream& operator<<(ostream& out, vector<GenotypeCombo>& combo);
ostream& operator<<(ostream& out, GenotypeCombo& g);


#endif
