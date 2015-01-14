#ifndef __GENOTYPE_H
#define __GENOTYPE_H


#include <iostream>
#include <vector>
#include <utility> // pair
#include <algorithm>
#include <numeric>
#include <vector>
#include <map>
#include <iterator>
#include <cmath>
#include <numeric>
#include <assert.h>
#include "Allele.h"
#include "Sample.h"
#include "Utility.h"
#include "Multinomial.h"
#include "CNV.h"
#include "Ewens.h"
#include "Bias.h"
#include "join.h"
#include "convert.h"

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
    long double permutationsln;  // aka, multinomialCoefficientLn(ploidy, counts())

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
        permutationsln = 0;

        if (!homozygous) {
            permutationsln = multinomialCoefficientLn(ploidy, counts());
        }

    }

    vector<Allele*> uniqueAlleles(void);
    int getPloidy(void);
    int alleleCount(const string& base);
    int alleleCount(Allele& allele);
    bool containsAllele(Allele& allele);
    bool containsAllele(const string& base);
    // returns true when the genotype is composed of a subset of the alleles
    bool matchesAlleles(vector<Allele>& alleles);
    vector<Allele> alternateAlleles(string& refbase);
    vector<string> alternateBases(string& refbase);
    vector<int> counts(void);
    // the probability of drawing each allele out of the genotype, ordered by allele
    vector<long double> alleleProbabilities(void);
    vector<long double> alleleProbabilities(Bias& observationBias);
    double alleleSamplingProb(const string& base);
    double alleleSamplingProb(Allele& allele);
    string str(void) const;
    string relativeGenotype(string& refbase, vector<Allele>& altbases);
    void relativeGenotype(vector<int>& spec, string& refbase, vector<Allele>& altbases);
    string relativeGenotype(string& refbase, string& altbase);
    void relativeGenotype(vector<int>& rg, vector<Allele>& alleles);
    bool isHeterozygous(void);
    bool isHomozygous(void);
    bool isHomozygousAlternate(void);
    bool isHomozygousReference(void);
    int containedAlleleTypes(void);
    vector<int> alleleObservationCounts(Sample& sample);
    int alleleObservationCount(Sample& sample);
    bool sampleHasSupportingObservations(Sample& sample);
    bool sampleHasSupportingObservationsForAllAlleles(Sample& sample);
    bool hasNullAllele(void);
    vector<Genotype*> nullMatchingGenotypes(vector<Genotype>& gts);

};


string IUPAC(Genotype& g);
string IUPAC2GenotypeStr(string iupac);

vector<Genotype> allPossibleGenotypes(int ploidy, vector<Allele>& potentialAlleles);

class SampleDataLikelihood {
public:
    string name;
    Genotype* genotype;
    long double prob;
    long double marginal;
    Sample* sample;
    bool hasObservations;
    int rank; // the rank of this data likelihood relative to others for the sample, 0 is best
    SampleDataLikelihood(string n, Sample* s, Genotype* g, long double p, int r)
        : name(n)
        , sample(s)
        , genotype(g)
        , prob(p)
        , rank(r)
        , marginal(0)
        , hasObservations(true)
    { }

    bool hasSupportingObservations(void) const {
        return genotype->sampleHasSupportingObservations(*sample);
    }

    int supportingObservationCount(void) const {
        return genotype->alleleObservationCount(*sample);
    }

};

class AlleleCounter {
public:
    int frequency;
    int observations;
    int forwardStrand; // supporting reads on the forward strand
    int reverseStrand; // supporting reads on the reverse strand
    int placedLeft;    // supporting reads placed to the left of the allele
    int placedRight;   // supporting reads placed to the right of the allele
    int placedStart;   // supporting reads for which the allele occurs in the first half of the read (5'-3')
    int placedEnd;     // supporting reads for which the allele occurs in the second half of the read (5'-3')
    AlleleCounter(void)
        : frequency(0)
        , observations(0)
        , forwardStrand(0)
        , reverseStrand(0)
        , placedLeft(0)
        , placedRight(0)
        , placedStart(0)
        , placedEnd(0)
    { }
};

// a combination of genotypes for the population of samples in the analysis
class GenotypeCombo : public vector<SampleDataLikelihood*> {
public:
    // GenotypeCombo::prob is equal to the sum of probs in the combo.  We
    // factor it out so that we can construct the probabilities efficiently as
    // we generate the genotype combinations
    long double probObsGivenGenotypes;  // aka data likelihood

    long double permutationsln;  // the number of perutations of unphased genotypes in the combo

    // these *must* be generated at construction time
    // for efficiency they can be updated as each genotype combo is generated
    //map<string, int> alleleCounts; // frequencies of each allele in the combo
    //map<string, pair<int, int> > alleleStrandCounts; // map from allele spec to (forword, reverse) counts
    //map<string, pair<int, int> > alleleReadPlacementCounts; // map from allele spec to (left, right) counts
    //map<string, pair<int, int> > alleleReadPositionCounts; // map from allele spec to (left, right) counts
    map<string, AlleleCounter> alleleCounters;
    map<Genotype*, int> genotypeCounts;

    GenotypeCombo(void)
        : probObsGivenGenotypes(0)
        , posteriorProb(0)
        , priorProb(0)
        , priorProbG_Af(0)
        , priorProbAf(0)
        , priorProbObservations(0)
        , permutationsln(0)
    { }

    void init(bool useObsExpectations);
    void addPriorAlleleCounts(map<string, int>& priorACs);

    // appends the other combo to this one,
    // updates the counts, and multiplies the probabilites,
    // assuming independence between the two combos
    void appendIndependentCombo(GenotypeCombo& other);

    int numberOfAlleles(void);
    vector<long double> alleleProbs(void);  // scales counts() by the total number of alleles
    int ploidy(void); // the number of copies of the locus in this combination
    int alleleCount(Allele& allele);
    int alleleCount(const string& allele);
    long double alleleFrequency(Allele& allele);
    long double alleleFrequency(const string& allele);
    long double genotypeFrequency(Genotype* genotype);
    void updateCachedCounts(Sample* sample, Genotype* oldGenotype, Genotype* newGenotype, bool useObsExpectations);
    map<string, int> countAlleles(void);
    map<int, int> countFrequencies(void);
    int hetCount(void);
    vector<int> counts(void); // the counts of frequencies of the alleles in the genotype combo
    vector<int> observationCounts(void); // the counts of observations of the alleles (in sorted order)
    int observationTotal(void);
    vector<string> alleles(void);  // the string representations of alleles in the genotype combo
    bool isHomozygous(void); // returns true if the combination is 100% homozygous across all individuals
                             // e.g. if there is no variation

    // posterior

    long double posteriorProb; // p(genotype combo) * p(observations | genotype combo)

    // priors

    long double priorProb; // p(genotype combo) = p(genotype combo | allele frequency) * p(allele frequency) * p(observations)
    long double priorProbG_Af; // p(genotype combo | allele frequency)
    long double priorProbAf; // p(allele frequency)
    long double priorProbObservations; // p(observations)
    long double priorProbGenotypesGivenHWE;

    //GenotypeCombo* combo,
    void calculatePosteriorProbability(
        long double theta,
        bool pooled,
        bool ewensPriors,
        bool permute,
        bool hwePriors,
        bool obsBinomialPriors,
        bool alleleBalancePriors,
        long double diffusionPriorScalarln);

    long double probabilityGivenAlleleFrequencyln(bool permute);

    long double hweExpectedFrequencyln(Genotype* genotype);
    long double hweProbGenotypeFrequencyln(Genotype* genotype);
    long double hweComboProb(void);

};

struct GenotypeComboResultSorter {
    bool operator()(const GenotypeCombo& gc1, const GenotypeCombo& gc2) {
        if (gc1.posteriorProb == gc2.posteriorProb) {
            return gc1 > gc2;
        } else {
            return gc1.posteriorProb > gc2.posteriorProb;
        }
    }
};

// for comparing GenotypeCombos which are empty
struct GenotypeComboResultEqual {
    bool operator()(const GenotypeCombo& gc1, const GenotypeCombo& gc2) {
	return gc1.posteriorProb == gc2.posteriorProb;
    }
};

// for sorting data likelihoods
struct SampleDataLikelihoodCompare {
    bool operator()(const SampleDataLikelihood& a,
            const SampleDataLikelihood& b) {
        return a.prob > b.prob;
    }
};

struct SampleMarginalCompare {
    bool operator()(const SampleDataLikelihood& a,
            const SampleDataLikelihood& b) {
        return a.marginal > b.marginal;
    }
};

struct SampleLikelihoodCompare {
    bool operator()(const SampleDataLikelihood& a,
            const SampleDataLikelihood& b) {
        return (a.marginal + a.prob) > (b.marginal + b.prob);
    }
};

struct SampleMarginalAndObsCompare {
    bool operator()(const SampleDataLikelihood& a,
            const SampleDataLikelihood& b) {
        int aObsCount = a.supportingObservationCount();
        int bObsCount = b.supportingObservationCount();
        if (aObsCount != bObsCount) {
            if (aObsCount == 0) {
                return false;
            } else if (bObsCount == 0) {
                return true;
            }
        }
        return (a.marginal + a.prob) > (b.marginal + b.prob);
    }
};

// a set of probabilities for a set of genotypes for a set of samples
typedef vector<vector<SampleDataLikelihood> > SampleDataLikelihoods;

void sortSampleDataLikelihoods(vector<SampleDataLikelihood>& likelihoods);
bool sortSampleDataLikelihoodsByMarginals(vector<SampleDataLikelihood>& likelihoods);
bool sortSampleDataLikelihoodsByMarginals(SampleDataLikelihoods& samplesLikelihoods);
bool sortSampleDataLikelihoodsByMarginalsAndObs(SampleDataLikelihoods& samplesLikelihoods);
bool sortSampleDataLikelihoodsScaledByMarginals(vector<SampleDataLikelihood>& likelihoods);
bool sortSampleDataLikelihoodsScaledByMarginals(SampleDataLikelihoods& samplesLikelihoods);

typedef map<string, SampleDataLikelihood*> GenotypeComboMap;

void genotypeCombo2Map(GenotypeCombo& gc, GenotypeComboMap& gcm);

void
orderedGenotypeCombo(
    GenotypeCombo& combo,
    GenotypeCombo& orderedCombo,
    SampleDataLikelihoods& sampleDataLikelihoods,
    long double theta,
    bool pooled,
    bool ewensPriors,
    bool permute,
    bool hwePriors,
    bool binomialObsPriors,
    bool alleleBalancePriors,
    long double diffusionPriorScalar);

void
makeComboByDatalLikelihoodRank(
    GenotypeCombo& combo,
    vector<int>& initialPosition,
    SampleDataLikelihoods& variantSampleDataLikelihoods,
    SampleDataLikelihoods& invariantSampleDataLikelihoods,
    map<string, int>& priorACs,
    long double theta,
    bool pooled,
    bool ewensPriors,
    bool permute,
    bool hwePriors,
    bool binomialObsPriors,
    bool alleleBalancePriors,
    long double diffusionPriorScalar);

void
dataLikelihoodMaxGenotypeCombo(
    GenotypeCombo& combo,
    SampleDataLikelihoods& sampleDataLikelihoods,
    long double theta,
    bool pooled,
    bool ewensPriors,
    bool permute,
    bool hwePriors,
    bool binomialObsPriors,
    bool alleleBalancePriors,
    long double diffusionPriorScalar);

bool
bandedGenotypeCombinations(
    list<GenotypeCombo>& combos,
    GenotypeCombo& comboKing,
    SampleDataLikelihoods& variantDataLikelihoods,
    SampleDataLikelihoods& invariantDataLikelihoods,
    Samples& samples,
    map<string, int>& priorACs,
    int bandwidth, int banddepth,
    long double theta,
    bool pooled,
    bool ewensPriors,
    bool permute,
    bool hwePriors,
    bool binomialObsPriors,
    bool alleleBalancePriors,
    long double diffusionPriorScalar);

void
allLocalGenotypeCombinations(
    list<GenotypeCombo>& combos,
    GenotypeCombo& comboKing,
    SampleDataLikelihoods& sampleDataLikelihoods,
    Samples& samples,
    map<string, int>& priorACs,
    long double theta,
    bool pooled,
    bool ewensPriors,
    bool permute,
    bool hwePriors,
    bool binomialObsPriors,
    bool alleleBalancePriors,
    long double diffusionPriorScalar,
    bool keepCombos);

void
convergentGenotypeComboSearch(
    list<GenotypeCombo>& combos,
    GenotypeCombo& comboKing,
    SampleDataLikelihoods& sampleDataLikelihoods,
    SampleDataLikelihoods& variantDataLikelihoods,
    SampleDataLikelihoods& invariantDataLikelihoods,
    Samples& samples,
    vector<Allele>& genotypeAlleles,
    map<string, int>& priorACs,
    int bandwidth, int banddepth,
    long double theta,
    bool pooled,
    bool ewensPriors,
    bool permute,
    bool hwePriors,
    bool binomialObsPriors,
    bool alleleBalancePriors,
    long double diffusionPriorScalar,
    int maxiterations,
    int& totaliterations,
    bool addHomozygousCombos);

void
addAllHomozygousCombos(
    list<GenotypeCombo>& combos,
    SampleDataLikelihoods& sampleDataLikelihoods,
    SampleDataLikelihoods& variantSampleDataLikelihoods,
    SampleDataLikelihoods& invariantSampleDataLikelihoods,
    Samples& samples,
    vector<Allele>& genotypeAlleles,
    long double theta,
    bool pooled,
    bool ewensPriors,
    bool permute,
    bool hwePriors,
    bool binomialObsPriors,
    bool alleleBalancePriors,
    long double diffusionPriorScalar);


vector<pair<Allele, int> > alternateAlleles(GenotypeCombo& combo, string referenceBase);

pair<int, int> alternateAndReferenceCount(vector<Allele*>& observations, string& refbase, string altbase);

ostream& operator<<(ostream& out, list<GenotypeCombo>& combo);
ostream& operator<<(ostream& out, GenotypeCombo& g);

map<int, vector<Genotype> > getGenotypesByPloidy(vector<int>& ploidies, vector<Allele>& genotypeAlleles);

void combinePopulationCombos(list<GenotypeCombo>& genotypeCombos,
                             map<string, list<GenotypeCombo> >& genotypeCombosByPopulation);

#endif
