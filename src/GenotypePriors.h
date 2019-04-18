#ifndef __GENOTYPE_PRIORS_H
#define __GENOTYPE_PRIORS_H

#include <map>
#include "Allele.h"
#include "Genotype.h"
#include "Multinomial.h"
#include "CNV.h"
#include "Utility.h"

using namespace std;

map<Allele, int> countAlleles(vector<Genotype*>& genotypeCombo);
map<int, int> countFrequencies(vector<Genotype*>& genotypeCombo);
long double alleleFrequencyProbability(const map<int, int>& alleleFrequencyCounts, long double theta);
long double alleleFrequencyProbabilityln(const map<int, int>& alleleFrequencyCounts, long double theta);
long double __alleleFrequencyProbabilityln(const map<int, int>& alleleFrequencyCounts, long double theta);
long double probabilityGenotypeComboGivenAlleleFrequencyln(GenotypeCombo& genotypeCombo, Allele& allele);

class AlleleFrequencyProbabilityCache : public map<map<int, int>, long double> {
public:
    long double alleleFrequencyProbabilityln(const map<int, int>& counts, long double theta) {
        map<map<int, int>, long double>::iterator p = find(counts);
        if (p == end()) {
            long double pln = __alleleFrequencyProbabilityln(counts, theta);
            insert(make_pair(counts, pln));
            return pln;
        } else {
            return p->second;
        }
    }
};

GenotypeComboResult
genotypeCombinationsPriorProbability(
        GenotypeCombo* combo,
        Allele& refAllele,
        long double theta,
        bool pooled,
        bool obsBinomialPriors,
        bool alleleBalancePriors,
        long double diffusionPriorScalarln);

void genotypeCombinationsPriorProbability(
        vector<GenotypeComboResult>& genotypeComboProbs,
        vector<GenotypeCombo>& bandedCombos,
        Allele& refAllele,
        long double theta,
        bool pooled,
        bool obsBinomialPriors,
        bool alleleBalancePriors,
        long double diffusionPriorScalarln);

#endif
