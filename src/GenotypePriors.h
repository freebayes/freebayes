#ifndef __GENOTYPE_PRIORS_H
#define __GENOTYPE_PRIORS_H

#include <map>
#include "Allele.h"
#include "Genotype.h"
#include "Multinomial.h"
#include "CNV.h"

using namespace std;

map<Allele, int> countAlleles(vector<Genotype*>& genotypeCombo);
map<int, int> countFrequencies(vector<Genotype*>& genotypeCombo);
long double alleleFrequencyProbability(map<int, int> alleleFrequencyCounts, long double theta);
long double alleleFrequencyProbabilityln(map<int, int> alleleFrequencyCounts, long double theta);
long double probabilityGenotypeComboGivenAlleleFrequencyln(GenotypeCombo& genotypeCombo, Allele& allele);

GenotypeComboResult
genotypeCombinationsPriorProbability(
        GenotypeCombo* combo,
        Allele& refAllele,
        long double theta,
        bool pooled,
        bool useBinomialProbs,
        long double diffusionPriorScalarln);

void genotypeCombinationsPriorProbability(
        vector<GenotypeComboResult>& genotypeComboProbs,
        vector<GenotypeCombo>& bandedCombos,
        Allele& refAllele,
        long double theta,
        bool pooled,
        bool useBinomialProbs,
        long double diffusionPriorScalarln);

#endif
