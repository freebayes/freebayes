#ifndef __GENOTYPE_PRIORS_H
#define __GENOTYPE_PRIORS_H

#include <map>
#include "Allele.h"
#include "Genotype.h"

using namespace std;

map<Allele, int> countAlleles(vector<Genotype*>& genotypeCombo);
map<int, int> countFrequencies(vector<Genotype*>& genotypeCombo);
long double alleleFrequencyProbability(map<int, int> alleleFrequencyCounts, long double theta);
long double alleleFrequencyProbabilityln(map<int, int> alleleFrequencyCounts, long double theta);
int alleleFrequency(GenotypeCombo& genotypeCombo, Allele& allele);
long double probabilityDiploidGenotypeComboGivenAlleleFrequencyln(GenotypeCombo& genotypeCombo, Allele& allele);

#endif
