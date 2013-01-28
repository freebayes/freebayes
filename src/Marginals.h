#ifndef __MARGINALS_H
#define __MARGINALS_H
#include <vector>
#include <map>
#include "Genotype.h"
#include "ResultData.h"
#include "Utility.h"

using namespace std;

//void marginalGenotypeLikelihoods(list<GenotypeCombo>& genotypeCombos, Results& results);
long double marginalGenotypeLikelihoods(list<GenotypeCombo>& genotypeCombos, SampleDataLikelihoods& likelihoods);
void bestMarginalGenotypeCombo(GenotypeCombo& combo,
        Results& results,
        SampleDataLikelihoods& samples,
        long double theta,
        bool pooled,
        bool permute,
        bool hwePriors,
        bool binomialObsPriors,
        bool alleleBalancePriors,
        long double diffusionPriorScalar);

long double balancedMarginalGenotypeLikelihoods(list<GenotypeCombo>& genotypeCombos, SampleDataLikelihoods& likelihoods);

#endif
