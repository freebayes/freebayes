#ifndef __MARGINALS_H
#define __MARGINALS_H
#include <vector>
#include <map>
#include "Genotype.h"
#include "ResultData.h"

void marginalGenotypeLikelihoods(long double posteriorNormalizer, list<GenotypeCombo>& genotypeCombos, Results& results);
void bestMarginalGenotypeCombo(GenotypeCombo& combo,
        Results& results,
        SampleDataLikelihoods& samples,
        long double theta,
        bool pooled,
        bool binomialObsPriors,
        bool alleleBalancePriors,
        long double diffusionPriorScalar);

#endif
