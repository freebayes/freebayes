#ifndef __DATALIKELIHOOD_H
#define __DATALIKELIHOOD_H

#include <iostream>
#include <vector>
#include <utility> // pair
#include <algorithm>
#include <numeric>
#include <vector>
#include <iterator>
#include "Allele.h"
#include "Genotype.h"
#include "Utility.h"
#include "Multinomial.h"


long double
probObservedAllelesGivenGenotype(
        vector<Allele*>& observedAlleles,
        Genotype& genotype);

vector<pair<Genotype, long double> >
probObservedAllelesGivenGenotypes(
        vector<Allele*>& observedAlleles,
        vector<Genotype>& genotype);

long double
likelihoodGivenTrueAlleles(
        const vector<Allele*>& observedAlleles,
        const vector<Allele*>& trueAlleles);

#endif
