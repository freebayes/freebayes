#ifndef __DATALIKELIHOOD_H
#define __DATALIKELIHOOD_H

#include <iostream>
#include <vector>
#include <utility> // pair
#include <algorithm>
#include <numeric>
#include <vector>
#include <iterator>
#include <cmath>
#include "Allele.h"
#include "Sample.h"
#include "Genotype.h"
#include "Utility.h"
#include "Multinomial.h"
#include "Dirichlet.h"
#include "Bias.h"
#include "Contamination.h"

using namespace std;

long double
probObservedAllelesGivenGenotype(
        Sample& sample,
        Genotype& genotype,
        double dependenceFactor,
        bool useMapQ,
        Bias& observationBias,
        bool standardGLs,
        vector<Allele>& genotypeAlleles,
        vector<Allele>& uniqueAllelesInGenotype,
        Contamination& contaminations,
        map<string, double>& freqs);

vector<pair<Genotype*, long double> >
probObservedAllelesGivenGenotypes(
        Sample& sample,
        vector<Genotype*>& genotypes,
        double dependenceFactor,
        bool useMapQ,
        Bias& observationBias,
        bool standardGLs,
        vector<Allele>& genotypeAlleles,
        Contamination& contaminations,
        map<string, double>& freqs);

#endif
