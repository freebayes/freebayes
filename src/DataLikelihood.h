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

using namespace std;

long double
probObservedAllelesGivenGenotype(
        Sample& sample,
        Genotype& genotype,
        long double dependenceFactor,
        bool useMapQ);

vector<pair<Genotype*, long double> >
probObservedAllelesGivenGenotypes(
        Sample& sample,
        vector<Genotype*>& genotypes,
        long double dependenceFactor,
        bool useMapQ);

#endif
