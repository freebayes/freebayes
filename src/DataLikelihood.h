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
#include "Sample.h"
#include "Genotype.h"
#include "Utility.h"
#include "Multinomial.h"
#include "Dirichlet.h"

using namespace std;

long double
probObservedAllelesGivenGenotype(
        Sample& sample,
        Genotype& genotype);

long double
bamBayesProbObservedAllelesGivenGenotype(
        Sample& sample,
        Genotype& genotype,
        long double dependenceFactor);

long double
approxProbObservedAllelesGivenGenotype(
        Sample& sample,
        Genotype& genotype,
        long double dependenceFactor);

vector<pair<Genotype*, long double> >
probObservedAllelesGivenGenotypes(
        Sample& sample,
        vector<Genotype>& genotype,
        long double dependenceFactor);

vector<pair<Genotype*, long double> >
bamBayesProbObservedAllelesGivenGenotypes(
        Sample& sample,
        vector<Genotype>& genotype,
        long double dependenceFactor);

long double
likelihoodGivenTrueAlleles(
        const vector<Allele*>& observedAlleles,
        const vector<Allele*>& trueAlleles);

#endif
