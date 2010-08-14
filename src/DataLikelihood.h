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
#include "Dirichlet.h"


long double
probObservedAllelesGivenGenotype(
        vector<Allele*>& observedAlleles,
        Genotype& genotype);

long double
approxProbObservedAllelesGivenGenotype(
        vector<Allele*>& observedAlleles,
        Genotype& genotype);

vector<pair<Genotype*, long double> >
probObservedAllelesGivenGenotypes(
        vector<Allele*>& observedAlleles,
        vector<Genotype>& genotype);


vector<pair<Genotype*, long double> >
exactProbObservedAllelesGivenGenotypes(
        vector<Allele*>& observedAlleles,
        vector<Genotype>& genotypes);

vector<pair<Genotype*, long double> >
approxProbObservedAllelesGivenGenotypes(
        vector<Allele*>& observedAlleles,
        vector<Genotype>& genotypes);

long double
likelihoodGivenTrueAlleles(
        const vector<Allele*>& observedAlleles,
        const vector<Allele*>& trueAlleles);

#endif
