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
#include "AlleleParser.h"
#include "ResultData.h"

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

void
calculateSampleDataLikelihoods(
    Samples& samples,
    Results& results,
    AlleleParser* parser,
    map<int, vector<Genotype> >& genotypesByPloidy,
    Parameters& parameters,
    bool usingNull,
    Bias& observationBias,
    vector<Allele>& genotypeAlleles,
    Contamination& contaminationEstimates,
    map<string, double>& estimatedAlleleFrequencies,
    map<string, vector<vector<SampleDataLikelihood> > >& sampleDataLikelihoodsByPopulation,
    map<string, vector<vector<SampleDataLikelihood> > >& variantSampleDataLikelihoodsByPopulation,
    map<string, vector<vector<SampleDataLikelihood> > >& invariantSampleDataLikelihoodsByPopulation);

#endif
