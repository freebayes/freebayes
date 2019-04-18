#ifndef __RESULT_H
#define __RESULT_H

#include <vector>
#include <string>
#include <algorithm>
#include <utility>
#include "Genotype.h"

using namespace std;

class Result : public vector<SampleDataLikelihood> {

public:

    string name;
    Sample* observations;

    void sortDataLikelihoods(void);

    //pair<Genotype*, long double> bestMarginalGenotype(void);

};

#endif
