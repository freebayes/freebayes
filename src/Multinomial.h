#ifndef __MULTINOMIAL_H
#define __MULTINOMIAL_H

#include "Utility.h"
#include <vector>

long double multinomialSamplingProb(const vector<long double>& probs, const vector<int>& obs);
long double multinomialSamplingProbLn(const vector<long double>& probs, const vector<int>& obs);
long double multinomialCoefficientLn(int n, const vector<int>& counts);

long double samplingProbLn(const vector<long double>& probs, const vector<int>& obs);

#endif
