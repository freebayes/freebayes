#ifndef __MULTINOMIAL_H
#define __MULTINOMIAL_H

#include "Utility.h"
#include <vector>

long double multinomialSamplingProb(vector<long double> probs, vector<int> obs);
long double multinomialSamplingProbLn(vector<long double> probs, vector<int> obs);
long double multinomialCoefficientLn(int n, vector<int> counts);

#endif
