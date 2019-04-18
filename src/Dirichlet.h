#include <vector>
#include "Utility.h"
#include "Sum.h"

long double dirichletMaximumLikelihoodRatio(const vector<long double>& probs, const vector<int>& obs, long double s = (long double) 1.0);
long double dirichlet(const vector<long double>& probs, const vector<int>& obs, long double s = (long double) 1.0);
long double dirichletMaximumLikelihoodRatioln(const vector<long double>& probs, const vector<int>& obs, long double s = (long double) 1.0);
long double dirichletln(const vector<long double>& probs, const vector<int>& obs, long double s = (long double) 1.0);
