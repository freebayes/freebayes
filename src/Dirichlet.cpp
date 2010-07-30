#include "Dirichlet.h"
#include "Sum.h"
#include "Product.h"


long double dirichlet(const vector<long double>& probs, 
        const vector<int>& obs, 
        long double s) {

    vector<long double> alphas;
    for (vector<int>::const_iterator o = obs.begin(); o != obs.end(); ++o)
        alphas.push_back(*o + 1 * s);

    vector<long double> obsProbs;
    vector<long double>::const_iterator a = alphas.begin();
    vector<long double>::const_iterator p = probs.begin();
    for (; p != probs.end() && a != alphas.end(); ++p, ++a) {
        obsProbs.push_back(pow(*p, *a - 1));
    }

    return 1.0 / beta(alphas) * product(obsProbs);

}

long double dirichletMaximumLikelihoodRatio(const vector<long double>& probs,
        const vector<int>& obs, 
        long double s) {
    int maximizingObs = obs.size() / sum(obs);
    vector<int> m(obs.size(), maximizingObs);
    return dirichlet(probs, obs, s) / dirichlet(probs, m, s);
}

