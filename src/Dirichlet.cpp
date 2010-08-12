#include "Dirichlet.h"
#include "Sum.h"
#include "Product.h"
#include <iostream>


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
    long double maximizingObs = obs.size() / sum(obs);
    vector<int> m(obs.size(), maximizingObs);
    return dirichlet(probs, obs, s) / dirichlet(probs, m, s);
}


// XXX the logspace versions are broken

long double dirichletln(const vector<long double>& probs, 
        const vector<int>& obs, 
        long double s) {

    vector<long double> alphas;
    for (vector<int>::const_iterator o = obs.begin(); o != obs.end(); ++o)
        alphas.push_back(*o + 1 * s);

    vector<long double> obsProbs;
    vector<long double>::const_iterator a = alphas.begin();
    vector<long double>::const_iterator p = probs.begin();
    for (; p != probs.end() && a != alphas.end(); ++p, ++a) {
        obsProbs.push_back(powln(log(*p), *a - 1));
    }
    
    return log(1.0) - (betaln(alphas) + sum(obsProbs));

}

long double dirichletMaximumLikelihoodRatioln(const vector<long double>& probs,
        const vector<int>& obs, 
        long double s) {
    long double maximizingObs = (long double) obs.size() / (long double) sum(obs);
    vector<int> m(obs.size(), maximizingObs);
    return dirichletln(probs, obs, s) - dirichletln(probs, m, s);
}
