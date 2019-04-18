#include "Multinomial.h"
#include "Sum.h"
#include "Product.h"


long double multinomialSamplingProb(const vector<long double>& probs, const vector<int>& obs) {
    vector<long double> factorials;
    vector<long double> probsPowObs;
    factorials.resize(obs.size());
    transform(obs.begin(), obs.end(), factorials.begin(), factorial);
    vector<long double>::const_iterator p = probs.begin();
    vector<int>::const_iterator o = obs.begin();
    for (; p != probs.end() && o != obs.end(); ++p, ++o) {
        probsPowObs.push_back(pow(*p, *o));
    }
    return factorial(sum(obs)) / product(factorials) * product(probsPowObs);
}

// TODO rename to reflect the fact that this is the multinomial sampling
// probability for obs counts given probs probabilities
long double multinomialSamplingProbLn(const vector<long double>& probs, const vector<int>& obs) {
    vector<long double> factorials;
    vector<long double> probsPowObs;
    factorials.resize(obs.size());
    transform(obs.begin(), obs.end(), factorials.begin(), factorialln);
    vector<long double>::const_iterator p = probs.begin();
    vector<int>::const_iterator o = obs.begin();
    for (; p != probs.end() && o != obs.end(); ++p, ++o) {
        probsPowObs.push_back(powln(log(*p), *o));
    }
    return factorialln(sum(obs)) - sum(factorials) + sum(probsPowObs);
}

long double multinomialCoefficientLn(int n, const vector<int>& counts) {
    vector<long double> count_factorials;
    count_factorials.resize(counts.size());
    transform(counts.begin(), counts.end(), count_factorials.begin(), factorialln);
    return factorialln(n) - sum(count_factorials);
}

long double samplingProbLn(const vector<long double>& probs, const vector<int>& obs) {
    vector<long double>::const_iterator p = probs.begin();
    vector<int>::const_iterator o = obs.begin();
    long double r = 0;
    for (; p != probs.end() && o != obs.end(); ++p, ++o) {
        r += powln(log(*p), *o);
    }
    return r;
}
