#include "Multinomial.h"
#include "Sum.h"
#include "Product.h"


long double multinomial(vector<long double> probs, vector<int> obs) {
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

long double multinomialln(vector<long double> probs, vector<int> obs) {
    vector<long double> factorials;
    vector<long double> probsPowObs;
    factorials.resize(obs.size());
    transform(obs.begin(), obs.end(), factorials.begin(), factorialln);
    vector<long double>::const_iterator p = probs.begin();
    vector<int>::const_iterator o = obs.begin();
    for (; p != probs.end() && o != obs.end(); ++p, ++o) {
        probsPowObs.push_back(log(pow(*p, *o)));
    }
    return factorialln(sum(obs)) - sum(factorials) + sum(probsPowObs);
}
