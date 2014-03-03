#include "Ewens.h"


long double alleleFrequencyProbability(const map<int, int>& alleleFrequencyCounts, long double theta) {

    int M = 0;
    long double p = 1;

    for (map<int, int>::const_iterator f = alleleFrequencyCounts.begin(); f != alleleFrequencyCounts.end(); ++f) {
        int frequency = f->first;
        int count = f->second;
        M += frequency * count;
        p *= (double) pow((double) theta, (double) count) / ((double) pow((double) frequency, (double) count) * factorial(count));
    }

    long double thetaH = 1;
    for (int h = 1; h < M; ++h)
        thetaH *= theta + h;

    return factorial(M) / (theta * thetaH) * p;

}

AlleleFrequencyProbabilityCache alleleFrequencyProbabilityCache;

long double alleleFrequencyProbabilityln(const map<int, int>& alleleFrequencyCounts, long double theta) {
    return alleleFrequencyProbabilityCache.alleleFrequencyProbabilityln(alleleFrequencyCounts, theta);
}

// Implements Ewens' Sampling Formula, which provides probability of a given
// partition of alleles in a sample from a population
long double __alleleFrequencyProbabilityln(const map<int, int>& alleleFrequencyCounts, long double theta) {

    int M = 0; // multiplicity of site
    long double p = 0;
    long double thetaln = log(theta);

    for (map<int, int>::const_iterator f = alleleFrequencyCounts.begin(); f != alleleFrequencyCounts.end(); ++f) {
        int frequency = f->first;
        int count = f->second;
        M += frequency * count;
        p += powln(thetaln, count) - (powln(log(frequency), count) + factorialln(count));
    }

    long double thetaH = 0;
    for (int h = 1; h < M; ++h)
        thetaH += log(theta + h);

    return factorialln(M) - (thetaln + thetaH) + p;

}
