#include <map>
#include <cmath>
#include "Utility.h"

using namespace std;

// genotype priors

long double alleleFrequencyProbability(const map<int, int>& alleleFrequencyCounts, long double theta);
long double alleleFrequencyProbabilityln(const map<int, int>& alleleFrequencyCounts, long double theta);
long double __alleleFrequencyProbabilityln(const map<int, int>& alleleFrequencyCounts, long double theta);

class AlleleFrequencyProbabilityCache : public map<map<int, int>, long double> {
public:
    long double alleleFrequencyProbabilityln(const map<int, int>& counts, long double theta) {
        map<map<int, int>, long double>::iterator p = find(counts);
        if (p == end()) {
            long double pln = __alleleFrequencyProbabilityln(counts, theta);
            insert(make_pair(counts, pln));
            return pln;
        } else {
            return p->second;
        }
    }
};
