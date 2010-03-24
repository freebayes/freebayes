// utility functions
//
#include "Utility.h"

short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}

double phred2float(int qual) {
    return pow(10,qual * -.1);
}

short float2phred(double prob) {
    return -10 * log10(prob);
}

// the probability that we have a completely true vector of qualities
short jointQuality(const std::vector<short>& quals) {
    std::vector<double> probs;
    for (int i = 0; i<quals.size(); ++i) {
        probs.push_back(phred2float(quals[i]));
    }
    // product of probability we don't have a true event for each element
    double prod = 1 - probs.front();
    for (int i = 1; i<probs.size(); ++i) {
        prod *= 1 - probs.at(i);
    }
    // and then invert it again to get probability of an event
    return float2phred(1 - prod);
}

