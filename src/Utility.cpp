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

// TODO do the following in phred (log) space for perf boost

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

short jointQuality(const std::string& qualstr) {

    std::vector<short> quals;
    for (int i=0; i<qualstr.size(); i++)
        quals.push_back(qualityChar2ShortInt(qualstr.at(i)));

    return jointQuality(quals);

}

short averageQuality(const std::string& qualstr) {

    double q = 0;
    for (int i=0; i<qualstr.size(); i++)
        q += phred2float(qualityChar2ShortInt(qualstr.at(i)));
    return float2phred(q / qualstr.size());

}

unsigned int factorial(int n) 
{
    int f = 1;
    for(int i=1; i<=n; i++) 
    {
        f *= i;
    }
    return f;
}
