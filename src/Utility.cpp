// utility functions
//
#include "Utility.h"

using namespace std;

short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}

long double factorial(int n) {
    long double f = 1;
    while ( n > 1) {
        f *= n;
        n -= 1;
    }
    return f;
}

long double phred2float(int qual) {
    return pow(10,qual * -.1);
}

int float2phred(long double prob) {
    return std::min(-10 * (long double) log10(prob), (long double) 99);
}

// TODO do the following in phred (log) space for perf boost

// the probability that we have a completely true vector of qualities
short jointQuality(const std::vector<short>& quals) {
    std::vector<long double> probs;
    for (int i = 0; i<quals.size(); ++i) {
        probs.push_back(phred2float(quals[i]));
    }
    // product of probability we don't have a true event for each element
    long double prod = 1 - probs.front();
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

    long double q = 0;
    for (int i=0; i<qualstr.size(); i++)
        q += phred2float(qualityChar2ShortInt(qualstr.at(i)));
    return float2phred(q / qualstr.size());

}

/*
unsigned int factorial(int n) 
{
    int f = 1;
    for(int i=1; i<=n; i++) 
    {
        f *= i;
    }
    return f;
}
*/

bool stringInVector(string item, vector<string> items) {
    for (vector<string>::iterator i = items.begin(); i != items.end(); ++i) {
        if (item == *i) {
            return true;
        }
    }
    return false;
}

